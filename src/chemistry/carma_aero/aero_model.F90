!===============================================================================
! CAMRA Aerosol Model
!===============================================================================
module aero_model
  use physics_buffer,    only: physics_buffer_desc, pbuf_get_index, pbuf_add_field, dtype_r8
  use shr_kind_mod,      only: r8 => shr_kind_r8
  use constituents,      only: pcnst, cnst_name, cnst_get_ind
  use perf_mod,          only: t_startf, t_stopf
  use ppgrid,            only: pcols, pver, pverp
  use phys_control,      only: phys_getopts, cam_physpkg_is
  use cam_abortutils,    only: endrun
  use cam_logfile,       only: iulog
  use physics_types,     only: physics_state, physics_ptend, physics_ptend_init
  use camsrfexch,        only: cam_in_t, cam_out_t
  use physics_buffer,    only: pbuf_get_field, pbuf_get_index, pbuf_set_field, dtype_r8
  use physconst,         only: gravit, rair, rhoh2o
  use spmd_utils,        only: masterproc
  use cam_history,       only: outfld
  use infnan,            only: nan, assignment(=)
  use rad_constituents,only: rad_cnst_get_info, rad_cnst_get_info_by_bin, &
                             rad_cnst_get_info_by_bin_spec, rad_cnst_get_bin_props_by_idx, &
                             rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_mmr, &
                             rad_cnst_get_bin_num

  implicit none
  private

  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch ! create, grow, change, and shrink aerosols.
  public :: aero_model_drydep     ! aerosol dry deposition and sediment
  public :: aero_model_wetdep     ! aerosol wet removal
  public :: aero_model_emissions  ! aerosol emissions
  public :: aero_model_surfarea    ! tropospheric aerosol wet surface area for chemistry
  public :: aero_model_strat_surfarea   ! stub

   ! Misc private data
  character(len=32), allocatable :: fieldname(:)    ! names for interstitial output fields
  character(len=32), allocatable :: fieldname_cw(:)    ! names for cloud_borne output fields

  ! number of modes
  integer :: wetdens_ap_idx      = 0
  integer :: qaerwat_idx         = 0

  integer :: fracis_idx          = 0
  integer :: prain_idx           = 0
  integer :: rprddp_idx          = 0
  integer :: rprdsh_idx          = 0
  integer :: nevapr_shcu_idx     = 0
  integer :: nevapr_dpcu_idx     = 0

  integer :: sulfeq_idx = -1

  integer :: nh3_ndx    = 0
  integer :: nh4_ndx    = 0

  ! variables for table lookup of aerosol impaction/interception scavenging rates
  integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12
  real(r8) :: dlndg_nimptblgrow
  real(r8),allocatable :: scavimptblnum(:,:)
  real(r8),allocatable :: scavimptblvol(:,:)


  ! description of bin aerosols
  integer, public, protected :: nspec_max = 0
  integer, public, protected :: nbins = 0
  integer, public, protected, allocatable :: nspec(:)

  ! local indexing for bins
  integer, allocatable :: bin_idx(:,:) ! table for local indexing of modal aero number and mmr
  integer :: ncnst_tot                  ! total number of mode number conc + mode species
  integer :: ncnst_extd                  ! twiece total number of mode number conc + mode species

  ! Indices for CARMA species in the ptend%q array.  Needed for prognostic aerosol case.
  logical, allocatable :: bin_cnst_lq(:,:)
  integer, allocatable :: bin_cnst_idx(:,:)


  ! ptr2d_t is used to create arrays of pointers to 2D fields
  type ptr2d_t
    real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
                               ! in the ptend object

    ! Namelist variables
  real(r8)          :: sol_facti_cloud_borne   = 1._r8
  real(r8)          :: sol_factb_interstitial  = 0.1_r8
  real(r8)          :: sol_factic_interstitial = 0.4_r8
  real(r8)          :: seasalt_emis_scale

 logical :: convproc_do_aer



contains

  !=============================================================================
  ! reads aerosol namelist options
  !=============================================================================
  subroutine aero_model_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use carma_aero_convproc,   only: ma_convproc_readnl
    !st use dust_model,      only: dust_readnl

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aero_model_readnl'

    ! Namelist variables
    !st character(len=16) :: aer_wetdep_list(pcnst) = ' '
    !st character(len=16) :: aer_drydep_list(pcnst) = ' '

    namelist /aerosol_nl/ sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aerosol_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aerosol_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    !st call mpibcast(aer_wetdep_list,   len(aer_wetdep_list(1))*pcnst, mpichar, 0, mpicom)
    !st call mpibcast(aer_drydep_list,   len(aer_drydep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast(sol_facti_cloud_borne, 1,                         mpir8,   0, mpicom)
    call mpibcast(sol_factb_interstitial, 1,                        mpir8,   0, mpicom)
    call mpibcast(sol_factic_interstitial, 1,                       mpir8,   0, mpicom)
    !st call mpibcast(modal_strat_sulfate,     1,                       mpilog,  0, mpicom)
    !st call mpibcast(seasalt_emis_scale, 1,                            mpir8,   0, mpicom)
    !st call mpibcast(modal_accum_coarse_exch, 1,                       mpilog,  0, mpicom)
#endif

    !st wetdep_list = aer_wetdep_list
    !st drydep_list = aer_drydep_list

    call ma_convproc_readnl(nlfile)
    !st call dust_readnl(nlfile)


  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register()


    integer :: m, l
    character(len=32) :: spec_name
    character(len=32) :: mmr_name

    integer :: idx

    call rad_cnst_get_info( 0, nbins=nbins)
    allocate( nspec(nbins) )

    ! add pbuf fields for interstitial (cloud borne) aerosols in CARMA
    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=nspec(m), mmr_name=mmr_name)
       call pbuf_add_field('CLD'//trim(mmr_name),'global',dtype_r8,(/pcols,pver/), idx)
       call pbuf_add_field('CLDNB'//trim(mmr_name),'global',dtype_r8,(/pcols,pver/), idx)
       do l = 1, nspec(m)
          call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=spec_name)
          call pbuf_add_field('CLD'//trim(spec_name),'global',dtype_r8,(/pcols,pver/),idx)
       enddo
    enddo

  end subroutine aero_model_register

  !=============================================================================
  !=============================================================================
  subroutine aero_model_init( pbuf2d )

    use mo_chem_utls,    only: get_inv_ndx
    use cam_history,     only: addfld, add_default, horiz_only
    use mo_chem_utls,    only: get_rxt_ndx, get_spc_ndx
    !st use modal_aero_data, only: cnst_name_cw
    !st use modal_aero_data, only: modal_aero_data_init
    !st use dust_model,      only: dust_init, dust_names, dust_active, dust_nbin, dust_nnum
    !st use seasalt_model,   only: seasalt_init, seasalt_names, seasalt_active,seasalt_nbin
    !st use drydep_mod,      only: inidrydep
    use wetdep,          only: wetdep_init
    !st use mo_setsox,       only: sox_inti

    !st use modal_aero_calcsize,   only: modal_aero_calcsize_init
    !st use modal_aero_coag,       only: modal_aero_coag_init
    !st use modal_aero_deposition, only: modal_aero_deposition_init
    !st use modal_aero_gasaerexch, only: modal_aero_gasaerexch_init
    !st use modal_aero_newnuc,     only: modal_aero_newnuc_init
    !st use modal_aero_rename,     only: modal_aero_rename_init

    use carma_aero_convproc,   only: ma_convproc_init

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    ! local vars
    character(len=*), parameter :: subrname = 'aero_model_init'
    integer :: m, n, id, ii, mm
    integer :: lptr      = -1
    integer :: idxtmp    = -1
    character(len=20) :: dummy

    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies
    logical  :: history_chemistry, history_cesm_forcing, history_dust

    integer :: l
    character(len=6) :: test_name
    character(len=64) :: errmes

    character(len=2)  :: unit_basename  ! Units 'kg' or '1'
    integer :: errcode
    !st character(len=fieldname_len) :: field_name

    ! aqueous chem initialization
    !st call sox_inti()

    fracis_idx      = pbuf_get_index('FRACIS')
    prain_idx       = pbuf_get_index('PRAIN')
    rprddp_idx      = pbuf_get_index('RPRDDP')
    rprdsh_idx      = pbuf_get_index('RPRDSH')
    nevapr_shcu_idx = pbuf_get_index('NEVAPR_SHCU')
    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
    !st sulfeq_idx      = pbuf_get_index('MAMH2SO4EQ',errcode)

    !st not sure if this is needed
    call phys_getopts(history_aerosol_out = history_aerosol, &
                      history_chemistry_out=history_chemistry, &
                      history_cesm_forcing_out=history_cesm_forcing, &
                      convproc_do_aer_out = convproc_do_aer)

    !st call modal_aero_data_init(pbuf2d)
    !st call modal_aero_bcscavcoef_init()

    !st  call modal_aero_rename_init( modal_accum_coarse_exch )
    !   calcsize call must follow rename call
    !st call modal_aero_calcsize_init( pbuf2d )
    !st call modal_aero_gasaerexch_init
    !   coag call must follow gasaerexch call
    !st call modal_aero_coag_init
    !st call modal_aero_newnuc_init

    ! call modal_aero_deposition_init only if the user has not specified
    ! prescribed aerosol deposition fluxes
    !st if (.not.aerodep_flx_prescribed()) then
    !st   call modal_aero_deposition_init
    !stendif

    if (convproc_do_aer) then
       call ma_convproc_init()
    endif

    !st call dust_init()
    !st call seasalt_init(seasalt_emis_scale)
    call wetdep_init()

    !st all CARMA species are deposited, therefore the following is not used
    !st nwetdep = 0
    !st ndrydep = 0

    !st count_species: do m = 1,pcnst
    !st    if ( len_trim(wetdep_list(m)) /= 0 ) then
    !st       nwetdep = nwetdep+1
    !st    endif
    !st    if ( len_trim(drydep_list(m)) /= 0 ) then
    !st       ndrydep = ndrydep+1
    !st    endif
    !st enddo count_species

    ! add plus one to include number, total mmr and nspec
    nspec_max = maxval(nspec) + 2

    ncnst_tot = nspec(1) + 2
    do m = 2, nbins
      ncnst_tot = ncnst_tot + nspec(m) + 2
    end do
    ncnst_extd = 2*ncnst_tot

    allocate( &
      bin_idx(nbins,nspec_max),      &
      bin_cnst_lq(nbins,nspec_max), &
      bin_cnst_idx(nbins,nspec_max), &
      fieldname_cw(ncnst_tot), &
      fieldname(ncnst_tot) )

    ii = 0
    do m = 1, nbins
      do l = 1, nspec(m) + 2    ! do through nspec plus mmr and number
         ii = ii + 1
         bin_idx(m,l) = ii

         if (l <= nspec(m) ) then   ! species
            call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=fieldname(ii), spec_name_cw=fieldname_cw(ii))
         else if (l == nspec(m) + 1) then   ! mmr
            call rad_cnst_get_info_by_bin(0, m,  mmr_name=fieldname(ii), mmr_name_cw=fieldname_cw(ii))
         else if (l == nspec(m) + 2) then   !number
            call rad_cnst_get_info_by_bin(0, m, num_name=fieldname(ii), num_name_cw=fieldname_cw(ii))
         end if

         call cnst_get_ind(fieldname(ii), idxtmp, abort=.false.)
          if (idxtmp.gt.0) then
             bin_cnst_lq(m,l) = .true.
             bin_cnst_idx(m,l) = idxtmp
             lq(idxtmp) = .true.
          else
             bin_cnst_lq(m,l) = .false.
             bin_cnst_idx(m,l) = 0
          end if

         mm = ii

         unit_basename = 'kg'
         if (l == nspec(m) + 2) then   ! number
          unit_basename = ' 1'
         end if

      if (l == nspec(m) + 1) then   ! mmr
       call addfld (trim(fieldname(mm))//'SFWET', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux at surface')
       call addfld (trim(fieldname(mm))//'SFSIC', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux (incloud, convective) at surface')
       call addfld (trim(fieldname(mm))//'SFSIS', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux (incloud, stratiform) at surface')
       call addfld (trim(fieldname(mm))//'SFSBC', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux (belowcloud, convective) at surface')
       call addfld (trim(fieldname(mm))//'SFSBS', &
            horiz_only,  'A',unit_basename//'/m2/s ','Wet deposition flux (belowcloud, stratiform) at surface')

       if (convproc_do_aer) then
          call addfld (trim(fieldname(mm))//'SFSES', &
               horiz_only,  'A','kg/m2/s','Wet deposition flux (precip evap, stratiform) at surface')
          call addfld (trim(fieldname(mm))//'SFSBD', &
               horiz_only,  'A','kg/m2/s','Wet deposition flux (belowcloud, deep convective) at surface')
       end if

       call addfld (trim(fieldname(mm))//'WET',(/ 'lev' /), 'A',unit_basename//'/kg/s ','wet deposition tendency')
       call addfld (trim(fieldname(mm))//'SIC',(/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(fieldname(mm))//' ic wet deposition')
       call addfld (trim(fieldname(mm))//'SIS',(/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(fieldname(mm))//' is wet deposition')
       call addfld (trim(fieldname(mm))//'SBC',(/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(fieldname(mm))//' bc wet deposition')
       call addfld (trim(fieldname(mm))//'SBS',(/ 'lev' /), 'A',unit_basename//'/kg/s ', &
            trim(fieldname(mm))//' bs wet deposition')

       if ( history_aerosol .or. history_chemistry ) then
          call add_default (trim(fieldname(mm))//'SFWET', 1, ' ')
       endif
       if ( history_aerosol ) then
          call add_default (trim(fieldname(mm))//'SFSIC', 1, ' ')
          call add_default (trim(fieldname(mm))//'SFSIS', 1, ' ')
          call add_default (trim(fieldname(mm))//'SFSBC', 1, ' ')
          call add_default (trim(fieldname(mm))//'SFSBS', 1, ' ')
          if (convproc_do_aer) then
             call add_default (trim(fieldname(mm))//'SFSES', 1, ' ')
             call add_default (trim(fieldname(mm))//'SFSBD', 1, ' ')
          end if
       endif

          call addfld( fieldname_cw(mm),                (/ 'lev' /), 'A', unit_basename//'/kg ',   &
               trim(fieldname_cw(mm))//' in cloud water')
          call addfld (trim(fieldname_cw(mm))//'SFWET', horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' wet deposition flux at surface')
          call addfld (trim(fieldname_cw(mm))//'SFSIC', horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' wet deposition flux (incloud, convective) at surface')
          call addfld (trim(fieldname_cw(mm))//'SFSIS', horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' wet deposition flux (incloud, stratiform) at surface')
          call addfld (trim(fieldname_cw(mm))//'SFSBC', horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' wet deposition flux (belowcloud, convective) at surface')
          call addfld (trim(fieldname_cw(mm))//'SFSBS', horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' wet deposition flux (belowcloud, stratiform) at surface')
          call addfld (trim(fieldname_cw(mm))//'DDF',   horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' dry deposition flux at bottom (grav + turb)')
          call addfld (trim(fieldname_cw(mm))//'TBF',   horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' turbulent dry deposition flux')
          call addfld (trim(fieldname_cw(mm))//'GVF',   horiz_only,  'A', unit_basename//'/m2/s ', &
               trim(fieldname_cw(mm))//' gravitational dry deposition flux')

          if (convproc_do_aer) then
             call addfld (trim(fieldname_cw(mm))//'SFSEC', &
                horiz_only,  'A','kg/m2/s','Wet deposition flux (precip evap, convective) at surface')
             call addfld (trim(fieldname_cw(mm))//'SFSES', &
                horiz_only,  'A','kg/m2/s','Wet deposition flux (precip evap, stratiform) at surface')
             call addfld (trim(fieldname_cw(mm))//'SFSBD', &
                horiz_only,  'A','kg/m2/s','Wet deposition flux (belowcloud, deep convective) at surface')
          end if


          if ( history_aerosol.or. history_chemistry ) then
             call add_default( fieldname_cw(mm), 1, ' ' )
             call add_default (trim(fieldname_cw(mm))//'SFWET', 1, ' ')
          endif
          if ( history_aerosol ) then
             call add_default (trim(fieldname_cw(mm))//'GVF', 1, ' ')
             call add_default (trim(fieldname_cw(mm))//'TBF', 1, ' ')
             call add_default (trim(fieldname_cw(mm))//'DDF', 1, ' ')
             call add_default (trim(fieldname_cw(mm))//'SFSBS', 1, ' ')
             call add_default (trim(fieldname_cw(mm))//'SFSIC', 1, ' ')
             call add_default (trim(fieldname_cw(mm))//'SFSBC', 1, ' ')
             call add_default (trim(fieldname_cw(mm))//'SFSIS', 1, ' ')
             if (convproc_do_aer) then
                call add_default (trim(fieldname_cw(mm))//'SFSEC', 1, ' ')
                call add_default (trim(fieldname_cw(mm))//'SFSES', 1, ' ')
                call add_default (trim(fieldname_cw(mm))//'SFSBD', 1, ' ')
             end if
          endif
         endif  ! only for nspec + 1
       enddo
    enddo



  end subroutine aero_model_init

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    ! args
    type(physics_state),    intent(in)    :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)
    real(r8),               intent(in)    :: ustar(:)  ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt        ! time step
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

  endsubroutine aero_model_drydep

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep( state, dt, dlf, cam_out, ptend, pbuf)

    use modal_aero_deposition, only: set_srf_wetdep
    use wetdep,                only: wetdepa_v2, wetdep_inputs_set, wetdep_inputs_t
    !st use modal_aero_wateruptake,only: modal_aero_wateruptake_dr
    use carma_aero_convproc,   only: deepconv_wetdep_history, ma_convproc_intr, convproc_do_evaprain_atonce


    ! args

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! local vars

    integer :: m ! tracer index

    integer :: lchnk ! chunk identifier
    integer :: ncol ! number of atmospheric columns

    real(r8) :: iscavt(pcols, pver)

    integer :: mm
    integer :: i,k


    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)
    !st real(r8) :: sol_factb, sol_facti
    real(r8) :: sol_factb(pcols,pver)
    real(r8) :: sol_facti(pcols,pver)
    real(r8) :: sol_factic(pcols,pver)

    real(r8) :: sflx(pcols) ! deposition flux

    integer :: jnv ! index for scavcoefnv 3rd dimension
    integer :: lphase ! index for interstitial / cloudborne aerosol
    integer :: strt_loop, end_loop, stride_loop !loop indices for the lphase loop
    integer :: l ! index for aerosol number / chem-mass / water-mass
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species
    real(r8) :: dqdt_tmp_sum(pcols,pver) ! temporary array to hold sum of tendency for all species
    real(r8) :: frac_dqdt(pcols,pver) ! temporary array to hold sum fraction between dqdt_tmp_sum and mmr per bin
    real(r8) :: f_act_conv(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: fracis_cw(pcols,pver)
    !st real(r8) :: dryr(pcols,pver)       ! dry radius from CARMA
    real(r8) :: stoke(pcols,pver)       ! dry radius from CARMA
    real(r8) :: J(pcols,pver)       ! dry radius from CARMA
    real(r8) :: prec(pcols) ! precipitation rate
    real(r8) :: q_tmp(pcols,pver) ! temporary array to hold "most current" mixing ratio for 1 species
    real(r8) :: scavcoefnv(pcols,pver,0:2) ! Dana and Hales coefficient (/mm) for
                                           ! cloud-borne num & vol (0),
                                           ! MAM: interstitial num (1), interstitial vol (2)
                                           ! CARMA does not use calculate scavcoefnv for num separately: interstitial num and vol (1)
    real(r8) :: tmpa, tmpb
    real(r8) :: tmpdust, tmpnacl
    real(r8) :: water_old, water_new ! temporary old/new aerosol water mix-rat
    logical  :: isprx(pcols,pver) ! true if precipation
    real(r8) :: aerdepwetis(pcols,ncnst_tot) ! aerosol wet deposition (interstitial)
    real(r8) :: aerdepwetcw(pcols,ncnst_tot) ! aerosol wet deposition (cloud water)

    ! For unified convection scheme
    !st logical, parameter :: do_aero_water_removal = .false. ! True if aerosol water reduction by wet removal is to be calculated
    !st                                                       ! (this has not been fully tested, so best to leave it off)
    logical :: do_lphase1, do_lphase2

    real(r8), pointer :: rprddp(:,:)     ! rain production, deep convection
    real(r8), pointer :: rprdsh(:,:)     ! rain production, shallow convection
    real(r8), pointer :: evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.
    real(r8), pointer :: evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.

    type(ptr2d_t), allocatable :: raer(:)     ! aerosol mass, number mixing ratios
    type(ptr2d_t), allocatable :: qqcw(:)


    real(r8) :: rprddpsum(pcols)
    real(r8) :: rprdshsum(pcols)
    real(r8) :: evapcdpsum(pcols)
    real(r8) :: evapcshsum(pcols)

    real(r8) :: tmp_resudp, tmp_resush
    real(r8) :: totalmmr
    real(r8) :: solmmr
    real(r8) :: solfact
    real(r8) :: rho_water
    real(r8) :: specdens              ! specie density from physprops files

    real(r8) :: sflxec(pcols), sflxecdp(pcols)  ! deposition flux
    real(r8) :: sflxic(pcols), sflxicdp(pcols)  ! deposition flux
    real(r8) :: sflxbc(pcols), sflxbcdp(pcols)  ! deposition flux
    real(r8) :: rcscavt(pcols, pver)
    real(r8) :: rsscavt(pcols, pver)
    real(r8) :: qqcw_in(pcols,pver), qqcw_sav(pcols,pver,0:nspec_max) ! temporary array to hold qqcw for the current mode
    real(r8) :: fldcw(pcols,pver)
    real(r8) :: rtscavt(pcols, pver, 0:nspec_max)

    integer, parameter :: nsrflx_mzaer2cnvpr = 2
    real(r8) :: qsrflx_mzaer2cnvpr(pcols,ncnst_tot,nsrflx_mzaer2cnvpr)
    ! End unified convection scheme


    !st real(r8), pointer :: dgnumwet(:,:,:)
    !st real(r8), pointer :: qaerwat(:,:,:)  ! aerosol water

    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble

    real(r8), pointer :: dryr(:,:)   ! CARMA dry radius in cm

    type(wetdep_inputs_t) :: dep_inputs

    real(r8) :: dcondt_resusp3d(ncnst_extd,pcols, pver)

    character(len=32) :: spectype
    character(len=32) :: bin_name
    character(len=2) :: tmpbin2

    lchnk = state%lchnk
    ncol  = state%ncol

    dcondt_resusp3d(:,:,:) = 0._r8

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_wetdep', lq=lq)

!st  CARMA is doing water uptake (in the CARMA code), we leave this out here
    ! Do calculations of mode radius and water uptake if:
    ! 1) modal aerosols are affecting the climate, or
    ! 2) prognostic modal aerosols are enabled

    !st call t_startf('calcsize')  not needed for CARMA
    ! for prognostic modal aerosols the transfer of mass between aitken and accumulation
    ! modes is done in conjunction with the dry radius calculation
    !st call modal_aero_calcsize_sub(state, ptend, dt, pbuf)
    !st call t_stopf('calcsize')

    !st call t_startf('wateruptake')
    !st call modal_aero_wateruptake_dr(state, pbuf)
    !st call t_stopf('wateruptake')

    !NOTE: currently all aerosols in CARMA are wet deposited, namelist setting or master table of ndep could be included later
    !st if (nwetdep<1) return

    call wetdep_inputs_set( state, pbuf, dep_inputs )

    allocate( &
      raer(ncnst_tot),                &
      qqcw(ncnst_tot)                 )


    !st needed to calcuate modal_aero_bcscavcoef_get scavcoefvol, scavcoefnum (need to get CARMA bin wet radius)

    !st call pbuf_get_field(pbuf, qaerwat_idx,        qaerwat,  start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
    call pbuf_get_field(pbuf, fracis_idx,         fracis, start=(/1,1,1/), kount=(/pcols, pver, ncnst_tot/) )

    prec(:ncol)=0._r8
    do k=1,pver
       where (prec(:ncol) >= 1.e-7_r8)
          isprx(:ncol,k) = .true.
       elsewhere
          isprx(:ncol,k) = .false.
       endwhere
       prec(:ncol) = prec(:ncol) + (dep_inputs%prain(:ncol,k) + dep_inputs%cmfdqr(:ncol,k) - dep_inputs%evapr(:ncol,k)) &
            *state%pdel(:ncol,k)/gravit
    end do

    if(convproc_do_aer) then
       qsrflx_mzaer2cnvpr(:,:,:) = 0.0_r8
       aerdepwetis(:,:)          = 0.0_r8
       aerdepwetcw(:,:)          = 0.0_r8
    else
       qsrflx_mzaer2cnvpr(:,:,:) = nan
       aerdepwetis(:,:)          = nan
       aerdepwetcw(:,:)          = nan
    endif

    scavcoefnv(:,:,0) = 0.0_r8 ! below-cloud scavcoef = 0.0 for cloud-borne species

    ! Counters for "without" unified convective treatment (i.e. default case)
    strt_loop   = 1
    end_loop    = 2
    stride_loop = 1
    if (convproc_do_aer) then
       !Do cloudborne first for unified convection scheme so that the resuspension of cloudborne
       !can be saved then applied to interstitial
       strt_loop   =  2
       end_loop    =  1
       stride_loop = -1
    endif

    do m = 1, nbins      ! main loop over aerosol bins
         ! r(m) is the dry bin radius
         ! taken here from CARMA pbuf field
         ! get bin info
         call rad_cnst_get_info_by_bin(0, m, bin_name=bin_name)
         call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_dryr"),dryr)
         !st write(iulog,*) 'bin_name ',bin_name

         ! Init pointers to mode number and specie mass mixing ratios in
         ! intersitial and cloud borne phases.
         do l = 1, nspec(m) + 2
            mm = bin_idx(m, l)
            if (l <= nspec(m)) then
              call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'a', state, pbuf, raer(mm)%fld)
              call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
            end if
            if (l == nspec(m)+1) then
              call rad_cnst_get_bin_mmr(0, m, 'a', state, pbuf, raer(mm)%fld)
              call rad_cnst_get_bin_mmr(0, m, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
            end if
            if (l == nspec(m)+2) then
              call rad_cnst_get_bin_num(0, m, 'a', state, pbuf, raer(mm)%fld)
              call rad_cnst_get_bin_num(0, m, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
            end if
         end do

         ! get  solfac
         ! note: solfact is for particle (group bin) not for element, need to get from CARMA, make namelist variable
         ! Francis, could  you add namelist variables that set solfact for pure sulates 'PR*" and mixed groups 'MX*"?
         tmpbin2  = trim(adjustl(bin_name(2:)))
         if (tmpbin2 == 'PR') then
               solfact = 0.3
         end if
         if (tmpbin2 == 'MX') then
              solfact = 0.2
         end if

         !--------------------------------------------------------------
         !sol_factb(:ncol,:)=solfact
         ! get mmr for each specie and multiply by the factors to derive sol_factb
         sol_factb(:ncol,:)=0._r8

           do k= 1,pver
             do i= 1,ncol
               totalmmr = 0._r8
               solmmr = 0._r8
               do l = 1, nspec(m)
                     mm = bin_idx(m, l)
                     totalmmr    = totalmmr + raer(mm)%fld(i,k)

                     call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)

                     if (trim(spectype) == 'sulfate') then
                        solmmr =  solmmr + raer(mm)%fld(i,k)*0.5
                     end if
                     if (trim(spectype) == 'black-c') then
                        solmmr =  solmmr + raer(mm)%fld(i,k)*0.1
                     end if
                     if (trim(spectype) == 'p-organic') then
                        solmmr =  solmmr + raer(mm)%fld(i,k)*0.2
                     end if
                     if (trim(spectype) == 'dust') then
                        solmmr =  solmmr + raer(mm)%fld(i,k)*0.1
                     end if
                     if (trim(spectype) == 'seasalt') then
                        solmmr =  solmmr + raer(mm)%fld(i,k)*0.1
                     end if
               !st write(iulog,*) 'spectype, totalmmr,solmmr   = ',spectype, totalmmr, solmmr
               end do
               if (totalmmr .gt. 0._r8) then
                 sol_factb(i,k) = solmmr/totalmmr
               !st  write(iulog,*) 'i,k, sol_factb first  = ', i,k, sol_factb(i,k)
               end if
               sol_factb(i,k) = min(max(sol_factb(i,k), 0.1_r8), 0.8_r8)
             end do
           end do

         ! the check has been put in the loop
         !st sol_factb = min(max(sol_factb, 0.1_r8), 0.8_r8)


       do lphase = strt_loop,end_loop, stride_loop ! loop over interstitial (1) and cloud-borne (2) forms

          dqdt_tmp_sum(:,:) = 0.0_r8
          frac_dqdt(:,:) = 0.0_r8

          ! sol_factb and sol_facti values
          ! sol_factb - currently this is basically a tuning factor
          ! sol_facti & sol_factic - currently has a physical basis, and reflects activation fraction
          !
          ! 2008-mar-07 rce - sol_factb (interstitial) changed from 0.3 to 0.1
          ! was done for MAM
          ! - sol_factic (interstitial, dust modes) changed from 1.0 to 0.5
          ! - sol_factic (cloud-borne, pcarb modes) no need to set it to 0.0
          ! because the cloud-borne pcarbon == 0 (no activation)
          !
          ! rce 2010/05/02
          ! prior to this date, sol_factic was used for convective in-cloud wet removal,
          ! and its value reflected a combination of an activation fraction (which varied between modes)
          ! and a tuning factor
          ! from this date forward, two parameters are used for convective in-cloud wet removal
          ! f_act_conv is the activation fraction
          ! note that "non-activation" of aerosol in air entrained into updrafts should
          ! be included here
          ! eventually we might use the activate routine (with w ~= 1 m/s) to calculate
          ! this, but there is still the entrainment issue
          ! sol_factic is strictly a tuning factor
          !

        do l = 1, nspec(m) + 2  ! loop over different npsec plus over total mmr and number
            mm = bin_idx(m,l)

         if (l <= nspec(m)) then  ! loop only over nspec
          call rad_cnst_get_bin_props_by_idx(0, m, l, density_aer=specdens)

          if (lphase == 1) then ! interstial aerosol

            ! This is calculation for scavcoefnv has been introducted by Pengfei Yu for CARMA and is required in particular for a sectional model
            ! The scavenging coefficients are calcuated as a function of the aerosol bin at each
            ! grid point (lat,lon,vertical). Refer to Dana and Hales (1975), however I use
            ! Marshall and Palmer Distribution instead log-normal for rain drop size distribution
            ! Author: Pengfei Yu
            ! Aug.25, 2011


              rho_water = 1.e3_r8                                               ! kg/m3
              sol_factic(:,:) = 0.0_r8
              sol_facti(:,:) = 0._r8

              !st z_scavcoef = 0.1_r8
              scavcoefnv(:,:,1) =  0.1_r8
              scavcoefnv(:,:,2) =  0.1_r8

              ! convert precipitation prain (kg/kg/s) to mm/h

              J = (dep_inputs%prain(:,:)+dep_inputs%cmfdqr(:,:)-dep_inputs%evapc(:,:)-dep_inputs%evapr(:,:))  &
                                            *state%pdel/gravit/rho_water*3.6e6_r8   ! m/s
              !st write(iulog,*) 'm, l, J   = ', m, l, J
              !st write(iulog,*) 'm, l, dryr   = ', m, l,  dryr

              do k= 1,pver
                 do i= 1,ncol
                     stoke(i,k) = 0.1_r8*(dryr(i,k)**2)*1.0e8_r8*specdens                    ! stoke parameter
                     if ((J(i,k) .le. 0._r8)  .or. (dryr(i,k) .le. 0._r8))  then
                           scavcoefnv(:,:,1) =  0.1_r8
                           scavcoefnv(:,:,2) =  0.1_r8
                     else
                       if (dryr(i,k) .lt. 1.e-4_r8) then
                          !z_scavcoef(i_col,k_ver) = 63.03750_r8*r(ibin)*(J(i_col,k_ver)**(-0.42_r8))
                          ! scavenging volume
                           scavcoefnv(i,k,2) = 252.15_r8*dryr(i,k)*(J(i,k)**(-0.42_r8))
                          ! scavenging number
                           !st z_scavcoef(i,k,0) = 252.15_r8*r(ibin)*(J(i_col,k_ver)**(-0.42_r8))
                       else
                          !z_scavcoef(i_col,k_ver) = 1.025_r8*(J(i_col,k_ver)**(-0.21_r8))*(((stoke - 0.0833_r8)/(stoke + 0.5833_r8))**1.5_r8)
                          ! scavenging volume
                           scavcoefnv(i,k,2) = 2.05_r8*(J(i,k)**(-0.21_r8))*(((stoke(i,k) - 0.0833_r8)/(stoke(i,k) + 0.5833_r8))**1.5_r8)
                          ! scavenging  number
                           !st z_scavcoef(i,k,1) = 2.05_r8*(J(i_col,k_ver)**(-0.21_r8))*(((stoke - 0.0833_r8)/(stoke + 0.5833_r8))**1.5_r8)
                       endif
                     endif
                        ! commented out in Pengfei's code
!                       ! let big particle sea salt falls on quicker
!                       if (cam_in%ocnfrac(i_col)>0._r8 .and. r(ibin) .gt. 1.e-4_r8 .and. k_ver .gt. 40) then
!                          sol_factic(i_col,k_ver) = 0.8_r8
!                       endif
                    !st write(iulog,*) 'i,k, scavcoefnv(i,k,2)   = ', i,k,  scavcoefnv(i,k,2)
                    enddo
              enddo

             sol_factic = sol_factic_interstitial

!             if (m == modeptr_pcarbon) then
!                ! sol_factic = 0.0_r8 ! conv in-cloud scav OFF (0.0 activation fraction)
!                f_act_conv = 0.0_r8 ! rce 2010/05/02
!             else if ((m == modeptr_finedust) .or. (m == modeptr_coardust)) then
!                ! sol_factic = 0.2_r8 ! conv in-cloud scav ON (0.5 activation fraction) ! tuned 1/4
!                f_act_conv = 0.4_r8 ! rce 2010/05/02
!             else
                ! sol_factic = 0.4_r8 ! conv in-cloud scav ON (1.0 activation fraction) ! tuned 1/4
                f_act_conv = 0.8_r8 ! rce 2010/05/02
!             end if

          else ! cloud-borne aerosol (borne by stratiform cloud drops)

             sol_factb  = 0.0_r8   ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
             sol_facti  = sol_facti_cloud_borne   ! strat  in-cloud scav cloud-borne tuning factor  ! 1.0 in Pengfei's model
             sol_factic = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
                                   !        that conv precip collects strat droplets)
             f_act_conv = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean

          end if
          if (convproc_do_aer .and. lphase == 1) then
             ! if modal aero convproc is turned on for aerosols, then
             !    turn off the convective in-cloud removal for interstitial aerosols
             !    (but leave the below-cloud on, as convproc only does in-cloud)
             !    and turn off the outfld SFWET, SFSIC, SFSID, SFSEC, and SFSED calls
             ! for (stratiform)-cloudborne aerosols, convective wet removal
             !    (all forms) is zero, so no action is needed
             sol_factic = 0.0_r8
          end if

         end if   !(l <= nspec(m)) then


            if (lphase == 1) then
             jnv = 2
            else
             jnv = 1
            end if

            ! cloud-borne num & vol (0),
            ! interstitial num & vol (1), interstitial vol (2)

             ! calculate removal for mmr
             ! lphase = 1 interstitial  (jnv = 1 for number and mmr)
             if (lphase == 1) then

              if (l <= nspec(m)) then
                dqdt_tmp(:,:) = 0.0_r8
                ! q_tmp reflects changes from modal_aero_calcsize and is the "most current" q
                ! q_tmp(1:ncol,:) = state%q(1:ncol,:,mm) + ptend%q(1:ncol,:,mm)*dt
                q_tmp(1:ncol,:) = raer(mm)%fld

                if(convproc_do_aer) then
                   !Feed in the saved cloudborne mixing ratios from phase 2
                   qqcw_in(:,:) = qqcw_sav(:,:,l)
                else
                   fldcw = qqcw(mm)%fld
                   qqcw_in(:,:) = fldcw(:,:)
                endif

               ! do not calculate for number
                call wetdepa_v2( state%pmid, state%q(:,:,1), state%pdel, &
                     dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                     dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                     dep_inputs%evapr, dep_inputs%totcond, q_tmp, dt, &
                     dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                     dlf, fracis(:,:,mm), sol_factb, ncol, &
                     scavcoefnv(:,:,jnv), &
                     is_strat_cloudborne=.false.,  &
                     qqcw=qqcw_in(:,:),  &
                     f_act_conv=f_act_conv, &
                     icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                     convproc_do_aer=convproc_do_aer, rcscavt=rcscavt, rsscavt=rsscavt,  &
                     sol_facti_in=sol_facti, sol_factic_in=sol_factic, &
                     convproc_do_evaprain_atonce_in=convproc_do_evaprain_atonce )

                if(convproc_do_aer) then
                   ! add resuspension of cloudborne species to dqdt of interstitial species
                   dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) + rtscavt(1:ncol,:,l)
                endif
                !st write(iulog,*) 'm, l, dqdt_tmp   = ', m, l,  dqdt_tmp

               ! sum up dqdt_tmp for all species (not total mmr and number)
               dqdt_tmp_sum(1:ncol,:) = dqdt_tmp_sum(1:ncol,:) + dqdt_tmp(1:ncol,:)

                ! some CARMA interstetial species are not advected, check:
               if (bin_cnst_lq(m,l)) then ! adveced species
                  ptend%q(1:ncol,:,mm) = ptend%q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)
               else
                 raer(mm)%fld(1:ncol,:) = raer(mm)%fld(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt
               end if


              else if (l == nspec(m) + 1) then   !mmr and number are not advected

               do k= 1,pver
                 do i= 1,ncol
                   if (raer(mm)%fld(i,k) .gt. 0._r8) then
                         frac_dqdt(i,k) = dqdt_tmp_sum(i,k) / raer(mm)%fld(i,k)
                   else
                         frac_dqdt(i,k) = 0._r8
                   end if
                   !st write(iulog,*) 'm, l, i, k,  frac_dqdt(i,k), raer(mm)%fld(i,k) = ', frac_dqdt(i,k), raer(mm)%fld(i,k)
                 end do
               end do
               dqdt_tmp(1:ncol,:)  = dqdt_tmp_sum(1:ncol,:)
               raer(mm)%fld(1:ncol,:) = raer(mm)%fld(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

               !st write(iulog,*) 'nspec+1: m, l, dqdt_tmp   = ', m, l,  dqdt_tmp

              else if (l == nspec(m) + 2) then   !num fraction is assuming dN/N = dM/M
                !st write(iulog,*) 'nspec+2: mm, m, l, i, k,  frac_dqdt, raer(mm)%fld = ',mm, m, l, frac_dqdt, raer(mm)%fld

               dqdt_tmp(1:ncol,:) = frac_dqdt(1:ncol,:) * raer(mm)%fld(1:ncol,:)
               raer(mm)%fld(1:ncol,:) = raer(mm)%fld(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

              end if   ! nspec

             !st for testing may keep all output fields
              if (l == nspec(m) + 1) then   ! write output only for mmr  to reduce output
                !st write(iulog,*) 'outfld   = ', m, l,  mm, fieldname(mm)
                call outfld( trim(fieldname(mm))//'WET', dqdt_tmp(:,:), pcols, lchnk)
                call outfld( trim(fieldname(mm))//'SIC', icscavt, pcols, lchnk)
                call outfld( trim(fieldname(mm))//'SIS', isscavt, pcols, lchnk)
                call outfld( trim(fieldname(mm))//'SBC', bcscavt, pcols, lchnk)
                call outfld( trim(fieldname(mm))//'SBS', bsscavt, pcols, lchnk)
              end if

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
              if (l == nspec(m) + 1) then   !mmr and number are not advected
                if (.not.convproc_do_aer) call outfld( trim(fieldname(mm))//'SFWET', sflx, pcols, lchnk)
              end if
                aerdepwetis(:ncol,mm) = sflx(:ncol)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
               if (l == nspec(m) + 1) then   !mmr only
                if (.not.convproc_do_aer) call outfld( trim(fieldname(mm))//'SFSIC', sflx, pcols, lchnk)
               end if

                if (convproc_do_aer) sflxic = sflx

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
               if (l == nspec(m) + 1) then   !mmr only
                call outfld( trim(fieldname(mm))//'SFSIS', sflx, pcols, lchnk)
               end if

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
               if (l == nspec(m) + 1) then   !mmr only
                call outfld( trim(fieldname(mm))//'SFSBC', sflx, pcols, lchnk)
               end if
                if (convproc_do_aer)sflxbc = sflx

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
               if (l == nspec(m) + 1) then   !mmr only
                call outfld( trim(fieldname(mm))//'SFSBS', sflx, pcols, lchnk)
               end if

                if (convproc_do_aer) then

                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+rcscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   sflxec = sflx

                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+rsscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
               if (l == nspec(m) + 1) then   !mmr only
                   call outfld( trim(fieldname(mm))//'SFSES', sflx, pcols, lchnk)
               end if  ! nspec + 1 for history fields only


                   ! apportion convective surface fluxes to deep and shallow conv
                   ! this could be done more accurately in subr wetdepa
                   ! since deep and shallow rarely occur simultaneously, and these
                   !    fields are just diagnostics, this approximate method is adequate
                   ! only do this for interstitial aerosol, because conv clouds to not
                   !    affect the stratiform-cloudborne aerosol

                  if ( deepconv_wetdep_history) then

                      call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
                      call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
                      call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
                      call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )

                      rprddpsum(:)  = 0.0_r8
                      rprdshsum(:)  = 0.0_r8
                      evapcdpsum(:) = 0.0_r8
                      evapcshsum(:) = 0.0_r8

                      do k = 1, pver
                         rprddpsum(:ncol)  = rprddpsum(:ncol)  +  rprddp(:ncol,k)*state%pdel(:ncol,k)/gravit
                         rprdshsum(:ncol)  = rprdshsum(:ncol)  +  rprdsh(:ncol,k)*state%pdel(:ncol,k)/gravit
                         evapcdpsum(:ncol) = evapcdpsum(:ncol) + evapcdp(:ncol,k)*state%pdel(:ncol,k)/gravit
                         evapcshsum(:ncol) = evapcshsum(:ncol) + evapcsh(:ncol,k)*state%pdel(:ncol,k)/gravit
                      end do

                      do i = 1, ncol
                         rprddpsum(i)  = max( rprddpsum(i),  1.0e-35_r8 )
                         rprdshsum(i)  = max( rprdshsum(i),  1.0e-35_r8 )
                         evapcdpsum(i) = max( evapcdpsum(i), 0.1e-35_r8 )
                         evapcshsum(i) = max( evapcshsum(i), 0.1e-35_r8 )

                         ! assume that in- and below-cloud removal are proportional to column precip production
                         tmpa = rprddpsum(i) / (rprddpsum(i) + rprdshsum(i))
                         tmpa = max( 0.0_r8, min( 1.0_r8, tmpa ) )
                         sflxicdp(i) = sflxic(i)*tmpa
                         sflxbcdp(i) = sflxbc(i)*tmpa

                         ! assume that resuspension is proportional to (wet removal)*[(precip evap)/(precip production)]
                         tmp_resudp =           tmpa  * min( (evapcdpsum(i)/rprddpsum(i)), 1.0_r8 )
                         tmp_resush = (1.0_r8 - tmpa) * min( (evapcshsum(i)/rprdshsum(i)), 1.0_r8 )
                         tmpb = max( tmp_resudp, 1.0e-35_r8 ) / max( (tmp_resudp+tmp_resush), 1.0e-35_r8 )
                         tmpb = max( 0.0_r8, min( 1.0_r8, tmpb ) )
                         sflxecdp(i) = sflxec(i)*tmpb
                      end do
                    if (l == nspec(m) + 1) then   !mmr and number are not advected
                      call outfld( trim(fieldname(mm))//'SFSBD', sflxbcdp, pcols, lchnk)
                    end if
                  else
                      sflxec(1:ncol)   = 0.0_r8
                      sflxecdp(1:ncol) = 0.0_r8
                  end if !  deepconv_wetdep_history

                   ! when ma_convproc_intr is used, convective in-cloud wet removal is done there
                   ! the convective (total and deep) precip-evap-resuspension includes in- and below-cloud
                   ! contributions
                   ! so pass the below-cloud contribution to ma_convproc_intr
                   qsrflx_mzaer2cnvpr(1:ncol,mm,1) = sflxec(  1:ncol)
                   qsrflx_mzaer2cnvpr(1:ncol,mm,2) = sflxecdp(1:ncol)

                endif ! convproc_do_aer



             elseif (lphase == 2) then
                do_lphase2 = .true.

                if (do_lphase2) then

                   dqdt_tmp(:,:) = 0.0_r8

                  if (l <= nspec(m) ) then   ! species

                   if (convproc_do_aer) then
                      fldcw = qqcw(mm)%fld
                      qqcw_sav(1:ncol,:,l) = fldcw(1:ncol,:)
                   else
                      fldcw = qqcw(mm)%fld
                   endif

                   call wetdepa_v2(state%pmid, state%q(:,:,1), state%pdel, &
                        dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                        dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                        dep_inputs%evapr, dep_inputs%totcond, fldcw, dt, &
                        dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                        dlf, fracis_cw, sol_factb, ncol, &
                        scavcoefnv(:,:,jnv), &
                        is_strat_cloudborne=.true.,  &
                        icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                        convproc_do_aer=convproc_do_aer, rcscavt=rcscavt, rsscavt=rsscavt,  &
                        sol_facti_in=sol_facti, sol_factic_in=sol_factic, &
                        convproc_do_evaprain_atonce_in=convproc_do_evaprain_atonce, &
                        bergso_in=dep_inputs%bergso )

                   if(convproc_do_aer) then
                      ! save resuspension of cloudborne species
                      rtscavt(1:ncol,:,l) = rcscavt(1:ncol,:) + rsscavt(1:ncol,:)
                      ! wetdepa_v2 adds the resuspension of cloudborne to the dqdt of cloudborne (as a source)
                      ! undo this, so the resuspension of cloudborne can be added to the dqdt of interstitial (above)
                      dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) - rtscavt(1:ncol,:,l)
                   endif

                   ! sum up dqdt_tmp for all species (not total mmr and number)
                   dqdt_tmp_sum(1:ncol,:) = dqdt_tmp_sum(1:ncol,:) + dqdt_tmp(1:ncol,:)

                   qqcw(mm)%fld(1:ncol,:) = qqcw(mm)%fld(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

                  else if (l == nspec(m) + 1) then     !mmr and number are not advected

                   do k= 1,pver
                      do i= 1,ncol
                        if (qqcw(mm)%fld(i,k) .gt. 0._r8) then
                              frac_dqdt(i,k) = dqdt_tmp_sum(i,k) / qqcw(mm)%fld(i,k)
                        else
                              frac_dqdt(i,k) = 0._r8
                        end if
                      end do
                   end do

                   dqdt_tmp(1:ncol,:)  = dqdt_tmp_sum(1:ncol,:)
                   qqcw(mm)%fld(1:ncol,:) = qqcw(mm)%fld(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

                  else if (l == nspec(m) + 2) then    !num fraction is assuming dN/N = dM/M
                   dqdt_tmp(1:ncol,:) = frac_dqdt(1:ncol,:) * qqcw(mm)%fld(1:ncol,:)
                   qqcw(mm)%fld(1:ncol,:) = qqcw(mm)%fld(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

                  end if   ! nspec

                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                  if (l == nspec(m) + 1) then   !for mmr only
                   call outfld( trim(fieldname_cw(mm))//'SFWET', sflx, pcols, lchnk)
                  end if
                   aerdepwetcw(:ncol,mm) = sflx(:ncol)

                  if (l == nspec(m) + 1) then   !for mmr only
                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(fieldname_cw(mm))//'SFSIC', sflx, pcols, lchnk)
                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(fieldname_cw(mm))//'SFSIS', sflx, pcols, lchnk)
                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(fieldname_cw(mm))//'SFSBC', sflx, pcols, lchnk)
                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(fieldname_cw(mm))//'SFSBS', sflx, pcols, lchnk)

                   if(convproc_do_aer) then
                      sflx(:)=0.0_r8
                      do k=1,pver
                         sflx(1:ncol)=sflx(1:ncol)+rcscavt(1:ncol,k)*state%pdel(1:ncol,k)/gravit
                      enddo
                      call outfld( trim(fieldname_cw(mm))//'SFSEC', sflx, pcols, lchnk)

                      sflx(:)=0.0_r8
                      do k=1,pver
                         sflx(1:ncol)=sflx(1:ncol)+rsscavt(1:ncol,k)*state%pdel(1:ncol,k)/gravit
                      enddo
                      call outfld( trim(fieldname_cw(mm))//'SFSES', sflx, pcols, lchnk)
                   endif
                  end if ! nspec +1 history fields only
                endif
             endif

          enddo ! l= 0, nspec(m)+2
       enddo ! lphase = 1, 2
    enddo ! m = 1, nbins


    if (convproc_do_aer) then
       call t_startf('ma_convproc')
       call ma_convproc_intr( state, ptend, pbuf, dt,                &
            nsrflx_mzaer2cnvpr, qsrflx_mzaer2cnvpr, aerdepwetis, &
            dcondt_resusp3d)

       if (convproc_do_evaprain_atonce) then
        ! st Francis can we ADD END RUN statement here: "convproc_do_evaprain_atonce does not work with CARMA, needs to be set to .false."
       end if

       call t_stopf('ma_convproc')
    endif

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    !st if (.not. aerodep_flx_prescribed()) then
    !st    call set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)
    !st endif

  endsubroutine aero_model_wetdep

  !-------------------------------------------------------------------------
  ! provides aerosol surface area info for sectional aerosols
  ! called from mo_usrrxt
  !-------------------------------------------------------------------------
  subroutine aero_model_surfarea( &
                  mmr, radmean, relhum, pmid, temp, strato_sad, sulfate,  m, ltrop, &
                  dlat, het1_ndx, pbuf, ncol, sfc, dm_aer, sad_total, reff_trop )

    ! dummy args
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: radmean      ! mean radii in cm
    real(r8), intent(in)    :: strato_sad(:,:)
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: ltrop(:)
    real(r8), intent(in)    :: dlat(:)                    ! degrees latitude
    integer,  intent(in)    :: het1_ndx
    real(r8), intent(in)    :: relhum(:,:)
    real(r8), intent(in)    :: m(:,:) ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_total(:,:)
    real(r8), intent(out)   :: reff_trop(:,:)

    reff_trop(:,:)=0._r8

  end subroutine aero_model_surfarea

  !-------------------------------------------------------------------------
  ! stub
  !-------------------------------------------------------------------------
  subroutine aero_model_strat_surfarea( ncol, mmr, pmid, temp, ltrop, pbuf, strato_sad, reff_strat )

    ! dummy args
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:) ! tropopause level indices
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out)   :: strato_sad(:,:)
    real(r8), intent(out)   :: reff_strat(:,:)

    strato_sad(:,:) = 0._r8
    reff_strat(:,:) = 0._r8

  end subroutine aero_model_strat_surfarea

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( loffset, ncol, lchnk, troplev, delt, reaction_rates, &
                                    tfld, pmid, pdel, mbar, relhum, &
                                    zm,  qh2o, cwat, cldfr, cldnum, &
                                    airdens, invariants, del_h2so4_gasprod,  &
                                    vmr0, vmr, pbuf )

    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
    integer,  intent(in) :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in) :: ncol                   ! number columns in chunk
    integer,  intent(in) :: lchnk                  ! chunk index
    integer,  intent(in) :: troplev(:)
    real(r8), intent(in) :: delt                   ! time step size (sec)
    real(r8), intent(in) :: reaction_rates(:,:,:)  ! reaction rates
    real(r8), intent(in) :: tfld(:,:)              ! temperature (K)
    real(r8), intent(in) :: pmid(:,:)              ! pressure at model levels (Pa)
    real(r8), intent(in) :: pdel(:,:)              ! pressure thickness of levels (Pa)
    real(r8), intent(in) :: mbar(:,:)              ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: relhum(:,:)            ! relative humidity
    real(r8), intent(in) :: airdens(:,:)           ! total atms density (molec/cm**3)
    real(r8), intent(in) :: invariants(:,:,:)
    real(r8), intent(in) :: del_h2so4_gasprod(:,:)
    real(r8), intent(in) :: zm(:,:)
    real(r8), intent(in) :: qh2o(:,:)
    real(r8), intent(in) :: cwat(:,:)          ! cloud liquid water content (kg/kg)
    real(r8), intent(in) :: cldfr(:,:)
    real(r8), intent(in) :: cldnum(:,:)       ! droplet number concentration (#/kg)
    real(r8), intent(in) :: vmr0(:,:,:)       ! initial mixing ratios (before gas-phase chem changes)
    real(r8), intent(inout) :: vmr(:,:,:)         ! mixing ratios ( vmr )

    type(physics_buffer_desc), pointer :: pbuf(:)

  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

  end subroutine aero_model_emissions


  !===============================================================================


end module aero_model
