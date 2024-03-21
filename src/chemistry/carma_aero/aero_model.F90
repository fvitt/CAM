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
  use physics_buffer,    only: pbuf_get_field, pbuf_set_field, dtype_r8
  use physconst,         only: gravit, rair, rhoh2o
  use spmd_utils,        only: masterproc
  use cam_history,       only: outfld
  use chem_mods,         only: gas_pcnst, adv_mass
  use mo_tracname,       only: solsym
  use infnan,            only: nan, assignment(=)
  use rad_constituents,  only: rad_cnst_get_info, rad_cnst_get_info_by_bin, &
                               rad_cnst_get_info_by_bin_spec, rad_cnst_get_bin_props_by_idx, &
                               rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_mmr, &
                               rad_cnst_get_bin_num
  use mo_setsox,         only: setsox, has_sox

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
  integer :: pblh_idx            = 0
  integer :: wetdens_ap_idx      = 0

  integer :: fracis_idx          = 0
  integer :: prain_idx           = 0
  integer :: rprddp_idx          = 0
  integer :: rprdsh_idx          = 0
  integer :: nevapr_shcu_idx     = 0
  integer :: nevapr_dpcu_idx     = 0

  integer :: sulfeq_idx = -1

  integer :: nh3_ndx    = 0
  integer :: nh4_ndx    = 0
  integer :: h2so4_ndx  = 0

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

    use carma_flags_mod,   only: carma_model

    integer :: m, l, i
    integer :: nsoa_vbs
    character(len=32) :: spectype
    character(len=32) :: spec_name
    character(len=32) :: mmr_name
    character(len=32) :: soag_name
    character(len=32) :: soa_name

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

  ! SOA information
  ! Define number of VBS bins (nsoa) based on number of SOAG chemistry species
    nsoa_vbs = 0
    do i = 1, pcnst
       if (cnst_name(i)(:4) == 'SOAG') then
          nsoa_vbs = nsoa_vbs + 1
       end if
    end do
    if (masterproc) then
       write(iulog,*) 'nsoa_vbs  = ', nsoa_vbs
    endif

   ! Define pbuf field for soa_fraction
    call pbuf_add_field('FRACVBS','global',dtype_r8,(/pcols,pver,nbins,nsoa_vbs/), idx)

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
    use mo_setsox,       only: sox_inti

    !st use modal_aero_calcsize,   only: modal_aero_calcsize_init
    !st use modal_aero_coag,       only: modal_aero_coag_init
    !st use modal_aero_deposition, only: modal_aero_deposition_init
    use carma_aero_gasaerexch, only: carma_aero_gasaerexch_init
    !st use modal_aero_newnuc,     only: modal_aero_newnuc_init
    !st use modal_aero_rename,     only: modal_aero_rename_init

    use carma_aero_convproc,   only: ma_convproc_init
    use time_manager,    only: is_first_step
    use constituents,    only: cnst_set_convtran2

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

    character(len=32) :: spectype
    character(len=32) :: spec_name
    character(len=32) :: mmr_name
    character(len=32) :: bin_name

    integer :: idx
    real(r8) :: nanval

    if (is_first_step()) then
       do m = 1, nbins
          call rad_cnst_get_info_by_bin(0, m, mmr_name=mmr_name)
          idx = pbuf_get_index('CLD'//trim(mmr_name))
          call pbuf_set_field(pbuf2d, idx, 0.0_r8)
          idx = pbuf_get_index('CLDNB'//trim(mmr_name))
          call pbuf_set_field(pbuf2d, idx, 0.0_r8)
          do l = 1, nspec(m)
             call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=spec_name)
             idx = pbuf_get_index('CLD'//trim(spec_name))
             call pbuf_set_field(pbuf2d, idx, 0.0_r8)
          enddo
       enddo
    endif

    ! define pbuf field for soa_fraction
    if (is_first_step()) then
       nanval = nan
       idx = pbuf_get_index('FRACVBS')
       call pbuf_set_field(pbuf2d, idx, nanval)
    end if

    ! aqueous chem initialization
    call sox_inti()

    h2so4_ndx = get_spc_ndx('H2SO4')
    nh3_ndx = get_spc_ndx('NH3')
    nh4_ndx = get_spc_ndx('NH4')



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

    call carma_aero_bcscavcoef_init(pbuf2d)

    !st  call modal_aero_rename_init( modal_accum_coarse_exch )
    !   calcsize call must follow rename call
    !st call modal_aero_calcsize_init( pbuf2d )
    call carma_aero_gasaerexch_init
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
             call cnst_set_convtran2(idxtmp, .not.convproc_do_aer)
          else
             bin_cnst_lq(m,l) = .false.
             bin_cnst_idx(m,l) = 0
          end if

         mm = ii

         unit_basename = 'kg'
         if (l == nspec(m) + 2) then   ! number
          unit_basename = ' 1'
         end if

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

       call addfld (trim(fieldname(mm))//'WETC',  (/ 'lev' /), 'A',unit_basename//'/kg/s ','wet deposition tendency??')
       call addfld (trim(fieldname(mm))//'CONU',  (/ 'lev' /), 'A',unit_basename//'/kg ','updraft mixing ratio??')
       call addfld (trim(fieldname(mm))//'QCONST',(/ 'lev' /), 'A',unit_basename//'/kg ','all mixing ratio??')

       call addfld (trim(fieldname_cw(mm))//'WETC',(/ 'lev' /), 'A',unit_basename//'/kg/s ','wet deposition tendency??')
       call addfld (trim(fieldname_cw(mm))//'CONU',(/ 'lev' /), 'A',unit_basename//'/kg ','updraft mixing ratio??')

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
          call addfld (trim(fieldname_cw(mm))//'WET',(/ 'lev' /), 'A',unit_basename//'/kg/s ','wet deposition tendency')
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
             call addfld( trim(fieldname_cw(mm))//'RSPTD', (/ 'lev' /), 'A', unit_basename//'/kg/s',   &
                trim(fieldname_cw(mm))//' resuspension tendency')
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
                call add_default (trim(fieldname_cw(mm))//'RSPTD', 1, ' ')
             end if
          endif
       enddo
    enddo

    do m = 1,gas_pcnst

       unit_basename = 'kg'  ! Units 'kg' or '1'

       call addfld( 'GS_'//trim(solsym(m)),horiz_only, 'A', unit_basename//'/m2/s ', &
                    trim(solsym(m))//' gas chemistry/wet removal (for gas species)')
       call addfld( 'AQ_'//trim(solsym(m)),horiz_only, 'A', unit_basename//'/m2/s ', &
                    trim(solsym(m))//' aqueous chemistry (for gas species)')
       if ( history_aerosol ) then
          call add_default( 'AQ_'//trim(solsym(m)), 1, ' ')
       endif

    enddo

    if (has_sox) then
       do n = 1, nbins
          do l = 1, nspec(n)   ! not for total mass or number
             mm = bin_idx(n, l)
             call addfld (&
                  trim(fieldname_cw(mm))//'AQSO4',horiz_only,  'A','kg/m2/s', &
                  trim(fieldname_cw(mm))//' aqueous phase chemistry')
             call addfld (&
                  trim(fieldname_cw(mm))//'AQH2SO4',horiz_only,  'A','kg/m2/s', &
                  trim(fieldname_cw(mm))//' aqueous phase chemistry')
             if ( history_aerosol ) then
                call add_default (trim(fieldname_cw(mm))//'AQSO4', 1, ' ')
                call add_default (trim(fieldname_cw(mm))//'AQH2SO4', 1, ' ')
             endif
          end do
       end do

       call addfld( 'XPH_LWC',    (/ 'lev' /), 'A','kg/kg',   'pH value multiplied by lwc')
       call addfld ('AQSO4_H2O2', horiz_only,  'A','kg/m2/s', 'SO4 aqueous phase chemistry due to H2O2')
       call addfld ('AQSO4_O3',   horiz_only,  'A','kg/m2/s', 'SO4 aqueous phase chemistry due to O3')

       if ( history_aerosol ) then
          call add_default ('XPH_LWC', 1, ' ')
          call add_default ('AQSO4_H2O2', 1, ' ')
          call add_default ('AQSO4_O3', 1, ' ')
       endif
    endif

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
    use time_manager,          only: is_first_step


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

    integer :: mm, lc, lpr
    integer :: i,k


    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)
    real(r8) :: totalmmr(pcols, pver)
    !st real(r8) :: sol_factb, sol_facti
    real(r8) :: sol_factb(pcols,pver)
    real(r8) :: sol_facti(pcols,pver)
    real(r8) :: sol_factic(pcols,pver)
    real(r8) :: sol_factb_c(pcols,pver)
    real(r8) :: sol_facti_c(pcols,pver)
    real(r8) :: sol_factic_c(pcols,pver)

    real(r8) :: sflx(pcols) ! deposition flux

    integer :: jnv ! index for scavcoefnv 3rd dimension
    integer :: lphase ! index for interstitial / cloudborne aerosol
    integer :: strt_loop, end_loop, stride_loop !loop indices for the lphase loop
    integer :: l ! index for aerosol number / chem-mass / water-mass
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species
    real(r8) :: icscavt_sum(pcols,pver) ! temporary array to hold sum of tendency for all species
    real(r8) :: isscavt_sum(pcols,pver) ! temporary array to hold sum of tendency for all species
    real(r8) :: bcscavt_sum(pcols,pver) ! temporary array to hold sum of tendency for all species
    real(r8) :: bsscavt_sum(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: rcscavt_sum(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: rsscavt_sum(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: f_act_conv(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: f_act_conv_c(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: fracis_cw(pcols,pver)
    real(r8) :: fracis_nadv(pcols,pver)
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
    real(r8) :: specmmr
    real(r8) :: solmmr
    real(r8) :: so4mmr
    real(r8) :: rho_water
    real(r8) :: specdens              ! specie density from physprops files

    real(r8) :: sflxec(pcols), sflxecdp(pcols)  ! deposition flux
    real(r8) :: sflxic(pcols), sflxicdp(pcols)  ! deposition flux
    real(r8) :: sflxbc(pcols), sflxbcdp(pcols)  ! deposition flux
    real(r8) :: rcscavt(pcols, pver)
    real(r8) :: rsscavt(pcols, pver)
    real(r8) :: qqcw_in(pcols,pver), qqcw_sav(pcols,pver,0:nspec_max) ! temporary array to hold qqcw for the current mode
    real(r8) :: fldcw(pcols,pver)
    real(r8) :: fldcw_sum(pcols,pver)
    real(r8) :: fld_sum(pcols,pver)
    real(r8) :: rtscavt(pcols, pver, 0:nspec_max)

    integer, parameter :: nsrflx_mzaer2cnvpr = 2
    real(r8) :: qsrflx_mzaer2cnvpr(pcols,ncnst_tot,nsrflx_mzaer2cnvpr)
    ! End unified convection scheme

    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble

    real(r8), pointer :: dryr(:,:)   ! CARMA dry radius in cm
    real(r8), pointer :: wetr(:,:)   ! CARMA wet radius in cm
    real(r8), pointer :: rmass_ptr(:)     ! CARMA rmass
    real(r8), allocatable :: rmass(:)     ! CARMA rmass

    type(wetdep_inputs_t) :: dep_inputs

    real(r8) :: dcondt_resusp3d(ncnst_extd,pcols, pver)

    character(len=32) :: spectype
    character(len=32) :: bin_name

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

    !NOTE: currently all aerosols in CARMA are wet deposited, namelist setting or master table of ndep could be included later
    !st if (nwetdep<1) return

    call wetdep_inputs_set( state, pbuf, dep_inputs )

    allocate( &
      raer(ncnst_tot),                &
      qqcw(ncnst_tot)                 )

    call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

    prec(:ncol)=0._r8
    do k=1,pver
       where (prec(:ncol) >= 1.e-7_r8)
          isprx(:ncol,k) = .true.
       elsewhere
          isprx(:ncol,k) = .false.
       endwhere
       prec(:ncol) = prec(:ncol) + (dep_inputs%prain(:ncol,k) + dep_inputs%cmfdqr(:ncol,k) - dep_inputs%evapr(:ncol,k)) &
                   * state%pdel(:ncol,k)/gravit
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

    allocate ( rmass(nbins) )
    do m = 1, nbins      ! main loop over aerosol bins
       ! r(m) is the dry bin radius
       ! taken here from CARMA pbuf field
       ! get bin info general for each bin
       call rad_cnst_get_info_by_bin(0, m, bin_name=bin_name)
       call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_dryr"),dryr)
       call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_wetr"),wetr)
       call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_rmass"),rmass_ptr)
       rmass(m) =  rmass_ptr(1)

       do lphase = strt_loop,end_loop, stride_loop ! loop over interstitial (1) and cloud-borne (2) forms
          ! Init pointers to mode number and specie mass mixing ratios in
          ! intersitial and cloud borne phases.
          !st note only raer is needed here to calcuted  solfact_b

          if (lphase == 2) then !  cloud-borne aerosol (borne by stratiform cloud drops)
             sol_factb_c  = 0.0_r8   ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
             sol_facti_c  = sol_facti_cloud_borne   ! strat  in-cloud scav cloud-borne tuning factor  ! 1.0 in Pengfei's model
             sol_factic_c = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
             f_act_conv_c = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
          end if

          fldcw_sum(:,:) = 0.0_r8
          fld_sum(:,:) = 0.0_r8
          icscavt_sum(:,:) = 0.0_r8
          isscavt_sum(:,:) = 0.0_r8
          bcscavt_sum(:,:) = 0.0_r8
          bsscavt_sum(:,:) = 0.0_r8
          rcscavt_sum(:,:) = 0.0_r8
          rsscavt_sum(:,:) = 0.0_r8

          if (lphase == 1) then ! interstial aerosol

             sol_facti(:,:) = 0._r8
             sol_factic = sol_factic_interstitial
             f_act_conv = 0.8_r8 ! rce 2010/05/02
             if (convproc_do_aer) then
                ! for (stratiform)-cloudborne aerosols, convective wet removal
                !    (all forms) is zero, so no action is needed
                sol_factic = 0.0_r8
             end if

             call carma_aero_bcscavcoef_get( m, ncol, isprx, wetr, dryr,&
                  scavcoefnv(:,:,1), scavcoefnv(:,:,2), pbuf )

          !alternative function to drive  scavcoefnv for CARMA  using Pengfei's method
          !load mmr and number


          !if  (dryr(i,k) .le. 0._r8)  then
                  !          scavcoefnv(i,k,1) =  0.0_r8
                  !          scavcoefnv(i,k,2) =  0.0_r8
                  !else


                   !follow MAMs approach for calculating the scavenging approach
                  ! call carma_aero_bcscavcoef_get( m, l, ncol, isprx, wetr, dryr,&
                  !      scavcoefnv(i,k,1), scavcoefnv(i,k,2) )
                  !end if

                   ! This is calculation for scavcoefnv has been introducted by Pengfei Yu for
                   ! CARMA and is required in particular for a sectional model
                   ! The scavenging coefficients are calcuated as a function of the aerosol bin at each
                   ! grid point (lat,lon,vertical). Refer to Dana and Hales (1975), however I use
                   ! Marshall and Palmer Distribution instead log-normal for rain drop size distribution
                   ! Author: Pengfei Yu
                   ! Aug.25, 2011



                   !st z_scavcoef = 0.1_r8
                   !scavcoefnv(:,:,1) =  0.1_r8
                   !scavcoefnv(:,:,2) =  0.1_r8

                   ! convert precipitation prain (kg/kg/s) to mm/h

                   !J(:ncol,:) = ( dep_inputs%prain(:ncol,:) + dep_inputs%cmfdqr(:ncol,:) &
                   !             - dep_inputs%evapc(:ncol,:) - dep_inputs%evapr(:ncol,:) ) &
                   !           * state%pdel(:ncol,:)/gravit/rho_water*3.6e6_r8   ! m/s
                   ! Note needs to be revisited, here scavenging function is dependent on each specie:

                   !do l = 1, nspec(m)   ! loop over different npsec plus over total mmr and number
                   ! call rad_cnst_get_bin_props_by_idx(0, m, l, density_aer=specdens)

                   !do k= 1,pver
                   !   do i= 1,ncol
                   !      stoke(i,k) = 0.1_r8*(dryr(i,k)**2)*1.0e8_r8*specdens                    ! stoke parameter
                   !      if ((J(i,k) .le. 0._r8)  .or. (dryr(i,k) .le. 0._r8))  then
                   !         scavcoefnv(i,k,1) =  0.0_r8
                   !         scavcoefnv(i,k,2) =  0.0_r8
                   !      else
                   !         if (dryr(i,k) .lt. 1.e-4_r8) then
                   !            !z_scavcoef(i_col,k_ver) = 63.03750_r8*r(ibin)*(J(i_col,k_ver)**(-0.42_r8))
                   !            ! scavenging volume
                   !            scavcoefnv(i,k,2) = 252.15_r8*dryr(i,k)*(J(i,k)**(-0.42_r8))
                   !            ! scavenging number
                   !            !st z_scavcoef(i,k,0) = 252.15_r8*r(ibin)*(J(i_col,k_ver)**(-0.42_r8))
                   !         else
                   !            !z_scavcoef(i_col,k_ver) = 1.025_r8*(J(i_col,k_ver)**(-0.21_r8))*(((stoke - 0.0833_r8)/(stoke + 0.5833_r8))**1.5_r8)
                   !            ! scavenging volume
                   !            scavcoefnv(i,k,2) = 2.05_r8*(J(i,k)**(-0.21_r8)) &
                   !                              * (((stoke(i,k) - 0.0833_r8)/(stoke(i,k) + 0.5833_r8))**1.5_r8)
                   !            ! scavenging  number
                   !            !st z_scavcoef(i,k,1) = 2.05_r8*(J(i_col,k_ver)**(-0.21_r8))*(((stoke - 0.0833_r8)/(stoke + 0.5833_r8))**1.5_r8)
                   !         endif
                   !      endif
                   !   enddo
                   !enddo


             ! derive sol_factb for CARMA  specific for bins and phase
             sol_factb(:ncol,:)=0._r8
             ! get total mass (per bin)
             l = nspec(m) + 1
             mm = bin_idx(m, l)
             call rad_cnst_get_bin_mmr(0, m, 'a', state, pbuf, raer(mm)%fld)
             totalmmr = raer(mm)%fld
             do l = 1, nspec(m)
                mm = bin_idx(m, l)
                call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'a', state, pbuf, raer(mm)%fld)
             end do

             do k= 1,pver
                do i= 1,ncol
                   specmmr  = 0._r8
                   solmmr = 0._r8
                   so4mmr = 0._r8

                   ! loop through species:
                   do l = 1, nspec(m)
                      mm = bin_idx(m, l)
                      call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)

                      if (trim(spectype) == 'sulfate') then
                         solmmr =  solmmr + raer(mm)%fld(i,k)*0.5_r8
                      end if
                      if (trim(spectype) == 'black-c') then
                         specmmr = specmmr + raer(mm)%fld(i,k)
                      end if
                      if (trim(spectype) == 'p-organic') then
                         solmmr =  solmmr + raer(mm)%fld(i,k)*0.2_r8
                         specmmr = specmmr + raer(mm)%fld(i,k)
                      end if
                      if (trim(spectype) == 's-organic') then
                         solmmr =  solmmr + raer(mm)%fld(i,k)*0.2_r8
                         specmmr = specmmr + raer(mm)%fld(i,k)
                      end if
                      if (trim(spectype) == 'dust') then
                         solmmr =  solmmr + raer(mm)%fld(i,k)*0.1_r8
                         specmmr = specmmr + raer(mm)%fld(i,k)
                      end if
                      if (trim(spectype) == 'seasalt') then
                         solmmr =  solmmr + raer(mm)%fld(i,k)*0.8_r8
                         specmmr = specmmr + raer(mm)%fld(i,k)
                      end if

                   end do   !nspec
                   if (totalmmr(i,k) .gt. 0._r8) then
                      sol_factb(i,k) = solmmr/totalmmr(i,k)
                   end if
                   sol_factb(i,k) = min(max(sol_factb(i,k), 0.1_r8), 0.8_r8)
                end do ! ncol
             end do !pver
          end if !phase 1
          ! for scavenging coefficient, for CARMA this is specie independent (different than for MAM)
          if (lphase == 1) then
             jnv = 2
          else
             jnv = 0
          end if

          !load mmr and number   specie loop starts
          do l = 1, nspec(m) + 2  ! loop over different npsec plus over total mmr and number
             mm = bin_idx(m,l)
             lpr = bin_cnst_idx(m,l)
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


             if (lphase == 1) then
                !ptend%lq(lpr) = .TRUE.
                if (l <= nspec(m)) then
                   dqdt_tmp(:,:) = 0.0_r8
                   q_tmp(:ncol,:) = raer(mm)%fld(:ncol,:)

                   if(convproc_do_aer) then
                      !Feed in the saved cloudborne mixing ratios from phase 2
                      qqcw_in(:ncol,:) = qqcw_sav(:ncol,:,l)
                   else
                      fldcw = qqcw(mm)%fld
                      qqcw_in(:ncol,:) = fldcw(:ncol,:)
                   endif

                   if (bin_cnst_lq(m,l)) then
                      ! advected species (differences between advected and non-advected is fracis
                      ! do not calculate for number
                      call wetdepa_v2( state%pmid, state%q(:,:,1), state%pdel, &
                           dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                           dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                           dep_inputs%evapr, dep_inputs%totcond, q_tmp, dt, &
                           dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                           dlf, fracis(:,:,lpr), sol_factb, ncol, &
                           scavcoefnv(:,:,jnv), &
                           is_strat_cloudborne=.false.,  &
                           qqcw=qqcw_in(:,:),  &
                           f_act_conv=f_act_conv, &
                           icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                           convproc_do_aer=convproc_do_aer, rcscavt=rcscavt, rsscavt=rsscavt,  &
                           sol_facti_in=sol_facti, sol_factic_in=sol_factic, &
                           convproc_do_evaprain_atonce_in=convproc_do_evaprain_atonce )
                   else
                      call wetdepa_v2( state%pmid, state%q(:,:,1), state%pdel, &
                           dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                           dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                           dep_inputs%evapr, dep_inputs%totcond, q_tmp, dt, &
                           dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                           dlf, fracis_nadv, sol_factb, ncol, &
                           scavcoefnv(:,:,jnv), &
                           is_strat_cloudborne=.false.,  &
                           qqcw=qqcw_in(:,:),  &
                           f_act_conv=f_act_conv, &
                           icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                           convproc_do_aer=convproc_do_aer, rcscavt=rcscavt, rsscavt=rsscavt,  &
                           sol_facti_in=sol_facti, sol_factic_in=sol_factic, &
                           convproc_do_evaprain_atonce_in=convproc_do_evaprain_atonce )
                   end if


                   if(convproc_do_aer) then
                      ! add resuspension of cloudborne species to dqdt of interstitial species
                      dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) + rtscavt(1:ncol,:,l)
                   endif

                   ! make sure raer values are not negative, adjust dqdt_tmp accordingly
                   do k= 1,pver
                      do i= 1,ncol

                         ! some CARMA interstitial species are not advected
                         if (bin_cnst_lq(m,l)) then
                            ! advected species, updated with ptend
                            !
                            ! NOTE: Need to take into account the tendency accumulated
                            ! so far as raer is the original state.
                            if ((raer(mm)%fld(i,k) + (ptend%q(i,k,lpr) + dqdt_tmp(i,k)) * dt) .lt. 0.0_r8) then
                               dqdt_tmp(i,k) = - raer(mm)%fld(i,k) / dt - ptend%q(i,k,lpr)
                            end if
                        else
                            ! pbuf species, updated directly
                            if (((raer(mm)%fld(i,k) + dqdt_tmp(i,k)) * dt) .lt. 0.0_r8) then
                               dqdt_tmp(i,k) = - raer(mm)%fld(i,k) / dt
                            end if
                        end if
                      end do
                   end do

                   ! some CARMA interstetial species are not advected, check:
                   if (bin_cnst_lq(m,l)) then
                      ! advected species, updated with ptend
                      ptend%q(1:ncol,:,lpr) = ptend%q(1:ncol,:,lpr) + dqdt_tmp(1:ncol,:)
                      fld_sum(1:ncol,:)  = fld_sum(1:ncol,:) + raer(mm)%fld(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt
                   else
                      ! pbuf species, updated directly
                      raer(mm)%fld(1:ncol,:) = max(raer(mm)%fld(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt, 0.0_r8)
                      fld_sum(1:ncol,:) = fld_sum(1:ncol,:) + raer(mm)%fld(1:ncol,:)
                   end if

                ! NOTE: For a flexible interface, you should not assume that mmr is advected
                ! and num is not advected. Removing the concentration element will cause both
                ! of these to be non-advected.
                else if (l == nspec(m) + 1) then   !mmr is advected and number is  not advected

                   if (bin_cnst_lq(m,l)) then ! adveced species
                      dqdt_tmp(1:ncol,:) = (fld_sum(1:ncol,:) - raer(mm)%fld(1:ncol,:)) / dt
                      ptend%q(1:ncol,:,lpr) = dqdt_tmp(1:ncol,:)
                   else
                      dqdt_tmp(1:ncol,:) = (raer(mm)%fld(1:ncol,:) - fld_sum(1:ncol,:)) / dt
                      raer(mm)%fld(1:ncol,:) = fld_sum(1:ncol,:)
                   end if

                else if (l == nspec(m) + 2) then   !num fraction is assuming dN/N = dM/M

                   if (bin_cnst_lq(m,l)) then ! adveced species
                      dqdt_tmp(1:ncol,:) = (fld_sum(1:ncol,:) / rmass(m) - raer(mm)%fld(1:ncol,:)) / dt
                      ptend%q(1:ncol,:,lpr) = dqdt_tmp(1:ncol,:)
                   else
                      dqdt_tmp(1:ncol,:) = (fld_sum(1:ncol,:) / rmass(m) - raer(mm)%fld(1:ncol,:)) / dt
                      raer(mm)%fld(1:ncol,:) = fld_sum(1:ncol,:) / rmass(m)
                   end if
                end if   ! nspec

                call outfld( trim(fieldname(mm))//'SIS', isscavt, pcols, lchnk)
                call outfld( trim(fieldname(mm))//'SBC', bcscavt, pcols, lchnk)
                call outfld( trim(fieldname(mm))//'SBS', bsscavt, pcols, lchnk)
                call outfld( trim(fieldname(mm))//'SIC', icscavt, pcols, lchnk)
                call outfld( trim(fieldname(mm))//'WET', dqdt_tmp(:,:), pcols, lchnk)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                if (.not.convproc_do_aer) call outfld( trim(fieldname(mm))//'SFWET', sflx, pcols, lchnk)
                aerdepwetis(:ncol,mm) = sflx(:ncol)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                if (.not.convproc_do_aer) call outfld( trim(fieldname(mm))//'SFSIC', sflx, pcols, lchnk)

                if (convproc_do_aer) sflxic = sflx

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(fieldname(mm))//'SFSIS', sflx, pcols, lchnk)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(fieldname(mm))//'SFSBC', sflx, pcols, lchnk)
                if (convproc_do_aer)sflxbc = sflx

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(fieldname(mm))//'SFSBS', sflx, pcols, lchnk)

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
                   call outfld( trim(fieldname(mm))//'SFSES', sflx, pcols, lchnk)


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
                      call outfld( trim(fieldname(mm))//'SFSBD', sflxbcdp, pcols, lchnk)
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
                ! for all cases in CARMA
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
                           dlf, fracis_cw, sol_factb_c, ncol, &
                           scavcoefnv(:,:,jnv), &
                           is_strat_cloudborne=.true.,  &
                           icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                           convproc_do_aer=convproc_do_aer, rcscavt=rcscavt, rsscavt=rsscavt,  &
                           sol_facti_in=sol_facti_c, sol_factic_in=sol_factic_c, &
                           convproc_do_evaprain_atonce_in=convproc_do_evaprain_atonce, &
                           bergso_in=dep_inputs%bergso )

                      if(convproc_do_aer) then
                         ! save resuspension of cloudborne species
                         rtscavt(1:ncol,:,l) = rcscavt(1:ncol,:) + rsscavt(1:ncol,:)
                         ! wetdepa_v2 adds the resuspension of cloudborne to the dqdt of cloudborne (as a source)
                         ! undo this, so the resuspension of cloudborne can be added to the dqdt of interstitial (above)
                         dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) - rtscavt(1:ncol,:,l)
                      endif

                      do k= 1,pver
                         do i= 1,ncol
                            if ( (qqcw(mm)%fld(i,k) + dqdt_tmp(i,k) * dt) .lt. 0.0_r8 )   then
                               dqdt_tmp(i,k) = - qqcw(mm)%fld(i,k) / dt
                            end if
                         end do
                      end do


                      ! sum up dqdt_tmp for all species (not total mmr and number)
                      do k= 1,pver
                         do i= 1,ncol
                            qqcw(mm)%fld(i,k) = max(qqcw(mm)%fld(i,k) + dqdt_tmp(i,k) * dt, 0.0_r8)
                            fldcw_sum(i,k) = fldcw_sum(i,k) + qqcw(mm)%fld(i,k)
                         end do
                      end do
                      icscavt_sum(:ncol,:) = icscavt_sum(:ncol,:) + icscavt(:ncol,:)
                      isscavt_sum(:ncol,:) = isscavt_sum(:ncol,:) + isscavt(:ncol,:)
                      bcscavt_sum(:ncol,:) = bcscavt_sum(:ncol,:) + bcscavt(:ncol,:)
                      bsscavt_sum(:ncol,:) = bsscavt_sum(:ncol,:) + bsscavt(:ncol,:)
                      rcscavt_sum(:ncol,:) = rcscavt_sum(:ncol,:) + rcscavt(:ncol,:)
                      rsscavt_sum(:ncol,:) = rsscavt_sum(:ncol,:) + rsscavt(:ncol,:)

                   else if (l == nspec(m) + 1) then     !mmr and number are not advected

                      dqdt_tmp(:ncol,:) = (fldcw_sum(:ncol,:) - qqcw(mm)%fld(:ncol,:)) / dt
                      qqcw(mm)%fld(:ncol,:) = fldcw_sum(:ncol,:)

                      icscavt(:ncol,:) = icscavt_sum(:ncol,:)
                      isscavt(:ncol,:) = isscavt_sum(:ncol,:)
                      bcscavt(:ncol,:) = bcscavt_sum(:ncol,:)
                      bsscavt(:ncol,:) = bsscavt_sum(:ncol,:)
                      rcscavt(:ncol,:) = rcscavt_sum(:ncol,:)
                      rsscavt(:ncol,:) = rsscavt_sum(:ncol,:)

                   else if (l == nspec(m) + 2) then    !num fraction is assuming dN/N = dM/M

                      dqdt_tmp(:ncol,:) = (fldcw_sum(:ncol,:) / rmass(m) - qqcw(mm)%fld(:ncol,:)) / dt
                      qqcw(mm)%fld(:ncol,:) = fldcw_sum(:ncol,:) / rmass(m)

                      icscavt(:ncol,:) = icscavt_sum(:ncol,:)
                      isscavt(:ncol,:) = isscavt_sum(:ncol,:)
                      bcscavt(:ncol,:) = bcscavt_sum(:ncol,:)
                      bsscavt(:ncol,:) = bsscavt_sum(:ncol,:)
                      rcscavt(:ncol,:) = rcscavt_sum(:ncol,:)
                      rsscavt(:ncol,:) = rsscavt_sum(:ncol,:)

                   end if   ! nspec

                   sflx(:)=0._r8
                   do k=1,pver
                      do i=1,ncol
                         sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                      enddo
                   enddo
                   call outfld( trim(fieldname_cw(mm))//'SFWET', sflx, pcols, lchnk)
                   aerdepwetcw(:ncol,mm) = sflx(:ncol)

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
          do m = 1, nbins ! loop over aerosol bins
             ! apply change to cloud-borne (lphase=2) aerosols
             do l = 1, nspec(m) + 2 ! loop over different npsec plus over total mmr and number
                mm = bin_idx(m, l)
                lc = mm + ncnst_tot
                qqcw(mm)%fld(:ncol,:) = qqcw(mm)%fld(:ncol,:) + dcondt_resusp3d(lc,:ncol,:)*dt
                call outfld( trim(fieldname_cw(mm))//'RSPTD', dcondt_resusp3d(lc,:ncol,:), ncol, lchnk )
             end do ! loop over number + mmr +  chem constituents
          end do ! m aerosol modes
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
  ! provides wet tropospheric aerosol surface area info for sectional aerosols
  ! called from mo_usrrxt
  !-------------------------------------------------------------------------
  subroutine aero_model_surfarea( &
                  state, mmr, radmean, relhum, pmid, temp, strato_sad, sulfate,  m, ltrop, &
                  dlat, het1_ndx, pbuf, ncol, sfc, dm_aer, sad_trop, reff_trop )

    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
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
    real(r8), intent(inout) :: sad_trop(:,:)  ! aerosol surface area density (cm2/cm3), zeroed above the tropopause
    real(r8), intent(out)   :: reff_trop(:,:) ! aerosol effective radius (cm), zeroed above the tropopause

    ! local vars
    integer :: beglev(ncol)
    integer :: endlev(ncol)

    beglev(:ncol)=ltrop(:ncol)+1
    endlev(:ncol)=pver
    call surf_area_dens( state, pbuf, ncol, mmr, beglev, endlev, sad_trop, reff_trop, sfc=sfc, dm_aer=dm_aer )

  end subroutine aero_model_surfarea

  !-------------------------------------------------------------------------
  ! provides wet stratospheric aerosol surface area info for sectional aerosols
  ! called from mo_gas_phase_chemdr.F90
  !-------------------------------------------------------------------------
  subroutine aero_model_strat_surfarea( state, ncol, mmr, pmid, temp, ltrop, pbuf, strato_sad, reff_strat )

    use ref_pres, only: clim_modal_aero_top_lev

    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:) ! tropopause level indices
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out)   :: strato_sad(:,:) ! aerosol surface area density (cm2/cm3), zeroed below the tropopause
    real(r8), intent(out)   :: reff_strat(:,:) ! aerosol effective radius (cm), zeroed below the tropopause

    ! local vars
    integer :: beglev(ncol)
    integer :: endlev(ncol)

    beglev(:ncol) = clim_modal_aero_top_lev
    endlev(:ncol) = ltrop(:ncol)

    call surf_area_dens( state, pbuf, ncol, mmr, beglev, endlev, strato_sad, reff_strat )

  end subroutine aero_model_strat_surfarea

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( state, loffset, ncol, lchnk, troplev, delt, reaction_rates, &
                                    tfld, pmid, pdel, mbar, relhum, &
                                    zm,  qh2o, cwat, cldfr, cldnum, &
                                    airdens, invariants, del_h2so4_gasprod,  &
                                    vmr0, vmr, pbuf )

    use carma_aero_gasaerexch, only : carma_aero_gasaerexch_sub
    use time_manager,          only : get_nstep
    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
    type(physics_state), intent(in)    :: state    ! Physics state variables
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

    ! local vars

    integer :: n, m, mm
    integer :: i,k,l
    integer :: nstep

    type(ptr2d_t), allocatable :: raer(:)     ! aerosol mass, number mixing ratios
    type(ptr2d_t), allocatable :: qqcw(:)

    real(r8) :: del_h2so4_aeruptk(ncol,pver)

    !st real(r8), pointer :: dgnum(:,:,:), dgnumwet(:,:,:), wetdens(:,:,:)
    real(r8), pointer :: pblh(:)                    ! pbl height (m)

    real(r8), dimension(ncol) :: wrk
    character(len=32)         :: name
    real(r8) :: dvmrcwdt(ncol,pver,ncnst_tot)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: delta_so4mass(ncol,pver,ncnst_tot)
    real(r8) :: wetr_n(pcols,pver,nbins)       ! wet radius from CARMA for different bin
    !st real(r8) :: vmrcw(ncol,pver,gas_pcnst)            ! cloud-borne aerosol (vmr)
    !st vmrcw is going only through CARMA aerosols (ncnst_tot)
    real(r8) :: vmrcw(ncol,pver,ncnst_tot)            ! cloud-borne aerosol (vmr)
    real(r8) :: mmrcw(ncol,pver,ncnst_tot)            ! cloud-borne aerosol (mmr)
    real(r8) :: raervmr(ncol,pver,ncnst_tot)            ! cloud-borne aerosol (vmr)

    real(r8) ::  aqso4(ncol,ncnst_tot)               ! aqueous phase chemistry
    real(r8) ::  aqh2so4(ncol,ncnst_tot)             ! aqueous phase chemistry
    real(r8) ::  aqso4_h2o2(ncol)                     ! SO4 aqueous phase chemistry due to H2O2
    real(r8) ::  aqso4_o3(ncol)                       ! SO4 aqueous phase chemistry due to O3
    real(r8) ::  xphlwc(ncol,pver)                    ! pH value multiplied by lwc
    real(r8) ::  nh3_beg(ncol,pver)
    real(r8) ::  mw_carma(ncnst_tot)
    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: sulfeq(:,:,:)
    real(r8), pointer :: dryr(:,:)   ! CARMA dry radius in cm
    real(r8), pointer :: wetr(:,:)   ! CARMA wet radius in cm
    real(r8), pointer :: rmass_ptr(:)     ! CARMA rmass
    real(r8), allocatable :: rmass(:)     ! CARMA rmass

    real(r8) :: old_total_mass
    real(r8) :: new_total_mass
    real(r8) :: old_total_number


    character(len=32) :: bin_name
    character(len=32) :: spectype
    logical :: is_spcam_m2005
!
! ... initialize nh3
!
    if ( nh3_ndx > 0 ) then
      nh3_beg = vmr(1:ncol,:,nh3_ndx)
    end if
!
    is_spcam_m2005   = cam_physpkg_is('spcam_m2005')

    !st call pbuf_get_field(pbuf, dgnum_idx,      dgnum)
    !st call pbuf_get_field(pbuf, dgnumwet_idx,   dgnumwet )
    !st call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens )
    !st call pbuf_get_field(pbuf, pblh_idx,       pblh)

  ! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    nstep = get_nstep()

    ! calculate tendency due to gas phase chemistry and processes
    dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
    do m = 1, gas_pcnst
      wrk(:) = 0.0_r8
      do k = 1,pver
        wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m)*adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
      end do
      name = 'GS_'//trim(solsym(m))
      call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

!
! Aerosol processes ...
!
    allocate( &
      rmass(nbins),                &
      raer(ncnst_tot),                &
      qqcw(ncnst_tot)                 )

    mw_carma(:) = 0.0_r8
    do m = 1, nbins      ! main loop over aerosol bins
       !st can we move this part to init???
       ! dryr is the dry bin radius
       ! wetr is the dry bin radius
       ! Note: taken here from CARMA pbuf field which may be not any more consistent with changed fields after carma was applied
       ! Need to add new code that recalcuates dryr and wetr
       ! get bin info
       call rad_cnst_get_info_by_bin(0, m, nspec=nspec(m), bin_name=bin_name)
       call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_wetr"),wetr)
       call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_rmass"),rmass_ptr)
       rmass(m) =  rmass_ptr(1)

       wetr_n(:,:,m) = wetr(:,:)

       ! Init pointers to mode number and specie mass mixing ratios in
       ! intersitial and cloud borne phases.
       do l = 1, nspec(m) + 2
          mm = bin_idx(m, l)
          if (l <= nspec(m)) then
             call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
             call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'a', state, pbuf, raer(mm)%fld)
             call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
             if (trim(spectype) == 'sulfate') then
                mw_carma(mm) = 96._r8
             end if
             if (trim(spectype) == 'black-c') then
                mw_carma(mm) = 12._r8
             end if
             if (trim(spectype) == 'p-organic') then
                mw_carma(mm) = 12._r8
             end if
             if (trim(spectype) == 's-organic') then
                mw_carma(mm) = 250._r8
             end if
             if (trim(spectype) == 'dust') then
                mw_carma(mm) = 12._r8
             end if
             if (trim(spectype) == 'seasalt') then
                mw_carma(mm) = 57._r8
             end if
          end if
          if (l == nspec(m)+1) then
             call rad_cnst_get_bin_mmr(0, m, 'a', state, pbuf, raer(mm)%fld)
             call rad_cnst_get_bin_mmr(0, m, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
          end if
          if (l == nspec(m)+2) then
             call rad_cnst_get_bin_num(0, m, 'a', state, pbuf, raer(mm)%fld)
             call rad_cnst_get_bin_num(0, m, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
          end if
          mmrcw(:ncol,:,mm) = qqcw(mm)%fld(:ncol,:)
          vmrcw(:ncol,:,mm) = qqcw(mm)%fld(:ncol,:)
          raervmr(:ncol,:,mm) = raer(mm)%fld(:ncol,:)
          !write(iulog,*) 'qqcw(mm)%fld) ', qqcw(mm)%fld(:ncol,:)
       end do
    end do
    !write(iulog,*) 'vmrcw(:,:,1) start', maxval(vmrcw(:ncol,:,1) )

    !write(iulog,*) 'mm start vmrcw, raervmr'
    ! qqcw2vrm is different from what is done in MAM, here we pass in the fields set by the qqcw and raer pointer
    ! for all the CARMA aerosols, species, mmr, and number, vmrcw (kg/kg) -> vmr
    call mmr2vmr_carma ( lchnk, vmrcw, mbar, mw_carma, ncol, loffset, rmass )
    !write(iulog,*) 'vmrcw(:,:,1) mmr', maxval(vmrcw(:,:,1))

    if (.not. is_spcam_m2005) then  ! regular CAM
       dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)   ! all adveced species no aerosols
       dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)  ! cloud borne carma aerosol species
    ! aqueous chemistry ...
    ! write(iulog,*) 'start has_sox'

    if( has_sox ) then
         call setsox(   &
              pbuf,     &
              ncol,     &
              lchnk,    &
              loffset,  &
              delt,     &
              pmid,     &
              pdel,     &
              tfld,     &
              mbar,     &
              cwat,     &
              cldfr,    &
              cldnum,   &
              airdens,  &
              invariants, &
              vmrcw,    &
              vmr,      &
              xphlwc,   &
              aqso4,    &
              aqh2so4,  &
              aqso4_h2o2, &
              aqso4_o3  &
              )

      !write(iulog,*) 'done with has_sox'
          do n = 1, nbins
            do l = 1, nspec(n)   ! not for total mass or number
                mm = bin_idx(n, l)
                call outfld( trim(fieldname_cw(mm))//'AQSO4',   aqso4(:ncol,mm),   ncol, lchnk)
                call outfld( trim(fieldname_cw(mm))//'AQH2SO4', aqh2so4(:ncol,mm), ncol, lchnk)
             end do
          end do

          call outfld( 'AQSO4_H2O2', aqso4_h2o2(:ncol), ncol, lchnk)
          call outfld( 'AQSO4_O3',   aqso4_o3(:ncol),   ncol, lchnk)
          call outfld( 'XPH_LWC',    xphlwc(:ncol,:),   ncol, lchnk )

    endif

!   Tendency due to aqueous chemistry
    dvmrdt = (vmr - dvmrdt) / delt
    dvmrcwdt = (vmrcw - dvmrcwdt) / delt

    do m = 1, gas_pcnst
       wrk(:) = 0.0_r8
       do k = 1,pver
          wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m) * adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
       end do
       name = 'AQ_'//trim(solsym(m))
       call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

    else if (is_spcam_m2005) then  ! SPCAM ECPP
! when ECPP is used, aqueous chemistry is done in ECPP,
! and not updated here.
! Minghuai Wang, 2010-02 (Minghuai.Wang@pnl.gov)

      dvmrdt = 0.0_r8
      dvmrcwdt = 0.0_r8
    endif

! do gas-aerosol exchange (h2so4, msa, nh3 condensation)

    if (h2so4_ndx > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,h2so4_ndx)
    else
       del_h2so4_aeruptk(:,:) = 0.0_r8
    endif


    !call t_startf('modal_gas-aer_exchng')

    !if ( sulfeq_idx>0 ) then
    !   call pbuf_get_field( pbuf, sulfeq_idx, sulfeq )
    !else
    !   nullify( sulfeq )
    !endif
    !write(iulog,*) 'start carma_aero_gasaerexch_sub'
    ! need to transform raer to raervmr from CARMA, routine requires vmr, note number wil not be changed here
    call mmr2vmr_carma ( lchnk, raervmr, mbar, mw_carma, ncol, loffset, rmass)
    !write(iulog,*) 'mm start raervmr done'

    call carma_aero_gasaerexch_sub(             &
          pbuf, lchnk,    ncol,     nstep,      &
          loffset,            delt, mbar ,      &
          tfld,     pmid,     pdel,             &
          qh2o,               troplev,          &
          vmr,                raervmr,          &
          wetr_n     )

    !if (h2so4_ndx > 0) then
    !   del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,h2so4_ndx) - del_h2so4_aeruptk(1:ncol,:)
    !endif

    !call t_stopf('modal_gas-aer_exchng')


    ! note vmr2qqcw does not change qqcw pointer (different than in MAM)
    call vmr2mmr_carma ( lchnk, vmrcw, mbar, mw_carma, ncol, loffset, rmass )

    !vmrcw in kg/kg
    ! change pointer value for total mmr and number. In order to do this correctly
    ! only mass has to be added to each bin (not number). This will require redistributing
    ! mass to different bins. Here, we change both mass and number until we have a better
    ! solution.
    delta_so4mass(:,:,:) = 0.0_r8
    do m = 1, nbins
       do l = 1, nspec(m)  ! for sulfate only
          mm = bin_idx(m, l)
         ! sulfate mass that needs to be added to the total mass
          call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
          if (trim(spectype) == 'sulfate') then
              ! only do loop if vmrcw has changed
              do k=1,pver
                 do i=1,ncol
                  if (vmrcw(i,k,mm) .gt. mmrcw(i,k,mm) .and. mmrcw(i,k,mm) /= 0.0_r8)  then
                   delta_so4mass(i,k,mm) = ( vmrcw(i,k,mm) - mmrcw(i,k,mm) )
                  else
                    delta_so4mass(i,k,mm) = 0.0_r8
                  end if
                 end do
              end do
         end if
       end do
    end do

    ! adjust total mass and number
    old_total_mass = 0.0_r8
    new_total_mass = 0.0_r8
    old_total_number = 0.0_r8
    do m = 1, nbins
       do k=1,pver
         do i=1,ncol
             l = nspec(m) + 1
             mm = bin_idx(m, l)
             if (delta_so4mass(i,k,mm) .gt. 0.0_r8) then
               call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
               old_total_mass = vmrcw(i,k,mm)
               vmrcw(i,k,mm) = vmrcw(i,k,mm) + delta_so4mass(i,k,mm)
               new_total_mass = vmrcw(i,k,mm)
               !  rmass is used to scale number
               l = nspec(m) + 2
               mm = bin_idx(m, l)
               old_total_number = vmrcw(i,k,mm)
               vmrcw(i,k,mm) = new_total_mass / rmass(m)
            end if
           enddo
        enddo
     enddo

    ! Is the loop here needed?
    do m = 1, nbins
       do l = 1, nspec(m) + 2 ! for sulfate only
          mm = bin_idx(m, l)
          qqcw(mm)%fld(:ncol,:) = vmrcw(:ncol,:,mm)
          call outfld( trim(fieldname_cw(mm)), qqcw(mm)%fld(:ncol,:), ncol, lchnk)
       end do
    end do


  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

  end subroutine aero_model_emissions


  !===============================================================================
  !===============================================================================
  ! private methods


  !=============================================================================
  !=============================================================================
  subroutine surf_area_dens( state, pbuf, ncol, mmr, beglev, endlev, sad, reff, sfc, dm_aer )
    use mo_constants, only: pi

    ! dummy args
    type(physics_state),    intent(in) :: state           ! Physics state variables
    type(physics_buffer_desc), pointer :: pbuf(:)
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: mmr(:,:,:)
    integer,  intent(in)  :: beglev(:)
    integer,  intent(in)  :: endlev(:)
    real(r8), intent(out) :: sad(:,:)    ! bulk surface area density in cm2/cm3 from beglev to endlev, zero elsewhere
    real(r8), intent(out) :: reff(:,:)   ! bulk effective radius in cm from beglev to endlev, zero elsewhere
    real(r8), optional, intent(out) :: sfc(:,:,:) ! surface area density per bin
    real(r8), optional, intent(out) :: dm_aer(:,:,:) ! diameter per bin

    ! local vars
    real(r8), pointer, dimension(:,:) :: reffaer ! bulk effective radius in cm
    real(r8), pointer, dimension(:,:) :: cmass,tmass ! carma element chemical and total mass
    real(r8) :: sad_bin(pcols,pver,nbins)
    integer  :: err, icol, ilev, ibin, ispec, reff_pbf_ndx
    real(r8) :: chm_mass, tot_mass
    character(len=32) :: bin_name   ! CARMA bin name
    character(len=32) :: spectype
    real(r8), pointer :: wetr(:,:)  ! CARMA bin wet radius in cm
    real(r8), pointer :: sad_carma(:,:)  ! CARMA bin wet surface area density in cm2/cm3
    real(r8), pointer :: aer_bin_mmr(:,:)

    !
    ! Compute surface aero for each bin.
    ! Total over all bins as the surface area for chemical reactions.
    !
    reff_pbf_ndx = pbuf_get_index("REFFAER",errcode=err) ! CARMA aerosol effective radius
    call pbuf_get_field(pbuf, reff_pbf_ndx, reffaer )

    sad = 0._r8
    sad_bin = 0._r8
    reff = 0._r8

    do ibin=1,nbins ! loop over aerosol bins
      call rad_cnst_get_info_by_bin(0, ibin, bin_name=bin_name)
      call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_sad"),sad_carma)
      call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_wetr"),wetr)
      if (present(dm_aer)) then
         dm_aer(:ncol,:,ibin) = 2._r8 * wetr(:ncol,:) ! convert wet radius (cm) to wet diameter (cm)
      endif
      sad_bin(:ncol,:,ibin) = sad_carma(:ncol,:) ! cm^2/cm^3
    end do

    do icol = 1,ncol
      do ilev = beglev(icol),endlev(icol)
        do ibin=1,nbins ! loop over aerosol bins
          !
          ! compute a mass weighting of the number
          !
          tot_mass = 0._r8
          chm_mass = 0._r8
          do ispec=1,nspec(ibin)

             call rad_cnst_get_bin_mmr_by_idx(0, ibin, ispec, 'a', state, pbuf, aer_bin_mmr)

             tot_mass = tot_mass + aer_bin_mmr(icol,ilev)

             call rad_cnst_get_bin_props_by_idx(0, ibin, ispec, spectype=spectype)

             if ( trim(spectype) == 'sulfate'   .or. &
                trim(spectype) == 's-organic' .or. &
                trim(spectype) == 'p-organic' .or. &
                trim(spectype) == 'black-c'   .or. &
                trim(spectype) == 'ammonium') then
                chm_mass = chm_mass + aer_bin_mmr(icol,ilev)
             end if

          end do
          if ( tot_mass > 0._r8 ) then
         ! surface area density
            sad_bin(icol,ilev,ibin) = chm_mass / tot_mass * sad_bin(icol,ilev,ibin) ! cm^2/cm^3
          else
            sad_bin(icol,ilev,ibin) = 0._r8
          end if
        end do
        sad(icol,ilev) = sum(sad_bin(icol,ilev,:))
        reff(icol,ilev) = reffaer(icol,ilev)

       end do
    end do

    if (present(sfc)) then
       sfc(:,:,:) = sad_bin(:,:,:)
    endif

  end subroutine surf_area_dens

  !===============================================================================
  !===============================================================================
  subroutine carma_aero_bcscavcoef_init ( pbuf2d )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Computes lookup table for aerosol impaction/interception scavenging rates
    !
    ! Authors: R. Easter
    ! Simone Tilmes Nov 2021
    ! added modifications for bin model, assuming sigma = 1.
    !
    !-----------------------------------------------------------------------

    use shr_kind_mod,    only: r8 => shr_kind_r8
    use cam_abortutils,  only: endrun
    use mo_constants, only:  pi
    use ppgrid,          only: begchunk

    implicit none

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !   local variables
    integer nnfit_maxd
    parameter (nnfit_maxd=27)

    integer m, i, l, jgrow, jdens, jpress, jtemp, nnfit
    integer lunerr

    character(len=32) :: bin_name
    character(len=32) :: spectype

    real(r8), pointer :: rmass_ptr(:)   ! CARMA rmass fixed per bin
    real(r8), allocatable :: rmass(:)     ! CARMA rmass

    real(r8) dg0, dg0_cgs, press, dg0_base, &
         rhodryaero, rhowetaero, rhowetaero_cgs, rmserr, &
         scavratenum, scavratevol, sigmag,                &
         temp, wetdiaratio, wetvolratio
    real(r8) :: specdens
    real(r8) aafitnum(1), xxfitnum(1,nnfit_maxd), yyfitnum(nnfit_maxd)
    real(r8) aafitvol(1), xxfitvol(1,nnfit_maxd), yyfitvol(nnfit_maxd)


    allocate(scavimptblnum(nimptblgrow_mind:nimptblgrow_maxd, nbins))
    allocate(scavimptblvol(nimptblgrow_mind:nimptblgrow_maxd, nbins))

    lunerr = iulog
    dlndg_nimptblgrow = log( 1.25_r8 )

    allocate ( rmass(nbins) )
    ! bin model: main loop over aerosol bins

    modeloop: do m = 1, nbins
       !write(*,*) 'mloop start ',m
       ! r(m) is the dry bin radius
       ! taken here from CARMA pbuf field
       ! get bin info
       call rad_cnst_get_info_by_bin(0, m, bin_name=bin_name)

       call pbuf_get_field(pbuf2d, begchunk, pbuf_get_index(trim(bin_name)//"_rmass"),rmass_ptr)
       rmass(m) = rmass_ptr(1)

       !   for setting up the lookup table, use the dry density of the first
      ! get specdens from sulfate (check)
       do l = 1, nspec(m)

        call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype, density_aer=specdens)

       !   chemical component of the aerosol type (which currently will be so4)
       !   For  CARMA, rmass per bin stays the same, while dry radius varies when the particle density varies
       !   rmass = 4/3 * Pi * density * dry radius
       !   We assume a fixed specie density
        if (trim(spectype) == 'sulfate') then
           rhodryaero = specdens
        end if
       end do
       dg0_base = 2._r8 * (0.75_r8*rmass(m) / pi  / (1.0e-3_r8*rhodryaero)) **(0.33_r8)    ! specdens kg/m3 to g/cm3, convert from radiust to diameter

           !sigmag = sigmag_amode(mode)
           !dg0_base = dcen_sect(m,n)*exp( -1.5*((log(sigmag))**2) )
           ! for bin approach sigma assumed to be 1., dg0_base equal dry radius
       sigmag = 1._r8


       !st rhodryaero = specdens_amode(1,mode)

       growloop: do jgrow = nimptblgrow_mind, nimptblgrow_maxd

          wetdiaratio = exp( jgrow*dlndg_nimptblgrow )
          !dg0 = dgnum_amode(mode)*wetdiaratio
          dg0 = dg0_base*wetdiaratio
          !st write(*,*) 'm,l,dg0 ',m,l,dg0

          wetvolratio = exp( jgrow*dlndg_nimptblgrow*3._r8 )
          rhowetaero = 1.0_r8 + (rhodryaero-1.0_r8)/wetvolratio
          rhowetaero = min( rhowetaero, rhodryaero )

          !
          !   compute impaction scavenging rates at 1 temp-press pair and save
          !
          nnfit = 0

          temp = 273.16_r8
          press = 0.75e6_r8   ! dynes/cm2
          rhowetaero = rhodryaero

          ! CARMA dry radius is in  cm
          !dg0_cgs = dg0*1.0e2_r8   ! m to cm
          dg0_cgs = dg0   ! CARMA  radius / diameter is already in cm

          rhowetaero_cgs = rhowetaero*1.0e-3_r8   ! kg/m3 to g/cm3


          call calc_1_impact_rate( &
               dg0_cgs, sigmag, rhowetaero_cgs, temp, press, &
               scavratenum, scavratevol, lunerr )


          nnfit = nnfit + 1
          if (nnfit > nnfit_maxd) then
             write(lunerr,9110)
             call endrun()
          end if
9110      format( '*** subr. carma_aero_bcscavcoef_init -- nnfit too big' )

          xxfitnum(1,nnfit) = 1._r8
          yyfitnum(nnfit) = log( scavratenum )

          xxfitvol(1,nnfit) = 1._r8
          yyfitvol(nnfit) = log( scavratevol )

          !
          ! skip mlinfit stuff because scav table no longer has dependencies on
          !    air temp, air press, and particle wet density
          ! just load the log( scavrate--- ) values
          !
          !!
          !!   do linear regression
          !!    log(scavrate) = a1 + a2*log(wetdens)
          !!
          !     call mlinft( xxfitnum, yyfitnum, aafitnum, nnfit, 1, 1, rmserr )
          !     call mlinft( xxfitvol, yyfitvol, aafitvol, nnfit, 1, 1, rmserr )
          !
          !     scavimptblnum(jgrow,mode) = aafitnum(1)
          !     scavimptblvol(jgrow,mode) = aafitvol(1)

         !depends on both bins and different species
          scavimptblnum(jgrow,m) = yyfitnum(1)
          scavimptblvol(jgrow,m) = yyfitvol(1)

       enddo growloop
    enddo modeloop

    return
  end subroutine carma_aero_bcscavcoef_init

  !===============================================================================
  !===============================================================================


  !===============================================================================
  subroutine carma_aero_bcscavcoef_get( m, ncol, isprx, wetr, dryr, scavcoefnum, scavcoefvol, pbuf )
    !  need to go through both bins and species
    ! need dry radius and wet radius

    !-----------------------------------------------------------------------

    use mo_constants, only:  pi

    implicit none

    integer,intent(in) :: m, ncol
    logical,intent(in):: isprx(pcols,pver)
    ! wet radius per bin dgn_awet -> wetr
    real(r8), intent(in) :: dryr(pcols,pver)
    real(r8), intent(in) :: wetr(pcols,pver)
    real(r8), intent(out) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer i, k, jgrow, l
    real(r8) dumdgratio, xgrow, dumfhi, dumflo, scavimpvol, scavimpnum, dg0_base, specdens, rhodryaero

    character(len=32) :: spectype
    character(len=32) :: bin_name

    real(r8), allocatable :: rmass(:)     ! CARMA rmass
    real(r8), pointer :: rmass_ptr(:)   ! CARMA rmass fixed per bin

    allocate ( rmass(nbins) )
    ! bin model: main loop over aerosol bins

    ! get bin info
     call rad_cnst_get_info_by_bin(0, m, bin_name=bin_name)
     call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_rmass"),rmass_ptr)
     rmass(m) = rmass_ptr(1)


    ! get rmass and specdens for sulfate
    do l = 1, nspec(m)
        call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype, density_aer=specdens)

       !   chemical component of the aerosol type (which currently will be so4)
       !   For  CARMA, rmass per bin stays the same, while dry radius varies when the particle density varies
       !   rmass = 4/3 * Pi * density * dry radius
       !   We assume a fixed specie density
        if (trim(spectype) == 'sulfate') then
           rhodryaero = specdens
        end if
    end do
    dg0_base = 2._r8 * (0.75_r8*rmass(m) / pi  / (1.0e-3_r8*rhodryaero)) **(0.33_r8)    ! specdens kg/m3 to g/cm3, convert from radiust to diameter
    !rg0_base = (0.75*rmass(m) / pi  / (1.0e-3_r8*specdens)) **(0.33_r8)    ! specdens kg/m3 to g/cm3

    do k = 1, pver
       do i = 1, ncol

          ! do only if no precip
          if ( isprx(i,k) .and. dryr(i,k).gt.0._r8) then
             !
             ! interpolate table values using log of (actual-wet-size)/(base-dry-size)

             ! dumdgratio = dgn_awet(i,k,m)/dgnum_amode(m)
             ! dgnum_amode(m) is the rg0_base radius.

             dumdgratio = wetr(i,k)/dg0_base

              if ((dumdgratio >= 0.99_r8) .and. (dumdgratio <= 1.01_r8)) then
                 scavimpvol = scavimptblvol(0,m)
                 scavimpnum = scavimptblnum(0,m)
              else
                xgrow = log( dumdgratio ) / dlndg_nimptblgrow
                jgrow = int( xgrow )
                if (xgrow < 0._r8) jgrow = jgrow - 1
                if (jgrow < nimptblgrow_mind) then
                   jgrow = nimptblgrow_mind
                   xgrow = jgrow
                else
                   jgrow = min( jgrow, nimptblgrow_maxd-1 )
                end if

                dumfhi = xgrow - jgrow
                dumflo = 1._r8 - dumfhi

                scavimpvol = dumflo*scavimptblvol(jgrow,m) + &
                     dumfhi*scavimptblvol(jgrow+1,m)
                scavimpnum = dumflo*scavimptblnum(jgrow,m) + &
                     dumfhi*scavimptblnum(jgrow+1,m)

             end if

             ! impaction scavenging removal amount for volume
             scavcoefvol(i,k) = exp( scavimpvol )
             ! impaction scavenging removal amount to number
             scavcoefnum(i,k) = exp( scavimpnum )

             ! scavcoef = impaction scav rate (1/h) for precip = 1 mm/h
             ! scavcoef = impaction scav rate (1/s) for precip = pfx_inrain
             ! (scavcoef/3600) = impaction scav rate (1/s) for precip = 1 mm/h
             ! (pfx_inrain*3600) = in-rain-area precip rate (mm/h)
             ! impactrate = (scavcoef/3600) * (pfx_inrain*3600)
          else
             scavcoefvol(i,k) = 0._r8
             scavcoefnum(i,k) = 0._r8
          end if

       end do
    end do

    return
  end subroutine carma_aero_bcscavcoef_get

  !===============================================================================
  subroutine calc_1_impact_rate(             &
       dg0, sigmag, rhoaero, temp, press, &
       scavratenum, scavratevol, lunerr )
    !
    !   routine computes a single impaction scavenging rate
    !    for precipitation rate of 1 mm/h
    !
    !   dg0 = geometric mean diameter of aerosol number size distrib. (for CARMA it is the dry radius) (cm)
    !   sigmag = geometric standard deviation of size distrib.
    !   rhoaero = density of aerosol particles (g/cm^3)
    !   temp = temperature (K)
    !   press = pressure (dyne/cm^2)
    !   scavratenum = number scavenging rate (1/h)
    !   scavratevol = volume or mass scavenging rate (1/h)
    !   lunerr = logical unit for error message
    !
    use shr_kind_mod, only: r8 => shr_kind_r8
    use mo_constants, only: boltz_cgs, pi, rhowater => rhoh2o_cgs, &
                           gravity => gravity_cgs, rgas => rgas_cgs

   implicit none

   !   subr. parameters
   integer lunerr
   real(r8) dg0, sigmag, rhoaero, temp, press, scavratenum, scavratevol

   !   local variables
   integer nrainsvmax
   parameter (nrainsvmax=50)
   real(r8) rrainsv(nrainsvmax), xnumrainsv(nrainsvmax),&
        vfallrainsv(nrainsvmax)

   integer naerosvmax
   parameter (naerosvmax=51)
   real(r8) aaerosv(naerosvmax), &
        ynumaerosv(naerosvmax), yvolaerosv(naerosvmax)

   integer i, ja, jr, na, nr
   real(r8) a, aerodiffus, aeromass, ag0, airdynvisc, airkinvisc
   real(r8) anumsum, avolsum, cair, chi
   real(r8) d, dr, dum, dumfuchs, dx
   real(r8) ebrown, eimpact, eintercept, etotal, freepath
   real(r8) precip, precipmmhr, precipsum
   real(r8) r, rainsweepout, reynolds, rhi, rhoair, rlo, rnumsum
   real(r8) scavsumnum, scavsumnumbb
   real(r8) scavsumvol, scavsumvolbb
   real(r8) schmidt, sqrtreynolds, sstar, stokes, sx
   real(r8) taurelax, vfall, vfallstp
   real(r8) x, xg0, xg3, xhi, xlo, xmuwaterair


   rlo = .005_r8
   rhi = .250_r8
   dr = 0.005_r8
   nr = 1 + nint( (rhi-rlo)/dr )
   if (nr > nrainsvmax) then
      write(lunerr,9110)
      call endrun()
   end if

9110 format( '*** subr. calc_1_impact_rate -- nr > nrainsvmax' )

   precipmmhr = 1.0_r8
   precip = precipmmhr/36000._r8

! if dg0 the diameter, than ag0 equals the radius
   ag0 = dg0/2._r8
  if (sigmag.ne.1._r8) then
   sx = log( sigmag )
   xg0 = log( ag0 )
   xg3 = xg0 + 3._r8*sx*sx

   xlo = xg3 - 4._r8*sx
   xhi = xg3 + 4._r8*sx
   dx = 0.2_r8*sx

   dx = max( 0.2_r8*sx, 0.01_r8 )
   xlo = xg3 - max( 4._r8*sx, 2._r8*dx )
   xhi = xg3 + max( 4._r8*sx, 2._r8*dx )

   na = 1 + nint( (xhi-xlo)/dx )
   if (na > naerosvmax) then
      write(lunerr,9120)
      call endrun()
   end if
  else
   na = 1
   a = ag0
  end if

9120 format( '*** subr. calc_1_impact_rate -- na > naerosvmax' )

   !   air molar density
   cair = press/(rgas*temp)
   !   air mass density
   rhoair = 28.966_r8*cair
   !   molecular freepath
   freepath = 2.8052e-10_r8/cair
   !   air dynamic viscosity
   airdynvisc = 1.8325e-4_r8 * (416.16_r8/(temp+120._r8)) *    &
        ((temp/296.16_r8)**1.5_r8)
   !   air kinemaic viscosity
   airkinvisc = airdynvisc/rhoair
   !   ratio of water viscosity to air viscosity (from Slinn)
   xmuwaterair = 60.0_r8

   !
   !   compute rain drop number concentrations
   !    rrainsv = raindrop radius (cm)
   !    xnumrainsv = raindrop number concentration (#/cm^3)
   !            (number in the bin, not number density)
   !    vfallrainsv = fall velocity (cm/s)
   !
   precipsum = 0._r8
   do i = 1, nr
      r = rlo + (i-1)*dr
      rrainsv(i) = r
      xnumrainsv(i) = exp( -r/2.7e-2_r8 )

      d = 2._r8*r
      if (d <= 0.007_r8) then
         vfallstp = 2.88e5_r8 * d**2._r8
      else if (d <= 0.025_r8) then
         vfallstp = 2.8008e4_r8 * d**1.528_r8
      else if (d <= 0.1_r8) then
         vfallstp = 4104.9_r8 * d**1.008_r8
      else if (d <= 0.25_r8) then
         vfallstp = 1812.1_r8 * d**0.638_r8
      else
         vfallstp = 1069.8_r8 * d**0.235_r8
      end if

      vfall = vfallstp * sqrt(1.204e-3_r8/rhoair)
      vfallrainsv(i) = vfall
      precipsum = precipsum + vfall*(r**3)*xnumrainsv(i)
   end do
   precipsum = precipsum*pi*1.333333_r8

   rnumsum = 0._r8
   do i = 1, nr
      xnumrainsv(i) = xnumrainsv(i)*(precip/precipsum)
      rnumsum = rnumsum + xnumrainsv(i)
   end do

   !
   !   compute aerosol concentrations
   !    aaerosv = particle radius (cm)
   !    fnumaerosv = fraction of total number in the bin (--)
   !    fvolaerosv = fraction of total volume in the bin (--)
   !


   anumsum = 0._r8
   avolsum = 0._r8
   ynumaerosv(:) = 1._r8
   yvolaerosv(:) = 1._r8
   aaerosv(:) = a
  if (na.ne.1) then
   do i = 1, na
      x = xlo + (i-1)*dx
      a = exp( x )
      aaerosv(i) = a
      dum = (x - xg0)/sx
      ynumaerosv(i) = exp( -0.5_r8*dum*dum )
      yvolaerosv(i) = ynumaerosv(i)*1.3333_r8*pi*a*a*a
      anumsum = anumsum + ynumaerosv(i)
      avolsum = avolsum + yvolaerosv(i)
   end do

   do i = 1, na
      ynumaerosv(i) = ynumaerosv(i)/anumsum
      yvolaerosv(i) = yvolaerosv(i)/avolsum
   end do
  end if


   !
   !   compute scavenging
   !
   scavsumnum = 0._r8
   scavsumvol = 0._r8
   !
   !   outer loop for rain drop radius
   !
   jr_loop: do jr = 1, nr

      r = rrainsv(jr)
      vfall = vfallrainsv(jr)

      reynolds = r * vfall / airkinvisc
      sqrtreynolds = sqrt( reynolds )

      !
      !   inner loop for aerosol particle radius
      !
      scavsumnumbb = 0._r8
      scavsumvolbb = 0._r8

      ja_loop: do ja = 1, na

         a = aaerosv(ja)

         chi = a/r

         dum = freepath/a
         dumfuchs = 1._r8 + 1.246_r8*dum + 0.42_r8*dum*exp(-0.87_r8/dum)
         taurelax = 2._r8*rhoaero*a*a*dumfuchs/(9._r8*rhoair*airkinvisc)


         aeromass = 4._r8*pi*a*a*a*rhoaero/3._r8
         aerodiffus = boltz_cgs*temp*taurelax/aeromass

         schmidt = airkinvisc/aerodiffus
         stokes = vfall*taurelax/r

         ebrown = 4._r8*(1._r8 + 0.4_r8*sqrtreynolds*(schmidt**0.3333333_r8)) /  &
              (reynolds*schmidt)

         dum = (1._r8 + 2._r8*xmuwaterair*chi) /         &
              (1._r8 + xmuwaterair/sqrtreynolds)
         eintercept = 4._r8*chi*(chi + dum)

         dum = log( 1._r8 + reynolds )
         sstar = (1.2_r8 + dum/12._r8) / (1._r8 + dum)
         eimpact = 0._r8
         if (stokes > sstar) then
            dum = stokes - sstar
            eimpact = (dum/(dum+0.6666667_r8)) ** 1.5_r8
         end if

         etotal = ebrown + eintercept + eimpact
         etotal = min( etotal, 1.0_r8 )

         rainsweepout = xnumrainsv(jr)*4._r8*pi*r*r*vfall

         scavsumnumbb = scavsumnumbb + rainsweepout*etotal*ynumaerosv(ja)
         scavsumvolbb = scavsumvolbb + rainsweepout*etotal*yvolaerosv(ja)

      enddo ja_loop

      scavsumnum = scavsumnum + scavsumnumbb
      scavsumvol = scavsumvol + scavsumvolbb

   enddo jr_loop

   scavratenum = scavsumnum*3600._r8
   scavratevol = scavsumvol*3600._r8

   return
 end subroutine calc_1_impact_rate

  !=============================================================================
  subroutine mmr2vmr_carma(lchnk, vmr, mbar, mw_carma, ncol, im, rmass)
    !-----------------------------------------------------------------
    !   ... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    !st use chem_mods, only : adv_mass, gas_pcnst

    implicit none

    !-----------------------------------------------------------------
    !   ... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: rmass(nbins)
    real(r8), intent(in)    :: mw_carma(ncnst_tot)
    real(r8), intent(inout) :: vmr(ncol,pver,ncnst_tot)
    real(r8)                :: vmr_total(ncol,pver)

    !-----------------------------------------------------------------
    !   ... Local variables
    !-----------------------------------------------------------------
    integer :: k, m, mm, l

    do m = 1, nbins
       vmr_total(:ncol,:) = 0._r8
       do l = 1, nspec(m)   ! for each species, not total mmr or number, information of mw are missing
             mm = bin_idx(m, l)
             do k=1,pver
                vmr(:ncol,k,mm) = mbar(:ncol,k) * vmr(:ncol,k,mm) / mw_carma(mm)
             end do
             vmr_total(:ncol,:) = vmr_total(:ncol,:) +  vmr(:ncol,:,mm)
       end do
       ! calculate mmr
       l =nspec(m) + 1
       mm = bin_idx(m, l)
       vmr(:ncol,:,mm) = vmr_total(:ncol,:)
       ! calculate number
       l =nspec(m) + 2
       mm = bin_idx(m, l)
       vmr(:ncol,:,mm) = vmr_total(:ncol,:) / rmass(m)
    end do
  end subroutine mmr2vmr_carma
    !=============================================================================

  !=============================================================================
  subroutine vmr2mmr_carma ( lchnk, vmr, mbar, mw_carma, ncol, im, rmass )
    !-----------------------------------------------------------------
    !   ... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------
    !   ... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: rmass(nbins)
    real(r8), intent(inout)    :: vmr(ncol,pver,ncnst_tot)
    real(r8), intent(in)    :: mw_carma(ncnst_tot)
    real(r8)                :: vmr_total(ncol,pver)

    !-----------------------------------------------------------------
    !   ... Local variables
    !-----------------------------------------------------------------
    integer :: k, m, mm, l
    !-----------------------------------------------------------------
    !   ... The non-group species
    !-----------------------------------------------------------------
    do m = 1, nbins
       vmr_total(:ncol,:) = 0._r8
       do l = 1, nspec(m)   ! for each species, not total mmr or number, information of mw are missing
             mm = bin_idx(m, l)
             do k=1,pver
                vmr(:ncol,k,mm) = mw_carma(mm) * vmr(:ncol,k,mm) / mbar(:ncol,k)
             end do
             vmr_total(:ncol,:) = vmr_total(:ncol,:) +  vmr(:ncol,:,mm)
       end do
       ! calculate mmr
       l =nspec(m) + 1
       mm = bin_idx(m, l)
       vmr(:ncol,:,mm) = vmr_total(:ncol,:)
       ! calculate number
       l =nspec(m) + 2
       mm = bin_idx(m, l)
       vmr(:ncol,:,mm) = vmr_total(:ncol,:) / rmass(m)
    end do

   end subroutine vmr2mmr_carma

end module aero_model
