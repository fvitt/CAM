module carma_aero_convproc
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to aerosol/trace-gas convective cloud processing scheme
!
! currently these routines assume stratiform and convective clouds only interact
! through the detrainment of convective cloudborne material into stratiform clouds
!
! thus the stratiform-cloudborne aerosols (in the qqcw array) are not processed
! by the convective up/downdrafts, but are affected by the detrainment
!
! Author: R. C. Easter
!
!---------------------------------------------------------------------------------

use shr_kind_mod,    only: r8=>shr_kind_r8, shr_kind_cs
use spmd_utils,      only: masterproc
use physconst,       only: gravit, rair
use ppgrid,          only: pver, pcols, pverp
use constituents,    only: pcnst, cnst_name, cnst_get_ind
use constituents,    only: cnst_species_class, cnst_spec_class_aerosol
use phys_control,    only: phys_getopts

use physics_types,   only: physics_state, physics_ptend
use physics_buffer,  only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

use time_manager,    only: get_nstep
use cam_history,     only: outfld, addfld, add_default, horiz_only
use cam_logfile,     only: iulog
use cam_abortutils,  only: endrun
use rad_constituents,only: rad_cnst_get_info, rad_cnst_get_info_by_bin, rad_cnst_get_info_by_bin_spec, &
                           rad_cnst_get_bin_props_by_idx, rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_mmr, &
                           rad_cnst_get_bin_num

use carma_aerosol_properties_mod, only: carma_aerosol_properties
use carma_aerosol_state_mod, only: carma_aerosol_state

implicit none
private
save

public :: &
   ma_convproc_register, &
   ma_convproc_init,     &
   ma_convproc_intr,     &
   ma_convproc_readnl

! description of bin aerosols
integer, public, protected :: nspec_max = 0
integer, public, protected :: nbins = 0
integer, public, protected, allocatable :: nspec(:)

! local indexing for bins
integer, allocatable :: bin_idx(:,:) ! table for local indexing of modal aero number and mmr
integer :: ncnst_tot                  ! total number of bin species, mmr, and number conc
integer :: ncnst_extd                 ! twice total number of bin species, mmr and number conc
integer :: pcnst_extd                 ! twice total number of bin species, mmr and number conc

! Indices for CARMA species in the ptend%q array.  Needed for prognostic aerosol case.
logical, allocatable :: bin_cnst_lq(:,:)
integer, allocatable :: bin_cnst_idx(:,:)

! ptr2d_t is used to create arrays of pointers to 2D fields
type ptr2d_t
   real(r8), pointer :: fld(:,:)
end type ptr2d_t

logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
                               ! in the ptend object

! namelist options
! NOTE: These are the defaults for CAM6.
logical, protected, public :: deepconv_wetdep_history = .true.
logical, protected, public :: convproc_do_deep = .true.
! NOTE: These are the defaults for the Eaton/Wang parameterization.
logical, protected, public :: convproc_do_evaprain_atonce = .false.
real(r8), protected, public :: convproc_pom_spechygro = -1._r8
real(r8), protected, public :: convproc_wup_max       = 4.0_r8

logical, parameter :: use_cwaer_for_activate_maxsat = .false.
logical, parameter :: apply_convproc_tend_to_ptend = .true.

real(r8) :: hund_ovr_g ! = 100.0_r8/gravit
!  used with zm_conv mass fluxes and delta-p
!     for mu = [mbar/s],   mu*hund_ovr_g = [kg/m2/s]
!     for dp = [mbar] and q = [kg/kg],   q*dp*hund_ovr_g = [kg/m2]

!  method1_activate_nlayers = number of layers (including cloud base) where activation is applied
integer, parameter  :: method1_activate_nlayers = 2
!  method2_activate_smaxmax = the uniform or peak supersat value (as 0-1 fraction = percent*0.01)
real(r8), parameter :: method2_activate_smaxmax = 0.003_r8

!  method_reduce_actfrac = 1 -- multiply activation fractions by factor_reduce_actfrac
!                               (this works ok with convproc_method_activate = 1 but not for ... = 2)
!                        = 2 -- do 2 iterations to get an overall reduction by factor_reduce_actfrac
!                               (this works ok with convproc_method_activate = 1 or 2)
!                        = other -- do nothing involving reduce_actfrac
integer, parameter  :: method_reduce_actfrac = 0
real(r8), parameter :: factor_reduce_actfrac = 0.5_r8

!  convproc_method_activate - 1=apply abdulrazzak-ghan to entrained aerosols for lowest nlayers
!                             2=do secondary activation with prescribed supersat
integer, parameter :: convproc_method_activate = 2

logical :: convproc_do_aer

! physics buffer indices
integer :: fracis_idx         = 0

integer :: rprddp_idx         = 0
integer :: rprdsh_idx         = 0
integer :: nevapr_shcu_idx    = 0
integer :: nevapr_dpcu_idx    = 0

integer :: icwmrdp_idx        = 0
integer :: icwmrsh_idx        = 0
integer :: sh_frac_idx        = 0
integer :: dp_frac_idx        = 0

integer :: zm_mu_idx          = 0
integer :: zm_eu_idx          = 0
integer :: zm_du_idx          = 0
integer :: zm_md_idx          = 0
integer :: zm_ed_idx          = 0
integer :: zm_dp_idx          = 0
integer :: zm_dsubcld_idx     = 0
integer :: zm_jt_idx          = 0
integer :: zm_maxg_idx        = 0
integer :: zm_ideep_idx       = 0

integer :: cmfmc_sh_idx       = 0
integer :: sh_e_ed_ratio_idx  = 0

integer :: istat

logical, parameter :: debug=.false.

character(len=32), allocatable :: fieldname(:)    ! names for interstitial output fields
character(len=32), allocatable :: fieldname_cw(:)    ! names for cloud_borne output fields
integer, allocatable :: ibl(:)    !  index that points to the larger bin for evaporation
integer, allocatable :: imx_bl(:)    ! index used to map pure sulfate bin to mixed sulfate bin
integer, allocatable :: imx_mmr_bl(:)    ! index used to map pure sulfate bin to mixed sulfate bin for mmr
integer, allocatable :: imx_num_bl(:)    ! index used to map pure sulfate bin to mixed sulfate bin for num

type(carma_aerosol_properties), pointer :: aero_props_obj => null()

!=========================================================================================
contains
!=========================================================================================

subroutine ma_convproc_register

end subroutine ma_convproc_register

!=========================================================================================
subroutine ma_convproc_readnl(nlfile)

  use namelist_utils, only: find_group_name
  use spmd_utils,     only: mpicom, masterprocid, mpi_real8, mpi_logical

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'ma_convproc_readnl'

  namelist /aerosol_convproc_opts/ deepconv_wetdep_history, convproc_do_deep, &
       convproc_do_evaprain_atonce, convproc_pom_spechygro, convproc_wup_max

  ! Read namelist
  if (masterproc) then
     open( newunit=unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'aerosol_convproc_opts', status=ierr)
     if (ierr == 0) then
        read(unitn, aerosol_convproc_opts, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
  end if

  ! Broadcast namelist variables
  call mpi_bcast( deepconv_wetdep_history,  1, mpi_logical, masterprocid, mpicom, ierr)
  call mpi_bcast( convproc_do_deep,  1, mpi_logical, masterprocid, mpicom, ierr)
  call mpi_bcast( convproc_do_evaprain_atonce,  1, mpi_logical, masterprocid, mpicom, ierr)
  call mpi_bcast( convproc_pom_spechygro,  1, mpi_real8, masterprocid, mpicom, ierr)
  call mpi_bcast( convproc_wup_max,  1, mpi_real8, masterprocid, mpicom, ierr)

  if (masterproc) then
     write(iulog,*) subname//': deepconv_wetdep_history = ',deepconv_wetdep_history
     write(iulog,*) subname//': convproc_do_deep = ',convproc_do_deep
     write(iulog,*) subname//': convproc_do_evaprain_atonce = ',convproc_do_evaprain_atonce
     write(iulog,*) subname//': convproc_pom_spechygro = ',convproc_pom_spechygro
     write(iulog,*) subname//': convproc_wup_max = ', convproc_wup_max
  end if

end subroutine ma_convproc_readnl

!=========================================================================================

subroutine ma_convproc_init

   integer :: n, l, ll, m, ii
   integer :: imx, imx_num, imx_mmr, ipr, ipr_num, ipr_mmr
   integer :: lptr      = -1
   integer :: idxtmp    = -1
   integer :: npass_calc_updraft
   logical :: history_aerosol
   character(len=32) :: spectype
   character(len=32) :: bin_name
   character(len=32) :: bin_name_l    ! bin name of the larger bin


   call phys_getopts( history_aerosol_out=history_aerosol, &
        convproc_do_aer_out = convproc_do_aer )

   call addfld('SH_MFUP_MAX', horiz_only, 'A', 'kg/m2', &
               'Shallow conv. column-max updraft mass flux' )
   call addfld('SH_WCLDBASE', horiz_only, 'A', 'm/s', &
               'Shallow conv. cloudbase vertical velocity' )
   call addfld('SH_KCLDBASE', horiz_only, 'A', '1', &
               'Shallow conv. cloudbase level index' )

   call addfld('DP_MFUP_MAX', horiz_only, 'A', 'kg/m2', &
               'Deep conv. column-max updraft mass flux' )
   call addfld('DP_WCLDBASE', horiz_only, 'A', 'm/s', &
               'Deep conv. cloudbase vertical velocity' )
   call addfld('DP_KCLDBASE', horiz_only, 'A', '1', &
               'Deep conv. cloudbase level index' )

   ! output wet deposition fields to history
   !    I = in-cloud removal;     E = precip-evap resuspension
   !    C = convective (total);   D = deep convective
   ! note that the precip-evap resuspension includes that resulting from
   !    below-cloud removal, calculated in mz_aero_wet_intr
   if (convproc_do_aer .and. apply_convproc_tend_to_ptend ) then

      call rad_cnst_get_info( 0, nbins=nbins)
      allocate( nspec(nbins) )

      do m = 1, nbins
         call rad_cnst_get_info_by_bin(0, m, nspec=nspec(m))
      end do
      ! add plus one to include number, total mmr and nspec
      nspec_max = maxval(nspec) + 2

      ncnst_tot = 0
      do m = 1, nbins
        ncnst_tot = ncnst_tot + nspec(m) + 2
      end do
      ncnst_extd = 2*ncnst_tot
      pcnst_extd = 2*ncnst_tot

      allocate( &
       bin_idx(nbins,nspec_max),      &
       bin_cnst_lq(nbins,nspec_max), &
       bin_cnst_idx(nbins,nspec_max), &
       ibl(ncnst_tot),              &
       imx_bl(nbins),              &
       imx_mmr_bl(nbins),              &
       imx_num_bl(nbins),              &
       fieldname_cw(ncnst_tot), &
       fieldname(ncnst_tot)  )

      ii = 0
      do m = 1, nbins
        do l = 1, nspec(m) + 2    ! do through nspec plus mmr and number
            ii = ii + 1
            bin_idx(m,l) = ii

            if (l <= nspec(m) ) then   ! species
               call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=fieldname(ii), spec_name_cw=fieldname_cw(ii))
            else if (l == nspec(m) + 1) then   ! bin mmr
               call rad_cnst_get_info_by_bin(0, m, mmr_name=fieldname(ii), mmr_name_cw=fieldname_cw(ii))
            else if (l == nspec(m) + 2) then   ! bin num dens
               call rad_cnst_get_info_by_bin(0, m, num_name=fieldname(ii), num_name_cw=fieldname_cw(ii))
            end if

            call addfld (trim(fieldname(ii))//'SFSEC', &
                 horiz_only,  'A','kg/m2/s','Wet deposition flux (precip evap, convective) at surface')
            if (history_aerosol) then
               call add_default(trim(fieldname(ii))//'SFSEC', 1, ' ')
            end if
            if ( deepconv_wetdep_history ) then
               call addfld (trim(fieldname(ii))//'SFSID', &
                    horiz_only,  'A','kg/m2/s','Wet deposition flux (incloud, deep convective) at surface')
               call addfld (trim(fieldname(ii))//'SFSED', &
                    horiz_only,  'A','kg/m2/s','Wet deposition flux (precip evap, deep convective) at surface')
               if (history_aerosol) then
                  call add_default(trim(fieldname(ii))//'SFSID', 1, ' ')
                  call add_default(trim(fieldname(ii))//'SFSED', 1, ' ')
               end if
            end if
            call cnst_get_ind(fieldname(ii), idxtmp, abort=.false.)
            if (idxtmp>0) then
               bin_cnst_lq(m,l) = .true.
               bin_cnst_idx(m,l) = idxtmp
               lq(idxtmp) = .true.
            else
               bin_cnst_lq(m,l) = .false.
               bin_cnst_idx(m,l) = 0
            end if
        end do
      end do

      imx = 0
      imx_mmr = 0
      imx_num = 0
      ipr = 0
      ipr_mmr = 0
      ipr_num = 0
      do m = 1, nbins
         call rad_cnst_get_info_by_bin(0, m, bin_name=bin_name)
         bin_name_l = ' '
         if (m.lt.nbins) then
            call rad_cnst_get_info_by_bin(0, m+1, bin_name=bin_name_l)
         end if
         do l = 1, nspec(m) + 2    ! do through nspec plus mmr and number
            ii = bin_idx(m,l)
            ibl(ii) = ii

            ! derive index  array for larger bin, for evaporation into larger bi
            if (l.le.nspec(m)) then !  only for elements
               call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
            else
               spectype = 'other'
            end if
            ! identification is required for pure and mixed aerosols, mixed aeroosols are moved to larger bin, pure aerosols are moved to mixed sulfate

            if (index(bin_name,'MXAER')>0 .and. index(bin_name_l,'MXAER')>0) then   ! for mixed aerosols
               ! find larger bin
               ibl(ii) =  bin_idx(m+1,l)
               ! define mixed aerosol sulfate index to be used for pure sulfate only
               if (trim(spectype) == 'sulfate') then
                  imx = imx + 1
                  imx_bl(imx) = ibl(ii)
                  !write(*,*) "MXSULF: m, l, imx, imx_bl(imx), fieldname(ii), spectype ",m, l, imx, imx_bl(imx), fieldname(ii),spectype
               end if
               if (l.eq.nspec(m)+1) then !  only for mmr
                  imx_mmr = imx_mmr + 1
                  ibl(ii) =  bin_idx(m+1,l)
                  imx_mmr_bl(imx_mmr) = ibl(ii)
               end if
               if (l.eq.nspec(m)+2) then !  only for num
                  imx_num = imx_num + 1
                  ibl(ii) =  bin_idx(m+1,l)
                  imx_num_bl(imx_num) = ibl(ii)
               end if
            end if ! MXAER

            if (index(bin_name,'PRSUL')>0 .and. index(bin_name_l,'PRSUL')>0) then  ! assuming  PRSULF and  MXSULF have the same number of bins
               if (trim(spectype) == 'sulfate') then
                  ipr = ipr +1
                  ibl(ii) = imx_bl(ipr)
               end if
               if (l.eq.nspec(m)+1) then !  only for mmr reset counter to only go from 1-20 bins
                  ipr_mmr = ipr_mmr +1
                  ibl(ii) = imx_mmr_bl(ipr_mmr)
               end if
               if (l.eq.nspec(m)+2) then !  only for num reset counter to only go from 1-20 bins
                  ipr_num= ipr_num +1
                  ibl(ii) = imx_num_bl(ipr_num)
               end if
            end if
            if (ibl(ii).eq.0) then
               ibl(ii) = ii
            end if
         end do
      end do

   end if

   if ( history_aerosol .and. convproc_do_aer ) then
      call add_default( 'SH_MFUP_MAX', 1, ' ' )
      call add_default( 'SH_WCLDBASE', 1, ' ' )
      call add_default( 'SH_KCLDBASE', 1, ' ' )
      call add_default( 'DP_MFUP_MAX', 1, ' ' )
      call add_default( 'DP_WCLDBASE', 1, ' ' )
      call add_default( 'DP_KCLDBASE', 1, ' ' )
   end if

   fracis_idx      = pbuf_get_index('FRACIS')

   rprddp_idx      = pbuf_get_index('RPRDDP')
   rprdsh_idx      = pbuf_get_index('RPRDSH')
   nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
   nevapr_shcu_idx = pbuf_get_index('NEVAPR_SHCU')

   icwmrdp_idx     = pbuf_get_index('ICWMRDP')
   icwmrsh_idx     = pbuf_get_index('ICWMRSH')
   dp_frac_idx     = pbuf_get_index('DP_FRAC')
   sh_frac_idx     = pbuf_get_index('SH_FRAC')

   zm_mu_idx       = pbuf_get_index('ZM_MU')
   zm_eu_idx       = pbuf_get_index('ZM_EU')
   zm_du_idx       = pbuf_get_index('ZM_DU')
   zm_md_idx       = pbuf_get_index('ZM_MD')
   zm_ed_idx       = pbuf_get_index('ZM_ED')
   zm_dp_idx       = pbuf_get_index('ZM_DP')
   zm_dsubcld_idx  = pbuf_get_index('ZM_DSUBCLD')
   zm_jt_idx       = pbuf_get_index('ZM_JT')
   zm_maxg_idx     = pbuf_get_index('ZM_MAXG')
   zm_ideep_idx    = pbuf_get_index('ZM_IDEEP')

   cmfmc_sh_idx    = pbuf_get_index('CMFMC_SH')
   sh_e_ed_ratio_idx = pbuf_get_index('SH_E_ED_RATIO', istat)

   if (masterproc ) then

      write(iulog,'(a,l12)')     'ma_convproc_init - convproc_do_aer               = ', &
         convproc_do_aer
      write(iulog,'(a,l12)')     'ma_convproc_init - use_cwaer_for_activate_maxsat = ', &
         use_cwaer_for_activate_maxsat
      write(iulog,'(a,l12)')     'ma_convproc_init - apply_convproc_tend_to_ptend  = ', &
         apply_convproc_tend_to_ptend
      write(iulog,'(a,i12)')     'ma_convproc_init - convproc_method_activate      = ', &
         convproc_method_activate
      write(iulog,'(a,i12)')     'ma_convproc_init - method1_activate_nlayers      = ', &
         method1_activate_nlayers
      write(iulog,'(a,1pe12.4)') 'ma_convproc_init - method2_activate_smaxmax      = ', &
         method2_activate_smaxmax
      write(iulog,'(a,i12)')     'ma_convproc_init - method_reduce_actfrac         = ', &
         method_reduce_actfrac
      write(iulog,'(a,1pe12.4)') 'ma_convproc_init - factor_reduce_actfrac         = ', &
         factor_reduce_actfrac

      npass_calc_updraft = 1
      if ( (method_reduce_actfrac == 2)      .and. &
         (factor_reduce_actfrac >= 0.0_r8) .and. &
         (factor_reduce_actfrac <= 1.0_r8) ) npass_calc_updraft = 2
      write(iulog,'(a,i12)')     'ma_convproc_init - npass_calc_updraft            = ', &
         npass_calc_updraft

   end if

   aero_props_obj => carma_aerosol_properties()

end subroutine ma_convproc_init

!=========================================================================================

subroutine ma_convproc_intr( state, ptend, pbuf, ztodt,             &
                           nsrflx_mzaer2cnvpr, qsrflx_mzaer2cnvpr,  &
                           aerdepwetis, dcondt_resusp3d )
!-----------------------------------------------------------------------
!
! Convective cloud processing (transport, activation/resuspension,
!    wet removal) of aerosols and trace gases.
!    (Currently no aqueous chemistry and no trace-gas wet removal)
! Does aerosols    when convproc_do_aer is .true.
!
! Does deep and shallow convection
! Uses mass fluxes, cloud water, precip production from the
!    convective cloud routines
!
! Author: R. Easter
!
!-----------------------------------------------------------------------


   ! Arguments
   type(physics_state),       intent(in )   :: state      ! Physics state variables
   type(physics_ptend),       intent(inout) :: ptend      ! %lq set in aero_model_wetdep
   type(physics_buffer_desc), pointer       :: pbuf(:)
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)

   type(ptr2d_t), allocatable :: raer(:)     ! aerosol mass, number mixing ratios


   integer,  intent(in)    :: nsrflx_mzaer2cnvpr
   real(r8), intent(in)    :: qsrflx_mzaer2cnvpr(pcols,ncnst_tot,nsrflx_mzaer2cnvpr)
   real(r8), intent(inout) :: aerdepwetis(pcols,ncnst_tot)  ! aerosol wet deposition (interstitial)
   real(r8), intent(inout) :: dcondt_resusp3d(pcnst_extd,pcols,pver)

   ! Local variables
   integer, parameter :: nsrflx = 5        ! last dimension of qsrflx
   integer  :: l, ll, lchnk
   integer  :: lpr
   integer  :: n, ncol

   real(r8) :: dqdt(pcols,pver,ncnst_tot)
   real(r8) :: dt
   real(r8) :: qa(pcols,pver,ncnst_tot), qb(pcols,pver,ncnst_tot)
   real(r8) :: qsrflx(pcols,ncnst_tot,nsrflx)
   real(r8) :: sflxic(pcols,ncnst_tot)
   real(r8) :: sflxid(pcols,ncnst_tot)
   real(r8) :: sflxec(pcols,ncnst_tot)
   real(r8) :: sflxed(pcols,ncnst_tot)

   logical  :: dotend(ncnst_tot)  ! will be all CARMA aerosols so will be set to .true. for all aerosols in ma_convproc_sh_intr
   !-------------------------------------------------------------------------------------------------

   ! Initialize
   lchnk = state%lchnk
   ncol  = state%ncol
   dt = ztodt

   hund_ovr_g = 100.0_r8/gravit
   !  used with zm_conv mass fluxes and delta-p
   !     for mu = [mbar/s],   mu*hund_ovr_g = [kg/m2/s]
   !     for dp = [mbar] and q = [kg/kg],   q*dp*hund_ovr_g = [kg/m2]

   sflxic(:,:) = 0.0_r8
   sflxid(:,:) = 0.0_r8
   sflxec(:,:) = 0.0_r8
   sflxed(:,:) = 0.0_r8

   allocate(raer(ncnst_tot))

   do n = 1, nbins      ! main loop over aerosol bins
      call rad_cnst_get_info_by_bin(0, n, nspec=nspec(n))
      do ll = 1, nspec(n) + 2
         l = bin_idx(n, ll)
         lpr = bin_cnst_idx(n,ll)

         sflxec(1:ncol,l) = qsrflx_mzaer2cnvpr(1:ncol,l,1)
         sflxed(1:ncol,l) = qsrflx_mzaer2cnvpr(1:ncol,l,2)

         if (ll <= nspec(n)) then
             call rad_cnst_get_bin_mmr_by_idx(0, n, ll, 'a', state, pbuf, raer(l)%fld)
         end if
         if (ll == nspec(n)+1) then ! mmr
              call rad_cnst_get_bin_mmr(0, n, 'a', state, pbuf, raer(l)%fld)
         end if
         if (ll == nspec(n)+2) then ! number
              call rad_cnst_get_bin_num(0, n, 'a', state, pbuf, raer(l)%fld)
         end if
        ! calc new raer(l) if from state, add ptend
         if (bin_cnst_lq(n,ll)) then ! advected species, non advective species have been updated
            qa(:ncol,:,l) = raer(l)%fld(1:ncol,:) + dt*ptend%q(1:ncol,:,lpr)
            qb(:ncol,:,l) = max( 0.0_r8, qa(1:ncol,:,l) )
         else
            qb(:ncol,:,l) = raer(l)%fld(:ncol,:)
         end if

      end do

      ! prepare for deep conv processing
      !st not done because calcaersize and mz_aero_wet_intr was not applied for CARMA

      !st do l = 1, pcnst
      !st    if ( ptend%lq(l) ) then
      !st       ! calc new q (after calcaersize and mz_aero_wet_intr)
      !st       qa(1:ncol,:,l) = state%q(1:ncol,:,l) + dt*ptend%q(1:ncol,:,l)
      !st       qb(1:ncol,:,l) = max( 0.0_r8, qa(1:ncol,:,l) )
      !st    else
      !st       ! use old q
      !st       qb(1:ncol,:,l) = state%q(1:ncol,:,l)
      !st    end if
      !st end do

   end do    !nbins
   dqdt(:,:,:) = 0.0_r8
   qsrflx(:,:,:) = 0.0_r8

   if (convproc_do_aer) then

      ! do deep conv processing
      if (convproc_do_deep) then
         call ma_convproc_dp_intr(                    &
            state, pbuf, dt,                          &
            qb, dqdt, dotend, nsrflx, qsrflx, dcondt_resusp3d )


         ! apply deep conv processing tendency and prepare for shallow conv processing
         !st do l = 1, pcnst  ! mm in aero_model
         do n = 1, nbins      ! main loop over aerosol bins
            do ll = 1, nspec(n) + 2
               l = bin_idx(n, ll)
               lpr = bin_cnst_idx(n,ll)

               if ( .not. dotend(l) ) cycle

               ! calc new q (after ma_convproc_dp_intr)
               qa(1:ncol,:,l) = qb(1:ncol,:,l) + dt*dqdt(1:ncol,:,l)
               qb(1:ncol,:,l) = max( 0.0_r8, qa(1:ncol,:,l) )

               if ( apply_convproc_tend_to_ptend ) then
                  ! add dqdt onto ptend%q and set ptend%lq   and set pbuf variable if not transported
                  if (bin_cnst_lq(n,ll)) then ! adveced species
                     ptend%q(1:ncol,:,lpr) = ptend%q(1:ncol,:,lpr) + dqdt(1:ncol,:,l)
                  else
                     raer(l)%fld(1:ncol,:) = max( 0.0_r8, raer(l)%fld(1:ncol,:) + dqdt(1:ncol,:,l) * dt )
                  end if

               end if

               !st done for all CARMA aerosols
               ! these used for history file wetdep diagnostics
               sflxic(1:ncol,l) = sflxic(1:ncol,l) + qsrflx(1:ncol,l,4)
               sflxid(1:ncol,l) = sflxid(1:ncol,l) + qsrflx(1:ncol,l,4)
               sflxec(1:ncol,l) = sflxec(1:ncol,l) + qsrflx(1:ncol,l,5)
               sflxed(1:ncol,l) = sflxed(1:ncol,l) + qsrflx(1:ncol,l,5)
               !st end if

               !st if (cnst_species_class(l) == cnst_spec_class_aerosol) then
               ! this used for surface coupling
               aerdepwetis(1:ncol,l) = aerdepwetis(1:ncol,l) &
                    + qsrflx(1:ncol,l,4) + qsrflx(1:ncol,l,5)
            end do
         end do
      end if

      dqdt(:,:,:) = 0.0_r8
      qsrflx(:,:,:) = 0.0_r8

   end if ! (convproc_do_aer) then

   if (convproc_do_aer .and. apply_convproc_tend_to_ptend ) then
      do n = 1, nbins      ! main loop over aerosol bins
         do ll = 1, nspec(n) + 2
            l = bin_idx(n, ll)

            call outfld( trim(fieldname(l))//'SFWET', aerdepwetis(:,l), pcols, lchnk )
            call outfld( trim(fieldname(l))//'SFSIC', sflxic(:,l), pcols, lchnk )
            call outfld( trim(fieldname(l))//'SFSEC', sflxec(:,l), pcols, lchnk )

            if ( deepconv_wetdep_history ) then
               call outfld( trim(fieldname(l))//'SFSID', sflxid(:,l), pcols, lchnk )
               call outfld( trim(fieldname(l))//'SFSED', sflxed(:,l), pcols, lchnk )
            end if
         end do
      end do
   end if

end subroutine ma_convproc_intr

!=========================================================================================

subroutine ma_convproc_dp_intr(                &
     state, pbuf, dt,                          &
     q, dqdt, dotend, nsrflx, qsrflx,  dcondt_resusp3d)
!-----------------------------------------------------------------------
!
! Convective cloud processing (transport, activation/resuspension,
!    wet removal) of aerosols and trace gases.
!    (Currently no aqueous chemistry and no trace-gas wet removal)
! Does aerosols    when convproc_do_aer is .true.
!
! This routine does deep convection
! Uses mass fluxes, cloud water, precip production from the
!    convective cloud routines
!
! Author: R. Easter
!
!-----------------------------------------------------------------------


   ! Arguments
   type(physics_state),       intent(in ) :: state          ! Physics state variables
   type(physics_buffer_desc), pointer     :: pbuf(:)

   real(r8), intent(in) :: dt                         ! delta t (model time increment)

   real(r8), intent(in)    :: q(pcols,pver,ncnst_tot)
   real(r8), intent(inout) :: dqdt(pcols,pver,ncnst_tot)
   logical,  intent(out)   :: dotend(ncnst_tot)
   integer,  intent(in)    :: nsrflx
   real(r8), intent(inout) :: qsrflx(pcols,ncnst_tot,nsrflx)
   real(r8), intent(inout) :: dcondt_resusp3d(pcnst_extd,pcols,pver)

   integer :: i, n, ll
   integer :: l, lchnk
   integer :: nstep
   integer :: ncol

   real(r8) :: dpdry(pcols,pver)     ! layer delta-p-dry (mb)
   real(r8) :: fracice(pcols,pver)   ! Ice fraction of cloud droplets
   real(r8) :: qaa(pcols,pver,ncnst_tot)
   real(r8) :: xx_mfup_max(pcols), xx_wcldbase(pcols), xx_kcldbase(pcols)

   ! physics buffer fields
   real(r8), pointer :: fracis(:,:,:)  ! fraction of transported species that are insoluble
   real(r8), pointer :: rprddp(:,:)    ! Deep conv precip production (kg/kg/s - grid avg)
   real(r8), pointer :: evapcdp(:,:)   ! Deep conv precip evaporation (kg/kg/s - grid avg)
   real(r8), pointer :: icwmrdp(:,:)   ! Deep conv cloud condensate (kg/kg - in cloud)
   real(r8), pointer :: dp_frac(:,:)   ! Deep conv cloud frac (0-1)
   ! mu, md, ..., ideep, lengath are all deep conv variables
   real(r8), pointer :: mu(:,:)        ! Updraft mass flux (positive) (pcols,pver)
   real(r8), pointer :: md(:,:)        ! Downdraft mass flux (negative) (pcols,pver)
   real(r8), pointer :: du(:,:)        ! Mass detrain rate from updraft (pcols,pver)
   real(r8), pointer :: eu(:,:)        ! Mass entrain rate into updraft (pcols,pver)
   real(r8), pointer :: ed(:,:)        ! Mass entrain rate into downdraft (pcols,pver)
                                       ! eu, ed, du are "d(massflux)/dp" and are all positive
   real(r8), pointer :: dp(:,:)        ! Delta pressure between interfaces (pcols,pver)
   real(r8), pointer :: dsubcld(:)     ! Delta pressure from cloud base to sfc (pcols)

   integer,  pointer :: jt(:)          ! Index of cloud top for each column (pcols)
   integer,  pointer :: maxg(:)        ! Index of cloud top for each column (pcols)
   integer,  pointer :: ideep(:)       ! Gathering array (pcols)
   integer           :: lengath        ! Gathered min lon indices over which to operate


   ! Initialize

   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ! Associate pointers with physics buffer fields
   call pbuf_get_field(pbuf, fracis_idx,      fracis)
   call pbuf_get_field(pbuf, rprddp_idx,      rprddp)
   call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp)
   call pbuf_get_field(pbuf, icwmrdp_idx,     icwmrdp)
   call pbuf_get_field(pbuf, dp_frac_idx,     dp_frac)
   call pbuf_get_field(pbuf, fracis_idx,      fracis)
   call pbuf_get_field(pbuf, zm_mu_idx,       mu)
   call pbuf_get_field(pbuf, zm_eu_idx,       eu)
   call pbuf_get_field(pbuf, zm_du_idx,       du)
   call pbuf_get_field(pbuf, zm_md_idx,       md)
   call pbuf_get_field(pbuf, zm_ed_idx,       ed)
   call pbuf_get_field(pbuf, zm_dp_idx,       dp)
   call pbuf_get_field(pbuf, zm_dsubcld_idx,  dsubcld)
   call pbuf_get_field(pbuf, zm_jt_idx,       jt)
   call pbuf_get_field(pbuf, zm_maxg_idx,     maxg)
   call pbuf_get_field(pbuf, zm_ideep_idx,    ideep)

   lengath = count(ideep > 0)

   fracice(:,:) = 0.0_r8

   ! initialize dpdry (units=mb), which is used for tracers of dry mixing ratio type
   dpdry = 0._r8
   do i = 1, lengath
      dpdry(i,:) = state%pdeldry(ideep(i),:)/100._r8
   end do

   qaa = q

   ! turn on/off calculations for aerosols and trace gases
   ! st CARMA all aerosols are turned on
   dotend(:) = .true.

   call ma_convproc_tend(                                            &
                     'deep',                                         &
                     state, pbuf,                                    &
                     lchnk,      pcnst,      nstep,      dt,         &
                     state%t,    state%pmid, state%pdel, qaa,        &
                     mu,         md,         du,         eu,         &
                     ed,         dp,         dpdry,      jt,         &
                     maxg,       ideep,      1,          lengath,    &
                     dp_frac,    icwmrdp,    rprddp,     evapcdp,    &
                     fracice,                                        &
                     dqdt,       dotend,     nsrflx,     qsrflx,     &
                     xx_mfup_max, xx_wcldbase, xx_kcldbase,          &
                     dcondt_resusp3d  )

   call outfld( 'DP_MFUP_MAX', xx_mfup_max, pcols, lchnk )
   call outfld( 'DP_WCLDBASE', xx_wcldbase, pcols, lchnk )
   call outfld( 'DP_KCLDBASE', xx_kcldbase, pcols, lchnk )

end subroutine ma_convproc_dp_intr

!=========================================================================================

subroutine ma_convproc_tend(                                           &
                     convtype,                                       &
                     state, pbuf,                                    &
                     lchnk,      ncnst,      nstep,      dt,         &
                     t,          pmid,       pdel,       q,          &
                     mu,         md,         du,         eu,         &
                     ed,         dp,         dpdry,      jt,         &
                     mx,         ideep,      il1g,       il2g,       &
                     cldfrac,    icwmr,      rprd,       evapc,      &
                     fracice,                                        &
                     dqdt,       doconvproc, nsrflx,     qsrflx,     &
                     xx_mfup_max, xx_wcldbase, xx_kcldbase,          &
                     dcondt_resusp3d )

!-----------------------------------------------------------------------
!
! Purpose:
! Convective transport of trace species.
! The trace species need not be conservative, and source/sink terms for
!    activation, resuspension, aqueous chemistry and gas uptake, and
!    wet removal are all applied.
! Currently this works with the ZM deep convection, but we should be able
!    to adapt it for both Hack and McCaa shallow convection
!
!
! Compare to subr convproc which does conservative trace species.
!
! A distinction between "moist" and "dry" mixing ratios is not currently made.
! (P. Rasch comment:  Note that we are still assuming that the tracers are
!  in a moist mixing ratio this will change soon)

!
! Method:
! Computes tracer mixing ratios in updraft and downdraft "cells" in a
! Lagrangian manner, with source/sinks applied in the updraft other.
! Then computes grid-cell-mean tendencies by considering
!    updraft and downdraft fluxes across layer boundaries
!    environment subsidence/lifting fluxes across layer boundaries
!    sources and sinks in the updraft
!    resuspension of activated species in the grid-cell as a whole
!
! Note1:  A better estimate or calculation of either the updraft velocity
!         or fractional area is needed.
! Note2:  If updraft area is a small fraction of over cloud area,
!         then aqueous chemistry is underestimated.  These are both
!         research areas.
!
! Authors: O. Seland and R. Easter, based on convtran by P. Rasch
!
!-----------------------------------------------------------------------

   implicit none

!-----------------------------------------------------------------------
!
! Input arguments
!
   type(physics_state), target, intent(in)    :: state

   type(physics_buffer_desc), pointer :: pbuf(:)

   character(len=*), intent(in) :: convtype  ! identifies the type of
                                             ! convection ("deep", "shcu")
   integer,  intent(in) :: lchnk             ! chunk identifier
   integer,  intent(in) :: ncnst             ! number of aerosol tracers to transport
   integer,  intent(in) :: nstep             ! Time step index
   real(r8), intent(in) :: dt                ! Model timestep
   real(r8), intent(in) :: t(pcols,pver)     ! Temperature
   real(r8), intent(in) :: pmid(pcols,pver)  ! Pressure at model levels
   real(r8), intent(in) :: pdel(pcols,pver)  ! Pressure thickness of levels
   real(r8), intent(in) :: q(pcols,pver,ncnst_tot) ! Tracer array including moisture

   real(r8), intent(in) :: mu(pcols,pver)    ! Updraft mass flux (positive)
   real(r8), intent(in) :: md(pcols,pver)    ! Downdraft mass flux (negative)
   real(r8), intent(in) :: du(pcols,pver)    ! Mass detrain rate from updraft
   real(r8), intent(in) :: eu(pcols,pver)    ! Mass entrain rate into updraft
   real(r8), intent(in) :: ed(pcols,pver)    ! Mass entrain rate into downdraft
! *** note1 - mu, md, eu, ed, du, dp, dpdry are GATHERED ARRAYS ***
! *** note2 - mu and md units are (mb/s), which is used in the zm_conv code
!           - eventually these should be changed to (kg/m2/s)
! *** note3 - eu, ed, du are "d(massflux)/dp" (with dp units = mb), and are all >= 0

   real(r8), intent(in) :: dp(pcols,pver)    ! Delta pressure between interfaces (mb)
   real(r8), intent(in) :: dpdry(pcols,pver) ! Delta dry-pressure (mb)
!  real(r8), intent(in) :: dsubcld(pcols)    ! Delta pressure from cloud base to sfc
   integer,  intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer,  intent(in) :: mx(pcols)         ! Index of cloud top for each column
   integer,  intent(in) :: ideep(pcols)      ! Gathering array indices
   integer,  intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer,  intent(in) :: il2g              ! Gathered max lon indices over which to operate
! *** note4 -- for il1g <= i <= il2g,  icol = ideep(i) is the "normal" chunk column index

   real(r8), intent(in) :: cldfrac(pcols,pver)  ! Convective cloud fractional area
   real(r8), intent(in) :: icwmr(pcols,pver)    ! Convective cloud water from zhang
   real(r8), intent(in) :: rprd(pcols,pver)     ! Convective precipitation formation rate
   real(r8), intent(in) :: evapc(pcols,pver)    ! Convective precipitation evaporation rate
   real(r8), intent(in) :: fracice(pcols,pver)  ! Ice fraction of cloud droplets

   real(r8), intent(out):: dqdt(pcols,pver,ncnst_tot)  ! Tracer tendency array
   logical,  intent(in) :: doconvproc(ncnst_tot) ! flag for doing convective transport
   integer,  intent(in) :: nsrflx            ! last dimension of qsrflx
   real(r8), intent(out):: qsrflx(pcols,ncnst_tot,nsrflx)
                              ! process-specific column tracer tendencies
                              ! (1=activation,  2=resuspension, 3=aqueous rxn,
                              !  4=wet removal, 5=renaming)
   real(r8), intent(out) :: xx_mfup_max(pcols)
   real(r8), intent(out) :: xx_wcldbase(pcols)
   real(r8), intent(out) :: xx_kcldbase(pcols)
   real(r8), intent(inout) :: dcondt_resusp3d(pcnst_extd,pcols,pver)

!--------------------------Local Variables------------------------------

! cloudborne aerosol, so the arrays are dimensioned with pcnst_extd = ncnst_tot*2

   integer :: i, icol         ! Work index
   integer :: iconvtype       ! 1=deep, 2=uw shallow
   integer :: iflux_method    ! 1=as in convtran (deep), 2=simpler
   integer :: ipass_calc_updraft
   integer :: itmpa, itmpb    ! Work variable
   integer :: j, jtsub        ! Work index
   integer :: k               ! Work index
   integer :: kactcnt         ! Counter for no. of levels having activation
   integer :: kactcntb        ! Counter for activation diagnostic output
   integer :: kactfirst       ! Lowest layer with activation (= cloudbase)
   integer :: kbot            ! Cloud-flux bottom layer for current i (=mx(i))
   integer :: kbot_prevap     ! Lowest layer for doing resuspension from evaporating precip
   integer :: ktop            ! Cloud-flux top    layer for current i (=jt(i))
                              ! Layers between kbot,ktop have mass fluxes
                              !    but not all have cloud water, because the
                              !    updraft starts below the cloud base
   integer :: km1, km1x       ! Work index
   integer :: kp1, kp1x       ! Work index
   integer :: l, ll, la, lc   ! Work index
   integer :: m, n, ii           ! Work index
   integer :: merr            ! number of errors (i.e., failed diagnostics)
                              ! for current column
   integer :: nerr            ! number of errors for entire run
   integer :: nerrmax         ! maximum number of errors to report
   integer :: ncnst_extd      !initialized above
   integer :: npass_calc_updraft
   integer :: ntsub           !

   logical  do_act_this_lev             ! flag for doing activation at current level
   logical  doconvproc_extd(pcnst_extd) ! flag for doing convective transport

   real(r8) aqfrac(pcnst_extd)       ! aqueous fraction of constituent in updraft
   real(r8) cldfrac_i(pver)          ! cldfrac at current i (with adjustments)

   real(r8) chat(pcnst_extd,pverp)   ! mix ratio in env at interfaces
   real(r8) cond(pcnst_extd,pverp)   ! mix ratio in downdraft at interfaces
   real(r8) const(pcnst_extd,pver)   ! gathered tracer array
   real(r8) conu(pcnst_extd,pverp)   ! mix ratio in updraft at interfaces

   real(r8) dcondt(pcnst_extd,pver)  ! grid-average TMR tendency for current column
   real(r8) dcondt_prevap(pcnst_extd,pver) ! portion of dcondt from precip evaporation
   real(r8) dcondt_resusp(pcnst_extd,pver) ! portion of dcondt from resuspension

   real(r8) dcondt_wetdep(pcnst_extd,pver) ! portion of dcondt from wet deposition
   real(r8) dconudt_activa(pcnst_extd,pverp) ! d(conu)/dt by activation
   real(r8) dconudt_aqchem(pcnst_extd,pverp) ! d(conu)/dt by aqueous chem
   real(r8) dconudt_wetdep(pcnst_extd,pverp) ! d(conu)/dt by wet removal

   real(r8) maxflux(pcnst_extd)      ! maximum (over layers) of fluxin and fluxout
   real(r8) maxflux2(pcnst_extd)     ! ditto but computed using method-2 fluxes
   real(r8) maxprevap(pcnst_extd)    ! maximum (over layers) of dcondt_prevap*dp
   real(r8) maxresusp(pcnst_extd)    ! maximum (over layers) of dcondt_resusp*dp
   real(r8) maxsrce(pcnst_extd)      ! maximum (over layers) of netsrce

   real(r8) sumflux(pcnst_extd)      ! sum (over layers) of netflux
   real(r8) sumflux2(pcnst_extd)     ! ditto but computed using method-2 fluxes
   real(r8) sumsrce(pcnst_extd)      ! sum (over layers) of dp*netsrce
   real(r8) sumchng(pcnst_extd)      ! sum (over layers) of dp*dcondt
   real(r8) sumchng3(pcnst_extd)     ! ditto but after call to resusp_conv
   real(r8) sumprevap(pcnst_extd)    ! sum (over layers) of dp*dcondt_prevap
   real(r8) sumwetdep(pcnst_extd)    ! sum (over layers) of dp*dconudt_wetdep

   real(r8) cabv                 ! mix ratio of constituent above
   real(r8) cbel                 ! mix ratio of constituent below
   real(r8) cdifr                ! normalized diff between cabv and cbel
   real(r8) cdt(pver)            ! (in-updraft first order wet removal rate) * dt
   real(r8) cdt2(pver)            ! (in-updraft first order wet removal rate) * dt
   real(r8) clw_cut              ! threshold clw value for doing updraft
                                 ! transformation and removal
   real(r8) courantmax           ! maximum courant no.
   real(r8) dddp(pver)           ! dd(i,k)*dp(i,k) at current i
   real(r8) dp_i(pver)           ! dp(i,k) at current i
   real(r8) dt_u(pver)           ! lagrangian transport time in the updraft
   real(r8) dudp(pver)           ! du(i,k)*dp(i,k) at current i
   real(r8) dqdt_i(pver,ncnst_tot)   ! dqdt(i,k,m) at current i
   real(r8) dtsub                ! dt/ntsub
   real(r8) dz                   ! working layer thickness (m)
   real(r8) eddp(pver)           ! ed(i,k)*dp(i,k) at current i
   real(r8) eudp(pver)           ! eu(i,k)*dp(i,k) at current i
   real(r8) expcdtm1             ! a work variable
   real(r8) expcdtm2             ! a work variable
   real(r8) fa_u(pver)           ! fractional area of in the updraft
   real(r8) fa_u_dp              ! current fa_u(k)*dp_i(k)
   real(r8) f_ent                ! fraction of the "before-detrainment" updraft
                                 ! massflux at k/k-1 interface resulting from
                                 ! entrainment of level k air
   real(r8) fluxin               ! a work variable
   real(r8) fluxout              ! a work variable
   real(r8) maxc                 ! a work variable
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) minc                 ! a work variable
   real(r8) md_m_eddp            ! a work variable
   real(r8) md_i(pverp)          ! md(i,k) at current i (note pverp dimension)
   real(r8) md_x(pverp)          ! md(i,k) at current i (note pverp dimension)
   real(r8) mu_i(pverp)          ! mu(i,k) at current i (note pverp dimension)
   real(r8) mu_x(pverp)          ! mu(i,k) at current i (note pverp dimension)
   ! md_i, md_x, mu_i, mu_x are all "dry" mass fluxes
   ! the mu_x/md_x are initially calculated from the incoming mu/md by applying dp/dpdry
   ! the mu_i/md_i are next calculated by applying the mbsth threshold
   real(r8) mu_p_eudp(pver)      ! = mu_i(kp1) + eudp(k)
   real(r8) netflux              ! a work variable
   real(r8) netsrce              ! a work variable
   real(r8) q_i(pver,ncnst_tot)      ! q(i,k,m) at current i
   real(r8) qsrflx_i(ncnst_tot,nsrflx) ! qsrflx(i,m,n) at current i
   real(r8) relerr_cut           ! relative error criterion for diagnostics
   real(r8) rhoair_i(pver)       ! air density at current i
   real(r8) small                ! a small number
   real(r8) tmpa, tmpb           ! work variables
   real(r8) tmpf                 ! work variables
   real(r8) tmpveca(pcnst_extd)  ! work variables
   real(r8) tmpmata(pcnst_extd,3) ! work variables
   real(r8) xinv_ntsub           ! 1.0/ntsub
   real(r8) wup(pver)            ! working updraft velocity (m/s)
   real(r8) zmagl(pver)          ! working height above surface (m)
   real(r8) zkm                  ! working height above surface (km)

   real(r8) conu2(pcols,pver,pcnst_extd)   !q(i,k,m) at current i
   real(r8) dcondt2(pcols,pver,pcnst_extd)   !q(i,k,m) at current i

   character(len=16) :: cnst_name_extd(pcnst_extd)
   !st character(len=32), allocatable :: fieldname(:)    ! names for interstitial output fields
   !st character(len=32), allocatable :: fieldname_cw(:)    ! names for cloud_borne output fields

   !Fractional area of ensemble mean updrafts in ZM scheme set to 0.01
   !Chosen to reproduce vertical vecocities in GATEIII GIGALES (Khairoutdinov etal 2009, JAMES)
   real(r8), parameter :: zm_areafrac = 0.01_r8
!-----------------------------------------------------------------------
!

!  if (nstep > 1) call endrun()

   if (convtype == 'deep') then
      iconvtype = 1
      iflux_method = 1
   else if (convtype == 'uwsh') then
      iconvtype = 2
      iflux_method = 2
   else
      call endrun( '*** ma_convproc_tend -- convtype is not |deep| or |uwsh|' )
   end if

   nerr = 0
   nerrmax = 99

   dcondt_resusp3d(:,:,:) = 0._r8

   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

   qsrflx(:,:,:) = 0.0_r8
   dqdt(:,:,:) = 0.0_r8
   xx_mfup_max(:) = 0.0_r8
   xx_wcldbase(:) = 0.0_r8
   xx_kcldbase(:) = 0.0_r8
   conu2(:,:,:) = 0.0_r8
   dcondt2(:,:,:) = 0.0_r8
   wup(:) = 0.0_r8

! set doconvproc_extd (extended array) values
! inititialize aqfrac to 1.0 for activated aerosol species, 0.0 otherwise
! st for CARMA all cloud borne aerosols are activate, and aqfrac is set to 1.0_r8 for those


   doconvproc_extd(:) = .false.
   doconvproc_extd(1:ncnst_tot) = doconvproc(1:ncnst_tot)
   aqfrac(:) = 0.0_r8

   do n = 1, nbins      ! main loop over aerosol bins
         do ll = 1, nspec(n) + 2
            l = bin_idx(n, ll)
            la = l
            lc = l + ncnst_tot
            cnst_name_extd(la) = fieldname(l)
            cnst_name_extd(lc) = fieldname_cw(l)
          end do
   end do

   !do l = 1, pcnst_extd
   !   if (l <= ncnst_tot) then
   !      cnst_name_extd(l) = fieldname(l)
   !   else
   !      cnst_name_extd(l) = fieldname_cw(l-ncnst_tot)
   !   end if
   !end do

   do n = 1, nbins
     do ll = 1, nspec(n) + 2
        l = bin_idx(n,ll)
        la = l
        lc = l + ncnst_tot
        if ( doconvproc(la) ) then
            doconvproc_extd(lc) = .true.
            aqfrac(lc) = 1.0_r8
            ! write(*,*) 'doconvproc_extd',la,lc
         end if
     end do ! ll
   end do ! n


! Loop ever each column that has convection
! *** i is index to gathered arrays; ideep(i) is index to "normal" chunk arrays
i_loop_main_aa: &
   do i = il1g, il2g
   icol = ideep(i)


   if ( (jt(i) <= 0) .and. (mx(i) <= 0) .and. (iconvtype /= 1) ) then
! shallow conv case with jt,mx <= 0, which means there is no shallow conv
! in this column -- skip this column
      cycle i_loop_main_aa

   else if ( (jt(i) < 1) .or. (mx(i) > pver) .or. (jt(i) > mx(i)) ) then
! invalid cloudtop and cloudbase indices -- skip this column
      write(*,9010) 'illegal jt, mx', convtype, lchnk, icol, i,    &
                                      jt(i), mx(i)
9010  format( '*** ma_convproc_tend error -- ', a, 5x, 'convtype = ', a /   &
              '*** lchnk, icol, il, jt, mx = ', 5(1x,i10) )
      cycle i_loop_main_aa

   else if (jt(i) == mx(i)) then
! cloudtop = cloudbase (1 layer cloud) -- skip this column
      write(*,9010) 'jt == mx', convtype, lchnk, icol, i, jt(i), mx(i)
      cycle i_loop_main_aa

   end if


!
! cloudtop and cloudbase indices are valid so proceed with calculations
!

! Load dp_i and cldfrac_i, and calc rhoair_i
      do k = 1, pver
         dp_i(k) = dpdry(i,k)
         cldfrac_i(k) = cldfrac(icol,k)
         rhoair_i(k) = pmid(icol,k)/(rair*t(icol,k))
      end do

! Calc dry mass fluxes
!    This is approximate because the updraft air is has different temp and qv than
!    the grid mean, but the whole convective parameterization is highly approximate
      mu_x(:) = 0.0_r8
      md_x(:) = 0.0_r8
! (eu-du) = d(mu)/dp -- integrate upwards, multiplying by dpdry
      do k = pver, 1, -1
         mu_x(k) = mu_x(k+1) + (eu(i,k)-du(i,k))*dp_i(k)
         xx_mfup_max(icol) = max( xx_mfup_max(icol), mu_x(k) )
      end do
! (ed) = d(md)/dp -- integrate downwards, multiplying by dpdry
      do k = 2, pver
         md_x(k) = md_x(k-1) - ed(i,k-1)*dp_i(k-1)
      end do

! Load mass fluxes over cloud layers
! (Note - use of arrays dimensioned k=1,pver+1 simplifies later coding)
! Zero out values below threshold
! Zero out values at "top of cloudtop", "base of cloudbase"
      ktop = jt(i)
      kbot = mx(i)
! usually the updraft ( & downdraft) start ( & end ) at kbot=pver, but sometimes kbot < pver
! transport, activation, resuspension, and wet removal only occur between kbot >= k >= ktop
! resuspension from evaporating precip can occur at k > kbot when kbot < pver
      kbot_prevap = pver
      mu_i(:) = 0.0_r8
      md_i(:) = 0.0_r8
      do k = ktop+1, kbot
         mu_i(k) = mu_x(k)
         if (mu_i(k) <= mbsth) mu_i(k) = 0.0_r8
         md_i(k) = md_x(k)
         if (md_i(k) >= -mbsth) md_i(k) = 0.0_r8
      end do
      mu_i(ktop) = 0.0_r8
      md_i(ktop) = 0.0_r8
      mu_i(kbot+1) = 0.0_r8
      md_i(kbot+1) = 0.0_r8

!  Compute updraft and downdraft "entrainment*dp" from eu and ed
!  Compute "detrainment*dp" from mass conservation
      eudp(:) = 0.0_r8
      dudp(:) = 0.0_r8
      eddp(:) = 0.0_r8
      dddp(:) = 0.0_r8
      courantmax = 0.0_r8
      do k = ktop, kbot
         if ((mu_i(k) > 0) .or. (mu_i(k+1) > 0)) then
            if (du(i,k) <= 0.0_r8) then
               eudp(k) = mu_i(k) - mu_i(k+1)
            else
               eudp(k) = max( eu(i,k)*dp_i(k), 0.0_r8 )
               dudp(k) = (mu_i(k+1) + eudp(k)) - mu_i(k)
               if (dudp(k) < 1.0e-12_r8*eudp(k)) then
                  eudp(k) = mu_i(k) - mu_i(k+1)
                  dudp(k) = 0.0_r8
               end if
            end if
         end if
         if ((md_i(k) < 0) .or. (md_i(k+1) < 0)) then
            eddp(k) = max( ed(i,k)*dp_i(k), 0.0_r8 )
            dddp(k) = (md_i(k+1) + eddp(k)) - md_i(k)
            if (dddp(k) < 1.0e-12_r8*eddp(k)) then
               eddp(k) = md_i(k) - md_i(k+1)
               dddp(k) = 0.0_r8
            end if
         end if
!        courantmax = max( courantmax, (eudp(k)+eddp(k))*dt/dp_i(k) )  ! old version - incorrect
         courantmax = max( courantmax, ( mu_i(k+1)+eudp(k)-md_i(k)+eddp(k) )*dt/dp_i(k) )
      end do ! k

! number of time substeps needed to maintain "courant number" <= 1
      ntsub = 1
      if (courantmax > (1.0_r8 + 1.0e-6_r8)) then
         ntsub = 1 + int( courantmax )
      end if
      xinv_ntsub = 1.0_r8/ntsub
      dtsub = dt*xinv_ntsub
      courantmax = courantmax*xinv_ntsub

! zmagl(k) = height above surface for middle of level k
      zmagl(pver) = 0.0_r8
      do k = pver, 1, -1
         if (k < pver) then
            zmagl(k) = zmagl(k+1) + 0.5_r8*dz
         end if
         dz = dp_i(k)*hund_ovr_g/rhoair_i(k)
         zmagl(k) = zmagl(k) + 0.5_r8*dz
      end do

!  load tracer mixing ratio array, which will be updated at the end of each jtsub interation
      q_i(1:pver,1:ncnst_tot) = q(icol,1:pver,1:ncnst_tot)
      conu2(icol,1:pver,1:ncnst_tot) = q(icol,1:pver,1:ncnst_tot)

!
!   when method_reduce_actfrac = 2, need to do the updraft calc twice
!   (1st to get non-adjusted activation amount, 2nd to apply reduction factor)
      npass_calc_updraft = 1
      if ( (method_reduce_actfrac == 2)      .and. &
           (factor_reduce_actfrac >= 0.0_r8) .and. &
           (factor_reduce_actfrac <= 1.0_r8) ) npass_calc_updraft = 2


jtsub_loop_main_aa: &
      do jtsub = 1, ntsub


ipass_calc_updraft_loop: &
      do ipass_calc_updraft = 1, npass_calc_updraft

      qsrflx_i(:,:) = 0.0_r8
      dqdt_i(:,:) = 0.0_r8

      const(:,:) = 0.0_r8 ! zero cloud-phase species
      chat(:,:) = 0.0_r8 ! zero cloud-phase species
      conu(:,:) = 0.0_r8
      cond(:,:) = 0.0_r8

      dcondt(:,:) = 0.0_r8
      dcondt_resusp(:,:) = 0.0_r8
      dcondt_wetdep(:,:) = 0.0_r8
      dcondt_prevap(:,:) = 0.0_r8
      dconudt_aqchem(:,:) = 0.0_r8
      dconudt_wetdep(:,:) = 0.0_r8
! only initialize the activation tendency on ipass=1
      if (ipass_calc_updraft == 1) dconudt_activa(:,:) = 0.0_r8

! initialize mixing ratio arrays (chat, const, conu, cond)
      do m = 1, ncnst_tot
      if ( doconvproc_extd(m) ) then

! Gather up the constituent
         do k = 1,pver
            const(m,k) = q_i(k,m)
         end do

! From now on work only with gathered data
! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            minc = min(const(m,km1),const(m,k))
            maxc = max(const(m,km1),const(m,k))
            if (minc < 0) then
               cdifr = 0._r8
            else
               cdifr = abs(const(m,k)-const(m,km1))/max(maxc,small)
            endif

! If the two layers differ significantly use a geometric averaging procedure
! But only do that for deep convection.  For shallow, use the simple
! averaging which is used in subr cmfmca
            if (iconvtype /= 1) then
               chat(m,k) = 0.5_r8* (const(m,k)+const(m,km1))
            else if (cdifr > 1.E-6_r8) then
               cabv = max(const(m,km1),maxc*1.e-12_r8)
               cbel = max(const(m,k),maxc*1.e-12_r8)
               chat(m,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel
            else             ! Small diff, so just arithmetic mean
               chat(m,k) = 0.5_r8* (const(m,k)+const(m,km1))
            end if

! Set provisional up and down draft values, and tendencies
            conu(m,k) = chat(m,k)
            cond(m,k) = chat(m,k)
         end do ! k

! Values at surface inferface == values in lowest layer
         chat(m,pver+1) = const(m,pver)
         conu(m,pver+1) = const(m,pver)
         cond(m,pver+1) = const(m,pver)
      end if
      end do ! m




! Compute updraft mixing ratios from cloudbase to cloudtop
! No special treatment is needed at k=pver because arrays
!    are dimensioned 1:pver+1
! A time-split approach is used.  First, entrainment is applied to produce
!    an initial conu(m,k) from conu(m,k+1).  Next, chemistry/physics are
!    applied to the initial conu(m,k) to produce a final conu(m,k).
!    Detrainment from the updraft uses this final conu(m,k).
! Note that different time-split approaches would give somewhat different
!    results
      kactcnt = 0 ; kactcntb = 0 ; kactfirst = 1
k_loop_main_bb: &
      do k = kbot, ktop, -1
         kp1 = k+1

! cldfrac = conv cloud fractional area.  This could represent anvil cirrus area,
!    and may not useful for aqueous chem and wet removal calculations
         cldfrac_i(k) = max( cldfrac_i(k), 0.005_r8 )
! mu_p_eudp(k) = updraft massflux at k, without detrainment between kp1,k
         mu_p_eudp(k) = mu_i(kp1) + eudp(k)

         fa_u(k) = 0.0_r8 !BSINGH(10/15/2014): Initialized so that it has a value if the following "if" check yeilds .false.
         if (mu_p_eudp(k) > mbsth) then
! if (mu_p_eudp(k) <= mbsth) the updraft mass flux is negligible at base and top
!    of current layer,
! so current layer is a "gap" between two unconnected updrafts,
! so essentially skip all the updraft calculations for this layer

! First apply changes from entrainment
            f_ent = eudp(k)/mu_p_eudp(k)
            f_ent = max( 0.0_r8, min( 1.0_r8, f_ent ) )
            tmpa = 1.0_r8 - f_ent
            do m = 1, pcnst_extd
              if (doconvproc_extd(m)) then
                  conu(m,k) = tmpa*conu(m,kp1) + f_ent*const(m,k)
              end if
            end do

! estimate updraft velocity (wup)
            if (iconvtype /= 1) then
! shallow - wup = (mup in kg/m2/s) / [rhoair * (updraft area)]
               wup(k) = (mu_i(kp1) + mu_i(k))*0.5_r8*hund_ovr_g &
                      / (rhoair_i(k) * (cldfrac_i(k)*0.5_r8))
            else
! deep - as in shallow, but assumed constant updraft_area with height zm_areafrac
               wup(k) = (mu_i(kp1) + mu_i(k))*0.5_r8*hund_ovr_g &
                      / (rhoair_i(k) * zm_areafrac)
            end if

! compute lagrangian transport time (dt_u) and updraft fractional area (fa_u)
! *** these must obey    dt_u(k)*mu_p_eudp(k) = dp_i(k)*fa_u(k)
            dt_u(k) = dz/wup(k)
            dt_u(k) = min( dt_u(k), dt )
            fa_u(k) = dt_u(k)*(mu_p_eudp(k)/dp_i(k))


! Now apply transformation and removal changes
!    Skip levels where icwmr(icol,k) <= clw_cut (= 1.0e-6) to eliminate
!    occasional very small icwmr values from the ZM module
            clw_cut = 1.0e-6_r8


            if (convproc_method_activate <= 1) then
! aerosol activation - method 1
!    skip levels that are completely glaciated (fracice(icol,k) == 1.0)
!    when kactcnt=1 (first/lowest layer with cloud water) apply
!       activatation to the entire updraft
!    when kactcnt>1 apply activatation to the amount entrained at this level
               if ((icwmr(icol,k) > clw_cut) .and. (fracice(icol,k) < 1.0_r8)) then
                  kactcnt = kactcnt + 1

                  if ((kactcnt == 1) .or. (f_ent > 0.0_r8)) then
                     kactcntb = kactcntb + 1
                  end if

                  if (kactcnt == 1) then
                     ! diagnostic fields
                     ! xx_wcldbase = w at first cloudy layer, estimated from mu and cldfrac
                     xx_wcldbase(icol) = (mu_i(kp1) + mu_i(k))*0.5_r8*hund_ovr_g &
                         / (rhoair_i(k) * (cldfrac_i(k)*0.5_r8))
                     xx_kcldbase(icol) = k

                     kactfirst = k
                     tmpa = 1.0_r8
                     call ma_activate_convproc(                            &
                        state, pbuf,                                       &
                        conu(:,k),  dconudt_activa(:,k), conu(:,k),        &
                        tmpa,       dt_u(k),             wup(k),           &
                        t(icol,k),  rhoair_i(k),         fracice(icol,k),  &
                        pcnst_extd, &
                        lchnk,      icol,                k,                &
                        ipass_calc_updraft                                 )
                  else if (f_ent > 0.0_r8) then
                     ! current layer is above cloud base (=first layer with activation)
                     !    only allow activation at k = kactfirst thru kactfirst-(method1_activate_nlayers-1)
                     if (k >= kactfirst-(method1_activate_nlayers-1)) then
                        call ma_activate_convproc(                            &
                           state, pbuf,                                       &
                           conu(:,k),  dconudt_activa(:,k), const(:,k),       &
                           f_ent,      dt_u(k),             wup(k),           &
                           t(icol,k),  rhoair_i(k),         fracice(icol,k),  &
                           pcnst_extd, &
                           lchnk,      icol,                k,                &
                           ipass_calc_updraft                                 )
                     end if
                  end if
! the following was for cam2 shallow convection (hack),
! but is not appropriate for cam5 (uwshcu)
!                 else if ((kactcnt > 0) .and. (iconvtype /= 1)) then
! !    for shallow conv, when you move from activation occuring to
! !       not occuring, reset kactcnt=0, because the hack scheme can
! !       produce multiple "1.5 layer clouds" separated by clear air
!                    kactcnt = 0
!                 end if
               end if ! ((icwmr(icol,k) > clw_cut) .and. (fracice(icol,k) < 1.0)) then

            else ! (convproc_method_activate >= 2)
! aerosol activation - method 2
!    skip levels that are completely glaciated (fracice(icol,k) == 1.0)
!    when kactcnt=1 (first/lowest layer with cloud water)
!       apply "primary" activatation to the entire updraft
!    when kactcnt>1
!       apply secondary activatation to the entire updraft
!       do this for all levels above cloud base (even if completely glaciated)
!          (this is something for sensitivity testing)
               do_act_this_lev = .false.
               if (kactcnt <= 0) then
                  if (icwmr(icol,k) > clw_cut) then
                     do_act_this_lev = .true.
                     kactcnt = 1
                     kactfirst = k
                     ! diagnostic fields
                     ! xx_wcldbase = w at first cloudy layer, estimated from mu and cldfrac
                     xx_wcldbase(icol) = (mu_i(kp1) + mu_i(k))*0.5_r8*hund_ovr_g &
                         / (rhoair_i(k) * (cldfrac_i(k)*0.5_r8))
                     xx_kcldbase(icol) = k
                  end if
               else
!                 if ((icwmr(icol,k) > clw_cut) .and. (fracice(icol,k) < 1.0)) then
                     do_act_this_lev = .true.
                     kactcnt = kactcnt + 1
!                 end if
               end if

               if ( do_act_this_lev ) then
                  kactcntb = kactcntb + 1

                  call ma_activate_convproc_method2(                    &
                     state, pbuf,                                       &
                     conu(:,k),  dconudt_activa(:,k),                   &
                     f_ent,      dt_u(k),             wup(k),           &
                     t(icol,k),  rhoair_i(k),         fracice(icol,k),  &
                     pcnst_extd, &
                     lchnk,      icol,                k,                &
                     kactfirst,  ipass_calc_updraft                     )
               end if
               conu2(icol,k,1:pver) = conu(1:pver,k)
            end if ! (convproc_method_activate <= 1)

! wet removal
!
! mirage2
!    rprd               = precip formation as a grid-cell average (kgW/kgA/s)
!    icwmr              = cloud water MR within updraft area (kgW/kgA)
!    fupdr              = updraft fractional area (--)
!    A = rprd/fupdr     = precip formation rate within updraft area (kgW/kgA/s)
!    B = A/icwmr = rprd/(icwmr*fupdr)
!                       = first-order removal rate (1/s)
!    C = dp/(mup/fupdr) = updraft air residence time in the layer (s)
!
!    fraction removed = (1.0 - exp(-cdt)) where
!                 cdt = B*C = (dp/mup)*rprd/icwmr
!
!    Note1:  fupdr cancels out in cdt, so need not be specified
!    Note2:  dp & mup units need only be consistent (e.g., mb & mb/s)
!    Note3:  for shallow conv, cdt = 1-beta (beta defined in Hack scheme)
!    Note4:  the "dp" in C above and code below should be the moist dp
!
! cam5
!    clw_preloss = cloud water MR before loss to precip
!                = icwmr + dt*(rprd/fupdr)
!    B = A/clw_preloss  = (rprd/fupdr)/(icwmr + dt*rprd/fupdr)
!                       = rprd/(fupdr*icwmr + dt*rprd)
!                       = first-order removal rate (1/s)
!
!    fraction removed = (1.0 - exp(-cdt)) where
!                 cdt = B*C = (fupdr*dp/mup)*[rprd/(fupdr*icwmr + dt*rprd)]
!
!    Note1:  *** cdt is now sensitive to fupdr, which we do not really know,
!                and is not the same as the convective cloud fraction
!    Note2:  dt is appropriate in the above cdt expression, not dtsub
!
!    Apply wet removal at levels where
!       icwmr(icol,k) > clw_cut  AND  rprd(icol,k) > 0.0
!    as wet removal occurs in both liquid and ice clouds


            cdt(k) = 0.0_r8
            if ((icwmr(icol,k) > clw_cut) .and. (rprd(icol,k) > 0.0_r8)) then
                  tmpf = 0.5_r8*cldfrac_i(k)
                  cdt(k) = (tmpf*dp(i,k)/mu_p_eudp(k)) * rprd(icol,k) / &
                        (tmpf*icwmr(icol,k) + dt*rprd(icol,k))
            end if

            if (cdt(k) > 0.0_r8) then
               expcdtm1 = exp(-cdt(k)) - 1.0_r8
               do m = 1, pcnst_extd
                  if (doconvproc_extd(m)) then
                     dconudt_wetdep(m,k) = conu(m,k)*aqfrac(m)*expcdtm1
                     conu(m,k) = conu(m,k) + dconudt_wetdep(m,k)
                     dconudt_wetdep(m,k) = dconudt_wetdep(m,k) / dt_u(k)
                     conu2(icol,k,m) = conu(m,k)
                  end if
               enddo
            end if

         end if    ! "(mu_p_eudp(k) > mbsth)"
      end do k_loop_main_bb ! "k = kbot, ktop, -1"

! when doing updraft calcs twice, only need to go this far on the first pass
      if ( (ipass_calc_updraft == 1) .and. &
           (npass_calc_updraft == 2) ) cycle ipass_calc_updraft_loop


! Compute downdraft mixing ratios from cloudtop to cloudbase
! No special treatment is needed at k=2
! No transformation or removal is applied in the downdraft
      do k = ktop, kbot
         kp1 = k + 1
! md_m_eddp = downdraft massflux at kp1, without detrainment between k,kp1
         md_m_eddp = md_i(k) - eddp(k)
         if (md_m_eddp < -mbsth) then
            do m = 1, pcnst_extd
               if (doconvproc_extd(m)) then
                  cond(m,kp1) = ( md_i(k)*cond(m,k)		&
                                - eddp(k)*const(m,k) ) / md_m_eddp
               endif
            end do
         end if
      end do ! k


! Now computes fluxes and tendencies
! NOTE:  The approach used in convtran applies to inert tracers and
!        must be modified to include source and sink terms
      sumflux(:) = 0.0_r8
      sumflux2(:) = 0.0_r8
      sumsrce(:) = 0.0_r8
      sumchng(:) = 0.0_r8
      sumchng3(:) = 0.0_r8
      sumwetdep(:) = 0.0_r8
      sumprevap(:) = 0.0_r8

      maxflux(:) = 0.0_r8
      maxflux2(:) = 0.0_r8
      maxresusp(:) = 0.0_r8
      maxsrce(:) = 0.0_r8
      maxprevap(:) = 0.0_r8

k_loop_main_cc: &
      do k = ktop, kbot
         kp1 = k+1
         km1 = k-1
         kp1x = min( kp1, pver )
         km1x = max( km1, 1 )
         fa_u_dp = fa_u(k)*dp_i(k)

         do m = 1, pcnst_extd
            if (doconvproc_extd(m)) then

! First compute fluxes using environment subsidence/lifting and
! entrainment/detrainment into up/downdrafts,
! to provide an additional mass balance check
! (this could be deleted after the code is well tested)
               fluxin  = mu_i(k)*min(chat(m,k),const(m,km1x))       &
                       - md_i(kp1)*min(chat(m,kp1),const(m,kp1x))   &
                       + dudp(k)*conu(m,k) + dddp(k)*cond(m,kp1)
               fluxout = mu_i(kp1)*min(chat(m,kp1),const(m,k))      &
                       - md_i(k)*min(chat(m,k),const(m,k))          &
                       + (eudp(k) + eddp(k))*const(m,k)

               netflux = fluxin - fluxout

               sumflux2(m) = sumflux2(m) + netflux
               maxflux2(m) = max( maxflux2(m), abs(fluxin), abs(fluxout) )

! Now compute fluxes as in convtran, and also source/sink terms
! (version 3 limit fluxes outside convection to mass in appropriate layer
! (these limiters are probably only safe for positive definite quantitities
! (it assumes that mu and md already satify a courant number limit of 1)
            if (iflux_method /= 2) then
               fluxin  =     mu_i(kp1)*conu(m,kp1)                     &
                           + mu_i(k  )*min(chat(m,k  ),const(m,km1x))  &
                         - ( md_i(k  )*cond(m,k)                       &
                           + md_i(kp1)*min(chat(m,kp1),const(m,kp1x)) )
               fluxout =     mu_i(k  )*conu(m,k)                       &
                           + mu_i(kp1)*min(chat(m,kp1),const(m,k   ))  &
                         - ( md_i(kp1)*cond(m,kp1)                     &
                           + md_i(k  )*min(chat(m,k  ),const(m,k   )) )
            else
               fluxin  =     mu_i(kp1)*conu(m,kp1)                     &
                         - ( md_i(k  )*cond(m,k) )
               fluxout =     mu_i(k  )*conu(m,k)                       &
                         - ( md_i(kp1)*cond(m,kp1) )
               tmpveca(1) = fluxin ; tmpveca(4) = -fluxout

               ! new method -- simple upstream method for the env subsidence
               ! tmpa = net env mass flux (positive up) at top of layer k
               tmpa = -( mu_i(k  ) + md_i(k  ) )
               if (tmpa <= 0.0_r8) then
                  fluxin  = fluxin  - tmpa*const(m,km1x)
               else
                  fluxout = fluxout + tmpa*const(m,k   )
               end if
               tmpveca(2) = fluxin ; tmpveca(5) = -fluxout
               ! tmpa = net env mass flux (positive up) at base of layer k
               tmpa = -( mu_i(kp1) + md_i(kp1) )
               if (tmpa >= 0.0_r8) then
                  fluxin  = fluxin  + tmpa*const(m,kp1x)
               else
                  fluxout = fluxout - tmpa*const(m,k   )
               end if
               tmpveca(3) = fluxin ; tmpveca(6) = -fluxout
            end if

            netflux = fluxin - fluxout
            netsrce = fa_u_dp*(dconudt_aqchem(m,k) + &
                        dconudt_activa(m,k) + dconudt_wetdep(m,k))

            dcondt(m,k) = (netflux+netsrce)/dp_i(k)

            dcondt_wetdep(m,k) = fa_u_dp*dconudt_wetdep(m,k)/dp_i(k)
            sumwetdep(m) = sumwetdep(m) + fa_u_dp*dconudt_wetdep(m,k)

            end if   ! "(doconvproc_extd(m))"

         end do      ! "m = 2,ncnst_extd"
      end do k_loop_main_cc ! "k = ktop, kbot"


! calculate effects of precipitation evaporation
      call ma_precpevap_convproc( dcondt, dcondt_wetdep,  dcondt_prevap,   &
                                  rprd,   evapc,          dp_i,            &
                                  icol,   ktop,           pcnst_extd,      &
                                  lchnk,           &
                                  doconvproc_extd                          )

! make adjustments to dcondt for activated & unactivated aerosol species
!    pairs to account any (or total) resuspension of convective-cloudborne aerosol
      call ma_resuspend_convproc( dcondt, dcondt_resusp,   &
                                  const, dp_i, ktop, kbot_prevap, pcnst_extd )

      ! Do resuspension of aerosols from rain only when the rain has
      ! totally evaporated.
      if (convproc_do_evaprain_atonce) then
         dcondt_resusp3d(ncnst_tot+1:pcnst_extd,icol,:) = dcondt_resusp(ncnst_tot+1:pcnst_extd,:)
         dcondt_resusp(ncnst_tot+1:pcnst_extd,:) = 0._r8
      end if

! calculate new column-tendency variables
      do m = 1, pcnst_extd
         if (doconvproc_extd(m)) then
            do k = ktop, kbot_prevap
               sumprevap(m) = sumprevap(m) + dcondt_prevap(m,k)*dp_i(k)
            end do
         end if
      end do ! m

!
! note again the ma_convproc_tend does not apply convective cloud processing
!    to the stratiform-cloudborne aerosol
! within this routine, cloudborne aerosols are convective-cloudborne
!
! before tendencies (dcondt, which is loaded into dqdt) are returned,
!    the convective-cloudborne aerosol tendencies must be combined
!    with the interstitial tendencies
! ma_resuspend_convproc has already done this for the dcondt
!
! the individual process column tendencies (sumwetdep, sumprevap, ...)
!    are just diagnostic fields that can be written to history
! tendencies for interstitial and convective-cloudborne aerosol could
!    both be passed back and output, if desired
! currently, however, the interstitial and convective-cloudborne tendencies
!    are combined (in the next code block) before being passed back (in qsrflx)
!
      do n = 1, nbins
        do ll = 1, nspec(n) + 2
           l = bin_idx(n,ll)
           la = l
           lc = l + ncnst_tot

           if (doconvproc(la)) then
               sumwetdep(la) = sumwetdep(la) + sumwetdep(lc)
               sumprevap(la) = sumprevap(la) + sumprevap(lc)
            end if
         enddo ! ll
      enddo ! n

!
! scatter overall tendency back to full array
!
      do m = 1, ncnst_tot
         if (doconvproc(m)) then
            do k = ktop, kbot_prevap
               dqdt_i(k,m) = dcondt(m,k)
               dqdt(icol,k,m) = dqdt(icol,k,m) + dqdt_i(k,m)*xinv_ntsub
            end do
!           dqdt_i(:,m) = 0.
         end if
      end do ! m

! scatter column burden tendencies for various processes to qsrflx
      do m = 1, ncnst_tot
         if (doconvproc(m)) then
            qsrflx_i(m,4) = sumwetdep(m)*hund_ovr_g
            qsrflx_i(m,5) = sumprevap(m)*hund_ovr_g
            qsrflx(icol,m,1:5) = qsrflx(icol,m,1:5) + qsrflx_i(m,1:5)*xinv_ntsub
         end if
      end do ! m

      if (jtsub < ntsub) then
         ! update the q_i for the next interation of the jtsub loop
         do m = 1, ncnst_tot
            if (doconvproc(m)) then
               do k = ktop, kbot_prevap
                  q_i(k,m) = max( (q_i(k,m) + dqdt_i(k,m)*dtsub), 0.0_r8 )
               end do
            end if
         end do ! m
      end if

      end do ipass_calc_updraft_loop

      end do jtsub_loop_main_aa  ! of the main "do jtsub = 1, ntsub" loop


   end do i_loop_main_aa  ! of the main "do i = il1g, il2g" loop

   do n = 1, nbins
      do ll = 1, nspec(n) + 2
         l = bin_idx(n,ll)
         la = l
         lc = l + ncnst_tot
         ! q(icol,1:pver,1:ncnst_tot)
         call outfld( trim(fieldname(l))//'QCONST', q(:,:,la), pcols, lchnk )
         call outfld( trim(fieldname(l))//'CONU', conu2(:,:,la), pcols, lchnk )
         call outfld( trim(fieldname(l))//'WETC', dcondt2(:,:,la), pcols, lchnk )
         call outfld( trim(fieldname_cw(l))//'CONU', conu2(:,:,lc), pcols, lchnk )
         call outfld( trim(fieldname_cw(l))//'WETC', dcondt2(:,:,lc), pcols, lchnk )
      end do
   end do

   return
end subroutine ma_convproc_tend



!=========================================================================================
   subroutine ma_precpevap_convproc(                           &
              dcondt,  dcondt_wetdep, dcondt_prevap,           &
              rprd,    evapc,         dp_i,                    &
              icol,    ktop,          pcnst_extd,              &
              lchnk,                   &
              doconvproc_extd                                  )
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate resuspension of wet-removed aerosol species resulting
!    precip evaporation
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

   implicit none

!-----------------------------------------------------------------------
! arguments
! (note:  TMR = tracer mixing ratio)
   integer,  intent(in)    :: pcnst_extd

   real(r8), intent(inout) :: dcondt(pcnst_extd,pver)
                              ! overall TMR tendency from convection
   real(r8), intent(in)    :: dcondt_wetdep(pcnst_extd,pver)
                              ! portion of TMR tendency due to wet removal
   real(r8), intent(inout) :: dcondt_prevap(pcnst_extd,pver)
                              ! portion of TMR tendency due to precip evaporation
                              ! (actually, due to the adjustments made here)
                              ! (on entry, this is 0.0)

   real(r8), intent(in)    :: rprd(pcols,pver)  ! conv precip production  rate (gathered)
   real(r8), intent(in)    :: evapc(pcols,pver)  ! conv precip evaporation rate (gathered)
   real(r8), intent(in)    :: dp_i(pver) ! pressure thickness of level (in mb)

   integer,  intent(in)    :: icol  ! normal (ungathered) i index for current column
   integer,  intent(in)    :: ktop  ! index of top cloud level for current column
   integer,  intent(in)    :: lchnk  ! chunk index

   logical,  intent(in)    :: doconvproc_extd(pcnst_extd)  ! indicates which species to process

!-----------------------------------------------------------------------
! local variables
   integer  :: k, l, ll, m, n
   real(r8) :: del_pr_flux_prod      ! change to precip flux from production  [(kg/kg/s)*mb]
   real(r8) :: del_pr_flux_evap      ! change to precip flux from evaporation [(kg/kg/s)*mb]
   real(r8) :: del_wd_flux_evap      ! change to wet deposition flux from evaporation [(kg/kg/s)*mb]
   real(r8) :: fdel_pr_flux_evap     ! fractional change to precip flux from evaporation
   real(r8) :: pr_flux               ! precip flux at base of current layer [(kg/kg/s)*mb]
   real(r8) :: pr_flux_old
   real(r8) :: tmpa, tmpb, tmpc, tmpd
   real(r8) :: tmpdp                 ! delta-pressure (mb)
   real(r8) :: wd_flux(pcnst_extd)   ! tracer wet deposition flux at base of current layer [(kg/kg/s)*mb]
   integer :: i
   character(len=4) :: spcstr
!-----------------------------------------------------------------------


   pr_flux = 0.0_r8
   wd_flux(:) = 0.0_r8

   do k = ktop, pver
      tmpdp = dp_i(k)

      pr_flux_old = pr_flux
      del_pr_flux_prod = tmpdp*max(0.0_r8, rprd(icol,k))
      pr_flux = pr_flux_old + del_pr_flux_prod

      del_pr_flux_evap = min( pr_flux, tmpdp*max(0.0_r8, evapc(icol,k)) )

      ! Do resuspension of aerosols from rain only when the rain has
      ! totally evaporated in one layer.
      if (convproc_do_evaprain_atonce .and. &
          (del_pr_flux_evap.ne.pr_flux)) del_pr_flux_evap = 0._r8

      fdel_pr_flux_evap = del_pr_flux_evap / max(pr_flux, 1.0e-35_r8)

      do m = 1, pcnst_extd
         if ( doconvproc_extd(m) ) then
            ! use -dcondt_wetdep(m,k) as it is negative (or zero)
            wd_flux(m) = wd_flux(m) + tmpdp*max(0.0_r8, -dcondt_wetdep(m,k))
            del_wd_flux_evap = wd_flux(m)*fdel_pr_flux_evap
            wd_flux(m) = max( 0.0_r8, wd_flux(m)-del_wd_flux_evap )

            dcondt_prevap(m,k) = del_wd_flux_evap/tmpdp
            dcondt(m,k) = dcondt(m,k) + dcondt_prevap(m,k)
         end if
      end do

     ! move dcondt_prevap to larger bin
      do m = 1, ncnst_tot
        if (ibl(m).ne.m) then
          dcondt_prevap(ibl(m),k) = dcondt_prevap(ibl(m),k) + dcondt_prevap(m,k)
          dcondt_prevap(m,k) = 0._r8
        end if
     end do
   end do ! k

     !   return
   end subroutine ma_precpevap_convproc

!=========================================================================================
   subroutine accumulate_to_larger_mode( spc_name, lptr, prevap )

     character(len=*), intent(in) :: spc_name
     integer,  intent(in) :: lptr(:)
     real(r8), intent(inout) :: prevap(:)

     integer :: m,n, nl,ns

     !     ! find constituent index of the largest mode for the species
     !     loop1: do m = 1,ntot_amode-1
     !        nl = lptr(mode_size_order(m))
     !        if (nl>0) exit loop1
     !     end do loop1
     !
     !     if (.not. nl>0) return
     !
     !     ! accumulate the smaller modes into the largest mode
     !     do n = m+1,ntot_amode
     !        ns = lptr(mode_size_order(n))
     !        if (ns>0) then
     !           prevap(nl) = prevap(nl) + prevap(ns)
     !           prevap(ns) = 0._r8
     !           if (masterproc .and. debug) then
     !              write(iulog,'(a,i3,a,i3)') trim(spc_name)//' mode number accumulate ',ns,'->',nl
     !           endif
     !        endif
     !     end do

   end subroutine accumulate_to_larger_mode

!=========================================================================================
   subroutine ma_activate_convproc(             &
              state, pbuf,                      &
              conu,       dconudt,   conent,    &
              f_ent,      dt_u,      wup,       &
              tair,       rhoair,    fracice,   &
              pcnst_extd, &
              lchnk,      i,         k,         &
              ipass_calc_updraft                )
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate activation of aerosol species in convective updraft
! for a single column and level
!
! Method:
! conu(l)    = Updraft TMR (tracer mixing ratio) at k/k-1 interface
! conent(l)  = TMR of air that is entrained into the updraft from level k
! f_ent      = Fraction of the "before-detrainment" updraft massflux at
!              k/k-1 interface" resulting from entrainment of level k air
!              (where k is the current level in subr ma_convproc_tend)
!
! On entry to this routine, the conu(l) represents the updraft TMR
! after entrainment, but before chemistry/physics and detrainment,
! and is equal to
!    conu(l) = f_ent*conent(l) + (1.0-f_ent)*conu_below(l)
! where
!    conu_below(l) = updraft TMR at the k+1/k interface, and
!    f_ent   = (eudp/mu_p_eudp) is the fraction of the updraft massflux
!              from level k entrainment
!
! This routine applies aerosol activation to the entrained tracer,
! then adjusts the conu so that on exit,
!   conu(la) = conu_incoming(la) - f_ent*conent(la)*f_act(la)
!   conu(lc) = conu_incoming(lc) + f_ent*conent(la)*f_act(la)
! where
!   la, lc   = indices for an unactivated/activated aerosol component pair
!   f_act    = fraction of conent(la) that is activated.  The f_act are
!              calculated with the Razzak-Ghan activation parameterization.
!              The f_act differ for each mode, and for number/surface/mass.
!
! Note:  At the lowest layer with cloud water, subr convproc calls this
! routine with conent==conu and f_ent==1.0, with the result that
! activation is applied to the entire updraft tracer flux
!
! *** The updraft velocity used for activation calculations is rather
!     uncertain and needs more work.  However, an updraft of 1-3 m/s
!     will activate essentially all of accumulation and coarse mode particles.
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

   use ndrop, only: activate_aerosol
   use ppgrid,          only: pcols

   implicit none

!-----------------------------------------------------------------------
! arguments  (note:  TMR = tracer mixing ratio)
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

   integer, intent(in)     :: pcnst_extd
   ! conu = tracer mixing ratios in updraft at top of this (current) level
   !        The conu are changed by activation
   real(r8), intent(inout) :: conu(pcnst_extd)
   ! conent = TMRs in the entrained air at this level
   real(r8), intent(in)    :: conent(pcnst_extd)
   real(r8), intent(inout) :: dconudt(pcnst_extd) ! TMR tendencies due to activation

   real(r8), intent(in)    :: f_ent  ! fraction of updraft massflux that was
                                     ! entrained across this layer == eudp/mu_p_eudp
   real(r8), intent(in)    :: dt_u   ! lagrangian transport time (s) in the
                                     ! updraft at current level
   real(r8), intent(in)    :: wup    ! mean updraft vertical velocity (m/s)
                                     ! at current level updraft

   real(r8), intent(in)    :: tair   ! Temperature in Kelvin
   real(r8), intent(in)    :: rhoair ! air density (kg/m3)

   real(r8), intent(in)    :: fracice ! Fraction of ice within the cloud
                                     ! used as in-cloud wet removal rate
   integer,  intent(in)    :: lchnk  ! chunk index
   integer,  intent(in)    :: i      ! column index
   integer,  intent(in)    :: k      ! level index
   integer,  intent(in)    :: ipass_calc_updraft

!-----------------------------------------------------------------------
! local variables
   integer  :: ll, la, lc, n, l
   integer  :: phase ! phase of aerosol

   real(r8) :: delact                ! working variable
   real(r8) :: dt_u_inv              ! 1.0/dt_u
   real(r8) :: fluxm(nbins)      ! to understand this, see subr activate_modal
   real(r8) :: fluxn(nbins)      ! to understand this, see subr activate_modal
   real(r8) :: flux_fullact           ! to understand this, see subr activate_modal
   real(r8) :: fm(nbins)         ! mass fraction of aerosols activated
   real(r8) :: fn(nbins)         ! number fraction of aerosols activated
   real(r8) :: hygro(nbins)      ! current hygroscopicity for int+act
   real(r8) :: maerosol(nbins)   ! interstitial+activated total mixing ratios conc (m3/m3)
   real(r8) :: naerosol(nbins)   ! interstitial+activated number conc (#/m3)
   real(r8) :: sigw                  ! standard deviation of updraft velocity (cm/s)
   real(r8) :: tmpa, tmpb, tmpc      ! working variable
   real(r8) :: tmp_fact              ! working variable
   real(r8) :: vaerosol(nbins)   ! int+act volume (m3/m3)
   real(r8) :: wbar                  ! mean updraft velocity (cm/s)
   real(r8) :: wdiab                 ! diabatic vertical velocity (cm/s)
   real(r8) :: wminf, wmaxf          ! limits for integration over updraft spectrum (cm/s)
   real(r8) :: specdens              ! specie density from physprops files
   real(r8) :: spechygro             ! hygroscopicity from physprops files

   type(carma_aerosol_state), pointer :: aero_state_obj

   real(r8) :: cs_a(pcols,pver)    ! air density (kg/m3)
   real(r8) :: naerosol_a(pcols)  ! number conc (1/m3)
   real(r8) :: vaerosol_a(pcols)  ! volume conc (m3/m3)
   real(r8) :: hygro_a(pcols)     ! bulk hygroscopicity of mode

   integer :: errnum
   character(len=shr_kind_cs) :: errstr
!-----------------------------------------------------------------------

   nullify(aero_state_obj)

! when ipass_calc_updraft == 2, apply the activation tendencies
!    from pass 1, but multiplied by factor_reduce_actfrac
! (can only have ipass_calc_updraft == 2 when method_reduce_actfrac = 2)
   if (ipass_calc_updraft == 2) then

   dt_u_inv = 1.0_r8/dt_u
   do n = 1, nbins
      do ll = 1, nspec(n) + 2
         l = bin_idx(n,ll)

         lc = l + ncnst_tot
         la = l

         delact = dconudt(lc)*dt_u * factor_reduce_actfrac
         delact = min( delact, conu(la) )
         delact = max( delact, 0.0_r8 )
         conu(la) = conu(la) - delact
         conu(lc) = conu(lc) + delact
         dconudt(la) = -delact*dt_u_inv
         dconudt(lc) =  delact*dt_u_inv
      end do
   end do   ! "n = 1, ntot_amode"
   return

   end if ! (ipass_calc_updraft == 2)

! check f_ent > 0
   if (f_ent <= 0.0_r8) return

! call Razzak-Ghan activation routine with single updraft
   wbar = max( wup, 0.5_r8 )  ! force wbar >= 0.5 m/s for now
   sigw = 0.0_r8
   wdiab = 0.0_r8
   wminf = wbar
   wmaxf = wbar

   phase = 1 ! interstitial
   aero_state_obj => carma_aerosol_state(state, pbuf)
   cs_a(i,k) = rhoair
   do n = 1, nbins
      call aero_state_obj%loadaer( aero_props_obj, i, i, k, n, cs_a, phase, naerosol_a, vaerosol_a, hygro_a, errnum, errstr)
      if (errnum/=0) then
         call endrun('ma_activate_convproc : '//trim(errstr))
      end if
      naerosol(n) = naerosol_a(i)
      vaerosol(n) = vaerosol_a(i)
      hygro(n) = hygro_a(i)
   end do
   deallocate(aero_state_obj)
   nullify(aero_state_obj)

   call activate_aerosol(                                                    &
         wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,                    &
         naerosol, nbins, vaerosol, hygro, aero_props_obj,                 &
         fn, fm, fluxn, fluxm, flux_fullact                                )

! apply the activation fractions to the updraft aerosol mixing ratios
   dt_u_inv = 1.0_r8/dt_u

   do n = 1, nbins
      do ll = 1, nspec(n) + 2
         l = bin_idx(n,ll)
         la = l
         lc = l + ncnst_tot

         if (ll == nspec(n) + 2) then
            tmp_fact = fn(n)
         else
            tmp_fact = fm(n)
         end if

         if ( (method_reduce_actfrac == 1)      .and. &
              (factor_reduce_actfrac >= 0.0_r8) .and. &
              (factor_reduce_actfrac <  1.0_r8) )     &
              tmp_fact = tmp_fact * factor_reduce_actfrac

         delact = min( conent(la)*tmp_fact*f_ent, conu(la) )
         delact = max( delact, 0.0_r8 )
         conu(la) = conu(la) - delact
         conu(lc) = conu(lc) + delact
         dconudt(la) = -delact*dt_u_inv
         dconudt(lc) =  delact*dt_u_inv
      end do
   end do   ! "n = 1, ntot_amode"

   return
   end subroutine ma_activate_convproc



!=========================================================================================
   subroutine ma_activate_convproc_method2(     &
              state, pbuf,                      &
              conu,       dconudt,              &
              f_ent,      dt_u,      wup,       &
              tair,       rhoair,    fracice,   &
              pcnst_extd, &
              lchnk,      i,         k,         &
              kactfirst,  ipass_calc_updraft    )
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate activation of aerosol species in convective updraft
! for a single column and level
!
! Method:
! conu(l)    = Updraft TMR (tracer mixing ratio) at k/k-1 interface
! f_ent      = Fraction of the "before-detrainment" updraft massflux at
!              k/k-1 interface" resulting from entrainment of level k air
!              (where k is the current level in subr ma_convproc_tend)
!
! On entry to this routine, the conu(l) represents the updraft TMR
! after entrainment, but before chemistry/physics and detrainment.
!
! This routine applies aerosol activation to the conu tracer mixing ratios,
! then adjusts the conu so that on exit,
!   conu(la) = conu_incoming(la) - conu(la)*f_act(la)
!   conu(lc) = conu_incoming(lc) + conu(la)*f_act(la)
! where
!   la, lc   = indices for an unactivated/activated aerosol component pair
!   f_act    = fraction of conu(la) that is activated.  The f_act are
!              calculated with the Razzak-Ghan activation parameterization.
!              The f_act differ for each mode, and for number/surface/mass.
!
! At cloud base (k==kactfirst), primary activation is done using the
! "standard" code in subr activate do diagnose maximum supersaturation.
! Above cloud base, secondary activation is done using a
! prescribed supersaturation.
!
! *** The updraft velocity used for activation calculations is rather
!     uncertain and needs more work.  However, an updraft of 1-3 m/s
!     will activate essentially all of accumulation and coarse mode particles.
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

   use ndrop, only: activate_aerosol
   use ppgrid,          only: pcols

   use rad_constituents,only: rad_cnst_get_info

   implicit none

!-----------------------------------------------------------------------
! arguments  (note:  TMR = tracer mixing ratio)
   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

   integer, intent(in)     :: pcnst_extd
   ! conu = tracer mixing ratios in updraft at top of this (current) level
   !        The conu are changed by activation
   real(r8), intent(inout) :: conu(pcnst_extd)
   real(r8), intent(inout) :: dconudt(pcnst_extd) ! TMR tendencies due to activation

   real(r8), intent(in)    :: f_ent  ! fraction of updraft massflux that was
                                     ! entrained across this layer == eudp/mu_p_eudp
   real(r8), intent(in)    :: dt_u   ! lagrangian transport time (s) in the
                                     ! updraft at current level
   real(r8), intent(in)    :: wup    ! mean updraft vertical velocity (m/s)
                                     ! at current level updraft

   real(r8), intent(in)    :: tair   ! Temperature in Kelvin
   real(r8), intent(in)    :: rhoair ! air density (kg/m3)

   real(r8), intent(in)    :: fracice ! Fraction of ice within the cloud
                                     ! used as in-cloud wet removal rate
   integer,  intent(in)    :: lchnk  ! chunk index
   integer,  intent(in)    :: i      ! column index
   integer,  intent(in)    :: k      ! level index
   integer,  intent(in)    :: kactfirst ! k at cloud base
   integer,  intent(in)    :: ipass_calc_updraft

!-----------------------------------------------------------------------
! local variables
   integer  :: ll, la, lc, n, l
   integer  :: phase ! phase of aerosol

   real(r8) :: delact                ! working variable
   real(r8) :: dt_u_inv              ! 1.0/dt_u
   real(r8) :: fluxm(nbins)      ! to understand this, see subr activate_modal
   real(r8) :: fluxn(nbins)      ! to understand this, see subr activate_modal
   real(r8) :: flux_fullact           ! to understand this, see subr activate_modal
   real(r8) :: fm(nbins)         ! mass fraction of aerosols activated
   real(r8) :: fn(nbins)         ! number fraction of aerosols activated
   real(r8) :: hygro(nbins)      ! current hygroscopicity for int+act
   real(r8) :: naerosol(nbins)   ! interstitial+activated number conc (#/m3)
   real(r8) :: maerosol(nbins)   ! interstitial+activated number conc (#/m3)
   real(r8) :: sigw                  ! standard deviation of updraft velocity (cm/s)
   real(r8) :: smax_prescribed       ! prescribed supersaturation for secondary activation (0-1 fraction)
   real(r8) :: tmpa, tmpb, tmpc      ! working variable
   real(r8) :: tmp_fact              ! working variable
   real(r8) :: vaerosol(nbins)   ! int+act volume (m3/m3)
   real(r8) :: wbar                  ! mean updraft velocity (cm/s)
   real(r8) :: wdiab                 ! diabatic vertical velocity (cm/s)
   real(r8) :: wminf, wmaxf          ! limits for integration over updraft spectrum (cm/s)
   real(r8) :: specdens              ! specie density from physprops files
   real(r8) :: spechygro             ! hygroscopicity from physprops files

   character(len=32) :: spectype

   type(carma_aerosol_state), pointer :: aero_state_obj

   real(r8) :: cs_a(pcols,pver)    ! air density (kg/m3)
   real(r8) :: naerosol_a(pcols)  ! number conc (1/m3)
   real(r8) :: vaerosol_a(pcols)  ! volume conc (m3/m3)
   real(r8) :: hygro_a(pcols)     ! bulk hygroscopicity of mode

   integer :: errnum
   character(len=shr_kind_cs) :: errstr
!-----------------------------------------------------------------------

   nullify(aero_state_obj)

! when ipass_calc_updraft == 2, apply the activation tendencies
!    from pass 1, but multiplied by factor_reduce_actfrac
! (can only have ipass_calc_updraft == 2 when method_reduce_actfrac = 2)
   if (ipass_calc_updraft == 2) then

   dt_u_inv = 1.0_r8/dt_u
   do n = 1, nbins
      do ll = 1, nspec(n) + 2
         l = bin_idx(n,ll)
         la = l
         lc = l + ncnst_tot

         delact = dconudt(lc)*dt_u * factor_reduce_actfrac
         delact = min( delact, conu(la) )
         delact = max( delact, 0.0_r8 )
         conu(la) = conu(la) - delact
         conu(lc) = conu(lc) + delact
         dconudt(la) = -delact*dt_u_inv
         dconudt(lc) =  delact*dt_u_inv

      end do
   end do   ! "n = 1, ntot_amode"
   return

   end if ! (ipass_calc_updraft == 2)


! check f_ent > 0
   if (f_ent <= 0.0_r8) return

   do n = 1, nbins
! compute a (or a+cw) volume and hygroscopicity
      tmpa = 0.0_r8
      tmpb = 0.0_r8
      do ll = 1, nspec(n)
         l = bin_idx(n,ll)

         call rad_cnst_get_bin_props_by_idx(0, n, ll, density_aer=specdens, hygro_aer=spechygro)

         tmpc = max( conu(l), 0.0_r8 )
         if ( use_cwaer_for_activate_maxsat ) &
          tmpc = tmpc + max( conu(l+ncnst_tot), 0.0_r8 )
          tmpc = tmpc / specdens
         tmpa = tmpa + tmpc

         ! Change the hygroscopicity of POM based on the discussion with Prof.
         ! Xiaohong Liu. Some observational studies found that the primary organic
         ! material from biomass burning emission shows very high hygroscopicity.
         ! Also, found that BC mass will be overestimated if all the aerosols in
         ! the primary mode are free to be removed. Therefore, set the hygroscopicity
         ! of POM here as 0.2 to enhance the wet scavenge of primary BC and POM.

         call rad_cnst_get_bin_props_by_idx(0, n, ll,spectype=spectype)
         if (spectype=='p-organic' .and. convproc_pom_spechygro>0._r8) then
             tmpb = tmpb + tmpc * convproc_pom_spechygro
         else
            tmpb = tmpb + tmpc * spechygro
         end if
      end do
      vaerosol(n) = tmpa * rhoair
      if (tmpa < 1.0e-35_r8) then
         hygro(n) = 0.2_r8
      else
         hygro(n) = tmpb/tmpa
      end if

! total mass not needed
! load a (or a+cw) total mmr  and bound it
!      l = bin_idx(n,nspec(n)+1)
!      tmpa = max( conu(l), 0.0_r8 )
!      if ( use_cwaer_for_activate_maxsat ) &
!      tmpa = tmpa + max( conu(l+ncnst_tot), 0.0_r8 )
!      maerosol(n) = tmpa * rhoair
!      maerosol(n) = max( vaerosol(n),maerosol(n))

! load a (or a+cw) number and bound it
      l = bin_idx(n,nspec(n)+2)
      tmpa = max( conu(l), 0.0_r8 )
      if ( use_cwaer_for_activate_maxsat ) &
      tmpa = tmpa + max( conu(l+ncnst_tot), 0.0_r8 )
!      naerosol(n) = tmpa * 1000.0_r8 * rhoair
      naerosol(n) = tmpa * rhoair
   end do

! call Razzak-Ghan activation routine with single updraft
   wbar = max( wup, 0.5_r8 )  ! force wbar >= 0.5 m/s for now
   sigw = 0.0_r8
   wdiab = 0.0_r8
   wminf = wbar
   wmaxf = wbar

   phase = 1 ! interstitial
   aero_state_obj => carma_aerosol_state(state, pbuf)
   cs_a(i,k) = rhoair
   do n = 1, nbins
      call aero_state_obj%loadaer( aero_props_obj, i, i, k, n, cs_a, phase, naerosol_a, vaerosol_a, hygro_a, errnum, errstr)
      if (errnum/=0) then
         call endrun('ma_activate_convproc_method2 : '//trim(errstr))
      end if
      naerosol(n) = naerosol_a(i)
      vaerosol(n) = vaerosol_a(i)
      hygro(n) = hygro_a(i)
   end do
   deallocate(aero_state_obj)
   nullify(aero_state_obj)

   if (k == kactfirst) then

      call activate_aerosol(                                                 &
         wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,                    &
         naerosol, nbins, vaerosol, hygro, aero_props_obj,                 &
         fn, fm, fluxn, fluxm, flux_fullact                                )


   else
! above cloud base - do secondary activation with prescribed supersat
! that is constant with height
      smax_prescribed = method2_activate_smaxmax
      call activate_aerosol(                                                 &
         wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,                    &
         naerosol, nbins, vaerosol, hygro, aero_props_obj,                 &
         fn, fm, fluxn, fluxm, flux_fullact, smax_prescribed               )
   end if

! apply the activation fractions to the updraft aerosol mixing ratios
   dt_u_inv = 1.0_r8/dt_u

   do n = 1, nbins
      do ll = 1, nspec(n) + 2
         l = bin_idx(n,ll)
         la = l
         lc = l + ncnst_tot

         if (ll == nspec(n) + 2) then
            tmp_fact = fn(n)
         else
            tmp_fact = fm(n)
         end if

         if ( (method_reduce_actfrac == 1)      .and. &
              (factor_reduce_actfrac >= 0.0_r8) .and. &
              (factor_reduce_actfrac <  1.0_r8) )     &
              tmp_fact = tmp_fact * factor_reduce_actfrac

         delact = min( conu(la)*tmp_fact, conu(la) )
         delact = max( delact, 0.0_r8 )
         conu(la) = conu(la) - delact
         conu(lc) = conu(lc) + delact
         dconudt(la) = -delact*dt_u_inv
         dconudt(lc) =  delact*dt_u_inv

      end do
   end do   ! "n = 1, nbin"

   return
   end subroutine ma_activate_convproc_method2



!=========================================================================================
   subroutine ma_resuspend_convproc(                           &
              dcondt,  dcondt_resusp,                          &
              const,   dp_i,          ktop,  kbot_prevap,  pcnst_extd )
!-----------------------------------------------------------------------
!
! Purpose:
! Calculate resuspension of activated aerosol species resulting from both
!    detrainment from updraft and downdraft into environment
!    subsidence and lifting of environment, which may move air from
!       levels with large-scale cloud to levels with no large-scale cloud
!
! Method:
! Three possible approaches were considered:
!
! 1. Ad-hoc #1 approach.  At each level, adjust dcondt for the activated
!    and unactivated portions of a particular aerosol species so that the
!    ratio of dcondt (activated/unactivate) is equal to the ratio of the
!    mixing ratios before convection.
!    THIS WAS IMPLEMENTED IN MIRAGE2
!
! 2. Ad-hoc #2 approach.  At each level, adjust dcondt for the activated
!    and unactivated portions of a particular aerosol species so that the
!    change to the activated portion is minimized (zero if possible).  The
!    would minimize effects of convection on the large-scale cloud.
!    THIS IS CURRENTLY IMPLEMENTED IN CAM5 where we assume that convective
!    clouds have no impact on the stratiform-cloudborne aerosol
!
! 3. Mechanistic approach that treats the details of interactions between
!    the large-scale and convective clouds.  (Something for the future.)
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

   implicit none

!-----------------------------------------------------------------------
! arguments
! (note:  TMR = tracer mixing ratio)
   integer,  intent(in)    :: pcnst_extd
   real(r8), intent(inout) :: dcondt(pcnst_extd,pver)
                              ! overall TMR tendency from convection
   real(r8), intent(inout) :: dcondt_resusp(pcnst_extd,pver)
                              ! portion of TMR tendency due to resuspension
                              ! (actually, due to the adjustments made here)
   real(r8), intent(in)    :: const(pcnst_extd,pver)  ! TMRs before convection

   real(r8), intent(in)    :: dp_i(pver) ! pressure thickness of level (in mb)
   integer,  intent(in)    :: ktop, kbot_prevap ! indices of top and bottom cloud levels

!-----------------------------------------------------------------------
! local variables
   integer  :: k, ll, la, lc, n, l
   real(r8) :: qa, qc, qac           ! working variables (mixing ratios)
   real(r8) :: qdota, qdotc, qdotac  ! working variables (MR tendencies)
!-----------------------------------------------------------------------

   do n = 1, nbins

      do ll = 1, nspec(n) + 2
         l = bin_idx(n,ll)
         la = l
         lc = l + ncnst_tot

         do k = ktop, kbot_prevap
            if (convproc_do_evaprain_atonce) then
               dcondt_resusp(la,k) = dcondt(la,k)
               dcondt_resusp(lc,k) = dcondt(lc,k)
            else
               qdota = dcondt(la,k)
               qdotc = dcondt(lc,k)
               qdotac = qdota + qdotc

               dcondt(la,k) = qdotac
               dcondt(lc,k) = 0.0_r8

               dcondt_resusp(la,k) = (dcondt(la,k) - qdota)
               dcondt_resusp(lc,k) = (dcondt(lc,k) - qdotc)
            end if
         end do

      end do   ! "ll = -1, nspec_amode(n)"
   end do      ! "n = 1, ntot_amode"

   return
   end subroutine ma_resuspend_convproc



!=========================================================================================



end module carma_aero_convproc
