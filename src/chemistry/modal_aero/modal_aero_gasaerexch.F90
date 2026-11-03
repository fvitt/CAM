! modal_aero_gasaerexch.F90


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_gasaerexch --- does modal aerosol gas-aerosol exchange
!
! !INTERFACE:
   module modal_aero_gasaerexch

! !USES:
  use shr_kind_mod,    only:  r8 => shr_kind_r8
  use chem_mods,       only:  gas_pcnst
  use ref_pres,        only:  top_lev => clim_modal_aero_top_lev
  use ppgrid,          only:  pcols, pver
  use aero_gasaerexch, only: soa_condevap_diag_t, col_tend_diag_t

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public modal_aero_gasaerexch_sub, modal_aero_gasaerexch_init

! !PUBLIC DATA MEMBERS:
  integer, parameter :: pcnstxx = gas_pcnst
  integer, protected, public :: maxspec_pcage != nspec_max

  integer, protected, public :: modefrm_pcage
  integer, protected, public :: nspecfrm_pcage
  integer :: modetoo_pcage

  integer, protected, allocatable, public :: lspecfrm_pcage(:)
  integer, protected, allocatable, public :: lspectoo_pcage(:)

  real(r8), parameter, public :: n_so4_monolayers_pcage = 8.0_r8

! number of so4(+nh4) monolayers needed to "age" a carbon particle

  real(r8), parameter, public :: &
              dr_so4_monolayers_pcage = n_so4_monolayers_pcage * 4.76e-10_r8
! thickness of the so4 monolayers (m)
! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3,
!    --> 1 mol so4(+nh4)  = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
! aging criterion is approximate so do not try to distinguish
!    sulfuric acid, bisulfate, ammonium sulfate

  real(r8), protected, allocatable, public :: soa_equivso4_factor(:)
! this factor converts an soa volume to a volume of so4(+nh4)
! having same hygroscopicity as the soa

  real (r8) :: fac_m2v_nh4, fac_m2v_so4
  real (r8), allocatable :: fac_m2v_soa(:)

  real (r8), allocatable :: fac_m2v_pcarbon(:)

! mode and species information, set from the aerosol properties object
  integer :: nmodes = 0                         ! number of modes
  integer :: nspec_max = 0                      ! maximum number of species in a mode
  integer, allocatable :: nspec_mode(:)         ! number of species in each mode
  integer, allocatable :: num_a_ndx(:)          ! constituent index of interstitial number in each mode
  integer, allocatable :: mass_a_ndx(:,:)       ! constituent index of interstitial species in each mode
  integer :: mode_pcarbon = -1                  ! primary carbon mode (<= 0 if none)
  integer :: mode_accum = -1                    ! accumulation mode (<= 0 if none)
  character(len=32), allocatable :: cw_cnst_name(:) ! cloud-borne name of each aerosol constituent

! numbers of soa and poa species (maximum over the modes),
! set from the aerosol properties object
  integer :: nsoa = 0
  integer :: npoa = 0

! constituent indices of interstitial aerosol species in each mode
! (set from the aerosol properties object; <= 0 if species is absent)
  integer, allocatable :: so4_a_ndx(:)      ! sulfate
  integer, allocatable :: nh4_a_ndx(:)      ! ammonium
  integer, allocatable :: soa_a_ndx(:,:)    ! secondary organic aerosol species
  integer, allocatable :: pom_a_ndx(:,:)    ! primary organic matter species

! constituent indices of the soa gas species (<= 0 if absent)
  integer, allocatable :: soag_cnst_ndx(:)

  real(r8), allocatable :: mw_soa_host(:)   ! molec wght of soa used in host code
  real(r8), allocatable :: mw_poa_host(:)   ! molec wght of poa used in host code

  type(soa_condevap_diag_t) :: soa_diag     ! SOA condensation/evaporation diagnostics

! column tendency diagnostics and their field index for each constituent (0 if none)
  type(col_tend_diag_t) :: sfgaex1_diag     ! gas-aerosol exchange ('_sfgaex1')
  type(col_tend_diag_t) :: sfgaex2_diag     ! renaming of interstitial species ('_sfgaex2')
  type(col_tend_diag_t) :: sfgaex2cw_diag   ! renaming of cloud-borne species ('_sfgaex2')
  integer, allocatable :: sfgaex1_fld(:)
  integer, allocatable :: sfgaex2_fld(:)
  integer, allocatable :: sfgaex2cw_fld(:)

! !DESCRIPTION: This module implements ...
!
! !REVISION HISTORY:
!
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! list private module data here

!EOC
!----------------------------------------------------------------------


  contains


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_gasaerexch_sub --- ...
!
! !INTERFACE:
subroutine modal_aero_gasaerexch_sub(                            &
                        aero_props,                              &
                        lchnk,    ncol,     nstep,               &
                        loffset,  deltat,                        &
                        t,        pmid,     pdel,                &
                        qh2o,               troplev,             &
                        q,                  qqcw,                &
                        dqdt_other,         dqqcwdt_other,       &
                        dgncur_a,           dgncur_awet,         &
                        sulfeq         )

! !USES:
use modal_aero_rename, only:  modal_aero_rename_sub
use aero_gasaerexch,   only:  aero_gasaerexch_uptkrates, aero_gasaerexch_soa_column, &
                              nh3_diffus_ratio, method_soa_reversible, soa_condevap_acc_t

use chem_mods,         only:  adv_mass
use constituents,      only:  pcnst, cnst_get_ind
use physconst,         only:  gravit, mwdry, rair
use cam_abortutils,    only:  endrun
use spmd_utils,        only:  iam, masterproc
use phys_control,      only:  cam_chempkg_is
use aerosol_properties_mod, only: aerosol_properties

implicit none

! !PARAMETERS:
   class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   integer,  intent(in)    :: troplev(pcols)       ! tropopause vertical index
   real(r8), intent(in)    :: deltat               ! time step (s)

   real(r8), intent(inout) :: q(ncol,pver,pcnstxx) ! tracer mixing ratio (TMR) array
                                                   ! *** MUST BE  #/kmol-air for number
                                                   ! *** MUST BE mol/mol-air for mass
                                                   ! *** NOTE ncol dimension
   real(r8), intent(inout) :: qqcw(ncol,pver,pcnstxx)
                                                   ! like q but for cloud-borner tracers
   real(r8), intent(in)    :: dqdt_other(ncol,pver,pcnstxx)
                                                   ! TMR tendency from other continuous
                                                   ! growth processes (aqchem, soa??)
                                                   ! *** NOTE ncol dimension
   real(r8), intent(in)    :: dqqcwdt_other(ncol,pver,pcnstxx)
                                                   ! like dqdt_other but for cloud-borner tracers
   real(r8), intent(in)    :: t(pcols,pver)        ! temperature at model levels (K)
   real(r8), intent(in)    :: pmid(pcols,pver)     ! pressure at model levels (Pa)
   real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: qh2o(pcols,pver)     ! water vapor mixing ratio (kg/kg)
   real(r8), intent(in)    :: dgncur_a(pcols,pver,nmodes)
   real(r8), intent(in)    :: dgncur_awet(pcols,pver,nmodes)
   real(r8), pointer       :: sulfeq(:,:,:)

                                 ! dry & wet geo. mean dia. (m) of number distrib.

! !DESCRIPTION:
! computes TMR (tracer mixing ratio) tendencies for gas condensation
!    onto aerosol particles
!
! this version does condensation of H2SO4, NH3, and MSA, both treated as
! completely non-volatile (gas --> aerosol, but no aerosol --> gas)
!    gas H2SO4 goes to aerosol SO4
!    gas MSA (if present) goes to aerosol SO4
!       aerosol MSA is not distinguished from aerosol SO4
!    gas NH3 (if present) goes to aerosol NH4
!       if gas NH3 is not present, then ????
!
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: jsrflx_gaexch = 1
   integer, parameter :: jsrflx_rename = 2
   integer, parameter :: method_soa = method_soa_reversible
!     method_soa=0 is no uptake
!     method_soa=1 is irreversible uptake done like h2so4 uptake
!     method_soa=2 is reversible uptake using subr aero_soaexch
!     (see aero_gasaerexch)

   integer :: i, iq
   integer :: ido_so4a(nmodes), ido_nh4a(nmodes)
   integer ::  ido_soaa(nmodes,nsoa)
   integer :: j, jsrf, jsoa
   integer :: k
   integer :: l, lb, lsfrm, lstoo
   integer :: l_so4g, l_nh4g, l_msag
   integer :: l_soag(nsoa)
   integer :: n, ntot_soamode

   logical :: is_dorename_atik, dorename_atik(ncol,pver)

   real (r8) :: avg_uprt_nh4, avg_uprt_so4
   real (r8) :: deltatxx
   real (r8) :: dqdt_nh4(nmodes), dqdt_so4(nmodes)
   real (r8) :: dqdt_soa(nmodes,nsoa)
   real (r8) :: fac_volsfc_pcarbon
   real (r8) :: fgain_nh4(nmodes), fgain_so4(nmodes)
   real (r8) :: pdel_fac
   real (r8) :: qmax_nh4, qnew_nh4, qnew_so4
   real (r8) :: qold_nh4(nmodes), qold_so4(nmodes)
   real (r8) :: soag_col(pcols,pver,nsoa)          ! soa gas mixing ratios
   real (r8) :: soa_col(pcols,pver,nmodes,nsoa)    ! soa aerosol mixing ratios
   real (r8) :: poa_col(pcols,pver,nmodes,npoa)    ! poa aerosol mixing ratios
   real (r8) :: dqdt_soa_col(pcols,pver,nmodes,nsoa) ! soa aerosol tendencies
   real (r8) :: sum_dqdt_soa_col(pcols,pver,nsoa)  ! soa gas uptake (negative of gas tendency)
   real (r8) :: sum_dqdt_msa, sum_dqdt_so4
   real (r8) :: sum_dqdt_soa(nsoa)
   real (r8) :: sum_dqdt_nh4, sum_dqdt_nh4_b
   real (r8) :: sum_uprt_nh4, sum_uprt_so4
   real (r8) :: tmp1, tmp2, tmpa
   real (r8) :: tmp_kxt, tmp_pxt
   real (r8) :: tmp_so4a_bgn, tmp_so4a_end
   real (r8) :: tmp_so4g_avg, tmp_so4g_bgn, tmp_so4g_equ
   real (r8) :: uptkrate(pcols,pver,nmodes) ! gas-to-aerosol mass transfer rates (1/s)
   real (r8) :: numdens(pcols,pver,nmodes)  ! mode number density (#/m3)
   real (r8) :: rhoair, aircon
   real (r8) :: uptkratebb(nmodes)
   real (r8) :: vol_core, vol_shell
   real (r8) :: xferfrac_pcage, xferfrac_max
   real (r8) :: xferrate

   logical  :: do_msag         ! true if msa gas is a species
   logical  :: do_nh4g         ! true if nh3 gas is a species
   logical  :: do_soag_any         ! true if soa gas is a species
   logical  :: do_soag(nsoa)       ! true if soa gas is a species
   logical  :: do_soamode(nmodes) ! true if mode takes part in soa exchange
   logical  :: do_soa(nmodes,nsoa) ! true if soa species takes up soa gas in mode
   logical  :: do_poa(nmodes)      ! true if poa is included in the soa exchange of a mode

   logical  :: dotend(pcnstxx)          ! identifies species directly involved in
                                        !    gas-aerosol exchange (gas condensation)
   logical  :: dotendqqcw(pcnstxx)      ! like dotend but for cloud-borner tracers
   logical  :: dotendrn(pcnstxx), dotendqqcwrn(pcnstxx)
                                        ! identifies species involved in renaming
                                        !    after "continuous growth"
                                        !    (gas-aerosol exchange and aqchem)

   integer, parameter :: nsrflx = 2     ! last dimension of qsrflx
   real(r8) :: dqdt(ncol,pver,pcnstxx)  ! TMR "delta q" array - NOTE dims
   real(r8) :: dqqcwdt(ncol,pver,pcnstxx) ! like dqdt but for cloud-borner tracers
   real(r8) :: qsrflx(pcols,pcnstxx,nsrflx)
                              ! process-specific column tracer tendencies
                              ! (1=renaming, 2=gas condensation)
   type(soa_condevap_acc_t) :: condevap   ! SOA condensation/evaporation diagnostics
   logical  :: condevap_mask(nmodes,nsoa)
   real(r8) :: condevap_fac(nmodes,nsoa)

   real(r8) :: qqcwsrflx(pcols,pcnstxx,nsrflx)

!----------------------------------------------------------------------

! set gas species indices
   call cnst_get_ind( 'H2SO4', l_so4g, .false. )
   call cnst_get_ind( 'NH3',   l_nh4g, .false. )
   if ( .not. cam_chempkg_is('geoschem_mam4') ) then
      call cnst_get_ind( 'MSA',   l_msag, .false. )
   else
      l_msag = 0
   endif
   l_so4g = l_so4g - loffset
   l_nh4g = l_nh4g - loffset
   l_msag = l_msag - loffset
   if ((l_so4g <= 0) .or. (l_so4g > pcnstxx)) then
      write( *, '(/a/a,2i7)' )   &
         '*** modal_aero_gasaerexch_sub -- cannot find H2SO4 species',   &
         '    l_so4g, loffset =', l_so4g, loffset
      call endrun( 'modal_aero_gasaerexch_sub error' )
   end if
   do_nh4g = .false.
   do_msag = .false.
   if ((l_nh4g > 0) .and. (l_nh4g <= pcnstxx)) do_nh4g = .true.
   if ((l_msag > 0) .and. (l_msag <= pcnstxx)) do_msag = .true.

   do_soag_any = .false.
   do_soag(:) = .false.
   do jsoa = 1, nsoa
      l_soag(jsoa) = soag_cnst_ndx(jsoa) - loffset
      if ((method_soa == 1) .or. (method_soa == 2)) then
         if ((l_soag(jsoa) > 0) .and. (l_soag(jsoa) <= pcnstxx)) then
            do_soag_any = .true.
            do_soag(jsoa) = .true.
         end if
      else if (method_soa /= 0) then
         write(*,'(/a,1x,i10)') '*** modal_aero_gasaerexch_sub - bad method_soa =', method_soa
         call endrun( 'modal_aero_gasaerexch_sub error' )
      end if
   end do ! jsoa

! set tendency flags
   dotend(:) = .false.
   dotendqqcw(:) = .false.
   ido_so4a(:) = 0
   ido_nh4a(:) = 0
   ido_soaa(:,:) = 0

   dotend(l_so4g) = .true.
   if ( do_nh4g ) dotend(l_nh4g) = .true.
   if ( do_msag ) dotend(l_msag) = .true.
   do jsoa = 1, nsoa
      if ( do_soag(jsoa) ) dotend(l_soag(jsoa)) = .true.
   end do

   ntot_soamode = 0
   do n = 1, nmodes
      l = so4_a_ndx(n)-loffset
      if ((l > 0) .and. (l <= pcnstxx)) then
         dotend(l) = .true.
         ido_so4a(n) = 1
         if ( do_nh4g ) then
            l = nh4_a_ndx(n)-loffset
            if ((l > 0) .and. (l <= pcnstxx)) then
               dotend(l) = .true.
               ido_nh4a(n) = 1
            end if
         end if
      end if

      do jsoa = 1, nsoa
         if ( do_soag(jsoa) ) then
            l = soa_a_ndx(n,jsoa)-loffset
            if ((l > 0) .and. (l <= pcnstxx)) then
               dotend(l) = .true.
               ido_soaa(n,jsoa) = 1
               ntot_soamode = n
            end if
         end if
      end do ! jsoa
   end do ! n


   if ( do_soag_any ) ntot_soamode = max( ntot_soamode, modefrm_pcage )
   do n = 1, nmodes
      do_soamode(n) = (n <= ntot_soamode)
   end do

   if (modefrm_pcage > 0) then
      ido_so4a(modefrm_pcage) = 2
      if (ido_nh4a(modetoo_pcage) == 1) ido_nh4a(modefrm_pcage) = 2
      do jsoa = 1, nsoa
         if (ido_soaa(modetoo_pcage,jsoa) == 1) ido_soaa(modefrm_pcage,jsoa) = 2
      end do
      do iq = 1, nspecfrm_pcage
         lsfrm = lspecfrm_pcage(iq)-loffset
         lstoo = lspectoo_pcage(iq)-loffset
         if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
            dotend(lsfrm) = .true.
            if ((lstoo > 0) .and. (lstoo <= pcnst)) then
               dotend(lstoo) = .true.
            end if
         end if
      end do


      n = mode_pcarbon
      fac_volsfc_pcarbon = exp( 2.5_r8*(aero_props%alogsig(n)**2) )
      xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps
   end if


! zero out tendencies and other
   dqdt(:,:,:) = 0.0_r8
   dqqcwdt(:,:,:) = 0.0_r8
   qsrflx(:,:,:) = 0.0_r8
   qqcwsrflx(:,:,:) = 0.0_r8

!-------Initialize evap/cond diagnostics (ncols x pver)-----------
   call soa_diag%reset(condevap)
   do jsoa = 1, nsoa
      do n = 1, nmodes
         condevap_mask(n,jsoa) = do_soag(jsoa) .and. (ido_soaa(n,jsoa) == 1)
         if (condevap_mask(n,jsoa)) then
            condevap_fac(n,jsoa) = adv_mass(soa_a_ndx(n,jsoa)-loffset)/mwdry
         else
            condevap_fac(n,jsoa) = 0.0_r8
         end if
      end do
   end do
!---------------------------------------------------

! compute gas-to-aerosol mass transfer rates
   numdens(:,:,:) = 0.0_r8
   do n = 1, nmodes
      do k = top_lev, pver
         do i = 1, ncol
            rhoair = pmid(i,k)/(rair*t(i,k))   ! (kg-air/m3)
!   number conc. (#/m3) -- note q(i,k,numptr) is (#/kmol-air)
!   so need aircon in (kmol-air/m3)
            aircon = rhoair/mwdry              ! (kmol-air/m3)
            numdens(i,k,n) = q(i,k,num_a_ndx(n)-loffset)*aircon
         end do
      end do
   end do
   call aero_gasaerexch_uptkrates( aero_props, ncol, pver, top_lev, numdens, dgncur_awet, &
                                   t, pmid, uptkrate )

! modes and species taking up soa gas
   do jsoa = 1, nsoa
      do n = 1, nmodes
         do_soa(n,jsoa) = (ido_soaa(n,jsoa) > 0)
      end do
   end do
! in previous code versions with nsoa=1,
!    qold_poa was non-zero (i.e., loaded from q) only when ido_soaa(n)=1
! thus qold_poa=0 for the primary carbon mode which has ido_soaa=2
! this is probably not how it should be
   do n = 1, nmodes
      do_poa(n) = any(ido_soaa(n,:) == 1)
   end do

! compute TMR tendencies for soag and soa interstial aerosol in all grid cells
   soag_col(:,:,:) = 0.0_r8
   soa_col(:,:,:,:) = 0.0_r8
   poa_col(:,:,:,:) = 0.0_r8
   if ( do_soag_any ) then
      do k = top_lev, pver
         do i = 1, ncol
            do jsoa = 1, nsoa
               if ( do_soag(jsoa) ) soag_col(i,k,jsoa) = q(i,k,l_soag(jsoa))
            end do
            do n = 1, nmodes
               do jsoa = 1, nsoa
                  if (ido_soaa(n,jsoa) == 1) then
                     soa_col(i,k,n,jsoa) = q(i,k,soa_a_ndx(n,jsoa)-loffset)
                  end if
               end do
               if ( do_poa(n) ) then
                  do j = 1, npoa
                     l = pom_a_ndx(n,j)-loffset
                     if (l > 0) poa_col(i,k,n,j) = q(i,k,l)
                  end do
               end if
            end do
         end do
      end do

      call aero_gasaerexch_soa_column( method_soa, ncol, pver, top_lev, deltat, t, pmid, &
           nmodes, nsoa, npoa, do_soamode, do_soa, do_soag, uptkrate, &
           mw_poa_host, mw_soa_host, soag_col, soa_col, poa_col, &
           dqdt_soa_col, sum_dqdt_soa_col )
   else ! no soa gas, no uptake
      dqdt_soa_col(:,:,:,:) = 0.0_r8
      sum_dqdt_soa_col(:,:,:) = 0.0_r8
   end if


! use this for tendency calcs to avoid generating very small negative values
   deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)


   jsrf = jsrflx_gaexch
   do k=top_lev,pver
      do i=1,ncol

!   fgain_so4(n) = fraction of total h2so4 uptake going to mode n
!   fgain_nh4(n) = fraction of total  nh3  uptake going to mode n
         sum_uprt_so4 = 0.0_r8
         sum_uprt_nh4 = 0.0_r8
         do n = 1, nmodes
            uptkratebb(n) = uptkrate(i,k,n)
            if (ido_so4a(n) > 0) then
               fgain_so4(n) = uptkratebb(n)
               sum_uprt_so4 = sum_uprt_so4 + fgain_so4(n)
               if (ido_so4a(n) == 1) then
                  qold_so4(n) = q(i,k,so4_a_ndx(n)-loffset)
               else
                  qold_so4(n) = 0.0_r8
               end if
            else
               fgain_so4(n) = 0.0_r8
               qold_so4(n) = 0.0_r8
            end if

            if (ido_nh4a(n) > 0) then
               !   nh3_diffus_ratio is for gas diffusivity (nh3/h2so4)
               !   differences in fuch-sutugin and accom coef ignored
               fgain_nh4(n) = uptkratebb(n)*nh3_diffus_ratio
               sum_uprt_nh4 = sum_uprt_nh4 + fgain_nh4(n)
               if (ido_nh4a(n) == 1) then
                  qold_nh4(n) = q(i,k,nh4_a_ndx(n)-loffset)
               else
                  qold_nh4(n) = 0.0_r8
               end if
            else
               fgain_nh4(n) = 0.0_r8
               qold_nh4(n) = 0.0_r8
            end if

         end do ! n

         if (sum_uprt_so4 > 0.0_r8) then
            do n = 1, nmodes
               fgain_so4(n) = fgain_so4(n) / sum_uprt_so4
            end do
         end if
!       at this point (sum_uprt_so4 <= 0.0) only when all the fgain_so4 are zero
         if (sum_uprt_nh4 > 0.0_r8) then
            do n = 1, nmodes
               fgain_nh4(n) = fgain_nh4(n) / sum_uprt_nh4
            end do
         end if

!   uptake amount (fraction of gas uptaken) over deltat
         avg_uprt_so4 = (1.0_r8 - exp(-deltatxx*sum_uprt_so4))/deltatxx
         avg_uprt_nh4 = (1.0_r8 - exp(-deltatxx*sum_uprt_nh4))/deltatxx

!   sum_dqdt_so4 = so4_a tendency from h2so4 gas uptake (mol/mol/s)
!   sum_dqdt_msa = msa_a tendency from msa   gas uptake (mol/mol/s)
!   sum_dqdt_nh4 = nh4_a tendency from nh3   gas uptake (mol/mol/s)
         sum_dqdt_so4 = q(i,k,l_so4g) * avg_uprt_so4
         if ( do_msag ) then
            sum_dqdt_msa = q(i,k,l_msag) * avg_uprt_so4
         else
            sum_dqdt_msa = 0.0_r8
         end if
         if ( do_nh4g ) then
            sum_dqdt_nh4 = q(i,k,l_nh4g) * avg_uprt_nh4
         else
            sum_dqdt_nh4 = 0.0_r8
         end if

         if ( associated(sulfeq) .and. (k <= troplev(i)) ) then
            !   compute TMR tendencies for so4 interstial aerosol due to reversible gas uptake
            !   only above the tropopause

            tmp_kxt = deltatxx*sum_uprt_so4  ! sum over modes of uptake_rate*deltat
            tmp_pxt = 0.0_r8
            do n = 1, nmodes
               if (ido_so4a(n) <= 0) cycle
               tmp_pxt = tmp_pxt + uptkratebb(n)*sulfeq(i,k,n)
            end do
            tmp_pxt = max( 0.0_r8, tmp_pxt*deltatxx )  ! sum over modes of uptake_rate*sulfeq*deltat
            tmp_so4g_bgn = q(i,k,l_so4g)
            ! calc avg h2so4(g) over deltat
            if (tmp_kxt >= 1.0e-5_r8) then
               ! exponential decay towards equilibrium value solution
               tmp_so4g_equ = tmp_pxt/tmp_kxt
               tmp_so4g_avg = tmp_so4g_equ + (tmp_so4g_bgn-tmp_so4g_equ)*(1.0_r8-exp(-tmp_kxt))/tmp_kxt
            else
               ! first order approx for tmp_kxt small
               tmp_so4g_avg = tmp_so4g_bgn*(1.0_r8-0.5_r8*tmp_kxt) + 0.5_r8*tmp_pxt
            end if
            sum_dqdt_so4 = 0.0_r8
            do n = 1, nmodes
               if (ido_so4a(n) <= 0) cycle
               ! calc change to so4(a) in mode n
               if (ido_so4a(n) == 1) then
                  l = so4_a_ndx(n)-loffset
                  tmp_so4a_bgn = q(i,k,l)
               else
                  tmp_so4a_bgn = 0.0_r8
               end if
               tmp_so4a_end = tmp_so4a_bgn + deltatxx*uptkratebb(n)*(tmp_so4g_avg-sulfeq(i,k,n))
               tmp_so4a_end = max( 0.0_r8, tmp_so4a_end )
               dqdt_so4(n) = (tmp_so4a_end - tmp_so4a_bgn)/deltatxx
               sum_dqdt_so4 = sum_dqdt_so4 + dqdt_so4(n)
            end do
            ! do not allow msa condensation in stratosphere
            ! ( Note that the code for msa has never been used.
            !   The plan was to simulate msa(g), treat it as non-volatile (like h2so4(g)),
            !   and treat condensed msa as sulfate, so just one additional tracer. )
            if ( do_msag ) sum_dqdt_msa = 0.0_r8

         else
            !   compute TMR tendencies for so4 interstial aerosol due to simple gas uptake
            do n = 1, nmodes
               dqdt_so4(n) = fgain_so4(n)*(sum_dqdt_so4 + sum_dqdt_msa)
            end do
         end if

         !   compute TMR tendencies for nh4 interstial aerosol due to simple gas uptake
         !   but force nh4/so4 molar ratio <= 2
         sum_dqdt_nh4_b = 0.0_r8
         dqdt_nh4(:) = 0._r8
         if ( do_nh4g ) then
            do n = 1, nmodes
               dqdt_nh4(n) = fgain_nh4(n)*sum_dqdt_nh4
               qnew_nh4 = qold_nh4(n) + dqdt_nh4(n)*deltat
               qnew_so4 = qold_so4(n) + dqdt_so4(n)*deltat
               qmax_nh4 = 2.0_r8*qnew_so4
               if (qnew_nh4 > qmax_nh4) then
                  dqdt_nh4(n) = (qmax_nh4 - qold_nh4(n))/deltatxx
               end if
               sum_dqdt_nh4_b = sum_dqdt_nh4_b + dqdt_nh4(n)
            end do
         end if

!   TMR tendencies for soag and soa interstial aerosol (computed above)
!   sum_dqdt_soa = soa_a tendency from soa   gas uptake (mol/mol/s)
         dqdt_soa(:,:) = dqdt_soa_col(i,k,:,:)
         sum_dqdt_soa(:) = sum_dqdt_soa_col(i,k,:)

         pdel_fac = pdel(i,k)/gravit
         do n = 1, nmodes
            if (ido_so4a(n) == 1) then
               l = so4_a_ndx(n)-loffset
               dqdt(i,k,l) = dqdt_so4(n)
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(n)*pdel_fac
            end if

            if ( do_nh4g ) then
               if (ido_nh4a(n) == 1) then
                  l = nh4_a_ndx(n)-loffset
                  dqdt(i,k,l) = dqdt_nh4(n)
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(n)*pdel_fac
               end if
            end if

            do jsoa = 1, nsoa
               if ( do_soag(jsoa) ) then
                  if (ido_soaa(n,jsoa) == 1) then
                     l = soa_a_ndx(n,jsoa)-loffset
                     dqdt(i,k,l) = dqdt_soa(n,jsoa) !calculated by  aero_soaexch for method_soa=2
                     qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(n,jsoa)*pdel_fac
                  end if
               end if
            end do
         end do ! n

!------- condensation/evaporation diagnostics
         call soa_diag%accum( condevap, i, k, nmodes, condevap_mask, condevap_fac, dqdt_soa )

!   compute TMR tendencies for h2so4, nh3, and msa gas
!   due to simple gas uptake
         l = l_so4g
         dqdt(i,k,l) = -sum_dqdt_so4
         qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac

         if ( do_msag ) then
            l = l_msag
            dqdt(i,k,l) = -sum_dqdt_msa
            qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
         end if

         if ( do_nh4g ) then
            l = l_nh4g
            dqdt(i,k,l) = -sum_dqdt_nh4_b
            qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
         end if

         do jsoa = 1, nsoa
            if ( do_soag(jsoa) ) then
               l = l_soag(jsoa)
               dqdt(i,k,l) = -sum_dqdt_soa(jsoa)
! dqdt for gas is negative of the sum of dqdt for aerosol soa species in each mode: Manish
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
            end if
         end do

!   compute TMR tendencies associated with primary carbon aging
         if (modefrm_pcage > 0) then
            n = mode_pcarbon
            tmpa = 0.0_r8
            do jsoa = 1, nsoa
               tmpa = tmpa + dqdt_soa(n,jsoa)*fac_m2v_soa(jsoa)*soa_equivso4_factor(jsoa)
            end do
            vol_shell = deltat *   &
                 ( dqdt_so4(n)*fac_m2v_so4 + dqdt_nh4(n)*fac_m2v_nh4 + tmpa )
            vol_core = 0.0_r8
            do l = 1, nspec_mode(n)
               vol_core = vol_core + &
                    q(i,k,mass_a_ndx(l,n)-loffset)*fac_m2v_pcarbon(l)
            end do
!   ratio1 = vol_shell/vol_core =
!      actual hygroscopic-shell-volume/carbon-core-volume after gas uptake
!   ratio2 = 6.0_r8*dr_so4_monolayers_pcage/(dgncur_a*fac_volsfc_pcarbon)
!      = (shell-volume corresponding to n_so4_monolayers_pcage)/core-volume
!      The 6.0/(dgncur_a*fac_volsfc_pcarbon) = (mode-surface-area/mode-volume)
!   Note that vol_shell includes both so4+nh4 AND soa as "equivalent so4",
!      The soa_equivso4_factor accounts for the lower hygroscopicity of soa.
!
!   Define xferfrac_pcage = min( 1.0, ratio1/ratio2)
!   But ratio1/ratio2 == tmp1/tmp2, and coding below avoids possible overflow
!
            tmp1 = vol_shell*dgncur_a(i,k,n)*fac_volsfc_pcarbon
            tmp2 = max( 6.0_r8*dr_so4_monolayers_pcage*vol_core, 0.0_r8 )
            if (tmp1 >= tmp2) then
               xferfrac_pcage = xferfrac_max
            else
               xferfrac_pcage = min( tmp1/tmp2, xferfrac_max )
            end if

            if (xferfrac_pcage > 0.0_r8) then
               do iq = 1, nspecfrm_pcage
                  lsfrm = lspecfrm_pcage(iq)-loffset
                  lstoo = lspectoo_pcage(iq)-loffset
                  xferrate = (xferfrac_pcage/deltat)*q(i,k,lsfrm)
                  dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xferrate
                  qsrflx(i,lsfrm,jsrf) = qsrflx(i,lsfrm,jsrf) - xferrate*pdel_fac
                  if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                     dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xferrate
                     qsrflx(i,lstoo,jsrf) = qsrflx(i,lstoo,jsrf) + xferrate*pdel_fac
                  end if
               end do

               if (ido_so4a(modetoo_pcage) > 0) then
                  l = so4_a_ndx(modetoo_pcage)-loffset
                  dqdt(i,k,l) = dqdt(i,k,l) + dqdt_so4(modefrm_pcage)
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(modefrm_pcage)*pdel_fac
               end if

               if (ido_nh4a(modetoo_pcage) > 0) then
                  l = nh4_a_ndx(modetoo_pcage)-loffset
                  dqdt(i,k,l) = dqdt(i,k,l) + dqdt_nh4(modefrm_pcage)
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(modefrm_pcage)*pdel_fac
               end if

               do jsoa = 1, nsoa
                  if (ido_soaa(modetoo_pcage,jsoa) > 0) then
                     l = soa_a_ndx(modetoo_pcage,jsoa)-loffset
                     dqdt(i,k,l) = dqdt(i,k,l) + dqdt_soa(modefrm_pcage,jsoa)
                     qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(modefrm_pcage,jsoa)*pdel_fac
                  end if
               end do

            end if

         end if

      end do   ! "i = 1, ncol"
   end do     ! "k = top_lev, pver"

!
! do renaming calcs
!
   dotendrn(:) = .false.
   dotendqqcwrn(:) = .false.
   dorename_atik(1:ncol,:) = .true.
   is_dorename_atik = .true.
   call modal_aero_rename_sub(                              &
        'modal_aero_gasaerexch_sub',            &
        lchnk,             ncol,      nstep,    &
        loffset,           deltat,              &
        pdel,              troplev,             &
        dotendrn,          q,                   &
        dqdt,              dqdt_other,          &
        dotendqqcwrn,      qqcw,                &
        dqqcwdt,           dqqcwdt_other,       &
        is_dorename_atik,  dorename_atik,       &
        jsrflx_rename,     nsrflx,              &
        qsrflx,            qqcwsrflx            )


!  This applies dqdt tendencies for all species
!  apply the dqdt to update q (and same for qqcw)
!
   do l = 1, pcnstxx
      if ( dotend(l) .or. dotendrn(l) ) then
         do k = top_lev, pver
            do i = 1, ncol
               q(i,k,l) = q(i,k,l) + dqdt(i,k,l)*deltat
            end do
         end do
      end if
      if ( dotendqqcw(l) .or. dotendqqcwrn(l) ) then
         do k = top_lev, pver
            do i = 1, ncol
               qqcw(i,k,l) = qqcw(i,k,l) + dqqcwdt(i,k,l)*deltat
            end do
         end do
      end if
   end do


!-----Outfld for condensation/evaporation------------------------------
   call soa_diag%outfld( condevap, lchnk )
!-----------------------------------------------------------------------

!   do history file column-tendency fields
!   (lb = l + loffset is a constituent index only for species with tendencies,
!    so the field index arrays are only accessed for those)
   do l = 1, pcnstxx
      lb = l + loffset
      if ( dotend(l) ) then
         if ( sfgaex1_fld(lb) > 0 ) then
            call sfgaex1_diag%outfld( sfgaex1_fld(lb), ncol, lchnk, adv_mass(l)/mwdry, &
                                      qsrflx(:,l,jsrflx_gaexch) )
         end if
      end if
      if ( dotendrn(l) ) then
         if ( sfgaex2_fld(lb) > 0 ) then
            call sfgaex2_diag%outfld( sfgaex2_fld(lb), ncol, lchnk, adv_mass(l)/mwdry, &
                                      qsrflx(:,l,jsrflx_rename) )
         end if
      end if
      if ( dotendqqcwrn(l) ) then
         if ( sfgaex2cw_fld(lb) > 0 ) then
            call sfgaex2cw_diag%outfld( sfgaex2cw_fld(lb), ncol, lchnk, adv_mass(l)/mwdry, &
                                        qqcwsrflx(:,l,jsrflx_rename) )
         end if
      end if
   end do ! l = ...

   return
   end subroutine modal_aero_gasaerexch_sub


!----------------------------------------------------------------------

      subroutine modal_aero_gasaerexch_init(aero_props)

!-----------------------------------------------------------------------
!
! Purpose:
!    set do_adjust and do_aitken flags
!    create history fields for column tendencies associated with
!       modal_aero_calcsize
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

use modal_aero_rename, only: npair_renamexf, nspecfrm_renamexf, lspecfrma_renamexf, lspecfrmc_renamexf, &
                             lspectooa_renamexf, lspectooc_renamexf

use cam_abortutils, only: endrun
use constituents,   only: pcnst, cnst_get_ind, cnst_name
use spmd_utils,     only: masterproc
use phys_control,   only: phys_getopts
use aerosol_properties_mod, only: aerosol_properties

implicit none

!-----------------------------------------------------------------------
! arguments
   class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object

!-----------------------------------------------------------------------
! local
   integer  :: ipair, iq, iqfrm, iqtoo
   integer  :: jsoa
   integer  :: l, l1, l2, lsfrm, lstoo, lunout
   integer  :: l_so4g, l_nh4g, l_msag
   integer  :: mfrm, mtoo
   integer  :: n
   integer  :: nchfrm, nchfrmskip, nchtoo, nchtooskip, nspec

   logical  :: do_msag, do_nh4g
   logical  :: do_soag_any
   logical, allocatable :: do_soag(:)
   logical  :: dotend(pcnst), dotendqqcw(pcnst)

   real(r8) :: tmp2
   real(r8) :: spec_mw, spec_dens, spec_hygro   ! species molecular weight, density, hygroscopicity

   character(128)                 :: msg

   logical                        :: history_aerosol      ! Output the MAM aerosol tendencies
   !-----------------------------------------------------------------------

        call phys_getopts( history_aerosol_out        = history_aerosol   )

        ! mode tables, species indices and numbers of soa/poa species (nsoa, npoa)
        call set_mode_tables(aero_props)
        call set_species_ndx(aero_props)

        maxspec_pcage = nspec_max
        allocate(lspecfrm_pcage(maxspec_pcage))
        allocate(lspectoo_pcage(maxspec_pcage))
        allocate(soa_equivso4_factor(nsoa))
        allocate(fac_m2v_soa(nsoa))
        allocate(fac_m2v_pcarbon(nspec_max))
        allocate(do_soag(nsoa))
	lunout = 6
!
!   define "from mode" and "to mode" for primary carbon aging
!
!   skip (turn off) aging if either is absent,
!      or if accum mode so4 is absent
!
	modefrm_pcage = -999888777
	modetoo_pcage = -999888777
        if ((mode_pcarbon <= 0) .or. (mode_accum <= 0)) goto 15000
        l = so4_a_ndx(mode_accum)
	if ((l < 1) .or. (l > pcnst)) goto 15000

        modefrm_pcage = mode_pcarbon
        modetoo_pcage = mode_accum

!
!   define species involved in each primary carbon aging pairing
!	(include aerosol water)
!
!
	mfrm = modefrm_pcage
	mtoo = modetoo_pcage

	if (mfrm < 10) then
	    nchfrmskip = 1
	else if (mfrm < 100) then
	    nchfrmskip = 2
	else
	    nchfrmskip = 3
	end if
	if (mtoo < 10) then
	    nchtooskip = 1
	else if (mtoo < 100) then
	    nchtooskip = 2
	else
	    nchtooskip = 3
	end if
	nspec = 0

aa_iqfrm: do iqfrm = -1, nspec_mode(mfrm)

	    if (iqfrm == -1) then
                lsfrm = num_a_ndx(mfrm)
                lstoo = num_a_ndx(mtoo)
	    else if (iqfrm == 0) then
!   bypass transfer of aerosol water due to primary-carbon aging
		cycle aa_iqfrm
!               lsfrm = lwaterptr_amode(mfrm)
!               lstoo = lwaterptr_amode(mtoo)
	    else
                lsfrm = mass_a_ndx(iqfrm,mfrm)
		lstoo = 0
	    end if
	    if ((lsfrm < 1) .or. (lsfrm > pcnst)) cycle aa_iqfrm

	    if (lsfrm>0 .and. iqfrm>0 ) then
		nchfrm = len( trim( cnst_name(lsfrm) ) ) - nchfrmskip

! find "too" species having same lspectype_amode as the "frm" species
! AND same cnst_name (except for last 1/2/3 characters which are the mode index)
                do iqtoo = 1, nspec_mode(mtoo)
!		    if ( lspectype_amode(iqtoo,mtoo) .eq.   &
!		         lspectype_amode(iqfrm,mfrm) ) then
                        lstoo = mass_a_ndx(iqtoo,mtoo)
			nchtoo = len( trim( cnst_name(lstoo) ) ) - nchtooskip
			if (cnst_name(lsfrm)(1:nchfrm) == cnst_name(lstoo)(1:nchtoo)) then
			    exit
			else
			    lstoo = 0
			end if
!		    end if
		end do
	    end if

	    if ((lstoo < 1) .or. (lstoo > pcnst)) lstoo = 0
	    nspec = nspec + 1
	    lspecfrm_pcage(nspec) = lsfrm
	    lspectoo_pcage(nspec) = lstoo
	end do aa_iqfrm

	nspecfrm_pcage = nspec

!
!   output results
!
	if ( masterproc ) then

	write(lunout,9310)

	  mfrm = modefrm_pcage
	  mtoo = modetoo_pcage
	  write(lunout,9320) 1, mfrm, mtoo

	  do iq = 1, nspecfrm_pcage
	    lsfrm = lspecfrm_pcage(iq)
	    lstoo = lspectoo_pcage(iq)
	    if (lstoo .gt. 0) then
		write(lunout,9330) lsfrm, cnst_name(lsfrm),   &
      			lstoo, cnst_name(lstoo)
	    else
		write(lunout,9340) lsfrm, cnst_name(lsfrm)
	    end if
	  end do

	write(lunout,*)

	end if ! ( masterproc )

9310	format( / 'subr. modal_aero_gasaerexch_init - primary carbon aging pointers' )
9320	format( 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3 )
9330	format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340	format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )


15000 continue

! set tendency flags and gas species indices and flags
      dotend(:) = .false.

      call cnst_get_ind( 'H2SO4', l_so4g, .false. )
      if ((l_so4g <= 0) .or. (l_so4g > pcnst)) then
         write( *, '(/a/a,2i7)' )   &
            '*** modal_aero_gasaerexch_init -- cannot find H2SO4 species',   &
            '    l_so4g=', l_so4g
         call endrun( 'modal_aero_gasaerexch_init error' )
      end if
      dotend(l_so4g) = .true.

      call cnst_get_ind( 'NH3',   l_nh4g, .false. )
      do_nh4g = .false.
      if ((l_nh4g > 0) .and. (l_nh4g <= pcnst)) then
         do_nh4g = .true.
         dotend(l_nh4g) = .true.
      end if

      call cnst_get_ind( 'MSA',   l_msag, .false. )
      do_msag = .false.
      if ((l_msag > 0) .and. (l_msag <= pcnst)) then
         do_msag = .true.
         dotend(l_msag) = .true.
      end if

      do_soag_any = .false.
      do_soag(:) = .false.
      do jsoa = 1, nsoa
         l = soag_cnst_ndx(jsoa)
         if ((l > 0) .and. (l <= pcnst)) then
            do_soag_any = .true.
            do_soag(jsoa) = .true.
            dotend(l) = .true.
         end if
      end do


      do n = 1, nmodes
         l = so4_a_ndx(n)
         if ((l > 0) .and. (l <= pcnst)) then
            dotend(l) = .true.
            if ( do_nh4g ) then
               l = nh4_a_ndx(n)
               if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
            end if
         end if
         do jsoa = 1, nsoa
            if ( do_soag(jsoa) ) then
               l = soa_a_ndx(n,jsoa)
               if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
            end if
         end do
      end do

      if (modefrm_pcage > 0) then
         do iq = 1, nspecfrm_pcage
            lsfrm = lspecfrm_pcage(iq)
            lstoo = lspectoo_pcage(iq)
            if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotend(lsfrm) = .true.
               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                  dotend(lstoo) = .true.
               end if
            end if
         end do
      end if

!---------define history fields for new cond/evap diagnostics----------------------------------------
!   nsoa=5:  totals over all SOA species
!   nsoa=15: fossil fuel (SOA species 1-5), biomass burning (6-10) and biogenic (11-15) groups
      if (nsoa == 5) then
         call soa_diag%init( nsoa, do_total=.true., do_spec=.false., history_aerosol=history_aerosol )
      else if (nsoa == 15) then
         call soa_diag%init( nsoa, do_total=.false., do_spec=.false., history_aerosol=history_aerosol, &
                             grp=(/ (1+(jsoa-1)/5, jsoa=1,nsoa) /), &
                             grp_sfx=(/ 'ff', 'bb', 'bg' /), &
                             grp_desc=(/ 'Fossil  ', 'Biomass ', 'Biogenic' /) )
      else
         call soa_diag%init( nsoa, do_total=.false., do_spec=.false., history_aerosol=history_aerosol )
      end if
!------------------------------------------------------------------------------

!  define history fields for basic gas-aer exchange
!  and primary carbon aging from that
      call init_col_tend( sfgaex1_diag, sfgaex1_fld, dotend, cnst_name, '_sfgaex1', &
                          'gas-aerosol-exchange primary column tendency', history_aerosol )
!  define history fields for aitken-->accum renaming
      dotend(:) = .false.
      dotendqqcw(:) = .false.
      do ipair = 1, npair_renamexf
         do iq = 1, nspecfrm_renamexf(ipair)
            lsfrm = lspecfrma_renamexf(iq,ipair)
            lstoo = lspectooa_renamexf(iq,ipair)
            if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotend(lsfrm) = .true.
               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                  dotend(lstoo) = .true.
               end if
            end if

            lsfrm = lspecfrmc_renamexf(iq,ipair)
            lstoo = lspectooc_renamexf(iq,ipair)
            if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotendqqcw(lsfrm) = .true.
               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                  dotendqqcw(lstoo) = .true.
               end if
            end if
         end do ! iq = ...
      end do ! ipair = ...

      call init_col_tend( sfgaex2_diag, sfgaex2_fld, dotend, cnst_name, '_sfgaex2', &
                          'gas-aerosol-exchange renaming column tendency', history_aerosol )
      call init_col_tend( sfgaex2cw_diag, sfgaex2cw_fld, dotendqqcw, cw_cnst_name, '_sfgaex2', &
                          'gas-aerosol-exchange renaming column tendency', history_aerosol )


! set for used in aging calcs:
!    fac_m2v_so4, fac_m2v_nh4, fac_m2v_soa(:)
!    soa_equivso4_factor(:)
      soa_equivso4_factor = 0.0_r8
      if (modefrm_pcage > 0) then
         n = mode_accum
         l = so4_a_ndx(n) ; l2 = -1
         if (l <= 0) call endrun( 'modal_aero_gasaerexch_init error a001 finding accum. so4' )
         do l1 = 1, nspec_mode(n)
            if (mass_a_ndx(l1,n) == l) then
!               l2 = lspectype_amode(l1,n)
               l2 = l1
!               fac_m2v_so4 = specmw_amode(l2) / specdens_amode(l2)
               call aero_props%get(n, l1, spec_mw=spec_mw, density=spec_dens, hygro=spec_hygro)
               fac_m2v_so4 = spec_mw / spec_dens
!               tmp2 = spechygro(l2)
               tmp2 = spec_hygro

            end if
         end do
         if (l2 <= 0) call endrun( 'modal_aero_gasaerexch_init error a002 finding accum. so4' )

         l = nh4_a_ndx(n) ; l2 = -1
         if (l > 0) then
            do l1 = 1, nspec_mode(n)
               if (mass_a_ndx(l1,n) == l) then
!                  l2 = lspectype_amode(l1,n)
                  l2 = l1
!                  fac_m2v_nh4 = specmw_amode(l2) / specdens_amode(l2)
                  call aero_props%get(n, l1, spec_mw=spec_mw, density=spec_dens)
                  fac_m2v_nh4 = spec_mw / spec_dens

               end if
            end do
            if (l2 <= 0) call endrun( 'modal_aero_gasaerexch_init error a002 finding accum. nh4' )
         else
            fac_m2v_nh4 = fac_m2v_so4
         end if

         do jsoa = 1, nsoa
            l = soa_a_ndx(n,jsoa) ; l2 = -1
            if (l <= 0) then
               write( msg, '(a,i4)') 'modal_aero_gasaerexch_init error a001 finding accum. jsoa =', jsoa
               call endrun( msg )
            end if
            do l1 = 1, nspec_mode(n)
               if (mass_a_ndx(l1,n) == l) then
!                  l2 = lspectype_amode(l1,n)
                  l2 = l1
!                  fac_m2v_soa(jsoa) = specmw_amode(l2) / specdens_amode(l2)
                  call aero_props%get(n, l1, spec_mw=spec_mw, density=spec_dens, hygro=spec_hygro)
                  fac_m2v_soa(jsoa) = spec_mw / spec_dens
!                  soa_equivso4_factor(jsoa) = spechygro(l2)/tmp2
                  soa_equivso4_factor(jsoa) = spec_hygro/tmp2
               end if
            end do
            if (l2 <= 0) then
               write( msg, '(a,i4)') 'modal_aero_gasaerexch_init error a002 finding accum. jsoa =', jsoa
               call endrun( msg )
            end if
         end do

         fac_m2v_pcarbon(:) = 0.0_r8
         n = mode_pcarbon
         do l = 1, nspec_mode(n)
!            l2 = lspectype_amode(l,n)
!      fac_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!           [m3-AP/kmol-AP]    = [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
!            fac_m2v_pcarbon(l) = specmw_amode(l2) / specdens_amode(l2)
            call aero_props%get(n, l, spec_mw=spec_mw, density=spec_dens)
            fac_m2v_pcarbon(l) = spec_mw / spec_dens
         end do
      end if


      return

      end subroutine modal_aero_gasaerexch_init

!----------------------------------------------------------------------
! registers the column tendency fields of the constituents selected by mask
! and sets the field index of each constituent (0 if none)
!----------------------------------------------------------------------
      subroutine init_col_tend( diag, fld_ndx, mask, cnames, suffix, desc, history_aerosol )

      type(col_tend_diag_t), intent(inout) :: diag
      integer, allocatable, intent(out) :: fld_ndx(:)  ! field index of each constituent
      logical, intent(in) :: mask(:)                   ! constituents with column tendency fields
      character(len=*), intent(in) :: cnames(:)        ! constituent names
      character(len=*), intent(in) :: suffix           ! field name suffix
      character(len=*), intent(in) :: desc             ! long name description
      logical, intent(in) :: history_aerosol           ! add fields to default history output

      integer :: l, nfld
      character(len=len(cnames)), allocatable :: names(:)

      allocate(fld_ndx(size(mask)))
      allocate(names(count(mask)))

      nfld = 0
      do l = 1, size(mask)
         if (mask(l)) then
            nfld = nfld + 1
            names(nfld) = cnames(l)
            fld_ndx(l) = nfld
         else
            fld_ndx(l) = 0
         end if
      end do

      call diag%init( names, suffix, desc, history_aerosol )

      end subroutine init_col_tend

!----------------------------------------------------------------------

!----------------------------------------------------------------------
! sets the numbers of modes and species, the constituent indices of mode
! number and species, the primary carbon and accumulation modes, and the
! cloud-borne constituent names from the aerosol properties object
!----------------------------------------------------------------------
      subroutine set_mode_tables(aero_props)

      use aerosol_properties_mod, only: aerosol_properties, aero_name_len
      use constituents,   only: pcnst, cnst_get_ind
      use cam_abortutils, only: endrun
      use modal_aero_data, only: ntot_amode_mad => ntot_amode, nspec_amode_mad => nspec_amode, &
                                 nspec_max_mad => nspec_max, numptr_amode_mad => numptr_amode, &
                                 lmassptr_amode_mad => lmassptr_amode, alnsg_amode_mad => alnsg_amode, &
                                 modeptr_pcarbon_mad => modeptr_pcarbon, modeptr_accum_mad => modeptr_accum, &
                                 specmw_amode_mad => specmw_amode, specdens_amode_mad => specdens_amode, &
                                 spechygro_mad => spechygro, cnst_name_cw_mad => cnst_name_cw

      class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object

      integer :: m, l, idx
      real(r8) :: spec_mw, spec_dens, spec_hygro
      character(len=aero_name_len) :: name_a, name_c, mode_type
      character(len=*), parameter :: subname = 'modal_aero_gasaerexch::set_mode_tables'

      nmodes = aero_props%nbins()

      allocate(nspec_mode(nmodes))
      do m = 1, nmodes
         nspec_mode(m) = aero_props%nspecies(m)
      end do
      nspec_max = maxval(nspec_mode)

      allocate(num_a_ndx(nmodes))
      allocate(mass_a_ndx(nspec_max,nmodes))
      allocate(cw_cnst_name(pcnst))
      mass_a_ndx(:,:) = -1
      cw_cnst_name(:) = ' '
      mode_pcarbon = -1
      mode_accum = -1

      do m = 1, nmodes
         call aero_props%num_names(m, name_a, name_c)
         call cnst_get_ind(name_a, idx, abort=.false.)
         num_a_ndx(m) = idx
         if (idx > 0) cw_cnst_name(idx) = name_c
         do l = 1, nspec_mode(m)
            call aero_props%mmr_names(m, l, name_a, name_c)
            call cnst_get_ind(name_a, idx, abort=.false.)
            mass_a_ndx(l,m) = idx
            if (idx > 0) cw_cnst_name(idx) = name_c
         end do
         ! the last mode of each type is used
         mode_type = aero_props%bin_name(m)
         if (trim(mode_type) == 'accum') mode_accum = m
         if (trim(mode_type) == 'primary_carbon') mode_pcarbon = m
      end do

      ! check consistency with modal_aero_data
      if (nmodes /= ntot_amode_mad .or. nspec_max /= nspec_max_mad) then
         call endrun(subname//': inconsistent number of modes or maximum number of species')
      end if
      if (.not. same_mode(mode_accum, modeptr_accum_mad) .or. &
          .not. same_mode(mode_pcarbon, modeptr_pcarbon_mad)) then
         call endrun(subname//': inconsistent accumulation or primary carbon mode')
      end if
      do m = 1, nmodes
         if (nspec_mode(m) /= nspec_amode_mad(m) .or. num_a_ndx(m) /= numptr_amode_mad(m)) then
            call endrun(subname//': inconsistent number of species or number index')
         end if
         if (aero_props%alogsig(m) /= alnsg_amode_mad(m)) then
            call endrun(subname//': inconsistent log of geometric standard deviation')
         end if
         do l = 1, nspec_mode(m)
            if (mass_a_ndx(l,m) /= lmassptr_amode_mad(l,m)) then
               call endrun(subname//': inconsistent species index')
            end if
            call aero_props%get(m, l, spec_mw=spec_mw, density=spec_dens, hygro=spec_hygro)
            if (spec_mw /= specmw_amode_mad(l,m) .or. spec_dens /= specdens_amode_mad(l,m) .or. &
                spec_hygro /= spechygro_mad(l,m)) then
               call endrun(subname//': inconsistent species properties')
            end if
         end do
      end do
      do idx = 1, pcnst
         if (trim(cw_cnst_name(idx)) /= trim(cnst_name_cw_mad(idx))) then
            call endrun(subname//': inconsistent cloud-borne name of '//trim(cnst_name_cw_mad(idx)))
         end if
      end do

      contains

         ! true if two mode indices are the same mode, or both indicate no mode (<= 0)
         pure logical function same_mode(m1, m2)
            integer, intent(in) :: m1, m2
            same_mode = (m1 == m2) .or. (m1 <= 0 .and. m2 <= 0)
         end function same_mode

      end subroutine set_mode_tables

!----------------------------------------------------------------------
! sets constituent indices of the interstitial aerosol species and soa gas
! species involved in gas-aerosol exchange, and the soa/poa molecular
! weights, from the aerosol properties object
!----------------------------------------------------------------------
      subroutine set_species_ndx(aero_props)

      use aerosol_properties_mod, only: aerosol_properties, aero_name_len
      use aero_gasaerexch, only: aero_gasaerexch_spec_ndx, aero_gasaerexch_soag_ndx
      use constituents,   only: pcnst, cnst_get_ind
      use cam_abortutils, only: endrun
      use modal_aero_data, only: lptr_so4_a_amode, lptr_nh4_a_amode, &
                                 lptr2_soa_a_amode, lptr2_pom_a_amode, lptr2_soa_g_amode, &
                                 nsoa_mad => nsoa, npoa_mad => npoa

      class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object

      integer :: m, j
      integer :: nspec_typ(nmodes)
      integer :: nspec_soa(nmodes), nspec_poa(nmodes)
      integer, allocatable :: spec_ndx(:,:)
      integer, allocatable :: spec_ndx_soa(:,:), spec_ndx_poa(:,:)
      integer, allocatable :: soag_ndx(:)
      real(r8) :: spec_mw
      logical :: mw_soa_set
      character(len=128) :: msg
      character(len=*), parameter :: subname = 'modal_aero_gasaerexch::set_species_ndx'

      ! numbers of soa and poa species
      call aero_gasaerexch_spec_ndx(aero_props, 's-organic', nspec_soa, spec_ndx_soa)
      call aero_gasaerexch_spec_ndx(aero_props, 'p-organic', nspec_poa, spec_ndx_poa)
      nsoa = maxval(nspec_soa)
      npoa = maxval(nspec_poa)

      allocate(so4_a_ndx(nmodes))
      allocate(nh4_a_ndx(nmodes))
      allocate(soa_a_ndx(nmodes,nsoa))
      allocate(pom_a_ndx(nmodes,npoa))
      allocate(soag_cnst_ndx(nsoa))
      allocate(mw_soa_host(nsoa))
      allocate(mw_poa_host(npoa))

      so4_a_ndx(:) = -1
      nh4_a_ndx(:) = -1
      soa_a_ndx(:,:) = -1
      pom_a_ndx(:,:) = -1
      soag_cnst_ndx(:) = -1
      mw_soa_host(:) = -huge(1.0_r8)
      mw_poa_host(:) = -huge(1.0_r8)
      mw_soa_set = .false.

      ! sulfate and ammonium -- the last species of the type in a mode is used
      call aero_gasaerexch_spec_ndx(aero_props, 'sulfate', nspec_typ, spec_ndx)
      do m = 1, nmodes
         if (nspec_typ(m) > 0) so4_a_ndx(m) = cnst_ndx(m, spec_ndx(m,nspec_typ(m)))
      end do

      call aero_gasaerexch_spec_ndx(aero_props, 'ammonium', nspec_typ, spec_ndx)
      do m = 1, nmodes
         if (nspec_typ(m) > 0) nh4_a_ndx(m) = cnst_ndx(m, spec_ndx(m,nspec_typ(m)))
      end do

      ! the last s-organic (p-organic) species found sets the soa (poa) molecular weight
      do m = 1, nmodes
         do j = 1, nspec_soa(m)
            soa_a_ndx(m,j) = cnst_ndx(m, spec_ndx_soa(m,j))
            call aero_props%get(m, spec_ndx_soa(m,j), spec_mw=spec_mw)
            mw_soa_host(:) = spec_mw
            mw_soa_set = .true.
         end do
      end do

      do m = 1, nmodes
         do j = 1, nspec_poa(m)
            pom_a_ndx(m,j) = cnst_ndx(m, spec_ndx_poa(m,j))
            call aero_props%get(m, spec_ndx_poa(m,j), spec_mw=spec_mw)
            mw_poa_host(:) = spec_mw
         end do
      end do

      if (nsoa > 0 .and. .not. mw_soa_set) then
         call endrun(subname//': no s-organic species found')
      end if

      ! soa gas species
      call aero_gasaerexch_soag_ndx(soag_ndx)
      do j = 1, min(nsoa, size(soag_ndx))
         soag_cnst_ndx(j) = soag_ndx(j)
      end do

      ! check consistency with the modal_aero_data species counts and pointers
      if (nsoa /= nsoa_mad .or. npoa /= npoa_mad) then
         write(msg,'(a,4i6)') subname//': inconsistent numbers of soa/poa species ', nsoa, nsoa_mad, npoa, npoa_mad
         call endrun(trim(msg))
      end if
      do m = 1, nmodes
         call check_ndx(so4_a_ndx(m), lptr_so4_a_amode(m), 'sulfate')
         call check_ndx(nh4_a_ndx(m), lptr_nh4_a_amode(m), 'ammonium')
         do j = 1, nsoa
            call check_ndx(soa_a_ndx(m,j), lptr2_soa_a_amode(m,j), 's-organic')
         end do
         do j = 1, npoa
            call check_ndx(pom_a_ndx(m,j), lptr2_pom_a_amode(m,j), 'p-organic')
         end do
      end do
      do j = 1, nsoa
         if (soag_cnst_ndx(j) /= lptr2_soa_g_amode(j)) then
            write(msg,'(a,3i8)') subname//': inconsistent soa gas index ', j, soag_cnst_ndx(j), lptr2_soa_g_amode(j)
            call endrun(trim(msg))
         end if
      end do

      contains

         ! constituent index of the interstitial species l of mode mm (<= 0 if not a constituent)
         ! (not pure: queries the aerosol properties object and the constituents module)
         integer function cnst_ndx(mm, l)
            integer, intent(in) :: mm  ! mode index
            integer, intent(in) :: l   ! species index within the mode

            character(len=aero_name_len) :: name_a, name_c

            call aero_props%mmr_names(mm, l, name_a, name_c)
            call cnst_get_ind(name_a, cnst_ndx, abort=.false.)

         end function cnst_ndx

         subroutine check_ndx(ndx, lptr, spectype)
            integer, intent(in) :: ndx              ! index from aerosol properties
            integer, intent(in) :: lptr             ! index from modal_aero_data
            character(len=*), intent(in) :: spectype ! species type

            if ((lptr > 0) .and. (lptr <= pcnst)) then
               if (ndx == lptr) return
            else
               if (ndx <= 0) return
            end if
            write(msg,'(a,2i8)') subname//': inconsistent '//trim(spectype)//' index for mode ', m, ndx
            call endrun(trim(msg))

         end subroutine check_ndx

      end subroutine set_species_ndx

!----------------------------------------------------------------------


!----------------------------------------------------------------------

end module modal_aero_gasaerexch
