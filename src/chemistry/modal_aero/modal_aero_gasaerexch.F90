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
  use modal_aero_data, only:  nspec_max, nsoa, npoa, soa_multi_species
  use ref_pres,        only:  top_lev => clim_modal_aero_top_lev
  use ppgrid,          only:  pcols, pver
  use modal_aero_data, only:  ntot_amode, numptr_amode, sigmag_amode
  use gas_aer_uptkrate_mod, only: gas_aer_uptkrates, nonvolatile_uptake_cell, h2so4_reversible_uptake_cell
  use modal_aero_data, only: lptr2_soa_g_amode, lptr2_soa_a_amode, lptr2_pom_a_amode
  use aero_soaexch_mod, only: soa_gasaerexch_column, method_soa
  use soa_gaex_diag_mod, only: soa_gaex_diag_t

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

  logical, allocatable :: do_soaexch(:)

! gases that condense on the modal aerosol, set in gas_table_init
  type :: condensing_gas_t
     integer  :: lgas = -1                         ! constituent index of the gas
     real(r8) :: diffus_ratio = -huge(1.0_r8)      ! gas diffusivity relative to h2so4
     integer, allocatable :: laer(:)               ! constituent index of the condensed species in each mode
                                                   !    (<= 0 if the mode does not have it)
     integer, allocatable :: ido(:)                ! uptake by each mode: 0 = none, 1 = uptake,
                                                   !    2 = uptake, then transfer to modetoo_pcage
                                                   !        by primary carbon aging
  end type condensing_gas_t

  integer :: ngas = 0                              ! number of condensing gases
  type(condensing_gas_t), allocatable :: gas(:)
  integer :: igas_h2so4 = -1                       ! index of h2so4 in gas(:)
  integer :: igas_nh3 = -1                         ! index of nh3 in gas(:) (-1 if nh3 is absent)
  integer, allocatable :: igas_soa(:)              ! index of each soa gas in gas(:) (-1 if absent)
  logical :: do_soag_any = .false.                 ! true if any soa gas is present
  logical, allocatable :: dotend_gaex(:)           ! true for constituents changed by gas uptake
                                                   !    and primary carbon aging (pcnst)
  real(r8) :: mw_soa_host = -huge(1.0_r8)          ! molec wght of soa used in host code (kg/kmol)
  type(soa_gaex_diag_t) :: soa_diag                ! soa condensation/evaporation diagnostics
  real(r8) :: fac_volsfc_pcarbon = -huge(1.0_r8)   ! (mode-surface-area/mode-volume) factor for the
                                                   !    primary carbon mode

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

      subroutine modal_aero_gasaerexch_init

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

use modal_aero_data
use modal_aero_rename

use cam_abortutils, only: endrun
use cam_history,    only: addfld, add_default, fieldname_len, horiz_only
use constituents,   only: pcnst, cnst_get_ind, cnst_name
use spmd_utils,     only: masterproc
use phys_control,   only: phys_getopts
use radiative_aerosol,  only: rad_aer_get_info

implicit none

!-----------------------------------------------------------------------
! arguments

!-----------------------------------------------------------------------
! local
   integer  :: ipair, iq, iqfrm, iqtoo
   integer  :: jac,jsoa
   integer  :: l, l1, l2, lsfrm, lstoo, lunout
   integer  :: mfrm, mtoo
   integer  :: n
   integer  :: nchfrm, nchfrmskip, nchtoo, nchtooskip, nspec

   logical  :: dotend(pcnst), dotendqqcw(pcnst)

   real(r8) :: tmp2

   character(len=fieldname_len)   :: tmpnamea
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(128)                 :: msg
   character(8)                   :: unit

   logical                        :: history_aerosol      ! Output the MAM aerosol tendencies
   character(len=32) :: spec_type
   character(len=len(cnst_name)) :: soag_names(nsoa)    ! soa gas species names
   !-----------------------------------------------------------------------

        call phys_getopts( history_aerosol_out        = history_aerosol   )

        maxspec_pcage = nspec_max
        allocate(lspecfrm_pcage(maxspec_pcage))
        allocate(lspectoo_pcage(maxspec_pcage))
        allocate(soa_equivso4_factor(nsoa))
        allocate(fac_m2v_soa(nsoa))
        allocate(fac_m2v_pcarbon(nspec_max))
        lunout = 6

        allocate(do_soaexch(ntot_amode))
        do_soaexch(:) = .false.

        do n = 1, ntot_amode
           !call rad_aer_get_info(0, n, mode_type=mode_type)

           do l = 1, nspec_amode(n)
              call rad_aer_get_info(0, n, l, spec_type=spec_type )
              if ( trim(spec_type) == 's-organic' .or. trim(spec_type) == 'p-organic' ) then
                 ! allow for exchange with primary organic mode
                 do_soaexch(n) = .true.
              end if
           end do
        end do
!print*,'FVDBG: do_soaexch: ',do_soaexch
!
!   define "from mode" and "to mode" for primary carbon aging
!
!   skip (turn off) aging if either is absent,
!      or if accum mode so4 is absent
!
	modefrm_pcage = -999888777
	modetoo_pcage = -999888777
	if ((modeptr_pcarbon <= 0) .or. (modeptr_accum <= 0)) goto 15000
	l = lptr_so4_a_amode(modeptr_accum)
	if ((l < 1) .or. (l > pcnst)) goto 15000

	modefrm_pcage = modeptr_pcarbon
	modetoo_pcage = modeptr_accum

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

aa_iqfrm: do iqfrm = -1, nspec_amode(mfrm)

	    if (iqfrm == -1) then
		lsfrm = numptr_amode(mfrm)
		lstoo = numptr_amode(mtoo)
	    else if (iqfrm == 0) then
!   bypass transfer of aerosol water due to primary-carbon aging
		cycle aa_iqfrm
!               lsfrm = lwaterptr_amode(mfrm)
!               lstoo = lwaterptr_amode(mtoo)
	    else
		lsfrm = lmassptr_amode(iqfrm,mfrm)
		lstoo = 0
	    end if
	    if ((lsfrm < 1) .or. (lsfrm > pcnst)) cycle aa_iqfrm

	    if (lsfrm>0 .and. iqfrm>0 ) then
		nchfrm = len( trim( cnst_name(lsfrm) ) ) - nchfrmskip

! find "too" species having same lspectype_amode as the "frm" species
! AND same cnst_name (except for last 1/2/3 characters which are the mode index)
		do iqtoo = 1, nspec_amode(mtoo)
!		    if ( lspectype_amode(iqtoo,mtoo) .eq.   &
!		         lspectype_amode(iqfrm,mfrm) ) then
			lstoo = lmassptr_amode(iqtoo,mtoo)
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

! set the table of condensing gases and the tendency flags
      call gas_table_init()

!---------define history fields for new cond/evap diagnostics----------------------------------------
      do jsoa = 1, nsoa
         soag_names(jsoa) = cnst_name(lptr2_soa_g_amode(jsoa))
      end do
      call soa_diag%init( soag_names, history_aerosol, add_default_species=.false. )
!------------------------------------------------------------------------------

!  define history fields for basic gas-aer exchange
!  and primary carbon aging from that
      do l = 1, pcnst
         if ( .not. dotend_gaex(l) ) cycle

         tmpnamea = cnst_name(l)
         fieldname = trim(tmpnamea) // '_sfgaex1'
         long_name = trim(tmpnamea) // ' gas-aerosol-exchange primary column tendency'
         unit = 'kg/m2/s'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if ( history_aerosol ) then
            call add_default( fieldname, 1, ' ' )
         endif
         if ( masterproc ) write(*,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit

      end do   ! l = ...
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

      do l = 1, pcnst
      do jac = 1, 2
         if (jac == 1) then
            if ( .not. dotend(l) ) cycle
            tmpnamea = cnst_name(l)
         else
            if ( .not. dotendqqcw(l) ) cycle
            tmpnamea = cnst_name_cw(l)
         end if

         fieldname = trim(tmpnamea) // '_sfgaex2'
         long_name = trim(tmpnamea) // ' gas-aerosol-exchange renaming column tendency'
         unit = 'kg/m2/s'
         if ((tmpnamea(1:3) == 'num') .or. &
             (tmpnamea(1:3) == 'NUM')) unit = '#/m2/s'
         call addfld( fieldname, horiz_only, 'A', unit, long_name )
         if ( history_aerosol ) then
            call add_default( fieldname, 1, ' ' )
         endif
         if ( masterproc ) write(*,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit
      end do   ! jac = ...
      end do   ! l = ...


! set for used in aging calcs:
!    fac_m2v_so4, fac_m2v_nh4, fac_m2v_soa(:)
!    soa_equivso4_factor(:)
      soa_equivso4_factor = 0.0_r8
      if (modefrm_pcage > 0) then
         n = modeptr_accum
         l = lptr_so4_a_amode(n) ; l2 = -1
         if (l <= 0) call endrun( 'modal_aero_gasaerexch_init error a001 finding accum. so4' )
         do l1 = 1, nspec_amode(n)
            if (lmassptr_amode(l1,n) == l) then
!               l2 = lspectype_amode(l1,n)
               l2 = l1
!               fac_m2v_so4 = specmw_amode(l2) / specdens_amode(l2)
               fac_m2v_so4 = specmw_amode(l1,n) / specdens_amode(l1,n)
!               tmp2 = spechygro(l2)
               tmp2 = spechygro(l1,n)

            end if
         end do
         if (l2 <= 0) call endrun( 'modal_aero_gasaerexch_init error a002 finding accum. so4' )

         l = lptr_nh4_a_amode(n) ; l2 = -1
         if (l > 0) then
            do l1 = 1, nspec_amode(n)
               if (lmassptr_amode(l1,n) == l) then
!                  l2 = lspectype_amode(l1,n)
                  l2 = l1
!                  fac_m2v_nh4 = specmw_amode(l2) / specdens_amode(l2)
                  fac_m2v_nh4 = specmw_amode(l1,n) / specdens_amode(l1,n)

               end if
            end do
            if (l2 <= 0) call endrun( 'modal_aero_gasaerexch_init error a002 finding accum. nh4' )
         else
            fac_m2v_nh4 = fac_m2v_so4
         end if

         do jsoa = 1, nsoa
            l = lptr2_soa_a_amode(n,jsoa) ; l2 = -1
            if (l <= 0) then
               write( msg, '(a,i4)') 'modal_aero_gasaerexch_init error a001 finding accum. jsoa =', jsoa
               call endrun( msg )
            end if
            do l1 = 1, nspec_amode(n)
               if (lmassptr_amode(l1,n) == l) then
!                  l2 = lspectype_amode(l1,n)
                  l2 = l1
!                  fac_m2v_soa(jsoa) = specmw_amode(l2) / specdens_amode(l2)
                  fac_m2v_soa(jsoa) = specmw_amode(l1,n) / specdens_amode(l1,n)
!                  soa_equivso4_factor(jsoa) = spechygro(l2)/tmp2
                  soa_equivso4_factor(jsoa) = spechygro(l1,n)/tmp2
               end if
            end do
            if (l2 <= 0) then
               write( msg, '(a,i4)') 'modal_aero_gasaerexch_init error a002 finding accum. jsoa =', jsoa
               call endrun( msg )
            end if
         end do

         fac_m2v_pcarbon(:) = 0.0_r8
         n = modeptr_pcarbon
         do l = 1, nspec_amode(n)
!            l2 = lspectype_amode(l,n)
!      fac_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!           [m3-AP/kmol-AP]    = [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
!            fac_m2v_pcarbon(l) = specmw_amode(l2) / specdens_amode(l2)
            fac_m2v_pcarbon(l) = specmw_amode(l,n) / specdens_amode(l,n)
         end do
      end if


      return

      end subroutine modal_aero_gasaerexch_init

!----------------------------------------------------------------------
subroutine gas_table_init

!-----------------------------------------------------------------------
!
! Purpose:
!    set the table of gases that condense on the modal aerosol
!    (h2so4, nh3 and soa gases), the modes that take them up,
!    and the constituents changed by gas uptake and primary carbon aging
!
!    called by modal_aero_gasaerexch_init after the primary carbon
!    aging pointers are set
!
!-----------------------------------------------------------------------

   use modal_aero_data, only: lptr_so4_a_amode, lptr_nh4_a_amode, lptr2_soa_a_amode, &
                              modeptr_pcarbon, alnsg_amode, nspec_amode, specmw_amode
   use radiative_aerosol, only: rad_aer_get_info
   use aero_soaexch_mod, only: soa_h2so4_diffus_ratio
   use constituents,    only: pcnst, cnst_get_ind
   use cam_abortutils,  only: endrun

   integer :: igas, iq, jsoa, l, lgas, lsfrm, lstoo, n
   character(len=32) :: spec_type

   if ((method_soa /= 0) .and. (method_soa /= 1) .and. (method_soa /= 2)) then
      write(*,'(/a,1x,i10)') '*** modal_aero_gasaerexch_init - bad method_soa =', method_soa
      call endrun( 'modal_aero_gasaerexch_init error' )
   end if

   allocate( gas(2+nsoa) )
   allocate( igas_soa(nsoa) )
   allocate( dotend_gaex(pcnst) )
   igas_soa(:) = -1
   dotend_gaex(:) = .false.
   ngas = 0

! h2so4 (required) condenses as so4
   call cnst_get_ind( 'H2SO4', lgas, .false. )
   if ((lgas <= 0) .or. (lgas > pcnst)) then
      write( *, '(/a/a,2i7)' )   &
         '*** modal_aero_gasaerexch_init -- cannot find H2SO4 species',   &
         '    l_so4g=', lgas
      call endrun( 'modal_aero_gasaerexch_init error' )
   end if
   ngas = ngas + 1
   igas_h2so4 = ngas
   gas(ngas)%lgas = lgas
   gas(ngas)%diffus_ratio = 1.0_r8
   allocate( gas(ngas)%laer(ntot_amode), gas(ngas)%ido(ntot_amode) )
   gas(ngas)%laer(:) = lptr_so4_a_amode(1:ntot_amode)
   gas(ngas)%ido(:) = 0
   do n = 1, ntot_amode
      l = gas(ngas)%laer(n)
      if ((l > 0) .and. (l <= pcnst)) gas(ngas)%ido(n) = 1
   end do

! nh3 (optional) condenses as nh4, only on modes that take up h2so4
   call cnst_get_ind( 'NH3', lgas, .false. )
   if ((lgas > 0) .and. (lgas <= pcnst)) then
      ngas = ngas + 1
      igas_nh3 = ngas
      gas(ngas)%lgas = lgas
      ! 2.08 factor is for gas diffusivity (nh3/h2so4)
      ! differences in fuch-sutugin and accom coef ignored
      gas(ngas)%diffus_ratio = 2.08_r8
      allocate( gas(ngas)%laer(ntot_amode), gas(ngas)%ido(ntot_amode) )
      gas(ngas)%laer(:) = lptr_nh4_a_amode(1:ntot_amode)
      gas(ngas)%ido(:) = 0
      do n = 1, ntot_amode
         l = gas(ngas)%laer(n)
         if ((gas(igas_h2so4)%ido(n) == 1) .and. (l > 0) .and. (l <= pcnst)) gas(ngas)%ido(n) = 1
      end do
   end if

! soa gases condense as soa (method_soa = 1 or 2)
   if (method_soa > 0) then
      do jsoa = 1, nsoa
         lgas = lptr2_soa_g_amode(jsoa)
         if ((lgas <= 0) .or. (lgas > pcnst)) cycle
         ngas = ngas + 1
         igas_soa(jsoa) = ngas
         gas(ngas)%lgas = lgas
         gas(ngas)%diffus_ratio = soa_h2so4_diffus_ratio
         allocate( gas(ngas)%laer(ntot_amode), gas(ngas)%ido(ntot_amode) )
         gas(ngas)%laer(:) = lptr2_soa_a_amode(1:ntot_amode,jsoa)
         gas(ngas)%ido(:) = 0
         do n = 1, ntot_amode
            l = gas(ngas)%laer(n)
            if ((l > 0) .and. (l <= pcnst)) gas(ngas)%ido(n) = 1
         end do
      end do
   end if
   do_soag_any = any( igas_soa(:) > 0 )

! molecular weight of soa used in the host code, from the species definitions
   do n = 1, ntot_amode
      do l = 1, nspec_amode(n)
         call rad_aer_get_info(0, n, l, spec_type=spec_type )
         select case( spec_type )
         case('s-organic')
            mw_soa_host = specmw_amode(l,n)
         end select
      end do
   end do

! constituents changed by gas uptake
! (set before the primary carbon aging changes ido below)
   do igas = 1, ngas
      dotend_gaex(gas(igas)%lgas) = .true.
      do n = 1, ntot_amode
         if (gas(igas)%ido(n) == 1) dotend_gaex(gas(igas)%laer(n)) = .true.
      end do
   end do

! primary carbon aging: the primary carbon mode takes up h2so4, and also takes up
! nh3 and soa gases that the accumulation mode takes up, with the condensed
! species transferred to the accumulation mode
   if (modefrm_pcage > 0) then
      gas(igas_h2so4)%ido(modefrm_pcage) = 2
      do igas = 1, ngas
         if (igas == igas_h2so4) cycle
         if (gas(igas)%ido(modetoo_pcage) == 1) gas(igas)%ido(modefrm_pcage) = 2
      end do

      do iq = 1, nspecfrm_pcage
         lsfrm = lspecfrm_pcage(iq)
         lstoo = lspectoo_pcage(iq)
         if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
            dotend_gaex(lsfrm) = .true.
            if ((lstoo > 0) .and. (lstoo <= pcnst)) then
               dotend_gaex(lstoo) = .true.
            end if
         end if
      end do

      fac_volsfc_pcarbon = exp( 2.5_r8*(alnsg_amode(modeptr_pcarbon)**2) )
   end if

end subroutine gas_table_init


!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_gasaerexch_sub --- ...
!
! !INTERFACE:
subroutine modal_aero_gasaerexch_sub(                            &
                        lchnk,    ncol,     nstep,               &
                        loffset,  deltat,                        &
                        t,        pmid,     pdel,                &
                        troplev,                                 &
                        q,                  qqcw,                &
                        dqdt_other,         dqqcwdt_other,       &
                        dgncur_a,           dgncur_awet,         &
                        sulfeq         )

! !USES:
use modal_aero_data,   only:  lmassptr_amode,cnst_name_cw
use modal_aero_data,   only:  lptr_so4_a_amode,lptr_nh4_a_amode
use modal_aero_data,   only:  modeptr_pcarbon,nspec_amode
use modal_aero_rename, only:  modal_aero_rename_sub

use cam_history,       only:  outfld, fieldname_len
use chem_mods,         only:  adv_mass
use constituents,      only:  pcnst, cnst_name
use physconst,         only:  gravit, mwdry, rair
use cam_abortutils,    only:  endrun

implicit none

! !PARAMETERS:
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
   real(r8), intent(in)    :: dgncur_a(pcols,pver,ntot_amode)
   real(r8), intent(in)    :: dgncur_awet(pcols,pver,ntot_amode)
   real(r8), pointer       :: sulfeq(:,:,:)

                                 ! dry & wet geo. mean dia. (m) of number distrib.

! !DESCRIPTION:
! computes TMR (tracer mixing ratio) tendencies for gas condensation
!    onto aerosol particles
!
! this version does condensation of H2SO4 and NH3, both treated as
! completely non-volatile (gas --> aerosol, but no aerosol --> gas)
!    gas H2SO4 goes to aerosol SO4
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

   integer :: i, iq, itmpa
   integer :: ido_so4a(ntot_amode), ido_nh4a(ntot_amode)
   integer ::  ido_soaa(ntot_amode,nsoa)
   integer :: j, jac, jsrf, jsoa
   integer :: k
   integer :: l, lb, lsfrm, lstoo
   integer :: l_so4g, l_nh4g
   integer :: l_soag(nsoa)
   integer :: n, nadv

   logical :: is_dorename_atik, dorename_atik(ncol,pver)

   character(len=fieldname_len+3) :: fieldname

   real (r8) :: deltatxx
   real (r8) :: dqdt_nh4(ntot_amode), dqdt_so4(ntot_amode)
   real (r8) :: dqdt_soa(ntot_amode,nsoa)
   real (r8) :: pdel_fac
   real (r8) :: qmax_nh4, qnew_nh4, qnew_so4
   real (r8) :: qold_nh4(ntot_amode), qold_so4(ntot_amode)
   real (r8), allocatable :: qold_soag_c(:,:,:)    ! (nsoa,pcols,pver) soa gas mixrat
   real (r8), allocatable :: qold_soa_c(:,:,:,:)   ! (ntot_amode,nsoa,pcols,pver) soa aerosol mixrat
   real (r8), allocatable :: qold_poa_c(:,:,:,:)   ! (ntot_amode,npoa,pcols,pver) poa aerosol mixrat
   real (r8), allocatable :: dqdt_soa_c(:,:,:,:)   ! (ntot_amode,nsoa,pcols,pver) soa aerosol tendency
   real (r8), allocatable :: sum_dqdt_soa_c(:,:,:) ! (nsoa,pcols,pver) soa tendency summed over modes
   real (r8) :: diag_fac(ntot_amode,nsoa)  ! converts soa tendency to kg/kg/s for the diagnostics
   logical   :: diag_on(ntot_amode,nsoa)   ! true if soa tendency is included in the diagnostics
   real (r8) :: sum_dqdt_so4
   real (r8) :: sum_dqdt_soa(nsoa)
   real (r8) :: sum_dqdt_nh4, sum_dqdt_nh4_b
   real (r8) :: sum_uprt_nh4, sum_uprt_so4
   real (r8) :: tmp1, tmp2, tmpa
   real (r8) :: uptkrate(ntot_amode,pcols,pver)
                ! gas-to-aerosol mass transfer rates (1/s)
   real (r8) :: vol_core, vol_shell
   real (r8) :: xferfrac_pcage, xferfrac_max
   real (r8) :: xferrate

   logical  :: do_nh4g         ! true if nh3 gas is a species
   logical  :: do_soag(nsoa)       ! true if soa gas is a species
   logical  :: so4_on(ntot_amode)      ! true if h2so4 gas condenses on mode
   logical  :: nh4_on(ntot_amode)      ! true if nh3 gas condenses on mode
   logical  :: soa_on(ntot_amode,nsoa) ! true if soa gas condenses on mode

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
   real(r8), allocatable :: qcon(:,:,:), qevap(:,:,:)  ! soa condensation/evaporation diagnostics
                                                       !    (pcols,pver,soa_diag%nfld())

   real(r8) :: qqcwsrflx(pcols,pcnstxx,nsrflx)

   real(r8) :: aircon(ncol,pver)
   real(r8) :: num_a(pcols,pver,ntot_amode)

!----------------------------------------------------------------------

! set gas species indices and tendency flags from the table of condensing gases
   l_so4g = gas(igas_h2so4)%lgas - loffset
   ido_so4a(:) = gas(igas_h2so4)%ido(:)

   do_nh4g = (igas_nh3 > 0)
   if ( do_nh4g ) then
      l_nh4g = gas(igas_nh3)%lgas - loffset
      ido_nh4a(:) = gas(igas_nh3)%ido(:)
   else
      l_nh4g = -1
      ido_nh4a(:) = 0
   end if

   do jsoa = 1, nsoa
      l_soag(jsoa) = lptr2_soa_g_amode(jsoa) - loffset
      do_soag(jsoa) = (igas_soa(jsoa) > 0)
      if ( do_soag(jsoa) ) then
         ido_soaa(:,jsoa) = gas(igas_soa(jsoa))%ido(:)
      else
         ido_soaa(:,jsoa) = 0
      end if
   end do ! jsoa

   ! local (chemistry) index l is constituent index l+loffset; chemistry species
   ! beyond pcnst-loffset (e.g., short-lived species) are not constituents
   nadv = min( pcnstxx, pcnst-loffset )
   dotend(:) = .false.
   dotend(1:nadv) = dotend_gaex(1+loffset:nadv+loffset)
   dotendqqcw(:) = .false.

   if (modefrm_pcage > 0) then
      xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps
   end if


! zero out tendencies and other
   dqdt(:,:,:) = 0.0_r8
   dqqcwdt(:,:,:) = 0.0_r8
   qsrflx(:,:,:) = 0.0_r8
   qqcwsrflx(:,:,:) = 0.0_r8

!-------Initialize evap/cond diagnostics (ncols x pver)-----------
   allocate( qcon(pcols,pver,soa_diag%nfld()), qevap(pcols,pver,soa_diag%nfld()) )
   qcon(:,:,:) = 0.0_r8
   qevap(:,:,:) = 0.0_r8
!---------------------------------------------------

!   number conc. (#/m3) -- note q(i,k,numptr) is (#/kmol-air)
!   so need aircon in (kmol-air/m3)
   aircon(:ncol,:) = pmid(:ncol,:)/(rair*t(:ncol,:))/mwdry              ! (kmol-air/m3)
   do n = 1,ntot_amode
      !"bounded" number conc. (#/m3)
      num_a(:ncol,:,n) = q(:ncol,:,numptr_amode(n)-loffset)*aircon(:ncol,:)
   end do

! compute gas-to-aerosol mass transfer rates
   call gas_aer_uptkrates( ncol, ntot_amode, num_a, dgncur_awet, t, pmid, uptkrate, &
                           sigmag=sigmag_amode )


! use this for tendency calcs to avoid generating very small negative values
   deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)


! gases condense on modes with ido > 0
! (ido = 2 for the primary carbon mode, which transfers the condensed species
!  to the accumulation mode by primary carbon aging)
   so4_on(:) = ido_so4a(:) > 0
   nh4_on(:) = ido_nh4a(:) > 0
   soa_on(:,:) = ido_soaa(:,:) > 0

!   compute TMR tendencies for soag and soa interstial aerosol, for all columns
   if ( do_soag_any ) then
      allocate( qold_soag_c(nsoa,pcols,pver), qold_soa_c(ntot_amode,nsoa,pcols,pver), &
                qold_poa_c(ntot_amode,npoa,pcols,pver), dqdt_soa_c(ntot_amode,nsoa,pcols,pver), &
                sum_dqdt_soa_c(nsoa,pcols,pver) )

      do n = 1, ntot_amode
         do jsoa = 1, nsoa
            diag_on(n,jsoa) = do_soag(jsoa) .and. (ido_soaa(n,jsoa) == 1)
            if ( diag_on(n,jsoa) ) then
               l = lptr2_soa_a_amode(n,jsoa)-loffset
               diag_fac(n,jsoa) = adv_mass(l)/mwdry
            else
               diag_fac(n,jsoa) = 0.0_r8
            end if
         end do
      end do

      do k=top_lev,pver
         do i=1,ncol
            do n = 1, ntot_amode
               do j = 1, npoa
                  l = lptr2_pom_a_amode(n,j)-loffset
                  if (l > 0) then
                     qold_poa_c(n,j,i,k) = q(i,k,l)
                  else
                     qold_poa_c(n,j,i,k) = 0.0_r8
                  end if
               end do

               itmpa = 0
               do jsoa = 1, nsoa
                  if (ido_soaa(n,jsoa) == 1) then
                     l = lptr2_soa_a_amode(n,jsoa)-loffset
                     qold_soa_c(n,jsoa,i,k) = q(i,k,l)
                     itmpa = itmpa + 1
                  else
                     qold_soa_c(n,jsoa,i,k) = 0.0_r8
                  end if
               end do ! jsoa
               ! in previous code versions with nsoa=1,
               !    qold_poa was non-zero (i.e., loaded from q) only when ido_soaa(n)=1
               ! thus qold_poa=0 for the primary carbon mode which has ido_soaa=2
               ! this is probably not how it should be
               if (itmpa == 0) qold_poa_c(n,:,i,k) = 0.0_r8
            end do ! n

            do jsoa = 1, nsoa
               if ( do_soag(jsoa) ) then
                  qold_soag_c(jsoa,i,k) = q(i,k,l_soag(jsoa))
               else
                  qold_soag_c(jsoa,i,k) = 0.0_r8
               end if
            end do
         end do
      end do

!   sum_dqdt_soa = soa_a tendency from soa   gas uptake (mol/mol/s)
      call soa_gasaerexch_column( ncol, ntot_amode, nsoa, npoa, deltat, t, pmid, &
           soa_on, do_soaexch, do_soag, uptkrate, mw_soa_host, &
           qold_soag_c, qold_soa_c, qold_poa_c, dqdt_soa_c, sum_dqdt_soa_c, &
           soa_diag, diag_on, diag_fac, qcon, qevap )
   end if

   jsrf = jsrflx_gaexch
   do k=top_lev,pver
      do i=1,ncol

         do n = 1, ntot_amode
            if (ido_so4a(n) == 1) then
               qold_so4(n) = q(i,k,lptr_so4_a_amode(n)-loffset)
            else
               qold_so4(n) = 0.0_r8
            end if

            if (ido_nh4a(n) == 1) then
               qold_nh4(n) = q(i,k,lptr_nh4_a_amode(n)-loffset)
            else
               qold_nh4(n) = 0.0_r8
            end if
         end do ! n

!   compute TMR tendencies for so4 interstial aerosol due to simple gas uptake
!   sum_dqdt_so4 = so4_a tendency from h2so4 gas uptake (mol/mol/s)
!   sum_uprt_so4 = h2so4 uptake rate summed over modes (1/s)
         call nonvolatile_uptake_cell( ntot_amode, deltat, uptkrate(:,i,k), gas(igas_h2so4)%diffus_ratio, &
              so4_on, q(i,k,l_so4g), dqdt_so4, sum_dqdt_so4, sum_uprt_so4 )

         if ( associated(sulfeq) .and. (k <= troplev(i)) ) then
!   compute TMR tendencies for so4 interstial aerosol due to reversible gas uptake
!   only above the tropopause (replaces the simple gas uptake tendencies)
            call h2so4_reversible_uptake_cell( ntot_amode, deltat, uptkrate(:,i,k), so4_on, &
                 q(i,k,l_so4g), qold_so4, sulfeq(i,k,:), dqdt_so4, sum_dqdt_so4 )
         end if

         !   compute TMR tendencies for nh4 interstial aerosol due to simple gas uptake
         !   but force nh4/so4 molar ratio <= 2
!   sum_dqdt_nh4 = nh4_a tendency from nh3   gas uptake (mol/mol/s)
         sum_dqdt_nh4_b = 0.0_r8
         dqdt_nh4(:) = 0._r8
         if ( do_nh4g ) then
            call nonvolatile_uptake_cell( ntot_amode, deltat, uptkrate(:,i,k), gas(igas_nh3)%diffus_ratio, &
                 nh4_on, q(i,k,l_nh4g), dqdt_nh4, sum_dqdt_nh4, sum_uprt_nh4 )
            do n = 1, ntot_amode
               qnew_nh4 = qold_nh4(n) + dqdt_nh4(n)*deltat
               qnew_so4 = qold_so4(n) + dqdt_so4(n)*deltat
               qmax_nh4 = 2.0_r8*qnew_so4
               if (qnew_nh4 > qmax_nh4) then
                  dqdt_nh4(n) = (qmax_nh4 - qold_nh4(n))/deltatxx
               end if
               sum_dqdt_nh4_b = sum_dqdt_nh4_b + dqdt_nh4(n)
            end do
         end if

         if ( do_soag_any ) then
!   soa tendencies computed above by soa_gasaerexch_column
            dqdt_soa(:,:) = dqdt_soa_c(:,:,i,k)
            sum_dqdt_soa(:) = sum_dqdt_soa_c(:,i,k)
         else ! no soa gas species, no uptake
            dqdt_soa(:,:) = 0.0_r8
            sum_dqdt_soa(:) = 0.0_r8
         end if

         pdel_fac = pdel(i,k)/gravit
         do n = 1, ntot_amode
            if (ido_so4a(n) == 1) then
               l = lptr_so4_a_amode(n)-loffset
               dqdt(i,k,l) = dqdt_so4(n)
               qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(n)*pdel_fac
            end if

            if ( do_nh4g ) then
               if (ido_nh4a(n) == 1) then
                  l = lptr_nh4_a_amode(n)-loffset
                  dqdt(i,k,l) = dqdt_nh4(n)
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(n)*pdel_fac
               end if
            end if

            do jsoa = 1, nsoa
               if ( do_soag(jsoa) ) then
                  if (ido_soaa(n,jsoa) == 1) then
                     l = lptr2_soa_a_amode(n,jsoa)-loffset
                     dqdt(i,k,l) = dqdt_soa(n,jsoa) !calculated by  modal_aero_soaexch for method_soa=2
                     qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(n,jsoa)*pdel_fac
                  end if
               end if
            end do
         end do ! n

!   compute TMR tendencies for h2so4 and nh3 gas
!   due to simple gas uptake
         l = l_so4g
         dqdt(i,k,l) = -sum_dqdt_so4
         qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac

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
            n = modeptr_pcarbon
            tmpa = 0.0_r8
            do jsoa = 1, nsoa
               tmpa = tmpa + dqdt_soa(n,jsoa)*fac_m2v_soa(jsoa)*soa_equivso4_factor(jsoa)
            end do
            vol_shell = deltat *   &
                 ( dqdt_so4(n)*fac_m2v_so4 + dqdt_nh4(n)*fac_m2v_nh4 + tmpa )
            vol_core = 0.0_r8
            do l = 1, nspec_amode(n)
               vol_core = vol_core + &
                    q(i,k,lmassptr_amode(l,n)-loffset)*fac_m2v_pcarbon(l)
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
                  l = lptr_so4_a_amode(modetoo_pcage)-loffset
                  dqdt(i,k,l) = dqdt(i,k,l) + dqdt_so4(modefrm_pcage)
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(modefrm_pcage)*pdel_fac
               end if

               if (ido_nh4a(modetoo_pcage) > 0) then
                  l = lptr_nh4_a_amode(modetoo_pcage)-loffset
                  dqdt(i,k,l) = dqdt(i,k,l) + dqdt_nh4(modefrm_pcage)
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(modefrm_pcage)*pdel_fac
               end if

               do jsoa = 1, nsoa
                  if (ido_soaa(modetoo_pcage,jsoa) > 0) then
                     l = lptr2_soa_a_amode(modetoo_pcage,jsoa)-loffset
                     dqdt(i,k,l) = dqdt(i,k,l) + dqdt_soa(modefrm_pcage,jsoa)
                     qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(modefrm_pcage,jsoa)*pdel_fac
                  end if
               end do

            end if

         end if

      end do   ! "i = 1, ncol"
   end do     ! "k = top_lev, pver"

   if ( do_soag_any ) then
      deallocate( qold_soag_c, qold_soa_c, qold_poa_c, dqdt_soa_c, sum_dqdt_soa_c )
   end if


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
   call soa_diag%write_fields( lchnk, qcon, qevap )
   deallocate( qcon, qevap )
!-----------------------------------------------------------------------

!   do history file column-tendency fields
   do l = 1, pcnstxx
      lb = l + loffset
      do jsrf = 1, 2
         do jac = 1, 2
	    if (jac == 1) then
               if (jsrf == jsrflx_gaexch) then
                  if ( .not. dotend(l) ) cycle
                  fieldname = trim(cnst_name(lb)) // '_sfgaex1'
               else if (jsrf == jsrflx_rename) then
                  if ( .not. dotendrn(l) ) cycle
                  fieldname = trim(cnst_name(lb)) // '_sfgaex2'
               else
                  cycle
               end if
               do i = 1, ncol
                  qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf)*(adv_mass(l)/mwdry)
               end do
               call outfld( fieldname, qsrflx(:,l,jsrf), pcols, lchnk )
	    else
               if (jsrf == jsrflx_gaexch) then
                  cycle
               else if (jsrf == jsrflx_rename) then
                  if ( .not. dotendqqcwrn(l) ) cycle
                  fieldname = trim(cnst_name_cw(lb)) // '_sfgaex2'
               else
                  cycle
               end if
               do i = 1, ncol
                  qqcwsrflx(i,l,jsrf) = qqcwsrflx(i,l,jsrf)*(adv_mass(l)/mwdry)
               end do
               call outfld( fieldname, qqcwsrflx(:,l,jsrf), pcols, lchnk )
	    end if
         end do ! jac = ...
      end do ! jsrf = ...
   end do ! l = ...

   return
   end subroutine modal_aero_gasaerexch_sub


!----------------------------------------------------------------------
!----------------------------------------------------------------------

!----------------------------------------------------------------------

!----------------------------------------------------------------------

!----------------------------------------------------------------------

end module modal_aero_gasaerexch
