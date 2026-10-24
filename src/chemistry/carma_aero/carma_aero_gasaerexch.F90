!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
!
! !MODULE: carma_aero_gasaerexch --- does carma aerosol gas-aerosol exchange for SOA
!
! !INTERFACE:
module carma_aero_gasaerexch

! !USES:
  use shr_kind_mod,    only:  r8 => shr_kind_r8
  use chem_mods,       only:  gas_pcnst
  use ref_pres,        only:  top_lev => clim_modal_aero_top_lev
  use ppgrid,          only:  pcols, pver
  use radiative_aerosol, only: rad_aer_get_info, rad_aer_get_info_by_bin, rad_aer_get_bin_props_by_idx, &
                              rad_aer_get_info_by_bin_spec
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use aerosol_properties_mod, only: aerosol_properties
  use aerosol_instances_mod, only: aerosol_instances_get_props, &
       aerosol_instances_get_state, aerosol_instances_get_num_models

  implicit none
  private
  public :: carma_aero_gasaerexch_sub
  public :: carma_aero_gasaerexch_init

  !PUBLIC DATA MEMBERS:

  ! description of bin aerosols
  integer, public, protected :: nspec_max = 0
  integer, public, protected :: nbins = 0
  integer, public, protected :: nsoa_vbs = 0
  integer, public, protected :: nsoa = 0
  integer, public, protected :: npoa = 0
  integer, public, protected, allocatable :: nspec(:)

  ! Misc private data
  character(len=32), allocatable :: fldname(:)    ! names for interstitial output fields
  character(len=32), allocatable :: fldname_cw(:)    ! names for cloud_borne output fields

  ! local indexing for bins
  integer :: ncnst_tot                  ! total number of mode number conc + mode species
  class(aerosol_properties), pointer :: aero_props =>null()
  integer :: iaermod_ = -1

  real(r8) :: mw_soa = 250._r8
  integer :: fracvbs_idx = -1
  integer, allocatable :: dqdtsoa_idx(:,:)
  integer, allocatable  :: cnsoa(:)         ! true if soa gas is a species and carma soa in bin
  integer, allocatable  :: cnpoa(:)         ! true if soa gas is a species and carma soa in bin
  integer, allocatable  :: l_soag(:)         ! true if soa gas is a species and carma soa in bin

  logical, allocatable  :: do_soag_any(:)         ! true if soa gas is a species and carma soa in bin
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

  subroutine carma_aero_gasaerexch_init

!-----------------------------------------------------------------------
!
! Purpose:
!      gas-aerosol exchange SOAG <-> soa
!
! Author: Simone Tilmes
!
!-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default, fieldname_len, horiz_only
    use constituents,   only: pcnst, cnst_name
    use phys_control,   only: phys_getopts
    use mo_chem_utls,   only: get_spc_ndx

!-----------------------------------------------------------------------
! arguments

!-----------------------------------------------------------------------
! local
    integer  :: j
    integer  :: i, ii
    integer  :: l
    integer  :: m
    integer  :: ns
    character(len=fieldname_len+3) :: fieldname
    character(len=32)              :: spectype
    character(len=32)              :: spec_name
    character(128)                 :: long_name
    character(8)                   :: unit
    character(len=2)               :: outsoa

    logical                        :: history_aerosol      ! Output the MAM aerosol tendencies
    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out        = history_aerosol   )

    !
    ! get info about the bin aerosols
    ! get nbins

    call rad_aer_get_info( 0, nbins=nbins)

    allocate( nspec(nbins) )
    allocate( cnsoa(nbins) )
    allocate( cnpoa(nbins) )

    do m = 1, nbins
       call rad_aer_get_info_by_bin(0, m, nspec=nspec(m))
    end do

    nspec_max = maxval(nspec)

    do iaermod_ = 1, aerosol_instances_get_num_models()
       aero_props => aerosol_instances_get_props(iaermod_, 0)
       if (aero_props%model_is('CARMA')) exit
    end do

    ncnst_tot = aero_props%ncnst_tot()

    allocate(  do_soag_any(nbins),       &
               fldname_cw(ncnst_tot),    &
               fldname(ncnst_tot) )

    ! SOAG / SOA / POM information
    ! Define number of VBS bins (nsoa) based on number of SOAG chemistry species

    nsoa_vbs = 0
    do i = 1, pcnst
       if (cnst_name(i)(:4) == 'SOAG') then
          nsoa_vbs = nsoa_vbs + 1
       end if
    end do
    allocate( l_soag(nsoa_vbs) )
    nsoa_vbs = 0
    do i = 1, pcnst
       if (cnst_name(i)(:4) == 'SOAG') then
          nsoa_vbs = nsoa_vbs + 1
          l_soag(nsoa_vbs) = get_spc_ndx(cnst_name(i))
       end if
    end do

    fracvbs_idx = pbuf_get_index('FRACVBS')

    ! identify number of SOA and POA in CARMA code (CARMA number cn)
    do m = 1, nbins
       cnsoa(m) = 0
       cnpoa(m) = 0
       do l = 1, nspec(m)
          call rad_aer_get_bin_props_by_idx(0, m, l,spectype=spectype)
          if (trim(spectype) == 's-organic') then
             cnsoa(m) = cnsoa(m) + 1
          end if
          if (trim(spectype) == 'p-organic') then
             cnpoa(m) = cnpoa(m) + 1
          end if
       end do
    end do
    ! some bins don't contain soa or poa
    nsoa= maxval(cnsoa)
    npoa= maxval(cnpoa)

    allocate( dqdtsoa_idx(nbins,nsoa)                       )
    do m = 1, nbins
       ns = 0
       do l = 1, nspec(m)
          call rad_aer_get_bin_props_by_idx(0, m, l,spectype=spectype)
          if (trim(spectype) == 's-organic') then
             call rad_aer_get_info_by_bin_spec(0, m, l, spec_name=spec_name)
             ns = ns + 1
             dqdtsoa_idx(m,ns) = pbuf_get_index('DQDT_'//trim(spec_name))
          end if
       end do
    end do

    do m = 1, nbins
       do_soag_any(m) = cnsoa(m)>0
    end do

!---------define history fields for new cond/evap diagnostics----------------------------------------

    fieldname=trim('qcon_gaex')
    long_name = trim('3D fields for SOA condensation')
    unit = 'kg/kg/s'
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
    if ( history_aerosol ) then
       call add_default( fieldname,  1, ' ' )
    endif

    do j = 1, nsoa_vbs
       write (outsoa, "(I2.2)") j
       fieldname = 'qcon_gaex'//outsoa
       long_name = '3D fields for SOA condensation for VBS bin'//outsoa
       call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif
       fieldname = 'qevap_gaex'//outsoa
       long_name = '3D fields for SOA evaporation for VBS bin'//outsoa
       call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif
    end do

    fieldname=trim('qevap_gaex')
    long_name = trim('3D fields for SOA evaporation')
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
    if ( history_aerosol ) then
       call add_default( fieldname,  1, ' ' )
    endif

!------------------------------------------------------------------------------

!  define history fields for basic gas-aer exchange
    do m = 1, nbins
       do l = 1, nspec(m)    ! do through nspec
          ii = aero_props%indexer(m,l)
          if (l <= nspec(m) ) then   ! species
             call rad_aer_get_info_by_bin_spec(0, m, l, spec_name=fldname(ii) )
             ! only write out SOA exchange here
             call rad_aer_get_bin_props_by_idx(0, m, l,spectype=spectype)
             if (trim(spectype) == 's-organic') then
                fieldname= trim(fldname(ii)) // '_sfgaex1'
                long_name = trim(fldname(ii)) // ' gas-aerosol-exchange primary column tendency'
                unit = 'kg/m2/s'
                call addfld( fieldname, horiz_only, 'A', unit, long_name )
                if ( history_aerosol ) then
                   call add_default( fieldname, 1, ' ' )
                endif
             end if
          end if
       end do

       write(fieldname,'("WETRAD_bin",I2.2)') m
       write(long_name,'("bin ",I2.2," wet radius in carma_aero_gasaerexch")') m

       call addfld(fieldname, (/'lev'/), 'A', 'cm', long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif

       write(fieldname,'("UPTKRATE_bin",I2.2)') m
       write(long_name,'("bin ",I2.2," up take rate in carma_aero_gasaerexch")') m

       call addfld(fieldname, (/'lev'/), 'A', 'sec-1', long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif

       write(fieldname,'("NUMDENS_bin",I2.2)') m
       write(long_name,'("bin ",I2.2," number density carma_aero_gasaerexch")') m

       call addfld(fieldname, (/'lev'/), 'A', 'm-3', long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif

    end do

    fieldname=trim('UPTKRATE')
    long_name = trim('total uptake rate in carma_aero_gasaerexch')
    call addfld(fieldname, (/'lev'/), 'A', 'sec-1', long_name )
    if ( history_aerosol ) then
       call add_default( fieldname,  1, ' ' )
    endif


  end subroutine carma_aero_gasaerexch_init


!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  carma_aero_gasaerexch_sub --- ...
!
! !INTERFACE:
subroutine carma_aero_gasaerexch_sub(  state, &
                        pbuf, lchnk,    ncol,     nstep,         &
                        loffset,  deltat,   mbar,                &
                        t,        pmid,     pdel,                &
                        qh2o,               troplev,             &
                        q,                  raervmr              )

  ! !USES:
  use cam_history,       only: outfld, fieldname_len
  use physconst,         only: gravit, mwdry
  use cam_abortutils,    only: endrun
  use time_manager,      only: is_first_step
  use aerosol_state_mod, only: aerosol_state
  use physics_types,     only: physics_state
  use physconst, only: mwdry, rair
  use aero_gasaerexch,   only: aero_gasaerexch_uptkrates, aero_gasaerexch_soa, &
                               method_soa_reversible

! !PARAMETERS:
  type(physics_state), target, intent(in) :: state    ! Physics state variables
  type(physics_buffer_desc), pointer :: pbuf(:)

  integer,  intent(in)    :: lchnk                ! chunk identifier
  integer,  intent(in)    :: ncol                 ! number of atmospheric column
  integer,  intent(in)    :: nstep                ! model time-step number
  integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
  integer,  intent(in)    :: troplev(pcols)       ! tropopause vertical index
  real(r8), intent(in)    :: deltat               ! time step (s)
  real(r8), intent(in)    :: mbar(ncol,pver)      ! mean wet atmospheric mass ( amu )

  real(r8), intent(inout) :: q(ncol,pver,gas_pcnst) ! tracer mixing ratio (TMR) array
  ! *** MUST BE  #/kmol-air for number
  ! *** MUST BE mol/mol-air for mass
  ! *** NOTE ncol dimension
  real(r8), intent(in)    :: raervmr (ncol,pver,ncnst_tot) ! aerosol mixing rations (vmr)
  real(r8), intent(in)    :: t(pcols,pver)        ! temperature at model levels (K)
  real(r8), intent(in)    :: pmid(pcols,pver)     ! pressure at model levels (Pa)
  real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
  real(r8), intent(in)    :: qh2o(pcols,pver)     ! water vapor mixing ratio (kg/kg)

! !DESCRIPTION:
! this version does only do condensation for SOA for CARMA
!     method_soa=0 is no uptake
!     method_soa=1 is irreversible uptake done like h2so4 uptake
!     method_soa=2 is reversible uptake using subr aero_soaexch
!     (see aero_gasaerexch)
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
  integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1, ldiag4=-1
  integer, parameter :: method_soa = method_soa_reversible

  real (r8) :: mw_poa_host(npoa)      ! molec wght of poa used in host code
  real (r8) :: mw_soa_host(nsoa_vbs)  ! molec wght of soa used in host code

  integer :: i
  integer :: j, jsoa
  integer :: k
  integer :: l
  integer :: mm, m, n, nn

  character(len=fieldname_len+3) :: fieldname
  character(len=32)              :: spectype
  character(len=2)               :: outsoa

  real (r8) :: dqdt_soa_vbs(nbins,nsoa_vbs)
  real (r8) :: dqdt_soa_all(pcols,pver,nbins,nsoa)
  real (r8) :: pdel_fac
  real (r8) :: num_bin(pcols,pver,nbins)
  real (r8) :: soa_vbs(pcols,pver,nbins,nsoa_vbs)
  real (r8) :: soa_c(pcols,pver,nbins,nsoa) ! SOA from CARMA
  real (r8) :: poa_c(pcols,pver,nbins,npoa) ! POA from CARMA
  real (r8) :: qold_poa(nbins,npoa)         ! POA from CARMA old
  real (r8) :: qold_soa(nbins,nsoa_vbs)     ! SOA on VBS bins  old
  real (r8) :: qnew_soa_vbs(nbins,nsoa_vbs) ! SOA on VBS bins  new
  real (r8) :: qnew_soa(nbins)              ! SOA new for combined VBS bin new for combined VBS binss
  real (r8) :: qold_soag(nsoa_vbs)
  real (r8) :: sum_dqdt_soa(nsoa_vbs)     !   sum_dqdt_soa = soa tendency from soa  gas uptake (mol/mol/s)
  real (r8) :: uptkrate(pcols,pver,nbins) ! gas-to-aerosol mass transfer rates (1/s)
  real (r8) :: wetdiam(pcols,pver,nbins)  ! wet diameter of each bin (m)
  real (r8) :: uptkrate_all(pcols,pver)
  real (r8) :: uptkratebb(nbins)
  logical   :: do_soa(nbins,nsoa_vbs)     ! soa species takes up gas in bin
  logical   :: do_soag(nsoa_vbs)          ! soa gas species is present

  integer, parameter :: nsrflx = 1 ! only one dimension of qsrflx, no renaming or changes in size for CARMA currently
  real(r8) :: dqdt(ncol,pver,gas_pcnst)  ! TMR "delta q" array - NOTE dims
  real(r8) :: qsrflx(pcols,nbins,nsoa)
  ! process-specific column tracer tendencies
  ! (1=gas condensation)
  real(r8) :: qcon_vbs(pcols,pver,nsoa_vbs)
  real(r8) :: qevap_vbs(pcols,pver,nsoa_vbs)
  real(r8) :: qcon(pcols,pver)
  real(r8) :: qevap(pcols,pver)
  real(r8) :: total_soag
  real(r8) :: soag(nsoa_vbs)

  real(r8), pointer :: frac_vbs(:,:,:,:)   ! fraction of vbs SOA bins to total SOA
  real(r8), pointer :: dqdt_soa(:,:)

  real(r8) :: rhoair(pcols,pver)
  real(r8), pointer :: nmr(:,:)
  class(aerosol_state), pointer :: aero_state

!----------------------------------------------------------------------
   aero_state => aerosol_instances_get_state(iaermod_, 0, lchnk)

  mw_poa_host(:) = 12.0_r8
  mw_soa_host(:) = 250.0_r8

!  map CARMA soa to working soa(nbins,nsoa)

  call pbuf_get_field(pbuf, fracvbs_idx, frac_vbs)

  num_bin(:,:,:) = 0._r8
  soa_c(:,:,:,:) = 0._r8
  poa_c(:,:,:,:) = 0._r8

  rhoair(:ncol,:) = pmid(:ncol,:)/(rair*t(:ncol,:))   ! (kg-air/m3)

  do m = 1, nbins      ! main loop over aerosol bins
     if (do_soag_any(m)) then  ! only bins that contain soa
        n = 0
        nn = 0
        do l = 1, nspec(m)
           mm = aero_props%indexer(m,l)
           call rad_aer_get_bin_props_by_idx(0, m, l, spectype=spectype)
           if (trim(spectype) == 's-organic') then
              n = n + 1
              soa_c(:ncol,:,m,n) = raervmr(:ncol,:,mm)
           end if
           if (trim(spectype) == 'p-organic') then
              nn = nn + 1
              poa_c(:ncol,:,m,nn) = raervmr(:ncol,:,mm)
           end if
        end do
        if (npoa .gt. 1) then
           call endrun( 'carma_aero_gasaerexch_sub error: CARMA currently only supports 1 POA element' )
        end if

        if (nsoa_vbs.eq.nsoa) then
           soa_vbs(:ncol,:,:,:) = soa_c(:ncol,:,:,:)
        else
           if (nsoa.eq.1) then
              if (is_first_step()) then
                 !first time step initialization only
                 do k=top_lev,pver
                    do i=1,ncol
                       total_soag = 0.0_r8
                       do j = 1, nsoa_vbs
                          soag(j) = q(i,k,l_soag(j))
                          total_soag = total_soag + soag(j)
                       end do
                       if (total_soag .gt. 0.0_r8) then
                          do j= 1, nsoa_vbs
                             frac_vbs(i,k,m,j) = soag(j)/total_soag
                          end do
                       end if
                    end do
                 end do
              end if
              ! end first time step, after that use fraction from previous time step
              do k=top_lev,pver
                 do i=1,ncol
                    do j= 1, nsoa_vbs
                       soa_vbs(i,k,m,j) = frac_vbs(i,k,m,j)*soa_c(i,k,m,nsoa)
                    end do
                 end do
              end do
           else
              ! error message this code only works if SOAG and SOA CARMA have the same number of species,
              ! or if SOA CARMA has only one species.
              call endrun( 'carma_aero_gasaerexch_sub error in number of SOA species' )
           end if

        end if

        ! get bin number densities for all aerosols
        call aero_state%get_ambient_num(m,nmr) !  #/kg
        num_bin(:ncol,:,m) = nmr(:ncol,:)*rhoair(:ncol,:) ! #/m3

     end if
  end do


! SOA will be updated in CARMA

! zero out tendencies and other
  dqdt(:,:,:) = 0.0_r8
  qsrflx(:,:,:) = 0.0_r8

!-------Initialize evap/cond diagnostics (ncols x pver)-----------
  qcon_vbs(:,:,:) = 0.0_r8
  qevap_vbs(:,:,:) = 0.0_r8
  qcon(:,:) = 0.0_r8
  qevap(:,:) = 0.0_r8
!---------------------------------------------------
! compute gas-to-aerosol mass transfer rates
  do m = 1, nbins
     wetdiam(:ncol,:,m) = aero_state%wet_diameter(m, ncol, pver)
  end do
  call aero_gasaerexch_uptkrates( aero_props, ncol, pver, top_lev, num_bin, wetdiam, &
                                  t, pmid, uptkrate )

  ! initialize to zero before summing
  uptkrate_all(:,:) = 0._r8

  do m = 1, nbins

     write(fieldname,'("NUMDENS_bin",I2.2)') m
     call outfld(fieldname, num_bin(:ncol,:,m), ncol, lchnk )

     write(fieldname,'("WETRAD_bin",I2.2)') m
     call outfld(fieldname, 0.5_r8*wetdiam(:ncol,:,m)*1.e2_r8, ncol, lchnk ) ! wet radius (cm)

     write(fieldname,'("UPTKRATE_bin",I2.2)') m
     call outfld(fieldname, uptkrate(:ncol,:,m), ncol, lchnk )

     uptkrate_all(:ncol,:) = uptkrate_all(:ncol,:) + uptkrate(:ncol,:,m)
  end do

  fieldname = trim('UPTKRATE')
  call outfld(fieldname, uptkrate_all(:ncol,:), ncol, lchnk )

  do n = 1, nbins
     do_soa(n,:) = do_soag_any(n)  ! all soa species in bins that contain soa
  end do
  do_soag(:) = .true.

  dqdt_soa_all(:,:,:,:) = 0.0_r8
  do k=top_lev,pver
     do i=1,ncol
        do n = 1, nbins
           uptkratebb(n) = uptkrate(i,k,n)
           if (do_soag_any(n)) then  ! only bins that contain soa
              do j = 1, npoa
                 qold_poa(n,j) = poa_c(i,k,n,j)
              end do
              do jsoa = 1, nsoa_vbs
                 qold_soa(n,jsoa) = soa_vbs(i,k,n,jsoa)
              end do
           else
              qold_poa(n,:) = 0.0_r8
              qold_soa(n,:) = 0.0_r8
           end if
        end do ! n

        do jsoa = 1, nsoa_vbs
           qold_soag(jsoa) = q(i,k,l_soag(jsoa))
        end do

!   compute tendencies for soag gas and soa interstitial aerosol
        call aero_gasaerexch_soa( method_soa, deltat, t(i,k), pmid(i,k), &
             nbins, nsoa_vbs, npoa, do_soag_any, do_soa, do_soag, uptkratebb, &
             mw_poa_host, mw_soa_host, &
             qold_soag, qold_soa, qold_poa, &
             dqdt_soa_vbs, sum_dqdt_soa )

        !  update soa to calcuate fractions (state variables and pbuf is not updated for SOA, will be done in CARMA)
        pdel_fac = pdel(i,k)/gravit
        qnew_soa(:) =0.0_r8
        qnew_soa_vbs(:,:) =0.0_r8

        do n = 1, nbins
           if ( do_soag_any(n) ) then
              if (nsoa.eq.nsoa_vbs) then
                 do jsoa = 1, nsoa_vbs
                    qsrflx(i,n,jsoa) = qsrflx(i,n,jsoa) + dqdt_soa_vbs(n,jsoa)*pdel_fac
                    dqdt_soa_all(i,k,n,jsoa) = dqdt_soa_vbs(n,jsoa) !  sum up for different volatility bins
                 end do
              else if (nsoa.eq.1) then
                 do jsoa = 1, nsoa_vbs
                    !  sum up for different volatility bins
                    dqdt_soa_all(i,k,n,nsoa) = dqdt_soa_all(i,k,n,nsoa) + dqdt_soa_vbs(n,jsoa)
                 end do
                 do jsoa = 1, nsoa_vbs
                    qsrflx(i,n,nsoa) = qsrflx(i,n,nsoa) + dqdt_soa_vbs(n,jsoa)*pdel_fac
                    qnew_soa_vbs(n,jsoa) = qold_soa(n,jsoa) + dqdt_soa_vbs(n,jsoa)*deltat
                    qnew_soa(n) = qnew_soa(n) + qnew_soa_vbs(n,jsoa) ! derive new fraction of SOA bin contributions
                 end do
                 do jsoa = 1, nsoa_vbs
                    if (qnew_soa(n) .gt. 0.0_r8) then
                       frac_vbs(i,k,n,jsoa) = qnew_soa_vbs(n,jsoa) / qnew_soa(n)
                    end if
                 end do
              else
                 call endrun( 'carma_aero_gasaerexch_sub error' )
              end if

!------- Add code for condensation/evaporation diagnostics sum of all bin---
              do jsoa = 1, nsoa_vbs
                 if (dqdt_soa_vbs(n,jsoa).ge.0.0_r8) then
                    qcon_vbs(i,k,jsoa)=qcon_vbs(i,k,jsoa) + dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                    qcon(i,k)=qcon(i,k)+dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                 else if (dqdt_soa_vbs(n,jsoa).lt.0.0_r8) then
                    qevap_vbs(i,k,jsoa)=qevap_vbs(i,k,jsoa) + dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                    qevap(i,k)=qevap(i,k)+dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                 endif
              end do
!---------------------------------------------------------------------------------------------------------------------
           end if
        end do ! n

!   compute TMR tendencies for SAOG gas
!   due to simple gas uptake
        do jsoa = 1, nsoa_vbs
           dqdt(i,k,l_soag(jsoa)) = -sum_dqdt_soa(jsoa)
        end do

     end do   ! "i = 1, ncol"
  end do     ! "k = top_lev, pver"

!  This applies dqdt tendencies for SOAG only , soa is done in CARMA
!  apply the dqdt to update q
!
  do jsoa = 1, nsoa_vbs
     do k = top_lev, pver
        do i = 1, ncol
           q(i,k,l_soag(jsoa)) = max (q(i,k,l_soag(jsoa)) + dqdt(i,k,l_soag(jsoa))*deltat, 1.0e-40_r8)
        end do
     end do
  end do


  !-----Outfld for condensation/evaporation------------------------------
  call outfld(trim('qcon_gaex'), qcon(:,:), pcols, lchnk )
  call outfld(trim('qevap_gaex'), qevap(:,:), pcols, lchnk )
  do jsoa = 1, nsoa_vbs
     write (outsoa, "(I2.2)") jsoa
     call outfld(trim('qcon_gaex')//outsoa, qcon_vbs(:,:,jsoa), pcols, lchnk )
     call outfld(trim('qevap_gaex')//outsoa, qevap_vbs(:,:,jsoa), pcols, lchnk )
  end do
  !-----------------------------------------------------------------------
  !   do history file of column-tendency fields over SOA fields (as defined in CARMA) and set pointer
  do m = 1, nbins
     if (do_soag_any(m)) then
        j  = 0
        do l = 1, nspec(m)
           mm = aero_props%indexer(m,l)
           call rad_aer_get_bin_props_by_idx(0, m, l,spectype=spectype)
           if (trim(spectype) == 's-organic') then
              j = j + 1
              fieldname= trim(fldname(mm)) // '_sfgaex1'
              do i = 1, ncol
                 qsrflx(i,m,j) = qsrflx(i,m,j)*(mw_soa/mwdry)
              end do
              call outfld( fieldname, qsrflx(:,m,j), pcols, lchnk )

              !set pointer field
              call pbuf_get_field(pbuf, dqdtsoa_idx(m,j),  dqdt_soa )

              dqdt_soa(:ncol,:) = dqdt_soa_all(:ncol,:,m,j) *(mw_soa/mbar(:ncol,:))
           end if
        end do ! l = ...
     end if
  end do ! m = ...

  nullify(aero_state)

end subroutine carma_aero_gasaerexch_sub

!----------------------------------------------------------------------

end module carma_aero_gasaerexch
