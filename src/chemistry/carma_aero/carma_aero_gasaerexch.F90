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
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use aerosol_properties_mod, only: aerosol_properties
  use aero_gasaerexch,  only: soa_condevap_diag_t, col_tend_diag_t

  implicit none
  private
  public :: carma_aero_gasaerexch_sub
  public :: carma_aero_gasaerexch_init

  ! description of bin aerosols
  integer :: nbins = 0                 ! number of bins
  integer :: nsoa = 0                  ! number of soa gas (SOAG) species, and of soa species
                                       ! in each bin that contains soa
  integer :: npoa = 0                  ! maximum number of poa species in a bin
  integer, allocatable :: nspec(:)     ! number of species in each bin

  ! Misc private data
  character(len=32), allocatable :: fldname(:)    ! names for interstitial output fields
  character(len=32), allocatable :: fldname_cw(:)    ! names for cloud_borne output fields

  ! local indexing for bins
  integer :: ncnst_tot                  ! total number of mode number conc + mode species

  real(r8), allocatable :: mw_soa_spec(:,:) ! molecular weight of each bin soa species (from aero_props)
  real(r8), allocatable :: mw_soa_host(:)   ! molec wght of soa used in host code
  real(r8), allocatable :: mw_poa_host(:)   ! molec wght of poa used in host code
  integer, allocatable :: dqdtsoa_idx(:,:)
  integer, allocatable  :: cnsoa(:)         ! number of soa species in each bin
  integer, allocatable  :: cnpoa(:)         ! number of poa species in each bin
  integer, allocatable  :: soa_ndx(:,:)     ! species indices within each bin of the soa species
  integer, allocatable  :: poa_ndx(:,:)     ! species indices within each bin of the poa species
  integer, allocatable  :: l_soag(:)         ! chemistry indices of the soa gas species

  logical, allocatable  :: do_soag_any(:)         ! true if soa gas is a species and carma soa in bin

  type(soa_condevap_diag_t) :: soa_diag   ! SOA condensation/evaporation diagnostics
  logical,  allocatable :: condevap_mask(:,:)  ! bin soa species included in the diagnostics
  real(r8), allocatable :: condevap_fac(:,:)   ! conversion factor to mass mixing ratio (mw/mwdry)
  type(col_tend_diag_t) :: sfgaex1_diag     ! SOA column tendency diagnostics ('_sfgaex1')
  integer, allocatable :: sfgaex1_fld(:,:)  ! field index of each bin soa species
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

  subroutine carma_aero_gasaerexch_init(aero_props)

!-----------------------------------------------------------------------
!
! Purpose:
!      gas-aerosol exchange SOAG <-> soa
!
! Author: Simone Tilmes
!
!-----------------------------------------------------------------------

    use constituents,   only: cnst_name
    use phys_control,   only: phys_getopts
    use physconst,      only: mwdry
    use mo_chem_utls,   only: get_spc_ndx
    use cam_abortutils, only: endrun
    use aero_gasaerexch, only: aero_gasaerexch_uptk_diag_init, aero_gasaerexch_spec_ndx, &
                               aero_gasaerexch_soag_ndx

!-----------------------------------------------------------------------
! arguments
    class(aerosol_properties), intent(in) :: aero_props ! CARMA aerosol properties object

!-----------------------------------------------------------------------
! local
    integer  :: ii
    integer  :: j
    integer  :: l
    integer  :: m
    integer  :: ns
    integer, allocatable :: soag_cnst_ndx(:)
    real(r8) :: spec_mw
    real(r8) :: mw_soa, mw_poa
    character(len=*), parameter :: subname = 'carma_aero_gasaerexch_init'
    character(len=128) :: msg
    integer  :: nfld
    character(len=32), allocatable :: sfgaex1_names(:)
    character(len=32)              :: spec_name, spec_name_cw

    logical                        :: history_aerosol      ! Output the MAM aerosol tendencies
    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out        = history_aerosol   )

    !
    ! get info about the bin aerosols
    ! get nbins

    nbins = aero_props%nbins()

    allocate( nspec(nbins) )
    allocate( cnsoa(nbins) )
    allocate( cnpoa(nbins) )

    do m = 1, nbins
       nspec(m) = aero_props%nspecies(m)
    end do

    ncnst_tot = aero_props%ncnst_tot()

    allocate(  do_soag_any(nbins),       &
               fldname_cw(ncnst_tot),    &
               fldname(ncnst_tot) )

    ! SOAG / SOA / POM information
    ! number of SOA species (nsoa) is the number of SOAG (VBS) chemistry species

    call aero_gasaerexch_soag_ndx(soag_cnst_ndx)
    nsoa = size(soag_cnst_ndx)
    allocate( l_soag(nsoa) )
    do j = 1, nsoa
       l_soag(j) = get_spc_ndx(cnst_name(soag_cnst_ndx(j)))
    end do

    ! identify SOA and POA species in CARMA bins (CARMA number cn)
    call aero_gasaerexch_spec_ndx(aero_props, 's-organic', cnsoa, soa_ndx)
    call aero_gasaerexch_spec_ndx(aero_props, 'p-organic', cnpoa, poa_ndx)
    ! some bins don't contain soa or poa
    npoa= maxval(cnpoa)

    ! each bin that contains soa has one soa species for each soa gas (SOAG) species,
    ! with soa species j taking up soa gas species j
    do m = 1, nbins
       if (cnsoa(m) > 0) then
          if (cnsoa(m) /= nsoa) then
             write(msg,'(a,i0,a,i0,a,i0)') ': number of SOA species in bin ', m, ' (', cnsoa(m), &
                  ') differs from the number of SOAG species ', nsoa
             call endrun(subname//trim(msg))
          end if
       end if
    end do

    allocate( dqdtsoa_idx(nbins,nsoa)                       )
    do m = 1, nbins
       do ns = 1, cnsoa(m)
          call aero_props%mmr_names(m, soa_ndx(m,ns), spec_name, spec_name_cw)
          dqdtsoa_idx(m,ns) = pbuf_get_index('DQDT_'//trim(spec_name))
       end do
    end do

    do m = 1, nbins
       do_soag_any(m) = cnsoa(m)>0
    end do

    ! SOA and POA molecular weights from the aerosol properties object;
    ! the SOA solver uses a single molecular weight for SOA (POA)
    allocate( mw_soa_spec(nbins,nsoa), mw_soa_host(nsoa), mw_poa_host(npoa) )
    mw_soa_spec(:,:) = -huge(1.0_r8)
    mw_soa = -huge(1.0_r8)
    mw_poa = -huge(1.0_r8)
    do m = 1, nbins
       do ns = 1, cnsoa(m)
          call aero_props%get(m, soa_ndx(m,ns), spec_mw=spec_mw)
          if (mw_soa > 0._r8) then
             if (spec_mw /= mw_soa) then
                call endrun(subname//': s-organic species molecular weights differ')
             end if
          end if
          mw_soa_spec(m,ns) = spec_mw
          mw_soa = spec_mw
       end do
       do ns = 1, cnpoa(m)
          call aero_props%get(m, poa_ndx(m,ns), spec_mw=spec_mw)
          if (mw_poa > 0._r8) then
             if (spec_mw /= mw_poa) then
                call endrun(subname//': p-organic species molecular weights differ')
             end if
          end if
          mw_poa = spec_mw
       end do
    end do
    mw_soa_host(:) = mw_soa
    mw_poa_host(:) = mw_poa

!---------define history fields for new cond/evap diagnostics----------------------------------------
!   totals over all SOA species and per SOAG (VBS) species

    call soa_diag%init( nsoa, do_total=.true., do_spec=.true., history_aerosol=history_aerosol )
    allocate( condevap_mask(nbins,nsoa), condevap_fac(nbins,nsoa) )
    do m = 1, nbins
       condevap_mask(m,:) = do_soag_any(m)
    end do
    condevap_fac(:,:) = mw_soa/mwdry

!------------------------------------------------------------------------------

!  define history fields for basic gas-aer exchange
    allocate( sfgaex1_fld(nbins,nsoa) )
    allocate( sfgaex1_names(sum(cnsoa)) )
    sfgaex1_fld(:,:) = 0
    nfld = 0
    do m = 1, nbins
       do l = 1, nspec(m)    ! do through nspec
          ii = aero_props%indexer(m,l)
          call aero_props%mmr_names(m, l, fldname(ii), fldname_cw(ii))
       end do
       ! only write out SOA exchange here
       do ns = 1, cnsoa(m)
          ii = aero_props%indexer(m,soa_ndx(m,ns))
          nfld = nfld + 1
          sfgaex1_names(nfld) = fldname(ii)
          sfgaex1_fld(m,ns) = nfld
       end do
    end do
    call sfgaex1_diag%init( sfgaex1_names, '_sfgaex1', 'gas-aerosol-exchange primary column tendency', &
                            history_aerosol )

!  uptake rate, number density and wet radius of each bin
    call aero_gasaerexch_uptk_diag_init( nbins, history_aerosol )


  end subroutine carma_aero_gasaerexch_init


!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  carma_aero_gasaerexch_sub --- ...
!
! !INTERFACE:
subroutine carma_aero_gasaerexch_sub(  aero_props, aero_state,  &
                        pbuf, lchnk,    ncol,                    &
                        deltat,   mbar,                          &
                        t,        pmid,     pdel,                &
                        q,                  raervmr              )

  ! !USES:
  use physconst,         only: gravit, mwdry
  use cam_abortutils,    only: endrun
  use aerosol_state_mod, only: aerosol_state
  use physconst, only: mwdry, rair
  use aero_gasaerexch,   only: aero_gasaerexch_uptkrates, aero_gasaerexch_soa_column, &
                               aero_gasaerexch_uptk_diag_out, soa_condevap_acc_t, &
                               method_soa_reversible

! !PARAMETERS:
  class(aerosol_properties), intent(in) :: aero_props ! CARMA aerosol properties object
  class(aerosol_state), intent(in) :: aero_state      ! CARMA aerosol state object
  type(physics_buffer_desc), pointer :: pbuf(:)

  integer,  intent(in)    :: lchnk                ! chunk identifier
  integer,  intent(in)    :: ncol                 ! number of atmospheric column
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
  integer, parameter :: method_soa = method_soa_reversible

  integer :: i
  integer :: j, jsoa
  integer :: k
  integer :: mm, m, n, nn


  real (r8) :: dqdt_soa_bin(nbins,nsoa)
  real (r8) :: dqdt_soa_all(pcols,pver,nbins,nsoa)
  real (r8) :: pdel_fac
  real (r8) :: num_bin(pcols,pver,nbins)
  real (r8) :: soa_c(pcols,pver,nbins,nsoa) ! SOA from CARMA
  real (r8) :: poa_c(pcols,pver,nbins,npoa) ! POA from CARMA
  real (r8) :: soag_col(ncol,pver,nsoa)               ! SOAG gas
  real (r8) :: dqdt_soa_col(pcols,pver,nbins,nsoa)    ! SOA tendencies
  real (r8) :: sum_dqdt_soa_col(pcols,pver,nsoa)      ! SOAG uptake (negative of gas tendency)
  real (r8) :: sum_dqdt_soa(nsoa)         !   sum_dqdt_soa = soa tendency from soa  gas uptake (mol/mol/s)
  real (r8) :: uptkrate(pcols,pver,nbins) ! gas-to-aerosol mass transfer rates (1/s)
  real (r8) :: wetdiam(pcols,pver,nbins)  ! wet diameter of each bin (m)
  logical   :: do_soa(nbins,nsoa)         ! soa species takes up gas in bin
  logical   :: do_soag(nsoa)              ! soa gas species is present

  real(r8) :: dqdt(ncol,pver,gas_pcnst)  ! TMR "delta q" array - NOTE dims
  real(r8) :: qsrflx(pcols,nbins,nsoa)
  ! process-specific column tracer tendencies
  ! (1=gas condensation)
  type(soa_condevap_acc_t) :: condevap   ! SOA condensation/evaporation diagnostics

  real(r8), pointer :: dqdt_soa(:,:)

  real(r8) :: rhoair(pcols,pver)
  real(r8), pointer :: nmr(:,:)

!----------------------------------------------------------------------
!  map CARMA soa to working soa(nbins,nsoa)

  num_bin(:,:,:) = 0._r8
  soa_c(:,:,:,:) = 0._r8
  poa_c(:,:,:,:) = 0._r8

  rhoair(:ncol,:) = pmid(:ncol,:)/(rair*t(:ncol,:))   ! (kg-air/m3)

  do m = 1, nbins      ! main loop over aerosol bins
     if (do_soag_any(m)) then  ! only bins that contain soa
        do n = 1, cnsoa(m)
           mm = aero_props%indexer(m,soa_ndx(m,n))
           soa_c(:ncol,:,m,n) = raervmr(:ncol,:,mm)
        end do
        do nn = 1, cnpoa(m)
           mm = aero_props%indexer(m,poa_ndx(m,nn))
           poa_c(:ncol,:,m,nn) = raervmr(:ncol,:,mm)
        end do
        if (npoa .gt. 1) then
           call endrun( 'carma_aero_gasaerexch_sub error: CARMA currently only supports 1 POA element' )
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
  call soa_diag%reset(condevap)
!---------------------------------------------------
! compute gas-to-aerosol mass transfer rates
  do m = 1, nbins
     wetdiam(:ncol,:,m) = aero_state%wet_diameter(m, ncol, pver)
  end do
  call aero_gasaerexch_uptkrates( aero_props, ncol, pver, top_lev, num_bin, wetdiam, &
                                  t, pmid, uptkrate )

  call aero_gasaerexch_uptk_diag_out( lchnk, ncol, nbins, num_bin, wetdiam, uptkrate )

  do n = 1, nbins
     do_soa(n,:) = do_soag_any(n)  ! all soa species in bins that contain soa
  end do
  do_soag(:) = .true.

!   compute tendencies for soag gas and soa interstitial aerosol in all grid cells
!   (soa_c is used only in bins that contain soa; poa_c is zero in other bins)
  do jsoa = 1, nsoa
     soag_col(:ncol,:,jsoa) = q(:ncol,:,l_soag(jsoa))
  end do
  call aero_gasaerexch_soa_column( method_soa, ncol, pver, top_lev, deltat, t, pmid, &
       nbins, nsoa, npoa, do_soag_any, do_soa, do_soag, uptkrate, &
       mw_poa_host, mw_soa_host, soag_col, soa_c, poa_c, &
       dqdt_soa_col, sum_dqdt_soa_col )

  dqdt_soa_all(:,:,:,:) = 0.0_r8
  do k=top_lev,pver
     do i=1,ncol
        dqdt_soa_bin(:,:) = dqdt_soa_col(i,k,:,:)
        sum_dqdt_soa(:) = sum_dqdt_soa_col(i,k,:)

        ! SOA is not updated here (state variables and pbuf), this is done in CARMA
        pdel_fac = pdel(i,k)/gravit

        do n = 1, nbins
           if ( do_soag_any(n) ) then
              do jsoa = 1, nsoa
                 qsrflx(i,n,jsoa) = qsrflx(i,n,jsoa) + dqdt_soa_bin(n,jsoa)*pdel_fac
                 dqdt_soa_all(i,k,n,jsoa) = dqdt_soa_bin(n,jsoa)
              end do
           end if
        end do ! n

!------- condensation/evaporation diagnostics sum of all bins
        call soa_diag%accum( condevap, i, k, nbins, condevap_mask, condevap_fac, dqdt_soa_bin )

!   compute TMR tendencies for SAOG gas
!   due to simple gas uptake
        do jsoa = 1, nsoa
           dqdt(i,k,l_soag(jsoa)) = -sum_dqdt_soa(jsoa)
        end do

     end do   ! "i = 1, ncol"
  end do     ! "k = top_lev, pver"

!  This applies dqdt tendencies for SOAG only , soa is done in CARMA
!  apply the dqdt to update q
!
  do jsoa = 1, nsoa
     do k = top_lev, pver
        do i = 1, ncol
           q(i,k,l_soag(jsoa)) = max (q(i,k,l_soag(jsoa)) + dqdt(i,k,l_soag(jsoa))*deltat, 1.0e-40_r8)
        end do
     end do
  end do


  !-----Outfld for condensation/evaporation------------------------------
  call soa_diag%outfld( condevap, lchnk )
  !-----------------------------------------------------------------------
  !   do history file of column-tendency fields over SOA fields (as defined in CARMA) and set pointer
  do m = 1, nbins
     if (do_soag_any(m)) then
        do j = 1, cnsoa(m)
           call sfgaex1_diag%outfld( sfgaex1_fld(m,j), ncol, lchnk, mw_soa_spec(m,j)/mwdry, qsrflx(:,m,j) )

           !set pointer field
           call pbuf_get_field(pbuf, dqdtsoa_idx(m,j),  dqdt_soa )

           dqdt_soa(:ncol,:) = dqdt_soa_all(:ncol,:,m,j) *(mw_soa_spec(m,j)/mbar(:ncol,:))
        end do ! j = ...
     end if
  end do ! m = ...

end subroutine carma_aero_gasaerexch_sub

!----------------------------------------------------------------------

end module carma_aero_gasaerexch
