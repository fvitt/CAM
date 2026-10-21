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
  use gas_aer_uptkrate_mod, only: gas_aer_uptkrates
  use aero_soaexch_mod, only: soa_gasaerexch_column
  use soa_gaex_diag_mod, only: soa_gaex_diag_t

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
  integer, allocatable  :: soa_cidx(:,:)    ! (nbins,nsoa) aerosol constituent index of the n-th soa species in bin
  integer, allocatable  :: poa_cidx(:,:)    ! (nbins,npoa) aerosol constituent index of the n-th poa species in bin
  integer, allocatable  :: l_soag(:)         ! true if soa gas is a species and carma soa in bin

  logical, allocatable  :: do_soaexch(:)

  type(soa_gaex_diag_t) :: soa_diag    ! soa condensation/evaporation diagnostics

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
    use cam_abortutils, only: endrun

!-----------------------------------------------------------------------
! arguments

!-----------------------------------------------------------------------
! local
    integer  :: i, ii
    integer  :: l
    integer  :: m
    integer  :: ns
    integer  :: np
    character(len=fieldname_len+3) :: fieldname
    character(len=32)              :: spectype
    character(len=32)              :: spec_name
    character(128)                 :: long_name
    character(8)                   :: unit
    character(len=len(cnst_name)), allocatable :: soag_names(:) ! soa gas species names

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

    allocate(  do_soaexch(nbins),       &
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
    allocate( l_soag(nsoa_vbs), soag_names(nsoa_vbs) )
    nsoa_vbs = 0
    do i = 1, pcnst
       if (cnst_name(i)(:4) == 'SOAG') then
          nsoa_vbs = nsoa_vbs + 1
          l_soag(nsoa_vbs) = get_spc_ndx(cnst_name(i))
          soag_names(nsoa_vbs) = cnst_name(i)
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

    if (npoa > 1) then
       call endrun( 'carma_aero_gasaerexch_init error: CARMA currently only supports 1 POA element' )
    end if
    ! this code only works if SOAG and SOA CARMA have the same number of species,
    ! or if SOA CARMA has only one species
    if ((nsoa > 0) .and. (nsoa_vbs /= nsoa) .and. (nsoa /= 1)) then
       call endrun( 'carma_aero_gasaerexch_init error in number of SOA species' )
    end if

    allocate( dqdtsoa_idx(nbins,nsoa)                       )
    allocate( soa_cidx(nbins,nsoa), poa_cidx(nbins,npoa)    )
    soa_cidx(:,:) = -1
    poa_cidx(:,:) = -1
    do m = 1, nbins
       ns = 0
       np = 0
       do l = 1, nspec(m)
          call rad_aer_get_bin_props_by_idx(0, m, l,spectype=spectype)
          if (trim(spectype) == 's-organic') then
             call rad_aer_get_info_by_bin_spec(0, m, l, spec_name=spec_name)
             ns = ns + 1
             dqdtsoa_idx(m,ns) = pbuf_get_index('DQDT_'//trim(spec_name))
             soa_cidx(m,ns) = aero_props%indexer(m,l)
          end if
          if (trim(spectype) == 'p-organic') then
             np = np + 1
             poa_cidx(m,np) = aero_props%indexer(m,l)
          end if
       end do
    end do

    do m = 1, nbins
       do_soaexch(m) = cnsoa(m)>0
    end do

!---------define history fields for new cond/evap diagnostics----------------------------------------

    call soa_diag%init( soag_names, history_aerosol, add_default_species=.true. )
    deallocate( soag_names )

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
subroutine carma_aero_gasaerexch_sub( &
                        pbuf, lchnk,    ncol,                    &
                        deltat,   mbar,                          &
                        t,        pmid,     pdel,                &
                        q,                  raervmr,             &
                        wetr_n                                   )

  ! !USES:
  use cam_history,       only: outfld, fieldname_len
  use physconst,         only: gravit, mwdry, rair
  use cam_abortutils,    only: endrun
  use time_manager,      only: is_first_step
  use aerosol_state_mod, only: aerosol_state

! !PARAMETERS:
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
  real(r8), intent(in)    :: wetr_n(pcols,pver,nbins) ! wet radius (cm) of each bin

! !DESCRIPTION:
! this version does only do condensation for SOA for CARMA
!     method_soa=0 is no uptake
!     method_soa=1 is irreversible uptake done like h2so4 uptake
!     method_soa=2 is reversible uptake using subr aero_soaexch
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
  real (r8), parameter :: mw_soa_host = 250.0_r8   ! molec wght of soa used in host code

  integer :: i
  integer :: j, jsoa
  integer :: k
  integer :: mm, m, n, nn

  character(len=fieldname_len+3) :: fieldname

  real (r8) :: dqdt_soa_vbs(nbins,nsoa_vbs)
  real (r8) :: dqdt_soa_all(pcols,pver,nbins,nsoa)
  real (r8) :: pdel_fac
  real (r8) :: num_bin(pcols,pver,nbins)
  real (r8) :: wetdiam(pcols,pver,nbins)    ! wet diameter (m)
  real (r8) :: soa_vbs(pcols,pver,nbins,nsoa_vbs)
  real (r8) :: soa_c(pcols,pver,nbins,nsoa) ! SOA from CARMA
  real (r8) :: poa_c(pcols,pver,nbins,npoa) ! POA from CARMA
  real (r8) :: qold_soa(nbins,nsoa_vbs)     ! SOA on VBS bins  old
  real (r8) :: qnew_soa_vbs(nbins,nsoa_vbs) ! SOA on VBS bins  new
  real (r8) :: qnew_soa(nbins)              ! SOA new for combined VBS bin new for combined VBS binss
  real (r8), allocatable :: qold_soag_c(:,:,:)    ! (nsoa_vbs,pcols,pver) soa gas mixrat
  real (r8), allocatable :: qold_soa_c(:,:,:,:)   ! (nbins,nsoa_vbs,pcols,pver) soa on VBS bins old
  real (r8), allocatable :: qold_poa_c(:,:,:,:)   ! (nbins,npoa,pcols,pver) poa from CARMA old
  real (r8), allocatable :: dqdt_soa_c(:,:,:,:)   ! (nbins,nsoa_vbs,pcols,pver) soa tendency
  real (r8), allocatable :: sum_dqdt_soa_c(:,:,:) ! (nsoa_vbs,pcols,pver) soa tendency summed over bins
  real (r8) :: diag_fac(nbins,nsoa_vbs)     ! converts soa tendency to kg/kg/s for the diagnostics
  real (r8) :: sum_dqdt_soa(nsoa_vbs)     !   sum_dqdt_soa = soa tendency from soa  gas uptake (mol/mol/s)
  real (r8) :: uptkrate(nbins,pcols,pver)
  real (r8) :: uptkrate_all(pcols,pver)
  ! gas-to-aerosol mass transfer rates (1/s)
  logical   :: soa_on(nbins,nsoa_vbs)   ! true if soa gas condenses on bin
  logical   :: do_soag(nsoa_vbs)        ! true if soa gas is a species
  logical   :: diag_on(nbins,nsoa_vbs)  ! true if soa tendency is included in the diagnostics

  real(r8) :: dqdt(ncol,pver,gas_pcnst)  ! TMR "delta q" array - NOTE dims
  real(r8) :: qsrflx(pcols,nbins,nsoa)
  ! process-specific column tracer tendencies
  ! (1=gas condensation)
  real(r8), allocatable :: qcon(:,:,:), qevap(:,:,:)  ! soa condensation/evaporation diagnostics
                                                      !    (pcols,pver,soa_diag%nfld())
  real(r8) :: total_soag
  real(r8) :: soag(nsoa_vbs)

  real(r8), pointer :: frac_vbs(:,:,:,:)   ! fraction of vbs SOA bins to total SOA
  real(r8), pointer :: dqdt_soa(:,:)

  real(r8) :: rhoair(pcols,pver)
  real(r8), pointer :: nmr(:,:)
  class(aerosol_state), pointer :: aero_state

!----------------------------------------------------------------------
   aero_state => aerosol_instances_get_state(iaermod_, 0, lchnk)

!  map CARMA soa to working soa(nbins,nsoa)

  call pbuf_get_field(pbuf, fracvbs_idx, frac_vbs)

  num_bin(:,:,:) = 0._r8
  soa_c(:,:,:,:) = 0._r8
  poa_c(:,:,:,:) = 0._r8

  rhoair(:ncol,:) = pmid(:ncol,:)/(rair*t(:ncol,:))   ! (kg-air/m3)

  do m = 1, nbins      ! main loop over aerosol bins
     if (do_soaexch(m)) then  ! only bins that contain soa
        do n = 1, cnsoa(m)
           soa_c(:ncol,:,m,n) = raervmr(:ncol,:,soa_cidx(m,n))
        end do
        do nn = 1, cnpoa(m)
           poa_c(:ncol,:,m,nn) = raervmr(:ncol,:,poa_cidx(m,nn))
        end do

        if (nsoa_vbs /= nsoa) then
           ! nsoa == 1 (checked in carma_aero_gasaerexch_init)
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
           end if

        end if

        ! get bin number densities for all aerosols
        call aero_state%get_ambient_num(m,nmr) !  #/kg
        num_bin(:ncol,:,m) = nmr(:ncol,:)*rhoair(:ncol,:) ! #/m3

     end if
  end do

  ! SOAG and SOA CARMA have the same number of species
  if (nsoa_vbs == nsoa) then
     soa_vbs(:ncol,:,:,:) = soa_c(:ncol,:,:,:)
  end if


! SOA will be updated in CARMA

! zero out tendencies and other
  dqdt(:,:,:) = 0.0_r8
  qsrflx(:,:,:) = 0.0_r8

!-------Initialize evap/cond diagnostics (ncols x pver)-----------
  allocate( qcon(pcols,pver,soa_diag%nfld()), qevap(pcols,pver,soa_diag%nfld()) )
  qcon(:,:,:) = 0.0_r8
  qevap(:,:,:) = 0.0_r8
!---------------------------------------------------
! compute gas-to-aerosol mass transfer rates
! check if only number is needed for this calculatuion!
  wetdiam(:ncol,:,:) = 2.0_r8 * wetr_n(:ncol,:,:) * 1.e-2_r8 ! wet radius (cm) -> diameter (m)
  call gas_aer_uptkrates( ncol, nbins, num_bin, wetdiam, t, pmid, uptkrate )

  ! initialize to zero before summing
  uptkrate_all(:,:) = 0._r8

  do m = 1, nbins

     write(fieldname,'("NUMDENS_bin",I2.2)') m
     call outfld(fieldname, num_bin(:ncol,:,m), ncol, lchnk )

     write(fieldname,'("WETRAD_bin",I2.2)') m
     call outfld(fieldname, wetr_n(:ncol,:,m), ncol, lchnk )

     write(fieldname,'("UPTKRATE_bin",I2.2)') m
     call outfld(fieldname, uptkrate(m,:ncol,:), ncol, lchnk )

     uptkrate_all(:ncol,:) = uptkrate_all(:ncol,:) + uptkrate(m,:ncol,:)
  end do

  fieldname = trim('UPTKRATE')
  call outfld(fieldname, uptkrate_all(:ncol,:), ncol, lchnk )

! every soa gas condenses on each bin that contains soa
  do jsoa = 1, nsoa_vbs
     soa_on(:,jsoa) = do_soaexch(:)
  end do
  do_soag(:) = .true.

! condensation/evaporation diagnostics (kg/kg/s), summed over all bins that contain soa
  do n = 1, nbins
     do jsoa = 1, nsoa_vbs
        diag_on(n,jsoa) = do_soaexch(n)
        diag_fac(n,jsoa) = mw_soa/mwdry
     end do
  end do

!   compute TMR tendencies for soag and soa interstial aerosol, for all columns
  allocate( qold_soag_c(nsoa_vbs,pcols,pver), qold_soa_c(nbins,nsoa_vbs,pcols,pver), &
            qold_poa_c(nbins,npoa,pcols,pver), dqdt_soa_c(nbins,nsoa_vbs,pcols,pver), &
            sum_dqdt_soa_c(nsoa_vbs,pcols,pver) )
  do k=top_lev,pver
     do i=1,ncol
        do n = 1, nbins
           if (do_soaexch(n)) then  ! only bins that contain soa
              do j = 1, npoa
                 qold_poa_c(n,j,i,k) = poa_c(i,k,n,j)
              end do
              do jsoa = 1, nsoa_vbs
                 qold_soa_c(n,jsoa,i,k) = soa_vbs(i,k,n,jsoa)
              end do
           else
              qold_poa_c(n,:,i,k) = 0.0_r8
              qold_soa_c(n,:,i,k) = 0.0_r8
           end if
        end do ! n

        do jsoa = 1, nsoa_vbs
           qold_soag_c(jsoa,i,k) = q(i,k,l_soag(jsoa))
        end do
     end do
  end do

!   sum_dqdt_soa = soa_a tendency from soa   gas uptake (mol/mol/s)
  call soa_gasaerexch_column( ncol, nbins, nsoa_vbs, npoa, deltat, t, pmid, &
       soa_on, do_soaexch, do_soag, uptkrate, mw_soa_host, &
       qold_soag_c, qold_soa_c, qold_poa_c, dqdt_soa_c, sum_dqdt_soa_c, &
       soa_diag, diag_on, diag_fac, qcon, qevap )

  dqdt_soa_all(:,:,:,:) = 0.0_r8
  do k=top_lev,pver
     do i=1,ncol
        qold_soa(:,:) = qold_soa_c(:,:,i,k)
        dqdt_soa_vbs(:,:) = dqdt_soa_c(:,:,i,k)
        sum_dqdt_soa(:) = sum_dqdt_soa_c(:,i,k)

        !  update soa to calcuate fractions (state variables and pbuf is not updated for SOA, will be done in CARMA)
        pdel_fac = pdel(i,k)/gravit
        qnew_soa(:) =0.0_r8
        qnew_soa_vbs(:,:) =0.0_r8

        do n = 1, nbins
           if ( do_soaexch(n) ) then
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
           end if
        end do ! n

!   compute TMR tendencies for SAOG gas
!   due to simple gas uptake
        do jsoa = 1, nsoa_vbs
           dqdt(i,k,l_soag(jsoa)) = -sum_dqdt_soa(jsoa)
        end do

     end do   ! "i = 1, ncol"
  end do     ! "k = top_lev, pver"

  deallocate( qold_soag_c, qold_soa_c, qold_poa_c, dqdt_soa_c, sum_dqdt_soa_c )

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
  call soa_diag%write_fields( lchnk, qcon, qevap )
  deallocate( qcon, qevap )
  !-----------------------------------------------------------------------
  !   do history file of column-tendency fields over SOA fields (as defined in CARMA) and set pointer
  do m = 1, nbins
     if (do_soaexch(m)) then
        do j = 1, cnsoa(m)
           mm = soa_cidx(m,j)
           fieldname= trim(fldname(mm)) // '_sfgaex1'
           do i = 1, ncol
              qsrflx(i,m,j) = qsrflx(i,m,j)*(mw_soa/mwdry)
           end do
           call outfld( fieldname, qsrflx(:,m,j), pcols, lchnk )

           !set pointer field
           call pbuf_get_field(pbuf, dqdtsoa_idx(m,j),  dqdt_soa )

           dqdt_soa(:ncol,:) = dqdt_soa_all(:ncol,:,m,j) *(mw_soa/mbar(:ncol,:))
        end do ! j = ...
     end if
  end do ! m = ...

  nullify(aero_state)

end subroutine carma_aero_gasaerexch_sub


!----------------------------------------------------------------------
!----------------------------------------------------------------------

!----------------------------------------------------------------------

!----------------------------------------------------------------------

end module carma_aero_gasaerexch
