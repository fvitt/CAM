!----------------------------------------------------------------------
! aero_soaexch_mod --- SOA gas-aerosol exchange shared by the modal (MAM)
!                      and sectional (CARMA) aerosol models
!----------------------------------------------------------------------
module aero_soaexch_mod

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: soa_gasaerexch_cell
  public :: aero_soaexch

  ! method_soa=0 is no uptake
  ! method_soa=1 is irreversible uptake done like h2so4 uptake
  ! method_soa=2 is reversible uptake using subr aero_soaexch
  integer, parameter, public :: method_soa = 2

  ! max allowed number of sub-step iterations in aero_soaexch
  integer, parameter, public :: soaexch_niter_max = 1000

  ! ratio of soa to h2so4 gas diffusivity applied to the h2so4 uptake rate
  ! (differences in fuchs-sutugin and accom coef ignored)
  real(r8), parameter, public :: soa_h2so4_diffus_ratio = 0.81_r8

contains

!----------------------------------------------------------------------
subroutine soa_gasaerexch_cell( deltat, temp, pres, nbins, nsoa, npoa, &
     soa_on, do_soaexch, do_soag, uptkrate, mw_poa_host, mw_soa_host, &
     qold_soag, qold_soa, qold_poa, dqdt_soa, sum_dqdt_soa )

!-----------------------------------------------------------------------
!
! Purpose:
!
! computes soa gas and soa aerosol mixing ratio tendencies due to
! gas-aerosol exchange in one grid cell, for all aerosol modes or bins
!
! the host model packs the soa gas, soa aerosol and poa aerosol mixing ratios
! (mol/mol-air at the host molecular weights) and applies the tendencies
!
!-----------------------------------------------------------------------

  real(r8), intent(in)  :: deltat                     ! time step (s)
  real(r8), intent(in)  :: temp                       ! air temperature (K)
  real(r8), intent(in)  :: pres                       ! air pressure (Pa)
  integer,  intent(in)  :: nbins                      ! number of aerosol modes or bins
  integer,  intent(in)  :: nsoa                       ! number of soa species
  integer,  intent(in)  :: npoa                       ! number of poa species
  logical,  intent(in)  :: soa_on(nbins,nsoa)         ! true if soa gas species condenses on mode/bin
  logical,  intent(in)  :: do_soaexch(nbins)          ! true if mode/bin participates in aero_soaexch
  logical,  intent(in)  :: do_soag(nsoa)              ! true if soa gas species is present
  real(r8), intent(in)  :: uptkrate(nbins)            ! h2so4 gas-to-aerosol mass transfer rate (1/s)
  real(r8), intent(in)  :: mw_poa_host                ! molec wght of poa used in host code
  real(r8), intent(in)  :: mw_soa_host                ! molec wght of soa used in host code
  real(r8), intent(in)  :: qold_soag(nsoa)            ! soa gas mixrat (mol/mol)
  real(r8), intent(in)  :: qold_soa(nbins,nsoa)       ! soa aerosol mixrat (mol/mol)
  real(r8), intent(in)  :: qold_poa(nbins,npoa)       ! poa aerosol mixrat (mol/mol)
  real(r8), intent(out) :: dqdt_soa(nbins,nsoa)       ! soa aerosol mixrat tendency (mol/mol/s)
  real(r8), intent(out) :: sum_dqdt_soa(nsoa)         ! soa aerosol tendency summed over modes/bins
                                                      ! = -(soa gas mixrat tendency) (mol/mol/s)

  integer  :: jsoa, n, niter
  real(r8) :: avg_uprt_soa(nsoa)
  real(r8) :: deltatxx
  real(r8) :: dqdt_soag(nsoa)
  real(r8) :: fgain_soa(nbins,nsoa)
  real(r8) :: sum_uprt_soa(nsoa)          ! total soa uptake rate over all modes/bins
  real(r8) :: uptkrate_soa(nbins,nsoa)    ! soa gas-to-aerosol mass transfer rates (1/s)

! use this for tendency calcs to avoid generating very small negative values
  deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)

!   fgain_soa(n,jsoa) = fraction of total soa uptake going to mode/bin n
  sum_uprt_soa(:) = 0.0_r8
  do n = 1, nbins
     do jsoa = 1, nsoa
        if (soa_on(n,jsoa)) then
           fgain_soa(n,jsoa) = uptkrate(n)*soa_h2so4_diffus_ratio
           sum_uprt_soa(jsoa) = sum_uprt_soa(jsoa) + fgain_soa(n,jsoa)
        else
           fgain_soa(n,jsoa) = 0.0_r8
        end if
        uptkrate_soa(n,jsoa) = fgain_soa(n,jsoa)
     end do
  end do

  do jsoa = 1, nsoa
     if (sum_uprt_soa(jsoa) > 0.0_r8) then
        do n = 1, nbins
           fgain_soa(n,jsoa) = fgain_soa(n,jsoa) / sum_uprt_soa(jsoa)
        end do
     end if
  end do

!   uptake amount (fraction of gas uptaken) over deltat
  do jsoa = 1, nsoa
     avg_uprt_soa(jsoa) = (1.0_r8 - exp(-deltatxx*sum_uprt_soa(jsoa)))/deltatxx
  end do

!   sum_dqdt_soa = soa_a tendency from soa gas uptake (mol/mol/s)
  do jsoa = 1, nsoa
     if (do_soag(jsoa)) then
        sum_dqdt_soa(jsoa) = qold_soag(jsoa) * avg_uprt_soa(jsoa)
     else
        sum_dqdt_soa(jsoa) = 0.0_r8
     end if
  end do

  if (method_soa > 1) then
!   compute TMR tendencies for soag and soa interstial aerosol
!   using soa parameterization
     dqdt_soa(:,:) = 0.0_r8
     dqdt_soag(:) = 0.0_r8

     call aero_soaexch( deltat, temp, pres, &
          niter, soaexch_niter_max, nbins, nsoa, npoa, &
          do_soaexch, &
          mw_poa_host, mw_soa_host, &
          qold_soag, qold_soa, qold_poa, uptkrate_soa, &
          dqdt_soag, dqdt_soa )

     sum_dqdt_soa(:) = -dqdt_soag(:)

  else if (method_soa == 1) then
!   compute TMR tendencies for soa interstial aerosol
!   due to simple gas uptake
     do jsoa = 1, nsoa
        do n = 1, nbins
           if (soa_on(n,jsoa)) then
              dqdt_soa(n,jsoa) = fgain_soa(n,jsoa)*sum_dqdt_soa(jsoa)
           else
              dqdt_soa(n,jsoa) = 0.0_r8
           end if
        end do
     end do

  else
!   no uptake
     dqdt_soa(:,:) = 0.0_r8
     sum_dqdt_soa(:) = 0.0_r8
  end if

end subroutine soa_gasaerexch_cell

!----------------------------------------------------------------------
subroutine aero_soaexch( dtfull, temp, pres, &
     niter, niter_max, nbins, ntot_soaspec, ntot_poaspec,  &
     do_soaexch, &
     mw_poa_host, mw_soa_host, &
     g_soa_in, a_soa_in, a_poa_in, xferrate_in, &
     g_soa_tend, a_soa_tend )

!-----------------------------------------------------------------------
!
! Purpose:
!
! calculates condensation/evaporation of "soa gas"
! to/from multiple aerosol modes in 1 grid cell
!
! key assumptions
! (1) ambient equilibrium vapor pressure of soa gas
!     is given by p0_soa_298 and delh_vap_soa
! (2) equilibrium vapor pressure of soa gas at aerosol
!     particle surface is given by raoults law in the form
!     g_star = g0_soa*[a_soa/(a_soa + a_opoa)]
! (3) (oxidized poa)/(total poa) is equal to frac_opoa (constant)
!
!
! Author: R. Easter and R. Zaveri
! Additions to run with multiple BC, SOA and POM's: Shrivastava et al., 2015
!-----------------------------------------------------------------------

  use mo_constants, only: rgas ! Gas constant (J/K/mol)

  real(r8), intent(in)  :: dtfull                          ! full integration time step (s)
  real(r8), intent(in)  :: temp                            ! air temperature (K)
  real(r8), intent(in)  :: pres                            ! air pressure (Pa)
  integer,  intent(out) :: niter                           ! number of iterations performed
  integer,  intent(in)  :: niter_max                       ! max allowed number of iterations
  integer,  intent(in)  :: nbins                           ! number of bins
  integer,  intent(in)  :: ntot_poaspec                    ! number of poa species
  integer,  intent(in)  :: ntot_soaspec                    ! number of soa species
  logical,  intent(in)  :: do_soaexch(nbins)               ! true if bin participates in soa exchange
  real(r8), intent(in)  :: mw_poa_host                     ! molec wght of poa used in host code
  real(r8), intent(in)  :: mw_soa_host                     ! molec wght of soa used in host code
  real(r8), intent(in)  :: g_soa_in(ntot_soaspec)          ! initial soa gas mixrat (mol/mol at host mw)
  real(r8), intent(in)  :: a_soa_in(nbins,ntot_soaspec)    ! initial soa aerosol mixrat (mol/mol at host mw)
  real(r8), intent(in)  :: a_poa_in(nbins,ntot_poaspec)    ! initial poa aerosol mixrat (mol/mol at host mw)
  real(r8), intent(in)  :: xferrate_in(nbins,ntot_soaspec) ! gas-aerosol mass transfer rate (1/s)
  real(r8), intent(out) :: g_soa_tend(ntot_soaspec)        ! soa gas mixrat tendency (mol/mol/s at host mw)
  real(r8), intent(out) :: a_soa_tend(nbins,ntot_soaspec)  ! soa aerosol mixrat tendency (mol/mol/s at host mw)

  integer :: ll
  integer :: m

  logical :: skip_soamode(nbins)   ! true if this bin does not have soa

  real(r8), parameter :: a_min1 = 1.0e-20_r8
  real(r8), parameter :: g_min1 = 1.0e-20_r8
  real(r8), parameter :: alpha = 0.05_r8        ! parameter used in calc of time step
  real(r8), parameter :: dtsub_fixed = -1.0_r8  ! fixed sub-step for time integration (s)

  real(r8) :: a_ooa_sum_tmp(nbins)          ! total ooa (=soa+opoa) in a bin
  real(r8) :: a_opoa(nbins)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
  real(r8) :: a_soa(nbins,ntot_soaspec)     ! soa aerosol mixrat (mol/mol at actual mw)
  real(r8) :: a_soa_tmp(nbins,ntot_soaspec) ! temporary soa aerosol mixrat (mol/mol)
  real(r8) :: beta(nbins,ntot_soaspec)      ! dtcur*xferrate
  real(r8) :: delh_vap_soa(ntot_soaspec)    ! delh_vap_soa = heat of vaporization for gas soa (J/mol)
  real(r8) :: del_g_soa_tmp(ntot_soaspec)
  real(r8) :: dtcur                         ! current time step (s)
  real(r8) :: dtmax                         ! = (dtfull-tcur)
  real(r8) :: g0_soa(ntot_soaspec)          ! ambient soa gas equilib mixrat (mol/mol at actual mw)
  real(r8) :: g_soa(ntot_soaspec)           ! soa gas mixrat (mol/mol at actual mw)
  real(r8) :: g_star(nbins,ntot_soaspec)    ! soa gas mixrat that is in equilib
  ! with each aerosol mode (mol/mol)
  real(r8) :: mw_poa                        ! actual molec wght of poa
  real(r8) :: mw_soa                        ! actual molec wght of soa
  real(r8) :: opoa_frac(ntot_poaspec)       ! fraction of poa that is opoa
  real(r8) :: phi(nbins,ntot_soaspec)       ! "relative driving force"
  real(r8) :: p0_soa(ntot_soaspec)          ! soa gas equilib vapor presssure (atm)
  real(r8) :: p0_soa_298(ntot_soaspec)      ! p0_soa_298 = soa gas equilib vapor presssure (atm) at 298 k
  real(r8) :: sat(nbins,ntot_soaspec)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
  !    used by the numerical integration scheme -- it is not a saturation rato!
  real(r8) :: tcur                          ! current integration time (from 0 s)
  real(r8) :: tmpa, tmpb, tmpf
  real(r8) :: tot_soa(ntot_soaspec)         ! g_soa + sum( a_soa(:) )
  real(r8) :: xferrate(nbins,ntot_soaspec)  ! gas-aerosol mass transfer rate (1/s)

! Changed by Manish Shrivastava
  opoa_frac(:) = 0.0_r8 !POA does not form solution with SOA for all runs; set opoa_frac=0.0_r8  by Manish Shrivastava
  mw_poa = 250.0_r8
  mw_soa = 250.0_r8

  ! New SOA properties added by Manish Shrivastava on 09/27/2012
  if (ntot_soaspec ==1) then
     p0_soa_298(:) = 9.7831E-11_r8
     delh_vap_soa(:) = 131.0e3_r8
     opoa_frac(:) = 0.0_r8
  elseif (ntot_soaspec ==2) then
     ! same for anthropogenic and biomass burning species
     p0_soa_298 (1) = 1.0e-10_r8
     p0_soa_298 (2) = 1.0e-10_r8
     delh_vap_soa(:) = 156.0e3_r8
  elseif(ntot_soaspec ==5) then
     ! 5 volatility bins for each of the a combined SOA classes ( including biomass burning, fossil fuel, biogenic)
     p0_soa_298 (1) = 9.7831E-13_r8 !soaff0 C*=0.01ug/m3
     p0_soa_298 (2) = 9.7831E-12_r8 !soaff1 C*=0.10ug/m3
     p0_soa_298 (3) = 9.7831E-11_r8 !soaff2 C*=1.0ug/m3
     p0_soa_298 (4) = 9.7831E-10_r8 !soaff3 C*=10.0ug/m3
     p0_soa_298 (5) = 9.7831E-9_r8  !soaff4 C*=100.0ug/m3

     delh_vap_soa(1) = 153.0e3_r8
     delh_vap_soa(2) = 142.0e3_r8
     delh_vap_soa(3) = 131.0e3_r8
     delh_vap_soa(4) = 120.0e3_r8
     delh_vap_soa(5) = 109.0e3_r8
  elseif(ntot_soaspec ==15) then
     !
     ! 5 volatility bins for each of the 3 SOA classes ( biomass burning, fossil fuel, biogenic)
     ! SOA species 1-5 are for anthropogenic while 6-10 are for biomass burning SOA
     ! SOA species 11-15 are for biogenic SOA, based on Cappa et al., Reference needs to be updated
     ! For MW=250.0
     p0_soa_298 (1) = 9.7831E-13_r8 !soaff0 C*=0.01ug/m3
     p0_soa_298 (2) = 9.7831E-12_r8 !soaff1 C*=0.10ug/m3
     p0_soa_298 (3) = 9.7831E-11_r8 !soaff2 C*=1.0ug/m3
     p0_soa_298 (4) = 9.7831E-10_r8 !soaff3 C*=10.0ug/m3
     p0_soa_298 (5) = 9.7831E-9_r8  !soaff4 C*=100.0ug/m3
     p0_soa_298 (6) = 9.7831E-13_r8 !soabb0 C*=0.01ug/m3
     p0_soa_298 (7) = 9.7831E-12_r8 !soabb1 C*=0.10ug/m3
     p0_soa_298 (8) = 9.7831E-11_r8 !soabb2 C*=1.0ug/m3
     p0_soa_298 (9) = 9.7831E-10_r8 !soabb3 C*=10.0ug/m3
     p0_soa_298 (10) = 9.7831E-9_r8  !soabb4 C*=100.0ug/m3
     p0_soa_298 (11) = 9.7831E-13_r8 !soabg0 C*=0.01ug/m3
     p0_soa_298 (12) = 9.7831E-12_r8 !soabg1 C*=0.1ug/m3
     p0_soa_298 (13) = 9.7831E-11_r8 !soabg2 C*=1.0ug/m3
     p0_soa_298 (14) = 9.7831E-10_r8 !soabg3 C*=10.0ug/m3
     p0_soa_298 (15) = 9.7831E-9_r8  !soabg4 C*=100.0ug/m3

     !
     ! have to be adjusted to 15 species, following the numbers by Epstein et al., 2012
     !
     delh_vap_soa(1) = 153.0e3_r8
     delh_vap_soa(2) = 142.0e3_r8
     delh_vap_soa(3) = 131.0e3_r8
     delh_vap_soa(4) = 120.0e3_r8
     delh_vap_soa(5) = 109.0e3_r8
     delh_vap_soa(6) = 153.0e3_r8
     delh_vap_soa(7) = 142.0e3_r8
     delh_vap_soa(8) = 131.0e3_r8
     delh_vap_soa(9) = 120.0e3_r8
     delh_vap_soa(10) = 109.0e3_r8
     delh_vap_soa(11) = 153.0e3_r8
     delh_vap_soa(12) = 142.0e3_r8
     delh_vap_soa(13) = 131.0e3_r8
     delh_vap_soa(14) = 120.0e3_r8
     delh_vap_soa(15) = 109.0e3_r8
  endif

  !BSINGH - Initialized g_soa_tend and a_soa_tend to circumvent the undefined behavior (04/16/12)
  g_soa_tend(:)   = 0.0_r8
  a_soa_tend(:,:) = 0.0_r8
  xferrate(:,:) = 0.0_r8

  skip_soamode(:) = .not.do_soaexch(:)

  ! determine which modes have non-zero transfer rates
  !    and are involved in the soa gas-aerosol transfer
  ! for diameter = 1 nm and number = 1 #/cm3, xferrate ~= 1e-9 s-1
  do m = 1, nbins
     if (do_soaexch(m)) then
        do ll = 1, ntot_soaspec
           xferrate(m,ll) = xferrate_in(m,ll)
        end do
     end if
  end do

  ! convert incoming mixing ratios from mol/mol at the "host-code" molec. weight (12.0 in cam5)
  !    to mol/mol at the "actual" molec. weight (currently assumed to be 250.0)
  ! also
  !    force things to be non-negative
  !    calc tot_soa(ll)
  !    calc a_opoa (always slightly >0)
  do ll = 1, ntot_soaspec
     tmpf = mw_soa_host/mw_soa
     g_soa(ll) = max( g_soa_in(ll), 0.0_r8 ) * tmpf
     tot_soa(ll) = g_soa(ll)
     do m = 1, nbins
        if ( skip_soamode(m) ) cycle
        a_soa(m,ll) = max( a_soa_in(m,ll), 0.0_r8 ) * tmpf
        tot_soa(ll) = tot_soa(ll) + a_soa(m,ll)
     end do
  end do


  do m = 1, nbins
     if ( skip_soamode(m) ) cycle
     a_opoa(m) = 0.0_r8
     do ll = 1, ntot_poaspec
        a_opoa(m) = a_opoa(m) + opoa_frac(ll)*a_poa_in(m,ll)
     end do
  end do

  ! calc ambient equilibrium soa gas
  do ll = 1, ntot_soaspec
     p0_soa(ll) = p0_soa_298(ll) * &
          exp( -(delh_vap_soa(ll)/rgas)*((1.0_r8/temp)-(1.0_r8/298.0_r8)) )
     g0_soa(ll) = 1.01325e5_r8*p0_soa(ll)/pres
  end do

  ! IF mw of soa EQ 12 (as in the MAM3 default case), this has to be in
  ! should actully talk the mw from the chemistry mechanism and substitute with 12.0

  niter = 0
  tcur = 0.0_r8
  dtcur = 0.0_r8
  phi(:,:) = 0.0_r8
  g_star(:,:) = 0.0_r8

! integration loop -- does multiple substeps to reach dtfull
  time_loop: do while (tcur < dtfull-1.0e-3_r8 )

     niter = niter + 1
     if (niter > niter_max) exit

     tmpa = 0.0_r8  ! time integration parameter for all soa species
     do m = 1, nbins
        if ( skip_soamode(m) ) cycle
        a_ooa_sum_tmp(m) = a_opoa(m) + sum( a_soa(m,1:ntot_soaspec) )
     end do
     do ll = 1, ntot_soaspec
        tmpb = 0.0_r8  ! time integration parameter for a single soa species
        do m = 1, nbins
           if ( skip_soamode(m) ) cycle
           sat(m,ll) = g0_soa(ll)/max( a_ooa_sum_tmp(m), a_min1 )
           g_star(m,ll) = sat(m,ll)*a_soa(m,ll)
           phi(m,ll) = (g_soa(ll) - g_star(m,ll))/max( g_soa(ll), g_star(m,ll), g_min1 )
           tmpb = tmpb + xferrate(m,ll)*abs(phi(m,ll))
        end do
        tmpa = max( tmpa, tmpb )
     end do

     if (dtsub_fixed > 0.0_r8) then
        dtcur = dtsub_fixed
        tcur = tcur + dtcur
     else
        dtmax = dtfull-tcur
        if (dtmax*tmpa <= alpha) then
! here alpha/tmpa >= dtmax, so this is final substep
           dtcur = dtmax
           tcur = dtfull
        else
           dtcur = alpha/tmpa
           tcur = tcur + dtcur
        end if
     end if

! step 1 - for modes where soa is condensing, estimate "new" a_soa(m,ll)
!    using an explicit calculation with "old" g_soa
!    and g_star(m,ll) calculated using "old" a_soa(m,ll)
! do this to get better estimate of "new" a_soa(m,ll) and sat(m,ll)
     do m = 1, nbins
        if ( skip_soamode(m) ) cycle
        do ll = 1, ntot_soaspec
           ! first ll loop calcs a_soa_tmp(m,ll) & a_ooa_sum_tmp
           a_soa_tmp(m,ll) = a_soa(m,ll)
           beta(m,ll) = dtcur*xferrate(m,ll)
           del_g_soa_tmp(ll) = g_soa(ll) - g_star(m,ll)
           if (del_g_soa_tmp(ll) > 0.0_r8) then
              a_soa_tmp(m,ll) = a_soa(m,ll) + beta(m,ll)*del_g_soa_tmp(ll)
           end if
        end do
        a_ooa_sum_tmp(m) = a_opoa(m) + sum( a_soa_tmp(m,1:ntot_soaspec) )
        do ll = 1, ntot_soaspec
           ! second ll loop calcs sat & g_star
           if (del_g_soa_tmp(ll) > 0.0_r8) then
              sat(m,ll) = g0_soa(ll)/max( a_ooa_sum_tmp(m), a_min1 )
              g_star(m,ll) = sat(m,ll)*a_soa_tmp(m,ll)   ! this just needed for diagnostics
           end if
        end do
     end do

! step 2 - implicit in g_soa and semi-implicit in a_soa,
!    with g_star(m,ll) calculated semi-implicitly
     do ll = 1, ntot_soaspec
        tmpa = 0.0_r8
        tmpb = 0.0_r8
        do m = 1, nbins
           if ( skip_soamode(m) ) cycle
           tmpa = tmpa + a_soa(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
           tmpb = tmpb + beta(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
        end do

        g_soa(ll) = (tot_soa(ll) - tmpa)/(1.0_r8 + tmpb)
        g_soa(ll) = max( 0.0_r8, g_soa(ll) )
        do m = 1, nbins
           if ( skip_soamode(m) ) cycle
           a_soa(m,ll) = (a_soa(m,ll) + beta(m,ll)*g_soa(ll))/   &
                (1.0_r8 + beta(m,ll)*sat(m,ll))
        end do
     end do

  end do time_loop

! calculate outgoing tendencies (at the host-code molec. weight)
! (a_soa & g_soa are at actual mw, but a_soa_in & g_soa_in are at host-code mw)
  do ll = 1, ntot_soaspec
     tmpf = mw_soa/mw_soa_host
     g_soa_tend(ll) = (g_soa(ll)*tmpf - g_soa_in(ll))/dtfull
     do m = 1, nbins
        if ( skip_soamode(m) ) cycle
        a_soa_tend(m,ll) = (a_soa(m,ll)*tmpf - a_soa_in(m,ll))/dtfull
     end do
  end do

end subroutine aero_soaexch

end module aero_soaexch_mod
