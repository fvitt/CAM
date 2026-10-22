!----------------------------------------------------------------------
! Reversible gas-aerosol partitioning of semi-volatile SOA species,
! shared by the aerosol models (MAM, CARMA) through their gas-aerosol
! exchange drivers.
!----------------------------------------------------------------------
module aero_soaexch_mod

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   private

   public :: aero_soaexch

contains

   !----------------------------------------------------------------------
   ! calculates condensation/evaporation of "soa gas"
   ! to/from multiple aerosol modes or bins in 1 grid cell
   !
   ! key assumptions
   ! (1) ambient equilibrium vapor pressure of soa gas
   !     is given by p0_soa_298 and delh_vap_soa
   ! (2) equilibrium vapor pressure of soa gas at aerosol
   !     particle surface is given by raoults law in the form
   !     g_star = g0_soa*[a_soa/(a_soa + a_opoa)]
   ! (3) (oxidized poa)/(total poa) is equal to frac_opoa (constant)
   !
   ! Author: R. Easter and R. Zaveri
   ! Additions to run with multiple BC, SOA and POM's: Shrivastava et al., 2015
   !----------------------------------------------------------------------
   subroutine aero_soaexch( dtfull, temp, pres, &
        niter, niter_max, nbins, do_bin, ntot_poaspec, ntot_soaspec, &
        mw_poa_host, mw_soa_host, &
        g_soa_in, a_soa_in, a_poa_in, xferrate_in, &
        g_soa_tend, a_soa_tend )

      use mo_constants,   only: rgas ! Gas constant (J/K/mol)
      use cam_abortutils, only: endrun

      real(r8), intent(in)  :: dtfull                         ! full integration time step (s)
      real(r8), intent(in)  :: temp                           ! air temperature (K)
      real(r8), intent(in)  :: pres                           ! air pressure (Pa)
      integer,  intent(out) :: niter                          ! number of iterations performed
      integer,  intent(in)  :: niter_max                      ! max allowed number of iterations
      integer,  intent(in)  :: nbins                          ! number of modes or bins
      logical,  intent(in)  :: do_bin(nbins)                  ! true if mode/bin takes part in soa exchange
      integer,  intent(in)  :: ntot_poaspec                   ! number of poa species
      integer,  intent(in)  :: ntot_soaspec                   ! number of soa species
      real(r8), intent(in)  :: mw_poa_host(ntot_poaspec)      ! molec wght of poa used in host code
      real(r8), intent(in)  :: mw_soa_host(ntot_soaspec)      ! molec wght of soa used in host code
      real(r8), intent(in)  :: g_soa_in(ntot_soaspec)         ! initial soa gas mixrat (mol/mol at host mw)
      real(r8), intent(in)  :: a_soa_in(nbins,ntot_soaspec)   ! initial soa aerosol mixrat (mol/mol at host mw)
      real(r8), intent(in)  :: a_poa_in(nbins,ntot_poaspec)   ! initial poa aerosol mixrat (mol/mol at host mw)
      real(r8), intent(in)  :: xferrate_in(nbins,ntot_soaspec) ! gas-aerosol mass transfer rate (1/s)
      real(r8), intent(out) :: g_soa_tend(ntot_soaspec)       ! soa gas mixrat tendency (mol/mol/s at host mw)
      real(r8), intent(out) :: a_soa_tend(nbins,ntot_soaspec) ! soa aerosol mixrat tendency (mol/mol/s at host mw)

      integer :: ll
      integer :: m

      real(r8), parameter :: a_min1 = 1.0e-20_r8
      real(r8), parameter :: g_min1 = 1.0e-20_r8
      real(r8), parameter :: alpha = 0.05_r8        ! parameter used in calc of time step
      real(r8), parameter :: dtsub_fixed = -1.0_r8  ! fixed sub-step for time integration (s)

      real(r8) :: a_ooa_sum_tmp(nbins)          ! total ooa (=soa+opoa) in a mode/bin
      real(r8) :: a_opoa(nbins)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa(nbins,ntot_soaspec)     ! soa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa_tmp(nbins,ntot_soaspec) ! temporary soa aerosol mixrat (mol/mol)
      real(r8) :: beta(nbins,ntot_soaspec)      ! dtcur*xferrate
      real(r8) :: delh_vap_soa(ntot_soaspec)    ! heat of vaporization for gas soa (J/mol)
      real(r8) :: del_g_soa_tmp(ntot_soaspec)
      real(r8) :: dtcur                         ! current time step (s)
      real(r8) :: dtmax                         ! = (dtfull-tcur)
      real(r8) :: g0_soa(ntot_soaspec)          ! ambient soa gas equilib mixrat (mol/mol at actual mw)
      real(r8) :: g_soa(ntot_soaspec)           ! soa gas mixrat (mol/mol at actual mw)
      real(r8) :: g_star(nbins,ntot_soaspec)    ! soa gas mixrat that is in equilib
                                                ! with each aerosol mode/bin (mol/mol)
      real(r8) :: mw_poa(ntot_poaspec)          ! actual molec wght of poa
      real(r8) :: mw_soa(ntot_soaspec)          ! actual molec wght of soa
      real(r8) :: opoa_frac(ntot_poaspec)       ! fraction of poa that is opoa
      real(r8) :: phi(nbins,ntot_soaspec)       ! "relative driving force"
      real(r8) :: p0_soa(ntot_soaspec)          ! soa gas equilib vapor presssure (atm)
      real(r8) :: p0_soa_298(ntot_soaspec)      ! soa gas equilib vapor presssure (atm) at 298 k
      real(r8) :: sat(nbins,ntot_soaspec)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
                                                !    used by the numerical integration scheme -- it is not a saturation rato!
      real(r8) :: tcur                          ! current integration time (from 0 s)
      real(r8) :: tmpa, tmpb, tmpf
      real(r8) :: tot_soa(ntot_soaspec)         ! g_soa + sum( a_soa(:) )
      real(r8) :: xferrate(nbins,ntot_soaspec)  ! gas-aerosol mass transfer rate (1/s)

      ! POA does not form solution with SOA for all runs; set opoa_frac=0.0_r8  by Manish Shrivastava
      opoa_frac(:) = 0.0_r8
      mw_poa(:) = 250.0_r8
      mw_soa(:) = 250.0_r8

      call soa_vbs_params( ntot_soaspec, p0_soa_298, delh_vap_soa )

      g_soa_tend(:)   = 0.0_r8
      a_soa_tend(:,:) = 0.0_r8
      xferrate(:,:)   = 0.0_r8

      ! modes/bins involved in the soa gas-aerosol transfer
      ! for diameter = 1 nm and number = 1 #/cm3, xferrate ~= 1e-9 s-1
      do m = 1, nbins
         if ( .not. do_bin(m) ) cycle
         do ll = 1, ntot_soaspec
            xferrate(m,ll) = xferrate_in(m,ll)
         end do
      end do

      ! convert incoming mixing ratios from mol/mol at the "host-code" molec. weight
      !    to mol/mol at the "actual" molec. weight (currently assumed to be 250.0)
      ! also
      !    force things to be non-negative
      !    calc tot_soa(ll)
      do ll = 1, ntot_soaspec
         tmpf = mw_soa_host(ll)/mw_soa(ll)
         g_soa(ll) = max( g_soa_in(ll), 0.0_r8 ) * tmpf
         tot_soa(ll) = g_soa(ll)
         do m = 1, nbins
            if ( .not. do_bin(m) ) cycle
            a_soa(m,ll) = max( a_soa_in(m,ll), 0.0_r8 ) * tmpf
            tot_soa(ll) = tot_soa(ll) + a_soa(m,ll)
         end do
      end do

      ! calc a_opoa
      ! note: a_opoa is overwritten for each poa species (not summed), and with
      ! opoa_frac=0 it is always zero; kept as is to preserve answers
      do m = 1, nbins
         if ( .not. do_bin(m) ) cycle
         a_opoa(m) = 0.0_r8
         do ll = 1, ntot_poaspec
            a_opoa(m) = opoa_frac(ll)*a_poa_in(m,ll)
         end do
      end do

      ! calc ambient equilibrium soa gas
      do ll = 1, ntot_soaspec
         p0_soa(ll) = p0_soa_298(ll) * &
              exp( -(delh_vap_soa(ll)/rgas)*((1.0_r8/temp)-(1.0_r8/298.0_r8)) )
         g0_soa(ll) = 1.01325e5_r8*p0_soa(ll)/pres
      end do

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
            if ( .not. do_bin(m) ) cycle
            a_ooa_sum_tmp(m) = a_opoa(m) + sum( a_soa(m,1:ntot_soaspec) )
         end do
         do ll = 1, ntot_soaspec
            tmpb = 0.0_r8  ! time integration parameter for a single soa species
            do m = 1, nbins
               if ( .not. do_bin(m) ) cycle
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
            if ( .not. do_bin(m) ) cycle
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
               if ( .not. do_bin(m) ) cycle
               tmpa = tmpa + a_soa(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
               tmpb = tmpb + beta(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
            end do

            g_soa(ll) = (tot_soa(ll) - tmpa)/(1.0_r8 + tmpb)
            g_soa(ll) = max( 0.0_r8, g_soa(ll) )
            do m = 1, nbins
               if ( .not. do_bin(m) ) cycle
               a_soa(m,ll) = (a_soa(m,ll) + beta(m,ll)*g_soa(ll))/   &
                    (1.0_r8 + beta(m,ll)*sat(m,ll))
            end do
         end do

      end do time_loop

      ! calculate outgoing tendencies (at the host-code molec. weight)
      ! (a_soa & g_soa are at actual mw, but a_soa_in & g_soa_in are at host-code mw)
      do ll = 1, ntot_soaspec
         tmpf = mw_soa(ll)/mw_soa_host(ll)
         g_soa_tend(ll) = (g_soa(ll)*tmpf - g_soa_in(ll))/dtfull
         do m = 1, nbins
            if ( .not. do_bin(m) ) cycle
            a_soa_tend(m,ll) = (a_soa(m,ll)*tmpf - a_soa_in(m,ll))/dtfull
         end do
      end do

   contains

      !----------------------------------------------------------------------
      ! volatility parameters of the soa species
      ! (New SOA properties added by Manish Shrivastava on 09/27/2012)
      !----------------------------------------------------------------------
      subroutine soa_vbs_params( nsoa, p0_soa_298, delh_vap_soa )

         integer,  intent(in)  :: nsoa                ! number of soa species
         real(r8), intent(out) :: p0_soa_298(nsoa)    ! soa gas equilib vapor presssure (atm) at 298 k
         real(r8), intent(out) :: delh_vap_soa(nsoa)  ! heat of vaporization for gas soa (J/mol)

         integer :: j

         select case (nsoa)
         case (1)
            p0_soa_298(:) = 9.7831E-11_r8
            delh_vap_soa(:) = 131.0e3_r8
         case (2)
            ! same for anthropogenic and biomass burning species
            p0_soa_298(:) = 1.0e-10_r8
            delh_vap_soa(:) = 156.0e3_r8
         case (5, 15)
            ! 5 volatility bins for each SOA class:
            !   nsoa=5  -- one combined class (biomass burning, fossil fuel, biogenic)
            !   nsoa=15 -- species 1-5 fossil fuel (anthropogenic), 6-10 biomass burning,
            !              11-15 biogenic SOA, based on Cappa et al., Reference needs to be updated
            ! delh_vap_soa following the numbers by Epstein et al., 2012
            ! For MW=250.0
            do j = 0, nsoa-5, 5
               p0_soa_298(j+1) = 9.7831E-13_r8 ! C*=0.01ug/m3
               p0_soa_298(j+2) = 9.7831E-12_r8 ! C*=0.10ug/m3
               p0_soa_298(j+3) = 9.7831E-11_r8 ! C*=1.0ug/m3
               p0_soa_298(j+4) = 9.7831E-10_r8 ! C*=10.0ug/m3
               p0_soa_298(j+5) = 9.7831E-9_r8  ! C*=100.0ug/m3

               delh_vap_soa(j+1) = 153.0e3_r8
               delh_vap_soa(j+2) = 142.0e3_r8
               delh_vap_soa(j+3) = 131.0e3_r8
               delh_vap_soa(j+4) = 120.0e3_r8
               delh_vap_soa(j+5) = 109.0e3_r8
            end do
         case default
            call endrun('aero_soaexch_mod::soa_vbs_params: unsupported number of SOA species')
         end select

      end subroutine soa_vbs_params

   end subroutine aero_soaexch

end module aero_soaexch_mod
