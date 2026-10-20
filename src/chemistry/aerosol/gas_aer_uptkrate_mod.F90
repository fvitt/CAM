!----------------------------------------------------------------------
! gas_aer_uptkrate_mod --- gas-to-aerosol mass transfer (uptake) rates
!                          shared by the modal (MAM) and sectional (CARMA)
!                          aerosol models
!----------------------------------------------------------------------
module gas_aer_uptkrate_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver

  implicit none
  private

  public :: gas_aer_uptkrates
  public :: nonvolatile_uptake_cell
  public :: h2so4_reversible_uptake_cell

contains

!----------------------------------------------------------------------
subroutine gas_aer_uptkrates( ncol, nbins, num, dp_in, t, pmid, uptkrate, sigmag )

!
! computes the gas-to-aerosol mass transfer rate of h2so4 gas
! for each aerosol mode or bin
!
!       gas_conden_rate(Dp) = 2 * pi * gasdiffus * Dp * F(Kn,ac)
!           F(Kn,ac) = Fuchs-Sutugin correction factor
!           Kn = Knudsen number
!           ac = accomodation coefficient
!
! log-normal modes (sigmag present):
!
!                         /
!              uptkrate = | dx  dN/dx  gas_conden_rate(Dp(x))
!                         /
!   using Gauss-Hermite quadrature of order nghq=2
!
!       Dp = particle diameter
!       x = ln(Dp)
!       dN/dx = log-normal particle number density distribution
!
!   the prefactor tworootpi = 2*sqrt(pi) combines the 2*pi of gas_conden_rate
!   with the 1/sqrt(pi) normalization of the Gauss-Hermite quadrature
!
! monodisperse bins (sigmag absent), e.g., sectional bins each treated as
! monodisperse at its wet diameter:
!
!       uptkrate = N * gas_conden_rate(Dp) = const * gasdiffus * F(Kn,ac),
!       const = 2 * pi * N * Dp
!       evaluated only where Dp > 0
!

   use ref_pres,      only: top_lev => clim_modal_aero_top_lev
   use shr_const_mod, only: pi => shr_const_pi

   integer,  intent(in)  :: ncol                   ! number of atmospheric columns
   integer,  intent(in)  :: nbins                  ! number of aerosol modes or bins
   real(r8), intent(in)  :: num(pcols,pver,nbins)  ! number concentration (#/m3)
   real(r8), intent(in)  :: dp_in(pcols,pver,nbins)! wet diameter (m) -- geometric mean diameter
                                                   ! of the number distribution for log-normal modes
   real(r8), intent(in)  :: t(pcols,pver)          ! temperature (K)
   real(r8), intent(in)  :: pmid(pcols,pver)       ! air pressure (Pa)
   real(r8), intent(out) :: uptkrate(nbins,pcols,pver) ! gas-to-aerosol mass transfer rates (1/s)
   real(r8), optional, intent(in) :: sigmag(nbins) ! geometric standard deviation of log-normal modes
                                                   ! (absent for monodisperse bins)

   integer, parameter :: nghq = 2
   integer :: i, iq, k, n

   ! Can use sqrt here once Lahey is gone.
   real(r8), parameter :: tworootpi = 3.5449077_r8
   real(r8), parameter :: twopi = 2.0_r8*pi
   real(r8), parameter :: root2 = 1.4142135_r8
   real(r8), parameter :: beta = 2.0_r8

   ! quadrature abscissae and weights
   real(r8), parameter :: xghq(nghq) = (/ 0.70710678_r8, -0.70710678_r8 /)
   real(r8), parameter :: wghq(nghq) = (/ 0.88622693_r8,  0.88622693_r8 /)

   real(r8) :: const
   real(r8) :: dp
   real(r8) :: gasdiffus, gasspeed
   real(r8) :: freepathx2
   real(r8) :: knudsen
   real(r8) :: lndp, lndpgn, lnsg
   real(r8) :: sumghq

   uptkrate(:,:,:) = 0.0_r8

   ! the branch is kept outside the loop nests so that each loop nest
   ! compiles (and rounds) as in the original model-specific routines

   if (present(sigmag)) then

      ! log-normal modes
      do n = 1, nbins
         do k = top_lev, pver
            do i = 1, ncol

               !   gasdiffus = h2so4 gas diffusivity from mosaic code (m^2/s)
               !               (pmid must be Pa)
               gasdiffus = 0.557e-4_r8 * (t(i,k)**1.75_r8) / pmid(i,k)
               !   gasspeed = h2so4 gas mean molecular speed from mosaic code (m/s)
               gasspeed  = 1.470e1_r8 * sqrt(t(i,k))
               !   freepathx2 = 2 * (h2so4 mean free path)  (m)
               freepathx2 = 6.0_r8*gasdiffus/gasspeed

               lnsg   = log( sigmag(n) )
               lndpgn = log( dp_in(i,k,n) )   ! (m)
               const  = tworootpi * num(i,k,n) * exp(beta*lndpgn + 0.5_r8*(beta*lnsg)**2)

               !   sum over gauss-hermite quadrature points
               sumghq = 0.0_r8
               do iq = 1, nghq
                  lndp = lndpgn + beta*lnsg**2 + root2*lnsg*xghq(iq)
                  dp = exp(lndp)

                  !   knudsen number
                  knudsen = freepathx2/dp
                  sumghq = sumghq + wghq(iq)*dp*fuchs_sutugin(knudsen)/(dp**beta)
               end do
               uptkrate(n,i,k) = const * gasdiffus * sumghq

            end do   ! "do i = 1, ncol"
         end do   ! "do k = top_lev, pver"
      end do   ! "do n = 1, nbins"

   else

      ! monodisperse bins
      do n = 1, nbins
         do k = top_lev, pver
            do i = 1, ncol
               if (dp_in(i,k,n) > 0.0_r8) then

                  !   gasdiffus = h2so4 gas diffusivity from mosaic code (m^2/s)
                  !               (pmid must be Pa)
                  gasdiffus = 0.557e-4_r8 * (t(i,k)**1.75_r8) / pmid(i,k)
                  !   gasspeed = h2so4 gas mean molecular speed from mosaic code (m/s)
                  gasspeed  = 1.470e1_r8 * sqrt(t(i,k))
                  !   freepathx2 = 2 * (h2so4 mean free path)  (m)
                  freepathx2 = 6.0_r8*gasdiffus/gasspeed

                  dp = dp_in(i,k,n)
                  const = twopi * num(i,k,n) * dp
                  !   knudsen number
                  knudsen = freepathx2/dp
                  uptkrate(n,i,k) = const * gasdiffus * fuchs_sutugin(knudsen)

               end if
            end do   ! "do i = 1, ncol"
         end do   ! "do k = top_lev, pver"
      end do   ! "do n = 1, nbins"

   end if

end subroutine gas_aer_uptkrates

!----------------------------------------------------------------------
subroutine nonvolatile_uptake_cell( nbins, deltat, uptkrate, diffus_ratio, mask, qgas, &
                                    dqdt_aer, sum_dqdt, sum_uprt )

!
! uptake of a non-volatile gas (gas --> aerosol, but no aerosol --> gas)
! by the aerosol modes or bins in one grid cell
!
! the gas uptake over deltat, assuming exponential decay of the gas,
! is distributed to the modes/bins in proportion to their uptake rates
!

   integer,  intent(in)  :: nbins            ! number of aerosol modes or bins
   real(r8), intent(in)  :: deltat           ! time step (s)
   real(r8), intent(in)  :: uptkrate(nbins)  ! h2so4 gas-to-aerosol mass transfer rate (1/s)
   real(r8), intent(in)  :: diffus_ratio     ! gas diffusivity relative to h2so4
   logical,  intent(in)  :: mask(nbins)      ! true if mode/bin takes up the gas
   real(r8), intent(in)  :: qgas             ! gas mixing ratio (mol/mol)
   real(r8), intent(out) :: dqdt_aer(nbins)  ! condensed species tendency in each mode/bin (mol/mol/s)
   real(r8), intent(out) :: sum_dqdt         ! condensed species tendency summed over modes/bins
                                             !    = -(gas tendency) (mol/mol/s)
   real(r8), intent(out) :: sum_uprt         ! gas uptake rate summed over modes/bins (1/s)

   integer  :: n
   real(r8) :: avg_uprt                      ! average uptake rate over deltat (1/s)
   real(r8) :: deltatxx
   real(r8) :: fgain(nbins)                  ! fraction of total gas uptake going to mode/bin n

! use this for tendency calcs to avoid generating very small negative values
   deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)

   sum_uprt = 0.0_r8
   do n = 1, nbins
      if (mask(n)) then
         fgain(n) = uptkrate(n)*diffus_ratio
         sum_uprt = sum_uprt + fgain(n)
      else
         fgain(n) = 0.0_r8
      end if
   end do

!  at this point (sum_uprt <= 0.0) only when all the fgain are zero
   if (sum_uprt > 0.0_r8) then
      do n = 1, nbins
         fgain(n) = fgain(n) / sum_uprt
      end do
   end if

!   uptake amount (fraction of gas uptaken) over deltat
   avg_uprt = (1.0_r8 - exp(-deltatxx*sum_uprt))/deltatxx

   sum_dqdt = qgas * avg_uprt

   do n = 1, nbins
      dqdt_aer(n) = fgain(n)*sum_dqdt
   end do

end subroutine nonvolatile_uptake_cell

!----------------------------------------------------------------------
subroutine h2so4_reversible_uptake_cell( nbins, deltat, uptkrate, mask, qgas, qaer, qgas_eq, &
                                         dqdt_aer, sum_dqdt )

!
! reversible uptake of h2so4 gas by the aerosol modes or bins in one grid cell
! (used above the tropopause, where so4 aerosol can evaporate)
!
! the h2so4 gas relaxes towards the uptake-rate weighted equilibrium value over deltat,
! and each mode/bin takes up (or loses) so4 according to the difference between the
! average h2so4 gas and its equilibrium h2so4 gas
!

   integer,  intent(in)  :: nbins            ! number of aerosol modes or bins
   real(r8), intent(in)  :: deltat           ! time step (s)
   real(r8), intent(in)  :: uptkrate(nbins)  ! h2so4 gas-to-aerosol mass transfer rate (1/s)
   logical,  intent(in)  :: mask(nbins)      ! true if mode/bin takes up h2so4
   real(r8), intent(in)  :: qgas             ! h2so4 gas mixing ratio (mol/mol)
   real(r8), intent(in)  :: qaer(nbins)      ! so4 aerosol mixing ratio in each mode/bin (mol/mol)
   real(r8), intent(in)  :: qgas_eq(nbins)   ! equilibrium h2so4 gas mixing ratio over each mode/bin (mol/mol)
   real(r8), intent(out) :: dqdt_aer(nbins)  ! so4 aerosol tendency in each mode/bin (mol/mol/s)
   real(r8), intent(out) :: sum_dqdt         ! so4 aerosol tendency summed over modes/bins
                                             !    = -(h2so4 gas tendency) (mol/mol/s)

   integer  :: n
   real(r8) :: deltatxx
   real(r8) :: sum_uprt                      ! h2so4 uptake rate summed over modes/bins (1/s)
   real(r8) :: tmp_kxt, tmp_pxt
   real(r8) :: tmp_so4a_bgn, tmp_so4a_end
   real(r8) :: tmp_so4g_avg, tmp_so4g_bgn, tmp_so4g_equ

! use this for tendency calcs to avoid generating very small negative values
   deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)

   sum_uprt = 0.0_r8
   do n = 1, nbins
      if (mask(n)) sum_uprt = sum_uprt + uptkrate(n)
   end do

   tmp_kxt = deltatxx*sum_uprt  ! sum over modes of uptake_rate*deltat
   tmp_pxt = 0.0_r8
   do n = 1, nbins
      if (.not. mask(n)) cycle
      tmp_pxt = tmp_pxt + uptkrate(n)*qgas_eq(n)
   end do
   tmp_pxt = max( 0.0_r8, tmp_pxt*deltatxx )  ! sum over modes of uptake_rate*qgas_eq*deltat
   tmp_so4g_bgn = qgas
   ! calc avg h2so4(g) over deltat
   if (tmp_kxt >= 1.0e-5_r8) then
      ! exponential decay towards equilibrium value solution
      tmp_so4g_equ = tmp_pxt/tmp_kxt
      tmp_so4g_avg = tmp_so4g_equ + (tmp_so4g_bgn-tmp_so4g_equ)*(1.0_r8-exp(-tmp_kxt))/tmp_kxt
   else
      ! first order approx for tmp_kxt small
      tmp_so4g_avg = tmp_so4g_bgn*(1.0_r8-0.5_r8*tmp_kxt) + 0.5_r8*tmp_pxt
   end if

   sum_dqdt = 0.0_r8
   do n = 1, nbins
      if (.not. mask(n)) then
         dqdt_aer(n) = 0.0_r8
         cycle
      end if
      ! calc change to so4(a) in mode n
      tmp_so4a_bgn = qaer(n)
      tmp_so4a_end = tmp_so4a_bgn + deltatxx*uptkrate(n)*(tmp_so4g_avg-qgas_eq(n))
      tmp_so4a_end = max( 0.0_r8, tmp_so4a_end )
      dqdt_aer(n) = (tmp_so4a_end - tmp_so4a_bgn)/deltatxx
      sum_dqdt = sum_dqdt + dqdt_aer(n)
   end do

end subroutine h2so4_reversible_uptake_cell

!----------------------------------------------------------------------
pure elemental function fuchs_sutugin(knudsen) result(fs)

!
! Fuchs-Sutugin correction factor F(Kn,ac) for accomodation coefficient ac = 1
!
! Changed by Manish Shrivastava on 7/17/2013 to use accom=1; because we do not know better
!   following assumes accomodation coefficient = ac = 1. instead 0.65 ! answer change needs to be tested
!   (Adams & Seinfeld, 2002, JGR, and references therein)
!           fuchs_sutugin = (0.75*ac*(1. + knudsen)) /
!                           (knudsen*(1.0 + knudsen + 0.283*ac) + 0.75*ac)
!

   real(r8), intent(in) :: knudsen ! Knudsen number
   real(r8) :: fs

   fs = (0.4875_r8*(1._r8 + knudsen)) /   &
        (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)

end function fuchs_sutugin

end module gas_aer_uptkrate_mod
