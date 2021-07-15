!------------------------------------------------------------------------------
! Utility routines for medium energy electron ionization base on Ap
! geomagnetic activity index
!------------------------------------------------------------------------------
module mee_ap_util_mod
  use shr_kind_mod,   only : r8 => shr_kind_r8

  implicit none

contains

  !------------------------------------------------------------------------------
  ! Calculating an energy spectrum for a specific L-shell and Ap calculation
  !
  ! A routine to return an energy spectrum for a specific L-shell and Ap
  ! Calcs based on:
  !
  ! van de Kamp, M., Seppala, A., Clilverd, M. A., Rodger, C. J., Verronen, P. T., and Whittaker, I. C. (2016),
  ! A model providing long‐term datasets of energetic electron precipitation during geomagnetic storms,
  ! J. Geophys. Res. Atmos., 121, 12,520– 12,540,
  ! [doi:10.1002/2015JD024212](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JD024212)
  !------------------------------------------------------------------------------
  function FluxSpectrum( energies, lshell, Ap ) result(flux)

    real(r8), intent(in) :: energies(:)
    real(r8), intent(in) :: lshell
    real(r8), intent(in) :: Ap

    real(r8) :: flux(size(energies))

    real(r8) :: lpp
    real(r8) :: Spp
    real(r8) :: A
    real(r8) :: b
    real(r8) :: c
    real(r8) :: s
    real(r8) :: d
    real(r8) :: F30
    real(r8) :: E
    real(r8) :: bk
    real(r8) :: sk
    real(r8) :: k
    real(r8) :: x

    lpp = -0.7430*log(Ap) + 6.5257
    Spp = lshell - lpp

    ! vdK2016 eqn.(8)

    A = 8.2091*Ap**(0.16255)
    b = 1.3754*Ap**(0.33042)
    c = 0.13334*Ap**(0.42616)
    s = 2.2833*Ap**(-0.22990)
    d = 2.7563e-4*Ap**(2.6116)

    F30 = exp(A) / (exp(-b*(Spp-s)) + exp(c*(Spp-s)) + d)

    ! vdK2016 eqn.(9)

    E = 3.3777*Ap**(-1.7038) + 0.15
    bk = 3.7632*Ap**(-0.16034)
    sk = 12.184*Ap**(-0.30111)

    k = -1.0 / (E*exp(-bk*Spp) + 0.30450*cosh(0.20098*(Spp-sk))) - 1._r8

    x=k+1
    c = F30*x/(1e3**x-30.**x)
    flux(:) = energies(:)**k*c

  end function FluxSpectrum

  !------------------------------------------------------------------------------
  ! Calculate how energy from top of atmosphere is deposited in rest of atmosphere
  !
  ! The function takes the energy spectrum at the top of the atmosphere and
  ! calculates how that energy is deposited in the atmosphere using the parameterization
  ! described in [Fang et al., (2010)](https://opensky.ucar.edu/islandora/object/articles:10653)
  !
  ! Fang, X., C. E. Randall, D. Lummerzheim, W. Wang, G. Lu, S. C. Solomon,
  ! and R. A. Frahm (2010), Parameterization of monoenergetic electron impact
  ! ionization, Geophys. Res. Lett., 37, L22106, [doi:10.1029/2010GL045406.]
  ! (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2010GL045406)
  !
  ! Application of the new parameterization requires the following steps:
  !
  ! 1. Calculate the Ci coefficients using equation (5) and Table 1.
  ! 2. Calculate the y values throughout the atmosphere using equation (1).
  ! 3. Calculate the normalized energy dissipation f values using equation (4).
  ! 4. Obtain the altitude profile of qtot by substituting the f values into equation (3).
  !
  !------------------------------------------------------------------------------
  function iprmono(e, flux, rho, scaleh) result(ipr)
    real(r8), intent(in) :: e(:)
    real(r8), intent(in) :: flux(:)
    real(r8), intent(in) :: rho(:)
    real(r8), intent(in) :: scaleh(:)

    real(r8) :: ipr(size(e),size(rho))
    integer :: nenergies, n
    integer :: nlayers, k

    real(r8) :: y(size(rho))
    real(r8) :: f(size(rho))
    real(r8) :: Qmono

    ! assign constants
    real(r8), parameter :: epsilon = 0.035 ! keV energy loss per ion pair produced

    ipr = 0._r8
    nenergies = size(e)
    nlayers = size(rho)

    do n = 1,nenergies

       ! step 1. (eq. 1)
       y(:) = (2/e(n))*(rho(:)*scaleh(:)/6.0e-6_r8)**0.7_r8

       do k = 1,nlayers
          f(k) = fang(y(k), e(n))
       end do

       ! calculate ipr (qtot) using eq. (3) for a specified flux at ea. energy
       Qmono = flux(n)*e(n) ! (keV cm−2 s−1)
       ipr(n,:) = f(:)*Qmono/(epsilon*scaleh(:))
    end do

  contains

    function fang(y, Emono) result(f)
      real(r8), intent(in) :: y, Emono
      real(r8) :: f
      ! Input:
      ! y - normalized atmospheric column mass as a function of vertical location (z)
      ! Emono - is incident electron energy (keV)
      ! Output:
      ! f - quanity calculated by eqn. (4)

      real(r8) :: lne, lne2, lne3
      real(r8) :: c(8)
      integer :: i

      ! Table 1.
      real(r8), parameter :: p1(8,4) = reshape( &
           (/ 1.24616E+0,  1.45903E+0, -2.42269E-1,  5.95459E-2, &
              2.23976E+0, -4.22918E-7,  1.36458E-2,  2.53332E-3, &
              1.41754E+0,  1.44597E-1,  1.70433E-2,  6.39717E-4, &
              2.48775E-1, -1.50890E-1,  6.30894E-9,  1.23707E-3, &
             -4.65119E-1, -1.05081E-1, -8.95701E-2,  1.22450E-2, &
              3.86019E-1,  1.75430E-3, -7.42960E-4,  4.60881E-4, &
             -6.45454E-1,  8.49555E-4, -4.28581E-2, -2.99302E-3, &
              9.48930E-1,  1.97385E-1, -2.50660E-3, -2.06938E-3 /), shape=(/8,4/),order=(/2,1/))


      ! terms in eq. (5)
      lne = log(Emono)
      lne2 = lne*lne
      lne3 = lne*lne2

      ! step 2. calculate the C array in (5)
      do i = 1,8
         c(i) = exp(p1(i,1) + p1(i,2)*lne + p1(i,3)*lne2 + p1(i,4)*lne3)
      end do

      ! eq. (4) - Normalized energy deposition
      f = c(1)*y**c(2)*exp(-c(3)*y**c(4)) + c(5)*y**c(6)*exp(-c(7)*y**c(8))

    end function fang

  end function iprmono


  !------------------------------------------------------------------------------
  ! Generate a grid of energy bins for the flux spectrum.
  ! The energy range of the spectrum is 30–1000 keV,
  ! with nbins of logarithmically spaced grid points.
  !------------------------------------------------------------------------------
  subroutine gen_energy_grid(nbins, energies, deltas)
    integer, intent(in) :: nbins
    real(r8),intent(out) :: energies(nbins)
    real(r8),intent(out) :: deltas(nbins)

    integer :: i
    real(r8) :: e1,e2
    real(r8) :: low,med,hig

    e1 = log10(30._r8)
    e2 = log10(1000._r8)

    do i = 1,nbins
       low = e1 + (e2-e1)*(i-1.0_r8)/nbins
       med = e1 + (e2-e1)*(i-0.5_r8)/nbins
       hig = e1 + (e2-e1)*(i)/nbins

       energies(i) = 10**med
       deltas(i) = (10**hig)-(10**low)
    end do

  end subroutine gen_energy_grid

  !------------------------------------------------------------------------------
  ! returns L-Shell number for a given magnetic latitude (radians)
  !------------------------------------------------------------------------------
  pure function maglat2lshell( rmaglat ) result(lshell)
    real(r8), intent(in) :: rmaglat ! mag latitude in radians
    real(r8) :: lshell  ! magnetosphere L-Shell number

    lshell = 1.01_r8/(cos(rmaglat)**2)

  end function maglat2lshell

end module mee_ap_util_mod
