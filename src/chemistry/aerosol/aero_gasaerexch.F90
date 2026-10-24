!----------------------------------------------------------------------
! Gas-aerosol exchange (gas condensation onto aerosol particles) common
! to the aerosol models (MAM, CARMA). The aerosol model is accessed only
! through the abstract aerosol_properties interface; model-specific
! bookkeeping (species indices, renaming, aging, VBS lumping, ...) stays
! with the calling aerosol model.
!----------------------------------------------------------------------
module aero_gasaerexch

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   private

   public :: aero_gasaerexch_uptkrates
   public :: aero_gasaerexch_soa

   ! gas diffusivity ratios relative to h2so4
   ! (differences in fuchs-sutugin and accom coef ignored)
   real(r8), parameter, public :: soa_diffus_ratio = 0.81_r8 ! soa/h2so4
   real(r8), parameter, public :: nh3_diffus_ratio = 2.08_r8 ! nh3/h2so4

   ! methods for soa gas uptake
   integer, parameter, public :: method_soa_none = 0         ! no uptake
   integer, parameter, public :: method_soa_irreversible = 1 ! irreversible uptake done like h2so4 uptake
   integer, parameter, public :: method_soa_reversible = 2   ! reversible uptake using aero_soaexch

contains

   !----------------------------------------------------------------------
   ! computes H2SO4 gas-to-aerosol mass transfer rates (1/s) for all
   ! aerosol modes/bins
   !----------------------------------------------------------------------
   subroutine aero_gasaerexch_uptkrates( aero_props, ncol, nlev, top_lev, numdens, wetdiam, &
                                         temp, pmid, uptkrate )

      use aerosol_properties_mod, only: aerosol_properties

      class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
      integer,  intent(in)  :: ncol              ! number of columns
      integer,  intent(in)  :: nlev              ! number of vertical levels
      integer,  intent(in)  :: top_lev           ! top level index for computation
      real(r8), intent(in)  :: numdens(:,:,:)    ! number density of each mode/bin (#/m3)
      real(r8), intent(in)  :: wetdiam(:,:,:)    ! wet diameter of each mode/bin (m)
      real(r8), intent(in)  :: temp(:,:)         ! temperature (K)
      real(r8), intent(in)  :: pmid(:,:)         ! air pressure (Pa)
      real(r8), intent(out) :: uptkrate(:,:,:)   ! gas-to-aerosol mass transfer rates (1/s)

      integer :: m

      do m = 1, aero_props%nbins()
         call aero_props%gas_uptkrate( m, ncol, nlev, top_lev, numdens(:,:,m), wetdiam(:,:,m), &
                                       temp, pmid, uptkrate(:,:,m) )
      end do

   end subroutine aero_gasaerexch_uptkrates

   !----------------------------------------------------------------------
   ! computes soa gas uptake tendencies in a single grid cell
   !
   !   method_soa_reversible   -- reversible uptake using aero_soaexch
   !   method_soa_irreversible -- irreversible uptake done like h2so4 uptake
   !   method_soa_none         -- no uptake
   !----------------------------------------------------------------------
   subroutine aero_gasaerexch_soa( method_soa, deltat, temp, pres, &
                                   nbins, nsoa, npoa, do_bin, do_soa, do_soag, uptkrate, &
                                   mw_poa_host, mw_soa_host, &
                                   qold_soag, qold_soa, qold_poa, &
                                   dqdt_soa, sum_dqdt_soa )

      use cam_abortutils,   only: endrun
      use aero_soaexch_mod, only: aero_soaexch

      integer,  intent(in)  :: method_soa                ! soa uptake method
      real(r8), intent(in)  :: deltat                    ! time step (s)
      real(r8), intent(in)  :: temp                      ! air temperature (K)
      real(r8), intent(in)  :: pres                      ! air pressure (Pa)
      integer,  intent(in)  :: nbins                     ! number of modes/bins
      integer,  intent(in)  :: nsoa                      ! number of soa species
      integer,  intent(in)  :: npoa                      ! number of poa species
      logical,  intent(in)  :: do_bin(nbins)             ! mode/bin takes part in reversible soa exchange
      logical,  intent(in)  :: do_soa(nbins,nsoa)        ! soa species takes up gas in mode/bin
      logical,  intent(in)  :: do_soag(nsoa)             ! soa gas species is present
      real(r8), intent(in)  :: uptkrate(nbins)           ! h2so4 gas-to-aerosol mass transfer rates (1/s)
      real(r8), intent(in)  :: mw_poa_host(npoa)         ! molec wght of poa used in host code
      real(r8), intent(in)  :: mw_soa_host(nsoa)         ! molec wght of soa used in host code
      real(r8), intent(in)  :: qold_soag(nsoa)           ! soa gas mixing ratios (mol/mol)
      real(r8), intent(in)  :: qold_soa(nbins,nsoa)      ! soa aerosol mixing ratios (mol/mol)
      real(r8), intent(in)  :: qold_poa(nbins,npoa)      ! poa aerosol mixing ratios (mol/mol)
      real(r8), intent(out) :: dqdt_soa(nbins,nsoa)      ! soa aerosol tendencies (mol/mol/s)
      real(r8), intent(out) :: sum_dqdt_soa(nsoa)        ! total soa gas uptake (mol/mol/s) -- negative of gas tendency

      integer, parameter :: niter_max = 1000

      integer  :: jsoa, n, niter
      real(r8) :: avg_uprt_soa(nsoa)
      real(r8) :: deltatxx
      real(r8) :: dqdt_soag(nsoa)
      real(r8) :: fgain_soa(nbins,nsoa)
      real(r8) :: sum_uprt_soa(nsoa)
      real(r8) :: uptkrate_soa(nbins,nsoa)

      ! use this for tendency calcs to avoid generating very small negative values
      deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)

      !   fgain_soa(n,jsoa) = fraction of total soa gas uptake going to mode/bin n
      sum_uprt_soa(:) = 0.0_r8
      do n = 1, nbins
         do jsoa = 1, nsoa
            if (do_soa(n,jsoa)) then
               fgain_soa(n,jsoa) = uptkrate(n)*soa_diffus_ratio
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
         if ( do_soag(jsoa) ) then
            sum_dqdt_soa(jsoa) = qold_soag(jsoa) * avg_uprt_soa(jsoa)
         else
            sum_dqdt_soa(jsoa) = 0.0_r8
         end if
      end do

      select case (method_soa)
      case (method_soa_reversible)
         !   compute tendencies for soa gas and soa interstitial aerosol
         !   using soa parameterization
         call aero_soaexch( deltat, temp, pres, &
              niter, niter_max, nbins, do_bin, npoa, nsoa, &
              mw_poa_host, mw_soa_host, &
              qold_soag, qold_soa, qold_poa, uptkrate_soa, &
              dqdt_soag, dqdt_soa )
         sum_dqdt_soa(:) = -dqdt_soag(:)

      case (method_soa_irreversible)
         !   compute tendencies for soa interstitial aerosol
         !   due to simple gas uptake
         do jsoa = 1, nsoa
            do n = 1, nbins
               dqdt_soa(n,jsoa) = fgain_soa(n,jsoa)*sum_dqdt_soa(jsoa)
            end do
         end do

      case (method_soa_none)
         dqdt_soa(:,:) = 0.0_r8
         sum_dqdt_soa(:) = 0.0_r8

      case default
         call endrun('aero_gasaerexch_soa: invalid method_soa')
      end select

   end subroutine aero_gasaerexch_soa

end module aero_gasaerexch
