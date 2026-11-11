!----------------------------------------------------------------------
! Gas-aerosol exchange (gas condensation onto aerosol particles) common
! to the aerosol models (MAM, CARMA). The aerosol model is accessed only
! through the abstract aerosol_properties interface; model-specific
! bookkeeping (constituent indices, renaming, aging, ...) stays with the
! calling aerosol model.
!----------------------------------------------------------------------
module aero_gasaerexch

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   private

   public :: aero_gasaerexch_uptkrates
   public :: aero_gasaerexch_soa
   public :: aero_gasaerexch_soa_column
   public :: aero_gasaerexch_partition
   public :: aero_gasaerexch_reversible
   public :: aero_gasaerexch_nh4_limit
   public :: aero_gasaerexch_uptk_diag_init
   public :: aero_gasaerexch_uptk_diag_out
   public :: aero_gasaerexch_spec_ndx
   public :: aero_gasaerexch_soag_ndx

   ! gas diffusivity ratios relative to h2so4
   ! (differences in fuchs-sutugin and accom coef ignored)
   real(r8), parameter, public :: h2so4_diffus_ratio = 1.0_r8 ! h2so4/h2so4
   real(r8), parameter, public :: soa_diffus_ratio = 0.81_r8 ! soa/h2so4
   real(r8), parameter, public :: nh3_diffus_ratio = 2.08_r8 ! nh3/h2so4

   ! methods for soa gas uptake
   integer, parameter, public :: method_soa_none = 0         ! no uptake
   integer, parameter, public :: method_soa_irreversible = 1 ! irreversible uptake done like h2so4 uptake
   integer, parameter, public :: method_soa_reversible = 2   ! reversible uptake using aero_soaexch

   !----------------------------------------------------------------------
   ! SOA condensation/evaporation history diagnostics (kg/kg/s)
   !
   ! Condensation (positive) and evaporation (negative) SOA tendencies,
   ! optionally output as totals over all SOA species ('qcon_gaex',
   ! 'qevap_gaex'), per SOA species ('qcon_gaexNN', 'qevap_gaexNN'), and/or
   ! per group of SOA species ('qcon<sfx>_gaex', 'qevap<sfx>_gaex').
   ! The configuration is set once at initialization; the sums for a call
   ! are kept in a soa_condevap_acc_t local to the caller (thread safe).
   !----------------------------------------------------------------------
   type, public :: soa_condevap_diag_t
      private
      integer :: nsoa = 0                        ! number of soa species
      logical :: do_total = .false.              ! output totals over all soa species
      logical :: do_spec = .false.               ! output per soa species
      integer :: ngrp = 0                        ! number of soa species groups output
      integer, allocatable :: grp(:)             ! group index of each soa species (0 = none)
      character(len=8), allocatable :: grp_sfx(:) ! field name suffix of each group
   contains
      procedure :: init => soa_condevap_init
      procedure :: reset => soa_condevap_reset
      procedure :: accum => soa_condevap_accum
      procedure :: accum_column => soa_condevap_accum_column
      procedure :: outfld => soa_condevap_outfld
   end type soa_condevap_diag_t

   ! sums of the SOA condensation/evaporation diagnostics for one chunk
   type, public :: soa_condevap_acc_t
      real(r8), allocatable :: con_tot(:,:)      ! (pcols,pver)
      real(r8), allocatable :: evap_tot(:,:)     ! (pcols,pver)
      real(r8), allocatable :: con_spc(:,:,:)    ! (pcols,pver,nsoa)
      real(r8), allocatable :: evap_spc(:,:,:)   ! (pcols,pver,nsoa)
      real(r8), allocatable :: con_grp(:,:,:)    ! (pcols,pver,ngrp)
      real(r8), allocatable :: evap_grp(:,:,:)   ! (pcols,pver,ngrp)
   end type soa_condevap_acc_t

   character(len=*), parameter :: condevap_units = 'kg/kg/s'

   !----------------------------------------------------------------------
   ! column tendency history diagnostics ('<species><suffix>', e.g.,
   ! 'so4_a1_sfgaex1') of species involved in gas-aerosol exchange
   ! (kg/m2/s, or #/m2/s for number)
   !
   ! The column tendencies are accumulated by the caller in mol/mol/s * kg/m2
   ! (#/kmol/s * kg/m2 for number) and are converted to the output units
   ! with a caller-supplied factor (molecular weight / mwdry) when output.
   !----------------------------------------------------------------------
   type, public :: col_tend_diag_t
      private
      integer :: nfld = 0                              ! number of fields
      character(len=64), allocatable :: fldname(:)     ! field names
   contains
      procedure :: init => col_tend_init
      procedure :: outfld => col_tend_outfld
   end type col_tend_diag_t

   !----------------------------------------------------------------------
   ! SOA species configuration of an aerosol model for the SOA gas-aerosol
   ! exchange, set from the aerosol properties object at initialization
   !
   ! SOA species j of a mode/bin takes up SOA gas (SOAG) species j.
   ! All SOA (POA) species have the same molecular weight.
   !----------------------------------------------------------------------
   type, public :: soa_exch_t
      integer :: nbins = 0                       ! number of modes/bins
      integer :: nsoa = 0                        ! number of soa species (maximum over the modes/bins)
      integer :: npoa = 0                        ! number of poa species (maximum over the modes/bins)
      integer :: nsoag = 0                       ! number of soa gas constituents
      integer, allocatable :: nsoa_bin(:)        ! number of soa species in each mode/bin
      integer, allocatable :: npoa_bin(:)        ! number of poa species in each mode/bin
      integer, allocatable :: soa_spec(:,:)      ! species index within the mode/bin of each soa species
                                                 ! (nbins,nsoa; <= 0 if absent)
      integer, allocatable :: poa_spec(:,:)      ! species index within the mode/bin of each poa species
                                                 ! (nbins,npoa; <= 0 if absent)
      integer, allocatable :: soag_cnst(:)       ! constituent index of the soa gas taken up by each soa
                                                 ! species (nsoa; <= 0 if absent)
      real(r8) :: mw_soa = -huge(1.0_r8)         ! molecular weight of the soa species (kg/kmol)
      real(r8) :: mw_poa = -huge(1.0_r8)         ! molecular weight of the poa species (kg/kmol)
      real(r8), allocatable :: mw_soa_host(:)    ! molec wght of soa used in host code (nsoa)
      real(r8), allocatable :: mw_poa_host(:)    ! molec wght of poa used in host code (npoa)
   contains
      procedure :: init => soa_exch_init
      procedure :: tend => soa_exch_tend
   end type soa_exch_t

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

   !----------------------------------------------------------------------
   ! computes soa gas uptake tendencies for all grid cells of a chunk
   ! (levels top_lev to nlev), using aero_gasaerexch_soa in each cell
   !
   ! In each cell, the soa aerosol mixing ratio of mode/bin n and soa species
   ! jsoa is soa(i,k,n,jsoa) where do_soa(n,jsoa) is true, and zero otherwise.
   !----------------------------------------------------------------------
   subroutine aero_gasaerexch_soa_column( method_soa, ncol, nlev, top_lev, deltat, temp, pmid, &
                                          nbins, nsoa, npoa, do_bin, do_soa, do_soag, uptkrate, &
                                          mw_poa_host, mw_soa_host, soag, soa, poa, &
                                          dqdt_soa, sum_dqdt_soa )

      integer,  intent(in)  :: method_soa                ! soa uptake method
      integer,  intent(in)  :: ncol                      ! number of columns
      integer,  intent(in)  :: nlev                      ! number of vertical levels
      integer,  intent(in)  :: top_lev                   ! top level index for computation
      real(r8), intent(in)  :: deltat                    ! time step (s)
      real(r8), intent(in)  :: temp(:,:)                 ! air temperature (K)
      real(r8), intent(in)  :: pmid(:,:)                 ! air pressure (Pa)
      integer,  intent(in)  :: nbins                     ! number of modes/bins
      integer,  intent(in)  :: nsoa                      ! number of soa species
      integer,  intent(in)  :: npoa                      ! number of poa species
      logical,  intent(in)  :: do_bin(nbins)             ! mode/bin takes part in reversible soa exchange
      logical,  intent(in)  :: do_soa(nbins,nsoa)        ! soa species takes up gas in mode/bin
      logical,  intent(in)  :: do_soag(nsoa)             ! soa gas species is present
      real(r8), intent(in)  :: uptkrate(:,:,:)           ! h2so4 gas-to-aerosol mass transfer rates (1/s) (i,k,n)
      real(r8), intent(in)  :: mw_poa_host(npoa)         ! molec wght of poa used in host code
      real(r8), intent(in)  :: mw_soa_host(nsoa)         ! molec wght of soa used in host code
      real(r8), intent(in)  :: soag(:,:,:)               ! soa gas mixing ratios (mol/mol) (i,k,jsoa)
      real(r8), intent(in)  :: soa(:,:,:,:)              ! soa aerosol mixing ratios (mol/mol) (i,k,n,jsoa)
      real(r8), intent(in)  :: poa(:,:,:,:)              ! poa aerosol mixing ratios (mol/mol) (i,k,n,j)
      real(r8), intent(out) :: dqdt_soa(:,:,:,:)         ! soa aerosol tendencies (mol/mol/s) (i,k,n,jsoa)
      real(r8), intent(out) :: sum_dqdt_soa(:,:,:)       ! total soa gas uptake (mol/mol/s) (i,k,jsoa)
                                                         ! -- negative of gas tendency

      integer  :: i, k, n, j
      real(r8) :: uptkratebb(nbins)
      real(r8) :: qold_soag(nsoa)
      real(r8) :: qold_soa(nbins,nsoa)
      real(r8) :: qold_poa(nbins,npoa)
      real(r8) :: dqdt_soa_cell(nbins,nsoa)
      real(r8) :: sum_dqdt_soa_cell(nsoa)

      dqdt_soa(:,:,:,:) = 0.0_r8
      sum_dqdt_soa(:,:,:) = 0.0_r8

      do k = top_lev, nlev
         do i = 1, ncol
            do n = 1, nbins
               uptkratebb(n) = uptkrate(i,k,n)
               do j = 1, nsoa
                  if (do_soa(n,j)) then
                     qold_soa(n,j) = soa(i,k,n,j)
                  else
                     qold_soa(n,j) = 0.0_r8
                  end if
               end do
               do j = 1, npoa
                  qold_poa(n,j) = poa(i,k,n,j)
               end do
            end do
            do j = 1, nsoa
               qold_soag(j) = soag(i,k,j)
            end do

            call aero_gasaerexch_soa( method_soa, deltat, temp(i,k), pmid(i,k), &
                                      nbins, nsoa, npoa, do_bin, do_soa, do_soag, uptkratebb, &
                                      mw_poa_host, mw_soa_host, &
                                      qold_soag, qold_soa, qold_poa, &
                                      dqdt_soa_cell, sum_dqdt_soa_cell )

            dqdt_soa(i,k,:,:) = dqdt_soa_cell(:,:)
            sum_dqdt_soa(i,k,:) = sum_dqdt_soa_cell(:)
         end do
      end do

   end subroutine aero_gasaerexch_soa_column

   !----------------------------------------------------------------------
   ! sets the SOA species configuration from the aerosol properties object:
   ! soa ('s-organic') and poa ('p-organic') species of each mode/bin, the
   ! soa gas (SOAG) constituents and the soa and poa molecular weights
   ! (endrun if the soa or poa species have different molecular weights)
   !----------------------------------------------------------------------
   subroutine soa_exch_init( self, aero_props )

      use aerosol_properties_mod, only: aerosol_properties
      use cam_abortutils, only: endrun

      class(soa_exch_t), intent(out) :: self
      class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object

      integer :: m, j
      integer, allocatable :: soag_ndx(:)
      real(r8) :: spec_mw
      character(len=*), parameter :: subname = 'aero_gasaerexch::soa_exch_init'

      self%nbins = aero_props%nbins()
      allocate( self%nsoa_bin(self%nbins), self%npoa_bin(self%nbins) )

      call aero_gasaerexch_spec_ndx( aero_props, 's-organic', self%nsoa_bin, self%soa_spec )
      call aero_gasaerexch_spec_ndx( aero_props, 'p-organic', self%npoa_bin, self%poa_spec )
      self%nsoa = maxval(self%nsoa_bin)
      self%npoa = maxval(self%npoa_bin)

      ! soa gas species j is taken up by soa species j
      call aero_gasaerexch_soag_ndx( soag_ndx )
      self%nsoag = size(soag_ndx)
      allocate( self%soag_cnst(self%nsoa) )
      self%soag_cnst(:) = -1
      do j = 1, min(self%nsoa, self%nsoag)
         self%soag_cnst(j) = soag_ndx(j)
      end do

      ! molecular weights
      do m = 1, self%nbins
         do j = 1, self%nsoa_bin(m)
            call aero_props%get( m, self%soa_spec(m,j), spec_mw=spec_mw )
            if (self%mw_soa > 0.0_r8) then
               if (spec_mw /= self%mw_soa) then
                  call endrun( subname//': s-organic species molecular weights differ' )
               end if
            end if
            self%mw_soa = spec_mw
         end do
         do j = 1, self%npoa_bin(m)
            call aero_props%get( m, self%poa_spec(m,j), spec_mw=spec_mw )
            if (self%mw_poa > 0.0_r8) then
               if (spec_mw /= self%mw_poa) then
                  call endrun( subname//': p-organic species molecular weights differ' )
               end if
            end if
            self%mw_poa = spec_mw
         end do
      end do

      allocate( self%mw_soa_host(self%nsoa), self%mw_poa_host(self%npoa) )
      self%mw_soa_host(:) = self%mw_soa
      self%mw_poa_host(:) = self%mw_poa

   end subroutine soa_exch_init

   !----------------------------------------------------------------------
   ! computes soa gas uptake tendencies for all grid cells of a chunk
   ! (levels top_lev to nlev) from the soa gas and soa/poa aerosol mixing
   ! ratios (mol/mol) in caller arrays:
   !   soa gas species j:            qgas(i,k,lgas_soa(j))       if do_soag(j)
   !   soa species j of mode/bin n:  qaer(i,k,laer_soa(n,j))     if laer_soa(n,j) > 0
   !   poa species j of mode/bin n:  qaer(i,k,laer_poa(n,j))     if laer_poa(n,j) > 0
   ! (zero otherwise). The tendencies are zero if no soa gas is present.
   !----------------------------------------------------------------------
   subroutine soa_exch_tend( self, method_soa, ncol, nlev, top_lev, deltat, temp, pmid, &
                             do_bin, do_soa, do_soag, uptkrate, &
                             qgas, lgas_soa, qaer, laer_soa, laer_poa, &
                             dqdt_soa, sum_dqdt_soa )

      class(soa_exch_t), intent(in) :: self
      integer,  intent(in)  :: method_soa                     ! soa uptake method
      integer,  intent(in)  :: ncol                           ! number of columns
      integer,  intent(in)  :: nlev                           ! number of vertical levels
      integer,  intent(in)  :: top_lev                        ! top level index for computation
      real(r8), intent(in)  :: deltat                         ! time step (s)
      real(r8), intent(in)  :: temp(:,:)                      ! air temperature (K)
      real(r8), intent(in)  :: pmid(:,:)                      ! air pressure (Pa)
      logical,  intent(in)  :: do_bin(self%nbins)             ! mode/bin takes part in reversible soa exchange
      logical,  intent(in)  :: do_soa(self%nbins,self%nsoa)   ! soa species takes up gas in mode/bin
      logical,  intent(in)  :: do_soag(self%nsoa)             ! soa gas species is present
      real(r8), intent(in)  :: uptkrate(:,:,:)                ! h2so4 gas-to-aerosol mass transfer rates (1/s)
      real(r8), intent(in)  :: qgas(:,:,:)                    ! gas mixing ratios (mol/mol)
      integer,  intent(in)  :: lgas_soa(self%nsoa)            ! qgas index of each soa gas species
      real(r8), intent(in)  :: qaer(:,:,:)                    ! aerosol mixing ratios (mol/mol)
      integer,  intent(in)  :: laer_soa(self%nbins,self%nsoa) ! qaer index of each soa species (<= 0 if none)
      integer,  intent(in)  :: laer_poa(self%nbins,self%npoa) ! qaer index of each poa species (<= 0 if none)
      real(r8), intent(out) :: dqdt_soa(:,:,:,:)              ! soa aerosol tendencies (mol/mol/s) (i,k,n,jsoa)
      real(r8), intent(out) :: sum_dqdt_soa(:,:,:)            ! total soa gas uptake (mol/mol/s) (i,k,jsoa)
                                                              ! -- negative of gas tendency

      integer  :: i, k, n, j, l
      real(r8) :: soag(ncol,nlev,self%nsoa)
      real(r8) :: soa(ncol,nlev,self%nbins,self%nsoa)
      real(r8) :: poa(ncol,nlev,self%nbins,self%npoa)

      if (.not. any(do_soag)) then ! no soa gas, no uptake
         dqdt_soa(:,:,:,:) = 0.0_r8
         sum_dqdt_soa(:,:,:) = 0.0_r8
         return
      end if

      soag(:,:,:) = 0.0_r8
      soa(:,:,:,:) = 0.0_r8
      poa(:,:,:,:) = 0.0_r8
      do k = top_lev, nlev
         do i = 1, ncol
            do j = 1, self%nsoa
               if (do_soag(j)) soag(i,k,j) = qgas(i,k,lgas_soa(j))
            end do
            do n = 1, self%nbins
               do j = 1, self%nsoa
                  l = laer_soa(n,j)
                  if (l > 0) soa(i,k,n,j) = qaer(i,k,l)
               end do
               do j = 1, self%npoa
                  l = laer_poa(n,j)
                  if (l > 0) poa(i,k,n,j) = qaer(i,k,l)
               end do
            end do
         end do
      end do

      call aero_gasaerexch_soa_column( method_soa, ncol, nlev, top_lev, deltat, temp, pmid, &
                                       self%nbins, self%nsoa, self%npoa, do_bin, do_soa, do_soag, uptkrate, &
                                       self%mw_poa_host, self%mw_soa_host, soag, soa, poa, &
                                       dqdt_soa, sum_dqdt_soa )

   end subroutine soa_exch_tend

   !----------------------------------------------------------------------
   ! partitions the uptake of a non-volatile gas among modes/bins in a grid
   ! cell:
   !   fgain(n) = fraction of the total gas uptake going to mode/bin n
   !   sum_uprt = total gas-to-aerosol mass transfer rate (1/s)
   !   avg_uprt = average uptake rate over the time step, such that
   !              (gas mixing ratio)*avg_uprt is the total uptake (mol/mol/s)
   !----------------------------------------------------------------------
   pure subroutine aero_gasaerexch_partition( nbins, do_bin, uptkrate, diffus_ratio, deltatxx, &
                                              fgain, sum_uprt, avg_uprt )

      integer,  intent(in)  :: nbins             ! number of modes/bins
      logical,  intent(in)  :: do_bin(nbins)     ! mode/bin takes up the gas
      real(r8), intent(in)  :: uptkrate(nbins)   ! h2so4 gas-to-aerosol mass transfer rates (1/s)
      real(r8), intent(in)  :: diffus_ratio      ! gas diffusivity ratio relative to h2so4
      real(r8), intent(in)  :: deltatxx          ! time step (s), slightly increased
      real(r8), intent(out) :: fgain(nbins)      ! fraction of the gas uptake going to each mode/bin
      real(r8), intent(out) :: sum_uprt          ! total gas-to-aerosol mass transfer rate (1/s)
      real(r8), intent(out) :: avg_uprt          ! average uptake rate over the time step (1/s)

      integer :: n

      sum_uprt = 0.0_r8
      do n = 1, nbins
         if (do_bin(n)) then
            fgain(n) = uptkrate(n)*diffus_ratio
            sum_uprt = sum_uprt + fgain(n)
         else
            fgain(n) = 0.0_r8
         end if
      end do

      if (sum_uprt > 0.0_r8) then
         do n = 1, nbins
            fgain(n) = fgain(n) / sum_uprt
         end do
      end if
      ! at this point (sum_uprt <= 0.0) only when all the fgain are zero

      ! uptake amount (fraction of gas uptaken) over deltat
      avg_uprt = (1.0_r8 - exp(-deltatxx*sum_uprt))/deltatxx

   end subroutine aero_gasaerexch_partition

   !----------------------------------------------------------------------
   ! computes aerosol tendencies of a gas (e.g., h2so4) with reversible uptake
   ! in a grid cell: the gas relaxes toward an equilibrium mixing ratio given
   ! for each mode/bin (e.g., sulfeq above the tropopause)
   !----------------------------------------------------------------------
   pure subroutine aero_gasaerexch_reversible( nbins, do_bin, uptkrate, gas_equil, deltatxx, &
                                               sum_uprt, gas_bgn, aer_bgn, dqdt_aer, sum_dqdt_aer )

      integer,  intent(in)  :: nbins              ! number of modes/bins
      logical,  intent(in)  :: do_bin(nbins)      ! mode/bin takes up the gas
      real(r8), intent(in)  :: uptkrate(nbins)    ! gas-to-aerosol mass transfer rates (1/s)
      real(r8), intent(in)  :: gas_equil(nbins)   ! equilibrium gas mixing ratio over each mode/bin (mol/mol)
      real(r8), intent(in)  :: deltatxx           ! time step (s), slightly increased
      real(r8), intent(in)  :: sum_uprt           ! total gas-to-aerosol mass transfer rate (1/s)
      real(r8), intent(in)  :: gas_bgn            ! initial gas mixing ratio (mol/mol)
      real(r8), intent(in)  :: aer_bgn(nbins)     ! initial aerosol mixing ratios (mol/mol)
      real(r8), intent(inout) :: dqdt_aer(nbins)  ! aerosol tendencies (mol/mol/s); left unchanged where .not.do_bin
      real(r8), intent(out) :: sum_dqdt_aer       ! total aerosol tendency (mol/mol/s)

      integer  :: n
      real(r8) :: tmp_kxt, tmp_pxt
      real(r8) :: tmp_aer_end
      real(r8) :: tmp_gas_avg, tmp_gas_equ

      tmp_kxt = deltatxx*sum_uprt  ! sum over modes of uptake_rate*deltat
      tmp_pxt = 0.0_r8
      do n = 1, nbins
         if (.not. do_bin(n)) cycle
         tmp_pxt = tmp_pxt + uptkrate(n)*gas_equil(n)
      end do
      tmp_pxt = max( 0.0_r8, tmp_pxt*deltatxx )  ! sum over modes of uptake_rate*gas_equil*deltat
      ! calc avg gas over deltat
      if (tmp_kxt >= 1.0e-5_r8) then
         ! exponential decay towards equilibrium value solution
         tmp_gas_equ = tmp_pxt/tmp_kxt
         tmp_gas_avg = tmp_gas_equ + (gas_bgn-tmp_gas_equ)*(1.0_r8-exp(-tmp_kxt))/tmp_kxt
      else
         ! first order approx for tmp_kxt small
         tmp_gas_avg = gas_bgn*(1.0_r8-0.5_r8*tmp_kxt) + 0.5_r8*tmp_pxt
      end if
      sum_dqdt_aer = 0.0_r8
      do n = 1, nbins
         if (.not. do_bin(n)) cycle
         ! calc change to aerosol in mode n
         tmp_aer_end = aer_bgn(n) + deltatxx*uptkrate(n)*(tmp_gas_avg-gas_equil(n))
         tmp_aer_end = max( 0.0_r8, tmp_aer_end )
         dqdt_aer(n) = (tmp_aer_end - aer_bgn(n))/deltatxx
         sum_dqdt_aer = sum_dqdt_aer + dqdt_aer(n)
      end do

   end subroutine aero_gasaerexch_reversible

   !----------------------------------------------------------------------
   ! computes nh4 aerosol tendencies from nh3 gas uptake in a grid cell, with
   ! the nh4/so4 molar ratio limited to 2
   !----------------------------------------------------------------------
   pure subroutine aero_gasaerexch_nh4_limit( nbins, fgain_nh4, sum_dqdt_nh4, deltat, deltatxx, &
                                              qold_nh4, qold_so4, dqdt_so4, dqdt_nh4, sum_dqdt_nh4_b )

      integer,  intent(in)  :: nbins              ! number of modes/bins
      real(r8), intent(in)  :: fgain_nh4(nbins)   ! fraction of the nh3 uptake going to each mode/bin
      real(r8), intent(in)  :: sum_dqdt_nh4       ! total nh3 uptake without the limit (mol/mol/s)
      real(r8), intent(in)  :: deltat             ! time step (s)
      real(r8), intent(in)  :: deltatxx           ! time step (s), slightly increased
      real(r8), intent(in)  :: qold_nh4(nbins)    ! initial nh4 aerosol mixing ratios (mol/mol)
      real(r8), intent(in)  :: qold_so4(nbins)    ! initial so4 aerosol mixing ratios (mol/mol)
      real(r8), intent(in)  :: dqdt_so4(nbins)    ! so4 aerosol tendencies (mol/mol/s)
      real(r8), intent(out) :: dqdt_nh4(nbins)    ! nh4 aerosol tendencies (mol/mol/s)
      real(r8), intent(out) :: sum_dqdt_nh4_b     ! total nh3 uptake with the limit (mol/mol/s)

      integer  :: n
      real(r8) :: qmax_nh4, qnew_nh4, qnew_so4

      sum_dqdt_nh4_b = 0.0_r8
      do n = 1, nbins
         dqdt_nh4(n) = fgain_nh4(n)*sum_dqdt_nh4
         qnew_nh4 = qold_nh4(n) + dqdt_nh4(n)*deltat
         qnew_so4 = qold_so4(n) + dqdt_so4(n)*deltat
         qmax_nh4 = 2.0_r8*qnew_so4
         if (qnew_nh4 > qmax_nh4) then
            dqdt_nh4(n) = (qmax_nh4 - qold_nh4(n))/deltatxx
         end if
         sum_dqdt_nh4_b = sum_dqdt_nh4_b + dqdt_nh4(n)
      end do

   end subroutine aero_gasaerexch_nh4_limit

   !----------------------------------------------------------------------
   ! finds the species of a given species type (e.g., 's-organic') in each
   ! aerosol mode/bin
   !----------------------------------------------------------------------
   subroutine aero_gasaerexch_spec_ndx( aero_props, spectype, nspec_typ, spec_ndx )

      use aerosol_properties_mod, only: aerosol_properties, aero_name_len

      class(aerosol_properties), intent(in) :: aero_props  ! aerosol properties object
      character(len=*), intent(in) :: spectype             ! species type
      integer, intent(out) :: nspec_typ(:)                 ! number of species of spectype in each mode/bin
      integer, allocatable, intent(out) :: spec_ndx(:,:)   ! species indices within the mode/bin (nbins, maxval(nspec_typ))
                                                           ! of the species of spectype, in species order

      integer :: m, l, nbins
      character(len=aero_name_len) :: spectype_ml

      nbins = aero_props%nbins()

      nspec_typ(:) = 0
      do m = 1, nbins
         do l = 1, aero_props%nspecies(m)
            call aero_props%species_type(m, l, spectype_ml)
            if (trim(spectype_ml) == trim(spectype)) nspec_typ(m) = nspec_typ(m) + 1
         end do
      end do

      allocate(spec_ndx(nbins, maxval(nspec_typ(1:nbins))))
      spec_ndx(:,:) = -1

      nspec_typ(:) = 0
      do m = 1, nbins
         do l = 1, aero_props%nspecies(m)
            call aero_props%species_type(m, l, spectype_ml)
            if (trim(spectype_ml) == trim(spectype)) then
               nspec_typ(m) = nspec_typ(m) + 1
               spec_ndx(m,nspec_typ(m)) = l
            end if
         end do
      end do

   end subroutine aero_gasaerexch_spec_ndx

   !----------------------------------------------------------------------
   ! finds the constituent indices of the SOA gas (SOAG*) species, in
   ! constituent order; SOAGX* species are excluded (as in modal_aero_data)
   !----------------------------------------------------------------------
   subroutine aero_gasaerexch_soag_ndx( cnst_ndx )

      use constituents, only: pcnst, cnst_name

      integer, allocatable, intent(out) :: cnst_ndx(:)   ! constituent indices of SOA gas species

      integer :: i, n

      n = 0
      do i = 1, pcnst
         if (is_soag(cnst_name(i))) n = n + 1
      end do

      allocate(cnst_ndx(n))

      n = 0
      do i = 1, pcnst
         if (is_soag(cnst_name(i))) then
            n = n + 1
            cnst_ndx(n) = i
         end if
      end do

   contains

      pure logical function is_soag(name)
         character(len=*), intent(in) :: name
         is_soag = (name(:4) == 'SOAG' .and. name(:5) /= 'SOAGX')
      end function is_soag

   end subroutine aero_gasaerexch_soag_ndx

   !----------------------------------------------------------------------
   ! registers history fields of the gas-to-aerosol mass transfer rates
   ! (uptake rates), number densities and wet radii of each mode/bin
   !----------------------------------------------------------------------
   subroutine aero_gasaerexch_uptk_diag_init( nbins, history_aerosol )

      use cam_history, only: addfld, add_default, fieldname_len

      integer, intent(in) :: nbins             ! number of modes/bins
      logical, intent(in) :: history_aerosol   ! add fields to default history output

      integer :: m
      character(len=fieldname_len) :: fieldname
      character(len=128) :: long_name

      do m = 1, nbins
         write(fieldname,'("WETRAD_bin",I2.2)') m
         write(long_name,'("bin ",I2.2," wet radius in gas-aerosol exchange")') m
         call addfld(fieldname, (/'lev'/), 'A', 'cm', long_name )
         if ( history_aerosol ) then
            call add_default( fieldname, 1, ' ' )
         end if

         write(fieldname,'("UPTKRATE_bin",I2.2)') m
         write(long_name,'("bin ",I2.2," uptake rate in gas-aerosol exchange")') m
         call addfld(fieldname, (/'lev'/), 'A', 'sec-1', long_name )
         if ( history_aerosol ) then
            call add_default( fieldname, 1, ' ' )
         end if

         write(fieldname,'("NUMDENS_bin",I2.2)') m
         write(long_name,'("bin ",I2.2," number density in gas-aerosol exchange")') m
         call addfld(fieldname, (/'lev'/), 'A', 'm-3', long_name )
         if ( history_aerosol ) then
            call add_default( fieldname, 1, ' ' )
         end if
      end do

      call addfld('UPTKRATE', (/'lev'/), 'A', 'sec-1', 'total uptake rate in gas-aerosol exchange' )
      if ( history_aerosol ) then
         call add_default( 'UPTKRATE', 1, ' ' )
      end if

   end subroutine aero_gasaerexch_uptk_diag_init

   !----------------------------------------------------------------------
   ! outputs the gas-to-aerosol mass transfer rates (uptake rates), number
   ! densities and wet radii of each mode/bin
   !----------------------------------------------------------------------
   subroutine aero_gasaerexch_uptk_diag_out( lchnk, ncol, nbins, numdens, wetdiam, uptkrate )

      use cam_history, only: outfld, fieldname_len

      integer,  intent(in) :: lchnk             ! chunk identifier
      integer,  intent(in) :: ncol              ! number of columns
      integer,  intent(in) :: nbins             ! number of modes/bins
      real(r8), intent(in) :: numdens(:,:,:)    ! number density of each mode/bin (#/m3)
      real(r8), intent(in) :: wetdiam(:,:,:)    ! wet diameter of each mode/bin (m)
      real(r8), intent(in) :: uptkrate(:,:,:)   ! gas-to-aerosol mass transfer rates (1/s)

      integer :: m
      character(len=fieldname_len) :: fieldname
      real(r8) :: uptkrate_all(ncol,size(uptkrate,2))

      ! initialize to zero before summing
      uptkrate_all(:,:) = 0._r8

      do m = 1, nbins

         write(fieldname,'("NUMDENS_bin",I2.2)') m
         call outfld(fieldname, numdens(:ncol,:,m), ncol, lchnk )

         write(fieldname,'("WETRAD_bin",I2.2)') m
         call outfld(fieldname, 0.5_r8*wetdiam(:ncol,:,m)*1.e2_r8, ncol, lchnk ) ! wet radius (cm)

         write(fieldname,'("UPTKRATE_bin",I2.2)') m
         call outfld(fieldname, uptkrate(:ncol,:,m), ncol, lchnk )

         uptkrate_all(:ncol,:) = uptkrate_all(:ncol,:) + uptkrate(:ncol,:,m)
      end do

      call outfld('UPTKRATE', uptkrate_all(:ncol,:), ncol, lchnk )

   end subroutine aero_gasaerexch_uptk_diag_out

   !----------------------------------------------------------------------
   ! configures and registers the SOA condensation/evaporation history fields
   !----------------------------------------------------------------------
   subroutine soa_condevap_init( self, nsoa, do_total, do_spec, history_aerosol, &
                                 grp, grp_sfx, grp_desc )

      use cam_history,    only: addfld, add_default, fieldname_len
      use cam_abortutils, only: endrun

      class(soa_condevap_diag_t), intent(inout) :: self
      integer, intent(in) :: nsoa                  ! number of soa species
      logical, intent(in) :: do_total              ! output totals over all soa species
      logical, intent(in) :: do_spec               ! output per soa species
      logical, intent(in) :: history_aerosol       ! add fields to default history output
      integer, optional, intent(in) :: grp(:)      ! group index of each soa species (0 = none)
      character(len=*), optional, intent(in) :: grp_sfx(:)  ! field name suffix of each group
      character(len=*), optional, intent(in) :: grp_desc(:) ! description of each group

      integer :: j
      character(len=2) :: outsoa
      character(len=fieldname_len) :: fieldname
      character(len=128) :: long_name
      character(len=*), parameter :: subname = 'aero_gasaerexch::soa_condevap_init'

      self%nsoa = nsoa
      self%do_total = do_total
      self%do_spec = do_spec

      allocate(self%grp(nsoa))
      self%grp(:) = 0
      self%ngrp = 0
      if (present(grp)) then
         if (.not.(present(grp_sfx) .and. present(grp_desc))) then
            call endrun(subname//': grp_sfx and grp_desc must be given with grp')
         end if
         if (size(grp) /= nsoa .or. size(grp_desc) /= size(grp_sfx)) then
            call endrun(subname//': inconsistent group arguments')
         end if
         self%ngrp = size(grp_sfx)
         if (any(grp < 0) .or. any(grp > self%ngrp)) then
            call endrun(subname//': invalid group index')
         end if
         self%grp(:) = grp(:)
         allocate(self%grp_sfx(self%ngrp))
         self%grp_sfx(:) = grp_sfx(:)
      else
         allocate(self%grp_sfx(0))
      end if

      if (self%do_total) then
         call add_condevap_fld('qcon_gaex', '3D fields for SOA condensation')
         call add_condevap_fld('qevap_gaex', '3D fields for SOA evaporation')
      end if

      if (self%do_spec) then
         do j = 1, nsoa
            write (outsoa, "(I2.2)") j
            call add_condevap_fld('qcon_gaex'//outsoa, '3D fields for SOA condensation for VBS bin'//outsoa)
            call add_condevap_fld('qevap_gaex'//outsoa, '3D fields for SOA evaporation for VBS bin'//outsoa)
         end do
      end if

      do j = 1, self%ngrp
         call add_condevap_fld('qcon'//trim(self%grp_sfx(j))//'_gaex', &
                               '3D fields for '//trim(grp_desc(j))//' SOA condensation')
         call add_condevap_fld('qevap'//trim(self%grp_sfx(j))//'_gaex', &
                               '3D fields for '//trim(grp_desc(j))//' SOA evaporation')
      end do

   contains

      subroutine add_condevap_fld(name, desc)
         character(len=*), intent(in) :: name
         character(len=*), intent(in) :: desc

         fieldname = name
         long_name = desc
         call addfld(fieldname, (/'lev'/), 'A', condevap_units, long_name )
         if ( history_aerosol ) then
            call add_default( fieldname, 1, ' ' )
         end if

      end subroutine add_condevap_fld

   end subroutine soa_condevap_init

   !----------------------------------------------------------------------
   ! allocates and zeros the SOA condensation/evaporation sums for a chunk
   !----------------------------------------------------------------------
   subroutine soa_condevap_reset( self, acc )

      use ppgrid, only: pcols, pver

      class(soa_condevap_diag_t), intent(in) :: self
      type(soa_condevap_acc_t), intent(out) :: acc   ! sums for a chunk

      allocate(acc%con_tot(pcols,pver), acc%evap_tot(pcols,pver))
      allocate(acc%con_spc(pcols,pver,self%nsoa), acc%evap_spc(pcols,pver,self%nsoa))
      allocate(acc%con_grp(pcols,pver,self%ngrp), acc%evap_grp(pcols,pver,self%ngrp))

      acc%con_tot(:,:) = 0.0_r8
      acc%evap_tot(:,:) = 0.0_r8
      acc%con_spc(:,:,:) = 0.0_r8
      acc%evap_spc(:,:,:) = 0.0_r8
      acc%con_grp(:,:,:) = 0.0_r8
      acc%evap_grp(:,:,:) = 0.0_r8

   end subroutine soa_condevap_reset

   !----------------------------------------------------------------------
   ! adds the SOA tendencies of one grid cell to the condensation
   ! (dqdt_soa >= 0) and evaporation (dqdt_soa < 0) sums
   !----------------------------------------------------------------------
   subroutine soa_condevap_accum( self, acc, i, k, nbins, mask, fac, dqdt_soa )

      class(soa_condevap_diag_t), intent(in) :: self
      type(soa_condevap_acc_t), intent(inout) :: acc     ! sums for a chunk
      integer,  intent(in) :: i                          ! column index
      integer,  intent(in) :: k                          ! level index
      integer,  intent(in) :: nbins                      ! number of modes/bins
      logical,  intent(in) :: mask(nbins,self%nsoa)      ! include mode/bin soa species
      real(r8), intent(in) :: fac(nbins,self%nsoa)       ! conversion factor to mass mixing ratio (mw/mwdry)
      real(r8), intent(in) :: dqdt_soa(nbins,self%nsoa)  ! soa tendencies (mol/mol/s)

      integer :: n, jsoa, jg

      ! each sum is updated as sum + dqdt_soa*fac (not via a shared product
      ! temporary) so that results do not depend on multiply-add contraction
      ! differently than the per-model code this replaces
      do n = 1, nbins
         do jsoa = 1, self%nsoa
            if (.not. mask(n,jsoa)) cycle
            jg = self%grp(jsoa)
            if (dqdt_soa(n,jsoa) >= 0.0_r8) then
               if (self%do_spec) acc%con_spc(i,k,jsoa) = acc%con_spc(i,k,jsoa) + dqdt_soa(n,jsoa)*fac(n,jsoa)
               if (self%do_total) acc%con_tot(i,k) = acc%con_tot(i,k) + dqdt_soa(n,jsoa)*fac(n,jsoa)
               if (jg > 0) acc%con_grp(i,k,jg) = acc%con_grp(i,k,jg) + dqdt_soa(n,jsoa)*fac(n,jsoa)
            else if (dqdt_soa(n,jsoa) < 0.0_r8) then
               if (self%do_spec) acc%evap_spc(i,k,jsoa) = acc%evap_spc(i,k,jsoa) + dqdt_soa(n,jsoa)*fac(n,jsoa)
               if (self%do_total) acc%evap_tot(i,k) = acc%evap_tot(i,k) + dqdt_soa(n,jsoa)*fac(n,jsoa)
               if (jg > 0) acc%evap_grp(i,k,jg) = acc%evap_grp(i,k,jg) + dqdt_soa(n,jsoa)*fac(n,jsoa)
            end if
         end do
      end do

   end subroutine soa_condevap_accum

   !----------------------------------------------------------------------
   ! adds the soa tendencies of all grid cells of a chunk (levels top_lev
   ! to nlev) to the condensation and evaporation sums (see accum)
   !----------------------------------------------------------------------
   subroutine soa_condevap_accum_column( self, acc, ncol, nlev, top_lev, nbins, mask, fac, dqdt_soa )

      class(soa_condevap_diag_t), intent(in) :: self
      type(soa_condevap_acc_t), intent(inout) :: acc     ! sums for a chunk
      integer,  intent(in) :: ncol                       ! number of columns
      integer,  intent(in) :: nlev                       ! number of vertical levels
      integer,  intent(in) :: top_lev                    ! top level index for computation
      integer,  intent(in) :: nbins                      ! number of modes/bins
      logical,  intent(in) :: mask(nbins,self%nsoa)      ! include mode/bin soa species
      real(r8), intent(in) :: fac(nbins,self%nsoa)       ! conversion factor to mass mixing ratio (mw/mwdry)
      real(r8), intent(in) :: dqdt_soa(:,:,:,:)          ! soa tendencies (mol/mol/s) (i,k,n,jsoa)

      integer :: i, k
      real(r8) :: dqdt_cell(nbins,self%nsoa)

      do k = top_lev, nlev
         do i = 1, ncol
            dqdt_cell(:,:) = dqdt_soa(i,k,:,:)
            call self%accum( acc, i, k, nbins, mask, fac, dqdt_cell )
         end do
      end do

   end subroutine soa_condevap_accum_column

   !----------------------------------------------------------------------
   ! outputs the SOA condensation/evaporation history fields
   !----------------------------------------------------------------------
   subroutine soa_condevap_outfld( self, acc, lchnk )

      use cam_history, only: outfld
      use ppgrid,      only: pcols

      class(soa_condevap_diag_t), intent(in) :: self
      type(soa_condevap_acc_t), intent(in) :: acc   ! sums for a chunk
      integer, intent(in) :: lchnk                  ! chunk identifier

      integer :: j
      character(len=2) :: outsoa

      if (self%do_total) then
         call outfld('qcon_gaex', acc%con_tot(:,:), pcols, lchnk )
         call outfld('qevap_gaex', acc%evap_tot(:,:), pcols, lchnk )
      end if

      if (self%do_spec) then
         do j = 1, self%nsoa
            write (outsoa, "(I2.2)") j
            call outfld('qcon_gaex'//outsoa, acc%con_spc(:,:,j), pcols, lchnk )
            call outfld('qevap_gaex'//outsoa, acc%evap_spc(:,:,j), pcols, lchnk )
         end do
      end if

      do j = 1, self%ngrp
         call outfld('qcon'//trim(self%grp_sfx(j))//'_gaex', acc%con_grp(:,:,j), pcols, lchnk )
         call outfld('qevap'//trim(self%grp_sfx(j))//'_gaex', acc%evap_grp(:,:,j), pcols, lchnk )
      end do

   end subroutine soa_condevap_outfld

   !----------------------------------------------------------------------
   ! registers the column tendency history fields '<spec_names(i)><suffix>'
   ! (field index i in calls to outfld)
   !----------------------------------------------------------------------
   subroutine col_tend_init( self, spec_names, suffix, desc, history_aerosol )

      use cam_history, only: addfld, add_default, horiz_only

      class(col_tend_diag_t), intent(inout) :: self
      character(len=*), intent(in) :: spec_names(:)  ! species names
      character(len=*), intent(in) :: suffix         ! field name suffix
      character(len=*), intent(in) :: desc           ! description appended to the species name in long name
      logical,          intent(in) :: history_aerosol ! add fields to default history output

      integer :: i
      character(len=8) :: unit
      character(len=128) :: long_name

      self%nfld = size(spec_names)
      allocate(self%fldname(self%nfld))

      do i = 1, self%nfld
         self%fldname(i) = trim(spec_names(i)) // suffix
         long_name = trim(spec_names(i)) // ' ' // desc
         unit = 'kg/m2/s'
         if ((spec_names(i)(1:3) == 'num') .or. (spec_names(i)(1:3) == 'NUM')) unit = '#/m2/s'
         call addfld( self%fldname(i), horiz_only, 'A', unit, long_name )
         if ( history_aerosol ) then
            call add_default( self%fldname(i), 1, ' ' )
         end if
      end do

   end subroutine col_tend_init

   !----------------------------------------------------------------------
   ! converts a column tendency to the output units and outputs it
   !----------------------------------------------------------------------
   subroutine col_tend_outfld( self, ifld, ncol, lchnk, fac, qsrflx )

      use cam_history, only: outfld
      use ppgrid,      only: pcols

      class(col_tend_diag_t), intent(in) :: self
      integer,  intent(in) :: ifld              ! field index (from init)
      integer,  intent(in) :: ncol              ! number of columns
      integer,  intent(in) :: lchnk             ! chunk identifier
      real(r8), intent(in) :: fac               ! conversion factor to output units (molecular weight / mwdry)
      real(r8), intent(inout) :: qsrflx(pcols)  ! column tendency; converted to output units on return

      integer :: i

      do i = 1, ncol
         qsrflx(i) = qsrflx(i)*fac
      end do
      call outfld( trim(self%fldname(ifld)), qsrflx, pcols, lchnk )

   end subroutine col_tend_outfld

end module aero_gasaerexch
