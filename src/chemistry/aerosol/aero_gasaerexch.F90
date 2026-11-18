!----------------------------------------------------------------------
! CAM-side gas-aerosol exchange (gas condensation onto aerosol particles)
! common to the aerosol models (MAM, CARMA): configuration from the
! abstract aerosol_properties interface, constituent lookups and history
! diagnostics. The portable science is in aero_gasaerexch_core
! (atmos_phys); model-specific bookkeeping (constituent indices, renaming,
! aging, ...) stays with the calling aerosol model.
!----------------------------------------------------------------------
module aero_gasaerexch

   use shr_kind_mod, only: r8 => shr_kind_r8
   use aero_gasaerexch_core, only: method_soa_none, method_soa_irreversible, method_soa_reversible

   implicit none
   private

   public :: aero_gasaerexch_uptkrates
   public :: aero_gasaerexch_uptk_diag_init
   public :: aero_gasaerexch_uptk_diag_out
   public :: aero_gasaerexch_spec_ndx
   public :: aero_gasaerexch_soag_ndx

   ! methods for soa gas uptake (see aero_gasaerexch_core)
   public :: method_soa_none
   public :: method_soa_irreversible
   public :: method_soa_reversible

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

      use aero_gasaerexch_core, only: aero_gasaerexch_soa_column
      use mo_constants,         only: rgas ! Gas constant (J/K/mol)
      use cam_abortutils,       only: endrun

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
      character(len=256) :: errmsg
      integer :: errflg

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

      call aero_gasaerexch_soa_column( method_soa, ncol, nlev, top_lev, deltat, temp, pmid, rgas, &
                                       self%nbins, self%nsoa, self%npoa, do_bin, do_soa, do_soag, uptkrate, &
                                       self%mw_poa_host, self%mw_soa_host, soag, soa, poa, &
                                       dqdt_soa, sum_dqdt_soa, errmsg, errflg )
      if (errflg /= 0) then
         call endrun('aero_gasaerexch::soa_exch_tend: '//trim(errmsg))
      end if

   end subroutine soa_exch_tend

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
