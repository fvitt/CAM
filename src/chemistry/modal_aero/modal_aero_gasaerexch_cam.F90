! CAM wrapper for modal_aero_gasaerexch.
! Resolves CAM-specific species indices and mode metadata from the aerosol
! properties object, calls the portable modal_aero_gasaerexch_init, and
! registers and outputs history fields through the shared gas-aerosol
! exchange diagnostics (aero_gasaerexch).
module modal_aero_gasaerexch_cam
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use aero_gasaerexch, only: soa_condevap_diag_t, col_tend_diag_t, soa_exch_t

   implicit none
   private

   public :: modal_aero_gasaerexch_cam_init
   public :: modal_aero_gasaerexch_cam_outfld

   ! mode and species information, set from the aerosol properties object
   integer :: nmodes = 0                         ! number of modes
   integer :: nspec_max = 0                      ! maximum number of species in a mode
   integer, allocatable :: nspec_mode(:)         ! number of species in each mode
   integer, allocatable :: num_a_ndx(:)          ! constituent index of interstitial number in each mode
   integer, allocatable :: mass_a_ndx(:,:)       ! constituent index of interstitial species in each mode
   integer :: mode_pcarbon = -1                  ! primary carbon mode (<= 0 if none)
   integer :: mode_accum = -1                    ! accumulation mode (<= 0 if none)
   character(len=32), allocatable :: cw_cnst_name(:) ! cloud-borne name of each aerosol constituent

   ! soa species configuration and numbers of soa and poa species (maximum
   ! over the modes), set from the aerosol properties object
   type(soa_exch_t) :: soa_exch
   integer :: nsoa = 0
   integer :: npoa = 0

   ! constituent indices of interstitial aerosol species in each mode
   ! (set from the aerosol properties object; <= 0 if species is absent)
   integer, allocatable :: so4_a_ndx(:)      ! sulfate
   integer, allocatable :: nh4_a_ndx(:)      ! ammonium
   integer, allocatable :: soa_a_ndx(:,:)    ! secondary organic aerosol species
   integer, allocatable :: pom_a_ndx(:,:)    ! primary organic matter species

   type(soa_condevap_diag_t) :: soa_diag     ! SOA condensation/evaporation diagnostics

   ! column tendency diagnostics and their field index for each constituent (0 if none)
   type(col_tend_diag_t) :: sfgaex1_diag     ! gas-aerosol exchange ('_sfgaex1')
   type(col_tend_diag_t) :: sfgaex2_diag     ! renaming of interstitial species ('_sfgaex2')
   type(col_tend_diag_t) :: sfgaex2cw_diag   ! renaming of cloud-borne species ('_sfgaex2')
   integer, allocatable :: sfgaex1_fld(:)
   integer, allocatable :: sfgaex2_fld(:)
   integer, allocatable :: sfgaex2cw_fld(:)

contains

   ! Resolve species indices and mode metadata from the aerosol properties
   ! object, call portable modal_aero_gasaerexch_init, and register history
   ! fields for gas-aerosol exchange diagnostics.
   subroutine modal_aero_gasaerexch_cam_init(aero_props)
      use modal_aero_gasaerexch, only: modal_aero_gasaerexch_init
      use modal_aero_rename_cam, only: npair_renamexf, nspecfrm_renamexf, &
                                       lspecfrma_renamexf, lspectooa_renamexf, &
                                       lspecfrmc_renamexf, lspectooc_renamexf
      use aerosol_properties_mod, only: aerosol_properties
      use constituents, only: pcnst, cnst_name, cnst_get_ind
      use physconst, only: rair, mwdry, r_universal
      use cam_abortutils, only: endrun
      use spmd_utils, only: masterproc
      use phys_control, only: phys_getopts, cam_chempkg_is
      use cam_logfile, only: iulog

      class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object

      ! indices are in pcnst (q) space
      ! run phase receives loffset for gas_pcnst (vmr) space

      ! local
      integer                          :: ipair, iq, iqfrm, iqtoo
      integer                          :: jsoa
      integer                          :: l, lsfrm, lstoo
      integer                          :: mfrm, mtoo
      integer                          :: n
      integer                          :: nchfrm, nchfrmskip, nchtoo, nchtooskip

      logical                          :: dotend(pcnst), dotendqqcw(pcnst)

      logical                          :: history_aerosol

      ! Resolved species indices (pcnst-space, for portable _init)
      integer                          :: idx_h2so4, idx_nh3, idx_msa
      integer, allocatable             :: idx_soag(:)

      ! pcage resolved arrays (pcnst-space)
      integer                          :: nspecfrm_pcage_resolved
      integer, allocatable             :: lspecfrm_pcage_resolved(:)
      integer, allocatable             :: lspectoo_pcage_resolved(:)

      ! mode and species properties from the aerosol properties object
      real(r8), allocatable            :: alnsg(:)
      real(r8), allocatable            :: specmw(:,:), specdens(:,:), spechygro(:,:)

      character(len=256) :: errmsg
      integer            :: errflg

      call phys_getopts(history_aerosol_out=history_aerosol)

      !-----------------------------------------------------------------------
      ! Part A: Resolve arguments and call portable _init
      !-----------------------------------------------------------------------

      ! mode tables, species indices and numbers of soa/poa species (nsoa, npoa)
      call set_mode_tables(aero_props)
      call set_species_ndx(aero_props)

      ! --- Gas-phase species indices (pcnst-space) ---
      call cnst_get_ind('H2SO4', idx_h2so4, abort=.false.)
      if ((idx_h2so4 <= 0) .or. (idx_h2so4 > pcnst)) then
         write (iulog, '(/a/a,i7)') &
            '*** modal_aero_gasaerexch_cam_init -- cannot find H2SO4 species', &
            '    idx_h2so4=', idx_h2so4
         call endrun('modal_aero_gasaerexch_cam_init error: H2SO4 not found')
      end if

      call cnst_get_ind('NH3', idx_nh3, abort=.false.)

      if (.not. cam_chempkg_is('geoschem_mam4')) then
         call cnst_get_ind('MSA', idx_msa, abort=.false.)
      else
         idx_msa = 0
      end if

      ! --- SOA gas-phase species indices (pcnst-space) ---
      allocate(idx_soag(nsoa))
      do jsoa = 1, nsoa
         l = soa_exch%soag_cnst(jsoa)
         if ((l > 0) .and. (l <= pcnst)) then
            idx_soag(jsoa) = l
         else
            idx_soag(jsoa) = 0
         end if
      end do

      ! --- mode and species properties ---
      allocate(alnsg(nmodes))
      allocate(specmw(nspec_max,nmodes), specdens(nspec_max,nmodes), spechygro(nspec_max,nmodes))
      specmw(:,:) = 0.0_r8
      specdens(:,:) = 0.0_r8
      spechygro(:,:) = 0.0_r8
      do n = 1, nmodes
         alnsg(n) = aero_props%alogsig(n)
         do l = 1, nspec_mode(n)
            call aero_props%get(n, l, spec_mw=specmw(l,n), density=specdens(l,n), hygro=spechygro(l,n))
         end do
      end do

      ! --- Resolve pcage species matching ---
      ! define "from mode" and "to mode" for primary carbon aging,
      ! skip (turn off) aging if either is absent, or if accum mode so4 is absent
      ! Indices stored in pcnst-space.
      allocate(lspecfrm_pcage_resolved(nspec_max))
      allocate(lspectoo_pcage_resolved(nspec_max))
      nspecfrm_pcage_resolved = 0
      lspecfrm_pcage_resolved(:) = 0
      lspectoo_pcage_resolved(:) = 0

      if ((mode_pcarbon > 0) .and. (mode_accum > 0)) then
         l = so4_a_ndx(mode_accum)
         if ((l >= 1) .and. (l <= pcnst)) then

            mfrm = mode_pcarbon
            mtoo = mode_accum

            if (mfrm < 10) then
               nchfrmskip = 1
            else if (mfrm < 100) then
               nchfrmskip = 2
            else
               nchfrmskip = 3
            end if
            if (mtoo < 10) then
               nchtooskip = 1
            else if (mtoo < 100) then
               nchtooskip = 2
            else
               nchtooskip = 3
            end if

            aa_iqfrm: do iqfrm = -1, nspec_mode(mfrm)

               if (iqfrm == -1) then
                  lsfrm = num_a_ndx(mfrm)
                  lstoo = num_a_ndx(mtoo)
               else if (iqfrm == 0) then
                  ! bypass transfer of aerosol water due to primary-carbon aging
                  cycle aa_iqfrm
               else
                  lsfrm = mass_a_ndx(iqfrm, mfrm)
                  lstoo = 0
               end if
               if ((lsfrm < 1) .or. (lsfrm > pcnst)) cycle aa_iqfrm

               if (lsfrm > 0 .and. iqfrm > 0) then
                  nchfrm = len(trim(cnst_name(lsfrm))) - nchfrmskip

                  ! find "too" species having same cnst_name (except for last 1/2/3 characters which are the mode index)
                  do iqtoo = 1, nspec_mode(mtoo)
                     lstoo = mass_a_ndx(iqtoo, mtoo)
                     nchtoo = len(trim(cnst_name(lstoo))) - nchtooskip
                     if (cnst_name(lsfrm) (1:nchfrm) == cnst_name(lstoo) (1:nchtoo)) then
                        exit
                     else
                        lstoo = 0
                     end if
                  end do
               end if

               if ((lstoo < 1) .or. (lstoo > pcnst)) lstoo = 0
               nspecfrm_pcage_resolved = nspecfrm_pcage_resolved + 1
               lspecfrm_pcage_resolved(nspecfrm_pcage_resolved) = lsfrm
               lspectoo_pcage_resolved(nspecfrm_pcage_resolved) = lstoo
            end do aa_iqfrm

            ! output results
            if (masterproc) then
               write (iulog, 9310)
               write (iulog, 9320) 1, mfrm, mtoo

               do iq = 1, nspecfrm_pcage_resolved
                  lsfrm = lspecfrm_pcage_resolved(iq)
                  lstoo = lspectoo_pcage_resolved(iq)
                  if (lstoo > 0) then
                     write (iulog, 9330) lsfrm, cnst_name(lsfrm), &
                        lstoo, cnst_name(lstoo)
                  else
                     write (iulog, 9340) lsfrm, cnst_name(lsfrm)
                  end if
               end do

               write (iulog, *)
            end if ! ( masterproc )

9310        format(/'subr. modal_aero_gasaerexch_cam_init - primary carbon aging pointers')
9320        format('pair', i3, 5x, 'mode', i3, ' ---> mode', i3)
9330        format(5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a)
9340        format(5x, 'spec', i3, '=', a, ' ---> LOSS')

         end if ! so4_a_ndx(mode_accum) valid
      end if ! mode_pcarbon > 0 .and. mode_accum > 0

      ! --- Call portable init ---
      call modal_aero_gasaerexch_init( &
         ntot_amode         = nmodes,                    &
         nsoa               = nsoa,                      &
         npoa               = npoa,                      &
         nspec_max          = nspec_max,                 &
         nspec_amode        = nspec_mode,                &
         modeptr_pcarbon    = mode_pcarbon,              &
         modeptr_accum      = mode_accum,                &
         alnsg_amode        = alnsg,                     &
         specmw_amode       = specmw,                    &
         specdens_amode     = specdens,                  &
         spechygro          = spechygro,                 &
         idx_h2so4          = idx_h2so4,                 &
         idx_nh3            = idx_nh3,                   &
         idx_msa            = idx_msa,                   &
         idx_soag           = idx_soag,                  &
         idx_so4_a          = so4_a_ndx,                 &
         idx_nh4_a          = nh4_a_ndx,                 &
         idx_soa_a          = soa_a_ndx,                 &
         idx_pom_a          = pom_a_ndx,                 &
         idx_num            = num_a_ndx,                 &
         idx_mass           = mass_a_ndx,                &
         pcnst_in           = pcnst,                     &
         nspecfrm_pcage_in  = nspecfrm_pcage_resolved,   &
         lspecfrm_pcage_in  = lspecfrm_pcage_resolved,   &
         lspectoo_pcage_in  = lspectoo_pcage_resolved,   &
         mw_soa_host        = soa_exch%mw_soa_host,      &
         mw_poa_host        = soa_exch%mw_poa_host,      &
         rair               = rair,                      &
         mwdry              = mwdry,                     &
         r_universal        = r_universal,               &
         errmsg             = errmsg,                    &
         errflg             = errflg)

      if (errflg /= 0) then
         call endrun('modal_aero_gasaerexch_cam_init: '//trim(errmsg))
      end if

      !-----------------------------------------------------------------------
      ! Part B: History field registration (addfld calls)
      !-----------------------------------------------------------------------

      ! --- Tendency flags for _sfgaex1 fields ---
      ! Determine which constituents get gas-aerosol exchange tendency output
      dotend(:) = .false.

      ! H2SO4 (required). below indices are in pcnst indexing:
      dotend(idx_h2so4) = .true.

      ! NH3 (optional)
      if ((idx_nh3 > 0) .and. (idx_nh3 <= pcnst)) dotend(idx_nh3) = .true.

      ! MSA (optional)
      if ((idx_msa > 0) .and. (idx_msa <= pcnst)) dotend(idx_msa) = .true.

      ! SOA gases
      do jsoa = 1, nsoa
         l = idx_soag(jsoa)
         if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
      end do

      ! Aerosol species in modes
      do n = 1, nmodes
         l = so4_a_ndx(n)
         if ((l > 0) .and. (l <= pcnst)) then
            dotend(l) = .true.
            if (idx_nh3 > 0) then
               l = nh4_a_ndx(n)
               if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
            end if
         end if
         do jsoa = 1, nsoa
            if (idx_soag(jsoa) > 0) then
               l = soa_a_ndx(n, jsoa)
               if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
            end if
         end do
      end do

      ! pcage species
      do iq = 1, nspecfrm_pcage_resolved
         lsfrm = lspecfrm_pcage_resolved(iq)
         lstoo = lspectoo_pcage_resolved(iq)
         if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
            dotend(lsfrm) = .true.
            if ((lstoo > 0) .and. (lstoo <= pcnst)) then
               dotend(lstoo) = .true.
            end if
         end if
      end do

      ! --- SOA condensation/evaporation diagnostics ---
      !   nsoa=5:  totals over all SOA species
      !   nsoa=15: fossil fuel (SOA species 1-5), biomass burning (6-10) and biogenic (11-15) groups
      if (nsoa == 5) then
         call soa_diag%init( nsoa, do_total=.true., do_spec=.false., history_aerosol=history_aerosol )
      else if (nsoa == 15) then
         call soa_diag%init( nsoa, do_total=.false., do_spec=.false., history_aerosol=history_aerosol, &
                             grp=(/ (1+(jsoa-1)/5, jsoa=1,nsoa) /), &
                             grp_sfx=(/ 'ff', 'bb', 'bg' /), &
                             grp_desc=(/ 'Fossil  ', 'Biomass ', 'Biogenic' /) )
      else
         call soa_diag%init( nsoa, do_total=.false., do_spec=.false., history_aerosol=history_aerosol )
      end if

      ! --- Per-species _sfgaex1 fields (gas-aerosol exchange and primary carbon aging) ---
      call init_col_tend( sfgaex1_diag, sfgaex1_fld, dotend, cnst_name, '_sfgaex1', &
                          'gas-aerosol-exchange primary column tendency', history_aerosol )

      ! --- Per-species _sfgaex2 fields (aitken-->accum renaming tendency) ---
      dotend(:) = .false.
      dotendqqcw(:) = .false.
      do ipair = 1, npair_renamexf
         do iq = 1, nspecfrm_renamexf(ipair)
            lsfrm = lspecfrma_renamexf(iq, ipair)
            lstoo = lspectooa_renamexf(iq, ipair)
            if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotend(lsfrm) = .true.
               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                  dotend(lstoo) = .true.
               end if
            end if

            lsfrm = lspecfrmc_renamexf(iq, ipair)
            lstoo = lspectooc_renamexf(iq, ipair)
            if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
               dotendqqcw(lsfrm) = .true.
               if ((lstoo > 0) .and. (lstoo <= pcnst)) then
                  dotendqqcw(lstoo) = .true.
               end if
            end if
         end do ! iq = ...
      end do ! ipair = ...

      call init_col_tend( sfgaex2_diag, sfgaex2_fld, dotend, cnst_name, '_sfgaex2', &
                          'gas-aerosol-exchange renaming column tendency', history_aerosol )
      call init_col_tend( sfgaex2cw_diag, sfgaex2cw_fld, dotendqqcw, cw_cnst_name, '_sfgaex2', &
                          'gas-aerosol-exchange renaming column tendency', history_aerosol )

   end subroutine modal_aero_gasaerexch_cam_init

   ! Output history fields for gas-aerosol exchange and renaming
   ! (column tendencies) and SOA condensation/evaporation.
   subroutine modal_aero_gasaerexch_cam_outfld(lchnk, ncol, loffset, dotend, dotendrn, dotendqqcwrn, &
                                               qsrflx_gaexch, qsrflx_rename, qqcwsrflx_rename, dqdt_soa_mode)
      use ppgrid,          only: pver
      use ref_pres,        only: top_lev => clim_modal_aero_top_lev
      use chem_mods,       only: adv_mass, gas_pcnst
      use physconst,       only: mwdry
      use aero_gasaerexch, only: soa_condevap_acc_t

      integer,  intent(in)    :: lchnk                    ! chunk identifier
      integer,  intent(in)    :: ncol                     ! number of columns
      integer,  intent(in)    :: loffset                  ! offset applied to modal aero "pointers"
      logical,  intent(in)    :: dotend(:)                ! (gas_pcnst) species with gas-aerosol exchange tendencies
      logical,  intent(in)    :: dotendrn(:)              ! (gas_pcnst) interstitial species with renaming tendencies
      logical,  intent(in)    :: dotendqqcwrn(:)          ! (gas_pcnst) cloud-borne species with renaming tendencies
      real(r8), intent(inout) :: qsrflx_gaexch(:,:)       ! (pcols,gas_pcnst) gas-aerosol exchange column tendencies
                                                          ! (converted to output units on return)
      real(r8), intent(inout) :: qsrflx_rename(:,:)       ! (pcols,gas_pcnst) interstitial renaming column tendencies
                                                          ! (converted to output units on return)
      real(r8), intent(inout) :: qqcwsrflx_rename(:,:)    ! (pcols,gas_pcnst) cloud-borne renaming column tendencies
                                                          ! (converted to output units on return)
      real(r8), intent(in)    :: dqdt_soa_mode(:,:,:,:)   ! (ncol,pver,nmodes,nsoa) soa aerosol tendencies of
                                                          ! each mode from soa gas uptake (mol/mol/s)

      integer  :: l, lb, n, jsoa
      logical  :: condevap_mask(nmodes,nsoa)
      real(r8) :: condevap_fac(nmodes,nsoa)
      type(soa_condevap_acc_t) :: condevap   ! SOA condensation/evaporation diagnostics

      ! SOA condensation/evaporation of soa species present in modes when the
      ! soa gas is present (excludes the primary carbon mode, which has no soa species)
      call soa_diag%reset(condevap)
      do jsoa = 1, nsoa
         do n = 1, nmodes
            l = soa_a_ndx(n,jsoa) - loffset
            condevap_mask(n,jsoa) = (soa_exch%soag_cnst(jsoa) - loffset > 0) .and. &
                                    (l > 0) .and. (l <= gas_pcnst)
            if (condevap_mask(n,jsoa)) then
               condevap_fac(n,jsoa) = adv_mass(l)/mwdry
            else
               condevap_fac(n,jsoa) = 0.0_r8
            end if
         end do
      end do
      call soa_diag%accum_column( condevap, ncol, pver, top_lev, nmodes, condevap_mask, condevap_fac, &
                                  dqdt_soa_mode )
      call soa_diag%outfld( condevap, lchnk )

      ! history file column-tendency fields
      ! (lb = l + loffset is a constituent index only for species with tendencies,
      !  so the field index arrays are only accessed for those)
      do l = 1, gas_pcnst
         lb = l + loffset
         if ( dotend(l) ) then
            if ( sfgaex1_fld(lb) > 0 ) then
               call sfgaex1_diag%outfld( sfgaex1_fld(lb), ncol, lchnk, adv_mass(l)/mwdry, &
                                         qsrflx_gaexch(:,l) )
            end if
         end if
         if ( dotendrn(l) ) then
            if ( sfgaex2_fld(lb) > 0 ) then
               call sfgaex2_diag%outfld( sfgaex2_fld(lb), ncol, lchnk, adv_mass(l)/mwdry, &
                                         qsrflx_rename(:,l) )
            end if
         end if
         if ( dotendqqcwrn(l) ) then
            if ( sfgaex2cw_fld(lb) > 0 ) then
               call sfgaex2cw_diag%outfld( sfgaex2cw_fld(lb), ncol, lchnk, adv_mass(l)/mwdry, &
                                           qqcwsrflx_rename(:,l) )
            end if
         end if
      end do

   end subroutine modal_aero_gasaerexch_cam_outfld

   ! registers the column tendency fields of the constituents selected by mask
   ! and sets the field index of each constituent (0 if none)
   subroutine init_col_tend( diag, fld_ndx, mask, cnames, suffix, desc, history_aerosol )

      type(col_tend_diag_t), intent(inout) :: diag
      integer, allocatable, intent(out) :: fld_ndx(:)  ! field index of each constituent
      logical, intent(in) :: mask(:)                   ! constituents with column tendency fields
      character(len=*), intent(in) :: cnames(:)        ! constituent names
      character(len=*), intent(in) :: suffix           ! field name suffix
      character(len=*), intent(in) :: desc             ! long name description
      logical, intent(in) :: history_aerosol           ! add fields to default history output

      integer :: l, nfld
      character(len=len(cnames)), allocatable :: names(:)

      allocate(fld_ndx(size(mask)))
      allocate(names(count(mask)))

      nfld = 0
      do l = 1, size(mask)
         if (mask(l)) then
            nfld = nfld + 1
            names(nfld) = cnames(l)
            fld_ndx(l) = nfld
         else
            fld_ndx(l) = 0
         end if
      end do

      call diag%init( names, suffix, desc, history_aerosol )

   end subroutine init_col_tend

   ! sets the numbers of modes and species, the constituent indices of mode
   ! number and species, the primary carbon and accumulation modes, and the
   ! cloud-borne constituent names from the aerosol properties object
   subroutine set_mode_tables(aero_props)

      use aerosol_properties_mod, only: aerosol_properties, aero_name_len
      use constituents,   only: pcnst, cnst_get_ind
      use cam_abortutils, only: endrun
      use modal_aero_data, only: ntot_amode_mad => ntot_amode, nspec_amode_mad => nspec_amode, &
                                 nspec_max_mad => nspec_max, numptr_amode_mad => numptr_amode, &
                                 lmassptr_amode_mad => lmassptr_amode, alnsg_amode_mad => alnsg_amode, &
                                 modeptr_pcarbon_mad => modeptr_pcarbon, modeptr_accum_mad => modeptr_accum, &
                                 specmw_amode_mad => specmw_amode, specdens_amode_mad => specdens_amode, &
                                 spechygro_mad => spechygro, cnst_name_cw_mad => cnst_name_cw

      class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object

      integer :: m, l, idx
      real(r8) :: spec_mw, spec_dens, spec_hygro
      character(len=aero_name_len) :: name_a, name_c, mode_type
      character(len=*), parameter :: subname = 'modal_aero_gasaerexch_cam::set_mode_tables'

      nmodes = aero_props%nbins()

      allocate(nspec_mode(nmodes))
      do m = 1, nmodes
         nspec_mode(m) = aero_props%nspecies(m)
      end do
      nspec_max = maxval(nspec_mode)

      allocate(num_a_ndx(nmodes))
      allocate(mass_a_ndx(nspec_max,nmodes))
      allocate(cw_cnst_name(pcnst))
      mass_a_ndx(:,:) = -1
      cw_cnst_name(:) = ' '
      mode_pcarbon = -1
      mode_accum = -1

      do m = 1, nmodes
         call aero_props%num_names(m, name_a, name_c)
         call cnst_get_ind(name_a, idx, abort=.false.)
         num_a_ndx(m) = idx
         if (idx > 0) cw_cnst_name(idx) = name_c
         do l = 1, nspec_mode(m)
            call aero_props%mmr_names(m, l, name_a, name_c)
            call cnst_get_ind(name_a, idx, abort=.false.)
            mass_a_ndx(l,m) = idx
            if (idx > 0) cw_cnst_name(idx) = name_c
         end do
         ! the last mode of each type is used
         mode_type = aero_props%bin_name(m)
         if (trim(mode_type) == 'accum') mode_accum = m
         if (trim(mode_type) == 'primary_carbon') mode_pcarbon = m
      end do

      ! check consistency with modal_aero_data
      if (nmodes /= ntot_amode_mad .or. nspec_max /= nspec_max_mad) then
         call endrun(subname//': inconsistent number of modes or maximum number of species')
      end if
      if (.not. same_mode(mode_accum, modeptr_accum_mad) .or. &
          .not. same_mode(mode_pcarbon, modeptr_pcarbon_mad)) then
         call endrun(subname//': inconsistent accumulation or primary carbon mode')
      end if
      do m = 1, nmodes
         if (nspec_mode(m) /= nspec_amode_mad(m) .or. num_a_ndx(m) /= numptr_amode_mad(m)) then
            call endrun(subname//': inconsistent number of species or number index')
         end if
         if (aero_props%alogsig(m) /= alnsg_amode_mad(m)) then
            call endrun(subname//': inconsistent log of geometric standard deviation')
         end if
         do l = 1, nspec_mode(m)
            if (mass_a_ndx(l,m) /= lmassptr_amode_mad(l,m)) then
               call endrun(subname//': inconsistent species index')
            end if
            call aero_props%get(m, l, spec_mw=spec_mw, density=spec_dens, hygro=spec_hygro)
            if (spec_mw /= specmw_amode_mad(l,m) .or. spec_dens /= specdens_amode_mad(l,m) .or. &
                spec_hygro /= spechygro_mad(l,m)) then
               call endrun(subname//': inconsistent species properties')
            end if
         end do
      end do
      do idx = 1, pcnst
         if (trim(cw_cnst_name(idx)) /= trim(cnst_name_cw_mad(idx))) then
            call endrun(subname//': inconsistent cloud-borne name of '//trim(cnst_name_cw_mad(idx)))
         end if
      end do

   contains

      ! true if two mode indices are the same mode, or both indicate no mode (<= 0)
      pure logical function same_mode(m1, m2)
         integer, intent(in) :: m1
         integer, intent(in) :: m2
         same_mode = (m1 == m2) .or. (m1 <= 0 .and. m2 <= 0)
      end function same_mode

   end subroutine set_mode_tables

   ! sets constituent indices of the interstitial aerosol species and soa gas
   ! species involved in gas-aerosol exchange, and the soa/poa molecular
   ! weights, from the aerosol properties object
   subroutine set_species_ndx(aero_props)

      use aerosol_properties_mod, only: aerosol_properties, aero_name_len
      use aero_gasaerexch, only: aero_gasaerexch_spec_ndx
      use constituents,   only: pcnst, cnst_get_ind
      use cam_abortutils, only: endrun
      use modal_aero_data, only: lptr_so4_a_amode, lptr_nh4_a_amode, &
                                 lptr2_soa_a_amode, lptr2_pom_a_amode, lptr2_soa_g_amode, &
                                 nsoa_mad => nsoa, npoa_mad => npoa

      class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object

      integer :: m, j
      integer :: nspec_typ(nmodes)
      integer, allocatable :: spec_ndx(:,:)
      character(len=128) :: msg
      character(len=*), parameter :: subname = 'modal_aero_gasaerexch_cam::set_species_ndx'

      ! soa and poa species, soa gas species and soa/poa molecular weights
      call soa_exch%init(aero_props)
      nsoa = soa_exch%nsoa
      npoa = soa_exch%npoa

      allocate(so4_a_ndx(nmodes))
      allocate(nh4_a_ndx(nmodes))
      allocate(soa_a_ndx(nmodes,nsoa))
      allocate(pom_a_ndx(nmodes,npoa))

      so4_a_ndx(:) = -1
      nh4_a_ndx(:) = -1
      soa_a_ndx(:,:) = -1
      pom_a_ndx(:,:) = -1

      ! sulfate and ammonium -- the last species of the type in a mode is used
      call aero_gasaerexch_spec_ndx(aero_props, 'sulfate', nspec_typ, spec_ndx)
      do m = 1, nmodes
         if (nspec_typ(m) > 0) so4_a_ndx(m) = cnst_ndx(m, spec_ndx(m,nspec_typ(m)))
      end do

      call aero_gasaerexch_spec_ndx(aero_props, 'ammonium', nspec_typ, spec_ndx)
      do m = 1, nmodes
         if (nspec_typ(m) > 0) nh4_a_ndx(m) = cnst_ndx(m, spec_ndx(m,nspec_typ(m)))
      end do

      do m = 1, nmodes
         do j = 1, soa_exch%nsoa_bin(m)
            soa_a_ndx(m,j) = cnst_ndx(m, soa_exch%soa_spec(m,j))
         end do
         do j = 1, soa_exch%npoa_bin(m)
            pom_a_ndx(m,j) = cnst_ndx(m, soa_exch%poa_spec(m,j))
         end do
      end do

      ! check consistency with the modal_aero_data species counts and pointers
      if (nsoa /= nsoa_mad .or. npoa /= npoa_mad) then
         write(msg,'(a,4i6)') subname//': inconsistent numbers of soa/poa species ', nsoa, nsoa_mad, npoa, npoa_mad
         call endrun(trim(msg))
      end if
      do m = 1, nmodes
         call check_ndx(so4_a_ndx(m), lptr_so4_a_amode(m), 'sulfate')
         call check_ndx(nh4_a_ndx(m), lptr_nh4_a_amode(m), 'ammonium')
         do j = 1, nsoa
            call check_ndx(soa_a_ndx(m,j), lptr2_soa_a_amode(m,j), 's-organic')
         end do
         do j = 1, npoa
            call check_ndx(pom_a_ndx(m,j), lptr2_pom_a_amode(m,j), 'p-organic')
         end do
      end do
      do j = 1, nsoa
         if (soa_exch%soag_cnst(j) /= lptr2_soa_g_amode(j)) then
            write(msg,'(a,3i8)') subname//': inconsistent soa gas index ', j, soa_exch%soag_cnst(j), &
                 lptr2_soa_g_amode(j)
            call endrun(trim(msg))
         end if
      end do

   contains

      ! constituent index of the interstitial species l of mode mm (<= 0 if not a constituent)
      ! (not pure: queries the aerosol properties object and the constituents module)
      integer function cnst_ndx(mm, l)
         integer, intent(in) :: mm  ! mode index
         integer, intent(in) :: l   ! species index within the mode

         character(len=aero_name_len) :: name_a, name_c

         call aero_props%mmr_names(mm, l, name_a, name_c)
         call cnst_get_ind(name_a, cnst_ndx, abort=.false.)

      end function cnst_ndx

      subroutine check_ndx(ndx, lptr, spectype)
         integer, intent(in) :: ndx              ! index from aerosol properties
         integer, intent(in) :: lptr             ! index from modal_aero_data
         character(len=*), intent(in) :: spectype ! species type

         if ((lptr > 0) .and. (lptr <= pcnst)) then
            if (ndx == lptr) return
         else
            if (ndx <= 0) return
         end if
         write(msg,'(a,2i8)') subname//': inconsistent '//trim(spectype)//' index for mode ', m, ndx
         call endrun(trim(msg))

      end subroutine check_ndx

   end subroutine set_species_ndx

end module modal_aero_gasaerexch_cam
