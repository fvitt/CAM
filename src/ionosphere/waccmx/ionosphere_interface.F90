module ionosphere_interface

   use shr_kind_mod,        only: r8 => shr_kind_r8, cl=>shr_kind_cl
   use cam_abortutils,      only: endrun
   use ppgrid,              only: begchunk, endchunk, pcols, pver
   use phys_grid,           only: get_ncols_p

   use dpie_coupling,       only: d_pie_init
   use dpie_coupling,       only: d_pie_epotent
   use dpie_coupling,       only: d_pie_coupling         ! WACCM-X ionosphere/electrodynamics coupling
   use short_lived_species, only: slvd_index, slvd_pbf_ndx => pbf_idx ! Routines to access short lived species

   use chem_mods,           only: adv_mass      ! Array holding mass values for short lived species
   use mo_chem_utls,        only: get_spc_ndx   ! Routine to get index of adv_mass array for short lived species
   use physics_buffer,      only: pbuf_get_chunk, pbuf_get_field
   use physics_buffer,      only: pbuf_get_index

   use constituents,        only: cnst_get_ind, cnst_mw
   use physconst,           only: rga
   use oplus,               only: oplus_init
   use edyn_init,           only: edynamo_init
   use pio,                 only: var_desc_t
   use perf_mod,            only: t_startf, t_stopf
   use epotential_params,   only: epot_active, epot_crit_colats
   use shr_const_mod,  only: SHR_CONST_REARTH ! meters

   implicit none

   private

   public :: ionosphere_readnl
   public :: ionosphere_init
   public :: ionosphere_run1
   public :: ionosphere_run2
   public :: ionosphere_init_restart
   public :: ionosphere_write_restart
   public :: ionosphere_read_restart
   public :: ionosphere_final

   ! private data

   ! opmmrtm1_phys is O+ at previous time step (phys grid decomposed)
   ! It needs to persist from time-step to time-step and across restarts
   ! On physics grid
   real(r8), allocatable :: Opmmrtm1_phys(:,:,:)
   type(var_desc_t)      :: Optm1_vdesc
   logical :: Opmmrtm1_initialized = .false.

   real(r8), allocatable :: NOpmmrtm1_phys(:,:,:)
   real(r8), allocatable :: O2pmmrtm1_phys(:,:,:)
   real(r8), allocatable :: Fepmmrtm1_phys(:,:,:)
   real(r8), allocatable :: Mgpmmrtm1_phys(:,:,:)
   real(r8), allocatable :: Napmmrtm1_phys(:,:,:)
   real(r8), allocatable :: Capmmrtm1_phys(:,:,:)
   real(r8), allocatable :: Kpmmrtm1_phys(:,:,:)
   real(r8), allocatable :: Sipmmrtm1_phys(:,:,:)
   logical :: NOpmmrtm1_initialized = .false.
   logical :: O2pmmrtm1_initialized = .false.
   logical :: Fepmmrtm1_initialized = .false.
   logical :: Mgpmmrtm1_initialized = .false.
   logical :: Napmmrtm1_initialized = .false.
   logical :: Capmmrtm1_initialized = .false.
   logical :: Kpmmrtm1_initialized = .false.
   logical :: Sipmmrtm1_initialized = .false.
   type(var_desc_t) :: NOptm1_vdesc
   type(var_desc_t) :: O2ptm1_vdesc
   type(var_desc_t) :: Feptm1_vdesc
   type(var_desc_t) :: Mgptm1_vdesc
   type(var_desc_t) :: Naptm1_vdesc
   type(var_desc_t) :: Captm1_vdesc
   type(var_desc_t) :: Kptm1_vdesc
   type(var_desc_t) :: Siptm1_vdesc

   integer :: index_ped, index_hall, index_te, index_ti
   integer :: index_ui, index_vi, index_wi

   integer :: ixo2=-1, ixo=-1, ixh=-1
   integer :: ixo2p=-1, ixnop=-1, ixn2p=-1, ixop=-1

   ! indices for accessing ions in pbuf when non-advected
   integer :: sIndxOp=-1, sIndxO2p=-1, sIndxNOp=-1, sIndxN2p=-1

   ! Metal Ions
   integer :: ixFep=-1  !Jianfei Wu added for Fe+ Na+ and K+
   integer :: ixMgp=-1  !Jianfei Wu added for Fe+ Na+ and K+
   integer :: ixNap=-1  !Jianfei Wu added for Fe+ Na+ and K+
   integer :: ixCap=-1  !Jianfei Wu added for Fe+ Na+ and K+
   integer :: ixKp=-1  !Jianfei Wu added for Fe+ Na+ and K+
   integer :: ixSip=-1  !Jianfei Wu added for Fe+ Na+ and K+

   integer :: sIndxFep=-1 !Jianfei Wu added for Fe+ Na+ and K+
   integer :: sIndxMgp=-1 !Jianfei Wu added for Fe+ Na+ and K+
   integer :: sIndxNap=-1 !Jianfei Wu added for Fe+ Na+ and K+
   integer :: sIndxCap=-1 !Jianfei Wu added for Fe+ Na+ and K+
   integer :: sIndxKp=-1 !Jianfei Wu added for Fe+ Na+ and K+
   integer :: sIndxSip=-1 !Jianfei Wu added for Fe+ Na+ and K+

   real(r8) :: rmassFep = 0._r8  ! Fe+ molecular weight kg/kmol
   real(r8) :: rmassMgp = 0._r8  ! Fe+ molecular weight kg/kmol
   real(r8) :: rmassNap = 0._r8  ! Fe+ molecular weight kg/kmol
   real(r8) :: rmassCap = 0._r8  ! Fe+ molecular weight kg/kmol
   real(r8) :: rmassKp  = 0._r8  ! Fe+ molecular weight kg/kmol
   real(r8) :: rmassSip = 0._r8  ! Fe+ molecular weight kg/kmol

   real(r8) :: rmassO2  = 0._r8  ! O2 molecular weight kg/kmol
   real(r8) :: rmassO1  = 0._r8  ! O atomic weight kg/kmol
   real(r8) :: rmassH   = 0._r8  ! H atomic weight kg/kmol
   real(r8) :: rmassN2  = 0._r8  ! N2 molecular weight kg/kmol
   real(r8) :: rmassO2p = 0._r8  ! O2+ molecular weight kg/kmol
   real(r8) :: rmassNOp = 0._r8  ! NO+ molecular weight kg/kmol
   real(r8) :: rmassN2p = 0._r8  ! N2+ molecular weight kg/kmol
   real(r8) :: rmassOp  = 0._r8  ! O+ molecular weight kg/kmol

   ! ionos_edyn_active == .true. will activate the edynamo which will
   !   generate ion drift velocities used in oplus transport, otherwise
   !   empirical ion drifts calculated in exbdrift (physics) will be used.
   logical, public,  protected :: ionos_edyn_active = .true.
   logical,          protected :: ionos_xport_active = .true.  ! if true, call d_pie_coupling
   !
   logical, public,  protected :: ionos_oplus_xport = .true.    ! if true, call sub oplus (based on tiegcm oplus.F)
   integer, public,  protected :: ionos_xport_nsplit = 5        ! number of substeps for O+ transport per model time step
   logical, public,  protected :: oplus_ring_polar_filter = .false. ! switch to apply ring polar filter
   integer, public,  protected :: vxb_xport_nsplit = 15        ! number of substeps for O+ transport per model time step

   real(r8) :: oplus_adiff_limiter = 1.5e+8_r8  ! limiter for ambipolar diffusion coefficient
   real(r8) :: oplus_shapiro_const = 0.03_r8    ! shapiro constant for spatial smoother
   logical  :: oplus_enforce_floor = .true.     ! switch to apply Stan's  floor

   integer, parameter :: max_num_files = 20
   character(len=cl) :: wei05_coefs_file = 'NONE' !'wei05sc.nc'
   character(len=cl) :: amienh_files(max_num_files) = 'NONE'
   character(len=cl) :: amiesh_files(max_num_files) = 'NONE'
   character(len=cl) :: ltr_files(max_num_files) = 'NONE'


   character(len=16) :: ionos_epotential_model = 'none'
   logical           :: ionos_epotential_amie = .false.
   logical           :: ionos_epotential_ltr = .false.
   integer           :: indxefx=-1, indxkev=-1

   integer           :: oplus_nlon, oplus_nlat   ! Oplus grid
   integer           :: ionos_npes = -1

   logical :: state_debug_checks = .false.
   logical :: ionos_debug_hist = .false.

   integer :: mag_nlon=0, mag_nlat=0, mag_nlev=0, mag_ngrid=0

   real(r8), parameter :: rearth_inv = 1._r8/SHR_CONST_REARTH ! /meters

 contains

   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine ionosphere_readnl( nlfile )

      use namelist_utils, only: find_group_name
      use units,          only: getunit, freeunit
      use spmd_utils,     only: masterproc, mpicom, masterprocid
      use spmd_utils,     only: mpi_real8, mpi_logical, mpi_integer, mpi_character
      use cam_logfile,    only: iulog

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr, ipos
      integer :: oplus_grid(2)
      character(len=8) :: edyn_grid
      integer :: total_pes
      character(len=*), parameter :: subname = 'ionosphere_readnl'

      namelist /ionosphere_nl/ ionos_xport_active, ionos_edyn_active, ionos_oplus_xport, ionos_xport_nsplit, vxb_xport_nsplit
      namelist /ionosphere_nl/ oplus_adiff_limiter, oplus_shapiro_const, oplus_enforce_floor, oplus_ring_polar_filter
      namelist /ionosphere_nl/ ionos_epotential_model, ionos_epotential_amie, ionos_epotential_ltr, wei05_coefs_file
      namelist /ionosphere_nl/ amienh_files, amiesh_files, wei05_coefs_file, ltr_files
      namelist /ionosphere_nl/ epot_crit_colats
      namelist /ionosphere_nl/ ionos_npes
      namelist /ionosphere_nl/ oplus_grid, edyn_grid
      namelist /ionosphere_nl/ ionos_debug_hist

      oplus_grid = 0

      ! Read namelist
      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'ionosphere_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, ionosphere_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

      ! Broadcast namelist variables
      call mpi_bcast(ionos_xport_active,  1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(ionos_edyn_active,   1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(ionos_oplus_xport,   1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(ionos_xport_nsplit,  1, mpi_integer, masterprocid, mpicom, ierr)
      call mpi_bcast(vxb_xport_nsplit,    1, mpi_integer, masterprocid, mpicom, ierr)
      call mpi_bcast(oplus_adiff_limiter, 1, mpi_real8,   masterprocid, mpicom, ierr)
      call mpi_bcast(ionos_epotential_model, len(ionos_epotential_model), mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(ionos_epotential_amie,1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(ionos_epotential_ltr,1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(wei05_coefs_file, len(wei05_coefs_file), mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(amienh_files, max_num_files*len(amienh_files(1)), mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(amiesh_files, max_num_files*len(amiesh_files(1)), mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(ltr_files, max_num_files*len(ltr_files(1)), mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(oplus_shapiro_const, 1, mpi_real8,   masterprocid, mpicom, ierr)
      call mpi_bcast(oplus_enforce_floor, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(oplus_ring_polar_filter,1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(epot_crit_colats,    2, mpi_real8,   masterprocid, mpicom, ierr)
      call mpi_bcast(ionos_npes,          1, mpi_integer, masterprocid, mpicom, ierr)
      call mpi_bcast(oplus_grid,          2, mpi_integer, masterprocid, mpicom, ierr)
      call mpi_bcast(edyn_grid,           8, mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(ionos_debug_hist,    1, mpi_logical, masterprocid, mpicom, ierr)

      ! Extract grid settings
      oplus_nlon = oplus_grid(1)
      oplus_nlat = oplus_grid(2)

      ipos = scan(edyn_grid,'x')
      read(edyn_grid(:ipos-1),*) mag_nlon
      read(edyn_grid(ipos+1:),*) mag_nlat

      mag_nlev = 5 + int(log(real(mag_nlon,r8)/80._r8)/log(2._r8))
      mag_ngrid = (mag_nlon/10)*2

      ! Set npes in case of default settings
      call mpi_comm_size(mpicom, total_pes, ierr)
      if (ionos_npes<1) then
         ionos_npes = total_pes
      else if (ionos_npes>total_pes) then
         call endrun('ionosphere_readnl: ionos_npes > total_pes')
      end if

      ! log the user settings
      if (masterproc) then
         write(iulog,*) 'ionosphere_readnl: ionos_xport_active     = ', ionos_xport_active
         write(iulog,*) 'ionosphere_readnl: ionos_edyn_active      = ', ionos_edyn_active
         write(iulog,*) 'ionosphere_readnl: ionos_oplus_xport      = ', ionos_oplus_xport
         write(iulog,*) 'ionosphere_readnl: ionos_xport_nsplit     = ', ionos_xport_nsplit
         write(iulog,*) 'ionosphere_readnl: vxb_xport_nsplit       = ', vxb_xport_nsplit
         write(iulog,*) 'ionosphere_readnl: ionos_epotential_model = ', trim(ionos_epotential_model)
         write(iulog,*) 'ionosphere_readnl: ionos_epotential_amie  = ', ionos_epotential_amie
         write(iulog,*) 'ionosphere_readnl: ionos_epotential_ltr   = ', ionos_epotential_ltr
         write(iulog,'(a,2(g12.4))') &
                        'ionosphere_readnl: epot_crit_colats       = ', epot_crit_colats
         write(iulog,'(a,i0)') 'ionosphere_readnl: ionos_npes = ',ionos_npes
         write(iulog,*) 'ionosphere_readnl: oplus_adiff_limiter    = ', oplus_adiff_limiter
         write(iulog,*) 'ionosphere_readnl: oplus_shapiro_const    = ', oplus_shapiro_const
         write(iulog,*) 'ionosphere_readnl: oplus_enforce_floor    = ', oplus_enforce_floor
         write(iulog,*) 'ionosphere_readnl: oplus_ring_polar_filter= ', oplus_ring_polar_filter
         if (ionos_xport_active) then
            write(iulog,'(a,i0)') 'ionosphere_readnl: oplus_nlon = ',oplus_nlon
            write(iulog,'(a,i0)') 'ionosphere_readnl: oplus_nlat = ',oplus_nlat
            write(iulog,'(a,i0)') 'ionosphere_readnl: edyn_grid = '//edyn_grid
            write(iulog,'(a,i0)') 'ionosphere_readnl: mag_nlon = ',mag_nlon
            write(iulog,'(a,i0)') 'ionosphere_readnl: mag_nlat = ',mag_nlat
            write(iulog,'(a,i0)') 'ionosphere_readnl: mag_nlev = ',mag_nlev
            write(iulog,'(a,i0)') 'ionosphere_readnl: mag_ngrid = ',mag_ngrid
         end if
      end if
      epot_active = .true.

   end subroutine ionosphere_readnl

   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine ionosphere_init()
      use spmd_utils,      only: mpicom, iam
      use physics_buffer,  only: pbuf_add_field, dtype_r8
      use cam_control_mod, only: initial_run
      use cam_history,     only: addfld, add_default, horiz_only
      use edyn_mpi,        only: mp_init
      use edyn_geogrid,    only: set_geogrid
      use edyn_maggrid,    only: alloc_maggrid
      use mo_apex,         only: mo_apex_init1
      ! Hybrid level definitions:
      use ref_pres,        only: pref_mid  ! target alev(pver) midpoint levels
      use ref_pres,        only: pref_edge ! target ailev(pverp) interface levels
      use amie_module,     only: init_amie
      use ltr_module,      only: init_ltr
      use wei05sc,         only: weimer05_init
      use phys_control,    only: phys_getopts
      use vxb, only: vxb_init

      ! local variables:
      integer :: sIndx
      character(len=*), parameter :: subname = 'ionosphere_init'

      call phys_getopts(state_debug_checks_out=state_debug_checks)

      if ( ionos_epotential_amie .or. ionos_epotential_ltr) then
         call pbuf_add_field('AUREFX', 'global', dtype_r8, (/pcols/), indxefx)  ! Prescribed Energy flux
         call pbuf_add_field('AURKEV', 'global', dtype_r8, (/pcols/), indxkev)  ! Prescribed Mean energy
      end if
      if (initial_run) then
         ! Read initial conditions (O+) on physics grid
         call ionosphere_read_ic()
      end if

      op_transport: if (ionos_xport_active) then

         index_ped  = pbuf_get_index('PedConduct')
         index_hall = pbuf_get_index('HallConduct')

         index_te   = pbuf_get_index('TElec')
         index_ti   = pbuf_get_index('TIon')
         !
         ! pbuf indices to empirical ion drifts, to be passed to oplus_xport,
         ! if ionos_edyn_active is false.
         !
         index_ui   = pbuf_get_index('UI')
         index_vi   = pbuf_get_index('VI')
         index_wi   = pbuf_get_index('WI')

         !---------------------------------------------------------------------
         ! Get indices for neutrals to get mixing ratios from state%q and masses
         !---------------------------------------------------------------------
         call cnst_get_ind('O2' ,ixo2 )
         call cnst_get_ind('O'  ,ixo )
         call cnst_get_ind('H'  ,ixh )
         !------------------------------------
         ! Get neutral molecular weights
         !------------------------------------
         rmassO2 = cnst_mw(ixo2)
         rmassO1 = cnst_mw(ixo)
         rmassH  = cnst_mw(ixh)
         rmassN2 = 28._r8

         ! set ion indexes and molecular masses
         call set_indices_mass('Op', ixop, sIndxOp, rmassOp)
         call set_indices_mass('O2p', ixo2p, sIndxO2p, rmassO2p)
         call set_indices_mass('NOp', ixnop, sIndxNOp, rmassNOp)
         call set_indices_mass('N2p', ixn2p, sIndxN2p, rmassN2p)

         call alloc_maggrid( mag_nlon, mag_nlat, mag_nlev, mag_ngrid )

         call mp_init(mpicom, ionos_npes, oplus_nlon, oplus_nlat, pver) ! set ntask,mytid

         ! set global geographic grid (sets coordinate distribution)
         ! lon0, lon1, etc. are set here
         call set_geogrid(oplus_nlon, oplus_nlat, pver, ionos_npes, iam, pref_mid, pref_edge)

         call edynamo_init(mpicom, ionos_debug_hist)

         call d_pie_init(ionos_edyn_active, ionos_oplus_xport, ionos_xport_nsplit,vxb_xport_nsplit, epot_crit_colats, &
                         ionos_debug_hist)

         call ionosphere_alloc()

         call oplus_init(oplus_adiff_limiter, oplus_shapiro_const, oplus_enforce_floor, &
                         oplus_ring_polar_filter, ionos_debug_hist)

         call vxb_init( oplus_adiff_limiter, oplus_shapiro_const, oplus_enforce_floor ) !Jianfei Wu

         call addfld('OpTM1&IC', (/ 'lev' /), 'I', 'kg/kg', 'O+ at time step minus 1', gridname='physgrid')
         call add_default ('OpTM1&IC',0, 'I')

         if (ixFep>0 .or. sIndxFep>0) then
            call addfld('NOpTM1&IC', (/ 'lev' /), 'I', 'kg/kg', 'NO+ at time step minus 1', gridname='physgrid')
            call add_default ('NOpTM1&IC',0, 'I')
            call addfld('O2pTM1&IC', (/ 'lev' /), 'I', 'kg/kg', 'O2+ at time step minus 1', gridname='physgrid')
            call add_default ('O2pTM1&IC',0, 'I')
            call addfld('FepTM1&IC', (/ 'lev' /), 'I', 'kg/kg', 'Fe+ at time step minus 1', gridname='physgrid')
            call add_default ('FepTM1&IC',0, 'I')
         end if
         if (ixMgp>0 .or. sIndxMGp>0) then
            call addfld('MgpTM1&IC', (/ 'lev' /), 'I', 'kg/kg', 'Mg+ at time step minus 1', gridname='physgrid')
            call add_default ('MgpTM1&IC',0, 'I')
         end if
         if (ixNap>0 .or. sIndxNap>0) then
            call addfld('NapTM1&IC', (/ 'lev' /), 'I', 'kg/kg', 'Na+ at time step minus 1', gridname='physgrid')
            call add_default ('NapTM1&IC',0, 'I')
         end if
         if (ixCap>0 .or. sIndxCap>0) then
            call addfld('CapTM1&IC', (/ 'lev' /), 'I', 'kg/kg', 'Ca+ at time step minus 1', gridname='physgrid')
            call add_default ('CapTM1&IC',0, 'I')
         end if
         if (ixKp>0 .or. sIndxKp>0) then
            call addfld('KpTM1&IC', (/ 'lev' /), 'I', 'kg/kg', 'K+ at time step minus 1', gridname='physgrid')
            call add_default ('KpTM1&IC',0, 'I')
         end if
         if (ixSip>0 .or. sIndxSip>0) then
            call addfld('SipTM1&IC', (/ 'lev' /), 'I', 'kg/kg', 'Si+ at time step minus 1', gridname='physgrid')
            call add_default ('SipTM1&IC',0, 'I')
         end if

         call addfld('Fep_phys0', (/ 'lev' /), 'I', 'kg/kg', 'Fep before d_pie_cpl ', gridname='physgrid')
         call addfld('Fep_phys1', (/ 'lev' /), 'I', 'kg/kg', 'Fep after d_pie_cpl ', gridname='physgrid')
         call addfld('Mgp_phys0', (/ 'lev' /), 'I', 'kg/kg', 'Mgp before d_pie_cpl ', gridname='physgrid')
         call addfld('Mgp_phys1', (/ 'lev' /), 'I', 'kg/kg', 'Mgp after d_pie_cpl ', gridname='physgrid')

      end if op_transport

      ! This has to be after edynamo_init (where maggrid is initialized)
      call mo_apex_init1()

      if (ionos_edyn_active) then
         call addfld ('UI',(/ 'lev' /),'I','m/s', 'UI Zonal ion drift from edynamo')
         call addfld ('VI',(/ 'lev' /),'I','m/s', 'VI Meridional ion drift from edynamo')
         call addfld ('WI',(/ 'lev' /),'I','m/s', 'WI Vertical ion drift from edynamo')
         call addfld ('UI&IC', (/ 'lev' /), 'I','m/s', 'Zonal ion drift velocity')
         call addfld ('VI&IC', (/ 'lev' /), 'I','m/s', 'Meridional ion drift velocity')
         call addfld ('WI&IC', (/ 'lev' /), 'I','m/s', 'Vertical ion drift velocity')
         call add_default ('UI&IC', 0, ' ')
         call add_default ('VI&IC', 0, ' ')
         call add_default ('WI&IC', 0, ' ')
      end if
      if ( ionos_epotential_amie ) then
         call init_amie(amienh_files,amiesh_files)
         call addfld ('amie_efx_phys', horiz_only, 'I', 'mW/m2', 'AMIE energy flux')
         call addfld ('amie_kev_phys', horiz_only, 'I', 'keV',   'AMIE mean energy')
      end if
      if ( ionos_epotential_ltr ) then
         call init_ltr(ltr_files)
         call addfld ('ltr_efx_phys', horiz_only, 'I', 'mW/m2', 'LTR energy flux')
         call addfld ('ltr_kev_phys', horiz_only, 'I', 'keV',  'LTR mean energy')
      end if
      if ( trim(ionos_epotential_model) == 'weimer' ) then
         call weimer05_init(wei05_coefs_file)
      end if

      ! d_pie_coupling diagnostics
      call addfld ('Z3GM',       (/ 'lev' /), 'I', 'm',                       &
           'Geometric height', gridname='physgrid')
      call addfld ('Z3GMI',      (/ 'lev' /), 'I', 'm',                       &
           'Geometric height (Interfaces)', gridname='physgrid')

   end subroutine ionosphere_init

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
   subroutine ionosphere_run1(pbuf2d)
      use physics_buffer, only: physics_buffer_desc
      use cam_history,    only: outfld, write_inithist

      ! args
      type(physics_buffer_desc), pointer :: pbuf2d(:,:)

      ! local vars
      integer :: i, j, lchnk,  blksize ! indices
      type(physics_buffer_desc), pointer :: pbuf_chnk(:)

      real(r8), pointer :: pbuf_efx(:) ! Pointer to prescribed energy flux in pbuf
      real(r8), pointer :: pbuf_kev(:) ! Pointer to prescribed mean energy in pbuf

      integer :: ncol
      real(r8), pointer :: prescr_efx(:) ! prescribed energy flux
      real(r8), pointer :: prescr_kev(:) ! prescribed characteristic mean energy

      if( write_inithist() .and. ionos_xport_active ) then
         do lchnk = begchunk, endchunk
            call outfld ('OpTM1&IC', opmmrtm1_phys(:,:,lchnk), pcols, lchnk)
         end do
         if (allocated(NOpmmrtm1_phys)) then
            do lchnk = begchunk, endchunk
               call outfld ('NOpTM1&IC', NOpmmrtm1_phys(:,:,lchnk), pcols, lchnk)
            end do
         end if
         if (allocated(O2pmmrtm1_phys)) then
            do lchnk = begchunk, endchunk
               call outfld ('O2pTM1&IC', O2pmmrtm1_phys(:,:,lchnk), pcols, lchnk)
            end do
         end if
         if (allocated(Fepmmrtm1_phys)) then
            do lchnk = begchunk, endchunk
               call outfld ('FepTM1&IC', Fepmmrtm1_phys(:,:,lchnk), pcols, lchnk)
            end do
         end if
         if (allocated(Mgpmmrtm1_phys)) then
            do lchnk = begchunk, endchunk
               call outfld ('MgpTM1&IC', Mgpmmrtm1_phys(:,:,lchnk), pcols, lchnk)
            end do
         end if
         if (allocated(Napmmrtm1_phys)) then
            do lchnk = begchunk, endchunk
               call outfld ('NapTM1&IC', Napmmrtm1_phys(:,:,lchnk), pcols, lchnk)
            end do
         end if
         if (allocated(Capmmrtm1_phys)) then
            do lchnk = begchunk, endchunk
               call outfld ('CapTM1&IC', Capmmrtm1_phys(:,:,lchnk), pcols, lchnk)
            end do
         end if
         if (allocated(Kpmmrtm1_phys)) then
            do lchnk = begchunk, endchunk
               call outfld ('KpTM1&IC', Kpmmrtm1_phys(:,:,lchnk), pcols, lchnk)
            end do
         end if
         if (allocated(Sipmmrtm1_phys)) then
            do lchnk = begchunk, endchunk
               call outfld ('SipTM1&IC', Sipmmrtm1_phys(:,:,lchnk), pcols, lchnk)
            end do
         end if

      end if

      nullify(prescr_efx)
      nullify(prescr_kev)
      prescribed_epot: if ( ionos_epotential_amie .or. ionos_epotential_ltr ) then
         blksize = 0
         do lchnk = begchunk, endchunk
            blksize = blksize + get_ncols_p(lchnk)
         end do

         allocate(prescr_efx(blksize))
         allocate(prescr_kev(blksize))

         ! data assimilated potential
         call d_pie_epotent(ionos_epotential_model, epot_crit_colats, &
              cols=1, cole=blksize, efx_phys=prescr_efx, kev_phys=prescr_kev, &
              amie_in=ionos_epotential_amie, ltr_in=ionos_epotential_ltr )

         ! transform to pbuf for aurora...

         j = 0
         chnk_loop1: do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
            call pbuf_get_field(pbuf_chnk, indxefx, pbuf_efx)
            call pbuf_get_field(pbuf_chnk, indxkev, pbuf_kev)

            do i = 1, ncol
               j = j + 1
               pbuf_efx(i) = prescr_efx(j)
               pbuf_kev(i) = prescr_kev(j)
            end do

            if ( ionos_epotential_amie ) then
               call outfld('amie_efx_phys', pbuf_efx, pcols, lchnk)
               call outfld('amie_kev_phys', pbuf_kev, pcols, lchnk)
            endif
            if ( ionos_epotential_ltr) then
               call outfld('ltr_efx_phys', pbuf_efx, pcols, lchnk )
               call outfld('ltr_kev_phys', pbuf_kev, pcols, lchnk )
            end if
         end do chnk_loop1

         deallocate(prescr_efx, prescr_kev)
         nullify(prescr_efx)
         nullify(prescr_kev)

      else

         ! set cross tail potential before physics --
         !   aurora uses weimer derived potential
         call d_pie_epotent( ionos_epotential_model, epot_crit_colats )

      end if prescribed_epot

   end subroutine ionosphere_run1

   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine ionosphere_run2(phys_state, pbuf2d)

      use physics_types,  only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use cam_history,    only: outfld, write_inithist, hist_fld_active
      use shr_assert_mod, only: shr_assert_in_domain

      ! - pull some fields from pbuf and dyn_in
      ! - invoke ionosphere/electro-dynamics coupling
      ! - push some fields back to physics via pbuf...

      ! args
      type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
      type(physics_buffer_desc), pointer       :: pbuf2d(:,:)

      ! local vars
      integer :: i,j,k, lchnk
      integer :: astat

      type(physics_buffer_desc), pointer :: pbuf_chnk(:)

      real(r8), pointer :: sigma_ped_phys(:,:) ! Pedersen Conductivity from pbuf
      real(r8), pointer :: sigma_hall_phys(:,:) ! Hall Conductivity from pbuf
      real(r8), pointer :: te_phys(:,:)      ! te from pbuf
      real(r8), pointer :: ti_phys(:,:)      ! ti from pbuf
      real(r8), pointer :: mmrPO2p_phys(:,:)=>null() ! O2+ from pbuf
      real(r8), pointer :: mmrPNOp_phys(:,:)=>null() ! NO+ from pbuf
      real(r8), pointer :: mmrPN2p_phys(:,:)=>null() ! N2+ from pbuf
      real(r8), pointer :: mmrPOp_phys(:,:) =>null() ! O+ from pbuf
      real(r8), pointer :: mmrPFep_phys(:,:)=>null() ! Fe+ from pbuf
      real(r8), pointer :: mmrPMgp_phys(:,:)=>null() ! Mg+ from pbuf
      real(r8), pointer :: mmrPNap_phys(:,:)=>null() ! Na+ from pbuf
      real(r8), pointer :: mmrPCap_phys(:,:)=>null() ! Ca+ from pbuf
      real(r8), pointer :: mmrPKp_phys(:,:) =>null() ! K+ from pbuf
      real(r8), pointer :: mmrPSip_phys(:,:)=>null() ! Si+ from pbuf
      !
      ! Empirical ion drifts from exbdrift (to be converted to blocked for dpie_coupling):
      real(r8), pointer :: ui_phys(:,:)       ! zonal ion drift from pbuf
      real(r8), pointer :: vi_phys(:,:)       ! meridional ion drift from pbuf
      real(r8), pointer :: wi_phys(:,:)       ! vertical ion drift from pbuf

      integer           :: ncol

      integer           :: blksize ! number of columns in 2D block

      real(r8), pointer :: sigma_ped_blck (:,:)
      real(r8), pointer :: sigma_hall_blck(:,:)
      real(r8), pointer :: ti_blck(:,:)
      real(r8), pointer :: te_blck(:,:)
      real(r8), pointer :: zi_blck(:,:) ! Geopotential on interfaces
      real(r8), pointer :: hi_blck(:,:) ! Geometric height on interfaces
      real(r8), pointer :: ui_blck(:,:)
      real(r8), pointer :: vi_blck(:,:)
      real(r8), pointer :: wi_blck(:,:)
      real(r8), pointer :: omega_blck(:,:)
      real(r8), pointer :: tn_blck(:,:)

      ! From physics state
      real(r8), pointer :: u_blck(:,:)
      real(r8), pointer :: v_blck(:,:)
      real(r8), pointer :: pmid_blck(:,:)
      real(r8), pointer :: phis(:)            ! surface geopotential
      ! Constituents
      real(r8), pointer :: n2mmr_blck(:,:)
      real(r8), pointer :: o2mmr_blck(:,:)
      real(r8), pointer :: o1mmr_blck(:,:)
      real(r8), pointer :: h1mmr_blck(:,:)
      real(r8), pointer :: o2pmmr_blck(:,:)=>null() ! O2+ (blocks)
      real(r8), pointer :: nopmmr_blck(:,:)=>null() ! NO+ (blocks)
      real(r8), pointer :: n2pmmr_blck(:,:)=>null() ! N2+ (blocks)
      real(r8), pointer :: opmmr_blck(:,:) =>null() ! O+ (blocks)
      real(r8), pointer :: opmmrtm1_blck(:,:)       ! O+ previous time step (blocks)
      real(r8), pointer :: mbar_blck(:,:)   ! mean molecular weight

      real(r8), pointer :: Fepmmr_blck(:,:)=>null() ! Fe+ (blocks)
      real(r8), pointer :: Mgpmmr_blck(:,:)=>null() ! Mg+ (blocks)
      real(r8), pointer :: Napmmr_blck(:,:)=>null() ! Na+ (blocks)
      real(r8), pointer :: Capmmr_blck(:,:)=>null() ! Ca+ (blocks)
      real(r8), pointer :: Kpmmr_blck(:,:) =>null() ! K+  (blocks)
      real(r8), pointer :: Sipmmr_blck(:,:)=>null() ! Si+ (blocks)

      real(r8), pointer :: O2pmmrtm1_blck(:,:)=>null() ! O2+ previous time step (blocks)
      real(r8), pointer :: NOpmmrtm1_blck(:,:)=>null() ! NO+ previous time step (blocks)
      real(r8), pointer :: Fepmmrtm1_blck(:,:)=>null() ! Fe+ previous time step (blocks)
      real(r8), pointer :: Mgpmmrtm1_blck(:,:)=>null() ! Mg+ previous time step (blocks)
      real(r8), pointer :: Napmmrtm1_blck(:,:)=>null() ! Na+ previous time step (blocks)
      real(r8), pointer :: Capmmrtm1_blck(:,:)=>null() ! Ca+ previous time step (blocks)
      real(r8), pointer :: Kpmmrtm1_blck(:,:) =>null() ! K+  previous time step (blocks)
      real(r8), pointer :: Sipmmrtm1_blck(:,:)=>null() ! Si+ previous time step (blocks)

     ! Temp fields for outfld
      real(r8), pointer :: mass_tmp(:,:)
      real(r8), pointer :: tempm(:,:) => null() ! Temp midpoint field for outfld
      real(r8), pointer :: tempi(:,:) => null() ! Temp interface field for outfld
      real(r8), parameter :: n2min = 1.e-6_r8  ! lower limit of N2 mixing ratios

      character(len=*), parameter :: subname = 'ionosphere_run2'

      ionos_cpl: if (ionos_xport_active) then

         blksize = 0
         do lchnk = begchunk, endchunk
            blksize = blksize + get_ncols_p(lchnk)
         end do

         allocate(phis(pcols), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate phis')
         end if
         allocate(u_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate u_blck')
         end if
         allocate(v_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate v_blck')
         end if
         allocate(sigma_ped_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate sigma_ped_blck')
         end if
         allocate(sigma_hall_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate sigma_hall_blck')
         end if
         allocate(ti_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate ti_blck')
         end if
         allocate(hi_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate hi_blck')
         end if
         allocate(te_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate te_blck')
         end if
         allocate(zi_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate zi_blck')
         end if
         allocate(ui_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate ui_blck')
         end if
         allocate(vi_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate vi_blck')
         end if
         allocate(wi_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate wi_blck')
         end if
         allocate(omega_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate omega_blck')
         end if
         allocate(tn_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate tn_blck')
         end if
         allocate(mass_tmp(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate mass_tmp')
         end if
         allocate(n2mmr_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate n2mmr_blck')
         end if
         allocate(o2mmr_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate o2mmr_blck')
         end if
         allocate(o1mmr_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate o1mmr_blck')
         end if
         allocate(h1mmr_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate h1mmr_blck')
         end if
         allocate(mbar_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate mbar_blck')
         end if
         allocate(pmid_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate pmid_blck')
         end if

         allocate(opmmrtm1_blck(pver, blksize), stat=astat)
         if (astat /= 0) then
            call endrun(subname//': failed to allocate opmmrtm1_blck')
         end if

         if (sIndxOp>0 .or. ixop>0) then
            allocate(opmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate opmmr_blck')
            end if
         end if
         if (sIndxO2p>0 .or. ixO2p>0) then
            allocate(o2pmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate o2pmmr_blck')
            end if
         end if
         if (sIndxNOp>0 .or. IxNOp>0) then
            allocate(nopmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate nopmmr_blck')
            end if
         end if
         if (sIndxN2p>0 .or. ixN2p>0) then
            allocate(n2pmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate n2pmmr_blck')
            end if
         end if

         if (sIndxFep>0 .or. ixFep>0) then
            allocate(Fepmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Fepmmr_blck')
            end if
            allocate(Fepmmrtm1_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Fepmmrtm1_blck')
            end if
            allocate(O2pmmrtm1_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate O2pmmrtm1_blck')
            end if
            allocate(NOpmmrtm1_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate NOpmmrtm1_blck')
            end if
         end if
         if (sIndxMgp>0 .or. ixMgp>0) then
            allocate(Mgpmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Mgpmmr_blck')
            end if
            allocate(Mgpmmrtm1_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Mgpmmrtm1_blck')
            end if
         end if
         if (sIndxNap>0 .or. ixNap>0) then
            allocate(Napmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Napmmr_blck')
            end if
            allocate(Napmmrtm1_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Napmmrtm1_blck')
            end if
         end if
         if (sIndxCap>0 .or. ixCap>0) then
            allocate(Capmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Capmmr_blck')
            end if
            allocate(Capmmrtm1_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Capmmrtm1_blck')
            end if
         end if
         if (sIndxKp>0 .or. ixKp>0) then
            allocate(Kpmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Kpmmr_blck')
            end if
            allocate(Kpmmrtm1_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Kpmmrtm1_blck')
            end if
         end if
         if (sIndxSip>0 .or. ixSip>0) then
            allocate(Sipmmr_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Sipmmr_blck')
            end if
            allocate(Sipmmrtm1_blck(pver, blksize), stat=astat)
            if (astat /= 0) then
               call endrun(subname//': failed to allocate Sipmmrtm1_blck')
            end if
         end if


         if (hist_fld_active('Z3GM')) then
            allocate(tempm(pcols, pver))
         end if

         if (hist_fld_active('Z3GMI')) then
            allocate(tempi(pcols, pver))
         end if

         if (.not.opmmrtm1_initialized) then
            call init_ionstm1(ixOp, sIndxOp, Opmmrtm1_phys)
            Opmmrtm1_initialized=.true.
         endif
         if (allocated(NOpmmrtm1_phys) .and. .not.NOpmmrtm1_initialized) then
            call init_ionstm1(ixNOp, sIndxNOp, NOpmmrtm1_phys)
            NOpmmrtm1_initialized=.true.
         endif
         if (allocated(O2pmmrtm1_phys) .and. .not.O2pmmrtm1_initialized) then
            call init_ionstm1(ixO2p, sIndxO2p, O2pmmrtm1_phys)
            O2pmmrtm1_initialized=.true.
         endif
         if (allocated(Fepmmrtm1_phys) .and. .not.Fepmmrtm1_initialized) then
            call init_ionstm1(ixFep, sIndxFep, Fepmmrtm1_phys)
            Fepmmrtm1_initialized=.true.
         endif
         if (allocated(Mgpmmrtm1_phys) .and. .not.Mgpmmrtm1_initialized) then
            call init_ionstm1(ixMgp, sIndxMgp, Mgpmmrtm1_phys)
            Mgpmmrtm1_initialized=.true.
         endif
         if (allocated(Napmmrtm1_phys) .and. .not.Napmmrtm1_initialized) then
            call init_ionstm1(ixNap, sIndxNap, Napmmrtm1_phys)
            Napmmrtm1_initialized=.true.
         endif
         if (allocated(Capmmrtm1_phys) .and. .not.Capmmrtm1_initialized) then
            call init_ionstm1(ixCap, sIndxCap, Capmmrtm1_phys)
            Capmmrtm1_initialized=.true.
         endif
         if (allocated(Kpmmrtm1_phys) .and. .not.Kpmmrtm1_initialized) then
            call init_ionstm1(ixKp, sIndxKp, Kpmmrtm1_phys)
            Kpmmrtm1_initialized=.true.
         endif
         if (allocated(Sipmmrtm1_phys) .and. .not.Sipmmrtm1_initialized) then
            call init_ionstm1(ixSip, sIndxSip, Sipmmrtm1_phys)
            Sipmmrtm1_initialized=.true.
         endif

         j = 0
         do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

            ! Gather data stored in pbuf and collect into blocked arrays
            ! Get Pedersen and Hall conductivities:
            call pbuf_get_field(pbuf_chnk, index_ped,  sigma_ped_phys)
            call pbuf_get_field(pbuf_chnk, index_hall, sigma_hall_phys)
            ! Get ion and electron temperatures
            call pbuf_get_field(pbuf_chnk, index_te, te_phys)
            call pbuf_get_field(pbuf_chnk, index_ti, ti_phys)
            ! Get components of ion drift velocities
            call pbuf_get_field(pbuf_chnk, index_ui, ui_phys)
            call pbuf_get_field(pbuf_chnk, index_vi, vi_phys)
            call pbuf_get_field(pbuf_chnk, index_wi, wi_phys)
            !--------------------------------------------------------
            ! Get ions from physics buffer if non-transported
            !--------------------------------------------------------
            if (sIndxO2p > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPO2p_phys,     &
                    start=(/1,1,sIndxO2p/), kount=(/pcols,pver,1/) )
            end if
            if (sIndxNOp > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPNOp_phys,     &
                    start=(/1,1,sIndxNOp/), kount=(/pcols,pver,1/) )
            end if
            if (sIndxN2p > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPN2p_phys,     &
                    start=(/1,1,sIndxN2p/), kount=(/pcols,pver,1/) )
            end if
            if (sIndxOp > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys,      &
                    start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
            end if

            if (sIndxFep > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPFep_phys,     &
                    start=(/1,1,sIndxFep/), kount=(/pcols,pver,1/) )
            end if
            if (sIndxMgp > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPMgp_phys,     &
                    start=(/1,1,sIndxMgp/), kount=(/pcols,pver,1/) )
            end if
            if (sIndxNap > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPNap_phys,     &
                    start=(/1,1,sIndxNap/), kount=(/pcols,pver,1/) )
            end if
            if (sIndxCap > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPCap_phys,     &
                    start=(/1,1,sIndxCap/), kount=(/pcols,pver,1/) )
            end if
            if (sIndxKp > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPKp_phys,     &
                    start=(/1,1,sIndxKp/), kount=(/pcols,pver,1/) )
            end if
            if (sIndxSip > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPSip_phys,     &
                    start=(/1,1,sIndxSip/), kount=(/pcols,pver,1/) )
            end if

            call outfld('Fep_phys0', phys_state(lchnk)%q(:,:,ixFep), pcols, lchnk)
            call outfld('Mgp_phys0', phys_state(lchnk)%q(:,:,ixMgp), pcols, lchnk)


            ! PHIS is from physics state
            phis(:ncol) = phys_state(lchnk)%phis(:ncol)
            do i = 1, ncol
               j = j + 1
               do k = 1, pver
                  ! physics state fields on levels
                  u_blck(k, j)     = phys_state(lchnk)%u(i, k)
                  v_blck(k, j)     = phys_state(lchnk)%v(i, k)
                  !------------------------------------------------------------
                  ! Might need geometric height on midpoints for output
                  !------------------------------------------------------------
                  if (hist_fld_active('Z3GM')) then
                     ! geometric altitude (meters above sea level)
                     tempm(i,k) = geometric_hgt(zgp=phys_state(lchnk)%zm(i,k), zsf=phis(i)*rga)
                  end if
                  ! physics state fields on interfaces (but only to pver)
                  zi_blck(k, j) = phys_state(lchnk)%zi(i, k) + phis(i)*rga
                  !------------------------------------------------------------
                  ! Convert geopotential to geometric height at interfaces:
                  !------------------------------------------------------------
                  ! Note: zht is pver instead of pverp because dynamo does not
                  !       use bottom interface
                  hi_blck(k,j) = geometric_hgt(zgp=phys_state(lchnk)%zi(i,k), zsf=phis(i)*rga)
                  if (hist_fld_active('Z3GMI')) then
                     tempi(i,k) = hi_blck(k, j)
                  end if
                  omega_blck(k, j) = phys_state(lchnk)%omega(i, k)
                  tn_blck(k, j)    = phys_state(lchnk)%t(i, k)
                  pmid_blck(k, j)  = phys_state(lchnk)%pmid(i, k)
                  ! Pedersen and Hall conductivities:
                  sigma_ped_blck(k, j) = sigma_ped_phys(i, k)
                  sigma_hall_blck(k, j) = sigma_hall_phys(i, k)
                  ! ion and electron temperatures
                  te_blck(k, j) = te_phys(i, k)
                  ti_blck(k, j) = ti_phys(i, k)
                  ! components of ion drift velocities
                  ui_blck(k, j) = ui_phys(i, k)
                  vi_blck(k, j) = vi_phys(i, k)
                  wi_blck(k, j) = wi_phys(i, k)
                  !------------------------------------------------------------
                  ! ions from physics state if transported, otherwise from pbuf
                  !------------------------------------------------------------
                  if (ixo2p > 0) then
                     o2pmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixo2p)
                  else if (sIndxO2p > 0) then
                     o2pmmr_blck(k, j) = mmrPO2p_phys(i, k)
                  else
                     call endrun(subname//': No source for O2p')
                  end if
                  if (ixnop > 0) then
                     nopmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixnop)
                  else if (sIndxNOp > 0) then
                     nopmmr_blck(k, j) = mmrPNOp_phys(i, k)
                  else
                     call endrun(subname//': No source for NOp')
                  end if
                  if (ixn2p > 0) then
                     n2pmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixn2p)
                  else if (sIndxN2p > 0) then
                     n2pmmr_blck(k, j) = mmrPN2p_phys(i, k)
                  else
                     call endrun(subname//': No source for N2p')
                  end if
                  if (ixop > 0) then
                     opmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixop)
                  else if (sIndxOp > 0) then
                     opmmr_blck(k, j) = mmrPOp_phys(i, k)
                  else
                     call endrun(subname//': No source for Op')
                  end if


                  if (ixFep > 0) then
                     Fepmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixFep)
                  else if (sIndxFep > 0) then
                     Fepmmr_blck(k, j) = mmrPFep_phys(i, k)
                  else
                     nullify(Fepmmr_blck)
                  end if
                  if (ixMgp > 0) then
                     Mgpmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixMgp)
                  else if (sIndxMgp > 0) then
                     Mgpmmr_blck(k, j) = mmrPMgp_phys(i, k)
                  else
                     nullify(Mgpmmr_blck)
                  end if
                  if (ixNap > 0) then
                     Napmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixNap)
                  else if (sIndxNap > 0) then
                     Napmmr_blck(k, j) = mmrPNap_phys(i, k)
                  else
                     nullify(Napmmr_blck)
                  end if
                  if (ixCap > 0) then
                     Capmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixCap)
                  else if (sIndxCap > 0) then
                     Capmmr_blck(k, j) = mmrPCap_phys(i, k)
                  else
                     nullify(Capmmr_blck)
                  end if
                  if (ixKp > 0) then
                     Kpmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixKp)
                  else if (sIndxKp > 0) then
                     Kpmmr_blck(k, j) = mmrPKp_phys(i, k)
                  else
                     nullify(Kpmmr_blck)
                  end if
                  if (ixSip > 0) then
                     Sipmmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixSip)
                  else if (sIndxSip > 0) then
                     Sipmmr_blck(k, j) = mmrPSip_phys(i, k)
                  else
                     nullify(Sipmmr_blck)
                  end if


                  opmmrtm1_blck(k, j) = opmmrtm1_phys(i, k, lchnk)

                  if (associated(O2pmmrtm1_blck)) then
                     O2pmmrtm1_blck(k, j) = O2pmmrtm1_phys(i, k, lchnk)
                  end if
                  if (associated(NOpmmrtm1_blck)) then
                     NOpmmrtm1_blck(k, j) = NOpmmrtm1_phys(i, k, lchnk)
                  end if
                  if (associated(Fepmmrtm1_blck)) then
                     Fepmmrtm1_blck(k, j) = Fepmmrtm1_phys(i, k, lchnk)
                  end if
                  if (associated(Mgpmmrtm1_blck)) then
                     Mgpmmrtm1_blck(k, j) = Mgpmmrtm1_phys(i, k, lchnk)
                  end if
                  if (associated(Napmmrtm1_blck)) then
                     Napmmrtm1_blck(k, j) = Napmmrtm1_phys(i, k, lchnk)
                  end if
                  if (associated(Capmmrtm1_blck)) then
                     Capmmrtm1_blck(k, j) = Capmmrtm1_phys(i, k, lchnk)
                  end if
                  if (associated(Kpmmrtm1_blck)) then
                     Kpmmrtm1_blck(k, j) = Kpmmrtm1_phys(i, k, lchnk)
                  end if
                  if (associated(Sipmmrtm1_blck)) then
                     Sipmmrtm1_blck(k, j) = Sipmmrtm1_phys(i, k, lchnk)
                  end if

                  !------------------------------------
                  ! neutrals from advected tracers array
                  !------------------------------------
                  o2mmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixo2)
                  o1mmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixo)
                  h1mmr_blck(k, j) = phys_state(lchnk)%q(i, k, ixh)
               end do
            end do ! do i = 1, ncol

            !------------------------------------------------------------------
            ! Save OMEGA and analytically derived geometric height
            !------------------------------------------------------------------
            if (hist_fld_active('Z3GM')) then
               tempm(ncol+1:, :) = 0.0_r8
               call outfld('Z3GM', tempm, pcols, lchnk)
            end if
            if (hist_fld_active('Z3GMI')) then
               tempi(ncol+1:, :) = 0.0_r8
               call outfld('Z3GMI', tempi, pcols, lchnk)
            end if
         end do ! do lchnk = begchunk, endchunk

         !---------------------------------------------------------------------
         ! Compute and save mean molecular weight:
         !---------------------------------------------------------------------
         mass_tmp(:,:) = o1mmr_blck(:,:) + o2mmr_blck(:,:) + h1mmr_blck(:,:)
         if (associated( Fepmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + 2.0_r8*Fepmmr_blck(:,:) ! why factor of 2 ???
         end if
         if (associated( Mgpmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + 2.0_r8*Mgpmmr_blck(:,:)
         end if
         if (associated( Napmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + 2.0_r8*Napmmr_blck(:,:)
         end if
         if (associated( Capmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + 2.0_r8*Capmmr_blck(:,:)
         end if
         if (associated( Kpmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + 2.0_r8*Kpmmr_blck(:,:)
         end if
         if (associated( Sipmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + 2.0_r8*Sipmmr_blck(:,:)
         end if

         n2mmr_blck(:,:) = max(1.0_r8 - mass_tmp(:,:) , n2min)

         mass_tmp(:,:) =                  o1mmr_blck(:,:) / rmassO1
         mass_tmp(:,:) = mass_tmp(:,:) + (o2mmr_blck(:,:) / rmassO2)
         mass_tmp(:,:) = mass_tmp(:,:) + (h1mmr_blck(:,:) / rmassH )
         mass_tmp(:,:) = mass_tmp(:,:) + (n2mmr_blck(:,:) / rmassN2)
         if (associated( Fepmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + (2._r8*Fepmmr_blck(:,:) / rmassFep)
         end if
         if (associated( Mgpmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + (2._r8*Mgpmmr_blck(:,:) / rmassMgp)
         end if
         if (associated( Napmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + (2._r8*Napmmr_blck(:,:) / rmassNap)
         end if
         if (associated( Capmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + (2._r8*Capmmr_blck(:,:) / rmassCap)
         end if
         if (associated( Kpmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + (2._r8*Kpmmr_blck(:,:) / rmassKp)
         end if
         if (associated( Sipmmr_blck )) then
            mass_tmp(:,:) = mass_tmp(:,:) + (2._r8*Sipmmr_blck(:,:) / rmassSip)
         end if

         mbar_blck(:,:) = 1.0_r8 / mass_tmp(:,:)

         if (state_debug_checks) then
            call shr_assert_in_domain(te_blck, is_nan=.false., varname="te_blck", msg="NaN found in te_blck in ionosphere_run2")
            call shr_assert_in_domain(ti_blck, is_nan=.false., varname="ti_blck", msg="NaN found in ti_blck in ionosphere_run2")
         end if

         call t_startf('d_pie_coupling')

         ! Compute geometric height and some diagnostic fields needed by
         ! the dynamo. Output some fields from physics grid
         ! This code is inside the timer as it is part of the coupling
         !
         ! waccmx ionosphere electro-dynamics -- transports O+ and
         !   provides updates to ion drift velocities (on physics grid)
         ! All fields are on physics mesh, (pver, blksize),
         !    where blksize is the total number of columns on this task

         call d_pie_coupling(omega_blck, pmid_blck, zi_blck, hi_blck,         &
              u_blck, v_blck, tn_blck, sigma_ped_blck, sigma_hall_blck,       &
              te_blck, ti_blck, mbar_blck, n2mmr_blck, o2mmr_blck,            &
              o1mmr_blck, o2pmmr_blck, nopmmr_blck, n2pmmr_blck,              &
              opmmr_blck, opmmrtm1_blck, ui_blck, vi_blck, wi_blck,           &
              rmassN2, rmassO2, rmassO1, &
              rmassO2p, rmassNOp, rmassN2p, rmassOp, 1, blksize, pver,        &
              rmassFep, rmassMgp, rmassNap, rmassCap, rmassKp, rmassSip,      &
              Fepmmr_blck, Mgpmmr_blck, Napmmr_blck, Capmmr_blck, Kpmmr_blck, Sipmmr_blck, &
              O2pmmrtm1_blck, NOpmmrtm1_blck, Fepmmrtm1_blck, Mgpmmrtm1_blck, Napmmrtm1_blck, &
              Capmmrtm1_blck, Kpmmrtm1_blck, Sipmmrtm1_blck )

         call t_stopf ('d_pie_coupling')

         if (state_debug_checks) then
            call shr_assert_in_domain(ui_blck, is_nan=.false., varname="ui_blck", msg="NaN found in ui_blck in ionosphere_run2")
            call shr_assert_in_domain(vi_blck, is_nan=.false., varname="vi_blck", msg="NaN found in vi_blck in ionosphere_run2")
            call shr_assert_in_domain(wi_blck, is_nan=.false., varname="wi_blck", msg="NaN found in wi_blck in ionosphere_run2")
            call shr_assert_in_domain(opmmr_blck, is_nan=.false., varname="opmmr_blck", msg="NaN found in opmmr_blck in ionosphere_run2")
            call shr_assert_in_domain(Fepmmr_blck, is_nan=.false., varname="Fepmmr_blck", msg="NaN found in Fepmmr_blck in ionosphere_run2")
            call shr_assert_in_domain(Mgpmmr_blck, is_nan=.false., varname="Mgpmmr_blck", msg="NaN found in Mgpmmr_blck in ionosphere_run2")
         end if

         !
         !----------------------------------------
         !  Put data back in to state or pbuf
         !----------------------------------------
         ! blocks --> physics chunks

         j = 0
         do lchnk = begchunk, endchunk
            ncol = phys_state(lchnk)%ncol
            pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

            call pbuf_get_field(pbuf_chnk, index_ui, ui_phys)
            call pbuf_get_field(pbuf_chnk, index_vi, vi_phys)
            call pbuf_get_field(pbuf_chnk, index_wi, wi_phys)
            if (sIndxOp > 0) then
               call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys,      &
                    start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/))
            end if
            do i = 1, ncol
               j = j + 1
               do k = 1, pver
                  ui_phys(i, k) = ui_blck(k, j)
                  vi_phys(i, k) = vi_blck(k, j)
                  wi_phys(i, k) = wi_blck(k, j)
                  if (ixop > 0) then
                     phys_state(lchnk)%q(i, k, ixop) = opmmr_blck(k, j)
                  else if (sIndxOp > 0) then
                     mmrPOp_phys(i, k) = opmmr_blck(k, j)
                  else
                     call endrun(subname//': No destination for Op')
                  end if
                  opmmrtm1_phys(i,k,lchnk) = opmmrtm1_blck(k,j)

                  ! Metal ions
                  if (ixFep > 0) then
                     phys_state(lchnk)%q(i, k, ixFep) = Fepmmr_blck(k, j)
                  else if (sIndxFep > 0) then
                     mmrPFep_phys(i, k) = Fepmmr_blck(k, j)
                  else
                     call endrun(subname//': No destination for Fep')
                  end if

                  if (ixMgp > 0) then
                     phys_state(lchnk)%q(i, k, ixMgp) = Mgpmmr_blck(k, j)
                  else if (sIndxMgp > 0) then
                     mmrPMgp_phys(i, k) = Mgpmmr_blck(k, j)
                  else
                     call endrun(subname//': No destination for Mgp')
                  end if
                  if (ixNap > 0) then
                     phys_state(lchnk)%q(i, k, ixNap) = Napmmr_blck(k, j)
                  else if (sIndxNap > 0) then
                     mmrPNap_phys(i, k) = Napmmr_blck(k, j)
                  end if
                  if (ixCap > 0) then
                     phys_state(lchnk)%q(i, k, ixCap) = Capmmr_blck(k, j)
                  else if (sIndxCap > 0) then
                     mmrPCap_phys(i, k) = Capmmr_blck(k, j)
                  end if
                  if (ixKp > 0) then
                     phys_state(lchnk)%q(i, k, ixKp) = Kpmmr_blck(k, j)
                  else if (sIndxKp > 0) then
                     mmrPKp_phys(i, k) = Kpmmr_blck(k, j)
                  end if
                  if (ixSip > 0) then
                     phys_state(lchnk)%q(i, k, ixSip) = Sipmmr_blck(k, j)
                  else if (sIndxSip > 0) then
                     mmrPSip_phys(i, k) = Sipmmr_blck(k, j)
                  end if

                  if (associated(O2pmmrtm1_blck)) then
                     O2pmmrtm1_phys(i,k,lchnk) = O2pmmrtm1_blck(k,j)
                  end if
                  if (associated(NOpmmrtm1_blck)) then
                     NOpmmrtm1_phys(i,k,lchnk) = NOpmmrtm1_blck(k,j)
                  end if
                  if (associated(Fepmmrtm1_blck)) then
                     Fepmmrtm1_phys(i,k,lchnk) = Fepmmrtm1_blck(k,j)
                  end if
                  if (associated(Mgpmmrtm1_blck)) then
                     Mgpmmrtm1_phys(i,k,lchnk) = Mgpmmrtm1_blck(k,j)
                  end if
                  if (associated(Napmmrtm1_blck)) then
                     Napmmrtm1_phys(i,k,lchnk) = Napmmrtm1_blck(k,j)
                  end if
                  if (associated(Capmmrtm1_blck)) then
                     Capmmrtm1_phys(i,k,lchnk) = Capmmrtm1_blck(k,j)
                  end if
                  if (associated(Kpmmrtm1_blck)) then
                     Kpmmrtm1_phys(i,k,lchnk) = Kpmmrtm1_blck(k,j)
                  end if
                  if (associated(Sipmmrtm1_blck)) then
                     Sipmmrtm1_phys(i,k,lchnk) = Sipmmrtm1_blck(k,j)
                  end if

               end do
            end do

            call outfld('Fep_phys1', phys_state(lchnk)%q(:,:,ixFep), pcols, lchnk)
            call outfld('Mgp_phys1', phys_state(lchnk)%q(:,:,ixMgp), pcols, lchnk)

            if (ionos_edyn_active) then
               call outfld('UI', ui_phys, pcols, lchnk)
               call outfld('VI', vi_phys, pcols, lchnk)
               call outfld('WI', wi_phys, pcols, lchnk)
               if (write_inithist()) then
                  call outfld('UI&IC', ui_phys, pcols, lchnk)
                  call outfld('VI&IC', vi_phys, pcols, lchnk)
                  call outfld('WI&IC', wi_phys, pcols, lchnk)
               end if
            end if
         end do

         if (associated(opmmr_blck)) then
            deallocate(opmmr_blck)
            nullify(opmmr_blck)
         end if
         if (associated(o2pmmr_blck)) then
            deallocate(o2pmmr_blck)
            nullify(o2pmmr_blck)
         end if
         if (associated(nopmmr_blck)) then
            deallocate(nopmmr_blck)
            nullify(nopmmr_blck)
         end if
         if (associated(n2pmmr_blck)) then
            deallocate(n2pmmr_blck)
            nullify(n2pmmr_blck)
         end if


         if (associated(Fepmmr_blck)) then
            deallocate(Fepmmr_blck)
            nullify(Fepmmr_blck)
         end if
         if (associated(Mgpmmr_blck)) then
            deallocate(Mgpmmr_blck)
            nullify(Mgpmmr_blck)
         end if
         if (associated(Napmmr_blck)) then
            deallocate(Napmmr_blck)
            nullify(Napmmr_blck)
         end if
         if (associated(Capmmr_blck)) then
            deallocate(Capmmr_blck)
            nullify(Capmmr_blck)
         end if
         if (associated(Kpmmr_blck)) then
            deallocate(Kpmmr_blck)
            nullify(Kpmmr_blck)
         end if
         if (associated(Sipmmr_blck)) then
            deallocate(Sipmmr_blck)
            nullify(Sipmmr_blck)
         end if


         if (associated(tempi)) then
            deallocate(tempi)
            nullify(tempi)
         end if
         if (associated(tempm)) then
            deallocate(tempm)
            nullify(tempm)
         end if
         deallocate(opmmrtm1_blck)
         nullify(opmmrtm1_blck)
         deallocate(phis)
         nullify(phis)
         deallocate(u_blck)
         nullify(u_blck)
         deallocate(v_blck)
         nullify(v_blck)
         deallocate(sigma_ped_blck)
         nullify(sigma_ped_blck)
         deallocate(sigma_hall_blck)
         nullify(sigma_hall_blck)
         deallocate(ti_blck)
         nullify(ti_blck)
         deallocate(hi_blck)
         nullify(hi_blck)
         deallocate(te_blck)
         nullify(te_blck)
         deallocate(zi_blck)
         nullify(zi_blck)
         deallocate(ui_blck)
         nullify(ui_blck)
         deallocate(vi_blck)
         nullify(vi_blck)
         deallocate(wi_blck)
         nullify(wi_blck)
         deallocate(omega_blck)
         nullify(omega_blck)
         deallocate(tn_blck)
         nullify(tn_blck)
         deallocate(n2mmr_blck)
         nullify(n2mmr_blck)
         deallocate(o2mmr_blck)
         nullify(o2mmr_blck)
         deallocate(o1mmr_blck)
         nullify(o1mmr_blck)
         deallocate(h1mmr_blck)
         nullify(h1mmr_blck)
         deallocate(mbar_blck)
         nullify(mbar_blck)
         deallocate(pmid_blck)
         nullify(pmid_blck)

      end if ionos_cpl

    contains

      subroutine init_ionstm1( state_ndx, slvd_ndx, ionstm1_phys )

        integer, intent(in) :: state_ndx
        integer, intent(in) :: slvd_ndx
        real(r8),intent(out):: ionstm1_phys(pcols, pver, begchunk:endchunk)

        real(r8), pointer :: slvd_ptr(:,:)=>null() ! pbuf fld
        type(physics_buffer_desc), pointer :: pbuf_chnk(:)

        integer :: lchnk, ncol

        do lchnk = begchunk, endchunk
           ncol = get_ncols_p(lchnk)
           if (slvd_ndx > 0) then
              pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
              call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, slvd_ptr,  &
                   start=(/1,1,slvd_ndx/), kount=(/pcols,pver,1/) )
              ionstm1_phys(:ncol,:pver,lchnk) = slvd_ptr(:ncol,:pver)
           else
              ionstm1_phys(:ncol,:pver,lchnk) = phys_state(lchnk)%q(:ncol,:pver,state_ndx)
           endif
        end do

      end subroutine init_ionstm1

    end subroutine ionosphere_run2

   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine ionosphere_init_restart(File)
      use pio,              only: file_desc_t, pio_double, pio_def_var
      use cam_pio_utils,    only: cam_pio_def_dim
      use cam_grid_support, only: cam_grid_id, cam_grid_write_attr
      use cam_grid_support, only: cam_grid_header_info_t

      type(File_desc_t), intent(inout) :: File

      integer                          :: grid_id
      integer                          :: hdimcnt, ierr, i
      integer                          :: dimids(3), ndims
      type(cam_grid_header_info_t)     :: info

      if (ionos_xport_active) then
         grid_id = cam_grid_id('physgrid')
         call cam_grid_write_attr(File, grid_id, info)
         hdimcnt = info%num_hdims()
         do i = 1, hdimcnt
            dimids(i) = info%get_hdimid(i)
         end do
         ndims = hdimcnt + 1

         call cam_pio_def_dim(File, 'lev',  pver,  dimids(ndims), existOK=.true.)

         ierr = pio_def_var(File, 'Optm1', pio_double, dimids(1:ndims), Optm1_vdesc)

         if (allocated(NOpmmrtm1_phys)) then
            ierr = pio_def_var(File, 'NOptm1', pio_double, dimids(1:ndims), NOptm1_vdesc)
         end if
         if (allocated(O2pmmrtm1_phys)) then
            ierr = pio_def_var(File, 'O2ptm1', pio_double, dimids(1:ndims), O2ptm1_vdesc)
         end if
         if (allocated(Fepmmrtm1_phys)) then
            ierr = pio_def_var(File, 'Feptm1', pio_double, dimids(1:ndims), Feptm1_vdesc)
         end if
         if (allocated(Mgpmmrtm1_phys)) then
            ierr = pio_def_var(File, 'Mgptm1', pio_double, dimids(1:ndims), Mgptm1_vdesc)
         end if
         if (allocated(Napmmrtm1_phys)) then
            ierr = pio_def_var(File, 'Naptm1', pio_double, dimids(1:ndims), Naptm1_vdesc)
         end if
         if (allocated(Capmmrtm1_phys)) then
            ierr = pio_def_var(File, 'Captm1', pio_double, dimids(1:ndims), Captm1_vdesc)
         end if
         if (allocated(Kpmmrtm1_phys)) then
            ierr = pio_def_var(File, 'Kptm1', pio_double, dimids(1:ndims), Kptm1_vdesc)
         end if
         if (allocated(Sipmmrtm1_phys)) then
            ierr = pio_def_var(File, 'Siptm1', pio_double, dimids(1:ndims), Siptm1_vdesc)
         end if

      end if
   end subroutine ionosphere_init_restart

   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine ionosphere_write_restart(File)
      use pio,              only: io_desc_t, file_desc_t, pio_write_darray
      use pio,              only: pio_double
      use cam_grid_support, only: cam_grid_id, cam_grid_write_var
      use cam_grid_support, only: cam_grid_get_decomp, cam_grid_dimensions
      use phys_grid,        only: phys_decomp

      type(file_desc_t), intent(inout) :: File

      integer                          :: ierr
      integer                          :: physgrid
      integer                          :: dims(3), gdims(3)
      integer                          :: nhdims
      type(io_desc_t), pointer         :: iodesc3d

      if (ionos_xport_active) then

         ! Write grid vars
         call cam_grid_write_var(File, phys_decomp)

         physgrid = cam_grid_id('physgrid')
         call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)
         nhdims = nhdims + 1
         gdims(nhdims) = pver
         dims(1) = pcols
         dims(2) = pver
         dims(3) = endchunk - begchunk + 1
         call cam_grid_get_decomp(physgrid, dims(1:3), gdims(1:nhdims),       &
              pio_double, iodesc3d)

         call pio_write_darray(File, Optm1_vdesc, iodesc3d, opmmrtm1_phys, ierr)

         if (allocated(NOpmmrtm1_phys)) then
            call pio_write_darray(File, NOptm1_vdesc, iodesc3d, NOpmmrtm1_phys, ierr)
         end if
         if (allocated(O2pmmrtm1_phys)) then
            call pio_write_darray(File, O2ptm1_vdesc, iodesc3d, O2pmmrtm1_phys, ierr)
         end if
         if (allocated(Fepmmrtm1_phys)) then
            call pio_write_darray(File, Feptm1_vdesc, iodesc3d, Fepmmrtm1_phys, ierr)
         end if
         if (allocated(Mgpmmrtm1_phys)) then
            call pio_write_darray(File, Mgptm1_vdesc, iodesc3d, Mgpmmrtm1_phys, ierr)
         end if
         if (allocated(Napmmrtm1_phys)) then
            call pio_write_darray(File, Naptm1_vdesc, iodesc3d, Napmmrtm1_phys, ierr)
         end if
         if (allocated(Capmmrtm1_phys)) then
            call pio_write_darray(File, Captm1_vdesc, iodesc3d, Capmmrtm1_phys, ierr)
         end if
         if (allocated(Kpmmrtm1_phys)) then
            call pio_write_darray(File, Kptm1_vdesc, iodesc3d, Kpmmrtm1_phys, ierr)
         end if
         if (allocated(Sipmmrtm1_phys)) then
            call pio_write_darray(File, Siptm1_vdesc, iodesc3d, Sipmmrtm1_phys, ierr)
         end if
      end if

   end subroutine ionosphere_write_restart

   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine ionosphere_read_restart(File)
      use pio,              only: io_desc_t, file_desc_t, pio_inq_varid
      use pio,              only: pio_read_darray, pio_double
      use cam_grid_support, only: cam_grid_id
      use cam_grid_support, only: cam_grid_get_decomp, cam_grid_dimensions

      type(file_desc_t), intent(inout) :: File

      integer                          :: ierr
      integer                          :: physgrid
      integer                          :: dims(3), gdims(3)
      integer                          :: nhdims
      type(io_desc_t), pointer         :: iodesc3d

      if (ionos_xport_active) then
         call ionosphere_alloc()

         physgrid = cam_grid_id('physgrid')
         call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)
         nhdims = nhdims + 1
         gdims(nhdims) = pver
         dims(1) = pcols
         dims(2) = pver
         dims(3) = endchunk - begchunk + 1
         call cam_grid_get_decomp(physgrid, dims(1:3), gdims(1:nhdims),       &
              pio_double, iodesc3d)

         ierr = pio_inq_varid(File, 'Optm1', Optm1_vdesc)
         call pio_read_darray(File, Optm1_vdesc, iodesc3d, opmmrtm1_phys, ierr)
         opmmrtm1_initialized = .true.

         if (allocated(NOpmmrtm1_phys)) then
            ierr = pio_inq_varid(File, 'NOptm1', NOptm1_vdesc)
            call pio_read_darray(File, NOptm1_vdesc, iodesc3d, NOpmmrtm1_phys, ierr)
            NOpmmrtm1_initialized = .true.
         end if
         if (allocated(O2pmmrtm1_phys)) then
            ierr = pio_inq_varid(File, 'O2ptm1', O2ptm1_vdesc)
            call pio_read_darray(File, O2ptm1_vdesc, iodesc3d, O2pmmrtm1_phys, ierr)
            O2pmmrtm1_initialized = .true.
         end if
         if (allocated(Fepmmrtm1_phys)) then
            ierr = pio_inq_varid(File, 'Feptm1', Feptm1_vdesc)
            call pio_read_darray(File, Feptm1_vdesc, iodesc3d, Fepmmrtm1_phys, ierr)
            Fepmmrtm1_initialized = .true.
         end if
         if (allocated(Mgpmmrtm1_phys)) then
            ierr = pio_inq_varid(File, 'Mgptm1', Mgptm1_vdesc)
            call pio_read_darray(File, Mgptm1_vdesc, iodesc3d, Mgpmmrtm1_phys, ierr)
            Mgpmmrtm1_initialized = .true.
         end if
         if (allocated(Napmmrtm1_phys)) then
            ierr = pio_inq_varid(File, 'Naptm1', Naptm1_vdesc)
            call pio_read_darray(File, Naptm1_vdesc, iodesc3d, Napmmrtm1_phys, ierr)
            Napmmrtm1_initialized = .true.
         end if
         if (allocated(Capmmrtm1_phys)) then
            ierr = pio_inq_varid(File, 'Captm1', Captm1_vdesc)
            call pio_read_darray(File, Captm1_vdesc, iodesc3d, Capmmrtm1_phys, ierr)
            Capmmrtm1_initialized = .true.
         end if
         if (allocated(Kpmmrtm1_phys)) then
            ierr = pio_inq_varid(File, 'Kptm1', Kptm1_vdesc)
            call pio_read_darray(File, Kptm1_vdesc, iodesc3d, Kpmmrtm1_phys, ierr)
            Kpmmrtm1_initialized = .true.
         end if
         if (allocated(Sipmmrtm1_phys)) then
            ierr = pio_inq_varid(File, 'Siptm1', Siptm1_vdesc)
            call pio_read_darray(File, Siptm1_vdesc, iodesc3d, Sipmmrtm1_phys, ierr)
            Sipmmrtm1_initialized = .true.
         end if

      end if

   end subroutine ionosphere_read_restart

   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine ionosphere_final

      use edyn_esmf, only: edyn_esmf_final

      call edyn_esmf_final()

      if (allocated(opmmrtm1_phys)) then
         deallocate(opmmrtm1_phys)
      end if

   end subroutine ionosphere_final

   !===========================================================================
   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine ionosphere_read_ic()

      use pio,              only: file_desc_t
      use ncdio_atm,        only: infld
      use cam_initfiles,    only: initial_file_get_id
      use cam_grid_support, only: cam_grid_check, cam_grid_id
      use cam_grid_support, only: cam_grid_get_dim_names

      type(file_desc_t), pointer :: fh_ini  ! PIO filehandle

      integer                    :: grid_id ! grid ID for data mapping
      character(len=8)           :: dim1name, dim2name
      character(len=*), parameter :: subname = 'ionosphere_read_ic'

      if ( ionos_xport_active ) then
         call ionosphere_alloc()

         fh_ini => initial_file_get_id()
         grid_id = cam_grid_id('physgrid')
         if (.not. cam_grid_check(grid_id)) then
            call endrun(trim(subname)//': Internal error, no "physgrid" grid')
         end if
         call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

         ! try reading in OpTM1 from the IC file
         call read_tm1_ion('Op', Opmmrtm1_phys, Opmmrtm1_initialized)

         if (allocated(NOpmmrtm1_phys)) then
            call read_tm1_ion('NOp', NOpmmrtm1_phys, NOpmmrtm1_initialized)
         end if
         if (allocated(O2pmmrtm1_phys)) then
            call read_tm1_ion('O2p', O2pmmrtm1_phys, O2pmmrtm1_initialized)
         end if
         if (allocated(Fepmmrtm1_phys)) then
            call read_tm1_ion('Fep', Fepmmrtm1_phys, Fepmmrtm1_initialized)
         end if
         if (allocated(Mgpmmrtm1_phys)) then
            call read_tm1_ion('Mgp', Mgpmmrtm1_phys, Mgpmmrtm1_initialized)
         end if
         if (allocated(Napmmrtm1_phys)) then
            call read_tm1_ion('Nap', Napmmrtm1_phys, Napmmrtm1_initialized)
         end if
         if (allocated(Capmmrtm1_phys)) then
            call read_tm1_ion('Cap', Capmmrtm1_phys, Capmmrtm1_initialized)
         end if
         if (allocated(Kpmmrtm1_phys)) then
            call read_tm1_ion('Kp', Kpmmrtm1_phys, Kpmmrtm1_initialized)
         end if
         if (allocated(Sipmmrtm1_phys)) then
            call read_tm1_ion('Sip', Sipmmrtm1_phys, Sipmmrtm1_initialized)
         end if

      end if

    contains

      subroutine read_tm1_ion(name, physfld, success)
        use infnan, only: nan, assignment(=)

        character(len=*), intent(in) :: name
        real(r8), intent(out) :: physfld(:,:,:)
        logical, intent(out) :: success

        success = .false.
        physfld = nan

        ! try reading in ion-TM1 from the IC file
        call infld(trim(name)//'TM1', fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, &
             begchunk, endchunk, physfld, success, gridname='physgrid')
        if (.not. success) then
           ! if ion-TM1 is not included in the IC file then try reading the ion
           call infld(trim(name), fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, &
                begchunk, endchunk, physfld, success, gridname='physgrid')
        end if

      end subroutine read_tm1_ion

   end subroutine ionosphere_read_ic

   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   subroutine ionosphere_alloc()
      use infnan, only: nan, assignment(=)
      integer :: astat

      if (.not. allocated(opmmrtm1_phys)) then
         allocate(opmmrtm1_phys(pcols, pver, begchunk:endchunk), stat=astat)
         if (astat /= 0) then
            call endrun('ionosphere_alloc: failed to allocate opmmrtm1_phys')
         end if
         opmmrtm1_phys = nan
         opmmrtm1_initialized = .false.

         ! set ion indexes and molecular masses
         call set_indices_mass('Fep', ixFep, sIndxFep, rmassFep)
         call set_indices_mass('Mgp', ixMgp, sIndxMgp, rmassMgp)
         call set_indices_mass('Nap', ixNap, sIndxNap, rmassNap)
         call set_indices_mass('Cap', ixCap, sIndxCap, rmassCap)
         call set_indices_mass('Kp',  ixKp,  sIndxKp,  rmassKp )
         call set_indices_mass('Sip', ixSip, sIndxSip, rmassSip)

         if (ixFep>0 .or. sIndxFep>0) then

            allocate(NOpmmrtm1_phys(pcols, pver, begchunk:endchunk), stat=astat)
            if (astat /= 0) then
               call endrun('ionosphere_alloc: failed to allocate NOpmmrtm1_phys')
            end if
            NOpmmrtm1_phys = nan
            NOpmmrtm1_initialized = .false.
            allocate(O2pmmrtm1_phys(pcols, pver, begchunk:endchunk), stat=astat)
            if (astat /= 0) then
               call endrun('ionosphere_alloc: failed to allocate O2pmmrtm1_phys')
            end if
            O2pmmrtm1_phys = nan
            O2pmmrtm1_initialized = .false.
            allocate(Fepmmrtm1_phys(pcols, pver, begchunk:endchunk), stat=astat)
            if (astat /= 0) then
               call endrun('ionosphere_alloc: failed to allocate Fepmmrtm1_phys')
            end if
            Fepmmrtm1_phys = nan
            Fepmmrtm1_initialized = .false.

         end if
         if (ixMgp>0 .or. sIndxMgp>0) then
            allocate(Mgpmmrtm1_phys(pcols, pver, begchunk:endchunk), stat=astat)
            if (astat /= 0) then
               call endrun('ionosphere_alloc: failed to allocate Mgpmmrtm1_phys')
            end if
            Mgpmmrtm1_phys = nan
            Mgpmmrtm1_initialized = .false.
         end if
         if (ixNap>0 .or. sIndxNap>0) then
            allocate(Napmmrtm1_phys(pcols, pver, begchunk:endchunk), stat=astat)
            if (astat /= 0) then
               call endrun('ionosphere_alloc: failed to allocate Napmmrtm1_phys')
            end if
            Napmmrtm1_phys = nan
            Napmmrtm1_initialized = .false.
         end if
         if (ixCap>0 .or. sIndxCap>0) then
            allocate(Capmmrtm1_phys(pcols, pver, begchunk:endchunk), stat=astat)
            if (astat /= 0) then
               call endrun('ionosphere_alloc: failed to allocate Capmmrtm1_phys')
            end if
            Capmmrtm1_phys = nan
            Capmmrtm1_initialized = .false.
         end if
         if (ixKp>0 .or. sIndxKp>0) then
            allocate(Kpmmrtm1_phys(pcols, pver, begchunk:endchunk), stat=astat)
            if (astat /= 0) then
               call endrun('ionosphere_alloc: failed to allocate Kpmmrtm1_phys')
            end if
            Kpmmrtm1_phys = nan
            Kpmmrtm1_initialized = .false.
         end if
         if (ixSip>0 .or. sIndxSip>0) then
            allocate(Sipmmrtm1_phys(pcols, pver, begchunk:endchunk), stat=astat)
            if (astat /= 0) then
               call endrun('ionosphere_alloc: failed to allocate Sipmmrtm1_phys')
            end if
            Sipmmrtm1_phys = nan
            Sipmmrtm1_initialized = .false.
         end if

      end if

   end subroutine ionosphere_alloc

   !==========================================================================

   ! utility routine to initialize indices and molec masses of ions
   subroutine set_indices_mass(name, cnstidx, slvdidx, molmass)

     character(len=*), intent(in) :: name ! cam constituent name
     integer, intent(out) :: cnstidx ! cam constituent index
     integer, intent(out) :: slvdidx ! short-lived species index
     real(r8),intent(out) :: molmass ! molecular masss

     integer :: idx

     cnstidx = -1
     slvdidx = -1
     molmass = 0._r8

     call cnst_get_ind(name, cnstidx, abort=.false.)
     if (cnstidx > 0) then
        molmass = cnst_mw(cnstidx)
     else
        slvdidx = slvd_index(name)
        if (slvdidx > 0) then
           idx = get_spc_ndx(name)
           molmass = adv_mass(idx)
        end if
     end if

   end subroutine set_indices_mass

   !==========================================================================

   ! calculates geometric height (meters above sea level)
   pure function geometric_hgt( zgp, zsf ) result(zgm)

     real(r8), intent(in) :: zgp ! geopotential height (m)
     real(r8), intent(in) :: zsf ! surface height above sea level (m)
     real(r8) :: zgm ! geometric height above sea level (m)

     real(r8) :: tmp

     ! Hanli's formulation:
     ! Z_gm = 1/(1 - (1+Zs/r) * Z_gp/r) * (Zs + (1+Zs/r) * Z_gp)
     ! Z_gm: geometric height
     ! Zs: Surface height
     ! Z_gp: model calculated geopotential height (zm and zi in the model)

     tmp = 1._r8+zsf*rearth_inv
     zgm = (zsf + tmp*zgp) / (1._r8 - tmp*zgp*rearth_inv)

   end function geometric_hgt

end module ionosphere_interface
