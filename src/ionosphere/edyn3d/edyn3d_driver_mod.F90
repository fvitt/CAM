module edyn3d_driver_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc, mpicom
  use mpi_module, only: mpi_size, mpi_rank
  use infnan, only: nan, assignment(=)
  use perf_mod, only: t_startf, t_stopf

  implicit none

  private
  public :: edyn3d_driver_init
  public :: edyn3d_driver_timestep
  public :: edyn3d_driver_final

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_driver_init( mpicom_atm, npes_edyn3D, edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt, hilat_pot_model, wei05_coefs_file )
    use mpi_module, only: mpi_init => init, setup_topology
    use mpi_module, only: mpi_rank, mpi_size, lat_size, lon_size, lat_rank, lon_rank
    use mpi_module, only: nmlon_task,mlon0_task,mlon1_task, nmlat_task,mlat0_task,mlat1_task
    use grid_module,only: generate_mag_grid
    use init_module,only: init_cons, init_fieldline, calculate_m, get_apex
    use params_module,only: nmlat_h,nmlon, nhgt_fix
    use alloc_module,only: alloc_fieldline, alloc_fieldline_lite
    use fieldline_module,only: &
         npts_p,npts_s1,npts_s2,npts_r, &
         jmax_p,jmax_s1,jmax_s2,jmax_r, &
         size_p,size_s1,size_s2,size_r, &
         qdlat_p,qdlat_s1,qdlat_s2,qdlat_r
    use fieldline_module,only: F_p,F_s1,F_s2,F_r,M3_p,M1_s1,M2_s2,M3_r
    use edyn3d_esmf_fields_rhandles, only: edyn3d_esmf_fields_rhandles_init
    use edyn3d_esmf_phys_mesh_mod, only:  edyn3d_esmf_phys_mesh_init
    use edyn3d_esmf_oplus_grid_mod, only:  edyn3d_esmf_oplus_grid_init
    use edyn3d_esmf_s1_mag_grid_mod, only: edyn3d_esmf_s1_mag_grid_init
    use edyn3d_esmf_s2_mag_grid_mod, only: edyn3d_esmf_s2_mag_grid_init

    use fieldline_module, only: glat_p, glon_p, glat_s1, glon_s1, glat_s2, glon_s2

    use prec, only: rp
    use cons_module, only: read_fac

    use mpi_module, only: mlon0, mlon1, mlat0, mlat1
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_grids_reg
    use cam_history, only: addfld, horiz_only
    use edyn3d_highlat_potential, only: edyn3d_highlat_potential_init

    integer, intent(in) :: mpicom_atm, npes_edyn3D
    integer, intent(in) :: edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt
    character(len=*),intent(in) :: hilat_pot_model
    character(len=*),intent(in) :: wei05_coefs_file

    integer :: ierror
    character(len=*), parameter :: prefix = 'edyn3d_driver_init: '

    if (r8 /= rp) then
       call endrun(prefix//'r8 /= rp')
    end if

    ! init mpi for 3D edynamo
    call mpi_init( mpicom_atm, npes_edyn3D )

    ! set up magnetic latitude and longitude grids
    call generate_mag_grid(edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt)

    ! log grid info:
    if (masterproc) then
       write(iulog,*) prefix,'3D Edyn grid params nmlat_h,nmlon,nhgt_fix: ',nmlat_h,nmlon,nhgt_fix
       write(iulog,*) prefix,'3D Edyn mpi_rank, mpi_size: ',mpi_rank,mpi_size
       write(iulog,*) prefix,'3D Edyn lon_size, lat_size: ',lon_size,lat_size
    end if

    ! set up constants
    call init_cons()

    call alloc_fieldline_lite(ierror)
    if (ierror/=0) then
       call endrun(prefix//'alloc_fieldline failed')
    end if

    ! set up field-line grids
    call init_fieldline(npts_p,npts_s1,npts_s2,npts_r, &
         jmax_p,jmax_s1,jmax_s2,jmax_r, &
         size_p,size_s1,size_s2,size_r, &
         qdlat_p,qdlat_s1,qdlat_s2,qdlat_r)


    ! set up MPI decomposition
    call setup_topology(nmlat_h,nmlon)

    call edyn3d_hist_mag_grids_reg()

    acitve_tasks: if (mpi_rank<mpi_size) then

       ! allocate memory for fieldline data
       call alloc_fieldline(ierror)
       if (ierror/=0) then
          call endrun(prefix//'alloc_fieldline failed')
       end if

       glat_p  = -huge(1._r8)
       glon_p  = -huge(1._r8)
       glat_s1 = -huge(1._r8)
       glon_s1 = -huge(1._r8)
       glat_s2 = -huge(1._r8)
       glon_s2 = -huge(1._r8)

       ! get apex coordinates and unit vectors
       call get_apex()

       ! calculate M coefficients - P,S1,S2,R
       call calculate_m(npts_p,npts_s1,npts_s2,npts_r, &
            F_p,F_s1,F_s2,F_r,M3_p,M1_s1,M2_s2,M3_r)

    end if acitve_tasks

    call edyn3d_esmf_s2_mag_grid_init()
    call edyn3d_esmf_s1_mag_grid_init()

    call edyn3d_esmf_phys_mesh_init()
    call edyn3d_esmf_oplus_grid_init()

    call edyn3d_esmf_fields_rhandles_init()

    call addfld ('sigma_ped_s1', horiz_only, 'I', 'K','Ped cond. on S1 mag field line grid', gridname='magfline_s1')
    call addfld ('sigma_hal_s1', horiz_only, 'I', 'K','Hal cond. on S1 mag field line grid', gridname='magfline_s1')
    call addfld ('sigma_ped_s2', horiz_only, 'I', 'K','Ped cond. on S2 mag field line grid', gridname='magfline_s2')
    call addfld ('sigma_hal_s2', horiz_only, 'I', 'K','Hal cond. on S2 mag field line grid', gridname='magfline_s2')

    call addfld ('un_s1', horiz_only, 'I', 'm/s','Zonal wind on S1 mag field line grid', gridname='magfline_s1')
    call addfld ('vn_s1', horiz_only, 'I', 'm/s','Meridional wind on S1 mag field line grid', gridname='magfline_s1')
    call addfld ('un_s2', horiz_only, 'I', 'm/s','Zonal wind on S2 mag field line grid', gridname='magfline_s2')
    call addfld ('vn_s2', horiz_only, 'I', 'm/s','Meridional wind on S2 mag field line grid', gridname='magfline_s2')

    call addfld ('IonU_s1', horiz_only, 'I', 'm/s','Zonal Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('IonV_s1', horiz_only, 'I', 'm/s','Meridional Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('IonW_s1', horiz_only, 'I', 'm/s','Verical Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('IonU_s2', horiz_only, 'I', 'm/s','Zonal Ion Drift Velocity on s1 grid', gridname='magfline_s2')
    call addfld ('IonV_s2', horiz_only, 'I', 'm/s','Meridional Ion Drift Velocity on s1 grid', gridname='magfline_s2')
    call addfld ('IonW_s2', horiz_only, 'I', 'm/s','Verical Ion Drift Velocity on s1 grid', gridname='magfline_s2')

    call addfld ('ELECPOTEN', horiz_only, 'I', 'Volts','Electric potential', gridname='geomag_grid')
    call addfld ('HILAT_POT', horiz_only, 'I', 'Volts','High-Latitude potential', gridname='geomag_grid')
    call addfld ('HILAT_FAC', horiz_only, 'I', '???','High-Latitude field-aligned current', gridname='geomag_grid')

    call addfld ('MLON_TEST', horiz_only, 'I', 'deg','Test fld', gridname='geomag_grid')
    call addfld ('MLAT_TEST', horiz_only, 'I', 'deg','Test fld', gridname='geomag_grid')

    read_fac = .false. ! prescribed high-lat potential (pot_hl) will be provided

    call edyn3d_highlat_potential_init(hilat_pot_model,wei05_coefs_file)

  end subroutine edyn3d_driver_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_driver_timestep( nphyscol, nphyslev, physalt, sigPed, sigHal, un, vn)
    use edyn3d_remap_mod, only: edyn3d_remap_phys2mag_s1
    use edyn3d_remap_mod, only: edyn3d_remap_phys2mag_s2
    use edyn3d_remap_mod, only: edyn3d_remap_mag2oplus, NOTSET
    use edyn3d_remap_mod, only: mag_fields_bundle_t, phys_fields_bundle_t, oplus_fields_bundle_t
    use edyn3d_esmf_fields_rhandles, only: magFieldDes_s1, rh_phys2mag_s1, phys2mag_nflds
    use edyn3D_esmf_fields_rhandles, only: mag2opls_nflds
    use mpi_module, only: mlat0, mlat1, mlon0, mlon1
    use mpi_module, only: mlond0, mlond1, mlatd0, mlatd1
    use params_module,only:  nhgt_fix
    use edyn_mpi, only: lon0,lon1,lat0,lat1,lev0,lev1
    use regridder, only: regrid_phys2geo_3d, regrid_geo2phys_3d
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_s1_out
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_s2_out
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mlonlat_out
    use mpi_module, only: sync_mlat_5d, sync_mlon_5d
    use calculate_terms_module, only: calculate_conductance
    use calculate_terms_module, only: calculate_n, calculate_je, calculate_s
    use calculate_terms_module, only: calculate_ed, calculate_ve, calculate_vxyz

    use stencil_module, only: calculate_coef, calculate_coef_ns2, calculate_coef_ns
    use stencil_module, only: calculate_bij
    use solver_module, only: linear_system
    use edyn3d_highlat_potential, only: edyn3d_highlat_potential_get

    !use fieldline_module,only: npts_s1,npts_s2, bmag_s1, bmag_s2, vmp_s1, vmp_s2
    !use fieldline_module,only: D_s1,M1_s1,d1d1_s1,d1d2_s1,d2d2_s1, D_s2,M2_s2,d1d2_s2,d2d2_s2
    !use fieldline_module,only: be3_s1, be3_s2, d1_s1, d2_s1, d1_s2, d2_s2
    use fieldline_module
    use params_module, only: ylonm, ylatm
    use cons_module, only: rtd

    integer,  intent(in) :: nphyscol, nphyslev
    real(r8), intent(in) :: physalt(nphyslev,nphyscol)

    real(r8), target, intent(in) :: sigPed(nphyslev,nphyscol)
    real(r8), target, intent(in) :: sigHal(nphyslev,nphyscol)
    real(r8), target, intent(in) :: un(nphyslev,nphyscol)
    real(r8), target, intent(in) :: vn(nphyslev,nphyscol)

    real(r8) :: opalt(lon0:lon1,lat0:lat1,lev0:lev1)

    type(phys_fields_bundle_t) :: phys_flds_bndl(phys2mag_nflds)
    type(mag_fields_bundle_t) :: mags1_flds_bndl(phys2mag_nflds)
    type(mag_fields_bundle_t) :: mags2_flds_bndl(phys2mag_nflds)
    type(mag_fields_bundle_t) :: magsrc_flds_bndl(mag2opls_nflds)
    type(oplus_fields_bundle_t) :: oplus_flds_bndl(mag2opls_nflds)

    real(r8), target :: sigped_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: sighal_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: un_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: vn_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8), target :: sigped_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: sighal_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: un_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: vn_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8), target :: ui_oplus(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8), target :: vi_oplus(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8), target :: wi_oplus(lon0:lon1,lat0:lat1,lev0:lev1)

    real(r8) :: tmp_ghost(4,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    integer :: i, rc , isn, j

    real(r8) :: sigP_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: zigP_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: sigH_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: zigH_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: sigP_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: zigP_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: sigH_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: zigH_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(r8) :: ntlU_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ntlV_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ntlU_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ntlV_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(r8) :: N1p_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: N1h_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: Je1D_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: N2p_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: N2h_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: Je2D_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: S_p(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(r8) :: coef(10,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: coef_ns2(10,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: coef_ns(10,2,mlatd0:mlatd1,mlond0:mlond1)

    real(r8) :: bij(mlatd0:mlatd1,mlond0:mlond1)

    real(r8) :: pot_hl_p(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: fac_hl_p(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: pot_p(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ed1_s1(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ed2_s1(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ve1_s1(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ve2_s1(2,mlatd0:mlatd1,mlond0:mlond1)

    real(r8) :: ed1_s2(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ed2_s2(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ve1_s2(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: ve2_s2(2,mlatd0:mlatd1,mlond0:mlond1)

    real(r8) :: vx_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: vy_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: vz_s1(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: vx_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: vy_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8) :: vz_s2(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1)

    real(r8) :: maglon(2,mlat0:mlat1,mlon0:mlon1)
    real(r8) :: maglat(2,mlat0:mlat1,mlon0:mlon1)

    logical,parameter :: setbij = .true.

    character(len=*), parameter :: subname = 'edyn3d_driver_timestep'

    call t_startf(subname)

    sigped_s1 = NOTSET
    sighal_s1 = NOTSET
    un_s1 = NOTSET
    vn_s1 = NOTSET

    sigped_s2 = NOTSET
    sighal_s2 = NOTSET
    un_s2 = NOTSET
    vn_s2 = NOTSET

    mags1_flds_bndl(1)%fld => sigped_s1
    mags1_flds_bndl(2)%fld => sighal_s1
    mags1_flds_bndl(3)%fld => un_s1
    mags1_flds_bndl(4)%fld => vn_s1

    mags2_flds_bndl(1)%fld => sigped_s2
    mags2_flds_bndl(2)%fld => sighal_s2
    mags2_flds_bndl(3)%fld => un_s2
    mags2_flds_bndl(4)%fld => vn_s2

    phys_flds_bndl(1)%fld => sigPed
    phys_flds_bndl(2)%fld => sigHal
    phys_flds_bndl(3)%fld => un
    phys_flds_bndl(4)%fld => vn

    call t_startf(subname//'->remap_phys2mag_s1')
    call edyn3d_remap_phys2mag_s1(nphyscol, nphyslev, physalt, phys_flds_bndl, mags1_flds_bndl)
    call t_stopf(subname//'->remap_phys2mag_s1')

    call edyn3d_hist_mag_s1_out('sigma_ped_s1',sigped_s1)
    call edyn3d_hist_mag_s1_out('sigma_hal_s1',sighal_s1)
    call edyn3d_hist_mag_s1_out('un_s1',un_s1)
    call edyn3d_hist_mag_s1_out('vn_s1',vn_s1)

    call t_startf(subname//'->remap_phys2mag_s2')
    call edyn3d_remap_phys2mag_s2(nphyscol, nphyslev, physalt, phys_flds_bndl, mags2_flds_bndl)
    call t_stopf(subname//'->remap_phys2mag_s2')

    call edyn3d_hist_mag_s2_out('sigma_ped_s2',sigped_s2)
    call edyn3d_hist_mag_s2_out('sigma_hal_s2',sighal_s2)
    call edyn3d_hist_mag_s2_out('un_s2',un_s2)
    call edyn3d_hist_mag_s2_out('vn_s2',vn_s2)

    if (mpi_rank<mpi_size) then


       sigP_s1 = nan
       zigP_s1 = nan
       sigH_s1 = nan
       zigH_s1 = nan

       sigP_s2 = nan
       zigP_s2 = nan
       sigH_s2 = nan
       zigH_s2 = nan

       ntlU_s1 = nan
       ntlV_s1 = nan
       ntlU_s2 = nan
       ntlV_s2 = nan

       N1p_s1 = nan
       N1h_s1 = nan
       Je1D_s1 = nan
       N2p_s2 = nan
       N2h_s2 = nan
       Je2D_s2 = nan

       S_p = nan
       coef = nan
       coef_ns = nan
       coef_ns2 = nan
       bij = nan

       call t_startf(subname//'->ghost_exchange')
       ! exchange S1 ghost points
       tmp_ghost = nan
       tmp_ghost(1,:,:,mlat0:mlat1,mlon0:mlon1) = sigped_s1(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(2,:,:,mlat0:mlat1,mlon0:mlon1) = sighal_s1(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(3,:,:,mlat0:mlat1,mlon0:mlon1) = un_s1(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(4,:,:,mlat0:mlat1,mlon0:mlon1) = vn_s1(:,:,mlat0:mlat1,mlon0:mlon1)
       call sync_mlat_5d(tmp_ghost(:,:,:,:,mlon0:mlon1), 4, nhgt_fix, 2)
       call sync_mlon_5d(tmp_ghost, 4, nhgt_fix, 2)

       sigP_s1(:,:,:,:) = tmp_ghost(1,:,:,:,:)
       sigH_s1(:,:,:,:) = tmp_ghost(2,:,:,:,:)
       ntlU_s1(:,:,:,:) = tmp_ghost(3,:,:,:,:)
       ntlV_s1(:,:,:,:) = tmp_ghost(4,:,:,:,:)

       ! exchange S2 ghost points
       tmp_ghost = nan
       tmp_ghost(1,:,:,mlat0:mlat1,mlon0:mlon1) = sigped_s2(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(2,:,:,mlat0:mlat1,mlon0:mlon1) = sighal_s2(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(3,:,:,mlat0:mlat1,mlon0:mlon1) = un_s2(:,:,mlat0:mlat1,mlon0:mlon1)
       tmp_ghost(4,:,:,mlat0:mlat1,mlon0:mlon1) = vn_s2(:,:,mlat0:mlat1,mlon0:mlon1)
       call sync_mlat_5d(tmp_ghost(:,:,:,:,mlon0:mlon1), 4, nhgt_fix, 2)
       call sync_mlon_5d(tmp_ghost, 4, nhgt_fix, 2)
       sigP_s2(:,:,:,:) = tmp_ghost(1,:,:,:,:)
       sigH_s2(:,:,:,:) = tmp_ghost(2,:,:,:,:)
       ntlU_s2(:,:,:,:) = tmp_ghost(3,:,:,:,:)
       ntlV_s2(:,:,:,:) = tmp_ghost(4,:,:,:,:)
       call t_stopf(subname//'->ghost_exchange')

       ! calculate field-line integrated conductance - S1,S2
       call calculate_conductance( &
            mlatd0,mlatd1,mlond0,mlond1, &
            npts_s1,npts_s2, &
            vmp_s1,bmag_s1,sigP_s1,sigH_s1, &
            vmp_s2,bmag_s2,sigP_s2,sigH_s2, &
            zigP_s1,zigH_s1,zigP_s2,zigH_s2 )

       ! calculate N coefficients - S1,S2
       call calculate_n( &
            mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
            D_s1,M1_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1, &
            D_s2,M2_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2, &
            N1p_s1,N1h_s1,N2p_s2,N2h_s2)

       ! calculate JeD-coefficients (right hand side) - S1,S2
       call calculate_je( &
            mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
            D_s1,be3_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1,ntlU_s1,ntlV_s1, &
            D_s2,be3_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2,ntlU_s2,ntlV_s2, &
            d1_s1,d2_s1,d1_s2,d2_s2,Je1D_s1,Je2D_s2)

       ! calculate S (right hand side)
       S_p = calculate_s(mlatd0,mlatd1,mlond0,mlond1, &
            npts_p,M1_s1,Je1D_s1,M2_s2,Je2D_s2,M3_r)

       ! calculate height-dependent matrix coefficients
       coef = calculate_coef(mlatd0,mlatd1,mlond0,mlond1, &
            npts_p,S_p,N1p_s1,N1h_s1,N2p_s2,N2h_s2)

       ! add the coefficients in height to get coefficients for each hemisphere
       coef_ns2 = calculate_coef_ns2(mlatd0,mlatd1,mlond0,mlond1,coef)

       ! set the coefficient matrix in both hemispheres
       coef_ns = calculate_coef_ns(mlatd0,mlatd1,mlond0,mlond1,coef_ns2)

       ! set field-aligned conductance (b) matrix
       if (setbij) then
          bij = calculate_bij(mlatd0,mlatd1,mlond0,mlond1,coef_ns2)
       else
          bij = 0._r8
       endif

       fac_hl_p = 0._r8

       pot_p = nan

       call edyn3d_highlat_potential_get(pot_hl_p)

       call edyn3d_hist_mlonlat_out('HILAT_POT', pot_hl_p(1:2,mlat0:mlat1,mlon0:mlon1))

       ! construct linear system and solve
       call t_startf(subname//'->linear_system_solve')
       call linear_system(mlatd0,mlatd1,mlond0,mlond1, bij,pot_hl_p,fac_hl_p,coef_ns,pot_p)
       call t_stopf(subname//'->linear_system_solve')

       call edyn3d_hist_mlonlat_out('HILAT_FAC',fac_hl_p(1:2,mlat0:mlat1,mlon0:mlon1))
       call edyn3d_hist_mlonlat_out('ELECPOTEN', pot_p(1:2,mlat0:mlat1,mlon0:mlon1))

       ! calculate electric fields
       call calculate_ed( &
            mlatd0,mlatd1,mlond0,mlond1, &
            pot_p,ed1_s1,ed2_s1,ed1_s2,ed2_s2)

       ! calculate drift velocity
       call calculate_ve( &
            mlatd0,mlatd1,mlond0,mlond1, &
            ed1_s1,ed2_s1,be3_s1(1,:,:,:), &
            ed1_s2,ed2_s2,be3_s2(1,:,:,:), &
            ve1_s1,ve2_s1,ve1_s2,ve2_s2)

       ! calculate drift velocity in geographic coordinates
       call calculate_vxyz( &
            mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
            ve1_s1,ve2_s1,e1_s1,e2_s1, &
            ve1_s2,ve2_s2,e1_s2,e2_s2, &
            vx_s1,vy_s1,vz_s1,vx_s2,vy_s2,vz_s2)

       call edyn3d_hist_mag_s1_out('IonU_s1',vx_s1(:,:,mlat0:mlat1,mlon0:mlon1))
       call edyn3d_hist_mag_s1_out('IonV_s1',vy_s1(:,:,mlat0:mlat1,mlon0:mlon1))
       call edyn3d_hist_mag_s1_out('IonW_s1',vz_s1(:,:,mlat0:mlat1,mlon0:mlon1))

       call edyn3d_hist_mag_s2_out('IonU_s2',vx_s2(:,:,mlat0:mlat1,mlon0:mlon1))
       call edyn3d_hist_mag_s2_out('IonV_s2',vy_s2(:,:,mlat0:mlat1,mlon0:mlon1))
       call edyn3d_hist_mag_s2_out('IonW_s2',vz_s2(:,:,mlat0:mlat1,mlon0:mlon1))

       do isn = 1,2
          do j = mlat0,mlat1
             do i = mlon0,mlon1

                maglat(isn,j,i) = ylatm(isn,j) * rtd
                maglon(isn,j,i) = ylonm(i) * rtd

             end do
          end do
       end do

       call edyn3d_hist_mlonlat_out('MLON_TEST', maglon )
       call edyn3d_hist_mlonlat_out('MLAT_TEST', maglat )

    end if

    call t_stopf(subname)

  end subroutine edyn3d_driver_timestep

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_driver_final()
    use edyn3d_esmf_fields_rhandles, only: edyn3d_esmf_fields_rhandles_destroy
    use edyn3d_esmf_oplus_grid_mod, only: edyn3d_esmf_oplus_grid_destroy
    use edyn3d_esmf_phys_mesh_mod, only: edyn3d_esmf_phys_mesh_destroy
    use edyn3d_esmf_s1_mag_grid_mod, only: edyn3d_esmf_s1_mag_grid_destroy
    use edyn3d_esmf_s2_mag_grid_mod, only: edyn3d_esmf_s2_mag_grid_destroy

    call edyn3d_esmf_fields_rhandles_destroy()
    call edyn3d_esmf_oplus_grid_destroy()
    call edyn3d_esmf_phys_mesh_destroy()
    call edyn3d_esmf_s1_mag_grid_destroy()
    call edyn3d_esmf_s2_mag_grid_destroy()

  end subroutine edyn3d_driver_final

end module edyn3d_driver_mod
