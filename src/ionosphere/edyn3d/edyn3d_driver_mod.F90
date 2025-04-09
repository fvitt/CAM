module edyn3d_driver_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc

  implicit none

  private
  public :: edyn3d_driver_init
  public :: edyn3d_driver_timestep

  real(r8), dimension(:,:,:), allocatable :: gmlat,gmlon

contains

  subroutine edyn3d_driver_init( mpicom_atm, npes_edyn3D, edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt )
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

    use infnan, only: nan, assignment(=)
    use mpi_module, only: mlon0, mlon1, mlat0, mlat1
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_grids_reg
    use cam_history, only: addfld, horiz_only

    integer, intent(in) :: mpicom_atm, npes_edyn3D
    integer, intent(in) :: edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt

    integer :: ierror
    character(len=*), parameter :: prefix = 'edyn3d_driver_init: '

    print*,'FVDBG.edyn3d_driver_init...0'
    if (r8 /= rp) then
       call endrun(prefix//'r8 /= rp')
    end if

    call mpi_init( mpicom_atm, npes_edyn3D )

    if (masterproc) then
       write(iulog,*) prefix,'3D Edyn grid params nmlat_h,nmlon,nhgt_fix: ',nmlat_h,nmlon,nhgt_fix
       write(iulog,*) prefix,'3D Edyn mpi_rank, mpi_size: ',mpi_rank,mpi_size
       write(iulog,*) prefix,'3D Edyn lon_size, lat_size: ',lon_size,lat_size
    end if

    ! set up magnetic latitude and longitude grids
    call generate_mag_grid(edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt)

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
    if (masterproc) then
       write(iulog,*) prefix,'3D Edyn nmlon_task: ',nmlon_task
       write(iulog,*) prefix,'3D Edyn mlon0_task: ',mlon0_task
       write(iulog,*) prefix,'3D Edyn mlon1_task: ',mlon1_task
       write(iulog,*) prefix,'3D Edyn nmlat_task: ',nmlat_task
       write(iulog,*) prefix,'3D Edyn mlat0_task: ',mlat0_task
       write(iulog,*) prefix,'3D Edyn mlat1_task: ',mlat1_task
    end if

    print*,'FVDBG.edyn3d_driver_init...mlon0, mlon1, mlat0, mlat1: ', mlon0, mlon1, mlat0, mlat1

    call edyn3d_hist_mag_grids_reg()

    acitve_tasks: if (mpi_rank<mpi_size) then

       print*,'FVDBG.edyn3d_driver_init...active mpi_rank: ',mpi_rank

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


    call addfld ('sigma_ped_s1', horiz_only, 'I', 'K','Ped cond. on S1 mag field line grid', &
                  gridname='magfline_s1')
    call addfld ('sigma_hal_s1', horiz_only, 'I', 'K','Hal cond. on S1 mag field line grid', &
                  gridname='magfline_s1')
    call addfld ('sigma_ped_s2', horiz_only, 'I', 'K','Ped cond. on S2 mag field line grid', &
                  gridname='magfline_s2')
    call addfld ('sigma_hal_s2', horiz_only, 'I', 'K','Hal cond. on S2 mag field line grid', &
                  gridname='magfline_s2')


    print*,'FVDBG.edyn3d_driver_init...END'

  end subroutine edyn3d_driver_init

  subroutine edyn3d_driver_timestep( nphyscol, nphyslev, physalt, sigPed, sigHal, un, vn)
    use edyn3d_remap_mod, only: edyn3d_remap_phys2mag_s1
    use edyn3d_remap_mod, only: edyn3d_remap_phys2mag_s2
    use edyn3d_remap_mod, only: edyn3d_remap_mag2oplus, NOTSET
    use edyn3d_remap_mod, only: mag_fields_bundle_t, phys_fields_bundle_t, oplus_fields_bundle_t
    use edyn3d_esmf_fields_rhandles, only: magFieldDes_s1, rh_phys2mag_s1, phys2mag_nflds
    use edyn3D_esmf_fields_rhandles, only: mag2opls_nflds
    use mpi_module, only: mlat0, mlat1, mlon0, mlon1
    use params_module,only:  nhgt_fix
    use edyn_mpi, only: lon0,lon1,lat0,lat1,lev0,lev1
    use regridder, only: regrid_phys2geo_3d, regrid_geo2phys_3d
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_s1_out
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_s2_out

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

    character(len=*), parameter :: subname = 'edyn3d_driver_timestep'

    print*,'FVDBG.edyn3D_driver_timestep... 0'

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

    call edyn3d_remap_phys2mag_s1(nphyscol, nphyslev, physalt, phys_flds_bndl, mags1_flds_bndl)

    call edyn3d_hist_mag_s1_out('sigma_ped_s1',sigped_s1)
    call edyn3d_hist_mag_s1_out('sigma_hal_s1',sighal_s1)

    call edyn3d_remap_phys2mag_s2(nphyscol, nphyslev, physalt, phys_flds_bndl, mags2_flds_bndl)

    call edyn3d_hist_mag_s2_out('sigma_ped_s2',sigped_s2)
    call edyn3d_hist_mag_s2_out('sigma_hal_s2',sighal_s2)

    magsrc_flds_bndl(1)%fld => un_s2
    magsrc_flds_bndl(2)%fld => vn_s2
    magsrc_flds_bndl(3)%fld => sighal_s2

    ui_oplus = NOTSET
    vi_oplus = NOTSET
    wi_oplus = NOTSET

    oplus_flds_bndl(1)%fld => ui_oplus
    oplus_flds_bndl(2)%fld => vi_oplus
    oplus_flds_bndl(3)%fld => wi_oplus

    call regrid_phys2geo_3d( physalt, opalt, nphyslev, 1, nphyscol )
    if (any(opalt==NOTSET)) then
       call endrun(subname//': regrid_phys2geo_3d physalt->opalt ERROR')
    end if

    call edyn3d_remap_mag2oplus( magsrc_flds_bndl, opalt, oplus_flds_bndl )
    if (any(ui_oplus==NOTSET)) then
       call endrun(subname//': edyn3d_remap_mag2oplus ERROR ui_oplus')
    end if
    if (any(vi_oplus==NOTSET)) then
       call endrun(subname//': edyn3d_remap_mag2oplus ERROR vi_oplus')
    end if
    if (any(wi_oplus==NOTSET)) then
       call endrun(subname//': edyn3d_remap_mag2oplus ERROR wi_oplus')
    end if

    print*,'FVDBG.edyn3D_driver_timestep... END'

  end subroutine edyn3d_driver_timestep

end module edyn3d_driver_mod
