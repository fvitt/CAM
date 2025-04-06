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

  subroutine edyn3d_driver_init( mpicom_atm, npes_edyn3D )
    use mpi_module, only: mpi_init => init, setup_topology
    use mpi_module, only: mpi_rank, mpi_size, lat_size, lon_size, lat_rank, lon_rank
    use mpi_module, only: nmlon_task,mlon0_task,mlon1_task, nmlat_task,mlat0_task,mlat1_task
    use grid_module,only: generate_mag_grid
    use init_module,only: init_cons, init_fieldline, calculate_m, get_apex, calc_magcoor_geogrid
    use params_module,only: nlat,nlon,nmlat_h,nmlon, nhgt_fix
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
    call generate_mag_grid()

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
    call setup_topology(nlat,nlon,nmlat_h,nmlon)
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

    real(r8), intent(in) :: sigPed(nphyslev,nphyscol)
    real(r8), intent(in) :: sigHal(nphyslev,nphyscol)
    real(r8), intent(in) :: un(nphyslev,nphyscol)
    real(r8), intent(in) :: vn(nphyslev,nphyscol)

    real(r8) :: physflds(nphyslev,nphyscol, phys2mag_nflds)
    real(r8) :: magflds(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1, phys2mag_nflds)
    real(r8) :: magflds2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1, mag2opls_nflds)

    real(r8) :: opalt(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8) :: oplusflds(lon0:lon1,lat0:lat1,lev0:lev1,mag2opls_nflds) ! field mapped to oplus grid

    character(len=*), parameter :: subname = 'edyn3d_driver_timestep'

    print*,'FVDBG.edyn3D_driver_timestep... 0'

    physflds = NOTSET
    magflds = NOTSET
    magflds2 = NOTSET
    opalt = NOTSET
    oplusflds = NOTSET

    physflds(:,:,1) = sigPed(:,:)
    physflds(:,:,2) = sigHal(:,:)
    physflds(:,:,3) = un(:,:)
    physflds(:,:,4) = vn(:,:)

    call edyn3d_remap_phys2mag_s1(nphyscol, nphyslev, physalt, physflds, magflds)

    call edyn3d_hist_mag_s1_out('sigma_ped_s1',magflds(:,:,:,:,1))
    call edyn3d_hist_mag_s1_out('sigma_hal_s1',magflds(:,:,:,:,2))

    call edyn3d_remap_phys2mag_s2(nphyscol, nphyslev, physalt, physflds, magflds)

    call edyn3d_hist_mag_s2_out('sigma_ped_s2',magflds(:,:,:,:,1))
    call edyn3d_hist_mag_s2_out('sigma_hal_s2',magflds(:,:,:,:,2))

    magflds2(:,:,:,:,1:3) = magflds(:,:,:,:,1:3)

    call regrid_phys2geo_3d( physalt, opalt, nphyslev, 1, nphyscol )
    if (any(opalt==NOTSET)) then
       call endrun(subname//': regrid_phys2geo_3d physalt->opalt ERROR')
    end if

    call edyn3d_remap_mag2oplus( magflds2, opalt, oplusflds )
    if (any(oplusflds==NOTSET)) then
       call endrun(subname//': edyn3d_remap_mag2oplus ERROR')
    end if

    print*,'FVDBG.edyn3D_driver_timestep... END'

  end subroutine edyn3d_driver_timestep

end module edyn3d_driver_mod
