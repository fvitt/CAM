!!XXgoldyXX: Memory leak every time edyn_esmf_update is called?
module edyn_esmf
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use ppgrid,         only: pdim1e => pcols
   use ppgrid,         only: pdim2s => begchunk, pdim2e => endchunk
   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use infnan,         only: nan, assignment(=)

#ifdef WACCMX_EDYN_ESMF
   use ESMF,           only: ESMF_Grid, ESMF_Mesh, ESMF_Field, ESMF_RouteHandle
   use ESMF,           only: ESMF_SUCCESS
   use ESMF,           only: ESMF_KIND_R8, ESMF_KIND_I4
   use ESMF,           only: ESMF_FieldGet, ESMF_GridWriteVTK
   use ESMF,           only: ESMF_STAGGERLOC_CENTER, ESMF_FieldRegridStore, ESMF_FieldRegrid
   use ESMF,           only: ESMF_STAGGERLOC_CORNER, ESMF_StaggerLoc
   use ESMF,           only: ESMF_REGRIDMETHOD_BILINEAR, ESMF_POLEMETHOD_ALLAVG
   use ESMF,           only: ESMF_REGRIDMETHOD_CONSERVE, ESMF_POLEMETHOD_NONE
   use ESMF,           only: ESMF_GridCreate1PeriDim, ESMF_INDEX_GLOBAL
   use ESMF,           only: ESMF_GridAddCoord, ESMF_GridGetCoord
   use ESMF,           only: ESMF_TYPEKIND_R8, ESMF_FieldCreate, ESMF_Array
   use ESMF,           only: ESMF_ArraySpec, ESMF_DistGrid, ESMF_DELayout
   use ESMF,           only: ESMF_GridGet, ESMF_ArraySpecSet
   use ESMF,           only: ESMF_ArrayCreate
   use ESMF,           only: ESMF_GridComp, ESMF_TERMORDER_SRCSEQ
   use ESMF,           only: ESMF_EXTRAPMETHOD_NEAREST_IDAVG
   use ESMF,           only: ESMF_UNMAPPEDACTION_IGNORE
   use edyn_mpi,       only: ntask, ntaski, ntaskj, tasks, lon0, lon1, lat0
   use edyn_mpi,       only: lat1, nmagtaski, nmagtaskj, mlon0, mlon1
   use edyn_mpi,       only: mlat0,mlat1
   use getapex,        only: gdlatdeg, gdlondeg
   use getapex,        only: gdlatdeg_corners, gdlondeg_corners
   ! dynamically allocated geo grid for Oplus transport model
   use edyn_geogrid,   only: nlon, nlev, glon, glat
   use edyn_maggrid,   only: gmlat, gmlon
   use spmd_utils,     only: masterproc
#endif

   implicit none
   save
   private

   public :: pdim1s, pdim1e, pdim2s, pdim2e
   public :: edyn_esmf_update

#ifdef WACCMX_EDYN_ESMF

   public :: edyn_esmf_final       ! Clean up any edyn usage of ESMF

   public :: edyn_esmf_regrid_phys2geo
   public :: edyn_esmf_regrid_geo2phys
   public :: edyn_esmf_regrid_phys2mag
   public :: edyn_esmf_regrid_mag2phys
   public :: edyn_esmf_regrid_geo2mag
   public :: edyn_esmf_regrid_mag2geo

   public :: edyn_esmf_get_1dfield
   public :: edyn_esmf_get_2dfield ! Retrieve a pointer to 2D ESMF field data
   public :: edyn_esmf_get_3dfield ! Retrieve a pointer to 3D ESMF field data
   public :: edyn_esmf_set3d_geo   ! Set ESMF field with 3D geo data
   public :: edyn_esmf_set2d_geo   ! Set ESMF field with 2D geo data
   public :: edyn_esmf_set3d_mag   ! Set ESMF field with 3D mag field data
   public :: edyn_esmf_set2d_mag   ! Set ESMF field with 2D mag field data
   public :: edyn_esmf_set3d_phys  ! Set ESMF field with 3D physics field data
   public :: edyn_esmf_set2d_phys  ! Set ESMF field with 2D physics field data
   public :: edyn_esmf_update_step ! .true. iff this is an update step
   public :: edyn_esmf_update_flag ! Set value of edyn_esmf_update_step
   public :: edyn_esmf_update_phys_mesh

   public :: phys_3dfld, phys_2dfld
   public :: geo_3dfld, geo_2dfld, geo2phys_3dfld
   public :: mag_des_3dfld, mag_des_2dfld
   public :: mag_src_3dfld, mag2phys_2dfld

   public :: edyn_esmf_chkerr

   !! Distribution information for grids
   type(ESMF_DELayout)  :: phys_delayout
   type(ESMF_DistGrid)  :: geo_distgrid
   type(ESMF_DistGrid)  :: mag_distgrid

   type(ESMF_Grid) :: &
        mag_src_grid, & ! source grid (will not have periodic pts)
        mag_des_grid, & ! destination grid (will have periodic pts)
        geo_grid        ! geographic grid for Oplus transport
   type(ESMF_Grid) :: geo2phys_grid
   type(ESMF_Grid) :: mag2phys_grid

   ! phys_mesh: Mesh representation of physics decomposition
   type(ESMF_Mesh), public, protected :: phys_mesh

   ! ESMF fields used for mapping between physics, oplus geographic, and geomagnetic grids
   type(ESMF_Field) :: phys_3dfld, phys_2dfld
   type(ESMF_Field) :: geo_3dfld, geo_2dfld, geo2phys_3dfld
   type(ESMF_Field) :: mag_des_3dfld, mag_des_2dfld, mag2phys_2dfld
   type(ESMF_Field) :: mag_src_3dfld


   type(ESMF_RouteHandle) ::     & ! ESMF route handles for regridding
        routehandle_phys2geo,    & ! for physics to geo 3-D regrid
        routehandle_geo2phys,    & ! for geo to physics 3-D regrid
        routehandle_phys2mag,    & ! for physics to mag 3-D regrid
!!$        routehandle_mag2phys,    & ! for mag to physics 3-D regrid
        routehandle_geo2mag,     & ! for geo to mag 3-D regrid
        routehandle_mag2geo,     & ! for geo to mag 3-D regrid
        routehandle_phys2mag_2d, & ! for 2d geo to phys
        routehandle_mag2phys_2d, & ! for 2d phys to geo for AMIE fields
!!$        routehandle_geo2phys_2d, & ! for 2d mag to phys
        routehandle_phys2geo_2d, & ! for 2d phys to geo
        routehandle_geo2mag_2d     ! for 2d geo to mag

   !
   real(r8), allocatable :: unitv(:)
   !

   logical, protected :: edyn_esmf_update_step = .true.
   logical, parameter :: debug = .false.
#endif
   integer, parameter :: pdim1s = 1

contains

#ifdef WACCMX_EDYN_ESMF
   subroutine edyn_esmf_chkerr(subname, routine, rc)
      use ESMF, only: ESMF_SUCCESS

      character(len=*), intent(in) :: subname
      character(len=*), intent(in) :: routine
      integer,          intent(in) :: rc

      character(len=512)           :: errmsg

      if (rc /= ESMF_SUCCESS) then
         write(errmsg, '(4a,i0)') subname, ': Error return from ',            &
              trim(routine), ', rc = ', rc
         if (masterproc) then
            write(iulog, '(2a)') 'ERROR: ', trim(errmsg)
         end if
         call endrun(trim(errmsg))
      end if
   end subroutine edyn_esmf_chkerr

   subroutine edyn_create_physmesh(mesh_out)
      use ESMF,         only: ESMF_DistGridCreate, ESMF_DistGrid
      use ESMF,         only: ESMF_FILEFORMAT_ESMFMESH, ESMF_MeshCreate, ESMF_MeshGet
      use phys_control, only: phys_getopts
      use phys_grid,    only: get_ncols_p, get_gcol_p, get_rlon_all_p, get_rlat_all_p
      use ppgrid,       only: pcols, begchunk, endchunk
      use shr_const_mod,only: shr_const_pi

      !
      ! Args:
      type(ESMF_Mesh), intent(out) :: mesh_out

      integer                       :: rc
      integer                       :: ncols
      integer                       :: lsize
      integer                       :: chnk, col, dindex
      integer,          allocatable :: decomp(:)
      character(len=256)            :: grid_file
      type(ESMF_DistGrid)           :: dist_grid
      character(len=*), parameter   :: subname = 'edyn_create_physmesh'


      integer                 :: spatialDim
      integer                 :: numOwnedElements
      real(r8), pointer       :: ownedElemCoords(:)
      real(r8), pointer       :: lat(:), latMesh(:)
      real(r8), pointer       :: lon(:), lonMesh(:)
      real(r8)                :: lats(pcols)                       ! array of chunk latitudes
      real(r8)                :: lons(pcols)                       ! array of chunk longitude
      integer :: i, c, n
      real(r8)        , parameter :: radtodeg = 180.0_r8/shr_const_pi
      character(len=80)       :: tempc1,tempc2
      character(len=300)      :: errstr

      ! Find the physics grid file
      call phys_getopts(physics_grid_out=grid_file)
      ! Compute the local decomp
      lsize = 0
      do chnk = begchunk, endchunk
         lsize = lsize + get_ncols_p(chnk)
      end do
      allocate(decomp(lsize))
      dindex = 0
      do chnk = begchunk, endchunk
         ncols = get_ncols_p(chnk)
         do col = 1, ncols
            dindex = dindex + 1
            decomp(dindex) = get_gcol_p(chnk, col)
         end do
      end do

      ! Create a DistGrid based on the physics decomp
      dist_grid = ESMF_DistGridCreate(arbSeqIndexList=decomp, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_DistGridCreate', rc)
      ! Create an ESMF_mesh for the physics decomposition
      mesh_out = ESMF_MeshCreate(trim(grid_file), ESMF_FILEFORMAT_ESMFMESH, &
           elementDistgrid=dist_grid, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_MeshCreateFromFile', rc)

      ! Check that the mesh coordinates are consistent with the model physics column coordinates

      ! obtain mesh lats and lons
      call ESMF_MeshGet(mesh_out, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_MeshGet', rc)

      if (numOwnedElements /= lsize) then
         write(tempc1,'(i10)') numOwnedElements
         write(tempc2,'(i10)') lsize
         call endrun(trim(subname)//": ERROR numOwnedElements "// &
                     trim(tempc1) //" not equal to local size "// trim(tempc2))
      end if

      allocate(ownedElemCoords(spatialDim*numOwnedElements))
      allocate(lonMesh(lsize), latMesh(lsize))
      call ESMF_MeshGet(mesh_out, ownedElemCoords=ownedElemCoords)

      do n = 1,lsize
         lonMesh(n) = ownedElemCoords(2*n-1)
         latMesh(n) = ownedElemCoords(2*n)
      end do

      ! obtain internally generated cam lats and lons
      allocate(lon(lsize)); lon(:) = 0._r8
      allocate(lat(lsize)); lat(:) = 0._r8
      n=0
      do c = begchunk, endchunk
         ncols = get_ncols_p(c)
         ! latitudes and longitudes returned in radians
         call get_rlat_all_p(c, ncols, lats)
         call get_rlon_all_p(c, ncols, lons)
         do i=1,ncols
            n = n+1
            lat(n) = lats(i)*radtodeg
            lon(n) = lons(i)*radtodeg
         end do
      end do

      errstr = ''
      ! error check differences between internally generated lons and those read in
      do n = 1,lsize
         if (abs(lonMesh(n) - lon(n)) > 0.000001_r8) then
            if ( (abs(lonMesh(n)-lon(n)) > 360.000001_r8) .or. (abs(lonMesh(n)-lon(n)) < 359.99999_r8) ) then
               write(errstr,100) n,lon(n),lonMesh(n), abs(lonMesh(n)-lon(n))
               write(*,*) trim(errstr)
            endif
         end if
         if (abs(latMesh(n) - lat(n)) > 0.000001_r8) then ! 1.e-12_r8
            if (.not.( (abs(lat(n))>88.0_r8) .and. (abs(latMesh(n))>88.0_r8) )) then
               write(errstr,101) n,lat(n),latMesh(n), abs(latMesh(n)-lat(n))
               write(*,*) trim(errstr)
            endif
         end if
      end do

      if ( len_trim(errstr) > 0 ) then
         call endrun(subname//': physics mesh coords do not match model coords')
      end if

      ! deallocate memory
      deallocate(ownedElemCoords)
      deallocate(lon, lonMesh)
      deallocate(lat, latMesh)
      deallocate(decomp)

100   format('edyn_create_physmesh: coord mismatch... n, lon(n), lonmesh(n), diff_lon = ',i6,2(f21.13,3x),d21.5)
101   format('edyn_create_physmesh: coord mismatch... n, lat(n), latmesh(n), diff_lat = ',i6,2(f21.13,3x),d21.5)

   end subroutine edyn_create_physmesh

   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_final

   end subroutine edyn_esmf_final

#endif

   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_update
      use getapex, only: get_apex, magfield, alonm
      use mo_apex, only: geomag_year_updated

#ifdef WACCMX_EDYN_ESMF
      ! Create ESMF grids for physics, geographic (ion transport), and
      ! magnetic grids, and create ESMF fields as necessary on each grid.
      ! Define the 2d coordinates for each grid, and save an ESMF
      ! routehandles for phys2mag, mag2phys, phys2geo, and geo2phys regridding.
      ! Note: In the future, mag2geo might be added
      !
      ! Local:
      integer :: rc ! return code for ESMF calls
      integer :: lbnd_destgeo(3), ubnd_destgeo(3) ! 3d bounds of dest geo grid
      integer :: lbnd_destmag(3), ubnd_destmag(3) ! 3d bounds of dest mag grid
      integer :: lbnd_srcgeo(3), ubnd_srcgeo(3)   ! 3d bounds of src geo grid
      integer :: lbnd_srcmag(3), ubnd_srcmag(3)   ! 3d bounds of src mag grid
      real(ESMF_KIND_R8),    pointer :: fptr(:,:,:)
      integer(ESMF_KIND_I4), pointer :: factorIndexList(:,:)
      real(ESMF_KIND_R8),    pointer :: factorList(:)
      integer                        :: smm_srctermproc,  smm_pipelinedep
#endif
      character(len=*), parameter :: subname = 'edyn_esmf_update'

      if (.not.geomag_year_updated .and. allocated(alonm)) then
         return
      end if
      !
      ! Get apex coordinates.
      !
      call get_apex()     ! get apex coordinates (allocates alonm)
      call magfield()     ! calculate magnetic field parameters

#ifdef WACCMX_EDYN_ESMF

      smm_srctermproc = 0
      smm_pipelinedep = 16
      !
      ! Set unit vector (this routine called once per run unless crossing
      !   year boundary):
      ! Handle year boundary by checking if field is allocated
      !
      if (.not.allocated(unitv)) then
         allocate(unitv(nlon))
      end if
      unitv(:) = 1._r8
      !
      ! Make geographic grid for phys2geo and geo2phys regridding:
      !
      call create_geo_grid(geo_grid)  ! geo (Oplus) grid
      call create_geo2phys_grid(geo2phys_grid)  ! geo (Oplus) grid
      !
      ! Make magnetic grid for phys2mag regridding:
      !
      call create_mag_grid(mag_des_grid, 'des')  ! mag destination grid
      !
      ! Make grid for mag2phys regridding:
      !
      call create_mag_grid(mag_src_grid, 'src')
      call create_mag2phys_grid(mag2phys_grid)
      !
      ! Make grid for mag2phys, phys2mag, geo2phys, and phys2geo regridding:
      !
      call edyn_create_physmesh(phys_mesh)
      !
      ! Create empty fields on geographic grid or phyiscs mesh that
      !   will be transformed to the magnetic grid and passed as input
      !   to the dynamo. This does not assign any values.
      !
      ! 3d fields (inputs to edynamo) on physics mesh for phys2mag:
      !
      call edyn_esmf_create_physfield(phys_2dfld, phys_mesh, 'PHYS_2DFLD', 0)
      call edyn_esmf_create_physfield(phys_3dfld, phys_mesh, 'PHYS_3DFLD', nlev)

      call edyn_esmf_create_geofield(geo_2dfld, geo_grid, 'GEO_2DFLD', 0)
      call edyn_esmf_create_geofield(geo_3dfld, geo_grid, 'GEO_3DFLD', nlev)
      call edyn_esmf_create_geofield(geo2phys_3dfld, geo2phys_grid, 'GEO2PHYS_3DFLD', nlev)

      call edyn_esmf_create_magfield(mag_des_2dfld, mag_des_grid, 'MAG_DES_2DFLD', 0)
      call edyn_esmf_create_magfield(mag_des_3dfld, mag_des_grid, 'MAG_DES_3DFLD', nlev)

      call edyn_esmf_create_magfield(mag2phys_2dfld, mag2phys_grid, 'MAG2PHYS_2DFLD', 0)
      call edyn_esmf_create_magfield(mag_src_3dfld, mag_src_grid, 'MAG_SRC_3DFLD', nlev)

      if (debug .and. masterproc) then
         !
         ! Get 3d bounds of source geo field:
         !
         call ESMF_FieldGet(geo_3dfld, localDe=0, farrayPtr=fptr,               &
              computationalLBound=lbnd_srcgeo,                                   &
              computationalUBound=ubnd_srcgeo, rc=rc)

         write(iulog,"(2a,i4,2(', ',i4),a,i4,2(', ',i4),a)") subname,         &
              ': Bounds of source geo field: lbnd_srcgeo = (',                &
              lbnd_srcgeo, '), ubnd_srcgeo = (', ubnd_srcgeo,')'
      end if

      if (debug .and. masterproc) then
         !
         ! Get 3d bounds of destination mag field:
         !
         call ESMF_FieldGet(mag_des_3dfld, localDe=0, farrayPtr=fptr,                  &
              computationalLBound=lbnd_destmag,                                  &
              computationalUBound=ubnd_destmag, rc=rc)

         write(iulog,"(2a,3i4,a,3i4,' gmlon=',2f9.3)") subname,               &
              ': Bounds of destination mag field: lbnd_destmag = ',           &
              lbnd_destmag, ' ubnd_destmag = ', ubnd_destmag
         write(iulog,"(a,': lon bnd_destmag =',2i4,' gmlon = ',2f9.3)")       &
              subname, lbnd_destmag(1), ubnd_destmag(1),                      &
              gmlon(lbnd_destmag(1)), gmlon(ubnd_destmag(1))
         write(iulog,"(a,': lat bnd_destmag = ',2i4,' gmlat = ',2f9.3)")      &
              subname, lbnd_destmag(2), ubnd_destmag(2),                      &
              gmlat(lbnd_destmag(2)), gmlat(ubnd_destmag(2))
         !
         ! Get 3d bounds of source mag field:
         !
         call ESMF_FieldGet(mag_src_3dfld, localDe=0, farrayPtr=fptr,                &
              computationalLBound=lbnd_srcmag,                                   &
              computationalUBound=ubnd_srcmag, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldGet, mag_phi3d', rc)

         write(iulog,"(a,2(a,i4),a)") subname, ': lon srcmag bounds = (',     &
              lbnd_srcmag(1), ', ', ubnd_srcmag(1), ')'
         write(iulog,"(a,2(a,i4),a)") subname, ': lat srcmag bounds = (',     &
              lbnd_srcmag(2), ', ', ubnd_srcmag(2), ')'
         !
         ! Get 3d bounds of destination geo field:
         !
         call ESMF_FieldGet(geo_3dfld, localDe=0, farrayPtr=fptr,               &
              computationalLBound=lbnd_destgeo,                                  &
              computationalUBound=ubnd_destgeo, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldGet, phys_phi3d', rc)

         write(iulog,"(a,': lon bnd_destgeo=',2i4)") subname,                 &
              lbnd_destgeo(1),ubnd_destgeo(1)
         write(iulog,"(a,': lat bnd_destgeo=',2i4)") subname,                 &
              lbnd_destgeo(2),ubnd_destgeo(2)
      end if
      !
      ! Save route handles for grid transformations in both directions
      ! phys2mag and mag2phys. FieldRegridStore needs to be called only
      ! once for each transformation before the timestep loop (src and
      ! dest fields are still required, so just use ped here). Once inside
      ! the timestep loop, the same routehandle can be used for all fields
      ! that are regridded in the given direction.
      !
      ! These calls will leave *.vtk info files in execdir:
      !   call ESMF_GridWriteVTK(geo_grid, &
      !     staggerloc=ESMF_STAGGERLOC_CENTER, filename="geoGrid",rc=rc)
      !   call ESMF_GridWriteVTK(mag_des_grid, &
      !     staggerloc=ESMF_STAGGERLOC_CENTER, filename="magGrid",rc=rc)

      !
      ! Compute and store route handle for phys2mag 2d fields:
      !
      call ESMF_FieldRegridStore(srcField=phys_2dfld, dstField=mag_des_2dfld, &
           regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
           polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
           extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
           routeHandle=routehandle_phys2mag_2d,                               &
           factorIndexList=factorIndexList,                                   &
           factorList=factorList, srcTermProcessing=smm_srctermproc,          &
           pipelineDepth=smm_pipelinedep, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldRegridStore for 2D phys2mag', rc)

      !
      ! Compute and store route handle for phys2mag 3d fields:
      !
      call ESMF_FieldRegridStore(srcField=phys_3dfld, dstField=mag_des_3dfld, &
           regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
           polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
           extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
           routeHandle=routehandle_phys2mag, factorIndexList=factorIndexList, &
           factorList=factorList, srcTermProcessing=smm_srctermproc,          &
           pipelineDepth=smm_pipelinedep, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldRegridStore for 3D phys2mag', rc)

      !
      ! Save route handle and get esmf indices and weights for mag2phys:
      ! (this overwrites factorIndexList and factorList from above)
      !
      ! These calls will leave *.vtk info files in execdir:
      !   call ESMF_GridWriteVTK(mag_src_grid, &
      !     staggerloc=ESMF_STAGGERLOC_CENTER, filename="magSrcGrid",rc=rc)
      !   call ESMF_GridWriteVTK(geo_des_grid, &
      !     staggerloc=ESMF_STAGGERLOC_CENTER, filename="geoDesGrid",rc=rc)
!!$
!!$      ! Compute and store route handle for mag2phys 3d fields:
!!$      !
!!$      call ESMF_FieldRegridStore(srcField=mag_src_3dfld, dstField=phys_3dfld,     &
!!$           regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
!!$           polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
!!$           routeHandle=routehandle_mag2phys, factorIndexList=factorIndexList, &
!!$           factorList=factorList, srcTermProcessing=smm_srctermproc,          &
!!$           pipelineDepth=smm_pipelinedep, rc=rc)
!!$      call edyn_esmf_chkerr(subname, 'ESMF_FieldRegridStore for 3D mag2phys', rc)

      ! Compute and store route handle for mag2phys 2d (amie) fields:
      call ESMF_FieldRegridStore(srcField=mag2phys_2dfld, dstField=phys_2dfld,&
           regridmethod=ESMF_REGRIDMETHOD_CONSERVE,                           &
           polemethod=ESMF_POLEMETHOD_NONE,                                   &
           routeHandle=routehandle_mag2phys_2d,                               &
           factorIndexList=factorIndexList,                                   &
           factorList=factorList, srcTermProcessing=smm_srctermproc,          &
           pipelineDepth=smm_pipelinedep, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldRegridStore for 2D mag2phys', rc)
      !
      ! Compute and store route handle for phys2geo 3d fields:
      !
      call ESMF_FieldRegridStore(srcField=phys_3dfld, dstField=geo_3dfld,     &
           regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
           polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
           extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
           routeHandle=routehandle_phys2geo, factorIndexList=factorIndexList, &
           factorList=factorList, srcTermProcessing=smm_srctermproc,          &
           pipelineDepth=smm_pipelinedep, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldRegridStore for 3D phys2geo', rc)

      !
      ! Compute and store route handle for geo2phys 3d fields:
      !
      call ESMF_FieldRegridStore(srcField=geo2phys_3dfld, dstField=phys_3dfld,&
           regridmethod=ESMF_REGRIDMETHOD_CONSERVE,                           &
           polemethod=ESMF_POLEMETHOD_NONE,                                   &
           routeHandle=routehandle_geo2phys, factorIndexList=factorIndexList, &
           factorList=factorList, srcTermProcessing=smm_srctermproc,          &
           pipelineDepth=smm_pipelinedep, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldRegridStore for 3D geo2phys', rc)

      !
      ! Compute and store route handle for geo2mag 3d fields:
      !
      call ESMF_FieldRegridStore(srcField=geo_3dfld, dstField=mag_des_3dfld,  &
           regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
           polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
           routeHandle=routehandle_geo2mag, factorIndexList=factorIndexList,  &
           factorList=factorList, srcTermProcessing=smm_srctermproc,          &
           pipelineDepth=smm_pipelinedep, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldRegridStore for 3D mag2geo', rc)
      !
      ! Compute and store route handle for geo2mag 2d fields:
      !
      call ESMF_FieldRegridStore(srcField=geo_2dfld, dstField=mag_des_2dfld,  &
           regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
           polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
           routeHandle=routehandle_geo2mag_2d,                                &
           factorIndexList=factorIndexList,                                   &
           factorList=factorList, srcTermProcessing=smm_srctermproc,          &
           pipelineDepth=smm_pipelinedep, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldRegridStore for 2D geo2mag', rc)

      !
      ! Compute and store route handle for mag2geo 3d fields:
      !
      call ESMF_FieldRegridStore(srcField=mag_src_3dfld, dstField=geo_3dfld,  &
           regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
           polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
           routeHandle=routehandle_mag2geo, factorIndexList=factorIndexList,  &
           factorList=factorList, srcTermProcessing=smm_srctermproc,          &
           pipelineDepth=smm_pipelinedep, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldRegridStore for 3D mag2geo', rc)

      edyn_esmf_update_step = .true.
#endif
   end subroutine edyn_esmf_update

#ifdef WACCMX_EDYN_ESMF
   !-----------------------------------------------------------------------
   subroutine create_mag_grid(grid_out, srcdes)
      !
      ! Create ESMF geomagnetic grid, w/ lon,lat coordinates.
      !
      ! Args:
      type(ESMF_Grid),  intent(out) :: grid_out
      character(len=*), intent(in)  :: srcdes
      !
      ! Local:
      integer                     :: i,j,n,rc
      real(ESMF_KIND_R8), pointer :: coordX(:,:),coordY(:,:)
      integer                     :: lbnd(2),ubnd(2)
      integer                     :: nmlons_task(ntaski) ! number of lons per task
      integer                     :: nmlats_task(ntaskj) ! number of lats per task
      character(len=*), parameter :: subname = 'create_mag_grid'
      !
      ! We are creating either a source grid or a destination grid:
      !
      if (srcdes /= 'src' .and. srcdes /= 'des') then
         write(iulog,"(a)") '>>> create_mag_grid: srcdes = ''',srcdes, &
              ''' but must be either ''src'' or ''des'''
         call endrun('create_mag_grid: srcdes')
      end if
      !
      ! nmlons_task(nmagtaski) = number of mag lons per task in lon dim
      !
      do i = 1, nmagtaski
         loop: do n = 0, ntask - 1
            if (tasks(n)%magtidi == i-1) then
               nmlons_task(i) = tasks(n)%nmaglons
               exit loop
            end if
         end do loop
      end do
      !
      ! Exclude periodic points (1 point fewer for mpi tasks at east end)
      ! for source grids (this overwrites above for eastern-most tasks):
      !
      if (srcdes == 'src') then
         do n = 0, ntask-1
            if (tasks(n)%magtidi == nmagtaski-1) then ! east edge of proc matrix
               nmlons_task(tasks(n)%magtidi+1) = tasks(n)%nmaglons-1
            end if
         end do
      end if
      !
      ! nmlats_task(nmagtaskj) = number of mag lats per task in lat dim
      !
      do j = 1, nmagtaskj
         loop1: do n = 0, ntask-1
            if (tasks(n)%magtidj == j-1) then
               nmlats_task(j) = tasks(n)%nmaglats
               exit loop1
            end if
         end do loop1
      end do
      !
      ! Create curvilinear magnetic grid (both coords depend
      ! on both dimensions, i.e., lon(i,j),lat(i,j)):
      !
      grid_out = ESMF_GridCreate1PeriDim(                  &
           countsPerDEDim1=nmlons_task, coordDep1=(/1,2/), &
           countsPerDEDim2=nmlats_task, coordDep2=(/1,2/), &
           indexflag=ESMF_INDEX_GLOBAL,rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCreate1PeriDim', rc)
      !
      ! Allocate coordinates:
      !
      call ESMF_GridAddCoord(grid_out,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridAddCoord', rc)
      !
      ! Get pointer and set mag grid longitude coordinates:
      !
      call ESMF_GridGetCoord(grid_out, coordDim=1, localDE=0,                 &
           computationalLBound=lbnd, computationalUBound=ubnd,                &
           farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridGetCoord for longitude coords', rc)

      do j = lbnd(2), ubnd(2)
         do i = lbnd(1), ubnd(1)
            coordX(i,j) = gdlondeg(i,j)
         end do
      end do
      !
      ! Get pointer and set mag grid latitude coordinates:
      !
      call ESMF_GridGetCoord(grid_out, coordDim=2, localDE=0,  &
           computationalLBound=lbnd, computationalUBound=ubnd, &
           farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridGetCoord for latitude coords', rc)

      do j = lbnd(2), ubnd(2)
         do i = lbnd(1), ubnd(1)
            coordY(i,j) = gdlatdeg(i,j)
         end do
      end do

      if (debug .and. masterproc) then
         write(iulog,"(3a,2i4,a,2i4,a,2i4,a,2i4)") 'Created ESMF ',srcdes,    &
              ' mag grid:  lbnd, ubnd_lon=', lbnd(1), ubnd(1),                &
              ' mlon0,1=', mlon0, mlon1, ' lbnd, ubnd_lat=',                  &
              lbnd(2), ubnd(2), ' mlat0,1=', mlat0, mlat1
      end if

   end subroutine create_mag_grid

   !-----------------------------------------------------------------------
   subroutine create_mag2phys_grid(grid_out)
     use edyn_params, only: rtd
     use edyn_maggrid,only: dlonm, nmlat
     !
     ! Create ESMF geomagnetic grid, w/ lon,lat coordinates.
     !
     ! Args:
     type(ESMF_Grid),  intent(out) :: grid_out
     !
     ! Local:
     integer                     :: i,j,n,rc
     real(ESMF_KIND_R8), pointer :: coordX(:,:),coordY(:,:)
     integer                     :: lbnd(2),ubnd(2)
     integer                     :: nmlons_task(ntaski) ! number of lons per task
     integer                     :: nmlats_task(ntaskj) ! number of lats per task
     character(len=*), parameter :: subname = 'create_mag2phys_grid'

     !
     ! nmlons_task(nmagtaski) = number of mag lons per task in lon dim
     !
     do i = 1, nmagtaski
        loop: do n = 0, ntask - 1
           if (tasks(n)%magtidi == i-1) then
              nmlons_task(i) = tasks(n)%nmaglons
              exit loop
           end if
        end do loop
     end do
     !
     ! Exclude periodic points (1 point fewer for mpi tasks at east end)
     ! for source grids (this overwrites above for eastern-most tasks):
     !
     do n = 0, ntask-1
        if (tasks(n)%magtidi == nmagtaski-1) then ! east edge of proc matrix
           nmlons_task(tasks(n)%magtidi+1) = tasks(n)%nmaglons-1
        end if
     end do

     !
     ! nmlats_task(nmagtaskj) = number of mag lats per task in lat dim
     !
     do j = 1, nmagtaskj
        loop1: do n = 0, ntask-1
           if (tasks(n)%magtidj == j-1) then
              nmlats_task(j) = tasks(n)%nmaglats
              exit loop1
           end if
        end do loop1
     end do

     !
     ! Create curvilinear magnetic grid (both coords depend
     ! on both dimensions, i.e., lon(i,j),lat(i,j)):
     !
     grid_out = ESMF_GridCreate1PeriDim(                  &
          countsPerDEDim1=nmlons_task, coordDep1=(/1,2/), &
          countsPerDEDim2=nmlats_task, coordDep2=(/1,2/), &
          indexflag=ESMF_INDEX_GLOBAL,rc=rc)
     call edyn_esmf_chkerr(subname, 'ESMF_GridCreate1PeriDim', rc)
     !
     ! Allocate coordinates:
     !
     call ESMF_GridAddCoord(grid_out,staggerloc=ESMF_STAGGERLOC_CORNER,rc=rc)
     call edyn_esmf_chkerr(subname, 'ESMF_GridAddCoord', rc)
     !
     ! Get pointer and set mag grid longitude coordinates:
     !
     call ESMF_GridGetCoord(grid_out, coordDim=1, localDE=0,                 &
          computationalLBound=lbnd, computationalUBound=ubnd,                &
          farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
     call edyn_esmf_chkerr(subname, 'ESMF_GridGetCoord for longitude coords', rc)

     do j = lbnd(2), ubnd(2)
        do i = lbnd(1), ubnd(1)
           coordX(i,j) = gdlondeg_corners(i,j)
        end do
     end do
     !
     ! Get pointer and set mag grid latitude coordinates:
     !
     call ESMF_GridGetCoord(grid_out, coordDim=2, localDE=0,  &
          computationalLBound=lbnd, computationalUBound=ubnd, &
          farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
     call edyn_esmf_chkerr(subname, 'ESMF_GridGetCoord for latitude coords', rc)

     do j = lbnd(2), ubnd(2)
        do i = lbnd(1), ubnd(1)
           coordY(i,j) = gdlatdeg_corners(i,j)
        end do
     end do

     if (debug .and. masterproc) then
        write(iulog,"(a,2i4,a,2i4,a,2i4,a,2i4)") &
             'Created ESMF mag grid: lbnd,ubnd_lon=', lbnd(1), ubnd(1), &
             ' mlon0,1=', mlon0, mlon1, ' lbnd, ubnd_lat=',&
             lbnd(2), ubnd(2), ' mlat0,1=', mlat0, mlat1
     end if

   end subroutine create_mag2phys_grid

   !-----------------------------------------------------------------------
   subroutine create_geo_grid(grid_out)
      !
      ! Args:
      type(ESMF_Grid),intent(out) :: grid_out
      !
      ! Local:
      integer                       :: i, j, n, rc
      integer                       :: lbnd_lat, ubnd_lat, lbnd_lon, ubnd_lon
      integer                       :: lbnd(1), ubnd(1)
      integer                       :: nlons_task(ntaski) ! # lons per task
      integer                       :: nlats_task(ntaskj) ! # lats per task
      real(ESMF_KIND_R8), pointer   :: coordX(:), coordY(:)
      character(len=*),   parameter :: subname = 'create_geo_grid'
      !
      ! nlons_task(ntaski) = number of geo lons per task.
      !
      do i = 1, ntaski
         loop: do n = 0, ntask-1
            if (tasks(n)%mytidi == i-1) then
               nlons_task(i) = tasks(n)%nlons
               exit loop
            end if
         end do loop
      end do
      !
      do j = 1, ntaskj
         loop1: do n = 0, ntask-1
            if (tasks(n)%mytidj == j-1) then
               nlats_task(j) = tasks(n)%nlats
               exit loop1
            end if
         end do loop1
      end do
      !
      !
      ! Create 2d geographic source grid (with poles)
      grid_out = ESMF_GridCreate1PeriDim(                       &
           countsPerDEDim1=nlons_task, coordDep1=(/1/),         &
           countsPerDEDim2=nlats_task, coordDep2=(/2/),         &
           indexflag=ESMF_INDEX_GLOBAL,minIndex=(/1,1/), rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCreate1PeriDim', rc)
      !
      ! Allocate coordinates:
      !
      call ESMF_GridAddCoord(grid_out, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridAddCoord', rc)
      !
      ! Get pointer and set geo grid longitude coordinates:
      !
      call ESMF_GridGetCoord(grid_out, coordDim=1, localDE=0,                 &
           computationalLBound=lbnd, computationalUBound=ubnd,                &
           farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridGetCoord for longitude coords', rc)
      !
      ! Note glon was shifted to +/-180 by sub set_geogrid (edyn_init.F90)
      !
      lbnd_lon = lbnd(1)
      ubnd_lon = ubnd(1)
      do i = lbnd_lon, ubnd_lon
         coordX(i) = glon(i)          ! 1 -> 72
      end do
      !
      ! Get pointer and set geo grid latitude coordinates, including poles:
      !
      call ESMF_GridGetCoord(grid_out, coordDim=2, localDE=0,                 &
           computationalLBound=lbnd, computationalUBound=ubnd,                &
           farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridGetCoord for latitude coords', rc)

      lbnd_lat = lbnd(1)
      ubnd_lat = ubnd(1)
      do i = lbnd_lat, ubnd_lat
         coordY(i) = glat(i)
      end do

      if (debug .and. masterproc) then
         write(iulog,"(2a,2i4,a,2i4,a,2i4,a,2i4)") 'Created ESMF geo_grid:',  &
              ' lbnd,ubnd_lon=', lbnd_lon, ubnd_lon, ' lon0,1=', lon0, lon1,  &
              ' lbnd,ubnd_lat=', lbnd_lat, ubnd_lat, ' lat0,1=', lat0, lat1
         write(iulog,"('coordX for geo grid = ',/,(8f10.4))") coordX
         write(iulog,"('coordY for geo grid = ',/,(8f10.4))") coordY
      end if

   end subroutine create_geo_grid

   !-----------------------------------------------------------------------
   subroutine create_geo2phys_grid(grid_out)
      use edyn_geogrid, only: glat_corner,glon_corner
      !
      ! Args:
      type(ESMF_Grid),intent(out) :: grid_out
      !
      ! Local:
      integer                       :: i, j, n, rc
      integer                       :: lbnd_lat, ubnd_lat, lbnd_lon, ubnd_lon
      integer                       :: lbnd(1), ubnd(1)
      integer                       :: nlons_task(ntaski) ! # lons per task
      integer                       :: nlats_task(ntaskj) ! # lats per task
      real(ESMF_KIND_R8), pointer   :: coordX(:), coordY(:)
      character(len=*),   parameter :: subname = 'create_geo_grid'
      !
      ! nlons_task(ntaski) = number of geo lons per task.
      !
      do i = 1, ntaski
         loop: do n = 0, ntask-1
            if (tasks(n)%mytidi == i-1) then
               nlons_task(i) = tasks(n)%nlons
               exit loop
            end if
         end do loop
      end do
      !
      do j = 1, ntaskj
         loop1: do n = 0, ntask-1
            if (tasks(n)%mytidj == j-1) then
               nlats_task(j) = tasks(n)%nlats
               exit loop1
            end if
         end do loop1
      end do
      !
      !
      ! Create 2d geographic source grid (with poles)
      grid_out = ESMF_GridCreate1PeriDim(                       &
           countsPerDEDim1=nlons_task, coordDep1=(/1/),         &
           countsPerDEDim2=nlats_task, coordDep2=(/2/),         &
           indexflag=ESMF_INDEX_GLOBAL,minIndex=(/1,1/), rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridCreate1PeriDim', rc)
      !
      ! Allocate coordinates:
      !
      call ESMF_GridAddCoord(grid_out, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridAddCoord', rc)
      !
      ! Get pointer and set geo grid longitude coordinates:
      !
      call ESMF_GridGetCoord(grid_out, coordDim=1, localDE=0,                 &
           computationalLBound=lbnd, computationalUBound=ubnd,                &
           farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridGetCoord for longitude coords', rc)
      !
      ! Note glon was shifted to +/-180 by sub set_geogrid (edyn_init.F90)
      !
      lbnd_lon = lbnd(1)
      ubnd_lon = ubnd(1)
      do i = lbnd_lon, ubnd_lon
         coordX(i) = glon_corner(i)
      end do
      !
      ! Get pointer and set geo grid latitude coordinates, including poles:
      !
      call ESMF_GridGetCoord(grid_out, coordDim=2, localDE=0,                 &
           computationalLBound=lbnd, computationalUBound=ubnd,                &
           farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridGetCoord for latitude coords', rc)

      lbnd_lat = lbnd(1)
      ubnd_lat = ubnd(1)
      do i = lbnd_lat, ubnd_lat
         coordY(i) = glat_corner(i)
      end do

      if (debug .and. masterproc) then
         write(iulog,"(2a,2i4,a,2i4,a,2i4,a,2i4)") 'Created ESMF geo2phys_grid:',  &
              ' lbnd,ubnd_lon=', lbnd_lon, ubnd_lon, ' lon0,1=', lon0, lon1,  &
              ' lbnd,ubnd_lat=', lbnd_lat, ubnd_lat, ' lat0,1=', lat0, lat1
         write(iulog,"('coordX for geo2phys grid = ',/,(8f10.4))") coordX
         write(iulog,"('coordY for geo2phys grid = ',/,(8f10.4))") coordY
      end if
   end subroutine create_geo2phys_grid
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_create_physfield(field, mesh, name, nlev)
      use ESMF, only: ESMF_FieldStatus_Flag, ESMF_FIELDSTATUS_EMPTY, ESMF_MESHLOC_ELEMENT, ESMF_SUCCESS
      use ESMF, only: operator(/=)
      use ESMF, only: operator(==)
      !
      ! Create ESMF field (2d or 3d) on physics mesh
      ! If nlev == 0, field is 2d (i,j), otherwise field is 3d,
      !   and 3rd dimension is ungridded
      !
      ! Args:
      integer,          intent(in)  :: nlev ! if nlev == 0, field is 2d (i,j)
      type(ESMF_Mesh),  intent(in)  :: mesh
      character(len=*), intent(in)  :: name
      type(ESMF_Field), intent(out) :: field
      !
      ! Local:
      integer                     :: rc
      type(ESMF_ArraySpec)        :: arrayspec
      character(len=*), parameter :: subname = 'edyn_esmf_create_physfield'
#if 0
      type(ESMF_FieldStatus_Flag) :: fstatus
      character(len=256)          :: errmsg
      !
      errmsg = ''
 ! this call to FieldGet is not successful when field is not instantiated
      call ESMF_FieldGet(field, status=fstatus, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)
      if (fstatus /= ESMF_FIELDSTATUS_EMPTY) then
         write(errmsg, '(a)') trim(subname) // ': Field, ' // trim(name) // ', is not empty'
         if (masterproc) then
            write(iulog, '(a)') 'ERROR: ' // trim(errmsg)
         end if
         call endrun(trim(errmsg))
      end if
#endif
      ! Create 3d field (i,j,k), with non-distributed vertical dimension:
      if (nlev > 0) then
         call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_ArraySpecSet 2D', rc)
         field = ESMF_FieldCreate(mesh, arrayspec, &
              gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT,  &
              ungriddedLBound=(/1/), ungriddedUBound=(/nlev/), rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldCreate 2D field', rc)
         !
         ! Create 2d field (i,j):
      else                ! create 2d field
         call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_ArraySpecSet 2D', rc)
         field = ESMF_FieldCreate(mesh, arrayspec,  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldCreate 2D field', rc)
      end if

    end subroutine edyn_esmf_create_physfield
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_create_geofield(field, grid, name, nlev)
      use ESMF, only: ESMF_FieldStatus_Flag, ESMF_FIELDSTATUS_EMPTY
      use ESMF, only: operator(/=)
      !
      ! Create ESMF field (2d or 3d) on geo grid (will exclude periodic points)
      ! If nlev == 0, field is 2d (i,j), otherwise field is 3d,
      !   and 3rd dimension is ungridded
      !
      ! Args:
      integer,          intent(in)  :: nlev ! if nlev == 0, field is 2d (i,j)
      type(ESMF_Grid),  intent(in)  :: grid
      character(len=*), intent(in)  :: name
      type(ESMF_Field), intent(out) :: field
      !
      ! Local:
      integer                     :: rc
      type(ESMF_ArraySpec)        :: arrayspec
      character(len=*), parameter :: subname = 'edyn_esmf_create_geofield'
      !
#if 0
      type(ESMF_FieldStatus_Flag) :: fstatus
      character(len=256)          :: errmsg
      errmsg = ''
      call ESMF_FieldGet(field, status=fstatus, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)
      if (fstatus /= ESMF_FIELDSTATUS_EMPTY) then
         write(errmsg, '(3a)') subname, ': Field, ',trim(name),', is not empty'
         if (masterproc) then
            write(iulog, '(2a)') 'ERROR: ', trim(errmsg)
         end if
         call endrun(trim(errmsg))
      end if
#endif
      ! Create 3d field (i,j,k), with non-distributed vertical dimension:
      if (nlev > 0) then
         call ESMF_ArraySpecSet(arrayspec,3,ESMF_TYPEKIND_R8,rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_ArraySpecSet 3D', rc)
         field = ESMF_FieldCreate(grid, arrayspec,ungriddedLBound=(/1/),      &
              ungriddedUBound=(/nlev/),staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldCreate 3D field', rc)
         !
         ! Create 2d field (i,j):
      else                ! create 2d field
         call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_ArraySpecSet 2D', rc)
         field = ESMF_FieldCreate(grid, arrayspec,&
              staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldCreate 2D field', rc)
      end if
   end subroutine edyn_esmf_create_geofield
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_create_magfield(field, grid, name, nlev)
      !
      ! Create ESMF field (2d or 3d) on mag grid. This will include the
      !   mag periodic point, which will be zero after regridding.
      ! If nlev == 0, field is 2d (i,j), otherwise field is 3d,
      !   and 3rd dimension is ungridded
      !
      ! Args:
      integer,          intent(in)  :: nlev ! if nlev == 0, field is 2d (i,j)
      type(ESMF_Grid),  intent(in)  :: grid
      character(len=*), intent(in)  :: name
      type(ESMF_Field), intent(out) :: field
      !
      ! Local:
      integer                     :: rc
      type(ESMF_ArraySpec)        :: arrayspec
      type(ESMF_Array)            :: array3d,array2d
      type(ESMF_DistGrid)         :: distgrid
      character(len=*), parameter :: subname = 'edyn_esmf_create_magfield'
      !
      ! Get necessary information from the mag grid:
      call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CENTER, &
           distgrid=distgrid,rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_GridGet', rc)
      !
      ! Create 3d mag field (i,j,k), with non-distributed vertical dimension:
      ! (add periodic point in longitude with computationalEdgeUWidth)
      !
      if (nlev > 0) then
         call ESMF_ArraySpecSet(arrayspec,3,ESMF_TYPEKIND_R8,rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_ArraySpecSet 3D field', rc)

         array3d = ESMF_ArrayCreate(arrayspec=arrayspec,                      &
              distgrid=distgrid,computationalEdgeUWidth=(/1,0/),              &
              undistLBound=(/1/),undistUBound=(/nlev/),                       &
              indexflag=ESMF_INDEX_GLOBAL,rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_ArrayCreate 3D field', rc)

         field = ESMF_FieldCreate(grid, array3d,                              &
              ungriddedLBound=(/1/), ungriddedUBound=(/nlev/),                &
              staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldCreate 3D field', rc)
         !
         ! Create 2d mag field (i,j):
         ! (add periodic point in longitude with computationalEdgeUWidth)
         !
      else                ! create 2d field
         call ESMF_ArraySpecSet(arrayspec,2,ESMF_TYPEKIND_R8,rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_ArraySpecSet 2D field', rc)

         array2d = ESMF_ArrayCreate(arrayspec=arrayspec,                      &
              distgrid=distgrid,computationalEdgeUWidth=(/1,0/),              &
              indexflag=ESMF_INDEX_GLOBAL,rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_ArrayCreate 2D field', rc)
         field = ESMF_FieldCreate(grid, array2d,                              &
              staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldCreate 2D field', rc)
      end if
   end subroutine edyn_esmf_create_magfield

   !-----------------------------------------------------------------------
   subroutine edyn_esmf_set3d_geo(field, fdata, ilon0, ilon1, ilat0, ilat1, ilev0, ilev1 )
      !
      ! Set values of a 3d ESMF field on geographic source grid, prior to
      !   geographic to physics grid transformation.
      ! Periodic points are excluded, geographic poles are at
      !   j==jspole and jnpole
      ! Note dimension order changes from input (k,i,j) to output (i,j,k).
      !
      ! Args:
      type(ESMF_Field), intent(in) :: field ! esmf fields on geo grid
      !
      ! field is input data on model subdomains (including periodic points)
      ! (note esmf source field excludes periodic points)
      !
      integer,  intent(in) :: ilev0, ilev1, ilon0, ilon1, ilat0, ilat1
      real(r8), intent(in) :: fdata(ilon0:ilon1,ilat0:ilat1,ilev0:ilev1)
      !
      ! Local:
      integer                         :: i, j, k, rc
      integer                         :: lbnd(3), ubnd(3) ! 3d field bounds
      !
      ! fptr is esmf      pointer (i,j,k) to 3d field, set by this subroutine
      real(ESMF_KIND_R8), pointer     :: fptr(:,:,:)
      character(len=*),   parameter   :: subname = 'edyn_esmf_set3d_geo'

      !
      ! Get and set pointer to the field:
      call ESMF_FieldGet(field, localDe=0, farrayPtr=fptr,        &
           computationalLBound=lbnd, computationalUBound=ubnd, rc=rc)
      call edyn_esmf_chkerr(subname, 'ESMF_FieldGet field', rc)

      fptr = nan
      do k = lbnd(3),ubnd(3)
         do j = lbnd(2),ubnd(2)
            do i = lbnd(1),ubnd(1)
               fptr(i,j,k) = fdata(i,j,k)
            end do
         end do
      end do

   end subroutine edyn_esmf_set3d_geo
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_set2d_geo(field, fdata, ilon0, ilon1, ilat0, ilat1)
     !
     ! Set values of a 2d ESMF field on geographic source grid, prior to
     !   geographic to physics grid transformation. (Essentially the same
     !   as esmf_set3d_geo, except for 2d fields instead of 3d)
     ! Periodic points are excluded, geographic poles are at j==jspole and jnpole
     !
     ! Args:
     type(ESMF_Field), intent(in) :: field
     integer,          intent(in) :: ilon0,ilon1,ilat0,ilat1
     real(r8),         intent(in) :: fdata(ilon0:ilon1,ilat0:ilat1)
     !
     ! Local:
     integer                       :: i, j, rc
     real(ESMF_KIND_R8), pointer   :: fptr(:,:)
     integer                       :: lbnd(2), ubnd(2)
     character(len=*),   parameter :: subname = 'edyn_esmf_set2d_geo'
     !
     ! Get pointer to the field:
     call ESMF_FieldGet(field, localDe=0, farrayPtr=fptr,                    &
          computationalLBound=lbnd, computationalUBound=ubnd, rc=rc)
     call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)
     !
     fptr(:,:) = nan
     !
     ! Set interior latitudes (excluding poles):
     do j = lbnd(2), ubnd(2)
        do i = lbnd(1), ubnd(1)
           fptr(i,j) = fdata(i,j)
        end do
     end do

   end subroutine edyn_esmf_set2d_geo

   !-----------------------------------------------------------------------
   subroutine edyn_esmf_set3d_mag(field, fdata, ilon0, ilon1, ilat0, ilat1, ilev0, ilev1 )
     !
     ! Set values of a 3d ESMF field on magnetic grid, prior to magnetic to
     ! physics grid transformation.
     !
     ! Args:
     type(ESMF_Field), intent(in) :: field ! esmf fields on mag grid
     !
     ! fdata is input data on model subdomains:
     !
     integer,          intent(in) :: ilev0, ilev1, ilon0, ilon1, ilat0, ilat1
     real(r8),         intent(in) :: fdata(ilon0:ilon1,ilat0:ilat1,ilev0:ilev1)
     !
     ! Local:
     integer                       :: i, j, k, rc
     integer                       :: lbnd(3), ubnd(3) ! 3d field bounds
     !
     ! fptr is esmf      pointer (i,j,k) to 3d field, set by this subroutine
     real(ESMF_KIND_R8), pointer   :: fptr(:,:,:)
     character(len=*),   parameter :: subname = 'edyn_esmf_set3d_mag'

     call ESMF_FieldGet(field, localDe=0, farrayPtr=fptr,              &
          computationalLBound=lbnd, computationalUBound=ubnd, rc=rc)
     call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)
     !
     fptr(:,:,:) = nan
     !
     ! Set ESMF pointer:
     !
     do k = lbnd(3), ubnd(3)
        do j = lbnd(2), ubnd(2)
           do i = lbnd(1), ubnd(1)
              fptr(i,j,k) = fdata(i,j,k)
           end do
        end do
     end do
   end subroutine edyn_esmf_set3d_mag
   !-----------------------------------------------------------------------

   subroutine edyn_esmf_set2d_mag(field, fdata, ilon0, ilon1, ilat0, ilat1)
     !
     ! Set values of a 2d ESMF field on magnetic grid, prior to magnetic to
     ! physics grid transformation.
     !
     ! Args:
     type(ESMF_Field), intent(in) :: field ! esmf field on mag grid
     !
     ! fdata is input data on model subdomains:
     !
     integer,          intent(in) :: ilon0, ilon1, ilat0, ilat1
     real(r8),         intent(in) :: fdata(ilon0:ilon1,ilat0:ilat1)
     !
     ! Local:
     integer                       :: i, j, rc
     integer                       :: lbnd(2), ubnd(2) ! 2d field bounds
     !
     ! fptr is esmf      pointer (i,j,k) to 2d field, set by this subroutine
     real(ESMF_KIND_R8), pointer   :: fptr(:,:)
     character(len=*),   parameter :: subname = 'edyn_esmf_set2d_mag'
     !
     call ESMF_FieldGet(field, localDe=0, farrayPtr=fptr,             &
          computationalLBound=lbnd, computationalUBound=ubnd, rc=rc)
     call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)
     !
     fptr(:,:) = nan
     !
     ! Set ESMF pointer:
     !
     do j = lbnd(2), ubnd(2)      ! lat
        do i = lbnd(1), ubnd(1)    ! lon
           fptr(i,j) = fdata(i,j)
        end do   ! mlon
     end do     ! mlat
     !
  end subroutine edyn_esmf_set2d_mag
  !-----------------------------------------------------------------------
  subroutine edyn_esmf_set3d_phys(field, fdata, ilev0, ilev1, icol0, icol1)
    !
    ! Set values ESMF field on physics mesh
    !
    ! Args:
    type(ESMF_Field), intent(in) :: field ! esmf field on phys grid
    !
    ! fdata is input data on model subdomains:
    !
    integer,          intent(in) :: ilev0, ilev1, icol0, icol1
    real(r8),         intent(in) :: fdata(ilev0:ilev1,icol0:icol1)
    !
    ! Local:
    integer                       :: i, k, rc
    integer                       :: lbnd(2), ubnd(2) ! 2d field bounds
    !
    ! fptr is esmf      pointer to 2d field, set by this subroutine
    real(ESMF_KIND_R8), pointer   :: fptr(:,:)
    character(len=*),   parameter :: subname = 'edyn_esmf_set2d_phys'
    !
    ! Fields loop:
    call ESMF_FieldGet(field, localDe=0, farrayPtr=fptr, &
         computationalLBound=lbnd, computationalUBound=ubnd, &
         rc=rc)
    call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)
    !
    fptr(:,:) = nan
    !
    ! Set ESMF pointer:
    !
    do i = lbnd(2), ubnd(2)
       do k = lbnd(1), ubnd(1)
          fptr(k,i) = fdata(k,i)
       end do
    end do
  end subroutine edyn_esmf_set3d_phys
  !-----------------------------------------------------------------------
  !
  subroutine edyn_esmf_set2d_phys(field, fdata, icol0, icol1)
    !
    ! Set values of a ESMF field on phys mesh
    !
    ! Args:
    type(ESMF_Field), intent(in) :: field  ! esmf field on mag grid
    !
    ! fdata is input data on model subdomains:
    !
    integer,          intent(in) :: icol0, icol1
    real(r8),         intent(in) :: fdata(icol0:icol1)
    !
    ! Local:
    integer                       :: i, rc
    integer                       :: lbnd(1), ubnd(1) ! 2d field bounds
    !
    ! fptr is esmf  1d field, set by this subroutine
    real(ESMF_KIND_R8), pointer   :: fptr(:)
    character(len=*),   parameter :: subname = 'edyn_esmf_set1d_phys'
    !
    call ESMF_FieldGet(field, localDe=0, farrayPtr=fptr,             &
         computationalLBound=lbnd, computationalUBound=ubnd, rc=rc)
    call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)
    !
    fptr(:) = nan
    !
    ! Set ESMF pointer:
    !
    do i = lbnd(1), ubnd(1)
       fptr(i) = fdata(i)
    end do
    !
  end subroutine edyn_esmf_set2d_phys
  !-----------------------------------------------------------------------
  subroutine edyn_esmf_get_3dfield(field, data, i0,i1,j0,j1,k0,k1 )
    !
    ! Get pointer to 3d esmf field (i,j,k):
    !
    ! Args:
    integer,           intent(in)  :: i0,i1,j0,j1,k0,k1
    type(ESMF_field),  intent(in)  :: field
    real(r8),          intent(out) :: data(i0:i1,j0:j1,k0:k1)
    !
    ! Local:
    real(r8), pointer           :: fptr(:,:,:)
    integer                     :: rc, lbnd(3), ubnd(3)
    character(len=*), parameter :: subname = 'edyn_esmf_get_3dfield'

    call ESMF_FieldGet(field, localDe=0, farrayPtr=fptr,                   &
         computationalLBound=lbnd, computationalUBound=ubnd, rc=rc)
    call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)

    data(:,:,:) = fptr(:,:,:)
  end subroutine edyn_esmf_get_3dfield
  !-----------------------------------------------------------------------
  subroutine edyn_esmf_get_2dfield(field, data, i0,i1,j0,j1 )
    !
    ! Get pointer to 2d esmf field (i,j):
    !
    ! Args:
    integer,           intent(in)  :: i0,i1,j0,j1
    type(ESMF_field),  intent(in)  :: field
    real(r8),          intent(out) :: data(i0:i1,j0:j1)
    !
    ! Local:
    real(r8), pointer :: fptr(:,:)
    integer :: rc
    character(len=*),   parameter :: subname = 'edyn_esmf_get_2dfield'

    call ESMF_FieldGet(field, localDe=0, farrayPtr=fptr, rc=rc)
    call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)

    data(:,:) = fptr(:,:)

  end subroutine edyn_esmf_get_2dfield
  !-----------------------------------------------------------------------
  subroutine edyn_esmf_get_1dfield(field, data, i0,i1 )
    !
    ! Get pointer to 2d esmf field (i,j):
    !
    ! Args:
    integer,           intent(in)  :: i0,i1
    type(ESMF_field),  intent(in)  :: field
    real(r8),          intent(out) :: data(i0:i1)
    !
    ! Local:
    real(r8), pointer :: fptr(:)
    integer :: rc
    character(len=*),   parameter :: subname = 'edyn_esmf_get_2dfield'

    call ESMF_FieldGet(field, localDe=0, farrayPtr=fptr, rc=rc)
    call edyn_esmf_chkerr(subname, 'ESMF_FieldGet', rc)

    data(:) = fptr(:)

  end subroutine edyn_esmf_get_1dfield
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_regrid_phys2mag(srcfield, dstfield, ndim)
      !
      ! Args:
      integer                         :: ndim
      type(ESMF_Field), intent(inout) :: srcfield, dstfield
      !
      ! Local:
      integer                     :: rc
      character(len=*), parameter :: subname = 'edyn_esmf_regrid_phys2mag'
      !
      if (ndim == 2) then
         !
         ! Do sparse matrix multiply for 2d phys2mag.
         !
         call ESMF_FieldRegrid(srcfield, dstfield, routehandle_phys2mag_2d,      &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid phys2mag 2D', rc)
      else ! 3d geo2mag
         !
         ! Do sparse matrix multiply for 3d geo2mag.
         !
         call ESMF_FieldRegrid(srcfield, dstfield, routehandle_phys2mag,         &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid phys2mag 3D', rc)
      end if
   end subroutine edyn_esmf_regrid_phys2mag
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_regrid_mag2phys(srcfield, dstfield, ndim)
      !
      ! Args:
      type(ESMF_Field), intent(inout) :: srcfield, dstfield
      integer                         :: ndim
      !
      ! Local:
      integer                     :: rc
      character(len=*), parameter :: subname = 'edyn_esmf_regrid_mag2phys'
      !
      if (ndim == 2) then
         call ESMF_FieldRegrid(srcfield, dstfield, routehandle_mag2phys_2d, &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid mag2phys 2D', rc)
      else
!!$         call ESMF_FieldRegrid(srcfield, dstfield, routehandle_mag2phys, &
!!$              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
!!$         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid mag2phys 3D', rc)
         call endrun('edyn_esmf_regrid_mag2phys: no 3D routehandle')
      end if
   end subroutine edyn_esmf_regrid_mag2phys
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_regrid_phys2geo(srcfield, dstfield, ndim)
      !
      ! Args:
      integer                         :: ndim
      type(ESMF_Field), intent(inout) :: srcfield, dstfield
      !
      ! Local:
      integer                     :: rc
      character(len=*), parameter :: subname = 'edyn_esmf_regrid_phys2geo'
      !
      if (ndim == 2) then
         !
         ! Do sparse matrix multiply for 2d phys2mag.
         !
!!$         call ESMF_FieldRegrid( srcfield, dstfield, routehandle_phys2geo_2d,         &
!!$              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
!!$         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid phys2geo 2D', rc)
         call endrun('edyn_esmf_regrid_phys2geo: 2D not working')
      else ! 3d phys2geo
         !
         ! Do sparse matrix multiply for 3d phys2geo.
         !
         call ESMF_FieldRegrid( srcfield, dstfield, routehandle_phys2geo,         &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid phys2geo 3D', rc)
      end if
   end subroutine edyn_esmf_regrid_phys2geo
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_regrid_geo2phys(srcfield, dstfield, ndim)
      !
      ! Args:
      integer                         :: ndim
      type(ESMF_Field), intent(inout) :: srcfield, dstfield
      !
      ! Local:
      integer                     :: rc
      character(len=*), parameter :: subname = 'edyn_esmf_regrid_geo2phys'
      !
      if (ndim == 2) then
!!$         call ESMF_FieldRegrid(srcfield, dstfield, routehandle_geo2phys_2d,      &
!!$           termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
!!$         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid geo2phys 2D', rc)
         call endrun('edyn_esmf_regrid_geo2phys: 2D not working')
      else
         call ESMF_FieldRegrid( srcfield, dstfield, routehandle_geo2phys,        &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid geo2phys 3D', rc)
      end if
   end subroutine edyn_esmf_regrid_geo2phys
   !-----------------------------------------------------------------------
   subroutine edyn_esmf_regrid_geo2mag(srcfield, dstfield, ndim)
      !
      ! Args:
      integer                         :: ndim
      type(ESMF_Field), intent(inout) :: srcfield, dstfield
      !
      ! Local:
      integer                     :: rc
      character(len=*), parameter :: subname = 'edyn_esmf_regrid_geo2mag'
      !
      if (ndim == 2) then
         !
         ! Do sparse matrix multiply for 2d geo2mag.
         !
         call ESMF_FieldRegrid(srcfield, dstfield, routehandle_geo2mag_2d,       &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid geo2mag 2D', rc)
      else ! 3d geo2mag
         !
         ! Do sparse matrix multiply for 3d geo2mag.
         !
         call ESMF_FieldRegrid(srcfield, dstfield, routehandle_geo2mag,          &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid geo2mag 3D', rc)
      end if
   end subroutine edyn_esmf_regrid_geo2mag
   !-----------------------------------------------------------------------

   subroutine edyn_esmf_regrid_mag2geo(srcfield, dstfield, ndim)
     !
     ! Args:
     integer                         :: ndim
     type(ESMF_Field), intent(inout) :: srcfield, dstfield
     !
     ! Local:
     integer                     :: rc
     character(len=*), parameter :: subname = 'edyn_esmf_regrid_mag2geo'
     !
     if (ndim == 2) then
#if 0
        !
        ! Do sparse matrix multiply for 2d geo2mag.
        !
        call ESMF_FieldRegrid(srcfield, dstfield, routehandle_mag2geo_2d,       &
             termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
        call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid geo2mag 2D', rc)
#else
        call endrun(subname//' need routehandle_mag2geo_2d ')
#endif
     else ! 3d geo2mag
        !
        ! Do sparse matrix multiply for 3d geo2mag.
        !
        call ESMF_FieldRegrid(srcfield, dstfield, routehandle_mag2geo,          &
             termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
        call edyn_esmf_chkerr(subname, 'ESMF_FieldRegrid geo2mag 3D', rc)
     end if
   end subroutine edyn_esmf_regrid_mag2geo

   !-----------------------------------------------------------------------

   subroutine edyn_esmf_update_flag( flag )
     logical, intent(in) :: flag
     edyn_esmf_update_step=flag
   end subroutine edyn_esmf_update_flag

   !-----------------------------------------------------------------------

   subroutine edyn_esmf_update_phys_mesh(new_phys_mesh)
      use ESMF, only: ESMF_Mesh, ESMF_MeshIsCreated, ESMF_MeshDestroy

      ! Dummy argument
      type(ESMF_Mesh), intent(in) :: new_phys_mesh

      ! Ignore return code here as all we need is an attempt to reclaim memory
      if (ESMF_MeshIsCreated(phys_mesh)) then
         call ESMF_MeshDestroy(phys_mesh)
      end if

      phys_mesh = new_phys_mesh

   end subroutine edyn_esmf_update_phys_mesh

#endif
end module edyn_esmf
