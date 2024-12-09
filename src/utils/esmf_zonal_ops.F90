module esmf_zonal_ops
  use shr_kind_mod, only: r8 => shr_kind_r8, cl=>SHR_KIND_CL
  use ppgrid, only: pcols, pver, begchunk, endchunk
  use phys_grid, only: get_ncols_p
  use cam_logfile,  only: iulog
  use cam_abortutils, only: endrun

  use spmd_utils, only: masterproc, mpicom
  use esmf_phys_mesh_mod, only: esmf_phys_mesh_init, physics_grid_mesh
  use shr_reprosum_mod,only: shr_reprosum_calc

  use ESMF, only: ESMF_Grid, ESMF_Field
  use ESMF, only: ESMF_SUCCESS, ESMF_KIND_R8, ESMF_TYPEKIND_R8
  use ESMF, only: ESMF_GridCreate1PeriDim, ESMF_INDEX_GLOBAL
  use ESMF, only: ESMF_GridAddCoord, ESMF_GridGetCoord, ESMF_STAGGERLOC_CENTER
  use ESMF, only: ESMF_ArraySpec, ESMF_MESHLOC_ELEMENT, ESMF_ArraySpecSet, ESMF_FieldCreate
  use ESMF, only: ESMF_FieldRegridStore, ESMF_REGRIDMETHOD_BILINEAR, ESMF_POLEMETHOD_ALLAVG
  use ESMF, only: ESMF_EXTRAPMETHOD_NEAREST_IDAVG, ESMF_RouteHandle
  use ESMF, only: ESMF_KIND_I4
  use ESMF, only: ESMF_FieldGet, ESMF_FieldRegrid, ESMF_TERMORDER_SRCSEQ

  implicit none

  integer :: nlats = -1
  integer :: nlons = -1

  real(r8), allocatable :: glats(:)
  real(r8), allocatable :: glons(:)

  integer, parameter :: minlats_per_pe = 2
  integer, parameter :: minlons_per_pe = 2
  integer :: ntasks_lat = -1
  integer :: ntasks_lon = -1
  integer :: npes = -1
  integer :: mytid = -1
  integer :: mytidi = -1
  integer :: mytidj = -1

  integer :: lon_beg = -1
  integer :: lon_end = -1
  integer :: lat_beg = -1
  integer :: lat_end = -1

  integer, allocatable :: itask_table(:,:)

  type(ESMF_RouteHandle) :: rh_phys2lonlat_2D
  type(ESMF_RouteHandle) :: rh_phys2lonlat_3D
  type(ESMF_Field) :: physfld_3d, physfld_2d
  type(ESMF_Field) :: lonlatfld_3d, lonlatfld_2d

  integer :: mynlats, mynlons

  integer :: rows_comm  ! communicators for each task row

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_zonal_ops_init(nlats_in)
    use phys_grid, only: get_grid_dims
    use mpi, only: mpi_comm_size, mpi_comm_rank, MPI_PROC_NULL, MPI_INTEGER

    integer, optional, intent(in) :: nlats_in

    integer :: hdim1_d,hdim2_d, nlat0, i,j, ndx, ierr, irank

    integer, parameter  :: ndelts = 8
    real(r8), parameter :: deltas(ndelts) = (/ 0.125_r8, 0.25_r8, 0.5_r8, 1.0_r8, 2.0_r8, 5.0_r8, 6.0_r8, 10._r8 /)
    real(r8) :: delt0, diff(ndelts), delx, dely

    integer :: n
    integer :: lons_per_task, lons_overflow, lats_per_task, lats_overflow
    integer :: task_cnt
    character(len=*), parameter :: subname  = 'esmf_zonal_ops_init'

    integer, allocatable :: petmap(:,:,:)
    integer :: petcnt

    integer, allocatable :: mytidi_send(:)
    integer, allocatable :: mytidj_send(:)
    integer, allocatable :: mytidi_recv(:)
    integer, allocatable :: mytidj_recv(:)

    integer, allocatable :: nlons_send(:)
    integer, allocatable :: nlats_send(:)
    integer, allocatable :: nlons_recv(:)
    integer, allocatable :: nlats_recv(:)

    integer, allocatable :: nlons_task(:)
    integer, allocatable :: nlats_task(:)

    type(ESMF_Grid) :: lonlat_grid

    integer                       :: lbnd_lat, ubnd_lat, lbnd_lon, ubnd_lon
    integer                       :: lbnd(1), ubnd(1)
    real(ESMF_KIND_R8), pointer   :: coordX(:), coordY(:)

    type(ESMF_ArraySpec) :: arrayspec
    integer(ESMF_KIND_I4), pointer :: factorIndexList(:,:)
    real(ESMF_KIND_R8),    pointer :: factorList(:)
    integer                        :: smm_srctermproc,  smm_pipelinedep

    ! create phys grid mesh
    call esmf_phys_mesh_init()

    ! create reg lat/lon grid for zonal mean calcs

    if (present(nlats_in)) then

       ! user specified resolution
       nlats = nlats_in
       nlons = 2*nlats

    else
       call get_grid_dims(hdim1_d,hdim2_d)

       if (hdim2_d>1) then

          ! on reg lat / lon FV grid
          nlons = hdim1_d
          nlats = hdim2_d
          delx = 360._r8/nlons
          dely = 180._r8/(nlats-1)

       else

          ! create grid close to the unstructured grid resolution
          nlat0 = sqrt(0.5_r8*hdim1_d) + 1
          delt0 = 180._r8/nlat0

          do i = 1, ndelts
             diff(i) = abs(deltas(i)-delt0)
          end do

          ndx = minloc(diff,1)

          nlats = 180._r8/deltas(ndx)

          nlons = 2*nlats
          nlats = nlats+1

          delx = deltas(ndx)
          dely = delx
       end if
    end if

    allocate(glons(nlons))
    allocate(glats(nlats))

    glons(1) = 0._r8
    glats(1) = -90._r8

    do i = 2,nlons
       glons(i) = glons(i-1) + delx
    end do
    do i = 2,nlats
       glats(i) = glats(i-1) + dely
    end do

    ! decompose the grid across mpi tasks ...

    call mpi_comm_size(mpicom, npes, ierr)
    call mpi_comm_rank(mpicom, mytid, ierr)

    do ntasks_lon = 1,nlons
       ntasks_lat = npes/ntasks_lon
       if ( (minlats_per_pe*ntasks_lat<nlats) .and. (ntasks_lat*ntasks_lon==npes) ) then
          exit
       endif
    end do

    if (masterproc) then
       write(iulog,'(a,3i6)') subname//': npes,nlons,nlats: ',npes,nlons,nlats
       write(iulog,'(a,2i6)') subname//': ntasks_lon,ntasks_lat : ',ntasks_lon,ntasks_lat
    endif

    if (ntasks_lat*ntasks_lon/=npes) then
       call endrun(subname//': ntasks_lat*ntasks_lon/=npes')
    endif


    ! figure out the starting and ending coordinates
    lons_per_task = nlons / ntasks_lon
    lons_overflow = MOD(nlons, ntasks_lon)
    lats_per_task = nlats / ntasks_lat
    lats_overflow = MOD(nlats, ntasks_lat)
    lon_beg = 1
    lon_end = 0
    lat_beg = 1
    lat_end = 0
    task_cnt= 0
    if (mytid<npes) then
       jloop: do j = 0,ntasks_lat-1
          lat_beg = lat_end + 1
          lat_end = lat_beg + lats_per_task - 1
          if (j<lats_overflow) then
             lat_end = lat_end + 1
          end if
          lon_end = 0
          do i = 0,ntasks_lon-1
             lon_beg = lon_end + 1
             lon_end = lon_beg + lons_per_task - 1
             if (i<lons_overflow) then
                lon_end = lon_end + 1
             end if
             task_cnt = task_cnt+1
             if (task_cnt>mytid) exit jloop
          end do
       enddo jloop
    endif

    mynlats = lat_end-lat_beg+1
    mynlons = lon_end-lon_beg+1

    if (mynlats<minlats_per_pe) then
       call endrun(subname//': mynlats < minlats_per_pe')
    end if
    if (mynlons<minlons_per_pe) then
       call endrun(subname//': mynlons < minlons_per_pe')
    end if

    !
    ! Allocate and set 2d table of tasks:
    !
    allocate(itask_table(-1:ntasks_lon,0:ntasks_lat-1),stat=ierr)
    if (ierr /= 0) then
       write(iulog,"(a,': Error allocating itask_table: ntaski,j=',2i4)") subname,ntasks_lon,ntasks_lat
       call endrun(subname//': Error allocating itask_table')
    endif

    ! setup MPI task table

    itask_table(:,:) = MPI_PROC_NULL

    irank = 0
    mytidi = -1
    mytidj = -1
    do j =  0, ntasks_lat-1
       do i = 0, ntasks_lon-1
          itask_table(i,j) = irank
          if (mytid == irank) then
             mytidi = i
             mytidj = j
          end if
          irank = irank+1
       end do
       !
       ! Tasks are periodic in longitude:
       !
       itask_table(-1,j) = itask_table(ntasks_lon-1,j)
       itask_table(ntasks_lon,j) = itask_table(0,j)

    end do ! j=0,ntaskj-1

    call mpi_comm_split(mpicom,mytidj,mytid,rows_comm,ierr)

    allocate(mytidi_send(npes))
    allocate(mytidj_send(npes))
    allocate(mytidi_recv(npes))
    allocate(mytidj_recv(npes))

    mytidi_send = mytidi
    mytidj_send = mytidj
    call mpi_alltoall(mytidi_send, 1, MPI_INTEGER, mytidi_recv, 1, MPI_INTEGER, mpicom, ierr)
    call mpi_alltoall(mytidj_send, 1, MPI_INTEGER, mytidj_recv, 1, MPI_INTEGER, mpicom, ierr)


! set up zm grid ...


    allocate(petmap(ntasks_lon,ntasks_lat,1))

    petcnt = 0
    do j = 1,ntasks_lat
       do i = 1,ntasks_lon
          petmap(i,j,1) = petcnt
          petcnt = petcnt+1
       end do
    end do

    allocate(nlons_send(npes))
    allocate(nlats_send(npes))
    allocate(nlons_recv(npes))
    allocate(nlats_recv(npes))

    nlons_send(:) = mynlons
    nlats_send(:) = mynlats

    call mpi_alltoall(nlons_send, 1, MPI_INTEGER, nlons_recv, 1, MPI_INTEGER, mpicom, ierr)
    call mpi_alltoall(nlats_send, 1, MPI_INTEGER, nlats_recv, 1, MPI_INTEGER, mpicom, ierr)

    deallocate(nlons_send)
    deallocate(nlats_send)


    allocate(nlons_task(ntasks_lon))
    allocate(nlats_task(ntasks_lat))


    do i = 1, ntasks_lon
       loop: do n = 1, npes
          if (mytidi_recv(n) == i-1) then
             nlons_task(i) = nlons_recv(n)
             exit loop
          end if
       end do loop
    end do

    do j = 1, ntasks_lat
       loop1: do n = 1, npes
          if (mytidj_recv(n) == j-1) then
             nlats_task(j) = nlats_recv(n)
             exit loop1
          end if
       end do loop1
    end do


    ! Create 2d geographic source grid (with poles)
    lonlat_grid = ESMF_GridCreate1PeriDim(                       &
         countsPerDEDim1=nlons_task, coordDep1=(/1/),         &
         countsPerDEDim2=nlats_task, coordDep2=(/2/), petmap=petmap, &
         indexflag=ESMF_INDEX_GLOBAL,minIndex=(/1,1/), rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_GridCreate1PeriDim ERROR')


    ! Set coordinates:

    call ESMF_GridAddCoord(lonlat_grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_GridAddCoord ERROR')

    if (mytid<npes) then
       call ESMF_GridGetCoord(lonlat_grid, coordDim=1, &
            computationalLBound=lbnd, computationalUBound=ubnd,  &
            farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=ierr)
       call check_esmf_error(ierr, subname//'ESMF_GridGetCoord for longitude coords ERROR')

       lbnd_lon = lbnd(1)
       ubnd_lon = ubnd(1)
       do i = lbnd_lon, ubnd_lon
          coordX(i) = glons(i)
       end do

       call ESMF_GridGetCoord(lonlat_grid, coordDim=2, &
            computationalLBound=lbnd, computationalUBound=ubnd, &
            farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=ierr)
       call check_esmf_error(ierr, subname//'ESMF_GridGetCoord for latitude coords ERROR')

       lbnd_lat = lbnd(1)
       ubnd_lat = ubnd(1)
       do i = lbnd_lat, ubnd_lat
          coordY(i) = glats(i)
       end do
    end if

    ! create ESMF fields

    ! 3D phys fld
    call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_ArraySpecSet 3D phys fld ERROR')

    physfld_3d = ESMF_FieldCreate(physics_grid_mesh, arrayspec, &
                                  gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, &
                                  ungriddedLBound=(/1/), ungriddedUBound=(/pver/), rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_FieldCreate 3D phys fld ERROR')

    ! 2D phys fld
    call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_ArraySpecSet 2D phys fld ERROR')

    physfld_2d = ESMF_FieldCreate(physics_grid_mesh, arrayspec, &
                                  meshloc=ESMF_MESHLOC_ELEMENT, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_FieldCreate 2D phys fld ERROR')

    ! 3D lon/lat grid
    call ESMF_ArraySpecSet(arrayspec, 3, ESMF_TYPEKIND_R8, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_ArraySpecSet 3D lonlat fld ERROR')

    lonlatfld_3d = ESMF_FieldCreate( lonlat_grid, arrayspec, staggerloc=ESMF_STAGGERLOC_CENTER, &
                                     ungriddedLBound=(/1/), ungriddedUBound=(/pver/), rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_FieldCreate 3D lonlat fld ERROR')

    ! 2D lon/lat grid
    call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_ArraySpecSet 2D lonlat fld ERROR')

    lonlatfld_2d = ESMF_FieldCreate( lonlat_grid, arrayspec, staggerloc=ESMF_STAGGERLOC_CENTER, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_FieldCreate 2D lonlat fld ERROR')


    ! route handles -- phys --> lonlat mapping

    smm_srctermproc = 0
    smm_pipelinedep = 16

    ! 3D

    call ESMF_FieldRegridStore(srcField=physfld_3d, dstField=lonlatfld_3d, &
         regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
         polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
         extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
         routeHandle=rh_phys2lonlat_3D, factorIndexList=factorIndexList, &
         factorList=factorList, srcTermProcessing=smm_srctermproc,          &
         pipelineDepth=smm_pipelinedep, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_FieldRegridStore 3D routehandle ERROR')

    ! 2D

    call ESMF_FieldRegridStore(srcField=physfld_2d, dstField=lonlatfld_2d, &
         regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
         polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
         extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
         routeHandle=rh_phys2lonlat_2D, factorIndexList=factorIndexList, &
         factorList=factorList, srcTermProcessing=smm_srctermproc,          &
         pipelineDepth=smm_pipelinedep, rc=ierr)
    call check_esmf_error(ierr, subname//'ESMF_FieldRegridStore 2D routehandle ERROR')

  end subroutine esmf_zonal_ops_init


  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function esmf_zonal_mean_2d(physfld) result(zmfld)

    real(r8),intent(in) :: physfld(pcols,begchunk:endchunk)

    real(r8) :: zmfld(lat_beg:lat_end)

    integer :: rc, i, ichnk, icol, ilon, ilat, ncol
    real(ESMF_KIND_R8), pointer :: physptr(:)
    real(ESMF_KIND_R8), pointer :: lonlatptr(:,:)

    character(len=*), parameter :: subname  = 'esmf_zm_calc_2d: '

    real(r8) :: arr(lon_beg:lon_end,1)
    real(r8) :: gsum(1)

    ! regrid to lat/lon

    call ESMF_FieldGet(physfld_2d, localDe=0, farrayPtr=physptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet physptr')

    i = 0
    do ichnk = begchunk, endchunk
       ncol = get_ncols_p(ichnk)
       do icol = 1,ncol
          i = i+1
          physptr(i) = physfld(icol,ichnk)
       end do
    end do

    call ESMF_FieldRegrid(physfld_2d, lonlatfld_2d, rh_phys2lonlat_2D, &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegrid physfld_2d->lonlatfld_2d')

    call ESMF_FieldGet(lonlatfld_2d, localDe=0, farrayPtr=lonlatptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet lonlatptr')

    ! zonal mean

    do ilat = lat_beg, lat_end
       arr(lon_beg:lon_end,1) = lonlatptr(lon_beg:lon_end,ilat)
       call shr_reprosum_calc(arr, gsum, mynlons, mynlons, 1, gbl_count=nlons, commid=rows_comm)
       zmfld(ilat) = gsum(1)/nlons
    end do

  end function esmf_zonal_mean_2d

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function esmf_zonal_mean_3d(physfld) result(zmfld)

    real(r8),intent(in) :: physfld(pver,pcols,begchunk:endchunk)

    real(r8) :: zmfld(lat_beg:lat_end,pver)


    integer :: rc, i, ichnk, icol, ilon, ilat, ilev, ncol
    real(ESMF_KIND_R8), pointer :: physptr(:,:)
    real(ESMF_KIND_R8), pointer :: lonlatptr(:,:,:)

    real(r8) :: arr(lon_beg:lon_end,pver)
    real(r8) :: gsum(pver)

    character(len=*), parameter :: subname  = 'esmf_zm_calc_3d: '

    integer :: k

    ! regrid to lat/lon

    call ESMF_FieldGet(physfld_3d, localDe=0, farrayPtr=physptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet physptr')

    i = 0
    do ichnk = begchunk, endchunk
       ncol = get_ncols_p(ichnk)
       do icol = 1,ncol
          i = i+1
          do ilev = 1,pver
             physptr(ilev,i) = physfld(ilev,icol,ichnk)
          end do
       end do
    end do

    call ESMF_FieldRegrid(physfld_3d, lonlatfld_3d, rh_phys2lonlat_3D, &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegrid physfld_3d->lonlatfld_3d')

    call ESMF_FieldGet(lonlatfld_3d, localDe=0, farrayPtr=lonlatptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet lonlatptr')

    ! zonal mean

    do ilat = lat_beg, lat_end
       arr(lon_beg:lon_end,:) = lonlatptr(lon_beg:lon_end,ilat,:)
       call shr_reprosum_calc(arr, gsum, mynlons, mynlons, pver, gbl_count=nlons, commid=rows_comm)
       zmfld(ilat,:) = gsum(:)/nlons
    end do

  end function esmf_zonal_mean_3d

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine check_esmf_error( rc, errmsg )
    integer, intent(in) :: rc
    character(len=*), intent(in) :: errmsg

    character(len=cl) :: errstr

    if (rc /= ESMF_SUCCESS) then
       write(errstr,'(a,i6)') 'esmf_zonal_ops::'//trim(errmsg)//' -- ESMF ERROR code: ',rc
       if (masterproc) write(iulog,*) trim(errstr)
       call endrun(trim(errstr))
    end if

  end subroutine check_esmf_error



end module esmf_zonal_ops
