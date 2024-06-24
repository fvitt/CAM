module edyn3D_oplus_grid
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  use ESMF

  implicit none

  type(ESMF_Grid) :: oplsGrid2D
  type(ESMF_Field) :: oplsField2D

  integer, allocatable :: petmap(:,:,:)

contains

  subroutine edyn3D_oplus_grid_init

    use edyn_mpi, only: mytid, ntask, ntaski, ntaskj, tasks

    use edyn_geogrid,   only: glon, glat

    integer :: i,j,n, localrc
    integer :: petcnt

    integer :: nlons_task(ntaski) ! # number of lons per task
    integer :: nlats_task(ntaskj) ! # number of lats per task

    integer :: lbnd_lat, ubnd_lat, lbnd_lon, ubnd_lon
    integer :: lbnd(1), ubnd(1)

    real(ESMF_KIND_R8), pointer   :: coordX(:), coordY(:)

    character(len=*), parameter :: subname = 'edyn3D_oplus_grid_init'

    real(r8) :: alt0, delz
    integer :: k

    type(ESMF_ArraySpec) :: arrayspec2D

    nlons_task = 0
    nlats_task = 0

    allocate(petmap(ntaski,ntaskj,1))

    petcnt = 0
    do j = 1,ntaskj
       do i = 1,ntaski
          petmap(i,j,1) = petcnt
          petcnt = petcnt+1
       end do
    end do

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

    oplsGrid2D = ESMF_GridCreate1PeriDim(               &
           countsPerDEDim1=nlons_task, coordDep1=(/1/), &
           countsPerDEDim2=nlats_task, coordDep2=(/2/), &
           petmap=petmap, &
           coordSys=ESMF_COORDSYS_SPH_DEG, &
           indexflag=ESMF_INDEX_GLOBAL,minIndex=(/1,1/), rc=localrc)

    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_FieldCreate oplsGrid2D ')
    end if

    ! set up coordinates:

    call ESMF_GridAddCoord(oplsGrid2D, staggerloc=ESMF_STAGGERLOC_CENTER, rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_GridAddCoord oplsGrid2D ')
    end if

    if (mytid<ntask) then

       ! Lon Coord
       call ESMF_GridGetCoord(oplsGrid2D, coordDim=1, localDE=0, &
            computationalLBound=lbnd, computationalUBound=ubnd, &
            farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=localrc)
       if (ESMF_LogFoundError(localrc)) then
          call endrun(subname//' ESMF_GridGetCoord longitude coord oplsGrid2D')
       end if

       lbnd_lon = lbnd(1)
       ubnd_lon = ubnd(1)
       do i = lbnd_lon, ubnd_lon
          coordX(i) = glon(i)
       end do

       ! Lat Coord
       call ESMF_GridGetCoord(oplsGrid2D, coordDim=2, localDE=0, &
            computationalLBound=lbnd, computationalUBound=ubnd, &
            farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=localrc)
       if (ESMF_LogFoundError(localrc)) then
          call endrun(subname//' ESMF_GridGetCoord longitude coord oplsGrid2D')
       end if

       lbnd_lat = lbnd(1)
       ubnd_lat = ubnd(1)
       do i = lbnd_lat, ubnd_lat
          coordY(i) = glat(i)
       end do

    end if

    ! Create 2D field
    call ESMF_ArraySpecSet(arrayspec2d, 2, ESMF_TYPEKIND_R8, rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_ArraySpecSet for oplsField2D')
    end if

    oplsField2D = ESMF_FieldCreate(oplsGrid2D, arrayspec2d, &
         staggerloc=ESMF_STAGGERLOC_CENTER, name="oplsField2D", rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_FieldCreate oplsField2D')
    end if

  end subroutine edyn3D_oplus_grid_init

end module edyn3D_oplus_grid
