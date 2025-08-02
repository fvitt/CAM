module edyn3d_esmf_oplus_grid_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc
  use edyn3d_esmf_error_mod, only: check_error

  use ESMF, only: ESMF_COORDSYS_SPH_DEG, ESMF_INDEX_GLOBAL, ESMF_KIND_R8, ESMF_Grid,  ESMF_GridCreate1PeriDim
  use ESMF, only: ESMF_GridAddCoord, ESMF_GridGetCoord, ESMF_STAGGERLOC_CENTER
  use ESMF, only: ESMF_GridDestroy

  private

  public :: edyn3d_esmf_oplus_grid_init
  public :: edyn3d_esmf_oplus_grid_destroy
  public :: oplus_grid

  type(ESMF_Grid), protected :: oplus_grid

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_esmf_oplus_grid_init()

    use edyn_geogrid, only: opluslon=>glon, opluslat=>glat ! oplus grid coordinates
    use edyn_mpi, only: edyn_ntask=>ntask, edyn_ntaski=>ntaski, edyn_ntaskj=>ntaskj, edyn_tasks=>tasks
    use edyn_mpi, only: mytid

    integer, allocatable :: petmap(:,:,:)
    integer :: petcnt
    integer :: nlons_task(edyn_ntaski) ! # number of lons per task
    integer :: nlats_task(edyn_ntaskj) ! # number of lats per task
    integer :: i,j, n, rc, astat
    real(ESMF_KIND_R8), pointer   :: coordX(:), coordY(:)
    integer :: lbnd(1), ubnd(1)

    character(len=*), parameter :: subname = 'dyn3d_esmf_oplus_grid_init'

    ! Oplus grid / field

    nlons_task = 0
    nlats_task = 0

    allocate(petmap(edyn_ntaski,edyn_ntaskj,1), stat=astat)
    if (astat/=0) then
       call endrun(subname//' : not able to allocate petmap array')
    end if

    petcnt = 0
    do j = 1,edyn_ntaskj
       do i = 1,edyn_ntaski
          petmap(i,j,1) = petcnt
          petcnt = petcnt+1
       end do
    end do

    do i = 1,edyn_ntaski
       loop: do n = 0, edyn_ntask-1
          if (edyn_tasks(n)%mytidi == i-1) then
             nlons_task(i) = edyn_tasks(n)%nlons
             exit loop
          end if
       end do loop
    end do
    !
    do j = 1, edyn_ntaskj
       loop1: do n = 0, edyn_ntask-1
          if (edyn_tasks(n)%mytidj == j-1) then
             nlats_task(j) = edyn_tasks(n)%nlats
             exit loop1
          end if
       end do loop1
    end do

    oplus_grid = ESMF_GridCreate1PeriDim(               &
           countsPerDEDim1=nlons_task, coordDep1=(/1/), &
           countsPerDEDim2=nlats_task, coordDep2=(/2/), &
           petmap=petmap, &
           coordSys=ESMF_COORDSYS_SPH_DEG, &
           indexflag=ESMF_INDEX_GLOBAL,minIndex=(/1,1/), rc=rc)
    call check_error(subname,'ESMF_GridCreate1PeriDim oplus_grid',rc)

    deallocate(petmap)

    ! set up coordinates:

    call ESMF_GridAddCoord(oplus_grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call check_error(subname,'ESMF_GridAddCoord  oplus_grid',rc)

    if (mytid<edyn_ntask) then

       ! Lon Coord
       call ESMF_GridGetCoord(oplus_grid, coordDim=1, localDE=0, &
            computationalLBound=lbnd, computationalUBound=ubnd, &
            farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
       call check_error(subname,'ESMF_GridGetCoord longitude coord oplus_grid',rc)

       do i = lbnd(1),ubnd(1)
          coordX(i) = opluslon(i)
       end do

       ! Lat Coord
       call ESMF_GridGetCoord(oplus_grid, coordDim=2, localDE=0, &
            computationalLBound=lbnd, computationalUBound=ubnd, &
            farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
       call check_error(subname,'ESMF_GridGetCoord latitude coord oplus_grid',rc)

       do j = lbnd(1),ubnd(1)
          coordY(j) = opluslat(j)
       end do

    end if

  end subroutine edyn3d_esmf_oplus_grid_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_esmf_oplus_grid_destroy()

    integer :: rc
    character(len=*), parameter :: subname = 'dyn3d_esmf_oplus_grid_destroy'

    call ESMF_GridDestroy(oplus_grid, rc=rc)
    call check_error(subname,'ESMF_GridDestroy oplus_grid', rc)

  end subroutine edyn3d_esmf_oplus_grid_destroy

end module edyn3d_esmf_oplus_grid_mod
