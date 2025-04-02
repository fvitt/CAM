module edyn3d_esmf_s1_mag_grid_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc, mpicom

  use ESMF

  implicit none

  type(ESMF_Grid), allocatable :: mag_s1_fdln_grid(:)

contains


  !%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>
  !       S1
  !     o--+--o P
  !  ^  |     |
  !  L  *     * S2
  !  A  |     |
  !  T  o--+--o
  !     LON -->
  !%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>%>

  subroutine edyn3d_esmf_s1_mag_grid_init
    use fieldline_module, only: npts_s1, glon_s1, glat_s1
    use params_module, only: nz=>nhgt_fix, nmlon, nmlat_h
    use mpi_module, only: mytid=>mpi_rank, ntask=>mpi_size, lon_size, lat_size
    use mpi_module, only: mlon0,mlon1,mlat0,mlat1, mlon0_task, mlon1_task, mlat0_task, mlat1_task
    use mpi_module, only: nmlon_task, nmlat_task

    integer :: rc
    !logical :: found_eq
    integer :: i,j,k,n, ii,jj, isn, nde
    integer, allocatable :: petmap(:,:,:)

    integer :: petcnt
    integer :: ncells_hlat
    integer :: klat_sz, j0, j1

    integer :: localDECount
    integer :: lonCellsPerDE(lon_size) ! number of lons per task
    integer, allocatable :: latCellsPerDE(:)  ! number of lons per task in each level

    integer :: lbnd(2),ubnd(2)
    real(kind=ESMF_KIND_R8), pointer :: loncoord(:,:), latcoord(:,:)

    character(len=*), parameter :: subname = 'edyn3d_esmf_s1_mag_grid_init'
    real(r8), parameter :: NOTSET = -huge(1._r8)

    allocate(mag_s1_fdln_grid(nz))

    do i = 1,lon_size
       lonCellsPerDE(i) = nmlon_task(i-1)
    end do

    print*,' '
    print*,'FVDBG.edyn3d_esmf_s1_mag_grid_init.. '

    vertloop: do k = 1, nz

       ! total number of grids cells per hemisphere
       !ncells_hlat = count(npts_s1>=k)
       ncells_hlat = nmlat_h - (k-1)

       print*,' ncells_hlat:',ncells_hlat  !,'   nmlat_h-(k-1):',nmlat_h - (k-1)

       ! find number of lat tasks that have grid cells in level k in 1 hemisphere
       i = 0
       n = nmlat_task(i)
       do while (n < ncells_hlat)
          i = i + 1
          n = n + nmlat_task(i)
       end do
       klat_sz = i+1

       print*,'FVDBG... k:',k,' number of lat tasks that have grid cells klat_sz:',klat_sz

       ! number of global lat grid cells for level k (1 DE straddles the equator)
       allocate(latCellsPerDE(klat_sz*2-1))
       latCellsPerDE = -huge(1)

       ! south pole to north pole
       ! southern hemisphere first
       do i = 1,klat_sz-1
          ii = (i-1)*lon_size + 1
          j0 = mlat0_task(ii-1)
          j1 = mlat1_task(ii-1)
          latCellsPerDE(i) = count( npts_s1(j0:j1)>=k )
       end do

       ! adjacent to the equator
       i = klat_sz
       ii = (i-1)*lon_size + 1
       j0 = mlat0_task(ii-1)
       j1 = min(mlat1_task(ii-1),ncells_hlat)
       latCellsPerDE(i) = 2*count( npts_s1(j0:j1)>=k ) -1 ! each side of the equator (only 1 at the equator)

       ! norther hemisphere
       do i = klat_sz+1,klat_sz*2-1
          ii = klat_sz*2-1 - i + 1
          latCellsPerDE(i) = latCellsPerDE(ii) ! mirror the southern hemisphere
       end do

       print*,'FVDBG... k:',k,' latCellsPerDE:',latCellsPerDE

       ! mpi task number for each DE (numLonDEs x numLatDEs)
       allocate(petmap(lon_size, klat_sz*2-1,1))
       petmap = -huge(1)

       petcnt = 0
       do j = 1,klat_sz
          do i = 1,lon_size
             petmap(i,j,1) = petcnt
             petcnt = petcnt+1
          end do
       end do
       petcnt = 0
       do j = 2*klat_sz-1,klat_sz,-1
          do i = 1,lon_size
             petmap(i,j,1) = petcnt
             petcnt = petcnt+1
          end do
       end do

       print*,'  petmap : '
       do j = 1, klat_sz*2-1
          print*,petmap(1:lon_size,j,1)
       end do

       ! 1 periodic dimension -- periodic logitude dim
       mag_s1_fdln_grid(k) = ESMF_GridCreate1PeriDim(  &
           countsPerDEDim1=lonCellsPerDE, coordDep1=(/1,2/), &
           countsPerDEDim2=latCellsPerDE, coordDep2=(/1,2/), petmap=petmap, &
           indexflag=ESMF_INDEX_GLOBAL,rc=rc)
       call check_error(subname, 'ESMF_GridCreate1PeriDim', rc)

       deallocate(petmap)
       deallocate(latCellsPerDE)

       ! get number of DEs for this MPI task
       call ESMF_GridGet(mag_s1_fdln_grid(k), localDECount=localDECount, rc=rc)
       call check_error(subname,'ESMF_GridGet localDECount',rc)

       ! S1 coordinates
       call ESMF_GridAddCoord(grid=mag_s1_fdln_grid(k),staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
       call check_error(subname,'ESMF_GridAddCoord mag_s1_fdln_grid',rc)

       do nde = 0,localDECount-1

          print*,'DE number nde: ',nde

          ! set S1 coordinates
          ! geographic longitudes
          call ESMF_GridGetCoord(mag_s1_fdln_grid(k), coordDim=1, localDE=nde, &
               computationalLBound=lbnd, computationalUBound=ubnd, &
               staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=loncoord, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord  S1',rc)
          loncoord = NOTSET

          ! geographic latitudes
          call ESMF_GridGetCoord(mag_s1_fdln_grid(k), coordDim=2, localDE=nde, &
               computationalLBound=lbnd, computationalUBound=ubnd, &
               staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=latcoord, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord  S1',rc)
          latcoord = NOTSET

          print*,' S1 grid lat lbnd : ',lbnd
          print*,' S1 grid lat ubnd : ',ubnd

          do i = lbnd(1),ubnd(1)
             do j = lbnd(2),ubnd(2)
                if (j>ncells_hlat) then
                   isn = 2
                   jj = 2*(ncells_hlat-1)+1 - j + 1
                else
                   isn = 1
                   jj = j
                end if
                loncoord(i,j) = glon_s1(k,isn,jj,i)
                latcoord(i,j) = glat_s1(k,isn,jj,i)
             end do
          end do

          if (any(loncoord==NOTSET)) then
             call endrun(subname//' loncoord not set correctly')
          end if
          if (any(latcoord==NOTSET)) then
             call endrun(subname//' latcoord not set correctly')
          end if

       end do

    end do vertloop

    print*,' '

!!$print*,'FVDBG.edyn3d_esmf_s1_mag_grid_init.. OK HERE'
!!$call mpi_barrier(mpicom, rc)
!!$call endrun('FVDBG.edyn3d_esmf_s1_mag_grid_init.. OK STOP HERE')

  end subroutine edyn3d_esmf_s1_mag_grid_init
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine check_error(subname, routine, rc)

    character(len=*), intent(in) :: subname
    character(len=*), intent(in) :: routine
    integer,          intent(in) :: rc

    character(len=cl) :: errmsg

    if (rc /= ESMF_SUCCESS) then
       write(errmsg, '(4a,i0)') trim(subname), ': Error return from ', trim(routine), ', rc = ', rc
       if (masterproc) then
          write(iulog, '(2a)') 'ERROR: ', trim(errmsg)
       end if
       call endrun(trim(errmsg))
    end if
  end subroutine check_error


end module edyn3d_esmf_s1_mag_grid_mod
