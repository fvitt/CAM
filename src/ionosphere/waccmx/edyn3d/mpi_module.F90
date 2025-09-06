module mpi_module

  use prec, only: rp
  use MPI
  use iso_fortran_env, only: real32,real64

  implicit none

  integer, protected :: dynamo_world=-huge(1), &
    mpi_rp=-huge(1), mpi_size=0, mpi_rank=-1, &
    lat_size=0, lon_size=0, lat_rank=-1, lon_rank=-1, &
    nmlat=0, maxmlat=-1, mlat0=1, mlat1=0, mlatd0=1, mlatd1=0, &
    nmlon=0, maxmlon=-1, mlon0=1, mlon1=0, mlond0=1, mlond1=0

  integer, protected, dimension(:), allocatable :: &
    nmlat_task, mlat0_task, mlat1_task, &
    nmlon_task, mlon0_task, mlon1_task

  interface gather_mag ! gather magnetic fields
    module procedure gather_mag_2d, gather_mag_3d, gather_mag_4d, gather_mag_5d
  endinterface

  interface bcast ! broadcast fields
    module procedure bcast_2d, bcast_3d
  endinterface bcast


  contains
!-----------------------------------------------------------------------
  subroutine init(mpi_comm_host, npes_edyn3d)

    integer, intent(in) :: mpi_comm_host
    integer, intent(in) :: npes_edyn3d

    integer :: ierror
    integer :: color, npes_host

    if (rp == real32) then
      mpi_rp = MPI_REAL4
    elseif (rp == real64) then
      mpi_rp = MPI_REAL8
    else
      stop 'unknown real precision'
    endif

    call mpi_comm_size(mpi_comm_host, npes_host, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_size', ierror)

    mpi_size = min(npes_host,npes_edyn3d)

    call MPI_Comm_rank(mpi_comm_host, mpi_rank, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Comm_rank', ierror)

    color = mpi_rank/mpi_size
    call mpi_comm_split(mpi_comm_host, color, mpi_rank, dynamo_world, ierror)


! factorize MPI process number to the nearest two numbers
    do lat_size = int(sqrt(real(mpi_size, kind=rp))), 1, -1
      lon_size = mpi_size / lat_size
      if (lon_size*lat_size == mpi_size) exit ! lon_size >= lat_size
    enddo

! 2D index increases faster in latitudes
! (stack along latitudes first then longitudes)
! (lat_size=3)
!  8  9 10 11
!  4  5  6  7
!  0  1  2  3 (lon_size=4)
    lat_rank = mpi_rank / lon_size
    lon_rank = modulo(mpi_rank, lon_size)

  endsubroutine init
!-----------------------------------------------------------------------
  subroutine setup_topology(nmlat_in, nmlon_in)
! setup MPI decompositions in mag coordinates and the connectivity matrix

    integer, intent(in) :: nmlat_in, nmlon_in

    integer :: i, j, rnk, rnki, rnkj

    allocate(nmlat_task(0:lat_size-1))
    allocate(nmlon_task(0:lon_size-1))

    nmlat_task = 0
    nmlon_task = 0

    allocate(mlat0_task(0:mpi_size-1))
    allocate(mlat1_task(0:mpi_size-1))
    allocate(mlon0_task(0:mpi_size-1))
    allocate(mlon1_task(0:mpi_size-1))

    mlat0_task = 1
    mlon0_task = 1
    mlat1_task = -1
    mlon1_task = -1

    nmlat = nmlat_in
    nmlon = nmlon_in

! setup magnetic grid decomposition
! each process can have unequal number of latitudes or longitudes

    nmlat_task = generate_minvar_list(nmlat, lat_size)
    maxmlat = maxval(nmlat_task)
    if (lat_rank<lat_size) then
       mlat0 = 1
       do j = 0, lat_rank-1
          mlat0 = mlat0 + nmlat_task(j)
       enddo
       mlat1 = mlat0 + nmlat_task(lat_rank) - 1
    endif

    nmlon_task = generate_minvar_list(nmlon, lon_size)
    maxmlon = maxval(nmlon_task)
    if (lat_rank<lat_size) then
       mlon0 = 1
       do i = 0, lon_rank-1
          mlon0 = mlon0 + nmlon_task(i)
       enddo
       mlon1 = mlon0 + nmlon_task(lon_rank) - 1
    endif

! each process keeps a record of the lat-lon decomposition
    do concurrent (rnk = 0:mpi_size-1)
      rnkj = rnk / lon_size
      rnki = modulo(rnk, lon_size)

      mlat0_task(rnk) = 1
      do j = 0, rnkj-1
        mlat0_task(rnk) = mlat0_task(rnk) + nmlat_task(j)
      enddo
      mlat1_task(rnk) = mlat0_task(rnk) + nmlat_task(rnkj) - 1

      mlon0_task(rnk) = 1
      do i = 0, rnki-1
        mlon0_task(rnk) = mlon0_task(rnk) + nmlon_task(i)
      enddo
      mlon1_task(rnk) = mlon0_task(rnk) + nmlon_task(rnki) - 1
    enddo

! halos
    mlatd0 = mlat0 - 1
    mlatd1 = mlat1 + 1
    mlond0 = mlon0 - 1
    mlond1 = mlon1 + 1

  endsubroutine setup_topology
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  subroutine sync_mlat_5d(var, l, m, n)
! longitude halo points are not included

    use MPI

    integer, intent(in) :: l, m, n
    real(kind=rp), dimension(l, m, n, mlatd0:mlatd1, mlon0:mlon1), intent(inout) :: var

    integer :: below, above, cnt, i, lc, mc, nc, ierror
    integer, dimension(4) :: request
    real(kind=rp), dimension(l, m, n, maxmlon) :: &
      send_to_below, send_to_above, recv_from_below, recv_from_above

! find the rank of adjacent processes
    if (lat_rank == 0) then
      below = MPI_PROC_NULL
    else
      below = mpi_rank - lon_size
    endif

    if (lat_rank == lat_size-1) then
      above = MPI_PROC_NULL
    else
      above = mpi_rank + lon_size
    endif

    cnt = l * m * n * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      send_to_below(lc, mc, nc, i) = var(lc, mc, nc, mlat0, i+mlon0-1)
      send_to_above(lc, mc, nc, i) = var(lc, mc, nc, mlat1, i+mlon0-1)
    enddo

! sync in latitude
    call MPI_Isend(send_to_below, cnt, mpi_rp, &
      below, 0, dynamo_world, request(1), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Isend(send_to_above, cnt, mpi_rp, &
      above, 1, dynamo_world, request(2), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Irecv(recv_from_above, cnt, mpi_rp, &
      above, 0, dynamo_world, request(3), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

    call MPI_Irecv(recv_from_below, cnt, mpi_rp, &
      below, 1, dynamo_world, request(4), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! unpack to model fields
    if (lat_rank /= 0) then
      do concurrent (i = mlon0:mlon1, nc = 1:n, mc = 1:m, lc = 1:l)
        var(lc, mc, nc, mlatd0, i) = recv_from_below(lc, mc, nc, i-mlon0+1)
      enddo
    endif

    if (lat_rank /= lat_size-1) then
      do concurrent (i = mlon0:mlon1, nc = 1:n, mc = 1:m, lc = 1:l)
        var(lc, mc, nc, mlatd1, i) = recv_from_above(lc, mc, nc, i-mlon0+1)
      enddo
    endif

  endsubroutine sync_mlat_5d
!-----------------------------------------------------------------------
  subroutine sync_mlon_5d(var, l, m, n)

    use MPI

    integer, intent(in) :: l, m, n
    real(kind=rp), dimension(l, m, n, mlatd0:mlatd1, mlond0:mlond1), intent(inout) :: var

    integer :: j, lc, mc, nc, left, right, cnt, ierror
    integer, dimension(4) :: request
    real(kind=rp), dimension(l, m, n, maxmlat+4) :: &
      send_to_left, send_to_right, recv_from_left, recv_from_right

! find the rank of adjacent processes
    if (lon_rank == 0) then
      left = mpi_rank - 1 + lon_size
    else
      left = mpi_rank - 1
    endif

    if (lon_rank == lon_size-1) then
      right = mpi_rank + 1 - lon_size
    else
      right = mpi_rank + 1
    endif

    cnt = l * m * n * (maxmlat + 4)

! load to work array
    do concurrent (j = 1:mlatd1-mlatd0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      send_to_left(lc, mc, nc, j) = var(lc, mc, nc, j+mlatd0-1, mlon0)
      send_to_right(lc, mc, nc, j) = var(lc, mc, nc, j+mlatd0-1, mlon1)
    enddo

! sync in longitude
    call MPI_Isend(send_to_left, cnt, mpi_rp, &
      left, 0, dynamo_world, request(1), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Isend(send_to_right, cnt, mpi_rp, &
      right, 1, dynamo_world, request(2), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

    call MPI_Irecv(recv_from_right, cnt, mpi_rp, &
      right, 0, dynamo_world, request(3), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

    call MPI_Irecv(recv_from_left, cnt, mpi_rp, &
      left, 1, dynamo_world, request(4), ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)

! wait for sync to complete
    call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! unpack to model fields
    do concurrent (j = mlatd0:mlatd1, nc = 1:n, mc = 1:m, lc = 1:l)
      var(lc, mc, nc, j, mlond0) = recv_from_left(lc, mc, nc, j-mlatd0+1)
      var(lc, mc, nc, j, mlond1) = recv_from_right(lc, mc, nc, j-mlatd0+1)
    enddo

  endsubroutine sync_mlon_5d
!-----------------------------------------------------------------------
  function gather_mlon_3d(varin, m, n) result(varout)

    use MPI

    integer, intent(in) :: m, n
    real(kind=rp), dimension(m, n, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(m, n, nmlon) :: varout

    integer :: i, mc, nc, cnt, rnki, i0, i1, ierror
    integer, dimension(0:lon_size*2-1) :: request
    real(kind=rp), dimension(m, n, maxmlon) :: sendbuf
    real(kind=rp), dimension(m, n, maxmlon, 0:lon_size-1) :: recvbuf

    cnt = m * n * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, nc = 1:n, mc = 1:m)
      sendbuf(mc, nc, i) = varin(mc, nc, i+mlon0-1)
    enddo

! gather longitudes
    do rnki = 0, lon_size-1
      call MPI_Isend(sendbuf, cnt, mpi_rp, &
        lat_rank*lon_size + rnki, 0, dynamo_world, &
        request(rnki), ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Isend', ierror)

      call MPI_Irecv(recvbuf(:, :, :, rnki), cnt, mpi_rp, &
        lat_rank*lon_size + rnki, 0, dynamo_world, &
        request(lon_size+rnki), ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Irecv', ierror)
    enddo

! wait for gather to complete
    call MPI_Waitall(lon_size*2, request, MPI_STATUSES_IGNORE, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Waitall', ierror)

! reconstruct longitude rings
    do concurrent (rnki = 0:lon_size-1)
      i0 = mlon0_task(lat_rank*lon_size + rnki)
      i1 = mlon1_task(lat_rank*lon_size + rnki)
      do concurrent (i = i0:i1, nc = 1:n, mc = 1:m)
        varout(mc, nc, i) = recvbuf(mc, nc, i-i0+1, rnki)
      enddo
    enddo

  endfunction gather_mlon_3d
!-----------------------------------------------------------------------
  function gather_mag_2d(varin, root) result(varout)

    use MPI

    integer, intent(in) :: root
    real(kind=rp), dimension(mlat0:mlat1, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(nmlat, nmlon) :: varout

    integer :: i, j, cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1)
      sendbuf(j, i) = varin(j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mpi_rank==root) then
      do concurrent (rnk = 0:mpi_size-1)
        i0 = mlon0_task(rnk)
        i1 = mlon1_task(rnk)
        j0 = mlat0_task(rnk)
        j1 = mlat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1)
          varout(j, i) = recvbuf(j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif

  endfunction gather_mag_2d
!-----------------------------------------------------------------------
  function gather_mag_3d(varin, n, root) result(varout)

    use MPI

    integer, intent(in) :: n, root
    real(kind=rp), dimension(n, mlat0:mlat1, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(n, nmlat, nmlon) :: varout

    integer :: i, j, nc, cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(n, maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(n, maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = n * maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1, nc = 1:n)
      sendbuf(nc, j, i) = varin(nc, j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mpi_rank==root) then
      do concurrent (rnk = 0:mpi_size-1)
        i0 = mlon0_task(rnk)
        i1 = mlon1_task(rnk)
        j0 = mlat0_task(rnk)
        j1 = mlat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1, nc = 1:n)
          varout(nc, j, i) = recvbuf(nc, j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif

  endfunction gather_mag_3d
!-----------------------------------------------------------------------
  function gather_mag_4d(varin, m, n, root) result(varout)

    use MPI

    integer, intent(in) :: m, n, root
    real(kind=rp), dimension(m, n, mlat0:mlat1, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(m, n, nmlat, nmlon) :: varout

    integer :: i, j, mc, nc, cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(m, n, maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(m, n, maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = m * n * maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1, nc = 1:n, mc = 1:m)
      sendbuf(mc, nc, j, i) = varin(mc, nc, j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mpi_rank==root) then
      do concurrent (rnk = 0:mpi_size-1)
        i0 = mlon0_task(rnk)
        i1 = mlon1_task(rnk)
        j0 = mlat0_task(rnk)
        j1 = mlat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1, nc = 1:n, mc = 1:m)
          varout(mc, nc, j, i) = recvbuf(mc, nc, j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif

  endfunction gather_mag_4d
!-----------------------------------------------------------------------
  function gather_mag_5d(varin, l, m, n, root) result(varout)

    use MPI

    integer, intent(in) :: l, m, n, root
    real(kind=rp), dimension(l, m, n, mlat0:mlat1, mlon0:mlon1), intent(in) :: varin
    real(kind=rp), dimension(l, m, n, nmlat, nmlon) :: varout

    integer :: i, j, lc, mc, nc, cnt, rnk, i0, i1, j0, j1, ierror
    real(kind=rp), dimension(l, m, n, maxmlat, maxmlon) :: sendbuf
    real(kind=rp), dimension(l, m, n, maxmlat, maxmlon, 0:mpi_size-1) :: recvbuf

    cnt = l * m * n * maxmlat * maxmlon

! load to work array
    do concurrent (i = 1:mlon1-mlon0+1, j = 1:mlat1-mlat0+1, nc = 1:n, mc = 1:m, lc = 1:l)
      sendbuf(lc, mc, nc, j, i) = varin(lc, mc, nc, j+mlat0-1, i+mlon0-1)
    enddo

    if (root < 0) then
      call MPI_Allgather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allgather', ierror)
    else
      call MPI_Gather(sendbuf, cnt, mpi_rp, &
        recvbuf, cnt, mpi_rp, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Gather', ierror)
    endif

! reconstruct based on lat-lon decomposition
    if (root<0 .or. mpi_rank==root) then
      do concurrent (rnk = 0:mpi_size-1)
        i0 = mlon0_task(rnk)
        i1 = mlon1_task(rnk)
        j0 = mlat0_task(rnk)
        j1 = mlat1_task(rnk)
        do concurrent (i = i0:i1, j = j0:j1, nc = 1:n, mc = 1:m, lc = 1:l)
          varout(lc, mc, nc, j, i) = recvbuf(lc, mc, nc, j-j0+1, i-i0+1, rnk)
        enddo
      enddo
    endif

  endfunction gather_mag_5d
!-----------------------------------------------------------------------
  subroutine bcast_2d(var, m, n, root)

    use MPI

    integer, intent(in) :: m, n, root
    real(kind=rp), dimension(m, n), intent(inout) :: var

    integer :: cnt, mc, nc, ierror
    real(kind=rp), dimension(m, n) :: buffer

    cnt = m * n

! load to work array
    do concurrent (nc = 1:n, mc = 1:m)
      buffer(mc, nc) = var(mc, nc)
    enddo

    call MPI_Bcast(buffer, cnt, mpi_rp, root, dynamo_world, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Bcast', ierror)

! unpack to model fields
    do concurrent (nc = 1:n, mc = 1:m)
      var(mc, nc) = buffer(mc, nc)
    enddo

  endsubroutine bcast_2d
!-----------------------------------------------------------------------
  subroutine bcast_3d(var, l, m, n, root)

    use MPI

    integer, intent(in) :: l, m, n, root
    real(kind=rp), dimension(l, m, n), intent(inout) :: var

    integer :: cnt, lc, mc, nc, ierror
    real(kind=rp), dimension(l, m, n) :: buffer

    cnt = l * m * n

! load to work array
    do concurrent (nc = 1:n, mc = 1:m, lc = 1:l)
      buffer(lc, mc, nc) = var(lc, mc, nc)
    enddo

    call MPI_Bcast(buffer, cnt, mpi_rp, root, dynamo_world, ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Bcast', ierror)

! unpack to model fields
    do concurrent (nc = 1:n, mc = 1:m, lc = 1:l)
      var(lc, mc, nc) = buffer(lc, mc, nc)
    enddo

  endsubroutine bcast_3d
!-----------------------------------------------------------------------
  function reduce_sum_1d(varin, n, root) result(varout)

    use MPI

    integer, intent(in) :: n, root
    real(kind=rp), dimension(n), intent(in) :: varin
    real(kind=rp), dimension(n) :: varout

    integer :: ierror

    if (root < 0) then
      call MPI_Allreduce(varin, varout, n, mpi_rp, &
        MPI_SUM, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Allreduce', ierror)
    else
      call MPI_Reduce(varin, varout, n, mpi_rp, &
        MPI_SUM, root, dynamo_world, ierror)
      if (ierror /= MPI_SUCCESS) call handle_error('MPI_Reduce', ierror)
    endif

  endfunction reduce_sum_1d
!-----------------------------------------------------------------------
  subroutine finalize

    use MPI

    integer :: ierror

    call MPI_Finalize(ierror)
    if (ierror /= MPI_SUCCESS) call handle_error('MPI_Finalize', ierror)

  endsubroutine finalize
!-----------------------------------------------------------------------
  subroutine handle_error(funcname, errorcode)

    use MPI

    character(len=*), intent(in) :: funcname
    integer, intent(in) :: errorcode

    character(len=MPI_MAX_ERROR_STRING) :: string
    integer :: resultlen, ierror

    call MPI_Error_string(errorcode, string, resultlen, ierror)
    write(6, "('MPI error encountered: ', a, ', when calling ', a, '. Finalizing...')") &
      trim(string), trim(funcname)
    call MPI_Finalize(ierror)

  endsubroutine handle_error
!-----------------------------------------------------------------------
  pure function generate_minvar_list(summation, nfactor) result(factor_list)
! generate a list of the given length summing up to the given number
! which has the minimum variances (all factors differ at most by 1)
! e.g., 5 out of 26 would be 5,5,5,5,6

    integer, intent(in) :: summation, nfactor
    integer, dimension(nfactor) :: factor_list

    integer :: minfactor, maxfactor, mincnt
    real(kind=rp) :: factor

    factor = real(summation, kind=rp) / nfactor
    minfactor = floor(factor)
    maxfactor = ceiling(factor)

    mincnt = maxfactor*nfactor - summation

    factor_list(1:mincnt) = minfactor
    factor_list(mincnt+1:nfactor) = maxfactor

  endfunction generate_minvar_list
!-----------------------------------------------------------------------
endmodule mpi_module
