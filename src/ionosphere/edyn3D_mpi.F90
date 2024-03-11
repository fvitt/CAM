module edyn3D_mpi
   use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun

   use spmd_utils,     only: masterproc
   use mpi,            only: mpi_comm_size, mpi_comm_rank, mpi_comm_split
   use mpi,            only: MPI_PROC_NULL, mpi_wait, MPI_STATUS_SIZE
   use mpi,            only: MPI_INTEGER, MPI_REAL8, MPI_SUCCESS, MPI_SUM

   implicit none
   private
   save

   ! Public data
   public :: mlon0_p, mlon1_p
   public :: mytid
   public :: ntask
   public :: nmagtasklon
   public :: tasks
   ! Public type
   public :: array_ptr_type
   ! Public interfaces
   public :: mp_init_edyn3D
!   public :: mp_scatter_phim
!   public :: mp_mag_foldhem
!   public :: mp_gather_edyn
   public :: ixfind
!   public :: mp_magpole_3d
!   public :: setpoles
!   public :: mp_gatherlons_f3d
!   public :: mp_scatterlons_f3d
   public :: mp_exchange_tasks_edyn3D
   public :: mp_distribute_mag_edyn3D
   public :: mlat0_p,mlat1_p
   public :: mlat0_r,mlat1_r
   public :: mlat0_s1,mlat1_s1
   public :: mlat0_s2,mlat1_s2
   public :: fldpts0_qdlat_p,fldpts1_qdlat_p
   public :: fldpts0_qdlat_r,fldpts1_qdlat_r
   public :: fldpts0_qdlat_s1,fldpts1_qdlat_s1
   public :: fldpts0_qdlat_s2,fldpts1_qdlat_s2
   public :: fldpts0_p,fldpts1_p
   public :: fldpts0_r,fldpts1_r
   public :: fldpts0_s1,fldpts1_s1
   public :: fldpts0_s2,fldpts1_s2

   !
   ! Number of MPI tasks and current task id:
   !
   integer :: &
        ntask,   & ! number of mpi tasks
        mytid      ! my task id
!   !
!   ! Magnetic computational subdomains for current task:
!   !
!   integer, protected :: &
!	 ncompmagtasklon_p,  & ! number of tasks in computation magnetic longitude dimension
!	 ncompmagtasklat_p,  & ! number of tasks in computation magnetic latitude dimension
!	 compmagtidlon_p,    & ! computation longitude coord for current task in task table
!	 compmagtidlat_p,    & ! computation latitude coord for current task in task table
!	 mclon_p0=1,mclon_p1=0,& ! first and last computation mag lons for each task
!	 mclat_p0=1,mclat_p1=0   ! first and last computation mag lats for each task
   !
   ! Magnetic output subdomains for current task:
   !
   integer, protected :: &
        nmagtasklon,   &       ! number of tasks in magnetic longitude dimension
        magtidlon,     &       ! longitude coord for current task in task table
        mlon0_p=1,mlon1_p=0    ! first and last mag lons for each task

   integer :: &
        mlat0_p,mlat1_p,     &
	mlat0_r,mlat1_r,     &
	mlat0_s1,mlat1_s1,   &
	mlat0_s2,mlat1_s2,   &
        fldpts0_qdlat_p,fldpts1_qdlat_p,   &
        fldpts0_qdlat_r,fldpts1_qdlat_r,   &
        fldpts0_qdlat_s1,fldpts1_qdlat_s1, &
        fldpts0_qdlat_s2,fldpts1_qdlat_s2, &
        fldpts0_p,fldpts1_p,	&
        fldpts0_r,fldpts1_r,	&
        fldpts0_s1,fldpts1_s1,  &
        fldpts0_s2,fldpts1_s2

!   integer :: &
!        mxcmaglon_p,  &   ! max number of computation mag subdomain lon points among all tasks
!        mxcmaglat_p,  &   ! max number of computation mag subdomain lat points among all tasks
!        mxhmaglon       ! max number of output mag subdomain lon points among all tasks

   integer, allocatable :: &
!        itask_table_comp_mag(:),  &  ! 1d table of tasks on mag computation grid (i)
        itask_table_mag(:)      ! 1d table of tasks on mag output grid (i)

!   integer :: cols_comm_comp_p_edyn3D  ! communicators for each computation task column
!   integer :: rows_comm_comp_p_edyn3D  ! communicators for each computation task row
   integer :: cols_comm_edyn3D  ! communicators for each output task column
   !
   ! Task type: subdomain information for all tasks, known by all tasks:
   !
   type task
      integer :: mytid       ! task id
      !
      ! Magnetic subdomains in task structure:
      !
      integer :: magtidlon = -1       ! task coord in output mag longitude dimension of task table
      integer :: nmaglons = 0            ! number of output mag longitudes calculated by this task
      integer :: mlon0 = 1,mlon1 = 0 ! first and last output longitude indices
   end type task
   !
   ! type(task) :: tasks(ntask) will be made available to all tasks
   ! (so each task has information about all tasks)
   !
   type(task), allocatable :: tasks(:)
   !
   ! Magnetic grid parameters
   !
!   integer :: ncmlon_p   ! number of computation longitudes
!   integer :: ncmlat_p   ! number of computation latitudes
   integer :: nmlon   ! number of output longitudes

!   type array_ptr_comp_type
!      real(r8),pointer :: ptr(:,:,:) ! (i,j,flp)
!   end type array_comp_ptr_type

   type array_ptr_type
      real(r8),pointer :: ptr(:,:) ! (i,flp)
   end type array_ptr_type

   integer, protected :: mpi_comm_edyn3D = -9999

   logical, parameter :: debug = .false.
!   logical, parameter :: debug = .true.

contains
   !-----------------------------------------------------------------------
   subroutine mp_init_edyn3D(mpi_comm, ionos_npes)
      !
      ! Initialize MPI, and allocate task table.
      !
      integer, intent(in) :: mpi_comm
      integer, intent(in) :: ionos_npes

      integer :: mpi_comm_edyn3D

      integer :: ierr
      integer :: color, npes
      character(len=cl) :: errmsg

      if (masterproc) then
         write(iulog,*) 'mpi_init_edyn3D... ionos_npes :', ionos_npes
      endif

      call mpi_comm_size(mpi_comm, npes, ierr)

      ntask = min(npes,ionos_npes)

      if (masterproc) then
         write(iulog,*) 'mpi_init_edyn3D: Initializing with ntask = ', ntask
      endif

      call mpi_comm_rank(mpi_comm, mytid, ierr)
      color = mytid/ntask !  ionos_npes
      call mpi_comm_split(mpi_comm, color, mytid, mpi_comm_edyn3D, ierr)
      !
      ! Allocate array of task structures:
      !
      allocate(tasks(0:npes-1), stat=ierr)
      if (ierr /= 0) then
         write(errmsg,"('>>> mp_init: error allocating tasks(',i3,')')") ntask
         write(iulog,*) trim(errmsg)
         call endrun(errmsg)
      endif
   end subroutine mp_init_edyn3D
   !-----------------------------------------------------------------------
   subroutine mp_distribute_mag_edyn3D(nmlon_in)
      !
      ! Args:
      integer, intent(in) :: nmlon_in ! number of longitudes
      !
      ! Local:
      integer                     :: i, n, irank, ier, nmaglon_p, tidcol, ncells
      character(len=cl)          :: errmsg
      character(len=*), parameter :: subname = 'mp_distribute_edyn3D_mag'
      !
      ! Number of tasks in mag lon:
      !
      nmagtasklon = ntask
      !
      ! Store magnetic grid parameters
      nmlon = nmlon_in

      if (mytid<ntask) then
         !
         ! Allocate and set 1d table of tasks:
	 !
         allocate(itask_table_mag(-1:nmagtasklon),stat=ier)
         if (ier /= 0) then
            write(errmsg, "(a,2(a,i0))") subname,                                &
                 '>>> Error allocating itable: nmagtasklon = ', nmagtasklon
            if (masterproc) then
               write(iulog, errmsg)
            end if
            call endrun(errmsg)
         endif

         itask_table_mag(:) = MPI_PROC_NULL

         irank = 0
         do i = 0, nmagtasklon-1
           itask_table_mag(i) = irank
           if (mytid == irank) then
              magtidlon = i
           endif
           irank = irank + 1
         end do

!         if (debug .and. masterproc) then
         if (debug) then
            !
            ! Print table to stdout:
            write(*,"(/,a,/a,i5,a,i5,a,i5,' Mag Task Table:')") subname,     &
                 'ntask=',ntask,' nmagtasklon=',nmagtasklon
            do i = 1, nmagtasklon
               write(iulog,"('i = ',i5,', itask_table_mag(i) = ',i5)") i,itask_table_mag(i)
            end do
         end if
         !
         ! Calculate start and end indices in mag longitude dimensions for each task:
         !
         call distribute_1d(1, nmlon, nmagtasklon, magtidlon, mlon0_p, mlon1_p)

         nmaglon_p = mlon1_p - mlon0_p + 1 ! number of mag longitudes for this task

         if (debug) then
            !
            ! Report my stats to stdout:
	    !
            write(*,"(/,a,i5,a,i5,a,2i5,a,i5)")            &
                 'mytid = ',mytid, ', magtidlon = ', magtidlon,                  &
                 'mlon0_p,1 = ', mlon0_p, mlon1_p, ' (', nmaglon_p, ') '
         end if
         !
         ! Define all task structures with current task values
         ! (redundant for alltoall):
         !
         do n=0,ntask-1
            tasks(n)%mytid     = mytid
            tasks(n)%magtidlon = magtidlon
            tasks(n)%nmaglons  = nmaglon_p
            tasks(n)%mlon0     = mlon0_p
            tasks(n)%mlon1     = mlon1_p
         enddo
         !
         ! All tasks must have at least 4 longitudes:
	 !
         do n = 0, ntask-1

            if (debug) then
               write(iulog,"('mp_distribute_mag_edyn3D: n=',i5,' tasks(n)%nmaglons=',i5)") &
                    n,tasks(n)%nmaglons
            endif

            if (tasks(n)%nmaglons < 4) then
               write(errmsg, "(3a,i0,', nmaglons = ',i5)") '>>> ', subname,      &
                    ': each task must carry at least 4 longitudes. task = ',     &
                    n, tasks(n)%nmaglons
               if (masterproc) then
                  write(iulog, errmsg)
               end if
               call endrun(errmsg)
            end if
         end do

      end if

   end subroutine mp_distribute_mag_edyn3D
   !-----------------------------------------------------------------------
   subroutine distribute_1d(n1,n2,nprocs,myrank,istart,iend)
      !
      ! Distribute work across a 1d vector(n1->n2) to nprocs.
      ! Return start and end indices for proc myrank.
      !
      ! Args:
      integer,intent(in)  :: n1,n2,nprocs,myrank
      integer,intent(out) :: istart,iend
      !
      ! Local:
      integer :: lenproc,iremain,n
      !
      n = n2-n1+1
      lenproc = n/nprocs
      iremain = mod(n,nprocs)
      istart = n1 + myrank*lenproc + min(myrank,iremain)
      iend = istart+lenproc-1
      if (iremain > myrank) iend = iend+1
   end subroutine distribute_1d
   !-----------------------------------------------------------------------
   subroutine mp_exchange_tasks_edyn3D(mpi_comm, iprint)
      !
      ! Args:
      integer,  intent(in) :: mpi_comm
      integer,  intent(in) :: iprint
      !
      ! Local:
      ! itasks_send(len_task_type,ntask) will be used to send tasks(:) info
      !   to all tasks (directly passing mpi derived data types is reportedly
      !   not stable, or not available until MPI 2.x).
      !
      integer :: n, ier
      integer, parameter :: len_task_type = 5 ! see type task above
      integer, allocatable :: &
           itasks_send(:,:), & ! send buffer
           itasks_recv(:,:)    ! send buffer
      integer :: npes

      call mpi_comm_size(mpi_comm, npes, ier)
      !
      ! Pack tasks(mytid) into itasks_send:
      !
      allocate(itasks_send(len_task_type,0:npes-1),stat=ier)
      if (ier /= 0) then
         write(iulog,"(i4,i4)") '>>> Error allocating itasks_send: len_task_type=',&
              len_task_type,' npes=',npes
         call endrun('mp_exchange_tasks: unable to allocate itasks_send')
      endif
      allocate(itasks_recv(len_task_type,0:npes-1),stat=ier)
      if (ier /= 0) then
         write(iulog,"(i4,i4)") '>>> Error allocating itasks_recv: len_task_type=',&
              len_task_type,' npes=',npes
         call endrun('mp_exchange_tasks: unable to allocate itasks_recv')
      endif
      do n=0,npes-1
         itasks_send(1,n) = tasks(mytid)%mytid

         itasks_send(2,n) = tasks(mytid)%magtidlon
         itasks_send(3,n) = tasks(mytid)%nmaglons
         itasks_send(4,n) = tasks(mytid)%mlon0
         itasks_send(5,n) = tasks(mytid)%mlon1
      enddo
      !
      ! Send itasks_send and receive itasks_recv:
      !
      call mpi_alltoall(itasks_send,len_task_type,MPI_INTEGER,&
           itasks_recv,len_task_type,MPI_INTEGER,&
           mpi_comm,ier)
      if (ier /= 0) &
           call handle_mpi_err(ier,'edyn_mpi: mpi_alltoall to send/recv itasks')
      !
      ! Unpack itasks_recv into tasks(n)
      !
      do n=0,npes-1
         tasks(n)%mytid  = itasks_recv(1,n)

         tasks(n)%magtidlon   = itasks_recv(2,n)
         tasks(n)%nmaglons  = itasks_recv(3,n)
         tasks(n)%mlon0   = itasks_recv(4,n)
         tasks(n)%mlon1   = itasks_recv(5,n)
         !
         ! Report to stdout:
         !
         if (n==mytid.and.iprint > 0) then
            write(*,"(/,'Task ',i3,':')") n
            write(*,"(/,'Subdomain on geomagnetic longitude grid:')")
            write(*,"('tasks(',i3,')%magtidlon=',i3)") n,tasks(n)%magtidlon
            write(*,"('tasks(',i3,')%nmaglons =',i3)") n,tasks(n)%nmaglons
            write(*,"('tasks(',i3,')%mlon0  =',i3)") n,tasks(n)%mlon0
            write(*,"('tasks(',i3,')%mlon1  =',i3)") n,tasks(n)%mlon1
            write(*,"('Number of mag subdomain grid points = ',i6)") &
                 tasks(n)%nmaglons
         endif
      enddo
      !
      ! Release locally allocated space:
      deallocate(itasks_send)
      deallocate(itasks_recv)

   end subroutine mp_exchange_tasks_edyn3D
   !-----------------------------------------------------------------------
   subroutine handle_mpi_err(ierrcode,string)
      !
      ! Args:
      integer,intent(in) :: ierrcode
      character(len=*) :: string
      !
      ! Local:
      character(len=80) :: errstring
      integer :: len_errstring, ierr
      !
      call mpi_error_string(ierrcode,errstring,len_errstring, ierr)
      write(iulog,"(/,'>>> mpi error: ',a)") trim(string)
      write(iulog,"('  ierrcode=',i3,': ',a)") trim(errstring)
   end subroutine handle_mpi_err
   !-----------------------------------------------------------------------
   integer function ixfind(iarray,idim,itarget,icount)
      !
      ! Search iarray(idim) for itarget, returning first index in iarray
      ! where iarray(idim)==target. Also return number of elements of
      ! iarray that == itarget in icount.
      !
      ! Args:
      integer,intent(in) :: idim,itarget
      integer,intent(in) :: iarray(idim)
      integer,intent(out) :: icount
      !
      ! Local:
      integer :: i
      !
      ixfind = 0
      icount = 0
      if (.not.any(iarray==itarget)) return
      icount = count(iarray==itarget)
      do i=1,idim
         if (iarray(i)==itarget) then
            ixfind = i
            exit
         endif
      enddo
   end function ixfind

end module edyn3D_mpi
