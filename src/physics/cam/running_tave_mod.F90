!================================================================================
!
! Manages running time average TTEND_DP from ZM deep convection scheme
!
!================================================================================
module running_tave_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
  use physics_buffer, only: physics_buffer_desc, pbuf_set_field, pbuf_add_field, pbuf_get_field, dtype_r8
  use ppgrid, only: pver, pcols, pverp
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc
  use time_manager, only: is_first_step, get_nstep

  implicit none

  integer :: accmax = 0 ! set by run-time namelist option

  integer :: accnum = 0
  integer :: data_idx = -1
  integer :: tave_idx = -1

contains

  !================================================================================
  ! read run-time namelist options
  !================================================================================
  subroutine running_tave_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use cam_abortutils, only: endrun
    use spmd_utils,     only: mpicom, masterprocid, mpi_integer, mpi_success

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'running_tave_readnl'

    integer :: running_tave_nsteps = 0

    namelist /running_tave_opts/ running_tave_nsteps

    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'running_tave_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, running_tave_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    call mpi_bcast(running_tave_nsteps, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: running_tave_nsteps')
    end if

    if (masterproc) then
       write(iulog, *) subname,': running_tave_nsteps : ',running_tave_nsteps
    end if

    accmax = running_tave_nsteps

  end subroutine running_tave_readnl

  !================================================================================
  ! register physics buffer fields
  !================================================================================
  subroutine running_tave_reg()

    ! holds the accumulated data
    call pbuf_add_field('TTEND_DP_DATA','global', dtype_r8, (/ pcols,pver,accmax /), data_idx )

    ! the running time average
    call pbuf_add_field('TTEND_DP_TAVE','global', dtype_r8, (/ pcols,pver /), tave_idx )

  end subroutine running_tave_reg

  !================================================================================
  ! initialize the physics buffer fields to zero at the beginning of the run
  !================================================================================
  subroutine running_tave_init(pbuf2d)

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (is_first_step()) then
       call pbuf_set_field(pbuf2d, data_idx, 0._r8)
       call pbuf_set_field(pbuf2d, tave_idx, 0._r8)
    end if

  end subroutine running_tave_init

  !================================================================================
  ! update the running time average which is stored in the physics buffer
  !================================================================================
  subroutine running_tave_update( ttend, ncol, pbuf )
    real(r8), intent(in) :: ttend(:,:)
    integer,  intent(in) :: ncol
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: n
    real(r8), pointer :: tave_ptr(:,:)
    real(r8), pointer :: data_ptr(:,:,:)
    integer :: nstep

    nstep = get_nstep()
    accnum = min(accmax,nstep+1)

    call pbuf_get_field(pbuf, data_idx, data_ptr)
    call pbuf_get_field(pbuf, tave_idx, tave_ptr)

    if (masterproc) then
       write(iulog,*) 'running_tave_update  nstep,accnum,accmax : ',nstep,accnum,accmax
    end if

    ! manage the storage of the accumulated ttend arrays
    do n = accmax, 2, -1
       ! shift the previous ttend entries
       data_ptr(:ncol,:,n) = data_ptr(:ncol,:,n-1)
    end do
    ! store current ttend at position 1
    data_ptr(:ncol,:,1) = ttend(:ncol,:)

    tave_ptr = 0.0_r8

    ! compute running time average
    do n = 1,accnum
       tave_ptr(:ncol,:) = tave_ptr(:ncol,:) + data_ptr(:ncol,:,n)
    end do
    tave_ptr(:ncol,:) = tave_ptr(:ncol,:)/accnum

  end subroutine running_tave_update

end module running_tave_mod
