!
! Manages running time average TTEND_DP from ZM deep convection scheme
!
module running_tave_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
  use physics_buffer, only: physics_buffer_desc, pbuf_set_field, pbuf_add_field, pbuf_get_field, dtype_r8
  use ppgrid, only: pver, pcols, pverp

  implicit none

  integer,parameter :: accmax = 10 ! could make this a run-time namelist option

  integer :: accnum = 0
  integer :: data_idx = -1
  integer :: tave_idx = -1

contains

  ! register physics buffer fields
  subroutine running_tave_reg()

    ! holds the accumulated data
    call pbuf_add_field('TTEND_DP_DATA','global', dtype_r8, (/ pcols,pver,accmax /), data_idx )

    ! the running time average
    call pbuf_add_field('TTEND_DP_TAVE','global', dtype_r8, (/ pcols,pver /), tave_idx )

  end subroutine running_tave_reg

  ! initialize the physics buffer fields to zero at the beginning of the run
  subroutine running_tave_init(pbuf2d)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    call pbuf_set_field(pbuf2d, data_idx, 0._r8)
    call pbuf_set_field(pbuf2d, tave_idx, 0._r8)

  end subroutine running_tave_init

  ! update the running time average which is stored in the physics buffer
  subroutine running_tave_update( ttend, ncol, pbuf )
    real(r8), intent(in) :: ttend(:,:)
    integer,  intent(in) :: ncol
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: n
    real(r8), pointer :: tave_ptr(:,:)
    real(r8), pointer :: data_ptr(:,:,:)

    call pbuf_get_field(pbuf, data_idx, data_ptr)
    call pbuf_get_field(pbuf, tave_idx, tave_ptr)

    accnum = accnum + 1
    accnum = min(accmax,accnum)

    tave_ptr = 0._r8

    ! manage the storage of the accumulated ttend arrays
    do n = 1,accmax-1
       ! shift the previous ttend entries
       data_ptr(:ncol,:,n+1) = data_ptr(:ncol,:,n)
    end do
    ! store current ttend at position 1
    data_ptr(:ncol,:,1) = ttend(:ncol,:)

    ! compute running time average
    do n = 1,accnum
       tave_ptr(:ncol,:) = tave_ptr(:ncol,:) + data_ptr(:ncol,:,n)
    end do
    tave_ptr(:ncol,:) = tave_ptr(:ncol,:)/accnum

  end subroutine running_tave_update

end module running_tave_mod
