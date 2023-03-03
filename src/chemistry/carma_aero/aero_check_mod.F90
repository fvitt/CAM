module aero_check_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8

  implicit none

  private
  public :: aero_check_routine
  public :: aero_check_errors

  !! NOTE: I would like to pass in state, but there is a circular
  !! dependency with physics types.
  abstract interface
     subroutine check_template(pdel, q, lchnk, ncol, nlev, packagename, fix, logmsg, abort)
       use shr_kind_mod, only: r8 => shr_kind_r8
       real(r8), intent(in) :: pdel(:,:)
       real(r8), intent(inout) :: q(:,:,:)
       integer, intent(in) :: lchnk, ncol, nlev
       character(len=*), intent(in) :: packagename
       logical, intent(in) :: fix
       logical, intent(in) :: logmsg
       logical, intent(in) :: abort
     end subroutine check_template
  end interface

  procedure(check_template), pointer :: aero_check_routine=>null()

contains

  subroutine aero_check_errors( pdel, q, lchnk, ncol, nlev, msg, fix, logmsg, abort)
    real(r8), intent(in) :: pdel(:,:)
    real(r8), intent(inout) :: q(:,:,:)
    integer, intent(in) :: lchnk, ncol, nlev
    character(len=*), intent(in) :: msg
    logical, intent(in) :: fix
    logical, intent(in) :: logmsg
    logical, intent(in) :: abort

    call aero_check_routine( pdel, q, lchnk, ncol, nlev, msg, fix, logmsg, abort )

  end subroutine aero_check_errors

end module aero_check_mod
