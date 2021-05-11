module aero_check_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none

  private
  public :: aero_check_routine
  public :: aero_check_errors

  abstract interface
     subroutine check_template(q, ncol, nlev, packagename, fix, logmsg, abort)
       use shr_kind_mod, only: r8 => shr_kind_r8
       real(r8), intent(inout) :: q(:,:,:)
       integer, intent(in) :: ncol, nlev
       character(len=*), intent(in) :: packagename
       logical, intent(in) :: fix
       logical, intent(in) :: logmsg
       logical, intent(in) :: abort
     end subroutine check_template
  end interface

  procedure(check_template), pointer :: aero_check_routine=>null()

contains

  subroutine aero_check_errors( q, ncol, nlev, msg, fix, logmsg, abort)
    real(r8), intent(inout) :: q(:,:,:)
    integer, intent(in) :: ncol, nlev
    character(len=*), intent(in) :: msg
    logical, intent(in) :: fix
    logical, intent(in) :: logmsg
    logical, intent(in) :: abort

    call aero_check_routine( q, ncol, nlev, msg, fix, logmsg, abort )

  end subroutine aero_check_errors

end module aero_check_mod
