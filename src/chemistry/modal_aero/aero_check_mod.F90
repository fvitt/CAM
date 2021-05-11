module aero_check_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

contains

  subroutine aero_check_errors( q, ncol, nlev, msg, fix, logmsg, abort)
    real(r8), intent(inout) :: q(:,:,:)
    integer, intent(in) :: ncol, nlev
    character(len=*), intent(in) :: msg
    logical, intent(in) :: fix
    logical, intent(in) :: logmsg
    logical, intent(in) :: abort

  end subroutine aero_check_errors

end module aero_check_mod
