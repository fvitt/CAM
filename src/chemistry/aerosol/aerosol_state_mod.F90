module aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: aerosol_state

  type, abstract :: aerosol_state
   contains
     procedure(aero_get_state_mmr), deferred :: get_ambient_mmr
     procedure(aero_get_state_mmr), deferred :: get_cldbrne_mmr
     procedure(aero_get_state_num), deferred :: get_ambient_num
     procedure(aero_get_state_num), deferred :: get_cldbrne_num
  end type aerosol_state

  interface
     function aero_get_state_mmr(self, l,m) result(x)
       import
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: l
       integer, intent(in) :: m
       real(r8), pointer :: x(:,:)
     end function aero_get_state_mmr

     function aero_get_state_num(self, m) result(x)
       import
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: m
       real(r8), pointer :: x(:,:)
     end function aero_get_state_num
  end interface

contains

end module aerosol_state_mod
