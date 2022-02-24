module aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: aerosol_properties

  type, abstract :: aerosol_properties
   contains
     procedure(aero_props_abdraz_f), deferred :: abdraz_f1
     procedure(aero_props_abdraz_f), deferred :: abdraz_f2
     procedure(aero_props_get), deferred :: get
  end type aerosol_properties

  interface
     subroutine aero_props_get(self, m,l, density,hygro)
       import
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: m,l
       real(r8), optional, intent(out) :: density
       real(r8), optional, intent(out) :: hygro
     end subroutine aero_props_get

     function aero_props_abdraz_f(self, m) result(f)
       import
       class(aerosol_properties), intent(in) :: self
       integer, intent(in) :: m
       real(r8) :: f
     end function aero_props_abdraz_f
  end interface

contains

end module aerosol_properties_mod
