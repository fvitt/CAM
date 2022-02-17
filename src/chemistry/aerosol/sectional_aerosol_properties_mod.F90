module sectional_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: sectional_aerosol_properties

  type, extends(aerosol_properties) :: sectional_aerosol_properties
     private
   contains
     procedure :: abdraz_f1
     procedure :: abdraz_f2
     final :: destructor
  end type sectional_aerosol_properties

  interface sectional_aerosol_properties
     procedure :: constructor
  end interface sectional_aerosol_properties

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor() result(newobj)
    type(sectional_aerosol_properties), pointer :: newobj
    allocate(newobj)
  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(sectional_aerosol_properties), intent(inout) :: self

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function abdraz_f1(self,m) result(f)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8) :: f

    f = 1._r8
  end function abdraz_f1

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function abdraz_f2(self,m) result(f)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8) :: f

    f = 1._r8
  end function abdraz_f2

end module sectional_aerosol_properties_mod
