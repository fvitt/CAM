module modal_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_properties_mod, only: aerosol_properties

  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_props

  implicit none

  private

  public :: modal_aerosol_properties

  type, extends(aerosol_properties) :: modal_aerosol_properties
     private
     integer :: nmodes
     real(r8), allocatable :: alogsig(:)
     real(r8), allocatable :: sigmag_amode(:)
     real(r8), allocatable :: f1(:)
     real(r8), allocatable :: f2(:)
   contains
     procedure :: abdraz_f1
     procedure :: abdraz_f2
     procedure :: get
     final :: destructor
  end type modal_aerosol_properties

  interface modal_aerosol_properties
     procedure :: constructor
  end interface modal_aerosol_properties

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor() result(newobj)
    type(modal_aerosol_properties), pointer :: newobj

    integer :: m

    allocate(newobj)

    call rad_cnst_get_info(0, nmodes=newobj%nmodes)

    allocate(newobj%alogsig(newobj%nmodes))
    allocate(newobj%sigmag_amode(newobj%nmodes))
    allocate(newobj%f1(newobj%nmodes))
    allocate(newobj%f2(newobj%nmodes))

    do m = 1, newobj%nmodes
       call rad_cnst_get_mode_props(0, m, sigmag=newobj%sigmag_amode(m))
       newobj%alogsig(m) = log(newobj%sigmag_amode(m))
       newobj%f1(m) = 0.5_r8*exp(2.5_r8*newobj%alogsig(m)*newobj%alogsig(m))
       newobj%f2(m) = 1._r8 + 0.25_r8*newobj%alogsig(m)
    end do

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_properties), intent(inout) :: self

    deallocate(self%alogsig)
    deallocate(self%sigmag_amode)
    deallocate(self%f1)
    deallocate(self%f2)

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get(self, m,l, density,hygro)
    use rad_constituents, only: rad_cnst_get_aer_props

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m,l
    real(r8), optional, intent(out) :: density
    real(r8), optional, intent(out) :: hygro

    call rad_cnst_get_aer_props(0, m, l, density_aer=density, hygro_aer=hygro)

  end subroutine get

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function abdraz_f1(self,m) result(f)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8) :: f

    f = self%f1(m)
  end function abdraz_f1

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function abdraz_f2(self,m) result(f)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8) :: f

    f = self%f2(m)
  end function abdraz_f2

end module modal_aerosol_properties_mod
