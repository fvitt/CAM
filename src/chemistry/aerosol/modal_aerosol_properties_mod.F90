module modal_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: modal_aerosol_properties

  type, extends(aerosol_properties) :: modal_aerosol_properties
     private
     real(r8), allocatable :: alogsig(:)
     real(r8), allocatable :: sigmag_amode(:)
     real(r8), allocatable :: f1(:)
     real(r8), allocatable :: f2(:)
     real(r8), allocatable :: voltonumblo_(:)
     real(r8), allocatable :: voltonumbhi_(:)
   contains
     procedure :: abdraz_f1
     procedure :: abdraz_f2
     procedure :: get
     procedure :: voltonumblo
     procedure :: voltonumbhi
     final :: destructor
  end type modal_aerosol_properties

  interface modal_aerosol_properties
     procedure :: constructor
  end interface modal_aerosol_properties

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor() result(newobj)
    use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_props
    use physconst,        only: pi

    type(modal_aerosol_properties), pointer :: newobj

    integer :: m, nmodes
    real(r8) :: dgnumlo
    real(r8) :: dgnumhi
    integer,allocatable :: nspecies(:)

    allocate(newobj)

    call rad_cnst_get_info(0, nmodes=nmodes)

    allocate(nspecies(nmodes))
    allocate(newobj%alogsig(nmodes))
    allocate(newobj%sigmag_amode(nmodes))
    allocate(newobj%f1(nmodes))
    allocate(newobj%f2(nmodes))
    allocate(newobj%voltonumblo_(nmodes))
    allocate(newobj%voltonumbhi_(nmodes))

    do m = 1, nmodes
       call rad_cnst_get_info(0, m, nspec=nspecies(m))
       call rad_cnst_get_mode_props(0, m, sigmag=newobj%sigmag_amode(m),  &
            dgnumhi=dgnumhi, dgnumlo=dgnumlo )
       newobj%alogsig(m) = log(newobj%sigmag_amode(m))
       newobj%f1(m) = 0.5_r8*exp(2.5_r8*newobj%alogsig(m)*newobj%alogsig(m))
       newobj%f2(m) = 1._r8 + 0.25_r8*newobj%alogsig(m)

       newobj%voltonumblo_(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumlo**3._r8)*exp(4.5_r8*newobj%alogsig(m)**2._r8) )
       newobj%voltonumbhi_(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumhi**3._r8)*exp(4.5_r8*newobj%alogsig(m)**2._r8) )

    end do

    call newobj%initialize(nmodes,nspecies)
    deallocate(nspecies)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_properties), intent(inout) :: self

    deallocate(self%alogsig)
    deallocate(self%sigmag_amode)
    deallocate(self%f1)
    deallocate(self%f2)

    call self%final()

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure real(r8) function voltonumblo(self,m)
    class(modal_aerosol_properties), intent(in) :: self
    integer,intent(in) :: m
    voltonumblo = self%voltonumblo_(m)
  end function voltonumblo

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure real(r8) function voltonumbhi(self,m)
    class(modal_aerosol_properties), intent(in) :: self
    integer,intent(in) :: m
    voltonumbhi = self%voltonumbhi_(m)
  end function voltonumbhi

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
