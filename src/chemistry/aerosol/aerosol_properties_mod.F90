module aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: aerosol_properties

  type, abstract :: aerosol_properties
     private
     integer :: nbins_ = 0
     integer, allocatable :: nspecies_(:)
   contains
     procedure :: initialize => aero_props_init
     procedure :: nbins
     procedure :: nspecies
     procedure(aero_props_abdraz_f), deferred :: abdraz_f1
     procedure(aero_props_abdraz_f), deferred :: abdraz_f2
     procedure(aero_props_get), deferred :: get
     procedure :: final=>aero_props_final
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

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_props_init(self, n, nspec )
    class(aerosol_properties), intent(inout) :: self
    integer :: n
    integer :: nspec(n)

    self%nbins_ = n
    allocate(self%nspecies_(n))
    self%nspecies_(:) = nspec(:)
  end subroutine aero_props_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine aero_props_final(self)
    class(aerosol_properties), intent(inout) :: self

    if (allocated(self%nspecies_)) then
       deallocate(self%nspecies_)
    end if

    self%nbins_ = 0
  end subroutine aero_props_final

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure integer function nspecies(self,m)
    class(aerosol_properties), intent(in) :: self
    integer, intent(in) :: m

    nspecies = self%nspecies_(m)
  end function nspecies

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure integer function nbins(self)
    class(aerosol_properties), intent(in) :: self

    nbins = self%nbins_
  end function nbins

end module aerosol_properties_mod
