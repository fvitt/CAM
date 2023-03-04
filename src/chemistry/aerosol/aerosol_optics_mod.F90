module aerosol_optics_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private
  public :: aerosol_optics

  type, abstract :: aerosol_optics

   contains

     procedure(aeropts_sw_props),deferred :: sw_props
     procedure(aeropts_lw_props),deferred :: lw_props

  end type aerosol_optics

  abstract interface

     subroutine aeropts_sw_props(self, ncol, ilev, iwav, pext, palb, pasm)
       import :: aerosol_optics, r8

       class(aerosol_optics), intent(in) :: self
       integer, intent(in) :: ncol
       integer, intent(in) :: ilev
       integer, intent(in) :: iwav
       real(r8),intent(out) :: pext(ncol)
       real(r8),intent(out) :: palb(ncol)
       real(r8),intent(out) :: pasm(ncol)

     end subroutine aeropts_sw_props

     subroutine aeropts_lw_props(self, ncol, ilev, iwav, pabs)
       import :: aerosol_optics, r8

       class(aerosol_optics), intent(in) :: self
       integer, intent(in) :: ncol
       integer, intent(in) :: ilev
       integer, intent(in) :: iwav
       real(r8),intent(out) :: pabs(ncol)

     end subroutine aeropts_lw_props

     subroutine aeropts_destroy(self)
       import :: aerosol_optics
       class(aerosol_optics), intent(inout) :: self
     end subroutine aeropts_destroy

  end interface

contains

end module aerosol_optics_mod
