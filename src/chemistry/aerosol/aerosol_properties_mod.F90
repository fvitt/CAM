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
     procedure :: maxsat     ! *** Does this belong with this class ??
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

  !------------------------------------------------------------------------------
  ! *** Does maxsat belong with this class ??
  !------------------------------------------------------------------------------
  function maxsat(self, zeta,eta,smc) result(smax)

    !-------------------------------------------------------------------------
    ! Calculates maximum supersaturation for multiple competing aerosol modes.
    !
    ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.
    !-------------------------------------------------------------------------

    class(aerosol_properties), intent(in) :: self
    real(r8), intent(in)  :: zeta(self%nbins_)
    real(r8), intent(in)  :: eta(self%nbins_)
    real(r8), intent(in)  :: smc(self%nbins_) ! critical supersaturation for number mode radius

    real(r8) :: smax ! maximum supersaturation

    integer  :: m
    integer  :: nbins
    real(r8) :: sum, g1, g2, g1sqrt, g2sqrt

    smax=0.0_r8
    nbins = self%nbins_

    check_loop: do m=1,nbins
       if((zeta(m) > 1.e5_r8*eta(m)) .or. (smc(m)*smc(m) > 1.e5_r8*eta(m))) then
          ! weak forcing -- essentially none activated
          smax=1.e-20_r8
       else
          ! significant activation of this mode -- calc activation all modes
          exit check_loop
       endif
       ! No significant activation in any mode.  Do nothing.
       if (m == nbins) return
    enddo check_loop

    sum=0.0_r8

    do m=1,nbins
       if(eta(m) > 1.e-20_r8)then
          g1=zeta(m)/eta(m)
          g1sqrt=sqrt(g1)
          g1=g1sqrt*g1
          g2=smc(m)/sqrt(eta(m)+3._r8*zeta(m))
          g2sqrt=sqrt(g2)
          g2=g2sqrt*g2
          sum=sum+(self%abdraz_f1(m)*g1+self%abdraz_f2(m)*g2)/(smc(m)*smc(m))
       else
          sum=1.e20_r8
       endif
    enddo

    smax=1._r8/sqrt(sum)

  end function maxsat

end module aerosol_properties_mod
