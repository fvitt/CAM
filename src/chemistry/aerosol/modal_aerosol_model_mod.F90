module modal_aerosol_model_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: pi => shr_const_pi

  use aerosol_model_mod, only: aerosol_model, twothird, sq2

  implicit none
  private
  public :: modal_aerosol_model

  type, extends(aerosol_model) :: modal_aerosol_model
   contains
     procedure :: model_init => modal_model_init
     procedure :: model_final => modal_model_final
  end type modal_aerosol_model

contains

  subroutine modal_model_init(self)
    class(modal_aerosol_model), intent(inout) :: self

    integer :: m, nspec_max, l,idx

    self%model_name = 'modal'

    allocate( &
         self%alogsig(self%aero_data%mtotal),      &
         self%exp45logsig(self%aero_data%mtotal),  &
         self%f1(self%aero_data%mtotal),           &
         self%f2(self%aero_data%mtotal) )

    allocate( self%amcubecoef(self%aero_data%mtotal) )
    allocate( self%argfactor(self%aero_data%mtotal) )

    do m = 1, self%aero_data%mtotal
       ! use only if width of size distribution is prescribed

       self%alogsig(m) = log(self%aero_data%sigmag_amode(m))

       self%exp45logsig(m) = exp(4.5_r8*self%alogsig(m)*self%alogsig(m))
       self%f1(m)          = 0.5_r8*exp(2.5_r8*self%alogsig(m)*self%alogsig(m))
       self%f2(m)          = 1._r8 + 0.25_r8*self%alogsig(m)
       self%amcubecoef(m)  = 3._r8/(4._r8*pi*self%exp45logsig(m))
       self%argfactor(m)   = twothird/(sq2*self%alogsig(m))
    end do

  end subroutine modal_model_init

  subroutine modal_model_final(self)
    class(modal_aerosol_model), intent(inout) :: self

    deallocate( &
         self%alogsig,      &
         self%exp45logsig,  &
         self%f1,           &
         self%f2 )
    deallocate( self%amcubecoef )
    deallocate( self%argfactor )

  end subroutine modal_model_final

end module modal_aerosol_model_mod
