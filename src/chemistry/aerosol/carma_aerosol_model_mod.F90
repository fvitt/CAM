module carma_aerosol_model_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: pi => shr_const_pi

  use aerosol_model_mod, only: aerosol_model, twothird, sq2
  use carma_cam_aerosol_data_mod, only: carma_cam_aerosol_data

  implicit none
  private
  public :: carma_aerosol_model

  type, extends(aerosol_model) :: carma_aerosol_model
   contains
     procedure :: model_init => carma_model_init
     procedure :: model_final => carma_model_final
     procedure :: err_funct => carma_err_funct ! override base routine
  end type carma_aerosol_model

contains

  subroutine carma_model_init(self)
    class(carma_aerosol_model), intent(inout) :: self

    integer :: l, m, nspec_max
    integer :: idxtmp

    allocate(carma_cam_aerosol_data::self%aero_data)
    call self%aero_data%initialize()

    self%model_name = 'carma'

    allocate( self%amcubecoef(self%aero_data%mtotal) )
    allocate( self%exp45logsig(self%aero_data%mtotal) )
    allocate( self%argfactor(self%aero_data%mtotal) )
    allocate( self%alogsig(self%aero_data%mtotal) )
    allocate( self%f1(self%aero_data%mtotal) )
    allocate( self%f2(self%aero_data%mtotal) )


    do m = 1, self%aero_data%mtotal
       self%amcubecoef(m)=3._r8/(4._r8*pi)
       self%argfactor(m)=twothird/(sq2*log(2._r8))
       self%exp45logsig(m)=1._r8
       self%alogsig(m)=1._r8
       self%f1(m)=1._r8
       self%f2(m)=1._r8
    end do

  end subroutine carma_model_init

  subroutine carma_model_final(self)
    class(carma_aerosol_model), intent(inout) :: self

    deallocate( self%amcubecoef )
    deallocate( self%argfactor )
    deallocate( self%alogsig )
    deallocate( self%f1 )
    deallocate( self%f2 )
    deallocate( self%indexer )

  end subroutine carma_model_final

  ! override the error function
  function carma_err_funct(self,x) result(err)
    class(carma_aerosol_model), intent(in) :: self
    real(r8), intent(in) :: x
    real(r8) :: err
    err = -1._r8
  end function carma_err_funct

end module carma_aerosol_model_mod
