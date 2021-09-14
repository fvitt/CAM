module carma_aerosol_model_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc
  use physconst, only: pi
  use aerosol_model_mod, only: aerosol_model, twothird, sq2
  use carma_cam_aerosol_data_mod, only: carma_cam_aerosol_data

  implicit none
  private
  public :: carma_aerosol_model

  type, extends(aerosol_model) :: carma_aerosol_model
     integer, allocatable :: nspec(:)
   contains
     procedure :: model_init => carma_model_init
     procedure :: model_final => carma_model_final
     procedure :: err_funct => carma_err_funct ! override base routine
  end type carma_aerosol_model

contains

  subroutine carma_model_init(self)
    class(carma_aerosol_model), intent(inout) :: self

    integer :: l, m, nbins, nspec_max
    integer :: idxtmp

    allocate(carma_cam_aerosol_data::self%aero_data)
    call self%aero_data%initialize()

    self%model_name = 'carma'

    ! get info about the modal aerosols
    ! get nbins

    select type (obj=>self%aero_data)
    type is (carma_cam_aerosol_data)
       nbins = obj%nbins
    end select

    self%mtotal = nbins

    allocate( self%nspec(nbins) )
    allocate( self%nmasses(nbins) )
    allocate( self%amcubecoef(nbins) )
    allocate( self%exp45logsig(nbins) )
    allocate( self%argfactor(nbins) )
    allocate( self%alogsig(nbins) )
    allocate( self%f1(nbins) )
    allocate( self%f2(nbins) )

    select type (obj=>self%aero_data)
    type is (carma_cam_aerosol_data)
       self%nspec(:) = obj%nspec(:)
    end select

    self%ncnst_tot = 0

    do m = 1, nbins
       self%nmasses(m) = self%nspec(m) + 1
       self%ncnst_tot =  self%ncnst_tot + self%nspec(m) + 2
       self%amcubecoef(m)=3._r8/(4._r8*pi)
       self%argfactor(m)=twothird/(sq2*log(2._r8))
       self%exp45logsig(m)=1._r8
       self%alogsig(m)=1._r8
       self%f1(m)=1._r8
       self%f2(m)=1._r8
    end do

    ! add plus one to include number, total mmr and nspec
    nspec_max = maxval(self%nspec) + 1
    allocate( self%indexer(nbins,0:nspec_max) )

  end subroutine carma_model_init

  subroutine carma_model_final(self)
    class(carma_aerosol_model), intent(inout) :: self

    deallocate( self%nspec )
    deallocate( self%nmasses )
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
