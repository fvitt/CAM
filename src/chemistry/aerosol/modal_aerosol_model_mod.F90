module modal_aerosol_model_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc
  use physconst,        only: pi
  use ppgrid,           only: pcols
  use cam_abortutils,   only: endrun
  use phys_control,     only: phys_getopts
  use aerosol_model_mod, only: aerosol_model, twothird, sq2
  use modal_cam_aerosol_data_mod, only: modal_cam_aerosol_data

  implicit none
  private
  public :: modal_aerosol_model

  type, extends(aerosol_model) :: modal_aerosol_model
     integer               :: ntot_amode     ! number of aerosol modes
     integer,  allocatable :: nspec_amode(:) ! number of chemical species in each aerosol mode
   contains
     procedure :: model_init => modal_model_init
     procedure :: model_final => modal_model_final
  end type modal_aerosol_model

contains

  subroutine modal_model_init(self)
    class(modal_aerosol_model), intent(inout) :: self

    integer :: m, nspec_max, l,idx

    allocate(modal_cam_aerosol_data::self%aero_data)
    call self%aero_data%initialize()

    self%model_name = 'modal'

    call phys_getopts(prog_modal_aero_out = self%prognostic)

    select type (obj=>self%aero_data)
    type is (modal_cam_aerosol_data)
       self%ntot_amode = obj%ntot_amode
    end select

    self%mtotal = self%ntot_amode

    allocate( &
         self%nspec_amode(self%ntot_amode),  &
         self%nmasses(self%ntot_amode),  &
         self%alogsig(self%ntot_amode),      &
         self%exp45logsig(self%ntot_amode),  &
         self%f1(self%ntot_amode),           &
         self%f2(self%ntot_amode) )

    allocate( self%amcubecoef(self%mtotal) )
    allocate( self%argfactor(self%mtotal) )

    select type (obj=>self%aero_data)
    type is (modal_cam_aerosol_data)
       self%nspec_amode(:) = obj%nspec_amode(:)
    end select

    self%ncnst_tot = 0

    do m = 1, self%ntot_amode
       ! use only if width of size distribution is prescribed

       self%nmasses(m) = self%nspec_amode(m)
       self%ncnst_tot =  self%ncnst_tot + self%nspec_amode(m) + 1

       select type (obj=>self%aero_data)
       type is (modal_cam_aerosol_data)
          self%alogsig(m) = log(obj%sigmag_amode(m))
       end select

       self%exp45logsig(m) = exp(4.5_r8*self%alogsig(m)*self%alogsig(m))
       self%f1(m)          = 0.5_r8*exp(2.5_r8*self%alogsig(m)*self%alogsig(m))
       self%f2(m)          = 1._r8 + 0.25_r8*self%alogsig(m)
       self%amcubecoef(m)=3._r8/(4._r8*pi*self%exp45logsig(m))
       self%argfactor(m)=twothird/(sq2*self%alogsig(m))
    end do

    ! Find max number of species in all the modes
    nspec_max = maxval(self%nspec_amode)
    allocate( self%indexer(self%ntot_amode,0:nspec_max) )

  end subroutine modal_model_init

  subroutine modal_model_final(self)
    class(modal_aerosol_model), intent(inout) :: self

    deallocate( &
         self%nspec_amode,  &
         self%nmasses, &
         self%alogsig,      &
         self%exp45logsig,  &
         self%f1,           &
         self%f2 )
    deallocate( self%amcubecoef )
    deallocate( self%argfactor )
    deallocate( self%indexer )

  end subroutine modal_model_final

end module modal_aerosol_model_mod
