module carma_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: carma_aerosol_properties

  type, extends(aerosol_properties) :: carma_aerosol_properties
     private
   contains
     procedure :: abdraz_f1
     procedure :: abdraz_f2
     procedure :: get
     final :: destructor
  end type carma_aerosol_properties

  interface carma_aerosol_properties
     procedure :: constructor
  end interface carma_aerosol_properties

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor() result(newobj)
    use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_info_by_bin
    use physconst,        only: pi

    type(carma_aerosol_properties), pointer :: newobj

    integer :: m, nbins
    integer,allocatable :: nspecies(:)
    real(r8),allocatable :: amcubecoefs(:)

    allocate(newobj)

    call rad_cnst_get_info( 0, nbins=nbins)
    allocate( nspecies(nbins) )
    allocate( amcubecoefs(nbins) )

    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=nspecies(m))
    end do

    amcubecoefs(:)=3._r8/(4._r8*pi)

    call newobj%initialize(nbins,nspecies,amcubecoefs)
    deallocate(nspecies)
    deallocate(amcubecoefs)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(carma_aerosol_properties), intent(inout) :: self

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get(self, m,l, density,hygro)
    use rad_constituents, only:  rad_cnst_get_bin_props_by_idx

    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m,l
    real(r8), optional, intent(out) :: density
    real(r8), optional, intent(out) :: hygro

    call rad_cnst_get_bin_props_by_idx(0, m, l, density_aer=density, hygro_aer=hygro)

  end subroutine get

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function abdraz_f1(self,m) result(f)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8) :: f

    f = 1._r8
  end function abdraz_f1

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function abdraz_f2(self,m) result(f)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8) :: f

    f = 1._r8
  end function abdraz_f2

end module carma_aerosol_properties_mod
