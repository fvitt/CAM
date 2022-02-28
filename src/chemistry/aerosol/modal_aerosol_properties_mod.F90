module modal_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi
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
     real(r8), allocatable :: exp45logsig_(:)
     real(r8), allocatable :: voltonumblo_(:)
     real(r8), allocatable :: voltonumbhi_(:)
   contains
     procedure :: abdraz_f1
     procedure :: abdraz_f2
     procedure :: get
     procedure :: exp45logsig
     procedure :: voltonumblo
     procedure :: voltonumbhi
     procedure :: amcube
     procedure :: actfracs
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

    type(modal_aerosol_properties), pointer :: newobj

    integer :: m, nmodes
    real(r8) :: dgnumlo
    real(r8) :: dgnumhi
    integer,allocatable :: nspecies(:)
    real(r8),allocatable :: amcubecoefs(:)

    allocate(newobj)

    call rad_cnst_get_info(0, nmodes=nmodes)

    allocate(nspecies(nmodes))
    allocate(amcubecoefs(nmodes))

    allocate(newobj%alogsig(nmodes))
    allocate(newobj%sigmag_amode(nmodes))
    allocate(newobj%f1(nmodes))
    allocate(newobj%f2(nmodes))
    allocate(newobj%exp45logsig_(nmodes))
    allocate(newobj%voltonumblo_(nmodes))
    allocate(newobj%voltonumbhi_(nmodes))

    do m = 1, nmodes
       call rad_cnst_get_info(0, m, nspec=nspecies(m))

       call rad_cnst_get_mode_props(0, m, sigmag=newobj%sigmag_amode(m), &
                                    dgnumhi=dgnumhi, dgnumlo=dgnumlo )

       newobj%alogsig(m) = log(newobj%sigmag_amode(m))

       newobj%exp45logsig_(m) = exp(4.5_r8*newobj%alogsig(m)*newobj%alogsig(m))

       amcubecoefs(m)=3._r8/(4._r8*pi*newobj%exp45logsig_(m))

       newobj%f1(m) = 0.5_r8*exp(2.5_r8*newobj%alogsig(m)*newobj%alogsig(m))
       newobj%f2(m) = 1._r8 + 0.25_r8*newobj%alogsig(m)

       newobj%voltonumblo_(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumlo**3._r8)*exp(4.5_r8*newobj%alogsig(m)**2._r8) )
       newobj%voltonumbhi_(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumhi**3._r8)*exp(4.5_r8*newobj%alogsig(m)**2._r8) )

    end do

    call newobj%initialize(nmodes,nspecies,amcubecoefs)
    deallocate(nspecies)
    deallocate(amcubecoefs)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_properties), intent(inout) :: self

    deallocate(self%alogsig)
    deallocate(self%sigmag_amode)
    deallocate(self%f1)
    deallocate(self%f2)
    deallocate(self%exp45logsig_)
    deallocate(self%voltonumblo_)
    deallocate(self%voltonumbhi_)

    call self%final()

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  pure real(r8) function exp45logsig(self,m)
    class(modal_aerosol_properties), intent(in) :: self
    integer,intent(in) :: m
    exp45logsig = self%exp45logsig_(m)
  end function exp45logsig

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

  !------------------------------------------------------------------------------
  ! amcube is overridden to keep MAM b4b
  !------------------------------------------------------------------------------
  pure real(r8) function amcube(self, m, volconc, numconc)

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8), intent(in) :: volconc
    real(r8), intent(in) :: numconc

    amcube = (3._r8*volconc/(4._r8*pi*self%exp45logsig_(m)*numconc))

  end function amcube

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine actfracs(self, m, smc, smax, fn, fm )
    use shr_spfn_mod, only: erf => shr_spfn_erf
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: m
    real(r8),intent(in) :: smc
    real(r8),intent(in) :: smax
    real(r8),intent(out) :: fn, fm

    real(r8) :: x,y
    real(r8), parameter :: twothird = 2._r8/3._r8
    real(r8), parameter :: sq2      = sqrt(2._r8)

    x=twothird*(log(smc)-log(smax))/(sq2*self%alogsig(m))
    y=x-1.5_r8*sq2*self%alogsig(m)

    fn = 0.5_r8*(1._r8-erf(x))
    fm = 0.5_r8*(1._r8-erf(y))

  end subroutine actfracs

end module modal_aerosol_properties_mod
