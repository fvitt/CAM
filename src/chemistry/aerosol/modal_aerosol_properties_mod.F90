module modal_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi
  use aerosol_properties_mod, only: aerosol_properties
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_props, rad_cnst_get_aer_props

  implicit none

  private

  public :: modal_aerosol_properties

  type, extends(aerosol_properties) :: modal_aerosol_properties
     private
     real(r8), allocatable :: exp45logsig_(:)
     real(r8), allocatable :: voltonumblo_(:)
     real(r8), allocatable :: voltonumbhi_(:)
   contains
     procedure :: get
     procedure :: voltonumblo
     procedure :: voltonumbhi
     procedure :: amcube
     procedure :: actfracs
     procedure :: num_names
     procedure :: mmr_names
     procedure :: amb_num_name
     procedure :: amb_mmr_name
     procedure :: species_type
     procedure :: icenuc_num
     procedure :: icenuc_mmr
     final :: destructor
  end type modal_aerosol_properties

  interface modal_aerosol_properties
     procedure :: constructor
  end interface modal_aerosol_properties

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor() result(newobj)

    type(modal_aerosol_properties), pointer :: newobj

    integer :: m, nmodes, ncnst_tot
    real(r8) :: dgnumlo
    real(r8) :: dgnumhi
    integer,allocatable :: nspecies(:)
    integer,allocatable :: nmasses(:)
    real(r8),allocatable :: sigmag(:)
    real(r8),allocatable :: alogsig(:)
    real(r8),allocatable :: f1(:)
    real(r8),allocatable :: f2(:)

    allocate(newobj)

    call rad_cnst_get_info(0, nmodes=nmodes)

    allocate(nspecies(nmodes))
    allocate(nmasses(nmodes))
    allocate(alogsig(nmodes))
    allocate( f1(nmodes) )
    allocate( f2(nmodes) )

    allocate(sigmag(nmodes))
    allocate(newobj%exp45logsig_(nmodes))
    allocate(newobj%voltonumblo_(nmodes))
    allocate(newobj%voltonumbhi_(nmodes))

    ncnst_tot = 0

    do m = 1, nmodes
       call rad_cnst_get_info(0, m, nspec=nspecies(m))

       ncnst_tot =  ncnst_tot + nspecies(m) + 1
       nmasses(m) = nspecies(m)

       call rad_cnst_get_mode_props(0, m, sigmag=sigmag(m), &
                                    dgnumhi=dgnumhi, dgnumlo=dgnumlo )

       alogsig(m) = log(sigmag(m))

       newobj%exp45logsig_(m) = exp(4.5_r8*alogsig(m)*alogsig(m))

       f1(m) = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
       f2(m) = 1._r8 + 0.25_r8*alogsig(m)

       newobj%voltonumblo_(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumlo**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )
       newobj%voltonumbhi_(m) = 1._r8 / ( (pi/6._r8)*                          &
            (dgnumhi**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )

    end do

    call newobj%initialize(nmodes,ncnst_tot,nspecies,nmasses,alogsig,f1,f2)
    deallocate(nspecies)
    deallocate(nmasses)
    deallocate(alogsig)
    deallocate(sigmag)
    deallocate(f1)
    deallocate(f2)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_properties), intent(inout) :: self

    deallocate(self%exp45logsig_)
    deallocate(self%voltonumblo_)
    deallocate(self%voltonumbhi_)

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
  subroutine get(self, bin_ndx, species_ndx, density,hygro)

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: species_ndx         ! species index
    real(r8), optional, intent(out) :: density
    real(r8), optional, intent(out) :: hygro

    call rad_cnst_get_aer_props(0, bin_ndx, species_ndx, density_aer=density, hygro_aer=hygro)

  end subroutine get

  !------------------------------------------------------------------------------
  ! returns radius^3 (m3) of a given bin number
  !------------------------------------------------------------------------------
  pure elemental real(r8) function amcube(self, bin_ndx, volconc, numconc)

    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    real(r8), intent(in) :: volconc ! volume conc (m3/m3)
    real(r8), intent(in) :: numconc ! number conc (1/m3)

    amcube = (3._r8*volconc/(4._r8*pi*self%exp45logsig_(bin_ndx)*numconc))

  end function amcube

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine actfracs(self, bin_ndx, smc, smax, fn, fm )
    use shr_spfn_mod, only: erf => shr_spfn_erf
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx
    real(r8),intent(in) :: smc
    real(r8),intent(in) :: smax
    real(r8),intent(out) :: fn, fm

    real(r8) :: x,y
    real(r8), parameter :: twothird = 2._r8/3._r8
    real(r8), parameter :: sq2      = sqrt(2._r8)

    x=twothird*(log(smc)-log(smax))/(sq2*self%alogsig(bin_ndx))
    y=x-1.5_r8*sq2*self%alogsig(bin_ndx)

    fn = 0.5_r8*(1._r8-erf(x))
    fm = 0.5_r8*(1._r8-erf(y))

  end subroutine actfracs

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine num_names(self, bin_ndx, name_a, name_c)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=32), intent(out) :: name_a, name_c

    call rad_cnst_get_info(0,bin_ndx, num_name=name_a, num_name_cw=name_c)
  end subroutine num_names

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=32), intent(out) :: name_a, name_c

    call rad_cnst_get_info(0, bin_ndx, species_ndx, spec_name=name_a, spec_name_cw=name_c)
  end subroutine mmr_names

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine amb_num_name(self, bin_ndx, name)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=32), intent(out) :: name   ! constituent name of ambient aerosol number dens

    call rad_cnst_get_info(0,bin_ndx, num_name=name)

  end subroutine amb_num_name
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine amb_mmr_name(self, bin_ndx, species_ndx, name)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=32), intent(out) :: name   ! constituent name of ambient aerosol MMR

    call rad_cnst_get_info(0, bin_ndx, species_ndx, spec_name=name)

  end subroutine amb_mmr_name

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine species_type(self, bin_ndx, species_ndx, spectype)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=32), intent(out) :: spectype ! species type

    call rad_cnst_get_info(0, bin_ndx, species_ndx, spec_type=spectype)

  end subroutine species_type

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function icenuc_num(self, bin_ndx) result(res)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    logical :: res

    character(len=32) :: spectype
    character(len=32) :: modetype
    integer :: spc_ndx

    res = .false.

    call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)
    if (.not.(modetype=='coarse' .or. modetype=='coarse_dust')) then
       return
    end if

    do spc_ndx = 1, self%nspecies(bin_ndx)
       call self%species_type( bin_ndx, spc_ndx, spectype)
       if (spectype=='dust') res = .true.
    end do

  end function icenuc_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function icenuc_mmr(self, bin_ndx, species_ndx) result(res)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    logical :: res

    character(len=32) :: spectype
    character(len=32) :: modetype

    res = .false.

    if (species_ndx>0) then

       call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)
       if (.not.(modetype=='coarse' .or. modetype=='coarse_dust')) then
          return
       end if

       call self%species_type( bin_ndx, species_ndx, spectype)
       if (spectype=='dust') res = .true.
    end if

  end function icenuc_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function icenuc_bin_wght(self, bin_ndx, species_type) result(wght)
    class(modal_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    character(len=*), intent(in) :: species_type  ! species type

    real(r8) :: wght
    character(len=32) :: modetype

    call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)

    wght = 0._r8

    if (species_type=='dust') then
       if (modetype=='coarse' .or. modetype=='coarse_dust') then
          wght = 1._r8
       end if
    end if

  end function icenuc_bin_wght

end module modal_aerosol_properties_mod
