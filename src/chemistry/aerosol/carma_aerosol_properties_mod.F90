module carma_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi
  use aerosol_properties_mod, only: aerosol_properties
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_bin_props_by_idx, &
                              rad_cnst_get_info_by_bin, rad_cnst_get_info_by_bin_spec

  implicit none

  private

  public :: carma_aerosol_properties

  type, extends(aerosol_properties) :: carma_aerosol_properties
     private
   contains
     procedure :: get
     procedure :: amcube
     procedure :: actfracs
     procedure :: num_names
     procedure :: mmr_names
     final :: destructor
  end type carma_aerosol_properties

  interface carma_aerosol_properties
     procedure :: constructor
  end interface carma_aerosol_properties

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor() result(newobj)

    type(carma_aerosol_properties), pointer :: newobj

    integer :: m, nbins, ncnst_tot
    integer,allocatable :: nspecies(:)
    integer,allocatable :: nmasses(:)
    real(r8),allocatable :: alogsig(:)
    real(r8),allocatable :: f1(:)
    real(r8),allocatable :: f2(:)

    allocate(newobj)

    call rad_cnst_get_info( 0, nbins=nbins)
    allocate( nspecies(nbins) )
    allocate( nmasses(nbins) )
    allocate( alogsig(nbins) )
    allocate( f1(nbins) )
    allocate( f2(nbins) )

    ncnst_tot = 0

    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=nspecies(m))
       ncnst_tot = ncnst_tot + nspecies(m) + 2
       nmasses(m) = nspecies(m) + 1
    end do

    alogsig(:) = log(2._r8)  !!!! ???? IS THIS RIGHT ???? !!!
    f1 = 1._r8
    f2 = 1._r8

    call newobj%initialize(nbins,ncnst_tot,nspecies,nmasses,alogsig,f1,f2)
    deallocate(nspecies)
    deallocate(nmasses)
    deallocate(alogsig)
    deallocate(f1)
    deallocate(f2)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(carma_aerosol_properties), intent(inout) :: self

    call self%final()

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get(self, bin_ndx, species_ndx, density,hygro)

    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: species_ndx         ! species index
    real(r8), optional, intent(out) :: density
    real(r8), optional, intent(out) :: hygro

    call rad_cnst_get_bin_props_by_idx(0, bin_ndx, species_ndx, density_aer=density, hygro_aer=hygro)

  end subroutine get

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine actfracs(self, bin_ndx, smc, smax, fn, fm )
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx
    real(r8),intent(in) :: smc
    real(r8),intent(in) :: smax
    real(r8),intent(out) :: fn, fm

    fn = 0._r8
    fm = 0._r8

    if (smc < smax) then
       fn = 1._r8
       fm = 1._r8
    end if

  end subroutine actfracs

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine num_names(self, bin_ndx, name_a, name_c)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=32), intent(out) :: name_a, name_c

    call rad_cnst_get_info_by_bin(0, bin_ndx, num_name=name_a, num_name_cw=name_c)

  end subroutine num_names

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  subroutine mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=32), intent(out) :: name_a, name_c

    if (species_ndx>1) then
       call rad_cnst_get_info_by_bin_spec(0, bin_ndx, species_ndx-1, spec_name=name_a, spec_name_cw=name_c)
    else
       call rad_cnst_get_info_by_bin(0, bin_ndx,  mmr_name=name_a, mmr_name_cw=name_c)
    end if

  end subroutine mmr_names

  !------------------------------------------------------------------------------
  ! returns radius^3 (m3) of a given bin number
  !------------------------------------------------------------------------------
  pure elemental real(r8) function amcube(self, bin_ndx, volconc, numconc)

    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    real(r8), intent(in) :: volconc ! volume conc (m3/m3)
    real(r8), intent(in) :: numconc ! number conc (1/m3)

    amcube = 3._r8/(4._r8*pi)*volconc/numconc

  end function amcube

end module carma_aerosol_properties_mod
