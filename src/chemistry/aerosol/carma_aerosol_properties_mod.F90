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
     procedure :: number_transported
     procedure :: get
     procedure :: amcube
     procedure :: actfracs
     procedure :: num_names
     procedure :: mmr_names
     procedure :: amb_num_name
     procedure :: amb_mmr_name
     procedure :: species_type
     procedure :: icenuc_updates_num
     procedure :: icenuc_updates_mmr
     procedure :: apply_number_limits
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
    integer :: ierr

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    call rad_cnst_get_info( 0, nbins=nbins)

    allocate( nspecies(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( nmasses(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( alogsig(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( f1(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( f2(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    ncnst_tot = 0

    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=nspecies(m))
       ncnst_tot = ncnst_tot + nspecies(m) + 2
       nmasses(m) = nspecies(m) + 1
    end do

    alogsig(:) = log(2._r8)  !!!! ???? IS THIS RIGHT ???? !!!
    f1 = 1._r8
    f2 = 1._r8

    call newobj%initialize(nbins,ncnst_tot,nspecies,nmasses,alogsig,f1,f2,ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
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
  ! returns number of transported aerosol constituents
  !------------------------------------------------------------------------------
  integer function number_transported(self)
    class(carma_aerosol_properties), intent(in) :: self
    ! to be implemented later
    number_transported = -1
  end function number_transported

  !------------------------------------------------------------------------
  ! returns aerosol properties:
  !  density
  !  hygroscopicity
  !------------------------------------------------------------------------
  subroutine get(self, bin_ndx, species_ndx, density,hygro)

    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: species_ndx         ! species index
    real(r8), optional, intent(out) :: density ! density (kg/m3)
    real(r8), optional, intent(out) :: hygro   ! hygroscopicity

    call rad_cnst_get_bin_props_by_idx(0, bin_ndx, species_ndx, density_aer=density, hygro_aer=hygro)

  end subroutine get

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

  !------------------------------------------------------------------------------
  ! returns mass and number activation fractions
  !------------------------------------------------------------------------------
  subroutine actfracs(self, bin_ndx, smc, smax, fn, fm )
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx   ! bin index
    real(r8),intent(in) :: smc       ! critical supersaturation for particles of bin radius
    real(r8),intent(in) :: smax      ! maximum supersaturation for multiple competing aerosols
    real(r8),intent(out) :: fn       ! activation fraction for aerosol number
    real(r8),intent(out) :: fm       ! activation fraction for aerosol mass

    fn = 0._r8
    fm = 0._r8

    if (smc < smax) then
       fn = 1._r8
       fm = 1._r8
    end if

  end subroutine actfracs

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine num_names(self, bin_ndx, name_a, name_c)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=32), intent(out) :: name_a ! constituent name of ambient aerosol number dens
    character(len=32), intent(out) :: name_c ! constituent name of cloud-borne aerosol number dens

    call rad_cnst_get_info_by_bin(0, bin_ndx, num_name=name_a, num_name_cw=name_c)

  end subroutine num_names

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=32), intent(out) :: name_a ! constituent name of ambient aerosol MMR
    character(len=32), intent(out) :: name_c ! constituent name of cloud-borne aerosol MMR

    if (species_ndx>1) then
       call rad_cnst_get_info_by_bin_spec(0, bin_ndx, species_ndx-1, spec_name=name_a, spec_name_cw=name_c)
    else
       call rad_cnst_get_info_by_bin(0, bin_ndx,  mmr_name=name_a, mmr_name_cw=name_c)
    end if

  end subroutine mmr_names

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_num_name(self, bin_ndx, name)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=32), intent(out) :: name   ! constituent name of ambient aerosol number dens

    call rad_cnst_get_info_by_bin(0, bin_ndx, num_name=name)

  end subroutine amb_num_name

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_mmr_name(self, bin_ndx, species_ndx, name)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=32), intent(out) :: name   ! constituent name of ambient aerosol MMR

    if (species_ndx>0) then
       call rad_cnst_get_info_by_bin_spec(0, bin_ndx, species_ndx, spec_name=name)
    else
       call rad_cnst_get_info_by_bin(0, bin_ndx,  mmr_name=name)
    end if

  end subroutine amb_mmr_name

  !------------------------------------------------------------------------
  ! returns species type
  !------------------------------------------------------------------------
  subroutine species_type(self, bin_ndx, species_ndx, spectype)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=32), intent(out) :: spectype ! species type

    call rad_cnst_get_info_by_bin_spec(0, bin_ndx, species_ndx, spec_type=spectype)

  end subroutine species_type

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to given aerosol bin number
  !------------------------------------------------------------------------------
  function icenuc_updates_num(self, bin_ndx) result(res)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    logical :: res

    character(len=32) :: spectype
    integer :: spc_ndx

    res = .false.

    do spc_ndx = 1, self%nspecies(bin_ndx)
       call self%species_type( bin_ndx, spc_ndx, spectype)
       if (trim(spectype)=='dust') res = .true.
       if (trim(spectype)=='sulfate') res = .true.
    end do

  end function icenuc_updates_num

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to a given species within a bin
  !------------------------------------------------------------------------------
  function icenuc_updates_mmr(self, bin_ndx, species_ndx) result(res)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    logical :: res

    character(len=32) :: spectype

    res = .false.

    if (species_ndx==0) then
       res = self%icenuc_updates_num(bin_ndx)
    else
       call self%species_type( bin_ndx, species_ndx, spectype)
       if (trim(spectype)=='dust') res = .true.
       if (trim(spectype)=='sulfate') res = .true.
    end if

  end function icenuc_updates_mmr

  !------------------------------------------------------------------------------
  ! apply max / min to number concentration
  !------------------------------------------------------------------------------
  subroutine apply_number_limits( self, naerosol, vaerosol, istart, istop, m )
    class(carma_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: naerosol(:)  ! number conc (1/m3)
    real(r8), intent(in)    :: vaerosol(:)  ! volume conc (m3/m3)
    integer,  intent(in) :: istart          ! start column index (1 <= istart <= istop <= pcols)
    integer,  intent(in) :: istop           ! stop column index
    integer,  intent(in) :: m               ! mode or bin index

  end subroutine apply_number_limits

end module carma_aerosol_properties_mod
