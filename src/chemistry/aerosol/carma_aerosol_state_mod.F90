module carma_aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_state_mod, only: aerosol_state, ptr2d_t

  use rad_constituents, only: rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_num, rad_cnst_get_bin_mmr
  use physics_buffer, only: physics_buffer_desc
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: carma_aerosol_state

  type, extends(aerosol_state) :: carma_aerosol_state
     private
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
   contains

     procedure :: get_ambient_mmr
     procedure :: get_cldbrne_mmr
     procedure :: get_ambient_num
     procedure :: get_cldbrne_num
     procedure :: get_states

     final :: destructor

  end type carma_aerosol_state

  interface carma_aerosol_state
     procedure :: constructor
  end interface carma_aerosol_state

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(state,pbuf,props) result(newobj)
    type(physics_state), target, optional :: state
    type(physics_buffer_desc), pointer, optional :: pbuf(:)
    class(aerosol_properties), intent(in) :: props

    type(carma_aerosol_state), pointer :: newobj

    allocate(newobj)

    newobj%state => state
    newobj%pbuf => pbuf

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(carma_aerosol_state), intent(inout) :: self

    nullify(self%state)
    nullify(self%pbuf)

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr(self, species_ndx, bin_ndx, mmr)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios

    if (species_ndx<1) then
       call rad_cnst_get_bin_mmr(0, bin_ndx, 'a', self%state, self%pbuf, mmr)
    else
       call rad_cnst_get_bin_mmr_by_idx(0, bin_ndx, species_ndx, 'a', self%state, self%pbuf, mmr)
    end if

  end subroutine get_ambient_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, mmr)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios

    if (species_ndx<1) then
       call rad_cnst_get_bin_mmr(0, bin_ndx, 'c', self%state, self%pbuf, mmr)
    else
       call rad_cnst_get_bin_mmr_by_idx(0, bin_ndx, species_ndx, 'c', self%state, self%pbuf, mmr)
    end if

  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, num)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    real(r8), pointer :: num(:,:)

    call rad_cnst_get_bin_num(0, bin_ndx, 'a', self%state, self%pbuf, num)
  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, num)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    real(r8), pointer :: num(:,:)

    call rad_cnst_get_bin_num(0, bin_ndx, 'c', self%state, self%pbuf, num)
  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(carma_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: m, mm, l

    do m = 1, aero_props%nbins()
       mm = aero_props%indexer(m, 0)
       call self%get_ambient_num(m, raer(mm)%fld)
       call self%get_cldbrne_num(m, qqcw(mm)%fld)
       mm = aero_props%indexer(m, 1)
       call self%get_ambient_mmr(0,m, raer(mm)%fld)
       call self%get_cldbrne_mmr(0,m, qqcw(mm)%fld)
       do l = 1, aero_props%nspecies(m)
          mm = aero_props%indexer(m, l+1)
          call self%get_ambient_mmr(l,m, raer(mm)%fld)
          call self%get_cldbrne_mmr(l,m, qqcw(mm)%fld)
       end do
    end do

  end subroutine get_states

end module carma_aerosol_state_mod
