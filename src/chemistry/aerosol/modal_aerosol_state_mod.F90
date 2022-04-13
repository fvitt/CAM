module modal_aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_state_mod, only: aerosol_state, ptr2d_t
  use rad_constituents, only: rad_cnst_get_aer_mmr, rad_cnst_get_mode_num
  use physics_buffer, only: physics_buffer_desc
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: modal_aerosol_state

  type, extends(aerosol_state) :: modal_aerosol_state
     private
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
   contains

     procedure :: get_ambient_mmr
     procedure :: get_cldbrne_mmr
     procedure :: get_ambient_num
     procedure :: get_cldbrne_num
     procedure :: get_states

     final :: destructor

  end type modal_aerosol_state

  interface modal_aerosol_state
     procedure :: constructor
  end interface modal_aerosol_state

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(pbuf,props) result(newobj)
    type(physics_buffer_desc), pointer :: pbuf(:)
    class(aerosol_properties), intent(in) :: props

    type(modal_aerosol_state), pointer :: newobj

    allocate(newobj)

    newobj%pbuf => pbuf

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_state), intent(inout) :: self

    nullify(self%pbuf)

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr(self, species_ndx, bin_ndx, cnst_array, mmr)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), target, intent(in) :: cnst_array(:,:,:) ! Constituent array
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios

    call rad_cnst_get_aer_mmr(0, bin_ndx, species_ndx, 'a', cnst_array, self%pbuf, mmr)
  end subroutine get_ambient_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, cnst_array, mmr)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), target, intent(in) :: cnst_array(:,:,:) ! Constituent array
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios

    call rad_cnst_get_aer_mmr(0, bin_ndx, species_ndx, 'c', cnst_array, self%pbuf, mmr)
  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, cnst_array, num)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    real(r8), target, intent(in) :: cnst_array(:,:,:) ! Constituent array
    real(r8), pointer :: num(:,:)

    call rad_cnst_get_mode_num(0, bin_ndx, 'a', cnst_array, self%pbuf, num)
  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, cnst_array, num)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    real(r8), target, intent(in) :: cnst_array(:,:,:) ! Constituent array
    real(r8), pointer :: num(:,:)

    call rad_cnst_get_mode_num(0, bin_ndx, 'c', cnst_array, self%pbuf, num)
  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, cnst_array, raer, qqcw )
    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    real(r8), target, intent(in) :: cnst_array(:,:,:) ! Constituent array
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: m,mm,l

    do m = 1, aero_props%nbins()
       mm = aero_props%indexer(m, 0)
       call self%get_ambient_num(m, cnst_array, raer(mm)%fld)
       call self%get_cldbrne_num(m, cnst_array, qqcw(mm)%fld)
       do l = 1, aero_props%nspecies(m)
          mm = aero_props%indexer(m, l)
          call self%get_ambient_mmr(l,m, cnst_array, raer(mm)%fld)
          call self%get_cldbrne_mmr(l,m, cnst_array, qqcw(mm)%fld)
       end do
    end do

  end subroutine get_states

end module modal_aerosol_state_mod
