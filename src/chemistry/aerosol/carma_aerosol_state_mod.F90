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
     procedure :: get_ambient_bin_mmr
     procedure :: get_cldbrne_bin_mmr
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
  function get_ambient_mmr(self,l,m) result(x)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'a', self%state, self%pbuf, x)
  end function get_ambient_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_cldbrne_mmr(self,l,m) result(x)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'c', self%state, self%pbuf, x)
  end function get_cldbrne_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_ambient_num(self,m) result(x)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_bin_num(0, m, 'a', self%state, self%pbuf, x)
  end function get_ambient_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_cldbrne_num(self,m) result(x)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_bin_num(0, m, 'c', self%state, self%pbuf, x)
  end function get_cldbrne_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_ambient_bin_mmr(self,m) result(x)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_bin_mmr(0, m, 'a', self%state, self%pbuf, x)
  end function get_ambient_bin_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_cldbrne_bin_mmr(self,m) result(x)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_bin_mmr(0, m, 'c', self%state, self%pbuf, x)
  end function get_cldbrne_bin_mmr

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
       raer(mm)%fld => self%get_ambient_num(m)
       qqcw(mm)%fld => self%get_cldbrne_num(m)
       mm = aero_props%indexer(m, 1)
       raer(mm)%fld => self%get_ambient_bin_mmr(m)
       qqcw(mm)%fld => self%get_cldbrne_bin_mmr(m)
       do l = 1, aero_props%nspecies(m)
          mm = aero_props%indexer(m, l+1)
          raer(mm)%fld => self%get_ambient_mmr(l,m)
          qqcw(mm)%fld => self%get_cldbrne_mmr(l,m)
       end do
    end do

  end subroutine get_states

end module carma_aerosol_state_mod
