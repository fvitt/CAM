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
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
     integer, pointer :: indexer_(:,:) => null()
   contains

     procedure :: get_ambient_mmr
     procedure :: get_cldbrne_mmr
     procedure :: get_ambient_num
     procedure :: get_cldbrne_num
     procedure :: indexer
     procedure :: get_states

     final :: destructor

  end type modal_aerosol_state

  interface modal_aerosol_state
     procedure :: constructor
  end interface modal_aerosol_state

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(state,pbuf,props) result(newobj)
    type(physics_state), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)
    class(aerosol_properties), intent(in) :: props

    type(modal_aerosol_state), pointer :: newobj

    integer :: l,m,mm

    allocate(newobj)

    newobj%state => state
    newobj%pbuf => pbuf

    allocate( newobj%indexer_(props%nbins(),0:props%nspec_max()) )

    newobj%indexer_ = -1
    mm = 0

    do m=1,props%nbins()
       do l = 0,props%nspecies(m) ! loop over number + chem constituents

          mm = mm+1
          newobj%indexer_(m,l) = mm

       end do
    end do

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_state), intent(inout) :: self

    nullify(self%state)
    nullify(self%pbuf)
    deallocate(self%indexer_)
    nullify(self%indexer_)

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_ambient_mmr(self,l,m) result(x)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_aer_mmr(0, m, l, 'a', self%state, self%pbuf, x)
  end function get_ambient_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_cldbrne_mmr(self,l,m) result(x)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: l
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_aer_mmr(0, m, l, 'c', self%state, self%pbuf, x)
  end function get_cldbrne_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_ambient_num(self,m) result(x)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_mode_num(0, m, 'a', self%state, self%pbuf, x)
  end function get_ambient_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_cldbrne_num(self,m) result(x)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: m
    real(r8), pointer :: x(:,:)

    call rad_cnst_get_mode_num(0, m, 'c', self%state, self%pbuf, x)
  end function get_cldbrne_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: m,mm,l

    do m = 1, aero_props%nbins()
       mm = self%indexer_(m, 0)
       raer(mm)%fld => self%get_ambient_num(m)
       qqcw(mm)%fld => self%get_cldbrne_num(m)
       do l = 1, aero_props%nspecies(m)
          mm = self%indexer_(m, l)
          raer(mm)%fld => self%get_ambient_mmr(l,m)
          qqcw(mm)%fld => self%get_cldbrne_mmr(l,m)
       end do
    end do

  end subroutine get_states

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  pure function indexer(self, m,l) result(ndx)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: m,l
    integer :: ndx
    ndx = self%indexer_(m,l)
  end function indexer

end module modal_aerosol_state_mod
