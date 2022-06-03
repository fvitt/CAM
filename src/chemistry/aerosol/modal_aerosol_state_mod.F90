module modal_aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_state_mod, only: aerosol_state, ptr2d_t
  use rad_constituents, only: rad_cnst_get_aer_mmr, rad_cnst_get_mode_num, rad_cnst_get_info
  use rad_constituents, only: rad_cnst_get_mode_props
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: modal_aerosol_state

  type, extends(aerosol_state) :: modal_aerosol_state
     private
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
   contains

     procedure :: get_ambient_mmr
     procedure :: get_cldbrne_mmr
     procedure :: get_ambient_num
     procedure :: get_cldbrne_num
     procedure :: get_states
     procedure :: icenuc_size_wght
     procedure :: icenuc_type_wght

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

    allocate(newobj)

    newobj%state => state
    newobj%pbuf => pbuf

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(modal_aerosol_state), intent(inout) :: self

    nullify(self%state)
    nullify(self%pbuf)

  end subroutine destructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr(self, species_ndx, bin_ndx, mmr)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios

    call rad_cnst_get_aer_mmr(0, bin_ndx, species_ndx, 'a', self%state, self%pbuf, mmr)
  end subroutine get_ambient_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, mmr)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios

    call rad_cnst_get_aer_mmr(0, bin_ndx, species_ndx, 'c', self%state, self%pbuf, mmr)
  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, num)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer   :: num(:,:)    ! number densities

    call rad_cnst_get_mode_num(0, bin_ndx, 'a', self%state, self%pbuf, num)
  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, num)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    real(r8), pointer :: num(:,:)

    call rad_cnst_get_mode_num(0, bin_ndx, 'c', self%state, self%pbuf, num)
  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(modal_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: m,mm,l

    do m = 1, aero_props%nbins()
       mm = aero_props%indexer(m, 0)
       call self%get_ambient_num(m, raer(mm)%fld)
       call self%get_cldbrne_num(m, qqcw(mm)%fld)
       do l = 1, aero_props%nspecies(m)
          mm = aero_props%indexer(m, l)
          call self%get_ambient_mmr(l,m, raer(mm)%fld)
          call self%get_cldbrne_mmr(l,m, qqcw(mm)%fld)
       end do
    end do

  end subroutine get_states

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght(self, bin_ndx, ncol, nlev, species_type, use_preexisting_ice, wght)
    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
    real(r8), intent(out) :: wght(:,:)

    character(len=32) :: modetype
    real(r8), pointer :: dgnum(:,:,:)    ! mode dry radius
    real(r8) :: sigmag_aitken
    integer :: i,k

    call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)

    wght = 0._r8

    select case ( trim(species_type) )
    case('dust')
       if (modetype=='coarse' .or. modetype=='coarse_dust') then
          wght(:ncol,:) = 1._r8
       end if
    case('sulfate')
       if (modetype=='aitken') then
          if ( use_preexisting_ice ) then
             wght(:ncol,:) = 1._r8
          else
             call rad_cnst_get_mode_props(0, bin_ndx, sigmag=sigmag_aitken)
             call pbuf_get_field(self%pbuf, pbuf_get_index('DGNUM' ), dgnum)
             do k = 1,nlev
                do i = 1,ncol
                   if (dgnum(i,k,bin_ndx) > 0._r8) then
                      ! only allow so4 with D>0.1 um in ice nucleation
                      wght(i,k) = max(0._r8,(0.5_r8 - 0.5_r8* &
                           erf(log(0.1e-6_r8/dgnum(i,k,bin_ndx))/ &
                           (2._r8**0.5_r8*log(sigmag_aitken)))  ))
                   end if
                end do
             end do
          endif
       endif
    case('black-c')
       if (modetype=='accum') then
          wght(:ncol,:) = 1._r8
       endif
    case('sulfate_strat')
       if (modetype=='accum' .or. modetype=='coarse' .or. modetype=='coarse_strat') then
          wght(:ncol,:) = 1._r8
       endif
    end select

  end subroutine icenuc_size_wght

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine icenuc_type_wght(self, bin_ndx, ncol, nlev, species_type, aero_props, wght)

    use aerosol_properties_mod, only: aerosol_properties

    class(modal_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol
    integer, intent(in) :: nlev
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props

    real(r8), intent(out) :: wght(:,:)

    character(len=32) :: modetype

    call rad_cnst_get_info(0, bin_ndx, mode_type=modetype)

    wght = 0._r8

    if (species_type == 'dust') then
       if (modetype=='coarse_dust') then
          wght(:ncol,:) = 1._r8
       else
          call self%icenuc_type_wght_base(bin_ndx, ncol, nlev, species_type, aero_props, wght)
       end if
    else if (species_type == 'sulfate_strat') then
       if (modetype=='accum') then
          wght(:ncol,:) = 1._r8
       elseif ( modetype=='coarse' .or. modetype=='coarse_strat') then
          call self%icenuc_type_wght_base(bin_ndx, ncol, nlev, species_type, aero_props, wght)
       endif
    else
       wght(:ncol,:) = 1._r8
    end if

  end subroutine icenuc_type_wght

end module modal_aerosol_state_mod
