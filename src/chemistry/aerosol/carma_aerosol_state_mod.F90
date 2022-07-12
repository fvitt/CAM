module carma_aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_state_mod, only: aerosol_state, ptr2d_t

  use rad_constituents, only: rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_num, rad_cnst_get_bin_mmr
  use rad_constituents, only: rad_cnst_get_info_by_bin
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties
  use cam_abortutils, only: endrun

  implicit none

  private

  public :: carma_aerosol_state

  type, extends(aerosol_state) :: carma_aerosol_state
     private
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
   contains

     procedure :: get_transported
     procedure :: set_transported
     procedure :: ambient_total_bin_mmr
     procedure :: get_ambient_mmr
     procedure :: get_cldbrne_mmr
     procedure :: get_ambient_num
     procedure :: get_cldbrne_num
     procedure :: get_states
     procedure :: icenuc_size_wght1
     procedure :: icenuc_size_wght2

     final :: destructor

  end type carma_aerosol_state

  interface carma_aerosol_state
     procedure :: constructor
  end interface carma_aerosol_state

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(state,pbuf) result(newobj)
    type(physics_state), target, optional :: state
    type(physics_buffer_desc), pointer, optional :: pbuf(:)

    type(carma_aerosol_state), pointer :: newobj

    character(len=*),parameter :: prefix = 'carma_aerosol_state::constructor: '
    integer :: ierr

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       call endrun(prefix//'error allocating newobj')
    end if

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
  ! sets transported components
  ! This aerosol model with the state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine set_transported( self, transported_array )
    class(carma_aerosol_state), intent(inout) :: self
    real(r8), intent(in) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine set_transported

  !------------------------------------------------------------------------------
  ! returns transported components
  ! This returns to current state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine get_transported( self, transported_array )
    class(carma_aerosol_state), intent(in) :: self
    real(r8), intent(out) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine get_transported

  !------------------------------------------------------------------------
  ! Total aerosol mass mixing ratio for a bin in a given grid box location (column and layer)
  !------------------------------------------------------------------------
  function ambient_total_bin_mmr(self, aero_props, bin_ndx, col_ndx, lyr_ndx) result(mmr_tot)
    class(carma_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    integer, intent(in) :: bin_ndx      ! bin index
    integer, intent(in) :: col_ndx      ! column index
    integer, intent(in) :: lyr_ndx      ! vertical layer index

    real(r8) :: mmr_tot                 ! mass mixing ratios totaled for all species

    real(r8),pointer :: mmrptr(:,:)

    call rad_cnst_get_bin_mmr(0, bin_ndx, 'a', self%state, self%pbuf, mmrptr)

    mmr_tot = mmrptr(col_ndx,lyr_ndx)

  end function ambient_total_bin_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given species index and bin index
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
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
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
  ! returns ambient aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, num)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer :: num(:,:)      ! number mixing ratios

    call rad_cnst_get_bin_num(0, bin_ndx, 'a', self%state, self%pbuf, num)
  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, num)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer :: num(:,:)      ! number mixing ratios

    call rad_cnst_get_bin_num(0, bin_ndx, 'c', self%state, self%pbuf, num)
  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  ! returns interstitial and cloud-borne aerosol states
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(carma_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: ibin,ispc, indx

    do ibin = 1, aero_props%nbins()
       indx = aero_props%indexer(ibin, 0)
       call self%get_ambient_num(ibin, raer(indx)%fld)
       call self%get_cldbrne_num(ibin, qqcw(indx)%fld)
       indx = aero_props%indexer(ibin, 1)
       call self%get_ambient_mmr(0,ibin, raer(indx)%fld)
       call self%get_cldbrne_mmr(0,ibin, qqcw(indx)%fld)
       do ispc = 1, aero_props%nspecies(ibin)
          indx = aero_props%indexer(ibin, ispc+1)
          call self%get_ambient_mmr(ispc,ibin, raer(indx)%fld)
          call self%get_cldbrne_mmr(ispc,ibin, qqcw(indx)%fld)
       end do
    end do

  end subroutine get_states

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght1(self, bin_ndx, ncol, nlev, species_type, use_preexisting_ice, wght)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
    real(r8), intent(out) :: wght(:,:)

    character(len=32) :: bin_name
    real(r8), pointer :: dryr(:,:)
    integer :: i,k
    real(r8) :: diamdry

    wght = 0._r8

    call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)
    call pbuf_get_field(self%pbuf, pbuf_get_index(trim(bin_name)//"_dryr"),dryr)

    do k = 1,nlev
       do i = 1,ncol
          diamdry = dryr(i,k) * 2.e4_r8  ! diameter in microns (from radius in cm)
          if (diamdry >= 0.1_r8) then ! size threashold
             wght(i,k) = 1._r8
          end if
       end do
    end do

  end subroutine icenuc_size_wght1

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin, column and verical layer
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght2(self, bin_ndx, col_ndx, lyr_ndx, species_type, use_preexisting_ice, wght)
    class(carma_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice    ! pre-existing ice flag
    real(r8), intent(out) :: wght

    character(len=32) :: bin_name
    real(r8), pointer :: dryr(:,:)
    real(r8) :: diamdry

    wght = 0._r8

    call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)
    call pbuf_get_field(self%pbuf, pbuf_get_index(trim(bin_name)//"_dryr"),dryr)

    diamdry = dryr(col_ndx,lyr_ndx) * 2.e4_r8  ! diameter in microns (from radius in cm)
    if (diamdry >= 0.1_r8) then ! size threashold
       wght = 1._r8
    end if

  end subroutine icenuc_size_wght2

end module carma_aerosol_state_mod
