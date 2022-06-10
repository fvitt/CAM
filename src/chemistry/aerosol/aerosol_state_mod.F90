module aerosol_state_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use cam_logfile, only: iulog
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

  public :: aerosol_state
  public :: ptr2d_t

  !> aerosol_state defines the interface to the time-varying aerosol state
  !! variables (e.g., mixing ratios, number concentrations). This includes the
  !! aerosol portion of the overall model state.
  !!
  !! Each aerosol package (e.g., MAM, CARMA, etc) must extend the aerosol_state
  !! class to allow access to the state information (transported and not transported)
  !! of the aerosol package. Any package must implement each of the deferred
  !! procedures of the abstract aerosol_state class, may include additional private
  !! data members and type-bound procedures, and may override functions of the
  !! abstract class.
  !!
  !! Please see the modal_aerosol_state module for an example of how the aerosol_state
  !! class can be extended for a specific aerosol package.
  type, abstract :: aerosol_state
   contains
     procedure :: number_transported
     procedure :: get_transported
     procedure :: set_transported
     procedure(aero_get_amb_total_bin_mmr), deferred :: ambient_total_bin_mmr
     procedure(aero_get_state_mmr), deferred :: get_ambient_mmr
     procedure(aero_get_state_mmr), deferred :: get_cldbrne_mmr
     procedure(aero_get_state_num), deferred :: get_ambient_num
     procedure(aero_get_state_num), deferred :: get_cldbrne_num
     procedure(aero_get_states), deferred :: get_states
     procedure, private :: loadaer1
     procedure, private :: loadaer2
     generic :: loadaer => loadaer1, loadaer2
     procedure(aero_icenuc_size_wght1), deferred :: icenuc_size_wght1
     procedure(aero_icenuc_size_wght2), deferred :: icenuc_size_wght2
     generic :: icenuc_size_wght => icenuc_size_wght1,icenuc_size_wght2
     procedure :: icenuc_type_wght_base
     procedure :: icenuc_type_wght => icenuc_type_wght_base
 end type aerosol_state

  ! for state fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  interface

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     function aero_get_amb_total_bin_mmr(self, aero_props, bin_ndx, col_ndx, lyr_ndx) result(mmr_tot)
       import
       class(aerosol_state), intent(in) :: self
       class(aerosol_properties), intent(in) :: aero_props
       integer, intent(in) :: bin_ndx      ! bin index
       integer, intent(in) :: col_ndx      ! column index
       integer, intent(in) :: lyr_ndx      ! vertical layer index

       real(r8) :: mmr_tot                 ! mass mixing ratios totaled for all species

     end function aero_get_amb_total_bin_mmr

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     subroutine aero_get_state_mmr(self, species_ndx, bin_ndx, mmr)
       import
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: species_ndx  ! species index
       integer, intent(in) :: bin_ndx      ! bin index
       real(r8), pointer :: mmr(:,:)       ! mass mixing ratios
     end subroutine aero_get_state_mmr

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     subroutine aero_get_state_num(self, bin_ndx, num)
       import
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx     ! bin index
       real(r8), pointer   :: num(:,:)    ! number densities

     end subroutine aero_get_state_num

     !------------------------------------------------------------------------
     !------------------------------------------------------------------------
     subroutine aero_get_states( self, aero_props, raer, qqcw )
       import

       class(aerosol_state), intent(in) :: self
       class(aerosol_properties), intent(in) :: aero_props ! properties of the aerosol model
       type(ptr2d_t), intent(out) :: raer(:) ! state of interstitual aerosols
       type(ptr2d_t), intent(out) :: qqcw(:) ! state of cloud-borne aerosols

     end subroutine aero_get_states

     !------------------------------------------------------------------------------
     !------------------------------------------------------------------------------
     subroutine aero_icenuc_size_wght1(self, bin_ndx, ncol, nlev, species_type, use_preexisting_ice, wght)
       import
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx           ! bin number
       integer, intent(in) :: ncol                ! number of columns
       integer, intent(in) :: nlev                ! number of vertical levels
       character(len=*), intent(in) :: species_type  ! species type
       logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
       real(r8), intent(out) :: wght(:,:)

     end subroutine aero_icenuc_size_wght1

     !------------------------------------------------------------------------------
     !------------------------------------------------------------------------------
     subroutine aero_icenuc_size_wght2(self, bin_ndx, col_ndx, lyr_ndx, species_type, use_preexisting_ice, wght)
       import
       class(aerosol_state), intent(in) :: self
       integer, intent(in) :: bin_ndx                ! bin number
       integer, intent(in) :: col_ndx                ! column index
       integer, intent(in) :: lyr_ndx                ! vertical layer index
       character(len=*), intent(in) :: species_type  ! species type
       logical, intent(in) :: use_preexisting_ice    ! pre-existing ice flag
       real(r8), intent(out) :: wght

     end subroutine aero_icenuc_size_wght2

  end interface

contains

  !------------------------------------------------------------------------------
  ! returns number of transported elements
  !------------------------------------------------------------------------------
  integer function number_transported(self)
    class(aerosol_state), intent(in) :: self
    ! to be implemented later
    number_transported = -1
  end function number_transported

  !------------------------------------------------------------------------------
  ! sets transported components
  !------------------------------------------------------------------------------
  subroutine set_transported( self, transported_array )
    class(aerosol_state), intent(inout) :: self
    real(r8), intent(in) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine set_transported

  !------------------------------------------------------------------------------
  ! returns transported components
  !------------------------------------------------------------------------------
  subroutine get_transported( self, transported_array )
    class(aerosol_state), intent(in) :: self
    real(r8), intent(out) :: transported_array(:,:,:)
    ! to be implemented later
  end subroutine get_transported

  !------------------------------------------------------------------------------
  ! returns aerosol number, volume concentrations, and bulk hygroscopicity
  !------------------------------------------------------------------------------
  subroutine loadaer1( self, aero_props, istart, istop, k,  m, cs, phase, &
                       naerosol, vaerosol, hygro)

    use aerosol_properties_mod, only: aerosol_properties
    use modal_aerosol_properties_mod, only: modal_aerosol_properties

    ! return aerosol number, volume concentrations, and bulk hygroscopicity

    ! input arguments
    class(aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer,  intent(in) :: istart      ! start column index (1 <= istart <= istop <= pcols)
    integer,  intent(in) :: istop       ! stop column index
    integer,  intent(in) :: m           ! mode or bin index
    integer,  intent(in) :: k           ! level index
    real(r8), intent(in) :: cs(:,:)     ! air density (kg/m3)
    integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

    ! output arguments
    real(r8), intent(out) :: naerosol(:)  ! number conc (1/m3)
    real(r8), intent(out) :: vaerosol(:)  ! volume conc (m3/m3)
    real(r8), intent(out) :: hygro(:)     ! bulk hygroscopicity of mode

    ! internal
    real(r8), pointer :: raer(:,:) ! interstitial aerosol mass, number mixing ratios
    real(r8), pointer :: qqcw(:,:) ! cloud-borne aerosol mass, number mixing ratios
    real(r8) :: specdens, spechygro

    real(r8) :: vol(istart:istop) ! aerosol volume mixing ratio
    integer  :: i, l
    !-------------------------------------------------------------------------------

    do i = istart, istop
       vaerosol(i) = 0._r8
       hygro(i)    = 0._r8
    end do

    do l = 1, aero_props%nspecies(m)

       call self%get_ambient_mmr(l,m, raer)
       call self%get_cldbrne_mmr(l,m, qqcw)
       call aero_props%get(m,l, density=specdens, hygro=spechygro)

       if (phase == 3) then
          do i = istart, istop
             vol(i) = max(raer(i,k) + qqcw(i,k), 0._r8)/specdens
          end do
       else if (phase == 2) then
          do i = istart, istop
             vol(i) = max(qqcw(i,k), 0._r8)/specdens
          end do
       else if (phase == 1) then
          do i = istart, istop
             vol(i) = max(raer(i,k), 0._r8)/specdens
          end do
       else
          write(iulog,*)'phase = ',phase,' in loadaer not recognized'
          call endrun('phase error in loadaer')
       end if

       do i = istart, istop
          vaerosol(i) = vaerosol(i) + vol(i)
          hygro(i)    = hygro(i) + vol(i)*spechygro
       end do

    end do

    do i = istart, istop
       if (vaerosol(i) > 1.0e-30_r8) then   ! +++xl add 8/2/2007
          hygro(i)    = hygro(i)/(vaerosol(i))
          vaerosol(i) = vaerosol(i)*cs(i,k)
       else
          hygro(i)    = 0.0_r8
          vaerosol(i) = 0.0_r8
       end if
    end do

    ! aerosol number
    call self%get_ambient_num(m, raer)
    call self%get_cldbrne_num(m, qqcw)
    if (phase == 3) then
       do i = istart, istop
          naerosol(i) = (raer(i,k) + qqcw(i,k))*cs(i,k)
       end do
    else if (phase == 2) then
       do i = istart, istop
          naerosol(i) = qqcw(i,k)*cs(i,k)
       end do
    else
       do i = istart, istop
          naerosol(i) = raer(i,k)*cs(i,k)
       end do
    end if

    select type(aero_props)
    type is (modal_aerosol_properties)

       ! adjust number so that dgnumlo < dgnum < dgnumhi not done for bins
       do i = istart, istop
          naerosol(i) = max(naerosol(i), vaerosol(i)*aero_props%voltonumbhi(m))
          naerosol(i) = min(naerosol(i), vaerosol(i)*aero_props%voltonumblo(m))
       end do
    end select

  end subroutine loadaer1

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine loadaer2( self, aero_props, i, k, m, cs, phase, &
                       naerosol, vaerosol, hygro )

    use aerosol_properties_mod, only: aerosol_properties
    use ppgrid, only: pcols, pver

    ! return aerosol number, volume concentrations, and bulk hygroscopicity

    ! input arguments
    class(aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer,  intent(in) :: i           ! column index
    integer,  intent(in) :: k           ! level index
    integer,  intent(in) :: m           ! mode or bin index
    real(r8), intent(in) :: cs          ! air density (kg/m3)
    integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

    ! output arguments
    real(r8), intent(out) :: naerosol  ! number conc (1/m3)
    real(r8), intent(out) :: vaerosol  ! volume conc (m3/m3)
    real(r8), intent(out) :: hygro     ! bulk hygroscopicity of mode

    real(r8) :: cs_a(pcols,pver)          ! air density (kg/m3)
    real(r8) :: naerosol_a(pcols)  ! number conc (1/m3)
    real(r8) :: vaerosol_a(pcols)  ! volume conc (m3/m3)
    real(r8) :: hygro_a(pcols)     ! bulk hygroscopicity of mode

    cs_a(i,k) = cs

    call self%loadaer(aero_props, i, i, k, m, cs_a, phase, naerosol_a, vaerosol_a, hygro_a)

    naerosol = naerosol_a(i)
    vaerosol = vaerosol_a(i)
    hygro = hygro_a(i)

  end subroutine loadaer2

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine icenuc_type_wght_base(self, bin_ndx, ncol, nlev, species_type, aero_props, wght)

    use aerosol_properties_mod, only: aerosol_properties

    class(aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol
    integer, intent(in) :: nlev
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props

    real(r8), intent(out) :: wght(:,:)

    real(r8) :: mmr(ncol,nlev)
    real(r8) :: totalmmr(ncol,nlev)
    real(r8), pointer :: aer_bin(:,:)

    character(len=32) :: spectype, sptype
    integer :: l

    wght(:,:) = 0._r8
    totalmmr(:,:) = 0._r8
    mmr(:,:)   = 0._r8

    if (species_type=='sulfate_strat') then
       sptype = 'sulfate'
    else
       sptype = species_type
    end if

    do l = 1, aero_props%nspecies(bin_ndx)

       call self%get_ambient_mmr(l, bin_ndx, aer_bin)
       call aero_props%species_type(bin_ndx, l, spectype=spectype)

       totalmmr(:ncol,:) = totalmmr(:ncol,:) + aer_bin(:ncol,:)

       if (trim(spectype) == trim(sptype)) then
          mmr(:ncol,:) = mmr(:ncol,:) + aer_bin(:ncol,:)
       end if

    end do

    where (totalmmr(:ncol,:) > 0._r8)
       wght(:ncol,:) = mmr(:ncol,:)/totalmmr(:ncol,:)
    end where

  end subroutine icenuc_type_wght_base

end module aerosol_state_mod
