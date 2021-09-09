module modal_cam_aerosol_data_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_aerosol_data_mod, only: cam_aerosol_data
  use aerosol_data_mod, only: ptr2d_t
  use cam_abortutils,   only: endrun
  use cam_logfile,      only: iulog
  use ppgrid,           only: pcols
  use rad_constituents, only: rad_cnst_get_aer_mmr, rad_cnst_get_aer_props
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_props, rad_cnst_get_mode_num
  use physconst,        only: pi

  implicit none

  type, extends(cam_aerosol_data) :: modal_cam_aerosol_data
     integer               :: ntot_amode     ! number of aerosol modes
     integer,  allocatable :: nspec_amode(:) ! number of chemical species in each aerosol mode
     real(r8), allocatable :: sigmag_amode(:)! geometric standard deviation for each aerosol mode
     real(r8), allocatable :: dgnumlo_amode(:)
     real(r8), allocatable :: dgnumhi_amode(:)
     real(r8), allocatable :: voltonumblo_amode(:)
     real(r8), allocatable :: voltonumbhi_amode(:)
  contains
     procedure :: initialize => modal_initialize
     procedure :: loadaer => modal_loadaer
     procedure :: set_ptrs => modal_set_ptrs
  end type modal_cam_aerosol_data

contains

  subroutine modal_initialize(self)
    class(modal_cam_aerosol_data), intent(inout) :: self

    real(r8),allocatable :: alogsig(:)
    integer :: m

    ! get info about the modal aerosols
    ! get ntot_amode
    call rad_cnst_get_info(0, nmodes=self%ntot_amode)

    allocate( &
         self%nspec_amode(self%ntot_amode),  &
         self%sigmag_amode(self%ntot_amode), &
         self%dgnumlo_amode(self%ntot_amode), &
         self%dgnumhi_amode(self%ntot_amode), &
         self%voltonumblo_amode(self%ntot_amode), &
         self%voltonumbhi_amode(self%ntot_amode)  )

    allocate( alogsig(self%ntot_amode) )

    do m = 1, self%ntot_amode
       ! use only if width of size distribution is prescribed

       ! get mode info
       call rad_cnst_get_info(0, m, nspec=self%nspec_amode(m))

       ! get mode properties
       call rad_cnst_get_mode_props(0, m, sigmag=self%sigmag_amode(m),  &
            dgnumhi=self%dgnumhi_amode(m), dgnumlo=self%dgnumlo_amode(m))

       alogsig(m) = log(self%sigmag_amode(m))

       self%voltonumblo_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
            (self%dgnumlo_amode(m)**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )
       self%voltonumbhi_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
            (self%dgnumhi_amode(m)**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )

    end do

    deallocate( alogsig )

  end subroutine modal_initialize

  subroutine modal_loadaer( self, istart, istop, k, m, cs, phase, naerosol, vaerosol, hygro)

    ! return aerosol number, volume concentrations, and bulk hygroscopicity

    class(modal_cam_aerosol_data), intent(in) :: self
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

    real(r8) :: vol(pcols) ! aerosol volume mixing ratio
    integer  :: i, l
    !-------------------------------------------------------------------------------

    do i = istart, istop
       vaerosol(i) = 0._r8
       hygro(i)    = 0._r8
    end do

    do l = 1, self%nspec_amode(m)

       call rad_cnst_get_aer_mmr(0, m, l, 'a', self%state, self%pbuf, raer)
       call rad_cnst_get_aer_mmr(0, m, l, 'c', self%state, self%pbuf, qqcw)
       call rad_cnst_get_aer_props(0, m, l, density_aer=specdens, hygro_aer=spechygro)

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
          write(iulog,*)'phase=',phase,' in loadaer'
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
    call rad_cnst_get_mode_num(0, m, 'a', self%state, self%pbuf, raer)
    call rad_cnst_get_mode_num(0, m, 'c', self%state, self%pbuf, qqcw)
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

    ! adjust number so that dgnumlo < dgnum < dgnumhi
    do i = istart, istop
       naerosol(i) = max(naerosol(i), vaerosol(i)*self%voltonumbhi_amode(m))
       naerosol(i) = min(naerosol(i), vaerosol(i)*self%voltonumblo_amode(m))
    end do

  end subroutine modal_loadaer

  subroutine modal_set_ptrs( self, indexer, raer, qqcw )
    class(modal_cam_aerosol_data), intent(in) :: self
    integer,       intent(in) :: indexer(:,0:)
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: m,mm,l

    do m = 1, self%ntot_amode
       mm = indexer(m, 0)
       call rad_cnst_get_mode_num(0, m, 'a', self%state, self%pbuf, raer(mm)%fld)
       call rad_cnst_get_mode_num(0, m, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
       do l = 1, self%nspec_amode(m)
          mm = indexer(m, l)
          call rad_cnst_get_aer_mmr(0, m, l, 'a', self%state, self%pbuf, raer(mm)%fld)
          call rad_cnst_get_aer_mmr(0, m, l, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
       end do
    end do

  end subroutine modal_set_ptrs

end module modal_cam_aerosol_data_mod
