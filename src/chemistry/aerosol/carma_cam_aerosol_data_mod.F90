module carma_cam_aerosol_data_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_aerosol_data_mod, only: cam_aerosol_data
  use aerosol_data_mod, only: ptr2d_t
  use cam_abortutils,   only: endrun
  use cam_logfile,      only: iulog
  use ppgrid,           only: pcols
  use rad_constituents, only: rad_cnst_get_bin_props_by_idx, rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_num
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_info_by_bin, rad_cnst_get_bin_mmr
  use rad_constituents, only: rad_cnst_get_info_by_bin_spec
  use constituents,     only: cnst_get_ind

  implicit none

  type, extends(cam_aerosol_data) :: carma_cam_aerosol_data
     integer, allocatable :: nspec(:)
     integer :: nbins
  contains
     procedure :: initialize => carma_initialize
     procedure :: loadaer => carma_loadaer
     procedure :: set_ptrs => carma_set_ptrs
  end type carma_cam_aerosol_data

contains

  subroutine carma_initialize(self)
    class(carma_cam_aerosol_data), intent(inout) :: self

    integer :: l, m, mm
    character(len=32) :: tmpname
    character(len=32) :: tmpname_cw
    integer :: nspec_max
    integer :: idxtmp

    ! get info about the modal aerosols
    ! get nbins

    call rad_cnst_get_info( 0, nbins=self%nbins)

    allocate( self%nspec(self%nbins) )

    self%ncnst_tot = 0

    do m = 1, self%nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=self%nspec(m))
       self%ncnst_tot =  self%ncnst_tot + self%nspec(m) + 2
    end do

    allocate( self%fieldname(self%ncnst_tot) )
    allocate( self%fieldname_cw(self%ncnst_tot) )

    nspec_max = maxval(self%nspec) + 1
    allocate( self%cnstndx(self%nbins,0:nspec_max) )

    self%cnstndx = -1
    mm = 0

    do m = 1, self%nbins
       do l = 0, self%nspec(m) + 1  ! loop over bin + aerosol constituents
          if (l == 0) then   ! number
             call rad_cnst_get_info_by_bin(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
          else if (l == 1) then
             call rad_cnst_get_info_by_bin(0, m,  mmr_name=tmpname, mmr_name_cw=tmpname_cw)
          else
             call rad_cnst_get_info_by_bin_spec(0, m, l-1, spec_name=tmpname, spec_name_cw=tmpname_cw)
          end if

          mm = mm+1

          self%fieldname(mm)    = trim(tmpname) // '_mixnuc1'
          self%fieldname_cw(mm) = trim(tmpname_cw) // '_mixnuc1'

          call cnst_get_ind(tmpname, idxtmp, abort=.false.)
          if (idxtmp>0) then
             self%cnstndx(m,l) = idxtmp
             self%lq(idxtmp) = .true.
          end if
       end do
    end do
  end subroutine carma_initialize


  subroutine carma_loadaer( self, istart, istop, k, m, cs, phase, naerosol, vaerosol, hygro)

    ! return aerosol number, volume concentrations, and bulk hygroscopicity

    class(carma_cam_aerosol_data), intent(in) :: self
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

    do l = 1, self%nspec(m)

       call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'a', self%state, self%pbuf, raer)
       call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'c', self%state, self%pbuf, qqcw)
       call rad_cnst_get_bin_props_by_idx(0, m, l, density_aer=specdens, hygro_aer=spechygro)

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
    call rad_cnst_get_bin_num(0, m, 'a', self%state, self%pbuf, raer)
    call rad_cnst_get_bin_num(0, m, 'c', self%state, self%pbuf, qqcw)  ! cloud-borne aerosol
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

  end subroutine carma_loadaer

  subroutine carma_set_ptrs( self, indexer, raer, qqcw )
    class(carma_cam_aerosol_data), intent(in) :: self
    integer,       intent(in) :: indexer(:,0:)
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: m, mm, l

    do m = 1, self%nbins
       mm = indexer(m, 0)
       call rad_cnst_get_bin_num(0, m, 'a', self%state, self%pbuf, raer(mm)%fld)
       call rad_cnst_get_bin_num(0, m, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
       mm = indexer(m, 1)
       call rad_cnst_get_bin_mmr(0, m, 'a', self%state, self%pbuf, raer(mm)%fld)
       call rad_cnst_get_bin_mmr(0, m, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
       do l = 2, self%nspec(m)+1
          mm = indexer(m, l)
          call rad_cnst_get_bin_mmr_by_idx(0, m, l-1, 'a', self%state, self%pbuf, raer(mm)%fld)
          call rad_cnst_get_bin_mmr_by_idx(0, m, l-1, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
       end do
    end do

  end subroutine carma_set_ptrs

end module carma_cam_aerosol_data_mod
