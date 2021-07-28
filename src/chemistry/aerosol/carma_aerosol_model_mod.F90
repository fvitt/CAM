module carma_aerosol_model_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_model_mod, only: aerosol_model, twothird, sq2
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_info_by_bin
  use rad_constituents, only: rad_cnst_get_bin_props_by_idx, rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_num
  use physconst, only: pi

  implicit none

  type, extends(aerosol_model) :: carma_aerosol_model
     integer, allocatable :: nspec(:)
   contains
     procedure :: loadaer => carma_loadaer
     procedure :: model_init => carma_model_init
     procedure :: model_final => carma_model_final
     procedure :: err_funct => carma_err_funct ! override base routine
  end type carma_aerosol_model

contains

  subroutine carma_model_init(self)
    class(carma_aerosol_model), intent(inout) :: self

    integer :: m, nbins

    if (masterproc) write(iulog,*) 'AERO_MODEL: Init CARMA data...'

    ! get info about the modal aerosols
    ! get nbins

    call rad_cnst_get_info( 0, nbins=nbins)

    self%mtotal = nbins

    allocate( self%nspec(nbins) )
    allocate( self%amcubecoef(nbins) )
    allocate( self%argfactor(nbins) )
    allocate( self%alogsig(nbins) )
    allocate( self%f1(nbins) )
    allocate( self%f2(nbins) )

    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=self%nspec(m))
       self%amcubecoef(m)=3._r8/(4._r8*pi)
       self%argfactor(m)=twothird/(sq2*log(2._r8))
       self%alogsig(m)=1._r8
       self%f1(m)=1._r8
       self%f2(m)=1._r8
    end do

  end subroutine carma_model_init

  subroutine carma_model_final(self)
    class(carma_aerosol_model), intent(inout) :: self

    deallocate( self%nspec )
    deallocate( self%amcubecoef )
    deallocate( self%argfactor )

  end subroutine carma_model_final

  subroutine carma_loadaer( self, istart, istop, k, m, cs, phase, naerosol, vaerosol, hygro)
    use cam_abortutils,   only: endrun
    use cam_logfile,      only: iulog
    use ppgrid,           only: pcols

    ! return aerosol number, volume concentrations, and bulk hygroscopicity

    class(carma_aerosol_model), intent(in) :: self
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

    if (masterproc) write(iulog,*) 'AERO_MODEL: CARMA load aerosol data..'

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

  ! override the error function
  function carma_err_funct(self,x) result(err)
    class(carma_aerosol_model), intent(in) :: self
    real(r8), intent(in) :: x
    real(r8) :: err
    err = -1._r8
  end function carma_err_funct

end module carma_aerosol_model_mod
