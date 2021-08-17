module carma_aerosol_model_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_model_mod, only: aerosol_model, twothird, sq2, ptr2d_t
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_info_by_bin, rad_cnst_get_info_by_bin_spec
  use rad_constituents, only: rad_cnst_get_bin_props_by_idx, rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_num
  use rad_constituents, only: rad_cnst_get_bin_mmr
  use physconst, only: pi
  use constituents,     only: cnst_get_ind

  implicit none

  type, extends(aerosol_model) :: carma_aerosol_model
     integer, allocatable :: nspec(:)
   contains
     procedure :: loadaer => carma_loadaer
     procedure :: set_ptrs => carma_set_ptrs
     procedure :: model_init => carma_model_init
     procedure :: model_final => carma_model_final
     procedure :: err_funct => carma_err_funct ! override base routine
  end type carma_aerosol_model

contains

  subroutine carma_model_init(self)
    class(carma_aerosol_model), intent(inout) :: self

    integer :: l, m, mm, nbins, nspec_max
    character(len=32) :: tmpname
    character(len=32) :: tmpname_cw
    integer :: idxtmp

    self%model_name = 'carma'

    ! get info about the modal aerosols
    ! get nbins

    call rad_cnst_get_info( 0, nbins=nbins)

    self%mtotal = nbins

    allocate( self%nspec(nbins) )
    allocate( self%nmasses(nbins) )
    allocate( self%amcubecoef(nbins) )
    allocate( self%exp45logsig(nbins) )
    allocate( self%argfactor(nbins) )
    allocate( self%alogsig(nbins) )
    allocate( self%f1(nbins) )
    allocate( self%f2(nbins) )

    self%ncnst_tot = 0

    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=self%nspec(m))
       self%nmasses(m) = self%nspec(m) + 1
       self%ncnst_tot =  self%ncnst_tot + self%nspec(m) + 2
       self%amcubecoef(m)=3._r8/(4._r8*pi)
       self%argfactor(m)=twothird/(sq2*log(2._r8))
       self%exp45logsig(m)=1._r8
       self%alogsig(m)=1._r8
       self%f1(m)=1._r8
       self%f2(m)=1._r8
    end do
    ! add plus one to include number, total mmr and nspec
    nspec_max = maxval(self%nspec) + 1
    allocate( self%indexer(nbins,0:nspec_max) )
    allocate( self%cnstndx(nbins,0:nspec_max) )
    allocate( self%fieldname(self%ncnst_tot) )
    allocate( self%fieldname_cw(self%ncnst_tot) )

    self%cnstndx = -1

    mm = 0

    do m = 1, nbins
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
  end subroutine carma_model_init

  subroutine carma_model_final(self)
    class(carma_aerosol_model), intent(inout) :: self

    deallocate( self%nspec )
    deallocate( self%nmasses )
    deallocate( self%amcubecoef )
    deallocate( self%argfactor )
    deallocate( self%alogsig )
    deallocate( self%f1 )
    deallocate( self%f2 )
    deallocate( self%indexer )
    deallocate( self%cnstndx )
    deallocate( self%fieldname_cw )

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

  subroutine carma_set_ptrs( self, raer, qqcw )
    class(carma_aerosol_model), intent(in) :: self
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: m, mm, l

    do m = 1, self%mtotal
       mm = self%indexer(m, 0)
       call rad_cnst_get_bin_num(0, m, 'a', self%state, self%pbuf, raer(mm)%fld)
       call rad_cnst_get_bin_num(0, m, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
       mm = self%indexer(m, 1)
       call rad_cnst_get_bin_mmr(0, m, 'a', self%state, self%pbuf, raer(mm)%fld)
       call rad_cnst_get_bin_mmr(0, m, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
       do l = 2, self%nspec(m)+1
          mm = self%indexer(m, l)
          call rad_cnst_get_bin_mmr_by_idx(0, m, l-1, 'a', self%state, self%pbuf, raer(mm)%fld)
          call rad_cnst_get_bin_mmr_by_idx(0, m, l-1, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
       end do
    end do

  end subroutine carma_set_ptrs

  ! override the error function
  function carma_err_funct(self,x) result(err)
    class(carma_aerosol_model), intent(in) :: self
    real(r8), intent(in) :: x
    real(r8) :: err
    err = -1._r8
  end function carma_err_funct

end module carma_aerosol_model_mod
