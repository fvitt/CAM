module modal_aerosol_model_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc
  use physconst,        only: pi
  use rad_constituents, only: rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, rad_cnst_get_mode_num
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_props
  use rad_constituents, only: rad_cnst_get_mam_mmr_idx, rad_cnst_get_mode_num_idx
  use ppgrid,           only: pcols
  use cam_abortutils,   only: endrun
  use phys_control,     only: phys_getopts
  use aerosol_model_mod, only: aerosol_model, twothird, sq2, ptr2d_t

  implicit none
  private
  public :: modal_aerosol_model

  type, extends(aerosol_model) :: modal_aerosol_model
     integer               :: ntot_amode     ! number of aerosol modes
     integer,  allocatable :: nspec_amode(:) ! number of chemical species in each aerosol mode
     real(r8), allocatable :: sigmag_amode(:)! geometric standard deviation for each aerosol mode
     real(r8), allocatable :: dgnumlo_amode(:)
     real(r8), allocatable :: dgnumhi_amode(:)
     real(r8), allocatable :: voltonumblo_amode(:)
     real(r8), allocatable :: voltonumbhi_amode(:)
   contains
     procedure :: loadaer => modal_loadaer
     procedure :: set_ptrs => modal_set_ptrs
     procedure :: model_init => modal_model_init
     procedure :: model_final => modal_model_final
  end type modal_aerosol_model

contains

  subroutine modal_model_init(self)
    class(modal_aerosol_model), intent(inout) :: self

    integer :: m, nspec_max, l,idx, mm
    character(len=32)   :: tmpname
    character(len=32)   :: tmpname_cw

    self%model_name = 'modal'

    call phys_getopts(prog_modal_aero_out = self%prognostic)

    ! get info about the modal aerosols
    ! get ntot_amode
    call rad_cnst_get_info(0, nmodes=self%ntot_amode)
    self%mtotal = self%ntot_amode

    allocate( &
         self%nspec_amode(self%ntot_amode),  &
         self%nmasses(self%ntot_amode),  &
         self%sigmag_amode(self%ntot_amode), &
         self%dgnumlo_amode(self%ntot_amode), &
         self%dgnumhi_amode(self%ntot_amode), &
         self%alogsig(self%ntot_amode),      &
         self%exp45logsig(self%ntot_amode),  &
         self%f1(self%ntot_amode),           &
         self%f2(self%ntot_amode),           &
         self%voltonumblo_amode(self%ntot_amode), &
         self%voltonumbhi_amode(self%ntot_amode)  )

    allocate( self%amcubecoef(self%mtotal) )
    allocate( self%argfactor(self%mtotal) )

    self%ncnst_tot = 0

    do m = 1, self%ntot_amode
       ! use only if width of size distribution is prescribed

       ! get mode info
       call rad_cnst_get_info(0, m, nspec=self%nspec_amode(m))
       self%nmasses(m) = self%nspec_amode(m)
       self%ncnst_tot =  self%ncnst_tot + self%nspec_amode(m) + 1

       ! get mode properties
       call rad_cnst_get_mode_props(0, m, sigmag=self%sigmag_amode(m),  &
            dgnumhi=self%dgnumhi_amode(m), dgnumlo=self%dgnumlo_amode(m))

       self%alogsig(m)     = log(self%sigmag_amode(m))
       self%exp45logsig(m) = exp(4.5_r8*self%alogsig(m)*self%alogsig(m))
       self%f1(m)          = 0.5_r8*exp(2.5_r8*self%alogsig(m)*self%alogsig(m))
       self%f2(m)          = 1._r8 + 0.25_r8*self%alogsig(m)

       self%voltonumblo_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
            (self%dgnumlo_amode(m)**3._r8)*exp(4.5_r8*self%alogsig(m)**2._r8) )
       self%voltonumbhi_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
            (self%dgnumhi_amode(m)**3._r8)*exp(4.5_r8*self%alogsig(m)**2._r8) )

       self%amcubecoef(m)=3._r8/(4._r8*pi*self%exp45logsig(m))
       self%argfactor(m)=twothird/(sq2*self%alogsig(m))
    end do

    ! Find max number of species in all the modes
    nspec_max = maxval(self%nspec_amode)
    allocate( self%indexer(self%ntot_amode,0:nspec_max) )
    allocate( self%cnstndx(self%ntot_amode,0:nspec_max) )
    allocate( self%fieldname(self%ncnst_tot) )
    allocate( self%fieldname_cw(self%ncnst_tot) )

    self%cnstndx = -1
    mm=0

    do m = 1, self%ntot_amode
       do l = 0, self%nspec_amode(m)   ! loop over number + chem constituents

          mm = mm+1

          if (l == 0) then   ! number
             call rad_cnst_get_info(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
          else
             call rad_cnst_get_info(0, m, l, spec_name=tmpname, spec_name_cw=tmpname_cw)
          end if

          self%fieldname(mm)    = trim(tmpname) // '_mixnuc1'
          self%fieldname_cw(mm) = trim(tmpname_cw) // '_mixnuc1'

          if (self%prognostic) then

             ! To set tendencies in the ptend object need to get the constituent indices
             ! for the prognostic species
             if (l == 0) then   ! number
                call rad_cnst_get_mode_num_idx(m, idx)
             else
                call rad_cnst_get_mam_mmr_idx(m, l, idx)
             end if
             self%cnstndx(m,l) = idx
             self%lq(idx) = .true.

          endif

       end do
    end do

  end subroutine modal_model_init

  subroutine modal_model_final(self)
    class(modal_aerosol_model), intent(inout) :: self

    deallocate( &
         self%nspec_amode,  &
         self%nmasses, &
         self%sigmag_amode, &
         self%dgnumlo_amode, &
         self%dgnumhi_amode, &
         self%alogsig,      &
         self%exp45logsig,  &
         self%f1,           &
         self%f2,           &
         self%voltonumblo_amode, &
         self%voltonumbhi_amode )
    deallocate( self%amcubecoef )
    deallocate( self%argfactor )
    deallocate( self%indexer )
    deallocate( self%cnstndx )
    deallocate( self%fieldname )
    deallocate( self%fieldname_cw )

  end subroutine modal_model_final


  subroutine modal_loadaer( self, istart, istop, k, m, cs, phase, naerosol, vaerosol, hygro)

    ! return aerosol number, volume concentrations, and bulk hygroscopicity

    class(modal_aerosol_model), intent(in) :: self
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

  subroutine modal_set_ptrs( self, raer, qqcw )
    class(modal_aerosol_model), intent(in) :: self
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: m,mm,l

   do m = 1, self%ntot_amode
      mm = self%indexer(m, 0)
      call rad_cnst_get_mode_num(0, m, 'a', self%state, self%pbuf, raer(mm)%fld)
      call rad_cnst_get_mode_num(0, m, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
      do l = 1, self%nspec_amode(m)
         mm = self%indexer(m, l)
         call rad_cnst_get_aer_mmr(0, m, l, 'a', self%state, self%pbuf, raer(mm)%fld)
         call rad_cnst_get_aer_mmr(0, m, l, 'c', self%state, self%pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
      end do
   end do
  end subroutine modal_set_ptrs

end module modal_aerosol_model_mod
