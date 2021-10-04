! Base Class for aerosol models

module aerosol_model_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_spfn_mod,   only: erf => shr_spfn_erf
  use shr_const_mod,  only: pi => shr_const_pi
  use shr_const_mod,  only: r_universal => shr_const_rgas
  use shr_const_mod,  only: latvap => shr_const_latvap

  ! not parameter
  use physconst,      only: gravit, rair, cpair
  use physconst,      only: mwh2o, rhoh2o, rh2o

  use wv_saturation,  only: qsat ! for activate

  use physics_types,  only: physics_ptend, physics_ptend_init
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use cam_history,    only: addfld, add_default, horiz_only, fieldname_len, outfld

  use aerosol_data_mod,only: aerosol_data, ptr2d_t
  use cam_aerosol_data_mod, only: cam_aerosol_data

  implicit none

  private

  public :: aerosol_model

  type, abstract :: aerosol_model
     class(aerosol_data), pointer :: aero_data
     character(len=16) :: model_name = 'base'
     logical :: prognostic = .true.
     real(r8), allocatable :: amcubecoef(:)
     real(r8), allocatable :: exp45logsig(:)
     real(r8), allocatable :: argfactor(:)
     real(r8), allocatable :: alogsig(:) ! natl log of geometric standard dev of aerosol
     real(r8), allocatable :: f1(:) ! abdul-razzak functions of width
     real(r8), allocatable :: f2(:) ! abdul-razzak functions of width
   contains
     procedure :: create => aero_create
     procedure :: destroy => aero_destroy
     procedure :: dropmixnuc => aero_dropmixnuc
     procedure :: activate => aero_activate
     procedure(aero_model_init), deferred :: model_init
     procedure(aero_model_final), deferred :: model_final
     procedure :: err_funct => aero_err_funct
     procedure, private :: ccncalc => aero_ccncalc
     procedure, private :: maxsat => aero_maxsat
     procedure, private :: explmix => aero_explmix
  end type aerosol_model

  interface

     subroutine aero_model_init( self )
       import
       class(aerosol_model), intent(inout) :: self
     end subroutine aero_model_init

     subroutine aero_model_final( self )
       import
       class(aerosol_model), intent(inout) :: self
     end subroutine aero_model_final

  end interface

  real(r8), parameter :: zero     = 0._r8
  real(r8), parameter :: third    = 1._r8/3._r8
  real(r8), public, parameter :: twothird = 2._r8*third
  real(r8), public, parameter :: sq2      = sqrt(2._r8)
  real(r8), parameter :: sixth    = 1._r8/6._r8
  real(r8), parameter :: sqpi     = sqrt(pi)
  real(r8), parameter :: alog2    = log(2._r8)
  real(r8), parameter :: alog3    = log(3._r8)

  real(r8), parameter :: t0       = 273._r8  ! reference temperature
  real(r8), parameter :: surften  = 0.076_r8 ! surface tension of water w/respect to air (N/m)

  integer,  parameter :: psat=6    ! number of supersaturations to calc ccn concentration
  real(r8), parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn concentration
                       (/ 0.02_r8, 0.05_r8, 0.1_r8, 0.2_r8, 0.5_r8, 1.0_r8 /)

  character(len=4),parameter :: ccn_name(psat) = &
                    (/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'/)
  real(r8) :: aten, alogaten
  real(r8) :: surften_coef

contains

  !===============================================================================

  subroutine aero_create( self )
    class(aerosol_model), intent(inout) :: self

    call self%model_init()

    aten     = 2._r8*mwh2o*surften/(r_universal*t0*rhoh2o)
    alogaten = log(aten)
    surften_coef = 2._r8*mwh2o*surften/(r_universal*rhoh2o)

  end subroutine aero_create

  !===============================================================================

  subroutine aero_destroy( self )
    class(aerosol_model), intent(inout) :: self

    call self%model_final()

  end subroutine aero_destroy

  !===============================================================================

  subroutine aero_dropmixnuc( self, ncol, nlev, top_lev, lchnk, &
       ptend, dtmicro, wsub, &
       ncldwtr, temp, pmid, pint, pdel, rpdel, zm, kvh, &
       cldn, cldo, cldliqf, tendnd, factnum, &
       from_spcam, gasndx, gasnames, rgas )

    ! vertical diffusion and nucleation of cloud droplets
    ! assume cloud presence controlled by cloud fraction
    ! doesn't distinguish between warm, cold clouds

    ! arguments
    class(aerosol_model), intent(inout) :: self
    integer, intent(in) :: ncol, nlev, top_lev, lchnk
    real(r8), intent(in) :: ncldwtr(:,:) ! droplet number concentration (#/kg)
    real(r8), intent(in) :: temp(:,:)    ! temperature (K)
    real(r8), intent(in) :: pmid(:,:)    ! mid-level pressure (Pa)
    real(r8), intent(in) :: pint(:,:)    ! pressure at layer interfaces (Pa)
    real(r8), intent(in) :: pdel(:,:)    ! pressure thickess of layer (Pa)
    real(r8), intent(in) :: rpdel(:,:)   ! inverse of pressure thickess of layer (/Pa)
    real(r8), intent(in) :: zm(:,:)      ! geopotential height of level (m)
    real(r8), intent(in) :: kvh(:,:)     ! vertical diffusivity (m2/s)

    type(physics_ptend),         intent(out)   :: ptend
    real(r8),                    intent(in)    :: dtmicro     ! time step for microphysics (s)

    ! arguments
    real(r8), intent(in) :: wsub(:,:)    ! subgrid vertical velocity
    real(r8), intent(in) :: cldn(:,:)    ! cloud fraction
    real(r8), intent(in) :: cldo(:,:)    ! cloud fraction on previous time step
    real(r8), intent(in) :: cldliqf(:,:) ! liquid cloud fraction (liquid / (liquid + ice))

    logical,  intent(in),optional :: from_spcam ! value insignificant - if variable present, is called from spcam
    integer,  intent(in),optional :: gasndx(:)
    character(len=*), intent(in),optional :: gasnames(:)
    real(r8), intent(in),optional :: rgas(:,:,:)

    ! output arguments
    real(r8), intent(out) :: tendnd(:,:)    ! change in droplet number concentration (#/kg/s)
    real(r8), intent(out) :: factnum(:,:,:) ! activation fraction for aerosol number

    !--------------------Local storage-------------------------------------

    type(ptr2d_t), allocatable :: raer(:)     ! aerosol mass, number mixing ratios
    type(ptr2d_t), allocatable :: qqcw(:)
    real(r8) :: raertend(nlev)  ! tendency of aerosol mass, number mixing ratios
    real(r8) :: qqcwtend(nlev)  ! tendency of cloudborne aerosol mass, number mixing ratios


    real(r8), parameter :: zkmin = 0.01_r8, zkmax = 100._r8
    real(r8), parameter :: wmixmin = 0.1_r8        ! minimum turbulence vertical velocity (m/s)
    real(r8), parameter :: sq2pi = sqrt(2._r8*pi)

    integer  :: i, k, l, m, mm, n, ngas
    integer  :: km1, kp1
    integer  :: nnew, nsav, ntemp
    integer  :: lptr
    integer  :: nsubmix, nsubmix_bnd
    integer, save :: count_submix(100)
    integer  :: phase ! phase of aerosol

    real(r8) :: arg
    real(r8) :: dtinv
    real(r8) :: dtmin, tinv, dtt
    real(r8) :: lcldn(ncol,nlev)
    real(r8) :: lcldo(ncol,nlev)

    real(r8) :: zs(nlev) ! inverse of distance between levels (m)
    real(r8) :: qcld(nlev) ! cloud droplet number mixing ratio (#/kg)
    real(r8) :: qncld(nlev)     ! droplet number nucleated on cloud boundaries
    real(r8) :: srcn(nlev)       ! droplet source rate (/s)
    real(r8) :: cs(ncol,nlev)      ! air density (kg/m3)
    real(r8) :: csbot(nlev)       ! air density at bottom (interface) of layer (kg/m3)
    real(r8) :: csbot_cscen(nlev) ! csbot(i)/cs(i,k)
    real(r8) :: dz(ncol,nlev)      ! geometric thickness of layers (m)

    real(r8) :: wtke(ncol,nlev)     ! turbulent vertical velocity at base of layer k (m/s)
    real(r8) :: wtke_cen(ncol,nlev) ! turbulent vertical velocity at center of layer k (m/s)
    real(r8) :: wbar, wmix, wmin, wmax

    real(r8) :: zn(nlev)   ! g/pdel (m2/g) for layer
    real(r8) :: flxconv    ! convergence of flux into lowest layer

    real(r8) :: wdiab           ! diabatic vertical velocity
    real(r8) :: ekd(nlev)       ! diffusivity for droplets (m2/s)
    real(r8) :: ekk(0:nlev)     ! density*diffusivity for droplets (kg/m3 m2/s)
    real(r8) :: ekkp(nlev)      ! zn*zs*density*diffusivity
    real(r8) :: ekkm(nlev)      ! zn*zs*density*diffusivity

    real(r8) :: dum, dumc
    real(r8) :: tmpa
    real(r8) :: dact
    real(r8) :: fluxntot         ! (#/cm2/s)
    real(r8) :: dtmix
    real(r8) :: alogarg
    real(r8) :: overlapp(nlev), overlapm(nlev) ! cloud overlap

    real(r8) :: nsource(ncol,nlev)            ! droplet number source (#/kg/s)
    real(r8) :: ndropmix(ncol,nlev)           ! droplet number mixing (#/kg/s)
    real(r8) :: ndropcol(ncol)                ! column droplet number (#/m2)
    real(r8) :: cldo_tmp, cldn_tmp
    real(r8) :: tau_cld_regenerate
    real(r8) :: taumix_internal_nlev_inv ! 1/(internal mixing time scale for k=nlev) (1/s)


    real(r8), allocatable :: nact(:,:)  ! fractional aero. number  activation rate (/s)
    real(r8), allocatable :: mact(:,:)  ! fractional aero. mass    activation rate (/s)

    real(r8), allocatable :: raercol(:,:,:)    ! single column of aerosol mass, number mixing ratios
    real(r8), allocatable :: raercol_cw(:,:,:) ! same as raercol but for cloud-borne phase


    real(r8) :: na(ncol), va(ncol), hy(ncol)
    real(r8), allocatable :: naermod(:)  ! (1/m3)
    real(r8), allocatable :: hygro(:)    ! hygroscopicity of aerosol mode
    real(r8), allocatable :: vaerosol(:) ! interstit+activated aerosol volume conc (cm3/cm3)

    real(r8) :: source(nlev)

    real(r8), allocatable :: fn(:)              ! activation fraction for aerosol number
    real(r8), allocatable :: fm(:)              ! activation fraction for aerosol mass

    real(r8), allocatable :: fluxn(:)           ! number  activation fraction flux (cm/s)
    real(r8), allocatable :: fluxm(:)           ! mass    activation fraction flux (cm/s)
    real(r8)              :: flux_fullact(nlev) ! 100%    activation fraction flux (cm/s)
    !     note:  activation fraction fluxes are defined as
    !     fluxn = [flux of activated aero. number into cloud (#/cm2/s)]
    !           / [aero. number conc. in updraft, just below cloudbase (#/cm3)]


    real(r8), allocatable :: coltend(:,:)       ! column tendency for diagnostic output
    real(r8), allocatable :: coltend_cw(:,:)    ! column tendency
    real(r8) :: ccn(ncol,nlev,psat)    ! number conc of aerosols activated at supersat

    !for gas species turbulent mixing
    real(r8), allocatable :: rgascol(:,:,:)
    real(r8), allocatable :: coltendgas(:)
    real(r8) :: zerogas(nlev)
    character(len=200) fieldnamegas
    character(len=*), parameter :: subname='%dropmixnuc'

    logical  :: called_from_spcam
    !-------------------------------------------------------------------------------

    ! Create the liquid weighted cloud fractions that were passsed in
    ! before. This doesn't seem like the best variable, since the cloud could
    ! have liquid condensate, but the part of it that is changing could be the
    ! ice portion; however, this is what was done before.
    lcldo(:ncol,:) = cldo(:ncol,:) * cldliqf(:ncol,:)
    lcldn(:ncol,:) = cldn(:ncol,:) * cldliqf(:ncol,:)


    arg = 1.0_r8
    if (abs(0.8427_r8 - erf(arg))/0.8427_r8 > 0.001_r8) then
       write(iulog,*) 'erf(1.0) = ',ERF(arg)
       call endrun(trim(self%model_name)//subname//': Error function error')
    endif
    arg = 0.0_r8
    if (erf(arg) /= 0.0_r8) then
       write(iulog,*) 'erf(0.0) = ',erf(arg)
       write(iulog,*) 'dropmixnuc: Error function error'
       call endrun(trim(self%model_name)//subname//': Error function error')
    endif

    dtinv = 1._r8/dtmicro

    allocate( &
         nact(nlev,self%aero_data%mtotal),          &
         mact(nlev,self%aero_data%mtotal),          &
         raer(self%aero_data%ncnst_tot),                &
         qqcw(self%aero_data%ncnst_tot),                &
         raercol(nlev,self%aero_data%ncnst_tot,2),      &
         raercol_cw(nlev,self%aero_data%ncnst_tot,2),   &
         coltend(ncol,self%aero_data%ncnst_tot),       &
         coltend_cw(ncol,self%aero_data%ncnst_tot),    &
         naermod(self%aero_data%mtotal),            &
         hygro(self%aero_data%mtotal),              &
         vaerosol(self%aero_data%mtotal),           &
         fn(self%aero_data%mtotal),                 &
         fm(self%aero_data%mtotal),                 &
         fluxn(self%aero_data%mtotal),              &
         fluxm(self%aero_data%mtotal)               )

    call self%aero_data%set_ptrs( raer, qqcw )

    if (present(from_spcam).and.present(gasndx).and.&
        present(gasnames).and.present(rgas)) then
       called_from_spcam = from_spcam
       ngas = size(gasndx)
    else
       called_from_spcam = .false.
       ngas = 0
    endif

    if (called_from_spcam) then
       allocate(rgascol(nlev, ngas, 2))
       allocate(coltendgas(ncol))
    endif

    factnum = 0._r8
    wtke = 0._r8
    nsource = 0._r8
    ndropmix = 0._r8
    ndropcol = 0._r8

    select type (obj=>self%aero_data)
    class is (cam_aerosol_data)
       if (self%prognostic) then
          ! aerosol tendencies
          call physics_ptend_init(ptend, obj%state%psetcols, 'ndrop_'//trim(self%model_name), lq=obj%lq)
       else
          ! no aerosol tendencies
          call physics_ptend_init(ptend, obj%state%psetcols, 'ndrop_'//trim(self%model_name))
       end if
    end select

    ! overall_main_i_loop
    do i = 1, ncol

       do k = top_lev, nlev-1
          zs(k) = 1._r8/(zm(i,k) - zm(i,k+1))
       end do
       zs(nlev) = zs(nlev-1)

       ! load number nucleated into qcld on cloud boundaries

       do k = top_lev, nlev

          qcld(k)  = ncldwtr(i,k)
          qncld(k) = 0._r8
          srcn(k)  = 0._r8
          cs(i,k)  = pmid(i,k)/(rair*temp(i,k))        ! air density (kg/m3)
          dz(i,k)  = 1._r8/(cs(i,k)*gravit*rpdel(i,k)) ! layer thickness in m

          do m = 1, self%aero_data%mtotal
             nact(k,m) = 0._r8
             mact(k,m) = 0._r8
          end do

          zn(k) = gravit*rpdel(i,k)

          if (k < nlev) then
             ekd(k)   = kvh(i,k+1)
             ekd(k)   = max(ekd(k), zkmin)
             ekd(k)   = min(ekd(k), zkmax)
             csbot(k) = 2.0_r8*pint(i,k+1)/(rair*(temp(i,k) + temp(i,k+1)))
             csbot_cscen(k) = csbot(k)/cs(i,k)
          else
             ekd(k)   = 0._r8
             csbot(k) = cs(i,k)
             csbot_cscen(k) = 1.0_r8
          end if

          ! rce-comment - define wtke at layer centers for new-cloud activation
          !    and at layer boundaries for old-cloud activation
          wtke_cen(i,k) = wsub(i,k)
          wtke(i,k)     = wsub(i,k)
          wtke_cen(i,k) = max(wtke_cen(i,k), wmixmin)
          wtke(i,k)     = max(wtke(i,k), wmixmin)

          nsource(i,k) = 0._r8

       end do

       nsav = 1
       nnew = 2
       do mm = 1,self%aero_data%ncnst_tot
          raercol_cw(:,mm,nsav) = 0.0_r8
          raercol(:,mm,nsav)    = 0.0_r8
          raercol_cw(top_lev:nlev,mm,nsav) = qqcw(mm)%fld(i,top_lev:nlev)
          raercol(top_lev:nlev,mm,nsav)    = raer(mm)%fld(i,top_lev:nlev)
       end do

       if (called_from_spcam) then
          !
          ! In the MMF model, turbulent mixing for tracer species are turned off.
          ! So the turbulent for gas species mixing are added here.
          ! (Previously, it had the turbulent mixing for aerosol species)
          !
          do m=1, ngas
             rgascol(:,m,nsav) = rgas(i,:,m)
          end do

       endif

       ! droplet nucleation/aerosol activation

       ! tau_cld_regenerate = time scale for regeneration of cloudy air
       !    by (horizontal) exchange with clear air
       tau_cld_regenerate = 3600.0_r8 * 3.0_r8

       if (called_from_spcam) then
          ! when this is called  in the MMF part, no cloud regeneration and decay.
          ! set the time scale be very long so that no cloud regeneration.
          tau_cld_regenerate = 3600.0_r8 * 24.0_r8 * 365.0_r8
       endif


       ! k-loop for growing/shrinking cloud calcs .............................
       ! grow_shrink_main_k_loop: &
       do k = top_lev, nlev

          ! This code was designed for liquid clouds, but the cloudbourne
          ! aerosol can be either from liquid or ice clouds. For the ice clouds,
          ! we do not do regeneration, but as cloud fraction decreases the
          ! aerosols should be returned interstitial. The lack of a liquid cloud
          ! should not mean that all of the aerosol is realease. Therefor a
          ! section has been added for shrinking ice clouds and checks were added
          ! to protect ice cloudbourne aerosols from being released when no
          ! liquid cloud is present.

          ! shrinking ice cloud ......................................................
          cldo_tmp = cldo(i,k) * (1._r8 - cldliqf(i,k))
          cldn_tmp = cldn(i,k) * (1._r8 - cldliqf(i,k))

          if (cldn_tmp < cldo_tmp) then

             ! convert activated aerosol to interstitial in decaying cloud

             dumc = (cldn_tmp - cldo_tmp)/cldo_tmp * (1._r8 - cldliqf(i,k))
             do mm = 1,self%aero_data%ncnst_tot
                dact   = raercol_cw(k,mm,nsav)*dumc
                raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact   ! cloud-borne aerosol
                raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
             end do
          end if

          ! shrinking liquid cloud ......................................................
          !    treat the reduction of cloud fraction from when cldn(i,k) < cldo(i,k)
          !    and also dissipate the portion of the cloud that will be regenerated
          cldo_tmp = lcldo(i,k)
          cldn_tmp = lcldn(i,k) * exp( -dtmicro/tau_cld_regenerate )
          !    alternate formulation
          !    cldn_tmp = cldn(i,k) * max( 0.0_r8, (1.0_r8-dtmicro/tau_cld_regenerate) )

          ! fraction is also provided.
          if (cldn_tmp < cldo_tmp) then
             !  droplet loss in decaying cloud
             !++ sungsup
             nsource(i,k) = nsource(i,k) + qcld(k)*(cldn_tmp - cldo_tmp)/cldo_tmp*cldliqf(i,k)*dtinv
             qcld(k)      = qcld(k)*(1._r8 + (cldn_tmp - cldo_tmp)/cldo_tmp)
             !-- sungsup

             ! convert activated aerosol to interstitial in decaying cloud

             dumc = (cldn_tmp - cldo_tmp)/cldo_tmp * cldliqf(i,k)
             do mm = 1,self%aero_data%ncnst_tot
                dact   = raercol_cw(k,mm,nsav)*dumc
                raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact   ! cloud-borne aerosol
                raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
             end do
          end if

          ! growing liquid cloud ......................................................
          !    treat the increase of cloud fraction from when cldn(i,k) > cldo(i,k)
          !    and also regenerate part of the cloud
          cldo_tmp = cldn_tmp
          cldn_tmp = lcldn(i,k)

          if (cldn_tmp-cldo_tmp > 0.01_r8) then

             ! rce-comment - use wtke at layer centers for new-cloud activation
             wbar  = wtke_cen(i,k)
             wmix  = 0._r8
             wmin  = 0._r8
             wmax  = 10._r8
             wdiab = 0._r8

             ! load aerosol properties, assuming external mixtures

             phase = 1 ! interstitial
             do m = 1, self%aero_data%mtotal
                call self%aero_data%loadaer( i, i, k, m, cs, phase, na, va, hy)
                naermod(m)  = na(i)
                vaerosol(m) = va(i)
                hygro(m)    = hy(i)
             end do

             call self%activate( &
                  wbar, wmix, wdiab, wmin, wmax,                       &
                  temp(i,k), cs(i,k), naermod, self%aero_data%mtotal,             &
                  vaerosol, hygro, fn, fm, fluxn,     &
                  fluxm,flux_fullact(k))

             factnum(i,k,:) = fn

             dumc = (cldn_tmp - cldo_tmp)
             do m = 1, self%aero_data%mtotal
                mm = self%aero_data%indexer(m,0)
                dact   = dumc*fn(m)*raer(mm)%fld(i,k) ! interstitial only
                qcld(k) = qcld(k) + dact
                nsource(i,k) = nsource(i,k) + dact*dtinv
                raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
                raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
                dum = dumc*fm(m)
                do l = 1, self%aero_data%nmasses(m)
                   mm = self%aero_data%indexer(m,l)
                   dact    = dum*raer(mm)%fld(i,k) ! interstitial only
                   raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
                   raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
                enddo
             enddo
          endif

       enddo  ! grow_shrink_main_k_loop
       ! end of k-loop for growing/shrinking cloud calcs ......................

       ! ......................................................................
       ! start of k-loop for calc of old cloud activation tendencies ..........
       !
       ! rce-comment
       !    changed this part of code to use current cloud fraction (cldn) exclusively
       !    consider case of cldo(:)=0, cldn(k)=1, cldn(k+1)=0
       !    previous code (which used cldo below here) would have no cloud-base activation
       !       into layer k.  however, activated particles in k mix out to k+1,
       !       so they are incorrectly depleted with no replacement

       ! old_cloud_main_k_loop
       do k = top_lev, nlev
          kp1 = min0(k+1, nlev)
          taumix_internal_nlev_inv = 0.0_r8

          if (lcldn(i,k) > 0.01_r8) then

             wdiab = 0._r8
             wmix  = 0._r8                       ! single updraft
             wbar  = wtke(i,k)                   ! single updraft
             if (k == nlev) wbar = wtke_cen(i,k) ! single updraft
             wmax  = 10._r8
             wmin  = 0._r8

             if (lcldn(i,k) - lcldn(i,kp1) > 0.01_r8 .or. k == nlev) then

                ! cloud base

                ! ekd(k) = wtke(i,k)*dz(i,k)/sq2pi
                ! rce-comments
                !   first, should probably have 1/zs(k) here rather than dz(i,k) because
                !      the turbulent flux is proportional to ekd(k)*zs(k),
                !      while the dz(i,k) is used to get flux divergences
                !      and mixing ratio tendency/change
                !   second and more importantly, using a single updraft velocity here
                !      means having monodisperse turbulent updraft and downdrafts.
                !      The sq2pi factor assumes a normal draft spectrum.
                !      The fluxn/fluxm from activate must be consistent with the
                !      fluxes calculated in explmix.
                ekd(k) = wbar/zs(k)

                alogarg = max(1.e-20_r8, 1._r8/lcldn(i,k) - 1._r8)
                wmin    = wbar + wmix*0.25_r8*sq2pi*log(alogarg)
                phase   = 1   ! interstitial

                do m = 1, self%aero_data%mtotal
                   ! rce-comment - use kp1 here as old-cloud activation involves
                   !   aerosol from layer below
                   call self%aero_data%loadaer( i, i, kp1, m, cs, phase, na, va, hy)
                   naermod(m)  = na(i)
                   vaerosol(m) = va(i)
                   hygro(m)    = hy(i)
                end do

                call self%activate( &
                     wbar, wmix, wdiab, wmin, wmax,                       &
                     temp(i,k), cs(i,k), naermod, self%aero_data%mtotal,             &
                     vaerosol, hygro,  fn, fm, fluxn,     &
                     fluxm, flux_fullact(k))

                factnum(i,k,:) = fn

                if (k < nlev) then
                   dumc = lcldn(i,k) - lcldn(i,kp1)
                else
                   dumc = lcldn(i,k)
                endif

                fluxntot = 0

                ! rce-comment 1
                !    flux of activated mass into layer k (in kg/m2/s)
                !       = "actmassflux" = dumc*fluxm*raercol(kp1,lmass)*csbot(k)
                !    source of activated mass (in kg/kg/s) = flux divergence
                !       = actmassflux/(cs(i,k)*dz(i,k))
                !    so need factor of csbot_cscen = csbot(k)/cs(i,k)
                !                   dum=1./(dz(i,k))
                dum=csbot_cscen(k)/(dz(i,k))

                ! rce-comment 2
                !    code for k=pver was changed to use the following conceptual model
                !    in k=pver, there can be no cloud-base activation unless one considers
                !       a scenario such as the layer being partially cloudy,
                !       with clear air at bottom and cloudy air at top
                !    assume this scenario, and that the clear/cloudy portions mix with
                !       a timescale taumix_internal = dz(i,pver)/wtke_cen(i,pver)
                !    in the absence of other sources/sinks, qact (the activated particle
                !       mixratio) attains a steady state value given by
                !          qact_ss = fcloud*fact*qtot
                !       where fcloud is cloud fraction, fact is activation fraction,
                !       qtot=qact+qint, qint is interstitial particle mixratio
                !    the activation rate (from mixing within the layer) can now be
                !       written as
                !          d(qact)/dt = (qact_ss - qact)/taumix_internal
                !                     = qtot*(fcloud*fact*wtke/dz) - qact*(wtke/dz)
                !    note that (fcloud*fact*wtke/dz) is equal to the nact/mact
                !    also, d(qact)/dt can be negative.  in the code below
                !       it is forced to be >= 0
                !
                ! steve --
                !    you will likely want to change this.  i did not really understand
                !       what was previously being done in k=pver
                !    in the cam3_5_3 code, wtke(i,pver) appears to be equal to the
                !       droplet deposition velocity which is quite small
                !    in the cam3_5_37 version, wtke is done differently and is much
                !       larger in k=pver, so the activation is stronger there
                !
                if (k == nlev) then
                   taumix_internal_nlev_inv = flux_fullact(k)/dz(i,k)
                end if

                do m = 1, self%aero_data%mtotal
                   mm = self%aero_data%indexer(m,0)
                   fluxn(m) = fluxn(m)*dumc
                   fluxm(m) = fluxm(m)*dumc
                   nact(k,m) = nact(k,m) + fluxn(m)*dum
                   mact(k,m) = mact(k,m) + fluxm(m)*dum
                   if (k < nlev) then
                      ! note that kp1 is used here
                      fluxntot = fluxntot &
                           + fluxn(m)*raercol(kp1,mm,nsav)*cs(i,k)
                   else
                      tmpa = raercol(kp1,mm,nsav)*fluxn(m) &
                           + raercol_cw(kp1,mm,nsav)*(fluxn(m) &
                           - taumix_internal_nlev_inv*dz(i,k))
                      fluxntot = fluxntot + max(0.0_r8, tmpa)*cs(i,k)
                   end if
                end do
                srcn(k)      = srcn(k) + fluxntot/(cs(i,k)*dz(i,k))
                nsource(i,k) = nsource(i,k) + fluxntot/(cs(i,k)*dz(i,k))

             endif  ! (cldn(i,k) - cldn(i,kp1) > 0.01 .or. k == nlev)

          else

             ! no liquid cloud
             nsource(i,k) = nsource(i,k) - qcld(k)*dtinv
             qcld(k)      = 0

             if (cldn(i,k) < 0.01_r8) then
                ! no ice cloud either

                ! convert activated aerosol to interstitial in decaying cloud

                do mm = 1,self%aero_data%ncnst_tot
                   raercol(k,mm,nsav)    = raercol(k,mm,nsav) + raercol_cw(k,mm,nsav)  ! cloud-borne aerosol
                   raercol_cw(k,mm,nsav) = 0._r8
                end do

             end if
          end if

       end do  ! old_cloud_main_k_loop

       ! switch nsav, nnew so that nnew is the updated aerosol
       ntemp = nsav
       nsav  = nnew
       nnew  = ntemp

       ! load new droplets in layers above, below clouds

       dtmin     = dtmicro
       ekk(top_lev-1)    = 0.0_r8
       ekk(nlev) = 0.0_r8
       do k = top_lev, nlev-1
          ! rce-comment -- ekd(k) is eddy-diffusivity at k/k+1 interface
          !   want ekk(k) = ekd(k) * (density at k/k+1 interface)
          !   so use pint(i,k+1) as pint is 1:nlevp
          !           ekk(k)=ekd(k)*2.*pint(i,k)/(rair*(temp(i,k)+temp(i,k+1)))
          !           ekk(k)=ekd(k)*2.*pint(i,k+1)/(rair*(temp(i,k)+temp(i,k+1)))
          ekk(k) = ekd(k)*csbot(k)
       end do

       do k = top_lev, nlev
          km1     = max0(k-1, top_lev)
          ekkp(k) = zn(k)*ekk(k)*zs(k)
          ekkm(k) = zn(k)*ekk(k-1)*zs(km1)
          tinv    = ekkp(k) + ekkm(k)

          ! rce-comment -- tinv is the sum of all first-order-loss-rates
          !    for the layer.  for most layers, the activation loss rate
          !    (for interstitial particles) is accounted for by the loss by
          !    turb-transfer to the layer above.
          !    k=nlev is special, and the loss rate for activation within
          !    the layer must be added to tinv.  if not, the time step
          !    can be too big, and explmix can produce negative values.
          !    the negative values are reset to zero, resulting in an
          !    artificial source.
          if (k == nlev) tinv = tinv + taumix_internal_nlev_inv

          if (tinv .gt. 1.e-6_r8) then
             dtt   = 1._r8/tinv
             dtmin = min(dtmin, dtt)
          end if
       end do

       dtmix   = 0.9_r8*dtmin
       nsubmix = dtmicro/dtmix + 1
       if (nsubmix > 100) then
          nsubmix_bnd = 100
       else
          nsubmix_bnd = nsubmix
       end if
       count_submix(nsubmix_bnd) = count_submix(nsubmix_bnd) + 1
       dtmix = dtmicro/nsubmix

       do k = top_lev, nlev
          kp1 = min(k+1, nlev)
          km1 = max(k-1, top_lev)
          ! maximum overlap assumption
          if (cldn(i,kp1) > 1.e-10_r8) then
             overlapp(k) = min(cldn(i,k)/cldn(i,kp1), 1._r8)
          else
             overlapp(k) = 1._r8
          end if
          if (cldn(i,km1) > 1.e-10_r8) then
             overlapm(k) = min(cldn(i,k)/cldn(i,km1), 1._r8)
          else
             overlapm(k) = 1._r8
          end if
       end do


       ! rce-comment
       !    the activation source(k) = mact(k,m)*raercol(kp1,lmass)
       !       should not exceed the rate of transfer of unactivated particles
       !       from kp1 to k which = ekkp(k)*raercol(kp1,lmass)
       !    however it might if things are not "just right" in subr activate
       !    the following is a safety measure to avoid negatives in explmix
       do k = top_lev, nlev-1
          do m = 1, self%aero_data%mtotal
             nact(k,m) = min( nact(k,m), ekkp(k) )
             mact(k,m) = min( mact(k,m), ekkp(k) )
          end do
       end do


       ! old_cloud_nsubmix_loop
       do n = 1, nsubmix
          qncld(:) = qcld(:)
          ! switch nsav, nnew so that nsav is the updated aerosol
          ntemp   = nsav
          nsav    = nnew
          nnew    = ntemp
          srcn(:) = 0.0_r8

          do m = 1, self%aero_data%mtotal
             mm = self%aero_data%indexer(m,0)

             ! update droplet source
             ! rce-comment- activation source in layer k involves particles from k+1
             !	       srcn(:)=srcn(:)+nact(:,m)*(raercol(:,mm,nsav))
             srcn(top_lev:nlev-1) = srcn(top_lev:nlev-1) + nact(top_lev:nlev-1,m)*(raercol(top_lev+1:nlev,mm,nsav))

             ! rce-comment- new formulation for k=pver
             !              srcn(  pver  )=srcn(  pver  )+nact(  pver  ,m)*(raercol(  pver,mm,nsav))
             tmpa = raercol(nlev,mm,nsav)*nact(nlev,m) &
                  + raercol_cw(nlev,mm,nsav)*(nact(nlev,m) - taumix_internal_nlev_inv)
             srcn(nlev) = srcn(nlev) + max(0.0_r8,tmpa)
          end do
          call self%explmix(  &
               qcld, srcn, ekkp, ekkm, overlapp,  &
               overlapm, qncld, zero, zero, nlev,top_lev, &
               dtmix, .false.)

          ! rce-comment
          !    the interstitial particle mixratio is different in clear/cloudy portions
          !    of a layer, and generally higher in the clear portion.  (we have/had
          !    a method for diagnosing the the clear/cloudy mixratios.)  the activation
          !    source terms involve clear air (from below) moving into cloudy air (above).
          !    in theory, the clear-portion mixratio should be used when calculating
          !    source terms
          do m = 1, self%aero_data%mtotal
             mm = self%aero_data%indexer(m,0)
             ! rce-comment -   activation source in layer k involves particles from k+1
             !	              source(:)= nact(:,m)*(raercol(:,mm,nsav))
             source(top_lev:nlev-1) = nact(top_lev:nlev-1,m)*(raercol(top_lev+1:nlev,mm,nsav))
             ! rce-comment - new formulation for k=pver
             !               source(  pver  )= nact(  pver,  m)*(raercol(  pver,mm,nsav))
             tmpa = raercol(nlev,mm,nsav)*nact(nlev,m) &
                  + raercol_cw(nlev,mm,nsav)*(nact(nlev,m) - taumix_internal_nlev_inv)
             source(nlev) = max(0.0_r8, tmpa)
             flxconv = 0._r8

             call self%explmix( &
                  raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
                  overlapm, raercol_cw(:,mm,nsav), zero, zero, nlev,top_lev,   &
                  dtmix, .false.)

             call self%explmix( &
                  raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
                  overlapm, raercol(:,mm,nsav), zero, flxconv, nlev,top_lev, &
                  dtmix, .true., raercol_cw(:,mm,nsav))

             do l = 1, self%aero_data%nmasses(m)
                mm = self%aero_data%indexer(m,l)
                ! rce-comment -   activation source in layer k involves particles from k+1
                !	          source(:)= mact(:,m)*(raercol(:,mm,nsav))
                source(top_lev:nlev-1) = mact(top_lev:nlev-1,m)*(raercol(top_lev+1:nlev,mm,nsav))
                ! rce-comment- new formulation for k=pver
                !                 source(  pver  )= mact(  pver  ,m)*(raercol(  pver,mm,nsav))
                tmpa = raercol(nlev,mm,nsav)*mact(nlev,m) &
                     + raercol_cw(nlev,mm,nsav)*(mact(nlev,m) - taumix_internal_nlev_inv)
                source(nlev) = max(0.0_r8, tmpa)
                flxconv = 0._r8

                call self%explmix( &
                     raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
                     overlapm, raercol_cw(:,mm,nsav), zero, zero, nlev,top_lev,   &
                     dtmix, .false.)

                call self%explmix( &
                     raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
                     overlapm, raercol(:,mm,nsav), zero, flxconv, nlev,top_lev, &
                     dtmix, .true., raercol_cw(:,mm,nsav))

             end do
          end do

          if (called_from_spcam) then
             !
             ! turbulent mixing for gas species .
             !
             do m=1, ngas
                flxconv = 0.0_r8
                zerogas(:) = 0.0_r8
                call self%explmix(rgascol(:,m,nnew),zerogas,ekkp,ekkm,overlapp,overlapm, &
                                  rgascol(:,m,nsav),zero, flxconv, nlev,top_lev, dtmix, &
                                  .true., zerogas)

             end do
          endif

       end do ! old_cloud_nsubmix_loop

       ! evaporate particles again if no cloud (either ice or liquid)

       do k = top_lev, nlev
          if (cldn(i,k) == 0._r8) then
             ! no ice or liquid cloud
             qcld(k)=0._r8

             ! convert activated aerosol to interstitial in decaying cloud
             do mm = 1,self%aero_data%ncnst_tot
                raercol(k,mm,nnew)    = raercol(k,mm,nnew) + raercol_cw(k,mm,nnew)
                raercol_cw(k,mm,nnew) = 0._r8
             end do

          end if
       end do

       ! droplet number

       ndropcol(i) = 0._r8
       do k = top_lev, nlev
          ndropmix(i,k) = (qcld(k) - ncldwtr(i,k))*dtinv - nsource(i,k)
          tendnd(i,k)   = (max(qcld(k), 1.e-6_r8) - ncldwtr(i,k))*dtinv
          ndropcol(i)   = ndropcol(i) + ncldwtr(i,k)*pdel(i,k)
       end do
       ndropcol(i) = ndropcol(i)/gravit

       if (self%prognostic) then

          call self%aero_data%update( pdel, raer, qqcw, raercol(:,:,nnew), raercol_cw(:,:,nnew), i,dtinv, &
                                      coltend(i,:),coltend_cw(i,:) )
       end if

       if (called_from_spcam) then
          !
          ! Gas tendency
          !
          do m=1, ngas
             mm = gasndx(m)
             ptend%lq(mm) = .true.
             ptend%q(i, :, mm) = (rgascol(:,m,nnew)-rgas(i,:,m)) * dtinv
          end do
       endif

    end do  ! overall_main_i_loop
    ! end of main loop over i/longitude ....................................

    call outfld('NDROPCOL', ndropcol, ncol, lchnk)
    call outfld('NDROPSRC', nsource,  ncol, lchnk)
    call outfld('NDROPMIX', ndropmix, ncol, lchnk)
    call outfld('WTKE    ', wtke,     ncol, lchnk)

    if(called_from_spcam) then
       call outfld('SPLCLOUD  ', cldn(:ncol,:), ncol, lchnk   )
       call outfld('SPKVH     ', kvh(:ncol,:), ncol, lchnk   )
    endif

    call self%ccncalc(ncol,nlev,top_lev, temp, cs, ccn)
    do l = 1, psat
       call outfld(ccn_name(l), ccn(:,:,l), ncol, lchnk)
    enddo

    ! do column tendencies
    if (self%prognostic) then
       do mm = 1,self%aero_data%ncnst_tot
          call outfld(self%aero_data%fieldname(mm),    coltend(:,mm),    ncol, lchnk)
          call outfld(self%aero_data%fieldname_cw(mm), coltend_cw(:,mm), ncol, lchnk)
       end do
    end if

    if(called_from_spcam) then
       !
       ! output column-integrated Gas tendency (this should be zero)
       !
       do m=1, ngas
          mm = gasndx(m)
          do i=1, ncol
             coltendgas(i) = sum( pdel(i,:)*ptend%q(i,:,mm) )/gravit
          end do
          fieldnamegas = trim(gasnames(m)) // '_mixnuc1sp'
          call outfld( trim(fieldnamegas), coltendgas, ncol, lchnk)
       end do
       deallocate(rgascol, coltendgas)
    end if

    deallocate( &
         nact,       &
         mact,       &
         raer,       &
         qqcw,       &
         raercol,    &
         raercol_cw, &
         coltend,    &
         coltend_cw, &
         naermod,    &
         hygro,      &
         vaerosol,   &
         fn,         &
         fm,         &
         fluxn,      &
         fluxm       )

  end subroutine aero_dropmixnuc

  !===============================================================================

  subroutine aero_ccncalc(self, ncol,nlev,toplev, tair, cs, ccn)

    ! calculates number concentration of aerosols activated as CCN at
    ! supersaturation supersat.
    ! assumes an internal mixture of a multiple externally-mixed aerosol modes
    ! cgs units

    ! Ghan et al., Atmos. Res., 1993, 198-221.

    ! arguments
    class(aerosol_model), intent(in) :: self
    integer,  intent(in)  :: ncol,nlev,toplev
    real(r8), intent(in)  :: tair(:,:)  ! air temperature (K)
    real(r8), intent(in)  :: cs(:,:)    ! air density (kg/m3)
    real(r8), intent(out) :: ccn(:,:,:) ! number conc of aerosols activated at supersat (#/m3)

    ! local

    real(r8) :: naerosol(ncol) ! interstit+activated aerosol number conc (/m3)
    real(r8) :: vaerosol(ncol) ! interstit+activated aerosol volume conc (m3/m3)

    real(r8) :: amcube(ncol)
    real(r8) :: a(ncol) ! surface tension parameter
    real(r8) :: hygro(ncol)  ! aerosol hygroscopicity
    real(r8) :: sm(ncol)  ! critical supersaturation at mode radius
    real(r8) :: arg(ncol)

    integer :: l,m,i,k
    real(r8) :: smcoef(ncol)
    integer :: phase ! phase of aerosol

    real(r8), parameter :: smcoefcoef=2._r8/sqrt(27._r8)
    real(r8), parameter :: super(psat)=supersat(:psat)*0.01_r8

    !-------------------------------------------------------------------------------

    ccn = 0._r8
    do k=toplev,nlev

       do i=1,ncol
          a(i)=surften_coef/tair(i,k)
          smcoef(i)=smcoefcoef*a(i)*sqrt(a(i))
       end do

       do m=1,self%aero_data%mtotal

          phase=3 ! interstitial+cloudborne

          call self%aero_data%loadaer( 1, ncol, k, m, cs, phase, naerosol, vaerosol, hygro)

          where(naerosol(:ncol)>1.e-3_r8 .and. hygro(:ncol)>0._r8)
             amcube(:ncol)=self%amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
             sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol)) ! critical supersaturation
          elsewhere
             sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
          endwhere
          do l=1,psat
             do i=1,ncol
                arg(i)=self%argfactor(m)*log(sm(i)/super(l))
                ccn(i,k,l)=ccn(i,k,l)+naerosol(i)*0.5_r8*(1._r8-erf(arg(i)))
             enddo
          enddo
       enddo
    enddo
    ccn(:ncol,:,:)=ccn(:ncol,:,:)*1.e-6_r8 ! convert from #/m3 to #/cm3

  end subroutine aero_ccncalc

  !===============================================================================

  subroutine aero_maxsat(self,zeta,eta,nmode,smc,smax)

    !      calculates maximum supersaturation for multiple
    !      competing aerosol modes.

    !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

    class(aerosol_model), intent(in) :: self
    integer,  intent(in)  :: nmode ! number of modes
    real(r8), intent(in)  :: smc(nmode) ! critical supersaturation for number mode radius
    real(r8), intent(in)  :: zeta(nmode)
    real(r8), intent(in)  :: eta(nmode)
    real(r8), intent(out) :: smax ! maximum supersaturation
    integer  :: m  ! mode index
    real(r8) :: sum, g1, g2, g1sqrt, g2sqrt

    do m=1,nmode
       if(zeta(m).gt.1.e5_r8*eta(m).or.smc(m)*smc(m).gt.1.e5_r8*eta(m))then
          ! weak forcing. essentially none activated
          smax=1.e-20_r8
       else
          ! significant activation of this mode. calc activation all modes.
          exit
       endif
       ! No significant activation in any mode.  Do nothing.
       if (m == nmode) return

    enddo

    sum=0.0_r8
    do m=1,nmode
       if(eta(m).gt.1.e-20_r8)then
          g1=zeta(m)/eta(m)
          g1sqrt=sqrt(g1)
          g1=g1sqrt*g1
          g2=smc(m)/sqrt(eta(m)+3._r8*zeta(m))
          g2sqrt=sqrt(g2)
          g2=g2sqrt*g2
          sum=sum+(self%f1(m)*g1+self%f2(m)*g2)/(smc(m)*smc(m))
       else
          sum=1.e20_r8
       endif
    enddo

    smax=1._r8/sqrt(sum)

  end subroutine aero_maxsat

  !===============================================================================

  subroutine aero_activate(self, wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,  &
       na, nmode, volume, hygro,  &
       fn, fm, fluxn, fluxm, flux_fullact, smax_prescribed, in_cloud_in, smax_f)

    !      calculates number, surface, and mass fraction of aerosols activated as CCN
    !      calculates flux of cloud droplets, surface area, and aerosol mass into cloud
    !      assumes an internal mixture within each of up to nbin multiple aerosol bins
    !      a gaussiam spectrum of updrafts can be treated. mode -> bins

    !      mks units

    !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.


    !      input
    class(aerosol_model), intent(in) :: self

    real(r8), intent(in) :: wbar          ! grid cell mean vertical velocity (m/s)
    real(r8), intent(in) :: sigw          ! subgrid standard deviation of vertical vel (m/s)
    real(r8), intent(in) :: wdiab         ! diabatic vertical velocity (0 if adiabatic)
    real(r8), intent(in) :: wminf         ! minimum updraft velocity for integration (m/s)
    real(r8), intent(in) :: wmaxf         ! maximum updraft velocity for integration (m/s)
    real(r8), intent(in) :: tair          ! air temperature (K)
    real(r8), intent(in) :: rhoair        ! air density (kg/m3)
    real(r8), intent(in) :: na(:)      ! aerosol number concentration (/m3)
    integer,  intent(in) :: nmode      ! number of aerosol modes or bins
    real(r8), intent(in) :: volume(:)  ! aerosol volume concentration (m3/m3)
    real(r8), intent(in) :: hygro(:)   ! hygroscopicity of aerosol mode

    !      output

    real(r8), intent(out) :: fn(:)      ! number fraction of aerosols activated
    real(r8), intent(out) :: fm(:)      ! mass fraction of aerosols activated
    real(r8), intent(out) :: fluxn(:)   ! flux of activated aerosol number fraction into cloud (cm/s)
    real(r8), intent(out) :: fluxm(:)   ! flux of activated aerosol mass fraction into cloud (cm/s)
    real(r8), intent(out) :: flux_fullact   ! flux of activated aerosol fraction assuming 100% activation (cm/s)
    !    rce-comment
    !    used for consistency check -- this should match (ekd(k)*zs(k))
    !    also, fluxm/flux_fullact gives fraction of aerosol mass flux
    !       that is activated

    !      optional
    real(r8), optional, intent(in) :: smax_prescribed  ! prescribed max. supersaturation for secondary activation
    logical,  optional, intent(in) :: in_cloud_in      ! switch to modify calculations when above cloud base
    real(r8), optional, intent(in) :: smax_f           ! droplet and rain size distr factor in the smax calculation
    ! used when in_cloud=.true.

    !      local

    integer, parameter:: nx=200
    real(r8) integ,integf
    real(r8), parameter :: p0 = 1013.25e2_r8    ! reference pressure (Pa)
    real(r8) pres ! pressure (Pa)
    real(r8) diff0,conduct0
    real(r8) es ! saturation vapor pressure
    real(r8) qs ! water vapor saturation mixing ratio
    real(r8) dqsdt ! change in qs with temperature
    real(r8) g ! thermodynamic function (m2/s)
    real(r8) zeta(nmode), eta(nmode)
    real(r8) lnsmax ! ln(smax)
    real(r8) alpha
    real(r8) gamma
    real(r8) beta
    real(r8) sqrtg
    real(r8) :: amcube(nmode) ! cube of dry mode radius (m)
    real(r8) :: lnsm(nmode)
    real(r8) smc(nmode) ! critical supersaturation for number mode radius
    real(r8) sumflx_fullact
    real(r8) sumflxn(nmode)
    real(r8) sumflxm(nmode)
    real(r8) sumfn(nmode)
    real(r8) sumfm(nmode)
    real(r8) fnold(nmode)   ! number fraction activated
    real(r8) fmold(nmode)   ! mass fraction activated
    real(r8) wold,gold
    real(r8) wmin,wmax,w,dw,dwmax,dwmin,wnuc,dwnew,wb
    real(r8) dfmin,dfmax,fnew,fold,fnmin,fnbar,fmbar
    real(r8) alw,sqrtalw
    real(r8) smax
    real(r8) x,arg
    real(r8) z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
    real(r8) etafactor1,etafactor2(nmode),etafactor2max
    real(r8) grow
    character(len=*), parameter :: subname='%activate'

    logical :: in_cloud
    integer m,n
    !      numerical integration parameters
    real(r8), parameter :: eps=0.3_r8,fmax=0.99_r8,sds=3._r8
    real(r8), parameter :: namin=1.e6_r8   ! minimum aerosol number concentration (/m3)

    integer ndist(nx)  ! accumulates frequency distribution of integration bins required
    data ndist/nx*0/
    save ndist

    if (present(in_cloud_in)) then
       if (.not. present(smax_f)) then
          call endrun(trim(self%model_name)//subname &
                      //' : smax_f must be supplied when in_cloud is used')
       end if
       in_cloud = in_cloud_in
    else
       in_cloud = .false.
    end if

    fn(:)=0._r8
    fm(:)=0._r8
    fluxn(:)=0._r8
    fluxm(:)=0._r8
    flux_fullact=0._r8

    if(nmode.eq.1.and.na(1).lt.1.e-20_r8)return

    if(sigw.le.1.e-5_r8.and.wbar.le.0._r8)return

    if ( present( smax_prescribed ) ) then
       if (smax_prescribed <= 0.0_r8) return
    end if


    pres=rair*rhoair*tair
    diff0=0.211e-4_r8*(p0/pres)*(tair/t0)**1.94_r8
    conduct0=(5.69_r8+0.017_r8*(tair-t0))*4.186e2_r8*1.e-5_r8 ! convert to J/m/s/deg
    call qsat(tair, pres, es, qs)
    dqsdt=latvap/(rh2o*tair*tair)*qs
    alpha=gravit*(latvap/(cpair*rh2o*tair*tair)-1._r8/(rair*tair))
    gamma=(1.0_r8+latvap/cpair*dqsdt)/(rhoair*qs)
    etafactor2max=1.e10_r8/(alpha*wmaxf)**1.5_r8 ! this should make eta big if na is very small.

    grow  = 1._r8/(rhoh2o/(diff0*rhoair*qs)  &
         + latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair) - 1._r8))
    sqrtg = sqrt(grow)
    beta  = 2._r8*pi*rhoh2o*grow*gamma

    do m=1,nmode

       if(volume(m).gt.1.e-39_r8.and.na(m).gt.1.e-39_r8)then

          amcube(m)=(3._r8*volume(m)/(4._r8*pi*self%exp45logsig(m)*na(m)))
          !           growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
          !           should depend on mean radius of mode to account for gas kinetic effects
          !           see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
          !           for approriate size to use for effective diffusivity.
          etafactor2(m)=1._r8/(na(m)*beta*sqrtg)
          if(hygro(m).gt.1.e-10_r8)then
             smc(m)=2._r8*aten*sqrt(aten/(27._r8*hygro(m)*amcube(m))) ! only if variable size dist
          else
             smc(m)=100._r8
          endif
       else
          smc(m)=1._r8
          etafactor2(m)=etafactor2max ! this should make eta big if na is very small.
       endif
       lnsm(m)=log(smc(m)) ! only if variable size dist
    enddo

    if(sigw.gt.1.e-5_r8)then ! spectrum of updrafts

       wmax=min(wmaxf,wbar+sds*sigw)
       wmin=max(wminf,-wdiab)
       wmin=max(wmin,wbar-sds*sigw)
       w=wmin
       dwmax=eps*sigw
       dw=dwmax
       dfmax=0.2_r8
       dfmin=0.1_r8
       if (wmax <= w) return
       do m=1,nmode
          sumflxn(m)=0._r8
          sumfn(m)=0._r8
          fnold(m)=0._r8
          sumflxm(m)=0._r8
          sumfm(m)=0._r8
          fmold(m)=0._r8
       enddo
       sumflx_fullact=0._r8

       fold=0._r8
       wold=0._r8
       gold=0._r8

       dwmin = min( dwmax, 0.01_r8 )
       do n = 1, nx

100       wnuc=w+wdiab
          !           write(iulog,*)'wnuc=',wnuc
          alw=alpha*wnuc
          sqrtalw=sqrt(alw)
          etafactor1=alw*sqrtalw

          do m=1,nmode
             eta(m)=etafactor1*etafactor2(m)
             zeta(m)=twothird*sqrtalw*aten/sqrtg
          enddo

          if ( present( smax_prescribed ) ) then
             smax = smax_prescribed
          else
             call self%maxsat(zeta,eta,nmode,smc,smax)
          endif
          !	      write(iulog,*)'w,smax=',w,smax

          lnsmax=log(smax)

          x=twothird*(lnsm(nmode)-lnsmax)/(sq2*self%alogsig(nmode))
          fnew=0.5_r8*(1._r8-self%err_funct(x))

          dwnew = dw
          if(fnew-fold.gt.dfmax.and.n.gt.1)then
             !              reduce updraft increment for greater accuracy in integration
             if (dw .gt. 1.01_r8*dwmin) then
                dw=0.7_r8*dw
                dw=max(dw,dwmin)
                w=wold+dw
                go to 100
             else
                dwnew = dwmin
             endif
          endif

          if(fnew-fold.lt.dfmin)then
             !              increase updraft increment to accelerate integration
             dwnew=min(1.5_r8*dw,dwmax)
          endif
          fold=fnew

          z=(w-wbar)/(sigw*sq2)
          g=exp(-z*z)
          fnmin=1._r8

          do m=1,nmode
             !              modal
             x=twothird*(lnsm(m)-lnsmax)/(sq2*self%alogsig(m))
             fn(m)=0.5_r8*(1._r8-self%err_funct(x))
             fnmin=min(fn(m),fnmin)
             !               integration is second order accurate
             !               assumes linear variation of f*g with w
             fnbar=(fn(m)*g+fnold(m)*gold)
             arg=x-1.5_r8*sq2*self%alogsig(m)
             fm(m)=0.5_r8*(1._r8-self%err_funct(arg))
             fmbar=(fm(m)*g+fmold(m)*gold)
             wb=(w+wold)
             if(w.gt.0._r8)then
                sumflxn(m)=sumflxn(m)+sixth*(wb*fnbar           &
                     +(fn(m)*g*w+fnold(m)*gold*wold))*dw
                sumflxm(m)=sumflxm(m)+sixth*(wb*fmbar           &
                     +(fm(m)*g*w+fmold(m)*gold*wold))*dw
             endif
             sumfn(m)=sumfn(m)+0.5_r8*fnbar*dw
             fnold(m)=fn(m)
             sumfm(m)=sumfm(m)+0.5_r8*fmbar*dw
             fmold(m)=fm(m)
          enddo
          !           same form as sumflxm but replace the fm with 1.0
          sumflx_fullact = sumflx_fullact &
               + sixth*(wb*(g+gold) + (g*w+gold*wold))*dw
          !            sumg=sumg+0.5_r8*(g+gold)*dw
          gold=g
          wold=w
          dw=dwnew
          if (n > 1 .and. (w > wmax .or. fnmin > fmax)) exit
          w=w+dw
          if (n == nx) then
             write(iulog,*)'do loop is too short in activate'
             write(iulog,*)'wmin=',wmin,' w=',w,' wmax=',wmax,' dw=',dw
             write(iulog,*)'wbar=',wbar,' sigw=',sigw,' wdiab=',wdiab
             write(iulog,*)'wnuc=',wnuc
             write(iulog,*)'na=',(na(m),m=1,nmode)
             write(iulog,*)'fn=',(fn(m),m=1,nmode)
             !   dump all subr parameters to allow testing with standalone code
             !   (build a driver that will read input and call activate)
             write(iulog,*)'wbar,sigw,wdiab,tair,rhoair,nmode='
             write(iulog,*) wbar,sigw,wdiab,tair,rhoair,nmode
             write(iulog,*)'na=',na
             write(iulog,*)'volume=', (volume(m),m=1,nmode)
             write(iulog,*)'hydro='
             write(iulog,*) hygro
             call endrun(trim(self%model_name)//subname)
          end if

       enddo

       ndist(n)=ndist(n)+1
       if(w.lt.wmaxf)then

          !            contribution from all updrafts stronger than wmax
          !            assuming constant f (close to fmax)
          wnuc=w+wdiab

          z1=(w-wbar)/(sigw*sq2)
          z2=(wmaxf-wbar)/(sigw*sq2)
          g=exp(-z1*z1)
          integ=sigw*0.5_r8*sq2*sqpi*(erf(z2)-erf(z1))
          !            consider only upward flow into cloud base when estimating flux
          wf1=max(w,zero)
          zf1=(wf1-wbar)/(sigw*sq2)
          gf1=exp(-zf1*zf1)
          wf2=max(wmaxf,zero)
          zf2=(wf2-wbar)/(sigw*sq2)
          gf2=exp(-zf2*zf2)
          gf=(gf1-gf2)
          integf=wbar*sigw*0.5_r8*sq2*sqpi*(erf(zf2)-erf(zf1))+sigw*sigw*gf

          do m=1,nmode
             sumflxn(m)=sumflxn(m)+integf*fn(m)
             sumfn(m)=sumfn(m)+fn(m)*integ
             sumflxm(m)=sumflxm(m)+integf*fm(m)
             sumfm(m)=sumfm(m)+fm(m)*integ
          enddo
          !           same form as sumflxm but replace the fm with 1.0
          sumflx_fullact = sumflx_fullact + integf
          !            sumg=sumg+integ
       endif


       do m=1,nmode
          fn(m)=sumfn(m)/(sq2*sqpi*sigw)
          !            fn(m)=sumfn(m)/(sumg)
          if(fn(m).gt.1.01_r8)then
             write(iulog,*)'fn=',fn(m),' > 1 in activate'
             write(iulog,*)'w,m,na,amcube=',w,m,na(m),amcube(m)
             write(iulog,*)'integ,sumfn,sigw=',integ,sumfn(m),sigw
             call endrun(trim(self%model_name)//subname)
          endif
          fluxn(m)=sumflxn(m)/(sq2*sqpi*sigw)
          fm(m)=sumfm(m)/(sq2*sqpi*sigw)
          !            fm(m)=sumfm(m)/(sumg)
          if(fm(m).gt.1.01_r8)then
             write(iulog,*)'fm=',fm(m),' > 1 in activate'
          endif
          fluxm(m)=sumflxm(m)/(sq2*sqpi*sigw)
       enddo
       !        same form as fluxm
       flux_fullact = sumflx_fullact/(sq2*sqpi*sigw)

    else

       !        single updraft
       wnuc=wbar+wdiab

       if(wnuc.gt.0._r8)then

          w=wbar

          if(in_cloud) then

             if (smax_f > 0._r8) then
                smax = alpha*w/(2.0_r8*pi*rhoh2o*grow*gamma*smax_f)
             else
                smax = 1.e-20_r8
             end if

          else ! at cloud base
             alw        = alpha*wnuc
             sqrtalw    = sqrt(alw)
             etafactor1 = alw*sqrtalw

             do m = 1, nmode
                eta(m)  = etafactor1*etafactor2(m)
                zeta(m) = twothird*sqrtalw*aten/sqrtg
             end do
             if ( present(smax_prescribed) ) then
                smax = smax_prescribed
             else
                call self%maxsat(zeta, eta, nmode, smc, smax)
             end if
          end if

          lnsmax=log(smax)

          do m=1,nmode
             !                 modal
             x=twothird*(lnsm(m)-lnsmax)/(sq2*self%alogsig(m))
             fn(m)=0.5_r8*(1._r8-self%err_funct(x))
             arg=x-1.5_r8*sq2*self%alogsig(m)
             fm(m)=0.5_r8*(1._r8-self%err_funct(arg))

             if(wbar.gt.0._r8)then
                fluxn(m)=fn(m)*w
                fluxm(m)=fm(m)*w
             endif
          enddo
          flux_fullact = w
       endif

    endif

  end subroutine aero_activate

  !===============================================================================

  subroutine aero_explmix( self, q, src, ekkp, ekkm, overlapp, overlapm, &
       qold, surfrate, flxconv, nlev, top_lev, dt, is_unact, qactold )

    !  explicit integration of droplet/aerosol mixing
    !     with source due to activation/nucleation


    class(aerosol_model), intent(in) :: self
    integer, intent(in) :: nlev ! number of levels
    integer, intent(in) :: top_lev !
    real(r8), intent(out) :: q(nlev) ! mixing ratio to be updated
    real(r8), intent(in) :: qold(nlev) ! mixing ratio from previous time step
    real(r8), intent(in) :: src(nlev) ! source due to activation/nucleation (/s)
    real(r8), intent(in) :: ekkp(nlev) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
    ! below layer k  (k,k+1 interface)
    real(r8), intent(in) :: ekkm(nlev) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
    ! above layer k  (k,k+1 interface)
    real(r8), intent(in) :: overlapp(nlev) ! cloud overlap below
    real(r8), intent(in) :: overlapm(nlev) ! cloud overlap above
    real(r8), intent(in) :: surfrate ! surface exchange rate (/s)
    real(r8), intent(in) :: flxconv ! convergence of flux from surface
    real(r8), intent(in) :: dt ! time step (s)
    logical, intent(in) :: is_unact ! true if this is an unactivated species
    real(r8), intent(in),optional :: qactold(nlev)
    ! mixing ratio of ACTIVATED species from previous step
    ! *** this should only be present
    !     if the current species is unactivated number/sfc/mass

    integer k,kp1,km1

    if ( is_unact ) then
       !     the qactold*(1-overlap) terms are resuspension of activated material
       do k=top_lev,nlev
          kp1=min(k+1,nlev)
          km1=max(k-1,top_lev)
          q(k) = qold(k) + dt*( - src(k) + ekkp(k)*(qold(kp1) - qold(k) +       &
               qactold(kp1)*(1.0_r8-overlapp(k)))               &
               + ekkm(k)*(qold(km1) - qold(k) +     &
               qactold(km1)*(1.0_r8-overlapm(k))) )
          !        force to non-negative
          !        if(q(k)<-1.e-30)then
          !           write(iulog,*)'q=',q(k),' in explmix'
          q(k)=max(q(k),0._r8)
          !        endif
       end do

       !     diffusion loss at base of lowest layer
       q(nlev)=q(nlev)-surfrate*qold(nlev)*dt+flxconv*dt
       !        force to non-negative
       !        if(q(nlev)<-1.e-30)then
       !           write(iulog,*)'q=',q(nlev),' in explmix'
       q(nlev)=max(q(nlev),0._r8)
       !        endif
    else
       do k=top_lev,nlev
          kp1=min(k+1,nlev)
          km1=max(k-1,top_lev)
          q(k) = qold(k) + dt*(src(k) + ekkp(k)*(overlapp(k)*qold(kp1)-qold(k)) +      &
               ekkm(k)*(overlapm(k)*qold(km1)-qold(k)) )
          !        force to non-negative
          !        if(q(k)<-1.e-30)then
          !           write(iulog,*)'q=',q(k),' in explmix'
          q(k)=max(q(k),0._r8)
          !        endif
       end do
       !     diffusion loss at base of lowest layer
       q(nlev)=q(nlev)-surfrate*qold(nlev)*dt+flxconv*dt
       !        force to non-negative
       !        if(q(nlev)<-1.e-30)then
       !           write(iulog,*)'q=',q(nlev),' in explmix'
       q(nlev)=max(q(nlev),0._r8)

    end if

  end subroutine aero_explmix

  !===============================================================================

  function aero_err_funct(self,x) result(err)
    class(aerosol_model), intent(in) :: self
    real(r8), intent(in) :: x
    real(r8) :: err
    err = erf(x)
  end function aero_err_funct

 end module aerosol_model_mod
