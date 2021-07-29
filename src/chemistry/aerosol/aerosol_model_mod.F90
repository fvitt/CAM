! Base Class for aerosol models

module aerosol_model_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physics_buffer, only: physics_buffer_desc
  use physics_types,  only: physics_state
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use ppgrid,         only: pcols, pver
  use ref_pres,       only: top_lev => trop_cloud_top_lev
  use physconst,      only: mwh2o, rhoh2o, r_universal, rh2o, pi
  use physconst,      only: gravit, latvap, cpair, rair

  implicit none

  type, abstract :: aerosol_model
     integer :: mtotal
     type(physics_state) :: state
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
     real(r8), allocatable :: amcubecoef(:)
     real(r8), allocatable :: argfactor(:)
     real(r8), allocatable :: alogsig(:) ! natl log of geometric standard dev of aerosol
     real(r8), allocatable :: f1(:) ! abdul-razzak functions of width
     real(r8), allocatable :: f2(:) ! abdul-razzak functions of width
   contains
     procedure :: create => aero_create
     procedure(aero_model_init), deferred :: model_init
     procedure(aero_model_final), deferred :: model_final
     procedure :: destroy => aero_destroy
     procedure(aero_loadaer), deferred :: loadaer
     procedure :: ccncalc => aero_ccncalc
     procedure :: maxsat => aero_maxsat
     procedure :: activate => aero_activate
     procedure :: err_funct => aero_err_funct
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

     subroutine aero_loadaer( self, istart, istop, k, m, cs, phase, naerosol, vaerosol, hygro)
       import

       class(aerosol_model), intent(in) :: self
       ! inputs
       integer,  intent(in) :: istart      ! start column index (1 <= istart <= istop <= pcols)
       integer,  intent(in) :: istop       ! stop column index
       integer,  intent(in) :: m           ! mode or bin index
       integer,  intent(in) :: k           ! level index
       real(r8), intent(in) :: cs(:,:)     ! air density (kg/m3)
       integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

       ! outputs
       real(r8), intent(out) :: naerosol(:)  ! number conc (1/m3)
       real(r8), intent(out) :: vaerosol(:)  ! volume conc (m3/m3)
       real(r8), intent(out) :: hygro(:)     ! bulk hygroscopicity of mode

     end subroutine aero_loadaer

  end interface

  integer,  parameter :: psat=6    ! number of supersaturations to calc ccn concentration
  real(r8), parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn concentration
                       (/ 0.02_r8, 0.05_r8, 0.1_r8, 0.2_r8, 0.5_r8, 1.0_r8 /)

  ! mathematical constants
  real(r8), parameter :: zero     = 0._r8
  real(r8), parameter :: third    = 1._r8/3._r8
  real(r8), parameter :: twothird = 2._r8*third
  real(r8), parameter :: sixth    = 1._r8/6._r8
  real(r8), parameter :: sq2      = sqrt(2._r8)
  real(r8), parameter :: sqpi     = sqrt(pi)
  real(r8), parameter :: t0       = 273._r8  ! reference temperature
  real(r8), parameter :: surften  = 0.076_r8 ! surface tension of water w/respect to air (N/m)
  real(r8), parameter :: alog2    = log(2._r8)
  real(r8), parameter :: alog3    = log(3._r8)

contains

  subroutine aero_create( self, state, pbuf )
    class(aerosol_model), intent(inout) :: self
    type(physics_state), intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    self%state = state
    self%pbuf => pbuf

    call self%model_init()

  end subroutine aero_create

  subroutine aero_destroy( self )
    class(aerosol_model), intent(inout) :: self

    nullify(self%pbuf)

    call self%model_final()

  end subroutine aero_destroy

  subroutine aero_ccncalc(self, state, cs, ccn)

    ! calculates number concentration of aerosols activated as CCN at
    ! supersaturation supersat.
    ! assumes an internal mixture of a multiple externally-mixed aerosol modes
    ! cgs units

    ! Ghan et al., Atmos. Res., 1993, 198-221.

    ! arguments
    class(aerosol_model), intent(in) :: self

    type(physics_state), target, intent(in)    :: state


    real(r8), intent(in)  :: cs(pcols,pver)       ! air density (kg/m3)
    real(r8), intent(out) :: ccn(pcols,pver,psat) ! number conc of aerosols activated at supersat (#/m3)

    ! local

    integer :: lchnk ! chunk index
    integer :: ncol  ! number of columns
    real(r8), pointer :: tair(:,:)     ! air temperature (K)

    real(r8) :: naerosol(pcols) ! interstit+activated aerosol number conc (/m3)
    real(r8) :: vaerosol(pcols) ! interstit+activated aerosol volume conc (m3/m3)

    real(r8) :: amcube(pcols)
    real(r8) :: a(pcols) ! surface tension parameter
    real(r8) :: hygro(pcols)  ! aerosol hygroscopicity
    real(r8) :: sm(pcols)  ! critical supersaturation at mode radius
    real(r8) :: arg(pcols)

    integer :: l,m,n,i,k
    real(r8) :: log,cc
    real(r8) :: smcoef(pcols)
    integer :: phase ! phase of aerosol

    real(r8), parameter :: smcoefcoef=2._r8/sqrt(27._r8)
    real(r8), parameter :: super(psat)=supersat(:psat)*0.01_r8
    real(r8) :: surften_coef

    !-------------------------------------------------------------------------------
    surften_coef = 2._r8*mwh2o*surften/(r_universal*rhoh2o)

    lchnk = state%lchnk
    ncol  = state%ncol
    tair  => state%t


    ccn = 0._r8
    do k=top_lev,pver

       do i=1,ncol
          a(i)=surften_coef/tair(i,k)
          smcoef(i)=smcoefcoef*a(i)*sqrt(a(i))
       end do

       do m=1,self%mtotal

          phase=3 ! interstitial+cloudborne

          call self%loadaer( 1, ncol, k, m, cs, phase, naerosol, vaerosol, hygro)

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

    use wv_saturation,only: qsat

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
    real(r8) rm ! number mode radius of aerosol at max supersat (cm)
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
    real(r8) xmincoeff
    real(r8) z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
    real(r8) etafactor1,etafactor2(nmode),etafactor2max
    real(r8) grow
    character(len=*), parameter :: subname='aero_model%activate'

    logical :: in_cloud
    integer m,n
    !      numerical integration parameters
    real(r8), parameter :: eps=0.3_r8,fmax=0.99_r8,sds=3._r8
    real(r8), parameter :: namin=1.e6_r8   ! minimum aerosol number concentration (/m3)
    real(r8) :: aten, alogaten

    integer ndist(nx)  ! accumulates frequency distribution of integration bins required
    data ndist/nx*0/
    save ndist

    aten     = 2._r8*mwh2o*surften/(r_universal*t0*rhoh2o)
    alogaten = log(aten)

    if (present(in_cloud_in)) then
       if (.not. present(smax_f)) call endrun(subname//' : smax_f must be supplied when in_cloud is used')
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

          amcube(m)= self%amcubecoef(m)*(volume(m)/na(m))
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
          xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3

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
             call endrun(subname)
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
             call endrun('activate')
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
          xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3


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


  function aero_err_funct(self,x) result(err)
    class(aerosol_model), intent(in) :: self
    real(r8), intent(in) :: x
    real(r8) :: err
    err = erf(x)
  end function aero_err_funct

 end module aerosol_model_mod
