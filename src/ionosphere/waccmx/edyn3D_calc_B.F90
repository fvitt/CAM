!---------------------------------------------------------------------- 
    module edyn3D_delB_module
! 
! nglon,nglat,glon,glat refer to grid for calculating delB.
! nggjlon,nggjlat,ggjlon,ggjlat refer to the grid for the currents,
!  which may be different from that for delB.
      use edyn3D_params,    only: nlat_qd,nmlon,h0,rtd,dtr,ylonm,nglon,nglat,glon,glat, &
                                  k_fix_ggjbot,nggjlon,nggjlat,nggjhgt,ggjhgt,ggjtop,delBsolution,&
				  h_LEO,mu0,l1ggj,l2ggj,lvhags,lshags,lwork,mdab
      use edyn3D_qd_module, only: qd,Keast,Ksouth,Jrtopgg,etagg
      use shr_kind_mod,     only: r8 => shr_kind_r8            ! 8-byte reals
      use shr_const_mod,    only: re => SHR_CONST_REARTH ! meters
      use physconst,        only: pi
! 
      implicit none
! 
      real(r8), dimension(mdab,nggjlat,nggjhgt) :: ar,ai,br,bi,cr,ci
      real(r8) :: wvhags(lvhags)
      real(r8) :: wshags(lshags)
      real(r8) :: work(lwork)

      real(r8) :: ggjbot(nggjhgt+1)
! SH coefficients for magnetic potential divided by radius (current layer base)
      real(r8) :: vorc(mdab,nggjlat,nggjhgt+1),vors(mdab,nggjlat,nggjhgt+1)
! SH coefficients for upward field within ionosphere (current layer base):
      real(r8) :: betac(mdab,nggjlat,nggjhgt+1),betas(mdab,nggjlat,nggjhgt+1)
! SH coefficients for external ground magnetic potential divided by radius
      real(r8) :: vextorc(mdab,nggjlat),vextors(mdab,nggjlat)
! SH coefficients for external upward field at ground:
      real(r8) :: betaextc(mdab,nggjlat),betaexts(mdab,nggjlat)
! SH coefficients for equivalent current function:
      real(r8) :: psicoefc(mdab,nggjlat),psicoefs(mdab,nggjlat)
      
      real(r8) :: bhat_g(nglon,nglat,nggjhgt+1,3), &   ! unit vector (east, north, up) along magnetc fieldline
           bhat_gLEO(nglon,nglat,3), &             ! unit vector (east, north, up) along magnetc fieldline
           qdlon_g(nglon,nglat),qdlat_g(nglon,nglat)
      
      real(r8),dimension(nglon,nglat) :: delbegrd, delbngrd ,delbugrd,psi
      real(r8),dimension(nglon,nglat) :: delbeLEO,delbnLEO,delbuLEO,delbsclLEO
      real(r8),dimension(nglon,nglat,nggjhgt+1) :: delbe, delbn ,delbu,delbscl

      real(r8),dimension(nmlon,nlat_qd-1) :: delbegrd_qd, delbngrd_qd,delbugrd_qd,psi_qd
      real(r8),dimension(nmlon,nlat_qd-1) :: delbeLEO_qd,delbnLEO_qd,delbuLEO_qd,delbscl_qd
      real(r8),dimension(nlat_qd-1,nggjhgt+1) :: delbe0ln_qd,delbn0ln_qd,delbu0ln_qd,delscllH_qd

     contains
!---------------------------------------------------------------------- 
      subroutine gen_geo_grid
!     
      implicit none
 !     
      integer :: ilateq,i
      real(r8) :: dlon,dlat
!      
! Set up geographic grid for outputting magnetic perturbations.
      dlon = 360./float(nglon-1)
      dlat = 180./float(nglat-1)
      ilateq = (nglat+1)/2
      do i=1,nglon
        glon(i) = 0.+(i-1)*dlon
      enddo
      do i=1,nglat
        glat(i) = (i-ilateq)*dlat
      enddo     
!      
      end subroutine gen_geo_grid
!---------------------------------------------------------------------- 
      subroutine calc_Bcoef
!     
! Uses SPHEREPACK 3.2 subroutines.
! shags is used to get SH coefficients of equivalent current function
!   etagg, and SH coefficients of radial current density.
! vhags is used to get SH coefficients of toroidal current. 
!     
      implicit none
      integer,parameter :: ldwork=max0(nggjlat*(nggjlat+4), &
                                (nggjlat*(3*nggjlat+9)+2)/2)
      double precision :: dwork(ldwork)
      integer :: ierror
      integer :: k,kk,mp,m,np,n
      real(r8) :: fac,fac2,fac3,bradius,rob,robnm1,bornp2,xnorm
!     
      ggjbot(1) = h0
      do k=1,nggjhgt
        ggjbot(k+1) = ggjtop(k)
      enddo
!
! Initialize whags and wvhags
!     subroutine shagsi(nlat,nlon,wshags,lshags,work,lwork,dwork,ldwork,
!    +                  ierror)
      call shagsi(nggjlat,nggjlon,wshags,lshags,work,lwork,dwork, &
                  ldwork,ierror)
      if (ierror.ne.0) then
        write (6,*) 'shagsi ierror= ',ierror
        stop
      endif
!      subroutine vhagsi(nlat,nlon,wvhags,lvhags,dwork,ldwork,ierror)
      call vhagsi(nggjlat,nggjlon,wvhags,lvhags,dwork,ldwork,ierror)
      if (ierror.ne.0) then
        write (6,*) 'vhagsi ierror= ',ierror
        stop
      endif
!     
! Compute delB coefficients for all layer bottoms and ground
!     
! First compute external coefficients associated with FAC
!     subroutine shags(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
!    1                    wshags,lshags,work,lwork,ierror)
      call shags(nggjlat,nggjlon,0,1,etagg,nggjlat,nggjlon,ar,ai &
                 ,mdab,nggjlat,wshags,lshags,work,lwork,ierror)
      if (ierror.ne.0) then
        write (6,*) 'shags (1) ierror= ',ierror
        stop
      endif
      bradius = re + ggjtop(nggjhgt)
! Contribution of FAC to external coefficients (divided by Re) and
!  upward field coefficients at Earth's surface
      rob = re/bradius
      robnm1 = 1./rob
      do np=2,nggjlat
        n = np - 1
        robnm1 = robnm1*rob
        fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
        fac  = n*fac2
        do mp=1,min0(np,mdab)
          vextorc(mp,np) = fac2*ar(mp,np,1)
          vextors(mp,np) = fac2*ai(mp,np,1)
          betaextc(mp,np) = -fac*ar(mp,np,1)
          betaexts(mp,np) = -fac*ai(mp,np,1)
        enddo ! mp
      enddo ! np
      if (delBsolution.ne.'quick_ground    ') then
! Contribution of FAC to external coefficients (divided by r) and
!  upward field coefficients at bottom of each current layer.
      do kk=1,nggjhgt+1
        rob = (re + ggjbot(kk))/bradius
        robnm1 = 1./rob
        do np=2,nggjlat
          n = np - 1
          robnm1 = robnm1*rob
          fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
          fac  = n*fac2
          do mp=1,min0(np,mdab)
            vorc(mp,np,kk) = fac2*ar(mp,np,1)
            vors(mp,np,kk) = fac2*ai(mp,np,1)
            betac(mp,np,kk) = -fac*ar(mp,np,1)
            betas(mp,np,kk) = -fac*ai(mp,np,1)
          enddo ! mp
        enddo ! np
      enddo ! kk
! Contribution of layer currents to external coefficients (divided by r) and
!  upward field coefficients at bottom of each current layer.
!      call vhags(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
!     +                 mdab,ndab,wvhags,lvhags,work,lwork,ierror)
! ityp=2 assumes horizontal divergence of current is zero, calculates only cr,ci
      call vhags(nggjlat,nggjlon,0,nggjhgt,Ksouth,Keast,nggjlat,nggjlon &
       ,br,bi,cr,ci,mdab,nggjlat,wvhags,lvhags,work,lwork,ierror)
      if (ierror.ne.0) then
        write (6,*) 'vhags ierror= ',ierror
        stop
      endif
! Rescale cr,ci to account for the facts that vhags uses a different
!  normalization of spherical harmonics for vectors than for scalars,
!  and that vhags uses a unit radius.
      do np=1,nggjlat
        xnorm = amax1(sqrt(float((np-1)*np)),1.)
        do kk=1,nggjhgt
          fac = (re+ggjhgt(kk))/xnorm
          do mp=1,mdab
            cr(mp,np,kk) = cr(mp,np,kk)*fac
            ci(mp,np,kk) = ci(mp,np,kk)*fac
          enddo
        enddo
      enddo

      do k=1,nggjhgt ! loop over all current layers
        bradius = re + ggjhgt(k)
! Calculate coefficients for layer bottoms lying above the kth current layer
        do kk=k+1,nggjhgt+1
          rob = (re + ggjbot(kk))/bradius
          bornp2 = 1./(rob*rob)
          do np=2,nggjlat
            n = np - 1
            bornp2 = bornp2/rob
            fac2 = mu0*bornp2*float(n)/(float(2*n+1)*bradius)
            fac  = (n+1)*fac2
            do mp=1,min0(np,mdab)
              vorc(mp,np,kk) = vorc(mp,np,kk) + fac2*cr(mp,np,k)
              vors(mp,np,kk) = vors(mp,np,kk) + fac2*ci(mp,np,k)
              betac(mp,np,kk) = betac(mp,np,kk) + fac*cr(mp,np,k)
              betas(mp,np,kk) = betas(mp,np,kk) + fac*ci(mp,np,k)
            enddo ! mp
          enddo ! np
        enddo ! kk
! Calculate coefficients for layer bottoms lying below the kth current layer
        do kk=1,k
          rob = (re + ggjbot(kk))/bradius
          robnm1 = 1./rob
          do np=2,nggjlat
            n = np - 1
            robnm1 = robnm1*rob
            fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
            fac  = n*fac2
            do mp=1,min0(np,mdab)
              vorc(mp,np,kk) = vorc(mp,np,kk) - fac2*cr(mp,np,k)
              vors(mp,np,kk) = vors(mp,np,kk) - fac2*ci(mp,np,k)
              betac(mp,np,kk) = betac(mp,np,kk) + fac*cr(mp,np,k)
              betas(mp,np,kk) = betas(mp,np,kk) + fac*ci(mp,np,k)
            enddo ! mp
          enddo ! np
        enddo ! kk
! Calculate coefficients for ground.
        rob = re/bradius
        robnm1 = 1./rob
        do np=2,nggjlat
          n = np - 1
          robnm1 = robnm1*rob
          fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
          fac  = n*fac2
          do mp=1,min0(np,mdab)
            vextorc(mp,np) = vextorc(mp,np) - fac2*cr(mp,np,k)
            vextors(mp,np) = vextors(mp,np) - fac2*ci(mp,np,k)
            betaextc(mp,np) = betaextc(mp,np) + fac*cr(mp,np,k)
            betaexts(mp,np) = betaexts(mp,np) + fac*ci(mp,np,k)
          enddo ! mp
        enddo ! np
      enddo ! k
! Compute coefficients of radial current density at tops of layers 
!     subroutine shags(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
!    1                    wshags,lshags,work,lwork,ierror)
      call shags(nggjlat,nggjlon,0,nggjhgt,Jrtopgg,nggjlat,nggjlon,ar,ai &
                 ,mdab,nggjlat,wshags,lshags,work,lwork,ierror)
      if (ierror.ne.0) then
        write (6,*) 'shags (2) ierror= ',ierror
        stop
      endif
! Convert ar,ai to be coefficients for tau/r at tops of layers
      do np=2,nggjlat
        n = np - 1
        do mp=1,min0(np,mdab)
          do k=1,nggjhgt
            fac3 = mu0*(re+ggjtop(k))/float(n*(n+1)) 
            ar(mp,np,k) = fac3*ar(mp,np,k) 
            ai(mp,np,k) = fac3*ai(mp,np,k) 
          enddo ! k
        enddo ! mp
      enddo ! np
      endif ! (delBsolution.ne.'quick_ground    ') 
!      
! Calculate equivalent current coefficients for 110 km current layer.
      bradius = re + 110.e3
      rob = re/bradius
      robnm1 = 1./rob
      do np=2,nggjlat
        n = np - 1
        robnm1 = robnm1*rob
        fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
        fac  = n*fac2
        do mp=1,min0(np,mdab)
          psicoefc(mp,np) = vextorc(mp,np)/fac2
          psicoefs(mp,np) = vextors(mp,np)/fac2
        enddo ! mp
      enddo ! np
!      
      end subroutine calc_Bcoef
!---------------------------------------------------------------------- 
      subroutine calc_B
! Examples that calculates the internal field using perfect conductor at
!  depth=600 km, and that calculates total perturbation fields over the
!  Earth's surface and at h_LEO, and globally at all nggjhgt+1 heights
!  of ggjbot.
! Uses SPHEREPACK 3.2 subroutines.  Note that "nlon" input to most
!  subroutines (except math2geo routines) excludes wraparound points,
!  and therefore equals nglon-1. nlon is assumed to be an even number.
!  nglat is assumed to be an odd number.
! Following assumes nglon includes 1 duplicated (wraparound) point.
      integer,parameter :: l1=min0(nglat,(nglon/2))
      integer,parameter :: l1s=min0(nglat,((nglon-1)/2)+1)
      integer,parameter :: l2=(nglat+1)/2
      integer,parameter :: l3=(l1*(2*nglat-l1+1))/2
      integer,parameter :: lvhses=2*l2*l3+(nglon-1)+15
      integer,parameter :: lshses=(l1s*l2*(2*nglat-l1s+1))/2+(nglon-1)+15
      integer,parameter :: &
        lworkvhsesi=3*(max0(l1s-2,0)*(nglat+nglat-l1s-1))/2+5*l2*nglat, &
        lworkgrades=2*(nggjhgt+1)*nglat*(nglon-1+l1)+nglat*nglon
      integer,parameter :: ldworkvhsesi=2*(nglat+1)
!      double precision :: dwork(ldworkvhsesi)
      real(r8) :: dwork(ldworkvhsesi)
      real(r8) :: wvhses(lvhses),work2(lworkgrades),wshses(lshses)
      real(r8), parameter :: depth = 6.e5  ! meters
      real(r8),dimension(mdab,nggjlat) :: &
        vintorc,betaintc, &
        vintors,betaints, &
        vortotc,betatotc, &
        vortots,betatots, &
        torc,tors
      real(r8),dimension(nglat,nglon) :: bsouth,beast,bup,gradts,gradte
      real(r8),dimension(nglat,nglon,nggjhgt+1) :: &
         bsouth3D,beast3D,bup3D,gradts3D,gradte3D
!      
      integer :: n,mp,np,lon,lat,k,kbelow,ierror,j,i
      real(r8) :: frac,fac,coa,coa2,coa2np1,aob,aobnp1,arbelow,aibelow
      real(r8) :: pertrb(nggjhgt)

! Find internal coefficients
      coa = (re-depth)/re
      coa2 = coa**2
      coa2np1 = coa
      do np=2,nggjlat
        n = np - 1
        coa2np1 = coa2np1*coa2
        fac = coa2np1*float(n)/float(n+1)
        do mp=1,min0(np,mdab)
          vintorc(mp,np)  =  fac*vextorc(mp,np)
          vintors(mp,np)  =  fac*vextors(mp,np)
          betaintc(mp,np) =  -coa2np1*betaextc(mp,np)
          betaints(mp,np) =  -coa2np1*betaexts(mp,np)
          vortotc(mp,np)  = vextorc(mp,np) + vintorc(mp,np) 
          vortots(mp,np)  = vextors(mp,np) + vintors(mp,np) 
          betatotc(mp,np) = betaextc(mp,np) + betaintc(mp,np) 
          betatots(mp,np) = betaexts(mp,np) + betaints(mp,np) 
        enddo ! mp
      enddo ! np
!
! Find k index for nearest coefficients below h_LEO
      do k=nggjhgt,1,-1
        if (ggjbot(k).lt.h_LEO) then
          kbelow = k
          exit 
        endif
      enddo ! k
! frac is fractional distance of h_LEO from ggjbot(kbelow+1) down towards
!   ggjbot(kbelow).
      frac = (ggjbot(kbelow+1) - h_LEO)/(ggjbot(kbelow+1)-ggjbot(kbelow))
      if (frac.lt.0.or.frac.gt.1.) then
        write (6,*) 'Stopped because frac =', frac
        stop
      endif
!
! Initialize wvhses and wshses
!     subroutine vhsesi(nlat,nlon,wvhses,lvhses,work,lwork,dwork,
!    +                  ldwork,ierror)
      call vhsesi(nglat,nglon-1,wvhses,lvhses,work,lworkvhsesi,dwork,ldworkvhsesi,ierror)
      if (ierror.ne.0) then
        write (6,*) 'vhsesi ierror= ',ierror
        stop
      endif
!     subroutine shsesi(nlat,nlon,wshses,lshses,work,lwork,dwork,
!    +                  ldwork,ierror)
      call shsesi(nglat,nglon-1,wshses,lshses,work,lworkvhsesi,dwork,ldworkvhsesi,ierror)
      if (ierror.ne.0) then
        write (6,*) 'shsesi ierror= ',ierror
        stop
      endif
!
! Calculate global fields at ground and h_LEO
!     subroutine grades(nlat,nlon,isym,nt,v,w,idvw,jdvw,a,b,mdab,ndab,
!    +                  wvhses,lvhses,work,lwork,ierror)
      call grades(nglat,nglon-1,0,1,bsouth,beast,nglat,nglon,vortotc,vortots &
         ,mdab,nggjlat,wvhses,lvhses,work2,lworkgrades,ierror)
      if (ierror.ne.0) then
        write (6,*) '(1) grades ierror= ',ierror
        stop
      endif
! Fill wrap-around points in longitude
      do j=1,nglat
        bsouth(j,nglon) = bsouth(j,1)
        beast(j,nglon) = beast(j,1)
      enddo
! Since grades calculates the positive gradient, whereas the magnetic
!   field is the negative gradient, reverse the signs of bsouth,beast
      bsouth = -bsouth
      beast = -beast
!     (4) subroutine math2geov(ig,nlat,nlon,vm,wm,ug,vg,work)
      call math2geov(0,nglat,nglon,bsouth,beast,delbegrd,delbngrd,work2)
      if (ierror.ne.0) then
        write (6,*) 'math2geov ierror= ',ierror
        stop
      endif

!     subroutine shses(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
!    +                 wshses,lshses,work,lwork,ierror)
      call shses(nglat,nglon-1,0,1,bup,nglat,nglon,betatotc,betatots,mdab,nggjlat, &
                      wshses,lshses,work2,lworkgrades,ierror)
      if (ierror.ne.0) then
        write (6,*) 'shses ierror= ',ierror
        stop
      endif

! Fill wrap-around points in longitude
      do j=1,nglat
        bup(j,nglon) = bup(j,1)
      enddo
!     (2) subroutine math2geos(ig,nlat,nlon,sm,sg,work)
      call math2geos(0,nglat,nglon,bup,delbugrd,work2)
      if (ierror.ne.0) then
        write (6,*) 'math2geos ierror= ',ierror
        stop
      endif
!
      aob = re/(re+h_LEO)
      aobnp1 = aob
      do np=2,nggjlat
        aobnp1 = aobnp1*aob
        n = np - 1
        do mp=1,min0(np,mdab)
          vortotc(mp,np) = vintorc(mp,np)*aobnp1 + frac*vorc(mp,np,kbelow) &
                + (1.-frac)*vorc(mp,np,kbelow+1)
          vortots(mp,np) = vintors(mp,np)*aobnp1 + frac*vors(mp,np,kbelow) &
                + (1.-frac)*vors(mp,np,kbelow+1)
          betatotc(mp,np) = betaintc(mp,np)*aobnp1*aob + frac*betac(mp,np,kbelow) &
                + (1.-frac)*betac(mp,np,kbelow+1)
          betatots(mp,np) = betaints(mp,np)*aobnp1*aob + frac*betas(mp,np,kbelow) &
                + (1.-frac)*betas(mp,np,kbelow+1)
! Whereas vortot and betatot refer to the bottom of the current layer,
!   ar,ai refer to the top of the current layer, so we shift the
!   height index by 1, and assume zero vertical current at the bottom
!   of the lowest layer:
          arbelow = 0.
          aibelow = 0.
          if (kbelow.gt.1) then
            arbelow = ar(mp,np,kbelow-1)
            aibelow = ai(mp,np,kbelow-1)
          endif
          torc(mp,np) = frac*arbelow + (1.-frac)*ar(mp,np,kbelow)
          tors(mp,np) = frac*aibelow + (1.-frac)*ai(mp,np,kbelow)
        enddo ! mp
      enddo ! np
!     subroutine grades(nlat,nlon,isym,nt,v,w,idvw,jdvw,a,b,mdab,ndab,
!    +                  wvhses,lvhses,work,lwork,ierror)
      call grades(nglat,nglon-1,0,1,bsouth,beast,nglat,nglon,vortotc,vortots &
         ,mdab,nggjlat,wvhses,lvhses,work2,lworkgrades,ierror)
      if (ierror.ne.0) then
        write (6,*) '(2) grades ierror= ',ierror
        stop
      endif

! Get toroidal magnetic field (divided by r*mu0)

!! Test
!      write(6,*) 'Reset torc, tors'
!      torc = 0.
!      tors = 0.
!      torc(2,3) = 1.
!      do np=2,nggjlat
!        n = np - 1
!        do mp=2,min0(np,mdab)
!!          torc(mp,np) = ar(mp,np,8)
!!          tors(mp,np) = ai(mp,np,8)
!          write(6,'(2i5,2e11.3)') n,mp-1,torc(mp,np),tors(mp,np)
!        enddo ! mp
!      enddo ! np

!     subroutine grades(nlat,nlon,isym,nt,v,w,idvw,jdvw,a,b,mdab,ndab,
!    +                  wvhses,lvhses,work,lwork,ierror)
      call grades(nglat,nglon-1,0,1,gradts,gradte,nglat,nglon,torc,tors &
         ,mdab,nggjlat,wvhses,lvhses,work2,lworkgrades,ierror)
      if (ierror.ne.0) then
        write (6,*) '(3) grades ierror= ',ierror
        stop
      endif

      do i=1,nglon-1
! Add contribution by local radial currents to bsouth,beast
        do j=1,nglat
! Since grades calculates the positive gradient, whereas the magnetic
!   field is the negative gradient, use negative of bsouth,beast on RHS
          bsouth(j,i) = -bsouth(j,i) + gradte(j,i)
          beast(j,i)  = -beast(j,i)  - gradts(j,i)
        enddo ! j
      enddo ! i
! Fill wrap-around points in longitude
      do j=1,nglat
        bsouth(j,nglon) = bsouth(j,1)
        beast(j,nglon) = beast(j,1)
      enddo
!     (4) subroutine math2geov(ig,nlat,nlon,vm,wm,ug,vg,work)
      call math2geov(0,nglat,nglon,bsouth,beast,delbeLEO,delbnLEO,work2)
      if (ierror.ne.0) then
        write (6,*) 'math2geov ierror= ',ierror
        stop
      endif

!     subroutine shses(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
!    +                 wshses,lshses,work,lwork,ierror)
      call shses(nglat,nglon-1,0,1,bup,nglat,nglon,betatotc,betatots,mdab,nggjlat, &
                      wshses,lshses,work2,lworkgrades,ierror)
      if (ierror.ne.0) then
        write (6,*) 'shses ierror= ',ierror
        stop
      endif

! Fill wrap-around points in longitude
      do j=1,nglat
        bup(j,nglon) = bup(j,1)
      enddo
!     (2) subroutine math2geos(ig,nlat,nlon,sm,sg,work)
      call math2geos(0,nglat,nglon,bup,delbuLEO,work2)
      if (ierror.ne.0) then
        write (6,*) 'math2geos ierror= ',ierror
        stop
      endif

      do lat=1,nglat
        do lon=1,nglon
	  delbsclLEO(lon,lat) = delbeLEO(lon,lat)*bhat_gLEO(lon,lat,1) + &
	    delbnLEO(lon,lat)*bhat_gLEO(lon,lat,2) + &
	    delbuLEO(lon,lat)*bhat_gLEO(lon,lat,3)
        enddo ! lon
      enddo ! lat
!
! Calculate equivalent current function psi.
!     subroutine shses(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
!    +                 wshses,lshses,work,lwork,ierror)
! Note: bup is being used as a dummy array below.
      call shses(nglat,nglon-1,0,1,bup,nglat,nglon,psicoefc,psicoefs,mdab,nggjlat, &
                      wshses,lshses,work2,lworkgrades,ierror)
      if (ierror.ne.0) then
        write (6,*) 'shses ierror= ',ierror
        stop
      endif
! Fill wrap-around points in longitude
      do j=1,nglat
        bup(j,nglon) = bup(j,1)
      enddo
!     (2) subroutine math2geos(ig,nlat,nlon,sm,sg,work)
      call math2geos(0,nglat,nglon,bup,psi,work2)
      if (ierror.ne.0) then
        write (6,*) 'math2geos ierror= ',ierror
        stop
      endif
!
! Calculate global fields at all heights
      do k=1,nggjhgt+1

      aob = re/(re+ggjbot(k))
      aobnp1 = aob
      do np=2,nggjlat
        aobnp1 = aobnp1*aob
        n = np - 1
        do mp=1,mdab
! Add contributions of internal potential and radial field to vorc,vors
!   and betac,betas at bottoms of all current layers.
          vorc(mp,np,k) = vintorc(mp,np)*aobnp1 + vorc(mp,np,k) 
          vors(mp,np,k) = vintors(mp,np)*aobnp1 + vors(mp,np,k) 
          betac(mp,np,k) = betaintc(mp,np)*aobnp1*aob + betac(mp,np,k) 
          betas(mp,np,k) = betaints(mp,np)*aobnp1*aob + betas(mp,np,k) 
        enddo ! mp
      enddo ! np

      enddo ! k

!     subroutine grades(nlat,nlon,isym,nt,v,w,idvw,jdvw,a,b,mdab,ndab,
!    +                  wvhses,lvhses,work,lwork,ierror)
      call grades(nglat,nglon-1,0,nggjhgt+1,bsouth3D,beast3D,nglat,nglon,vorc,vors &
         ,mdab,nggjlat,wvhses,lvhses,work2,lworkgrades,ierror)
! Note: at this point, bsouth3D,beast3D have the wrong sign, which is
!   taken into account later.
      if (ierror.ne.0) then
        write (6,*) '(4) grades ierror= ',ierror
        stop
      endif
!
! Get radial magnetic field 
!     subroutine shses(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
!    +                 wshses,lshses,work,lwork,ierror)
      call shses(nglat,nglon-1,0,nggjhgt+1,bup3D,nglat,nglon,betac,betas,mdab,nggjlat, &
                      wshses,lshses,work2,lworkgrades,ierror)
      if (ierror.ne.0) then
        write (6,*) 'shses ierror= ',ierror
        stop
      endif
      do k=1,nggjhgt+1
! Fill wrap-around points in longitude
        do j=1,nglat
          bup3D(j,nglon,k) = bup3D(j,1,k)
        enddo
!     (2) subroutine math2geos(ig,nlat,nlon,sm,sg,work)
        call math2geos(0,nglat,nglon,bup3D(1,1,k),delbu(1,1,k),work2)
        if (ierror.ne.0) then
          write (6,*) 'math2geos ierror= ',ierror,'  k=',k
          stop
        endif
      enddo
!
      gradts3D = 0.
      gradte3D = 0.
! Get gradient of tau at bottoms of layers
!     subroutine grades(nlat,nlon,isym,nt,v,w,idvw,jdvw,a,b,mdab,ndab,
!    +                  wvhses,lvhses,work,lwork,ierror)
! Whereas gradts3D and gradte3D refer to the bottom of the current layer,
!   ar,ai refer to the top of the current layer, so we shift the
!   height index by 1.
      call grades(nglat,nglon-1,0,nggjhgt,gradts3D(1,1,2) &
         ,gradte3D(1,1,2),nglat,nglon,ar,ai &
         ,mdab,nggjlat,wvhses,lvhses,work2,lworkgrades,ierror)
      if (ierror.ne.0) then
        write (6,*) '(5) grades ierror= ',ierror
        stop
      endif
!
      do k=1,nggjhgt+1
        do i=1,nglon-1
          do j=1,nglat

! Since grades calculates the positive gradient, whereas the magnetic
!   field is the negative gradient, reverse the signs of bsouth3D,beast3D
!   and add the toroidal field.
            bsouth(j,i) = -bsouth3D(j,i,k) + gradte3D(j,i,k)
            beast(j,i)  = -beast3D(j,i,k)  - gradts3D(j,i,k)
          enddo ! j
        enddo ! i
! Fill wrap-around points in longitude
        do j=1,nglat
          bsouth(j,nglon) = bsouth(j,1)
          beast(j,nglon) = beast(j,1)
        enddo
!     (4) subroutine math2geov(ig,nlat,nlon,vm,wm,ug,vg,work)
        call math2geov(0,nglat,nglon,bsouth,beast,delbe(1,1,k),delbn(1,1,k),work2)
        if (ierror.ne.0) then
          write (6,*) 'math2geov ierror= ',ierror,'  k=',k
          stop
        endif

        do i=1,nglon
          do j=1,nglat
            delbscl(i,j,k) = delbe(i,j,k)*bhat_g(i,j,k,1) + &
      	      delbn(i,j,k)*bhat_g(i,j,k,2) + &
      	      delbu(i,j,k)*bhat_g(i,j,k,3)
          enddo ! j
        enddo ! i
      enddo ! k
!      
      return
      end subroutine calc_B
!---------------------------------------------------------------------- 
    end module edyn3D_delB_module
!---------------------------------------------------------------------- 
