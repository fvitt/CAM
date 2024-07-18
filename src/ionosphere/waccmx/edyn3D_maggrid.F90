module edyn3D_maggrid

   use shr_kind_mod,   only : r8 => shr_kind_r8            ! 8-byte reals
   use cam_logfile,    only: iulog
!   use edyn_params,    only: finit

!   use edyn3D_params, only: nmlat_h,nmlatS2_h,nmlat_T1,nmlat_T2,nmlon, &
!        nmlonp1,ylatm,ylonm,ylatm_s,ylonm_s,pi,rho,rho_s,rtd,dtr,re=>rearth_m,r0,h0, &
!        m2km,km2m,nhgt_fix,nhgt_fix_r,hgt_fix,hgt_fix_r,ha,ha_s

   implicit none

   !
   ! Global geomagnetic grid:
   !
!   integer, protected ::       &
!        nmlat, &   ! number of mag latitudes
!        nmlath, &  ! index of magnetic equator
!        nmlon, &   ! number of mag longitudes
!        nmlonp1    ! number of longitudes plus periodic point

   !
   ! geomagnetic grid resolution parameters:
   !
!   integer, protected :: res_nlev
!   integer, protected :: res_ngrid

   !
   ! Mag grid coordinates:
   !
!   real(r8), allocatable, protected :: &
!        ylatm(:),   & ! magnetic latitudes (radians)
!        ylonm(:),   & ! magnetic longitudes (radians)
!        gmlat(:),   & ! magnetic latitudes (degrees)
!        gmlon(:)      ! magnetic longitudes (degrees)
!   real(r8), protected :: dlonm,dlatm
   !
   ! Level coordinates will be same as geographic levels:
   !
!   integer, protected :: nmlev ! number of levels (same as nlev in geographic)

!   real(r8), allocatable, protected :: &
!        rcos0s(:),    & ! cos(theta0)/cos(thetas)
!        dt0dts(:),    & ! d(theta0)/d(thetas)
!        dt1dts(:)       ! dt0dts/abs(sinim) (non-zero at equator)


!   real(r8), protected :: table(91,2) = finit

!    logical, private :: debug = .false. ! set true for prints to stdout at each call
   logical, private :: debug = .true. ! set true for prints to stdout at each call


 contains

!-------------------------------------------------------------------------------------------

      subroutine gen_highres_grid
!
! Generate the field line, latitude, and longitude grid for dynamo with P, R, S1, S2 points
!
      use edyn3D_params, only: nmlat_h,nmlatS2_h,nmlat_T1,nmlat_T2,nmlon, &
        nmlonp1,ylatm,ylonm,ylatm_s,ylonm_s,pi,rho,rho_s,rtd,dtr,re=>rearth_m,r0,h0, &
        m2km,km2m,nhgt_fix,nhgt_fix_r,hgt_fix,hgt_fix_r,ha,ha_s

      implicit none
      !
      ! Local variables:
      !
      integer :: j,jns,k,isn                ! Loop indexing variables
      real(r8) :: b,c,d,e,f,g,h,i           ! Variables used in calculating grid regions
      real(r8) :: d1,d2,d3,d4,d5            ! Apex latitudes for defining grid latitude regions (Degrees)
      real(r8) :: lam1,lam2,lam3,lam4,lam5  ! Same as d latitudes above but in radians (Radians)
      real(r8) :: pio2                      ! pi/2
      real(r8) :: y,y1,y2,y3,y4,y5          ! Definition variables for eight latitude regions
      real(r8) :: lam                       ! Latitude variable used in calculating ranges for latitude regions (Radians)
      real(r8) :: th0,lamb,yb,d0,hb         ! Initial coordinate variables used to set latitude region grid points
      real(r8) :: lamc,yc,dc,hs,dhdy        ! Initial coordinate variables used to set latitude region grid points
      real(r8) :: rhob,rhoc,rho1,rho_loc    ! Cosine of magnetic latitudes
      real(r8) :: ymax                      ! y value for pole
      real(r8) :: ha_loc                    ! Local variables for apex height and cosine of latitude
      real(r8) :: hc,h1,r0km,rekm,h0km      ! Variables related to earth radius
      real(r8) :: cos_avg                   ! Average of cosine of latitude for S2-grid latitude calculation
      real(r8) :: dlonm                     ! Delta used to calculate magnetic longitude
      !
      ! Set relative densities of grid points (dy/dlam) at midlatitudes (d), auroral zone (f), and polar cap (h)
      !
      parameter (d=6._r8, f=22._r8, h=20._r8)
      !
      ! Set region boundary values y at the upper boundary of Region I (0<y<yb) and Region II (yb<y<yc)
      !
      parameter (yb = 4._r8, yc = 18._r8)
      !
      ! Set values of magnetic latitudes (degrees) separating Regions III-IV (d1), IV-V (d2), V-VI (d3), VI-VII (d4), and VII-VIII (d5)
      !
      parameter (d1=30._r8, d2=50._r8, d3=60._r8, d4=72._r8, d5=82._r8)
      !
      ! Set values of height hs and height increment dhdy, where ha=h0+hs*y in Region I and ha=h0+hs*y+.5*dhdz*(y-yb)^2 in Region II.
      !
      parameter (hs=8, dhdy=8._r8)    ! [km]
      !
      !      logical, parameter :: debug=.false.    ! [km]
      if(debug) write(iulog,*) 'edyn3D Mag Grid settings:'
      !
      ! Calculate bottom height of dynamo grid from values of the radius of the earth and factors to extend beyond
      !
      r0km = r0*m2km
      rekm = re*m2km
      h0km = h0*m2km
      pio2 = pi*0.5_r8               ! pi/2
      !
      ! Calculate apex height hb of field line for which y=yb which is between Region I and Region II of grid.
      !
      hb   = h0km + hs*yb            ! ~ 110 km
      !
      ! Calculate values of rho and lam corresponding to apex height hb.
      !
      rhob = sqrt(r0km/(hb+rekm))    ! 0.9976829
      lamb = acos(rhob)          ! 0.0680087 [rad]
      !
      ! Calculate values of h, rho and lam corresponding to y=yc.
      !
      hc   = h0km + hs*yc + .5_r8*dhdy*(yc-yb)**2
      rhoc = sqrt(r0km/(hc+rekm))
      lamc = acos(rhoc)
      !
      ! Convert latitude region values from degrees to radians
      !
      lam1 = d1*dtr
      lam2 = d2*dtr
      lam3 = d3*dtr
      lam4 = d4*dtr
      lam5 = d5*dtr
      !
      ! Calculate apex height h1 corresponding to y=y1
      !
      h1 = r0km/cos(lam1)**2 - rekm
      if(debug) write(iulog,*) 'hb =',hb,'    hc =',hc,'    h1 =',h1
      !
      ! Convert lamb and lamc to radians
      !
      d0 = lamb*rtd
      dc = lamc*rtd
      if(debug) write(iulog,*) 'd0 =',d0,'    dc =',dc
      !
      ! Calculate rho1, corresponding to lam1 (and h1 and y1)
      !
      rho1 = cos(lam1)
      if(debug) write(iulog,*) 'rhob =',rhob,'    rhoc =',rhoc,'    rho1 =',rho1
      !
      ! Determine the value of c that makes d(ha)/d(y) continuous at y=yc
      ! and the value of b that makes d(y)/d(lam) continuous at y=y1
      !
      c = hs + dhdy*(yc-yb)
!      b = (2*r0km*sin(lam1)/(d*rho1**3)+c)/(r0km/rho1**2-rekm-hc)
      !
      !ADR230827 Above formula appears to be in error, and should be below with a -c
      !
      b = (2*r0km*sin(lam1)/(d*rho1**3)-c)/(r0km/rho1**2-rekm-hc)
      !
      if(debug) write(iulog,*) 'b =',b,'    c =',c
      !
      ! Determine y values y1-y5 corresponding to latitude boundaries lam1-lam5
      !
      !
!      y1 = yc + alog(b*(h1-hc)/c+1)/b
      y1 = yc + log(b*(h1-hc)/c+1)/b
      if(debug) write(iulog,*) 'y1 =',y1
      y2 = y1 + d*(lam2-lam1)
      if(debug) write(iulog,*) 'y2 =',y2
      e = (f-d)/(2.e0_r8*(lam3-lam2))
      y3 = y2 + d*(lam3-lam2) + e*(lam3-lam2)**2
      if(debug) write(iulog,*) 'y3 =',y3
      y4 = y3 + f*(lam4-lam3)
      if(debug) write(iulog,*) 'y4 =',y4
      g = (f-h)/(2.e0_r8*(lam5-lam4))
      y5 = y4 + f*(lam5-lam4) - g*(lam5-lam4)**2
      if(debug) write(iulog,*) 'y5 =',y5
      ymax = y5 + h*(pio2-lam5)
      if(debug) write(iulog,*) 'ymax =',ymax
      !
      if(debug) write(iulog,*) 'j    rho         ha     lam(deg)  lam(rad) '
      !
      ! Step in latitude along the lower boundary at h=h0 from the equator to the north magnetic pole, using
      !   a constant increment of y, to calculate ha, lam, and rho for the P and S1 grid points.
      !   Different formulas are used for Regions I-VIII.
      !   Note that j in this do-loop steps from equator to pole, opposite to j in the notes that
      !   steps from the pole to the equator.  [Changing j->k would remove this confusion.]
      !
      do j=1,nmlat_h   ! goes from equator to pole S1, P points
	y = (j-1)*ymax/float(nmlat_h-1)
        !
        ! Next line added april 2015 but not used here since above line matches ADR230827 notes
        !
!	y = (float(j)-.5_r8)*ymax/(float(nmlat_h)-.5_r8)
        !
        ! Region I grid points
        !
        ha_loc = h0km + hs*y
        rho_loc = sqrt(r0km/(ha_loc+rekm))
        lam = acos(rho_loc)
        !
        ! Region II grid points
        !
	if (y.gt.yb.and.y.le.yc) then
          ha_loc = h0km + hs*y +.5_r8*dhdy*(y-yb)**2
          rho_loc = sqrt(r0km/(ha_loc+rekm))
          lam = acos(rho_loc)
        endif
        !
        ! Region III grid points
        !
	if (y.gt.yc.and.y.le.y1) then
          ha_loc = hc + c*(exp(b*(y-yc))-1)/b
          rho_loc = sqrt(r0km/(ha_loc+rekm))
          lam = acos(rho_loc)
        endif
        !
        ! Region IV grid points
        !
	if (y.gt.y1.and.y.le.y2) then
       	  lam = (y-y1)/d + lam1
          rho_loc = cos(lam)
          ha_loc = r0km/rho_loc**2 - rekm
        endif
        !
        ! Region V grid points
        !
	if (y.gt.y2.and.y.le.y3) then
      	  lam = lam2 + (sqrt(d**2+4._r8*e*(y-y2))-d)/(2._r8*e)
          rho_loc = cos(lam)
          ha_loc = r0km/rho_loc**2 - rekm
        endif
        !
        ! Region VI grid points
        !
	if (y.gt.y3.and.y.le.y4) then
      	  lam = (y-y3)/f + lam3
          rho_loc = cos(lam)
          ha_loc = r0km/rho_loc**2 - rekm
        endif
        !
        ! Region VII grid points
        !
	if (y.gt.y4.and.y.le.y5) then
      	  lam = lam4 + (f-sqrt(f**2-4._r8*g*(y-y4)))/(2._r8*g)
          rho_loc = cos(lam)
          ha_loc = r0km/rho_loc**2 - rekm
        endif
        !
        ! Region VIII grid points
        !
	if (y.gt.y5) then
      	  lam = (y-y5)/h + lam5
          rho_loc = cos(lam)
          ha_loc = 9999999._r8
          if (j.ne.nmlat_h) ha_loc = r0km/rho_loc**2 - rekm
        endif
        !
        ! Copy grid points for northern hemisphere/southern hemisphere symmetry
        ! jns corresponds to j of notes, increasing from 1 at pole to J=nmlat_h at equator.
        !
	jns = nmlat_h-j+1     ! jns: nmlat_h (eq) to 1 (pole)
        ylatm(jns,2) =  lam   ! northern hemisphere
        ylatm(jns,1) = -ylatm(jns,2) ! southern hemisphere
        rho(jns,1)   = cos(ylatm(jns,1))  ! cos(ylatm)
        rho(jns,2)   = cos(ylatm(jns,2))  ! cos(ylatm)
        !
        ! Overwrite numerical inaccuracy in calculating rho for j=nmlat_h / jns = 1:
        !
	if (jns.eq.1) then
          rho(jns,1)   = 0._r8
          rho(jns,2)   = 0._r8
	endif
	ha(jns) = ha_loc*km2m
        !
        ! Set layer heights hgt_fix to correspond to apexes of P field lines.
        !
        if(j.le.nhgt_fix) hgt_fix(j) = ha_loc*km2m
        !if(debug .and. j.le.nhgt_fix) write(iulog,'(i3,1(x,f15.2))') j,hgt_fix(j)
	if(debug) write(iulog,20) jns,rho(jns,2),ha_loc*km2m,90._r8*ylatm(jns,2)/pio2,ylatm(jns,2)
   20   format(i4,f10.6,2f15.2,f10.4)
!   20   format(i2,f12.8,f11.2,f11.3,f11.5,f11.3,f11.5)
      enddo ! goes from equator to pole S1, P points
      !
      ! Above is for P-grid points for potential.  Need the same for S2-grid point latitude, cosine of latitude, and apex height values.
      ! S2 points are for conductivities, wind, and equatorward/downward current
      !
      if(debug) write(iulog,*)   " "
      !
      ! In the following do-loop j corresponds to jns of the previous loop
      !
      do j=1,nmlatS2_h
         cos_avg = cos(ylatm(j,2))+cos(ylatm(j+1,2))
         cos_avg = cos_avg*0.5_r8
         ylatm_s(j,2) = acos(cos_avg)
         ylatm_s(j,1) = -ylatm_s(j,2)
         rho_s(j,1) = cos_avg  ! cos(ylatm_s)
         rho_s(j,2) = cos_avg  ! cos(ylatm_s)
         ha_s(j)    = (r0km/cos_avg**2 - rekm)*km2m
	 if(debug) write(iulog,'(i4,3(1x,f10.6))') j,rho_s(j,2),90._r8*ylatm_s(j,2)/pio2,ylatm_s(j,2)
!	 if(debug) write(iulog,*) j,ylatm(j,2),ylatm(j+1,2),ylatm_s(j,2)
      enddo
      !
      ! Also need R-grid point values.  Set up height levels [m] for R-grid on which Je3 is defined. These are apex heights of ylatm_s S2 points
      !
      isn=2  ! need to specify only one hemisphere
	if(debug) write(iulog,*)   " "
      do k=2,nhgt_fix_r
        j = nmlatS2_h-k+2
        hgt_fix_r(k) = apex_height(rho_s(j,isn))
	if(debug)  write(iulog,'(i3,1(1x,f15.2))') k,hgt_fix_r(k)
      enddo
      k=1
      hgt_fix_r(1) = h0
      if(debug)  write(iulog,'(i4,3(1x,f15.2))') k,hgt_fix_r(k)
      !
      ! Calculate magnetic longitudes [radians]
      !
      dlonm = 2._r8*pi/float(nmlon)
      do j=1,nmlon
        ylonm(j) = -pi+float(j-1)*dlonm
      enddo ! i=1,nmlon
      !
      ! Calculate magnetic longitudes [radians] for S1-grid points which are between the P-grid points
      !
      do j=1,nmlon
        ylonm_s(j) =ylonm(j)+0.5_r8*dlonm
      enddo ! i=1,nmlon
!
     contains

      function apex_height(rho)

      use edyn3D_params, only:  re=>rearth_m,r0
      ! Calculate apex height of each fieldline
      ! eq.(3.3) [Richmond 1995]
      ! ha = R/cos^2lam_m - R_E  with lam_m modified apex latitude
      !                               R_E mean Earth radius (code re)
      !                               h_R reference height  (code h0)
      !                               R = R_E + h_R         (code r0)

      implicit none
      real(r8), intent(in):: rho
      real(r8) :: apex_height
      !
      apex_height = r0/(rho)**2 -re ! [m]
      !
      end function apex_height

    end subroutine gen_highres_grid
!
!---------------------------------------------------------------------- 

      subroutine edyn3D_gen_ggj_grid
     
  	! Sets up geographic grid (in radians and m) for calculating currents
  	!  to use for magnetic perturbations.
  	! In this version, colatitude points are Gaussian, which exclude the poles.
  	! i=1 is 0 geographic longitude. 
  	! j=1 is near, but not at, the North Pole (pi/2).
  	! j=nggjlat is near, but not at, the South Pole (-pi/2).
   
  	use edyn3D_params, only: pi,nhgt_fix,nhgt_fix_r,hgt_fix,hgt_fix_r, &
  	  nggjlon,nggjlat,nggjhgt,ggjlon,ggjclat,wts,ggjhgt,ggjtop, &
  	  ktop,k_fix_ggjbot,delBsolution

  	implicit none
  !	 
  	integer :: i,lwork,ierror,k,kk,ktop_tmp,kend
  	real(r8) :: dggjlon,k_tmp(nhgt_fix_r),w
      
  	  k_tmp = hgt_fix_r - 400000._r8
  	  ktop_tmp = minloc(abs(k_tmp),dim=1)
  	  write(iulog,'(a4,i4)') 'ktop',minloc(abs(k_tmp))
  	  
  	dggjlon = 2._r8*pi/float(nggjlon)
  	do i=1,nggjlon
  	  ggjlon(i) = float(i-1)*dggjlon
  	enddo ! i
  ! Gaussian colatitude points (radians) and weights
  	call gaqd(nggjlat,ggjclat,wts,w,lwork,ierror)
  	if (ierror.ne.0) then
  	  write (iulog,*) 'gaqd error code ',ierror
  	  stop
  	endif
  	k_fix_ggjbot(1) = 1
  	if (delBsolution.eq.'full_hgt_delB  ') then
  ! In this version, heights for current layers used for magnetic
  !   perturbation calculations are set to heights of rho and QD grids,
  	  ktop = nhgt_fix_r
  	  do k=1,nggjhgt
  	    ggjhgt(k) = hgt_fix(k)
  ! Note that ggjtop refers to the height of the top of the ggj layer,
  !   while hgt_fix_r refers to the height of the bottom of the hgt_fix
  !   layer.
  	    ggjtop(k) = hgt_fix_r(k+1)
  	    k_fix_ggjbot(k+1) = k+1
  	 enddo
  	elseif (delBsolution.eq.'quick_ground	 ') then
  ! In this version, all horizontal currents up to 397140 m are combined
  !   to a single layer at 110 km, and delB is calculated only at the
  !   ground.
  	  k_tmp = hgt_fix_r - 400000._r8
  	  ktop = minloc(abs(k_tmp),dim=1)
  	  !ktop = 39  ! hgt_fix_r(39) = 397140. m
  	  ggjhgt(1) = 110.e3_r8  ! m
  	  ggjtop(nggjhgt) = hgt_fix_r(ktop)  
  	  k_fix_ggjbot(nggjhgt+1) = ktop
  	elseif (delBsolution.eq.'quick_ground,LEO') then
  ! In this version, all horizontal currents up to 397140 m are combined
  !   to a single layer at 110 km, and delB is calculated only at the
  !   ground and at 397139 m.
  	  k_tmp = hgt_fix_r - 400000._r8
  	  ktop = minloc(abs(k_tmp),dim=1)
  	  !ktop = 39  ! hgt_fix_r(39) = 397140. m
  	  ggjhgt(1) = 110.e3_r8  ! m
  	  ggjtop(nggjhgt) = hgt_fix_r(ktop)  
  	  k_fix_ggjbot(nggjhgt+1) = ktop
  	elseif (delBsolution.eq.'ground,LEO	') then
  ! In this version, horizontal currents are combined in thick layers,
  !   and delB is calculated at the ground and at heights between 299702 m
  !   and 540303 m, and between 764409 m and 948105 m.
  	  k_tmp = hgt_fix_r - 998000._r8
  	  ktop = minloc(abs(k_tmp),dim=1)			 ! ktop = 54  ! hgt_fix_r(54) = 948105 m 
  	  
  	  kk = minloc(abs(hgt_fix-109000._r8),dim=1)
  	  ggjhgt(1) = hgt_fix(kk)				 ! hgt_fix(14)      ! 109486 m
  	  k_fix_ggjbot(2) = minloc(abs(hgt_fix_r-140000._r8),dim=1) ! k_fix_ggjbot(2) = 21
  	  ggjtop(1) = hgt_fix_r(k_fix_ggjbot(2))		 ! hgt_fix_r(21)    ! 139381 m
  	   
  	  kk = minloc(abs(hgt_fix-140000._r8),dim=1)
  	  ggjhgt(2) = hgt_fix(kk)				 ! hgt_fix(25)      ! 179576 m
  	  k_fix_ggjbot(3) = minloc(abs(hgt_fix_r-22000._r8),dim=1)  ! k_fix_ggjbot(3) = 29
  	  ggjtop(2) = hgt_fix_r(k_fix_ggjbot(3))		 ! hgt_fix_r(29)    ! 222139 m
  	  
  	  kk = minloc(abs(hgt_fix-225000._r8),dim=1)
  	  ggjhgt(3) = hgt_fix(kk)				 ! hgt_fix(31)      ! 258343 m
  	  k_fix_ggjbot(4) = minloc(abs(hgt_fix_r-300000._r8),dim=1) ! k_fix_ggjbot(4) = 34
  	  ggjtop(3) = hgt_fix_r(k_fix_ggjbot(4))		 ! hgt_fix_r(34)  ! 299702 m
  	  
  	  k = k_fix_ggjbot(4)					 ! k= 34
  	  kend = minloc(abs(hgt_fix_r - 540000._r8),dim=1) 	 ! ggjtop = 540303 m  
  	  kend = kend-k+1					 ! how many levels between	     
  	  do kk=4,kend  					 ! ggjtop = 317599 m to 540303 m
  	    ggjhgt(kk) = hgt_fix(k)
  	    ggjtop(kk) = hgt_fix_r(k+1)
  	    k_fix_ggjbot(kk+1) = k+1
  	    k = k+1
  	  enddo
  	  kk = minloc(abs(hgt_fix-580000._r8),dim=1)
  	  ggjhgt(kend+1) = hgt_fix(kk)  				   ! ggjhgt(15) = hgt_fix(46)	   ! 580485 m
  	  k_fix_ggjbot(kend+2) = minloc(abs(hgt_fix_r-630000._r8),dim=1)	   ! k_fix_ggjbot(16) = 48
  	  ggjtop(kend+1) = hgt_fix_r(k_fix_ggjbot(kend+2))		   ! ggjtop(15) = hgt_fix_r(48)    ! 622616 m
  	  kk = minloc(abs(hgt_fix-670000._r8),dim=1)
  	  ggjhgt(kend+2) = hgt_fix(kk)  				   ! ggjhgt(16) = hgt_fix(50)	   ! 696619 m
  	  k_fix_ggjbot(kend+3) =  minloc(abs(hgt_fix_r-760000._r8),dim=1)     ! k_fix_ggjbot(17) = 52
  	  ggjtop(kend+2) = hgt_fix_r(k_fix_ggjbot(kend+3))		   ! ggjtop(16) = hgt_fix_r(52)    ! 764409 m
  	  kk = minloc(abs(hgt_fix-790000._r8),dim=1)
  	  ggjhgt(kend+3) = hgt_fix(kk)  				   ! ggjhgt(17) = hgt_fix(52)	   ! 793704 m
  	  k_fix_ggjbot(kend+4) = minloc(abs(hgt_fix_r-840000._r8),dim=1)	   ! k_fix_ggjbot(18) = 53
  	  ggjtop(kend+3) = hgt_fix_r(k_fix_ggjbot(kend+4))		   ! ggjtop(17) = hgt_fix_r(53)    ! 837615 m
  	  kk = minloc(abs(hgt_fix-880000._r8),dim=1)
  	  ggjhgt(kend+4) = hgt_fix(kk)  				   ! ggjhgt(18) = hgt_fix(53)	   ! 881931 m
  	  k_fix_ggjbot(nggjhgt+1) = ktop				   ! ???? should this be kend+5 instead of nggjhgt+1
  	  ggjtop(kend+4) = hgt_fix_r(k_fix_ggjbot(nggjhgt+1))		   ! ggjtop(18) = hgt_fix_r(54)    ! 948105 m
  	else
  	  write (iulog,*) 'delBsolution must be specified in params'
  	  stop
  	endif
  !	  
       end subroutine edyn3D_gen_ggj_grid
!
!-----------------------------------------------------------------------

      subroutine edyn3D_gen_qd_grid
!     
!  set up quasi dipole grid (longitude, latitude and height)   
!     
      use edyn3D_params, only: nmlon,ylonm,ylonm_s,nlat_qd,nlat_qd_h,pi,rtd, &
                               nhgt_fix,nhgt_fix_r,hgt_fix,hgt_fix_r
!     
      implicit none
!
      integer :: i,l,k
      real(r8) :: dlatm,fac_r      

      real(r8) :: lat_qd_ed(nlat_qd)    ! quasi latitude of edge of volume l-.5
      real(r8) :: lat_qd_mp(nlat_qd-1)  ! quasi latitude of midpoint of volume l
      real(r8) :: lon_qd_ed(nmlon)      ! quasi longitude of edge of volume     
      real(r8) :: lon_qd_mp(nmlon)      ! quasi longitude of midpoint of volume     
      real(r8) :: hgt_qd_mp(nhgt_fix)   ! height of quasi dipole grid = p height level
      real(r8) :: hgt_qd_ed(nhgt_fix_r) ! height of quasi dipole grid = r height level
!    
      lon_qd_mp = ylonm    ! [rad] same as p-points
      lon_qd_ed = ylonm_s  ! [rad] same as s1-points
!      
      dlatm = pi/float(nlat_qd-1)
      do l =1,nlat_qd
        lat_qd_ed(l) = pi*float(l-nlat_qd_h)/float(nlat_qd-1) ! equally distributed
        !write(6,*) 'lat_qd_ed',i,lat_qd_ed(i)*rtd
      end do
!      
! midpoints are in the middle of the volume l
      do l =1,nlat_qd-1
        lat_qd_mp(l) = 0.5_r8*(lat_qd_ed(l)+lat_qd_ed(l+1)) ! equally distributed
        !write(6,*) 'lat_qd_mp',i,lat_qd_mp(i)*rtd
      end do
!       
      hgt_qd_mp = hgt_fix    ! height of p,s1,s2 points
      hgt_qd_ed = hgt_fix_r  ! height of r points
        
!         
      end subroutine edyn3D_gen_qd_grid

!---------------------------------------------------------------------- 

      subroutine edyn3D_gen_geo_grid
!     
      use edyn3D_params, only: nglon,nglat,glon,glat
!
      implicit none
!     
      integer :: ilateq,i
      real(r8) :: dlon,dlat
!      
! Set up geographic grid for outputting magnetic perturbations.  Also used for 
! interpolation between physics grid and 3D dynamo grid
!
      dlon = 360._r8/float(nglon-1)
      dlat = 180._r8/float(nglat-1)
      ilateq = (nglat+1)/2
      do i=1,nglon
        glon(i) = 0._r8+(i-1)*dlon
      enddo
      do i=1,nglat
        glat(i) = (i-ilateq)*dlat
      enddo     
!      
      end subroutine edyn3D_gen_geo_grid

!----------------------------------------------------------------------------- 
      subroutine edyn3D_qcoef
      
!  Calculates coefficients needed to calculate functions Q of Richmond
!  (1974), including functions with odd n-m, used for symmetric
!  high-latitude FAC.  These functions are applied to the top of the
!  region of resolved 3D currents, above which current is assumed to
!  flow only along dipolar field lines out to the apex, and then
!  radially to infinity (for symmetric FAC).
!     
      implicit none
!     
! parameters and arrays for calculating equivalent current function
!   associated with field-aligned currents above ggjtop(nggjhgt)
      integer, parameter :: nmax = 6
      integer :: m,n,p,po2,num,k
      real(r8) :: PMOPMMO(nmax+1), R(0:nmax,0:nmax), SQ2,pi,sq4pi
      real(r8) :: a(0:nmax,0:nmax,0:nmax), mc(0:nmax,0:nmax,0:nmax), &
        md(0:nmax,0:nmax)

      real(r8) :: x

      pi = 4._r8*atan(1._r8)
      SQ2 = sqrt(2.e0_r8)
      sq4pi = sqrt(4._r8*pi)
      R = 0._r8
      
      do m=0,nmax
        if (m.ne.0) PMOPMMO(m) = sqrt(1._r8 + .5_r8/float(m))
        do n=max0(1,m),nmax
          R(n,m) = sqrt(float(n*n-m*m))/sqrt(float(4*n*n-1))
        enddo
      enddo
      
      a = 0._r8
      a(0,0,0) = 1._r8

      do m=1,nmax
        a(m,m,0) = sqrt(float(2*m+1)/float(2*m))*a(m-1,m-1,0)
        n = m + 1
        if (n.gt.nmax) cycle
        a(n,m,0) = sqrt(float(2*m+3))*a(m-1,m-1,0)
        n = m + 2
        if (n.gt.nmax) cycle
        a(n,m,0) = (1._r8 - R(n-1,m)**2 - R(n-2,m)**2)*a(n-2,m,0)/ &
          (R(n,m)*R(n-1,m))
        a(n,m,1) = ((1._r8 - R(n-1,m)**2 - R(n-2,m)**2)*a(n-2,m,1) &
           - a(n-2,m,0))/(R(n,m)*R(n-1,m))
        n = m + 3
        if (n.gt.nmax) cycle
        a(n,m,0) = (1._r8 - R(n-1,m)**2 - R(n-2,m)**2)*a(n-2,m,0)/ &
          (R(n,m)*R(n-1,m))
        a(n,m,1) = ((1._r8 - R(n-1,m)**2 - R(n-2,m)**2)*a(n-2,m,1) &
           - a(n-2,m,0))/(R(n,m)*R(n-1,m))
        if (m+4.gt.nmax) cycle
        do n=m+4,nmax
          do p=0,n-m,2
            po2 = p/2
            x = 0._r8
            if (p.ge.2) x = a(n-2,m,po2-1)
            a(n,m,po2) = ((1._r8-R(n-1,m)**2-R(n-2,m)**2)*a(n-2,m,po2) - &
              x - R(n-2,m)*R(n-3,m)*a(n-4,m,po2))/(R(n,m)*R(n-1,m))
          enddo
        enddo
      enddo
! Now multiply a_nmp of notes by m and divide by (2n-m-p).
      do m=1,nmax
        do n=m,nmax
          do p=0,n-m,2
            po2 = p/2
            a(n,m,po2) = m*a(n,m,po2)/float(2*n-m-p)
! a contains a_nmp of notes multiplied by m and divided by (2n-m-p).
          enddo
        enddo
      enddo
      mc = 0._r8
! Case 1: m,n are both odd
      do m=1,nmax,2
        do n=m,nmax,2
          mc(n,m,0) = a(n,m,0)
          if (2.gt.2*n-m-1) cycle
          do p=2,2*n-m-1,2
            po2 = p/2
            mc(n,m,po2) = a(n,m,po2) + &
              mc(n,m,po2-1)*float(2*n-m-p+1)/float(2*n-m-p)
! mc contains c_nm of notes multiplied by m.
          enddo
        enddo
      enddo
! Case 2: m,n are both even
      md = 0._r8
      do m=2,nmax,2
        do n=m,nmax,2
          mc(n,m,0) = a(n,m,0)
          do p=0,n-m,2
            po2 = p/2
            num = 2*n-m-p
            x = a(n,m,po2)*(num-1)
            do k=1,n-m/2
              if (num.eq.2) exit
              num = num-2
              x = x*(num-1)/float(num)
            enddo
            md(n,m) = md(n,m) + x
! md contains d_nm of notes multiplied by m.
          enddo
          if (2.gt.2*n-m-2) cycle
          do p=2,2*n-m-2,2
            po2 = p/2
            mc(n,m,po2) = a(n,m,po2) + &
              mc(n,m,po2-1)*float(2*n-m-p+1)/float(2*n-m-p)
! mc contains c_nm of notes multiplied by m.
          enddo
        enddo
      enddo

! Test
!      do m=0,nmax
!        do n=m,nmax
!          write (6,'(a1,2i3,7e10.3)') 'a',n,m,(a(n,m,po2),po2=0,nmax)
!        enddo
!      enddo
!      do m=0,nmax
!        do n=m,nmax
!          write (6,'(a2,2i3,7e10.3)') 'mc',n,m,(mc(n,m,po2),po2=0,nmax)
!        enddo
!      enddo
!      do m=0,nmax
!        do n=m,nmax
!          write (6,'(a2,2i3,7e10.3)') 'md',n,m,md(n,m)
!        enddo
!      enddo

! Case 3: If n-m is odd use array a, which contains 
!   a_nmp of notes multiplied by m and divided by (2n-m-p).
!
      end subroutine edyn3D_qcoef

!-----------------------------------------------------------------------------

  subroutine edyn3D_calculate_mf

! Calculation of normalized integrated areas from the pole to an s2
!  surface in the meridional (a1) and horizontal (a3) planes, and of
!  the factors m1f,m2f,m3f, which are independent of magnetic longitude.
!  These latter factors give M1,M2,M3 when divided by F.

       use edyn3D_params, only: nmlatS2_h,hgt_fix,nmlat_h,nhgt_fix,nhgt_fix_r, &
                                hgt_fix_r,r0,ylonm_s,rho_s,re=>rearth_m,pi, &
                                m1f,m2f,m3f

       implicit none
!       
! local variables         
       integer ::  j,k,isn,jmax      
       real(r8) ::  fac1,fac2,fac3,dlonm,rm,rp,ra,rbar
! a1,a3 vary from 0 at the magnetic pole to 1 at the magnetic equator.
       real(r8) ::  a1(nmlat_h+1,nhgt_fix)  ! Normalized integral from pole of M1*F
       real(r8) ::  a3(nmlat_h+1,nhgt_fix_r)! Normalized integral from pole of M3*F
!                   
       dlonm = ylonm_s(2)-ylonm_s(1) ! assumes equidistant longitudinal gridpoints
! Assumes hemispherical symmetry
        isn = 1 !Should get the same result for isn = 2, since rho_s(j,1) = rho_s(j,2)
!  
! In order to include the pole, the first index of a1,a3 represents the
!  location j-0.5, not j+0.5. However, the first index of m1f,m2f,m3f
!  represents the location j+0.5, i.e., the s2 points.
! 
        ! Set a1,a3 to 1 for equator and beyond (points before equator will be
        !  overwritten later).
        a1 = 1.
        a3 = 1.
        ! Set m1f,m2f,m3f to 0 beyond equator (points before equator will be
        !  overwritten later).
        m1f = 0.
        m2f = 0.
        m3f = 0.
        do k=1,nhgt_fix 
           a1(1,k) = 0.
           a3(1,k) = 0.
           jmax = nmlatS2_h - k + 1 ! number of s2 points at level k 
           if (jmax.lt.1) then
! Error trap; this condition should never occur.
              write(6,*) 'Stopped in calc_mf because jmax=',jmax
              stop
           endif
! Normalized radii of the top and bottom of layer k.
           rp = (hgt_fix_r(k+1)+re)/r0           ! r_k+0.5 /R
           rm = (hgt_fix_r(k)+re)/r0             ! r_k-0.5 /R
!
           fac2 = 0.5*pi*(rp - rm)*r0   ! Pi/2*((r_k+0.5)-(r_k-0.5)) 
           fac3 = 2.*dlonm*(hgt_fix(k)+re)**3/r0 ! 2*dlon*r_k^3/R
           do j=1,jmax     ! loop over all s2 points
              fac1 = sqrt(rm)*rho_s(j,isn)      ! sqrt[r_k-0.5/R]*rho(j+0.5) 
! First index of a1,a3 is j+1 because this corresponds to position j+0.5.
              a3(j+1,k) = 1. - sqrt(1.-fac1**2)
              m3f(j,k) = (hgt_fix_r(k)+re)**2 *dlonm*(a3(j+1,k)-a3(j,k))
! Calculate the normalized radius within the layer, rbar, that gives the most
!  accurate value of a1 when a1 is computed by the approximation below.
!  This calculation of rbar assumes the radius of field lines near the
!  equator is parabolic with respect to magnetic latitude.
              ra = 1./rho_s(j,isn)**2
! Prevent ra from getting large enough to affect the numerical accuracy
!  of rbar, which rapidly asymptotes to .5*(rp+rm) as ra increases.
              ra = min(ra,rp+16.*(rp-rm))
              ra = max(ra,rp)  ! ra was less than rp when it should have been equal
              !
              rbar = ra - (2.*(sqrt(ra-rm)**3-sqrt(ra-rp)**3)/(3.*(rp-rm)))**2  ! eq. (45') page 4c r*/R
              a1(j+1,k) = 2.*asin(sqrt(rbar)*rho_s(j,isn))/pi                   ! eq. (47')
              m1f(j,k) = ((hgt_fix(k)+re)/r0)**2.5 *fac2*(a1(j+1,k)-a1(j,k))*r0
! Make sure fac1 does not numerically exceed 1,
!  so that sqrt(1-fac1**2) can be computed.
              fac1 = sqrt(rp)*rho_s(j,isn)        ! sqrt[r_k+0.5/R]*rho(j+0.5) 
              fac1 = min(fac1,1.)
              m2f(j,k) = fac3*(sqrt(1.-rm*rho_s(j,isn)**2) &
                - sqrt(1.-fac1**2))*sqrt(1.-.75*rho_s(j,isn)**2)/rho_s(j,isn)
           enddo !j
            m1f(jmax+1,k) = ((hgt_fix(k)+re)/r0)**2.5 *fac2*(1-a1(jmax+1,k))*r0
           m3f(jmax+1,k) = (hgt_fix_r(k)+re)**2 *dlonm*(1.-a3(jmax+1,k))
        enddo !k
!
! Now do a3,m3f for top level
        k = nhgt_fix_r
        a3(1,k) = 0.
        jmax = nmlatS2_h - k + 1 ! number of s2 points at k-0.5 
        if (jmax.ge.1) then
           rm = (hgt_fix_r(k)+re)/r0
           do j=1,jmax     ! loop over all s2 points (level k has nmlatS2_h-k+1 s2 points)
              fac1 = sqrt(rm)*rho_s(j,isn)      ! sqrt[r_k-0.5/R]*rho(j+0.5) 
              a3(j+1,k) = 1. - sqrt(1.-fac1**2)
              m3f(j,k) = (hgt_fix_r(k)+re)**2 *dlonm*(a3(j+1,k)-a3(j,k))
            enddo !j
            m3f(jmax+1,k) = (hgt_fix_r(k)+re)**2 *dlonm*(1.-a3(jmax+1,k))
        endif
!       
  end subroutine edyn3D_calculate_mf

!---------------------------------------------------------------------- 

end module edyn3D_maggrid
