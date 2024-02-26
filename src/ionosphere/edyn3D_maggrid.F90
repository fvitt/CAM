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

!!-------------------------------------------------------------------------------------------
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
      e = (f-d)/(2.e0*(lam3-lam2))
      y3 = y2 + d*(lam3-lam2) + e*(lam3-lam2)**2
      if(debug) write(iulog,*) 'y3 =',y3
      y4 = y3 + f*(lam4-lam3)
      if(debug) write(iulog,*) 'y4 =',y4
      g = (f-h)/(2.e0*(lam5-lam4))
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
      	  lam = lam2 + (sqrt(d**2+4.*e*(y-y2))-d)/(2._r8*e)
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
	 if(debug) write(iulog,'(i4,3(x,f10.6))') j,rho_s(j,2),90.*ylatm_s(j,2)/pio2,ylatm_s(j,2)
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
	if(debug)  write(iulog,'(i3,1(x,f15.2))') k,hgt_fix_r(k)
      enddo
      k=1
      hgt_fix_r(1) = h0
      if(debug)  write(iulog,'(i4,3(x,f15.2))') k,hgt_fix_r(k)
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

end module edyn3D_maggrid
