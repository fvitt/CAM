    module edyn3D_qd_module
!     
! sets up quasi dipole grid: same longitudes as the p-points,
! but different, regular latitude grid
! 
      use edyn3D_params, only: nlat_qd,nlat_qd_h,nmlon,ylonm,ylonm_s,pi,nhgt_fix,&
                               hgt_fix,rtd,nhgt_fix_r,hgt_fix_r,nggjlon,nggjlat,nggjhgt
      use shr_kind_mod,  only: r8 => shr_kind_r8            ! 8-byte reals
!     
      implicit none
!    
      type qd_grid 
        real(r8) :: F1(nhgt_fix)      ! for I1qd points
        real(r8) :: M1q(nhgt_fix)
        real(r8) :: M3q(nhgt_fix_r)
        real(r8) :: F3(nhgt_fix_r)
        real(r8) :: glat(nhgt_fix_r),glon(nhgt_fix_r)
        real(r8) :: Fqd(nhgt_fix)         ! for I1qd points
        real(r8) :: f11(nhgt_fix),f12(nhgt_fix),f21(nhgt_fix),f22(nhgt_fix)
        real(r8) :: f1d2(nhgt_fix),f2d2(nhgt_fix),kd2(nhgt_fix)
	real(r8) :: cosIm(nhgt_fix),DsinI(nhgt_fix)
	real(r8) :: bhat_mp(3,nhgt_fix)
        real(r8) :: bhat(3,nhgt_fix_r)    ! unit vector in magnetic field direction for qd midpoints in lon/lat
      end type qd_grid
!      
      type (qd_grid) :: qd(nmlon,nlat_qd-1) ! at midpoints points in latitude -> nlat_qd-1
!       
      real(r8) :: lat_qd_ed(nlat_qd)    ! quasi latitude of edge of volume l-.5
      real(r8) :: lat_qd_mp(nlat_qd-1)  ! quasi latitude of midpoint of volume l
      real(r8) :: lon_qd_ed(nmlon)      ! quasi longitude of edge of volume     
      real(r8) :: lon_qd_mp(nmlon)      ! quasi longitude of midpoint of volume     
      real(r8) :: hgt_qd_mp(nhgt_fix)   ! height of quasi dipole grid = p height level
      real(r8) :: hgt_qd_ed(nhgt_fix_r) ! height of quasi dipole grid = r height level
!      
! I1(qd) is at i+/-0.5 ,  l+/-0.5, k
! AM 2022-08-08 mapping of LI1
!      real :: a1q(nlat_qd_h,nhgt_fix)             ! normalized area invariant in longitude; edges l+0.5
      real(r8) :: a1q(nlat_qd_h,nhgt_fix,4)            ! normalized area invariant in longitude; edges l+0.5
! AM 2022-08-08 mapping of LI1 
!      real :: wgt1(nlat_qd_h,nhgt_fix)             ! weighting function invariant in longitude; edges l+0.5
      real(r8) :: wgt1(nlat_qd_h,nhgt_fix,4)             ! weighting function invariant in longitude; edges l+0.5
      real(r8) :: q1(nlat_qd_h-1,nhgt_fix)             ! areas factors for eastward current, midpoints
      real(r8) :: LI1qd(nmlon,nlat_qd_h,nhgt_fix,2)    ! integrated I1qd ; edges l+0.5
      real(r8) :: I1qd(nmlon,nlat_qd_h-1,nhgt_fix,2)   ! I1qd ; i+/-0.5; midpoints l
      real(r8) :: Jf1qd(nmlon,nlat_qd_h-1,nhgt_fix,2)  ! Jf1qd ; i+/-0.5; midpoints l
      integer :: jl_qd(nlat_qd_h,nhgt_fix)         ! index in which each element lies l+0.5
!
      real(r8) :: q2(nlat_qd_h,nhgt_fix)               ! areas factors for northward current, midpoints
!
! I3(qd) is at k+/-0.5 i,  l+/-0.5     
      real(r8) :: a3q(nlat_qd_h,nhgt_fix_r)               ! normalized area invariant in longitude edges l+0.5
! 2022-7-26 ADR Weights for cubic spline interpolation
      real(r8) :: wgt3(nlat_qd_h,nhgt_fix_r,4)              ! weighting function invariant in longitude edges l+0.5
!      real :: wgt3(nlat_qd_h,nhgt_fix_r)              ! weighting function invariant in longitude edges l+0.5
      real(r8) :: q3(nlat_qd_h-1,nhgt_fix_r)              ! areas factors for upward current, midpoints
      real(r8) :: LI3qd(nmlon,nlat_qd_h,nhgt_fix_r,2)     ! integrated I3qd; edges l+0.5
      real(r8) :: I3qd(nmlon,nlat_qd_h-1,nhgt_fix_r,2)    ! I3qd ; i k+/-0.5; midpoints l
      real(r8) :: Jrqd(nmlon,nlat_qd_h-1,nhgt_fix_r,2)    ! Jrqd ; i k+/-0.5; midpoints l
      integer :: jl3_qd(nlat_qd_h,nhgt_fix_r)         ! index in which each element lies l+0.5
!
! I2(qd) is at l+0.5 , i, k
      real(r8) :: I2qd(nmlon,nlat_qd_h,nhgt_fix,2)      ! I2qd ;  at edges l+0.5 see eq. 247 (page 7; 2015 April 20)
      real(r8) :: Jf2qd(nmlon,nlat_qd_h,nhgt_fix,2)      ! Jf2qd; at edges l+0.5 (see eq. 255)
!
! Jhor &Jr at point i,j,k 
      real(r8), dimension(nmlon,nlat_qd_h-1,nhgt_fix,2) :: &
         Jf1,&	        ! Jf1
         Jf2,&	        ! Jf2
         Jr,&	        ! Jr
         Jf1hor,&	! Jf1hor
         Jf2hor,&	! Jf2hor
         Jeej,&	        ! total Jeastward in f1 dir with Jf2 mapped into it too
         Jfac_qd,&	! field aligned current
         Je2J_qd,&	! Je2 component of total current
         Je2J_qA,&	! Je2 component of total current; alternative calculation
         Je2J_qA1,&	! Je2 component of total current; alternative calculation
         Je2J_qA2,&	! Je2 component of total current; alternative calculation
         g13,&  	! k^*g1
         g23		! k^*g2
      real(r8), dimension(nggjlat,nggjlon,nggjhgt) :: &
         Keast,Ksouth,Jrtopgg
      real(r8), dimension(nggjlat,nggjlon) :: etagg    ! equivalent current
!          function at height of ggtop(nggjhgt), representing the
!          effects of currents above this height on magnetic
!          perturbations below and at this height.

!     
      contains
!-----------------------------------------------------------------------
      subroutine gen_qd_grid
!     
!  set up quasi dipole grid (longitude, latitude and height)   
!     
      implicit none
!
      integer :: i,l,k
      real(r8) :: dlatm,fac_r      
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
        lat_qd_mp(l) = 0.5*(lat_qd_ed(l)+lat_qd_ed(l+1)) ! equally distributed
	!write(6,*) 'lat_qd_mp',i,lat_qd_mp(i)*rtd
      end do
!       
      hgt_qd_mp = hgt_fix    ! height of p,s1,s2 points
      hgt_qd_ed = hgt_fix_r  ! height of r points
	
!         
      end subroutine gen_qd_grid
!-----------------------------------------------------------------------
      end module edyn3D_qd_module
!-----------------------------------------------------------------------
