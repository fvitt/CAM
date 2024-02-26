module edyn3D_params
   !
   ! Constants for edynamo.
   !
   use shr_kind_mod,  only: r8 => shr_kind_r8            ! 8-byte reals
   use shr_const_mod, only: rearth_m => SHR_CONST_REARTH ! meters
   use physconst,     only: pi

   implicit none

   private

   public :: nmlat_h,nmlatS2_h,nmlat_T1,nmlat_T2,nmlon
   public :: nmlonp1,ylatm,ylonm,ylatm_s,ylonm_s,pi,rho,rho_s,rtd,dtr,rearth_m,r0,h0
   public :: m2km,km2m,nhgt_fix,nhgt_fix_r,hgt_fix,hgt_fix_r,ha,ha_s
   public :: nlat_qd,nlat_qd_h,nptsp_total,nptsr_total,nptss1_total,nptss2_total
   public :: nptsp_max,nptsr_max,nptss1_max,nptss2_max

   integer,parameter :: &
!       nmlat_h  = 91,	&	   ! org. 81 pts; a number of magnetic latitudes in one hemisphere P,S1, and R points
! am 10/31/2014 test higher resolution case
       nmlat_h  = 161,  &          ! a test see grid.f90 / number of magnetic latitudes in one hemisphere P,S1, and R points
!*****
!ADR230826: nmlat_h = J of notes.  How are R points defined?
!*****
       nmlatS2_h= nmlat_h-1,	&  ! number of magnetic latitudes in one hemisphere S2 points
       nmlat_T1 = (2*nmlat_h-1),&  ! total number of magnetic latitudes P/S1 points- equator value is double
       nmlat_T2 = (2*nmlatS2_h),&  ! total number of magnetic latitudes S2 points
       nmlat_T3 = (2*nmlat_h)  ,&  ! total number of magnetic latitudes R points- no equator value
       nlat_qd = 321,  &  ! 361; 221; 1181 ; 21,  number of quasi dipole latitudes, edge points j-0.5
!       nlat_qd = 17,  &  !1181 ; 21,           &  ! number of quasi dipole latitudes, edge points j-0.5
!       nlat_qd = 181,  &  !1181 ; 21,          &  ! number of quasi dipole latitudes, edge points j-0.5
       nlat_qd_h=(nlat_qd+1)*0.5, &! half of the hemisphere (assumes point at equator) edge points j-0.5

       nmlon	= 180,  	&  ! number of magnetic longitudes P,S1,S2,R points
       nmlonp1  = nmlon+1,	&
       nhgt_fix   = 82 , 	&  ! Number of height levels on which P points lie.
       nhgt_fix_r = nhgt_fix+1     ! Number of height levels encompassing the lower and upper faces of elemental volumes

   real(r8), parameter ::	      &
       dtr = pi/180._r8,	      & ! Conversion factor when going from degrees to radians
       rtd = 180._r8/pi	                ! Conversion factor when going from radians to degrees

   real(r8) :: ylonm(nmlon),	      & ! magnetic longitudes of p and s2 grid; same for both hemispheres
	       ylatm(nmlat_h,2),      & ! ylatm_s(j) is the magnetic latitude of the equatorward/lower face of an elemental volume for which the P point is at latitude ylatm(j).
	       ylatm_s(nmlatS2_h,2),  & ! magnetic latitudes of s2-points; 2 index for hemisphere
	       ylonm_s(nmlon),        & ! magnetic longitude of the eastern face of an elemental volume for which the P point is at longitude ylonm(i)
	       rho(nmlat_h,2),        & ! cos of magnetic latitudes of p and s1 points; 2 index for hemisphere
	       rho_s(nmlatS2_h,2),    & ! cos of magnetic latitudes of s2-points; 2 index for hemisphere
	       hgt_fix(nhgt_fix),     & ! array with fixed heights for s&p-grids
	       hgt_fix_r(nhgt_fix_r), & ! height of the lower face of an elemental volume for which the P point is at height hgt_fix(k)
	       ha(nmlat_h),	      & ! apex height calculated in grid.f90 for ylatm points; same for both hemispheres
	       ha_s(nmlatS2_h)          ! apex height calculated in grid.f90 for ylatm_s points; same for both hemispheres
!
! global constants
   real, parameter ::	          &
	   h0 = 8.0e4_r8,         &   ! Initial value for bottom height of dynamo grid
	   r0 =rearth_m+h0,	  &   ! Mean Earth radius plus height of bottom of dynamo region (h0) [m]
	   m2km=1.e-3,            &   ! Conversion factor when going from meters to kilometers
	   km2m=1.e3                  ! Conversion factor when going from kilometers to meters
   !
   ! Field point parameters
   !
   integer :: nptsp_max  = 0  !Total number of field points on p grid
   integer :: nptsr_max  = 0  !Total number of field points on 3(r) grid
   integer :: nptss1_max = 0  !Total number of field points on s1 grid
   integer :: nptss2_max = 0  !Total number of field points on s2 grid

   integer :: nptsp_total  = 0  !Total number of field points on p grid
   integer :: nptsr_total  = 0  !Total number of field points on 3(r) grid
   integer :: nptss1_total = 0  !Total number of field points on s1 grid
   integer :: nptss2_total = 0  !Total number of field points on s2 grid


end module edyn3D_params
