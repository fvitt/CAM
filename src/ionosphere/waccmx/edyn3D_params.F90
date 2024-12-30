module edyn3D_params
   !
   ! Constants for edyn3D
   !
   use shr_kind_mod,  only: r8 => shr_kind_r8            ! 8-byte reals
   use shr_const_mod, only: rearth_m => SHR_CONST_REARTH ! meters
   use physconst,     only: pi

   implicit none

   private

   public :: nmlat_h,nmlatS2_h,nmlat_T1,nmlat_T2,nmlon,nlonlat
   public :: nmlonp1,ylatm,ylonm,ylatm_s,ylonm_s,pi,rho,rho_s,rtd,dtr,rearth_m,r0,h0
   public :: m2km,km2m,nhgt_fix,nhgt_fix_r,hgt_fix,hgt_fix_r,ha,ha_s
   public :: nlat_qd,nlat_qd_h,nptsp_total,nptsr_total,nptss1_total,nptss2_total
   public :: nptsp_max,nptsr_max,nptss1_max,nptss2_max
   public :: nglon,nglat,glon,glat
   public :: m1f,m2f,m3f
   public :: nggjlon,nggjlat,nggjhgt,ggjlon,ggjclat,ggjhgt,ggjtop
   public :: wts,ktop,k_fix_ggjbot,delBsolution
   public :: J3LB,test_pot,jpg_add,no_wind,use_lbJ,val_fill,use_stabil
   public :: h_leo,mu0,l1ggj,l2ggj,lvhags,lshags,lwork,mdab
   public :: iw,fw,ls,fs,f11gg,f12gg,f21gg,f22gg,f31oFgg,f32oFgg
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
       nlat_qd_h=(nlat_qd+1)/2, &! half of the hemisphere (assumes point at equator) edge points j-0.5

       nmlon	= 180,  	&  ! number of magnetic longitudes P,S1,S2,R points
       nmlonp1  = nmlon+1,	&
       nlonlat = nmlon*nmlat_h, &
       nhgt_fix   = 82 , 	&  ! Number of height levels on which P points lie.
       nhgt_fix_r = nhgt_fix+1, &   ! Number of height levels encompassing the lower and upper faces of elemental volumes
       ! nglon=73, nglat=91 gives 5 x 2 deg lon-lat grid for magnetic perturbations
       nglon =179, nglat = 319     ! dimensions of geographic grid
!       nglon = 141, nglat = 91      ! dimensions of geographic grid
!       nglon = 5, nglat = 5      ! dimensions of geographic grid
!        for calculating delB (cannot be larger than nggjlat)

       !
       ! Select speed and accuracy of delB calculations by the choice of
       !   delBsolution and nggjhgt below.  nggjhgt is the number of current
       !   layers used to get spherical-harmonic coefficients for magnetic
       !   perturbations.
       !
       ! Single current layer combining horizontal currents up to 397.14 km,
       !   with delB calculated only at the ground.
       !     character(len=16), parameter :: delBsolution = 'quick_ground    '
       !     integer,parameter :: nggjhgt=1
       !     real,parameter :: h_LEO=397139.  ! Nominal height for calculating LEO delB
       !
       ! Single current layer combining horizontal currents up to 397.14 km,
       !   with delB calculated at the ground and and at h_LEO (which may not
       !   exceed ggjtop(1)=397139 m for this case).
       !     character(len=16), parameter :: delBsolution = 'quick_ground,LEO'
       !     integer,parameter :: nggjhgt=1
       !     real,parameter :: h_LEO=397139.  ! Nominal height for calculating LEO delB
       !
       ! Horizontal currents are combined into 18 thick layers,
       !   and delB is calculated at the ground and at heights between 299702 m
       !   and 540303 m, and between 764409 m and 948105 m.
       !     character(len=16), parameter :: delBsolution = 'ground,LEO     '
       !     integer,parameter :: nggjhgt=18
 	    real,parameter :: h_LEO=397139.  ! Nominal height for calculating LEO delB
       !      (h_LEO may not exceed ggjtop(nggjhgt)=948105.)
       !
       ! Heights for current layers are set to heights of rho and QD grids.
       !   This gives the most accurate but time-consuming solution.
       character(len=16), parameter :: delBsolution = 'full_hgt_delB  '
       integer,parameter :: nggjhgt=nhgt_fix
!       real,parameter :: h_LEO=400000.  ! Nominal height for calculating LEO delB
!        (h_LEO may not exceed ggjtop(nggjhgt)=1114811.)
!
       integer, parameter :: &
               nggjlon = 360, & ! 120 number of geographic longitude points
!               nggjlon = 12, & ! number of geographic longitude points
!                        for calculating currents, to get spherical-harmonic
!                        coefficients for magnetic perturbations (delB).
!               nggjlat = 120, & ! number of geographic latitude points
!               nggjlat = 16, & ! number of geographic latitude points
               nggjlat = 360, & ! 180  number of geographic latitude points
!                        for calculating currents, to get spherical-harmonic
!                        coefficients for magnetic perturbations.
       ! Following parameters l1ggj through mdab are for SPHEREPACK 3.2
               l1ggj=min0(nggjlat,(nggjlon+2)/2), &
               l2ggj=(nggjlat+1)/2, &
               lvhags=l2ggj*nggjlat*(2*nggjlat-l1ggj+1)+nggjlon+15,  &
               lshags=nggjlat*(3*(l1ggj+l2ggj)-2)+(l1ggj-1)*       &
                      (l2ggj*(2*nggjlat-l1ggj)-3*l1ggj)/2+nggjlon+15, &
               lwork=max0(4*nggjlat*(nggjlat+2)+2, &
                      (2*nggjhgt+1)*nggjlat*nggjlon), &
               mdab=min0(nggjlat,(nggjlon+2)/2)

       real(r8) :: ggjlon(nggjlon)    ! geographic longitude grid for J (radians)
       real(r8) ::  &
           ggjclat(nggjlat), &  ! geographic Gaussian colatitude grid for J (rad)
           wts(nggjlat)         ! Gaussian weights
!
       real(r8) :: ggjhgt(nggjhgt), & ! geographic height grid for J (m)
               ggjtop(nggjhgt)    ! top heights of J layers (m)

       integer ::  ktop         ! k index for the highest hgt_fix layer
!        used for calculations of currents, plus 1.
       integer :: k_fix_ggjbot(nhgt_fix_r) ! array that gives the appropriate
!        index of hgt_fix_r corresponding to the k index of ggjbot(k) (defined
!        in calc_B.f90 such that ggjbot(k+1)=ggjtop(k)).

   real(r8), parameter ::	      &
       dtr = pi/180._r8,	      & ! Conversion factor when going from degrees to radians
       rtd = 180._r8/pi	                ! Conversion factor when going from radians to degrees

   real(r8) :: ylonm(nmlon),      & ! magnetic longitudes of p and s2 grid; same for both hemispheres
	       ylatm(nmlat_h,2),      & ! ylatm_s(j) is the magnetic latitude of the equatorward/lower face of an elemental volume for which the P point is at latitude ylatm(j).
	       ylatm_s(nmlatS2_h,2),  & ! magnetic latitudes of s2-points; 2 index for hemisphere
	       ylonm_s(nmlon),    & ! magnetic longitude of the eastern face of an elemental volume for which the P point is at longitude ylonm(i)
	       rho(nmlat_h,2),        & ! cos of magnetic latitudes of p and s1 points; 2 index for hemisphere
	       rho_s(nmlatS2_h,2),    & ! cos of magnetic latitudes of s2-points; 2 index for hemisphere
	       hgt_fix(nhgt_fix),     & ! array with fixed heights for s&p-grids
	       hgt_fix_r(nhgt_fix_r), & ! height of the lower face of an elemental volume for which the P point is at height hgt_fix(k)
	       ha(nmlat_h),	      & ! apex height calculated in grid.f90 for ylatm points; same for both hemispheres
	       ha_s(nmlatS2_h),       & ! apex height calculated in grid.f90 for ylatm_s points; same for both hemispheres
               glon(nglon),           & ! geographic longitude grid for delB (degrees)
               glat(nglat)              ! geographic latitude grid for delB (degrees)
   !
   ! M1*F, M2*F, M3*F
   !
   real(r8),dimension(nmlat_h,nhgt_fix) :: m1f,m2f
   real(r8),dimension(nmlat_h,nhgt_fix_r) :: m3f
!
! Dimensions below are flipped so use above
!   !
!   ! M1*F, M2*F, M3*F
!   !
!   real(r8),dimension(nhgt_fix,nmlat_h) :: m1f,m2f
!   real(r8),dimension(nhgt_fix_r,nmlat_h) :: m3f
!
! Each gg current point falls in a QD lat-lon box.
! The western edge of the box lies on QD longitude lon_qd_mp(iw(i,j,k)),
!  with iw=1 corresponding to QD longitude = -pi.
! The southern edge of the box lies on QD latitude lat_qd_ed(ls(i,j,k)),
!  with ls=1 corresponding to QD latitude = -pi/2 (south QD pole).
! The fractional distance of the gg point from the western box edge
!  is fw(i,j,k).
! The fractional distance of the gg point from the southern box edge
!  is fs(i,j,k).
!
     real(r8) :: iw(nggjlon,nggjlat,nggjhgt), &
                 fw(nggjlon,nggjlat,nggjhgt), &
                 ls(nggjlon,nggjlat,nggjhgt), &
                 fs(nggjlon,nggjlat,nggjhgt), &
                 f11gg(nggjlon,nggjlat,nggjhgt), &
                 f12gg(nggjlon,nggjlat,nggjhgt), &
                 f21gg(nggjlon,nggjlat,nggjhgt), &
                 f22gg(nggjlon,nggjlat,nggjhgt), &
                 f31oFgg(nggjlon,nggjlat,nggjhgt), &
                 f32oFgg(nggjlon,nggjlat,nggjhgt)
   !
   ! global constants
   !
   real(r8), parameter ::                 &
                   h0 = 8.0e4_r8,         &   ! Initial value for bottom height of dynamo grid
                   r0 =rearth_m+h0,       &   ! Mean Earth radius plus height of bottom of dynamo region (h0) [m]
                   mu0 = 4.e-7*pi,        &
                   m2km=1.e-3_r8,         &   ! Conversion factor when going from meters to kilometers
                   km2m=1.e3_r8,          &   ! Conversion factor when going from kilometers to meters
                   val_fill=999999._r8        ! fill value
    !
    ! Boundary conditions
    ! lower boundary Je2LB defined by lower atmosphere model
    !
    logical, parameter :: use_lbJ = .false.
    real  :: J3LB(nmlon,nmlat_h,2)    ! r-points current from lower atmosphere [A/m2]
    !
    ! Boundary conditions
    ! lower boundary Je2LB defined by lower atmosphere model
    logical, parameter :: test_pot =.false.
    !
    ! No wind forcing
    !
    logical, parameter :: no_wind =.false.
    !
    ! Jpg
    ! Ionospheric current
    !
    logical, parameter :: Jpg     =.false.
    logical, parameter :: Jpg_add =.false.  ! flag can be reomved once Jpg is tested

    logical, parameter :: use_stabil =.false. ! DO NOT CHANGE should be false
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
