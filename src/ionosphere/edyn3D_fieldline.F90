module edyn3D_fieldline
!
!-------------------------------------------------------------------------------------
! Purpose:
! Define field line p, r, s1, and s2 structures and allocate and populate variables
!-------------------------------------------------------------------------------------
!
  use edyn3d_params, only: dtr,ylatm,ylatm_s,hgt_fix,nhgt_fix,hgt_fix_r,nhgt_fix_r, &
    ylonm,ylonm_s,rho,rho_s,ha,ha_s,nmlon,nmlat_T1
  use shr_kind_mod,  only: r8 => shr_kind_r8            ! 8-byte reals
  use cam_logfile,    only: iulog
  use spmd_utils,     only: masterproc
!
  implicit none
  save
  private
!
   ! Public data
   public :: fline_p
   public :: fline_r
   public :: fline_s1
   public :: fline_s2
   public :: fieldline_init
   public :: fieldline_getapex

!   public :: tasks
   ! Public type
!   public :: array_ptr_type
   ! Public interfaces
!   public :: mp_init_edyn3D
!   public :: mp_scatter_phim
!   public :: mp_mag_foldhem
!   public :: mp_gather_edyn
!   public :: ixfind
!   public :: mp_magpole_3d
!   public :: setpoles
!   public :: mp_gatherlons_f3d
!   public :: mp_scatterlons_f3d
!   public :: mp_exchange_tasks_edyn3D
!   public :: mp_distribute_mag_edyn3D
  !
  ! Define p field line structure
  !
  type fieldline_p
     integer :: npts

     real(r8) :: ha 	! apex height
     real(r8) :: mlat_m	! modified apex latitude
     real(r8) :: mlon_m	! modified apex longitude

     real(r8) :: pot       ! electric potential
     real(r8) :: pot_test  ! electric potential am 1/2015 for testing
     real(r8) :: fac_hl    ! fieldaligned current high latitude

     real(r8), allocatable :: mlat_qd(:)    ! quasi dipole latitude
     real(r8), allocatable :: mlon_qd(:)    ! quasi dipole longitude
     real(r8), allocatable :: hgt_pt(:)     ! height of point

     real(r8), allocatable :: glat(:)	 ! geog. latitude
     real(r8), allocatable :: glon(:)	 ! geog. longitude

     real(r8), allocatable :: D(:)	 !  D
     real(r8), allocatable :: F(:)	 !  F
     real(r8), allocatable :: sinI(:)	 !  sinI
     real(r8), allocatable :: d1k(:)	 !  d1 dot k vector
     real(r8), allocatable :: d2k(:)	 !  d2 dot k vector
     real(r8), allocatable :: M3(:)	 !  M3
     real(r8), allocatable :: S(:)	 !  S
     real(r8), allocatable :: Jr(:)	 !  Jr
     real(r8), allocatable :: I1hor(:)   !  I1horizontal
     real(r8), allocatable :: I2hor(:)   !  I2horizontal

     integer, allocatable :: ngh_pts(:,:)  ! neighboring points lat_index
  end type fieldline_p

  type magfld_t
     real(r8), allocatable :: fld(:)
  end type magfld_t



  !
  ! Define r field line structure
  !
  type fieldline_r
     integer :: npts

     real(r8) :: ha 	! apex height
     real(r8) :: mlat_m	! modified apex latitude
     real(r8) :: mlon_m	! modified apex longitude
     real(r8) :: pot	! electric potential

     real(r8), allocatable :: mlat_qd(:)    ! quasi dipole latitude
     real(r8), allocatable :: mlon_qd(:)    ! quasi dipole longitude
     real(r8), allocatable :: hgt_pt(:)     ! height of point

     real(r8), allocatable :: glat(:)	 ! geog. latitude
     real(r8), allocatable :: glon(:)	 ! geog. longitude
     real(r8), allocatable :: sinI(:)	 ! sinI coefficient
     real(r8), allocatable :: D(:)	 ! D coefficient
     real(r8), allocatable :: F(:)	 ! F factor
     real(r8), allocatable :: M3(:)	 ! M3 coefficient
     real(r8), allocatable :: I3(:)	 ! I3 current

     real(r8), allocatable :: a3(:)	 ! for mapping: A3(l)/A3(nlat_max)
     real(r8), allocatable :: aa3(:)	 ! for mapping from mod.apex to quasi dipole grid A3=Sum_pole_j M3(j)
     real(r8), allocatable :: LI3(:)	 ! for mapping from mod.apex to quasi dipole grid latitudinal integrated I3

     real(r8), allocatable :: je3(:)	 ! je3 current
     real(r8), allocatable :: Jr(:)	 ! Jr current (diagnostic)

     integer, allocatable :: ngh_pts(:,:)  ! neighboring points lat_index
  end type fieldline_r
  !
  ! Define s1 field line structure
  !
  ! fieldline_s1 are the points on which Je1 will be calculated (these are in longitudinal direction from the p point)
  !
  type fieldline_s1
     integer :: npts

     real(r8) :: ha	    ! apex height
     real(r8) :: mlat_m     ! modified apex latitude
     real(r8) :: mlon_m     ! modified apex longitude
     real(r8) :: zigP	    ! Pedersen Conductance
     real(r8) :: zigH	    ! Hall Conductance
     real(r8) :: Ed1	    ! Ed1 electric field
     real(r8) :: Ed2	    ! Ed2 electric field
     real(r8) :: ve1	    ! ExB ion velocity field
     real(r8) :: ve2	    ! ExB ion velocity field

     real(r8), allocatable :: mlat_qd(:)	! quasi dipole latitude
     real(r8), allocatable :: mlon_qd(:)	! quasi dipole longitude
     real(r8), allocatable :: hgt_pt(:)	! height of point

     real(r8), allocatable :: glat(:)    ! geog. latitude
     real(r8), allocatable :: glon(:)    ! geog. longitude

     real(r8), allocatable :: Vmp(:)     ! magnetic fpotential [TM]
     real(r8), allocatable :: Bmag(:)    ! magnetic field magnitude []
     real(r8), allocatable :: sinI(:)    ! local inclination sin I
     real(r8), allocatable :: bo(:,:)    ! magnetic field component up positive
     real(r8), allocatable :: be3(:)     ! Be3
     real(r8), allocatable :: D(:)       !  D
     real(r8), allocatable :: F(:)       !  D
     real(r8), allocatable :: ds(:)      ! distance between points on fieldline
     real(r8), allocatable :: d1(:,:)    ! d1 vector
     real(r8), allocatable :: d2(:,:)    ! d2 vector
     real(r8), allocatable :: d3(:,:)    ! d3 vector
     real(r8), allocatable :: d1d1(:)    ! d1 dot d1 vector
     real(r8), allocatable :: d1d2(:)    ! d1 dot d2 vector
     real(r8), allocatable :: d2d2(:)    ! d2 dot d2 vector
     real(r8), allocatable :: e1g2(:)    ! e1 dot g2 vector
     real(r8), allocatable :: e2g2(:)    ! e2 dot g2 vector
     real(r8), allocatable :: e1k(:)     ! e1 dot k vector
     real(r8), allocatable :: e2k(:)     ! e2 dot k vector
     real(r8), allocatable :: e3(:,:)    ! e3 vector
     real(r8), allocatable :: M1(:)      ! M1
     real(r8), allocatable :: N1p(:)     ! N1p
     real(r8), allocatable :: N1h(:)     ! N1h
     real(r8), allocatable :: Je1D(:)    ! Je1D
     real(r8), allocatable :: Je1Ion(:)  ! Je1_ionosphere
     real(r8), allocatable :: Je2Ion(:)  ! Je2_ionosphere
     real(r8), allocatable :: M1Je1D(:)  ! M1Je1D addition to stabilize close to the equator with Cowling conductivity
     real(r8), allocatable :: sigH(:)    ! Hall conductivity [S/m]
     real(r8), allocatable :: sigP(:)    ! Pedersen conductivity [S/m]
     real(r8), allocatable :: un(:)      ! zonal neutral wind [m/s] (pos. eastward)
     real(r8), allocatable :: vn(:)      ! meridonal neutral wind [m/s] (pos. northward)

     ! diagnostic
     real(r8), allocatable :: Ne(:)      ! electron density [#/m3]
     real(r8), allocatable :: Tei(:)     ! Te+Ti [K]

     real(r8), allocatable :: je1(:)      ! je1 current
     real(r8), allocatable :: I1(:)       ! I1 current eastward current integrated over the lat/hgt surface
     ! diagnostic
     real(r8), allocatable :: I13d_1(:)  ! I1
     real(r8), allocatable :: I13d_2(:)  ! I1
     real(r8), allocatable :: I13d_3(:)  ! I1

     integer, allocatable :: ngh_pts(:,:)  ! neighboring points lat_index
  end type fieldline_s1
  !
  ! Define s2 field line structure
  !
  ! fieldline_s2 are the points on which Je2 will be calculated (these are
  !   in latitudinal direction from the p point
  !
  type fieldline_s2
     integer :: npts

     real(r8) :: ha	    ! apex height
     real(r8) :: mlat_m     ! modified apex latitude
     real(r8) :: mlon_m     ! modified apex longitude
     real(r8) :: zigP	    ! Pedersen Conductance
     real(r8) :: zigH	    ! Hall Conductance
     real(r8) :: Ed1	    ! Ed1 electric field
     real(r8) :: Ed2	    ! Ed2 electric field

     real(r8), allocatable :: mlat_qd(:)	! quasi dipole latitude
     real(r8), allocatable :: mlon_qd(:)	! quasi dipole longitude
     real(r8), allocatable :: hgt_pt(:)	! height of point

     real(r8), allocatable :: glat(:)    ! geog. latitude
     real(r8), allocatable :: glon(:)    ! geog. longitude

     real(r8), allocatable :: Vmp(:)     ! magnetic fpotential [TM]
     real(r8), allocatable :: Bmag(:)    ! magnetic field magnitude [T]
     real(r8), allocatable :: sinI(:)    ! local inclination sin I
     real(r8), allocatable :: bo(:,:)    ! magnetic field component up positive
     real(r8), allocatable :: be3(:)     ! Be3
     real(r8), allocatable :: D(:)       !  D
     real(r8), allocatable :: F(:)       !  F
     real(r8), allocatable :: ds(:)      ! distance between points on fieldline
     real(r8), allocatable :: d1(:,:)    ! d1 vector
     real(r8), allocatable :: d2(:,:)    ! d2 vector
     real(r8), allocatable :: d1d2(:)    ! d1 dot d2 vector
     real(r8), allocatable :: d2d2(:)    ! d2 dot d2 vector
     real(r8), allocatable :: d3(:,:)    ! d3 vector
     real(r8), allocatable :: e1g2(:)    ! e1 dot g2 vector
     real(r8), allocatable :: e2g2(:)    ! e2 dot g2 vector
     real(r8), allocatable :: e1k(:)     ! e1 dot k vector
     real(r8), allocatable :: e2k(:)     ! e2 dot k vector
     real(r8), allocatable :: e3(:,:)    ! e3 vector
     real(r8), allocatable :: M2(:)      ! M1
     real(r8), allocatable :: N2p(:)     ! N1p
     real(r8), allocatable :: N2h(:)     ! N1h
     real(r8), allocatable :: Je2D(:)    ! Je2D
     real(r8), allocatable :: Je2Ion(:)  ! Je2_ionosphere
     real(r8), allocatable :: Je1Ion(:)  ! Je1_ionosphere
     ! di(agnostic
     real(r8), allocatable :: I23d_1(:)  ! I2
     real(r8), allocatable :: I23d_2(:)  ! I2
     real(r8), allocatable :: I23d_3(:)  ! I2
     real(r8), allocatable :: I2oM2(:)   ! I2/M2
     ! diagnostic
     real(r8), allocatable :: d1d1(:)  ! dot (d1,d1)
     real(r8), allocatable :: e1g1(:)  ! dot (e1,g1)
     real(r8), allocatable :: e2g1(:)  ! dot (e2,g1)
     real(r8), allocatable :: bg1(:)   ! dot (bhat,g1)
     real(r8), allocatable :: bg2(:)   ! dot (bhat,g2)
     real(r8), allocatable :: Jf1Dyn(:)! Jf1(wind driven)
     real(r8), allocatable :: Jf2Dyn(:)! Jf2(wind driven)
     real(r8), allocatable :: Jf1Ion2(:)! Jpg2 in f1 direction
     real(r8), allocatable :: Jf2Ion2(:)! (wind driven)

     real(r8), allocatable :: sigH(:)    ! Hall conductivity [S/m]
     real(r8), allocatable :: sigP(:)    ! Pedersen conductivity [S/m]
     real(r8), allocatable :: un(:)      ! zonal neutral wind [m/s] (pos. eastward)
     real(r8), allocatable :: vn(:)      ! meridonal neutral wind [m/s] (pos. northward)
     ! diagnostic
     real(r8), allocatable :: Ne(:)      ! electron density [#/m3]
     real(r8), allocatable :: Tei(:)     ! Te+Ti [K]

     real(r8), allocatable :: je2(:)      ! je2 current
     real(r8), allocatable :: I2(:)       ! I2 current meridional current integrated over the lon/hgt surface

     integer, allocatable :: ngh_pts(:,:)  ! neighboring points lat_index
  end type fieldline_s2
  !
  ! Declare field line structures
  !
  type (fieldline_p), allocatable  :: fline_p(:,:,:)
  type (fieldline_r), allocatable  :: fline_r(:,:,:)
  type (fieldline_s1), allocatable :: fline_s1(:,:,:)
  type (fieldline_s2), allocatable :: fline_s2(:,:,:)

  type hgt_fld
      integer :: npts		   ! # of latitudes still have points at hgt(k)
      integer, allocatable :: ilat(:)  ! latitude index
  end type hgt_fld

  type (hgt_fld), allocatable :: hgt_fl(:)

  real(r8) :: poten_hl(nmlon,nmlat_T1)  ! high latitude potential r-points
  real(r8) :: poten_hl3(nmlon,nmlat_T1) ! high latitude potential r-points

  real,parameter ::   Je2Ion_eq(nmlon)=0.      ! at S1 points at k=1 and j=nmlat_h
                                               ! otherwise could be interpolated?
                                               ! read in or set later

  type(magfld_t), allocatable :: magfld(:,:,:)

  public :: magfld, magfld_t

  contains
!-----------------------------------------------------------------------------
    subroutine fieldline_init

      use edyn3D_params, only: nmlat_h, &       ! For p and r points
                               nmlatS2_h,rtd    ! For s1 and s2 points

      !,hgt_fix_r,ha,ylatm,ylonm, & ! For r points

      integer :: i,j,k,nlat_k,lat_k(nmlat_h),isn,ilon,is,jns,ier
!      integer :: npt_fldline          ! function
!      integer :: npt_fldline_r        ! function
!      real(r8) :: lamqd_from_apex_coord   ! function
!      real(r8) :: apex_height             ! function

      if (masterproc) then
         write(iulog,"(' edyn3D_fieldline, fieldline_init: Allocating fline structures')")
      endif
      if (masterproc) then
         write(iulog,*) ' edyn3Dmpi, mpi_init_edyn3D: Initializing with nmlon, nmlat_h', nmlon, nmlat_h
      endif
       if (masterproc) then
         write(iulog,*) ' edyn3Dmpi, mpi_init_edyn3D: Longitude p grid ylonm ', ylonm
      endif

      allocate(fline_p(nmlon,nmlat_h,2),stat=ier)
      allocate(magfld(nmlon,nmlat_h,2),stat=ier)
!      allocate(fline_r(nmlon,nmlat_h,2),stat=ier)
!      allocate(fline_s1(nmlon,nmlat_h,2),stat=ier)     ! note same points as p-fieldlines
!      allocate(fline_s2(nmlon,nmlatS2_h,2),stat=ier)   ! note one point less than p-fieldlines
      !
      ! Fieldlines p, r, and s1 can be done together since dimensions are the same
      !
      if (masterproc) then
         write(iulog,"(' edyn3D_fieldline, fieldline_init: Allocating p grid structure ')")
      endif
      do isn = 1,2  ! loop over hemisphere
        do j=1,nmlat_h  ! loop over latitudes (pole to equator)
          fline_p(:,j,isn)%ha     = ha(j)                             ! apex_height
          fline_p(:,j,isn)%npts   = npt_fldline(fline_p(1,j,isn)%ha)  ! points on fieldline
          fline_p(:,j,isn)%mlat_m = ylatm(j,isn)                      ! same as magnetic grid
!	   !
!	   ! r points are
!	   !
!	   fline_r(:,j,isn)%ha     = ha(j)				! apex_height on same fieldline as p points
!	   fline_r(:,j,isn)%npts   = npt_fldline_r(fline_r(1,j,isn)%ha) ! points on fieldline
!	   fline_r(:,j,isn)%mlat_m = ylatm(j,isn)			! same as magnetic grid
!	   !
!	   ! s1 points are in between p-points with respect to longitude, but ylatm is same as p points
!	   !
!	   fline_s1(:,j,isn)%ha     = fline_p(:,j,isn)%ha		! apex_height from p-grid
!	   fline_s1(:,j,isn)%npts   = fline_p(:,j,isn)%npts		! points on fieldline from p-grid
!	   fline_s1(:,j,isn)%mlat_m = ylatm(j,isn)			! magnetic latitude  from p-grid
	   !
	   ! Allocate p, r, s1, s2 field line structure variables
	   !
	   do i=1,nmlon ! loop over longitude

	     fline_p(i,j,isn)%mlon_m = ylonm(i)  !

	     allocate(fline_p(i,j,isn)%hgt_pt(fline_p(i,j,isn)%npts))  ! should be independent of longitude
	     allocate(fline_p(i,j,isn)%mlat_qd(fline_p(i,j,isn)%npts)) ! should be independent of longitude
	    allocate(fline_p(i,j,isn)%mlon_qd(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%glon(fline_p(i,j,isn)%npts))
            allocate(fline_p(i,j,isn)%glat(fline_p(i,j,isn)%npts))
            fline_p(i,j,isn)%glon = -huge(1._r8)
            fline_p(i,j,isn)%glat = -huge(1._r8)
	    allocate(fline_p(i,j,isn)%ngh_pts(2,fline_p(i,j,isn)%npts)) ! lat_ind of neighboring point
	    allocate(fline_p(i,j,isn)%D(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%F(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%sinI(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%d1k(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%d2k(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%M3(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%S(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%Jr(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%I1hor(fline_p(i,j,isn)%npts))
            allocate(fline_p(i,j,isn)%I2hor(fline_p(i,j,isn)%npts))

            allocate(magfld(i,j,isn)%fld(fline_p(i,j,isn)%npts))

!	   allocate(fline_p(i,j,isn)%pot(fline_p(i,j,isn)%npts))
!	   allocate(fline_p(i,j,isn)%pot_test(fline_p(i,j,isn)%npts)) ! am 1/2015 for testing

!	     fline_r(i,j,isn)%mlon_m = ylonm(i)  !
!
!	     allocate(fline_r(i,j,isn)%hgt_pt(fline_r(i,j,isn)%npts))  ! should be independent of longitude
!	     allocate(fline_r(i,j,isn)%mlat_qd(fline_r(i,j,isn)%npts)) ! should be independent of longitude
!	     allocate(fline_r(i,j,isn)%mlon_qd(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%glon(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%glat(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%ngh_pts(2,fline_r(i,j,isn)%npts)) ! lat_ind of neighboring point
!	     allocate(fline_r(i,j,isn)%D(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%F(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%sinI(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%M3(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%I3(fline_r(i,j,isn)%npts))
!
!	     allocate(fline_r(i,j,isn)%a3(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%aa3(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%LI3(fline_r(i,j,isn)%npts))
!
!	     allocate(fline_r(i,j,isn)%je3(fline_r(i,j,isn)%npts))
!	     allocate(fline_r(i,j,isn)%Jr(fline_r(i,j,isn)%npts))
!	     !
!	     !  Relationship between P and S1 points for the same index (i,j)
!	     !  P(i,j) then is really S1(i+0.5,j) with j increasing equatorward
!	     !  coefficient is calculated at P points
!	     !
!	     fline_s1(i,j,isn)%mlon_m = ylonm_s(i)  !
!
!	     allocate(fline_s1(i,j,isn)%hgt_pt(fline_s1(i,j,isn)%npts))  ! should be independent of longitude
!	     allocate(fline_s1(i,j,isn)%mlat_qd(fline_s1(i,j,isn)%npts)) ! should be independent of longitude
!	     allocate(fline_s1(i,j,isn)%mlon_qd(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%glon(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%glat(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%Vmp(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%Bmag(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%sinI(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%bo(3,fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%be3(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%D(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%F(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%d1d1(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%d1d2(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%d2d2(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%d1(3,fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%d2(3,fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%d3(3,fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%e1g2(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%e2g2(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%e1k(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%e2k(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%e3(3,fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%M1(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%N1p(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%N1h(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%Je1D(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%M1Je1D(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%Je1Ion(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%Je2Ion(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%sigH(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%sigP(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%un(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%vn(fline_s1(i,j,isn)%npts))
!	     !
!	     ! diagnostic
!	     allocate(fline_s1(i,j,isn)%Ne(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%Tei(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%I13d_1(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%I13d_2(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%I13d_3(fline_s1(i,j,isn)%npts))
!	     !
!	     allocate(fline_s1(i,j,isn)%ngh_pts(2,fline_s1(i,j,isn)%npts)) ! lat_ind of neighboring point
!	     allocate(fline_s1(i,j,isn)%je1(fline_s1(i,j,isn)%npts))
!	     allocate(fline_s1(i,j,isn)%I1(fline_s1(i,j,isn)%npts))
!
	     do k=1,fline_p(i,j,isn)%npts

	      fline_p(i,j,isn)%hgt_pt(k) = hgt_fix(k)  ! [m] assumes ordering goes from bottom of fieldline to top
						      ! fix heights go also from bottome to top
	      fline_p(i,j,isn)%mlon_qd(k) = ylonm(i)   ! independent of latitude and height
	      fline_p(i,j,isn)%mlat_qd(k) = lamqd_from_apex_coord(fline_p(i,j,isn)%mlat_m,hgt_fix(k))	! quasi dipole latitude

	     enddo
!
!	     do k=1,fline_r(i,j,isn)%npts
!
!	       fline_r(i,j,isn)%hgt_pt(k) = hgt_fix_r(k)  ! [m] assumes ordering goes from bottom of fieldline to top
!							! fix heights go also from bottome to top
!	       fline_r(i,j,isn)%mlon_qd(k) = ylonm(i)	! independent of latitude and height
!!		fline_r(i,j,isn)%mlat_qd(k) = lamqd_from_apex_coord(fline_r(i,j,isn)%mlat_m,hgt_fix_r(k))   ! quasi dipole latitude
!!
!	     enddo
!!
!	     do k=1,fline_s1(i,j,isn)%npts
!!
!	       fline_s1(i,j,isn)%hgt_pt(k) = hgt_fix(k)  ! [m] assumes ordering goes from bottom of fieldline to top
!							 ! fix heights go also from bottome to top
!	       fline_s1(i,j,isn)%mlon_qd(k) = ylonm_s(i)   ! independent of latitude and height
!	       fline_s1(i,j,isn)%mlat_qd(k) = fline_p(i,j,isn)%mlat_qd(k)   ! quasi dipole latitude from p-grid
!
!	     enddo
!!
	    enddo   ! end loop longitudes
	  enddo  ! end loop latitude
if (masterproc) then
   write(iulog,"('edyn3d_fieldline, fieldline_init: Done with allocation of r fieldline structure ')")
endif
        !
        !
        !  Relationship between P and S2 points for the same index (i,j)
        !  P(i,j) then is really 2(i,j+0.5) with j increasing equatorward
        !  coefficient is calculated at P points
        !
!	 do j=1,nmlatS2_h    ! loop over latitudes (direction pole to equator)
!			   ! s2 point are inbetween p-points with respect to rho=cos(ylatm),but same mlon as P points
!	   fline_s2(:,j,isn)%ha     = ha_s(j)				! apex_height
!	   fline_s2(:,j,isn)%npts   = npt_fldline(fline_s2(1,j,isn)%ha) ! points on fieldline
!	   fline_s2(:,j,isn)%mlat_m = ylatm_s(j,isn)			! magnetic latitude
!!
!	   do i=1,nmlon ! loop over longitude
!!	if (masterproc) then
!!	   write(iulog,*) 'Looping longitudes allocate s2 fieldline structure fieldline_init:', isn, j, i
!!	endif
!
!	     fline_s2(i,j,isn)%mlon_m = ylonm(i)  !
!
!	     allocate(fline_s2(i,j,isn)%hgt_pt(fline_s2(i,j,isn)%npts))  ! should be independent of longitude
!	     allocate(fline_s2(i,j,isn)%mlat_qd(fline_s2(i,j,isn)%npts)) ! should be independent of longitude
!	     allocate(fline_s2(i,j,isn)%mlon_qd(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%glon(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%glat(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Vmp(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Bmag(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%sinI(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%bo(3,fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%be3(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%D(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%F(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%d1d2(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%d2d2(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%d1(3,fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%d2(3,fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%d3(3,fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%e1g2(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%e2g2(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%e1k(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%e2k(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%e3(3,fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%M2(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%N2p(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%N2h(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Je2D(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Je1Ion(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Je2Ion(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%sigH(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%sigP(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%un(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%vn(fline_s2(i,j,isn)%npts))
!	     ! diagnostic
!	     allocate(fline_s2(i,j,isn)%Ne(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Tei(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%I23d_1(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%I23d_2(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%I23d_3(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%I2oM2(fline_s2(i,j,isn)%npts))
!	     ! diagnostic
!	     allocate(fline_s2(i,j,isn)%d1d1(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%e1g1(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%e2g1(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%bg1(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%bg2(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Jf1Dyn(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Jf2Dyn(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Jf1Ion2(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%Jf2Ion2(fline_s2(i,j,isn)%npts))
!
!	     allocate(fline_s2(i,j,isn)%ngh_pts(2,fline_s2(i,j,isn)%npts)) ! lat_ind of neighboring point
!	     allocate(fline_s2(i,j,isn)%je2(fline_s2(i,j,isn)%npts))
!	     allocate(fline_s2(i,j,isn)%I2(fline_s2(i,j,isn)%npts))
!!
!	     do k=1,fline_s2(i,j,isn)%npts
!!
!	       fline_s2(i,j,isn)%hgt_pt(k)  = hgt_fix(k)  ! [m] assumes ordering goes from bottom of fieldline to top
!							 ! fix heights go also from bottome to top
!	       fline_s2(i,j,isn)%mlon_qd(k) = ylonm(i)   ! independent of latitude and height
!!	if (masterproc .and. isn == 1 .and. j == 12 .and. i == 133) then
!!	   write(iulog,*) 'Before lamqd call for s2 fieldline structure fieldline_init:', isn, j, i, k, fline_s2(i,j,isn)%npts
!!	  write(iulog,*) 'Before lamqd call for s2 fieldline structure fieldline_init:', fline_s2(i,j,isn)%mlat_m, hgt_fix(k)
!!	endif
!!		fline_s2(i,j,isn)%mlat_qd(k) = lamqd_from_apex_coord(fline_s2(i,j,isn)%mlat_m,hgt_fix(k))   ! quasi dipole latitude
!!
!	     enddo
!
!	   enddo   ! end loop longitudes
!	 enddo  ! end loop latitude
      enddo  ! end loop hemisphere
      !
      ! Create list with lat at each fixed height
      !
      allocate(hgt_fl(nhgt_fix))
      i=1
      isn = 1 ! southern hemisphere it will be the same in the northern hemisphere
      do k =1, nhgt_fix  ! assumes each longitide is the same (no loop over longitude)
	  nlat_k = 0
	  do j=1,nmlat_h  ! latitude loop from pole to equator
	    if(fline_p(i,j,isn)%npts >= k) then 	  ! check if #of pts on fldline is => height => intersects
	      nlat_k = nlat_k + 1   ! increase number of latitudinal points at that height k
	      lat_k(nlat_k) = j     ! get latitudinal index
	    endif
	  enddo

	  allocate(hgt_fl(k)%ilat(nlat_k))
	  hgt_fl(k)%npts= nlat_k		     ! number of fieldlines intersecting with that height k
	  hgt_fl(k)%ilat(1:nlat_k)= lat_k(1:nlat_k)  ! latitudinal index of fldline intersecting with that height k
	  !
	  ! Now use the list of latitudes at each height to set the neighboring points for each fieldline point
	  !
	  do j=1,hgt_fl(k)%npts ! set neighboring points for fldlne
	   do is = 1,2
	    do ilon = 1,nmlon
	     if(j==1) then
		fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(1,k)  = -99
		fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(2,k)  = hgt_fl(k)%ilat(j+1)
	     elseif(j ==  hgt_fl(k)%npts) then
		fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(1,k)  = hgt_fl(k)%ilat(j-1)
		fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(2,k)  = -99
	     else
		fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(1,k)  = hgt_fl(k)%ilat(j-1)
		fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(2,k)  = hgt_fl(k)%ilat(j+1)
	    endif
	   enddo  ! end lon loop
	  enddo  ! end is loop
	 enddo  ! end loop field line points
      enddo  ! end loop field line apex heights
!
     contains
!-----------------------------------------------------------------------
      integer function npt_fldline(apex_height)
! calculates number of points along a fieldline
! uses apex_height of that fieldline and fixed height grid
!
      use edyn3D_params, only: nhgt_fix,hgt_fix
      use shr_kind_mod,  only: r8 => shr_kind_r8	    ! 8-byte reals

      real(r8),intent(in) :: apex_height
      integer :: i
!
      i = 1
      npt_fldline =0
      do while (i <= nhgt_fix.and.hgt_fix(i) <= apex_height)  ! round to the nearest number since hgt_fix
        npt_fldline =npt_fldline + 1                                 ! was slightly larger in the last decimal than 90 km
        i = i + 1
        if(i > nhgt_fix) exit
      end do
!
      end function  npt_fldline
!--------------------------------------------------------------------------------------------
      integer function npt_fldline_r(apex_height)
! calculates number of points along a fieldline
! uses apex_height of that fieldline and fixed height grid
!
      use edyn3D_params, only: nhgt_fix_r,hgt_fix_r
      use shr_kind_mod,  only: r8 => shr_kind_r8            ! 8-byte reals

      real(r8),intent(in) :: apex_height
      integer :: i
!
      i = 1
      npt_fldline_r =0
      do while (i <= nhgt_fix_r.and.hgt_fix_r(i) <= apex_height)  ! round to the nearest number since hgt_fix
        npt_fldline_r =npt_fldline_r + 1                                 ! was slightly larger in the last decimal than 90 km
        i = i + 1
        if(i > nhgt_fix_r) exit
      end do
!
      end function  npt_fldline_r
!--------------------------------------------------------------------------------
      function lamqd_from_apex_coord(lat,h) result(lamqd)

      ! calculate quasi-dipole latitude lamqd
      ! from modified apex latitude and height of point
      ! lat modified apex latitude of the fieldline
      ! lamqd = +/- acos([(Re+h)/(Re+hr)]^0.5*cos(lam_m)) eq. (6.2) [Richmond, 1995]

	use edyn3D_params, only: rearth_m,r0
	use shr_kind_mod,  only: r8 => shr_kind_r8	     ! 8-byte reals

	implicit none

	real(r8),intent(in) :: lat,h
	real(r8) :: lamqd

	real(r8) :: fac
!      if (masterproc .and. h > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord rearth_m,r0,h,lat fieldline_init:', rearth_m,r0,h,lat
!      endif

	fac = sqrt((rearth_m+h)/r0)*cos(lat)
!      if (masterproc .and. h > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord fac1 fieldline_init:', fac
!      endif
	!
	! ensure fac was below 1 and will not cause problem with acos
	!
	if (abs(fac) > 1.0_r8) fac = sign(1.0_r8,fac)
!      if (masterproc .and. h > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord fac2 fieldline_init:', fac
!      endif

	lamqd = sign(acos(fac),lat)
!      if (masterproc .and. h > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord lamqd fieldline_init:', lamqd
!      endif

      end function lamqd_from_apex_coord
!
!
!
!
!
!!
!    real function lamqd_from_apex_coord(latm,hgt_fix)
! calculate quasi dipole latitude lamq from modified apex latitude/mod.apex latitude
!	h height of point
!      lamm mod. apex latitude of the fieldline
! lamq= +/- acos([(Re+h)/(Re+hr)]^0.5*cos(lam_m)) eq. (6.2) [Richmond, 1995]
!      use edyn3D_params, only: r0,rearth_m
!      use shr_kind_mod,  only: r8 => shr_kind_r8	    ! 8-byte reals
!!
!      implicit none
!
!      real(r8),intent(in)     :: latm
!      real(r8),intent(in)     :: hgt_fix
!      real(r8) :: fac,tmp
!      real(r8),parameter :: eps = 1.e-6 ! am 10/2014 had to increase for running on glade system/geyser
!
!      if (masterproc .and. hgt_fix > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord rearth_m,r0,hgt_fix fieldline_init:', rearth_m,r0,hgt_fix,latm
!      endif
!      fac = (rearth_m*100+hgt_fix)/r0  ! all units [cm]
!      fac = (rearth_m+hgt_fix)/r0  ! all units [cm]
!      if (masterproc .and. hgt_fix > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord fac1 fieldline_init:', fac
!      endif
!      fac = sqrt(fac)
!      if (masterproc .and. hgt_fix > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord fac2 fieldline_init:', fac
!      endif
!      fac = fac*cos(latm)
!      if (masterproc .and. hgt_fix > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord fac3 fieldline_init:', fac
!      endif
!
!      if (abs(abs(fac)-1.0).lt.eps) fac = 1.0 ! fac was 1.e-15 larger than 1 and cuased problem with acos
! lamq needs to be same sign as lam_m
!      if (masterproc .and. hgt_fix > 1000000) then
!	 write(iulog,*) 'Inside lamqd_from_apex_coord fac4 fieldline_init:', fac
!      endif
!      lamqd_from_apex_coord = sign(acos(fac),latm)
!      if (masterproc .and. hgt_fix > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord acos(fac),lamqd_from_apex_coord fieldline_init:', acos(fac),lamqd_from_apex_coord
!      endif
! set minimum value to eps
!     if(abs(lamqd_from_apex_coord) <= eps)	lamqd_from_apex_coord =  sign(eps,latm)
!
!      if (masterproc .and. hgt_fix > 1000000) then
!	  write(iulog,*) 'Inside lamqd_from_apex_coord eps,lamqd_from_apex_coord fieldline_init:', eps,lamqd_from_apex_coord
!      endif
!      end function  lamqd_from_apex_coord
!-----------------------------------------------------------------------
    end subroutine fieldline_init

    subroutine fieldline_getapex
      use apex, only: apex_q2g, apex_mall
      use edyn3d_params, only: h0, nmlat_h
      use edyn3d_mpi, only: mlon0_p,mlon1_p
      use physconst, only: pi

      integer :: i,j,k,isn, ierr
      real(r8) :: qdlat,qdlon,alt, gdlat,gdlon
      real(r8) :: bmag,si,alon,xlatm,vmp,w,d,be3,sim,xlatqd,f
      real(r8), dimension(3) :: b,bhat,d1,d2,d3,e1,e2,e3,f1,f2

      real(r8), parameter :: href = h0*1e-3_r8
      real(r8), parameter :: r2d = 180._r8/pi

      do isn = 1,2
         do j = 1,nmlat_h
            do i = mlon0_p,mlon1_p
               do k = 1,fline_p(i,j,isn)%npts
                  qdlat = fline_p(i,j,isn)%mlat_qd(k)*r2d ! get quasi-dipole latitude
                  qdlon = fline_p(i,j,isn)%mlon_qd(k)*r2d ! get quasi-dipole longitude
                  alt = fline_p(i,j,isn)%hgt_pt(k)*1e-3 ! convert height from [m] to [km]

                  call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon, ierr)

                  fline_p(i,j,isn)%glon(k) = gdlon
                  fline_p(i,j,isn)%glat(k) = gdlat

                  call apex_mall(gdlat,gdlon,alt,href,b,bhat,bmag,si, &
                       alon,xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, &
                       xlatqd,f,f1,f2, ierr)
                  fline_p(i,j,isn)%sinI(k) = si ! sin(I)
                  fline_p(i,j,isn)%D(k) = d
                  fline_p(i,j,isn)%F(k) = f
                  fline_p(i,j,isn)%d1k(k) = d1(3) ! k: unit vector upward
                  fline_p(i,j,isn)%d2k(k) = d2(3) ! k: unit vector upward
               enddo
            end do
         end do
      end do

    end subroutine fieldline_getapex

end module edyn3D_fieldline
