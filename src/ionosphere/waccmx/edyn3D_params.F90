module edyn3D_params
  use shr_kind_mod,  only: r8 => shr_kind_r8            ! 8-byte reals
  use shr_const_mod, only: rearth_m => SHR_CONST_REARTH ! meters
  use physconst,     only: pi

  implicit none

  private

  public :: edyn3D_params_init
  public :: edyn3D_params_alloc

  public :: nmlat_h,nmlatS2_h,nmlat_T1,nmlat_T2,nmlon,nlonlat,ylatm_JT,jlatm_JT,phi_pol
  public :: nmlonp1,ylatm,ylonm,ylatm_s,ylonm_s,pi,rho,rho_s,rtd,dtr,rearth_m,r0,h0
  public :: m2km,km2m,nhgt_fix,nhgt_fix_r,hgt_fix,hgt_fix_r,ha,ha_s
  public :: nptsp_total,nptsr_total,nptss1_total,nptss2_total
  public :: nptsp_max,nptsr_max,nptss1_max,nptss2_max
  public :: m1f,m2f,m3f
  public :: J3LB,test_pot,jpg_add,no_wind,use_lbJ,val_fill,use_stabil
  public :: nlonlat_T1
  public :: Je2Ion_eq
  public :: poten_hl

  !
  ! global constants
  !
  real(r8), parameter ::  &
       dtr = pi/180._r8, & ! Conversion factor when going from degrees to radians
       rtd = 180._r8/pi, & ! Conversion factor when going from radians to degrees
       h0 = 8.0e4_r8,    & ! Initial value for bottom height of dynamo grid
       r0 = rearth_m+h0, & ! Mean Earth radius plus height of bottom of dynamo region (h0) [m]
       m2km = 1.e-3_r8,  & ! Conversion factor when going from meters to kilometers
       km2m = 1.e3_r8,   & ! Conversion factor when going from kilometers to meters
       val_fill = 999999._r8  ! fill value

  ! run-time specified resolution parameters
  integer, protected :: nmlat_h = 0
  integer, protected :: nmlats2_h = 0
  integer, protected :: nmlat_T1 = 0
  integer, protected :: nmlat_T2 = 0
  integer, protected :: nmlat_T3 = 0
  integer, protected :: nmlon = 0
  integer, protected :: nmlonp1 = 0
  integer, protected :: nlonlat = 0
  integer, protected :: nlonlat_T1 = 0
  integer, protected :: nhgt_fix = 0
  integer, protected :: nhgt_fix_r = 0

  real(r8) :: &
       ylatm_JT = 45*dtr,     & ! transition latitude where potential becomes symmetric/asymmetric
       jlatm_JT,              & ! latitude index corresponding to the transition latitude
       phi_pol = 0              ! north pole potential

  real(r8), allocatable :: ylonm(:)
  real(r8), allocatable :: ylatm(:,:)
  real(r8), allocatable :: ylonm_s(:)
  real(r8), allocatable :: ylatm_s(:,:)

  real(r8), allocatable :: hgt_fix(:)
  real(r8), allocatable :: hgt_fix_r(:)
  real(r8), allocatable :: ha(:)
  real(r8), allocatable :: ha_s(:)
  real(r8), allocatable :: rho(:,:)
  real(r8), allocatable :: rho_s(:,:)
  !
  ! M1*F, M2*F, M3*F
  !
  real(r8), allocatable :: m1f(:,:)
  real(r8), allocatable :: m2f(:,:)
  real(r8), allocatable :: m3f(:,:)

  real(r8), allocatable :: poten_hl(:,:)
  real(r8), allocatable :: Je2Ion_eq(:)

  !
  ! Boundary conditions
  ! lower boundary Je2LB defined by lower atmosphere model
  !
  logical, parameter :: use_lbJ = .false.
  real(r8), allocatable :: J3LB(:,:,:)

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

contains

  subroutine edyn3D_params_init(edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt)

    integer, intent(in) :: edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt

    nmlat_h = edyn3d_nmlat_h  !  number of magnetic latitudes in one hemisphere P,S1, and R points
    nmlon = edyn3d_nmlon ! number of magnetic longitudes P,S1,S2,R points
    nhgt_fix = edyn3d_nhgt ! Number of height levels on which P points

    nmlatS2_h = nmlat_h-1  ! number of magnetic latitudes in one hemisphere S2 points
    nmlat_T1 = (2*nmlat_h-1) ! total number of magnetic latitudes P/S1 points- equator value is double
    nmlat_T2 = (2*nmlatS2_h) ! total number of magnetic latitudes S2 points
    nmlat_T3 = (2*nmlat_h) ! total number of magnetic latitudes R points- no equator value

    nmlonp1  = nmlon+1
    nlonlat = nmlon*nmlat_h
    nhgt_fix_r = nhgt_fix+1 ! Number of height levels encompassing the lower and upper faces of elemental volumes

    nlonlat_T1 = nmlat_T1*nmlon ! solve the whole globe

    ! light-weight global arrays
    allocate(ylonm(nmlon))
    allocate(ylatm(nmlat_h,2))
    allocate(ylonm_s(nmlon))
    allocate(ylatm_s(nmlatS2_h,2))
    allocate(hgt_fix(nhgt_fix))
    allocate(hgt_fix_r(nhgt_fix_r))
    allocate(ha(nmlat_h))
    allocate(ha_s(nmlatS2_h))
    allocate(rho(nmlat_h,2))
    allocate(rho_s(nmlatS2_h,2))

  end subroutine edyn3D_params_init

  subroutine edyn3D_params_alloc

    ! allocate only on edyn3D active tasks
    allocate(m1f(nmlat_h,nhgt_fix))
    allocate(m2f(nmlat_h,nhgt_fix))
    allocate(m3f(nmlat_h,nhgt_fix_r))

    allocate(poten_hl(0:nmlon+1,nmlat_T1))

  end subroutine edyn3D_params_alloc

end module edyn3D_params
