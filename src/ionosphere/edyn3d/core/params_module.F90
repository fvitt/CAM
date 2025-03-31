module params_module
! grid parameters (geo/mag)

  use prec,only:rp

  implicit none

! geographic grid parameters:
  integer :: nlon,nlat,nlev,nlevp1
  real(kind=rp),dimension(:),allocatable :: glon,glat,zpint

! new grid since 2015/02 (see Art's notes 2015/02/12), updated 2015/04
! reference height at k=0.5, grids start at 80 km
! with closer latitude spacing at low latitudes and in the auroral region

! magnetic grid parameters:
  integer,parameter :: &
    nmlon = 180, & ! number of magnetic longitudes - P,S1,S2,R

! number of magnetic latitudes in one hemisphere
    nmlat_h   = 91, &        ! P,S1,R
    nmlatS2_h = nmlat_h-1, & ! S2

! total number of magnetic latitudes
    nmlat_T1 = 2*nmlat_h-1, & ! P,S1,R
    nmlat_T2 = 2*nmlatS2_h, & ! S2

! number of fixed heights
    nhgt_fix   = 54, &      ! P,S1,S2
    nhgt_fix_r = nhgt_fix+1 ! R

! magnetic longitudes including halo points
  real(kind=rp),dimension(0:nmlon+1) :: &
    ylonm, & ! P,S2,R
    ylonm_s  ! S1

! magnetic latitudes, 2 for both hemispheres
  real(kind=rp),dimension(2,nmlat_h  ) :: ylatm   ! P,S1,R
  real(kind=rp),dimension(2,nmlatS2_h) :: ylatm_s ! S2

! cos(lambda), same for both hemispheres - P,S1,R
  real(kind=rp),dimension(nmlat_h) :: &
    rho, & ! cosine of magnetic latitudes
    ha     ! apex heights

! cos(lambda), same for both hemispheres - S2
  real(kind=rp),dimension(nmlatS2_h) :: &
    rho_s, & ! cosine of magnetic latitudes
    ha_s     ! apex heights

! fixed heights
  real(kind=rp),dimension(nhgt_fix  ) :: hgt_fix   ! P,S1,S2
  real(kind=rp),dimension(nhgt_fix_r) :: hgt_fix_r ! R

endmodule params_module
