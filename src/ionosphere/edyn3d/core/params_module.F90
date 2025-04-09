module params_module
! grid parameters (geo/mag)

  use prec,only:rp

  implicit none

! geographic grid parameters:
  integer :: nlon,nlat,nlev,nlevp1
  real(kind=rp),dimension(:),allocatable :: glon,glat,zpint

! magnetic grid parameters:
  integer,parameter :: nmlon = 180 ! number of magnetic longitudes - P,S1,S2,R

  integer :: & ! number of magnetic latitudes and fixed heights
    nmlat_h,nmlat_T1, &   ! P,S1,R
    nmlatS2_h,nmlat_T2, & ! S2
    nhgt_fix, &           ! P,S1,S2
    nhgt_fix_r            ! R

! magnetic longitudes including halo points
  real(kind=rp),dimension(0:nmlon+1) :: &
    ylonm, & ! P,S2,R
    ylonm_s  ! S1

! magnetic latitudes separated in two hemispheres
  real(kind=rp),dimension(:,:),allocatable :: &
    ylatm, & ! P,S1,R
    ylatm_s  ! S2

  real(kind=rp),dimension(:),allocatable :: &
    rho,ha, &     ! cos(ylatm) and apex heights - P,S1,R
    rho_s,ha_s, & ! cos(ylatm_s) and apex heights - S2
    hgt_fix, &    ! fixed heights - P,S1,S2
    hgt_fix_r     ! fixed heights - R

endmodule params_module
