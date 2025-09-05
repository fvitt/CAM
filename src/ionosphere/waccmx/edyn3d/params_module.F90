module params_module
! grid parameters (geo/mag)

  use prec, only: rp

  implicit none

! magnetic grid parameters:
  integer :: nmlon = 0 ! number of magnetic longitudes - P,S1,S2,R

  integer :: & ! number of magnetic latitudes and fixed heights
    nmlat_h = 0, nmlat_T1 = 0, &   ! P,S1,R
    nmlatS2_h = 0, nmlat_T2 = 0, & ! S2
    nhgt_fix = 0, &                ! P,S1,S2
    nhgt_fix_r = 0                 ! R

! magnetic longitudes including halo points
  real(kind=rp),dimension(:),allocatable :: &
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

  logical :: &
    read_pot = .false., & ! whether potential is used at high latitude
    read_fac = .false.

end module params_module
