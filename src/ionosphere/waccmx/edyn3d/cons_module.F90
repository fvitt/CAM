module cons_module

  use prec,only:rp

  implicit none

  real(kind=rp),parameter :: &
    re = 6.37122e6_rp, &   ! earth radius (m)
    pi = 4*atan(1.0_rp), &
    rtd = 180/pi, &        ! radians to degrees
    dtr = pi/180, &        ! degrees to radians
    h0 = 8e4_rp, &         ! reference height (m) for dynamo calculations
    r0 = re+h0, &
    ylatm_JT = 45*dtr, &   ! transition latitude where potential becomes symmetric/asymmetric
    fill_value = huge(0.0) ! filling value for uninitialized fields

  logical,parameter :: &
    read_pot = .true., & ! whether potential is used at high latitude
    read_fac = .false., & ! whether FAC is used at high latitude
    diagnostics = .true., & ! whether QD currents are calculated
    esmf = & ! whether to use ESMF regridding or native Fortran interface
#ifdef USE_ESMF
    .true.
#else
    .false.
#endif

  integer :: &
    jlatm_JT, & ! latitude index corresponding to the transition latitude
    ndays, &    ! number of days in this year
    ntime       ! number of time steps in this run
  integer,dimension(12) :: days_in_month ! number of days in each month

! lower boundary condition Je2LB is defined by lower atmosphere model
  real(kind=rp),dimension(:,:,:),allocatable :: J3LB ! R current [A/m2]

endmodule cons_module
