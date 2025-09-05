module cons_module

  use prec, only: rp, sp

  implicit none

  real(kind=rp),parameter :: &
    re = 6.37122e6_rp, &   ! earth radius (m)
    pi = 4*atan(1.0_rp), &
    rtd = 180._rp/pi, &        ! radians to degrees
    dtr = pi/180._rp, &        ! degrees to radians
    h0 = 8.e4_rp, &         ! reference height (m) for dynamo calculations
    r0 = re+h0, &
    ylatm_JT = 45._rp*dtr, &   ! transition latitude where potential becomes symmetric/asymmetric
    fill_value = huge(1._sp) ! filling value for uninitialized fields

  integer :: &
    jlatm_JT ! latitude index corresponding to the transition latitude

! lower boundary condition Je2LB is defined by lower atmosphere model
  real(kind=rp),dimension(:,:,:),allocatable :: J3LB ! R current [A/m2]

endmodule cons_module
