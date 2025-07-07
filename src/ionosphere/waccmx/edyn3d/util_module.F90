module util_module

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  pure function find_index(group,instance) result(idx)
! find the index of a string in another string array

    character(len=*),dimension(:),intent(in) :: group
    character(len=*),intent(in) :: instance
    integer :: idx

    do idx = 1,size(group)
      if (trim(instance) == trim(group(idx))) return
    enddo
    idx = 0

  endfunction find_index
!-----------------------------------------------------------------------
  pure function to_month_day(days_in_month,doy) result(month_day)
! convert doy to month,day pair

    integer,dimension(12),intent(in) :: days_in_month
    integer,intent(in) :: doy
    integer,dimension(2) :: month_day

    integer :: month,day

    day = doy
    do month = 1,12
      if (day > days_in_month(month)) then
        day = day - days_in_month(month)
      else
        exit
      endif
    enddo

    month_day(1) = month
    month_day(2) = day

  endfunction to_month_day
!-----------------------------------------------------------------------
  elemental function isclose(a,b,rtol,atol) result(flag)
! check whether a and b are numerically close

    real(kind=rp),intent(in) :: a,b
    real(kind=rp),intent(in),optional :: rtol,atol
    logical :: flag

    real(kind=rp) :: reltol,abstol

    if (present(rtol)) then
      reltol = rtol
    else
      reltol = 1e-9_rp
    endif
    if (present(atol)) then
      abstol = atol
    else
      abstol = 0
    endif

    flag = abs(a-b) <= max(reltol * max(abs(a), abs(b)), abstol)

  endfunction isclose
!-----------------------------------------------------------------------
  pure function inside(x,y,x0,y0) result(in)
! given a convex polygon, check if the test point is inside

    real(kind=rp),dimension(:),intent(in) :: x
    real(kind=rp),dimension(size(x)),intent(in) :: y
    real(kind=rp),intent(in) :: x0,y0
    logical :: in

    integer :: n,i,i0,i1
    real(kind=rp),dimension(size(x)) :: orient

    n = size(x)

    do concurrent (i = 1:n)
      i0 = i
      if (i == n) then
        i1 = 1
      else
        i1 = i+1
      endif
      orient(i) = (x(i1)-x(i0))*(y0   -y(i0))- &
                  (x0   -x(i0))*(y(i1)-y(i0))
    enddo

    if (all(orient >= 0) .or. all(orient <= 0)) then
      in = .true.
    else
      in = .false.
    endif

  endfunction inside
!-----------------------------------------------------------------------
  elemental function lamqd_from_apex_coord(re,r0,lam_m,h) result(lam_qd)
! calculate quasi-dipole latitude lam_qd
! from modified apex latitude and height of point
! h: height of point
! lam_m: modified apex latitude of the field line
! lam_qd = +/- acos([(Re+h)/(Re+hr)]^0.5*cos(lam_m)) Richmond 1995 Eq (6.2)

    real(kind=rp),intent(in) :: re,r0,lam_m,h
    real(kind=rp) :: lam_qd

    real(kind=rp) :: fac

    fac = sqrt((re+h)/r0)*cos(lam_m)

! ensure fac was below 1 and will not cause problem with acos
    if (abs(fac) > 1) fac = sign(1.0_rp,fac)

! lam_qd needs to have the same sign as lam_m
    lam_qd = sign(acos(fac),lam_m)

  endfunction lamqd_from_apex_coord
!-----------------------------------------------------------------------
  pure function find(x,x0) result(i)
! 1D binary search
! x (input): value array, required to be monotonic
! x0 (input): value to be found
! i (output): the index of x0 in x
!   satisfying x(i)<=x0<x(i+1) (if x is increasing)
!   or         x(i)>=x0>x(i+1) (if x is decreasing)
!   0 and size(x) indicate x0 is out of the range of x

! The monotonicity of x is not enforced in this function.
! Caller should ensure x is in order (increasing or decreasing),
! otherwise the result is meaningless.

    real(kind=rp),dimension(:),intent(in) :: x
    real(kind=rp),intent(in) :: x0
    integer :: i

    integer :: nx,i0,i1

    nx = size(x)

! x is in increasing order
    if (x(1) < x(nx)) then
      if (x0 < x(1)) then
        i = 0
      elseif (x0 >= x(nx)) then
        i = nx
      else
        i0 = 1
        i1 = nx
        do while (i0+1 < i1)
          i = (i0+i1)/2
          if (x(i) <= x0) then
            i0 = i
          else
            i1 = i
          endif
        enddo
        i = (i0+i1)/2
      endif

! x is in decreasing order
    else
      if (x0 > x(1)) then
        i = 0
      elseif (x0 <= x(nx)) then
        i = nx
      else
        i0 = 1
        i1 = nx
        do while (i0+1 < i1)
          i = (i0+i1)/2
          if (x(i) >= x0) then
            i0 = i
          else
            i1 = i
          endif
        enddo
        i = (i0+i1)/2
      endif
    endif

  endfunction find
!-----------------------------------------------------------------------
endmodule util_module
