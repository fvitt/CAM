!--------------------------------------------------------------------------
! interpolate 2-d, 3-d, and 4-d look up tables
!--------------------------------------------------------------------------
module table_interp_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none

  private
  public :: table_interp

  interface table_interp
     module procedure table_interp_2d
     module procedure table_interp_3d
     module procedure table_interp_4d
  end interface table_interp

contains

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  function table_interp_4d( ws,xs,ys,zs, f, w,x,y,z ) result(g)

    real(r8), intent(in) :: ws(:), xs(:), ys(:), zs(:)
    real(r8), intent(in) :: f(:,:,:,:)
    real(r8), intent(in) :: w,x,y,z

    real(r8) :: g

    real(r8) :: g1
    real(r8) :: g2

    integer :: iz

    iz = find_index( zs, z, dimname='Z' )

    g1 = table_interp_3d( ws,xs,ys, f(:,:,:,iz-1), w,x,y )
    g2 = table_interp_3d( ws,xs,ys, f(:,:,:,iz  ), w,x,y )

    g =  interp1d( zs(iz-1),zs(iz), g1, g2, z )

  end function table_interp_4d

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  function table_interp_3d( ws,xs,ys, f, w,x,y ) result(g)

    real(r8), intent(in) :: ws(:), xs(:), ys(:)
    real(r8), intent(in) :: f(:,:,:)
    real(r8), intent(in) :: w,x,y

    real(r8) :: g
    real(r8) :: g1
    real(r8) :: g2
    integer :: iy

    iy = find_index( ys, y, dimname='Y' )

    g1 = table_interp_2d( ws,xs, f(:,:,iy-1), w,x )
    g2 = table_interp_2d( ws,xs, f(:,:,iy  ), w,x )

    g =  interp1d( ys(iy-1),ys(iy), g1, g2, y )

  end function table_interp_3d

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  function table_interp_2d( ws,xs, f, w,x ) result(g)

    real(r8), intent(in) :: ws(:), xs(:)
    real(r8), intent(in) :: f(:,:)
    real(r8), intent(in) :: w,x

    real(r8) :: g

    real(r8) :: g1
    real(r8) :: g2
    integer :: iw,ix

    character(len=*), parameter :: subname = 'table_interp_2d'

    iw = find_index( ws, w, dimname='W' )
    ix = find_index( xs, x, dimname='X' )

    g1 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix-1),f(iw,ix-1), w)
    g2 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix  ),f(iw,ix  ), w)

    g =  interp1d( xs(ix-1),xs(ix), g1, g2, x )

  end function table_interp_2d

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  pure function interp1d ( x1, x2, y1, y2, x ) result(y)
    real(r8), intent(in) :: x1, x2, y1, y2
    real(r8), intent(in) :: x

    real(r8) :: y

    y = y1 + (x-x1)*(y2-y1)/(x2-x1)

  end function interp1d

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  function find_index( values, vx, dimname ) result(ndx)
    real(r8), intent(in) :: values(:)
    real(r8), intent(in) :: vx
    character(len=*),intent(in) :: dimname

    integer :: ndx

    integer  :: nn, i

    ndx = -huge(1)

    nn = size(values)

    do i = 2,nn
       if (.not.(values(i) > values(i-1))) then
          call endrun('table_interp_mod: '//dimname//' values must be monotonically increasing')
       end if
    end do

    if ((vx < values(1)) .or. (vx > values(nn))) then
       call endrun('table_interp_mod: extrapolation is not allowed in the '//dimname//' dimension')
    end if

    find_ndx: do ndx = 1, nn-1
       if (values(ndx)>vx) exit find_ndx
    end do find_ndx

  end function find_index

end module table_interp_mod
