module table_interp_mod

  use shr_kind_mod, only: r8 => shr_kind_r8, shr_kind_cl
  use cam_abortutils, only: endrun

  implicit none

  private
  public :: table_interp_4d

contains

  function table_interp_4d( ws,xs,ys,zs, f, w,x,y,z) result(g)

    real(r8), intent(in) :: ws(:), xs(:), ys(:), zs(:)
    real(r8), intent(in) :: f(:,:,:,:)
    real(r8), intent(in) :: w,x,y,z

    real(r8) :: g

    real(r8) :: g1
    real(r8) :: g2
    real(r8) :: g11
    real(r8) :: g22
    real(r8) :: g111
    real(r8) :: g222
    integer :: iw,ix,iy,iz

    character(len=*), parameter :: subname = 'table_interp_4d'

    iw = find_index( ws, w, dimnum=1, dimname='W' )
    ix = find_index( xs, x, dimnum=2, dimname='X' )
    iy = find_index( ys, y, dimnum=3, dimname='Y' )
    iz = find_index( zs, z, dimnum=4, dimname='Z' )

    g1 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix-1,iy-1,iz-1),f(iw,ix-1,iy-1,iz-1), w)
    g2 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix  ,iy-1,iz-1),f(iw,ix  ,iy-1,iz-1), w)
    g11 = interp1d( xs(ix-1),xs(ix), g1, g2, x )

    g1 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix-1,iy  ,iz-1),f(iw,ix-1,iy  ,iz-1), w)
    g2 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix  ,iy  ,iz-1),f(iw,ix  ,iy  ,iz-1), w)
    g22 = interp1d( xs(ix-1),xs(ix), g1, g2, x )

    g111 = interp1d( ys(iy-1),ys(iy), g11, g22, y )

    g1 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix-1,iy-1,iz  ),f(iw,ix-1,iy-1,iz  ), w)
    g2 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix  ,iy-1,iz  ),f(iw,ix  ,iy-1,iz  ), w)
    g11 = interp1d( xs(ix-1),xs(ix), g1, g2, x )

    g1 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix-1,iy  ,iz  ),f(iw,ix-1,iy  ,iz  ), w)
    g2 = interp1d( ws(iw-1),ws(iw), f(iw-1,ix  ,iy  ,iz  ),f(iw,ix  ,iy  ,iz  ), w)
    g22 = interp1d( xs(ix-1),xs(ix), g1, g2, x )

    g222 = interp1d( ys(iy-1),ys(iy), g11, g22, y )

    g = interp1d( zs(iz-1),zs(iz), g111, g222, z )

  contains

    function find_index( values, vx, dimnum, dimname ) result(ndx)
      real(r8), intent(in) :: values(:)
      real(r8), intent(in) :: vx
      integer , intent(in) :: dimnum
      character(len=*),intent(in) :: dimname

      integer :: ndx

      integer  :: nn, i

      ndx = -huge(1)

      nn = size(values)
      if (size(f,dim=dimnum) /= nn) then
         call endrun(subname//': '//dimname//' dimension size of f is incorrect')
      end if

      do i = 2,nn
         if (.not.(values(i) > values(i-1))) then
            call endrun(subname//': '//dimname//' values must be monotonically increasing')
         end if
      end do

      if ((vx < values(1)) .or. (vx >= values(nn))) then
         call endrun(subname//': extrapolation is not allowed in the '//dimname//' dimension')
      end if

      find_ndx: do ndx = 1, nn-1
         if (values(ndx)>vx) exit find_ndx
      end do find_ndx

    end function find_index

    function interp1d ( x1, x2, y1, y2, x ) result(g)
      real(r8), intent(in) :: x1, x2, y1, y2
      real(r8), intent(in) :: x

      real(r8) :: g

      g = y1 + (x-x1)*(y2-y1)/(x2-x1)

    end function interp1d

  end function table_interp_4d


end module table_interp_mod
