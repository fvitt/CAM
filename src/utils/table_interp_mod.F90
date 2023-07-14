module table_interp_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none

  private
  public :: table_interp
  public :: table_interp_wghts
  public :: table_interp_setwghts
  public :: table_interp_delwghts

  interface table_interp
     module procedure interp2d
  end interface table_interp

  type :: table_interp_wghts
     real(r8), pointer :: w1(:) => null()
     real(r8), pointer :: w2(:) => null()
     integer, pointer :: ix1(:) => null()
     integer, pointer :: ix2(:) => null()
  end type table_interp_wghts

contains

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  function interp2d( ncoef,ncol,nxs,nys, xwghts,ywghts, tbl ) result(res)

    integer, intent(in)  :: ncoef,ncol,nxs,nys
    real(r8), intent(in) :: tbl(ncoef,nxs,nys)
    type(table_interp_wghts), intent(in) :: xwghts
    type(table_interp_wghts), intent(in) :: ywghts

    real(r8) :: res(ncoef,ncol)

    real(r8) :: fx(ncoef,2)

    integer :: i

    do i = 1,ncol

       fx(:,1) = xwghts%w1(i)*tbl(:,xwghts%ix1(i),ywghts%ix1(i)) &
               + xwghts%w2(i)*tbl(:,xwghts%ix2(i),ywghts%ix1(i))
       fx(:,2) = xwghts%w1(i)*tbl(:,xwghts%ix1(i),ywghts%ix2(i)) &
               + xwghts%w2(i)*tbl(:,xwghts%ix2(i),ywghts%ix2(i))

       res(:,i) = ywghts%w1(i)*fx(:,1) + ywghts%w2(i)*fx(:,2)

    end do


  end function interp2d

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  subroutine table_interp_setwghts( ngrid, xgrid, ncols, xcols, wghts )
    integer,  intent(in) :: ngrid
    real(r8), intent(in) :: xgrid(ngrid)
    integer,  intent(in) :: ncols
    real(r8), intent(in) :: xcols(ncols)
    type(table_interp_wghts), intent(out) :: wghts

    integer :: i, ierr

    character(len=*), parameter :: prefix = 'table_interp_setwghts: '

    allocate(wghts%w1(ncols), stat=ierr)
    if( ierr /= 0 ) then
       call endrun(prefix//'failed to allocate wghts%w1')
    end if
    allocate(wghts%w2(ncols), stat=ierr)
    if( ierr /= 0 ) then
       call endrun(prefix//'failed to allocate wghts%w2')
    end if
    allocate(wghts%ix1(ncols), stat=ierr)
    if( ierr /= 0 ) then
       call endrun(prefix//'failed to allocate wghts%ix1')
    end if
    allocate(wghts%ix2(ncols), stat=ierr)
    if( ierr /= 0 ) then
       call endrun(prefix//'failed to allocate wghts%ix2')
    end if


    do i = 1,ncols
       wghts%ix2(i) = find_index(ngrid,xgrid,xcols(i))
       wghts%ix1(i) = wghts%ix2(i) - 1
       wghts%w1(i) = (xgrid(wghts%ix2(i))-xcols(i)) &
                    /(xgrid(wghts%ix2(i))-xgrid(wghts%ix1(i)))
       wghts%w2(i) = (xcols(i)-xgrid(wghts%ix1(i))) &
                    /(xgrid(wghts%ix2(i))-xgrid(wghts%ix1(i)))
    end do

  end subroutine table_interp_setwghts

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  subroutine table_interp_delwghts( wghts )
    type(table_interp_wghts), intent(inout) :: wghts

    deallocate(wghts%w1)
    nullify(wghts%w1)
    deallocate(wghts%w2)
    nullify(wghts%w2)
    deallocate(wghts%ix1)
    nullify(wghts%ix1)
    deallocate(wghts%ix2)
    nullify(wghts%ix2)

  end subroutine table_interp_delwghts

  ! private methods
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  pure function find_index( nvals, vals, vx ) result(ndx)
    integer,  intent(in) :: nvals
    real(r8), intent(in) :: vals(nvals)
    real(r8), intent(in) :: vx

    integer :: ndx

    find_ndx: do ndx = 1, nvals-1
       if (vals(ndx)>vx) exit find_ndx
    end do find_ndx

  end function find_index

end module table_interp_mod
