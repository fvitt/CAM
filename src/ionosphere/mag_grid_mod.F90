module mag_grid_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst,    only: pi

  implicit none

  integer :: mlon0,mlon1
  integer, parameter :: nmlon = 18

  real(r8) :: gmlon(nmlon) = 0._r8

  integer, parameter :: nmlat = 10
  real(r8) :: gmlat(nmlat) = 0._r8

  real(r8), parameter :: r2d = 180._r8/pi
  real(r8), parameter :: d2r = pi/180._r8

contains

  subroutine mag_grid_mod_reg
    use cam_history,  only: addfld, horiz_only
    use spmd_utils, only: mpicom, iam

    integer :: i

    if (iam==0) then
       mlon0 = 1
       mlon1 = nmlon/2
    elseif (iam==1) then
       mlon0 = (nmlon/2) + 1
       mlon1 = nmlon
    end if

    do i = 1,nmlon
       gmlon(i) = -180._r8 + 360._r8*dble(i-1)/dble(nmlon)
    end do

    do i = 1,nmlat
       gmlat(i) = -90._r8 + 180._r8*dble(i-1)/dble(nmlat-1)
    end do

    call reg_hist_grid()

    call addfld ('MAGTEST', horiz_only, 'I', 'units','Test field' ,gridname='geomag_grid')

  end subroutine mag_grid_mod_reg


  subroutine reg_hist_grid

    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
    use cam_grid_support, only: cam_grid_register


    type(horiz_coord_t), pointer :: lat_coord => null()
    type(horiz_coord_t), pointer :: lon_coord => null()
    integer(iMap),       pointer :: grid_map(:,:) => null()
    integer(iMap),       pointer :: coord_map(:) => null()
    integer                      :: i, j, ind

    integer, parameter :: mag_decomp = 121 ! Must be unique within CAM

    nullify(lat_coord)
    nullify(lon_coord)
    nullify(grid_map)
    nullify(coord_map)

    allocate(grid_map(4, ((mlon1-mlon0+1) * nmlat)))
    grid_map = -huge(1_iMap)

    ind = 0
    do i = 1,nmlat
       do j = mlon0,mlon1
          ind = ind + 1
          grid_map(1, ind) = j
          grid_map(2, ind) = i
          grid_map(3, ind) = j
          grid_map(4, ind) = i
       end do
    end do

    allocate(coord_map(nmlat))
    if (mlon0==1) then
       coord_map = (/ (i, i = 1, nmlat) /)
    else
       coord_map = 0
    end if

    lat_coord => horiz_coord_create('mlat', '', nmlat, 'latitude', &
         'degrees_north', 1,nmlat, gmlat(1:nmlat), map=coord_map)
    nullify(coord_map)

    allocate(coord_map(mlon1 - mlon0 + 1))

    coord_map = (/ (i, i = mlon0, mlon1) /)

    lon_coord => horiz_coord_create('mlon', '', nmlon, 'longitude', &
         'degrees_east', mlon0, mlon1, gmlon(mlon0:mlon1), map=coord_map)
    nullify(coord_map)

    call cam_grid_register('geomag_grid', mag_decomp, lat_coord, lon_coord, grid_map, unstruct=.false.)

    nullify(grid_map)

  end subroutine reg_hist_grid

  subroutine mag_grid_timestep
    use cam_history,  only: outfld

    real(r8) :: testvals(mlon0:mlon1,nmlat)
    integer :: i,j

    testvals=0._r8

    do j = 1,nmlat
       do i = mlon0,mlon1
          testvals(i,j) = sin(d2r*gmlat(j))*cos(d2r*gmlon(i))
       end do
    end do


    do j = 1,nmlat
       call outfld('MAGTEST', testvals(mlon0:mlon1,j), mlon1-mlon0+1, j)
    end do

  end subroutine mag_grid_timestep
end module mag_grid_mod
