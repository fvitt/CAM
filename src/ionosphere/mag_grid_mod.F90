module mag_grid_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils, only: masterproc


  use edyn3D_maggrid, only: gen_highres_grid

  use edyn3D_mpi, only: mp_init_edyn3D, mp_distribute_mag_edyn3D, mp_exchange_tasks_edyn3D
  use edyn3D_fieldline, only: fieldline_init, fieldline_getapex
  use edyn3D_params, only: nmlon, nmlat_h, nptsp_total ,ylonm

  use edyn3D_esmf_regrid

  use physconst,     only: pi

  implicit none

  real(r8), parameter :: r2d = 180./pi

contains

  subroutine mag_grid_mod_readnl(nlfile)
    use mo_apex, only: mo_apex_readnl

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    call mo_apex_readnl(nlfile)

  end subroutine mag_grid_mod_readnl

  subroutine mag_grid_mod_reg
    use cam_history,  only: addfld, horiz_only
    use spmd_utils, only: mpicom, iam
    use mo_apex, only: mo_apex_init1

    call mo_apex_init1()

    call mp_init_edyn3D(mpicom)

    call gen_highres_grid()

    call mp_distribute_mag_edyn3D(nmlon)

    call mp_exchange_tasks_edyn3D(mpicom, iprint=0)
    call fieldline_init()  ! Allocate and populate the p, r, s1, and s2 field line structures for computations
    call fieldline_getapex()

    call reg_hist_grid()

    call edyn3D_esmf_regrid_init()

    call addfld ('GEOGALT', horiz_only, 'I', 'km','magnetic field line point geographic (geodetic) altitude', &
                  gridname='mag_fldpnts')
    call addfld ('GEOGLAT', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) latitude', &
                  gridname='mag_fldpnts')
    call addfld ('GEOGLON', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) longitude', &
                  gridname='mag_fldpnts')

  end subroutine mag_grid_mod_reg


  subroutine reg_hist_grid

    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
    use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
    use edyn3d_mpi, only: mlon0_p,mlon1_p
    use edyn3D_fieldline, only: fline_p
    use edyn3d_params, only: nhgt_fix, nmlat_T1

    type(horiz_coord_t), pointer :: flp_coord
    type(horiz_coord_t), pointer :: lon_coord
    integer(iMap),       pointer :: grid_map(:,:)
    integer(iMap),       pointer :: coord_map(:)
    integer                      :: h, i, j, ind

    real(r8), pointer :: latvals(:)
    real(r8), pointer :: altvals(:)

    integer, parameter :: mag_decomp = 703 ! Must be unique within CAM

    real(r8) :: xdel

    integer :: isn, jj, k, ncnt
    integer :: dk,k0,k1


    nullify(flp_coord)
    nullify(lon_coord)
    nullify(grid_map)
    nullify(coord_map)
    nullify(latvals)
    nullify(altvals)

    ncnt = 0
    do j = 1,nmlat_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = fline_p(mlon0_p,j,isn)%npts
             dk = 1
          else
             k0 = fline_p(mlon0_p,j,isn)%npts
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
          end do
       end do
    end do

    nptsp_total = ncnt
    allocate(latvals(nptsp_total))
    allocate(altvals(nptsp_total))

    ncnt = 0
    do j = 1,nmlat_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = fline_p(mlon0_p,j,isn)%npts
             dk = 1
          else
             k0 = fline_p(mlon0_p,j,isn)%npts
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
             latvals(ncnt) = fline_p(mlon0_p,j,isn)%mlat_qd(k)*r2d ! degrees
             altvals(ncnt) = fline_p(mlon0_p,j,isn)%hgt_pt(k)*1.e-3_r8 ! km
          end do
       end do
    end do

    allocate(grid_map(4, ((mlon1_p-mlon0_p+1) * nptsp_total)))
    grid_map = -huge(1_iMap)
    ind = 0
    do i = 1,nptsp_total
       do j = mlon0_p,mlon1_p
          ind = ind + 1
          grid_map(1, ind) = j
          grid_map(2, ind) = i
          grid_map(3, ind) = j
          grid_map(4, ind) = i
       end do
    end do


    allocate(coord_map(nptsp_total))
    if (mlon0_p==1) then
       coord_map = (/ (i, i = 1,nptsp_total ) /)
    else
       coord_map = 0
    end if

    flp_coord => horiz_coord_create('pmlat', 'pflpt', nptsp_total, 'magnetic latitude', &
                                    'degrees_north', 1, nptsp_total, latvals, map=coord_map)
    nullify(coord_map)


    allocate(coord_map(mlon1_p-mlon0_p+1))
    coord_map = (/ (i, i = mlon0_p,mlon1_p ) /)

    lon_coord => horiz_coord_create('pmlon', '', nmlon, 'magnetic longitude', &
                                    'degrees_east', mlon0_p,mlon1_p, r2d*ylonm(mlon0_p:mlon1_p), map=coord_map)
    nullify(coord_map)


    call cam_grid_register('mag_fldpnts', mag_decomp, flp_coord, lon_coord, grid_map, unstruct=.false.)

    nullify(grid_map)
    nullify(flp_coord)
    nullify(lon_coord)

    call cam_grid_attribute_register('mag_fldpnts', 'pmalt', 'magnetic field line point altitude (km)', &
                                     'pflpt', altvals)

    nullify(latvals)
    nullify(altvals)

  end subroutine reg_hist_grid

  subroutine mag_grid_timestep
    use edyn3d_mpi, only: mlon0_p,mlon1_p
    use cam_history,  only: outfld
    use edyn3D_fieldline, only: fline_p

    real(r8) :: geogalts(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglats(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglons(mlon0_p:mlon1_p, nptsp_total)

    integer :: i,j,k,isn,ncnt
    integer :: dk,k0,k1

    geogalts=-huge(1._r8)
    geoglats=-huge(1._r8)
    geoglons=-huge(1._r8)

    do i = mlon0_p,mlon1_p
       ncnt = 0
       do j = 1,nmlat_h
          do isn = 1,2

             if (isn==1) then
                k0 = 1
                k1 = fline_p(i,j,isn)%npts
                dk = 1
             else
                k0 = fline_p(i,j,isn)%npts
                k1 = 1
                dk = -1
             endif

             do k = k0,k1,dk
                ncnt = ncnt + 1
                geogalts(i,ncnt) = fline_p(i,j,isn)%hgt_pt(k)
                geoglats(i,ncnt) = fline_p(i,j,isn)%glat(k)
                geoglons(i,ncnt) = fline_p(i,j,isn)%glon(k)
             end do
          end do
       end do
    end do

    do j = 1,nptsp_total
       call outfld('GEOGALT', geogalts(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       call outfld('GEOGLAT', geoglats(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       call outfld('GEOGLON', geoglons(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
    end do

  end subroutine mag_grid_timestep
end module mag_grid_mod
