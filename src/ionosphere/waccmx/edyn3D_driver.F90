module edyn3D_driver
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils, only: masterproc

  use edyn3D_maggrid, only: gen_highres_grid

  use edyn3D_mpi, only: mp_init_edyn3D, mp_distribute_mag_edyn3D, mp_exchange_tasks_edyn3D
  use edyn3D_mpi, only: mytid, ntask
  use edyn3D_fieldline, only: fieldline_init, fieldline_getapex
  use edyn3D_params, only: nmlon, nmlat_h, nptsp_total,nptss1_total,nptss2_total, ylonm,ylonm_s

  use edyn3D_fline_fields
  use edyn3D_esmf_regrid
  use edyn3D_regridder

  use physconst, only: pi

  implicit none

  private
  public :: edyn3D_driver_reg
  public :: edyn3D_driver_timestep

  real(r8), parameter :: r2d = 180./pi

contains
  subroutine edyn3D_driver_reg(mpicom, npes)
    use cam_history,  only: addfld, horiz_only
    use mo_apex, only: mo_apex_init1

    integer, intent(in) :: mpicom, npes

    call mo_apex_init1()

    call mp_init_edyn3D(mpicom, npes)

    call gen_highres_grid()

    call mp_distribute_mag_edyn3D(nmlon)

    call mp_exchange_tasks_edyn3D(mpicom, iprint=0)
    call fieldline_init()  ! Allocate and populate the p, r, s1, and s2 field line structures for computations
    call fieldline_getapex()

    call reg_hist_grid()

    call edyn3D_esmf_regrid_init()

    call addfld ('Tn_mag', horiz_only, 'I', 'K','Neutral Temperature on geo-magnetic field line grid', &
                  gridname='magfline_p')

    call addfld ('GEOGALTp', horiz_only, 'I', 'km','magnetic field line point geographic (geodetic) altitude', &
                  gridname='magfline_p')
    call addfld ('GEOGLATp', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) latitude', &
                  gridname='magfline_p')
    call addfld ('GEOGLONp', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) longitude', &
                  gridname='magfline_p')

    call addfld ('GEOGALTs1', horiz_only, 'I', 'km','magnetic field line point geographic (geodetic) altitude', &
                  gridname='magfline_s1')
    call addfld ('GEOGLATs1', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) latitude', &
                  gridname='magfline_s1')
    call addfld ('GEOGLONs1', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) longitude', &
                  gridname='magfline_s1')

    call addfld ('GEOGALTs2', horiz_only, 'I', 'km','magnetic field line point geographic (geodetic) altitude', &
                  gridname='magfline_s2')
    call addfld ('GEOGLATs2', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) latitude', &
                  gridname='magfline_s2')
    call addfld ('GEOGLONs2', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) longitude', &
                  gridname='magfline_s2')


    call addfld ('simga_ped_s1', horiz_only, 'I', 'K','Ped cond. on S1 mag field line grid', &
                  gridname='magfline_s1')
    call addfld ('simga_hal_s1', horiz_only, 'I', 'K','Hal cond. on S1 mag field line grid', &
                  gridname='magfline_s1')
    call addfld ('simga_ped_s2', horiz_only, 'I', 'K','Ped cond. on S2 mag field line grid', &
                  gridname='magfline_s2')
    call addfld ('simga_hal_s2', horiz_only, 'I', 'K','Hal cond. on S2 mag field line grid', &
                  gridname='magfline_s2')

    call edyn3D_fline_fields_alloc()

  end subroutine edyn3D_driver_reg

  subroutine edyn3D_driver_timestep( nphyscol, nphyslev, physalt, tn, sigPed, sigHal, tn_out )
    use edyn3d_mpi, only: mlon0_p,mlon1_p
    use cam_history,  only: outfld
    use edyn3D_fieldline, only: fline_p, fline_s1, fline_s2

    integer,  intent(in) :: nphyscol, nphyslev
    real(r8), intent(in) :: physalt(nphyslev,nphyscol)
    real(r8), intent(in) :: tn(nphyslev,nphyscol)
    real(r8), intent(in) :: sigPed(nphyslev,nphyscol)
    real(r8), intent(in) :: sigHal(nphyslev,nphyscol)
    real(r8), intent(out) :: tn_out(nphyslev,nphyscol)

    real(r8) :: geogaltp(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglatp(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglonp(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geogalts1(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglats1(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglons1(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geogalts2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: geoglats2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: geoglons2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: Tn_tmp(mlon0_p:mlon1_p, nptsp_total)

    integer :: i,j,k,isn,ncnt, ier
    integer :: dk,k0,k1

    call edyn3D_regridder_phys2mag(sigPed,physalt,nphyscol,nphyslev,sigma_ped_s1)
    call edyn3D_regridder_phys2mag(sigPed,physalt,nphyscol,nphyslev,sigma_ped_s2)

    call output_fline_field('simga_ped_s1',sigma_ped_s1)
    call output_fline_field('simga_ped_s2',sigma_ped_s2)

    call edyn3D_regridder_phys2mag(sigHal,physalt,nphyscol,nphyslev,sigma_hal_s1)
    call edyn3D_regridder_phys2mag(sigHal,physalt,nphyscol,nphyslev,sigma_hal_s2)

    call output_fline_field('simga_hal_s1',sigma_hal_s1)
    call output_fline_field('simga_hal_s2',sigma_hal_s2)

    if (mytid<ntask) then
       geogaltp=-huge(1._r8)
       geoglatp=-huge(1._r8)
       geoglonp=-huge(1._r8)

       geogalts1=-huge(1._r8)
       geoglats1=-huge(1._r8)
       geoglons1=-huge(1._r8)

       geogalts2=-huge(1._r8)
       geoglats2=-huge(1._r8)
       geoglons2=-huge(1._r8)

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
                   geogaltp(i,ncnt) = fline_p(i,j,isn)%hgt_pt(k)
                   geoglatp(i,ncnt) = fline_p(i,j,isn)%glat(k)
                   geoglonp(i,ncnt) = fline_p(i,j,isn)%glon(k)
                   geogalts1(i,ncnt) = fline_s1(i,j,isn)%hgt_pt(k)
                   geoglats1(i,ncnt) = fline_s1(i,j,isn)%glat(k)
                   geoglons1(i,ncnt) = fline_s1(i,j,isn)%glon(k)
               end do
             end do
          end do
          ncnt = 0
          do j = 1,nmlats2_h
             do isn = 1,2

                if (isn==1) then
                   k0 = 1
                   k1 = fline_s2(i,j,isn)%npts
                   dk = 1
                else
                   k0 = fline_s2(i,j,isn)%npts
                   k1 = 1
                   dk = -1
                endif

                do k = k0,k1,dk
                   ncnt = ncnt + 1
                   geogalts2(i,ncnt) = fline_s2(i,j,isn)%hgt_pt(k)
                   geoglats2(i,ncnt) = fline_s2(i,j,isn)%glat(k)
                   geoglons2(i,ncnt) = fline_s2(i,j,isn)%glon(k)
               end do
             end do
          end do
       end do

       do j = 1,nptsp_total
          call outfld('GEOGALTp', geogaltp(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLATp', geoglatp(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLONp', geoglonp(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGALTs1', geogalts1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLATs1', geoglats1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLONs1', geoglons1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

       do j = 1,nptss2_total
          call outfld('GEOGALTs2', geogalts2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLATs2', geoglats2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLONs2', geoglons2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

    end if

    call edyn3D_regridder_phys2mag(tn,physalt,nphyscol,nphyslev,Tn_p)

    call output_fline_field('Tn_mag', Tn_p)

    call edyn3D_regridder_mag2phys(Tn_p, physalt, nphyscol,nphyslev, tn_out)

  end subroutine edyn3D_driver_timestep

  subroutine reg_hist_grid

    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
    use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
    use edyn3d_mpi, only: mlon0_p,mlon1_p
    use edyn3D_fieldline, only: fline_p, fline_s2
    use edyn3d_params, only: nhgt_fix, nmlat_T1

    type(horiz_coord_t), pointer :: flp_coord
    type(horiz_coord_t), pointer :: fls1_coord
    type(horiz_coord_t), pointer :: fls2_coord
    type(horiz_coord_t), pointer :: lonp_coord
    type(horiz_coord_t), pointer :: lons1_coord
    type(horiz_coord_t), pointer :: lons2_coord
    integer(iMap),       pointer :: grid_map(:,:)
    integer(iMap),       pointer :: coord_map(:)
    integer                      :: h, i, j, ind

    real(r8), pointer :: latvals(:)
    real(r8), pointer :: altvals(:)
    real(r8), pointer :: latvalss2(:)
    real(r8), pointer :: altvalss2(:)

    integer, parameter :: magp_decomp = 701 ! Must be unique within CAM
    integer, parameter :: mags1_decomp = 702 ! Must be unique within CAM
    integer, parameter :: mags2_decomp = 703 ! Must be unique within CAM

    real(r8) :: xdel

    integer :: isn, jj, k, ncnt
    integer :: dk,k0,k1

    if (mytid>=ntask) then
       if (mlon0_p/=1) then
          call endrun('register_grids: mlat0_p needs to be 1 on inactive PEs')
       end if
       if (mlon1_p/=0) then
          call endrun('register_grids: mlat1_p needs to be 0 on inactive PEs')
       end if
    end if

    nullify(flp_coord)
    nullify(fls1_coord)
    nullify(fls2_coord)
    nullify(lonp_coord)
    nullify(lons1_coord)
    nullify(lons2_coord)
    nullify(grid_map)
    nullify(coord_map)
    nullify(latvals)
    nullify(altvals)
    nullify(latvalss2)
    nullify(altvalss2)

 ! p and S1 -grids

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
    latvals = huge(1._r8)
    altvals = -huge(1._r8)

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

    flp_coord => horiz_coord_create('lat_p', 'pflpt', nptsp_total, 'magnetic latitude', &
                                    'degrees_north', 1, nptsp_total, latvals)
    fls1_coord => horiz_coord_create('lat_s1', 'pflpt', nptsp_total, 'magnetic latitude', &
                                    'degrees_north', 1, nptsp_total, latvals)


    allocate(coord_map(mlon1_p-mlon0_p+1))
    coord_map = (/ (i, i = mlon0_p,mlon1_p ) /)

    lonp_coord => horiz_coord_create('lon_p', '', nmlon, 'magnetic longitude', &
                                    'degrees_east', mlon0_p,mlon1_p, r2d*ylonm(mlon0_p:mlon1_p), map=coord_map)
    lons2_coord => horiz_coord_create('lon_s2', '', nmlon, 'magnetic longitude', &
                                    'degrees_east', mlon0_p,mlon1_p, r2d*ylonm(mlon0_p:mlon1_p), map=coord_map)
    lons1_coord => horiz_coord_create('lon_s1', '', nmlon, 'magnetic longitude', &
                                    'degrees_east', mlon0_p,mlon1_p, r2d*ylonm_s(mlon0_p:mlon1_p), map=coord_map)

    nullify(coord_map)

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


    call cam_grid_register('magfline_p', magp_decomp, flp_coord, lonp_coord, grid_map, unstruct=.false.)
    call cam_grid_register('magfline_s1', mags1_decomp, fls1_coord, lons1_coord, grid_map, unstruct=.false.)
    nullify(flp_coord)
    nullify(fls1_coord)
    nullify(lons1_coord)
    nullify(grid_map)


 ! S2-grid


    ncnt = 0
    do j = 1,nmlats2_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = fline_s2(mlon0_p,j,isn)%npts
             dk = 1
          else
             k0 = fline_s2(mlon0_p,j,isn)%npts
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
          end do
       end do
    end do

    nptss2_total = ncnt
    allocate(latvalss2(nptss2_total))
    allocate(altvalss2(nptss2_total))
    latvalss2 = huge(1._r8)
    altvalss2 = -huge(1._r8)

    ncnt = 0
    do j = 1,nmlats2_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = fline_s2(mlon0_p,j,isn)%npts
             dk = 1
          else
             k0 = fline_s2(mlon0_p,j,isn)%npts
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
             latvalss2(ncnt) = fline_s2(mlon0_p,j,isn)%mlat_qd(k)*r2d ! degrees
             altvalss2(ncnt) = fline_s2(mlon0_p,j,isn)%hgt_pt(k)*1.e-3_r8 ! km
          end do
       end do
    end do

    fls2_coord => horiz_coord_create('lat_s2', 's2flpt', nptss2_total, 'magnetic latitude', &
                                    'degrees_north', 1, nptss2_total, latvalss2)

    allocate(grid_map(4, ((mlon1_p-mlon0_p+1) * nptss2_total)))
    grid_map = -huge(1_iMap)
    ind = 0
    do i = 1,nptss2_total
       do j = mlon0_p,mlon1_p
          ind = ind + 1
          grid_map(1, ind) = j
          grid_map(2, ind) = i
          grid_map(3, ind) = j
          grid_map(4, ind) = i
       end do
    end do

    call cam_grid_register('magfline_s2', mags2_decomp, fls2_coord, lons2_coord, grid_map, unstruct=.false.)
    nullify(fls2_coord)
    nullify(lonp_coord)
    nullify(lons1_coord)
    nullify(lons2_coord)
    nullify(grid_map)

    call cam_grid_attribute_register('magfline_p', 'alt_p', 'magnetic field line p-grid altitude (km)', &
                                     'pflpt', altvals)
    call cam_grid_attribute_register('magfline_s1', 'alt_s1', 'magnetic field line s1-grid altitude (km)', &
                                     'pflpt', altvals)

    call cam_grid_attribute_register('magfline_s2', 'alt_s2', 'magnetic field line s1-grid altitude (km)', &
                                     's2flpt', altvalss2)

    nullify(latvals)
    nullify(altvals)

    nullify(latvalss2)
    nullify(altvalss2)
  end subroutine reg_hist_grid


  subroutine output_fline_field( hfld_name, magfld )
    use cam_history,  only: outfld
    character(len=*), intent(in) :: hfld_name
    type(magfield_t), intent(in) :: magfld

    integer :: i,j,k,isn,ncnt
    integer :: dk,k0,k1
    real(r8) :: fld_tmp(magfld%mlon0:magfld%mlon1,magfld%nptstot)

    if (mytid<ntask) then

       do i = magfld%mlon0,magfld%mlon1
          ncnt = 0
          do j = 1,magfld%nmlat_h
             do isn = 1,2

                if (isn==1) then
                   k0 = 1
                   k1 = magfld%flines(i,j,isn)%npts
                   dk = 1
                else
                   k0 = magfld%flines(i,j,isn)%npts
                   k1 = 1
                   dk = -1
                endif

                do k = k0,k1,dk
                   ncnt = ncnt + 1
                   fld_tmp(i,ncnt) = Tn_p%flines(i,j,isn)%fld(k)
                end do
             end do
          end do
       end do

       do j = 1,magfld%nptstot
          call outfld(hfld_name, fld_tmp(magfld%mlon0:magfld%mlon1,j), magfld%mlon1-magfld%mlon0+1, j)
       end do
    end if

  end subroutine output_fline_field

end module edyn3D_driver
