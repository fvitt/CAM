module edyn3d_hist_mag_grids_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use fieldline_module, only: npts_s1, npts_s2
  use params_module, only: hgt_fix, nhgt_fix, nmlon, nmlat_h, nmlatS2_h
  use mpi_module, only: mlon0, mlon1, mlat0, mlat1
  use mpi_module, only: mytid=>mpi_rank, ntask=>mpi_size
  use cam_abortutils, only: endrun
  use infnan, only: nan, assignment(=)
  use spmd_utils,      only: masterproc
  use cam_logfile,         only: iulog

  implicit none

  private
  public :: edyn3d_hist_mag_grids_reg
  public :: edyn3d_hist_mag_s1_out
  public :: edyn3d_hist_mag_s2_out

  integer :: s1flpt0=1, s1flpt1=0
  integer :: s2flpt0=1, s2flpt1=0

  integer :: flpts1ndx(nhgt_fix,2,nmlat_h) = -1
  integer :: flpts2ndx(nhgt_fix,2,nmlatS2_h) = -1

  real(r8), parameter :: NOTSET = -huge(1._r8)

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_hist_mag_grids_reg
    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
    use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
    use fieldline_module, only: qdlat_s1, qdlat_s2
    use params_module, only: ylonm, ylonm_s
    use cons_module, only: rtd

    integer, parameter :: mags1_decomp = 719 ! Must be unique within CAM
    integer, parameter :: mags2_decomp = 694
    integer, parameter :: geomag_decomp = 762

    type(horiz_coord_t), pointer :: flns1_coord => null()
    type(horiz_coord_t), pointer :: flns2_coord => null()
    type(horiz_coord_t), pointer :: lons1_coord => null()
    type(horiz_coord_t), pointer :: lons2_coord => null()
    integer(iMap),       pointer :: grid_map(:,:) => null()
    integer(iMap),       pointer :: coord_map(:) => null()

    integer :: ncnt, i, j, k, isn, k0, k1, dk, ind
    integer :: npts1_tot, npts2_tot

    real(r8), pointer :: latvals1(:) => null()
    real(r8), pointer :: altvals1(:) => null()
    real(r8), pointer :: latvals2(:) => null()
    real(r8), pointer :: altvals2(:) => null()
    real(r8) :: lonvals1(nmlon)
    real(r8) :: lonvals2(nmlon)

    character(len=*), parameter :: subname = 'edyn3d_hist_mag_grids_reg'

    if (masterproc) then
       write(iulog,*) subname,'Reg mag fieldline history grid START'
    end if
    if (mytid>=ntask) then
       if (mlon0/=1) then
          call endrun(subname//': mlon0 needs to be 1 on inactive PEs')
       end if
       if (mlat0/=1) then
          call endrun(subname//': mlat0 needs to be 1 on inactive PEs')
       end if
    end if

    ncnt = 0
    do j = 1,nmlat_h
       do k = 1,npts_s1(j)
          ncnt = ncnt + 1
       end do
    end do
    npts1_tot = 2*ncnt

    allocate(latvals1(npts1_tot))
    allocate(altvals1(npts1_tot))
    latvals1 = huge(1._r8)
    altvals1 = -huge(1._r8)

    ncnt = 0
    do j = 1,nmlat_h
       if (j==mlat0) s1flpt0 = ncnt + 1
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = npts_s1(j)
             dk = 1
          else
             k0 = npts_s1(j)
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
             latvals1(ncnt) = qdlat_s1(k,isn,j)*rtd ! degrees
             altvals1(ncnt) = hgt_fix(k)*1.e-3_r8 ! km
             flpts1ndx(k,isn,j) = ncnt
          end do

       end do
       if (j==min(mlat1,nmlat_h)) s1flpt1 = ncnt
    end do

    ncnt = 0
    do j = 1,nmlatS2_h
       do k = 1,npts_s2(j)
          ncnt = ncnt + 1
       end do
    end do
    npts2_tot = 2*ncnt

    allocate(latvals2(npts2_tot))
    allocate(altvals2(npts2_tot))
    latvals2 = huge(1._r8)
    altvals2 = -huge(1._r8)

    ncnt = 0
    do j = 1,nmlatS2_h
       if (j==mlat0) s2flpt0 = ncnt + 1
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = npts_s2(j)
             dk = 1
          else
             k0 = npts_s2(j)
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
             latvals2(ncnt) = qdlat_s2(k,isn,j)*rtd ! degrees
             altvals2(ncnt) = hgt_fix(k)*1.e-3_r8 ! km
             flpts2ndx(k,isn,j) = ncnt
          end do

       end do
       if (j==min(mlat1,nmlatS2_h)) s2flpt1 = ncnt
    end do

    allocate(coord_map(s1flpt1-s1flpt0+1))
    if (mlon0==1) then
       coord_map = (/ (i, i = s1flpt0,s1flpt1 ) /)
    else
       coord_map = 0
    end if
    flns1_coord => horiz_coord_create('lat_s1', 'pflpts1', npts1_tot, 'magnetic latitude', &
                                      'degrees_north', s1flpt0,s1flpt1, latvals1(s1flpt0:s1flpt1), map=coord_map)
    nullify(coord_map)

    allocate(coord_map(s2flpt1-s2flpt0+1))
    if (mlon0==1) then
       coord_map = (/ (i, i = s2flpt0,s2flpt1 ) /)
    else
       coord_map = 0
    end if

    flns2_coord => horiz_coord_create('lat_s2', 'pflpts2', npts2_tot, 'magnetic latitude', &
                                      'degrees_north', s2flpt0,s2flpt1, latvals2(s2flpt0:s2flpt1), map=coord_map)
    nullify(coord_map)

    lonvals1(1:nmlon) = rtd*ylonm_s(1:nmlon)
    lonvals2(1:nmlon) = rtd*ylonm(1:nmlon)

    allocate(coord_map(mlon1-mlon0+1))
    coord_map = (/ (i, i = mlon0,mlon1 ) /)

    lons1_coord => horiz_coord_create('lon_s1', '', nmlon, 'magnetic longitude', &
                                      'degrees_east', mlon0,mlon1, lonvals1(mlon0:mlon1), map=coord_map)
    lons2_coord => horiz_coord_create('lon_s2', '', nmlon, 'magnetic longitude', &
                                      'degrees_east', mlon0,mlon1, lonvals2(mlon0:mlon1), map=coord_map)
    nullify(coord_map)

    allocate(grid_map(4, ((mlon1 - mlon0 + 1) * (s1flpt1 - s1flpt0 + 1))))
    grid_map = -huge(1_iMap)
    ind = 0
    do i = s1flpt0, s1flpt1
       do j = mlon0, mlon1
          ind = ind + 1
          grid_map(1, ind) = j
          grid_map(2, ind) = i
          grid_map(3, ind) = j
          grid_map(4, ind) = i
       end do
    end do

    call cam_grid_register('magfline_s1', mags1_decomp, flns1_coord, lons1_coord, grid_map, unstruct=.false.)

    nullify(grid_map)

    allocate(grid_map(4, ((mlon1 - mlon0 + 1) * (s2flpt1 - s2flpt0 + 1))))
    grid_map = -huge(1_iMap)
    ind = 0
    do i = s2flpt0, s2flpt1
       do j = mlon0, mlon1
          ind = ind + 1
          grid_map(1, ind) = j
          grid_map(2, ind) = i
          grid_map(3, ind) = j
          grid_map(4, ind) = i
       end do
    end do

    call cam_grid_register('magfline_s2', mags2_decomp, flns2_coord, lons2_coord, grid_map, unstruct=.false.)

    nullify(grid_map)

    call cam_grid_attribute_register('magfline_s1', 'alt_s1', 'magnetic field line s1-grid altitude (km)', 'pflpts1', altvals1)
    call cam_grid_attribute_register('magfline_s2', 'alt_s2', 'magnetic field line s2-grid altitude (km)', 'pflpts2', altvals2)

    nullify(latvals1)
    nullify(altvals1)

    nullify(latvals2)
    nullify(altvals2)

    if (masterproc) then
       write(iulog,*) subname,'Reg mag fieldline history grid FINISHED'
    end if

  end subroutine edyn3d_hist_mag_grids_reg

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_hist_mag_s1_out( fldname, fldarray )
    use cam_history, only: outfld

    character(len=*), intent(in) :: fldname
    real(r8), intent(in) :: fldarray(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: tmparray(mlon0:mlon1,s1flpt0:s1flpt1)

    integer :: n, i, j, k, isn, k0, k1, dk

    tmparray = NOTSET

    do j = mlat0,min(mlat1,nmlat_h)
       do isn = 1,2
          if (isn==1) then
             k0 = 1
             k1 = npts_s1(j)
             dk = 1
          else
             k0 = npts_s1(j)
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             n = flpts1ndx(k,isn,j)
             tmparray(:,n)  = fldarray(k,isn,j,:)
          end do
       end do
    end do

    if (any(tmparray==NOTSET)) then
       call endrun('edyn3d_hist_mag_s1_out: tmparray not set correctly')
    end if

    do j = s1flpt0,s1flpt1
       call outfld(fldname, tmparray(mlon0:mlon1,j), mlon1-mlon0+1, j)
    end do

  end subroutine edyn3d_hist_mag_s1_out


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_hist_mag_s2_out( fldname, fldarray )
    use cam_history, only: outfld

    character(len=*), intent(in) :: fldname
    real(r8), intent(in) :: fldarray(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: tmparray(mlon0:mlon1,s2flpt0:s2flpt1)

    integer :: n, i, j, k, isn, k0, k1, dk

    tmparray = NOTSET

    do j = mlat0,min(mlat1,nmlatS2_h)
       do isn = 1,2
          if (isn==1) then
             k0 = 1
             k1 = npts_s2(j)
             dk = 1
          else
             k0 = npts_s2(j)
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             n = flpts2ndx(k,isn,j)
             tmparray(:,n)  = fldarray(k,isn,j,:)
          end do
       end do
    end do

    if (any(tmparray==NOTSET)) then
       call endrun('edyn3d_hist_mag_s2_out: tmparray not set correctly')
    end if

    do j = s2flpt0,s2flpt1
       call outfld(fldname, tmparray(mlon0:mlon1,j), mlon1-mlon0+1, j)
    end do

  end subroutine edyn3d_hist_mag_s2_out


end module edyn3d_hist_mag_grids_mod
