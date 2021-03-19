module land_types

  use shr_kind_mod,   only : r8 => shr_kind_r8, shr_kind_cl
  use cam_logfile,    only : iulog
  use spmd_utils,     only : masterproc
  use dyn_grid,       only : get_dyn_grid_parm, get_horiz_grid_d
  use scamMod,        only : single_column
  use ppgrid,         only : pcols, begchunk, endchunk
  use cam_abortutils, only : endrun
  use ioFileMod,      only : getfil
  use cam_pio_utils,  only : cam_pio_openfile, cam_pio_closefile
  use pio,            only : file_desc_t,var_desc_t, pio_nowrite, pio_inq_dimid
  use pio,            only : pio_inq_dimlen, pio_inq_varid, pio_get_var
  use infnan,         only : nan, assignment(=)

  implicit none
  private

  public :: land_types_init, n_land_type, fraction_landuse

  real(r8), protected, allocatable  :: fraction_landuse(:,:,:)

  integer, parameter :: n_land_type = 11

contains

  subroutine land_types_init( depvel_lnd_file, drydep_srf_file )
    use dycore, only : dycore_is

    character(len=*), intent(in) :: depvel_lnd_file
    character(len=*), intent(in) :: drydep_srf_file

    integer :: nlon_veg, nlat_veg, npft_veg
    integer :: dimid
    integer :: i
    integer :: astat
    integer :: plon, plat
    integer :: ierr

    real(r8), allocatable :: vegetation_map(:,:,:)
    real(r8), allocatable :: work(:,:)
    real(r8), allocatable :: landmask(:,:)
    real(r8), allocatable :: urban(:,:)
    real(r8), allocatable :: lake(:,:)
    real(r8), allocatable :: wetland(:,:)
    real(r8), allocatable :: lon_veg_edge(:)
    real(r8), allocatable :: lat_veg_edge(:)

    type(file_desc_t) :: piofile
    type(var_desc_t) :: vid
    character(len=shr_kind_cl) :: locfn

    allocate( fraction_landuse(pcols,n_land_type, begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate fraction_landuse; error = ',astat
       call endrun
    end if
    fraction_landuse = nan

    if(dycore_is('UNSTRUCTURED') ) then
       call get_landuse_and_soilw_from_file(drydep_srf_file)
    else
       plon = get_dyn_grid_parm('plon')
       plat = get_dyn_grid_parm('plat')

       !---------------------------------------------------------------------------
       ! 	... read landuse map
       !---------------------------------------------------------------------------
       call getfil (depvel_lnd_file, locfn, 0)
       call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
       !---------------------------------------------------------------------------
       ! 	... get the dimensions
       !---------------------------------------------------------------------------
       ierr = pio_inq_dimid( piofile, 'lon', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlon_veg )
       ierr = pio_inq_dimid( piofile, 'lat', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlat_veg )
       ierr = pio_inq_dimid( piofile, 'pft', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, npft_veg )
       !---------------------------------------------------------------------------
       ! 	... allocate arrays
       !---------------------------------------------------------------------------
       allocate( vegetation_map(nlon_veg,nlat_veg,npft_veg), work(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
          call endrun
       end if
       allocate( urban(nlon_veg,nlat_veg), lake(nlon_veg,nlat_veg), &
            landmask(nlon_veg,nlat_veg), wetland(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
          call endrun
       end if
       allocate( lon_veg_edge(nlon_veg+1), lat_veg_edge(nlat_veg+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation lon, lat arrays; error = ',astat
          call endrun
       end if
       !---------------------------------------------------------------------------
       ! 	... read the vegetation map and landmask
       !---------------------------------------------------------------------------
       ierr = pio_inq_varid( piofile, 'PCT_PFT', vid )
       ierr = pio_get_var( piofile, vid, vegetation_map )

       ierr = pio_inq_varid( piofile, 'LANDMASK', vid )
       ierr = pio_get_var( piofile, vid, landmask )

       ierr = pio_inq_varid( piofile, 'PCT_URBAN', vid )
       ierr = pio_get_var( piofile, vid, urban )

       ierr = pio_inq_varid( piofile, 'PCT_LAKE', vid )
       ierr = pio_get_var( piofile, vid, lake )

       ierr = pio_inq_varid( piofile, 'PCT_WETLAND', vid )
       ierr = pio_get_var( piofile, vid, wetland )

       call cam_pio_closefile( piofile )

       !---------------------------------------------------------------------------
       ! scale vegetation, urban, lake, and wetland to fraction
       !---------------------------------------------------------------------------
       vegetation_map(:,:,:) = .01_r8 * vegetation_map(:,:,:)
       wetland(:,:)          = .01_r8 * wetland(:,:)
       lake(:,:)             = .01_r8 * lake(:,:)
       urban(:,:)            = .01_r8 * urban(:,:)
#ifdef DEBUG
       if(masterproc) then
          write(iulog,*) 'minmax vegetation_map ',minval(vegetation_map),maxval(vegetation_map)
          write(iulog,*) 'minmax wetland        ',minval(wetland),maxval(wetland)
          write(iulog,*) 'minmax landmask       ',minval(landmask),maxval(landmask)
       end if
#endif
       !---------------------------------------------------------------------------
       ! 	... define lat-lon of vegetation map (1x1)
       !---------------------------------------------------------------------------
       lat_veg_edge(:) = (/ (-90.0_r8 + (i-1),i=1,nlat_veg+1) /)
       lon_veg_edge(:) = (/ (  0.0_r8 + (i-1),i=1,nlon_veg+1) /)

       !---------------------------------------------------------------------------
       ! 	... regrid to model grid
       !---------------------------------------------------------------------------
       call interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg_edge, &
            lon_veg_edge, landmask, urban, lake, &
            wetland, vegetation_map )

       deallocate( vegetation_map, work, stat=astat )
       deallocate( lon_veg_edge, lat_veg_edge, stat=astat )
       deallocate( landmask, urban, lake, wetland, stat=astat )
    endif  ! Unstructured grid
  end subroutine land_types_init


  !-------------------------------------------------------------------------------------
  subroutine get_landuse_and_soilw_from_file(drydep_srf_file)
    use ncdio_atm, only : infld

    character(len=*), intent(in) :: drydep_srf_file

    logical :: readvar

    type(file_desc_t) :: piofile
    character(len=shr_kind_cl) :: locfn
    logical :: lexist

    call getfil (drydep_srf_file, locfn, 1, lexist)
    if(lexist) then
       call cam_pio_openfile(piofile, locfn, PIO_NOWRITE)

       call infld('fraction_landuse', piofile, 'ncol','class',1,pcols,1,n_land_type, begchunk,endchunk, &
            fraction_landuse, readvar, gridname='physgrid')
       if (.not. readvar) then
          write(iulog,*)'**************************************'
          write(iulog,*)'get_landuse_and_soilw_from_file: INFO:'
          write(iulog,*)' fraction_landuse not read from file: '
          write(iulog,*)' ', trim(locfn)
          write(iulog,*)' setting all values to zero'
          write(iulog,*)'**************************************'
          fraction_landuse = 0._r8
       end if

       call cam_pio_closefile(piofile)
    else
       call endrun('Unstructured grids require drydep_srf_file ')
    end if


  end subroutine get_landuse_and_soilw_from_file

  !-------------------------------------------------------------------------------------
  subroutine interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg_edge, &
                         lon_veg_edge, landmask, urban, lake, &
                         wetland, vegetation_map )

    use mo_constants, only : r2d
    use scamMod, only : latiop,loniop,scmlat,scmlon,scm_cambfb_mode
    use shr_scam_mod  , only: shr_scam_getCloseLatLon  ! Standardized system subroutines
    use cam_initfiles, only: initial_file_get_id
    use dycore, only : dycore_is
    use phys_grid,     only : get_rlat_all_p, get_rlon_all_p, get_ncols_p

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer,  intent(in)         :: plon, plat, nlon_veg, nlat_veg, npft_veg
    real(r8), intent(in)         :: landmask(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: urban(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: lake(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: wetland(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: vegetation_map(nlon_veg,nlat_veg,npft_veg)
    real(r8), intent(in)         :: lon_veg_edge(nlon_veg+1)
    real(r8), intent(in)         :: lat_veg_edge(nlat_veg+1)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8) :: closelat,closelon
    integer :: latidx,lonidx

    integer, parameter              :: veg_ext = 20
    type(file_desc_t), pointer      :: piofile
    integer                         :: i, j, ii, jj, i_ndx, n
    integer, dimension(plon+1)      :: ind_lon
    integer, dimension(plat+1)      :: ind_lat
    real(r8)                        :: total_land
    real(r8), dimension(plon+1)     :: lon_edge
    real(r8), dimension(plat+1)     :: lat_edge
    real(r8)                        :: lat1, lon1
    real(r8)                        :: x1, x2, y1, y2, dx, dy
    real(r8)                        :: area, total_area
    real(r8), dimension(npft_veg+3) :: fraction
    real(r8), dimension(-veg_ext:nlon_veg+veg_ext) :: lon_veg_edge_ext
    integer, dimension(-veg_ext:nlon_veg+veg_ext) :: mapping_ext

    real(r8), allocatable :: lam(:), phi(:)

    logical, parameter :: has_npole = .true.
    integer :: ploniop,platiop
    real(r8) :: tmp_frac_lu(plon,n_land_type,plat)

    real(r8):: rlats(pcols), rlons(pcols)
    integer :: lchnk, ncol, icol
    logical :: found

    if(dycore_is('UNSTRUCTURED') ) then
       call endrun('mo_drydep::interp_map called for UNSTRUCTURED grid')
    endif

    allocate(lam(plon), phi(plat))
    call get_horiz_grid_d(plat, clat_d_out=phi)
    call get_horiz_grid_d(plon, clon_d_out=lam)

    if (single_column) then
       if (scm_cambfb_mode) then
          piofile => initial_file_get_id()
          call shr_scam_getCloseLatLon(piofile%fh,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          ploniop=size(loniop)
          platiop=size(latiop)
       else
          latidx=1
          lonidx=1
          ploniop=1
          platiop=1
       end if

       lon_edge(1) = loniop(lonidx) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d

       if (lonidx.lt.ploniop) then
          lon_edge(2) = loniop(lonidx+1) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d
       else
          lon_edge(2) = lon_edge(1) + (loniop(2) - loniop(1)) * r2d
       end if

       lat_edge(1) = latiop(latidx) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d

       if (latidx.lt.platiop) then
          lat_edge(2) = latiop(latidx+1) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d
       else
          lat_edge(2) = lat_edge(1) + (latiop(2) - latiop(1)) * r2d
       end if
    else
       do i = 1,plon
          lon_edge(i) = lam(i) * r2d - .5_r8*(lam(2) - lam(1)) * r2d
       end do
       lon_edge(plon+1) = lon_edge(plon) + (lam(2) - lam(1)) * r2d
       if( .not. has_npole ) then
          do j = 1,plat+1
             lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
          end do
       else
          do j = 1,plat
             lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
          end do
          lat_edge(plat+1) = lat_edge(plat) + (phi(2) - phi(1)) * r2d
       end if
    end if
    do j = 1,plat+1
       lat_edge(j) = min( lat_edge(j), 90._r8 )
       lat_edge(j) = max( lat_edge(j),-90._r8 )
    end do

    !-------------------------------------------------------------------------------------
    ! wrap around the longitudes
    !-------------------------------------------------------------------------------------
    do i = -veg_ext,0
       lon_veg_edge_ext(i) = lon_veg_edge(nlon_veg+i) - 360._r8
       mapping_ext     (i) =              nlon_veg+i
    end do
    do i = 1,nlon_veg
       lon_veg_edge_ext(i) = lon_veg_edge(i)
       mapping_ext     (i) =              i
    end do
    do i = nlon_veg+1,nlon_veg+veg_ext
       lon_veg_edge_ext(i) = lon_veg_edge(i-nlon_veg) + 360._r8
       mapping_ext     (i) =              i-nlon_veg
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : lon_edge ',lon_edge
    write(iulog,*) 'interp_map : lat_edge ',lat_edge
    write(iulog,*) 'interp_map : mapping_ext ',mapping_ext
#endif
    do j = 1,plon+1
       lon1 = lon_edge(j)
       do i = -veg_ext,nlon_veg+veg_ext
          dx = lon_veg_edge_ext(i  ) - lon1
          dy = lon_veg_edge_ext(i+1) - lon1
          if( dx*dy <= 0._r8 ) then
             ind_lon(j) = i
             exit
          end if
       end do
    end do

    do j = 1,plat+1
       lat1 = lat_edge(j)
       do i = 1,nlat_veg
          dx = lat_veg_edge(i  ) - lat1
          dy = lat_veg_edge(i+1) - lat1
          if( dx*dy <= 0._r8 ) then
             ind_lat(j) = i
             exit
          end if
       end do
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : ind_lon ',ind_lon
    write(iulog,*) 'interp_map : ind_lat ',ind_lat
#endif
    lat_loop : do j = 1,plat
       lon_loop : do i = 1,plon
          total_area       = 0._r8
          fraction         = 0._r8
          do jj = ind_lat(j),ind_lat(j+1)
             y1 = max( lat_edge(j),lat_veg_edge(jj) )
             y2 = min( lat_edge(j+1),lat_veg_edge(jj+1) )
             dy = (y2 - y1)/(lat_veg_edge(jj+1) - lat_veg_edge(jj))
             do ii =ind_lon(i),ind_lon(i+1)
                i_ndx = mapping_ext(ii)
                x1 = max( lon_edge(i),lon_veg_edge_ext(ii) )
                x2 = min( lon_edge(i+1),lon_veg_edge_ext(ii+1) )
                dx = (x2 - x1)/(lon_veg_edge_ext(ii+1) - lon_veg_edge_ext(ii))
                area = dx * dy
                total_area = total_area + area
                !-----------------------------------------------------------------
                ! 	... special case for ocean grid point
                !-----------------------------------------------------------------
                if( nint(landmask(i_ndx,jj)) == 0 ) then
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area
                else
                   do n = 1,npft_veg
                      fraction(n) = fraction(n) + vegetation_map(i_ndx,jj,n) * area
                   end do
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area * lake   (i_ndx,jj)
                   fraction(npft_veg+2) = fraction(npft_veg+2) + area * wetland(i_ndx,jj)
                   fraction(npft_veg+3) = fraction(npft_veg+3) + area * urban  (i_ndx,jj)
                   !-----------------------------------------------------------------
                   ! 	... check if land accounts for the whole area.
                   !           If not, the remaining area is in the ocean
                   !-----------------------------------------------------------------
                   total_land = sum(vegetation_map(i_ndx,jj,:)) &
                              + urban  (i_ndx,jj) &
                              + lake   (i_ndx,jj) &
                              + wetland(i_ndx,jj)
                   if( total_land < 1._r8 ) then
                      fraction(npft_veg+1) = fraction(npft_veg+1) + (1._r8 - total_land) * area
                   end if
                end if
             end do
          end do
          !-------------------------------------------------------------------------------------
          ! 	... divide by total area of grid box
          !-------------------------------------------------------------------------------------
          fraction(:) = fraction(:)/total_area
          !-------------------------------------------------------------------------------------
          ! 	... make sure we don't have too much or too little
          !-------------------------------------------------------------------------------------
          if( abs( sum(fraction) - 1._r8) > .001_r8 ) then
             fraction(:) = fraction(:)/sum(fraction)
          end if
          !-------------------------------------------------------------------------------------
          ! 	... map to Wesely land classification
          !-------------------------------------------------------------------------------------
          tmp_frac_lu(i, 1, j) =     fraction(20)
          tmp_frac_lu(i, 2, j) = sum(fraction(16:17))
          tmp_frac_lu(i, 3, j) = sum(fraction(13:15))
          tmp_frac_lu(i, 4, j) = sum(fraction( 5: 9))
          tmp_frac_lu(i, 5, j) = sum(fraction( 2: 4))
          tmp_frac_lu(i, 6, j) =     fraction(19)
          tmp_frac_lu(i, 7, j) =     fraction(18)
          tmp_frac_lu(i, 8, j) =     fraction( 1)
          tmp_frac_lu(i, 9, j) = 0._r8
          tmp_frac_lu(i,10, j) = 0._r8
          tmp_frac_lu(i,11, j) = sum(fraction(10:12))
       end do lon_loop
    end do lat_loop

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       call get_rlat_all_p(lchnk, ncol, rlats(:ncol))
       call get_rlon_all_p(lchnk, ncol, rlons(:ncol))
       do icol= 1,ncol
          found=.false.
          find_col: do j = 1,plat
             do i = 1,plon
                if (rlats(icol)==phi(j) .and. rlons(icol)==lam(i)) then
                   found=.true.
                   exit find_col
                endif
             enddo
          enddo find_col

          if (.not.found) call endrun('mo_drydep::interp_map not able find physics column coordinate')
          fraction_landuse(icol,1:n_land_type,lchnk) =  tmp_frac_lu(i,1:n_land_type,j)

       end do

       !-------------------------------------------------------------------------------------
       ! 	... make sure there are no out of range values
       !-------------------------------------------------------------------------------------
       where (fraction_landuse(:ncol,:n_land_type,lchnk) < 0._r8) fraction_landuse(:ncol,:n_land_type,lchnk) = 0._r8
       where (fraction_landuse(:ncol,:n_land_type,lchnk) > 1._r8) fraction_landuse(:ncol,:n_land_type,lchnk) = 1._r8
    end do

  end subroutine interp_map

end module land_types
