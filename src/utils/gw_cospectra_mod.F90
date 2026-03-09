module gw_cospectra_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver, begchunk, endchunk
  use phys_grid, only: get_ncols_p
  use physics_types, only: physics_state
  use esmf_lonlat_grid_mod, only: nlon, nlat, glats, lon_beg, lon_end, lat_beg, lat_end
  use esmf_zonal_fft_mod, only : esmf_zonal_fft_3d, esmf_zonal_fft_init
  use esmf_zonal_mean_mod, only: esmf_zonal_mean_calc
  use esmf_phys2lonlat_mod, only: esmf_phys2lonlat_regrid
  use esmf_lonlat2phys_mod, only: esmf_lonlat2phys_init

  use, intrinsic :: iso_c_binding

  use spmd_utils, only: masterproc
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun

  use pio

  implicit none

  integer :: nftnum = 0
  integer, parameter :: ntime = 8

  real(r8), allocatable :: accum_cospectra_u(:,:,:,:)
  real(r8), allocatable :: accum_cospectra_v(:,:,:,:)

  real(r8), pointer, protected, public :: frcxu_phys(:,:,:) => null() !(pcols,pver,begchunk:endchunk)
  real(r8), pointer, protected, public :: frcyu_phys(:,:,:) => null() !(pcols,pver,begchunk:endchunk)

  type(var_desc_t) :: rest_accum_u_desc, rest_accum_v_desc

contains

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_reg
    use cam_history_support, only: add_hist_coord

    nftnum = nlon/2+1

    ! add history coordinate for FT number
    call add_hist_coord('fft_num', nftnum, 'Fourier Transform Number')

    allocate(accum_cospectra_u( nftnum, lat_beg:lat_end, pver, ntime ))
    accum_cospectra_u = 0._r8
    allocate(accum_cospectra_v( nftnum, lat_beg:lat_end, pver, ntime ))
    accum_cospectra_v = 0._r8

  end subroutine gw_cospectra_reg

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_init()
    use cam_history, only: addfld

    call esmf_zonal_fft_init()

    call addfld('U_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of FT U', gridname='ctem_zm')
    call addfld('U_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of FT U', gridname='ctem_zm')

    call addfld('V_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of FT U', gridname='ctem_zm')
    call addfld('V_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of FT U', gridname='ctem_zm')

    call addfld('OMEGA_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of FT omega', gridname='ctem_zm')
    call addfld('OMEGA_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of FT omega', gridname='ctem_zm')
    call addfld('THETA_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of FT theta', gridname='ctem_zm')
    call addfld('THETA_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of FT theta', gridname='ctem_zm')
    call addfld('WSTAR_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of the complex conjugat of FT omega', gridname='ctem_zm')
    call addfld('WSTAR_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of the complex conjugat of FT omega', gridname='ctem_zm')


    call addfld('UWSTAR_cosp', (/'fft_num','lev    '/), 'I', '1', 'U*WSTAR cospectra', gridname='ctem_zm')
    call addfld('VWSTAR_cosp', (/'fft_num','lev    '/), 'I', '1', 'V*WSTAR cospectra', gridname='ctem_zm')
    call addfld('TWSTAR_cosp', (/'fft_num','lev    '/), 'I', '1', 'THETA*WSTAR cospectra', gridname='ctem_zm')

    call addfld('RHOBAR',  (/'lev'/), 'I', 'kg/m3', 'Zonal mean air density', gridname='ctem_zm')

    call addfld('MFLXXUP', (/'lev'/), 'I', '1', 'Positive unresolved zonal momentum flux', gridname='ctem_zm')
    call addfld('MFLXXUN', (/'lev'/), 'I', '1', 'Negative unresolved zonal momentum flux', gridname='ctem_zm')
    call addfld('MFLXYUP', (/'lev'/), 'I', '1', 'Positive unresolved meridianal momentum flux', gridname='ctem_zm')
    call addfld('MFLXYUN', (/'lev'/), 'I', '1', 'Negative unresolved meridianal momentum flux', gridname='ctem_zm')

    call addfld('FRCXR', (/'lev'/), 'I', '1', 'Resolved zonal momentum flux', gridname='ctem_zm')
    call addfld('FRCXU', (/'lev'/), 'I', '1', 'Unresolved zonal momentum flux', gridname='ctem_zm')
    call addfld('FRCYR', (/'lev'/), 'I', '1', 'Resolved meridianal momentum flux', gridname='ctem_zm')
    call addfld('FRCYU', (/'lev'/), 'I', '1', 'Unresolved meridianal momentum flux', gridname='ctem_zm')

    call addfld('SFRCXR', (/'lev'/), 'I', '1', 'Smoothed Resolved zonal momentum flux', gridname='ctem_zm')
    call addfld('SFRCXU', (/'lev'/), 'I', '1', 'Smoothed Unresolved zonal momentum flux', gridname='ctem_zm')
    call addfld('SFRCYR', (/'lev'/), 'I', '1', 'Smoothed Resolved meridianal momentum flux', gridname='ctem_zm')
    call addfld('SFRCYU', (/'lev'/), 'I', '1', 'Smoothed Unresolved meridianal momentum flux', gridname='ctem_zm')

    call addfld('SFRCXU_phys', (/'lev'/), 'I', 'meters/sec2', 'Smoothed Unresolved zonal forcing', gridname='physgrid')
    call addfld('SFRCYU_phys', (/'lev'/), 'I', 'meters/sec2', 'Smoothed Unresolved meridianal zonal forcing', gridname='physgrid')

    call esmf_lonlat2phys_init()

    allocate(frcxu_phys(pcols,pver,begchunk:endchunk))
    frcxu_phys = 0._r8
    allocate(frcyu_phys(pcols,pver,begchunk:endchunk))
    frcyu_phys = 0._r8

  end subroutine gw_cospectra_init

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_calc(phys_state)
    use cam_history, only: outfld
    use perf_mod, only: t_startf, t_stopf
    use cospext_mod, only: cospext
    use air_composition, only: rairv  ! composition dependent gas constant (J/K/kg)
    use ref_pres, only: pref_mid
    use esmf_phys2lonlat_mod, only: p2l_bdl=>fields_bundle_t, phys2lonlat_nflds=>nflds
    use esmf_lonlat2phys_mod, only: l2p_bdl=>fields_bundle_t, lonlat2phys_nflds=>nflds, esmf_lonlat2phys_regrid

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    complex(C_DOUBLE_COMPLEX) :: u_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: v_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: w_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: t_fft(nftnum, lat_beg:lat_end, pver)

    complex(C_DOUBLE_COMPLEX) :: wstar(nftnum, lat_beg:lat_end, pver)

    real(r8),target :: ufld(pver,pcols,begchunk:endchunk)
    real(r8),target :: vfld(pver,pcols,begchunk:endchunk)
    real(r8),target :: wfld(pver,pcols,begchunk:endchunk)
    real(r8),target :: tfld(pver,pcols,begchunk:endchunk)
    integer :: lchnk, ncol, icol

    real(r8),target :: rho_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: u_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: v_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: w_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: t_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)

    type(p2l_bdl) :: physflds(phys2lonlat_nflds)
    type(p2l_bdl) :: lonlatflds(phys2lonlat_nflds)

    type(l2p_bdl) :: physfrcs(lonlat2phys_nflds)
    type(l2p_bdl) :: lonlatfrcs(lonlat2phys_nflds)

    integer :: i,j,k, n

    complex(r8) :: tmpfld(nftnum, lat_beg:lat_end, pver)
    real(r8) :: cospectra(nftnum, lat_beg:lat_end, pver)
    real(r8) :: latrad(lat_beg:lat_end)
    real(r8),target :: rho(pver,pcols,begchunk:endchunk) ! air mass density
    real(r8) :: rhobar(lat_beg:lat_end, pver)

    real(r8) :: mflxxup(lat_beg:lat_end,pver), mflxxun(lat_beg:lat_end,pver)
    real(r8) :: mflxyup(lat_beg:lat_end,pver), mflxyun(lat_beg:lat_end,pver)
    real(r8) :: frcxr(lat_beg:lat_end,pver), frcxu(lat_beg:lat_end,pver)
    real(r8) :: frcyr(lat_beg:lat_end,pver), frcyu(lat_beg:lat_end,pver)

    real(r8),target :: frcxu_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: frcyu_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)

    real(r8) :: wvlxbeg, wvlxend
    real(r8) :: mflux_glb(nlat,pver)

    real(r8), parameter :: pi = 4._r8*atan(1._r8)
    real(r8), parameter :: deg2rad = pi/180._r8
    character(len=*), parameter :: subname  = 'gw_cospectra_calc'

    wvlxbeg = 200.e3_r8  ! 200 km
    wvlxend = 20.e3_r8   ! 20 km

    call t_startf ('gw_cospectra_calc')

    latrad(lat_beg:lat_end) = glats(lat_beg:lat_end)*deg2rad

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       do icol = 1,ncol
          ufld(:pver,icol,lchnk) = phys_state(lchnk)%u(icol,:pver)
          vfld(:pver,icol,lchnk) = phys_state(lchnk)%v(icol,:pver)
          wfld(:pver,icol,lchnk) = phys_state(lchnk)%omega(icol,:pver)
          tfld(:pver,icol,lchnk) = phys_state(lchnk)%t(icol,:pver) * phys_state(lchnk)%exner(icol,:pver)
          rho (:pver,icol,lchnk) = phys_state(lchnk)%pmid(icol,:pver)/(rairv(icol,:pver,lchnk)*phys_state(lchnk)%t(icol,:pver)) ! kg/m3
       end do
    end do

    physflds(1)%fld => ufld
    physflds(2)%fld => vfld
    physflds(3)%fld => wfld
    physflds(4)%fld => tfld
    physflds(5)%fld => rho

    lonlatflds(1)%fld => u_lonlat
    lonlatflds(2)%fld => v_lonlat
    lonlatflds(3)%fld => w_lonlat
    lonlatflds(4)%fld => t_lonlat
    lonlatflds(5)%fld => rho_lonlat

    call esmf_phys2lonlat_regrid(physflds, lonlatflds)

    u_fft = esmf_zonal_fft_3d(u_lonlat)
    call output_fld(u_fft, name='U')

    v_fft = esmf_zonal_fft_3d(v_lonlat)
    call output_fld(v_fft, name='V')

    w_fft = esmf_zonal_fft_3d(w_lonlat)
    call output_fld(w_fft, name='OMEGA')

    t_fft = esmf_zonal_fft_3d(t_lonlat)
    call output_fld(t_fft, name='THETA')

    wstar = conjg(w_fft)
    call output_fld(wstar, name='WSTAR')

    do i = 1,ntime-1
       accum_cospectra_u(:,:,:,i) = accum_cospectra_u(:,:,:,i+1)
       accum_cospectra_v(:,:,:,i) = accum_cospectra_v(:,:,:,i+1)
    end do

    tmpfld = u_fft * wstar   ! times 2 to account for the other half of the spectrum
    cospectra = tmpfld%re
    cospectra(2:,:,:) = 2._r8 * cospectra(2:,:,:)
    call output_cosp(cospectra,'U')

    accum_cospectra_u(:,:,:,ntime) = cospectra(:,:,:)

    tmpfld = v_fft * wstar
    cospectra = tmpfld%re
    cospectra(2:,:,:) = 2._r8 * cospectra(2:,:,:)
    call output_cosp(cospectra,'V')

    accum_cospectra_v(:,:,:,ntime) = cospectra(:,:,:)

    tmpfld = t_fft * wstar
    cospectra = tmpfld%re
    cospectra(2:,:,:) = 2._r8 * cospectra(2:,:,:)
    call output_cosp(cospectra,'T')

    call esmf_zonal_mean_calc(rho_lonlat, rhobar)

    do icol = lat_beg, lat_end
       call outfld('RHOBAR', rhobar(icol,:),1,icol)
    end do

    ! zonal component
    call cospext(nftnum, lat_beg,lat_end, pver,ntime, latrad, accum_cospectra_u, wvlxbeg,wvlxend, pref_mid, rhobar, mflxxup,mflxxun, frcxr,frcxu)

    ! meridianal component
    call cospext(nftnum, lat_beg,lat_end, pver,ntime, latrad, accum_cospectra_v, wvlxbeg,wvlxend, pref_mid, rhobar, mflxyup,mflxyun, frcyr,frcyu)

    do icol = lat_beg, lat_end
       call outfld('MFLXXUP', mflxxup(icol,:),1,icol)
       call outfld('MFLXXUN', mflxxun(icol,:),1,icol)
       call outfld('MFLXYUP', mflxyup(icol,:),1,icol)
       call outfld('MFLXYUN', mflxyun(icol,:),1,icol)

       call outfld('FRCXR', frcxr(icol,:),1,icol)
       call outfld('FRCXU', frcxu(icol,:),1,icol)
       call outfld('FRCYR', frcyr(icol,:),1,icol)
       call outfld('FRCYU', frcyu(icol,:),1,icol)

    end do

    ! gather and smooth

    frcxr = gather_and_smooth(frcxr)
    frcxu = gather_and_smooth(frcxu)
    frcyr = gather_and_smooth(frcyr)
    frcyu = gather_and_smooth(frcyu)

    do icol = lat_beg, lat_end

       call outfld('SFRCXR', frcxr(icol,:),1,icol)
       call outfld('SFRCXU', frcxu(icol,:),1,icol)
       call outfld('SFRCYR', frcyr(icol,:),1,icol)
       call outfld('SFRCYU', frcyu(icol,:),1,icol)

    end do

    frcxu_lonlat = -huge(1._r8)
    frcyu_lonlat = -huge(1._r8)

    do k = 1,pver
       do j = lat_beg,lat_end
          frcxu_lonlat(lon_beg:lon_end,j,k) = frcxu(j,k)
          frcyu_lonlat(lon_beg:lon_end,j,k) = frcyu(j,k)
       end do
    end do

    frcxu_phys = -huge(1._r8)
    frcyu_phys = -huge(1._r8)

    lonlatfrcs(1)%fld => frcxu_lonlat
    lonlatfrcs(2)%fld => frcyu_lonlat
    physfrcs(1)%fld => frcxu_phys
    physfrcs(2)%fld => frcyu_phys


    call esmf_lonlat2phys_regrid( lonlatfrcs, physfrcs )

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       call outfld('SFRCXU_phys',frcxu_phys(:,:,lchnk), pcols, lchnk)
       call outfld('SFRCYU_phys',frcyu_phys(:,:,lchnk), pcols, lchnk)
    end do

    call t_stopf ('gw_cospectra_calc')

  contains

    !==========================================================================
    function gather_and_smooth( flux ) result (sflx)
      real(r8),intent(in) :: flux(lat_beg:lat_end,1:pver)
      real(r8) :: sflx(lat_beg:lat_end,1:pver)

      real(r8) :: flxglb1(nlat,pver)
      real(r8) :: flxglb2(nlat,pver)

      integer :: k

      sflx = 0._r8

      flxglb1 = gather_fluxes(flux)

      do k = 1,pver
         flxglb2(1:nlat,k) = smooth(flxglb1(1:nlat,k),nlat,11)
      end do
      sflx(lat_beg:lat_end,:) = flxglb2(lat_beg:lat_end,:)

    end function gather_and_smooth

    !==========================================================================
    function smooth(data, n_points, window_size) result(smoothed_data)

      real(r8), intent(in) :: data(n_points)
      integer, intent(in) :: n_points, window_size

      real(r8) :: smoothed_data(n_points)

      integer :: ilat, j, start_idx, end_idx, count
      real(r8) :: sum_val

      smoothed_data = data ! Initialize with original data

      do ilat = 1, n_points
         sum_val = 0.0_r8
         count = 0
         start_idx = max(1,      ilat - ((window_size-1)/2))
         end_idx = min(n_points, ilat + ((window_size-1)/2))

         do j = start_idx, end_idx
            sum_val = sum_val + data(j)
            count = count + 1
         end do
         smoothed_data(ilat) = sum_val / real(count,kind=r8)
      end do

    end function smooth

    !==========================================================================
    function gather_fluxes( flx_loc ) result(flxglb)
      use mpi, only: MPI_REAL8, MPI_SUCCESS, MPI_SUM
      use esmf_lonlat_grid_mod, only: merid_comm

      real(r8),intent(in) :: flx_loc(lat_beg:lat_end,1:pver)

      real(r8) :: flxglb(nlat,pver)
      real(r8) :: sndbuf(nlat,pver)
      integer :: rc, len

      len = nlat*pver

      flxglb = 0._r8
      sndbuf = 0._r8
      sndbuf(lat_beg:lat_end,1:pver) = flx_loc(lat_beg:lat_end,1:pver)

      call mpi_allreduce(sndbuf,flxglb,len,MPI_REAL8,MPI_SUM, merid_comm, rc)
      if ( rc /= MPI_SUCCESS ) then
         call endrun('gw_cospectra_mod::gather_fluxes: mpi_allreduce FAILED')
      end if

    end function gather_fluxes

    !==========================================================================
    subroutine output_fld(out_fft, name)

      complex(C_DOUBLE_COMPLEX), intent(in) :: out_fft(nftnum, lat_beg:lat_end, pver)
      character(len=*), intent(in) :: name

      real(r8) ::  tmpr(nftnum, pver)
      real(r8) ::  tmpi(nftnum, pver)

      do icol = lat_beg, lat_end
         do n = 1,nftnum
            do k = 1,pver
               tmpr(n,k) = out_fft(n,icol,k)%re
               tmpi(n,k) = out_fft(n,icol,k)%im
            end do
         end do
         call outfld(trim(name)//'_FFT_real', tmpr, 1, icol)
         call outfld(trim(name)//'_FFT_imag', tmpi, 1, icol)
      end do

    end subroutine output_fld

    !==========================================================================
    subroutine output_cosp(out_fld, name)

      real(r8), intent(in) :: out_fld(nftnum, lat_beg:lat_end, pver)
      character(len=*), intent(in) :: name

      real(r8) ::  tmp(nftnum, pver)

      do icol = lat_beg, lat_end
         do n = 1,nftnum
            do k = 1,pver
               tmp(n,k) = out_fld(n,icol,k)
            end do
         end do
         call outfld(trim(name)//'WSTAR_cosp', tmp, 1, icol)
      end do

    end subroutine output_cosp

  end subroutine gw_cospectra_calc

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_restart_init(file)

    type(file_desc_t),  intent(inout) :: file

    integer :: ierr
    integer :: nftnum_dimid, nlat_dimid, ntime_dimid, nlev_dimid

    ierr = pio_def_dim(file, 'gw_csp_nftnum', nftnum, nftnum_dimid)
    ierr = pio_def_dim(file, 'gw_csp_nlat', nlat, nlat_dimid)
    ierr = pio_def_dim(file, 'gw_csp_ntime', ntime, ntime_dimid)
    ierr = pio_inq_dimid(file, 'lev', nlev_dimid)

    ierr = pio_def_var(file, 'gw_csp_accum_u', pio_double, &
         (/nftnum_dimid, nlat_dimid, nlev_dimid, ntime_dimid/),rest_accum_u_desc)
    ierr = pio_def_var(file, 'gw_csp_accum_v', pio_double, &
         (/nftnum_dimid, nlat_dimid, nlev_dimid, ntime_dimid/),rest_accum_v_desc)

  end subroutine gw_cospectra_restart_init

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_restart_write(file)
    use cam_pio_utils, only: pio_subsystem

    type(file_desc_t), intent(inout) :: file

    integer :: ierr
    type(io_desc_t) :: iodesc
    integer(PIO_OFFSET_KIND), pointer :: ldof(:)

    ldof => get_restart_decomp()
    call pio_initdecomp(pio_subsystem, pio_double, (/nftnum, nlat, pver, ntime/), ldof, iodesc)
    deallocate(ldof)

    call pio_write_darray(file, rest_accum_u_desc, iodesc, accum_cospectra_u, ierr)
    call pio_write_darray(file, rest_accum_v_desc, iodesc, accum_cospectra_v, ierr)

    call pio_freedecomp(file, iodesc)

  end subroutine gw_cospectra_restart_write

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_restart_read(file)
    use cam_pio_utils, only: pio_subsystem

    type(file_desc_t), intent(inout) :: file

    integer :: ierr
    type(io_desc_t) :: iodesc
    integer(PIO_OFFSET_KIND), pointer :: ldof(:)

    ierr = pio_inq_varid(file, 'gw_csp_accum_u', rest_accum_u_desc)
    ierr = pio_inq_varid(file, 'gw_csp_accum_v', rest_accum_v_desc)

    ldof => get_restart_decomp()
    call pio_initdecomp(pio_subsystem, pio_double, (/nftnum, nlat, pver, ntime/), ldof, iodesc)
    deallocate(ldof)

    call pio_read_darray(file, rest_accum_u_desc, iodesc, accum_cospectra_u, ierr)

    call pio_read_darray(file, rest_accum_v_desc, iodesc, accum_cospectra_v, ierr)

    call pio_freedecomp(file, iodesc)

  end subroutine gw_cospectra_restart_read

  ! utility routines
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_restart_decomp() result(ldof)

    integer(PIO_OFFSET_KIND), pointer :: ldof(:)

    ! local variables
    integer :: i, k, j, t
    integer :: lcnt

    lcnt = ntime*pver*(lat_end-lat_beg+1)*nftnum
    allocate(ldof(lcnt))
    ldof(:) = 0

    lcnt = 0

    do t = 1, ntime
       do k = 1,pver
          do j = lat_beg,lat_end
             do i = 1,nftnum

                lcnt = lcnt + 1
                ldof(lcnt) = i + (j-1)*nftnum + (k-1)*nftnum*nlat + (t-1)*nftnum*nlat*pver

             end do
          end do
       end do
    end do

  end function get_restart_decomp

end module gw_cospectra_mod
