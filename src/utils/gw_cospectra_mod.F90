module gw_cospectra_mod

  use shr_kind_mod, only: r8 => shr_kind_r8, cx => SHR_KIND_CX
  use ppgrid, only: pcols, pver, begchunk, endchunk
  use phys_grid, only: get_ncols_p
  use physics_types, only: physics_state
  use esmf_lonlat_grid_mod, only: nlon, nlat, glats, lon_beg, lon_end, lat_beg, lat_end
  use esmf_zonal_fft_mod, only : esmf_zonal_fft_3d, esmf_zonal_ifft_3d, esmf_zonal_fft_init
  use esmf_zonal_mean_mod, only: esmf_zonal_mean_calc, esmf_zonal_mean_wsums
  use esmf_phys2lonlat_mod, only: esmf_phys2lonlat_regrid
  use esmf_lonlat2phys_mod, only: esmf_lonlat2phys_init

  use, intrinsic :: iso_c_binding

  use spmd_utils, only: masterproc
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun

  use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field, pbuf_get_index

  use pio

  implicit none

  private
  public :: gw_cospectra_readnl
  public :: gw_cospectra_reg
  public :: gw_cospectra_init
  public :: gw_cospectra_calc
  public :: gw_cospectra_adj_tends
  public :: gw_cospectra_restart_init
  public :: gw_cospectra_restart_write
  public :: gw_cospectra_restart_read
  public :: gw_cospectra_final

  logical, protected, public :: gw_cospectra_active = .false.

  real(r8), pointer :: frcxu_phys(:,:,:) => null() !(pcols,pver,begchunk:endchunk)
  real(r8), pointer :: frcyu_phys(:,:,:) => null() !(pcols,pver,begchunk:endchunk)
  real(r8), pointer :: frcxyu_phys(:,:,:) => null() !(pcols,pver,begchunk:endchunk)

  integer :: nftnum = 0
  integer :: ntime = 0

  real(r8), allocatable :: accum_cospectra_u(:,:,:,:)
  real(r8), allocatable :: accum_cospectra_v(:,:,:,:)
  real(r8), allocatable :: accum_cospectra_uv(:,:,:,:)

  integer :: gw_cospectra_smooth_nlats = 0
  real(r8):: gw_cospectra_wavelen_begin = 0.0_r8 ! meters

contains

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils, only : mpicom, masterprocid, mpi_integer, mpi_success, mpi_real8

    character(len=*), intent(in) :: nlfile
    integer :: unitn, ierr
    character(len=cx) :: iomsg

    integer :: gw_cospectra_accum_ntimes

    character(len=*), parameter :: prefix = 'gw_cospectra_readnl: '

    namelist /gw_cospectra_nl/ gw_cospectra_accum_ntimes, gw_cospectra_smooth_nlats, &
         gw_cospectra_wavelen_begin

    if (masterproc) then
       ! read namelist
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'gw_cospectra_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, gw_cospectra_nl, iostat=ierr, iomsg=iomsg)
          if (ierr /= 0) then
             call endrun(prefix//'gw_cospectra_nl: ERROR reading namelist: '//trim(iomsg))
          end if
       else
          gw_cospectra_accum_ntimes = 0
       end if
       close(unitn)
    end if

    call mpi_bcast(gw_cospectra_accum_ntimes, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(prefix//'mpi_bcast error : gw_cospectra_accum_ntimes')
    call mpi_bcast(gw_cospectra_smooth_nlats, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(prefix//'mpi_bcast error : gw_cospectra_smooth_nlats')
    call mpi_bcast(gw_cospectra_wavelen_begin, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(prefix//'mpi_bcast error : gw_cospectra_wavelen_begin')

    ntime = gw_cospectra_accum_ntimes
    gw_cospectra_active = ntime > 0

    if (masterproc) then
       write(iulog,*) prefix//'gw_cospectra_active : ', gw_cospectra_active
       write(iulog,*) prefix//'gw_cospectra_accum_ntimes: ', ntime
       write(iulog,*) prefix//'gw_cospectra_smooth_nlats: ',gw_cospectra_smooth_nlats
       write(iulog,*) prefix//'gw_cospectra_wavelen_begin (meters): ',gw_cospectra_wavelen_begin
    end if

  end subroutine gw_cospectra_readnl

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_reg
    use cam_history_support, only: add_hist_coord

    if (.not.gw_cospectra_active) return

    nftnum = nlon/2+1

    ! add history coordinate for FT number
    call add_hist_coord('fft_num', nftnum, 'Fourier Transform Number')

    allocate(accum_cospectra_u( nftnum, lat_beg:lat_end, pver, ntime ))
    accum_cospectra_u = 0._r8
    allocate(accum_cospectra_v( nftnum, lat_beg:lat_end, pver, ntime ))
    accum_cospectra_v = 0._r8
    allocate(accum_cospectra_uv( nftnum, lat_beg:lat_end, pver, ntime ))
    accum_cospectra_uv = 0._r8

  end subroutine gw_cospectra_reg

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_init()
    use cam_history, only: addfld

    if (.not.gw_cospectra_active) return

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

    call addfld('FRCXR', (/'lev'/), 'I', '1', 'Resolved zonal forcing', gridname='ctem_zm')
    call addfld('FRCXU', (/'lev'/), 'I', '1', 'Unresolved zonal forcing', gridname='ctem_zm')
    call addfld('FRCYR', (/'lev'/), 'I', '1', 'Resolved meridianal forcing', gridname='ctem_zm')
    call addfld('FRCYU', (/'lev'/), 'I', '1', 'Unresolved meridianal forcing', gridname='ctem_zm')

    call addfld('SFRCXR', (/'lev'/), 'I', '1', 'Smoothed Resolved zonal forcing', gridname='ctem_zm')
    call addfld('SFRCXU', (/'lev'/), 'I', '1', 'Smoothed Unresolved zonal forcing', gridname='ctem_zm')
    call addfld('SFRCYR', (/'lev'/), 'I', '1', 'Smoothed Resolved meridianal forcing', gridname='ctem_zm')
    call addfld('SFRCYU', (/'lev'/), 'I', '1', 'Smoothed Unresolved meridianal forcing', gridname='ctem_zm')
    call addfld('SFRCXYR', (/'lev'/), 'I', '1', 'Smoothed Resolved zonal forcing by uv', gridname='ctem_zm')
    call addfld('SFRCXYU', (/'lev'/), 'I', '1', 'Smoothed Unresolved zonal forcing by uv', gridname='ctem_zm')

    call addfld('SFRCXU_phys', (/'lev'/), 'I', 'meters/sec2', 'Smoothed Unresolved zonal forcing', gridname='physgrid')
    call addfld('SFRCYU_phys', (/'lev'/), 'I', 'meters/sec2', 'Smoothed Unresolved meridianal zonal forcing', gridname='physgrid')
    call addfld('SFRCXYU_phys', (/'lev'/), 'I', 'meters/sec2', 'Smoothed Unresolved zonal forcing by uv', gridname='physgrid')

    call esmf_lonlat2phys_init()

    allocate(frcxu_phys(pcols,pver,begchunk:endchunk))
    frcxu_phys = 0._r8
    allocate(frcyu_phys(pcols,pver,begchunk:endchunk))
    frcyu_phys = 0._r8
    allocate(frcxyu_phys(pcols,pver,begchunk:endchunk))
    frcxyu_phys = 0._r8

  end subroutine gw_cospectra_init

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_calc(phys_state, pbuf2d)
    use cam_history, only: outfld
    use perf_mod, only: t_startf, t_stopf
    use cospext_mod, only: cospext, NOTSET, kxmin
    use air_composition, only: rairv  ! composition dependent gas constant (J/K/kg)
    use ref_pres, only: pref_mid, do_molec_diff, nbot_molec
    use esmf_phys2lonlat_mod, only: p2l_bdl=>fields_bundle_t, phys2lonlat_nflds=>nflds
    use esmf_lonlat2phys_mod, only: l2p_bdl=>fields_bundle_t, lonlat2phys_nflds=>nflds, esmf_lonlat2phys_regrid
    use physconst, only: cpair, gravit, rair
    use molec_diff, only: km_fac, km_pwr=>pwr  ! molecular viscosity constants (non-WACCM-X)

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    complex(C_DOUBLE_COMPLEX) :: u_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: v_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: w_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: t_fft(nftnum, lat_beg:lat_end, pver)

    complex(C_DOUBLE_COMPLEX) :: wstar(nftnum, lat_beg:lat_end, pver)

    real(r8),target :: ufld(pver,pcols,begchunk:endchunk)
    real(r8),target :: vfld(pver,pcols,begchunk:endchunk)
    real(r8),target :: wfld(pver,pcols,begchunk:endchunk)
    real(r8),target :: tfld(pver,pcols,begchunk:endchunk)
    real(r8),target :: kvmfld(pver,pcols,begchunk:endchunk) ! kvm at midpoints
    real(r8),target :: nmfld(pver,pcols,begchunk:endchunk)  ! Brunt-Vaisala frequency at midpoints
    integer :: lchnk, ncol, icol

    real(r8),target :: rho_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: u_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: v_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: w_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: t_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: kvm_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: nm_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)

    real(r8) :: kvm_zm(lat_beg:lat_end,pver)
    real(r8) :: nm_zm(lat_beg:lat_end,pver)
    real(r8) :: wvlx_visc(lat_beg:lat_end,pver)
    real(r8) :: wvlxend(lat_beg:lat_end,pver)
    real(r8), pointer :: kvm_pbuf(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)
    integer, save :: kvm_pbuf_idx = -1
    logical, save :: kvm_idx_set = .false.
    real(r8), parameter :: wvlx_visc_floor = 4.e4_r8 ! lower bound on wvlxend (m)
    real(r8), parameter :: n2min = 5.e-5_r8          ! floor on N^2 (matches gw prepare_profiles)
    real(r8) :: ti_col(pver+1), rhoi_col(pver+1), ni_col(pver+1)
    real(r8) :: mkvisc_i(pver+1)
    real(r8) :: dtdp, n2

    ! band-pass filtered spectra (wavenumbers outside [kxbeg,kxend] zeroed)
    complex(C_DOUBLE_COMPLEX) :: u_fft_f(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: v_fft_f(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: w_fft_f(nftnum, lat_beg:lat_end, pver)

    ! filtered fields in physical (lon,lat) space
    real(r8) :: u_filt(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8) :: v_filt(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8) :: w_filt(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8) :: t_filt(lon_beg:lon_end,lat_beg:lat_end,pver)

    ! element-wise products of filtered fields
    real(r8) :: uw_filt(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8) :: vw_filt(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8) :: uv_filt(lon_beg:lon_end,lat_beg:lat_end,pver)

    ! wavenumber-band edges and scaling factors (positive/negative branches)
    integer :: kxbeg(lat_beg:lat_end), kxl(lat_beg:lat_end), kx
    integer :: kxend(lat_beg:lat_end,pver)
    integer :: lat0, lat1
    real(r8) :: fp_uw(lat_beg:lat_end,pver), fn_uw(lat_beg:lat_end,pver)
    real(r8) :: fp_vw(lat_beg:lat_end,pver), fn_vw(lat_beg:lat_end,pver)
    real(r8) :: fp_uv(lat_beg:lat_end,pver), fn_uv(lat_beg:lat_end,pver)
    ! zonal sums of absolute filtered products, used as weighting norms
    real(r8) :: uw_abs_sum(lat_beg:lat_end,pver)
    real(r8) :: vw_abs_sum(lat_beg:lat_end,pver)
    real(r8) :: uv_abs_sum(lat_beg:lat_end,pver)
    real(r8), parameter :: sum_tiny = 1.e-30_r8
    real(r8), parameter :: rearth = 6.371e6_r8
    real(r8), parameter :: krat = 1.6_r8
    real(r8), parameter :: lat_scl_max = 85._r8   ! only scale equatorward of this |lat| (deg)

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

    real(r8) :: mflxxrp(lat_beg:lat_end,pver), mflxxrn(lat_beg:lat_end,pver)
    real(r8) :: mflxyrp(lat_beg:lat_end,pver), mflxyrn(lat_beg:lat_end,pver)
    real(r8) :: mflxxyrp(lat_beg:lat_end,pver), mflxxyrn(lat_beg:lat_end,pver)
    real(r8) :: mflxxup(lat_beg:lat_end,pver), mflxxun(lat_beg:lat_end,pver)
    real(r8) :: mflxyup(lat_beg:lat_end,pver), mflxyun(lat_beg:lat_end,pver)
    real(r8) :: mflxxyup(lat_beg:lat_end,pver), mflxxyun(lat_beg:lat_end,pver)
    real(r8) :: frcxr(lat_beg:lat_end,pver), frcxu(lat_beg:lat_end,pver)
    real(r8) :: frcyr(lat_beg:lat_end,pver), frcyu(lat_beg:lat_end,pver)
    real(r8) :: frcxyr(lat_beg:lat_end,pver), frcxyu(lat_beg:lat_end,pver)

    ! spectral slopes for positive and negative branches (zonal, meridional, and uv components)
    real(r8) :: slppx(lat_beg:lat_end,pver),  slpnx(lat_beg:lat_end,pver)
    real(r8) :: slppy(lat_beg:lat_end,pver),  slpny(lat_beg:lat_end,pver)
    real(r8) :: slppxy(lat_beg:lat_end,pver), slpnxy(lat_beg:lat_end,pver)

    real(r8),target :: frcxu_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: frcyu_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8),target :: frcxyu_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)

    real(r8) :: wvlxbeg
    real(r8) :: mflux_glb(nlat,pver)
    integer :: kk

    real(r8), parameter :: pi = 4._r8*atan(1._r8)
    real(r8), parameter :: deg2rad = pi/180._r8
    character(len=*), parameter :: subname  = 'gw_cospectra_calc'

    if (.not.gw_cospectra_active) return

    wvlxbeg = gw_cospectra_wavelen_begin

    call t_startf ('gw_cospectra_calc')

    if (.not. kvm_idx_set) then
       kvm_pbuf_idx = pbuf_get_index('kvm')
       kvm_idx_set = .true.
    end if

    latrad(lat_beg:lat_end) = glats(lat_beg:lat_end)*deg2rad

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)

       ! Pull previous-step kvm (interfaces, pverp) from pbuf, average to midpoints
       pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
       call pbuf_get_field(pbuf_chnk, kvm_pbuf_idx, kvm_pbuf)

       do icol = 1,ncol
          ufld(:pver,icol,lchnk) = phys_state(lchnk)%u(icol,:pver)
          vfld(:pver,icol,lchnk) = phys_state(lchnk)%v(icol,:pver)
          wfld(:pver,icol,lchnk) = phys_state(lchnk)%omega(icol,:pver)
          tfld(:pver,icol,lchnk) = phys_state(lchnk)%t(icol,:pver) * phys_state(lchnk)%exner(icol,:pver)
          rho (:pver,icol,lchnk) = phys_state(lchnk)%pmid(icol,:pver)/(rairv(icol,:pver,lchnk)*phys_state(lchnk)%t(icol,:pver)) ! kg/m3

          ! Brunt-Vaisala frequency on midpoints (same formulation as
          ! gravity_wave_drag_prepare_profiles_run in atmos_phys). Build interface
          ! T, density, and N^2 from phys_state, then average to midpoints.
          ti_col(1) = phys_state(lchnk)%t(icol,1)
          rhoi_col(1) = phys_state(lchnk)%pint(icol,1)/(rair*ti_col(1))
          ni_col(1) = sqrt(gravit*gravit/(cpair*ti_col(1)))
          do kk = 2,pver
             ti_col(kk) = 0.5_r8*(phys_state(lchnk)%t(icol,kk-1) + phys_state(lchnk)%t(icol,kk))
             rhoi_col(kk) = phys_state(lchnk)%pint(icol,kk)/(rair*ti_col(kk))
             dtdp = (phys_state(lchnk)%t(icol,kk) - phys_state(lchnk)%t(icol,kk-1)) &
                  / (phys_state(lchnk)%pmid(icol,kk) - phys_state(lchnk)%pmid(icol,kk-1))
             n2 = gravit*gravit/ti_col(kk)*(1._r8/cpair - rhoi_col(kk)*dtdp)
             ni_col(kk) = sqrt(max(n2min, n2))
          end do
          ti_col(pver+1) = phys_state(lchnk)%t(icol,pver)
          rhoi_col(pver+1) = phys_state(lchnk)%pint(icol,pver+1)/(rair*ti_col(pver+1))
          ni_col(pver+1) = ni_col(pver)
          do kk = 1,pver
             nmfld(kk,icol,lchnk) = 0.5_r8*(ni_col(kk) + ni_col(kk+1))
          end do

          ! kvm total at interfaces = eddy kvm (from pbuf, previous step) +
          ! molecular kinematic viscosity (non-WACCM-X branch of molec_diff.F90:
          ! mkvisc = km_fac * tint^pwr / rhoi, applied for k = 1..nbot_molec).
          ! Then average to midpoints.
          mkvisc_i(:) = 0._r8
          if (do_molec_diff) then
             do kk = 1,nbot_molec
                mkvisc_i(kk) = km_fac * ti_col(kk)**km_pwr / rhoi_col(kk)
             end do
          end if
          do kk = 1,pver
             kvmfld(kk,icol,lchnk) = 0.5_r8*((kvm_pbuf(icol,kk)   + mkvisc_i(kk)) &
                                           + (kvm_pbuf(icol,kk+1) + mkvisc_i(kk+1)))
          end do
       end do
    end do

    physflds(1)%fld => ufld
    physflds(2)%fld => vfld
    physflds(3)%fld => wfld
    physflds(4)%fld => tfld
    physflds(5)%fld => rho
    physflds(6)%fld => kvmfld
    physflds(7)%fld => nmfld

    lonlatflds(1)%fld => u_lonlat
    lonlatflds(2)%fld => v_lonlat
    lonlatflds(3)%fld => w_lonlat
    lonlatflds(4)%fld => t_lonlat
    lonlatflds(5)%fld => rho_lonlat
    lonlatflds(6)%fld => kvm_lonlat
    lonlatflds(7)%fld => nm_lonlat

    call esmf_phys2lonlat_regrid(physflds, lonlatflds)

    ! Zonal-mean kvm and nm, then viscous horizontal scale and wvlxend
    call esmf_zonal_mean_calc(kvm_lonlat, kvm_zm)
    call esmf_zonal_mean_calc(nm_lonlat, nm_zm)
    do k = 1, pver
       do j = lat_beg, lat_end
          if (nm_zm(j,k) > 0._r8 .and. kvm_zm(j,k) > 0._r8) then
             wvlx_visc(j,k) = 2._r8*pi*sqrt(kvm_zm(j,k)/nm_zm(j,k))
          else
             wvlx_visc(j,k) = 0._r8
          end if
          if (wvlx_visc(j,k) < wvlx_visc_floor) then
             wvlxend(j,k) = wvlx_visc_floor
          else if (wvlx_visc(j,k) <= wvlxbeg) then
             wvlxend(j,k) = wvlx_visc(j,k)
          else
             wvlxend(j,k) = wvlxbeg
          end if
       end do
    end do

    u_fft = esmf_zonal_fft_3d(u_lonlat)
    call output_fld(u_fft, name='U')

    v_fft = esmf_zonal_fft_3d(v_lonlat)
    call output_fld(v_fft, name='V')

    w_fft = esmf_zonal_fft_3d(w_lonlat)
    call output_fld(w_fft, name='OMEGA')

    t_fft = esmf_zonal_fft_3d(t_lonlat)
    call output_fld(t_fft, name='THETA')

    ! zonal band-pass filter: keep only wavenumbers in [kxbeg,kxend] per latitude
    ! and level, then inverse FFT back to physical space.
    ! kxl matches cospext's scale-invariant lower bound (kxbeg/krat/krat).
    do j = lat_beg, lat_end
       kxbeg(j) = nint(2._r8*pi*rearth*cos(latrad(j))/wvlxbeg)
       kxl(j)   = nint(nint(kxbeg(j)/krat)/krat)
    end do
    do k = 1, pver
       do j = lat_beg, lat_end
          kxend(j,k) = nint(2._r8*pi*rearth*cos(latrad(j))/wvlxend(j,k))
       end do
    end do

    u_fft_f = u_fft
    v_fft_f = v_fft
    w_fft_f = w_fft
    do k = 1, pver
       do j = lat_beg, lat_end
          do i = 1, nftnum
             kx = i - 1
             if (kx < kxbeg(j) .or. kx > kxend(j,k)) then
                u_fft_f(i,j,k) = (0._r8, 0._r8)
                v_fft_f(i,j,k) = (0._r8, 0._r8)
                w_fft_f(i,j,k) = (0._r8, 0._r8)
             end if
          end do
       end do
    end do

    u_filt = esmf_zonal_ifft_3d(u_fft_f)
    v_filt = esmf_zonal_ifft_3d(v_fft_f)
    w_filt = esmf_zonal_ifft_3d(w_fft_f)

    ! element-wise products of the filtered fields
    uw_filt = u_filt * w_filt
    vw_filt = v_filt * w_filt
    uv_filt = u_filt * v_filt

    ! Longitudinal smoothing (moving average) of the filtered products.
    ! Window ~= wvlxbeg in physical distance, with wrap-around at 0/360 deg.
    call lon_smooth_3d(uw_filt, wvlxbeg)
    call lon_smooth_3d(vw_filt, wvlxbeg)
    call lon_smooth_3d(uv_filt, wvlxbeg)

    wstar = conjg(w_fft)
    call output_fld(wstar, name='WSTAR')

    do i = 1,ntime-1
       accum_cospectra_u(:,:,:,i) = accum_cospectra_u(:,:,:,i+1)
       accum_cospectra_v(:,:,:,i) = accum_cospectra_v(:,:,:,i+1)
       accum_cospectra_uv(:,:,:,i) = accum_cospectra_uv(:,:,:,i+1)
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

    tmpfld = u_fft * conjg(v_fft)
    cospectra = tmpfld%re
    cospectra(2:,:,:) = 2._r8 * cospectra(2:,:,:)

    accum_cospectra_uv(:,:,:,ntime) = cospectra(:,:,:)

    tmpfld = t_fft * wstar
    cospectra = tmpfld%re
    cospectra(2:,:,:) = 2._r8 * cospectra(2:,:,:)
    call output_cosp(cospectra,'T')

    call esmf_zonal_mean_calc(rho_lonlat, rhobar)

    do icol = lat_beg, lat_end
       call outfld('RHOBAR', rhobar(icol,:),1,icol)
    end do

    ! zonal component
    call cospext(nftnum, lat_beg,lat_end, pver,ntime, latrad, accum_cospectra_u, wvlxbeg,wvlxend, pref_mid, rhobar, mflxxup,mflxxun,mflxxrp,mflxxrn,force_r=frcxr,force_u=frcxu)

    ! meridional component
    call cospext(nftnum, lat_beg,lat_end, pver,ntime, latrad, accum_cospectra_v, wvlxbeg,wvlxend, pref_mid, rhobar, mflxyup,mflxyun,mflxyrp,mflxyrn,force_r=frcyr,force_u=frcyu)

    ! meridional momentum flux
    call cospext(nftnum, lat_beg,lat_end, pver,ntime, latrad, accum_cospectra_uv, wvlxbeg,wvlxend, pref_mid, rhobar, mflxxyup,mflxxyun,mflxxyrp,mflxxyrn)

    ! Scale filtered products by branch-dependent factors derived from the spectral
    ! slopes (see cospext_mod::momentum_fluxes for the fp/fn formulation).
    ! Scaling is applied only equatorward of lat_scl_max, matching cospext.
    lat0 = 0
    lat1 = -1
    do j = lat_beg, lat_end
       if ( abs(latrad(j))*180._r8/pi < lat_scl_max ) then
          if (lat0 < 1) lat0 = j
          lat1 = j
       end if
    end do

!    call compute_scaling(slppx,  slpnx,  fp_uw, fn_uw)
!    call compute_scaling(slppy,  slpny,  fp_vw, fn_vw)
!    call compute_scaling(slppxy, slpnxy, fp_uv, fn_uv)

    frcxyu = gather_and_latderiv(mflxxyup,mflxxyun,rhobar)
    frcxyr = gather_and_latderiv(mflxxyrp,mflxxyrn,rhobar)

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

    ! gather and smooth in latitude

    frcxr = gather_and_smooth(frcxr)
    frcxu = gather_and_smooth(frcxu)
    frcyr = gather_and_smooth(frcyr)
    frcyu = gather_and_smooth(frcyu)
    frcxyr = gather_and_smooth(frcxyr)
    frcxyu = gather_and_smooth(frcxyu)

    ! Where the viscous horizontal scale exceeds wvlxbeg the unresolved window
    ! collapses, so the unresolved forcings are set to zero.
    do k = 1, pver
       do j = lat_beg, lat_end
          if (wvlx_visc(j,k) > wvlxbeg) then
             frcxu(j,k)  = 0._r8
             frcyu(j,k)  = 0._r8
             frcxyu(j,k) = 0._r8
          end if
       end do
    end do

    do icol = lat_beg, lat_end

       call outfld('SFRCXR', frcxr(icol,:),1,icol)
       call outfld('SFRCXU', frcxu(icol,:),1,icol)
       call outfld('SFRCYR', frcyr(icol,:),1,icol)
       call outfld('SFRCYU', frcyu(icol,:),1,icol)
       call outfld('SFRCXYR', frcxyr(icol,:),1,icol)
       call outfld('SFRCXYU', frcxyu(icol,:),1,icol)

    end do

    ! Distribute the zonal-mean forcings across longitudes using the amplitude
    ! (absolute value) of the filtered momentum-flux products as the weight:
    !     frcxu_lonlat(i,j,k) = frcxu(j,k) * nlon * |uw_filt(i,j,k)| / sum_i(|uw_filt|)
    ! and analogously for frcyu_lonlat (|vw_filt|) and frcxyu_lonlat (|uv_filt|).
    ! By construction the zonal mean of frcxu_lonlat equals frcxu, so the
    ! zonal-mean forcing is preserved while the longitudinal structure follows
    ! the gravity-wave amplitude. Only applied equatorward of lat_scl_max;
    ! latitudes outside stay at zero to add no gravity-wave forcing.
    uw_abs_sum = esmf_zonal_mean_wsums(abs(uw_filt))
    vw_abs_sum = esmf_zonal_mean_wsums(abs(vw_filt))
    uv_abs_sum = esmf_zonal_mean_wsums(abs(uv_filt))

    frcxu_lonlat  = 0._r8
    frcyu_lonlat  = 0._r8
    frcxyu_lonlat = 0._r8

    do k = 1, pver
       do j = lat0, lat1
          if (wvlx_visc(j,k) > wvlxbeg) cycle
          if (uw_abs_sum(j,k) > sum_tiny) then
             do i = lon_beg, lon_end
                frcxu_lonlat(i,j,k) = frcxu(j,k) * real(nlon,r8) * abs(uw_filt(i,j,k)) / uw_abs_sum(j,k)
             end do
          end if
          if (vw_abs_sum(j,k) > sum_tiny) then
             do i = lon_beg, lon_end
                frcyu_lonlat(i,j,k) = frcyu(j,k) * real(nlon,r8) * abs(vw_filt(i,j,k)) / vw_abs_sum(j,k)
             end do
          end if
          if (uv_abs_sum(j,k) > sum_tiny) then
             do i = lon_beg, lon_end
                frcxyu_lonlat(i,j,k) = frcxyu(j,k) * real(nlon,r8) * abs(uv_filt(i,j,k)) / uv_abs_sum(j,k)
             end do
          end if
       end do
    end do

    frcxu_phys = -huge(1._r8)
    frcyu_phys = -huge(1._r8)
    frcxyu_phys = -huge(1._r8)

    lonlatfrcs(1)%fld => frcxu_lonlat
    lonlatfrcs(2)%fld => frcyu_lonlat
    lonlatfrcs(3)%fld => frcxyu_lonlat
    physfrcs(1)%fld => frcxu_phys
    physfrcs(2)%fld => frcyu_phys
    physfrcs(3)%fld => frcxyu_phys

    ! map forcings to physics grid
    call esmf_lonlat2phys_regrid( lonlatfrcs, physfrcs )

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       call outfld('SFRCXU_phys',frcxu_phys(:,:,lchnk), pcols, lchnk)
       call outfld('SFRCYU_phys',frcyu_phys(:,:,lchnk), pcols, lchnk)
       call outfld('SFRCXYU_phys',frcxyu_phys(:,:,lchnk), pcols, lchnk)
    end do

    call t_stopf ('gw_cospectra_calc')

  contains

    !==========================================================================
    ! Vertical derivative dy/dz on the model level grid (one-sided at boundaries,
    ! centered in the interior). Same as cospext_mod::vertdiv.
    function vertdiv( z, y ) result(dydz)
      real(r8), intent(in) :: z(pver)
      real(r8), intent(in) :: y(pver)
      real(r8) :: dydz(pver)

      integer :: kk

      do kk = 2, pver-1
         dydz(kk) = (y(kk+1) - y(kk-1)) / (z(kk+1) - z(kk-1))
      end do
      dydz(1)    = (y(2)    - y(1))      / (z(2)    - z(1))
      dydz(pver) = (y(pver) - y(pver-1)) / (z(pver) - z(pver-1))

    end function vertdiv

    !==========================================================================
    ! Compute per-latitude/level positive and negative branch scaling factors
    ! from spectral slopes, using the same formulation as cospext_mod's
    ! momentum_fluxes routine. fp/fn are set to 0 where the corresponding
    ! slope equals NOTSET or where kxbeg is at/below kxmin.
    subroutine compute_scaling(slpp, slpn, fp, fn)
      real(r8), intent(in)  :: slpp(lat_beg:lat_end,pver), slpn(lat_beg:lat_end,pver)
      real(r8), intent(out) :: fp(lat_beg:lat_end,pver),   fn(lat_beg:lat_end,pver)

      real(r8) :: bp, bn
      integer  :: jj, kk

      fp = 0._r8
      fn = 0._r8

      do kk = 1, pver
         do jj = lat0, lat1
            if (kxbeg(jj) > kxmin) then

               if (slpp(jj,kk) /= NOTSET) then
                  if (slpp(jj,kk) == 1._r8) then
                     fp(jj,kk) = log(real(kxend(jj,kk),r8)/real(kxbeg(jj),r8)) &
                               / log(real(kxbeg(jj),r8)/real(kxl(jj),r8))
                  else
                     bp = 1._r8 - slpp(jj,kk)
                     fp(jj,kk) = (real(kxend(jj,kk),r8)**bp - real(kxbeg(jj),r8)**bp) &
                               / (real(kxbeg(jj),r8)**bp - real(kxl(jj),r8)**bp)
                  end if
               end if

               if (slpn(jj,kk) /= NOTSET) then
                  if (slpn(jj,kk) == 1._r8) then
                     fn(jj,kk) = log(real(kxend(jj,kk),r8)/real(kxbeg(jj),r8)) &
                               / log(real(kxbeg(jj),r8)/real(kxl(jj),r8))
                  else
                     bn = 1._r8 - slpn(jj,kk)
                     fn(jj,kk) = (real(kxend(jj,kk),r8)**bn - real(kxbeg(jj),r8)**bn) &
                               / (real(kxbeg(jj),r8)**bn - real(kxl(jj),r8)**bn)
                  end if
               end if

            end if
         end do
      end do

    end subroutine compute_scaling

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
         flxglb2(1:nlat,k) = smooth(flxglb1(1:nlat,k),nlat,gw_cospectra_smooth_nlats)
      end do
      sflx(lat_beg:lat_end,:) = flxglb2(lat_beg:lat_end,:)

    end function gather_and_smooth

    !==========================================================================
    function gather_and_latderiv( fluxp, fluxn, rho) result (frc)
      real(r8),intent(in) :: fluxp(lat_beg:lat_end,1:pver),fluxn(lat_beg:lat_end,1:pver),rho(lat_beg:lat_end,1:pver)
      real(r8) :: frc(lat_beg:lat_end,1:pver)

      real(r8) :: flxglb1(nlat,pver)
      real(r8) :: flxglb2(nlat,pver)
      real(r8) :: rhoglb(nlat,pver)
      real(r8) :: frcglb(nlat,pver)
      real(r8), parameter :: pi = 4._r8*atan(1._r8)
      real(r8), parameter :: re = 6.371e6_r8

      real(r8) :: dy

      integer :: k,j

      frc = 0._r8

      flxglb1 = gather_fluxes(fluxp)
      flxglb2 = gather_fluxes(fluxn)
      rhoglb = gather_fluxes(rho)

      dy = pi*re/(nlat-1)

      do k = 1,pver
         do j = 2,nlat-1
            frcglb(j,k) = -((flxglb1(j+1,k)+flxglb2(j+1,k))*rhoglb(j+1,k)-(flxglb1(j-1,k)+flxglb2(j-1,k))*rhoglb(j-1,k))/2._r8/dy/rhoglb(j,k)
         end do
         frcglb(1,k) = 0._r8
         frcglb(nlat,k) = 0._r8
      end do
      frc(lat_beg:lat_end,:) = frcglb(lat_beg:lat_end,:)

    end function gather_and_latderiv

    !==========================================================================
    function smooth(data, n_points, window_size) result(smoothed_data)

      integer, intent(in) :: n_points, window_size
      real(r8), intent(in) :: data(n_points)

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
    ! Longitudinal moving-average smoothing of a 3D lon-lat-lev field, with
    ! wrap-around at 0/360 deg. The window in physical distance is wvlx_win
    ! (meters); the number of grid points varies with latitude via
    ! nlon / kxbeg(j) = wvlx_win * nlon / (2*pi*rearth*cos(lat)) and is capped
    ! to avoid overlapping the wrap.
    subroutine lon_smooth_3d(fld, wvlx_win)
      use mpi, only: MPI_REAL8, MPI_SUCCESS, MPI_SUM
      use esmf_lonlat_grid_mod, only: zonal_comm

      real(r8), intent(inout) :: fld(lon_beg:lon_end, lat_beg:lat_end, pver)
      real(r8), intent(in)    :: wvlx_win

      real(r8) :: sndbf(nlon, lat_beg:lat_end, pver)
      real(r8) :: rcvbf(nlon, lat_beg:lat_end, pver)
      integer  :: rc, len
      integer  :: ii, jj, kk, mm, idx, half, win
      real(r8) :: acc, coslat

      len = nlon * (lat_end-lat_beg+1) * pver

      sndbf = 0._r8
      sndbf(lon_beg:lon_end, lat_beg:lat_end, 1:pver) = &
           fld(lon_beg:lon_end, lat_beg:lat_end, 1:pver)

      call mpi_allreduce(sndbf, rcvbf, len, MPI_REAL8, MPI_SUM, zonal_comm, rc)
      if (rc /= MPI_SUCCESS) then
         call endrun('gw_cospectra_mod::lon_smooth_3d: mpi_allreduce FAILED')
      end if

      do jj = lat_beg, lat_end
         coslat = max(cos(latrad(jj)), 1.e-6_r8)
         win  = max(1, nint(wvlx_win * real(nlon,r8) / (2._r8*pi*rearth*coslat)))
         half = min(win/2, (nlon-1)/2)
         do kk = 1, pver
            do ii = lon_beg, lon_end
               acc = 0._r8
               do mm = -half, half
                  idx = ii + mm
                  if (idx < 1)    idx = idx + nlon
                  if (idx > nlon) idx = idx - nlon
                  acc = acc + rcvbf(idx, jj, kk)
               end do
               fld(ii, jj, kk) = acc / real(2*half+1, r8)
            end do
         end do
      end do

    end subroutine lon_smooth_3d

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
  subroutine gw_cospectra_adj_tends( ncol, lchnk, utend, vtend )
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(inout) :: utend(:,:)
    real(r8), intent(inout) :: vtend(:,:)

    integer :: i,k
    real(r8) :: utgw_cosp   ! temporary array for deduced unresolved forcing
    real(r8) :: vtgw_cosp
    real(r8), parameter :: gwtnd_cosp_max = 400._r8/86400._r8

    ! apply vertical smoothing
    do k = 2, pver-1
       do i = 1, ncol
          utgw_cosp = sum(frcxu_phys(i,k-1:k+1,lchnk)+frcxyu_phys(i,k-1:k+1,lchnk))/3._r8
          utgw_cosp = utend(i,k) + utgw_cosp
          utend(i,k) = SIGN(MIN(ABS(utgw_cosp),gwtnd_cosp_max),utgw_cosp)

          vtgw_cosp = sum(frcyu_phys(i,k-1:k+1,lchnk))/3._r8
          vtgw_cosp = vtend(i,k) + vtgw_cosp
          vtend(i,k) = SIGN(MIN(ABS(vtgw_cosp),gwtnd_cosp_max),vtgw_cosp)
       end do
    end do
    do i = 1, ncol
       utgw_cosp = frcxu_phys(i,1,lchnk)+frcxyu_phys(i,1,lchnk)
       utgw_cosp = utend(i,1) + utgw_cosp
       utend(i,1) = SIGN(MIN(ABS(utgw_cosp),gwtnd_cosp_max),utgw_cosp)

       vtgw_cosp = frcyu_phys(i,1,lchnk)
       vtgw_cosp = vtend(i,1) + vtgw_cosp
       vtend(i,1) = SIGN(MIN(ABS(vtgw_cosp),gwtnd_cosp_max),vtgw_cosp)
    end do

    do i = 1, ncol
       utgw_cosp = frcxu_phys(i,pver,lchnk)+frcxyu_phys(i,pver,lchnk)
       utgw_cosp = utend(i,pver) + utgw_cosp
       utend(i,pver) = SIGN(MIN(ABS(utgw_cosp),gwtnd_cosp_max),utgw_cosp)

       vtgw_cosp = frcyu_phys(i,pver,lchnk)
       vtgw_cosp = vtend(i,pver) + vtgw_cosp
       vtend(i,pver) = SIGN(MIN(ABS(vtgw_cosp),gwtnd_cosp_max),vtgw_cosp)
    end do

  end subroutine gw_cospectra_adj_tends

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_final()

    if (.not.gw_cospectra_active) return

    ! free up memory
    deallocate(accum_cospectra_u)
    deallocate(accum_cospectra_v)
    deallocate(accum_cospectra_uv)
    deallocate(frcxu_phys)
    nullify(frcxu_phys)
    deallocate(frcyu_phys)
    nullify(frcyu_phys)
    deallocate(frcxyu_phys)
    nullify(frcxyu_phys)
  end subroutine gw_cospectra_final

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_restart_init(file)

    type(file_desc_t),  intent(inout) :: file

    integer :: ierr, t
    integer :: nftnum_dimid, nlat_dimid, nlev_dimid
    type(var_desc_t) :: rest_accum_u_desc, rest_accum_v_desc, rest_accum_uv_desc
    character(len=2) :: numstr

    if (.not.gw_cospectra_active) return

    ierr = pio_def_dim(file, 'gw_csp_nftnum', nftnum, nftnum_dimid)
    ierr = pio_def_dim(file, 'gw_csp_nlat', nlat, nlat_dimid)
    ierr = pio_inq_dimid(file, 'lev', nlev_dimid)

    do t = 1,ntime
       write(numstr,'(I2.2)') t
       ierr = pio_def_var(file, 'gw_csp_accum_u_t'//numstr, pio_double, &
            (/nftnum_dimid, nlat_dimid, nlev_dimid/),rest_accum_u_desc)
       ierr = pio_def_var(file, 'gw_csp_accum_v_t'//numstr, pio_double, &
            (/nftnum_dimid, nlat_dimid, nlev_dimid/),rest_accum_v_desc)
       ierr = pio_def_var(file, 'gw_csp_accum_uv_t'//numstr, pio_double, &
            (/nftnum_dimid, nlat_dimid, nlev_dimid/),rest_accum_uv_desc)
    end do

  end subroutine gw_cospectra_restart_init

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_restart_write(file)
    use cam_pio_utils, only: pio_subsystem

    type(file_desc_t), intent(inout) :: file

    integer :: ierr, t
    type(io_desc_t) :: iodesc
    integer(PIO_OFFSET_KIND), pointer :: ldof(:)
    character(len=2) :: numstr
    type(var_desc_t) :: u_desc, v_desc, uv_desc

    if (.not.gw_cospectra_active) return

    ldof => get_restart_decomp(active=lon_beg==1)
    call pio_initdecomp(pio_subsystem, pio_double, (/nftnum, nlat, pver/), ldof, iodesc)
    deallocate(ldof)

    do t = 1,ntime
       write(numstr,'(I2.2)') t
       ierr = pio_inq_varid(file, 'gw_csp_accum_u_t'//numstr, u_desc)
       ierr = pio_inq_varid(file, 'gw_csp_accum_v_t'//numstr, v_desc)
       ierr = pio_inq_varid(file, 'gw_csp_accum_uv_t'//numstr, uv_desc)

       call pio_write_darray(file, u_desc, iodesc, accum_cospectra_u(:,:,:,t), ierr)
       call pio_write_darray(file, v_desc, iodesc, accum_cospectra_v(:,:,:,t), ierr)
       call pio_write_darray(file, uv_desc, iodesc, accum_cospectra_uv(:,:,:,t), ierr)
    end do

    call pio_freedecomp(file, iodesc)

  end subroutine gw_cospectra_restart_write

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine gw_cospectra_restart_read(file)
    use cam_pio_utils, only: pio_subsystem

    type(file_desc_t), intent(inout) :: file

    integer :: ierr, t
    type(io_desc_t) :: iodesc
    integer(PIO_OFFSET_KIND), pointer :: ldof(:)
    character(len=2) :: numstr
    type(var_desc_t) :: u_desc, v_desc, uv_desc

    if (.not.gw_cospectra_active) return

    ldof => get_restart_decomp(active=.true.)
    call pio_initdecomp(pio_subsystem, pio_double, (/nftnum, nlat, pver/), ldof, iodesc)
    deallocate(ldof)

    do t = 1,ntime
       write(numstr,'(I2.2)') t
       ierr = pio_inq_varid(file, 'gw_csp_accum_u_t'//numstr, u_desc)
       ierr = pio_inq_varid(file, 'gw_csp_accum_v_t'//numstr, v_desc)
       ierr = pio_inq_varid(file, 'gw_csp_accum_uv_t'//numstr, uv_desc)

       call pio_read_darray(file, u_desc, iodesc, accum_cospectra_u(:,:,:,t), ierr)
       call pio_read_darray(file, v_desc, iodesc, accum_cospectra_v(:,:,:,t), ierr)
       call pio_read_darray(file, uv_desc, iodesc, accum_cospectra_uv(:,:,:,t), ierr)
    end do

    call pio_freedecomp(file, iodesc)

  end subroutine gw_cospectra_restart_read

  ! utility routines
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function get_restart_decomp(active) result(ldof)

    logical, intent(in) :: active

    integer(PIO_OFFSET_KIND), pointer :: ldof(:)

    ! local variables
    integer :: i, k, j
    integer :: lcnt

    lcnt = pver*(lat_end-lat_beg+1)*nftnum
    allocate(ldof(lcnt))
    ldof(:) = 0

    if (active) then
       lcnt = 0
       do k = 1,pver
          do j = lat_beg,lat_end
             do i = 1,nftnum
                lcnt = lcnt + 1
                ldof(lcnt) = i + (j-1)*nftnum + (k-1)*nftnum*nlat
             end do
          end do
       end do
    endif

  end function get_restart_decomp

end module gw_cospectra_mod
