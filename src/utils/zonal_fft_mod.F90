module zonal_fft_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver, begchunk, endchunk
  use phys_grid, only: get_ncols_p
  use physics_types, only: physics_state
  use esmf_zonal_ops, only : lat_beg,lat_end, lon_beg,lon_end, nlons, glats, nlats, zonal_comm
  use esmf_zonal_ops, only : esmf_zonal_fft_3d, esmf_zonal_mean_3d
  use, intrinsic :: iso_c_binding

  use spmd_utils, only: masterproc
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun

  implicit none

  integer :: nftnum = 0
  integer :: ntime = 0

  real(r8), allocatable :: accum_cospectra_u(:,:,:,:)
  real(r8), allocatable :: accum_cospectra_v(:,:,:,:)

contains

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine zonal_fft_reg
    use cam_history_support, only: add_hist_coord

    nftnum = nlons/2+1

    ! add history coordinate for FT number
    call add_hist_coord('fft_num', nftnum, 'Fourier Transform Number')

  end subroutine zonal_fft_reg

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine zonal_fft_init(ntime_in)
    use cam_history, only: addfld

    integer, intent(in) :: ntime_in

    call addfld('U_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of FT U', gridname='esmf_zonal_mean')
    call addfld('U_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of FT U', gridname='esmf_zonal_mean')

    call addfld('V_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of FT U', gridname='esmf_zonal_mean')
    call addfld('V_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of FT U', gridname='esmf_zonal_mean')

    call addfld('OMEGA_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of FT omega', gridname='esmf_zonal_mean')
    call addfld('OMEGA_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of FT omega', gridname='esmf_zonal_mean')
    call addfld('THETA_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of FT theta', gridname='esmf_zonal_mean')
    call addfld('THETA_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of FT theta', gridname='esmf_zonal_mean')
    call addfld('WSTAR_FFT_real', (/'fft_num','lev    '/), 'I', '1', 'Real part of the complex conjugat of FT omega', gridname='esmf_zonal_mean')
    call addfld('WSTAR_FFT_imag', (/'fft_num','lev    '/), 'I', '1', 'Imaginary part of the complex conjugat of FT omega', gridname='esmf_zonal_mean')


    call addfld('UWSTAR_cosp', (/'fft_num','lev    '/), 'I', '1', 'U*WSTAR cospectra', gridname='esmf_zonal_mean')
    call addfld('VWSTAR_cosp', (/'fft_num','lev    '/), 'I', '1', 'V*WSTAR cospectra', gridname='esmf_zonal_mean')
    call addfld('TWSTAR_cosp', (/'fft_num','lev    '/), 'I', '1', 'THETA*WSTAR cospectra', gridname='esmf_zonal_mean')

    call addfld('RHOBAR',  (/'lev'/), 'I', 'kg/m3', 'Zonal mean air density', gridname='esmf_zonal_mean')

    call addfld('MFLXXUP', (/'lev'/), 'I', '1', 'Positive unresolved zonal momentum flux', gridname='esmf_zonal_mean')
    call addfld('MFLXXUN', (/'lev'/), 'I', '1', 'Negative unresolved zonal momentum flux', gridname='esmf_zonal_mean')
    call addfld('MFLXYUP', (/'lev'/), 'I', '1', 'Positive unresolved meridianal momentum flux', gridname='esmf_zonal_mean')
    call addfld('MFLXYUN', (/'lev'/), 'I', '1', 'Negative unresolved meridianal momentum flux', gridname='esmf_zonal_mean')

    call addfld('FLXXR', (/'lev'/), 'I', '1', 'Resolved zonal momentum flux', gridname='esmf_zonal_mean')
    call addfld('FLXXU', (/'lev'/), 'I', '1', 'Unresolved zonal momentum flux', gridname='esmf_zonal_mean')
    call addfld('FLXYR', (/'lev'/), 'I', '1', 'Resolved meridianal momentum flux', gridname='esmf_zonal_mean')
    call addfld('FLXYU', (/'lev'/), 'I', '1', 'Unresolved meridianal momentum flux', gridname='esmf_zonal_mean')

    ntime = ntime_in
    allocate(accum_cospectra_u( nftnum, lat_beg:lat_end, pver, ntime ))
    accum_cospectra_u = 0._r8
    allocate(accum_cospectra_v( nftnum, lat_beg:lat_end, pver, ntime ))
    accum_cospectra_v = 0._r8

  end subroutine zonal_fft_init

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine zonal_fft_calc(phys_state)
    use cam_history, only: outfld
    use perf_mod, only: t_startf, t_stopf
    use cospext_mod, only: cospext
    use esmf_zonal_ops, only: glats
    use air_composition, only: rairv  ! composition dependent gas constant (J/K/kg)
    use ref_pres, only: pref_mid

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    complex(C_DOUBLE_COMPLEX) :: u_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: v_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: w_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: t_fft(nftnum, lat_beg:lat_end, pver)

    complex(C_DOUBLE_COMPLEX) :: wstar(nftnum, lat_beg:lat_end, pver)

    real(r8) :: ufld(pver,pcols,begchunk:endchunk)
    real(r8) :: vfld(pver,pcols,begchunk:endchunk)
    real(r8) :: wfld(pver,pcols,begchunk:endchunk)
    real(r8) :: tfld(pver,pcols,begchunk:endchunk)
    integer :: lchnk, ncol, icol

    integer :: i,k, n

    complex(r8) :: tmpfld(nftnum, lat_beg:lat_end, pver)
    real(r8) :: cospectra(nftnum, lat_beg:lat_end, pver)
    real(r8) :: latrad(lat_beg:lat_end)
    real(r8) :: rho(pver,pcols,begchunk:endchunk) ! air mass density
    real(r8) :: rhobar(lat_beg:lat_end, pver)

    real(r8) :: mflxxup(lat_beg:lat_end,pver), mflxxun(lat_beg:lat_end,pver)
    real(r8) :: mflxyup(lat_beg:lat_end,pver), mflxyun(lat_beg:lat_end,pver)
    real(r8) :: flxxr(lat_beg:lat_end,pver), flxxu(lat_beg:lat_end,pver)
    real(r8) :: flxyr(lat_beg:lat_end,pver), flxyu(lat_beg:lat_end,pver)

    real(r8) :: wvlxbeg, wvlxend
    real(r8) :: mflux_glb(nlats,pver)

    real(r8), parameter :: pi = 4._r8*atan(1._r8)
    real(r8), parameter :: deg2rad = pi/180._r8
    character(len=*), parameter :: subname  = 'zonal_fft_calc'

    wvlxbeg = 200.e3_r8  ! 200 km
    wvlxend = 20.e3_r8   ! 20 km

    call t_startf ('zonal_fft_calc')

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

    rhobar = esmf_zonal_mean_3d(rho)

    u_fft = esmf_zonal_fft_3d(ufld)
    call output_fld(u_fft, name='U')

    v_fft = esmf_zonal_fft_3d(vfld)
    call output_fld(v_fft, name='V')

    w_fft = esmf_zonal_fft_3d(wfld)
    call output_fld(w_fft, name='OMEGA')

    t_fft = esmf_zonal_fft_3d(tfld)
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

    do icol = lat_beg, lat_end
       call outfld('RHOBAR', rhobar(icol,:),1,icol)
    end do

    ! zonal component
    call cospext(nftnum, lat_beg,lat_end, pver,ntime, latrad, accum_cospectra_u, wvlxbeg,wvlxend, pref_mid, rhobar, mflxxup,mflxxun, flxxr,flxxu)

    ! meridianal component
    call cospext(nftnum, lat_beg,lat_end, pver,ntime, latrad, accum_cospectra_v, wvlxbeg,wvlxend, pref_mid, rhobar, mflxyup,mflxyun, flxyr,flxyu)

    do icol = lat_beg, lat_end
       call outfld('MFLXXUP', mflxxup(icol,:),1,icol)
       call outfld('MFLXXUN', mflxxun(icol,:),1,icol)
       call outfld('MFLXYUP', mflxyup(icol,:),1,icol)
       call outfld('MFLXYUN', mflxyun(icol,:),1,icol)

       call outfld('FLXXR', flxxr(icol,:),1,icol)
       call outfld('FLXXU', flxxu(icol,:),1,icol)
       call outfld('FLXYR', flxyr(icol,:),1,icol)
       call outfld('FLXYU', flxyu(icol,:),1,icol)

    end do

    ! gather and smooth

    call t_stopf ('zonal_fft_calc')

  contains

    function gather_fluxes( flx_loc ) result(flxglb)
      use mpi, only: MPI_REAL8, MPI_SUCCESS, MPI_SUM

      real(r8),intent(in) :: flx_loc(lat_beg:lat_end,1:pver)

      real(r8) :: flxglb(nlats,pver)
      real(r8) :: sndbuf(nlats,pver)
      integer :: rc, len

      len = nlats*pver

      flxglb = 0._r8
      sndbuf = 0._r8
      sndbuf(lat_beg:lat_end,1:pver) = flx_loc(lat_beg:lat_end,1:pver)

      call mpi_allreduce(sndbuf,flxglb,len,MPI_REAL8,MPI_SUM,zonal_comm,rc)
      if ( rc /= MPI_SUCCESS ) then
         call endrun('zonal_fft_mod::gather_fluxes: mpi_allreduce FAILED')
      end if

    end function gather_fluxes


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

  end subroutine zonal_fft_calc


end module zonal_fft_mod
