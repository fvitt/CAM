module zonal_fft_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver, begchunk, endchunk
  use phys_grid, only: get_ncols_p
  use physics_types, only: physics_state
  use esmf_zonal_ops, only : lat_beg,lat_end, lon_beg,lon_end, nlons, glats, nlats
  use esmf_zonal_ops, only : esmf_zonal_fft_3d
  use, intrinsic :: iso_c_binding

  implicit none

  integer :: nftnum = 0

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
  subroutine zonal_fft_init
    use cam_history, only: addfld

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

  end subroutine zonal_fft_init

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine zonal_fft_calc(phys_state)
    use cam_history, only: outfld

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

    integer :: n,k

    complex(r8) :: tmpfld(nftnum, lat_beg:lat_end, pver)
    real(r8) :: cospectra(nftnum, lat_beg:lat_end, pver)

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       do icol = 1,ncol
          ufld(:pver,icol,lchnk) = phys_state(lchnk)%u(icol,:pver)
          vfld(:pver,icol,lchnk) = phys_state(lchnk)%v(icol,:pver)
          wfld(:pver,icol,lchnk) = phys_state(lchnk)%omega(icol,:pver)
          ! = -sheight(:ncol,:) *  phys_state(lchnk)%omega(:ncol,:) / phys_state(lchnk)%pmid(:ncol,:)
          tfld(:pver,icol,lchnk) = phys_state(lchnk)%t(icol,:pver) * phys_state(lchnk)%exner(icol,:pver)
       end do
    end do

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

    tmpfld = u_fft * wstar
    cospectra = wstar%re
    call output_cosp(cospectra,'U')

    tmpfld = v_fft * wstar
    cospectra = wstar%re
    call output_cosp(cospectra,'V')

    tmpfld = t_fft * wstar
    cospectra = wstar%re
    call output_cosp(cospectra,'T')

  contains

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
