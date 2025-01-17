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

  end subroutine zonal_fft_init

  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  subroutine zonal_fft_calc(phys_state)
    use cam_history, only: outfld

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    complex(C_DOUBLE_COMPLEX) :: u_fft(nftnum, lat_beg:lat_end, pver)
    complex(C_DOUBLE_COMPLEX) :: v_fft(nftnum, lat_beg:lat_end, pver)
    real(r8) :: ufld(pver,pcols,begchunk:endchunk)
    real(r8) :: vfld(pver,pcols,begchunk:endchunk)
    integer :: lchnk, ncol, icol
    real(r8) ::  tmpr(nftnum, pver)
    real(r8) ::  tmpi(nftnum, pver)

    integer :: n,k

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       do icol = 1,ncol
          ufld(:pver,icol,lchnk) = phys_state(lchnk)%u(icol,:pver)
          vfld(:pver,icol,lchnk) = phys_state(lchnk)%v(icol,:pver)
       end do
    end do

    u_fft = esmf_zonal_fft_3d(ufld)
    v_fft = esmf_zonal_fft_3d(vfld)

    do icol = lat_beg, lat_end
       do n = 1,nftnum
          do k = 1,pver
             tmpr(n,k) = u_fft(n,icol,k)%re
             tmpi(n,k) = u_fft(n,icol,k)%im
          end do
       end do

       call outfld('U_FFT_real', tmpr, 1, icol)
       call outfld('U_FFT_imag', tmpi, 1, icol)

       do n = 1,nftnum
          do k = 1,pver
             tmpr(n,k) = v_fft(n,icol,k)%re
             tmpi(n,k) = v_fft(n,icol,k)%im
          end do
       end do

       call outfld('V_FFT_real', tmpr, 1, icol)
       call outfld('V_FFT_imag', tmpi, 1, icol)
    end do


  end subroutine zonal_fft_calc


end module zonal_fft_mod
