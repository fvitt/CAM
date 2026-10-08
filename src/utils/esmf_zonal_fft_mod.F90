module esmf_zonal_fft_mod

  use shr_kind_mod, only: r8 => shr_kind_r8, cl=>SHR_KIND_CL
  use ppgrid, only: pcols, pver, begchunk, endchunk
  use esmf_lonlat_grid_mod, only: nlon, lon_beg, lon_end, lat_beg, lat_end
  use esmf_lonlat_grid_mod, only: zonal_comm

  use mpi, only: MPI_REAL8, MPI_SUCCESS, MPI_SUM
  use perf_mod, only: t_startf, t_stopf
  use cam_abortutils, only: endrun

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  type(C_PTR) :: fftw_plan
  real(C_DOUBLE), allocatable :: fftw_in(:)
  complex(C_DOUBLE_COMPLEX), allocatable :: fftw_out(:)

  ! inverse (complex -> real) plan and scratch arrays
  type(C_PTR) :: fftw_iplan
  complex(C_DOUBLE_COMPLEX), allocatable :: fftw_iin(:)
  real(C_DOUBLE), allocatable :: fftw_iout(:)

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_zonal_fft_init

    ! FFTW3 variables
    integer(C_INTPTR_T) :: fftw_n
    type(C_PTR) :: plan
    integer(C_INTPTR_T) :: local_n, local_start
    integer(C_INTPTR_T) :: local_ni, local_i_start, local_no, local_o_start

    allocate(fftw_in(nlon))
    allocate(fftw_out(nlon/2+1))

    ! Create a forward FFT plan -- real to complex 1-dimensional
    fftw_plan = fftw_plan_dft_r2c_1d(nlon, fftw_in, fftw_out, FFTW_ESTIMATE)

    ! Create an inverse FFT plan -- complex to real 1-dimensional
    allocate(fftw_iin(nlon/2+1))
    allocate(fftw_iout(nlon))
    fftw_iplan = fftw_plan_dft_c2r_1d(nlon, fftw_iin, fftw_iout, FFTW_ESTIMATE)

  end subroutine esmf_zonal_fft_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function esmf_zonal_fft_3d(fld_lonlat) result(zfft)

    real(r8),intent(in) :: fld_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)

    complex(C_DOUBLE_COMPLEX) :: zfft(nlon/2+1, lat_beg:lat_end, pver)

    integer :: rc, i, ichnk, icol, ilon, ilat, ilev, ncol, len

    real(C_DOUBLE) :: fld(nlon)

    real(r8) :: sndbf(nlon)
    real(r8) :: rcvbf(nlon)

    real(r8) :: sndbf2(nlon,lat_beg:lat_end,1:pver)
    real(r8) :: rcvbf2(nlon,lat_beg:lat_end,1:pver)

    character(len=*), parameter :: subname = ': esmf_zonal_fft_3d'

    call t_startf('esmf_zonal_fft_3d')

    ! zonal FFT


    ! gather all longitudes ....

    len = pver*nlon*(lat_end-lat_beg+1)

    sndbf2(:,:,:) = 0._r8
    sndbf2(lon_beg:lon_end,lat_beg:lat_end,1:pver) = fld_lonlat(lon_beg:lon_end,lat_beg:lat_end,1:pver)
    call mpi_allreduce( sndbf2, rcvbf2, len, MPI_REAL8, MPI_SUM, zonal_comm, rc )
    if ( rc /= MPI_SUCCESS ) then
       call endrun(subname//'mpi_allreduce failed 2')
    end if

    do ilev = 1, pver
       do ilat = lat_beg, lat_end
          fftw_in(:) = real(rcvbf2(:,ilat,ilev),kind=C_DOUBLE)
          call fftw_execute_dft_r2c(fftw_plan, fftw_in, fftw_out)
          zfft(:,ilat,ilev) = fftw_out(:)/nlon ! normalize
       end do
    end do

    call t_stopf('esmf_zonal_fft_3d')

  end function esmf_zonal_fft_3d

  !------------------------------------------------------------------------------
  ! Inverse zonal FFT: complex spectrum (nlon/2+1, lat_beg:lat_end, pver) ->
  ! real physical field (lon_beg:lon_end, lat_beg:lat_end, pver).
  ! Assumes the input spectrum is normalized as produced by esmf_zonal_fft_3d.
  ! Each rank holds full-longitude spectra for its local latitudes, so no
  ! interprocessor communication is required.
  !------------------------------------------------------------------------------
  function esmf_zonal_ifft_3d(zfft) result(fld_lonlat)

    complex(C_DOUBLE_COMPLEX), intent(in) :: zfft(nlon/2+1, lat_beg:lat_end, pver)

    real(r8) :: fld_lonlat(lon_beg:lon_end, lat_beg:lat_end, pver)

    integer :: ilev, ilat
    character(len=*), parameter :: subname = ': esmf_zonal_ifft_3d'

    call t_startf('esmf_zonal_ifft_3d')

    do ilev = 1, pver
       do ilat = lat_beg, lat_end
          ! c2r transform destroys its input; copy into the plan's buffer
          fftw_iin(:) = zfft(:,ilat,ilev)
          call fftw_execute_dft_c2r(fftw_iplan, fftw_iin, fftw_iout)
          fld_lonlat(lon_beg:lon_end, ilat, ilev) = real(fftw_iout(lon_beg:lon_end), kind=r8)
       end do
    end do

    call t_stopf('esmf_zonal_ifft_3d')

  end function esmf_zonal_ifft_3d

end module esmf_zonal_fft_mod
