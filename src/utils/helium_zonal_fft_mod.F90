!-------------------------------------------------------------------------------
! Utility module for Helium fluxes at the upper boundary of WACCM-X
!-------------------------------------------------------------------------------
module helium_zonal_fft_mod

  use shr_kind_mod, only: r8 => shr_kind_r8, cl=>SHR_KIND_CL
  use ppgrid, only: pcols, pver, begchunk, endchunk
  use esmf_lonlat_grid_mod, only: esmf_lonlat_grid

  use mpi, only: MPI_REAL8, MPI_SUCCESS, MPI_SUM
  use perf_mod, only: t_startf, t_stopf
  use cam_abortutils, only: endrun

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  private
  public :: helium_zonal_fft_init
  public :: helium_zonal_fft_forward
  public :: helium_zonal_fft_backward

  type(C_PTR) :: plan_forward
  type(C_PTR) :: plan_backward

  type(esmf_lonlat_grid), pointer :: he_grid => null()

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine helium_zonal_fft_init(grid)

    class(esmf_lonlat_grid), pointer, intent(in) :: grid

    complex(C_DOUBLE_COMPLEX),dimension(:),allocatable :: zin,zout

    he_grid => grid

    allocate(zin(grid%nlon))
    allocate(zout(grid%nlon))

    ! Create a forward FFT plan -- real to complex 1-dimensional
    plan_forward = fftw_plan_dft_1d(grid%nlon,zin,zout,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_backward = fftw_plan_dft_1d(grid%nlon,zout,zin,FFTW_BACKWARD,FFTW_ESTIMATE)

    deallocate(zin)
    deallocate(zout)

  end subroutine helium_zonal_fft_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function helium_zonal_fft_forward(fld_lonlat) result(fld_out)

    real(r8),intent(in) :: fld_lonlat(he_grid%lon_beg:he_grid%lon_end,he_grid%lat_beg:he_grid%lat_end)

    complex(r8) :: fld_out(he_grid%nlon, he_grid%lat_beg:he_grid%lat_end)

    complex(C_DOUBLE_COMPLEX) :: zout(he_grid%nlon)
    complex(C_DOUBLE_COMPLEX) :: zin(he_grid%nlon)

    integer :: rc, ilat, len

    real(r8) :: sndbf(he_grid%nlon,he_grid%lat_beg:he_grid%lat_end)
    real(r8) :: rcvbf(he_grid%nlon,he_grid%lat_beg:he_grid%lat_end)

    character(len=*), parameter :: subname = 'helium_zonal_fft_forward: '

    call t_startf('helium_zonal_fft_forward')

    ! zonal FFT

    ! gather all longitudes ....

    associate( nlon=>he_grid%nlon, &
         lat_beg=>he_grid%lat_beg, &
         lat_end=>he_grid%lat_end, &
         lon_beg=>he_grid%lon_beg, &
         lon_end=>he_grid%lon_end, &
         zonal_comm=>he_grid%zonal_comm )

    len = nlon*(lat_end-lat_beg+1)

    rcvbf(:,:) = 0._r8
    sndbf(:,:) = 0._r8
    sndbf(lon_beg:lon_end,lat_beg:lat_end) = fld_lonlat(lon_beg:lon_end,lat_beg:lat_end)
    call mpi_allreduce( sndbf, rcvbf, len, MPI_REAL8, MPI_SUM, zonal_comm, rc )
    if ( rc /= MPI_SUCCESS ) then
       call endrun(subname//'mpi_allreduce failed 2')
    end if

    do ilat = lat_beg, lat_end
       zin(:) = rcvbf(:,ilat)
       call fftw_execute_dft(plan_forward, zin, zout)
       fld_out(:,ilat) = zout(:)
    end do

    end associate

    call t_stopf('helium_zonal_fft_forward')

  end function helium_zonal_fft_forward

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function helium_zonal_fft_backward(fld_in) result(fld_lonlat)

    complex(r8), intent(in) :: fld_in(he_grid%nlon, he_grid%lat_beg:he_grid%lat_end)

    complex(r8) :: fld_lonlat(he_grid%nlon,he_grid%lat_beg:he_grid%lat_end)

    complex(C_DOUBLE_COMPLEX) :: zin(he_grid%nlon)
    complex(C_DOUBLE_COMPLEX) :: zout(he_grid%nlon)

    integer :: ilat

    associate( &
         lat_beg=>he_grid%lat_beg, &
         lat_end=>he_grid%lat_end )

    do ilat = lat_beg, lat_end
       zin(:) = fld_in(:,ilat)
       call fftw_execute_dft(plan_backward, zin, zout)
       fld_lonlat(:,ilat) = zout(:)
    end do

    end associate

  end function helium_zonal_fft_backward


end module helium_zonal_fft_mod
