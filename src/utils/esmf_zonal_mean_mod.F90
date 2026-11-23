!------------------------------------------------------------------------------
! Provides methods for calculating zonal means on the ESMF regular longitude
! / latitude grid
!------------------------------------------------------------------------------
module esmf_zonal_mean_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils, only: masterproc
  use cam_history_support, only : fillvalue

  implicit none

  private

  public :: esmf_zonal_mean_calc
  public :: esmf_zonal_mean_masked
  public :: esmf_zonal_mean_wsums

  interface esmf_zonal_mean_calc
     module procedure esmf_zonal_mean_calc_2d
     module procedure esmf_zonal_mean_calc_3d
  end interface esmf_zonal_mean_calc

contains

  !------------------------------------------------------------------------------
  ! Calculates zonal means of 3D fields. The wght option can be used to mask out
  ! regions such as mountains.
  !------------------------------------------------------------------------------
  subroutine esmf_zonal_mean_calc_3d(lonlatarr, zmarr, wght)
    use ppgrid, only: pver
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end, nlon
    use esmf_lonlat_grid_mod, only: zonal_comm
    use shr_reprosum_mod,only: shr_reprosum_calc

    real(r8), intent(in) :: lonlatarr(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8), intent(out) :: zmarr(lat_beg:lat_end,pver)

    real(r8), optional, intent(in) :: wght(lon_beg:lon_end,lat_beg:lat_end,pver)

    real(r8) :: wsums(lat_beg:lat_end,pver)

    ! zonal mean
    if (present(wght)) then

       wsums = esmf_zonal_mean_wsums(wght)
       call esmf_zonal_mean_masked(lonlatarr, wght, wsums, zmarr)

    else

       ! zonal sums of all (lat,lev) columns in one reproducible-sum call
       zmarr = zonal_sums_3d(lonlatarr)/nlon

    end if

  end subroutine esmf_zonal_mean_calc_3d

  !------------------------------------------------------------------------------
  ! Computes zonal mean for 2D lon / lat fields.
  !------------------------------------------------------------------------------
  subroutine esmf_zonal_mean_calc_2d(lonlatarr, zmarr)
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end, nlon
    use esmf_lonlat_grid_mod, only: zonal_comm
    use shr_reprosum_mod,only: shr_reprosum_calc

    real(r8), intent(in) :: lonlatarr(lon_beg:lon_end,lat_beg:lat_end)
    real(r8), intent(out) :: zmarr(lat_beg:lat_end)

    real(r8) :: gsum(lat_beg:lat_end)

    integer :: numlons, numlats

    numlons = lon_end-lon_beg+1
    numlats = lat_end-lat_beg+1

    ! zonal mean

    call shr_reprosum_calc(lonlatarr, gsum, numlons, numlons, numlats, gbl_count=nlon, commid=zonal_comm)
    zmarr(:) = gsum(:)/nlon

  end subroutine esmf_zonal_mean_calc_2d

  !------------------------------------------------------------------------------
  ! Computes longitude sums of grid cell weights.
  !------------------------------------------------------------------------------
  function esmf_zonal_mean_wsums(wght) result(wsums)
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end, nlon
    use esmf_lonlat_grid_mod, only: zonal_comm
    use shr_reprosum_mod,only: shr_reprosum_calc
    use ppgrid, only: pver

    real(r8), intent(in) :: wght(lon_beg:lon_end,lat_beg:lat_end,pver)

    real(r8) :: wsums(lat_beg:lat_end,pver)

    wsums = zonal_sums_3d(wght)

  end function esmf_zonal_mean_wsums

  !------------------------------------------------------------------------------
  ! Longitude sums of a 3D lon/lat/lev field for every (lat,lev) column, done
  ! with a single reproducible-sum call over zonal_comm (each (lat,lev) column
  ! is a separate field). This replaces one collective call per latitude; the
  ! reproducible sum of each field is independent of the other fields in the
  ! call, so the results are unchanged.
  !------------------------------------------------------------------------------
  function zonal_sums_3d(arr) result(zsums)
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end, nlon
    use esmf_lonlat_grid_mod, only: zonal_comm
    use shr_reprosum_mod,only: shr_reprosum_calc
    use ppgrid, only: pver

    real(r8), intent(in) :: arr(lon_beg:lon_end,lat_beg:lat_end,pver)

    real(r8) :: zsums(lat_beg:lat_end,pver)

    integer :: numlons, nflds
    real(r8) :: arr2d(lon_end-lon_beg+1,(lat_end-lat_beg+1)*pver)
    real(r8) :: gsum((lat_end-lat_beg+1)*pver)

    numlons = lon_end-lon_beg+1
    nflds = (lat_end-lat_beg+1)*pver

    if (nflds < 1) return

    ! (lon,lat,lev) -> (lon, field) with field = ilat-lat_beg+1 + (ilev-1)*nlat_local
    arr2d = reshape(arr, shape(arr2d))

    call shr_reprosum_calc(arr2d, gsum, numlons, numlons, nflds, gbl_count=nlon, commid=zonal_comm)

    zsums = reshape(gsum, shape(zsums))

  end function zonal_sums_3d

  !------------------------------------------------------------------------------
  ! Masks out regions (e.g. mountains) from zonal mean calculation.
  !------------------------------------------------------------------------------
  subroutine esmf_zonal_mean_masked(lonlatarr, wght, wsums, zmarr)
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end, nlon
    use esmf_lonlat_grid_mod, only: zonal_comm
    use shr_reprosum_mod,only: shr_reprosum_calc
    use ppgrid, only: pver

    real(r8), intent(in) :: lonlatarr(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8), intent(in) :: wght(lon_beg:lon_end,lat_beg:lat_end,pver) ! grid cell weights
    real(r8), intent(in) :: wsums(lat_beg:lat_end,pver) ! pre-computed sums of grid cell weights
    real(r8), intent(out) :: zmarr(lat_beg:lat_end,pver) ! zonal means

    real(r8) :: tmparr(lon_beg:lon_end,pver)
    integer :: numlons, ilat, ilev
    real(r8) :: gsum(pver)

    numlons = lon_end-lon_beg+1

    do ilat = lat_beg, lat_end

       tmparr(lon_beg:lon_end,:) = wght(lon_beg:lon_end,ilat,:)*lonlatarr(lon_beg:lon_end,ilat,:)
       call shr_reprosum_calc(tmparr, gsum, numlons, numlons, pver, gbl_count=nlon, commid=zonal_comm)

       do ilev = 1,pver
          if (wsums(ilat,ilev)>0._r8) then
             zmarr(ilat,ilev) = gsum(ilev)/wsums(ilat,ilev)
          else
             zmarr(ilat,ilev) = fillvalue
          end if
       end do

    end do

  end subroutine esmf_zonal_mean_masked

end module esmf_zonal_mean_mod
