module esmf_zonal_mean_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils, only: masterproc

  implicit none

  integer :: zonal_mean_nlats = 0

contains

  subroutine esmf_zonal_mean_reg

    use esmf_lonlat_grid_mod, only: esmf_lonlat_grid_init
    use esmf_phys_mesh_mod, only: esmf_phys_mesh_init
    use esmf_phys2lonlat_mod, only: esmf_phys2lonlat_init

    zonal_mean_nlats = 90

    call esmf_lonlat_grid_init(zonal_mean_nlats)
    call esmf_phys_mesh_init()
    call esmf_phys2lonlat_init()

  end subroutine esmf_zonal_mean_reg

  subroutine esmf_zonal_mean_calc(lonlatarr, zmarr)
    use ppgrid, only: pver
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end, nlon
    use esmf_lonlat_grid_mod, only: zonal_comm
    use shr_reprosum_mod,only: shr_reprosum_calc

    real(r8), intent(in) :: lonlatarr(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8), intent(out) :: zmarr(lat_beg:lat_end,pver)

    real(r8) :: tmparr(lon_beg:lon_end,pver)
    real(r8) :: gsum(pver)

    integer :: numlons, ilat

    numlons = lon_end-lon_beg+1

    ! zonal mean

    do ilat = lat_beg, lat_end
       tmparr(lon_beg:lon_end,:) = lonlatarr(lon_beg:lon_end,ilat,:)
       call shr_reprosum_calc(tmparr, gsum, numlons, numlons, pver, gbl_count=nlon, commid=zonal_comm)
       zmarr(ilat,:) = gsum(:)/nlon
    end do

  end subroutine esmf_zonal_mean_calc

end module esmf_zonal_mean_mod
