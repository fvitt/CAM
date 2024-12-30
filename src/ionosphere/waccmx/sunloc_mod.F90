module sunloc_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi

  use time_manager, only: get_curr_date
  use getapex,      only: alonm ! (nlonp1,0:nlatp1)
  use edyn_geogrid, only: nlon, nlat, dphi, dlamda ! oplus grid

  implicit none

contains
  !-----------------------------------------------------------------------
  !
  ! Return sun's longitude in dipole coordinates
  !
  subroutine sunloc_calc(sunlon)

    !
    ! Args:
    real(r8),intent(out) :: sunlon ! magnetic longitude
    !
    ! Local:
    integer :: j, i, ii, isun, jsun
    real(r8) :: glats, glons, pisun, pjsun, sndlons, csdlons
    real(r8) :: rlonm(nlon+4, nlat) ! (nlon+4,nlat)
    real(r8) :: r8_isun, r8_jsun

    integer :: yr, mon, day, tod

    call get_curr_date(yr, mon, day, tod)

    ! Sun's geographic coordinates:
    glats = asin(.398749_r8*sin(2._r8 * pi * real(day-80, r8) / 365._r8))
    glons = pi * (1._r8 - (2._r8 * real(tod, r8) / 86400._r8))

    do j = 1, nlat
       do i = 1, nlon
          ii = i + 2
          rlonm(ii, j) = alonm(i, j)
       end do
       do i = 1, 2
          rlonm(i, j) = rlonm(i+nlon, j)
          rlonm(i+nlon+2, j) = rlonm(i+2, j)
       end do
    end do

    pisun = ((glons + pi) / dlamda) + 1._r8
    pjsun = ((glats + (.5_r8 * (pi - dphi))) / dphi) + 1._r8
    isun = int(pisun)
    jsun = int(pjsun)
    r8_isun = real(isun, r8)
    r8_jsun = real(jsun, r8)
    pisun = pisun - r8_isun
    pjsun = pjsun - r8_jsun

    sndlons = &
         (1._r8-pisun) * (1._r8-pjsun) * sin(rlonm(isun+2,jsun  )) + &
                pisun  * (1._r8-pjsun) * sin(rlonm(isun+3,jsun  )) + &
                pisun  *        pjsun  * sin(rlonm(isun+3,jsun+1)) + &
         (1._r8-pisun) *        pjsun  * sin(rlonm(isun+2,jsun+1))

    csdlons = &
         (1._r8-pisun) * (1._r8-pjsun) * cos(rlonm(isun+2,jsun)) +       &
                pisun  * (1._r8-pjsun) * cos(rlonm(isun+3,jsun))+        &
                pisun  *        pjsun  * cos(rlonm(isun+3,jsun+1))+      &
         (1._r8-pisun) *        pjsun  * cos(rlonm(isun+2,jsun+1))

    sunlon = atan2(sndlons, csdlons)

  end subroutine sunloc_calc

end module sunloc_mod
