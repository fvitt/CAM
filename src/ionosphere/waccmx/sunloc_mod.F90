module sunloc_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi

  use time_manager, only: get_curr_date

  implicit none

contains
  !-----------------------------------------------------------------------  !
  ! Return sun's longitude in dipole coordinates
  !
  subroutine sunloc_calc(sunlon)
    use getapex,      only: alonm ! (nlonp1,0:nlatp1)
    use edyn_geogrid, only: nlon, nlat, dphi, dlamda ! oplus grid

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


  !-----------------------------------------------------------------------
  subroutine sunloc_calc2(sunlon)
    use apex,only:apex_mall
    use cons_module, only: rtd, dtr

    real(r8),intent(out) :: sunlon ! magnetic longitude

    real(r8) :: glats, glons
    real(r8) :: alt, hr
    real(r8) :: &
         b(3)             ,& ! Magnetic field components (east, north, up), in nT
         bhat(3)          ,& ! components (east, north, up) of unit vector along
                             ! geomagnetic field direction
         bmag             ,& ! Magnitude of magnetic field (nT)
         si               ,& ! sin(i)
         alon             ,& ! Apex longitude = modified apex longitude =
                             ! quasi-dipole longitude (deg)
         xlatm            ,& ! Modified Apex latitude (deg)
         vmp              ,& ! Magnetic potential (T.m)
         w                ,& ! W of Richmond reference above, in km**2 /nT (i.e., 10**15 m**2 /T)
         d                ,& ! D of Richmond reference above
         be3              ,& ! B_e3 of reference above (= Bmag/D), in nT
         sim              ,& ! sin(I_m) described in Richmond reference above
         xlatqd           ,& ! Quasi-dipole latitude (deg)
         f                   ! F described in ref above for quasi-dipole coordinates

    real(r8),dimension(3) :: d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3 ! Components of base vectors
    integer :: ier ! error return
    integer :: yr, mon, day, tod

    call get_curr_date(yr, mon, day, tod)

    ! Sun's geographic coordinates:
    glats = asin(.398749_r8*sin(2._r8 * pi * real(day-80, r8) / 365._r8))
    glons = pi * (1._r8 - (2._r8 * real(tod, r8) / 86400._r8))

    alt = 100._r8
    hr = 100._r8

    glats = glats * rtd
    glons = glons * rtd

    call apex_mall(glats,glons,alt,hr, b,bhat,bmag,si,alon,xlatm,vmp,w,&
         d,be3,sim,d1,d2,d3,e1,e2,e3,xlatqd,f,f1,f2,f3,g1,g2,g3,ier)

    sunlon = alon * dtr
  end subroutine sunloc_calc2

end module sunloc_mod
