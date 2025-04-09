module grid_module

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  subroutine generate_mag_grid( edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt )
! new grid since 2015/02 (see Art's notes 2015/02/12), updated 2015/04
! reference height at k=0.5, grids start at 80 km
! with closer latitude spacing at low latitudes and in the auroral region

    use params_module,only:nmlon,nmlat_h,nmlatS2_h,nmlat_T1,nmlat_T2,nhgt_fix,nhgt_fix_r, &
      ylonm,ylonm_s,ylatm,ylatm_s,rho,rho_s,ha,ha_s,hgt_fix,hgt_fix_r
    use cons_module,only:re,h0,r0,pi,dtr,fill_value

    integer, intent(in) :: edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt

    real(kind=rp),parameter :: rekm = re*1e-3_rp, h0km = h0*1e-3_rp, r0km = r0*1e-3_rp, &
      d=10, f=15, h=10, yb=5, yc=18, d1=30, d2=50, d3=55, d4=75, d5=82, hs=6, dhdz=6
    integer :: i,j,k,jns
    real(kind=rp) :: lam1,lam2,lam3,lam4,lam5,hc,h1,rho1,b,c,e,g, &
      y1,y2,y3,y4,y5,ymax,y,lam,rho_loc,ha_loc,dlonm

    ! set run-time grid resolution parameters
    nmlon = edyn3d_nmlon
    nmlat_h = edyn3d_nmlat_h
    nmlatS2_h = nmlat_h-1
    nhgt_fix = edyn3d_nhgt
    nhgt_fix_r = nhgt_fix+1

    allocate(ylonm(0:nmlon+1))
    allocate(ylonm_s(0:nmlon+1))
    allocate(ylatm(2,nmlat_h))
    allocate(ylatm_s(2,nmlatS2_h))
    allocate(rho(nmlat_h))
    allocate(rho_s(nmlatS2_h))
    allocate(ha(nmlat_h))
    allocate(ha_s(nmlatS2_h))
    allocate(hgt_fix(nhgt_fix))
    allocate(hgt_fix_r(nhgt_fix_r))

    lam1 = d1*dtr
    lam2 = d2*dtr
    lam3 = d3*dtr
    lam4 = d4*dtr
    lam5 = d5*dtr

    hc = h0km + hs*yc + dhdz*(yc-yb)**2/2
    h1 = r0km/cos(lam1)**2 - rekm
    rho1 = cos(lam1)

    c = hs + dhdz*(yc-yb)
    b = (2*r0km*sin(lam1)/(d*rho1**3)-c)/(r0km/rho1**2-rekm-hc)
    e = (f-d)/(2*(lam3-lam2))
    g = (f-h)/(2*(lam5-lam4))

    y1 = yc + log(b*(h1-hc)/c+1)/b
    y2 = y1 + d*(lam2-lam1)
    y3 = y2 + d*(lam3-lam2) + e*(lam3-lam2)**2
    y4 = y3 + f*(lam4-lam3)
    y5 = y4 + f*(lam5-lam4) - g*(lam5-lam4)**2
    ymax = y5 + h*(pi/2-lam5)

! ordering goes from equator to pole - P,S1
    do concurrent (j = 1:nmlat_h)
      y = (j-1)*ymax/(nmlat_h-1)

      if (y <= yb) then
        ha_loc = h0km + hs*y
        rho_loc = sqrt(r0km/(ha_loc+rekm))
        lam = acos(rho_loc)
      elseif (y <= yc) then
        ha_loc = h0km + hs*y + dhdz*(y-yb)**2/2
        rho_loc = sqrt(r0km/(ha_loc+rekm))
        lam = acos(rho_loc)
      elseif (y <= y1) then
        ha_loc = hc + c*(exp(b*(y-yc))-1)/b
        rho_loc = sqrt(r0km/(ha_loc+rekm))
        lam = acos(rho_loc)
      elseif (y <= y2) then
        lam = (y-y1)/d + lam1
        rho_loc = cos(lam)
        ha_loc = r0km/rho_loc**2 - rekm
      elseif (y <= y3) then
        lam = lam2 + (sqrt(d**2+4*e*(y-y2))-d)/(2*e)
        rho_loc = cos(lam)
        ha_loc = r0km/rho_loc**2 - rekm
      elseif (y <= y4) then
        lam = (y-y3)/f + lam3
        rho_loc = cos(lam)
        ha_loc = r0km/rho_loc**2 - rekm
      elseif (y <= y5) then
        lam = lam4 + (f-sqrt(f**2-4*g*(y-y4)))/(2*g)
        rho_loc = cos(lam)
        ha_loc = r0km/rho_loc**2 - rekm
      else
        lam = (y-y5)/h + lam5
        rho_loc = cos(lam)
        if (j == nmlat_h) then
          ha_loc = fill_value
        else
          ha_loc = r0km/rho_loc**2 - rekm
        endif
      endif

! jns: nmlat_h (equator) to 1 (pole)
      jns = nmlat_h-j+1

      ylatm(2,jns) = lam
      ylatm(1,jns) = -lam

! added 2015/04
! overwrite numerical inaccuracy in calculating rho for j=nmlat_h / jns=1
      if (jns == 1) then
        rho(jns) = 0
        ha(jns) = fill_value
      else
        rho(jns) = rho_loc
        ha(jns) = ha_loc
      endif
    enddo

! set up height levels - P,S1,R
! reference height at k=0.5
    do concurrent (k = 1:nhgt_fix)
      j = nmlat_h-k+1
      hgt_fix(k) = ha(j)
    enddo

! set S2 latitude and rho values
    do concurrent (j = 1:nmlatS2_h)
      rho_s(j) = (rho(j)+rho(j+1))/2
      ylatm_s(2,j) = acos(rho_s(j))
      ylatm_s(1,j) = -ylatm_s(2,j)

! calculate apex height of each field line, Richmond 1995 Eq (3.3)
! ha=r0/cos^2lambda-re with lambda modified apex latitude
      ha_s(j) = r0km/rho_s(j)**2 - rekm
    enddo

! set up height levels - R
! these are apex heights of ylatm_s
    hgt_fix_r(1) = h0km
    do concurrent (k = 2:nhgt_fix_r)
      j = nmlatS2_h-k+2
      hgt_fix_r(k) = ha_s(j)
    enddo

    nmlat_T1 = 2*nmlat_h-1
    nmlat_T2 = 2*nmlatS2_h

! convert from km to m
    do k = 1,nhgt_fix
      hgt_fix(k) = hgt_fix(k)*1000
    enddo
    do k = 1,nhgt_fix_r
      hgt_fix_r(k) = hgt_fix_r(k)*1000
    enddo
    do j = 1,nmlat_h
      ha(j) = ha(j)*1000
    enddo
    do j = 1,nmlatS2_h
      ha_s(j) = ha_s(j)*1000
    enddo

! magnetic longitudes are equidistant
    dlonm = 2*pi/nmlon
    do concurrent (i = 0:nmlon+1)
      ylonm(i) = -pi+(i-1)*dlonm    ! P,S2,R
      ylonm_s(i) = ylonm(i)+dlonm/2 ! S1
    enddo

  endsubroutine generate_mag_grid
!-----------------------------------------------------------------------
  subroutine generate_simplified_grid
! 5 parameter setup:
!   1. top height
!   2. height spacing (in scale height)
!   3. number of transition latitudes
!   4. high-latitude spacing
!   5. longitude spacing

    use params_module,only:nmlon,nmlat_h,nmlatS2_h,nmlat_T1,nmlat_T2,nhgt_fix,nhgt_fix_r, &
      ylonm,ylonm_s,ylatm,ylatm_s,rho,rho_s,ha,ha_s,hgt_fix,hgt_fix_r
    use cons_module,only:re,h0,r0,pi,dtr,rtd,fill_value

    real(kind=rp),parameter :: rekm = re*1e-3_rp, h0km = h0*1e-3_rp, r0km = r0*1e-3_rp, &
      a = 2, b = 1, htop = 1200, dh = 1.0_rp/3.0_rp, dlat_high = 1
    integer,parameter :: ntrans = 5
    integer :: i,j,k,nhigh
    real(kind=rp) :: dlat_low,ddlat,lat0_high,dlonm
    real(kind=rp),dimension(ntrans) :: dlat_trans,lat_trans
    real(kind=rp),dimension(:),allocatable :: lat_low,lat_high,lat

! given top height and height spacing, calculate the number of height grids
! heights are assumed to be parabolically increasing with log-pressure levels (scale height)
    nhgt_fix = 1 + (sqrt(b**2+4*a*(htop-h0km))-b)/(2*a*dh)

! init height grids and the footpoint latitudes (coupled height-latitude grid)
    allocate(hgt_fix(nhgt_fix))
    allocate(lat_low(nhgt_fix))
    do concurrent (k = 1:nhgt_fix)
      hgt_fix(k) = a*((k-1)*dh)**2 + b*(k-1)*dh + h0km
      lat_low(k) = acos(sqrt(r0km/(hgt_fix(k)+rekm)))*rtd
    enddo

! the latitude spacing at the low-latitude boundary
    dlat_low = lat_low(nhgt_fix)-lat_low(nhgt_fix-1)

! we want a smooth transition from low to high latitudes within ntrans grids

! given grid spacing is increasing linearly, find the spacing of each transition grid
    ddlat = (dlat_high-dlat_low)/(ntrans+1)
    do concurrent (j = 1:ntrans)
      dlat_trans(j) = dlat_low+j*ddlat
    enddo

! given the spacing of transition grids, find the latitude of transition grids
    lat_trans(1) = lat_low(nhgt_fix)+dlat_trans(1)
    do j = 2,ntrans
      lat_trans(j) = lat_trans(j-1)+dlat_trans(j)
    enddo

! the starting latitude of high-latitude grids
    lat0_high = lat_low(nhgt_fix) + ntrans*(dlat_low+dlat_high)/2

! the number of high-latitude grids
    nhigh = (90-lat0_high)/dlat_high

! the latitude of high-latitude grids
    allocate(lat_high(nhigh))
    do concurrent (j = 1:nhigh)
      lat_high(j) = 90 - (nhigh-j)*dlat_high
    enddo

! total number of latitude grids is the sum of
! low-latitude grids, transition grids, high-latitude grids
    nmlat_h = nhgt_fix+ntrans+nhigh

    allocate(lat(nmlat_h))
    lat(1:nhgt_fix) = lat_low
    lat(nhgt_fix+1:nhgt_fix+ntrans) = lat_trans
    lat(nhgt_fix+ntrans+1:nmlat_h) = lat_high

! given latitude grids, find other grid parameters
    allocate(ylatm(2,nmlat_h))
    allocate(rho(nmlat_h))
    allocate(ha(nmlat_h))
    do concurrent (j = 1:nmlat_h)
      k = nmlat_h-j+1
      ylatm(2,j) = lat(k)*dtr
      ylatm(1,j) = -ylatm(2,j)
      rho(j) = cos(ylatm(2,j))
      if (j == 1) then ! the pole does not have an associated apex height
        ha(j) = fill_value
      elseif (j <= nmlat_h-nhgt_fix) then
        ha(j) = r0km/rho(j)**2 - rekm
      else ! directly take from calculated values to avoid round-off errors
        ha(j) = hgt_fix(k)
      endif
    enddo

! set S2 latitude and rho values
    nmlatS2_h = nmlat_h-1
    allocate(ylatm_s(2,nmlatS2_h))
    allocate(rho_s(nmlatS2_h))
    allocate(ha_s(nmlatS2_h))

    do concurrent (j = 1:nmlatS2_h)
      rho_s(j) = (rho(j)+rho(j+1))/2
      ylatm_s(2,j) = acos(rho_s(j))
      ylatm_s(1,j) = -ylatm_s(2,j)

! calculate apex height of each field line, Richmond 1995 Eq (3.3)
! ha=r0/cos^2lambda-re with lambda modified apex latitude
      ha_s(j) = r0km/rho_s(j)**2 - rekm
    enddo

! set up height levels - R
! these are apex heights of ylatm_s
    nhgt_fix_r = nhgt_fix+1
    allocate(hgt_fix_r(nhgt_fix_r))
    hgt_fix_r(1) = h0km
    do concurrent (k = 2:nhgt_fix_r)
      j = nmlatS2_h-k+2
      hgt_fix_r(k) = ha_s(j)
    enddo

    nmlat_T1 = 2*nmlat_h-1
    nmlat_T2 = 2*nmlatS2_h

! convert from km to m
    do k = 1,nhgt_fix
      hgt_fix(k) = hgt_fix(k)*1000
    enddo
    do k = 1,nhgt_fix_r
      hgt_fix_r(k) = hgt_fix_r(k)*1000
    enddo
    do j = 1,nmlat_h
      ha(j) = ha(j)*1000
    enddo
    do j = 1,nmlatS2_h
      ha_s(j) = ha_s(j)*1000
    enddo

! magnetic longitudes are equidistant
    dlonm = 2*pi/nmlon
    do concurrent (i = 0:nmlon+1)
      ylonm(i) = -pi+(i-1)*dlonm    ! P,S2,R
      ylonm_s(i) = ylonm(i)+dlonm/2 ! S1
    enddo

    deallocate(lat_low)
    deallocate(lat_high)
    deallocate(lat)

  endsubroutine generate_simplified_grid
!-----------------------------------------------------------------------
endmodule grid_module
