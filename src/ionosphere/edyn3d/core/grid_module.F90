module grid_module

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  subroutine generate_mag_grid

    use params_module,only:nmlon,nmlat_h,nmlatS2_h,nhgt_fix,nhgt_fix_r, &
      ylonm,ylonm_s,ylatm,ylatm_s,rho,rho_s,ha,ha_s,hgt_fix,hgt_fix_r
    use cons_module,only:re,h0,r0,pi,dtr,fill_value

    real(kind=rp),parameter :: &
      d=10, f=15, h=10, yb=5, yc=18, d1=30, d2=50, d3=55, d4=75, d5=82, hs=6, dhdz=6
    integer :: i,j,k,jns
    real(kind=rp) :: rekm,h0km,r0km,lam1,lam2,lam3,lam4,lam5,hc,h1,rho1,b,c,e,g, &
      y1,y2,y3,y4,y5,ymax,y,lam,rho_loc,ha_loc,dlonm

    rekm = re/1000
    h0km = h0/1000
    r0km = r0/1000

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
        ha(jns) = ha_loc*1e3_rp
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
      ha_s(j) = r0/rho_s(j)**2 - re
    enddo

! set up height levels - R
! these are apex heights of ylatm_s
    hgt_fix_r(1) = h0
    do concurrent (k = 2:nhgt_fix_r)
      j = nmlatS2_h-k+2
      hgt_fix_r(k) = ha_s(j)
    enddo

    dlonm = 2*pi/nmlon
    do concurrent (i = 0:nmlon+1)
! magnetic longitudes - P,S2,R
      ylonm(i) = -pi+(i-1)*dlonm

! magnetic longitudes - S1
      ylonm_s(i) = ylonm(i)+dlonm/2
    enddo

  endsubroutine generate_mag_grid
!-----------------------------------------------------------------------
endmodule grid_module
