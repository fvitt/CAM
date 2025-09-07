module init_module

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  subroutine init_cons
! Set derived constants

    use params_module,only:nmlat_h,nmlon,ylatm
    use cons_module,only:ylatm_JT,jlatm_JT,J3LB
    use solver_module,only:nlonlat

    integer,dimension(1) :: idx

! get the index of ylatm_JT (index is counted from the pole)
! which is the transition latitude between symmetric and asymmetric potentials
    idx = minloc(abs(ylatm(2,:)-ylatm_JT))
    jlatm_JT = idx(1)

! total number of latitudes to be solved (high-lat*2 + low-lat)
! (jlatm_JT-1-2+1)*2 + (nmlat_h-jlatm_JT+1) = nmlat_h+jlatm_JT-3
    nlonlat = (nmlat_h+jlatm_JT-3)*nmlon+1

    allocate(J3LB(2,nmlat_h,0:nmlon+1))
    J3LB = 0

  endsubroutine init_cons
!-----------------------------------------------------------------------
  pure subroutine init_fieldline(npts_p,npts_s1,npts_s2,npts_r, &
    jmax_p,jmax_s1,jmax_s2,jmax_r,size_p,size_s1,size_s2,size_r, &
    qdlat_p,qdlat_s1,qdlat_s2,qdlat_r)

    use params_module,only:nhgt_fix,nhgt_fix_r,nmlat_h,nmlatS2_h, &
      hgt_fix,hgt_fix_r,ha,ha_s,ylatm,ylatm_s
    use cons_module,only:re,r0,fill_value
    use util_module,only:find
    use util_module,only:lamqd_from_apex_coord

    integer,dimension(nmlat_h),intent(out) :: npts_p,npts_s1,npts_r
    integer,dimension(nmlatS2_h),intent(out) :: npts_s2
    integer,dimension(nhgt_fix),intent(out) :: &
      jmax_p,jmax_s1,jmax_s2,size_p,size_s1,size_s2
    integer,dimension(nhgt_fix_r),intent(out) :: jmax_r,size_r
    real(kind=rp),dimension(nhgt_fix,2,nmlat_h),intent(out) :: qdlat_p,qdlat_s1
    real(kind=rp),dimension(nhgt_fix,2,nmlatS2_h),intent(out) :: qdlat_s2
    real(kind=rp),dimension(nhgt_fix_r,2,nmlat_h),intent(out) :: qdlat_r

    integer :: j,isn,k,n

    do concurrent (j = 1:nmlat_h)
      npts_p(j) = find(hgt_fix,ha(j))
      npts_s1(j) = find(hgt_fix,ha(j))
      npts_r(j) = find(hgt_fix_r,ha(j))
    enddo

    do concurrent (j = 1:nmlatS2_h)
      npts_s2(j) = find(hgt_fix,ha_s(j))
    enddo

! the number of points at a fixed height
    do k = 1,nhgt_fix
      do j = 1,nmlat_h ! from open to closed field line
        if (ha(j) < hgt_fix(k)) exit ! find the latitude whose apex height is below this height
      enddo
      n = j-1 ! the apex height of open field lines will be higher than this height
      jmax_p(k) = n
      jmax_s1(k) = n
      if (n == nmlat_h) then
        size_p(k) = jmax_p(k)*2-1
        size_s1(k) = jmax_s1(k)*2-1
      else
        size_p(k) = jmax_p(k)*2
        size_s1(k) = jmax_s1(k)*2
      endif

      do j = 1,nmlatS2_h
        if (ha_s(j) < hgt_fix(k)) exit
      enddo
      n = j-1
      jmax_s2(k) = n
      size_s2(k) = jmax_s2(k)*2
    enddo

    do k = 1,nhgt_fix_r
      do j = 1,nmlat_h
        if (ha(j) < hgt_fix_r(k)) exit
      enddo
      n = j-1
      jmax_r(k) = n
      if (n == nmlat_h) then
        size_r(k) = jmax_r(k)*2-1
      else
        size_r(k) = jmax_r(k)*2
      endif
    enddo

! relationship between P,S1,S2 points for the same index (i,j):
! P(i,j) is S1(i+0.5,j) and S2(i,j+0.5) with j increasing equatorward

    qdlat_p = fill_value
    qdlat_s1 = fill_value
    qdlat_s2 = fill_value
    qdlat_r = fill_value

    do concurrent (j = 1:nmlat_h, isn = 1:2)
      do concurrent (k = 1:npts_p(j))
        qdlat_p(k,isn,j) = lamqd_from_apex_coord(re,r0,ylatm(isn,j),hgt_fix(k))
      enddo

! S1 points are in between P points with respect to mlon, but same mlat
      do concurrent (k = 1:npts_s1(j))
        qdlat_s1(k,isn,j) = lamqd_from_apex_coord(re,r0,ylatm(isn,j),hgt_fix(k))
      enddo

      do concurrent (k = 1:npts_r(j))
        qdlat_r(k,isn,j) = lamqd_from_apex_coord(re,r0,ylatm(isn,j),hgt_fix_r(k))
      enddo
    enddo

! S2 points are in between P points with respect to rho=cos(ylatm), but same mlon
    do concurrent (j = 1:nmlatS2_h, isn = 1:2)
      do concurrent (k = 1:npts_s2(j))
        qdlat_s2(k,isn,j) = lamqd_from_apex_coord(re,r0,ylatm_s(isn,j),hgt_fix(k))
      enddo
    enddo

  endsubroutine init_fieldline
!-----------------------------------------------------------------------
  subroutine get_apex( npts_p,npts_s1,npts_s2,npts_r, &
       qdlat_p,qdlat_s1,qdlat_s2,qdlat_r, &
       glat_p, glon_p, glat_s1, glon_s1, glat_s2, glon_s2, glat_r, glon_r )

    use params_module,only:nhgt_fix,nhgt_fix_r,nmlat_h,nmlatS2_h, &
      hgt_fix,hgt_fix_r,ha,ha_s,ylatm,ylatm_s
    use params_module,only:ylonm,ylonm_s
    use cons_module,only:h0,rtd
    use fieldline_module
    use apex,only: apex_mall,apex_q2g
    use mpi_module,only:mlond0,mlond1,mlatd0,mlatd1

    integer,dimension(nmlat_h),intent(in) :: npts_p,npts_s1,npts_r
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(nhgt_fix,2,nmlat_h),intent(in) :: qdlat_p,qdlat_s1
    real(kind=rp),dimension(nhgt_fix,2,nmlatS2_h),intent(in) :: qdlat_s2
    real(kind=rp),dimension(nhgt_fix_r,2,nmlat_h),intent(in) :: qdlat_r
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1), intent(out) :: &
         glat_p, glon_p, glat_s1, glon_s1, glat_s2, glon_s2
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1), intent(out) :: &
         glat_r, glon_r

    real(kind=rp),parameter :: hr = h0*1e-3_rp
    integer :: i,j,isn,k,icomp, &
      ist ! for generating interpolation grids in apex
    real(kind=rp) :: qdlat,qdlon,alt,gdlat,gdlon, &
! scalar arguments returned by APXMALL
      bmag,si,alon,xlatm,vmp,w,d,be3,sim,xlatqd,f
! non-scalar arguments returned by APXMALL
    real(kind=rp),dimension(3) :: b,bhat,d1,d2,d3,e1,e2,e3,f1,f2,f3,g1,g2,g3

! convert from quasi-dipole to geodetic coordinates
! APXQ2G (input magnetic, output geodetic) is the functional inverse
! of APXALL or APXMALL (input geodetic, output magnetic)

! now the coordinate system is [phi_m,rho=coslam_m,+/-h]
! the base vectors are d' and e'
! so there are new base vectors di' and ei'
! with d1' = d1; d2' = d2; e3' = e3; D; B0;
! what is different is d3'; e1'; e2'
! d3' = -k^/(DsinI)
! e1' = d2'xd3' = (R k^ x grad(rho))/(sqrt(1-3/4rho^2)D sinI)
! e2' = d3'xd1' = (R rho grad(rho) x k^)/(D sinI)
! volume and area vectors:
! Wm'  = |k^ dot grad(phi_m) x grad(rho)|^-1 = R^2 rho/(D |sinI| sqrt(1-3/4rho^2))
! am1' = Wm' grad(phi_m) = R d1'/(D |sinI| sqrt(1-3/4rho^2))
! am2' = Wm' grad(rho_m) = R rho d2'/(D |sinI|)
! am3' = -/+Wm' k^       = R^2 rho d3'/(sqrt(1-3/4rho^2))
! with R = r0 = Re + h0 = Re + hr

! wrap longitude if it is out of bound (-180,180)
! this only happens for halo points in the left and right most mpi task

    do i = mlond0,mlond1
      do j = mlatd0,mlatd1
        if (j>=1 .and. j<=nmlat_h) then
          do isn = 1,2
            do k = 1,npts_p(j)
              qdlat = qdlat_p(k,isn,j)*rtd ! get quasi-dipole latitude
              qdlon = ylonm(i)*rtd ! get quasi-dipole longitude
              if (qdlon < -180) qdlon = qdlon+360
              if (qdlon > 180) qdlon = qdlon-360
              alt = hgt_fix(k)*1.e-3_rp ! convert height from [m] to [km]

              call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ist)
              call apex_mall(gdlat,gdlon,alt,hr,b,bhat,bmag,si, &
                alon,xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, &
                xlatqd,f,f1,f2,f3,g1,g2,g3,ist)

              glat_p(k,isn,j,i) = gdlat
              glon_p(k,isn,j,i) = gdlon

              D_p(k,isn,j,i) = d
              F_p(k,isn,j,i) = f
              vmp_p(k,isn,j,i) = vmp ! magnitude potential Tm (diagnostic for ds calculation)
              bmag_p(k,isn,j,i) = bmag*1e-9_rp ! magnitude of magnetic field, convert from [nT] to [T]
            enddo

            do k = 1,npts_s1(j)
              qdlat = qdlat_s1(k,isn,j)*rtd ! get quasi-dipole latitude
              qdlon = ylonm_s(i)*rtd ! get quasi-dipole longitude
              if (qdlon < -180) qdlon = qdlon+360
              if (qdlon > 180) qdlon = qdlon-360
              alt = hgt_fix(k)/1000 ! convert height from [m] to [km]

              call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ist)
              call apex_mall(gdlat,gdlon,alt,hr,b,bhat,bmag,si, &
                alon,xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, &
                xlatqd,f,f1,f2,f3,g1,g2,g3,ist)

              glat_s1(k,isn,j,i) = gdlat
              glon_s1(k,isn,j,i) = gdlon

! these are the same using the "new" coordinate system and the one from the paper
              D_s1(k,isn,j,i) = d
              F_s1(k,isn,j,i) = f
              be3_s1(k,isn,j,i) = be3*1e-9_rp ! B0=Be3*e3, convert from [nT] to [T]
              d1d1_s1(k,isn,j,i) = dot_product(d1,d1)
              d1d2_s1(k,isn,j,i) = dot_product(d1,d2)
              d2d2_s1(k,isn,j,i) = dot_product(d2,d2)
              do concurrent (icomp = 1:3) ! components (east, north, up) of base vectors
                d1_s1(icomp,k,isn,j,i) = d1(icomp)
                d2_s1(icomp,k,isn,j,i) = d2(icomp)
                d3_s1(icomp,k,isn,j,i) = d3(icomp)
                e1_s1(icomp,k,isn,j,i) = e1(icomp)
                e2_s1(icomp,k,isn,j,i) = e2(icomp)
                e3_s1(icomp,k,isn,j,i) = e3(icomp)
              enddo
            enddo

            do k = 1,npts_r(j)
              qdlat = qdlat_r(k,isn,j)*rtd ! get quasi-dipole latitude
              qdlon = ylonm(i)*rtd ! get quasi-dipole longitude
              if (qdlon < -180) qdlon = qdlon+360
              if (qdlon > 180) qdlon = qdlon-360
              alt = hgt_fix_r(k)/1000 ! convert height from [m] to [km]

              call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ist)
              call apex_mall(gdlat,gdlon,alt,hr,b,bhat,bmag,si, &
                alon,xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, &
                xlatqd,f,f1,f2,f3,g1,g2,g3,ist)

              glat_r(k,isn,j,i) = gdlat
              glon_r(k,isn,j,i) = gdlon

              D_r(k,isn,j,i) = d
              F_r(k,isn,j,i) = f
            enddo
          enddo
        endif

        if (j>=1 .and. j<=nmlatS2_h) then
          do isn = 1,2
            do k = 1,npts_s2(j)
              qdlat = qdlat_s2(k,isn,j)*rtd ! get quasi-dipole latitude
              qdlon = ylonm(i)*rtd ! get quasi-dipole longitude
              if (qdlon < -180) qdlon = qdlon+360
              if (qdlon > 180) qdlon = qdlon-360
              alt = hgt_fix(k)/1000 ! convert height from [m] to [km]

              call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ist)
              call apex_mall(gdlat,gdlon,alt,hr,b,bhat,bmag,si, &
                alon,xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, &
                xlatqd,f,f1,f2,f3,g1,g2,g3,ist)

              glat_s2(k,isn,j,i) = gdlat
              glon_s2(k,isn,j,i) = gdlon

! these are the same using the "new" coordinate system and the one from the paper
              D_s2(k,isn,j,i) = d
              F_s2(k,isn,j,i) = f
              be3_s2(k,isn,j,i) = be3*1e-9_rp ! B0=Be3*e3, convert from [nT] to [T]
              d1d1_s2(k,isn,j,i) = dot_product(d1,d1)
              d1d2_s2(k,isn,j,i) = dot_product(d1,d2)
              d2d2_s2(k,isn,j,i) = dot_product(d2,d2)
              do concurrent (icomp = 1:3) ! components (east, north, up) of base vectors
                d1_s2(icomp,k,isn,j,i) = d1(icomp)
                d2_s2(icomp,k,isn,j,i) = d2(icomp)
                d3_s2(icomp,k,isn,j,i) = d3(icomp)
                e1_s2(icomp,k,isn,j,i) = e1(icomp)
                e2_s2(icomp,k,isn,j,i) = e2(icomp)
                e3_s2(icomp,k,isn,j,i) = e3(icomp)
              enddo
            enddo
          enddo
        endif
      enddo
    enddo

  endsubroutine get_apex
!-----------------------------------------------------------------------
  pure subroutine calculate_a(a1,a3)
! a1,a3 vary from 0 at the magnetic pole to 1 at the magnetic equator

    use params_module,only:nhgt_fix,nhgt_fix_r,nmlat_h,rho_s,hgt_fix_r
    use cons_module,only:pi,re,r0

    real(kind=rp),dimension(nhgt_fix  ,nmlat_h+1),intent(out) :: a1
    real(kind=rp),dimension(nhgt_fix_r,nmlat_h+1),intent(out) :: a3

    integer :: j,k
    real(kind=rp) :: rm1,rp1,dr,tmp,ra,rbar

! in order to include the pole, the first index of a1,a3 represents the location j-0.5, not j+0.5
! however, the first index of m1f,m2f,m3f represents the location j+0.5, which are the S2 points

! set a1,a3 to 1 for equator and beyond
! (points before equator will be overwritten later)
    a1 = 1
    a3 = 1

    do concurrent (k = 1:nhgt_fix)
      a1(k,1) = 0
      a3(k,1) = 0

! normalized radii of the top and bottom of layer k
      rp1 = (hgt_fix_r(k+1)+re)/r0 ! r_k+0.5/R
      rm1 = (hgt_fix_r(k  )+re)/r0 ! r_k-0.5/R
      dr = rp1-rm1

      do concurrent (j = 1:nmlat_h-k)

! first index of a1,a3 is j+1 because this corresponds to position j+0.5
        a3(k,j+1) = 1-sqrt(max(1-rm1*rho_s(j)**2,0.0_rp))

! calculate the normalized radius within the layer, rbar,
! that gives the most accurate value of a1 when a1 is computed by the approximation below
! this calculation of rbar assumes the radius of field lines near the equator
! is parabolic with respect to magnetic latitude
        ra = 1/rho_s(j)**2

! prevent ra from getting too large to affect the numerical accuracy of rbar
! which rapidly asymptotes to 0.5*(rp1+rm1) as ra increases
        if (ra > rp1+16*dr) ra = rp1+16*dr

! Eq (45') page 4c r*/R
        tmp = sqrt(max(ra-rm1,0.0_rp))**3-sqrt(max(ra-rp1,0.0_rp))**3
        rbar = ra-(2*tmp/(3*dr))**2

        tmp = sqrt(max(rbar,0.0_rp))*rho_s(j)
        if (tmp > 1) tmp = 1
        if (tmp < -1) tmp = -1
        a1(k,j+1) = 2*asin(tmp)/pi ! Eq (47')
      enddo
    enddo

! now do a3 for top level
    k = nhgt_fix_r
    a3(k,1) = 0
    rm1 = (hgt_fix_r(k)+re)/r0
    do concurrent (j = 1:nmlat_h-k)
      a3(k,j+1) = 1-sqrt(max(1-rm1*rho_s(j)**2,0.0_rp))
    enddo

  end subroutine calculate_a
!-----------------------------------------------------------------------
  pure subroutine calculate_m( npts_p,npts_s1,npts_s2,npts_r, &
    F_p,F_s1,F_s2,F_r,M3_p,M1_s1,M2_s2,M3_r)
! calculate integrated areas for each surface of a volume (Chapter 5)

    use params_module,only:nhgt_fix,nhgt_fix_r, &
      nmlat_h,nmlatS2_h,ylonm,rho_s,hgt_fix,hgt_fix_r
    use cons_module,only:re,r0,fill_value
    use mpi_module,only:mlond0,mlond1,mlatd0,mlatd1

    integer,dimension(nmlat_h),intent(in) :: npts_p,npts_s1,npts_r
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: F_p,F_s1,F_s2
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: F_r
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: M3_p,M1_s1,M2_s2
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: M3_r

    integer :: i,j,isn,k
    real(kind=rp) :: dlonm,rp1,rm1,rp1s,rm1s,rap1,ram1,tmp

! M1*F, M2*F, M3*F (independent of magnetic longitude)
! these factors later give M1,M2,M3 when divided by F
    real(kind=rp),dimension(nhgt_fix  ,nmlat_h  ) :: m1f
    real(kind=rp),dimension(nhgt_fix  ,nmlatS2_h) :: m2f
    real(kind=rp),dimension(nhgt_fix_r,nmlat_h  ) :: m3f

! assume equidistant longitudinal grid points
    dlonm = ylonm(2)-ylonm(1)

! M1*F (this is different from the note)
    do concurrent (j = 1:nmlat_h, k = 1:nhgt_fix)
      rp1 = (hgt_fix_r(k+1)+re)/r0
      rm1 = (hgt_fix_r(k  )+re)/r0
      rp1s = sqrt(rp1)
      rm1s = sqrt(rm1)

      if (j == 1) then ! let rho(0.5) = rho_s(0) = 0
        rap1 = 1/rho_s(j)**2
        tmp = rp1s*asin(min(rp1s*rho_s(j),1.0_rp)) - rm1s*asin(min(rm1s*rho_s(j),1.0_rp)) + &
              sqrt(max(rap1-rp1,0.0_rp)) - sqrt(max(rap1-rm1,0.0_rp))
      elseif (j == nmlat_h) then ! let rho(J+0.5) = rho_s(J) = 1
        rap1 = 1
        ram1 = 1/rho_s(j-1)**2
        tmp = (rp1s*asin(min(rp1s           ,1.0_rp)) - rm1s*asin(min(rm1s           ,1.0_rp))) - &
              (rp1s*asin(min(rp1s*rho_s(j-1),1.0_rp)) - rm1s*asin(min(rm1s*rho_s(j-1),1.0_rp))) + &
              (sqrt(max(rap1-rp1,0.0_rp)) - sqrt(max(rap1-rm1,0.0_rp))) - &
              (sqrt(max(ram1-rp1,0.0_rp)) - sqrt(max(ram1-rm1,0.0_rp)))
      else
        rap1 = 1/rho_s(j  )**2
        ram1 = 1/rho_s(j-1)**2
        tmp = (rp1s*asin(min(rp1s*rho_s(j  ),1.0_rp)) - rm1s*asin(min(rm1s*rho_s(j  ),1.0_rp))) - &
              (rp1s*asin(min(rp1s*rho_s(j-1),1.0_rp)) - rm1s*asin(min(rm1s*rho_s(j-1),1.0_rp))) + &
              (sqrt(max(rap1-rp1,0.0_rp)) - sqrt(max(rap1-rm1,0.0_rp))) - &
              (sqrt(max(ram1-rp1,0.0_rp)) - sqrt(max(ram1-rm1,0.0_rp)))
      endif
      m1f(k,j) = 2*(hgt_fix(k)+re)**3/r0*max(tmp,0.0_rp)
    enddo

! M2*F
    do concurrent (j = 1:nmlatS2_h, k = 1:nhgt_fix)
      rp1 = (hgt_fix_r(k+1)+re)/r0
      rm1 = (hgt_fix_r(k  )+re)/r0
      rap1 = 1/rho_s(j)**2

      tmp = sqrt(max(rap1-rm1,0.0_rp))-sqrt(max(rap1-rp1,0.0_rp))
      m2f(k,j) = 2*(hgt_fix(k)+re)**3/r0*sqrt(1-3*rho_s(j)**2/4)*dlonm*max(tmp,0.0_rp)
    enddo

! M3*F
    do concurrent (j = 1:nmlat_h, k = 1:nhgt_fix_r)
      rm1 = (hgt_fix_r(k)+re)/r0

      if (j == 1) then ! let rho(0.5) = rho_s(0) = 0
        tmp = 1-sqrt(max(1-rm1*rho_s(j)**2,0.0_rp))
      elseif (j == nmlat_h) then ! let rho(J+0.5) = rho_s(J) = 1
        tmp = sqrt(max(1-rm1*rho_s(j-1)**2,0.0_rp))-sqrt(max(1-rm1,0.0_rp))
      else
        tmp = sqrt(max(1-rm1*rho_s(j-1)**2,0.0_rp))-sqrt(max(1-rm1*rho_s(j)**2,0.0_rp))
      endif
      m3f(k,j) = (hgt_fix_r(k)+re)**2*dlonm*max(tmp,0.0_rp)
    enddo

    M3_p = fill_value
    M1_s1 = fill_value
    M2_s2 = fill_value
    M3_r = fill_value

    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
      do concurrent (k = 1:npts_p(j))
        M3_p(k,isn,j,i) = m3f(k,j)/F_p(k,isn,j,i) ! M3 is used in Jr calculation
      enddo

      do concurrent (k = 1:npts_s1(j))
        M1_s1(k,isn,j,i) = m1f(k,j)/F_s1(k,isn,j,i)
      enddo

      do concurrent (k = 1:npts_r(j))
        M3_r(k,isn,j,i) = m3f(k,j)/F_r(k,isn,j,i)
      enddo
    enddo

    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlatS2_h)
      do concurrent (k = 1:npts_s2(j))
        M2_s2(k,isn,j,i) = m2f(k,j)/F_s2(k,isn,j,i)
      enddo
    enddo

  endsubroutine calculate_m
!-----------------------------------------------------------------------
endmodule init_module
