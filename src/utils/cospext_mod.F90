!===============================================================================
!  !! This can be after each cospectra is calculated in zonal_fft_mod
!  !! After u_fft*wstar
!  call cospext(cspr,wvlong,wvlxbeg,wvlxend,ntime,mflxxup,mflxxun)
!  !! After v_fft*wstar
!  call cospext(cspr,wvlong,wvlxbeg,wvlxend,ntime,mflxyup,mflxyun)
!
!  !! Calculate vertical divergence to get zonal mean zonal and meridonal forcing
!
!  !! rbar is the zonal mean mass density ( u = unresolved scale , x = zonal direction, y = meridianal)
!
!  fzonal(:,k) = -((mflxxup(:,k-1)+mfluxxun(:,k-1))*rbar(:,k-1)-(mflxxup(:,k+1)+mflxxun(:,k+1))*rbar(:,k+1))/
!                 (pmid(k-1)-pmid(k+1))/rbar(:,k)
!  fmerid(:,k) = -((mflxyup(:,k-1)+mfluxyun(:,k-1))*rbar(:,k-1)-(mflxyup(:,k+1)+mflxyun(:,k+1))*rbar(:,k+1))/
!                 (pmid(k-1)-pmid(k+1))/rbar(:,k)
!
!  !! Check to make sure it's pmid
!
!  !! Then this should be scattered to physics mesh at all longitudes at the correponding latitudes.
!
!  Will need cospectra module to compute and acculate to cospectras
!  - write / read acculated coespectras to / from IC and restart files
!  - compute forcings as described above on "zonal" grid
!  - scatter forcings to physics grid
!
!===============================================================================

module cospext_mod
  use shr_kind_mod, only: r8 => SHR_KIND_R8
  use shr_const_mod,only: pi => SHR_CONST_PI
  use shr_const_mod,only: rearth => SHR_CONST_REARTH ! meters

  implicit none

  real(r8), parameter :: NOTSET = -huge(1._r8)
  integer, parameter :: kxmin = 5

contains

  !=============================================================================
  !! This subroutine should be called once per day. The cospectra should be collected either every hour (if the calculation is
  !! is expensive), or every time step (if the calculation is not so expensive), or somewhere in between.
  !! Will need to save th daily forcing on IC (.i.) file and restart file (.r*), so GW forcing is available for startup and
  !! restart runs, including short (less one day) runs. The daily forcing should also be made optional for the IC file.
  !! The run can still start without the forcing.
  !=============================================================================
  subroutine cospext(nftnum,lat_beg,lat_end,pver,ntime,latrad,cspr,wvlxbeg,wvlxend,press,rhozm,mflxup,mflxun,mflxrp,mflxrn,force_r,force_u,slpp,slpn,ioldest)

    integer, intent(in) :: nftnum,lat_beg,lat_end,pver,ntime
    real(r8), intent(in) :: latrad(lat_beg:lat_end)
    real(r8), intent(in) :: cspr(nftnum, lat_beg:lat_end, pver, ntime)  !! co-spectra resolved by the model

    !! wvlong: longer wavelength side for power index calculation. (2000x10^3 m used in current calculation)
    !! wvlxbeg: longer wavelength of the unresolved range/short end of the resolved range. This is grid size depedent
    !! wvlxend: short wavelength cutoff of the unresolved range -- 2D (lat,lev) field set from a viscous horizontal scale

    real(r8), intent(in) :: wvlxbeg
    real(r8), intent(in) :: wvlxend(lat_beg:lat_end,pver)
    real(r8), intent(in) :: press(pver)
    real(r8), intent(in) :: rhozm(lat_beg:lat_end,pver)

    ! resolved and unresolved forces
    real(r8), optional, intent(out) :: force_r(lat_beg:lat_end,pver), force_u(lat_beg:lat_end,pver)

    !! mflxup: total momentum flux in the positive direction -- resolved and unresolved
    !! mflxun: total momentum flux in the negative direction
    real(r8), intent(out) :: mflxup(lat_beg:lat_end,pver), mflxun(lat_beg:lat_end,pver)
    real(r8), intent(out) :: mflxrp(lat_beg:lat_end,pver), mflxrn(lat_beg:lat_end,pver)

    !! slpp, slpn: spectral slopes for positive and negative branches (optional output)
    real(r8), optional, intent(out) :: slpp(lat_beg:lat_end,pver), slpn(lat_beg:lat_end,pver)

    !! ioldest: time slot of cspr holding the oldest entry when cspr is a ring
    !! buffer over time (default 1, i.e. slots in chronological order). The time
    !! sums are always taken oldest->newest, so results do not depend on where
    !! the ring buffer happens to start (needed for exact restarts).
    integer, optional, intent(in) :: ioldest

    !! kxl: zonal wavenumber corresponding to wvlong
    !! kxm: 2*kxl, kxr: 4*kxl (wavenumbers used to calculate spectral slope (Liu, 2019)
    !! kxbeg: zonal wavenumber corresponding to wvlexbeg
    !! kxend: zonal wavenumber corresponding to wvlexend

    real(r8) :: circlat(lat_beg:lat_end)   !! circumference at a specific latitude
    real(r8) :: slpp_wrk(lat_beg:lat_end,pver), slpn_wrk(lat_beg:lat_end,pver) !! spectral slopes of cospectra for each latitude and level

    integer :: kxl(lat_beg:lat_end), kxm(lat_beg:lat_end), kxr(lat_beg:lat_end), kxbeg(lat_beg:lat_end)
    integer :: kxend(lat_beg:lat_end,pver)
    real(r8), parameter :: krat = 1.6_r8
    integer :: j
    real(r8) :: csprp(nftnum, lat_beg:lat_end, pver),csprn(nftnum, lat_beg:lat_end, pver)
    real(r8) :: csprpi(nftnum, lat_beg:lat_end, pver),csprni(nftnum, lat_beg:lat_end, pver)

    real(r8) :: r2d = 180._r8/pi
    integer :: lat0, lat1
    integer :: t0

    t0 = 1
    if (present(ioldest)) t0 = ioldest

    mflxup = 0._r8
    mflxun = 0._r8

    ! initialize slopes to NOTSET so latitudes outside the calculation band
    ! (|lat| >= 85 deg) have a defined value when exported via optional args
    slpp_wrk = NOTSET
    slpn_wrk = NOTSET

    do j = lat_beg,lat_end
       circlat(j) = 2._r8 * pi * rearth * cos(latrad(j))     !! rearth is Earth radius in m
    enddo

    kxbeg(:)  = nint(circlat(:)/wvlxbeg) ! wvl_r maximum resolved zonal wavenumber
    do j = lat_beg,lat_end
       kxend(j,:) = nint(circlat(j)/wvlxend(j,:)) ! wvl_c zonal wavenumber of cutoff wavelength
    end do

    !! These following 3 lines calculate the scale invariance range according to the short wavelength of the resolved range
    kxr(:) = kxbeg(:)
    kxm(:) = nint(kxr(:)/krat)
    kxl(:) = nint(kxm(:)/krat)

    lat0 = 0
    lat1 = -1

    !! calculations are done equatorward of 85 deg latitude
    do j = lat_beg,lat_end
       if ( abs(latrad(j)*r2d) < 85._r8 ) then
          if (lat0<1) lat0 = j
          lat1 = j
       end if
    end do

    !! separate the cospectra into positive and negative branches from the accumulated spectra
    call spectral_separate(cspr,csprp,csprn,csprpi,csprni)
    call spectral_slope( csprpi,kxl,slpp_wrk)
    call spectral_slope(-csprni,kxl,slpn_wrk)
    call momentum_fluxes(kxl,kxbeg,kxend,csprp,csprn,slpp_wrk,slpn_wrk,mflxup,mflxun, mflxrp,mflxrn )

    if (present(force_r) .and. present(force_u)) then
       ! Calculate the vertical divergence of the resolved and unresolved fluxes (tendencies or forces)
       do j = lat_beg,lat_end
          force_r(j,:) = -vertdiv( press(:), rhozm(j,:)*(mflxrp(j,:) + mflxrn(j,:)) )/rhozm(j,:)
          force_u(j,:) = -vertdiv( press(:), rhozm(j,:)*(mflxup(j,:) + mflxun(j,:)) )/rhozm(j,:)
       end do
    end if

    if (present(slpp)) slpp = slpp_wrk
    if (present(slpn)) slpn = slpn_wrk

  contains

    !===========================================================================
    !===========================================================================
    subroutine spectral_separate(cspr,csprp,csprn,csprpi,csprni)

      !! co-spectra resolved by the model
      real(r8), intent(in) :: cspr(nftnum, lat_beg:lat_end, pver, ntime)

      !! csprp: cospectra that is positive and averaged over the accumulated time period
      !! csprn: cospectra that is negative and averaged over the accumulated time period
      real(r8), intent(out) :: csprp(nftnum, lat_beg:lat_end, pver)
      real(r8), intent(out) :: csprn(nftnum, lat_beg:lat_end, pver)

      ! interpolated over the zeros
      real(r8), intent(out) :: csprpi(nftnum, lat_beg:lat_end, pver)
      real(r8), intent(out) :: csprni(nftnum, lat_beg:lat_end, pver)

      !! running sums and counts of the positive / negative components over time
      real(r8) :: sump(nftnum), sumn(nftnum)
      integer  :: npos(nftnum), nneg(nftnum)
      real(r8) :: val
      integer :: i,j,k,t,tt

      csprp = 0._r8
      csprn = 0._r8
      csprpi = 0._r8
      csprni = 0._r8

      !! Find the time-average positive and negative cospectral components.
      !! Single pass over cspr with unit-stride inner loop; no 4-D temporaries.
      !! Time slices are summed in chronological order starting at slot t0
      !! (t0..ntime, then 1..t0-1) so the result is independent of the
      !! ring-buffer offset.

      do k = 1, pver
         do j = lat0,lat1

            sump = 0._r8
            sumn = 0._r8
            npos = 0
            nneg = 0
            do tt = 0, ntime-1
               t = t0 + tt
               if (t > ntime) t = t - ntime
               do i = 1, nftnum
                  val = cspr(i,j,k,t)
                  if (val > 0._r8) then
                     sump(i) = sump(i) + val
                     npos(i) = npos(i) + 1
                  else if (val < 0._r8) then
                     sumn(i) = sumn(i) + val
                     nneg(i) = nneg(i) + 1
                  end if
               end do
            end do

            do i = 1, nftnum
               if (npos(i) /= 0) csprp(i,j,k) = sump(i) / real(npos(i),kind=r8)
               if (nneg(i) /= 0) csprn(i,j,k) = sumn(i) / real(nneg(i),kind=r8)
            end do

            ! interpolate to replace the zeros -- the interpolated spectra are
            ! only used by spectral_slope, which skips kxbeg <= kxmin
            if (kxbeg(j) > kxmin) then
               csprpi(:,j,k) = interp_cospec(csprp(:,j,k))
               csprni(:,j,k) = interp_cospec(csprn(:,j,k))
            end if

         end do
      end do

    end subroutine spectral_separate

    !===========================================================================
    subroutine spectral_slope(spct,kxl,slp)

      real(r8), intent(in) :: spct(nftnum,lat_beg:lat_end,pver)
      integer, intent(in) :: kxl(lat_beg:lat_end)
      real(r8), intent(out) :: slp(lat_beg:lat_end,pver)  ! spectral slope of power law over kxl to 4*kxl

      real(r8) :: silm, simr  ! integration of spct over kxl to 2*kxl, and over 2*kxl to 4*kxl
      integer :: j,k
      real(r8), parameter :: nominal_slope = 7._r8/6._r8

      do k = 1,pver
         do j=lat0,lat1
            if (kxbeg(j) > kxmin) then
!               silm = sum(spct(kxl(j):2*kxl(j),j,k))
!               simr = sum(spct(2*kxl(j):4*kxl(j),j,k))
               silm = sum(spct(kxl(j):kxm(j),j,k))
               simr = sum(spct(kxm(j):kxr(j),j,k))
               if (silm>0.0_r8 .and. simr>0.0_r8) then
!                  slp(j,k) = MAX(1._r8-log(simr/silm)/log(2._r8),-0.3_r8)
!                  slp(j,k) = MAX(1._r8-log(simr/silm)/log(krat),-0.3_r8)
                  slp(j,k) = MAX(1._r8-log(simr/silm)/log(krat),0.667_r8)
               else
                  slp(j,k) = NOTSET
               end if
            else
               slp(j,k) = nominal_slope
            end if
         enddo
      enddo

    end subroutine spectral_slope

    !===========================================================================
    subroutine momentum_fluxes(kxl,kxbeg,kxend, csprp,csprn,slpp,slpn, mflxup,mflxun, mflxrp,mflxrn )

      real(r8), intent(in)  :: csprp(nftnum,lat_beg:lat_end,pver)
      real(r8), intent(in)  :: csprn(nftnum,lat_beg:lat_end,pver)
      integer,  intent(in)  :: kxl(lat_beg:lat_end)
      integer,  intent(in)  :: kxbeg(lat_beg:lat_end)
      integer,  intent(in)  :: kxend(lat_beg:lat_end,pver)
      real(r8), intent(in)  :: slpp(lat_beg:lat_end,pver)
      real(r8), intent(in)  :: slpn(lat_beg:lat_end,pver)

      real(r8), intent(out) :: mflxup(lat_beg:lat_end,pver)
      real(r8), intent(out) :: mflxun(lat_beg:lat_end,pver)
      real(r8), intent(out) :: mflxrp(lat_beg:lat_end,pver)
      real(r8), intent(out) :: mflxrn(lat_beg:lat_end,pver)

      real(r8) :: bp, bn, fp, fn
      integer :: j,k

      real(r8), parameter :: kxlc = 10._r8
      mflxrp = 0._r8
      mflxrn = 0._r8
      mflxup = 0._r8
      mflxun = 0._r8

      do k=1,pver
         do j=lat0,lat1
            if (kxbeg(j) > kxmin) then
               if (slpp(j,k)/=NOTSET.and.slpp(j,k)/=1._r8) then
                  mflxrp(j,k) = sum(csprp(kxl(j):kxbeg(j),j,k),1)
                  bp = 1._r8-slpp(j,k)
                  fp = (real(kxend(j,k),r8)**bp-real(kxbeg(j),r8)**bp)/(real(kxbeg(j),r8)**bp-kxl(j)**bp)
                  mflxup(j,k) = mflxrp(j,k)*fp
!                  mflxrp(j,k) = sum(csprp(nint(kxlc):kxbeg(j),j,k),1)
!                  bp = 1._r8-slpp(j,k)
!                  fp = (real(kxend(j,k),r8)**bp-real(kxbeg(j),r8)**bp)/(real(kxbeg(j),r8)**bp-kxlc**bp)
!                  mflxup(j,k) = mflxrp(j,k)*fp
               end if
               if (slpp(j,k)/=NOTSET.and.slpp(j,k)==1._r8) then
                  mflxrp(j,k) = sum(csprp(kxl(j):kxbeg(j),j,k),1)
                  fp = log(real(kxend(j,k),r8)/real(kxbeg(j),r8))/log(real(kxbeg(j),r8)/real(kxl(j),r8))
                  mflxup(j,k) = mflxrp(j,k)*fp
!                  mflxrp(j,k) = sum(csprp(nint(kxlc):kxbeg(j),j,k),1)
!                  fp = log(real(kxend(j,k),r8)/real(kxbeg(j),r8))/log(real(kxbeg(j),r8)/kxlc)
!                  mflxup(j,k) = mflxrp(j,k)*fp
               end if
               if (slpp(j,k)==NOTSET) then
                  mflxrp(j,k) = 0._r8
                  mflxup(j,k) = 0._r8
               end if

               if (slpn(j,k)/=NOTSET.and.slpn(j,k)/=1._r8) then
                  mflxrn(j,k) = sum(csprn(kxl(j):kxbeg(j),j,k),1)
                  bn = 1._r8-slpn(j,k)
                  fn = (real(kxend(j,k),r8)**bn-real(kxbeg(j),r8)**bn)/(real(kxbeg(j),r8)**bn-kxl(j)**bn)
                  mflxun(j,k) = mflxrn(j,k)*fn
!                  mflxrn(j,k) = sum(csprn(nint(kxlc):kxbeg(j),j,k),1)
!                  bn = 1._r8-slpn(j,k)
!                  fn = (real(kxend(j,k),r8)**bn-real(kxbeg(j),r8)**bn)/(real(kxbeg(j),r8)**bn-kxlc**bn)
!                  mflxun(j,k) = mflxrn(j,k)*fn
               end if
               if (slpn(j,k)/=NOTSET.and.slpn(j,k)==1._r8) then
                  mflxrn(j,k) = sum(csprn(kxl(j):kxbeg(j),j,k),1)
                  fn = log(real(kxend(j,k),r8)/real(kxbeg(j),r8))/log(real(kxbeg(j),r8)/real(kxl(j),r8))
                  mflxun(j,k) = mflxrn(j,k)*fn
!                  mflxrn(j,k) = sum(csprn(nint(kxlc):kxbeg(j),j,k),1)
!                  fn = log(real(kxend(j,k),r8)/real(kxbeg(j),r8))/log(real(kxbeg(j),r8)/kxlc)
!                  mflxun(j,k) = mflxrn(j,k)*fn
               end if
               if (slpn(j,k)==NOTSET) then
                  mflxrn(j,k) = 0._r8
                  mflxun(j,k) = 0._r8
               end if
            endif
         enddo
      enddo

    end subroutine momentum_fluxes

    !=============================================================================
    ! interpolate to replace zeros in the cospectra
    ! Zeros are filled by linear interpolation in log(|cosp|) vs. wavenumber
    ! index between the nearest nonzero neighbors; leading/trailing zeros take
    ! the nearest nonzero value. If fewer than two nonzero values exist the
    ! result is all zeros. Done in a single O(nftnum) sweep (equivalent to the
    ! former lininterp-based version, which had an O(nzeros*nonzeros) search
    ! plus heap allocations on every call).
    function interp_cospec(cosp) result(cospi)

      real(r8), intent(in) :: cosp(nftnum)
      real(r8) :: cospi(nftnum)

      integer :: i, ii, ilo, ihi, n_nonzero
      real(r8) :: sgn, loglo, loghi, dx, wgts, wgtn

      cospi = cosp
      n_nonzero = count(cosp/=0._r8)

      if (n_nonzero == nftnum) return   ! no zeros to fill

      if (n_nonzero < 2) then
         cospi = 0._r8
         return
      end if

      if (any(cosp<0._r8)) then
         sgn = -1._r8
      else
         sgn = 1._r8
      end if

      ilo = 0   ! index of the most recent nonzero value (0 = none yet)
      i = 1
      do while (i <= nftnum)
         if (cosp(i) /= 0._r8) then
            ilo = i
            i = i + 1
            cycle
         end if

         ! cosp(i) starts a run of zeros; find the next nonzero value
         ihi = i + 1
         do while (ihi <= nftnum)
            if (cosp(ihi) /= 0._r8) exit
            ihi = ihi + 1
         end do

         if (ilo == 0) then
            ! leading zeros: copy the first nonzero value
            loghi = log(sgn*cosp(ihi))
            cospi(i:ihi-1) = sgn*exp(loghi)
         else if (ihi > nftnum) then
            ! trailing zeros: copy the last nonzero value
            loglo = log(sgn*cosp(ilo))
            cospi(i:nftnum) = sgn*exp(loglo)
         else
            loglo = log(sgn*cosp(ilo))
            loghi = log(sgn*cosp(ihi))
            dx = real(ihi,r8) - real(ilo,r8)
            do ii = i, ihi-1
               wgts = (real(ihi,r8) - real(ii,r8))/dx
               wgtn = (real(ii,r8) - real(ilo,r8))/dx
               cospi(ii) = sgn*exp(loglo*wgts + loghi*wgtn)
            end do
         end if

         i = ihi
      end do

    end function interp_cospec

    !===========================================================================
    ! vertical divergence
    function vertdiv( z, y ) result(dydz)

      real(r8),intent(in) :: z(pver)
      real(r8),intent(in) :: y(pver)
      real(r8) :: dydz(pver)

      integer :: k

      ! compute vertical derivative dy/dz (vertical divergence)
      do k = 2,pver-1
         dydz(k) = (y(k+1)-y(k-1)) / (z(k+1)-z(k-1))
      end do

      ! boundary points: use one-sided differences
      dydz(   1) = (y(2) - y(1)) / (z(2) - z(1))
      dydz(pver) = (y(pver) - y(pver-1)) / (z(pver) - z(pver-1))

    end function vertdiv

  end subroutine cospext

end module cospext_mod
