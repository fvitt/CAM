!  !!Psudo code that calls cospext. This can be after each cospectra is calculated in zonal_fft_mod
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
!  ...

module cospext_mod
  use shr_kind_mod, only: r8 => SHR_KIND_R8
  use shr_const_mod,only: pi => SHR_CONST_PI
  use shr_const_mod,only: rearth => SHR_CONST_REARTH ! meters

  implicit none

  real(r8), parameter :: NOTSET = -huge(1._r8)

contains

  !! This subroutine should be called once per day. The cospectra should be collected either every hour (if the calculation is
  !! is expensive), or every time step (if the calculation is not so expensive), or somewhere in between.
  !! Will need to save th daily forcing on IC (.i.) file and restart file (.r*), so GW forcing is available for startup and
  !! restart runs, including short (less one day) runs. The daily forcing should also be made optional for the IC file.
  !! The run can still start without the forcing.

  subroutine cospext(nftnum,lat_beg,lat_end, pver,ntime, latrad, cspr, wvlxbeg,wvlxend, mflxup,mflxun)

    integer, intent(in) :: nftnum,lat_beg,lat_end,pver,ntime
    real(r8), intent(in) :: latrad(lat_beg:lat_end)
    real(r8), intent(in) :: cspr(nftnum, lat_beg:lat_end, pver, ntime)  !! co-spectra resolved by the model

    !! wvlong: longer wavelength side for power index calculation. (2000x10^3 m used in current calculation)
    !! wvlxbeg: longer wavelength of the unresolved range. This is grid size depedent
    !! wvlxend: short wavelength cutoff of the unresolved range (20x10^2 m assumed in current calculation)

    real(r8), intent(in) :: wvlxbeg,wvlxend

    !! mflxup: total momentum flux in the positive direction over unresolved scales
    !! mflxun: total momentum flux in the negative direction over unresolved scales
    real(r8), intent(out) :: mflxup(lat_beg:lat_end,pver), mflxun(lat_beg:lat_end,pver)

    !! kxl: zonal wavenumber corresponding to wvlong
    !! kxm: 2*kxl, kxr: 4*kxl (wavenumbers used to calculate spectral slope (Liu, 2019)
    !! kxbeg: zonal wavenumber corresponding to wvlexbeg
    !! kxend: zonal wavenumber corresponding to wvlexend

    real(r8) :: circlat(lat_beg:lat_end)   !! circumference at a specific latitude
    real(r8) :: slpp(lat_beg:lat_end,pver), slpn(lat_beg:lat_end,pver) !! spectral slopes of cospectra for each latitude and level

    integer :: kxl(lat_beg:lat_end), kxm(lat_beg:lat_end), kxr(lat_beg:lat_end), kxbeg(lat_beg:lat_end), kxend(lat_beg:lat_end)
    integer :: j
    real(r8) :: csprp(nftnum, lat_beg:lat_end, pver),csprn(nftnum, lat_beg:lat_end, pver)

    do j = lat_beg,lat_end
       circlat(j) = 2._r8 * pi * rearth * cos(latrad(j))     !! rearth is Earth radius in m
    enddo

    kxbeg(:)  = nint(circlat(:)/wvlxbeg)
    kxend(:)  = nint(circlat(:)/wvlxend)

!! These following 3 lines calculate the scale invariance range according to the short wavelength of the resolved range
    kxr(:) = kxbeg(:)
    kxm(:) = kxr(:)/2
    kxl(:) = kxm(:)/2

    !! The following calculation should be done equatorward of 85 deg latitude. A conditional should be set in the subroutine

    call spectral_separate(nftnum,lat_beg,lat_end,pver,ntime, cspr,csprp,csprn)  !! separate the cospectra into positive and negative branches from the accumulated spectra
    call spectral_slope(nftnum,lat_beg,lat_end,pver,csprp,kxl,slpp)
    call spectral_slope(nftnum,lat_beg,lat_end,pver,-csprn,kxl,slpn)
    call scale_unres(nftnum,lat_beg,lat_end,pver,kxl,kxbeg,kxend,csprp,csprn,slpp,slpn,mflxup,mflxun)

    return
  end subroutine cospext


  subroutine spectral_separate(nftnum,lat_beg,lat_end,pver,ntime,cspr,csprp,csprn)

    integer, intent(in) :: nftnum,lat_beg,lat_end,pver,ntime
    real(r8), intent(in) :: cspr(nftnum, lat_beg:lat_end, pver, ntime)  !! co-spectra resolved by the model

    !! csprp: cospectra that is positive and averaged over the accumulated time period
    !! csprn: cospectra that is negative and averaged over the accumulated time period
    real(r8), intent(out) :: csprp(nftnum, lat_beg:lat_end, pver),csprn(nftnum, lat_beg:lat_end, pver)

    !! csprps: cospectra that is positive at each time slice
    !! csprns: cospectra that is negative at each time slice
    real(r8) :: csprps(nftnum,lat_beg:lat_end,pver,ntime),csprns(nftnum,lat_beg:lat_end,pver,ntime)
    real(r8) :: wrk1(ntime)
    integer :: i,j,k
    integer :: npos, nneg

    csprps = 0._r8
    csprns = 0._r8
    csprp = 0._r8
    csprn = 0._r8

    !! separate out the positve and negative spectral component at each time step

    where(cspr>0._r8)
       csprps = cspr
    end where
    where(cspr<0._r8)
       csprns = cspr
    end where

    !! Find the average positive and negative cospectral components

    do k = 1, pver
       do j = lat_beg, lat_end
          do i = 1, nftnum

             wrk1 = csprps(i, j, k, :)
             npos = count(wrk1/=0._r8)
             if (npos /= 0) then
                csprp(i, j, k) = sum(wrk1) / real(npos,kind=r8)
             end if

             wrk1 = csprns(i, j, k, :)
             nneg = count(wrk1/=0._r8)
             if (nneg /= 0) then
                csprn(i, j, k) = sum(wrk1) / real(nneg,kind=r8)
             end if

          end do
       end do
    end do
    return
  end subroutine spectral_separate

  subroutine spectral_slope(nftnum,lat_beg,lat_end,pver,spct,kxl,slp)
    integer, intent(in) :: nftnum,lat_beg,lat_end,pver
    real(r8), intent(in) :: spct(nftnum,lat_beg:lat_end,pver)
    integer, intent(in) :: kxl(lat_beg:lat_end)
    real(r8), intent(out) :: slp(lat_beg:lat_end,pver)  ! spectral slope of power law over kxl to 4*kxl

    real(r8) :: silm, simr  ! integration of spct over kxl to 2*kxl, and over 2*kxl to 4*kxl
    integer :: j,k

    do k = 1,pver
       do j=lat_beg,lat_end
          if (kxl(j)>0) then
             silm = sum(spct(kxl(j):2*kxl(j),j,k))
             simr = sum(spct(2*kxl(j):4*kxl(j),j,k))
             if (silm>0.0_r8 .and. simr>0.0_r8) then
                slp(j,k) = MAX(1._r8-log(simr/silm)/log(2._r8),-1._r8)
             else
                slp(j,k) = NOTSET
             end if
          end if
       enddo
    enddo

    return
  end subroutine spectral_slope

  subroutine scale_unres(nftnum,lat_beg,lat_end,pver,kxl,kxbeg,kxend,csprp,csprn,slpp,slpn,mflxup,mflxun)

    integer,  intent(in)  :: nftnum,lat_beg,lat_end,pver
    real(r8), intent(in)  :: csprp(nftnum,lat_beg:lat_end,pver)
    real(r8), intent(in)  :: csprn(nftnum,lat_beg:lat_end,pver)
    integer,  intent(in)  :: kxl(lat_beg:lat_end)
    integer,  intent(in)  :: kxbeg(lat_beg:lat_end)
    integer,  intent(in)  :: kxend(lat_beg:lat_end)
    real(r8), intent(in)  :: slpp(lat_beg:lat_end,pver)
    real(r8), intent(in)  :: slpn(lat_beg:lat_end,pver)

    real(r8), intent(out) :: mflxup(lat_beg:lat_end,pver)
    real(r8), intent(out) :: mflxun(lat_beg:lat_end,pver)

    real(r8) :: siresp, siresn, bp, bn, fp, fn
    integer :: j,k

    mflxup = 0._r8
    mflxun = 0._r8

    do k=1,pver
       do j=lat_beg,lat_end
          if (kxl(j) > 2) then
             if (slpp(j,k)/=NOTSET.and.slpp(j,k)/=1._r8) then
                siresp = sum(csprp(kxl(j):kxbeg(j),j,k),1)
                bp = 1._r8-slpp(j,k)
                fp = (real(kxend(j),r8)**bp-real(kxbeg(j),r8)**bp)/(real(kxbeg(j),r8)**bp-real(kxl(j),r8)**bp)
                mflxup(j,k) = siresp*fp
             end if
             if (slpp(j,k)/=NOTSET.and.slpp(j,k)==1._r8) then
                siresp = sum(csprp(kxl(j):kxbeg(j),j,k),1)
                fp = log(real(kxend(j),r8)/real(kxbeg(j),r8))/log(real(kxbeg(j),r8)/real(kxl(j),r8))
                mflxup(j,k) = siresp*fp
             end if
             if (slpp(j,k)==NOTSET) then
                mflxup(j,k) = 0._r8
             end if
                
             if (slpn(j,k)/=NOTSET.and.slpn(j,k)/=1._r8) then
                siresn = sum(csprn(kxl(j):kxbeg(j),j,k),1)
                bn = 1._r8-slpn(j,k)
                fn = (real(kxend(j),r8)**bn-real(kxbeg(j),r8)**bn)/(real(kxbeg(j),r8)**bn-real(kxl(j),r8)**bn)
                mflxun(j,k) = siresn*fn
             end if
             if (slpn(j,k)/=NOTSET.and.slpn(j,k)==1._r8) then
                siresp = sum(csprn(kxl(j):kxbeg(j),j,k),1)
                fn = log(real(kxend(j),r8)/real(kxbeg(j),r8))/log(real(kxbeg(j),r8)/real(kxl(j),r8))
                mflxun(j,k) = siresn*fn
             end if
             if (slpn(j,k)==NOTSET) then
                mflxun(j,k) = 0._r8
             end if
             
          endif
       enddo
    enddo

    return
  end subroutine scale_unres


end module cospext_mod
