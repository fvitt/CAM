module major_mod

  use iso_fortran_env, only: rp=>real64

  implicit none

  private
  public :: comp

! exponent factor for diff_fac
  real(kind=rp),dimension(3),parameter :: ss = (/1.710_rp,1.749_rp,1.718_rp/)

! mutual thermal diffusion coefficients among major species
  real(kind=rp),dimension(3,4),parameter :: &
    psi = reshape( &
      (/0.0_rp ,0.673_rp,0.270_rp, &
        1.35_rp,0.0_rp  ,0.404_rp, &
        2.16_rp,1.616_rp,0.0_rp  , &
        1.11_rp,0.769_rp,0.322_rp/),(/3,4/))

  real(kind=rp),parameter :: &
    rmass_o2 = 32, rmass_o1 = 16, rmass_he = 4 , rmass_n2  = 28, &
    rmass_n1 = 14, rmass_no = 30, rmass_ar = 40, & ! molar mass
    grav = 870  ! acceleration due to gravity (dependent on lower boundary)

  real(kind=rp),parameter :: p0 = 5e-4_rp         ! standard pressure

contains
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  pure subroutine comp(nlevp1,dz,expzm,expzmid, step,dfactor,tlbc,bo2,bo1,bhe,he_ubc, &
    difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
    o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm, &
    o2_nm_hd,o1_nm_hd,he_nm_hd, &
    prod,loss,o2_upd,o1_upd,he_upd)
! advance major species O2, O, He and N2

    use lbc_mod,only:fb,b
    use matutil_mod,only:matinv3

    integer,intent(in) :: nlevp1

    real(kind=rp),intent(in) :: step,dfactor,tlbc,bo2,bo1,bhe,he_ubc
    real(kind=rp),dimension(nlevp1),intent(in) :: &
      difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
      o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm, &
      o2_nm_hd,o1_nm_hd,he_nm_hd
    real(rp),intent(in) :: dz(nlevp1)
    real(rp),intent(in) :: expzm(nlevp1)
    real(rp),intent(in) :: expzmid
    real(kind=rp),dimension(3,nlevp1),intent(in) :: prod
    real(kind=rp),dimension(3,3,nlevp1),intent(in) :: loss
    real(kind=rp),dimension(nlevp1),intent(out) :: o2_upd,o1_upd,he_upd

    real(kind=rp),parameter :: tau = 1.86e3_rp, t00 = 273, &
      thdiffalpha = -0.38_rp ! thermal diffusion coefficient (alpha) for Helium
    integer,dimension(3,3),parameter :: delta = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
    integer :: k,m,n
    real(kind=rp) :: &
      bn2,bmbar, &      ! at midpoint level 0 (not interface level 1)
      flx00,o1_ub,he_ub ! Helium Mass Flux at upper boundary
    real(kind=rp),dimension(3) :: epep, &
      diff_fac ! correction factor for diffusion coefficients between He and O2, O, N2
    real(kind=rp),dimension(3,3) :: invalpha
    real(kind=rp),dimension(nlevp1) :: dtdz,dmdz,wks1, &
      eddyp,eddyq,eddyr,eddyppart,eddyrpart,eddyp1part,eddyr1part
    real(kind=rp),dimension(3,nlevp1) :: &
      ep,fk,upd,dpdt,eddydif,veradv,loss_out,moldif
    real(kind=rp),dimension(3,3,nlevp1) :: &
      alpha,molp,molq,molr,molp1,molr1,pk,qk,rk

! N2, mbar at midpoint level 0 (not interface level 1)
    bn2 = max(1-bo2-bo1-bhe,0.0_rp)
    bmbar = 1/(bo2/rmass_o2+ &
               bo1/rmass_o1+ &
               bhe/rmass_he+ &
               bn2/rmass_n2)

    dtdz(1) = (tn(1)-tlbc)*2/dz(1)
    dmdz(1) = (mbar(1)-bmbar)/dz(1)
    do k = 2,nlevp1
      dtdz(k) = (tn(k)-tn(k-1))/dz(k)
      dmdz(k) = (mbar(k)-mbar(k-1))/dz(k)
    enddo

! WKS1 = MBAR/M4*(T00/T)**0.25/TAU
    wks1 = barm*(t00/tni)**0.25_rp/(tau*rmass_n2)

! EP = 1-(M+DMBAR/DZ)/MBAR
    ep(1,:) = 1-(rmass_o2+dmdz)/barm
    ep(2,:) = 1-(rmass_o1+dmdz)/barm
    ep(3,:) = 1-(rmass_he+dmdz)/barm-thdiffalpha*dtdz/tni

    do k = 1,nlevp1

! correction factors for mutual diffusion between He and O2, O, N2
      do n = 1,3
        diff_fac(n) = (tni(k)/t00)**(1.75_rp-ss(n))
      enddo

! alpha matrix
      alpha(1,1,k) = -psi(1,4)- &
        (psi(1,2)-psi(1,4))*o1i(k)- &
        (diff_fac(1)*psi(1,3)-psi(1,4))*hei(k)
      alpha(2,2,k) = -psi(2,4)- &
        (psi(2,1)-psi(2,4))*o2i(k)- &
        (diff_fac(2)*psi(2,3)-psi(2,4))*hei(k)
      alpha(3,3,k) = -diff_fac(3)*psi(3,4)- &
        (diff_fac(1)*psi(3,1)-diff_fac(3)*psi(3,4))*o2i(k)- &
        (diff_fac(2)*psi(3,2)-diff_fac(3)*psi(3,4))*o1i(k)
      alpha(1,2,k) = (psi(1,2)-psi(1,4))*o2i(k)
      alpha(1,3,k) = (diff_fac(1)*psi(1,3)-psi(1,4))*o2i(k)
      alpha(2,1,k) = (psi(2,1)-psi(2,4))*o1i(k)
      alpha(2,3,k) = (diff_fac(2)*psi(2,3)-psi(2,4))*o1i(k)
      alpha(3,1,k) = (diff_fac(1)*psi(3,1)-diff_fac(3)*psi(3,4))*hei(k)
      alpha(3,2,k) = (diff_fac(2)*psi(3,2)-diff_fac(3)*psi(3,4))*hei(k)

! molecular diffusion coefficients of O2, O, He
      invalpha = matinv3(alpha(:,:,k))
      do n = 1,3
        do m = 1,3
          molp (m,n,k) = invalpha(m,n)*wks1(k)*(1/dz(k)+ep(n,k)/2)
          molr1(m,n,k) = invalpha(m,n)*wks1(k)*(1/dz(k)-ep(n,k)/2)
        enddo
      enddo
    enddo

! eddy diffusion coefficients (part)
    eddyppart  = difk*(1/dz-dmdz/(barm*2))
    eddyr1part = difk*(1/dz+dmdz/(barm*2))

    do k = 1,nlevp1-1
      molp1(:,:,k) = molp (:,:,k+1)
      molr (:,:,k) = molr1(:,:,k+1)
      eddyp1part(k) = eddyppart (k+1)
      eddyrpart (k) = eddyr1part(k+1)
    enddo
    molp1(:,:,nlevp1) = 2*molp (:,:,nlevp1)-molp (:,:,nlevp1-1)
    molr (:,:,nlevp1) = 2*molr1(:,:,nlevp1)-molr1(:,:,nlevp1-1)
    eddyp1part(nlevp1) =     2*eddyppart (nlevp1)-eddyppart (nlevp1-1)
    eddyrpart (nlevp1) = max(2*eddyr1part(nlevp1)-eddyr1part(nlevp1-1),0.0_rp)

    molq = molp1+molr1

! finish the remaining part of eddy diffusion coefficients
    eddyp = dfactor*eddyppart/expzmid
    eddyr = dfactor*eddyrpart*expzmid
    eddyq = dfactor*(eddyp1part*expzmid+eddyr1part/expzmid)

    do n = 1,3
      do m = 1,3
        pk(m,n,:) = (molp(m,n,:)-expzm*delta(m,n)*(eddyp+wmid/2))/dz(k)
        rk(m,n,:) = (molr(m,n,:)-expzm*delta(m,n)*(eddyr-wmid/2))/dz(k)
        qk(m,n,:) = -molq(m,n,:)/dz(k)+ &
          expzm*(delta(m,n)*(eddyq/dz(k)+1/(2*step))-loss(m,n,:))
      enddo
    enddo

! add explicit source terms to fk
    fk(1,:) = expzm*(prod(1,:)+o2_nm/(2*step)+o2_nm_hd-o2_hadv)
    fk(2,:) = expzm*(prod(2,:)+o1_nm/(2*step)+o1_nm_hd-o1_hadv)
    fk(3,:) = expzm*(prod(3,:)+he_nm/(2*step)+he_nm_hd-he_hadv)

! lower boundaries
    qk(:,:,1) = qk(:,:,1)+matmul(pk(:,:,1),b)
    do n = 1,3
      fk(n,1) = fk(n,1)-dot_product(pk(n,:,1),fb)
    enddo
    pk(:,:,1) = 0

! upper boundary
    epep = (2+ep(:,nlevp1)*dz(nlevp1))/(2-ep(:,nlevp1)*dz(nlevp1))
    do n = 1,3
      do m = 1,3
        qk(m,n,nlevp1-1) = qk(m,n,nlevp1-1)+epep(n)*rk(m,n,nlevp1-1)
      enddo
    enddo

! Eric Sutton: calculate Helium lateral exospheric transport mass flux at upper boundary
    flx00 = wks1(nlevp1)*p0/grav
    o1_ub = he_ubc*(alpha(2,3,nlevp1)-alpha(2,2,nlevp1))/(flx00*(1/dz(nlevp1)-ep(2,nlevp1)/2))
    he_ub = he_ubc*(alpha(3,3,nlevp1)-alpha(3,2,nlevp1))/(flx00*(1/dz(nlevp1)-ep(3,nlevp1)/2))
    fk(:,nlevp1-1) = fk(:,nlevp1-1)-rk(:,2,nlevp1-1)*o1_ub-rk(:,3,nlevp1-1)*he_ub
    rk(:,:,nlevp1-1) = 0

    upd = blktri(pk,qk,rk,fk,nlevp1)

! upper boundaries
    upd(:,nlevp1) = epep*upd(:,nlevp1-1)
    upd(2,nlevp1) = upd(2,nlevp1)+o1_ub
    upd(3,nlevp1) = upd(3,nlevp1)+he_ub

    dpdt(1,:) = (upd(1,:)-o2_nm)/(2*step)
    dpdt(2,:) = (upd(2,:)-o1_nm)/(2*step)
    dpdt(3,:) = (upd(3,:)-he_nm)/(2*step)
    do n = 1,3
      do k = 1,nlevp1
        loss_out(n,k) = dot_product(loss(n,:,k),upd(:,k))
      enddo
      do k = 2,nlevp1-1
        moldif(n,k) = &
          (dot_product(molp(n,:,k),upd(:,k-1))- &
           dot_product(molq(n,:,k),upd(:,k  ))+ &
           dot_product(molr(n,:,k),upd(:,k+1)))/dz(k)/expzm(k)
        eddydif(n,k) = &
          (eddyp(k)*upd(n,k-1)- &
           eddyq(k)*upd(n,k  )+ &
           eddyr(k)*upd(n,k+1))/dz(k)
        veradv(n,k) = wmid(k)*(upd(n,k+1)-upd(n,k-1))/(2*dz(k))
      enddo
    enddo

! ensure non-negative O2, O, He
    o2_upd = max(upd(1,:),0.0_rp)
    o1_upd = max(upd(2,:),0.0_rp)
    he_upd = max(upd(3,:),0.0_rp)

  endsubroutine comp


!-----------------------------------------------------------------------
  pure function blktri(pk,qk,rk,fk,nk) result(upd)

    use matutil_mod,only:matinv3

    integer,intent(in) :: nk
    real(kind=rp),dimension(3,3,nk),intent(in) :: pk,qk,rk
    real(kind=rp),dimension(3,nk),intent(in) :: fk
    real(kind=rp),dimension(3,nk) :: upd

    integer :: n,k
    real(kind=rp),dimension(3) :: wkv1
    real(kind=rp),dimension(3,3) :: wkm1
    real(kind=rp),dimension(3,nk) :: zz
    real(kind=rp),dimension(3,3,nk) :: gama

    zz(:,1) = 0
    gama(:,:,1) = 0

    do k = 1,nk-1
! ALFA = Q(K)-P(K)*GAMA(K)
! ALFA refers to the block diagonal matrices,
!   and GAMA to the upper block diagonal matrices
!   in the Thomas algorithm solution
!   to the block tridiagonal system of equations
! WKM1 = INV(ALFA)
      wkm1 = matinv3(qk(:,:,k)-matmul(pk(:,:,k),gama(:,:,k)))

! WKV1 = F(K)-P(K)*Z(K)
      do n = 1,3
        wkv1(n) = fk(n,k)-dot_product(pk(n,:,k),zz(:,k))
      enddo

! GAMA(K+1) = WKM1*R(K)
      gama(:,:,k+1) = matmul(wkm1,rk(:,:,k))

! Z(K+1) = WKM1*WKV1
      do n = 1,3
        zz(n,k+1) = dot_product(wkm1(n,:),wkv1)
      enddo
    enddo

! set upper boundary to zero
    upd(:,nk) = 0

! downward sweep
    do k = nk-1,1,-1
      do n = 1,3
        upd(n,k) = zz(n,k+1)-dot_product(gama(n,:,k+1),upd(:,k+1))
      enddo
    enddo

  endfunction blktri
!-----------------------------------------------------------------------
endmodule major_mod
