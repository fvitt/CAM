module calculate_terms_module

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  pure subroutine calculate_conductance( &
    mlatd0,mlatd1,mlond0,mlond1, &
    npts_s1,npts_s2, &
    vmp_s1,bmag_s1,sigP_s1,sigH_s1, &
    vmp_s2,bmag_s2,sigP_s2,sigH_s2, &
    zigP_s1,zigH_s1,zigP_s2,zigH_s2, &
    npts_p, vmp_p, bmag_p, sigP_p, zigP_p)
! calculate field-line integrated conductance

    use params_module,only:nhgt_fix,nmlat_h,nmlatS2_h
    use cons_module,only:fill_value

    ! Dummy Args
    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_s1
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
         vmp_s1,bmag_s1,sigP_s1,sigH_s1,vmp_s2,bmag_s2,sigP_s2,sigH_s2
    real(rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
         zigP_s1,zigH_s1,zigP_s2,zigH_s2

    ! Optional Args
    integer,dimension(nmlat_h), optional ,intent(in) :: npts_p
    real(rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1), optional ,intent(in) :: &
         vmp_p, bmag_p, sigP_p
    real(rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1), optional ,intent(out) :: &
         zigP_p

    integer :: i,j,isn,k
    real(kind=rp) :: ds,sumP,sumH

    zigP_s1 = fill_value
    zigH_s1 = fill_value
    zigP_s2 = fill_value
    zigH_s2 = fill_value

    if (present(npts_p) .and. present(vmp_p) .and. present(bmag_p) .and. present(sigP_p) .and. present(zigP_p)) then
       zigP_p = fill_value
       do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
          sumP = 0
          do k = 1,npts_p(j)-1
             ds = 2*abs(vmp_p(k+1,isn,j,i)-vmp_p(k,isn,j,i))/ &
                  (bmag_p(k+1,isn,j,i)+bmag_p(k,isn,j,i))
             sumP = sumP+sigP_p(k,isn,j,i)*ds
          enddo
          zigP_p(isn,j,i) = sumP
       end do
    end if

    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
      sumP = 0
      sumH = 0
      do k = 1,npts_s1(j)-1
        ds = 2*abs(vmp_s1(k+1,isn,j,i)-vmp_s1(k,isn,j,i))/ &
          (bmag_s1(k+1,isn,j,i)+bmag_s1(k,isn,j,i))
        sumP = sumP+sigP_s1(k,isn,j,i)*ds
        sumH = sumH+sigH_s1(k,isn,j,i)*ds
      enddo
      zigP_s1(isn,j,i) = sumP
      zigH_s1(isn,j,i) = sumH
    enddo

    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlatS2_h)
      sumP = 0
      sumH = 0
      do k = 1,npts_s2(j)-1
        ds = 2*abs(vmp_s2(k+1,isn,j,i)-vmp_s2(k,isn,j,i))/ &
          (bmag_s2(k+1,isn,j,i)+bmag_s2(k,isn,j,i))
        sumP = sumP+sigP_s2(k,isn,j,i)*ds
        sumH = sumH+sigH_s2(k,isn,j,i)*ds
      enddo
      zigP_s2(isn,j,i) = sumP
      zigH_s2(isn,j,i) = sumH
    enddo

  endsubroutine calculate_conductance
!-----------------------------------------------------------------------
  pure subroutine calculate_n( &
    mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
    D_s1,M1_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1, &
    D_s2,M2_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2, &
    N1p_s1,N1h_s1,N2p_s2,N2h_s2)
! Updated 2015/08/23

    use params_module,only:nhgt_fix,nmlat_h,nmlatS2_h,ylonm,rho,rho_s
    use cons_module,only:r0,fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_s1
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      D_s1,M1_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1, &
      D_s2,M2_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
      N1p_s1,N1h_s1,N2p_s2,N2h_s2

    integer :: i,j,isn,k
    real(kind=rp) :: dlonm,drho,sigC

    N1p_s1 = fill_value
    N1h_s1 = fill_value
    N2p_s2 = fill_value
    N2h_s2 = fill_value

! assume equidistant longitudinal grid points
    dlonm = ylonm(2)-ylonm(1)

! calculate N coefficients for S1 points, these are at (i+0.5,j,k)
! do not calculate N1p & N1h for the pole (not used)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h)
      if (j == nmlat_h) then
        drho = 2*(rho(j)-rho(j-1))
      else
        drho = rho(j+1)-rho(j-1)
      endif

      do concurrent (k = 1:npts_s1(j))

! N1P(i+0.5) = M1(i+0.5)*[sigP*d1^2](i+0.5)/R/rho(j)/(phi(i+1)-phi(i))
        N1p_s1(k,isn,j,i) = M1_s1(k,isn,j,i)* &
          sigP_s1(k,isn,j,i)*d1d1_s1(k,isn,j,i)/r0/rho(j)/dlonm

! N1H(i+0.5) = M1(i+0.5)*[sigH*D-sigP*d1*d2](i+0.5)*sqrt(1-0.75*rho^2(j))/2/R/(rho(j+1)-rho(j-1))
        N1h_s1(k,isn,j,i) = M1_s1(k,isn,j,i)* &
          (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)- &
          sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))* &
          sqrt(1-3*rho(j)**2/4)/2/r0/drho
      enddo
    enddo

! lowest equatorial volume for i+0.5,j,k
! overwrite values from above N1P and N1H (page 12 2014/01/30 Art's notes)
! N1H = 0
! N1P -> N1C
! N1C(i+0.5,j,k) = M1(i+0.5,j,k)*sigC(i+0.5,j,k)/R/rho(j)/(phi(i+1)-phi(i))
! sigC = sigP*d1^2+(sigH*D-sigP*d1*d2)*(sigH*D+sigP*d1*d2)/sigP/(d2*d2)
    j = nmlat_h
    k = 1
    if (j>=mlatd0 .and. j<=mlatd1) then
      do concurrent (i = mlond0:mlond1, isn = 1:2)
        N1h_s1(k,isn,j,i) = 0
        sigC = sigP_s1(k,isn,j,i)*d1d1_s1(k,isn,j,i)+ &
          (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)- &
          sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))* &
          (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)+ &
          sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))/ &
          sigP_s1(k,isn,j,i)/d2d2_s1(k,isn,j,i)
        N1p_s1(k,isn,j,i) = M1_s1(k,isn,j,i)*sigC/r0/rho(j)/dlonm
      enddo
    endif

! calculate N coefficients for S2 points, these are at (i,j+0.5,k)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlatS2_h)
      drho = rho(j+1)-rho(j)

      do concurrent (k = 1:npts_s2(j))

! N2H(j+0.5) = M2(j+0.5)*[sigH*D+sigP*d1*d2](j+0.5)/4/R/rho(j+0.5)/(phi(i+1)-phi(i)))
        N2h_s2(k,isn,j,i) = M2_s2(k,isn,j,i)* &
          (sigH_s2(k,isn,j,i)*D_s2(k,isn,j,i)+ &
          sigP_s2(k,isn,j,i)*d1d2_s2(k,isn,j,i))/ &
          4/r0/rho_s(j)/dlonm

! N2P(j+0.5) = M2(j+0.5)[sigP*d2^2](j+0.5)*sqrt(1-0.75*rho^2(j+0.5))/R/(rho(j+1)-rho(j))
        N2p_s2(k,isn,j,i) = M2_s2(k,isn,j,i)* &
          sigP_s2(k,isn,j,i)*d2d2_s2(k,isn,j,i)* &
          sqrt(1-3*rho_s(j)**2/4)/r0/drho
      enddo
    enddo

  endsubroutine calculate_n
!-----------------------------------------------------------------------
  pure subroutine calculate_je( &
    mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
    D_s1,be3_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1,un_s1,vn_s1, &
    D_s2,be3_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2,un_s2,vn_s2, &
    d1_s1,d2_s1,d1_s2,d2_s2,Je1D_s1,Je2D_s2)

    use params_module,only:nhgt_fix,nmlat_h,nmlatS2_h
    use cons_module,only:J3LB,fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_s1
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      D_s1,be3_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1,un_s1,vn_s1, &
      D_s2,be3_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2,un_s2,vn_s2
    real(kind=rp),dimension(3,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      d1_s1,d2_s1,d1_s2,d2_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
      Je1D_s1,Je2D_s2

    integer :: i,j,isn,k
    real(kind=rp) :: ue1,ue2,fac

    Je1D_s1 = fill_value
    Je2D_s2 = fill_value

! at i+0.5,j,k calculate Je1D = sigP*d1^2*ue2*Be3+(sigH*D-sigP*d1*d2)*ue1*Be3+Je1^Ion
! calculate values for S1 points, these are at (i+0.5,j,k)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h) ! no pole
      do concurrent (k = 1:npts_s1(j))
        ue1 = un_s1(k,isn,j,i)*d1_s1(1,k,isn,j,i)+ &
              vn_s1(k,isn,j,i)*d1_s1(2,k,isn,j,i)
        ue2 = un_s1(k,isn,j,i)*d2_s1(1,k,isn,j,i)+ &
              vn_s1(k,isn,j,i)*d2_s1(2,k,isn,j,i)
        Je1D_s1(k,isn,j,i) = &
          sigP_s1(k,isn,j,i)*d1d1_s1(k,isn,j,i)* &
          ue2*be3_s1(k,isn,j,i)+ &
          (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)- &
          sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))* &
          ue1*be3_s1(k,isn,j,i)
      enddo
    enddo

! lowest equatorial volume at i+0.5,j,k
! calculate Je2 = (sigH*D+sigP*d1*d2)*ue2*Be3-sigP*d2^2*ue1*Be3+Je1^Ion
! Je1D -> Je1S
! with Je1S = Je1D - (sigH*D-sigP*d1*d2)/sigP/d2^2*(Je2LB-Je2D) at (i+0.5,j,k)
! Je2LB is the Je2 given through coupling with the lower atmosphere
! we assume that this is Je2LB(i) = -J3LB(i,nmlat_h)
! not exactly since there is half a height level in between but should be close
! could do for only one hemisphere (since the same point) and then copy into other hemisphere
    j = nmlat_h
    k = 1
    if (j>=mlatd0 .and. j<=mlatd1) then
      do i = mlond0,mlond1
        do isn = 1,2
          ue1 = un_s1(k,isn,j,i)*d1_s1(1,k,isn,j,i)+ &
                vn_s1(k,isn,j,i)*d1_s1(2,k,isn,j,i)
          ue2 = un_s1(k,isn,j,i)*d2_s1(1,k,isn,j,i)+ &
                vn_s1(k,isn,j,i)*d2_s1(2,k,isn,j,i)

! Je2D at the equator (there is no S2 point therefore needs to be calculated)
          fac = (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)+ &
            sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))* &
            ue2*be3_s1(k,isn,j,i)- &
            sigP_s1(k,isn,j,i)*d2d2_s1(k,isn,j,i)* &
            ue1*be3_s1(k,isn,j,i)

! Je1S = Je1D - (sigH*D-sigP*d1*d2)/sigP/d2^2*(Je2LB-Je2D) at (i+0.5,j,k)
! we assume that this is Je2LB(i)= -J3LB(i,nmlat_h)
! Correction from H. Wu: change calculation based on comments
! sigP_s1(k,isn,j,i)*d2d2_s1(k,isn,j,i) is now corrected as
! sigP_s1(k,isn,j,i)/d2d2_s1(k,isn,j,i)
          Je1D_s1(k,isn,j,i) = Je1D_s1(k,isn,j,i)- &
            (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)- &
            sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))/ &
            sigP_s1(k,isn,j,i)/d2d2_s1(k,isn,j,i)* &
            (-J3LB(isn,j,i)-fac)
        enddo
      enddo
    endif

! at i,j+0.5,k calculate Je2 = (sigH*D+sigP*d1*d2)*ue2*Be3-sigP*d2^2*ue1*Be3+Je1^Ion
! calculate values for S2 points, these are at (i,j+0.5,k)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlatS2_h)
      do concurrent (k = 1:npts_s2(j))
        ue1 = un_s2(k,isn,j,i)*d1_s2(1,k,isn,j,i)+ &
              vn_s2(k,isn,j,i)*d1_s2(2,k,isn,j,i)
        ue2 = un_s2(k,isn,j,i)*d2_s2(1,k,isn,j,i)+ &
              vn_s2(k,isn,j,i)*d2_s2(2,k,isn,j,i)
        Je2D_s2(k,isn,j,i) = &
          (sigH_s2(k,isn,j,i)*D_s2(k,isn,j,i)+ &
          sigP_s2(k,isn,j,i)*d1d2_s2(k,isn,j,i))* &
          ue2*be3_s2(k,isn,j,i)- &
          sigP_s2(k,isn,j,i)*d2d2_s2(k,isn,j,i)* &
          ue1*be3_s2(k,isn,j,i)
      enddo
    enddo

  endsubroutine calculate_je
!-----------------------------------------------------------------------
  pure function calculate_s(mlatd0,mlatd1,mlond0,mlond1,npts_p, &
    M1_s1,Je1D_s1,M2_s2,Je2D_s2,M3_r) result(S_p)
! S is the wind driven and ionospheric current sources (89)
! I3^E are the external current sources from the magnetosphere I3^M and from the lower atmosphere
! see Eq (108) I3^E = -I3^M(NH+SH) + I3_lb(NH) + I3_lb(SH)
! Eq (115) I3_total = sum_k[S(NH) + S(SH)]  + I3^E

    use params_module,only:nhgt_fix,nhgt_fix_r,nmlat_h
    use cons_module,only:J3LB,fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_p
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      M1_s1,Je1D_s1,M2_s2,Je2D_s2
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: M3_r
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1) :: S_p

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,k

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    S_p = fill_value

! page 12 Eq (89) Art's script
! -S(i,j,k) = M1(i-0.5,j,k)*Je1D(i-0.5,j,k)-M1(i+0.5,j,k)*Je1D(i+0.5,j,k)+
!             M2(i,j-0.5,k)*Je2D(i,j-0.5,k)-M2(i,j+0.5,k)*Je2D(i,j+0.5,k)
! the -S(i,j,k) denotes that it is on the left hand side,
! but for the dynamo equation we need it on the right hand side,
! therefore multiply by -1
! for equatorial volumes at different heights M2(i,j+0.5,k) = 0

! the following is taken care of in calculate_je
! for lowest equatorial volume
! S = M1(i-0.5,j,k)*Je1S(i-0.5,j,k) - M1(i+0.5,j,k)*Je1S(i+0.5,j,k) + M2(i,j-0.5,k)*Je2D(i,j-0.5,k)
! with Je1S = Je1D - (sigH*D-sigP*d1*d2)/sigP/d2^2*(Je2LB-Je2D) at (i+0.5,j,k)
! Je2LB is the Je2 given through coupling with the lower atmosphere
! at the moment it is set to zero

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)
      do concurrent (k = 1:npts_p(j))

! A. Maute 2023/01: set pole value (no value set before)
        if (j == 1) then

! S(i,1,k) = -M2(i,3/2,k)*Je2D(i,3/2,k) Art Jan 16 2023 page 16 eq(322)
! change the sign of S since in the code S is on the RHS but in the write up it is on the LHS
! north pole value are not really used but south pole
          S_p(k,isn,j,i) = M2_s2(k,isn,j,i)*Je2D_s2(k,isn,j,i)

! j >= 2, pole seperate
        else

! top volume at equator M2(i,j+0.5,k) = 0
          S_p(k,isn,j,i) = &
            -(M1_s1(k,isn,j  ,i-1)*Je1D_s1(k,isn,j  ,i-1)- &
              M1_s1(k,isn,j  ,i  )*Je1D_s1(k,isn,j  ,i  )+ &
              M2_s2(k,isn,j-1,i  )*Je2D_s2(k,isn,j-1,i  ))
          if (k /= nmlat_h-j+1) &
            S_p(k,isn,j,i) = S_p(k,isn,j,i)+ &
            M2_s2(k,isn,j,i)*Je2D_s2(k,isn,j,i)
        endif
      enddo

! add the current from the lower atmosphere J3_lb*M3 at R point for k=0.5 index=1
! I3S_lb = J3LB*M3 -> change sign here as well since S is on RHS
! H. Wu: this is not consistent with below, comment it out for now
!      S_p(1,isn,j,i) = S_p(1,isn,j,i)-J3LB(isn,j,i)*M3_r(1,isn,j,i)

! add the current from the lower atmosphere J3_lb*M3 at R point for k=0.5 index=1
! I3S_lb = J3LB*M3
      S_p(1,isn,j,i) = S_p(1,isn,j,i)+J3LB(isn,j,i)*M3_r(1,isn,j,i)
    enddo

  endfunction calculate_s
!-----------------------------------------------------------------------
  function correct_fac_hl(mlatd0,mlatd1,mlond0,mlond1, &
    zigP_p,M3_p,fac_hl_in_p) result(fac_hl_out_p)
! correct the input high latitude FAC to make sure it is zero
! when integrated in each hemisphere

! this is called when direct FAC is read in (read_fac==.true.)
! and needs to be corrected (direct FAC input may be unbalanced)
! M3_p is only the bottom level

    use params_module,only:nmlat_h
    use mpi_module,only:reduce_sum_1d

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      zigP_p,M3_p,fac_hl_in_p
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1) :: fac_hl_out_p

    real(kind=rp),parameter :: thres = 1.5_rp
    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn
    real(kind=rp) :: facArea
    real(kind=rp),dimension(2) :: sumfac,sumzigP,corr
    real(kind=rp),dimension(4) :: tmp_sub,tmp_full

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    fac_hl_out_p = fac_hl_in_p

! corr = - zigP*abs(Jmr)/sinI * [sum_i^N Jmr*area] / [sum_i^N abs(Jmr)*zigP/sinI*area]
! sinI is ignored by assuming it is close to 1 at high latitudes
    sumfac = 0
    sumzigP = 0

! exclude halo points to avoid double counting
    do i = mlon0,mlon1
      do j = mlat0,mlat1
        if (j>=2 .and. j<=nmlat_h) then ! no pole
          do isn = 1,2
            if (zigP_p(isn,j,i) > thres) then
              facArea = fac_hl_out_p(isn,j,i)*M3_p(isn,j,i)
              sumfac(isn) = sumfac(isn)+facArea
              sumzigP(isn) = sumzigP(isn)+zigP_p(isn,j,i)*abs(facArea)
            else
              fac_hl_out_p(isn,j,i) = 0
            endif
          enddo
        endif
      enddo
    enddo

    tmp_sub(1) = sumfac(1)
    tmp_sub(2) = sumfac(2)
    tmp_sub(3) = sumzigP(1)
    tmp_sub(4) = sumzigP(2)
    tmp_full = reduce_sum_1d(tmp_sub,4,-1)
    sumfac(1) = tmp_full(1)
    sumfac(2) = tmp_full(2)
    sumzigP(1) = tmp_full(3)
    sumzigP(2) = tmp_full(4)

    do concurrent (isn = 1:2)
      corr(isn) = sumfac(isn)/sumzigP(isn)
    enddo

! correct fac_hl
    do i = mlond0,mlond1
      do j = mlatd0,mlatd1
        if (j>=2 .and. j<=nmlat_h) then ! no pole (potential set later)
          do isn = 1,2
            fac_hl_out_p(isn,j,i) = fac_hl_out_p(isn,j,i)- &
              zigP_p(isn,j,i)*abs(fac_hl_out_p(isn,j,i))*corr(isn)
          enddo
        endif
      enddo
    enddo

  endfunction correct_fac_hl
!-----------------------------------------------------------------------
  pure subroutine calculate_ed(mlatd0,mlatd1,mlond0,mlond1, &
    pot_p,ed1_s1,ed2_s1,ed1_s2,ed2_s2)
! calculates electric field Ed1,Ed2 at S1 and S2 points

    use params_module,only:nmlat_h,nmlatS2_h,ylonm,rho,rho_s
    use cons_module,only:r0,fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_p
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
      ed1_s1,ed2_s1,ed1_s2,ed2_s2

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn
    real(kind=rp) :: fac,facj

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    ed1_s1 = fill_value
    ed2_s1 = fill_value
    ed1_s2 = fill_value
    ed2_s2 = fill_value

    fac = 1/(r0*(ylonm(2)-ylonm(1))) ! regular spaced in longitude

! S1 loop
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)
      if (j == 1) then ! pole
        facj = sqrt(1-3*rho(j)**2/4)/r0/2/(rho(j+1)-rho(j))
        ed1_s1(isn,j,i) = (pot_p(isn,j,i)-pot_p(isn,j,i+1))*fac/rho(j+1)
        ed2_s1(isn,j,i) = facj* &
          (pot_p(isn,j  ,i)+pot_p(isn,j  ,i+1)- &
           pot_p(isn,j+1,i)-pot_p(isn,j+1,i+1))
      elseif (j == nmlat_h) then ! equator
        facj = sqrt(1-3*rho(j)**2/4)/r0/2/(rho(j)-rho(j-1))
        ed1_s1(isn,j,i) = (pot_p(isn,j,i)-pot_p(isn,j,i+1))*fac/rho(j)
        ed2_s1(isn,j,i) = facj* &
          (pot_p(isn,j-1,i)+pot_p(isn,j-1,i+1)- &
           pot_p(isn,j  ,i)-pot_p(isn,j  ,i+1))
      else ! not the pole or equator
        facj = sqrt(1-3*rho(j)**2/4)/r0/2/(rho(j+1)-rho(j-1))
        ed1_s1(isn,j,i) = (pot_p(isn,j,i)-pot_p(isn,j,i+1))*fac/rho(j)
        ed2_s1(isn,j,i) = facj* &
          (pot_p(isn,j-1,i)+pot_p(isn,j-1,i+1)- &
           pot_p(isn,j+1,i)-pot_p(isn,j+1,i+1))
      endif
    enddo

! S2 loop
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j <= nmlatS2_h)
      facj = sqrt(1-3*rho_s(j)**2/4)/r0/(rho(j+1)-rho(j))
      ed1_s2(isn,j,i) = fac/4/rho_s(j)* &
        (pot_p(isn,j,i-1)+pot_p(isn,j+1,i-1)- &
         pot_p(isn,j,i+1)-pot_p(isn,j+1,i+1))
      ed2_s2(isn,j,i) = facj*(pot_p(isn,j,i)-pot_p(isn,j+1,i))
    enddo

  endsubroutine calculate_ed
!-----------------------------------------------------------------------
  pure subroutine calculate_ve( &
    mlatd0,mlatd1,mlond0,mlond1, &
    ed1_s1,ed2_s1,be3_s1, &
    ed1_s2,ed2_s2,be3_s2, &
    ve1_s1,ve2_s1,ve1_s2,ve2_s2)
! calculates drift velocity ve1,ve2 at S1 and S2 points
! ve1 = Ed2/Be3   &   ve2 = -Ed1/Be3
! be3 is only the bottom level

    use params_module,only:nmlatS2_h
    use cons_module,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      ed1_s1,ed2_s1,be3_s1,ed1_s2,ed2_s2,be3_s2
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
      ve1_s1,ve2_s1,ve1_s2,ve2_s2

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    ve1_s1 = fill_value
    ve2_s1 = fill_value
    ve1_s2 = fill_value
    ve2_s2 = fill_value

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)
      ve1_s1(isn,j,i) =  ed2_s1(isn,j,i)/be3_s1(isn,j,i)
      ve2_s1(isn,j,i) = -ed1_s1(isn,j,i)/be3_s1(isn,j,i)
    enddo

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j <= nmlatS2_h)
      ve1_s2(isn,j,i) =  ed2_s2(isn,j,i)/be3_s2(isn,j,i)
      ve2_s2(isn,j,i) = -ed1_s2(isn,j,i)/be3_s2(isn,j,i)
    enddo

  endsubroutine calculate_ve
!-----------------------------------------------------------------------
  pure subroutine calculate_exyz( &
    mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
    ed1_s1,ed2_s1,d1_s1,d2_s1, &
    ed1_s2,ed2_s2,d1_s2,d2_s2, &
    ex_s1,ey_s1,ez_s1,ex_s2,ey_s2,ez_s2)
! get electric fields in geographic coordinates (Ed1,2 -> Ex,y,z)

    use params_module,only:nhgt_fix,nmlat_h,nmlatS2_h
    use cons_module,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_s1
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      ed1_s1,ed2_s1,ed1_s2,ed2_s2
    real(kind=rp),dimension(3,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      d1_s1,d2_s1,d1_s2,d2_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
      ex_s1,ey_s1,ez_s1,ex_s2,ey_s2,ez_s2

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,k

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    ex_s1 = fill_value
    ey_s1 = fill_value
    ez_s1 = fill_value
    ex_s2 = fill_value
    ey_s2 = fill_value
    ez_s2 = fill_value

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)
      do concurrent (k = 1:npts_s1(j))
        ex_s1(k,isn,j,i) = &
          ed1_s1(isn,j,i)*d1_s1(1,k,isn,j,i)+ &
          ed2_s1(isn,j,i)*d2_s1(1,k,isn,j,i)
        ey_s1(k,isn,j,i) = &
          ed1_s1(isn,j,i)*d1_s1(2,k,isn,j,i)+ &
          ed2_s1(isn,j,i)*d2_s1(2,k,isn,j,i)
        ez_s1(k,isn,j,i) = &
          ed1_s1(isn,j,i)*d1_s1(3,k,isn,j,i)+ &
          ed2_s1(isn,j,i)*d2_s1(3,k,isn,j,i)
      enddo
    enddo

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j <= nmlatS2_h)
      do concurrent (k = 1:npts_s2(j))
        ex_s2(k,isn,j,i) = &
          ed1_s2(isn,j,i)*d1_s2(1,k,isn,j,i)+ &
          ed2_s2(isn,j,i)*d2_s2(1,k,isn,j,i)
        ey_s2(k,isn,j,i) = &
          ed1_s2(isn,j,i)*d1_s2(2,k,isn,j,i)+ &
          ed2_s2(isn,j,i)*d2_s2(2,k,isn,j,i)
        ez_s2(k,isn,j,i) = &
          ed1_s2(isn,j,i)*d1_s2(3,k,isn,j,i)+ &
          ed2_s2(isn,j,i)*d2_s2(3,k,isn,j,i)
      enddo
    enddo

  endsubroutine calculate_exyz
!-----------------------------------------------------------------------
  pure subroutine calculate_vxyz( &
    mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
    ve1_s1,ve2_s1,e1_s1,e2_s1, &
    ve1_s2,ve2_s2,e1_s2,e2_s2, &
    vx_s1,vy_s1,vz_s1,vx_s2,vy_s2,vz_s2)
! get drift velocities in geographic coordinates (Ve1,2 -> Vx,y,z)

    use params_module,only:nhgt_fix,nmlat_h,nmlatS2_h
    use cons_module,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_s1
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      ve1_s1,ve2_s1,ve1_s2,ve2_s2
    real(kind=rp),dimension(3,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      e1_s1,e2_s1,e1_s2,e2_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
      vx_s1,vy_s1,vz_s1,vx_s2,vy_s2,vz_s2

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,k

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    vx_s1 = fill_value
    vy_s1 = fill_value
    vz_s1 = fill_value
    vx_s2 = fill_value
    vy_s2 = fill_value
    vz_s2 = fill_value

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)
      do concurrent (k = 1:npts_s1(j))
        vx_s1(k,isn,j,i) = &
          ve1_s1(isn,j,i)*e1_s1(1,k,isn,j,i)+ &
          ve2_s1(isn,j,i)*e2_s1(1,k,isn,j,i)
        vy_s1(k,isn,j,i) = &
          ve1_s1(isn,j,i)*e1_s1(2,k,isn,j,i)+ &
          ve2_s1(isn,j,i)*e2_s1(2,k,isn,j,i)
        vz_s1(k,isn,j,i) = &
          ve1_s1(isn,j,i)*e1_s1(3,k,isn,j,i)+ &
          ve2_s1(isn,j,i)*e2_s1(3,k,isn,j,i)
      enddo
    enddo

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j <= nmlatS2_h)
      do concurrent (k = 1:npts_s2(j))
        vx_s2(k,isn,j,i) = &
          ve1_s2(isn,j,i)*e1_s2(1,k,isn,j,i)+ &
          ve2_s2(isn,j,i)*e2_s2(1,k,isn,j,i)
        vy_s2(k,isn,j,i) = &
          ve1_s2(isn,j,i)*e1_s2(2,k,isn,j,i)+ &
          ve2_s2(isn,j,i)*e2_s2(2,k,isn,j,i)
        vz_s2(k,isn,j,i) = &
          ve1_s2(isn,j,i)*e1_s2(3,k,isn,j,i)+ &
          ve2_s2(isn,j,i)*e2_s2(3,k,isn,j,i)
      enddo
    enddo

  endsubroutine calculate_vxyz
!-----------------------------------------------------------------------
  subroutine calculate_current( &
    mlatd0,mlatd1,mlond0,mlond1, &
    npts_p,npts_s1,npts_s2,npts_r, &
    pot_p,M3_p,M1_s1,N1p_s1,N1h_s1,Je1D_s1, &
    M2_s2,N2p_s2,N2h_s2,Je2D_s2,M3_r, &
    I1_s1,I2_s2,I3_r,Jr_p,Jr_r)
! Eq (83) page 10 Art's script
! I1(i+0.5,j,k) = N1p(i+0.5,j,k)*[Phi(i,j)-Phi(i+1,j)]
!                -N1h(i+0.5,j,k)*[Phi(i,j-1)+Phi(i+1,j-1)-Phi(i,j+1)-Phi(i+1,j+1)]
!                +M1 (i+0.5,j,k)*Je1D(i+0.5,j,k)

! Eq (85) page 10 Art's script
! I2(i,j+0.5,k) = N2h(i,j+0.5,k)*[Phi(i-1,j)+Phi(i-1,j+1)-Phi(i+1,j)-Phi(i+1,j+1)]
!                +N2p(i,j+0.5,k)*[Phi(i,j)-Phi(i,j+1)]
!                +M2 (i,j+0.5,k)*Je2D(i,j+0.5,k)

! Eq (63') page 7 Art's script
! for all i and j=2,nmlat_h
! I3(i,j,k+0.5) = I3(i,j,k-0.5) + I1(i-0.5,j,k)-I1(i+0.5,j,k)+I2(i,j-0.5,k)-I2(i,j+0.5,k)
! for k=1, the lowest level I3(i,j,0.5) is given
! in variable J3LB (from lower atmosphere) in parms

! how does the index in notes relate to the index in code
!    points       notes             code       quantities
! P  points   i    ,j    ,k       i  ,j  ,k   potential, S
! S1 points   i-0.5,j    ,k       i-1,j  ,k
! S2 points   i    ,j-0.5,k       i  ,j-1,k
! R  points   i    ,j    ,k-0.5   i  ,j  ,k

!    points       notes           code         quantities
! P  points   i    ,j    ,k       i,j,k       potential, S
! S1 points   i+0.5,j    ,k       i,j,k
! S2 points   i    ,j+0.5,k       i,j,k
! R  points   i    ,j    ,k-0.5   i,j,k

    use params_module,only:nhgt_fix,nhgt_fix_r,nmlat_h,nmlatS2_h,nmlon
    use cons_module,only:J3LB,fill_value
    use mpi_module,only:gather_mlon_3d,sync_mlat_5d,sync_mlon_5d

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_p,npts_s1,npts_r
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_p
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      M3_p,M1_s1,N1p_s1,N1h_s1,Je1D_s1,M2_s2,N2p_s2,N2h_s2,Je2D_s2
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: M3_r
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: I1_s1,I2_s2,Jr_p
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: I3_r,Jr_r

    integer :: mlat0,mlat1,mlon0,mlon1,isn,i,j,k,iconj
    real(kind=rp),dimension(nhgt_fix,2,mlond0+1:mlond1-1) :: Je1_sub
    real(kind=rp),dimension(nhgt_fix,2,nmlon) :: Je1_full
    real(kind=rp),dimension(nhgt_fix,2,mlatd0+1:mlatd1-1,mlond0+1:mlond1-1) :: &
      I1_1,I1_2,I1_3,I2_1,I2_2,I2_3
    real(kind=rp),dimension(2,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1) :: tmpI

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    I1_s1 = fill_value
    I2_s2 = fill_value
    I3_r = fill_value
    Jr_p = fill_value
    Jr_r = fill_value

! for the lowest equatorial volume, j=nmlat_h and k=1, Eq (93) Art's notes
! I(i-0.5,j,k) = N_1^C(i-0.5,j,k)[Phi(i-1,j)-Phi(i,j)]+M1(i-0.5,j,k)Je1^S(i-0.5,j,k)
! N1^c(i-0.5,j,k) = M1(i-0.5,j,k)[sig_c](i-0.5,j,k)/R/rhoj(phi_i-phi_i-1)
! N1^c saved in N1^p
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j >= 2) ! no pole
      do concurrent (k = 1:npts_s1(j))
        if (j == nmlat_h) then
          I1_1(k,isn,j,i) = 0
        else
          I1_1(k,isn,j,i) = -N1h_s1(k,isn,j,i)* &
            (pot_p(isn,j-1,i)+pot_p(isn,j-1,i+1)- &
             pot_p(isn,j+1,i)-pot_p(isn,j+1,i+1))
        endif
        I1_2(k,isn,j,i) = N1p_s1(k,isn,j,i)*(pot_p(isn,j,i)-pot_p(isn,j,i+1))
        I1_3(k,isn,j,i) = M1_s1(k,isn,j,i)*Je1D_s1(k,isn,j,i)
        I1_s1(k,isn,j,i) = I1_1(k,isn,j,i)+I1_2(k,isn,j,i)+I1_3(k,isn,j,i)
      enddo
    enddo

! pole value see Eq (216') page 1 2015/04/19
! Je1(i-0.5,1,k) = 0.5*[Je1(i-0.5,2,k) - Je1(iconj-0.5,2,k)]
! iconj is the conjugate longitude
! I1(i-0.5,1,k) = Je1(i-0.5,1,k)*M1(i-0.5,1,k)
! we do not have iconj in the current process, so first gather all longitudes
    j = mlat0
    Je1_sub = fill_value
    do concurrent (i = mlon0:mlon1, isn = 1:2, k = 1:npts_s1(j))
      Je1_sub(k,isn,i) = I1_s1(k,isn,j+1,i)/M1_s1(k,isn,j+1,i) ! Je1(i-0.5,2,k)
    enddo
    Je1_full = gather_mlon_3d(Je1_sub,nhgt_fix,2)
    if (j == 1) then
      do concurrent (i = mlon0:mlon1, isn = 1:2, k = 1:npts_s1(j))
        if (i > nmlon/2) then
          iconj = i-nmlon/2
        else
          iconj = i+nmlon/2
        endif
        I1_s1(k,isn,j,i) = (Je1_full(k,isn,i)-Je1_full(k,isn,iconj))/2*M1_s1(k,isn,j,i)
      enddo
    endif

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j <= nmlatS2_h)
      do concurrent (k = 1:npts_s2(j))
        I2_1(k,isn,j,i) = N2h_s2(k,isn,j,i)* &
          (pot_p(isn,j,i-1)+pot_p(isn,j+1,i-1)- &
           pot_p(isn,j,i+1)-pot_p(isn,j+1,i+1))
        I2_2(k,isn,j,i) = N2p_s2(k,isn,j,i)*(pot_p(isn,j,i)-pot_p(isn,j+1,i))
        I2_3(k,isn,j,i) = M2_s2(k,isn,j,i)*Je2D_s2(k,isn,j,i)
        I2_s2(k,isn,j,i) = I2_1(k,isn,j,i)+I2_2(k,isn,j,i)+I2_3(k,isn,j,i)
      enddo
    enddo

! I1 at i-1 and I2 at j-1 are used in I3 calculation, so sync before used
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, k = 1:nhgt_fix)
      tmpI(1,k,isn,j,i) = I1_s1(k,isn,j,i)
      tmpI(2,k,isn,j,i) = I2_s2(k,isn,j,i)
    enddo
    call sync_mlat_5d(tmpI(:,:,:,:,mlon0:mlon1),2,nhgt_fix,2)
    call sync_mlon_5d(tmpI,2,nhgt_fix,2)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, k = 1:nhgt_fix)
      I1_s1(k,isn,j,i) = tmpI(1,k,isn,j,i)
      I2_s2(k,isn,j,i) = tmpI(2,k,isn,j,i)
    enddo

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)

! lowest level given by lower atmosphere coupling: J3LB [A/m2], M3 [m2]
      I3_r(1,isn,j,i) = J3LB(isn,j,i)*M3_r(1,isn,j,i)

! 2015/10/14 Make calculation of I3 at pole consistent with I1 values
      do k = 2,npts_r(j)
        I3_r(k,isn,j,i) = I3_r(k-1,isn,j,i)+ &
          I1_s1(k-1,isn,j,i-1)-I1_s1(k-1,isn,j,i)- &
          I2_s2(k-1,isn,j,i)
        if (j /= 1) I3_r(k,isn,j,i) = I3_r(k,isn,j,i)+I2_s2(k-1,isn,j-1,i)
      enddo
    enddo

! Jr(i,j,k) = I3(i,j,k)/M3(i,j,k) (122)
! but top volume at equator with k=km
! Jr(i,j,k) = sqrt(2) Ir(i,j,k)/M3(i,j,k-0.5) (124)
!   and Ir = (0.5-0.5^1.5)[I1(i-0.5,j,k)-I1(i+0.5,j,k)] -(0.5^0.5-0.5)I2(i,j-0.5,k)+0.5*I3(i,j,k-0.5) (123)
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j >= 2)
      k = npts_p(j)
      Jr_p(k,isn,j,i) = I3_r(k,isn,j,i)/(sqrt(2.0_rp)*M3_r(k,isn,j,i))
      do concurrent (k = 1:npts_p(j)-1)
        Jr_p(k,isn,j,i) = (I3_r(k,isn,j,i)+I3_r(k+1,isn,j,i))/(2*M3_p(k,isn,j,i))
      enddo
      do concurrent (k = 1:npts_r(j))
        Jr_r(k,isn,j,i) = I3_r(k,isn,j,i)/M3_r(k,isn,j,i)
      enddo
    enddo

  endsubroutine calculate_current
!-----------------------------------------------------------------------
endmodule calculate_terms_module
