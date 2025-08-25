module calculate_terms_module

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  pure subroutine calculate_conductance( &
    mlatd0,mlatd1,mlond0,mlond1, &
    npts_p,vmp_p,bmag_p,sigP_p,zigP_p)
! calculate field-line integrated conductance

    use params_module,only:nhgt_fix,nmlat_h
    use cons_module,only:fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_p
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: vmp_p,bmag_p,sigP_p
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: zigP_p

    integer :: i,j,isn,k
    real(kind=rp) :: sumP

    zigP_p = fill_value

    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
      sumP = 0
      do k = 1,npts_p(j)-1
        sumP = sumP+ &
          (sigP_p(k+1,isn,j,i)+sigP_p(k,isn,j,i))* &
          abs(vmp_p(k+1,isn,j,i)-vmp_p(k,isn,j,i))/ &
          (bmag_p(k+1,isn,j,i)+bmag_p(k,isn,j,i))
      enddo
      zigP_p(isn,j,i) = sumP
    enddo

  endsubroutine calculate_conductance
!-----------------------------------------------------------------------
  pure subroutine calculate_n( &
    mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
    D_s1,M1_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1, &
    D_s2,M2_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2, &
    N1p_s1,N1h_s1,N2p_s2,N2h_s2)

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
    real(kind=rp) :: dlonm,sigC

    N1p_s1 = fill_value
    N1h_s1 = fill_value
    N2p_s2 = fill_value
    N2h_s2 = fill_value

! assume equidistant longitudinal grid points
    dlonm = ylonm(2)-ylonm(1)

! calculate N coefficients at (i+0.5,j,k) for S1 points (no pole or equator)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h-1)
      do concurrent (k = 1:npts_s1(j))

! Equation (6.11)
! N1P(i+0.5) = M1(i+0.5)*[sigP*d1^2](i+0.5)/R/rho(j)/(phi(i+1)-phi(i))
        N1p_s1(k,isn,j,i) = M1_s1(k,isn,j,i)* &
          sigP_s1(k,isn,j,i)*d1d1_s1(k,isn,j,i)/r0/rho(j)/dlonm

! Equation (6.12)
! N1H(i+0.5) = M1(i+0.5)*[sigH*D-sigP*d1*d2](i+0.5)*sqrt(1-0.75*rho(j)^2)/2/R/(rho(j+1)-rho(j-1))
        N1h_s1(k,isn,j,i) = M1_s1(k,isn,j,i)* &
          (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)- &
          sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))* &
          sqrt(1-3*rho(j)**2/4)/2/r0/(rho(j+1)-rho(j-1))
      enddo
    enddo

! Equations (6.16) and (6.17) for the equatorial volume at (i+0.5,j=J,k=1)
! N1P <- M1(i+0.5)*sigC(i+0.5)/R/rho(j)/(phi(i+1)-phi(i))
! sigC = sigP*d1^2+((sigH*D)^2-(sigP*d1*d2)^2)/(sigP*d2^2)
! N1H <- 0
    j = nmlat_h
    k = 1
    if (j>=mlatd0 .and. j<=mlatd1) then
      do concurrent (i = mlond0:mlond1, isn = 1:2)
        sigC = sigP_s1(k,isn,j,i)*d1d1_s1(k,isn,j,i)+ &
          ((sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i))**2- &
          (sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))**2)/ &
          (sigP_s1(k,isn,j,i)*d2d2_s1(k,isn,j,i))
        N1p_s1(k,isn,j,i) = M1_s1(k,isn,j,i)*sigC/r0/rho(j)/dlonm
        N1h_s1(k,isn,j,i) = 0
      enddo
    endif

! calculate N coefficients at (i,j+0.5,k) for S2 points
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlatS2_h)
      do concurrent (k = 1:npts_s2(j))

! Equation (6.20)
! N2P(j+0.5) = M2(j+0.5)*[sigP*d2^2](j+0.5)*sqrt(1-0.75*rho(j+0.5)^2)/R/(rho(j+1)-rho(j))
        N2p_s2(k,isn,j,i) = M2_s2(k,isn,j,i)* &
          sigP_s2(k,isn,j,i)*d2d2_s2(k,isn,j,i)* &
          sqrt(1-3*rho_s(j)**2/4)/r0/(rho(j+1)-rho(j))

! Equation (6.21)
! N2H(j+0.5) = M2(j+0.5)*[sigH*D+sigP*d1*d2](j+0.5)/2/R/rho(j+0.5)/(phi(i+1)-phi(i-1)))
        N2h_s2(k,isn,j,i) = M2_s2(k,isn,j,i)* &
          (sigH_s2(k,isn,j,i)*D_s2(k,isn,j,i)+ &
          sigP_s2(k,isn,j,i)*d1d2_s2(k,isn,j,i))/ &
          4/r0/rho_s(j)/dlonm
      enddo
    enddo

  endsubroutine calculate_n
!-----------------------------------------------------------------------
  pure subroutine calculate_je( &
    mlatd0,mlatd1,mlond0,mlond1,npts_s1,npts_s2, &
    D_s1,be3_s1,d1d1_s1,d1d2_s1,d2d2_s1,sigP_s1,sigH_s1,un_s1,vn_s1, &
    D_s2,be3_s2,d1d2_s2,d2d2_s2,sigP_s2,sigH_s2,un_s2,vn_s2, &
    d1_s1,d2_s1,d1_s2,d2_s2,Je1D_s1,Je2D_s2)
! calculate wind driven currents

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
    real(kind=rp) :: ue1,ue2,je2d

    Je1D_s1 = fill_value
    Je2D_s2 = fill_value

! Equation (3.8)
! calculate Je1D = sigP*d1^2*ue2*Be3+(sigH*D-sigP*d1*d2)*ue1*Be3 at (i+0.5,j,k) for S1 points
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
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

! Equation (3.16) for the equatorial volume at (i+0.5,j=J,k=1)
! Je1D <- Je1D - (sigH*D-sigP*d1*d2)/(sigP*d2^2)*(Je2LB-Je2D)
    j = nmlat_h
    k = 1
    if (j>=mlatd0 .and. j<=mlatd1) then
      do i = mlond0,mlond1
        do isn = 1,2
          ue1 = un_s1(k,isn,j,i)*d1_s1(1,k,isn,j,i)+ &
                vn_s1(k,isn,j,i)*d1_s1(2,k,isn,j,i)
          ue2 = un_s1(k,isn,j,i)*d2_s1(1,k,isn,j,i)+ &
                vn_s1(k,isn,j,i)*d2_s1(2,k,isn,j,i)

! calculate Je2D = (sigH*D+sigP*d1*d2)*ue2*Be3-sigP*d2^2*ue1*Be3 at the equator
! there is no S2 point at the equator therefore it needs to be calculated
          je2d = (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)+ &
            sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))* &
            ue2*be3_s1(k,isn,j,i)- &
            sigP_s1(k,isn,j,i)*d2d2_s1(k,isn,j,i)* &
            ue1*be3_s1(k,isn,j,i)

! Je2LB is the current from lower atmosphere, we assume Je2LB = -J3LB
! not exactly since they are half height level apart but should be close
          Je1D_s1(k,isn,j,i) = Je1D_s1(k,isn,j,i)- &
            (sigH_s1(k,isn,j,i)*D_s1(k,isn,j,i)- &
            sigP_s1(k,isn,j,i)*d1d2_s1(k,isn,j,i))/ &
            (sigP_s1(k,isn,j,i)*d2d2_s1(k,isn,j,i))* &
            (-J3LB(isn,j,i)-je2d)
        enddo
      enddo
    endif

! Equation (3.9)
! calculate Je2D = (sigH*D+sigP*d1*d2)*ue2*Be3-sigP*d2^2*ue1*Be3 at (i,j+0.5,k) for S2 points
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
  function balance_fac_hl(mlatd0,mlatd1,mlond0,mlond1, &
    M3_p,fac_hl_in_p) result(fac_hl_out_p)
! balance the input high latitude FAC to make sure it is zero
! when integrated over the globe

! this is called when direct FAC is read in (read_fac==.true.)
! and needs to be corrected (direct FAC input may be unbalanced)
! M3_p is only the bottom level

    use params_module,only:nmlat_h
    use mpi_module,only:reduce_sum_1d

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      M3_p,fac_hl_in_p
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1) :: fac_hl_out_p

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn
    real(kind=rp) :: downfac
    real(kind=rp),dimension(2) :: subsum,fullsum

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    fac_hl_out_p = 0 ! not fill value

! calculate the sum of upward and downward FACs
    subsum = 0

! exclude halo points to avoid double counting
    do i = mlon0,mlon1
      do j = mlat0,mlat1
        if (j>=2 .and. j<=nmlat_h) then ! no pole
          do isn = 1,2
            if (fac_hl_in_p(isn,j,i) > 0) &
              subsum(1) = subsum(1)+fac_hl_in_p(isn,j,i)*M3_p(isn,j,i)
            if (fac_hl_in_p(isn,j,i) < 0) &
              subsum(2) = subsum(2)-fac_hl_in_p(isn,j,i)*M3_p(isn,j,i)
          enddo
        endif
      enddo
    enddo

    fullsum = reduce_sum_1d(subsum,2,-1)

! if the integrated upward current is stronger than the downward current
! scale the upward current down globally
    if (fullsum(1) > fullsum(2)) then
      downfac = fullsum(2)/fullsum(1)
      do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h)
        if (fac_hl_in_p(isn,j,i) > 0) fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)*downfac
        if (fac_hl_in_p(isn,j,i) < 0) fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)
      enddo

! if the integrated upward current is weaker than the downward current
! scale the downward current down globally
    else
      downfac = fullsum(1)/fullsum(2)
      do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h)
        if (fac_hl_in_p(isn,j,i) < 0) fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)*downfac
        if (fac_hl_in_p(isn,j,i) > 0) fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)
      enddo
    endif

  endfunction balance_fac_hl
!-----------------------------------------------------------------------
  function balance_fac_hl_2(mlatd0,mlatd1,mlond0,mlond1, &
    zigP_p,M3_p,fac_hl_in_p) result(fac_hl_out_p)
! balance the input high latitude FAC to make sure it is zero
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

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn
    real(kind=rp) :: fac
    real(kind=rp),dimension(2) :: sumfac,sumzigP,corr
    real(kind=rp),dimension(4) :: tmp_sub,tmp_full

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    fac_hl_out_p = 0 ! not fill value

! corr = - zigP*abs(Jmr)/sinI * [sum_i^N Jmr*area] / [sum_i^N abs(Jmr)*zigP/sinI*area]
! sinI is ignored by assuming it is close to 1 at high latitudes
    sumfac = 0
    sumzigP = 0

! exclude halo points to avoid double counting
    do i = mlon0,mlon1
      do j = mlat0,mlat1
        if (j>=2 .and. j<=nmlat_h) then ! no pole
          do isn = 1,2
            fac = fac_hl_in_p(isn,j,i)*M3_p(isn,j,i)
            sumfac(isn) = sumfac(isn)+fac
            sumzigP(isn) = sumzigP(isn)+zigP_p(isn,j,i)*abs(fac)
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
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=2 .and. j<=nmlat_h)
      fac_hl_out_p(isn,j,i) = fac_hl_in_p(isn,j,i)- &
        zigP_p(isn,j,i)*abs(fac_hl_in_p(isn,j,i))*corr(isn)
    enddo

  endfunction balance_fac_hl_2
!-----------------------------------------------------------------------
  pure subroutine calculate_ed(mlatd0,mlatd1,mlond0,mlond1, &
    pot_p,ed1_s1,ed2_s1,ed1_s2,ed2_s2)
! calculates electric field Ed1, Ed2 at S1 and S2 points

    use params_module,only:nmlat_h,nmlatS2_h,ylonm,rho,rho_s
    use cons_module,only:r0,fill_value

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_p
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: &
      ed1_s1,ed2_s1,ed1_s2,ed2_s2

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn
    real(kind=rp) :: dlonm

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    ed1_s1 = fill_value
    ed2_s1 = fill_value
    ed1_s2 = fill_value
    ed2_s2 = fill_value

! assume equidistant longitudinal grid points
    dlonm = ylonm(2)-ylonm(1)

! Equations (6.4) and (6.5) (no pole or equator)
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j>=2 .and. j<=nmlat_h-1)
      ed1_s1(isn,j,i) = (pot_p(isn,j,i)-pot_p(isn,j,i+1))/r0/rho(j)/dlonm
      ed2_s1(isn,j,i) = &
        (pot_p(isn,j-1,i)+pot_p(isn,j-1,i+1)- &
         pot_p(isn,j+1,i)-pot_p(isn,j+1,i+1))* &
        sqrt(1-3*rho(j)**2/4)/r0/2/(rho(j+1)-rho(j-1))
    enddo

! Equations (6.8) and (6.9)
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j <= nmlatS2_h)
      ed1_s2(isn,j,i) = &
        (pot_p(isn,j,i-1)+pot_p(isn,j+1,i-1)- &
         pot_p(isn,j,i+1)-pot_p(isn,j+1,i+1))/ &
        4/r0/rho_s(j)/dlonm
      ed2_s2(isn,j,i) = (pot_p(isn,j,i)-pot_p(isn,j+1,i))* &
        sqrt(1-3*rho_s(j)**2/4)/r0/(rho(j+1)-rho(j))
    enddo

  endsubroutine calculate_ed
!-----------------------------------------------------------------------
  pure subroutine calculate_ve( &
    mlatd0,mlatd1,mlond0,mlond1, &
    ed1_s1,ed2_s1,be3_s1, &
    ed1_s2,ed2_s2,be3_s2, &
    ve1_s1,ve2_s1,ve1_s2,ve2_s2)
! calculates drift velocity ve1, ve2 at S1 and S2 points
! ve1 =  Ed2/Be3
! ve2 = -Ed1/Be3
! be3 is only the bottom level

    use params_module,only:nmlat_h,nmlatS2_h
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

! no pole or equator
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j>=2 .and. j<=nmlat_h-1)
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

! no pole or equator
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j>=2 .and. j<=nmlat_h-1)
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

! no pole or equator
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j>=2 .and. j<=nmlat_h-1)
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
    npts_s1,npts_s2,npts_r,pot_p, &
    M1_s1,N1p_s1,N1h_s1,Je1D_s1, &
    M2_s2,N2p_s2,N2h_s2,Je2D_s2, &
    M3_r,I1_s1,I2_s2,I3_r)

    use params_module,only:nhgt_fix,nhgt_fix_r,nmlat_h,nmlatS2_h,nmlon
    use cons_module,only:J3LB,fill_value
    use mpi_module,only:gather_mlon_3d,sync_mlat_5d,sync_mlon_5d

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_s1,npts_r
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_p
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      M1_s1,N1p_s1,N1h_s1,Je1D_s1,M2_s2,N2p_s2,N2h_s2,Je2D_s2
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: M3_r
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: I1_s1,I2_s2
    real(kind=rp),dimension(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: I3_r

    integer :: mlat0,mlat1,mlon0,mlon1,isn,i,j,k,iconj
    real(kind=rp),dimension(nhgt_fix,2,mlond0:mlond1) :: Je1_sub
    real(kind=rp),dimension(nhgt_fix,2,nmlon) :: Je1_full
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1) :: &
      I1_1,I1_2,I1_3,I2_1,I2_2,I2_3
    real(kind=rp),dimension(2,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1) :: tmpI

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1
    I1_s1 = fill_value
    I2_s2 = fill_value
    I3_r = fill_value

! Equation (6.10)
! Equations (6.16) to (6.18) for the equatorial volume at (i+0.5,j=J,k=1)
! no need for separate cases since they have already been considered in N1P and Je1D
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j >= 2) ! no pole
      do concurrent (k = 1:npts_s1(j))
        I1_1(k,isn,j,i) = N1p_s1(k,isn,j,i)*(pot_p(isn,j,i)-pot_p(isn,j,i+1))
        if (j == nmlat_h) then
          I1_2(k,isn,j,i) = 0
        else
          I1_2(k,isn,j,i) = -N1h_s1(k,isn,j,i)* &
            (pot_p(isn,j-1,i)+pot_p(isn,j-1,i+1)- &
             pot_p(isn,j+1,i)-pot_p(isn,j+1,i+1))
        endif
        I1_3(k,isn,j,i) = M1_s1(k,isn,j,i)*Je1D_s1(k,isn,j,i)
        I1_s1(k,isn,j,i) = I1_1(k,isn,j,i)+I1_2(k,isn,j,i)+I1_3(k,isn,j,i)
      enddo
    enddo

! pole value see Equation (6.15)
! Je1(i+0.5,j=1,k) = 0.5*[Je1(i+0.5,j=2,k) - Je1(i'+0.5,j=2,k)]
! i' is the conjugate longitude
! I1(i+0.5,j=1,k) = Je1(i+0.5,j=1,k)*M1(i+0.5,j=1,k)
! we do not have i' in the current process, so first gather all longitudes
    j = mlat0
    Je1_sub = fill_value
    do concurrent (i = mlon0:mlon1, isn = 1:2, k = 1:npts_s1(j))
      Je1_sub(k,isn,i) = I1_s1(k,isn,j+1,i)/M1_s1(k,isn,j+1,i) ! Je1(i+0.5,2,k)
    enddo
    Je1_full = gather_mlon_3d(Je1_sub(:,:,mlon0:mlon1),nhgt_fix,2)
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

! Equation (6.19)
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, j <= nmlatS2_h)
      do concurrent (k = 1:npts_s2(j))
        I2_1(k,isn,j,i) = N2p_s2(k,isn,j,i)*(pot_p(isn,j,i)-pot_p(isn,j+1,i))
        I2_2(k,isn,j,i) = N2h_s2(k,isn,j,i)* &
          (pot_p(isn,j,i-1)+pot_p(isn,j+1,i-1)- &
           pot_p(isn,j,i+1)-pot_p(isn,j+1,i+1))
        I2_3(k,isn,j,i) = M2_s2(k,isn,j,i)*Je2D_s2(k,isn,j,i)
        I2_s2(k,isn,j,i) = I2_1(k,isn,j,i)+I2_2(k,isn,j,i)+I2_3(k,isn,j,i)
      enddo
    enddo

! I1 at i-1 and I2 at j-1 are used in I3 calculation, so sync before used
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, k = 1:nhgt_fix)
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

! lower boundary is the current from lower atmosphere J3LB*M3
      k = 1
      I3_r(k,isn,j,i) = J3LB(isn,j,i)*M3_r(k,isn,j,i)

! Equation (6.27)
      do k = 2,npts_r(j)
        I3_r(k,isn,j,i) = I3_r(k-1,isn,j,i)+ &
          I1_s1(k-1,isn,j,i-1)-I1_s1(k-1,isn,j,i)- &
          I2_s2(k-1,isn,j,i)
        if (j /= 1) I3_r(k,isn,j,i) = I3_r(k,isn,j,i)+I2_s2(k-1,isn,j-1,i)
      enddo
    enddo

  endsubroutine calculate_current
!-----------------------------------------------------------------------
endmodule calculate_terms_module
