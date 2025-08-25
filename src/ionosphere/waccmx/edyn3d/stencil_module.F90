module stencil_module

  use prec,only:rp

  implicit none

! coefficients are calculated at P points
! coef ordering
! coef(4) (i-1,j+1)  coef(3) (i,j+1)  coef(2) (i+1,j+1)
! coef(5) (i-1,j  )  coef(9) (i,j  )  coef(1) (i+1,j  )
! coef(6) (i-1,j-1)  coef(7) (i,j-1)  coef(8) (i+1,j-1)

  contains
!-----------------------------------------------------------------------
  pure function calculate_coef3d(mlatd0,mlatd1,mlond0,mlond1, &
    npts_p,npts_s2,N1p_s1,N1h_s1,N2p_s2,N2h_s2) result(coef3d)
! calculate height-dependent matrix coefficients

    use params_module,only:nhgt_fix,nmlat_h,nmlatS2_h

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_p
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      N1p_s1,N1h_s1,N2p_s2,N2h_s2
    real(kind=rp),dimension(9,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1) :: coef3d

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,k
    real(kind=rp) :: N2p_p,N2h_p,N2h_pp,N2h_pm

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    coef3d = 0

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)
      do concurrent (k = 1:npts_p(j))
        N2p_p = 0
        N2h_p = 0
        N2h_pp = 0
        N2h_pm = 0

! S2 grid does not have this point
        if (j <= nmlatS2_h) then
          if (k <= npts_s2(j)) then
            N2p_p = N2p_s2(k,isn,j,i)
            N2h_p = N2h_s2(k,isn,j,i)
            N2h_pp = N2h_s2(k,isn,j,i+1)
            N2h_pm = N2h_s2(k,isn,j,i-1)
          endif
        endif

! the volume at the south pole is a combined cell of all longitudes (a prism with nmlon edges)
! Equation (6.50) will be combined with the b matrix later to give Equation (7.23)
! Equation (7.23) is the actual equation to be solved in matrix form

! summation over longitudes is not done here (no gathering)
! it will be carried out when constructing the matrix
        if (j == 1) then

! Equation (6.51)
          coef3d(3,k,isn,j,i) = N2p_p-N2h_pp+N2h_pm

! Equation (6.52)
          coef3d(9,k,isn,j,i) = -N2p_p

        else

! Equations (6.33) to (6.42)
          coef3d(1,k,isn,j,i) = &
              N1p_s1(k,isn,j  ,i  ) &
            - N2h_s2(k,isn,j-1,i  ) + N2h_p
          coef3d(2,k,isn,j,i) = &
            - N1h_s1(k,isn,j  ,i  ) &
            + N2h_p
          coef3d(3,k,isn,j,i) = &
              N1h_s1(k,isn,j  ,i-1) - N1h_s1(k,isn,j  ,i  ) &
            + N2p_p
          coef3d(4,k,isn,j,i) = &
              N1h_s1(k,isn,j  ,i-1) &
            - N2h_p
          coef3d(5,k,isn,j,i) = &
              N1p_s1(k,isn,j  ,i-1) &
            + N2h_s2(k,isn,j-1,i  ) - N2h_p
          coef3d(6,k,isn,j,i) = &
            - N1h_s1(k,isn,j  ,i-1) &
            + N2h_s2(k,isn,j-1,i  )
          coef3d(7,k,isn,j,i) = &
            - N1h_s1(k,isn,j  ,i-1) + N1h_s1(k,isn,j  ,i  ) &
            + N2p_s2(k,isn,j-1,i  )
          coef3d(8,k,isn,j,i) = &
              N1h_s1(k,isn,j  ,i  ) &
            - N2h_s2(k,isn,j-1,i  )
          coef3d(9,k,isn,j,i) = &
            - N1p_s1(k,isn,j  ,i-1) - N1p_s1(k,isn,j  ,i  ) &
            - N2p_s2(k,isn,j-1,i  ) - N2p_p
        endif
      enddo
    enddo

  endfunction calculate_coef3d
!-----------------------------------------------------------------------
  pure function calculate_coef2d(mlatd0,mlatd1,mlond0,mlond1,coef3d) result(coef2d)
! add the coefficients in height to get the coefficients for each hemisphere

    use params_module,only:nhgt_fix

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(9,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef3d
    real(kind=rp),dimension(9,2,mlatd0:mlatd1,mlond0:mlond1) :: coef2d

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,ic

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

! include poles altogether (north pole is not used)
! only ic=3,9 are actually used in south pole, all other coefficients are zero
! coef3d is zero above the topmost level at low latitudes (no contribution to summation)
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2, ic = 1:9)
      coef2d(ic,isn,j,i) = sum(coef3d(ic,:,isn,j,i))
    enddo

  endfunction calculate_coef2d
!-----------------------------------------------------------------------
  pure function calculate_src3d(mlatd0,mlatd1,mlond0,mlond1, &
    npts_p,npts_s2,M1_s1,Je1D_s1,M2_s2,Je2D_s2,M3_r) result(src3d)
! calculate wind driven ionospheric current sources

! M3_r is only the bottom level

    use params_module,only:nhgt_fix,nmlat_h,nmlatS2_h
    use cons_module,only:J3LB

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_p
    integer,dimension(nmlatS2_h),intent(in) :: npts_s2
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      M1_s1,Je1D_s1,M2_s2,Je2D_s2
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: M3_r
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1) :: src3d

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,k
    real(kind=rp) :: je2d_p

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    src3d = 0

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)
      do concurrent (k = 1:npts_p(j))
        je2d_p = 0

! S2 grid does not have this point
        if (j <= nmlatS2_h) then
          if (k <= npts_s2(j)) je2d_p = M2_s2(k,isn,j,i)*Je2D_s2(k,isn,j,i)
        endif

! Equation (6.50)
! S(i,1,k) = M2(i,1.5,k)*Je2D(i,1.5,k)
! north pole values are not really used but south pole
        if (j == 1) then
          src3d(k,isn,j,i) = -je2d_p

! Equation (6.42)
! S(i,j,k) = M1(i-0.5,j,k)*Je1D(i-0.5,j,k)-M1(i+0.5,j,k)*Je1D(i+0.5,j,k)+
!            M2(i,j-0.5,k)*Je2D(i,j-0.5,k)-M2(i,j+0.5,k)*Je2D(i,j+0.5,k)
        else
          src3d(k,isn,j,i) = &
            M1_s1(k,isn,j  ,i-1)*Je1D_s1(k,isn,j  ,i-1)- &
            M1_s1(k,isn,j  ,i  )*Je1D_s1(k,isn,j  ,i  )+ &
            M2_s2(k,isn,j-1,i  )*Je2D_s2(k,isn,j-1,i  )- &
            je2d_p
        endif
      enddo

! add the current from lower atmosphere J3LB*M3
      src3d(1,isn,j,i) = src3d(1,isn,j,i)+J3LB(isn,j,i)*M3_r(isn,j,i)
    enddo

  endfunction calculate_src3d
!-----------------------------------------------------------------------
  pure function calculate_src2d(mlatd0,mlatd1,mlond0,mlond1,src3d) result(src2d)
! add the source in height to get the source for each hemisphere

    use params_module,only:nhgt_fix

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: src3d
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1) :: src2d

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)
      src2d(isn,j,i) = sum(src3d(:,isn,j,i))
    enddo

  endfunction calculate_src2d
!-----------------------------------------------------------------------
  pure function calculate_bij(mlatd0,mlatd1,mlond0,mlond1,coef2d) result(bij)
! set field-aligned conductance (b) matrix

    use params_module,only:rho,rho_s
    use cons_module,only:dtr,jlatm_JT

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(9,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef2d
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1) :: bij

    real(kind=rp),parameter :: &

! b_mult is [|Phi|/Delta(Phi)]*(R/L)^2, where Phi is a characteristic potential value,
! Delta(Phi) is a characteristic allowed interhemispheric potential difference,
! R is Earth radius, and L is a characteristic N-S length scale for Phi.
! It is assumed that b_mult is similar for middle and auroral latitudes.
      b_mult = 1e3_rp, &

! pccolat is the polar cap colatitude, which for now is fixed.
! But it can be made variable w.r.t. time and magnetic longitude in the future.
      pccolat = 14, & ! we can move this further equatorward
      rho_pc = sin(pccolat*dtr)

    integer :: mlat0,mlat1,mlon0,mlon1,i,j
    real(kind=rp) :: fac3

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

! initialize b to zero
! this also sets b in low latitudes where fieldlines are assumed to be equipotential
! and it should not be used there (zero would actually mean two hemispheres are uncoupled)
    bij = 0

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, j>=2 .and. j<=jlatm_JT-1)
      bij(j,i) = b_mult*(rho_s(j)-rho_s(j-1))**2/ &
        (1/(coef2d(3,1,j,i)+coef2d(7,1,j,i))+ &
         1/(coef2d(3,2,j,i)+coef2d(7,2,j,i)))
    enddo

! set b to zero within polar caps, transitioning linearly
! to the full original value over a distance of about (1/3) pccolat
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, j <= jlatm_JT-1)
      fac3 = 3*(rho(j)/rho_pc-1)
      if (fac3 <= 0) bij(j,i) = 0
      if (fac3>0 .and. fac3<1) bij(j,i) = fac3*bij(j,i)

! if fac3>=1, b remains unmodified
    enddo

  endfunction calculate_bij
!-----------------------------------------------------------------------
endmodule stencil_module
