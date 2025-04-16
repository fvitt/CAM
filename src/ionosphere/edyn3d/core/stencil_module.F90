module stencil_module

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  pure function calculate_coef(mlatd0,mlatd1,mlond0,mlond1,npts_p, &
    S_p,N1p_s1,N1h_s1,N2p_s2,N2h_s2) result(coef)
! calculate height-dependent matrix coefficients

    use params_module,only:nhgt_fix,nmlat_h

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    integer,dimension(nmlat_h),intent(in) :: npts_p
    real(kind=rp),dimension(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: &
      S_p,N1p_s1,N1h_s1,N2p_s2,N2h_s2
    real(kind=rp),dimension(10,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1) :: coef

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,k
    real(kind=rp) :: N2p_p,N2h_p

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

! coef ordering from TIEGCM
! ^   equatorward
! coef(4) (i-1,j+1)      coef(3) (i,j+1)    coef(2) (i+1,j+1)
! coef(5) (i-1,j)        coef(9) (i,j)      coef(1) (i+1,j)
! coef(6) (i-1,j-1)      coef(7) (i,j-1)    coef(8) (i+1,j-1)
! v   poleward

! relationship between P,S1,S2 points for the same index (i,j):
! P(i,j) is S1(i+0.5,j) and S2(i,j+0.5) with j increasing equatorward
! coefficients are calculated at P points

    coef = 0._rp

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, isn = 1:2)
      do concurrent (k = 1:npts_p(j))

! south pole Eq 338 in Art 2023/01/17 notes
! Sum_i C3(i,1) * Phi(i,j+1) + C9P * Phi^P + beta * Phi^P* + I3^TP = 0
! note for I3^TP the sign will be changed since in the code it is on the RHS

! C3(i,1) * Phi(i,j+1)
!   -> Phi(i,j+1) * Sum_k=1^K [ N2P(i,3/2,k)-N2H(i+1,3/2,k)+N2H(i-1,3/2,k) ]
!   -> coef(3) south pole
! -Sum_k=1^K [ Sum_i=1^nmlon N2P(i,3/2,k) ] - Sum_i^nmlon b(i,1)
!   -> coef(9) south pole
! -Sum_k=1^K [ Sum_i=1^nmlon S(i,1,k) ] + [ Sum_i^nmlon b(i,1) ] Phi^NP
!   -> coef(10) south pole
! beta = Sum_i b(i,1) -> precaluclated
! I3^TP = Sum_i [ Sum_k^K S(i,1,k) - I3^R(i,1) ] + I3^P1/2
!   with I3^P(1/2) = Sum_i I3(i,1,1/2)
!        I3^R(i,1) = M3(i,1,K+1/2) * Jr^R(i,1)
!          -> upper boundary is calculated later in calculate_fac_hl
!        S(i,1,k)  = -M2(i,3/2,k) * Je2^D(i,3/2,k)
!          -> calculated in calculate_s and already included the minus sign on the RHS
!             includes I3^P(1/2) at S(i,1,1)
        if (j == 1) then

! N2P(i,3/2,k)-N2H(i+1,3/2,k)+N2H(i-1,3/2,k)
          coef(3,k,isn,j,i) = &
              N2p_s2(k,isn,j,i  ) &
            - N2h_s2(k,isn,j,i+1) + N2h_s2(k,isn,j,i-1)

! -N2P(i,3/2,k)
          coef(9,k,isn,j,i) = -N2p_s2(k,isn,j,i)

          coef(10,k,isn,j,i) = S_p(k,isn,j,i)

        else
          if (k == nmlat_h-j+1) then ! top volume at equator
            N2p_p = 0._rp
            N2h_p = 0._rp
          else ! i,j+0.5
            N2p_p = N2p_s2(k,isn,j,i)
            N2h_p = N2h_s2(k,isn,j,i)
          endif

          coef(1,k,isn,j,i) = &
              N1p_s1(k,isn,j  ,i  ) &
            - N2h_s2(k,isn,j-1,i  ) + N2h_p
          coef(2,k,isn,j,i) = &
            - N1h_s1(k,isn,j  ,i  ) &
            + N2h_p
          coef(3,k,isn,j,i) = &
              N1h_s1(k,isn,j  ,i-1) - N1h_s1(k,isn,j  ,i  ) &
            + N2p_p
          coef(4,k,isn,j,i) = &
              N1h_s1(k,isn,j  ,i-1) &
            - N2h_p
          coef(5,k,isn,j,i) = &
              N1p_s1(k,isn,j  ,i-1) &
            + N2h_s2(k,isn,j-1,i  ) - N2h_p
          coef(6,k,isn,j,i) = &
            - N1h_s1(k,isn,j  ,i-1) &
            + N2h_s2(k,isn,j-1,i  )
          coef(7,k,isn,j,i) = &
            - N1h_s1(k,isn,j  ,i-1) + N1h_s1(k,isn,j  ,i  ) &
            + N2p_s2(k,isn,j-1,i  )
          coef(8,k,isn,j,i) = &
              N1h_s1(k,isn,j  ,i  ) &
            - N2h_s2(k,isn,j-1,i  )
          coef(9,k,isn,j,i) = &
            - N1p_s1(k,isn,j  ,i-1) - N1p_s1(k,isn,j  ,i  ) &
            - N2p_s2(k,isn,j-1,i  ) - N2p_p

          coef(10,k,isn,j,i) = S_p(k,isn,j,i)
        endif
      enddo
    enddo

  endfunction calculate_coef
!-----------------------------------------------------------------------
  pure function calculate_coef_ns2(mlatd0,mlatd1,mlond0,mlond1,coef) result(coef_ns2)
! add the coefficients in height to get coefficients for each hemisphere
! needed for calculating high latitude FAC

    use params_module,only:nhgt_fix
    use cons_module,only:phi_pol

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(10,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef
    real(kind=rp),dimension(10,2,mlatd0:mlatd1,mlond0:mlond1) :: coef_ns2

    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,k,ic

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    coef_ns2 = 0._rp

    do i = mlon0,mlon1
      do j = mlat0,mlat1

! Phi^SP(i=1,j=1) = Phi(i,j=1) for i = 2,nmlon
!   -> C9N(i,1) = 1 C*(1,1) = -1 C1S-C8S = 0 C10S = 0
! in NH Phi(i,1) = Phi^NP
!   -> C9N(i,1) = 1 C10N(i,1) = Phi^NP C1N-C8N = 0
        if (j == 1) then

! A. Maute 2023/01: south pole
          isn = 1
          do k = 1,nhgt_fix

! need to move C3(i,j) to the appropriate place on the LHS
            coef_ns2(3,isn,j,i) = coef_ns2(3,isn,j,i)+coef(3,k,isn,j,i)

! -Sum_k=1^K N2P(i,3/2,k) -> put into coef_ns2(i,j,isn,9)
            coef_ns2(9,isn,j,i) = coef_ns2(9,isn,j,i)+coef(9,k,isn,j,i)

! Sum_k=1^K S(i,1,k) -> put into coef_ns2(i,j,isn,1-)
            coef_ns2(10,isn,j,i) = coef_ns2(10,isn,j,i)+coef(10,k,isn,j,i)
          enddo

! A. Richmond 2023/06/20: north pole
          isn = 2
          coef_ns2(10,isn,j,i) = phi_pol ! north pole for each i Phi^N(i,1) = Phi^NP
          coef_ns2(9,isn,j,i) = 1
          do concurrent (ic = 1:8)
            coef_ns2(ic,isn,j,i) = 0._rp
          enddo

! no pole (done above)
        else
          do isn = 1,2
            do k = 1,nhgt_fix
              do ic = 1,10
                coef_ns2(ic,isn,j,i) = coef_ns2(ic,isn,j,i)+coef(ic,k,isn,j,i)
              enddo
            enddo
          enddo
        endif
      enddo
    enddo

  endfunction calculate_coef_ns2
!-----------------------------------------------------------------------
  pure function calculate_coef_ns(mlatd0,mlatd1,mlond0,mlond1,coef_ns2) result(coef_ns)
! set the coefficient matrix in both hemispheres
! decide where to add the SH and NH stencil and where not
! change the direction of NH stencil from coef_ns2 to coef_ns
! LHS+RHS for each P-point
! A. Maute 2023/02: solve two hemispheres

    use params_module,only:nmlat_h
    use cons_module,only:jlatm_JT

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(10,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_ns2
    real(kind=rp),dimension(10,2,mlatd0:mlatd1,mlond0:mlond1) :: coef_ns

    integer :: i,j,isn,ic
    real(kind=rp) :: s

! from pole to latm_JT, set coefficients separately in two hemispheres
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, ic = 1:10, j>=1 .and. j<=jlatm_JT-1)
      coef_ns(ic,isn,j,i) = coef_ns2(ic,isn,j,i)
    enddo

    j = jlatm_JT
    if (j>=mlatd0 .and. j<=mlatd1) then
      do concurrent (i = mlond0:mlond1, ic = 1:10)
        if (ic>=6 .and. ic<=8) then ! don't add values from the other hemisphere
          do concurrent (isn = 1:2)
            coef_ns(ic,isn,j,i) = coef_ns2(ic,isn,j,i)
          enddo
        else ! add values from both hemispheres
          s = coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i)
          coef_ns(ic,1,j,i) = s
          coef_ns(ic,2,j,i) = s
        endif
      enddo
    endif

! from latm_JT to equator, add values from both hemispheres
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, ic = 1:10, j>=jlatm_JT+1 .and. j<=nmlat_h)
      s = coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i)
      coef_ns(ic,1,j,i) = s
      coef_ns(ic,2,j,i) = s
    enddo

! set equatorial boundary condition (page 14 Art's notes)
! there should be just one equator value
    j = nmlat_h
    if (j>=mlatd0 .and. j<=mlatd1) then
      do i = mlond0,mlond1
        do isn = 1,2
          do ic = 2,4
            s = (coef_ns(ic,isn,j,i)+coef_ns(10-ic,isn,j,i))/2
            coef_ns(ic,isn,j,i) = s
            coef_ns(10-ic,isn,j,i) = s
          enddo
        enddo
      enddo
    endif

  endfunction calculate_coef_ns
!-----------------------------------------------------------------------
  pure function calculate_bij(mlatd0,mlatd1,mlond0,mlond1,coef_ns2) result(bij)
! set field-aligned conductance (b) matrix

    use params_module,only:rho,rho_s
    use cons_module,only:jlatm_JT

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(10,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_ns2
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1) :: bij

    real(kind=rp),parameter :: &

! b_mult is [|Phi|/Delta(Phi)]*(R/L)^2, where Phi is a characteristic potential value,
! Delta(Phi) is a characteristic allowed interhemispheric potential difference,
! R is Earth radius, and L is a characteristic N-S length scale for Phi.
! It is assumed that b_mult is similar for middle and auroral latitudes.
      b_mult = 1e3_rp, &

! pccolatrad is the polar cap colatitude in radians, which for now is fixed.
! But it can be made variable w.r.t. time and magnetic longitude in the future.
      pccolatrad = 0.25_rp, & ! 14 degree
      rho_pc = sin(pccolatrad)

    integer :: mlat0,mlat1,mlon0,mlon1,i,j
    real(kind=rp) :: fac3

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

! initialize bij to zero
! this also sets bij in low latitudes where fieldlines are assumed to be equipotential
! and it should not be used there (zero would actually mean two hemispheres are uncoupled)
    bij = 0

    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, j>=2 .and. j<=jlatm_JT)
      bij(j,i) = b_mult*(rho_s(j)-rho_s(j-1))**2/ &
        (1/(coef_ns2(3,1,j,i)+coef_ns2(7,1,j,i))+ &
         1/(coef_ns2(3,2,j,i)+coef_ns2(7,2,j,i)))
    enddo

! set bij to zero within polar caps, transitioning linearly
! to the full original value over a distance of about (1/3) pccolatrad
    do concurrent (i = mlon0:mlon1, j = mlat0:mlat1, j <= jlatm_JT)
      fac3 = 3*(rho(j)/rho_pc-1)

      if (fac3 <= 0) bij(j,i) = 0

! the same at conjugate points since bij is the same
      if (fac3>0 .and. fac3<1) bij(j,i) = fac3*bij(j,i)

! if fac3>=1 bij remains unmodified
    enddo

  endfunction calculate_bij
!-----------------------------------------------------------------------
endmodule stencil_module
