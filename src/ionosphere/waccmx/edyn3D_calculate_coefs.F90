   module edyn3D_calculate_coefs
     !
     !  This module calculates coefficients and right hand side and left hand side
     !  to solve for potential.  Based on standalone module calc_coef with gather and
     !  scatter added
     !
     use edyn3D_params,  only: nmlon,nmlat_h,nhgt_fix,nlonlat,nmlatS2_h,nmlat_T1
     use edyn3d_mpi,     only: mlon0_p,mlon1_p,mp_poten_halos_edyn3D
     use shr_kind_mod,   only: r8 => shr_kind_r8            ! 8-byte reals
     use cam_logfile,    only: iulog
     use spmd_utils,     only: masterproc
     use cam_abortutils, only: endrun

     implicit none
     private

     public :: edyn3D_calculate_coef,edyn3D_calculate_coef_ns2,edyn3D_calculate_coef_ns,edyn3D_calculate_bij

!     real(r8), allocatable :: coef(:,:,:,:,:)       ! for each P-point each hemisphere and height
!     real(r8), allocatable :: coef_ns(:,:,:)        ! lhs+rhs: for each P-point
!     real(r8), allocatable :: coef_ns_glb(:,:,:)    ! lhs+rhs:  globally
!     real(r8), allocatable :: coef_ns2(:,:,:,:)     ! lhs+rhs: for each P-point
!     real(r8), allocatable :: bij(:,:)              ! Field aligned conductance
!     real(r8), allocatable :: rhs_ns(:)             ! forcing: for each P-point
!     real(r8), allocatable :: lhs_ns(:,:)           ! lhs: for each P-point
!     real(r8), allocatable :: poten_glb(:,:,:)      ! Potential: globally

     contains

!-----------------------------------------------------------------------------
     subroutine edyn3D_calculate_coef(fline_p,fline_s1,fline_s2,coef)
     !
     ! Calculate matrix coefficients for left hand and right hand side inputs to global potential solver
     !
     use edyn3D_fieldline, only: fieldline_p,fieldline_s1,fieldline_s2

     implicit none

     type(fieldline_p),dimension(mlon0_p-1:mlon1_p+1,nmlat_h,2),intent(in) :: fline_p
     type(fieldline_s1),dimension(mlon0_p-1:mlon1_p+1,nmlat_h,2),intent(in) :: fline_s1
     type(fieldline_s2),dimension(mlon0_p-1:mlon1_p+1,nmlatS2_h,2),intent(in) :: fline_s2
     real(r8),dimension(mlon0_p:mlon1_p,nmlat_h,nhgt_fix,10,2),intent(out) :: coef

     !
     ! coef ordering from TIEGCM
     ! ^   equatorward
     ! coef(4) (i-1,j+1)      coef(3) (i,j+1)	 coef(2) (i+1,j+1)
     ! coef(5) (i-1,j)        coef(9) (i,j)	 coef(1) (i+1,j)
     ! coef(6) (i-1,j-1)      coef(7) (i,j-1)	 coef(8) (i+1,j-1)
     ! v   poleward
     !
     ! relationship between P,S1, and S2 point for the same index (i,j)
     !  P(i,j) then is really S1(i+0.5,j) and S2(i,j+0.5) with j increasing equatorward
     !  coefficient is calculated at P points
     !

     integer :: isn,i,j,k,im,nmax,status,ic
     real(r8) :: N2p_p,N2h_p

     coef = 0._r8

     !
     ! Calculate initial coefficients from s1 and s2 field line variables
     !
     do i=mlon0_p,mlon1_p ! loop over task longitudes

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
     !  	-> upper boundary is calculated later in calculate_fac_hl
     !        S(i,1,k)  = -M2(i,3/2,k) * Je2^D(i,3/2,k)
     !  	-> calculated in calculate_s and already included the minus sign on the RHS
     !  	   includes I3^P(1/2) at S(i,1,1)

       j = 1
       do isn = 1,2
 	 do k = 1,fline_p(i,j,isn)%npts

 	   ! N2P(i,3/2,k)-N2H(i+1,3/2,k)+N2H(i-1,3/2,k)
 	   coef(i,j,k,3,isn) = &
 	       fline_s2(i,j,isn  )%N2p(k) &
 	     - fline_s2(i+1,j,isn)%N2h(k) + fline_s2(i-1,j,isn)%N2h(k)

 	   ! -N2P(i,3/2,k)
 	   coef(i,j,k,9,isn) = -fline_s2(i,j,isn)%N2p(k)

 	   coef(i,j,k,10,isn) = fline_p(i,j,isn)%S(k)
 	 enddo
       enddo

       do isn = 1,2
	 do j = 2,nmlat_h

  	   do k = 1,fline_p(i,j,isn)%npts
	     if (k == nmlat_h-j+1) then ! top volume at equator
	       N2p_p = 0
	       N2h_p = 0
	     else ! i,j+0.5
!	       N2p_p = fline_s2(i,j-1,isn)%N2p(k)
!	       N2h_p = fline_s2(i,j-1,isn)%N2h(k)
	       N2p_p = fline_s2(i,j,isn)%N2p(k)
	       N2h_p = fline_s2(i,j,isn)%N2h(k)
	     endif

	     coef(i,j,k,1,isn) = &
		 fline_s1(i,j,isn)%N1p(k) &
	       - fline_s2(i,j-1,isn)%N2h(k) + N2h_p
	     coef(i,j,k,2,isn) = &
	       - fline_s1(i,j,isn)%N1h(k) &
	       + N2h_p
	     coef(i,j,k,3,isn) = &
		 fline_s1(i-1,j,isn)%N1h(k) - fline_s1(i,j,isn)%N1h(k) &
	       + N2p_p
	     coef(i,j,k,4,isn) = &
		 fline_s1(i-1,j,isn)%N1h(k) &
	       - N2h_p
	     coef(i,j,k,5,isn) = &
		 fline_s1(i-1,j,isn)%N1p(k) &
	       + fline_s2(i,j-1,isn)%N2h(k) - N2h_p
	     coef(i,j,k,6,isn) = &
	       - fline_s1(i-1,j,isn)%N1h(k) &
	       + fline_s2(i,j-1,isn)%N2h(k)
	     coef(i,j,k,7,isn) = &
	       - fline_s1(i-1,j,isn)%N1h(k) + fline_s1(i,j,isn)%N1h(k) &
	       + fline_s2(i,j-1,isn)%N2p(k)
	     coef(i,j,k,8,isn) = &
		 fline_s1(i,j,isn)%N1h(k) &
	       - fline_s2(i,j-1,isn)%N2h(k)
	     coef(i,j,k,9,isn) = &
	       - fline_s1(i-1,j,isn)%N1p(k) - fline_s1(i,j,isn)%N1p(k) &
	       - fline_s2(i,j-1,isn)%N2p(k) - N2p_p

	     coef(i,j,k,10,isn) = fline_p(i,j,isn)%S(k)
	   enddo ! end lat/fieldline loop
	 enddo ! end height loop
       enddo ! end hemisphere loop
     enddo ! end longitude loop

     end subroutine edyn3d_calculate_coef
!-----------------------------------------------------------------------
     subroutine edyn3D_calculate_coef_ns2(coef,coef_ns2)

     ! add the coefficients in height to get coefficients for each hemisphere
     ! needed for calculating high latitude FAC

     use edyn3D_params, only: phi_pol

     real(r8),dimension(mlon0_p:mlon1_p,nmlat_h,nhgt_fix,10,2),intent(in) :: coef
     real(r8),dimension(mlon0_p:mlon1_p,nmlat_h,2,10),intent(out) :: coef_ns2

     integer :: i,j,isn,k,ic

     coef_ns2(:,:,:,:) = 0._r8

     ! Phi^SP(i=1,j=1) = Phi(i,j=1) for i = 2,nmlon
     !   -> C9N(i,1) = 1 C*(1,1) = -1 C1S-C8S = 0 C10S = 0
     ! in NH Phi(i,1) = Phi^NP
     !   -> C9N(i,1) = 1 C10N(i,1) = Phi^NP C1N-C8N = 0

     do i = mlon0_p,mlon1_p

   ! A. Maute 2023/01: south pole

       j = 1
       isn = 1

       do k = 1,nhgt_fix

   ! need to move C3(i,j) to the appropriate place on the LHS
   	 coef_ns2(i,j,isn,3) = coef_ns2(i,j,isn,3)+coef(i,j,k,3,isn)

   ! -Sum_k=1^K N2P(i,3/2,k) -> put into coef_ns2(i,j,isn,9)
   	 coef_ns2(i,j,isn,9) = coef_ns2(i,j,isn,9)+coef(i,j,k,9,isn)

   ! Sum_k=1^K S(i,1,k) -> put into coef_ns2(i,j,isn,1-)
   	 coef_ns2(i,j,isn,10) = coef_ns2(i,j,isn,10)+coef(i,j,k,10,isn)

       enddo

   ! A. Richmond 2023/06/20: north pole

       j = 1
       isn = 2
       coef_ns2(i,j,isn,10) = phi_pol ! north pole for each i Phi^N(i,1) = Phi^NP
       coef_ns2(i,j,isn,9) = 1._r8
       do ic = 1,8
 	 coef_ns2(i,j,isn,ic) = 0._r8
       enddo

       do j = 2,nmlat_h ! no pole
 	 do isn = 1,2
 	   do k = 1,nhgt_fix
 	     do ic = 1,10
 	       coef_ns2(i,j,isn,ic) = coef_ns2(i,j,isn,ic)+coef(i,j,k,ic,isn)
 	     enddo
 	   enddo
 	 enddo
       enddo
     enddo

     end subroutine edyn3D_calculate_coef_ns2

!-----------------------------------------------------------------------------
     subroutine edyn3D_calculate_coef_ns(coef_ns2,coef_ns)
   ! set the the coefficient matrix in both hemispheres
   ! decide where to add the SH and NH stencil and where not
   ! change the direction of NH stencil from coef_ns2 to coef_ns
   ! LHS+RHS for each P-point
   ! A. Maute 2023/02: solve two hemispheres

       use edyn3D_params, only:nmlat_h,nmlat_T1,jlatm_JT,ylatm,ylatm_JT
       use edyn3D_mpi,only:mlon0_p,mlon1_p

       real(r8),dimension(mlon0_p:mlon1_p,nmlat_h,2,10),intent(in) :: coef_ns2
       real(r8),dimension(mlon0_p:mlon1_p,nmlat_T1,10),intent(out) :: coef_ns

       integer :: i,j,ic, &
   	 jS,jN ! overall index from pole to equator

       do i = mlon0_p,mlon1_p

   ! from pole to latm_JT, set coefficients separately in two hemispheres
   	 do j = 1,jlatm_JT-1
   	   jS = j
   	   jN = nmlat_T1-j+1

   	   do ic = 1,10
   	     coef_ns(i,jS,ic) = coef_ns2(i,j,1,ic)
   	     coef_ns(i,jN,ic) = coef_ns2(i,j,2,ic)
   	   enddo

   	 enddo

   	 j = jlatm_JT
   	 jS = j
   	 jN = nmlat_T1-j+1

   ! add values from both hemispheres
   	 do ic = 1,5
   	   coef_ns(i,jS,ic) = coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic)
   	   coef_ns(i,jN,ic) = coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic)
   	 enddo
      do ic = 9,10
        coef_ns(i,jS,ic) = coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic)
        coef_ns(i,jN,ic) = coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic)
      enddo

   ! don't add values from the other hemisphere
   	 do ic = 6,8
   	   coef_ns(i,jS,ic) = coef_ns2(i,j,1,ic)
   	   coef_ns(i,jN,ic) = coef_ns2(i,j,2,ic)
   	 enddo

   ! from latm_JT to equator, add values from both hemispheres
   	 do j = jlatm_JT+1,nmlat_h-1
   	   jS = j
   	   jN = nmlat_T1-j+1

   	   do ic = 1,10
   	     coef_ns(i,jS,ic) = coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic)
   	     coef_ns(i,jN,ic) = coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic)
   	   enddo
   	 enddo

   ! set equatorial boundary condition (page 14 Art's notes)
   ! there should be just one equator value
   	 j = nmlat_h
   	 coef_ns(i,j,1) = coef_ns2(i,j,1,1)+coef_ns2(i,j,2,1)
   	 coef_ns(i,j,5) = coef_ns2(i,j,1,5)+coef_ns2(i,j,2,5)
   	 coef_ns(i,j,9) = coef_ns2(i,j,1,9)+coef_ns2(i,j,2,9)
   	 coef_ns(i,j,10) = coef_ns2(i,j,1,10)+coef_ns2(i,j,2,10)
   	 do ic = 6,8
   	   coef_ns(i,j,ic) = (coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic))/2._r8
   	 enddo

   ! Q from Wu: Why isn't it (coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic))/2?
   	 do ic = 2,4
   	   coef_ns(i,j,ic) = coef_ns(i,j,10-ic) ! 2,3,4 <- 8,7,6
   	 enddo
       enddo

     end subroutine edyn3D_calculate_coef_ns

!-----------------------------------------------------------------------
    subroutine edyn3D_calculate_bij(coef_ns2,bij)
    !     subroutine calculate_bij(coef_ns2,bij)
    ! set field-aligned conductance (b) matrix

       use edyn3D_params, only:nmlat_h,rho,rho_s,jlatm_JT

       real(r8),dimension(mlon0_p:mlon1_p,nmlat_h),intent(out) :: bij
       real(r8),dimension(mlon0_p:mlon1_p,nmlat_h,2,10),intent(in) :: coef_ns2

       real(r8),parameter :: &

   ! b_mult is [|Phi|/Delta(Phi)]*(R/L)^2, where Phi is a characteristic potential value,
   ! Delta(Phi) is a characteristic allowed interhemispheric potential difference,
   ! R is Earth radius, and L is a characteristic N-S length scale for Phi.
   ! It is assumed that b_mult is similar for middle and auroral latitudes.
   	 b_mult = 1e3_r8, &

   ! pccolatrad is the polar cap colatitude in radians, which for now is fixed.
   ! But it can be made variable w.r.t. time and magnetic longitude in the future.
   	 pccolatrad = 0.25_r8, & ! 14 degree
   	 rho_pc = sin(pccolatrad)

       integer :: i,j
       real(r8) :: fac3

   ! initialize bij to zero
   ! this also sets bij in low latitudes where fieldlines are assumed to be equipotential
   ! and it should not be used there (zero would actually mean two hemispheres are uncoupled)
       bij = 0

       do i = mlon0_p,mlon1_p
   	 do j = 2,jlatm_JT
   	   bij(i,j) = b_mult*(rho_s(j,1)-rho_s(j-1,1))**2/ &
   	     (1/(coef_ns2(i,j,1,3)+coef_ns2(i,j,1,7))+ &
   	      1/(coef_ns2(i,j,2,3)+coef_ns2(i,j,2,7)))
!   	   bij(i,j) = b_mult*(rho_s(j,1)-rho_s(j-1,1))**2/ &
!   	     (1/(coef_ns2(3,1,j,i)+coef_ns2(7,1,j,i))+ &
!   	      1/(coef_ns2(3,2,j,i)+coef_ns2(7,2,j,i)))
   	 enddo
       enddo

   ! set bij to zero within polar caps, transitioning linearly
   ! to the full original value over a distance of about (1/3) pccolatrad
       do j = 1,jlatm_JT
   	 fac3 = 3*(rho(j,1)/rho_pc-1)

   	 if (fac3 <= 0) then
   	   do i = mlon0_p,mlon1_p
   	     bij(i,j) = 0._r8
   	   enddo
   	 endif

   ! the same at conjugate points since bij is the same
   	 if (fac3>0 .and. fac3<1) then
   	   do i = mlon0_p,mlon1_p
   	     bij(i,j) = fac3*bij(i,j)
   	   enddo
   	 endif

   ! if fac3>=1 bij remains unmodified
       enddo

     end subroutine edyn3D_calculate_bij

!-----------------------------------------------------------------------------

   end module edyn3D_calculate_coefs
