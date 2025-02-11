module edyn3D_serial_solver

  use shr_kind_mod,   only: r8 => shr_kind_r8            ! 8-byte reals
  use edyn3D_params,  only: nmlat_T1,nmlon

  implicit none

  integer,parameter :: nlonlat = nmlat_T1*nmlon ! solve the whole globe

  contains
!-----------------------------------------------------------------------
  subroutine linear_system(bij_full,pot_hl_full,fac_hl_full,coef_ns_full,pot_full)
!
! construct linear system based on bij and coef_ns and solve for potential
!
! if FAC is read in, pot_hl is not used, only fac_hl is used
! if potential is read in, pot_hl is used, fac_hl is output

    use edyn3D_params, only: nmlat_h,nmlat_T1
    use edyn3D_mpi,    only: mlon0_p,mlon1_p,mp_gather_edyn3D,mp_scatter_edyn3D,mytid
     use cam_logfile,    only: iulog

    real(r8),dimension(nmlat_h,nmlon),intent(in) :: bij_full
    real(r8),dimension(2,nmlat_h,nmlon),intent(in) :: pot_hl_full
    real(r8),dimension(2,nmlat_h,nmlon),intent(out) :: fac_hl_full
    real(r8),dimension(10,nmlat_T1,nmlon),intent(inout) :: coef_ns_full
    real(r8),dimension(2,nmlat_h,nmlon),intent(out) :: pot_full

    logical,parameter :: use_mkl = .false.

!    logical,parameter :: use_mkl = &
!#ifdef MKL
!      .true.
!#else
!      .false.
!#endif

    integer,parameter :: root = 0
    integer :: i,j,jj,isn,nnz,ncnt1
    integer,dimension(10*nlonlat) :: irow,jcol
    real(r8),dimension(10*nlonlat) :: values
    real(r8),dimension(nlonlat) :: rhs,z,pot_hl_f,sol

! Q from Wu: why coefficients at the south pole are gathered to i=1?
      j = 1
      coef_ns_full(9,j,1) = sum(coef_ns_full(9,j,:))
      coef_ns_full(10,j,1) = sum(coef_ns_full(10,j,:))
      do i = 2,nmlon
        coef_ns_full(9,j,i) = 0
        coef_ns_full(10,j,i) = 0
      enddo

!!$ write(iulog,*) "linear_system: min/max pot_hl_full ", MINVAL(pot_hl_full), MAXVAL(pot_hl_full)
!!$
!!$ write(iulog,*) "linear_system: min/max fac_hl_full ", MINVAL(fac_hl_full), MAXVAL(fac_hl_full)
!!$
!!$ write(iulog,*) "linear_system: min/max coef_ns_full ", MINVAL(coef_ns_full), MAXVAL(coef_ns_full)
!!$
!!$ write(iulog,*) "linear_system: coef_ns_full(9,1,1), coef_ns_full(10,1,1) ", coef_ns_full(9,1,1), coef_ns_full(10,1,1)
!!$
! construct LHS matrix in COO format
      call construct_lhs(bij_full,coef_ns_full(1:9,:,:),nnz,irow,jcol,values)

! RHS is dense
      call construct_rhs(coef_ns_full(10,:,:),rhs)

! determine FAC forcing (dense)
!      if (read_fac) then ! input is corrected fac_hl, pot_hl is not used
!        z = flatten(fac_hl_full)
!
!      else ! input is pot_hl, fac_hl is to be calculated (output)

! A. Maute 2023/11/21: put the high latitude potential in X
! and then use LHS to calculate the RHS FAC
        pot_hl_f = flatten(pot_hl_full)

! z = matmul(lhs, pot_hl)
        z = 0
        do i = 1,nnz
          z(irow(i)) = z(irow(i))+values(i)*pot_hl_f(jcol(i))
        enddo

! no need for correction since it is from the divergence of horizontal current

! reconstruct 2D distribution of FAC based on z
        fac_hl_full(:,:,1:nmlon) = unravel(z)
!        fac_hl_full(:,:,1:nmlon) = 0.0_r8     ! test: settin HL fac to 0

!! add periodic points
!        do j = 1,nmlat_h
!          do isn = 1,2
!            fac_hl_full(isn,j,0) = fac_hl_full(isn,j,nmlon)
!            fac_hl_full(isn,j,nmlon+1) = fac_hl_full(isn,j,1)
!          enddo
!        enddo
!      endif

! add FAC forcing to RHS
      do i = 1,nlonlat
        rhs(i) = rhs(i)+z(i)
      enddo

!!$ write(iulog,*) "linear_system: nlonlat,nnz,min/max irow,min/max jcol,min/max values,min/max rhs ", nlonlat,nnz,MINVAL(irow), MAXVAL(irow), MINVAL(jcol), MAXVAL(jcol), MINVAL(values), MAXVAL(values), MINVAL(rhs), MAXVAL(rhs)

      if (use_mkl) then
        call solve_mkl(nlonlat,nnz,irow(1:nnz),jcol(1:nnz),values(1:nnz),rhs,sol)
      else
        call solve_superlu(nlonlat,nnz,irow(1:nnz),jcol(1:nnz),values(1:nnz),rhs,sol)
      endif

! reconstruct 2D distribution of potential based on the solution
      pot_full(:,:,1:nmlon) = unravel(sol)

!!$ write(iulog,*) "linear_system: min/max pot_full ", MINVAL(pot_full), MAXVAL(pot_full)

!    call mp_scatter_edyn3D(pot_2r,mlon0_p,mlon1_p,pot,nmlon+2,nmlat_h,2)
!
!    pot_full(mlon0_p-1,:,1) = pot_2r(mlon0_p-1,:,1)
!    pot(mlon0_p-1,:,2) = pot_2r(mlon0_p-1,:,2)
!    pot(mlon1_p+1,:,1) = pot_2r(mlon1_p+1,:,1)
!    pot(mlon1_p+1,:,2) = pot_2r(mlon1_p+1,:,2)
!
!    call mp_scatter_edyn3D(fac_hl_2r,mlon0_p,mlon1_p,fac_hl,nmlon+2,nmlat_h,2)
!
!    fac_hl(mlon0_p-1,:,1) = fac_hl_2r(mlon0_p-1,:,1)
!    fac_hl(mlon0_p-1,:,2) = fac_hl_2r(mlon0_p-1,:,2)
!    fac_hl(mlon1_p+1,:,1) = fac_hl_2r(mlon1_p+1,:,1)
!    fac_hl(mlon1_p+1,:,2) = fac_hl_2r(mlon1_p+1,:,2)

  endsubroutine linear_system
!-----------------------------------------------------------------------
  subroutine construct_lhs(bij,coef,nnz,irow,jcol,values)
! construct LHS matrix (sparse, coordinate form)

! need to set where the two hemispheres are connected
! This is not in the 9-point stencil but needs to be done manually

    use edyn3D_params, only: nmlat_h,jlatm_JT
!    use cons_module,only:jlatm_JT

    real(r8),dimension(nmlat_h,nmlon),intent(in) :: bij
    real(r8),dimension(9,nmlat_T1,nmlon),intent(in) :: coef
    integer,intent(out) :: nnz
    integer,dimension(10*nlonlat),intent(out) :: irow,jcol
    real(r8),dimension(10*nlonlat),intent(out) :: values

! if two hemispheres are uncoupled at high latitudes, set bijSum to zero
    real(r8),parameter :: bijSum = 0
    integer :: i,j,jS,jN,ij,it,im,ip,itcon

    nnz = 0

! set up poles, there are no c6,c7,c8 values
    j = 1
    jS = j
    jN = nmlat_T1-j+1

! for equation i=1 at the south pole
    i = 1
    ij = (i-1)*nmlat_T1+jS

! A. Richmond 2023/06/20: add bij effects
    it = (i-1)*nmlat_T1+jS
!    lhs(ij,it) = coef(9,jS,i)-bijSum
    nnz = nnz+1
    irow(nnz) = ij
    jcol(nnz) = it
    values(nnz) = coef(9,jS,i)-bijSum

! conjugate point (north pole at i=1)
    it = (i-1)*nmlat_T1+jN
!    lhs(ij,it) = bijSum
    nnz = nnz+1
    irow(nnz) = ij
    jcol(nnz) = it
    values(nnz) = bijSum

    do i = 1,nmlon ! Sum_i=1^nmlon C3(i,1) Phi(i,2)
      it = (i-1)*nmlat_T1+jS+1
!      lhs(ij,it) = coef(3,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(3,jS,i)
    enddo

! for other equations (i/=1) at the south pole
    do i = 2,nmlon
      ij = (i-1)*nmlat_T1+jS

      it = (i-1)*nmlat_T1+jS
!      lhs(ij,it) = 1
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = 1

      it = (1-1)*nmlat_T1+jS ! Phi(i,1) = Phi(1,1)
!      lhs(ij,it) = -1
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = -1
    enddo

! for each i, set Phi(i,nmlat_T1) = phi_pol at the north pole
    do i = 1,nmlon
      ij = (i-1)*nmlat_T1+jN

      it = (i-1)*nmlat_T1+jN
!      lhs(ij,it) = 1
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = 1
    enddo

    do i = 1,nmlon
      if (i == 1) then
        im = nmlon
      else
        im = i-1
      endif
      if (i == nmlon) then
        ip = 1
      else
        ip = i+1
      endif

! from pole to latm_JT, two hemispheres are uncoupled
      do j = 2,jlatm_JT-1
        jS = j
        jN = nmlat_T1-j+1

        ij = (i-1)*nmlat_T1+jS

        it = (ip-1)*nmlat_T1+jS
!        lhs(ij,it) = coef(1,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(1,jS,i)

        it = (ip-1)*nmlat_T1+jS+1
!        lhs(ij,it) = coef(2,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(2,jS,i)

        it = (i-1)*nmlat_T1+jS+1
!        lhs(ij,it) = coef(3,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(3,jS,i)

        it = (im-1)*nmlat_T1+jS+1
!        lhs(ij,it) = coef(4,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(4,jS,i)

        it = (im-1)*nmlat_T1+jS
!        lhs(ij,it) = coef(5,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(5,jS,i)

        it = (im-1)*nmlat_T1+jS-1
!        lhs(ij,it) = coef(6,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(6,jS,i)

        it = (i-1)*nmlat_T1+jS-1
!        lhs(ij,it) = coef(7,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(7,jS,i)

        it = (ip-1)*nmlat_T1+jS-1
!        lhs(ij,it) = coef(8,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(8,jS,i)

        it = (i-1)*nmlat_T1+jS
!        lhs(ij,it) = coef(9,jS,i)-bij(j,i) ! should be 1
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(9,jS,i)-bij(j,i)

        itcon = (i-1)*nmlat_T1+jN ! conjugate point
!        lhs(ij,itcon) = bij(j,i) ! b(i,j) Phi*(i,j)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = itcon
        values(nnz) = bij(j,i)

! note that coef has the direction switched in the NH
        ij = (i-1)*nmlat_T1+jN

        it = (ip-1)*nmlat_T1+jN
!        lhs(ij,it) = coef(1,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(1,jN,i)

        it = (ip-1)*nmlat_T1+jN+1
!        lhs(ij,it) = coef(8,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(8,jN,i)

        it = (i-1)*nmlat_T1+jN+1
!        lhs(ij,it) = coef(7,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(7,jN,i)

        it = (im-1)*nmlat_T1+jN+1
!        lhs(ij,it) = coef(6,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(6,jN,i)

        it = (im-1)*nmlat_T1+jN
!        lhs(ij,it) = coef(5,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(5,jN,i)

        it = (im-1)*nmlat_T1+jN-1
!        lhs(ij,it) = coef(4,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(4,jN,i)

        it = (i-1)*nmlat_T1+jN-1
!        lhs(ij,it) = coef(3,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(3,jN,i)

        it = (ip-1)*nmlat_T1+jN-1
!        lhs(ij,it) = coef(2,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(2,jN,i)

        it = (i-1)*nmlat_T1+jN
!        lhs(ij,it) = coef(9,jN,i)-bij(j,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(9,jN,i)-bij(j,i)

        itcon = (i-1)*nmlat_T1+jS ! conjugate point
!        lhs(ij,itcon) = bij(j,i) ! b(i,j) Phi*(i,j)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = itcon
        values(nnz) = bij(j,i)
      enddo

! at latm_JT, two hemispheres are coupled at j-1
      j = jlatm_JT
      jS = j
      jN = nmlat_T1-j+1

      ij = (i-1)*nmlat_T1+jS

      it = (ip-1)*nmlat_T1+jS
!      lhs(ij,it) = coef(1,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(1,jS,i)

      it = (ip-1)*nmlat_T1+jS+1
!      lhs(ij,it) = coef(2,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(2,jS,i)

      it = (i-1)*nmlat_T1+jS+1
!      lhs(ij,it) = coef(3,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(3,jS,i)

      it = (im-1)*nmlat_T1+jS+1
!      lhs(ij,it) = coef(4,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(4,jS,i)

      it = (im-1)*nmlat_T1+jS
!      lhs(ij,it) = coef(5,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(5,jS,i)

      it = (im-1)*nmlat_T1+jS-1
!      lhs(ij,it) = coef(6,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(6,jS,i)

      it = (i-1)*nmlat_T1+jS-1
!      lhs(ij,it) = coef(7,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(7,jS,i)

      it = (ip-1)*nmlat_T1+jS-1
!      lhs(ij,it) = coef(8,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(8,jS,i)

      it = (i-1)*nmlat_T1+jS
!      lhs(ij,it) = coef(9,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(9,jS,i)

! note that coef has the direction switched in the NH
! therefore the original coef_ns2 c6,c7,c8 at j-1 poleward
! become coef_ns2 c4,c3,c2 at j'+1 poleward
      itcon = (im-1)*nmlat_T1+jN+1
!      lhs(ij,itcon) = coef(6,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = itcon
      values(nnz) = coef(6,jN,i)

      itcon = (i-1)*nmlat_T1+jN+1
!      lhs(ij,itcon) = coef(7,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = itcon
      values(nnz) = coef(7,jN,i)

      itcon = (ip-1)*nmlat_T1+jN+1
!      lhs(ij,itcon) = coef(8,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = itcon
      values(nnz) = coef(8,jN,i)

! note that coef has the direction switched in the NH
      ij = (i-1)*nmlat_T1+jN

      it = (ip-1)*nmlat_T1+jN
!      lhs(ij,it) = coef(1,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(1,jN,i)

      it = (ip-1)*nmlat_T1+jN+1
!      lhs(ij,it) = coef(8,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(8,jN,i)

      it = (i-1)*nmlat_T1+jN+1
!      lhs(ij,it) = coef(7,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(7,jN,i)

      it = (im-1)*nmlat_T1+jN+1
!      lhs(ij,it) = coef(6,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(6,jN,i)

      it = (im-1)*nmlat_T1+jN
!      lhs(ij,it) = coef(5,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(5,jN,i)

      it = (im-1)*nmlat_T1+jN-1
!      lhs(ij,it) = coef(4,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(4,jN,i)

      it = (i-1)*nmlat_T1+jN-1
!      lhs(ij,it) = coef(3,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(3,jN,i)

      it = (ip-1)*nmlat_T1+jN-1
!      lhs(ij,it) = coef(2,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(2,jN,i)

      it = (i-1)*nmlat_T1+jN
!      lhs(ij,it) = coef(9,jN,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(9,jN,i)

! SH coef does not have a direction switch
! therefore coef at j-1 is still coef at j'-1
      itcon = (im-1)*nmlat_T1+jS-1
!      lhs(ij,itcon) = coef(6,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = itcon
      values(nnz) = coef(6,jS,i)

      itcon = (i-1)*nmlat_T1+jS-1
!      lhs(ij,itcon) = coef(7,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = itcon
      values(nnz) = coef(7,jS,i)

      itcon = (ip-1)*nmlat_T1+jS-1
!      lhs(ij,itcon) = coef(8,jS,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = itcon
      values(nnz) = coef(8,jS,i)

! from latm_JT to equator, symmetric solution
      do j = jlatm_JT+1,nmlat_h-1
        jS = j
        jN = nmlat_T1-j+1

        ij = (i-1)*nmlat_T1+jS

        it = (ip-1)*nmlat_T1+jS
!        lhs(ij,it) = coef(1,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(1,jS,i)

        it = (ip-1)*nmlat_T1+jS+1
!        lhs(ij,it) = coef(2,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(2,jS,i)

        it = (i-1)*nmlat_T1+jS+1
!        lhs(ij,it) = coef(3,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(3,jS,i)

        it = (im-1)*nmlat_T1+jS+1
!        lhs(ij,it) = coef(4,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(4,jS,i)

        it = (im-1)*nmlat_T1+jS
!        lhs(ij,it) = coef(5,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(5,jS,i)

        it = (im-1)*nmlat_T1+jS-1
!        lhs(ij,it) = coef(6,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(6,jS,i)

        it = (i-1)*nmlat_T1+jS-1
!        lhs(ij,it) = coef(7,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(7,jS,i)

        it = (ip-1)*nmlat_T1+jS-1
!        lhs(ij,it) = coef(8,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(8,jS,i)

        it = (i-1)*nmlat_T1+jS
!        lhs(ij,it) = coef(9,jS,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(9,jS,i)

! note that coef has the direction switched in the NH
        ij = (i-1)*nmlat_T1+jN

        it = (ip-1)*nmlat_T1+jN
!        lhs(ij,it) = coef(1,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(1,jN,i)

        it = (ip-1)*nmlat_T1+jN+1
!        lhs(ij,it) = coef(8,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(8,jN,i)

        it = (i-1)*nmlat_T1+jN+1
!        lhs(ij,it) = coef(7,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(7,jN,i)

        it = (im-1)*nmlat_T1+jN+1
!        lhs(ij,it) = coef(6,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(6,jN,i)

        it = (im-1)*nmlat_T1+jN
!        lhs(ij,it) = coef(5,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(5,jN,i)

        it = (im-1)*nmlat_T1+jN-1
!        lhs(ij,it) = coef(4,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(4,jN,i)

        it = (i-1)*nmlat_T1+jN-1
!        lhs(ij,it) = coef(3,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(3,jN,i)

        it = (ip-1)*nmlat_T1+jN-1
!        lhs(ij,it) = coef(2,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(2,jN,i)

        it = (i-1)*nmlat_T1+jN
!        lhs(ij,it) = coef(9,jN,i)
        nnz = nnz+1
        irow(nnz) = ij
        jcol(nnz) = it
        values(nnz) = coef(9,jN,i)
      enddo

! equator
      j = nmlat_h

      ij = (i-1)*nmlat_T1+j

      it = (ip-1)*nmlat_T1+j
!      lhs(ij,it) = coef(1,j,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(1,j,i)

      it = (ip-1)*nmlat_T1+j+1
!      lhs(ij,it) = coef(2,j,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(2,j,i)

      it = (i-1)*nmlat_T1+j+1
!      lhs(ij,it) = coef(3,j,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(3,j,i)

      it = (im-1)*nmlat_T1+j+1
!      lhs(ij,it) = coef(4,j,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(4,j,i)

      it = (im-1)*nmlat_T1+j
!      lhs(ij,it) = coef(5,j,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(5,j,i)

      it = (im-1)*nmlat_T1+j-1
!      lhs(ij,it) = coef(6,j,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(6,j,i)

      it = (i-1)*nmlat_T1+j-1
!      lhs(ij,it) = coef(7,j,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(7,j,i)

      it = (ip-1)*nmlat_T1+j-1
!      lhs(ij,it) = coef(8,j,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(8,j,i)

      it = (i-1)*nmlat_T1+j
!      lhs(ij,it) = coef(9,j,i)
      nnz = nnz+1
      irow(nnz) = ij
      jcol(nnz) = it
      values(nnz) = coef(9,j,i)
    enddo

  endsubroutine construct_lhs
!-----------------------------------------------------------------------
  subroutine construct_rhs(coef_10,rhs)
! construct vector RHS
! this is different from flatten

    use edyn3D_params, only: nmlat_h,phi_pol
!    use cons_module,only:phi_pol

    real(r8),dimension(nmlat_T1,nmlon),intent(in) :: coef_10
    real(r8),dimension(nlonlat),intent(out) :: rhs

    integer :: i,j,ij

    rhs = 0

! set up poles, there are no c6,c7,c8 values
    j = 1

! for equation i=1 at the south pole
    i = 1
    ij = (i-1)*nmlat_T1+j
    rhs(ij) = coef_10(j,i)

! for each i, set Phi(i,nmlat_T1) = phi_pol at the north pole
    do i = 1,nmlon
      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      rhs(ij) = phi_pol
    enddo

    do i = 1,nmlon
      do j = 2,nmlat_h-1
        ij = (i-1)*nmlat_T1+j
        rhs(ij) = coef_10(j,i)

        ij = (i-1)*nmlat_T1+nmlat_T1-j+1
        rhs(ij) = coef_10(nmlat_T1-j+1,i)
      enddo

! equator
      j = nmlat_h

      ij = (i-1)*nmlat_T1+j
      rhs(ij) = coef_10(j,i)
    enddo

  endsubroutine construct_rhs
!-----------------------------------------------------------------------
  pure function flatten(fin) result(fout)
! reorder 2D fields (lat-lon) into 1D vector (RHS)
! northern/southern hemispheres are either separate or averaged
! based on their latitude ranges (high-lat, transition, low-lat, equator)

    use edyn3D_params, only: nmlat_h,jlatm_JT
!    use cons_module,only:jlatm_JT

    real(r8),dimension(2,nmlat_h,nmlon),intent(in) :: fin
    real(r8),dimension(nlonlat) :: fout

    integer :: i,j,ij

    do i = 1,nmlon

! from pole to latm_JT, two hemispheres are uncoupled
      do j = 1,jlatm_JT
        ij = (i-1)*nmlat_T1+j
        fout(ij) = fin(1,j,i)

        ij = (i-1)*nmlat_T1+nmlat_T1-j+1
        fout(ij) = fin(2,j,i)
      enddo

! from latm_JT to equator, symmetric solution
      do j = jlatm_JT+1,nmlat_h-1
        ij = (i-1)*nmlat_T1+j
        fout(ij) = (fin(1,j,i)+fin(2,j,i))/2

        ij = (i-1)*nmlat_T1+nmlat_T1-j+1
        fout(ij) = (fin(1,j,i)+fin(2,j,i))/2
      enddo

! equator
      j = nmlat_h
      ij = (i-1)*nmlat_T1+j
      fout(ij) = fin(1,j,i)
    enddo

  endfunction flatten
!-----------------------------------------------------------------------
  pure function unravel(fin) result(fout)
! reorder 1D vector (RHS) into 2D fields (lat-lon)

    use edyn3D_params, only: nmlat_h

    real(r8),dimension(nlonlat),intent(in) :: fin
    real(r8),dimension(2,nmlat_h,nmlon) :: fout

    integer :: i,j,isn,ij

    do i = 1,nmlon
      do j = 1,nmlat_h
        do isn = 1,2
          if (isn == 1) then
            ij = (i-1)*nmlat_T1+j
          else
            ij = (i-1)*nmlat_T1+nmlat_T1-j+1
          endif
          fout(isn,j,i) = fin(ij)
        enddo
      enddo
    enddo

  endfunction unravel
!-----------------------------------------------------------------------
  subroutine solve_mkl(n,nnz,irow,jcol,values,rhs,sol)

#ifdef MKL
    include 'mkl_pardiso.fi'
#endif

    integer,intent(in) :: n,nnz
    integer,dimension(nnz),intent(in) :: irow,jcol
    real(r8),dimension(nnz),intent(in) :: values
    real(r8),dimension(n),intent(in) :: rhs
    real(r8),dimension(n),intent(out) :: sol

#ifdef MKL
! for PARDISO sparse matrix solver (CSR format)
    logical,parameter :: debug = .false.
    integer,parameter :: maxfct = 1, mnum = 1, nrhs = 1, &
      mtype = 11, & ! real and non-symmetric matrix
      msglvl = 1 ! print statistical information
    integer :: i,phase,error
    integer,dimension(64) :: iparm
    integer,dimension(n) :: perm
    integer,dimension(n+1) :: rowptr
    integer,dimension(nnz) :: colind
    real(r8),dimension(n) :: rhs_cp
    real(r8),dimension(nnz) :: nzval
    type(MKL_PARDISO_HANDLE),dimension(64) :: pt

! MKL needs CSR format
    call coo_to_csr(n,n,nnz,irow,jcol,values,rowptr,colind,nzval)

! initialize PARDISO with default parameters in accordance with the matrix type
    call pardisoinit(pt,mtype,iparm)

! set some nonzero iparm elements
    if (debug) then

! 0: iparm(2) - iparm(64) are filled with default values
      iparm(1) = 1

! report the number of non-zero elements in the factors
      iparm(18) = -1

! report number of floating point operations (in 10^6 floating point operations)
! that are necessary to factor the matrix A
      iparm(19) = -1

! matrix checker, 1 checks integer arrays rowptr and colind
      iparm(27) = 1
    endif

! fill permutation vector with zero (not used)
    perm = 0

! analysis, numerical factorization, solve
! this is the same as doing calls to 11, 22, 33 in order
    phase = 13

! pardiso doesn't have intent(in) attribute for rhs, so make a copy
    do i = 1,n
      rhs_cp(i) = rhs(i)
    enddo

! use 32-bit integer version
! if the number of non-zero elements is on the order of 500 million or more
! then use pardiso_64 (64-bit integer version)
    call pardiso(pt, maxfct, mnum, mtype, phase, n, &
      nzval, rowptr, colind, perm, nrhs, iparm, msglvl, rhs_cp, sol, error)

    write(6,"('phase ',i4,' error ',i4)") phase,error
#endif

  endsubroutine solve_mkl
!-----------------------------------------------------------------------
  subroutine solve_superlu(n,nnz,irow,jcol,values,rhs,sol)
    use iso_c_binding, only: c_long_long

    integer,intent(in) :: n,nnz
    integer,dimension(nnz),intent(in) :: irow,jcol
    real(r8),dimension(nnz),intent(in) :: values
    real(r8),dimension(n),intent(in) :: rhs
    real(r8),dimension(n),intent(out) :: sol

! for SuperLU sparse matrix solver (CSC format)
    integer,parameter :: nrhs = 1
    integer :: i,iopt,info
    integer(kind=c_long_long) :: f_factors
    integer,dimension(n+1) :: colptr
    integer,dimension(nnz) :: rowind
    real(r8),dimension(nnz) :: nzval

    interface
      subroutine c_fortran_dgssv_(iopt,n,nnz,nrhs,values,rowind,colptr,b,ldb,f_factors,info) bind(c,name='c_fortran_dgssv_')
        use iso_c_binding,only:c_int,c_long_long,c_double
        integer(kind=c_int) :: iopt,n,nnz,nrhs,ldb,info
        real(kind=c_double),dimension(nnz) :: values
        integer(kind=c_int),dimension(nnz) :: rowind
        integer(kind=c_int),dimension(n+1) :: colptr
        real(kind=c_double),dimension(ldb) :: b
        integer(kind=c_long_long) :: f_factors
      endsubroutine c_fortran_dgssv_
    endinterface

! SuperLU needs CSC format
    call coo_to_csc(n,n,nnz,irow,jcol,values,colptr,rowind,nzval)

    do i = 1,n
      sol(i) = rhs(i)
    enddo

! first, factorize the matrix, the factors are stored in *f_factors* handle
    iopt = 1
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      nzval, rowind, colptr, sol, n, f_factors, info)
    write(6,"('INFO from LU decomposition = ',i4)") info

! second, solve the system using the existing factors
    iopt = 2
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      nzval, rowind, colptr, sol, n, f_factors, info)
    write(6,"('INFO from triangular solve = ',i4)") info

! last, free the storage allocated inside SuperLU
    iopt = 3
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      nzval, rowind, colptr, sol, n, f_factors, info)

  endsubroutine solve_superlu
!-----------------------------------------------------------------------
  subroutine coo_to_csr(nrow,ncol,nnz,irow,jcol,values,rowptr,colind,nzval)
! convert sparse matrix format from COO to CSR

! Compressed Sparse Row (CSR) format:
! rowptr is the cumulative sum of the number of non-zero elements in each row
! i.e., rowptr(i+1)-rowptr(i) is the number of non-zero elements in the ith row
! colind(rowptr(i):rowptr(i+1)-1) are the column indices of non-zero elements in the ith column

    integer,intent(in) :: nrow,ncol,nnz
    integer,dimension(nnz),intent(in) :: irow,jcol
    real(r8),dimension(nnz),intent(in) :: values
    integer,dimension(nrow+1),intent(out) :: rowptr
    integer,dimension(nnz),intent(out) :: colind
    real(r8),dimension(nnz),intent(out) :: nzval

    integer :: i,j,last_colind,nnz_i
    integer,dimension(ncol) :: jcol_i,idx
    real(r8),dimension(ncol) :: values_i

! start with the first row
    last_colind = 1
    rowptr(1) = last_colind

! loop through rows
    do i = 1,nrow

! squeeze column indices only in the ith row
      nnz_i = 0
      do j = 1,nnz
        if (irow(j) == i) then
          nnz_i = nnz_i+1
          jcol_i(nnz_i) = jcol(j)
          values_i(nnz_i) = values(j)
        endif
      enddo

! if nnz_i is zero, no element is in the ith row, the matrix is degenerate

! sort values based on column indices
      if (nnz_i > 0) then
        idx(1:nnz_i) = argsort(jcol_i(1:nnz_i),nnz_i,1,ncol)
        do j = 1,nnz_i
          colind(last_colind+j-1) = jcol_i(idx(j))
          nzval(last_colind+j-1) = values_i(idx(j))
        enddo
      endif

! update the pointer to the new row
      last_colind = last_colind+nnz_i
      rowptr(i+1) = last_colind
    enddo

! last_colind should equal nnz+1 at this point

  endsubroutine coo_to_csr
!-----------------------------------------------------------------------
  subroutine coo_to_csc(nrow,ncol,nnz,irow,jcol,values,colptr,rowind,nzval)
! convert sparse matrix format from COO to CSC

! Compressed Sparse Column (CSC) format:
! colptr is the cumulative sum of the number of non-zero elements in each column
! i.e., colptr(j+1)-colptr(j) is the number of non-zero elements in the jth column
! rowind(colptr(j):colptr(j+1)-1) are the row indices of non-zero elements in the jth row

    integer,intent(in) :: nrow,ncol,nnz
    integer,dimension(nnz),intent(in) :: irow,jcol
    real(r8),dimension(nnz),intent(in) :: values
    integer,dimension(ncol+1),intent(out) :: colptr
    integer,dimension(nnz),intent(out) :: rowind
    real(r8),dimension(nnz),intent(out) :: nzval

    integer :: i,j,last_rowind,nnz_j
    integer,dimension(nrow) :: irow_j,idx
    real(r8),dimension(nrow) :: values_j

! start with the first column
    last_rowind = 1
    colptr(1) = last_rowind

! loop through columns
    do j = 1,ncol

! squeeze row indices only in the jth column
      nnz_j = 0
      do i = 1,nnz
        if (jcol(i) == j) then
          nnz_j = nnz_j+1
          irow_j(nnz_j) = irow(i)
          values_j(nnz_j) = values(i)
        endif
      enddo

! if nnz_j is zero, no element is in the jth column, the matrix is degenerate

! sort values based on row indices
      if (nnz_j > 0) then
        idx(1:nnz_j) = argsort(irow_j(1:nnz_j),nnz_j,1,nrow)
        do i = 1,nnz_j
          rowind(last_rowind+i-1) = irow_j(idx(i))
          nzval(last_rowind+i-1) = values_j(idx(i))
        enddo
      endif

! update the pointer to the new column
      last_rowind = last_rowind+nnz_j
      colptr(j+1) = last_rowind
    enddo

! last_rowind should equal nnz+1 at this point

  endsubroutine coo_to_csc
!-----------------------------------------------------------------------
  pure function argsort(a,n,amin,amax) result(idx)
! a simplified radix sort for row or column indices

    integer,intent(in) :: n,amin,amax
    integer,dimension(n),intent(in) :: a
    integer,dimension(n) :: idx

    integer :: i,cnt
    integer,dimension(amin:amax) :: b

! put all elements in place
    b = amin-1
    do i = 1,n
      b(a(i)) = i
    enddo

! squeeze non-existent indices
    cnt = 1
    do i = amin,amax
      if (b(i) /= amin-1) then
        idx(cnt) = b(i)
        cnt = cnt+1
        if (cnt > n) exit
      endif
    enddo

  endfunction argsort
!-----------------------------------------------------------------------
end module edyn3D_serial_solver
