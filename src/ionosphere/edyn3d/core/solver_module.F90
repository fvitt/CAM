module solver_module
  use perf_mod, only: t_startf, t_stopf

  use prec,only:rp

  implicit none

  contains
!-----------------------------------------------------------------------
  subroutine linear_system(mlatd0,mlatd1,mlond0,mlond1, &
    bij,pot_hl,fac_hl,coef_ns,pot)
! construct linear system based on bij and coef and solve in pot

! if FAC is read in, pot_hl is not used, only fac_hl is used
! if potential is read in, pot_hl is used, fac_hl is output

    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:read_fac
    use mpi_module,only:gather_mag,bcast_3d,mpi_rank, dynamo_world

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_hl
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(inout) :: fac_hl
    real(kind=rp),dimension(10,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_ns
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: pot

    logical,parameter :: use_mkl = &
#ifdef MKL
      .true.
#else
      .false.
#endif
    integer,parameter :: root = 0
    integer :: nlonlat,mlat0,mlat1,mlon0,mlon1,i,j,isn,ic,nnz
    integer,dimension(nmlat_T1*nmlon+1) :: rowptr,colptr
    integer,dimension(12*nmlat_T1*nmlon) :: colind,rowind
    real(kind=rp),dimension(12*nmlat_T1*nmlon) :: values_csr,values_csc
    real(kind=rp),dimension(nmlat_h,nmlon) :: bij_full
    real(kind=rp),dimension(2,nmlat_h,nmlon) :: pot_hl_full,fac_hl_full
    real(kind=rp),dimension(10,2,nmlat_h,nmlon) :: coef_ns_full
    real(kind=rp),dimension(10,nmlat_T1,nmlon) :: coef_full
    real(kind=rp),dimension(2,nmlat_h,0:nmlon+1) :: fac_hl_2,pot_2
    real(kind=rp),dimension(nmlat_T1*nmlon) :: rhs,z,pot_hl_f,sol

    integer :: ier

    call t_startf('linear_system')

    nlonlat = nmlat_T1*nmlon ! solve the whole globe
    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    call t_startf('linear_system->gather_mag')

    bij_full = gather_mag(bij(mlat0:mlat1,mlon0:mlon1),root)
    if (read_fac) then
      fac_hl_full = gather_mag(fac_hl(:,mlat0:mlat1,mlon0:mlon1),2,root)
    else
      pot_hl_full = gather_mag(pot_hl(:,mlat0:mlat1,mlon0:mlon1),2,root)
    endif
    coef_ns_full = gather_mag(coef_ns(:,:,mlat0:mlat1,mlon0:mlon1),10,2,root)

    call t_stopf('linear_system->gather_mag')

    if (mpi_rank == root) then

! concatenate two hemispheres
      do concurrent (i = 1:nmlon, j = 1:nmlat_h, ic = 1:10)
        coef_full(ic,j,i) = coef_ns_full(ic,1,j,i)
        coef_full(ic,nmlat_T1-j+1,i) = coef_ns_full(ic,2,j,i)
      enddo

! construct LHS matrix in CSR format
      call construct_lhs(bij_full,coef_full(1:9,:,:),rowptr,colind,values_csr)
      nnz = rowptr(nlonlat+1)-1

! RHS is dense
      rhs = construct_rhs(coef_full(10,:,:))

! determine FAC forcing (dense)
      if (read_fac) then ! input is corrected fac_hl, pot_hl is not used
        z = flatten(fac_hl_full)

      else ! input is pot_hl, fac_hl is to be calculated (output)

! A. Maute 2023/11/21: put the high latitude potential in X
! and then use LHS to calculate the RHS FAC
        pot_hl_f = flatten(pot_hl_full)

! z = matmul(lhs, pot_hl)
        z = 0
        do i = 1,nlonlat
          do j = rowptr(i),rowptr(i+1)-1
            z(i) = z(i)+values_csr(j)*pot_hl_f(colind(j))
          enddo
        enddo

! no need for correction since it is from the divergence of horizontal current

! reconstruct 2D distribution of FAC based on z
        fac_hl_2(:,:,1:nmlon) = unravel(z)

! add periodic points
        do j = 1,nmlat_h
          do isn = 1,2
            fac_hl_2(isn,j,0) = fac_hl_2(isn,j,nmlon)
            fac_hl_2(isn,j,nmlon+1) = fac_hl_2(isn,j,1)
          enddo
        enddo
      endif

! add FAC forcing to RHS
      do i = 1,nlonlat
        rhs(i) = rhs(i)+z(i)
      enddo

      if (use_mkl) then
        sol = solve_mkl(nlonlat,nnz,rowptr,colind(1:nnz),values_csr(1:nnz),rhs)
      else
        call csr_to_csc(nlonlat,nlonlat,nnz, &
          rowptr,colind(1:nnz),values_csr(1:nnz), &
          colptr,rowind(1:nnz),values_csc(1:nnz))

        call t_startf('linear_system->solve_superlu')
        sol = solve_superlu(nlonlat,nnz,colptr,rowind(1:nnz),values_csc(1:nnz),rhs)
        call t_stopf('linear_system->solve_superlu')

      endif

! reconstruct 2D distribution of potential based on the solution
      pot_2(:,:,1:nmlon) = unravel(sol)

! periodic points
      do j = 1,nmlat_h
        do isn = 1,2
          pot_2(isn,j,0) = pot_2(isn,j,nmlon)
          pot_2(isn,j,nmlon+1) = pot_2(isn,j,1)
        enddo
      enddo
    endif

    call mpi_barrier (dynamo_world, ier)

    call t_startf('linear_system->bcast_3d')
    if (.not. read_fac) then
      call bcast_3d(fac_hl_2,2,nmlat_h,nmlon+2,root)

      do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
        fac_hl(isn,j,i) = fac_hl_2(isn,j,i)
      enddo
    endif

    call bcast_3d(pot_2,2,nmlat_h,nmlon+2,root)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
      pot(isn,j,i) = pot_2(isn,j,i)
    enddo
    call t_stopf('linear_system->bcast_3d')

    call t_stopf('linear_system')

  endsubroutine linear_system
!-----------------------------------------------------------------------
  pure subroutine construct_lhs(bij,coef,rowptr,colind,values)
! construct LHS matrix (CSR format)

! need to set where the two hemispheres are connected
! This is not in the 9-point stencil but needs to be done manually

    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(nmlat_h,nmlon),intent(in) :: bij
    real(kind=rp),dimension(9,nmlat_T1,nmlon),intent(in) :: coef
    integer,dimension(nmlat_T1*nmlon+1),intent(out) :: rowptr
    integer,dimension(12*nmlat_T1*nmlon),intent(out) :: colind
    real(kind=rp),dimension(12*nmlat_T1*nmlon),intent(out) :: values

! if two hemispheres are uncoupled at high latitudes, set bijSum to zero
    real(kind=rp),parameter :: bijSum = 0
    integer :: nlonlat,i,j,jS,jN,ij,isub,im,ip
    integer,dimension(nmlat_T1*nmlon) :: rowcnt

! the first row has most elements (nmlon+2)
    integer,dimension(nmlon+2) :: jcol1
    real(kind=rp),dimension(nmlon+2) :: nzval1

! other rows have at most 12 elements
    integer,dimension(12,nmlat_T1*nmlon) :: jcol
    real(kind=rp),dimension(12,nmlat_T1*nmlon) :: nzval

    nlonlat = nmlat_T1*nmlon
    rowcnt = 0

! set up poles
    j = 1

! there are no c6,c7,c8 values at the south pole

! for longitude i=1 at the south pole
    i = 1
!    ij = (i-1)*nmlat_T1+j
!    lhs(ij,(i-1)*nmlat_T1+j           ) = sum(coef(9,j,:))-bijSum ! A. Richmond 2023/06/20: add bij effects
!    lhs(ij,(i-1)*nmlat_T1+j+1         ) = coef(3,j,i)             ! Sum_i=1^nmlon C3(i,1) Phi(i,2)
!    lhs(ij,(i-1)*nmlat_T1+nmlat_T1-j+1) = bijSum                  ! conjugate point (north pole at i=1)
    jcol1(1:3) = (/(i-1)*nmlat_T1+j,(i-1)*nmlat_T1+j+1,(i-1)*nmlat_T1+nmlat_T1-j+1/)
    nzval1(1:3) = (/sum(coef(9,j,:))-bijSum,coef(3,j,i),bijSum/)
    do isub = 2,nmlon ! Sum_i=1^nmlon C3(i,1) Phi(i,2)
!      lhs(ij,(isub-1)*nmlat_T1+j+1) = coef(3,j,isub)
      jcol1(isub+2) = (isub-1)*nmlat_T1+j+1
      nzval1(isub+2) = coef(3,j,isub)
    enddo
    rowcnt(1) = nmlon+2

! for other longitudes (i/=1) at the south pole
    do i = 2,nmlon
      ij = (i-1)*nmlat_T1+j
!      lhs(ij,(1-1)*nmlat_T1+j) = -1 ! Phi(i,1) = Phi(1,1)
!      lhs(ij,(i-1)*nmlat_T1+j) = 1
      jcol(rowcnt(ij)+1:rowcnt(ij)+2,ij) = (/(1-1)*nmlat_T1+j,(i-1)*nmlat_T1+j/)
      nzval(rowcnt(ij)+1:rowcnt(ij)+2,ij) = (/-1,1/)
      rowcnt(ij) = rowcnt(ij)+2
    enddo

! for each i, set Phi(i,nmlat_T1) = phi_pol at the north pole
    do i = 1,nmlon
      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
!      lhs(ij,(i-1)*nmlat_T1+nmlat_T1-j+1) = 1
      rowcnt(ij) = rowcnt(ij)+1
      jcol(rowcnt(ij),ij) = (i-1)*nmlat_T1+nmlat_T1-j+1
      nzval(rowcnt(ij),ij) = 1
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
!        lhs(ij,(im-1)*nmlat_T1+jS-1) = coef(6,jS,i)
!        lhs(ij,(im-1)*nmlat_T1+jS  ) = coef(5,jS,i)
!        lhs(ij,(im-1)*nmlat_T1+jS+1) = coef(4,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS-1) = coef(7,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS  ) = coef(9,jS,i)-bij(j,i) ! should be 1
!        lhs(ij, (i-1)*nmlat_T1+jS+1) = coef(3,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jN  ) = bij(j,i)              ! conjugate point, b(i,j) Phi*(i,j)
!        lhs(ij,(ip-1)*nmlat_T1+jS-1) = coef(8,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS  ) = coef(1,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS+1) = coef(2,jS,i)
        if (i == 1) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN, &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1, &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(7,jS,i),coef(9,jS,i)-bij(j,i),coef(3,jS,i),bij(j,i), &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i), &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i)/)
        elseif (i == nmlon) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1, &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1, &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i), &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i), &
            coef(7,jS,i),coef(9,jS,i)-bij(j,i),coef(3,jS,i),bij(j,i)/)
        else
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1, &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN, &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i), &
            coef(7,jS,i),coef(9,jS,i)-bij(j,i),coef(3,jS,i),bij(j,i), &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i)/)
        endif
        rowcnt(ij) = rowcnt(ij)+10

! note that coef has the direction switched in the NH
        ij = (i-1)*nmlat_T1+jN
!        lhs(ij,(im-1)*nmlat_T1+jN-1) = coef(4,jN,i)
!        lhs(ij,(im-1)*nmlat_T1+jN  ) = coef(5,jN,i)
!        lhs(ij,(im-1)*nmlat_T1+jN+1) = coef(6,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jS  ) = bij(j,i)              ! conjugate point, b(i,j) Phi*(i,j)
!        lhs(ij, (i-1)*nmlat_T1+jN-1) = coef(3,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jN  ) = coef(9,jN,i)-bij(j,i)
!        lhs(ij, (i-1)*nmlat_T1+jN+1) = coef(7,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN-1) = coef(2,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN  ) = coef(1,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN+1) = coef(8,jN,i)
        if (i == 1) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            bij(j,i),coef(3,jN,i),coef(9,jN,i)-bij(j,i),coef(7,jN,i), &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i)/)
        elseif (i == nmlon) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
            (i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
            bij(j,i),coef(3,jN,i),coef(9,jN,i)-bij(j,i),coef(7,jN,i)/)
        else
          jcol(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
            (i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+10,ij) = (/ &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
            bij(j,i),coef(3,jN,i),coef(9,jN,i)-bij(j,i),coef(7,jN,i), &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i)/)
        endif
        rowcnt(ij) = rowcnt(ij)+10
      enddo

! at latm_JT, two hemispheres are coupled at j-1
      j = jlatm_JT
      jS = j
      jN = nmlat_T1-j+1

! note that coef has the direction switched in the NH
! therefore the original coef_ns2 c6,c7,c8 at j-1 poleward
! become coef_ns2 c4,c3,c2 at j'+1 poleward
      ij = (i-1)*nmlat_T1+jS
!      lhs(ij,(im-1)*nmlat_T1+jS-1) = coef(6,jS,i)
!      lhs(ij,(im-1)*nmlat_T1+jS  ) = coef(5,jS,i)
!      lhs(ij,(im-1)*nmlat_T1+jS+1) = coef(4,jS,i)
!      lhs(ij,(im-1)*nmlat_T1+jN+1) = coef(6,jN,i)
!      lhs(ij, (i-1)*nmlat_T1+jS-1) = coef(7,jS,i)
!      lhs(ij, (i-1)*nmlat_T1+jS  ) = coef(9,jS,i)
!      lhs(ij, (i-1)*nmlat_T1+jS+1) = coef(3,jS,i)
!      lhs(ij, (i-1)*nmlat_T1+jN+1) = coef(7,jN,i)
!      lhs(ij,(ip-1)*nmlat_T1+jS-1) = coef(8,jS,i)
!      lhs(ij,(ip-1)*nmlat_T1+jS  ) = coef(1,jS,i)
!      lhs(ij,(ip-1)*nmlat_T1+jS+1) = coef(2,jS,i)
!      lhs(ij,(ip-1)*nmlat_T1+jN+1) = coef(8,jN,i)
      if (i == 1) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN+1, &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1,(ip-1)*nmlat_T1+jN+1, &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1,(im-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(7,jS,i),coef(9,jS,i),coef(3,jS,i),coef(7,jN,i), &
          coef(8,jS,i),coef(1,jS,i),coef(2,jS,i),coef(8,jN,i), &
          coef(6,jS,i),coef(5,jS,i),coef(4,jS,i),coef(6,jN,i)/)
      elseif (i == nmlon) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1,(ip-1)*nmlat_T1+jN+1, &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1,(im-1)*nmlat_T1+jN+1, &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(8,jS,i),coef(1,jS,i),coef(2,jS,i),coef(8,jN,i), &
          coef(6,jS,i),coef(5,jS,i),coef(4,jS,i),coef(6,jN,i), &
          coef(7,jS,i),coef(9,jS,i),coef(3,jS,i),coef(7,jN,i)/)
      else
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1,(im-1)*nmlat_T1+jN+1, &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1,(i-1)*nmlat_T1+jN+1, &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1,(ip-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(6,jS,i),coef(5,jS,i),coef(4,jS,i),coef(6,jN,i), &
          coef(7,jS,i),coef(9,jS,i),coef(3,jS,i),coef(7,jN,i), &
          coef(8,jS,i),coef(1,jS,i),coef(2,jS,i),coef(8,jN,i)/)
      endif
      rowcnt(ij) = rowcnt(ij)+12

! note that coef has the direction switched in the NH
! SH coef does not have a direction switch
! therefore coef at j-1 is still coef at j'-1
      ij = (i-1)*nmlat_T1+jN
!      lhs(ij,(im-1)*nmlat_T1+jS-1) = coef(6,jS,i)
!      lhs(ij,(im-1)*nmlat_T1+jN-1) = coef(4,jN,i)
!      lhs(ij,(im-1)*nmlat_T1+jN  ) = coef(5,jN,i)
!      lhs(ij,(im-1)*nmlat_T1+jN+1) = coef(6,jN,i)
!      lhs(ij, (i-1)*nmlat_T1+jS-1) = coef(7,jS,i)
!      lhs(ij, (i-1)*nmlat_T1+jN-1) = coef(3,jN,i)
!      lhs(ij, (i-1)*nmlat_T1+jN  ) = coef(9,jN,i)
!      lhs(ij, (i-1)*nmlat_T1+jN+1) = coef(7,jN,i)
!      lhs(ij,(ip-1)*nmlat_T1+jS-1) = coef(8,jS,i)
!      lhs(ij,(ip-1)*nmlat_T1+jN-1) = coef(2,jN,i)
!      lhs(ij,(ip-1)*nmlat_T1+jN  ) = coef(1,jN,i)
!      lhs(ij,(ip-1)*nmlat_T1+jN+1) = coef(8,jN,i)
      if (i == 1) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(7,jS,i),coef(3,jN,i),coef(9,jN,i),coef(7,jN,i), &
          coef(8,jS,i),coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
          coef(6,jS,i),coef(4,jN,i),coef(5,jN,i),coef(6,jN,i)/)
      elseif (i == nmlon) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(8,jS,i),coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
          coef(6,jS,i),coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
          coef(7,jS,i),coef(3,jN,i),coef(9,jN,i),coef(7,jN,i)/)
      else
        jcol(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
          (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
          (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+12,ij) = (/ &
          coef(6,jS,i),coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
          coef(7,jS,i),coef(3,jN,i),coef(9,jN,i),coef(7,jN,i), &
          coef(8,jS,i),coef(2,jN,i),coef(1,jN,i),coef(8,jN,i)/)
      endif
      rowcnt(ij) = rowcnt(ij)+12

! from latm_JT to equator, symmetric solution
      do j = jlatm_JT+1,nmlat_h-1
        jS = j
        jN = nmlat_T1-j+1

        ij = (i-1)*nmlat_T1+jS
!        lhs(ij,(im-1)*nmlat_T1+jS-1) = coef(6,jS,i)
!        lhs(ij,(im-1)*nmlat_T1+jS  ) = coef(5,jS,i)
!        lhs(ij,(im-1)*nmlat_T1+jS+1) = coef(4,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS-1) = coef(7,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS  ) = coef(9,jS,i)
!        lhs(ij, (i-1)*nmlat_T1+jS+1) = coef(3,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS-1) = coef(8,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS  ) = coef(1,jS,i)
!        lhs(ij,(ip-1)*nmlat_T1+jS+1) = coef(2,jS,i)
        if (i == 1) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1, &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1, &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(7,jS,i),coef(9,jS,i),coef(3,jS,i), &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i), &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i)/)
        elseif (i == nmlon) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1, &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1, &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i), &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i), &
            coef(7,jS,i),coef(9,jS,i),coef(3,jS,i)/)
        else
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (im-1)*nmlat_T1+jS-1,(im-1)*nmlat_T1+jS,(im-1)*nmlat_T1+jS+1, &
            (i-1)*nmlat_T1+jS-1,(i-1)*nmlat_T1+jS,(i-1)*nmlat_T1+jS+1, &
            (ip-1)*nmlat_T1+jS-1,(ip-1)*nmlat_T1+jS,(ip-1)*nmlat_T1+jS+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(6,jS,i),coef(5,jS,i),coef(4,jS,i), &
            coef(7,jS,i),coef(9,jS,i),coef(3,jS,i), &
            coef(8,jS,i),coef(1,jS,i),coef(2,jS,i)/)
        endif
        rowcnt(ij) = rowcnt(ij)+9

! note that coef has the direction switched in the NH
        ij = (i-1)*nmlat_T1+jN
!        lhs(ij,(im-1)*nmlat_T1+jN-1) = coef(4,jN,i)
!        lhs(ij,(im-1)*nmlat_T1+jN  ) = coef(5,jN,i)
!        lhs(ij,(im-1)*nmlat_T1+jN+1) = coef(6,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jN-1) = coef(3,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jN  ) = coef(9,jN,i)
!        lhs(ij, (i-1)*nmlat_T1+jN+1) = coef(7,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN-1) = coef(2,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN  ) = coef(1,jN,i)
!        lhs(ij,(ip-1)*nmlat_T1+jN+1) = coef(8,jN,i)
        if (i == 1) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(3,jN,i),coef(9,jN,i),coef(7,jN,i), &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i)/)
        elseif (i == nmlon) then
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1, &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
            (i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i), &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
            coef(3,jN,i),coef(9,jN,i),coef(7,jN,i)/)
        else
          jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            (im-1)*nmlat_T1+jN-1,(im-1)*nmlat_T1+jN,(im-1)*nmlat_T1+jN+1, &
            (i-1)*nmlat_T1+jN-1,(i-1)*nmlat_T1+jN,(i-1)*nmlat_T1+jN+1, &
            (ip-1)*nmlat_T1+jN-1,(ip-1)*nmlat_T1+jN,(ip-1)*nmlat_T1+jN+1/)
          nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
            coef(4,jN,i),coef(5,jN,i),coef(6,jN,i), &
            coef(3,jN,i),coef(9,jN,i),coef(7,jN,i), &
            coef(2,jN,i),coef(1,jN,i),coef(8,jN,i)/)
        endif
        rowcnt(ij) = rowcnt(ij)+9
      enddo

! equator
      j = nmlat_h

      ij = (i-1)*nmlat_T1+j
!      lhs(ij,(im-1)*nmlat_T1+j-1) = coef(6,j,i)
!      lhs(ij,(im-1)*nmlat_T1+j  ) = coef(5,j,i)
!      lhs(ij,(im-1)*nmlat_T1+j+1) = coef(4,j,i)
!      lhs(ij, (i-1)*nmlat_T1+j-1) = coef(7,j,i)
!      lhs(ij, (i-1)*nmlat_T1+j  ) = coef(9,j,i)
!      lhs(ij, (i-1)*nmlat_T1+j+1) = coef(3,j,i)
!      lhs(ij,(ip-1)*nmlat_T1+j-1) = coef(8,j,i)
!      lhs(ij,(ip-1)*nmlat_T1+j  ) = coef(1,j,i)
!      lhs(ij,(ip-1)*nmlat_T1+j+1) = coef(2,j,i)
      if (i == 1) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          (i-1)*nmlat_T1+j-1,(i-1)*nmlat_T1+j,(i-1)*nmlat_T1+j+1, &
          (ip-1)*nmlat_T1+j-1,(ip-1)*nmlat_T1+j,(ip-1)*nmlat_T1+j+1, &
          (im-1)*nmlat_T1+j-1,(im-1)*nmlat_T1+j,(im-1)*nmlat_T1+j+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          coef(7,j,i),coef(9,j,i),coef(3,j,i), &
          coef(8,j,i),coef(1,j,i),coef(2,j,i), &
          coef(6,j,i),coef(5,j,i),coef(4,j,i)/)
      elseif (i == nmlon) then
        jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          (ip-1)*nmlat_T1+j-1,(ip-1)*nmlat_T1+j,(ip-1)*nmlat_T1+j+1, &
          (im-1)*nmlat_T1+j-1,(im-1)*nmlat_T1+j,(im-1)*nmlat_T1+j+1, &
          (i-1)*nmlat_T1+j-1,(i-1)*nmlat_T1+j,(i-1)*nmlat_T1+j+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          coef(8,j,i),coef(1,j,i),coef(2,j,i), &
          coef(6,j,i),coef(5,j,i),coef(4,j,i), &
          coef(7,j,i),coef(9,j,i),coef(3,j,i)/)
      else
        jcol(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          (im-1)*nmlat_T1+j-1,(im-1)*nmlat_T1+j,(im-1)*nmlat_T1+j+1, &
          (i-1)*nmlat_T1+j-1,(i-1)*nmlat_T1+j,(i-1)*nmlat_T1+j+1, &
          (ip-1)*nmlat_T1+j-1,(ip-1)*nmlat_T1+j,(ip-1)*nmlat_T1+j+1/)
        nzval(rowcnt(ij)+1:rowcnt(ij)+9,ij) = (/ &
          coef(6,j,i),coef(5,j,i),coef(4,j,i), &
          coef(7,j,i),coef(9,j,i),coef(3,j,i), &
          coef(8,j,i),coef(1,j,i),coef(2,j,i)/)
      endif
      rowcnt(ij) = rowcnt(ij)+9
    enddo

    rowptr(1) = 1
    do i = 2,nlonlat+1
      rowptr(i) = rowptr(i-1)+rowcnt(i-1)
    enddo

    do j = 1,nmlon+2
      colind(j) = jcol1(j)
      values(j) = nzval1(j)
    enddo

    do i = 2,nlonlat
      do j = 1,rowcnt(i)
        colind(rowptr(i)+j-1) = jcol(j,i)
        values(rowptr(i)+j-1) = nzval(j,i)
      enddo
    enddo

  endsubroutine construct_lhs
!-----------------------------------------------------------------------
  pure function construct_rhs(coef_10) result(rhs)
! construct vector RHS
! this is different from flatten

    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:phi_pol

    real(kind=rp),dimension(nmlat_T1,nmlon),intent(in) :: coef_10
    real(kind=rp),dimension(nmlat_T1*nmlon) :: rhs

    integer :: i,j,ij

    rhs = 0

! set up poles, there are no c6,c7,c8 values
    j = 1

! for longitude i=1 at the south pole
    i = 1
    ij = (i-1)*nmlat_T1+j
    rhs(ij) = sum(coef_10(j,:))

! for each i, set Phi(i,nmlat_T1) = phi_pol at the north pole
    do concurrent (i = 1:nmlon)
      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      rhs(ij) = phi_pol
    enddo

    do concurrent (i = 1:nmlon, j = 2:nmlat_h)
      ij = (i-1)*nmlat_T1+j
      rhs(ij) = coef_10(j,i)

      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      rhs(ij) = coef_10(nmlat_T1-j+1,i)
    enddo

  endfunction construct_rhs
!-----------------------------------------------------------------------
  function solve_mkl(n,nnz,rowptr,colind,values,rhs) result(sol)

#ifdef MKL
    include 'mkl_pardiso.fi'
#endif

    integer,intent(in) :: n,nnz
    integer,dimension(n+1),intent(in) :: rowptr
    integer,dimension(nnz),intent(in) :: colind
    real(kind=rp),dimension(nnz),intent(in) :: values
    real(kind=rp),dimension(n),intent(in) :: rhs
    real(kind=rp),dimension(n) :: sol

#ifdef MKL
! for PARDISO sparse matrix solver
    logical,parameter :: debug = .true.
    integer,parameter :: maxfct = 1, mnum = 1, nrhs = 1, &
      mtype = 11, & ! real and non-symmetric matrix
      msglvl = 1 ! print statistical information
    integer :: i,phase,error
    integer,dimension(64) :: iparm
    integer,dimension(n) :: perm
    real(kind=rp),dimension(n) :: rhs_cp
    type(MKL_PARDISO_HANDLE),dimension(64) :: pt

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
    do concurrent (i = 1:n)
      rhs_cp(i) = rhs(i)
    enddo

! use 32-bit integer version
! if the number of non-zero elements is on the order of 500 million or more
! then use pardiso_64 (64-bit integer version)
    call pardiso(pt, maxfct, mnum, mtype, phase, n, &
      values, rowptr, colind, perm, nrhs, iparm, msglvl, rhs_cp, sol, error)

    write(6,"('phase ',i4,' error ',i4)") phase,error
#else
    sol = 0 ! to suppress the "return value not set" warning
#endif

  endfunction solve_mkl
!-----------------------------------------------------------------------
  function solve_superlu(n,nnz,colptr,rowind,values,rhs) result(sol)
    use iso_c_binding,only:c_int,c_long_long,c_double

    integer,intent(in) :: n,nnz
    integer(kind=c_int),dimension(n+1),intent(in) :: colptr
    integer(kind=c_int),dimension(nnz),intent(in) :: rowind
    real(kind=c_double),dimension(nnz),intent(in) :: values
    real(kind=rp),dimension(n),intent(in) :: rhs
    real(kind=rp),dimension(n) :: sol

! for SuperLU sparse matrix solver
    integer,parameter :: nrhs = 1
    integer :: i,iopt,info
    integer(kind=c_long_long) :: f_factors

    interface
      subroutine c_fortran_dgssv(iopt,n,nnz,nrhs, &
        values,rowind,colptr,b,ldb,f_factors,info) &
        bind(c,name='c_fortran_dgssv_')
        use iso_c_binding,only:c_int,c_long_long,c_double
        integer(kind=c_int) :: iopt,n,nnz,nrhs,ldb,info
        real(kind=c_double),dimension(nnz) :: values
        integer(kind=c_int),dimension(nnz) :: rowind
        integer(kind=c_int),dimension(n+1) :: colptr
        real(kind=c_double),dimension(ldb) :: b
        integer(kind=c_long_long) :: f_factors
      endsubroutine c_fortran_dgssv
    endinterface

    do concurrent (i = 1:n)
      sol(i) = rhs(i)
    enddo

! first, factorize the matrix, the factors are stored in *f_factors* handle
    iopt = 1
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      values, rowind, colptr, sol, n, f_factors, info)
    write(6,"('INFO from LU decomposition = ',i4)") info

! second, solve the system using the existing factors
    iopt = 2
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      values, rowind, colptr, sol, n, f_factors, info)
    write(6,"('INFO from triangular solve = ',i4)") info

! last, free the storage allocated inside SuperLU
    iopt = 3
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      values, rowind, colptr, sol, n, f_factors, info)

  endfunction solve_superlu
!-----------------------------------------------------------------------
  pure subroutine csr_to_csc(nrow,ncol,nnz, &
    rowptr,colind,values_csr,colptr,rowind,values_csc)

    integer,intent(in) :: nrow,ncol,nnz
    integer,dimension(nrow+1),intent(in) :: rowptr
    integer,dimension(nnz),intent(in) :: colind
    real(kind=rp),dimension(nnz),intent(in) :: values_csr
    integer,dimension(ncol+1),intent(out) :: colptr
    integer,dimension(nnz),intent(out) :: rowind
    real(kind=rp),dimension(nnz),intent(out) :: values_csc

    integer :: n,i,j,newidx
    integer,dimension(ncol) :: colcnt,cnt

    colcnt = 0
    do n = 1,nnz
      j = colind(n)
      colcnt(j) = colcnt(j)+1
    enddo

    colptr(1) = 1
    do j = 2,ncol+1
      colptr(j) = colptr(j-1)+colcnt(j-1)
    enddo

    cnt = 0
    do i = 1,nrow
      do n = rowptr(i),rowptr(i+1)-1
        j = colind(n)
        newidx = colptr(j)+cnt(j)
        rowind(newidx) = i
        values_csc(newidx) = values_csr(n)
        cnt(j) = cnt(j)+1
      enddo
    enddo

  endsubroutine csr_to_csc
!-----------------------------------------------------------------------
  pure function flatten(fin) result(fout)
! reorder 2D fields (lat-lon) into 1D vector (RHS)
! northern/southern hemispheres are either separate or averaged
! based on their latitude ranges (high-lat, transition, low-lat, equator)

    use params_module,only:nmlat_h,nmlat_T1,nmlon
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(2,nmlat_h,nmlon),intent(in) :: fin
    real(kind=rp),dimension(nmlat_T1*nmlon) :: fout

    integer :: i,j,ij
    real(kind=rp) :: avg

! from pole to latm_JT, two hemispheres are uncoupled
    do concurrent (i = 1:nmlon, j = 1:jlatm_JT)
      ij = (i-1)*nmlat_T1+j
      fout(ij) = fin(1,j,i)

      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      fout(ij) = fin(2,j,i)
    enddo

! from latm_JT to equator, symmetric solution
    do concurrent (i = 1:nmlon, j = jlatm_JT+1:nmlat_h)
      avg = (fin(1,j,i)+fin(2,j,i))/2

      ij = (i-1)*nmlat_T1+j
      fout(ij) = avg

      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      fout(ij) = avg
    enddo

  endfunction flatten
!-----------------------------------------------------------------------
  pure function unravel(fin) result(fout)
! reorder 1D vector (RHS) into 2D fields (lat-lon)

    use params_module,only:nmlat_h,nmlat_T1,nmlon

    real(kind=rp),dimension(nmlat_T1*nmlon),intent(in) :: fin
    real(kind=rp),dimension(2,nmlat_h,nmlon) :: fout

    integer :: i,j,isn,ij

    do concurrent (i = 1:nmlon, j = 1:nmlat_h, isn = 1:2)
      if (isn == 1) then
        ij = (i-1)*nmlat_T1+j
      else
        ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      endif
      fout(isn,j,i) = fin(ij)
    enddo

  endfunction unravel
!-----------------------------------------------------------------------
endmodule solver_module
