module edyn3D_solver_module

  use shr_kind_mod,  only: r8 => shr_kind_r8            ! 8-byte reals
  use edyn3D_params, only: nmlat_h, nmlat_T1, nmlon, jlatm_JT, phi_pol, nlonlat=>nlonlat_T1
  use perf_mod, only: t_startf, t_stopf

  implicit none

  private
  public :: linear_system

  contains
!-----------------------------------------------------------------------
  subroutine linear_system(bij_full,pot_hl_full,fac_hl_full,coef_ns_full,pot_full)
! construct linear system based on bij and coef and solve in pot

! if FAC is read in, pot_hl is not used, only fac_hl is used
! if potential is read in, pot_hl is used, fac_hl is output

    real(r8),dimension(nmlat_h,nmlon),intent(in) :: bij_full
    real(r8),dimension(2,nmlat_h,nmlon),intent(in) :: pot_hl_full
    real(r8),dimension(2,nmlat_h,nmlon),intent(out) :: fac_hl_full
    real(r8),dimension(10,nmlat_T1,nmlon),intent(inout) :: coef_ns_full
    real(r8),dimension(2,nmlat_h,nmlon),intent(out) :: pot_full

    logical,parameter :: use_mkl = &
#ifdef MKL
      .true.
#else
      .false.
#endif
    integer,parameter :: root = 0
    integer :: i,j,isn,nnz
    integer,dimension(nlonlat+1) :: rowptr,colptr
    integer,dimension(12*nlonlat) :: colind,rowind

    real(r8),dimension(12*nlonlat) :: values_csr,values_csc

    real(r8),dimension(nlonlat) :: rhs,z,pot_hl_f,sol

    call t_startf('linear_system')

    ! Q from Wu: why coefficients at the south pole are gathered to i=1?
    j = 1
    coef_ns_full(9,j,1) = sum(coef_ns_full(9,j,:))
    coef_ns_full(10,j,1) = sum(coef_ns_full(10,j,:))
    do concurrent (i = 2:nmlon)
       coef_ns_full(9,j,i) = 0
       coef_ns_full(10,j,i) = 0
    enddo

    call t_startf('linear_system->construct_lhs')
    ! construct LHS matrix in CSR format
    call construct_lhs(bij_full,coef_ns_full(1:9,:,:),rowptr,colind,values_csr)
    nnz = rowptr(nlonlat+1)-1
    call t_stopf('linear_system->construct_lhs')

    call t_startf('linear_system->construct_rhs')
    ! RHS is dense
    rhs = construct_rhs(coef_ns_full(10,:,:))
    call t_stopf('linear_system->construct_rhs')

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
    fac_hl_full(:,:,1:nmlon) = unravel(z)

    ! add FAC forcing to RHS
    do i = 1,nlonlat
       rhs(i) = rhs(i)+z(i)
    enddo

    call t_startf('linear_system->solve')
    if (use_mkl) then
       sol = solve_mkl(nlonlat,nnz,rowptr,colind(1:nnz),values_csr(1:nnz),rhs)
    else
       call t_startf('linear_system->csr_to_csc')
       call csr_to_csc(nlonlat,nlonlat,nnz, &
            rowptr,colind(1:nnz),values_csr(1:nnz), &
            colptr,rowind(1:nnz),values_csc(1:nnz))
       call t_stopf('linear_system->csr_to_csc')
       sol = solve_superlu(nlonlat,nnz,colptr,rowind(1:nnz),values_csc(1:nnz),rhs)
    endif
    call t_stopf('linear_system->solve')

    ! reconstruct 2D distribution of potential based on the solution
    pot_full(:,:,1:nmlon) = unravel(sol)

    call t_stopf('linear_system')

  end subroutine linear_system

!-----------------------------------------------------------------------
  subroutine construct_lhs(bij,coef,rowptr,colind,values)
! construct LHS matrix (CSR format)

! need to set where the two hemispheres are connected
! This is not in the 9-point stencil but needs to be done manually

    real(r8),dimension(nmlat_h,nmlon),intent(in) :: bij
    real(r8),dimension(9,nmlat_T1,nmlon),intent(in) :: coef
    integer,dimension(nlonlat+1),intent(out) :: rowptr
    integer,dimension(12*nlonlat),intent(out) :: colind
    real(r8),dimension(12*nlonlat),intent(out) :: values

! if two hemispheres are uncoupled at high latitudes, set bijSum to zero
    real(r8),parameter :: bijSum = 0
    integer :: i,j,jS,jN,ij,isub,im,ip
    integer,dimension(nlonlat) :: rowcnt

! the first row has most elements (nmlon+2)
    integer,dimension(nmlon+2) :: jcol1,sorted,idx
    real(r8),dimension(nmlon+2) :: nzval1

! other rows have at most 12 elements
    integer,dimension(12,nlonlat) :: jcol
    real(r8),dimension(12,nlonlat) :: nzval

    rowcnt = 0

! set up poles
    j = 1

! there are no c6,c7,c8 values at the south pole

! for longitude i=1 at the south pole
    i = 1
!    ij = (i-1)*nmlat_T1+j
!    lhs(ij,(i-1)*nmlat_T1+j           ) = coef(9,j,i)-bijSum ! A. Richmond 2023/06/20: add bij effects
!    lhs(ij,(i-1)*nmlat_T1+j+1         ) = coef(3,j,i)        ! Sum_i=1^nmlon C3(i,1) Phi(i,2)
!    lhs(ij,(i-1)*nmlat_T1+nmlat_T1-j+1) = bijSum             ! conjugate point (north pole at i=1)
    jcol1(1:3) = (/(i-1)*nmlat_T1+j,(i-1)*nmlat_T1+j+1,(i-1)*nmlat_T1+nmlat_T1-j+1/)
    nzval1(1:3) = (/coef(9,j,i)-bijSum,coef(3,j,i),bijSum/)
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

    if (inorder(jcol1)) then
      do j = 1,nmlon+2
        colind(j) = jcol1(j)
        values(j) = nzval1(j)
      enddo
    else
      !write(6,"('Row 1 is not in order, sorting')")
      call argsort(nmlon+2,nlonlat,jcol1,sorted,idx)
      do j = 1,nmlon+2
        colind(j) = sorted(j)
        values(j) = nzval1(idx(j))
      enddo
    endif

    do i = 2,nlonlat
      if (inorder(jcol(1:rowcnt(i),i))) then
        do j = 1,rowcnt(i)
          colind(rowptr(i)+j-1) = jcol(j,i)
          values(rowptr(i)+j-1) = nzval(j,i)
        enddo
      else
        !write(6,"('Row ',i6,' is not in order, sorting')") i
        call argsort(rowcnt(i),nlonlat,jcol(1:rowcnt(i),i), &
          sorted(1:rowcnt(i)),idx(1:rowcnt(i)))
        do j = 1,rowcnt(i)
          colind(rowptr(i)+j-1) = sorted(j)
          values(rowptr(i)+j-1) = nzval(idx(j),i)
        enddo
      endif
    enddo

  end subroutine construct_lhs
!-----------------------------------------------------------------------
  pure function construct_rhs(coef_10) result(rhs)
! construct vector RHS
! this is different from flatten

    real(r8),dimension(nmlat_T1,nmlon),intent(in) :: coef_10
    real(r8),dimension(nlonlat) :: rhs

    integer :: i,j,ij

    rhs = 0

! set up poles, there are no c6,c7,c8 values
    j = 1

! for longitude i=1 at the south pole
    i = 1
    ij = (i-1)*nmlat_T1+j
    rhs(ij) = coef_10(j,i)

! for each i, set Phi(i,nmlat_T1) = phi_pol at the north pole
    do concurrent (i = 1:nmlon)
      ij = (i-1)*nmlat_T1+nmlat_T1-j+1
      rhs(ij) = phi_pol
    enddo

    do concurrent (i = 1:nmlon, j = 2:nmlat_h)
      ij = (i-1)*nmlat_T1+j
      rhs(ij) = coef_10(j,i)

      if (j /= nmlat_h) then
        ij = (i-1)*nmlat_T1+nmlat_T1-j+1
        rhs(ij) = coef_10(nmlat_T1-j+1,i)
      endif
    enddo

  end function construct_rhs
!-----------------------------------------------------------------------
  function solve_mkl(n,nnz,rowptr,colind,values,rhs) result(sol)

#ifdef MKL
    include 'mkl_pardiso.fi'
#endif

    integer,intent(in) :: n,nnz
    integer,dimension(n+1),intent(in) :: rowptr
    integer,dimension(nnz),intent(in) :: colind
    real(r8),dimension(nnz),intent(in) :: values
    real(r8),dimension(n),intent(in) :: rhs
    real(r8),dimension(n) :: sol

#ifdef MKL
! for PARDISO sparse matrix solver
    logical,parameter :: debug = .true.
    integer,parameter :: maxfct = 1, mnum = 1, nrhs = 1, &
      mtype = 11, & ! real and non-symmetric matrix
      msglvl = 1 ! print statistical information
    integer :: i,phase,error
    integer,dimension(64) :: iparm
    integer,dimension(n) :: perm
    real(r8),dimension(n) :: rhs_cp
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

  end function solve_mkl
!-----------------------------------------------------------------------
  function solve_superlu(n,nnz,colptr,rowind,values,rhs) result(sol)
    use iso_c_binding, only: c_long_long

    integer,intent(in) :: n,nnz
    integer,dimension(n+1),intent(in) :: colptr
    integer,dimension(nnz),intent(in) :: rowind
    real(r8),dimension(nnz),intent(in) :: values
    real(r8),dimension(n),intent(in) :: rhs
    real(r8),dimension(n) :: sol

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
    !write(6,"('INFO from LU decomposition = ',i4)") info

! second, solve the system using the existing factors
    iopt = 2
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      values, rowind, colptr, sol, n, f_factors, info)
    !write(6,"('INFO from triangular solve = ',i4)") info

! last, free the storage allocated inside SuperLU
    iopt = 3
    call c_fortran_dgssv(iopt, n, nnz, nrhs, &
      values, rowind, colptr, sol, n, f_factors, info)

  end function solve_superlu
!-----------------------------------------------------------------------
  subroutine csr_to_csc(nrow,ncol,nnz, &
    rowptr,colind,values_csr,colptr,rowind,values_csc)

    integer,intent(in) :: nrow,ncol,nnz
    integer,dimension(nrow+1),intent(in) :: rowptr
    integer,dimension(nnz),intent(in) :: colind
    real(r8),dimension(nnz),intent(in) :: values_csr
    integer,dimension(ncol+1),intent(out) :: colptr
    integer,dimension(nnz),intent(out) :: rowind
    real(r8),dimension(nnz),intent(out) :: values_csc

    integer :: n,i,j, ndx
    integer,dimension(ncol) :: colcnt,cnt
    integer,dimension(nrow) :: sorted,idx
    real(r8),dimension(nrow) :: values

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
        ndx = colptr(j)+cnt(j)
        rowind(ndx) = i
        values_csc(ndx) = values_csr(n)
        cnt(j) = cnt(j)+1
      enddo
    enddo

    do j = 1,ncol
      if (.not. inorder(rowind(colptr(j):colptr(j+1)-1))) then
        !write(6,"('Column ',i6,' is not in order, sorting')")
        n = colptr(j+1)-colptr(j)
        call argsort(n,nrow,rowind(colptr(j):colptr(j+1)-1),sorted(1:n),idx(1:n))
        values(1:n) = values_csc(colptr(j):colptr(j+1)-1)
        do i = 1,n
          rowind(colptr(j)+i-1) = sorted(i)
          values_csc(colptr(j)+i-1) = values(idx(i))
        enddo
      endif
    enddo

  end subroutine csr_to_csc
!-----------------------------------------------------------------------
  pure function flatten(fin) result(fout)
! reorder 2D fields (lat-lon) into 1D vector (RHS)
! northern/southern hemispheres are either separate or averaged
! based on their latitude ranges (high-lat, transition, low-lat, equator)

    real(r8),dimension(2,nmlat_h,nmlon),intent(in) :: fin
    real(r8),dimension(nlonlat) :: fout

    integer :: i,j,ij

    do concurrent (i = 1:nmlon, j = 1:nmlat_h)
! j should take one of the following three situations

! from pole to latm_JT, two hemispheres are uncoupled
      if (j <= jlatm_JT) then
        ij = (i-1)*nmlat_T1+j
        fout(ij) = fin(1,j,i)

        ij = (i-1)*nmlat_T1+nmlat_T1-j+1
        fout(ij) = fin(2,j,i)
      endif

! from latm_JT to equator, symmetric solution
      if (j>=jlatm_JT+1 .and. j<=nmlat_h-1) then
        ij = (i-1)*nmlat_T1+j
        fout(ij) = (fin(1,j,i)+fin(2,j,i))/2

        ij = (i-1)*nmlat_T1+nmlat_T1-j+1
        fout(ij) = (fin(1,j,i)+fin(2,j,i))/2
      endif

! equator
      if (j == nmlat_h) then
        ij = (i-1)*nmlat_T1+j
        fout(ij) = fin(1,j,i)
!        fout(ij) = (fin(1,j,i)+fin(2,j,i))/2
      endif
    enddo

  end function flatten
!-----------------------------------------------------------------------
  pure function unravel(fin) result(fout)
! reorder 1D vector (RHS) into 2D fields (lat-lon)

    real(r8),dimension(nlonlat),intent(in) :: fin
    real(r8),dimension(2,nmlat_h,nmlon) :: fout

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
  pure function inorder(a) result(flag)

    integer,dimension(:),intent(in) :: a
    logical :: flag

    integer :: i

    flag = .true.
    do i = 1,size(a)-1
      if (a(i) > a(i+1)) then
        flag = .false.
        return
      endif
    enddo

  end function inorder
!-----------------------------------------------------------------------
  pure subroutine argsort(n,amax,a,sorted,idx)
! radix sort for row or column indices

    integer,intent(in) :: n,amax
    integer,dimension(n),intent(in) :: a
    integer,dimension(n),intent(out) :: sorted,idx

    integer :: i,e

    do i = 1,n
      sorted(i) = a(i)
      idx(i) = i
    enddo

    e = 1
    do while (e < amax)
      call counting_sort(n,e,sorted,idx)
      e = e*10
    enddo

  end subroutine argsort
!-----------------------------------------------------------------------
  pure subroutine counting_sort(n,e,a,idx)
! counting sort for individual digits

    integer,intent(in) :: n,e
    integer,dimension(n),intent(inout) :: a,idx

    integer :: i,digit,ind
    integer,dimension(0:9) :: cnt
    integer,dimension(n) :: sorted,sort_idx

    cnt = 0

    do i = 1,n
      digit = modulo(a(i)/e,10)
      cnt(digit) = cnt(digit)+1
    enddo

    do i = 1,9
      cnt(i) = cnt(i)+cnt(i-1)
    enddo

    do i = n,1,-1
      digit = modulo(a(i)/e,10)
      ind = cnt(digit)
      sorted(ind) = a(i)
      sort_idx(ind) = idx(i)
      cnt(digit) = cnt(digit)-1
    enddo

    do i = 1,n
      a(i) = sorted(i)
      idx(i) = sort_idx(i)
    enddo

  end subroutine counting_sort
!-----------------------------------------------------------------------
end module edyn3D_solver_module
