module solver_module

  use prec,only:rp

  implicit none

  real(kind=rp),parameter :: &
    phi_np = 0, & ! north pole potential
    beta = 0 ! if two hemispheres are uncoupled at high latitudes, set beta to zero
  integer :: nlonlat ! total number of grids to be solved

! convention for grid numbers: longitude increases first then latitude
! so external loop is on latitudes, internal loop is on longitudes

! south pole:
! isn = 1, j = 1, i = 1; then ij = 1

! southern latitudes (no pole):
! isn = 1, 2 <= j <= nmlat_h; then ij = (j-2)*nmlon+i+1

! northern high latitudes (no pole):
! isn = 2, 2 <= j <= jlatm_JT-1; then ij = (nmlat_h+jlatm_JT-2-j)*nmlon+i+1

  contains
!-----------------------------------------------------------------------
  subroutine linear_system(mlatd0,mlatd1,mlond0,mlond1,bij,pot_hl,fac_hl,src,coef,pot)
! construct linear system based on src and coef and solve in pot

! if potential is read in, pot_hl is used, fac_hl is output
! if FAC is read in, pot_hl is not used, only fac_hl is used

    use params_module,only:nmlat_h,nmlon
    use params_module,only:read_pot,read_fac
    use mpi_module,only:gather_mag,bcast,mpi_rank

    integer,intent(in) :: mlatd0,mlatd1,mlond0,mlond1
    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_hl,src
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(inout) :: fac_hl
    real(kind=rp),dimension(9,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: pot

    integer,parameter :: root = 0
    integer :: mlat0,mlat1,mlon0,mlon1,i,j,isn,nnz
    integer,dimension(nlonlat+1) :: rowptr,colptr
    integer,dimension(12*nlonlat) :: colind,rowind
    real(kind=rp) :: offset
    real(kind=rp),dimension(nlonlat) :: rhs,z,pot_hl_f,sol
    real(kind=rp),dimension(12*nlonlat) :: values_csr,values_csc
    real(kind=rp),dimension(nmlat_h,nmlon) :: bij_full
    real(kind=rp),dimension(2,nmlat_h,nmlon) :: pot_hl_full,fac_hl_full,src_full
    real(kind=rp),dimension(9,2,nmlat_h,nmlon) :: coef_full
    real(kind=rp),dimension(2,nmlat_h,0:nmlon+1) :: fac_hl_2,pot_2

    mlat0 = mlatd0+1
    mlat1 = mlatd1-1
    mlon0 = mlond0+1
    mlon1 = mlond1-1

    bij_full = gather_mag(bij(mlat0:mlat1,mlon0:mlon1),root)
    if (read_pot) pot_hl_full = gather_mag(pot_hl(:,mlat0:mlat1,mlon0:mlon1),2,root)
    if (read_fac) fac_hl_full = gather_mag(fac_hl(:,mlat0:mlat1,mlon0:mlon1),2,root)
    src_full = gather_mag(src(:,mlat0:mlat1,mlon0:mlon1),2,root)
    coef_full = gather_mag(coef(:,:,mlat0:mlat1,mlon0:mlon1),9,2,root)

    if (mpi_rank == root) then

! construct LHS matrix in CSR format
      call construct_lhs(bij_full,coef_full,rowptr,colind,values_csr)
      nnz = rowptr(nlonlat+1)-1

! RHS is vector (dense)
      rhs = construct_rhs(src_full, &
        coef_full(6,2,2,:)+coef_full(7,2,2,:)+coef_full(8,2,2,:))

! determine FAC forcing (dense)

! input is pot_hl, fac_hl is to be calculated (output)
      if (read_pot) then

! move north pole potential to phi_np (this might not be necessary)
        offset = phi_np-pot_hl_full(2,1,1)
        do i = 1,nmlon
          do j = 1,nmlat_h
            do isn = 1,2
              pot_hl_full(isn,j,i) = pot_hl_full(isn,j,i)+offset
            enddo
          enddo
        enddo

! A. Maute 2023/11/21: put the high latitude potential in X
! and then use LHS to calculate the RHS FAC
        pot_hl_f = ravel(pot_hl_full)

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
        fac_hl_2(2,1,:) = 0 ! north pole is not set in unravel

! add periodic points
        do j = 1,nmlat_h
          do isn = 1,2
            fac_hl_2(isn,j,0) = fac_hl_2(isn,j,nmlon)
            fac_hl_2(isn,j,nmlon+1) = fac_hl_2(isn,j,1)
          enddo
        enddo
      endif

! input is corrected fac_hl, pot_hl is not used
      if (read_fac) z = ravel(fac_hl_full)

! add FAC forcing to RHS
      do i = 1,nlonlat
        rhs(i) = rhs(i)+z(i)
      enddo

#ifdef USE_MKL
      sol = solve_mkl(nlonlat,nnz,rowptr,colind(1:nnz),values_csr(1:nnz),rhs)
#else
      call csr_to_csc(nlonlat,nlonlat,nnz, &
        rowptr,colind(1:nnz),values_csr(1:nnz), &
        colptr,rowind(1:nnz),values_csc(1:nnz))

      sol = solve_superlu(nlonlat,nnz,colptr,rowind(1:nnz),values_csc(1:nnz),rhs)

#endif

! reconstruct 2D distribution of potential based on the solution
      pot_2(:,:,1:nmlon) = unravel(sol)
      pot_2(2,1,:) = phi_np ! north pole is not set in unravel

! periodic points
      do j = 1,nmlat_h
        do isn = 1,2
          pot_2(isn,j,0) = pot_2(isn,j,nmlon)
          pot_2(isn,j,nmlon+1) = pot_2(isn,j,1)
        enddo
      enddo
    endif

    if (read_pot) then
      call bcast(fac_hl_2,2,nmlat_h,nmlon+2,root)

      do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
        fac_hl(isn,j,i) = fac_hl_2(isn,j,i)
      enddo
    endif

    call bcast(pot_2,2,nmlat_h,nmlon+2,root)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
      pot(isn,j,i) = pot_2(isn,j,i)
    enddo

  endsubroutine linear_system
!-----------------------------------------------------------------------
  pure subroutine construct_lhs(bij,coef,rowptr,colind,values)
! construct LHS matrix (CSR format)
! rearrange each element in row and then concatenate

! this is not a 9-point stencil

    use params_module,only:nmlat_h,nmlon
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(nmlat_h,nmlon),intent(in) :: bij
    real(kind=rp),dimension(9,2,nmlat_h,nmlon),intent(in) :: coef
    integer,dimension(nlonlat+1),intent(out) :: rowptr
    integer,dimension(12*nlonlat),intent(out) :: colind
    real(kind=rp),dimension(12*nlonlat),intent(out) :: values

    integer :: i,j,isn,ij,isub,im,ip
    real(kind=rp) :: c1,c2,c3,c4,c5,c6,c7,c8,c9

! each row corresponds to one grid in the lat-lon decomposition
! each grid has rowcnt(ij) non-zero elements
    integer,dimension(nlonlat) :: rowcnt

! the first row is longitude i=1 at south pole
! it has nmlon+1 elements
    integer,dimension(nmlon+1) :: jcol1
    real(kind=rp),dimension(nmlon+1) :: nzval1

! other rows are remaining grids except for longitude i=1 at south pole
! they have at most 12 elements
    integer,dimension(12,2:nlonlat) :: jcol
    real(kind=rp),dimension(12,2:nlonlat) :: nzval

! start with southern hemisphere
    isn = 1

! south pole has only one point (longitudes are squeezed)
! the volume is a combined cell of all longitudes (a prism with nmlon edges)
    j = 1
    i = 1 ! Equation (7.23) is applied to i=1
    ij = index_3to1(isn,j,i)
    rowcnt(ij) = nmlon+1

! A. Richmond 2023/06/20: add b effects
! (Sum_i C9S(i) - beta) PhiS(1,1)
    jcol1(1) = index_3to1(isn,j,i)
    nzval1(1) = sum(coef(9,isn,j,:))-beta

! Sum_i C3S(i) PhiS(2,i)
    do concurrent (isub = 1:nmlon)
      jcol1(isub+1) = index_3to1(isn,j+1,isub)
      nzval1(isub+1) = coef(3,isn,j,isub)
    enddo

! j=2 in the southern hemisphere is different from others
! south pole (j=1) has only one point, so coefficients are summed up
    j = 2
    do concurrent (i = 1:nmlon)
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
      ij = index_3to1(isn,j,i)
      rowcnt(ij) = 8

! C1S     PhiS(j  ,i+1)
! C2S     PhiS(j+1,i+1)
! C3S     PhiS(j+1,i  )
! C4S     PhiS(j+1,i-1)
! C5S     PhiS(j  ,i-1)
! (C9S-b) PhiS(j  ,i  )
! b       PhiN(j  ,i  )

! PhiS(j-1,i-1), PhiS(j-1,i), PhiS(j-1,i+1) are south pole
! so their coefficients C6S, C7S, C8S are summed up
! which becomes (C6S+C7S+C8S) PhiS(j-1,i)
      c1 = coef(1,isn,j,i)
      c2 = coef(2,isn,j,i)
      c3 = coef(3,isn,j,i)
      c4 = coef(4,isn,j,i)
      c5 = coef(5,isn,j,i)
      c7 = coef(6,isn,j,i)+coef(7,isn,j,i)+coef(8,isn,j,i)
      c9 = coef(9,isn,j,i)-bij(j,i)

      if (i == 1) then
        jcol(1:8,ij) = &
          (/index_3to1(isn,j-1,i), &
            index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im), &
            index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im), &
            index_3to1(2  ,j  ,i)/)
        nzval(1:8,ij) = &
          (/c7, &
            c9,c1,c5, &
            c3,c2,c4, &
            bij(j,i)/)
      elseif (i == nmlon) then
        jcol(1:8,ij) = &
          (/                                              index_3to1(isn,j-1,i), &
            index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im),index_3to1(isn,j  ,i), &
            index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im),index_3to1(isn,j+1,i), &
                                                          index_3to1(2  ,j  ,i)/)
        nzval(1:8,ij) = &
          (/      c7, &
            c1,c5,c9, &
            c2,c4,c3, &
            bij(j,i)/)
      else
        jcol(1:8,ij) = &
          (/                       index_3to1(isn,j-1,i), &
            index_3to1(isn,j  ,im),index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip), &
            index_3to1(isn,j+1,im),index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip), &
                                   index_3to1(2  ,j  ,i)/)
        nzval(1:8,ij) = &
          (/   c7, &
            c5,c9,c1, &
            c4,c3,c2, &
            bij(j,i)/)
      endif
    enddo

! southern high latitudes, Equation (7.22)
    do concurrent (i = 1:nmlon, j = 3:jlatm_JT-1)
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
      ij = index_3to1(isn,j,i)
      rowcnt(ij) = 10

! C1S     PhiS(j  ,i+1)
! C2S     PhiS(j+1,i+1)
! C3S     PhiS(j+1,i  )
! C4S     PhiS(j+1,i-1)
! C5S     PhiS(j  ,i-1)
! C6S     PhiS(j-1,i-1)
! C7S     PhiS(j-1,i  )
! C8S     PhiS(j-1,i+1)
! (C9S-b) PhiS(j  ,i  )
! b       PhiN(j  ,i  )
      c1 = coef(1,isn,j,i)
      c2 = coef(2,isn,j,i)
      c3 = coef(3,isn,j,i)
      c4 = coef(4,isn,j,i)
      c5 = coef(5,isn,j,i)
      c6 = coef(6,isn,j,i)
      c7 = coef(7,isn,j,i)
      c8 = coef(8,isn,j,i)
      c9 = coef(9,isn,j,i)-bij(j,i)

      if (i == 1) then
        jcol(1:10,ij) = &
          (/index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im), &
            index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im), &
            index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im), &
            index_3to1(2  ,j  ,i)/)
        nzval(1:10,ij) = &
          (/c7,c8,c6, &
            c9,c1,c5, &
            c3,c2,c4, &
            bij(j,i)/)
      elseif (i == nmlon) then
        jcol(1:10,ij) = &
          (/index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im),index_3to1(isn,j-1,i), &
            index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im),index_3to1(isn,j  ,i), &
            index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im),index_3to1(isn,j+1,i), &
                                                          index_3to1(2  ,j  ,i)/)
        nzval(1:10,ij) = &
          (/c8,c6,c7, &
            c1,c5,c9, &
            c2,c4,c3, &
            bij(j,i)/)
      else
        jcol(1:10,ij) = &
          (/index_3to1(isn,j-1,im),index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip), &
            index_3to1(isn,j  ,im),index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip), &
            index_3to1(isn,j+1,im),index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip), &
                                   index_3to1(2  ,j  ,i)/)
        nzval(1:10,ij) = &
          (/c6,c7,c8, &
            c5,c9,c1, &
            c4,c3,c2, &
            bij(j,i)/)
      endif
    enddo

! southern transition latitude, Equation (7.13)
    j = jlatm_JT
    do concurrent (i = 1:nmlon)
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
      ij = index_3to1(isn,j,i)
      rowcnt(ij) = 12

! (C1S+C1N) PhiS(j  ,i+1)
! (C2S+C2N) PhiS(j+1,i+1)
! (C3S+C3N) PhiS(j+1,i  )
! (C4S+C4N) PhiS(j+1,i-1)
! (C5S+C5N) PhiS(j  ,i-1)
! C6S       PhiS(j-1,i-1)
! C7S       PhiS(j-1,i  )
! C8S       PhiS(j-1,i+1)
! (C9S+C9N) PhiS(j  ,i  )
! C6N       PhiN(j-1,i-1)
! C7N       PhiN(j-1,i  )
! C8N       PhiN(j-1,i+1)
      c1 = coef(1,1,j,i)+coef(1,2,j,i)
      c2 = coef(2,1,j,i)+coef(2,2,j,i)
      c3 = coef(3,1,j,i)+coef(3,2,j,i)
      c4 = coef(4,1,j,i)+coef(4,2,j,i)
      c5 = coef(5,1,j,i)+coef(5,2,j,i)
      c6 = coef(6,isn,j,i)
      c7 = coef(7,isn,j,i)
      c8 = coef(8,isn,j,i)
      c9 = coef(9,1,j,i)+coef(9,2,j,i)

      if (i == 1) then
        jcol(1:12,ij) = &
          (/index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im), &
            index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im), &
            index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im), &
            index_3to1(2  ,j-1,i),index_3to1(2  ,j-1,ip),index_3to1(2  ,j-1,im)/)
        nzval(1:12,ij) = &
          (/c7,c8,c6, &
            c9,c1,c5, &
            c3,c2,c4, &
            coef(7,2,j,i),coef(8,2,j,i),coef(6,2,j,i)/)
      elseif (i == nmlon) then
        jcol(1:12,ij) = &
          (/index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im),index_3to1(isn,j-1,i), &
            index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im),index_3to1(isn,j  ,i), &
            index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im),index_3to1(isn,j+1,i), &
            index_3to1(2  ,j-1,ip),index_3to1(2  ,j-1,im),index_3to1(2  ,j-1,i)/)
        nzval(1:12,ij) = &
          (/c8,c6,c7, &
            c1,c5,c9, &
            c2,c4,c3, &
            coef(8,2,j,i),coef(6,2,j,i),coef(7,2,j,i)/)
      else
        jcol(1:12,ij) = &
          (/index_3to1(isn,j-1,im),index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip), &
            index_3to1(isn,j  ,im),index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip), &
            index_3to1(isn,j+1,im),index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip), &
            index_3to1(2  ,j-1,im),index_3to1(2  ,j-1,i),index_3to1(2  ,j-1,ip)/)
        nzval(1:12,ij) = &
          (/c6,c7,c8, &
            c5,c9,c1, &
            c4,c3,c2, &
            coef(6,2,j,i),coef(7,2,j,i),coef(8,2,j,i)/)
      endif
    enddo

! southern low latitudes, Equation (7.14)
    do concurrent (i = 1:nmlon, j = jlatm_JT+1:nmlat_h-1)
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
      ij = index_3to1(isn,j,i)
      rowcnt(ij) = 9

! (C1S+C1N) PhiS(j  ,i+1)
! (C2S+C2N) PhiS(j+1,i+1)
! (C3S+C3N) PhiS(j+1,i  )
! (C4S+C4N) PhiS(j+1,i-1)
! (C5S+C5N) PhiS(j  ,i-1)
! (C6S+C6N) PhiS(j-1,i-1)
! (C7S+C7N) PhiS(j-1,i  )
! (C8S+C8N) PhiS(j-1,i+1)
! (C9S+C9N) PhiS(j  ,i  )
      c1 = coef(1,1,j,i)+coef(1,2,j,i)
      c2 = coef(2,1,j,i)+coef(2,2,j,i)
      c3 = coef(3,1,j,i)+coef(3,2,j,i)
      c4 = coef(4,1,j,i)+coef(4,2,j,i)
      c5 = coef(5,1,j,i)+coef(5,2,j,i)
      c6 = coef(6,1,j,i)+coef(6,2,j,i)
      c7 = coef(7,1,j,i)+coef(7,2,j,i)
      c8 = coef(8,1,j,i)+coef(8,2,j,i)
      c9 = coef(9,1,j,i)+coef(9,2,j,i)

      if (i == 1) then
        jcol(1:9,ij) = &
          (/index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im), &
            index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im), &
            index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im)/)
        nzval(1:9,ij) = &
          (/c7,c8,c6, &
            c9,c1,c5, &
            c3,c2,c4/)
      elseif (i == nmlon) then
        jcol(1:9,ij) = &
          (/index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im),index_3to1(isn,j-1,i), &
            index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im),index_3to1(isn,j  ,i), &
            index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im),index_3to1(isn,j+1,i)/)
        nzval(1:9,ij) = &
          (/c8,c6,c7, &
            c1,c5,c9, &
            c2,c4,c3/)
      else
        jcol(1:9,ij) = &
          (/index_3to1(isn,j-1,im),index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip), &
            index_3to1(isn,j  ,im),index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip), &
            index_3to1(isn,j+1,im),index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip)/)
        nzval(1:9,ij) = &
          (/c6,c7,c8, &
            c5,c9,c1, &
            c4,c3,c2/)
      endif
    enddo

! equator, Equation (7.14)
    j = nmlat_h
    do concurrent (i = 1:nmlon)
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
      ij = index_3to1(isn,j,i)
      rowcnt(ij) = 6

! (C1S+C1N) PhiS(j  ,i+1)
! (C5S+C5N) PhiS(j  ,i-1)
! (C6S+C6N) PhiS(j-1,i-1)
! (C7S+C7N) PhiS(j-1,i  )
! (C8S+C8N) PhiS(j-1,i+1)
! (C9S+C9N) PhiS(j  ,i  )

! there are no PhiS(j+1,i-1), PhiS(j+1,i), PhiS(j+1,i+1)
      c1 = coef(1,1,j,i)+coef(1,2,j,i)
      c5 = coef(5,1,j,i)+coef(5,2,j,i)
      c6 = coef(6,1,j,i)+coef(6,2,j,i)
      c7 = coef(7,1,j,i)+coef(7,2,j,i)
      c8 = coef(8,1,j,i)+coef(8,2,j,i)
      c9 = coef(9,1,j,i)+coef(9,2,j,i)

      if (i == 1) then
        jcol(1:6,ij) = &
          (/index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im), &
            index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im)/)
        nzval(1:6,ij) = &
          (/c7,c8,c6, &
            c9,c1,c5/)
      elseif (i == nmlon) then
        jcol(1:6,ij) = &
          (/index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im),index_3to1(isn,j-1,i), &
            index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im),index_3to1(isn,j  ,i)/)
        nzval(1:6,ij) = &
          (/c8,c6,c7, &
            c1,c5,c9/)
      else
        jcol(1:6,ij) = &
          (/index_3to1(isn,j-1,im),index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip), &
            index_3to1(isn,j  ,im),index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip)/)
        nzval(1:6,ij) = &
          (/c6,c7,c8, &
            c5,c9,c1/)
      endif
    enddo

! continue with northern hemisphere
    isn = 2

! skip northern low latitudes (duplicated from southern low latitudes)

! northern high latitudes, Equation (7.22)
    do concurrent (i = 1:nmlon, j = 3:jlatm_JT-1)
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
      ij = index_3to1(isn,j,i)
      rowcnt(ij) = 10

! C1N     PhiN(j  ,i+1)
! C2N     PhiN(j+1,i+1)
! C3N     PhiN(j+1,i  )
! C4N     PhiN(j+1,i-1)
! C5N     PhiN(j  ,i-1)
! C6N     PhiN(j-1,i-1)
! C7N     PhiN(j-1,i  )
! C8N     PhiN(j-1,i+1)
! (C9N-b) PhiN(j  ,i  )
! b       PhiS(j  ,i  )
      c1 = coef(1,isn,j,i)
      c2 = coef(2,isn,j,i)
      c3 = coef(3,isn,j,i)
      c4 = coef(4,isn,j,i)
      c5 = coef(5,isn,j,i)
      c6 = coef(6,isn,j,i)
      c7 = coef(7,isn,j,i)
      c8 = coef(8,isn,j,i)
      c9 = coef(9,isn,j,i)-bij(j,i)

      if (i == 1) then
        jcol(1:10,ij) = &
          (/index_3to1(1  ,j  ,i), &
            index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im), &
            index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im), &
            index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im)/)
        nzval(1:10,ij) = &
          (/bij(j,i), &
            c3,c2,c4, &
            c9,c1,c5, &
            c7,c8,c6/)
      elseif (i == nmlon) then
        jcol(1:10,ij) = &
          (/                                              index_3to1(1  ,j  ,i), &
            index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im),index_3to1(isn,j+1,i), &
            index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im),index_3to1(isn,j  ,i), &
            index_3to1(isn,j-1,ip),index_3to1(isn,j-1,im),index_3to1(isn,j-1,i)/)
        nzval(1:10,ij) = &
          (/bij(j,i), &
            c2,c4,c3, &
            c1,c5,c9, &
            c8,c6,c7/)
      else
        jcol(1:10,ij) = &
          (/                       index_3to1(1  ,j  ,i), &
            index_3to1(isn,j+1,im),index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip), &
            index_3to1(isn,j  ,im),index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip), &
            index_3to1(isn,j-1,im),index_3to1(isn,j-1,i),index_3to1(isn,j-1,ip)/)
        nzval(1:10,ij) = &
          (/bij(j,i), &
            c4,c3,c2, &
            c5,c9,c1, &
            c6,c7,c8/)
      endif
    enddo

! j=2 in the northern hemisphere is different from others
! north pole (j=1) is set to phi_np, so it is on the right hand side
    j = 2
    do concurrent (i = 1:nmlon)
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
      ij = index_3to1(isn,j,i)
      rowcnt(ij) = 7

! C1N     PhiN(j  ,i+1)
! C2N     PhiN(j+1,i+1)
! C3N     PhiN(j+1,i  )
! C4N     PhiN(j+1,i-1)
! C5N     PhiN(j  ,i-1)
! (C9N-b) PhiN(j  ,i  )
! b       PhiS(j  ,i  )

! PhiN(j-1,i-1), PhiN(j-1,i), PhiN(j-1,i+1) are north pole
! so they are moved to right hand side
      c1 = coef(1,isn,j,i)
      c2 = coef(2,isn,j,i)
      c3 = coef(3,isn,j,i)
      c4 = coef(4,isn,j,i)
      c5 = coef(5,isn,j,i)
      c9 = coef(9,isn,j,i)-bij(j,i)

      if (i == 1) then
        jcol(1:7,ij) = &
          (/index_3to1(1  ,j  ,i), &
            index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im), &
            index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im)/)
        nzval(1:7,ij) = &
          (/bij(j,i), &
            c3,c2,c4, &
            c9,c1,c5/)
      elseif (i == nmlon) then
        jcol(1:7,ij) = &
          (/                                              index_3to1(1  ,j  ,i), &
            index_3to1(isn,j+1,ip),index_3to1(isn,j+1,im),index_3to1(isn,j+1,i), &
            index_3to1(isn,j  ,ip),index_3to1(isn,j  ,im),index_3to1(isn,j  ,i)/)
        nzval(1:7,ij) = &
          (/bij(j,i), &
            c2,c4,c3, &
            c1,c5,c9/)
      else
        jcol(1:7,ij) = &
          (/                       index_3to1(1  ,j  ,i), &
            index_3to1(isn,j+1,im),index_3to1(isn,j+1,i),index_3to1(isn,j+1,ip), &
            index_3to1(isn,j  ,im),index_3to1(isn,j  ,i),index_3to1(isn,j  ,ip)/)
        nzval(1:7,ij) = &
          (/bij(j,i), &
            c4,c3,c2, &
            c5,c9,c1/)
      endif
    enddo

    rowptr(1) = 1
    do i = 2,nlonlat+1
      rowptr(i) = rowptr(i-1)+rowcnt(i-1)
    enddo

    do concurrent (j = 1:rowcnt(1))
      colind(j) = jcol1(j)
      values(j) = nzval1(j)
    enddo

    do concurrent (i = 2:nlonlat)
      do concurrent (j = 1:rowcnt(i))
        colind(rowptr(i)+j-1) = jcol(j,i)
        values(rowptr(i)+j-1) = nzval(j,i)
      enddo
    enddo

  endsubroutine construct_lhs
!-----------------------------------------------------------------------
  pure function construct_rhs(src,coef678) result(rhs)
! construct vector RHS
! this is different from ravel

    use params_module,only:nmlat_h,nmlon
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(2,nmlat_h,nmlon),intent(in) :: src
    real(kind=rp),dimension(nmlon),intent(in) :: coef678
    real(kind=rp),dimension(nlonlat) :: rhs

    integer :: i,j,isn,ij

! start with southern hemisphere
    isn = 1

! south pole has only one point (longitudes are squeezed)
! the volume is a combined cell of all longitudes (a prism with nmlon edges)
    j = 1
    i = 1 ! Equation (7.23) is applied to i=1
    ij = index_3to1(isn,j,i)

! Sum_i S(1,i) + beta PhiN(1,1)
! north pole is set to phi_np, so it on the right hand side
    rhs(ij) = -sum(src(isn,j,:))-beta*phi_np

! other latitudes are the same as ravel

! southern high latitudes, Equation (7.22)
    do concurrent (i = 1:nmlon, j = 2:jlatm_JT-1)
      ij = index_3to1(isn,j,i)
      rhs(ij) = -src(isn,j,i)
    enddo

! southern low latitudes, Equation (7.14)
    do concurrent (i = 1:nmlon, j = jlatm_JT:nmlat_h)
      ij = index_3to1(isn,j,i)
      rhs(ij) = -(src(1,j,i)+src(2,j,i))
    enddo

! continue with northern hemisphere
    isn = 2

! skip northern low latitudes (duplicated from southern low latitudes)

! northern high latitudes, Equation (7.22)
    do concurrent (i = 1:nmlon, j = 3:jlatm_JT-1)
      ij = index_3to1(isn,j,i)
      rhs(ij) = -src(isn,j,i)
    enddo

! j=2 in the northern hemisphere is different from others
! north pole (j=1) is set to phi_np, so it is on the right hand side
! PhiN(j-1,i-1), PhiN(j-1,i), PhiN(j-1,i+1) are north pole
! so they are moved to right hand side
    j = 2
    do concurrent (i = 1:nmlon)
      ij = index_3to1(isn,j,i)
      rhs(ij) = -src(isn,j,i)-coef678(i)*phi_np
    enddo

! skip north pole (not used)

  endfunction construct_rhs
!-----------------------------------------------------------------------
  pure function index_3to1(isn,j,i) result(ij)
! convert lon-lat-hemi triplet (isn,j,i) to linear index ij
! this is the inverse of index_1to3

! isn,j,i need to satisfy 1<=isn<=2, 1<=j<=nmlat_h, 1<=i<=nmlon
! but those conditions are not enforced

    use params_module,only:nmlat_h,nmlon
    use cons_module,only:jlatm_JT

    integer,intent(in) :: isn,j,i
    integer :: ij

    if (isn == 1) then
      if (j == 1) then
        ij = 1
      else
        ij = (j-2)*nmlon+i+1
      endif
    else
      if (j == 1) then
        ij = 0
      elseif (j <= jlatm_JT-1) then
        ij = (nmlat_h+jlatm_JT-2-j)*nmlon+i+1
      else
        ij = (j-2)*nmlon+i+1
      endif
    endif

  endfunction index_3to1
!-----------------------------------------------------------------------
  pure subroutine index_1to3(ij,isn,j,i)
! convert linear index ij to lon-lat-hemi triplet (isn,j,i)
! this is the inverse of index_3to1

! ij needs to be within [1, nlonlat]
! but it is not enforced

    use params_module,only:nmlat_h,nmlon
    use cons_module,only:jlatm_JT

    integer,intent(in) :: ij
    integer,intent(out) :: isn,j,i

    if (ij == 1) then
      isn = 1
      j = 1
      i = 1
    elseif (ij <= (nmlat_h-1)*nmlon+1) then
      isn = 1
      j = (ij-1)/nmlon + 2
      i = modulo(ij-1,nmlon)
      if (i == 0) then
        j = j-1
        i = nmlon
      endif
    else
      isn = 2
      j = nmlat_h+jlatm_JT-2 - (ij-1)/nmlon
      i = modulo(ij-1,nmlon)
      if (i == 0) then
        j = j+1
        i = nmlon
      endif
    endif

  endsubroutine index_1to3
!-----------------------------------------------------------------------
  pure function ravel(fin) result(fout)
! reorder 2D fields (lat-lon) into 1D vector (RHS)
! this is the inverse of unravel (except for the north pole)

    use params_module,only:nmlat_h,nmlon
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(2,nmlat_h,nmlon),intent(in) :: fin
    real(kind=rp),dimension(nlonlat) :: fout

    integer :: i,j,isn,ij

! start with southern hemisphere
    isn = 1

! south pole has only one point (longitudes are squeezed)
    j = 1
    i = 1
    ij = index_3to1(isn,j,i)
    fout(ij) = fin(isn,j,i)

! southern latitudes (no pole)
    do concurrent (i = 1:nmlon, j = 2:nmlat_h)
      ij = index_3to1(isn,j,i)
      fout(ij) = fin(isn,j,i)
    enddo

! continue with northern hemisphere
    isn = 2

! skip northern low latitudes (duplicated from southern low latitudes)

! northern high latitudes (no pole)
    do concurrent (i = 1:nmlon, j = 2:jlatm_JT-1)
      ij = index_3to1(isn,j,i)
      fout(ij) = fin(isn,j,i)
    enddo

! skip north pole (not used)

  endfunction ravel
!-----------------------------------------------------------------------
  pure function unravel(fin) result(fout)
! reorder 1D vector (RHS) into 2D fields (lat-lon)
! this is the inverse of ravel (except for the north pole)

    use params_module,only:nmlat_h,nmlon
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(nlonlat),intent(in) :: fin
    real(kind=rp),dimension(2,nmlat_h,nmlon) :: fout

    integer :: i,j,isn,ij

! start with southern hemisphere
    isn = 1

! south pole has only one point (longitudes are squeezed)
    j = 1
    i = 1
    ij = index_3to1(isn,j,i)
    do concurrent (i = 1:nmlon)
      fout(isn,j,i) = fin(ij)
    enddo

! southern latitudes (no pole)
    do concurrent (i = 1:nmlon, j = 2:nmlat_h)
      ij = index_3to1(isn,j,i)
      fout(isn,j,i) = fin(ij)
    enddo

! continue with northern hemisphere
    isn = 2

! northern low latitudes (duplicated from southern low latitudes)
    do concurrent (i = 1:nmlon, j = jlatm_JT:nmlat_h)
      ij = index_3to1(1,j,i)
      fout(isn,j,i) = fin(ij)
    enddo

! northern high latitudes (no pole)
    do concurrent (i = 1:nmlon, j = 2:jlatm_JT-1)
      ij = index_3to1(isn,j,i)
      fout(isn,j,i) = fin(ij)
    enddo

! north pole is left unset (not in fin)

  endfunction unravel
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
#ifdef USE_MKL
  function solve_mkl(n,nnz,rowptr,colind,values,rhs) result(sol)

    include 'mkl_pardiso.fi'

    integer,intent(in) :: n,nnz
    integer,dimension(n+1),intent(in) :: rowptr
    integer,dimension(nnz),intent(in) :: colind
    real(kind=rp),dimension(nnz),intent(in) :: values
    real(kind=rp),dimension(n),intent(in) :: rhs
    real(kind=rp),dimension(n) :: sol

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

  endfunction solve_mkl
#endif
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
endmodule solver_module
