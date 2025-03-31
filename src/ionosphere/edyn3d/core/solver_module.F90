module solver_module

  use prec,only:rp
  use params_module,only:nmlat_T1,nmlon

  implicit none

  integer,parameter :: nlonlat = nmlat_T1*nmlon ! solve the whole globe

  contains
!-----------------------------------------------------------------------
  subroutine linear_system(bij,pot_hl,fac_hl,coef_ns2,pot, ierr, errmsg)
! construct linear system based on bij and coef and solve in pot

! if FAC is read in, pot_hl is not used, only fac_hl is used
! if potential is read in, pot_hl is used, fac_hl is output

    use params_module,only:nmlat_h
    use cons_module,only:read_fac
    use mpi_module,only:gather_mag,bcast_3d,mpi_rank, &
      mlon0,mlon1,mlond0,mlond1,mlat0,mlat1,mlatd0,mlatd1

    real(kind=rp),dimension(mlatd0:mlatd1,mlond0:mlond1),intent(in) :: bij
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: pot_hl
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(inout) :: fac_hl
    real(kind=rp),dimension(10,2,mlatd0:mlatd1,mlond0:mlond1),intent(in) :: coef_ns2
    real(kind=rp),dimension(2,mlatd0:mlatd1,mlond0:mlond1),intent(out) :: pot

    integer, intent(out) :: ierr
    character(len=*), intent(out) :: errmsg

#ifdef MKL
    logical,parameter :: use_mkl = .true.
#else
    logical,parameter :: use_mkl = .false.
#endif
    integer,parameter :: root = 0
    integer :: i,j,isn,nnz

    integer,dimension(:),allocatable :: irow,jcol
    real(kind=rp),dimension(:),allocatable :: values
    real(kind=rp),dimension(:,:),allocatable :: bij_full
    real(kind=rp),dimension(:,:,:),allocatable :: pot_hl_full,fac_hl_full
    real(kind=rp),dimension(:,:,:,:),allocatable :: coef_ns2_full
    real(kind=rp),dimension(:,:,:),allocatable :: coef_full
    real(kind=rp),dimension(:,:,:),allocatable :: fac_hl_2,pot_2
    real(kind=rp),dimension(:),allocatable :: rhs,z,pot_hl_f,sol

    character(len=*), parameter :: prefix = 'solver_module::linear_system: '

    errmsg = ' '
    ierr = 0

    allocate(bij_full(nmlat_h,nmlon), stat=ierr)
    if (ierr/=0) then
       errmsg = prefix//'not able to allocate bij_full'
       return
    end if
    allocate(pot_hl_full(2,nmlat_h,nmlon), stat=ierr)
    if (ierr/=0) then
       errmsg = prefix//'not able to allocate pot_hl_full'
       return
    end if
    allocate(fac_hl_full(2,nmlat_h,nmlon), stat=ierr)
    if (ierr/=0) then
       errmsg = prefix//'not able to allocate fac_hl_full'
       return
    end if
    allocate(coef_ns2_full(10,2,nmlat_h,nmlon), stat=ierr)
    if (ierr/=0) then
       errmsg = prefix//'not able to allocate coef_ns2_full'
       return
    end if

    bij_full = gather_mag(bij(mlat0:mlat1,mlon0:mlon1),root)
    if (read_fac) then
      fac_hl_full = gather_mag(fac_hl(:,mlat0:mlat1,mlon0:mlon1),2,root)
    else
      pot_hl_full = gather_mag(pot_hl(:,mlat0:mlat1,mlon0:mlon1),2,root)
    endif
    coef_ns2_full = gather_mag(coef_ns2(:,:,mlat0:mlat1,mlon0:mlon1),10,2,root)

    allocate(fac_hl_2(2,nmlat_h,0:nmlon+1), stat=ierr)
    if (ierr/=0) then
       errmsg = prefix//'not able to allocate fac_hl_2'
       return
    end if
    allocate(pot_2(2,nmlat_h,0:nmlon+1), stat=ierr)
    if (ierr/=0) then
       errmsg = prefix//'not able to allocate pot_2'
       return
    end if

    root_task: if (mpi_rank == root) then

      allocate(coef_full(10,nmlat_T1,nmlon), stat=ierr)
      if (ierr/=0) then
         errmsg = prefix//'not able to allocate coef_full'
         return
      end if

! set the coefficient matrix in both hemispheres
      coef_full = calculate_coef_ns(coef_ns2_full)

! Q from Wu: why coefficients at the south pole are gathered to i=1?
      j = 1
      coef_full(9,j,1) = sum(coef_full(9,j,:))
      coef_full(10,j,1) = sum(coef_full(10,j,:))
      do concurrent (i = 2:nmlon)
        coef_full(9,j,i) = 0
        coef_full(10,j,i) = 0
      enddo

      allocate(irow(10*nlonlat), stat=ierr)
      if (ierr/=0) then
         errmsg = prefix//'not able to allocate irow'
         return
      end if
      allocate(jcol(10*nlonlat), stat=ierr)
      if (ierr/=0) then
         errmsg = prefix//'not able to allocate jcol'
         return
      end if
      allocate(values(10*nlonlat), stat=ierr)
      if (ierr/=0) then
         errmsg = prefix//'not able to allocate values'
         return
      end if

! construct LHS matrix in COO format
      call construct_lhs(bij_full,coef_full(1:9,:,:),nnz,irow,jcol,values)

      allocate(rhs(nlonlat), stat=ierr)
      if (ierr/=0) then
         errmsg = prefix//'not able to allocate rhs'
         return
      end if

! RHS is dense
      rhs = construct_rhs(coef_full(10,:,:))

      deallocate(coef_full)

      allocate(z(nlonlat), stat=ierr)
      if (ierr/=0) then
         errmsg = prefix//'not able to allocate z'
         return
      end if

! determine FAC forcing (dense)
      if (read_fac) then ! input is corrected fac_hl, pot_hl is not used
        z = flatten(fac_hl_full)

      else ! input is pot_hl, fac_hl is to be calculated (output)

        allocate(pot_hl_f(nlonlat), stat=ierr)
        if (ierr/=0) then
           errmsg = prefix//'not able to allocate pot_hl_f'
           return
        end if

! A. Maute 2023/11/21: put the high latitude potential in X
! and then use LHS to calculate the RHS FAC
        pot_hl_f = flatten(pot_hl_full)

! z = matmul(lhs, pot_hl)
        z = 0
        do i = 1,nnz
          z(irow(i)) = z(irow(i))+values(i)*pot_hl_f(jcol(i))
        enddo

        deallocate(pot_hl_f)

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

      deallocate(z)
      allocate(sol(nlonlat), stat=ierr)
      if (ierr/=0) then
         errmsg = prefix//'not able to allocate sol'
         return
      end if

      if (use_mkl) then
        sol = solve_mkl(nlonlat,nnz,irow(1:nnz),jcol(1:nnz),values(1:nnz),rhs)
      else
        sol = solve_superlu(nlonlat,nnz,irow(1:nnz),jcol(1:nnz),values(1:nnz),rhs)
      endif

      deallocate(irow)
      deallocate(jcol)
      deallocate(values)
      deallocate(rhs)

! reconstruct 2D distribution of potential based on the solution
      pot_2(:,:,1:nmlon) = unravel(sol)

      deallocate(sol)

! periodic points
      do j = 1,nmlat_h
        do isn = 1,2
          pot_2(isn,j,0) = pot_2(isn,j,nmlon)
          pot_2(isn,j,nmlon+1) = pot_2(isn,j,1)
        enddo
      enddo

    endif root_task

    deallocate(bij_full)
    deallocate(pot_hl_full)
    deallocate(fac_hl_full)
    deallocate(coef_ns2_full)

    if (.not. read_fac) then
      call bcast_3d(fac_hl_2,2,nmlat_h,nmlon+2,root)

      do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
        fac_hl(isn,j,i) = fac_hl_2(isn,j,i)
      enddo
    endif

    deallocate(fac_hl_2)

    call bcast_3d(pot_2,2,nmlat_h,nmlon+2,root)
    do concurrent (i = mlond0:mlond1, j = mlatd0:mlatd1, isn = 1:2, j>=1 .and. j<=nmlat_h)
      pot(isn,j,i) = pot_2(isn,j,i)
    enddo

    deallocate(pot_2)

  endsubroutine linear_system
!-----------------------------------------------------------------------
  pure function calculate_coef_ns(coef_ns2) result(coef_ns)
! set the coefficient matrix in both hemispheres
! decide where to add the SH and NH stencil and where not
! change the direction of NH stencil from coef_ns2 to coef_ns
! LHS+RHS for each P-point
! A. Maute 2023/02: solve two hemispheres

    use params_module,only:nmlat_h
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(10,2,nmlat_h,nmlon),intent(in) :: coef_ns2
    real(kind=rp),dimension(10,nmlat_T1,nmlon) :: coef_ns

    integer :: i,j,ic,jS,jN ! overall index from pole to equator

    do concurrent (i = 1:nmlon, j = 1:nmlat_h)
      jS = j
      jN = nmlat_T1-j+1

! j should take one of the following four situations

! from pole to latm_JT, set coefficients separately in two hemispheres
      if (j <= jlatm_JT-1) then
        do concurrent (ic = 1:10)
          coef_ns(ic,jS,i) = coef_ns2(ic,1,j,i)
          coef_ns(ic,jN,i) = coef_ns2(ic,2,j,i)
        enddo
      endif

      if (j == jlatm_JT) then

! add values from both hemispheres
        do concurrent (ic = 1:5)
          coef_ns(ic,jS,i) = coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i)
          coef_ns(ic,jN,i) = coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i)
        enddo
        do concurrent (ic = 9:10)
          coef_ns(ic,jS,i) = coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i)
          coef_ns(ic,jN,i) = coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i)
        enddo

! don't add values from the other hemisphere
        do concurrent (ic = 6:8)
          coef_ns(ic,jS,i) = coef_ns2(ic,1,j,i)
          coef_ns(ic,jN,i) = coef_ns2(ic,2,j,i)
        enddo
      endif

! from latm_JT to equator, add values from both hemispheres
      if (j>=jlatm_JT+1 .and. j<=nmlat_h-1) then
        do concurrent (ic = 1:10)
          coef_ns(ic,jS,i) = coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i)
          coef_ns(ic,jN,i) = coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i)
        enddo
      endif

! set equatorial boundary condition (page 14 Art's notes)
! there should be just one equator value
      if (j == nmlat_h) then
        coef_ns(1,j,i) = coef_ns2(1,1,j,i)+coef_ns2(1,2,j,i)
        coef_ns(5,j,i) = coef_ns2(5,1,j,i)+coef_ns2(5,2,j,i)
        coef_ns(9,j,i) = coef_ns2(9,1,j,i)+coef_ns2(9,2,j,i)
        coef_ns(10,j,i) = coef_ns2(10,1,j,i)+coef_ns2(10,2,j,i)
        do concurrent (ic = 6:8)
          coef_ns(ic,j,i) = (coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i))/2
        enddo

! Q from Wu: Why isn't it (coef_ns2(ic,1,j,i)+coef_ns2(ic,2,j,i))/2?
        do ic = 2,4
          coef_ns(ic,j,i) = coef_ns(10-ic,j,i) ! 2,3,4 <- 8,7,6
        enddo
      endif
    enddo

  endfunction calculate_coef_ns
!-----------------------------------------------------------------------
  subroutine construct_lhs(bij,coef,nnz,irow,jcol,values)
! construct LHS matrix (sparse, coordinate form)

! need to set where the two hemispheres are connected
! This is not in the 9-point stencil but needs to be done manually

    use params_module,only:nmlat_h
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(nmlat_h,nmlon),intent(in) :: bij
    real(kind=rp),dimension(9,nmlat_T1,nmlon),intent(in) :: coef
    integer,intent(out) :: nnz
    integer,dimension(10*nlonlat),intent(out) :: irow,jcol
    real(kind=rp),dimension(10*nlonlat),intent(out) :: values

! if two hemispheres are uncoupled at high latitudes, set bijSum to zero
    real(kind=rp),parameter :: bijSum = 0
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
  pure function construct_rhs(coef_10) result(rhs)
! construct vector RHS
! this is different from flatten

    use params_module,only:nmlat_h
    use cons_module,only:phi_pol

    real(kind=rp),dimension(nmlat_T1,nmlon),intent(in) :: coef_10
    real(kind=rp),dimension(nlonlat) :: rhs

    integer :: i,j,ij

    rhs = 0

! set up poles, there are no c6,c7,c8 values
    j = 1

! for equation i=1 at the south pole
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

  endfunction construct_rhs
!-----------------------------------------------------------------------
  function solve_mkl(n,nnz,irow,jcol,values,rhs) result(sol)

#ifdef MKL
    include 'mkl_pardiso.fi'
#endif

    integer,intent(in) :: n,nnz
    integer,dimension(nnz),intent(in) :: irow,jcol
    real(kind=rp),dimension(nnz),intent(in) :: values
    real(kind=rp),dimension(n),intent(in) :: rhs
    real(kind=rp),dimension(n) :: sol

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
    real(kind=rp),dimension(n) :: rhs_cp
    real(kind=rp),dimension(nnz) :: nzval
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
    do concurrent (i = 1:n)
      rhs_cp(i) = rhs(i)
    enddo

! use 32-bit integer version
! if the number of non-zero elements is on the order of 500 million or more
! then use pardiso_64 (64-bit integer version)
    call pardiso(pt, maxfct, mnum, mtype, phase, n, &
      nzval, rowptr, colind, perm, nrhs, iparm, msglvl, rhs_cp, sol, error)

    write(6,"('phase ',i4,' error ',i4)") phase,error
#else
    sol = 0 ! to suppress the "return value not set" warning
#endif

  endfunction solve_mkl
!-----------------------------------------------------------------------
  function solve_superlu(n,nnz,irow,jcol,values,rhs) result(sol)

    integer,intent(in) :: n,nnz
    integer,dimension(nnz),intent(in) :: irow,jcol
    real(kind=rp),dimension(nnz),intent(in) :: values
    real(kind=rp),dimension(n),intent(in) :: rhs
    real(kind=rp),dimension(n) :: sol

! for SuperLU sparse matrix solver (CSC format)
    integer,parameter :: nrhs = 1
    integer :: i,iopt,info
    integer(kind=8) :: f_factors
    integer,dimension(n+1) :: colptr
    integer,dimension(nnz) :: rowind
    real(kind=rp),dimension(nnz) :: nzval

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

! SuperLU needs CSC format
    call coo_to_csc(n,n,nnz,irow,jcol,values,colptr,rowind,nzval)

    do concurrent (i = 1:n)
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

  endfunction solve_superlu
!-----------------------------------------------------------------------
  subroutine coo_to_csr(nrow,ncol,nnz,irow,jcol,values,rowptr,colind,nzval)
! convert sparse matrix format from COO to CSR

! Compressed Sparse Row (CSR) format:
! rowptr is the cumulative sum of the number of non-zero elements in each row
! i.e., rowptr(i+1)-rowptr(i) is the number of non-zero elements in the ith row
! colind(rowptr(i):rowptr(i+1)-1) are the column indices of non-zero elements in the ith column

    integer,intent(in) :: nrow,ncol,nnz
    integer,dimension(nnz),intent(in) :: irow,jcol
    real(kind=rp),dimension(nnz),intent(in) :: values
    integer,dimension(nrow+1),intent(out) :: rowptr
    integer,dimension(nnz),intent(out) :: colind
    real(kind=rp),dimension(nnz),intent(out) :: nzval

    integer :: i,j,last_colind,nnz_i
    integer,dimension(ncol) :: jcol_i,idx
    real(kind=rp),dimension(ncol) :: values_i

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
    real(kind=rp),dimension(nnz),intent(in) :: values
    integer,dimension(ncol+1),intent(out) :: colptr
    integer,dimension(nnz),intent(out) :: rowind
    real(kind=rp),dimension(nnz),intent(out) :: nzval

    integer :: i,j,last_rowind,nnz_j
    integer,dimension(nrow) :: irow_j,idx
    real(kind=rp),dimension(nrow) :: values_j

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
  pure function flatten(fin) result(fout)
! reorder 2D fields (lat-lon) into 1D vector (RHS)
! northern/southern hemispheres are either separate or averaged
! based on their latitude ranges (high-lat, transition, low-lat, equator)

    use params_module,only:nmlat_h
    use cons_module,only:jlatm_JT

    real(kind=rp),dimension(2,nmlat_h,nmlon),intent(in) :: fin
    real(kind=rp),dimension(nlonlat) :: fout

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

  endfunction flatten
!-----------------------------------------------------------------------
  pure function unravel(fin) result(fout)
! reorder 1D vector (RHS) into 2D fields (lat-lon)

    use params_module,only:nmlat_h

    real(kind=rp),dimension(nlonlat),intent(in) :: fin
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
  pure function argsort(a,n,amin,amax) result(idx)
! a simplified radix sort for row or column indices

    integer,intent(in) :: n,amin,amax
    integer,dimension(n),intent(in) :: a
    integer,dimension(n) :: idx

    integer :: i,cnt
    integer,dimension(amin:amax) :: b

! put all elements in place
    b = amin-1
    do concurrent (i = 1:n)
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
endmodule solver_module
