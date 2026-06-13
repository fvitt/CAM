module lbc_mod

  use iso_fortran_env, only: rp=>real64

  implicit none

! lower boundary condition
  real(kind=rp),dimension(3) :: fb = -huge(1._rp)
  real(kind=rp),dimension(3,3) :: b = -huge(1._rp)

  contains
!-----------------------------------------------------------------------
  subroutine init(dz)

    use matutil_mod,only:matinv3

    real(kind=rp),intent(in) :: dz

    real(kind=rp),parameter :: &
      alfa = 0.234_rp, &    ! lower boundary for O2+O (0.22+0.14)
      pshelb = 0.1154e-5_rp ! lower boundary for Helium (mmr)
    real(kind=rp),dimension(3),parameter :: &
      g = -(/alfa,0.0_rp,pshelb/) ! g = -(O2+O 0 He)
    real(kind=rp),dimension(3,3),parameter :: &
!     |0 0 0|
! e = |0 1 0|
!     |0 0 0|
      e = reshape((/0,0,0,0,1,0,0,0,0/),(/3,3/)), &
!     |1  1  0|
! f = |0 -1  0|
!     |0  0  1|
      f = reshape((/1,0,0,1,-1,0,0,0,1/),(/3,3/))
    integer :: n
    real(kind=rp),dimension(3,3) :: wm1,wm2,wm3

! calculate matrix b(3,3) and vector fb(3)
! representing the lower boundary condition in major,
! where psi = (O2 O He) are calculated as
! psi(k=-1/2) = b * psi(k=1/2) + fb

! first define 3x3 matrices e, f and length-3 vector g
! in the general lower boundary condition
! e * d(psi)/ds + f * psi + g = 0

! then evaluates b and fb from:
! b = (e/ds - f/2)**(-1) * (e/ds + f/2)
! fb = (e/ds - f/2)**(-1) * g

! wm1 = (e/ds - f/2)
! wm2 = (e/ds + f/2)
    wm1 = e/dz - f/2
    wm2 = e/dz + f/2

! now invert wm1 in wm3
    wm3 = matinv3(wm1)

! b = wm3 * wm2
    b = matmul(wm3,wm2)

! fb = wm3 * g
    do n = 1,3
      fb(n) = dot_product(wm3(n,:),g)
    enddo

  endsubroutine init
!-----------------------------------------------------------------------
  elemental subroutine major_lbc(o2,o1,he,bo2,bo1,bhe)
! O2, O, He at midpoint level 0 (one level below the midpoint level 1)
! this is half level below the lower boundary of the model

    real(kind=rp),intent(in) :: o2,o1,he
    real(kind=rp),intent(out) :: bo2,bo1,bhe

    bo2 = fb(1)+b(1,1)*o2+b(1,2)*o1+b(1,3)*he
    bo1 = fb(2)+b(2,1)*o2+b(2,2)*o1+b(2,3)*he
    bhe = fb(3)+b(3,1)*o2+b(3,2)*o1+b(3,3)*he

  endsubroutine major_lbc
!-----------------------------------------------------------------------
  elemental subroutine set_lbc(tlbc,ulbc,vlbc,zlbc)

    real(kind=rp),intent(out) :: tlbc,ulbc,vlbc,zlbc

    real(kind=rp),parameter :: &
      tbound = 181, &                    ! background T at lower boundary (K)
      zbound = 136.291e5_rp/sqrt(2.0_rp) ! background Z at lower boundary (cm)

    tlbc = tbound
    ulbc = 0
    vlbc = 0
    zlbc = zbound

  endsubroutine set_lbc
!-----------------------------------------------------------------------
endmodule lbc_mod
