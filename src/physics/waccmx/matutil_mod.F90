module matutil_mod
! Perform direct calculations of small matrices (rank <= 4).
! Copied from https://fortranwiki.org/fortran/show/Matrix+inversion

  use iso_fortran_env, only: rp=>real64

  implicit none

  contains
!-------------------------------------------------------------------
  pure function matdet2(A) result(d)
! Calculate the determinant of the matrix

    real(kind=rp),dimension(2,2),intent(in) :: A
    real(kind=rp) :: d

    d = A(1,1)*A(2,2) - A(1,2)*A(2,1)

  endfunction matdet2
!-------------------------------------------------------------------
  pure function matdet3(A) result(d)
! Calculate the determinant of the matrix

    real(kind=rp),dimension(3,3),intent(in) :: A
    real(kind=rp) :: d

    d = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
      - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
      + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

  endfunction matdet3
!-------------------------------------------------------------------
  pure function matdet4(A) result(d)
! Calculate the determinant of the matrix

    real(kind=rp),dimension(4,4),intent(in) :: A
    real(kind=rp) :: d

    d = A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))) &
      - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))) &
      + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))) &
      - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

  endfunction matdet4
!-------------------------------------------------------------------
  pure function matadj2(A) result(B)
! Calculate the adjugate of the matrix

    real(kind=rp),dimension(2,2),intent(in) :: A
    real(kind=rp),dimension(2,2) :: B

    B(1,1) =  A(2,2)
    B(2,1) = -A(2,1)
    B(1,2) = -A(1,2)
    B(2,2) =  A(1,1)

  endfunction matadj2
!-------------------------------------------------------------------
  pure function matadj3(A) result(B)
! Calculate the adjugate of the matrix

    real(kind=rp),dimension(3,3),intent(in) :: A
    real(kind=rp),dimension(3,3) :: B

    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))

  endfunction matadj3
!-------------------------------------------------------------------
  pure function matadj4(A) result(B)
! Calculate the adjugate of the matrix

    real(kind=rp),dimension(4,4),intent(in) :: A
    real(kind=rp),dimension(4,4) :: B

    B(1,1) = A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
    B(2,1) = A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
    B(3,1) = A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    B(4,1) = A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    B(1,2) = A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
    B(2,2) = A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
    B(3,2) = A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
    B(4,2) = A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
    B(1,3) = A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
    B(2,3) = A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
    B(3,3) = A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
    B(4,3) = A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
    B(1,4) = A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
    B(2,4) = A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    B(3,4) = A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
    B(4,4) = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

  endfunction matadj4
!-------------------------------------------------------------------
  pure function matinv2(A) result(B)
! Calculate the inverse of the matrix

    real(kind=rp),dimension(2,2),intent(in) :: A
    real(kind=rp),dimension(2,2) :: B

    B = matadj2(A)/matdet2(A)

  endfunction matinv2
!-------------------------------------------------------------------
  pure function matinv3(A) result(B)
! Calculate the inverse of the matrix

    real(kind=rp),dimension(3,3),intent(in) :: A
    real(kind=rp),dimension(3,3) :: B

    B = matadj3(A)/matdet3(A)

  endfunction matinv3
!-------------------------------------------------------------------
  pure function matinv4(A) result(B)
! Calculate the inverse of the matrix

    real(kind=rp),dimension(4,4),intent(in) :: A
    real(kind=rp),dimension(4,4) :: B

    B = matadj4(A)/matdet4(A)

  endfunction matinv4
!-------------------------------------------------------------------
endmodule matutil_mod
