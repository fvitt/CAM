!-----------------------------------------------------------------------
!  Heelis high-latitude empirical electric potential for 3D dynamo
!-----------------------------------------------------------------------
module edyn3D_heelis

  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi

  implicit none

  private
  public :: edyn3D_heelis_set_hlat_pot

contains

  ! set poten_hl for edyn solver input
  subroutine edyn3D_heelis_set_hlat_pot(sunlon)

    use edyn3D_fieldline, only: poten_hl ! (0:nmlon+1,nmlat_T1)
    use edyn3D_params, only: nmlon, nmlat_h, nmlat_T1, ylatm, ylonm
    use heelis_mod,only: heelis_update, heelis_flwv32

    real(r8),intent(in) :: sunlon ! mag longitude sun location

    real(r8) :: xlat(nmlon)
    real(r8) :: xlon(nmlon)
    real(r8) :: ratio(nmlon)
    integer :: iflag(nmlon)
    integer :: isn, j, jj

    call heelis_update()

    ratio(:) = 1._r8
    poten_hl(:,:) = 0._r8

    do isn = 1,2
       do j = 1,nmlat_h
          if (abs(ylatm(j,isn)) > pi/6._r8) then
             if (isn==1) then
                jj = j
             else
                jj = nmlat_T1 - j + 1
             end if

             xlat(:) = ylatm(j,isn)
             xlon(:) = ylonm(1:nmlon)-sunlon
             iflag(:) = 1 ! must be updated at each j

             call heelis_flwv32(xlat,xlon,ratio,pi,iflag,nmlon,poten_hl(1:nmlon+1,jj))
          end if
       end do
    end do

    ! wrap around longitude point
    poten_hl(0,:) = poten_hl(nmlon,:)

  end subroutine edyn3D_heelis_set_hlat_pot

end module edyn3D_heelis
