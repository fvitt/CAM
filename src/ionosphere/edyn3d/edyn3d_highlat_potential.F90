module edyn3d_highlat_potential
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use mpi_module, only: mlond0, mlond1, mlatd0, mlatd1, mlon0, mlon1

  implicit none

  real(r8), protected, allocatable :: hilat_potential(:,:,:)

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_highlat_potential_alloc()

    integer :: astat

    character, parameter :: prefix = 'edyn3d_highlat_potential_alloc: '

    allocate(hilat_potential(2,mlatd0:mlatd1,mlond0:mlond1), stat=astat)
    if (astat /= 0) then
       call endrun(prefix//'failed to allocate hilat_potential')
    end if

    hilat_potential = 0._r8

  end subroutine edyn3d_highlat_potential_alloc

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_heelis_update()
    use params_module, only: nmlat_h, nmlat_T1, nmlon
    use params_module, only: ylonm, ylatm
    use heelis_mod,only: heelis_update, heelis_flwv32
    use sunloc_mod, only: sunloc_calc
    use physconst, only: pi

    real(r8) :: sunlon ! mag longitude sun location
    real(r8) :: xlat(nmlon)
    real(r8) :: xlon(nmlon)
    real(r8) :: pot(nmlon+1)
    real(r8) :: ratio(nmlon)
    integer :: iflag(nmlon)
    integer :: isn, j, jj

    call heelis_update()

    ratio(:) = 1._r8
    hilat_potential = -huge(1._r8)

    call sunloc_calc(sunlon)

    do isn = 1,2
       do j = max(1,mlatd0),min(mlatd1,nmlat_h)

          if (isn==1) then
             jj = j
          else
             jj = nmlat_T1 - j + 1
          end if

          xlat(:) = ylatm(isn,j)
          xlon(:) = ylonm(1:nmlon)-sunlon
          iflag(:) = 1 ! must be updated at each j

          call heelis_flwv32(xlat,xlon,ratio,pi,iflag,nmlon,pot)
          hilat_potential(isn,j,max(1,mlond0):mlond1) = pot(max(1,mlond0):mlond1)

          ! wrap around longitude points
          if (mlon0 == 1) then
             hilat_potential(isn,j,0) = pot(nmlon)
          else if(mlon1 == nmlon) then
             hilat_potential(isn,j,nmlon+1) = pot(1)
          end if

       end do
    end do

  end subroutine edyn3d_heelis_update

end module edyn3d_highlat_potential
