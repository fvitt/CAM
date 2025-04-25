module edyn3d_highlat_potential
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use sunloc_mod, only: sunloc_calc, sunloc_calc2
  use physconst, only: pi
  use mpi_module, only: mlond0, mlond1, mlatd0, mlatd1, mlon0, mlon1
  use params_module, only: nmlat_h, nmlat_T1, nmlon
  use params_module, only: ylonm, ylatm
  use infnan, only: nan, assignment(=)

  implicit none

  private

  public :: edyn3d_highlat_potential_init
  public :: edyn3d_highlat_potential_update
  public :: edyn3d_highlat_potential_get

  character(len=6) :: hilat_pot_model = 'NONE'

  real(r8), allocatable :: phihm(:,:)
  real(r8), allocatable :: maglon(:)
  real(r8), allocatable :: maglat(:)

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_highlat_potential_init(method,wei05_ncfile)
    use weimer_highlat_potential, only: weimer_highlat_potential_init

    character(len=*),intent(in) :: method
    character(len=*),intent(in) :: wei05_ncfile

    character(len=*), parameter :: prefix = 'edyn3d_highlat_potential_init: '

    integer :: astat

    hilat_pot_model = method

    allocate(phihm(nmlon,nmlat_T1), stat=astat)
    if (astat /= 0) then
       call endrun(prefix//'not able to allocate phihm')
    end if
    phihm = nan

    allocate(maglon(nmlon), stat=astat)
    if (astat /= 0) then
       call endrun(prefix//'not able to allocate maglon')
    end if
    allocate(maglat(nmlat_T1), stat=astat)
    if (astat /= 0) then
       call endrun(prefix//'not able to allocate maglat')
    end if

    maglon(1:nmlon) = ylonm(1:nmlon)
    maglat(1:nmlat_h) = ylatm(1,1:nmlat_h)
    maglat(nmlat_h:)  = ylatm(2,nmlat_h:1:-1)

    if (hilat_pot_model == 'heelis') then

    else if (hilat_pot_model == 'weimer') then
       call weimer_highlat_potential_init(wei05_ncfile)
    else
       call endrun(prefix//' hilat_pot_model '//trim(hilat_pot_model)//' not recognized')
    end if

  end subroutine edyn3d_highlat_potential_init


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_highlat_potential_update

    character(len=*), parameter :: prefix = 'edyn3d_highlat_potential_update: '

    if (hilat_pot_model == 'heelis') then
       call edyn3d_heelis_update()
    else if (hilat_pot_model == 'weimer') then
       call edyn3d_weimer_update()
    else
       call endrun(prefix//' hilat_pot_model '//trim(hilat_pot_model)//' not recognized')
    end if

  end subroutine edyn3d_highlat_potential_update


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_highlat_potential_get(hl_pot)
    real(r8), intent(out) :: hl_pot(2,mlatd0:mlatd1,mlond0:mlond1)

    integer :: h,i,j,jj

    do h = 1,2
       do i = max(1,mlond0),min(mlond1,nmlon)
          do j = max(1,mlatd0),min(mlatd1,nmlat_h)
             if (h==1) then
                jj = j
             else
                jj = nmlat_T1 - j + 1
             end if
             hl_pot(h,j,i) = phihm(i,jj)

             ! wrap around longitude points
             if (i==1) then
                hl_pot(h,j,0) = phihm(nmlon,jj)
             else if (i==nmlon) then
                hl_pot(h,j,nmlon+1) = phihm(1,jj)
             end if

          end do

       end do
    end do

  end subroutine edyn3d_highlat_potential_get

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_heelis_update()
    use heelis_mod,only: heelis_update, heelis_flwv32

    real(r8) :: sunlon ! mag longitude sun location
    real(r8) :: xlat(nmlon)
    real(r8) :: xlon(nmlon)
    real(r8) :: pot(nmlon+1)
    real(r8) :: ratio(nmlon)
    integer :: iflag(nmlon)
    integer :: j

    call heelis_update()

    ratio(:) = 1._r8

    call sunloc_calc2(sunlon)

    do j = 1,nmlat_T1

       xlat(:) = maglat(j)
       xlon(:) = maglon(:)-sunlon
       iflag(:) = 1 ! must be updated at each j

       call heelis_flwv32(xlat,xlon,ratio,pi,iflag,nmlon,pot)
       phihm(1:nmlon,j) = pot(1:nmlon)

    end do

  end subroutine edyn3d_heelis_update

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_weimer_update()

    use weimer_highlat_potential, only: weimer_highlat_potential_update
    use solar_wind_data,  only: solar_wind_advance
    use solar_wind_data,  only: bzimf=>solar_wind_bzimf
    use solar_wind_data,  only: byimf=>solar_wind_byimf
    use solar_wind_data,  only: swvel=>solar_wind_swvel
    use solar_wind_data,  only: swden=>solar_wind_swden

    real(r8) :: sunlon ! mag longitude sun location

    ! update solar wind data (IMF, etc.)
    call solar_wind_advance()

    call sunloc_calc2(sunlon)

    call weimer_highlat_potential_update(byimf, bzimf, swvel, swden, sunlon, nmlon, nmlat_T1, &
                                         maglon, maglat, phihm)

  end subroutine edyn3d_weimer_update

end module edyn3d_highlat_potential
