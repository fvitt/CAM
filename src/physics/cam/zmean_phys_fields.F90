module zmean_phys_fields
  use shr_kind_mod, only: r8=>SHR_KIND_R8
  use ppgrid,       only: begchunk, endchunk, pcols, pver
  use physics_types,only: physics_state
  use cam_history,  only: addfld, outfld, horiz_only
  use infnan,       only: nan, assignment(=)
  use physconst,    only: pi

  use zmean_fields, only: zmean_fields_init, zmean_3d, zmean_2d

  use Zonal_Mean,   only: ZonalAverage_t

  implicit none

  integer, parameter :: nzalat = 16
  type(ZonalAverage_t) :: ZA

  real(r8) :: zalats(nzalat)


contains

  subroutine zmean_phys_fields_init


    real(r8) :: area(nzalat)

    call zmean_fields_init

    call addfld ( 'Tfld', (/ 'lev' /), 'A', 'K', 'T field')
    call addfld ( 'Tzmn', (/ 'lev' /), 'A', 'K', 'T zonal mean')
    call addfld ( 'Tzmn2', (/ 'lev' /), 'A', 'K', 'T zonal mean')
    call addfld ( 'Tzavg', (/ 'lev' /), 'A', 'K', 'T zonal mean')
    call addfld ( 'Tzavg2', (/ 'lev' /), 'A', 'K', 'T zonal mean')
    call addfld ( 'Ufld', (/ 'lev' /), 'A', 'K', 'T field')
    call addfld ( 'Uzmn', (/ 'lev' /), 'A', 'K', 'T zonal mean')
    call addfld ( 'Vfld', (/ 'lev' /), 'A', 'K', 'T field')
    call addfld ( 'Vzmn', (/ 'lev' /), 'A', 'K', 'T zonal mean')

    call addfld ( 'PSfld', horiz_only, 'A', 'Pa', 'T zonal mean')
    call addfld ( 'PSzmn', horiz_only, 'A', 'Pa', 'T zonal mean')
    call addfld ( 'PSzavg', horiz_only, 'A','Pa', 'T zonal mean')

    call ZA%init(zalats,area,nzalat,GEN_GAUSSLATS=.true.)

  end subroutine zmean_phys_fields_init

  subroutine zmean_phys_fields_timestep_init(phys_state)
    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    real(r8) :: Tfld(pcols,pver,begchunk:endchunk)
    real(r8) :: Tzmfld(pcols,pver,begchunk:endchunk)
    real(r8) :: Tzafld(pcols,pver,begchunk:endchunk)
    real(r8) :: Tzafld2(pcols,pver,begchunk:endchunk)
    real(r8) :: Tzmfld2(pcols,pver,begchunk:endchunk)
    real(r8) :: Ufld(pcols,pver,begchunk:endchunk)
    real(r8) :: Uzmfld(pcols,pver,begchunk:endchunk)
    real(r8) :: Vfld(pcols,pver,begchunk:endchunk)
    real(r8) :: Vzmfld(pcols,pver,begchunk:endchunk)

    real(r8) :: PSfld(pcols,begchunk:endchunk)
    real(r8) :: PSzmfld(pcols,begchunk:endchunk)

    real(r8) :: fld_tmp(pcols,pver)
    integer :: lchnk,ncol, i, k

    real(r8) :: Tzavg(nzalat,pver)
    real(r8) :: Tzavg2(nzalat,pver)
    real(r8) :: PSzavg(nzalat)
    real(r8) :: PSzafld(pcols,begchunk:endchunk)
    real(r8) :: dlat, clat
    integer :: j, jzm

    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol
       do i = 1,ncol
          Tfld(i,:,lchnk) = phys_state(lchnk)%t(i,:)
          Ufld(i,:,lchnk) = phys_state(lchnk)%u(i,:)
          Vfld(i,:,lchnk) = phys_state(lchnk)%v(i,:)
          PSfld(i,lchnk) = phys_state(lchnk)%ps(i)
       end do
    end do

    Tzmfld = zmean_3d( Tfld, pver )
    Uzmfld = zmean_3d( Ufld, pver )
    Vzmfld = zmean_3d( Vfld, pver )
    PSzmfld = zmean_2d( PSfld )

    do k = 1,pver
       Tzmfld2(:,k,:) = zmean_2d(Tfld(:,k,:))
    end do


    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol

       fld_tmp(:ncol,:) = Tfld(:ncol,:,lchnk)
       call outfld( 'Tfld', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = Tzmfld(:ncol,:,lchnk)
       call outfld( 'Tzmn', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = Tzmfld2(:ncol,:,lchnk)
       call outfld( 'Tzmn2', fld_tmp(:ncol,:), ncol, lchnk)

       fld_tmp(:ncol,:) = Ufld(:ncol,:,lchnk)
       call outfld( 'Ufld', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = Uzmfld(:ncol,:,lchnk)
       call outfld( 'Uzmn', fld_tmp(:ncol,:), ncol, lchnk)

       fld_tmp(:ncol,:) = Vfld(:ncol,:,lchnk)
       call outfld( 'Vfld', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = Vzmfld(:ncol,:,lchnk)
       call outfld( 'Vzmn', fld_tmp(:ncol,:), ncol, lchnk)


       fld_tmp(:ncol,1) = PSfld(:ncol,lchnk)
       call outfld( 'PSfld', fld_tmp(:ncol,1), ncol, lchnk)
       fld_tmp(:ncol,1) = PSzmfld(:ncol,lchnk)
       call outfld( 'PSzmn', fld_tmp(:ncol,1), ncol, lchnk)

    end do

    call ZA%binAvg(PSfld,PSzavg)
    call ZA%binAvg(Tfld,Tzavg)

    do k=1,pver
       call ZA%binAvg(Tfld(:,k,:),Tzavg2(:,k))
    end do

    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol
       do i = 1,ncol
          clat = phys_state(lchnk)%lat(i)*180._r8/pi !degrees

          dlat = huge(1._r8)
          jzm = -huge(1)
          do j = 1,nzalat
             if ( abs( clat - zalats(j) ) < dlat ) then
                dlat = abs( clat - zalats(j) )
                jzm = j
             end if
          end do
          PSzafld(i,lchnk) = PSzavg(jzm)
          do k=1,pver
             Tzafld(i,k,lchnk) = Tzavg(jzm,k)
             Tzafld2(i,k,lchnk) = Tzavg2(jzm,k)
          end do
       end do

       fld_tmp(:ncol,1) = PSzafld(:ncol,lchnk)
       call outfld('PSzavg', fld_tmp(:ncol,1), ncol, lchnk)

       fld_tmp(:ncol,:) = Tzafld(:ncol,:,lchnk)
       call outfld('Tzavg', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = Tzafld2(:ncol,:,lchnk)
       call outfld('Tzavg2', fld_tmp(:ncol,:), ncol, lchnk)

    end do

  end subroutine zmean_phys_fields_timestep_init


end module zmean_phys_fields
