module zmean_phys_fields
  use shr_kind_mod, only: r8=>SHR_KIND_R8
  use ppgrid,       only: begchunk, endchunk, pcols, pver
  use physics_types,only: physics_state
  use cam_history,  only: addfld, outfld, horiz_only
  use infnan,       only: nan, assignment(=)
  use physconst,    only: pi

  use zmean_fields, only: zmean_fields_init, zmean_3d, zmean_2d

  use Zonal_Mean,   only: ZonalAverage_t
  use spmd_utils,   only: masterproc

  implicit none

  integer, parameter :: nzalat = 16
  type(ZonalAverage_t) :: ZA

  real(r8) :: zalats(nzalat)

contains

  subroutine zmean_phys_fields_reg

    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap, cam_grid_register

    type(horiz_coord_t), pointer :: zalon_coord
    type(horiz_coord_t), pointer :: zalat_coord
    integer(iMap),       pointer :: grid_map(:,:)

    real(r8) :: area(nzalat)
    real(r8) :: zalons(1)
    integer :: j

    integer, parameter :: zmean_phys_decomp = 201 ! Must be unique within CAM

    nullify(zalat_coord)
    nullify(zalon_coord)
    nullify(grid_map)

    zalons(1) = 0._r8

    call ZA%init(zalats,area,nzalat,GEN_GAUSSLATS=.true.)

    ! Zonal average grid

    zalat_coord => horiz_coord_create('zalat', '', nzalat, 'latitude',                &
         'degrees_north', 1, nzalat, zalats)
    zalon_coord => horiz_coord_create('zalon', '', 1, 'longitude',                &
         'degrees_east', 1, 1, zalons)

    ! grid decomposition map
    allocate(grid_map(4,nzalat))

    do j = 1,nzalat
       grid_map(1,j) = 1
       grid_map(2,j) = j
       if (masterproc) then
          grid_map(3,j) = 1
          grid_map(4,j) = j
       else
          grid_map(3,j) = 0
          grid_map(4,j) = 0
       end if
    end do

    ! register the zonal average grid
    call cam_grid_register('zavg_phys', zmean_phys_decomp, zalat_coord, &
         zalon_coord, grid_map, unstruct=.false., zonal_grid=.true.)

  end subroutine zmean_phys_fields_reg

  subroutine zmean_phys_fields_init

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

    call addfld ( 'PS_ZA', horiz_only, 'A','Pa', 'PS zonal mean',gridname='zavg_phys')
    call addfld ( 'T_ZA', (/ 'lev' /), 'A', 'K', 'T zonal mean', gridname='zavg_phys')
    call addfld ( 'T_ZA2',(/ 'lev' /), 'A', 'K', 'T zonal mean', gridname='zavg_phys')

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

    do j = 1,nzalat
       call outfld('PS_ZA',PSzavg(j),1,j)
       call outfld('T_ZA',Tzavg(j,:),1,j)
       call outfld('T_ZA2',Tzavg2(j,:),1,j)
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
