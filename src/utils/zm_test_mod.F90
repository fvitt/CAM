module zm_test_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: begchunk, endchunk, pcols, pver
  use physics_types, only: physics_state
  use phys_grid, only: get_ncols_p
  use spmd_utils, only: masterproc
  use cam_history, only: horiz_only, addfld, outfld

  use esmf_zm_mod, only: lon_beg, lon_end, lat_beg, lat_end, esmf_zm_calc_2d, esmf_zm_calc_3d
  use esmf_zm_mod, only: nlats, glats

  implicit none

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine zm_test_reg()

    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap, cam_grid_register

    integer :: ind, j, astat

    integer, parameter :: esmf_zonal_mean_decomp = 334 ! Must be unique within CAM
    type(horiz_coord_t), pointer :: zmlon_coord
    type(horiz_coord_t), pointer :: zmlat_coord
    integer(iMap),       pointer :: grid_map(:,:)
    real(r8) :: zmlons(1)

    ! Zonal mean grid for history fields
    zmlons = 0._r8

    zmlat_coord => horiz_coord_create('zmlat', '', nlats, 'latitude', 'degrees_north', 1, nlats, glats)
    zmlon_coord => horiz_coord_create('zmlon', '', 1,    'longitude', 'degrees_east',  1, 1,   zmlons)

    ! grid decomposition map
    allocate(grid_map(4,lat_end-lat_beg+1), stat=astat)

    ind = 0
    do j = lat_beg, lat_end
       ind = ind + 1
       grid_map(1,ind) = 1
       grid_map(2,ind) = j
       if (lon_beg==1) then
          grid_map(3,ind) = 1
          grid_map(4,ind) = j
       else
          grid_map(3,ind) = 0
          grid_map(4,ind) = 0
       end if
    end do

    ! register the zonal average grid
    call cam_grid_register('esmf_zonal_mean', esmf_zonal_mean_decomp, zmlat_coord, zmlon_coord, grid_map, &
                           unstruct=.false., zonal_grid=.true.)


  end subroutine zm_test_reg

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine zm_test_init()

    call addfld ('T_emsfzm',  (/'lev'/), 'A', 'K',  'Zonal-Mean temperature', gridname='esmf_zonal_mean' )
    call addfld ('PS_emsfzm', horiz_only, 'A','Pa', 'Zonal-Mean srf press',   gridname='esmf_zonal_mean' )

  end subroutine zm_test_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine zm_test_run(phys_state)

    use time_manager, only: get_nstep

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    real(r8) :: physfld(pver,pcols,begchunk:endchunk)
    real(r8) :: psfld(pcols,begchunk:endchunk)

    integer :: lchnk, ncol, icol, nstep, j

    real(r8) :: temp_zm(lat_beg:lat_end,pver)
    real(r8) :: ps_zm(lat_beg:lat_end)

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       do icol = 1,ncol
          physfld(:pver,icol,lchnk) = phys_state(lchnk)%t(icol,:pver)
          psfld(icol,lchnk) = phys_state(lchnk)%ps(icol)
       end do
    end do

    temp_zm = esmf_zm_calc_3d(physfld)

    nstep = get_nstep()

    ps_zm = esmf_zm_calc_2d(psfld)

    do j = lat_beg, lat_end
       call outfld('PS_emsfzm',ps_zm(j),1,j)
       call outfld('T_emsfzm',temp_zm(j,:),1,j)
    end do

  end subroutine zm_test_run


end module zm_test_mod
