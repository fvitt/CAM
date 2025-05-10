module zm_test_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: begchunk, endchunk, pcols, pver
  use physics_types, only: physics_state
  use phys_grid, only: get_ncols_p
  use spmd_utils, only: masterproc
  use cam_history, only: horiz_only, addfld, outfld


  use esmf_lonlat_grid_mod, only: lon_beg, lon_end, lat_beg, lat_end, nlat, glats
  use esmf_phys2lonlat_mod, only: esmf_phys2lonlat_regrid
  use esmf_zonal_mean_mod, only: esmf_zonal_mean_calc,  esmf_zonal_mean_reg

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

    call esmf_zonal_mean_reg()

    ! Zonal mean grid for history fields
    zmlons = 0._r8

    zmlat_coord => horiz_coord_create('zmlat', '', nlat, 'latitude', 'degrees_north', 1, nlat, glats)
    zmlon_coord => horiz_coord_create('zmlon', '', 1,    'longitude', 'degrees_east', 1, 1,   zmlons)

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

    call addfld ('T_test',  (/'lev'/), 'A', 'K',  'Zonal-Mean temperature', gridname='physgrid' )
    call addfld ('U_test',  (/'lev'/), 'A', 'm s-1', 'Zonal-Mean zonal wind', gridname='physgrid' )
    call addfld ('V_test',  (/'lev'/), 'A', 'm s-1', 'Zonal-Mean meridianal wind', gridname='physgrid' )

    call addfld ('T_testzm',  (/'lev'/), 'A', 'K',  'Zonal-Mean temperature', gridname='esmf_zonal_mean' )
    call addfld ('U_testzm',  (/'lev'/), 'A', 'm s-1', 'Zonal-Mean zonal wind', gridname='esmf_zonal_mean' )
    call addfld ('V_testzm',  (/'lev'/), 'A', 'm s-1', 'Zonal-Mean meridianal wind', gridname='esmf_zonal_mean' )

  end subroutine zm_test_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine zm_test_run(phys_state)

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    real(r8) :: tfld(pver,pcols,begchunk:endchunk)
    real(r8) :: ufld(pver,pcols,begchunk:endchunk)
    real(r8) :: vfld(pver,pcols,begchunk:endchunk)
    real(r8) :: psfld(pcols,begchunk:endchunk)

    integer :: lchnk, ncol, icol

    real(r8) :: t_zm(lat_beg:lat_end,pver)
    real(r8) :: u_zm(lat_beg:lat_end,pver)
    real(r8) :: v_zm(lat_beg:lat_end,pver)

    real(r8) :: t_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8) :: u_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)
    real(r8) :: v_lonlat(lon_beg:lon_end,lat_beg:lat_end,pver)

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       do icol = 1,ncol
          tfld(:pver,icol,lchnk) = phys_state(lchnk)%t(icol,:pver)
          ufld(:pver,icol,lchnk) = phys_state(lchnk)%u(icol,:pver)
          vfld(:pver,icol,lchnk) = phys_state(lchnk)%v(icol,:pver)
       end do

       call outfld('T_test', phys_state(lchnk)%t, pcols, lchnk)
       call outfld('U_test', phys_state(lchnk)%u, pcols, lchnk)
       call outfld('V_test', phys_state(lchnk)%v, pcols, lchnk)

    end do

    call esmf_phys2lonlat_regrid(tfld, t_lonlat)
    call esmf_zonal_mean_calc(t_lonlat, t_zm)

    call esmf_phys2lonlat_regrid(ufld, u_lonlat)
    call esmf_zonal_mean_calc(u_lonlat, u_zm)

    call esmf_phys2lonlat_regrid(vfld, v_lonlat)
    call esmf_zonal_mean_calc(v_lonlat, v_zm)

    do icol = lat_beg, lat_end
       call outfld('T_testzm',t_zm(icol,:),1,icol)
       call outfld('U_testzm',u_zm(icol,:),1,icol)
       call outfld('V_testzm',v_zm(icol,:),1,icol)
    end do

  end subroutine zm_test_run


end module zm_test_mod
