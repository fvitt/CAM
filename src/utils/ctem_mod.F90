module ctem_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: begchunk, endchunk, pcols, pver, pverp
  use physics_types, only: physics_state
  use phys_grid, only: get_ncols_p
  use spmd_utils, only: masterproc
  use ref_pres, only: pref_mid
  use esmf_lonlat_grid_mod, only: beglon=>lon_beg, endlon=>lon_end, beglat=>lat_beg, endlat=>lat_end
  use cam_history,   only: addfld, outfld
  use perf_mod, only: t_startf, t_stopf
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun

  implicit none

  integer :: ctem_diags_numlats = 0
  logical :: ctem_diags_active = .false.

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine ctem_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils, only : mpicom, masterprocid, mpi_integer, mpi_success
    use string_utils, only : to_lower

    character(len=*), intent(in) :: nlfile
    integer :: unitn, ierr, m, n, ndx_co, ndx_ar

    character(len=*), parameter :: prefix = 'ctem_readnl: '

    namelist /ctem_diags_nl/ ctem_diags_numlats

    if (masterproc) then
       ! read namelist
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'ctem_diags_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, ctem_diags_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(prefix//'ctem_diags_nl: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    call mpi_bcast(ctem_diags_numlats, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(prefix//'mpi_bcast error : ctem_diags_numlats')

    ctem_diags_active = ctem_diags_numlats > 0

    if (masterproc) then
       write(iulog,*) prefix//'ctem_diags_numlats: ', ctem_diags_numlats
       write(iulog,*) prefix//'ctem_diags_active : ', ctem_diags_active
    end if

  end subroutine ctem_readnl

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine ctem_reg()

    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap, cam_grid_register
    use esmf_lonlat_grid_mod, only: glats, nlat, glons, nlon
    use esmf_lonlat_grid_mod, only: esmf_lonlat_grid_init
    use esmf_phys_mesh_mod, only: esmf_phys_mesh_init
    use esmf_phys2lonlat_mod, only: esmf_phys2lonlat_init

    integer :: ind, j, astat

    integer, parameter :: zm_decomp  = 854 ! Must be unique within CAM
    integer, parameter :: reg_decomp = 654

    type(horiz_coord_t), pointer :: zmlon_coord
    type(horiz_coord_t), pointer :: zmlat_coord

    integer(iMap),       pointer :: grid_map(:,:)
    real(r8) :: zmlons(1)

    integer(iMap),       pointer :: coord_map(:) => null()
    type(horiz_coord_t), pointer :: lon_coord
    type(horiz_coord_t), pointer :: lat_coord
    integer :: i

    if (.not.ctem_diags_active) return

    ! initialize grids and mapping
    call esmf_lonlat_grid_init(ctem_diags_numlats)
    call esmf_phys_mesh_init()
    call esmf_phys2lonlat_init()

    ! Zonal mean grid for history fields
    zmlons = 0._r8

    zmlat_coord => horiz_coord_create('zmlat', '', nlat, 'latitude', 'degrees_north', 1, nlat, glats)
    zmlon_coord => horiz_coord_create('zmlon', '', 1,    'longitude', 'degrees_east', 1, 1,   zmlons)

    ! grid decomposition map
    allocate(grid_map(4,endlat-beglat+1), stat=astat)

    ind = 0
    do j = beglat, endlat
       ind = ind + 1
       grid_map(1,ind) = 1
       grid_map(2,ind) = j
       if (beglon==1) then
          grid_map(3,ind) = 1
          grid_map(4,ind) = j
       else
          grid_map(3,ind) = 0
          grid_map(4,ind) = 0
       end if
    end do

    ! register the zonal average grid
    call cam_grid_register('ctem_zm', zm_decomp, zmlat_coord, zmlon_coord, grid_map, &
                           unstruct=.false., zonal_grid=.true.)

    nullify(grid_map)

    ! for the lon-lat grid
    allocate(grid_map(4, ((endlon - beglon + 1) * (endlat - beglat + 1))))
    ind = 0
    do i = beglat, endlat
       do j = beglon, endlon
          ind = ind + 1
          grid_map(1, ind) = j
          grid_map(2, ind) = i
          grid_map(3, ind) = j
          grid_map(4, ind) = i
       end do
    end do

    allocate(coord_map(endlat - beglat + 1))
    if (beglon==1) then
       coord_map = (/ (i, i = beglat, endlat) /)
    else
       coord_map = 0
    end if
    lat_coord => horiz_coord_create('reglat', '', nlat, 'latitude',  'degrees_north', beglat, endlat, &
                                    glats(beglat:endlat),  map=coord_map)

    nullify(coord_map)

    allocate(coord_map(endlon - beglon + 1))
    if (beglat==1) then
       coord_map = (/ (i, i = beglon, endlon) /)
    else
       coord_map = 0
    end if

    lon_coord => horiz_coord_create('reglon', '', nlon, 'longitude',  'degrees_east', beglon, endlon, &
                                    glons(beglon:endlon),  map=coord_map)

    nullify(coord_map)

    call cam_grid_register('ctem_reg', reg_decomp, lat_coord, lon_coord, grid_map, unstruct=.false.)

    nullify(grid_map)

  end subroutine ctem_reg

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine ctem_init()

    if (.not.ctem_diags_active) return

    ! fields on reg lon lat grid
    call addfld ('THreg', (/'lev'/), 'A','K',      'Potential temp', gridname='ctem_reg' )
    call addfld ('Ureg',  (/'lev'/), 'A','m s-1',  'Zonal-Mean zonal wind', gridname='ctem_reg' )
    call addfld ('Vreg',  (/'lev'/), 'A','m s-1',  'Zonal-Mean meridional wind', gridname='ctem_reg' )
    call addfld ('Wreg',  (/'lev'/), 'A','m s-1',  'Zonal-Mean vertical wind', gridname='ctem_reg' )

    call addfld ('VTHreg',(/'lev'/), 'A','K m s-1','Meridional Heat Flux:', gridname='ctem_reg')
    call addfld ('WTHreg',(/'lev'/), 'A','K m s-1','Vertical Heat Flux:', gridname='ctem_reg')
    call addfld ('UVreg', (/'lev'/), 'A','m2 s-2', 'Meridional Flux of Zonal Momentum', gridname='ctem_reg')
    call addfld ('UWreg', (/'lev'/), 'A','m2 s-2', 'Vertical Flux of Zonal Momentum', gridname='ctem_reg')

    ! fields on zonal mean grid
    call addfld ('Uzm',  (/'lev'/), 'A','m s-1',  'Zonal-Mean zonal wind', gridname='ctem_zm' )
    call addfld ('Vzm',  (/'lev'/), 'A','m s-1',  'Zonal-Mean meridional wind', gridname='ctem_zm' )
    call addfld ('Wzm',  (/'lev'/), 'A','m s-1',  'Zonal-Mean vertical wind', gridname='ctem_zm' )
    call addfld ('THzm', (/'lev'/), 'A','K',      'Zonal-Mean potential temp', gridname='ctem_zm' )
    call addfld ('VTHzm',(/'lev'/), 'A','K m s-1','Meridional Heat Flux:', gridname='ctem_zm')
    call addfld ('WTHzm',(/'lev'/), 'A','K m s-1','Vertical Heat Flux:', gridname='ctem_zm')
    call addfld ('UVzm', (/'lev'/), 'A','m2 s-2', 'Meridional Flux of Zonal Momentum', gridname='ctem_zm')
    call addfld ('UWzm', (/'lev'/), 'A','m2 s-2', 'Vertical Flux of Zonal Momentum', gridname='ctem_zm')

  end subroutine ctem_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine ctem_calc(phys_state)
    use air_composition, only: mbarv ! g/mole
    use shr_const_mod, only: rgas => shr_const_rgas ! J/K/kmole
    use shr_const_mod, only: grav => shr_const_g ! m/s2
    use esmf_phys2lonlat_mod, only: esmf_phys2lonlat_regrid
    use esmf_zonal_mean_mod, only: esmf_zonal_mean_calc
    use interpolate_data, only: lininterp

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    real(r8) :: u_phys(pver,pcols,begchunk:endchunk)
    real(r8) :: v_phys(pver,pcols,begchunk:endchunk)
    real(r8) :: w_phys(pver,pcols,begchunk:endchunk)
    real(r8) :: t_phys(pver,pcols,begchunk:endchunk)
    real(r8) :: p_phys(pver,pcols,begchunk:endchunk)

    real(r8) :: u_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: v_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: w_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: t_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: p_lonlat(beglon:endlon,beglat:endlat,pver)

    real(r8) :: ui_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: vi_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: wi_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: ti_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: pi_lonlat(beglon:endlon,beglat:endlat,pver)

    real(r8) :: u_zm(beglat:endlat,pver)
    real(r8) :: v_zm(beglat:endlat,pver)
    real(r8) :: w_zm(beglat:endlat,pver)
    real(r8) :: t_zm(beglat:endlat,pver)

    real(r8) :: ud_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: vd_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: wd_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: td_lonlat(beglon:endlon,beglat:endlat,pver)

    real(r8) :: vtp(beglon:endlon,beglat:endlat,pver)
    real(r8) :: wtp(beglon:endlon,beglat:endlat,pver)
    real(r8) :: uwp(beglon:endlon,beglat:endlat,pver)
    real(r8) :: uvp(beglon:endlon,beglat:endlat,pver)

    real(r8) :: vtp_zm(beglat:endlat,pver)
    real(r8) :: wtp_zm(beglat:endlat,pver)
    real(r8) :: uwp_zm(beglat:endlat,pver)
    real(r8) :: uvp_zm(beglat:endlat,pver)

    integer  :: lchnk, ncol, i, j, k
    real(r8) :: sheight(pver) ! pressure scale height (m)

    real(r8) :: outtmp(beglon:endlon,pver)
    integer :: outcnt

    if (.not.ctem_diags_active) return

    call t_startf('ctem_calc')

    call t_startf('ctem_calc-setarrs')

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       do i = 1,ncol
          ! wind components
          u_phys(:,i,lchnk) =  phys_state(lchnk)%u(i,:)
          v_phys(:,i,lchnk) =  phys_state(lchnk)%v(i,:)

          ! scale height
          sheight(:) = phys_state(lchnk)%t(i,:) * rgas / ( mbarv(i,:,lchnk) * grav ) ! meters

          ! vertical velocity
          w_phys(:,i,lchnk) = -sheight(:) *  phys_state(lchnk)%omega(i,:) / phys_state(lchnk)%pmid(i,:)

          ! potential temperature
          t_phys(:,i,lchnk) = phys_state(lchnk)%t(i,:) * phys_state(lchnk)%exner(i,:)

          ! mid point press
          p_phys(:,i,lchnk) = phys_state(lchnk)%pmid(i,:)
       end do
    end do

    call t_stopf('ctem_calc-setarrs')

    call t_startf('ctem_calc-regrid')

    ! regrid to lon/lat grid
    call esmf_phys2lonlat_regrid(u_phys, u_lonlat)
    call esmf_phys2lonlat_regrid(v_phys, v_lonlat)
    call esmf_phys2lonlat_regrid(w_phys, w_lonlat)
    call esmf_phys2lonlat_regrid(t_phys, t_lonlat)
    call esmf_phys2lonlat_regrid(p_phys, p_lonlat)

    call t_stopf('ctem_calc-regrid')

    outcnt = endlon-beglon+1

    call t_startf('ctem_calc-interp')

    ! vertically intepolate to ref press
    do i = beglon,endlon
       do j = beglat,endlat

          call lininterp( u_lonlat(i,j,:), p_lonlat(i,j,:), pver, &
                          ui_lonlat(i,j,:), pref_mid(:), pver )

          call lininterp( v_lonlat(i,j,:), p_lonlat(i,j,:), pver, &
                          vi_lonlat(i,j,:), pref_mid(:), pver )

          call lininterp( w_lonlat(i,j,:), p_lonlat(i,j,:), pver, &
                          wi_lonlat(i,j,:), pref_mid(:), pver )

          call lininterp( t_lonlat(i,j,:), p_lonlat(i,j,:), pver, &
                          ti_lonlat(i,j,:), pref_mid(:), pver )

       end do
    end do

    call t_stopf('ctem_calc-interp')

    call t_startf('ctem_calc-zonal_mean-uvwt')

    ! calculate zonal means from interpolated fields
    call esmf_zonal_mean_calc(ui_lonlat, u_zm)
    call esmf_zonal_mean_calc(vi_lonlat, v_zm)
    call esmf_zonal_mean_calc(wi_lonlat, w_zm)
    call esmf_zonal_mean_calc(ti_lonlat, t_zm)

    call t_stopf('ctem_calc-zonal_mean-uvwt')

    call t_startf('ctem_calc-calc_deviations')

    ! Calculate zonal deviations from zonal means
    do j = beglat,endlat
       do i = beglon,endlon
          ud_lonlat(i,j,:) = ui_lonlat(i,j,:) - u_zm(j,:)
          vd_lonlat(i,j,:) = vi_lonlat(i,j,:) - v_zm(j,:)
          wd_lonlat(i,j,:) = wi_lonlat(i,j,:) - w_zm(j,:)
          td_lonlat(i,j,:) = ti_lonlat(i,j,:) - t_zm(j,:)
       end do
    end do

    call t_stopf('ctem_calc-calc_deviations')

    call t_startf('ctem_calc-calc_fluxes')

    ! Calculate fluxes
    vtp(:,:,:) = vd_lonlat(:,:,:) * td_lonlat(:,:,:)
    wtp(:,:,:) = wd_lonlat(:,:,:) * td_lonlat(:,:,:)
    uwp(:,:,:) = ud_lonlat(:,:,:) * wd_lonlat(:,:,:)
    uvp(:,:,:) = ud_lonlat(:,:,:) * vd_lonlat(:,:,:)

    call t_stopf('ctem_calc-calc_fluxes')

    call t_startf('ctem_calc-zonal_mean-p')

    call esmf_zonal_mean_calc(vtp, vtp_zm)
    call esmf_zonal_mean_calc(wtp, wtp_zm)
    call esmf_zonal_mean_calc(uwp, uwp_zm)
    call esmf_zonal_mean_calc(uvp, uvp_zm)

    call t_stopf('ctem_calc-zonal_mean-p')

    call t_startf('ctem_calc-output')

    ! output diagnostics
    do j = beglat,endlat
       outtmp(beglon:endlon,1:pver) = ti_lonlat(beglon:endlon,j,1:pver)
       call outfld('THreg',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = ui_lonlat(beglon:endlon,j,1:pver)
       call outfld('Ureg',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = vi_lonlat(beglon:endlon,j,1:pver)
       call outfld('Vreg',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = wi_lonlat(beglon:endlon,j,1:pver)
       call outfld('Wreg',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = vtp(beglon:endlon,j,1:pver)
       call outfld('VTHreg',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = wtp(beglon:endlon,j,1:pver)
       call outfld('WTHreg',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = uvp(beglon:endlon,j,1:pver)
       call outfld('UVreg',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = uwp(beglon:endlon,j,1:pver)
       call outfld('UWreg',outtmp, outcnt, j)

       call outfld('Uzm',  u_zm(j,:), 1,j)
       call outfld('Vzm',  v_zm(j,:), 1,j)
       call outfld('Wzm',  w_zm(j,:), 1,j)
       call outfld('THzm', t_zm(j,:), 1,j)

       call outfld('VTHzm',vtp_zm(j,:),1,j)
       call outfld('WTHzm',wtp_zm(j,:),1,j)
       call outfld('UVzm', uvp_zm(j,:),1,j)
       call outfld('UWzm', uwp_zm(j,:),1,j)
    end do

    call t_stopf('ctem_calc-output')

    call t_stopf('ctem_calc')

  end subroutine ctem_calc

end module ctem_mod
