module phys_grid_ctem
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: begchunk, endchunk, pcols, pver, pverp
  use ref_pres,     only: pref_edge
  use interpolate_data, only: vertinterp
  use physics_types,only: physics_state
  use cam_history,  only: addfld, outfld

  use zmean_fields, only: zmean_fields_init, zmean_3d

  implicit none

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine phys_grid_ctem_init

    call addfld ('VTHzmean',(/ 'ilev' /),'A','MK/S', 'Meridional Heat Flux: 3D zon. mean', gridname='physgrid' )
    call addfld ('WTHzmean',(/ 'ilev' /),'A','MK/S', 'Vertical Heat Flux: 3D zon. mean', gridname='physgrid' )
    call addfld ('UVzmean', (/ 'ilev' /),'A','M2/S2','Meridional Flux of Zonal Momentum: 3D zon. mean', gridname='physgrid' )
    call addfld ('UWzmean', (/ 'ilev' /),'A','M2/S2','Vertical Flux of Zonal Momentum: 3D zon. mean', gridname='physgrid' )

    call addfld ('Uzmean',  (/ 'ilev' /),'A','M/S',  'Zonal-Mean zonal wind - defined on ilev', gridname='physgrid')
    call addfld ('Vzmean',  (/ 'ilev' /),'A','M/S',  'Zonal-Mean meridional wind - defined on ilev', gridname='physgrid' )
    call addfld ('Wzmean',  (/ 'ilev' /),'A','M/S',  'Zonal-Mean vertical wind - defined on ilev', gridname='physgrid' )
    call addfld ('THzmean', (/ 'ilev' /),'A',  'K',  'Zonal-Mean potential temp - defined on ilev', gridname='physgrid' )
    call addfld ('THphys', (/ 'ilev' /),'A',  'K',  'Zonal-Mean potential temp - defined on ilev', gridname='physgrid' )

  end subroutine phys_grid_ctem_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine phys_grid_ctem_diags(phys_state)
    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    real(r8) :: ui(pcols,pverp,begchunk:endchunk)
    real(r8) :: vi(pcols,pverp,begchunk:endchunk)
    real(r8) :: wi(pcols,pverp,begchunk:endchunk)

    real(r8) :: uzm(pcols,pverp,begchunk:endchunk)
    real(r8) :: vzm(pcols,pverp,begchunk:endchunk)
    real(r8) :: wzm(pcols,pverp,begchunk:endchunk)

    real(r8) :: ud(pcols,pverp,begchunk:endchunk)
    real(r8) :: vd(pcols,pverp,begchunk:endchunk)
    real(r8) :: wd(pcols,pverp,begchunk:endchunk)
    real(r8) :: thd(pcols,pverp,begchunk:endchunk)

    real(r8) :: uvp(pcols,pverp,begchunk:endchunk)
    real(r8) :: uwp(pcols,pverp,begchunk:endchunk)
    real(r8) :: vthp(pcols,pverp,begchunk:endchunk)
    real(r8) :: wthp(pcols,pverp,begchunk:endchunk)

    real(r8) :: uv(pcols,pverp,begchunk:endchunk)
    real(r8) :: uw(pcols,pverp,begchunk:endchunk)
    real(r8) :: vth(pcols,pverp,begchunk:endchunk)
    real(r8) :: wth(pcols,pverp,begchunk:endchunk)

    integer :: lchnk, ncol, k
    real(r8) :: fld_tmp(pcols,pverp)

    real(r8) :: theta(pcols,pver,begchunk:endchunk) ! potential temperature
    real(r8) :: thi(pcols,pverp,begchunk:endchunk)
    real(r8) :: thzm(pcols,pverp,begchunk:endchunk)

    real(r8) :: w(pcols,pver,begchunk:endchunk)

    real(r8), parameter :: hscale = 7000._r8          ! pressure scale height

    ui(:,:,:) = 0._r8
    vi(:,:,:) = 0._r8
    wi(:,:,:) = 0._r8
    thi(:,:,:) = 0._r8

    uzm(:,:,:) = 0._r8
    vzm(:,:,:) = 0._r8
    wzm(:,:,:) = 0._r8
    thzm(:,:,:) = 0._r8

    ud(:,:,:) = 0._r8
    vd(:,:,:) = 0._r8
    uvp(:,:,:) = 0._r8

    do lchnk = begchunk,endchunk

       ncol = phys_state(lchnk)%ncol

       theta(:ncol,:,lchnk) = phys_state(lchnk)%t(:ncol,:) * phys_state(lchnk)%exner(:ncol,:)
       w(:ncol,:,lchnk)  = - hscale *  phys_state(lchnk)%omega(:ncol,:) / phys_state(lchnk)%pmid(:ncol,:)

       do k = 1,pverp
          call vertinterp( ncol, pcols, pver, phys_state(lchnk)%pmid(:,:), pref_edge(k), phys_state(lchnk)%u(:,:), ui(:,k,lchnk) )
          call vertinterp( ncol, pcols, pver, phys_state(lchnk)%pmid(:,:), pref_edge(k), phys_state(lchnk)%v(:,:), vi(:,k,lchnk) )
          call vertinterp( ncol, pcols, pver, phys_state(lchnk)%pmid(:,:), pref_edge(k), theta(:,:,lchnk), thi(:,k,lchnk) )
          call vertinterp( ncol, pcols, pver, phys_state(lchnk)%pmid(:,:), pref_edge(k), w(:,:,lchnk), wi(:,k,lchnk) )
       end do

    end do

    ! these need to be evaluated on the physics grid (3D)
    ! to be used in the deviations calculation below
    uzm(:,:,:) = zmean_3d(ui(:,:,:), nlev=pverp)
    vzm(:,:,:) = zmean_3d(vi(:,:,:), nlev=pverp)
    wzm(:,:,:) = zmean_3d(wi(:,:,:), nlev=pverp)
    thzm(:,:,:) = zmean_3d(thi(:,:,:), nlev=pverp)

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       do k = 1,pverp
          ! zonal deviations
          thd(:ncol,k,lchnk) = thi(:ncol,k,lchnk) - thzm(:ncol,k,lchnk)
          ud(:ncol,k,lchnk) = ui(:ncol,k,lchnk) - uzm(:ncol,k,lchnk)
          vd(:ncol,k,lchnk) = vi(:ncol,k,lchnk) - vzm(:ncol,k,lchnk)
          wd(:ncol,k,lchnk) = wi(:ncol,k,lchnk) - wzm(:ncol,k,lchnk)
          ! fluxes
          uvp(:ncol,k,lchnk) = ud(:ncol,k,lchnk) * vd(:ncol,k,lchnk)
          uwp(:ncol,k,lchnk) = ud(:ncol,k,lchnk) * wd(:ncol,k,lchnk)
          vthp(:ncol,k,lchnk) = vd(:ncol,k,lchnk) * thd(:ncol,k,lchnk)
          wthp(:ncol,k,lchnk) = wd(:ncol,k,lchnk) * thd(:ncol,k,lchnk)
       end do
    end do

    ! idealy want to evailuate and output thsese on a zonal-mean grid (2D rather than 3D)
    uv(:,:,:) = zmean_3d(uvp(:,:,:), nlev=pverp)
    uw(:,:,:) = zmean_3d(uwp(:,:,:), nlev=pverp)
    vth(:,:,:) = zmean_3d(vthp(:,:,:), nlev=pverp)
    wth(:,:,:) = zmean_3d(wthp(:,:,:), nlev=pverp)

    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol

       fld_tmp(:ncol,:) = thi(:ncol,:,lchnk)
       call outfld( 'THphys', fld_tmp(:ncol,:), ncol, lchnk)

       fld_tmp(:ncol,:) = thzm(:ncol,:,lchnk)
       call outfld( 'THzmean', fld_tmp(:ncol,:), ncol, lchnk)


       fld_tmp(:ncol,:) = vth(:ncol,:,lchnk)
       call outfld( 'VTHzmean', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = wth(:ncol,:,lchnk)
       call outfld( 'WTHzmean', fld_tmp(:ncol,:), ncol, lchnk)

       fld_tmp(:ncol,:) = uv(:ncol,:,lchnk)
       call outfld( 'UVzmean', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = uw(:ncol,:,lchnk)
       call outfld( 'UWzmean', fld_tmp(:ncol,:), ncol, lchnk)

       fld_tmp(:ncol,:) = uzm(:ncol,:,lchnk)
       call outfld( 'Uzmean', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = vzm(:ncol,:,lchnk)
       call outfld( 'Vzmean', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = wzm(:ncol,:,lchnk)
       call outfld( 'Wzmean', fld_tmp(:ncol,:), ncol, lchnk)
    end do

  end subroutine phys_grid_ctem_diags

end module phys_grid_ctem
