module edyn3D_driver
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils, only: masterproc
  use cam_abortutils, only: endrun

  use edyn3D_maggrid, only: gen_highres_grid, edyn3D_gen_ggj_grid, edyn3D_gen_qd_grid, &
                            edyn3D_gen_geo_grid, edyn3D_qcoef, edyn3D_calculate_mf

  use edyn3D_mpi, only: mp_init_edyn3D, mp_distribute_mag_edyn3D, mp_exchange_tasks_edyn3D
  use edyn3D_mpi, only: mytid, ntask
  use edyn3D_fieldline, only: fieldline_init, fieldline_getapex
  use edyn3D_params, only: nmlon,nmlonp1,nmlat_h,nptsp_total,nptss1_total,nptss2_total, &
                           ylonm,ylonm_s,nhgt_fix,nmlat_T1,nmlatS2_h

  use edyn3D_esmf_regrid
  use edyn3D_regridder

  use perf_mod, only: t_startf, t_stopf

  use physconst, only: pi

  implicit none

  private
  public :: edyn3D_driver_reg
  public :: edyn3D_driver_timestep

  real(r8), parameter :: r2d = 180._r8/pi

contains
  subroutine edyn3D_driver_reg(mpicom, npes)
    use cam_history,         only: addfld, horiz_only
    use mo_apex,             only: mo_apex_init1,geomag_year
    use edyn3D_fline_fields, only: edyn3D_fline_fields_alloc

    integer, intent(in) :: mpicom, npes

    call mo_apex_init1()

    call mp_init_edyn3D(mpicom, npes)

    call gen_highres_grid()

    call edyn3D_gen_ggj_grid()

    call edyn3D_gen_qd_grid()

    call edyn3D_gen_geo_grid()

    call edyn3D_qcoef()

    call edyn3D_calculate_mf()

    call mp_distribute_mag_edyn3D(nmlon)

    call mp_exchange_tasks_edyn3D(mpicom, iprint=0)

    call fieldline_init()  ! Allocate and populate the p, r, s1, and s2 field line structures for computations

    call fieldline_getapex()

    call reg_hist_grid()

    call edyn3D_esmf_regrid_init()

!    call apxparm(geomag_year)

    call addfld ('Tn_mag', horiz_only, 'I', 'K','Neutral Temperature on geo-magnetic field line grid', &
                  gridname='magfline_p')

    call addfld ('GEOGALTp', horiz_only, 'I', 'km','magnetic field line point geographic (geodetic) altitude', &
                  gridname='magfline_p')
    call addfld ('GEOGLATp', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) latitude', &
                  gridname='magfline_p')
    call addfld ('GEOGLONp', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) longitude', &
                  gridname='magfline_p')

    call addfld ('GEOGALTs1', horiz_only, 'I', 'km','magnetic field line point geographic (geodetic) altitude', &
                  gridname='magfline_s1')
    call addfld ('GEOGLATs1', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) latitude', &
                  gridname='magfline_s1')
    call addfld ('GEOGLONs1', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) longitude', &
                  gridname='magfline_s1')

    call addfld ('GEOGALTs2', horiz_only, 'I', 'km','magnetic field line point geographic (geodetic) altitude', &
                  gridname='magfline_s2')
    call addfld ('GEOGLATs2', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) latitude', &
                  gridname='magfline_s2')
    call addfld ('GEOGLONs2', horiz_only, 'I', 'Degrees','magnetic field line point geographic (geodetic) longitude', &
                  gridname='magfline_s2')

    call addfld ('height_s1', horiz_only, 'I', 'm','altitude', &
                  gridname='magfline_s1')
    call addfld ('height_s2', horiz_only, 'I', 'm','altitude', &
                  gridname='magfline_s2')

    call addfld ('sigma_ped_s1', horiz_only, 'I', 'K','Ped cond. on S1 mag field line grid', &
                  gridname='magfline_s1')
    call addfld ('sigma_hal_s1', horiz_only, 'I', 'K','Hal cond. on S1 mag field line grid', &
                  gridname='magfline_s1')
    call addfld ('sigma_ped_s2', horiz_only, 'I', 'K','Ped cond. on S2 mag field line grid', &
                  gridname='magfline_s2')
    call addfld ('sigma_hal_s2', horiz_only, 'I', 'K','Hal cond. on S2 mag field line grid', &
                  gridname='magfline_s2')

    call addfld ('un_s1', horiz_only, 'I', 'm/s','Zonal wind on S1 mag field line grid', &
                  gridname='magfline_s1')
    call addfld ('vn_s1', horiz_only, 'I', 'm/s','Meridional wind on S1 mag field line grid', &
                  gridname='magfline_s1')
    call addfld ('un_s2', horiz_only, 'I', 'm/s','Zonal wind on S2 mag field line grid', &
                  gridname='magfline_s2')
    call addfld ('vn_s2', horiz_only, 'I', 'm/s','Meridional wind on S2 mag field line grid', &
                  gridname='magfline_s2')

    call addfld ('Tn_opg0', (/ 'lev' /), 'I', 'K','Tn_opg0 test field' , gridname='geo_grid')
    call addfld ('Tn_opg1', (/ 'lev' /), 'I', 'K','Tn_opg1 test field' , gridname='geo_grid')

    call addfld ('POTENp', horiz_only, 'I', 'Volts','magnetic field line point electric potential', gridname='magfline_p')

    call edyn3D_fline_fields_alloc()

    call addfld ('ED1s1', horiz_only, 'I', 'V/m','Eastward electric field on s1 grid', gridname='magfline_s1')
    call addfld ('ED2s1', horiz_only, 'I', 'V/m','Equatorward electric field on s1 grid', gridname='magfline_s1')

    call addfld ('ED1s2', horiz_only, 'I', 'V/m','Eastward electric field on s2 grid', gridname='magfline_s2')
    call addfld ('ED2s2', horiz_only, 'I', 'V/m','Equatorward electric field on s2 grid', gridname='magfline_s2')

    call addfld ('Ve1s1', horiz_only, 'I', 'm/s','Ion Drift Velocity 1 on s1 grid', gridname='magfline_s1')
    call addfld ('Ve2s1', horiz_only, 'I', 'm/s','Ion Drift Velocity 2 on s1 grid', gridname='magfline_s1')

    call addfld ('Ve1s2', horiz_only, 'I', 'm/s','Ion Drift Velocity 1 on s2 grid', gridname='magfline_s2')
    call addfld ('Ve2s2', horiz_only, 'I', 'm/s','Ion Drift Velocity 2 on s2 grid', gridname='magfline_s2')

    call addfld ('UI_s1', horiz_only, 'I', 'm/s','Zonal Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('VI_s1', horiz_only, 'I', 'm/s','Meridional Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('WI_s1', horiz_only, 'I', 'm/s','Vertical Ion Drift Velocity on s1 grid', gridname='magfline_s1')

    call addfld ('IonU_s1', horiz_only, 'I', 'm/s','Zonal Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('IonV_s1', horiz_only, 'I', 'm/s','Meridional Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('IonW_s1', horiz_only, 'I', 'm/s','Verical Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('IonU_opg', (/ 'lev' /), 'I', 'm/s','Zonal Ion Drift Velocity on oplus grid' , gridname='geo_grid')
    call addfld ('IonV_opg', (/ 'lev' /), 'I', 'm/s','Meridional Ion Drift Velocity on oplus grid' , gridname='geo_grid')
    call addfld ('IonW_opg', (/ 'lev' /), 'I', 'm/s','Vertical Ion Drift Velocity on oplus grid' , gridname='geo_grid')

  end subroutine edyn3D_driver_reg

  subroutine edyn3D_driver_timestep( nphyscol, nphyslev, physalt, tn, sigPed, sigHal, un, vn, tn_out, tn_out2, ui_out, vi_out, wi_out )

    use edyn3d_mpi, only: mlon0_p,mlon1_p
    use cam_history,  only: outfld
    use edyn3D_fieldline, only: fline_p, fline_s1, fline_s2
    use edyn3D_fline_fields, only: Tn_p, height_s1, height_s2, IonV_s1, IonU_s1, IonW_s1
    use edyn3D_fline_fields, only: sigma_ped_s1,sigma_hal_s1,sigma_ped_s2,sigma_hal_s2,un_s1,vn_s1,un_s2,vn_s2
    use edyn_mpi, only: lon0,lon1,lat0,lat1,lev0,lev1
    use regridder, only: regrid_phys2geo_3d, regrid_geo2phys_3d
    use edyn3D_calc_coef_fac_const_rhs, only: edyn3D_calc_coef,edyn3D_calc_FAC,edyn3D_add_coef_ns, &
                                              edyn3D_gather_coef_ns,edyn3D_const_rhs,edyn3D_scatter_poten

    integer,  intent(in) :: nphyscol, nphyslev
    real(r8), intent(in) :: physalt(nphyslev,nphyscol)
    real(r8), intent(in) :: tn(nphyslev,nphyscol)
    real(r8), intent(in) :: sigPed(nphyslev,nphyscol)
    real(r8), intent(in) :: sigHal(nphyslev,nphyscol)
    real(r8), intent(in) :: un(nphyslev,nphyscol)
    real(r8), intent(in) :: vn(nphyslev,nphyscol)
    real(r8), intent(out) :: tn_out(nphyslev,nphyscol)
    real(r8), intent(out) :: tn_out2(nphyslev,nphyscol)
    real(r8), intent(out) :: ui_out(nphyslev,nphyscol)
    real(r8), intent(out) :: vi_out(nphyslev,nphyscol)
    real(r8), intent(out) :: wi_out(nphyslev,nphyscol)

    real(r8) :: geogaltp(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglatp(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglonp(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geogalts1(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglats1(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geoglons1(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: geogalts2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: geoglats2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: geoglons2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: Tn_tmp(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: potential(mlon0_p:mlon1_p, nptsp_total)
    real(r8) :: ed1_s1(mlon0_p:mlon1_p, nptss1_total)
    real(r8) :: ed2_s1(mlon0_p:mlon1_p, nptss1_total)
    real(r8) :: ed1_s2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: ed2_s2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: ve1_s1(mlon0_p:mlon1_p, nptss1_total)
    real(r8) :: ve2_s1(mlon0_p:mlon1_p, nptss1_total)
    real(r8) :: ve1_s2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: ve2_s2(mlon0_p:mlon1_p, nptss2_total)
    real(r8) :: ui_s1(mlon0_p:mlon1_p, nptss1_total)
    real(r8) :: vi_s1(mlon0_p:mlon1_p, nptss1_total)
    real(r8) :: wi_s1(mlon0_p:mlon1_p, nptss1_total)

    integer,parameter :: ndyn = 7
    real(r8) :: efield_fline(ndyn,nhgt_fix,nmlat_T1,mlon0_p:mlon1_p)

    integer :: i,j,jj,k,isn,ncnt,ncnt1,ncnt2,ncnt3,ier
    integer :: dk,k0,k1
    integer :: npts_s1

    real(r8) :: opalt (lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8) :: Tn_oplus0(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8) :: Tn_oplus1(lon0:lon1,lat0:lat1,lev0:lev1)

    real(r8) :: IonU_oplus(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8) :: IonV_oplus(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8) :: IonW_oplus(lon0:lon1,lat0:lat1,lev0:lev1)

    call t_startf('edyn3D_driver_timestep.1')
    call edyn3D_regridder_phys2mag(physalt,physalt,nphyscol,nphyslev,height_s1)
    call edyn3D_regridder_phys2mag(physalt,physalt,nphyscol,nphyslev,height_s2)

    call output_fline_field(height_s1)
    call output_fline_field(height_s2)

    call edyn3D_regridder_phys2mag(sigPed,physalt,nphyscol,nphyslev,sigma_ped_s1)
    call edyn3D_regridder_phys2mag(sigPed,physalt,nphyscol,nphyslev,sigma_ped_s2)

    call output_fline_field(sigma_ped_s1)
    call output_fline_field(sigma_ped_s2)

    call edyn3D_regridder_phys2mag(sigHal,physalt,nphyscol,nphyslev,sigma_hal_s1)
    call edyn3D_regridder_phys2mag(sigHal,physalt,nphyscol,nphyslev,sigma_hal_s2)

    call output_fline_field(sigma_hal_s1)
    call output_fline_field(sigma_hal_s2)

    call edyn3D_regridder_phys2mag(un,physalt,nphyscol,nphyslev,un_s1)
    call edyn3D_regridder_phys2mag(un,physalt,nphyscol,nphyslev,un_s2)

    call output_fline_field(un_s1)
    call output_fline_field(un_s2)

    call edyn3D_regridder_phys2mag(vn,physalt,nphyscol,nphyslev,vn_s1)
    call edyn3D_regridder_phys2mag(vn,physalt,nphyscol,nphyslev,vn_s2)

    call output_fline_field(vn_s1)
    call output_fline_field(vn_s2)

    if (mytid<ntask) then
       geogaltp=-huge(1._r8)
       geoglatp=-huge(1._r8)
       geoglonp=-huge(1._r8)

       geogalts1=-huge(1._r8)
       geoglats1=-huge(1._r8)
       geoglons1=-huge(1._r8)

       geogalts2=-huge(1._r8)
       geoglats2=-huge(1._r8)
       geoglons2=-huge(1._r8)

       do i = mlon0_p,mlon1_p
          ncnt = 0
          do j = 1,nmlat_h
             do isn = 1,2

                fline_s1(i,j,isn)%sigP(:) = sigma_ped_s1%flines(i,j,isn)%fld(:)
                fline_s1(i,j,isn)%sigH(:) = sigma_hal_s1%flines(i,j,isn)%fld(:)

                fline_s1(i,j,isn)%un(:) = un_s1%flines(i,j,isn)%fld(:)
                fline_s1(i,j,isn)%vn(:) = vn_s1%flines(i,j,isn)%fld(:)

                if (isn==1) then
                   k0 = 1
                   k1 = fline_p(i,j,isn)%npts
                   dk = 1
                else
                   k0 = fline_p(i,j,isn)%npts
                   k1 = 1
                   dk = -1
                endif

                do k = k0,k1,dk
                   ncnt = ncnt + 1
                   geogaltp(i,ncnt) = fline_p(i,j,isn)%hgt_pt(k)
                   geoglatp(i,ncnt) = fline_p(i,j,isn)%glat(k)
                   geoglonp(i,ncnt) = fline_p(i,j,isn)%glon(k)
                   geogalts1(i,ncnt) = fline_s1(i,j,isn)%hgt_pt(k)
                   geoglats1(i,ncnt) = fline_s1(i,j,isn)%glat(k)
                   geoglons1(i,ncnt) = fline_s1(i,j,isn)%glon(k)
                end do
             end do
          end do
          ncnt = 0
          do j = 1,nmlats2_h
             do isn = 1,2

                fline_s2(i,j,isn)%sigP(:) = sigma_ped_s2%flines(i,j,isn)%fld(:)
                fline_s2(i,j,isn)%sigH(:) = sigma_hal_s2%flines(i,j,isn)%fld(:)
                fline_s2(i,j,isn)%un(:) = un_s2%flines(i,j,isn)%fld(:)
                fline_s2(i,j,isn)%vn(:) = vn_s2%flines(i,j,isn)%fld(:)

                if (isn==1) then
                   k0 = 1
                   k1 = fline_s2(i,j,isn)%npts
                   dk = 1
                else
                   k0 = fline_s2(i,j,isn)%npts
                   k1 = 1
                   dk = -1
                endif

                do k = k0,k1,dk
                   ncnt = ncnt + 1
                   geogalts2(i,ncnt) = fline_s2(i,j,isn)%hgt_pt(k)
                   geoglats2(i,ncnt) = fline_s2(i,j,isn)%glat(k)
                   geoglons2(i,ncnt) = fline_s2(i,j,isn)%glon(k)
               end do
             end do
          end do
       end do

       do j = 1,nptsp_total
          call outfld('GEOGALTp', geogaltp(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLATp', geoglatp(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLONp', geoglonp(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGALTs1', geogalts1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLATs1', geoglats1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLONs1', geoglons1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

       do j = 1,nptss2_total
          call outfld('GEOGALTs2', geogalts2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLATs2', geoglats2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('GEOGLONs2', geoglons2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

    end if

    call edyn3D_regridder_phys2mag(tn,physalt,nphyscol,nphyslev,Tn_p)

    call output_fline_field(Tn_p)

    call t_stopf('edyn3D_driver_timestep.1')

    !
    ! Call 3D dynamo routine for solving
    !
    if (mytid<ntask) then

      call t_startf('edyn3D_driver_timestep.2')

      call edyn3D_get_conduct     ! - Get conductivities for edyn3D_calc_FAC

      call edyn3D_calc_mn_s1s2    ! - (calc conductivities, each timestep)

      call edyn3D_calc_je_s1s2    ! - (must be calculated each timestep)

      call edyn3D_coef_halos      ! - Get halo points required in following routines

      call edyn3D_calc_S          ! - (must be calculated each timestep)

      call edyn3D_calc_coef       ! - calc LHS & RHS

      call edyn3D_calc_FAC        ! - calc high latitude

      call edyn3D_add_coef_ns     ! - add North & South coef

      call t_stopf('edyn3D_driver_timestep.2')

      call t_startf('edyn3D_driver_timestep.3.gather')

      call edyn3D_gather_coef_ns  ! - gather coef_ns for solver

      call t_stopf('edyn3D_driver_timestep.3.gather')

      call t_startf('edyn3D_driver_timestep.4.solve')
      if (mytid == 0) then
        call edyn3D_const_rhs     ! - solver - solve for rhs (electric potential)
      endif
      call t_stopf('edyn3D_driver_timestep.4.solve')

      call t_startf('edyn3D_driver_timestep.5.scatter')
      call edyn3D_scatter_poten   ! - Send global potential to each task
      call t_stopf('edyn3D_driver_timestep.5.scatter')

      call t_startf('edyn3D_driver_timestep.6')
      call edyn3D_poten_halos     ! - Get potential halo points required in next call

      call edyn3D_calc_efield    ! - Calculate the electric field and ion drift velocities
      call t_stopf('edyn3D_driver_timestep.6')

    endif

    call t_startf('edyn3D_driver_timestep.7')

    call edyn3D_regridder_mag2phys(Tn_p, physalt, nphyscol,nphyslev, tn_out)

    call regrid_phys2geo_3d( Tn, Tn_oplus0, nphyslev, 1, nphyscol )
    do j = lat0,lat1
       call outfld( 'Tn_opg0', Tn_oplus0(lon0:lon1,j,lev0:lev1), lon1-lon0+1, j )
    end do

    call regrid_phys2geo_3d( physalt, opalt, nphyslev, 1, nphyscol )
    call edyn3D_regridder_mag2oplus( opalt, Tn_p, Tn_oplus1 )
    do j = lat0,lat1
       call outfld( 'Tn_opg1', Tn_oplus1(lon0:lon1,j,lev0:lev1), lon1-lon0+1, j )
    end do

    call regrid_geo2phys_3d( Tn_oplus1, Tn_out2, nphyslev, 1, nphyscol )

    !  diagnostics ...

    proc_tasks: if (mytid<ntask) then

       do i = mlon0_p,mlon1_p
          ncnt1 = 0
          do j = 1,nmlat_h
             do isn = 1,2

                if (isn==1) then
                   k0 = 1
                   k1 = fline_p(i,j,isn)%npts
                   dk = 1
                else
                   k0 = fline_p(i,j,isn)%npts
                   k1 = 1
                   dk = -1
                endif

                do k = k0,k1,dk
                   ncnt1 = ncnt1 + 1
                   potential(i,ncnt1) = fline_p(i,j,isn)%pot
                end do

             end do
          end do
       end do

       do i = mlon0_p,mlon1_p
          ncnt2 = 0
          do j = 1,nmlat_h
             do isn = 1,2

                if (isn==1) then
                   k0 = 1
                   k1 = fline_s1(i,j,isn)%npts
                   dk = 1
                else
                   k0 = fline_s1(i,j,isn)%npts
                   k1 = 1
                   dk = -1
                endif

                do k = k0,k1,dk
                   ncnt2 = ncnt2 + 1
                   ed1_s1(i,ncnt2) = fline_s1(i,j,isn)%ed1
                   ed2_s1(i,ncnt2) = fline_s1(i,j,isn)%ed2
                   ve1_s1(i,ncnt2) = fline_s1(i,j,isn)%ve1
                   ve2_s1(i,ncnt2) = fline_s1(i,j,isn)%ve2
                   ui_s1(i,ncnt2)  = fline_s1(i,j,isn)%ve1*fline_s1(i,j,isn)%e1(1,k)+ &
                                     fline_s1(i,j,isn)%ve2*fline_s1(i,j,isn)%e2(1,k)
                   vi_s1(i,ncnt2)  = fline_s1(i,j,isn)%ve1*fline_s1(i,j,isn)%e1(2,k)+ &
                                     fline_s1(i,j,isn)%ve2*fline_s1(i,j,isn)%e2(2,k)
                   wi_s1(i,ncnt2)  = fline_s1(i,j,isn)%ve1*fline_s1(i,j,isn)%e1(3,k)+ &
                                     fline_s1(i,j,isn)%ve2*fline_s1(i,j,isn)%e2(3,k)

                   IonU_s1%flines(i,j,isn)%fld(k) = ui_s1(i,ncnt2)
                   IonV_s1%flines(i,j,isn)%fld(k) = vi_s1(i,ncnt2)
                   IonW_s1%flines(i,j,isn)%fld(k) = wi_s1(i,ncnt2)

                end do

             end do
          end do
          if (i == mlon0_p) nptss1_total = ncnt2
       end do

       call output_fline_field(IonU_s1)
       call output_fline_field(IonV_s1)
       call output_fline_field(IonW_s1)

       call edyn3D_regridder_mag2oplus( opalt, IonU_s1, IonU_oplus )
       call edyn3D_regridder_mag2oplus( opalt, IonV_s1, IonV_oplus )
       call edyn3D_regridder_mag2oplus( opalt, IonW_s1, IonW_oplus )

       do j = lat0,lat1
          call outfld( 'IonU_opg', IonU_oplus(lon0:lon1,j,lev0:lev1), lon1-lon0+1, j )
          call outfld( 'IonV_opg', IonV_oplus(lon0:lon1,j,lev0:lev1), lon1-lon0+1, j )
          call outfld( 'IonW_opg', IonW_oplus(lon0:lon1,j,lev0:lev1), lon1-lon0+1, j )
       end do

       do i = mlon0_p,mlon1_p
          ncnt3 = 0
          do j = 1,nmlatS2_h
             do isn = 1,2

		if (isn==1) then
		   k0 = 1
		   k1 = fline_s2(i,j,isn)%npts
		   dk = 1
		else
		   k0 = fline_s2(i,j,isn)%npts
		   k1 = 1
		   dk = -1
		endif

		do k = k0,k1,dk
		   ncnt3 = ncnt3 + 1
		   ed1_s2(i,ncnt3) = fline_s2(i,j,isn)%ed1
		   ed2_s2(i,ncnt3) = fline_s2(i,j,isn)%ed2
		   ve1_s2(i,ncnt3) = fline_s2(i,j,isn)%ve1
		   ve2_s2(i,ncnt3) = fline_s2(i,j,isn)%ve2

	        end do

             end do
          end do
       end do

       do j = 1,nptsp_total
          call outfld('POTENp',  potential(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

       do j = 1,nptss1_total
          call outfld('ED1s1',  ed1_s1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('ED2s1',  ed2_s1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

       do j = 1,nptss2_total
          call outfld('ED1s2',  ed1_s2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('ED2s2',  ed2_s2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

       do j = 1,nptss1_total
          call outfld('Ve1s1',  ve1_s1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('Ve2s1',  ve2_s1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

       do j = 1,nptss1_total
          call outfld('UI_s1',  ui_s1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('VI_s1',  vi_s1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('WI_s1',  wi_s1(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

       do j = 1,nptss2_total
          call outfld('Ve1s2',  ve1_s2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
          call outfld('Ve2s2',  ve2_s2(mlon0_p:mlon1_p,j), mlon1_p-mlon0_p+1, j)
       end do

    end if proc_tasks

    call edyn3D_regridder_mag2phys(IonU_s1, physalt, nphyscol,nphyslev, ui_out)
    call edyn3D_regridder_mag2phys(IonV_s1, physalt, nphyscol,nphyslev, vi_out)
    call edyn3D_regridder_mag2phys(IonW_s1, physalt, nphyscol,nphyslev, wi_out)

    call t_stopf('edyn3D_driver_timestep.7')

  end subroutine edyn3D_driver_timestep

  subroutine reg_hist_grid

    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
    use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
    use edyn3d_mpi, only: mlon0_p,mlon1_p
    use edyn3D_fieldline, only: gmapex_p, gmapex_s
    use edyn3d_params, only: nhgt_fix, nmlat_T1

    type(horiz_coord_t), pointer :: flp_coord
    type(horiz_coord_t), pointer :: fls1_coord
    type(horiz_coord_t), pointer :: fls2_coord
    type(horiz_coord_t), pointer :: lonp_coord
    type(horiz_coord_t), pointer :: lons1_coord
    type(horiz_coord_t), pointer :: lons2_coord
    integer(iMap),       pointer :: grid_map(:,:)
    integer(iMap),       pointer :: coord_map(:)
    integer                      :: h, i, j, ind

    real(r8), pointer :: latvals(:)
    real(r8), pointer :: altvals(:)
    real(r8), pointer :: latvalss2(:)
    real(r8), pointer :: altvalss2(:)

    integer, parameter :: magp_decomp = 701 ! Must be unique within CAM
    integer, parameter :: mags1_decomp = 702 ! Must be unique within CAM
    integer, parameter :: mags2_decomp = 703 ! Must be unique within CAM

    real(r8) :: xdel

    integer :: isn, jj, k, ncnt
    integer :: dk,k0,k1

    if (mytid>=ntask) then
       if (mlon0_p/=1) then
          call endrun('register_grids: mlat0_p needs to be 1 on inactive PEs')
       end if
       if (mlon1_p/=0) then
          call endrun('register_grids: mlat1_p needs to be 0 on inactive PEs')
       end if
    end if

    nullify(flp_coord)
    nullify(fls1_coord)
    nullify(fls2_coord)
    nullify(lonp_coord)
    nullify(lons1_coord)
    nullify(lons2_coord)
    nullify(grid_map)
    nullify(coord_map)
    nullify(latvals)
    nullify(altvals)
    nullify(latvalss2)
    nullify(altvalss2)

 ! p and S1 -grids

    ncnt = 0
    do j = 1,nmlat_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = gmapex_p(j,isn)%npts
             dk = 1
          else
             k0 = gmapex_p(j,isn)%npts
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
          end do
       end do
    end do

    nptsp_total = ncnt
    allocate(latvals(nptsp_total))
    allocate(altvals(nptsp_total))
    latvals = huge(1._r8)
    altvals = -huge(1._r8)

    nptss1_total = nptsp_total  ! s1 grid points same as p grid points

    ncnt = 0
    do j = 1,nmlat_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = gmapex_p(j,isn)%npts
             dk = 1
          else
             k0 = gmapex_p(j,isn)%npts
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
             latvals(ncnt) = gmapex_p(j,isn)%mlat_qd(k)*r2d ! degrees
             altvals(ncnt) = gmapex_p(j,isn)%hgt_pt(k)*1.e-3_r8 ! km
          end do
       end do
    end do

    flp_coord => horiz_coord_create('lat_p', 'pflpt', nptsp_total, 'magnetic latitude', &
                                    'degrees_north', 1, nptsp_total, latvals)
    fls1_coord => horiz_coord_create('lat_s1', 'pflpt', nptsp_total, 'magnetic latitude', &
                                    'degrees_north', 1, nptsp_total, latvals)


    allocate(coord_map(mlon1_p-mlon0_p+1))
    coord_map = (/ (i, i = mlon0_p,mlon1_p ) /)

    lonp_coord => horiz_coord_create('lon_p', '', nmlon, 'magnetic longitude', &
                                    'degrees_east', mlon0_p,mlon1_p, r2d*ylonm(mlon0_p:mlon1_p), map=coord_map)
    lons2_coord => horiz_coord_create('lon_s2', '', nmlon, 'magnetic longitude', &
                                    'degrees_east', mlon0_p,mlon1_p, r2d*ylonm(mlon0_p:mlon1_p), map=coord_map)
    lons1_coord => horiz_coord_create('lon_s1', '', nmlon, 'magnetic longitude', &
                                    'degrees_east', mlon0_p,mlon1_p, r2d*ylonm_s(mlon0_p:mlon1_p), map=coord_map)

    nullify(coord_map)

    allocate(grid_map(4, ((mlon1_p-mlon0_p+1) * nptsp_total)))
    grid_map = -huge(1_iMap)
    ind = 0
    do i = 1,nptsp_total
       do j = mlon0_p,mlon1_p
          ind = ind + 1
          grid_map(1, ind) = j
          grid_map(2, ind) = i
          grid_map(3, ind) = j
          grid_map(4, ind) = i
       end do
    end do


    call cam_grid_register('magfline_p', magp_decomp, flp_coord, lonp_coord, grid_map, unstruct=.false.)
    call cam_grid_register('magfline_s1', mags1_decomp, fls1_coord, lons1_coord, grid_map, unstruct=.false.)
    nullify(flp_coord)
    nullify(fls1_coord)
    nullify(lons1_coord)
    nullify(grid_map)


 ! S2-grid


    ncnt = 0
    do j = 1,nmlats2_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = gmapex_s(j,isn)%npts
             dk = 1
          else
             k0 = gmapex_s(j,isn)%npts
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
          end do
       end do
    end do

    nptss2_total = ncnt
    allocate(latvalss2(nptss2_total))
    allocate(altvalss2(nptss2_total))
    latvalss2 = huge(1._r8)
    altvalss2 = -huge(1._r8)

    ncnt = 0
    do j = 1,nmlats2_h
       do isn = 1,2

          if (isn==1) then
             k0 = 1
             k1 = gmapex_s(j,isn)%npts
             dk = 1
          else
             k0 = gmapex_s(j,isn)%npts
             k1 = 1
             dk = -1
          endif

          do k = k0,k1,dk
             ncnt = ncnt + 1
             latvalss2(ncnt) = gmapex_s(j,isn)%mlat_qd(k)*r2d ! degrees
             altvalss2(ncnt) = gmapex_s(j,isn)%hgt_pt(k)*1.e-3_r8 ! km
          end do
       end do
    end do

    fls2_coord => horiz_coord_create('lat_s2', 's2flpt', nptss2_total, 'magnetic latitude', &
                                    'degrees_north', 1, nptss2_total, latvalss2)

    allocate(grid_map(4, ((mlon1_p-mlon0_p+1) * nptss2_total)))
    grid_map = -huge(1_iMap)
    ind = 0
    do i = 1,nptss2_total
       do j = mlon0_p,mlon1_p
          ind = ind + 1
          grid_map(1, ind) = j
          grid_map(2, ind) = i
          grid_map(3, ind) = j
          grid_map(4, ind) = i
       end do
    end do

    call cam_grid_register('magfline_s2', mags2_decomp, fls2_coord, lons2_coord, grid_map, unstruct=.false.)
    nullify(fls2_coord)
    nullify(lonp_coord)
    nullify(lons1_coord)
    nullify(lons2_coord)
    nullify(grid_map)

    call cam_grid_attribute_register('magfline_p', 'alt_p', 'magnetic field line p-grid altitude (km)', &
                                     'pflpt', altvals)
    call cam_grid_attribute_register('magfline_s1', 'alt_s1', 'magnetic field line s1-grid altitude (km)', &
                                     'pflpt', altvals)

    call cam_grid_attribute_register('magfline_s2', 'alt_s2', 'magnetic field line s1-grid altitude (km)', &
                                     's2flpt', altvalss2)

    nullify(latvals)
    nullify(altvals)

    nullify(latvalss2)
    nullify(altvalss2)
  end subroutine reg_hist_grid

!--------------------------------------------------------------------------------
  subroutine edyn3D_coef_halos
!
!    Calculate halo points in require variables
!
       use edyn3D_fieldline,only: fline_s1
       use edyn3D_mpi,      only: mytid,ntask,mp_mag_halos_edyn3D,mlon0_p,mlon1_p
       use edyn3D_params,   only: nptss1_max
!
!    Local:
!    5 fields to get halo points:
!
       integer, parameter :: nf = 5
       real(r8), allocatable :: fmsub(:,:,:,:)
!            real(r8) :: fmsub(mlon0_p-1:mlon1_p+1,nmlat_h,nptss1_max,nf)
!       real(r8) :: fmsub(mlon0_p-1:mlon1_p+1,nptss1_max,nf)
       real(r8) :: tempfield
       integer :: isn,i,j,k,status,nlon_p

       nlon_p = mlon1_p-mlon0_p+3

       allocate(fmsub(mlon0_p-1:mlon1_p+1,nmlat_h,nptss1_max,nf),STAT=status)
       if (status /= 0 ) then
         write(iulog,*) 'fmsub allocation failed'
	 call endrun('edyn3D_coef_halos')
       endif

       do isn = 1,2
	 !
	 ! Reset input to halos routine for each hemisphere
	 !
         fmsub(:,:,:,:) = 0._r8

         do j = 1, nmlat_h

           !
	   ! Reset input to halos routine for each latitude
	   !
!           fmsub(:,:,:) = 0._r8

           do k = 1,fline_s1(mlon0_p,j,isn)%npts ! Every longitude has same number of field line points

             do i = mlon0_p,mlon1_p

               fmsub(i,j,k,1) = fline_s1(i,j,isn)%M1(k)
               fmsub(i,j,k,2) = fline_s1(i,j,isn)%Je1D(k)
               fmsub(i,j,k,3) = fline_s1(i,j,isn)%M1Je1D(k)
               fmsub(i,j,k,4) = fline_s1(i,j,isn)%N1h(k)
               fmsub(i,j,k,5) = fline_s1(i,j,isn)%N1p(k)

             enddo ! longitude

           enddo ! field line points

         enddo ! latitude

	   call mp_mag_halos_edyn3D(fmsub,mlon0_p,mlon1_p,nmlat_h,nptss1_max,nf)

          do j = 1, nmlat_h
	    do k = 1,fline_s1(mlon0_p,j,isn)%npts ! Field lines at every longitude have same number of points

		fline_s1(mlon0_p-1,j,isn)%M1(k) = fmsub(mlon0_p-1,j,k,1)
		fline_s1(mlon1_p+1,j,isn)%M1(k) = fmsub(mlon1_p+1,j,k,1)
                fline_s1(mlon0_p-1,j,isn)%Je1D(k) = fmsub(mlon0_p-1,j,k,2)
                fline_s1(mlon1_p+1,j,isn)%Je1D(k) = fmsub(mlon1_p+1,j,k,2)
                fline_s1(mlon0_p-1,j,isn)%M1Je1D(k) = fmsub(mlon0_p-1,j,k,3)
                fline_s1(mlon1_p+1,j,isn)%M1Je1D(k) = fmsub(mlon1_p+1,j,k,3)
                fline_s1(mlon0_p-1,j,isn)%N1h(k) = fmsub(mlon0_p-1,j,k,4)
                fline_s1(mlon1_p+1,j,isn)%N1h(k) = fmsub(mlon1_p+1,j,k,4)
                fline_s1(mlon0_p-1,j,isn)%N1p(k) = fmsub(mlon0_p-1,j,k,5)
                fline_s1(mlon1_p+1,j,isn)%N1p(k) = fmsub(mlon1_p+1,j,k,5)

           enddo ! field line points
	 enddo ! latitude
       enddo ! N/S hemisphere

       deallocate(fmsub,STAT=status)
       if(status /= 0) then
	  write(iulog,*) 'deallocation of fmsub not successful'
	  call endrun('edyn3D_coef_halos')
       endif

  end subroutine edyn3D_coef_halos

!--------------------------------------------------------------------------------
  subroutine edyn3D_poten_halos
    !
    !    Calculate potential halo points for efield and ion drift velocity calculations
    !
    use edyn3D_fieldline,only: fline_p
    use edyn3D_mpi,	 only: mytid,ntask,mp_poten_halos_edyn3D,mlon0_p,mlon1_p
!    use edyn3D_params,   only: nptss1_max
!
! Local:
! Potential field to get halo points:
!
    real(r8), allocatable :: potensub(:,:)
!	  real(r8) :: fmsub(mlon0_p-1:mlon1_p+1,nmlat_h,nptss1_max,nf)
!    real(r8) :: fmsub(mlon0_p-1:mlon1_p+1,nptss1_max,nf)
    integer :: isn,i,j,k,status,nlon_p

    nlon_p = mlon1_p-mlon0_p+3

    allocate(potensub(mlon0_p-1:mlon1_p+1,nmlat_h),STAT=status)
    if (status /= 0 ) then
      write(iulog,*) 'potensub allocation failed'
      call endrun('edyn3D_poten_halos')
    endif

    do isn = 1,2
      !
      ! Reset input to halos routine for each hemisphere
      !
      potensub(:,:) = 0._r8

      do j = 1, nmlat_h
 	do i = mlon0_p,mlon1_p

 	  potensub(i,j) = fline_p(i,j,isn)%pot

 	enddo ! longitude
      enddo ! latitude

      call mp_poten_halos_edyn3D(potensub,mlon0_p,mlon1_p,nmlat_h)

      do j = 1, nmlat_h

     	fline_p(mlon0_p-1,j,isn)%pot = potensub(mlon0_p-1,j)
     	fline_p(mlon1_p+1,j,isn)%pot = potensub(mlon1_p+1,j)

      enddo ! latitude
    enddo ! N/S hemisphere

    deallocate(potensub,STAT=status)
    if (status /= 0) then
       write(iulog,*) 'deallocation of potensub not successful'
       call endrun('edyn3D_poten_halos')
    endif

  end subroutine edyn3D_poten_halos

!--------------------------------------------------------------------------------
  subroutine output_fline_field( magfld )
    use cam_history,  only: outfld
    use edyn3D_fline_fields, only: magfield_t

    type(magfield_t), intent(in) :: magfld

    integer :: i,j,k,isn,ncnt
    integer :: dk,k0,k1
    real(r8) :: fld_tmp(magfld%mlon0:magfld%mlon1,magfld%nptstot)

    if (mytid<ntask) then

       do i = magfld%mlon0,magfld%mlon1
          ncnt = 0
          do j = 1,magfld%nmlat_h
             do isn = 1,2

                if (isn==1) then
                   k0 = 1
                   k1 = magfld%flines(i,j,isn)%npts
                   dk = 1
                else
                   k0 = magfld%flines(i,j,isn)%npts
                   k1 = 1
                   dk = -1
                endif

                do k = k0,k1,dk
                   ncnt = ncnt + 1
                   fld_tmp(i,ncnt) = magfld%flines(i,j,isn)%fld(k)
                end do
             end do
          end do
       end do

       do j = 1,magfld%nptstot
          call outfld(magfld%name, fld_tmp(magfld%mlon0:magfld%mlon1,j), magfld%mlon1-magfld%mlon0+1, j)
       end do
    end if

  end subroutine output_fline_field

end module edyn3D_driver
