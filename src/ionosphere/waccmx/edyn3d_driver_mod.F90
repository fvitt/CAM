module edyn3d_driver_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc, mpicom
  use mpi_module, only: mpi_size, mpi_rank
  use infnan, only: nan, assignment(=)
  use perf_mod, only: t_startf, t_stopf

  implicit none

  private
  public :: edyn3d_driver_init
  public :: edyn3d_driver_timestep
  public :: edyn3d_driver_final

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_driver_init( mpicom_atm, npes_edyn3D, edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt, hilat_pot_model, wei05_coefs_file )
    use mpi_module, only: mpi_rank, mpi_size, lat_size, lon_size
    use params_module,only: nmlat_h, nmlon, nhgt_fix, hgt_fix_r

    use edyn3d_esmf_fields_rhandles, only: edyn3d_esmf_fields_rhandles_init
    use edyn3d_esmf_phys_mesh_mod, only:  edyn3d_esmf_phys_mesh_init
    use edyn3d_esmf_oplus_grid_mod, only:  edyn3d_esmf_oplus_grid_init
    use edyn3d_esmf_s1_mag_grid_mod, only: edyn3d_esmf_s1_mag_grid_init
    use edyn3d_esmf_s2_mag_grid_mod, only: edyn3d_esmf_s2_mag_grid_init
    use edyn3d_esmf_mag_ref_p_grid_mod, only: edyn3d_esmf_mag_ref_p_grid_init

    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_grids_reg
    use cam_history, only: addfld, horiz_only
    use edyn3d_highlat_potential, only: edyn3d_highlat_potential_init

    use mo_apex, only: mo_apex_init1
    use dynamo_interface_mod, only: dynamo_init1, dynamo_init2

    integer, intent(in) :: mpicom_atm, npes_edyn3D
    integer, intent(in) :: edyn3d_nmlat_h, edyn3d_nmlon, edyn3d_nhgt
    character(len=*),intent(in) :: hilat_pot_model
    character(len=*),intent(in) :: wei05_coefs_file

    character(len=*), parameter :: prefix = 'edyn3d_driver_init: '

    ! before apex init
    call dynamo_init1( set_hilat_pot_in=.true., set_hilat_fac_in=.false., &
         mpicom_atm=mpicom_atm, npes_edyn3D=npes_edyn3D, &
         edyn3d_nmlat_h=edyn3d_nmlat_h, edyn3d_nmlon=edyn3d_nmlon, edyn3d_nhgt=edyn3d_nhgt, &
         real_kind=r8 )

    ! log grid info:
    if (masterproc) then
       write(iulog,*) prefix,'3D Edyn grid params nmlat_h,nmlon,nhgt_fix: ',nmlat_h,nmlon,nhgt_fix
       write(iulog,*) prefix,'3D Edyn mpi_rank, mpi_size: ',mpi_rank,mpi_size
       write(iulog,*) prefix,'3D Edyn lon_size, lat_size: ',lon_size,lat_size
    end if

    ! initialize APEX -- after dynamo_init1 and before dynamo_init2
    call mo_apex_init1( alts_in=hgt_fix_r*1.e-3_r8 ) ! m --> km

    ! after apex init
    call dynamo_init2()

    ! setup cam history grids
    call edyn3d_hist_mag_grids_reg()

    ! setup ESMF grids
    call edyn3d_esmf_s2_mag_grid_init()
    call edyn3d_esmf_s1_mag_grid_init()
    call edyn3d_esmf_mag_ref_p_grid_init()

    call edyn3d_esmf_phys_mesh_init()
    call edyn3d_esmf_oplus_grid_init()

    ! setup ESMF regrid weights
    call edyn3d_esmf_fields_rhandles_init()

    ! initialize high-latitude electric potential inputs
    call edyn3d_highlat_potential_init(hilat_pot_model,wei05_coefs_file)

    ! cam history fields
    call addfld ('sigma_ped_s1', horiz_only, 'I', 'K','Ped cond. on S1 mag field line grid', gridname='magfline_s1')
    call addfld ('sigma_hal_s1', horiz_only, 'I', 'K','Hal cond. on S1 mag field line grid', gridname='magfline_s1')
    call addfld ('sigma_ped_s2', horiz_only, 'I', 'K','Ped cond. on S2 mag field line grid', gridname='magfline_s2')
    call addfld ('sigma_hal_s2', horiz_only, 'I', 'K','Hal cond. on S2 mag field line grid', gridname='magfline_s2')

    call addfld ('un_s1', horiz_only, 'I', 'm/s','Zonal wind on S1 mag field line grid', gridname='magfline_s1')
    call addfld ('vn_s1', horiz_only, 'I', 'm/s','Meridional wind on S1 mag field line grid', gridname='magfline_s1')
    call addfld ('un_s2', horiz_only, 'I', 'm/s','Zonal wind on S2 mag field line grid', gridname='magfline_s2')
    call addfld ('vn_s2', horiz_only, 'I', 'm/s','Meridional wind on S2 mag field line grid', gridname='magfline_s2')

    call addfld ('IonU_s1', horiz_only, 'I', 'm/s','Zonal Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('IonV_s1', horiz_only, 'I', 'm/s','Meridional Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('IonW_s1', horiz_only, 'I', 'm/s','Verical Ion Drift Velocity on s1 grid', gridname='magfline_s1')
    call addfld ('IonU_s2', horiz_only, 'I', 'm/s','Zonal Ion Drift Velocity on s1 grid', gridname='magfline_s2')
    call addfld ('IonV_s2', horiz_only, 'I', 'm/s','Meridional Ion Drift Velocity on s1 grid', gridname='magfline_s2')
    call addfld ('IonW_s2', horiz_only, 'I', 'm/s','Verical Ion Drift Velocity on s1 grid', gridname='magfline_s2')

    call addfld ('IonU_opg', (/ 'lev' /), 'I', 'm/s','Zonal Ion Drift Velocity on oplus grid' , gridname='geo_grid')
    call addfld ('IonV_opg', (/ 'lev' /), 'I', 'm/s','Meridional Ion Drift Velocity on oplus grid' , gridname='geo_grid')
    call addfld ('IonW_opg', (/ 'lev' /), 'I', 'm/s','Vertical Ion Drift Velocity on oplus grid' , gridname='geo_grid')

    call addfld ('ELECPOTEN', horiz_only, 'I', 'Volts','Electric potential', gridname='geomag_p')
    call addfld ('HILAT_POT', horiz_only, 'I', 'Volts','High-Latitude potential', gridname='geomag_p')
    call addfld ('HILAT_FAC', horiz_only, 'I', '???','High-Latitude field-aligned current', gridname='geomag_p')

    call addfld ('POTEN_opg', horiz_only, 'I', 'Volts', 'Electric potential', gridname='geo_grid')
    call addfld ('HLPOT_opg', horiz_only, 'I', 'Volts', 'High-latitude potential', gridname='geo_grid')
    call addfld ('HLFAC_opg', horiz_only, 'I', '???', 'High-Latitude field-aligned current', gridname='geo_grid')

    call addfld ('ED1s1', horiz_only, 'I', 'V/m','Electric field component', gridname='geomag_s1')
    call addfld ('ED2s1', horiz_only, 'I', 'V/m','Electric field component', gridname='geomag_s1')
    call addfld ('ED1s2', horiz_only, 'I', 'V/m','Electric field component', gridname='geomag_s2')
    call addfld ('ED2s2', horiz_only, 'I', 'V/m','Electric field component', gridname='geomag_s2')

    call addfld ('Ve1s1', horiz_only, 'I', 'm/s','Ion Drift Velocity', gridname='geomag_s1')
    call addfld ('Ve2s1', horiz_only, 'I', 'm/s','Ion Drift Velocity', gridname='geomag_s1')
    call addfld ('Ve1s2', horiz_only, 'I', 'm/s','Ion Drift Velocity', gridname='geomag_s2')
    call addfld ('Ve2s2', horiz_only, 'I', 'm/s','Ion Drift Velocity', gridname='geomag_s2')

  end subroutine edyn3d_driver_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_driver_timestep( nphyscol, nphyslev, physalt, sigPed, sigHal, un, vn, ui_oplus, vi_oplus, wi_oplus )
    use edyn3d_remap_mod, only: edyn3d_remap_phys2mag_s1
    use edyn3d_remap_mod, only: edyn3d_remap_phys2mag_s2
    use edyn3d_remap_mod, only: edyn3d_remap_refp_mag2oplus
    use edyn3d_remap_mod, only: edyn3d_remap_mag2oplus, NOTSET
    use edyn3d_remap_mod, only: mag_fields_bundle_t, phys_fields_bundle_t, oplus_fields_bundle_t
    use edyn3d_remap_mod, only: mag_2d_fields_bundle_t, oplus_2d_fields_bundle_t
    use edyn3d_esmf_fields_rhandles, only: phys2mag_nflds
    use edyn3D_esmf_fields_rhandles, only: mag2opls_nflds
    use mpi_module, only: mlat0, mlat1, mlon0, mlon1
    use mpi_module, only: mlond0, mlond1, mlatd0, mlatd1
    use params_module,only:  nhgt_fix
    use edyn_mpi, only: lon0,lon1,lat0,lat1,lev0,lev1
    use regridder, only: regrid_phys2geo_3d, regrid_geo2phys_3d
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_s1_out
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_s2_out
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mlonlat_out
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mlonlat_s_out

    use dynamo_interface_mod, only: dynamo_calc
    use edyn3d_highlat_potential, only: edyn3d_highlat_potential_get

    use cam_history,  only: outfld

    integer,  intent(in) :: nphyscol, nphyslev
    real(r8), intent(in) :: physalt(nphyslev,nphyscol)

    real(r8), target, intent(in) :: sigPed(nphyslev,nphyscol)
    real(r8), target, intent(in) :: sigHal(nphyslev,nphyscol)
    real(r8), target, intent(in) :: un(nphyslev,nphyscol)
    real(r8), target, intent(in) :: vn(nphyslev,nphyscol)
    real(r8), target, intent(out) :: ui_oplus(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8), target, intent(out) :: vi_oplus(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8), target, intent(out) :: wi_oplus(lon0:lon1,lat0:lat1,lev0:lev1)

    real(r8) :: opalt(lon0:lon1,lat0:lat1,lev0:lev1)

    type(phys_fields_bundle_t) :: phys_flds_bndl(phys2mag_nflds)
    type(mag_fields_bundle_t) :: mags1_flds_bndl(phys2mag_nflds)
    type(mag_fields_bundle_t) :: mags2_flds_bndl(phys2mag_nflds)
    type(mag_fields_bundle_t) :: magsrc_flds_bndl(mag2opls_nflds)
    type(oplus_fields_bundle_t) :: oplus_flds_bndl(mag2opls_nflds)

    type(mag_2d_fields_bundle_t) :: magsrc_2d_flds_bndl(mag2opls_nflds)
    type(oplus_2d_fields_bundle_t) :: oplus_2d_flds_bndl(mag2opls_nflds)

    real(r8), target :: sigped_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: sighal_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: un_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: vn_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8), target :: sigped_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: sighal_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: un_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: vn_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8), target :: pot_hl_p(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8), target :: fac_hl_p(2,mlatd0:mlatd1,mlond0:mlond1)
    real(r8), target :: hlfac_op(lon0:lon1,lat0:lat1)
    real(r8), target :: hlpot_op(lon0:lon1,lat0:lat1)
    real(r8), target :: poten_op(lon0:lon1,lat0:lat1)

    real(r8) :: ui_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8) :: vi_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8) :: wi_s1(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8), target :: ui_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: vi_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)
    real(r8), target :: wi_s2(nhgt_fix,2,mlat0:mlat1,mlon0:mlon1)

    real(r8), target :: elec_pot_p(2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: efld1_s1(2,mlat0:mlat1,mlon0:mlon1)
    real(r8) :: efld2_s1(2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: efld1_s2(2,mlat0:mlat1,mlon0:mlon1)
    real(r8) :: efld2_s2(2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: ionvel1_s1(2,mlat0:mlat1,mlon0:mlon1)
    real(r8) :: ionvel2_s1(2,mlat0:mlat1,mlon0:mlon1)

    real(r8) :: ionvel1_s2(2,mlat0:mlat1,mlon0:mlon1)
    real(r8) :: ionvel2_s2(2,mlat0:mlat1,mlon0:mlon1)

    integer :: j

    character(len=*), parameter :: subname = 'edyn3d_driver_timestep'

    call t_startf(subname)

    sigped_s1 = NOTSET
    sighal_s1 = NOTSET
    un_s1 = NOTSET
    vn_s1 = NOTSET

    sigped_s2 = NOTSET
    sighal_s2 = NOTSET
    un_s2 = NOTSET
    vn_s2 = NOTSET

    mags1_flds_bndl(1)%fld => sigped_s1
    mags1_flds_bndl(2)%fld => sighal_s1
    mags1_flds_bndl(3)%fld => un_s1
    mags1_flds_bndl(4)%fld => vn_s1

    mags2_flds_bndl(1)%fld => sigped_s2
    mags2_flds_bndl(2)%fld => sighal_s2
    mags2_flds_bndl(3)%fld => un_s2
    mags2_flds_bndl(4)%fld => vn_s2

    phys_flds_bndl(1)%fld => sigPed
    phys_flds_bndl(2)%fld => sigHal
    phys_flds_bndl(3)%fld => un
    phys_flds_bndl(4)%fld => vn

    call t_startf(subname//'->remap_phys2mag_s1')
    call edyn3d_remap_phys2mag_s1(nphyscol, nphyslev, physalt, phys_flds_bndl, mags1_flds_bndl)
    call t_stopf(subname//'->remap_phys2mag_s1')

    call edyn3d_hist_mag_s1_out('sigma_ped_s1',sigped_s1)
    call edyn3d_hist_mag_s1_out('sigma_hal_s1',sighal_s1)
    call edyn3d_hist_mag_s1_out('un_s1',un_s1)
    call edyn3d_hist_mag_s1_out('vn_s1',vn_s1)

    call t_startf(subname//'->remap_phys2mag_s2')
    call edyn3d_remap_phys2mag_s2(nphyscol, nphyslev, physalt, phys_flds_bndl, mags2_flds_bndl)
    call t_stopf(subname//'->remap_phys2mag_s2')

    call edyn3d_hist_mag_s2_out('sigma_ped_s2',sigped_s2)
    call edyn3d_hist_mag_s2_out('sigma_hal_s2',sighal_s2)
    call edyn3d_hist_mag_s2_out('un_s2',un_s2)
    call edyn3d_hist_mag_s2_out('vn_s2',vn_s2)

    if (mpi_rank<mpi_size) then

       call edyn3d_highlat_potential_get(pot_hl_p)

       call dynamo_calc( &
            sigped_s1, sighal_s1, un_s1, vn_s1, &
            sigped_s2, sighal_s2, un_s2, vn_s2, &
            pot_hl_p, fac_hl_p, &
            ui_s1, vi_s1, wi_s1, ui_s2, vi_s2, wi_s2, &
            elec_pot_p, &
            efld1_s1, efld2_s1, efld1_s2, efld2_s2, &
            ionvel1_s1, ionvel2_s1, ionvel1_s2, ionvel2_s2)

       call edyn3d_hist_mlonlat_out('HILAT_POT',   pot_hl_p(1:2,mlat0:mlat1,mlon0:mlon1))
       call edyn3d_hist_mlonlat_out('HILAT_FAC',   fac_hl_p(1:2,mlat0:mlat1,mlon0:mlon1))
       call edyn3d_hist_mlonlat_out('ELECPOTEN', elec_pot_p(1:2,mlat0:mlat1,mlon0:mlon1))

       call edyn3d_hist_mlonlat_out('ED1s1', efld1_s1)
       call edyn3d_hist_mlonlat_out('ED2s1', efld2_s1)
       call edyn3d_hist_mlonlat_s_out('ED1s2', efld1_s2)
       call edyn3d_hist_mlonlat_s_out('ED2s2', efld2_s2)

       call edyn3d_hist_mlonlat_out('Ve1s1',  ionvel1_s1)
       call edyn3d_hist_mlonlat_out('Ve2s1',  ionvel2_s1)
       call edyn3d_hist_mlonlat_s_out('Ve1s2', ionvel1_s2)
       call edyn3d_hist_mlonlat_s_out('Ve2s2', ionvel2_s2)

       call edyn3d_hist_mag_s1_out('IonU_s1', ui_s1)
       call edyn3d_hist_mag_s1_out('IonV_s1', vi_s1)
       call edyn3d_hist_mag_s1_out('IonW_s1', wi_s1)

       call edyn3d_hist_mag_s2_out('IonU_s2', ui_s2)
       call edyn3d_hist_mag_s2_out('IonV_s2', vi_s2)
       call edyn3d_hist_mag_s2_out('IonW_s2', wi_s2)

    end if

    ! map to oplus geographic grid and output diagnostics
    magsrc_2d_flds_bndl(1)%fld => fac_hl_p
    magsrc_2d_flds_bndl(2)%fld => pot_hl_p
    magsrc_2d_flds_bndl(3)%fld => elec_pot_p
    oplus_2d_flds_bndl(1)%fld => hlfac_op
    oplus_2d_flds_bndl(2)%fld => hlpot_op
    oplus_2d_flds_bndl(3)%fld => poten_op

    call edyn3d_remap_refp_mag2oplus( magsrc_2d_flds_bndl, oplus_2d_flds_bndl )

    do j = lat0,lat1
       call outfld( 'HLFAC_opg', hlfac_op(lon0:lon1,j), lon1-lon0+1, j )
       call outfld( 'HLPOT_opg', hlpot_op(lon0:lon1,j), lon1-lon0+1, j )
       call outfld( 'POTEN_opg', poten_op(lon0:lon1,j), lon1-lon0+1, j )
    end do

    ! map ion vels to oplus xport grid (geographic)
    magsrc_flds_bndl(1)%fld => ui_s2
    magsrc_flds_bndl(2)%fld => vi_s2
    magsrc_flds_bndl(3)%fld => wi_s2

    oplus_flds_bndl(1)%fld => ui_oplus
    oplus_flds_bndl(2)%fld => vi_oplus
    oplus_flds_bndl(3)%fld => wi_oplus

    call t_startf(subname//'->regrid_phys2geo_3d')
    call regrid_phys2geo_3d( physalt, opalt, nphyslev, 1, nphyscol )
    call t_stopf(subname//'->regrid_phys2geo_3d')

    call t_startf(subname//'->remap_mag2oplus')
    ! invert the oplus altitudes (need bottom-up order)
    call edyn3d_remap_mag2oplus( magsrc_flds_bndl, opalt(:,:,lev1:lev0:-1), oplus_flds_bndl )
    call t_stopf(subname//'->remap_mag2oplus')

    do j = lat0,lat1
       call outfld( 'IonU_opg', ui_oplus(lon0:lon1,j,lev1:lev0:-1), lon1-lon0+1, j )
       call outfld( 'IonV_opg', vi_oplus(lon0:lon1,j,lev1:lev0:-1), lon1-lon0+1, j )
       call outfld( 'IonW_opg', wi_oplus(lon0:lon1,j,lev1:lev0:-1), lon1-lon0+1, j )
    end do

    call t_stopf(subname)

  end subroutine edyn3d_driver_timestep

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_driver_final()
    use edyn3d_esmf_fields_rhandles, only: edyn3d_esmf_fields_rhandles_destroy
    use edyn3d_esmf_oplus_grid_mod, only: edyn3d_esmf_oplus_grid_destroy
    use edyn3d_esmf_phys_mesh_mod, only: edyn3d_esmf_phys_mesh_destroy
    use edyn3d_esmf_s1_mag_grid_mod, only: edyn3d_esmf_s1_mag_grid_destroy
    use edyn3d_esmf_s2_mag_grid_mod, only: edyn3d_esmf_s2_mag_grid_destroy
    use edyn3d_esmf_mag_ref_p_grid_mod, only: edyn3d_esmf_mag_ref_p_grid_destroy
    use edyn3d_hist_mag_grids_mod, only: edyn3d_hist_mag_grids_final

    call edyn3d_esmf_fields_rhandles_destroy()
    call edyn3d_esmf_oplus_grid_destroy()
    call edyn3d_esmf_phys_mesh_destroy()
    call edyn3d_esmf_s1_mag_grid_destroy()
    call edyn3d_esmf_s2_mag_grid_destroy()
    call edyn3d_esmf_mag_ref_p_grid_destroy()

    call edyn3d_hist_mag_grids_final()

  end subroutine edyn3d_driver_final

end module edyn3d_driver_mod
