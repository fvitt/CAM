module dpie_coupling
  !
  ! Dynamics/Physics Ionosphere/Electrodynamics coupler.
  !
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_logfile,         only: iulog
  use cam_history,         only: hist_fld_active, outfld
  use cam_history,         only: addfld, horiz_only
  use cam_history_support, only: fillvalue
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc, mpi_logical, mpicom, masterprocid
  use edyn_mpi,            only: array_ptr_type
  use perf_mod,            only: t_startf, t_stopf
  use amie_module,         only: getamie
  use ltr_module,          only: getltr
  use edyn_solve,          only: phihm
  use edyn_params,         only: dtr, rtd
  use aurora_params,       only: prescribed_period ! turns on overwrite of energy fields in aurora phys

  implicit none

  private
  public :: d_pie_init
  public :: d_pie_epotent  ! sets electric potential
  public :: d_pie_coupling ! handles coupling with edynamo and ion transport

  logical  :: ionos_edyn_active ! if true, call oplus_xport for O+ transport
  logical  :: ionos_oplus_xport ! if true, call oplus_xport for O+ transport
  integer  :: nspltop           ! nsplit for oplus_xport
  integer  :: nsplt_vxb=-1      ! nsplit for oplus_xport

  logical  :: debug = .false.

  real(r8) :: crad(2), crit(2)
  logical  :: crit_user_set = .false.
  real(r8), parameter :: amie_default_crit(2) = (/ 35._r8, 40._r8 /)

  logical :: debug_hist

contains
  !----------------------------------------------------------------------
  subroutine d_pie_init( edyn_active_in, oplus_xport_in, oplus_nsplit_in, vxb_nsplit_in, crit_colats_deg, ionos_debug_hist )

    logical, intent(in) :: edyn_active_in, oplus_xport_in
    integer, intent(in) :: oplus_nsplit_in
    integer, intent(in) :: vxb_nsplit_in
    real(r8),intent(in) :: crit_colats_deg(:)
    logical, intent(in) :: ionos_debug_hist

    debug_hist = ionos_debug_hist

    ionos_edyn_active = edyn_active_in
    ionos_oplus_xport = oplus_xport_in
    nspltop = oplus_nsplit_in
    nsplt_vxb = vxb_nsplit_in

    crit_user_set = all( crit_colats_deg(:) > 0._r8 )
    if (crit_user_set) then
       crit(:) = crit_colats_deg(:)*dtr
    end if

    call addfld ('HMF2'       , horiz_only , 'I', 'km'  ,'Height of the F2 Layer'      , gridname='physgrid')
    call addfld ('NMF2'       , horiz_only , 'I', 'cm-3','Peak Density of the F2 Layer', gridname='physgrid')

    call addfld ('OpDens'  ,(/ 'lev' /), 'I', 'cm^-3','O+  Number Density', gridname='physgrid')
    call addfld ('NOpDens' ,(/ 'lev' /), 'I', 'cm^-3','NO+ Number Density', gridname='physgrid')
    call addfld ('O2pDens' ,(/ 'lev' /), 'I', 'cm^-3','O2+ Number Density', gridname='physgrid')
    call addfld ('FepDens' ,(/ 'lev' /), 'I', 'cm^-3','Fe+ Number Density', gridname='physgrid')
    call addfld ('MgpDens' ,(/ 'lev' /), 'I', 'cm^-3','Mg+ Number Density', gridname='physgrid')
    call addfld ('NapDens' ,(/ 'lev' /), 'I', 'cm^-3','Na+ Number Density', gridname='physgrid')
    call addfld ('CapDens' ,(/ 'lev' /), 'I', 'cm^-3','Ca+ Number Density', gridname='physgrid')
    call addfld ('KpDens'  ,(/ 'lev' /), 'I', 'cm^-3','K+  Number Density', gridname='physgrid')
    call addfld ('SipDens' ,(/ 'lev' /), 'I', 'cm^-3','Si+ Number Density', gridname='physgrid')

    call addfld ('EDens'  ,(/ 'lev' /), 'I', 'cm^-3','e Number Density (sum of O2+,NO+,N2+,O+)', gridname='physgrid')

    call addfld ('prescr_efxp'  , horiz_only, 'I','mW/m2','Prescribed energy flux on geo grid'     ,gridname='physgrid')
    call addfld ('prescr_kevp'  , horiz_only, 'I','keV  ','Prescribed mean energy on geo grid'     ,gridname='physgrid')

    call addfld ('prescr_phihm' , horiz_only, 'I','VOLTS','Prescribed Electric Potential-mag grid' ,gridname='gmag_grid')
    call addfld ('prescr_efxm'  , horiz_only, 'I','mW/m2','Prescribed energy flux on mag grid'     ,gridname='gmag_grid')
    call addfld ('prescr_kevm'  , horiz_only, 'I','keV  ','Prescribed mean energy on mag grid'     ,gridname='gmag_grid')

    if (debug_hist) then
       ! Dynamo inputs (called from dpie_coupling. Fields are in waccm format, in CGS units):
       call addfld ('DPIE_OMEGA',(/ 'lev' /), 'I', 'Pa/s    ','OMEGA input to DPIE coupling', gridname='physgrid')
       call addfld ('DPIE_MBAR' ,(/ 'lev' /), 'I', 'kg/kmole','MBAR Mean Mass from dpie_coupling', gridname='physgrid')
       call addfld ('DPIE_TN   ',(/ 'lev' /), 'I', 'deg K   ','DPIE_TN'   , gridname='physgrid')
       call addfld ('DPIE_UN   ',(/ 'lev' /), 'I', 'cm/s    ','DPIE_UN'   , gridname='physgrid')
       call addfld ('DPIE_VN   ',(/ 'lev' /), 'I', 'cm/s    ','DPIE_VN'   , gridname='physgrid')
       call addfld ('DPIE_WN   ',(/ 'lev' /), 'I', 'cm/s    ','DPIE_WN'   , gridname='physgrid')
       call addfld ('DPIE_OM   ',(/ 'lev' /), 'I', 's-1     ','DPIE_OM'   , gridname='physgrid')
       call addfld ('DPIE_ZHT  ',(/ 'lev' /), 'I', 'cm      ','DPIE_ZHT (geometric height,simple)', gridname='physgrid')
       call addfld ('DPIE_ZGI  ',(/ 'lev' /), 'I', 'cm      ','DPIE_ZGI (geopotential height on interfaces)', gridname='physgrid')
       call addfld ('DPIE_O2   ',(/ 'lev' /), 'I', 'mmr     ','DPIE_O2'   , gridname='physgrid')
       call addfld ('DPIE_O    ',(/ 'lev' /), 'I', 'mmr     ','DPIE_O'    , gridname='physgrid')
       call addfld ('DPIE_N2   ',(/ 'lev' /), 'I', 'mmr     ','DPIE_N2'   , gridname='physgrid')
       call addfld ('DPIE_TE   ',(/ 'lev' /), 'I', 'deg K   ','DPIE_TE'   , gridname='physgrid')
       call addfld ('DPIE_TI   ',(/ 'lev' /), 'I', 'deg K   ','DPIE_TI'   , gridname='physgrid')

       call addfld ('PED_phys',(/ 'lev' /), 'I', 'S/m','Pedersen Conductivity'  , gridname='physgrid')
       call addfld ('HAL_phys',(/ 'lev' /), 'I', 'S/m','Hall Conductivity'   , gridname='physgrid')

       call addfld ('DPIE_OPMMR' ,(/ 'lev' /), 'I', 'mmr'  ,'DPIE_OPMMR'  , gridname='physgrid')
       call addfld ('DPIE_FePMMR' ,(/ 'lev' /), 'I', 'mmr'  ,'DPIE_FeMMR'  , gridname='physgrid')
       call addfld ('DPIE_MgPMMR' ,(/ 'lev' /), 'I', 'mmr'  ,'DPIE_MgMMR'  , gridname='physgrid')
       call addfld ('DPIE_O2P',(/ 'lev' /), 'I', 'm^-3','DPIE_O2P(dpie input)', gridname='physgrid')
       call addfld ('DPIE_NOP',(/ 'lev' /), 'I', 'm^-3','DPIE_NOP(dpie input)', gridname='physgrid')
       call addfld ('DPIE_N2P',(/ 'lev' /), 'I', 'm^-3','DPIE_N2P(dpie input)', gridname='physgrid')

       call addfld ('DPIE_FEP',(/ 'lev' /), 'I', 'm^-3','DPIE Fep(dpie input)', gridname='physgrid')
       call addfld ('DPIE_MGP',(/ 'lev' /), 'I', 'm^-3','DPIE Mgp(dpie input)', gridname='physgrid')
       call addfld ('DPIE_NAP',(/ 'lev' /), 'I', 'm^-3','DPIE Nap(dpie input)', gridname='physgrid')
       call addfld ('DPIE_CAP',(/ 'lev' /), 'I', 'm^-3','DPIE Cap(dpie input)', gridname='physgrid')
       call addfld ('DPIE_KP', (/ 'lev' /), 'I', 'm^-3','DPIE Kp (dpie input)', gridname='physgrid')
       call addfld ('DPIE_SIP',(/ 'lev' /), 'I', 'm^-3','DPIE Sip(dpie input)', gridname='physgrid')

       call addfld ('WACCM_UI'   ,(/ 'lev' /), 'I', 'm/s'  ,'WACCM_UI (dpie output)', gridname='physgrid')
       call addfld ('WACCM_VI'   ,(/ 'lev' /), 'I', 'm/s'  ,'WACCM_VI (dpie output)', gridname='physgrid')
       call addfld ('WACCM_WI'   ,(/ 'lev' /), 'I', 'm/s'  ,'WACCM_WI (dpie output)', gridname='physgrid')
       call addfld ('WACCM_OP'   ,(/ 'lev' /), 'I', 'kg/kg'  ,'WACCM_OP (dpie output)', gridname='physgrid')
       call addfld ('WACCM_FeP'  ,(/ 'lev' /), 'I', 'kg/kg'  ,'WACCM_FeP (dpie output)', gridname='physgrid')
       call addfld ('WACCM_MgP'  ,(/ 'lev' /), 'I', 'kg/kg'  ,'WACCM_MgP (dpie output)', gridname='physgrid')

       call addfld ('EDYN_ADOTV1 ', (/ 'lev' /), 'I', '        ','EDYN_ADOTV1' , gridname='geo_grid')
       call addfld ('EDYN_ADOTV2 ', (/ 'lev' /), 'I', '        ','EDYN_ADOTV2' , gridname='geo_grid')
       call addfld ('EDYN_ADOTA1 ', horiz_only , 'I', '        ','EDYN_ADOTA1' , gridname='geo_grid')
       call addfld ('EDYN_ADOTA2 ', horiz_only , 'I', '        ','EDYN_ADOTA2' , gridname='geo_grid')
       call addfld ('EDYN_A1DTA2 ', horiz_only , 'I', '        ','EDYN_A1DTA2' , gridname='geo_grid')

       call addfld ('EDYN_SINI   ', horiz_only , 'I', '        ','EDYN_SINI'   , gridname='geo_grid')
       call addfld ('EDYN_BE3    ', horiz_only , 'I', '        ','EDYN_BE3'    , gridname='geo_grid')

       call addfld ('ADOTA1_MAG', horiz_only , 'I', ' ','ADOTA1 in geo-mag coords' , gridname='gmag_grid')
       call addfld ('SINI_MAG',   horiz_only , 'I', ' ','sini in geo-mag coords' , gridname='gmag_grid')

       call addfld ('OPLUS', (/ 'lev' /), 'I', 'cm^3','O+ (oplus_xport output)',    gridname='geo_grid')
       call addfld ('OPtm1i',(/ 'lev' /), 'I', 'cm^3','O+ (oplus_xport output)',    gridname='geo_grid')
       call addfld ('OPtm1o',(/ 'lev' /), 'I', 'cm^3','O+ (oplus_xport output)',    gridname='geo_grid')

       call addfld ('FePtm1i',(/ 'lev' /), 'I', 'cm^3','Fe+ (oplus_xport output)',    gridname='geo_grid')
       call addfld ('FePtm1o',(/ 'lev' /), 'I', 'cm^3','Fe+ (oplus_xport output)',    gridname='geo_grid')

       call addfld ('MgPtm1i',(/ 'lev' /), 'I', 'cm^3','Mg+ (oplus_xport output)',    gridname='geo_grid')
       call addfld ('MgPtm1o',(/ 'lev' /), 'I', 'cm^3','Mg+ (oplus_xport output)',    gridname='geo_grid')

    endif


    call addfld ('OPi',(/ 'lev' /), 'I', 'cm^3','O+ (oplus_xport output)',    gridname='geo_grid')
    call addfld ('OPo',(/ 'lev' /), 'I', 'cm^3','O+ (oplus_xport output)',    gridname='geo_grid')

    call addfld ('FePi',(/ 'lev' /), 'I', 'cm^3','Fe+ (oplus_xport output)',    gridname='geo_grid')
    call addfld ('FePo',(/ 'lev' /), 'I', 'cm^3','Fe+ (oplus_xport output)',    gridname='geo_grid')

    call addfld ('MgPi',(/ 'lev' /), 'I', 'cm^3','Mg+ (oplus_xport output)',    gridname='geo_grid')
    call addfld ('MgPo',(/ 'lev' /), 'I', 'cm^3','Mg+ (oplus_xport output)',    gridname='geo_grid')


  end subroutine d_pie_init

  !-----------------------------------------------------------------------
  subroutine d_pie_epotent( highlat_potential_model, crit_out, cols, cole, efx_phys, kev_phys, amie_in, ltr_in )
    use edyn_solve,       only: pfrac    ! NH fraction of potential (nmlonp1,nmlat0)
    use time_manager,     only: get_curr_date
    use heelis,           only: heelis_model
    use wei05sc,          only: weimer05  ! driver for weimer high-lat convection model
    use edyn_esmf,        only: edyn_esmf_update
    use solar_parms_data, only: solar_parms_advance
    use solar_wind_data,  only: solar_wind_advance
    use solar_wind_data,  only: bzimf=>solar_wind_bzimf
    use solar_wind_data,  only: byimf=>solar_wind_byimf
    use solar_wind_data,  only: swvel=>solar_wind_swvel
    use solar_wind_data,  only: swden=>solar_wind_swden
    use edyn_mpi,         only: mlat0, mlat1, mlon0, mlon1, omlon1, ntask, mytid
    use edyn_maggrid,     only: nmlonp1, nmlat
    use regridder,  only: regrid_mag2phys_2d

    ! Args:
    !
    character(len=*),  intent(in)  :: highlat_potential_model
    real(r8),          intent(out) :: crit_out(2) ! critical colatitudes (degrees)
    integer, optional, intent(in)  :: cols, cole
    logical, optional,intent(in) :: amie_in
    logical, optional,intent(in) :: ltr_in

    ! Prescribed energy flux
    real(r8), optional, intent(out) :: efx_phys(:)
    ! Prescribed characteristic mean energy
    real(r8), optional, intent(out) :: kev_phys(:)

    !
    ! local vars
    !
    logical :: amie_inputs, ltr_inputs

    real(r8)             :: secs               ! time of day in seconds
    integer              :: iyear,imo,iday,tod ! tod is time-of-day in seconds
    real(r8)             :: sunlon

    integer              :: iprint
    integer              :: j, iamie, iltr, ierr
    !
    ! AMIE fields (extra dimension added for longitude switch)
    !
    real(r8) :: prescr_efxm(nmlonp1,nmlat), prescr_kevm(nmlonp1,nmlat)
    real(r8) :: prescr_phihm(nmlonp1,nmlat)

    call edyn_esmf_update()

    call get_curr_date(iyear, imo,iday, tod)
    ! tod is integer time-of-day in seconds
    secs = real(tod, r8)

    ! update solar wind data (IMF, etc.)
    call solar_wind_advance()

    ! update kp -- phys timestep init happens later ...
    call solar_parms_advance()
    if ( mytid<ntask ) then
       !
       ! Get sun's longitude at latitudes (geographic):
       !
       call sunloc(iday, secs, sunlon) ! sunlon is returned
       !
       ! Get high-latitude convection from empirical model (heelis or weimer).
       ! High-latitude potential phihm (edyn_solve) is defined for edynamo.
       !
       if (trim(highlat_potential_model) == 'heelis') then
          call heelis_model(sunlon) ! heelis.F90
       elseif (trim(highlat_potential_model) == 'weimer') then
          !
          call weimer05(byimf, bzimf, swvel, swden, sunlon)
          if (debug .and. masterproc) then
             write(iulog, "(a,2f8.2,a,2f8.2)")                                   &
                  'dpie_coupling call weimer05: byimf,bzimf=',                   &
                  byimf, bzimf, ' swvel,swden=', swvel, swden
          end if
       else
          call endrun('d_pie_epotent: Unknown highlat_potential_model')
       end if
    endif

    amie_inputs=.false.
    ltr_inputs=.false.
    if (present(amie_in)) amie_inputs=amie_in
    if (present(ltr_in))   ltr_inputs= ltr_in

    prescribed_inputs: if (amie_inputs .or. ltr_inputs) then

       if (.not. (present(kev_phys).and.present(efx_phys)) ) then
          call endrun('d_pie_epotent: kev_phys and efx_phys must be present')
       end if

       iprint = 1
       if (amie_inputs) then
          if (masterproc) then
             write(iulog,*) 'Calling getamie >>> '
          end if

          call getamie(iyear, imo, iday, tod, sunlon, iprint, iamie, &
               prescr_phihm, prescr_efxm, prescr_kevm, crad)

          if (masterproc) then
             write(iulog,"('After Calling getamie >>> iamie = ', i2)") iamie
          end if
          prescribed_period = iamie == 1
       else
          if (masterproc) then
             write(iulog,*) 'Calling getltr >>> '
          end if

          call getltr(iyear, imo, iday, tod,sunlon, iprint, iltr, &
               prescr_phihm, prescr_efxm, prescr_kevm )

          if (masterproc) then
             write(iulog,"('After Calling getltr >>> iltr = ', i2)") iltr
          end if
          prescribed_period = iltr == 1
       end if

       do j = mlat0, mlat1
          call outfld('prescr_phihm',prescr_phihm(mlon0:omlon1,j),omlon1-mlon0+1,j)
          call outfld('prescr_efxm', prescr_efxm(mlon0:omlon1,j), omlon1-mlon0+1,j)
          call outfld('prescr_kevm', prescr_kevm(mlon0:omlon1,j), omlon1-mlon0+1,j)
       end do

       if (prescribed_period) then
          phihm = prescr_phihm
       end if

       call mpi_bcast(prescribed_period, 1, mpi_logical, masterprocid, mpicom, ierr)

       call regrid_mag2phys_2d(prescr_kevm(mlon0:mlon1,mlat0:mlat1), kev_phys, cols, cole)
       call regrid_mag2phys_2d(prescr_efxm(mlon0:mlon1,mlat0:mlat1), efx_phys, cols, cole)

       call outfld_phys1d( 'prescr_efxp', efx_phys )
       call outfld_phys1d( 'prescr_kevp', kev_phys )

    end if prescribed_inputs

    if ( mytid<ntask ) then

       call calc_pfrac(sunlon, pfrac) ! returns pfrac for dynamo (edyn_solve)

       crit_out(:) = crit(:)*rtd ! degrees
    endif
  end subroutine d_pie_epotent

  !-----------------------------------------------------------------------
  subroutine d_pie_coupling(omega, pmid, zgi, zht, u, v, tn,                  &
       sigma_ped, sigma_hall, te, ti, mbar, n2mmr, o2mmr, o1mmr, o2pmmr,      &
       nopmmr, n2pmmr, opmmr, opmmrtm1, ui, vi, wi,                           &
       rmassN2, rmassO2, rmassO1, &
       rmassO2p, rmassNOp, rmassN2p, rmassOp, cols, cole, plev, &
       rmassFep, rmassMgp, rmassNap, rmassCap, rmassKp, rmassSip, &
       Fepmmr, Mgpmmr, Napmmr, Capmmr, Kpmmr, Sipmmr, &
       NOpmmrtm1, O2pmmrtm1, Fepmmrtm1, Mgpmmrtm1, Napmmrtm1, Capmmrtm1, Kpmmrtm1, Sipmmrtm1)
     !
     ! Call dynamo to calculate electric potential and electric field
     ! Note: dynamo calculates ion drifts.
     ! Then call oplus_xport to transport O+, which is passed back to physics.
     !
     use edyn_geogrid,  only: nlev
     use shr_const_mod, only: grav => shr_const_g ! gravitational accel. (m/s^2)
     use shr_const_mod, only: kboltz => shr_const_boltz ! Boltzmann const. (J/K/molecule)
     use time_manager,  only: get_nstep, get_curr_date
     use edynamo,       only: dynamo
     use edyn_mpi,      only: mp_geo_halos, mp_pole_halos
     ! Mag grid distribution info
     use edyn_mpi,      only: mlon0, mlon1, mlat0, mlat1, mlev0, mlev1
     use edyn_mpi,      only: lon0, lon1, lat0, lat1, lev0, lev1, ntask, mytid
     use oplus,         only: oplus_xport
     use ref_pres,      only: pref_mid
     use regridder,  only: regrid_phys2geo_3d, regrid_phys2mag_3d, regrid_geo2phys_3d
     use regridder,  only: regrid_geo2mag_3d, regrid_geo2mag_2d
     use adotv_mod,  only: calc_adotv
     use vxb, only: vxb_xport

     !
     ! Args:
     !
     integer,intent(in) :: &
          cols,            & ! First column (usually 1)
          cole,            & ! Last column
          plev               ! Number of levels in physics midpoint fields

     ! Input arguments on physics grid
     real(r8), intent(in)    :: omega(plev, cols:cole)      ! pressure velocity on midpoints (Pa/s) (i,k,j)
     real(r8), intent(in)    :: pmid(plev, cols:cole)       ! pressure at midpoints (Pa)
     real(r8), intent(in)    :: zgi(plev, cols:cole)        ! geopotential height (on interfaces) (m)
     real(r8), intent(in)    :: zht(plev, cols:cole)        ! geometric height (m) (Simple method - interfaces)
     real(r8), intent(in)    :: u(plev, cols:cole)          ! U-wind (m/s)
     real(r8), intent(in)    :: v(plev, cols:cole)          ! V-wind (m/s)
     real(r8), intent(in)    :: tn(plev, cols:cole)         ! neutral temperature (K)
     real(r8), intent(in)    :: sigma_ped(plev, cols:cole)  ! Pedersen conductivity
     real(r8), intent(in)    :: sigma_hall(plev, cols:cole) ! Hall conductivity
     real(r8), intent(in)    :: te(plev, cols:cole)         ! electron temperature
     real(r8), intent(in)    :: ti(plev, cols:cole)         ! ion temperature
     real(r8), intent(in)    :: mbar(plev, cols:cole)       ! mean molecular weight kg/kmole
     real(r8), intent(in)    :: n2mmr(plev, cols:cole)      ! N2 mass mixing ratio (for oplus)
     real(r8), intent(in)    :: o2mmr(plev, cols:cole)      ! O2 mass mixing ratio (for oplus)
     real(r8), intent(in)    :: o1mmr(plev, cols:cole)      ! O mass mixing ratio (for oplus)
     real(r8), intent(inout) :: o2pmmr(plev, cols:cole)     ! O2+ mass mixing ratio (for oplus)
     real(r8), intent(inout) :: nopmmr(plev, cols:cole)     ! NO+ mass mixing ratio (for oplus)
     real(r8), intent(in)    :: n2pmmr(plev, cols:cole)     ! N2+ mass mixing ratio (for oplus)
     real(r8), intent(inout) :: opmmr(plev, cols:cole)      ! O+ mass mixing ratio (oplus_xport output)
     real(r8), intent(inout) :: opmmrtm1(plev, cols:cole)   ! O+ previous time step (oplus_xport output)
     real(r8), intent(inout) :: ui(plev, cols:cole)         ! zonal ion drift (edynamo or empirical)
     real(r8), intent(inout) :: vi(plev, cols:cole)         ! meridional ion drift (edynamo or empirical)
     real(r8), intent(inout) :: wi(plev, cols:cole)         ! vertical ion drift (edynamo or empirical)
     real(r8), intent(in)    :: rmassN2, rmassO2, rmassO1
     real(r8), intent(in)    :: rmassO2p                    ! O2+ molecular weight kg/kmol
     real(r8), intent(in)    :: rmassNOp                    ! NO+ molecular weight kg/kmol
     real(r8), intent(in)    :: rmassN2p                    ! N2+ molecular weight kg/kmol
     real(r8), intent(in)    :: rmassOp                     ! O+ molecular weight kg/kmol

     real(r8), intent(in)    :: rmassFep                    ! Fe+ molecular weight kg/kmol
     real(r8), intent(in)    :: rmassMgp                    ! Mg+ molecular weight kg/kmol
     real(r8), intent(in)    :: rmassNap                    ! Na+ molecular weight kg/kmol
     real(r8), intent(in)    :: rmassCap                    ! Ca+ molecular weight kg/kmol
     real(r8), intent(in)    :: rmassKp                     ! K+ molecular weight kg/kmol
     real(r8), intent(in)    :: rmassSip                    ! Si+ molecular weight kg/kmol
     real(r8), pointer :: Fepmmr(:,:)
     real(r8), pointer :: Mgpmmr(:,:)
     real(r8), pointer :: Napmmr(:,:)
     real(r8), pointer :: Capmmr(:,:)
     real(r8), pointer :: Kpmmr(:,:)
     real(r8), pointer :: Sipmmr(:,:)
     real(r8), pointer :: NOpmmrtm1(:,:)
     real(r8), pointer :: O2pmmrtm1(:,:)
     real(r8), pointer :: Fepmmrtm1(:,:)
     real(r8), pointer :: Mgpmmrtm1(:,:)
     real(r8), pointer :: Napmmrtm1(:,:)
     real(r8), pointer :: Capmmrtm1(:,:)
     real(r8), pointer :: Kpmmrtm1(:,:)
     real(r8), pointer :: Sipmmrtm1(:,:)

     !
     ! Local:
     !
     integer :: i, j, k, kk
     integer :: kx                  ! Vertical index at peak of F2 layer electron density
     integer :: nstep
     integer :: nfields             ! Number of fields for multi-field calls
     integer :: iyear,imo,iday,tod  ! tod is time-of-day in seconds
     integer :: isplit              ! loop index

     real(r8) :: secs           ! time of day in seconds

     real(r8), parameter :: small = 1.e-25_r8 ! for fields not currently available

     real(r8) :: wn(plev, cols:cole)     ! vertical velocity (from omega)
     real(r8) :: pmid_inv(nlev)                       ! inverted reference pressure at midpoints (Pa)

     real(r8),dimension(plev, cols:cole) :: & ! ion number densities (m^3)
          o2p, nop, n2p, op, ne, optm1
     real(r8),dimension(plev, cols:cole) :: Fep, Mgp, Nap, Cap, Kp, Sip

     real(r8) :: op_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: optm1_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: zpot_in(nlev,lon0:lon1,lat0:lat1) ! geopotential height (cm)
     real(r8) :: ui_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: vi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: wi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: wn_in(nlev,lon0:lon1,lat0:lat1) ! vertical velocity (cm/s)
     real(r8) :: eobx(nlev,lon0:lon1,lat0:lat1) ! zonal  E-field/B-field term
     real(r8) :: eoby(nlev,lon0:lon1,lat0:lat1) ! meridional  E-field/B-field term
     real(r8) :: eobz(nlev,lon0:lon1,lat0:lat1) ! vertical  E-field/B-field term

     real(r8) :: op_out(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: optm1_out(nlev,lon0:lon1,lat0:lat1)

     real(r8),dimension(plev, cols:cole) :: NOptm1, O2ptm1
     real(r8),dimension(plev, cols:cole) :: Feptm1, Mgptm1, Naptm1, Captm1, Kptm1, Siptm1

     real(r8) :: O2pxi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: O2p_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: O2ptm1_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: O2p_out(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: O2ptm1_out(nlev,lon0:lon1,lat0:lat1)

     real(r8) :: NOpxi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: NOp_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: NOptm1_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: NOp_out(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: NOptm1_out(nlev,lon0:lon1,lat0:lat1)

     real(r8) :: Fepxi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Fep_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Feptm1_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Fep_out(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Feptm1_out(nlev,lon0:lon1,lat0:lat1)

     real(r8) :: Mgpxi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Mgp_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Mgptm1_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Mgp_out(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Mgptm1_out(nlev,lon0:lon1,lat0:lat1)

     real(r8) :: Napxi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Nap_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Naptm1_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Nap_out(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Naptm1_out(nlev,lon0:lon1,lat0:lat1)

     real(r8) :: Capxi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Cap_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Captm1_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Cap_out(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Captm1_out(nlev,lon0:lon1,lat0:lat1)

     real(r8) :: Kpxi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Kp_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Kptm1_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Kp_out(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Kptm1_out(nlev,lon0:lon1,lat0:lat1)

     real(r8) :: Sipxi_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Sip_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Siptm1_in(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Sip_out(nlev,lon0:lon1,lat0:lat1)
     real(r8) :: Siptm1_out(nlev,lon0:lon1,lat0:lat1)

     real(r8),target :: halo_tn(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! neutral temperature (deg K)
     real(r8),target :: halo_te(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! electron temperature (deg K)
     real(r8),target :: halo_ti(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! ion temperature (deg K)
     real(r8),target :: halo_un(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! neutral zonal wind (cm/s)
     real(r8),target :: halo_vn(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! neutral meridional wind (cm/s)
     real(r8),target :: halo_om(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! omega (1/s)
     real(r8),target :: halo_o2(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! o2 (mmr)
     real(r8),target :: halo_o1(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! o (mmr)
     real(r8),target :: halo_n2(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! n2 (mmr)
     real(r8),target :: halo_ne(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! n2 (mmr)
     real(r8),target :: halo_op(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! o (mmr)
     real(r8),target :: halo_mbar(nlev,lon0-2:lon1+2,lat0-2:lat1+2) ! mean molecular weight
     real(r8), allocatable :: polesign(:)
    !
     real(r8) :: nmf2(cols:cole) ! Electron number density at F2 peak (m-3 converted to cm-3)
     real(r8) :: hmf2(cols:cole) ! Height of electron number density F2 peak (m converted to km)
     real(r8) :: &
          height(3),    &  ! Surrounding heights when locating electron density F2 peak
          nde(3)           ! Surround densities when locating electron density F2 peak
     real(r8) h12,h22,h32,deltx,atx,ax,btx,bx,ctx,cx ! Variables used for weighting when locating F2 peak
     !
     logical :: do_integrals
!
! Pointers for multiple-field calls:
    type(array_ptr_type),allocatable :: ptrs(:)

     character(len=*), parameter :: subname = 'd_pie_coupling'

     real(r8), dimension(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1) :: & ! 3d fields on mag grid
            ped_mag,     & ! pedersen conductivity on magnetic grid
            hal_mag,     & ! hall conductivity on magnetic grid
            zpot_mag       ! geopotential on magnetic grid

     real(r8), dimension(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1) :: & ! 3d fields on mag grid
            ped_mag_in,     & ! pedersen conductivity on magnetic grid
            hal_mag_in,     & ! hall conductivity on magnetic grid
            zpot_mag_in       ! geopotential on magnetic grid

     real(r8), dimension(lon0:lon1,lat0:lat1,lev0:lev1) :: & ! 3d fields on geo grid
            zpot_geo, &       ! geopotential on magnetic grid
            tn_geo, &
            te_geo, &
            ti_geo, &
            un_geo, &
            vn_geo, &
            wn_geo, &
            ui_geo, &
            vi_geo, &
            wi_geo, &
            omega_geo, &
            o2_geo, &
            o_geo, &
            n2_geo, &
            ne_geo, &
            op_geo, &
            optm1_geo, &
            pmid_geo, &
            mbar_geo

     real(r8), dimension(lon0:lon1,lat0:lat1,lev0:lev1) :: & ! 3d fields on geo grid
          NOp_geo, &
          O2p_geo

     real(r8), dimension(lon0:lon1,lat0:lat1,lev0:lev1) :: & ! 3d fields on geo grid
          NOptm1_geo, &
          O2ptm1_geo, &
          Feptm1_geo, &
          Mgptm1_geo, &
          Naptm1_geo, &
          Captm1_geo, &
          Kptm1_geo, &
          Siptm1_geo, &
          Fep_geo, &
          Mgp_geo, &
          Nap_geo, &
          Cap_geo, &
          Kp_geo, &
          Sip_geo

     real(r8), dimension(lon0:lon1,lat0:lat1,lev0:lev1) :: &
          adotv1_in, adotv2_in
     real(r8), dimension(lon0:lon1,lat0:lat1) :: &
          adota1_in, adota2_in, a1dta2_in, be3_in, sini_in

     real(r8), dimension(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1) :: &
          adotv1_mag, adotv2_mag
     real(r8), dimension(mlon0:mlon1,mlat0:mlat1) :: &
          adota1_mag, adota2_mag, a1dta2_mag, be3_mag, sini_mag

     call t_startf(subname)

     if (debug .and. masterproc) then

        nstep = get_nstep()
        call get_curr_date(iyear, imo, iday, tod)
        secs = real(tod, r8)

        write(iulog,"(3a,i8,a,3i5,a,f6.2)")                                   &
             'Enter ',subname,': nstep = ',nstep, ', iyear,imo,iday = ',      &
             iyear, imo, iday, ' ut (hrs) = ', secs/3600._r8

        write(iulog,"(a,': nspltop = ',i3)") subname, nspltop
     end if

     !---------------------------------------------------------------
     ! Calculate vertical neutral wind velocity wn(i,j,k).
     ! (omega (Pa/s), tn (K), and mbar (kg/kmole) are inputs, grav is m/s^2)
     !---------------------------------------------------------------
     call calc_wn(tn, omega, pmid, mbar, grav, wn, cols, cole, nlev)

     !---------------------------------------------------------------
     ! Convert from mmr to number densities (m^3):
     !---------------------------------------------------------------
     do i = cols, cole
        do k = 1, nlev
           ! O2+, NO+, N2+, O+:
           o2p(k, i) = o2pmmr(k, i) * mbar(k, i) / rmassO2p * &
                pmid(k,i) / (kboltz * tn(k, i))
           nop(k, i) = nopmmr(k, i) * mbar(k, i) / rmassNOp * &
                pmid(k,i) / (kboltz * tn(k, i))
           n2p(k, i) = n2pmmr(k, i) * mbar(k, i) / rmassN2p * &
                pmid(k,i) / (kboltz * tn(k, i))
           op(k, i)  = opmmr(k, i)  * mbar(k, i) / rmassOp  * &
                pmid(k,i) / (kboltz * tn(k, i))
           optm1(k, i)  = opmmrtm1(k, i)  * mbar(k, i) / rmassOp  * &
                pmid(k,i) / (kboltz * tn(k, i))
           ne(k, i) = o2p(k,i)+nop(k,i)+n2p(k,i)+op(k,i)
        end do
     end do ! k=1,nlev

     Fep=fillvalue; Mgp=fillvalue; Nap=fillvalue; Cap=fillvalue; Kp=fillvalue; Sip=fillvalue;

     if (associated(Fepmmr)) then
        Fep(1:nlev,cols:cole) = Fepmmr(1:nlev,cols:cole) * mbar(1:nlev,cols:cole) / rmassFep * &
             pmid(1:nlev,cols:cole) / (kboltz * tn(1:nlev,cols:cole))
        ne(1:nlev,cols:cole) = ne(1:nlev,cols:cole) + Fep(1:nlev,cols:cole)
     end if
     if (associated(Mgpmmr)) then
        Mgp(1:nlev,cols:cole) = Mgpmmr(1:nlev,cols:cole) * mbar(1:nlev,cols:cole) / rmassMgp * &
             pmid(1:nlev,cols:cole) / (kboltz * tn(1:nlev,cols:cole))
        ne(1:nlev,cols:cole) = ne(1:nlev,cols:cole) + Mgp(1:nlev,cols:cole)
     end if
     if (associated(Napmmr)) then
        Nap(1:nlev,cols:cole) = Napmmr(1:nlev,cols:cole) * mbar(1:nlev,cols:cole) / rmassNap * &
             pmid(1:nlev,cols:cole) / (kboltz * tn(1:nlev,cols:cole))
        ne(1:nlev,cols:cole) = ne(1:nlev,cols:cole) + Nap(1:nlev,cols:cole)
     end if
     if (associated(Capmmr)) then
        Cap(1:nlev,cols:cole) = Capmmr(1:nlev,cols:cole) * mbar(1:nlev,cols:cole) / rmassCap * &
             pmid(1:nlev,cols:cole) / (kboltz * tn(1:nlev,cols:cole))
        ne(1:nlev,cols:cole) = ne(1:nlev,cols:cole) + Cap(1:nlev,cols:cole)
     end if
     if (associated(Kpmmr)) then
        Kp(1:nlev,cols:cole) = Kpmmr(1:nlev,cols:cole) * mbar(1:nlev,cols:cole) / rmassKp * &
             pmid(1:nlev,cols:cole) / (kboltz * tn(1:nlev,cols:cole))
        ne(1:nlev,cols:cole) = ne(1:nlev,cols:cole) + Kp(1:nlev,cols:cole)
     end if
     if (associated(Sipmmr)) then
        Sip(1:nlev,cols:cole) = Sipmmr(1:nlev,cols:cole) * mbar(1:nlev,cols:cole) / rmassSip * &
             pmid(1:nlev,cols:cole) / (kboltz * tn(1:nlev,cols:cole))
        ne(1:nlev,cols:cole) = ne(1:nlev,cols:cole) + Sip(1:nlev,cols:cole)
     end if

     if (associated(NOpmmrtm1)) then
        NOptm1(:,:) = NOpmmrtm1(:,:) * mbar(:,:) / rmassOp / pmid(:,:) / (kboltz * tn(:,:))
     end if
     if (associated(O2pmmrtm1)) then
        O2ptm1(:,:) = O2pmmrtm1(:,:) * mbar(:,:) / rmassOp / pmid(:,:) / (kboltz * tn(:,:))
     end if
     if (associated(Fepmmrtm1)) then
        Feptm1(:,:) = Fepmmrtm1(:,:) * mbar(:,:) / rmassOp / pmid(:,:) / (kboltz * tn(:,:))
     end if
     if (associated(Mgpmmrtm1)) then
        Mgptm1(:,:) = Mgpmmrtm1(:,:) * mbar(:,:) / rmassOp / pmid(:,:) / (kboltz * tn(:,:))
     end if
     if (associated(Napmmrtm1)) then
        Naptm1(:,:) = Napmmrtm1(:,:) * mbar(:,:) / rmassOp / pmid(:,:) / (kboltz * tn(:,:))
     end if
     if (associated(Capmmrtm1)) then
        Captm1(:,:) = Capmmrtm1(:,:) * mbar(:,:) / rmassOp / pmid(:,:) / (kboltz * tn(:,:))
     end if
     if (associated(Kpmmrtm1)) then
        Kptm1(:,:) = Kpmmrtm1(:,:) * mbar(:,:) / rmassOp / pmid(:,:) / (kboltz * tn(:,:))
     end if
     if (associated(Sipmmrtm1)) then
        Siptm1(:,:) = Sipmmrtm1(:,:) * mbar(:,:) / rmassOp / pmid(:,:) / (kboltz * tn(:,:))
     end if

     if (debug_hist) then
        call outfld_phys('DPIE_TN',tn)
        call outfld_phys('DPIE_UN',u* 100._r8)
        call outfld_phys('DPIE_VN',v* 100._r8)
        call outfld_phys('DPIE_WN',wn* 100._r8)
        call outfld_phys('DPIE_ZHT',zht* 100._r8)
        call outfld_phys('DPIE_ZGI',zgi* 100._r8)
        call outfld_phys('DPIE_MBAR',mbar)

        call outfld_phys('DPIE_N2',n2mmr)
        call outfld_phys('DPIE_O2',o2mmr)
        call outfld_phys('DPIE_O',o1mmr)

        call outfld_phys('DPIE_OMEGA',omega)
        call outfld_phys('DPIE_OM',-omega/pmid)

        call outfld_phys('DPIE_TE',te)
        call outfld_phys('DPIE_TI',ti)

        call outfld_phys('DPIE_O2P',o2p)
        call outfld_phys('DPIE_NOP',nop)
        call outfld_phys('DPIE_N2P',n2p)

        call outfld_phys('DPIE_FEP',Fep)
        call outfld_phys('DPIE_MGP',Mgp)
        call outfld_phys('DPIE_NAP',Nap)
        call outfld_phys('DPIE_CAP',Cap)
        call outfld_phys('DPIE_KP', Kp)
        call outfld_phys('DPIE_SIP',Sip)

     endif
     call outfld_phys('EDens',ne/1.E6_r8)
     call outfld_phys('OpDens',op/1.E6_r8)
     call outfld_phys('NOpDens',nop/1.E6_r8)
     call outfld_phys('O2pDens',o2p/1.E6_r8)

     call outfld_phys('FepDens',Fep/1.E6_r8)
     call outfld_phys('MgpDens',Mgp/1.E6_r8)
     call outfld_phys('NapDens',Nap/1.E6_r8)
     call outfld_phys('CapDens',Cap/1.E6_r8)
     call outfld_phys('KpDens' ,Kp/1.E6_r8)
     call outfld_phys('SipDens',Sip/1.E6_r8)

     !-------------------------------------------------------------------------
     !  Derive diagnostics nmF2 and hmF2 for output based on TIE-GCM algorithm
     !-------------------------------------------------------------------------
     if (hist_fld_active('HMF2') .or. hist_fld_active('NMF2')) then
        iloop: do i = cols, cole
           kx = 0
           kloop: do k= 2, nlev
              if (ne(k,i) >= ne(k-1,i) .and. ne(k,i) >= ne(k+1,i)) then
                 kx = k
                 exit kloop
              end if
           end do kloop

           if (kx==0) then
              hmf2(i) = fillvalue
              nmf2(i) = fillvalue
              exit iloop
           end if

           height = (/zht(kx+1,i),zht(kx,i),zht(kx-1,i)/)
           nde = (/ne(kx+1,i),ne(kx,i),ne(kx-1,i)/)

           h12 = height(1)*height(1)
           h22 = height(2)*height(2)
           h32 = height(3)*height(3)

           deltx=h12*height(2)+h22*height(3)+h32*height(1)-h32*height(2)-h12*height(3)-h22*height(1)

           atx=nde(1)*height(2)+nde(2)*height(3)+nde(3)*height(1)-height(2)*nde(3)-height(3)*nde(1)-height(1)*nde(2)
           ax=atx/deltx

           btx=h12*nde(2)+h22*nde(3)+h32*nde(1)-h32*nde(2)-h12*nde(3)-h22*nde(1)
           bx=btx/deltx
           ctx=h12*height(2)*nde(3)+h22*height(3)*nde(1)+h32*height(1)*nde(2)-h32*height(2)*nde(1)- &
                h12*height(3)*nde(2)-h22*height(1)*nde(3)
           cx=ctx/deltx

           hmf2(i)=-(bx/(2._r8*ax)) * 1.E-03_r8
           nmf2(i)=-((bx*bx-4._r8*ax*cx)/(4._r8*ax)) * 1.E-06_r8

        end do iloop ! i=cols, cole

        call outfld_phys1d('HMF2',hmf2)
        call outfld_phys1d('NMF2',nmf2)
     end if
     if (debug_hist) then
        call outfld_phys('DPIE_OPMMR', opmmr)
        call outfld_phys('DPIE_FePMMR', fepmmr)
        call outfld_phys('DPIE_MgPMMR', mgpmmr)
        call outfld_phys('PED_phys', sigma_ped )
        call outfld_phys('HAL_phys', sigma_hall )
     endif
     if (ionos_edyn_active .or. ionos_oplus_xport) then

        call regrid_phys2geo_3d( zgi,zpot_geo, plev, cols, cole )
        call regrid_phys2geo_3d( u, un_geo, plev, cols, cole )
        call regrid_phys2geo_3d( v, vn_geo, plev, cols, cole )
        call regrid_phys2geo_3d( wn,wn_geo, plev, cols, cole )
        call regrid_phys2geo_3d( ui, ui_geo, plev, cols, cole )
        call regrid_phys2geo_3d( vi, vi_geo, plev, cols, cole )
        call regrid_phys2geo_3d( wi, wi_geo, plev, cols, cole )

        do k = 1, nlev
           kk = nlev-k+1
           do j = lat0, lat1
              do i = lon0, lon1
                 zpot_in(kk,i,j) = zpot_geo(i,j,k) * 100._r8 ! m -> cm
                 halo_un(kk,i,j) = un_geo(i,j,k) * 100._r8 ! m/s -> cm/s
                 halo_vn(kk,i,j) = vn_geo(i,j,k) * 100._r8 ! m/s -> cm/s
                 wn_in(kk,i,j)   = wn_geo(i,j,k) * 100._r8 ! m/s -> cm/s
                 ui_in(kk,i,j)   = ui_geo(i,j,k) * 100._r8 ! zonal ion drift (m/s -> cm/s)
                 vi_in(kk,i,j)   = vi_geo(i,j,k) * 100._r8 ! meridional ion drift (m/s -> cm/s)
                 wi_in(kk,i,j)   = wi_geo(i,j,k) * 100._r8 ! vertical ion drift (m/s -> cm/s)
              end do
           end do
        end do

     end if

    !
    !
    ! Call electrodynamo (edynamo.F90)
    ! If using time3d conductances, tell dynamo to *not* do fieldline
    ! integrations (i.e., do_integrals == false). In this case, edynamo
    ! conductances zigmxx,rim1,2 from time3d will be set by subroutine
    ! transform_glbin in time3d module.
    !
    do_integrals = .true.
    !
    ! If ionos_edyn_active=false, then empirical ion drifts were passed in from physics,
    ! otherwise dynamo calculates them here, and they will be passed to physics.
    !
    if (ionos_edyn_active) then

       call t_startf('dpie_ionos_dynamo')

       call calc_adotv( zpot_in(lev0:lev1,lon0:lon1,lat0:lat1), &
            halo_un(lev0:lev1,lon0:lon1,lat0:lat1), &
            halo_vn(lev0:lev1,lon0:lon1,lat0:lat1), &
            wn_in(lev0:lev1,lon0:lon1,lat0:lat1), &
            adotv1_in, adotv2_in, adota1_in, adota2_in, &
            a1dta2_in, be3_in, sini_in, lev0, lev1, lon0, lon1, lat0, lat1)

       call regrid_geo2mag_3d( adotv1_in, adotv1_mag )
       call regrid_geo2mag_3d( adotv2_in, adotv2_mag )
       if (debug_hist) then
          call outfld_geo('EDYN_ADOTV1', adotv1_in(:,:,lev1:lev0:-1) )
          call outfld_geo('EDYN_ADOTV2', adotv2_in(:,:,lev1:lev0:-1) )

          call outfld_geo2d( 'EDYN_ADOTA1', adota1_in )
          call outfld_geo2d( 'EDYN_ADOTA2', adota2_in )
          call outfld_geo2d( 'EDYN_A1DTA2', a1dta2_in )
          call outfld_geo2d( 'EDYN_BE3' , be3_in )
          call outfld_geo2d( 'EDYN_SINI', sini_in )
       endif
       call regrid_geo2mag_2d( adota1_in, adota1_mag )
       call regrid_geo2mag_2d( adota2_in, adota2_mag )
       call regrid_geo2mag_2d( a1dta2_in, a1dta2_mag )
       call regrid_geo2mag_2d( be3_in, be3_mag )
       call regrid_geo2mag_2d( sini_in, sini_mag )
       if (debug_hist) then
          call outfld_mag2d('ADOTA1_MAG', adota1_mag )
          call outfld_mag2d('SINI_MAG', sini_mag )
       endif
       call regrid_phys2mag_3d( sigma_ped, ped_mag, plev, cols, cole )
       call regrid_phys2mag_3d( sigma_hall, hal_mag, plev, cols, cole )
       call regrid_phys2mag_3d( zgi, zpot_mag, plev, cols, cole )

       if (mytid<ntask) then
          zpot_mag_in(:,:,mlev0:mlev1) = zpot_mag(:,:,mlev1:mlev0:-1) * 100._r8 ! m -> cm
          ped_mag_in(:,:,mlev0:mlev1) = ped_mag(:,:,mlev1:mlev0:-1)
          hal_mag_in(:,:,mlev0:mlev1) = hal_mag(:,:,mlev1:mlev0:-1)

          call  dynamo( zpot_mag_in, ped_mag_in, hal_mag_in, adotv1_mag, adotv2_mag, adota1_mag, &
               adota2_mag, a1dta2_mag, be3_mag, sini_mag,  &
               zpot_in, ui_in, vi_in, wi_in, eobx,eoby,eobz, &
               lon0,lon1, lat0,lat1, lev0,lev1, do_integrals )
       endif

       call t_stopf ('dpie_ionos_dynamo')

    else
       if (debug .and. masterproc) then
          write(iulog,"('dpie_coupling (dynamo NOT called): nstep=',i8)") nstep
          write(iulog,"('  empirical ExB ui min,max (cm/s)=',2es12.4)")       &
               minval(ui),maxval(ui)
          write(iulog,"('  empirical ExB vi min,max (cm/s)=',2es12.4)")       &
               minval(vi),maxval(vi)
          write(iulog,"('  empirical ExB wi min,max (cm/s)=',2es12.4)")       &
               minval(wi),maxval(wi)
       end if
    end if

    !
    ! Call O+ transport routine.  Now all inputs to oplus_xport should be in
    ! tiegcm-format wrt longitude (-180->180), vertical (bot2top), and units (CGS).
    ! (Composition is mmr, ne is cm^3, winds are cm/s)
    ! Output op_out and opnm_out will be in cm^3, converted to mmr below.
    !
    if (ionos_oplus_xport) then
       pmid_inv(1:nlev) = pref_mid(nlev:1:-1) ! invert ref pressure (Pa) as in tiegcm


       !
       ! Transport O+ (all args in 'TIEGCM format')
       ! Subcycle oplus_xport nspltop times.
       !
       if (debug .and. masterproc) then
          write(iulog,"(a,i8,a,i3)")                                         &
               'dpie_coupling before subcycling oplus_xport: nstep = ',      &
               nstep, ' nspltop = ', nspltop
       end if

       call regrid_phys2geo_3d( tn, tn_geo, plev, cols, cole )
       call regrid_phys2geo_3d( te, te_geo, plev, cols, cole )
       call regrid_phys2geo_3d( ti, ti_geo, plev, cols, cole )
       call regrid_phys2geo_3d( omega, omega_geo, plev, cols, cole )
       call regrid_phys2geo_3d( o2mmr, o2_geo, plev, cols, cole )
       call regrid_phys2geo_3d( n2mmr, n2_geo, plev, cols, cole )
       call regrid_phys2geo_3d( o1mmr, o_geo, plev, cols, cole )
       call regrid_phys2geo_3d( op, op_geo, plev, cols, cole )
       call regrid_phys2geo_3d( optm1, optm1_geo, plev, cols, cole )
       call regrid_phys2geo_3d( pmid, pmid_geo, plev, cols, cole )
       call regrid_phys2geo_3d( mbar, mbar_geo, plev, cols, cole )
       call regrid_phys2geo_3d( ne, ne_geo, plev, cols, cole )

       if (associated(NOpmmrtm1)) then
          call regrid_phys2geo_3d( NOp, NOp_geo, plev, cols, cole )
          call regrid_phys2geo_3d( NOptm1, NOptm1_geo, plev, cols, cole )
       end if
       if (associated(O2pmmrtm1)) then
          call regrid_phys2geo_3d( O2p, O2p_geo, plev, cols, cole )
          call regrid_phys2geo_3d( O2ptm1, O2ptm1_geo, plev, cols, cole )
       end if

       if (associated(Fepmmrtm1)) then
          call regrid_phys2geo_3d( Feptm1, Feptm1_geo, plev, cols, cole )
       end if
       if (associated(Mgpmmrtm1)) then
          call regrid_phys2geo_3d( Mgptm1, Mgptm1_geo, plev, cols, cole )
       end if
       if (associated(Napmmrtm1)) then
          call regrid_phys2geo_3d( Naptm1, Naptm1_geo, plev, cols, cole )
       end if
       if (associated(Capmmrtm1)) then
          call regrid_phys2geo_3d( Captm1, Captm1_geo, plev, cols, cole )
       end if
       if (associated(Kpmmrtm1)) then
          call regrid_phys2geo_3d( Kptm1, Kptm1_geo, plev, cols, cole )
       end if
       if (associated(Sipmmrtm1)) then
          call regrid_phys2geo_3d( Siptm1, Siptm1_geo, plev, cols, cole )
       end if

       if (associated(Fepmmr)) then
          call regrid_phys2geo_3d( Fep, Fep_geo, plev, cols, cole )
       end if
       if (associated(Mgpmmr)) then
          call regrid_phys2geo_3d( Mgp, Mgp_geo, plev, cols, cole )
       end if
       if (associated(Napmmr)) then
          call regrid_phys2geo_3d( Nap, Nap_geo, plev, cols, cole )
       end if
       if (associated(Capmmr)) then
          call regrid_phys2geo_3d( Cap, Cap_geo, plev, cols, cole )
       end if
       if (associated(Kpmmr)) then
          call regrid_phys2geo_3d( Kp, Kp_geo, plev, cols, cole )
       end if
       if (associated(Sipmmr)) then
          call regrid_phys2geo_3d( Sip, Sip_geo, plev, cols, cole )
       end if


       if (mytid<ntask) then

          call t_startf('dpie_halo')
          do k = 1, nlev
             kk = nlev-k+1
             do j = lat0, lat1
                do i = lon0, lon1
                   halo_tn(kk,i,j)   = tn_geo(i,j,k)
                   halo_te(kk,i,j)   = te_geo(i,j,k)
                   halo_ti(kk,i,j)   = ti_geo(i,j,k)
                   halo_om(kk,i,j)   = -omega_geo(i,j,k) / pmid_geo(i,j,k) ! Pa/s -> 1/s
                   halo_o2(kk,i,j)   = o2_geo(i,j,k)
                   halo_o1(kk,i,j)   = o_geo(i,j,k)
                   halo_n2(kk,i,j)   = n2_geo(i,j,k)
                   halo_ne(kk,i,j)   = ne_geo(i,j,k)
                   halo_op(kk,i,j)   = op_geo(i,j,k)
                   halo_mbar(kk,i,j) = mbar_geo(i,j,k)
                   op_in(kk,i,j)     = op_geo(i,j,k) / 1.e6_r8  ! m^3 -> cm^3
                   optm1_in(kk,i,j)  = optm1_geo(i,j,k) / 1.e6_r8  ! m^3 -> cm^3
                end do
             end do
             if (associated(NOpmmrtm1)) then
                NOp_in(kk,:,:)    = NOp_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
                NOptm1_in(kk,:,:) = NOptm1_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
             endif
             if (associated(O2pmmrtm1)) then
                O2p_in(kk,:,:)    = O2p_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
                O2ptm1_in(kk,:,:) = O2ptm1_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
             endif
             if (associated(Fepmmrtm1)) then
                Fep_in(kk,:,:)    = Fep_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
                Feptm1_in(kk,:,:) = Feptm1_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
             endif
             if (associated(Mgpmmrtm1)) then
                Mgp_in(kk,:,:)    = Mgp_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
                Mgptm1_in(kk,:,:) = Mgptm1_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
             endif
             if (associated(Napmmrtm1)) then
                Nap_in(kk,:,:)    = Nap_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
                Naptm1_in(kk,:,:) = Naptm1_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
             endif
             if (associated(Capmmrtm1)) then
                Cap_in(kk,:,:)    = Cap_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
                Captm1_in(kk,:,:) = Captm1_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
             endif
             if (associated(Kpmmrtm1)) then
                Kp_in(kk,:,:)    = Kp_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
                Kptm1_in(kk,:,:) = Kptm1_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
             endif
             if (associated(Sipmmrtm1)) then
                Sip_in(kk,:,:)    = Sip_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
                Siptm1_in(kk,:,:) = Siptm1_geo(:,:,k) / 1.e6_r8  ! m^3 -> cm^3
             endif

          end do
          !
          ! Define halo points on inputs:
          ! WACCM has global longitude values at the poles (j=1,j=nlev)
          ! (they are constant for most, except the winds.)
          !
          ! Set two halo points in lat,lon:
          !
          nfields = 12
          allocate(ptrs(nfields), polesign(nfields))
          ptrs(1)%ptr => halo_tn ; ptrs(2)%ptr => halo_te ; ptrs(3)%ptr => halo_ti
          ptrs(4)%ptr => halo_un ; ptrs(5)%ptr => halo_vn ; ptrs(6)%ptr => halo_om
          ptrs(7)%ptr => halo_o2 ; ptrs(8)%ptr => halo_o1 ; ptrs(9)%ptr => halo_n2
          ptrs(10)%ptr => halo_mbar
          ptrs(11)%ptr => halo_ne
          ptrs(12)%ptr => halo_op

          polesign = 1._r8
          polesign(4:5) = -1._r8 ! un,vn

          call mp_geo_halos(ptrs,1,nlev,lon0,lon1,lat0,lat1,nfields)
          !
          ! Set latitude halo points over the poles (this does not change the poles).
          ! (the 2nd halo over the poles will not actually be used (assuming lat loops
          !  are lat=2,plat-1), because jp1,jm1 will be the pole itself, and jp2,jm2
          !  will be the first halo over the pole)
          !
          call mp_pole_halos(ptrs,1,nlev,lon0,lon1,lat0,lat1,nfields,polesign)
          deallocate(ptrs,polesign)
          call t_stopf('dpie_halo')
          if (debug_hist) then
             call outfld_geokij( 'OPtm1i',optm1_in, lev0,lev1, lon0,lon1, lat0,lat1 )
          endif
          call t_startf('dpie_oplus_xport')

          call outfld_geokij( 'OPi',op_in, lev0,lev1, lon0,lon1, lat0,lat1 )

          do isplit = 1, nspltop

             if (isplit > 1) then
                op_in = op_out
                optm1_in = optm1_out
             end if

             call oplus_xport(halo_tn, halo_te, halo_ti, halo_un, halo_vn, halo_om, &
                  zpot_in, halo_o2, halo_o1, halo_n2, op_in, optm1_in, &
                  halo_mbar, ui_in, vi_in, wi_in, pmid_inv, &
                  op_out, optm1_out, &
                  lon0, lon1, lat0, lat1, nspltop, isplit)

          end do ! isplit=1,nspltop

          call outfld_geokij( 'OPo',op_out, lev0,lev1, lon0,lon1, lat0,lat1 )

          if (associated(O2pmmrtm1)) then
             call calc_xi(halo_o2(:,lon0:lon1,lat0:lat1),halo_o1(:,lon0:lon1,lat0:lat1),halo_n2(:,lon0:lon1,lat0:lat1),&
                          rmassO2p, rmassN2, rmassO2, rmassO1, O2pxi_in, lon0, lon1, lat0, lat1, nlev)
             do isplit = 1,nsplt_vxb
                if (isplit > 1) then
                   O2p_in = O2p_out
                   O2ptm1_in = O2ptm1_out
                end if
                call vxb_xport(halo_tn, halo_te, halo_ti, halo_ne, halo_op, halo_un, halo_vn, halo_om,  &
                     zpot_in, halo_o2, halo_o1, halo_n2, O2p_in, O2ptm1_in, O2pxi_in,       &
                     halo_mbar, ui_in, vi_in, wi_in,  eobx, eoby, eobz, pmid_inv,    &
                     O2p_out, O2ptm1_out, rmassO2p, lon0, lon1, lat0, lat1, nsplt_vxb, isplit )
             end do
          endif

          if (associated(NOpmmrtm1)) then
             call calc_xi(halo_o2(:,lon0:lon1,lat0:lat1),halo_o1(:,lon0:lon1,lat0:lat1),halo_n2(:,lon0:lon1,lat0:lat1),&
                          rmassNOp, rmassN2, rmassO2, rmassO1, NOpxi_in, lon0, lon1, lat0, lat1, nlev)
             do isplit = 1,nsplt_vxb
                if (isplit > 1) then
                   NOp_in = NOp_out
                   NOptm1_in = NOptm1_out
                end if
                call vxb_xport(halo_tn, halo_te, halo_ti, halo_ne, halo_op, halo_un, halo_vn, halo_om,  &
                     zpot_in, halo_o2, halo_o1, halo_n2, NOp_in, NOptm1_in, NOpxi_in,       &
                     halo_mbar, ui_in, vi_in, wi_in,  eobx, eoby, eobz, pmid_inv,    &
                     NOp_out, NOptm1_out, rmassNOp, lon0, lon1, lat0, lat1, nsplt_vxb, isplit )
             end do
          endif


          if (associated(Fepmmrtm1)) then
             if (debug_hist) then
                call outfld_geokij( 'FePtm1i',Feptm1_in, lev0,lev1, lon0,lon1, lat0,lat1 )
             endif
             call outfld_geokij( 'FePi',Fep_in, lev0,lev1, lon0,lon1, lat0,lat1 )

             call calc_xi(halo_o2(:,lon0:lon1,lat0:lat1),halo_o1(:,lon0:lon1,lat0:lat1),halo_n2(:,lon0:lon1,lat0:lat1),&
                          rmassFep, rmassN2, rmassO2, rmassO1, Fepxi_in, lon0, lon1, lat0, lat1, nlev)
             do isplit = 1,nsplt_vxb
                if (isplit > 1) then
                   Fep_in = Fep_out
                   Feptm1_in = Feptm1_out
                end if

                call vxb_xport(halo_tn, halo_te, halo_ti, halo_ne, halo_op, halo_un, halo_vn, halo_om,  &
                     zpot_in, halo_o2, halo_o1, halo_n2, Fep_in, Feptm1_in, Fepxi_in,       &
                     halo_mbar, ui_in, vi_in, wi_in,  eobx, eoby, eobz, pmid_inv,    &
                     Fep_out, Feptm1_out, rmassFep, lon0, lon1, lat0, lat1, nsplt_vxb, isplit )
             end do

             call outfld_geokij( 'FePo',Fep_out, lev0,lev1, lon0,lon1, lat0,lat1 )


          endif

          if (associated(Mgpmmrtm1)) then
             if (debug_hist) then
                call outfld_geokij( 'MgPtm1i',Mgptm1_in, lev0,lev1, lon0,lon1, lat0,lat1 )
             endif

             call outfld_geokij( 'MgPi',Mgp_in, lev0,lev1, lon0,lon1, lat0,lat1 )

             call calc_xi(halo_o2(:,lon0:lon1,lat0:lat1),halo_o1(:,lon0:lon1,lat0:lat1),halo_n2(:,lon0:lon1,lat0:lat1),&
                          rmassMgp, rmassN2, rmassO2, rmassO1, Mgpxi_in, lon0, lon1, lat0, lat1, nlev)
             do isplit = 1,nsplt_vxb
                if (isplit > 1) then
                   Mgp_in = Mgp_out
                   Mgptm1_in = Mgptm1_out
                end if
                call vxb_xport(halo_tn, halo_te, halo_ti, halo_ne, halo_op, halo_un, halo_vn, halo_om,  &
                     zpot_in, halo_o2, halo_o1, halo_n2, Mgp_in, Mgptm1_in, Mgpxi_in,       &
                     halo_mbar, ui_in, vi_in, wi_in,  eobx, eoby, eobz, pmid_inv,    &
                     Mgp_out, Mgptm1_out, rmassMgp, lon0, lon1, lat0, lat1, nsplt_vxb, isplit )
             end do


             call outfld_geokij( 'MgPo',Mgp_out, lev0,lev1, lon0,lon1, lat0,lat1 )

          endif
          if (associated(Napmmrtm1)) then
             call calc_xi(halo_o2(:,lon0:lon1,lat0:lat1),halo_o1(:,lon0:lon1,lat0:lat1),halo_n2(:,lon0:lon1,lat0:lat1),&
                          rmassNap, rmassN2, rmassO2, rmassO1, Napxi_in, lon0, lon1, lat0, lat1, nlev)
             do isplit = 1,nsplt_vxb
                if (isplit > 1) then
                   Nap_in = Nap_out
                   Naptm1_in = Naptm1_out
                end if
                call vxb_xport(halo_tn, halo_te, halo_ti, halo_ne, halo_op, halo_un, halo_vn, halo_om,  &
                     zpot_in, halo_o2, halo_o1, halo_n2, Nap_in, Naptm1_in, Napxi_in,       &
                     halo_mbar, ui_in, vi_in, wi_in,  eobx, eoby, eobz, pmid_inv,    &
                     Nap_out, Naptm1_out, rmassNap, lon0, lon1, lat0, lat1, nsplt_vxb, isplit )
             end do
          endif
          if (associated(Capmmrtm1)) then
             call calc_xi(halo_o2(:,lon0:lon1,lat0:lat1),halo_o1(:,lon0:lon1,lat0:lat1),halo_n2(:,lon0:lon1,lat0:lat1),&
                          rmassCap, rmassN2, rmassO2, rmassO1, Capxi_in, lon0, lon1, lat0, lat1, nlev)
             do isplit = 1,nsplt_vxb
                if (isplit > 1) then
                   Cap_in = Cap_out
                   Captm1_in = Captm1_out
                end if
                call vxb_xport(halo_tn, halo_te, halo_ti, halo_ne, halo_op, halo_un, halo_vn, halo_om,  &
                     zpot_in, halo_o2, halo_o1, halo_n2, Cap_in, Captm1_in, Capxi_in,       &
                     halo_mbar, ui_in, vi_in, wi_in,  eobx, eoby, eobz, pmid_inv,    &
                     Cap_out, Captm1_out, rmassCap, lon0, lon1, lat0, lat1, nsplt_vxb, isplit )
             end do
          endif
          if (associated(Kpmmrtm1)) then
             call calc_xi(halo_o2(:,lon0:lon1,lat0:lat1),halo_o1(:,lon0:lon1,lat0:lat1),halo_n2(:,lon0:lon1,lat0:lat1),&
                          rmassKp, rmassN2, rmassO2, rmassO1, Kpxi_in, lon0, lon1, lat0, lat1, nlev)
             do isplit = 1,nsplt_vxb
                if (isplit > 1) then
                   Kp_in = Kp_out
                   Kptm1_in = Kptm1_out
                end if
                call vxb_xport(halo_tn, halo_te, halo_ti, halo_ne, halo_op, halo_un, halo_vn, halo_om,  &
                     zpot_in, halo_o2, halo_o1, halo_n2, Kp_in, Kptm1_in, Kpxi_in,       &
                     halo_mbar, ui_in, vi_in, wi_in,  eobx, eoby, eobz, pmid_inv,    &
                     Kp_out, Kptm1_out, rmassKp, lon0, lon1, lat0, lat1, nsplt_vxb, isplit )
             end do
          endif
          if (associated(Sipmmrtm1)) then
             call calc_xi(halo_o2(:,lon0:lon1,lat0:lat1),halo_o1(:,lon0:lon1,lat0:lat1),halo_n2(:,lon0:lon1,lat0:lat1),&
                          rmassSip, rmassN2, rmassO2, rmassO1, Sipxi_in, lon0, lon1, lat0, lat1, nlev)
             do isplit = 1,nsplt_vxb
                if (isplit > 1) then
                   Sip_in = Sip_out
                   Siptm1_in = Siptm1_out
                end if
                call vxb_xport(halo_tn, halo_te, halo_ti, halo_ne, halo_op, halo_un, halo_vn, halo_om,  &
                     zpot_in, halo_o2, halo_o1, halo_n2, Sip_in, Siptm1_in, Sipxi_in,       &
                     halo_mbar, ui_in, vi_in, wi_in,  eobx, eoby, eobz, pmid_inv,    &
                     Sip_out, Siptm1_out, rmassSip, lon0, lon1, lat0, lat1, nsplt_vxb, isplit )
             end do
          endif

          call t_stopf ('dpie_oplus_xport')
          if (debug.and.masterproc) then
             write(iulog,"('dpie_coupling after subcycling oplus_xport: nstep=',i8,' nspltop=',i3)") &
                  nstep,nspltop
             write(iulog,"('  op_out   min,max (cm^3)=',2es12.4)") minval(op_out)   ,maxval(op_out)
             write(iulog,"(' optm1_out min,max (cm^3)=',2es12.4)") minval(optm1_out),maxval(optm1_out)
          end if
          if (debug_hist) then
             call outfld_geokij( 'OPLUS', op_out, lev0,lev1, lon0,lon1, lat0,lat1 )
             call outfld_geokij( 'OPtm1o',optm1_out, lev0,lev1, lon0,lon1, lat0,lat1 )
             call outfld_geokij( 'FePtm1o',Feptm1_out, lev0,lev1, lon0,lon1, lat0,lat1 )
             call outfld_geokij( 'MgPtm1o',Mgptm1_out, lev0,lev1, lon0,lon1, lat0,lat1 )
          endif
          !
          ! Pass new O+ for current and previous time step back to physics (convert from cm^3 to m^3 and back to mmr).
          !
          do k=1,nlev
             kk = nlev-k+1
             do j = lat0,lat1
                do i = lon0,lon1
                   op_geo(i,j,k) = op_out(kk,i,j)*1.e6_r8 * rmassOp / mbar_geo(i,j,k) * &
                        (kboltz * tn_geo(i,j,k)) / pmid_geo(i,j,k)
                   optm1_geo(i,j,k) = optm1_out(kk,i,j)*1.e6_r8 * rmassOp / mbar_geo(i,j,k) * &
                        (kboltz * tn_geo(i,j,k)) / pmid_geo(i,j,k)
                   ui_geo(i,j,k) = ui_in(kk,i,j)/100._r8 ! cm/s -> m/s
                   vi_geo(i,j,k) = vi_in(kk,i,j)/100._r8 ! cm/s -> m/s
                   wi_geo(i,j,k) = wi_in(kk,i,j)/100._r8 ! cm/s -> m/s
                end do
             end do
             if (associated(NOpmmrtm1)) then
                NOp_geo(:,:,k) = NOp_out(kk,:,:)*1.e6_r8 * rmassNOp / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
                NOptm1_geo(:,:,k) = NOptm1_out(kk,:,:)*1.e6_r8 * rmassNOp / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
             end if
             if (associated(O2pmmrtm1)) then
                O2p_geo(:,:,k) = O2p_out(kk,:,:)*1.e6_r8 * rmassO2p / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
                O2ptm1_geo(:,:,k) = O2ptm1_out(kk,:,:)*1.e6_r8 * rmassO2p / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
             end if

             if (associated(Fepmmrtm1)) then
                Fep_geo(:,:,k) = Fep_out(kk,:,:)*1.e6_r8 * rmassFep / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
                Feptm1_geo(:,:,k) = Feptm1_out(kk,:,:)*1.e6_r8 * rmassFep / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
             end if
             if (associated(Mgpmmrtm1)) then
                Mgp_geo(:,:,k) = Mgp_out(kk,:,:)*1.e6_r8 * rmassMgp / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
                Mgptm1_geo(:,:,k) = Mgptm1_out(kk,:,:)*1.e6_r8 * rmassMgp / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
             end if
             if (associated(Napmmrtm1)) then
                Nap_geo(:,:,k) = Nap_out(kk,:,:)*1.e6_r8 * rmassNap / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
                Naptm1_geo(:,:,k) = Naptm1_out(kk,:,:)*1.e6_r8 * rmassNap / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
             end if
             if (associated(Capmmrtm1)) then
                Cap_geo(:,:,k) = Cap_out(kk,:,:)*1.e6_r8 * rmassCap / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
                Captm1_geo(:,:,k) = Captm1_out(kk,:,:)*1.e6_r8 * rmassCap / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
             end if
             if (associated(Kpmmrtm1)) then
                Kp_geo(:,:,k) = Kp_out(kk,:,:)*1.e6_r8 * rmassKp / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
                Kptm1_geo(:,:,k) = Kptm1_out(kk,:,:)*1.e6_r8 * rmassKp / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
             end if
             if (associated(Sipmmrtm1)) then
                Sip_geo(:,:,k) = Sip_out(kk,:,:)*1.e6_r8 * rmassSip / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
                Siptm1_geo(:,:,k) = Siptm1_out(kk,:,:)*1.e6_r8 * rmassSip / mbar_geo(:,:,k) * &
                                 (kboltz * tn_geo(:,:,k)) / pmid_geo(:,:,k)
             end if

          end do

       endif

       call regrid_geo2phys_3d( op_geo, opmmr, plev, cols, cole )
       call regrid_geo2phys_3d( optm1_geo, opmmrtm1, plev, cols, cole )
       call regrid_geo2phys_3d( ui_geo, ui, plev, cols, cole )
       call regrid_geo2phys_3d( vi_geo, vi, plev, cols, cole )
       call regrid_geo2phys_3d( wi_geo, wi, plev, cols, cole )

       if (associated(NOpmmrtm1)) then
          call regrid_geo2phys_3d( NOp_geo, NOpmmr, plev, cols, cole )
          call regrid_geo2phys_3d( NOptm1_geo, NOpmmrtm1, plev, cols, cole )
       end if
       if (associated(O2pmmrtm1)) then
          call regrid_geo2phys_3d( O2p_geo, O2pmmr, plev, cols, cole )
          call regrid_geo2phys_3d( O2ptm1_geo, O2pmmrtm1, plev, cols, cole )
       end if

       if (associated(Fepmmrtm1)) then
          call regrid_geo2phys_3d( Fep_geo, Fepmmr, plev, cols, cole )
          call regrid_geo2phys_3d( Feptm1_geo, Fepmmrtm1, plev, cols, cole )
       end if
       if (associated(Mgpmmrtm1)) then
          call regrid_geo2phys_3d( Mgp_geo, Mgpmmr, plev, cols, cole )
          call regrid_geo2phys_3d( Mgptm1_geo, Mgpmmrtm1, plev, cols, cole )
       end if
       if (associated(Napmmrtm1)) then
          call regrid_geo2phys_3d( Nap_geo, Napmmr, plev, cols, cole )
          call regrid_geo2phys_3d( Naptm1_geo, Napmmrtm1, plev, cols, cole )
       end if
       if (associated(Capmmrtm1)) then
          call regrid_geo2phys_3d( Cap_geo, Capmmr, plev, cols, cole )
          call regrid_geo2phys_3d( Captm1_geo, Capmmrtm1, plev, cols, cole )
       end if
       if (associated(Kpmmrtm1)) then
          call regrid_geo2phys_3d( Kp_geo, Kpmmr, plev, cols, cole )
          call regrid_geo2phys_3d( Kptm1_geo, Kpmmrtm1, plev, cols, cole )
       end if
       if (associated(Sipmmrtm1)) then
          call regrid_geo2phys_3d( Sip_geo, Sipmmr, plev, cols, cole )
          call regrid_geo2phys_3d( Siptm1_geo, Sipmmrtm1, plev, cols, cole )
       end if

    end if ! ionos_oplus_xport

    if (debug_hist) then
       call outfld_phys('WACCM_UI',ui)
       call outfld_phys('WACCM_VI',vi)
       call outfld_phys('WACCM_WI',wi)
       call outfld_phys('WACCM_OP',opmmr)
       call outfld_phys('WACCM_FeP',fepmmr)
       call outfld_phys('WACCM_MgP',mgpmmr)
    endif
    call t_stopf('d_pie_coupling')

  end subroutine d_pie_coupling
  !-----------------------------------------------------------------------
  subroutine calc_wn(tn,omega,pmid,mbar,grav,wn,cols,cole,nlev)
     use shr_const_mod,only : shr_const_rgas ! Universal gas constant
     !
     ! Calculate neutral vertical wind on midpoints (m/s)
     !
     ! Inputs:
     integer,intent(in) :: cols, cole, nlev
     real(r8),dimension(nlev, cols:cole),intent(in) :: &
          tn,   &  ! neutral temperature (deg K)
          omega,&  ! pressure velocity (Pa/s)
          mbar     ! mean molecular weight
     real(r8),dimension(nlev,cols:cole),intent(in) :: &
          pmid     ! pressure at midpoints (Pa)
     real(r8),intent(in) :: grav ! m/s^2
     !
     ! Output:
     real(r8),intent(out) :: wn(nlev, cols:cole) ! vertical velocity output (m/s)
     !
     ! Local:
     integer :: i,k
     real(r8) :: scheight(nlev, cols:cole) ! dimensioned for vectorization

     do i = cols, cole
        do k = 1, nlev
           scheight(k,i) = shr_const_rgas*tn(k,i)/(mbar(k,i)*grav)
           wn(k,i) = -omega(k,i)*scheight(k,i)/pmid(k,i)
        end do
     end do
  end subroutine calc_wn
  !-----------------------------------------------------------------------
  subroutine calc_pfrac(sunlon,pfrac)
    !
    ! Calculate pfrac fractional presence of dynamo equation using critical
    !  convection colatitudes crit(2).
    !
    use edyn_maggrid  ,only: nmlonp1,ylonm,ylatm
    use edyn_solve    ,only: nmlat0
    use aurora_params ,only: offc, dskofc, theta0, aurora_params_set

    implicit none
    !
    ! Args:
    real(r8),intent(in) :: sunlon  ! Sun's longitude in dipole coordinates
    !
    ! Output: fractional presence of dynamo equation using critical colatitudes
    !
    real(r8),intent(out) :: pfrac(nmlonp1,nmlat0) ! NH fraction of potential
    !
    ! Local:
    integer :: j,i
    real(r8),dimension(nmlonp1,nmlat0) :: colatc
    real(r8) :: sinlat,coslat,aslonc,ofdc,cosofc,sinofc,crit1deg

    if (.not. crit_user_set) then
       if (prescribed_period) then
          crit(:) = amie_default_crit(:)*dtr
       else
          crit1deg = max(15._r8,0.5_r8*(theta0(1)+theta0(2))*rtd + 5._r8)
          crit1deg = min(30._r8,crit1deg)
          !
          ! Critical colatitudes:
          crit(1) = crit1deg*dtr
          crit(2) = crit(1) + 15._r8*dtr
       end if
    end if

    if (.not.aurora_params_set) then
       offc(:)   =  1._r8*dtr
       dskofc(:) =  0._r8
    end if

    !
    ! offc(2), dskofc(2) are for northern hemisphere aurora
    !
    ofdc = sqrt(offc(2)**2 + dskofc(2)**2)
    cosofc = cos(ofdc)
    sinofc = sin(ofdc)
    aslonc = asin(dskofc(2)/ofdc)
    !
    ! Define colatc with northern convection circle coordinates
    !
    do j=1,nmlat0
       sinlat = sin(abs(ylatm(j+nmlat0-1)))
       coslat = cos(    ylatm(j+nmlat0-1))
       do i=1,nmlonp1
          colatc(i,j) = cos(ylonm(i)-sunlon+aslonc)
          colatc(i,j) = acos(cosofc*sinlat-sinofc*coslat*colatc(i,j))
       end do ! i=1,nmlonp1
       !
       ! Calculate fractional presence of dynamo equation at each northern
       ! hemisphere geomagnetic grid point. Output in pfrac(nmlonp1,nmlat0)
       !
       do i = 1 , nmlonp1
          pfrac(i,j) = (colatc(i,j)-crit(1)) / (crit(2)-crit(1))
          if (pfrac(i,j) < 0._r8) then
             pfrac(i,j)  = 0._r8
          end if
          if (pfrac(i,j) >= 1._r8) then
             pfrac(i,j) = 1._r8
          end if
       end do ! i=1,nmlonp1
    end do ! j=1,nmlat0
    !
  end subroutine calc_pfrac
  !-----------------------------------------------------------------------
  subroutine calc_xi(o2,o1,n2,rmass,rmassN2, rmassO2, rmassO1, ans,i0,i1,j0,j1,nlev)
    use physconst, only: pi
    use getapex  , only: bmod2
    ! Calculate ratio between collision frequency and gyrofrequency
    !
    ! Inputs:
    integer,intent(in) :: i0,i1,j0,j1,nlev
    real(r8),dimension(nlev,i0:i1,j0:j1),intent(in) :: &
         o2,   &  ! m^3
         o1,&  !
         n2     !
    real(r8),intent(in) :: rmass,rmassN2,rmassO2,rmassO1 ! m/s^2
    !
    ! Output:
    real(r8),intent(out) :: ans(nlev,i0:i1,j0:j1)    ! xi/(1+xi^2) output
    real(r8):: collision(nlev,i0:i1,j0:j1) ! dimensioned for vectorization
    !
    ! Local:
    integer :: i,j,k
    real(r8) :: fixpart
    real(r8) :: gyrofrequency(nlev,i0:i1,j0:j1) ! dimensioned for vectorization

    real(r8),parameter :: na=6.02e26_r8   !kg-mol^-1
    real(r8),parameter :: ele=1.6e-19_r8   !C
    real(r8),parameter :: eps=8.85e-12_r8   !epsilon F/m
    real(r8),parameter :: alpha_n2=1.74e-30_r8   !polarizability of N2 m^3
    real(r8),parameter :: alpha_o2=1.6e-30_r8   !polarizability of O2 m^3
    real(r8),parameter :: alpha_o1=8.02e-31_r8   !polarizability of O m^3

    fixpart=2._r8*pi*ele*sqrt(1/(4._r8*pi*eps))
    do k=1,nlev
       do j=j0,j1
          do i=i0,i1
             collision(k,i,j)=(sqrt(alpha_n2/(rmass*rmassN2/((rmass+rmassN2)*na)))*(rmassN2/(rmass+rmassN2))*n2(k,i,j)+ &
                               sqrt(alpha_o2/(rmass*rmassO2/((rmass+rmassO2)*na)))*(rmassO2/(rmass+rmassO2))*o2(k,i,j)+ &
                               sqrt(alpha_o1/(rmass*rmassO1/((rmass+rmassO1)*na)))*(rmassO1/(rmass+rmassO1))*o1(k,i,j))*fixpart
             gyrofrequency(k,i,j)=ele*bmod2(i,j)*1.0e-4_r8/(rmass/na)   !1e-4 gauss-->T
             ans(k,i,j)=collision(k,i,j)/gyrofrequency(k,i,j)
          enddo
       enddo
    enddo
  end subroutine calc_xi
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine sunloc(iday, secs, sunlon)
    !
    ! Given day of year and ut, return sun's longitude in dipole coordinates
    ! in sunlon
    !
    use getapex,      only: alonm ! (nlonp1,0:nlatp1)
    use edyn_geogrid, only: nlon, nlat, dphi, dlamda
    use edyn_params,  only: pi
    !
    ! Args:
    integer,intent(in)   :: iday          ! day of year
    real(r8),intent(in)  :: secs          ! ut in seconds
    real(r8),intent(out) :: sunlon        ! output
    !
    ! Local:
    integer :: j, i, ii, isun, jsun
    real(r8) :: glats, glons, pisun, pjsun, sndlons, csdlons
    real(r8) :: rlonm(nlon+4, nlat) ! (nlon+4,nlat)
    real(r8) :: r8_isun, r8_jsun

    !
    ! Sun's geographic coordinates:
    glats   = asin(.398749_r8*sin(2._r8 * pi * real(iday-80, r8) / 365._r8))
    glons   = pi * (1._r8 - (2._r8 * secs / 86400._r8))

    do j = 1, nlat
       do i = 1, nlon
          ii = i + 2
          rlonm(ii, j) = alonm(i, j)
       end do
       do i = 1, 2
          rlonm(i, j) = rlonm(i+nlon, j)
          rlonm(i+nlon+2, j) = rlonm(i+2, j)
       end do
    end do

    pisun = ((glons + pi) / dlamda) + 1._r8
    pjsun = ((glats + (.5_r8 * (pi - dphi))) / dphi) + 1._r8
    isun = int(pisun)
    jsun = int(pjsun)
    r8_isun = real(isun, r8)
    r8_jsun = real(jsun, r8)
    pisun = pisun - r8_isun
    pjsun = pjsun - r8_jsun

    sndlons = (1._r8-pisun) * (1._r8-pjsun) * sin(rlonm(isun+2, jsun)) +      &
         pisun*(1._r8-pjsun)       *          sin(rlonm(isun+3,jsun)) +       &
         pisun*pjsun               *          sin(rlonm(isun+3,jsun+1)) +     &
         (1._r8-pisun)*pjsun       *          sin(rlonm(isun+2,jsun+1))
    csdlons = (1._r8-pisun) * (1._r8-pjsun) * cos(rlonm(isun+2,jsun)) +       &
         pisun*(1._r8-pjsun)       *          cos(rlonm(isun+3,jsun))+        &
         pisun*pjsun               *          cos(rlonm(isun+3,jsun+1))+      &
         (1._r8-pisun)*pjsun       *          cos(rlonm(isun+2,jsun+1))
    sunlon = atan2(sndlons, csdlons)

  end subroutine sunloc

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine outfld_phys1d( fldname, array )
    use ppgrid,   only: pcols, begchunk, endchunk
    use phys_grid,only: get_ncols_p

    character(len=*), intent(in) :: fldname
    real(r8), intent(in) :: array(:)

    integer :: i,j, lchnk,ncol
    real(r8) :: tmparr(pcols)

    if (hist_fld_active(fldname)) then
       j = 0
       do lchnk = begchunk, endchunk
          ncol = get_ncols_p(lchnk)
          do i = 1, ncol
             j = j + 1
             tmparr(i) = array(j)
          enddo
          call outfld(fldname,tmparr(:ncol),ncol,lchnk)
       enddo
    end if

  end subroutine outfld_phys1d
  !-----------------------------------------------------------------------
  subroutine outfld_phys( fldname, array )
    use ppgrid,   only: pcols, pver, begchunk, endchunk
    use phys_grid,only: get_ncols_p

    character(len=*), intent(in) :: fldname
    real(r8), intent(in) :: array(:,:)

    integer :: i,j,k, lchnk,ncol
    real(r8) :: tmparr(pcols, pver)

    if (hist_fld_active(fldname)) then
       j = 0
       do lchnk = begchunk, endchunk
          ncol = get_ncols_p(lchnk)
          do i = 1, ncol
             j = j + 1
             do k = 1, pver
                tmparr(i,k) = array(k,j)
             enddo
          enddo
          call outfld(fldname,tmparr(:ncol,:),ncol,lchnk)
       enddo
    end if

  end subroutine outfld_phys
  !-----------------------------------------------------------------------
  subroutine outfld_geokij( name, array, ilev0,ilev1, ilon0,ilon1, ilat0,ilat1 )

    character(len=*), intent(in) :: name
    integer,  intent(in) :: ilev0,ilev1, ilon0,ilon1, ilat0,ilat1
    real(r8), intent(in) :: array(ilev0:ilev1, ilon0:ilon1, ilat0:ilat1)

    integer :: j,k
    real(r8) :: tmpout(ilon0:ilon1,ilev0:ilev1)

    do j = ilat0,ilat1
       do k = ilev0,ilev1
          tmpout(ilon0:ilon1,k) = array(ilev1-k+1,ilon0:ilon1,j)
       end do
       call outfld( name, tmpout, ilon1-ilon0+1, j )
    end do
  end subroutine outfld_geokij
  !-----------------------------------------------------------------------
  subroutine outfld_geo( fldname, array )
    use edyn_mpi, only: lon0, lon1, lat0, lat1, lev0, lev1

    character(len=*), intent(in) :: fldname
    real(r8), intent(in) :: array(lon0:lon1,lat0:lat1,lev0:lev1)

    integer :: j

    do j = lat0,lat1
       call outfld( fldname, array(lon0:lon1,j,lev0:lev1), lon1-lon0+1, j )
    end do

  end subroutine outfld_geo
  !-----------------------------------------------------------------------
  subroutine outfld_geo2d( fldname, array )
    use edyn_mpi, only: lon0, lon1, lat0, lat1

    character(len=*), intent(in) :: fldname
    real(r8), intent(in) :: array(lon0:lon1,lat0:lat1)

    integer :: j

    do j = lat0,lat1
       call outfld( fldname, array(lon0:lon1,j), lon1-lon0+1, j )
    end do

  end subroutine outfld_geo2d
  !-----------------------------------------------------------------------
  subroutine outfld_mag( fldname, array )
    use edyn_mpi, only: omlon1, mlon0, mlon1, mlat0, mlat1, mlev0, mlev1

    character(len=*), intent(in) :: fldname
    real(r8), intent(in) :: array(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)

    integer :: j

    do j = mlat0,mlat1
       call outfld( fldname, array(mlon0:omlon1,j,mlev0:mlev1),omlon1-mlon0+1,j)
    end do

  end subroutine outfld_mag
  !-----------------------------------------------------------------------
  subroutine outfld_mag2d( fldname, array )
    use edyn_mpi, only: mlon0, mlon1, mlat0, mlat1

    character(len=*), intent(in) :: fldname
    real(r8), intent(in) :: array(mlon0:mlon1,mlat0:mlat1)

    integer :: j

    do j = mlat0,mlat1
       call outfld( fldname, array(mlon0:mlon1,j), mlon1-mlon0+1, j )
    end do

  end subroutine outfld_mag2d

end module dpie_coupling
