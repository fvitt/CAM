module micm_mod
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use cam_logfile,     only: iulog
  use spmd_utils,      only: masterproc
  use cam_abortutils,  only: endrun
  use const_props_mod, only: const_props_type
  use ppgrid,          only: pcols, pver
  use json_loader,     only: json_loader_read

  implicit none
  private

  public :: micm_readnl
  public :: micm_register
  public :: micm_initialize
  public :: micm_timestep_tend
  public :: micm_timestep_init

  type(const_props_type), allocatable :: cnst_props(:)
  character(len=16), allocatable :: jnames(:)
  integer, allocatable :: micm_cnst_map(:)

  integer :: nspecies
  integer :: njrxns
  integer :: nkrxns

  integer :: nwave
  integer :: h2o_ndx = -1
  integer :: o2_ndx = -1
  integer :: o3_ndx = -1
  integer :: so2_ndx = -1
  integer :: no2_ndx = -1
  integer :: ice_ndx = -1
  integer :: liq_ndx = -1
  integer :: ndx_cldfr

  real(r8) :: dtime

  character(len=16), allocatable :: knames(:)
  character(len=24), allocatable :: photname(:)
  
contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine micm_readnl(nlfile)
    use tuv_photolysis, only: tuv_photolysis_readnl

    character(len=*), intent(in) :: nlfile

    integer :: errflg
    character(len=1200) :: errmsg

    call tuv_photolysis_readnl('./photolysis_options_in', errmsg, errflg)
    if (errflg/=0) then
       call endrun('micm_readnl: '//trim(errmsg))
    end if
    
  end subroutine micm_readnl
  
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine micm_register

    use constituents, only: cnst_add

    real(r8), parameter :: cptmp = -huge(1._r8) ! does not seem to be used in CAM
    real(r8), parameter :: qmin = 0._r8
    real(r8) :: wght
    character(len=16) :: name
    integer :: m,n

    call initialize_chemistry_information()

    allocate( micm_cnst_map(nspecies) )
    micm_cnst_map(:) = -1
    do m = 1,nspecies
       name = trim(cnst_props(m)%get_name())//'_micm'
       wght = cnst_props(m)%get_wght()
       call cnst_add( name, wght, cptmp, qmin, n, cam_outfld=.false. )
       micm_cnst_map(m) = n
    enddo

  end subroutine micm_register

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine micm_initialize
    use cam_history,            only: addfld, add_default, horiz_only
    use molec_ox_xsect,         only: molec_ox_xsect_init
    use tuv_photolysis,         only: tuv_photolysis_init
    use tuv_radiation_transfer, only: tuv_radiation_transfer_init
    use chemistry_driver,       only: chemistry_driver_init
    use time_manager,           only: get_step_size
    use constituents,           only: cnst_get_ind 
    use physics_buffer,         only: pbuf_get_index

    integer :: errflg
    character(len=1200) :: errmsg

    character(len=16) :: name
    integer :: m
    character(len=128) :: reaction_names(nkrxns+njrxns)
    character(len=3) :: str

    do m = 1,nspecies
       name = trim(cnst_props(m)%get_name())//'_micm'
       call addfld( name, (/ 'lev' /), 'A', 'mole/mole', 'MICM '//trim(cnst_props(m)%get_name()) )
       call add_default( name, 1, ' ' )
       call add_default( name, 10, ' ' )
    end do

    allocate(photname(njrxns))
    do m = 1,njrxns
       write(str,'(i2.2)') m
       photname(m) = trim(jnames(m))//'_micm'//trim(str)
       call addfld( photname(m), (/ 'lev' /), 'A', '1/sec', 'MICM '//trim(jnames(m)) )
       call add_default( photname(m), 10, ' ' )
    end do

    allocate(knames(nkrxns))

    do m = 1, nkrxns
       write(knames(m),fmt='(a1,i3.3,a5)') 'k',m,'_micm'
       call addfld( knames(m), (/ 'lev' /), 'A', '1/sec', 'MICM '//trim(knames(m)) )
       call add_default( knames(m), 10, ' ' )
    end do
    
    call molec_ox_xsect_init( errmsg, errflg )
    if (errflg/=0) then
       call endrun('micm_initialize: '//trim(errmsg))
    end if

    call tuv_radiation_transfer_init( r8, errmsg, errflg )
    if (errflg/=0) then
       call endrun('micm_initialize: '//trim(errmsg))
    end if

    call tuv_photolysis_init( r8, pver, jnames, nwave, errmsg, errflg )
    if (errflg/=0) then
       call endrun('micm_initialize: '//trim(errmsg))
    end if

    dtime = get_step_size()

    call chemistry_driver_init(nspecies, nkrxns, njrxns, reaction_names, TimeStart=0._r8, TimeEnd=dtime, dt=dtime, &
         options_filepath='./solver_options_in', print_log_message=masterproc, errmsg=errmsg, errflg=errflg)
    if (errflg/=0) then
       call endrun('micm_initialize: '//trim(errmsg))
    end if

    call cnst_get_ind ('Q',   h2o_ndx)
    call cnst_get_ind ('O2',   o2_ndx)
    call cnst_get_ind ('O3',   o3_ndx)
    call cnst_get_ind ('SO2', so2_ndx)
    call cnst_get_ind ('NO2', no2_ndx)
    call cnst_get_ind( 'CLDLIQ', liq_ndx )
    call cnst_get_ind( 'CLDICE', ice_ndx )

    ndx_cldfr = pbuf_get_index('CLD')

  end subroutine micm_initialize

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine micm_timestep_init

    use wavelength_grid,  only: tuv_nw=>nwave, tuv_we=>wl
    use module_prates_tuv,only: update_etf
    use solar_irrad_data, only: data_nbins=>nbins, data_we=>we, data_etf=>sol_irrad ! W/m2/nm
    use mo_util,          only: rebin

    real(r8) :: etf(tuv_nw) 

    call rebin( data_nbins, tuv_nw, data_we, tuv_we, data_etf, etf )

    call update_etf(etf)
  end subroutine micm_timestep_init
  
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine micm_timestep_tend( cam_in, state, ptend, pbuf )
    use camsrfexch,    only: cam_in_t
    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use physics_types, only: physics_ptend_sum
    use constituents,  only: cnst_mw, pcnst
    use cam_history,   only: outfld
    use physconst,     only: mbarv, rairv
    use mo_constants,  only: boltz_cgs
    use wv_saturation, only: qsat
    use k_rateConst,   only: k_rateConst_run
    use orbit,         only: zenith
    use shr_const_mod, only: shr_const_pi
    use ref_pres,      only: ptop_ref
    use phys_grid,     only: get_rlat_all_p, get_rlon_all_p
    use time_manager,  only: get_curr_calday, get_step_size
    use molec_ox_xsect,only: molec_ox_xsect_run
    use tuv_radiation_transfer, only: tuv_radiation_transfer_run
    use tuv_photolysis,only: tuv_photolysis_run
    use chemistry_driver, only: chemistry_driver_run
    use physics_buffer,only: physics_buffer_desc, pbuf_get_field

    type(cam_in_t),      intent(in)  :: cam_in
    type(physics_state), intent(in)  :: state           ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend           ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: i,k,m,n
    logical :: lq(pcnst)

    real(r8) :: vmr(state%ncol,pver,nspecies)
    real(r8) :: vmr_out(state%ncol,pver,nspecies)
    real(r8) :: k_rateConsts(pcols,pver,nkrxns)
    real(r8) :: j_rateConsts(pcols,pver,njrxns)
    real(r8) :: reaction_rates(njrxns+nkrxns)
    real(r8) :: reaction_rate_constants(njrxns+nkrxns)
    real(r8) :: sad(4)
    real(r8) :: adiam(4)
    integer :: lchnk, ncol
    character(len=16) :: name

    integer :: errflg
    character(len=1200) :: errmsg
    
    real(r8) :: airdens(pcols,pver)          ! molecules/cm3
    real(r8) :: h2o_vmr(pcols,pver)          ! mole/mole
    real(r8) :: o2_vmr(pcols,pver)           ! mole/mole
    real(r8) :: o3_vmr(pcols,pver)           ! mole/mole
    real(r8) :: so2_vmr(pcols,pver)          ! mole/mole
    real(r8) :: no2_vmr(pcols,pver)          ! mole/mole
    real(r8) :: cld_mmr(pcols,pver)          ! mole/mole
    real(r8) :: satv(pcols,pver)             ! wrk array for relative humidity
    real(r8) :: satq(pcols,pver)             ! wrk array for relative humidity
    real(r8) :: relhum(pcols,pver)           ! relative humidity (fraction)
    real(r8), parameter :: Pa_xfac = 10._r8  ! Pascals to dyne/cm^2
    real(r8) :: zen_angle(pcols), sza(pcols) ! solar zenith angle
    real(r8), parameter :: rad2deg = 180._r8/shr_const_pi ! radians to degrees conversion factor
    real(r8) :: rlats(pcols), rlons(pcols)
    real(r8) :: calday
    real(r8) :: watdens(pver) ! g/m3 <-- cldw (kg/kg)
    real(r8) :: cloudfr(pver) ! 

    real(r8) :: dto2(pver,nwave)
    real(r8) :: srb_o2_xs(nwave,pver)
    real(r8) :: radfld(nwave,pver)
    real(r8), pointer :: cldfr(:,:)
    type(physics_ptend) :: ptend_loc           ! Local parameterization tendencies
   
    lchnk = state%lchnk
    ncol = state%ncol
    
    lq(:) = .false.
    do m = 1,nspecies
       n = micm_cnst_map(m)
       lq(n) = .true.
       do k = 1,pver
          vmr(:ncol,k,m) = mbarv(:ncol,k,lchnk) * state%q(:ncol,k,n) / cnst_mw(n)
       end do
    end do
    call physics_ptend_init(ptend_loc, state%psetcols, 'MICM', lq=lq)

    ptend_loc%q = 0.0_r8

    calday = get_curr_calday( )
    call get_rlat_all_p( lchnk, ncol, rlats(:ncol) )
    call get_rlon_all_p( lchnk, ncol, rlons(:ncol) )
    call zenith( calday, rlats(:ncol), rlons(:ncol), zen_angle(:ncol), ncol )
    sza(:ncol) = acos( zen_angle(:ncol) ) * rad2deg

    !-----------------------------------------------------------------
    !	... compute the relative humidity
    !-----------------------------------------------------------------
    call qsat(state%t(:ncol,:), state%pmid(:ncol,:), satv(:ncol,:), satq(:ncol,:))

    do k = 1,pver
       h2o_vmr(:ncol,k) = mbarv(:ncol,k,lchnk) * state%q(:ncol,k,h2o_ndx) / cnst_mw(h2o_ndx)
       o2_vmr(:ncol,k) = mbarv(:ncol,k,lchnk) * state%q(:ncol,k,o2_ndx) / cnst_mw(o2_ndx)
       o3_vmr(:ncol,k) = mbarv(:ncol,k,lchnk) * state%q(:ncol,k,o3_ndx) / cnst_mw(o3_ndx)
       so2_vmr(:ncol,k) = mbarv(:ncol,k,lchnk) * state%q(:ncol,k,so2_ndx) / cnst_mw(o2_ndx)
       no2_vmr(:ncol,k) = mbarv(:ncol,k,lchnk) * state%q(:ncol,k,no2_ndx) / cnst_mw(o2_ndx)
       airdens(:ncol,k) = Pa_xfac * state%pmid(:ncol,k) / (boltz_cgs*state%t(:ncol,k))
       relhum(:ncol,k) = .622_r8 * h2o_vmr(:ncol,k) / satq(:ncol,k)
       relhum(:ncol,k) = max( 0._r8,min( 1._r8,relhum(:ncol,k) ) )
       cld_mmr(:ncol,k) = state%q(:ncol,k,liq_ndx) + state%q(:ncol,k,ice_ndx)
    end do

    call pbuf_get_field(pbuf, ndx_cldfr, cldfr )

    do i = 1,ncol
       call molec_ox_xsect_run( pver, sza(i), state%zm(i,:), state%t(i,:), state%pmid(i,:), ptop_ref, o2_vmr(i,:), dto2, srb_o2_xs, errmsg, errflg )
       if (errflg/=0) then
          call endrun('tuv_cam_run: '//trim(errmsg))
       end if

       !            1.e3_r8 * kg/kg * Pa / ((J/K/kg) * K) --> g/m3
       watdens(:) = 1.e3_r8 * cld_mmr(i,:)*state%pmid(i,:)/(rairv(i,:,lchnk)*state%t(i,:))
       cloudfr(:) = cldfr(i,:)

       call tuv_radiation_transfer_run( pver, nwave, sza(i), cam_in%asdir(i), state%pmid(i,:), ptop_ref, state%zm(i,:), state%t(i,:), &
                                        o3_vmr(i,:), so2_vmr(i,:), no2_vmr(i,:), cloudfr, watdens, dto2, radfld, errmsg, errflg )
       if (errflg/=0) then
          call endrun('micm_timestep_tend: '//trim(errmsg))
       end if

       call tuv_photolysis_run( pver, state%t(i,:), state%pmid(i,:), radfld, srb_o2_xs, j_rateConsts(i,:,:), errmsg, errflg )
       if (errflg/=0) then
          call endrun('micm_timestep_tend: '//trim(errmsg))
       end if
   end do

   do k = 1,pver
      do i = 1,ncol
         call k_rateConst_run(nkrxns, njrxns, k_rateConsts(i,k,:), airdens(i,k), relhum(i,k), h2o_vmr(i,k), o2_vmr(i,k), &
                              state%t(i,k), state%pmid(i,k),  sad, adiam,  errmsg, errflg)
         if (errflg/=0) then
            call endrun('micm_timestep_tend: '//trim(errmsg))
         end if

         vmr_out(i,k,:nspecies) = vmr(i,k,:nspecies)
         call chemistry_driver_run(vmr_out(i,k,:nspecies), 0._r8, dtime, j_rateConsts(i,k,:), k_rateConsts(i,k,:), airdens(i,k), &
                                   reaction_rates, reaction_rate_constants, errmsg, errflg)
         if (errflg/=0) then
            call endrun('micm_timestep_tend: '//trim(errmsg))
         end if
         do m = 1,nspecies
            n = micm_cnst_map(m)

            ptend_loc%q(i,k,n) =  (vmr_out(i,k,m) - vmr(i,k,m)) * cnst_mw(n) / mbarv(i,k,lchnk) / dtime

         end do

      end do
   end do

   do m = 1, nkrxns
      call outfld( knames(m), k_rateConsts(:ncol,:,m), ncol ,lchnk )
   end do

   do m = 1,njrxns
      call outfld( photname(m), j_rateConsts(:ncol,:,m), ncol ,lchnk )
   end do

   do m = 1,nspecies
      name = trim(cnst_props(m)%get_name())//'_micm'
      call outfld( name, vmr_out(:ncol,:,m), ncol ,lchnk )
   end do

   call physics_ptend_sum(ptend_loc, ptend, ncol)
    
  end subroutine micm_timestep_tend

  ! prviate methods

  subroutine initialize_chemistry_information

    character(len=*), parameter :: jsonfile = '/glade/u/home/fvitt/camdev/micm_trop1/src/chemistry/micm/trop1_kinetics/mechanism.json'
    integer :: i

    call json_loader_read( jsonfile, cnst_props, nspecies, nkrxns, njrxns, jnames )

    if (masterproc) then
       write(iulog,*) 'micm_mod: nspecies = ',nspecies 
       do i=1,nspecies
          write(iulog,'(a,f6.2)') 'micm_mod: '//trim(cnst_props(i)%get_name())//' wght: ',cnst_props(i)%get_wght()
       enddo
       write(iulog,*) 'micm_mod: njrxns = ',njrxns
       do i=1,njrxns
          write(iulog,'(a,i4)') 'micm_mod: '//trim(jnames(i)),i
       enddo
    endif

  end subroutine initialize_chemistry_information

end module micm_mod
