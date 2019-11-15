module radxfr_cam
  use shr_kind_mod,   only : r8 => shr_kind_r8
  use ppgrid,         only : pcols, pver
  use ppgrid,         only : begchunk, endchunk
  use cam_logfile,    only : iulog
  use cam_abortutils, only : endrun
  use spmd_utils,     only : masterproc
  use infnan, only : nan, assignment(=)

  use tuv_radiation_transfer, only: tuv_radiation_transfer_init
  use tuv_radiation_transfer, only: tuv_radiation_transfer_run
  use molec_ox_xsect, only: molec_ox_xsect_init
  use molec_ox_xsect, only: molec_ox_xsect_run
  use wavelength_grid, only: nwave
  use interpolate_data, only: lininterp_init, lininterp, interp_type
  use radconstants, only: nswbands

  implicit none

  logical :: do_radxfr = .false.
  logical :: has_aer_ra_feedback = .true.
  real(r8), protected, allocatable :: actinic_fluxes(:,:,:,:) ! (nwave, pver, pcols, begchunk:endchunk )
  
  integer :: swaertau_idx   = -1
  integer :: swaertauw_idx  = -1
  integer :: swaertauwg_idx = -1
  
  real(r8) :: rrtmg_wavelength(nswbands-1)
  
  type (interp_type) :: interp_wgts
  
  integer, parameter  :: nintegrals = 6
  real(r8), protected, allocatable :: swgt(:,:) ! (nintegrals, nwave)
  character(len=16) :: sname(nintegrals)

  integer, parameter  :: nfluxes = 9
  real(r8), parameter :: fwc(nfluxes) = (/ 0., 121., 180., 210., 240., 308., 400., 600., 999. /)
  integer             :: ituv600 = -1
  integer             :: fwv(nfluxes) = -1
  character(len=16)   :: fname(nfluxes)
  
  
contains
  subroutine radxfr_cam_readnl(nlfile)
    use spmd_utils, only : mpicom, masterprocid, mpi_character, mpi_logical
    use units,      only : getunit, freeunit
    use namelist_utils, only : find_group_name  
    use wavelength_grid, only: wavelength_grid_init

    use params_mod, only: input_data_root

    character(len=*), intent(in) :: nlfile

    character(len=64)  :: radxfr_wavelength_grid_file = 'NONE'
    character(len=265) :: radxfr_input_data_root = 'NONE'
    character(len=512) :: wavelen_grid_filepath
    character(len=512) :: errmsg
    logical            :: radxfr_has_aer_ra_feedback
    integer :: errflg
    integer :: unitn, ierr

    namelist /photo_radxfr_opts/ radxfr_wavelength_grid_file, radxfr_input_data_root, radxfr_has_aer_ra_feedback

    errflg=0
    errmsg=' '

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       write(*,*) 'read : '//trim(nlfile)
       call find_group_name(unitn, 'photo_radxfr_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, photo_radxfr_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun('radxfr_cam_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    call mpi_bcast(radxfr_wavelength_grid_file, len(radxfr_wavelength_grid_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_input_data_root, len(radxfr_input_data_root), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(radxfr_has_aer_ra_feedback, 1, mpi_logical, masterprocid, mpicom, ierr)

    do_radxfr = radxfr_wavelength_grid_file /= 'NONE'

    if (.not.do_radxfr) return

    if (masterproc) then    
       write(iulog,*) 'radxfr_cam_readnl: radxfr_input_data_root = '//trim(radxfr_input_data_root)
       write(iulog,*) 'radxfr_cam_readnl: radxfr_wavelength_grid_file = '//trim(radxfr_wavelength_grid_file)
       write(iulog,*) 'radxfr_cam_readnl: radxfr_has_aer_ra_feedback = ', radxfr_has_aer_ra_feedback
    end if

    input_data_root = trim(radxfr_input_data_root)
    wavelen_grid_filepath = trim(input_data_root)//'/'//trim(radxfr_wavelength_grid_file)

    ! call this here since nwave needs to be known earlier than the init phase
    call wavelength_grid_init(wavelen_grid_filepath, errmsg, errflg)
    if(errflg/=0) then
       call endrun('radxfr_cam_readnl: '//trim(errmsg))
    end if

    has_aer_ra_feedback  = radxfr_has_aer_ra_feedback

  end subroutine radxfr_cam_readnl

  subroutine radxfr_cam_init
    use physics_buffer,  only: pbuf_get_index
    use radconstants,    only: get_sw_spectral_boundaries
    use cam_history,     only: addfld, horiz_only
    use wavelength_grid, only: nwave, wc

    integer :: errflg
    integer :: i, iw, j
    character(len=1200) :: errmsg
    real(r8) :: wvl_low(nswbands)
    real(r8) :: wvl_high(nswbands)
    real(r8) :: wdiff
    
    if (.not.do_radxfr) return

    ! physic buffer fields for aerosol optical properties.
    swaertau_idx   = pbuf_get_index('SWAERTAU')
    swaertauw_idx  = pbuf_get_index('SWAERTAUW')
    swaertauwg_idx = pbuf_get_index('SWAERTAUWG')

    allocate( actinic_fluxes(nwave, pver, pcols, begchunk:endchunk ) )
    actinic_fluxes = nan

    allocate( swgt(nintegrals, nwave ) )
    
    errflg=0
    errmsg=' '

    call molec_ox_xsect_init( errmsg, errflg )
    if (errflg/=0) then
       call endrun('radxfr_cam_init: '//trim(errmsg))
    end if

    call tuv_radiation_transfer_init( r8, errmsg, errflg )
    if (errflg/=0) then
       call endrun('radxfr_cam_init: '//trim(errmsg))
    end if
    
    ! Get the RRTMG wavenumber edges and convert to a wavelength center.
    !
    ! NOTE: Last band is a broadband that overlaps the other bands, so skip it.
    call get_sw_spectral_boundaries(wvl_low, wvl_high, "nm")
    rrtmg_wavelength = (wvl_low(1:nswbands-1) + wvl_high(1:nswbands-1)) / 2._r8

    ! Calculate weights needed to interpolate from the RRTMG wavelengths to the
    ! radxfr wavelengths.
    call lininterp_init(rrtmg_wavelength, nswbands-1, wc, nwave, 1, interp_wgts)

    write(iulog, *) 'radxfr: rrtmg wavelengths: ', rrtmg_wavelength
    write(iulog, *) 'radxfr: tuv wavelengths: ', wc 
    
    ! Find the point closest to 600 nm on the TUV wavelength grid
    fname(1) = 'TUV_AFLXLOW'
    fwv(1) = 1
    
    fname(nfluxes) = 'TUV_AFLXHIGH'
    fwv(nfluxes) = nwave
    
    fname(2) = 'TUV_AFLX121'
    fname(3) = 'TUV_AFLX180'
    fname(4) = 'TUV_AFLX210'
    fname(5) = 'TUV_AFLX240'
    fname(6) = 'TUV_AFLX308'
    fname(7) = 'TUV_AFLX400'
    fname(8) = 'TUV_AFLX600'
    
    do j = 2, nfluxes-1
      wdiff = 1000._r8
      fwv(j) = -1
      do i = 1, nwave
        if (abs(wc(i) - fwc(j)) .lt. wdiff) then
          fwv(j) = i
          wdiff = abs(wc(i) - fwc(j))
        end if
      end do

      write(iulog, *) 'radxfr: closest band to ', fwc(j), ' nm ', fwv(j), wc(fwv(j))
    end do
    
    ituv600 = fwv(8)

       ! Add output fields for the tuv aerosol optics in each of the bands: 300, 400, 600, 999 nm.
     ! Add output fields for the tuv aerosol optics in each of the bands: 300, 400, 600, 999 nm.
    call addfld(trim(fname(1)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, shortest', flag_xyfill=.true.)
    call addfld('TUV_TAULOW', (/ 'lev' /), 'A', '1', 'TUV aerosol optical depth, shortest', flag_xyfill=.true.)
    call addfld('TUV_WLOW',   (/ 'lev' /), 'A', '1', 'TUV single scatter albedo, shortest',flag_xyfill=.true.)
    call addfld('TUV_GLOW',   (/ 'lev' /), 'A', '1', 'TUV asymmetry factor, shortest', flag_xyfill=.true.)
    call addfld(trim(fname(8)),  (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 600 nm', flag_xyfill=.true.)
    call addfld('TUV_TAU600', (/ 'lev' /), 'A', '1', 'TUV aerosol optical depth, 600 nm', flag_xyfill=.true.)
    call addfld('TUV_W600',   (/ 'lev' /), 'A', '1', 'TUV single scatter albedo, 600 nm', flag_xyfill=.true.)
    call addfld('TUV_G600',   (/ 'lev' /), 'A', '1', 'TUV asymmetry factor, 600 nm', flag_xyfill=.true.)
    call addfld(trim(fname(nfluxes)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, longest', flag_xyfill=.true.)
    call addfld('TUV_TAUHIGH',(/ 'lev' /), 'A', '1', 'TUV aerosol optical depth, longest', flag_xyfill=.true.)
    call addfld('TUV_WHIGH',  (/ 'lev' /), 'A', '1', 'TUV single scatter albedo, longest', flag_xyfill=.true.)
    call addfld('TUV_GHIGH',  (/ 'lev' /), 'A', '1', 'TUV asymmetry factor, longest', flag_xyfill=.true.)
    call addfld(trim(fname(2)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 121 nm', flag_xyfill=.true.)
    call addfld(trim(fname(3)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 180 nm', flag_xyfill=.true.)
    call addfld(trim(fname(4)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 210 nm', flag_xyfill=.true.)
    call addfld(trim(fname(5)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 240 nm', flag_xyfill=.true.)
    call addfld(trim(fname(6)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 308 nm', flag_xyfill=.true.)
    call addfld(trim(fname(7)), (/ 'lev' /), 'A', 'photons cm-2 s-1 nm-1', 'TUV actinic flux, 400 nm', flag_xyfill=.true.)
    
    ! Create the weights for the spectral integrals and add the fields

    ! UV-A (315-400 nm)
    j=1
    sname(j) = 'TUV_UVA'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'UV-A, 315-400 nm', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 315._r8 .and. wc(iw) .lt. 400._r8) then
          swgt(j,iw) = 1._r8
       else
          swgt(j,iw) = 0._r8
       end if
    end do

    ! UV-B (280-330 nm)
    j=2
    sname(j) = 'TUV_UVB'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'UV-B (280-315 nm)', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 280._r8 .and. wc(iw) .lt. 315._r8) then
          swgt(j,iw) = 1._r8
       else
          swgt(j,iw) = 0._r8
       end if
    end do

    ! UV-B* (280-315 nm)
    j=3
    sname(j) = 'TUV_UVBSTAR'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'UV-B* (280-320 nm)', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 280._r8 .and. wc(iw) .lt. 320._r8) then
          swgt(j,iw) = 1._r8
       else
          swgt(j,iw) = 0._r8
       end if
    end do

    ! UV index (Canadian - WMO/WHO)
    ! Report of the WMO Meeting of experts on UV-B measurements, data quality 
    ! and standardization of UV indices, World Meteorological Organization
    ! (WMO), report No. 95, Geneva, 1994.
    ! based on the CIE erythema weighting, multiplied by 40.
    j=4
    sname(j) = 'TUV_PAR'
    call addfld(trim(sname(j)), horiz_only, 'A', 'umol m-2 s-1', 'Photosynthetic Active Radiation (400 < PAR < 700 nm)', flag_xyfill=.true.)
    do iw = 1, nwave
       if (wc(iw) .gt. 400._r8 .and. wc(iw) .lt. 700._r8) then
          swgt(j,iw) = 8.36e-3_r8 * wc(iw)
       else
          swgt(j,iw) = 0._r8
       end if
    end do

    ! UV index (Canadian - WMO/WHO)
    ! Report of the WMO Meeting of experts on UV-B measurements, data quality 
    ! and standardization of UV indices, World Meteorological Organization
    ! (WMO), report No. 95, Geneva, 1994.
    ! based on the CIE erythema weighting, multiplied by 40.
    j=5
    sname(j) = 'TUV_UVINDEX'
    call addfld(trim(sname(j)), horiz_only, 'A', 'W m-2', 'UV index (Canadian - WMO/WHO)', flag_xyfill=.true.)
    do iw = 1, nwave
       swgt(j,iw) = 40._r8 * fery(wc(iw))
    end do

    ! UV index (Canadian - WMO/WHO)
    ! Report of the WMO Meeting of experts on UV-B measurements, data quality 
    ! and standardization of UV indices, World Meteorological Organization
    ! (WMO), report No. 95, Geneva, 1994.
    ! based on the CIE erythema weighting, multiplied by 40.
    j=6
    sname(j) = 'TUV_UVINDEXMX'
    call addfld(trim(sname(j)), horiz_only, 'X', '1', 'UV index (Canadian - WMO/WHO), Maximum', flag_xyfill=.true.)
    do iw = 1, nwave
       swgt(j,iw) = 40._r8 * fery(wc(iw))
    end do

  end subroutine radxfr_cam_init

  subroutine radxfr_cam_update(ncol, lchnk, esfact, zenith, albedo, press_mid, alt, temp, o2vmr, o3vmr, so2vmr, no2vmr, novmr, cldfrac, cldw, pbuf)
    use physconst,       only : rairv
    use ref_pres,        only : press_top=>ptop_ref
    use physics_buffer,  only : pbuf_get_field, physics_buffer_desc
    use cam_history,     only : outfld
    use wavelength_grid, only : nwave, wl

    integer,  intent(in) :: ncol, lchnk
    real(r8), intent(in) :: esfact
    real(r8), intent(in) :: zenith(:)
    real(r8), intent(in) :: albedo(:)
    real(r8), intent(in) :: press_mid(:,:)
    real(r8), intent(in) :: alt(:,:) !  kilometers
    real(r8), intent(in) :: temp(:,:)
    real(r8), intent(in) :: o2vmr(:,:)
    real(r8), intent(in) :: o3vmr(:,:)
    real(r8), intent(in) :: so2vmr(:,:)
    real(r8), intent(in) :: no2vmr(:,:)
    real(r8), intent(in) :: novmr(:,:)
    real(r8), intent(in) :: cldfrac(:,:)
    real(r8), intent(in) :: cldw(:,:) ! kg/kg
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: i, j, k, iwv
    integer :: errflg
    character(len=512) :: errmsg

    real(r8) :: alt_meters(pver)
    real(r8) :: watdens(pver) ! g/m3 <-- cldw (kg/kg)
    real(r8) :: cloudfr(pver) ! 

    real(r8) :: dto2(pver,nwave)
    real(r8) :: srb_o2_xs(nwave,pver)
    
    real(r8), pointer, dimension(:,:,:) :: swaertau   ! shortwave aerosol tau
    real(r8), pointer, dimension(:,:,:) :: swaertauw  ! shortwave aerosol tau * w
    real(r8), pointer, dimension(:,:,:) :: swaertauwg ! shortwave aerosol tau * w * g

    real(r8) :: swaerw(pcols, pver, nswbands)
    real(r8) :: swaerg(pcols, pver, nswbands)

    real(r8) :: tauaer(pcols, pver, nwave)
    real(r8) :: waer(pcols, pver, nwave)
    real(r8) :: gaer(pcols, pver, nwave)
    
    real(r8) :: spectral_irradiance(pver, nwave, pcols)
    real(r8) :: sintegral(nintegrals, pcols)
    real(r8) :: tmp(pcols, pver)
    
    if (.not.do_radxfr) return

    errflg=0
    errmsg=' '
    
    ! Get the aerosol optical properties.
    if (has_aer_ra_feedback) then
      call pbuf_get_field(pbuf, swaertau_idx,   swaertau)
      call pbuf_get_field(pbuf, swaertauw_idx,  swaertauw)
      call pbuf_get_field(pbuf, swaertauwg_idx, swaertauwg)
    
      ! Need to convert tau*w to w and tau*w*g to g for the radiation code.
      where(swaertau .ne. 0._r8)
        swaerw = swaertauw / swaertau
      elsewhere
        swaerw = 1._r8
      end where

      where(swaertauw .ne. 0._r8)
        swaerg = swaertauwg / swaertauw
      elsewhere
        swaerg = 0._r8
      end where
    
      ! The CESM wavelengths to the wavelength grid used by TUV.
      do i = 1, ncol
        do k = 1, pver
          call lininterp(swaertau(i,k,1:nswbands-1), nswbands-1, tauaer(i,k,:), nwave, interp_wgts)
          call lininterp(swaerw(i,k,1:nswbands-1), nswbands-1, waer(i,k,:), nwave, interp_wgts)
          call lininterp(swaerg(i,k,1:nswbands-1), nswbands-1, gaer(i,k,:), nwave, interp_wgts)
        end do
      end do
    else
      tauaer(:,:,:) = 0._r8
      waer(:,:,:) = 1._r8
      gaer(:,:,:)   = 0._r8
    end if
    
    ! For DEBUG, output the aerosol optical properties for the tuv bands.
    call outfld('TUV_TAULOW', tauaer(:ncol,:,1), ncol, lchnk)
    call outfld('TUV_WLOW',   waer(:ncol,:,1), ncol, lchnk)
    call outfld('TUV_GLOW',   gaer(:ncol,:,1), ncol, lchnk)
    call outfld('TUV_TAU600', tauaer(:ncol,:,ituv600), ncol, lchnk)
    call outfld('TUV_W600',   waer(:ncol,:,ituv600), ncol, lchnk)
    call outfld('TUV_G600',   gaer(:ncol,:,ituv600), ncol, lchnk)
    call outfld('TUV_TAUHIGH', tauaer(:ncol,:,nwave), ncol, lchnk)
    call outfld('TUV_WHIGH',   waer(:ncol,:,nwave), ncol, lchnk)
    call outfld('TUV_GHIGH',   gaer(:ncol,:,nwave), ncol, lchnk)    
    
    do i = 1,ncol    

       alt_meters(:) = alt(i,:)*1.e3_r8 ! km --> m

       call molec_ox_xsect_run( pver, zenith(i), alt_meters, temp(i,:), press_mid(i,:), press_top, o2vmr(i,:), dto2, srb_o2_xs, errmsg, errflg )
       if (errflg/=0) then
          call endrun('radxfr_cam_update: '//trim(errmsg))
       end if

       !            1.e3_r8 * kg/kg * Pa / ((J/K/kg) * K) --> g/m3
       watdens(:) = 1.e3_r8 * cldw(i,:)*press_mid(i,:)/(rairv(i,:,lchnk)*temp(i,:)) 
       cloudfr(:) = cldfrac(i,:)

       call tuv_radiation_transfer_run( pver, nwave, &
            zenith(i), albedo(i), press_mid(i,:), press_top, alt_meters, temp(i,:), &
            o3vmr(i,:), so2vmr(i,:), no2vmr(i,:), cloudfr, watdens, dto2, &
            has_aer_ra_feedback, tauaer(i,:,:), waer(i,:,:), gaer(i,:,:), &
            actinic_fluxes(:,:,i,lchnk), spectral_irradiance(:,:,i), errmsg, errflg )
       if (errflg/=0) then
          call endrun('radxfr_cam_update: '//trim(errmsg))
       end if
       
               
      ! Check for small negative values. Actinic flux should be positive.
      do k = 1, pver
        do iwv = 1, nwave
          if (actinic_fluxes(iwv, k, i, lchnk) .lt. 0._r8) then
            actinic_fluxes(iwv, k, i, lchnk) = 0._r8
          end if

          if (spectral_irradiance(k, iwv, i) .lt. 0._r8) then
            spectral_irradiance(k, iwv, i) = 0._r8
          end if          
        end do
      end do
    end do
    
    ! Only need this at the surface for now.
    !
    ! NOTE: The actinic flux and the spectral irradiance are defined with k=1 as the bottom,
    ! not the top as is normal in CAM.
    sintegral(:,:) = 0._r8
    do j = 1, nintegrals
      do iwv = 1, nwave-1
        sintegral(j, :ncol) = sintegral(j, :ncol) + &
          esfact * swgt(j, iwv) * spectral_irradiance(1, iwv, :ncol) * abs(wl(iwv+1) - wl(iwv))
      end do
      
      ! Output the diagnostic.
      call outfld(trim(sname(j)), sintegral(j, :ncol), ncol, lchnk)
    end do
    
    ! output actinic flux at certain wavelengths import to stratospheric chemistry
    !
    ! NOTE: Need to flip the fluxes for output and also reverse in the vertical.
    do j = 1, nfluxes
      do i = 1, ncol
        do k = 1, pver
          tmp(i,k) = esfact * actinic_fluxes(fwv(j),pver-k+1,i,lchnk)
        end do
      end do
      call outfld(trim(fname(j)), tmp(:ncol, :), ncol, lchnk)
    end do   

  end subroutine radxfr_cam_update

!---------------------------------------------------------------------------
!  PURPOSE:                                                                
!  Calculate the action spectrum value for erythema at a given wavelength   
!  according to: McKinlay, A.F and B.L.Diffey, A reference action spectrum  
!  for ultraviolet induced erythema in human skin, CIE Journal, vol 6,      
!  pp 17-22, 1987.                                                          
!  Value at 300 nm = 0.6486                                                
!---------------------------------------------------------------------------
!  PARAMETERS:                                                              
!  W - REAL, wavelength (nm)                                             
!---------------------------------------------------------------------------

  function fery(w)

    implicit none

    ! input:
    real(r8), intent(in) :: w 

    ! function value:
    real(r8) :: fery

    if (w .lt. 250._r8) then
       fery = 1._r8
       ! outside the ery spectrum range
    elseif ((w .ge. 250._r8) .and. (w .lt. 298._r8)) then
       fery = 1._r8
    elseif ((w .ge. 298._r8) .and. (w .lt. 328._r8)) then
       fery = 10._r8**( 0.094_r8*(298._r8 - w) )
    elseif ((w .ge. 328._r8) .and. (w .lt. 400._r8)) then
       fery = 10._r8**( 0.015_r8*(139._r8 - w) )
    else
       fery = 1.e-36_r8
       ! outside the ery spectrum range
    endif

  end function fery

end module radxfr_cam
