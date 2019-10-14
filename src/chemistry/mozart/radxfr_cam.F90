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

  implicit none

  logical :: do_radxfr = .false.
  real(r8), protected, allocatable :: actinic_fluxes(:,:,:,:) ! (nwave, pver, pcols, begchunk:endchunk )
  
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
    integer :: errflg
    integer :: unitn, ierr

    namelist /photo_radxfr_opts/ radxfr_wavelength_grid_file, radxfr_input_data_root

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

    do_radxfr = radxfr_wavelength_grid_file /= 'NONE'

    if (.not.do_radxfr) return

    if (masterproc) then    
       write(iulog,*) 'radxfr_cam_readnl: radxfr_input_data_root = '//trim(radxfr_input_data_root)
       write(iulog,*) 'radxfr_cam_readnl: radxfr_wavelength_grid_file = '//trim(radxfr_wavelength_grid_file)
    end if

    input_data_root = trim(radxfr_input_data_root)
    wavelen_grid_filepath = trim(input_data_root)//'/'//trim(radxfr_wavelength_grid_file)

    ! call this here since nwave needs to be known earlier than the init phase
    call wavelength_grid_init(wavelen_grid_filepath, errmsg, errflg)
    if(errflg/=0) then
       call endrun('radxfr_cam_readnl: '//trim(errmsg))
    end if

  end subroutine radxfr_cam_readnl

  subroutine radxfr_cam_init
    integer :: errflg
    character(len=1200) :: errmsg
    
    if (.not.do_radxfr) return

    allocate( actinic_fluxes(nwave, pver, pcols, begchunk:endchunk ) )
    actinic_fluxes = nan
    
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

  end subroutine radxfr_cam_init

  subroutine radxfr_cam_update(ncol, lchnk, zenith, albedo, press_mid, alt, temp, o2vmr, o3vmr, so2vmr, no2vmr, novmr, cldfrac, cldw)
    use physconst,       only : rairv
    use ref_pres, only: press_top=>ptop_ref

    integer,  intent(in) :: ncol, lchnk
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

    integer :: i
    integer :: errflg
    character(len=512) :: errmsg

    real(r8) :: alt_meters(pver)
    real(r8) :: watdens(pver) ! g/m3 <-- cldw (kg/kg)
    real(r8) :: cloudfr(pver) ! 

    real(r8) :: dto2(pver,nwave)
    real(r8) :: srb_o2_xs(nwave,pver)

    if (.not.do_radxfr) return

    errflg=0
    errmsg=' '

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
            actinic_fluxes(:,:,i,lchnk) , errmsg, errflg )
       if (errflg/=0) then
          call endrun('radxfr_cam_update: '//trim(errmsg))
       end if

    end do

  end subroutine radxfr_cam_update
  
  end module radxfr_cam
