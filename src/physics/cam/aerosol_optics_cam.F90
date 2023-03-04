module aerosol_optics_cam
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_kind_mod, only: cl => shr_kind_cl
  use cam_logfile,  only: iulog
  use radconstants, only: nswbands, nlwbands
  use physics_types,only: physics_state
  use physics_buffer,only: physics_buffer_desc
  use ppgrid, only: pcols, pver
  use physconst, only: rga
  use cam_abortutils, only: endrun
  use spmd_utils, only : masterproc
  use radconstants, only: nswbands, nlwbands
  use wv_saturation, only: qsat

  use aerosol_properties_mod, only: aerosol_properties
  use modal_aerosol_properties_mod, only: modal_aerosol_properties
  use carma_aerosol_properties_mod, only: carma_aerosol_properties

  use aerosol_state_mod,      only: aerosol_state
  use modal_aerosol_state_mod,only: modal_aerosol_state
  use carma_aerosol_state_mod,only: carma_aerosol_state

  use aerosol_optics_mod,     only: aerosol_optics
  use refactive_aerosol_optics_mod, only: refactive_aerosol_optics
  use hygrocoreshell_aerosol_optics_mod, only: hygrocoreshell_aerosol_optics
  use hygrowghtpct_aerosol_optics_mod, only: hygrowghtpct_aerosol_optics

  implicit none

  private

  public :: aerosol_optics_cam_readnl
  public :: aerosol_optics_cam_init
  public :: aerosol_optics_cam_final
  public :: aerosol_optics_cam_sw
  public :: aerosol_optics_cam_lw

  class(aerosol_properties), pointer :: aero_props => null()

  ! refractive index for water read in read_water_refindex
  complex(r8) :: crefwsw(nswbands) = -huge(1._r8) ! complex refractive index for water visible
  complex(r8) :: crefwlw(nlwbands) = -huge(1._r8) ! complex refractive index for water infrared
  character(len=cl) :: aerwat_refindex_file = 'NONE' ! full pathname for water refractive index dataset

  logical :: carma_active = .false.
  logical :: modal_active = .false.

contains

  !===============================================================================
  subroutine aerosol_optics_cam_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils,     only : mpicom, masterprocid, mpi_character, mpi_success

    character(len=*), intent(in)  :: nlfile  ! filepath for file containing namelist input

    integer                       :: unitn, ierr
    character(len=*), parameter   :: subname = 'aerosol_optics_cam_readnl'

    ! ===================
    ! Namelist definition
    ! ===================
    namelist /aerosol_optics_nl/ aerwat_refindex_file

    ! =============
    ! Read namelist
    ! =============
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aerosol_optics_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aerosol_optics_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! ============================
    ! Broadcast namelist variables
    ! ============================
    call mpi_bcast(aerwat_refindex_file, len(aerwat_refindex_file),  mpi_character, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname // ':: ERROR mpi_bcast '//trim(aerwat_refindex_file))
    end if

    if (masterproc) then
       write(iulog,*) subname,': aerwat_refindex_file = ',trim(aerwat_refindex_file)
    end if

  end subroutine aerosol_optics_cam_readnl

  !===============================================================================
  subroutine aerosol_optics_cam_init
    use rad_constituents, only: rad_cnst_get_info

    integer :: nmodes=0, nbins=0

    call rad_cnst_get_info(0, nmodes=nmodes, nbins=nbins)
    carma_active = nbins>0
    modal_active = nmodes>0

    if (modal_active) then
       aero_props => modal_aerosol_properties()
    else if (carma_active) then
       aero_props => carma_aerosol_properties()
    end if

    if (aerwat_refindex_file/='NONE') then
       call read_water_refindex(aerwat_refindex_file)
    end if

  end subroutine aerosol_optics_cam_init

  !===============================================================================
  subroutine aerosol_optics_cam_final
    deallocate(aero_props)
    nullify(aero_props)
  end subroutine aerosol_optics_cam_final

  !===============================================================================
  subroutine aerosol_optics_cam_sw(list_idx, state, pbuf, nnite, idxnite, &
       tauxar, wa, ga, fa)

    ! calculates aerosol sw radiative properties

    integer,             intent(in) :: list_idx       ! index of the climate or a diagnostic list
    type(physics_state), intent(in), target :: state          ! state variables

    type(physics_buffer_desc), pointer :: pbuf(:)
    integer,             intent(in) :: nnite          ! number of night columns
    integer,             intent(in) :: idxnite(nnite) ! local column indices of night columns

    real(r8), intent(inout) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
    real(r8), intent(inout) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
    real(r8), intent(inout) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
    real(r8), intent(inout) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

    character(len=*), parameter :: prefix = 'aerosol_optics_cam_sw: '

    integer :: errcode
    character(len=cl) :: errmsg

    integer :: ibin, nbins
    integer :: iwav, ilev
    integer :: ncol, icol, istat

    class(aerosol_state), pointer :: aero_state
    class(aerosol_optics), pointer :: aero_optics

    real(r8) :: dopaer(pcols)
    real(r8) :: mass(pcols,pver)

    real(r8), allocatable :: pext(:)
    real(r8), allocatable :: palb(:)
    real(r8), allocatable :: pasm(:)

    real(r8) :: relh(pcols,pver)
    real(r8) :: sate(pcols,pver)     ! saturation vapor pressure
    real(r8) :: satq(pcols,pver)     ! saturation specific humidity

    character(len=32) :: opticstype

    nullify(aero_state)
    nullify(aero_optics)

    if (.not.associated(aero_props)) return

    if (modal_active) then
       aero_state => modal_aerosol_state( state, pbuf )
    else if (carma_active) then
       aero_state => carma_aerosol_state( state, pbuf )
    end if

    if (.not.associated(aero_state)) then
       call endrun(prefix//'aero_state not associated')
    end if

    ncol = state%ncol

    mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

    allocate(pext(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pext')
    end if
    allocate(palb(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: palb')
    end if
    allocate(pasm(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pasm')
    end if

    nbins=aero_props%nbins()


    do ibin = 1, nbins

       call aero_props%optics_params(list_idx, ibin, opticstype=opticstype)

       select case (trim(opticstype))
       case('modal') ! refactive method
          aero_optics=>refactive_aerosol_optics(aero_props, aero_state, list_idx, ibin, ncol, pver, nswbands, nlwbands, crefwsw, crefwlw)
       case('hygroscopic_coreshell')
          ! calculate relative humidity for table lookup into rh grid
          call qsat(state%t(:ncol,:), state%pmid(:ncol,:), sate(:ncol,:), satq(:ncol,:), ncol, pver)
          relh(:ncol,:) = state%q(1:ncol,:,1) / satq(:ncol,:)
          relh(:ncol,:) = max(1.e-20_r8,relh(:ncol,:))
          aero_optics=>hygrocoreshell_aerosol_optics(aero_props, aero_state, list_idx, ibin, ncol, pver, relh(:ncol,:))
       case('hygroscopic_wtp')
          aero_optics=>hygrowghtpct_aerosol_optics(aero_props, aero_state, list_idx, ibin, ncol, pver)
       case default
          call endrun(prefix//'optics method not recognized')
       end select

       if (associated(aero_optics)) then
          do iwav = 1, nswbands

             do ilev = 1, pver

                call aero_optics%sw_props(ncol, ilev, iwav, pext, palb, pasm )

                do icol = 1,ncol
                   dopaer(icol) = pext(icol)*mass(icol,ilev)
                   tauxar(icol,ilev,iwav) = tauxar(icol,ilev,iwav) + dopaer(icol)
                   wa(icol,ilev,iwav) = wa(icol,ilev,iwav) + dopaer(icol)*palb(icol)
                   ga(icol,ilev,iwav) = ga(icol,ilev,iwav) + dopaer(icol)*palb(icol)*pasm(icol)
                   fa(icol,ilev,iwav) = fa(icol,ilev,iwav) + dopaer(icol)*palb(icol)*pasm(icol)*pasm(icol)
                end do

             end do
          end do
       else
          call endrun(prefix//'aero_optics object pointer not associated')
       end if

       deallocate(aero_optics)
       nullify(aero_optics)

    end do

    deallocate(pext)
    deallocate(palb)
    deallocate(pasm)

    deallocate(aero_state)
    nullify(aero_state)

  end subroutine aerosol_optics_cam_sw

  !===============================================================================
  subroutine aerosol_optics_cam_lw(list_idx, state, pbuf, tauxar)

    ! calculates aerosol lw radiative properties

    integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
    type(physics_state), intent(in), target :: state    ! state variables

    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

    class(aerosol_state), pointer :: aero_state_obj

    real(r8) :: dopaer(pcols)
    real(r8) :: mass(pcols,pver)

    character(len=*), parameter :: prefix = 'aerosol_optics_cam_lw: '

    integer :: errcode
    character(len=cl) :: errmsg

    integer :: ibin, nbins
    integer :: iwav, ilev
    integer :: ncol, icol, istat

    class(aerosol_state), pointer :: aero_state
    class(aerosol_optics), pointer :: aero_optics

    real(r8), allocatable :: pabs(:)

    real(r8) :: relh(pcols,pver)
    real(r8) :: sate(pcols,pver)     ! saturation vapor pressure
    real(r8) :: satq(pcols,pver)     ! saturation specific humidity

    character(len=32) :: opticstype

    nullify(aero_state)
    nullify(aero_optics)

    if (.not.associated(aero_props)) return

    if (modal_active) then
       aero_state => modal_aerosol_state( state, pbuf )
    else if (carma_active) then
       aero_state => carma_aerosol_state( state, pbuf )
    end if

    if (.not.associated(aero_state)) then
       call endrun(prefix//'aero_state not associated')
    end if

    ncol = state%ncol

    mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

    allocate(pabs(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pabs')
    end if

    nbins=aero_props%nbins()

    do ibin = 1, nbins

       call aero_props%optics_params(list_idx, ibin, opticstype=opticstype)

       select case (trim(opticstype))
       case('modal') ! refactive method
          aero_optics=>refactive_aerosol_optics(aero_props, aero_state, list_idx, ibin, ncol, pver, nswbands, nlwbands, crefwsw, crefwlw)
       case('hygroscopic_coreshell')
          ! calculate relative humidity for table lookup into rh grid
          call qsat(state%t(:ncol,:), state%pmid(:ncol,:), sate(:ncol,:), satq(:ncol,:), ncol, pver)
          relh(:ncol,:) = state%q(1:ncol,:,1) / satq(:ncol,:)
          relh(:ncol,:) = max(1.e-20_r8,relh(:ncol,:))
          aero_optics=>hygrocoreshell_aerosol_optics(aero_props, aero_state, list_idx, ibin, ncol, pver, relh(:ncol,:))
       case('hygroscopic_wtp')
          aero_optics=>hygrowghtpct_aerosol_optics(aero_props, aero_state, list_idx, ibin, ncol, pver)
       case default
          call endrun(prefix//'optics method not recognized')
       end select

       if (associated(aero_optics)) then

          do iwav = 1, nswbands

             do ilev = 1, pver
                call aero_optics%lw_props(ncol, ilev, iwav, pabs )

                do icol = 1, ncol
                   dopaer(icol) = pabs(icol)*mass(icol,ilev)
                   tauxar(icol,ilev,iwav) = tauxar(icol,ilev,iwav) + dopaer(icol)
                end do

             end do

          end do

       else
          call endrun(prefix//'aero_optics object pointer not associated')
       end if

       nullify(aero_optics)

    end do

    deallocate(pabs)

  end subroutine aerosol_optics_cam_lw

  !===============================================================================
  ! Private routines
  !===============================================================================

  subroutine read_water_refindex(infilename)
    use cam_pio_utils, only: cam_pio_openfile
    use pio, only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
                   pio_get_var, PIO_NOWRITE, pio_closefile


    ! read water refractive index file and set module data

    character*(*), intent(in) :: infilename   ! modal optics filename

    ! Local variables

    integer            :: i, ierr
    type(file_desc_t)  :: ncid              ! pio file handle
    integer            :: did               ! dimension ids
    integer            :: dimlen            ! dimension lengths
    type(var_desc_t)   :: vid               ! variable ids
    real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
    real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared
    !----------------------------------------------------------------------------

    ! open file
    call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

    ! inquire dimensions.  Check that file values match parameter values.

    ierr = pio_inq_dimid(ncid, 'lw_band', did)
    ierr = pio_inq_dimlen(ncid, did, dimlen)
    if (dimlen .ne. nlwbands) then
       write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
       call endrun('read_modal_optics: bad lw_band value')
    endif

    ierr = pio_inq_dimid(ncid, 'sw_band', did)
    ierr = pio_inq_dimlen(ncid, did, dimlen)
    if (dimlen .ne. nswbands) then
       write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
       call endrun('read_modal_optics: bad sw_band value')
    endif

    ! read variables
    ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
    ierr = pio_get_var(ncid, vid, refrwsw)

    ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
    ierr = pio_get_var(ncid, vid, refiwsw)

    ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
    ierr = pio_get_var(ncid, vid, refrwlw)

    ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
    ierr = pio_get_var(ncid, vid, refiwlw)

    ! set complex representation of refractive indices as module data
    do i = 1, nswbands
       crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)),kind=r8)
    end do
    do i = 1, nlwbands
       crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)),kind=r8)
    end do

    call pio_closefile(ncid)

  end subroutine read_water_refindex

end module aerosol_optics_cam
