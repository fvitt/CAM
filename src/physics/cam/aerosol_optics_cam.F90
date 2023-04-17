module aerosol_optics_cam
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_kind_mod, only: cl => shr_kind_cl
  use cam_logfile,  only: iulog
  use radconstants, only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
  use radconstants, only: ot_length, get_lw_spectral_boundaries
  use physics_types,only: physics_state
  use physics_buffer,only: physics_buffer_desc
  use ppgrid, only: pcols, pver
  use physconst, only: rga, rair
  use cam_abortutils, only: endrun
  use spmd_utils, only : masterproc
  use wv_saturation, only: qsat
  use rad_constituents,  only: n_diag, rad_cnst_get_call_list
  use cam_history,       only: addfld, add_default, outfld, horiz_only
  use cam_history_support, only: fillvalue

  use aerosol_properties_mod, only: aerosol_properties
  use modal_aerosol_properties_mod, only: modal_aerosol_properties
  use carma_aerosol_properties_mod, only: carma_aerosol_properties

  use aerosol_state_mod,      only: aerosol_state
  use modal_aerosol_state_mod,only: modal_aerosol_state
  use carma_aerosol_state_mod,only: carma_aerosol_state

  use aerosol_optics_mod,     only: aerosol_optics
  use refractive_aerosol_optics_mod, only: refractive_aerosol_optics
  use hygrocoreshell_aerosol_optics_mod, only: hygrocoreshell_aerosol_optics
  use hygrowghtpct_aerosol_optics_mod, only: hygrowghtpct_aerosol_optics

  implicit none

  private

  public :: aerosol_optics_cam_readnl
  public :: aerosol_optics_cam_init
  public :: aerosol_optics_cam_final
  public :: aerosol_optics_cam_sw
  public :: aerosol_optics_cam_lw

  type aero_props_t
     class(aerosol_properties), pointer :: obj => null()
  end type aero_props_t
  type aero_state_t
     class(aerosol_state), pointer :: obj => null()
  end type aero_state_t

  type(aero_props_t), allocatable :: aero_props(:)

  ! refractive index for water read in read_water_refindex
  complex(r8) :: crefwsw(nswbands) = -huge(1._r8) ! complex refractive index for water visible
  complex(r8) :: crefwlw(nlwbands) = -huge(1._r8) ! complex refractive index for water infrared
  character(len=cl) :: aerwat_refindex_file = 'NONE' ! full pathname for water refractive index dataset

  logical :: carma_active = .false.
  logical :: modal_active = .false.
  integer :: num_aero_models = 0
  integer :: lw10um_indx = -1

  character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

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

    character(len=*), parameter :: prefix = 'aerosol_optics_cam_sw: '
    integer :: nmodes=0, nbins=0, iaermod, istat, ilist, i

    logical :: call_list(0:n_diag)
    real(r8) :: lwavlen_lo(nlwbands), lwavlen_hi(nlwbands)

    num_aero_models = 0

    call rad_cnst_get_info(0, nmodes=nmodes, nbins=nbins)
    modal_active = nmodes>0
    carma_active = nbins>0

    if (modal_active) then
       num_aero_models = num_aero_models+1
    end if
    if (carma_active) then
       num_aero_models = num_aero_models+1
    end if

    if (num_aero_models>0) then
       allocate(aero_props(num_aero_models), stat=istat)
       if (istat/=0) then
          call endrun(prefix//'array allocation error: aero_props')
       end if
    end if

    iaermod = 0

    if (modal_active) then
       iaermod = iaermod+1
       aero_props(iaermod)%obj => modal_aerosol_properties()
    else if (carma_active) then
       iaermod = iaermod+1
       aero_props(iaermod)%obj => carma_aerosol_properties()
    end if

    if (aerwat_refindex_file/='NONE') then
       call read_water_refindex(aerwat_refindex_file)
    end if

    call get_lw_spectral_boundaries(lwavlen_lo, lwavlen_hi, units='um')
    do i = 1,nlwbands
       if ((lwavlen_lo(i)<=10._r8) .and. (lwavlen_hi(i)>=10._r8)) then
          lw10um_indx = i
       end if
    end do
    call rad_cnst_get_call_list(call_list)

    do ilist = 0, n_diag
       if (call_list(ilist)) then
          call addfld ('EXTINCTtest'//diag(ilist),    (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 550 nm, day only', flag_xyfill=.true.)
          call addfld ('EXTINCTUVtest'//diag(ilist),  (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 350 nm, day only', flag_xyfill=.true.)
          call addfld ('EXTINCTNIRtest'//diag(ilist), (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 1020 nm, day only', flag_xyfill=.true.)
          call addfld ('ABSORBtest'//diag(ilist),     (/ 'lev' /), 'A','/m',&
               'Aerosol absorption, day only', flag_xyfill=.true.)
          call addfld ('AODVIStest'//diag(ilist),   horiz_only,  'A','  ', &
               'Aerosol optical depth 550 nm', flag_xyfill=.true.)
          call addfld ('AODUVtest'//diag(ilist),      horiz_only,  'A','  ', &
               'Aerosol optical depth 350 nm, day only', flag_xyfill=.true.)
          call addfld ('AODNIRtest'//diag(ilist),     horiz_only,  'A','  ', &
               'Aerosol optical depth 1020 nm, day only',flag_xyfill=.true.)
          call addfld ('AODABStest'//diag(ilist),     horiz_only,  'A','  ', &
               'Aerosol absorption optical depth 550 nm, day only', flag_xyfill=.true.)
          call addfld ('AODxASYMtest'//diag(ilist),   horiz_only,  'A','  ', &
               'Aerosol optical depth 550 * asymmetry factor, day only', flag_xyfill=.true.)
          call addfld ('EXTxASYMtest'//diag(ilist),   (/ 'lev' /), 'A','  ', &
               'extinction 550 nm * asymmetry factor, day only',  flag_xyfill=.true.)
          call addfld ('AODTOTtest'//diag(ilist), horiz_only, 'A','1',&
               'Aerosol optical depth summed over all sw wavelenghts', flag_xyfill=.true.)

          call addfld ('EXTINCTdntest'//diag(ilist),    (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 550 nm, day only')
          call addfld ('EXTINCTUVdntest'//diag(ilist),  (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 350 nm, day only')
          call addfld ('EXTINCTNIRdntest'//diag(ilist), (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 1020 nm, day only')
          call addfld ('ABSORBdntest'//diag(ilist),     (/ 'lev' /), 'A','/m',&
               'Aerosol absorption, day only')
          call addfld ('AODVISdntest'//diag(ilist),   horiz_only,  'A','  ', &
               'Aerosol optical depth 550 nm')
          call addfld ('AODUVdntest'//diag(ilist),      horiz_only,  'A','  ', &
               'Aerosol optical depth 350 nm, day only')
          call addfld ('AODNIRdntest'//diag(ilist),     horiz_only,  'A','  ', &
               'Aerosol optical depth 1020 nm, day only',flag_xyfill=.true.)
          call addfld ('AODABSdntest'//diag(ilist),     horiz_only,  'A','  ', &
               'Aerosol absorption optical depth 550 nm, day only')
          call addfld ('AODxASYMdntest'//diag(ilist),   horiz_only,  'A','  ', &
               'Aerosol optical depth 550 * asymmetry factor, day only')
          call addfld ('EXTxASYMdntest'//diag(ilist),   (/ 'lev' /), 'A','  ', &
               'extinction 550 nm * asymmetry factor, day only',  flag_xyfill=.true.)
          call addfld ('AODTOTdntest'//diag(ilist), horiz_only, 'A','1',&
               'Aerosol optical depth summed over all sw wavelenghts')

          if (lw10um_indx>0) then
             call addfld('AODABSLWtest'//diag(ilist), (/ 'lev' /), 'A','/m',&
                  'Aerosol long-wave absorption optical depth at 10 microns')
          end if
          call addfld ('TOTABSLWtest'//diag(ilist), (/ 'lev' /), 'A',' ', &
               'LW Aero total abs')
       end if
    end do

  end subroutine aerosol_optics_cam_init

  !===============================================================================
  subroutine aerosol_optics_cam_final

    integer :: iaermod

    do iaermod = 1,num_aero_models
       deallocate(aero_props(iaermod)%obj)
       nullify(aero_props(iaermod)%obj)
    end do

    if (allocated(aero_props)) then
       deallocate(aero_props)
    endif

  end subroutine aerosol_optics_cam_final

  !===============================================================================
  subroutine aerosol_optics_cam_sw(list_idx, state, pbuf, nnite, idxnite, tauxar, wa, ga, fa)

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

    integer :: ibin, nbins
    integer :: iwav, ilev
    integer :: icol, istat
    integer :: lchnk, ncol

    type(aero_state_t), allocatable :: aero_state(:)

    class(aerosol_optics), pointer :: aero_optics

    real(r8) :: dopaer(pcols)
    real(r8) :: mass(pcols,pver)
    real(r8) :: air_density(pcols,pver)

    real(r8), allocatable :: pext(:)
    real(r8), allocatable :: pabs(:)
    real(r8), allocatable :: palb(:)
    real(r8), allocatable :: pasm(:)

    real(r8) :: relh(pcols,pver)
    real(r8) :: sate(pcols,pver)     ! saturation vapor pressure
    real(r8) :: satq(pcols,pver)     ! saturation specific humidity

    character(len=ot_length) :: opticstype
    integer :: iaermod

    real(r8) :: aodvis(pcols)              ! extinction optical depth in vis
    real(r8) :: aoduv(pcols)               ! extinction optical depth in uv
    real(r8) :: aodnir(pcols)              ! extinction optical depth in nir
    real(r8) :: absorb(pcols,pver)
    real(r8) :: aodabs(pcols)              ! absorption optical depth

    real(r8) :: aodtot(pcols)

    real(r8) :: extinct(pcols,pver)
    real(r8) :: extinctnir(pcols,pver)
    real(r8) :: extinctuv(pcols,pver)

    real(r8) :: asymvis(pcols)              ! asymmetry factor * optical depth
    real(r8) :: asymext(pcols,pver)         ! asymmetry factor * extinction

    lchnk = state%lchnk
    ncol  = state%ncol

    nullify(aero_optics)

    mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
    air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

    aodvis = 0._r8
    aodnir = 0._r8
    aoduv = 0._r8
    aodabs = 0._r8
    absorb = 0._r8
    aodtot = 0._r8
    tauxar = 0._r8
    extinct = 0._r8
    extinctnir = 0._r8
    extinctuv = 0._r8
    asymvis = 0.0_r8
    asymext = 0.0_r8

    if (num_aero_models<1) return

    allocate(aero_state(num_aero_models), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: aero_state')
    end if

    iaermod = 0
    if (modal_active) then
       iaermod = iaermod+1
       aero_state(iaermod)%obj => modal_aerosol_state( state, pbuf )
    else if (carma_active) then
       iaermod = iaermod+1
       aero_state(iaermod)%obj => carma_aerosol_state( state, pbuf )
    end if

    allocate(pext(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pext')
    end if
    allocate(pabs(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pabs')
    end if
    allocate(palb(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: palb')
    end if
    allocate(pasm(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pasm')
    end if

    do iaermod = 1,num_aero_models

       nbins=aero_props(iaermod)%obj%nbins()

       do ibin = 1, nbins

          call aero_props(iaermod)%obj%optics_params(list_idx, ibin, opticstype=opticstype)

          select case (trim(opticstype))
          case('modal') ! refractive method
             aero_optics=>refractive_aerosol_optics(aero_props(iaermod)%obj, aero_state(iaermod)%obj, list_idx, ibin, &
                                                    ncol, pver, nswbands, nlwbands, crefwsw, crefwlw)
          case('hygroscopic_coreshell')
             ! calculate relative humidity for table lookup into rh grid
             call qsat(state%t(:ncol,:), state%pmid(:ncol,:), sate(:ncol,:), satq(:ncol,:), ncol, pver)
             relh(:ncol,:) = state%q(1:ncol,:,1) / satq(:ncol,:)
             relh(:ncol,:) = max(1.e-20_r8,relh(:ncol,:))
             aero_optics=>hygrocoreshell_aerosol_optics(aero_props(iaermod)%obj, aero_state(iaermod)%obj, list_idx, &
                                                        ibin, ncol, pver, relh(:ncol,:))
          case('hygroscopic_wtp')
             aero_optics=>hygrowghtpct_aerosol_optics(aero_props(iaermod)%obj, aero_state(iaermod)%obj, list_idx, &
                                                      ibin, ncol, pver)
          case default
             call endrun(prefix//'optics method not recognized')
          end select

          if (associated(aero_optics)) then
             do iwav = 1, nswbands

                do ilev = 1, pver

                   call aero_optics%sw_props(ncol, ilev, iwav, pext, pabs, palb, pasm )

                   do icol = 1,ncol
                      dopaer(icol) = pext(icol)*mass(icol,ilev)
                      tauxar(icol,ilev,iwav) = tauxar(icol,ilev,iwav) + dopaer(icol)
                      wa(icol,ilev,iwav) = wa(icol,ilev,iwav) + dopaer(icol)*palb(icol)
                      ga(icol,ilev,iwav) = ga(icol,ilev,iwav) + dopaer(icol)*palb(icol)*pasm(icol)
                      fa(icol,ilev,iwav) = fa(icol,ilev,iwav) + dopaer(icol)*palb(icol)*pasm(icol)*pasm(icol)

                      if (iwav==idx_uv_diag) then
                         aoduv(icol) = aoduv(icol) + dopaer(icol)
                         extinctuv(icol,ilev) = extinctuv(icol,ilev) + dopaer(icol)*air_density(icol,ilev)/mass(icol,ilev)
                      else if (iwav==idx_sw_diag) then ! vis
                         aodvis(icol) = aodvis(icol) + dopaer(icol)
                         aodabs(icol) = aodabs(icol) + pabs(icol)*mass(icol,ilev)
                         extinct(icol,ilev) = extinct(icol,ilev) + dopaer(icol)*air_density(icol,ilev)/mass(icol,ilev)
                         absorb(icol,ilev)  = absorb(icol,ilev) + pabs(icol)*air_density(icol,ilev)
                         asymvis(icol)      = asymvis(icol) + dopaer(icol)*pasm(icol)
                         asymext(icol,ilev) = asymext(icol,ilev) + dopaer(icol)*pasm(icol)*air_density(icol,ilev)/mass(icol,ilev)
                      else if (iwav==idx_nir_diag) then
                         aodnir(icol) = aodnir(icol) + dopaer(icol)
                         extinctnir(icol,ilev) = extinctnir(icol,ilev) + dopaer(icol)*air_density(icol,ilev)/mass(icol,ilev)
                      end if

                      aodtot(icol) = aodtot(icol) + dopaer(icol)

                   end do

                end do
             end do
          else
             call endrun(prefix//'aero_optics object pointer not associated')
          end if

          deallocate(aero_optics)
          nullify(aero_optics)

       end do
    end do

    call outfld('AODUVdntest'//diag(list_idx),  aoduv,  pcols, lchnk)
    call outfld('AODVISdntest'//diag(list_idx), aodvis, pcols, lchnk)
    call outfld('AODNIRdntest'//diag(list_idx), aodnir, pcols, lchnk)
    call outfld('AODABSdntest'//diag(list_idx), aodabs, pcols, lchnk)
    call outfld('AODTOTdntest'//diag(list_idx), aodtot, pcols, lchnk)
    call outfld('EXTINCTUVdntest'//diag(list_idx),  extinctuv,  pcols, lchnk)
    call outfld('EXTINCTNIRdntest'//diag(list_idx), extinctnir, pcols, lchnk)
    call outfld('EXTINCTdntest'//diag(list_idx),  extinct,  pcols, lchnk)
    call outfld('ABSORBdntest'//diag(list_idx),   absorb,  pcols, lchnk)
    call outfld('EXTxASYMdntest'//diag(list_idx), asymext, pcols, lchnk)
    call outfld('AODxASYMdntest'//diag(list_idx), asymvis, pcols, lchnk)

    do icol = 1, nnite
       aodvis(idxnite(icol)) = fillvalue
       aodnir(idxnite(icol)) = fillvalue
       aoduv(idxnite(icol)) = fillvalue
       aodabs(idxnite(icol)) = fillvalue
       aodtot(idxnite(icol)) = fillvalue
       extinct(idxnite(icol),:) = fillvalue
       extinctnir(idxnite(icol),:) = fillvalue
       extinctuv(idxnite(icol),:) = fillvalue
       absorb(idxnite(icol),:)  = fillvalue
       asymext(idxnite(icol),:) = fillvalue
       asymvis(idxnite(icol)) = fillvalue
    end do

    call outfld('AODUVtest'//diag(list_idx),  aoduv,  pcols, lchnk)
    call outfld('AODVIStest'//diag(list_idx), aodvis, pcols, lchnk)
    call outfld('AODNIRtest'//diag(list_idx), aodnir, pcols, lchnk)
    call outfld('AODABStest'//diag(list_idx), aodabs, pcols, lchnk)
    call outfld('AODTOTtest'//diag(list_idx), aodtot, pcols, lchnk)
    call outfld('EXTINCTUVtest'//diag(list_idx),  extinctuv,  pcols, lchnk)
    call outfld('EXTINCTNIRtest'//diag(list_idx), extinctnir, pcols, lchnk)
    call outfld('EXTINCTtest'//diag(list_idx),  extinct,  pcols, lchnk)
    call outfld('ABSORBtest'//diag(list_idx),   absorb,  pcols, lchnk)
    call outfld('EXTxASYMtest'//diag(list_idx), asymext, pcols, lchnk)
    call outfld('AODxASYMtest'//diag(list_idx), asymvis, pcols, lchnk)

    deallocate(pext)
    deallocate(pabs)
    deallocate(palb)
    deallocate(pasm)

    do iaermod = 1,num_aero_models
       deallocate(aero_state(iaermod)%obj)
       nullify(aero_state(iaermod)%obj)
    end do

    deallocate(aero_state)

  end subroutine aerosol_optics_cam_sw

  !===============================================================================
  subroutine aerosol_optics_cam_lw(list_idx, state, pbuf, tauxar)

    ! calculates aerosol lw radiative properties

    integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
    type(physics_state), intent(in), target :: state    ! state variables

    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth


    real(r8) :: dopaer(pcols)
    real(r8) :: mass(pcols,pver)

    character(len=*), parameter :: prefix = 'aerosol_optics_cam_lw: '

    integer :: ibin, nbins
    integer :: iwav, ilev
    integer :: ncol, icol, istat

    type(aero_state_t), allocatable :: aero_state(:)

    class(aerosol_optics), pointer :: aero_optics

    real(r8), allocatable :: pabs(:)

    real(r8) :: relh(pcols,pver)
    real(r8) :: sate(pcols,pver)     ! saturation vapor pressure
    real(r8) :: satq(pcols,pver)     ! saturation specific humidity

    character(len=32) :: opticstype
    integer :: iaermod

    real(r8) :: lwabs(pcols,pver)
    lwabs = 0._r8
    tauxar = 0._r8

    nullify(aero_optics)

    allocate(aero_state(num_aero_models), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: aero_state')
    end if

    iaermod = 0
    if (modal_active) then
       iaermod = iaermod+1
       aero_state(iaermod)%obj => modal_aerosol_state( state, pbuf )
    else if (carma_active) then
       iaermod = iaermod+1
       aero_state(iaermod)%obj => carma_aerosol_state( state, pbuf )
    end if

    ncol = state%ncol

    mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

    allocate(pabs(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pabs')
    end if

    do iaermod = 1,num_aero_models

       nbins=aero_props(iaermod)%obj%nbins()

       do ibin = 1, nbins

          call aero_props(iaermod)%obj%optics_params(list_idx, ibin, opticstype=opticstype)

          select case (trim(opticstype))
          case('modal') ! refractive method
             aero_optics=>refractive_aerosol_optics(aero_props(iaermod)%obj, aero_state(iaermod)%obj, list_idx, ibin, ncol, pver, nswbands, nlwbands, crefwsw, crefwlw)
          case('hygroscopic_coreshell')
             ! calculate relative humidity for table lookup into rh grid
             call qsat(state%t(:ncol,:), state%pmid(:ncol,:), sate(:ncol,:), satq(:ncol,:), ncol, pver)
             relh(:ncol,:) = state%q(1:ncol,:,1) / satq(:ncol,:)
             relh(:ncol,:) = max(1.e-20_r8,relh(:ncol,:))
             aero_optics=>hygrocoreshell_aerosol_optics(aero_props(iaermod)%obj, aero_state(iaermod)%obj, list_idx, ibin, ncol, pver, relh(:ncol,:))
          case('hygroscopic_wtp')
             aero_optics=>hygrowghtpct_aerosol_optics(aero_props(iaermod)%obj, aero_state(iaermod)%obj, list_idx, ibin, ncol, pver)
          case default
             call endrun(prefix//'optics method not recognized')
          end select

          if (associated(aero_optics)) then

             do iwav = 1, nlwbands

                do ilev = 1, pver
                   call aero_optics%lw_props(ncol, ilev, iwav, pabs )

                   do icol = 1, ncol
                      dopaer(icol) = pabs(icol)*mass(icol,ilev)
                      tauxar(icol,ilev,iwav) = tauxar(icol,ilev,iwav) + dopaer(icol)
                      lwabs(icol,ilev) = lwabs(icol,ilev) + pabs(icol)
                   end do

                end do

             end do

          else
             call endrun(prefix//'aero_optics object pointer not associated')
          end if

          deallocate(aero_optics)
          nullify(aero_optics)

       end do
    end do

    call outfld('TOTABSLWtest'//diag(list_idx),  lwabs(:,:), pcols, state%lchnk)

    if (lw10um_indx>0) then
       call outfld('AODABSLWtest'//diag(list_idx), tauxar(:,:,lw10um_indx), pcols, state%lchnk)
    end if

    deallocate(pabs)

    do iaermod = 1,num_aero_models
       deallocate(aero_state(iaermod)%obj)
       nullify(aero_state(iaermod)%obj)
    end do

    deallocate(aero_state)

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
