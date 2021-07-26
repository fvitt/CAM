module mee_ap_ionization
  use shr_kind_mod, only: r8 => shr_kind_r8
  use solar_parms_data, only: Ap=>solar_parms_ap ! geomag activity index
  use mo_apex, only: alatm !(icol,lchnk) apex mag latitude at each geographic grid point (radians)
  use ppgrid, only: pcols, pver
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc
  use cam_abortutils, only: endrun
  use mee_ap_util_mod,only: mee_ap_init, mee_ap_error, mee_ap_iprs

  implicit none

  private
  public :: mee_ap_ion_readnl
  public :: mee_ap_ion_init
  public :: mee_ap_ionpairs

  logical :: mee_ap_ion_inline = .false.
  real(r8) :: mee_ap_ion_blc = -huge(1._r8) ! bounce cone angle (degrees)

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ap_ion_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use spmd_utils, only: mpicom, mpi_logical, mpi_real8, masterprocid

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'mee_ap_ion_readnl'

    namelist /mee_ap_ion_nl/ mee_ap_ion_inline, mee_ap_ion_blc


    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'mee_ap_ion_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, mee_ap_ion_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(mee_ap_ion_inline, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(mee_ap_ion_blc, 1, mpi_real8, masterprocid, mpicom, ierr)
    if ( masterproc ) then
       write(iulog,*) subname//':: mee_ap_ion_inline = ', mee_ap_ion_inline
       write(iulog,*) subname//':: mee_ap_ion_blc = ', mee_ap_ion_blc
    endif

  end subroutine mee_ap_ion_readnl

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ap_ion_init()
    use cam_history, only: addfld

    integer :: err

    if (.not.mee_ap_ion_inline) return

    call mee_ap_init(mee_ap_ion_blc,err)
    if (err==mee_ap_error) then
       call endrun('mee_ap_ion_init: not able to initialize Ap based MEE ionization')
    endif

    call addfld( 'APMEEionprs', (/ 'lev' /), 'A', 'pairs/cm3/sec', 'Ap generated MEE ionization rate' )
  end subroutine mee_ap_ion_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ap_ionpairs(ncol,lchnk, pmid, temp, ionpairs )

    use physconst, only: mbarv  ! kg/kmole
    use physconst, only: gravit
    use physconst, only: rairv  ! composition dependent gas constant (J/K/kg)
    use physconst, only: boltz  ! Boltzman's constant (J/K/molecule)
    use physconst, only: avogad ! Avogadro's number (molecules/kmole)

    use cam_history, only : outfld

    integer,  intent(in) :: ncol,lchnk
    real(r8), intent(in) :: pmid(:,:)
    real(r8), intent(in) :: temp(:,:)
    real(r8), intent(out) :: ionpairs(:,:)

    real(r8) :: rho(pcols,pver)
    real(r8) :: scaleh(pcols,pver)
    integer :: err

    if (.not.mee_ap_ion_inline) then
       ionpairs(:,:) = 0._r8
       return
    endif

    rho(:ncol,:) = pmid(:ncol,:)/(rairv(:ncol,:,lchnk)*temp(:ncol,:)) ! kg/m3
    rho(:ncol,:) = rho(:ncol,:)*1.0e-3_r8 ! kg/m3 --> g/cm3

    scaleh(:ncol,:) = avogad * boltz*temp(:ncol,:)/(mbarv(:ncol,:,lchnk)*gravit) ! m
    scaleh(:ncol,:) = scaleh(:ncol,:) * 1.0e2_r8 ! m -> cm

    ionpairs(:ncol,:) = mee_ap_iprs(ncol, pver, rho(:ncol,:), scaleh(:ncol,:), Ap, status=err, maglat=alatm(:ncol,lchnk))
    if (err==mee_ap_error) then
       call endrun('mee_ap_ionpairs: error in Ap based MEE ionization calculation')
    endif

    call outfld( 'APMEEionprs', ionpairs(:ncol,:), ncol, lchnk )

  end subroutine mee_ap_ionpairs


end module mee_ap_ionization
