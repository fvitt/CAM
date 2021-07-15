module mee_ap_ionization
  use shr_kind_mod, only: r8 => shr_kind_r8
  use solar_parms_data, only: Ap=>solar_parms_ap ! geomag activity index
  use mo_apex, only: alatm !(icol,lchnk) apex mag latitude at each geographic grid point (radians)
  use ppgrid, only: pcols, pver
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc
  use cam_abortutils, only : endrun

  implicit none

  private
  public :: mee_ap_ion_readnl
  public :: mee_ap_ion_init
  public :: mee_ap_ionpairs

  integer, parameter :: nbins=100

  real(r8) :: energies(nbins)
  real(r8) :: denergies(nbins) ! width of each energy bin

  logical :: mee_ap_ion_inline = .false.

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ap_ion_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use spmd_utils,     only: mpicom, mpi_logical, masterprocid

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'mee_ap_ion_readnl'

    namelist /mee_ap_ion_nl/ mee_ap_ion_inline


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
    if ( masterproc ) then
       write(iulog,*) subname//':: mee_ap_ion_inline = ',mee_ap_ion_inline
    endif

  end subroutine mee_ap_ion_readnl

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ap_ion_init()
    use cam_history, only: addfld
    use mee_ap_util_mod, only: gen_energy_grid

    if (.not.mee_ap_ion_inline) return

    call gen_energy_grid(nbins, energies, denergies)

    call addfld( 'APMEEionprs', (/ 'lev' /), 'A', 'pairs/cm3/sec', 'Ap generated MEE ionization rate' )
  end subroutine mee_ap_ion_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine mee_ap_ionpairs(ncol,lchnk, pmid, temp, ionpairs )

    use physconst, only: pi
    use physconst, only: mbarv  ! kg/kmole
    use physconst, only: gravit
    use physconst, only: rairv  ! composition dependent gas constant (J/K/kg)
    use physconst, only: boltz  ! Boltzman's constant (J/K/molecule)
    use physconst, only: avogad ! Avogadro's number (molecules/kmole)

    use cam_history, only : outfld
    use mee_ap_util_mod, only: maglat2lshell, FluxSpectrum, iprmono

    integer,  intent(in) :: ncol,lchnk
    real(r8), intent(in) :: pmid(:,:)
    real(r8), intent(in) :: temp(:,:)
    real(r8), intent(out) :: ionpairs(:,:)

    real(r8) :: rho(pcols,pver)
    real(r8) :: scaleh(pcols,pver)

    real(r8) :: lshell
    real(r8) :: flux_sd(nbins)
    real(r8) :: flux(nbins)
    real(r8) :: ipr(nbins,pver)

    integer :: i,k

    if (.not.mee_ap_ion_inline) return

    rho(:ncol,:) = pmid(:ncol,:)/(rairv(:ncol,:,lchnk)*temp(:ncol,:)) ! kg/m3
    rho(:ncol,:) = rho(:ncol,:)*1.0e-3_r8 !  kg/m3 --> g/cm3

    scaleh(:ncol,:) =  avogad * boltz*temp(:ncol,:)/(mbarv(:ncol,:,lchnk)*gravit ) ! m
    scaleh(:ncol,:) = scaleh(:ncol,:) * 1.0e2_r8 ! m -> cm

    ionpairs(:,:) = 0._r8

    do i = 1,ncol

       if ( abs(alatm(i,lchnk)) < 85._r8*pi/180._r8 ) then

          ! get L-shell value corresponeding to the column geo-mag latitude
          lshell = maglat2lshell( alatm(i,lchnk) )

          ! calculate the top of the atmosphere energetic electron energy spectrum
          flux_sd(:) = FluxSpectrum(energies, lshell, Ap)

          ! van de Kamp is per steradian (electrons / (cm2 sr s keV))
          ! assume flux is isotropic inside a nominal bounce loss cone (BLC) angle
          ! of 66.3˚. The area of the BLC in sr is 2pi(1-cosd(66.3))
          flux(:) = 2._r8*pi*(1._r8-cos(66.3*pi/180._r8)) * flux_sd(:)

          ! calculate the IPR as a function f height and energy
          ipr(:,:) = iprmono(energies, flux, rho(i,:), scaleh(i,:))

          ! integrate across the energy range to get total IPR
          do k=1,pver
             ionpairs(i,k) = ionpairs(i,k) + sum(ipr(:,k)*denergies(:))
          end do

       end if

    end do

    call outfld( 'APMEEionprs', ionpairs(:ncol,:), ncol, lchnk )

  end subroutine mee_ap_ionpairs


end module mee_ap_ionization
