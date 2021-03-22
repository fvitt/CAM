module mo_drydep

  !---------------------------------------------------------------------
  !       ... Dry deposition
  !---------------------------------------------------------------------

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use chem_mods,        only : gas_pcnst
  use pmgrid,           only : plev
  use spmd_utils,       only : masterproc
  use ppgrid,           only : pcols, begchunk, endchunk
  use mo_tracname,      only : solsym
  use cam_abortutils,   only : endrun
  use cam_logfile,      only : iulog
  use seq_drydep_mod,   only : nddvels => n_drydep, drydep_list, mapping
  use seq_drydep_mod,   only : n_land_type => NLUse
  use physconst,        only : karman

  implicit none

  private

  public :: drydep_init, drydep, has_drydep
  public :: drydep_update

  integer :: pan_ndx, mpan_ndx, o3_ndx, ch4_ndx, co_ndx, h2_ndx, ch3cooh_ndx

  integer :: so2_ndx, ch3cn_ndx, hcn_ndx, hcooh_ndx

  integer :: o3a_ndx,xpan_ndx,xmpan_ndx

  integer :: cohc_ndx=-1, come_ndx=-1
  integer, parameter :: NTAGS = 50
  integer :: cotag_ndx(NTAGS)
  integer :: tag_cnt

  real(r8), parameter    :: small_value = 1.e-36_r8
  real(r8), parameter    :: large_value = 1.e36_r8
  real(r8), parameter    :: diffm       = 1.789e-5_r8
  real(r8), parameter    :: diffk       = 1.461e-5_r8
  real(r8), parameter    :: difft       = 2.060e-5_r8
  real(r8), parameter    :: vonkar      = karman
  real(r8), parameter    :: ric         = 0.2_r8
  real(r8), parameter    :: r           = 287.04_r8
  real(r8), parameter    :: cp          = 1004._r8
  real(r8), parameter    :: grav        = 9.81_r8
  real(r8), parameter    :: p00         = 100000._r8
  real(r8), parameter    :: wh2o        = 18.0153_r8
  real(r8), parameter    :: ph          = 1.e-5_r8
  real(r8), parameter    :: ph_inv      = 1._r8/ph
  real(r8), parameter    :: rovcp       = r/cp
  real(r8), parameter    :: crb         = (difft/diffm)**(2._r8/3._r8)

  logical :: has_dvel(gas_pcnst) = .false.
  integer :: map_dvel(gas_pcnst) = 0

  real(r8), allocatable, dimension(:,:,:) :: dep_ra ! [s/m] aerodynamic resistance
  real(r8), allocatable, dimension(:,:,:) :: dep_rb ! [s/m] resistance across sublayer

  integer, allocatable :: spc_ndx(:) ! nddvels

  type lnd_dvel_type
     real(r8), pointer :: dvel(:,:)   ! deposition velocity over land (cm/s)
  end type lnd_dvel_type

  type(lnd_dvel_type), allocatable :: lnd(:)

  integer, parameter :: ocean_type = 7
  integer, parameter :: seaice_type = 8
  integer, parameter :: beglandtype = 7
  integer, parameter :: endlandtype = 8

contains

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_update( state, cam_in )
    use physics_types,   only : physics_state
    use camsrfexch,      only : cam_in_t

    type(physics_state), intent(in) :: state           ! Physics state variables
    type(cam_in_t),  intent(in) :: cam_in

    if (nddvels<1) return

    lnd(state%lchnk)%dvel => cam_in%depvel

  end subroutine drydep_update

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep( ocnfrac, icefrac, sfc_temp, pressure_sfc,  &
                     wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                     snow, solar_flux, dvelocity, dflx, mmr, &
                     tv, ncol, lchnk )

    !-------------------------------------------------------------------------------------
    ! combines the deposition velocities provided by the land model with deposition
    ! velocities over ocean and sea ice
    !-------------------------------------------------------------------------------------

    use ppgrid,         only : pcols
    use chem_mods,      only : gas_pcnst

#if (defined OFFLINE_DYN)
    use metdata, only: get_met_fields
#endif

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------

    real(r8), target, intent(in) :: icefrac(pcols)
    real(r8), target, intent(in) :: ocnfrac(pcols)
    integer,  intent(in)      :: ncol
    integer,  intent(in)      :: lchnk                    ! chunk number
    real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
    real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
    real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
    real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
    real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
    real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
    real(r8), intent(in)      :: rain(pcols)
    real(r8), intent(in)      :: snow(pcols)              ! snow height (m)
    real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
    real(r8), intent(in)      :: tv(pcols)                ! potential temperature
    real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst) ! constituent concentration (kg/kg)
    real(r8), intent(out)     :: dvelocity(ncol,gas_pcnst) ! deposition velocity (cm/s)
    real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)     ! deposition flux (/cm^2/s)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8) :: ocnice_dvel(ncol,gas_pcnst)
    real(r8) :: ocnice_dflx(pcols,gas_pcnst)

    real(r8), dimension(ncol) :: term    ! work array
    integer  :: ispec
    real(r8) :: lndfrac(pcols)
#if (defined OFFLINE_DYN)
    real(r8), target :: met_ocnfrac(pcols)
    real(r8), target :: met_icefrac(pcols)
#endif
    integer :: i
    real(r8), pointer :: ocnfrac_ptr(:)
    real(r8), pointer :: icefrac_ptr(:)

    lndfrac(:ncol) = 1._r8 - ocnfrac(:ncol) - icefrac(:ncol)

    where( lndfrac(:ncol) < 0._r8 )
       lndfrac(:ncol) = 0._r8
    endwhere

#if (defined OFFLINE_DYN)
    call get_met_fields(lndfrac, met_ocnfrac, met_icefrac, lchnk, ncol)
    ocnfrac_ptr=>met_ocnfrac
    icefrac_ptr=>met_icefrac
#else
    ocnfrac_ptr=>ocnfrac
    icefrac_ptr=>icefrac
#endif

    !-------------------------------------------------------------------------------------
    !   ... initialize
    !-------------------------------------------------------------------------------------
    dvelocity(:,:) = 0._r8

    !-------------------------------------------------------------------------------------
    !   ... compute the dep velocities over ocean and sea ice
    !-------------------------------------------------------------------------------------
    call drydep_xactive( sfc_temp, pressure_sfc,  &
                         wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                         snow, solar_flux, ocnice_dvel, ocnice_dflx, mmr, &
                         tv, ncol, lchnk, &
                         ocnfrc=ocnfrac_ptr, icefrc=icefrac_ptr )

    term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))

    do ispec = 1,nddvels
       !-------------------------------------------------------------------------------------
       !        ... merge the land component with the non-land component
       !            ocn and ice already have fractions factored in
       !-------------------------------------------------------------------------------------
       dvelocity(:ncol,spc_ndx(ispec)) = lnd(lchnk)%dvel(:ncol,ispec)*lndfrac(:ncol) &
                                       + ocnice_dvel(:ncol,spc_ndx(ispec))
    enddo

    !-------------------------------------------------------------------------------------
    !        ... special adjustments
    !-------------------------------------------------------------------------------------
    if( mpan_ndx>0 ) then
       dvelocity(:ncol,mpan_ndx) = dvelocity(:ncol,mpan_ndx)/3._r8
    endif
    if( xmpan_ndx>0 ) then
       dvelocity(:ncol,xmpan_ndx) = dvelocity(:ncol,xmpan_ndx)/3._r8
    endif
    if( hcn_ndx>0 ) then
       dvelocity(:ncol,hcn_ndx) = ocnice_dvel(:ncol,hcn_ndx) ! should be zero over land
    endif
    if( ch3cn_ndx>0 ) then
       dvelocity(:ncol,ch3cn_ndx) = ocnice_dvel(:ncol,ch3cn_ndx) ! should be zero over land
    endif

    ! HCOOH, use CH3COOH dep.vel
    if( hcooh_ndx > 0 .and. ch3cooh_ndx > 0 ) then
       if( has_dvel(hcooh_ndx) ) then
          dvelocity(:ncol,hcooh_ndx) = dvelocity(:ncol,ch3cooh_ndx)
       end if
    end if

    !-------------------------------------------------------------------------------------
    !        ... assign CO tags to CO
    ! put this kludge in for now ...
    !  -- should be able to set all these via the table mapping in seq_drydep_mod
    !-------------------------------------------------------------------------------------
    if( cohc_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,cohc_ndx) = dvelocity(:ncol,co_ndx)
       dflx(:ncol,cohc_ndx) = dvelocity(:ncol,co_ndx) * term(:ncol) * mmr(:ncol,plev,cohc_ndx)
    endif
    if( come_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,come_ndx) = dvelocity(:ncol,co_ndx)
       dflx(:ncol,come_ndx) = dvelocity(:ncol,co_ndx) * term(:ncol) * mmr(:ncol,plev,come_ndx)
    endif

    if ( co_ndx>0 ) then
       do i=1,tag_cnt
          dvelocity(:ncol,cotag_ndx(i)) = dvelocity(:ncol,co_ndx)
          dflx(:ncol,cotag_ndx(i)) = dvelocity(:ncol,co_ndx) * term(:ncol) * mmr(:ncol,plev,cotag_ndx(i))
       enddo
    endif

    do ispec = 1,nddvels
       !-------------------------------------------------------------------------------------
       !        ... compute the deposition flux
       !-------------------------------------------------------------------------------------
       dflx(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,spc_ndx(ispec)) * term(:ncol) * mmr(:ncol,plev,spc_ndx(ispec))
    end do

  end subroutine drydep

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_init( )
    !-------------------------------------------------------------------------------------
    ! 	... intialize interactive drydep
    !-------------------------------------------------------------------------------------
    use mo_chem_utls,  only : get_spc_ndx

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    integer :: i
    integer :: m
    integer :: astat
    integer :: ndx
    integer :: ispc

    character(len=32) :: test_name
    character(len=4) :: tag_name

    allocate(spc_ndx(nddvels))
    allocate( lnd(begchunk:endchunk) )

    do ispc = 1,nddvels

       spc_ndx(ispc) = get_spc_ndx(drydep_list(ispc))
       if (spc_ndx(ispc) < 1) then
          write(*,*) 'drydep_inti: '//trim(drydep_list(ispc))//' is not included in species set'
          call endrun('drydep_init: invalid dry deposition species')
       endif

    enddo

    if( masterproc ) then
       write(iulog,*) 'drydep_inti: following species have dry deposition'
       do i=1,nddvels
          if( len_trim(drydep_list(i)) > 0 ) then
             write(iulog,*) 'drydep_inti: '//trim(drydep_list(i))//' is requested to have dry dep'
          endif
       enddo
       write(iulog,*) 'drydep_inti:'
    endif

    !-------------------------------------------------------------------------------------
    ! 	... get species indices
    !-------------------------------------------------------------------------------------
    xpan_ndx    = get_spc_ndx( 'XPAN' )
    xmpan_ndx   = get_spc_ndx( 'XMPAN' )
    o3a_ndx     = get_spc_ndx( 'O3A' )

    ch4_ndx     = get_spc_ndx( 'CH4' )
    h2_ndx      = get_spc_ndx( 'H2' )
    co_ndx      = get_spc_ndx( 'CO' )
    pan_ndx     = get_spc_ndx( 'PAN' )
    mpan_ndx    = get_spc_ndx( 'MPAN' )
    o3_ndx      = get_spc_ndx( 'OX' )
    if( o3_ndx < 0 ) then
       o3_ndx   = get_spc_ndx( 'O3' )
    end if
    so2_ndx     = get_spc_ndx( 'SO2' )
    ch3cooh_ndx = get_spc_ndx( 'CH3COOH')

    hcn_ndx     = get_spc_ndx( 'HCN')
    ch3cn_ndx   = get_spc_ndx( 'CH3CN')

    cohc_ndx    = get_spc_ndx( 'COhc' )
    come_ndx    = get_spc_ndx( 'COme' )

    tag_cnt=0
    cotag_ndx(:)=-1
    do i = 1,NTAGS
       write(tag_name,'(a2,i2.2)') 'CO',i
       ndx = get_spc_ndx(tag_name)
       if (ndx>0) then
          tag_cnt = tag_cnt+1
          cotag_ndx(tag_cnt) = ndx
       endif
    enddo

    do i=1,nddvels
       if ( mapping(i) > 0 ) then
          test_name = drydep_list(i)
          m = get_spc_ndx( test_name )
          has_dvel(m) = .true.
          map_dvel(m) = i
       endif
    enddo

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    !---------------------------------------------------------------------------
    ! 	... allocate module variables
    !---------------------------------------------------------------------------
    allocate( dep_ra(pcols,beglandtype:endlandtype,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_ra; error = ',astat
       call endrun
    end if
    allocate( dep_rb(pcols,beglandtype:endlandtype,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_rb; error = ',astat
       call endrun
    end if

  end subroutine drydep_init

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_xactive( sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, dvel, dflx, mmr, &
                             tv, ncol, lchnk, ocnfrc, icefrc )

    !-------------------------------------------------------------------------------------
    !   code based on wesely (atmospheric environment, 1989, vol 23, p. 1293-1304) for
    !   calculation of r_c, and on walcek et. al. (atmospheric enviroment, 1986,
    !   vol. 20, p. 949-964) for calculation of r_a and r_b
    !
    !   as suggested in walcek (u_i)(u*_i) = (u_a)(u*_a)
    !   is kept constant where i represents a subgrid environment and a the
    !   grid average environment. thus the calculation proceeds as follows:
    !   va the grid averaged wind is calculated on dots
    !   z0(i) the grid averaged roughness coefficient is calculated
    !   ri(i) the grid averaged richardson number is calculated
    !   --> the grid averaged (u_a)(u*_a) is calculated
    !   --> subgrid scale u*_i is calculated assuming (u_i) given as above
    !   --> final deposotion velocity is weighted average of subgrid scale velocities
    !
    ! code written by P. Hess, rewritten in fortran 90 by JFL (August 2000)
    ! modified by JFL to be used in MOZART-2 (October 2002)
    !-------------------------------------------------------------------------------------

    use seq_drydep_mod, only: z0, rgso, rgss, ri, rclo, rcls, rlu, rac
    use seq_drydep_mod, only: seq_drydep_setHCoeff, foxd, drat
    use physconst,      only: tmelt

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer, intent(in)   :: ncol
    real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
    real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
    real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
    real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
    real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
    real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
    real(r8), intent(in)      :: rain(pcols)
    real(r8), intent(in)      :: snow(pcols)              ! snow height (m)

    real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
    real(r8), intent(in)      :: tv(pcols)                ! potential temperature
    real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst)    ! constituent concentration (kg/kg)
    real(r8), intent(out)     :: dvel(ncol,gas_pcnst)        ! deposition velocity (cm/s)
    real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)        ! deposition flux (/cm^2/s)

    integer, intent(in)     ::   lchnk                   ! chunk number

    real(r8), intent(in), optional      :: ocnfrc(pcols)
    real(r8), intent(in), optional      :: icefrc(pcols)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8), parameter :: scaling_to_cm_per_s = 100._r8
    real(r8), parameter :: rain_threshold      = 1.e-7_r8  ! of the order of 1cm/day expressed in m/s

    integer :: i, ispec, lt, m
    integer :: sndx

    real(r8) :: slope = 0._r8
    real(r8) :: z0water ! revised z0 over water
    real(r8) :: p       ! pressure at midpoint first layer
    real(r8) :: pg      ! surface pressure
    real(r8) :: es      ! saturation vapor pressure
    real(r8) :: ws      ! saturation mixing ratio
    real(r8) :: hvar    ! constant to compute xmol
    real(r8) :: h       ! constant to compute xmol
    real(r8) :: psih    ! stability correction factor
    real(r8) :: rs      ! constant for calculating rsmx
    real(r8) :: rmx     ! resistance by vegetation
    real(r8) :: zovl    ! ratio of z to  m-o length
    real(r8) :: cvarb   ! cvar averaged over landtypes
    real(r8) :: bb      ! b averaged over landtypes
    real(r8) :: ustarb  ! ustar averaged over landtypes
    real(r8) :: tc(ncol)  ! temperature in celsius
    real(r8) :: cts(ncol) ! correction to rlu rcl and rgs for frost

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,nddvels) :: heff

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location only
    !-------------------------------------------------------------------------------------
    integer                :: index_season(ncol,beglandtype:endlandtype)
    real(r8), dimension(ncol) :: tha     ! atmospheric virtual potential temperature
    real(r8), dimension(ncol) :: thg     ! ground virtual potential temperature
    real(r8), dimension(ncol) :: z       ! height of lowest level
    real(r8), dimension(ncol) :: va      ! magnitude of v on cross points
    real(r8), dimension(ncol) :: ribn    ! richardson number
    real(r8), dimension(ncol) :: qs      ! saturation specific humidity
    real(r8), dimension(ncol) :: crs     ! multiplier to calculate crs
    real(r8), dimension(ncol) :: rdc     ! part of lower canopy resistance
    real(r8), dimension(ncol) :: uustar  ! u*ustar (assumed constant over grid)
    real(r8), dimension(ncol) :: z0b     ! average roughness length over grid
    real(r8), dimension(ncol) :: wrk     ! work array
    real(r8), dimension(ncol) :: term    ! work array
    real(r8), dimension(ncol) :: resc    ! work array
    real(r8), dimension(ncol) :: lnd_frc ! work array
    logical,  dimension(ncol) :: unstable
    logical,  dimension(ncol) :: has_rain
    logical,  dimension(ncol) :: has_dew

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and landtype
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,beglandtype:endlandtype) :: rds   ! resistance for deposition of sulfate
    real(r8), dimension(ncol,beglandtype:endlandtype) :: b     ! buoyancy parameter for unstable conditions
    real(r8), dimension(ncol,beglandtype:endlandtype) :: cvar  ! height parameter
    real(r8), dimension(ncol,beglandtype:endlandtype) :: ustar ! friction velocity
    real(r8), dimension(ncol,beglandtype:endlandtype) :: xmol  ! monin-obukhov length

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location, landtype and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,beglandtype:endlandtype,gas_pcnst) :: rsmx  ! vegetative resistance (plant mesophyll)
    real(r8), dimension(ncol,beglandtype:endlandtype,gas_pcnst) :: rclx  ! lower canopy resistance
    real(r8), dimension(ncol,beglandtype:endlandtype,gas_pcnst) :: rlux  ! vegetative resistance (upper canopy)
    real(r8), dimension(ncol,beglandtype:endlandtype)           :: rlux_o3 ! vegetative resistance (upper canopy)
    real(r8), dimension(ncol,beglandtype:endlandtype,gas_pcnst) :: rgsx  ! ground resistance
    real(r8) :: vds
    logical  :: fr_lnduse(ncol,n_land_type)           ! wrking array
    real(r8) :: dewm                                  ! multiplier for rs when dew occurs

    real(r8) :: lcl_frc_landuse(ncol,n_land_type)

    !-------------------------------------------------------------------------------------
    ! jfl : mods for PAN
    !-------------------------------------------------------------------------------------
    real(r8) :: dv_pan
    real(r8) :: c0_pan(11) = (/ 0.000_r8, 0.006_r8, 0.002_r8, 0.009_r8, 0.015_r8, &
                                0.006_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.002_r8, 0.002_r8 /)
    real(r8) :: k_pan (11) = (/ 0.000_r8, 0.010_r8, 0.005_r8, 0.004_r8, 0.003_r8, &
                                0.005_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.075_r8, 0.002_r8 /)

    !-------------------------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------------------------
    do m = 1,gas_pcnst
       dvel(:,m) = 0._r8
    end do

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    !-------------------------------------------------------------------------------------
    ! define species-dependent parameters (temperature dependent)
    !-------------------------------------------------------------------------------------
    call seq_drydep_setHCoeff( ncol, sfc_temp, heff )

    do lt = beglandtype,endlandtype
       dep_ra (:,lt,lchnk)   = 0._r8
       dep_rb (:,lt,lchnk)   = 0._r8
       rds(:,lt)   = 0._r8
    end do

    !-------------------------------------------------------------------------------------
    ! season index only for ocn and sea ice
    !-------------------------------------------------------------------------------------
    index_season(:,:) = 4

    !-------------------------------------------------------------------------------------
    ! scale rain and define logical arrays
    !-------------------------------------------------------------------------------------
    has_rain(:ncol) = rain(:ncol) > rain_threshold

    !-------------------------------------------------------------------------------------
    ! loop over longitude points
    !-------------------------------------------------------------------------------------
    col_loop :  do i = 1,ncol
       p   = pressure_10m(i)
       pg  = pressure_sfc(i)
       !-------------------------------------------------------------------------------------
       ! potential temperature
       !-------------------------------------------------------------------------------------
       tha(i) = air_temp(i) * (p00/p )**rovcp * (1._r8 + .61_r8*spec_hum(i))
       thg(i) = sfc_temp(i) * (p00/pg)**rovcp * (1._r8 + .61_r8*spec_hum(i))
       !-------------------------------------------------------------------------------------
       ! height of 1st level
       !-------------------------------------------------------------------------------------
       z(i) = - r/grav * air_temp(i) * (1._r8 + .61_r8*spec_hum(i)) * log(p/pg)
       !-------------------------------------------------------------------------------------
       ! wind speed
       !-------------------------------------------------------------------------------------
       va(i) = max( .01_r8,wind_speed(i) )
       !-------------------------------------------------------------------------------------
       ! Richardson number
       !-------------------------------------------------------------------------------------
       ribn(i) = z(i) * grav * (tha(i) - thg(i))/thg(i) / (va(i)*va(i))
       ribn(i) = min( ribn(i),ric )
       unstable(i) = ribn(i) < 0._r8
       !-------------------------------------------------------------------------------------
       ! saturation vapor pressure (Pascals)
       ! saturation mixing ratio
       ! saturation specific humidity
       !-------------------------------------------------------------------------------------
       es    = 611._r8*exp( 5414.77_r8*(sfc_temp(i) - tmelt)/(tmelt*sfc_temp(i)) )
       ws    = .622_r8*es/(pg - es)
       qs(i) = ws/(1._r8 + ws)
       has_dew(i) = .false.
       if( qs(i) <= spec_hum(i) ) then
          has_dew(i) = .true.
       end if
       if( sfc_temp(i) < tmelt ) then
          has_dew(i) = .false.
       end if
       !-------------------------------------------------------------------------------------
       ! constant in determining rs
       !-------------------------------------------------------------------------------------
       tc(i) = sfc_temp(i) - tmelt
       if( sfc_temp(i) > tmelt .and. sfc_temp(i) < 313.15_r8 ) then
          crs(i) = (1._r8 + (200._r8/(solar_flux(i) + .1_r8))**2) * (400._r8/(tc(i)*(40._r8 - tc(i))))
       else
          crs(i) = large_value
       end if
       !-------------------------------------------------------------------------------------
       ! rdc (lower canopy res)
       !-------------------------------------------------------------------------------------
       rdc(i) = 100._r8*(1._r8 + 1000._r8/(solar_flux(i) + 10._r8))/(1._r8 + 1000._r8*slope)
    end do col_loop

    !-------------------------------------------------------------------------------------
    ! 	... form working arrays
    !-------------------------------------------------------------------------------------
    lcl_frc_landuse(:,:) = 0._r8

    if ( present(ocnfrc) .and. present(icefrc) ) then
       do i=1,ncol
          ! land type 7 is used for ocean
          ! land type 8 is used for sea ice
          lcl_frc_landuse(i,ocean_type) = ocnfrc(i)
          lcl_frc_landuse(i,seaice_type) = icefrc(i)
       enddo
    endif
    do lt = beglandtype,endlandtype
       do i=1,ncol
          fr_lnduse(i,lt) = lcl_frc_landuse(i,lt) > 0._r8
       enddo
    end do

    !-------------------------------------------------------------------------------------
    ! find grid averaged z0: z0bar (the roughness length) z_o=exp[S(f_i*ln(z_oi))]
    ! this is calculated so as to find u_i, assuming u*u=u_i*u_i
    !-------------------------------------------------------------------------------------
    z0b(:) = 0._r8
    do lt = beglandtype,endlandtype
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             z0b(i) = z0b(i) + lcl_frc_landuse(i,lt) * log( z0(index_season(i,lt),lt) )
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! find the constant velocity uu*=(u_i)(u*_i)
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       z0b(i) = exp( z0b(i) )
       cvarb  = vonkar/log( z(i)/z0b(i) )
       !-------------------------------------------------------------------------------------
       ! unstable and stable cases
       !-------------------------------------------------------------------------------------
       if( unstable(i) ) then
          bb = 9.4_r8*(cvarb**2)*sqrt( abs(ribn(i))*z(i)/z0b(i) )
          ustarb = cvarb * va(i) * sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*bb)) )
       else
          ustarb = cvarb * va(i)/(1._r8 + 4.7_r8*ribn(i))
       end if
       uustar(i) = va(i)*ustarb
    end do

    !-------------------------------------------------------------------------------------
    ! calculate the friction velocity for each land type u_i=uustar/u*_i
    !-------------------------------------------------------------------------------------
    do lt = beglandtype,endlandtype
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
             call ustar_calc( uustar(i), cvar(i,lt), ribn(i), z(i), z0(index_season(i,lt),lt), unstable(i), b(i,lt), ustar(i,lt) )
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! revise calculation of friction velocity and z0 over water
    !-------------------------------------------------------------------------------------
    lt = ocean_type
    do i = 1,ncol
       if( fr_lnduse(i,lt) ) then
          z0water    = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
          cvar(i,lt) = vonkar/(log( z(i)/z0water ))
          call ustar_calc( uustar(i), cvar(i,lt), ribn(i), z(i), z0water, unstable(i), b(i,lt), ustar(i,lt) )
       end if
    end do

    !-------------------------------------------------------------------------------------
    ! compute monin-obukhov length for unstable and stable conditions/ sublayer resistance
    !-------------------------------------------------------------------------------------
    do lt = beglandtype,endlandtype
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             hvar = (va(i)/0.74_r8) * (tha(i) - thg(i)) * (cvar(i,lt)**2)
             if( unstable(i) ) then
                h = hvar*(1._r8 - (9.4_r8*ribn(i)/(1._r8 + 5.3_r8*b(i,lt))))
             else
                h = hvar/((1._r8+4.7_r8*ribn(i))**2)
             end if
             xmol(i,lt) = thg(i) * ustar(i,lt) * ustar(i,lt) / (vonkar * grav * h)
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! psih
    !-------------------------------------------------------------------------------------
    do lt = beglandtype,endlandtype
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( xmol(i,lt) < 0._r8 ) then
                zovl = z(i)/xmol(i,lt)
                zovl = max( -1._r8,zovl )
                psih = exp( .598_r8 + .39_r8*log( -zovl ) - .09_r8*(log( -zovl ))**2 )
                vds  = 2.e-3_r8*ustar(i,lt) * (1._r8 + (300/(-xmol(i,lt)))**0.666_r8)
             else
                zovl = z(i)/xmol(i,lt)
                zovl = min( 1._r8,zovl )
                psih = -5._r8 * zovl
                vds  = 2.e-3_r8*ustar(i,lt)
             end if
             dep_ra (i,lt,lchnk) = (vonkar - psih*cvar(i,lt))/(ustar(i,lt)*vonkar*cvar(i,lt))
             dep_rb (i,lt,lchnk) = (2._r8/(vonkar*ustar(i,lt))) * crb
             rds(i,lt) = 1._r8/vds
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! surface resistance : depends on both land type and species
    ! land types are computed seperately, then resistance is computed as average of values
    ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1
    !
    ! compute rsmx = 1/(rs+rm) : multiply by 3 if surface is wet
    !-------------------------------------------------------------------------------------
    species_loop1 :  do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          m = map_dvel(ispec)
          do lt = beglandtype,endlandtype
             do i = 1,ncol
                if( fr_lnduse(i,lt) ) then
                   sndx = index_season(i,lt)
                   if( ispec == o3_ndx .or. ispec == o3a_ndx .or. ispec == so2_ndx ) then
                      rmx = 0._r8
                   else
                      rmx = 1._r8/(heff(i,m)/3000._r8 + 100._r8*foxd(m))
                   end if
                   cts(i) = 1000._r8*exp( - tc(i) - 4._r8 ) ! correction for frost
                   rgsx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rgss(sndx,lt))) + (foxd(m)/rgso(sndx,lt)))
                   !-------------------------------------------------------------------------------------
                   ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
                   !-------------------------------------------------------------------------------------
                   if( ispec == h2_ndx .or. ispec == co_ndx .or. ispec == ch4_ndx ) then
                      !-------------------------------------------------------------------------------------
                      ! no deposition on snow, ice, desert, and water
                      !-------------------------------------------------------------------------------------
                      rgsx(i,lt,ispec) = large_value
                   end if
                   if( lt == ocean_type ) then
                      rclx(i,lt,ispec) = large_value
                      rsmx(i,lt,ispec) = large_value
                      rlux(i,lt,ispec) = large_value
                   else
                      rs = ri(sndx,lt)*crs(i)
                      if ( has_dew(i) .or. has_rain(i) ) then
                         dewm = 3._r8
                      else
                         dewm = 1._r8
                      end if
                      rsmx(i,lt,ispec) = (dewm*rs*drat(m) + rmx)
                      !-------------------------------------------------------------------------------------
                      ! jfl : special case for PAN
                      !-------------------------------------------------------------------------------------
                      if( ispec == pan_ndx .or. ispec == xpan_ndx ) then
                         dv_pan =  c0_pan(lt) * (1._r8 - exp( -k_pan(lt)*(dewm*rs*drat(m))*1.e-2_r8 ))
                         if( dv_pan > 0._r8 .and. sndx /= 4 ) then
                            rsmx(i,lt,ispec) = ( 1._r8/dv_pan )
                         end if
                      end if
                      rclx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rcls(sndx,lt))) + (foxd(m)/rclo(sndx,lt)))
                      rlux(i,lt,ispec) = cts(i) + rlu(sndx,lt)/(1.e-5_r8*heff(i,m) + foxd(m))
                   end if
                end if
             end do
          end do
       end if
    end do species_loop1

    lt = seaice_type
    do i = 1,ncol
       if( fr_lnduse(i,lt) ) then
          sndx = index_season(i,lt)
          !-------------------------------------------------------------------------------------
          ! 	... no effect if sfc_temp < O C
          !-------------------------------------------------------------------------------------
          if( sfc_temp(i) > tmelt ) then
             if( has_dew(i) ) then
                rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + rlu(sndx,lt))
                if( o3_ndx > 0 ) then
                   rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                endif
                if( o3a_ndx > 0 ) then
                   rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                endif
             end if
             if( has_rain(i) ) then
                rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + 3._r8*rlu(sndx,lt))
                if( o3_ndx > 0 ) then
                   rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                endif
                if( o3a_ndx > 0 ) then
                   rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                endif
             end if
          end if

          if ( o3_ndx > 0 ) then
             rclx(i,lt,o3_ndx) = cts(i) + rclo(index_season(i,lt),lt)
             rlux(i,lt,o3_ndx) = cts(i) + rlux(i,lt,o3_ndx)
          end if
          if ( o3a_ndx > 0 ) then
             rclx(i,lt,o3a_ndx) = cts(i) + rclo(index_season(i,lt),lt)
             rlux(i,lt,o3a_ndx) = cts(i) + rlux(i,lt,o3a_ndx)
          end if

       end if
    end do

    lt = seaice_type
    species_loop2 : do ispec = 1,gas_pcnst
       m = map_dvel(ispec)
       if( has_dvel(ispec) ) then
          if( ispec /= o3_ndx .and. ispec /= o3a_ndx .and. ispec /= so2_ndx ) then
             do i = 1,ncol
                if( fr_lnduse(i,lt) ) then
                   !-------------------------------------------------------------------------------------
                   ! no effect if sfc_temp < O C
                   !-------------------------------------------------------------------------------------
                   if( sfc_temp(i) > tmelt ) then
                      if( has_dew(i) ) then
                         rlux(i,lt,ispec) = 1._r8/((1._r8/(3._r8*rlux(i,lt,ispec))) &
                              + 1.e-7_r8*heff(i,m) + foxd(m)/rlux_o3(i,lt))
                      end if
                   end if

                end if
             end do
          else if( ispec == so2_ndx ) then
             do i = 1,ncol
                if( fr_lnduse(i,lt) ) then
                   !-------------------------------------------------------------------------------------
                   ! no effect if sfc_temp < O C
                   !-------------------------------------------------------------------------------------
                   if( sfc_temp(i) > tmelt ) then
                      if( qs(i) <= spec_hum(i) ) then
                         rlux(i,lt,ispec) = 100._r8
                      end if
                      if( has_rain(i) ) then
                         rlux(i,lt,ispec) = 15._r8*rlu(index_season(i,lt),lt)/(5._r8 + 3.e-3_r8*rlu(index_season(i,lt),lt))
                      end if
                   end if
                   rclx(i,lt,ispec) = cts(i) + rcls(index_season(i,lt),lt)
                   rlux(i,lt,ispec) = cts(i) + rlux(i,lt,ispec)

                end if
             end do
          end if
       end if
    end do species_loop2

    !-------------------------------------------------------------------------------------
    ! compute rc
    !-------------------------------------------------------------------------------------
    term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))
    species_loop3 : do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          wrk(:) = 0._r8
          lt_loop: do lt = beglandtype,endlandtype
             do i = 1,ncol
                if (fr_lnduse(i,lt)) then
                   resc(i) = 1._r8/( 1._r8/rsmx(i,lt,ispec) + 1._r8/rlux(i,lt,ispec) &
                                   + 1._r8/(rdc(i) + rclx(i,lt,ispec)) &
                                   + 1._r8/(rac(index_season(i,lt),lt) + rgsx(i,lt,ispec)))

                   resc(i) = max( 10._r8,resc(i) )

                   lnd_frc(i) = lcl_frc_landuse(i,lt)
                endif
             enddo
             !-------------------------------------------------------------------------------------
             ! 	... compute average deposition velocity
             !-------------------------------------------------------------------------------------
             select case( solsym(ispec) )
             case( 'SO2' )
                if( lt == ocean_type ) then
                   where( fr_lnduse(:ncol,lt) )
                      ! assume no surface resistance for SO2 over water`
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk))
                   endwhere
                else
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:))
                   endwhere
                end if

                !  JFL - increase in dry deposition of SO2 to improve bias over US/Europe
                wrk(:) = wrk(:) * 2._r8

             case( 'SO4' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + rds(:,lt))
                endwhere
             case( 'NH4', 'NH4NO3', 'XNH4NO3' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + 0.5_r8*rds(:,lt))
                endwhere

             !-------------------------------------------------------------------------------------
             !  ... special case for Pb (for consistency with offline code)
             !-------------------------------------------------------------------------------------
             case( 'Pb' )
                if( lt == ocean_type ) then
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:) = wrk(:) + lnd_frc(:) * 0.05e-2_r8
                   endwhere
                else
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_r8
                   endwhere
                end if

             !-------------------------------------------------------------------------------------
             !  ... special case for carbon aerosols
             !-------------------------------------------------------------------------------------
             case( 'CB1', 'CB2', 'OC1', 'OC2', 'SOAM', 'SOAI', 'SOAT', 'SOAB','SOAX' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.10e-2_r8
                endwhere

             !-------------------------------------------------------------------------------------
             ! deposition over ocean for HCN, CH3CN
             !    velocity estimated from aircraft measurements (E.Apel, INTEX-B)
             !-------------------------------------------------------------------------------------
             case( 'HCN','CH3CN' )
                if( lt == ocean_type ) then ! over ocean only
                   where( fr_lnduse(:ncol,lt) .and. snow(:ncol) < 0.01_r8  )
                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_r8
                   endwhere
                end if
             case default
                where( fr_lnduse(:ncol,lt) )
                   wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:ncol))
                endwhere
             end select
          end do lt_loop
          dvel(:ncol,ispec) = wrk(:ncol) * scaling_to_cm_per_s
          dflx(:ncol,ispec) = term(:ncol) * dvel(:ncol,ispec) * mmr(:ncol,plev,ispec)
       end if

    end do species_loop3

  contains

    subroutine ustar_calc( uustar, cvar, ribn, z, z0, unstable, b, ustr )
      real(r8), intent(in) :: uustar, cvar, ribn, z, z0
      logical, intent(in) :: unstable

      real(r8), intent(out) :: b, ustr

      if( unstable ) then
         b = 9.4_r8*(cvar      **2)* sqrt( abs(ribn   )*z   /z0 )
         ustr = sqrt( cvar      *uustar*   sqrt( 1._r8 - (9.4_r8*ribn   /(1._r8 + 7.4_r8*b      )) ) )
      else
         ustr = sqrt( cvar      *uustar   /(1._r8 + 4.7_r8*ribn   ) )
      end if

    end subroutine ustar_calc

  end subroutine drydep_xactive

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  function has_drydep( name )

    character(len=*), intent(in) :: name

    logical :: has_drydep
    integer :: i

    has_drydep = .false.

    do i=1,nddvels
       if ( trim(name) == trim(drydep_list(i)) ) then
         has_drydep = .true.
         exit
       endif
    enddo

  endfunction has_drydep

end module mo_drydep
