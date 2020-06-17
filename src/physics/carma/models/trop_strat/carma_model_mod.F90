!! This CARMA model is for dust aerosols and is based upon Su & Toon, JGR, 2009;
!! Su & Toon, ACP 2011.
!!
!! These dust are not currently radiatively active and do not replace the dust
!! in CAM; however, this is something that could be done in the future.
!!
!! This module defines several constants needed by CARMA, extends a couple of CARMA
!! interface methods:
!!
!!   - CARMA_DefineModel()
!!   - CARMA_EmitParticle()
!!
!! and adds some local functions used to do sea salt emission:
!!
!!   - CARMA_SurfaceWind()
!!   - CARMA_SurfaceWind()
!!   - WeibullWind()
!!
!! @version April-2020
!! @author  Simone Tilmes, Lin Su, Pengfei Yu, Chuck Bardeen 
module carma_model_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagas_mod
  use carmagroup_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use carma_flags_mod
  use carma_model_flags_mod

  use spmd_utils,     only: masterproc
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use radconstants,   only: nswbands, nlwbands
  use cam_abortutils, only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc

  implicit none

  private

  ! Declare the public methods.
  public CARMA_DefineModel
  public CARMA_Detrain
  public CARMA_DiagnoseBins
  public CARMA_DiagnoseBulk
  public CARMA_EmitParticle
  public CARMA_InitializeModel
  public CARMA_InitializeParticle
  public CARMA_WetDeposition
  public CARMA_SaltFlux

  ! Declare public constants
  integer, public, parameter      :: NGROUP   = 1               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 5               !! Number of particle elements
  integer, public, parameter      :: NBIN     = 20              !! Number of particle bins
  integer, public, parameter      :: NSOLUTE  = 0               !! Number of particle solutes
  integer, public, parameter      :: NGAS     = 2               !! Number of gases

  ! These need to be defined, but are only used when the particles are radiatively active.
  integer, public, parameter      :: NMIE_RH  = 10              !! Number of relative humidities for mie calculations
  real(kind=f), public            :: mie_rh(NMIE_RH)

  ! Defines whether the groups should undergo deep convection in phase 1 or phase 2.
  ! Water vapor and cloud particles are convected in phase 1, while all other constituents
  ! are done in phase 2.
  logical, public                 :: is_convtran1(NGROUP) = .false.  !! Should the group be transported in the first phase?

  ! Define any particle compositions that are used. Each composition type
  ! should have a unique number.
  integer, public, parameter      :: I_H2SO4          = 1       !! H2SO4 coposition
  integer, public, parameter      :: I_OC             = 2       !! OC composition
  integer, public, parameter      :: I_BC             = 3       !! BC composition
  integer, public, parameter      :: I_DUST           = 4       !! dust composition
  integer, public, parameter      :: I_SALT           = 2       !! sea salt composition

  integer, public, parameter      :: I_GRP_MIXAER     = 1       !! mixed aerosol

  integer, public, parameter      :: I_ELEM_CRMIX     = 1       !! aerosol
  integer, public, parameter      :: I_ELEM_CROC      = 2       !! organics aerosol
  integer, public, parameter      :: I_ELEM_CRBC      = 3       !! black carbon
  integer, public, parameter      :: I_ELEM_CRDUST    = 4       !! dust aerosol
  integer, public, parameter      :: I_ELEM_CRSALT    = 5       !! sea salt aerosol

  integer, public, parameter      :: I_GAS_H2O        = 1              !! water vapor
  integer, public, parameter      :: I_GAS_H2SO4      = 2              !! sulphuric acid  

  real(kind=f), public, parameter         :: Kappa_OC = 0.5_f      !! hygroscopicity of OC
  real(kind=f), public, parameter         :: Kappa_BC = 0.1_f
  real(kind=f), public, parameter         :: Kappa_DUST = 0.2_f
  real(kind=f), public, parameter         :: Kappa_SALT = 1.0_f
  real(kind=f), public, parameter         :: Kappa_SULF = 0.5_f

  real(kind=f), public, parameter         :: RHO_obc  = 1.35_f                  !! dry density of smoke aerosol
  real(kind=f), public, parameter         :: RHO_DUST = 2.65_f                  !! dry density of dust particles (g/cm^3) -Lin Su
  real(kind=f), public, parameter         :: RHO_SALT = 2.65_f                  !! dry density of sea salt particles (g/cm)
  real(kind=f), public, parameter         :: RHO_SULFATE  = 1.923_f     !! dry density of sulfate particles (g/cm3)

  ! for sea salt flux calculation
  real(r8), parameter             :: uth_salt = 4._r8                !! threshold wind velocity  


  ! for dust calculation
  real(kind=f), parameter         :: rClay = 1e-4_f         !! silt/clay particle radius boundary (cm)

  integer                         :: nClay                  !! Number of clay bins (r < 1 um)
  integer                         :: nSilt                  !! Number of silt bins
  real(kind=f)                    :: clay_mf(NBIN)          !! clay mass fraction (fraction)  
  real(kind=f)                    :: kappa(NBIN)            !! Kappa for Aerosols   
  real(kind=f)                    :: hygro(NBIN)            !! Hygroscopicity for Aerosols  
  real(kind=f), allocatable, dimension(:,:) :: soil_factor  !! Soil Erosion Factor (fraction) 
  real(kind=f), public, parameter :: WTMOL_H2SO4    = 98.078479_f    !! molecular weight of sulphuric acid
  real(kind=f), allocatable, dimension(:,:,:)   :: PCT_LeafArea                 !! Soil Erosion Factor (fraction) 



! NOTE: The WeibullK distribution is not currently supported, since the coefficients are not
! generated. This can be added later.
  real(r8), allocatable, dimension(:,:) :: Weibull_k            ! Weibull K(nlat,nlon
  real(kind=f), public, parameter     :: rmin_MIXAER     = 5e-6_f     ! minimum radius (cm)
  real(kind=f), public, parameter     :: vmrat_MIXAER    = 2.2588_f    !2.4610_f        ! volume ratio
! sea-salt
!  real(kind=f), parameter            :: rmin_MIXAER     = 1e-6_f    ! minimum radius (cm)
! real(kind=f), parameter            :: vmrat_MIXAER    = 4.32_f    ! volume ratio


contains


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DefineModel(carma, rc)
    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure
    
    ! Local variables
    integer                            :: LUNOPRT              ! logical unit number for output
    logical                            :: do_print             ! do print output?
!   real(kind=f), parameter            :: RHO_DUST = 2.65_f    ! dry density of dust particles (g/cm^3) -Lin Su 
!   real(kind=f), parameter            :: RHO_SALT = 2.65_f    ! dry density of sea salt particles (g/cm)
!dust    real(kind=f), parameter            :: rmin     = 1.19e-5_f ! minimum radius (cm)
!dust    real(kind=f), parameter            :: vmrat    = 2.371_f   ! volume ratio
    complex(kind=f)                    :: refidx(NWAVE)        ! refractice indices
    
    ! Default return code.
    rc = RC_OK    
    
    ! Report model specific namelist configuration parameters.
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_DefineModel: CARMA_Get failed.")
    
      if (do_print) write(LUNOPRT,*) ''
      if (do_print) write(LUNOPRT,*) 'CARMA ', trim(carma_model), ' specific settings :'
      if (do_print) write(LUNOPRT,*) '  carma_soilerosion_file = ', carma_soilerosion_file
      if (do_print) write(LUNOPRT,*) '  carma_seasalt_emis = ', trim(carma_seasalt_emis)
    end if

    ! Define the Groups
    !
    ! NOTE: For CAM, the optional do_wetdep and do_drydep flags should be
    ! defined. If wetdep is defined, then the optional solubility factor
    ! should also be defined.
    !call CARMAGROUP_Create(carma, I_GRP_MIXAER, "mixed aerosol", rmin_MIXAER, vmrat_MIXAER, I_SPHERE, 1._f, .false., &
    !                       rc, do_wetdep=.true., do_drydep=.true., solfac=0.2_f, &
    !                       scavcoef=0.1_f, shortname="CRMIX", refidx=refidx, &
    !                       refidxS=refidxS, refidxC=refidxC, do_mie=.true., &
    !                       irhswell=I_MIX, irhswcomp=I_SWG_URBAN, icoreshell=1,imiertn=I_MIERTN_TOON1981)
    !if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    call CARMAGROUP_Create(carma, I_GRP_MIXAER, "mixed aerosol", rmin_MIXAER, vmrat_MIXAER, I_SPHERE, 1._f, .false., &
                           rc, do_wetdep=.true., do_drydep=.true., solfac=0.2_f, &
                           scavcoef=0.1_f, shortname="CRMIX", refidx=refidx, &
                           irhswell=I_MIX, imiertn=I_MIERTN_TOON1981)
   !                       irhswell=I_MIX, irhswcomp=I_SWG_URBAN,imiertn=I_MIERTN_TOON1981)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    
    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    !call CARMAELEMENT_Create(carma, I_ELEM_PURSUL, I_GRP_PURSUL, "Sulfate", RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="SULF")
    !if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRMIX,  I_GRP_MIXAER, "Sulfate in mixed sulfate", RHO_SULFATE, I_INVOLATILE, I_H2SO4, rc, kappa=Kappa_SULF, shortname="CRMIX")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CROC,   I_GRP_MIXAER, "organic carbon", RHO_obc, I_COREMASS, I_OC, rc, kappa=Kappa_OC, shortname="CROC")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRBC,   I_GRP_MIXAER, "black carbon", RHO_obc, I_COREMASS, I_BC, rc, kappa=Kappa_BC,  shortname="CRBC")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRDUST, I_GRP_MIXAER, "dust", RHO_DUST, I_COREMASS, I_DUST, rc,  kappa=Kappa_DUST, shortname="CRDUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_CRSALT, I_GRP_MIXAER, "SALT in mixed sulfate", RHO_SALT, I_COREMASS, I_SALT, rc,  kappa=Kappa_SALT, shortname="CRSALT")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    
    ! Define the Solutes
    


    ! Define the Gases
    call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, &
                         rc, shortname = "Q", ds_threshold=-0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    call CARMAGAS_Create(carma, I_GAS_H2SO4, "Sulfuric Acid", WTMOL_H2SO4, I_VAPRTN_H2SO4_AYERS1980, &
                         I_GCOMP_H2SO4, rc, shortname = "H2SO4", ds_threshold=-0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    
    ! Define the Processes

!   call CARMA_AddGrowth(carma, I_ELEM_PURSUL, I_GAS_H2SO4, rc)
!   if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddGrowth(carma, I_ELEM_CRMIX, I_GAS_H2SO4, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

!   call CARMA_AddNucleation(carma, I_ELEM_PURSUL, I_ELEM_PURSUL, I_HOMNUC, 0._f, rc, igas=I_GAS_H2SO4)
!   if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

!   call CARMA_AddCoagulation(carma, I_GRP_PURSUl, I_GRP_PURSUL, I_GRP_PURSUL, I_COLLEC_FUCHS, rc)
!   if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

!   call CARMA_AddCoagulation(carma, I_GRP_PURSUl, I_GRP_MIXAER, I_GRP_MIXAER, I_COLLEC_DATA, rc)
!   if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_MIXAER, I_GRP_MIXAER, I_GRP_MIXAER, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

!   call pbuf_add_field('SADSULF', 'global', dtype_r8, (/pcols, pver/), ipbuf4sad)

    
    return
  end subroutine CARMA_DefineModel


  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009 
  !!  @author  Chuck Bardeen 
  !!
  !!  @see CARMASTATE_SetDetrain
  subroutine CARMA_Detrain(carma, cstate, cam_in, dlf, state, icol, dt, rc, rliq, prec_str, snow_str, &
     tnd_qsnow, tnd_nsnow)
    use camsrfexch,         only: cam_in_t
    use physconst,          only: latice, latvap, cpair

    implicit none

    type(carma_type), intent(in)         :: carma            !! the carma object
    type(carmastate_type), intent(inout) :: cstate           !! the carma state object
    type(cam_in_t),  intent(in)          :: cam_in           !! surface input
    real(r8), intent(in)                 :: dlf(pcols, pver) !! Detraining cld H20 from convection (kg/kg/s)
    type(physics_state), intent(in)      :: state            !! physics state variables
    integer, intent(in)                  :: icol             !! column index
    real(r8), intent(in)                 :: dt               !! time step (s)
    integer, intent(out)                 :: rc               !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(out), optional      :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(out), optional      :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)

    ! Default return code.
    rc = RC_OK
    
    return
  end subroutine CARMA_Detrain


  !! For diagnostic groups, sets up up the CARMA bins based upon the CAM state.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBins(carma, cstate, state, pbuf, icol, dt, rc, rliq, prec_str, snow_str)
    use time_manager,     only: is_first_step

    implicit none

    type(carma_type), intent(in)          :: carma        !! the carma object
    type(carmastate_type), intent(inout)  :: cstate       !! the carma state object
    type(physics_state), intent(in)       :: state        !! physics state variables
    type(physics_buffer_desc), pointer    :: pbuf(:)      !! physics buffer
    integer, intent(in)                   :: icol         !! column index
    real(r8), intent(in)                  :: dt           !! time step
    integer, intent(out)                  :: rc           !! return code, negative indicates failure
    real(r8), intent(in), optional        :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional     :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional     :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    
    real(r8)                             :: mmr(pver) !! elements mass mixing ratio
    integer                              :: ibin      !! bin index
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the mass in each bin from the CAM state.
    
    return
  end subroutine CARMA_DiagnoseBins


  !! For diagnostic groups, determines the tendencies on the CAM state from the CARMA bins.
  !!
  !!  @version July-2009 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, tnd_qsnow, tnd_nsnow, re_ice)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma     !! the carma object
    type(carmastate_type), intent(inout) :: cstate    !! the carma state object
    type(cam_out_t),      intent(inout)  :: cam_out   !! cam output to surface models
    type(physics_state), intent(in)      :: state     !! physics state variables
    type(physics_buffer_desc), pointer   :: pbuf(:)   !! physics buffer
    type(physics_ptend), intent(inout)   :: ptend     !! constituent tendencies
    integer, intent(in)                  :: icol      !! column index
    real(r8), intent(in)                 :: dt        !! time step
    integer, intent(out)                 :: rc        !! return code, negative indicates failure
    real(r8), intent(inout), optional    :: rliq(pcols)      !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), intent(inout), optional    :: prec_str(pcols)  !! [Total] sfc flux of precip from stratiform (m/s) 
    real(r8), intent(inout), optional    :: snow_str(pcols)  !! [Total] sfc flux of snow from stratiform (m/s)
    real(r8), intent(inout), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), intent(inout), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(inout), optional    :: tnd_qsnow(pcols,pver) !! snow mass tendency (kg/kg/s)
    real(r8), intent(inout), optional    :: tnd_nsnow(pcols,pver) !! snow number tendency (#/kg/s)
    real(r8), intent(out), optional      :: re_ice(pcols,pver)    !! ice effective radius (m)
    
    ! Default return code.
    rc = RC_OK
    
    ! By default, do nothing. If diagnosed groups exist, this needs to be replaced by
    ! code to determine the bulk mass from the CARMA state.
    
    return
  end subroutine CARMA_DiagnoseBulk


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Lin Su, Pengfei Yu, Chuck Bardeen
  !! @version Dec-2010
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use phys_grid,     only: get_lon_all_p, get_lat_all_p
    use camsrfexch,    only: cam_in_t
    use cam_history,   only: outfld
    
    implicit none
    
    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure
    
    integer      :: ilat(pcols)             ! latitude index 
    integer      :: ilon(pcols)             ! longitude index
    integer      :: lchnk                   ! chunk identifier
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    integer      :: igroup                  ! the index of the carma aerosol group
    character(len=32) :: shortname          ! the shortname of the group
    
    ! -------- local variables added for dust model ------------
    real(r8), parameter :: ch = 0.5e-9_r8                     ! dimensional factor & tuning number,
                                                              ! as it's model resolution dependent (kgs^2/m^5)!!!
    real(r8)            :: r(NBIN)                            ! bin center (cm)
    real(r8)            :: uth                                ! threshold wind velocity (m/s)
    real(r8)            :: uv10                               ! 10 m wind speed (m/s)
    real(r8)            :: cd10                               ! 10-m drag coefficient ()
    real(r8)            :: wwd                                ! raw wind speed (m/s) 
    real(r8)            :: sp                                 ! mass fraction for soil factor
    integer             :: idustbin                           ! ibin to use for dust production, smallest silt bin for clay
    real(r8)            :: soilfact(pcols)                    ! soil erosion factor (for debug)

    ! ------------ local variables added for organics model ----------------------
    real(r8)     :: dr(carma%f_NBIN)
    real(r8)     :: rmass(carma%f_NBIN)
    real(r8)     :: SaltFlux(pcols)             ! sea salt flux to calculate marine POA

    ! Default return code.
    rc = RC_OK

    ! Determine the latitude and longitude of each column.
    lchnk = state%lchnk
    ncol = state%ncol

    call get_lat_all_p(lchnk, ncol, ilat)
    call get_lon_all_p(lchnk, ncol, ilon)

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8
    
    ! For emissions into the atmosphere, put the emission here.
    !
    ! NOTE: Do not set tendency to be the surface flux. Surface source is put in to
    ! the bottom layer by vertical diffusion. See vertical_solver module, line 355.            
    tendency(:ncol, :pver) = 0.0_r8
        
    call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname)
    if (RC < RC_ERROR) return
    
    call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname, r=r)
    if (RC < RC_ERROR) return
    
    if (ielem == I_ELEM_CRDUST) then
    ! st if (shortname .eq. "CRDUST") then  ! done by Pengfei
    
      ! Is this clay or silt?
      !
      ! NOTE: It is assumed that 90% of the mass will be silt and 10% will
      ! be clay.
      !
      ! NOTE: For clay bins, use the smallest silt bin to calculate the
      ! mass and then scale that into each clay bin based upon interpolation of
      ! Tegen and Lacis [1996].
      if (r(ibin) >= rClay) then
        sp         = 0.9_r8 / nSilt
        idustbin   = ibin
      else
        sp         = 0.1_r8 / nClay
        idustbin   = nClay + 1 
      end if

      ! Process each column.
      do icol = 1,ncol
      
        call CARMA_SurfaceWind(carma, state, icol, ilat(icol), ilon(icol), ielem, igroup, idustbin, cam_in, uv10, wwd, uth, rc) 

        ! Is the wind above the threshold for dust production?
        if (uv10 > uth) then
          surfaceFlux(icol) = ch * soil_factor(ilat(icol), ilon(icol)) * sp * &
                              wwd * (uv10 - uth)           
        endif
        
        ! Scale the clay bins based upon the smallest silt bin.   
        surfaceFlux(icol) = clay_mf(ibin) * surfaceFlux(icol)
        
        ! Save off the soil erosion factor so it can be output.
        soilfact(icol) = soil_factor(ilat(icol), ilon(icol))
      end do

      ! For debug purposes, output the soil erosion factor.
      call outfld('CRSLERFC', soilfact, pcols, lchnk)

    ! st elseif (sname_elem .eq. "CRSALT") then ! this is how Pengfei did
    elseif (ielem == I_ELEM_CRSALT) then
       !! introduce marine POA emission, use ChlorophyII-dependent mass contribution of OC
       !! see Gantt et al., 2009
       !! for sub-micron, I use sea salt flux instead of sub-micron marine particles
       !! needed to verify later
       !! Added by Pengfei Yu
       !! Oct.6.2012
       ! get sea spray aerosol flux first (for ibin; SaltFlux(:ncol) unit:kg/m2/s)

      call CARMA_SaltFlux(carma, ibin, state, r, dr, rmass, cam_in, SaltFlux, rc)

!st  not used currently  but done by Pengfei
       ! get [Chl-a] data
  !!   do icol = 1, ncol
  !!       if (Chla(lat(icol), lon(icol)) .lt. 0._r8) then
  !!          Fsub(icol) = 0._r8
  !!       else
  !!          Fsub(icol) = Chla(lat(icol), lon(icol)) * 0.63_r8 + 0.1_r8
  !!       endif
  !!       Fsub(icol) = min(Fsub(icol), 1._r8)
  !!   enddo
  !!   surfaceFlux(:ncol) = SaltFlux(:ncol)
  !!   ! sea salt (NaCl) flux should exclude marine organics and marine sulfate
  !!   if (carma%f_group(igroup)%f_r(ibin) .le. 0.5e-4_r8) then
  !!       !surfaceFlux(:ncol) = SaltFlux(:ncol)*(1._r8-0.0983_r8) - SaltFlux(:ncol) * Fsub(:ncol)
  !!        surfaceFlux(:ncol) = (SaltFlux(:ncol) - SaltFlux(:ncol)*Fsub(:ncol))/1.0983_r8
  !!   else
  !!       !surfaceFlux(:ncol) = SaltFlux(:ncol)*(1._r8-0.0983_r8) - SaltFlux(:ncol) * (Fsub(:ncol)*0.03_r8)
  !!        surfaceFlux(:ncol) = (SaltFlux(:ncol) - SaltFlux(:ncol)*Fsub(:ncol)*0.03_r8)/1.0983_r8
  !!   endif
       surfaceFlux(:ncol) = SaltFlux(:ncol)
    end if        
    
    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, rc)
    use cam_history,  only: addfld, add_default, horiz_only
    use constituents, only: pcnst

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    ! -------- local variables ----------
    integer            :: ibin                                ! CARMA bin index
    real(r8)           :: r(carma%f_NBIN)                     ! bin center (cm)
    integer            :: count_Silt                          ! count number for Silt
    integer            :: igroup                              ! the index of the carma aerosol group
    integer            :: ielem                               ! the index of the carma aerosol element
    character(len=32)  :: shortname                           ! the shortname of the element
    integer            :: LUNOPRT                             ! logical unit number for output
    logical            :: do_print                            ! do print output?
    
        
    ! Default return code.
    rc = RC_OK

    ! Determine how many clay and how many silt bins there are, based
    ! upon the bin definitions and rClay.
    !
    ! TBD: This should use the radii rather than being hard coded.
    ! nClay = 8
    ! nSilt = NBIN - nClay   
    do ielem = 1, NELEM    
       ! To get particle radius
       call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname)
       if (RC < RC_ERROR) return
       
       call CARMAGROUP_GET(carma, igroup, rc, r=r)
       if (RC < RC_ERROR) return
       
       if (shortname .eq. "CRDUST") then
          count_Silt = 0
          do ibin = 1, NBIN
             if (r(ibin) >= rclay) then
                count_Silt = count_Silt + 1
             else
             end if
          end do       
          nSilt = count_Silt
          nClay = NBIN - nSilt     
       end if       
    end do
    
    ! Read in the soil factors.
    call CARMA_ReadSoilErosionFactor(carma, rc)
    if (RC < RC_ERROR) return
    
    ! To determine Clay Mass Fraction
    do ielem = 1, NELEM    
       ! To get particle radius
       call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname)
       if (RC < RC_ERROR) return

       if (shortname .eq. "CRDUST") then
          call CARMA_ClayMassFraction(carma, igroup, rc) 
       end if       
    end do
    
    if (masterproc) then
      call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
      if (rc < 0) call endrun("CARMA_InitializeModel: CARMA_Get failed.")

      if (do_print) then
        write(carma%f_LUNOPRT,*) 'Initializing CARMA dust model ...'
        write(carma%f_LUNOPRT,*) 'nClay = ', nClay, ' nSilt = ', nSilt
        write(carma%f_LUNOPRT,*) 'clay_mf = ', clay_mf    
        write(carma%f_LUNOPRT,*) 'soil_factor = ', soil_factor
        
        write(carma%f_LUNOPRT,*) 'CARMA dust initialization complete'
      end if
    end if
    
    call addfld('CRSLERFC', horiz_only, 'A', 'fraction', 'CARMA soil erosion factor')
    
    return
  end subroutine CARMA_InitializeModel


  !! Sets the initial condition for CARMA aerosol particles. By default, there are no
  !! particles, but this routine can be overridden for models that wish to have an
  !! initial value.
  !!
  !! NOTE: If CARMA constituents appear in the initial condition file, then those
  !! values will override anything set here.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeParticle(carma, ielem, ibin, latvals, lonvals, mask, q, rc)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plat, plev, plon

    implicit none

    type(carma_type), intent(in)  :: carma      !! the carma object
    integer,          intent(in)  :: ielem      !! element index
    integer,          intent(in)  :: ibin       !! bin index
    real(r8),         intent(in)  :: latvals(:) !! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) !! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    !! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     !! mass mixing ratio (gcol, lev)
    integer,          intent(out) :: rc         !! return code, negative indicates failure

    ! Default return code.
    rc = RC_OK

    ! Add initial condition here.
    !
    ! NOTE: Initialized to 0. by the caller, so nothing needs to be done.

    return
  end subroutine CARMA_InitializeParticle
  

  !!  Called after wet deposition has been performed. Allows the specific model to add
  !!  wet deposition of CARMA aerosols to the aerosols being communicated to the surface.
  !!
  !!  @version July-2011 
  !!  @author  Chuck Bardeen 
  subroutine CARMA_WetDeposition(carma, ielem, ibin, sflx, cam_out, state, rc)
    use camsrfexch,       only: cam_out_t

    implicit none
    
    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: ielem       !! element index
    integer, intent(in)                  :: ibin        !! bin index
    real(r8), intent(in)                 :: sflx(pcols) !! surface flux (kg/m2/s)
    type(cam_out_t), intent(inout)       :: cam_out     !! cam output to surface models
    type(physics_state), intent(in)      :: state       !! physics state variables
    integer, intent(out)                 :: rc          !! return code, negative indicates failure
    
    integer    :: icol
 
    ! Default return code.
    rc = RC_OK
    
    return
  end subroutine CARMA_WetDeposition 


  !! Calculates the emissions for CARMA sea salt aerosol particles. 
  !!
  !! @author  Tianyi Fan, Chuck Bardeen, Pengfei Yu
  !! @version Dec-2010
  !! originally calculate sea salt flux in EmitParticle, Pengfei Yu make
  !! it a separate subroutine since multiple aerosol types need salt flux
  !! e.g. sea salt, sea salt sulfate, marine organics
  subroutine CARMA_SaltFlux(carma, ibin, state, r, dr, rmass, cam_in, SaltFlux, rc)
    use shr_kind_mod,  only: r8 => shr_kind_r8
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use phys_grid,     only: get_lon_all_p, get_lat_all_p, get_rlat_all_p
    use camsrfexch,    only: cam_in_t
    use cam_history,   only: outfld

    implicit none

    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ibin                  !! bin index
    type(physics_state), intent(in)    :: state                 !! physics state
    real(r8), intent(in)               :: r(NBIN)               !! bin center (cm)
    real(r8), intent(in)               :: dr(NBIN)              !! bin width (cm)
    real(r8), intent(in)               :: rmass(NBIN)           !! bin mass (g)    
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: SaltFlux(pcols)       !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer      :: ilat(pcols)              ! latitude index 
    integer      :: ilon(pcols)              ! longitude index
    integer      :: lchnk                   ! chunk identifier
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index


    ! -------- local variables added for sea salt model ------------
    real(r8)            :: rdrycm, rdry                       ! dry radius [cm], [um]     
    real(r8)            :: r80cm, r80                         ! wet radius at relatige humidity of 80% [cm]
    real(r8)            :: ncflx                              ! dF/dr [#/m2/s/um]
    real(r8)            :: Monahan, Clarke, Smith             ! dF/dr [#/m2/s/um]
    real(r8)            :: A_para, B_para, sita_para          ! A, B, and sita parameters in Gong
    real(r8)            :: B_mona                             ! the parameter used in Monahan
    real(r8)            :: W_Caff                             ! Correction factor in Caffrey
    real(r8)            :: u14, ustar_smith, cd_smith         ! 14m wind velocity, friction velocity, and drag coefficient as desired by Andreas source function              
    real(r8)            :: wcap                               ! whitecap coverage
    real(r8)            :: fref                               ! correction factor suggested by Hoppe2005
    real(r8), parameter :: xkar = 0.4_r8                      ! Von Karman constant
    real(r8)            :: u10in                              ! 10 meter wind speed use in the emission rate

    ! ------------------------------------------------------------------------------------------------
    ! -- Martensson source function. Coefficients for the parameterization of Ak(c4-c0) and Bk(d4-d0)
    ! -------------------------------------------------------------------------------------------------
    real(r8), parameter :: c41 = -2.576e35_r8
    real(r8), parameter :: c42 = -2.452e33_r8
    real(r8), parameter :: c43 =  1.085e29_r8
    real(r8), parameter :: c31 =  5.932e28_r8
    real(r8), parameter :: c32 =  2.404e27_r8
    real(r8), parameter :: c33 = -9.841e23_r8
    real(r8), parameter :: c21 = -2.867e21_r8
    real(r8), parameter :: c22 = -8.148e20_r8
    real(r8), parameter :: c23 =  3.132e18_r8
    real(r8), parameter :: c11 = -3.003e13_r8
    real(r8), parameter :: c12 =  1.183e14_r8
    real(r8), parameter :: c13 = -4.165e12_r8
    real(r8), parameter :: c01 = -2.881e6_r8
    real(r8), parameter :: c02 = -6.743e6_r8
    real(r8), parameter :: c03 =  2.181e6_r8
    real(r8), parameter :: d41 = 7.188e37_r8
    real(r8), parameter :: d42 = 7.368e35_r8
    real(r8), parameter :: d43 = -2.859e31_r8
    real(r8), parameter :: d31 =-1.616e31_r8
    real(r8), parameter :: d32 =-7.310e29_r8
    real(r8), parameter :: d33 = 2.601e26_r8
    real(r8), parameter :: d21 = 6.791e23_r8
    real(r8), parameter :: d22 = 2.528e23_r8
    real(r8), parameter :: d23 =-8.297e20_r8
    real(r8), parameter :: d11 = 1.829e16_r8
    real(r8), parameter :: d12 =-3.787e16_r8
    real(r8), parameter :: d13 = 1.105e15_r8
    real(r8), parameter :: d01 = 7.609e8_r8
    real(r8), parameter :: d02 = 2.279e9_r8
    real(r8), parameter :: d03 =-5.800e8_r8

    real(r8)            :: rpdry                              ! dry radius 
    real(r8)            :: Ak1                                ! Coefficient Ak in Martensson's source function
    real(r8)            :: Ak2
    real(r8)            :: Ak3
    real(r8)            :: Bk1                                ! Coefficient Bk in Martensson's source function
    real(r8)            :: Bk2
    real(r8)            :: Bk3
    Ak1(rpdry)= c41*(2._r8*rpdry)**4 + c31*(2._r8*rpdry) ** 3 + c21*(2._r8*rpdry)**2 + c11*(2._r8*rpdry)+ c01
    Ak2(rpdry)= c42*(2._r8*rpdry)**4 + c32*(2._r8*rpdry) ** 3 + c22*(2._r8*rpdry)**2 + c12*(2._r8*rpdry)+ c02
    Ak3(rpdry)= c43*(2._r8*rpdry)**4 + c33*(2._r8*rpdry) ** 3 + c23*(2._r8*rpdry)**2 + c13*(2._r8*rpdry)+ c03
    Bk1(rpdry)= d41*(2._r8*rpdry)**4 + d31*(2._r8*rpdry) ** 3 + d21*(2._r8*rpdry)**2 + d11*(2._r8*rpdry)+ d01
    Bk2(rpdry)= d42*(2._r8*rpdry)**4 + d32*(2._r8*rpdry) ** 3 + d22*(2._r8*rpdry)**2 + d12*(2._r8*rpdry)+ d02
    Bk3(rpdry)= d43*(2._r8*rpdry)**4 + d33*(2._r8*rpdry) ** 3 + d23*(2._r8*rpdry)**2 + d13*(2._r8*rpdry)+ d03

   ! ------------------------------------------------------------
    ! ----  Clarke Source Function. Coefficients for Ai    -------
    ! ------------------------------------------------------------
    real(r8), parameter :: beta01 =-5.001e3_r8
    real(r8), parameter :: beta11 = 0.808e6_r8
    real(r8), parameter :: beta21 =-1.980e7_r8
    real(r8), parameter :: beta31 = 2.188e8_r8
    real(r8), parameter :: beta41 =-1.144e9_r8
    real(r8), parameter :: beta51 = 2.290e9_r8
    real(r8), parameter :: beta02 = 3.854e3_r8
    real(r8), parameter :: beta12 = 1.168e4_r8
    real(r8), parameter :: beta22 =-6.572e4_r8
    real(r8), parameter :: beta32 = 1.003e5_r8
    real(r8), parameter :: beta42 =-6.407e4_r8
    real(r8), parameter :: beta52 = 1.493e4_r8
    real(r8), parameter :: beta03 = 4.498e2_r8
    real(r8), parameter :: beta13 = 0.839e3_r8
    real(r8), parameter :: beta23 =-5.394e2_r8
    real(r8), parameter :: beta33 = 1.218e2_r8
    real(r8), parameter :: beta43 =-1.213e1_r8
    real(r8), parameter :: beta53 = 4.514e-1_r8
    real(r8)            :: A1                                ! Coefficient Ak in Clarkes's source function
    real(r8)            :: A2
    real(r8)            :: A3
    A1(rpdry) = beta01 + beta11*(2._r8*rpdry) + beta21*(2._r8*rpdry)**2 + beta31*(2._r8*rpdry)**3 + beta41*(2._r8*rpdry)**4 + beta51*(2._r8*rpdry)**5
    A2(rpdry) = beta02 + beta12*(2._r8*rpdry) + beta22*(2._r8*rpdry)**2 + beta32*(2._r8*rpdry)**3 + beta42*(2._r8*rpdry)**4 + beta52*(2._r8*rpdry)**5
    A3(rpdry) = beta03 + beta13*(2._r8*rpdry) + beta23*(2._r8*rpdry)**2 + beta33*(2._r8*rpdry)**3 + beta43*(2._r8*rpdry)**4 + beta53*(2._r8*rpdry)**5

    ! ---------------------------------------------
    ! coefficient A1, A2 in Andreas's Source funcion
    ! ---------------------------------------------
    real(r8)            ::A1A92
    real(r8)            ::A2A92

    ! ---------------------------------------------
    ! coefficient in Smith's Source funcion
    ! --------------------------------------------- 
    real(r8), parameter ::  f1 = 3.1_r8
    real(r8), parameter ::  f2 = 3.3_r8
    real(r8), parameter ::  r1 = 2.1_r8
    real(r8), parameter ::  r2 = 9.2_r8
    real(r8), parameter ::  delta = 10._r8

    ! --------------------------------------------------------------------
    ! ---- constants in calculating the particle wet radius [Gerber, 1985]   
    ! --------------------------------------------------------------------        
    real(r8), parameter :: c1   = 0.7674_r8        ! .
    real(r8), parameter :: c2   = 3.079_r8         ! .
    real(r8), parameter :: c3   = 2.573e-11_r8     ! .
    real(r8), parameter :: c4   = -1.424_r8        ! constants in calculating the particel wet radius

    ! Default return code.
    rc = RC_OK

    ! Determine the latitude and longitude of each column.
    lchnk = state%lchnk
    ncol = state%ncol

    call get_lat_all_p(lchnk, ncol, ilat)
    call get_lon_all_p(lchnk, ncol, ilon)

    ! Add any surface flux here.
    SaltFlux(:ncol) = 0.0_r8

   ! Are we configured for one of the known emission schemes?
      if(carma_seasalt_emis .ne. "Gong"       .and. &
         carma_seasalt_emis .ne. "Martensson" .and. &
         carma_seasalt_emis .ne. "Clarke"     .and. &
         carma_seasalt_emis .ne. "Andreas"    .and. &
         carma_seasalt_emis .ne. "Caffrey"    .and. &
         carma_seasalt_emis .ne. "CMS"        .and. &
         carma_seasalt_emis .ne. "NONE"       .and. &
         carma_seasalt_emis .ne. "CONST"        ) then

        call endrun('carma_EmitParticle:: Invalid sea salt emission scheme.')
      end if

      !**********************************
      ! wet sea salt radius at RH = 80%  
      !**********************************    
      r80cm   = (c1 *  (r(ibin)) ** c2 / (c3 * r(ibin) ** c4 - log10(0.8_r8)) + (r(ibin))**3) ** (1._r8/3._r8) ! [cm]
      rdrycm  = r(ibin)  ! [cm]   
      r80     = r80cm *1.e4_r8    ! [um]
      rdry    = rdrycm*1.e4_r8  ! [um]
     
       do icol = 1,ncol

        ! Only generate sea salt over the ocean.
        if (cam_in%ocnfrac(icol) > 0._r8) then

          !**********************************
          !    WIND for seasalt production
          !**********************************    
          call CARMA_SurfaceWind_salt(carma, state, icol, ilat(icol), ilon(icol), cam_in, u10in, rc)

          ! Add any surface flux here.       
          ncflx       = 0.0_r8
          Monahan     = 0.0_r8
          Clarke      = 0.0_r8
          Smith       = 0.0_r8

          !**********************************
          !        Whitecap Coverage
          !**********************************
          wcap = 3.84e-6_r8 * u10in ** 3.41_r8      ! in percent, ie., 75%, wcap = 0.75

          !****************************************
          !        Hoppel correction factor
          !        Smith drag coefficients and etc
          !****************************************
          if (u10in .le. 10._r8) then
            cd_smith = 1.14e-3_r8
          else
            cd_smith = (0.49_r8 + 0.065_r8 * u10in) * 1.e-3_r8
          end if

          ustar_smith = cd_smith **0.5_r8 * u10in

          ! We don't have vg yet, since that is calculated by CARMA. That will require
          ! a different interface for the emissions, storing vg in the physics buffer,
          ! and/or doing some duplicate calculations for vg assuming 80% RH.
!          fref = (delta/state%zm(icol, pver))**(vg(icol, ibin, igelem(i))/(xkar*ustar_smith))
          fref = 1.0_r8

          !**********************************
          !        Source Functions
          !**********************************
          if (carma_seasalt_emis .eq. 'NONE') then
            ncflx = 0._r8
          end if

          if (carma_seasalt_emis .eq. 'CONST') then
            ncflx = 1.e-5_r8
          end if
     
         !-------Gong source function------
          if (carma_seasalt_emis == "Gong") then
            sita_para = 30
            A_para = - 4.7_r8 * (1+ sita_para * r80) ** (- 0.017_r8 * r80** (-1.44_r8))
            B_para = (0.433_r8 - log10(r80)) / 0.433_r8
            ncflx = 1.373_r8* u10in ** 3.41_r8 * r80 ** A_para * (1._r8 + 0.057_r8 * r80**3.45_r8) * 10._r8 ** (1.607_r8 * exp(- B_para **2))
!            if (do_print) write(LUNOPRT, *) "Gong: ncflx = ", ncflx, ", u10n = ", u10in
          end if

          !------Martensson source function-----
          if (carma_seasalt_emis == "Martensson") then
            if (rdry .le. 0.0725_r8) then
              ncflx = (Ak1(rdry*1.0e-6_r8)* (25._r8+273._r8) + Bk1(rdry*1.0e-6_r8)) * wcap      ! dF/dlogr [#/s/m2]                                     
              ncflx = ncflx / (2.30258509_r8 * rdry)                                            ! dF/dr    [#/s/m2/um]        
            elseif (rdry .gt. 0.0725_r8 .and. rdry .le. 0.2095_r8) then
              ncflx = (Ak2(rdry*1.0e-6_r8)* (25._r8+273._r8) + Bk2(rdry*1.0e-6_r8)) * wcap      ! dF/dlogr [#/s/m2]
              ncflx = ncflx / (2.30258509_r8 * rdry)                                            ! dF/dr    [#/s/m2/um]                                
            elseif (rdry .gt. 0.2095_r8 .and. rdry .le. 1.4_r8) then
              ncflx = (Ak3(rdry*1.0e-6_r8)* (25._r8+273._r8) + Bk3(rdry*1.0e-6_r8)) * wcap      ! dF/dlogr [#/s/m2]
              ncflx = ncflx / (2.30258509_r8 * rdry)                                            ! dF/dr    [#/s/m2/um] 
            else
              ncflx = 0._r8
            end if
          end if

          !-------Clarke source function------- 
          if (carma_seasalt_emis == "Clarke")then
            if (rdry .lt. 0.066_r8) then
             ncflx = A1(rdry) * 1.e4_r8 * wcap                              ! dF/dlogr [#/s/m2] 
              ncflx = ncflx / (2.30258509_r8 * rdry)                        ! dF/dr    [#/s/m2/um]
            elseif (rdry .ge. 0.066_r8 .and. rdry .lt. 0.6_r8) then
              ncflx = A2(rdry) * 1.e4_r8 * wcap                             ! dF/dlogr [#/s/m2] 
              ncflx = ncflx / (2.30258509_r8 * rdry)                        ! dF/dr    [#/s/m2/um]
            elseif (rdry .ge. 0.6_r8 .and. rdry .lt. 4.0_r8) then
              ncflx = A3(rdry) * 1.e4_r8 * wcap                             ! dF/dlogr [#/s/m2] 
              ncflx= ncflx / (2.30258509_r8 * rdry)                         ! dF/dr    [#/s/m2/um]
            else
              ncflx = 0._r8
            end if
          end if

          !-----------Caffrey source function------------
          if (carma_seasalt_emis == "Caffrey") then

            !Monahan        
            B_mona = (0.38_r8 - log10(r80)) / 0.65_r8
            Monahan = 1.373_r8 * (u10in**3.41_r8) * r80**(-3._r8) * (1._r8 + 0.057 *r80**1.05_r8)  * 10._r8 ** (1.19_r8 * exp(-1. * B_mona**2)) ! dF/dr

            !Smith
            u14 = u10in * (1._r8 + cd_smith**0.5_r8 / xkar * log(14._r8 / 10._r8))  ! 14 meter wind
            A1A92 = 10._r8 ** (0.0676_r8 * u14 + 2.430_r8)
            A2A92 = 10._r8 ** (0.9590_r8 * u14**0.5_r8 - 1.476_r8)
            Smith = A1A92*exp(-f1 *(log(r80/r1))**2) + A2A92*exp(-f2 * (log(r80/r2))**2)     ! dF/dr   [#/m2/s/um]                    

            !Caffrey based on Monahan and Smith
            W_Caff = 1.136_r8 **(-1._r8 * rdry ** (-0.855_r8))*(1._r8 + 0.2_r8/rdry)
            if (rdry .lt. 0.15_r8) then
              ncflx = Monahan
            else
              if (u10in .le. 9._r8) then
                ncflx = Monahan
              else
                if(Monahan .ge. Smith) then
                  ncflx = Monahan
                else
                  ncflx = Smith
                end if
              end if
            end if

            ncflx = ncflx * W_Caff

            !%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Apply Hoppel correction
            !%%%%%%%%%%%%%%%%%%%%%%%%%
            ncflx = ncflx * fref
          end if

          !--------CMS (Clarke, Monahan, and Smith source function)-------
          if (carma_seasalt_emis == "CMS") then

            !Clarke  
            if (rdry .lt. 0.066_r8) then
              Clarke = A1(rdry) * 1.e4_r8 * wcap                     ! dF/dlogr [#/s/m2] 
              Clarke = Clarke / (2.30258509_r8 * rdry)               ! dF/dr    [#/s/m2/um]
            elseif ((rdry .ge. 0.066_r8) .and. (rdry .lt. 0.6_r8)) then
              Clarke = A2(rdry) * 1.e4_r8 * wcap                     ! dF/dlogr [#/s/m2] 
              Clarke = Clarke / (2.30258509_r8 * rdry)               ! dF/dr    [#/s/m2/um]
            elseif ((rdry .ge. 0.6_r8) .and. (rdry .lt. 4.0_r8)) then
              Clarke = A3(rdry) * 1.e4_r8 * wcap                      ! dF/dlogr [#/s/m2] 
              Clarke= Clarke / (2.30258509_r8 * rdry)                 ! dF/dr    [#/s/m2/um]
            end if

          !Monahan   
            B_Mona = (0.38_r8 - log10(r80)) / 0.65_r8
            Monahan = 1.373_r8 * u10in ** 3.41_r8 * r80 ** (-3._r8) * (1._r8 + 0.057_r8 * r80**1.05_r8) * 10._r8 ** (1.19_r8 * exp(- B_Mona **2))

            !Smith
            u14 = u10in * (1._r8 + cd_smith**0.5_r8 / xkar*log(14._r8 / 10._r8))  ! 14 meter wind
            A1A92 = 10._r8 ** (0.0676_r8 * u14 + 2.430_r8)
            A2A92 = 10._r8 ** (0.9590_r8 * u14**0.5_r8 - 1.476_r8)
            Smith = A1A92*exp(-f1 *(log(r80 / r1))**2) + A2A92*exp(-f2 * (log(r80 / r2))**2)     ! dF/dr   [#/m2/s/um]                    

            !%%%%%%%%%%%%%%%%%%%%%%%%%
            !     CMS1 or CMS2
            !%%%%%%%%%%%%%%%%%%%%%%%%%
  !          if (rdry .lt. 0.1_r8) then   ! originally cut at 0.1 um
            ! ***CMS1*****
            if (rdry .lt. 1._r8) then    ! cut at 1.0 um
            ! ***CMS2*****
  !          if (rdry .lt. 2._r8) then    ! cut at 2.0 um
              ncflx = Clarke
            else
              if (u10in .lt. 9._r8) then
                ncflx = Monahan
              else
                if (Monahan .gt. Smith) then
                  ncflx = Monahan
                else
                  ncflx = Smith
                end if
              end if
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Apply Hoppel correction
            !%%%%%%%%%%%%%%%%%%%%%%%%%       
            ncflx = ncflx * fref
          end if
    
          ! convert ncflx [#/m^2/s/um] to surfaceFlx [kg/m^2/s]              
          ! SaltFlux(icol) = ncflx * dr(ibin) * rmass(ibin) * 10._r8      ! *1e4[um/cm] * 1.e-3[kg/g]
                  SaltFlux(icol) = ncflx * dr(ibin) * rmass(ibin) * 10._r8
                        
!          if (do_print) write(LUNOPRT, *) "ibin = ", ibin, ", igroup = ", igroup
!          if (do_print) write(LUNOPRT, *) "dr = ", dr(ibin), ", rmass = ", rmass(ibin)
!          if (do_print) write(LUNOPRT, *) "ncflx = " , ncflx, ", SaltFlux = ", SaltFlux(icol)

          ! weighted by the ocean fraction
          SaltFlux(icol) = SaltFlux(icol) * cam_in%ocnfrac(icol)
        end if
      end do

    return
  end subroutine CARMA_SaltFlux
 
  
  !! Calculate the sea surface wind with a Weibull distribution.
  !!
  !! @author  Tianyi Fan
  !! @version August-2010
  subroutine CARMA_SurfaceWind_salt(carma, state, icol, ilat, ilon, cam_in, u10in, rc)
    use ppgrid,           only: pcols, pver
    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t

    implicit none

    ! in and out field
    type(carma_type), intent(in)        :: carma                 !! the carma object
    type(physics_state), intent(in)     :: state                 !! physics state
    integer, intent(in)                 :: icol                  !! column index
    integer, intent(in)                 :: ilat                  !! latitude index 
    integer, intent(in)                 :: ilon                  !! longitude index
    type(cam_in_t), intent(in)          :: cam_in                !! surface inputs
    real(r8), intent(out)               :: u10in                 !! the 10m wind speed put into the source function
    integer, intent(out)                :: rc                    !! return code, negative indicates failure

    ! local variables
    real(r8) :: uWB341              ! the nth mean wind with integration using Weibull Distribution(integrate from threshold wind velocity)      

    rc = RC_OK

    uWB341 = 0._r8

    ! calc. the Weibull wind distribution         
    u10in = cam_in%u10(icol)

    ! Use Weibull with prescribed coefficients?
    if (carma_do_WeibullK) then
      call WeibullWind(u10in, uth_salt, 3.41_r8, uWB341, Weibull_k(ilat, ilon))
    else
      call WeibullWind(u10in, uth_salt, 3.41_r8, uWB341)
    end if

    u10in = uWB341 ** (1._r8 / 3.41_r8)

!    if (do_print) write(LUNOPRT, *) 'CARMA_SurfaceWind: icol ',icol, ', u10 =', cam_in%u10(icol), ', u10in =', u10in

    return
  end subroutine CARMA_SurfaceWind_salt



  !!  Determines the mass fraction for the clay (submicron) bins based upon
  !!  Tegen and Lacis [1996]. The total fraction for all clay bins should
  !!  add up to 1.
  !!
  !!  NOTE: WOuld it be better to interpolate this into the bins rather than
  !!  assigning all CARMA bins within a Tegen & Lacis bin the same value?
  !!
  !!  NOTE: Should any mass go to bins smaller than the smallest one used by
  !!  Tegen and Lacis?
  !!
  !!  @version July-2012 
  !!  @author  Lin Su, Pengfei Yu, Chuck Bardeen 
  subroutine CARMA_ClayMassFraction(carma, igroup, rc)
    implicit none
    
    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: igroup      !! the carma group index
    integer, intent(inout)               :: rc          !! return code, negative indicates failure

    ! Bins and mass fraction from Tegen and Lacis.
    integer, parameter  :: NBIN_TEGEN = 4    
    real(r8)            :: tl_rmin(NBIN_TEGEN) = (/ 1.e-5_r8,  1.8e-5_r8, 3.e-5_r8, 6.e-5_r8 /)
    real(r8)            :: tl_rmax(NBIN_TEGEN) = (/ 1.8e-5_r8, 3.e-5_r8,  6.e-5_r8, 1.e-4_r8 /)
    real(r8)            :: tl_mf(NBIN_TEGEN)   = (/ 0.009_r8,  0.081_r8,  0.234_r8, 0.676_r8 /)

    ! Local Variables
    integer, parameter  :: IBELOW = 1    
    integer, parameter  :: IABOVE = 6    
    integer             :: tl_count(NBIN_TEGEN+2)  ! count number in Tegen and Lacis ranges
    integer             :: ind_up(NBIN_TEGEN+2)
    integer             :: ind_low(NBIN_TEGEN+2)
    integer             :: j                    ! local index number
    integer             :: ibin                 ! carma bin index
    real(r8)            :: r(carma%f_NBIN)      ! CARMA bin center (cm)
     
    ! Default return code.
    rc = RC_OK
    
    ! Interpolate from Tegen and Lacis.
    call CARMAGROUP_GET(carma, igroup, rc, r=r)
    if (RC < RC_ERROR) return
    
    ! Figure out how many of the CARMA bins are in each of the Tegen and Lacis
    ! ranges.
    tl_count(:) = 0
    
    do ibin = 1, NBIN
    
      ! Smaller than the range.
      if (r(ibin) < tl_rmin(1)) then
        tl_count(IBELOW) = tl_count(IBELOW) + 1
      end if
      
      ! In the range
      do j = 1, NBIN_TEGEN
        if (r(ibin) < tl_rmax(j) .and. r(ibin) >= tl_rmin(j)) then
          tl_count(j+1) = tl_count(j+1) + 1
        end if
      end do

      ! Bigger than the range.
      if (r(ibin) >= tl_rmax(NBIN_TEGEN)) then
        tl_count(IABOVE) = tl_count(IABOVE) + 1
      end if       
    end do

    ! Determine where the boundaries are between the TEGEN bins and
    ! the CARMA bin structure.
    ind_up(:)   = 0
    ind_low(:)  = 0
    ind_up (IBELOW)  = tl_count(IBELOW)
    ind_low(IBELOW)  = min(1, tl_count(IBELOW))
    
    do j = 1, 5
      ind_up (j+1) = ind_up(j) + tl_count(j+1)
      ind_low(j+1) = ind_up(j) + min(tl_count(j+1), 1)
    end do
    
    ! No mass to bins smaller than the smallest size.
    clay_mf(:) = 0._r8
    
    ! NOTE: This won't work right if the dust bins are coarser than
    ! the Tegen and Lacis bins. In this case mass fraction would need
    ! to be combined from the Tegen & Lacis bins into a CARMA bin. 
    do j = 1, NBIN_TEGEN
      if (tl_count(j+1) > 0) then
        clay_mf(ind_low(j+1):ind_up(j+1)) = tl_mf(j) / tl_count(j+1)
      end if
    end do 
                                       
    clay_mf(ind_low(IABOVE):) = 1._r8

    return
  end subroutine CARMA_ClayMassFraction

                                              
  !! Calculate the sea surface wind with a Weibull distribution.
  !!
  !! NOTE: This should be combined with a similar routine in the sea salt
  !! model, and any differences should be control by parameters into this
  !! routine (and perhaps namelist variables).
  !!
  !! @author  Lin Su, Pengfei Yu, Chuck Bardeen
  !! @version July-2012
  subroutine CARMA_SurfaceWind(carma, state, icol, ilat, ilon, ielem, igroup, ibin, cam_in, uv10, wwd, uth, rc)
    use ppgrid,           only: pcols, pver
    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t
 
    implicit none

    ! in and out field
    type(carma_type), intent(in)        :: carma                 !! the carma object
    type(physics_state), intent(in)     :: state                 !! physics state   
    integer, intent(in)                 :: icol                  !! column index
    integer, intent(in)                 :: ilat                  !! latitude index 
    integer, intent(in)                 :: ilon                  !! longitude index
    integer, intent(in)                 :: ielem                 !! element index
    integer, intent(in)                 :: igroup                !! group index
    integer, intent(in)                 :: ibin                  !! bin index
    type(cam_in_t), intent(in)          :: cam_in                !! surface inputs
    real(r8), intent(out)               :: uv10                  !! the 10m wind speed (m/s)
    real(r8), intent(out)               :: wwd                   !! the 10m wind speed  with Weibull applied (m/s)
    real(r8), intent(out)               :: uth                   !! the 10m wind threshold (m/s)
    integer,  intent(inout)             :: rc                    !! return code, negative indicates failure

    real(r8), parameter                 :: vk = 0.4_r8           ! von Karman constant
    real(r8)                            :: r(NBIN)               ! CARMA bin center (cm)
    real(r8)                            :: rhop(NBIN)            ! CARMA partile element density (g/cm3)
    real(r8)                            :: uthfact               !     
    integer                             :: iepart                ! element in group containing the particle concentration
    real(r8), parameter                 :: rhoa = 1.25e-3_r8     ! Air density at surface
    
    rc = RC_OK
    
    ! Get the 10 meter wind speed
    uv10 = cam_in%u10(icol)

    ! Calculate the threshold wind speed of each bin [Marticorena and Bergametti,1995]
    ! note that in cgs units --> m/s
    call CARMAGROUP_GET(carma, igroup, rc, r=r)
    if (RC < RC_ERROR) return
    
    ! Define particle # concentration element index for current group
    call CARMAELEMENT_Get(carma, ielem, rc, rho=rhop)
    if (RC < RC_ERROR) return
        
    if (cam_in%soilw(icol) > 0._r8 .AND. cam_in%soilw(icol) < 0.5_r8) then
       uthfact = 1.2_r8 + 0.2_r8*log10(cam_in%soilw(icol))
       if (r(ibin) > 2.825e-5_r8) then  ! r(4) = 2.825e-5 cm
           uth = uthfact * 1.e-2_r8 * 0.13_r8 * sqrt(rhop(ibin)*GRAV*r(ibin)*2._r8/rhoa) &
                       * sqrt(1._r8 + .006_r8/rhop(ibin)/GRAV/(r(ibin)*2._r8)**2.5_r8) &
                       / sqrt(1.928_r8*(1331._r8*(r(ibin)*2._r8)**1.56_r8 + .38_r8)**.092_r8 - 1._r8)
       else
           uth = uthfact*1.e-2_r8* 0.13_r8 * sqrt(rhop(ibin)*GRAV*(.75e-4_r8)*2./rhoa)   &
                       * sqrt(1._r8 + .006_r8/rhop(ibin)/GRAV/((.75e-4_r8)*2._r8)**2.5_r8) &
                       / sqrt(1.928_r8*(1331._r8*((.75e-4_r8)*2._r8)**1.56_r8 + .38_r8)**.092_r8 - 1._r8)
       endif
    else
       uth = uv10
    endif    

    ! Use Weibull with Lansing's estimate for shape.
    call WeibullWind(uv10, uth, 2._r8, wwd)

    return
  end subroutine CARMA_SurfaceWind


  !! Read in the dust source (soil) erodibility factor from a NETCDF file. In this
  !! processes, the data is regridded from the source size to the size needed by the
  !! model.
  !!
  !! NOTE: This is currently doing 2-D interpolation, but it really should be doing
  !! regridding.
  !!
  !! @author  Pengfei Yu
  !! @version July-2012

!! st
!! could use /components/cam/src/chemistry/aerosol/soil_erod_mod.F90 here insted of this routine?
  subroutine CARMA_ReadSoilErosionFactor(carma, rc)
    use pmgrid,        only: plat, plon
    use ioFileMod,     only: getfil
    use wrap_nf
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish    
    
    implicit none

    type(carma_type), intent(in)              :: carma                 !! the carma object
    integer, intent(out)                      :: rc                    !! return code, negative indicates failure

    ! local variables
    integer                                   :: idvar, f_nlon, f_nlat, idlat, idlon
    integer                                   :: fid, fid_lon, fid_lat
    real(r8), allocatable, dimension(:,:)     :: ero_factor, ero_factor1
    character(len=256)                        :: ero_file
    real(r8), allocatable, dimension(:)       :: ero_lat               ! latitude dimension
    real(r8), allocatable, dimension(:)       :: ero_lon               ! latitude dimension
    type (interp_type)                        :: wgt1, wgt2
    real(r8)                                  :: lat(plat), lon(plon)
    integer                                   :: i

    rc = RC_OK

    ! Open the netcdf file (read only)
    call getfil(carma_soilerosion_file, ero_file, 0)
    call wrap_open(ero_file, 0, fid)
  
    ! Get file dimensions
    call wrap_inq_dimid(fid, 'plon', fid_lon)
    call wrap_inq_dimid(fid, 'plat', fid_lat)
    call wrap_inq_dimlen(fid, fid_lon, f_nlon)
    call wrap_inq_dimlen(fid, fid_lat, f_nlat)
  
    allocate(ero_lat(f_nlat))
    allocate(ero_lon(f_nlon))
    allocate(ero_factor (f_nlon, f_nlat))
    allocate(ero_factor1(plon, plat))
    allocate(soil_factor(plat, plon))
    
    ! Read in the tables.
    call wrap_inq_varid(fid, 'new_source', idvar)
    i = nf90_get_var (fid, idvar, ero_factor)
    if (i/=NF90_NOERR) then
       write(iulog,*)'CARMA_ReadSoilErosionFactor: error reading varid =', idvar
       call handle_error (i)
    end if
    call wrap_inq_varid(fid, 'plat', idlat)
    call wrap_get_var_realx(fid, idlat,  ero_lat)
    call wrap_inq_varid(fid, 'plon', idlon)
    call wrap_get_var_realx(fid, idlon,  ero_lon)
            
    ! Close the file.
    call wrap_close(fid)
    
    ! NOTE: Is there a better way to get all of the dimensions
    ! needed for the model grid? Seems like it shouldn't be hard
    ! coded here.
    do i = 1, plat
       lat(i) = 180._r8 / (plat-1) * (i-1) - 90._r8
    end do
    
    do i = 1, plon
       lon(i) = 360._r8 / plon * (i-1)
    end do
    
    call lininterp_init(ero_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(ero_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(ero_factor, f_nlon, f_nlat, ero_factor1, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)
    
    soil_factor(:plat, :plon) = transpose(ero_factor1(:plon, :plat))
    
    deallocate(ero_lat)
    deallocate(ero_lon)
    deallocate(ero_factor)
    deallocate(ero_factor1)
    
    return
  end subroutine CARMA_ReadSoilErosionFactor


  !! Calculate the nth mean of u using Weibull wind distribution
  !! considering the threshold wind velocity. This algorithm
  !! integrates from uth to infinite (u^n P(u)du )
  !!  
  !! @author  Tianyi Fan
  !! @version August-2010
   subroutine WeibullWind(u, uth, n, uwb, wbk)
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use shr_spfn_mod, only: gamma =>  shr_spfn_gamma, &
         igamma => shr_spfn_igamma

    implicit none
  
    real(r8), intent(in)  :: u      ! mean wind speed
    real(r8), intent(in)  :: uth    ! threshold velocity
    real(r8), intent(in)  :: n      ! the rank of u in the integration
    real(r8), intent(out) :: uwb    ! the Weibull distribution
    real(r8), intent(in), optional ::  wbk    ! the shape parameter
  
    ! local variable
    real(r8)  :: k                  ! the shape parameter in Weibull distribution
    real(r8)  :: c                  ! the scale parameter in Weibull distribution
  
    if (present(wbk)) then
      k = wbk
    else
      k = 0.94*u**0.5_r8            ! follow Grini and Zender, 2004JGR
 !    k = 2.5_r8                   ! Lansing's estimate
    end if 
  
    ! If u is 0, then k can be 0, which makes a lot of this undefined.
    ! Just return 0. in this case.
    if (u == 0._r8) then
      uwb = 0._r8
    else 
      c   = u * (gamma(1._r8 + 1._r8 / k))**(-1._r8)  
      uwb = c**n * igamma(n / k + 1._r8, (uth / c)**k)
    end if

  end subroutine WeibullWind
  
end module
