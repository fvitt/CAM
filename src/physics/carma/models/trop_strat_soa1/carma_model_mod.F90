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
!!   - WeibullWind()
!!
!! @version April-2020
!! @author  Simone Tilmes, Lin Su, Pengfei Yu, Chuck Bardeen
!!  changes to pervious version: rename PURSULF to PRSULF to be easier read in in CAM

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
  use cam_abortutils, only: endrun
  use physics_types,  only: physics_state, physics_ptend
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc, pbuf_set_field, pbuf_get_field, pbuf_get_index
  use time_manager,   only: is_first_step
  use cam_logfile,    only: iulog

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
  integer, public, parameter      :: NGROUP   = 2               !! Number of particle groups
  integer, public, parameter      :: NELEM    = 7               !! Number of particle elements
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
  integer, public, parameter      :: I_SOA            = 3       !! SOA composition
  integer, public, parameter      :: I_BC             = 4       !! BC composition
  integer, public, parameter      :: I_DUST           = 5       !! dust composition
  integer, public, parameter      :: I_SALT           = 6       !! sea salt composition

  integer, public, parameter      :: I_GRP_PRSUL     = 1        !! sulfate aerosol
  integer, public, parameter      :: I_GRP_MXAER     = 2        !! mixed aerosol

  integer, public, parameter      :: I_ELEM_PRSUL     = 1       !! sulfate aerosol;  nameing needs to only have 2 charaters  before the element name to work with
                                                                !! partsof the code reading different elements
  integer, public, parameter      :: I_ELEM_MXAER     = 2       !! aerosol
  integer, public, parameter      :: I_ELEM_MXOC      = 3       !! organics aerosol
  integer, public, parameter      :: I_ELEM_MXSOA     = 4       !! secondary organic aerosol
  integer, public, parameter      :: I_ELEM_MXBC      = 5       !! black carbon
  integer, public, parameter      :: I_ELEM_MXDUST    = 6       !! dust aerosol
  integer, public, parameter      :: I_ELEM_MXSALT    = 7       !! sea salt aerosol

  integer, public, parameter      :: I_GAS_H2O        = 1       !! water vapor
  integer, public, parameter      :: I_GAS_H2SO4      = 2       !! sulphuric acid

  real(kind=f), public, parameter         :: Kappa_OC = 0.5_f      !! hygroscopicity of OC
  real(kind=f), public, parameter         :: Kappa_SOA = 0.5_f     !! hygroscopicity of SOA
  real(kind=f), public, parameter         :: Kappa_BC = 0.1_f
  real(kind=f), public, parameter         :: Kappa_DUST = 0.2_f
  real(kind=f), public, parameter         :: Kappa_SALT = 1.0_f
  real(kind=f), public, parameter         :: Kappa_SULF = 0.5_f

  real(kind=f), public, parameter         :: RHO_obc  = 1.35_f          !! dry density of smoke aerosol
  real(kind=f), public, parameter         :: RHO_DUST = 2.65_f          !! dry density of dust particles (g/cm^3) -Lin Su
  real(kind=f), public, parameter         :: RHO_SALT = 2.65_f          !! dry density of sea salt particles (g/cm)
  real(kind=f), public, parameter         :: RHO_SULFATE  = 1.923_f     !! dry density of sulfate particles (g/cm3)

 ! see CARMA_SmokeEmissionRead
! real(kind=f), allocatable, dimension(:,:)     ::   Chla                                       ! Chlorophy11 data (mg/m3)
  real(r8), allocatable, dimension(:,:,:)       ::   BCnew                              ! #/cm2/s
  real(r8), allocatable, dimension(:,:,:)       ::   OCnew


  ! for sea salt flux calculation
  real(r8), parameter             :: uth_salt = 4._r8                !! threshold wind velocity


  ! for dust calculation
  real(kind=f), parameter         :: rClay = 1e-4_f         !! silt/clay particle radius boundary (cm)

  integer                         :: nClay                  !! Number of clay bins (r < 1 um)
  integer                         :: nSilt                  !! Number of silt bins
  real(kind=f)                    :: clay_mf(NBIN)=-huge(1._f) !! clay mass fraction (fraction)
  real(kind=f), allocatable, dimension(:,:) :: soil_factor  !! Soil Erosion Factor (fraction)
  real(kind=f), public, parameter :: WTMOL_H2SO4    = 98.078479_f    !! molecular weight of sulphuric acid

! NOTE: The WeibullK distribution is not currently supported, since the coefficients are not
! generated. This can be added later.
!  real(r8), allocatable, dimension(:,:) :: Weibull_k            ! Weibull K(nlat,nlon
  real(kind=f), public, parameter     :: rmin_PRSUL     = 3.43e-8_f  ! minimum radius (cm)
  real(kind=f), public, parameter     :: vmrat_PRSUL    = 3.67_f     ! volume ratio
  real(kind=f), public, parameter     :: rmin_MXAER     = 5e-6_f     ! minimum radius (cm)
  real(kind=f), public, parameter     :: vmrat_MXAER    = 2.2588_f    !2.4610_f        ! volume ratio

! Physics buffer index for sulfate surface area density
  integer      :: ipbuf4sadsulf = -1 ! total aerosol surface area density over all groups (cm2/cm3)
  integer      :: ipbuf4wtp = -1 ! sulfate weight percent H2SO4 composition
  integer      :: ipbuf4reffaer = -1 ! total aerosol effective radius over all groups (cm)
  integer      :: ipbuf4reff(NGROUP) = -1 ! aerosol effective radius per group (m)
  integer      :: ipbuf4numnkg(NGROUP) = -1 ! aerosol number density (#/kg air)
  integer      :: ipbuf4sad(NBIN,NGROUP) = -1 ! aerosol surface area density per bin (cm2/cm3)
  integer      :: ipbuf4elem1mr(NBIN,NGROUP) = -1
  integer      :: ipbuf4binnkg(NBIN,NGROUP) = -1
  integer      :: ipbuf4kappa(NBIN,NGROUP) = -1 ! hygroscopicity factor per bin
  integer      :: ipbuf4wetr(NBIN,NGROUP) = -1 ! aerosol wet radius per bin
  integer      :: ipbuf4dryr(NBIN,NGROUP) = -1 ! aerosol dry radius per bin
  integer      :: ipbuf4rmass(NBIN,NGROUP) = -1 ! aerosol mass per bin
  integer      :: ipbuf4soa(NBIN) = -1
  integer      :: ipbuf4jno2 = -1
  real(kind=f) :: aeronet_fraction(NBIN)  !! fraction of BC dV/dlnr in each bin (100%)

  integer :: bc_srfemis_ndx=-1, oc_srfemis_ndx=-1

contains

  !! Defines all the CARMA components (groups, elements, solutes and gases) and process
  !! (coagulation, growth, nucleation) that will be part of the microphysical model.
  !!
  !!  @version May-2009
  !!  @author  Chuck Bardeen
  subroutine CARMA_DefineModel(carma, rc)

    use physics_buffer, only: pbuf_add_field, dtype_r8

    type(carma_type), intent(inout)    :: carma     !! the carma object
    integer, intent(out)               :: rc        !! return code, negative indicates failure


    ! Local variables
    integer                            :: LUNOPRT              ! logical unit number for output
    character(len=2)                   :: outputname,outputbin
    logical                            :: do_print             ! do print output?
    complex(kind=f)                    :: refidx(NWAVE)        ! refractice indices
    complex(kind=f)                    :: refidxS(NWAVE)       ! refractice indices for Shell
    complex(kind=f)                    :: refidxC(NWAVE)       ! refractice indices for Core

    integer                            :: igroup,ibin
    character(len=8)                   :: sname                ! short (CAM) name

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

    !call CARMAGROUP_Create(carma, I_GRP_PURSUL, "sulfate", rmin_PRSUL, vmrat_PRSUL, I_SPHERE, 1._f, .false., &
    !                       rc, irhswell=I_WTPCT_H2SO4, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
    !                       scavcoef=0.1_f, is_sulfate=.true., shortname="PRSULF", icoreshell=0, &
    !                       refidx = refidx, refidxS = refidx, refidxC = refidx, do_mie=.true.,imiertn=I_MIERTN_TOON1981)
    !if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    call CARMAGROUP_Create(carma, I_GRP_PRSUL, "sulfate", rmin_PRSUL, vmrat_PRSUL, I_SPHERE, 1._f, .false., &
                           rc, irhswell=I_WTPCT_H2SO4, do_wetdep=.false., do_drydep=.true., solfac=0.3_f, &
                           scavcoef=0.1_f, is_sulfate=.true., shortname="PRSUL", &
                           imiertn=I_MIERTN_TOON1981)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')


    !call CARMAGROUP_Create(carma, I_GRP_MIXAER, "mixed aerosol", rmin_MIXAER, vmrat_MIXAER, I_SPHERE, 1._f, .false., &
    !                       rc, do_wetdep=.true., do_drydep=.true., solfac=0.2_f, &
    !                       scavcoef=0.1_f, shortname="CRMIX", refidx=refidx, &
    !                       refidxS=refidxS, refidxC=refidxC, do_mie=.true., &
    !                       irhswell=I_MIX, irhswcomp=I_SWG_URBAN, icoreshell=1,imiertn=I_MIERTN_TOON1981)
    !if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')

    call CARMAGROUP_Create(carma, I_GRP_MXAER, "mixed aerosol", rmin_MXAER, vmrat_MXAER, I_SPHERE, 1._f, .false., &
                           rc, do_wetdep=.false., do_drydep=.true., solfac=0.2_f, &
                           scavcoef=0.1_f, shortname="MXAER", refidx=refidx, irhswell=I_PETTERS, imiertn=I_MIERTN_TOON1981,neutral_volfrc=-1._f)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddGroup failed.')


    ! Define the Elements
    !
    ! NOTE: For CAM, the optional shortname needs to be provided for the group. These names
    ! should be 6 characters or less and without spaces.
    call CARMAELEMENT_Create(carma, I_ELEM_PRSUL, I_GRP_PRSUL, "Sulfate", &
                             RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="PRSUL")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXAER,  I_GRP_MXAER, "Sulfate in mixed sulfate", &
                             RHO_SULFATE, I_VOLATILE, I_H2SO4, rc,  kappa=Kappa_SULF, shortname="MXAER")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXOC,   I_GRP_MXAER, "organic carbon", &
                             RHO_obc, I_COREMASS, I_OC, rc, kappa=Kappa_OC, shortname="MXOC")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXSOA,   I_GRP_MXAER, "secondary organic aerosol", &
                             RHO_obc, I_COREMASS, I_SOA, rc, kappa=Kappa_SOA, shortname="MXSOA")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXBC,   I_GRP_MXAER, "black carbon", &
                             RHO_obc, I_COREMASS, I_BC, rc, kappa=Kappa_BC, shortname="MXBC")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXDUST, I_GRP_MXAER, "dust", &
                             RHO_DUST, I_COREMASS, I_DUST, rc,  kappa=Kappa_DUST, shortname="MXDUST")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')

    call CARMAELEMENT_Create(carma, I_ELEM_MXSALT, I_GRP_MXAER, "SALT in mixed sulfate", &
                             RHO_SALT, I_COREMASS, I_SALT, rc, kappa=Kappa_SALT, shortname="MXSALT")
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddElement failed.')


    ! Define the Solutes



    ! Define the Gases
    call CARMAGAS_Create(carma, I_GAS_H2O, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, &
                         rc, shortname = "Q", ds_threshold=-0.2_f)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')

    call CARMAGAS_Create(carma, I_GAS_H2SO4, "Sulfuric Acid", WTMOL_H2SO4, I_VAPRTN_H2SO4_AYERS1980, &
                          I_GCOMP_H2SO4, rc, shortname = "H2SO4")
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMAGAS_Create failed.')


    ! Define the Processes

    call CARMA_AddGrowth(carma, I_ELEM_PRSUL, I_GAS_H2SO4, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddGrowth(carma, I_ELEM_MXAER, I_GAS_H2SO4, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddGrowth failed.')

    call CARMA_AddNucleation(carma, I_ELEM_PRSUL, I_ELEM_PRSUL, I_HOMNUC, 0._f, rc, igas=I_GAS_H2SO4)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddNucleation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_PRSUL, I_GRP_PRSUL, I_COLLEC_FUCHS, rc)
    if (rc < RC_OK) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_PRSUL, I_GRP_MXAER, I_GRP_MXAER, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    call CARMA_AddCoagulation(carma, I_GRP_MXAER, I_GRP_MXAER, I_GRP_MXAER, I_COLLEC_DATA, rc)
    if (rc < 0) call endrun('CARMA_DefineModel::CARMA_AddCoagulation failed.')

    !----------------- add pbuf ------------------
    do igroup = 1, NGROUP

      call CARMAGROUP_Get(carma, igroup, rc, shortname=sname)
      if (rc < 0) call endrun('carma_register::CARMAGROUP_Get failed.')
      !write(*,*) "igroup",igroup,"sname",sname

      ! effective radius. ERMIXAER
      call pbuf_add_field('ER'//trim(sname), 'global', dtype_r8, (/pcols, pver/), ipbuf4reff(igroup))
      !write(*,*) "'ER'//trim(sname)",'ER'//trim(sname)

      ! number density #/kg  NBMXAER
      call pbuf_add_field('NB'//trim(sname), 'global', dtype_r8, (/pcols, pver/), ipbuf4numnkg(igroup))
      !write(*,*) "'NB'//trim(sname)",'NB'//trim(sname)

      ! sulfate mass and number density for each bin
      ! e.g. CRSULF01 first element mass mixing ratio; NBMXAER01 #/kg
      do ibin=1,NBIN
         write (outputbin, "(I2.2)") ibin
         write (outputname,"(A2)") trim(sname)
         if (trim(sname)//outputbin .ne. outputname//"SULF"//outputbin) then
            !write(*,*) "sname//outputbin",trim(sname)//outputbin
            call pbuf_add_field(outputname//"SULF"//outputbin,'global', dtype_r8, (/pcols, pver/), ipbuf4elem1mr(ibin,igroup))
            !write(*,*) outputname//"SULF"//outputbin,"   ", trim(sname)//outputbin
         end if
         call pbuf_add_field("NB"//trim(sname)//outputbin,'global', dtype_r8, (/pcols, pver/), ipbuf4binnkg(ibin,igroup))
         call pbuf_add_field(trim(sname)//outputbin//"_kappa",'global', dtype_r8, (/pcols, pver/), ipbuf4kappa(ibin,igroup))
         !write(*,*) trim(sname)//outputbin//"_kappa"
         call pbuf_add_field(trim(sname)//outputbin//"_wetr",'global', dtype_r8, (/pcols, pver/), ipbuf4wetr(ibin,igroup))
         call pbuf_add_field(trim(sname)//outputbin//"_dryr",'global', dtype_r8, (/pcols, pver/), ipbuf4dryr(ibin,igroup))
         call pbuf_add_field(trim(sname)//outputbin//"_rmass",'global', dtype_r8, (/pcols/), ipbuf4rmass(ibin,igroup))
         call pbuf_add_field(trim(sname)//outputbin//"_sad", 'global', dtype_r8, (/pcols, pver/), ipbuf4sad(ibin,igroup))
         if (igroup==I_GRP_MXAER) then
           call pbuf_add_field("DQDT_MXSOA"//outputbin,'global',dtype_r8,(/pcols,pver/), ipbuf4soa(ibin))
         end if
      end do
   end do

    ! total surface area density (cm2/cm3)
    call pbuf_add_field('SADSULF', 'global', dtype_r8, (/pcols, pver/), ipbuf4sadsulf)
    ! total effective radius (cm) for REFF_AERO history output
    call pbuf_add_field('REFFAER', 'global', dtype_r8, (/pcols, pver/), ipbuf4reffaer)
    ! weight percent H2SO4
    call pbuf_add_field('WTP','global', dtype_r8, (/pcols, pver/), ipbuf4wtp)
    ! no2 photolysis rate constant (/sec)
    call pbuf_add_field('JNO2', 'global', dtype_r8, (/pcols,pver/), ipbuf4jno2)

    !---------------------------------------------

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
    use camsrfexch, only: cam_in_t

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

    ! local variables
    character(len=8)                     :: snamecore                !! short (CAM) name
    character(len=8)                     :: shortname                !! short (CAM) name
    real(r8), pointer, dimension(:,:)    :: dqdt_soa              !! soa tendency due to gas-aerosol exchange  kg/kg/s
    real(r8), pointer, dimension(:,:)    :: jno2_rate             !! jno2 tendency due to gas-aerosol exchange  kg/kg/s
    real(r8)                             :: mmr_total(cstate%f_NZ)!! mass mixing ratio of a group (kg/kg)
    real(r8)                             :: mmr_core(cstate%f_NZ)!! mass mixing ratio of the core (kg/kg)
    real(r8)                             :: mmr_soa(cstate%f_NZ)  !! mass mixing ratio of soa element (kg/kg)
    real(r8)                             :: mmr(cstate%f_NZ)      !! mass mixing ratio per bin (kg/kg)
    real(r8)                             :: delta_soa(cstate%f_NZ)     !! mass mixing ratio differences from soa gas-aerosol-exchange
    integer                              :: icorelem(NELEM), ncore,ienconc,icore, ielem, ielem_soa, igroup, ibin, icomposition, n, err

    ! Default return code.
    rc = RC_OK

    ! get no2 photolysis rates if they exist
    call pbuf_get_field(pbuf, ipbuf4jno2, jno2_rate)     ! surface area density

    ! get SOA tendency pbuf field for the mixed group and every bin

    igroup = I_GRP_MXAER

    call CARMAGROUP_Get(carma, igroup, rc, ienconc=ienconc,ncore=ncore,icorelem=icorelem)

    do ibin = 1, NBIN

       call CARMASTATE_GetBin(cstate, ienconc, ibin, mmr_total(:), rc)
       if (rc /= RC_OK) call endrun('CARMA_DiagnoseBins::CARMASTATE_GetBin failed.')

       mmr_soa(:) = 0.0_r8
       mmr_core(:) = 0.0_r8

       do ielem = 1, ncore
          delta_soa(:) = 0.0_r8

          call CARMASTATE_GetBin(cstate, icorelem(ielem), ibin, mmr(:), rc)

          call CARMAELEMENT_GET(carma, icorelem(ielem), rc, igroup=igroup, shortname=snamecore, icomposition=icomposition)

          if (icomposition==I_SOA) then
             call pbuf_get_field(pbuf, ipbuf4soa(ibin), dqdt_soa)     ! surface area density
             ielem_soa = ielem
             mmr_soa = mmr

             !add soa tendency to mmr_soa  ; dqdt = kg/kg/s

             mmr_soa(:) = mmr_soa(:) + dqdt_soa(icol,:) * dt

             ! substract photolysis rates
             mmr_soa(:) = mmr_soa(:) - 0.0004_r8*jno2_rate(icol,:)*mmr_soa(:) * dt

             mmr_soa(:) = max(mmr_soa(:),0.0_r8)

             delta_soa(:) = mmr_soa(:) - mmr(:)

             ! set mmr to new mmr

             mmr(:) = mmr_soa(:)

          end if  !mxsoa
          mmr_core(:) = mmr_core(:) + mmr(:)

       end do  !ielem

       !update mmr_total and check that not smaller than core mass
       do n = 1, cstate%f_NZ
          mmr_total(n) = mmr_total(n) + delta_soa(n)
          if (mmr_total(n) .lt.  mmr_core(n))  then
             mmr_total(n) = mmr_core(n)
          end if
       end do

       call CARMASTATE_SetBin(cstate, icorelem(ielem_soa), ibin, mmr_soa, rc)
       if (rc /= RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_SetBin failed.')

       call CARMASTATE_SetBin(cstate, ienconc, ibin, mmr_total, rc)
       if (rc /= RC_OK) call endrun('CARMA_DiagnoseBins::CARMAGROUP_SetBin failed.')

    end do  !ibin

  end subroutine CARMA_DiagnoseBins


  !! For diagnostic groups, determines the tendencies on the CAM state from the CARMA bins.
  !!
  !!  @version July-2009
  !!  @author  Chuck Bardeen
  subroutine CARMA_DiagnoseBulk(carma, cstate, cam_out, state, pbuf, ptend, icol, dt, rc, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, tnd_qsnow, tnd_nsnow, re_ice)
    use camsrfexch, only: cam_out_t

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

   ! Local variables
    real(r8)                             :: numberDensity(cstate%f_NZ)
    real(r8)                             :: totad(cstate%f_NZ)
    real(r8)                             :: ad(cstate%f_NZ)       !! aerosol wet surface area density (cm2/cm3)
    real(r8)                             :: totreff(cstate%f_NZ)  !! total volume density, used to calculate total effective radius (cm) for history output
    real(r8)                             :: reff(cstate%f_NZ)     !! wet effective radius (m)
    real(r8)                             :: mmr(cstate%f_NZ)      !! mass mixing ratio per bin (kg/kg)
    real(r8)                             :: mmr_total(cstate%f_NZ)!! mass mixing ratio of a group (kg/kg)
    real(r8)                             :: coremmr(cstate%f_NZ)  !! mmr of all the core
    real(r8)                             :: mmr_gas(cstate%f_NZ)  !! gas mass mixing ratio (kg/kg)
    real(r8)                             :: numnkg(cstate%f_NZ)   !! total number density (#/kg)
    real(r8)                             :: r_wet(cstate%f_NZ)    !! Sulfate aerosol bin wet radius (cm)
    real(r8)                             :: elem1mr(cstate%f_NZ)  !! First element mass mixing ratio (kg/kg)
    real(r8)                             :: binnkg(cstate%f_NZ)   !! number density per bin (#/kg)
    real(r8)                             :: kappa(cstate%f_NZ)    !! hygroscopicity parameter (Petters & Kreidenweis, ACP, 2007)
    real(r8)                             :: rhoa_wet(cstate%f_NZ) !! wet air density (kg/m3)
    real(r8)                             :: wtpct(cstate%f_NZ)    !! sulfate weight percent
    real(r8)                             :: rmass(NBIN)           !! dry mass
    real(r8)                             :: rhop_dry(cstate%f_NZ) !! dry particle density [g/cm3]

    integer                              :: ibin, igroup, igas, icomposition
    integer                              :: icorelem(NELEM), ncore,ienconc,icore
    character(len=8)                     :: sname                 !! short (CAM) name

    real(r8), pointer, dimension(:,:)    :: sadsulf_ptr           !! Total surface area density pointer (cm2/cm3)
    real(r8), pointer, dimension(:,:)    :: reffaer_ptr           !! Total effective radius pointer (cm) for history output
    real(r8), pointer, dimension(:,:)    :: wtp_ptr		  !! weight percent pointer
    real(r8), pointer, dimension(:,:)    :: sad_ptr               !! Surface area density pointer
    real(r8), pointer, dimension(:,:)    :: reff_ptr              !! Effective radius pointer
    real(r8), pointer, dimension(:,:)    :: numnkg_ptr            !! Each group number density pointer
    real(r8), pointer, dimension(:,:)    :: binnkg_ptr            !! Each bin number density pointer
    real(r8), pointer, dimension(:,:)    :: elem1mr_ptr           !! First element mmr pointer
    real(r8), pointer, dimension(:,:)    :: kappa_ptr		  !! kappa pointer
    real(r8), pointer, dimension(:,:)    :: wetr_ptr           !! wet radius pointer
    real(r8), pointer, dimension(:,:)    :: dryr_ptr             !! dry radius

    ! Default return code.
    rc = RC_OK

    ! Get the air density
    call CARMASTATE_GetState(cstate, rc, rhoa_wet=rhoa_wet)

    totad(:) = 0.0_r8   ! total aerosol surface area density (cm2/cm3)
    totreff(:) = 0.0_r8 ! total volume density, used to calculate total effective radius (cm) for REFF_AERO history output

    ! calculated SAD, RE, and number density (#/kg) for each group
    do igroup = 1, NGROUP

      ad(:)  = 0.0_r8     ! wet aerosol surface area density (cm2/cm3)
      reff(:)  = 0.0_r8    ! effective radius (m)
      numnkg(:) = 0.0_r8    ! number density (#/kg)

      call CARMAGROUP_Get(carma, igroup, rc, ienconc=ienconc,ncore=ncore,icorelem=icorelem, rmass=rmass)
      do ibin = 1, NBIN

        call CARMASTATE_GetBin(cstate, ienconc, ibin, mmr_total(:), rc, &
                             numberDensity=numberDensity, r_wet=r_wet)
        if (rc < 0) call endrun('CARMA_DiagnoseBulk::CARMASTATE_GetBin failed.')

	if (ncore .ne. 0)then
          do icore = 1,ncore
            call CARMASTATE_GetBin(cstate, icorelem(icore), ibin, mmr(:), rc)
	    call CARMAELEMENT_Get(carma, icorelem(icore), rc, icomposition=icomposition)
          end do
        end if

        ! Calculate the total densities.
        ! NOTE: Calculate AD in cm2/cm3.
        if (numberDensity(1) /= CAM_FILL) then
          ad(:)  = ad(:)  + numberDensity(:) * (r_wet(:)**2)
          reff(:) = reff(:) + numberDensity(:) * (r_wet(:)**3)
        end if
      end do

      totreff(:) = totreff(:) + reff(:) ! volume density in cm3
      reff(:) = reff(:) / ad(:) ! wet effective radius in cm
      reff(:) = reff(:) / 100.0_r8 ! cm -> m

      ad(:) = ad(:) * 4.0_r8 * PI ! surface area density in cm2/cm3

      ! airdensity from carma state
      ! convert the number density from #/cm3 to #/kg
      ! number Density #/cm3; rhoa_wet kg/m3
      numnkg(:) = 1.e6_r8*numberDensity(:)/rhoa_wet(:)    !#/kg

      call pbuf_get_field(pbuf, ipbuf4reff(igroup), reff_ptr)   ! effective radius m
      call pbuf_get_field(pbuf, ipbuf4numnkg(igroup), numnkg_ptr) ! number density #/kg

      ! put variables in physics buffer
      reff_ptr(icol, :cstate%f_NZ) = reff(:cstate%f_NZ)
      numnkg_ptr(icol, :cstate%f_NZ)= numnkg(:cstate%f_NZ)

      !calculate the surface area density including all bins (cm2/cm3)
      totad(:) = totad(:)+ad(:)

    end do

    totreff(:) = 4.0_r8 * PI * totreff(:) / totad(:) ! cm

    call pbuf_get_field(pbuf, ipbuf4sadsulf, sadsulf_ptr)     ! surface area density (cm2/cm3)
    sadsulf_ptr(icol, :cstate%f_NZ)  = totad(:cstate%f_NZ)

    call pbuf_get_field(pbuf, ipbuf4reffaer, reffaer_ptr)     ! effective radius (cm)
    reffaer_ptr(icol, :cstate%f_NZ)  = totreff(:cstate%f_NZ)

    do igas = 1,NGAS
      if(igas .eq. I_GAS_H2SO4)then ! only output the sulfate weight percent
        call CARMASTATE_GetGas(cstate, igas, mmr_gas(:), rc, wtpct=wtpct)
        call pbuf_get_field(pbuf, ipbuf4wtp, wtp_ptr)
        wtp_ptr(icol, :cstate%f_NZ)  = wtpct(:cstate%f_NZ)
      end if
    end do

    ! add new CRSULF element that is only the sulfur part from MIXAER:
    ! calculate CRSULF mass mixing ratio : MXAER minus CROC+CRBC+CRBC+CRDUST+CRSALT (+SUM(CRSOA*)

    do igroup = 1, NGROUP
      call CARMAGROUP_Get(carma, igroup, rc, shortname=sname,ienconc=ienconc,ncore=ncore,icorelem=icorelem, rmass=rmass)
      !write(*,*) "igroup",igroup,"ncore",ncore,"ienconc",ienconc,"icorelem",icorelem
      do ibin = 1, NBIN

        elem1mr(:) = 0._r8
        binnkg(:) = 0._r8
	coremmr(:) = 0._r8

        call CARMASTATE_GetBin(cstate, ienconc, ibin, mmr_total(:), rc, numberDensity=numberDensity, kappa=kappa, r_wet=r_wet,rhop_dry=rhop_dry)

        if (ncore .ne. 0)then
          do icore = 1,ncore
            call CARMASTATE_GetBin(cstate, icorelem(icore), ibin, mmr(:), rc)
            coremmr(:) = coremmr(:) + mmr(:)
          end do

          if (any(coremmr(:) .gt. mmr_total(:))) then
            where(coremmr(:) .gt. mmr_total(:)) mmr_total = coremmr
	    call CARMASTATE_SetBin(cstate, ienconc, ibin, mmr_total(:), rc)
            call CARMASTATE_GetBin(cstate, ienconc, ibin, mmr_total(:), rc, numberDensity=numberDensity, r_wet=r_wet,rhop_dry=rhop_dry)
          end if
	  elem1mr(:) = mmr_total(:)-coremmr(:)
          elem1mr(:) = max(elem1mr(:),0._f)

          call pbuf_get_field(pbuf, ipbuf4elem1mr(ibin,igroup), elem1mr_ptr)
	  elem1mr_ptr(icol, :cstate%f_NZ) = elem1mr(:cstate%f_NZ)
	else
	  elem1mr(:) = mmr_total(:)
	  call pbuf_get_field(pbuf, ipbuf4elem1mr(ibin,igroup), elem1mr_ptr)
          elem1mr_ptr(icol, :cstate%f_NZ) = elem1mr(:cstate%f_NZ)
	end if
        binnkg(:) = 1.e6_r8*numberDensity(:)/rhoa_wet(:)     !#/kg

        call pbuf_get_field(pbuf, ipbuf4binnkg(ibin,igroup), binnkg_ptr)
	call pbuf_get_field(pbuf, ipbuf4kappa(ibin,igroup), kappa_ptr)

        binnkg_ptr(icol, :cstate%f_NZ) = binnkg(:cstate%f_NZ)
	kappa_ptr(icol, :cstate%f_NZ) = kappa(:cstate%f_NZ)

	call pbuf_get_field(pbuf, ipbuf4wetr(ibin,igroup), wetr_ptr)
	call pbuf_get_field(pbuf, ipbuf4dryr(ibin,igroup), dryr_ptr)
	call pbuf_get_field(pbuf, ipbuf4sad(ibin,igroup), sad_ptr)

	wetr_ptr(icol, :cstate%f_NZ) = r_wet(:cstate%f_NZ)
        !dryr_ptr(icol, :cstate%f_NZ) = (rmass(ibin)/(4._f/3._f*PI*rhop_dry(:cstate%f_NZ)))**(1._f/3._f) !cm
        ! r = (mass / rho / 4 * 3 / PI)^(1/3)
        dryr_ptr(icol, :cstate%f_NZ) = (rmass(ibin)/rhop_dry(:cstate%f_NZ) / 4._f * 3._f / PI)**(1._f/3._f) !cm

        sad_ptr(icol, :cstate%f_NZ) = 4.0_r8 * PI * numberDensity(:cstate%f_NZ) * (r_wet(:cstate%f_NZ)**2) !cm2/cm3

        !write(*,*) 'CARMA igroup, ibin, rmass(ibin), rhop_dry',igroup, ibin, rmass(ibin), rhop_dry(:cstate%f_NZ)
        !write(*,*) 'CARMA dryr igroup, ibin',igroup, ibin, dryr_ptr(icol, :cstate%f_NZ)

      end do !NBIN
    end do !NGROUP

    return
  end subroutine CARMA_DiagnoseBulk


  !! Calculates the emissions for CARMA aerosol particles. By default, there is no
  !! emission, but this routine can be overridden for models that wish to have
  !! an aerosol emission.
  !!
  !! @author  Lin Su, Pengfei Yu, Chuck Bardeen
  !! @version Dec-2010
  subroutine CARMA_EmitParticle(carma, ielem, ibin, icnst, dt, state, cam_in, tendency, surfaceFlux, pbuf, rc)
    use ppgrid,        only: pcols, pver
    use physics_types, only: physics_state
    use phys_grid,     only: get_lon_all_p, get_lat_all_p
    use time_manager,  only: get_curr_date, get_perp_date, get_curr_calday, &
                             is_perpetual
    use camsrfexch,    only: cam_in_t
    use cam_history,   only: outfld

    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ielem                 !! element index
    integer, intent(in)                :: ibin                  !! bin index
    integer, intent(in)                :: icnst                 !! consituent index
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_state), intent(in)    :: state                 !! physics state
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: tendency(pcols, pver) !! constituent tendency (kg/kg/s)
    real(r8), intent(out)              :: surfaceFlux(pcols)    !! constituent surface flux (kg/m^2/s)
    type(physics_buffer_desc), pointer :: pbuf(:)               !! physics buffer
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    integer      :: ilat(pcols)             ! latitude index
    integer      :: ilon(pcols)             ! longitude index
    real(r8)     :: clat(pcols)             ! latitude
    integer      :: lchnk                   ! chunk identifier
    integer      :: ncol                    ! number of columns in chunk
    integer      :: icol                    ! column index
    real(r8)     :: calday                  ! current calendar day
    integer      :: yr                      ! year
    integer      :: mon                     ! month
    integer      :: day                     ! day of month
    integer      :: ncsec                   ! time of day (seconds)
    real(r8)     :: smoke(pcols)            ! smoke emission flux (molecues/cm2/s)
    integer      :: igroup                  ! the index of the carma aerosol group
    character(len=32) :: shortname          ! the shortname of the group



    ! -------- local variables added for dust and sea-salt model ------------
    real(r8), parameter :: ch = 0.5e-9_r8                     ! dimensional factor & tuning number,
    real(r8)            :: rmass(NBIN)                        ! bin mass (g)
    real(r8)            :: r                                  ! bin center (cm)
    real(r8)            :: rdust                              ! dust bin center (cm)
    real(r8)            :: dustFlux                           ! dust flux (kg/m2/s)
    real(r8)            :: rsalt                              ! salt bin center (cm)
    real(r8)            :: drsalt                             ! salt bin width (cm)
    real(r8)            :: rhop(NBIN)                         ! element density (g/cm3)
    real(r8)            :: vrfact
    real(r8)            :: uth                                ! threshold wind velocity (m/s)
    real(r8)            :: uv10                               ! 10 m wind speed (m/s)
    real(r8)            :: cd10                               ! 10-m drag coefficient ()
    real(r8)            :: wwd                                ! raw wind speed (m/s)
    real(r8)            :: sp                                 ! mass fraction for soil factor
    integer             :: idustbin                           ! ibin to use for dust production, smallest silt bin for clay

! ------------ local variables added for organics model ----------------------
    real(r8)     :: dr
    real(r8)     :: aeronet(NBIN)                       ! AERONET DATA, Sep.20, 2002, Jaru Reserve, Brazil (refer to MATICHUK et al., 2008)
    real(r8)     :: saltFlux(pcols)                     ! sea salt flux to calculate marine POA
    integer      :: LUNOPRT                             ! logical unit number for output
    logical      :: do_print                            ! do print output?

    real(r8),parameter :: OMtoOCratio = 1.8_r8           ! Need better names and doc
    real(r8),parameter :: SmoketoSufaceFlux = 1.9934e-22_r8 ! SmoketoSufaceFlux = BC molecular weight
                                                            ! (12 g/mol)/avocadro constant (6e-23 #/mol) *10
    real(r8), pointer :: BCemis_ptr(:), OCemis_ptr(:)

!   currently not used
!   real(r8)     :: MPOAFlux(pcols)             ! marine POA flux
!   real(r8)     :: Fsub(pcols)                 ! marine Chlorophy11-dependent mass contribution of sub-micron organics
!   real(r8)     :: sub_micron(pcols)                   ! total sub-micron sea spray particles
!   real(r8)     :: PBAPFlux(pcols)                             ! Primary biological aerosol particles
!   real(r8)     :: spor_mon_N(12) = (/0.5_r8,0.5_r8,0.5_r8,1.5_r8,1.5_r8,1.5_r8,1.5_r8,1.5_r8,1.5_r8,0.5_r8,0.5_r8,0.5_r8/)
!   real(r8)     :: spor_mon_S(12) = (/1.5_r8,1.5_r8,1.5_r8,0.5_r8,0.5_r8,0.5_r8,0.5_r8,0.5_r8,0.5_r8,1.5_r8,1.5_r8,1.5_r8/)
!   real(r8)     :: Spbin                                               ! fraction of emission for each bin
!   real(r8)     :: Rrh(pcols)                                  ! RH effect for spore emission
!   real(r8)     :: sporemass                   ! spore mass per particle

    ! Default return code.
    rc = RC_OK
    smoke(:) = -huge(1._r8)

    ! Determine the day of year.
    calday = get_curr_calday()
    if ( is_perpetual() ) then
      call get_perp_date(yr, mon, day, ncsec)
    else
      call get_curr_date(yr, mon, day, ncsec)
    end if

    ! Determine the latitude and longitude of each column.
    lchnk = state%lchnk
    ncol = state%ncol

    ! Add any surface flux here.
    surfaceFlux(:ncol) = 0.0_r8

    ! For emissions into the atmosphere, put the emission here.
    !
    ! NOTE: Do not set tendency to be the surface flux. Surface source is put in to
    ! the bottom layer by vertical diffusion. See vertical_solver module, line 355.
    tendency(:ncol, :pver) = 0.0_r8

     ! Add Emission (surfaceFlux) here.

    !!*******************************************************************************************************

    !! add an element, first element is total number with emission from both OC and BC;
    !! second element is BC mass
    !! by Pengfei Yu
    !! Feb.22 2012
    !!*******************************************************************************************************


    call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname)
    if (RC < RC_ERROR) return

    call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname, rmass=rmass)
    if (RC < RC_ERROR) return

     !!*******************************************************************************************************

    !if (masterproc) then
    !  call CARMA_Get(carma, rc, do_print=do_print, LUNOPRT=LUNOPRT)
    !
    ! if (do_print) then
    !   write(carma%f_LUNOPRT,*) 'AERONET', aeronet
    !   write(carma%f_LUNOPRT,*) 'dr', dr
    !   write(carma%f_LUNOPRT,*) 'r', r
    ! end if
    !end if

    !!*******************************************************************************************************

    if(carma_BCOCemissions == 'Specified')then
      call pbuf_get_field(pbuf, bc_srfemis_ndx, BCemis_ptr)
      call pbuf_get_field(pbuf, oc_srfemis_ndx, OCemis_ptr)
    end if

    ! Organic carbon emssions
    if ((ielem == I_ELEM_MXOC) .or. (ielem == I_ELEM_MXAER)) then
       if (carma_BCOCemissions == 'Yu2015') then
          call get_lat_all_p(lchnk, ncol, ilat)
          call get_lon_all_p(lchnk, ncol, ilon)
          do icol = 1,ncol
             smoke(icol) = OCnew(ilat(icol), ilon(icol), mon)*OMtoOCratio
          end do
       elseif(carma_BCOCemissions == 'Specified')then
          smoke(:ncol) = OCemis_ptr(:ncol)
       end if

!  st  scip Fsub PBAFlux etcfor now
       surfaceFlux(:ncol) = surfaceFlux(:ncol) + smoke(:ncol)*aeronet_fraction(ibin)*SmoketoSufaceFlux
    end if

    ! Black carbon emissions
    if ((ielem == I_ELEM_MXBC) .or. (ielem == I_ELEM_MXAER)) then
       if (carma_BCOCemissions == 'Yu2015') then
          do icol = 1,ncol
             smoke(icol) = BCnew(ilat(icol), ilon(icol), mon)
          end do
       elseif(carma_BCOCemissions == 'Specified') then
          smoke(:ncol) = BCemis_ptr(:ncol)
       end if

       surfaceFlux(:ncol) = surfaceFlux(:ncol) + smoke(:ncol)*aeronet_fraction(ibin)*SmoketoSufaceFlux
    end if

    ! Dust emissions
    if ((ielem == I_ELEM_MXDUST) .or. (ielem == I_ELEM_MXAER)) then

      ! The radius should be determined by the dust density not the group
      ! density
      call CARMAELEMENT_Get(carma, I_ELEM_MXDUST, rc, rho=rhop)
      if (RC < RC_ERROR) return

      ! Calculate the radius assuming that all the mass will be emitted as this
      ! element.
      rdust = (3._r8 * rmass(ibin) / 4._r8 / PI / rhop(ibin)) ** (1._r8 / 3._r8)

      ! Is this clay or silt?
      !
      ! NOTE: It is assumed that 90% of the mass will be silt and 10% will
      ! be clay.
      !
      ! NOTE: For clay bins, use the smallest silt bin to calculate the
      ! mass and then scale that into each clay bin based upon interpolation of
      ! Tegen and Lacis [1996].
      if (rdust >= rClay) then
        sp         = 0.9_r8 / nSilt
        idustbin   = ibin
      else
        sp         = 0.1_r8 / nClay
        idustbin   = nClay + 1
      end if

      ! Process each column.
      do icol = 1,ncol

        call CARMA_SurfaceWind(carma, icol, ielem, igroup, idustbin, cam_in, uv10, wwd, uth, rc)

        ! Is the wind above the threshold for dust production?
        if (sqrt(wwd) > uth) then
          dustFlux = ch * soil_factor(icol, lchnk) * sp * &
                              wwd * (sqrt(wwd) - uth)
        else
          dustFlux = 0._r8
        endif

        ! Scale the clay bins based upon the smallest silt bin.
        dustFlux = clay_mf(ibin) * dustFlux

        ! Add the dust flux to the accumulated emissions (important for I_ELEM_MXAER)
        surfaceFlux(icol) = surfaceFlux(icol) + dustFlux
      end do

      ! For debug purposes, output the soil erosion factor.
      call outfld('CRSLERFC', soil_factor(:ncol, lchnk), ncol, lchnk)
    end if


    ! Sea salt emissions
    if ((ielem == I_ELEM_MXSALT) .or. (ielem == I_ELEM_MXAER)) then

      ! The radius should be determined by the dust density not the group
      ! density
      call CARMAELEMENT_Get(carma, I_ELEM_MXSALT, rc, rho=rhop)
      if (RC < RC_ERROR) return

      ! Calculate the radius assuming that all the mass will be emitted as sea
      ! salt.
      vrfact = ((3._r8/2._r8 / PI / (vmrat_MXAER + 1._r8))**(1._r8 / 3._r8)) * ((vmrat_MXAER**(1._r8 / 3._r8)) - 1._r8)
      rsalt = (3._r8 * rmass(ibin) / 4._r8 / PI / rhop(ibin))**(1._r8 / 3._r8)
      drsalt = vrfact * ((rmass(ibin)/rhop(ibin))**(1._r8 / 3._r8))

      ! get sea spray aerosol flux first (for ibin; SaltFlux(:ncol) unit:kg/m2/s)
      call CARMA_SaltFlux(carma, ibin, state, rsalt, drsalt, rmass(ibin), cam_in, saltFlux, rc)

!st  not used currently  but done by Pengfei
       !! introduce marine POA emission, use ChlorophyII-dependent mass contribution of OC
       !! see Gantt et al., 2009
       !! for sub-micron, I use sea salt flux instead of sub-micron marine particles
       !! needed to verify later
       !! Added by Pengfei Yu
       !! Oct.6.2012
       ! get [Chl-a] data
  !!   do icol = 1, ncol
  !!       if (Chla(ilat(icol), ilon(icol)) .lt. 0._r8) then
  !!          Fsub(icol) = 0._r8
  !!       else
  !!          Fsub(icol) = Chla(ilat(icol), ilon(icol)) * 0.63_r8 + 0.1_r8
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
      surfaceFlux(:ncol) = surfaceFlux(:ncol) + saltFlux(:ncol)
    end if

    return
  end subroutine CARMA_EmitParticle


  !! Allows the model to perform its own initialization in addition to what is done
  !! by default in CARMA_init.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  subroutine CARMA_InitializeModel(carma, lq_carma, pbuf2d, rc)
    use cam_history,  only: addfld,  horiz_only
    use constituents, only: pcnst

    type(carma_type), intent(in)       :: carma                 !! the carma object
    logical, intent(inout)             :: lq_carma(pcnst)       !! flags to indicate whether the constituent
                                                                !! could have a CARMA tendency
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

    ! -------- local variables ----------
    integer            :: ibin                                ! CARMA bin index
    real(r8)           :: rdust(NBIN),robc(NBIN),drobc(NBIN),rm(NBIN),rhop(NBIN)       ! bin center (cm)
    integer            :: count_Silt                          ! count number for Silt
    integer            :: igroup                              ! the index of the carma aerosol group
    integer            :: ielem                               ! the index of the carma aerosol element
    character(len=32)  :: shortname                           ! the shortname of the element
    integer            :: LUNOPRT                             ! logical unit number for output
    logical            :: do_print                            ! do print output?

    integer :: idata,isizebin,ibin_local
    integer,parameter :: aeronet_dim1 = 22
    integer,parameter :: aeronet_dim2 = 4
    real(r8),dimension(aeronet_dim1,aeronet_dim2) :: sizedist_aeronet
    real(r8),dimension(aeronet_dim1) :: sizedist_avg
    real(r8),dimension(NBIN) :: sizedist_carmabin
    real(r8) :: rmass(NBIN) !! dry mass
    real(r8) :: vrfact

    real(r8),parameter :: size_aeronet(aeronet_dim1) = (/0.050000_r8,0.065604_r8,0.086077_r8,0.112939_r8,0.148184_r8, &
         0.194429_r8,0.255105_r8,0.334716_r8,0.439173_r8,0.576227_r8,0.756052_r8,0.991996_r8,1.301571_r8,1.707757_r8, &
         2.240702_r8,2.939966_r8,3.857452_r8,5.061260_r8,6.640745_r8,8.713145_r8,11.432287_r8,15.000000_r8/)*1.e-4_r8 !um to cm

    ! Default return code.
    rc = RC_OK

    ! Determine how many clay and how many silt bins there are, based
    ! upon the bin definitions and rClay.
    !
    ! TBD: This should use the radii rather than being hard coded.
    ! nClay = 8
    ! nSilt = NBIN - nClay
    do ielem = 1, NELEM

       ! To get particle radius, need to derive from rmass and density of dust.
       call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname, rho=rhop)
       if (RC < RC_ERROR) return

       call CARMAGROUP_GET(carma, igroup, rc, rmass=rmass)
       if (RC < RC_ERROR) return

       if (shortname .eq. "MXDUST") then

          count_Silt = 0
          do ibin = 1, NBIN

             ! Calculate the radius assuming that all the mass will be emitted as this
             ! element.
             rdust(ibin) = (3._r8 * rmass(ibin) / 4._r8 / PI / rhop(ibin))**(1._r8 / 3._r8)

             if (rdust(ibin) >= rclay) then
                count_Silt = count_Silt + 1
             else
             end if
          end do
          nSilt = count_Silt
          nClay = NBIN - nSilt
       end if
    end do

    ! Read in the soil factors.
    call CARMA_ReadSoilErosionFactor(rc)
    if (RC < RC_ERROR) return

    ! To determine Clay Mass Fraction
    do ielem = 1, NELEM
       ! To get particle radius
       call CARMAELEMENT_GET(carma, ielem, rc, igroup=igroup, shortname=shortname)
       if (RC < RC_ERROR) return

       if (shortname .eq. "MXDUST") then
          call CARMA_ClayMassFraction(carma, igroup, rdust, rc)
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

    if (carma_BCOCemissions == 'Yu2015')then
       ! Added by Pengfei Yu to read smoke emission data
       call CARMA_BCOCread(rc)
    end if
    if(carma_BCOCemissions == 'Specified')then
       bc_srfemis_ndx = pbuf_get_index("BC_srfemis")
       oc_srfemis_ndx = pbuf_get_index("OC_srfemis")
    end if

    if (is_first_step()) then
       ! Initialize physics buffer fields
       do igroup = 1, NGROUP
          call pbuf_set_field(pbuf2d, ipbuf4reff(igroup), 0.0_r8 )
          call pbuf_set_field(pbuf2d, ipbuf4numnkg(igroup), 0.0_r8 )

          call CARMAGROUP_Get(carma, igroup, rc, rmass=rmass)
          if (RC /= RC_OK) return

          do ibin=1,NBIN
             if (ipbuf4elem1mr(ibin,igroup)>0) then
                call pbuf_set_field(pbuf2d, ipbuf4elem1mr(ibin,igroup), 0.0_r8 )
             endif
             call pbuf_set_field(pbuf2d, ipbuf4binnkg(ibin,igroup), 0.0_r8 )
             call pbuf_set_field(pbuf2d, ipbuf4kappa(ibin,igroup), 0.0_r8 )
             call pbuf_set_field(pbuf2d, ipbuf4wetr(ibin,igroup), 0.0_r8 )
             call pbuf_set_field(pbuf2d, ipbuf4dryr(ibin,igroup), 0.0_r8 )
             call pbuf_set_field(pbuf2d, ipbuf4sad(ibin,igroup), 0.0_r8 )
             call pbuf_set_field(pbuf2d, ipbuf4rmass(ibin,igroup), 1.0e-3_r8*rmass(ibin) ) ! convert rmass from g to kg
             if (igroup==I_GRP_MXAER) then
                call pbuf_set_field(pbuf2d, ipbuf4soa(ibin), 0.0_r8 )
             end if
          end do
       end do

       call pbuf_set_field(pbuf2d, ipbuf4sadsulf, 0.0_r8 )
       call pbuf_set_field(pbuf2d, ipbuf4reffaer, 0.0_r8 )
       call pbuf_set_field(pbuf2d, ipbuf4wtp, 0.0_r8 )
       call pbuf_set_field(pbuf2d, ipbuf4jno2, 0.0_r8 )
    endif

    sizedist_aeronet(:aeronet_dim1,1) = (/0.000585_r8,0.006080_r8,0.025113_r8,0.052255_r8,0.079131_r8,0.081938_r8, &
         0.035791_r8,0.010982_r8,0.005904_r8,0.007106_r8,0.011088_r8,0.012340_r8,0.010812_r8,0.010423_r8, &
         0.011892_r8,0.016529_r8,0.023967_r8,0.026854_r8,0.017901_r8,0.007226_r8,0.002161_r8,0.000544_r8/)
    sizedist_aeronet(:aeronet_dim1,2) = (/0.000541_r8,0.006524_r8,0.026103_r8,0.050825_r8,0.077730_r8,0.080545_r8, &
         0.035400_r8,0.011143_r8,0.005753_r8,0.006095_r8,0.008730_r8,0.010794_r8,0.011517_r8,0.012051_r8, &
         0.012362_r8,0.014710_r8,0.019738_r8,0.022156_r8,0.014892_r8,0.005976_r8,0.001891_r8,0.000573_r8/)
    sizedist_aeronet(:aeronet_dim1,3) = (/0.000747_r8,0.009291_r8,0.043556_r8,0.099216_r8,0.142377_r8,0.108606_r8, &
         0.043723_r8,0.016385_r8,0.008318_r8,0.005597_r8,0.004431_r8,0.004131_r8,0.004980_r8,0.007484_r8, &
         0.011795_r8,0.017235_r8,0.022404_r8,0.025216_r8,0.022521_r8,0.013752_r8,0.005051_r8,0.001057_r8/)
    sizedist_aeronet(:aeronet_dim1,4) = (/0.000979_r8,0.007724_r8,0.034451_r8,0.090410_r8,0.135893_r8,0.103115_r8, &
         0.046047_r8,0.018989_r8,0.009149_r8,0.005034_r8,0.003199_r8,0.002680_r8,0.003249_r8,0.005105_r8, &
         0.008370_r8,0.012542_r8,0.016973_r8,0.021107_r8,0.022077_r8,0.015639_r8,0.006001_r8,0.001115_r8/)

    sizedist_avg(:) = 0._r8
    do idata = 1,aeronet_dim2
       sizedist_avg(:) = sizedist_avg(:) + sizedist_aeronet(:,idata)
    end do
    sizedist_avg(:) = sizedist_avg(:)*0.25_r8

    do igroup = 1,NGROUP
      call CARMAGROUP_GET(carma, igroup, rc, shortname=shortname, rmass=rmass)


      if (shortname .eq. "MXAER") then

        !interpolate into carma bin
        sizedist_carmabin = 0._r8

        do ibin_local = 1, NBIN
          ! Calculate the radius assuming that all the mass will be emitted as this
          ! element.
          vrfact = ((3._r8/2._r8 / PI / (vmrat_MXAER + 1._r8))**(1._r8 / 3._r8)) * ((vmrat_MXAER**(1._r8 / 3._r8)) - 1._r8)
          robc(ibin_local) = (3._r8 * rmass(ibin_local) / 4._r8 / PI / rho_obc)**(1._r8 / 3._r8)
          drobc(ibin_local) = vrfact * ((rmass(ibin_local)/rho_obc) **(1._r8 / 3._r8))

          if(robc(ibin_local) .lt. size_aeronet(1)) then
            sizedist_carmabin(ibin_local) = sizedist_avg(1)
          end if
          if(robc(ibin_local) .ge. size_aeronet(aeronet_dim1)) then
            sizedist_carmabin(ibin_local) = sizedist_avg(aeronet_dim1)
          end if
          do isizebin= 1,aeronet_dim1-1
            if( robc(ibin_local) .ge. size_aeronet(isizebin) .and.  robc(ibin_local) .lt. size_aeronet(isizebin+1))then
              sizedist_carmabin(ibin_local) = sizedist_avg(isizebin)*(size_aeronet(isizebin+1)-robc(ibin_local))/&
                  (size_aeronet(isizebin+1)-size_aeronet(isizebin))&
                  +sizedist_avg(isizebin+1)*(robc(ibin_local)-size_aeronet(isizebin))&
                  /(size_aeronet(isizebin+1)-size_aeronet(isizebin))
            end if
          end do
        end do

        rm(:) = 0._r8
        do ibin_local = 1, NBIN
          rm(ibin_local) = sizedist_carmabin(ibin_local)*drobc(ibin_local)/robc(ibin_local)*RHO_obc*1.e-15_r8         ! kg
        enddo

        do ibin_local = 1, NBIN
          aeronet_fraction(ibin_local) = rm(ibin_local)/sum(rm(:))
        end do

      end if
    end do

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

    type(carma_type), intent(in)  :: carma      !! the carma object
    integer,          intent(in)  :: ielem      !! element index
    integer,          intent(in)  :: ibin       !! bin index
    real(r8),         intent(in)  :: latvals(:) !! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) !! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    !! Only initialize where .true.
    real(r8),         intent(inout) :: q(:,:)     !! mass mixing ratio (gcol, lev)
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
    use camsrfexch, only: cam_out_t

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
    use ppgrid,        only: pcols
    use physics_types, only: physics_state
    use phys_grid,     only: get_lon_all_p, get_lat_all_p
    use camsrfexch,    only: cam_in_t

    type(carma_type), intent(in)       :: carma                 !! the carma object
    integer, intent(in)                :: ibin                  !! bin index
    type(physics_state), intent(in)    :: state                 !! physics state
    real(r8), intent(in)               :: r                     !! bin center (cm)
    real(r8), intent(in)               :: dr                    !! bin width (cm)
    real(r8), intent(in)               :: rmass                 !! bin mass (g)
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    real(r8), intent(out)              :: SaltFlux(pcols)       !! constituent surface flux (kg/m^2/s)
    integer, intent(out)               :: rc                    !! return code, negative indicates failure

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
    real(r8), parameter :: c4   = -1.424_r8        ! constants in calculating the particle wet radius

    ! Default return code.
    rc = RC_OK

    ! Determine the latitude and longitude of each column.
    lchnk = state%lchnk
    ncol = state%ncol

    ! Add any surface flux here.
    SaltFlux(:ncol) = 0.0_r8

    ! Are we configured for one of the known emission schemes?
    if( carma_seasalt_emis .ne. "Gong"       .and. &
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
    r80cm   = (c1 *  (r) ** c2 / (c3 * r ** c4 - log10(0.8_r8)) + (r)**3) ** (1._r8/3._r8) ! [cm]
    rdrycm  = r  ! [cm]
    r80     = r80cm *1.e4_r8    ! [um]
    rdry    = rdrycm*1.e4_r8  ! [um]

    do icol = 1,ncol

       ! Only generate sea salt over the ocean.
       if (cam_in%ocnfrac(icol) > 0._r8) then

          !**********************************
          !    WIND for seasalt production
          !**********************************
          call CARMA_SurfaceWind_salt(icol, cam_in, u10in, rc)

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

          ! ustar_smith = cd_smith **0.5_r8 * u10in
          !
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
             Monahan = 1.373_r8 * (u10in**3.41_r8) * r80**(-3._r8) * (1._r8 + 0.057_r8 *r80**1.05_r8)  * 10._r8 ** (1.19_r8 * exp(-1._r8 * B_mona**2)) ! dF/dr

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
          SaltFlux(icol) = ncflx * dr * rmass * 10._r8      ! *1e4[um/cm] * 1.e-3[kg/g]

          !          if (do_print) write(LUNOPRT, *) "ibin = ", ibin, ", igroup = ", igroup
          !          if (do_print) write(LUNOPRT, *) "dr = ", dr, ", rmass = ", rmass
          !          if (do_print) write(LUNOPRT, *) "ncflx = " , ncflx, ", SaltFlux = ", SaltFlux(icol)

          ! weighted by the ocean fraction
          SaltFlux(icol) = SaltFlux(icol) * cam_in%ocnfrac(icol)
       end if
    end do

  contains

    ! Coefficient Ak in Martensson's source functions
    pure real(r8) function Ak1(rpdry)
      real(r8),intent(in) :: rpdry
      Ak1 = c41*(2._r8*rpdry)**4 + c31*(2._r8*rpdry) ** 3 + c21*(2._r8*rpdry)**2 + c11*(2._r8*rpdry)+ c01
    end function Ak1

    pure real(r8) function Ak2(rpdry)
      real(r8),intent(in) :: rpdry
      Ak2 = c42*(2._r8*rpdry)**4 + c32*(2._r8*rpdry) ** 3 + c22*(2._r8*rpdry)**2 + c12*(2._r8*rpdry)+ c02
    end function Ak2

    pure real(r8) function Ak3(rpdry)
      real(r8),intent(in) :: rpdry
      Ak3 = c43*(2._r8*rpdry)**4 + c33*(2._r8*rpdry) ** 3 + c23*(2._r8*rpdry)**2 + c13*(2._r8*rpdry)+ c03
    end function Ak3

    ! Coefficient Bk in Martensson's source functions
    pure real(r8) function Bk1(rpdry)
      real(r8),intent(in) :: rpdry
      Bk1= d41*(2._r8*rpdry)**4 + d31*(2._r8*rpdry) ** 3 + d21*(2._r8*rpdry)**2 + d11*(2._r8*rpdry)+ d01
    end function Bk1

    pure real(r8) function Bk2(rpdry)
      real(r8),intent(in) :: rpdry
      Bk2 = d42*(2._r8*rpdry)**4 + d32*(2._r8*rpdry) ** 3 + d22*(2._r8*rpdry)**2 + d12*(2._r8*rpdry)+ d02
    end function Bk2

    pure real(r8) function Bk3(rpdry)
      real(r8),intent(in) :: rpdry
      Bk3 = d43*(2._r8*rpdry)**4 + d33*(2._r8*rpdry) ** 3 + d23*(2._r8*rpdry)**2 + d13*(2._r8*rpdry)+ d03
    end function Bk3

    ! Coefficient Ak in Clarkes's source function
    pure real(r8) function A1(rpdry)
      real(r8),intent(in) :: rpdry
      A1 = beta01 + beta11*(2._r8*rpdry) + beta21*(2._r8*rpdry)**2 + beta31*(2._r8*rpdry)**3 &
           + beta41*(2._r8*rpdry)**4 + beta51*(2._r8*rpdry)**5
    end function A1

    pure real(r8) function A2(rpdry)
      real(r8),intent(in) :: rpdry
      A2 = beta02 + beta12*(2._r8*rpdry) + beta22*(2._r8*rpdry)**2 + beta32*(2._r8*rpdry)**3 &
           + beta42*(2._r8*rpdry)**4 + beta52*(2._r8*rpdry)**5
    end function A2

    pure real(r8) function A3(rpdry)
      real(r8),intent(in) :: rpdry
      A3 = beta03 + beta13*(2._r8*rpdry) + beta23*(2._r8*rpdry)**2 + beta33*(2._r8*rpdry)**3 &
           + beta43*(2._r8*rpdry)**4 + beta53*(2._r8*rpdry)**5
    end function A3

  end subroutine CARMA_SaltFlux


  !! Calculate the sea surface wind with a Weibull distribution.
  !!
  !! @author  Tianyi Fan
  !! @version August-2010
  subroutine CARMA_SurfaceWind_salt(icol, cam_in, u10in, rc)
    use camsrfexch, only: cam_in_t

    ! in and out field
    integer, intent(in)                 :: icol                  !! column index
    type(cam_in_t), intent(in)          :: cam_in                !! surface inputs
    real(r8), intent(out)               :: u10in                 !! the 10m wind speed put into the source function
    integer, intent(out)                :: rc                    !! return code, negative indicates failure

    ! local variables
    real(r8) :: uWB341              ! the nth mean wind with integration using Weibull Distribution(integrate from threshold wind velocity)

    rc = RC_OK

    uWB341 = 0._r8

    ! calc. the Weibull wind distribution
    u10in = cam_in%u10(icol)

    call WeibullWind(u10in, uth_salt, 3.41_r8, uWB341)

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
  subroutine CARMA_ClayMassFraction(carma, igroup, rdust, rc)

    type(carma_type), intent(in)         :: carma       !! the carma object
    integer, intent(in)                  :: igroup      !! the carma group index
    real(r8), intent(in)                 :: rdust(NBIN) !! radius assuming entire particle is dust
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

    ! Figure out how many of the CARMA bins are in each of the Tegen and Lacis
    ! ranges.
    tl_count(:) = 0

    do ibin = 1, NBIN

      ! Smaller than the range.
      if (rdust(ibin) < tl_rmin(1)) then
        tl_count(IBELOW) = tl_count(IBELOW) + 1
      end if

      ! In the range
      do j = 1, NBIN_TEGEN
        if (rdust(ibin) < tl_rmax(j) .and. rdust(ibin) >= tl_rmin(j)) then
          tl_count(j+1) = tl_count(j+1) + 1
        end if
      end do

      ! Bigger than the range.
      if (rdust(ibin) >= tl_rmax(NBIN_TEGEN)) then
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
  subroutine CARMA_SurfaceWind(carma, icol, ielem, igroup, ibin, cam_in, uv10, wwd, uth, rc)
    use camsrfexch, only: cam_in_t

    ! in and out field
    type(carma_type), intent(in)        :: carma                 !! the carma object
    integer, intent(in)                 :: icol                  !! column index
    integer, intent(in)                 :: ielem                 !! element index
    integer, intent(in)                 :: igroup                !! group index
    integer, intent(in)                 :: ibin                  !! bin index
    type(cam_in_t), intent(in)          :: cam_in                !! surface inputs
    real(r8), intent(out)               :: uv10                  !! the 10m wind speed (m/s)
    real(r8), intent(out)               :: wwd                   !! the 10m wind speed  with Weibull applied (m/s)
    real(r8), intent(out)               :: uth                   !! the 10m wind threshold (m/s)
    integer,  intent(inout)             :: rc                    !! return code, negative indicates failure

    real(r8), parameter                 :: vk = 0.4_r8           ! von Karman constant
    real(r8)                            :: rmass(NBIN)           ! CARMA bin mass (g)
    real(r8)                            :: r                     ! CARMA bin center (cm)
    real(r8)                            :: rhop(NBIN)            ! CARMA partile element density (g/cm3)
    real(r8)                            :: uthfact               !
    real(r8), parameter                 :: rhoa = 1.25e-3_r8     ! Air density at surface

    rc = RC_OK

    ! Get the 10 meter wind speed
    uv10 = cam_in%u10(icol)

    ! Calculate the threshold wind speed of each bin [Marticorena and Bergametti,1995]
    ! note that in cgs units --> m/s
    call CARMAGROUP_GET(carma, igroup, rc, rmass=rmass)
    if (RC < RC_ERROR) return

    ! Define particle # concentration element index for current group
    call CARMAELEMENT_Get(carma, ielem, rc, rho=rhop)
    if (RC < RC_ERROR) return

    ! Calculate the radius assuming that all the mass will be emitted as this
    ! element.
    r = (3._r8 * rmass(ibin) / 4._r8 / PI / rhop(ibin))**(1._r8 / 3._r8)

    if (cam_in%soilw(icol) >= 0._r8 .AND. cam_in%soilw(icol) < 0.5_r8) then

       ! Prevent small values of soilw from driving uthfact negative, but allow
       ! for dust emissions even when soilw is 0.
       uthfact = 1.2_r8 + 0.2_r8*log10(max(0.001_r8, cam_in%soilw(icol)))

       if (r > 2.825e-5_r8) then  ! r(4) = 2.825e-5 cm
           uth = uthfact * 1.e-2_r8 * 0.13_r8 * sqrt(rhop(ibin)*GRAV*r*2._r8/rhoa) &
                       * sqrt(1._r8 + .006_r8/rhop(ibin)/GRAV/(r*2._r8)**2.5_r8) &
                       / sqrt(1.928_r8*(1331._r8*(r*2._r8)**1.56_r8 + .38_r8)**.092_r8 - 1._r8)
       else
           uth = uthfact*1.e-2_r8* 0.13_r8 * sqrt(rhop(ibin)*GRAV*(.75e-4_r8)*2._r8/rhoa)   &
                       * sqrt(1._r8 + .006_r8/rhop(ibin)/GRAV/((.75e-4_r8)*2._r8)**2.5_r8) &
                       / sqrt(1.928_r8*(1331._r8*((.75e-4_r8)*2._r8)**1.56_r8 + .38_r8)**.092_r8 - 1._r8)
       endif
    else
       uth = uv10
    endif

    ! Use Weibull with Lansing's estimate for shape.
    call WeibullWind(uv10, uth, 2._r8, wwd)

    ! Set the threshold to the weibull wind value if sol moisture >= 0.5,
    ! to turn off emissions.
    if (cam_in%soilw(icol) >= 0.5_r8) then
      uth = sqrt(wwd)
    end if

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
  subroutine CARMA_ReadSoilErosionFactor(rc)
    use ppgrid,             only: begchunk, endchunk, pcols
    use ioFileMod,     	    only: getfil
    use interpolate_data,   only: lininterp_init, lininterp, interp_type, lininterp_finish
    use phys_grid,          only: get_rlon_all_p, get_rlat_all_p, get_ncols_p
    use wrap_nf

    integer, intent(out)                      :: rc                    !! return code, negative indicates failure

    ! local variables
    integer                                   :: idvar, f_nlon, f_nlat, idlat, idlon
    integer                                   :: fid, fid_lon, fid_lat
    real(r8), allocatable, dimension(:,:)     :: ero_factor
    character(len=256)                        :: ero_file
    real(r8), allocatable, dimension(:)       :: ero_lat               ! latitude dimension
    real(r8), allocatable, dimension(:)       :: ero_lon               ! latitude dimension
    type (interp_type)                        :: lat_wght, lon_wght
    real(r8)                                  :: lat(pcols)            ! latitude index
    real(r8)                                  :: lon(pcols)            ! longitude index
    integer                                   :: i, ii
    integer                                   :: lchnk                 ! chunk identifier
    integer                                   :: ncol                  ! number of columns in chunk

    real(r8), parameter   :: zero=0_r8, twopi=2_r8*pi, degs2rads = pi/180._r8

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
    allocate(soil_factor(pcols, begchunk:endchunk))

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

    ero_lat(:) = ero_lat(:)*degs2rads
    ero_lon(:) = ero_lon(:)*degs2rads

    ! Close the file.
    call wrap_close(fid)

    do lchnk=begchunk, endchunk
       ncol = get_ncols_p(lchnk)

       call get_rlat_all_p(lchnk, pcols, lat)
       call get_rlon_all_p(lchnk, pcols, lon)

       call lininterp_init(ero_lon, f_nlon, lon, ncol, 2, lon_wght, zero, twopi)
       call lininterp_init(ero_lat, f_nlat, lat, ncol, 1, lat_wght)

       call lininterp(ero_factor, f_nlon, f_nlat, soil_factor(1:ncol,lchnk), ncol, lon_wght, lat_wght)

       call lininterp_finish(lon_wght)
       call lininterp_finish(lat_wght)
    end do

    deallocate(ero_lat)
    deallocate(ero_lon)
    deallocate(ero_factor)

  end subroutine CARMA_ReadSoilErosionFactor

  !! Calculate the nth mean of u using Weibull wind distribution
  !! considering the threshold wind velocity. This algorithm
  !! integrates from uth to infinite (u^n P(u)du )
  !!
  !! @author  Tianyi Fan
  !! @version August-2010
   subroutine WeibullWind(u, uth, n, uwb, wbk)
    use shr_spfn_mod, only: gamma =>  shr_spfn_gamma, igamma => shr_spfn_igamma

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
      k = 0.94_r8*u**0.5_r8        ! follow Grini and Zender, 2004JGR
 !    k = 2.5_r8                   ! Lansing's estimate
    end if

    ! If u is 0, then k can be 0, which makes a lot of this undefined.
    ! Just return 0. in this case.
    if (u < 0.35_r8) then
      uwb = 0._r8
    else
      c   = u * (gamma(1._r8 + 1._r8 / k))**(-1._r8)
      uwb = c**n * igamma(n / k + 1._r8, (uth / c)**k)
    end if

  end subroutine WeibullWind

  !! Read BC data from three components:
  !! 1. GAINS anthropogenic; 2. Ship Emission; 3. GFEDv3; 4. Aircraft
  !! GAINS unit: kt/year; 2D; lon:-180-180
  !! Ship Emission unit: kg/m2/s; 3D (month,lat,lon); lon:0-360
  !! GFEDv3 unit: g/m2/month; 3D (month,lat,lon); lon:-180-180
  !!
  !! @author  Pengfei Yu
  !! @version May-2013
  subroutine CARMA_BCOCRead(rc)
    use pmgrid,        only: plat, plon
    use ioFileMod,     only: getfil
    use cam_pio_utils, only: cam_pio_openfile
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish
    use pio,            only : file_desc_t, var_desc_t, &
                               pio_inq_dimid, pio_inq_varid, &
                               pio_get_var, pio_nowrite, pio_inq_dimlen, &
                               pio_inq_dimlen, pio_closefile
    use dycore,        only: dycore_is

    integer, intent(out)                      :: rc                    !! return code, negative indicates failure

    ! local variables
    integer                                   :: f_nlon, f_nlat, f_ntime
    integer                                   :: fid_lon, fid_lat, fid_time
    real(r8), allocatable, dimension(:,:)     :: BC_f2d, BC2d, OC_f2d, OC2d
    real(r8), allocatable, dimension(:,:,:)   :: BC_f3d, BC3d, OC_f3d, OC3d
!
    character(len=256)                        :: BC_GAINS_file
    character(len=256)                        :: OC_GAINS_file
    character(len=256)                        :: BC_GFEDv3_file
    character(len=256)                        :: OC_GFEDv3_file
    character(len=256)                        :: BC_ship_file
    character(len=256)                        :: OC_ship_file
!
    real(r8), allocatable, dimension(:,:,:)       :: BC_anthro_GAINS
    real(r8), allocatable, dimension(:,:,:)       :: OC_anthro_GAINS
    real(r8), allocatable, dimension(:,:,:)       :: BC_GFEDv3
    real(r8), allocatable, dimension(:,:,:)       :: OC_GFEDv3
    real(r8), allocatable, dimension(:,:,:)       :: BC_ship_GAINS
    real(r8), allocatable, dimension(:,:,:)       :: OC_ship_GAINS
!
    real(r8), allocatable, dimension(:)       :: BC_lat, OC_lat       ! latitude dimension
    real(r8), allocatable, dimension(:)       :: BC_lon, OC_lon       ! latitude dimension
    type (interp_type)                        :: wgt1, wgt2
    real(r8)                                  :: lat(plat), lon(plon)
    integer                                   :: i, itime
    real(r8)                                                              :: rearth, gridarea
    integer                                                                       :: nmonth
    real(r8)                                                              :: tempor(plon,plat)
    real(r8), allocatable, dimension(:,:,:)       :: tempor3d
    real(r8), allocatable, dimension(:,:)         :: tempor2d
    real(r8), allocatable, dimension(:)           :: tempor1d
    integer                                                                       :: mid_idx
    real(r8), allocatable, dimension(:,:)         :: BC_dom_f2d, OC_dom_f2d
    real(r8), allocatable, dimension(:,:,:)       :: BC_dom_f3d, OC_dom_f3d
    real(r8), allocatable, dimension(:,:,:)       :: BC_awb_f3d, OC_awb_f3d
    real(r8), allocatable, dimension(:,:)         :: BC2d_dom, OC2d_dom
    real(r8), allocatable, dimension(:)           :: facH, facL
    integer                                                                       :: ind_15N, ind_45N, ierr
    type(file_desc_t) :: fid
    type(var_desc_t) :: idvar, idlat, idlon, idvar_dom, idvar_awb

    real(r8) :: nlats

    rc = RC_OK

    if(dycore_is('UNSTRUCTURED') ) then
       call endrun('CARMA_InitializeModel: Yu2015 emissions not implemented for unstructured grids' )
    end if

    ! get model lat and lon
    nlats = plat-1 ! gnu compiler workaround
    do i = 1, plat
       lat(i) = 180._r8/(nlats)*(i-1)-90._r8
    end do
    do i = 1, plon
       lon(i) = 360._r8/plon*(i-1)
    end do

!
    nmonth = 12

    if(carma_BCOCemissions == 'Yu2015')then
       ! allocate BCnew and OCnew, unit is #/cm2/s
       allocate(BCnew(plat, plon, nmonth))
       allocate(OCnew(plat, plon, nmonth))
       BCnew = -huge(1._r8)
       OCnew = -huge(1._r8)
    endif

! monthly fraction of domestic emission
    allocate(facH(nmonth))
    allocate(facL(nmonth))
    facH = (/0.18_r8,0.14_r8,0.13_r8,0.08_r8,0.04_r8,0.02_r8,0.01_r8,&
                0.02_r8,0.03_r8,0.07_r8,0.11_r8,0.17_r8/)
    facL = (/0.17_r8,0.14_r8,0.11_r8,0.06_r8,0.04_r8,0.04_r8,0.04_r8,&
                0.04_r8,0.04_r8,0.06_r8,0.10_r8,0.15_r8/)

! find index for 15N and 45N
    do i = 1, plat
       if (lat(i) .gt. 15._r8) then
          ind_15N = i
          exit
       endif
    end do
!
    do i = 1, plat
       if (lat(i) .gt. 45._r8) then
          ind_45N = i
          exit
       endif
    end do

    ! Part 1a: BC anthropogenic from GAINS
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(BC_GAINS_filename, BC_GAINS_file, 0)
    call cam_pio_openfile( fid, BC_GAINS_file, PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'time', fid_time)
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_time,f_ntime)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(BC_lat(f_nlat))
    allocate(BC_lon(f_nlon))
    allocate(BC_f3d(f_nlon, f_nlat, f_ntime))
    allocate(BC_f2d(f_nlon, f_nlat))
    allocate(BC_dom_f2d(f_nlon, f_nlat))
    allocate(BC_dom_f3d(f_nlon, f_nlat, f_ntime))
    allocate(BC_awb_f3d(f_nlon, f_nlat, f_ntime))
    allocate(BC2d (plon, plat))
    allocate(BC2d_dom (plon, plat))
    allocate(BC_anthro_GAINS(nmonth, plat, plon))

    ! Read in the tables.
    ierr = pio_inq_varid(fid, 'emis_all', idvar)
    ierr = pio_get_var(fid, idvar, BC_f3d )
    ierr = pio_inq_varid(fid, 'emis_dom', idvar_dom)
    ierr = pio_get_var(fid, idvar, BC_dom_f3d )
    ierr = pio_inq_varid(fid, 'emis_awb', idvar_awb)
    ierr = pio_get_var(fid, idvar, BC_awb_f3d )
    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, BC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, BC_lon )

    ! Close the file.
    call pio_closefile(fid)
    ! get emission excluding domestic and agriculture waste buring
    BC_f2d = BC_f3d(:,:,1) - BC_dom_f3d(:,:,1) - BC_awb_f3d(:,:,1)
    BC_dom_f2d = BC_dom_f3d(:,:,1)

    ! make sure file longitude range from 0-360
    if (BC_lon(1) < -160._r8) then
       allocate(tempor2d(f_nlon, f_nlat))
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       ! emission excluding dom
       tempor2d(1:mid_idx,:f_nlat) = BC_f2d(mid_idx+1:f_nlon,:f_nlat)
       tempor1d(1:mid_idx) = BC_lon(mid_idx+1:f_nlon)
       tempor2d(mid_idx+1:f_nlon,:f_nlat) = BC_f2d(1:mid_idx,:f_nlat)
       tempor1d(mid_idx+1:f_nlon) = BC_lon(1:mid_idx)+360._r8
       BC_f2d = tempor2d
       ! dom emission
       tempor2d(1:mid_idx,:f_nlat) = BC_dom_f2d(mid_idx+1:f_nlon,:f_nlat)
       tempor2d(mid_idx+1:f_nlon,:f_nlat) = BC_dom_f2d(1:mid_idx,:f_nlat)
       BC_dom_f2d = tempor2d
       !
       BC_lon = tempor1d
       deallocate(tempor2d)
       deallocate(tempor1d)
    else
       BC_lon = BC_lon
    endif

    ! Convert kt/year ----> #/cm2/s
    rearth = 6.371e6_r8 ! m
    do i = 1, f_nlat
       gridarea = 2.0_r8*3.14159_r8*rearth/f_nlat * &
                          2.0_r8*3.14159_r8*rearth/f_nlon*cos(BC_lat(i)/180._r8*3.14159_r8)
       !
       BC_f2d(:f_nlon,i) = BC_f2d(:f_nlon,i)/365._r8/86400._r8*1.e9_r8/ &       ! g/s
                                        12._r8*6.02e23_r8/gridarea*1.e-4_r8                     ! #/cm2/s
       !
       BC_dom_f2d(:f_nlon,i) = BC_dom_f2d(:f_nlon,i)/365._r8/86400._r8*1.e9_r8/ &       ! g/s
                                        12._r8*6.02e23_r8/gridarea*1.e-4_r8                     ! #/cm2/s
    end do

    call lininterp_init(BC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(BC_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(BC_f2d, f_nlon, f_nlat, BC2d, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    call lininterp_init(BC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(BC_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(BC_dom_f2d, f_nlon, f_nlat, BC2d_dom, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    ! To implement Monthly data for dom emssion
    ! methods from Stohl et al., 2013
    ! facH works for high latitudes: 45-90N
    ! facL works for low latitudes: 15-45N
    ! below 15N, no seasonal variation
    !
    do itime = 1, nmonth
       ! 45N-90N
       BC2d(:plon, ind_45N:plat) = BC2d(:plon, ind_45N:plat) + &
                                   BC2d_dom(:plon, ind_45N:plat)*facH(itime)*12._r8
       ! 15N-45N
       BC2d(:plon, ind_15N:ind_45N-1) = BC2d(:plon, ind_15N:ind_45N-1) + &
                                        BC2d_dom(:plon, ind_15N:ind_45N-1)*facL(itime)*12._r8
       ! 90S-15N
       BC2d(:plon, 1:ind_15N-1) = BC2d(:plon, 1:ind_15N-1) + &
                                  BC2d_dom(:plon, 1:ind_15N-1)

       BC_anthro_GAINS(itime, :plat, :plon) = transpose(BC2d(:plon, :plat))
    end do

    deallocate(BC_lat)
    deallocate(BC_lon)
    deallocate(BC_f2d)
    deallocate(BC_f3d)
    deallocate(BC_dom_f2d)
    deallocate(BC_dom_f3d)
    deallocate(BC_awb_f3d)
    deallocate(BC2d)
    deallocate(BC2d_dom)

    ! Part 1b: OC anthropogenic from GAINS
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(OC_GAINS_filename, OC_GAINS_file, 0)
    call cam_pio_openfile(fid, trim(OC_GAINS_file), PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'time', fid_time)
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_time,f_ntime)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(OC_lat(f_nlat))
    allocate(OC_lon(f_nlon))
    allocate(OC_f2d(f_nlon, f_nlat))
    allocate(OC_f3d(f_nlon, f_nlat, f_ntime))
    allocate(OC_dom_f2d(f_nlon, f_nlat))
    allocate(OC_dom_f3d(f_nlon, f_nlat, f_ntime))
    allocate(OC_awb_f3d(f_nlon, f_nlat, f_ntime))
    allocate(OC2d (plon, plat))
    allocate(OC2d_dom (plon, plat))
    allocate(OC_anthro_GAINS(nmonth, plat, plon))

    ! Read in the tables.
    ierr = pio_inq_varid(fid, 'emis_all', idvar)
    ierr = pio_get_var(fid, idvar, OC_f3d )
    ierr = pio_inq_varid(fid, 'emis_dom', idvar_dom)
    ierr = pio_get_var(fid, idvar, OC_dom_f3d )
    ierr = pio_inq_varid(fid, 'emis_awb', idvar_awb)
    ierr = pio_get_var(fid, idvar, OC_awb_f3d )
    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, OC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, OC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! get emission excluding domestic and agriculture waste burning
    OC_f2d(:,:) = OC_f3d(:,:,1) - OC_dom_f3d(:,:,1) - OC_awb_f3d(:,:,1)
    OC_dom_f2d = OC_dom_f3d(:,:,1)

    ! make sure file longitude range from -180-180 to 0-360
    if (OC_lon(1) < -160._r8) then
       allocate(tempor2d(f_nlon, f_nlat))
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       ! emission excluding dom
       tempor2d(1:mid_idx,:f_nlat) = OC_f2d(mid_idx+1:f_nlon,:f_nlat)
       tempor1d(1:mid_idx) = OC_lon(mid_idx+1:f_nlon)
       tempor2d(mid_idx+1:f_nlon,:f_nlat) = OC_f2d(1:mid_idx,:f_nlat)
       tempor1d(mid_idx+1:f_nlon) = OC_lon(1:mid_idx)+360._r8
       OC_f2d = tempor2d
       ! dom emission
       tempor2d(1:mid_idx,:f_nlat) = OC_dom_f2d(mid_idx+1:f_nlon,:f_nlat)
       tempor2d(mid_idx+1:f_nlon,:f_nlat) = OC_dom_f2d(1:mid_idx,:f_nlat)
       OC_dom_f2d = tempor2d
       !
       OC_lon = tempor1d
       deallocate(tempor2d)
       deallocate(tempor1d)
    else
       OC_lon = OC_lon
    endif

    ! Convert kt/year ----> #/cm2/s
    rearth = 6.371e6_r8 ! m
    do i = 1, f_nlat
       gridarea = 2.0_r8*3.14159_r8*rearth/f_nlat * &
                  2.0_r8*3.14159_r8*rearth/f_nlon*cos(OC_lat(i)/180._r8*3.14159_r8)
       !
       OC_f2d(:f_nlon,i) = OC_f2d(:f_nlon,i)/365._r8/86400._r8*1.e9_r8/ &       ! g/s
                           12._r8*6.02e23_r8/gridarea*1.e-4_r8                     ! #/cm2/s
       !
       OC_dom_f2d(:f_nlon,i) = OC_dom_f2d(:f_nlon,i)/365._r8/86400._r8*1.e9_r8/ &       ! g/s
                               12._r8*6.02e23_r8/gridarea*1.e-4_r8                     ! #/cm2/s
    end do

    call lininterp_init(OC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(OC_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(OC_f2d, f_nlon, f_nlat, OC2d, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    call lininterp_init(OC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(OC_lon, f_nlon, lon, plon, 1, wgt2)
    call lininterp(OC_dom_f2d, f_nlon, f_nlat, OC2d_dom, plon, plat, wgt2, wgt1)
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    ! To implement Monthly data for dom emssion
    ! methods from Stohl et al., 2013
    ! facH works for high latitudes: 45-90N
    ! facL works for low latitudes: 15-45N
    ! below 15N, no seasonal variation
    !
    do itime = 1, nmonth
       ! 45N-90N
       OC2d(:plon, ind_45N:plat) = OC2d(:plon, ind_45N:plat) + &
                                   OC2d_dom(:plon, ind_45N:plat)*facH(itime)*12._r8
       ! 15N-45N
       OC2d(:plon, ind_15N:ind_45N-1) = OC2d(:plon, ind_15N:ind_45N-1) + &
                                        OC2d_dom(:plon, ind_15N:ind_45N-1)*facL(itime)*12._r8
       ! 90S-15N
       OC2d(:plon, 1:ind_15N-1) = OC2d(:plon, 1:ind_15N-1) + &
                                  OC2d_dom(:plon, 1:ind_15N-1)

       OC_anthro_GAINS(itime, :plat, :plon) = transpose(OC2d(:plon, :plat))
    end do

    deallocate(OC_lat)
    deallocate(OC_lon)
    deallocate(OC_f2d)
    deallocate(OC_f3d)
    deallocate(OC_dom_f2d)
    deallocate(OC_dom_f3d)
    deallocate(OC_awb_f3d)
    deallocate(OC2d)
    deallocate(OC2d_dom)

    ! Part 2a: BC ship
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(BC_ship_filename, BC_ship_file, 0)
    call cam_pio_openfile(fid, trim(BC_ship_file), PIO_NOWRITE)
    !call wrap_open(BC_ship_file, 0, fid)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(BC_lat(f_nlat))
    allocate(BC_lon(f_nlon))
    allocate(BC_f3d(f_nlon, f_nlat, nmonth))
    allocate(BC3d (plon, plat, nmonth))
    allocate(BC_ship_GAINS(nmonth, plat, plon))

   ! Read in the tables.
    ierr = pio_inq_varid(fid, 'emiss_shp', idvar)
    ierr = pio_get_var(fid, idvar, BC_f3d )
    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, BC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, BC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! make sure file longitude range from -180-180 to 0-360
    if (BC_lon(1) < -160._r8) then
       allocate(tempor3d(f_nlon, f_nlat, nmonth))
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       tempor3d(1:mid_idx,:f_nlat,:nmonth) = BC_f3d(mid_idx+1:f_nlon,:f_nlat,:nmonth)
       tempor1d(1:mid_idx) = BC_lon(mid_idx+1:f_nlon)
       tempor3d(mid_idx+1:f_nlon,:f_nlat,:nmonth) = BC_f3d(1:mid_idx,:f_nlat,:nmonth)
       tempor1d(mid_idx+1:f_nlon) = BC_lon(1:mid_idx)+360._r8
       BC_f3d = tempor3d
       BC_lon = tempor1d
       deallocate(tempor3d)
       deallocate(tempor1d)
    else
       BC_lon = BC_lon
    endif

    ! convert unit from kg/m2/s to #/cm2/s
    BC_f3d = BC_f3d*1.e3_r8/1.e4_r8/12._r8*6.02e23_r8

    call lininterp_init(BC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(BC_lon, f_nlon, lon, plon, 1, wgt2)
    do itime = 1, nmonth
       call lininterp(BC_f3d(:,:,itime), f_nlon, f_nlat, tempor(:,:), plon, plat, wgt2, wgt1)
       BC3d(:,:,itime) = tempor(:,:)
    end do
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    do itime = 1, nmonth
       BC_ship_GAINS(itime, :plat, :plon) = transpose(BC3d(:plon, :plat, itime))
    end do

    deallocate(BC_lat)
    deallocate(BC_lon)
    deallocate(BC_f3d)
    deallocate(BC3d)

    ! Part 2b: OC Ship
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(OC_ship_filename, OC_ship_file, 0)
    call cam_pio_openfile(fid, trim(OC_ship_file), PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(OC_lat(f_nlat))
    allocate(OC_lon(f_nlon))
    allocate(OC_f3d(f_nlon, f_nlat, nmonth))
    allocate(OC3d (plon, plat, nmonth))
    allocate(OC_ship_GAINS(nmonth, plat, plon))

    ! Read in the tables.
    ierr = pio_inq_varid(fid, 'emiss_shp', idvar)
    ierr = pio_get_var(fid, idvar, OC_f3d )
    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, OC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, OC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! make sure file longitude range from -180-180 to 0-360
    if (OC_lon(1) < -160._r8) then
       allocate(tempor3d(f_nlon, f_nlat, nmonth))
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       tempor3d(1:mid_idx,:f_nlat,:nmonth) = OC_f3d(mid_idx+1:f_nlon,:f_nlat,:nmonth)
       tempor1d(1:mid_idx) = OC_lon(mid_idx+1:f_nlon)
       tempor3d(mid_idx+1:f_nlon,:f_nlat,:nmonth) = OC_f3d(1:mid_idx,:f_nlat,:nmonth)
       tempor1d(mid_idx+1:f_nlon) = OC_lon(1:mid_idx)+360._r8
       OC_f3d = tempor3d
       OC_lon = tempor1d
       deallocate(tempor3d)
       deallocate(tempor1d)
    else
       OC_lon = OC_lon
    endif

    ! convert unit from kg/m2/s to #/cm2/s
    OC_f3d = OC_f3d*1.e3_r8/1.e4_r8/12._r8*6.02e23_r8

    call lininterp_init(OC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(OC_lon, f_nlon, lon, plon, 1, wgt2)
    do itime = 1, nmonth
       call lininterp(OC_f3d(:,:,itime), f_nlon, f_nlat, tempor(:,:), plon, plat, wgt2, wgt1)
       OC3d(:,:,itime) = tempor(:,:)
    end do
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    do itime = 1, nmonth
       OC_ship_GAINS(itime, :plat, :plon) = transpose(OC3d(:plon, :plat, itime))
    end do

    deallocate(OC_lat)
    deallocate(OC_lon)
    deallocate(OC_f3d)
    deallocate(OC3d)

    ! Part 3a: BC GFEDv3
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(BC_GFEDv3_filename, BC_GFEDv3_file, 0)
    call cam_pio_openfile(fid, trim(BC_GFEDv3_file), PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    allocate(BC_lat(f_nlat))
    allocate(BC_lon(f_nlon))
    allocate(BC_f3d(f_nlon, f_nlat, nmonth))
    allocate(tempor3d(f_nlon, f_nlat, nmonth))
    allocate(BC3d (plon, plat, nmonth))
    allocate(BC_GFEDv3(nmonth, plat, plon))

    ! Read in the tables.
    BC_f3d = 0._r8
    ierr = pio_inq_varid(fid, 'emis', idvar)
    ierr = pio_get_var(fid, idvar, tempor3d )
    !call wrap_inq_varid(fid, 'emis', idvar)
    !call wrap_get_var_realx(fid, idvar,  tempor3d)
    BC_f3d = BC_f3d + tempor3d
    ! excluding non-real values
    where (BC_f3d(:,:,:) .ge. 1.e10_r8)
        BC_f3d(:,:,:) = 1.e-30_r8
    end where

    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, BC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, BC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! make sure file longitude range from -180-180 to 0-360
    if (BC_lon(1) < -160._r8) then
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       tempor3d(1:mid_idx,:f_nlat,:nmonth) = BC_f3d(mid_idx+1:f_nlon,:f_nlat,:nmonth)
       tempor1d(1:mid_idx) = BC_lon(mid_idx+1:f_nlon)
       tempor3d(mid_idx+1:f_nlon,:f_nlat,:nmonth) = BC_f3d(1:mid_idx,:f_nlat,:nmonth)
       tempor1d(mid_idx+1:f_nlon) = BC_lon(1:mid_idx)+360._r8
       BC_f3d = tempor3d
       BC_lon = tempor1d
       deallocate(tempor1d)
    else
       BC_lon = BC_lon
    endif

    ! convert unit from g/m2/month to #/cm2/s
    BC_f3d = BC_f3d/1.e4_r8/30._r8/86400._r8/12._r8*6.02e23_r8

    call lininterp_init(BC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(BC_lon, f_nlon, lon, plon, 1, wgt2)
    do itime = 1, nmonth
       call lininterp(BC_f3d(:,:,itime), f_nlon, f_nlat, tempor(:,:), plon, plat, wgt2, wgt1)
       BC3d(:,:,itime) = tempor(:,:)
    end do
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    do itime = 1, nmonth
       BC_GFEDv3(itime, :plat, :plon) = transpose(BC3d(:plon, :plat, itime))
    end do

    deallocate(BC_lat)
    deallocate(BC_lon)
    deallocate(BC_f3d)
    deallocate(BC3d)
    deallocate(tempor3d)

    ! Part 3b: OC GFEDv3
    ! -------------------------------------------------
    ! Open the netcdf file (read only)
    call getfil(OC_GFEDv3_filename, OC_GFEDv3_file, 0)
    call cam_pio_openfile(fid, trim(OC_GFEDv3_file), PIO_NOWRITE)

    ! Get file dimensions
    ierr = pio_inq_dimid(fid, 'lon', fid_lon)
    ierr = pio_inq_dimid(fid, 'lat', fid_lat)
    ierr = pio_inq_dimlen(fid, fid_lon, f_nlon)
    ierr = pio_inq_dimlen(fid, fid_lat, f_nlat)

    ! write(carma%f_LUNOPRT,*) ''
    ! write(carma%f_LUNOPRT,*) 'f_lon = ', f_nlon
    ! write(carma%f_LUNOPRT,*) 'f_lat = ', f_nlat
    ! write(carma%f_LUNOPRT,*) ''

    allocate(OC_lat(f_nlat))
    allocate(OC_lon(f_nlon))
    allocate(OC_f3d(f_nlon, f_nlat, nmonth))
    allocate(tempor3d(f_nlon, f_nlat, nmonth))
    allocate(OC3d (plon, plat, nmonth))
    allocate(OC_GFEDv3(nmonth, plat, plon))

    ! Read in the tables.
     OC_f3d = 0._r8
    ierr = pio_inq_varid(fid, 'emis', idvar)
    ierr = pio_get_var(fid, idvar, tempor3d )
    !call wrap_inq_varid(fid, 'emis', idvar)
    !call wrap_get_var_realx(fid, idvar,  tempor3d)
    OC_f3d = OC_f3d + tempor3d
    ! excluding non-real values
    where (OC_f3d(:,:,:) .ge. 1.e10_r8)
        OC_f3d(:,:,:) = 1.e-30_r8
    end where

    ierr = pio_inq_varid(fid, 'lat', idlat)
    ierr = pio_get_var(fid, idlat, OC_lat )
    ierr = pio_inq_varid(fid, 'lon ', idlon)
    ierr = pio_get_var(fid, idlon, OC_lon )

    ! Close the file.
    call pio_closefile(fid)

    ! make sure file longitude range from -180-180 to 0-360
    if (OC_lon(1) < -160._r8) then
       allocate(tempor1d(f_nlon))
       mid_idx = floor(f_nlon/2._r8)
       tempor3d(1:mid_idx,:f_nlat,:nmonth) = OC_f3d(mid_idx+1:f_nlon,:f_nlat,:nmonth)
       tempor1d(1:mid_idx) = OC_lon(mid_idx+1:f_nlon)
       tempor3d(mid_idx+1:f_nlon,:f_nlat,:nmonth) = OC_f3d(1:mid_idx,:f_nlat,:nmonth)
       tempor1d(mid_idx+1:f_nlon) = OC_lon(1:mid_idx)+360._r8
       OC_f3d = tempor3d
       OC_lon = tempor1d
       deallocate(tempor1d)
    else
       OC_lon = OC_lon
    endif
    call lininterp_init(OC_lat, f_nlat, lat, plat, 1, wgt1)
    call lininterp_init(OC_lon, f_nlon, lon, plon, 1, wgt2)
    do itime = 1, nmonth
       call lininterp(OC_f3d(:,:,itime), f_nlon, f_nlat, tempor(:,:), plon, plat, wgt2, wgt1)
       OC3d(:,:,itime) = tempor(:,:)
    end do
    call lininterp_finish(wgt1)
    call lininterp_finish(wgt2)

    do itime = 1, nmonth
       OC_GFEDv3(itime, :plat, :plon) = transpose(OC3d(:plon, :plat, itime))
    end do

    deallocate(OC_lat)
    deallocate(OC_lon)
    deallocate(OC_f3d)
    deallocate(OC3d)
    deallocate(tempor3d)

! Sum
    do itime = 1, nmonth
       BCnew(:plat, :plon, itime) = BC_anthro_GAINS(itime, :plat, :plon) +  &
             BC_ship_GAINS(itime, :plat, :plon) +  BC_GFEDv3(itime, :plat, :plon)
!
       OCnew(:plat, :plon, itime) = OC_anthro_GAINS(itime, :plat, :plon) +  &
             OC_ship_GAINS(itime, :plat, :plon) +  OC_GFEDv3(itime, :plat, :plon)
    end do
!
    deallocate(BC_anthro_GAINS)
    deallocate(OC_anthro_GAINS)
    deallocate(BC_ship_GAINS)
    deallocate(OC_ship_GAINS)
    deallocate(BC_GFEDv3)
    deallocate(OC_GFEDv3)
    deallocate(facH)
    deallocate(facL)
!
    return
  end subroutine CARMA_BCOCRead

end module carma_model_mod
