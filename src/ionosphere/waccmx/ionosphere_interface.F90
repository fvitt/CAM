module ionosphere_interface

  use shr_kind_mod,        only: r8 => shr_kind_r8
  use phys_grid,           only: begchunk, endchunk, get_ncols_p
  use pmgrid,              only: plat, plon, plev
  use ppgrid,              only: pcols, pver

  use dpie_coupling,       only: d_pie_init
  use dpie_coupling,       only: d_pie_epotent
  use dpie_coupling,       only: d_pie_coupling         ! WACCM-X ionosphere/electrodynamics coupling
  use short_lived_species, only: slvd_index,slvd_pbf_ndx => pbf_idx ! Routines to access short lived species 

  use chem_mods,           only: adv_mass      ! Array holding mass values for short lived species
  use mo_chem_utls,        only: get_spc_ndx   ! Routine to get index of adv_mass array for short lived species
  use physics_buffer,      only: pbuf_get_chunk, pbuf_get_field, pbuf_get_index

  use cam_abortutils,      only: endrun
  use constituents,        only: cnst_get_ind, cnst_mw  !Needed to access constituent molecular weights
  use phys_grid,           only: get_lon_all_p, get_lat_all_p, transpose_block_to_chunk, transpose_chunk_to_block
  use phys_grid,           only: chunk_to_block_send_pters, chunk_to_block_recv_pters, block_to_chunk_send_pters, &
                                 block_to_chunk_recv_pters
  use physconst,           only: gravit
  use oplus,               only: oplus_init
  use vxb,                 only: vxb_init   !vxB jianfei Wu added
  use edyn_init,           only: edynamo_init
  use pio,                 only: var_desc_t
  use spmd_dyn,            only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
  use dyn_internal_state,  only: get_dyn_state_grid
  use dynamics_vars,       only: t_fvdycore_grid
  use perf_mod
  use epotential_params,   only: epot_active, epot_crit_colats
  use spmd_utils,          only: masterproc
  use cam_logfile,         only: iulog

  implicit none

  private

  public :: ionosphere_readnl
  public :: ionosphere_init
  public :: ionosphere_run1
  public :: ionosphere_run2
  public :: ionosphere_init_restart
  public :: ionosphere_write_restart
  public :: ionosphere_read_restart
  public :: ionosphere_final

  ! private data 

  ! this needs to persist from time-step to time-step and across restarts
  real(r8), allocatable :: opmmrtm1_blck(:,:,:)   ! O+ at previous time step(blocks)
  real(r8), allocatable :: nopmmrtm1_blck(:,:,:)   ! O+ at previous time step(blocks)
  real(r8), allocatable :: o2pmmrtm1_blck(:,:,:)   ! O+ at previous time step(blocks)
  !-------------Jianfei Wu----------------------------------------------
  real(r8), allocatable :: fepmmrtm1_blck(:,:,:)   ! Fe+ at previous time step(blocks)
  real(r8), allocatable :: mgpmmrtm1_blck(:,:,:)   ! Fe+ at previous time step(blocks)
  real(r8), allocatable :: napmmrtm1_blck(:,:,:)   ! Fe+ at previous time step(blocks)

  type(var_desc_t) :: Feptm1_vdesc 
  type(var_desc_t) :: Mgptm1_vdesc 
  type(var_desc_t) :: Naptm1_vdesc 
!Add more metal ions, Wuhu Feng, 11 April 2022

    real(r8), allocatable :: capmmrtm1_blck(:,:,:)   ! Ca+ at previous time step(blocks)
    real(r8), allocatable :: kpmmrtm1_blck(:,:,:)   ! K+ at previous time step(blocks)
    real(r8), allocatable :: sipmmrtm1_blck(:,:,:)   ! Si+ at previous time step(blocks)
 
    type(var_desc_t) :: captm1_vdesc
    type(var_desc_t) :: kptm1_vdesc
    type(var_desc_t) :: siptm1_vdesc
!-------------------------------------------------------
  type(var_desc_t) :: Optm1_vdesc 
  type(var_desc_t) :: NOptm1_vdesc 
  type(var_desc_t) :: O2ptm1_vdesc 
  integer :: index_ped, index_hall, index_te, index_ti
  integer :: index_ui, index_vi, index_wi

  integer :: ixo2=-1, ixo=-1, ixh=-1
  integer :: ixo2p=-1, ixnop=-1, ixn2p=-1, ixop=-1
  integer :: ixfep=-1  !Jianfei Wu added for Fe+ Na+ and K+
  integer :: ixmgp=-1  !Jianfei Wu added for Fe+ Na+ and K+
  integer :: ixnap=-1  !Jianfei Wu added for Fe+ Na+ and K+

  ! indices for accessing ions in pbuf when non-advected
  integer :: sIndxOp=-1, sIndxO2p=-1, sIndxNOp=-1, sIndxN2p=-1  
  integer :: sIndxFep=-1 !Jianfei Wu added for Fe+ Na+ and K+ 
  integer :: sIndxMgp=-1 !Jianfei Wu added for Fe+ Na+ and K+ 
  integer :: sIndxNap=-1 !Jianfei Wu added for Fe+ Na+ and K+ 
!Add more metal ions, Wuhu Feng, 11 April 2022

    integer :: ixcap=-1  !Jianfei Wu added for Fe+ Na+ and K+
    integer :: ixkp=-1  !Jianfei Wu added for Fe+ Na+ and K+
    integer :: ixsip=-1  !Jianfei Wu added for Fe+ Na+ and K+
    integer :: sIndxcap=-1 !Jianfei Wu added for Fe+ Na+ and K+
    integer :: sIndxkp=-1 !Jianfei Wu added for Fe+ Na+ and K+
    integer :: sIndxsip=-1 !Jianfei Wu added for Fe+ Na+ and K+
    real(r8) :: rmasscap    ! Fe+ molecular weight kg/kmol
    real(r8) :: rmasskp    ! Fe+ molecular weight kg/kmol
    real(r8) :: rmasssip    ! Fe+ molecular weight kg/kmol



  real(r8) :: rmassO2    ! O2 molecular weight kg/kmol
  real(r8) :: rmassO1    ! O atomic weight kg/kmol
  real(r8) :: rmassH     ! H atomic weight kg/kmol
  real(r8) :: rmassN2    ! N2 molecular weight kg/kmol
  real(r8) :: rmassO2p   ! O2+ molecular weight kg/kmol
  real(r8) :: rmassNOp   ! NO+ molecular weight kg/kmol
  real(r8) :: rmassN2p   ! N2+ molecular weight kg/kmol
  real(r8) :: rmassOp    ! O+ molecular weight kg/kmol
  !-------------Jianfei Wu added---------------------
  real(r8) :: rmassFep    ! Fe+ molecular weight kg/kmol
  real(r8) :: rmassMgp    ! Fe+ molecular weight kg/kmol
  real(r8) :: rmassNap    ! Fe+ molecular weight kg/kmol
  !-------------end------------------------------------

  logical, public,  protected :: ionos_edyn_active = .true.   ! if true, edynamo will generate ion drifts
  logical, public,  protected :: ionos_xport_active = .true.  ! if true, call d_pie_coupling from dp_coupling
  !
  ! ionos_edyn_active = .true. will activate the edynamo which will generate ion drift velocities 
  !  used in oplus transport, otherwise empirical ion drifts calculated in exbdrift (physics) will be used.
  !
  logical, public,  protected :: ionos_oplus_xport = .true.    ! if true, call sub oplus (based on tiegcm oplus.F)
  integer, public,  protected :: ionos_xport_nsplit = 5        ! number of substeps for O+ transport per model time step
  integer, public,  protected :: vxb_xport_nsplit = 15       ! number of substeps for O+ transport per model time step

  real(r8), public, protected :: oplus_adiff_limiter = 1.5e+8_r8  ! limiter for ambipolar diffusion coefficient
  real(r8), public, protected :: oplus_shapiro_const = 0.03_r8    ! shapiro constant for spatial smoother
  logical,  public, protected :: oplus_enforce_floor = .true.     ! switch to apply Stan's  floor

  character(len=256) :: wei05_coefs_file = 'NONE' !'wei05sc.nc'
  character(len=256) :: amienh_file  = 'NONE'
  character(len=256) :: amiesh_file  = 'NONE'

  character(len=16), public, protected :: ionos_epotential_model = 'none'
  logical,           public, protected :: ionos_epotential_amie = .false.
  integer ::  indxAMIEefxg=-1, indxAMIEkevg=-1

contains

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_readnl( nlfile )

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: mpicom, masterprocid, mpi_real8, mpi_logical, mpi_integer, mpi_character
    use cam_logfile,    only: iulog
    use spmd_utils,     only: masterproc

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'ionosphere_readnl'

    namelist /ionosphere_nl/ ionos_xport_active, ionos_edyn_active, ionos_oplus_xport, ionos_xport_nsplit, vxb_xport_nsplit
    namelist /ionosphere_nl/ oplus_adiff_limiter, oplus_shapiro_const, oplus_enforce_floor
    namelist /ionosphere_nl/ ionos_epotential_model, ionos_epotential_amie, wei05_coefs_file
    namelist /ionosphere_nl/ amienh_file, amiesh_file, wei05_coefs_file
    namelist /ionosphere_nl/ epot_crit_colats

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'ionosphere_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, ionosphere_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(ionos_xport_active,  1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_edyn_active,   1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_oplus_xport,   1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_xport_nsplit,  1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(vxb_xport_nsplit,  1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(oplus_adiff_limiter, 1, mpi_real8,   masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_epotential_model, len(ionos_epotential_model), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_epotential_amie,1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(wei05_coefs_file, len(wei05_coefs_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(amienh_file, len(amienh_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(amiesh_file, len(amiesh_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(oplus_shapiro_const, 1, mpi_real8,   masterprocid, mpicom, ierr)
    call mpi_bcast(oplus_enforce_floor, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(epot_crit_colats,    2, mpi_real8,   masterprocid, mpicom, ierr)

    ! log the user settings
    if (masterproc) then
       write(iulog,*) 'ionosphere_readnl: ionos_xport_active  = ', ionos_xport_active
       write(iulog,*) 'ionosphere_readnl: ionos_edyn_active   = ', ionos_edyn_active
       write(iulog,*) 'ionosphere_readnl: ionos_oplus_xport   = ', ionos_oplus_xport
       write(iulog,*) 'ionosphere_readnl: ionos_xport_nsplit  = ', ionos_xport_nsplit
       write(iulog,*) 'ionosphere_readnl: vxb_xport_nsplit  = ', vxb_xport_nsplit
       write(iulog,*) 'ionosphere_readnl: ionos_epotential_model = ', trim(ionos_epotential_model)
       write(iulog,*) 'ionosphere_readnl: ionos_epotential_amie  = ', ionos_epotential_amie
       write(iulog,'(a,2(g12.4))') &
                     ' ionosphere_readnl: epot_crit_colats       = ', epot_crit_colats
       write(iulog,*) 'ionosphere_readnl: oplus_adiff_limiter = ', oplus_adiff_limiter
       write(iulog,*) 'ionosphere_readnl: oplus_shapiro_const = ', oplus_shapiro_const
       write(iulog,*) 'ionosphere_readnl: oplus_enforce_floor = ', oplus_enforce_floor
    endif
    epot_active = .true.

  end subroutine ionosphere_readnl

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_init()
    use physics_buffer, only: pbuf_add_field, dtype_r8
    use cam_history,    only: addfld, add_default, horiz_only
    use mo_apex,        only: mo_apex_init1
    use cam_control_mod,only: initial_run
    use dyn_grid,       only: get_horiz_grid_d
    use ref_pres,  only : & ! Hybrid level definitions:
      pref_mid,           & ! target alev(plev) midpoint levels coord
      pref_edge             ! target ailev(plevp) interface levels coord
    use amie_module,    only: init_amie
    use wei05sc,        only: weimer05_init

    ! local variables:
    type (t_fvdycore_grid), pointer :: grid
    integer :: sIndx

    integer :: mpicomm         ! MPI communicator
    integer :: ntaski, ntaskj  ! number of MPI tasks in lon,lat dimensions
    integer :: lat0,lat1       ! first and last latitude  indices
    integer :: lon0,lon1       ! first and last longitude indices
    integer :: lev0,lev1       ! first and last pressure indices
    real(r8), allocatable :: glon(:) ! global geo-graphic longitudes (degrees)
    real(r8), allocatable :: glat(:) ! global geo-graphic latitudes (degrees)

    if ( ionos_epotential_amie ) then
       call pbuf_add_field('AMIE_efxg', 'global', dtype_r8, (/pcols/), indxAMIEefxg)  ! Energy flux from AMIE
       call pbuf_add_field('AMIE_kevg', 'global', dtype_r8, (/pcols/), indxAMIEkevg)  ! Mean energy from AMIE  
    endif
    if (initial_run) then
       call ionosphere_read_ic()
    endif

    call mo_apex_init1()

    op_transport: if (ionos_xport_active) then

       grid => get_dyn_state_grid()

       index_ped  = pbuf_get_index('PedConduct')
       index_hall = pbuf_get_index('HallConduct')

       index_te   = pbuf_get_index('TElec')
       index_ti   = pbuf_get_index('TIon')
       !
       ! pbuf indices to empirical ion drifts, to be passed to oplus_xport,
       ! if ionos_edyn_active is false.
       !
       index_ui   = pbuf_get_index('UI')
       index_vi   = pbuf_get_index('VI')
       index_wi   = pbuf_get_index('WI')

       !-----------------------------------------------------------------------
       !  Get indices for neutrals to get mixing ratios from state%q and masses
       !-----------------------------------------------------------------------
       call cnst_get_ind('O2' ,ixo2 )
       call cnst_get_ind('O'  ,ixo )
       call cnst_get_ind('H'  ,ixh )
       !------------------------------------
       ! Get neutral molecular weights
       !------------------------------------
       rmassO2 = cnst_mw(ixo2)
       rmassO1 = cnst_mw(ixo)
       rmassH  = cnst_mw(ixh)
       rmassN2 = 28._r8

       call cnst_get_ind('Op',ixop, abort=.false.)
       if (ixop > 0) then
          rMassOp = cnst_mw(ixop)
          if (masterproc) write(iulog,"('ixop:',i4,' rMassOp',f8.3)") ixop,rMassOp
       else
          sIndxOp  = slvd_index( 'Op' )
          if (sIndxOp > 0) then
             sIndx = get_spc_ndx( 'Op' )
             rmassOp = adv_mass(sIndx)
             if (masterproc) write(iulog,"('SIndxOp:',i4,' rMassOp',f8.3)") SIndxOp,rMassOp
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for Op')
          endif
       endif
      !-------------------------------------------------
      !-------------Jianfei Wu added--------------------
      !--------------------------------------------------
       call cnst_get_ind('Fep',ixfep, abort=.false.)
       if (ixfep > 0) then
          rMassFep = cnst_mw(ixfep)
          if (masterproc) write(iulog,"('ixfep:',i4,' rMassFep',f8.3)") ixfep,rMassFep
       else
          sIndxFep  = slvd_index( 'Fep' )
          if (sIndxFep > 0) then
             sIndx = get_spc_ndx( 'Fep' )
             rmassFep = adv_mass(sIndx)
             if (masterproc) write(iulog,"('SIndxFep:',i4' rMassFep',f8.3)") SIndxFep,rMassFep
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for Fep')
          endif
       endif
       if (masterproc) write(iulog,"('SIndxFep:',i4' rMassFep',f8.3)") SIndxFep,rMassFep

       call cnst_get_ind('Mgp',ixmgp, abort=.false.)
       if (ixmgp > 0) then
          rMassMgp = cnst_mw(ixmgp)
          if (masterproc) write(iulog,"('ixmgp:',i4,' rMassMgp',f8.3)") ixmgp,rMassMgp
       else
          sIndxMgp  = slvd_index( 'Mgp' )
          if (sIndxMgp > 0) then
             sIndx = get_spc_ndx( 'Mgp' )
             rmassMgp = adv_mass(sIndx)
             if (masterproc) write(iulog,"('SIndxMgp:',i4' rMassMgp',f8.3)") SIndxMgp,rMassMgp
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for Mgp')
          endif
       endif
       if (masterproc) write(iulog,"('SIndxMgp:',i4' rMassMgp',f8.3)") SIndxMgp,rMassMgp

       call cnst_get_ind('Nap',ixnap, abort=.false.)
       if (ixnap > 0) then
          rMassNap = cnst_mw(ixnap)
          if (masterproc) write(iulog,"('ixnap:',i4,' rMassNap',f8.3)") ixnap,rMassNap
       else
          sIndxNap  = slvd_index( 'Nap' )
          if (sIndxNap > 0) then
             sIndx = get_spc_ndx( 'Nap' )
             rmassNap = adv_mass(sIndx)
             if (masterproc) write(iulog,"('SIndxNap:',i4' rMassNap',f8.3)") SIndxNap,rMassNap
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for Nap')
          endif
       endif
       if (masterproc) write(iulog,"('SIndxNap:',i4' rMassNap',f8.3)") SIndxNap,rMassNap

       !------------------Jianfei Wu end---------------------------------------------
!Add more metal ions, Wuhu Feng, 11 April 2022
         call cnst_get_ind('Cap',ixcap, abort=.false.)
         if (ixcap > 0) then
            rMassCap = cnst_mw(ixcap)
            if (masterproc) write(iulog,"('ixcap:',i4,' rMassFep',f8.3)") ixcap,rMassCap
         else
            sIndxCap  = slvd_index( 'Cap' )
            if (sIndxCap > 0) then
               sIndx = get_spc_ndx( 'Cap' )
               rmassCap = adv_mass(sIndx)
               if (masterproc) write(iulog,"('SIndxCap:',i4' rMassFep',f8.3)") SIndxCap,rMassCap
            else
               call endrun('ionosphere_init: Cannot find state or pbuf index for Cap')
            endif
         endif
         if (masterproc) write(iulog,"('SIndxCap:',i4' rMassCap',f8.3)") SIndxCap,rMassCap
 
         call cnst_get_ind('Kp',ixkp, abort=.false.)
         if (ixkp > 0) then
            rMassKp = cnst_mw(ixkp)
            if (masterproc) write(iulog,"('ixkp:',i4,' rMassKp',f8.3)") ixkp,rMassKp
         else
            sIndxKp  = slvd_index( 'Kp' )
            if (sIndxKp > 0) then
               sIndx = get_spc_ndx( 'Kp' )
               rmassKp = adv_mass(sIndx)
               if (masterproc) write(iulog,"('SIndxKp:',i4' rMassKp',f8.3)") SIndxKp,rMassKp
            else
               call endrun('ionosphere_init: Cannot find state or pbuf index for Kp')
            endif
         endif
         if (masterproc) write(iulog,"('SIndxKp:',i4' rMassKp',f8.3)") SIndxKp,rMassKp
 
         call cnst_get_ind('Sip',ixsip, abort=.false.)
         if (ixsip > 0) then
            rMassSip = cnst_mw(ixsip)
            if (masterproc) write(iulog,"('ixsip:',i4,' rMassSip',f8.3)") ixsip,rMassSip
         else
            sIndxSip  = slvd_index( 'Sip' )
            if (sIndxSip > 0) then
               sIndx = get_spc_ndx( 'Sip' )
               rmassSip = adv_mass(sIndx)
               if (masterproc) write(iulog,"('SIndxSip:',i4' rMassSip',f8.3)") SIndxSip,rMassSip
            else
               call endrun('ionosphere_init: Cannot find state or pbuf index for Sip')
            endif
         endif
         if (masterproc) write(iulog,"('SIndxSip:',i4' rMassSip',f8.3)") SIndxSip,rMassSip


       call cnst_get_ind('O2p',ixo2p, abort=.false.)
       if (ixo2p > 0) then
          rMassO2p = cnst_mw(ixo2p)
       else
          sIndxO2p  = slvd_index( 'O2p' )
          if (sIndxO2p > 0) then
             sIndx = get_spc_ndx( 'O2p' )
             rmassO2p = adv_mass(sIndx)
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for O2p')
          endif
       endif

       call cnst_get_ind('NOp',ixnop, abort=.false.)
       if (ixnop > 0) then
          rMassNOp = cnst_mw(ixnop)
       else
          sIndxNOp  = slvd_index( 'NOp' )
          if (sIndxNOp > 0) then
             sIndx = get_spc_ndx( 'NOp' )
             rmassNOp = adv_mass(sIndx)
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for NOp')
          endif
       endif

       call cnst_get_ind('N2p',ixn2p, abort=.false.)
       if (ixn2p > 0) then
          rMassN2p = cnst_mw(ixn2p)
       else
          sIndxN2p  = slvd_index( 'N2p' )
          if (sIndxN2p > 0) then
             sIndx = get_spc_ndx( 'N2p' )
             rmassN2p = adv_mass(sIndx)
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for N2p')
          endif
       endif

       call d_pie_init( ionos_edyn_active, ionos_oplus_xport, ionos_xport_nsplit , vxb_xport_nsplit, epot_crit_colats)
       if ( grid%iam < grid%npes_xy ) then
          
          allocate(glon(plon))
          allocate(glat(plat))
          call get_horiz_grid_d( plon, lon_d_out=glon )
          call get_horiz_grid_d( plat, lat_d_out=glat )

          mpicomm = grid%commxy
          lon0 = grid%ifirstxy ; lon1 = grid%ilastxy
          lat0 = grid%jfirstxy ; lat1 = grid%jlastxy
          lev0 = 1             ; lev1 = grid%km
          ntaski = grid%nprxy_x
          ntaskj = grid%nprxy_y

          call edynamo_init( mpicomm, plon, plat, plev, lon0,lon1,lat0,lat1,lev0,lev1, ntaski,ntaskj, &
                             glon, glat, pref_mid,pref_edge)
          call ionosphere_alloc()
          call oplus_init( oplus_adiff_limiter, oplus_shapiro_const, oplus_enforce_floor )
          call vxb_init( oplus_adiff_limiter, oplus_shapiro_const, oplus_enforce_floor ) !Jianfei Wu

          deallocate(glon,glat)
       endif

       call addfld ('OpTM1&IC', (/ 'lev' /),'I','kg/kg','O+ at time step minus 1',gridname='fv_centers')
       call add_default ('OpTM1&IC',0, 'I')
       !-----------------Jianfei Wu----------------
       call addfld ('NOpTM1&IC', (/ 'lev' /),'I','kg/kg','NO+ at time step minus 1',gridname='fv_centers')
       call add_default ('NOpTM1&IC',0, 'I')
       call addfld ('O2pTM1&IC', (/ 'lev' /),'I','kg/kg','O2+ at time step minus 1',gridname='fv_centers')
       call add_default ('O2pTM1&IC',0, 'I')

       call addfld ('MgpTM1&IC', (/ 'lev' /),'I','kg/kg','Mg+ at time step minus 1',gridname='fv_centers')
       call add_default ('MgpTM1&IC',0, 'I')
       !------------------------------------------------------------------------------
       call addfld ('NapTM1&IC', (/ 'lev' /),'I','kg/kg','Mg+ at time step minus 1',gridname='fv_centers')
       call add_default ('NapTM1&IC',0, 'I')
       !------------------------------------------------------------------------------

    endif op_transport

    if (ionos_edyn_active) then
       call addfld ('UI',(/ 'lev' /),'I','m/s', 'UI Zonal ion drift from edynamo') 
       call addfld ('VI',(/ 'lev' /),'I','m/s', 'VI Meridional ion drift from edynamo')
       call addfld ('WI',(/ 'lev' /),'I','m/s', 'WI Vertical ion drift from edynamo')
       call addfld ('UI&IC', (/ 'lev' /), 'I','m/s', 'Zonal ion drift velocity')
       call addfld ('VI&IC', (/ 'lev' /), 'I','m/s', 'Meridional ion drift velocity')
       call addfld ('WI&IC', (/ 'lev' /), 'I','m/s', 'Vertical ion drift velocity')
       call add_default ('UI&IC', 0, ' ')
       call add_default ('VI&IC', 0, ' ')
       call add_default ('WI&IC', 0, ' ')
    endif
    if ( ionos_epotential_amie ) then
       call init_amie(amienh_file,amiesh_file)
       call addfld ('amie_efx_phys',horiz_only,'I','mW/m2', 'AMIE energy flux') 
       call addfld ('amie_kev_phys',horiz_only,'I','keV'  , 'AMIE mean energy')
    end if
    if ( trim(ionos_epotential_model) == 'weimer' ) then
       call weimer05_init(wei05_coefs_file)
    endif

  end subroutine ionosphere_init

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_run1(pbuf2d)
    use physics_buffer, only: physics_buffer_desc
    use cam_history  , only: outfld, write_inithist
    use phys_grid,      only: get_ncols_p

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local vars
    integer :: i, j, k, lchnk  ! indices
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy, km, idim
    real(r8), allocatable :: tmp(:,:), tmp1(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    type(t_fvdycore_grid), pointer :: grid

    real(r8), pointer :: pbuf_amie_efxg(:)     ! Pointer to access AMIE energy flux in pbuf
    real(r8), pointer :: pbuf_amie_kevg(:)     ! Pointer to access AMIE mean energy in pbuf
    
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude in
    integer :: blksiz                ! number of columns in 2D block
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: iam, astat
    integer :: ib, ic, jc,ncol
    integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    real(r8), allocatable, dimension(:) :: bbuffer, cbuffer
    real(r8), allocatable :: amie_efxg(:,:) ! energy flux from AMIE
    real(r8), allocatable :: amie_kevg(:,:) ! characteristic mean energy from AMIE

    grid => get_dyn_state_grid()
    iam = grid%iam

    ifirstxy     = grid%ifirstxy
    ilastxy      = grid%ilastxy
    jfirstxy     = grid%jfirstxy
    jlastxy      = grid%jlastxy
    km           = grid%km

    if( write_inithist() .and. ionos_xport_active ) then

       allocate( tmp(ifirstxy:ilastxy,km) )
       !---------------Jinafei Wu added-------------------
       allocate( tmp1(ifirstxy:ilastxy,km) )
       !---------------------Jianfei end------------------

       idim = ilastxy - ifirstxy + 1
       do j = jfirstxy, jlastxy
          do k = 1, km
             do i = ifirstxy, ilastxy
                tmp(i,k) = opmmrtm1_blck(i,j,k)
             enddo
          enddo
          call outfld ('OpTM1&IC', tmp, idim, j) 
       enddo

       deallocate( tmp )
       deallocate( tmp1 )
    endif

    amie_active: if ( ionos_epotential_amie ) then
       allocate(amie_efxg(ifirstxy:ilastxy,jfirstxy:jlastxy))
       allocate(amie_kevg(ifirstxy:ilastxy,jfirstxy:jlastxy))

       ! data assimilated potential
       call d_pie_epotent( ionos_epotential_model, epot_crit_colats, &
                           i0=ifirstxy,i1=ilastxy,j0=jfirstxy,j1=jlastxy, &
                           efxg=amie_efxg,kevg=amie_kevg )

       ! transform to physics grid for aurora...

       ! blocks --> physics chunks

       blcks2phys_local: if (local_dp_map) then

          chnk_loop1 : do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)
             call get_lon_all_p(lchnk, ncol, lons)
             call get_lat_all_p(lchnk, ncol, lats)

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
             call pbuf_get_field(pbuf_chnk, indxAMIEefxg, pbuf_amie_efxg)
             call pbuf_get_field(pbuf_chnk, indxAMIEkevg, pbuf_amie_kevg)

             do i=1,ncol
                ic = lons(i)
                jc = lats(i)
                pbuf_amie_efxg(i) = amie_efxg(ic,jc)
                pbuf_amie_kevg(i) = amie_kevg(ic,jc)
             end do
             call outfld ( 'amie_efx_phys', pbuf_amie_efxg, pcols, lchnk )
             call outfld ( 'amie_kev_phys', pbuf_amie_kevg, pcols, lchnk )
          end do chnk_loop1

       else ! blcks2phys_local

          tsize = 2
          blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
          allocate( bpter(blksiz,0:km),stat=astat )
          allocate( bbuffer(tsize*block_buf_nrecs),stat=astat )
          allocate( cbuffer(tsize*chunk_buf_nrecs),stat=astat )

          if (iam < grid%npes_xy) then 
             call block_to_chunk_send_pters(iam+1,blksiz,pver+1,tsize,bpter)
          endif

          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)
                bbuffer(bpter(ib,0)+0) = amie_efxg(i,j)
                bbuffer(bpter(ib,0)+1) = amie_kevg(i,j)
             end do
          end do

          call transpose_block_to_chunk(tsize, bbuffer, cbuffer)

          chnk_loop2: do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
             call pbuf_get_field(pbuf_chnk, indxAMIEefxg, pbuf_amie_efxg)
             call pbuf_get_field(pbuf_chnk, indxAMIEkevg, pbuf_amie_kevg)
             call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)
             do i=1,ncol
                pbuf_amie_efxg(i) = cbuffer(cpter(i,0)+0)
                pbuf_amie_kevg(i) = cbuffer(cpter(i,0)+1)
             end do
             call outfld ( 'amie_efx_phys', pbuf_amie_efxg, pcols, lchnk )
             call outfld ( 'amie_kev_phys', pbuf_amie_kevg, pcols, lchnk )
          end do chnk_loop2

          deallocate(bpter)
          deallocate(bbuffer)
          deallocate(cbuffer)


       end if blcks2phys_local

       deallocate(amie_efxg,amie_kevg)

    else
       
       ! set cross tail potential before physics -- aurora uses weimer derived potential
       call d_pie_epotent( ionos_epotential_model, epot_crit_colats )

    end if amie_active

  end subroutine ionosphere_run1

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_run2( phys_state, dyn_in, pbuf2d )

    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    use dyn_comp,       only: dyn_import_t
    use cam_history,    only: outfld, write_inithist

    ! - pull some fields from pbuf and dyn_in
    ! - invoke ionosphere/electro-dynamics coupling
    ! - push some fields back to physics via pbuf...

    ! args
    type(physics_state),    intent(inout), target :: phys_state(begchunk:endchunk)
    type(dyn_import_t),  intent(inout) :: dyn_in  ! dynamics inputs
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local vars
    integer :: i,j,k, lchnk
    integer :: astat

    integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    real(r8), allocatable, dimension(:) :: bbuffer, cbuffer

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    real(r8), pointer :: sigma_ped_phys(:,:)  ! physics pointer to Pedersen Conductivity
    real(r8), pointer :: sigma_hall_phys(:,:) ! physics pointer fo Hall Conductivity
    real(r8), pointer :: te_phys(:,:)         ! te from pbuf
    real(r8), pointer :: ti_phys(:,:)         ! ti from pbuf
    real(r8), pointer :: mmrPO2p_phys(:,:)    ! Pointer to access O2+ in pbuf
    real(r8), pointer :: mmrPNOp_phys(:,:)    ! Pointer to access NO+ in pbuf
    real(r8), pointer :: mmrPN2p_phys(:,:)    ! Pointer to access N2+ in pbuf
    real(r8), pointer :: mmrPOp_phys(:,:)     ! Pointer to access O+ in pbuf
    !------------Jianfei Wu-------------------------------------------------
    !----------------------------------------------------------------------
!
! Empirical ion drifts from exbdrift (to be converted to blocked for dpie_coupling):
    real(r8), pointer :: ui_phys(:,:)         ! zonal ion drift from pbuf
    real(r8), pointer :: vi_phys(:,:)         ! meridional ion drift from pbuf
    real(r8), pointer :: wi_phys(:,:)         ! vertical ion drift from pbuf
    real(r8), dimension(:,:,:), pointer :: q

    real(r8), pointer :: n2pmmr_blck(:,:,:) => null()     ! N2+ (blocks)
    real(r8), pointer :: opmmr_blck(:,:,:)  => null()     ! O+ (blocks)

    real(r8), pointer :: tracer(:,:,:,:)
    real(r8), pointer :: u3s(:,:,:)
    real(r8), pointer :: v3s(:,:,:)
    real(r8), pointer :: pexy(:,:,:)

    real(r8), pointer :: phis(:,:)            ! surface geopotential

    real(r8), pointer :: o2mmr_blck(:,:,:)
    real(r8), pointer :: o1mmr_blck(:,:,:)
    real(r8), pointer :: h1mmr_blck(:,:,:)

    integer :: ib, ic, jc, ifirstxy, ilastxy, jfirstxy, jlastxy, km, ncol

    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: nSIons                        ! number of ions set to non-advected
    integer :: ibuffOp,ibuffO2p,ibuffNOp, ibuffN2p ! Buffer indices for non-advected ions
    integer :: ibuffFep ! Buffer indices for non-advected ions
    integer :: ibuffMgp ! Buffer indices for non-advected ions
    integer :: ibuffNap ! Buffer indices for non-advected ions

    integer :: blksiz                 ! number of columns in 2D block
    integer :: tsize                  ! amount of data per grid point passed to physics
    integer :: dsize                  ! amount of data per grid point passed to physics  Jianfei Wu
    integer :: iam

    real(r8), allocatable :: wuxy(:,:,:)
    real(r8), allocatable :: wvxy(:,:,:)
    real(r8), allocatable :: sigma_ped_blck (:,:,:)
    real(r8), allocatable :: sigma_hall_blck(:,:,:)
    real(r8), allocatable :: ti_blck(:,:,:)
    real(r8), allocatable :: te_blck(:,:,:)
    real(r8), allocatable :: zi_blck(:,:,:)
    real(r8), allocatable :: zm_blck(:,:,:)
    real(r8), allocatable :: ui_blck(:,:,:)
    real(r8), allocatable :: vi_blck(:,:,:)
    real(r8), allocatable :: wi_blck(:,:,:)
    real(r8), allocatable :: omega_blck(:,:,:)
    real(r8), allocatable :: tn_blck(:,:,:)
    !----------------Jianfei Wu added---------------------------------
    real(r8), allocatable :: fepmmr_blck(:,:,:)    ! Fe+ (blocks)
    real(r8), allocatable :: mgpmmr_blck(:,:,:)    ! Fe+ (blocks)
    real(r8), allocatable :: napmmr_blck(:,:,:)    ! Fe+ (blocks)
!Add more metal ions, Wuhu Feng, 11 April 2022
    real(r8), allocatable :: capmmr_blck(:,:,:)    ! Fe+ (blocks)
    real(r8), allocatable :: kpmmr_blck(:,:,:)    ! Fe+ (blocks)
    real(r8), allocatable :: sipmmr_blck(:,:,:)    ! Fe+ (blocks)
    real(r8), allocatable :: o2pmmr_blck(:,:,:)     ! O2+ (blocks)
    real(r8), allocatable :: nopmmr_blck(:,:,:)     ! NO+ (blocks)
    !--------------------------------------------------------------------

    type (t_fvdycore_grid), pointer :: grid

    ionos_cpl: if (ionos_xport_active) then 

       grid => get_dyn_state_grid()
       iam = grid%iam

       allocate( wuxy(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( wvxy(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( sigma_ped_blck (grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( sigma_hall_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( ti_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( te_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( zi_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( zm_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( ui_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( vi_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( wi_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( omega_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( tn_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( fepmmr_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( mgpmmr_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( napmmr_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
!Add more metal ions, Wuhu Feng, 11 April 2022
       allocate( capmmr_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( kpmmr_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( sipmmr_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( nopmmr_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( o2pmmr_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )

       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       phis   => dyn_in%phis

       tracer => dyn_in%tracer
       pexy   => dyn_in%pe

       u3s    => dyn_in%u3s
       v3s    => dyn_in%v3s

       if (iam < grid%npes_xy) then 
          call d2a3dijk( grid, u3s, v3s, wuxy, wvxy )
       endif

       if (sIndxOp>0) then
          allocate(opmmr_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
          if (astat /= 0) call endrun('ionos_intr_d_p_cplng: failed to allocate opmmr_blck')
       endif
       if (sIndxN2p>0) then
          allocate(n2pmmr_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
          if (astat /= 0) call endrun('ionos_intr_d_p_cplng: failed to allocate n2pmmr_blck')
       endif
       !-------------------------Jianfei Wu------------------------------------------
       !----------------------------------------------------------------------------------

       phys2blcks_local: if (local_dp_map) then

          do lchnk = begchunk,endchunk

             ncol = get_ncols_p(lchnk)
             call get_lon_all_p(lchnk, ncol, lons)
             call get_lat_all_p(lchnk, ncol, lats)
             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

             ! Get Pedersen and Hall conductivities:
             call pbuf_get_field(pbuf_chnk, index_ped,  sigma_ped_phys)
             call pbuf_get_field(pbuf_chnk, index_hall, sigma_hall_phys)
             do k=1,km
                do i=1,ncol
                   sigma_ped_blck(lons(i),lats(i),k) = sigma_ped_phys(i,k)
                   sigma_hall_blck(lons(i),lats(i),k) = sigma_hall_phys(i,k)
                end do
             enddo

             ! Get ion and electron temperatures 
             call pbuf_get_field(pbuf_chnk, index_te, te_phys)
             call pbuf_get_field(pbuf_chnk, index_ti, ti_phys)
             do k=1,km
                do i=1,ncol
                   te_blck(lons(i),lats(i),k) = te_phys(i,k)
                   ti_blck(lons(i),lats(i),k) = ti_phys(i,k)
                end do
             enddo

             ! Get components of ion drift velocities
             call pbuf_get_field(pbuf_chnk, index_ui, ui_phys)
             call pbuf_get_field(pbuf_chnk, index_vi, vi_phys)
             call pbuf_get_field(pbuf_chnk, index_wi, wi_phys)
             do k=1,km
                do i=1,ncol
                   ui_blck(lons(i),lats(i),k) = ui_phys(i,k)
                   vi_blck(lons(i),lats(i),k) = vi_phys(i,k)
                   wi_blck(lons(i),lats(i),k) = wi_phys(i,k)
                   zi_blck(lons(i),lats(i),k)    = phys_state(lchnk)%zi(i,k)
                   zm_blck(lons(i),lats(i),k)    = phys_state(lchnk)%zm(i,k)
                   omega_blck(lons(i),lats(i),k) = phys_state(lchnk)%omega(i,k)
                   tn_blck(lons(i),lats(i),k)    = phys_state(lchnk)%t(i,k)
                   fepmmr_blck(lons(i),lats(i),k)    = phys_state(lchnk)%q(i,k,ixfep)
                   mgpmmr_blck(lons(i),lats(i),k)    = phys_state(lchnk)%q(i,k,ixmgp)
                   napmmr_blck(lons(i),lats(i),k)    = phys_state(lchnk)%q(i,k,ixnap)
!add more metal ions, WUhu Feng 12/04/2022
                   capmmr_blck(lons(i),lats(i),k)    = phys_state(lchnk)%q(i,k,ixcap)
                   kpmmr_blck(lons(i),lats(i),k)    = phys_state(lchnk)%q(i,k,ixkp)
                   sipmmr_blck(lons(i),lats(i),k)    = phys_state(lchnk)%q(i,k,ixsip)
                enddo
             enddo

             !--------------------------------------------------------
             ! Get ions from physics buffer if non-transported
             !--------------------------------------------------------
             if (sIndxN2p > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPN2p_phys, &
                     start=(/1,1,sIndxN2p/), kount=(/pcols,pver,1/) )
                do k=1,km
                   do i=1,ncol
                      n2pmmr_blck(lons(i),lats(i),k) = mmrPN2p_phys(i,k)
                   end do
                enddo
             endif
             if (sIndxOp > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys, &
                     start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
                do k=1,km
                   do i=1,ncol
                      opmmr_blck(lons(i),lats(i),k) = mmrPOp_phys(i,k)
                   end do
                enddo
             endif
             !------------------------------Jianfei Wu------------------
             !---------------------------------------------------------------

          enddo ! do lchnk = begchunk,endchunk

       else ! phys2blcks_local

!Add more metal ions, WUhu Feng, 11/04/2022
!         tsize = 16
          tsize = 19

          nSIons = 0
          if (sIndxOp > 0)  then 
             ibuffOp = tsize + nSIons
             nSIons = nSIons + 1
          endif
          if (sIndxN2p > 0) then
             ibuffN2p = tsize + nSIons
             nSIons = nSIons + 1
          endif
          !-------------Jianfei Wu-------------
          !----------------------------------
          tsize = tsize + nSIons

          blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
          allocate(bpter(blksiz,0:km))
          allocate(bbuffer(tsize*block_buf_nrecs))
          allocate(cbuffer(tsize*chunk_buf_nrecs))

          do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

             ! Get Pedersen and Hall conductivities:
             call pbuf_get_field(pbuf_chnk, index_ped,  sigma_ped_phys)
             call pbuf_get_field(pbuf_chnk, index_hall, sigma_hall_phys)

             ! Get ion and electron temperatures 
             call pbuf_get_field(pbuf_chnk, index_te,  te_phys)
             call pbuf_get_field(pbuf_chnk, index_ti,  ti_phys)

             ! Get components of ion drift velocities
             call pbuf_get_field(pbuf_chnk, index_ui,  ui_phys)
             call pbuf_get_field(pbuf_chnk, index_vi,  vi_phys)
             call pbuf_get_field(pbuf_chnk, index_wi,  wi_phys)
 
             !--------------------------------------------------------
             ! Get ions from physics buffer if non-transported
             !--------------------------------------------------------

             if (sIndxOp > 0)  call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys,  &
                  start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
             if (sIndxN2p > 0) call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPN2p_phys, &
                  start=(/1,1,sIndxN2p/), kount=(/pcols,pver,1/) )
              !------------------------------Jianfei Wu-------------------------------

             call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

             do i=1,ncol
                cbuffer(cpter(i,0):cpter(i,0)+tsize-1) = 0.0_r8
             end do

             do k=1,km
                do i=1,ncol

                   cbuffer(cpter(i,k)+0) = sigma_ped_phys(i,k)
                   cbuffer(cpter(i,k)+1) = sigma_hall_phys(i,k)
                   cbuffer(cpter(i,k)+2) = te_phys(i,k)
                   cbuffer(cpter(i,k)+3) = ti_phys(i,k)
                   cbuffer(cpter(i,k)+4) = phys_state(lchnk)%zi(i,k)
                   cbuffer(cpter(i,k)+5) = phys_state(lchnk)%zm(i,k)
                   cbuffer(cpter(i,k)+6) = ui_phys(i,k)
                   cbuffer(cpter(i,k)+7) = vi_phys(i,k)
                   cbuffer(cpter(i,k)+8) = wi_phys(i,k)
                   cbuffer(cpter(i,k)+9) = phys_state(lchnk)%omega(i,k)
                   cbuffer(cpter(i,k)+10) = phys_state(lchnk)%t(i,k)
                   cbuffer(cpter(i,k)+11) = phys_state(lchnk)%q(i,k,ixfep)
                   cbuffer(cpter(i,k)+12) = phys_state(lchnk)%q(i,k,ixmgp)
                   cbuffer(cpter(i,k)+13) = phys_state(lchnk)%q(i,k,ixnap)

                   cbuffer(cpter(i,k)+14) = phys_state(lchnk)%q(i,k,ixnop)
                   cbuffer(cpter(i,k)+15) = phys_state(lchnk)%q(i,k,ixo2p)
!Add more metal ions, WUhu Feng, 11/04/2022
                   cbuffer(cpter(i,k)+16) = phys_state(lchnk)%q(i,k,ixcap)
                   cbuffer(cpter(i,k)+17) = phys_state(lchnk)%q(i,k,ixkp)
                   cbuffer(cpter(i,k)+18) = phys_state(lchnk)%q(i,k,ixsip)

                   if (sIndxN2p > 0)cbuffer(cpter(i,k)+ibuffN2p) = mmrPN2p_phys(i,k)
                   if (sIndxOp > 0) cbuffer(cpter(i,k)+ibuffOp)  = mmrPOp_phys(i,k)
                   !-----------------Jianfei Wu-----------------------------------

                end do

             end do

          end do

          call t_barrierf('sync_chk_to_blk', grid%commxy)
          call t_startf ('chunk_to_block')
          call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
          call t_stopf  ('chunk_to_block')

          if (iam < grid%npes_xy) then 
             call chunk_to_block_recv_pters(iam+1,blksiz,pver+1,tsize,bpter)
          endif

          do j=jfirstxy,jlastxy
             do k=1,km
                do i=ifirstxy,ilastxy
                   ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)

                   sigma_ped_blck(i,j,k)  = bbuffer(bpter(ib,k)+0)
                   sigma_hall_blck(i,j,k) = bbuffer(bpter(ib,k)+1)
                   te_blck(i,j,k)         = bbuffer(bpter(ib,k)+2)
                   ti_blck(i,j,k)         = bbuffer(bpter(ib,k)+3)
                   zi_blck(i,j,k)         = bbuffer(bpter(ib,k)+4)
                   zm_blck(i,j,k)         = bbuffer(bpter(ib,k)+5)
                   ui_blck(i,j,k)         = bbuffer(bpter(ib,k)+6)
                   vi_blck(i,j,k)         = bbuffer(bpter(ib,k)+7)
                   wi_blck(i,j,k)         = bbuffer(bpter(ib,k)+8)
                   omega_blck(i,j,k)      = bbuffer(bpter(ib,k)+9)
                   tn_blck(i,j,k)         = bbuffer(bpter(ib,k)+10)
                   fepmmr_blck(i,j,k)     = bbuffer(bpter(ib,k)+11)
                   mgpmmr_blck(i,j,k)     = bbuffer(bpter(ib,k)+12)
                   napmmr_blck(i,j,k)     = bbuffer(bpter(ib,k)+13)

                   nopmmr_blck(i,j,k) = bbuffer(bpter(ib,k)+14)
                   o2pmmr_blck(i,j,k) = bbuffer(bpter(ib,k)+15)

!WUHU FENG., extra metal ions
                   capmmr_blck(i,j,k)     = bbuffer(bpter(ib,k)+16)
                   kpmmr_blck(i,j,k)     = bbuffer(bpter(ib,k)+17)
                   sipmmr_blck(i,j,k)     = bbuffer(bpter(ib,k)+18)

                   if (sIndxN2p > 0) n2pmmr_blck(i,j,k) = bbuffer(bpter(ib,k)+ibuffN2p)
                   if (sIndxOp > 0)  opmmr_blck(i,j,k)  = bbuffer(bpter(ib,k)+ibuffOp)
                   !-----------------Jianfei Wu--------------------------------

                enddo
             enddo
          enddo

          deallocate(bpter)
          deallocate(bbuffer)
          deallocate(cbuffer)

       endif phys2blcks_local

       !-------------------------------------------------------------------------------------------
       !  Set dpie_coupling input ions if they are advected ...
       !-------------------------------------------------------------------------------------------
       if (ixn2p > 0) then
          n2pmmr_blck => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixn2p)
       endif
       if (ixop > 0) then
          opmmr_blck => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixop)
       endif
       !----------------------Jianfei Wu-----------------------------------------------

       !------------------------------------
       ! Get neutrals from advected tracers array
       !------------------------------------

       o2mmr_blck  => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixo2)
       o1mmr_blck  => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixo)
       h1mmr_blck  => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixh)

       !
       !   Make geopotential height (m) for d_pie_coupling. 
       !
       do k=1,km
          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                zi_blck(i,j,k) = zi_blck(i,j,k)+phis(i,j)/gravit ! phis is redundant in k
                zm_blck(i,j,k) = zm_blck(i,j,k)+phis(i,j)/gravit ! phis is redundant in k
             enddo
          enddo
       enddo

       call t_startf('d_pie_coupling')

       if (iam < grid%npes_xy) then 
          ! waccmx ionosphere electro-dynamics -- transports O+ and provides updates to ion drift velocities
          call d_pie_coupling(omega_blck,pexy,zi_blck,zm_blck,wuxy,wvxy,tn_blck,                        &
               sigma_ped_blck,sigma_hall_blck,te_blck,ti_blck,                      &
               o2mmr_blck,o1mmr_blck,h1mmr_blck,o2pmmr_blck,o2pmmrtm1_blck,nopmmr_blck,&
               nopmmrtm1_blck,n2pmmr_blck, opmmr_blck,opmmrtm1_blck,                   &
               fepmmr_blck,fepmmrtm1_blck,mgpmmr_blck,mgpmmrtm1_blck, napmmr_blck,napmmrtm1_blck,&
!WUHU FENG, add more metal ions,
               capmmr_blck,capmmrtm1_blck,kpmmr_blck,kpmmrtm1_blck, sipmmr_blck,sipmmrtm1_blck,&
               ui_blck,vi_blck,wi_blck,    & !Jianfei Wu
               rmassO2,rmassO1,rmassH,rmassN2,rmassO2p,rmassNOp,rmassN2p, rmassOp,rmassFep,rmassMgp,rmassNap, &
!WUHU FENG, add more metal ions,
               rmassCap,rmassKp,rmassSip, &
               ifirstxy,ilastxy, jfirstxy,jlastxy)
       endif

       call t_stopf ('d_pie_coupling')

       !
       !----------------------------------------
       !  Put data back in to state%q or pbuf
       !----------------------------------------
       if (ixop > 0) then
          tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixop) = opmmr_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km)           
       endif
       !-------------------Jianfei Wu-----------------------------------------------------------------------------
       !-------------------Jianfei Wu-----------------------------------------------------------------------------
       !-------------------Jianfei Wu-----------------------------------------------------------------------------
       !----------------------------------------------------------------------------------------------------------

       ! blocks --> physics chunks

       blcks2phys_local: if (local_dp_map) then

          chnk_loop1 : do lchnk = begchunk,endchunk
             ncol = phys_state(lchnk)%ncol
             call get_lon_all_p(lchnk, ncol, lons)
             call get_lat_all_p(lchnk, ncol, lats)

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

             call pbuf_get_field(pbuf_chnk, index_ui, ui_phys)
             call pbuf_get_field(pbuf_chnk, index_vi, vi_phys)
             call pbuf_get_field(pbuf_chnk, index_wi, wi_phys)
             if (sIndxOp > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys, &
                     start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
             endif
             !----------------------Jianfei Wu--------------------------
             !-----------------------------------------------------------
             q => phys_state(lchnk)%q
             do k=1,km
                do i=1,ncol
                   ic = lons(i)
                   jc = lats(i)
                   ui_phys(i,k) = ui_blck(ic,jc,k)
                   vi_phys(i,k) = vi_blck(ic,jc,k)
                   wi_phys(i,k) = wi_blck(ic,jc,k)
                   q(i,k,ixfep)=fepmmr_blck(ic,jc,k)
                   q(i,k,ixmgp)=mgpmmr_blck(ic,jc,k)
                   q(i,k,ixnap)=napmmr_blck(ic,jc,k)
!WUHU FENG, more metal ions
                     q(i,k,ixcap)=capmmr_blck(ic,jc,k)
                     q(i,k,ixkp)=kpmmr_blck(ic,jc,k)
                     q(i,k,ixsip)=sipmmr_blck(ic,jc,k)

                   q(i,k,ixnop) = nopmmr_blck(ic,jc,k)
                   q(i,k,ixo2p) = o2pmmr_blck(ic,jc,k)
                   if (sIndxOp > 0) mmrPOp_phys(i,k) = opmmr_blck(ic,jc,k)
                end do
             end do

             if (ionos_edyn_active) then
                call outfld ( 'UI', ui_phys, pcols, lchnk )
                call outfld ( 'VI', vi_phys, pcols, lchnk )
                call outfld ( 'WI', wi_phys, pcols, lchnk )
                if (write_inithist()) then
                   call outfld ( 'UI&IC', ui_phys, pcols, lchnk )
                   call outfld ( 'VI&IC', vi_phys, pcols, lchnk )
                   call outfld ( 'WI&IC', wi_phys, pcols, lchnk )
                endif
             endif

          end do chnk_loop1

       else ! blcks2phys_local

          if (sIndxOp > 0) then
!WUHU FENG, more metal ions
!            tsize = 9 ! for ui,vi,wi,op
             tsize = 12 ! for ui,vi,wi,cap,kp,sip,op
          else
!            tsize = 8 ! for ui,vi,wi
             tsize = 11 ! for ui,vi,wi
          endif
          tsize=tsize+1

          blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
          allocate( bpter(blksiz,0:km),stat=astat )
          allocate( bbuffer(tsize*block_buf_nrecs),stat=astat )
          allocate( cbuffer(tsize*chunk_buf_nrecs),stat=astat )

          if (iam < grid%npes_xy) then 
             call block_to_chunk_send_pters(iam+1,blksiz,km+1,tsize,bpter)
          endif

          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)

                do k=1,km

                   bbuffer(bpter(ib,k)) = ui_blck(i,j,k)
                   bbuffer(bpter(ib,k)+1) = vi_blck(i,j,k)
                   bbuffer(bpter(ib,k)+2) = wi_blck(i,j,k)
                   bbuffer(bpter(ib,k)+3) = fepmmr_blck(i,j,k)
                   bbuffer(bpter(ib,k)+4) = mgpmmr_blck(i,j,k)
                   bbuffer(bpter(ib,k)+5) = napmmr_blck(i,j,k)
                   bbuffer(bpter(ib,k)+6) = nopmmr_blck(i,j,k)
                   bbuffer(bpter(ib,k)+7) = o2pmmr_blck(i,j,k)
!WUHU FENG, more metal ions
                   bbuffer(bpter(ib,k)+8) = capmmr_blck(i,j,k)
                   bbuffer(bpter(ib,k)+9) = kpmmr_blck(i,j,k)
                   bbuffer(bpter(ib,k)+10) = sipmmr_blck(i,j,k)
                   !------------------Jianfei Wu-------------------------
                   if (sIndxOp > 0) then
!                      bbuffer(bpter(ib,k)+8) = opmmr_blck(i,j,k)
!WUHU FENG, more metal ions
                       bbuffer(bpter(ib,k)+11) = opmmr_blck(i,j,k)
                   endif
                   !-------------------------------------------------------------

                end do
             end do
          end do

          call t_barrierf('sync_ionos_blk_to_chk', grid%commxy)
          call t_startf ('ionos_block_to_chunk')
          call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
          call t_stopf  ('ionos_block_to_chunk')

          chnk_loop2: do lchnk = begchunk,endchunk
             ncol = phys_state(lchnk)%ncol

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

             call pbuf_get_field(pbuf_chnk, index_ui, ui_phys)
             call pbuf_get_field(pbuf_chnk, index_vi, vi_phys)
             call pbuf_get_field(pbuf_chnk, index_wi, wi_phys)
             if (sIndxOp > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys, &
                     start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
             endif

             call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)
             q => phys_state(lchnk)%q

             do i=1,ncol

                do k=1,km
                   ui_phys(i,k) = cbuffer(cpter(i,k))
                   vi_phys(i,k) = cbuffer(cpter(i,k)+1)
                   wi_phys(i,k) = cbuffer(cpter(i,k)+2)
                   q(i,k,ixfep)=cbuffer(cpter(i,k)+3)
                   q(i,k,ixmgp)=cbuffer(cpter(i,k)+4)
                   q(i,k,ixnap)=cbuffer(cpter(i,k)+5)
                   q(i,k,ixnop)=cbuffer(cpter(i,k)+6)
                   q(i,k,ixo2p)=cbuffer(cpter(i,k)+7)
                   q(i,k,ixcap)=cbuffer(cpter(i,k)+8)
                   q(i,k,ixkp)=cbuffer(cpter(i,k)+9)
                   q(i,k,ixsip)=cbuffer(cpter(i,k)+10)
!Wuhu Feng, more metal ions
                   if (sIndxOp > 0) then
!Wuhu Feng, more metal ions
!                     mmrPOp_phys(i,k) = cbuffer(cpter(i,k)+8)
                      mmrPOp_phys(i,k) = cbuffer(cpter(i,k)+11)
                   endif
                end do ! k=1,km
             end do ! i=1,ncol

             if (ionos_edyn_active) then
                call outfld ( 'UI', ui_phys, pcols, lchnk )
                call outfld ( 'VI', vi_phys, pcols, lchnk )
                call outfld ( 'WI', wi_phys, pcols, lchnk )
                if (write_inithist()) then
                   call outfld ( 'UI&IC', ui_phys, pcols, lchnk )
                   call outfld ( 'VI&IC', vi_phys, pcols, lchnk )
                   call outfld ( 'WI&IC', wi_phys, pcols, lchnk )
                endif
             endif

          end do chnk_loop2

          deallocate(bpter)
          deallocate(bbuffer)
          deallocate(cbuffer)

       endif blcks2phys_local

       if (sIndxOp>0) then
          deallocate(opmmr_blck)
          nullify(opmmr_blck)
       endif
       if (sIndxN2p>0) then
          deallocate(n2pmmr_blck)
          nullify(n2pmmr_blck)
       endif

       deallocate( wuxy )
       deallocate( wvxy )
       deallocate( sigma_ped_blck )
       deallocate( sigma_hall_blck )
       deallocate( ti_blck )
       deallocate( te_blck )
       deallocate( zi_blck )
       deallocate( ui_blck )
       deallocate( vi_blck )
       deallocate( wi_blck )
       deallocate( omega_blck )
       deallocate( tn_blck )
       deallocate( fepmmr_blck )
       deallocate( mgpmmr_blck )
       deallocate( napmmr_blck )
       deallocate( nopmmr_blck )
       deallocate( o2pmmr_blck )
!WUHU FENG, more metal ions
         deallocate( capmmr_blck )
         deallocate( kpmmr_blck )
         deallocate( sipmmr_blck )


    endif ionos_cpl

  end subroutine ionosphere_run2

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
 subroutine ionosphere_init_restart(File)
    use pio, only: file_desc_t, pio_double, pio_def_var
    use cam_pio_utils, only: cam_pio_def_dim
    use dyn_grid,      only: get_horiz_grid_dim_d

    type(File_desc_t),  intent(inout) :: File

    integer :: ierr,hdim1,hdim2, dimids(3)

    call get_horiz_grid_dim_d(hdim1, hdim2)

    call cam_pio_def_dim(File, 'lon',  hdim1, dimids(1),  existOK=.true.)
    call cam_pio_def_dim(File, 'lat',  hdim2, dimids(2),  existOK=.true.)
    call cam_pio_def_dim(File, 'lev',  pver,  dimids(3),  existOK=.true.)

    if (ionos_xport_active) then
       ierr = PIO_Def_Var(File, 'Optm1', pio_double, dimids, Optm1_vdesc)
       ierr = PIO_Def_Var(File, 'NOptm1', pio_double, dimids, NOptm1_vdesc)
       ierr = PIO_Def_Var(File, 'O2ptm1', pio_double, dimids, O2ptm1_vdesc)
       !-------------------Jianfei Wu------------------------------------
       ierr = PIO_Def_Var(File, 'Feptm1', pio_double, dimids, Feptm1_vdesc)
       ierr = PIO_Def_Var(File, 'Mgptm1', pio_double, dimids, Mgptm1_vdesc)
       ierr = PIO_Def_Var(File, 'Naptm1', pio_double, dimids, Naptm1_vdesc)
!Wuhu Feng, more metal ions
        ierr = PIO_Def_Var(File, 'Captm1', pio_double, dimids, Captm1_vdesc)
         ierr = PIO_Def_Var(File, 'Kptm1', pio_double, dimids, Kptm1_vdesc)
         ierr = PIO_Def_Var(File, 'Siptm1', pio_double, dimids, Siptm1_vdesc)

    endif
  end subroutine ionosphere_init_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_write_restart(File)
    use pio, only: io_desc_t, file_desc_t, pio_write_darray, pio_initdecomp, pio_double
    use cam_pio_utils, only: pio_subsystem
    use dyn_grid,      only: get_horiz_grid_dim_d

    type(File_desc_t), intent(inout) :: File

    type(io_desc_t) :: iodesc3d
    integer :: hdim1, hdim2
    integer, pointer :: ldof(:)
    integer :: ierr

    if (ionos_xport_active) then
       call get_horiz_grid_dim_d(hdim1, hdim2)
       ldof => get_restart_decomp(hdim1, hdim2, pver)
       call pio_initdecomp(pio_subsystem, pio_double, (/hdim1, hdim2, pver/), ldof, iodesc3d)
       deallocate(ldof)

       call pio_write_darray(File, Optm1_vdesc, iodesc3d, opmmrtm1_blck, ierr)
       call pio_write_darray(File, NOptm1_vdesc, iodesc3d, nopmmrtm1_blck, ierr)
       call pio_write_darray(File, O2ptm1_vdesc, iodesc3d, o2pmmrtm1_blck, ierr)
       !-----------------------Jianfei Wu------------------------------------
       call pio_write_darray(File, Feptm1_vdesc, iodesc3d, fepmmrtm1_blck, ierr)
       call pio_write_darray(File, Mgptm1_vdesc, iodesc3d, mgpmmrtm1_blck, ierr)
       call pio_write_darray(File, Naptm1_vdesc, iodesc3d, napmmrtm1_blck, ierr)
!Wuhu Feng, more metal ions
       call pio_write_darray(File, Captm1_vdesc, iodesc3d, capmmrtm1_blck, ierr)
       call pio_write_darray(File, Kptm1_vdesc, iodesc3d, kpmmrtm1_blck, ierr)
       call pio_write_darray(File, Siptm1_vdesc, iodesc3d, sipmmrtm1_blck, ierr)
    endif

  end subroutine ionosphere_write_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_read_restart(File)
    use pio, only: io_desc_t, file_desc_t, pio_inq_varid, pio_read_darray, pio_initdecomp, pio_double
    use cam_pio_utils, only: pio_subsystem
    use dyn_grid,      only: get_horiz_grid_dim_d

    type(file_desc_t), intent(inout) :: File

    integer :: ierr
    type(io_desc_t) :: iodesc3d
    integer :: hdim1, hdim2
    integer, pointer :: ldof(:)

    if (ionos_xport_active) then
       call ionosphere_alloc

       call get_horiz_grid_dim_d(hdim1, hdim2)
       ldof => get_restart_decomp(hdim1, hdim2, pver)
       call pio_initdecomp(pio_subsystem, pio_double, (/hdim1, hdim2, pver/), ldof, iodesc3d)
       deallocate(ldof)

       ierr = pio_inq_varid(File, 'Optm1', Optm1_vdesc)
       call pio_read_darray(File, Optm1_vdesc, iodesc3d, opmmrtm1_blck, ierr)
       ierr = pio_inq_varid(File, 'NOptm1', NOptm1_vdesc)
       call pio_read_darray(File, NOptm1_vdesc, iodesc3d, nopmmrtm1_blck, ierr)
       ierr = pio_inq_varid(File, 'O2ptm1', O2ptm1_vdesc)
       call pio_read_darray(File, O2ptm1_vdesc, iodesc3d, o2pmmrtm1_blck, ierr)
       !---------------------Jianfei Wu-------------------------------------
       ierr = pio_inq_varid(File, 'Feptm1', Feptm1_vdesc)
       call pio_read_darray(File, Feptm1_vdesc, iodesc3d, fepmmrtm1_blck, ierr)
       ierr = pio_inq_varid(File, 'Mgptm1', Mgptm1_vdesc)
       call pio_read_darray(File, Mgptm1_vdesc, iodesc3d, mgpmmrtm1_blck, ierr)
       ierr = pio_inq_varid(File, 'Naptm1', Naptm1_vdesc)
       call pio_read_darray(File, Naptm1_vdesc, iodesc3d, napmmrtm1_blck, ierr)
!Wuhu Feng, more metal ions
         ierr = pio_inq_varid(File, 'Captm1', Captm1_vdesc)
         call pio_read_darray(File, Captm1_vdesc, iodesc3d, capmmrtm1_blck, ierr)
         ierr = pio_inq_varid(File, 'Kptm1', Kptm1_vdesc)
         call pio_read_darray(File, Kptm1_vdesc, iodesc3d, kpmmrtm1_blck, ierr)
         ierr = pio_inq_varid(File, 'Siptm1', Siptm1_vdesc)
         call pio_read_darray(File, Siptm1_vdesc, iodesc3d, sipmmrtm1_blck, ierr)


    endif

  end subroutine ionosphere_read_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_final

#ifdef WACCMX_EDYN_ESMF
    use edyn_esmf, only: edyn_esmf_final

    call edyn_esmf_final()
#endif

    if (allocated(opmmrtm1_blck)) deallocate(opmmrtm1_blck)
    if (allocated(nopmmrtm1_blck)) deallocate(nopmmrtm1_blck)
    if (allocated(o2pmmrtm1_blck)) deallocate(o2pmmrtm1_blck)
    !--------------Jianfei Wu-----------------------------
    if (allocated(fepmmrtm1_blck)) deallocate(fepmmrtm1_blck)
    if (allocated(mgpmmrtm1_blck)) deallocate(mgpmmrtm1_blck)
    if (allocated(napmmrtm1_blck)) deallocate(napmmrtm1_blck)
!Wuhu Feng, more metal ions
    if (allocated(capmmrtm1_blck)) deallocate(capmmrtm1_blck)
    if (allocated(kpmmrtm1_blck)) deallocate(kpmmrtm1_blck)
    if (allocated(sipmmrtm1_blck)) deallocate(sipmmrtm1_blck)

  end subroutine ionosphere_final

!=========================================================================================
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_read_ic()

    use pio,          only: file_desc_t
    use ncdio_atm,    only: infld
    use cam_initfiles,      only: initial_file_get_id

    type(file_desc_t), pointer :: fh_ini    ! PIO filehandle

    type (t_fvdycore_grid), pointer :: grid
    integer :: ifirstxy,ilastxy,jfirstxy,jlastxy,km
    logical :: readvar

    if ( ionos_xport_active ) then
       call ionosphere_alloc()

       fh_ini   => initial_file_get_id()
       grid     => get_dyn_state_grid()
       ifirstxy =  grid%ifirstxy
       ilastxy  =  grid%ilastxy
       jfirstxy =  grid%jfirstxy
       jlastxy  =  grid%jlastxy
       km       =  grid%km

       ! try reading in OpTM1 from the IC file
       call infld('OpTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
            1, km, opmmrtm1_blck, readvar, gridname='fv_centers')

       if (.not.readvar) then
          ! if OpTM1 is not included in the IC file then try using O+
          call infld('Op', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1, km, opmmrtm1_blck, readvar, gridname='fv_centers')
       endif
       ! try reading in NOpTM1 from the IC file
       call infld('NOpTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
            1, km, nopmmrtm1_blck, readvar, gridname='fv_centers')

       if (.not.readvar) then
          ! if NOpTM1 is not included in the IC file then try using O+
          call infld('NOp', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1, km, nopmmrtm1_blck, readvar, gridname='fv_centers')
       endif
       ! try reading in O2pTM1 from the IC file
       call infld('O2pTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
            1, km, o2pmmrtm1_blck, readvar, gridname='fv_centers')

       if (.not.readvar) then
          ! if O2pTM1 is not included in the IC file then try using O+
          call infld('O2p', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1, km, o2pmmrtm1_blck, readvar, gridname='fv_centers')
       endif
       !-----------------------Jianfei Wu------------------------------------------------------
       ! try reading in FepTM1 from the IC file
       call infld('FepTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
            1, km, fepmmrtm1_blck, readvar, gridname='fv_centers')

       if (.not.readvar) then
          ! if FepTM1 is not included in the IC file then try using Fe+
          call infld('Fep', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1, km, fepmmrtm1_blck, readvar, gridname='fv_centers')
       endif
       ! try reading in MgpTM1 from the IC file
       call infld('MgpTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
            1, km, mgpmmrtm1_blck, readvar, gridname='fv_centers')

       if (.not.readvar) then
          ! if MgpTM1 is not included in the IC file then try using Fe+
          call infld('Mgp', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1, km, mgpmmrtm1_blck, readvar, gridname='fv_centers')
       endif
       ! try reading in NapTM1 from the IC file
       call infld('NapTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
            1, km, napmmrtm1_blck, readvar, gridname='fv_centers')

       if (.not.readvar) then
          ! if NapTM1 is not included in the IC file then try using Fe+
          call infld('Nap', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1, km, napmmrtm1_blck, readvar, gridname='fv_centers')
       endif
       !--------------------------------------------------------------------------------------
!Wuhu Feng, more metal ions
         call infld('CapTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
              1, km, capmmrtm1_blck, readvar, gridname='fv_centers')
 
         if (.not.readvar) then
            call infld('Cap', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
                 1, km, capmmrtm1_blck, readvar, gridname='fv_centers')
         endif
         call infld('KpTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
              1, km, kpmmrtm1_blck, readvar, gridname='fv_centers')
 
         if (.not.readvar) then
            call infld('Kp', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
                 1, km, kpmmrtm1_blck, readvar, gridname='fv_centers')
         endif
         call infld('SipTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
              1, km, sipmmrtm1_blck, readvar, gridname='fv_centers')
 
         if (.not.readvar) then
            call infld('Sip', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
                 1, km, sipmmrtm1_blck, readvar, gridname='fv_centers')
         endif
 
!--------------------------------------------------------------------------------------

    endif

  end subroutine ionosphere_read_ic

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_alloc

    type(T_FVDYCORE_GRID),pointer :: grid ! FV Dynamics grid
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy, km
    integer :: astat

    if (.not. allocated(opmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(opmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate opmmrtm1_blck')
       opmmrtm1_blck = 0._r8

    endif
    if (.not. allocated(nopmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(nopmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate nopmmrtm1_blck')
       nopmmrtm1_blck = 0._r8

    endif
    if (.not. allocated(o2pmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(o2pmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate o2pmmrtm1_blck')
       o2pmmrtm1_blck = 0._r8

    endif
    !----------------------Jianfei Wu--------------------------------------------
    if (.not. allocated(fepmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(fepmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate fepmmrtm1_blck')
       fepmmrtm1_blck = 0._r8

    endif
    if (.not. allocated(mgpmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(mgpmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate mgpmmrtm1_blck')
       mgpmmrtm1_blck = 0._r8

    endif
    if (.not. allocated(napmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(napmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate napmmrtm1_blck')
       napmmrtm1_blck = 0._r8

    endif
!WUHU FENG, more metal ions
    if (.not. allocated(capmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(capmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate napmmrtm1_blck')
       capmmrtm1_blck = 0._r8

    endif
    if (.not. allocated(kpmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(kpmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate napmmrtm1_blck')
       kpmmrtm1_blck = 0._r8

    endif
    if (.not. allocated(sipmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(sipmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate napmmrtm1_blck')
       sipmmrtm1_blck = 0._r8

    endif
    !-----------------------------------------------------------------------------------

  end subroutine ionosphere_alloc


  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
function get_restart_decomp(hdim1, hdim2, nlev) result(ldof)
   use dyn_grid, only: get_dyn_grid_parm

   ! Get the integer mapping of a variable in the dynamics decomp in memory.  
   ! The canonical ordering is as on the file. A 0 value indicates that the
   ! variable is not on the file (eg halo or boundary values)

   ! arguments
   integer, intent(in) :: hdim1, hdim2, nlev
   integer, pointer :: ldof(:)

   ! local variables
   integer :: i, k, j
   integer :: lcnt
   integer :: beglatxy, beglonxy, endlatxy, endlonxy
   !----------------------------------------------------------------------------

   beglonxy = get_dyn_grid_parm('beglonxy')
   endlonxy = get_dyn_grid_parm('endlonxy')
   beglatxy = get_dyn_grid_parm('beglatxy')
   endlatxy = get_dyn_grid_parm('endlatxy')

   lcnt = (endlatxy-beglatxy+1)*nlev*(endlonxy-beglonxy+1)
   allocate(ldof(lcnt))
   ldof(:) = 0	

   lcnt = 0
   do k = 1, nlev
      do j = beglatxy, endlatxy
         do i = beglonxy, endlonxy
            lcnt = lcnt + 1
            ldof(lcnt) = i + (j-(plat-hdim2+1))*hdim1+(k-1)*hdim1*hdim2
         end do
      end do
   end do

end function get_restart_decomp

!=========================================================================================


end module ionosphere_interface
