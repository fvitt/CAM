module upper_bc

!---------------------------------------------------------------------------------
! Module to compute the upper boundary condition for temperature (dry static energy)
! and trace gases. Uses the MSIS model, and SNOE and TIME GCM data.
!
! original code by Stacy Walters
! adapted by B. A. Boville
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod,only: grav   => shr_const_g,     &   ! gravitational constant (m/s^2)
                          kboltz => shr_const_boltz, &   ! Boltzmann constant
                          pi => shr_const_pi,        &   ! pi
                          rEarth => shr_const_rearth     ! Earth radius
  use ppgrid,       only: pcols, pver, pverp
  use constituents, only: pcnst
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc
  use ref_pres,     only: do_molec_diff, ptop_ref
  use shr_kind_mod, only: cx=>SHR_KIND_CX
  use cam_abortutils,only: endrun
    use infnan,       only : nan, assignment(=)

  implicit none
  private
  save
!
! Public interfaces
!
  public :: ubc_readnl         ! read namelist options for UBCs
  public :: ubc_init           ! global initialization
  public :: ubc_timestep_init  ! time step initialization
  public :: ubc_get_vals       ! get ubc values for this step
  public :: ubc_get_flxs       ! get ub fluxes for this step
  public :: ubc_fixed_conc

  character(len=64) :: ubc_fields(pcnst) = 'NOTSET'
  character(len=cx) :: ubc_filepath = 'NONE'

  character(len=16) :: ubc_flds(pcnst) = 'NOTSET'
  real(r8) :: fixed_mmr(pcnst) = -huge(1._r8)
  real(r8) :: fixed_vmr(pcnst) = -huge(1._r8)

  integer :: num_fixed = 0
  character(len=2), parameter :: msis_spc(4) = &
       (/ 'H ','N ','O ','O2' /)
  character(len=2), parameter :: tgcm_spc(1) = &
       (/ 'H2' /)
  character(len=2), parameter :: snoe_spc(1) = &
       (/ 'NO' /)

  logical, public, protected :: msis_active =.false.
  logical :: tgcm_active =.false.
  logical :: snoe_active =.false.

! Namelist variables
  character(len=256) :: snoe_ubc_file = ' '
  real(r8)           :: t_pert_ubc  = 0._r8
  real(r8)           :: no_xfac_ubc = 1._r8

  integer :: h_ndx=-1

  character(len=256) :: tgcm_ubc_file = ' '
  integer            :: tgcm_ubc_cycle_yr = 0
  integer            :: tgcm_ubc_fixed_ymd = 0
  integer            :: tgcm_ubc_fixed_tod = 0
  character(len=32)  :: tgcm_ubc_data_type = 'CYCLICAL'

  logical :: apply_upper_bc = .false.
  logical :: reported = .false.

!================================================================================================
contains
!================================================================================================

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine ubc_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils, only : mpicom, masterprocid, mpi_character, mpi_integer, mpi_real8

    character(len=*), intent(in) :: nlfile
    integer :: unitn, ierr, n, ndx

    character(len=*), parameter :: prefix = 'ubc_readnl: '

    namelist /upper_bc_opts/ tgcm_ubc_file,tgcm_ubc_data_type,tgcm_ubc_cycle_yr,tgcm_ubc_fixed_ymd, &
                             tgcm_ubc_fixed_tod, snoe_ubc_file, no_xfac_ubc, t_pert_ubc
    namelist /upper_bc_opts/ ubc_fields, ubc_filepath

    if (masterproc) then
       ! read namelist
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'upper_bc_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, upper_bc_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(prefix//'upper_bc_opts: ERROR reading namelist')
          end if
       end if
       close(unitn)

       ! log the UBC options
       write(iulog,*) prefix//'tgcm_ubc_file = '//trim(tgcm_ubc_file)
       write(iulog,*) prefix//'tgcm_ubc_data_type = '//trim(tgcm_ubc_data_type)
       write(iulog,*) prefix//'tgcm_ubc_cycle_yr = ', tgcm_ubc_cycle_yr
       write(iulog,*) prefix//'tgcm_ubc_fixed_ymd = ', tgcm_ubc_fixed_ymd
       write(iulog,*) prefix//'tgcm_ubc_fixed_tod = ', tgcm_ubc_fixed_tod
       write(iulog,*) prefix//'snoe_ubc_file = '//trim(snoe_ubc_file)
       write(iulog,*) prefix//'t_pert_ubc = ', t_pert_ubc
       write(iulog,*) prefix//'no_xfac_ubc = ', no_xfac_ubc

       write(iulog,*) prefix//'ubc_filepath = '//trim(ubc_filepath)

       write(iulog,*) prefix//'ubc_fields : '
       n=1
       do while(ubc_fields(n)/='NOTSET')
          write(iulog,'(i4,a)') n,'  '//trim(ubc_fields(n))
          ndx = index(ubc_fields(n),'=')
          if (ndx>0) then
             ubc_flds(n) = ubc_fields(n)(:ndx-1)
          else
             ubc_flds(n) = ubc_fields(n)
          end if
          n=n+1
       end do
       num_fixed=n-1
    end if


    ! broadcast to all MPI tasks
    call mpi_bcast(num_fixed, 1, mpi_integer, masterprocid, mpicom, ierr)

    call mpi_bcast (tgcm_ubc_file,      len(tgcm_ubc_file), mpi_character,masterprocid, mpicom, ierr)
    call mpi_bcast (tgcm_ubc_data_type, len(tgcm_ubc_data_type),mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast (tgcm_ubc_cycle_yr,  1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast (tgcm_ubc_fixed_ymd, 1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast (tgcm_ubc_fixed_tod, 1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast (snoe_ubc_file, len(snoe_ubc_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast (t_pert_ubc,    1,  mpi_real8, masterprocid, mpicom, ierr)
    call mpi_bcast (no_xfac_ubc,   1,  mpi_real8, masterprocid, mpicom, ierr)

    call mpi_bcast(ubc_filepath, 1, mpi_character,masterprocid, mpicom, ierr)
    call mpi_bcast(ubc_fields, pcnst*len(ubc_fields(1)), mpi_character,masterprocid, mpicom, ierr)
    call mpi_bcast(ubc_flds, pcnst*len(ubc_flds(1)), mpi_character,masterprocid, mpicom, ierr)

  end subroutine ubc_readnl

!===============================================================================

  subroutine ubc_init()
!-----------------------------------------------------------------------
! Initialization of time independent fields for the upper boundary condition
! Calls initialization routine for MSIS, TGCM and SNOE
!-----------------------------------------------------------------------
    use mo_tgcm_ubc, only: tgcm_ubc_inti
    use mo_snoe,     only: snoe_inti
    use mo_msis_ubc, only: msis_ubc_inti
    use constituents,only: cnst_get_ind

!---------------------------Local workspace-----------------------------
    logical, parameter :: zonal_avg = .false.
    integer :: m, spc_ndx
    integer :: ndx, mmrndx, vmrndx

    real(r8) :: val
    character(len=32) :: fld, str, valstr
    character(len=*), parameter :: prefix = 'ubc_init: '


    !-----------------------------------------------------------------------

    if (masterproc)   write(iulog,*) 'FVDBG.ubc_init ptop_ref = ', ptop_ref

    apply_upper_bc = .true. !do_molec_diff ! ptop_ref<1._r8 ! Pa

    if (.not.apply_upper_bc) return

    call cnst_get_ind('H', h_ndx, abort=.false.)

    if (ubc_filepath == 'NONE') then

       do m = 1,num_fixed
          if (ubc_fields(m)=='T') then
             msis_active = .true.
          else

             ndx = index(ubc_fields(m),'=')
             if (ndx>0) then
                fld = trim(ubc_fields(m)(:ndx-1))
                valstr = trim(ubc_fields(m)(ndx+1:))
                call cnst_get_ind(fld, spc_ndx, abort=.false.)
                if (spc_ndx>0) then
                   mmrndx = index(valstr,'mmr')
                   vmrndx = index(valstr,'vmr')
                   if (masterproc) write(iulog,*) 'FVDBG.ubc_init field, m,mmrndx,vmrndx: '//trim(ubc_fields(m)),&
                                                     m,mmrndx,vmrndx
                   if (mmrndx>0) then
                      str = valstr(:mmrndx-1)
                      read(str,*) val
                      fixed_mmr(spc_ndx) = val
                      if (masterproc) write(iulog,*) 'FVDBG.ubc_init field: '//trim(ubc_fields(m))//' MMR',val
                   else if (vmrndx>0) then
                      str = valstr(:vmrndx-1)
                      if (masterproc) write(iulog,*) 'FVDBG.ubc_init valstr: >|'//trim(valstr)//'|< vmrndx=',vmrndx
                      if (masterproc) write(iulog,*) 'FVDBG.ubc_init    str: >|'//trim(valstr)//'|< '
                      read(str,*) val
                      fixed_vmr(spc_ndx) = val
                      if (masterproc) write(iulog,*) 'FVDBG.ubc_init field: '//trim(ubc_fields(m))//' VMR',val
                   else
                      call endrun(prefix//'ubc value must include mmr or vmr units')
                   end if

                end if
             else
                fld = trim(ubc_fields(m))
                call cnst_get_ind(fld, spc_ndx, abort=.false.)
             endif

             if (spc_ndx>0) then
                if (.not.msis_active) then
                   msis_active = (any(msis_spc == ubc_flds(m)))
                endif
                if (.not.snoe_active) then
                   snoe_active = (any(snoe_spc == ubc_flds(m)))
                endif
                if (.not.tgcm_active) then
                   tgcm_active = (any(tgcm_spc == ubc_flds(m)))
                endif
             endif
          end if
       end do
       if (masterproc) then
          write(iulog,*) 'FVDBG.ubc_init msis_active = ',msis_active
          write(iulog,*) 'FVDBG.ubc_init snoe_active = ',snoe_active
          write(iulog,*) 'FVDBG.ubc_init tgcm_active = ',tgcm_active
       end if

       if (tgcm_active) then
          !-----------------------------------------------------------------------
          !       ... initialize the tgcm upper boundary module
          !-----------------------------------------------------------------------
          call tgcm_ubc_inti( tgcm_ubc_file, tgcm_ubc_data_type, tgcm_ubc_cycle_yr, &
               tgcm_ubc_fixed_ymd, tgcm_ubc_fixed_tod)
          if (masterproc) write(iulog,*) 'ubc_init: after tgcm_ubc_inti'
       endif

       if (snoe_active) then
          !-----------------------------------------------------------------------
          !       ... initialize the snoe module
          !-----------------------------------------------------------------------
          call snoe_inti(snoe_ubc_file)
          if (masterproc) write(iulog,*) 'ubc_init: after snoe_inti'
       endif

       if (msis_active) then
          !-----------------------------------------------------------------------
          !       ... initialize the msis module
          !-----------------------------------------------------------------------
          call msis_ubc_inti( zonal_avg )
          if (masterproc) write(iulog,*) 'ubc_init: after msis_ubc_inti'
       endif

    else

    endif

  end subroutine ubc_init

!===============================================================================
!===============================================================================

  logical function ubc_fixed_conc(name)

    character(len=*), intent(in) :: name

    integer :: m

    ubc_fixed_conc = .false.

    do m = 1,num_fixed
       if ( trim(ubc_flds(m)) == trim(name) ) then
          ubc_fixed_conc = .true.
          return
       endif
    end do

  end function ubc_fixed_conc

!===============================================================================

  subroutine ubc_timestep_init(pbuf2d, state)
!-----------------------------------------------------------------------
! timestep dependent setting
!-----------------------------------------------------------------------

    use solar_parms_data, only: kp=>solar_parms_kp, ap=>solar_parms_ap, f107=>solar_parms_f107
    use solar_parms_data, only: f107a=>solar_parms_f107a, f107p=>solar_parms_f107p
    use mo_msis_ubc,      only: msis_timestep_init
    use mo_tgcm_ubc,      only: tgcm_timestep_init
    use mo_snoe,          only: snoe_timestep_init
    use physics_types,    only: physics_state
    use ppgrid,           only: begchunk, endchunk
    use physics_buffer,   only: physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (.not.apply_upper_bc) return

    if (msis_active) &
         call msis_timestep_init( ap, f107p, f107a )
    if (tgcm_active) &
         call tgcm_timestep_init( pbuf2d, state )
    if (snoe_active) &
         call snoe_timestep_init( kp, f107 )

  end subroutine ubc_timestep_init

!===============================================================================

  subroutine ubc_get_vals (lchnk, ncol, pint, zi, msis_temp, ubc_mmr)

!-----------------------------------------------------------------------
! interface routine for vertical diffusion and pbl scheme
!-----------------------------------------------------------------------
    use mo_msis_ubc,      only: get_msis_ubc
    use mo_snoe,          only: set_no_ubc, ndx_no
    use mo_tgcm_ubc,      only: set_tgcm_ubc
    use cam_abortutils,   only: endrun
    use physconst,        only: rairv, mbarv
    use constituents,     only: cnst_mw, cnst_fixed_ubc, cnst_name  ! Needed for ubc_flux

!------------------------------Arguments--------------------------------
    integer,  intent(in)  :: lchnk                 ! chunk identifier
    integer,  intent(in)  :: ncol                  ! number of atmospheric columns
    real(r8), intent(in)  :: pint(pcols,pverp)     ! interface pressures
    real(r8), intent(in)  :: zi(pcols,pverp)       ! interface geoptl height above sfc

    real(r8), intent(out) :: msis_temp(pcols)      ! upper bndy temperature (K)
    real(r8), intent(out) :: ubc_mmr(pcols,pcnst)  ! upper bndy mixing ratios (kg/kg)

!---------------------------Local storage-------------------------------
    integer :: m                                   ! constituent index
    real(r8) :: rho_top(pcols)                     ! density at top interface
    real(r8) :: z_top(pcols)                       ! height of top interface (km)

    real(r8), parameter :: m2km = 1.e-3_r8         ! meter to km

    !-----------------------------------------------------------------------

    character(len=16) :: source(pcnst) = 'NONE'

    real(r8), parameter :: NOTSET = -huge(1._r8)

    ubc_mmr(:,:) = nan
    msis_temp(:) = nan
    ubc_mmr(:,:) = NOTSET ! nan

    if (.not. apply_upper_bc) return

    if (ubc_filepath == 'NONE') then

       if (msis_active) then
          call get_msis_ubc( lchnk, ncol, msis_temp, ubc_mmr )
          if( t_pert_ubc /= 0._r8 ) then
             msis_temp(:ncol) = msis_temp(:ncol) + t_pert_ubc
             if( any( msis_temp(:ncol) < 0._r8 ) ) then
                write(iulog,*) 'ubc_get_vals: msis temp < 0 after applying offset = ',t_pert_ubc
                call endrun
             end if
          end if
       end if
       do m = 1,pcnst
          if (source(m)=='NONE' .and. any(ubc_mmr(:ncol,m)/=NOTSET) )  &
               source(m)='MSIS'
       end do

       if (snoe_active) then

          rho_top(:ncol) = pint(:ncol,1) / (rairv(:ncol,1,lchnk)*msis_temp(:ncol))
          z_top(:ncol)   = m2km * zi(:ncol,1)

          call set_no_ubc  ( lchnk, ncol, z_top, ubc_mmr, rho_top )
          if( ndx_no > 0 .and. no_xfac_ubc /= 1._r8 ) then
             ubc_mmr(:ncol,ndx_no) = no_xfac_ubc * ubc_mmr(:ncol,ndx_no)
          end if

       endif

       do m = 1,pcnst
          if (source(m)=='NONE' .and. any(ubc_mmr(:ncol,m)/=NOTSET) )  &
               source(m)='SNOE'
       end do


       if (tgcm_active) then
          call set_tgcm_ubc( lchnk, ncol, ubc_mmr )
       endif

       do m = 1,pcnst
          if (source(m)=='NONE' .and. any(ubc_mmr(:ncol,m)/=NOTSET) )  &
               source(m)='TGCM'
       end do

    else

    endif

    ! fixed values
    do m = 1,pcnst
       if (fixed_mmr(m) > 0._r8) then
          ubc_mmr(:ncol, m) = fixed_mmr(m)
          if (masterproc.and..not.reported) write(iulog,*) 'FVDBG cnst_name(m),fixed_mmr(m):',cnst_name(m),fixed_mmr(m)
       elseif (fixed_vmr(m) > 0._r8) then
          ubc_mmr(:ncol, m) = cnst_mw(m)*fixed_vmr(m)/mbarv(:ncol,1,lchnk)
          if (masterproc.and..not.reported) write(iulog,*) 'FVDBG cnst_name(m),fixed_vmr(m):',cnst_name(m),fixed_vmr(m)
       end if
    end do

    do m = 1,pcnst
       if (source(m)=='NONE' .and. any(ubc_mmr(:ncol,m)/=NOTSET) )  &
            source(m)='CONST_MMR'
    end do

    ! Zero out constituent ubc's that are not used.
    do m = 1, pcnst
       if (.not. cnst_fixed_ubc(m)) then
          ubc_mmr(:,m) = 0.0_r8
       end if
    end do


    if (masterproc.and..not.reported) then
       do m = 1,pcnst
          if ( cnst_fixed_ubc(m) ) &
             write(iulog,*) 'FVDBG.ubc_get_vals  cnst_name,source,max mmr: '//cnst_name(m)//source(m), maxval(ubc_mmr(:ncol,m)), cnst_fixed_ubc(m)
       end do
       reported=.true.
       write(iulog,*) 'FVDBG.ubc_get_vals  msis_temp: ',msis_temp(1) !maxval( msis_temp(:ncol))
    endif

  end subroutine ubc_get_vals

!===============================================================================

  subroutine ubc_get_flxs (lchnk, ncol, pint, zi, t, q, phis, ubc_flux)

    use physconst, only: avogad, rairv, rga
    use constituents, only: cnst_mw
!------------------------------Arguments--------------------------------
    integer,  intent(in)  :: lchnk                 ! chunk identifier
    integer,  intent(in)  :: ncol                  ! number of atmospheric columns
    real(r8), intent(in)  :: pint(pcols,pverp)     ! interface pressures
    real(r8), intent(in)  :: zi(pcols,pverp)       ! interface geoptl height above sfc
    real(r8), intent(in)  :: t(pcols,pver)         ! midpoint temperature
    real(r8), intent(in),target :: q(pcols,pver,pcnst)   ! contituent mixing ratios (kg/kg)
    real(r8), intent(in)  :: phis(pcols)           ! Surface geopotential (m2/s2)

    real(r8), intent(out) :: ubc_flux(pcols,pcnst) ! upper bndy flux (kg/s/m^2)

!---------------------------Local storage-------------------------------
    integer :: iCol                                ! column loop counter

    real(r8), parameter :: hfluxlimitfac = 0.72_r8 ! Hydrogen upper boundary flux limiting factor

    real(r8) :: nmbartop                           ! Top level density (rho)
    real(r8) :: zkt                                ! Factor for H Jean's escape flux calculation

    real(r8), pointer :: qh_top(:)         ! Top level hydrogen mixing ratio (kg/kg)

    ubc_flux(:,:) = nan

    if (.not. apply_upper_bc) return
!!$    if (.not. ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) ) return

    qh_top => q(:,1,h_ndx)

    do iCol = 1, ncol
       !--------------------------------------------------
       ! Get total density (rho) at top level
       !--------------------------------------------------
       nmbartop = 0.5_r8 * (pint(iCol,1) + pint(iCol,2)) / ( rairv(iCol,1,lchnk) * t(iCol,1) )

       !---------------------------------------------------------------------
       ! Calculate factor for Jean's escape flux once here, used twice below
       !---------------------------------------------------------------------
       zkt = (rEarth + ( 0.5_r8 * ( zi(iCol,1) + zi(iCol,2) ) + rga * phis(iCol) ) ) * &
            cnst_mw(h_ndx) / avogad * grav / ( kboltz * t(iCol,1) )

       ubc_flux(iCol,h_ndx) = hfluxlimitfac * SQRT(kboltz/(2.0_r8 * pi * cnst_mw(h_ndx) / avogad)) * &
            qh_top(iCol) * nmbartop * &
            SQRT(t(iCol,1)) * (1._r8 + zkt) * EXP(-zkt)

       ubc_flux(iCol,h_ndx) = ubc_flux(iCol,h_ndx) * &
            (2.03E-13_r8 * qh_top(iCol) * nmbartop / (cnst_mw(h_ndx) / avogad) * t(iCol,1))

    enddo

  end subroutine ubc_get_flxs

end module upper_bc
