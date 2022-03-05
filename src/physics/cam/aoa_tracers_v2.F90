!---------------------------------------------------------------------------
! Age of air (AOA) diagnostic tracers -- Version 2
!---------------------------------------------------------------------------
module aoa_tracers_v2
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils, only: masterproc

  implicit none

  private
  public :: aoa_tracers_v2_readnl
  public :: aoa_tracers_v2_reg
  public :: aoa_tracers_v2_impl_cnst
  public :: aoa_tracers_v2_init_cnst
  public :: aoa_tracers_v2_init
  public :: aoa_tracers_v2_tend

  integer, parameter :: ncnst = 1 ! number of constituents implemented by this module

  ! AOA tracer cnstituents names
  character(len=16), parameter :: aoa_names(ncnst) = (/'AOA1mf'/)

  ! AOA tracer long names
  character(len=64), parameter :: lng_names(ncnst) = (/'Age-of_air tracer 1'/)

  ! AOA tracer cnstituents indices
  integer :: aoa_cnst_ndx(ncnst) = -1

  logical :: aoa_v2_active = .false.
  logical :: aoa_v2_readic = .true.

contains

  !---------------------------------------------------------------------------
  ! reads namelist options included in aoa_tracers_v2_opts group (in atm_in)
  !---------------------------------------------------------------------------
  subroutine aoa_tracers_v2_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils, only : mpicom, masterprocid, mpi_logical
    use cam_abortutils, only: endrun

    character(len=*), intent(in) :: nlfile
    integer :: unitn, ierr

    namelist /aoa_tracers_v2_opts/ aoa_v2_active, aoa_v2_readic

    ! read namelist only by masterpoc MPI tasks
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aoa_tracers_v2_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, aoa_tracers_v2_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun('aoa_tracers_v2_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! broadcast to all MPI tasks
    call mpi_bcast(aoa_v2_active, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(aoa_v2_readic, 1, mpi_logical, masterprocid, mpicom, ierr)

  end subroutine aoa_tracers_v2_readnl

  !---------------------------------------------------------------------------
  ! register the AOA tracers to be advected by CAM
  !---------------------------------------------------------------------------
  subroutine aoa_tracers_v2_reg()
    use physconst, only: cpair, mwdry
    use constituents, only: cnst_add

    integer :: i

    if (.not.aoa_v2_active) return

    do i=1,ncnst
       ! register the tracers to be advected by CAM
       call cnst_add( aoa_names(i), mwdry, cpair, 0._r8, aoa_cnst_ndx(i), &
                      readiv=aoa_v2_readic, longname=lng_names(i) )
    end do

  end subroutine aoa_tracers_v2_reg

  !---------------------------------------------------------------------------
  ! used by cam's read initial conditions code
  !---------------------------------------------------------------------------
  function aoa_tracers_v2_impl_cnst(name)
    character(len=*), intent(in) :: name       ! constituent name
    logical :: aoa_tracers_v2_impl_cnst  ! return value

    integer :: i

    aoa_tracers_v2_impl_cnst = .false.

    if (.not.aoa_v2_active) return

    do i = 1, ncnst
       if (trim(name) == trim(aoa_names(i))) then
          aoa_tracers_v2_impl_cnst = .true.
          return
       end if
    end do
  end function aoa_tracers_v2_impl_cnst

  !---------------------------------------------------------------------------
  ! used by cam's read initial conditions code
  !---------------------------------------------------------------------------
  subroutine aoa_tracers_v2_init_cnst(name, latvals, lonvals, mask, q)
    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize AOA tracers mixing ratio fields
    !  This subroutine is called at the beginning of an initial run ONLY
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in)  :: name       ! name of the tracer
    real(r8),         intent(in)  :: latvals(:) ! lat in degrees
    real(r8),         intent(in)  :: lonvals(:) ! lon in degrees
    logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     ! MMR (kg/kg) (gcol, plev)

    integer :: m
    !-----------------------------------------------------------------------

    do m = 1, ncnst
       if (name ==  aoa_names(m))  then
          q(:,:) = 0._r8
       endif
    end do

  end subroutine aoa_tracers_v2_init_cnst

  !---------------------------------------------------------------------------
  ! initializer -- add the AOA tracers to cam history
  !---------------------------------------------------------------------------
  subroutine aoa_tracers_v2_init()
    use cam_history, only: addfld, add_default

    integer :: i

    if (.not.aoa_v2_active) return

    ! add the AOA tracers to cam history
    do i=1,ncnst
       call addfld(aoa_names(i), (/ 'lev' /), 'A', 'kg/kg', lng_names(i))
       call add_default(aoa_names(i), 1, ' ')
    end do

  end subroutine aoa_tracers_v2_init

  !---------------------------------------------------------------------------
  ! return tendencies of the AOA tracers
  !---------------------------------------------------------------------------
  subroutine aoa_tracers_v2_tend(dt,state,ptend)
    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use ppgrid, only: pver
    use constituents, only: pcnst

    real(r8), intent(in) :: dt
    type(physics_state), intent(in) :: state
    type(physics_ptend), intent(out) :: ptend

    logical :: lq(pcnst)
    integer :: i, ncol, lchnk
    real(r8) :: dtinv
    real(r8) :: fixed_lbc_mmr

    if (.not.aoa_v2_active) return

    ncol = state%ncol
    lchnk = state%lchnk

    lq(:) = .false.

    do i=1,ncnst
       lq(aoa_cnst_ndx(i))=.true.
    end do

    call physics_ptend_init(ptend, state%psetcols, 'aoa_tracers_v2', lq=lq)

    dtinv = 1._r8/dt

    fixed_lbc_mmr = 1.e-6_r8 ! some tempary number for testing -- will be set
                             ! via some time-dependent algorithm

    do i=1,ncnst
       ptend%q(:ncol,:,aoa_cnst_ndx(i)) = 0._r8

       ! calculate the tendency of the MMR in the lowest layer corresponding
       ! to the fixed concentration at the lower boundary
       ptend%q(:ncol,pver,aoa_cnst_ndx(i)) = &
            (fixed_lbc_mmr-state%q(:ncol,pver,aoa_cnst_ndx(i)))*dtinv
    end do

  end subroutine aoa_tracers_v2_tend

end module aoa_tracers_v2
