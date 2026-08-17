!-------------------------------------------------------------------------------
! solar energetic proton data -- mean energy (eo), energy flux (fe), and number flux (fn)
!-------------------------------------------------------------------------------
module solar_proton_data

  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use input_data_utils, only : time_coordinate
  use infnan,           only : nan, assignment(=)
  use spmd_utils,       only: masterproc

  implicit none

  private

  ! public interface

  public :: solar_proton_init
  public :: solar_proton_advance

  logical, public :: solar_proton_on = .false.

  ! time-interpolated quantities

  real(r8), public, protected :: solar_proton_eo = 0._r8 ! mean energy
  real(r8), public, protected :: solar_proton_fe = 0._r8 ! energy flux
  real(r8), public, protected :: solar_proton_fn = 0._r8 ! number flux

  ! private data

  real(r8), allocatable :: eo_in(:)
  real(r8), allocatable :: fe_in(:)
  real(r8), allocatable :: fn_in(:)

  type(time_coordinate) :: time_coord

contains

  !=============================================================================
  subroutine solar_proton_init(filepath, fixed, fixed_ymd, fixed_tod)
    !---------------------------------------------------------------
    !	... initialize solar parmaters
    !---------------------------------------------------------------

    use ioFileMod
    use error_messages, only: alloc_err
    use cam_pio_utils,  only: cam_pio_openfile
    use pio,            only: file_desc_t, var_desc_t, pio_get_var, &
         pio_inq_varid, pio_closefile, pio_nowrite

    !---------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------
    character(len=*), intent(in) :: filepath
    logical, intent(in) :: fixed
    integer, intent(in) :: fixed_ymd
    integer, intent(in) :: fixed_tod

    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    type(file_desc_t)  :: ncid
    type(var_desc_t)  :: varid
    integer  :: astat
    character(len=shr_kind_cl) :: locfn
    integer :: ierr

    solar_proton_eo = 0._r8
    solar_proton_fe = 0._r8
    solar_proton_fn = 0._r8

    solar_proton_on = (trim(filepath).ne.'NONE' .and. len_trim(filepath)>0)

    if (.not.solar_proton_on) return

    !-----------------------------------------------------------------------
    !	... readin the solar parms dataset
    !-----------------------------------------------------------------------

    call getfil(filepath,  locfn, 0)
    call cam_pio_openfile ( ncid, locfn, PIO_NOWRITE)

    call time_coord%initialize( filepath, fixed=fixed, fixed_ymd=fixed_ymd, &
         fixed_tod=fixed_tod, force_time_interp=.true. )

    !---------------------------------------------------------------
    !	... allocate and read spe parms
    !---------------------------------------------------------------
    allocate( eo_in(time_coord%ntimes), fe_in(time_coord%ntimes), &
              fn_in(time_coord%ntimes), stat=astat )
    if( astat /= 0 ) then
       call alloc_err( astat, 'solar_proton_init', 'eo_in ... fe_in ', &
            time_coord%ntimes )
    end if

    ierr = pio_inq_varid( ncid, 'eo', varid )
    ierr = pio_get_var( ncid, varid, eo_in )
    ierr = pio_inq_varid( ncid, 'fe', varid )
    ierr = pio_get_var( ncid, varid, fe_in )
    ierr = pio_inq_varid( ncid, 'fn', varid )
    ierr = pio_get_var( ncid, varid, fn_in )

    call pio_closefile( ncid )

  end subroutine solar_proton_init

  !=============================================================================
  subroutine solar_proton_advance
    !---------------------------------------------------------------
    ! time interpolate space wx indices
    !---------------------------------------------------------------

    if (solar_proton_on) then
       call time_coord%advance()
       !  time interpolate
       solar_proton_eo = time_coord%wghts(1)*eo_in(time_coord%indxs(1)) &
            + time_coord%wghts(2)*eo_in(time_coord%indxs(2))
       solar_proton_fe = time_coord%wghts(1)*fe_in(time_coord%indxs(1)) &
            + time_coord%wghts(2)*fe_in(time_coord%indxs(2))
       solar_proton_fn = time_coord%wghts(1)*fn_in(time_coord%indxs(1)) &
            + time_coord%wghts(2)*fn_in(time_coord%indxs(2))

       if (masterproc) then
          write(6,"('solar_proton_data:solar_proton_eo,solar_proton_fe,solar_proton_fn = ',3g10.3)") &
               solar_proton_eo,solar_proton_fe,solar_proton_fn
       end if
    endif

  end subroutine solar_proton_advance

end module solar_proton_data
