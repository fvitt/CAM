!-------------------------------------------------------------------------------
! solar energetic proton data -- mean energy (eo), energy flux (fe), and number flux (fn)
!-------------------------------------------------------------------------------
module spe_data

  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use input_data_utils, only : time_coordinate
  use infnan,           only : nan, assignment(=)
  use spmd_utils,       only: masterproc

  implicit none

  private
  save

 ! public interface

  public :: spe_init
  public :: spe_advance

  logical, public :: spe_on = .false.

 ! time-interpolated quantities

  real(r8), public, protected :: spe_eo = 0._r8
  real(r8), public, protected :: spe_fe = 0._r8
  real(r8), public, protected :: spe_fn = 0._r8

 ! private data

  real(r8), allocatable :: eo_in(:)
  real(r8), allocatable :: fe_in(:)
  real(r8), allocatable :: fn_in(:)

  type(time_coordinate) :: time_coord

contains

  subroutine spe_init(filepath, fixed, fixed_ymd, fixed_tod)
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

    spe_eo = 0._r8
    spe_fe = 0._r8
    spe_fn = 0._r8

 !   if (masterproc) write(6,"('spe_init:intitial fixed_ymd, fixed_tod = ',2i10)")fixed_ymd, fixed_tod

    spe_on = (trim(filepath).ne.'NONE' .and. len_trim(filepath)>0)

    if (.not.spe_on) return

    !-----------------------------------------------------------------------
    !	... readin the solar parms dataset
    !-----------------------------------------------------------------------

    call getfil(filepath,  locfn, 0)
    call cam_pio_openfile ( ncid, locfn, PIO_NOWRITE)

    call time_coord%initialize( filepath, fixed=fixed, fixed_ymd=fixed_ymd, &
                                fixed_tod=fixed_tod, force_time_interp=.true. )

 !   if (masterproc) write(6,"('spe_init:after fixed_ymd, fixed_tod = ',2i10)")fixed_ymd, fixed_tod

    !---------------------------------------------------------------
    !	... allocate and read spe parms
    !---------------------------------------------------------------
    allocate( eo_in(time_coord%ntimes), fe_in(time_coord%ntimes), &
              fn_in(time_coord%ntimes), stat=astat )
    if( astat /= 0 ) then
       call alloc_err( astat, 'spe_init', 'eo_in ... fe_in ', &
                       time_coord%ntimes )
    end if

    ierr = pio_inq_varid( ncid, 'eo', varid )
    ierr = pio_get_var( ncid, varid, eo_in )
    ierr = pio_inq_varid( ncid, 'fe', varid )
    ierr = pio_get_var( ncid, varid, fe_in )
    ierr = pio_inq_varid( ncid, 'fn', varid )
    ierr = pio_get_var( ncid, varid, fn_in )

!   write(6,"('spe_init: eo_in = ',/,(10f10.3))") eo_in(:)

    call pio_closefile( ncid )

end subroutine spe_init

subroutine spe_advance
  !---------------------------------------------------------------
  ! time interpolate space wx indices
  !---------------------------------------------------------------

  if (spe_on) then
     call time_coord%advance()
    !  time interpolate
     spe_eo = time_coord%wghts(1)*eo_in(time_coord%indxs(1)) &
                      + time_coord%wghts(2)*eo_in(time_coord%indxs(2))
     spe_fe = time_coord%wghts(1)*fe_in(time_coord%indxs(1)) &
                      + time_coord%wghts(2)*fe_in(time_coord%indxs(2))
     spe_fn = time_coord%wghts(1)*fn_in(time_coord%indxs(1)) &
                      + time_coord%wghts(2)*fn_in(time_coord%indxs(2))

     if (masterproc) &
       write(6,"('spe_data:spe_eo,spe_fe,spe_fn = ',3g10.3)")spe_eo,spe_fe,spe_fn
  endif

end subroutine spe_advance

end module spe_data
