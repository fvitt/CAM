module met_input_stream
  use shr_kind_mod, only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use spmd_utils, only : mpicom, masterproc, iam
  use cam_logfile       , only : iulog
  use cam_abortutils    , only : endrun
  use dshr_strdata_mod  , only : shr_strdata_type
  use ppgrid,       only: pcols, pver

  use ESMF, only: ESMF_Finalize, ESMF_LogFoundError
  use ESMF, only: ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT

  implicit none

  logical :: met_input_stream_initialized = .false.

  type(shr_strdata_type) :: strmdata
  character(len=CL) :: input_mesh_filename
  character(len=CL) :: input_data_filename
  integer :: input_data_year_first ! first year in stream to use
  integer :: input_data_year_last  ! last year in stream to use
  integer :: input_data_year_align ! align stream_year_firstndep with

  character(len=CS) :: stream_varlist(1)

contains

  subroutine met_input_stream_datainit(model_mesh, model_clock)

    use dshr_strdata_mod, only: shr_strdata_init_from_inline
    use ESMF, only : ESMF_Clock, ESMF_Mesh

    ! input/output variables
    type(ESMF_CLock), intent(in)  :: model_clock
    type(ESMF_Mesh) , intent(in)  :: model_mesh

    integer :: rc


    stream_varlist = (/ 'T' /)

    input_mesh_filename = '/terminator-data1/home/fvitt/cdeps_test_inputs/MERRA2_orig_res_ESMF_mesh.nc'

    input_data_filename = '/terminator-data1/home/fvitt/cdeps_test_inputs/MERRA2_orig_res_20240101.nc'
    input_data_year_first = 2024
    input_data_year_last = 2024
    input_data_year_align = 1

    !return

    ! Initialize the cdeps data type sdat_ndep
    call shr_strdata_init_from_inline(strmdata,                    &
         my_task             = iam,                                 &
         logunit             = iulog,                               &
         compname            = 'ATM',                               &
         model_clock         = model_clock,                         &
         model_mesh          = model_mesh,                          &
         stream_meshfile     = trim(input_mesh_filename),     &
         stream_filenames    = (/trim(input_data_filename)/), &
         stream_yearFirst    = input_data_year_first,              &
         stream_yearLast     = input_data_year_last,               &
         stream_yearAlign    = input_data_year_align,              &
         stream_fldlistFile  = stream_varlist,                 &
         stream_fldListModel = stream_varlist,                 &
         stream_lev_dimname  = 'lev',                              &
!!$         stream_mapalgo      = 'bilinear',                          &
         stream_mapalgo      = 'nn',                          &
         stream_offset       = 0,                                   &
         stream_taxmode      = 'cycle',                             &
         stream_dtlimit      = 1.0e3_r8,                           &
         stream_tintalgo     = 'linear',                            &
         stream_name         = 'Met Inputs Stream',         &
         rc                  = rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    if (rc/=ESMF_SUCCESS) then
       call endrun('met_input_stream_init ERROR')
    end if

    met_input_stream_initialized = .true.

  end subroutine met_input_stream_datainit


  subroutine met_input_stream_histinit

    use cam_history, only: addfld, horiz_only

    call addfld ('T_input', (/ 'lev' /), 'A','K','T input test')
    call addfld ('T_input_srf', horiz_only, 'A','K','Surf T input test')

  end subroutine met_input_stream_histinit


  subroutine met_input_stream_adv(phys_state,pbuf2d)
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance

    use ppgrid           , only : begchunk, endchunk
    use time_manager     , only : get_curr_date
    use phys_grid        , only : get_ncols_p
    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    use interpolate_data, only: lininterp

    use cam_history_support, only: fillvalue
    use cam_history, only: outfld

    type(physics_state), intent(in) :: phys_state(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local variables
    integer :: i,c,g,k,  rc
    integer :: year    ! year (0, ...) for nstep+1
    integer :: mon     ! month (1, ..., 12) for nstep+1
    integer :: day     ! day of month (1, ..., 31) for nstep+1
    integer :: sec     ! seconds into current date for nstep+1
    integer :: mcdate  ! Current model date (yyyymmdd)
    real(r8), pointer :: dataptr(:,:)
    integer, parameter :: nlevs_met = 72

    real(r8) :: tmp_out(pcols,pver)
    real(r8) :: tmp_in(pcols,nlevs_met)
    !-----------------------------------------------------------------------


    ! Advance sdat stream
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day
    call shr_strdata_advance(strmdata, ymd=mcdate, tod=sec, logunit=iulog, istr='ndepdyn', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    if (rc/=ESMF_SUCCESS) then
       call endrun('met_input_stream_adv ERROR')
    end if


    ! Get pointer for stream data that is time and spatially interpolated to model time and grid
    call dshr_fldbun_getFldPtr(strmdata%pstrm(1)%fldbun_model, stream_varlist(1), fldptr2=dataptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    if (rc/=ESMF_SUCCESS) then
       call endrun('met_input_stream_adv ERROR')
    end if

    g = 1
    do c = begchunk,endchunk
       tmp_in = fillvalue
       do i = 1,get_ncols_p(c)
          do k = 1,nlevs_met
             tmp_in(i,k) = dataptr(k,g)
          end do
          g = g + 1
          call outfld('T_input_srf', tmp_in(:,nlevs_met), pcols, c)
       end do
    end do

  end subroutine met_input_stream_adv

end module met_input_stream
