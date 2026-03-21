module meped_input_stream

  use shr_kind_mod, only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use dshr_strdata_mod, only : shr_strdata_type
  use spmd_utils, only : iam
  use cam_logfile, only : iulog

  use ESMF, only : ESMF_Clock, ESMF_Mesh
  use ESMF, only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
  use ESMF, only : ESMF_Finalize, ESMF_LogFoundError

  use ppgrid, only : pcols, begchunk,endchunk
  use cam_history, only: addfld, add_default, horiz_only, outfld

  implicit none

  type(shr_strdata_type) :: stream_north ! input data stream for northern hemisphere
  type(shr_strdata_type) :: stream_south ! input data stream for southern hemisphere

  character(len=*), parameter :: meped_filepath = &
       '/data/terminator-data1/home/fvitt/camdev/ganglu_epp/inputs/med_ring_exponential_oct03.ncf5.time2.nc'
  character(len=*), parameter :: meped_north_mesh = &
       '/data/terminator-data1/home/fvitt/camdev/ganglu_epp/inputs/ESMFmesh_MEDEP_north_nomask_c260320.cdf5.nc'
  character(len=*), parameter :: meped_south_mesh = &
       '/data/terminator-data1/home/fvitt/camdev/ganglu_epp/inputs/ESMFmesh_MEDEP_south_nomask_c260320.cdf5.nc'

  character(len=*), parameter :: meped_north_varlist(*) = (/'e_flux_nh','e_ekev_nh','p_flux_nh','p_ekev_nh'/)
  character(len=*), parameter :: meped_south_varlist(*) = (/'e_flux_sh','e_ekev_sh','p_flux_sh','p_ekev_sh'/)

  real(r8),protected, pointer :: meped_e_flux(:,:) => null()
  real(r8),protected, pointer :: meped_p_flux(:,:) => null()
  real(r8),protected, pointer :: meped_e_ekev(:,:) => null()
  real(r8),protected, pointer :: meped_p_ekev(:,:) => null()

  integer, parameter :: stream_meped_year_first = 2003 ! first year in stream to use
  integer, parameter :: stream_meped_year_last  = 2003 ! last year in stream to use
  integer, parameter :: stream_meped_year_align = 2003 ! align stream_year_firstndep with

  logical :: stream_meped_is_initialized = .false.
  logical,protected :: stream_meped_is_active = .true.

contains


  subroutine meped_input_stream_init(model_mesh, model_clock, rc)
    use dshr_strdata_mod, only: shr_strdata_init_from_inline

    ! input/output variables
    type(ESMF_CLock), intent(in)  :: model_clock
    type(ESMF_Mesh) , intent(in)  :: model_mesh
    integer         , intent(out) :: rc

    rc = ESMF_SUCCESS

    ! Initialize northern hemisphere input stream
    call shr_strdata_init_from_inline(stream_north,                 &
         my_task             = iam,                                 &
         logunit             = iulog,                               &
         compname            = 'ATM',                               &
         model_clock         = model_clock,                         &
         model_mesh          = model_mesh,                          &
         stream_meshfile     = trim(meped_north_mesh),              &
         stream_filenames    = (/trim(meped_filepath)/),            &
         stream_yearFirst    = stream_meped_year_first,             &
         stream_yearLast     = stream_meped_year_last,              &
         stream_yearAlign    = stream_meped_year_align,             &
         stream_fldlistFile  = meped_north_varlist,                 &
         stream_fldListModel = meped_north_varlist,                 &
         stream_lev_dimname  = 'null',                              &
         stream_mapalgo      = 'consd',                             &
         stream_offset       = 0,                                   &
         stream_taxmode      = 'limit',                             &
         stream_dtlimit      = 1.0e30_r8,                           &
         stream_tintalgo     = 'linear',                            &
         stream_name         = 'Northern hemisphere MEPED inputs ', &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Initialize southern hemisphere input stream
    call shr_strdata_init_from_inline(stream_south,                 &
         my_task             = iam,                                 &
         logunit             = iulog,                               &
         compname            = 'ATM',                               &
         model_clock         = model_clock,                         &
         model_mesh          = model_mesh,                          &
         stream_meshfile     = trim(meped_south_mesh),              &
         stream_filenames    = (/trim(meped_filepath)/),            &
         stream_yearFirst    = stream_meped_year_first,             &
         stream_yearLast     = stream_meped_year_last,              &
         stream_yearAlign    = stream_meped_year_align,             &
         stream_fldlistFile  = meped_south_varlist,                 &
         stream_fldListModel = meped_south_varlist,                 &
         stream_lev_dimname  = 'null',                              &
         stream_mapalgo      = 'consd',                             &
         stream_offset       = 0,                                   &
         stream_taxmode      = 'limit',                             &
         stream_dtlimit      = 1.0e30_r8,                           &
         stream_tintalgo     = 'linear',                            &
         stream_name         = 'Southern hemisphere MEPED inputs ', &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    allocate(meped_e_flux(pcols,begchunk:endchunk))
    meped_e_flux = 0._r8
    allocate(meped_e_ekev(pcols,begchunk:endchunk))
    meped_e_ekev = 0._r8
    allocate(meped_p_flux(pcols,begchunk:endchunk))
    meped_p_flux = 0._r8
    allocate(meped_p_ekev(pcols,begchunk:endchunk))
    meped_p_ekev = 0._r8

  end subroutine meped_input_stream_init

  subroutine meped_input_stream_advance(rc)
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance
    use time_manager     , only : get_curr_date
    use phys_grid        , only : get_ncols_p

    integer, intent(out)  :: rc

    integer :: year
    integer :: mon
    integer :: day
    integer :: sec     ! seconds into current day
    integer :: mcdate  ! Current model date (yyyymmdd)

    real(r8), pointer :: dataptr1n(:)
    real(r8), pointer :: dataptr2n(:)
    real(r8), pointer :: dataptr3n(:)
    real(r8), pointer :: dataptr4n(:)

    real(r8), pointer :: dataptr1s(:)
    real(r8), pointer :: dataptr2s(:)
    real(r8), pointer :: dataptr3s(:)
    real(r8), pointer :: dataptr4s(:)

    integer :: i, c, n, ncol

    real(r8) :: tmp1(pcols), tmp2(pcols), tmp3(pcols), tmp4(pcols)

    ! Advance input streams
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day

    call shr_strdata_advance(stream_north, ymd=mcdate, tod=sec, logunit=iulog, istr='north-meped', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    call shr_strdata_advance(stream_south, ymd=mcdate, tod=sec, logunit=iulog, istr='south-meped', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Get pointer for stream data that is time and spatially interpolated to model time and grid
    call dshr_fldbun_getFldPtr(stream_north%pstrm(1)%fldbun_model, meped_north_varlist(1), fldptr1=dataptr1n, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_north%pstrm(1)%fldbun_model, meped_north_varlist(2), fldptr1=dataptr2n, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_north%pstrm(1)%fldbun_model, meped_north_varlist(3), fldptr1=dataptr3n, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_north%pstrm(1)%fldbun_model, meped_north_varlist(4), fldptr1=dataptr4n, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    n = 0
    do c = begchunk,endchunk
       ncol = get_ncols_p(c)
       do i = 1,ncol
          n = n+1
          tmp1(i) = dataptr1n(n)
          tmp2(i) = dataptr2n(n)
          tmp3(i) = dataptr3n(n)
          tmp4(i) = dataptr4n(n)
       end do
       call outfld(meped_north_varlist(1), tmp1, pcols, c)
       call outfld(meped_north_varlist(2), tmp2, pcols, c)
       call outfld(meped_north_varlist(3), tmp3, pcols, c)
       call outfld(meped_north_varlist(4), tmp4, pcols, c)
    end do


    call dshr_fldbun_getFldPtr(stream_south%pstrm(1)%fldbun_model, meped_south_varlist(1), fldptr1=dataptr1s, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_south%pstrm(1)%fldbun_model, meped_south_varlist(2), fldptr1=dataptr2s, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_south%pstrm(1)%fldbun_model, meped_south_varlist(3), fldptr1=dataptr3s, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_south%pstrm(1)%fldbun_model, meped_south_varlist(4), fldptr1=dataptr4s, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    n = 0
    do c = begchunk,endchunk
       ncol = get_ncols_p(c)
       do i = 1,ncol
          n = n+1
          tmp1(i) = dataptr1s(n)
          tmp2(i) = dataptr2s(n)
          tmp3(i) = dataptr3s(n)
          tmp4(i) = dataptr4s(n)
       end do
       call outfld(meped_south_varlist(1), tmp1, pcols, c)
       call outfld(meped_south_varlist(2), tmp2, pcols, c)
       call outfld(meped_south_varlist(3), tmp3, pcols, c)
       call outfld(meped_south_varlist(4), tmp4, pcols, c)
    end do

    n = 0
    do c = begchunk,endchunk
       ncol = get_ncols_p(c)
       do i = 1,ncol
          n = n+1
          meped_e_flux(i,c) = dataptr1n(n) + dataptr1s(n)
          meped_e_ekev(i,c) = dataptr2n(n) + dataptr2s(n)
          meped_p_flux(i,c) = dataptr3n(n) + dataptr3s(n)
          meped_p_ekev(i,c) = dataptr4n(n) + dataptr4s(n)
       end do
    end do

    rc = ESMF_SUCCESS

  end subroutine meped_input_stream_advance

end module meped_input_stream
