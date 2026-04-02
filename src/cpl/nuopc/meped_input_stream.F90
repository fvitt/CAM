!%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
! Manages input data streams of PEMED mapped fluxes and mean energies
!
! Uses CDEP's input data stream utility which handles time interpolation
! and mapping to the physics grid
!%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%
module meped_input_stream

  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
  use shr_log_mod, only: errMsg => shr_log_errMsg
  use dshr_strdata_mod, only: shr_strdata_type
  use spmd_utils, only: iam, masterproc
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun
  use ppgrid, only: pcols, begchunk,endchunk

  use ESMF, only: ESMF_Clock, ESMF_Mesh
  use ESMF, only: ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT
  use ESMF, only: ESMF_Finalize, ESMF_LogFoundError

  implicit none

  private
  public :: meped_input_stream_readnl
  public :: meped_input_stream_init
  public :: meped_input_stream_advance
  public :: meped_input_stream_final
  public :: meped_input_is_active
  public :: meped_e_flux
  public :: meped_p_flux
  public :: meped_e_ekev
  public :: meped_p_ekev

  logical,protected :: meped_input_is_active = .false.

  real(r8),protected, pointer :: meped_e_flux(:,:) => null()
  real(r8),protected, pointer :: meped_p_flux(:,:) => null()
  real(r8),protected, pointer :: meped_e_ekev(:,:) => null()
  real(r8),protected, pointer :: meped_p_ekev(:,:) => null()

  integer :: stream_meped_year_first = -huge(1) ! 2003 ! first year in stream to use
  integer :: stream_meped_year_last  = -huge(1) ! 2003 ! last year in stream to use
  integer :: stream_meped_year_align = -huge(1) ! 2003 ! align stream_meped_year_first with

  character(len=CL) :: meped_filepath = 'NONE'
  character(len=CL) :: meped_north_mesh = 'NONE'
  character(len=CL) :: meped_south_mesh = 'NONE'

  type(shr_strdata_type) :: stream_north ! input data stream for northern hemisphere
  type(shr_strdata_type) :: stream_south ! input data stream for southern hemisphere

  character(len=*), parameter :: north_varlist(*) = (/'e_flux_nh','e_ekev_nh','p_flux_nh','p_ekev_nh'/)
  character(len=*), parameter :: south_varlist(*) = (/'e_flux_sh','e_ekev_sh','p_flux_sh','p_ekev_sh'/)

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine meped_input_stream_readnl(nlfile)
    use shr_nl_mod, only: shr_nl_find_group_name
    use spmd_utils, only: mpicom, mpi_character, mpi_integer, mpi_success

    ! input/output variables
    character(len=*), intent(in) :: nlfile

    ! local variables
    integer                 :: nu_nml                 ! unit for namelist file
    integer                 :: nml_error              ! namelist i/o error flag
    integer                 :: ierr
    character(*), parameter :: subname = 'meped_input_stream_readnl: '

    namelist /meped_stream_nl/ &
         meped_filepath, &
         meped_north_mesh, &
         meped_south_mesh, &
         stream_meped_year_first, &
         stream_meped_year_last, &
         stream_meped_year_align

    ! Read meped_stream_nl namelist
    if (masterproc) then
       open( newunit=nu_nml, file=trim(nlfile), status='old', iostat=nml_error )
       if (nml_error /= 0) then
          call endrun(subname//'ERROR opening '//trim(nlfile)//errMsg(__FILE__, __LINE__))
       end if
       call shr_nl_find_group_name(nu_nml, 'meped_stream_nl', status=nml_error)
       if (nml_error == 0) then
          read(nu_nml, nml=meped_stream_nl, iostat=nml_error)
          if (nml_error /= 0) then
             call endrun(subname//'ERROR reading meped_stream_nl namelist'//errMsg(__FILE__, __LINE__))
          end if
       end if
       close(nu_nml)
    endif

    call mpi_bcast(meped_filepath, len(meped_filepath), mpi_character, 0, mpicom, ierr)
    if (ierr/=mpi_success) call endrun(trim(subname)//"FATAL: mpi_bcast: meped_filepath")
    call mpi_bcast(meped_north_mesh, len(meped_north_mesh), mpi_character, 0, mpicom, ierr)
    if (ierr/=mpi_success) call endrun(trim(subname)//"FATAL: mpi_bcast: meped_north_mesh")
    call mpi_bcast(meped_south_mesh, len(meped_south_mesh), mpi_character, 0, mpicom, ierr)
    if (ierr/=mpi_success) call endrun(trim(subname)//"FATAL: mpi_bcast: meped_south_mesh")

    call mpi_bcast(stream_meped_year_first, 1, mpi_integer, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//"FATAL: mpi_bcast: stream_meped_year_first")
    call mpi_bcast(stream_meped_year_last, 1, mpi_integer, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//"FATAL: mpi_bcast: stream_meped_year_last")
    call mpi_bcast(stream_meped_year_align, 1, mpi_integer, 0, mpicom, ierr)
    if (ierr /= 0) call endrun(trim(subname)//"FATAL: mpi_bcast: stream_meped_year_align")

    meped_input_is_active = (meped_filepath/='NONE') .and. (len_trim(meped_filepath)>0)

    if (masterproc) then
       write(iulog,*) subname,'meped_input_is_active = ',meped_input_is_active
       if (meped_input_is_active) then
          write(iulog,*) subname,'meped_filepath = ',trim(meped_filepath)
          write(iulog,*) subname,'meped_north_mesh = ',trim(meped_north_mesh)
          write(iulog,*) subname,'meped_south_mesh = ',trim(meped_south_mesh)
          write(iulog,*) subname,'stream_meped_year_first = ',stream_meped_year_first
          write(iulog,*) subname,'stream_meped_year_last = ',stream_meped_year_last
          write(iulog,*) subname,'stream_meped_year_align = ',stream_meped_year_align
       end if
    end if

  end subroutine meped_input_stream_readnl

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine meped_input_stream_init(model_mesh, model_clock, rc)
    use dshr_strdata_mod, only: shr_strdata_init_from_inline

    ! input/output variables
    type(ESMF_CLock), intent(in)  :: model_clock
    type(ESMF_Mesh) , intent(in)  :: model_mesh
    integer         , intent(out) :: rc

    integer :: astat ! allocate status
    character(*), parameter :: subname = 'meped_input_stream_init: '

    rc = ESMF_SUCCESS
    if (.not.meped_input_is_active) return

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
         stream_fldlistFile  = north_varlist,                       &
         stream_fldListModel = north_varlist,                       &
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
         stream_fldlistFile  = south_varlist,                       &
         stream_fldListModel = south_varlist,                       &
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

    allocate(meped_e_flux(pcols,begchunk:endchunk), stat=astat)
    if( astat /= 0 ) call endrun(subname//'Failed to allocate meped_e_flux')
    meped_e_flux(:,:) = 0._r8

    allocate(meped_e_ekev(pcols,begchunk:endchunk), stat=astat)
    if( astat /= 0 ) call endrun(subname//'Failed to allocate meped_e_ekev')
    meped_e_ekev(:,:) = 0._r8

    allocate(meped_p_flux(pcols,begchunk:endchunk), stat=astat)
    if( astat /= 0 ) call endrun(subname//'Failed to allocate meped_p_flux')
    meped_p_flux(:,:) = 0._r8

    allocate(meped_p_ekev(pcols,begchunk:endchunk), stat=astat)
    if( astat /= 0 ) call endrun(subname//'Failed to allocate meped_p_ekev')
    meped_p_ekev(:,:) = 0._r8

  end subroutine meped_input_stream_init

  !-----------------------------------------------------------------------------
  ! Updates MEPED input streams for current model date and time
  !-----------------------------------------------------------------------------
  subroutine meped_input_stream_advance()
    use dshr_methods_mod , only : dshr_fldbun_getfldptr
    use dshr_strdata_mod , only : shr_strdata_advance
    use time_manager     , only : get_curr_date
    use phys_grid        , only : get_ncols_p

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

    integer :: rc
    integer :: i, c, n, ncol

    real(r8) :: tmp1(pcols), tmp2(pcols), tmp3(pcols), tmp4(pcols)

    if (.not.meped_input_is_active) return

    ! get current model date, time
    call get_curr_date(year, mon, day, sec)
    mcdate = year*10000 + mon*100 + day

    ! update input streams
    call shr_strdata_advance(stream_north, ymd=mcdate, tod=sec, logunit=iulog, istr='north-meped', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call shr_strdata_advance(stream_south, ymd=mcdate, tod=sec, logunit=iulog, istr='south-meped', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! Get pointers for stream data that is time and spatially interpolated to model time and grid
    ! southern hemisphere
    call dshr_fldbun_getFldPtr(stream_north%pstrm(1)%fldbun_model, north_varlist(1), fldptr1=dataptr1n, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_north%pstrm(1)%fldbun_model, north_varlist(2), fldptr1=dataptr2n, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_north%pstrm(1)%fldbun_model, north_varlist(3), fldptr1=dataptr3n, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_north%pstrm(1)%fldbun_model, north_varlist(4), fldptr1=dataptr4n, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! southern hemisphere
    call dshr_fldbun_getFldPtr(stream_south%pstrm(1)%fldbun_model, south_varlist(1), fldptr1=dataptr1s, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_south%pstrm(1)%fldbun_model, south_varlist(2), fldptr1=dataptr2s, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_south%pstrm(1)%fldbun_model, south_varlist(3), fldptr1=dataptr3s, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    call dshr_fldbun_getFldPtr(stream_south%pstrm(1)%fldbun_model, south_varlist(4), fldptr1=dataptr4s, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    ! set gridded MEPED arrays
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

  end subroutine meped_input_stream_advance

  !-----------------------------------------------------------------------------
  ! free up allocated memory, etc.
  !-----------------------------------------------------------------------------
  subroutine meped_input_stream_final()
    deallocate(meped_e_flux)
    deallocate(meped_p_flux)
    deallocate(meped_e_ekev)
    deallocate(meped_p_ekev)
    nullify(meped_e_flux)
    nullify(meped_p_flux)
    nullify(meped_e_ekev)
    nullify(meped_p_ekev)
  end subroutine meped_input_stream_final

end module meped_input_stream
