module meped_module
!
! Import electron data from MEPED (Medium Energy Proton and Electron
! Detector) instrument on NOAA POES (TIROS) satellite.
!
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl => shr_kind_cl
  use cam_logfile,    only: iulog
  use spmd_utils,     only: masterproc
  use edyn_maggrid,   only: nmlat, nmlonp1
  use edyn_maggrid,   only: ylonm     ! magnetic latitudes (nmlat) (radians)
  use edyn_maggrid,   only: ylatm     ! magnetic longtitudes (nmlonp1) (radians)
  use cam_pio_utils,  only: cam_pio_openfile, cam_pio_closefile
  use pio,            only: pio_inq_dimid, pio_inquire_dimension
  use pio,            only: pio_inquire, pio_inq_varid
  use pio,            only: file_desc_t, pio_noerr, pio_nowrite, pio_get_var
  use utils_mod,      only: check_ncerr, check_alloc, boxcar_ave
  use edyn_mpi,       only: ntask, mytid
  use edyn_params,    only: pi, dtr, rtd
!  use input_data_utils, only: time_coordinate
  use time_manager,   only : get_curr_calday

  implicit none

  private
  public :: init_meped
  public :: getmeped, imeped, prescribed_meped_period
  logical  :: prescribed_meped_period= .false.
  integer  :: imeped

!
! Define MEPED output fields
  integer,allocatable, dimension(:), save :: year,jday,ut
  real(r8),allocatable, dimension(:), save :: glon,glat
  real(r8),allocatable, dimension(:), save :: time_meped
  real(r8),allocatable,dimension(:,:,:), save :: & ! ntimes,gnlon,gnlat
          e_flux_nh,e_flux_sh,e_ekev_nh,e_ekev_sh,      &
          p_flux_nh,p_flux_sh,p_ekev_nh,p_ekev_sh
!
  type(file_desc_t) :: ncid
  character(len=cl), allocatable :: meped_files(:)
  integer :: num_files, file_ndx
 ! type(time_coordinate) :: time_coord
  integer                      :: ntimes, gnlat,gnlon

      contains
!-----------------------------------------------------------------------
      subroutine init_meped(meped_list)
!
      character(len=*),intent(in) :: meped_list(:)

      integer :: n, nfiles

      nfiles = min( size(meped_list), size(meped_list) )
      num_files = 0

      count_files: do n = 1,nfiles
         if (len_trim(meped_list(n))<1 .or. trim(meped_list(n))=='NONE') then
            exit count_files
         else
            num_files = num_files + 1
         end if
      end do count_files

      allocate(meped_files(num_files))
      meped_files(:num_files) = meped_list(:num_files)
      file_ndx = 1
      call open_files()

      end subroutine init_meped
!-----------------------------------------------------------------------
      subroutine rdmeped(meped)
!
! Read MEPED data file (this is called once per run from init, if
! namelist read meped_file is non-blank).
!
   character(len=*), intent(in) :: meped
!
!  Local variables:
    integer                      :: istat, ndims, nvars, ngatts
    integer                      :: idunlim, ier
    integer                      :: id_time,idv_year,idv_day,idv_ut
    integer                      :: idv_glat,idv_glon,id_glat,id_glon
    character(len=*), parameter  :: subname = 'rdmeped'

    integer :: idv_e_flux_nh,idv_e_flux_sh,idv_p_flux_nh,     &
               idv_p_flux_sh,idv_e_ekev_nh,idv_e_ekev_sh,     &
               idv_p_ekev_nh,idv_p_ekev_sh
!
    if (masterproc) then
       write(iulog, "(/,72('-'))")
       write(iulog, "(a,': read MEPED data:')") subname
    end if
    !
!
! From ncdump (only electron vars):
!
! netcdf meped {
! dimensions:
!         time = 6992 ;
!variables:
! float GLON(structure_elements, dim1_GLON) ;
! float GLAT(structure_elements, dim1_GLAT) ;
! short YEAR(structure_elements, dim1_YEAR) ;
!         YEAR:lon_name = "year" ;
! short JDAY(structure_elements, dim1_JDAY) ;
!         JDAY:lon_name = "Julian day of year" ;
! float UT(structure_elements, dim1_UT) ;
!         UT:long_name = "universal time (hr)" ;
! float E_FLUX_NH(structure_elements, dim3_E_FLUX_NH, dim2_E_FLUX_NH, dim1_E_FLUX_NH) ;
!         E_FLUX_NH:lon_name = "Fitted Northern NOAA-MED Electron Energy Flux" ;
!         E_FLUX_NH:units = "mW m^-2" ;
! float E_EKEV_NH(structure_elements, dim3_E_EKEV_NH, dim2_E_EKEV_NH, dim1_E_EKEV_NH) ;
!         E_EKEV_NH:lon_name = "Fitted Northern NOAA-MED Electron Mean Energy" ;
!         E_EKEV_NH:units = "keV" ;
! float E_FLUX_SH(structure_elements, dim3_E_FLUX_SH, dim2_E_FLUX_SH, dim1_E_FLUX_SH) ;
!         E_FLUX_SH:lon_name = "Fitted Southern NOAA-MED Electron Energy Flux" ;
!         E_FLUX_SH:units = "mW m^-2" ;
! float E_EKEV_SH(structure_elements, dim3_E_EKEV_SH, dim2_E_EKEV_SH, dim1_E_EKEV_SH) ;
!         E_EKEV_SH:lon_name = "Fitted Southern NOAA-MED Electron Mean Energy" ;
!         E_EKEV_SH:units = "keV" ;
! float P_FLUX_NH(structure_elements, dim3_P_FLUX_NH, dim2_P_FLUX_NH, dim1_P_FLUX_NH) ;
!         P_FLUX_NH:lon_name = "Fitted Northern NOAA-MED Proton Energy Flux" ;
!         P_FLUX_NH:units = "mW m^-2" ;
! float P_EKEV_NH(structure_elements, dim3_P_EKEV_NH, dim2_P_EKEV_NH, dim1_P_EKEV_NH) ;
!         P_EKEV_NH:lon_name = "Fitted Northern NOAA-MED Proton Mean Energy" ;
!         P_EKEV_NH:units = "keV" ;
! float P_FLUX_SH(structure_elements, dim3_P_FLUX_SH, dim2_P_FLUX_SH, dim1_P_FLUX_SH) ;
!         P_FLUX_SH:lon_name = "Fitted Southern NOAA-MED Proton Energy Flux" ;
!         P_FLUX_SH:units = "mW m^-2" ;
! float P_EKEV_SH(structure_elements, dim3_P_EKEV_SH, dim2_P_EKEV_SH, dim1_P_EKEV_SH) ;
!         P_EKEV_SH:lon_name = "Fitted Southern NOAA-MED Proton Mean Energy" ;
!         el_char:units = "KeV" ;
!
    ! Open netcdf file:
    call cam_pio_openfile(ncid, meped, pio_nowrite)
    !
    ! Get MEPED grid dimension:
    istat = pio_inq_dimid(ncid, 'nlon', id_glon)
    istat = pio_inquire_dimension(ncid, id_glon, len=gnlon)
    call check_ncerr(istat, subname, 'MEPED longitude ndimension')

    istat = pio_inq_dimid(ncid, 'nlat', id_glat)
    istat = pio_inquire_dimension(ncid, id_glat, len=gnlat)
    call check_ncerr(istat, subname, 'MEPED latitude dimension')

!    call time_coord%initialize( meped, set_weights=.false. )

    !
    ! Get time dimension:
    istat = pio_inq_dimid(ncid, 'ntime', id_time)
    istat = pio_inquire_dimension(ncid, id_time, len=ntimes)
    call check_ncerr(istat, subname, 'MEPED time dimension')

    if (masterproc) then
       write(iulog,"('rdmeped: ntime,gnlon,gnlat =',3i7)") ntimes,gnlon,gnlat
    endif

    if (.not. allocated(time_meped)) then
       allocate(time_meped(ntimes), stat=ier)
       call check_alloc(ier, subname, 'time_meped', ntimes=ntimes)
    end if

    if (.not. allocated(YEAR)) then
       allocate(year(ntimes), stat=ier)
       call check_alloc(ier, subname, 'year', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid, 'year', idv_year)
    call check_ncerr(istat, subname, 'MEPED year id')
    istat = pio_get_var(ncid, idv_year, year)
    call check_ncerr(istat, subname, 'MEPED year')
    if (.not. allocated(JDAY)) then
       allocate(jday(ntimes), stat=ier)
       call check_alloc(ier, subname, 'jday', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid, 'jday', idv_day)
    call check_ncerr(istat, subname, 'MEPED jday id')
    istat = pio_get_var(ncid, idv_day, jday)
    call check_ncerr(istat, subname, 'MEPED jday')
    if (.not. allocated(ut)) then
       allocate(ut(ntimes), stat=ier)
       call check_alloc(ier, subname, 'ut', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid, 'ut', idv_ut)
    call check_ncerr(istat, subname, 'MEPED ut id')
    istat = pio_get_var(ncid, idv_ut, ut)
    call check_ncerr(istat, subname, 'MEPED ut')

    !
    ! Allocate 3-d fields:
    if (.not. allocated(e_flux_nh)) then
       allocate(e_flux_nh(ntimes,gnlon,gnlat), stat=ier)
       call check_alloc(ier, subname, 'e_flux_nh id')
    end if
    if (.not. allocated(e_flux_sh)) then
       allocate(e_flux_sh(ntimes,gnlon,gnlat), stat=ier)
       call check_alloc(ier, subname, 'e_flux_sh id')
    end if
    if (.not. allocated(e_ekev_nh)) then
       allocate(e_ekev_nh(ntimes,gnlon,gnlat), stat=ier)
       call check_alloc(ier, subname, 'e_ekev_nh id')
    end if
    if (.not. allocated(e_ekev_sh)) then
       allocate(e_ekev_sh(ntimes,gnlon,gnlat), stat=ier)
       call check_alloc(ier, subname, 'e_ekev_sh id')
    end if
    if (.not. allocated(p_flux_nh)) then
       allocate(p_flux_nh(ntimes,gnlon,gnlat), stat=ier)
       call check_alloc(ier, subname, 'p_flux_nh id')
    end if
    if (.not. allocated(p_flux_sh)) then
       allocate(p_flux_sh(ntimes,gnlon,gnlat), stat=ier)
       call check_alloc(ier, subname, 'p_flux_sh id')
    end if
    if (.not. allocated(p_ekev_nh)) then
       allocate(p_ekev_nh(ntimes,gnlon,gnlat), stat=ier)
       call check_alloc(ier, subname, 'p_ekev_nh id')
    end if
    if (.not. allocated(p_ekev_sh)) then
       allocate(p_ekev_sh(ntimes,gnlon,gnlat), stat=ier)
       call check_alloc(ier, subname, 'p_ekev_sh id')
    end if

!
! Get MEPED mean energy and energy flux for electrons and protons

    ! ELECTRONS
    istat = pio_inq_varid(ncid, 'e_flux_nh', idv_e_flux_nh)
    call check_ncerr(istat, subname, 'inq_varid idv_e_flux_nh')
    istat = pio_get_var(ncid, idv_e_flux_nh, e_flux_nh)
    call check_ncerr(istat, subname, 'get_var e_flux_nh')

    istat = pio_inq_varid(ncid, 'e_flux_sh', idv_e_flux_sh)
    call check_ncerr(istat, subname, 'inq_varid idv_e_flux_sh')
    istat = pio_get_var(ncid, idv_e_flux_nh, e_flux_sh)
    call check_ncerr(istat, subname, 'get_var e_flux_sh')

    istat = pio_inq_varid(ncid, 'e_ekev_nh', idv_e_ekev_nh)
    call check_ncerr(istat, subname, 'inq_varid idv_e_ekev_nh')
    istat = pio_get_var(ncid, idv_e_ekev_nh, e_ekev_nh)
    call check_ncerr(istat, subname, 'get_var e_ekev_nh')

    istat = pio_inq_varid(ncid, 'e_ekev_sh', idv_e_ekev_sh)
    call check_ncerr(istat, subname, 'inq_varid idv_e_ekev_sh')
    istat = pio_get_var(ncid, idv_e_ekev_nh, e_ekev_sh)
    call check_ncerr(istat, subname, 'get_var e_ekev_sh')

    ! PROTONS
    istat = pio_inq_varid(ncid, 'p_flux_nh', idv_p_flux_nh)
    call check_ncerr(istat, subname, 'inq_varid idv_p_flux_nh')
    istat = pio_get_var(ncid, idv_p_flux_nh, p_flux_nh)
    call check_ncerr(istat, subname, 'get_var p_flux_nh')

    istat = pio_inq_varid(ncid, 'p_flux_sh', idv_p_flux_sh)
    call check_ncerr(istat, subname, 'inq_varid idv_p_flux_sh')
    istat = pio_get_var(ncid, idv_p_flux_nh, p_flux_sh)
    call check_ncerr(istat, subname, 'get_var p_flux_sh')

    istat = pio_inq_varid(ncid, 'p_ekev_nh', idv_p_ekev_nh)
    call check_ncerr(istat, subname, 'inq_varid idv_p_ekev_nh')
    istat = pio_get_var(ncid, idv_p_ekev_nh, p_ekev_nh)
    call check_ncerr(istat, subname, 'get_var p_ekev_nh')

    istat = pio_inq_varid(ncid, 'p_ekev_sh', idv_p_ekev_sh)
    call check_ncerr(istat, subname, 'inq_varid idv_p_ekev_sh')
    istat = pio_get_var(ncid, idv_p_ekev_nh, p_ekev_sh)
    call check_ncerr(istat, subname, 'get_var p_ekev_sh')

  end subroutine rdmeped
!-----------------------------------------------------------------------
  subroutine getmeped(model_year,model_jday,rmodel_secs,iprint,imeped, &
                   E_meped_flux,E_meped_ekev,P_meped_flux,P_meped_ekev)
!
! Interpolate MEPED data to current model time.
! This is called at every timestep from aurora.F, if namelist read
!   meped_file is non-blank (data was read at run initialization
!   by sub rdmeped).
!
   use cam_history_support, only: fillvalue
   use rgrd_mod,            only: rgrd2
   use edyn_geogrid,        only: nlat, nlon, nlonp1, nlatp1, ylatg, ylong
! note: ylatg & ylong are in radians, ylatg from -pi/2 to pi/2, ylong frpm -pi to pi
!
! Args:
      integer,intent(in) :: model_year,model_jday,iprint
      integer,intent(out) :: imeped
      real(r8),intent(in) :: rmodel_secs
      real(r8), intent(out) ::       & ! on geographic grid
        E_meped_flux(nlonp1,nlat),P_meped_flux(nlonp1,nlat), &
        E_meped_ekev(nlonp1,nlat),P_meped_ekev(nlonp1,nlat)

!
! Local:
      integer :: i,j,i0,i1,ind,jnd,ndpy,nspy,meped_debug
      integer,parameter :: ngglat = 60,ngglon=180
      integer,parameter :: nxglat = 180,nxglon=ngglon
      real(r8) :: gglat(ngglat),gglon(ngglon),ylatgg(nlatp1+1)
      real(r8) :: xglat(nxglat),xglon(ngglon)
      real(r8) :: total_model_secs,f0,f1
      real(r8),dimension(gnlon,gnlat) ::                        &
       ave_e_flux_nh,ave_e_flux_sh,ave_e_ekev_nh,ave_e_ekev_sh, &
       ave_p_flux_nh,ave_p_flux_sh,ave_p_ekev_nh,ave_p_ekev_sh
      real(r8),dimension(nxglon,nxglat) ::                      &
       ave_e_flux,ave_p_flux,ave_e_ekev,ave_p_ekev
      integer(kind=8) :: meped_total_sec0,meped_total_sec1,     &
        model_total_secs
      integer :: ier,lw,liw,isign,intpol(2)
      integer,allocatable :: iw(:)
      real(r8),allocatable :: w(:)
      real(r8) ::                                                    &
        E_fit_flux(nlonp1,nlat),P_fit_flux(nlonp1,nlat), &
        E_fit_ekev(nlonp1,nlat),P_fit_ekev(nlonp1,nlat), &
        yout(nlonp1,nlat)
      character(len=*), parameter :: subname = 'getmeped'
!
      E_meped_flux(:,:) = 0._r8
      E_meped_ekev(:,:) = 30._r8
      P_meped_flux(:,:) = 0._r8
      P_meped_ekev(:,:) = 30._r8
!
! Check times:
!
      if (year(1) /= model_year) then
        write(6,"('>>> getmeped: wrong year,model_year=',2i4)")         &
           year(1),model_year
        stop 'getmeped'
      endif
! Total model seconds (since Jan 1, year 0000)
      ndpy = 365
      if (mod(model_year,4)==0) ndpy = 366 ! leap year
      nspy = 3600*24*ndpy ! number of seconds per year
      total_model_secs = float(nspy)*float(model_year)+         &
        3600.*24.*float(model_jday)+rmodel_secs
      time_meped(:) = float(nspy)*float(year(:))+               &
        3600.*24.*float(jday(:))+3600.*float(ut(:))
!
! Bracket model time:
      i0 = 0
      do i=2,ntimes
        if (total_model_secs <= time_meped(i) .and.             &
            total_model_secs >= time_meped(i-1)) then
          i0 = i-1
          i1 = i
        endif
      enddo ! i=1,ntimes
      if (total_model_secs >= time_meped(ntimes)) then
        i0=ntimes
        i1=i0
      endif
      if (i0==0 .or. i0==ntimes) then
         if (masterproc) then
            write(iulog,"(/,'>>> getmeped: could not bracket model ',   &
                 'time: year=',i4,' day=',i20,' secs=',f10.2)")         &
                 model_year,model_jday,rmodel_secs
            write(iulog,"('MEPED data file contains data for ',         &
                 'year,day ',i4,', ',i3,' to ',i4,', ',i3)")           &
                 year(1),jday(1),year(ntimes),jday(ntimes)
            !       call shutdown('getmeped')
            write(iulog,"('>>> getmeped: Am turning meped data off ','(setting imeped=0)')")
         end if
        imeped = 0
        return
      else
         imeped = 1
         if (masterproc) then
            write(iulog,"('getmeped: bracket model time: i0=',i5,' i1=',i5, &
                 ' time_meped(i0)=',f15.2,' time_meped(i1)=',f15.2)")      &
                 i0,i1,time_meped(i0),time_meped(i1)
         end if
      endif
!
! Linear interpolation (8-byte ints for total seconds):
!
      meped_total_sec0 = int(time_meped(i0),8)
      meped_total_sec1 = int(time_meped(i1),8)
      model_total_secs = int(total_model_secs,8)
      f0 = model_total_secs-meped_total_sec0
      f1 = meped_total_sec1-meped_total_sec0
      if (total_model_secs >= time_meped(ntimes)) f1=1
!       write(6,"('getmeped: bracket model time: i0=',i5,' i1=',i5,     &
!         ' time_meped(i0)=',f15.2,' time_meped(i1)=',f15.2,            &
!         ' f0=',f15.2,' f1=',f15.2)") i0,i1,time_meped(i0),time_meped(i1),f0,f1
      ave_e_flux_nh(:,:) = e_flux_nh(i0,:,:)+(e_flux_nh(i1,:,:)-e_flux_nh(i0,:,:))*f0/f1
      ave_e_flux_sh(:,:) = e_flux_sh(i0,:,:)+(e_flux_sh(i1,:,:)-e_flux_sh(i0,:,:))*f0/f1
      ave_e_ekev_nh(:,:) = e_ekev_nh(i0,:,:)+(e_ekev_nh(i1,:,:)-e_ekev_nh(i0,:,:))*f0/f1
      ave_e_ekev_sh(:,:) = e_ekev_sh(i0,:,:)+(e_ekev_sh(i1,:,:)-e_ekev_sh(i0,:,:))*f0/f1
      ave_p_flux_nh(:,:) = p_flux_nh(i0,:,:)+(p_flux_nh(i1,:,:)-p_flux_nh(i0,:,:))*f0/f1
      ave_p_flux_sh(:,:) = p_flux_sh(i0,:,:)+(p_flux_sh(i1,:,:)-p_flux_sh(i0,:,:))*f0/f1
      ave_p_ekev_nh(:,:) = p_ekev_nh(i0,:,:)+(p_ekev_nh(i1,:,:)-p_ekev_nh(i0,:,:))*f0/f1
      ave_p_ekev_sh(:,:) = p_ekev_sh(i0,:,:)+(p_ekev_sh(i1,:,:)-p_ekev_sh(i0,:,:))*f0/f1

      if (total_model_secs >= time_meped(ntimes)) then
         if (masterproc) then
            write(iulog,"('getmeped: model year,day,secs=',i4,', ',i3,', ',i6,    &
                 ' model total secs=',i12)") model_year,model_jday,              &
                 int(rmodel_secs),model_total_secs
            write(iulog,"('getmeped: meped year0,day0,ut0=',i4,', ',i3,', ',      &
                 i6,' meped total secs0=',i15)") year(i0),jday(i0),              &
                 ut(i0),meped_total_sec0
            write(iulog,"('getmeped: meped year1,day1,ut1=',i4,', ',i3,', ',      &
                 i6,' meped total secs1=',i15)") year(i1),jday(i1),              &
                 ut(i1),meped_total_sec1
            !     write(iulog,"('elec_pw0,1  =',2e12.4,' epower=',e12.4)") elec_pw(i0), &
            !       elec_pw(i1),epower
            !     write(iulog,"('sec0, sec1, secs, sec1-sec0, secs-sec0 = ',5i15)")     &
            !      meped_total_sec0,meped_total_sec1,model_total_secs,              &
            !      meped_total_sec1-meped_total_sec0,model_total_secs-meped_total_sec0
            write(iulog,"('max e_flux_nh=',2e12.4,' max ave_e_flux_nh =',e12.4)") &
                 maxval(e_flux_nh(i0,:,:)),maxval(e_flux_nh(i1,:,:)),            &
                 maxval(ave_e_flux_nh(:,:))
            write(iulog,"('max p_flux_nh=',2e12.4,' max ave_p_flux_nh =',e12.4)") &
                 maxval(p_flux_nh(i0,:,:)),maxval(p_flux_nh(i1,:,:)),            &
                 maxval(ave_p_flux_nh(:,:))
         endif
      endif

! Convert to WACCMX 2-D geographic grid
!Note: the MEPED longitude from (1.,3., ...179,-179,....,-1), needs to be reorganized to
! to WACCMX glon which extands from -180 to 180.
      do i=1,nxglon
        xglon(i) = -181._r8 +i*2._r8
      enddo
      do j=1,nxglat
        xglat(j) = j*1._r8 - 90.5_r8
        ave_e_flux(:,j) = 0._r8
        ave_e_ekev(:,j) = 30._r8
        ave_p_flux(:,j) = 0._r8
        ave_p_ekev(:,j) = 30._r8
        do i=1,nxglon
          ind = int(xglon(i)/2._r8)
          if (xglon(i) < 0.) ind = int((360._r8+xglon(i))/2._r8)
          if (xglat(j) < -30._r8) then
              jnd = int((-29.5_r8-xglat(j))/1._r8)
              ave_e_flux(i,j) = ave_e_flux_sh(ind,jnd)
              ave_e_ekev(i,j) = ave_e_ekev_sh(ind,jnd)
              ave_p_flux(i,j) = ave_p_flux_sh(ind,jnd)
              ave_p_ekev(i,j) = ave_p_ekev_sh(ind,jnd)
          endif
          if (xglat(j) > 30._r8) then
              jnd = int((xglat(j) - 29.5_r8)/1._r8)
              ave_e_flux(i,j) = ave_e_flux_nh(ind,jnd)
              ave_e_ekev(i,j) = ave_e_ekev_nh(ind,jnd)
              ave_p_flux(i,j) = ave_p_flux_nh(ind,jnd)
              ave_p_ekev(i,j) = ave_p_ekev_nh(ind,jnd)
          endif
        enddo
      enddo

     !write(6,"('nxglat,nxglon,nlat,nlonp1 = ',4i10))") nxglat,nxglon,nlat,nlonp1
     !write(6,"('xglat = ',/,(10f10.4))") xglat(:)
     !write(6,"('xglon = ',/,(10f10.4))") xglon(:)
     !write(6,"('ylatg = ',/,(10f10.4))") ylatg(:)/dtr
     !write(6,"('ylong = ',/,(10f10.4))") ylong(:)/dtr

      xglat(:) = xglat(:)*dtr
      xglon(:) = xglon(:)*dtr
      xglat(1) = ylatg(1)
      xglat(nxglat) = ylatg(nlat)
      xglon(1) = ylong(1)
      xglon(nxglon) = ylong(nlonp1)

      lw = nlonp1+nlat+2*nlonp1
      liw = nlonp1 + nlat
      if (.not. allocated(w)) allocate(w(lw),stat=ier)
      if (ier /= 0) write(6,"('>>> horizontal_interp: error allocating',        &
        ' w(lw): lw=',i6,' ier=',i4)") lw,ier
      if (.not. allocated(iw)) allocate(iw(liw),stat=ier)
      if (ier /= 0) then
        write(6,"('>>> horzontal_interp: error allocating',     &
        ' iw(liw): liw=',i6,' ier=',i4)") liw,ier
       stop 'allocating liw'
      endif
      intpol(:) = 1 ! linear (not cubic) interp in both dimensions

      call rgrd2(nxglon,nxglat,xglon,xglat,ave_e_flux,nlonp1,nlat,      &
                 ylong,ylatg,E_fit_flux,intpol,w,lw,iw,liw,ier)
      if (ier /= 0) then
        write(6,"('>>> rgrd2: error in interpolating ier=',i4)") ier
       stop 'rgrd2'
      endif
!     write(6,"('max,min E_meped_flux = ',2g9.2)")      &
!         maxval(E_meped_flux(:,:)),minval(E_meped_flux(:,:))
      call rgrd2(nxglon,nxglat,xglon,xglat,ave_e_ekev,nlonp1,nlat,      &
                 ylong,ylatg,E_fit_ekev,intpol,w,lw,iw,liw,ier)
      call rgrd2(nxglon,nxglat,xglon,xglat,ave_p_flux,nlonp1,nlat,      &
                 ylong,ylatg,P_fit_flux,intpol,w,lw,iw,liw,ier)
      call rgrd2(nxglon,nxglat,xglon,xglat,ave_p_ekev,nlonp1,nlat,      &
                 ylong,ylatg,P_fit_ekev,intpol,w,lw,iw,liw,ier)

      call bin2d(E_fit_flux,yout,nlonp1,nlat)
      E_meped_flux(:,:) = yout(:,:)
      call bin2d(E_fit_ekev,yout,nlonp1,nlat)
      E_meped_ekev(:,:) = yout(:,:)
      call bin2d(P_fit_flux,yout,nlonp1,nlat)
      P_meped_flux(:,:) = yout(:,:)
      call bin2d(P_fit_ekev,yout,nlonp1,nlat)
      P_meped_ekev(:,:) = yout(:,:)

!
      if (iprint > 0) then
         if (masterproc) then
            write(iulog,"('getmeped: MEPED data interpolated to date and time')")
            write(iulog,"(72('-'),/)")
         endif
      endif

      end subroutine getmeped
!-----------------------------------------------------------------------
      subroutine bin2d(x,y,lon,lat)
!
! perform spacial average over the nearest 4 points
!
! Args:
      integer,intent(in) :: lon,lat
      real(r8),dimension(lon,lat), intent(in) :: x
      real(r8),dimension(lon,lat), intent(out) :: y
! Local:
      integer :: i, j
!
      y(:,:) = 0.
      do j=2,lat-1
        do i=3,lon-2
          y(i,j) = (x(i,j) + 0.5*(x(i-1,j)+x(i+1,j))+   &
                    0.25*(x(i-2,j)+x(i+2,j))+           &
                    0.25*(x(i,j-1)+x(i,j+1)))/3.
          y(i,1) = (x(i,1) + 0.5*(x(i-1,1)+x(i+1,1))+   &
                    0.25*(x(i-2,1)+x(i+2,1))+           &
                    0.25*(x(i,1)+x(i,2)))/3.
          y(i,lat) = (x(i,lat) + 0.5*(x(i-1,lat)+x(i+1,lat))+   &
                    0.25*(x(i-2,lat)+x(i+2,lat))+               &
                    0.25*(x(i,lat-1)+x(i,lat)))/3.
        enddo
        y(1,j) = (x(1,j) + 0.5*(x(lon,j)+x(2,j))+       &
                    0.25*(x(lon-1,j)+x(3,j))+           &
                    0.25*(x(1,j-1)+x(1,j+1)))/3.
        y(2,j) = (x(2,j) + 0.5*(x(1,j)+x(3,j))+         &
                    0.25*(x(lon,j)+x(4,j))+             &
                    0.25*(x(2,j-1)+x(2,j+1)))/3.
        y(lon,j) = (x(lon,j) + 0.5*(x(lon-1,j)+x(1,j))+ &
                    0.25*(x(lon-2,j)+x(2,j))+           &
                    0.25*(x(lon,j-1)+x(lon,j+1)))/3.
        y(lon-1,j) = (x(lon-1,j) + 0.5*(x(lon-2,j)+x(lon,j))+   &
                    0.25*(x(lon-3,j)+x(1,j))+                   &
                    0.25*(x(lon-1,j-1)+x(lon-1,j+1)))/3.
      enddo
!
      end subroutine bin2d
!-------------------------------------------------------------------
  subroutine close_files

    call cam_pio_closefile(ncid)

  end subroutine close_files
  !-----------------------------------------------------------------------
  subroutine open_files()

    call rdmeped(meped_files(file_ndx))

  end subroutine open_files
!--------------------------------------------------------------------
end module meped_module
