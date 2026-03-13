module helium_ubc_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_kind_mod, only: cl => shr_kind_cl
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils, only: masterproc, host_mpicom=>mpicom, masterprocid, mpi_character, mpi_success
  use ppgrid, only: begchunk, endchunk, pcols, pver, pverp
  use physics_types, only: physics_state

  use esmf_phys_mesh_mod, only: esmf_phys_mesh_init
  use esmf_lonlat_grid_mod, only: esmf_lonlat_grid_init, lon_beg,lon_end,lat_beg,lat_end
!  use esmf_zonal_fft_mod, only : esmf_zonal_fft_3d, esmf_zonal_fft_init
  use esmf_lonlatphys_regrid_mod, only: esmf_lonlatphys_regrid_init, regrid_lonlat2phys, regrid_phys2lonlat

  implicit none

  private
  public :: helium_ubc_readnl
  public :: helium_ubc_init
  public :: helium_ubc_calc
  public :: helium_ubc_fluxes

  real(r8), protected, allocatable :: helium_ubc_fluxes(:,:)

  character(len=cl) :: helium_ubc_coefs = 'NONE'

  integer :: mytid = -huge(1)

  real(r8), parameter :: gask = 8.314e7_r8      ! gas constant
  real(r8), parameter :: rmass_he = 4._r8
  real(r8), parameter :: grav = 870._r8
  real(r8), parameter :: pi = 4._r8*atan(1.0_r8)
  real(r8), parameter :: re = 6.37122e8_r8      ! earth radius (cm)
  real(r8), parameter :: p0 = 5.e-4_r8           ! standard pressure

  integer :: he_cnst_ndx = -1
  integer :: nlat_he = -1

  real(r8),dimension(:,:,:),allocatable :: pmn,zmn

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine helium_ubc_readnl(nlfile)
    use namelist_utils, only: find_group_name

    character(len=*), intent(in) :: nlfile

    integer :: unitn, ierr
    character(len=*), parameter :: prefix = 'helium_ubc_readnl: '

    namelist /helium_ubc_opts/ helium_ubc_coefs

    if (masterproc) then
       ! read namelist
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'helium_ubc_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, helium_ubc_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(prefix//'helium_ubc_opts: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! broadcast to all MPI tasks
    call mpi_bcast(helium_ubc_coefs, len(helium_ubc_coefs), mpi_character, masterprocid, host_mpicom, ierr)
    if (ierr /= mpi_success) call endrun(prefix//'mpi_bcast error : helium_ubc_coefs')

    ! write params to atm log
    if (masterproc) then
       write(iulog,*) prefix, 'helium_ubc_coefs file path: ',helium_ubc_coefs
    end if

  end subroutine helium_ubc_readnl

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine helium_ubc_init()

    use mpi, only: mpi_comm_size, mpi_comm_rank, mpi_comm_split
    use constituents, only: cnst_mw, cnst_get_ind
    use cam_history, only: addfld, horiz_only

    integer :: host_npes, he_npes
    integer :: ierr
    character(len=*), parameter :: prefix = 'helium_ubc_init: '

    call mpi_comm_size(host_mpicom, host_npes, ierr)
    if (ierr /= mpi_success) then
       call endrun(prefix//'MPI ERROR -- mpi_comm_size')
    end if

    call read_coefs_file()

    he_npes = min(host_npes, nlat_he/2)

    ! write to atm log
    if (masterproc) then
       write(iulog,*) prefix, 'host model npes : ', host_npes
       write(iulog,*) prefix, 'helium_ubc_npes : ', he_npes
       write(iulog,*) prefix, 'helium_ubc_nlats: ', nlat_he
       write(iulog,*) prefix, 'host model mpicom: ',host_mpicom
    end if

    call esmf_phys_mesh_init()
    call esmf_lonlat_grid_init(nlat_he)

    call esmf_lonlatphys_regrid_init()

!    call esmf_zonal_fft_init()

    allocate(helium_ubc_fluxes(pcols,begchunk:endchunk))
    helium_ubc_fluxes = 0._r8

    call cnst_get_ind( 'H', he_cnst_ndx )

    call addfld('HEFLUX_TST1', horiz_only,  'A', ' ', 'He Flux Test fld1' )
    call addfld('HEFLUX_TST2', horiz_only,  'A', ' ', 'He Flux Test fld2' )

  end subroutine helium_ubc_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine helium_ubc_calc(phys_state)
    use cam_history, only: outfld

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)


    real(r8) :: flx_phys(pcols,begchunk:endchunk)
    real(r8) :: tmp_phys(pcols,begchunk:endchunk)

    real(r8) :: flx_lonlat(lon_beg:lon_end,lat_beg:lat_end)
    real(r8) :: tn, he_mmr
    integer :: lchnk, ncol, i

!    flx_arg(:ncol) = -4._r8*p0*sqrt((gask*tni(:ncol)/(rmass_he*grav))**3)* &
!         barm(:ncol)*(1._r8+tni(:ncol)/3330._r8)*hei(:ncol)/(re**2*sqrt(2._r8*pi*grav)*rmass_he)

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       do i = 1,ncol
          tn = phys_state(lchnk)%t(i,1) ! top layer temperature
          he_mmr = phys_state(lchnk)%q(i,1, he_cnst_ndx)
          flx_phys(i,lchnk) = -4._r8*p0*sqrt((gask*tn/(rmass_he*grav))**3) ! ...
       end do
       call outfld('HEFLUX_TST1', flx_phys(:ncol,lchnk), ncol, lchnk)
    end do

    call regrid_phys2lonlat(flx_phys,flx_lonlat)

    call regrid_lonlat2phys(flx_lonlat,tmp_phys)

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       call outfld('HEFLUX_TST2', tmp_phys(:ncol,lchnk), ncol, lchnk)
    end do

  end subroutine helium_ubc_calc

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine read_coefs_file()

    use netcdf,only:nf90_open,nf90_inq_dimid,nf90_inquire_dimension, &
         nf90_inq_varid,nf90_get_var,nf90_close,nf90_nowrite,nf90_noerr

    character(len=*), parameter :: prefix = 'helium_ubc_mod->read_ceofs_file : '
    integer :: ncid, stat
    integer :: dimid, varid, length

    stat = nf90_open(trim(helium_ubc_coefs),nf90_nowrite,ncid)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_open',stat)

    stat = nf90_inq_dimid(ncid,'lat1',dimid)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=nlat_he)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_inquire_dimension',stat)

    stat = nf90_inq_dimid(ncid,'lat2',dimid)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=length)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_inquire_dimension',stat)

    if (length /= nlat_he-1) then
      call endrun(prefix//'Helium coefficient file dimension lat2 does not conform')
    end if

    stat = nf90_inq_dimid(ncid,'lat3',dimid)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_inq_dimid',stat)

    stat = nf90_inquire_dimension(ncid,dimid,len=length)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_inquire_dimension',stat)

    if (length /= nlat_he) then
      call endrun(prefix//'Helium coefficient file dimension lat3 does not conform')
    end if

    allocate(pmn(nlat_he,0:nlat_he-2,0:nlat_he-1))
    allocate(zmn(nlat_he,0:nlat_he-2,0:nlat_he-1))

    stat = nf90_inq_varid(ncid,'pmn',varid)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,pmn)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_get_var',stat)

    stat = nf90_inq_varid(ncid,'zmn',varid)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_inq_varid',stat)

    stat = nf90_get_var(ncid,varid,zmn)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_get_var',stat)

    stat = nf90_close(ncid)
    if (stat /= nf90_noerr) call handle_nc_error('nf90_close',stat)

  contains

    subroutine handle_nc_error(funcname,ncerr)

      use netcdf, only: nf90_strerror

      character(len=*),intent(in) :: funcname
      integer,intent(in) :: ncerr

      character(len=cl) :: errstr

      write(errstr,"('NetCDF error encountered: ',a,', when calling ',a)") &
           trim(nf90_strerror(ncerr)), funcname

      call endrun(prefix//trim(errstr))

    endsubroutine handle_nc_error

  end subroutine read_coefs_file


end module helium_ubc_mod
