module helium_ubc_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_kind_mod, only: cl => shr_kind_cl
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils, only: masterproc, host_mpicom=>mpicom, masterprocid, mpi_character, mpi_success
  use ppgrid, only: begchunk, endchunk, pcols, pver, pverp
  use physics_types, only: physics_state

  use esmf_phys_mesh_mod, only: esmf_phys_mesh_init
  use esmf_lonlat_grid_mod, only: esmf_lonlat_grid
  use esmf_lonlatphys_regrid_mod, only: esmf_lonlatphys_regrid_init, regrid_lonlat2phys, regrid_phys2lonlat
  use helium_zonal_fft_mod, only: helium_zonal_fft_init, helium_zonal_fft_forward, helium_zonal_fft_backward

  implicit none

  private
  public :: helium_ubc_readnl
  public :: helium_ubc_init
  public :: helium_ubc_calc
  public :: helium_ubc_fluxes

  real(r8), protected, pointer :: helium_ubc_fluxes(:,:) => null()

  character(len=cl) :: helium_ubc_coefs = 'NONE'

  integer :: mytid = -huge(1)

  real(r8), parameter :: gask = 8.314e7_r8      ! gas constant (erg / K / mole)
  real(r8), parameter :: rmass_he = 4._r8       ! (grams/mole)
  real(r8), parameter :: grav = 870._r8         ! cm/sec/sec
  real(r8), parameter :: pi = 4._r8*atan(1.0_r8)!
  real(r8), parameter :: re = 6.37122e8_r8      ! earth radius (cm)
  real(r8), parameter :: p0 = 5.e-4_r8          ! standard pressure Pa ?

  integer :: he_cnst_ndx = -1
  integer :: nlat_he = -1
  integer :: nlon_he = -1

  real(r8),dimension(:,:,:),allocatable :: pmn,zmn

  logical :: he_ubc_active = .false.

  class(esmf_lonlat_grid), pointer :: he_grid => null()

  integer :: beglat = 0
  integer :: endlat = 0
  integer :: beglon = 0
  integer :: endlon = 0

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

    he_ubc_active = (helium_ubc_coefs /= 'NONE') .and. (len_trim(helium_ubc_coefs)>0)

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

    integer :: ierr, npes, npes_he
    character(len=*), parameter :: prefix = 'helium_ubc_init: '

    if (.not.he_ubc_active) return

    call read_coefs_file()

    call mpi_comm_size(host_mpicom, npes, ierr)

    npes_he = min( npes, nlat_he/2 )

    ! write to atm log
    if (masterproc) then
       write(iulog,*) prefix, 'helium_ubc_nlats: ', nlat_he
       write(iulog,*) prefix, 'helium_ubc  npes: ', npes_he
       write(iulog,*) prefix, 'host model mpicom: ',host_mpicom
    end if

    call esmf_phys_mesh_init()
    he_grid => esmf_lonlat_grid(nlat_he, npes_he)

    call esmf_lonlatphys_regrid_init(he_grid)

    call helium_zonal_fft_init(he_grid)

    nlon_he = he_grid%nlon
    beglon = he_grid%lon_beg
    endlon = he_grid%lon_end
    beglat = he_grid%lat_beg
    endlat = he_grid%lat_end

    allocate(helium_ubc_fluxes(pcols,begchunk:endchunk))
    helium_ubc_fluxes = 0._r8

    call cnst_get_ind( 'HE', he_cnst_ndx )

    call addfld('HE_UBC_FLXi', horiz_only,  'A', '#/cm4/sec', 'Initial phyics grid He flux / Re^2 ' )
    call addfld('HE_UBC_FLUX', horiz_only,  'A', '#/cm4/sec', 'Upper boundary He flux / Re^2' )

  end subroutine helium_ubc_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine helium_ubc_calc(phys_state)
    use cam_history, only: outfld
    use mo_mean_mass, only: set_mean_mass
    use chemistry, only: imozart
    use ref_pres, only: ptop_ref

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    integer,parameter :: nmax = 35
    integer,parameter :: truncdeg = 8

    real(r8) :: flx_phys(pcols,begchunk:endchunk)

    real(r8) :: flx_lonlat(beglon:endlon,beglat:endlat)
    real(r8) :: tn, he_mmr
    real(r8) :: barm(pcols,pver) ! mean molecular weight (g/mole)

    integer :: lchnk, ncol, i, j, m,n
    complex(r8) :: zin (nlon_he, beglat:endlat)
    complex(r8) :: zout(nlon_he, beglat:endlat)

    ! Fourier coefficients ordered in real/image pairs
    real(kind=r8),dimension(nlat_he,nlon_he+2) :: fx_f

    real(kind=r8),dimension(0:nmax-1,0:nmax) :: amn ! a(m,n) spectral coefficient
    real(kind=r8),dimension(1:nmax-1,0:nmax) :: bmn ! b(m,n) spectral coefficient

    real(r8) :: ptop ! top of model pressure in cgs units

    if (.not.he_ubc_active) return

    ptop = ptop_ref * 10._r8 ! Pa --> Ba (dyne/cm2 or g/cm2/sec2)

    associate( nlon_he=>he_grid%nlon, &
               lon_beg=>he_grid%lon_beg, &
               lon_end=>he_grid%lon_end, &
               lat_beg=>he_grid%lat_beg, &
               lat_end=>he_grid%lat_end )

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       call set_mean_mass( ncol, lchnk, phys_state(lchnk)%q(:,:,imozart:), barm )
       do i = 1,ncol
          tn = phys_state(lchnk)%t(i,1) ! top layer temperature ! nmbr dens flux / Re^2 (#/cm4/sec)
          he_mmr = phys_state(lchnk)%q(i,1, he_cnst_ndx)
          flx_phys(i,lchnk) = -4._r8*ptop*sqrt((gask*tn/(rmass_he*grav))**3)*barm(i,1) & ! p0 --> p_top interface
               *(1._r8+tn/3330._r8)*he_mmr/(re**2*sqrt(2._r8*pi*grav)*rmass_he)
       end do
       call outfld('HE_UBC_FLXi', flx_phys(:ncol,lchnk), ncol, lchnk)
    end do

    call regrid_phys2lonlat(flx_phys,flx_lonlat)

! Forward transform from gridpoint to Fourier space:
    zout = helium_zonal_fft_forward(flx_lonlat)

    fx_f(:,:) = 0._r8

! reorder complex Fourier coefficients to real/image pairs
    do j = lat_beg,lat_end
       do i = 1,nlon_he/2+1
          fx_f(j,i*2-1) = zout(i,j)%re ! real(zout(i),kind=r8)
          fx_f(j,i*2)   = zout(i,j)%im ! aimag(zout(i))
       enddo
    enddo

! FFTW leaves an extra N after forward Fourier transform
    do j = lat_beg,lat_end
       do i = 1,nlon_he+2
          fx_f(j,i) = fx_f(j,i)/nlon_he
       enddo
    enddo

    ! global fluxes are needed for A and B accumulations
    fx_f(:,:) = allgather_flx(fx_f)

! order of coefficients:
! A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
! where B(0)=B(N/2)=0

! Fit (now on global grid)
    amn = 0
    bmn = 0
    do n = 0,truncdeg ! for no truncation, sum to nmax-1
       do j = 1,nlat_he
          amn(0,n) = amn(0,n)+zmn(j,0,n)*fx_f(j,1)
       enddo
       do m = 1,n
          do j = 1,nlat_he
             amn(m,n) = amn(m,n)+zmn(j,m,n)*fx_f(j,2*m+1)
             bmn(m,n) = bmn(m,n)+zmn(j,m,n)*fx_f(j,2*m+2)
          enddo
       enddo
    enddo

! Synthesis
    fx_f = 0
    do n = 0,truncdeg ! for no truncation, sum to nmax-1
       do j = lat_beg,lat_end
          fx_f(j,1) = fx_f(j,1)-n*(n+1)*amn(0,n)*pmn(j,0,n) ! A(0)
       enddo
       do m = 1,n
          do j = lat_beg,lat_end
             fx_f(j,2*m+1) = fx_f(j,2*m+1)-n*(n+1)*amn(m,n)*pmn(j,m,n) ! A(1),A(2),...
             fx_f(j,2*m+2) = fx_f(j,2*m+2)-n*(n+1)*bmn(m,n)*pmn(j,m,n) ! B(1),B(2),...
          enddo
       enddo
    enddo

! reconstruct complex Fourier coefficients from real/image pairs
    do j = lat_beg,lat_end
       zin(1,j) = fx_f(j,1)
       zin(nlon_he/2+1,j) = fx_f(j,nlon_he+1)
       do i = 2,nlon_he/2
          zin(i,j)           = cmplx(fx_f(j,i*2-1), fx_f(j,i*2),kind=r8)
          zin(nlon_he+2-i,j) = cmplx(fx_f(j,i*2-1),-fx_f(j,i*2),kind=r8)
       enddo
    end do

! Inverse transform from Fourier space to gridpoint:
    zout = helium_zonal_fft_backward(zin)

    do j = lat_beg,lat_end
       do i = lon_beg,lon_end
          flx_lonlat(i,j) = zout(i,j)%re
       end do
    end do

! regrid lon-lat to phyics column grid
    call regrid_lonlat2phys( flx_lonlat, helium_ubc_fluxes )

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       call outfld('HE_UBC_FLUX', helium_ubc_fluxes(:ncol,lchnk), ncol, lchnk)
    end do

    end associate

  contains

    function allgather_flx(fx_in) result(fx_glb)
      use mpi, only: MPI_REAL8, MPI_SUCCESS, MPI_SUM

      real(r8), intent(in) :: fx_in(nlat_he,nlon_he+2)

      real(r8) :: fx_glb(nlat_he,nlon_he+2)
      real(r8) :: sndbf(nlat_he,nlon_he+2)
      real(r8) :: rcvbf(nlat_he,nlon_he+2)
      integer :: len, rc

      character(len=*),parameter :: subname = 'helium_ubc_calc.allgather_flx: '

      len = nlat_he*(nlon_he+2)
      rcvbf(:,:) = 0._r8
      sndbf(:,:) = 0._r8
      sndbf(beglat:endlat,:) = fx_f(beglat:endlat,:)
      call mpi_allreduce( sndbf, rcvbf, len, MPI_REAL8, MPI_SUM, he_grid%merid_comm, rc )
      if ( rc /= MPI_SUCCESS ) then
         call endrun(subname//'mpi_allreduce failed')
      end if

      fx_glb(:,:) = rcvbf(:,:)

    end function allgather_flx

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
