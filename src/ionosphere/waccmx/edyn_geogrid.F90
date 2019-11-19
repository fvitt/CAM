module edyn_geogrid
!
! Global geographic grid.
! See sub set_geogrid (edyn_init.F90)
!
  use shr_kind_mod, only: r8 => shr_kind_r8 ! 8-byte reals
  use infnan,       only: nan, assignment(=)
  use cam_logfile,  only: iulog
  use cam_abortutils, only: endrun

  implicit none
  private
  save

  integer, public, protected :: & ! dimensions
    nlat,    & ! number of latitudes
    nlon,    & ! number of longitudes
    nlev,    & ! number of midpoint levels
    nilev      ! number of interface latitudes

  real(r8), public, protected, allocatable , dimension(:) :: & ! coordinate vars
    glat,    & ! latitude coordinates (degrees)
    glon,    & ! longitude coordinates (degrees)
    ylatg,   & ! latitudes (radians)
    ylong,   & ! longitudes (radians)
    zlev,    & ! midpoint vertical coordinates
    zilev      ! interface vertical coordinates

  real(r8), public, allocatable, protected :: cs(:)   ! cos(phi) (0:nlat+1)
  real(r8), public, allocatable            :: expz(:) ! exp(-zp)
  real(r8), public, allocatable            :: zp(:)   ! log pressure (as in tiegcm lev(nlev))

  integer, public, protected :: & ! model independent (set by sub get_geogrid)
    nlonp1,  & ! nlon+1
    nlonp2,  & ! nlon+2
    nlatp1     ! nlat+1

  real(r8) :: dlatg, dlong
  public   :: dlong
  real(r8), public, protected :: dphi
  real(r8), public, protected :: dlamda
!
! Using p0 in microbars, as in TIEGCM.
  real(r8), parameter, public :: p0 = 5.0e-4_r8  ! standard pressure (microbars)

  integer, public, protected :: & ! model dependent (set by subs read_tgcm, read_waccm)
    jspole,  & ! latitude index to geographic south pole
    jnpole     ! latitude index to geographic north pole
!
! lev_sequence is a string indicating ordering of the vertical
! coordinates lev and ilev, and of the field arrays along the
! vertical dimension. lev_sequence can have 1 of 2 values:
!
!  'bottom2top' means lev(1) is the bottom boundary, lev(nlev) is the top boundary
!  'top2bottom' means lev(1) is the top boundary, lev(nlev) is the bottom boundary
!
! For example, TIMEGCM history files are bottom2top, whereas
! WACCM files are top2bottom. The edynamo code assumes bottom2top,
! so WACCM input fields are reversed to be bottom2top for the edynamo
! calculations, then reversed back to the native WACCM sequence
! (top2bottom) before writing to the edynamo output file.
!
  character(len=10) :: lev_sequence
!
! lon_sequence is a string indicating ordering of the longitude
! coordinate lon, and of the field arrays along this dimension.
! lon_sequece can have 1 of 2 values:
!
!   '-180to180' means lon(1) is -180 deg west longitude, lon(nlon) is +180 east
!   'zeroto360' means lon(1) is 0 deg west longitude, lon(nlon) is 360 deg east
!
! Note that TIMEGCM convention is '-180to180' and WACCM convention is 'zeroto360'
! (this is treating similarly to lev_sequence above)
!
  character(len=9) :: lon_sequence

  ! distribute_geogrid lays out the oplus grid
!!XXgoldyXX: v debug only
!  public :: distribute_geo_grid
!!XXgoldyXX: ^ debug only
  ! set_geogrid sets up a distributed finite-volume lat/lon grid
  public :: set_geogrid

  logical :: debug = .false. ! set true for prints to stdout at each call

contains

!!XXgoldyXX: v debug only
   ! !-----------------------------------------------------------------------
   ! subroutine distribute_geo_grid(mpicom, oplus_npes, oplus_nlon, oplus_nlat, &
   !            lon0, lon1, lat0, lat1, lev0, lev1, ntaski, ntaskj, glon, glat)
   !    ! distribute_geo_grid imitates the action of the FV dycore's
   !    ! decomposition routine. This allows us to separete the geogrid from
   !    ! the FV dycore

   !    ! Dummy arguments
   !    integer, intent(in) :: mpicom
   !    integer, intent(in) :: oplus_npes
   !    integer, intent(in) :: oplus_nlon
   !    integer, intent(in) :: oplus_nat
   !    integer, intent(in) :: lat0, lat1       ! first and last latitude  indices
   !    integer, intent(in) :: lon0, lon1       ! first and last longitude indices
   !    integer, intent(in) :: lev0, lev1       ! first and last pressure indices
   !    integer, intent(in) :: ntaski  ! number of MPI tasks in lon dimensions
   !    integer, intent(in) :: ntaskj  ! number of MPI tasks in lat dimensions
   !    real(r8), allocatable, intent(out) :: glon(:) ! global geo-graphic longitudes (degrees)
   !    real(r8), allocatable, intent(out) :: glat(:) ! global geo-graphic latitudes (degrees)

   !    ! Local variables
   !    integer :: iam, ierr       ! MPI task number

   !    if (iam < oplus_npes, ionos_dynamo_npes)) then

   !          allocate(glon(oplus_nlon))
   !          allocate(glat(oplus_nlat))

   !          lon0 = grid%ifirstxy ; lon1 = grid%ilastxy
   !          lat0 = grid%jfirstxy ; lat1 = grid%jlastxy
   !          lev0 = 1             ; lev1 = grid%km
   !          ntaski = grid%nprxy_x
   !          ntaskj = grid%nprxy_y
!!XXgoldyXX: ^ debug only

   !-----------------------------------------------------------------------
   subroutine set_geogrid(nlon_g, nlat_g, nlev_in, npes_in, atm_mpicom, pres_mid_in, pres_edge_in, min_lat_pe_in)
      use shr_const_mod,  only: pi => shr_const_pi
      use edyn_params,    only: kbotdyn, pbotdyn
      use edyn_mpi,       only: mp_distribute_geo
      use spmd_utils,     only: masterproc

      ! Dummy Args
      integer,            intent(in) :: nlon_g        ! Global num longitudes
      integer,            intent(in) :: nlat_g        ! Global num latitudes
      integer,            intent(in) :: nlev_in       ! Num levels
      integer,            intent(in) :: npes_in
      integer,            intent(in) :: atm_mpicom
      real(r8),           intent(in) :: pres_mid_in(:)
      real(r8),           intent(in) :: pres_edge_in(:)
      integer,  optional, intent(in) :: min_lat_pe_in ! Min # lats / PE
      !
      ! Local:
      integer                        :: latind, lonind, js, k
      integer                        :: lon_beg, lon_end, lat_beg, lat_end
      integer                        :: lons_per_task, lats_per_task
      integer                        :: lons_overflow, lats_overflow
      integer                        :: ntasks_lat, ntasks_lon, ntasks_atm
      integer                        :: task_cnt, i,j
      integer                        :: minlats_per_pe, npes
      integer                        :: geogrid_mpicom
      integer                        :: iam, color
      integer                        :: ierr
      real(r8)                       :: phi
      real(r8)                       :: delta ! Coordinate spacing
      real(r8), parameter            :: eps = 1.e-6_r8

      real(r8)                       :: pmid(nlev_in)

      nlon = nlon_g
      nlat = nlat_g
      nlev = nlev_in
      npes = npes_in

      nilev = nlev+1

      nlonp1 = nlon + 1
      nlonp2 = nlon + 2
      nlatp1 = nlat + 1

      jspole = 1
      jnpole = nlat

      if (present(min_lat_pe_in)) then
         minlats_per_pe = min_lat_pe_in
      else
         minlats_per_pe = 3
      end if

      dphi   = pi / real(nlat,r8)
      dlamda = 2._r8*pi / real(nlon,r8)

      !
      ! Allocate coordinate variables:
      !
      allocate(glon(nlon))
      allocate(glat(nlat))
      !
      ! Create a finite-volume coordinate grid (in degrees)
      !
      delta = 360.0_r8 / real(nlon, r8)
      do lonind = 1, nlon
         glon(lonind) = -180.0_r8 + ((lonind - 1) * delta)
      end do
      delta = 180.0_r8 / real((nlat - 1), r8)
      ! Set the poles exactly (they might be checked later)
      glat(1) = -90.0_r8
      glat(nlat) = 90.0_r8
      do latind = 2, nlat - 1
         glat(latind) = -90.0_r8 + ((latind - 1) * delta)
      end do

      allocate(zlev(nlev))
      allocate(zilev(nilev))
      !
      ! zp and expz are not set until oplus is called from dpie_coupling.
      allocate(zp(nlev))      ! log pressure (as in TIEGCM)
      allocate(expz(nlev))    ! exp(-zp)
      zp = nan
      expz = nan
      !
      ! Hybrid-sigma levels from ref_pres module:
      !
      zlev(:nlev)  = pres_mid_in(:)  ! midpoints vertical coord (top down)
      zilev(:nilev) = pres_edge_in(:nilev)  ! interfaces vertical coord

      ! do bottom up search for kbotdyn
      pmid(:nlev) = zlev(nlev:1:-1)
      kloop: do k = 1, nlev
         if ( pmid(k) <= pbotdyn) then
            kbotdyn = k
            exit kloop
         end if
      end do kloop
      if ( kbotdyn < 1 ) then
         call endrun('set_geogrid: kbotdyn is not set')
      endif
      if (debug) then
         write(iulog,"('set_geogrid: kbotdyn=',i4,' pmid(kbotdyn)=',es12.4)") kbotdyn,pmid(kbotdyn)
      endif

      !
      ! Setup a decomposition for the geogrid
      !
      ! First, try using a 1-D latitude decomposition
      call MPI_comm_size(atm_mpicom, ntasks_atm, ierr)
      if (npes > ntasks_atm) then
         call endrun('set_geogrid: npes too large for ntasks_atm')
      else if (npes <= 0) then
         npes = ntasks_atm ! Use all tasks by default
      end if
      call MPI_comm_rank(atm_mpicom, iam, ierr)
      if (nlat / npes < minlats_per_pe) then
         ! We need to work up a 2-D decomposition
         ntasks_lat = floor(real(nlat, r8) / real(minlats_per_pe, r8))
         ntasks_lon = 1
         do lonind = 2, npes / ntasks_lat
            if (nlon / lonind < minlats_per_pe) then
               exit
            else
               ntasks_lon = lonind
            end if
            if (ntasks_lon * ntasks_lat >= npes) then
               exit
            end if
         end do
      else
         ntasks_lat = npes
         ntasks_lon = 1
      end if

      ! Create the geogrid communicator
      color = iam / (ntasks_lat * ntasks_lon)
      call MPI_comm_split(atm_mpicom, color, iam, geogrid_mpicom, ierr)
      ! A quick sanity check
      call MPI_comm_size(geogrid_mpicom, ntasks_atm, ierr)
      color = ntasks_lat * ntasks_lon
      if ((iam < color) .and. (ntasks_atm /= color)) then
         call endrun('set_geogrid: Incorrect size for geogrid_mpicom')
      end if

      ! Now, figure the starting and ending coordinates
      lons_per_task = nlon / ntasks_lon
      lons_overflow = MOD(nlon, ntasks_lon)
      lats_per_task = nlat / ntasks_lat
      lats_overflow = MOD(nlat, ntasks_lat)
      lon_beg = 1
      lon_end = 0
      lat_beg = 1
      lat_end = 0
      task_cnt= 0
      jloop: do j = 0,ntasks_lat-1
         lat_beg = lat_end + 1
         lat_end = lat_beg + lats_per_task - 1
         if (j<lats_overflow) then
            lat_end = lat_end + 1
         end if
         lon_end = 0
         do i = 0,ntasks_lon-1
            lon_beg = lon_end + 1
            lon_end = lon_beg + lons_per_task - 1
            if (i<lons_overflow) then
               lon_end = lon_end + 1
            end if
            task_cnt = task_cnt+1
            if (task_cnt>iam) exit jloop
         end do
      enddo jloop

      if (masterproc) then
         write(iulog,'(a,)') '****************************************************** '
         write(iulog,'(a,3i)') '*** set_geogrid: npes, ntasks_lon, ntasks_lat : ', npes, ntasks_lon, ntasks_lat
         write(iulog,'(a,)') '****************************************************** '
      endif
      
      call mp_distribute_geo(lon_beg, lon_end, lat_beg, lat_end, 1, nlev, ntasks_lon, ntasks_lat)
      !!XXgoldyXX: Do we need the geogrid_mpicom for anything?
      call MPI_comm_free(geogrid_mpicom, ierr)
      !
      ! Set horizontal geographic grid in radians (for apex code):
      !
      allocate(ylatg(nlat))   ! waccm grid includes poles
      allocate(ylong(nlonp1)) ! single periodic point
      dlatg = pi / real(nlat,r8)
      dlong = 2._r8*pi / real(nlon,r8)
      ylatg(1)    = -pi/2._r8+eps ! south pole
      ylatg(nlat) =  pi/2._r8-eps ! north pole
      do latind = 2, nlat-1
         ylatg(latind) = -0.5_r8*(pi-dlatg)+real(latind-1,r8)*dlatg
      end do
      do lonind = 1, nlonp1
         ylong(lonind) = -pi+real(lonind-1,r8)*dlong
      end do
      !
      ! Calculate cosine of latitude
      !
      allocate(cs(0:nlat+1))
      js = -(nlat/2)
      do latind = 1, nlat
         phi = (latind + js - .5_r8) * dphi
         cs(latind) = cos(phi)
      end do
      cs(0) = -cs(1)
      cs(nlat+1) = -cs(nlat)

   end subroutine set_geogrid

end module edyn_geogrid
