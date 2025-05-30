module vxb
!
! Horizontally transport the O+ ion, adapted for WACCM-X from TIEGCM. 
! Input O+ is received from WACCM physics/chemistry, transported O+ 
! (op_out and opnm_out) are passed back to chemistry.
!
! B. Foster (foster@ucar.edu), May, 2015.
!
  use shr_kind_mod   ,only: r8 => shr_kind_r8
  use cam_abortutils ,only: endrun
  use cam_logfile    ,only: iulog
  use savefield_waccm,only: savefld_waccm, savefld_waccm_switch     ! save field to waccm history
  use edyn_geogrid   ,only: dphi,dlamda,cs,zp,expz,p0 !, nlon, nlat, nlev
  use getapex        ,only: bx,by,bz,bmod2       ! (0:nlonp1,jspole-1:jnpole+1)
  use edyn_params    ,only: re
  use time_manager   ,only: get_step_size,is_first_step,is_first_restart_step
  use edyn_mpi       ,only: array_ptr_type
  use shr_const_mod  ,only: shr_const_g         ! gravitational constant (m/s^2)
  use spmd_utils     ,only: masterproc

  implicit none
  private
  public :: vxb_xport,  vxb_init
  public :: kbot
  public :: ktop

  real(r8) :: pi,rtd
!
! Constants in CGS:
!
  real(r8),parameter :: boltz = 1.38E-16_r8 ! boltzman's constant (erg/kelvin)
  real(r8),parameter :: gask  = 8.314e7_r8  ! gas constant (erg/mol)
!
! Collision factor (tuneable) (see also local colfac in iondrag.F90)
! FIX: Collision factor colfac is set locally in iondrag.F90 and here.
!      It should be in one location and shared between ionosphere and
!      dpie_coupling.
!
    real(r8),parameter :: colfac = 1.5_r8   ! see also iondrag.F90
!
! Reciprocal of molecular mass (multiply is cheaper than divide)
    real(r8),parameter :: rmassinv_o2=1._r8/32._r8, rmassinv_o1=1._r8/16._r8, &
                          rmassinv_n2=1._r8/28._r8
!    real(r8),parameter :: rmass_fep=24.3_r8

    real(r8) :: dzp                       ! delta zp (typically 0.5 from kbot to top)
    real(r8) :: grav_cm                   ! gravitational constant (cm/s^2)
    integer, protected :: kbot = -999     ! k-index corresponding to ~pbot
    real(r8),parameter :: pbot = 10._r8 ! Pa -- bottom of metal ions transport (near 80 km) need check
    integer, protected :: ktop = -999     ! k-index corresponding to ~ptop
    real(r8),parameter :: ptop = 1.0e-7_r8 ! Pa -- top of metal ions transport (near 180 km) need check
    integer, protected :: kminor = -999     ! k-index corresponding to ~ptop
    real(r8),parameter :: pminor = 1.0e-3_r8 ! Pa -- top of metal ions transport (near 180 km) need check
!
! The shapiro constant .03 is used for spatial smoothing of oplus,
! (shapiro is tuneable, and maybe should be a function of timestep size).
! dtsmooth and dtsmooth_div2 are used in the time smoothing.
! To turn off all smoothing here, set shapiro=0. and dtsmooth = 1. 
!

    real(r8),parameter ::    & 
      dtsmooth = 0.95_r8,    &                ! for time smoother
      dtsmooth_div2 = 0.5_r8*(1._r8-dtsmooth)
 
    real(r8) :: adiff_limiter
    real(r8) :: shapiro_const
    logical  :: enforce_floor
    logical, parameter :: debug = .false.

  contains

!-----------------------------------------------------------------------
  subroutine vxb_init( adiff_limiter_in, shapiro_const_in, enforce_floor_in )

    use cam_history,  only : addfld, horiz_only
    use vxbfilter_module,only : vxbfilter_init

    real(r8), intent(in) :: adiff_limiter_in
    real(r8), intent(in) :: shapiro_const_in
    logical , intent(in) :: enforce_floor_in

    shapiro_const = shapiro_const_in
    enforce_floor = enforce_floor_in
    adiff_limiter = adiff_limiter_in
    
    call vxbfilter_init

    !
    ! Save fields from oplus module:
    !

  end subroutine vxb_init

!-----------------------------------------------------------------------
  subroutine vxb_xport(tn,te,ti,ne,nop,un,vn,om,zg,o2,o1,n2,op_in,opnm_in,xi, &
                         mbar,ui,vi,wi,eeex,eeey,eeez,pmid,op_out,opnm_out,rmass_fep, &
                         i0,i1,j0,j1,nspltop,ispltop )
!
! All input fields from dpie_coupling are in "TIEGCM" format, i.e., 
! longitude (-180->180), vertical (bot2top), and units (CGS).
!
    use edyn_mpi,only:  mp_geo_halos,mp_pole_halos,setpoles
    use edyn_geogrid,only : glat, nlat, nlev
!
! Transport O+ ion.
! March-May, 2015 B.Foster: Adapted from TIEGCM (oplus.F) for WACCM-X.
!
! Notes:
! - waccmx_opt='ionosphere' must be set in user_nl_cam for te,ti inputs to have values
!
! Args:
!
    integer,intent(in) :: &
      i0,                 & ! grid%ifirstxy
      i1,                 & ! grid%ilastxy
      j0,                 & ! grid%jfirstxy
      j1                    ! grid%jlastxy
    integer,intent(in) :: nspltop,ispltop
!
! Input fields without halo points (lon +/-180, vertical bot2top, CGS units):
!
    real(r8),intent(in) :: tn   (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral temperature (deg K)
    real(r8),intent(in) :: te   (nlev,i0-2:i1+2,j0-2:j1+2) ! electron temperature (deg K)
    real(r8),intent(in) :: ti   (nlev,i0-2:i1+2,j0-2:j1+2) ! ion temperature (deg K)
    real(r8),intent(in) :: ne   (nlev,i0-2:i1+2,j0-2:j1+2) ! electron number density ion (/cm3)
    real(r8),intent(in) :: nop   (nlev,i0-2:i1+2,j0-2:j1+2) ! op number density ion (/cm3)
    real(r8),intent(in) :: un   (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral zonal wind (cm/s)
    real(r8),intent(in) :: vn   (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral meridional wind (cm/s)
    real(r8),intent(in) :: om   (nlev,i0-2:i1+2,j0-2:j1+2) ! omega (1/s)
    real(r8),intent(in) :: o2   (nlev,i0-2:i1+2,j0-2:j1+2) ! o2 (mmr)
    real(r8),intent(in) :: o1   (nlev,i0-2:i1+2,j0-2:j1+2) ! o (mmr)
    real(r8),intent(in) :: n2   (nlev,i0-2:i1+2,j0-2:j1+2) ! n2 (mmr)
    real(r8),intent(in) :: mbar (nlev,i0-2:i1+2,j0-2:j1+2) ! mean molecular weight

    real(r8),intent(in) :: op_in(nlev,i0:i1,j0:j1) ! O+ density (cm^3)
    real(r8),intent(in) :: opnm_in(nlev,i0:i1,j0:j1) ! O+ density (cm^3) at time-1
    real(r8),intent(in) :: rmass_fep
    real(r8),intent(in) :: ui(nlev,i0:i1,j0:j1)   ! zonal ion drift
    real(r8),intent(in) :: vi(nlev,i0:i1,j0:j1)   ! meridional ion drift
    real(r8),intent(in) :: wi(nlev,i0:i1,j0:j1)   ! vertical ion drift
    real(r8),intent(in) :: eeex(nlev,i0:i1,j0:j1)   ! zonal ion drift
    real(r8),intent(in) :: eeey(nlev,i0:i1,j0:j1)   ! meridional ion drift
    real(r8),intent(in) :: eeez(nlev,i0:i1,j0:j1)   ! vertical ion drift
    real(r8),intent(in) :: xi(nlev,i0:i1,j0:j1) ! Jianfei Wu added xi/(1+xi^2)
    real(r8),intent(in) :: zg   (nlev,i0:i1,j0:j1) ! geopotential height (cm)
!
! Ion drifts from edynamo (also in tiegcm-format):
!
    real(r8),intent(in) :: pmid(nlev)             ! pressure at midpoints (Pa)
!
! Output:
!
    real(r8),intent(out) ::       &
      op_out  (nlev,i0:i1,j0:j1), & ! O+ output
      opnm_out(nlev,i0:i1,j0:j1)    ! O+ output at time n-1
!
! Local:
!
    integer :: i,j,k,lat,jm1,jp1,jm2,jp2,lat0,lat1
    real(r8),dimension(i0:i1,j0:j1) :: &
      dvb         ! divergence of B-field
!
! Local inputs with added halo points in lat,lon:
! 
    real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: op, opnm,xiop,xiopnm, &
        xxiop,xxiopnm,xxxiop,xxxiopnm,xxxxiop,xxxxiopnm,eeeex,eeeey,eeeez,uuii,vvii,wwii,uui,vvi,wwi,wn,&
        ri,rri,rrri,rrrri,bdotu


    real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: &
      tr          ,&   ! Reduced temperature (.5*(tn+ti))
      tp          ,&   ! Plasma temperature N(O+)*(te+ti)
      tpp          ,&   ! Plasma temperature N(O+)*(te+ti)
      dj          ,&   ! diffusion coefficients
      dji          ,&   ! diffusion coefficients
      bvel        ,&   ! bvel @ j   = (B.U)*N(O+)
      diffj       ,&   ! (D/(H*DZ)*2.*TP+M*G/R)*N(O+)
      diffji       ,&   ! (D/(H*DZ)*2.*TP+M*G/R)*N(O+)
      diffjj       ,&   ! (D/(H*DZ)*2.*TP+M*G/R)*N(O+)
      bdotdh_op   ,&   ! (b(h)*del(h))*phi
      bdotdh_opi   ,&   ! (b(h)*del(h))*phi
      ambvel   ,&   ! (b(h)*del(h))*phi
      bdotdh_opj  ,&   ! (b(h)*del(h))*phi
      bdotdh_opj1  ,&   ! (b(h)*del(h))*phi
      bdotdh_opj2  ,&   ! (b(h)*del(h))*phi
      bdotdh_diff ,&   ! (b(h)*del(h))*phi
      opnm_smooth      ! O+ at time-1, smoothed
    real(r8),dimension(i0-2:i1+2,j0-2:j1+2),target :: &
      abk        !upper boundry codition 


    real(r8),dimension(nlev,i0:i1,j0:j1) :: &
        rwii,uii,vii,wii,  &! metal ion velocity from subroutine vxbvel  Jianfei Wu
        hadudth,epsh,scrh, &
        nulh, mulh, flxh, diff, LORHOT, LNRHOT, RHOT, RHOTD, FABS, FSGN, TERM, TERP,aa,va,rva
    real(r8),dimension(nlev,i0:i1,j0-1:j1+1) :: hj ! scale height
    real(r8) :: gmr,dtime,dtx2,dtx2inv, ONE6TH,ONE3RD
    real(r8),dimension(nlev,i0:i1) :: &
      bdzdvb_op,   &
      hdz,         &
      tp1,         &
      tphdz0,      &
      tphdz1,      &
      djint,       &
      divbz,       &
      hdzmbz,      &
      hdzpbz
!
!
! Arguments for tridiagonal solver trsolv (no halos):
    real(r8),dimension(nlev,i0:i1,j0:j1) :: &
      explicit1,explicit2,&
      explicit,explicit_a1,explicit_b1,&
      explicit3,explicit_a2,explicit_b2,&
      explicit4,explicit_a3,explicit_b3,&
      explicit5,explicit_a4,explicit_b4,&
      explicit_a5,explicit_b5,&
      explicit6,explicit7
    real(r8),dimension(i0:i1) :: ubca, ubcb
    real(r8),parameter :: one=1._r8
    logical :: calltrsolv
!
! Pointers for multiple-field calls (e.g., mp_geo_halos)
    integer :: nfields
    real(r8),allocatable :: polesign(:)
    type(array_ptr_type),allocatable :: ptrs(:)

    real(r8) :: zpmid(nlev), opfloor
    real(r8),parameter :: opmin=3000.0_r8
    Real(r8),Parameter :: BIGNUM = Huge(1._r8)
!
! Execute:
!
    dtime = get_step_size() ! step size in seconds
    dtime = dtime / dble(nspltop)
    dtx2 = 2._r8*dtime
    dtx2inv = 1._r8/dtx2

    if ((is_first_step().or.is_first_restart_step()).and.ispltop==1) then
      if (masterproc) write(iulog,"('mgplus: shapiro=',es12.4,' dtsmooth=',es12.4,' dtsmooth_div2=',es12.4)") &
        shapiro_const,dtsmooth,dtsmooth_div2
      if (masterproc) write(iulog,"('mgplus: shr_const_g=',f8.3)") shr_const_g
    endif

    !
    ! zp,expz are declared in edyn_geogrid.F90, and allocated in sub 
    ! set_geogrid (edyn_init.F90). pmid was passed in here (bot2top)
    ! from dpie_coupling.
    !
    ! kbot is the k-index at the bottom of O+ transport calculations,
    !   corresponding to pressure pbot. 
    !
    if ((is_first_step().or.is_first_restart_step()).and.ispltop==1) then
       kloop: do k=1,nlev
          if ( pmid(k) <= pbot) then
             kbot = k
             exit kloop
          end if
       enddo kloop
       kloop2: do k=1,nlev
          if ( pmid(k) <= ptop) then
             ktop = k-1
             exit kloop2
          end if
       enddo kloop2
       kloop3: do k=1,nlev
          if ( pmid(k) <= pminor) then
             kminor = k-1
             exit kloop3
          end if
       enddo kloop3
       do k=1,nlev
          zp(k) = -log(pmid(k)*10._r8/p0)
          expz(k) = exp(-zp(k))
       enddo
       if (debug.and.masterproc) then
          write(iulog,"('mgplus: kbot=',i4,'ktop=',i4,' pmid(kbot)=',es12.4,' zp(kbot)=',es12.4)") &
               kbot,ktop,pmid(kbot),zp(kbot)
       endif
    endif

    if (kbot < 1.or.ktop < kbot) then
       call endrun('vxb_xport: kbot or ktop is not set')
    endif

    dzp = zp(nlev)-zp(nlev-1)  ! use top 2 levels (typically dzp=0.5)

    if (debug.and.masterproc) then
      write(iulog,"('mgplus: nlev=',i3,' zp (bot2top)   =',/,(6es12.3))") nlev,zp
      write(iulog,"('mgplus: nlev=',i3,' expz (bot2top) =',/,(6es12.3))") nlev,expz
      write(iulog,"('mgplus: nlev=',i3,' dzp  =',/,(6es12.3))") nlev,dzp
    endif
!
! Set subdomain blocks from input (composition is in mmr):
!
!$omp parallel do private(i, j, k)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          ri(k,i,j) = xi(k,i,j)/(1._r8+xi(k,i,j)**2)
          rrri(k,i,j) = xi(k,i,j)**2/(1._r8+xi(k,i,j)**2)
          rrrri(k,i,j) = 1._r8/(1._r8+xi(k,i,j)**2)
          !rrrri(k,i,j) = (xi(k,i,j)**2)/(1._r8+xi(k,i,j)**2)
          rri(k,i,j) = 1._r8/(1._r8+xi(k,i,j)**2)
          op(k,i,j)   = op_in(k,i,j)
          opnm(k,i,j) = opnm_in(k,i,j)
          xiop(k,i,j)   = op_in(k,i,j)*ri(k,i,j)
          xiopnm(k,i,j) = opnm_in(k,i,j)*ri(k,i,j)
          xxiop(k,i,j)   = op_in(k,i,j)*rri(k,i,j)
          xxiopnm(k,i,j) = opnm_in(k,i,j)*rri(k,i,j)
          xxxiop(k,i,j)   = op_in(k,i,j)*rrri(k,i,j)
          xxxiopnm(k,i,j) = opnm_in(k,i,j)*rrri(k,i,j)
          xxxxiop(k,i,j)   = op_in(k,i,j)*rrrri(k,i,j)
          xxxxiopnm(k,i,j) = opnm_in(k,i,j)*rrrri(k,i,j)
          uui(k,i,j)=ui(k,i,j)
          vvi(k,i,j)=vi(k,i,j)
          wwi(k,i,j)=wi(k,i,j)
          eeeex(k,i,j)=eeex(k,i,j)!*bmod2(i,j)
          eeeey(k,i,j)=eeey(k,i,j)!*bmod2(i,j)
          eeeez(k,i,j)=eeez(k,i,j)!*bmod2(i,j)
          
        enddo
      enddo
    enddo

!
! Define halo points on inputs:
! WACCM has global longitude values at the poles (j=1,j=nlev)
! (they are constant for most, except the winds.)
!
! Set two halo points in lat,lon:
!   real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: tn,te,etc.
!
    nfields = 20
    allocate(ptrs(nfields),polesign(nfields))

    ptrs(1)%ptr => op ; ptrs(2)%ptr => opnm
    ptrs(3)%ptr =>xiop ; ptrs(4)%ptr => xiopnm
    ptrs(5)%ptr =>xxiop ; ptrs(6)%ptr => xxiopnm
    ptrs(7)%ptr =>eeeex ; ptrs(8)%ptr => eeeey
    ptrs(9)%ptr =>eeeez 
    ptrs(10)%ptr =>uui ; ptrs(11)%ptr => vvi
    ptrs(12)%ptr =>wwi 
    ptrs(13)%ptr =>xxxiop ; ptrs(14)%ptr => xxxiopnm
    ptrs(15)%ptr =>xxxxiop ; ptrs(16)%ptr => xxxxiopnm
    ptrs(17)%ptr =>ri ; ptrs(18)%ptr => rri
    ptrs(19)%ptr =>rrri ; ptrs(20)%ptr => rrrri
    polesign = 1._r8
    polesign(7:8) = -1._r8
    polesign(10:11) = -1._r8
!
! mp_geo_halos first arg:
!     type(array_ptr_type) :: fmsub(nf) ! (lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
!
    call mp_geo_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields)
!
! Set latitude halo points over the poles (this does not change the poles).
! (the 2nd halo over the poles will not actually be used (assuming lat loops
!  are lat=2,nlat-1), because jp1,jm1 will be the pole itself, and jp2,jm2 
!  will be the first halo over the pole)
!
! mp_pole_halos first arg:
!   type(array_ptr_type) :: f(nf) ! (nlev,i0-2:i1+2,j0-2:j1+2)

    call mp_pole_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields,polesign)
    deallocate(ptrs,polesign)

!
! Use below to exclude the poles (lat=2,nlat-1) from latitude scans.
!
    lat0 = j0
    lat1 = j1
    if (j0 == 1)    lat0 = 2
    if (j1 == nlat) lat1 = nlat-1
!
! Save input fields to WACCM histories. Sub savefld_waccm_switch converts
! fields from tiegcm-format to waccm-format before saving to waccm histories.
!
!
! Initialize output op_out with input op at 1:kbot-1, to retain values from 
! bottom of column up to kbot. This routine will change (transport) these 
! outputs only from kbot to the top (nlev).
!
    op_out   = 0._r8 
    opnm_out = 0._r8 
    op_out  (1:kbot-1,i0:i1,j0:j1) = op  (1:kbot-1,i0:i1,j0:j1)
    opnm_out(1:kbot-1,i0:i1,j0:j1) = opnm(1:kbot-1,i0:i1,j0:j1)
    op_out  (ktop+1:nlev,i0:i1,j0:j1) = op  (ktop+1:nlev,i0:i1,j0:j1)
    opnm_out(ktop+1:nlev,i0:i1,j0:j1) = opnm(ktop+1:nlev,i0:i1,j0:j1)
!
!
! Divergence of B (mag field) is returned by divb in dvb(i0:i1,j0:j1)
!
    call divb(dvb,i0,i1,j0,j1)
! The solver will be called only if calltrsolv=true. It is sometimes
! set false when skipping parts of the code for debug purposes.
! 
    calltrsolv = .true.
    tr  = 0._r8
    tp  = 0._r8
    tpp  = 0._r8
    dj  = 0._r8
    dji  = 0._r8
    hj  = 0._r8
    bvel= 0._r8
    diffj = 0._r8
    diffji = 0._r8
    diffjj = 0._r8
    abk = 1._r8
    opnm_smooth = 0._r8
    grav_cm = shr_const_g * 100._r8 ! m/s^2 -> cm/s^2
    hadudth= 0._r8
    epsh=0._r8
    scrh=0._r8
    nulh=0._r8
    mulh=0._r8
    flxh=0._r8
    diff=0._r8
    aa=0._r8
    va=0._r8
    rva=0._r8
    bdotu       = 0._r8
    !!!!!
    ktop=ktop+1
    kbot=kbot-1
!
!----------------------- Begin first latitude scan ---------------------
    do lat=lat0,lat1
      jm2 = lat-2
      jm1 = lat-1
      jp1 = lat+1
      jp2 = lat+2
!
! as of April, 2015, TIEGCM incorrectly uses te+ti instead of tn+ti
! This has not been fixed in TIEGCM, because fixing it causes a tuning 
! problem (ask Hanli and Wenbin). For WACCM, it is correct as below.
!
!$omp parallel do private(i,k)
      do i=i0,i1
!
! Reduced temperature (tpj in tiegcm):
! 'OPLUS_TR' (has constants at poles)
!
        do k=kbot,ktop
          tr(k,i,jm1) = 0.5_r8*(tn(k,i,jm1)+ti(k,i,jm1))
          tr(k,i,lat) = 0.5_r8*(tn(k,i,lat)+ti(k,i,lat))
          tr(k,i,jp1) = 0.5_r8*(tn(k,i,jp1)+ti(k,i,jp1))
        enddo
      enddo ! i=i0,i1
!
! rrk returns ambipolar diffusion coefficients in d(jm1),dj(lat),djp1(jp1):
! 'OPLUS_DJ' (has constants at poles)
!
      call rrk(                                            &
        tn(kbot:ktop,i0:i1,jm1),ti(kbot:ktop,i0:i1,jm1),mbar(kbot:ktop,i0:i1,jm1), &
        o2(kbot:ktop,i0:i1,jm1),o1  (kbot:ktop,i0:i1,jm1), &
        n2(kbot:ktop,i0:i1,jm1),nop(kbot:ktop,i0:i1,jm1),op(kbot:ktop,i0:i1,jm1), &
        rmass_fep,tr(kbot:ktop,i0:i1,jm1), &
        dj(kbot:ktop,i0:i1,jm1),dji(kbot:ktop,i0:i1,jm1),i0,i1,kbot,ktop,kminor)

      call rrk(                                            &
        tn(kbot:ktop,i0:i1,lat),ti(kbot:ktop,i0:i1,lat),mbar(kbot:ktop,i0:i1,lat), &
        o2(kbot:ktop,i0:i1,lat),o1  (kbot:ktop,i0:i1,lat), &
        n2(kbot:ktop,i0:i1,lat),nop(kbot:ktop,i0:i1,lat),op(kbot:ktop,i0:i1,lat), &
        rmass_fep,tr  (kbot:ktop,i0:i1,lat),&
        dj(kbot:ktop,i0:i1,lat),dji(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,kminor)

      call rrk(                                            &
        tn(kbot:ktop,i0:i1,jp1),ti(kbot:ktop,i0:i1,jp1),mbar(kbot:ktop,i0:i1,jp1), &
        o2(kbot:ktop,i0:i1,jp1),o1  (kbot:ktop,i0:i1,jp1), &
        n2(kbot:ktop,i0:i1,jp1),nop(kbot:ktop,i0:i1,jp1),op(kbot:ktop,i0:i1,jp1), &
        rmass_fep,tr(kbot:ktop,i0:i1,jp1), &
        dj(kbot:ktop,i0:i1,jp1),dji(kbot:ktop,i0:i1,jp1),i0,i1,kbot,ktop,kminor)
!
! Plasma temperature:
! 'OPLUS_TP0' (tp will get poles from jm1 and jp1)
!
!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,ktop
          tp(k,i,jm1) = te(k,i,jm1)+ti(k,i,jm1)
          tp(k,i,lat) = te(k,i,lat)+ti(k,i,lat)
          tp(k,i,jp1) = te(k,i,jp1)+ti(k,i,jp1)
        enddo
      enddo
!
! Neutral scale height:
! 'OPLUS_HJ' (has constants at poles)
!
!$omp parallel do private(i,k)

      do i=i0,i1
        do k=1,nlev
          hj(k,i,jm1) = gask * tn(k,i,jm1) / (mbar(k,i,jm1) * grav_cm)
          hj(k,i,lat) = gask * tn(k,i,lat) / (mbar(k,i,lat) * grav_cm)
          hj(k,i,jp1) = gask * tn(k,i,jp1) / (mbar(k,i,jp1) * grav_cm)
        enddo
      enddo
! bvel @ jm1 = (B.U)*N(O+)    (J-1)
! bvel @ j   = (B.U)*N(O+)      (J)
! bvel @ jp1 = (B.U)*N(O+)    (J+1)
! 'OPLUS_BVEL' (has constants at poles)
!
! Note bx,by,bz were set globally for all tasks by sub magfield
! (getapex.F90)
!
!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,ktop
          bvel(k,i,jm1) = &
            rri(k,i,jm1)*(bx(i,jm1)*un(k,i,jm1)+by(i,jm1)*vn(k,i,jm1)+   &
            hj(k,i,jm1)*bz(i,jm1)*om(k,i,jm1))*op(k,i,jm1)
          bvel(k,i,lat) = &
           rri(k,i,lat)*(bx(i,lat)*un(k,i,lat)+by(i,lat)*vn(k,i,lat)+   &
            hj(k,i,lat)*bz(i,lat)*om(k,i,lat))*op(k,i,lat)
          bvel(k,i,jp1) = &
            rri(k,i,jp1)*(bx(i,jp1)*un(k,i,jp1)+by(i,jp1)*vn(k,i,jp1)+   &
            hj(k,i,jp1)*bz(i,jp1)*om(k,i,jp1))*op(k,i,jp1)
        enddo ! k=kbot,ktop
      enddo ! i=lon0,lon1
!
! Ambipolar diffusion is returned in diffj:
! 'OPLUS_DIFFJ' (will have constants at poles after this lat scan)
!
      call diffus(tp(kbot:ktop,i0:i1,jm1),ti(kbot:ktop,i0:i1,jm1),te(kbot:ktop,i0:i1,jm1),&
          nop(kbot:ktop,i0:i1,jm1),ne(kbot:ktop,i0:i1,jm1),&
          op(kbot:ktop,i0:i1,jm1),&
          xxxxiop(kbot:ktop,i0:i1,jm1),rmass_fep,hj(kbot:ktop,:,jm1), &
          diffj(kbot:ktop,i0:i1,jm1),diffji(kbot:ktop,i0:i1,jm1),i0,i1,kbot,ktop,lat)
      call diffus(tp(kbot:ktop,i0:i1,lat),ti(kbot:ktop,i0:i1,lat),te(kbot:ktop,i0:i1,lat),&
          nop(kbot:ktop,i0:i1,lat),ne(kbot:ktop,i0:i1,lat),&
          op(kbot:ktop,i0:i1,lat),&
          xxxxiop(kbot:ktop,i0:i1,lat),rmass_fep,hj(kbot:ktop,:,lat), &
          diffj(kbot:ktop,i0:i1,lat),diffji(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat)
      call diffus(tp(kbot:ktop,i0:i1,jp1),ti(kbot:ktop,i0:i1,jp1),te(kbot:ktop,i0:i1,jp1),&
          nop(kbot:ktop,i0:i1,jp1),ne(kbot:ktop,i0:i1,jp1),&
          op(kbot:ktop,i0:i1,jp1),&
          xxxxiop(kbot:ktop,i0:i1,jp1),rmass_fep,hj(kbot:ktop,:,jp1), &
           diffj(kbot:ktop,i0:i1,jp1),diffji(kbot:ktop,i0:i1,jp1),i0,i1,kbot,ktop,lat)
      call diffusj(tp(kbot:ktop,i0:i1,jm1),ti(kbot:ktop,i0:i1,jm1),te(kbot:ktop,i0:i1,jm1),&
          nop(kbot:ktop,i0:i1,jm1),ne(kbot:ktop,i0:i1,jm1),&
          op(kbot:ktop,i0:i1,jm1),rmass_fep,hj(kbot:ktop,:,jm1),dj(kbot:ktop,i0:i1,jm1), &
            dji(kbot:ktop,i0:i1,jm1),diffjj(kbot:ktop,i0:i1,jm1),i0,i1,kbot,ktop,lat)
      call diffusj(tp(kbot:ktop,i0:i1,lat),ti(kbot:ktop,i0:i1,lat),te(kbot:ktop,i0:i1,lat),&
          nop(kbot:ktop,i0:i1,lat),ne(kbot:ktop,i0:i1,lat),&
          op(kbot:ktop,i0:i1,lat),rmass_fep,hj(kbot:ktop,:,lat),dj(kbot:ktop,i0:i1,lat), &
           dji(kbot:ktop,i0:i1,lat),diffjj(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat)
      call diffusj(tp(kbot:ktop,i0:i1,jp1),ti(kbot:ktop,i0:i1,jp1),te(kbot:ktop,i0:i1,jp1),&
          nop(kbot:ktop,i0:i1,jp1),ne(kbot:ktop,i0:i1,jp1),&
          op(kbot:ktop,i0:i1,jp1),rmass_fep,hj(kbot:ktop,:,jp1),dj(kbot:ktop,i0:i1,jp1), &
           dji(kbot:ktop,i0:i1,jp1), diffjj(kbot:ktop,i0:i1,jp1),i0,i1,kbot,ktop,lat)
!
! 'OPLUS_TP1' (constants at the poles)
!
!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,ktop
          tpp(k,i,jm2) = xxxxiop(k,i,jm2)*ti(k,i,jm2)
          tpp(k,i,jm1) = ti(k,i,jm1)*xxxxiop(k,i,jm1)
          tpp(k,i,lat) = ti(k,i,lat)*xxxxiop(k,i,lat)
          tpp(k,i,jp1) = ti(k,i,jp1)*xxxxiop(k,i,jp1)
          tpp(k,i,jp2) = xxxxiop(k,i,jp2)*ti(k,i,jp2)
        enddo
      enddo
!
!
! Latidinal shapiro smoother: opnm is O+ at time n-1.
! opnm_smooth will be used in explicit terms below.
! Smooth in latitude:
! 'OPNM_SMOOTH' (zero at poles)
!
!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,ktop
          opnm_smooth(k,i,lat) = opnm(k,i,lat)-shapiro_const* &
                                (opnm(k,i,jp2)+opnm(k,i,jm2)-4._r8*         &
                                (opnm(k,i,jp1)+opnm(k,i,jm1))+6._r8*        &
                                 opnm(k,i,lat))
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
    enddo ! end first latitude scan (lat=lat0,lat1)
!
!------------------------- End first latitude scan ---------------------
!
! Set pole values for opnm_smooth. Do this before savefld calls, so plots will
! include the poles. All other fields in 1st lat scan got values at the poles 
! via jm1,jp1 above.
!
    call setpoles(opnm_smooth(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)
!
! Save to history file (exclude halo points)
!

    call vxbvel(un,vn,om,hj,uii,vii,wii,1,nlev,i0,i1,j0,j1)
!$omp parallel do private(k, j, i)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          uuii(k,i,j)=uii(k,i,j)
          vvii(k,i,j)=vii(k,i,j)
          wwii(k,i,j)=wii(k,i,j)
          wn(k,i,j)=om(k,i,j)*hj(k,i,j) 
          bdotu(k,i,j) = rri(k,i,j)*(bx(i,j)*un(k,i,j)+by(i,j)*vn(k,i,j)+ &
            wn(k,i,j)*bz(i,j))*bz(i,j)
        enddo
      enddo
    enddo
!
! Set halo points where needed.
!
    nfields = 13
    allocate(ptrs(nfields),polesign(nfields))
    ptrs(1)%ptr => dj ; ptrs(2)%ptr => bvel ; ptrs(3)%ptr => diffj
    ptrs(4)%ptr => tp ; ptrs(5)%ptr => opnm_smooth
    ptrs(6)%ptr =>uuii ; ptrs(7)%ptr => vvii
    ptrs(8)%ptr =>wwii ;ptrs(9)%ptr =>wn 
    ptrs(10)%ptr => diffjj
    ptrs(11)%ptr => diffji
    ptrs(12)%ptr => tpp
    ptrs(13)%ptr => bdotu
    polesign = 1._r8
    polesign(6:7) = -1._r8


    call mp_geo_halos (ptrs,1,nlev,i0,i1,j0,j1,nfields)
    call mp_pole_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields,polesign)

    deallocate(ptrs,polesign)
!----------------------- Begin second latitude scan --------------------
    bdotdh_op  = 0._r8
    bdotdh_opi  = 0._r8
    ambvel  = 0._r8
    bdotdh_opj = 0._r8
    bdotdh_opj1 = 0._r8
    bdotdh_opj2 = 0._r8

    do lat=lat0,lat1
      jm2 = lat-2
      jm1 = lat-1
      jp1 = lat+1
      jp2 = lat+2
!
! bdotdh_op = (B(H).DEL(H))*(D/(H*DZ)*TP+M*G/R)*N(O+)
! then bdotdh_op = d*bz*bdotdh_op
! real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2) :: diffj
! real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2) :: bdotdh_op
! 'BDOTDH_OP' (zero at the poles)
!
      call bdotdh(                    &
        diffj(kbot:ktop,i0:i1,jm1),   &
        diffj(kbot:ktop,:,lat    ),   & ! includes longitude halos
        diffj(kbot:ktop,i0:i1,jp1),   &
        bdotdh_op(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat)
      call bdotdh(                    &
        diffji(kbot:ktop,i0:i1,jm1),   &
        diffji(kbot:ktop,:,lat    ),   & ! includes longitude halos
        diffji(kbot:ktop,i0:i1,jp1),   &
        bdotdh_opi(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat)
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,ktop
          bdotdh_op(k,i,lat) = bz(i,lat)*(dj(k,i,lat)*bdotdh_op(k,i,lat)+&
              dji(k,i,lat)*bdotdh_opi(k,i,lat)) ! BDOTDH_OP
          ambvel(k,i,lat) = (bz(i,lat)**2)*diffjj(k,i,lat) ! BDOTDH_OP
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
!
! bdotdh_opjm1 = (B(H).DEL(H))*2.*TP*N(O+)    (J-1)
! bdotdh_opj   = (B(H).DEL(H))*2.*TP*N(O+)      (J)
! bdotdh_opjp1 = (B(H).DEL(H))*2.*TP*N(O+)    (J+1)
! 'BDOTDH_OPJ' (has reasonable non-constant values at poles)
!
      call bdotdh(                &
        tpp(kbot:ktop,i0:i1,jm2),  &
        tpp(kbot:ktop,:,jm1),      &
        tpp(kbot:ktop,i0:i1,lat),  &
        bdotdh_opj(kbot:ktop,i0:i1,jm1),i0,i1,kbot,ktop,jm1)
      call bdotdh(                &
        tpp(kbot:ktop,i0:i1,jm1),  &
        tpp(kbot:ktop,:,lat),      &
        tpp(kbot:ktop,i0:i1,jp1),  &
        bdotdh_opj(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat)
      call bdotdh(                &
        tpp(kbot:ktop,i0:i1,lat),  &
        tpp(kbot:ktop,:,jp1),      &
        tpp(kbot:ktop,i0:i1,jp2),  &
        bdotdh_opj(kbot:ktop,i0:i1,jp1),i0,i1,kbot,ktop,jp1)

      call firstbdotdh(                &
        tp(kbot:ktop,i0:i1,jm2),  &
        tp(kbot:ktop,:,jm1),      &
        tp(kbot:ktop,i0:i1,lat),  &
        nop(kbot:ktop,i0:i1,jm2),  &
        nop(kbot:ktop,:,jm1),      &
        nop(kbot:ktop,i0:i1,lat),  &
        xxxxiop(kbot:ktop,i0:i1,jm1),  &
        bdotdh_opj1(kbot:ktop,i0:i1,jm1),i0,i1,kbot,ktop,jm1)
      call firstbdotdh(                &
        tp(kbot:ktop,i0:i1,jm1),  &
        tp(kbot:ktop,:,lat),      &
        tp(kbot:ktop,i0:i1,jp1),  &
        nop(kbot:ktop,i0:i1,jm1),  &
        nop(kbot:ktop,:,lat),      &
        nop(kbot:ktop,i0:i1,jp1),  &
        xxxxiop(kbot:ktop,i0:i1,lat),  &
        bdotdh_opj1(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat)
      call firstbdotdh(                &
        tp(kbot:ktop,i0:i1,lat),  &
        tp(kbot:ktop,:,jp1),      &
        tp(kbot:ktop,i0:i1,jp2),  &
        nop(kbot:ktop,i0:i1,lat),  &
        nop(kbot:ktop,:,jp1),      &
        nop(kbot:ktop,i0:i1,jp2),  &
        xxxxiop(kbot:ktop,i0:i1,jp1),  &
        bdotdh_opj1(kbot:ktop,i0:i1,jp1),i0,i1,kbot,ktop,jp1)

      call firstbdotdh(                &
        te(kbot:ktop,i0:i1,jm2),  &
        te(kbot:ktop,:,jm1),      &
        te(kbot:ktop,i0:i1,lat),  &
        ne(kbot:ktop,i0:i1,jm2),  &
        ne(kbot:ktop,:,jm1),      &
        ne(kbot:ktop,i0:i1,lat),  &
        xxxxiop(kbot:ktop,i0:i1,jm1),  &
        bdotdh_opj2(kbot:ktop,i0:i1,jm1),i0,i1,kbot,ktop,jm1)
      call firstbdotdh(                &
        te(kbot:ktop,i0:i1,jm1),  &
        te(kbot:ktop,:,lat),      &
        te(kbot:ktop,i0:i1,jp1),  &
        ne(kbot:ktop,i0:i1,jm1),  &
        ne(kbot:ktop,:,lat),      &
        ne(kbot:ktop,i0:i1,jp1),  &
        xxxxiop(kbot:ktop,i0:i1,lat),  &
        bdotdh_opj2(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat)
      call firstbdotdh(                &
        te(kbot:ktop,i0:i1,lat),  &
        te(kbot:ktop,:,jp1),      &
        te(kbot:ktop,i0:i1,jp2),  &
        ne(kbot:ktop,i0:i1,lat),  &
        ne(kbot:ktop,:,jp1),      &
        ne(kbot:ktop,i0:i1,jp2),  &
        xxxxiop(kbot:ktop,i0:i1,jp1),  &
        bdotdh_opj2(kbot:ktop,i0:i1,jp1),i0,i1,kbot,ktop,jp1)
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,ktop
          bdotdh_opj(k,i,jm1) = bdotdh_opj(k,i,jm1)*dji(k,i,jm1)+&
              bdotdh_opj1(k,i,jm1)*dj(k,i,jm1)+&
              bdotdh_opj2(k,i,jm1)*dji(k,i,jm1)
          bdotdh_opj(k,i,lat) = bdotdh_opj(k,i,lat)*dji(k,i,lat)+&
              bdotdh_opj1(k,i,lat)*dj(k,i,lat)+&
              bdotdh_opj2(k,i,lat)*dji(k,i,lat)
          bdotdh_opj(k,i,jp1) = bdotdh_opj(k,i,jp1)*dji(k,i,jp1)+&
              bdotdh_opj1(k,i,jp1)*dj(k,i,jp1)+&
              bdotdh_opj2(k,i,jp1)*dji(k,i,jp1)
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
    enddo ! lat=j0,j1 (end second lat scan)
!
!------------------------ End second latitude scan ---------------------
!
! bdotdh_opj already has non-constant polar values, but bdotdh_op poles are zero.
! Sub setpoles will set poles to the zonal average of the latitude below each pole.
!
! This may not be necessary, but do it for plotting:
    call setpoles(bdotdh_op(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)
    call setpoles(ambvel(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)

!
! Note mp_geo_halos will overwrite jm1,jp1 that was set above.
! bdotdh_opj needs longitude halos for the bdotdh call below.
!
! real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: bdotdh_op,opj
!
    allocate(ptrs(2),polesign(2))
    ptrs(1)%ptr => bdotdh_opj
    ptrs(2)%ptr => ambvel
    polesign = 1._r8
    call mp_geo_halos (ptrs,1,nlev,i0,i1,j0,j1,2)
    call mp_pole_halos(ptrs,1,nlev,i0,i1,j0,j1,2,polesign)
    deallocate(ptrs,polesign)
!
    !-----------------------------------------------------------------
    !  calculate metal ion velocity 


!----------------------- Begin third latitude scan ---------------------
!
!
    bdotdh_diff = 0._r8
    bdzdvb_op   = 0._r8
    explicit(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a1(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b1(1:nlev,i0:i1,j0:j1)=0._r8
    explicit3(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a2(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b2(1:nlev,i0:i1,j0:j1)=0._r8
    explicit4(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a3(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b3(1:nlev,i0:i1,j0:j1)=0._r8
    explicit5(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a4(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b4(1:nlev,i0:i1,j0:j1)=0._r8
    explicit7(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a5(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b5(1:nlev,i0:i1,j0:j1)=0._r8
    explicit1(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit2(1:nlev,i0:i1,j0:j1)=0._r8; explicit6(1:nlev,i0:i1,j0:j1)=0._r8
    hdz         = 0._r8
    tphdz0      = 0._r8
    tphdz1      = 0._r8
    djint       = 0._r8
    divbz       = 0._r8
    hdzmbz      = 0._r8
    hdzpbz      = 0._r8
    LORHOT= 0._r8
    LNRHOT= 0._r8
    RHOT= 0._r8
    RHOTD= 0._r8
    FABS= 0._r8
    FSGN= 0._r8
    TERM= 0._r8
    TERP= 0._r8
    rwii= 0._r8


!
! gmr = G*M(O+)/(2.*R)
!
    gmr = grav_cm*rmass_fep/(2._r8*gask)
!
! Globally, this loop is lat=2,nlat-1 (i.e., skipping the poles)
!
!    call filter2_op(wwi(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)
    do lat=lat0,lat1
      jm2 = lat-2
      jm1 = lat-1 ! this will be south pole for southern pes (j==1)
      jp1 = lat+1 ! this will be north pole for northern pes (j==nlat)
      jp2 = lat+2
!
!
! bdotdh_opj = (B(H).DEL(H))*D*(B(H).DEL(H))*2.*TP*N(O+)   (J)
! 'BDOTDH_DIFF' (zero at the poles)
!
      call bdotdh(                       &
        bdotdh_opj(kbot:ktop,i0:i1,jm1), &
        bdotdh_opj(kbot:ktop,:,lat),     & ! includes longitude halos
        bdotdh_opj(kbot:ktop,i0:i1,jp1), &
        bdotdh_diff(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat) ! BDOTDH_DIFF
!
! bdzdvb_op = (BZ*D/(H*DZ)+DIV(*B))*S2
! bdzdvb returns bdzdvb_op(k,i).
! 'BDZDVB_OP' (zero at the poles)
!
! real(r8),dimension(i0:i1,j0:j1) :: dvb
! real(r8),dimension(nlev,i0:i1,j0-1:j1+1) :: hj ! scale height
! real(r8),dimension(nlev,i0:i1) :: bdzdvb_op

! subroutine bdzdvb(phi,dvb,h,ans,lev0,lev1,lon0,lon1,lat)
!   real(r8),intent(in) :: dvb(lon0:lon1)
!   real(r8),dimension(lev0:lev1,lon0:lon1),intent(in)    :: phi,h
!   real(r8),dimension(lev0:lev1,lon0:lon1),intent(out)   :: ans
!
      call bdzdvb(bdotdh_opj(kbot:ktop,i0:i1,lat),dvb(:,lat),hj(kbot:ktop,i0:i1,lat), &
        bdzdvb_op(kbot:ktop,i0:i1),kbot,ktop,i0,i1,lat)
!
!
! Collect explicit terms:
! 'EXPLICIT0' (this will have poles set after third lat scan, before
!              plotting. The poles will be constant in longitude, and
!              may differ structurally from adjacent latitudes.
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,ktop
          explicit(k,i,lat) = -(1._r8)*(bdzdvb_op(k,i)+bdotdh_diff(k,i,lat)+ &
            bdotdh_op(k,i,lat))
!            bdotdh_op(k,i,lat)+ddh_op(k,i,lat))
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
!
! Ion drifts are interpolated to midpoints (is this necessary in WACCM?).
!
! Need lon,lat halos for op, bvel, and bmod2
! op,bvel halos were set above, bmod2 was set in magfield (getapex.F90)
! (ui,vi,wi halos are not used here.)
!
! bmod2 halos are set in sub magfield (getapex.F90), including nlat-1,nlat,nlat+1,
!   and 1 halo point in longitude. Note bmod2 is global in lon and lat for all pe's.
! use getapex,only: bmod2  ! (0:nlonp1,jspole-1:jnpole+1)
!
! When looping lat=2,nlat-1, this explicit has zero pole values,
! but there are still problems at processor longitude boundaries,
! especially near the south pole:
! 'EXPLICIT1' (zero at the poles)


!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,ktop
!
! Original TIEGCM statement:
!         explicit(k,i) = explicit(k,i)+1._r8/(2._r8*re)*        &
!           (1._r8/(cs(lat)*dlamda)*(bx(i,lat)*                  &
!           (bvel(k,i+1,lat)-bvel(k,i-1,lat))+                   &
!           0.5_r8*(ui(k,i,lat)+ui(k+1,i,lat))*bmod2(i,lat)**2*  &
!           (op(k,i+1,lat)/bmod2(i+1,lat)**2-                    &
!            op(k,i-1,lat)/bmod2(i-1,lat)**2))+                  &
!
!           1._r8/dphi*(by(i,lat)*(bvel(k,i,jp1)-bvel(k,i,jm1))+ &
!           0.5_r8*(vi(k,i,lat)+vi(k+1,i,lat))*bmod2(i,lat)**2*  &
!           (op(k,i,jp1)/bmod2(i,jp1)**2-                        &
!            op(k,i,jm1)/bmod2(i,jm1)**2)))
!
! Break it into two pieces and put together for debug:
! 'EXPLICITa'
         explicit_a1(k,i,lat) = (bx(i+1,lat)*                       &
            bvel(k,i+1,lat)-bx(i-1,lat)*bvel(k,i-1,lat))
         explicit_a2(k,i,lat) =                        &
            (uui(k,i+1,lat)*xxiop(k,i+1,lat)-                    &
             uui(k,i-1,lat)*xxiop(k,i-1,lat))                
         explicit_a3(k,i,lat) =                       &
            (uuii(k,i+1,lat)*xiop(k,i+1,lat)-                    &
             uuii(k,i-1,lat)*xiop(k,i-1,lat))
         explicit_a4(k,i,lat) =                    &
            (eeeex(k,i+1,lat)*xiop(k,i+1,lat)-   &
            eeeex(k,i-1,lat)*xiop(k,i-1,lat))
         explicit_a5(k,i,lat) =                    &
            (un(k,i+1,lat)*xxxiop(k,i+1,lat)-                    &
             un(k,i-1,lat)*xxxiop(k,i-1,lat))
!
! 'EXPLICITb'
!
         explicit_b1(k,i,lat) = &
            (by(i,jp1)*bvel(k,i,jp1)-by(i,jm1)*bvel(k,i,jm1))
         explicit_b2(k,i,lat) = &
            (vvi(k,i,jp1)*xxiop(k,i,jp1)-                        &
             vvi(k,i,jm1)*xxiop(k,i,jm1))
         explicit_b3(k,i,lat) = &
            (vvii(k,i,jp1)*xiop(k,i,jp1)-                        &
             vvii(k,i,jm1)*xiop(k,i,jm1))
         explicit_b4(k,i,lat) = &
             (eeeey(k,i,jp1)*xiop(k,i,jp1)-   &
             eeeey(k,i,jm1)*xiop(k,i,jm1))
         explicit_b5(k,i,lat) = &
            (vn(k,i,jp1)*xxxiop(k,i,jp1)-                        &
             vn(k,i,jm1)*xxxiop(k,i,jm1))
!            (xiop(k,i,jp1)/bmod2(i,jp1)-                        &
!             xiop(k,i,jm1)/bmod2(i,jm1)))
!
! 'EXPLICIT1'
! explicit will receive polar values after this latitude scan.
         explicit3(k,i,lat) = 1._r8/(2._r8*re)* &
           (1._r8/(cs(lat)*dlamda)*explicit_a1(k,i,lat)+          &
           1._r8/dphi*explicit_b1(k,i,lat))
         explicit4(k,i,lat) = 1._r8/(2._r8*re)* &
           (1._r8/(cs(lat)*dlamda)*explicit_a2(k,i,lat)+          &
           1._r8/dphi*explicit_b2(k,i,lat))
         explicit5(k,i,lat) = 1._r8/(2._r8*re)* &
           (1._r8/(cs(lat)*dlamda)*explicit_a3(k,i,lat)+          &
           1._r8/dphi*explicit_b3(k,i,lat))
         explicit6(k,i,lat) = 1._r8/(2._r8*re)* &
           (1._r8/(cs(lat)*dlamda)*explicit_a4(k,i,lat)+          &
           1._r8/dphi*explicit_b4(k,i,lat))
         explicit7(k,i,lat) = 1._r8/(2._r8*re)* &
           (1._r8/(cs(lat)*dlamda)*explicit_a5(k,i,lat)+          &
           1._r8/dphi*explicit_b5(k,i,lat))
!
         explicit(k,i,lat) = explicit(k,i,lat)+ &
             explicit3(k,i,lat)+  &
             explicit4(k,i,lat)+  &
             explicit5(k,i,lat)+  &
             explicit6(k,i,lat)!+  &
 !            explicit7(k,i,lat)

!
!
! explicit is bad at i=1,72,73,144 near south pole (npole appears to be ok)
! This does not appear to adversely affect the final O+ output, and TIEGCM
! has the same high magnitudes, so am ignoring this for now. The high magnitudes
! are near the south pole, at processor longitude boundaries (implicating an error
! with longitude halo points).
!
         if (debug) then
            if (explicit(k,i,lat) < -300._r8 .or. explicit(k,i,lat) > 300._r8) then
               write(iulog,"('>>> bad explicit: k,i,lat=',3i4,' explicit=',es12.4)") &
                    k,i,lat,explicit(k,i,lat)
               write(iulog,"('  cs(lat)	       =',3es12.4)") cs(lat)
               write(iulog,"('  mgp(k,i-1:i+1,lat)  =',3es12.4)") op(k,i-1:i+1,lat)
               write(iulog,"('  mgp(k,i,jm1:jp1)    =',3es12.4)") op(k,i,jm1:jp1)
               write(iulog,"('  bmod2(i-1:i+1,lat) =',3es12.4)") bmod2(i-1:i+1,lat)
               write(iulog,"('  bmod2(i,jm1:jp1)   =',3es12.4)") bmod2(i,jm1:jp1)
               write(iulog,"('  vxb_ui(k:k+1,i,lat)    =',2es12.4)") ui(k:k+1,i,lat)
               write(iulog,"('  vxb_vi(k:k+1,i,lat)    =',2es12.4)") vi(k:k+1,i,lat)
               write(iulog,"('  bx,by(i,lat)       =',2es12.4)") bx(i,lat),by(i,lat)
            endif
         endif

        enddo ! k=kbot,nlev-1
      enddo ! i=i0,i1
!$omp parallel do private( k )

!      do i=i0,i1
!        do k=kbot,ktop
!          explicit(k,i,lat)=explicit(k,i,lat)+0.5_r8*(xiop(k,i,lat)+xiop(k+1,i,lat))*dvvxb(k,i,lat)+ &
!              0.5_r8*(xiop(k,i,lat)+xiop(k+1,i,lat))/bmod2(i,lat)*dve(k,i,lat)+&
!              0.5_r8*(xxxiop(k,i,lat)+xxxiop(k+1,i,lat))*dvv(k,i,lat)+&!              0.5_r8*(xxiop(k,i,lat)+xxiop(k+1,i,lati))*dvexb(k,i,lat)
 !             0.5_r8*(xxiop(k,i,lat)+xxiop(k+1,i,lat))*dvexb(k,i,lat)
!            if (explicit(k,i,lat) < -0.1_r8) then
!                explicit(k,i,lat)=-0.1_r8
!            endif
!            if (explicit(k,i,lat) > 0.1_r8) then
!                explicit(k,i,lat)=0.1_r8
!            endif
 !       enddo ! k=kbot,nlev
 !     enddo ! i=i0,i1


!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,ktop
          rwii(k,i,lat)=0.5_r8*(wwii(k,i,lat)*ri(k,i,lat)+wwii(k+1,i,lat)*ri(k+1,i,lat))&
!              +0.5_r8*(wn(k,i,lat)*rrri(k,i,lat)+wn(k+1,i,lat)*rrri(k+1,i,lat))&
              +0.5_r8*(wwi(k,i,lat)*rri(k,i,lat)+wwi(k+1,i,lat)*rri(k+1,i,lat))&
              +0.5_r8*(eeeez(k,i,lat)*ri(k,i,lat)+eeeez(k+1,i,lat)*ri(k+1,i,lat))&
              +0.5_r8*(bdotu(k,i,lat)+bdotu(k+1,i,lat))&
              -0.5_r8*(rrrri(k,i,lat)*(ambvel(k,i,lat))+&
              rrrri(k+1,i,lat)*(ambvel(k+1,i,lat)))
!              -0.5_r8*(rrri(k,i,lat)*(ambvel1(k,i,lat))+&
!              rrri(k+1,i,lat)*(ambvel1(k+1,i,lat)))
!              -0.5_r8*(rrrri(k,i,lat)*(ambvel(k,i,lat)+mambvel(k,i,lat))+&
!              rrrri(k+1,i,lat)*(ambvel(k+1,i,lat)+mambvel(k+1,i,lat)))
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
!      do i=i0,i1
!        do k=kbot,ktop
!           rwii(k,i,lat)=-1._r8*rwii(k,i,lat)
!          if (rwii (k,i,lat) > 1.e4_r8) rwii  (k,i,lat) = 1.e4_r8
!          if (rwii (k,i,lat) < -1.e4_r8) rwii  (k,i,lat) = -1.e4_r8
!        enddo ! k=lev0,lev1-1
!      enddo ! i=lon0,lon1
!          p_coeff(k,i,lat) =   hdzmbz(k,i)*djint(k  ,i)*tphdz0(k  ,i)

! Begin coefficients p_coeff, q_coeff, r_coeff
!
!
! Continue coefficients with vertical ion drift:
! wi is converted from interfaces to midpoints (first use of wi).
! The p,q,r coeffs are still zero at top boundary k=nlev, and at poles.
!
! Continue coefficients with vertical ion drift:
! wi is converted from interfaces to midpoints (first use of wi).
! The p,q,r coeffs are still zero at top boundary k=nlev, and at poles.
!
! Then TIEGCM "Add source term to RHS (explicit terms)", and calculates
! lower boundary condition N(O+) = Q/L (q_coeff, explicit, p_coeff), and
! finally calls trsolv.
!
 300 continue
    enddo ! end third latitude scan (lat=lat0,lat1)
!    call filter1_op(rwii(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)
!    call filter2_op(rwii(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)
!    call filter1_op(explicit(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)
!    call filter2_op(explicit(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)
!
!------------------------ End third latitude scan ---------------------
    ktop=ktop-1
    kbot=kbot+1
    ONE6TH = 1.0_r8/6.0_r8
    ONE3RD = 1.0_r8/3.0_r8
!$omp parallel do private(k, j, i)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          aa(k,i,j)=re*re*dphi*dlamda*cs(j)
!          aa(k,i,j)=1.0_r8
          va(k,i,j)=aa(k,i,j)*hj(k,i,j)*dzp
          rva(k,i,j)=1._r8/va(k,i,j)
        enddo
      enddo
    enddo
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          hadudth(k,i,j)=dtime*rwii(k,i,j)*aa(k,i,j)
          epsh(k,i,j) =hadudth(k,i,j)*rva(k,i,j)
          scrh(k,i,j) = min ( ONE6TH, ABS(epsh(k,i,j)) )
          scrh(k,i,j) = ONE3RD*scrh(k,i,j)**2
          hadudth(k,i,j)=0.5_r8*hadudth(k,i,j)
          nulh(k,i,j) = ONE6TH + ONE3RD*(epsh(k,i,j) + scrh(k,i,j))* &
                                        (epsh(k,i,j) - scrh(k,i,j))
          mulh(k,i,j) =  0.25_r8 - 0.5_r8*nulh(k,i,j)
          nulh(k,i,j) = va(k,i,j)*(nulh(k,i,j) + scrh(k,i,j))
          mulh(k,i,j) = va(k,i,j)*(mulh(k,i,j) + scrh(k,i,j))
        enddo
      enddo
    enddo
!$omp parallel do private(k, j, i)
      do j=j0,j1
        do i=i0,i1
          op(ktop+1,i,j) = abk(i,j)*op(ktop,i,j)
        enddo
      enddo
    do k=2,ktop+1
      do j=j0,j1
        do i=i0,i1
          flxh(k,i,j) = hadudth(k,i,j)* ( op(k,i,j) + op(k-1,i,j) )
          diff(k,i,j) = nulh(k,i,j) * ( op(k,i,j) - op(k-1,i,j) )
        enddo
      enddo
    enddo
!$omp parallel do private(lat, i, k )
    do lat=lat0,lat1
      do i=i0,i1
        do k=kbot,ktop
             LORHOT(k,i,lat) = va(k,i,lat)*op(k,i,lat) - &
                 explicit(k,i,lat)*dtime*va(k,i,lat) +&
                 (flxh(k,i,lat)-flxh(k+1,i,lat))
             LNRHOT(k,i,lat) = LORHOT(k,i,lat) + (diff(k+1,i,lat) - diff(k,i,lat))
             RHOT(k,i,lat)  = LORHOT(k,i,lat)*rva(k,i,lat)
             RHOTD(k,i,lat) = LNRHOT(k,i,lat)*rva(k,i,lat)
        enddo ! k=kbot,nlev
        RHOT(kbot-1,i,lat)=RHOT(kbot,i,lat)
        RHOTD(kbot-1,i,lat)=RHOTD(kbot,i,lat)
        RHOT(ktop+1,i,lat)=RHOT(ktop,i,lat)*abk(i,lat)
        RHOTD(ktop+1,i,lat)=RHOTD(ktop,i,lat)*abk(i,lat)
      enddo ! i=i0,i1
      do i=i0,i1
        do k=kbot,ktop
             flxh(k,i,lat) = mulh(k,i,lat) * ( RHOT(k,i,lat) - RHOT(k-1,i,lat) )
             diff(k,i,lat) = RHOTD(k,i,lat) - RHOTD(k-1,i,lat)
        enddo ! k=kbot,nlev
        flxh(ktop+1,i,lat)=mulh(ktop+1,i,lat) * ( RHOT(ktop+1,i,lat) - RHOT(ktop,i,lat) )
        diff(ktop+1,i,lat)=RHOTD(ktop+1,i,lat) - RHOTD(ktop,i,lat)
      enddo ! i=i0,i1
      do i=i0,i1
        do k=kbot,ktop
             FABS(k+1,i,lat) = ABS ( flxh (k+1,i,lat) )
             FSGN(k+1,i,lat) = SIGN ( 0.999_r8, diff(k+1,i,lat) )
             TERM(k+1,i,lat) = FSGN(k+1,i,lat)*va(k,i,lat)*diff(k,i,lat)
             TERP(k,i,lat) = FSGN(k,i,lat)*va(k,i,lat)*DIFF(k+1,i,lat)
        enddo ! k=kbot,nlev
        TERP(ktop+1,i,lat)=BIGNUM
        TERM(kbot,i,lat)=BIGNUM
      enddo ! i=i0,i1
      do i=i0,i1
        flxh(kbot,i,lat) = FSGN(kbot,i,lat) * max ( 0.0_r8, &
                      min ( TERM(kbot,i,lat), FABS(kbot,i,lat), TERP(kbot,i,lat) ) )
        do k=kbot,ktop
             flxh(k+1,i,lat) = FSGN(k+1,i,lat) * MAX ( 0.0_r8, &
                     MIN ( TERM(k+1,i,lat), FABS(k+1,i,lat), TERP(k+1,i,lat) ) )
             op_out(k,i,lat) =  ( LNRHOT(k,i,lat) + (flxh(k,i,lat) - flxh(k+1,i,lat)) )*rva(k,i,lat)
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
    enddo ! end third latitude scan (lat=lat0,lat1)
!
! Set poles for selected diagnostics:
!
!
! All tasks have global 2d bmod2.
! bmod2 was set by sub magfield (getapex.F90)
!   allocate(bmod2(0:nlonp1,jspole-1:jnpole+1))
! Copy bmod2 poles to diagnostic array.
!
!
! Set poles for selected diagnostics:
!
!
    call setpoles(explicit(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)
!
! Call solver, defining O+ output op_out:
!
! Its best not to call this unless the coefficients and explicit terms
! have been properly set in the third latitude scan above (e.g., during 
! "goto 300" debugging above, where the coeffs may not have been calculated).
!
!
! Write fields from third latitude scan to waccm history:
!
!------------------------------------------------------------------------
!
! Filter O+ output from solver:
! (TIMEGCM calls both filters, whereas TIEGCM calls only filter2)
!
!   call filter2_op(op_out(kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1)
!
!   call filter2_op(op_out(kbot:ktop,i0:i1,j0:j1),kbot,ktop,i0,i1,j0,j1)
!
!----------------------- Begin fourth latitude scan ---------------------
!
!$omp parallel do private(lat, i, k)
    do lat=j0,j1
!      do i=i0,i1
!          op_out(ktop,i,lat)=op_out(ktop+1,i,lat)
!          op_out(ktop-1,i,lat)=op_out(ktop+1,i,lat)
!          op_out(ktop-2,i,lat)=op_out(ktop+1,i,lat)
!          op_out(kbot,i,lat)=op_out(kbot-1,i,lat)
!          op_out(kbot+2,i,lat)=op_out(kbot-1,i,lat)
!      enddo

      do i=i0,i1
        do k=kbot,ktop
          opnm_out(k,i,lat) = dtsmooth*op(k,i,lat)+dtsmooth_div2* &
            (opnm(k,i,lat)+op_out(k,i,lat))
        enddo
      enddo
!
! Insure non-negative O+ output:
      do i=i0,i1
        do k=kbot,ktop
          if (op_out  (k,i,lat) < 1.e-7_r8) op_out  (k,i,lat) = 1.e-7_r8
          if (opnm_out(k,i,lat) < 1.e-7_r8) opnm_out(k,i,lat) = 1.e-7_r8
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!      do i=i0,i1
!        do k=kbot,ktop
!          if (op_out  (k,i,lat) > 1.e4_r8) op_out  (k,i,lat) = 1.e4_r8
!          if (opnm_out(k,i,lat) > 1.e4_r8) opnm_out(k,i,lat) = 1.e4_r8
!        enddo ! k=lev0,lev1-1
!      enddo ! i=lon0,lon1
!
! Enforce O+ minimum if enforce_opfloor is true.
! Opfloor is Stan's "smooth floor" (product of two Gaussians, 
!   dependent on latitude and pressure level) (opmin=3000.0):
!
!      if (enforce_floor) then
!         zpmid(kbot:ktop) = log(50.e-6_r8/pmid(kbot:ktop)) ! tgcm levs -- maybe done once at init time
!         do k=kbot,ktop
!            opfloor = opmin*exp(-(glat(lat)/90.0_r8)**2/0.3_r8) &
!                 *exp(-((zpmid(k)-4.25_r8)/zpmid(ktop))**2/0.1_r8)
!            do i=i0,i1
!               if (op_out(k,i,lat) < opfloor) then
!                  op_out(k,i,lat) = opfloor
!               endif ! opout < opfloor
!            enddo ! i=lon0,lon1
!         enddo ! k=lev0,lev1-1
!      endif ! enforce_opfloor

    enddo ! lat=lat0,lat1
!
! Save O+ output to WACCM history (cm^3):
  end subroutine vxb_xport
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine divb(dvb,i0,i1,j0,j1)
!
! Evaluate divergence of B, the unit magnetic field vector.
! (all processors have the full global 2d field)
!
! Args:
    integer,intent(in) :: i0,i1,j0,j1
    real(r8),intent(out) :: dvb(i0:i1,j0:j1)
!
! Local:
    integer :: i,j,jm1,jp1
    real(r8),parameter :: re = 6.37122e8_r8  ! earth radius (cm)

    dvb = 0._r8

!
! Note re is in cm.
! (bx,by,bz are set by sub magfield (getapex.F90))
! (dphi,dlamda, and cs are set by sub set_geogrid (edyn_init.F90))
!
    do j=j0,j1
      jm1 = j-1
      jp1 = j+1
      do i=i0,i1
        dvb(i,j) = (((bx(i+1,j)-bx(i-1,j))/(2._r8*dlamda)+      &
          (cs(jp1)*by(i,jp1)-cs(jm1)*by(i,jm1))/(2._r8*dphi))/  &
          cs(j)+2._r8*bz(i,j))/re
      enddo ! i=i0,i1
    enddo ! j=j0,j1
  end subroutine divb
!-----------------------------------------------------------------------
  subroutine rrk(t,ti,rms,ps1,ps2,n2,nop,en,rmass_fep,tr,ans,ans1,lon0,lon1,lev0,lev1,levn)
!
! Returns ambipolar diffusion coefficient in ans.
!
! Args:
    integer,intent(in) :: lon0,lon1,lev0,lev1,levn
    real(r8),intent(in) :: rmass_fep
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(in) :: &
      t,ti,rms,ps1,ps2,n2,tr,nop,en
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out) :: ans, ans1
    real(r8),dimension(lev0:lev1,lon0:lon1) :: nopn,enn
!
! Local:
!
    integer :: k,i
      do i=lon0,lon1
        do k=lev0,lev1
        nopn(k,i)=nop(k,i)
        enn(k,i)=en(k,i)
          if (nopn  (k,i) < 1.e-5_r8) nopn  (k,i) = 1.e-5_r8
          if (enn  (k,i) < 1.e-5_r8) enn  (k,i) = 1.e-5_r8
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!$omp parallel do private(i,k)
    do i=lon0,lon1
      do k=lev0,lev1-1
        ans(k,i) = 1.42e17_r8*boltz*t(k,i)/(p0*expz(k)*.5_r8*(rms(k,i)+    &
          rms(k+1,i))*(ps2(k,i)*rmassinv_o1*sqrt(tr(k,i))*(1._r8-0.064_r8* &
          log10(tr(k,i)))**2*colfac+18.6_r8*n2(k,i)*rmassinv_n2+18.1_r8*   &
          ps1(k,i)*rmassinv_o2))

!        ans(k,i) = 0._r8
!        ans1(k,i) = 0._r8
        ans1(k,i) = 0._r8
!        ans(k,i) = 8.31e7_r8/coo(k,i)/rmass_fep

      enddo ! k=lev0,lev1
      ans(lev1,i) = ans(lev1-1,i) ! should not need to do this
      ans1(lev1,i) = ans1(lev1-1,i) ! should not need to do this

    enddo ! i=lon0,lon1
    do i=lon0,lon1
      do k=levn,lev1-1
        ans1(k,i) = gask/(rmass_fep*0.0851_r8*sqrt((rmass_fep+16._r8)/(rmass_fep*16._r8))* &
            ti(k,i)**(-1.5_r8)*nopn(k,i)*log(10518._r8*ti(k,i)**1.5_r8/sqrt(nopn(k,i))))
      enddo ! k=lev0,lev1
      ans1(lev1,i) = ans1(lev1-1,i) ! should not need to do this

    enddo ! i=lon0,lon1
      do i=lon0,lon1
        do k=lev0,lev1
          if (nopn  (k,i) < 2.e4_r8) then
!          if (nopn  (k,i) < 5.e2_r8 .or. nopn  (k,i) < 100._r8*enn(k,i)) then
!          if (nopn  (k,i) < 5.e2_r8 .or. nopn  (k,i) < 100._r8*enn(k,i)) then
              ans1  (k,i) = 0._r8
          endif
          if (nopn  (k,i) < 1.e3_r8) then
              ans  (k,i) = 0._r8
          endif
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! Cap ambipolar diffusion coefficient in ans.
!
    ! acceptable range for limiter 1.e8 to 1.e9 ...
    where( ans(:,:) > adiff_limiter )
      ans(:,:) = adiff_limiter
    endwhere
    where( ans1(:,:) > adiff_limiter )
      ans1(:,:) = adiff_limiter
    endwhere

  end subroutine rrk
!-----------------------------------------------------------------------
  subroutine diffus(tp,ti,te,nop,ne,en,een,rmass_fep,hj,ans,ans1,i0,i1,lev0,lev1,lat)
!                                      kbot,nlev
! Evaluates ans = (d/(h*dz)*tp+m*g/r)*en
! Remember: "bot2top": lev0=kbot=bottom, lev1=nlev=top
!
! Args:
    integer :: i0,i1,lev0,lev1,lat
    real(r8),intent(in) :: rmass_fep
    real(r8),dimension(lev0:lev1,i0:i1),intent(in) :: tp,ti,te,nop,ne,en,een,hj
    real(r8),dimension(lev0:lev1,i0:i1),intent(out) :: ans,ans1
    real(r8),dimension(lev0:lev1,i0:i1) :: nopn,enn
!
! Local:
    integer :: i,k
    real(r8) :: mgr,mgro
      do i=i0,i1
        do k=lev0,lev1
        nopn (k,i)=nop  (k,i)
        enn (k,i)=en  (k,i)
          if (nopn  (k,i) < 1.e-5_r8) nopn  (k,i) = 1.e-5_r8
          if (enn(k,i) < 1.e-5_r8) enn(k,i) = 1.e-5_r8
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

    mgr = rmass_fep*grav_cm/gask
    mgro = 16._r8*grav_cm/gask

!$omp parallel do private(i,k)
    do i=i0,i1
      do k=lev0,lev1-2
        ans(k+1,i) =1.0_r8/(2._r8*hj(k+1,i)*dzp)*&
        (nopn(k+2,i)*tp(k+2,i)-nopn(k,i)*tp(k,i))/nopn(k+1,i)*een(k+1,i)&
!        ans(k+1,i) = 1.0_r8/(2._r8*hj(k+1,i)*dzp)*((en(k+2,i)*ti(k+2,i)- &
!          en(k,i)*ti(k,i))+(en(k+2,i)*te(k+2,i)-en(k,i)*te(k,i)))&
          +mgro*een(k+1,i)
        ans1(k+1,i) =1.0_r8/(2._r8*hj(k+1,i)*dzp)*&
        (((enn(k+2,i)*ti(k+2,i)-enn(k,i)*ti(k,i))/enn(k+1,i)&
        +(ne(k+2,i)*te(k+2,i)-ne(k,i)*te(k,i))/ne(k+1,i))*een(k+1,i))&
          +mgr*een(k+1,i)
      enddo
      if (debug) then
        write(iulog,"('diffus: lat=',i4,' i=',i4,' ans(lev0:lev1-1,i)=',2es12.4)") &
          lat,i,minval(ans(lev0:lev1-1,i)),maxval(ans(lev0:lev1-1,i))
      endif
    enddo
!
! Upper and lower boundaries:
!
!$omp parallel do private(i)
    do i=i0,i1
!
! Upper boundary:
!      ans(lev1,i) = 1._r8/(hj(lev1,i)*dzp)*((ti(lev1,i)*en(lev1,i)- &
!        ti(lev1-1,i)*en(lev1-1,i))+((ne(lev1,i)*te(lev1,i)-ne(lev1-1,i)&
!        *te(lev1-1,i))/(ne(lev1,i))*en(lev1,i)))+mgr*en(lev1,i)
      ans(lev1,i) = 0._r8
      ans1(lev1,i) = 0._r8
!
! Lower boundary:
      ans(lev0,i) = 1._r8/(hj(lev0,i)*dzp)*(tp(lev0+1,i)*nopn(lev0+1,i)- &
        tp(lev0,i)*nopn(lev0,i))/nopn(lev0,i)*een(lev0,i)+mgro*een(lev0,i)
      ans1(lev0,i) = 1._r8/(hj(lev0,i)*dzp)*&
          (((ti(lev0+1,i)*enn(lev0+1,i)-ti(lev0,i)*enn(lev0,i))/enn(lev0,i)+&
          (te(lev0+1,i)*ne(lev0+1,i)-te(lev0,i)*ne(lev0,i))/ne(lev0,i))*een(lev0,i))+mgr*een(lev0,i)
!      ans(lev0,i) = 1._r8/(hj(lev0,i)*dzp)*((ti(lev0+1,i)*en(lev0+1,i)- &
!        ti(lev0,i)*en(lev0,i))+(en(lev0+1,i)*te(lev0+1,i)-en(lev0,i)&
!        *te(lev0-1,i)))+mgr*en(lev0,i)
    enddo
      do i=i0,i1
        do k=lev0,lev1
          if (nopn  (k,i) < 2.e4_r8) then
!          if (nopn  (k,i) < 5.e2_r8 .or. nopn  (k,i) < 100._r8*enn(k,i)) then
              ans1  (k,i) = 0._r8
          endif
          if (nopn  (k,i) < 1.e3_r8) then
              ans  (k,i) = 0._r8
          endif
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
  end subroutine diffus
!-----------------------------------------------------------------------
  subroutine diffusj(tp,ti,te,nop,ne,en,rmass_fep,hj,dj,dji,ans,i0,i1,lev0,lev1,lat)
!                                      kbot,nlev
! Evaluates ans = (d/(h*dz)*tp+m*g/r)*en
! Remember: "bot2top": lev0=kbot=bottom, lev1=nlev=top
!
! Args:
    integer :: i0,i1,lev0,lev1,lat
    real(r8),intent(in) :: rmass_fep
    real(r8),dimension(lev0:lev1,i0:i1),intent(in) :: tp,ti,te,nop,ne,en,hj,dj,dji
    real(r8),dimension(lev0:lev1,i0:i1),intent(out) :: ans
    real(r8),dimension(lev0:lev1,i0:i1) :: ans1,ans2
    real(r8),dimension(lev0:lev1,i0:i1) :: nopn,enn
!
! Local:
    integer :: i,k
    real(r8) :: mgr,mgro
      do i=i0,i1
        do k=lev0,lev1
        nopn (k,i)=nop  (k,i)
        enn (k,i)=en  (k,i)
          if (nopn  (k,i) < 1.e-5_r8) nopn  (k,i) = 1.e-5_r8
          if (enn(k,i) < 1.e-5_r8) enn(k,i) = 1.e-5_r8
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

    mgr = rmass_fep*grav_cm/gask
    mgro = 16._r8*grav_cm/gask
!    where( coo(:,:) <1.e-1_r8  )
!      coo(:,:) = 1.e-1_r8
!    endwhere

!$omp parallel do private(i,k)
    do i=i0,i1
      do k=lev0,lev1-2
        ans1(k+1,i) =1.0_r8/(2._r8*hj(k+1,i)*dzp)*&
            ((nopn(k+2,i)*tp(k+2,i)- nopn(k,i)*tp(k,i))/nopn(k+1,i))&
          *dj(k+1,i)+mgro*dj(k+1,i)
!        ans1(k+1,i) =mgro*dj(k+1,i)
          if (nopn  (k,i) < 1.e3_r8) then
              ans1  (k,i) = 0._r8
          endif
 !       ans1(k+1,i) = 0.0_r8
 !      ans2(k+1,i) = 0.0_r8
        ans2(k+1,i) = 1.0_r8/(2._r8*hj(k+1,i)*dzp)*((enn(k+2,i)*ti(k+2,i)- &
          enn(k,i)*ti(k,i))/(enn(k+1,i))+(ne(k+2,i)*te(k+2,i)-ne(k,i)*te(k,i))/(ne(k+1,i)))&
          *dji(k+1,i)+mgr*dji(k+1,i)
!        ans2(k+1,i) = mgr*dji(k+1,i)
          if (enn  (k,i) < 1.e-4_r8) then
              ans2  (k,i) = 0._r8
          endif
        ans(k+1,i)=ans1(k+1,i)+ans2(k+1,i)
!        ans(k+1,i) = 1.0_r8/(2._r8*hj(k+1,i)*dzp)*((ti(k+2,i)-ti(k,i)))&
!          *boltz*6.02e23_r8/(rmass_fep*coo(k+1,i))
!        ans(k+1,i) = grav_cm/coo(k+1,i)
      enddo
      if (debug) then
        write(iulog,"('diffus: lat=',i4,' i=',i4,' ans(lev0:lev1-1,i)=',2es12.4)") &
          lat,i,minval(ans(lev0:lev1-1,i)),maxval(ans(lev0:lev1-1,i))
      endif
    enddo

! Upper and lower boundaries:
!
!$omp parallel do private(i)
    do i=i0,i1
!
! Upper boundary:
!        ans(lev1,i) = 1.0_r8/(hj(lev1,i)*dzp)*((en(lev1,i)*ti(lev1,i)- &
!          en(lev1-1,i)*ti(lev1-1,i))/(en(lev1,i))+(ne(lev1,i)*te(lev1,i)-ne(lev1-1,i)*te(lev1-1,i))/(ne(lev1,i)))&
 !         *boltz*6.02e23_r8/(rmass_fep*coo(lev1,i))+grav_cm/coo(lev1,i)
        ans(lev1,i) = 0.0_r8
!      ans(lev1,i) = grav_cm/coo(lev1,i)
!

!        ans1(i) = 1._r8-((ti(lev1,i)- &
!          ti(lev1-1,i))+(nop(lev1,i)*te(lev1,i)-nop(lev1-1,i)*te(lev1-1,i))/nop(lev1,i)&
!          +rmass_fep*grav_cm/gask*hj(lev1,i)*dzp)/ti(lev1,i)
!          ti(lev1-1,i))+(te(lev1,i)-te(lev1-1,i))&
!          +rmass_fep*grav_cm/gask*hj(lev1,i)*dzp)/(ti(lev1,i)+te(lev1,i))
! Lower boundary:
!        ans1(lev0,i) =1.0_r8/(hj(lev0,i)*dzp)*&
!            ((ne(lev0+1,i)*tp(lev0+1,i)-ne(lev0,i)*tp(lev0,i))/ne(lev0,i))&
!          *dj(lev0,i)+mgro*dj(lev0,i)
        ans1(lev0,i) =1.0_r8/(hj(lev0,i)*dzp)*&
            ((nopn(lev0+1,i)*tp(lev0+1,i)-nopn(lev0,i)*tp(lev0,i))/nopn(lev0,i))&
          *dj(lev0,i)+mgro*dj(lev0,i)
          if (nopn  (lev0,i) < 1.e3_r8) then
              ans1  (lev0,i) = 0._r8
          endif
 !       ans1(k+1,i) = 0.0_r8
!        ans1(lev0,i) = 0._r8
!        ans2(lev0,i) = 0._r8
        ans2(lev0,i) = 1.0_r8/(hj(lev0,i)*dzp)*((enn(lev0+1,i)*ti(lev0+1,i)- &
          enn(lev0,i)*ti(lev0,i))/(enn(lev0,i))+(ne(lev0+1,i)*te(lev0+1,i)-ne(lev0,i)*te(lev0,i))/(ne(lev0,i)))&
          *dji(lev0,i)+mgr*dji(lev0,i)
!        ans2(lev0,i) = mgr*dji(lev0,i)
          if (enn  (lev0,i) < 1.e-4_r8) then
              ans2  (lev0,i) = 0._r8
          endif
        ans(lev0,i)=ans1(lev0,i)+ans2(lev0,i)
!        ans(lev0,i) = 1.0_r8/(hj(lev0,i)*dzp)*( &
!          (ti(lev0+1,i)-ti(lev0,i)))&
!          *boltz*6.02e23_r8/(rmass_fep*coo(lev0,i))
!      ans(lev0,i) = grav_cm/coo(lev0,i)
    enddo
    where( ans(:,:) >2.0e4_r8  )
      ans(:,:) = 2.0e4_r8
    endwhere
    where( ans(:,:) <-2.0e4_r8  )
      ans(:,:) = -2.0e4_r8
    endwhere
  end subroutine diffusj
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine bdotdh(phijm1,phij,phijp1,ans,lon0,lon1,lev0,lev1,lat)
!
! Evaluates ans = (b(h)*del(h))*phi
!
! Args:
    integer,intent(in) :: lon0,lon1,lev0,lev1,lat
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(in) :: phijm1,phijp1
    real(r8),dimension(lev0:lev1,lon0-2:lon1+2),intent(inout) :: phij ! why intent(inout)?
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out) :: ans
!
! Local:
    integer :: k,i
!
! Note phij longitude dimension is lon0-2:lon1+2 (only i-1 and i+1 are used).
! Halo longitudes i-1 and i+1 must have been set before this routine is
! called. ('by' is use-associated above)
!
!$omp parallel do private( i, k )
    do i=lon0,lon1
      do k=lev0,lev1
        ans(k,i) = 1._r8/re*(bx(i,lat)/(cs(lat)*2._r8*dlamda)* &
          (phij(k,i+1)-phij(k,i-1))+by(i,lat)*                 &
          (phijp1(k,i)-phijm1(k,i))/(2._r8*dphi))
      enddo ! k=lev0,lev1
    enddo ! i=lon0,lon1
!
  end subroutine bdotdh
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine firstbdotdh(phijm1,phij,phijp1,nijm1,nij,nijp1,enn,ans,lon0,lon1,lev0,lev1,lat)
!
! Evaluates ans = (b(h)*del(h))*phi
!
! Args:
    integer,intent(in) :: lon0,lon1,lev0,lev1,lat
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(in) :: phijm1,phijp1,nijm1,nijp1,enn
    real(r8),dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: phij, nij! why intent(inout)?
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out) :: ans
!
! Local:
    integer :: k,i
!
! Note phij longitude dimension is lon0-2:lon1+2 (only i-1 and i+1 are used).
! Halo longitudes i-1 and i+1 must have been set before this routine is
! called. ('by' is use-associated above)
!
!$omp parallel do private( i, k )
    ans = 0._r8
    do i=lon0,lon1
      do k=lev0,lev1
          if (nij  (k,i) > 1.e3_r8) then
        ans(k,i) = 1._r8/re*(bx(i,lat)/(cs(lat)*2._r8*dlamda)* &
          (phij(k,i+1)*nij(k,i+1)-phij(k,i-1)*nij(k,i-1))+by(i,lat)*                 &
          (phijp1(k,i)*nijp1(k,i)-phijm1(k,i)*nijm1(k,i))/(2._r8*dphi))/ &
          nij(k,i)*enn(k,i)
          endif
      enddo ! k=lev0,lev1
    enddo ! i=lon0,lon1
!
  end subroutine firstbdotdh
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine bdzdvb(phi,dvb,h,ans,lev0,lev1,lon0,lon1,lat)
!
! Evaluates  ans = (bz*d/(h*dz)+divb)*phi
!
! Args:
    integer,intent(in) :: lev0,lev1,lon0,lon1,lat
    real(r8),intent(in) :: dvb(lon0:lon1)
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(in)    :: phi,h
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out)   :: ans
!
! Local:
    integer :: k,i
!
!$omp parallel do private( i, k )
    do i=lon0,lon1
      do k=lev0+1,lev1-1
        ans(k,i) = bz(i,lat)/(2._r8*h(k,i)*dzp)*(phi(k+1,i)-phi(k-1,i))+ &
          dvb(i)*phi(k,i)
      enddo ! k=lev0+1,lev1-1
    enddo ! i=lon0,lon1
!
! Upper and lower boundaries:
!$omp parallel do private( i )
    do i=lon0,lon1
      ans(lev1,i) = bz(i,lat)/(h(lev1,i)*dzp)*(phi(lev1,i)- &
        phi(lev1-1,i))+dvb(i)*phi(lev1,i)
      ans(lev0,i) = bz(i,lat)/(h(lev0,i)*dzp)* &
        (phi(lev0+1,i)-phi(lev0,i))+dvb(i)*phi(lev0,i)
    enddo ! i=lon0,lon1
  end subroutine bdzdvb
!-----------------------------------------------------------------------
  subroutine vxbvel(un,vn,om,hj,umi,vmi,wmi,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Calculate 3d VnxB metal ion drifts
! on geographic grid.
!
! Args:
    integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
    real(r8),intent(in),dimension(lev0:lev1,lon0:lon1,lat0-1:lat1+1) :: &
      hj  ! geopotential from input (cm)
    real(r8),intent(in),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) :: &
      un,vn,om  ! geopotential from input (cm)
    real(r8),intent(out),dimension(lev0:lev1,lon0:lon1,lat0:lat1) :: &
      umi,vmi,wmi
!
! Local:
    integer :: i,ii,k,j
!
! Scan geographic latitude subdomain:
!
    do j=lat0,lat1
!
! umi = zonal, vmi = meridional, wmi = vertical
!
      do k=lev0,lev1
        do i=lon0,lon1
          ii = i
          umi(k,i,j) = vn(k,i,j)*bz(ii,j)-hj(k,i,j)*om(k,i,j)*by(ii,j)
          vmi(k,i,j) = hj(k,i,j)*om(k,i,j)*bx(ii,j)-un(k,i,j)*bz(ii,j)
          wmi(k,i,j) = un(k,i,j)*by(ii,j)-vn(k,i,j)*bx(ii,j)
        enddo ! i=lon0,lon1
      enddo ! k=lev0,lev1
    enddo ! j=lat0,lat1

  end subroutine vxbvel
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1)
!
! Tri-diagonal solver.
!   a(k,i)*x(k-1,i) + b(k,i)*x(k,i) + c(k,i)*x(k+1,i) = f(k,i)
!
    implicit none
!
! Args:
    integer,intent(in) :: lev0,lev1,k1,k2,lon0,lon1
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(in) :: &
      a, & ! input coefficients
      b, & ! input coefficients
      c, & ! input coefficients
      f    ! input RHS
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out) :: &
      x  ! output
!
! Local:
    integer :: k,kk,i
    real(r8),dimension(lev0:lev1,lon0:lon1) :: w1,w2,w3  ! work arrays

!
! Lower boundary (W(K1)=B(K1):
    do i=lon0,lon1
      w1(lev0,i) = b(lev0,i) 
    enddo
!
! Set up work arrays:
    do i=lon0,lon1
      do k=k1+1,k2
!
! W(KF+K-1)=C(K-1)/W(K-1):
        w2(k-1,i) = c(k-1,i) / w1(k-1,i)
!
! W(K)=A(K)*W(KF+K-1)
        w1(k,i) = a(k,i) * w2(k-1,i)
!
! W(K)=B(K)-W(K)
        w1(k,i) = b(k,i) - w1(k,i)
      enddo ! k=k1+1,k2
    enddo ! i=lon0,lon1
!
! Lower boundary (W(2*KF+K1)=F(K1)/W(K1)):
    do i=lon0,lon1
      w3(k1,i) = f(k1,i) / w1(k1,i)
    enddo
!
    do i=lon0,lon1
      do k=k1+1,k2
!
! W(2*KF+K)=A(K)*W(2*KF+K-1)
        w3(k,i) = a(k,i) * w3(k-1,i)
!
! W(2*KF+K)=F(K)-W(2*KF+K)
        w3(k,i) = f(k,i) - w3(k,i)         
!
! W(2*KF+K)=W(2*KF+K)/W(K)
        w3(k,i) = w3(k,i) / w1(k,i)
      enddo ! k=k1+1,k2
    enddo ! i=lon0,lon1
!
! Upper boundary (X(K2)=W(2*KF+K2)):
    do i=lon0,lon1
      x(k2,i) = w3(k2,i)       
    enddo
!

! Back substitution:
    do i=lon0,lon1
      do kk=k1+1,k2          
        k = k1+k2-kk ! k2-1,k1,-1
!
! X(K)=W(KF+K)*X(K+1)
        x(k,i) = w2(k,i) * x(k+1,i)
!
! X(K)=W(2*KF+K)-X(K)
        x(k,i) = w3(k,i) - x(k,i)
      enddo ! k=k1+1,k2
    enddo
  end subroutine trsolv
!-----------------------------------------------------------------------
  subroutine printpoles(f,klev,k0,k1,i0,i1,j0,j1,name)
    use edyn_geogrid,only : nlat
!
! Args:
    integer,intent(in) :: klev,k0,k1,i0,i1,j0,j1
    real(r8),intent(in) :: f(k0:k1,i0:i1,j0:j1)
    character(len=*),intent(in) :: name
!
! Print values at the poles at klev:
    if (j0==1) then
      if (debug.and.masterproc) write(iulog,"(/,'printpoles ',a,' spole: klev=',i4,' f(klev,i0:i1,j0)=',/,(8es12.4))") &
        name,klev,f(klev,i0:i1,j0)
    endif
    if (j1==nlat) then
      if (debug.and.masterproc) write(iulog,"(/,'printpoles ',a,' npole: klev=',i4,' f(klev,i0:i1,j1)=',/,(8es12.4))") &
        name,klev,f(klev,i0:i1,j1)
    endif

  end subroutine printpoles
!-----------------------------------------------------------------------
  subroutine filter1_op(f,k0,k1,i0,i1,j0,j1)
!
! Polar fft filter, option 1 (see filter.F90).
!
    use vxbfilter_module,only: vxbfilter1,vxbkut1
    use edyn_mpi     ,only: mp_gatherlons_f3d,mytidi
    use edyn_mpi     ,only: mp_scatterlons_f3d
    use edyn_geogrid ,only: nlon
!
! Args:
    integer,intent(in) :: k0,k1,i0,i1,j0,j1
    real(r8),intent(inout) :: f(k0:k1,i0:i1,j0:j1)
!
! Local:
    integer :: i,j,k,nlevs
    real(r8) :: fik(nlon,k1-k0+1)
    type(array_ptr_type) :: fkij(1) ! fkij(1)%ptr(k1-k0+1,nlon,j0:j1)

    nlevs = k1-k0+1
!
! Define lons in fkij from current task subdomain:
!
    allocate(fkij(1)%ptr(nlevs,nlon,j0:j1))
    do j=j0,j1
      do i=i0,i1
        do k=k0,k1 ! kbot,nlev
          fkij(1)%ptr(k-k0+1,i,j) = f(k,i,j)
        enddo
      enddo
    enddo
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0
!   gather lons from other tasks in that row). This includes all latitudes.
!
    call mp_gatherlons_f3d(fkij,1,nlevs,i0,i1,j0,j1,1)
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
!
    if (mytidi==0) then
!
! Define 2d array with all longitudes for filter at each latitude:
!
      latscan: do j=j0,j1
        if (vxbkut1(j) >= nlon/2) cycle latscan
        do i=1,nlon
          do k=k0,k1
            fik(i,k-k0+1) = fkij(1)%ptr(k-k0+1,i,j)
          enddo
        enddo 
!
! Remove wave numbers > kut(lat):
!
        call vxbfilter1(fik,1,nlevs,j)
!
! Return filtered array to fkij:
!
        do i=1,nlon
          do k=k0,k1
            fkij(1)%ptr(k-k0+1,i,j) = fik(i,k-k0+1)
          enddo
        enddo ! i=1,nlon
      enddo latscan ! j=j0,j1
    endif ! mytidi==0
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
    call mp_scatterlons_f3d(fkij,1,nlevs,i0,i1,j0,j1,1)
!
! Return filtered array to inout field at task subdomain:
!
    do j=j0,j1
      do i=i0,i1
        do k=k0,k1
          f(k,i,j) = fkij(1)%ptr(k-k0+1,i,j)
        enddo
      enddo
    enddo
    deallocate(fkij(1)%ptr)
  end subroutine filter1_op
!-----------------------------------------------------------------------
  subroutine filter2_op(f,k0,k1,i0,i1,j0,j1)
    use vxbfilter_module,only: vxbfilter2
    use edyn_mpi     ,only: mp_gatherlons_f3d,mytidi
    use edyn_mpi     ,only: mp_scatterlons_f3d
    use edyn_geogrid ,only: nlon
!
! Args:
    integer,intent(in) :: k0,k1,i0,i1,j0,j1
    real(r8),intent(inout) :: f(k0:k1,i0:i1,j0:j1)
!
! Local:
    integer :: i,j,k,nlevs
    real(r8) :: fik(nlon,k1-k0+1)
    type(array_ptr_type) :: fkij(1) ! fkij(1)%ptr(k1-k0+1,nlon,j0:j1)

    nlevs = k1-k0+1
!
! Define lons in fkij from current task subdomain:
!
    allocate(fkij(1)%ptr(nlevs,nlon,j0:j1))
!$omp parallel do private( i,j,k )
    do j=j0,j1
      do i=i0,i1
        do k=k0,k1
          fkij(1)%ptr(k-k0+1,i,j) = f(k,i,j)
        enddo
      enddo
    enddo
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0
!   gather lons from other tasks in that row). This includes all latitudes.
!
    call mp_gatherlons_f3d(fkij,1,nlevs,i0,i1,j0,j1,1)
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
!
    if (mytidi==0) then
!
! Define 2d array with all longitudes for filter at each latitude:
!
      do j=j0,j1
        do i=1,nlon
          do k=k0,k1
            fik(i,k-k0+1) = fkij(1)%ptr(k-k0+1,i,j)
          enddo
        enddo 
!
! Remove wave numbers > kut(lat):
!
        call vxbfilter2(fik,1,nlevs,j)
!
! Return filtered array to fkij:
!
        do i=1,nlon
          do k=k0,k1
            fkij(1)%ptr(k-k0+1,i,j) = fik(i,k-k0+1)
          enddo
        enddo ! i=1,nlon
      enddo ! j=j0,j1
    endif ! mytidi==0
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
    call mp_scatterlons_f3d(fkij,1,nlevs,i0,i1,j0,j1,1)
!
! Return filtered array to inout field at task subdomain:
    do j=j0,j1
      do i=i0,i1
        do k=k0,k1
          f(k,i,j) = fkij(1)%ptr(k-k0+1,i,j)
        enddo
      enddo
    enddo
    deallocate(fkij(1)%ptr)
  end subroutine filter2_op
!-----------------------------------------------------------------------
end module vxb
