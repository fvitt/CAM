module vxb
!
! Transport ions (metal ions, O2+, NO+) by the VxB drifts and ambipolar
! diffusion on the WACCM-X geographic grid. Adapted from the O+ transport
! in oplus.F90 (B. Foster, May 2015), which was adapted from TIEGCM.
! Input ion densities come from dpie_coupling; transported densities
! (op_out and opnm_out in vxb_xport) are passed back to dpie_coupling.
!
  use shr_kind_mod   ,only: r8 => shr_kind_r8
  use cam_abortutils ,only: endrun
  use cam_logfile    ,only: iulog
  use edyn_geogrid   ,only: dphi,dlamda,cs, p0, nlev !, nlon, nlat, nlev
  use getapex        ,only: bx,by,bz,bmod2       ! (0:nlonp1,jspole-1:jnpole+1)
  use edyn_params    ,only: re=>Rearth ! Radius of Earth (cm)
  use time_manager   ,only: get_step_size
  use edyn_mpi       ,only: array_ptr_type
  use shr_const_mod  ,only: shr_const_g         ! gravitational constant (m/s^2)
  use spmd_utils     ,only: masterproc

  implicit none
  private
  public :: vxb_xport,  vxb_init
  public :: vxb_prep
  public :: vxb_ion_prep
  public :: kbot
  public :: ktop

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
    logical, parameter :: debug = .false.

    real(r8), allocatable :: expz(:) ! exp(-zp)
    real(r8), allocatable :: zp(:)   ! log pressure (as in tiegcm lev(nlev))
!
! Fields that do not depend on the transported ion species or the substep.
! They are set once per time step by vxb_prep and used by every vxb_xport call.
! Arrays with halo points are (nlev,i0-2:i1+2,j0-2:j1+2).
!
    real(r8), allocatable, target :: uui(:,:,:), vvi(:,:,:), wwi(:,:,:)       ! ion drifts (cm/s)
    real(r8), allocatable, target :: eeeex(:,:,:), eeeey(:,:,:), eeeez(:,:,:) ! E/B drift terms
    real(r8), allocatable, target :: uuii(:,:,:), vvii(:,:,:), wwii(:,:,:)    ! VnxB drifts from vxbvel
    real(r8), allocatable, target :: wn(:,:,:)  ! neutral vertical velocity, om*hj (cm/s)
    real(r8), allocatable, target :: tr(:,:,:)  ! Reduced temperature (.5*(tn+ti))
    real(r8), allocatable, target :: tp(:,:,:)  ! Plasma temperature (te+ti)
    real(r8), allocatable :: hj(:,:,:)   ! scale height (nlev,i0:i1,j0-1:j1+1)
    real(r8), allocatable :: dvb(:,:)    ! divergence of B-field (i0:i1,j0:j1)
    real(r8), allocatable :: aa(:,:,:)   ! cell horizontal area (nlev,i0:i1,j0:j1)
    real(r8), allocatable :: va(:,:,:)   ! cell volume (nlev,i0:i1,j0:j1)
    real(r8), allocatable :: rva(:,:,:)  ! 1/va (nlev,i0:i1,j0:j1)
!
! Fields that depend on the transported ion (through xi) but not on the substep.
! They are set once per ion by vxb_ion_prep. (nlev,i0-2:i1+2,j0-2:j1+2)
!
    real(r8), allocatable, target :: ri(:,:,:)     ! xi/(1+xi**2)
    real(r8), allocatable, target :: rri(:,:,:)    ! 1/(1+xi**2)
    real(r8), allocatable, target :: rrri(:,:,:)   ! xi**2/(1+xi**2)
    real(r8), allocatable, target :: rrrri(:,:,:)  ! 1/(1+xi**2)
    real(r8), allocatable, target :: bdotu(:,:,:)  ! rri*(B.U)*bz

  contains

!-----------------------------------------------------------------------
  subroutine vxb_init( adiff_limiter_in, shapiro_const_in )
!
! Set parameters, and the level indices and log-pressure coordinate of the
! transport region from the reference pressures. Stops if the vertical grid
! does not cover the transport region (pbot, pminor, ptop).
!
    use ref_pres, only: pref_mid

    real(r8), intent(in) :: adiff_limiter_in
    real(r8), intent(in) :: shapiro_const_in

    integer :: k
    real(r8) :: pmid(nlev)        ! reference pressure at midpoints, bottom to top (Pa)
    character(len=256) :: errmsg

    shapiro_const = shapiro_const_in
    adiff_limiter = adiff_limiter_in

    allocate(zp(nlev))      ! log pressure (as in TIEGCM)
    allocate(expz(nlev))    ! exp(-zp)

    pmid(1:nlev) = pref_mid(nlev:1:-1) ! invert ref pressure (Pa) as in tiegcm
!
! kbot is the k-index at the bottom of ion transport calculations,
!   corresponding to pressure pbot. ktop is the top, below ptop.
!   kminor is the bottom of the Coulomb term in rrk, below pminor.
!
    kbot = -999
    ktop = -999
    kminor = -999
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
!
! vxb_xport uses levels kbot-1 to ktop+2, and rrk needs kbot <= kminor <= ktop.
!
    if (kbot < 2) then
       write(errmsg,"('vxb_init: need a level with pressure <= pbot=',es10.3,' Pa above the lowest level; kbot=',i5)") &
            pbot, kbot
       call endrun(trim(errmsg))
    end if
    if (ktop < kbot) then
       write(errmsg,"('vxb_init: need a level with pressure <= ptop=',es10.3,' Pa above kbot=',i4,'; ktop=',i5)") &
            ptop, kbot, ktop
       call endrun(trim(errmsg))
    end if
    if (ktop+2 > nlev) then
       write(errmsg,"('vxb_init: need two levels above ktop=',i4,' (pressure <= ptop=',es10.3,' Pa); nlev=',i4)") &
            ktop, ptop, nlev
       call endrun(trim(errmsg))
    end if
    if (kminor < kbot .or. kminor > ktop) then
       write(errmsg,"('vxb_init: kminor=',i5,' (pressure <= pminor=',es10.3,' Pa) is outside kbot=',i4,' to ktop=',i4)") &
            kminor, pminor, kbot, ktop
       call endrun(trim(errmsg))
    end if

    do k=1,nlev
       zp(k) = -log(pmid(k)*10._r8/p0)
       expz(k) = exp(-zp(k))
    enddo
    dzp = zp(nlev)-zp(nlev-1)  ! use top 2 levels (typically dzp=0.5)

    if (masterproc) then
       write(iulog,"('vxb_init: shapiro=',es12.4,' dtsmooth=',es12.4,' dtsmooth_div2=',es12.4)") &
            shapiro_const,dtsmooth,dtsmooth_div2
       write(iulog,"('vxb_init: shr_const_g=',f8.3)") shr_const_g
       write(iulog,"('vxb_init: kbot=',i4,' (',es10.3,' Pa)  kminor=',i4,' (',es10.3,' Pa)  ktop=',i4,' (',es10.3,' Pa)')") &
            kbot, pmid(kbot), kminor, pmid(kminor), ktop, pmid(ktop)
       if (debug) then
          write(iulog,"('vxb_init: nlev=',i3,' zp (bot2top)   =',/,(6es12.3))") nlev,zp
          write(iulog,"('vxb_init: nlev=',i3,' expz (bot2top) =',/,(6es12.3))") nlev,expz
          write(iulog,"('vxb_init: nlev=',i3,' dzp  =',/,(6es12.3))") nlev,dzp
       end if
    end if

  end subroutine vxb_init

!-----------------------------------------------------------------------
  subroutine vxb_prep(tn,te,ti,un,vn,om,mbar,ui,vi,wi,eeex,eeey,eeez, &
                      i0,i1,j0,j1)
!
! Set the fields used by vxb_xport that do not depend on the transported
! ion species or the substep. Call once per time step, before the
! vxb_xport calls. Inputs are in "TIEGCM" format, as for vxb_xport.
!
    use edyn_mpi,     only: mp_geo_halos, mp_pole_halos
    use edyn_geogrid, only: nlat
!
! Args:
!
    integer, intent(in) :: i0  ! grid%ifirstxy
    integer, intent(in) :: i1  ! grid%ilastxy
    integer, intent(in) :: j0  ! grid%jfirstxy
    integer, intent(in) :: j1  ! grid%jlastxy
    real(r8), intent(in) :: tn  (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral temperature (deg K)
    real(r8), intent(in) :: te  (nlev,i0-2:i1+2,j0-2:j1+2) ! electron temperature (deg K)
    real(r8), intent(in) :: ti  (nlev,i0-2:i1+2,j0-2:j1+2) ! ion temperature (deg K)
    real(r8), intent(in) :: un  (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral zonal wind (cm/s)
    real(r8), intent(in) :: vn  (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral meridional wind (cm/s)
    real(r8), intent(in) :: om  (nlev,i0-2:i1+2,j0-2:j1+2) ! omega (1/s)
    real(r8), intent(in) :: mbar(nlev,i0-2:i1+2,j0-2:j1+2) ! mean molecular weight
    real(r8), intent(in) :: ui  (nlev,i0:i1,j0:j1)         ! zonal ion drift
    real(r8), intent(in) :: vi  (nlev,i0:i1,j0:j1)         ! meridional ion drift
    real(r8), intent(in) :: wi  (nlev,i0:i1,j0:j1)         ! vertical ion drift
    real(r8), intent(in) :: eeex(nlev,i0:i1,j0:j1)         ! zonal E/B drift term
    real(r8), intent(in) :: eeey(nlev,i0:i1,j0:j1)         ! meridional E/B drift term
    real(r8), intent(in) :: eeez(nlev,i0:i1,j0:j1)         ! vertical E/B drift term
!
! Local:
!
    integer :: i, j, k, lat0, lat1
    integer :: kb, kt   ! level range of tr and tp (one level beyond kbot,ktop)
    integer :: nfields
    real(r8), allocatable :: polesign(:)
    type(array_ptr_type), allocatable :: ptrs(:)
    real(r8), dimension(nlev,i0:i1,j0:j1) :: uii, vii, wii ! VnxB drifts from vxbvel

    if (.not. allocated(zp)) then
       call endrun('vxb_prep: vxb_init must be called before vxb_prep')
    end if

    if (.not. allocated(tr)) then
       allocate(uui  (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(vvi  (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(wwi  (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(eeeex(nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(eeeey(nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(eeeez(nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(uuii (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(vvii (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(wwii (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(wn   (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(tr   (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(tp   (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(hj   (nlev,i0:i1,j0-1:j1+1))
       allocate(dvb  (i0:i1,j0:j1))
       allocate(aa   (nlev,i0:i1,j0:j1))
       allocate(va   (nlev,i0:i1,j0:j1))
       allocate(rva  (nlev,i0:i1,j0:j1))
       allocate(ri   (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(rri  (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(rrri (nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(rrrri(nlev,i0-2:i1+2,j0-2:j1+2))
       allocate(bdotu(nlev,i0-2:i1+2,j0-2:j1+2))
    end if

    grav_cm = shr_const_g * 100._r8 ! m/s^2 -> cm/s^2
!
! Ion drifts and E/B terms with halo points:
!
!$omp parallel do private(i, j, k)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          uui(k,i,j)=ui(k,i,j)
          vvi(k,i,j)=vi(k,i,j)
          wwi(k,i,j)=wi(k,i,j)
          eeeex(k,i,j)=eeex(k,i,j)
          eeeey(k,i,j)=eeey(k,i,j)
          eeeez(k,i,j)=eeez(k,i,j)
        enddo
      enddo
    enddo

    nfields = 6
    allocate(ptrs(nfields),polesign(nfields))
    ptrs(1)%ptr => uui ; ptrs(2)%ptr => vvi ; ptrs(3)%ptr => wwi
    ptrs(4)%ptr => eeeex ; ptrs(5)%ptr => eeeey ; ptrs(6)%ptr => eeeez
    polesign = 1._r8
    polesign(1:2) = -1._r8 ! uui,vvi
    polesign(4:5) = -1._r8 ! eeeex,eeeey

    call mp_geo_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields)
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
! Divergence of B (mag field) is returned by divb in dvb(i0:i1,j0:j1)
!
    call divb(dvb,i0,i1,j0,j1)

    tr  = 0._r8
    tp  = 0._r8
    hj  = 0._r8
    kb = kbot-1
    kt = ktop+1
!
! The latitude scans in vxb_xport use rows lat-1:lat+1 of tr and hj, and
! lat-2:lat+2 of tp (the outer rows come from the halo exchange below).
! Each row j=lat0-1,lat1+1 is set once, so the rows are independent.
!
    if (lat0 <= lat1) then
!$omp parallel do private(j, i, k)
      do j=lat0-1,lat1+1
        do i=i0,i1
!
! Reduced temperature (tpj in tiegcm):
! 'OPLUS_TR' (has constants at poles)
!
          do k=kb,kt
            tr(k,i,j) = 0.5_r8*(tn(k,i,j)+ti(k,i,j))
          enddo
!
! Plasma temperature:
! 'OPLUS_TP0' (tp gets poles from rows lat0-1 and lat1+1)
!
          do k=kb,kt
            tp(k,i,j) = te(k,i,j)+ti(k,i,j)
          enddo
!
! Neutral scale height:
! 'OPLUS_HJ' (has constants at poles)
!
          do k=1,nlev
            hj(k,i,j) = gask * tn(k,i,j) / (mbar(k,i,j) * grav_cm)
          enddo
        enddo ! i=i0,i1
      enddo ! j=lat0-1,lat1+1
    endif

    call vxbvel(un,vn,om,hj,uii,vii,wii,1,nlev,i0,i1,j0,j1)
!$omp parallel do private(k, j, i)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          uuii(k,i,j)=uii(k,i,j)
          vvii(k,i,j)=vii(k,i,j)
          wwii(k,i,j)=wii(k,i,j)
          wn(k,i,j)=om(k,i,j)*hj(k,i,j)
        enddo
      enddo
    enddo
!
! Set halo points where needed.
!
    nfields = 5
    allocate(ptrs(nfields),polesign(nfields))
    ptrs(1)%ptr => tp
    ptrs(2)%ptr => uuii ; ptrs(3)%ptr => vvii
    ptrs(4)%ptr => wwii ; ptrs(5)%ptr => wn
    polesign = 1._r8
    polesign(2:3) = -1._r8 ! uuii,vvii

    call mp_geo_halos (ptrs,1,nlev,i0,i1,j0,j1,nfields)
    call mp_pole_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields,polesign)
    deallocate(ptrs,polesign)
!
! Cell area and volume for the vertical advection scheme:
!
!$omp parallel do private(k, j, i)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          aa(k,i,j)=re*re*dphi*dlamda*cs(j)
          va(k,i,j)=aa(k,i,j)*hj(k,i,j)*dzp
          rva(k,i,j)=1._r8/va(k,i,j)
        enddo
      enddo
    enddo

  end subroutine vxb_prep

!-----------------------------------------------------------------------
  subroutine vxb_ion_prep(xi,un,vn,i0,i1,j0,j1)
!
! Set the fields used by vxb_xport that depend on the transported ion
! (through xi) but not on the substep. Call once per ion, after vxb_prep
! and before that ion's vxb_xport calls.
!
    use edyn_mpi, only: mp_geo_halos, mp_pole_halos
!
! Args:
!
    integer, intent(in) :: i0  ! grid%ifirstxy
    integer, intent(in) :: i1  ! grid%ilastxy
    integer, intent(in) :: j0  ! grid%jfirstxy
    integer, intent(in) :: j1  ! grid%jlastxy
    real(r8), intent(in) :: xi(nlev,i0:i1,j0:j1)          ! ion-neutral collision freq / ion gyrofreq
    real(r8), intent(in) :: un(nlev,i0-2:i1+2,j0-2:j1+2)  ! neutral zonal wind (cm/s)
    real(r8), intent(in) :: vn(nlev,i0-2:i1+2,j0-2:j1+2)  ! neutral meridional wind (cm/s)
!
! Local:
!
    integer :: i, j, k
    integer :: nfields
    real(r8), allocatable :: polesign(:)
    type(array_ptr_type), allocatable :: ptrs(:)

    if (.not. allocated(tr)) then
       call endrun('vxb_ion_prep: vxb_prep must be called before vxb_ion_prep')
    endif

!$omp parallel do private(i, j, k)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          ri(k,i,j) = xi(k,i,j)/(1._r8+xi(k,i,j)**2)
          rrri(k,i,j) = xi(k,i,j)**2/(1._r8+xi(k,i,j)**2)
          rrrri(k,i,j) = 1._r8/(1._r8+xi(k,i,j)**2)
          rri(k,i,j) = 1._r8/(1._r8+xi(k,i,j)**2)
        enddo
      enddo
    enddo

    bdotu = 0._r8
!$omp parallel do private(k, j, i)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          bdotu(k,i,j) = rri(k,i,j)*(bx(i,j)*un(k,i,j)+by(i,j)*vn(k,i,j)+ &
            wn(k,i,j)*bz(i,j))*bz(i,j)
        enddo
      enddo
    enddo

    nfields = 5
    allocate(ptrs(nfields),polesign(nfields))
    ptrs(1)%ptr => ri ; ptrs(2)%ptr => rri
    ptrs(3)%ptr => rrri ; ptrs(4)%ptr => rrrri
    ptrs(5)%ptr => bdotu
    polesign = 1._r8

    call mp_geo_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields)
    call mp_pole_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields,polesign)
    deallocate(ptrs,polesign)

  end subroutine vxb_ion_prep

!-----------------------------------------------------------------------
  subroutine vxb_xport(tn,te,ti,ne,nop,un,vn,om,o2,o1,n2,op_in,opnm_in, &
                         mbar,op_out,opnm_out,rmass_ion, &
                         i0,i1,j0,j1,nspltop )
!
! All input fields from dpie_coupling are in "TIEGCM" format, i.e.,
! longitude (-180->180), vertical (bot2top), and units (CGS).
!
! vxb_prep must be called earlier in the time step to set the fields
! that do not depend on the ion species (ion drifts, tp, hj, ...), and
! vxb_ion_prep before the first substep of each ion (ri, rri, bdotu, ...).
!
    use edyn_mpi,only:  mp_geo_halos,mp_pole_halos,setpoles
    use edyn_geogrid,only : nlat, nlev
!
! Transport one ion species over one substep. Adapted from the O+ transport
! in oplus.F90 (March-May 2015, B. Foster), which was adapted from TIEGCM (oplus.F).
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
    integer,intent(in) :: nspltop
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

    real(r8),intent(in) :: op_in(nlev,i0:i1,j0:j1) ! ion density (cm^3)
    real(r8),intent(in) :: opnm_in(nlev,i0:i1,j0:j1) ! ion density (cm^3) at time-1
    real(r8),intent(in) :: rmass_ion  ! ion molecular weight (g/mol)
!
! Output:
!
    real(r8),intent(out) ::       &
      op_out  (nlev,i0:i1,j0:j1), & ! ion density output (cm^3)
      opnm_out(nlev,i0:i1,j0:j1)    ! ion density output at time n-1 (cm^3)
!
! Local:
!
    integer :: i,j,k,lat,jm1,jp1,lat0,lat1
!
! Local inputs with added halo points in lat,lon:
! (ri, rri, rrri, rrrri and bdotu are module arrays set by vxb_ion_prep.)
!
    real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: op, opnm,xiop,xiopnm, &
        xxiop,xxiopnm,xxxiop,xxxiopnm,xxxxiop,xxxxiopnm


    real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: &
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
        rwii,              &! metal ion vertical velocity  Jianfei Wu
        hadudth,epsh,scrh, &
        nulh, mulh, flxh, diff, LORHOT, LNRHOT, RHOT, RHOTD, FABS, FSGN, TERM, TERP
    real(r8) :: dtime, ONE6TH,ONE3RD
    real(r8),dimension(nlev,i0:i1) :: &
      bdzdvb_op
!
! Explicit (horizontal transport) terms (no halos):
    real(r8),dimension(nlev,i0:i1,j0:j1) :: &
      explicit,explicit_a1,explicit_b1,&
      explicit3,explicit_a2,explicit_b2,&
      explicit4,explicit_a3,explicit_b3,&
      explicit5,explicit_a4,explicit_b4,&
      explicit6
!
! Pointers for multiple-field calls (e.g., mp_geo_halos)
    integer :: nfields
    real(r8),allocatable :: polesign(:)
    type(array_ptr_type),allocatable :: ptrs(:)

    Real(r8),Parameter :: BIGNUM = Huge(1._r8)
!
! Execute:
!
    dtime = get_step_size() ! step size in seconds
    dtime = dtime / dble(nspltop)

    if (.not. allocated(tr)) then
       call endrun('vxb_xport: vxb_prep must be called before vxb_xport')
    endif
!
! Set subdomain blocks from input (composition is in mmr):
!
!$omp parallel do private(i, j, k)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
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
    nfields = 10
    allocate(ptrs(nfields),polesign(nfields))

    ptrs(1)%ptr => op ; ptrs(2)%ptr => opnm
    ptrs(3)%ptr =>xiop ; ptrs(4)%ptr => xiopnm
    ptrs(5)%ptr =>xxiop ; ptrs(6)%ptr => xxiopnm
    ptrs(7)%ptr =>xxxiop ; ptrs(8)%ptr => xxxiopnm
    ptrs(9)%ptr =>xxxxiop ; ptrs(10)%ptr => xxxxiopnm
    polesign = 1._r8
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
    tpp  = 0._r8
    dj  = 0._r8
    dji  = 0._r8
    bvel= 0._r8
    diffj = 0._r8
    diffji = 0._r8
    diffjj = 0._r8
    abk = 1._r8
    opnm_smooth = 0._r8
    hadudth= 0._r8
    epsh=0._r8
    scrh=0._r8
    nulh=0._r8
    mulh=0._r8
    flxh=0._r8
    diff=0._r8
    ktop=ktop+1
    kbot=kbot-1
!
!----------------------- Begin first latitude scan ---------------------
!
! The scans at lat=lat0,lat1 use dj, dji, bvel, diffj, diffji and diffjj at
! rows lat-1:lat+1, and tpp at rows lat-2:lat+2. Each row is computed once
! here, so the rows are independent and can be done in parallel.
!
! as of April, 2015, TIEGCM incorrectly uses te+ti instead of tn+ti
! This has not been fixed in TIEGCM, because fixing it causes a tuning
! problem (ask Hanli and Wenbin). For WACCM, it is correct as below.
! (tr, tp and hj are set in vxb_prep.)
!
    if (lat0 <= lat1) then
!$omp parallel do private(j, i, k)
      do j=lat0-1,lat1+1
!
! rrk returns ambipolar diffusion coefficients in dj and dji:
! 'OPLUS_DJ' (has constants at poles)
!
        call rrk(                                            &
          tn(kbot:ktop,i0:i1,j),ti(kbot:ktop,i0:i1,j),mbar(kbot:ktop,i0:i1,j), &
          o2(kbot:ktop,i0:i1,j),o1  (kbot:ktop,i0:i1,j), &
          n2(kbot:ktop,i0:i1,j),nop(kbot:ktop,i0:i1,j),op(kbot:ktop,i0:i1,j), &
          rmass_ion,tr(kbot:ktop,i0:i1,j), &
          dj(kbot:ktop,i0:i1,j),dji(kbot:ktop,i0:i1,j),i0,i1,kbot,ktop,kminor)
!
! bvel = (B.U)*N(O+)
! 'OPLUS_BVEL' (has constants at poles)
!
! Note bx,by,bz were set globally for all tasks by sub magfield
! (getapex.F90)
!
        do i=i0,i1
          do k=kbot,ktop
            bvel(k,i,j) = &
              rri(k,i,j)*(bx(i,j)*un(k,i,j)+by(i,j)*vn(k,i,j)+   &
              hj(k,i,j)*bz(i,j)*om(k,i,j))*op(k,i,j)
          enddo ! k=kbot,ktop
        enddo ! i=lon0,lon1
!
! Ambipolar diffusion is returned in diffj:
! 'OPLUS_DIFFJ' (has constants at poles)
!
        call diffus(tp(kbot:ktop,i0:i1,j),ti(kbot:ktop,i0:i1,j),te(kbot:ktop,i0:i1,j),&
            nop(kbot:ktop,i0:i1,j),ne(kbot:ktop,i0:i1,j),&
            op(kbot:ktop,i0:i1,j),&
            xxxxiop(kbot:ktop,i0:i1,j),rmass_ion,hj(kbot:ktop,:,j), &
            diffj(kbot:ktop,i0:i1,j),diffji(kbot:ktop,i0:i1,j),i0,i1,kbot,ktop,j)
        call diffusj(tp(kbot:ktop,i0:i1,j),ti(kbot:ktop,i0:i1,j),te(kbot:ktop,i0:i1,j),&
            nop(kbot:ktop,i0:i1,j),ne(kbot:ktop,i0:i1,j),&
            op(kbot:ktop,i0:i1,j),rmass_ion,hj(kbot:ktop,:,j),dj(kbot:ktop,i0:i1,j), &
            dji(kbot:ktop,i0:i1,j),diffjj(kbot:ktop,i0:i1,j),i0,i1,kbot,ktop,j)
      enddo ! j=lat0-1,lat1+1
!
! 'OPLUS_TP1' (constants at the poles)
!
!$omp parallel do private(j, i, k)
      do j=lat0-2,lat1+2
        do i=i0,i1
          do k=kbot,ktop
            tpp(k,i,j) = ti(k,i,j)*xxxxiop(k,i,j)
          enddo
        enddo
      enddo ! j=lat0-2,lat1+2
    endif
!
! Latidinal shapiro smoother: opnm is O+ at time n-1.
! opnm_smooth will be used in explicit terms below.
! Smooth in latitude:
! 'OPNM_SMOOTH' (zero at poles)
!
!$omp parallel do private(lat, i, k)
    do lat=lat0,lat1
      do i=i0,i1
        do k=kbot,ktop
          opnm_smooth(k,i,lat) = opnm(k,i,lat)-shapiro_const* &
                                (opnm(k,i,lat+2)+opnm(k,i,lat-2)-4._r8*     &
                                (opnm(k,i,lat+1)+opnm(k,i,lat-1))+6._r8*    &
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

!
! Set halo points where needed (tp,uuii,vvii,wwii,wn halos are set in vxb_prep,
! bdotu halos in vxb_ion_prep).
!
    nfields = 7
    allocate(ptrs(nfields),polesign(nfields))
    ptrs(1)%ptr => dj ; ptrs(2)%ptr => bvel ; ptrs(3)%ptr => diffj
    ptrs(4)%ptr => opnm_smooth
    ptrs(5)%ptr => diffjj
    ptrs(6)%ptr => diffji
    ptrs(7)%ptr => tpp
    polesign = 1._r8


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

!$omp parallel do private(lat, i, k)
    do lat=lat0,lat1
!
! bdotdh_op = (B(H).DEL(H))*(D/(H*DZ)*TP+M*G/R)*N(O+)
! then bdotdh_op = d*bz*bdotdh_op
! real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2) :: diffj
! real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2) :: bdotdh_op
! 'BDOTDH_OP' (zero at the poles)
!
      call bdotdh(                    &
        diffj(kbot:ktop,i0:i1,lat-1), &
        diffj(kbot:ktop,:,lat    ),   & ! includes longitude halos
        diffj(kbot:ktop,i0:i1,lat+1), &
        bdotdh_op(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat)
      call bdotdh(                    &
        diffji(kbot:ktop,i0:i1,lat-1), &
        diffji(kbot:ktop,:,lat    ),   & ! includes longitude halos
        diffji(kbot:ktop,i0:i1,lat+1), &
        bdotdh_opi(kbot:ktop,i0:i1,lat),i0,i1,kbot,ktop,lat)
!
      do i=i0,i1
        do k=kbot,ktop
          bdotdh_op(k,i,lat) = bz(i,lat)*(dj(k,i,lat)*bdotdh_op(k,i,lat)+&
              dji(k,i,lat)*bdotdh_opi(k,i,lat)) ! BDOTDH_OP
          ambvel(k,i,lat) = (bz(i,lat)**2)*diffjj(k,i,lat) ! BDOTDH_OP
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
    enddo ! lat=lat0,lat1
!
! bdotdh_opj = (B(H).DEL(H))*2.*TP*N(O+) is needed at rows lat-1:lat+1.
! Each row j=lat0-1,lat1+1 is computed once, so the rows are independent.
! 'BDOTDH_OPJ' (has reasonable non-constant values at poles)
!
    if (lat0 <= lat1) then
!$omp parallel do private(j, i, k)
      do j=lat0-1,lat1+1
        call bdotdh(                &
          tpp(kbot:ktop,i0:i1,j-1),  &
          tpp(kbot:ktop,:,j),        &
          tpp(kbot:ktop,i0:i1,j+1),  &
          bdotdh_opj(kbot:ktop,i0:i1,j),i0,i1,kbot,ktop,j)

        call firstbdotdh(                &
          tp(kbot:ktop,i0:i1,j-1),  &
          tp(kbot:ktop,:,j),        &
          tp(kbot:ktop,i0:i1,j+1),  &
          nop(kbot:ktop,i0:i1,j-1),  &
          nop(kbot:ktop,:,j),        &
          nop(kbot:ktop,i0:i1,j+1),  &
          xxxxiop(kbot:ktop,i0:i1,j),  &
          bdotdh_opj1(kbot:ktop,i0:i1,j),i0,i1,kbot,ktop,j)

        call firstbdotdh(                &
          te(kbot:ktop,i0:i1,j-1),  &
          te(kbot:ktop,:,j),        &
          te(kbot:ktop,i0:i1,j+1),  &
          ne(kbot:ktop,i0:i1,j-1),  &
          ne(kbot:ktop,:,j),        &
          ne(kbot:ktop,i0:i1,j+1),  &
          xxxxiop(kbot:ktop,i0:i1,j),  &
          bdotdh_opj2(kbot:ktop,i0:i1,j),i0,i1,kbot,ktop,j)
!
        do i=i0,i1
          do k=kbot,ktop
            bdotdh_opj(k,i,j) = bdotdh_opj(k,i,j)*dji(k,i,j)+&
                bdotdh_opj1(k,i,j)*dj(k,i,j)+&
                bdotdh_opj2(k,i,j)*dji(k,i,j)
          enddo ! k=kbot,nlev
        enddo ! i=i0,i1
      enddo ! j=lat0-1,lat1+1 (end second lat scan)
    endif
!
!------------------------ End second latitude scan ---------------------
!
! bdotdh_opj already has non-constant polar values. The poles of bdotdh_op
! and ambvel are zero, but they are only used at lat=lat0,lat1.
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
!----------------------- Begin third latitude scan ---------------------
!
    bdotdh_diff = 0._r8
    explicit(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a1(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b1(1:nlev,i0:i1,j0:j1)=0._r8
    explicit3(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a2(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b2(1:nlev,i0:i1,j0:j1)=0._r8
    explicit4(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a3(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b3(1:nlev,i0:i1,j0:j1)=0._r8
    explicit5(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a4(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b4(1:nlev,i0:i1,j0:j1)=0._r8
    explicit6(1:nlev,i0:i1,j0:j1)    = 0._r8
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
! Globally, this loop is lat=2,nlat-1 (i.e., skipping the poles)
! Each latitude writes only row lat (bdzdvb_op is a per-latitude work array).
!
!$omp parallel do private(lat, jm1, jp1, i, k, bdzdvb_op)
    do lat=lat0,lat1
      jm1 = lat-1 ! this will be south pole for southern pes (j==1)
      jp1 = lat+1 ! this will be north pole for northern pes (j==nlat)
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
      do i=i0,i1
        do k=kbot,ktop
          explicit(k,i,lat) = -(1._r8)*(bdzdvb_op(k,i)+bdotdh_diff(k,i,lat)+ &
            bdotdh_op(k,i,lat))
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
!
         explicit(k,i,lat) = explicit(k,i,lat)+ &
             explicit3(k,i,lat)+  &
             explicit4(k,i,lat)+  &
             explicit5(k,i,lat)+  &
             explicit6(k,i,lat)
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
               write(iulog,"('  op(k,i-1:i+1,lat)   =',3es12.4)") op(k,i-1:i+1,lat)
               write(iulog,"('  op(k,i,jm1:jp1)     =',3es12.4)") op(k,i,jm1:jp1)
               write(iulog,"('  bmod2(i-1:i+1,lat) =',3es12.4)") bmod2(i-1:i+1,lat)
               write(iulog,"('  bmod2(i,jm1:jp1)   =',3es12.4)") bmod2(i,jm1:jp1)
               write(iulog,"('  vxb_ui(k:k+1,i,lat)    =',2es12.4)") uui(k:k+1,i,lat)
               write(iulog,"('  vxb_vi(k:k+1,i,lat)    =',2es12.4)") vvi(k:k+1,i,lat)
               write(iulog,"('  bx,by(i,lat)       =',2es12.4)") bx(i,lat),by(i,lat)
            endif
         endif

        enddo ! k=kbot,nlev-1
      enddo ! i=i0,i1

      do i=i0,i1
        do k=kbot,ktop
          rwii(k,i,lat)=0.5_r8*(wwii(k,i,lat)*ri(k,i,lat)+wwii(k+1,i,lat)*ri(k+1,i,lat))&
              +0.5_r8*(wwi(k,i,lat)*rri(k,i,lat)+wwi(k+1,i,lat)*rri(k+1,i,lat))&
              +0.5_r8*(eeeez(k,i,lat)*ri(k,i,lat)+eeeez(k+1,i,lat)*ri(k+1,i,lat))&
              +0.5_r8*(bdotu(k,i,lat)+bdotu(k+1,i,lat))&
              -0.5_r8*(rrrri(k,i,lat)*(ambvel(k,i,lat))+&
              rrrri(k+1,i,lat)*(ambvel(k+1,i,lat)))
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
    enddo ! end third latitude scan (lat=lat0,lat1)
!
!------------------------ End third latitude scan ---------------------
    ktop=ktop-1
    kbot=kbot+1
    ONE6TH = 1.0_r8/6.0_r8
    ONE3RD = 1.0_r8/3.0_r8
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
!----------------------- Begin fourth latitude scan ---------------------
!
!$omp parallel do private(lat, i, k)
    do lat=j0,j1
      do i=i0,i1
        do k=kbot,ktop
          opnm_out(k,i,lat) = dtsmooth*op(k,i,lat)+dtsmooth_div2* &
            (opnm(k,i,lat)+op_out(k,i,lat))
        enddo
      enddo
!
! Insure non-negative ion output:
      do i=i0,i1
        do k=kbot,ktop
          if (op_out  (k,i,lat) < 1.e-7_r8) op_out  (k,i,lat) = 1.e-7_r8
          if (opnm_out(k,i,lat) < 1.e-7_r8) opnm_out(k,i,lat) = 1.e-7_r8
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
    enddo ! lat=j0,j1
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
  subroutine rrk(t,ti,rms,ps1,ps2,n2,nop,en,rmass_ion,tr,ans,ans1,lon0,lon1,lev0,lev1,levn)
!
! Returns ambipolar diffusion coefficient in ans.
!
! Args:
    integer,intent(in) :: lon0,lon1,lev0,lev1,levn
    real(r8),intent(in) :: rmass_ion
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
    do i=lon0,lon1
      do k=lev0,lev1-1
        ans(k,i) = 1.42e17_r8*boltz*t(k,i)/(p0*expz(k)*.5_r8*(rms(k,i)+    &
          rms(k+1,i))*(ps2(k,i)*rmassinv_o1*sqrt(tr(k,i))*(1._r8-0.064_r8* &
          log10(tr(k,i)))**2*colfac+18.6_r8*n2(k,i)*rmassinv_n2+18.1_r8*   &
          ps1(k,i)*rmassinv_o2))
        ans1(k,i) = 0._r8
      enddo ! k=lev0,lev1
      ans(lev1,i) = ans(lev1-1,i) ! should not need to do this
      ans1(lev1,i) = ans1(lev1-1,i) ! should not need to do this

    enddo ! i=lon0,lon1
    do i=lon0,lon1
      do k=levn,lev1-1
        ans1(k,i) = gask/(rmass_ion*0.0851_r8*sqrt((rmass_ion+16._r8)/(rmass_ion*16._r8))* &
            ti(k,i)**(-1.5_r8)*nopn(k,i)*log(10518._r8*ti(k,i)**1.5_r8/sqrt(nopn(k,i))))
      enddo ! k=lev0,lev1
      ans1(lev1,i) = ans1(lev1-1,i) ! should not need to do this

    enddo ! i=lon0,lon1
      do i=lon0,lon1
        do k=lev0,lev1
          if (nopn  (k,i) < 2.e4_r8) then
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
  subroutine diffus(tp,ti,te,nop,ne,en,een,rmass_ion,hj,ans,ans1,i0,i1,lev0,lev1,lat)
!                                      kbot,nlev
! Evaluates ans = (d/(h*dz)*tp+m*g/r)*en
! Remember: "bot2top": lev0=kbot=bottom, lev1=nlev=top
!
! Args:
    integer :: i0,i1,lev0,lev1,lat
    real(r8),intent(in) :: rmass_ion
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

    mgr = rmass_ion*grav_cm/gask
    mgro = 16._r8*grav_cm/gask

    do i=i0,i1
      do k=lev0,lev1-2
        ans(k+1,i) =1.0_r8/(2._r8*hj(k+1,i)*dzp)*&
        (nopn(k+2,i)*tp(k+2,i)-nopn(k,i)*tp(k,i))/nopn(k+1,i)*een(k+1,i)&
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
    do i=i0,i1
!
! Upper boundary:
      ans(lev1,i) = 0._r8
      ans1(lev1,i) = 0._r8
!
! Lower boundary:
      ans(lev0,i) = 1._r8/(hj(lev0,i)*dzp)*(tp(lev0+1,i)*nopn(lev0+1,i)- &
        tp(lev0,i)*nopn(lev0,i))/nopn(lev0,i)*een(lev0,i)+mgro*een(lev0,i)
      ans1(lev0,i) = 1._r8/(hj(lev0,i)*dzp)*&
          (((ti(lev0+1,i)*enn(lev0+1,i)-ti(lev0,i)*enn(lev0,i))/enn(lev0,i)+&
          (te(lev0+1,i)*ne(lev0+1,i)-te(lev0,i)*ne(lev0,i))/ne(lev0,i))*een(lev0,i))+mgr*een(lev0,i)
    enddo
      do i=i0,i1
        do k=lev0,lev1
          if (nopn  (k,i) < 2.e4_r8) then
              ans1  (k,i) = 0._r8
          endif
          if (nopn  (k,i) < 1.e3_r8) then
              ans  (k,i) = 0._r8
          endif
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
  end subroutine diffus
!-----------------------------------------------------------------------
  subroutine diffusj(tp,ti,te,nop,ne,en,rmass_ion,hj,dj,dji,ans,i0,i1,lev0,lev1,lat)
!                                      kbot,nlev
! Evaluates ans = (d/(h*dz)*tp+m*g/r)*en
! Remember: "bot2top": lev0=kbot=bottom, lev1=nlev=top
!
! Args:
    integer :: i0,i1,lev0,lev1,lat
    real(r8),intent(in) :: rmass_ion
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

    mgr = rmass_ion*grav_cm/gask
    mgro = 16._r8*grav_cm/gask
!
! diffusj is called from inside the threaded loop in vxb_xport, so this
! nested region normally runs on one thread. The directive is kept because
! removing it changes which multiply-add the compiler fuses into an FMA
! below, which changes answers at roundoff level (seen with gfortran -O3).
!
!$omp parallel do private(i,k)
    do i=i0,i1
      do k=lev0,lev1-2
        ans1(k+1,i) =1.0_r8/(2._r8*hj(k+1,i)*dzp)*&
            ((nopn(k+2,i)*tp(k+2,i)- nopn(k,i)*tp(k,i))/nopn(k+1,i))&
          *dj(k+1,i)+mgro*dj(k+1,i)
          if (nopn  (k,i) < 1.e3_r8) then
              ans1  (k,i) = 0._r8
          endif
        ans2(k+1,i) = 1.0_r8/(2._r8*hj(k+1,i)*dzp)*((enn(k+2,i)*ti(k+2,i)- &
          enn(k,i)*ti(k,i))/(enn(k+1,i))+(ne(k+2,i)*te(k+2,i)-ne(k,i)*te(k,i))/(ne(k+1,i)))&
          *dji(k+1,i)+mgr*dji(k+1,i)
          if (enn  (k,i) < 1.e-4_r8) then
              ans2  (k,i) = 0._r8
          endif
        ans(k+1,i)=ans1(k+1,i)+ans2(k+1,i)
      enddo
      if (debug) then
        write(iulog,"('diffus: lat=',i4,' i=',i4,' ans(lev0:lev1-1,i)=',2es12.4)") &
          lat,i,minval(ans(lev0:lev1-1,i)),maxval(ans(lev0:lev1-1,i))
      endif
    enddo

! Upper and lower boundaries:
!
    do i=i0,i1
!
! Upper boundary:
        ans(lev1,i) = 0.0_r8
!
! Lower boundary:
        ans1(lev0,i) =1.0_r8/(hj(lev0,i)*dzp)*&
            ((nopn(lev0+1,i)*tp(lev0+1,i)-nopn(lev0,i)*tp(lev0,i))/nopn(lev0,i))&
          *dj(lev0,i)+mgro*dj(lev0,i)
          if (nopn  (lev0,i) < 1.e3_r8) then
              ans1  (lev0,i) = 0._r8
          endif
        ans2(lev0,i) = 1.0_r8/(hj(lev0,i)*dzp)*((enn(lev0+1,i)*ti(lev0+1,i)- &
          enn(lev0,i)*ti(lev0,i))/(enn(lev0,i))+(ne(lev0+1,i)*te(lev0+1,i)-ne(lev0,i)*te(lev0,i))/(ne(lev0,i)))&
          *dji(lev0,i)+mgr*dji(lev0,i)
          if (enn  (lev0,i) < 1.e-4_r8) then
              ans2  (lev0,i) = 0._r8
          endif
        ans(lev0,i)=ans1(lev0,i)+ans2(lev0,i)
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
    real(r8),dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: phij
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out) :: ans
!
! Local:
    integer :: k,i
!
! Note phij longitude dimension is lon0-2:lon1+2 (only i-1 and i+1 are used).
! Halo longitudes i-1 and i+1 must have been set before this routine is
! called. ('by' is use-associated above)
!
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
    do i=lon0,lon1
      do k=lev0+1,lev1-1
        ans(k,i) = bz(i,lat)/(2._r8*h(k,i)*dzp)*(phi(k+1,i)-phi(k-1,i))+ &
          dvb(i)*phi(k,i)
      enddo ! k=lev0+1,lev1-1
    enddo ! i=lon0,lon1
!
! Upper and lower boundaries:
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
end module vxb
