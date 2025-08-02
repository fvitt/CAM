module edyn3d_esmf_fields_rhandles
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc, mpicom
  use edyn3d_esmf_error_mod, only: check_error

  use ESMF, only: ESMF_KIND_I4, ESMF_KIND_R8, ESMF_TYPEKIND_R8, ESMF_MESHLOC_ELEMENT, ESMF_STAGGERLOC_CENTER
  use ESMF, only: ESMF_Field, ESMF_RouteHandle, ESMF_ArraySpec, ESMF_ArraySpecSet, ESMF_FieldCreate
  use ESMF, only: ESMF_FieldRegridStore, ESMF_REGRIDMETHOD_BILINEAR, ESMF_POLEMETHOD_ALLAVG, ESMF_EXTRAPMETHOD_NEAREST_IDAVG
  use ESMF, only: ESMF_FieldDestroy, ESMF_RouteHandleDestroy

  implicit none

  private

  public :: physFieldSrc
  public :: magFieldDes_s1
  public :: magFieldDes_s2
  public :: magFieldSrc_s2
  public :: oplusFieldDes

  public :: magFieldSrc_ref_p
  public :: oplusFieldDes_ref_p

  public :: rh_mag2oplus_s2
  public :: rh_phys2mag_s1
  public :: rh_phys2mag_s2
  public :: rh_mag2oplus_ref_p

  public :: edyn3d_esmf_fields_rhandles_init
  public :: edyn3d_esmf_fields_rhandles_destroy

  type(ESMF_Field) :: physFieldSrc
  type(ESMF_Field) :: oplusFieldDes

  type(ESMF_Field) :: magFieldSrc_ref_p
  type(ESMF_Field) :: oplusFieldDes_ref_p

  type(ESMF_Field), allocatable :: magFieldDes_s1(:)
  type(ESMF_Field), allocatable :: magFieldDes_s2(:)
  type(ESMF_Field), allocatable :: magFieldSrc_s2(:)

  type(ESMF_RouteHandle), allocatable :: rh_phys2mag_s1(:)
  type(ESMF_RouteHandle), allocatable :: rh_phys2mag_s2(:)
  type(ESMF_RouteHandle), allocatable :: rh_mag2oplus_s2(:)
  type(ESMF_RouteHandle) :: rh_mag2oplus_ref_p

  integer, public, parameter :: phys2mag_nflds = 4
  integer, public, parameter :: mag2opls_nflds = 3
  integer, public, parameter :: mag2opls_ref_p_nflds = 3

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_esmf_fields_rhandles_init()
    use params_module, only: nz=>nhgt_fix
    use edyn3d_esmf_phys_mesh_mod, only: phys_mesh
    use edyn3d_esmf_oplus_grid_mod, only: oplus_grid
    use edyn3d_esmf_s1_mag_grid_mod, only: mag_s1_fdln_grid
    use edyn3d_esmf_s2_mag_grid_mod, only: mag_s2_fdln_grid
    use edyn3d_esmf_mag_ref_p_grid_mod, only: mag_ref_p_fdln_grid

    integer :: k, rc, astat
    type(ESMF_ArraySpec) :: arrayspec

    integer :: smm_srctermproc, smm_pipelinedep
    integer(ESMF_KIND_I4), pointer :: factorIndexList(:,:)
    real(ESMF_KIND_R8),    pointer :: factorList(:)

    character(len=*), parameter :: subname = 'edyn3d_esmf_fields_rhandles_init'

    ! Physics Field

    call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=rc)
    call check_error(subname,'ESMF_ArraySpecSet',rc)

    physFieldSrc = ESMF_FieldCreate(phys_mesh, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, &
                                    ungriddedLBound=(/1/), ungriddedUBound=(/phys2mag_nflds/), rc=rc)
    call check_error(subname,'ESMF_FieldCreate physFieldSrc',rc)

    ! Oplus field

    call ESMF_ArraySpecSet(arrayspec, 3, ESMF_TYPEKIND_R8, rc=rc)
    call check_error(subname,'ESMF_ArraySpecSet for oplusFeildDes',rc)

    oplusFieldDes = ESMF_FieldCreate(oplus_grid, arrayspec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         name="oplusFieldDes", ungriddedLBound=(/1/), ungriddedUBound=(/mag2opls_nflds/), rc=rc)
    call check_error(subname,'ESMF_FieldCreate oplusFieldDes',rc)

    oplusFieldDes_ref_p = ESMF_FieldCreate(oplus_grid, arrayspec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         name="oplusFieldDes", ungriddedLBound=(/1/), ungriddedUBound=(/mag2opls_ref_p_nflds/), rc=rc)
    call check_error(subname,'ESMF_FieldCreate oplusFieldDes',rc)

    ! nz vertical mag field-line levels

    allocate(magFieldDes_s1(nz), stat=astat)
    if (astat/=0) then
       call endrun(subname//' : not able to allocate magFieldDes_s1')
    end if
    allocate(magFieldDes_s2(nz), stat=astat)
    if (astat/=0) then
       call endrun(subname//' : not able to allocate magFieldDes_s2')
    end if
    allocate(magFieldSrc_s2(nz), stat=astat)
    if (astat/=0) then
       call endrun(subname//' : not able to allocate magFieldSrc_s2')
    end if

    allocate(rh_phys2mag_s1(nz), stat=astat)
    if (astat/=0) then
       call endrun(subname//' : not able to allocate rh_phys2mag_s1')
    end if
    allocate(rh_phys2mag_s2(nz), stat=astat)
    if (astat/=0) then
       call endrun(subname//' : not able to allocate rh_phys2mag_s2')
    end if
    allocate(rh_mag2oplus_s2(nz), stat=astat)
    if (astat/=0) then
       call endrun(subname//' : not able to allocate rh_mag2oplus_s2')
    end if

    ! Create route handles

    smm_srctermproc = 0
    smm_pipelinedep = 16

    vertloop: do k = 1, nz

       ! Mag fields

       magFieldDes_s2(k) = ESMF_FieldCreate( grid=mag_s2_fdln_grid(k), &
            staggerloc=ESMF_STAGGERLOC_CENTER, typekind=ESMF_TYPEKIND_R8, &
            ungriddedLBound=(/1/), ungriddedUBound=(/phys2mag_nflds/), rc=rc)
       call check_error(subname,'ESMF_FieldCreate magFieldDes_s2',rc)

       magFieldDes_s1(k) = ESMF_FieldCreate( grid=mag_s1_fdln_grid(k), &
            staggerloc=ESMF_STAGGERLOC_CENTER, typekind=ESMF_TYPEKIND_R8, &
            ungriddedLBound=(/1/), ungriddedUBound=(/phys2mag_nflds/), rc=rc)
       call check_error(subname,'ESMF_FieldCreate magFieldDes_s1',rc)

       magFieldSrc_s2(k) = ESMF_FieldCreate( grid=mag_s2_fdln_grid(k), &
            staggerloc=ESMF_STAGGERLOC_CENTER, typekind=ESMF_TYPEKIND_R8, &
            ungriddedLBound=(/1/), ungriddedUBound=(/mag2opls_nflds/), rc=rc)
       call check_error(subname,'ESMF_FieldCreate magFieldSrc_s2',rc)

       ! phys->mag s1
       call ESMF_FieldRegridStore( &
            srcField=physFieldSrc, dstField=magFieldDes_s1(k), &
            routehandle=rh_phys2mag_s1(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call check_error(subname,'FieldRegridStore rh_phys2mag_s1 route handle',rc)

       ! phys->mag s2
       call ESMF_FieldRegridStore( &
            srcField=physFieldSrc, dstField=magFieldDes_s2(k), &
            routehandle=rh_phys2mag_s2(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call check_error(subname,'FieldRegridStore rh_phys2mag_s2 route handle',rc)

       ! mag s2 -> oplus
       call ESMF_FieldRegridStore( &
            srcField=magFieldSrc_s2(k), dstField=oplusFieldDes, &
            routehandle=rh_mag2oplus_s2(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call check_error(subname,'FieldRegridStore rh_mag2oplus_s2(k) route handle',rc)

    end do vertloop

    magFieldSrc_ref_p = ESMF_FieldCreate( grid=mag_ref_p_fdln_grid, &
         staggerloc=ESMF_STAGGERLOC_CENTER, typekind=ESMF_TYPEKIND_R8, &
         ungriddedLBound=(/1/), ungriddedUBound=(/mag2opls_ref_p_nflds/), rc=rc)
    call check_error(subname,'ESMF_FieldCreate magFieldSrc_s2',rc)

    ! mag ref p -> oplus
    call ESMF_FieldRegridStore( &
         srcField=magFieldSrc_ref_p, dstField=oplusFieldDes_ref_p, &
         routehandle=rh_mag2oplus_ref_p, &
         regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
         polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
         extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
         factorIndexList=factorIndexList,                                   &
         factorList=factorList, srcTermProcessing=smm_srctermproc,          &
         pipelineDepth=smm_pipelinedep, rc=rc)
    call check_error(subname,'FieldRegridStore rh_mag2oplus_s2(k) route handle',rc)


  end subroutine edyn3d_esmf_fields_rhandles_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine edyn3d_esmf_fields_rhandles_destroy()
    use params_module, only: nz=>nhgt_fix

    integer :: k, rc
    character(len=*), parameter :: subname = 'edyn3d_esmf_fields_rhandles_destroy'

    call ESMF_FieldDestroy(physFieldSrc, rc=rc)
    call check_error(subname,'ESMF_FieldDestroy physFieldSrc', rc)

    call ESMF_FieldDestroy(oplusFieldDes, rc=rc)
    call check_error(subname,'ESMF_FieldDestroy oplusFieldDes', rc)

    call ESMF_FieldDestroy(magFieldSrc_ref_p, rc=rc)
    call check_error(subname,'ESMF_FieldDestroy magFieldSrc_ref_p', rc)

    call ESMF_FieldDestroy(oplusFieldDes_ref_p, rc=rc)
    call check_error(subname,'ESMF_FieldDestroy oplusFieldDes_ref_p', rc)

    do k = 1, nz
       call ESMF_FieldDestroy(magFieldDes_s2(k), rc=rc)
       call check_error(subname,'ESMF_FieldDestroy magFieldDes_s2', rc)

       call ESMF_FieldDestroy(magFieldDes_s1(k), rc=rc)
       call check_error(subname,'ESMF_FieldDestroy magFieldDes_s1', rc)

       call ESMF_FieldDestroy(magFieldSrc_s2(k), rc=rc)
       call check_error(subname,'ESMF_FieldDestroy magFieldSrc_s2', rc)

       call ESMF_RouteHandleDestroy(rh_phys2mag_s1(k), rc=rc)
       call check_error(subname,'ESMF_RouteHandleDestroy rh_phys2mag_s1', rc)

       call ESMF_RouteHandleDestroy(rh_phys2mag_s2(k), rc=rc)
       call check_error(subname,'ESMF_RouteHandleDestroy rh_phys2mag_s2', rc)

       call ESMF_RouteHandleDestroy(rh_mag2oplus_s2(k), rc=rc)
       call check_error(subname,'ESMF_RouteHandleDestroy rh_mag2oplus_s2', rc)
    end do

    call ESMF_RouteHandleDestroy(rh_mag2oplus_ref_p, rc=rc)
    call check_error(subname,'ESMF_RouteHandleDestroy rh_mag2oplus_ref_p', rc)

    deallocate(magFieldDes_s2)
    deallocate(magFieldDes_s1)
    deallocate(magFieldSrc_s2)
    deallocate(rh_phys2mag_s1)
    deallocate(rh_phys2mag_s2)
    deallocate(rh_mag2oplus_s2)

  end subroutine edyn3d_esmf_fields_rhandles_destroy

end module edyn3D_esmf_fields_rhandles
