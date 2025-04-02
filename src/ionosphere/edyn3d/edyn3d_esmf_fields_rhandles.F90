module edyn3d_esmf_fields_rhandles
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc, mpicom

  use ESMF

  implicit none

  private

  public :: physFieldSrc
  public :: magFieldDes_s1
  public :: magFieldDes_s2
  public :: magFieldSrc_s2
  public :: oplusFieldDes

  public :: rh_mag2oplus_s2
  public :: rh_phys2mag_s1
  public :: rh_phys2mag_s2

  public :: edyn3d_esmf_fields_rhandles_init

  type(ESMF_Field) :: physFieldSrc
  type(ESMF_Field) :: oplusFieldDes

  type(ESMF_Field), allocatable :: magFieldDes_s1(:)
  type(ESMF_Field), allocatable :: magFieldDes_s2(:)
  type(ESMF_Field), allocatable :: magFieldSrc_s2(:)

  type(ESMF_RouteHandle) :: rh_phys2oplus
  type(ESMF_RouteHandle) :: rh_oplus2phys

  type(ESMF_RouteHandle), allocatable :: rh_phys2mag_s1(:)
  type(ESMF_RouteHandle), allocatable :: rh_phys2mag_s2(:)

  type(ESMF_RouteHandle), allocatable :: rh_mag2oplus_s2(:)

  integer, public, parameter :: phys2mag_nflds = 4
  integer, public, parameter :: mag2opls_nflds = 3

contains

  subroutine edyn3d_esmf_fields_rhandles_init
    use params_module, only: nz=>nhgt_fix
    use edyn3d_esmf_phys_mesh_mod, only: phys_mesh
    use edyn3d_esmf_oplus_grid_mod, only: oplus_grid
    use edyn3d_esmf_s1_mag_grid_mod, only: mag_s1_fdln_grid
    use edyn3d_esmf_s2_mag_grid_mod, only: mag_s2_fdln_grid

    integer :: k, rc
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

    ! nz vertical mag field-line levels

    allocate(magFieldDes_s1(nz))
    allocate(magFieldDes_s2(nz))
    allocate(magFieldSrc_s2(nz))

    allocate(rh_phys2mag_s1(nz))
    allocate(rh_phys2mag_s2(nz))
    allocate(rh_mag2oplus_s2(nz))

    ! Create route handles

    smm_srctermproc = 0
    smm_pipelinedep = 16

    vertloop: do k = 1, nz

       print*,' FVDBG... create rhandles lev k : ',k
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

!!$
!!$print*,'FVDBG.edyn3d_esmf_fields_rhandles_init.. OK HERE'
!!$call mpi_barrier(mpicom, rc)
!!$call endrun('FVDBG.edyn3d_esmf_fields_rhandles_init.. OK STOP HERE')
  end subroutine edyn3d_esmf_fields_rhandles_init

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine check_error(subname, routine, rc)

    character(len=*), intent(in) :: subname
    character(len=*), intent(in) :: routine
    integer,          intent(in) :: rc

    character(len=cl) :: errmsg

    if (rc /= ESMF_SUCCESS) then
       write(errmsg, '(4a,i0)') trim(subname), ': Error return from ', trim(routine), ', rc = ', rc
       if (masterproc) then
          write(iulog, '(2a)') 'ERROR: ', trim(errmsg)
       end if
       call endrun(trim(errmsg))
    end if
  end subroutine check_error

end module edyn3D_esmf_fields_rhandles
