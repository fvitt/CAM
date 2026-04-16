!------------------------------------------------------------------------------
! Provides methods for mapping between regular longitude / latitude grid and
! physics grid via ESMF regridding capabilities
!------------------------------------------------------------------------------
module esmf_lonlatphys_regrid_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile,  only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils, only: masterproc
  use ppgrid, only: pver

  use esmf_lonlat_grid_mod, only: mytid, lonlat_npes

  use ESMF, only: ESMF_RouteHandle, ESMF_Field, ESMF_ArraySpec, ESMF_ArraySpecSet
  use ESMF, only: ESMF_FieldCreate, ESMF_FieldRegridStore
  use ESMF, only: ESMF_FieldGet, ESMF_FieldRegrid
  use ESMF, only: ESMF_KIND_I4, ESMF_KIND_R8, ESMF_TYPEKIND_R8
  use ESMF, only: ESMF_REGRIDMETHOD_CONSERVE
  use ESMF, only: ESMF_TERMORDER_SRCSEQ, ESMF_MESHLOC_ELEMENT, ESMF_STAGGERLOC_CENTER
  use ESMF, only: ESMF_FieldDestroy, ESMF_RouteHandleDestroy
  use esmf_check_error_mod, only: check_esmf_error

  implicit none

  private

  public :: esmf_lonlatphys_regrid_init
  public :: esmf_lonlatphys_regrid_destroy
  public :: regrid_phys2lonlat
  public :: regrid_lonlat2phys

  type(ESMF_RouteHandle) :: rh_phys2lonlat
  type(ESMF_RouteHandle) :: rh_lonlat2phys

  type(ESMF_Field) :: physfld
  type(ESMF_Field) :: lonlatfld

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_lonlatphys_regrid_init()
    use esmf_phys_mesh_mod, only: physics_grid_mesh
    use esmf_lonlat_grid_mod, only: lonlat_grid

    type(ESMF_ArraySpec) :: arrayspec
    integer                        :: smm_srctermproc,  smm_pipelinedep, rc

    character(len=*), parameter :: subname  = 'esmf_lonlatphys_regrid_init: '

    smm_srctermproc = 0
    smm_pipelinedep = 16

    ! create ESMF fields 2-D ...

    ! phys fld
    call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_ArraySpecSet 2D phys fld ERROR')

    physfld = ESMF_FieldCreate(physics_grid_mesh, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldCreate 2D phys fld ERROR')

    ! lon lat fld
    call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_ArraySpecSet 2D lonlat fld ERROR')

    lonlatfld = ESMF_FieldCreate( lonlat_grid, arrayspec, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldCreate 2D lonlat fld ERROR')

    call ESMF_FieldRegridStore(srcField=physfld, dstField=lonlatfld, &
         regridMethod=ESMF_REGRIDMETHOD_CONSERVE, &
         routeHandle=rh_phys2lonlat, &
         srcTermProcessing=smm_srctermproc, &
         pipelineDepth=smm_pipelinedep, &
         rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegridStore phys2lonlat routehandle ERROR')

    call ESMF_FieldRegridStore(srcField=lonlatfld, dstField=physfld, &
         regridMethod=ESMF_REGRIDMETHOD_CONSERVE, &
         routeHandle=rh_lonlat2phys, &
         srcTermProcessing=smm_srctermproc, &
         pipelineDepth=smm_pipelinedep, &
         rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegridStore lonlat2phys routehandle ERROR')

  end subroutine esmf_lonlatphys_regrid_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine regrid_phys2lonlat(physarr, lonlatarr)
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end
    use ppgrid, only: pcols, pver, begchunk, endchunk
    use phys_grid, only: get_ncols_p

    real(r8), intent(in) :: physarr(pcols,begchunk:endchunk)
    real(r8), intent(out) :: lonlatarr(lon_beg:lon_end,lat_beg:lat_end)

    integer :: i, ichnk, ncol, icol, rc
    real(ESMF_KIND_R8), pointer :: physptr(:)
    real(ESMF_KIND_R8), pointer :: lonlatptr(:,:)

    character(len=*), parameter :: subname = 'regrid_phys2lonlat: '

    call ESMF_FieldGet(physfld, localDe=0, farrayPtr=physptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet physptr')

    i = 0
    do ichnk = begchunk, endchunk
      ncol = get_ncols_p(ichnk)
      do icol = 1,ncol
        i = i+1
        physptr(i) = physarr(icol,ichnk)
      end do
    end do

    call ESMF_FieldRegrid(physfld, lonlatfld, rh_phys2lonlat, &
                          termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegrid physfld->lonlatfld')

    if (mytid<lonlat_npes) then
       call ESMF_FieldGet(lonlatfld, localDe=0, farrayPtr=lonlatptr, rc=rc)
       call check_esmf_error(rc, subname//'ESMF_FieldGet lonlatptr')
       lonlatarr(lon_beg:lon_end,lat_beg:lat_end) = lonlatptr(lon_beg:lon_end,lat_beg:lat_end)
    endif

  end subroutine regrid_phys2lonlat


  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine regrid_lonlat2phys(lonlatarr, physarr)
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end
    use ppgrid, only: pcols, pver, begchunk, endchunk
    use phys_grid, only: get_ncols_p

    real(r8), intent(in) :: lonlatarr(lon_beg:lon_end,lat_beg:lat_end)
    real(r8), intent(out) :: physarr(pcols,begchunk:endchunk)

    integer :: i, ichnk, ncol, icol, rc
    real(ESMF_KIND_R8), pointer :: physptr(:)
    real(ESMF_KIND_R8), pointer :: lonlatptr(:,:)

    character(len=*), parameter :: subname = 'regrid_lonlat2phys: '

    if (mytid<lonlat_npes) then
       call ESMF_FieldGet(lonlatfld, localDe=0, farrayPtr=lonlatptr, rc=rc)
       call check_esmf_error(rc, subname//'ESMF_FieldGet lonlatptr')
       lonlatptr(lon_beg:lon_end,lat_beg:lat_end) = lonlatarr(lon_beg:lon_end,lat_beg:lat_end)
    end if

    call ESMF_FieldRegrid(lonlatfld, physfld, rh_lonlat2phys, &
                          termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegrid lonlatfld->physfld')

    call ESMF_FieldGet(physfld, localDe=0, farrayPtr=physptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet physptr')

    i = 0
    do ichnk = begchunk, endchunk
      ncol = get_ncols_p(ichnk)
      do icol = 1,ncol
        i = i+1
        physarr(icol,ichnk) = physptr(i)
      end do
    end do

  end subroutine regrid_lonlat2phys

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_lonlatphys_regrid_destroy()

    integer :: rc
    character(len=*), parameter :: subname = 'esmf_lonlatphys_regrid_destroy: '

    call ESMF_RouteHandleDestroy(rh_phys2lonlat, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy rh_phys2lonlat')

    call ESMF_RouteHandleDestroy(rh_lonlat2phys, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy rh_lonlat2phys')

    call ESMF_FieldDestroy(lonlatfld, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy lonlatfld')

    call ESMF_FieldDestroy(physfld, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy physfld')

  end subroutine esmf_lonlatphys_regrid_destroy

end module esmf_lonlatphys_regrid_mod
