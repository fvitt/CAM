module esmf_phys2lonlat_mod
  use shr_kind_mod, only: r8 => shr_kind_r8, cl=>SHR_KIND_CL
  use cam_logfile,  only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils, only: masterproc
  use ppgrid, only: pver

  use ESMF

  implicit none

  private

  public :: esmf_phys2lonlat_init
  public :: esmf_phys2lonlat_regrid

  type(ESMF_RouteHandle) :: rh_phys2lonlat

  type(ESMF_Field) :: physfld
  type(ESMF_Field) :: lonlatfld

contains

  subroutine esmf_phys2lonlat_init()
    use esmf_phys_mesh_mod, only: physics_grid_mesh
    use esmf_lonlat_grid_mod, only: lonlat_grid

    type(ESMF_ArraySpec) :: arrayspec
    integer(ESMF_KIND_I4), pointer :: factorIndexList(:,:)
    real(ESMF_KIND_R8),    pointer :: factorList(:)
    integer                        :: smm_srctermproc,  smm_pipelinedep, rc

    character(len=*), parameter :: subname  = 'esmf_phys2lonlat_init: '

    ! create ESMF fields

    ! 3D phys fld
    call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_ArraySpecSet 3D phys fld ERROR')

    physfld = ESMF_FieldCreate(physics_grid_mesh, arrayspec, &
                               gridToFieldMap=(/2/), meshloc=ESMF_MESHLOC_ELEMENT, &
                               ungriddedLBound=(/1/), ungriddedUBound=(/pver/), rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldCreate 3D phys fld ERROR')

    ! 3D lon lat grid
    call ESMF_ArraySpecSet(arrayspec, 3, ESMF_TYPEKIND_R8, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_ArraySpecSet 3D lonlat fld ERROR')

    lonlatfld = ESMF_FieldCreate( lonlat_grid, arrayspec, staggerloc=ESMF_STAGGERLOC_CENTER, &
                                     ungriddedLBound=(/1/), ungriddedUBound=(/pver/), rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldCreate 3D lonlat fld ERROR')


    call ESMF_FieldRegridStore(srcField=physfld, dstField=lonlatfld, &
         regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
         polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
         extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
         routeHandle=rh_phys2lonlat, factorIndexList=factorIndexList, &
         factorList=factorList, srcTermProcessing=smm_srctermproc,          &
         pipelineDepth=smm_pipelinedep, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegridStore 3D routehandle ERROR')

  end subroutine esmf_phys2lonlat_init

 !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_phys2lonlat_regrid(physarr, lonlatarr)
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end
    use ppgrid, only: pcols, pver, begchunk, endchunk
    use phys_grid, only: get_ncols_p

    real(r8),intent(in) :: physarr(pver,pcols,begchunk:endchunk)
    real(r8),intent(out) :: lonlatarr(lon_beg:lon_end,lat_beg:lat_end,pver)

    integer :: i, ichnk, ncol, ilev, icol, rc
    real(ESMF_KIND_R8), pointer :: physptr(:,:)
    real(ESMF_KIND_R8), pointer :: lonlatptr(:,:,:)

    character(len=*), parameter :: subname = 'esmf_phys2lonlat_regrid: '

    call ESMF_FieldGet(physfld, localDe=0, farrayPtr=physptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet physptr')

    i = 0
    do ichnk = begchunk, endchunk
       ncol = get_ncols_p(ichnk)
       do icol = 1,ncol
          i = i+1
          do ilev = 1,pver
             physptr(ilev,i) = physarr(ilev,icol,ichnk)
          end do
       end do
    end do

    call ESMF_FieldRegrid(physfld, lonlatfld, rh_phys2lonlat, &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegrid physfld_3d->lonlatfld_3d')

    call ESMF_FieldGet(lonlatfld, localDe=0, farrayPtr=lonlatptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet lonlatptr')

    lonlatarr(lon_beg:lon_end,lat_beg:lat_end,1:pver) = lonlatptr(lon_beg:lon_end,lat_beg:lat_end,1:pver)

  end subroutine esmf_phys2lonlat_regrid

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine check_esmf_error( rc, errmsg )
    integer, intent(in) :: rc
    character(len=*), intent(in) :: errmsg

    character(len=cl) :: errstr

    if (rc /= ESMF_SUCCESS) then
       write(errstr,'(a,i6)') 'esmf_zonal_mod::'//trim(errmsg)//' -- ESMF ERROR code: ',rc
       if (masterproc) write(iulog,*) trim(errstr)
       call endrun(trim(errstr))
    end if

  end subroutine check_esmf_error

end module esmf_phys2lonlat_mod
