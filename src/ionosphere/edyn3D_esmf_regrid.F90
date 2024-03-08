module edyn3D_esmf_regrid
  use shr_kind_mod, only: r8 => shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc

  use ESMF

  use ppgrid,       only: begchunk, endchunk
  use phys_grid,    only: get_ncols_p, get_gcol_p
  use edyn3d_params, only: nz=>nhgt_fix, nmlat_h

  implicit none

  type(ESMF_Field) :: physField
  type(ESMF_Field), allocatable :: magField(:)
  type(ESMF_RouteHandle), allocatable :: rh_phys2mag(:), rh_mag2phys(:)

  ! dist_grid_2d: DistGrid for 2D fields
  type(ESMF_DistGrid) :: dist_grid_2d

contains

  subroutine edyn3D_esmf_regrid_init
    use phys_control, only: phys_getopts
    use edyn3D_mpi, only: ntask, tasks
    use edyn3D_fieldline, only: fline_p
    use edyn3d_mpi, only: mlon0_p,mlon1_p

    character(len=cl) :: mesh_file

    integer :: total_cols, rc, i,j,k,isn,jj,ii,  nmlat
    integer :: ncols, chnk, col, dindex
    integer,allocatable :: decomp(:)

    type(ESMF_Grid) :: mag_grid
    type(ESMF_Mesh) :: phys_mesh
    type(ESMF_ArraySpec) :: arrayspec

    integer :: nmlons_task(ntask) ! number of lons per task
    integer :: nmlats_task(1) ! number of lats per task

    integer :: sp
    real(kind=ESMF_KIND_R8), pointer :: fptr2d(:,:)

    character(len=*), parameter :: subname = 'edyn3D_esmf_regrid_init'
    integer(ESMF_KIND_I4), pointer :: factorIndexList(:,:)
    real(ESMF_KIND_R8),    pointer :: factorList(:)
    integer                        :: smm_srctermproc,  smm_pipelinedep

    smm_srctermproc = 0
    smm_pipelinedep = 16

    call phys_getopts(physics_grid_out=mesh_file)

    ! Compute the local decomp
    total_cols = 0
    do chnk = begchunk, endchunk
       total_cols = total_cols + get_ncols_p(chnk)
    end do
    allocate(decomp(total_cols))
     dindex = 0
    do chnk = begchunk, endchunk
       ncols = get_ncols_p(chnk)
       do col = 1, ncols
          dindex = dindex + 1
          decomp(dindex) = get_gcol_p(chnk, col)
       end do
    end do

    ! Create a DistGrid based on the physics decomp
    dist_grid_2d = ESMF_DistGridCreate(arbSeqIndexList=decomp, rc=rc)
    call edyn3D_esmf_chkerr(subname,'ESMF_DistGridCreate phys decomp',rc)

    phys_mesh = ESMF_MeshCreate(trim(mesh_file), ESMF_FILEFORMAT_ESMFMESH,  &
                                elementDistgrid=dist_grid_2d, rc=rc)
    call edyn3D_esmf_chkerr(subname,'ESMF_MeshCreate phys_mesh',rc)

    call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=rc)
    call edyn3D_esmf_chkerr(subname,'ESMF_ArraySpecSet',rc)

    physField = ESMF_FieldCreate(phys_mesh, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    call edyn3D_esmf_chkerr(subname,'ESMF_FieldCreate physField',rc)

    ! 2 for fline_s1 and fline_s2
    allocate(magField(nz))
    allocate(rh_phys2mag(nz))
    allocate(rh_mag2phys(nz))

    do i = 1,ntask
       nmlons_task = tasks(i-1)%nmaglons
    end do

    do k = 1,nz

       nmlat = (nmlat_h - (k-1))*2
       nmlats_task(1) = nmlat
       mag_grid = ESMF_GridCreate1PeriDim( &
            countsPerDEDim1=nmlons_task, countsPerDEDim2=nmlats_task, &
            coordDep1=(/1, 2/), coordDep2=(/1, 2/), rc=rc)

       call ESMF_GridAddCoord(grid=mag_grid, rc=rc)
       call edyn3D_esmf_chkerr(subname,'ESMF_GridAddCoord mag_grid',rc)

       !
       ! Get pointer and set mag grid geographic longitude coordinates:
       !
       call ESMF_GridGetCoord(grid=mag_grid, coordDim=1, farrayPtr=fptr2d, rc=rc)
       call edyn3D_esmf_chkerr(subname,'ESMF_GridGetCoord glon',rc)

       do j = 1, nmlat
          if (j<=nmlat_h-(k-1)) then
             isn = 1
             jj = j
          else
             isn = 2
             jj = nmlat-j+1
          end if
          do i = mlon0_p,mlon1_p
             ii = i-mlon0_p+1
             fptr2d(ii,j) = fline_p(i,jj,isn)%glon(k)
          end do
       enddo

       !
       ! Get pointer and set mag grid geographic latitdue coordinates:
       !
       call ESMF_GridGetCoord(grid=mag_grid, coordDim=2, farrayPtr=fptr2d, rc=rc)
       call edyn3D_esmf_chkerr(subname,'ESMF_GridGetCoord glat',rc)

       do j = 1, nmlat
          if (j<=nmlat_h-(k-1)) then
             isn = 1
             jj = j
          else
             isn = 2
             jj = nmlat-j+1
          end if
          do i = mlon0_p,mlon1_p
             ii = i-mlon0_p+1
             fptr2d(ii,j) = fline_p(i,jj,isn)%glat(k)
          end do
       enddo

       magField(k) = ESMF_FieldCreate( grid=mag_grid, typekind=ESMF_TYPEKIND_R8, rc=rc)
       call edyn3D_esmf_chkerr(subname,'ESMF_FieldCreate magField',rc)

       ! phys->mag
       call ESMF_FieldRegridStore( &
            srcField=physField, dstField=magField(k), &
            routehandle=rh_phys2mag(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call edyn3D_esmf_chkerr(subname,'FieldRegridStore phys2mag route handle',rc)

       ! mag->phys
       call ESMF_FieldRegridStore( &
            srcField=magField(k), dstField=physField, &
            routehandle=rh_mag2phys(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call edyn3D_esmf_chkerr(subname,'FieldRegridStore mag2phys route handle',rc)

    end do

  end subroutine edyn3D_esmf_regrid_init


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine edyn3D_esmf_regrid_phys2mag(physfld,physalt,nphyscol,nphyslev,magfld)
    use edyn3d_params, only: hgt_fix,nhgt_fix
    use edyn3D_fieldline, only: magfld_t
    use edyn3d_mpi, only: mlon0_p,mlon1_p
    use interpolate_data, only: lininterp

    real(r8), intent(in) :: physfld(nphyscol,nphyslev)
    real(r8), intent(in) :: physalt(nphyslev)
    integer,  intent(in) :: nphyscol,nphyslev
    type(magfld_t), intent(out) :: magfld(mlon0_p:mlon1_p,nmlat_h,2)

    real(r8) :: physfld_tmp(nphyscol,nhgt_fix)

    integer :: rc, i,j,k,isn,jj, nmlat
    character(len=*), parameter :: subname = 'edyn3D_esmf_regrid_phys2mag'

    integer :: lbnd1d(1), ubnd1d(1) ! 1d field bounds
    integer :: lbnd2d(2), ubnd2d(2) ! 2d field bounds

    real(ESMF_KIND_R8), pointer :: fptr2d(:,:)
    real(ESMF_KIND_R8), pointer :: fptr1d(:)

    do i = 1,nphyscol
       call lininterp(physfld(i,nphyslev:1:-1),physalt(nphyslev:1:-1),nphyslev,&
                      physfld_tmp(i,:),hgt_fix(:),nhgt_fix)
    end do

    do k = 1,nhgt_fix

       call ESMF_FieldGet(field=physField, localDe=0, farrayPtr=fptr1d, &
            computationalLBound=lbnd1d, computationalUBound=ubnd1d, rc=rc)
       call edyn3D_esmf_chkerr(subname,'ESMF_FieldGet physField',rc)

       do i = lbnd1d(1), ubnd1d(1)
          fptr1d(i) = physfld_tmp(i,k)
       end do

       call ESMF_FieldRegrid(physField, magField(k), rh_phys2mag(k), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
       call edyn3D_esmf_chkerr(subname,'ESMF_FieldRegrid phys2mag',rc)

       call ESMF_FieldGet(magField(k), localDe=0, farrayPtr=fptr2d,             &
            computationalLBound=lbnd2d, computationalUBound=ubnd2d, rc=rc)
       call edyn3D_esmf_chkerr(subname,'ESMF_FieldGet physField',rc)

       nmlat = (nmlat_h - (k-1))*2

       do j = lbnd2d(2), ubnd2d(2)
          if (j<=nmlat_h-(k-1)) then
             isn = 1
             jj = j
          else
             isn = 2
             jj = nmlat-j+1
          end if
          do i = lbnd2d(1), ubnd2d(1)
             magfld(i,jj,isn)%fld(k) = fptr2d(i,j)
          end do
       end do

    end do

  end subroutine edyn3D_esmf_regrid_phys2mag


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine edyn3D_esmf_regrid_mag2phys
  end subroutine edyn3D_esmf_regrid_mag2phys

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine edyn3D_esmf_chkerr(subname, routine, rc)
    use shr_kind_mod, only: shr_kind_cl

    character(len=*), intent(in) :: subname
    character(len=*), intent(in) :: routine
    integer,          intent(in) :: rc

    character(len=shr_kind_cl) :: errmsg

    if (rc /= ESMF_SUCCESS) then
       write(errmsg, '(4a,i0)') trim(subname), ': Error return from ', trim(routine), ', rc = ', rc
       if (masterproc) then
          write(iulog, '(2a)') 'ERROR: ', trim(errmsg)
       end if
       call endrun(trim(errmsg))
    end if
  end subroutine edyn3D_esmf_chkerr

end module edyn3D_esmf_regrid
