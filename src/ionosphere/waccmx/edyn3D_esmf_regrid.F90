module edyn3D_esmf_regrid
  use shr_kind_mod, only: r8 => shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc

  use ESMF

  use ppgrid,       only: begchunk, endchunk
  use phys_grid,    only: get_ncols_p, get_gcol_p
  use edyn3d_params, only: nz=>nhgt_fix, nmlat_h, nmlats2_h
  use edyn3D_mpi, only: ntask, tasks, mytid
  use edyn_mpi, only: edyn_mytid=>mytid, edyn_ntask=>ntask, edyn_ntaski=>ntaski, edyn_ntaskj=>ntaskj, edyn_tasks=>tasks
  use edyn_geogrid, only: glon, glat

  implicit none

  type(ESMF_Grid),  pointer :: mag_grid_p_src(:)
  type(ESMF_Grid),  pointer :: mag_grid_p_des(:)

  type(ESMF_Field), pointer :: magField_p_src(:)
  type(ESMF_Field), pointer :: magField_p_des(:)

  type(ESMF_Grid),  pointer :: mag_grid_s1(:) ! destination only
  type(ESMF_Field), pointer :: magField_s1(:)
  type(ESMF_Grid),  pointer :: mag_grid_s2(:)
  type(ESMF_Field), pointer :: magField_s2(:)

  type(ESMF_Grid)  :: oplusGrid
  type(ESMF_Field) :: oplusField

  type(ESMF_RouteHandle), pointer :: rh_phys2mag_p(:)
  type(ESMF_RouteHandle), pointer :: rh_phys2mag_s1(:), rh_phys2mag_s2(:)
  type(ESMF_RouteHandle), pointer :: rh_mag_p2phys(:), rh_mag_p2oplus(:)

  ! dist_grid_2d: DistGrid for 2D physics fields
  type(ESMF_DistGrid) :: dist_grid_2d
  type(ESMF_Mesh)  :: phys_mesh
  type(ESMF_Field) :: physField

  integer, allocatable :: petmap(:,:,:)

contains

  subroutine edyn3D_esmf_regrid_init
    use phys_control, only: phys_getopts
    use edyn3D_fieldline, only: fline_p, fline_s1, fline_s2
    use edyn3d_mpi, only: mlon0_p,mlon1_p

    character(len=cl) :: mesh_file

    integer :: total_cols, rc, i,j,k,isn,jj,nmlat, nmlat_s2
    integer :: ncols, chnk, col, dindex
    integer,allocatable :: decomp(:)

    type(ESMF_ArraySpec) :: arrayspec

    integer :: nmlons_task(ntask) ! number of lons per task
    integer :: nmlats_task(1) ! number of lats per task

    integer :: sp
    real(kind=ESMF_KIND_R8), pointer :: fptr2d(:,:)

    character(len=*), parameter :: subname = 'edyn3D_esmf_regrid_init'
    integer(ESMF_KIND_I4), pointer :: factorIndexList(:,:)
    real(ESMF_KIND_R8),    pointer :: factorList(:)
    integer                        :: smm_srctermproc,  smm_pipelinedep

    integer :: petcnt
    integer :: nlons_task(edyn_ntaski) ! # number of lons per task
    integer :: nlats_task(edyn_ntaskj) ! # number of lats per task
    real(ESMF_KIND_R8), pointer   :: coordX(:), coordY(:)
    integer :: lbnd(1), ubnd(1), n

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
    call check_error(subname,'ESMF_DistGridCreate phys decomp',rc)

    phys_mesh = ESMF_MeshCreate(trim(mesh_file), ESMF_FILEFORMAT_ESMFMESH,  &
                                elementDistgrid=dist_grid_2d, rc=rc)
    call check_error(subname,'ESMF_MeshCreate phys_mesh',rc)

    call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=rc)
    call check_error(subname,'ESMF_ArraySpecSet',rc)

    physField = ESMF_FieldCreate(phys_mesh, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    call check_error(subname,'ESMF_FieldCreate physField',rc)


    ! Oplus grid / field

    nlons_task = 0
    nlats_task = 0

    allocate(petmap(edyn_ntaski,edyn_ntaskj,1))

    petcnt = 0
    do j = 1,edyn_ntaskj
       do i = 1,edyn_ntaski
          petmap(i,j,1) = petcnt
          petcnt = petcnt+1
       end do
    end do

    do i = 1,edyn_ntaski
       loop: do n = 0, edyn_ntask-1
          if (edyn_tasks(n)%mytidi == i-1) then
             nlons_task(i) = edyn_tasks(n)%nlons
             exit loop
          end if
       end do loop
    end do
    !
    do j = 1, edyn_ntaskj
       loop1: do n = 0, edyn_ntask-1
          if (edyn_tasks(n)%mytidj == j-1) then
             nlats_task(j) = edyn_tasks(n)%nlats
             exit loop1
          end if
       end do loop1
    end do

    oplusGrid = ESMF_GridCreate1PeriDim(               &
           countsPerDEDim1=nlons_task, coordDep1=(/1/), &
           countsPerDEDim2=nlats_task, coordDep2=(/2/), &
           petmap=petmap, &
           coordSys=ESMF_COORDSYS_SPH_DEG, &
           indexflag=ESMF_INDEX_GLOBAL,minIndex=(/1,1/), rc=rc)
    call check_error(subname,'ESMF_GridCreate1PeriDim oplusGrid',rc)

    ! set up coordinates:

    call ESMF_GridAddCoord(oplusGrid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call check_error(subname,'ESMF_GridAddCoord  oplusGrid',rc)

    if (mytid<edyn_ntask) then

       ! Lon Coord
       call ESMF_GridGetCoord(oplusGrid, coordDim=1, localDE=0, &
            computationalLBound=lbnd, computationalUBound=ubnd, &
            farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
       call check_error(subname,'ESMF_GridGetCoord longitude coord oplusGrid',rc)

       do i = lbnd(1),ubnd(1)
          coordX(i) = glon(i)
       end do

       ! Lat Coord
       call ESMF_GridGetCoord(oplusGrid, coordDim=2, localDE=0, &
            computationalLBound=lbnd, computationalUBound=ubnd, &
            farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
       call check_error(subname,'ESMF_GridGetCoord latitude coord oplusGrid',rc)

       do j = lbnd(1),ubnd(1)
          coordY(j) = glat(j)
       end do

    end if

    ! Create 2D field
    call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=rc)
    call check_error(subname,'ESMF_ArraySpecSet for oplusFeild',rc)

    oplusField = ESMF_FieldCreate(oplusGrid, arrayspec, &
         staggerloc=ESMF_STAGGERLOC_CENTER, name="oplusField", rc=rc)
    call check_error(subname,'ESMF_FieldCreate oplusFeild',rc)

    deallocate(petmap)

    ! 2 for fline_s1 and fline_s2
    allocate(mag_grid_p_src(nz))
    allocate(mag_grid_p_des(nz))

    allocate(magField_p_src(nz))
    allocate(magField_p_des(nz))

    allocate(mag_grid_s1(nz))
    allocate(magField_s1(nz))

    allocate(mag_grid_s2(nz))
    allocate(magField_s2(nz))

    allocate(rh_phys2mag_p(nz))
    allocate(rh_phys2mag_s1(nz))
    allocate(rh_phys2mag_s2(nz))

    allocate(rh_mag_p2phys(nz))
    allocate(rh_mag_p2oplus(nz))

    allocate(petmap(ntask,1,1))

    do i = 1,ntask
       nmlons_task(i) = tasks(i-1)%nmaglons
       petmap(i,1,1) = i-1
    end do

    do k = 1,nz

       ! des mag p grid
       nmlat = (nmlat_h - (k-1))*2
       nmlats_task(1) = nmlat
       mag_grid_p_des(k) = ESMF_GridCreate1PeriDim( &
            countsPerDEDim1=nmlons_task, countsPerDEDim2=nmlats_task, &
            coordDep1=(/1, 2/), coordDep2=(/1, 2/), &
            petmap=petmap, indexflag=ESMF_INDEX_GLOBAL, rc=rc)
       call check_error(subname,'ESMF_GridCreate1PeriDim mag_grid_p_des',rc)

       call ESMF_GridAddCoord(grid=mag_grid_p_des(k),staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
       call check_error(subname,'ESMF_GridAddCoord mag_grid',rc)

       if (mytid<ntask) then
          !
          ! Get pointer and set mag grid geographic longitude coordinates:
          !
          call ESMF_GridGetCoord(grid=mag_grid_p_des(k), coordDim=1, farrayPtr=fptr2d, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord mag_grid_p_des glon',rc)

          do j = 1, nmlat
             if (j>nmlat/2) then
                isn = 2
                jj = nmlat-j+1
             else
                isn = 1
                jj = j
             endif
             do i = mlon0_p,mlon1_p
                fptr2d(i,j) = fline_p(i,jj,isn)%glon(k)
             end do
          enddo

          !
          ! Get pointer and set mag grid geographic latitdue coordinates:
          !
          call ESMF_GridGetCoord(grid=mag_grid_p_des(k), coordDim=2, farrayPtr=fptr2d, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord mag_grid_p_des glat',rc)

          do j = 1, nmlat
             if (j>nmlat/2) then
                isn = 2
                jj = nmlat-j+1
             else
                isn = 1
                jj = j
             endif
             do i = mlon0_p,mlon1_p
                fptr2d(i,j) = fline_p(i,jj,isn)%glat(k)
             end do
          enddo
       end if

       ! src mag p grid
       nmlat = (nmlat_h - (k-1))*2 - 1
       nmlats_task(1) = nmlat
       mag_grid_p_src(k) = ESMF_GridCreate1PeriDim( &
            countsPerDEDim1=nmlons_task, countsPerDEDim2=nmlats_task, &
            coordDep1=(/1, 2/), coordDep2=(/1, 2/), &
            petmap=petmap, indexflag=ESMF_INDEX_GLOBAL, rc=rc)
       call check_error(subname,'ESMF_GridCreate1PeriDim mag_grid_p_src',rc)

       call ESMF_GridAddCoord(grid=mag_grid_p_src(k),staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
       call check_error(subname,'ESMF_GridAddCoord mag_grid_p_src',rc)

       if (mytid<ntask) then
          !
          ! Get pointer and set mag grid geographic longitude coordinates:
          !
          call ESMF_GridGetCoord(grid=mag_grid_p_src(k), coordDim=1, farrayPtr=fptr2d, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord mag_grid_p_src glon',rc)

          do j = 1, nmlat
             if (j>nmlat/2) then
                isn = 2
                jj = nmlat-j+1
             else
                isn = 1
                jj = j
             endif
             do i = mlon0_p,mlon1_p
                fptr2d(i,j) = fline_p(i,jj,isn)%glon(k)
             end do
          enddo

          !
          ! Get pointer and set mag grid geographic latitdue coordinates:
          !
          call ESMF_GridGetCoord(grid=mag_grid_p_src(k), coordDim=2, farrayPtr=fptr2d, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord mag_grid_p_src glat',rc)

          do j = 1, nmlat
             if (j>nmlat/2) then
                isn = 2
                jj = nmlat-j+1
             else
                isn = 1
                jj = j
             endif
             do i = mlon0_p,mlon1_p
                fptr2d(i,j) = fline_p(i,jj,isn)%glat(k)
             end do
          enddo
       end if

       ! des and src mag p fields

       magField_p_des(k) = ESMF_FieldCreate( grid=mag_grid_p_des(k), typekind=ESMF_TYPEKIND_R8, rc=rc)
       call check_error(subname,'ESMF_FieldCreate magField_p_des',rc)

       magField_p_src(k) = ESMF_FieldCreate( grid=mag_grid_p_src(k), typekind=ESMF_TYPEKIND_R8, rc=rc)
       call check_error(subname,'ESMF_FieldCreate magField_p_src',rc)

       ! phys->mag
       call ESMF_FieldRegridStore( &
            srcField=physField, dstField=magField_p_des(k), &
            routehandle=rh_phys2mag_p(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call check_error(subname,'FieldRegridStore phys2mag route handle',rc)

       ! mag->phys
       call ESMF_FieldRegridStore( &
            srcField=magField_p_src(k), dstField=physField, &
            routehandle=rh_mag_p2phys(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call check_error(subname,'FieldRegridStore mag2phys route handle',rc)

       ! mag->oplus
       call ESMF_FieldRegridStore( &
            srcField=magField_p_src(k), dstField=oplusField, &
            routehandle=rh_mag_p2oplus(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call check_error(subname,'FieldRegridStore rh_mag2oplus route handle',rc)

       ! des mag s1 grid
       nmlat = (nmlat_h - (k-1))*2
       nmlats_task(1) = nmlat
       mag_grid_s1(k) = ESMF_GridCreate1PeriDim( &
            countsPerDEDim1=nmlons_task, countsPerDEDim2=nmlats_task, &
            coordDep1=(/1, 2/), coordDep2=(/1, 2/), &
            petmap=petmap, indexflag=ESMF_INDEX_GLOBAL, rc=rc)
       call check_error(subname,'ESMF_GridCreate1PeriDim mag_grid_s1',rc)

       call ESMF_GridAddCoord(grid=mag_grid_s1(k),staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
       call check_error(subname,'ESMF_GridAddCoord mag_grid_s1',rc)

       if (mytid<ntask) then
          !
          ! Get pointer and set mag grid geographic longitude coordinates:
          !
          call ESMF_GridGetCoord(grid=mag_grid_s1(k), coordDim=1, farrayPtr=fptr2d, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord glon',rc)

          do j = 1, nmlat
             if (j>nmlat/2) then
                isn = 2
                jj = nmlat-j+1
             else
                isn = 1
                jj = j
             endif
             do i = mlon0_p,mlon1_p
                fptr2d(i,j) = fline_p(i,jj,isn)%glon(k)
             end do
          enddo

          !
          ! Get pointer and set mag grid geographic latitdue coordinates:
          !
          call ESMF_GridGetCoord(grid=mag_grid_s1(k), coordDim=2, farrayPtr=fptr2d, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord glat',rc)

          do j = 1, nmlat
             if (j>nmlat/2) then
                isn = 2
                jj = nmlat-j+1
             else
                isn = 1
                jj = j
             endif
             do i = mlon0_p,mlon1_p
                fptr2d(i,j) = fline_s1(i,jj,isn)%glat(k)
             end do
          enddo
       end if

       magField_s1(k) = ESMF_FieldCreate( grid=mag_grid_s1(k), typekind=ESMF_TYPEKIND_R8, rc=rc)
       call check_error(subname,'ESMF_FieldCreate magField',rc)

       ! phys->mag
       call ESMF_FieldRegridStore( &
            srcField=physField, dstField=magField_s1(k), &
            routehandle=rh_phys2mag_s1(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call check_error(subname,'FieldRegridStore phys2mag route handle rh_phys2mag_s1',rc)

       ! s2 grid
       nmlat_s2 = (nmlatS2_h - (k-1))*2
       nmlats_task(1) = nmlat_s2
       mag_grid_s2(k) = ESMF_GridCreate1PeriDim( &
            countsPerDEDim1=nmlons_task, countsPerDEDim2=nmlats_task, &
            coordDep1=(/1, 2/), coordDep2=(/1, 2/), &
            petmap=petmap, indexflag=ESMF_INDEX_GLOBAL, rc=rc)
       call check_error(subname,'ESMF_GridCreate1PeriDim mag_grid_s2',rc)

       call ESMF_GridAddCoord(grid=mag_grid_s2(k),staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
       call check_error(subname,'ESMF_GridAddCoord mag_grid_s2',rc)

       if (mytid<ntask) then
          !
          ! Get pointer and set mag grid geographic longitude coordinates:
          !
          call ESMF_GridGetCoord(grid=mag_grid_s2(k), coordDim=1, farrayPtr=fptr2d, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord glon',rc)

          do j = 1, nmlat_s2
             if (j>nmlat_s2/2) then
                isn = 2
                jj = nmlat_s2-j+1
             else
                isn = 1
                jj = j
             endif
             do i = mlon0_p,mlon1_p
                fptr2d(i,j) = fline_s2(i,jj,isn)%glon(k)
             end do
          enddo

          !
          ! Get pointer and set mag grid geographic latitdue coordinates:
          !
          call ESMF_GridGetCoord(grid=mag_grid_s2(k), coordDim=2, farrayPtr=fptr2d, rc=rc)
          call check_error(subname,'ESMF_GridGetCoord glat',rc)

          do j = 1, nmlat_s2
             if (j>nmlat_s2/2) then
                isn = 2
                jj = nmlat_s2-j+1
             else
                isn = 1
                jj = j
             endif
             do i = mlon0_p,mlon1_p
                fptr2d(i,j) = fline_s2(i,jj,isn)%glat(k)
             end do
          enddo
       end if

       magField_s2(k) = ESMF_FieldCreate( grid=mag_grid_s2(k), typekind=ESMF_TYPEKIND_R8, rc=rc)
       call check_error(subname,'ESMF_FieldCreate magField',rc)

       ! phys->mag
       call ESMF_FieldRegridStore( &
            srcField=physField, dstField=magField_s2(k), &
            routehandle=rh_phys2mag_s2(k), &
            regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                           &
            polemethod=ESMF_POLEMETHOD_ALLAVG,                                 &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                      &
            factorIndexList=factorIndexList,                                   &
            factorList=factorList, srcTermProcessing=smm_srctermproc,          &
            pipelineDepth=smm_pipelinedep, rc=rc)
       call check_error(subname,'FieldRegridStore phys2mag route handle rh_phys2mag_s2',rc)

    end do

  end subroutine edyn3D_esmf_regrid_init

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine check_error(subname, routine, rc)
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
  end subroutine check_error

end module edyn3D_esmf_regrid
