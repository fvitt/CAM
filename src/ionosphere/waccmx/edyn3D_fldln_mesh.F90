module edyn3D_fldln_mesh
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  use edyn3d_params, only: nz=>nhgt_fix, nmlat_h, nmlon
  use edyn3d_mpi, only: mlon0_p,mlon1_p
  use edyn3D_fieldline, only: gmapex_p, fline_p
  use edyn3D_fieldline, only: glbgrid_p

  use edyn3D_fline_fields, only: magfield_t

  use ESMF

  implicit none

  type(ESMF_Mesh) :: fldln_mesh
  type(ESMF_Field) :: fldln_field

contains

  subroutine edyn3D_fldln_mesh_init

    integer, pointer :: nodeIds(:),nodeOwners(:)
    real(ESMF_KIND_R8), pointer :: nodeCoords(:)
    integer, pointer :: elemIds(:),elemTypes(:),elemConn(:)
    integer :: numNodes, numElems

    integer :: id, ndx, owner, ii, nElems, id0

    integer :: isn, i,j,k, nmlat
    integer :: nhexa, ntetr

    integer :: localrc, rc
    integer :: localPet, petCount

    type(ESMF_VM) :: vm

    type(ESMF_ArraySpec) :: arrayspec

    character(len=*), parameter :: subname = 'edyn3D_fldln_mesh_init'

    integer :: crdndx(mlon0_p:mlon1_p+1,nmlat_h,2,nz)

    numNodes = 0
    numElems = 0
    nElems = 0

    do k = 1,nz
       if (k==1) then
          nmlat = nmlat_h - 1
       else
          nmlat = nmlat_h - (k-1)
       end if

       do i = mlon0_p,mlon1_p+1

          do isn = 1,2

             numNodes = numNodes + nmlat
             if (k>1 .and. i<=mlon1_p) then
                numElems = numElems + (nmlat-1)
             end if
             if (k>1 .and. i==mlon0_p) then
                nElems = nElems + (nmlat-1)
             end if

          end do
       end do
    end do

    ! get pet info
    call ESMF_VMGetGlobal(vm, rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_VMGetGlobal vm error')
    end if


    call ESMF_VMGet(vm, petCount=petCount, localPet=localpet, rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_VMGet petCount, localPet error')
    end if

    ! Allocate and fill the node arrays.
    allocate(nodeIds(numNodes))
    allocate(nodeCoords(3*numNodes))
    allocate(nodeOwners(numNodes))

    ndx = 0

    do k = 1,nz
       if (k==1) then
          nmlat = nmlat_h - 1
       else
          nmlat = nmlat_h - (k-1)
       end if

       do i = mlon0_p,mlon1_p+1

          ii = i
          owner = localpet
          if (i==mlon1_p+1) then
             if (mlon1_p==nmlon) then
                ii = 1
                owner = 0
             else
                owner = owner+1
             end if
          end if

          do isn = 1,2

             do j = 1,nmlat
                ndx = ndx+1
                nodeIds(ndx) = glbgrid_p(ii,j,isn)%ID(k)
                nodeCoords(3*(ndx-1)+1) = glbgrid_p(ii,j,isn)%geolon(k)
                nodeCoords(3*(ndx-1)+2) = glbgrid_p(ii,j,isn)%geolat(k)
                nodeCoords(3*(ndx-1)+3) = glbgrid_p(ii,j,isn)%geoalt(k)
                nodeOwners(ndx) = owner
                crdndx(i,j,isn,k) = ndx
             end do
          end do
       end do
    end do

    ! Allocate and fill the element arrays.
    allocate(elemIds(numElems))
    allocate(elemTypes(numElems))
    allocate(elemConn(8*numElems))

    elemIds = -huge(1)

    ndx = 0
    id0 = nElems*(mlon0_p-1)

    do k = 2,nz
       nmlat = nmlat_h - (k-1) -1
       do i = mlon0_p,mlon1_p
          do isn = 1,2
             do j = 1,nmlat
                ndx = ndx+1
                elemIds(ndx) = id0+ndx
                elemTypes(ndx) = ESMF_MESHELEMTYPE_HEX
                elemConn(8*(ndx-1)+1) = crdndx(i,j,isn,k-1)
                elemConn(8*(ndx-1)+2) = crdndx(i,j+1,isn,k-1)
                elemConn(8*(ndx-1)+3) = crdndx(i+1,j,isn,k-1)
                elemConn(8*(ndx-1)+4) = crdndx(i+1,j+1,isn,k-1)
                elemConn(8*(ndx-1)+5) = crdndx(i,j,isn,k)
                elemConn(8*(ndx-1)+6) = crdndx(i,j+1,isn,k)
                elemConn(8*(ndx-1)+7) = crdndx(i+1,j,isn,k)
                elemConn(8*(ndx-1)+8) = crdndx(i+1,j+1,isn,k)
             end do
          end do
       end do
    end do

    ! Create 3D Mesh structure
    fldln_mesh = ESMF_MeshCreate(parametricDim=3,spatialDim=3, &
                                 nodeIds=nodeIds, nodeCoords=nodeCoords, &
                                 nodeOwners=nodeOwners, elementIds=elemIds,&
                                 elementTypes=elemTypes, elementConn=elemConn, &
                                 coordSys=ESMF_COORDSYS_SPH_DEG, rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_MeshCreate fldln_mesh')
    end if

    ! Create field
    call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_ArraySpecSet arrayspec ')
    end if
    fldln_field = ESMF_FieldCreate(fldln_mesh, arrayspec, name="MagFieldLine", rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_FieldCreate fldln_field ')
    end if


  end subroutine edyn3D_fldln_mesh_init

  subroutine edyn3D_fldln_mesh_setfld( magfld )

    type(magfield_t), intent(in) :: magfld


    integer :: localrc, ndx, nmlat, i,j,k, isn

    real(ESMF_KIND_R8), pointer :: farrayPtr1D(:)
    character(len=*), parameter :: subname = 'edyn3D_fldln_mesh_setfld'

    ! Load test data into the 3D mag field-line Field
    ! Should only be 1 localDE
    call ESMF_FieldGet(fldln_field, 0, farrayPtr1D,  rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_FieldGet fldln_field error')
    end if

    farrayPtr1D = -huge(1._r8)

    ndx = 0

    do k = 1,nz
       if (k==1) then
          nmlat = nmlat_h - 1
       else
          nmlat = nmlat_h - (k-1)
       end if

       do i = mlon0_p,mlon1_p
           do isn = 1,2
             do j = 1,nmlat
                ndx = ndx+1

                farrayPtr1D(ndx) = magfld%flines(i,j,isn)%fld(k)

             end do
          end do
       end do
    end do

  end subroutine edyn3D_fldln_mesh_setfld

  subroutine edyn3D_fldln_mesh_getfld( magfld )

    type(magfield_t), intent(inout) :: magfld

    integer :: localrc, ndx, nmlat, i,j,k, isn

    real(ESMF_KIND_R8), pointer :: farrayPtr1D(:)
    character(len=*), parameter :: subname = 'edyn3D_fldln_mesh_setfld'

    ! Load test data into the 3D mag field-line Field
    ! Should only be 1 localDE
    call ESMF_FieldGet(fldln_field, 0, farrayPtr1D,  rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_FieldGet fldln_field error')
    end if

    ndx = 0

    do k = 1,nz
       if (k==1) then
          nmlat = nmlat_h - 1
       else
          nmlat = nmlat_h - (k-1)
       end if

       do i = mlon0_p,mlon1_p
           do isn = 1,2
             do j = 1,nmlat
                ndx = ndx+1
                magfld%flines(i,j,isn)%fld(k) = farrayPtr1D(ndx)
             end do
          end do
       end do
    end do

    do i = mlon0_p,mlon1_p
       do isn = 1,2
          magfld%flines(i,nmlat_h,isn)%fld(1) = 0.5_r8 &
               *(magfld%flines(i,nmlat_h-1,1)%fld(1) + magfld%flines(i,nmlat_h-1,2)%fld(1))
       end do
    end do

  end subroutine edyn3D_fldln_mesh_getfld

end module edyn3D_fldln_mesh
