module edyn3d_remap_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use spmd_utils, only: masterproc, mpicom
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use mpi_module, only: ntask3D=>mpi_size, mytid=>mpi_rank
  use edyn_mpi,   only: ntaskOp=>ntask
  use params_module, only: hgt_fix, nhgt_fix, nmlat_h, nmlatS2_h
  use interpolate_data, only: lininterp

  use mpi_module, only: mlat0, mlat1, mlon0, mlon1

  use ESMF

  implicit none

  private

  public :: edyn3d_remap_phys2mag_s1
  public :: edyn3d_remap_phys2mag_s2
  public :: edyn3d_remap_mag2oplus
  public :: phys_fields_bundle_t
  public :: mag_fields_bundle_t
  public :: oplus_fields_bundle_t
  public :: NOTSET

  type :: phys_fields_bundle_t
     real(r8),  pointer :: fld(:,:)
  end type phys_fields_bundle_t

  type :: mag_fields_bundle_t
     real(r8),  pointer :: fld(:,:,:,:)
  end type mag_fields_bundle_t

  type :: oplus_fields_bundle_t
     real(r8),  pointer :: fld(:,:,:)
  end type oplus_fields_bundle_t

  real(r8), parameter :: NOTSET = -huge(1._r8)

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine edyn3d_remap_phys2mag_s1(nphyscol, nphyslev, physalt, physflds, magflds)

    use edyn3d_esmf_s1_mag_grid_mod, only: mag_s1_fdln_grid
    use fieldline_module, only: npts_s1
    use edyn3d_esmf_fields_rhandles, only: magFieldDes_s1, rh_phys2mag_s1, nflds=>phys2mag_nflds
    use edyn3D_esmf_fields_rhandles, only: physFieldSrc

    integer,  intent(in) :: nphyscol, nphyslev
    real(r8), intent(in) :: physalt(nphyslev,nphyscol)
    type(phys_fields_bundle_t), intent(in) :: physflds(nflds)
    type(mag_fields_bundle_t), intent(out) :: magflds(nflds)

    real(r8) :: physflds_tmp(nphyscol,nhgt_fix,nflds)

    integer :: n, i, j, k, isn, jj, nmlat, rc
    integer :: ncells_hlat, localDECount, nde

    real(kind=ESMF_KIND_R8), pointer :: fptr2d(:,:)
    real(kind=ESMF_KIND_R8), pointer :: fptr3d(:,:,:)

    integer :: lbnd2d(2), ubnd2d(2) !
    integer :: lbnd3d(3), ubnd3d(3) !

    character(len=*), parameter :: subname = 'edyn3d_remap_phys2mag'


    do n = 1,nflds
       magflds(n)%fld = NOTSET
       do i = 1,nphyscol
          call lininterp(physflds(n)%fld(nphyslev:1:-1,i),physalt(nphyslev:1:-1,i),nphyslev, &
                         physflds_tmp(i,:,n),hgt_fix(:),nhgt_fix)
       end do
    end do

    vertloop: do k = 1,nhgt_fix

       call ESMF_FieldGet(field=physFieldSrc, localDe=0, farrayPtr=fptr2d, &
                          computationalLBound=lbnd2d, computationalUBound=ubnd2d, rc=rc)
       call check_error(subname,'ESMF_FieldGet physFieldSrc',rc)
       fptr2d = NOTSET

       do n = lbnd2d(2), ubnd2d(2) ! 1,nflds
          do i = lbnd2d(1), ubnd2d(1)
             fptr2d(i,n) = physflds_tmp(i,k,n)
          end do
       end do

       call ESMF_FieldRegrid(physFieldSrc, magFieldDes_s1(k), rh_phys2mag_s1(k), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
       call check_error(subname,'ESMF_FieldRegrid phys2mag',rc)

       call ESMF_GridGet(mag_s1_fdln_grid(k), localDECount=localDECount, rc=rc)
       call check_error(subname,'ESMF_GridGet localDECount',rc)

       ! total number of grids cells per hemisphere
       ncells_hlat = nmlat_h - (k-1)

       DE_num: do nde = 0,localDECount-1

          call ESMF_FieldGet(magFieldDes_s1(k), localDe=nde, farrayPtr=fptr3d, &
               computationalLBound=lbnd3d, computationalUBound=ubnd3d, rc=rc)
          call check_error(subname,'ESMF_FieldGet magFieldDes_s1',rc)

          do n = lbnd3d(3), ubnd3d(3) ! 1,nflds
             do j = lbnd3d(2), ubnd3d(2)
                if (j>ncells_hlat) then
                   isn = 2
                   jj = 2*(ncells_hlat-1)+1 - j + 1
                else
                   isn = 1
                   jj = j
                end if
                do i = lbnd3d(1), ubnd3d(1)
                   magflds(n)%fld(k,isn,jj,i) = fptr3d(i,j,n)
                end do
                if (j==ncells_hlat) then ! at equator set point north to south
                   magflds(n)%fld(k,2,jj,:) = magflds(n)%fld(k,1,jj,:)
                end if
             end do
          end do

       end do DE_num

    end do vertloop

    do isn = 1,2
       do j = mlat0,mlat1
          do k = 1,npts_s1(j)
             do i = mlon0,mlon1
                do n = 1,nflds
                   if (magflds(n)%fld(k,isn,j,i)==NOTSET) then
                      write(*,*) subname,': magflds not set correctly at k,isn,j,i,n ',k,isn,j,i,n
                      call endrun(subname//': magflds not set correctly')
                   end if
                end do
             end do
          end do
       end do
    end do

  end subroutine edyn3d_remap_phys2mag_s1

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine edyn3d_remap_phys2mag_s2(nphyscol, nphyslev, physalt, physflds, magflds)

    use edyn3d_esmf_s2_mag_grid_mod, only: mag_s2_fdln_grid
    use fieldline_module, only: npts_s2
    use edyn3d_esmf_fields_rhandles, only: magFieldDes_s2, rh_phys2mag_s2, nflds=>phys2mag_nflds
    use edyn3D_esmf_fields_rhandles, only: physFieldSrc

    integer,  intent(in) :: nphyscol, nphyslev
    real(r8), intent(in) :: physalt(nphyslev,nphyscol)
    type(phys_fields_bundle_t), intent(in) :: physflds(nflds)
    type(mag_fields_bundle_t), intent(out) :: magflds(nflds)

    real(r8) :: physflds_tmp(nphyscol,nhgt_fix,nflds)

    integer :: n, i, j, k, isn, jj, nmlat, rc
    integer :: ncells_hlat, localDECount, nde

    real(kind=ESMF_KIND_R8), pointer :: fptr2d(:,:)
    real(kind=ESMF_KIND_R8), pointer :: fptr3d(:,:,:)

    integer :: lbnd2d(2), ubnd2d(2) !
    integer :: lbnd3d(3), ubnd3d(3) !

    character(len=*), parameter :: subname = 'edyn3d_remap_phys2mag'

    do n = 1,nflds
       magflds(n)%fld = NOTSET
       do i = 1,nphyscol
          call lininterp(physflds(n)%fld(nphyslev:1:-1,i),physalt(nphyslev:1:-1,i),nphyslev, &
                         physflds_tmp(i,:,n),hgt_fix(:),nhgt_fix)
       end do
    end do

    vertloop: do k = 1,nhgt_fix

       call ESMF_FieldGet(field=physFieldSrc, localDe=0, farrayPtr=fptr2d, &
                          computationalLBound=lbnd2d, computationalUBound=ubnd2d, rc=rc)
       call check_error(subname,'ESMF_FieldGet physFieldSrc',rc)
       fptr2d = NOTSET

       do n = lbnd2d(2), ubnd2d(2) ! 1,nflds
          do i = lbnd2d(1), ubnd2d(1)
             fptr2d(i,n) = physflds_tmp(i,k,n)
          end do
       end do

       call ESMF_FieldRegrid(physFieldSrc, magFieldDes_s2(k), rh_phys2mag_s2(k), &
                             termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
       call check_error(subname,'ESMF_FieldRegrid phys2mag',rc)

       call ESMF_GridGet(mag_s2_fdln_grid(k), localDECount=localDECount, rc=rc)
       call check_error(subname,'ESMF_GridGet localDECount',rc)

       ! total number of grids cells per hemisphere
       ncells_hlat = nmlatS2_h - (k-1)

       DE_num: do nde = 0,localDECount-1

          call ESMF_FieldGet(magFieldDes_s2(k), localDe=nde, farrayPtr=fptr3d, &
                             computationalLBound=lbnd3d, computationalUBound=ubnd3d, rc=rc)
          call check_error(subname,'ESMF_FieldGet magFieldDes_s2',rc)

          do n = lbnd3d(3), ubnd3d(3) ! 1,nflds
             do j = lbnd3d(2), ubnd3d(2)
                if (j>ncells_hlat) then
                   isn = 2
                   jj = 2*ncells_hlat - j + 1
                else
                   isn = 1
                   jj = j
                end if
                do i = lbnd3d(1), ubnd3d(1)
                   magflds(n)%fld(k,isn,jj,i) = fptr3d(i,j,n)
                end do
             end do
          end do

       end do DE_num

    end do vertloop

    do isn = 1,2
       do j = mlat0,min(mlat1,nmlatS2_h)
          do k = 1,npts_s2(j)
             do i = mlon0,mlon1
                do n = 1,nflds
                   if (magflds(n)%fld(k,isn,j,i)==NOTSET) then
                      write(*,*) subname,': magflds not set correctly at k,isn,j,i,n ',k,isn,j,i,n
                      call endrun(subname//': magflds not set correctly')
                   end if
                end do
             end do
          end do
       end do
    end do

  end subroutine edyn3d_remap_phys2mag_s2

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine edyn3d_remap_mag2oplus( magflds, opalt, oplusflds )

    use edyn3D_esmf_fields_rhandles, only: magFieldSrc_s2, oplusFieldDes, rh_mag2oplus_s2
    use edyn3D_esmf_fields_rhandles, only: nflds => mag2opls_nflds
    use edyn_mpi, only: lon0,lon1,lat0,lat1,lev0,lev1
    use edyn_geogrid, only: nlevo=>nlev
    use edyn3d_esmf_s2_mag_grid_mod, only: mag_s2_fdln_grid

    type(mag_fields_bundle_t), intent(in) :: magflds(nflds)
    real(r8), intent(in) :: opalt(lon0:lon1,lat0:lat1,lev0:lev1) ! oplus grid altitudes
    type(oplus_fields_bundle_t), intent(out) :: oplusflds(nflds) ! field mapped to oplus grid

    real(r8) :: f_tmp(lon0:lon1,lat0:lat1,nflds,nhgt_fix)
    integer :: lbnd3d(3), ubnd3d(3) ! field bounds
    real(ESMF_KIND_R8), pointer :: fptr3d(:,:,:)

    integer :: i,j,k,jj,isn,n, rc
    integer :: localDECount, nde, ncells_hlat

    character(len=*), parameter :: subname = 'edyn3d_remap_mag2oplus'

    do n = 1,nflds
       oplusflds(n)%fld = NOTSET
    end do

    vertloop: do k = 1,nhgt_fix

       call ESMF_GridGet(mag_s2_fdln_grid(k), localDECount=localDECount, rc=rc)
       call check_error(subname,'ESMF_GridGet localDECount',rc)

       ! total number of grids cells per hemisphere
       ncells_hlat = nmlatS2_h - (k-1)

       DE_num: do nde = 0,localDECount-1

          call ESMF_FieldGet(magFieldSrc_s2(k), localDe=nde, farrayPtr=fptr3d, &
               computationalLBound=lbnd3d, computationalUBound=ubnd3d, rc=rc)
          call check_error(subname,'ESMF_FieldGet magFieldSrc_s2(k)',rc)
          fptr3d = NOTSET

          do n = lbnd3d(3), ubnd3d(3) ! 1,nflds
             do j = lbnd3d(2), ubnd3d(2)
                if (j>ncells_hlat) then
                   isn = 2
                   jj = 2*ncells_hlat - j + 1
                else
                   isn = 1
                   jj = j
                end if
                do i = lbnd3d(1), ubnd3d(1)
                   fptr3d(i,j,n) = magflds(n)%fld(k,isn,jj,i)
                end do
             end do
          end do

          if (any(fptr3d==NOTSET)) then
             call endrun(subname//': fptr3d not set correctly')
          end if

       end do DE_num

       call ESMF_FieldRegrid(magFieldSrc_s2(k), oplusFieldDes, rh_mag2oplus_s2(k), &
            termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
       call check_error(subname,'ESMF_FieldRegrid mag2oplus',rc)

       if (mytid<ntaskOp) then

          call ESMF_FieldGet(field=oplusFieldDes, localDe=0, farrayPtr=fptr3d, &
               computationalLBound=lbnd3d, computationalUBound=ubnd3d, rc=rc)
          call check_error(subname,'ESMF_FieldGet oplusFieldDes',rc)

          do n = lbnd3d(3), ubnd3d(3) ! 1,nflds
             do j = lbnd3d(2), ubnd3d(2)
                do i = lbnd3d(1), ubnd3d(1)
                   f_tmp(i,j,n,k) = fptr3d(i,j,n)
                end do
             end do
          end do

       endif

    enddo vertloop

    if (any(f_tmp==NOTSET)) then
       call endrun(subname//': f_tmp not set correctly')
    end if

    do n = 1,nflds
       do i = lon0,lon1
          do j = lat0,lat1
             !vert interpolate...
             call lininterp(f_tmp(i,j,n,:), hgt_fix(:), nhgt_fix, oplusflds(n)%fld(i,j,:), opalt(i,j,:), nlevo )
          end do
       end do
       if (any(oplusflds(n)%fld==NOTSET)) then
          call endrun(subname//': oplusflds(n)%fld not set correctly')
       end if
    end do

  end subroutine edyn3d_remap_mag2oplus

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

end module edyn3d_remap_mod
