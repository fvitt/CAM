module edyn3D_remap_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use spmd_utils, only: masterproc
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use edyn3D_mpi, only: ntask, mytid
  use edyn3d_params, only: hgt_fix,nhgt_fix
  use edyn3D_fline_fields, only: magfield_t
  use interpolate_data, only: lininterp

  use ESMF

  implicit none

  private

  public :: edyn3D_remap_phys2mag
  public :: edyn3D_remap_mag2oplus

contains

  subroutine edyn3D_remap_phys2mag(physflds, physalt, nphyscol, nphyslev, nflds, desfields, routehandles, magflds)

    use edyn3D_esmf_fields_rhandles, only: physFieldSrc

    integer,  intent(in) :: nphyscol,nphyslev, nflds
    real(r8), intent(in) :: physflds(nphyslev,nphyscol,nflds)
    real(r8), intent(in) :: physalt(nphyslev,nphyscol)

    character(len=*), parameter :: subname = 'edyn3D_remap_phys2mag'

    type(ESMF_Field), intent(inout) :: desfields(:)
    type(ESMF_RouteHandle), intent(inout) :: routehandles(:)
    type(magfield_t) , intent(inout) :: magflds(nflds)

    real(r8) :: physflds_tmp(nphyscol,nhgt_fix,nflds)

    integer :: n, i, j, k, isn, jj, nmlat, rc

    real(kind=ESMF_KIND_R8), pointer :: fptr2d(:,:)
    real(kind=ESMF_KIND_R8), pointer :: fptr3d(:,:,:)

    integer :: lbnd2d(2), ubnd2d(2) !
    integer :: lbnd3d(3), ubnd3d(3) !

    do n = 1,nflds
      do i = 1,nphyscol
         call lininterp(physflds(nphyslev:1:-1,i,n),physalt(nphyslev:1:-1,i),nphyslev,&
                        physflds_tmp(i,:,n),hgt_fix(:),nhgt_fix)
      end do
    end do

    do k = 1,nhgt_fix

       call ESMF_FieldGet(field=physFieldSrc, localDe=0, farrayPtr=fptr2d, &
                          computationalLBound=lbnd2d, computationalUBound=ubnd2d, rc=rc)
       call check_error(subname,'ESMF_FieldGet physFieldSrc',rc)

       do n = lbnd2d(2), ubnd2d(2) ! 1,nflds
          do i = lbnd2d(1), ubnd2d(1)
             fptr2d(i,n) = physflds_tmp(i,k,n)
          end do
       end do

       call ESMF_FieldRegrid(physFieldSrc, desfields(k), routehandles(k), &
                             termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
       call check_error(subname,'ESMF_FieldRegrid phys2mag',rc)

       if (mytid<ntask) then

          call ESMF_FieldGet(desfields(k), localDe=0, farrayPtr=fptr3d, &
                             computationalLBound=lbnd3d, computationalUBound=ubnd3d, rc=rc)
          call check_error(subname,'ESMF_FieldGet desfields',rc)

          nmlat = (magflds(1)%nmlat_h - (k-1))*2

          do n = lbnd3d(3), ubnd3d(3) ! 1,nflds
             do j = lbnd3d(2), ubnd3d(2)
                if (j>nmlat/2) then
                   isn = 2
                   jj = nmlat-j+1
                else
                   isn = 1
                   jj = j
                end if
                do i = lbnd3d(1), ubnd3d(1)
                   magflds(n)%flines(i,jj,isn)%fld(k) = fptr3d(i,j,n)
                end do
             end do
          end do

       end if

    end do

  end subroutine edyn3D_remap_phys2mag



  subroutine edyn3D_remap_mag2oplus( magflds, opalt, oplusflds )

    use edyn3D_esmf_fields_rhandles, only: magFieldSrc_s1, oplusFieldDes, rh_mag2plus_s1
    use edyn3D_esmf_fields_rhandles, only: nflds => mag2opls_nflds
    use edyn_mpi, only: lon0,lon1,lat0,lat1
    use edyn_geogrid, only: nlevo=>nlev

    character(len=*), parameter :: subname = 'edyn3D_remap_mag2oplus'

    type(magfield_t) , intent(inout) :: magflds(nflds)
    real(r8), intent(in) :: opalt(lon0:lon1,lat0:lat1,nlevo) ! oplus grid altitudes
    real(r8), intent(out) :: oplusflds(lon0:lon1,lat0:lat1,nlevo,nflds) ! field mapped to oplus grid

    real(r8) :: f_tmp(lon0:lon1,lat0:lat1,nflds,nhgt_fix)
    integer :: lbnd3d(3), ubnd3d(3) ! field bounds
    real(ESMF_KIND_R8), pointer :: fptr3d(:,:,:)

    integer :: i,j,k,jj,isn,n, nmlat, rc

    if (mytid<ntask) then

       do k = 1,nhgt_fix

          call ESMF_FieldGet(magFieldSrc_s1(k), localDe=0, farrayPtr=fptr3d, &
               computationalLBound=lbnd3d, computationalUBound=ubnd3d, rc=rc)
          call check_error(subname,'ESMF_FieldGet magFieldSrc_s1(k)',rc)

          nmlat = (magflds(1)%nmlat_h - (k-1))*2

          do n = lbnd3d(3), ubnd3d(3) ! 1,nflds
             do j = lbnd3d(2), ubnd3d(2)
                if (j>nmlat/2) then
                   isn = 2
                   jj = nmlat-j+1
                else
                   isn = 1
                   jj = j
                end if
                do i = lbnd3d(1), ubnd3d(1)
                   fptr3d(i,j,n) = magflds(n)%flines(i,jj,isn)%fld(k)
                end do
             end do
          end do

          call ESMF_FieldRegrid(magFieldSrc_s1(k), oplusFieldDes, rh_mag2plus_s1(k), &
               termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
          call check_error(subname,'ESMF_FieldRegrid mag2oplus',rc)

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

       end do

       do n = 1,nflds
          do i = lon0,lon1
             do j = lat0,lat1
                !vert interpolate...
                call lininterp(f_tmp(i,j,n,:), hgt_fix(:), nhgt_fix, &
                               oplusflds(i,j,:,n), opalt(i,j,:), nlevo )

             end do
          end do
       end do

    end if


  end subroutine edyn3D_remap_mag2oplus


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

end module edyn3D_remap_mod
