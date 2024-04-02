module edyn3D_regridder
  use shr_kind_mod, only: r8 => shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl
  use edyn3D_esmf_regrid, only: physField
  use edyn3D_mpi, only: mytid, ntask

  use ESMF

  implicit none

  contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine edyn3D_regridder_phys2mag(physfld,physalt,nphyscol,nphyslev,magfld)
    use edyn3d_params, only: hgt_fix,nhgt_fix
    use edyn3D_fline_fields, only: magfield_t
    use interpolate_data, only: lininterp

    integer,  intent(in) :: nphyscol,nphyslev
    real(r8), intent(in) :: physfld(nphyslev,nphyscol)
    real(r8), intent(in) :: physalt(nphyslev,nphyscol)
    type(magfield_t) , intent(inout) :: magfld

    real(r8) :: physfld_tmp(nphyscol,nhgt_fix)

    integer :: rc, i,j,k,isn,jj, nmlat
    character(len=*), parameter :: subname = 'edyn3D_regridder_phys2mag'

    integer :: lbnd1d(1), ubnd1d(1) ! 1d field bounds
    integer :: lbnd2d(2), ubnd2d(2) ! 2d field bounds

    real(ESMF_KIND_R8), pointer :: fptr2d(:,:)
    real(ESMF_KIND_R8), pointer :: fptr1d(:)

    do i = 1,nphyscol
       call lininterp(physfld(nphyslev:1:-1,i),physalt(nphyslev:1:-1,i),nphyslev,&
                      physfld_tmp(i,:),hgt_fix(:),nhgt_fix)
    end do

    do k = 1,nhgt_fix

       call ESMF_FieldGet(field=physField, localDe=0, farrayPtr=fptr1d, &
                          computationalLBound=lbnd1d, computationalUBound=ubnd1d, rc=rc)
       call check_errror(subname,'ESMF_FieldGet physField',rc)

       do i = lbnd1d(1), ubnd1d(1)
          fptr1d(i) = physfld_tmp(i,k)
       end do

       call ESMF_FieldRegrid(physField, magfld%esmf_fld(k), magfld%rhandle_phys2mag(k), &
                             termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
       call check_errror(subname,'ESMF_FieldRegrid phys2mag',rc)

       if (mytid<ntask) then

          call ESMF_FieldGet(magfld%esmf_fld(k), localDe=0, farrayPtr=fptr2d, &
                             computationalLBound=lbnd2d, computationalUBound=ubnd2d, rc=rc)
          call check_errror(subname,'ESMF_FieldGet physField',rc)

          nmlat = (magfld%nmlat_h - (k-1))*2

          do j = lbnd2d(2), ubnd2d(2)
             if (j<=magfld%nmlat_h-(k-1)) then
                isn = 1
                jj = j
             else
                isn = 2
                jj = nmlat-j+1
             end if
             do i = lbnd2d(1), ubnd2d(1)
                magfld%flines(i,jj,isn)%fld(k) = fptr2d(i,j)
             end do
          end do

       end if

    end do

  end subroutine edyn3D_regridder_phys2mag

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine edyn3D_regridder_mag2phys(magfld, physalt, nphyscol,nphyslev, physfld)
    use edyn3d_params, only: hgt_fix,nhgt_fix
    use edyn3D_fline_fields, only: magfield_t
    use interpolate_data, only: lininterp

    integer,  intent(in) :: nphyscol,nphyslev
    type(magfield_t) , intent(in) :: magfld
    real(r8), intent(in) :: physalt(nphyslev,nphyscol)
    real(r8), intent(inout) :: physfld(nphyslev,nphyscol)

    real(r8) :: physfld_tmp(nphyscol,nhgt_fix)

    integer :: rc, i,j,k,isn,jj, nmlat
    character(len=*), parameter :: subname = 'edyn3D_regridder_mag2phys'

    integer :: lbnd1d(1), ubnd1d(1) ! 1d field bounds
    integer :: lbnd2d(2), ubnd2d(2) ! 2d field bounds

    real(ESMF_KIND_R8), pointer :: fptr2d(:,:)
    real(ESMF_KIND_R8), pointer :: fptr1d(:)

    do k = 1,nhgt_fix

       if (mytid<ntask) then

          call ESMF_FieldGet(magfld%esmf_fld(k), localDe=0, farrayPtr=fptr2d, &
                             computationalLBound=lbnd2d, computationalUBound=ubnd2d, rc=rc)
          call check_errror(subname,'ESMF_FieldGet physField',rc)

          nmlat = (magfld%nmlat_h - (k-1))*2

          do j = lbnd2d(2), ubnd2d(2)
             if (j<=magfld%nmlat_h-(k-1)) then
                isn = 1
                jj = j
             else
                isn = 2
                jj = nmlat-j+1
             end if
             do i = lbnd2d(1), ubnd2d(1)
                fptr2d(i,j) = magfld%flines(i,jj,isn)%fld(k)
             end do
          end do

       end if

       call ESMF_FieldRegrid(magfld%esmf_fld(k), physField, magfld%rhandle_mag2phys(k), &
                             termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
       call check_errror(subname,'ESMF_FieldRegrid mag2phys',rc)

       call ESMF_FieldGet(field=physField, localDe=0, farrayPtr=fptr1d, &
                          computationalLBound=lbnd1d, computationalUBound=ubnd1d, rc=rc)
       call check_errror(subname,'ESMF_FieldGet physField',rc)

       do i = lbnd1d(1), ubnd1d(1)
          physfld_tmp(i,k) = fptr1d(i)
       end do

    end do

    do i = 1,nphyscol
       call lininterp(physfld_tmp(i,:),hgt_fix(:),nhgt_fix, &
                      physfld(nphyslev:1:-1,i),physalt(nphyslev:1:-1,i),nphyslev )
    end do

  end subroutine edyn3D_regridder_mag2phys

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine check_errror(subname, routine, rc)
    use shr_kind_mod, only: shr_kind_cl
    use spmd_utils, only: masterproc
    use cam_logfile, only: iulog
    use cam_abortutils, only: endrun

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
  end subroutine check_errror

end module edyn3D_regridder
