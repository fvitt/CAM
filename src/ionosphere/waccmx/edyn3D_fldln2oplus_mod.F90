module edyn3D_fldln2oplus_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use infnan,         only: nan, assignment(=)
  use edyn3D_mpi, only: mytid, ntask

  use edyn3D_fline_fields, only: magfield_t
  use edyn3D_oplus_grid, only: oplus_field
  use edyn3D_fldln_mesh, only: fldln_field
  use edyn3D_fldln_mesh, only: edyn3D_fldln_mesh_setfld
  use edyn3D_oplus_grid, only: edyn3D_oplus_grid_getfld

  use ESMF

  implicit none

  type(ESMF_RouteHandle) :: routeHandle3D

contains

  subroutine edyn3D_fldln2oplus_init

    integer :: localrc, smm_srctermproc,  smm_pipelinedep
    character(len=*), parameter :: subname = 'edyn3D_fldln2oplus_init'

    smm_srctermproc = 0
    smm_pipelinedep = 16

    ! Regrid store
    call ESMF_FieldRegridStore( &
         fldln_field, &
         dstField=oplus_field, &
         routeHandle=routeHandle3D, &
         regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
         polemethod=ESMF_POLEMETHOD_ALLAVG, &
         extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG, &
         srcTermProcessing=smm_srctermproc, &
         pipelineDepth=smm_pipelinedep, rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_FieldRegridStore')
    end if

  end subroutine edyn3D_fldln2oplus_init

  subroutine edyn3D_fldln2oplus_flg2opg( opalt, magfld, opfld )
    use edyn_mpi, only: lon0,lon1,lat0,lat1,lev0,lev1
    use edyn_geogrid, only: nlevo=>nlev
    use edyn3D_oplus_grid, only: nlevg=>nlev, galt
    use interpolate_data, only: lininterp

    real(r8), intent(in) :: opalt(lon0:lon1,lat0:lat1,nlevo)
    type(magfield_t), intent(in) :: magfld
    real(r8), intent(out) :: opfld(lon0:lon1,lat0:lat1,nlevo)

    integer :: i,j
    integer :: localrc
    character(len=*), parameter :: subname = 'edyn3D_fldln2oplus_flg2opg'

    real(r8) :: fdata(lon0:lon1,lat0:lat1,lev0:lev1)

    call edyn3D_fldln_mesh_setfld( magfld )

    call ESMF_FieldRegrid(fldln_field, oplus_field, routeHandle3D, rc=localrc)
    if (ESMF_LogFoundError(localrc)) then
       call endrun(subname//' ESMF_FieldRegrid')
    end if

    call edyn3D_oplus_grid_getfld(fdata)

    do i = lon0,lon1
       do j = lat0,lat1
          !vert interpolate...
          call lininterp(fdata(i,j,:), galt(:), nlevg, &
                         opfld(i,j,:), opalt(i,j,:), nlevo )

       end do
    end do

  end subroutine edyn3D_fldln2oplus_flg2opg

  subroutine edyn3D_fldln2oplus_flg2opg_v2( opalt, magfld, opfld )
    use edyn3D_esmf_regrid, only: rh_mag2opls
    use edyn3D_oplus_grid, only: oplsField2D
    use edyn_mpi, only: lon0,lon1,lat0,lat1,lev0,lev1
    use edyn_geogrid, only: nlevo=>nlev
    use edyn3d_params, only: hgt_fix,nhgt_fix
    use interpolate_data, only: lininterp

    real(r8), intent(in) :: opalt(lon0:lon1,lat0:lat1,nlevo)
    type(magfield_t), intent(in) :: magfld
    real(r8), intent(out) :: opfld(lon0:lon1,lat0:lat1,nlevo)

    real(r8) :: f_tmp(lon0:lon1,lat0:lat1,nhgt_fix)
    integer :: lbnd1d(1), ubnd1d(1) ! 1d field bounds
    integer :: lbnd2d(2), ubnd2d(2) ! 2d field bounds
    real(ESMF_KIND_R8), pointer :: fptr2d(:,:)

    integer :: rc
    integer :: i,j,k, isn, jj, nmlat
    character(len=*), parameter :: subname = 'edyn3D_fldln2oplus_flg2opg_v2'

    if (mytid<ntask) then
       do k = 1,nhgt_fix

          if (mytid<ntask) then

             call ESMF_FieldGet(magfld%esmf_fld(k), localDe=0, farrayPtr=fptr2d, &
                  computationalLBound=lbnd2d, computationalUBound=ubnd2d, rc=rc)
             if (ESMF_LogFoundError(rc)) then
                call endrun(subname//' ESMF_FieldGet magfld%esmf_fld(k)')
             end if

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

          call ESMF_FieldRegrid(magfld%esmf_fld(k), oplsField2D, rh_mag2opls(k), &
               termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
          if (ESMF_LogFoundError(rc)) then
             call endrun(subname//' ESMF_FieldRegrid')
          end if

          call ESMF_FieldGet(field=oplsField2D, localDe=0, farrayPtr=fptr2d, &
               computationalLBound=lbnd2d, computationalUBound=ubnd2d, rc=rc)
          if (ESMF_LogFoundError(rc)) then
             call endrun(subname//' ESMF_FieldGet moplsField2D')
          end if

          do j = lbnd2d(2), ubnd2d(2)
             do i = lbnd2d(1), ubnd2d(1)
                f_tmp(i,j,k) = fptr2d(i,j)
             end do
          end do

       end do

       do i = lon0,lon1
          do j = lat0,lat1
             !vert interpolate...
             call lininterp(f_tmp(i,j,:), hgt_fix(:), nhgt_fix, &
                  opfld(i,j,:), opalt(i,j,:), nlevo )

          end do
       end do
    end if

  end subroutine edyn3D_fldln2oplus_flg2opg_v2

end module edyn3D_fldln2oplus_mod
