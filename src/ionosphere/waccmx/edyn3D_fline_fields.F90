module edyn3D_fline_fields
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ESMF, only: ESMF_RouteHandle, ESMF_Field

  implicit none

  type fieldline_t
     integer :: npts
     real(r8), allocatable :: fld(:) ! k dim
  end type fieldline_t

  type magfield_t

     integer :: mlon0 = -huge(1)
     integer :: mlon1 = -huge(1)
     integer :: nmlat_h = -huge(1)
     integer :: nptstot = 0
     type(ESMF_Field),       pointer :: esmf_fld_src(:) => null()
     type(ESMF_Field),       pointer :: esmf_fld_des(:) => null()
     type(ESMF_RouteHandle), pointer :: rhandle_phys2mag(:) => null()
     type(ESMF_RouteHandle), pointer :: rhandle_mag2phys(:) => null()
     type(fieldline_t),      pointer :: flines(:,:,:) => null()
     character(len=32) :: name = ' '

  end type magfield_t

  type(magfield_t) :: Tn_p

  type(magfield_t) :: sigma_hal_s1
  type(magfield_t) :: sigma_ped_s1
  type(magfield_t) :: sigma_hal_s2
  type(magfield_t) :: sigma_ped_s2

  type(magfield_t) :: height_s1
  type(magfield_t) :: height_s2

  type(magfield_t) :: un_s1
  type(magfield_t) :: vn_s1
  type(magfield_t) :: un_s2
  type(magfield_t) :: vn_s2

  type(magfield_t) :: IonU_s1
  type(magfield_t) :: IonV_s1
  type(magfield_t) :: IonW_s1

contains

  subroutine edyn3D_fline_fields_alloc()

    use edyn3D_params, only: nmlat_h,nmlatS2_h, nptsp_total,nptss2_total, nmlat_T1, nmlat_T2
    use edyn3d_mpi, only: mlon0_p,mlon1_p
    use edyn3D_fieldline, only: fline_p,fline_s1,fline_s2

    use edyn3D_esmf_regrid, only: magField_p_src
    use edyn3D_esmf_regrid, only: magField_p_des
    use edyn3D_esmf_regrid, only: magField_s1_src
    use edyn3D_esmf_regrid, only: magField_s1_des
    use edyn3D_esmf_regrid, only: magField_s2

    use edyn3D_esmf_regrid, only: rh_phys2mag_p
    use edyn3D_esmf_regrid, only: rh_phys2mag_s1
    use edyn3D_esmf_regrid, only: rh_phys2mag_s2
    use edyn3D_esmf_regrid, only: rh_mag_p2phys, rh_mag_p2oplus
    use edyn3D_esmf_regrid, only: rh_s1mag2phys

    use infnan, only: nan, assignment(=)

    integer :: h,i,j,k

    Tn_p%name = 'Tn_mag'
    Tn_p%mlon0 = mlon0_p
    Tn_p%mlon1 = mlon1_p
    Tn_p%nmlat_h = nmlat_h
    Tn_p%nptstot = nptsp_total
    Tn_p%rhandle_phys2mag => rh_phys2mag_p
    Tn_p%rhandle_mag2phys => rh_mag_p2phys
    Tn_p%esmf_fld_src => magField_p_src
    Tn_p%esmf_fld_des => magField_p_des
    allocate(Tn_p%flines(mlon0_p:mlon1_p,nmlat_h,2))

    do h = 1,2
       do j = 1,nmlat_h
          do i = mlon0_p,mlon1_p

             Tn_p%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(Tn_p%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             Tn_p%flines(i,j,h)%fld = nan

          end do
       end do
    end do

    height_s1%name = 'height_s1'
    height_s1%mlon0 = mlon0_p
    height_s1%mlon1 = mlon1_p
    height_s1%nmlat_h = nmlat_h
    height_s1%nptstot = nptsp_total
    height_s1%rhandle_phys2mag => rh_phys2mag_s1
    height_s1%esmf_fld_des => magField_s1_des
    height_s1%esmf_fld_src => null()
    allocate(height_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    sigma_hal_s1%name = 'sigma_hal_s1'
    sigma_hal_s1%mlon0 = mlon0_p
    sigma_hal_s1%mlon1 = mlon1_p
    sigma_hal_s1%nmlat_h = nmlat_h
    sigma_hal_s1%nptstot = nptsp_total
    sigma_hal_s1%rhandle_phys2mag => rh_phys2mag_s1
    sigma_hal_s1%esmf_fld_des => magField_s1_des
    sigma_hal_s1%esmf_fld_src => null()
    allocate(sigma_hal_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    sigma_ped_s1%name = 'sigma_ped_s1'
    sigma_ped_s1%mlon0 = mlon0_p
    sigma_ped_s1%mlon1 = mlon1_p
    sigma_ped_s1%nmlat_h = nmlat_h
    sigma_ped_s1%nptstot = nptsp_total
    sigma_ped_s1%rhandle_phys2mag => rh_phys2mag_s1
    sigma_ped_s1%esmf_fld_des => magField_s1_des
    sigma_ped_s1%esmf_fld_src => null()
    allocate(sigma_ped_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    un_s1%name = 'un_s1'
    un_s1%mlon0 = mlon0_p
    un_s1%mlon1 = mlon1_p
    un_s1%nmlat_h = nmlat_h
    un_s1%nptstot = nptsp_total
    un_s1%rhandle_phys2mag => rh_phys2mag_s1
    un_s1%esmf_fld_des => magField_s1_des
    un_s1%esmf_fld_src => null()
    allocate(un_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    vn_s1%name = 'vn_s1'
    vn_s1%mlon0 = mlon0_p
    vn_s1%mlon1 = mlon1_p
    vn_s1%nmlat_h = nmlat_h
    vn_s1%nptstot = nptsp_total
    vn_s1%rhandle_phys2mag => rh_phys2mag_s1
    vn_s1%esmf_fld_des => magField_s1_des
    vn_s1%esmf_fld_src => null()
    allocate(vn_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    IonU_s1%name = 'IonU_s1'
    IonU_s1%mlon0 = mlon0_p
    IonU_s1%mlon1 = mlon1_p
    IonU_s1%nmlat_h = nmlat_h
    IonU_s1%nptstot = nptsp_total
    IonU_s1%rhandle_mag2phys => rh_s1mag2phys
    IonU_s1%rhandle_phys2mag => rh_phys2mag_s1
    IonU_s1%esmf_fld_des => magField_s1_des
    IonU_s1%esmf_fld_src => magField_s1_src
    allocate(IonU_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    IonV_s1%name = 'IonV_s1'
    IonV_s1%mlon0 = mlon0_p
    IonV_s1%mlon1 = mlon1_p
    IonV_s1%nmlat_h = nmlat_h
    IonV_s1%nptstot = nptsp_total
    IonV_s1%rhandle_mag2phys => rh_s1mag2phys
    IonV_s1%rhandle_phys2mag => rh_phys2mag_s1
    IonV_s1%esmf_fld_des => magField_s1_des
    IonV_s1%esmf_fld_src => magField_s1_src
    allocate(IonV_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    IonW_s1%name = 'IonW_s1'
    IonW_s1%mlon0 = mlon0_p
    IonW_s1%mlon1 = mlon1_p
    IonW_s1%nmlat_h = nmlat_h
    IonW_s1%nptstot = nptsp_total
    IonW_s1%rhandle_mag2phys => rh_s1mag2phys
    IonW_s1%rhandle_phys2mag => rh_phys2mag_s1
    IonW_s1%esmf_fld_des => magField_s1_des
    IonW_s1%esmf_fld_src => magField_s1_src
    allocate(IonW_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    do h = 1,2
       do j = 1,nmlat_h
          do i = mlon0_p,mlon1_p

             height_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(height_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             height_s1%flines(i,j,h)%fld = nan

             sigma_hal_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(sigma_hal_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             sigma_hal_s1%flines(i,j,h)%fld = nan

             sigma_ped_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(sigma_ped_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             sigma_ped_s1%flines(i,j,h)%fld = nan

             un_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(un_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             un_s1%flines(i,j,h)%fld = nan

             vn_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(vn_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             vn_s1%flines(i,j,h)%fld = nan

             IonU_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(IonU_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             IonU_s1%flines(i,j,h)%fld = nan

             IonV_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(IonV_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             IonV_s1%flines(i,j,h)%fld = nan

             IonW_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(IonW_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             IonW_s1%flines(i,j,h)%fld = nan

          end do
       end do
    end do

    height_s2%name = 'height_s2'
    height_s2%mlon0 = mlon0_p
    height_s2%mlon1 = mlon1_p
    height_s2%nmlat_h = nmlatS2_h
    height_s2%nptstot = nptss2_total
    height_s2%rhandle_phys2mag => rh_phys2mag_s2
    height_s2%esmf_fld_des => magField_s2
    allocate(height_s2%flines(mlon0_p:mlon1_p,nmlatS2_h,2))

    sigma_hal_s2%name = 'sigma_hal_s2'
    sigma_hal_s2%mlon0 = mlon0_p
    sigma_hal_s2%mlon1 = mlon1_p
    sigma_hal_s2%nmlat_h = nmlatS2_h
    sigma_hal_s2%nptstot = nptss2_total
    sigma_hal_s2%rhandle_phys2mag => rh_phys2mag_s2
    sigma_hal_s2%esmf_fld_des => magField_s2
    allocate(sigma_hal_s2%flines(mlon0_p:mlon1_p,nmlatS2_h,2))

    sigma_ped_s2%name = 'sigma_ped_s2'
    sigma_ped_s2%mlon0 = mlon0_p
    sigma_ped_s2%mlon1 = mlon1_p
    sigma_ped_s2%nmlat_h = nmlatS2_h
    sigma_ped_s2%nptstot = nptss2_total
    sigma_ped_s2%rhandle_phys2mag => rh_phys2mag_s2
    sigma_ped_s2%esmf_fld_des => magField_s2
    allocate(sigma_ped_s2%flines(mlon0_p:mlon1_p,nmlatS2_h,2))

    un_s2%name = 'un_s2'
    un_s2%mlon0 = mlon0_p
    un_s2%mlon1 = mlon1_p
    un_s2%nmlat_h = nmlatS2_h
    un_s2%nptstot = nptss2_total
    un_s2%rhandle_phys2mag => rh_phys2mag_s2
    un_s2%esmf_fld_des => magField_s2
    allocate(un_s2%flines(mlon0_p:mlon1_p,nmlatS2_h,2))

    vn_s2%name = 'vn_s2'
    vn_s2%mlon0 = mlon0_p
    vn_s2%mlon1 = mlon1_p
    vn_s2%nmlat_h = nmlatS2_h
    vn_s2%nptstot = nptss2_total
    vn_s2%rhandle_phys2mag => rh_phys2mag_s2
    vn_s2%esmf_fld_des => magField_s2
    allocate(vn_s2%flines(mlon0_p:mlon1_p,nmlatS2_h,2))

    do h = 1,2
       do j = 1,nmlatS2_h
          do i = mlon0_p,mlon1_p

             height_s2%flines(i,j,h)%npts = fline_s2(i,j,h)%npts
             allocate(height_s2%flines(i,j,h)%fld(fline_s2(i,j,h)%npts))
             height_s2%flines(i,j,h)%fld = nan

             sigma_hal_s2%flines(i,j,h)%npts = fline_s2(i,j,h)%npts
             allocate(sigma_hal_s2%flines(i,j,h)%fld(fline_s2(i,j,h)%npts))
             sigma_hal_s2%flines(i,j,h)%fld = nan

             sigma_ped_s2%flines(i,j,h)%npts = fline_s2(i,j,h)%npts
             allocate(sigma_ped_s2%flines(i,j,h)%fld(fline_s2(i,j,h)%npts))
             sigma_ped_s2%flines(i,j,h)%fld = nan

             un_s2%flines(i,j,h)%npts = fline_s2(i,j,h)%npts
             allocate(un_s2%flines(i,j,h)%fld(fline_s2(i,j,h)%npts))
             un_s2%flines(i,j,h)%fld = nan

             vn_s2%flines(i,j,h)%npts = fline_s2(i,j,h)%npts
             allocate(vn_s2%flines(i,j,h)%fld(fline_s2(i,j,h)%npts))
             vn_s2%flines(i,j,h)%fld = nan

          end do
       end do
    end do

  end subroutine edyn3D_fline_fields_alloc

end module edyn3D_fline_fields
