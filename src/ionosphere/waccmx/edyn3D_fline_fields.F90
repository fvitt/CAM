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
     type(ESMF_Field),       pointer :: esmf_fld(:) => null()
     type(ESMF_RouteHandle), pointer :: rhandle_phys2mag(:) => null()
     type(ESMF_RouteHandle), pointer :: rhandle_mag2phys(:) => null()
     type(fieldline_t),      pointer :: flines(:,:,:) => null()

  end type magfield_t

  type(magfield_t) :: Tn_p

  type(magfield_t) :: sigma_hal_s1
  type(magfield_t) :: sigma_ped_s1
  type(magfield_t) :: sigma_hal_s2
  type(magfield_t) :: sigma_ped_s2

contains

  subroutine edyn3D_fline_fields_alloc()

    use edyn3D_params, only: nmlat_h,nmlatS2_h, nptsp_total,nptss2_total
    use edyn3d_mpi, only: mlon0_p,mlon1_p
    use edyn3D_fieldline, only: fline_p,fline_s1,fline_s2
    use edyn3D_esmf_regrid, only: magField, rh_phys2mag, rh_mag2phys
    use edyn3D_esmf_regrid, only: magField_s1, rh_phys2mag_s1
    use edyn3D_esmf_regrid, only: magField_s2, rh_phys2mag_s2
    use infnan, only: nan, assignment(=)

    integer :: h,i,j,k

    Tn_p%mlon0 = mlon0_p
    Tn_p%mlon1 = mlon1_p
    Tn_p%nmlat_h = nmlat_h
    Tn_p%nptstot = nptsp_total
    Tn_p%rhandle_phys2mag => rh_phys2mag
    Tn_p%rhandle_mag2phys => rh_mag2phys
    Tn_p%esmf_fld => magField
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

    sigma_hal_s1%mlon0 = mlon0_p
    sigma_hal_s1%mlon1 = mlon1_p
    sigma_hal_s1%nmlat_h = nmlat_h
    sigma_hal_s1%nptstot = nptsp_total
    sigma_hal_s1%rhandle_phys2mag => rh_phys2mag_s1
    sigma_hal_s1%esmf_fld => magField_s1
    allocate(sigma_hal_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    sigma_ped_s1%mlon0 = mlon0_p
    sigma_ped_s1%mlon1 = mlon1_p
    sigma_ped_s1%nmlat_h = nmlat_h
    sigma_ped_s1%nptstot = nptsp_total
    sigma_ped_s1%rhandle_phys2mag => rh_phys2mag_s1
    sigma_ped_s1%esmf_fld => magField_s1
    allocate(sigma_ped_s1%flines(mlon0_p:mlon1_p,nmlat_h,2))

    do h = 1,2
       do j = 1,nmlat_h
          do i = mlon0_p,mlon1_p

             sigma_hal_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(sigma_hal_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             sigma_hal_s1%flines(i,j,h)%fld = nan

             sigma_ped_s1%flines(i,j,h)%npts = fline_s1(i,j,h)%npts
             allocate(sigma_ped_s1%flines(i,j,h)%fld(fline_s1(i,j,h)%npts))
             sigma_ped_s1%flines(i,j,h)%fld = nan

          end do
       end do
    end do

    sigma_hal_s2%mlon0 = mlon0_p
    sigma_hal_s2%mlon1 = mlon1_p
    sigma_hal_s2%nmlat_h = nmlatS2_h
    sigma_hal_s2%nptstot = nptss2_total
    sigma_hal_s2%rhandle_phys2mag => rh_phys2mag_s2
    sigma_hal_s2%esmf_fld => magField_s2
    allocate(sigma_hal_s2%flines(mlon0_p:mlon1_p,nmlat_h,2))

    sigma_ped_s2%mlon0 = mlon0_p
    sigma_ped_s2%mlon1 = mlon1_p
    sigma_ped_s2%nmlat_h = nmlatS2_h
    sigma_ped_s2%nptstot = nptss2_total
    sigma_ped_s2%rhandle_phys2mag => rh_phys2mag_s2
    sigma_ped_s2%esmf_fld => magField_s2
    allocate(sigma_ped_s2%flines(mlon0_p:mlon1_p,nmlat_h,2))

    do h = 1,2
       do j = 1,nmlatS2_h
          do i = mlon0_p,mlon1_p

             sigma_hal_s2%flines(i,j,h)%npts = fline_s2(i,j,h)%npts
             allocate(sigma_hal_s2%flines(i,j,h)%fld(fline_s2(i,j,h)%npts))
             sigma_hal_s2%flines(i,j,h)%fld = nan

             sigma_ped_s2%flines(i,j,h)%npts = fline_s2(i,j,h)%npts
             allocate(sigma_ped_s2%flines(i,j,h)%fld(fline_s2(i,j,h)%npts))
             sigma_ped_s2%flines(i,j,h)%fld = nan

          end do
       end do
    end do

  end subroutine edyn3D_fline_fields_alloc

end module edyn3D_fline_fields
