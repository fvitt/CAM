module edyn3D_glblslv_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none

  private
  public :: edyn3D_glblslv_poten

contains

  subroutine edyn3D_glblslv_poten(fline_p, fline_s1, fline_s2)

    use edyn3D_fieldline, only: fieldline_p,fieldline_s1,fieldline_s2
    use edyn3d_mpi, only: mytid, mlon0_p, mlon1_p, mp_gather_edyn3D, mp_scatter_edyn3D
    use edyn3D_params, only: nmlon, nmlat_T1,nmlat_h,nmlatS2_h,nhgt_fix, poten_hl
    use edyn3D_calculate_coefs, only: edyn3D_calculate_coef, edyn3D_calculate_coef_ns2
    use edyn3D_calculate_coefs, only: edyn3D_calculate_coef_ns, edyn3D_calculate_bij
    use edyn3D_serial_solver, only: linear_system

    ! args

    type(fieldline_p), intent(inout) :: fline_p(mlon0_p-1:mlon1_p+1,nmlat_h,2)
    type(fieldline_s1),intent(in) :: fline_s1(mlon0_p-1:mlon1_p+1,nmlat_h,2)
    type(fieldline_s2),intent(in) :: fline_s2(mlon0_p-1:mlon1_p+1,nmlatS2_h,2)

    ! local vars

    integer :: i,j,jj, astat
    integer :: isn, ncnt
    real(r8) :: coef(mlon0_p:mlon1_p,nmlat_h,nhgt_fix,10,2)
    real(r8) :: coef_ns(mlon0_p:mlon1_p,nmlat_T1,10)
    real(r8) :: coef_ns2(mlon0_p:mlon1_p,nmlat_h,2,10)
    real(r8) :: bij(mlon0_p:mlon1_p,nmlat_h)
    real(r8) :: fmsub(mlon0_p:mlon1_p,nmlat_h,1)
    real(r8) :: fmsub4(mlon0_p:mlon1_p,nmlat_h,4)

    real(r8),allocatable :: fmglb(:,:,:)
    real(r8),allocatable :: fmglb_T1(:,:,:)
    real(r8),allocatable :: fmglb4(:,:,:)

    real(r8),allocatable :: bij_glb(:,:)
    real(r8),allocatable :: pot_glb(:,:,:)
    real(r8),allocatable :: fac_hl_glb(:,:,:)
    real(r8),allocatable :: coef_ns_glb(:,:,:)
    real(r8),allocatable :: pot_hl_glb(:,:,:)

    character(len=*), parameter :: prefix = 'edyn3D_glblslv_poten: '

    call edyn3D_calculate_coef(fline_p,fline_s1,fline_s2,coef)

    call edyn3D_calculate_coef_ns2(coef,coef_ns2)

    call edyn3D_calculate_coef_ns(coef_ns2,coef_ns)

    call edyn3D_calculate_bij(coef_ns2,bij)

    allocate(fmglb(nmlon,nmlat_h,1), stat=astat)
    if (astat/=0) then
       call endrun(prefix//'fmglb array allocation failed')
    end if

    fmsub(mlon0_p:mlon1_p,:,1) = bij(mlon0_p:mlon1_p,:)
    call mp_gather_edyn3D(fmsub,mlon0_p,mlon1_p,fmglb,nmlon,nmlat_h,1)

    if (mytid==0) then
       allocate(bij_glb(nmlat_h,nmlon), stat=astat)
       if (astat/=0) then
          call endrun(prefix//'bij_glb array allocation failed')
       end if
       do i = 1,nmlon
          do j = 1,nmlat_h
             bij_glb(j,i) = fmglb(i,j,1)
          enddo
       end do
    end if

    deallocate(fmglb)


    allocate(fmglb_T1(nmlon,nmlat_T1,10), stat=astat)
    if (astat/=0) then
       call endrun(prefix//'fmglb_T1 array allocation failed')
    end if

    call mp_gather_edyn3D(coef_ns,mlon0_p,mlon1_p,fmglb_T1,nmlon,nmlat_T1,10)

    if (mytid==0) then
       allocate(coef_ns_glb(10,nmlat_T1,nmlon), stat=astat)
       if (astat/=0) then
          call endrun(prefix//'coef_ns_glb array allocation failed')
       end if

       do i = 1,nmlon
          do j = 1,nmlat_T1
             coef_ns_glb(:,j,i) = fmglb_T1(i,j,:)
          enddo
       enddo
    end if

    deallocate(fmglb_T1)

    if (mytid==0) then

       allocate(pot_hl_glb(2,nmlat_h,nmlon), stat=astat)
       if (astat/=0) then
          call endrun(prefix//'pot_hl_glb array allocation failed')
       end if

       do i = 1,nmlon
          do j = 1,nmlat_h
             do isn = 1,2

                if (isn==1) then
                   jj = j
                else
                   jj = nmlat_T1 - j + 1
                endif

                pot_hl_glb(isn,j,i) = poten_hl(i,jj)

             end do
          end do
       end do

       allocate(pot_glb(2,nmlat_h,nmlon), stat=astat)
       if (astat/=0) then
          call endrun(prefix//'pot_glb array allocation failed')
       end if
       allocate(fac_hl_glb(2,nmlat_h,nmlon), stat=astat)
       if (astat/=0) then
          call endrun(prefix//'fac_hl_glb array allocation failed')
       end if

       call linear_system(bij_glb,pot_hl_glb,fac_hl_glb,coef_ns_glb,pot_glb)

       deallocate(bij_glb,pot_hl_glb,coef_ns_glb)

    end if

    allocate(fmglb4(nmlon,nmlat_h,4), stat=astat)
    if (astat/=0) then
       call endrun(prefix//'fmglb4 array allocation failed')
    end if

    if (mytid==0) then

       do i = 1,nmlon
          do j = 1,nmlat_h
             fmglb4(i,j,1) = pot_glb(1,j,i)
             fmglb4(i,j,2) = pot_glb(2,j,i)
             fmglb4(i,j,3) = fac_hl_glb(1,j,i)
             fmglb4(i,j,4) = fac_hl_glb(2,j,i)
          enddo
       enddo

       deallocate(pot_glb,fac_hl_glb)

    end if ! task 0

    call mp_scatter_edyn3D(fmglb4,mlon0_p,mlon1_p,fmsub4,nmlon,nmlat_h,4)

    deallocate(fmglb4)

    do i=mlon0_p,mlon1_p ! loop over task longitudes
       do j=1,nmlat_h ! loop over all latitudes in one hemisphere
          fline_p(i,j,1)%pot    = fmsub4(i,j,1)
          fline_p(i,j,2)%pot    = fmsub4(i,j,2)
          fline_p(i,j,1)%fac_hl = fmsub4(i,j,3)
          fline_p(i,j,2)%fac_hl = fmsub4(i,j,4)
       end do
    end do

  end subroutine edyn3D_glblslv_poten


end module edyn3D_glblslv_mod
