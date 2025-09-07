module alloc_module
  use prec, only: rp
  use params_module, only: nmlat_h, nmlatS2_h, nhgt_fix, nhgt_fix_r
  use mpi_module, only: mlond0,mlond1,mlatd0,mlatd1
  use fieldline_module

  implicit none

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine alloc_fieldline(ierr)

! All grids include one latitude halo point on each side,
! but only internal grids (1<=j<=nmlat_h) are defined,
! the halo points at j==0 and j==nmlat_h+1 are not defined.

! Although S2 grids have equal latitudes with P/S1/R grids,
! the grid at j==nmlat_h is not defined.

    integer, intent(out) :: ierr

    ierr = 0

    allocate(f3d(nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1,26), stat=ierr)
    if (ierr /= 0) return
    allocate(f3d_r(nhgt_fix_r,2,mlatd0:mlatd1,mlond0:mlond1,4), stat=ierr)
    if (ierr /= 0) return
    allocate(uvec(3,nhgt_fix,2,mlatd0:mlatd1,mlond0:mlond1,12), stat=ierr)
    if (ierr /= 0) return

    sinI_p(1:,1:,mlatd0:,mlond0:)  => f3d(:,:,:,:,1)
    D_p(1:,1:,mlatd0:,mlond0:)     => f3d(:,:,:,:,2)
    F_p(1:,1:,mlatd0:,mlond0:)     => f3d(:,:,:,:,3)
    vmp_p(1:,1:,mlatd0:,mlond0:)   => f3d(:,:,:,:,4)
    bmag_p(1:,1:,mlatd0:,mlond0:)  => f3d(:,:,:,:,5)
    M3_p(1:,1:,mlatd0:,mlond0:)    => f3d(:,:,:,:,6)

    sinI_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,7)
    D_s1(1:,1:,mlatd0:,mlond0:)    => f3d(:,:,:,:,8)
    F_s1(1:,1:,mlatd0:,mlond0:)    => f3d(:,:,:,:,9)
    vmp_s1(1:,1:,mlatd0:,mlond0:)  => f3d(:,:,:,:,10)
    bmag_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,11)
    be3_s1(1:,1:,mlatd0:,mlond0:)  => f3d(:,:,:,:,12)
    M1_s1(1:,1:,mlatd0:,mlond0:)   => f3d(:,:,:,:,13)
    d1d1_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,14)
    d1d2_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,15)
    d2d2_s1(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,16)

! S2 grid only goes to nmlatS2_h
    sinI_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,17)
    D_s2(1:,1:,mlatd0:,mlond0:)    => f3d(:,:,:,:,18)
    F_s2(1:,1:,mlatd0:,mlond0:)    => f3d(:,:,:,:,19)
    vmp_s2(1:,1:,mlatd0:,mlond0:)  => f3d(:,:,:,:,20)
    bmag_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,21)
    be3_s2(1:,1:,mlatd0:,mlond0:)  => f3d(:,:,:,:,22)
    M2_s2(1:,1:,mlatd0:,mlond0:)   => f3d(:,:,:,:,23)
    d1d1_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,24)
    d1d2_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,25)
    d2d2_s2(1:,1:,mlatd0:,mlond0:) => f3d(:,:,:,:,26)

    sinI_r(1:,1:,mlatd0:,mlond0:)  => f3d_r(:,:,:,:,1)
    D_r(1:,1:,mlatd0:,mlond0:)     => f3d_r(:,:,:,:,2)
    F_r(1:,1:,mlatd0:,mlond0:)     => f3d_r(:,:,:,:,3)
    M3_r(1:,1:,mlatd0:,mlond0:)    => f3d_r(:,:,:,:,4)

    d1_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,1)
    d2_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,2)
    d3_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,3)
    e1_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,4)
    e2_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,5)
    e3_s1(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,6)

! S2 grid only goes to nmlatS2_h
    d1_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,7)
    d2_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,8)
    d3_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,9)
    e1_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,10)
    e2_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,11)
    e3_s2(1:,1:,1:,mlatd0:,mlond0:) => uvec(:,:,:,:,:,12)

  endsubroutine alloc_fieldline
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine dealloc_fieldline()

    if (allocated(f3d)) deallocate(f3d)
    if (allocated(f3d_r)) deallocate(f3d_r)
    if (allocated(uvec)) deallocate(uvec)

    nullify(sinI_p)
    nullify(D_p)
    nullify(F_p)
    nullify(vmp_p)
    nullify(bmag_p)
    nullify(M3_p)

    nullify(sinI_s1)
    nullify(D_s1)
    nullify(F_s1)
    nullify(vmp_s1)
    nullify(bmag_s1)
    nullify(be3_s1)
    nullify(M1_s1)
    nullify(d1d1_s1)
    nullify(d1d2_s1)
    nullify(d2d2_s1)


    nullify(sinI_s2)
    nullify(D_s2)
    nullify(F_s2)
    nullify(vmp_s2)
    nullify(bmag_s2)
    nullify(be3_s2)
    nullify(M2_s2)
    nullify(d1d1_s2)
    nullify(d1d2_s2)
    nullify(d2d2_s2)

    nullify(sinI_r)
    nullify(D_r)
    nullify(F_r)
    nullify(M3_r)

    nullify(d1_s1)
    nullify(d2_s1)
    nullify(d3_s1)
    nullify(e1_s1)
    nullify(e2_s1)
    nullify(e3_s1)


    nullify(d1_s2)
    nullify(d2_s2)
    nullify(d3_s2)
    nullify(e1_s2)
    nullify(e2_s2)
    nullify(e3_s2)

  end subroutine dealloc_fieldline
!-----------------------------------------------------------------------

endmodule alloc_module
