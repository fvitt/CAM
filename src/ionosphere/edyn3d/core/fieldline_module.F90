module fieldline_module

  use prec,only:rp
  use params_module,only:nhgt_fix,nhgt_fix_r,nmlat_h,nmlatS2_h

  implicit none

  integer,dimension(:),allocatable :: npts_p,npts_s1,npts_r
  integer,dimension(:),allocatable :: npts_s2
  integer,dimension(:),allocatable :: jmax_p,jmax_s1,jmax_s2,size_p,size_s1,size_s2
  integer,dimension(:),allocatable :: jmax_r,size_r

  real(kind=rp),dimension(:,:,:),allocatable :: qdlat_p,qdlat_s1
  real(kind=rp),dimension(:,:,:),allocatable :: qdlat_s2
  real(kind=rp),dimension(:,:,:),allocatable :: qdlat_r

  real(kind=rp),dimension(:,:,:,:,:),allocatable,target :: f3d,f3d_r
  real(kind=rp),dimension(:,:,:,:,:,:),allocatable,target :: uvec

  real(kind=rp),dimension(:,:,:,:),pointer :: &
    glat_p,glon_p,sinI_p,D_p,F_p,vmp_p,bmag_p,M3_p, &
    glat_s1,glon_s1,sinI_s1,D_s1,F_s1, &
    vmp_s1,bmag_s1,be3_s1,M1_s1, &
    d1d1_s1,d1d2_s1,d2d2_s1, &
    glat_s2,glon_s2,sinI_s2,D_s2,F_s2, &
    vmp_s2,bmag_s2,be3_s2,M2_s2, &
    d1d1_s2,d1d2_s2,d2d2_s2, &
    glat_r,glon_r,sinI_r,D_r,F_r,M3_r

  real(kind=rp),dimension(:,:,:,:,:),pointer :: &
    d1_s1,d2_s1,d3_s1,e1_s1,e2_s1,e3_s1, &
    d1_s2,d2_s2,d3_s2,e1_s2,e2_s2,e3_s2

endmodule fieldline_module
