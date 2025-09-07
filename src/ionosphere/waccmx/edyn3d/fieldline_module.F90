module fieldline_module

  use prec,only:rp

  implicit none

  real(kind=rp),dimension(:,:,:,:,:),allocatable,target :: f3d,f3d_r
  real(kind=rp),dimension(:,:,:,:,:,:),allocatable,target :: uvec

  real(kind=rp),dimension(:,:,:,:),pointer :: &
       sinI_p=>null(), &
       D_p=>null(), &
       F_p=>null(), &
       vmp_p=>null(), &
       bmag_p=>null(), &
       M3_p=>null(), &
       sinI_s1=>null(), &
       D_s1=>null(), &
       F_s1=>null(), &
       vmp_s1=>null(), &
       bmag_s1=>null(), &
       be3_s1=>null(), &
       M1_s1=>null(), &
       d1d1_s1=>null(), &
       d1d2_s1=>null(), &
       d2d2_s1=>null(), &
       sinI_s2=>null(), &
       D_s2=>null(), &
       F_s2=>null(), &
       vmp_s2=>null(), &
       bmag_s2=>null(), &
       be3_s2=>null(), &
       M2_s2=>null(), &
       d1d1_s2=>null(), &
       d1d2_s2=>null(), &
       d2d2_s2=>null(), &
       sinI_r=>null(), &
       D_r=>null(), &
       F_r=>null(), &
       M3_r=>null()

  real(kind=rp),dimension(:,:,:,:,:),pointer :: &
       d1_s1=>null(), &
       d2_s1=>null(), &
       d3_s1=>null(), &
       e1_s1=>null(), &
       e2_s1=>null(), &
       e3_s1=>null(), &
       d1_s2=>null(), &
       d2_s2=>null(), &
       d3_s2=>null(), &
       e1_s2=>null(), &
       e2_s2=>null(), &
       e3_s2=>null()

endmodule fieldline_module
