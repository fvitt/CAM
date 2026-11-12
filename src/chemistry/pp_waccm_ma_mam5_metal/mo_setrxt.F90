
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,105) = 0.000258_r8
      rate(:,:,106) = 0.085_r8
      rate(:,:,107) = 1.2e-10_r8
      rate(:,:,112) = 1.2e-10_r8
      rate(:,:,113) = 1e-20_r8
      rate(:,:,114) = 1.3e-16_r8
      rate(:,:,116) = 4.2e-13_r8
      rate(:,:,118) = 8e-14_r8
      rate(:,:,119) = 3.9e-17_r8
      rate(:,:,126) = 6.9e-12_r8
      rate(:,:,127) = 7.2e-11_r8
      rate(:,:,128) = 1.6e-12_r8
      rate(:,:,134) = 1.8e-12_r8
      rate(:,:,138) = 1.8e-12_r8
      rate(:,:,141) = 1.06e-05_r8
      rate(:,:,143) = 7e-11_r8
      rate(:,:,144) = 7e-13_r8
      rate(:,:,152) = 3.5e-12_r8
      rate(:,:,154) = 1.3e-11_r8
      rate(:,:,155) = 2.2e-11_r8
      rate(:,:,156) = 5e-11_r8
      rate(:,:,194) = 1.7e-13_r8
      rate(:,:,196) = 2.607e-10_r8
      rate(:,:,197) = 9.75e-11_r8
      rate(:,:,198) = 2.07e-10_r8
      rate(:,:,199) = 2.088e-10_r8
      rate(:,:,200) = 1.17e-10_r8
      rate(:,:,201) = 4.644e-11_r8
      rate(:,:,202) = 1.204e-10_r8
      rate(:,:,203) = 9.9e-11_r8
      rate(:,:,204) = 3.3e-12_r8
      rate(:,:,223) = 4.5e-11_r8
      rate(:,:,224) = 4.62e-10_r8
      rate(:,:,225) = 1.2e-10_r8
      rate(:,:,226) = 9e-11_r8
      rate(:,:,227) = 3e-11_r8
      rate(:,:,232) = 2.14e-11_r8
      rate(:,:,233) = 1.9e-10_r8
      rate(:,:,246) = 2.57e-10_r8
      rate(:,:,247) = 1.8e-10_r8
      rate(:,:,248) = 1.794e-10_r8
      rate(:,:,249) = 1.3e-10_r8
      rate(:,:,250) = 7.65e-11_r8
      rate(:,:,258) = 1.31e-10_r8
      rate(:,:,259) = 3.5e-11_r8
      rate(:,:,260) = 9e-12_r8
      rate(:,:,266) = 2.3e-12_r8
      rate(:,:,268) = 1.2e-11_r8
      rate(:,:,269) = 5.7e-11_r8
      rate(:,:,270) = 2.8e-11_r8
      rate(:,:,271) = 6.6e-11_r8
      rate(:,:,272) = 1.4e-11_r8
      rate(:,:,275) = 1.9e-12_r8
      rate(:,:,300) = 0.047_r8
      rate(:,:,301) = 7.7e-05_r8
      rate(:,:,302) = 0.171_r8
      rate(:,:,306) = 6e-11_r8
      rate(:,:,309) = 1e-12_r8
      rate(:,:,310) = 4e-10_r8
      rate(:,:,311) = 2e-10_r8
      rate(:,:,312) = 1e-10_r8
      rate(:,:,313) = 5e-16_r8
      rate(:,:,314) = 4.4e-10_r8
      rate(:,:,315) = 9e-10_r8
      rate(:,:,317) = 1.3e-10_r8
      rate(:,:,320) = 8e-10_r8
      rate(:,:,321) = 5e-12_r8
      rate(:,:,322) = 7e-10_r8
      rate(:,:,325) = 4.8e-10_r8
      rate(:,:,326) = 1e-10_r8
      rate(:,:,327) = 4e-10_r8
      rate(:,:,331) = 4.9e-12_r8
      rate(:,:,333) = 3.9e-11_r8
      rate(:,:,339) = 2.7e-9_r8
      rate(:,:,340) = 8.0e-10_r8
      rate(:,:,343) = 6.e-10_r8
      rate(:,:,344) = 6.e-10_r8
      rate(:,:,345) = 4.e-10_r8
      rate(:,:,346) = 1.e-11_r8
      rate(:,:,347) = 1.e-12_r8
      rate(:,:,348) = 5.e-12_r8
      rate(:,:,360) = 5.0e-12_r8
      rate(:,:,365) = 3e-11_r8
      rate(:,:,366) = 5.0e-11_r8
      rate(:,:,367) = 5.0e-11_r8
      rate(:,:,368) = 1.1e-9_r8
      rate(:,:,369) = 9.2e-10_r8
      rate(:,:,370) = 9.0e-10_r8
      rate(:,:,385) = 1.0e-12_r8
      rate(:,:,388) = 6.7e-12_r8
      rate(:,:,390) = 8.0e-14_r8
      rate(:,:,392) = 2.0e-10_r8
      rate(:,:,393) = 2.0e-12_r8
      rate(:,:,394) = 9.0e-10_r8
      rate(:,:,395) = 1.2e-9_r8
      rate(:,:,396) = 8.2e-10_r8
      rate(:,:,397) = 1.17e-9_r8
      rate(:,:,398) = 5.9e-10_r8
      rate(:,:,400) = 3.5e-12_r8
      rate(:,:,402) = 6.5e-10_r8
      rate(:,:,403) = 1.8e-10_r8
      rate(:,:,404) = 3.3e-10_r8
      rate(:,:,413) = 4.0e-10_r8
      rate(:,:,414) = 5.0e-13_r8
      rate(:,:,419) = 3.3e-10_r8
      rate(:,:,420) = 3.2e-10_r8
      rate(:,:,424) = 6.0e-10_r8
      rate(:,:,426) = 3.2e-10_r8
      rate(:,:,427) = 2.0e-10_r8
      rate(:,:,428) = 2.0e-10_r8
      rate(:,:,433) = 2.0e-10_r8
      rate(:,:,436) = 2.7e-12_r8
      rate(:,:,442) = 5.0e-12_r8
      rate(:,:,443) = 5.0e-12_r8
      rate(:,:,447) = 0.0e-11_r8
      rate(:,:,448) = 2.0e-11_r8
      rate(:,:,449) = 1.7e-11_r8
      rate(:,:,451) = 0.0e-12_r8
      rate(:,:,452) = 6.7e-12_r8
      rate(:,:,453) = 1.2e-11_r8
      rate(:,:,454) = 1.0e-10_r8
      rate(:,:,456) = 2.0e-10_r8
      rate(:,:,457) = 1.5e-10_r8
      rate(:,:,458) = 9.0e-8_r8
      rate(:,:,459) = 9.0e-8_r8
      rate(:,:,460) = 9.0e-8_r8
      rate(:,:,461) = 9.0e-8_r8
      rate(:,:,462) = 9.0e-8_r8
      rate(:,:,463) = 1.8e-9_r8
      rate(:,:,464) = 4.0e-9_r8
      rate(:,:,465) = 3.9e-10_r8
      rate(:,:,466) = 4.2e-11_r8
      rate(:,:,468) = 1.0e-10_r8
      rate(:,:,471) = 1.2e-10_r8
      rate(:,:,472) = 1.3e-9_r8
      rate(:,:,473) = 4.e-10_r8
      rate(:,:,479) = 3e-7_r8
      rate(:,:,481) = 3e-10_r8
      rate(:,:,493) = 2.7e-7_r8
      rate(:,:,494) = 9.4e-10_r8
      rate(:,:,495) = 3.2e-9_r8
      rate(:,:,524) = 1.0e-12_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,108) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      rate(:,:,109) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,110) = 2.64e-11_r8 * exp_fac(:,:)
      rate(:,:,111) = 6.6e-12_r8 * exp_fac(:,:)
      rate(:,:,115) = 3.6e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,117) = 1.8e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,120) = 3.5e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,121) = 8e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,124) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      rate(:,:,125) = 1.4e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,130) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,221) = 5.5e-12_r8 * exp_fac(:,:)
      rate(:,:,256) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,131) = 1e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,132) = 1.4e-10_r8 * exp( -470._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1800._r8 * itemp(:,:) )
      rate(:,:,133) = 2.8e-12_r8 * exp_fac(:,:)
      rate(:,:,501) = 2.6e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,135) = 4.8e-11_r8 * exp_fac(:,:)
      rate(:,:,219) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,136) = 1.8e-11_r8 * exp( 180._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -940._r8 * itemp(:,:) )
      rate(:,:,137) = 1.7e-12_r8 * exp_fac(:,:)
      rate(:,:,337) = 5.e-10_r8 * exp_fac(:,:)
      rate(:,:,142) = 4.5e-13_r8 * exp( 610._r8 * itemp(:,:) )
      rate(:,:,145) = 2.1e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,146) = 2.9e-12_r8 * exp_fac(:,:)
      rate(:,:,147) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,148) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,149) = 5.1e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,150) = 1.2e-13_r8 * exp_fac(:,:)
      rate(:,:,179) = 3e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 125._r8 * itemp(:,:) )
      rate(:,:,153) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,253) = 5.5e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,157) = 3.44e-12_r8 * exp_fac(:,:)
      rate(:,:,212) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,215) = 8.8e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,158) = 3e-12_r8 * exp_fac(:,:)
      rate(:,:,220) = 5.8e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,160) = 7.26e-11_r8 * exp_fac(:,:)
      rate(:,:,161) = 4.64e-11_r8 * exp_fac(:,:)
      rate(:,:,171) = 8.1e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,172) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:,:) )
      rate(:,:,173) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,174) = 1.1e-11_r8 * exp( -980._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,175) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,189) = 7.4e-12_r8 * exp_fac(:,:)
      rate(:,:,176) = 3.6e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,177) = 2.3e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,178) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,180) = 1e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,181) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,182) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,183) = 6.4e-12_r8 * exp_fac(:,:)
      rate(:,:,213) = 4.1e-13_r8 * exp_fac(:,:)
      rate(:,:,184) = 6.5e-12_r8 * exp( 135._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,186) = 3.6e-12_r8 * exp_fac(:,:)
      rate(:,:,235) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,187) = 1.2e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,188) = 2.8e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,190) = 6e-13_r8 * exp_fac(:,:)
      rate(:,:,210) = 1.5e-12_r8 * exp_fac(:,:)
      rate(:,:,218) = 1.9e-11_r8 * exp_fac(:,:)
      rate(:,:,191) = 1e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,192) = 1.8e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,193) = 3.4e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,195) = 3e-12_r8 * exp_fac(:,:)
      rate(:,:,229) = 1.4e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,207) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,234) = 6.3e-12_r8 * exp_fac(:,:)
      rate(:,:,208) = 4.8e-12_r8 * exp( -310._r8 * itemp(:,:) )
      rate(:,:,209) = 1.6e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,211) = 9.5e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,214) = 4.5e-12_r8 * exp( 460._r8 * itemp(:,:) )
      rate(:,:,217) = 1.9e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,222) = 1.2e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,228) = 1.6e-10_r8 * exp( -260._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,230) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,232) = 2.14e-11_r8 * exp_fac(:,:)
      rate(:,:,233) = 1.9e-10_r8 * exp_fac(:,:)
      rate(:,:,246) = 2.57e-10_r8 * exp_fac(:,:)
      rate(:,:,247) = 1.8e-10_r8 * exp_fac(:,:)
      rate(:,:,248) = 1.794e-10_r8 * exp_fac(:,:)
      rate(:,:,249) = 1.3e-10_r8 * exp_fac(:,:)
      rate(:,:,250) = 7.65e-11_r8 * exp_fac(:,:)
      rate(:,:,258) = 1.31e-10_r8 * exp_fac(:,:)
      rate(:,:,259) = 3.5e-11_r8 * exp_fac(:,:)
      rate(:,:,260) = 9e-12_r8 * exp_fac(:,:)
      rate(:,:,266) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,268) = 1.2e-11_r8 * exp_fac(:,:)
      rate(:,:,269) = 5.7e-11_r8 * exp_fac(:,:)
      rate(:,:,270) = 2.8e-11_r8 * exp_fac(:,:)
      rate(:,:,271) = 6.6e-11_r8 * exp_fac(:,:)
      rate(:,:,272) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,275) = 1.9e-12_r8 * exp_fac(:,:)
      rate(:,:,300) = 0.047_r8 * exp_fac(:,:)
      rate(:,:,301) = 7.7e-05_r8 * exp_fac(:,:)
      rate(:,:,302) = 0.171_r8 * exp_fac(:,:)
      rate(:,:,306) = 6e-11_r8 * exp_fac(:,:)
      rate(:,:,309) = 1e-12_r8 * exp_fac(:,:)
      rate(:,:,310) = 4e-10_r8 * exp_fac(:,:)
      rate(:,:,311) = 2e-10_r8 * exp_fac(:,:)
      rate(:,:,312) = 1e-10_r8 * exp_fac(:,:)
      rate(:,:,313) = 5e-16_r8 * exp_fac(:,:)
      rate(:,:,314) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,315) = 9e-10_r8 * exp_fac(:,:)
      rate(:,:,317) = 1.3e-10_r8 * exp_fac(:,:)
      rate(:,:,320) = 8e-10_r8 * exp_fac(:,:)
      rate(:,:,321) = 5e-12_r8 * exp_fac(:,:)
      rate(:,:,322) = 7e-10_r8 * exp_fac(:,:)
      rate(:,:,325) = 4.8e-10_r8 * exp_fac(:,:)
      rate(:,:,326) = 1e-10_r8 * exp_fac(:,:)
      rate(:,:,327) = 4e-10_r8 * exp_fac(:,:)
      rate(:,:,331) = 4.9e-12_r8 * exp_fac(:,:)
      rate(:,:,333) = 3.9e-11_r8 * exp_fac(:,:)
      rate(:,:,339) = 2.7e-9_r8 * exp_fac(:,:)
      rate(:,:,340) = 8.0e-10_r8 * exp_fac(:,:)
      rate(:,:,343) = 6.e-10_r8 * exp_fac(:,:)
      rate(:,:,344) = 6.e-10_r8 * exp_fac(:,:)
      rate(:,:,345) = 4.e-10_r8 * exp_fac(:,:)
      rate(:,:,346) = 1.e-11_r8 * exp_fac(:,:)
      rate(:,:,347) = 1.e-12_r8 * exp_fac(:,:)
      rate(:,:,348) = 5.e-12_r8 * exp_fac(:,:)
      rate(:,:,360) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,365) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,366) = 5.0e-11_r8 * exp_fac(:,:)
      rate(:,:,367) = 5.0e-11_r8 * exp_fac(:,:)
      rate(:,:,368) = 1.1e-9_r8 * exp_fac(:,:)
      rate(:,:,369) = 9.2e-10_r8 * exp_fac(:,:)
      rate(:,:,370) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,385) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,388) = 6.7e-12_r8 * exp_fac(:,:)
      rate(:,:,390) = 8.0e-14_r8 * exp_fac(:,:)
      rate(:,:,392) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,393) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,394) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,395) = 1.2e-9_r8 * exp_fac(:,:)
      rate(:,:,396) = 8.2e-10_r8 * exp_fac(:,:)
      rate(:,:,397) = 1.17e-9_r8 * exp_fac(:,:)
      rate(:,:,398) = 5.9e-10_r8 * exp_fac(:,:)
      rate(:,:,400) = 3.5e-12_r8 * exp_fac(:,:)
      rate(:,:,402) = 6.5e-10_r8 * exp_fac(:,:)
      rate(:,:,403) = 1.8e-10_r8 * exp_fac(:,:)
      rate(:,:,404) = 3.3e-10_r8 * exp_fac(:,:)
      rate(:,:,413) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,414) = 5.0e-13_r8 * exp_fac(:,:)
      rate(:,:,419) = 3.3e-10_r8 * exp_fac(:,:)
      rate(:,:,420) = 3.2e-10_r8 * exp_fac(:,:)
      rate(:,:,424) = 6.0e-10_r8 * exp_fac(:,:)
      rate(:,:,426) = 3.2e-10_r8 * exp_fac(:,:)
      rate(:,:,427) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,428) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,433) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,436) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,442) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,443) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,447) = 0.0e-11_r8 * exp_fac(:,:)
      rate(:,:,448) = 2.0e-11_r8 * exp_fac(:,:)
      rate(:,:,449) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,451) = 0.0e-12_r8 * exp_fac(:,:)
      rate(:,:,452) = 6.7e-12_r8 * exp_fac(:,:)
      rate(:,:,453) = 1.2e-11_r8 * exp_fac(:,:)
      rate(:,:,454) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,456) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,457) = 1.5e-10_r8 * exp_fac(:,:)
      rate(:,:,458) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,459) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,460) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,461) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,462) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,463) = 1.8e-9_r8 * exp_fac(:,:)
      rate(:,:,464) = 4.0e-9_r8 * exp_fac(:,:)
      rate(:,:,465) = 3.9e-10_r8 * exp_fac(:,:)
      rate(:,:,466) = 4.2e-11_r8 * exp_fac(:,:)
      rate(:,:,468) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,471) = 1.2e-10_r8 * exp_fac(:,:)
      rate(:,:,472) = 1.3e-9_r8 * exp_fac(:,:)
      rate(:,:,473) = 4.e-10_r8 * exp_fac(:,:)
      rate(:,:,479) = 3e-7_r8 * exp_fac(:,:)
      rate(:,:,481) = 3e-10_r8 * exp_fac(:,:)
      rate(:,:,493) = 2.7e-7_r8 * exp_fac(:,:)
      rate(:,:,494) = 9.4e-10_r8 * exp_fac(:,:)
      rate(:,:,495) = 3.2e-9_r8 * exp_fac(:,:)
      rate(:,:,524) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,231) = 6e-12_r8 * exp( 400._r8 * itemp(:,:) )
      rate(:,:,236) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:,:) )
      rate(:,:,237) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:,:) )
      rate(:,:,238) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      rate(:,:,239) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1200._r8 * itemp(:,:) )
      rate(:,:,240) = 1.96e-12_r8 * exp_fac(:,:)
      rate(:,:,379) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,241) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,242) = 9e-13_r8 * exp( -360._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,243) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,252) = 3.4e-11_r8 * exp_fac(:,:)
      rate(:,:,244) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,245) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:,:) )
      rate(:,:,251) = 6e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,254) = 4.1e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,255) = 2.8e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,257) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,262) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,263) = 1.1e-11_r8 * exp( -280._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2200._r8 * itemp(:,:) )
      rate(:,:,264) = 2.1e-11_r8 * exp_fac(:,:)
      rate(:,:,511) = 1.4e-9_r8 * exp_fac(:,:)
      rate(:,:,265) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:,:) )
      rate(:,:,273) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:,:) )
      rate(:,:,274) = 3.4e-12_r8 * exp( -1100._r8 * itemp(:,:) )
      rate(:,:,276) = 2.6e-11_r8 * exp( 330._r8 * itemp(:,:) )
      rate(:,:,328) = 1.1e-9_r8 * exp( -116._r8 * itemp(:,:) )
      rate(:,:,330) = 3.2e-10_r8 * exp( -550._r8 * itemp(:,:) )
      rate(:,:,332) = 5.06e-10_r8 * exp( -240._r8 * itemp(:,:) )
      rate(:,:,354) = 2.94e-10_r8 * exp( -174._r8 * itemp(:,:) )
      rate(:,:,355) = 4.6e-10_r8 * exp( -350._r8 * itemp(:,:) )
      rate(:,:,356) = 3.0e-10_r8 * exp( -177._r8 * itemp(:,:) )
      rate(:,:,357) = 1.4e-10_r8 * exp( -580._r8 * itemp(:,:) )
      rate(:,:,358) = 4.4e-10_r8 * exp( -170._r8 * itemp(:,:) )
      rate(:,:,359) = 2.3e-10_r8 * exp( -2310._r8 * itemp(:,:) )
      rate(:,:,361) = 3.3e-10_r8 * exp( -302._r8 * itemp(:,:) )
      rate(:,:,362) = 3.0e-10_r8 * exp( -796._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -600._r8 * itemp(:,:) )
      rate(:,:,363) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,391) = 1.0e-11_r8 * exp_fac(:,:)
      rate(:,:,364) = 7.6e-10_r8 * exp( -240.5_r8 * itemp(:,:) )
      rate(:,:,380) = 2.28e-10_r8 * exp( -139._r8 * itemp(:,:) )
      rate(:,:,382) = 2.19e-10_r8 * exp( -548._r8 * itemp(:,:) )
      rate(:,:,435) = 8.23e-10_r8 * exp( -192._r8 * itemp(:,:) )
      rate(:,:,437) = 1.1e-9_r8 * exp( -421._r8 * itemp(:,:) )
      rate(:,:,438) = 5.7e-10_r8 * exp( -267._r8 * itemp(:,:) )
      rate(:,:,440) = 4.4e-11_r8 * exp( -202._r8 * itemp(:,:) )
      rate(:,:,446) = 7.0e-10_r8 * exp( -4017._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -120._r8 * itemp(:,:) )
      rate(:,:,482) = 1.15e-9_r8 * exp_fac(:,:)
      rate(:,:,483) = 2.00e-10_r8 * exp_fac(:,:)
      rate(:,:,486) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,487) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,488) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,489) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,490) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,485) = 6.9e-10_r8 * exp( -385._r8 * itemp(:,:) )
      rate(:,:,492) = 4.5e-11_r8 * exp( -3590._r8 * itemp(:,:) )
      rate(:,:,497) = 2.8e-8_r8 * exp( -1680._r8 * itemp(:,:) )
      rate(:,:,499) = 1.5e-9_r8 * exp( -820._r8 * itemp(:,:) )
      rate(:,:,507) = 2.8e-10_r8 * exp( -2220._r8 * itemp(:,:) )
      rate(:,:,509) = 5.0e-10_r8 * exp( -752._r8 * itemp(:,:) )
      rate(:,:,514) = 1.0e-9_r8 * exp( -873._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 5.3e-32_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 9.5e-11_r8 * itemp(:,:)**(-0.4_r8)
      call jpl( rate(1,1,129), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.9e-31_r8 * itemp(:,:)**1._r8
      kinf(:,:) = 2.6e-11_r8
      call jpl( rate(1,1,139), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.5e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,151), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3e-11_r8
      call jpl( rate(1,1,159), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.9e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 4e-12_r8 * itemp(:,:)**0.3_r8
      call jpl( rate(1,1,162), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.4e-30_r8 * itemp(:,:)**3._r8
      kinf(:,:) = 1.6e-12_r8 * itemp(:,:)**(-0.1_r8)
      call jpl( rate(1,1,163), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.8e-30_r8 * itemp(:,:)**3._r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,164), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.8e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,185), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.9e-32_r8 * itemp(:,:)**3.6_r8
      kinf(:,:) = 3.7e-12_r8 * itemp(:,:)**1.6_r8
      call jpl( rate(1,1,205), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.2e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,216), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.9e-31_r8 * itemp(:,:)**4.1_r8
      kinf(:,:) = 1.7e-12_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,267), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)

      rate(:,:kbot,113) = 1e-20_r8
      rate(:,:kbot,114) = 1.3e-16_r8
      rate(:,:kbot,118) = 8e-14_r8
      rate(:,:kbot,119) = 3.9e-17_r8
      rate(:,:kbot,126) = 6.9e-12_r8
      rate(:,:kbot,143) = 7e-11_r8
      rate(:,:kbot,144) = 7e-13_r8
      rate(:,:kbot,300) = 0.047_r8
      rate(:,:kbot,301) = 7.7e-05_r8
      rate(:,:kbot,302) = 0.171_r8
      rate(:,:kbot,306) = 6e-11_r8
      rate(:,:kbot,309) = 1e-12_r8
      rate(:,:kbot,310) = 4e-10_r8
      rate(:,:kbot,311) = 2e-10_r8
      rate(:,:kbot,312) = 1e-10_r8
      rate(:,:kbot,314) = 4.4e-10_r8
      rate(:,:kbot,317) = 1.3e-10_r8
      rate(:,:kbot,320) = 8e-10_r8
      rate(:,:kbot,321) = 5e-12_r8
      rate(:,:kbot,322) = 7e-10_r8
      rate(:,:kbot,325) = 4.8e-10_r8
      rate(:,:kbot,326) = 1e-10_r8
      rate(:,:kbot,327) = 4e-10_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,109) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,110) = 2.64e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,111) = 6.6e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,115) = 3.6e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,117) = 1.8e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,120) = 3.5e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,121) = 8e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,130) = 3e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,131) = 1e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,132) = 1.4e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,135) = 4.8e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,136) = 1.8e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,137) = 1.7e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,145) = 2.1e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,149) = 5.1e-12_r8 * exp( 210._r8 * itemp(:,:) )
      rate(:,:kbot,157) = 3.44e-12_r8 * exp( 260._r8 * itemp(:,:) )
      rate(:,:kbot,158) = 3e-12_r8 * exp( -1500._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 5.3e-32_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 9.5e-11_r8 * itemp(:,:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,129) = wrk(:,:)











      end subroutine setrxt_hrates

      end module mo_setrxt
