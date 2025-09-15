
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
      rate(:,:,142) = 7e-13_r8
      rate(:,:,143) = 5e-12_r8
      rate(:,:,152) = 3.5e-12_r8
      rate(:,:,154) = 1.3e-11_r8
      rate(:,:,155) = 2.2e-11_r8
      rate(:,:,156) = 5e-11_r8
      rate(:,:,191) = 1.7e-13_r8
      rate(:,:,193) = 2.607e-10_r8
      rate(:,:,194) = 9.75e-11_r8
      rate(:,:,195) = 2.07e-10_r8
      rate(:,:,196) = 2.088e-10_r8
      rate(:,:,197) = 1.17e-10_r8
      rate(:,:,198) = 4.644e-11_r8
      rate(:,:,199) = 1.204e-10_r8
      rate(:,:,200) = 9.9e-11_r8
      rate(:,:,201) = 3.3e-12_r8
      rate(:,:,220) = 4.5e-11_r8
      rate(:,:,221) = 4.62e-10_r8
      rate(:,:,222) = 1.2e-10_r8
      rate(:,:,223) = 9e-11_r8
      rate(:,:,224) = 3e-11_r8
      rate(:,:,229) = 2.14e-11_r8
      rate(:,:,230) = 1.9e-10_r8
      rate(:,:,243) = 2.57e-10_r8
      rate(:,:,244) = 1.8e-10_r8
      rate(:,:,245) = 1.794e-10_r8
      rate(:,:,246) = 1.3e-10_r8
      rate(:,:,247) = 7.65e-11_r8
      rate(:,:,255) = 1.31e-10_r8
      rate(:,:,256) = 3.5e-11_r8
      rate(:,:,257) = 9e-12_r8
      rate(:,:,263) = 2.3e-12_r8
      rate(:,:,265) = 1.2e-11_r8
      rate(:,:,266) = 5.7e-11_r8
      rate(:,:,267) = 2.8e-11_r8
      rate(:,:,268) = 6.6e-11_r8
      rate(:,:,269) = 1.4e-11_r8
      rate(:,:,272) = 1.9e-12_r8
      rate(:,:,297) = 0.047_r8
      rate(:,:,298) = 7.7e-05_r8
      rate(:,:,299) = 0.171_r8
      rate(:,:,303) = 6e-11_r8
      rate(:,:,306) = 1e-12_r8
      rate(:,:,307) = 4e-10_r8
      rate(:,:,308) = 2e-10_r8
      rate(:,:,309) = 1e-10_r8
      rate(:,:,310) = 5e-16_r8
      rate(:,:,311) = 4.4e-10_r8
      rate(:,:,312) = 9e-10_r8
      rate(:,:,314) = 1.3e-10_r8
      rate(:,:,317) = 8e-10_r8
      rate(:,:,318) = 5e-12_r8
      rate(:,:,319) = 7e-10_r8
      rate(:,:,322) = 4.8e-10_r8
      rate(:,:,323) = 1e-10_r8
      rate(:,:,324) = 4e-10_r8
      rate(:,:,328) = 4.9e-12_r8
      rate(:,:,330) = 3.9e-11_r8
      rate(:,:,336) = 2.7e-9_r8
      rate(:,:,337) = 8.0e-10_r8
      rate(:,:,340) = 6.e-10_r8
      rate(:,:,341) = 6.e-10_r8
      rate(:,:,342) = 4.e-10_r8
      rate(:,:,343) = 1.e-11_r8
      rate(:,:,344) = 1.e-12_r8
      rate(:,:,345) = 5.e-12_r8
      rate(:,:,357) = 5.0e-12_r8
      rate(:,:,362) = 3e-11_r8
      rate(:,:,363) = 5.0e-11_r8
      rate(:,:,364) = 5.0e-11_r8
      rate(:,:,365) = 1.1e-9_r8
      rate(:,:,366) = 9.2e-10_r8
      rate(:,:,367) = 9.0e-10_r8
      rate(:,:,382) = 1.0e-12_r8
      rate(:,:,385) = 6.7e-12_r8
      rate(:,:,387) = 8.0e-14_r8
      rate(:,:,389) = 2.0e-10_r8
      rate(:,:,390) = 2.0e-12_r8
      rate(:,:,391) = 9.0e-10_r8
      rate(:,:,392) = 1.2e-9_r8
      rate(:,:,393) = 8.2e-10_r8
      rate(:,:,394) = 1.17e-9_r8
      rate(:,:,395) = 5.9e-10_r8
      rate(:,:,397) = 3.5e-12_r8
      rate(:,:,399) = 6.5e-10_r8
      rate(:,:,400) = 1.8e-10_r8
      rate(:,:,401) = 3.3e-10_r8
      rate(:,:,410) = 4.0e-10_r8
      rate(:,:,411) = 5.0e-13_r8
      rate(:,:,416) = 3.3e-10_r8
      rate(:,:,417) = 3.2e-10_r8
      rate(:,:,421) = 6.0e-10_r8
      rate(:,:,423) = 3.2e-10_r8
      rate(:,:,424) = 2.0e-10_r8
      rate(:,:,425) = 2.0e-10_r8
      rate(:,:,430) = 2.0e-10_r8
      rate(:,:,433) = 2.7e-12_r8
      rate(:,:,439) = 5.0e-12_r8
      rate(:,:,440) = 5.0e-12_r8
      rate(:,:,444) = 0.0e-11_r8
      rate(:,:,445) = 2.0e-11_r8
      rate(:,:,446) = 1.7e-11_r8
      rate(:,:,448) = 0.0e-12_r8
      rate(:,:,449) = 6.7e-12_r8
      rate(:,:,450) = 1.2e-11_r8
      rate(:,:,451) = 1.0e-10_r8
      rate(:,:,453) = 2.0e-10_r8
      rate(:,:,454) = 1.5e-10_r8
      rate(:,:,455) = 9.0e-8_r8
      rate(:,:,456) = 9.0e-8_r8
      rate(:,:,457) = 9.0e-8_r8
      rate(:,:,458) = 9.0e-8_r8
      rate(:,:,459) = 9.0e-8_r8
      rate(:,:,460) = 1.8e-9_r8
      rate(:,:,461) = 4.0e-9_r8
      rate(:,:,462) = 3.9e-10_r8
      rate(:,:,463) = 4.2e-11_r8
      rate(:,:,465) = 1.0e-10_r8
      rate(:,:,468) = 1.2e-10_r8
      rate(:,:,469) = 1.3e-9_r8
      rate(:,:,470) = 4.e-10_r8
      rate(:,:,476) = 3e-7_r8
      rate(:,:,478) = 3e-10_r8
      rate(:,:,490) = 2.7e-7_r8
      rate(:,:,491) = 9.4e-10_r8
      rate(:,:,492) = 3.2e-9_r8
      rate(:,:,521) = 1.0e-12_r8
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
      rate(:,:,218) = 5.5e-12_r8 * exp_fac(:,:)
      rate(:,:,253) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,131) = 1e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,132) = 1.4e-10_r8 * exp( -470._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1800._r8 * itemp(:,:) )
      rate(:,:,133) = 2.8e-12_r8 * exp_fac(:,:)
      rate(:,:,498) = 2.6e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,135) = 4.8e-11_r8 * exp_fac(:,:)
      rate(:,:,216) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,136) = 1.8e-11_r8 * exp( 180._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -940._r8 * itemp(:,:) )
      rate(:,:,137) = 1.7e-12_r8 * exp_fac(:,:)
      rate(:,:,334) = 5.e-10_r8 * exp_fac(:,:)
      rate(:,:,141) = 4.5e-13_r8 * exp( 610._r8 * itemp(:,:) )
      rate(:,:,144) = 2.1e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,145) = 2.9e-12_r8 * exp_fac(:,:)
      rate(:,:,146) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,147) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,148) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:,:) )
      rate(:,:,149) = 5.1e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,150) = 1.2e-13_r8 * exp_fac(:,:)
      rate(:,:,176) = 3e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 125._r8 * itemp(:,:) )
      rate(:,:,153) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,250) = 5.5e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,157) = 3.44e-12_r8 * exp_fac(:,:)
      rate(:,:,209) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,212) = 8.8e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,158) = 3e-12_r8 * exp_fac(:,:)
      rate(:,:,217) = 5.8e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,160) = 7.26e-11_r8 * exp_fac(:,:)
      rate(:,:,161) = 4.64e-11_r8 * exp_fac(:,:)
      rate(:,:,168) = 8.1e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,169) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:,:) )
      rate(:,:,170) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,171) = 1.1e-11_r8 * exp( -980._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,172) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,186) = 7.4e-12_r8 * exp_fac(:,:)
      rate(:,:,173) = 3.6e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,174) = 2.3e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,175) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,177) = 1e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,178) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,179) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,180) = 6.4e-12_r8 * exp_fac(:,:)
      rate(:,:,210) = 4.1e-13_r8 * exp_fac(:,:)
      rate(:,:,181) = 6.5e-12_r8 * exp( 135._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,183) = 3.6e-12_r8 * exp_fac(:,:)
      rate(:,:,232) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,184) = 1.2e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,185) = 2.8e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,187) = 6e-13_r8 * exp_fac(:,:)
      rate(:,:,207) = 1.5e-12_r8 * exp_fac(:,:)
      rate(:,:,215) = 1.9e-11_r8 * exp_fac(:,:)
      rate(:,:,188) = 1e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,189) = 1.8e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,190) = 3.4e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,192) = 3e-12_r8 * exp_fac(:,:)
      rate(:,:,226) = 1.4e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,204) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,231) = 6.3e-12_r8 * exp_fac(:,:)
      rate(:,:,205) = 4.8e-12_r8 * exp( -310._r8 * itemp(:,:) )
      rate(:,:,206) = 1.6e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,208) = 9.5e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,211) = 4.5e-12_r8 * exp( 460._r8 * itemp(:,:) )
      rate(:,:,214) = 1.9e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,219) = 1.2e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,225) = 1.6e-10_r8 * exp( -260._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,227) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,229) = 2.14e-11_r8 * exp_fac(:,:)
      rate(:,:,230) = 1.9e-10_r8 * exp_fac(:,:)
      rate(:,:,243) = 2.57e-10_r8 * exp_fac(:,:)
      rate(:,:,244) = 1.8e-10_r8 * exp_fac(:,:)
      rate(:,:,245) = 1.794e-10_r8 * exp_fac(:,:)
      rate(:,:,246) = 1.3e-10_r8 * exp_fac(:,:)
      rate(:,:,247) = 7.65e-11_r8 * exp_fac(:,:)
      rate(:,:,255) = 1.31e-10_r8 * exp_fac(:,:)
      rate(:,:,256) = 3.5e-11_r8 * exp_fac(:,:)
      rate(:,:,257) = 9e-12_r8 * exp_fac(:,:)
      rate(:,:,263) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,265) = 1.2e-11_r8 * exp_fac(:,:)
      rate(:,:,266) = 5.7e-11_r8 * exp_fac(:,:)
      rate(:,:,267) = 2.8e-11_r8 * exp_fac(:,:)
      rate(:,:,268) = 6.6e-11_r8 * exp_fac(:,:)
      rate(:,:,269) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,272) = 1.9e-12_r8 * exp_fac(:,:)
      rate(:,:,297) = 0.047_r8 * exp_fac(:,:)
      rate(:,:,298) = 7.7e-05_r8 * exp_fac(:,:)
      rate(:,:,299) = 0.171_r8 * exp_fac(:,:)
      rate(:,:,303) = 6e-11_r8 * exp_fac(:,:)
      rate(:,:,306) = 1e-12_r8 * exp_fac(:,:)
      rate(:,:,307) = 4e-10_r8 * exp_fac(:,:)
      rate(:,:,308) = 2e-10_r8 * exp_fac(:,:)
      rate(:,:,309) = 1e-10_r8 * exp_fac(:,:)
      rate(:,:,310) = 5e-16_r8 * exp_fac(:,:)
      rate(:,:,311) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,312) = 9e-10_r8 * exp_fac(:,:)
      rate(:,:,314) = 1.3e-10_r8 * exp_fac(:,:)
      rate(:,:,317) = 8e-10_r8 * exp_fac(:,:)
      rate(:,:,318) = 5e-12_r8 * exp_fac(:,:)
      rate(:,:,319) = 7e-10_r8 * exp_fac(:,:)
      rate(:,:,322) = 4.8e-10_r8 * exp_fac(:,:)
      rate(:,:,323) = 1e-10_r8 * exp_fac(:,:)
      rate(:,:,324) = 4e-10_r8 * exp_fac(:,:)
      rate(:,:,328) = 4.9e-12_r8 * exp_fac(:,:)
      rate(:,:,330) = 3.9e-11_r8 * exp_fac(:,:)
      rate(:,:,336) = 2.7e-9_r8 * exp_fac(:,:)
      rate(:,:,337) = 8.0e-10_r8 * exp_fac(:,:)
      rate(:,:,340) = 6.e-10_r8 * exp_fac(:,:)
      rate(:,:,341) = 6.e-10_r8 * exp_fac(:,:)
      rate(:,:,342) = 4.e-10_r8 * exp_fac(:,:)
      rate(:,:,343) = 1.e-11_r8 * exp_fac(:,:)
      rate(:,:,344) = 1.e-12_r8 * exp_fac(:,:)
      rate(:,:,345) = 5.e-12_r8 * exp_fac(:,:)
      rate(:,:,357) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,362) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,363) = 5.0e-11_r8 * exp_fac(:,:)
      rate(:,:,364) = 5.0e-11_r8 * exp_fac(:,:)
      rate(:,:,365) = 1.1e-9_r8 * exp_fac(:,:)
      rate(:,:,366) = 9.2e-10_r8 * exp_fac(:,:)
      rate(:,:,367) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,382) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,385) = 6.7e-12_r8 * exp_fac(:,:)
      rate(:,:,387) = 8.0e-14_r8 * exp_fac(:,:)
      rate(:,:,389) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,390) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,391) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,392) = 1.2e-9_r8 * exp_fac(:,:)
      rate(:,:,393) = 8.2e-10_r8 * exp_fac(:,:)
      rate(:,:,394) = 1.17e-9_r8 * exp_fac(:,:)
      rate(:,:,395) = 5.9e-10_r8 * exp_fac(:,:)
      rate(:,:,397) = 3.5e-12_r8 * exp_fac(:,:)
      rate(:,:,399) = 6.5e-10_r8 * exp_fac(:,:)
      rate(:,:,400) = 1.8e-10_r8 * exp_fac(:,:)
      rate(:,:,401) = 3.3e-10_r8 * exp_fac(:,:)
      rate(:,:,410) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,411) = 5.0e-13_r8 * exp_fac(:,:)
      rate(:,:,416) = 3.3e-10_r8 * exp_fac(:,:)
      rate(:,:,417) = 3.2e-10_r8 * exp_fac(:,:)
      rate(:,:,421) = 6.0e-10_r8 * exp_fac(:,:)
      rate(:,:,423) = 3.2e-10_r8 * exp_fac(:,:)
      rate(:,:,424) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,425) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,430) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,433) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,439) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,440) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,444) = 0.0e-11_r8 * exp_fac(:,:)
      rate(:,:,445) = 2.0e-11_r8 * exp_fac(:,:)
      rate(:,:,446) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,448) = 0.0e-12_r8 * exp_fac(:,:)
      rate(:,:,449) = 6.7e-12_r8 * exp_fac(:,:)
      rate(:,:,450) = 1.2e-11_r8 * exp_fac(:,:)
      rate(:,:,451) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,453) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,454) = 1.5e-10_r8 * exp_fac(:,:)
      rate(:,:,455) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,456) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,457) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,458) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,459) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,460) = 1.8e-9_r8 * exp_fac(:,:)
      rate(:,:,461) = 4.0e-9_r8 * exp_fac(:,:)
      rate(:,:,462) = 3.9e-10_r8 * exp_fac(:,:)
      rate(:,:,463) = 4.2e-11_r8 * exp_fac(:,:)
      rate(:,:,465) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,468) = 1.2e-10_r8 * exp_fac(:,:)
      rate(:,:,469) = 1.3e-9_r8 * exp_fac(:,:)
      rate(:,:,470) = 4.e-10_r8 * exp_fac(:,:)
      rate(:,:,476) = 3e-7_r8 * exp_fac(:,:)
      rate(:,:,478) = 3e-10_r8 * exp_fac(:,:)
      rate(:,:,490) = 2.7e-7_r8 * exp_fac(:,:)
      rate(:,:,491) = 9.4e-10_r8 * exp_fac(:,:)
      rate(:,:,492) = 3.2e-9_r8 * exp_fac(:,:)
      rate(:,:,521) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,228) = 6e-12_r8 * exp( 400._r8 * itemp(:,:) )
      rate(:,:,233) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:,:) )
      rate(:,:,234) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:,:) )
      rate(:,:,235) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      rate(:,:,236) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1200._r8 * itemp(:,:) )
      rate(:,:,237) = 1.96e-12_r8 * exp_fac(:,:)
      rate(:,:,376) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,238) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,239) = 9e-13_r8 * exp( -360._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,240) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,249) = 3.4e-11_r8 * exp_fac(:,:)
      rate(:,:,241) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,242) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:,:) )
      rate(:,:,248) = 6e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,251) = 4.1e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,252) = 2.8e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,254) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,259) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,260) = 1.1e-11_r8 * exp( -280._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2200._r8 * itemp(:,:) )
      rate(:,:,261) = 2.1e-11_r8 * exp_fac(:,:)
      rate(:,:,508) = 1.4e-9_r8 * exp_fac(:,:)
      rate(:,:,262) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:,:) )
      rate(:,:,270) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:,:) )
      rate(:,:,271) = 3.4e-12_r8 * exp( -1100._r8 * itemp(:,:) )
      rate(:,:,273) = 2.6e-11_r8 * exp( 330._r8 * itemp(:,:) )
      rate(:,:,325) = 1.1e-9_r8 * exp( -116._r8 * itemp(:,:) )
      rate(:,:,327) = 3.2e-10_r8 * exp( -550._r8 * itemp(:,:) )
      rate(:,:,329) = 5.06e-10_r8 * exp( -240._r8 * itemp(:,:) )
      rate(:,:,351) = 2.94e-10_r8 * exp( -174._r8 * itemp(:,:) )
      rate(:,:,352) = 4.6e-10_r8 * exp( -350._r8 * itemp(:,:) )
      rate(:,:,353) = 3.0e-10_r8 * exp( -177._r8 * itemp(:,:) )
      rate(:,:,354) = 1.4e-10_r8 * exp( -580._r8 * itemp(:,:) )
      rate(:,:,355) = 4.4e-10_r8 * exp( -170._r8 * itemp(:,:) )
      rate(:,:,356) = 2.3e-10_r8 * exp( -2310._r8 * itemp(:,:) )
      rate(:,:,358) = 3.3e-10_r8 * exp( -302._r8 * itemp(:,:) )
      rate(:,:,359) = 3.0e-10_r8 * exp( -796._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -600._r8 * itemp(:,:) )
      rate(:,:,360) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,388) = 1.0e-11_r8 * exp_fac(:,:)
      rate(:,:,361) = 7.6e-10_r8 * exp( -240.5_r8 * itemp(:,:) )
      rate(:,:,377) = 2.28e-10_r8 * exp( -139._r8 * itemp(:,:) )
      rate(:,:,379) = 2.19e-10_r8 * exp( -548._r8 * itemp(:,:) )
      rate(:,:,432) = 8.23e-10_r8 * exp( -192._r8 * itemp(:,:) )
      rate(:,:,434) = 1.1e-9_r8 * exp( -421._r8 * itemp(:,:) )
      rate(:,:,435) = 5.7e-10_r8 * exp( -267._r8 * itemp(:,:) )
      rate(:,:,437) = 4.4e-11_r8 * exp( -202._r8 * itemp(:,:) )
      rate(:,:,443) = 7.0e-10_r8 * exp( -4017._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -120._r8 * itemp(:,:) )
      rate(:,:,479) = 1.15e-9_r8 * exp_fac(:,:)
      rate(:,:,480) = 2.00e-10_r8 * exp_fac(:,:)
      rate(:,:,483) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,484) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,485) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,486) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,487) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,482) = 6.9e-10_r8 * exp( -385._r8 * itemp(:,:) )
      rate(:,:,489) = 4.5e-11_r8 * exp( -3590._r8 * itemp(:,:) )
      rate(:,:,494) = 2.8e-8_r8 * exp( -1680._r8 * itemp(:,:) )
      rate(:,:,496) = 1.5e-9_r8 * exp( -820._r8 * itemp(:,:) )
      rate(:,:,504) = 2.8e-10_r8 * exp( -2220._r8 * itemp(:,:) )
      rate(:,:,506) = 5.0e-10_r8 * exp( -752._r8 * itemp(:,:) )
      rate(:,:,511) = 1.0e-9_r8 * exp( -873._r8 * itemp(:,:) )

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
      call jpl( rate(1,1,182), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.9e-32_r8 * itemp(:,:)**3.6_r8
      kinf(:,:) = 3.7e-12_r8 * itemp(:,:)**1.6_r8
      call jpl( rate(1,1,202), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.2e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,213), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.9e-31_r8 * itemp(:,:)**4.1_r8
      kinf(:,:) = 1.7e-12_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,264), m, 0.6_r8, ko, kinf, n )

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
      rate(:,:kbot,142) = 7e-13_r8
      rate(:,:kbot,143) = 5e-12_r8
      rate(:,:kbot,297) = 0.047_r8
      rate(:,:kbot,298) = 7.7e-05_r8
      rate(:,:kbot,299) = 0.171_r8
      rate(:,:kbot,303) = 6e-11_r8
      rate(:,:kbot,306) = 1e-12_r8
      rate(:,:kbot,307) = 4e-10_r8
      rate(:,:kbot,308) = 2e-10_r8
      rate(:,:kbot,309) = 1e-10_r8
      rate(:,:kbot,311) = 4.4e-10_r8
      rate(:,:kbot,314) = 1.3e-10_r8
      rate(:,:kbot,317) = 8e-10_r8
      rate(:,:kbot,318) = 5e-12_r8
      rate(:,:kbot,319) = 7e-10_r8
      rate(:,:kbot,322) = 4.8e-10_r8
      rate(:,:kbot,323) = 1e-10_r8
      rate(:,:kbot,324) = 4e-10_r8
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
      rate(:,:kbot,144) = 2.1e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,148) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:,:) )
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
