
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      real(r8)  :: itemp(ncol*pver)
      real(r8)  :: exp_fac(ncol*pver)
      real(r8)  :: ko(ncol*pver)
      real(r8)  :: kinf(ncol*pver)

      rate(:,124) = 1.2e-10_r8
      rate(:,128) = 1.2e-10_r8
      rate(:,134) = 6.9e-12_r8
      rate(:,135) = 7.2e-11_r8
      rate(:,136) = 1.6e-12_r8
      rate(:,142) = 1.8e-12_r8
      rate(:,146) = 1.8e-12_r8
      rate(:,158) = 3.5e-12_r8
      rate(:,160) = 1e-11_r8
      rate(:,161) = 2.2e-11_r8
      rate(:,162) = 5e-11_r8
      rate(:,197) = 1.7e-13_r8
      rate(:,199) = 2.607e-10_r8
      rate(:,200) = 9.75e-11_r8
      rate(:,201) = 2.07e-10_r8
      rate(:,202) = 2.088e-10_r8
      rate(:,203) = 1.17e-10_r8
      rate(:,204) = 4.644e-11_r8
      rate(:,205) = 1.204e-10_r8
      rate(:,206) = 9.9e-11_r8
      rate(:,207) = 3.3e-12_r8
      rate(:,226) = 4.5e-11_r8
      rate(:,227) = 4.62e-10_r8
      rate(:,228) = 1.2e-10_r8
      rate(:,229) = 9e-11_r8
      rate(:,230) = 3e-11_r8
      rate(:,235) = 2.14e-11_r8
      rate(:,236) = 1.9e-10_r8
      rate(:,249) = 2.57e-10_r8
      rate(:,250) = 1.8e-10_r8
      rate(:,251) = 1.794e-10_r8
      rate(:,252) = 1.3e-10_r8
      rate(:,253) = 7.65e-11_r8
      rate(:,267) = 4e-13_r8
      rate(:,271) = 1.31e-10_r8
      rate(:,272) = 3.5e-11_r8
      rate(:,273) = 9e-12_r8
      rate(:,306) = 6.8e-14_r8
      rate(:,307) = 2e-13_r8
      rate(:,321) = 7e-13_r8
      rate(:,322) = 1e-12_r8
      rate(:,326) = 1e-14_r8
      rate(:,327) = 1e-11_r8
      rate(:,328) = 1.15e-11_r8
      rate(:,329) = 4e-14_r8
      rate(:,342) = 3e-12_r8
      rate(:,343) = 6.7e-13_r8
      rate(:,353) = 3.5e-13_r8
      rate(:,354) = 5.4e-11_r8
      rate(:,357) = 2e-12_r8
      rate(:,358) = 1.4e-11_r8
      rate(:,361) = 2.4e-12_r8
      rate(:,372) = 5e-12_r8
      rate(:,382) = 1.6e-12_r8
      rate(:,384) = 6.7e-12_r8
      rate(:,387) = 3.5e-12_r8
      rate(:,390) = 1.3e-11_r8
      rate(:,391) = 1.4e-11_r8
      rate(:,395) = 2.4e-12_r8
      rate(:,396) = 1.4e-11_r8
      rate(:,401) = 2.4e-12_r8
      rate(:,402) = 4e-11_r8
      rate(:,403) = 4e-11_r8
      rate(:,405) = 1.4e-11_r8
      rate(:,409) = 2.4e-12_r8
      rate(:,410) = 4e-11_r8
      rate(:,414) = 7e-11_r8
      rate(:,415) = 1e-10_r8
      rate(:,420) = 2.4e-12_r8
      rate(:,435) = 4.7e-11_r8
      rate(:,448) = 2.1e-12_r8
      rate(:,449) = 2.8e-13_r8
      rate(:,457) = 1.7e-11_r8
      rate(:,463) = 8.4e-11_r8
      rate(:,465) = 1.9e-11_r8
      rate(:,466) = 1.2e-14_r8
      rate(:,467) = 2e-10_r8
      rate(:,474) = 2.4e-12_r8
      rate(:,475) = 2e-11_r8
      rate(:,479) = 2.3e-11_r8
      rate(:,480) = 2e-11_r8
      rate(:,484) = 3.3e-11_r8
      rate(:,485) = 1e-12_r8
      rate(:,486) = 5.7e-11_r8
      rate(:,487) = 3.4e-11_r8
      rate(:,492) = 2.3e-12_r8
      rate(:,493) = 1.2e-11_r8
      rate(:,494) = 5.7e-11_r8
      rate(:,495) = 2.8e-11_r8
      rate(:,496) = 6.6e-11_r8
      rate(:,497) = 1.4e-11_r8
      rate(:,500) = 1.9e-12_r8
      rate(:,514) = 6.34e-08_r8
      rate(:,520) = 1.9e-11_r8
      rate(:,523) = 1.2e-14_r8
      rate(:,524) = 2e-10_r8
      rate(:,535) = 1.34e-11_r8
      rate(:,541) = 1.34e-11_r8
      rate(:,545) = 1.7e-11_r8
      rate(:,565) = 1.29e-07_r8
      rate(:,566) = 2.31e-07_r8
      rate(:,567) = 2.31e-06_r8
      rate(:,568) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,125) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,126) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,127) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,129) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,132) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,133) = 1.4e-12_r8 * exp_fac(:)
      rate(:,411) = 1.05e-14_r8 * exp_fac(:)
      rate(:,531) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,138) = 3e-11_r8 * exp_fac(:)
      rate(:,224) = 5.5e-12_r8 * exp_fac(:)
      rate(:,263) = 3.8e-12_r8 * exp_fac(:)
      rate(:,311) = 3.8e-12_r8 * exp_fac(:)
      rate(:,338) = 3.8e-12_r8 * exp_fac(:)
      rate(:,346) = 3.8e-12_r8 * exp_fac(:)
      rate(:,350) = 3.8e-12_r8 * exp_fac(:)
      rate(:,366) = 2.3e-11_r8 * exp_fac(:)
      rate(:,376) = 3.8e-12_r8 * exp_fac(:)
      rate(:,386) = 3.8e-12_r8 * exp_fac(:)
      rate(:,413) = 1.52e-11_r8 * exp_fac(:)
      rate(:,421) = 1.52e-12_r8 * exp_fac(:)
      rate(:,427) = 3.8e-12_r8 * exp_fac(:)
      rate(:,430) = 3.8e-12_r8 * exp_fac(:)
      rate(:,434) = 3.8e-12_r8 * exp_fac(:)
      rate(:,450) = 3.8e-12_r8 * exp_fac(:)
      rate(:,454) = 3.8e-12_r8 * exp_fac(:)
      rate(:,460) = 3.8e-12_r8 * exp_fac(:)
      rate(:,464) = 3.8e-12_r8 * exp_fac(:)
      rate(:,139) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,140) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,141) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,143) = 4.8e-11_r8 * exp_fac(:)
      rate(:,222) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,144) = 1.8e-11_r8 * exp_fac(:)
      rate(:,324) = 4.2e-12_r8 * exp_fac(:)
      rate(:,337) = 4.2e-12_r8 * exp_fac(:)
      rate(:,345) = 4.2e-12_r8 * exp_fac(:)
      rate(:,374) = 4.2e-12_r8 * exp_fac(:)
      rate(:,394) = 4.4e-12_r8 * exp_fac(:)
      rate(:,400) = 4.4e-12_r8 * exp_fac(:)
      rate(:,473) = 4.2e-12_r8 * exp_fac(:)
      rate(:,478) = 4.2e-12_r8 * exp_fac(:)
      rate(:,483) = 4.2e-12_r8 * exp_fac(:)
      rate(:,145) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,149) = 1.3e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,150) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,151) = 2.9e-12_r8 * exp_fac(:)
      rate(:,152) = 1.45e-12_r8 * exp_fac(:)
      rate(:,153) = 1.45e-12_r8 * exp_fac(:)
      rate(:,154) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,155) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,156) = 1.2e-13_r8 * exp_fac(:)
      rate(:,182) = 3e-11_r8 * exp_fac(:)
      rate(:,159) = 1.5e-11_r8 * exp( 170._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,163) = 3.3e-12_r8 * exp_fac(:)
      rate(:,178) = 1.4e-11_r8 * exp_fac(:)
      rate(:,192) = 7.4e-12_r8 * exp_fac(:)
      rate(:,320) = 8.1e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,164) = 3e-12_r8 * exp_fac(:)
      rate(:,223) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,166) = 7.26e-11_r8 * exp_fac(:)
      rate(:,167) = 4.64e-11_r8 * exp_fac(:)
      rate(:,174) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,175) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,176) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,177) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,179) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,180) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,181) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,183) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,184) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,185) = 2.6e-12_r8 * exp_fac(:)
      rate(:,186) = 6.4e-12_r8 * exp_fac(:)
      rate(:,216) = 4.1e-13_r8 * exp_fac(:)
      rate(:,423) = 7.5e-12_r8 * exp_fac(:)
      rate(:,437) = 7.5e-12_r8 * exp_fac(:)
      rate(:,440) = 7.5e-12_r8 * exp_fac(:)
      rate(:,443) = 7.5e-12_r8 * exp_fac(:)
      rate(:,187) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,189) = 3.6e-12_r8 * exp_fac(:)
      rate(:,238) = 2e-12_r8 * exp_fac(:)
      rate(:,190) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,191) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,193) = 6e-13_r8 * exp_fac(:)
      rate(:,213) = 1.5e-12_r8 * exp_fac(:)
      rate(:,221) = 1.9e-11_r8 * exp_fac(:)
      rate(:,194) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,195) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,196) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,198) = 3e-12_r8 * exp_fac(:)
      rate(:,232) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,210) = 1.7e-11_r8 * exp_fac(:)
      rate(:,237) = 6.3e-12_r8 * exp_fac(:)
      rate(:,211) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,212) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,214) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,215) = 2.3e-12_r8 * exp_fac(:)
      rate(:,218) = 8.8e-12_r8 * exp_fac(:)
      rate(:,217) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,220) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,225) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,231) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,233) = 1.4e-11_r8 * exp_fac(:)
      rate(:,235) = 2.14e-11_r8 * exp_fac(:)
      rate(:,236) = 1.9e-10_r8 * exp_fac(:)
      rate(:,249) = 2.57e-10_r8 * exp_fac(:)
      rate(:,250) = 1.8e-10_r8 * exp_fac(:)
      rate(:,251) = 1.794e-10_r8 * exp_fac(:)
      rate(:,252) = 1.3e-10_r8 * exp_fac(:)
      rate(:,253) = 7.65e-11_r8 * exp_fac(:)
      rate(:,267) = 4e-13_r8 * exp_fac(:)
      rate(:,271) = 1.31e-10_r8 * exp_fac(:)
      rate(:,272) = 3.5e-11_r8 * exp_fac(:)
      rate(:,273) = 9e-12_r8 * exp_fac(:)
      rate(:,306) = 6.8e-14_r8 * exp_fac(:)
      rate(:,307) = 2e-13_r8 * exp_fac(:)
      rate(:,321) = 7e-13_r8 * exp_fac(:)
      rate(:,322) = 1e-12_r8 * exp_fac(:)
      rate(:,326) = 1e-14_r8 * exp_fac(:)
      rate(:,327) = 1e-11_r8 * exp_fac(:)
      rate(:,328) = 1.15e-11_r8 * exp_fac(:)
      rate(:,329) = 4e-14_r8 * exp_fac(:)
      rate(:,342) = 3e-12_r8 * exp_fac(:)
      rate(:,343) = 6.7e-13_r8 * exp_fac(:)
      rate(:,353) = 3.5e-13_r8 * exp_fac(:)
      rate(:,354) = 5.4e-11_r8 * exp_fac(:)
      rate(:,357) = 2e-12_r8 * exp_fac(:)
      rate(:,358) = 1.4e-11_r8 * exp_fac(:)
      rate(:,361) = 2.4e-12_r8 * exp_fac(:)
      rate(:,372) = 5e-12_r8 * exp_fac(:)
      rate(:,382) = 1.6e-12_r8 * exp_fac(:)
      rate(:,384) = 6.7e-12_r8 * exp_fac(:)
      rate(:,387) = 3.5e-12_r8 * exp_fac(:)
      rate(:,390) = 1.3e-11_r8 * exp_fac(:)
      rate(:,391) = 1.4e-11_r8 * exp_fac(:)
      rate(:,395) = 2.4e-12_r8 * exp_fac(:)
      rate(:,396) = 1.4e-11_r8 * exp_fac(:)
      rate(:,401) = 2.4e-12_r8 * exp_fac(:)
      rate(:,402) = 4e-11_r8 * exp_fac(:)
      rate(:,403) = 4e-11_r8 * exp_fac(:)
      rate(:,405) = 1.4e-11_r8 * exp_fac(:)
      rate(:,409) = 2.4e-12_r8 * exp_fac(:)
      rate(:,410) = 4e-11_r8 * exp_fac(:)
      rate(:,414) = 7e-11_r8 * exp_fac(:)
      rate(:,415) = 1e-10_r8 * exp_fac(:)
      rate(:,420) = 2.4e-12_r8 * exp_fac(:)
      rate(:,435) = 4.7e-11_r8 * exp_fac(:)
      rate(:,448) = 2.1e-12_r8 * exp_fac(:)
      rate(:,449) = 2.8e-13_r8 * exp_fac(:)
      rate(:,457) = 1.7e-11_r8 * exp_fac(:)
      rate(:,463) = 8.4e-11_r8 * exp_fac(:)
      rate(:,465) = 1.9e-11_r8 * exp_fac(:)
      rate(:,466) = 1.2e-14_r8 * exp_fac(:)
      rate(:,467) = 2e-10_r8 * exp_fac(:)
      rate(:,474) = 2.4e-12_r8 * exp_fac(:)
      rate(:,475) = 2e-11_r8 * exp_fac(:)
      rate(:,479) = 2.3e-11_r8 * exp_fac(:)
      rate(:,480) = 2e-11_r8 * exp_fac(:)
      rate(:,484) = 3.3e-11_r8 * exp_fac(:)
      rate(:,485) = 1e-12_r8 * exp_fac(:)
      rate(:,486) = 5.7e-11_r8 * exp_fac(:)
      rate(:,487) = 3.4e-11_r8 * exp_fac(:)
      rate(:,492) = 2.3e-12_r8 * exp_fac(:)
      rate(:,493) = 1.2e-11_r8 * exp_fac(:)
      rate(:,494) = 5.7e-11_r8 * exp_fac(:)
      rate(:,495) = 2.8e-11_r8 * exp_fac(:)
      rate(:,496) = 6.6e-11_r8 * exp_fac(:)
      rate(:,497) = 1.4e-11_r8 * exp_fac(:)
      rate(:,500) = 1.9e-12_r8 * exp_fac(:)
      rate(:,514) = 6.34e-08_r8 * exp_fac(:)
      rate(:,520) = 1.9e-11_r8 * exp_fac(:)
      rate(:,523) = 1.2e-14_r8 * exp_fac(:)
      rate(:,524) = 2e-10_r8 * exp_fac(:)
      rate(:,535) = 1.34e-11_r8 * exp_fac(:)
      rate(:,541) = 1.34e-11_r8 * exp_fac(:)
      rate(:,545) = 1.7e-11_r8 * exp_fac(:)
      rate(:,565) = 1.29e-07_r8 * exp_fac(:)
      rate(:,566) = 2.31e-07_r8 * exp_fac(:)
      rate(:,567) = 2.31e-06_r8 * exp_fac(:)
      rate(:,568) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,234) = 6e-12_r8 * exp_fac(:)
      rate(:,359) = 5e-13_r8 * exp_fac(:)
      rate(:,392) = 5e-13_r8 * exp_fac(:)
      rate(:,397) = 5e-13_r8 * exp_fac(:)
      rate(:,406) = 5e-13_r8 * exp_fac(:)
      rate(:,417) = 5e-13_r8 * exp_fac(:)
      rate(:,239) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,240) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,241) = 1.64e-12_r8 * exp_fac(:)
      rate(:,378) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,242) = 2.03e-11_r8 * exp_fac(:)
      rate(:,499) = 3.4e-12_r8 * exp_fac(:)
      rate(:,243) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,244) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,245) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,246) = 1.25e-12_r8 * exp_fac(:)
      rate(:,256) = 3.4e-11_r8 * exp_fac(:)
      rate(:,247) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,248) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,254) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,255) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,257) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,258) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,259) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,260) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,261) = 2.8e-12_r8 * exp_fac(:)
      rate(:,349) = 2.9e-12_r8 * exp_fac(:)
      rate(:,262) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,264) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,268) = 7.5e-13_r8 * exp_fac(:)
      rate(:,308) = 7.5e-13_r8 * exp_fac(:)
      rate(:,323) = 7.5e-13_r8 * exp_fac(:)
      rate(:,336) = 7.5e-13_r8 * exp_fac(:)
      rate(:,344) = 7.5e-13_r8 * exp_fac(:)
      rate(:,348) = 8.6e-13_r8 * exp_fac(:)
      rate(:,360) = 8e-13_r8 * exp_fac(:)
      rate(:,373) = 7.5e-13_r8 * exp_fac(:)
      rate(:,383) = 7.5e-13_r8 * exp_fac(:)
      rate(:,393) = 8e-13_r8 * exp_fac(:)
      rate(:,398) = 8e-13_r8 * exp_fac(:)
      rate(:,407) = 8e-13_r8 * exp_fac(:)
      rate(:,418) = 8e-13_r8 * exp_fac(:)
      rate(:,425) = 7.5e-13_r8 * exp_fac(:)
      rate(:,429) = 7.5e-13_r8 * exp_fac(:)
      rate(:,432) = 7.5e-13_r8 * exp_fac(:)
      rate(:,445) = 7.5e-13_r8 * exp_fac(:)
      rate(:,452) = 7.5e-13_r8 * exp_fac(:)
      rate(:,458) = 7.5e-13_r8 * exp_fac(:)
      rate(:,461) = 7.5e-13_r8 * exp_fac(:)
      rate(:,472) = 7.5e-13_r8 * exp_fac(:)
      rate(:,477) = 7.5e-13_r8 * exp_fac(:)
      rate(:,482) = 7.5e-13_r8 * exp_fac(:)
      rate(:,526) = 7.5e-13_r8 * exp_fac(:)
      rate(:,533) = 7.5e-13_r8 * exp_fac(:)
      rate(:,543) = 7.5e-13_r8 * exp_fac(:)
      rate(:,546) = 7.5e-13_r8 * exp_fac(:)
      rate(:,269) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,270) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,274) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,305) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,309) = 2.6e-12_r8 * exp_fac(:)
      rate(:,426) = 2.6e-12_r8 * exp_fac(:)
      rate(:,431) = 2.6e-12_r8 * exp_fac(:)
      rate(:,433) = 2.6e-12_r8 * exp_fac(:)
      rate(:,446) = 2.6e-12_r8 * exp_fac(:)
      rate(:,453) = 2.6e-12_r8 * exp_fac(:)
      rate(:,459) = 2.6e-12_r8 * exp_fac(:)
      rate(:,462) = 2.6e-12_r8 * exp_fac(:)
      rate(:,527) = 2.6e-12_r8 * exp_fac(:)
      rate(:,534) = 2.6e-12_r8 * exp_fac(:)
      rate(:,544) = 2.6e-12_r8 * exp_fac(:)
      rate(:,547) = 2.6e-12_r8 * exp_fac(:)
      rate(:,310) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,312) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,313) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,314) = 1.4e-12_r8 * exp_fac(:)
      rate(:,334) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,315) = 4.63e-12_r8 * exp_fac(:)
      rate(:,530) = 2.7e-12_r8 * exp_fac(:)
      rate(:,316) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,317) = 2.9e-12_r8 * exp_fac(:)
      rate(:,318) = 2e-12_r8 * exp_fac(:)
      rate(:,347) = 7.1e-13_r8 * exp_fac(:)
      rate(:,368) = 2e-12_r8 * exp_fac(:)
      rate(:,471) = 2e-12_r8 * exp_fac(:)
      rate(:,476) = 2e-12_r8 * exp_fac(:)
      rate(:,481) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,319) = 4.3e-13_r8 * exp_fac(:)
      rate(:,369) = 4.3e-13_r8 * exp_fac(:)
      rate(:,422) = 4.3e-13_r8 * exp_fac(:)
      rate(:,436) = 4.3e-13_r8 * exp_fac(:)
      rate(:,439) = 4.3e-13_r8 * exp_fac(:)
      rate(:,442) = 4.3e-13_r8 * exp_fac(:)
      rate(:,325) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,333) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,335) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,339) = 8.7e-12_r8 * exp( -615._r8 * itemp(:) )
      rate(:,340) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,341) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,355) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,356) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,362) = 2.7e-12_r8 * exp_fac(:)
      rate(:,363) = 1.3e-13_r8 * exp_fac(:)
      rate(:,365) = 9.6e-12_r8 * exp_fac(:)
      rate(:,371) = 5.3e-12_r8 * exp_fac(:)
      rate(:,408) = 2.7e-12_r8 * exp_fac(:)
      rate(:,419) = 2.7e-12_r8 * exp_fac(:)
      rate(:,522) = 2.7e-12_r8 * exp_fac(:)
      rate(:,538) = 2.7e-12_r8 * exp_fac(:)
      rate(:,364) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,367) = 4.6e-12_r8 * exp_fac(:)
      rate(:,370) = 2.3e-12_r8 * exp_fac(:)
      rate(:,375) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,379) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,385) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,388) = 1.86e-11_r8 * exp_fac(:)
      rate(:,389) = 1.86e-11_r8 * exp_fac(:)
      rate(:,399) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,404) = 3.03e-12_r8 * exp_fac(:)
      rate(:,528) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,412) = 2.54e-11_r8 * exp_fac(:)
      rate(:,532) = 2.54e-11_r8 * exp_fac(:)
      rate(:,416) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,424) = 2.3e-12_r8 * exp_fac(:)
      rate(:,525) = 2.3e-12_r8 * exp_fac(:)
      rate(:,428) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,447) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,455) = 1.7e-12_r8 * exp_fac(:)
      rate(:,542) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,468) = 1.2e-12_r8 * exp_fac(:)
      rate(:,536) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,469) = 6.3e-16_r8 * exp_fac(:)
      rate(:,539) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,470) = 1.2e-11_r8 * exp_fac(:)
      rate(:,540) = 1.2e-11_r8 * exp_fac(:)
      rate(:,488) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,489) = 9.6e-12_r8 * exp( -234._r8 * itemp(:) )
      rate(:,490) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,491) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,498) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,501) = 2.7e-11_r8 * exp( 335._r8 * itemp(:) )
      rate(:,505) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,521) = 2.75e-13_r8 * exp_fac(:)
      rate(:,529) = 2.12e-13_r8 * exp_fac(:)
      rate(:,537) = 2.6e-13_r8 * exp_fac(:)

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,137), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,147), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,157), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,165), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,168), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,169), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,170), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,188), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,208), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,219), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,265), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.28e-33_r8
      kinf(:) = 9.3e-15_r8 * itemp(:)**(-4.42_r8)
      call jpl( rate(:,266), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,277), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,279), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,281), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,283), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,285), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,287), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,289), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,291), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,293), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,295), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,297), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,299), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,301), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,302), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,303), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,304), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,330), m, 0.48_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,331), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,351), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,377), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,380), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,438), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,441), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,444), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,451), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      integer   ::  k
      real(r8)  :: itemp(ncol*kbot)
      real(r8)  :: exp_fac(ncol*kbot)
      real(r8)  :: ko(ncol*kbot)
      real(r8)  :: kinf(ncol*kbot)
      real(r8)  :: wrk(ncol*kbot)
 
      n = ncol*kbot

      rate(:n,134) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,126) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,129) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,138) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,139) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,140) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,143) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,144) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,145) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,150) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,154) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,155) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,163) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,164) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,137) = wrk(:)





































      end subroutine setrxt_hrates

      end module mo_setrxt
