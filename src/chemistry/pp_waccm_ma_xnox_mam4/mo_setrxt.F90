
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

      rate(:,109) = 0.000258_r8
      rate(:,110) = 0.085_r8
      rate(:,111) = 1.2e-10_r8
      rate(:,116) = 1.2e-10_r8
      rate(:,117) = 1.2e-10_r8
      rate(:,118) = 1e-20_r8
      rate(:,119) = 1.3e-16_r8
      rate(:,121) = 4.2e-13_r8
      rate(:,123) = 8e-14_r8
      rate(:,124) = 3.9e-17_r8
      rate(:,132) = 1.2e-10_r8
      rate(:,137) = 1.2e-10_r8
      rate(:,143) = 6.9e-12_r8
      rate(:,144) = 7.2e-11_r8
      rate(:,145) = 1.6e-12_r8
      rate(:,154) = 1.8e-12_r8
      rate(:,158) = 1.8e-12_r8
      rate(:,164) = 7e-13_r8
      rate(:,165) = 5e-12_r8
      rate(:,177) = 3.5e-12_r8
      rate(:,179) = 1e-11_r8
      rate(:,180) = 2.2e-11_r8
      rate(:,182) = 1e-11_r8
      rate(:,183) = 5e-11_r8
      rate(:,212) = 7e-13_r8
      rate(:,213) = 5e-12_r8
      rate(:,222) = 3.5e-12_r8
      rate(:,224) = 1e-11_r8
      rate(:,225) = 2.2e-11_r8
      rate(:,226) = 5e-11_r8
      rate(:,261) = 1.7e-13_r8
      rate(:,263) = 1.7e-13_r8
      rate(:,264) = 2.607e-10_r8
      rate(:,265) = 9.75e-11_r8
      rate(:,266) = 2.07e-10_r8
      rate(:,267) = 2.088e-10_r8
      rate(:,268) = 1.17e-10_r8
      rate(:,269) = 4.644e-11_r8
      rate(:,270) = 1.204e-10_r8
      rate(:,271) = 9.9e-11_r8
      rate(:,272) = 3.3e-12_r8
      rate(:,278) = 2.607e-10_r8
      rate(:,279) = 9.75e-11_r8
      rate(:,280) = 2.07e-10_r8
      rate(:,281) = 2.088e-10_r8
      rate(:,282) = 1.17e-10_r8
      rate(:,283) = 4.644e-11_r8
      rate(:,284) = 1.204e-10_r8
      rate(:,285) = 9.9e-11_r8
      rate(:,286) = 3.3e-12_r8
      rate(:,310) = 4.5e-11_r8
      rate(:,311) = 4.62e-10_r8
      rate(:,312) = 1.2e-10_r8
      rate(:,313) = 9e-11_r8
      rate(:,314) = 3e-11_r8
      rate(:,316) = 4.5e-11_r8
      rate(:,317) = 4.62e-10_r8
      rate(:,318) = 1.2e-10_r8
      rate(:,319) = 9e-11_r8
      rate(:,320) = 3e-11_r8
      rate(:,326) = 2.14e-11_r8
      rate(:,327) = 1.9e-10_r8
      rate(:,328) = 2.14e-11_r8
      rate(:,329) = 1.9e-10_r8
      rate(:,342) = 2.57e-10_r8
      rate(:,343) = 1.8e-10_r8
      rate(:,344) = 1.794e-10_r8
      rate(:,345) = 1.3e-10_r8
      rate(:,346) = 7.65e-11_r8
      rate(:,347) = 2.57e-10_r8
      rate(:,348) = 1.8e-10_r8
      rate(:,349) = 1.794e-10_r8
      rate(:,350) = 1.3e-10_r8
      rate(:,351) = 7.65e-11_r8
      rate(:,363) = 1.31e-10_r8
      rate(:,364) = 3.5e-11_r8
      rate(:,365) = 9e-12_r8
      rate(:,367) = 1.31e-10_r8
      rate(:,368) = 3.5e-11_r8
      rate(:,369) = 9e-12_r8
      rate(:,375) = 2.3e-12_r8
      rate(:,376) = 1.2e-11_r8
      rate(:,377) = 5.7e-11_r8
      rate(:,378) = 2.8e-11_r8
      rate(:,379) = 6.6e-11_r8
      rate(:,380) = 1.4e-11_r8
      rate(:,383) = 1.9e-12_r8
      rate(:,385) = 1.4e-11_r8
      rate(:,387) = 1.2e-11_r8
      rate(:,412) = 1._r8
      rate(:,413) = 1._r8
      rate(:,414) = 1._r8
      rate(:,415) = 1._r8
      rate(:,416) = 1._r8
      rate(:,417) = 1._r8
      rate(:,418) = 1._r8
      rate(:,419) = 1._r8
      rate(:,420) = 1._r8
      rate(:,421) = 1._r8
      rate(:,422) = 1._r8
      rate(:,423) = 1._r8
      rate(:,424) = 1._r8
      rate(:,425) = 1._r8
      rate(:,426) = 1._r8
      rate(:,427) = 1._r8
      rate(:,428) = 1._r8
      rate(:,429) = 1._r8
      rate(:,433) = 6e-11_r8
      rate(:,436) = 1e-12_r8
      rate(:,437) = 4e-10_r8
      rate(:,438) = 2e-10_r8
      rate(:,439) = 1e-10_r8
      rate(:,440) = 5e-16_r8
      rate(:,441) = 4.4e-10_r8
      rate(:,442) = 9e-10_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      exp_fac(:) = exp( 60._r8 * itemp(:) )
      rate(:,112) = 1.63e-10_r8 * exp_fac(:)
      rate(:,133) = 1.63e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( 110._r8 * itemp(:) )
      rate(:,113) = 2.15e-11_r8 * exp_fac(:)
      rate(:,134) = 2.15e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,114) = 2.64e-11_r8 * exp_fac(:)
      rate(:,115) = 6.6e-12_r8 * exp_fac(:)
      rate(:,135) = 2.64e-11_r8 * exp_fac(:)
      rate(:,136) = 6.6e-12_r8 * exp_fac(:)
      rate(:,120) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,122) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,125) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      exp_fac(:) = exp( -2060._r8 * itemp(:) )
      rate(:,126) = 8e-12_r8 * exp_fac(:)
      rate(:,127) = 8e-12_r8 * exp_fac(:)
      rate(:,138) = 8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -4570._r8 * itemp(:) )
      rate(:,139) = 1.6e-11_r8 * exp_fac(:)
      rate(:,142) = 1.6e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,140) = 1.4e-12_r8 * exp_fac(:)
      rate(:,141) = 1.4e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,147) = 3e-11_r8 * exp_fac(:)
      rate(:,149) = 3e-11_r8 * exp_fac(:)
      rate(:,306) = 5.5e-12_r8 * exp_fac(:)
      rate(:,360) = 3.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -490._r8 * itemp(:) )
      rate(:,148) = 1e-14_r8 * exp_fac(:)
      rate(:,150) = 1e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( -470._r8 * itemp(:) )
      rate(:,151) = 1.4e-10_r8 * exp_fac(:)
      rate(:,152) = 1.4e-10_r8 * exp_fac(:)
      rate(:,153) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,155) = 4.8e-11_r8 * exp_fac(:)
      rate(:,300) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,156) = 1.8e-11_r8 * exp_fac(:)
      rate(:,160) = 1.8e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -940._r8 * itemp(:) )
      rate(:,157) = 1.7e-12_r8 * exp_fac(:)
      rate(:,161) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 380._r8 * itemp(:) )
      rate(:,163) = 1.3e-12_r8 * exp_fac(:)
      rate(:,211) = 1.3e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 100._r8 * itemp(:) )
      rate(:,166) = 2.1e-11_r8 * exp_fac(:)
      rate(:,189) = 2.1e-11_r8 * exp_fac(:)
      rate(:,214) = 2.1e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,167) = 2.9e-12_r8 * exp_fac(:)
      rate(:,168) = 1.45e-12_r8 * exp_fac(:)
      rate(:,169) = 1.45e-12_r8 * exp_fac(:)
      rate(:,190) = 2.9e-12_r8 * exp_fac(:)
      rate(:,191) = 1.45e-12_r8 * exp_fac(:)
      rate(:,192) = 1.45e-12_r8 * exp_fac(:)
      rate(:,215) = 2.9e-12_r8 * exp_fac(:)
      rate(:,216) = 1.45e-12_r8 * exp_fac(:)
      rate(:,217) = 1.45e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -3600._r8 * itemp(:) )
      rate(:,170) = 1.5e-11_r8 * exp_fac(:)
      rate(:,218) = 1.5e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 210._r8 * itemp(:) )
      rate(:,171) = 5.1e-12_r8 * exp_fac(:)
      rate(:,174) = 5.1e-12_r8 * exp_fac(:)
      rate(:,219) = 5.1e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,172) = 1.2e-13_r8 * exp_fac(:)
      rate(:,175) = 1.2e-13_r8 * exp_fac(:)
      rate(:,220) = 1.2e-13_r8 * exp_fac(:)
      rate(:,240) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 170._r8 * itemp(:) )
      rate(:,178) = 1.5e-11_r8 * exp_fac(:)
      rate(:,181) = 1.5e-11_r8 * exp_fac(:)
      rate(:,223) = 1.5e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,184) = 3.3e-12_r8 * exp_fac(:)
      rate(:,227) = 3.3e-12_r8 * exp_fac(:)
      rate(:,236) = 1.4e-11_r8 * exp_fac(:)
      rate(:,251) = 7.4e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,185) = 3e-12_r8 * exp_fac(:)
      rate(:,187) = 3e-12_r8 * exp_fac(:)
      rate(:,228) = 3e-12_r8 * exp_fac(:)
      rate(:,305) = 5.8e-12_r8 * exp_fac(:)
      rate(:,307) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,193) = 7.26e-11_r8 * exp_fac(:)
      rate(:,194) = 4.64e-11_r8 * exp_fac(:)
      rate(:,230) = 7.26e-11_r8 * exp_fac(:)
      rate(:,231) = 4.64e-11_r8 * exp_fac(:)
      rate(:,232) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,233) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,234) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,235) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,237) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      exp_fac(:) = exp( -200._r8 * itemp(:) )
      rate(:,238) = 2.3e-11_r8 * exp_fac(:)
      rate(:,256) = 2.3e-11_r8 * exp_fac(:)
      rate(:,239) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,241) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,242) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,243) = 2.6e-12_r8 * exp_fac(:)
      rate(:,244) = 6.4e-12_r8 * exp_fac(:)
      rate(:,253) = 6.4e-12_r8 * exp_fac(:)
      rate(:,293) = 4.1e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 135._r8 * itemp(:) )
      rate(:,245) = 6.5e-12_r8 * exp_fac(:)
      rate(:,275) = 6.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,247) = 3.6e-12_r8 * exp_fac(:)
      rate(:,249) = 3.6e-12_r8 * exp_fac(:)
      rate(:,276) = 3.6e-12_r8 * exp_fac(:)
      rate(:,331) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -330._r8 * itemp(:) )
      rate(:,248) = 1.2e-12_r8 * exp_fac(:)
      rate(:,277) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 85._r8 * itemp(:) )
      rate(:,250) = 2.8e-11_r8 * exp_fac(:)
      rate(:,255) = 2.8e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,252) = 6e-13_r8 * exp_fac(:)
      rate(:,290) = 1.5e-12_r8 * exp_fac(:)
      rate(:,299) = 1.9e-11_r8 * exp_fac(:)
      rate(:,303) = 1.9e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -3300._r8 * itemp(:) )
      rate(:,257) = 1e-11_r8 * exp_fac(:)
      rate(:,259) = 1e-11_r8 * exp_fac(:)
      rate(:,258) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,260) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,262) = 3e-12_r8 * exp_fac(:)
      rate(:,322) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,287) = 1.7e-11_r8 * exp_fac(:)
      rate(:,330) = 6.3e-12_r8 * exp_fac(:)
      rate(:,288) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      exp_fac(:) = exp( -780._r8 * itemp(:) )
      rate(:,289) = 1.6e-11_r8 * exp_fac(:)
      rate(:,304) = 1.6e-11_r8 * exp_fac(:)
      rate(:,291) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,292) = 2.3e-12_r8 * exp_fac(:)
      rate(:,295) = 8.8e-12_r8 * exp_fac(:)
      rate(:,301) = 8.8e-12_r8 * exp_fac(:)
      rate(:,294) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      exp_fac(:) = exp( 215._r8 * itemp(:) )
      rate(:,297) = 1.9e-11_r8 * exp_fac(:)
      rate(:,298) = 1.9e-11_r8 * exp_fac(:)
      rate(:,315) = 1.9e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -430._r8 * itemp(:) )
      rate(:,308) = 1.2e-10_r8 * exp_fac(:)
      rate(:,309) = 1.2e-10_r8 * exp_fac(:)
      rate(:,321) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,323) = 1.4e-11_r8 * exp_fac(:)
      rate(:,326) = 2.14e-11_r8 * exp_fac(:)
      rate(:,327) = 1.9e-10_r8 * exp_fac(:)
      rate(:,328) = 2.14e-11_r8 * exp_fac(:)
      rate(:,329) = 1.9e-10_r8 * exp_fac(:)
      rate(:,342) = 2.57e-10_r8 * exp_fac(:)
      rate(:,343) = 1.8e-10_r8 * exp_fac(:)
      rate(:,344) = 1.794e-10_r8 * exp_fac(:)
      rate(:,345) = 1.3e-10_r8 * exp_fac(:)
      rate(:,346) = 7.65e-11_r8 * exp_fac(:)
      rate(:,347) = 2.57e-10_r8 * exp_fac(:)
      rate(:,348) = 1.8e-10_r8 * exp_fac(:)
      rate(:,349) = 1.794e-10_r8 * exp_fac(:)
      rate(:,350) = 1.3e-10_r8 * exp_fac(:)
      rate(:,351) = 7.65e-11_r8 * exp_fac(:)
      rate(:,363) = 1.31e-10_r8 * exp_fac(:)
      rate(:,364) = 3.5e-11_r8 * exp_fac(:)
      rate(:,365) = 9e-12_r8 * exp_fac(:)
      rate(:,367) = 1.31e-10_r8 * exp_fac(:)
      rate(:,368) = 3.5e-11_r8 * exp_fac(:)
      rate(:,369) = 9e-12_r8 * exp_fac(:)
      rate(:,375) = 2.3e-12_r8 * exp_fac(:)
      rate(:,376) = 1.2e-11_r8 * exp_fac(:)
      rate(:,377) = 5.7e-11_r8 * exp_fac(:)
      rate(:,378) = 2.8e-11_r8 * exp_fac(:)
      rate(:,379) = 6.6e-11_r8 * exp_fac(:)
      rate(:,380) = 1.4e-11_r8 * exp_fac(:)
      rate(:,383) = 1.9e-12_r8 * exp_fac(:)
      rate(:,385) = 1.4e-11_r8 * exp_fac(:)
      rate(:,387) = 1.2e-11_r8 * exp_fac(:)
      rate(:,412) = 1._r8 * exp_fac(:)
      rate(:,413) = 1._r8 * exp_fac(:)
      rate(:,414) = 1._r8 * exp_fac(:)
      rate(:,415) = 1._r8 * exp_fac(:)
      rate(:,416) = 1._r8 * exp_fac(:)
      rate(:,417) = 1._r8 * exp_fac(:)
      rate(:,418) = 1._r8 * exp_fac(:)
      rate(:,419) = 1._r8 * exp_fac(:)
      rate(:,420) = 1._r8 * exp_fac(:)
      rate(:,421) = 1._r8 * exp_fac(:)
      rate(:,422) = 1._r8 * exp_fac(:)
      rate(:,423) = 1._r8 * exp_fac(:)
      rate(:,424) = 1._r8 * exp_fac(:)
      rate(:,425) = 1._r8 * exp_fac(:)
      rate(:,426) = 1._r8 * exp_fac(:)
      rate(:,427) = 1._r8 * exp_fac(:)
      rate(:,428) = 1._r8 * exp_fac(:)
      rate(:,429) = 1._r8 * exp_fac(:)
      rate(:,433) = 6e-11_r8 * exp_fac(:)
      rate(:,436) = 1e-12_r8 * exp_fac(:)
      rate(:,437) = 4e-10_r8 * exp_fac(:)
      rate(:,438) = 2e-10_r8 * exp_fac(:)
      rate(:,439) = 1e-10_r8 * exp_fac(:)
      rate(:,440) = 5e-16_r8 * exp_fac(:)
      rate(:,441) = 4.4e-10_r8 * exp_fac(:)
      rate(:,442) = 9e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,324) = 6e-12_r8 * exp_fac(:)
      rate(:,325) = 6e-12_r8 * exp_fac(:)
      rate(:,332) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,333) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      rate(:,334) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:) )
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,335) = 2.03e-11_r8 * exp_fac(:)
      rate(:,382) = 3.4e-12_r8 * exp_fac(:)
      rate(:,386) = 3.4e-12_r8 * exp_fac(:)
      rate(:,336) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,337) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,338) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,339) = 1.25e-12_r8 * exp_fac(:)
      rate(:,353) = 3.4e-11_r8 * exp_fac(:)
      rate(:,356) = 3.4e-11_r8 * exp_fac(:)
      rate(:,340) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,341) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      exp_fac(:) = exp( -2058._r8 * itemp(:) )
      rate(:,352) = 6e-13_r8 * exp_fac(:)
      rate(:,355) = 6e-13_r8 * exp_fac(:)
      rate(:,354) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,357) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,358) = 2.8e-12_r8 * exp_fac(:)
      rate(:,359) = 2.8e-12_r8 * exp_fac(:)
      rate(:,361) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 520._r8 * itemp(:) )
      rate(:,370) = 1.9e-13_r8 * exp_fac(:)
      rate(:,372) = 1.9e-13_r8 * exp_fac(:)
      rate(:,371) = 9.6e-12_r8 * exp( -234._r8 * itemp(:) )
      rate(:,373) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,374) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,381) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,384) = 2.7e-11_r8 * exp( 335._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,146), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,159), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,173), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,176), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,186), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,188), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,195), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,196), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,197), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,198), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,199), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,200), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,201), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,202), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,221), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,229), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,246), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,254), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,273), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,296), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,302), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,362), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,118) = 1e-20_r8
      rate(:n,119) = 1.3e-16_r8
      rate(:n,123) = 8e-14_r8
      rate(:n,124) = 3.9e-17_r8
      rate(:n,143) = 6.9e-12_r8
      rate(:n,164) = 7e-13_r8
      rate(:n,165) = 5e-12_r8
      rate(:n,433) = 6e-11_r8
      rate(:n,436) = 1e-12_r8
      rate(:n,437) = 4e-10_r8
      rate(:n,438) = 2e-10_r8
      rate(:n,439) = 1e-10_r8
      rate(:n,441) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,113) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,114) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,115) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,120) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,122) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,125) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,126) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,147) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,148) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,151) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,155) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,156) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,157) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,166) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,170) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,171) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,184) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,185) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,146) = wrk(:)






















      end subroutine setrxt_hrates

      end module mo_setrxt
