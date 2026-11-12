
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

      rate(:,:,104) = 8.00e-14_r8
      rate(:,:,105) = 3.90e-17_r8
      rate(:,:,108) = 4.20e-13_r8
      rate(:,:,109) = 8.50e-2_r8
      rate(:,:,110) = 1.30e-16_r8
      rate(:,:,112) = 1.00e-20_r8
      rate(:,:,113) = 2.58e-04_r8
      rate(:,:,120) = 1.20e-10_r8
      rate(:,:,121) = 2.02e-10_r8
      rate(:,:,122) = 1.204e-10_r8
      rate(:,:,123) = 1.50e-10_r8
      rate(:,:,124) = 9.75e-11_r8
      rate(:,:,125) = 1.50e-11_r8
      rate(:,:,126) = 7.20e-11_r8
      rate(:,:,127) = 1.794e-10_r8
      rate(:,:,128) = 1.628e-10_r8
      rate(:,:,129) = 2.84e-10_r8
      rate(:,:,130) = 1.674e-10_r8
      rate(:,:,131) = 9.60e-11_r8
      rate(:,:,132) = 4.10e-11_r8
      rate(:,:,133) = 1.012e-10_r8
      rate(:,:,134) = 1.20e-10_r8
      rate(:,:,135) = 4.49e-10_r8
      rate(:,:,136) = 2.57e-10_r8
      rate(:,:,137) = 2.14e-11_r8
      rate(:,:,138) = 1.90e-10_r8
      rate(:,:,139) = 1.31e-10_r8
      rate(:,:,140) = 3.50e-11_r8
      rate(:,:,141) = 9.00e-12_r8
      rate(:,:,142) = 1.20e-10_r8
      rate(:,:,143) = 1.50e-10_r8
      rate(:,:,144) = 1.20e-10_r8
      rate(:,:,147) = 7.20e-11_r8
      rate(:,:,148) = 6.90e-12_r8
      rate(:,:,149) = 1.60e-12_r8
      rate(:,:,153) = 1.80e-12_r8
      rate(:,:,156) = 1.80e-12_r8
      rate(:,:,162) = 5.00e-12_r8
      rate(:,:,163) = 7.00e-13_r8
      rate(:,:,164) = 5.00e-11_r8
      rate(:,:,181) = 1.00e-11_r8
      rate(:,:,182) = 2.20e-11_r8
      rate(:,:,183) = 3.50e-12_r8
      rate(:,:,208) = 1.70e-13_r8
      rate(:,:,280) = 9.0e-10_r8
      rate(:,:,281) = 1.0e-10_r8
      rate(:,:,282) = 4.4e-10_r8
      rate(:,:,283) = 4.0e-10_r8
      rate(:,:,284) = 2.0e-10_r8
      rate(:,:,285) = 1.0e-12_r8
      rate(:,:,286) = 6.0e-11_r8
      rate(:,:,287) = 5.0e-16_r8
      rate(:,:,291) = 4.8e-10_r8
      rate(:,:,292) = 1.0e-10_r8
      rate(:,:,293) = 4.0e-10_r8
      rate(:,:,296) = 5.0e-12_r8
      rate(:,:,297) = 7.0e-10_r8
      rate(:,:,298) = 8.0e-10_r8
      rate(:,:,300) = 4.7e-2_r8
      rate(:,:,301) = 1.71e-1_r8
      rate(:,:,302) = 7.7e-5_r8
      rate(:,:,306) = 4.9e-12_r8
      rate(:,:,308) = 3.9e-11_r8
      rate(:,:,314) = 2.7e-9_r8
      rate(:,:,315) = 8.0e-10_r8
      rate(:,:,318) = 6.e-10_r8
      rate(:,:,319) = 6.e-10_r8
      rate(:,:,320) = 4.e-10_r8
      rate(:,:,321) = 1.e-11_r8
      rate(:,:,322) = 1.e-12_r8
      rate(:,:,323) = 5.e-12_r8
      rate(:,:,335) = 5.0e-12_r8
      rate(:,:,340) = 3e-11_r8
      rate(:,:,341) = 5.0e-11_r8
      rate(:,:,342) = 5.0e-11_r8
      rate(:,:,343) = 1.1e-9_r8
      rate(:,:,344) = 9.2e-10_r8
      rate(:,:,345) = 9.0e-10_r8
      rate(:,:,360) = 1.0e-12_r8
      rate(:,:,363) = 6.7e-12_r8
      rate(:,:,365) = 8.0e-14_r8
      rate(:,:,367) = 2.0e-10_r8
      rate(:,:,368) = 2.0e-12_r8
      rate(:,:,369) = 9.0e-10_r8
      rate(:,:,370) = 1.2e-9_r8
      rate(:,:,371) = 8.2e-10_r8
      rate(:,:,372) = 1.17e-9_r8
      rate(:,:,373) = 5.9e-10_r8
      rate(:,:,375) = 3.5e-12_r8
      rate(:,:,377) = 6.5e-10_r8
      rate(:,:,378) = 1.8e-10_r8
      rate(:,:,379) = 3.3e-10_r8
      rate(:,:,388) = 4.0e-10_r8
      rate(:,:,389) = 5.0e-13_r8
      rate(:,:,394) = 3.3e-10_r8
      rate(:,:,395) = 3.2e-10_r8
      rate(:,:,399) = 6.0e-10_r8
      rate(:,:,401) = 3.2e-10_r8
      rate(:,:,402) = 2.0e-10_r8
      rate(:,:,403) = 2.0e-10_r8
      rate(:,:,408) = 2.0e-10_r8
      rate(:,:,411) = 2.7e-12_r8
      rate(:,:,417) = 5.0e-12_r8
      rate(:,:,418) = 5.0e-12_r8
      rate(:,:,422) = 0.0e-11_r8
      rate(:,:,423) = 2.0e-11_r8
      rate(:,:,424) = 1.7e-11_r8
      rate(:,:,426) = 0.0e-12_r8
      rate(:,:,427) = 6.7e-12_r8
      rate(:,:,428) = 1.2e-11_r8
      rate(:,:,429) = 1.0e-10_r8
      rate(:,:,431) = 2.0e-10_r8
      rate(:,:,432) = 1.5e-10_r8
      rate(:,:,433) = 9.0e-8_r8
      rate(:,:,434) = 9.0e-8_r8
      rate(:,:,435) = 9.0e-8_r8
      rate(:,:,436) = 9.0e-8_r8
      rate(:,:,437) = 9.0e-8_r8
      rate(:,:,438) = 1.8e-9_r8
      rate(:,:,439) = 4.0e-9_r8
      rate(:,:,440) = 3.9e-10_r8
      rate(:,:,441) = 4.2e-11_r8
      rate(:,:,443) = 1.0e-10_r8
      rate(:,:,446) = 1.2e-10_r8
      rate(:,:,447) = 1.3e-9_r8
      rate(:,:,448) = 4.e-10_r8
      rate(:,:,454) = 3e-7_r8
      rate(:,:,456) = 3e-10_r8
      rate(:,:,468) = 2.7e-7_r8
      rate(:,:,469) = 9.4e-10_r8
      rate(:,:,470) = 3.2e-9_r8
      rate(:,:,499) = 1.0e-12_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,102) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,106) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,107) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,111) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,114) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,115) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,116) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,117) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,118) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,119) = 4.63e-11_r8 * exp_fac(:,:)
      rate(:,:,146) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,150) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -940._r8 * itemp(:,:) )
      rate(:,:,151) = 1.70e-12_r8 * exp_fac(:,:)
      rate(:,:,312) = 5.e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,152) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,218) = 1.70e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1800._r8 * itemp(:,:) )
      rate(:,:,155) = 2.80e-12_r8 * exp_fac(:,:)
      rate(:,:,476) = 2.6e-10_r8 * exp_fac(:,:)
      rate(:,:,157) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,158) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,226) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,254) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,159) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,161) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      rate(:,:,165) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,166) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,167) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,168) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,169) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,171) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,190) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,195) = 7.40e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,172) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,227) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,173) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,175) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,201) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,180) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,185) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,187) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,188) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,189) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,191) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,192) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,193) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,194) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,196) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,217) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,225) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,197) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,199) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,224) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,198) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,202) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,203) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,206) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,207) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,209) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,210) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,231) = 1.40e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,211) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,242) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,212) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,213) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,214) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,215) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,216) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,244) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,219) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,220) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,223) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,222) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,228) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,229) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,230) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,280) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,281) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,282) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,283) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,284) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,285) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,286) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,287) = 5.0e-16_r8 * exp_fac(:,:)
      rate(:,:,291) = 4.8e-10_r8 * exp_fac(:,:)
      rate(:,:,292) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,293) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,296) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,297) = 7.0e-10_r8 * exp_fac(:,:)
      rate(:,:,298) = 8.0e-10_r8 * exp_fac(:,:)
      rate(:,:,300) = 4.7e-2_r8 * exp_fac(:,:)
      rate(:,:,301) = 1.71e-1_r8 * exp_fac(:,:)
      rate(:,:,302) = 7.7e-5_r8 * exp_fac(:,:)
      rate(:,:,306) = 4.9e-12_r8 * exp_fac(:,:)
      rate(:,:,308) = 3.9e-11_r8 * exp_fac(:,:)
      rate(:,:,314) = 2.7e-9_r8 * exp_fac(:,:)
      rate(:,:,315) = 8.0e-10_r8 * exp_fac(:,:)
      rate(:,:,318) = 6.e-10_r8 * exp_fac(:,:)
      rate(:,:,319) = 6.e-10_r8 * exp_fac(:,:)
      rate(:,:,320) = 4.e-10_r8 * exp_fac(:,:)
      rate(:,:,321) = 1.e-11_r8 * exp_fac(:,:)
      rate(:,:,322) = 1.e-12_r8 * exp_fac(:,:)
      rate(:,:,323) = 5.e-12_r8 * exp_fac(:,:)
      rate(:,:,335) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,340) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,341) = 5.0e-11_r8 * exp_fac(:,:)
      rate(:,:,342) = 5.0e-11_r8 * exp_fac(:,:)
      rate(:,:,343) = 1.1e-9_r8 * exp_fac(:,:)
      rate(:,:,344) = 9.2e-10_r8 * exp_fac(:,:)
      rate(:,:,345) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,360) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,363) = 6.7e-12_r8 * exp_fac(:,:)
      rate(:,:,365) = 8.0e-14_r8 * exp_fac(:,:)
      rate(:,:,367) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,368) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,369) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,370) = 1.2e-9_r8 * exp_fac(:,:)
      rate(:,:,371) = 8.2e-10_r8 * exp_fac(:,:)
      rate(:,:,372) = 1.17e-9_r8 * exp_fac(:,:)
      rate(:,:,373) = 5.9e-10_r8 * exp_fac(:,:)
      rate(:,:,375) = 3.5e-12_r8 * exp_fac(:,:)
      rate(:,:,377) = 6.5e-10_r8 * exp_fac(:,:)
      rate(:,:,378) = 1.8e-10_r8 * exp_fac(:,:)
      rate(:,:,379) = 3.3e-10_r8 * exp_fac(:,:)
      rate(:,:,388) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,389) = 5.0e-13_r8 * exp_fac(:,:)
      rate(:,:,394) = 3.3e-10_r8 * exp_fac(:,:)
      rate(:,:,395) = 3.2e-10_r8 * exp_fac(:,:)
      rate(:,:,399) = 6.0e-10_r8 * exp_fac(:,:)
      rate(:,:,401) = 3.2e-10_r8 * exp_fac(:,:)
      rate(:,:,402) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,403) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,408) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,411) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,417) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,418) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,422) = 0.0e-11_r8 * exp_fac(:,:)
      rate(:,:,423) = 2.0e-11_r8 * exp_fac(:,:)
      rate(:,:,424) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,426) = 0.0e-12_r8 * exp_fac(:,:)
      rate(:,:,427) = 6.7e-12_r8 * exp_fac(:,:)
      rate(:,:,428) = 1.2e-11_r8 * exp_fac(:,:)
      rate(:,:,429) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,431) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,432) = 1.5e-10_r8 * exp_fac(:,:)
      rate(:,:,433) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,434) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,435) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,436) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,437) = 9.0e-8_r8 * exp_fac(:,:)
      rate(:,:,438) = 1.8e-9_r8 * exp_fac(:,:)
      rate(:,:,439) = 4.0e-9_r8 * exp_fac(:,:)
      rate(:,:,440) = 3.9e-10_r8 * exp_fac(:,:)
      rate(:,:,441) = 4.2e-11_r8 * exp_fac(:,:)
      rate(:,:,443) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,446) = 1.2e-10_r8 * exp_fac(:,:)
      rate(:,:,447) = 1.3e-9_r8 * exp_fac(:,:)
      rate(:,:,448) = 4.e-10_r8 * exp_fac(:,:)
      rate(:,:,454) = 3e-7_r8 * exp_fac(:,:)
      rate(:,:,456) = 3e-10_r8 * exp_fac(:,:)
      rate(:,:,468) = 2.7e-7_r8 * exp_fac(:,:)
      rate(:,:,469) = 9.4e-10_r8 * exp_fac(:,:)
      rate(:,:,470) = 3.2e-9_r8 * exp_fac(:,:)
      rate(:,:,499) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,232) = 1.60e-10_r8 * exp( -260._r8 * itemp(:,:) )
      rate(:,:,233) = 6.00e-12_r8 * exp( 400._r8 * itemp(:,:) )
      rate(:,:,234) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,235) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,236) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,237) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,240) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,251) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,238) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,239) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,241) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -600._r8 * itemp(:,:) )
      rate(:,:,243) = 1.35e-12_r8 * exp_fac(:,:)
      rate(:,:,338) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,366) = 1.0e-11_r8 * exp_fac(:,:)
      rate(:,:,245) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,246) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,249) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,250) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,252) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,253) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,303) = 1.1e-9_r8 * exp( -116._r8 * itemp(:,:) )
      rate(:,:,305) = 3.2e-10_r8 * exp( -550._r8 * itemp(:,:) )
      rate(:,:,307) = 5.06e-10_r8 * exp( -240._r8 * itemp(:,:) )
      rate(:,:,329) = 2.94e-10_r8 * exp( -174._r8 * itemp(:,:) )
      rate(:,:,330) = 4.6e-10_r8 * exp( -350._r8 * itemp(:,:) )
      rate(:,:,331) = 3.0e-10_r8 * exp( -177._r8 * itemp(:,:) )
      rate(:,:,332) = 1.4e-10_r8 * exp( -580._r8 * itemp(:,:) )
      rate(:,:,333) = 4.4e-10_r8 * exp( -170._r8 * itemp(:,:) )
      rate(:,:,334) = 2.3e-10_r8 * exp( -2310._r8 * itemp(:,:) )
      rate(:,:,336) = 3.3e-10_r8 * exp( -302._r8 * itemp(:,:) )
      rate(:,:,337) = 3.0e-10_r8 * exp( -796._r8 * itemp(:,:) )
      rate(:,:,339) = 7.6e-10_r8 * exp( -240.5_r8 * itemp(:,:) )
      rate(:,:,354) = 6.0e-11_r8 * exp( -1200._r8 * itemp(:,:) )
      rate(:,:,355) = 2.28e-10_r8 * exp( -139._r8 * itemp(:,:) )
      rate(:,:,357) = 2.19e-10_r8 * exp( -548._r8 * itemp(:,:) )
      rate(:,:,410) = 8.23e-10_r8 * exp( -192._r8 * itemp(:,:) )
      rate(:,:,412) = 1.1e-9_r8 * exp( -421._r8 * itemp(:,:) )
      rate(:,:,413) = 5.7e-10_r8 * exp( -267._r8 * itemp(:,:) )
      rate(:,:,415) = 4.4e-11_r8 * exp( -202._r8 * itemp(:,:) )
      rate(:,:,421) = 7.0e-10_r8 * exp( -4017._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -120._r8 * itemp(:,:) )
      rate(:,:,457) = 1.15e-9_r8 * exp_fac(:,:)
      rate(:,:,458) = 2.00e-10_r8 * exp_fac(:,:)
      rate(:,:,461) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,462) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,463) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,464) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,465) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,460) = 6.9e-10_r8 * exp( -385._r8 * itemp(:,:) )
      rate(:,:,467) = 4.5e-11_r8 * exp( -3590._r8 * itemp(:,:) )
      rate(:,:,472) = 2.8e-8_r8 * exp( -1680._r8 * itemp(:,:) )
      rate(:,:,474) = 1.5e-9_r8 * exp( -820._r8 * itemp(:,:) )
      rate(:,:,482) = 2.8e-10_r8 * exp( -2220._r8 * itemp(:,:) )
      rate(:,:,484) = 5.0e-10_r8 * exp( -752._r8 * itemp(:,:) )
      rate(:,:,486) = 1.4e-9_r8 * exp( -2200._r8 * itemp(:,:) )
      rate(:,:,489) = 1.0e-9_r8 * exp( -873._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,145), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,154), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,170), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,174), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,176), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,178), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,184), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,200), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,204), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,221), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,248), m, 0.6_r8, ko, kinf, n )

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

      rate(:,:kbot,104) = 8.00e-14_r8
      rate(:,:kbot,105) = 3.90e-17_r8
      rate(:,:kbot,110) = 1.30e-16_r8
      rate(:,:kbot,112) = 1.00e-20_r8
      rate(:,:kbot,148) = 6.90e-12_r8
      rate(:,:kbot,162) = 5.00e-12_r8
      rate(:,:kbot,163) = 7.00e-13_r8
      rate(:,:kbot,281) = 1.0e-10_r8
      rate(:,:kbot,282) = 4.4e-10_r8
      rate(:,:kbot,283) = 4.0e-10_r8
      rate(:,:kbot,284) = 2.0e-10_r8
      rate(:,:kbot,285) = 1.0e-12_r8
      rate(:,:kbot,286) = 6.0e-11_r8
      rate(:,:kbot,291) = 4.8e-10_r8
      rate(:,:kbot,292) = 1.0e-10_r8
      rate(:,:kbot,293) = 4.0e-10_r8
      rate(:,:kbot,296) = 5.0e-12_r8
      rate(:,:kbot,297) = 7.0e-10_r8
      rate(:,:kbot,298) = 8.0e-10_r8
      rate(:,:kbot,300) = 4.7e-2_r8
      rate(:,:kbot,301) = 1.71e-1_r8
      rate(:,:kbot,302) = 7.7e-5_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,102) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,106) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,107) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,111) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,114) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,115) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,116) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,146) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,150) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,151) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,152) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,158) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,159) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,165) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,166) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,171) = 3.30e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,172) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,173) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,145) = wrk(:,:)











      end subroutine setrxt_hrates

      end module mo_setrxt
