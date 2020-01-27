
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

      rate(:,180) = 0.000258_r8
      rate(:,181) = 0.085_r8
      rate(:,182) = 1.2e-10_r8
      rate(:,187) = 1.2e-10_r8
      rate(:,188) = 1.2e-10_r8
      rate(:,189) = 1e-20_r8
      rate(:,190) = 1.3e-16_r8
      rate(:,192) = 4.2e-13_r8
      rate(:,194) = 8e-14_r8
      rate(:,195) = 3.9e-17_r8
      rate(:,203) = 1.2e-10_r8
      rate(:,208) = 1.2e-10_r8
      rate(:,214) = 6.9e-12_r8
      rate(:,215) = 7.2e-11_r8
      rate(:,216) = 1.6e-12_r8
      rate(:,225) = 1.8e-12_r8
      rate(:,229) = 1.8e-12_r8
      rate(:,235) = 7e-13_r8
      rate(:,236) = 5e-12_r8
      rate(:,248) = 3.5e-12_r8
      rate(:,250) = 1e-11_r8
      rate(:,251) = 2.2e-11_r8
      rate(:,253) = 1e-11_r8
      rate(:,254) = 5e-11_r8
      rate(:,283) = 7e-13_r8
      rate(:,284) = 5e-12_r8
      rate(:,293) = 3.5e-12_r8
      rate(:,295) = 1e-11_r8
      rate(:,296) = 2.2e-11_r8
      rate(:,297) = 5e-11_r8
      rate(:,332) = 1.7e-13_r8
      rate(:,334) = 1.7e-13_r8
      rate(:,335) = 2.607e-10_r8
      rate(:,336) = 9.75e-11_r8
      rate(:,337) = 2.07e-10_r8
      rate(:,338) = 2.088e-10_r8
      rate(:,339) = 1.17e-10_r8
      rate(:,340) = 4.644e-11_r8
      rate(:,341) = 1.204e-10_r8
      rate(:,342) = 9.9e-11_r8
      rate(:,343) = 3.3e-12_r8
      rate(:,349) = 2.607e-10_r8
      rate(:,350) = 9.75e-11_r8
      rate(:,351) = 2.07e-10_r8
      rate(:,352) = 2.088e-10_r8
      rate(:,353) = 1.17e-10_r8
      rate(:,354) = 4.644e-11_r8
      rate(:,355) = 1.204e-10_r8
      rate(:,356) = 9.9e-11_r8
      rate(:,357) = 3.3e-12_r8
      rate(:,381) = 4.5e-11_r8
      rate(:,382) = 4.62e-10_r8
      rate(:,383) = 1.2e-10_r8
      rate(:,384) = 9e-11_r8
      rate(:,385) = 3e-11_r8
      rate(:,387) = 4.5e-11_r8
      rate(:,388) = 4.62e-10_r8
      rate(:,389) = 1.2e-10_r8
      rate(:,390) = 9e-11_r8
      rate(:,391) = 3e-11_r8
      rate(:,397) = 2.14e-11_r8
      rate(:,398) = 1.9e-10_r8
      rate(:,399) = 2.14e-11_r8
      rate(:,400) = 1.9e-10_r8
      rate(:,413) = 2.57e-10_r8
      rate(:,414) = 1.8e-10_r8
      rate(:,415) = 1.794e-10_r8
      rate(:,416) = 1.3e-10_r8
      rate(:,417) = 7.65e-11_r8
      rate(:,418) = 2.57e-10_r8
      rate(:,419) = 1.8e-10_r8
      rate(:,420) = 1.794e-10_r8
      rate(:,421) = 1.3e-10_r8
      rate(:,422) = 7.65e-11_r8
      rate(:,439) = 4e-13_r8
      rate(:,444) = 1.31e-10_r8
      rate(:,445) = 3.5e-11_r8
      rate(:,446) = 9e-12_r8
      rate(:,449) = 1.31e-10_r8
      rate(:,450) = 3.5e-11_r8
      rate(:,451) = 9e-12_r8
      rate(:,458) = 6.8e-14_r8
      rate(:,459) = 2e-13_r8
      rate(:,476) = 7e-13_r8
      rate(:,477) = 1e-12_r8
      rate(:,482) = 1e-14_r8
      rate(:,483) = 1e-11_r8
      rate(:,484) = 1.15e-11_r8
      rate(:,485) = 4e-14_r8
      rate(:,491) = 4e-14_r8
      rate(:,505) = 3e-12_r8
      rate(:,506) = 6.7e-13_r8
      rate(:,518) = 6.7e-13_r8
      rate(:,519) = 3.5e-13_r8
      rate(:,520) = 5.4e-11_r8
      rate(:,521) = 3.5e-13_r8
      rate(:,526) = 2e-12_r8
      rate(:,527) = 1.4e-11_r8
      rate(:,530) = 2.4e-12_r8
      rate(:,533) = 2.4e-12_r8
      rate(:,545) = 5e-12_r8
      rate(:,547) = 5e-12_r8
      rate(:,561) = 2e-12_r8
      rate(:,563) = 1.6e-12_r8
      rate(:,565) = 6.7e-12_r8
      rate(:,567) = 6.7e-12_r8
      rate(:,570) = 3.5e-12_r8
      rate(:,573) = 1.3e-11_r8
      rate(:,574) = 1.4e-11_r8
      rate(:,578) = 2.4e-12_r8
      rate(:,580) = 2.4e-12_r8
      rate(:,581) = 1.4e-11_r8
      rate(:,586) = 2.4e-12_r8
      rate(:,588) = 2.4e-12_r8
      rate(:,589) = 4e-11_r8
      rate(:,590) = 4e-11_r8
      rate(:,592) = 1.4e-11_r8
      rate(:,596) = 2.4e-12_r8
      rate(:,598) = 2.4e-12_r8
      rate(:,599) = 4e-11_r8
      rate(:,605) = 7e-11_r8
      rate(:,606) = 1e-10_r8
      rate(:,607) = 1.6e-12_r8
      rate(:,608) = 4e-11_r8
      rate(:,609) = 4e-11_r8
      rate(:,610) = 1.4e-11_r8
      rate(:,614) = 2.4e-12_r8
      rate(:,615) = 4e-11_r8
      rate(:,616) = 7e-11_r8
      rate(:,617) = 1e-10_r8
      rate(:,622) = 2.4e-12_r8
      rate(:,624) = 2.4e-12_r8
      rate(:,643) = 4.7e-11_r8
      rate(:,663) = 2.1e-12_r8
      rate(:,664) = 2.8e-13_r8
      rate(:,666) = 2.1e-12_r8
      rate(:,667) = 2.8e-13_r8
      rate(:,677) = 1.7e-11_r8
      rate(:,685) = 8.4e-11_r8
      rate(:,687) = 1.9e-11_r8
      rate(:,688) = 1.2e-14_r8
      rate(:,689) = 2e-10_r8
      rate(:,690) = 1.9e-11_r8
      rate(:,691) = 1.2e-14_r8
      rate(:,700) = 2.4e-12_r8
      rate(:,702) = 2.4e-12_r8
      rate(:,703) = 2e-11_r8
      rate(:,708) = 2.3e-11_r8
      rate(:,709) = 2e-11_r8
      rate(:,714) = 3.3e-11_r8
      rate(:,715) = 1e-12_r8
      rate(:,716) = 5.7e-11_r8
      rate(:,717) = 1e-12_r8
      rate(:,718) = 3.4e-11_r8
      rate(:,722) = 2.4e-12_r8
      rate(:,723) = 2e-11_r8
      rate(:,724) = 2e-11_r8
      rate(:,727) = 2.3e-12_r8
      rate(:,728) = 1.2e-11_r8
      rate(:,729) = 5.7e-11_r8
      rate(:,730) = 2.8e-11_r8
      rate(:,731) = 6.6e-11_r8
      rate(:,732) = 1.4e-11_r8
      rate(:,735) = 1.9e-12_r8
      rate(:,737) = 1.4e-11_r8
      rate(:,739) = 1.2e-11_r8
      rate(:,755) = 6.34e-08_r8
      rate(:,774) = 1.9e-11_r8
      rate(:,775) = 1.2e-14_r8
      rate(:,776) = 2e-10_r8
      rate(:,781) = 1.34e-11_r8
      rate(:,785) = 1.34e-11_r8
      rate(:,787) = 1.7e-11_r8
      rate(:,826) = 6e-11_r8
      rate(:,829) = 1e-12_r8
      rate(:,830) = 4e-10_r8
      rate(:,831) = 2e-10_r8
      rate(:,832) = 1e-10_r8
      rate(:,833) = 5e-16_r8
      rate(:,834) = 4.4e-10_r8
      rate(:,835) = 9e-10_r8
      rate(:,838) = 1.29e-07_r8
      rate(:,839) = 2.31e-07_r8
      rate(:,840) = 2.31e-06_r8
      rate(:,841) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      exp_fac(:) = exp( 60._r8 * itemp(:) )
      rate(:,183) = 1.63e-10_r8 * exp_fac(:)
      rate(:,204) = 1.63e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( 110._r8 * itemp(:) )
      rate(:,184) = 2.15e-11_r8 * exp_fac(:)
      rate(:,205) = 2.15e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,185) = 2.64e-11_r8 * exp_fac(:)
      rate(:,186) = 6.6e-12_r8 * exp_fac(:)
      rate(:,206) = 2.64e-11_r8 * exp_fac(:)
      rate(:,207) = 6.6e-12_r8 * exp_fac(:)
      rate(:,191) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,193) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,196) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      exp_fac(:) = exp( -2060._r8 * itemp(:) )
      rate(:,197) = 8e-12_r8 * exp_fac(:)
      rate(:,198) = 8e-12_r8 * exp_fac(:)
      rate(:,209) = 8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -4570._r8 * itemp(:) )
      rate(:,210) = 1.6e-11_r8 * exp_fac(:)
      rate(:,213) = 1.6e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,211) = 1.4e-12_r8 * exp_fac(:)
      rate(:,212) = 1.4e-12_r8 * exp_fac(:)
      rate(:,600) = 1.05e-14_r8 * exp_fac(:)
      rate(:,604) = 1.05e-14_r8 * exp_fac(:)
      rate(:,779) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,218) = 3e-11_r8 * exp_fac(:)
      rate(:,220) = 3e-11_r8 * exp_fac(:)
      rate(:,377) = 5.5e-12_r8 * exp_fac(:)
      rate(:,435) = 3.8e-12_r8 * exp_fac(:)
      rate(:,464) = 3.8e-12_r8 * exp_fac(:)
      rate(:,500) = 3.8e-12_r8 * exp_fac(:)
      rate(:,510) = 3.8e-12_r8 * exp_fac(:)
      rate(:,515) = 3.8e-12_r8 * exp_fac(:)
      rate(:,538) = 2.3e-11_r8 * exp_fac(:)
      rate(:,552) = 3.8e-12_r8 * exp_fac(:)
      rate(:,569) = 3.8e-12_r8 * exp_fac(:)
      rate(:,602) = 1.52e-11_r8 * exp_fac(:)
      rate(:,625) = 1.52e-12_r8 * exp_fac(:)
      rate(:,633) = 3.8e-12_r8 * exp_fac(:)
      rate(:,636) = 3.8e-12_r8 * exp_fac(:)
      rate(:,642) = 3.8e-12_r8 * exp_fac(:)
      rate(:,665) = 3.8e-12_r8 * exp_fac(:)
      rate(:,673) = 3.8e-12_r8 * exp_fac(:)
      rate(:,681) = 3.8e-12_r8 * exp_fac(:)
      rate(:,686) = 3.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -490._r8 * itemp(:) )
      rate(:,219) = 1e-14_r8 * exp_fac(:)
      rate(:,221) = 1e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( -470._r8 * itemp(:) )
      rate(:,222) = 1.4e-10_r8 * exp_fac(:)
      rate(:,223) = 1.4e-10_r8 * exp_fac(:)
      rate(:,224) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,226) = 4.8e-11_r8 * exp_fac(:)
      rate(:,372) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,227) = 1.8e-11_r8 * exp_fac(:)
      rate(:,231) = 1.8e-11_r8 * exp_fac(:)
      rate(:,479) = 4.2e-12_r8 * exp_fac(:)
      rate(:,480) = 4.2e-12_r8 * exp_fac(:)
      rate(:,498) = 4.2e-12_r8 * exp_fac(:)
      rate(:,499) = 4.2e-12_r8 * exp_fac(:)
      rate(:,508) = 4.2e-12_r8 * exp_fac(:)
      rate(:,509) = 4.2e-12_r8 * exp_fac(:)
      rate(:,549) = 4.2e-12_r8 * exp_fac(:)
      rate(:,550) = 4.2e-12_r8 * exp_fac(:)
      rate(:,577) = 4.4e-12_r8 * exp_fac(:)
      rate(:,579) = 4.4e-12_r8 * exp_fac(:)
      rate(:,585) = 4.4e-12_r8 * exp_fac(:)
      rate(:,587) = 4.4e-12_r8 * exp_fac(:)
      rate(:,699) = 4.2e-12_r8 * exp_fac(:)
      rate(:,701) = 4.2e-12_r8 * exp_fac(:)
      rate(:,706) = 4.2e-12_r8 * exp_fac(:)
      rate(:,707) = 4.2e-12_r8 * exp_fac(:)
      rate(:,712) = 4.2e-12_r8 * exp_fac(:)
      rate(:,713) = 4.2e-12_r8 * exp_fac(:)
      rate(:,721) = 4.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -940._r8 * itemp(:) )
      rate(:,228) = 1.7e-12_r8 * exp_fac(:)
      rate(:,232) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 380._r8 * itemp(:) )
      rate(:,234) = 1.3e-12_r8 * exp_fac(:)
      rate(:,282) = 1.3e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 100._r8 * itemp(:) )
      rate(:,237) = 2.1e-11_r8 * exp_fac(:)
      rate(:,260) = 2.1e-11_r8 * exp_fac(:)
      rate(:,285) = 2.1e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,238) = 2.9e-12_r8 * exp_fac(:)
      rate(:,239) = 1.45e-12_r8 * exp_fac(:)
      rate(:,240) = 1.45e-12_r8 * exp_fac(:)
      rate(:,261) = 2.9e-12_r8 * exp_fac(:)
      rate(:,262) = 1.45e-12_r8 * exp_fac(:)
      rate(:,263) = 1.45e-12_r8 * exp_fac(:)
      rate(:,286) = 2.9e-12_r8 * exp_fac(:)
      rate(:,287) = 1.45e-12_r8 * exp_fac(:)
      rate(:,288) = 1.45e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -3600._r8 * itemp(:) )
      rate(:,241) = 1.5e-11_r8 * exp_fac(:)
      rate(:,289) = 1.5e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 210._r8 * itemp(:) )
      rate(:,242) = 5.1e-12_r8 * exp_fac(:)
      rate(:,245) = 5.1e-12_r8 * exp_fac(:)
      rate(:,290) = 5.1e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,243) = 1.2e-13_r8 * exp_fac(:)
      rate(:,246) = 1.2e-13_r8 * exp_fac(:)
      rate(:,291) = 1.2e-13_r8 * exp_fac(:)
      rate(:,311) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 170._r8 * itemp(:) )
      rate(:,249) = 1.5e-11_r8 * exp_fac(:)
      rate(:,252) = 1.5e-11_r8 * exp_fac(:)
      rate(:,294) = 1.5e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,255) = 3.3e-12_r8 * exp_fac(:)
      rate(:,298) = 3.3e-12_r8 * exp_fac(:)
      rate(:,307) = 1.4e-11_r8 * exp_fac(:)
      rate(:,322) = 7.4e-12_r8 * exp_fac(:)
      rate(:,474) = 8.1e-12_r8 * exp_fac(:)
      rate(:,475) = 8.1e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,256) = 3e-12_r8 * exp_fac(:)
      rate(:,258) = 3e-12_r8 * exp_fac(:)
      rate(:,299) = 3e-12_r8 * exp_fac(:)
      rate(:,376) = 5.8e-12_r8 * exp_fac(:)
      rate(:,378) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,264) = 7.26e-11_r8 * exp_fac(:)
      rate(:,265) = 4.64e-11_r8 * exp_fac(:)
      rate(:,301) = 7.26e-11_r8 * exp_fac(:)
      rate(:,302) = 4.64e-11_r8 * exp_fac(:)
      rate(:,303) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,304) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,305) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,306) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,308) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      exp_fac(:) = exp( -200._r8 * itemp(:) )
      rate(:,309) = 2.3e-11_r8 * exp_fac(:)
      rate(:,327) = 2.3e-11_r8 * exp_fac(:)
      rate(:,310) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,312) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,313) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,314) = 2.6e-12_r8 * exp_fac(:)
      rate(:,315) = 6.4e-12_r8 * exp_fac(:)
      rate(:,324) = 6.4e-12_r8 * exp_fac(:)
      rate(:,364) = 4.1e-13_r8 * exp_fac(:)
      rate(:,627) = 7.5e-12_r8 * exp_fac(:)
      rate(:,628) = 7.5e-12_r8 * exp_fac(:)
      rate(:,645) = 7.5e-12_r8 * exp_fac(:)
      rate(:,647) = 7.5e-12_r8 * exp_fac(:)
      rate(:,650) = 7.5e-12_r8 * exp_fac(:)
      rate(:,652) = 7.5e-12_r8 * exp_fac(:)
      rate(:,655) = 7.5e-12_r8 * exp_fac(:)
      rate(:,657) = 7.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 135._r8 * itemp(:) )
      rate(:,316) = 6.5e-12_r8 * exp_fac(:)
      rate(:,346) = 6.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,318) = 3.6e-12_r8 * exp_fac(:)
      rate(:,320) = 3.6e-12_r8 * exp_fac(:)
      rate(:,347) = 3.6e-12_r8 * exp_fac(:)
      rate(:,402) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -330._r8 * itemp(:) )
      rate(:,319) = 1.2e-12_r8 * exp_fac(:)
      rate(:,348) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 85._r8 * itemp(:) )
      rate(:,321) = 2.8e-11_r8 * exp_fac(:)
      rate(:,326) = 2.8e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,323) = 6e-13_r8 * exp_fac(:)
      rate(:,361) = 1.5e-12_r8 * exp_fac(:)
      rate(:,371) = 1.9e-11_r8 * exp_fac(:)
      rate(:,374) = 1.9e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -3300._r8 * itemp(:) )
      rate(:,328) = 1e-11_r8 * exp_fac(:)
      rate(:,330) = 1e-11_r8 * exp_fac(:)
      rate(:,329) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,331) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,333) = 3e-12_r8 * exp_fac(:)
      rate(:,393) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,358) = 1.7e-11_r8 * exp_fac(:)
      rate(:,401) = 6.3e-12_r8 * exp_fac(:)
      rate(:,359) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      exp_fac(:) = exp( -780._r8 * itemp(:) )
      rate(:,360) = 1.6e-11_r8 * exp_fac(:)
      rate(:,375) = 1.6e-11_r8 * exp_fac(:)
      rate(:,362) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,363) = 2.3e-12_r8 * exp_fac(:)
      rate(:,366) = 8.8e-12_r8 * exp_fac(:)
      rate(:,367) = 8.8e-12_r8 * exp_fac(:)
      rate(:,365) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      exp_fac(:) = exp( 215._r8 * itemp(:) )
      rate(:,369) = 1.9e-11_r8 * exp_fac(:)
      rate(:,370) = 1.9e-11_r8 * exp_fac(:)
      rate(:,386) = 1.9e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -430._r8 * itemp(:) )
      rate(:,379) = 1.2e-10_r8 * exp_fac(:)
      rate(:,380) = 1.2e-10_r8 * exp_fac(:)
      rate(:,392) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,394) = 1.4e-11_r8 * exp_fac(:)
      rate(:,397) = 2.14e-11_r8 * exp_fac(:)
      rate(:,398) = 1.9e-10_r8 * exp_fac(:)
      rate(:,399) = 2.14e-11_r8 * exp_fac(:)
      rate(:,400) = 1.9e-10_r8 * exp_fac(:)
      rate(:,413) = 2.57e-10_r8 * exp_fac(:)
      rate(:,414) = 1.8e-10_r8 * exp_fac(:)
      rate(:,415) = 1.794e-10_r8 * exp_fac(:)
      rate(:,416) = 1.3e-10_r8 * exp_fac(:)
      rate(:,417) = 7.65e-11_r8 * exp_fac(:)
      rate(:,418) = 2.57e-10_r8 * exp_fac(:)
      rate(:,419) = 1.8e-10_r8 * exp_fac(:)
      rate(:,420) = 1.794e-10_r8 * exp_fac(:)
      rate(:,421) = 1.3e-10_r8 * exp_fac(:)
      rate(:,422) = 7.65e-11_r8 * exp_fac(:)
      rate(:,439) = 4e-13_r8 * exp_fac(:)
      rate(:,444) = 1.31e-10_r8 * exp_fac(:)
      rate(:,445) = 3.5e-11_r8 * exp_fac(:)
      rate(:,446) = 9e-12_r8 * exp_fac(:)
      rate(:,449) = 1.31e-10_r8 * exp_fac(:)
      rate(:,450) = 3.5e-11_r8 * exp_fac(:)
      rate(:,451) = 9e-12_r8 * exp_fac(:)
      rate(:,458) = 6.8e-14_r8 * exp_fac(:)
      rate(:,459) = 2e-13_r8 * exp_fac(:)
      rate(:,476) = 7e-13_r8 * exp_fac(:)
      rate(:,477) = 1e-12_r8 * exp_fac(:)
      rate(:,482) = 1e-14_r8 * exp_fac(:)
      rate(:,483) = 1e-11_r8 * exp_fac(:)
      rate(:,484) = 1.15e-11_r8 * exp_fac(:)
      rate(:,485) = 4e-14_r8 * exp_fac(:)
      rate(:,491) = 4e-14_r8 * exp_fac(:)
      rate(:,505) = 3e-12_r8 * exp_fac(:)
      rate(:,506) = 6.7e-13_r8 * exp_fac(:)
      rate(:,518) = 6.7e-13_r8 * exp_fac(:)
      rate(:,519) = 3.5e-13_r8 * exp_fac(:)
      rate(:,520) = 5.4e-11_r8 * exp_fac(:)
      rate(:,521) = 3.5e-13_r8 * exp_fac(:)
      rate(:,526) = 2e-12_r8 * exp_fac(:)
      rate(:,527) = 1.4e-11_r8 * exp_fac(:)
      rate(:,530) = 2.4e-12_r8 * exp_fac(:)
      rate(:,533) = 2.4e-12_r8 * exp_fac(:)
      rate(:,545) = 5e-12_r8 * exp_fac(:)
      rate(:,547) = 5e-12_r8 * exp_fac(:)
      rate(:,561) = 2e-12_r8 * exp_fac(:)
      rate(:,563) = 1.6e-12_r8 * exp_fac(:)
      rate(:,565) = 6.7e-12_r8 * exp_fac(:)
      rate(:,567) = 6.7e-12_r8 * exp_fac(:)
      rate(:,570) = 3.5e-12_r8 * exp_fac(:)
      rate(:,573) = 1.3e-11_r8 * exp_fac(:)
      rate(:,574) = 1.4e-11_r8 * exp_fac(:)
      rate(:,578) = 2.4e-12_r8 * exp_fac(:)
      rate(:,580) = 2.4e-12_r8 * exp_fac(:)
      rate(:,581) = 1.4e-11_r8 * exp_fac(:)
      rate(:,586) = 2.4e-12_r8 * exp_fac(:)
      rate(:,588) = 2.4e-12_r8 * exp_fac(:)
      rate(:,589) = 4e-11_r8 * exp_fac(:)
      rate(:,590) = 4e-11_r8 * exp_fac(:)
      rate(:,592) = 1.4e-11_r8 * exp_fac(:)
      rate(:,596) = 2.4e-12_r8 * exp_fac(:)
      rate(:,598) = 2.4e-12_r8 * exp_fac(:)
      rate(:,599) = 4e-11_r8 * exp_fac(:)
      rate(:,605) = 7e-11_r8 * exp_fac(:)
      rate(:,606) = 1e-10_r8 * exp_fac(:)
      rate(:,607) = 1.6e-12_r8 * exp_fac(:)
      rate(:,608) = 4e-11_r8 * exp_fac(:)
      rate(:,609) = 4e-11_r8 * exp_fac(:)
      rate(:,610) = 1.4e-11_r8 * exp_fac(:)
      rate(:,614) = 2.4e-12_r8 * exp_fac(:)
      rate(:,615) = 4e-11_r8 * exp_fac(:)
      rate(:,616) = 7e-11_r8 * exp_fac(:)
      rate(:,617) = 1e-10_r8 * exp_fac(:)
      rate(:,622) = 2.4e-12_r8 * exp_fac(:)
      rate(:,624) = 2.4e-12_r8 * exp_fac(:)
      rate(:,643) = 4.7e-11_r8 * exp_fac(:)
      rate(:,663) = 2.1e-12_r8 * exp_fac(:)
      rate(:,664) = 2.8e-13_r8 * exp_fac(:)
      rate(:,666) = 2.1e-12_r8 * exp_fac(:)
      rate(:,667) = 2.8e-13_r8 * exp_fac(:)
      rate(:,677) = 1.7e-11_r8 * exp_fac(:)
      rate(:,685) = 8.4e-11_r8 * exp_fac(:)
      rate(:,687) = 1.9e-11_r8 * exp_fac(:)
      rate(:,688) = 1.2e-14_r8 * exp_fac(:)
      rate(:,689) = 2e-10_r8 * exp_fac(:)
      rate(:,690) = 1.9e-11_r8 * exp_fac(:)
      rate(:,691) = 1.2e-14_r8 * exp_fac(:)
      rate(:,700) = 2.4e-12_r8 * exp_fac(:)
      rate(:,702) = 2.4e-12_r8 * exp_fac(:)
      rate(:,703) = 2e-11_r8 * exp_fac(:)
      rate(:,708) = 2.3e-11_r8 * exp_fac(:)
      rate(:,709) = 2e-11_r8 * exp_fac(:)
      rate(:,714) = 3.3e-11_r8 * exp_fac(:)
      rate(:,715) = 1e-12_r8 * exp_fac(:)
      rate(:,716) = 5.7e-11_r8 * exp_fac(:)
      rate(:,717) = 1e-12_r8 * exp_fac(:)
      rate(:,718) = 3.4e-11_r8 * exp_fac(:)
      rate(:,722) = 2.4e-12_r8 * exp_fac(:)
      rate(:,723) = 2e-11_r8 * exp_fac(:)
      rate(:,724) = 2e-11_r8 * exp_fac(:)
      rate(:,727) = 2.3e-12_r8 * exp_fac(:)
      rate(:,728) = 1.2e-11_r8 * exp_fac(:)
      rate(:,729) = 5.7e-11_r8 * exp_fac(:)
      rate(:,730) = 2.8e-11_r8 * exp_fac(:)
      rate(:,731) = 6.6e-11_r8 * exp_fac(:)
      rate(:,732) = 1.4e-11_r8 * exp_fac(:)
      rate(:,735) = 1.9e-12_r8 * exp_fac(:)
      rate(:,737) = 1.4e-11_r8 * exp_fac(:)
      rate(:,739) = 1.2e-11_r8 * exp_fac(:)
      rate(:,755) = 6.34e-08_r8 * exp_fac(:)
      rate(:,774) = 1.9e-11_r8 * exp_fac(:)
      rate(:,775) = 1.2e-14_r8 * exp_fac(:)
      rate(:,776) = 2e-10_r8 * exp_fac(:)
      rate(:,781) = 1.34e-11_r8 * exp_fac(:)
      rate(:,785) = 1.34e-11_r8 * exp_fac(:)
      rate(:,787) = 1.7e-11_r8 * exp_fac(:)
      rate(:,826) = 6e-11_r8 * exp_fac(:)
      rate(:,829) = 1e-12_r8 * exp_fac(:)
      rate(:,830) = 4e-10_r8 * exp_fac(:)
      rate(:,831) = 2e-10_r8 * exp_fac(:)
      rate(:,832) = 1e-10_r8 * exp_fac(:)
      rate(:,833) = 5e-16_r8 * exp_fac(:)
      rate(:,834) = 4.4e-10_r8 * exp_fac(:)
      rate(:,835) = 9e-10_r8 * exp_fac(:)
      rate(:,838) = 1.29e-07_r8 * exp_fac(:)
      rate(:,839) = 2.31e-07_r8 * exp_fac(:)
      rate(:,840) = 2.31e-06_r8 * exp_fac(:)
      rate(:,841) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,395) = 6e-12_r8 * exp_fac(:)
      rate(:,396) = 6e-12_r8 * exp_fac(:)
      rate(:,528) = 5e-13_r8 * exp_fac(:)
      rate(:,575) = 5e-13_r8 * exp_fac(:)
      rate(:,582) = 5e-13_r8 * exp_fac(:)
      rate(:,593) = 5e-13_r8 * exp_fac(:)
      rate(:,611) = 5e-13_r8 * exp_fac(:)
      rate(:,619) = 5e-13_r8 * exp_fac(:)
      rate(:,403) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,404) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,405) = 1.64e-12_r8 * exp_fac(:)
      rate(:,554) = 8.5e-16_r8 * exp_fac(:)
      rate(:,556) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,406) = 2.03e-11_r8 * exp_fac(:)
      rate(:,734) = 3.4e-12_r8 * exp_fac(:)
      rate(:,738) = 3.4e-12_r8 * exp_fac(:)
      rate(:,407) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,408) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,409) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,410) = 1.25e-12_r8 * exp_fac(:)
      rate(:,425) = 3.4e-11_r8 * exp_fac(:)
      rate(:,428) = 3.4e-11_r8 * exp_fac(:)
      rate(:,411) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,412) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,423) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      exp_fac(:) = exp( -2058._r8 * itemp(:) )
      rate(:,424) = 6e-13_r8 * exp_fac(:)
      rate(:,427) = 6e-13_r8 * exp_fac(:)
      rate(:,426) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,429) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,430) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,431) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,432) = 2.8e-12_r8 * exp_fac(:)
      rate(:,433) = 2.8e-12_r8 * exp_fac(:)
      rate(:,513) = 2.9e-12_r8 * exp_fac(:)
      rate(:,514) = 2.9e-12_r8 * exp_fac(:)
      rate(:,434) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,436) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,440) = 7.5e-13_r8 * exp_fac(:)
      rate(:,460) = 7.5e-13_r8 * exp_fac(:)
      rate(:,478) = 7.5e-13_r8 * exp_fac(:)
      rate(:,497) = 7.5e-13_r8 * exp_fac(:)
      rate(:,507) = 7.5e-13_r8 * exp_fac(:)
      rate(:,512) = 8.6e-13_r8 * exp_fac(:)
      rate(:,529) = 8e-13_r8 * exp_fac(:)
      rate(:,548) = 7.5e-13_r8 * exp_fac(:)
      rate(:,564) = 7.5e-13_r8 * exp_fac(:)
      rate(:,576) = 8e-13_r8 * exp_fac(:)
      rate(:,583) = 8e-13_r8 * exp_fac(:)
      rate(:,594) = 8e-13_r8 * exp_fac(:)
      rate(:,612) = 8e-13_r8 * exp_fac(:)
      rate(:,620) = 8e-13_r8 * exp_fac(:)
      rate(:,630) = 7.5e-13_r8 * exp_fac(:)
      rate(:,635) = 7.5e-13_r8 * exp_fac(:)
      rate(:,639) = 7.5e-13_r8 * exp_fac(:)
      rate(:,659) = 7.5e-13_r8 * exp_fac(:)
      rate(:,670) = 7.5e-13_r8 * exp_fac(:)
      rate(:,678) = 7.5e-13_r8 * exp_fac(:)
      rate(:,682) = 7.5e-13_r8 * exp_fac(:)
      rate(:,698) = 7.5e-13_r8 * exp_fac(:)
      rate(:,705) = 7.5e-13_r8 * exp_fac(:)
      rate(:,711) = 7.5e-13_r8 * exp_fac(:)
      rate(:,720) = 7.5e-13_r8 * exp_fac(:)
      rate(:,441) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      exp_fac(:) = exp( 265._r8 * itemp(:) )
      rate(:,442) = 2.6e-12_r8 * exp_fac(:)
      rate(:,443) = 2.6e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 105._r8 * itemp(:) )
      rate(:,447) = 1.08e-10_r8 * exp_fac(:)
      rate(:,452) = 1.08e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -2630._r8 * itemp(:) )
      rate(:,456) = 1.2e-14_r8 * exp_fac(:)
      rate(:,457) = 1.2e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,461) = 2.6e-12_r8 * exp_fac(:)
      rate(:,462) = 2.6e-12_r8 * exp_fac(:)
      rate(:,631) = 2.6e-12_r8 * exp_fac(:)
      rate(:,632) = 2.6e-12_r8 * exp_fac(:)
      rate(:,637) = 2.6e-12_r8 * exp_fac(:)
      rate(:,638) = 2.6e-12_r8 * exp_fac(:)
      rate(:,640) = 2.6e-12_r8 * exp_fac(:)
      rate(:,641) = 2.6e-12_r8 * exp_fac(:)
      rate(:,660) = 2.6e-12_r8 * exp_fac(:)
      rate(:,661) = 2.6e-12_r8 * exp_fac(:)
      rate(:,671) = 2.6e-12_r8 * exp_fac(:)
      rate(:,672) = 2.6e-12_r8 * exp_fac(:)
      rate(:,679) = 2.6e-12_r8 * exp_fac(:)
      rate(:,680) = 2.6e-12_r8 * exp_fac(:)
      rate(:,683) = 2.6e-12_r8 * exp_fac(:)
      rate(:,684) = 2.6e-12_r8 * exp_fac(:)
      rate(:,463) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,465) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,466) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,467) = 1.4e-12_r8 * exp_fac(:)
      rate(:,469) = 1.4e-12_r8 * exp_fac(:)
      rate(:,493) = 6.5e-15_r8 * exp_fac(:)
      rate(:,495) = 6.5e-15_r8 * exp_fac(:)
      rate(:,468) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      rate(:,470) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,471) = 2.9e-12_r8 * exp_fac(:)
      rate(:,472) = 2e-12_r8 * exp_fac(:)
      rate(:,511) = 7.1e-13_r8 * exp_fac(:)
      rate(:,541) = 2e-12_r8 * exp_fac(:)
      rate(:,697) = 2e-12_r8 * exp_fac(:)
      rate(:,704) = 2e-12_r8 * exp_fac(:)
      rate(:,710) = 2e-12_r8 * exp_fac(:)
      rate(:,719) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,473) = 4.3e-13_r8 * exp_fac(:)
      rate(:,542) = 4.3e-13_r8 * exp_fac(:)
      rate(:,626) = 4.3e-13_r8 * exp_fac(:)
      rate(:,644) = 4.3e-13_r8 * exp_fac(:)
      rate(:,649) = 4.3e-13_r8 * exp_fac(:)
      rate(:,654) = 4.3e-13_r8 * exp_fac(:)
      rate(:,481) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      exp_fac(:) = exp( -1156._r8 * itemp(:) )
      rate(:,492) = 4.6e-13_r8 * exp_fac(:)
      rate(:,494) = 4.6e-13_r8 * exp_fac(:)
      rate(:,496) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,501) = 8.7e-12_r8 * exp( -615._r8 * itemp(:) )
      exp_fac(:) = exp( -1860._r8 * itemp(:) )
      rate(:,502) = 1.4e-12_r8 * exp_fac(:)
      rate(:,504) = 1.4e-12_r8 * exp_fac(:)
      rate(:,503) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      exp_fac(:) = exp( 120._r8 * itemp(:) )
      rate(:,522) = 4.8e-12_r8 * exp_fac(:)
      rate(:,524) = 4.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 693._r8 * itemp(:) )
      rate(:,523) = 5.1e-14_r8 * exp_fac(:)
      rate(:,525) = 5.1e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,531) = 2.7e-12_r8 * exp_fac(:)
      rate(:,532) = 1.3e-13_r8 * exp_fac(:)
      rate(:,534) = 2.7e-12_r8 * exp_fac(:)
      rate(:,535) = 1.3e-13_r8 * exp_fac(:)
      rate(:,537) = 9.6e-12_r8 * exp_fac(:)
      rate(:,544) = 5.3e-12_r8 * exp_fac(:)
      rate(:,546) = 5.3e-12_r8 * exp_fac(:)
      rate(:,595) = 2.7e-12_r8 * exp_fac(:)
      rate(:,597) = 2.7e-12_r8 * exp_fac(:)
      rate(:,613) = 2.7e-12_r8 * exp_fac(:)
      rate(:,621) = 2.7e-12_r8 * exp_fac(:)
      rate(:,623) = 2.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -2100._r8 * itemp(:) )
      rate(:,536) = 1.5e-15_r8 * exp_fac(:)
      rate(:,539) = 1.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,540) = 4.6e-12_r8 * exp_fac(:)
      rate(:,543) = 2.3e-12_r8 * exp_fac(:)
      rate(:,551) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,555) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      exp_fac(:) = exp( 870._r8 * itemp(:) )
      rate(:,566) = 5.4e-14_r8 * exp_fac(:)
      rate(:,568) = 5.4e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,571) = 1.86e-11_r8 * exp_fac(:)
      rate(:,572) = 1.86e-11_r8 * exp_fac(:)
      rate(:,584) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,591) = 3.03e-12_r8 * exp_fac(:)
      rate(:,603) = 3.03e-12_r8 * exp_fac(:)
      rate(:,778) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,601) = 2.54e-11_r8 * exp_fac(:)
      rate(:,780) = 2.54e-11_r8 * exp_fac(:)
      rate(:,618) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,629) = 2.3e-12_r8 * exp_fac(:)
      rate(:,777) = 2.3e-12_r8 * exp_fac(:)
      rate(:,634) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,662) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,674) = 1.7e-12_r8 * exp_fac(:)
      rate(:,786) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,692) = 1.2e-12_r8 * exp_fac(:)
      rate(:,695) = 1.2e-12_r8 * exp_fac(:)
      rate(:,782) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,693) = 6.3e-16_r8 * exp_fac(:)
      rate(:,696) = 6.3e-16_r8 * exp_fac(:)
      rate(:,783) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,694) = 1.2e-11_r8 * exp_fac(:)
      rate(:,784) = 1.2e-11_r8 * exp_fac(:)
      rate(:,725) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,726) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,733) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,736) = 2.7e-11_r8 * exp( 335._r8 * itemp(:) )
      exp_fac(:) = exp( 520._r8 * itemp(:) )
      rate(:,742) = 1.9e-13_r8 * exp_fac(:)
      rate(:,744) = 1.9e-13_r8 * exp_fac(:)
      rate(:,743) = 9.6e-12_r8 * exp( -234._r8 * itemp(:) )
      rate(:,745) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,217), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,230), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,244), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,247), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,257), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,259), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,266), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,267), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,268), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,269), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,270), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,271), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,272), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,273), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,292), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,300), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,317), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,325), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,344), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,368), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,373), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,437), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.28e-33_r8
      kinf(:) = 9.3e-15_r8 * itemp(:)**(-4.42_r8)
      call jpl( rate(:,438), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,453), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,454), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,455), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,486), m, 0.48_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,487), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,488), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,516), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,553), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,562), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,646), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,648), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,651), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,653), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,656), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,658), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,668), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,669), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,189) = 1e-20_r8
      rate(:n,190) = 1.3e-16_r8
      rate(:n,194) = 8e-14_r8
      rate(:n,195) = 3.9e-17_r8
      rate(:n,214) = 6.9e-12_r8
      rate(:n,235) = 7e-13_r8
      rate(:n,236) = 5e-12_r8
      rate(:n,826) = 6e-11_r8
      rate(:n,829) = 1e-12_r8
      rate(:n,830) = 4e-10_r8
      rate(:n,831) = 2e-10_r8
      rate(:n,832) = 1e-10_r8
      rate(:n,834) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,184) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,185) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,186) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,191) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,193) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,196) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,197) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,218) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,219) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,222) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,226) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,227) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,228) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,237) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,241) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,242) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,255) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,256) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,217) = wrk(:)








































      end subroutine setrxt_hrates

      end module mo_setrxt
