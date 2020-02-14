      module mo_nln_matrix
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only: veclen
      private
      public :: nlnmat
      contains
      subroutine nlnmat01( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,631) = -(rxt(k,563)*y(k,250))
         mat(k,1884) = -rxt(k,563)*y(k,1)
         mat(k,2390) = rxt(k,566)*y(k,219)
         mat(k,1059) = rxt(k,566)*y(k,124)
         mat(k,654) = -(rxt(k,569)*y(k,250))
         mat(k,1885) = -rxt(k,569)*y(k,2)
         mat(k,1060) = rxt(k,564)*y(k,232)
         mat(k,2130) = rxt(k,564)*y(k,219)
         mat(k,1005) = -(rxt(k,687)*y(k,126) + rxt(k,688)*y(k,135) + rxt(k,689) &
                      *y(k,250))
         mat(k,2314) = -rxt(k,687)*y(k,6)
         mat(k,2218) = -rxt(k,688)*y(k,6)
         mat(k,1912) = -rxt(k,689)*y(k,6)
         mat(k,91) = -(rxt(k,629)*y(k,250))
         mat(k,1805) = -rxt(k,629)*y(k,7)
         mat(k,362) = -(rxt(k,633)*y(k,250))
         mat(k,1852) = -rxt(k,633)*y(k,8)
         mat(k,589) = rxt(k,630)*y(k,232)
         mat(k,2105) = rxt(k,630)*y(k,220)
         mat(k,92) = .120_r8*rxt(k,629)*y(k,250)
         mat(k,1806) = .120_r8*rxt(k,629)*y(k,7)
         mat(k,999) = .100_r8*rxt(k,688)*y(k,135)
         mat(k,1031) = .100_r8*rxt(k,693)*y(k,135)
         mat(k,2206) = .100_r8*rxt(k,688)*y(k,6) + .100_r8*rxt(k,693)*y(k,110)
         mat(k,2374) = .500_r8*rxt(k,631)*y(k,220) + .200_r8*rxt(k,671)*y(k,257) &
                      + .060_r8*rxt(k,679)*y(k,261)
         mat(k,590) = .500_r8*rxt(k,631)*y(k,124)
         mat(k,784) = .200_r8*rxt(k,671)*y(k,124)
         mat(k,810) = .060_r8*rxt(k,679)*y(k,124)
         mat(k,2368) = .200_r8*rxt(k,671)*y(k,257) + .200_r8*rxt(k,679)*y(k,261)
         mat(k,783) = .200_r8*rxt(k,671)*y(k,124)
         mat(k,808) = .200_r8*rxt(k,679)*y(k,124)
         mat(k,2384) = .200_r8*rxt(k,671)*y(k,257) + .150_r8*rxt(k,679)*y(k,261)
         mat(k,786) = .200_r8*rxt(k,671)*y(k,124)
         mat(k,811) = .150_r8*rxt(k,679)*y(k,124)
         mat(k,2369) = .210_r8*rxt(k,679)*y(k,261)
         mat(k,809) = .210_r8*rxt(k,679)*y(k,124)
         mat(k,180) = -(rxt(k,570)*y(k,250))
         mat(k,1823) = -rxt(k,570)*y(k,15)
         mat(k,998) = .050_r8*rxt(k,688)*y(k,135)
         mat(k,1030) = .050_r8*rxt(k,693)*y(k,135)
         mat(k,2204) = .050_r8*rxt(k,688)*y(k,6) + .050_r8*rxt(k,693)*y(k,110)
         mat(k,406) = -(rxt(k,519)*y(k,126) + rxt(k,520)*y(k,250))
         mat(k,2306) = -rxt(k,519)*y(k,16)
         mat(k,1858) = -rxt(k,520)*y(k,16)
         mat(k,1657) = -(rxt(k,358)*y(k,42) + rxt(k,359)*y(k,232) + rxt(k,360) &
                      *y(k,135))
         mat(k,1775) = -rxt(k,358)*y(k,17)
         mat(k,2179) = -rxt(k,359)*y(k,17)
         mat(k,2248) = -rxt(k,360)*y(k,17)
         mat(k,2568) = 4.000_r8*rxt(k,361)*y(k,19) + (rxt(k,362)+rxt(k,363))*y(k,59) &
                      + rxt(k,366)*y(k,124) + rxt(k,371)*y(k,133) + rxt(k,729) &
                      *y(k,150) + rxt(k,372)*y(k,250)
         mat(k,1980) = (rxt(k,362)+rxt(k,363))*y(k,19)
         mat(k,849) = rxt(k,376)*y(k,133) + rxt(k,384)*y(k,246) + rxt(k,377)*y(k,250)
         mat(k,2439) = rxt(k,366)*y(k,19)
         mat(k,2490) = rxt(k,371)*y(k,19) + rxt(k,376)*y(k,81)
         mat(k,1616) = rxt(k,729)*y(k,19)
         mat(k,2280) = rxt(k,384)*y(k,81)
         mat(k,1947) = rxt(k,372)*y(k,19) + rxt(k,377)*y(k,81)
         mat(k,2555) = rxt(k,364)*y(k,59)
         mat(k,1967) = rxt(k,364)*y(k,19)
         mat(k,1630) = (rxt(k,796)+rxt(k,801))*y(k,91)
         mat(k,728) = (rxt(k,796)+rxt(k,801))*y(k,85)
         mat(k,2583) = -(4._r8*rxt(k,361)*y(k,19) + (rxt(k,362) + rxt(k,363) + rxt(k,364) &
                      ) * y(k,59) + rxt(k,365)*y(k,232) + rxt(k,366)*y(k,124) &
                      + rxt(k,368)*y(k,125) + rxt(k,371)*y(k,133) + rxt(k,372) &
                      *y(k,250) + rxt(k,729)*y(k,150))
         mat(k,1996) = -(rxt(k,362) + rxt(k,363) + rxt(k,364)) * y(k,19)
         mat(k,2195) = -rxt(k,365)*y(k,19)
         mat(k,2455) = -rxt(k,366)*y(k,19)
         mat(k,2046) = -rxt(k,368)*y(k,19)
         mat(k,2506) = -rxt(k,371)*y(k,19)
         mat(k,1963) = -rxt(k,372)*y(k,19)
         mat(k,1627) = -rxt(k,729)*y(k,19)
         mat(k,1666) = rxt(k,360)*y(k,135)
         mat(k,512) = rxt(k,369)*y(k,133)
         mat(k,854) = rxt(k,385)*y(k,246)
         mat(k,736) = rxt(k,379)*y(k,133)
         mat(k,2506) = mat(k,2506) + rxt(k,369)*y(k,20) + rxt(k,379)*y(k,91)
         mat(k,2264) = rxt(k,360)*y(k,17)
         mat(k,2296) = rxt(k,385)*y(k,81)
         mat(k,505) = -(rxt(k,369)*y(k,133))
         mat(k,2468) = -rxt(k,369)*y(k,20)
         mat(k,2559) = rxt(k,368)*y(k,125)
         mat(k,2007) = rxt(k,368)*y(k,19)
         mat(k,186) = -(rxt(k,634)*y(k,250))
         mat(k,1824) = -rxt(k,634)*y(k,22)
         mat(k,2367) = rxt(k,637)*y(k,221)
         mat(k,513) = rxt(k,637)*y(k,124)
         mat(k,287) = -(rxt(k,636)*y(k,250))
         mat(k,1839) = -rxt(k,636)*y(k,23)
         mat(k,514) = rxt(k,635)*y(k,232)
         mat(k,2098) = rxt(k,635)*y(k,221)
         mat(k,253) = -(rxt(k,453)*y(k,56) + rxt(k,454)*y(k,250))
         mat(k,2049) = -rxt(k,453)*y(k,24)
         mat(k,1836) = -rxt(k,454)*y(k,24)
         mat(k,496) = -(rxt(k,455)*y(k,56) + rxt(k,456)*y(k,135) + rxt(k,486)*y(k,250))
         mat(k,2054) = -rxt(k,455)*y(k,25)
         mat(k,2209) = -rxt(k,456)*y(k,25)
         mat(k,1871) = -rxt(k,486)*y(k,25)
         mat(k,217) = -(rxt(k,463)*y(k,250))
         mat(k,1830) = -rxt(k,463)*y(k,26)
         mat(k,1118) = .800_r8*rxt(k,458)*y(k,222) + .200_r8*rxt(k,459)*y(k,226)
         mat(k,1670) = .200_r8*rxt(k,459)*y(k,222)
         mat(k,292) = -(rxt(k,464)*y(k,250))
         mat(k,1840) = -rxt(k,464)*y(k,27)
         mat(k,1119) = rxt(k,460)*y(k,232)
         mat(k,2099) = rxt(k,460)*y(k,222)
         mat(k,259) = -(rxt(k,465)*y(k,56) + rxt(k,466)*y(k,250))
         mat(k,2050) = -rxt(k,465)*y(k,28)
         mat(k,1837) = -rxt(k,466)*y(k,28)
         mat(k,1155) = -(rxt(k,492)*y(k,126) + rxt(k,493)*y(k,135) + rxt(k,516) &
                      *y(k,250))
         mat(k,2324) = -rxt(k,492)*y(k,29)
         mat(k,2228) = -rxt(k,493)*y(k,29)
         mat(k,1922) = -rxt(k,516)*y(k,29)
         mat(k,967) = .130_r8*rxt(k,600)*y(k,135)
         mat(k,2228) = mat(k,2228) + .130_r8*rxt(k,600)*y(k,98)
         mat(k,374) = -(rxt(k,500)*y(k,250))
         mat(k,1854) = -rxt(k,500)*y(k,30)
         mat(k,898) = rxt(k,497)*y(k,232)
         mat(k,2107) = rxt(k,497)*y(k,223)
         mat(k,57) = -(rxt(k,501)*y(k,250))
         mat(k,1802) = -rxt(k,501)*y(k,31)
         mat(k,221) = -(rxt(k,642)*y(k,250))
         mat(k,1831) = -rxt(k,642)*y(k,32)
         mat(k,738) = rxt(k,639)*y(k,232)
         mat(k,2095) = rxt(k,639)*y(k,224)
         mat(k,1778) = -(rxt(k,303)*y(k,56) + rxt(k,358)*y(k,17) + rxt(k,423)*y(k,232) &
                      + rxt(k,424)*y(k,126) + rxt(k,425)*y(k,133) + rxt(k,426) &
                      *y(k,250))
         mat(k,2075) = -rxt(k,303)*y(k,42)
         mat(k,1659) = -rxt(k,358)*y(k,42)
         mat(k,2182) = -rxt(k,423)*y(k,42)
         mat(k,2350) = -rxt(k,424)*y(k,42)
         mat(k,2493) = -rxt(k,425)*y(k,42)
         mat(k,1950) = -rxt(k,426)*y(k,42)
         mat(k,637) = .400_r8*rxt(k,563)*y(k,250)
         mat(k,1019) = .340_r8*rxt(k,688)*y(k,135)
         mat(k,412) = .500_r8*rxt(k,519)*y(k,126)
         mat(k,500) = rxt(k,456)*y(k,135)
         mat(k,1163) = .500_r8*rxt(k,493)*y(k,135)
         mat(k,456) = .500_r8*rxt(k,477)*y(k,250)
         mat(k,780) = rxt(k,434)*y(k,250)
         mat(k,388) = .300_r8*rxt(k,435)*y(k,250)
         mat(k,1983) = rxt(k,310)*y(k,226)
         mat(k,1192) = .800_r8*rxt(k,483)*y(k,250)
         mat(k,976) = .910_r8*rxt(k,600)*y(k,135)
         mat(k,584) = .300_r8*rxt(k,589)*y(k,250)
         mat(k,1390) = .800_r8*rxt(k,593)*y(k,226)
         mat(k,1407) = .120_r8*rxt(k,536)*y(k,135)
         mat(k,476) = .500_r8*rxt(k,553)*y(k,250)
         mat(k,1051) = .340_r8*rxt(k,693)*y(k,135)
         mat(k,1537) = .600_r8*rxt(k,554)*y(k,135)
         mat(k,2442) = .100_r8*rxt(k,565)*y(k,219) + rxt(k,432)*y(k,226) &
                      + .500_r8*rxt(k,522)*y(k,229) + .500_r8*rxt(k,479)*y(k,231) &
                      + .920_r8*rxt(k,577)*y(k,234) + .250_r8*rxt(k,531)*y(k,236) &
                      + rxt(k,544)*y(k,238) + rxt(k,508)*y(k,253) + rxt(k,513) &
                      *y(k,254) + .340_r8*rxt(k,706)*y(k,255) + .320_r8*rxt(k,712) &
                      *y(k,256) + .250_r8*rxt(k,621)*y(k,260)
         mat(k,2350) = mat(k,2350) + .500_r8*rxt(k,519)*y(k,16) + rxt(k,578)*y(k,234) &
                      + .250_r8*rxt(k,530)*y(k,236) + rxt(k,545)*y(k,238)
         mat(k,2251) = .340_r8*rxt(k,688)*y(k,6) + rxt(k,456)*y(k,25) &
                      + .500_r8*rxt(k,493)*y(k,29) + .910_r8*rxt(k,600)*y(k,98) &
                      + .120_r8*rxt(k,536)*y(k,105) + .340_r8*rxt(k,693)*y(k,110) &
                      + .600_r8*rxt(k,554)*y(k,111)
         mat(k,442) = rxt(k,485)*y(k,250)
         mat(k,1230) = .680_r8*rxt(k,718)*y(k,250)
         mat(k,1069) = .100_r8*rxt(k,565)*y(k,124)
         mat(k,1128) = .700_r8*rxt(k,459)*y(k,226)
         mat(k,906) = rxt(k,496)*y(k,226)
         mat(k,1597) = rxt(k,472)*y(k,226) + rxt(k,574)*y(k,234) + .250_r8*rxt(k,527) &
                      *y(k,236) + rxt(k,540)*y(k,238) + .250_r8*rxt(k,618)*y(k,260)
         mat(k,1713) = rxt(k,310)*y(k,59) + .800_r8*rxt(k,593)*y(k,101) + rxt(k,432) &
                      *y(k,124) + .700_r8*rxt(k,459)*y(k,222) + rxt(k,496)*y(k,223) &
                      + rxt(k,472)*y(k,225) + (4.000_r8*rxt(k,429)+2.000_r8*rxt(k,430)) &
                      *y(k,226) + 1.500_r8*rxt(k,575)*y(k,234) + .750_r8*rxt(k,582) &
                      *y(k,235) + .880_r8*rxt(k,528)*y(k,236) + 2.000_r8*rxt(k,541) &
                      *y(k,238) + .750_r8*rxt(k,697)*y(k,245) + .800_r8*rxt(k,511) &
                      *y(k,254) + .930_r8*rxt(k,704)*y(k,255) + .950_r8*rxt(k,710) &
                      *y(k,256) + .800_r8*rxt(k,619)*y(k,260)
         mat(k,705) = .500_r8*rxt(k,522)*y(k,124)
         mat(k,887) = .500_r8*rxt(k,479)*y(k,124)
         mat(k,2182) = mat(k,2182) + .450_r8*rxt(k,542)*y(k,238) + .150_r8*rxt(k,512) &
                      *y(k,254)
         mat(k,1487) = .920_r8*rxt(k,577)*y(k,124) + rxt(k,578)*y(k,126) + rxt(k,574) &
                      *y(k,225) + 1.500_r8*rxt(k,575)*y(k,226)
         mat(k,1454) = .750_r8*rxt(k,582)*y(k,226)
         mat(k,1514) = .250_r8*rxt(k,531)*y(k,124) + .250_r8*rxt(k,530)*y(k,126) &
                      + .250_r8*rxt(k,527)*y(k,225) + .880_r8*rxt(k,528)*y(k,226)
         mat(k,1559) = rxt(k,544)*y(k,124) + rxt(k,545)*y(k,126) + rxt(k,540)*y(k,225) &
                      + 2.000_r8*rxt(k,541)*y(k,226) + .450_r8*rxt(k,542)*y(k,232) &
                      + 4.000_r8*rxt(k,543)*y(k,238)
         mat(k,1309) = .750_r8*rxt(k,697)*y(k,226)
         mat(k,1950) = mat(k,1950) + .400_r8*rxt(k,563)*y(k,1) + .500_r8*rxt(k,477) &
                      *y(k,51) + rxt(k,434)*y(k,52) + .300_r8*rxt(k,435)*y(k,53) &
                      + .800_r8*rxt(k,483)*y(k,74) + .300_r8*rxt(k,589)*y(k,99) &
                      + .500_r8*rxt(k,553)*y(k,109) + rxt(k,485)*y(k,139) &
                      + .680_r8*rxt(k,718)*y(k,178)
         mat(k,863) = rxt(k,508)*y(k,124)
         mat(k,1328) = rxt(k,513)*y(k,124) + .800_r8*rxt(k,511)*y(k,226) &
                      + .150_r8*rxt(k,512)*y(k,232)
         mat(k,1287) = .340_r8*rxt(k,706)*y(k,124) + .930_r8*rxt(k,704)*y(k,226)
         mat(k,1262) = .320_r8*rxt(k,712)*y(k,124) + .950_r8*rxt(k,710)*y(k,226)
         mat(k,1365) = .250_r8*rxt(k,621)*y(k,124) + .250_r8*rxt(k,618)*y(k,225) &
                      + .800_r8*rxt(k,619)*y(k,226)
      end do
      end subroutine nlnmat01
      subroutine nlnmat02( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,1198) = -(rxt(k,467)*y(k,126) + rxt(k,468)*y(k,250))
         mat(k,2327) = -rxt(k,467)*y(k,45)
         mat(k,1925) = -rxt(k,468)*y(k,45)
         mat(k,635) = .800_r8*rxt(k,563)*y(k,250)
         mat(k,411) = rxt(k,519)*y(k,126)
         mat(k,218) = rxt(k,463)*y(k,250)
         mat(k,294) = .500_r8*rxt(k,464)*y(k,250)
         mat(k,1156) = .500_r8*rxt(k,493)*y(k,135)
         mat(k,1527) = .100_r8*rxt(k,554)*y(k,135)
         mat(k,2420) = .400_r8*rxt(k,565)*y(k,219) + rxt(k,461)*y(k,222) &
                      + .270_r8*rxt(k,498)*y(k,223) + rxt(k,522)*y(k,229) + rxt(k,549) &
                      *y(k,240) + rxt(k,508)*y(k,253)
         mat(k,2327) = mat(k,2327) + rxt(k,519)*y(k,16)
         mat(k,2230) = .500_r8*rxt(k,493)*y(k,29) + .100_r8*rxt(k,554)*y(k,111)
         mat(k,1067) = .400_r8*rxt(k,565)*y(k,124)
         mat(k,1125) = rxt(k,461)*y(k,124) + 3.200_r8*rxt(k,458)*y(k,222) &
                      + .800_r8*rxt(k,459)*y(k,226)
         mat(k,903) = .270_r8*rxt(k,498)*y(k,124)
         mat(k,1693) = .800_r8*rxt(k,459)*y(k,222)
         mat(k,703) = rxt(k,522)*y(k,124)
         mat(k,2159) = .200_r8*rxt(k,548)*y(k,240)
         mat(k,722) = rxt(k,549)*y(k,124) + .200_r8*rxt(k,548)*y(k,232)
         mat(k,1925) = mat(k,1925) + .800_r8*rxt(k,563)*y(k,1) + rxt(k,463)*y(k,26) &
                      + .500_r8*rxt(k,464)*y(k,27)
         mat(k,860) = rxt(k,508)*y(k,124)
         mat(k,52) = -(rxt(k,470)*y(k,250))
         mat(k,1801) = -rxt(k,470)*y(k,47)
         mat(k,1074) = -(rxt(k,517)*y(k,250))
         mat(k,1915) = -rxt(k,517)*y(k,48)
         mat(k,633) = .800_r8*rxt(k,563)*y(k,250)
         mat(k,1007) = .520_r8*rxt(k,688)*y(k,135)
         mat(k,408) = .500_r8*rxt(k,519)*y(k,126)
         mat(k,1039) = .520_r8*rxt(k,693)*y(k,135)
         mat(k,2411) = .250_r8*rxt(k,565)*y(k,219) + .820_r8*rxt(k,498)*y(k,223) &
                      + .500_r8*rxt(k,522)*y(k,229) + .270_r8*rxt(k,706)*y(k,255) &
                      + .040_r8*rxt(k,712)*y(k,256)
         mat(k,2317) = .500_r8*rxt(k,519)*y(k,16)
         mat(k,2221) = .520_r8*rxt(k,688)*y(k,6) + .520_r8*rxt(k,693)*y(k,110)
         mat(k,1224) = .500_r8*rxt(k,718)*y(k,250)
         mat(k,1063) = .250_r8*rxt(k,565)*y(k,124)
         mat(k,900) = .820_r8*rxt(k,498)*y(k,124) + .820_r8*rxt(k,496)*y(k,226)
         mat(k,1684) = .820_r8*rxt(k,496)*y(k,223) + .150_r8*rxt(k,704)*y(k,255) &
                      + .025_r8*rxt(k,710)*y(k,256)
         mat(k,699) = .500_r8*rxt(k,522)*y(k,124)
         mat(k,1915) = mat(k,1915) + .800_r8*rxt(k,563)*y(k,1) + .500_r8*rxt(k,718) &
                      *y(k,178)
         mat(k,1275) = .270_r8*rxt(k,706)*y(k,124) + .150_r8*rxt(k,704)*y(k,226)
         mat(k,1248) = .040_r8*rxt(k,712)*y(k,124) + .025_r8*rxt(k,710)*y(k,226)
         mat(k,1417) = -(rxt(k,502)*y(k,126) + rxt(k,503)*y(k,250))
         mat(k,2339) = -rxt(k,502)*y(k,49)
         mat(k,1938) = -rxt(k,503)*y(k,49)
         mat(k,1239) = rxt(k,505)*y(k,250)
         mat(k,1403) = .880_r8*rxt(k,536)*y(k,135)
         mat(k,1530) = .500_r8*rxt(k,554)*y(k,135)
         mat(k,2432) = .170_r8*rxt(k,645)*y(k,227) + .050_r8*rxt(k,585)*y(k,235) &
                      + .250_r8*rxt(k,531)*y(k,236) + .170_r8*rxt(k,655)*y(k,239) &
                      + .400_r8*rxt(k,671)*y(k,257) + .250_r8*rxt(k,621)*y(k,260) &
                      + .540_r8*rxt(k,679)*y(k,261) + .510_r8*rxt(k,683)*y(k,262)
         mat(k,2339) = mat(k,2339) + .050_r8*rxt(k,586)*y(k,235) + .250_r8*rxt(k,530) &
                      *y(k,236) + .250_r8*rxt(k,622)*y(k,260)
         mat(k,893) = rxt(k,506)*y(k,250)
         mat(k,2240) = .880_r8*rxt(k,536)*y(k,105) + .500_r8*rxt(k,554)*y(k,111)
         mat(k,1588) = .250_r8*rxt(k,527)*y(k,236) + .250_r8*rxt(k,618)*y(k,260)
         mat(k,1704) = .240_r8*rxt(k,528)*y(k,236) + .500_r8*rxt(k,511)*y(k,254) &
                      + .100_r8*rxt(k,619)*y(k,260)
         mat(k,831) = .170_r8*rxt(k,645)*y(k,124) + .070_r8*rxt(k,644)*y(k,232)
         mat(k,2171) = .070_r8*rxt(k,644)*y(k,227) + .070_r8*rxt(k,654)*y(k,239)
         mat(k,1447) = .050_r8*rxt(k,585)*y(k,124) + .050_r8*rxt(k,586)*y(k,126)
         mat(k,1509) = .250_r8*rxt(k,531)*y(k,124) + .250_r8*rxt(k,530)*y(k,126) &
                      + .250_r8*rxt(k,527)*y(k,225) + .240_r8*rxt(k,528)*y(k,226)
         mat(k,987) = .170_r8*rxt(k,655)*y(k,124) + .070_r8*rxt(k,654)*y(k,232)
         mat(k,1938) = mat(k,1938) + rxt(k,505)*y(k,95) + rxt(k,506)*y(k,127)
         mat(k,1325) = .500_r8*rxt(k,511)*y(k,226)
         mat(k,795) = .400_r8*rxt(k,671)*y(k,124)
         mat(k,1362) = .250_r8*rxt(k,621)*y(k,124) + .250_r8*rxt(k,622)*y(k,126) &
                      + .250_r8*rxt(k,618)*y(k,225) + .100_r8*rxt(k,619)*y(k,226)
         mat(k,821) = .540_r8*rxt(k,679)*y(k,124)
         mat(k,619) = .510_r8*rxt(k,683)*y(k,124)
         mat(k,469) = -(rxt(k,476)*y(k,250))
         mat(k,1866) = -rxt(k,476)*y(k,50)
         mat(k,1146) = .120_r8*rxt(k,493)*y(k,135)
         mat(k,2207) = .120_r8*rxt(k,493)*y(k,29)
         mat(k,1575) = .100_r8*rxt(k,472)*y(k,226) + .150_r8*rxt(k,473)*y(k,232)
         mat(k,1678) = .100_r8*rxt(k,472)*y(k,225)
         mat(k,2116) = .150_r8*rxt(k,473)*y(k,225) + .150_r8*rxt(k,542)*y(k,238)
         mat(k,1549) = .150_r8*rxt(k,542)*y(k,232)
         mat(k,453) = -(rxt(k,477)*y(k,250))
         mat(k,1864) = -rxt(k,477)*y(k,51)
         mat(k,1574) = .400_r8*rxt(k,473)*y(k,232)
         mat(k,2115) = .400_r8*rxt(k,473)*y(k,225) + .400_r8*rxt(k,542)*y(k,238)
         mat(k,1548) = .400_r8*rxt(k,542)*y(k,232)
         mat(k,779) = -(rxt(k,434)*y(k,250))
         mat(k,1894) = -rxt(k,434)*y(k,52)
         mat(k,1375) = .200_r8*rxt(k,593)*y(k,226)
         mat(k,1120) = .300_r8*rxt(k,459)*y(k,226)
         mat(k,1680) = .200_r8*rxt(k,593)*y(k,101) + .300_r8*rxt(k,459)*y(k,222) &
                      + 2.000_r8*rxt(k,430)*y(k,226) + .250_r8*rxt(k,575)*y(k,234) &
                      + .250_r8*rxt(k,582)*y(k,235) + .250_r8*rxt(k,528)*y(k,236) &
                      + .250_r8*rxt(k,697)*y(k,245) + .500_r8*rxt(k,511)*y(k,254) &
                      + .250_r8*rxt(k,704)*y(k,255) + .250_r8*rxt(k,710)*y(k,256) &
                      + .300_r8*rxt(k,619)*y(k,260)
         mat(k,1468) = .250_r8*rxt(k,575)*y(k,226)
         mat(k,1432) = .250_r8*rxt(k,582)*y(k,226)
         mat(k,1499) = .250_r8*rxt(k,528)*y(k,226)
         mat(k,1297) = .250_r8*rxt(k,697)*y(k,226)
         mat(k,1319) = .500_r8*rxt(k,511)*y(k,226)
         mat(k,1274) = .250_r8*rxt(k,704)*y(k,226)
         mat(k,1247) = .250_r8*rxt(k,710)*y(k,226)
         mat(k,1353) = .300_r8*rxt(k,619)*y(k,226)
         mat(k,386) = -(rxt(k,435)*y(k,250))
         mat(k,1855) = -rxt(k,435)*y(k,53)
         mat(k,1676) = rxt(k,431)*y(k,232)
         mat(k,2108) = rxt(k,431)*y(k,226)
         mat(k,2079) = -(rxt(k,303)*y(k,42) + rxt(k,305)*y(k,77) + rxt(k,306)*y(k,79) &
                      + (rxt(k,307) + rxt(k,308)) * y(k,232) + rxt(k,309)*y(k,135) &
                      + rxt(k,316)*y(k,60) + rxt(k,331)*y(k,92) + rxt(k,465)*y(k,28))
         mat(k,1782) = -rxt(k,303)*y(k,56)
         mat(k,1345) = -rxt(k,305)*y(k,56)
         mat(k,576) = -rxt(k,306)*y(k,56)
         mat(k,2186) = -(rxt(k,307) + rxt(k,308)) * y(k,56)
         mat(k,2255) = -rxt(k,309)*y(k,56)
         mat(k,949) = -rxt(k,316)*y(k,56)
         mat(k,843) = -rxt(k,331)*y(k,56)
         mat(k,263) = -rxt(k,465)*y(k,56)
         mat(k,2574) = rxt(k,363)*y(k,59)
         mat(k,1987) = rxt(k,363)*y(k,19) + (4.000_r8*rxt(k,311)+2.000_r8*rxt(k,313)) &
                      *y(k,59) + rxt(k,315)*y(k,124) + rxt(k,321)*y(k,133) &
                      + rxt(k,730)*y(k,150) + rxt(k,310)*y(k,226) + rxt(k,322) &
                      *y(k,250)
         mat(k,210) = rxt(k,398)*y(k,246)
         mat(k,1645) = rxt(k,328)*y(k,133) + rxt(k,342)*y(k,246) + rxt(k,329)*y(k,250)
         mat(k,2446) = rxt(k,315)*y(k,59)
         mat(k,2497) = rxt(k,321)*y(k,59) + rxt(k,328)*y(k,85)
         mat(k,1621) = rxt(k,730)*y(k,59)
         mat(k,1717) = rxt(k,310)*y(k,59)
         mat(k,2287) = rxt(k,398)*y(k,65) + rxt(k,342)*y(k,85)
         mat(k,1954) = rxt(k,322)*y(k,59) + rxt(k,329)*y(k,85)
         mat(k,2048) = rxt(k,316)*y(k,60)
         mat(k,1966) = 2.000_r8*rxt(k,312)*y(k,59)
         mat(k,941) = rxt(k,316)*y(k,56) + (rxt(k,794)+rxt(k,799)+rxt(k,804))*y(k,85)
         mat(k,1629) = (rxt(k,794)+rxt(k,799)+rxt(k,804))*y(k,60) + (rxt(k,789) &
                       +rxt(k,795)+rxt(k,800))*y(k,92)
         mat(k,837) = (rxt(k,789)+rxt(k,795)+rxt(k,800))*y(k,85)
         mat(k,1965) = 2.000_r8*rxt(k,344)*y(k,59)
         mat(k,1985) = -(rxt(k,310)*y(k,226) + (4._r8*rxt(k,311) + 4._r8*rxt(k,312) &
                      + 4._r8*rxt(k,313) + 4._r8*rxt(k,344)) * y(k,59) + rxt(k,314) &
                      *y(k,232) + rxt(k,315)*y(k,124) + rxt(k,317)*y(k,125) + rxt(k,321) &
                      *y(k,133) + (rxt(k,322) + rxt(k,323)) * y(k,250) + (rxt(k,362) &
                      + rxt(k,363) + rxt(k,364)) * y(k,19) + rxt(k,730)*y(k,150))
         mat(k,1715) = -rxt(k,310)*y(k,59)
         mat(k,2184) = -rxt(k,314)*y(k,59)
         mat(k,2444) = -rxt(k,315)*y(k,59)
         mat(k,2035) = -rxt(k,317)*y(k,59)
         mat(k,2495) = -rxt(k,321)*y(k,59)
         mat(k,1952) = -(rxt(k,322) + rxt(k,323)) * y(k,59)
         mat(k,2572) = -(rxt(k,362) + rxt(k,363) + rxt(k,364)) * y(k,59)
         mat(k,1619) = -rxt(k,730)*y(k,59)
         mat(k,2077) = rxt(k,331)*y(k,92) + rxt(k,309)*y(k,135) + rxt(k,308)*y(k,232)
         mat(k,947) = rxt(k,318)*y(k,133)
         mat(k,1643) = rxt(k,343)*y(k,246)
         mat(k,842) = rxt(k,331)*y(k,56) + rxt(k,332)*y(k,133) + rxt(k,333)*y(k,250)
         mat(k,2495) = mat(k,2495) + rxt(k,318)*y(k,60) + rxt(k,332)*y(k,92)
         mat(k,2253) = rxt(k,309)*y(k,56)
         mat(k,279) = rxt(k,735)*y(k,150)
         mat(k,1619) = mat(k,1619) + rxt(k,735)*y(k,136)
         mat(k,2184) = mat(k,2184) + rxt(k,308)*y(k,56)
         mat(k,2285) = rxt(k,343)*y(k,85)
         mat(k,1952) = mat(k,1952) + rxt(k,333)*y(k,92)
         mat(k,944) = -(rxt(k,316)*y(k,56) + rxt(k,318)*y(k,133) + rxt(k,319)*y(k,250) &
                      + (rxt(k,794) + rxt(k,799) + rxt(k,804)) * y(k,85))
         mat(k,2060) = -rxt(k,316)*y(k,60)
         mat(k,2481) = -rxt(k,318)*y(k,60)
         mat(k,1909) = -rxt(k,319)*y(k,60)
         mat(k,1636) = -(rxt(k,794) + rxt(k,799) + rxt(k,804)) * y(k,60)
         mat(k,1973) = rxt(k,317)*y(k,125)
         mat(k,2017) = rxt(k,317)*y(k,59)
         mat(k,1234) = -((rxt(k,437) + rxt(k,448)) * y(k,250))
         mat(k,1928) = -(rxt(k,437) + rxt(k,448)) * y(k,62)
         mat(k,1013) = .230_r8*rxt(k,688)*y(k,135)
         mat(k,1656) = rxt(k,358)*y(k,42)
         mat(k,256) = .350_r8*rxt(k,454)*y(k,250)
         mat(k,499) = .630_r8*rxt(k,456)*y(k,135)
         mat(k,1157) = .560_r8*rxt(k,493)*y(k,135)
         mat(k,1772) = rxt(k,358)*y(k,17) + rxt(k,303)*y(k,56) + rxt(k,424)*y(k,126) &
                      + rxt(k,425)*y(k,133) + rxt(k,426)*y(k,250)
         mat(k,1416) = rxt(k,502)*y(k,126) + rxt(k,503)*y(k,250)
         mat(k,2068) = rxt(k,303)*y(k,42)
         mat(k,914) = rxt(k,484)*y(k,250)
         mat(k,968) = .620_r8*rxt(k,600)*y(k,135)
         mat(k,1401) = .650_r8*rxt(k,536)*y(k,135)
         mat(k,1045) = .230_r8*rxt(k,693)*y(k,135)
         mat(k,1528) = .560_r8*rxt(k,554)*y(k,135)
         mat(k,2423) = .170_r8*rxt(k,645)*y(k,227) + .220_r8*rxt(k,531)*y(k,236) &
                      + .400_r8*rxt(k,650)*y(k,237) + .350_r8*rxt(k,655)*y(k,239) &
                      + .225_r8*rxt(k,706)*y(k,255) + .250_r8*rxt(k,621)*y(k,260)
         mat(k,2330) = rxt(k,424)*y(k,42) + rxt(k,502)*y(k,49) + .220_r8*rxt(k,530) &
                      *y(k,236) + .500_r8*rxt(k,622)*y(k,260)
         mat(k,2486) = rxt(k,425)*y(k,42) + rxt(k,725)*y(k,137)
         mat(k,2232) = .230_r8*rxt(k,688)*y(k,6) + .630_r8*rxt(k,456)*y(k,25) &
                      + .560_r8*rxt(k,493)*y(k,29) + .620_r8*rxt(k,600)*y(k,98) &
                      + .650_r8*rxt(k,536)*y(k,105) + .230_r8*rxt(k,693)*y(k,110) &
                      + .560_r8*rxt(k,554)*y(k,111)
         mat(k,324) = rxt(k,725)*y(k,133) + rxt(k,726)*y(k,250)
         mat(k,1226) = .700_r8*rxt(k,718)*y(k,250)
         mat(k,1583) = .220_r8*rxt(k,527)*y(k,236) + .250_r8*rxt(k,618)*y(k,260)
         mat(k,1695) = .110_r8*rxt(k,528)*y(k,236) + .125_r8*rxt(k,704)*y(k,255) &
                      + .200_r8*rxt(k,619)*y(k,260)
         mat(k,830) = .170_r8*rxt(k,645)*y(k,124) + .070_r8*rxt(k,644)*y(k,232)
         mat(k,2161) = .070_r8*rxt(k,644)*y(k,227) + .160_r8*rxt(k,649)*y(k,237) &
                      + .140_r8*rxt(k,654)*y(k,239)
         mat(k,1506) = .220_r8*rxt(k,531)*y(k,124) + .220_r8*rxt(k,530)*y(k,126) &
                      + .220_r8*rxt(k,527)*y(k,225) + .110_r8*rxt(k,528)*y(k,226)
         mat(k,804) = .400_r8*rxt(k,650)*y(k,124) + .160_r8*rxt(k,649)*y(k,232)
         mat(k,986) = .350_r8*rxt(k,655)*y(k,124) + .140_r8*rxt(k,654)*y(k,232)
         mat(k,1928) = mat(k,1928) + .350_r8*rxt(k,454)*y(k,24) + rxt(k,426)*y(k,42) &
                      + rxt(k,503)*y(k,49) + rxt(k,484)*y(k,75) + rxt(k,726)*y(k,137) &
                      + .700_r8*rxt(k,718)*y(k,178)
         mat(k,1282) = .225_r8*rxt(k,706)*y(k,124) + .125_r8*rxt(k,704)*y(k,226)
         mat(k,1359) = .250_r8*rxt(k,621)*y(k,124) + .500_r8*rxt(k,622)*y(k,126) &
                      + .250_r8*rxt(k,618)*y(k,225) + .200_r8*rxt(k,619)*y(k,226)
      end do
      end subroutine nlnmat02
      subroutine nlnmat03( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,143) = -(rxt(k,397)*y(k,246))
         mat(k,2267) = -rxt(k,397)*y(k,64)
         mat(k,207) = -(rxt(k,398)*y(k,246))
         mat(k,2269) = -rxt(k,398)*y(k,65)
         mat(k,120) = -(rxt(k,643)*y(k,250))
         mat(k,1809) = -rxt(k,643)*y(k,66)
         mat(k,114) = .180_r8*rxt(k,674)*y(k,250)
         mat(k,1809) = mat(k,1809) + .180_r8*rxt(k,674)*y(k,180)
         mat(k,336) = -(rxt(k,742)*y(k,126) + (rxt(k,743) + rxt(k,746)) * y(k,250))
         mat(k,2305) = -rxt(k,742)*y(k,67)
         mat(k,1848) = -(rxt(k,743) + rxt(k,746)) * y(k,67)
         mat(k,880) = rxt(k,478)*y(k,232)
         mat(k,2090) = rxt(k,478)*y(k,231)
         mat(k,870) = -(rxt(k,393)*y(k,77) + rxt(k,394)*y(k,263) + rxt(k,395)*y(k,89))
         mat(k,1337) = -rxt(k,393)*y(k,73)
         mat(k,2590) = -rxt(k,394)*y(k,73)
         mat(k,2509) = -rxt(k,395)*y(k,73)
         mat(k,145) = 2.000_r8*rxt(k,397)*y(k,246)
         mat(k,209) = rxt(k,398)*y(k,246)
         mat(k,2274) = 2.000_r8*rxt(k,397)*y(k,64) + rxt(k,398)*y(k,65)
         mat(k,1190) = -(rxt(k,483)*y(k,250))
         mat(k,1924) = -rxt(k,483)*y(k,74)
         mat(k,581) = .700_r8*rxt(k,589)*y(k,250)
         mat(k,533) = .500_r8*rxt(k,590)*y(k,250)
         mat(k,346) = rxt(k,605)*y(k,250)
         mat(k,2419) = .050_r8*rxt(k,585)*y(k,235) + .530_r8*rxt(k,531)*y(k,236) &
                      + .225_r8*rxt(k,706)*y(k,255) + .250_r8*rxt(k,621)*y(k,260)
         mat(k,2326) = .050_r8*rxt(k,586)*y(k,235) + .530_r8*rxt(k,530)*y(k,236) &
                      + .250_r8*rxt(k,622)*y(k,260)
         mat(k,1746) = rxt(k,482)*y(k,230)
         mat(k,1582) = .530_r8*rxt(k,527)*y(k,236) + .250_r8*rxt(k,618)*y(k,260)
         mat(k,1692) = .260_r8*rxt(k,528)*y(k,236) + .125_r8*rxt(k,704)*y(k,255) &
                      + .100_r8*rxt(k,619)*y(k,260)
         mat(k,402) = rxt(k,482)*y(k,134)
         mat(k,1440) = .050_r8*rxt(k,585)*y(k,124) + .050_r8*rxt(k,586)*y(k,126)
         mat(k,1504) = .530_r8*rxt(k,531)*y(k,124) + .530_r8*rxt(k,530)*y(k,126) &
                      + .530_r8*rxt(k,527)*y(k,225) + .260_r8*rxt(k,528)*y(k,226)
         mat(k,1924) = mat(k,1924) + .700_r8*rxt(k,589)*y(k,99) + .500_r8*rxt(k,590) &
                      *y(k,100) + rxt(k,605)*y(k,115)
         mat(k,1280) = .225_r8*rxt(k,706)*y(k,124) + .125_r8*rxt(k,704)*y(k,226)
         mat(k,1358) = .250_r8*rxt(k,621)*y(k,124) + .250_r8*rxt(k,622)*y(k,126) &
                      + .250_r8*rxt(k,618)*y(k,225) + .100_r8*rxt(k,619)*y(k,226)
         mat(k,913) = -(rxt(k,484)*y(k,250))
         mat(k,1906) = -rxt(k,484)*y(k,75)
         mat(k,255) = .650_r8*rxt(k,454)*y(k,250)
         mat(k,1189) = .200_r8*rxt(k,483)*y(k,250)
         mat(k,1173) = rxt(k,606)*y(k,250)
         mat(k,2407) = rxt(k,631)*y(k,220) + .050_r8*rxt(k,585)*y(k,235) &
                      + .400_r8*rxt(k,650)*y(k,237) + .170_r8*rxt(k,655)*y(k,239) &
                      + .700_r8*rxt(k,660)*y(k,252) + .600_r8*rxt(k,671)*y(k,257) &
                      + .250_r8*rxt(k,621)*y(k,260) + .340_r8*rxt(k,679)*y(k,261) &
                      + .170_r8*rxt(k,683)*y(k,262)
         mat(k,2311) = .050_r8*rxt(k,586)*y(k,235) + .250_r8*rxt(k,622)*y(k,260)
         mat(k,593) = rxt(k,631)*y(k,124)
         mat(k,1576) = .250_r8*rxt(k,618)*y(k,260)
         mat(k,1683) = .100_r8*rxt(k,619)*y(k,260)
         mat(k,2148) = .160_r8*rxt(k,649)*y(k,237) + .070_r8*rxt(k,654)*y(k,239)
         mat(k,1434) = .050_r8*rxt(k,585)*y(k,124) + .050_r8*rxt(k,586)*y(k,126)
         mat(k,801) = .400_r8*rxt(k,650)*y(k,124) + .160_r8*rxt(k,649)*y(k,232)
         mat(k,982) = .170_r8*rxt(k,655)*y(k,124) + .070_r8*rxt(k,654)*y(k,232)
         mat(k,1906) = mat(k,1906) + .650_r8*rxt(k,454)*y(k,24) + .200_r8*rxt(k,483) &
                      *y(k,74) + rxt(k,606)*y(k,116)
         mat(k,541) = .700_r8*rxt(k,660)*y(k,124)
         mat(k,790) = .600_r8*rxt(k,671)*y(k,124)
         mat(k,1354) = .250_r8*rxt(k,621)*y(k,124) + .250_r8*rxt(k,622)*y(k,126) &
                      + .250_r8*rxt(k,618)*y(k,225) + .100_r8*rxt(k,619)*y(k,226)
         mat(k,816) = .340_r8*rxt(k,679)*y(k,124)
         mat(k,616) = .170_r8*rxt(k,683)*y(k,124)
         mat(k,2552) = -((rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,232) + rxt(k,217) &
                      *y(k,134) + rxt(k,222)*y(k,135))
         mat(k,2194) = -(rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,76)
         mat(k,1761) = -rxt(k,217)*y(k,76)
         mat(k,2263) = -rxt(k,222)*y(k,76)
         mat(k,1790) = rxt(k,426)*y(k,250)
         mat(k,2087) = rxt(k,305)*y(k,77)
         mat(k,1237) = rxt(k,448)*y(k,250)
         mat(k,878) = rxt(k,393)*y(k,77)
         mat(k,1350) = rxt(k,305)*y(k,56) + rxt(k,393)*y(k,73) + rxt(k,210)*y(k,133) &
                      + rxt(k,182)*y(k,246) + rxt(k,224)*y(k,250)
         mat(k,853) = rxt(k,385)*y(k,246)
         mat(k,1650) = rxt(k,343)*y(k,246)
         mat(k,940) = rxt(k,254)*y(k,250)
         mat(k,2505) = rxt(k,210)*y(k,77) + rxt(k,227)*y(k,250)
         mat(k,328) = rxt(k,726)*y(k,250)
         mat(k,717) = rxt(k,731)*y(k,250)
         mat(k,1626) = rxt(k,736)*y(k,250)
         mat(k,2295) = rxt(k,182)*y(k,77) + rxt(k,385)*y(k,81) + rxt(k,343)*y(k,85)
         mat(k,1962) = rxt(k,426)*y(k,42) + rxt(k,448)*y(k,62) + rxt(k,224)*y(k,77) &
                      + rxt(k,254)*y(k,112) + rxt(k,227)*y(k,133) + rxt(k,726) &
                      *y(k,137) + rxt(k,731)*y(k,148) + rxt(k,736)*y(k,150)
         mat(k,1341) = -(rxt(k,182)*y(k,246) + rxt(k,210)*y(k,133) + rxt(k,224) &
                      *y(k,250) + rxt(k,305)*y(k,56) + rxt(k,393)*y(k,73))
         mat(k,2278) = -rxt(k,182)*y(k,77)
         mat(k,2487) = -rxt(k,210)*y(k,77)
         mat(k,1934) = -rxt(k,224)*y(k,77)
         mat(k,2069) = -rxt(k,305)*y(k,77)
         mat(k,873) = -rxt(k,393)*y(k,77)
         mat(k,2535) = rxt(k,214)*y(k,232)
         mat(k,2167) = rxt(k,214)*y(k,76)
         mat(k,573) = -(rxt(k,211)*y(k,133) + rxt(k,225)*y(k,250) + rxt(k,306)*y(k,56))
         mat(k,2469) = -rxt(k,211)*y(k,79)
         mat(k,1878) = -rxt(k,225)*y(k,79)
         mat(k,2055) = -rxt(k,306)*y(k,79)
         mat(k,2126) = 2.000_r8*rxt(k,233)*y(k,232)
         mat(k,1878) = mat(k,1878) + 2.000_r8*rxt(k,230)*y(k,250)
         mat(k,212) = rxt(k,741)*y(k,263)
         mat(k,2585) = rxt(k,741)*y(k,152)
         mat(k,848) = -(rxt(k,376)*y(k,133) + rxt(k,377)*y(k,250) + (rxt(k,384) &
                      + rxt(k,385)) * y(k,246))
         mat(k,2478) = -rxt(k,376)*y(k,81)
         mat(k,1900) = -rxt(k,377)*y(k,81)
         mat(k,2273) = -(rxt(k,384) + rxt(k,385)) * y(k,81)
         mat(k,1655) = rxt(k,358)*y(k,42) + rxt(k,359)*y(k,232)
         mat(k,1767) = rxt(k,358)*y(k,17)
         mat(k,2143) = rxt(k,359)*y(k,17)
         mat(k,1640) = -(rxt(k,328)*y(k,133) + rxt(k,329)*y(k,250) + (rxt(k,342) &
                      + rxt(k,343)) * y(k,246) + (rxt(k,789) + rxt(k,795) + rxt(k,800) &
                      ) * y(k,92) + (rxt(k,794) + rxt(k,799) + rxt(k,804)) * y(k,60) &
                      + (rxt(k,796) + rxt(k,801)) * y(k,91))
         mat(k,2489) = -rxt(k,328)*y(k,85)
         mat(k,1946) = -rxt(k,329)*y(k,85)
         mat(k,2279) = -(rxt(k,342) + rxt(k,343)) * y(k,85)
         mat(k,840) = -(rxt(k,789) + rxt(k,795) + rxt(k,800)) * y(k,85)
         mat(k,945) = -(rxt(k,794) + rxt(k,799) + rxt(k,804)) * y(k,85)
         mat(k,731) = -(rxt(k,796) + rxt(k,801)) * y(k,85)
         mat(k,261) = rxt(k,465)*y(k,56)
         mat(k,1774) = rxt(k,303)*y(k,56)
         mat(k,2071) = rxt(k,465)*y(k,28) + rxt(k,303)*y(k,42) + rxt(k,305)*y(k,77) &
                      + rxt(k,306)*y(k,79) + rxt(k,331)*y(k,92) + rxt(k,307)*y(k,232)
         mat(k,1979) = rxt(k,323)*y(k,250)
         mat(k,1342) = rxt(k,305)*y(k,56)
         mat(k,574) = rxt(k,306)*y(k,56)
         mat(k,840) = mat(k,840) + rxt(k,331)*y(k,56)
         mat(k,2178) = rxt(k,307)*y(k,56)
         mat(k,1946) = mat(k,1946) + rxt(k,323)*y(k,59)
         mat(k,199) = -(rxt(k,438)*y(k,250) + rxt(k,447)*y(k,246))
         mat(k,1827) = -rxt(k,438)*y(k,86)
         mat(k,2268) = -rxt(k,447)*y(k,86)
         mat(k,775) = -(rxt(k,439)*y(k,250))
         mat(k,1893) = -rxt(k,439)*y(k,87)
         mat(k,1003) = .050_r8*rxt(k,688)*y(k,135)
         mat(k,254) = .350_r8*rxt(k,454)*y(k,250)
         mat(k,497) = .370_r8*rxt(k,456)*y(k,135)
         mat(k,1149) = .120_r8*rxt(k,493)*y(k,135)
         mat(k,962) = .110_r8*rxt(k,600)*y(k,135)
         mat(k,1400) = .330_r8*rxt(k,536)*y(k,135)
         mat(k,1035) = .050_r8*rxt(k,693)*y(k,135)
         mat(k,1525) = .120_r8*rxt(k,554)*y(k,135)
         mat(k,2398) = rxt(k,442)*y(k,233)
         mat(k,2213) = .050_r8*rxt(k,688)*y(k,6) + .370_r8*rxt(k,456)*y(k,25) &
                      + .120_r8*rxt(k,493)*y(k,29) + .110_r8*rxt(k,600)*y(k,98) &
                      + .330_r8*rxt(k,536)*y(k,105) + .050_r8*rxt(k,693)*y(k,110) &
                      + .120_r8*rxt(k,554)*y(k,111)
         mat(k,2137) = rxt(k,440)*y(k,233)
         mat(k,524) = rxt(k,442)*y(k,124) + rxt(k,440)*y(k,232)
         mat(k,1893) = mat(k,1893) + .350_r8*rxt(k,454)*y(k,24)
         mat(k,869) = rxt(k,393)*y(k,77) + rxt(k,395)*y(k,89) + rxt(k,394)*y(k,263)
         mat(k,1334) = rxt(k,393)*y(k,73)
         mat(k,2508) = rxt(k,395)*y(k,73)
         mat(k,2586) = rxt(k,394)*y(k,73)
         mat(k,2529) = -(rxt(k,274)*y(k,250) + rxt(k,395)*y(k,73))
         mat(k,1961) = -rxt(k,274)*y(k,89)
         mat(k,877) = -rxt(k,395)*y(k,89)
         mat(k,1789) = rxt(k,424)*y(k,126)
         mat(k,1205) = rxt(k,467)*y(k,126)
         mat(k,1422) = rxt(k,502)*y(k,126)
         mat(k,952) = (rxt(k,794)+rxt(k,799)+rxt(k,804))*y(k,85)
         mat(k,343) = rxt(k,742)*y(k,126)
         mat(k,1649) = (rxt(k,794)+rxt(k,799)+rxt(k,804))*y(k,60)
         mat(k,2044) = rxt(k,268)*y(k,250)
         mat(k,2361) = rxt(k,424)*y(k,42) + rxt(k,467)*y(k,45) + rxt(k,502)*y(k,49) &
                      + rxt(k,742)*y(k,67)
         mat(k,1961) = mat(k,1961) + rxt(k,268)*y(k,125)
         mat(k,416) = -(rxt(k,234)*y(k,250))
         mat(k,1859) = -rxt(k,234)*y(k,90)
         mat(k,2002) = rxt(k,266)*y(k,232)
         mat(k,2112) = rxt(k,266)*y(k,125)
         mat(k,730) = -(rxt(k,379)*y(k,133) + (rxt(k,796) + rxt(k,801)) * y(k,85))
         mat(k,2473) = -rxt(k,379)*y(k,91)
         mat(k,1634) = -(rxt(k,796) + rxt(k,801)) * y(k,91)
         mat(k,2560) = rxt(k,365)*y(k,232)
         mat(k,2135) = rxt(k,365)*y(k,19)
         mat(k,839) = -(rxt(k,331)*y(k,56) + rxt(k,332)*y(k,133) + rxt(k,333)*y(k,250) &
                      + (rxt(k,789) + rxt(k,795) + rxt(k,800)) * y(k,85))
         mat(k,2057) = -rxt(k,331)*y(k,92)
         mat(k,2477) = -rxt(k,332)*y(k,92)
         mat(k,1899) = -rxt(k,333)*y(k,92)
         mat(k,1635) = -(rxt(k,789) + rxt(k,795) + rxt(k,800)) * y(k,92)
         mat(k,1971) = rxt(k,314)*y(k,232)
         mat(k,943) = rxt(k,319)*y(k,250)
         mat(k,2142) = rxt(k,314)*y(k,59)
         mat(k,1899) = mat(k,1899) + rxt(k,319)*y(k,60)
         mat(k,1211) = -(rxt(k,526)*y(k,250))
         mat(k,1926) = -rxt(k,526)*y(k,93)
         mat(k,582) = .300_r8*rxt(k,589)*y(k,250)
         mat(k,534) = .500_r8*rxt(k,590)*y(k,250)
         mat(k,2421) = rxt(k,523)*y(k,229) + rxt(k,532)*y(k,236)
         mat(k,704) = rxt(k,523)*y(k,124)
         mat(k,1505) = rxt(k,532)*y(k,124)
         mat(k,1926) = mat(k,1926) + .300_r8*rxt(k,589)*y(k,99) + .500_r8*rxt(k,590) &
                      *y(k,100)
      end do
      end subroutine nlnmat03
      subroutine nlnmat04( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,172) = -(rxt(k,571)*y(k,250))
         mat(k,1821) = -rxt(k,571)*y(k,94)
         mat(k,1238) = -(rxt(k,505)*y(k,250))
         mat(k,1929) = -rxt(k,505)*y(k,95)
         mat(k,583) = .700_r8*rxt(k,589)*y(k,250)
         mat(k,535) = .500_r8*rxt(k,590)*y(k,250)
         mat(k,474) = .500_r8*rxt(k,553)*y(k,250)
         mat(k,2424) = .050_r8*rxt(k,585)*y(k,235) + .220_r8*rxt(k,531)*y(k,236) &
                      + .250_r8*rxt(k,621)*y(k,260)
         mat(k,2331) = .050_r8*rxt(k,586)*y(k,235) + .220_r8*rxt(k,530)*y(k,236) &
                      + .250_r8*rxt(k,622)*y(k,260)
         mat(k,484) = .500_r8*rxt(k,510)*y(k,250)
         mat(k,1584) = .220_r8*rxt(k,527)*y(k,236) + .250_r8*rxt(k,618)*y(k,260)
         mat(k,1696) = .230_r8*rxt(k,528)*y(k,236) + .200_r8*rxt(k,511)*y(k,254) &
                      + .100_r8*rxt(k,619)*y(k,260)
         mat(k,1443) = .050_r8*rxt(k,585)*y(k,124) + .050_r8*rxt(k,586)*y(k,126)
         mat(k,1507) = .220_r8*rxt(k,531)*y(k,124) + .220_r8*rxt(k,530)*y(k,126) &
                      + .220_r8*rxt(k,527)*y(k,225) + .230_r8*rxt(k,528)*y(k,226)
         mat(k,1929) = mat(k,1929) + .700_r8*rxt(k,589)*y(k,99) + .500_r8*rxt(k,590) &
                      *y(k,100) + .500_r8*rxt(k,553)*y(k,109) + .500_r8*rxt(k,510) &
                      *y(k,146)
         mat(k,1323) = .200_r8*rxt(k,511)*y(k,226)
         mat(k,1360) = .250_r8*rxt(k,621)*y(k,124) + .250_r8*rxt(k,622)*y(k,126) &
                      + .250_r8*rxt(k,618)*y(k,225) + .100_r8*rxt(k,619)*y(k,226)
         mat(k,310) = -(rxt(k,572)*y(k,250))
         mat(k,1844) = -rxt(k,572)*y(k,96)
         mat(k,2371) = .870_r8*rxt(k,585)*y(k,235)
         mat(k,2303) = .950_r8*rxt(k,586)*y(k,235)
         mat(k,1571) = rxt(k,581)*y(k,235)
         mat(k,1672) = .750_r8*rxt(k,582)*y(k,235)
         mat(k,1428) = .870_r8*rxt(k,585)*y(k,124) + .950_r8*rxt(k,586)*y(k,126) &
                      + rxt(k,581)*y(k,225) + .750_r8*rxt(k,582)*y(k,226)
         mat(k,69) = -(rxt(k,573)*y(k,250))
         mat(k,1804) = -rxt(k,573)*y(k,97)
         mat(k,674) = .600_r8*rxt(k,602)*y(k,250)
         mat(k,1804) = mat(k,1804) + .600_r8*rxt(k,602)*y(k,103)
         mat(k,963) = -(rxt(k,591)*y(k,126) + rxt(k,600)*y(k,135) + rxt(k,601) &
                      *y(k,250))
         mat(k,2313) = -rxt(k,591)*y(k,98)
         mat(k,2217) = -rxt(k,600)*y(k,98)
         mat(k,1910) = -rxt(k,601)*y(k,98)
         mat(k,580) = -(rxt(k,589)*y(k,250))
         mat(k,1879) = -rxt(k,589)*y(k,99)
         mat(k,2385) = .080_r8*rxt(k,577)*y(k,234)
         mat(k,1466) = .080_r8*rxt(k,577)*y(k,124)
         mat(k,531) = -(rxt(k,590)*y(k,250))
         mat(k,1873) = -rxt(k,590)*y(k,100)
         mat(k,2382) = .080_r8*rxt(k,585)*y(k,235)
         mat(k,1429) = .080_r8*rxt(k,585)*y(k,124)
         mat(k,1384) = -(rxt(k,592)*y(k,225) + rxt(k,593)*y(k,226) + rxt(k,594) &
                      *y(k,232) + rxt(k,595)*y(k,124) + rxt(k,596)*y(k,126))
         mat(k,1586) = -rxt(k,592)*y(k,101)
         mat(k,1702) = -rxt(k,593)*y(k,101)
         mat(k,2169) = -rxt(k,594)*y(k,101)
         mat(k,2430) = -rxt(k,595)*y(k,101)
         mat(k,2337) = -rxt(k,596)*y(k,101)
         mat(k,969) = rxt(k,591)*y(k,126)
         mat(k,2337) = mat(k,2337) + rxt(k,591)*y(k,98)
         mat(k,368) = -(rxt(k,599)*y(k,250))
         mat(k,1853) = -rxt(k,599)*y(k,102)
         mat(k,1373) = rxt(k,594)*y(k,232)
         mat(k,2106) = rxt(k,594)*y(k,101)
         mat(k,675) = -(rxt(k,602)*y(k,250))
         mat(k,1887) = -rxt(k,602)*y(k,103)
         mat(k,2132) = rxt(k,576)*y(k,234) + rxt(k,583)*y(k,235)
         mat(k,1467) = rxt(k,576)*y(k,232)
         mat(k,1431) = rxt(k,583)*y(k,232)
         mat(k,39) = -(rxt(k,781)*y(k,250))
         mat(k,1798) = -rxt(k,781)*y(k,104)
         mat(k,1402) = -(rxt(k,536)*y(k,135) + rxt(k,537)*y(k,250))
         mat(k,2239) = -rxt(k,536)*y(k,105)
         mat(k,1937) = -rxt(k,537)*y(k,105)
         mat(k,970) = .300_r8*rxt(k,600)*y(k,135)
         mat(k,2431) = .360_r8*rxt(k,577)*y(k,234)
         mat(k,2338) = .400_r8*rxt(k,578)*y(k,234)
         mat(k,2239) = mat(k,2239) + .300_r8*rxt(k,600)*y(k,98)
         mat(k,1587) = .390_r8*rxt(k,574)*y(k,234)
         mat(k,1703) = .310_r8*rxt(k,575)*y(k,234)
         mat(k,1479) = .360_r8*rxt(k,577)*y(k,124) + .400_r8*rxt(k,578)*y(k,126) &
                      + .390_r8*rxt(k,574)*y(k,225) + .310_r8*rxt(k,575)*y(k,226)
         mat(k,301) = -(rxt(k,538)*y(k,250))
         mat(k,1842) = -rxt(k,538)*y(k,106)
         mat(k,2100) = rxt(k,529)*y(k,236)
         mat(k,1498) = rxt(k,529)*y(k,232)
         mat(k,459) = -(rxt(k,551)*y(k,250))
         mat(k,1865) = -rxt(k,551)*y(k,107)
         mat(k,2377) = .800_r8*rxt(k,565)*y(k,219)
         mat(k,1058) = .800_r8*rxt(k,565)*y(k,124)
         mat(k,316) = -(rxt(k,552)*y(k,250))
         mat(k,1846) = -rxt(k,552)*y(k,108)
         mat(k,2101) = .800_r8*rxt(k,548)*y(k,240)
         mat(k,718) = .800_r8*rxt(k,548)*y(k,232)
         mat(k,473) = -(rxt(k,553)*y(k,250))
         mat(k,1867) = -rxt(k,553)*y(k,109)
         mat(k,2004) = rxt(k,557)*y(k,238)
         mat(k,1550) = rxt(k,557)*y(k,125)
         mat(k,1037) = -(rxt(k,692)*y(k,126) + rxt(k,693)*y(k,135) + rxt(k,694) &
                      *y(k,250))
         mat(k,2315) = -rxt(k,692)*y(k,110)
         mat(k,2219) = -rxt(k,693)*y(k,110)
         mat(k,1913) = -rxt(k,694)*y(k,110)
         mat(k,1532) = -(rxt(k,554)*y(k,135) + rxt(k,555)*y(k,250))
         mat(k,2244) = -rxt(k,554)*y(k,111)
         mat(k,1942) = -rxt(k,555)*y(k,111)
         mat(k,973) = .200_r8*rxt(k,600)*y(k,135)
         mat(k,2436) = .560_r8*rxt(k,577)*y(k,234)
         mat(k,2343) = .600_r8*rxt(k,578)*y(k,234)
         mat(k,2244) = mat(k,2244) + .200_r8*rxt(k,600)*y(k,98)
         mat(k,1592) = .610_r8*rxt(k,574)*y(k,234)
         mat(k,1708) = .440_r8*rxt(k,575)*y(k,234)
         mat(k,1483) = .560_r8*rxt(k,577)*y(k,124) + .600_r8*rxt(k,578)*y(k,126) &
                      + .610_r8*rxt(k,574)*y(k,225) + .440_r8*rxt(k,575)*y(k,226)
         mat(k,930) = -(rxt(k,237)*y(k,124) + (rxt(k,238) + rxt(k,239) + rxt(k,240) &
                      ) * y(k,125) + rxt(k,241)*y(k,134) + rxt(k,254)*y(k,250) &
                      + rxt(k,832)*y(k,249))
         mat(k,2408) = -rxt(k,237)*y(k,112)
         mat(k,2016) = -(rxt(k,238) + rxt(k,239) + rxt(k,240)) * y(k,112)
         mat(k,1742) = -rxt(k,241)*y(k,112)
         mat(k,1908) = -rxt(k,254)*y(k,112)
         mat(k,751) = -rxt(k,832)*y(k,112)
         mat(k,2480) = rxt(k,235)*y(k,241) + rxt(k,829)*y(k,244)
         mat(k,1742) = mat(k,1742) + rxt(k,830)*y(k,244)
         mat(k,769) = 1.100_r8*rxt(k,825)*y(k,242) + .200_r8*rxt(k,823)*y(k,243)
         mat(k,465) = rxt(k,235)*y(k,133)
         mat(k,647) = 1.100_r8*rxt(k,825)*y(k,228)
         mat(k,759) = .200_r8*rxt(k,823)*y(k,228)
         mat(k,450) = rxt(k,829)*y(k,133) + rxt(k,830)*y(k,134)
         mat(k,2001) = rxt(k,267)*y(k,126)
         mat(k,2301) = rxt(k,267)*y(k,125)
         mat(k,344) = -(rxt(k,605)*y(k,250))
         mat(k,1849) = -rxt(k,605)*y(k,115)
         mat(k,1372) = .200_r8*rxt(k,593)*y(k,226)
         mat(k,1675) = .200_r8*rxt(k,593)*y(k,101)
         mat(k,1178) = -(rxt(k,606)*y(k,250))
         mat(k,1923) = -rxt(k,606)*y(k,116)
         mat(k,1380) = rxt(k,595)*y(k,124) + rxt(k,596)*y(k,126) + rxt(k,592)*y(k,225) &
                      + .800_r8*rxt(k,593)*y(k,226)
         mat(k,2418) = rxt(k,595)*y(k,101)
         mat(k,2325) = rxt(k,596)*y(k,101)
         mat(k,1581) = rxt(k,592)*y(k,101)
         mat(k,1691) = .800_r8*rxt(k,593)*y(k,101)
         mat(k,49) = -(rxt(k,745)*y(k,250))
         mat(k,1800) = -rxt(k,745)*y(k,120)
         mat(k,2451) = -(rxt(k,237)*y(k,112) + rxt(k,249)*y(k,126) + rxt(k,255) &
                      *y(k,232) + rxt(k,256)*y(k,135) + rxt(k,257)*y(k,133) + rxt(k,315) &
                      *y(k,59) + rxt(k,366)*y(k,19) + rxt(k,432)*y(k,226) + rxt(k,442) &
                      *y(k,233) + rxt(k,461)*y(k,222) + rxt(k,474)*y(k,225) + rxt(k,479) &
                      *y(k,231) + rxt(k,498)*y(k,223) + rxt(k,508)*y(k,253) + rxt(k,513) &
                      *y(k,254) + (rxt(k,522) + rxt(k,523)) * y(k,229) + (rxt(k,531) &
                      + rxt(k,532)) * y(k,236) + rxt(k,544)*y(k,238) + rxt(k,549) &
                      *y(k,240) + (rxt(k,565) + rxt(k,566)) * y(k,219) + rxt(k,577) &
                      *y(k,234) + rxt(k,585)*y(k,235) + rxt(k,595)*y(k,101) + rxt(k,621) &
                      *y(k,260) + rxt(k,627)*y(k,218) + rxt(k,631)*y(k,220) + rxt(k,637) &
                      *y(k,221) + rxt(k,640)*y(k,224) + rxt(k,645)*y(k,227) + rxt(k,650) &
                      *y(k,237) + rxt(k,655)*y(k,239) + rxt(k,660)*y(k,252) + rxt(k,671) &
                      *y(k,257) + rxt(k,679)*y(k,261) + rxt(k,683)*y(k,262) + rxt(k,699) &
                      *y(k,245) + rxt(k,706)*y(k,255) + rxt(k,712)*y(k,256) + rxt(k,834) &
                      *y(k,249))
         mat(k,938) = -rxt(k,237)*y(k,124)
         mat(k,2359) = -rxt(k,249)*y(k,124)
         mat(k,2191) = -rxt(k,255)*y(k,124)
         mat(k,2260) = -rxt(k,256)*y(k,124)
         mat(k,2502) = -rxt(k,257)*y(k,124)
         mat(k,1992) = -rxt(k,315)*y(k,124)
         mat(k,2579) = -rxt(k,366)*y(k,124)
         mat(k,1721) = -rxt(k,432)*y(k,124)
         mat(k,530) = -rxt(k,442)*y(k,124)
         mat(k,1132) = -rxt(k,461)*y(k,124)
         mat(k,1603) = -rxt(k,474)*y(k,124)
         mat(k,891) = -rxt(k,479)*y(k,124)
         mat(k,910) = -rxt(k,498)*y(k,124)
         mat(k,867) = -rxt(k,508)*y(k,124)
         mat(k,1332) = -rxt(k,513)*y(k,124)
         mat(k,708) = -(rxt(k,522) + rxt(k,523)) * y(k,124)
         mat(k,1519) = -(rxt(k,531) + rxt(k,532)) * y(k,124)
         mat(k,1565) = -rxt(k,544)*y(k,124)
         mat(k,727) = -rxt(k,549)*y(k,124)
         mat(k,1073) = -(rxt(k,565) + rxt(k,566)) * y(k,124)
         mat(k,1493) = -rxt(k,577)*y(k,124)
         mat(k,1460) = -rxt(k,585)*y(k,124)
         mat(k,1395) = -rxt(k,595)*y(k,124)
         mat(k,1370) = -rxt(k,621)*y(k,124)
         mat(k,673) = -rxt(k,627)*y(k,124)
         mat(k,599) = -rxt(k,631)*y(k,124)
         mat(k,522) = -rxt(k,637)*y(k,124)
         mat(k,747) = -rxt(k,640)*y(k,124)
         mat(k,836) = -rxt(k,645)*y(k,124)
         mat(k,807) = -rxt(k,650)*y(k,124)
         mat(k,992) = -rxt(k,655)*y(k,124)
         mat(k,547) = -rxt(k,660)*y(k,124)
         mat(k,799) = -rxt(k,671)*y(k,124)
         mat(k,826) = -rxt(k,679)*y(k,124)
         mat(k,623) = -rxt(k,683)*y(k,124)
         mat(k,1314) = -rxt(k,699)*y(k,124)
         mat(k,1291) = -rxt(k,706)*y(k,124)
         mat(k,1267) = -rxt(k,712)*y(k,124)
         mat(k,754) = -rxt(k,834)*y(k,124)
         mat(k,938) = mat(k,938) + 2.000_r8*rxt(k,239)*y(k,125) + rxt(k,241)*y(k,134) &
                      + rxt(k,254)*y(k,250)
         mat(k,2042) = 2.000_r8*rxt(k,239)*y(k,112) + rxt(k,242)*y(k,133) + rxt(k,732) &
                      *y(k,150)
         mat(k,2502) = mat(k,2502) + rxt(k,242)*y(k,125)
         mat(k,1759) = rxt(k,241)*y(k,112) + rxt(k,236)*y(k,241)
         mat(k,1624) = rxt(k,732)*y(k,125)
         mat(k,468) = rxt(k,236)*y(k,134)
         mat(k,1959) = rxt(k,254)*y(k,112)
      end do
      end subroutine nlnmat04
      subroutine nlnmat05( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,2036) = -((rxt(k,238) + rxt(k,239) + rxt(k,240)) * y(k,112) + (rxt(k,242) &
                      + rxt(k,244)) * y(k,133) + rxt(k,243)*y(k,135) + rxt(k,266) &
                      *y(k,232) + rxt(k,267)*y(k,126) + rxt(k,268)*y(k,250) + rxt(k,317) &
                      *y(k,59) + rxt(k,368)*y(k,19) + rxt(k,487)*y(k,225) + rxt(k,557) &
                      *y(k,238) + rxt(k,646)*y(k,227) + rxt(k,651)*y(k,237) + rxt(k,656) &
                      *y(k,239) + rxt(k,663)*y(k,141) + rxt(k,668)*y(k,218) + rxt(k,732) &
                      *y(k,150))
         mat(k,936) = -(rxt(k,238) + rxt(k,239) + rxt(k,240)) * y(k,125)
         mat(k,2496) = -(rxt(k,242) + rxt(k,244)) * y(k,125)
         mat(k,2254) = -rxt(k,243)*y(k,125)
         mat(k,2185) = -rxt(k,266)*y(k,125)
         mat(k,2353) = -rxt(k,267)*y(k,125)
         mat(k,1953) = -rxt(k,268)*y(k,125)
         mat(k,1986) = -rxt(k,317)*y(k,125)
         mat(k,2573) = -rxt(k,368)*y(k,125)
         mat(k,1599) = -rxt(k,487)*y(k,125)
         mat(k,1561) = -rxt(k,557)*y(k,125)
         mat(k,834) = -rxt(k,646)*y(k,125)
         mat(k,805) = -rxt(k,651)*y(k,125)
         mat(k,990) = -rxt(k,656)*y(k,125)
         mat(k,605) = -rxt(k,663)*y(k,125)
         mat(k,671) = -rxt(k,668)*y(k,125)
         mat(k,1620) = -rxt(k,732)*y(k,125)
         mat(k,639) = rxt(k,563)*y(k,250)
         mat(k,414) = rxt(k,519)*y(k,126)
         mat(k,2573) = mat(k,2573) + rxt(k,366)*y(k,124)
         mat(k,1986) = mat(k,1986) + rxt(k,315)*y(k,124)
         mat(k,419) = rxt(k,234)*y(k,250)
         mat(k,586) = .700_r8*rxt(k,589)*y(k,250)
         mat(k,1392) = rxt(k,595)*y(k,124) + rxt(k,596)*y(k,126)
         mat(k,2445) = rxt(k,366)*y(k,19) + rxt(k,315)*y(k,59) + rxt(k,595)*y(k,101) &
                      + 2.000_r8*rxt(k,249)*y(k,126) + rxt(k,257)*y(k,133) &
                      + rxt(k,256)*y(k,135) + rxt(k,627)*y(k,218) + rxt(k,565) &
                      *y(k,219) + rxt(k,631)*y(k,220) + rxt(k,637)*y(k,221) &
                      + rxt(k,461)*y(k,222) + rxt(k,498)*y(k,223) + rxt(k,640) &
                      *y(k,224) + rxt(k,474)*y(k,225) + rxt(k,432)*y(k,226) &
                      + rxt(k,645)*y(k,227) + rxt(k,522)*y(k,229) + rxt(k,479) &
                      *y(k,231) + rxt(k,255)*y(k,232) + rxt(k,442)*y(k,233) &
                      + .920_r8*rxt(k,577)*y(k,234) + .920_r8*rxt(k,585)*y(k,235) &
                      + rxt(k,531)*y(k,236) + rxt(k,650)*y(k,237) + rxt(k,544) &
                      *y(k,238) + rxt(k,655)*y(k,239) + rxt(k,549)*y(k,240) &
                      + 1.600_r8*rxt(k,699)*y(k,245) + rxt(k,660)*y(k,252) &
                      + rxt(k,508)*y(k,253) + rxt(k,513)*y(k,254) + .900_r8*rxt(k,706) &
                      *y(k,255) + .800_r8*rxt(k,712)*y(k,256) + rxt(k,671)*y(k,257) &
                      + rxt(k,621)*y(k,260) + rxt(k,679)*y(k,261) + rxt(k,683) &
                      *y(k,262)
         mat(k,2353) = mat(k,2353) + rxt(k,519)*y(k,16) + rxt(k,596)*y(k,101) &
                      + 2.000_r8*rxt(k,249)*y(k,124) + rxt(k,250)*y(k,133) &
                      + rxt(k,248)*y(k,232) + rxt(k,578)*y(k,234) + rxt(k,586) &
                      *y(k,235) + rxt(k,530)*y(k,236) + rxt(k,545)*y(k,238) &
                      + 2.000_r8*rxt(k,700)*y(k,245) + rxt(k,251)*y(k,250) &
                      + rxt(k,622)*y(k,260)
         mat(k,897) = rxt(k,506)*y(k,250)
         mat(k,2496) = mat(k,2496) + rxt(k,257)*y(k,124) + rxt(k,250)*y(k,126)
         mat(k,2254) = mat(k,2254) + rxt(k,256)*y(k,124)
         mat(k,627) = rxt(k,709)*y(k,250)
         mat(k,671) = mat(k,671) + rxt(k,627)*y(k,124)
         mat(k,1071) = rxt(k,565)*y(k,124)
         mat(k,597) = rxt(k,631)*y(k,124)
         mat(k,520) = rxt(k,637)*y(k,124)
         mat(k,1130) = rxt(k,461)*y(k,124)
         mat(k,908) = rxt(k,498)*y(k,124)
         mat(k,744) = rxt(k,640)*y(k,124)
         mat(k,1599) = mat(k,1599) + rxt(k,474)*y(k,124)
         mat(k,1716) = rxt(k,432)*y(k,124) + .500_r8*rxt(k,697)*y(k,245)
         mat(k,834) = mat(k,834) + rxt(k,645)*y(k,124)
         mat(k,706) = rxt(k,522)*y(k,124)
         mat(k,889) = rxt(k,479)*y(k,124)
         mat(k,2185) = mat(k,2185) + rxt(k,255)*y(k,124) + rxt(k,248)*y(k,126)
         mat(k,528) = rxt(k,442)*y(k,124)
         mat(k,1489) = .920_r8*rxt(k,577)*y(k,124) + rxt(k,578)*y(k,126)
         mat(k,1456) = .920_r8*rxt(k,585)*y(k,124) + rxt(k,586)*y(k,126)
         mat(k,1516) = rxt(k,531)*y(k,124) + rxt(k,530)*y(k,126)
         mat(k,805) = mat(k,805) + rxt(k,650)*y(k,124)
         mat(k,1561) = mat(k,1561) + rxt(k,544)*y(k,124) + rxt(k,545)*y(k,126)
         mat(k,990) = mat(k,990) + rxt(k,655)*y(k,124)
         mat(k,725) = rxt(k,549)*y(k,124)
         mat(k,1311) = 1.600_r8*rxt(k,699)*y(k,124) + 2.000_r8*rxt(k,700)*y(k,126) &
                      + .500_r8*rxt(k,697)*y(k,226)
         mat(k,1953) = mat(k,1953) + rxt(k,563)*y(k,1) + rxt(k,234)*y(k,90) &
                      + .700_r8*rxt(k,589)*y(k,99) + rxt(k,251)*y(k,126) + rxt(k,506) &
                      *y(k,127) + rxt(k,709)*y(k,175)
         mat(k,545) = rxt(k,660)*y(k,124)
         mat(k,865) = rxt(k,508)*y(k,124)
         mat(k,1330) = rxt(k,513)*y(k,124)
         mat(k,1289) = .900_r8*rxt(k,706)*y(k,124)
         mat(k,1264) = .800_r8*rxt(k,712)*y(k,124)
         mat(k,797) = rxt(k,671)*y(k,124)
         mat(k,1367) = rxt(k,621)*y(k,124) + rxt(k,622)*y(k,126)
         mat(k,824) = rxt(k,679)*y(k,124)
         mat(k,621) = rxt(k,683)*y(k,124)
         mat(k,2358) = -(rxt(k,248)*y(k,232) + rxt(k,249)*y(k,124) + rxt(k,250) &
                      *y(k,133) + rxt(k,251)*y(k,250) + rxt(k,267)*y(k,125) + rxt(k,424) &
                      *y(k,42) + rxt(k,467)*y(k,45) + rxt(k,492)*y(k,29) + rxt(k,502) &
                      *y(k,49) + rxt(k,519)*y(k,16) + rxt(k,530)*y(k,236) + rxt(k,545) &
                      *y(k,238) + rxt(k,578)*y(k,234) + rxt(k,586)*y(k,235) + rxt(k,591) &
                      *y(k,98) + rxt(k,596)*y(k,101) + rxt(k,622)*y(k,260) + rxt(k,687) &
                      *y(k,6) + rxt(k,692)*y(k,110) + rxt(k,700)*y(k,245) + rxt(k,715) &
                      *y(k,177) + rxt(k,742)*y(k,67))
         mat(k,2190) = -rxt(k,248)*y(k,126)
         mat(k,2450) = -rxt(k,249)*y(k,126)
         mat(k,2501) = -rxt(k,250)*y(k,126)
         mat(k,1958) = -rxt(k,251)*y(k,126)
         mat(k,2041) = -rxt(k,267)*y(k,126)
         mat(k,1786) = -rxt(k,424)*y(k,126)
         mat(k,1204) = -rxt(k,467)*y(k,126)
         mat(k,1168) = -rxt(k,492)*y(k,126)
         mat(k,1421) = -rxt(k,502)*y(k,126)
         mat(k,415) = -rxt(k,519)*y(k,126)
         mat(k,1518) = -rxt(k,530)*y(k,126)
         mat(k,1564) = -rxt(k,545)*y(k,126)
         mat(k,1492) = -rxt(k,578)*y(k,126)
         mat(k,1459) = -rxt(k,586)*y(k,126)
         mat(k,980) = -rxt(k,591)*y(k,126)
         mat(k,1394) = -rxt(k,596)*y(k,126)
         mat(k,1369) = -rxt(k,622)*y(k,126)
         mat(k,1023) = -rxt(k,687)*y(k,126)
         mat(k,1055) = -rxt(k,692)*y(k,126)
         mat(k,1313) = -rxt(k,700)*y(k,126)
         mat(k,1144) = -rxt(k,715)*y(k,126)
         mat(k,342) = -rxt(k,742)*y(k,126)
         mat(k,509) = rxt(k,369)*y(k,133)
         mat(k,2083) = rxt(k,316)*y(k,60)
         mat(k,950) = rxt(k,316)*y(k,56) + rxt(k,318)*y(k,133) + rxt(k,319)*y(k,250)
         mat(k,876) = rxt(k,395)*y(k,89)
         mat(k,2526) = rxt(k,395)*y(k,73) + rxt(k,274)*y(k,250)
         mat(k,480) = .500_r8*rxt(k,553)*y(k,250)
         mat(k,2041) = mat(k,2041) + rxt(k,244)*y(k,133) + rxt(k,243)*y(k,135)
         mat(k,2501) = mat(k,2501) + rxt(k,369)*y(k,20) + rxt(k,318)*y(k,60) &
                      + rxt(k,244)*y(k,125)
         mat(k,2259) = rxt(k,243)*y(k,125)
         mat(k,445) = rxt(k,485)*y(k,250)
         mat(k,1958) = mat(k,1958) + rxt(k,319)*y(k,60) + rxt(k,274)*y(k,89) &
                      + .500_r8*rxt(k,553)*y(k,109) + rxt(k,485)*y(k,139)
         mat(k,892) = -(rxt(k,506)*y(k,250))
         mat(k,1904) = -rxt(k,506)*y(k,127)
         mat(k,1151) = rxt(k,492)*y(k,126)
         mat(k,532) = .500_r8*rxt(k,590)*y(k,250)
         mat(k,370) = rxt(k,599)*y(k,250)
         mat(k,345) = rxt(k,605)*y(k,250)
         mat(k,1172) = rxt(k,606)*y(k,250)
         mat(k,2310) = rxt(k,492)*y(k,29)
         mat(k,1904) = mat(k,1904) + .500_r8*rxt(k,590)*y(k,100) + rxt(k,599)*y(k,102) &
                      + rxt(k,605)*y(k,115) + rxt(k,606)*y(k,116)
         mat(k,350) = -(rxt(k,703)*y(k,250))
         mat(k,1850) = -rxt(k,703)*y(k,128)
         mat(k,2103) = rxt(k,698)*y(k,245)
         mat(k,1295) = rxt(k,698)*y(k,232)
         mat(k,2503) = -(rxt(k,197)*y(k,135) + 4._r8*rxt(k,199)*y(k,133) + rxt(k,200) &
                      *y(k,134) + rxt(k,210)*y(k,77) + rxt(k,211)*y(k,79) + rxt(k,218) &
                      *y(k,232) + rxt(k,227)*y(k,250) + (rxt(k,242) + rxt(k,244) &
                      ) * y(k,125) + rxt(k,250)*y(k,126) + rxt(k,257)*y(k,124) &
                      + rxt(k,318)*y(k,60) + rxt(k,321)*y(k,59) + rxt(k,328)*y(k,85) &
                      + rxt(k,332)*y(k,92) + rxt(k,369)*y(k,20) + rxt(k,371)*y(k,19) &
                      + rxt(k,376)*y(k,81) + rxt(k,379)*y(k,91) + rxt(k,425)*y(k,42) &
                      + rxt(k,725)*y(k,137) + (rxt(k,827) + rxt(k,828)) * y(k,242) &
                      + rxt(k,829)*y(k,244))
         mat(k,2261) = -rxt(k,197)*y(k,133)
         mat(k,1760) = -rxt(k,200)*y(k,133)
         mat(k,1348) = -rxt(k,210)*y(k,133)
         mat(k,578) = -rxt(k,211)*y(k,133)
         mat(k,2192) = -rxt(k,218)*y(k,133)
         mat(k,1960) = -rxt(k,227)*y(k,133)
         mat(k,2043) = -(rxt(k,242) + rxt(k,244)) * y(k,133)
         mat(k,2360) = -rxt(k,250)*y(k,133)
         mat(k,2452) = -rxt(k,257)*y(k,133)
         mat(k,951) = -rxt(k,318)*y(k,133)
         mat(k,1993) = -rxt(k,321)*y(k,133)
         mat(k,1648) = -rxt(k,328)*y(k,133)
         mat(k,844) = -rxt(k,332)*y(k,133)
         mat(k,510) = -rxt(k,369)*y(k,133)
         mat(k,2580) = -rxt(k,371)*y(k,133)
         mat(k,852) = -rxt(k,376)*y(k,133)
         mat(k,735) = -rxt(k,379)*y(k,133)
         mat(k,1788) = -rxt(k,425)*y(k,133)
         mat(k,327) = -rxt(k,725)*y(k,133)
         mat(k,652) = -(rxt(k,827) + rxt(k,828)) * y(k,133)
         mat(k,452) = -rxt(k,829)*y(k,133)
         mat(k,2550) = rxt(k,216)*y(k,232)
         mat(k,939) = rxt(k,237)*y(k,124) + rxt(k,238)*y(k,125) + rxt(k,241)*y(k,134) &
                      + rxt(k,832)*y(k,249)
         mat(k,2452) = mat(k,2452) + rxt(k,237)*y(k,112)
         mat(k,2043) = mat(k,2043) + rxt(k,238)*y(k,112)
         mat(k,1760) = mat(k,1760) + rxt(k,241)*y(k,112) + rxt(k,727)*y(k,148) &
                      + rxt(k,733)*y(k,150) + rxt(k,831)*y(k,244) + (rxt(k,185) &
                       +rxt(k,186))*y(k,246) + rxt(k,837)*y(k,251)
         mat(k,716) = rxt(k,727)*y(k,134)
         mat(k,1625) = rxt(k,733)*y(k,134)
         mat(k,774) = rxt(k,823)*y(k,243) + 1.150_r8*rxt(k,824)*y(k,249)
         mat(k,2192) = mat(k,2192) + rxt(k,216)*y(k,76)
         mat(k,763) = rxt(k,823)*y(k,228)
         mat(k,452) = mat(k,452) + rxt(k,831)*y(k,134)
         mat(k,2293) = (rxt(k,185)+rxt(k,186))*y(k,134)
         mat(k,755) = rxt(k,832)*y(k,112) + 1.150_r8*rxt(k,824)*y(k,228)
         mat(k,1960) = mat(k,1960) + 2.000_r8*rxt(k,229)*y(k,250)
         mat(k,613) = rxt(k,837)*y(k,134)
         mat(k,1750) = -(rxt(k,185)*y(k,246) + rxt(k,191)*y(k,247) + rxt(k,200) &
                      *y(k,133) + rxt(k,217)*y(k,76) + rxt(k,236)*y(k,241) + rxt(k,241) &
                      *y(k,112) + rxt(k,482)*y(k,230) + rxt(k,727)*y(k,148) + rxt(k,733) &
                      *y(k,150) + rxt(k,826)*y(k,242) + (rxt(k,830) + rxt(k,831) &
                      ) * y(k,244) + rxt(k,837)*y(k,251))
         mat(k,2282) = -rxt(k,185)*y(k,134)
         mat(k,79) = -rxt(k,191)*y(k,134)
         mat(k,2492) = -rxt(k,200)*y(k,134)
         mat(k,2539) = -rxt(k,217)*y(k,134)
         mat(k,466) = -rxt(k,236)*y(k,134)
         mat(k,934) = -rxt(k,241)*y(k,134)
         mat(k,403) = -rxt(k,482)*y(k,134)
         mat(k,713) = -rxt(k,727)*y(k,134)
         mat(k,1617) = -rxt(k,733)*y(k,134)
         mat(k,649) = -rxt(k,826)*y(k,134)
         mat(k,451) = -(rxt(k,830) + rxt(k,831)) * y(k,134)
         mat(k,612) = -rxt(k,837)*y(k,134)
         mat(k,1658) = rxt(k,360)*y(k,135) + rxt(k,359)*y(k,232)
         mat(k,2569) = 2.000_r8*rxt(k,361)*y(k,19) + (rxt(k,363)+rxt(k,364))*y(k,59) &
                      + rxt(k,371)*y(k,133) + rxt(k,365)*y(k,232)
         mat(k,2074) = rxt(k,309)*y(k,135) + rxt(k,307)*y(k,232)
         mat(k,1982) = (rxt(k,363)+rxt(k,364))*y(k,19) + (2.000_r8*rxt(k,311) &
                       +2.000_r8*rxt(k,312))*y(k,59) + rxt(k,321)*y(k,133) &
                      + rxt(k,314)*y(k,232) + rxt(k,323)*y(k,250)
         mat(k,2539) = mat(k,2539) + rxt(k,222)*y(k,135) + rxt(k,214)*y(k,232)
         mat(k,417) = rxt(k,234)*y(k,250)
         mat(k,934) = mat(k,934) + rxt(k,240)*y(k,125)
         mat(k,2441) = rxt(k,256)*y(k,135) + rxt(k,834)*y(k,249)
         mat(k,2032) = rxt(k,240)*y(k,112) + rxt(k,242)*y(k,133) + rxt(k,243)*y(k,135)
         mat(k,2349) = rxt(k,250)*y(k,133) + rxt(k,248)*y(k,232)
         mat(k,2492) = mat(k,2492) + rxt(k,371)*y(k,19) + rxt(k,321)*y(k,59) &
                      + rxt(k,242)*y(k,125) + rxt(k,250)*y(k,126) &
                      + 2.000_r8*rxt(k,199)*y(k,133) + 2.000_r8*rxt(k,197)*y(k,135) &
                      + rxt(k,218)*y(k,232) + rxt(k,190)*y(k,247) + rxt(k,227) &
                      *y(k,250)
         mat(k,1750) = mat(k,1750) + 2.000_r8*rxt(k,191)*y(k,247)
         mat(k,2250) = rxt(k,360)*y(k,17) + rxt(k,309)*y(k,56) + rxt(k,222)*y(k,76) &
                      + rxt(k,256)*y(k,124) + rxt(k,243)*y(k,125) &
                      + 2.000_r8*rxt(k,197)*y(k,133) + rxt(k,728)*y(k,148) &
                      + rxt(k,734)*y(k,150) + 2.000_r8*rxt(k,219)*y(k,232) &
                      + 2.000_r8*rxt(k,187)*y(k,246) + rxt(k,228)*y(k,250)
         mat(k,713) = mat(k,713) + rxt(k,728)*y(k,135)
         mat(k,1617) = mat(k,1617) + rxt(k,734)*y(k,135)
         mat(k,1127) = rxt(k,460)*y(k,232)
         mat(k,905) = rxt(k,497)*y(k,232)
         mat(k,1712) = rxt(k,431)*y(k,232)
         mat(k,2181) = rxt(k,359)*y(k,17) + rxt(k,365)*y(k,19) + rxt(k,307)*y(k,56) &
                      + rxt(k,314)*y(k,59) + rxt(k,214)*y(k,76) + rxt(k,248)*y(k,126) &
                      + rxt(k,218)*y(k,133) + 2.000_r8*rxt(k,219)*y(k,135) &
                      + rxt(k,460)*y(k,222) + rxt(k,497)*y(k,223) + rxt(k,431) &
                      *y(k,226) + 2.000_r8*rxt(k,233)*y(k,232) + rxt(k,226)*y(k,250) &
                      + rxt(k,507)*y(k,253)
         mat(k,2282) = mat(k,2282) + 2.000_r8*rxt(k,187)*y(k,135)
         mat(k,79) = mat(k,79) + rxt(k,190)*y(k,133) + 2.000_r8*rxt(k,191)*y(k,134)
         mat(k,752) = rxt(k,834)*y(k,124)
         mat(k,1949) = rxt(k,323)*y(k,59) + rxt(k,234)*y(k,90) + rxt(k,227)*y(k,133) &
                      + rxt(k,228)*y(k,135) + rxt(k,226)*y(k,232)
         mat(k,862) = rxt(k,507)*y(k,232)
      end do
      end subroutine nlnmat05
      subroutine nlnmat06( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,2257) = -(rxt(k,187)*y(k,246) + rxt(k,197)*y(k,133) + rxt(k,219) &
                      *y(k,232) + rxt(k,222)*y(k,76) + rxt(k,228)*y(k,250) + rxt(k,243) &
                      *y(k,125) + rxt(k,256)*y(k,124) + rxt(k,309)*y(k,56) + rxt(k,360) &
                      *y(k,17) + rxt(k,456)*y(k,25) + rxt(k,493)*y(k,29) + rxt(k,536) &
                      *y(k,105) + rxt(k,554)*y(k,111) + rxt(k,600)*y(k,98) + rxt(k,664) &
                      *y(k,141) + rxt(k,688)*y(k,6) + rxt(k,693)*y(k,110) + rxt(k,728) &
                      *y(k,148) + rxt(k,734)*y(k,150))
         mat(k,2289) = -rxt(k,187)*y(k,135)
         mat(k,2499) = -rxt(k,197)*y(k,135)
         mat(k,2188) = -rxt(k,219)*y(k,135)
         mat(k,2546) = -rxt(k,222)*y(k,135)
         mat(k,1956) = -rxt(k,228)*y(k,135)
         mat(k,2039) = -rxt(k,243)*y(k,135)
         mat(k,2448) = -rxt(k,256)*y(k,135)
         mat(k,2081) = -rxt(k,309)*y(k,135)
         mat(k,1662) = -rxt(k,360)*y(k,135)
         mat(k,503) = -rxt(k,456)*y(k,135)
         mat(k,1167) = -rxt(k,493)*y(k,135)
         mat(k,1410) = -rxt(k,536)*y(k,135)
         mat(k,1541) = -rxt(k,554)*y(k,135)
         mat(k,979) = -rxt(k,600)*y(k,135)
         mat(k,606) = -rxt(k,664)*y(k,135)
         mat(k,1022) = -rxt(k,688)*y(k,135)
         mat(k,1054) = -rxt(k,693)*y(k,135)
         mat(k,715) = -rxt(k,728)*y(k,135)
         mat(k,1623) = -rxt(k,734)*y(k,135)
         mat(k,2499) = mat(k,2499) + rxt(k,200)*y(k,134)
         mat(k,1757) = rxt(k,200)*y(k,133)
         mat(k,1601) = .150_r8*rxt(k,473)*y(k,232)
         mat(k,2188) = mat(k,2188) + .150_r8*rxt(k,473)*y(k,225) + .150_r8*rxt(k,542) &
                      *y(k,238)
         mat(k,1563) = .150_r8*rxt(k,542)*y(k,232)
         mat(k,276) = -(rxt(k,735)*y(k,150))
         mat(k,1607) = -rxt(k,735)*y(k,136)
         mat(k,2557) = rxt(k,362)*y(k,59)
         mat(k,1968) = rxt(k,362)*y(k,19) + 2.000_r8*rxt(k,313)*y(k,59)
         mat(k,321) = -(rxt(k,725)*y(k,133) + rxt(k,726)*y(k,250))
         mat(k,2462) = -rxt(k,725)*y(k,137)
         mat(k,1847) = -rxt(k,726)*y(k,137)
         mat(k,1207) = rxt(k,526)*y(k,250)
         mat(k,2366) = .100_r8*rxt(k,706)*y(k,255)
         mat(k,1822) = rxt(k,526)*y(k,93)
         mat(k,1271) = .100_r8*rxt(k,706)*y(k,124)
         mat(k,439) = -(rxt(k,485)*y(k,250))
         mat(k,1863) = -rxt(k,485)*y(k,139)
         mat(k,2003) = rxt(k,487)*y(k,225)
         mat(k,1573) = rxt(k,487)*y(k,125)
         mat(k,1998) = rxt(k,668)*y(k,218)
         mat(k,665) = rxt(k,668)*y(k,125)
         mat(k,602) = -(rxt(k,663)*y(k,125) + rxt(k,664)*y(k,135))
         mat(k,2008) = -rxt(k,663)*y(k,141)
         mat(k,2210) = -rxt(k,664)*y(k,141)
         mat(k,122) = .070_r8*rxt(k,643)*y(k,250)
         mat(k,2387) = rxt(k,640)*y(k,224)
         mat(k,103) = .060_r8*rxt(k,662)*y(k,250)
         mat(k,168) = .070_r8*rxt(k,685)*y(k,250)
         mat(k,739) = rxt(k,640)*y(k,124)
         mat(k,1881) = .070_r8*rxt(k,643)*y(k,66) + .060_r8*rxt(k,662)*y(k,142) &
                      + .070_r8*rxt(k,685)*y(k,213)
         mat(k,101) = -(rxt(k,662)*y(k,250))
         mat(k,1807) = -rxt(k,662)*y(k,142)
         mat(k,93) = .530_r8*rxt(k,629)*y(k,250)
         mat(k,1807) = mat(k,1807) + .530_r8*rxt(k,629)*y(k,7)
         mat(k,281) = -(rxt(k,665)*y(k,250))
         mat(k,1838) = -rxt(k,665)*y(k,143)
         mat(k,2097) = rxt(k,659)*y(k,252)
         mat(k,539) = rxt(k,659)*y(k,232)
         mat(k,481) = -(rxt(k,510)*y(k,250))
         mat(k,1868) = -rxt(k,510)*y(k,146)
         mat(k,2117) = rxt(k,507)*y(k,253)
         mat(k,856) = rxt(k,507)*y(k,232)
         mat(k,356) = -(rxt(k,515)*y(k,250))
         mat(k,1851) = -rxt(k,515)*y(k,147)
         mat(k,2104) = .850_r8*rxt(k,512)*y(k,254)
         mat(k,1318) = .850_r8*rxt(k,512)*y(k,232)
         mat(k,711) = -(rxt(k,727)*y(k,134) + rxt(k,728)*y(k,135) + rxt(k,731) &
                      *y(k,250))
         mat(k,1737) = -rxt(k,727)*y(k,148)
         mat(k,2211) = -rxt(k,728)*y(k,148)
         mat(k,1890) = -rxt(k,731)*y(k,148)
         mat(k,1615) = -(rxt(k,729)*y(k,19) + rxt(k,730)*y(k,59) + rxt(k,732)*y(k,125) &
                      + rxt(k,733)*y(k,134) + rxt(k,734)*y(k,135) + rxt(k,735) &
                      *y(k,136) + rxt(k,736)*y(k,250))
         mat(k,2566) = -rxt(k,729)*y(k,150)
         mat(k,1978) = -rxt(k,730)*y(k,150)
         mat(k,2028) = -rxt(k,732)*y(k,150)
         mat(k,1748) = -rxt(k,733)*y(k,150)
         mat(k,2247) = -rxt(k,734)*y(k,150)
         mat(k,278) = -rxt(k,735)*y(k,150)
         mat(k,1945) = -rxt(k,736)*y(k,150)
         mat(k,2488) = rxt(k,725)*y(k,137)
         mat(k,1748) = mat(k,1748) + rxt(k,727)*y(k,148)
         mat(k,2247) = mat(k,2247) + rxt(k,728)*y(k,148)
         mat(k,325) = rxt(k,725)*y(k,133)
         mat(k,712) = rxt(k,727)*y(k,134) + rxt(k,728)*y(k,135) + rxt(k,731)*y(k,250)
         mat(k,1945) = mat(k,1945) + rxt(k,731)*y(k,148)
         mat(k,918) = -(rxt(k,740)*y(k,250))
         mat(k,1907) = -rxt(k,740)*y(k,151)
         mat(k,2561) = rxt(k,729)*y(k,150)
         mat(k,1972) = rxt(k,730)*y(k,150)
         mat(k,337) = rxt(k,742)*y(k,126) + (rxt(k,743)+.500_r8*rxt(k,746))*y(k,250)
         mat(k,2015) = rxt(k,732)*y(k,150)
         mat(k,2312) = rxt(k,742)*y(k,67)
         mat(k,1741) = rxt(k,733)*y(k,150)
         mat(k,2216) = rxt(k,734)*y(k,150)
         mat(k,277) = rxt(k,735)*y(k,150)
         mat(k,323) = rxt(k,726)*y(k,250)
         mat(k,1611) = rxt(k,729)*y(k,19) + rxt(k,730)*y(k,59) + rxt(k,732)*y(k,125) &
                      + rxt(k,733)*y(k,134) + rxt(k,734)*y(k,135) + rxt(k,735) &
                      *y(k,136) + rxt(k,736)*y(k,250)
         mat(k,1907) = mat(k,1907) + (rxt(k,743)+.500_r8*rxt(k,746))*y(k,67) &
                      + rxt(k,726)*y(k,137) + rxt(k,736)*y(k,150)
         mat(k,213) = -(rxt(k,741)*y(k,263))
         mat(k,2587) = -rxt(k,741)*y(k,152)
         mat(k,917) = rxt(k,740)*y(k,250)
         mat(k,1829) = rxt(k,740)*y(k,151)
         mat(k,993) = .2202005_r8*rxt(k,775)*y(k,135) + .2202005_r8*rxt(k,776) &
                      *y(k,250)
         mat(k,86) = .0023005_r8*rxt(k,777)*y(k,250)
         mat(k,954) = .0031005_r8*rxt(k,780)*y(k,250)
         mat(k,34) = .2381005_r8*rxt(k,781)*y(k,250)
         mat(k,1025) = .0508005_r8*rxt(k,783)*y(k,135) + .0508005_r8*rxt(k,784) &
                      *y(k,250)
         mat(k,2197) = .2202005_r8*rxt(k,775)*y(k,6) + .0508005_r8*rxt(k,783)*y(k,110)
         mat(k,40) = .5931005_r8*rxt(k,785)*y(k,250)
         mat(k,108) = .1364005_r8*rxt(k,786)*y(k,250)
         mat(k,153) = .1677005_r8*rxt(k,787)*y(k,250)
         mat(k,1793) = .2202005_r8*rxt(k,776)*y(k,6) + .0023005_r8*rxt(k,777)*y(k,7) &
                      + .0031005_r8*rxt(k,780)*y(k,98) + .2381005_r8*rxt(k,781) &
                      *y(k,104) + .0508005_r8*rxt(k,784)*y(k,110) &
                      + .5931005_r8*rxt(k,785)*y(k,172) + .1364005_r8*rxt(k,786) &
                      *y(k,180) + .1677005_r8*rxt(k,787)*y(k,211)
         mat(k,994) = .2067005_r8*rxt(k,775)*y(k,135) + .2067005_r8*rxt(k,776) &
                      *y(k,250)
         mat(k,87) = .0008005_r8*rxt(k,777)*y(k,250)
         mat(k,955) = .0035005_r8*rxt(k,780)*y(k,250)
         mat(k,35) = .1308005_r8*rxt(k,781)*y(k,250)
         mat(k,1026) = .1149005_r8*rxt(k,783)*y(k,135) + .1149005_r8*rxt(k,784) &
                      *y(k,250)
         mat(k,2198) = .2067005_r8*rxt(k,775)*y(k,6) + .1149005_r8*rxt(k,783)*y(k,110)
         mat(k,41) = .1534005_r8*rxt(k,785)*y(k,250)
         mat(k,109) = .0101005_r8*rxt(k,786)*y(k,250)
         mat(k,154) = .0174005_r8*rxt(k,787)*y(k,250)
         mat(k,1794) = .2067005_r8*rxt(k,776)*y(k,6) + .0008005_r8*rxt(k,777)*y(k,7) &
                      + .0035005_r8*rxt(k,780)*y(k,98) + .1308005_r8*rxt(k,781) &
                      *y(k,104) + .1149005_r8*rxt(k,784)*y(k,110) &
                      + .1534005_r8*rxt(k,785)*y(k,172) + .0101005_r8*rxt(k,786) &
                      *y(k,180) + .0174005_r8*rxt(k,787)*y(k,211)
         mat(k,995) = .0653005_r8*rxt(k,775)*y(k,135) + .0653005_r8*rxt(k,776) &
                      *y(k,250)
         mat(k,88) = .0843005_r8*rxt(k,777)*y(k,250)
         mat(k,956) = .0003005_r8*rxt(k,780)*y(k,250)
         mat(k,36) = .0348005_r8*rxt(k,781)*y(k,250)
         mat(k,1027) = .0348005_r8*rxt(k,783)*y(k,135) + .0348005_r8*rxt(k,784) &
                      *y(k,250)
         mat(k,2199) = .0653005_r8*rxt(k,775)*y(k,6) + .0348005_r8*rxt(k,783)*y(k,110)
         mat(k,42) = .0459005_r8*rxt(k,785)*y(k,250)
         mat(k,110) = .0763005_r8*rxt(k,786)*y(k,250)
         mat(k,155) = .086_r8*rxt(k,787)*y(k,250)
         mat(k,1795) = .0653005_r8*rxt(k,776)*y(k,6) + .0843005_r8*rxt(k,777)*y(k,7) &
                      + .0003005_r8*rxt(k,780)*y(k,98) + .0348005_r8*rxt(k,781) &
                      *y(k,104) + .0348005_r8*rxt(k,784)*y(k,110) &
                      + .0459005_r8*rxt(k,785)*y(k,172) + .0763005_r8*rxt(k,786) &
                      *y(k,180) + .086_r8*rxt(k,787)*y(k,211)
         mat(k,996) = .1749305_r8*rxt(k,774)*y(k,126) + .1284005_r8*rxt(k,775) &
                      *y(k,135) + .1284005_r8*rxt(k,776)*y(k,250)
         mat(k,89) = .0443005_r8*rxt(k,777)*y(k,250)
         mat(k,957) = .0590245_r8*rxt(k,778)*y(k,126) + .0033005_r8*rxt(k,779) &
                      *y(k,135) + .0271005_r8*rxt(k,780)*y(k,250)
         mat(k,37) = .0076005_r8*rxt(k,781)*y(k,250)
         mat(k,1028) = .1749305_r8*rxt(k,782)*y(k,126) + .0554005_r8*rxt(k,783) &
                      *y(k,135) + .0554005_r8*rxt(k,784)*y(k,250)
         mat(k,2298) = .1749305_r8*rxt(k,774)*y(k,6) + .0590245_r8*rxt(k,778)*y(k,98) &
                      + .1749305_r8*rxt(k,782)*y(k,110)
         mat(k,2200) = .1284005_r8*rxt(k,775)*y(k,6) + .0033005_r8*rxt(k,779)*y(k,98) &
                      + .0554005_r8*rxt(k,783)*y(k,110)
         mat(k,43) = .0085005_r8*rxt(k,785)*y(k,250)
         mat(k,111) = .2157005_r8*rxt(k,786)*y(k,250)
         mat(k,156) = .0512005_r8*rxt(k,787)*y(k,250)
         mat(k,1796) = .1284005_r8*rxt(k,776)*y(k,6) + .0443005_r8*rxt(k,777)*y(k,7) &
                      + .0271005_r8*rxt(k,780)*y(k,98) + .0076005_r8*rxt(k,781) &
                      *y(k,104) + .0554005_r8*rxt(k,784)*y(k,110) &
                      + .0085005_r8*rxt(k,785)*y(k,172) + .2157005_r8*rxt(k,786) &
                      *y(k,180) + .0512005_r8*rxt(k,787)*y(k,211)
         mat(k,997) = .5901905_r8*rxt(k,774)*y(k,126) + .114_r8*rxt(k,775)*y(k,135) &
                      + .114_r8*rxt(k,776)*y(k,250)
         mat(k,90) = .1621005_r8*rxt(k,777)*y(k,250)
         mat(k,958) = .0250245_r8*rxt(k,778)*y(k,126) + .0474005_r8*rxt(k,780) &
                      *y(k,250)
         mat(k,38) = .0113005_r8*rxt(k,781)*y(k,250)
         mat(k,1029) = .5901905_r8*rxt(k,782)*y(k,126) + .1278005_r8*rxt(k,783) &
                      *y(k,135) + .1278005_r8*rxt(k,784)*y(k,250)
         mat(k,2299) = .5901905_r8*rxt(k,774)*y(k,6) + .0250245_r8*rxt(k,778)*y(k,98) &
                      + .5901905_r8*rxt(k,782)*y(k,110)
         mat(k,2201) = .114_r8*rxt(k,775)*y(k,6) + .1278005_r8*rxt(k,783)*y(k,110)
         mat(k,44) = .0128005_r8*rxt(k,785)*y(k,250)
         mat(k,112) = .0232005_r8*rxt(k,786)*y(k,250)
         mat(k,157) = .1598005_r8*rxt(k,787)*y(k,250)
         mat(k,1797) = .114_r8*rxt(k,776)*y(k,6) + .1621005_r8*rxt(k,777)*y(k,7) &
                      + .0474005_r8*rxt(k,780)*y(k,98) + .0113005_r8*rxt(k,781) &
                      *y(k,104) + .1278005_r8*rxt(k,784)*y(k,110) &
                      + .0128005_r8*rxt(k,785)*y(k,172) + .0232005_r8*rxt(k,786) &
                      *y(k,180) + .1598005_r8*rxt(k,787)*y(k,211)
      end do
      end subroutine nlnmat06
      subroutine nlnmat07( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,45) = -(rxt(k,785)*y(k,250))
         mat(k,1799) = -rxt(k,785)*y(k,172)
         mat(k,115) = .100_r8*rxt(k,674)*y(k,250)
         mat(k,158) = .230_r8*rxt(k,677)*y(k,250)
         mat(k,1810) = .100_r8*rxt(k,674)*y(k,180) + .230_r8*rxt(k,677)*y(k,211)
         mat(k,548) = -(rxt(k,708)*y(k,250))
         mat(k,1875) = -rxt(k,708)*y(k,174)
         mat(k,2123) = rxt(k,705)*y(k,255)
         mat(k,1273) = rxt(k,705)*y(k,232)
         mat(k,624) = -(rxt(k,709)*y(k,250))
         mat(k,1883) = -rxt(k,709)*y(k,175)
         mat(k,2389) = .200_r8*rxt(k,699)*y(k,245) + .200_r8*rxt(k,712)*y(k,256)
         mat(k,1679) = .500_r8*rxt(k,697)*y(k,245)
         mat(k,1296) = .200_r8*rxt(k,699)*y(k,124) + .500_r8*rxt(k,697)*y(k,226)
         mat(k,1246) = .200_r8*rxt(k,712)*y(k,124)
         mat(k,423) = -(rxt(k,714)*y(k,250))
         mat(k,1860) = -rxt(k,714)*y(k,176)
         mat(k,2113) = rxt(k,711)*y(k,256)
         mat(k,1245) = rxt(k,711)*y(k,232)
         mat(k,1137) = -(rxt(k,715)*y(k,126) + rxt(k,716)*y(k,250))
         mat(k,2323) = -rxt(k,715)*y(k,177)
         mat(k,1921) = -rxt(k,716)*y(k,177)
         mat(k,1011) = .330_r8*rxt(k,688)*y(k,135)
         mat(k,1043) = .330_r8*rxt(k,693)*y(k,135)
         mat(k,2417) = .800_r8*rxt(k,699)*y(k,245) + .800_r8*rxt(k,712)*y(k,256)
         mat(k,2323) = mat(k,2323) + rxt(k,700)*y(k,245)
         mat(k,2227) = .330_r8*rxt(k,688)*y(k,6) + .330_r8*rxt(k,693)*y(k,110)
         mat(k,625) = rxt(k,709)*y(k,250)
         mat(k,1690) = .500_r8*rxt(k,697)*y(k,245) + rxt(k,710)*y(k,256)
         mat(k,1301) = .800_r8*rxt(k,699)*y(k,124) + rxt(k,700)*y(k,126) &
                      + .500_r8*rxt(k,697)*y(k,226)
         mat(k,1921) = mat(k,1921) + rxt(k,709)*y(k,175)
         mat(k,1253) = .800_r8*rxt(k,712)*y(k,124) + rxt(k,710)*y(k,226)
         mat(k,1225) = -(rxt(k,718)*y(k,250))
         mat(k,1927) = -rxt(k,718)*y(k,178)
         mat(k,1012) = .300_r8*rxt(k,688)*y(k,135)
         mat(k,1044) = .300_r8*rxt(k,693)*y(k,135)
         mat(k,2422) = .900_r8*rxt(k,706)*y(k,255)
         mat(k,2231) = .300_r8*rxt(k,688)*y(k,6) + .300_r8*rxt(k,693)*y(k,110)
         mat(k,1694) = rxt(k,704)*y(k,255)
         mat(k,1281) = .900_r8*rxt(k,706)*y(k,124) + rxt(k,704)*y(k,226)
         mat(k,559) = -(rxt(k,673)*y(k,250))
         mat(k,1876) = -rxt(k,673)*y(k,179)
         mat(k,2124) = rxt(k,670)*y(k,257)
         mat(k,785) = rxt(k,670)*y(k,232)
         mat(k,113) = -(rxt(k,674)*y(k,250))
         mat(k,1808) = -rxt(k,674)*y(k,180)
         mat(k,61) = -(rxt(k,607)*y(k,250))
         mat(k,1803) = -rxt(k,607)*y(k,181)
         mat(k,1083) = rxt(k,568)*y(k,219)
         mat(k,1057) = rxt(k,568)*y(k,195)
         mat(k,225) = -(rxt(k,386)*y(k,133))
         mat(k,2461) = -rxt(k,386)*y(k,182)
         mat(k,2556) = rxt(k,373)*y(k,196)
         mat(k,1110) = rxt(k,373)*y(k,19)
         mat(k,430) = -(rxt(k,346)*y(k,56) + rxt(k,347)*y(k,133) + rxt(k,348)*y(k,250) &
                      + (rxt(k,811) + rxt(k,817) + rxt(k,822)) * y(k,85))
         mat(k,2051) = -rxt(k,346)*y(k,183)
         mat(k,2463) = -rxt(k,347)*y(k,183)
         mat(k,1861) = -rxt(k,348)*y(k,183)
         mat(k,1632) = -(rxt(k,811) + rxt(k,817) + rxt(k,822)) * y(k,183)
         mat(k,1969) = rxt(k,325)*y(k,196)
         mat(k,1112) = rxt(k,325)*y(k,59)
         mat(k,1080) = -(rxt(k,277)*y(k,250) + rxt(k,396)*y(k,73))
         mat(k,1916) = -rxt(k,277)*y(k,184)
         mat(k,871) = -rxt(k,396)*y(k,184)
         mat(k,1768) = rxt(k,427)*y(k,199)
         mat(k,1195) = rxt(k,469)*y(k,199)
         mat(k,1413) = rxt(k,504)*y(k,199)
         mat(k,338) = rxt(k,744)*y(k,199)
         mat(k,1637) = (rxt(k,811)+rxt(k,817)+rxt(k,822))*y(k,183)
         mat(k,431) = (rxt(k,811)+rxt(k,817)+rxt(k,822))*y(k,85)
         mat(k,1114) = rxt(k,272)*y(k,250)
         mat(k,1101) = rxt(k,427)*y(k,42) + rxt(k,469)*y(k,45) + rxt(k,504)*y(k,49) &
                      + rxt(k,744)*y(k,67)
         mat(k,1916) = mat(k,1916) + rxt(k,272)*y(k,196)
         mat(k,129) = -(rxt(k,282)*y(k,250))
         mat(k,1811) = -rxt(k,282)*y(k,185)
         mat(k,1106) = rxt(k,270)*y(k,232)
         mat(k,2091) = rxt(k,270)*y(k,196)
         mat(k,298) = -(rxt(k,561)*y(k,250))
         mat(k,1841) = -rxt(k,561)*y(k,186)
         mat(k,230) = .300_r8*rxt(k,608)*y(k,250)
         mat(k,235) = .500_r8*rxt(k,609)*y(k,250)
         mat(k,1088) = rxt(k,525)*y(k,229) + rxt(k,535)*y(k,236)
         mat(k,697) = rxt(k,525)*y(k,195)
         mat(k,1497) = rxt(k,535)*y(k,195)
         mat(k,1841) = mat(k,1841) + .300_r8*rxt(k,608)*y(k,187) + .500_r8*rxt(k,609) &
                      *y(k,188)
         mat(k,229) = -(rxt(k,608)*y(k,250))
         mat(k,1832) = -rxt(k,608)*y(k,187)
         mat(k,1086) = .080_r8*rxt(k,579)*y(k,234)
         mat(k,1464) = .080_r8*rxt(k,579)*y(k,195)
         mat(k,234) = -(rxt(k,609)*y(k,250))
         mat(k,1833) = -rxt(k,609)*y(k,188)
         mat(k,1087) = .080_r8*rxt(k,587)*y(k,235)
         mat(k,1426) = .080_r8*rxt(k,587)*y(k,195)
         mat(k,333) = -(rxt(k,610)*y(k,225) + rxt(k,611)*y(k,226) + rxt(k,612) &
                      *y(k,232) + rxt(k,613)*y(k,124) + rxt(k,614)*y(k,126))
         mat(k,1572) = -rxt(k,610)*y(k,189)
         mat(k,1674) = -rxt(k,611)*y(k,189)
         mat(k,2102) = -rxt(k,612)*y(k,189)
         mat(k,2373) = -rxt(k,613)*y(k,189)
         mat(k,2304) = -rxt(k,614)*y(k,189)
         mat(k,959) = rxt(k,603)*y(k,199)
         mat(k,1098) = rxt(k,603)*y(k,98)
         mat(k,132) = -(rxt(k,615)*y(k,250))
         mat(k,1812) = -rxt(k,615)*y(k,190)
         mat(k,329) = rxt(k,612)*y(k,232)
         mat(k,2092) = rxt(k,612)*y(k,189)
         mat(k,135) = -(rxt(k,562)*y(k,250))
         mat(k,1813) = -rxt(k,562)*y(k,191)
         mat(k,1107) = rxt(k,558)*y(k,238)
         mat(k,1547) = rxt(k,558)*y(k,196)
         mat(k,138) = -(rxt(k,285)*y(k,124) + (rxt(k,286) + rxt(k,287) + rxt(k,288) &
                      ) * y(k,125) + rxt(k,289)*y(k,134) + rxt(k,297)*y(k,250))
         mat(k,2365) = -rxt(k,285)*y(k,192)
         mat(k,2000) = -(rxt(k,286) + rxt(k,287) + rxt(k,288)) * y(k,192)
         mat(k,1728) = -rxt(k,289)*y(k,192)
         mat(k,1814) = -rxt(k,297)*y(k,192)
         mat(k,2460) = rxt(k,283)*y(k,258)
         mat(k,84) = rxt(k,283)*y(k,133)
         mat(k,140) = -(rxt(k,616)*y(k,250))
         mat(k,1815) = -rxt(k,616)*y(k,193)
         mat(k,330) = .200_r8*rxt(k,611)*y(k,226)
         mat(k,1668) = .200_r8*rxt(k,611)*y(k,189)
         mat(k,307) = -(rxt(k,617)*y(k,250))
         mat(k,1843) = -rxt(k,617)*y(k,194)
         mat(k,2370) = rxt(k,613)*y(k,189)
         mat(k,2302) = rxt(k,614)*y(k,189)
         mat(k,332) = rxt(k,613)*y(k,124) + rxt(k,614)*y(k,126) + rxt(k,610)*y(k,225) &
                      + .800_r8*rxt(k,611)*y(k,226)
         mat(k,1570) = rxt(k,610)*y(k,189)
         mat(k,1671) = .800_r8*rxt(k,611)*y(k,189)
         mat(k,1092) = -(rxt(k,252)*y(k,126) + rxt(k,260)*y(k,112) + rxt(k,298) &
                      *y(k,232) + rxt(k,299)*y(k,135) + rxt(k,300)*y(k,133) + rxt(k,324) &
                      *y(k,59) + rxt(k,367)*y(k,19) + rxt(k,433)*y(k,226) + rxt(k,443) &
                      *y(k,233) + rxt(k,462)*y(k,222) + rxt(k,475)*y(k,225) + rxt(k,480) &
                      *y(k,231) + rxt(k,499)*y(k,223) + rxt(k,509)*y(k,253) + rxt(k,514) &
                      *y(k,254) + (rxt(k,524) + rxt(k,525)) * y(k,229) + (rxt(k,534) &
                      + rxt(k,535)) * y(k,236) + rxt(k,546)*y(k,238) + rxt(k,550) &
                      *y(k,240) + (rxt(k,567) + rxt(k,568)) * y(k,219) + rxt(k,579) &
                      *y(k,234) + rxt(k,587)*y(k,235) + rxt(k,597)*y(k,101) + rxt(k,623) &
                      *y(k,260) + rxt(k,628)*y(k,218) + rxt(k,632)*y(k,220) + rxt(k,638) &
                      *y(k,221) + rxt(k,641)*y(k,224) + rxt(k,647)*y(k,227) + rxt(k,652) &
                      *y(k,237) + rxt(k,657)*y(k,239) + rxt(k,661)*y(k,252) + rxt(k,672) &
                      *y(k,257) + rxt(k,680)*y(k,261) + rxt(k,684)*y(k,262) + rxt(k,701) &
                      *y(k,245) + rxt(k,707)*y(k,255) + rxt(k,713)*y(k,256))
         mat(k,2319) = -rxt(k,252)*y(k,195)
         mat(k,931) = -rxt(k,260)*y(k,195)
         mat(k,2153) = -rxt(k,298)*y(k,195)
         mat(k,2223) = -rxt(k,299)*y(k,195)
         mat(k,2483) = -rxt(k,300)*y(k,195)
         mat(k,1975) = -rxt(k,324)*y(k,195)
         mat(k,2563) = -rxt(k,367)*y(k,195)
         mat(k,1686) = -rxt(k,433)*y(k,195)
         mat(k,525) = -rxt(k,443)*y(k,195)
         mat(k,1121) = -rxt(k,462)*y(k,195)
         mat(k,1578) = -rxt(k,475)*y(k,195)
         mat(k,883) = -rxt(k,480)*y(k,195)
         mat(k,901) = -rxt(k,499)*y(k,195)
         mat(k,858) = -rxt(k,509)*y(k,195)
         mat(k,1320) = -rxt(k,514)*y(k,195)
         mat(k,701) = -(rxt(k,524) + rxt(k,525)) * y(k,195)
         mat(k,1501) = -(rxt(k,534) + rxt(k,535)) * y(k,195)
         mat(k,1551) = -rxt(k,546)*y(k,195)
         mat(k,720) = -rxt(k,550)*y(k,195)
         mat(k,1064) = -(rxt(k,567) + rxt(k,568)) * y(k,195)
         mat(k,1470) = -rxt(k,579)*y(k,195)
         mat(k,1437) = -rxt(k,587)*y(k,195)
         mat(k,1377) = -rxt(k,597)*y(k,195)
         mat(k,1355) = -rxt(k,623)*y(k,195)
         mat(k,668) = -rxt(k,628)*y(k,195)
         mat(k,594) = -rxt(k,632)*y(k,195)
         mat(k,517) = -rxt(k,638)*y(k,195)
         mat(k,741) = -rxt(k,641)*y(k,195)
         mat(k,828) = -rxt(k,647)*y(k,195)
         mat(k,802) = -rxt(k,652)*y(k,195)
         mat(k,984) = -rxt(k,657)*y(k,195)
         mat(k,542) = -rxt(k,661)*y(k,195)
         mat(k,792) = -rxt(k,672)*y(k,195)
         mat(k,818) = -rxt(k,680)*y(k,195)
         mat(k,617) = -rxt(k,684)*y(k,195)
         mat(k,1298) = -rxt(k,701)*y(k,195)
         mat(k,1277) = -rxt(k,707)*y(k,195)
         mat(k,1250) = -rxt(k,713)*y(k,195)
         mat(k,931) = mat(k,931) + rxt(k,262)*y(k,196)
         mat(k,2020) = rxt(k,287)*y(k,192)
         mat(k,2483) = mat(k,2483) + rxt(k,290)*y(k,196)
         mat(k,1743) = rxt(k,289)*y(k,192) + rxt(k,284)*y(k,258)
         mat(k,1612) = rxt(k,737)*y(k,196)
         mat(k,139) = rxt(k,287)*y(k,125) + rxt(k,289)*y(k,134) + rxt(k,297)*y(k,250)
         mat(k,1115) = rxt(k,262)*y(k,112) + rxt(k,290)*y(k,133) + rxt(k,737)*y(k,150)
         mat(k,1917) = rxt(k,297)*y(k,192)
         mat(k,85) = rxt(k,284)*y(k,134)
         mat(k,1117) = -((rxt(k,261) + rxt(k,262) + rxt(k,263)) * y(k,112) + rxt(k,270) &
                      *y(k,232) + rxt(k,271)*y(k,126) + rxt(k,272)*y(k,250) + rxt(k,273) &
                      *y(k,199) + (rxt(k,290) + rxt(k,292)) * y(k,133) + rxt(k,291) &
                      *y(k,135) + rxt(k,325)*y(k,59) + rxt(k,373)*y(k,19) + rxt(k,488) &
                      *y(k,225) + rxt(k,558)*y(k,238) + rxt(k,648)*y(k,227) + rxt(k,653) &
                      *y(k,237) + rxt(k,658)*y(k,239) + rxt(k,666)*y(k,141) + rxt(k,669) &
                      *y(k,218) + rxt(k,737)*y(k,150))
         mat(k,932) = -(rxt(k,261) + rxt(k,262) + rxt(k,263)) * y(k,196)
         mat(k,2155) = -rxt(k,270)*y(k,196)
         mat(k,2321) = -rxt(k,271)*y(k,196)
         mat(k,1919) = -rxt(k,272)*y(k,196)
         mat(k,1104) = -rxt(k,273)*y(k,196)
         mat(k,2485) = -(rxt(k,290) + rxt(k,292)) * y(k,196)
         mat(k,2225) = -rxt(k,291)*y(k,196)
         mat(k,1977) = -rxt(k,325)*y(k,196)
         mat(k,2565) = -rxt(k,373)*y(k,196)
         mat(k,1580) = -rxt(k,488)*y(k,196)
         mat(k,1553) = -rxt(k,558)*y(k,196)
         mat(k,829) = -rxt(k,648)*y(k,196)
         mat(k,803) = -rxt(k,653)*y(k,196)
         mat(k,985) = -rxt(k,658)*y(k,196)
         mat(k,604) = -rxt(k,666)*y(k,196)
         mat(k,669) = -rxt(k,669)*y(k,196)
         mat(k,1614) = -rxt(k,737)*y(k,196)
         mat(k,410) = rxt(k,521)*y(k,199)
         mat(k,2565) = mat(k,2565) + rxt(k,367)*y(k,195)
         mat(k,1977) = mat(k,1977) + rxt(k,324)*y(k,195)
         mat(k,1379) = rxt(k,597)*y(k,195) + rxt(k,598)*y(k,199)
         mat(k,2415) = rxt(k,294)*y(k,199) + .600_r8*rxt(k,721)*y(k,202)
         mat(k,2321) = mat(k,2321) + rxt(k,252)*y(k,195) + rxt(k,722)*y(k,202)
         mat(k,2485) = mat(k,2485) + rxt(k,300)*y(k,195) + rxt(k,295)*y(k,199)
         mat(k,2225) = mat(k,2225) + rxt(k,299)*y(k,195)
         mat(k,62) = rxt(k,607)*y(k,250)
         mat(k,131) = rxt(k,282)*y(k,250)
         mat(k,232) = .700_r8*rxt(k,608)*y(k,250)
         mat(k,1094) = rxt(k,367)*y(k,19) + rxt(k,324)*y(k,59) + rxt(k,597)*y(k,101) &
                      + rxt(k,252)*y(k,126) + rxt(k,300)*y(k,133) + rxt(k,299) &
                      *y(k,135) + rxt(k,628)*y(k,218) + rxt(k,567)*y(k,219) &
                      + rxt(k,632)*y(k,220) + rxt(k,638)*y(k,221) + rxt(k,462) &
                      *y(k,222) + rxt(k,499)*y(k,223) + rxt(k,641)*y(k,224) &
                      + rxt(k,475)*y(k,225) + rxt(k,433)*y(k,226) + rxt(k,647) &
                      *y(k,227) + rxt(k,524)*y(k,229) + rxt(k,480)*y(k,231) &
                      + rxt(k,298)*y(k,232) + rxt(k,443)*y(k,233) + .920_r8*rxt(k,579) &
                      *y(k,234) + .920_r8*rxt(k,587)*y(k,235) + rxt(k,534)*y(k,236) &
                      + rxt(k,652)*y(k,237) + rxt(k,546)*y(k,238) + rxt(k,657) &
                      *y(k,239) + rxt(k,550)*y(k,240) + rxt(k,701)*y(k,245) &
                      + rxt(k,661)*y(k,252) + rxt(k,509)*y(k,253) + rxt(k,514) &
                      *y(k,254) + .900_r8*rxt(k,707)*y(k,255) + .800_r8*rxt(k,713) &
                      *y(k,256) + rxt(k,672)*y(k,257) + rxt(k,623)*y(k,260) &
                      + rxt(k,680)*y(k,261) + rxt(k,684)*y(k,262)
         mat(k,1104) = mat(k,1104) + rxt(k,521)*y(k,16) + rxt(k,598)*y(k,101) &
                      + rxt(k,294)*y(k,124) + rxt(k,295)*y(k,133) + rxt(k,293) &
                      *y(k,232) + rxt(k,580)*y(k,234) + rxt(k,588)*y(k,235) &
                      + rxt(k,533)*y(k,236) + rxt(k,547)*y(k,238) + rxt(k,702) &
                      *y(k,245) + rxt(k,296)*y(k,250) + rxt(k,624)*y(k,260)
         mat(k,195) = rxt(k,518)*y(k,250)
         mat(k,438) = .600_r8*rxt(k,721)*y(k,124) + rxt(k,722)*y(k,126) &
                      + .500_r8*rxt(k,719)*y(k,226)
         mat(k,315) = rxt(k,724)*y(k,250)
         mat(k,669) = mat(k,669) + rxt(k,628)*y(k,195)
         mat(k,1065) = rxt(k,567)*y(k,195)
         mat(k,595) = rxt(k,632)*y(k,195)
         mat(k,518) = rxt(k,638)*y(k,195)
         mat(k,1123) = rxt(k,462)*y(k,195)
         mat(k,902) = rxt(k,499)*y(k,195)
         mat(k,742) = rxt(k,641)*y(k,195)
         mat(k,1580) = mat(k,1580) + rxt(k,475)*y(k,195)
         mat(k,1688) = rxt(k,433)*y(k,195) + .500_r8*rxt(k,719)*y(k,202)
         mat(k,829) = mat(k,829) + rxt(k,647)*y(k,195)
         mat(k,702) = rxt(k,524)*y(k,195)
         mat(k,884) = rxt(k,480)*y(k,195)
         mat(k,2155) = mat(k,2155) + rxt(k,298)*y(k,195) + rxt(k,293)*y(k,199)
         mat(k,526) = rxt(k,443)*y(k,195)
         mat(k,1472) = .920_r8*rxt(k,579)*y(k,195) + rxt(k,580)*y(k,199)
         mat(k,1439) = .920_r8*rxt(k,587)*y(k,195) + rxt(k,588)*y(k,199)
         mat(k,1503) = rxt(k,534)*y(k,195) + rxt(k,533)*y(k,199)
         mat(k,803) = mat(k,803) + rxt(k,652)*y(k,195)
         mat(k,1553) = mat(k,1553) + rxt(k,546)*y(k,195) + rxt(k,547)*y(k,199)
         mat(k,985) = mat(k,985) + rxt(k,657)*y(k,195)
         mat(k,721) = rxt(k,550)*y(k,195)
         mat(k,1300) = rxt(k,701)*y(k,195) + rxt(k,702)*y(k,199)
         mat(k,1919) = mat(k,1919) + rxt(k,607)*y(k,181) + rxt(k,282)*y(k,185) &
                      + .700_r8*rxt(k,608)*y(k,187) + rxt(k,296)*y(k,199) + rxt(k,518) &
                      *y(k,201) + rxt(k,724)*y(k,210)
         mat(k,543) = rxt(k,661)*y(k,195)
         mat(k,859) = rxt(k,509)*y(k,195)
         mat(k,1322) = rxt(k,514)*y(k,195)
         mat(k,1279) = .900_r8*rxt(k,707)*y(k,195)
         mat(k,1252) = .800_r8*rxt(k,713)*y(k,195)
         mat(k,793) = rxt(k,672)*y(k,195)
         mat(k,1357) = rxt(k,623)*y(k,195) + rxt(k,624)*y(k,199)
         mat(k,819) = rxt(k,680)*y(k,195)
         mat(k,618) = rxt(k,684)*y(k,195)
      end do
      end subroutine nlnmat07
      subroutine nlnmat08( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,2300) = rxt(k,271)*y(k,196)
         mat(k,1109) = rxt(k,271)*y(k,126)
         mat(k,1111) = rxt(k,273)*y(k,199)
         mat(k,1097) = rxt(k,273)*y(k,196)
         mat(k,1103) = -(rxt(k,269)*y(k,125) + rxt(k,273)*y(k,196) + rxt(k,293) &
                      *y(k,232) + rxt(k,294)*y(k,124) + rxt(k,295)*y(k,133) + rxt(k,296) &
                      *y(k,250) + rxt(k,427)*y(k,42) + rxt(k,469)*y(k,45) + rxt(k,494) &
                      *y(k,29) + rxt(k,504)*y(k,49) + rxt(k,521)*y(k,16) + rxt(k,533) &
                      *y(k,236) + rxt(k,547)*y(k,238) + rxt(k,580)*y(k,234) + rxt(k,588) &
                      *y(k,235) + rxt(k,598)*y(k,101) + rxt(k,603)*y(k,98) + rxt(k,624) &
                      *y(k,260) + rxt(k,690)*y(k,6) + rxt(k,695)*y(k,110) + rxt(k,702) &
                      *y(k,245) + rxt(k,717)*y(k,177) + rxt(k,744)*y(k,67))
         mat(k,2021) = -rxt(k,269)*y(k,199)
         mat(k,1116) = -rxt(k,273)*y(k,199)
         mat(k,2154) = -rxt(k,293)*y(k,199)
         mat(k,2414) = -rxt(k,294)*y(k,199)
         mat(k,2484) = -rxt(k,295)*y(k,199)
         mat(k,1918) = -rxt(k,296)*y(k,199)
         mat(k,1770) = -rxt(k,427)*y(k,199)
         mat(k,1196) = -rxt(k,469)*y(k,199)
         mat(k,1153) = -rxt(k,494)*y(k,199)
         mat(k,1414) = -rxt(k,504)*y(k,199)
         mat(k,409) = -rxt(k,521)*y(k,199)
         mat(k,1502) = -rxt(k,533)*y(k,199)
         mat(k,1552) = -rxt(k,547)*y(k,199)
         mat(k,1471) = -rxt(k,580)*y(k,199)
         mat(k,1438) = -rxt(k,588)*y(k,199)
         mat(k,1378) = -rxt(k,598)*y(k,199)
         mat(k,965) = -rxt(k,603)*y(k,199)
         mat(k,1356) = -rxt(k,624)*y(k,199)
         mat(k,1009) = -rxt(k,690)*y(k,199)
         mat(k,1041) = -rxt(k,695)*y(k,199)
         mat(k,1299) = -rxt(k,702)*y(k,199)
         mat(k,1135) = -rxt(k,717)*y(k,199)
         mat(k,339) = -rxt(k,744)*y(k,199)
         mat(k,2063) = rxt(k,346)*y(k,183)
         mat(k,872) = rxt(k,396)*y(k,184)
         mat(k,2484) = mat(k,2484) + rxt(k,386)*y(k,182) + rxt(k,347)*y(k,183) &
                      + rxt(k,292)*y(k,196)
         mat(k,2224) = rxt(k,291)*y(k,196)
         mat(k,227) = rxt(k,386)*y(k,133)
         mat(k,432) = rxt(k,346)*y(k,56) + rxt(k,347)*y(k,133) + rxt(k,348)*y(k,250)
         mat(k,1081) = rxt(k,396)*y(k,73) + rxt(k,277)*y(k,250)
         mat(k,136) = .500_r8*rxt(k,562)*y(k,250)
         mat(k,1116) = mat(k,1116) + rxt(k,292)*y(k,133) + rxt(k,291)*y(k,135)
         mat(k,151) = rxt(k,491)*y(k,250)
         mat(k,1918) = mat(k,1918) + rxt(k,348)*y(k,183) + rxt(k,277)*y(k,184) &
                      + .500_r8*rxt(k,562)*y(k,191) + rxt(k,491)*y(k,208)
         mat(k,1999) = rxt(k,269)*y(k,199)
         mat(k,1095) = rxt(k,269)*y(k,125)
         mat(k,194) = -(rxt(k,518)*y(k,250))
         mat(k,1825) = -rxt(k,518)*y(k,201)
         mat(k,1145) = rxt(k,494)*y(k,199)
         mat(k,233) = .500_r8*rxt(k,609)*y(k,250)
         mat(k,133) = rxt(k,615)*y(k,250)
         mat(k,141) = rxt(k,616)*y(k,250)
         mat(k,306) = rxt(k,617)*y(k,250)
         mat(k,1096) = rxt(k,494)*y(k,29)
         mat(k,1825) = mat(k,1825) + .500_r8*rxt(k,609)*y(k,188) + rxt(k,615)*y(k,190) &
                      + rxt(k,616)*y(k,193) + rxt(k,617)*y(k,194)
         mat(k,436) = -(rxt(k,719)*y(k,226) + rxt(k,720)*y(k,232) + rxt(k,721) &
                      *y(k,124) + rxt(k,722)*y(k,126))
         mat(k,1677) = -rxt(k,719)*y(k,202)
         mat(k,2114) = -rxt(k,720)*y(k,202)
         mat(k,2376) = -rxt(k,721)*y(k,202)
         mat(k,2307) = -rxt(k,722)*y(k,202)
         mat(k,1000) = rxt(k,690)*y(k,199)
         mat(k,1032) = rxt(k,695)*y(k,199)
         mat(k,1133) = .500_r8*rxt(k,717)*y(k,199)
         mat(k,1099) = rxt(k,690)*y(k,6) + rxt(k,695)*y(k,110) + .500_r8*rxt(k,717) &
                      *y(k,177)
         mat(k,239) = rxt(k,723)*y(k,250)
         mat(k,1862) = rxt(k,723)*y(k,203)
         mat(k,238) = -(rxt(k,723)*y(k,250))
         mat(k,1834) = -rxt(k,723)*y(k,203)
         mat(k,434) = rxt(k,720)*y(k,232)
         mat(k,2096) = rxt(k,720)*y(k,202)
         mat(k,493) = -(rxt(k,201)*y(k,133) + rxt(k,202)*y(k,134) + rxt(k,209) &
                      *y(k,135) + rxt(k,212)*y(k,79) + rxt(k,213)*y(k,77) + rxt(k,220) &
                      *y(k,232) + rxt(k,231)*y(k,250) + (rxt(k,245) + rxt(k,247) &
                      ) * y(k,125) + rxt(k,253)*y(k,126) + rxt(k,259)*y(k,124) &
                      + rxt(k,320)*y(k,60) + rxt(k,326)*y(k,59) + rxt(k,330)*y(k,85) &
                      + rxt(k,334)*y(k,92) + rxt(k,370)*y(k,20) + rxt(k,374)*y(k,19) &
                      + rxt(k,378)*y(k,81) + rxt(k,380)*y(k,91) + rxt(k,428)*y(k,42))
         mat(k,2467) = -rxt(k,201)*y(k,204)
         mat(k,1734) = -rxt(k,202)*y(k,204)
         mat(k,2208) = -rxt(k,209)*y(k,204)
         mat(k,572) = -rxt(k,212)*y(k,204)
         mat(k,1336) = -rxt(k,213)*y(k,204)
         mat(k,2119) = -rxt(k,220)*y(k,204)
         mat(k,1870) = -rxt(k,231)*y(k,204)
         mat(k,2006) = -(rxt(k,245) + rxt(k,247)) * y(k,204)
         mat(k,2308) = -rxt(k,253)*y(k,204)
         mat(k,2379) = -rxt(k,259)*y(k,204)
         mat(k,942) = -rxt(k,320)*y(k,204)
         mat(k,1970) = -rxt(k,326)*y(k,204)
         mat(k,1633) = -rxt(k,330)*y(k,204)
         mat(k,838) = -rxt(k,334)*y(k,204)
         mat(k,504) = -rxt(k,370)*y(k,204)
         mat(k,2558) = -rxt(k,374)*y(k,204)
         mat(k,847) = -rxt(k,378)*y(k,204)
         mat(k,729) = -rxt(k,380)*y(k,204)
         mat(k,1764) = -rxt(k,428)*y(k,204)
         mat(k,925) = rxt(k,261)*y(k,196)
         mat(k,1734) = mat(k,1734) + (rxt(k,206)+rxt(k,207))*y(k,259)
         mat(k,1113) = rxt(k,261)*y(k,112)
         mat(k,243) = (rxt(k,206)+rxt(k,207))*y(k,134)
         mat(k,490) = -(rxt(k,188)*y(k,246) + rxt(k,198)*y(k,133) + rxt(k,221) &
                      *y(k,232) + rxt(k,223)*y(k,76) + rxt(k,232)*y(k,250) + rxt(k,246) &
                      *y(k,125) + rxt(k,258)*y(k,124) + rxt(k,327)*y(k,56) + rxt(k,375) &
                      *y(k,17) + rxt(k,457)*y(k,25) + rxt(k,495)*y(k,29) + rxt(k,539) &
                      *y(k,105) + rxt(k,556)*y(k,111) + rxt(k,604)*y(k,98) + rxt(k,667) &
                      *y(k,141) + rxt(k,691)*y(k,6) + rxt(k,696)*y(k,110) + rxt(k,738) &
                      *y(k,150) + rxt(k,739)*y(k,148))
         mat(k,2271) = -rxt(k,188)*y(k,205)
         mat(k,2466) = -rxt(k,198)*y(k,205)
         mat(k,2118) = -rxt(k,221)*y(k,205)
         mat(k,2533) = -rxt(k,223)*y(k,205)
         mat(k,1869) = -rxt(k,232)*y(k,205)
         mat(k,2005) = -rxt(k,246)*y(k,205)
         mat(k,2378) = -rxt(k,258)*y(k,205)
         mat(k,2052) = -rxt(k,327)*y(k,205)
         mat(k,1653) = -rxt(k,375)*y(k,205)
         mat(k,494) = -rxt(k,457)*y(k,205)
         mat(k,1147) = -rxt(k,495)*y(k,205)
         mat(k,1398) = -rxt(k,539)*y(k,205)
         mat(k,1523) = -rxt(k,556)*y(k,205)
         mat(k,960) = -rxt(k,604)*y(k,205)
         mat(k,600) = -rxt(k,667)*y(k,205)
         mat(k,1001) = -rxt(k,691)*y(k,205)
         mat(k,1033) = -rxt(k,696)*y(k,205)
         mat(k,1608) = -rxt(k,738)*y(k,205)
         mat(k,709) = -rxt(k,739)*y(k,205)
         mat(k,1733) = rxt(k,202)*y(k,204)
         mat(k,492) = rxt(k,202)*y(k,134)
         mat(k,297) = rxt(k,561)*y(k,250)
         mat(k,1085) = .100_r8*rxt(k,707)*y(k,255)
         mat(k,1826) = rxt(k,561)*y(k,186)
         mat(k,1272) = .100_r8*rxt(k,707)*y(k,195)
         mat(k,147) = -(rxt(k,625)*y(k,250))
         mat(k,1817) = -rxt(k,625)*y(k,207)
         mat(k,2093) = rxt(k,620)*y(k,260)
         mat(k,1352) = rxt(k,620)*y(k,232)
         mat(k,150) = -(rxt(k,491)*y(k,250))
         mat(k,1818) = -rxt(k,491)*y(k,208)
         mat(k,1108) = rxt(k,488)*y(k,225)
         mat(k,1569) = rxt(k,488)*y(k,196)
         mat(k,1105) = rxt(k,669)*y(k,218)
         mat(k,664) = rxt(k,669)*y(k,196)
         mat(k,313) = -(rxt(k,724)*y(k,250))
         mat(k,1845) = -rxt(k,724)*y(k,210)
         mat(k,2372) = .200_r8*rxt(k,721)*y(k,202)
         mat(k,1089) = .200_r8*rxt(k,713)*y(k,256)
         mat(k,435) = .200_r8*rxt(k,721)*y(k,124) + .500_r8*rxt(k,719)*y(k,226)
         mat(k,1673) = .500_r8*rxt(k,719)*y(k,202)
         mat(k,1244) = .200_r8*rxt(k,713)*y(k,195)
         mat(k,159) = -(rxt(k,677)*y(k,250))
         mat(k,1819) = -rxt(k,677)*y(k,211)
         mat(k,686) = -(rxt(k,681)*y(k,250))
         mat(k,1888) = -rxt(k,681)*y(k,212)
         mat(k,2133) = rxt(k,678)*y(k,261)
         mat(k,812) = rxt(k,678)*y(k,232)
         mat(k,167) = -(rxt(k,685)*y(k,250))
         mat(k,1820) = -rxt(k,685)*y(k,213)
         mat(k,160) = .150_r8*rxt(k,677)*y(k,250)
         mat(k,1820) = mat(k,1820) + .150_r8*rxt(k,677)*y(k,211)
         mat(k,392) = -(rxt(k,686)*y(k,250))
         mat(k,1856) = -rxt(k,686)*y(k,214)
         mat(k,2109) = rxt(k,682)*y(k,262)
         mat(k,614) = rxt(k,682)*y(k,232)
         mat(k,666) = -(rxt(k,626)*y(k,232) + rxt(k,627)*y(k,124) + rxt(k,668) &
                      *y(k,125))
         mat(k,2131) = -rxt(k,626)*y(k,218)
         mat(k,2391) = -rxt(k,627)*y(k,218)
         mat(k,2009) = -rxt(k,668)*y(k,218)
         mat(k,187) = rxt(k,634)*y(k,250)
         mat(k,1886) = rxt(k,634)*y(k,22)
         mat(k,1062) = -(rxt(k,564)*y(k,232) + (rxt(k,565) + rxt(k,566)) * y(k,124))
         mat(k,2150) = -rxt(k,564)*y(k,219)
         mat(k,2410) = -(rxt(k,565) + rxt(k,566)) * y(k,219)
         mat(k,656) = rxt(k,569)*y(k,250)
         mat(k,181) = rxt(k,570)*y(k,250)
         mat(k,1914) = rxt(k,569)*y(k,2) + rxt(k,570)*y(k,15)
         mat(k,591) = -(rxt(k,630)*y(k,232) + rxt(k,631)*y(k,124))
         mat(k,2127) = -rxt(k,630)*y(k,220)
         mat(k,2386) = -rxt(k,631)*y(k,220)
         mat(k,94) = .350_r8*rxt(k,629)*y(k,250)
         mat(k,364) = rxt(k,633)*y(k,250)
         mat(k,1880) = .350_r8*rxt(k,629)*y(k,7) + rxt(k,633)*y(k,8)
         mat(k,515) = -(rxt(k,635)*y(k,232) + rxt(k,637)*y(k,124))
         mat(k,2120) = -rxt(k,635)*y(k,221)
         mat(k,2380) = -rxt(k,637)*y(k,221)
         mat(k,288) = rxt(k,636)*y(k,250)
         mat(k,116) = .070_r8*rxt(k,674)*y(k,250)
         mat(k,161) = .060_r8*rxt(k,677)*y(k,250)
         mat(k,1872) = rxt(k,636)*y(k,23) + .070_r8*rxt(k,674)*y(k,180) &
                      + .060_r8*rxt(k,677)*y(k,211)
         mat(k,1124) = -(4._r8*rxt(k,458)*y(k,222) + rxt(k,459)*y(k,226) + rxt(k,460) &
                      *y(k,232) + rxt(k,461)*y(k,124))
         mat(k,1689) = -rxt(k,459)*y(k,222)
         mat(k,2156) = -rxt(k,460)*y(k,222)
         mat(k,2416) = -rxt(k,461)*y(k,222)
         mat(k,293) = .500_r8*rxt(k,464)*y(k,250)
         mat(k,260) = rxt(k,465)*y(k,56) + rxt(k,466)*y(k,250)
         mat(k,2065) = rxt(k,465)*y(k,28)
         mat(k,1920) = .500_r8*rxt(k,464)*y(k,27) + rxt(k,466)*y(k,28)
      end do
      end subroutine nlnmat08
      subroutine nlnmat09( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,899) = -(rxt(k,496)*y(k,226) + rxt(k,497)*y(k,232) + rxt(k,498) &
                      *y(k,124))
         mat(k,1682) = -rxt(k,496)*y(k,223)
         mat(k,2147) = -rxt(k,497)*y(k,223)
         mat(k,2406) = -rxt(k,498)*y(k,223)
         mat(k,375) = rxt(k,500)*y(k,250)
         mat(k,58) = rxt(k,501)*y(k,250)
         mat(k,1905) = rxt(k,500)*y(k,30) + rxt(k,501)*y(k,31)
         mat(k,740) = -(rxt(k,639)*y(k,232) + rxt(k,640)*y(k,124))
         mat(k,2136) = -rxt(k,639)*y(k,224)
         mat(k,2394) = -rxt(k,640)*y(k,224)
         mat(k,223) = rxt(k,642)*y(k,250)
         mat(k,2394) = mat(k,2394) + rxt(k,627)*y(k,218)
         mat(k,2212) = rxt(k,664)*y(k,141)
         mat(k,603) = rxt(k,664)*y(k,135)
         mat(k,667) = rxt(k,627)*y(k,124) + .400_r8*rxt(k,626)*y(k,232)
         mat(k,2136) = mat(k,2136) + .400_r8*rxt(k,626)*y(k,218)
         mat(k,1892) = rxt(k,642)*y(k,32)
         mat(k,1594) = -(4._r8*rxt(k,471)*y(k,225) + rxt(k,472)*y(k,226) + rxt(k,473) &
                      *y(k,232) + rxt(k,474)*y(k,124) + rxt(k,487)*y(k,125) + rxt(k,527) &
                      *y(k,236) + rxt(k,574)*y(k,234) + rxt(k,581)*y(k,235) + rxt(k,592) &
                      *y(k,101) + rxt(k,618)*y(k,260))
         mat(k,1710) = -rxt(k,472)*y(k,225)
         mat(k,2177) = -rxt(k,473)*y(k,225)
         mat(k,2438) = -rxt(k,474)*y(k,225)
         mat(k,2027) = -rxt(k,487)*y(k,225)
         mat(k,1512) = -rxt(k,527)*y(k,225)
         mat(k,1485) = -rxt(k,574)*y(k,225)
         mat(k,1452) = -rxt(k,581)*y(k,225)
         mat(k,1388) = -rxt(k,592)*y(k,225)
         mat(k,1363) = -rxt(k,618)*y(k,225)
         mat(k,1018) = .060_r8*rxt(k,688)*y(k,135)
         mat(k,1200) = rxt(k,467)*y(k,126) + rxt(k,468)*y(k,250)
         mat(k,1418) = rxt(k,502)*y(k,126) + rxt(k,503)*y(k,250)
         mat(k,454) = .500_r8*rxt(k,477)*y(k,250)
         mat(k,974) = .080_r8*rxt(k,600)*y(k,135)
         mat(k,1406) = .100_r8*rxt(k,536)*y(k,135)
         mat(k,1050) = .060_r8*rxt(k,693)*y(k,135)
         mat(k,1534) = .280_r8*rxt(k,554)*y(k,135)
         mat(k,2438) = mat(k,2438) + .530_r8*rxt(k,531)*y(k,236) + rxt(k,544)*y(k,238) &
                      + rxt(k,549)*y(k,240) + rxt(k,513)*y(k,254)
         mat(k,2345) = rxt(k,467)*y(k,45) + rxt(k,502)*y(k,49) + .530_r8*rxt(k,530) &
                      *y(k,236) + rxt(k,545)*y(k,238)
         mat(k,2246) = .060_r8*rxt(k,688)*y(k,6) + .080_r8*rxt(k,600)*y(k,98) &
                      + .100_r8*rxt(k,536)*y(k,105) + .060_r8*rxt(k,693)*y(k,110) &
                      + .280_r8*rxt(k,554)*y(k,111)
         mat(k,1228) = .650_r8*rxt(k,718)*y(k,250)
         mat(k,1594) = mat(k,1594) + .530_r8*rxt(k,527)*y(k,236)
         mat(k,1710) = mat(k,1710) + .260_r8*rxt(k,528)*y(k,236) + rxt(k,541)*y(k,238) &
                      + .300_r8*rxt(k,511)*y(k,254)
         mat(k,2177) = mat(k,2177) + .450_r8*rxt(k,542)*y(k,238) + .200_r8*rxt(k,548) &
                      *y(k,240) + .150_r8*rxt(k,512)*y(k,254)
         mat(k,1512) = mat(k,1512) + .530_r8*rxt(k,531)*y(k,124) + .530_r8*rxt(k,530) &
                      *y(k,126) + .530_r8*rxt(k,527)*y(k,225) + .260_r8*rxt(k,528) &
                      *y(k,226)
         mat(k,1557) = rxt(k,544)*y(k,124) + rxt(k,545)*y(k,126) + rxt(k,541)*y(k,226) &
                      + .450_r8*rxt(k,542)*y(k,232) + 4.000_r8*rxt(k,543)*y(k,238)
         mat(k,723) = rxt(k,549)*y(k,124) + .200_r8*rxt(k,548)*y(k,232)
         mat(k,1944) = rxt(k,468)*y(k,45) + rxt(k,503)*y(k,49) + .500_r8*rxt(k,477) &
                      *y(k,51) + .650_r8*rxt(k,718)*y(k,178)
         mat(k,1326) = rxt(k,513)*y(k,124) + .300_r8*rxt(k,511)*y(k,226) &
                      + .150_r8*rxt(k,512)*y(k,232)
         mat(k,1711) = -(rxt(k,310)*y(k,59) + (4._r8*rxt(k,429) + 4._r8*rxt(k,430) &
                      ) * y(k,226) + rxt(k,431)*y(k,232) + rxt(k,432)*y(k,124) &
                      + rxt(k,459)*y(k,222) + rxt(k,472)*y(k,225) + rxt(k,496) &
                      *y(k,223) + rxt(k,511)*y(k,254) + rxt(k,528)*y(k,236) + rxt(k,541) &
                      *y(k,238) + rxt(k,575)*y(k,234) + rxt(k,582)*y(k,235) + rxt(k,593) &
                      *y(k,101) + rxt(k,619)*y(k,260) + rxt(k,697)*y(k,245) + rxt(k,704) &
                      *y(k,255) + rxt(k,710)*y(k,256))
         mat(k,1981) = -rxt(k,310)*y(k,226)
         mat(k,2180) = -rxt(k,431)*y(k,226)
         mat(k,2440) = -rxt(k,432)*y(k,226)
         mat(k,1126) = -rxt(k,459)*y(k,226)
         mat(k,1595) = -rxt(k,472)*y(k,226)
         mat(k,904) = -rxt(k,496)*y(k,226)
         mat(k,1327) = -rxt(k,511)*y(k,226)
         mat(k,1513) = -rxt(k,528)*y(k,226)
         mat(k,1558) = -rxt(k,541)*y(k,226)
         mat(k,1486) = -rxt(k,575)*y(k,226)
         mat(k,1453) = -rxt(k,582)*y(k,226)
         mat(k,1389) = -rxt(k,593)*y(k,226)
         mat(k,1364) = -rxt(k,619)*y(k,226)
         mat(k,1308) = -rxt(k,697)*y(k,226)
         mat(k,1286) = -rxt(k,704)*y(k,226)
         mat(k,1261) = -rxt(k,710)*y(k,226)
         mat(k,1161) = .280_r8*rxt(k,493)*y(k,135)
         mat(k,470) = rxt(k,476)*y(k,250)
         mat(k,387) = .700_r8*rxt(k,435)*y(k,250)
         mat(k,975) = .050_r8*rxt(k,600)*y(k,135)
         mat(k,1389) = mat(k,1389) + rxt(k,592)*y(k,225)
         mat(k,2440) = mat(k,2440) + rxt(k,474)*y(k,225) + .830_r8*rxt(k,645)*y(k,227) &
                      + .170_r8*rxt(k,655)*y(k,239)
         mat(k,2249) = .280_r8*rxt(k,493)*y(k,29) + .050_r8*rxt(k,600)*y(k,98)
         mat(k,1595) = mat(k,1595) + rxt(k,592)*y(k,101) + rxt(k,474)*y(k,124) &
                      + 4.000_r8*rxt(k,471)*y(k,225) + .900_r8*rxt(k,472)*y(k,226) &
                      + .450_r8*rxt(k,473)*y(k,232) + rxt(k,574)*y(k,234) + rxt(k,581) &
                      *y(k,235) + rxt(k,527)*y(k,236) + rxt(k,540)*y(k,238) &
                      + rxt(k,618)*y(k,260)
         mat(k,1711) = mat(k,1711) + .900_r8*rxt(k,472)*y(k,225)
         mat(k,832) = .830_r8*rxt(k,645)*y(k,124) + .330_r8*rxt(k,644)*y(k,232)
         mat(k,2180) = mat(k,2180) + .450_r8*rxt(k,473)*y(k,225) + .330_r8*rxt(k,644) &
                      *y(k,227) + .070_r8*rxt(k,654)*y(k,239)
         mat(k,1486) = mat(k,1486) + rxt(k,574)*y(k,225)
         mat(k,1453) = mat(k,1453) + rxt(k,581)*y(k,225)
         mat(k,1513) = mat(k,1513) + rxt(k,527)*y(k,225)
         mat(k,1558) = mat(k,1558) + rxt(k,540)*y(k,225)
         mat(k,988) = .170_r8*rxt(k,655)*y(k,124) + .070_r8*rxt(k,654)*y(k,232)
         mat(k,1948) = rxt(k,476)*y(k,50) + .700_r8*rxt(k,435)*y(k,53)
         mat(k,1364) = mat(k,1364) + rxt(k,618)*y(k,225)
         mat(k,827) = -(rxt(k,644)*y(k,232) + rxt(k,645)*y(k,124) + rxt(k,646) &
                      *y(k,125))
         mat(k,2141) = -rxt(k,644)*y(k,227)
         mat(k,2402) = -rxt(k,645)*y(k,227)
         mat(k,2013) = -rxt(k,646)*y(k,227)
         mat(k,768) = -(rxt(k,823)*y(k,243) + rxt(k,824)*y(k,249) + rxt(k,825) &
                      *y(k,242))
         mat(k,758) = -rxt(k,823)*y(k,228)
         mat(k,750) = -rxt(k,824)*y(k,228)
         mat(k,646) = -rxt(k,825)*y(k,228)
         mat(k,698) = -((rxt(k,522) + rxt(k,523)) * y(k,124))
         mat(k,2392) = -(rxt(k,522) + rxt(k,523)) * y(k,229)
         mat(k,407) = rxt(k,520)*y(k,250)
         mat(k,1889) = rxt(k,520)*y(k,16)
         mat(k,401) = -(rxt(k,482)*y(k,134))
         mat(k,1730) = -rxt(k,482)*y(k,230)
         mat(k,2375) = .750_r8*rxt(k,479)*y(k,231)
         mat(k,881) = .750_r8*rxt(k,479)*y(k,124)
         mat(k,882) = -(rxt(k,478)*y(k,232) + rxt(k,479)*y(k,124))
         mat(k,2145) = -rxt(k,478)*y(k,231)
         mat(k,2404) = -rxt(k,479)*y(k,231)
         mat(k,498) = rxt(k,486)*y(k,250)
         mat(k,1903) = rxt(k,486)*y(k,25)
         mat(k,2187) = -((rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,76) + rxt(k,218) &
                      *y(k,133) + rxt(k,219)*y(k,135) + rxt(k,226)*y(k,250) &
                      + 4._r8*rxt(k,233)*y(k,232) + rxt(k,248)*y(k,126) + rxt(k,255) &
                      *y(k,124) + rxt(k,266)*y(k,125) + (rxt(k,307) + rxt(k,308) &
                      ) * y(k,56) + rxt(k,314)*y(k,59) + rxt(k,359)*y(k,17) + rxt(k,365) &
                      *y(k,19) + rxt(k,423)*y(k,42) + rxt(k,431)*y(k,226) + rxt(k,440) &
                      *y(k,233) + rxt(k,460)*y(k,222) + rxt(k,473)*y(k,225) + rxt(k,478) &
                      *y(k,231) + rxt(k,497)*y(k,223) + rxt(k,507)*y(k,253) + rxt(k,512) &
                      *y(k,254) + rxt(k,529)*y(k,236) + rxt(k,542)*y(k,238) + rxt(k,548) &
                      *y(k,240) + rxt(k,564)*y(k,219) + rxt(k,576)*y(k,234) + rxt(k,583) &
                      *y(k,235) + rxt(k,594)*y(k,101) + rxt(k,620)*y(k,260) + rxt(k,626) &
                      *y(k,218) + rxt(k,630)*y(k,220) + rxt(k,635)*y(k,221) + rxt(k,639) &
                      *y(k,224) + rxt(k,644)*y(k,227) + rxt(k,649)*y(k,237) + rxt(k,654) &
                      *y(k,239) + rxt(k,659)*y(k,252) + rxt(k,670)*y(k,257) + rxt(k,678) &
                      *y(k,261) + rxt(k,682)*y(k,262) + rxt(k,698)*y(k,245) + rxt(k,705) &
                      *y(k,255) + rxt(k,711)*y(k,256))
         mat(k,2545) = -(rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,232)
         mat(k,2498) = -rxt(k,218)*y(k,232)
         mat(k,2256) = -rxt(k,219)*y(k,232)
         mat(k,1955) = -rxt(k,226)*y(k,232)
         mat(k,2355) = -rxt(k,248)*y(k,232)
         mat(k,2447) = -rxt(k,255)*y(k,232)
         mat(k,2038) = -rxt(k,266)*y(k,232)
         mat(k,2080) = -(rxt(k,307) + rxt(k,308)) * y(k,232)
         mat(k,1988) = -rxt(k,314)*y(k,232)
         mat(k,1661) = -rxt(k,359)*y(k,232)
         mat(k,2575) = -rxt(k,365)*y(k,232)
         mat(k,1783) = -rxt(k,423)*y(k,232)
         mat(k,1718) = -rxt(k,431)*y(k,232)
         mat(k,529) = -rxt(k,440)*y(k,232)
         mat(k,1131) = -rxt(k,460)*y(k,232)
         mat(k,1600) = -rxt(k,473)*y(k,232)
         mat(k,890) = -rxt(k,478)*y(k,232)
         mat(k,909) = -rxt(k,497)*y(k,232)
         mat(k,866) = -rxt(k,507)*y(k,232)
         mat(k,1331) = -rxt(k,512)*y(k,232)
         mat(k,1517) = -rxt(k,529)*y(k,232)
         mat(k,1562) = -rxt(k,542)*y(k,232)
         mat(k,726) = -rxt(k,548)*y(k,232)
         mat(k,1072) = -rxt(k,564)*y(k,232)
         mat(k,1490) = -rxt(k,576)*y(k,232)
         mat(k,1457) = -rxt(k,583)*y(k,232)
         mat(k,1393) = -rxt(k,594)*y(k,232)
         mat(k,1368) = -rxt(k,620)*y(k,232)
         mat(k,672) = -rxt(k,626)*y(k,232)
         mat(k,598) = -rxt(k,630)*y(k,232)
         mat(k,521) = -rxt(k,635)*y(k,232)
         mat(k,745) = -rxt(k,639)*y(k,232)
         mat(k,835) = -rxt(k,644)*y(k,232)
         mat(k,806) = -rxt(k,649)*y(k,232)
         mat(k,991) = -rxt(k,654)*y(k,232)
         mat(k,546) = -rxt(k,659)*y(k,232)
         mat(k,798) = -rxt(k,670)*y(k,232)
         mat(k,825) = -rxt(k,678)*y(k,232)
         mat(k,622) = -rxt(k,682)*y(k,232)
         mat(k,1312) = -rxt(k,698)*y(k,232)
         mat(k,1290) = -rxt(k,705)*y(k,232)
         mat(k,1265) = -rxt(k,711)*y(k,232)
         mat(k,1021) = .570_r8*rxt(k,688)*y(k,135)
         mat(k,96) = .650_r8*rxt(k,629)*y(k,250)
         mat(k,1661) = mat(k,1661) + rxt(k,358)*y(k,42)
         mat(k,2575) = mat(k,2575) + rxt(k,372)*y(k,250)
         mat(k,258) = .350_r8*rxt(k,454)*y(k,250)
         mat(k,502) = .130_r8*rxt(k,456)*y(k,135)
         mat(k,220) = rxt(k,463)*y(k,250)
         mat(k,1166) = .280_r8*rxt(k,493)*y(k,135)
         mat(k,1783) = mat(k,1783) + rxt(k,358)*y(k,17) + rxt(k,303)*y(k,56) &
                      + rxt(k,424)*y(k,126) + rxt(k,425)*y(k,133)
         mat(k,54) = rxt(k,470)*y(k,250)
         mat(k,782) = rxt(k,434)*y(k,250)
         mat(k,2080) = mat(k,2080) + rxt(k,303)*y(k,42) + rxt(k,306)*y(k,79)
         mat(k,1988) = mat(k,1988) + rxt(k,310)*y(k,226) + rxt(k,322)*y(k,250)
         mat(k,1236) = rxt(k,437)*y(k,250)
         mat(k,124) = .730_r8*rxt(k,643)*y(k,250)
         mat(k,341) = .500_r8*rxt(k,746)*y(k,250)
         mat(k,1194) = rxt(k,483)*y(k,250)
         mat(k,916) = rxt(k,484)*y(k,250)
         mat(k,2545) = mat(k,2545) + rxt(k,217)*y(k,134)
         mat(k,577) = rxt(k,306)*y(k,56) + rxt(k,211)*y(k,133) + rxt(k,225)*y(k,250)
         mat(k,202) = rxt(k,438)*y(k,250)
         mat(k,777) = rxt(k,439)*y(k,250)
         mat(k,1220) = rxt(k,526)*y(k,250)
         mat(k,1243) = rxt(k,505)*y(k,250)
         mat(k,978) = .370_r8*rxt(k,600)*y(k,135)
         mat(k,587) = .300_r8*rxt(k,589)*y(k,250)
         mat(k,537) = rxt(k,590)*y(k,250)
         mat(k,1393) = mat(k,1393) + rxt(k,595)*y(k,124) + rxt(k,596)*y(k,126) &
                      + rxt(k,592)*y(k,225) + 1.200_r8*rxt(k,593)*y(k,226)
         mat(k,373) = rxt(k,599)*y(k,250)
         mat(k,1409) = .140_r8*rxt(k,536)*y(k,135)
         mat(k,305) = .200_r8*rxt(k,538)*y(k,250)
         mat(k,479) = .500_r8*rxt(k,553)*y(k,250)
         mat(k,1053) = .570_r8*rxt(k,693)*y(k,135)
         mat(k,1540) = .280_r8*rxt(k,554)*y(k,135)
         mat(k,348) = rxt(k,605)*y(k,250)
         mat(k,1186) = rxt(k,606)*y(k,250)
         mat(k,2447) = mat(k,2447) + rxt(k,595)*y(k,101) + rxt(k,565)*y(k,219) &
                      + rxt(k,631)*y(k,220) + rxt(k,637)*y(k,221) + rxt(k,461) &
                      *y(k,222) + rxt(k,498)*y(k,223) + rxt(k,432)*y(k,226) &
                      + .170_r8*rxt(k,645)*y(k,227) + rxt(k,522)*y(k,229) &
                      + .250_r8*rxt(k,479)*y(k,231) + rxt(k,442)*y(k,233) &
                      + .920_r8*rxt(k,577)*y(k,234) + .920_r8*rxt(k,585)*y(k,235) &
                      + .470_r8*rxt(k,531)*y(k,236) + .400_r8*rxt(k,650)*y(k,237) &
                      + .830_r8*rxt(k,655)*y(k,239) + rxt(k,660)*y(k,252) + rxt(k,508) &
                      *y(k,253) + .900_r8*rxt(k,706)*y(k,255) + .800_r8*rxt(k,712) &
                      *y(k,256) + rxt(k,671)*y(k,257) + rxt(k,621)*y(k,260) &
                      + rxt(k,679)*y(k,261) + rxt(k,683)*y(k,262)
         mat(k,2355) = mat(k,2355) + rxt(k,424)*y(k,42) + rxt(k,596)*y(k,101) &
                      + rxt(k,578)*y(k,234) + rxt(k,586)*y(k,235) + .470_r8*rxt(k,530) &
                      *y(k,236) + rxt(k,251)*y(k,250) + rxt(k,622)*y(k,260)
         mat(k,2498) = mat(k,2498) + rxt(k,425)*y(k,42) + rxt(k,211)*y(k,79)
         mat(k,1756) = rxt(k,217)*y(k,76) + rxt(k,482)*y(k,230)
         mat(k,2256) = mat(k,2256) + .570_r8*rxt(k,688)*y(k,6) + .130_r8*rxt(k,456) &
                      *y(k,25) + .280_r8*rxt(k,493)*y(k,29) + .370_r8*rxt(k,600) &
                      *y(k,98) + .140_r8*rxt(k,536)*y(k,105) + .570_r8*rxt(k,693) &
                      *y(k,110) + .280_r8*rxt(k,554)*y(k,111) + rxt(k,228)*y(k,250)
         mat(k,105) = .800_r8*rxt(k,662)*y(k,250)
         mat(k,921) = rxt(k,740)*y(k,250)
         mat(k,1232) = .200_r8*rxt(k,718)*y(k,250)
         mat(k,119) = .280_r8*rxt(k,674)*y(k,250)
         mat(k,166) = .380_r8*rxt(k,677)*y(k,250)
         mat(k,171) = .630_r8*rxt(k,685)*y(k,250)
         mat(k,1072) = mat(k,1072) + rxt(k,565)*y(k,124)
         mat(k,598) = mat(k,598) + rxt(k,631)*y(k,124)
         mat(k,521) = mat(k,521) + rxt(k,637)*y(k,124)
         mat(k,1131) = mat(k,1131) + rxt(k,461)*y(k,124) + 2.400_r8*rxt(k,458) &
                      *y(k,222) + rxt(k,459)*y(k,226)
         mat(k,909) = mat(k,909) + rxt(k,498)*y(k,124) + rxt(k,496)*y(k,226)
         mat(k,1600) = mat(k,1600) + rxt(k,592)*y(k,101) + .900_r8*rxt(k,472)*y(k,226) &
                      + rxt(k,574)*y(k,234) + rxt(k,581)*y(k,235) + .470_r8*rxt(k,527) &
                      *y(k,236) + rxt(k,618)*y(k,260)
         mat(k,1718) = mat(k,1718) + rxt(k,310)*y(k,59) + 1.200_r8*rxt(k,593)*y(k,101) &
                      + rxt(k,432)*y(k,124) + rxt(k,459)*y(k,222) + rxt(k,496) &
                      *y(k,223) + .900_r8*rxt(k,472)*y(k,225) + 4.000_r8*rxt(k,429) &
                      *y(k,226) + rxt(k,575)*y(k,234) + rxt(k,582)*y(k,235) &
                      + .730_r8*rxt(k,528)*y(k,236) + rxt(k,541)*y(k,238) &
                      + .500_r8*rxt(k,697)*y(k,245) + .300_r8*rxt(k,511)*y(k,254) &
                      + rxt(k,704)*y(k,255) + rxt(k,710)*y(k,256) + .800_r8*rxt(k,619) &
                      *y(k,260)
         mat(k,835) = mat(k,835) + .170_r8*rxt(k,645)*y(k,124) + .070_r8*rxt(k,644) &
                      *y(k,232)
         mat(k,707) = rxt(k,522)*y(k,124)
         mat(k,405) = rxt(k,482)*y(k,134)
         mat(k,890) = mat(k,890) + .250_r8*rxt(k,479)*y(k,124)
         mat(k,2187) = mat(k,2187) + .070_r8*rxt(k,644)*y(k,227) + .160_r8*rxt(k,649) &
                      *y(k,237) + .330_r8*rxt(k,654)*y(k,239)
         mat(k,529) = mat(k,529) + rxt(k,442)*y(k,124)
         mat(k,1490) = mat(k,1490) + .920_r8*rxt(k,577)*y(k,124) + rxt(k,578)*y(k,126) &
                      + rxt(k,574)*y(k,225) + rxt(k,575)*y(k,226)
         mat(k,1457) = mat(k,1457) + .920_r8*rxt(k,585)*y(k,124) + rxt(k,586)*y(k,126) &
                      + rxt(k,581)*y(k,225) + rxt(k,582)*y(k,226)
         mat(k,1517) = mat(k,1517) + .470_r8*rxt(k,531)*y(k,124) + .470_r8*rxt(k,530) &
                      *y(k,126) + .470_r8*rxt(k,527)*y(k,225) + .730_r8*rxt(k,528) &
                      *y(k,226)
         mat(k,806) = mat(k,806) + .400_r8*rxt(k,650)*y(k,124) + .160_r8*rxt(k,649) &
                      *y(k,232)
         mat(k,1562) = mat(k,1562) + rxt(k,541)*y(k,226)
         mat(k,991) = mat(k,991) + .830_r8*rxt(k,655)*y(k,124) + .330_r8*rxt(k,654) &
                      *y(k,232)
         mat(k,1312) = mat(k,1312) + .500_r8*rxt(k,697)*y(k,226)
         mat(k,1955) = mat(k,1955) + .650_r8*rxt(k,629)*y(k,7) + rxt(k,372)*y(k,19) &
                      + .350_r8*rxt(k,454)*y(k,24) + rxt(k,463)*y(k,26) + rxt(k,470) &
                      *y(k,47) + rxt(k,434)*y(k,52) + rxt(k,322)*y(k,59) + rxt(k,437) &
                      *y(k,62) + .730_r8*rxt(k,643)*y(k,66) + .500_r8*rxt(k,746) &
                      *y(k,67) + rxt(k,483)*y(k,74) + rxt(k,484)*y(k,75) + rxt(k,225) &
                      *y(k,79) + rxt(k,438)*y(k,86) + rxt(k,439)*y(k,87) + rxt(k,526) &
                      *y(k,93) + rxt(k,505)*y(k,95) + .300_r8*rxt(k,589)*y(k,99) &
                      + rxt(k,590)*y(k,100) + rxt(k,599)*y(k,102) + .200_r8*rxt(k,538) &
                      *y(k,106) + .500_r8*rxt(k,553)*y(k,109) + rxt(k,605)*y(k,115) &
                      + rxt(k,606)*y(k,116) + rxt(k,251)*y(k,126) + rxt(k,228) &
                      *y(k,135) + .800_r8*rxt(k,662)*y(k,142) + rxt(k,740)*y(k,151) &
                      + .200_r8*rxt(k,718)*y(k,178) + .280_r8*rxt(k,674)*y(k,180) &
                      + .380_r8*rxt(k,677)*y(k,211) + .630_r8*rxt(k,685)*y(k,213)
         mat(k,546) = mat(k,546) + rxt(k,660)*y(k,124)
         mat(k,866) = mat(k,866) + rxt(k,508)*y(k,124)
         mat(k,1331) = mat(k,1331) + .300_r8*rxt(k,511)*y(k,226)
         mat(k,1290) = mat(k,1290) + .900_r8*rxt(k,706)*y(k,124) + rxt(k,704)*y(k,226)
         mat(k,1265) = mat(k,1265) + .800_r8*rxt(k,712)*y(k,124) + rxt(k,710)*y(k,226)
         mat(k,798) = mat(k,798) + rxt(k,671)*y(k,124)
         mat(k,1368) = mat(k,1368) + rxt(k,621)*y(k,124) + rxt(k,622)*y(k,126) &
                      + rxt(k,618)*y(k,225) + .800_r8*rxt(k,619)*y(k,226)
         mat(k,825) = mat(k,825) + rxt(k,679)*y(k,124)
         mat(k,622) = mat(k,622) + rxt(k,683)*y(k,124)
      end do
      end subroutine nlnmat09
      subroutine nlnmat10( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,523) = -(rxt(k,440)*y(k,232) + rxt(k,442)*y(k,124))
         mat(k,2121) = -rxt(k,440)*y(k,233)
         mat(k,2381) = -rxt(k,442)*y(k,233)
         mat(k,1765) = rxt(k,423)*y(k,232)
         mat(k,2121) = mat(k,2121) + rxt(k,423)*y(k,42)
         mat(k,1481) = -(rxt(k,574)*y(k,225) + rxt(k,575)*y(k,226) + rxt(k,576) &
                      *y(k,232) + rxt(k,577)*y(k,124) + rxt(k,578)*y(k,126))
         mat(k,1590) = -rxt(k,574)*y(k,234)
         mat(k,1706) = -rxt(k,575)*y(k,234)
         mat(k,2173) = -rxt(k,576)*y(k,234)
         mat(k,2434) = -rxt(k,577)*y(k,234)
         mat(k,2341) = -rxt(k,578)*y(k,234)
         mat(k,972) = .600_r8*rxt(k,601)*y(k,250)
         mat(k,1940) = .600_r8*rxt(k,601)*y(k,98)
         mat(k,1448) = -(rxt(k,581)*y(k,225) + rxt(k,582)*y(k,226) + rxt(k,583) &
                      *y(k,232) + rxt(k,585)*y(k,124) + rxt(k,586)*y(k,126))
         mat(k,1589) = -rxt(k,581)*y(k,235)
         mat(k,1705) = -rxt(k,582)*y(k,235)
         mat(k,2172) = -rxt(k,583)*y(k,235)
         mat(k,2433) = -rxt(k,585)*y(k,235)
         mat(k,2340) = -rxt(k,586)*y(k,235)
         mat(k,971) = .400_r8*rxt(k,601)*y(k,250)
         mat(k,1939) = .400_r8*rxt(k,601)*y(k,98)
         mat(k,1510) = -(rxt(k,527)*y(k,225) + rxt(k,528)*y(k,226) + rxt(k,529) &
                      *y(k,232) + rxt(k,530)*y(k,126) + (rxt(k,531) + rxt(k,532) &
                      ) * y(k,124))
         mat(k,1591) = -rxt(k,527)*y(k,236)
         mat(k,1707) = -rxt(k,528)*y(k,236)
         mat(k,2174) = -rxt(k,529)*y(k,236)
         mat(k,2342) = -rxt(k,530)*y(k,236)
         mat(k,2435) = -(rxt(k,531) + rxt(k,532)) * y(k,236)
         mat(k,1404) = .500_r8*rxt(k,537)*y(k,250)
         mat(k,302) = .200_r8*rxt(k,538)*y(k,250)
         mat(k,1531) = rxt(k,555)*y(k,250)
         mat(k,1941) = .500_r8*rxt(k,537)*y(k,105) + .200_r8*rxt(k,538)*y(k,106) &
                      + rxt(k,555)*y(k,111)
         mat(k,800) = -(rxt(k,649)*y(k,232) + rxt(k,650)*y(k,124) + rxt(k,651) &
                      *y(k,125))
         mat(k,2139) = -rxt(k,649)*y(k,237)
         mat(k,2400) = -rxt(k,650)*y(k,237)
         mat(k,2012) = -rxt(k,651)*y(k,237)
         mat(k,1556) = -(rxt(k,540)*y(k,225) + rxt(k,541)*y(k,226) + rxt(k,542) &
                      *y(k,232) + 4._r8*rxt(k,543)*y(k,238) + rxt(k,544)*y(k,124) &
                      + rxt(k,545)*y(k,126) + rxt(k,557)*y(k,125))
         mat(k,1593) = -rxt(k,540)*y(k,238)
         mat(k,1709) = -rxt(k,541)*y(k,238)
         mat(k,2176) = -rxt(k,542)*y(k,238)
         mat(k,2437) = -rxt(k,544)*y(k,238)
         mat(k,2344) = -rxt(k,545)*y(k,238)
         mat(k,2026) = -rxt(k,557)*y(k,238)
         mat(k,1405) = .500_r8*rxt(k,537)*y(k,250)
         mat(k,303) = .500_r8*rxt(k,538)*y(k,250)
         mat(k,1943) = .500_r8*rxt(k,537)*y(k,105) + .500_r8*rxt(k,538)*y(k,106)
         mat(k,983) = -(rxt(k,654)*y(k,232) + rxt(k,655)*y(k,124) + rxt(k,656) &
                      *y(k,125))
         mat(k,2149) = -rxt(k,654)*y(k,239)
         mat(k,2409) = -rxt(k,655)*y(k,239)
         mat(k,2018) = -rxt(k,656)*y(k,239)
         mat(k,719) = -(rxt(k,548)*y(k,232) + rxt(k,549)*y(k,124))
         mat(k,2134) = -rxt(k,548)*y(k,240)
         mat(k,2393) = -rxt(k,549)*y(k,240)
         mat(k,460) = rxt(k,551)*y(k,250)
         mat(k,317) = rxt(k,552)*y(k,250)
         mat(k,1891) = rxt(k,551)*y(k,107) + rxt(k,552)*y(k,108)
         mat(k,464) = -(rxt(k,235)*y(k,133) + rxt(k,236)*y(k,134))
         mat(k,2465) = -rxt(k,235)*y(k,241)
         mat(k,1732) = -rxt(k,236)*y(k,241)
         mat(k,2465) = mat(k,2465) + rxt(k,827)*y(k,242)
         mat(k,764) = .900_r8*rxt(k,825)*y(k,242) + .800_r8*rxt(k,823)*y(k,243)
         mat(k,641) = rxt(k,827)*y(k,133) + .900_r8*rxt(k,825)*y(k,228)
         mat(k,756) = .800_r8*rxt(k,823)*y(k,228)
         mat(k,643) = -(rxt(k,825)*y(k,228) + rxt(k,826)*y(k,134) + (rxt(k,827) &
                      + rxt(k,828)) * y(k,133))
         mat(k,765) = -rxt(k,825)*y(k,242)
         mat(k,1736) = -rxt(k,826)*y(k,242)
         mat(k,2471) = -(rxt(k,827) + rxt(k,828)) * y(k,242)
         mat(k,757) = -(rxt(k,823)*y(k,228))
         mat(k,767) = -rxt(k,823)*y(k,243)
         mat(k,928) = rxt(k,832)*y(k,249)
         mat(k,2396) = rxt(k,834)*y(k,249)
         mat(k,2475) = rxt(k,827)*y(k,242)
         mat(k,1739) = rxt(k,831)*y(k,244)
         mat(k,645) = rxt(k,827)*y(k,133)
         mat(k,449) = rxt(k,831)*y(k,134)
         mat(k,749) = rxt(k,832)*y(k,112) + rxt(k,834)*y(k,124)
         mat(k,446) = -(rxt(k,829)*y(k,133) + (rxt(k,830) + rxt(k,831)) * y(k,134))
         mat(k,2464) = -rxt(k,829)*y(k,244)
         mat(k,1731) = -(rxt(k,830) + rxt(k,831)) * y(k,244)
         mat(k,1305) = -(rxt(k,697)*y(k,226) + rxt(k,698)*y(k,232) + rxt(k,699) &
                      *y(k,124) + rxt(k,700)*y(k,126))
         mat(k,1699) = -rxt(k,697)*y(k,245)
         mat(k,2165) = -rxt(k,698)*y(k,245)
         mat(k,2427) = -rxt(k,699)*y(k,245)
         mat(k,2334) = -rxt(k,700)*y(k,245)
         mat(k,1015) = rxt(k,687)*y(k,126)
         mat(k,1047) = rxt(k,692)*y(k,126)
         mat(k,2334) = mat(k,2334) + rxt(k,687)*y(k,6) + rxt(k,692)*y(k,110) &
                      + .500_r8*rxt(k,715)*y(k,177)
         mat(k,352) = rxt(k,703)*y(k,250)
         mat(k,1141) = .500_r8*rxt(k,715)*y(k,126)
         mat(k,1932) = rxt(k,703)*y(k,128)
         mat(k,2290) = -(rxt(k,182)*y(k,77) + rxt(k,183)*y(k,263) + (rxt(k,185) &
                      + rxt(k,186)) * y(k,134) + rxt(k,187)*y(k,135) + (rxt(k,342) &
                      + rxt(k,343)) * y(k,85) + (rxt(k,384) + rxt(k,385)) * y(k,81) &
                      + rxt(k,397)*y(k,64) + rxt(k,398)*y(k,65) + rxt(k,447)*y(k,86))
         mat(k,1346) = -rxt(k,182)*y(k,246)
         mat(k,2608) = -rxt(k,183)*y(k,246)
         mat(k,1758) = -(rxt(k,185) + rxt(k,186)) * y(k,246)
         mat(k,2258) = -rxt(k,187)*y(k,246)
         mat(k,1646) = -(rxt(k,342) + rxt(k,343)) * y(k,246)
         mat(k,851) = -(rxt(k,384) + rxt(k,385)) * y(k,246)
         mat(k,146) = -rxt(k,397)*y(k,246)
         mat(k,211) = -rxt(k,398)*y(k,246)
         mat(k,203) = -rxt(k,447)*y(k,246)
         mat(k,1758) = mat(k,1758) + rxt(k,236)*y(k,241)
         mat(k,772) = .850_r8*rxt(k,824)*y(k,249)
         mat(k,467) = rxt(k,236)*y(k,134)
         mat(k,753) = .850_r8*rxt(k,824)*y(k,228)
         mat(k,78) = -(rxt(k,190)*y(k,133) + rxt(k,191)*y(k,134))
         mat(k,2457) = -rxt(k,190)*y(k,247)
         mat(k,1725) = -rxt(k,191)*y(k,247)
         mat(k,2457) = mat(k,2457) + rxt(k,194)*y(k,248)
         mat(k,1725) = mat(k,1725) + rxt(k,195)*y(k,248)
         mat(k,2202) = rxt(k,196)*y(k,248)
         mat(k,80) = rxt(k,194)*y(k,133) + rxt(k,195)*y(k,134) + rxt(k,196)*y(k,135)
         mat(k,81) = -(rxt(k,194)*y(k,133) + rxt(k,195)*y(k,134) + rxt(k,196)*y(k,135))
         mat(k,2458) = -rxt(k,194)*y(k,248)
         mat(k,1726) = -rxt(k,195)*y(k,248)
         mat(k,2203) = -rxt(k,196)*y(k,248)
         mat(k,1726) = mat(k,1726) + rxt(k,185)*y(k,246)
         mat(k,2266) = rxt(k,185)*y(k,134)
         mat(k,748) = -(rxt(k,824)*y(k,228) + rxt(k,832)*y(k,112) + rxt(k,834) &
                      *y(k,124))
         mat(k,766) = -rxt(k,824)*y(k,249)
         mat(k,927) = -rxt(k,832)*y(k,249)
         mat(k,2395) = -rxt(k,834)*y(k,249)
         mat(k,1738) = rxt(k,826)*y(k,242) + rxt(k,830)*y(k,244) + rxt(k,837)*y(k,251)
         mat(k,644) = rxt(k,826)*y(k,134)
         mat(k,448) = rxt(k,830)*y(k,134)
         mat(k,608) = rxt(k,837)*y(k,134)
         mat(k,1951) = -(rxt(k,224)*y(k,77) + rxt(k,225)*y(k,79) + rxt(k,226)*y(k,232) &
                      + rxt(k,227)*y(k,133) + rxt(k,228)*y(k,135) + (4._r8*rxt(k,229) &
                      + 4._r8*rxt(k,230)) * y(k,250) + rxt(k,234)*y(k,90) + rxt(k,251) &
                      *y(k,126) + rxt(k,254)*y(k,112) + rxt(k,268)*y(k,125) + rxt(k,274) &
                      *y(k,89) + rxt(k,319)*y(k,60) + (rxt(k,322) + rxt(k,323) &
                      ) * y(k,59) + rxt(k,329)*y(k,85) + rxt(k,333)*y(k,92) + rxt(k,372) &
                      *y(k,19) + rxt(k,377)*y(k,81) + rxt(k,426)*y(k,42) + rxt(k,434) &
                      *y(k,52) + rxt(k,435)*y(k,53) + (rxt(k,437) + rxt(k,448) &
                      ) * y(k,62) + rxt(k,438)*y(k,86) + rxt(k,439)*y(k,87) + rxt(k,454) &
                      *y(k,24) + rxt(k,463)*y(k,26) + rxt(k,464)*y(k,27) + rxt(k,466) &
                      *y(k,28) + rxt(k,468)*y(k,45) + rxt(k,470)*y(k,47) + rxt(k,476) &
                      *y(k,50) + rxt(k,477)*y(k,51) + rxt(k,483)*y(k,74) + rxt(k,484) &
                      *y(k,75) + rxt(k,485)*y(k,139) + rxt(k,486)*y(k,25) + rxt(k,500) &
                      *y(k,30) + rxt(k,501)*y(k,31) + rxt(k,503)*y(k,49) + rxt(k,505) &
                      *y(k,95) + rxt(k,506)*y(k,127) + rxt(k,510)*y(k,146) + rxt(k,515) &
                      *y(k,147) + rxt(k,516)*y(k,29) + rxt(k,517)*y(k,48) + rxt(k,520) &
                      *y(k,16) + rxt(k,526)*y(k,93) + rxt(k,537)*y(k,105) + rxt(k,538) &
                      *y(k,106) + rxt(k,551)*y(k,107) + rxt(k,552)*y(k,108) + rxt(k,553) &
                      *y(k,109) + rxt(k,555)*y(k,111) + rxt(k,563)*y(k,1) + rxt(k,569) &
                      *y(k,2) + rxt(k,570)*y(k,15) + rxt(k,571)*y(k,94) + rxt(k,572) &
                      *y(k,96) + rxt(k,573)*y(k,97) + rxt(k,589)*y(k,99) + rxt(k,590) &
                      *y(k,100) + rxt(k,599)*y(k,102) + rxt(k,601)*y(k,98) + rxt(k,602) &
                      *y(k,103) + rxt(k,605)*y(k,115) + rxt(k,606)*y(k,116) + rxt(k,625) &
                      *y(k,207) + rxt(k,629)*y(k,7) + rxt(k,633)*y(k,8) + rxt(k,634) &
                      *y(k,22) + rxt(k,636)*y(k,23) + rxt(k,642)*y(k,32) + rxt(k,643) &
                      *y(k,66) + rxt(k,662)*y(k,142) + rxt(k,665)*y(k,143) + rxt(k,673) &
                      *y(k,179) + rxt(k,674)*y(k,180) + rxt(k,677)*y(k,211) + rxt(k,681) &
                      *y(k,212) + rxt(k,685)*y(k,213) + rxt(k,686)*y(k,214) + rxt(k,689) &
                      *y(k,6) + rxt(k,694)*y(k,110) + rxt(k,703)*y(k,128) + rxt(k,708) &
                      *y(k,174) + rxt(k,709)*y(k,175) + rxt(k,714)*y(k,176) + rxt(k,716) &
                      *y(k,177) + rxt(k,718)*y(k,178) + rxt(k,726)*y(k,137) + rxt(k,731) &
                      *y(k,148) + rxt(k,736)*y(k,150) + rxt(k,740)*y(k,151) + (rxt(k,743) &
                      + rxt(k,746)) * y(k,67) + rxt(k,745)*y(k,120))
         mat(k,1344) = -rxt(k,224)*y(k,250)
         mat(k,575) = -rxt(k,225)*y(k,250)
         mat(k,2183) = -rxt(k,226)*y(k,250)
         mat(k,2494) = -rxt(k,227)*y(k,250)
         mat(k,2252) = -rxt(k,228)*y(k,250)
         mat(k,418) = -rxt(k,234)*y(k,250)
         mat(k,2351) = -rxt(k,251)*y(k,250)
         mat(k,935) = -rxt(k,254)*y(k,250)
         mat(k,2034) = -rxt(k,268)*y(k,250)
         mat(k,2519) = -rxt(k,274)*y(k,250)
         mat(k,946) = -rxt(k,319)*y(k,250)
         mat(k,1984) = -(rxt(k,322) + rxt(k,323)) * y(k,250)
         mat(k,1642) = -rxt(k,329)*y(k,250)
         mat(k,841) = -rxt(k,333)*y(k,250)
         mat(k,2571) = -rxt(k,372)*y(k,250)
         mat(k,850) = -rxt(k,377)*y(k,250)
         mat(k,1779) = -rxt(k,426)*y(k,250)
         mat(k,781) = -rxt(k,434)*y(k,250)
         mat(k,389) = -rxt(k,435)*y(k,250)
         mat(k,1235) = -(rxt(k,437) + rxt(k,448)) * y(k,250)
         mat(k,201) = -rxt(k,438)*y(k,250)
         mat(k,776) = -rxt(k,439)*y(k,250)
         mat(k,257) = -rxt(k,454)*y(k,250)
         mat(k,219) = -rxt(k,463)*y(k,250)
         mat(k,295) = -rxt(k,464)*y(k,250)
         mat(k,262) = -rxt(k,466)*y(k,250)
         mat(k,1202) = -rxt(k,468)*y(k,250)
         mat(k,53) = -rxt(k,470)*y(k,250)
         mat(k,471) = -rxt(k,476)*y(k,250)
         mat(k,457) = -rxt(k,477)*y(k,250)
         mat(k,1193) = -rxt(k,483)*y(k,250)
         mat(k,915) = -rxt(k,484)*y(k,250)
         mat(k,443) = -rxt(k,485)*y(k,250)
         mat(k,501) = -rxt(k,486)*y(k,250)
         mat(k,377) = -rxt(k,500)*y(k,250)
         mat(k,59) = -rxt(k,501)*y(k,250)
         mat(k,1419) = -rxt(k,503)*y(k,250)
         mat(k,1242) = -rxt(k,505)*y(k,250)
         mat(k,896) = -rxt(k,506)*y(k,250)
         mat(k,486) = -rxt(k,510)*y(k,250)
         mat(k,360) = -rxt(k,515)*y(k,250)
         mat(k,1164) = -rxt(k,516)*y(k,250)
         mat(k,1078) = -rxt(k,517)*y(k,250)
         mat(k,413) = -rxt(k,520)*y(k,250)
         mat(k,1218) = -rxt(k,526)*y(k,250)
         mat(k,1408) = -rxt(k,537)*y(k,250)
         mat(k,304) = -rxt(k,538)*y(k,250)
         mat(k,463) = -rxt(k,551)*y(k,250)
         mat(k,320) = -rxt(k,552)*y(k,250)
         mat(k,477) = -rxt(k,553)*y(k,250)
         mat(k,1538) = -rxt(k,555)*y(k,250)
         mat(k,638) = -rxt(k,563)*y(k,250)
         mat(k,662) = -rxt(k,569)*y(k,250)
         mat(k,182) = -rxt(k,570)*y(k,250)
         mat(k,175) = -rxt(k,571)*y(k,250)
         mat(k,312) = -rxt(k,572)*y(k,250)
         mat(k,71) = -rxt(k,573)*y(k,250)
         mat(k,585) = -rxt(k,589)*y(k,250)
         mat(k,536) = -rxt(k,590)*y(k,250)
         mat(k,371) = -rxt(k,599)*y(k,250)
         mat(k,977) = -rxt(k,601)*y(k,250)
         mat(k,680) = -rxt(k,602)*y(k,250)
         mat(k,347) = -rxt(k,605)*y(k,250)
         mat(k,1184) = -rxt(k,606)*y(k,250)
         mat(k,149) = -rxt(k,625)*y(k,250)
         mat(k,95) = -rxt(k,629)*y(k,250)
         mat(k,366) = -rxt(k,633)*y(k,250)
         mat(k,188) = -rxt(k,634)*y(k,250)
         mat(k,290) = -rxt(k,636)*y(k,250)
         mat(k,224) = -rxt(k,642)*y(k,250)
         mat(k,123) = -rxt(k,643)*y(k,250)
         mat(k,104) = -rxt(k,662)*y(k,250)
         mat(k,284) = -rxt(k,665)*y(k,250)
         mat(k,566) = -rxt(k,673)*y(k,250)
         mat(k,118) = -rxt(k,674)*y(k,250)
         mat(k,165) = -rxt(k,677)*y(k,250)
         mat(k,695) = -rxt(k,681)*y(k,250)
         mat(k,170) = -rxt(k,685)*y(k,250)
         mat(k,396) = -rxt(k,686)*y(k,250)
         mat(k,1020) = -rxt(k,689)*y(k,250)
         mat(k,1052) = -rxt(k,694)*y(k,250)
         mat(k,353) = -rxt(k,703)*y(k,250)
         mat(k,555) = -rxt(k,708)*y(k,250)
         mat(k,626) = -rxt(k,709)*y(k,250)
         mat(k,428) = -rxt(k,714)*y(k,250)
         mat(k,1142) = -rxt(k,716)*y(k,250)
         mat(k,1231) = -rxt(k,718)*y(k,250)
         mat(k,326) = -rxt(k,726)*y(k,250)
         mat(k,714) = -rxt(k,731)*y(k,250)
         mat(k,1618) = -rxt(k,736)*y(k,250)
         mat(k,920) = -rxt(k,740)*y(k,250)
         mat(k,340) = -(rxt(k,743) + rxt(k,746)) * y(k,250)
         mat(k,50) = -rxt(k,745)*y(k,250)
         mat(k,1020) = mat(k,1020) + .630_r8*rxt(k,688)*y(k,135)
         mat(k,257) = mat(k,257) + .650_r8*rxt(k,454)*y(k,250)
         mat(k,501) = mat(k,501) + .130_r8*rxt(k,456)*y(k,135)
         mat(k,295) = mat(k,295) + .500_r8*rxt(k,464)*y(k,250)
         mat(k,1164) = mat(k,1164) + .360_r8*rxt(k,493)*y(k,135)
         mat(k,1779) = mat(k,1779) + rxt(k,425)*y(k,133)
         mat(k,389) = mat(k,389) + .300_r8*rxt(k,435)*y(k,250)
         mat(k,2076) = rxt(k,308)*y(k,232)
         mat(k,875) = rxt(k,394)*y(k,263)
         mat(k,2541) = rxt(k,222)*y(k,135) + 2.000_r8*rxt(k,215)*y(k,232)
         mat(k,1344) = mat(k,1344) + rxt(k,210)*y(k,133) + rxt(k,182)*y(k,246)
         mat(k,575) = mat(k,575) + rxt(k,211)*y(k,133)
         mat(k,850) = mat(k,850) + rxt(k,376)*y(k,133) + rxt(k,384)*y(k,246)
         mat(k,1642) = mat(k,1642) + rxt(k,328)*y(k,133) + rxt(k,342)*y(k,246)
         mat(k,201) = mat(k,201) + rxt(k,447)*y(k,246)
         mat(k,733) = rxt(k,379)*y(k,133)
         mat(k,841) = mat(k,841) + rxt(k,332)*y(k,133)
         mat(k,977) = mat(k,977) + .320_r8*rxt(k,600)*y(k,135)
         mat(k,680) = mat(k,680) + .600_r8*rxt(k,602)*y(k,250)
         mat(k,1408) = mat(k,1408) + .240_r8*rxt(k,536)*y(k,135)
         mat(k,304) = mat(k,304) + .100_r8*rxt(k,538)*y(k,250)
         mat(k,1052) = mat(k,1052) + .630_r8*rxt(k,693)*y(k,135)
         mat(k,1538) = mat(k,1538) + .360_r8*rxt(k,554)*y(k,135)
         mat(k,2443) = rxt(k,255)*y(k,232)
         mat(k,2351) = mat(k,2351) + rxt(k,248)*y(k,232)
         mat(k,2494) = mat(k,2494) + rxt(k,425)*y(k,42) + rxt(k,210)*y(k,77) &
                      + rxt(k,211)*y(k,79) + rxt(k,376)*y(k,81) + rxt(k,328)*y(k,85) &
                      + rxt(k,379)*y(k,91) + rxt(k,332)*y(k,92) + rxt(k,218)*y(k,232)
         mat(k,2252) = mat(k,2252) + .630_r8*rxt(k,688)*y(k,6) + .130_r8*rxt(k,456) &
                      *y(k,25) + .360_r8*rxt(k,493)*y(k,29) + rxt(k,222)*y(k,76) &
                      + .320_r8*rxt(k,600)*y(k,98) + .240_r8*rxt(k,536)*y(k,105) &
                      + .630_r8*rxt(k,693)*y(k,110) + .360_r8*rxt(k,554)*y(k,111) &
                      + rxt(k,219)*y(k,232)
         mat(k,486) = mat(k,486) + .500_r8*rxt(k,510)*y(k,250)
         mat(k,149) = mat(k,149) + .500_r8*rxt(k,625)*y(k,250)
         mat(k,670) = .400_r8*rxt(k,626)*y(k,232)
         mat(k,1598) = .450_r8*rxt(k,473)*y(k,232)
         mat(k,833) = .400_r8*rxt(k,644)*y(k,232)
         mat(k,2183) = mat(k,2183) + rxt(k,308)*y(k,56) + 2.000_r8*rxt(k,215)*y(k,76) &
                      + rxt(k,255)*y(k,124) + rxt(k,248)*y(k,126) + rxt(k,218) &
                      *y(k,133) + rxt(k,219)*y(k,135) + .400_r8*rxt(k,626)*y(k,218) &
                      + .450_r8*rxt(k,473)*y(k,225) + .400_r8*rxt(k,644)*y(k,227) &
                      + .450_r8*rxt(k,542)*y(k,238) + .400_r8*rxt(k,654)*y(k,239) &
                      + .200_r8*rxt(k,548)*y(k,240) + .150_r8*rxt(k,512)*y(k,254)
         mat(k,1560) = .450_r8*rxt(k,542)*y(k,232)
         mat(k,989) = .400_r8*rxt(k,654)*y(k,232)
         mat(k,724) = .200_r8*rxt(k,548)*y(k,232)
         mat(k,2284) = rxt(k,182)*y(k,77) + rxt(k,384)*y(k,81) + rxt(k,342)*y(k,85) &
                      + rxt(k,447)*y(k,86) + 2.000_r8*rxt(k,183)*y(k,263)
         mat(k,1951) = mat(k,1951) + .650_r8*rxt(k,454)*y(k,24) + .500_r8*rxt(k,464) &
                      *y(k,27) + .300_r8*rxt(k,435)*y(k,53) + .600_r8*rxt(k,602) &
                      *y(k,103) + .100_r8*rxt(k,538)*y(k,106) + .500_r8*rxt(k,510) &
                      *y(k,146) + .500_r8*rxt(k,625)*y(k,207)
         mat(k,1329) = .150_r8*rxt(k,512)*y(k,232)
         mat(k,2602) = rxt(k,394)*y(k,73) + 2.000_r8*rxt(k,183)*y(k,246)
      end do
      end subroutine nlnmat10
      subroutine nlnmat11( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,607) = -(rxt(k,837)*y(k,134))
         mat(k,1735) = -rxt(k,837)*y(k,251)
         mat(k,2470) = rxt(k,828)*y(k,242) + rxt(k,829)*y(k,244)
         mat(k,642) = rxt(k,828)*y(k,133)
         mat(k,447) = rxt(k,829)*y(k,133)
         mat(k,540) = -(rxt(k,659)*y(k,232) + rxt(k,660)*y(k,124))
         mat(k,2122) = -rxt(k,659)*y(k,252)
         mat(k,2383) = -rxt(k,660)*y(k,252)
         mat(k,121) = .200_r8*rxt(k,643)*y(k,250)
         mat(k,102) = .140_r8*rxt(k,662)*y(k,250)
         mat(k,282) = rxt(k,665)*y(k,250)
         mat(k,1874) = .200_r8*rxt(k,643)*y(k,66) + .140_r8*rxt(k,662)*y(k,142) &
                      + rxt(k,665)*y(k,143)
         mat(k,857) = -(rxt(k,507)*y(k,232) + rxt(k,508)*y(k,124))
         mat(k,2144) = -rxt(k,507)*y(k,253)
         mat(k,2403) = -rxt(k,508)*y(k,253)
         mat(k,1150) = rxt(k,516)*y(k,250)
         mat(k,482) = .500_r8*rxt(k,510)*y(k,250)
         mat(k,1901) = rxt(k,516)*y(k,29) + .500_r8*rxt(k,510)*y(k,146)
         mat(k,1324) = -(rxt(k,511)*y(k,226) + rxt(k,512)*y(k,232) + rxt(k,513) &
                      *y(k,124))
         mat(k,1700) = -rxt(k,511)*y(k,254)
         mat(k,2166) = -rxt(k,512)*y(k,254)
         mat(k,2428) = -rxt(k,513)*y(k,254)
         mat(k,1016) = .060_r8*rxt(k,688)*y(k,135)
         mat(k,1075) = rxt(k,517)*y(k,250)
         mat(k,1048) = .060_r8*rxt(k,693)*y(k,135)
         mat(k,2237) = .060_r8*rxt(k,688)*y(k,6) + .060_r8*rxt(k,693)*y(k,110)
         mat(k,357) = rxt(k,515)*y(k,250)
         mat(k,1227) = .150_r8*rxt(k,718)*y(k,250)
         mat(k,1933) = rxt(k,517)*y(k,48) + rxt(k,515)*y(k,147) + .150_r8*rxt(k,718) &
                      *y(k,178)
         mat(k,1283) = -(rxt(k,704)*y(k,226) + rxt(k,705)*y(k,232) + rxt(k,706) &
                      *y(k,124))
         mat(k,1698) = -rxt(k,704)*y(k,255)
         mat(k,2164) = -rxt(k,705)*y(k,255)
         mat(k,2426) = -rxt(k,706)*y(k,255)
         mat(k,2333) = .500_r8*rxt(k,715)*y(k,177)
         mat(k,553) = rxt(k,708)*y(k,250)
         mat(k,1140) = .500_r8*rxt(k,715)*y(k,126) + rxt(k,716)*y(k,250)
         mat(k,1931) = rxt(k,708)*y(k,174) + rxt(k,716)*y(k,177)
         mat(k,1256) = -(rxt(k,710)*y(k,226) + rxt(k,711)*y(k,232) + rxt(k,712) &
                      *y(k,124))
         mat(k,1697) = -rxt(k,710)*y(k,256)
         mat(k,2163) = -rxt(k,711)*y(k,256)
         mat(k,2425) = -rxt(k,712)*y(k,256)
         mat(k,1014) = rxt(k,689)*y(k,250)
         mat(k,1046) = rxt(k,694)*y(k,250)
         mat(k,426) = rxt(k,714)*y(k,250)
         mat(k,1930) = rxt(k,689)*y(k,6) + rxt(k,694)*y(k,110) + rxt(k,714)*y(k,176)
         mat(k,787) = -(rxt(k,670)*y(k,232) + rxt(k,671)*y(k,124))
         mat(k,2138) = -rxt(k,670)*y(k,257)
         mat(k,2399) = -rxt(k,671)*y(k,257)
         mat(k,561) = rxt(k,673)*y(k,250)
         mat(k,117) = .650_r8*rxt(k,674)*y(k,250)
         mat(k,1895) = rxt(k,673)*y(k,179) + .650_r8*rxt(k,674)*y(k,180)
         mat(k,83) = -(rxt(k,283)*y(k,133) + rxt(k,284)*y(k,134))
         mat(k,2459) = -rxt(k,283)*y(k,258)
         mat(k,1727) = -rxt(k,284)*y(k,258)
         mat(k,242) = -(rxt(k,203)*y(k,77) + rxt(k,204)*y(k,263) + (rxt(k,206) &
                      + rxt(k,207)) * y(k,134) + rxt(k,208)*y(k,135) + (rxt(k,356) &
                      + rxt(k,357)) * y(k,85) + (rxt(k,390) + rxt(k,391)) * y(k,81) &
                      + rxt(k,399)*y(k,64) + rxt(k,400)*y(k,65) + rxt(k,452)*y(k,86))
         mat(k,1335) = -rxt(k,203)*y(k,259)
         mat(k,2588) = -rxt(k,204)*y(k,259)
         mat(k,1729) = -(rxt(k,206) + rxt(k,207)) * y(k,259)
         mat(k,2205) = -rxt(k,208)*y(k,259)
         mat(k,1631) = -(rxt(k,356) + rxt(k,357)) * y(k,259)
         mat(k,846) = -(rxt(k,390) + rxt(k,391)) * y(k,259)
         mat(k,144) = -rxt(k,399)*y(k,259)
         mat(k,208) = -rxt(k,400)*y(k,259)
         mat(k,200) = -rxt(k,452)*y(k,259)
         mat(k,1361) = -(rxt(k,618)*y(k,225) + rxt(k,619)*y(k,226) + rxt(k,620) &
                      *y(k,232) + rxt(k,621)*y(k,124) + rxt(k,622)*y(k,126))
         mat(k,1585) = -rxt(k,618)*y(k,260)
         mat(k,1701) = -rxt(k,619)*y(k,260)
         mat(k,2168) = -rxt(k,620)*y(k,260)
         mat(k,2429) = -rxt(k,621)*y(k,260)
         mat(k,2336) = -rxt(k,622)*y(k,260)
         mat(k,174) = rxt(k,571)*y(k,250)
         mat(k,311) = rxt(k,572)*y(k,250)
         mat(k,70) = rxt(k,573)*y(k,250)
         mat(k,676) = .400_r8*rxt(k,602)*y(k,250)
         mat(k,148) = .500_r8*rxt(k,625)*y(k,250)
         mat(k,1935) = rxt(k,571)*y(k,94) + rxt(k,572)*y(k,96) + rxt(k,573)*y(k,97) &
                      + .400_r8*rxt(k,602)*y(k,103) + .500_r8*rxt(k,625)*y(k,207)
         mat(k,814) = -(rxt(k,678)*y(k,232) + rxt(k,679)*y(k,124))
         mat(k,2140) = -rxt(k,678)*y(k,261)
         mat(k,2401) = -rxt(k,679)*y(k,261)
         mat(k,162) = .560_r8*rxt(k,677)*y(k,250)
         mat(k,688) = rxt(k,681)*y(k,250)
         mat(k,1897) = .560_r8*rxt(k,677)*y(k,211) + rxt(k,681)*y(k,212)
         mat(k,615) = -(rxt(k,682)*y(k,232) + rxt(k,683)*y(k,124))
         mat(k,2129) = -rxt(k,682)*y(k,262)
         mat(k,2388) = -rxt(k,683)*y(k,262)
         mat(k,169) = .300_r8*rxt(k,685)*y(k,250)
         mat(k,393) = rxt(k,686)*y(k,250)
         mat(k,1882) = .300_r8*rxt(k,685)*y(k,213) + rxt(k,686)*y(k,214)
         mat(k,2615) = -(rxt(k,183)*y(k,246) + rxt(k,394)*y(k,73) + rxt(k,741) &
                      *y(k,152))
         mat(k,2297) = -rxt(k,183)*y(k,263)
         mat(k,879) = -rxt(k,394)*y(k,263)
         mat(k,216) = -rxt(k,741)*y(k,263)
         mat(k,264) = rxt(k,466)*y(k,250)
         mat(k,379) = rxt(k,500)*y(k,250)
         mat(k,60) = rxt(k,501)*y(k,250)
         mat(k,1792) = rxt(k,426)*y(k,250)
         mat(k,1206) = rxt(k,468)*y(k,250)
         mat(k,1079) = rxt(k,517)*y(k,250)
         mat(k,1424) = rxt(k,503)*y(k,250)
         mat(k,472) = rxt(k,476)*y(k,250)
         mat(k,458) = rxt(k,477)*y(k,250)
         mat(k,391) = rxt(k,435)*y(k,250)
         mat(k,2554) = rxt(k,216)*y(k,232)
         mat(k,1351) = rxt(k,224)*y(k,250)
         mat(k,579) = rxt(k,225)*y(k,250)
         mat(k,855) = rxt(k,377)*y(k,250)
         mat(k,1652) = (rxt(k,796)+rxt(k,801))*y(k,91) + (rxt(k,789)+rxt(k,795) &
                       +rxt(k,800))*y(k,92) + rxt(k,329)*y(k,250)
         mat(k,778) = rxt(k,439)*y(k,250)
         mat(k,2532) = rxt(k,274)*y(k,250)
         mat(k,422) = rxt(k,234)*y(k,250)
         mat(k,737) = (rxt(k,796)+rxt(k,801))*y(k,85)
         mat(k,845) = (rxt(k,789)+rxt(k,795)+rxt(k,800))*y(k,85) + rxt(k,333)*y(k,250)
         mat(k,1412) = .500_r8*rxt(k,537)*y(k,250)
         mat(k,51) = rxt(k,745)*y(k,250)
         mat(k,488) = rxt(k,510)*y(k,250)
         mat(k,361) = rxt(k,515)*y(k,250)
         mat(k,2196) = rxt(k,216)*y(k,76) + rxt(k,226)*y(k,250)
         mat(k,1964) = rxt(k,466)*y(k,28) + rxt(k,500)*y(k,30) + rxt(k,501)*y(k,31) &
                      + rxt(k,426)*y(k,42) + rxt(k,468)*y(k,45) + rxt(k,517)*y(k,48) &
                      + rxt(k,503)*y(k,49) + rxt(k,476)*y(k,50) + rxt(k,477)*y(k,51) &
                      + rxt(k,435)*y(k,53) + rxt(k,224)*y(k,77) + rxt(k,225)*y(k,79) &
                      + rxt(k,377)*y(k,81) + rxt(k,329)*y(k,85) + rxt(k,439)*y(k,87) &
                      + rxt(k,274)*y(k,89) + rxt(k,234)*y(k,90) + rxt(k,333)*y(k,92) &
                      + .500_r8*rxt(k,537)*y(k,105) + rxt(k,745)*y(k,120) + rxt(k,510) &
                      *y(k,146) + rxt(k,515)*y(k,147) + rxt(k,226)*y(k,232) &
                      + 2.000_r8*rxt(k,229)*y(k,250)
      end do
      end subroutine nlnmat11
      subroutine nlnmat_finit( avec_len, mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: dti(veclen)
      real(r8), intent(in) :: lmat(veclen,nzcnt)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k, 1) = lmat(k, 1)
         mat(k, 2) = lmat(k, 2)
         mat(k, 3) = lmat(k, 3)
         mat(k, 4) = lmat(k, 4)
         mat(k, 5) = lmat(k, 5)
         mat(k, 6) = lmat(k, 6)
         mat(k, 7) = lmat(k, 7)
         mat(k, 8) = lmat(k, 8)
         mat(k, 9) = lmat(k, 9)
         mat(k, 10) = lmat(k, 10)
         mat(k, 11) = lmat(k, 11)
         mat(k, 12) = lmat(k, 12)
         mat(k, 13) = lmat(k, 13)
         mat(k, 14) = lmat(k, 14)
         mat(k, 15) = lmat(k, 15)
         mat(k, 16) = lmat(k, 16)
         mat(k, 17) = lmat(k, 17)
         mat(k, 18) = lmat(k, 18)
         mat(k, 19) = lmat(k, 19)
         mat(k, 20) = lmat(k, 20)
         mat(k, 21) = lmat(k, 21)
         mat(k, 22) = lmat(k, 22)
         mat(k, 23) = lmat(k, 23)
         mat(k, 24) = lmat(k, 24)
         mat(k, 25) = lmat(k, 25)
         mat(k, 26) = lmat(k, 26)
         mat(k, 27) = lmat(k, 27)
         mat(k, 28) = lmat(k, 28)
         mat(k, 29) = lmat(k, 29)
         mat(k, 30) = lmat(k, 30)
         mat(k, 31) = lmat(k, 31)
         mat(k, 32) = lmat(k, 32)
         mat(k, 33) = lmat(k, 33)
         mat(k, 39) = mat(k, 39) + lmat(k, 39)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 46) = lmat(k, 46)
         mat(k, 47) = lmat(k, 47)
         mat(k, 48) = lmat(k, 48)
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 52) = mat(k, 52) + lmat(k, 52)
         mat(k, 55) = lmat(k, 55)
         mat(k, 56) = lmat(k, 56)
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 63) = lmat(k, 63)
         mat(k, 64) = lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = lmat(k, 68)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 72) = lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 82) = lmat(k, 82)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 91) = mat(k, 91) + lmat(k, 91)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = lmat(k, 98)
         mat(k, 99) = lmat(k, 99)
         mat(k, 100) = lmat(k, 100)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 106) = lmat(k, 106)
         mat(k, 107) = lmat(k, 107)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 125) = lmat(k, 125)
         mat(k, 126) = lmat(k, 126)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = lmat(k, 128)
         mat(k, 129) = mat(k, 129) + lmat(k, 129)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 132) = mat(k, 132) + lmat(k, 132)
         mat(k, 134) = lmat(k, 134)
         mat(k, 135) = mat(k, 135) + lmat(k, 135)
         mat(k, 137) = lmat(k, 137)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 140) = mat(k, 140) + lmat(k, 140)
         mat(k, 142) = lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = lmat(k, 152)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 173) = lmat(k, 173)
         mat(k, 175) = mat(k, 175) + lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = lmat(k, 178)
         mat(k, 179) = lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = lmat(k, 184)
         mat(k, 185) = lmat(k, 185)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 189) = lmat(k, 189)
         mat(k, 190) = lmat(k, 190)
         mat(k, 191) = lmat(k, 191)
         mat(k, 192) = lmat(k, 192)
         mat(k, 193) = lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 196) = lmat(k, 196)
         mat(k, 197) = lmat(k, 197)
         mat(k, 198) = lmat(k, 198)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 204) = lmat(k, 204)
         mat(k, 205) = lmat(k, 205)
         mat(k, 206) = lmat(k, 206)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 213) = mat(k, 213) + lmat(k, 213)
         mat(k, 214) = lmat(k, 214)
         mat(k, 215) = lmat(k, 215)
         mat(k, 217) = mat(k, 217) + lmat(k, 217)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 222) = lmat(k, 222)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 226) = lmat(k, 226)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 228) = lmat(k, 228)
         mat(k, 229) = mat(k, 229) + lmat(k, 229)
         mat(k, 231) = lmat(k, 231)
         mat(k, 234) = mat(k, 234) + lmat(k, 234)
         mat(k, 236) = lmat(k, 236)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 240) = lmat(k, 240)
         mat(k, 241) = lmat(k, 241)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 243) = mat(k, 243) + lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 247) = lmat(k, 247)
         mat(k, 248) = lmat(k, 248)
         mat(k, 249) = lmat(k, 249)
         mat(k, 250) = lmat(k, 250)
         mat(k, 251) = lmat(k, 251)
         mat(k, 252) = lmat(k, 252)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 265) = lmat(k, 265)
         mat(k, 266) = lmat(k, 266)
         mat(k, 267) = lmat(k, 267)
         mat(k, 268) = lmat(k, 268)
         mat(k, 269) = lmat(k, 269)
         mat(k, 270) = lmat(k, 270)
         mat(k, 271) = lmat(k, 271)
         mat(k, 272) = lmat(k, 272)
         mat(k, 273) = lmat(k, 273)
         mat(k, 274) = lmat(k, 274)
         mat(k, 275) = lmat(k, 275)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 279) = mat(k, 279) + lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 285) = lmat(k, 285)
         mat(k, 286) = lmat(k, 286)
         mat(k, 287) = mat(k, 287) + lmat(k, 287)
         mat(k, 290) = mat(k, 290) + lmat(k, 290)
         mat(k, 291) = lmat(k, 291)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 295) = mat(k, 295) + lmat(k, 295)
         mat(k, 296) = lmat(k, 296)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 299) = lmat(k, 299)
         mat(k, 300) = lmat(k, 300)
         mat(k, 301) = mat(k, 301) + lmat(k, 301)
         mat(k, 307) = mat(k, 307) + lmat(k, 307)
         mat(k, 308) = lmat(k, 308)
         mat(k, 309) = lmat(k, 309)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 314) = lmat(k, 314)
         mat(k, 315) = mat(k, 315) + lmat(k, 315)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 318) = lmat(k, 318)
         mat(k, 319) = lmat(k, 319)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 322) = lmat(k, 322)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 336) = mat(k, 336) + lmat(k, 336)
         mat(k, 344) = mat(k, 344) + lmat(k, 344)
         mat(k, 349) = lmat(k, 349)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 351) = lmat(k, 351)
         mat(k, 353) = mat(k, 353) + lmat(k, 353)
         mat(k, 354) = lmat(k, 354)
         mat(k, 355) = lmat(k, 355)
         mat(k, 356) = mat(k, 356) + lmat(k, 356)
         mat(k, 358) = lmat(k, 358)
         mat(k, 359) = lmat(k, 359)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 362) = mat(k, 362) + lmat(k, 362)
         mat(k, 363) = lmat(k, 363)
         mat(k, 365) = lmat(k, 365)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 367) = lmat(k, 367)
         mat(k, 368) = mat(k, 368) + lmat(k, 368)
         mat(k, 369) = lmat(k, 369)
         mat(k, 372) = lmat(k, 372)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 374) = mat(k, 374) + lmat(k, 374)
         mat(k, 376) = lmat(k, 376)
         mat(k, 377) = mat(k, 377) + lmat(k, 377)
         mat(k, 378) = lmat(k, 378)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = lmat(k, 382)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = lmat(k, 384)
         mat(k, 385) = lmat(k, 385)
         mat(k, 386) = mat(k, 386) + lmat(k, 386)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 390) = lmat(k, 390)
         mat(k, 392) = mat(k, 392) + lmat(k, 392)
         mat(k, 394) = lmat(k, 394)
         mat(k, 395) = lmat(k, 395)
         mat(k, 396) = mat(k, 396) + lmat(k, 396)
         mat(k, 397) = lmat(k, 397)
         mat(k, 398) = lmat(k, 398)
         mat(k, 399) = lmat(k, 399)
         mat(k, 400) = lmat(k, 400)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 404) = lmat(k, 404)
         mat(k, 405) = mat(k, 405) + lmat(k, 405)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 420) = lmat(k, 420)
         mat(k, 421) = lmat(k, 421)
         mat(k, 423) = mat(k, 423) + lmat(k, 423)
         mat(k, 424) = lmat(k, 424)
         mat(k, 425) = lmat(k, 425)
         mat(k, 427) = lmat(k, 427)
         mat(k, 428) = mat(k, 428) + lmat(k, 428)
         mat(k, 429) = lmat(k, 429)
         mat(k, 430) = mat(k, 430) + lmat(k, 430)
         mat(k, 431) = mat(k, 431) + lmat(k, 431)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 433) = lmat(k, 433)
         mat(k, 436) = mat(k, 436) + lmat(k, 436)
         mat(k, 439) = mat(k, 439) + lmat(k, 439)
         mat(k, 440) = lmat(k, 440)
         mat(k, 441) = lmat(k, 441)
         mat(k, 444) = lmat(k, 444)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 455) = lmat(k, 455)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 459) = mat(k, 459) + lmat(k, 459)
         mat(k, 461) = lmat(k, 461)
         mat(k, 462) = lmat(k, 462)
         mat(k, 464) = mat(k, 464) + lmat(k, 464)
         mat(k, 469) = mat(k, 469) + lmat(k, 469)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 475) = lmat(k, 475)
         mat(k, 478) = lmat(k, 478)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 483) = lmat(k, 483)
         mat(k, 485) = lmat(k, 485)
         mat(k, 486) = mat(k, 486) + lmat(k, 486)
         mat(k, 487) = lmat(k, 487)
         mat(k, 489) = lmat(k, 489)
         mat(k, 490) = mat(k, 490) + lmat(k, 490)
         mat(k, 491) = lmat(k, 491)
         mat(k, 493) = mat(k, 493) + lmat(k, 493)
         mat(k, 496) = mat(k, 496) + lmat(k, 496)
         mat(k, 505) = mat(k, 505) + lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = lmat(k, 507)
         mat(k, 508) = lmat(k, 508)
         mat(k, 509) = mat(k, 509) + lmat(k, 509)
         mat(k, 511) = lmat(k, 511)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 515) = mat(k, 515) + lmat(k, 515)
         mat(k, 523) = mat(k, 523) + lmat(k, 523)
         mat(k, 527) = lmat(k, 527)
         mat(k, 529) = mat(k, 529) + lmat(k, 529)
         mat(k, 531) = mat(k, 531) + lmat(k, 531)
         mat(k, 538) = lmat(k, 538)
         mat(k, 540) = mat(k, 540) + lmat(k, 540)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 549) = lmat(k, 549)
         mat(k, 550) = lmat(k, 550)
         mat(k, 551) = lmat(k, 551)
         mat(k, 552) = lmat(k, 552)
         mat(k, 554) = lmat(k, 554)
         mat(k, 555) = mat(k, 555) + lmat(k, 555)
         mat(k, 556) = lmat(k, 556)
         mat(k, 557) = lmat(k, 557)
         mat(k, 558) = lmat(k, 558)
         mat(k, 559) = mat(k, 559) + lmat(k, 559)
         mat(k, 560) = lmat(k, 560)
         mat(k, 564) = lmat(k, 564)
         mat(k, 565) = lmat(k, 565)
         mat(k, 566) = mat(k, 566) + lmat(k, 566)
         mat(k, 567) = lmat(k, 567)
         mat(k, 568) = lmat(k, 568)
         mat(k, 569) = lmat(k, 569)
         mat(k, 570) = lmat(k, 570)
         mat(k, 571) = lmat(k, 571)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 575) = mat(k, 575) + lmat(k, 575)
         mat(k, 580) = mat(k, 580) + lmat(k, 580)
         mat(k, 588) = lmat(k, 588)
         mat(k, 591) = mat(k, 591) + lmat(k, 591)
         mat(k, 602) = mat(k, 602) + lmat(k, 602)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 608) = mat(k, 608) + lmat(k, 608)
         mat(k, 609) = lmat(k, 609)
         mat(k, 610) = lmat(k, 610)
         mat(k, 611) = lmat(k, 611)
         mat(k, 615) = mat(k, 615) + lmat(k, 615)
         mat(k, 624) = mat(k, 624) + lmat(k, 624)
         mat(k, 625) = mat(k, 625) + lmat(k, 625)
         mat(k, 627) = mat(k, 627) + lmat(k, 627)
         mat(k, 628) = lmat(k, 628)
         mat(k, 629) = lmat(k, 629)
         mat(k, 630) = lmat(k, 630)
         mat(k, 631) = mat(k, 631) + lmat(k, 631)
         mat(k, 633) = mat(k, 633) + lmat(k, 633)
         mat(k, 635) = mat(k, 635) + lmat(k, 635)
         mat(k, 637) = mat(k, 637) + lmat(k, 637)
         mat(k, 639) = mat(k, 639) + lmat(k, 639)
         mat(k, 640) = lmat(k, 640)
         mat(k, 643) = mat(k, 643) + lmat(k, 643)
         mat(k, 653) = lmat(k, 653)
         mat(k, 654) = mat(k, 654) + lmat(k, 654)
         mat(k, 657) = lmat(k, 657)
         mat(k, 659) = lmat(k, 659)
         mat(k, 661) = lmat(k, 661)
         mat(k, 662) = mat(k, 662) + lmat(k, 662)
         mat(k, 663) = lmat(k, 663)
         mat(k, 666) = mat(k, 666) + lmat(k, 666)
         mat(k, 675) = mat(k, 675) + lmat(k, 675)
         mat(k, 677) = lmat(k, 677)
         mat(k, 678) = lmat(k, 678)
         mat(k, 679) = lmat(k, 679)
         mat(k, 681) = lmat(k, 681)
         mat(k, 682) = lmat(k, 682)
         mat(k, 683) = lmat(k, 683)
         mat(k, 684) = lmat(k, 684)
         mat(k, 685) = lmat(k, 685)
         mat(k, 686) = mat(k, 686) + lmat(k, 686)
         mat(k, 690) = lmat(k, 690)
         mat(k, 693) = lmat(k, 693)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 696) = lmat(k, 696)
         mat(k, 698) = mat(k, 698) + lmat(k, 698)
         mat(k, 711) = mat(k, 711) + lmat(k, 711)
         mat(k, 719) = mat(k, 719) + lmat(k, 719)
         mat(k, 730) = mat(k, 730) + lmat(k, 730)
         mat(k, 732) = lmat(k, 732)
         mat(k, 733) = mat(k, 733) + lmat(k, 733)
         mat(k, 740) = mat(k, 740) + lmat(k, 740)
         mat(k, 748) = mat(k, 748) + lmat(k, 748)
         mat(k, 749) = mat(k, 749) + lmat(k, 749)
         mat(k, 754) = mat(k, 754) + lmat(k, 754)
         mat(k, 757) = mat(k, 757) + lmat(k, 757)
         mat(k, 768) = mat(k, 768) + lmat(k, 768)
         mat(k, 775) = mat(k, 775) + lmat(k, 775)
         mat(k, 779) = mat(k, 779) + lmat(k, 779)
         mat(k, 787) = mat(k, 787) + lmat(k, 787)
         mat(k, 800) = mat(k, 800) + lmat(k, 800)
         mat(k, 814) = mat(k, 814) + lmat(k, 814)
         mat(k, 827) = mat(k, 827) + lmat(k, 827)
         mat(k, 839) = mat(k, 839) + lmat(k, 839)
         mat(k, 841) = mat(k, 841) + lmat(k, 841)
         mat(k, 843) = mat(k, 843) + lmat(k, 843)
         mat(k, 848) = mat(k, 848) + lmat(k, 848)
         mat(k, 849) = mat(k, 849) + lmat(k, 849)
         mat(k, 853) = mat(k, 853) + lmat(k, 853)
         mat(k, 857) = mat(k, 857) + lmat(k, 857)
         mat(k, 869) = mat(k, 869) + lmat(k, 869)
         mat(k, 870) = mat(k, 870) + lmat(k, 870)
         mat(k, 874) = lmat(k, 874)
         mat(k, 882) = mat(k, 882) + lmat(k, 882)
         mat(k, 892) = mat(k, 892) + lmat(k, 892)
         mat(k, 894) = lmat(k, 894)
         mat(k, 895) = lmat(k, 895)
         mat(k, 897) = mat(k, 897) + lmat(k, 897)
         mat(k, 899) = mat(k, 899) + lmat(k, 899)
         mat(k, 912) = lmat(k, 912)
         mat(k, 913) = mat(k, 913) + lmat(k, 913)
         mat(k, 914) = mat(k, 914) + lmat(k, 914)
         mat(k, 916) = mat(k, 916) + lmat(k, 916)
         mat(k, 918) = mat(k, 918) + lmat(k, 918)
         mat(k, 919) = lmat(k, 919)
         mat(k, 922) = lmat(k, 922)
         mat(k, 924) = lmat(k, 924)
         mat(k, 929) = lmat(k, 929)
         mat(k, 930) = mat(k, 930) + lmat(k, 930)
         mat(k, 943) = mat(k, 943) + lmat(k, 943)
         mat(k, 944) = mat(k, 944) + lmat(k, 944)
         mat(k, 947) = mat(k, 947) + lmat(k, 947)
         mat(k, 948) = lmat(k, 948)
         mat(k, 949) = mat(k, 949) + lmat(k, 949)
         mat(k, 950) = mat(k, 950) + lmat(k, 950)
         mat(k, 952) = mat(k, 952) + lmat(k, 952)
         mat(k, 963) = mat(k, 963) + lmat(k, 963)
         mat(k, 983) = mat(k, 983) + lmat(k, 983)
         mat(k,1005) = mat(k,1005) + lmat(k,1005)
         mat(k,1037) = mat(k,1037) + lmat(k,1037)
         mat(k,1062) = mat(k,1062) + lmat(k,1062)
         mat(k,1074) = mat(k,1074) + lmat(k,1074)
         mat(k,1076) = lmat(k,1076)
         mat(k,1077) = lmat(k,1077)
         mat(k,1080) = mat(k,1080) + lmat(k,1080)
         mat(k,1082) = lmat(k,1082)
         mat(k,1084) = lmat(k,1084)
         mat(k,1090) = lmat(k,1090)
         mat(k,1092) = mat(k,1092) + lmat(k,1092)
         mat(k,1100) = lmat(k,1100)
         mat(k,1101) = mat(k,1101) + lmat(k,1101)
         mat(k,1102) = lmat(k,1102)
         mat(k,1103) = mat(k,1103) + lmat(k,1103)
         mat(k,1104) = mat(k,1104) + lmat(k,1104)
         mat(k,1113) = mat(k,1113) + lmat(k,1113)
         mat(k,1114) = mat(k,1114) + lmat(k,1114)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1117) = mat(k,1117) + lmat(k,1117)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1137) = mat(k,1137) + lmat(k,1137)
         mat(k,1138) = lmat(k,1138)
         mat(k,1139) = lmat(k,1139)
         mat(k,1143) = lmat(k,1143)
         mat(k,1155) = mat(k,1155) + lmat(k,1155)
         mat(k,1171) = lmat(k,1171)
         mat(k,1178) = mat(k,1178) + lmat(k,1178)
         mat(k,1185) = lmat(k,1185)
         mat(k,1186) = mat(k,1186) + lmat(k,1186)
         mat(k,1188) = lmat(k,1188)
         mat(k,1190) = mat(k,1190) + lmat(k,1190)
         mat(k,1191) = lmat(k,1191)
         mat(k,1192) = mat(k,1192) + lmat(k,1192)
         mat(k,1194) = mat(k,1194) + lmat(k,1194)
         mat(k,1198) = mat(k,1198) + lmat(k,1198)
         mat(k,1199) = lmat(k,1199)
         mat(k,1201) = lmat(k,1201)
         mat(k,1203) = lmat(k,1203)
         mat(k,1208) = lmat(k,1208)
         mat(k,1209) = lmat(k,1209)
         mat(k,1210) = lmat(k,1210)
         mat(k,1211) = mat(k,1211) + lmat(k,1211)
         mat(k,1212) = lmat(k,1212)
         mat(k,1213) = lmat(k,1213)
         mat(k,1215) = lmat(k,1215)
         mat(k,1217) = lmat(k,1217)
         mat(k,1219) = lmat(k,1219)
         mat(k,1220) = mat(k,1220) + lmat(k,1220)
         mat(k,1222) = lmat(k,1222)
         mat(k,1224) = mat(k,1224) + lmat(k,1224)
         mat(k,1225) = mat(k,1225) + lmat(k,1225)
         mat(k,1226) = mat(k,1226) + lmat(k,1226)
         mat(k,1227) = mat(k,1227) + lmat(k,1227)
         mat(k,1228) = mat(k,1228) + lmat(k,1228)
         mat(k,1230) = mat(k,1230) + lmat(k,1230)
         mat(k,1232) = mat(k,1232) + lmat(k,1232)
         mat(k,1234) = mat(k,1234) + lmat(k,1234)
         mat(k,1238) = mat(k,1238) + lmat(k,1238)
         mat(k,1240) = lmat(k,1240)
         mat(k,1241) = lmat(k,1241)
         mat(k,1243) = mat(k,1243) + lmat(k,1243)
         mat(k,1256) = mat(k,1256) + lmat(k,1256)
         mat(k,1283) = mat(k,1283) + lmat(k,1283)
         mat(k,1305) = mat(k,1305) + lmat(k,1305)
         mat(k,1324) = mat(k,1324) + lmat(k,1324)
         mat(k,1341) = mat(k,1341) + lmat(k,1341)
         mat(k,1361) = mat(k,1361) + lmat(k,1361)
         mat(k,1384) = mat(k,1384) + lmat(k,1384)
         mat(k,1401) = mat(k,1401) + lmat(k,1401)
         mat(k,1402) = mat(k,1402) + lmat(k,1402)
         mat(k,1405) = mat(k,1405) + lmat(k,1405)
         mat(k,1406) = mat(k,1406) + lmat(k,1406)
         mat(k,1407) = mat(k,1407) + lmat(k,1407)
         mat(k,1409) = mat(k,1409) + lmat(k,1409)
         mat(k,1416) = mat(k,1416) + lmat(k,1416)
         mat(k,1417) = mat(k,1417) + lmat(k,1417)
         mat(k,1418) = mat(k,1418) + lmat(k,1418)
         mat(k,1420) = lmat(k,1420)
         mat(k,1425) = lmat(k,1425)
         mat(k,1448) = mat(k,1448) + lmat(k,1448)
         mat(k,1457) = mat(k,1457) + lmat(k,1457)
         mat(k,1481) = mat(k,1481) + lmat(k,1481)
         mat(k,1510) = mat(k,1510) + lmat(k,1510)
         mat(k,1526) = lmat(k,1526)
         mat(k,1528) = mat(k,1528) + lmat(k,1528)
         mat(k,1532) = mat(k,1532) + lmat(k,1532)
         mat(k,1534) = mat(k,1534) + lmat(k,1534)
         mat(k,1535) = lmat(k,1535)
         mat(k,1556) = mat(k,1556) + lmat(k,1556)
         mat(k,1594) = mat(k,1594) + lmat(k,1594)
         mat(k,1610) = lmat(k,1610)
         mat(k,1615) = mat(k,1615) + lmat(k,1615)
         mat(k,1625) = mat(k,1625) + lmat(k,1625)
         mat(k,1640) = mat(k,1640) + lmat(k,1640)
         mat(k,1645) = mat(k,1645) + lmat(k,1645)
         mat(k,1650) = mat(k,1650) + lmat(k,1650)
         mat(k,1657) = mat(k,1657) + lmat(k,1657)
         mat(k,1711) = mat(k,1711) + lmat(k,1711)
         mat(k,1735) = mat(k,1735) + lmat(k,1735)
         mat(k,1738) = mat(k,1738) + lmat(k,1738)
         mat(k,1740) = lmat(k,1740)
         mat(k,1750) = mat(k,1750) + lmat(k,1750)
         mat(k,1758) = mat(k,1758) + lmat(k,1758)
         mat(k,1760) = mat(k,1760) + lmat(k,1760)
         mat(k,1772) = mat(k,1772) + lmat(k,1772)
         mat(k,1773) = lmat(k,1773)
         mat(k,1778) = mat(k,1778) + lmat(k,1778)
         mat(k,1790) = mat(k,1790) + lmat(k,1790)
         mat(k,1816) = lmat(k,1816)
         mat(k,1828) = lmat(k,1828)
         mat(k,1947) = mat(k,1947) + lmat(k,1947)
         mat(k,1948) = mat(k,1948) + lmat(k,1948)
         mat(k,1951) = mat(k,1951) + lmat(k,1951)
         mat(k,1954) = mat(k,1954) + lmat(k,1954)
         mat(k,1955) = mat(k,1955) + lmat(k,1955)
         mat(k,1964) = mat(k,1964) + lmat(k,1964)
         mat(k,1985) = mat(k,1985) + lmat(k,1985)
         mat(k,1987) = mat(k,1987) + lmat(k,1987)
         mat(k,1993) = mat(k,1993) + lmat(k,1993)
         mat(k,2034) = mat(k,2034) + lmat(k,2034)
         mat(k,2036) = mat(k,2036) + lmat(k,2036)
         mat(k,2042) = mat(k,2042) + lmat(k,2042)
         mat(k,2043) = mat(k,2043) + lmat(k,2043)
         mat(k,2044) = mat(k,2044) + lmat(k,2044)
         mat(k,2068) = mat(k,2068) + lmat(k,2068)
         mat(k,2071) = mat(k,2071) + lmat(k,2071)
         mat(k,2072) = lmat(k,2072)
         mat(k,2073) = lmat(k,2073)
         mat(k,2079) = mat(k,2079) + lmat(k,2079)
         mat(k,2080) = mat(k,2080) + lmat(k,2080)
         mat(k,2126) = mat(k,2126) + lmat(k,2126)
         mat(k,2187) = mat(k,2187) + lmat(k,2187)
         mat(k,2202) = mat(k,2202) + lmat(k,2202)
         mat(k,2250) = mat(k,2250) + lmat(k,2250)
         mat(k,2257) = mat(k,2257) + lmat(k,2257)
         mat(k,2258) = mat(k,2258) + lmat(k,2258)
         mat(k,2261) = mat(k,2261) + lmat(k,2261)
         mat(k,2267) = mat(k,2267) + lmat(k,2267)
         mat(k,2269) = mat(k,2269) + lmat(k,2269)
         mat(k,2274) = mat(k,2274) + lmat(k,2274)
         mat(k,2278) = mat(k,2278) + lmat(k,2278)
         mat(k,2280) = mat(k,2280) + lmat(k,2280)
         mat(k,2281) = lmat(k,2281)
         mat(k,2282) = mat(k,2282) + lmat(k,2282)
         mat(k,2283) = lmat(k,2283)
         mat(k,2284) = mat(k,2284) + lmat(k,2284)
         mat(k,2287) = mat(k,2287) + lmat(k,2287)
         mat(k,2288) = lmat(k,2288)
         mat(k,2290) = mat(k,2290) + lmat(k,2290)
         mat(k,2292) = lmat(k,2292)
         mat(k,2293) = mat(k,2293) + lmat(k,2293)
         mat(k,2295) = mat(k,2295) + lmat(k,2295)
         mat(k,2349) = mat(k,2349) + lmat(k,2349)
         mat(k,2353) = mat(k,2353) + lmat(k,2353)
         mat(k,2358) = mat(k,2358) + lmat(k,2358)
         mat(k,2359) = mat(k,2359) + lmat(k,2359)
         mat(k,2360) = mat(k,2360) + lmat(k,2360)
         mat(k,2361) = mat(k,2361) + lmat(k,2361)
         mat(k,2396) = mat(k,2396) + lmat(k,2396)
         mat(k,2397) = lmat(k,2397)
         mat(k,2408) = mat(k,2408) + lmat(k,2408)
         mat(k,2451) = mat(k,2451) + lmat(k,2451)
         mat(k,2452) = mat(k,2452) + lmat(k,2452)
         mat(k,2470) = mat(k,2470) + lmat(k,2470)
         mat(k,2476) = lmat(k,2476)
         mat(k,2503) = mat(k,2503) + lmat(k,2503)
         mat(k,2519) = mat(k,2519) + lmat(k,2519)
         mat(k,2521) = lmat(k,2521)
         mat(k,2529) = mat(k,2529) + lmat(k,2529)
         mat(k,2552) = mat(k,2552) + lmat(k,2552)
         mat(k,2568) = mat(k,2568) + lmat(k,2568)
         mat(k,2580) = mat(k,2580) + lmat(k,2580)
         mat(k,2583) = mat(k,2583) + lmat(k,2583)
         mat(k,2595) = lmat(k,2595)
         mat(k,2602) = mat(k,2602) + lmat(k,2602)
         mat(k,2608) = mat(k,2608) + lmat(k,2608)
         mat(k,2611) = lmat(k,2611)
         mat(k,2613) = lmat(k,2613)
         mat(k,2615) = mat(k,2615) + lmat(k,2615)
         mat(k, 163) = 0._r8
         mat(k, 164) = 0._r8
         mat(k, 237) = 0._r8
         mat(k, 289) = 0._r8
         mat(k, 331) = 0._r8
         mat(k, 334) = 0._r8
         mat(k, 335) = 0._r8
         mat(k, 437) = 0._r8
         mat(k, 495) = 0._r8
         mat(k, 516) = 0._r8
         mat(k, 519) = 0._r8
         mat(k, 544) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 563) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 596) = 0._r8
         mat(k, 601) = 0._r8
         mat(k, 620) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 636) = 0._r8
         mat(k, 648) = 0._r8
         mat(k, 650) = 0._r8
         mat(k, 651) = 0._r8
         mat(k, 655) = 0._r8
         mat(k, 658) = 0._r8
         mat(k, 660) = 0._r8
         mat(k, 687) = 0._r8
         mat(k, 689) = 0._r8
         mat(k, 691) = 0._r8
         mat(k, 692) = 0._r8
         mat(k, 694) = 0._r8
         mat(k, 700) = 0._r8
         mat(k, 710) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 743) = 0._r8
         mat(k, 746) = 0._r8
         mat(k, 760) = 0._r8
         mat(k, 761) = 0._r8
         mat(k, 762) = 0._r8
         mat(k, 770) = 0._r8
         mat(k, 771) = 0._r8
         mat(k, 773) = 0._r8
         mat(k, 788) = 0._r8
         mat(k, 789) = 0._r8
         mat(k, 791) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 813) = 0._r8
         mat(k, 815) = 0._r8
         mat(k, 817) = 0._r8
         mat(k, 820) = 0._r8
         mat(k, 822) = 0._r8
         mat(k, 823) = 0._r8
         mat(k, 861) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 868) = 0._r8
         mat(k, 885) = 0._r8
         mat(k, 886) = 0._r8
         mat(k, 888) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 911) = 0._r8
         mat(k, 923) = 0._r8
         mat(k, 926) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 961) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 981) = 0._r8
         mat(k,1002) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1006) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1010) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1024) = 0._r8
         mat(k,1034) = 0._r8
         mat(k,1036) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1042) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1056) = 0._r8
         mat(k,1061) = 0._r8
         mat(k,1066) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1070) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1122) = 0._r8
         mat(k,1129) = 0._r8
         mat(k,1134) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1158) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1174) = 0._r8
         mat(k,1175) = 0._r8
         mat(k,1176) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1197) = 0._r8
         mat(k,1214) = 0._r8
         mat(k,1216) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1223) = 0._r8
         mat(k,1229) = 0._r8
         mat(k,1233) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1254) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1258) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1260) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1266) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1269) = 0._r8
         mat(k,1270) = 0._r8
         mat(k,1276) = 0._r8
         mat(k,1278) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1285) = 0._r8
         mat(k,1288) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1307) = 0._r8
         mat(k,1310) = 0._r8
         mat(k,1315) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1317) = 0._r8
         mat(k,1321) = 0._r8
         mat(k,1333) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1349) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1371) = 0._r8
         mat(k,1374) = 0._r8
         mat(k,1376) = 0._r8
         mat(k,1381) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1385) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1391) = 0._r8
         mat(k,1396) = 0._r8
         mat(k,1397) = 0._r8
         mat(k,1399) = 0._r8
         mat(k,1411) = 0._r8
         mat(k,1415) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1430) = 0._r8
         mat(k,1433) = 0._r8
         mat(k,1435) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1441) = 0._r8
         mat(k,1442) = 0._r8
         mat(k,1444) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1446) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1451) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1461) = 0._r8
         mat(k,1462) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1465) = 0._r8
         mat(k,1469) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1474) = 0._r8
         mat(k,1475) = 0._r8
         mat(k,1476) = 0._r8
         mat(k,1477) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1480) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1484) = 0._r8
         mat(k,1488) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1494) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1500) = 0._r8
         mat(k,1508) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1515) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1522) = 0._r8
         mat(k,1524) = 0._r8
         mat(k,1529) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1536) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1542) = 0._r8
         mat(k,1543) = 0._r8
         mat(k,1544) = 0._r8
         mat(k,1545) = 0._r8
         mat(k,1546) = 0._r8
         mat(k,1554) = 0._r8
         mat(k,1555) = 0._r8
         mat(k,1566) = 0._r8
         mat(k,1567) = 0._r8
         mat(k,1568) = 0._r8
         mat(k,1577) = 0._r8
         mat(k,1579) = 0._r8
         mat(k,1596) = 0._r8
         mat(k,1602) = 0._r8
         mat(k,1604) = 0._r8
         mat(k,1605) = 0._r8
         mat(k,1606) = 0._r8
         mat(k,1609) = 0._r8
         mat(k,1613) = 0._r8
         mat(k,1622) = 0._r8
         mat(k,1628) = 0._r8
         mat(k,1638) = 0._r8
         mat(k,1639) = 0._r8
         mat(k,1641) = 0._r8
         mat(k,1644) = 0._r8
         mat(k,1647) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1654) = 0._r8
         mat(k,1660) = 0._r8
         mat(k,1663) = 0._r8
         mat(k,1664) = 0._r8
         mat(k,1665) = 0._r8
         mat(k,1667) = 0._r8
         mat(k,1669) = 0._r8
         mat(k,1681) = 0._r8
         mat(k,1685) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1719) = 0._r8
         mat(k,1720) = 0._r8
         mat(k,1722) = 0._r8
         mat(k,1723) = 0._r8
         mat(k,1724) = 0._r8
         mat(k,1744) = 0._r8
         mat(k,1745) = 0._r8
         mat(k,1747) = 0._r8
         mat(k,1749) = 0._r8
         mat(k,1751) = 0._r8
         mat(k,1752) = 0._r8
         mat(k,1753) = 0._r8
         mat(k,1754) = 0._r8
         mat(k,1755) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1763) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1769) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1776) = 0._r8
         mat(k,1777) = 0._r8
         mat(k,1780) = 0._r8
         mat(k,1781) = 0._r8
         mat(k,1784) = 0._r8
         mat(k,1785) = 0._r8
         mat(k,1787) = 0._r8
         mat(k,1791) = 0._r8
         mat(k,1835) = 0._r8
         mat(k,1857) = 0._r8
         mat(k,1877) = 0._r8
         mat(k,1896) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1902) = 0._r8
         mat(k,1911) = 0._r8
         mat(k,1936) = 0._r8
         mat(k,1957) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1976) = 0._r8
         mat(k,1989) = 0._r8
         mat(k,1990) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,1997) = 0._r8
         mat(k,2010) = 0._r8
         mat(k,2011) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2019) = 0._r8
         mat(k,2022) = 0._r8
         mat(k,2023) = 0._r8
         mat(k,2024) = 0._r8
         mat(k,2025) = 0._r8
         mat(k,2029) = 0._r8
         mat(k,2030) = 0._r8
         mat(k,2031) = 0._r8
         mat(k,2033) = 0._r8
         mat(k,2037) = 0._r8
         mat(k,2040) = 0._r8
         mat(k,2045) = 0._r8
         mat(k,2047) = 0._r8
         mat(k,2053) = 0._r8
         mat(k,2056) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2059) = 0._r8
         mat(k,2061) = 0._r8
         mat(k,2062) = 0._r8
         mat(k,2064) = 0._r8
         mat(k,2066) = 0._r8
         mat(k,2067) = 0._r8
         mat(k,2070) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2082) = 0._r8
         mat(k,2084) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2088) = 0._r8
         mat(k,2089) = 0._r8
         mat(k,2094) = 0._r8
         mat(k,2110) = 0._r8
         mat(k,2111) = 0._r8
         mat(k,2125) = 0._r8
         mat(k,2128) = 0._r8
         mat(k,2146) = 0._r8
         mat(k,2151) = 0._r8
         mat(k,2152) = 0._r8
         mat(k,2157) = 0._r8
         mat(k,2158) = 0._r8
         mat(k,2160) = 0._r8
         mat(k,2162) = 0._r8
         mat(k,2170) = 0._r8
         mat(k,2175) = 0._r8
         mat(k,2189) = 0._r8
         mat(k,2193) = 0._r8
         mat(k,2214) = 0._r8
         mat(k,2215) = 0._r8
         mat(k,2220) = 0._r8
         mat(k,2222) = 0._r8
         mat(k,2226) = 0._r8
         mat(k,2229) = 0._r8
         mat(k,2233) = 0._r8
         mat(k,2234) = 0._r8
         mat(k,2235) = 0._r8
         mat(k,2236) = 0._r8
         mat(k,2238) = 0._r8
         mat(k,2241) = 0._r8
         mat(k,2242) = 0._r8
         mat(k,2243) = 0._r8
         mat(k,2245) = 0._r8
         mat(k,2262) = 0._r8
         mat(k,2265) = 0._r8
         mat(k,2270) = 0._r8
         mat(k,2272) = 0._r8
         mat(k,2275) = 0._r8
         mat(k,2276) = 0._r8
         mat(k,2277) = 0._r8
         mat(k,2286) = 0._r8
         mat(k,2291) = 0._r8
         mat(k,2294) = 0._r8
         mat(k,2309) = 0._r8
         mat(k,2316) = 0._r8
         mat(k,2318) = 0._r8
         mat(k,2320) = 0._r8
         mat(k,2322) = 0._r8
         mat(k,2328) = 0._r8
         mat(k,2329) = 0._r8
         mat(k,2332) = 0._r8
         mat(k,2335) = 0._r8
         mat(k,2346) = 0._r8
         mat(k,2347) = 0._r8
         mat(k,2348) = 0._r8
         mat(k,2352) = 0._r8
         mat(k,2354) = 0._r8
         mat(k,2356) = 0._r8
         mat(k,2357) = 0._r8
         mat(k,2362) = 0._r8
         mat(k,2363) = 0._r8
         mat(k,2364) = 0._r8
         mat(k,2405) = 0._r8
         mat(k,2412) = 0._r8
         mat(k,2413) = 0._r8
         mat(k,2449) = 0._r8
         mat(k,2453) = 0._r8
         mat(k,2454) = 0._r8
         mat(k,2456) = 0._r8
         mat(k,2472) = 0._r8
         mat(k,2474) = 0._r8
         mat(k,2479) = 0._r8
         mat(k,2482) = 0._r8
         mat(k,2491) = 0._r8
         mat(k,2500) = 0._r8
         mat(k,2504) = 0._r8
         mat(k,2507) = 0._r8
         mat(k,2510) = 0._r8
         mat(k,2511) = 0._r8
         mat(k,2512) = 0._r8
         mat(k,2513) = 0._r8
         mat(k,2514) = 0._r8
         mat(k,2515) = 0._r8
         mat(k,2516) = 0._r8
         mat(k,2517) = 0._r8
         mat(k,2518) = 0._r8
         mat(k,2520) = 0._r8
         mat(k,2522) = 0._r8
         mat(k,2523) = 0._r8
         mat(k,2524) = 0._r8
         mat(k,2525) = 0._r8
         mat(k,2527) = 0._r8
         mat(k,2528) = 0._r8
         mat(k,2530) = 0._r8
         mat(k,2531) = 0._r8
         mat(k,2534) = 0._r8
         mat(k,2536) = 0._r8
         mat(k,2537) = 0._r8
         mat(k,2538) = 0._r8
         mat(k,2540) = 0._r8
         mat(k,2542) = 0._r8
         mat(k,2543) = 0._r8
         mat(k,2544) = 0._r8
         mat(k,2547) = 0._r8
         mat(k,2548) = 0._r8
         mat(k,2549) = 0._r8
         mat(k,2551) = 0._r8
         mat(k,2553) = 0._r8
         mat(k,2562) = 0._r8
         mat(k,2564) = 0._r8
         mat(k,2567) = 0._r8
         mat(k,2570) = 0._r8
         mat(k,2576) = 0._r8
         mat(k,2577) = 0._r8
         mat(k,2578) = 0._r8
         mat(k,2581) = 0._r8
         mat(k,2582) = 0._r8
         mat(k,2584) = 0._r8
         mat(k,2589) = 0._r8
         mat(k,2591) = 0._r8
         mat(k,2592) = 0._r8
         mat(k,2593) = 0._r8
         mat(k,2594) = 0._r8
         mat(k,2596) = 0._r8
         mat(k,2597) = 0._r8
         mat(k,2598) = 0._r8
         mat(k,2599) = 0._r8
         mat(k,2600) = 0._r8
         mat(k,2601) = 0._r8
         mat(k,2603) = 0._r8
         mat(k,2604) = 0._r8
         mat(k,2605) = 0._r8
         mat(k,2606) = 0._r8
         mat(k,2607) = 0._r8
         mat(k,2609) = 0._r8
         mat(k,2610) = 0._r8
         mat(k,2612) = 0._r8
         mat(k,2614) = 0._r8
         mat(k, 1) = mat(k, 1) - dti(k)
         mat(k, 2) = mat(k, 2) - dti(k)
         mat(k, 3) = mat(k, 3) - dti(k)
         mat(k, 4) = mat(k, 4) - dti(k)
         mat(k, 5) = mat(k, 5) - dti(k)
         mat(k, 6) = mat(k, 6) - dti(k)
         mat(k, 7) = mat(k, 7) - dti(k)
         mat(k, 8) = mat(k, 8) - dti(k)
         mat(k, 9) = mat(k, 9) - dti(k)
         mat(k, 10) = mat(k, 10) - dti(k)
         mat(k, 11) = mat(k, 11) - dti(k)
         mat(k, 12) = mat(k, 12) - dti(k)
         mat(k, 13) = mat(k, 13) - dti(k)
         mat(k, 14) = mat(k, 14) - dti(k)
         mat(k, 15) = mat(k, 15) - dti(k)
         mat(k, 16) = mat(k, 16) - dti(k)
         mat(k, 17) = mat(k, 17) - dti(k)
         mat(k, 18) = mat(k, 18) - dti(k)
         mat(k, 19) = mat(k, 19) - dti(k)
         mat(k, 20) = mat(k, 20) - dti(k)
         mat(k, 21) = mat(k, 21) - dti(k)
         mat(k, 22) = mat(k, 22) - dti(k)
         mat(k, 23) = mat(k, 23) - dti(k)
         mat(k, 24) = mat(k, 24) - dti(k)
         mat(k, 25) = mat(k, 25) - dti(k)
         mat(k, 26) = mat(k, 26) - dti(k)
         mat(k, 27) = mat(k, 27) - dti(k)
         mat(k, 28) = mat(k, 28) - dti(k)
         mat(k, 29) = mat(k, 29) - dti(k)
         mat(k, 30) = mat(k, 30) - dti(k)
         mat(k, 31) = mat(k, 31) - dti(k)
         mat(k, 32) = mat(k, 32) - dti(k)
         mat(k, 33) = mat(k, 33) - dti(k)
         mat(k, 39) = mat(k, 39) - dti(k)
         mat(k, 45) = mat(k, 45) - dti(k)
         mat(k, 46) = mat(k, 46) - dti(k)
         mat(k, 49) = mat(k, 49) - dti(k)
         mat(k, 52) = mat(k, 52) - dti(k)
         mat(k, 55) = mat(k, 55) - dti(k)
         mat(k, 57) = mat(k, 57) - dti(k)
         mat(k, 61) = mat(k, 61) - dti(k)
         mat(k, 63) = mat(k, 63) - dti(k)
         mat(k, 66) = mat(k, 66) - dti(k)
         mat(k, 69) = mat(k, 69) - dti(k)
         mat(k, 72) = mat(k, 72) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 78) = mat(k, 78) - dti(k)
         mat(k, 81) = mat(k, 81) - dti(k)
         mat(k, 83) = mat(k, 83) - dti(k)
         mat(k, 91) = mat(k, 91) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 101) = mat(k, 101) - dti(k)
         mat(k, 106) = mat(k, 106) - dti(k)
         mat(k, 113) = mat(k, 113) - dti(k)
         mat(k, 120) = mat(k, 120) - dti(k)
         mat(k, 125) = mat(k, 125) - dti(k)
         mat(k, 129) = mat(k, 129) - dti(k)
         mat(k, 132) = mat(k, 132) - dti(k)
         mat(k, 135) = mat(k, 135) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 140) = mat(k, 140) - dti(k)
         mat(k, 143) = mat(k, 143) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 150) = mat(k, 150) - dti(k)
         mat(k, 159) = mat(k, 159) - dti(k)
         mat(k, 167) = mat(k, 167) - dti(k)
         mat(k, 172) = mat(k, 172) - dti(k)
         mat(k, 177) = mat(k, 177) - dti(k)
         mat(k, 180) = mat(k, 180) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 186) = mat(k, 186) - dti(k)
         mat(k, 189) = mat(k, 189) - dti(k)
         mat(k, 194) = mat(k, 194) - dti(k)
         mat(k, 196) = mat(k, 196) - dti(k)
         mat(k, 199) = mat(k, 199) - dti(k)
         mat(k, 204) = mat(k, 204) - dti(k)
         mat(k, 207) = mat(k, 207) - dti(k)
         mat(k, 213) = mat(k, 213) - dti(k)
         mat(k, 217) = mat(k, 217) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 225) = mat(k, 225) - dti(k)
         mat(k, 229) = mat(k, 229) - dti(k)
         mat(k, 234) = mat(k, 234) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 242) = mat(k, 242) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 250) = mat(k, 250) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 259) = mat(k, 259) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 271) = mat(k, 271) - dti(k)
         mat(k, 276) = mat(k, 276) - dti(k)
         mat(k, 281) = mat(k, 281) - dti(k)
         mat(k, 287) = mat(k, 287) - dti(k)
         mat(k, 292) = mat(k, 292) - dti(k)
         mat(k, 298) = mat(k, 298) - dti(k)
         mat(k, 301) = mat(k, 301) - dti(k)
         mat(k, 307) = mat(k, 307) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 313) = mat(k, 313) - dti(k)
         mat(k, 316) = mat(k, 316) - dti(k)
         mat(k, 321) = mat(k, 321) - dti(k)
         mat(k, 333) = mat(k, 333) - dti(k)
         mat(k, 336) = mat(k, 336) - dti(k)
         mat(k, 344) = mat(k, 344) - dti(k)
         mat(k, 350) = mat(k, 350) - dti(k)
         mat(k, 356) = mat(k, 356) - dti(k)
         mat(k, 362) = mat(k, 362) - dti(k)
         mat(k, 368) = mat(k, 368) - dti(k)
         mat(k, 374) = mat(k, 374) - dti(k)
         mat(k, 380) = mat(k, 380) - dti(k)
         mat(k, 386) = mat(k, 386) - dti(k)
         mat(k, 392) = mat(k, 392) - dti(k)
         mat(k, 398) = mat(k, 398) - dti(k)
         mat(k, 401) = mat(k, 401) - dti(k)
         mat(k, 406) = mat(k, 406) - dti(k)
         mat(k, 416) = mat(k, 416) - dti(k)
         mat(k, 423) = mat(k, 423) - dti(k)
         mat(k, 430) = mat(k, 430) - dti(k)
         mat(k, 436) = mat(k, 436) - dti(k)
         mat(k, 439) = mat(k, 439) - dti(k)
         mat(k, 446) = mat(k, 446) - dti(k)
         mat(k, 453) = mat(k, 453) - dti(k)
         mat(k, 459) = mat(k, 459) - dti(k)
         mat(k, 464) = mat(k, 464) - dti(k)
         mat(k, 469) = mat(k, 469) - dti(k)
         mat(k, 473) = mat(k, 473) - dti(k)
         mat(k, 481) = mat(k, 481) - dti(k)
         mat(k, 490) = mat(k, 490) - dti(k)
         mat(k, 493) = mat(k, 493) - dti(k)
         mat(k, 496) = mat(k, 496) - dti(k)
         mat(k, 505) = mat(k, 505) - dti(k)
         mat(k, 515) = mat(k, 515) - dti(k)
         mat(k, 523) = mat(k, 523) - dti(k)
         mat(k, 531) = mat(k, 531) - dti(k)
         mat(k, 540) = mat(k, 540) - dti(k)
         mat(k, 548) = mat(k, 548) - dti(k)
         mat(k, 559) = mat(k, 559) - dti(k)
         mat(k, 568) = mat(k, 568) - dti(k)
         mat(k, 573) = mat(k, 573) - dti(k)
         mat(k, 580) = mat(k, 580) - dti(k)
         mat(k, 591) = mat(k, 591) - dti(k)
         mat(k, 602) = mat(k, 602) - dti(k)
         mat(k, 607) = mat(k, 607) - dti(k)
         mat(k, 615) = mat(k, 615) - dti(k)
         mat(k, 624) = mat(k, 624) - dti(k)
         mat(k, 631) = mat(k, 631) - dti(k)
         mat(k, 643) = mat(k, 643) - dti(k)
         mat(k, 654) = mat(k, 654) - dti(k)
         mat(k, 666) = mat(k, 666) - dti(k)
         mat(k, 675) = mat(k, 675) - dti(k)
         mat(k, 686) = mat(k, 686) - dti(k)
         mat(k, 698) = mat(k, 698) - dti(k)
         mat(k, 711) = mat(k, 711) - dti(k)
         mat(k, 719) = mat(k, 719) - dti(k)
         mat(k, 730) = mat(k, 730) - dti(k)
         mat(k, 740) = mat(k, 740) - dti(k)
         mat(k, 748) = mat(k, 748) - dti(k)
         mat(k, 757) = mat(k, 757) - dti(k)
         mat(k, 768) = mat(k, 768) - dti(k)
         mat(k, 775) = mat(k, 775) - dti(k)
         mat(k, 779) = mat(k, 779) - dti(k)
         mat(k, 787) = mat(k, 787) - dti(k)
         mat(k, 800) = mat(k, 800) - dti(k)
         mat(k, 814) = mat(k, 814) - dti(k)
         mat(k, 827) = mat(k, 827) - dti(k)
         mat(k, 839) = mat(k, 839) - dti(k)
         mat(k, 848) = mat(k, 848) - dti(k)
         mat(k, 857) = mat(k, 857) - dti(k)
         mat(k, 870) = mat(k, 870) - dti(k)
         mat(k, 882) = mat(k, 882) - dti(k)
         mat(k, 892) = mat(k, 892) - dti(k)
         mat(k, 899) = mat(k, 899) - dti(k)
         mat(k, 913) = mat(k, 913) - dti(k)
         mat(k, 918) = mat(k, 918) - dti(k)
         mat(k, 930) = mat(k, 930) - dti(k)
         mat(k, 944) = mat(k, 944) - dti(k)
         mat(k, 963) = mat(k, 963) - dti(k)
         mat(k, 983) = mat(k, 983) - dti(k)
         mat(k,1005) = mat(k,1005) - dti(k)
         mat(k,1037) = mat(k,1037) - dti(k)
         mat(k,1062) = mat(k,1062) - dti(k)
         mat(k,1074) = mat(k,1074) - dti(k)
         mat(k,1080) = mat(k,1080) - dti(k)
         mat(k,1092) = mat(k,1092) - dti(k)
         mat(k,1103) = mat(k,1103) - dti(k)
         mat(k,1117) = mat(k,1117) - dti(k)
         mat(k,1124) = mat(k,1124) - dti(k)
         mat(k,1137) = mat(k,1137) - dti(k)
         mat(k,1155) = mat(k,1155) - dti(k)
         mat(k,1178) = mat(k,1178) - dti(k)
         mat(k,1190) = mat(k,1190) - dti(k)
         mat(k,1198) = mat(k,1198) - dti(k)
         mat(k,1211) = mat(k,1211) - dti(k)
         mat(k,1225) = mat(k,1225) - dti(k)
         mat(k,1234) = mat(k,1234) - dti(k)
         mat(k,1238) = mat(k,1238) - dti(k)
         mat(k,1256) = mat(k,1256) - dti(k)
         mat(k,1283) = mat(k,1283) - dti(k)
         mat(k,1305) = mat(k,1305) - dti(k)
         mat(k,1324) = mat(k,1324) - dti(k)
         mat(k,1341) = mat(k,1341) - dti(k)
         mat(k,1361) = mat(k,1361) - dti(k)
         mat(k,1384) = mat(k,1384) - dti(k)
         mat(k,1402) = mat(k,1402) - dti(k)
         mat(k,1417) = mat(k,1417) - dti(k)
         mat(k,1448) = mat(k,1448) - dti(k)
         mat(k,1481) = mat(k,1481) - dti(k)
         mat(k,1510) = mat(k,1510) - dti(k)
         mat(k,1532) = mat(k,1532) - dti(k)
         mat(k,1556) = mat(k,1556) - dti(k)
         mat(k,1594) = mat(k,1594) - dti(k)
         mat(k,1615) = mat(k,1615) - dti(k)
         mat(k,1640) = mat(k,1640) - dti(k)
         mat(k,1657) = mat(k,1657) - dti(k)
         mat(k,1711) = mat(k,1711) - dti(k)
         mat(k,1750) = mat(k,1750) - dti(k)
         mat(k,1778) = mat(k,1778) - dti(k)
         mat(k,1951) = mat(k,1951) - dti(k)
         mat(k,1985) = mat(k,1985) - dti(k)
         mat(k,2036) = mat(k,2036) - dti(k)
         mat(k,2079) = mat(k,2079) - dti(k)
         mat(k,2187) = mat(k,2187) - dti(k)
         mat(k,2257) = mat(k,2257) - dti(k)
         mat(k,2290) = mat(k,2290) - dti(k)
         mat(k,2358) = mat(k,2358) - dti(k)
         mat(k,2451) = mat(k,2451) - dti(k)
         mat(k,2503) = mat(k,2503) - dti(k)
         mat(k,2529) = mat(k,2529) - dti(k)
         mat(k,2552) = mat(k,2552) - dti(k)
         mat(k,2583) = mat(k,2583) - dti(k)
         mat(k,2615) = mat(k,2615) - dti(k)
      end do
      end subroutine nlnmat_finit
      subroutine nlnmat( avec_len, mat, y, rxt, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: dti(veclen)
      real(r8), intent(in) :: lmat(veclen,nzcnt)
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
      call nlnmat01( avec_len, mat, y, rxt )
      call nlnmat02( avec_len, mat, y, rxt )
      call nlnmat03( avec_len, mat, y, rxt )
      call nlnmat04( avec_len, mat, y, rxt )
      call nlnmat05( avec_len, mat, y, rxt )
      call nlnmat06( avec_len, mat, y, rxt )
      call nlnmat07( avec_len, mat, y, rxt )
      call nlnmat08( avec_len, mat, y, rxt )
      call nlnmat09( avec_len, mat, y, rxt )
      call nlnmat10( avec_len, mat, y, rxt )
      call nlnmat11( avec_len, mat, y, rxt )
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
