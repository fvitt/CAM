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
         mat(k,985) = -(rxt(k,204)*y(k,17) + rxt(k,205)*y(k,87) + rxt(k,206)*y(k,70))
         mat(k,1599) = -rxt(k,204)*y(k,3)
         mat(k,1374) = -rxt(k,205)*y(k,3)
         mat(k,1762) = -rxt(k,206)*y(k,3)
         mat(k,1686) = 4.000_r8*rxt(k,207)*y(k,5) + (rxt(k,208)+rxt(k,209))*y(k,28) &
                      + rxt(k,212)*y(k,61) + rxt(k,215)*y(k,68) + rxt(k,266)*y(k,77) &
                      + rxt(k,216)*y(k,96)
         mat(k,66) = rxt(k,194)*y(k,92)
         mat(k,72) = rxt(k,220)*y(k,92)
         mat(k,293) = 2.000_r8*rxt(k,231)*y(k,25) + 2.000_r8*rxt(k,243)*y(k,92) &
                      + 2.000_r8*rxt(k,232)*y(k,96)
         mat(k,352) = rxt(k,233)*y(k,25) + rxt(k,244)*y(k,92) + rxt(k,234)*y(k,96)
         mat(k,248) = 3.000_r8*rxt(k,238)*y(k,25) + 3.000_r8*rxt(k,221)*y(k,92) &
                      + 3.000_r8*rxt(k,239)*y(k,96)
         mat(k,1488) = 2.000_r8*rxt(k,231)*y(k,16) + rxt(k,233)*y(k,18) &
                      + 3.000_r8*rxt(k,238)*y(k,24)
         mat(k,1793) = (rxt(k,208)+rxt(k,209))*y(k,5)
         mat(k,43) = 2.000_r8*rxt(k,222)*y(k,92)
         mat(k,541) = rxt(k,217)*y(k,68) + rxt(k,223)*y(k,92) + rxt(k,218)*y(k,96)
         mat(k,1707) = rxt(k,212)*y(k,5)
         mat(k,1574) = rxt(k,215)*y(k,5) + rxt(k,217)*y(k,45)
         mat(k,881) = rxt(k,266)*y(k,5)
         mat(k,1454) = rxt(k,194)*y(k,9) + rxt(k,220)*y(k,10) + 2.000_r8*rxt(k,243) &
                      *y(k,16) + rxt(k,244)*y(k,18) + 3.000_r8*rxt(k,221)*y(k,24) &
                      + 2.000_r8*rxt(k,222)*y(k,42) + rxt(k,223)*y(k,45)
         mat(k,1344) = rxt(k,216)*y(k,5) + 2.000_r8*rxt(k,232)*y(k,16) + rxt(k,234) &
                      *y(k,18) + 3.000_r8*rxt(k,239)*y(k,24) + rxt(k,218)*y(k,45)
         mat(k,1679) = rxt(k,210)*y(k,28)
         mat(k,1785) = rxt(k,210)*y(k,5)
         mat(k,947) = (rxt(k,288)+rxt(k,293))*y(k,53)
         mat(k,441) = (rxt(k,288)+rxt(k,293))*y(k,49)
         mat(k,1700) = -(4._r8*rxt(k,207)*y(k,5) + (rxt(k,208) + rxt(k,209) + rxt(k,210) &
                      ) * y(k,28) + rxt(k,211)*y(k,87) + rxt(k,212)*y(k,61) + rxt(k,213) &
                      *y(k,62) + rxt(k,215)*y(k,68) + rxt(k,216)*y(k,96) + rxt(k,266) &
                      *y(k,77))
         mat(k,1807) = -(rxt(k,208) + rxt(k,209) + rxt(k,210)) * y(k,5)
         mat(k,1388) = -rxt(k,211)*y(k,5)
         mat(k,1723) = -rxt(k,212)*y(k,5)
         mat(k,1420) = -rxt(k,213)*y(k,5)
         mat(k,1591) = -rxt(k,215)*y(k,5)
         mat(k,1361) = -rxt(k,216)*y(k,5)
         mat(k,889) = -rxt(k,266)*y(k,5)
         mat(k,994) = rxt(k,206)*y(k,70)
         mat(k,327) = rxt(k,214)*y(k,68)
         mat(k,546) = rxt(k,224)*y(k,92)
         mat(k,448) = rxt(k,219)*y(k,68)
         mat(k,1591) = mat(k,1591) + rxt(k,214)*y(k,6) + rxt(k,219)*y(k,53)
         mat(k,1778) = rxt(k,206)*y(k,3)
         mat(k,1471) = rxt(k,224)*y(k,45)
         mat(k,320) = -(rxt(k,214)*y(k,68))
         mat(k,1522) = -rxt(k,214)*y(k,6)
         mat(k,1681) = rxt(k,213)*y(k,62)
         mat(k,1396) = rxt(k,213)*y(k,5)
         mat(k,32) = -(rxt(k,193)*y(k,92))
         mat(k,1425) = -rxt(k,193)*y(k,8)
         mat(k,64) = -(rxt(k,194)*y(k,92))
         mat(k,1430) = -rxt(k,194)*y(k,9)
         mat(k,69) = -(rxt(k,220)*y(k,92))
         mat(k,1431) = -rxt(k,220)*y(k,10)
         mat(k,45) = -(rxt(k,195)*y(k,92))
         mat(k,1427) = -rxt(k,195)*y(k,11)
         mat(k,74) = -(rxt(k,196)*y(k,92))
         mat(k,1432) = -rxt(k,196)*y(k,12)
         mat(k,49) = -(rxt(k,197)*y(k,92))
         mat(k,1428) = -rxt(k,197)*y(k,13)
         mat(k,79) = -(rxt(k,198)*y(k,92))
         mat(k,1433) = -rxt(k,198)*y(k,14)
         mat(k,53) = -(rxt(k,199)*y(k,92))
         mat(k,1429) = -rxt(k,199)*y(k,15)
         mat(k,291) = -(rxt(k,231)*y(k,25) + rxt(k,232)*y(k,96) + rxt(k,243)*y(k,92))
         mat(k,1479) = -rxt(k,231)*y(k,16)
         mat(k,1327) = -rxt(k,232)*y(k,16)
         mat(k,1442) = -rxt(k,243)*y(k,16)
         mat(k,1610) = -(rxt(k,168)*y(k,25) + rxt(k,204)*y(k,3) + rxt(k,248)*y(k,63) &
                      + rxt(k,249)*y(k,68) + rxt(k,250)*y(k,96))
         mat(k,1499) = -rxt(k,168)*y(k,17)
         mat(k,993) = -rxt(k,204)*y(k,17)
         mat(k,1671) = -rxt(k,248)*y(k,17)
         mat(k,1588) = -rxt(k,249)*y(k,17)
         mat(k,1358) = -rxt(k,250)*y(k,17)
         mat(k,826) = rxt(k,175)*y(k,28) + rxt(k,252)*y(k,61)
         mat(k,238) = .300_r8*rxt(k,253)*y(k,96)
         mat(k,800) = (rxt(k,256)+rxt(k,257))*y(k,92)
         mat(k,1804) = rxt(k,175)*y(k,21)
         mat(k,1720) = rxt(k,252)*y(k,21)
         mat(k,1468) = (rxt(k,256)+rxt(k,257))*y(k,23)
         mat(k,1358) = mat(k,1358) + .300_r8*rxt(k,253)*y(k,22)
         mat(k,349) = -(rxt(k,233)*y(k,25) + rxt(k,234)*y(k,96) + rxt(k,244)*y(k,92))
         mat(k,1480) = -rxt(k,233)*y(k,18)
         mat(k,1328) = -rxt(k,234)*y(k,18)
         mat(k,1443) = -rxt(k,244)*y(k,18)
         mat(k,57) = -(rxt(k,235)*y(k,96))
         mat(k,1313) = -rxt(k,235)*y(k,19)
         mat(k,218) = -(rxt(k,236)*y(k,25) + rxt(k,237)*y(k,96))
         mat(k,1477) = -rxt(k,236)*y(k,20)
         mat(k,1322) = -rxt(k,237)*y(k,20)
         mat(k,819) = -(rxt(k,175)*y(k,28) + rxt(k,251)*y(k,87) + rxt(k,252)*y(k,61))
         mat(k,1790) = -rxt(k,175)*y(k,21)
         mat(k,1372) = -rxt(k,251)*y(k,21)
         mat(k,1705) = -rxt(k,252)*y(k,21)
         mat(k,235) = .700_r8*rxt(k,253)*y(k,96)
         mat(k,790) = rxt(k,169)*y(k,25) + rxt(k,225)*y(k,39) + rxt(k,255)*y(k,92) &
                      + rxt(k,254)*y(k,96)
         mat(k,1486) = rxt(k,169)*y(k,23)
         mat(k,566) = rxt(k,225)*y(k,23)
         mat(k,1449) = rxt(k,255)*y(k,23)
         mat(k,1339) = .700_r8*rxt(k,253)*y(k,22) + rxt(k,254)*y(k,23)
         mat(k,234) = -(rxt(k,253)*y(k,96))
         mat(k,1324) = -rxt(k,253)*y(k,22)
         mat(k,818) = rxt(k,251)*y(k,87)
         mat(k,1366) = rxt(k,251)*y(k,21)
         mat(k,789) = -(rxt(k,169)*y(k,25) + rxt(k,225)*y(k,39) + rxt(k,254)*y(k,96) &
                      + (rxt(k,255) + rxt(k,256) + rxt(k,257)) * y(k,92))
         mat(k,1485) = -rxt(k,169)*y(k,23)
         mat(k,565) = -rxt(k,225)*y(k,23)
         mat(k,1338) = -rxt(k,254)*y(k,23)
         mat(k,1448) = -(rxt(k,255) + rxt(k,256) + rxt(k,257)) * y(k,23)
         mat(k,246) = -(rxt(k,221)*y(k,92) + rxt(k,238)*y(k,25) + rxt(k,239)*y(k,96))
         mat(k,1441) = -rxt(k,221)*y(k,24)
         mat(k,1478) = -rxt(k,238)*y(k,24)
         mat(k,1325) = -rxt(k,239)*y(k,24)
         mat(k,1497) = -(rxt(k,168)*y(k,17) + rxt(k,169)*y(k,23) + rxt(k,170)*y(k,41) &
                      + rxt(k,171)*y(k,43) + (rxt(k,172) + rxt(k,173)) * y(k,87) &
                      + rxt(k,174)*y(k,70) + rxt(k,181)*y(k,29) + rxt(k,190)*y(k,54) &
                      + rxt(k,231)*y(k,16) + rxt(k,233)*y(k,18) + rxt(k,236)*y(k,20) &
                      + rxt(k,238)*y(k,24))
         mat(k,1608) = -rxt(k,168)*y(k,25)
         mat(k,798) = -rxt(k,169)*y(k,25)
         mat(k,1645) = -rxt(k,170)*y(k,25)
         mat(k,362) = -rxt(k,171)*y(k,25)
         mat(k,1383) = -(rxt(k,172) + rxt(k,173)) * y(k,25)
         mat(k,1773) = -rxt(k,174)*y(k,25)
         mat(k,587) = -rxt(k,181)*y(k,25)
         mat(k,552) = -rxt(k,190)*y(k,25)
         mat(k,296) = -rxt(k,231)*y(k,25)
         mat(k,356) = -rxt(k,233)*y(k,25)
         mat(k,224) = -rxt(k,236)*y(k,25)
         mat(k,251) = -rxt(k,238)*y(k,25)
         mat(k,1695) = rxt(k,209)*y(k,28)
         mat(k,34) = 4.000_r8*rxt(k,193)*y(k,92)
         mat(k,68) = rxt(k,194)*y(k,92)
         mat(k,48) = 2.000_r8*rxt(k,195)*y(k,92)
         mat(k,78) = 2.000_r8*rxt(k,196)*y(k,92)
         mat(k,52) = 2.000_r8*rxt(k,197)*y(k,92)
         mat(k,83) = rxt(k,198)*y(k,92)
         mat(k,56) = 2.000_r8*rxt(k,199)*y(k,92)
         mat(k,59) = 3.000_r8*rxt(k,235)*y(k,96)
         mat(k,224) = mat(k,224) + rxt(k,237)*y(k,96)
         mat(k,825) = rxt(k,175)*y(k,28)
         mat(k,1802) = rxt(k,209)*y(k,5) + rxt(k,175)*y(k,21) + (4.000_r8*rxt(k,176) &
                       +2.000_r8*rxt(k,178))*y(k,28) + rxt(k,180)*y(k,61) + rxt(k,185) &
                      *y(k,68) + rxt(k,267)*y(k,77) + rxt(k,186)*y(k,96)
         mat(k,107) = rxt(k,230)*y(k,92)
         mat(k,103) = rxt(k,245)*y(k,92) + rxt(k,240)*y(k,96)
         mat(k,95) = rxt(k,246)*y(k,92) + rxt(k,241)*y(k,96)
         mat(k,168) = rxt(k,247)*y(k,92) + rxt(k,242)*y(k,96)
         mat(k,958) = rxt(k,188)*y(k,68) + rxt(k,200)*y(k,92) + rxt(k,189)*y(k,96)
         mat(k,1718) = rxt(k,180)*y(k,28)
         mat(k,1586) = rxt(k,185)*y(k,28) + rxt(k,188)*y(k,49)
         mat(k,887) = rxt(k,267)*y(k,28)
         mat(k,1466) = 4.000_r8*rxt(k,193)*y(k,8) + rxt(k,194)*y(k,9) &
                      + 2.000_r8*rxt(k,195)*y(k,11) + 2.000_r8*rxt(k,196)*y(k,12) &
                      + 2.000_r8*rxt(k,197)*y(k,13) + rxt(k,198)*y(k,14) &
                      + 2.000_r8*rxt(k,199)*y(k,15) + rxt(k,230)*y(k,34) + rxt(k,245) &
                      *y(k,46) + rxt(k,246)*y(k,47) + rxt(k,247)*y(k,48) + rxt(k,200) &
                      *y(k,49)
         mat(k,1356) = 3.000_r8*rxt(k,235)*y(k,19) + rxt(k,237)*y(k,20) + rxt(k,186) &
                      *y(k,28) + rxt(k,240)*y(k,46) + rxt(k,241)*y(k,47) + rxt(k,242) &
                      *y(k,48) + rxt(k,189)*y(k,49)
         mat(k,1476) = rxt(k,181)*y(k,29)
         mat(k,1784) = 2.000_r8*rxt(k,177)*y(k,28)
         mat(k,580) = rxt(k,181)*y(k,25) + (rxt(k,286)+rxt(k,291)+rxt(k,296))*y(k,49)
         mat(k,946) = (rxt(k,286)+rxt(k,291)+rxt(k,296))*y(k,29) + (rxt(k,281) &
                       +rxt(k,287)+rxt(k,292))*y(k,54)
         mat(k,548) = (rxt(k,281)+rxt(k,287)+rxt(k,292))*y(k,49)
         mat(k,1783) = 2.000_r8*rxt(k,202)*y(k,28)
         mat(k,1810) = -(rxt(k,175)*y(k,21) + (4._r8*rxt(k,176) + 4._r8*rxt(k,177) &
                      + 4._r8*rxt(k,178) + 4._r8*rxt(k,202)) * y(k,28) + rxt(k,179) &
                      *y(k,87) + rxt(k,180)*y(k,61) + rxt(k,182)*y(k,62) + rxt(k,185) &
                      *y(k,68) + (rxt(k,186) + rxt(k,187)) * y(k,96) + (rxt(k,208) &
                      + rxt(k,209) + rxt(k,210)) * y(k,5) + rxt(k,267)*y(k,77))
         mat(k,828) = -rxt(k,175)*y(k,28)
         mat(k,1391) = -rxt(k,179)*y(k,28)
         mat(k,1726) = -rxt(k,180)*y(k,28)
         mat(k,1423) = -rxt(k,182)*y(k,28)
         mat(k,1594) = -rxt(k,185)*y(k,28)
         mat(k,1364) = -(rxt(k,186) + rxt(k,187)) * y(k,28)
         mat(k,1703) = -(rxt(k,208) + rxt(k,209) + rxt(k,210)) * y(k,28)
         mat(k,892) = -rxt(k,267)*y(k,28)
         mat(k,1505) = rxt(k,190)*y(k,54) + rxt(k,174)*y(k,70) + rxt(k,173)*y(k,87)
         mat(k,590) = rxt(k,183)*y(k,68)
         mat(k,962) = rxt(k,201)*y(k,92)
         mat(k,554) = rxt(k,190)*y(k,25) + rxt(k,191)*y(k,68) + rxt(k,192)*y(k,96)
         mat(k,1594) = mat(k,1594) + rxt(k,183)*y(k,29) + rxt(k,191)*y(k,54)
         mat(k,1781) = rxt(k,174)*y(k,25)
         mat(k,183) = rxt(k,272)*y(k,77)
         mat(k,892) = mat(k,892) + rxt(k,272)*y(k,71)
         mat(k,1391) = mat(k,1391) + rxt(k,173)*y(k,25)
         mat(k,1474) = rxt(k,201)*y(k,49)
         mat(k,1364) = mat(k,1364) + rxt(k,192)*y(k,54)
         mat(k,582) = -(rxt(k,181)*y(k,25) + rxt(k,183)*y(k,68) + rxt(k,184)*y(k,96) &
                      + (rxt(k,286) + rxt(k,291) + rxt(k,296)) * y(k,49))
         mat(k,1484) = -rxt(k,181)*y(k,29)
         mat(k,1547) = -rxt(k,183)*y(k,29)
         mat(k,1336) = -rxt(k,184)*y(k,29)
         mat(k,950) = -(rxt(k,286) + rxt(k,291) + rxt(k,296)) * y(k,29)
         mat(k,1789) = rxt(k,182)*y(k,62)
         mat(k,1399) = rxt(k,182)*y(k,28)
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
         mat(k,385) = -(rxt(k,258)*y(k,96))
         mat(k,1330) = -rxt(k,258)*y(k,31)
         mat(k,983) = rxt(k,204)*y(k,17)
         mat(k,1596) = rxt(k,204)*y(k,3) + rxt(k,168)*y(k,25) + rxt(k,248)*y(k,63) &
                      + rxt(k,249)*y(k,68) + rxt(k,250)*y(k,96)
         mat(k,219) = rxt(k,236)*y(k,25)
         mat(k,1482) = rxt(k,168)*y(k,17) + rxt(k,236)*y(k,20)
         mat(k,1126) = rxt(k,312)*y(k,97)
         mat(k,1657) = rxt(k,248)*y(k,17)
         mat(k,1526) = rxt(k,249)*y(k,17) + rxt(k,261)*y(k,72)
         mat(k,227) = rxt(k,261)*y(k,68) + rxt(k,262)*y(k,96)
         mat(k,1330) = mat(k,1330) + rxt(k,250)*y(k,17) + rxt(k,262)*y(k,72)
         mat(k,652) = rxt(k,312)*y(k,32)
         mat(k,1156) = -(rxt(k,312)*y(k,97) + rxt(k,331)*y(k,103) + rxt(k,339) &
                      *y(k,105) + rxt(k,340)*y(k,107) + rxt(k,384)*y(k,122) + rxt(k,403) &
                      *y(k,130) + rxt(k,440)*y(k,148) + rxt(k,442)*y(k,146) + rxt(k,467) &
                      *y(k,155) + rxt(k,488)*y(k,168) + rxt(k,499)*y(k,170) + rxt(k,503) &
                      *y(k,171) + rxt(k,507)*y(k,173) + rxt(k,512)*y(k,172))
         mat(k,659) = -rxt(k,312)*y(k,32)
         mat(k,368) = -rxt(k,331)*y(k,32)
         mat(k,432) = -rxt(k,339)*y(k,32)
         mat(k,636) = -rxt(k,340)*y(k,32)
         mat(k,534) = -rxt(k,384)*y(k,32)
         mat(k,674) = -rxt(k,403)*y(k,32)
         mat(k,647) = -rxt(k,440)*y(k,32)
         mat(k,710) = -rxt(k,442)*y(k,32)
         mat(k,687) = -rxt(k,467)*y(k,32)
         mat(k,381) = -rxt(k,488)*y(k,32)
         mat(k,773) = -rxt(k,499)*y(k,32)
         mat(k,737) = -rxt(k,503)*y(k,32)
         mat(k,761) = -rxt(k,507)*y(k,32)
         mat(k,749) = -rxt(k,512)*y(k,32)
         mat(k,386) = rxt(k,258)*y(k,96)
         mat(k,1200) = rxt(k,332)*y(k,104) + rxt(k,444)*y(k,151) + rxt(k,489)*y(k,169)
         mat(k,1578) = rxt(k,385)*y(k,127) + rxt(k,443)*y(k,151)
         mat(k,1297) = rxt(k,468)*y(k,163)
         mat(k,1045) = rxt(k,347)*y(k,108) + rxt(k,407)*y(k,135) + rxt(k,471)*y(k,163) &
                      + rxt(k,518)*y(k,174)
         mat(k,1348) = rxt(k,258)*y(k,31)
         mat(k,242) = rxt(k,332)*y(k,40)
         mat(k,192) = rxt(k,347)*y(k,86)
         mat(k,126) = rxt(k,385)*y(k,68)
         mat(k,143) = rxt(k,407)*y(k,86)
         mat(k,509) = rxt(k,444)*y(k,40) + rxt(k,443)*y(k,68) + rxt(k,448)*y(k,176)
         mat(k,481) = rxt(k,468)*y(k,69) + rxt(k,471)*y(k,86) + rxt(k,469)*y(k,176)
         mat(k,259) = rxt(k,489)*y(k,40)
         mat(k,490) = rxt(k,518)*y(k,86) + rxt(k,515)*y(k,176)
         mat(k,1869) = rxt(k,448)*y(k,151) + rxt(k,469)*y(k,163) + rxt(k,515)*y(k,174)
         mat(k,96) = -(rxt(k,229)*y(k,92))
         mat(k,1435) = -rxt(k,229)*y(k,33)
         mat(k,65) = rxt(k,194)*y(k,92)
         mat(k,70) = rxt(k,220)*y(k,92)
         mat(k,75) = rxt(k,196)*y(k,92)
         mat(k,50) = 2.000_r8*rxt(k,197)*y(k,92)
         mat(k,80) = 2.000_r8*rxt(k,198)*y(k,92)
         mat(k,54) = rxt(k,199)*y(k,92)
         mat(k,42) = 2.000_r8*rxt(k,222)*y(k,92)
         mat(k,92) = rxt(k,246)*y(k,92) + rxt(k,241)*y(k,96)
         mat(k,163) = rxt(k,247)*y(k,92) + rxt(k,242)*y(k,96)
         mat(k,1435) = mat(k,1435) + rxt(k,194)*y(k,9) + rxt(k,220)*y(k,10) &
                      + rxt(k,196)*y(k,12) + 2.000_r8*rxt(k,197)*y(k,13) &
                      + 2.000_r8*rxt(k,198)*y(k,14) + rxt(k,199)*y(k,15) &
                      + 2.000_r8*rxt(k,222)*y(k,42) + rxt(k,246)*y(k,47) + rxt(k,247) &
                      *y(k,48)
         mat(k,1315) = rxt(k,241)*y(k,47) + rxt(k,242)*y(k,48)
         mat(k,104) = -(rxt(k,230)*y(k,92))
         mat(k,1437) = -rxt(k,230)*y(k,34)
         mat(k,46) = rxt(k,195)*y(k,92)
         mat(k,76) = rxt(k,196)*y(k,92)
         mat(k,100) = rxt(k,245)*y(k,92) + rxt(k,240)*y(k,96)
         mat(k,1437) = mat(k,1437) + rxt(k,195)*y(k,11) + rxt(k,196)*y(k,12) &
                      + rxt(k,245)*y(k,46)
         mat(k,1317) = rxt(k,240)*y(k,46)
         mat(k,170) = -(rxt(k,259)*y(k,63) + (rxt(k,260) + rxt(k,274)) * y(k,96))
         mat(k,1656) = -rxt(k,259)*y(k,35)
         mat(k,1321) = -(rxt(k,260) + rxt(k,274)) * y(k,35)
         mat(k,564) = -(rxt(k,225)*y(k,23) + rxt(k,226)*y(k,41) + rxt(k,227)*y(k,176) &
                      + rxt(k,228)*y(k,51))
         mat(k,788) = -rxt(k,225)*y(k,39)
         mat(k,1623) = -rxt(k,226)*y(k,39)
         mat(k,1837) = -rxt(k,227)*y(k,39)
         mat(k,1213) = -rxt(k,228)*y(k,39)
         mat(k,71) = rxt(k,220)*y(k,92)
         mat(k,81) = rxt(k,198)*y(k,92)
         mat(k,97) = 2.000_r8*rxt(k,229)*y(k,92)
         mat(k,105) = rxt(k,230)*y(k,92)
         mat(k,1446) = rxt(k,220)*y(k,10) + rxt(k,198)*y(k,14) + 2.000_r8*rxt(k,229) &
                      *y(k,33) + rxt(k,230)*y(k,34)
         mat(k,1201) = -((rxt(k,126) + rxt(k,127) + rxt(k,128)) * y(k,87) + rxt(k,129) &
                      *y(k,69) + rxt(k,132)*y(k,70) + rxt(k,330)*y(k,103) + rxt(k,332) &
                      *y(k,104) + rxt(k,358)*y(k,114) + rxt(k,359)*y(k,113) + (rxt(k,360) &
                      + rxt(k,376)) * y(k,115) + rxt(k,388)*y(k,124) + rxt(k,389) &
                      *y(k,129) + rxt(k,390)*y(k,126) + rxt(k,425)*y(k,144) + rxt(k,444) &
                      *y(k,151) + rxt(k,446)*y(k,148) + rxt(k,449)*y(k,149) + rxt(k,450) &
                      *y(k,147) + rxt(k,451)*y(k,150) + rxt(k,486)*y(k,167) + rxt(k,487) &
                      *y(k,168) + rxt(k,489)*y(k,169))
         mat(k,1376) = -(rxt(k,126) + rxt(k,127) + rxt(k,128)) * y(k,40)
         mat(k,1298) = -rxt(k,129)*y(k,40)
         mat(k,1767) = -rxt(k,132)*y(k,40)
         mat(k,369) = -rxt(k,330)*y(k,40)
         mat(k,243) = -rxt(k,332)*y(k,40)
         mat(k,195) = -rxt(k,358)*y(k,40)
         mat(k,609) = -rxt(k,359)*y(k,40)
         mat(k,336) = -(rxt(k,360) + rxt(k,376)) * y(k,40)
         mat(k,206) = -rxt(k,388)*y(k,40)
         mat(k,134) = -rxt(k,389)*y(k,40)
         mat(k,401) = -rxt(k,390)*y(k,40)
         mat(k,464) = -rxt(k,425)*y(k,40)
         mat(k,510) = -rxt(k,444)*y(k,40)
         mat(k,648) = -rxt(k,446)*y(k,40)
         mat(k,255) = -rxt(k,449)*y(k,40)
         mat(k,698) = -rxt(k,450)*y(k,40)
         mat(k,524) = -rxt(k,451)*y(k,40)
         mat(k,516) = -rxt(k,486)*y(k,40)
         mat(k,382) = -rxt(k,487)*y(k,40)
         mat(k,260) = -rxt(k,489)*y(k,40)
         mat(k,1601) = rxt(k,250)*y(k,96)
         mat(k,793) = rxt(k,256)*y(k,92)
         mat(k,1490) = rxt(k,170)*y(k,41)
         mat(k,567) = rxt(k,226)*y(k,41)
         mat(k,1638) = rxt(k,170)*y(k,25) + rxt(k,226)*y(k,39) + rxt(k,124)*y(k,68) &
                      + rxt(k,107)*y(k,92) + rxt(k,133)*y(k,96) + rxt(k,328)*y(k,101) &
                      + rxt(k,423)*y(k,142) + rxt(k,485)*y(k,166)
         mat(k,542) = rxt(k,224)*y(k,92)
         mat(k,953) = rxt(k,201)*y(k,92)
         mat(k,869) = rxt(k,156)*y(k,96)
         mat(k,1579) = rxt(k,124)*y(k,41) + rxt(k,136)*y(k,96)
         mat(k,231) = rxt(k,262)*y(k,96)
         mat(k,406) = rxt(k,268)*y(k,96)
         mat(k,882) = rxt(k,273)*y(k,96)
         mat(k,1459) = rxt(k,256)*y(k,23) + rxt(k,107)*y(k,41) + rxt(k,224)*y(k,45) &
                      + rxt(k,201)*y(k,49)
         mat(k,1349) = rxt(k,250)*y(k,17) + rxt(k,133)*y(k,41) + rxt(k,156)*y(k,55) &
                      + rxt(k,136)*y(k,68) + rxt(k,262)*y(k,72) + rxt(k,268)*y(k,75) &
                      + rxt(k,273)*y(k,77) + rxt(k,414)*y(k,137)
         mat(k,596) = rxt(k,328)*y(k,41)
         mat(k,781) = rxt(k,414)*y(k,96)
         mat(k,855) = rxt(k,418)*y(k,176)
         mat(k,975) = rxt(k,423)*y(k,41)
         mat(k,721) = rxt(k,485)*y(k,41)
         mat(k,1870) = rxt(k,418)*y(k,141)
         mat(k,1648) = -(rxt(k,107)*y(k,92) + rxt(k,124)*y(k,68) + rxt(k,133)*y(k,96) &
                      + rxt(k,170)*y(k,25) + rxt(k,226)*y(k,39) + rxt(k,328)*y(k,101) &
                      + rxt(k,423)*y(k,142) + rxt(k,485)*y(k,166))
         mat(k,1469) = -rxt(k,107)*y(k,41)
         mat(k,1589) = -rxt(k,124)*y(k,41)
         mat(k,1359) = -rxt(k,133)*y(k,41)
         mat(k,1500) = -rxt(k,170)*y(k,41)
         mat(k,570) = -rxt(k,226)*y(k,41)
         mat(k,600) = -rxt(k,328)*y(k,41)
         mat(k,979) = -rxt(k,423)*y(k,41)
         mat(k,726) = -rxt(k,485)*y(k,41)
         mat(k,801) = rxt(k,257)*y(k,92)
         mat(k,1208) = rxt(k,126)*y(k,87) + rxt(k,376)*y(k,115)
         mat(k,1386) = rxt(k,126)*y(k,40)
         mat(k,1469) = mat(k,1469) + rxt(k,257)*y(k,23)
         mat(k,338) = rxt(k,376)*y(k,40)
         mat(k,41) = -(rxt(k,222)*y(k,92))
         mat(k,1426) = -rxt(k,222)*y(k,42)
         mat(k,358) = -(rxt(k,125)*y(k,68) + rxt(k,134)*y(k,96) + rxt(k,171)*y(k,25))
         mat(k,1524) = -rxt(k,125)*y(k,43)
         mat(k,1329) = -rxt(k,134)*y(k,43)
         mat(k,1481) = -rxt(k,171)*y(k,43)
         mat(k,1368) = 2.000_r8*rxt(k,140)*y(k,87)
         mat(k,1329) = mat(k,1329) + 2.000_r8*rxt(k,139)*y(k,96)
         mat(k,115) = rxt(k,275)*y(k,176)
         mat(k,1814) = rxt(k,275)*y(k,79)
         mat(k,540) = -(rxt(k,217)*y(k,68) + rxt(k,218)*y(k,96) + (rxt(k,223) &
                      + rxt(k,224)) * y(k,92))
         mat(k,1543) = -rxt(k,217)*y(k,45)
         mat(k,1332) = -rxt(k,218)*y(k,45)
         mat(k,1444) = -(rxt(k,223) + rxt(k,224)) * y(k,45)
         mat(k,984) = rxt(k,204)*y(k,17) + rxt(k,205)*y(k,87)
         mat(k,1597) = rxt(k,204)*y(k,3)
         mat(k,1370) = rxt(k,205)*y(k,3)
         mat(k,99) = -(rxt(k,240)*y(k,96) + rxt(k,245)*y(k,92))
         mat(k,1316) = -rxt(k,240)*y(k,46)
         mat(k,1436) = -rxt(k,245)*y(k,46)
         mat(k,91) = -(rxt(k,241)*y(k,96) + rxt(k,246)*y(k,92))
         mat(k,1314) = -rxt(k,241)*y(k,47)
         mat(k,1434) = -rxt(k,246)*y(k,47)
         mat(k,164) = -(rxt(k,242)*y(k,96) + rxt(k,247)*y(k,92))
         mat(k,1320) = -rxt(k,242)*y(k,48)
         mat(k,1440) = -rxt(k,247)*y(k,48)
         mat(k,951) = -(rxt(k,188)*y(k,68) + rxt(k,189)*y(k,96) + (rxt(k,200) &
                      + rxt(k,201)) * y(k,92) + (rxt(k,281) + rxt(k,287) + rxt(k,292) &
                      ) * y(k,54) + (rxt(k,286) + rxt(k,291) + rxt(k,296)) * y(k,29) &
                      + (rxt(k,288) + rxt(k,293)) * y(k,53))
         mat(k,1572) = -rxt(k,188)*y(k,49)
         mat(k,1342) = -rxt(k,189)*y(k,49)
         mat(k,1452) = -(rxt(k,200) + rxt(k,201)) * y(k,49)
         mat(k,550) = -(rxt(k,281) + rxt(k,287) + rxt(k,292)) * y(k,49)
         mat(k,583) = -(rxt(k,286) + rxt(k,291) + rxt(k,296)) * y(k,49)
         mat(k,443) = -(rxt(k,288) + rxt(k,293)) * y(k,49)
         mat(k,292) = rxt(k,231)*y(k,25)
         mat(k,1598) = rxt(k,168)*y(k,25)
         mat(k,351) = rxt(k,233)*y(k,25)
         mat(k,221) = 2.000_r8*rxt(k,236)*y(k,25)
         mat(k,791) = rxt(k,169)*y(k,25)
         mat(k,247) = rxt(k,238)*y(k,25)
         mat(k,1487) = rxt(k,231)*y(k,16) + rxt(k,168)*y(k,17) + rxt(k,233)*y(k,18) &
                      + 2.000_r8*rxt(k,236)*y(k,20) + rxt(k,169)*y(k,23) + rxt(k,238) &
                      *y(k,24) + rxt(k,170)*y(k,41) + rxt(k,171)*y(k,43) + rxt(k,190) &
                      *y(k,54) + rxt(k,172)*y(k,87)
         mat(k,1792) = rxt(k,187)*y(k,96)
         mat(k,1631) = rxt(k,170)*y(k,25)
         mat(k,359) = rxt(k,171)*y(k,25)
         mat(k,550) = mat(k,550) + rxt(k,190)*y(k,25)
         mat(k,1373) = rxt(k,172)*y(k,25)
         mat(k,1342) = mat(k,1342) + rxt(k,187)*y(k,28)
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
         mat(k,786) = rxt(k,225)*y(k,39)
         mat(k,563) = rxt(k,225)*y(k,23) + rxt(k,226)*y(k,41) + rxt(k,228)*y(k,51) &
                      + rxt(k,227)*y(k,176)
         mat(k,1618) = rxt(k,226)*y(k,39)
         mat(k,1212) = rxt(k,228)*y(k,39)
         mat(k,1819) = rxt(k,227)*y(k,39)
         mat(k,1220) = -(rxt(k,165)*y(k,96) + rxt(k,228)*y(k,39))
         mat(k,1350) = -rxt(k,165)*y(k,51)
         mat(k,568) = -rxt(k,228)*y(k,51)
         mat(k,1602) = rxt(k,248)*y(k,63)
         mat(k,584) = (rxt(k,286)+rxt(k,291)+rxt(k,296))*y(k,49)
         mat(k,172) = rxt(k,259)*y(k,63)
         mat(k,954) = (rxt(k,286)+rxt(k,291)+rxt(k,296))*y(k,29)
         mat(k,1409) = rxt(k,164)*y(k,96)
         mat(k,1663) = rxt(k,248)*y(k,17) + rxt(k,259)*y(k,35)
         mat(k,1350) = mat(k,1350) + rxt(k,164)*y(k,62)
         mat(k,276) = -(rxt(k,141)*y(k,96))
         mat(k,1326) = -rxt(k,141)*y(k,52)
         mat(k,1395) = rxt(k,162)*y(k,87)
         mat(k,1367) = rxt(k,162)*y(k,62)
         mat(k,442) = -(rxt(k,219)*y(k,68) + (rxt(k,288) + rxt(k,293)) * y(k,49))
         mat(k,1533) = -rxt(k,219)*y(k,53)
         mat(k,948) = -(rxt(k,288) + rxt(k,293)) * y(k,53)
         mat(k,1682) = rxt(k,211)*y(k,87)
         mat(k,1369) = rxt(k,211)*y(k,5)
         mat(k,549) = -(rxt(k,190)*y(k,25) + rxt(k,191)*y(k,68) + rxt(k,192)*y(k,96) &
                      + (rxt(k,281) + rxt(k,287) + rxt(k,292)) * y(k,49))
         mat(k,1483) = -rxt(k,190)*y(k,54)
         mat(k,1544) = -rxt(k,191)*y(k,54)
         mat(k,1333) = -rxt(k,192)*y(k,54)
         mat(k,949) = -(rxt(k,281) + rxt(k,287) + rxt(k,292)) * y(k,54)
         mat(k,1787) = rxt(k,179)*y(k,87)
         mat(k,581) = rxt(k,184)*y(k,96)
         mat(k,1371) = rxt(k,179)*y(k,28)
         mat(k,1333) = mat(k,1333) + rxt(k,184)*y(k,29)
         mat(k,864) = -(rxt(k,144)*y(k,61) + (rxt(k,145) + rxt(k,146) + rxt(k,147) &
                      ) * y(k,62) + rxt(k,148)*y(k,69) + rxt(k,156)*y(k,96) + rxt(k,309) &
                      *y(k,95))
         mat(k,1706) = -rxt(k,144)*y(k,55)
         mat(k,1400) = -(rxt(k,145) + rxt(k,146) + rxt(k,147)) * y(k,55)
         mat(k,1287) = -rxt(k,148)*y(k,55)
         mat(k,1340) = -rxt(k,156)*y(k,55)
         mat(k,1066) = -rxt(k,309)*y(k,55)
         mat(k,1567) = rxt(k,142)*y(k,88) + rxt(k,306)*y(k,91)
         mat(k,1287) = mat(k,1287) + rxt(k,307)*y(k,91)
         mat(k,1037) = 1.100_r8*rxt(k,302)*y(k,89) + .200_r8*rxt(k,300)*y(k,90)
         mat(k,436) = rxt(k,142)*y(k,68)
         mat(k,343) = 1.100_r8*rxt(k,302)*y(k,86)
         mat(k,1097) = .200_r8*rxt(k,300)*y(k,86)
         mat(k,499) = rxt(k,306)*y(k,68) + rxt(k,307)*y(k,69)
         mat(k,111) = -((rxt(k,160) + rxt(k,161)) * y(k,92))
         mat(k,1438) = -(rxt(k,160) + rxt(k,161)) * y(k,56)
         mat(k,861) = rxt(k,145)*y(k,62)
         mat(k,1393) = rxt(k,145)*y(k,55)
         mat(k,1394) = rxt(k,163)*y(k,63)
         mat(k,1655) = rxt(k,163)*y(k,62)
         mat(k,1724) = -(rxt(k,144)*y(k,55) + rxt(k,153)*y(k,63) + rxt(k,157)*y(k,87) &
                      + rxt(k,158)*y(k,70) + rxt(k,159)*y(k,68) + rxt(k,180)*y(k,28) &
                      + rxt(k,212)*y(k,5) + rxt(k,252)*y(k,21) + rxt(k,311)*y(k,95))
         mat(k,876) = -rxt(k,144)*y(k,61)
         mat(k,1675) = -rxt(k,153)*y(k,61)
         mat(k,1389) = -rxt(k,157)*y(k,61)
         mat(k,1779) = -rxt(k,158)*y(k,61)
         mat(k,1592) = -rxt(k,159)*y(k,61)
         mat(k,1808) = -rxt(k,180)*y(k,61)
         mat(k,1701) = -rxt(k,212)*y(k,61)
         mat(k,827) = -rxt(k,252)*y(k,61)
         mat(k,1083) = -rxt(k,311)*y(k,61)
         mat(k,876) = mat(k,876) + 2.000_r8*rxt(k,146)*y(k,62) + rxt(k,148)*y(k,69) &
                      + rxt(k,156)*y(k,96)
         mat(k,114) = 2.000_r8*rxt(k,160)*y(k,92)
         mat(k,1421) = 2.000_r8*rxt(k,146)*y(k,55) + rxt(k,149)*y(k,68) + rxt(k,269) &
                      *y(k,77)
         mat(k,1592) = mat(k,1592) + rxt(k,149)*y(k,62)
         mat(k,1309) = rxt(k,148)*y(k,55) + rxt(k,143)*y(k,88)
         mat(k,890) = rxt(k,269)*y(k,62)
         mat(k,440) = rxt(k,143)*y(k,69)
         mat(k,1114) = rxt(k,337)*y(k,100) + rxt(k,366)*y(k,110) + rxt(k,393)*y(k,121) &
                      + rxt(k,461)*y(k,145) + rxt(k,491)*y(k,165)
         mat(k,1472) = 2.000_r8*rxt(k,160)*y(k,56)
         mat(k,1362) = rxt(k,156)*y(k,55)
         mat(k,908) = rxt(k,337)*y(k,90)
         mat(k,815) = rxt(k,366)*y(k,90)
         mat(k,840) = rxt(k,393)*y(k,90)
         mat(k,943) = rxt(k,461)*y(k,90)
         mat(k,925) = rxt(k,491)*y(k,90)
         mat(k,1413) = -((rxt(k,145) + rxt(k,146) + rxt(k,147)) * y(k,55) + (rxt(k,149) &
                      + rxt(k,151)) * y(k,68) + rxt(k,150)*y(k,70) + rxt(k,162) &
                      *y(k,87) + rxt(k,163)*y(k,63) + rxt(k,164)*y(k,96) + rxt(k,182) &
                      *y(k,28) + rxt(k,213)*y(k,5) + rxt(k,269)*y(k,77))
         mat(k,873) = -(rxt(k,145) + rxt(k,146) + rxt(k,147)) * y(k,62)
         mat(k,1584) = -(rxt(k,149) + rxt(k,151)) * y(k,62)
         mat(k,1771) = -rxt(k,150)*y(k,62)
         mat(k,1381) = -rxt(k,162)*y(k,62)
         mat(k,1667) = -rxt(k,163)*y(k,62)
         mat(k,1354) = -rxt(k,164)*y(k,62)
         mat(k,1800) = -rxt(k,182)*y(k,62)
         mat(k,1693) = -rxt(k,213)*y(k,62)
         mat(k,886) = -rxt(k,269)*y(k,62)
         mat(k,1693) = mat(k,1693) + rxt(k,212)*y(k,61)
         mat(k,824) = rxt(k,252)*y(k,61)
         mat(k,1800) = mat(k,1800) + rxt(k,180)*y(k,61)
         mat(k,280) = rxt(k,141)*y(k,96)
         mat(k,1716) = rxt(k,212)*y(k,5) + rxt(k,252)*y(k,21) + rxt(k,180)*y(k,28) &
                      + 2.000_r8*rxt(k,153)*y(k,63) + rxt(k,159)*y(k,68) + rxt(k,158) &
                      *y(k,70) + rxt(k,157)*y(k,87)
         mat(k,1667) = mat(k,1667) + 2.000_r8*rxt(k,153)*y(k,61) + rxt(k,154)*y(k,68) &
                      + rxt(k,152)*y(k,87) + rxt(k,155)*y(k,96)
         mat(k,1584) = mat(k,1584) + rxt(k,159)*y(k,61) + rxt(k,154)*y(k,63)
         mat(k,1771) = mat(k,1771) + rxt(k,158)*y(k,61)
         mat(k,1381) = mat(k,1381) + rxt(k,157)*y(k,61) + rxt(k,152)*y(k,63)
         mat(k,1354) = mat(k,1354) + rxt(k,141)*y(k,52) + rxt(k,155)*y(k,63)
         mat(k,1673) = -(rxt(k,152)*y(k,87) + rxt(k,153)*y(k,61) + rxt(k,154)*y(k,68) &
                      + rxt(k,155)*y(k,96) + rxt(k,163)*y(k,62) + rxt(k,248)*y(k,17) &
                      + rxt(k,259)*y(k,35))
         mat(k,1387) = -rxt(k,152)*y(k,63)
         mat(k,1722) = -rxt(k,153)*y(k,63)
         mat(k,1590) = -rxt(k,154)*y(k,63)
         mat(k,1360) = -rxt(k,155)*y(k,63)
         mat(k,1419) = -rxt(k,163)*y(k,63)
         mat(k,1612) = -rxt(k,248)*y(k,63)
         mat(k,175) = -rxt(k,259)*y(k,63)
         mat(k,326) = rxt(k,214)*y(k,68)
         mat(k,1501) = rxt(k,181)*y(k,29)
         mat(k,589) = rxt(k,181)*y(k,25) + rxt(k,183)*y(k,68) + rxt(k,184)*y(k,96)
         mat(k,571) = rxt(k,228)*y(k,51)
         mat(k,1230) = rxt(k,228)*y(k,39) + rxt(k,165)*y(k,96)
         mat(k,1419) = mat(k,1419) + rxt(k,151)*y(k,68) + rxt(k,150)*y(k,70)
         mat(k,1590) = mat(k,1590) + rxt(k,214)*y(k,6) + rxt(k,183)*y(k,29) &
                      + rxt(k,151)*y(k,62)
         mat(k,1777) = rxt(k,150)*y(k,62)
         mat(k,1360) = mat(k,1360) + rxt(k,184)*y(k,29) + rxt(k,165)*y(k,51)
         mat(k,1587) = -(rxt(k,121)*y(k,70) + 4._r8*rxt(k,122)*y(k,68) + rxt(k,123) &
                      *y(k,69) + rxt(k,124)*y(k,41) + rxt(k,125)*y(k,43) + rxt(k,130) &
                      *y(k,87) + rxt(k,136)*y(k,96) + (rxt(k,149) + rxt(k,151) &
                      ) * y(k,62) + rxt(k,154)*y(k,63) + rxt(k,159)*y(k,61) + rxt(k,183) &
                      *y(k,29) + rxt(k,185)*y(k,28) + rxt(k,188)*y(k,49) + rxt(k,191) &
                      *y(k,54) + rxt(k,214)*y(k,6) + rxt(k,215)*y(k,5) + rxt(k,217) &
                      *y(k,45) + rxt(k,219)*y(k,53) + rxt(k,249)*y(k,17) + rxt(k,261) &
                      *y(k,72) + (rxt(k,304) + rxt(k,305)) * y(k,89) + rxt(k,306) &
                      *y(k,91) + rxt(k,326)*y(k,101) + rxt(k,334)*y(k,102) + rxt(k,342) &
                      *y(k,107) + rxt(k,343)*y(k,106) + rxt(k,352)*y(k,111) + rxt(k,354) &
                      *y(k,112) + rxt(k,356)*y(k,113) + rxt(k,362)*y(k,118) + rxt(k,363) &
                      *y(k,120) + rxt(k,364)*y(k,119) + rxt(k,378)*y(k,122) + rxt(k,380) &
                      *y(k,123) + rxt(k,385)*y(k,127) + rxt(k,387)*y(k,128) + rxt(k,395) &
                      *y(k,131) + rxt(k,399)*y(k,132) + rxt(k,424)*y(k,144) + rxt(k,430) &
                      *y(k,142) + rxt(k,434)*y(k,146) + rxt(k,437)*y(k,147) + rxt(k,443) &
                      *y(k,151) + rxt(k,445)*y(k,148) + rxt(k,453)*y(k,153) + rxt(k,454) &
                      *y(k,152) + rxt(k,463)*y(k,160) + rxt(k,465)*y(k,161) + rxt(k,480) &
                      *y(k,166) + rxt(k,483)*y(k,167) + rxt(k,497)*y(k,170) + rxt(k,501) &
                      *y(k,171) + rxt(k,505)*y(k,172))
         mat(k,1774) = -rxt(k,121)*y(k,68)
         mat(k,1305) = -rxt(k,123)*y(k,68)
         mat(k,1646) = -rxt(k,124)*y(k,68)
         mat(k,363) = -rxt(k,125)*y(k,68)
         mat(k,1384) = -rxt(k,130)*y(k,68)
         mat(k,1357) = -rxt(k,136)*y(k,68)
         mat(k,1416) = -(rxt(k,149) + rxt(k,151)) * y(k,68)
         mat(k,1670) = -rxt(k,154)*y(k,68)
         mat(k,1719) = -rxt(k,159)*y(k,68)
         mat(k,588) = -rxt(k,183)*y(k,68)
         mat(k,1803) = -rxt(k,185)*y(k,68)
         mat(k,959) = -rxt(k,188)*y(k,68)
         mat(k,553) = -rxt(k,191)*y(k,68)
         mat(k,325) = -rxt(k,214)*y(k,68)
         mat(k,1696) = -rxt(k,215)*y(k,68)
         mat(k,545) = -rxt(k,217)*y(k,68)
         mat(k,447) = -rxt(k,219)*y(k,68)
         mat(k,1609) = -rxt(k,249)*y(k,68)
         mat(k,233) = -rxt(k,261)*y(k,68)
         mat(k,348) = -(rxt(k,304) + rxt(k,305)) * y(k,68)
         mat(k,503) = -rxt(k,306)*y(k,68)
         mat(k,599) = -rxt(k,326)*y(k,68)
         mat(k,188) = -rxt(k,334)*y(k,68)
         mat(k,638) = -rxt(k,342)*y(k,68)
         mat(k,395) = -rxt(k,343)*y(k,68)
         mat(k,625) = -rxt(k,352)*y(k,68)
         mat(k,302) = -rxt(k,354)*y(k,68)
         mat(k,612) = -rxt(k,356)*y(k,68)
         mat(k,377) = -rxt(k,362)*y(k,68)
         mat(k,201) = -rxt(k,363)*y(k,68)
         mat(k,309) = -rxt(k,364)*y(k,68)
         mat(k,537) = -rxt(k,378)*y(k,68)
         mat(k,315) = -rxt(k,380)*y(k,68)
         mat(k,127) = -rxt(k,385)*y(k,68)
         mat(k,131) = -rxt(k,387)*y(k,68)
         mat(k,417) = -rxt(k,395)*y(k,68)
         mat(k,286) = -rxt(k,399)*y(k,68)
         mat(k,467) = -rxt(k,424)*y(k,68)
         mat(k,978) = -rxt(k,430)*y(k,68)
         mat(k,713) = -rxt(k,434)*y(k,68)
         mat(k,700) = -rxt(k,437)*y(k,68)
         mat(k,511) = -rxt(k,443)*y(k,68)
         mat(k,650) = -rxt(k,445)*y(k,68)
         mat(k,149) = -rxt(k,453)*y(k,68)
         mat(k,90) = -rxt(k,454)*y(k,68)
         mat(k,473) = -rxt(k,463)*y(k,68)
         mat(k,290) = -rxt(k,465)*y(k,68)
         mat(k,725) = -rxt(k,480)*y(k,68)
         mat(k,519) = -rxt(k,483)*y(k,68)
         mat(k,775) = -rxt(k,497)*y(k,68)
         mat(k,739) = -rxt(k,501)*y(k,68)
         mat(k,751) = -rxt(k,505)*y(k,68)
         mat(k,1163) = rxt(k,507)*y(k,173)
         mat(k,1207) = rxt(k,128)*y(k,87) + rxt(k,425)*y(k,144) + rxt(k,450)*y(k,147)
         mat(k,875) = rxt(k,144)*y(k,61) + rxt(k,145)*y(k,62) + rxt(k,148)*y(k,69) &
                      + rxt(k,309)*y(k,95)
         mat(k,1719) = mat(k,1719) + rxt(k,144)*y(k,55)
         mat(k,1416) = mat(k,1416) + rxt(k,145)*y(k,55)
         mat(k,1305) = mat(k,1305) + rxt(k,148)*y(k,55) + rxt(k,263)*y(k,75) &
                      + rxt(k,270)*y(k,77) + rxt(k,308)*y(k,91) + (rxt(k,110) &
                       +rxt(k,111))*y(k,92) + rxt(k,315)*y(k,97) + rxt(k,319)*y(k,98) &
                      + rxt(k,409)*y(k,136) + rxt(k,506)*y(k,173)
         mat(k,409) = rxt(k,263)*y(k,69)
         mat(k,888) = rxt(k,270)*y(k,69)
         mat(k,1052) = rxt(k,300)*y(k,90) + 1.150_r8*rxt(k,301)*y(k,95) + rxt(k,349) &
                      *y(k,106) + rxt(k,373)*y(k,118) + rxt(k,405)*y(k,131) &
                      + rxt(k,426)*y(k,142) + rxt(k,473)*y(k,160) + rxt(k,520) &
                      *y(k,173)
         mat(k,1384) = mat(k,1384) + rxt(k,128)*y(k,40)
         mat(k,439) = rxt(k,314)*y(k,97)
         mat(k,1112) = rxt(k,300)*y(k,86)
         mat(k,503) = mat(k,503) + rxt(k,308)*y(k,69)
         mat(k,1467) = (rxt(k,110)+rxt(k,111))*y(k,69)
         mat(k,562) = rxt(k,419)*y(k,141) + rxt(k,433)*y(k,145)
         mat(k,1081) = rxt(k,309)*y(k,55) + 1.150_r8*rxt(k,301)*y(k,86)
         mat(k,1357) = mat(k,1357) + 2.000_r8*rxt(k,138)*y(k,96)
         mat(k,664) = rxt(k,315)*y(k,69) + rxt(k,314)*y(k,88)
         mat(k,275) = rxt(k,319)*y(k,69)
         mat(k,395) = mat(k,395) + rxt(k,349)*y(k,86)
         mat(k,377) = mat(k,377) + rxt(k,373)*y(k,86)
         mat(k,417) = mat(k,417) + rxt(k,405)*y(k,86)
         mat(k,495) = rxt(k,409)*y(k,69)
         mat(k,858) = rxt(k,419)*y(k,93)
         mat(k,978) = mat(k,978) + rxt(k,426)*y(k,86)
         mat(k,467) = mat(k,467) + rxt(k,425)*y(k,40)
         mat(k,942) = rxt(k,433)*y(k,93)
         mat(k,700) = mat(k,700) + rxt(k,450)*y(k,40)
         mat(k,473) = mat(k,473) + rxt(k,473)*y(k,86)
         mat(k,763) = rxt(k,507)*y(k,32) + rxt(k,506)*y(k,69) + rxt(k,520)*y(k,86) &
                      + rxt(k,509)*y(k,176)
         mat(k,1878) = rxt(k,509)*y(k,173)
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
         mat(k,1299) = -(rxt(k,110)*y(k,92) + rxt(k,115)*y(k,93) + rxt(k,123)*y(k,68) &
                      + rxt(k,129)*y(k,40) + rxt(k,143)*y(k,88) + rxt(k,148)*y(k,55) &
                      + rxt(k,263)*y(k,75) + rxt(k,270)*y(k,77) + rxt(k,303)*y(k,89) &
                      + (rxt(k,307) + rxt(k,308)) * y(k,91) + rxt(k,315)*y(k,97) &
                      + rxt(k,319)*y(k,98) + rxt(k,333)*y(k,100) + rxt(k,345)*y(k,106) &
                      + rxt(k,368)*y(k,111) + rxt(k,371)*y(k,117) + rxt(k,383) &
                      *y(k,122) + rxt(k,386)*y(k,123) + rxt(k,397)*y(k,134) + rxt(k,398) &
                      *y(k,130) + rxt(k,409)*y(k,136) + rxt(k,420)*y(k,141) + rxt(k,436) &
                      *y(k,145) + rxt(k,441)*y(k,146) + rxt(k,452)*y(k,150) + rxt(k,464) &
                      *y(k,155) + rxt(k,468)*y(k,163) + rxt(k,470)*y(k,162) + rxt(k,478) &
                      *y(k,164) + rxt(k,481)*y(k,165) + rxt(k,495)*y(k,170) + rxt(k,506) &
                      *y(k,173) + rxt(k,511)*y(k,171))
         mat(k,1461) = -rxt(k,110)*y(k,69)
         mat(k,561) = -rxt(k,115)*y(k,69)
         mat(k,1581) = -rxt(k,123)*y(k,69)
         mat(k,1202) = -rxt(k,129)*y(k,69)
         mat(k,437) = -rxt(k,143)*y(k,69)
         mat(k,870) = -rxt(k,148)*y(k,69)
         mat(k,407) = -rxt(k,263)*y(k,69)
         mat(k,883) = -rxt(k,270)*y(k,69)
         mat(k,347) = -rxt(k,303)*y(k,69)
         mat(k,502) = -(rxt(k,307) + rxt(k,308)) * y(k,69)
         mat(k,660) = -rxt(k,315)*y(k,69)
         mat(k,274) = -rxt(k,319)*y(k,69)
         mat(k,904) = -rxt(k,333)*y(k,69)
         mat(k,394) = -rxt(k,345)*y(k,69)
         mat(k,623) = -rxt(k,368)*y(k,69)
         mat(k,456) = -rxt(k,371)*y(k,69)
         mat(k,536) = -rxt(k,383)*y(k,69)
         mat(k,314) = -rxt(k,386)*y(k,69)
         mat(k,212) = -rxt(k,397)*y(k,69)
         mat(k,675) = -rxt(k,398)*y(k,69)
         mat(k,494) = -rxt(k,409)*y(k,69)
         mat(k,856) = -rxt(k,420)*y(k,69)
         mat(k,940) = -rxt(k,436)*y(k,69)
         mat(k,712) = -rxt(k,441)*y(k,69)
         mat(k,525) = -rxt(k,452)*y(k,69)
         mat(k,688) = -rxt(k,464)*y(k,69)
         mat(k,482) = -rxt(k,468)*y(k,69)
         mat(k,424) = -rxt(k,470)*y(k,69)
         mat(k,217) = -rxt(k,478)*y(k,69)
         mat(k,920) = -rxt(k,481)*y(k,69)
         mat(k,774) = -rxt(k,495)*y(k,69)
         mat(k,762) = -rxt(k,506)*y(k,69)
         mat(k,738) = -rxt(k,511)*y(k,69)
         mat(k,988) = rxt(k,206)*y(k,70) + rxt(k,205)*y(k,87)
         mat(k,1690) = 2.000_r8*rxt(k,207)*y(k,5) + (rxt(k,209)+rxt(k,210))*y(k,28) &
                      + rxt(k,215)*y(k,68) + rxt(k,211)*y(k,87)
         mat(k,821) = rxt(k,251)*y(k,87)
         mat(k,1492) = rxt(k,174)*y(k,70) + rxt(k,172)*y(k,87)
         mat(k,1797) = (rxt(k,209)+rxt(k,210))*y(k,5) + (2.000_r8*rxt(k,176) &
                       +2.000_r8*rxt(k,177))*y(k,28) + rxt(k,185)*y(k,68) + rxt(k,179) &
                      *y(k,87) + rxt(k,187)*y(k,96)
         mat(k,1158) = rxt(k,440)*y(k,148) + rxt(k,512)*y(k,172)
         mat(k,1202) = mat(k,1202) + rxt(k,132)*y(k,70) + rxt(k,126)*y(k,87) &
                      + rxt(k,359)*y(k,113) + rxt(k,390)*y(k,126) + rxt(k,446) &
                      *y(k,148)
         mat(k,277) = rxt(k,141)*y(k,96)
         mat(k,870) = mat(k,870) + rxt(k,147)*y(k,62)
         mat(k,112) = rxt(k,161)*y(k,92)
         mat(k,1713) = rxt(k,158)*y(k,70) + rxt(k,311)*y(k,95)
         mat(k,1410) = rxt(k,147)*y(k,55) + rxt(k,149)*y(k,68) + rxt(k,150)*y(k,70)
         mat(k,1664) = rxt(k,154)*y(k,68) + rxt(k,152)*y(k,87)
         mat(k,1581) = mat(k,1581) + rxt(k,215)*y(k,5) + rxt(k,185)*y(k,28) &
                      + rxt(k,149)*y(k,62) + rxt(k,154)*y(k,63) + 2.000_r8*rxt(k,122) &
                      *y(k,68) + 2.000_r8*rxt(k,121)*y(k,70) + rxt(k,130)*y(k,87) &
                      + rxt(k,114)*y(k,93) + rxt(k,136)*y(k,96) + rxt(k,326)*y(k,101) &
                      + rxt(k,334)*y(k,102) + rxt(k,343)*y(k,106) + rxt(k,352) &
                      *y(k,111) + rxt(k,354)*y(k,112) + rxt(k,356)*y(k,113) &
                      + rxt(k,362)*y(k,118) + rxt(k,364)*y(k,119) + rxt(k,378) &
                      *y(k,122) + rxt(k,380)*y(k,123) + rxt(k,387)*y(k,128) &
                      + rxt(k,395)*y(k,131) + rxt(k,399)*y(k,132) + rxt(k,430) &
                      *y(k,142) + rxt(k,424)*y(k,144) + rxt(k,434)*y(k,146) &
                      + rxt(k,437)*y(k,147) + rxt(k,445)*y(k,148) + rxt(k,454) &
                      *y(k,152) + rxt(k,453)*y(k,153) + rxt(k,463)*y(k,160) &
                      + rxt(k,465)*y(k,161) + rxt(k,480)*y(k,166) + rxt(k,483) &
                      *y(k,167) + rxt(k,505)*y(k,172)
         mat(k,1299) = mat(k,1299) + 2.000_r8*rxt(k,115)*y(k,93)
         mat(k,1768) = rxt(k,206)*y(k,3) + rxt(k,174)*y(k,25) + rxt(k,132)*y(k,40) &
                      + rxt(k,158)*y(k,61) + rxt(k,150)*y(k,62) + 2.000_r8*rxt(k,121) &
                      *y(k,68) + rxt(k,265)*y(k,75) + rxt(k,271)*y(k,77) &
                      + 2.000_r8*rxt(k,131)*y(k,87) + 2.000_r8*rxt(k,112)*y(k,92) &
                      + rxt(k,137)*y(k,96) + rxt(k,325)*y(k,100) + 2.000_r8*rxt(k,327) &
                      *y(k,101) + rxt(k,351)*y(k,110) + rxt(k,353)*y(k,111) &
                      + rxt(k,355)*y(k,112) + rxt(k,361)*y(k,117) + rxt(k,377) &
                      *y(k,121) + rxt(k,379)*y(k,122) + rxt(k,394)*y(k,130) + ( &
                      + 2.000_r8*rxt(k,400)+rxt(k,401))*y(k,131) + rxt(k,410)*y(k,136) &
                      + rxt(k,411)*y(k,137) + rxt(k,416)*y(k,141) + rxt(k,432) &
                      *y(k,145) + rxt(k,435)*y(k,146) + rxt(k,447)*y(k,147) &
                      + rxt(k,462)*y(k,155) + rxt(k,479)*y(k,165) + rxt(k,482) &
                      *y(k,166)
         mat(k,407) = mat(k,407) + rxt(k,265)*y(k,70)
         mat(k,883) = mat(k,883) + rxt(k,271)*y(k,70)
         mat(k,1047) = rxt(k,374)*y(k,119) + rxt(k,428)*y(k,144) + rxt(k,517)*y(k,172)
         mat(k,1378) = rxt(k,205)*y(k,3) + rxt(k,211)*y(k,5) + rxt(k,251)*y(k,21) &
                      + rxt(k,172)*y(k,25) + rxt(k,179)*y(k,28) + rxt(k,126)*y(k,40) &
                      + rxt(k,152)*y(k,63) + rxt(k,130)*y(k,68) + 2.000_r8*rxt(k,131) &
                      *y(k,70) + 2.000_r8*rxt(k,140)*y(k,87) + rxt(k,135)*y(k,96)
         mat(k,1461) = mat(k,1461) + rxt(k,161)*y(k,56) + 2.000_r8*rxt(k,112)*y(k,70)
         mat(k,561) = mat(k,561) + rxt(k,114)*y(k,68) + 2.000_r8*rxt(k,115)*y(k,69)
         mat(k,1076) = rxt(k,311)*y(k,61) + rxt(k,336)*y(k,100) + rxt(k,365)*y(k,110) &
                      + rxt(k,392)*y(k,121) + rxt(k,429)*y(k,137) + rxt(k,460) &
                      *y(k,145) + rxt(k,492)*y(k,165)
         mat(k,1351) = rxt(k,187)*y(k,28) + rxt(k,141)*y(k,52) + rxt(k,136)*y(k,68) &
                      + rxt(k,137)*y(k,70) + rxt(k,135)*y(k,87)
         mat(k,904) = mat(k,904) + rxt(k,325)*y(k,70) + rxt(k,336)*y(k,95)
         mat(k,597) = rxt(k,326)*y(k,68) + 2.000_r8*rxt(k,327)*y(k,70)
         mat(k,187) = rxt(k,334)*y(k,68)
         mat(k,394) = mat(k,394) + rxt(k,343)*y(k,68)
         mat(k,811) = rxt(k,351)*y(k,70) + rxt(k,365)*y(k,95)
         mat(k,623) = mat(k,623) + rxt(k,352)*y(k,68) + rxt(k,353)*y(k,70)
         mat(k,301) = rxt(k,354)*y(k,68) + rxt(k,355)*y(k,70)
         mat(k,610) = rxt(k,359)*y(k,40) + rxt(k,356)*y(k,68) + rxt(k,357)*y(k,176)
         mat(k,456) = mat(k,456) + rxt(k,361)*y(k,70)
         mat(k,376) = rxt(k,362)*y(k,68)
         mat(k,308) = rxt(k,364)*y(k,68) + rxt(k,374)*y(k,86)
         mat(k,838) = rxt(k,377)*y(k,70) + rxt(k,392)*y(k,95)
         mat(k,536) = mat(k,536) + rxt(k,378)*y(k,68) + rxt(k,379)*y(k,70)
         mat(k,314) = mat(k,314) + rxt(k,380)*y(k,68)
         mat(k,402) = rxt(k,390)*y(k,40) + rxt(k,382)*y(k,176)
         mat(k,130) = rxt(k,387)*y(k,68)
         mat(k,675) = mat(k,675) + rxt(k,394)*y(k,70)
         mat(k,416) = rxt(k,395)*y(k,68) + (2.000_r8*rxt(k,400)+rxt(k,401))*y(k,70)
         mat(k,285) = rxt(k,399)*y(k,68)
         mat(k,494) = mat(k,494) + rxt(k,410)*y(k,70)
         mat(k,782) = rxt(k,411)*y(k,70) + rxt(k,429)*y(k,95)
         mat(k,856) = mat(k,856) + rxt(k,416)*y(k,70)
         mat(k,976) = rxt(k,430)*y(k,68)
         mat(k,465) = rxt(k,424)*y(k,68) + rxt(k,428)*y(k,86)
         mat(k,940) = mat(k,940) + rxt(k,432)*y(k,70) + rxt(k,460)*y(k,95)
         mat(k,712) = mat(k,712) + rxt(k,434)*y(k,68) + rxt(k,435)*y(k,70)
         mat(k,699) = rxt(k,437)*y(k,68) + rxt(k,447)*y(k,70)
         mat(k,649) = rxt(k,440)*y(k,32) + rxt(k,446)*y(k,40) + rxt(k,445)*y(k,68) &
                      + rxt(k,439)*y(k,176)
         mat(k,89) = rxt(k,454)*y(k,68)
         mat(k,148) = rxt(k,453)*y(k,68)
         mat(k,688) = mat(k,688) + rxt(k,462)*y(k,70)
         mat(k,472) = rxt(k,463)*y(k,68)
         mat(k,289) = rxt(k,465)*y(k,68)
         mat(k,920) = mat(k,920) + rxt(k,479)*y(k,70) + rxt(k,492)*y(k,95)
         mat(k,722) = rxt(k,480)*y(k,68) + rxt(k,482)*y(k,70)
         mat(k,517) = rxt(k,483)*y(k,68)
         mat(k,750) = rxt(k,512)*y(k,32) + rxt(k,505)*y(k,68) + rxt(k,517)*y(k,86) &
                      + rxt(k,513)*y(k,176)
         mat(k,1872) = rxt(k,357)*y(k,113) + rxt(k,382)*y(k,126) + rxt(k,439)*y(k,148) &
                      + rxt(k,513)*y(k,172)
         mat(k,1780) = -(rxt(k,112)*y(k,92) + rxt(k,121)*y(k,68) + rxt(k,131)*y(k,87) &
                      + rxt(k,132)*y(k,40) + rxt(k,137)*y(k,96) + rxt(k,150)*y(k,62) &
                      + rxt(k,158)*y(k,61) + rxt(k,174)*y(k,25) + rxt(k,206)*y(k,3) &
                      + rxt(k,265)*y(k,75) + rxt(k,271)*y(k,77) + rxt(k,325)*y(k,100) &
                      + rxt(k,327)*y(k,101) + rxt(k,351)*y(k,110) + rxt(k,353) &
                      *y(k,111) + rxt(k,355)*y(k,112) + rxt(k,361)*y(k,117) + rxt(k,377) &
                      *y(k,121) + rxt(k,379)*y(k,122) + rxt(k,394)*y(k,130) + (rxt(k,400) &
                      + rxt(k,401)) * y(k,131) + rxt(k,410)*y(k,136) + rxt(k,411) &
                      *y(k,137) + (rxt(k,416) + rxt(k,417)) * y(k,141) + rxt(k,421) &
                      *y(k,142) + rxt(k,432)*y(k,145) + rxt(k,435)*y(k,146) + rxt(k,447) &
                      *y(k,147) + rxt(k,462)*y(k,155) + rxt(k,479)*y(k,165) + rxt(k,482) &
                      *y(k,166))
         mat(k,1473) = -rxt(k,112)*y(k,70)
         mat(k,1593) = -rxt(k,121)*y(k,70)
         mat(k,1390) = -rxt(k,131)*y(k,70)
         mat(k,1210) = -rxt(k,132)*y(k,70)
         mat(k,1363) = -rxt(k,137)*y(k,70)
         mat(k,1422) = -rxt(k,150)*y(k,70)
         mat(k,1725) = -rxt(k,158)*y(k,70)
         mat(k,1504) = -rxt(k,174)*y(k,70)
         mat(k,995) = -rxt(k,206)*y(k,70)
         mat(k,410) = -rxt(k,265)*y(k,70)
         mat(k,891) = -rxt(k,271)*y(k,70)
         mat(k,909) = -rxt(k,325)*y(k,70)
         mat(k,601) = -rxt(k,327)*y(k,70)
         mat(k,816) = -rxt(k,351)*y(k,70)
         mat(k,627) = -rxt(k,353)*y(k,70)
         mat(k,303) = -rxt(k,355)*y(k,70)
         mat(k,458) = -rxt(k,361)*y(k,70)
         mat(k,841) = -rxt(k,377)*y(k,70)
         mat(k,538) = -rxt(k,379)*y(k,70)
         mat(k,677) = -rxt(k,394)*y(k,70)
         mat(k,418) = -(rxt(k,400) + rxt(k,401)) * y(k,70)
         mat(k,496) = -rxt(k,410)*y(k,70)
         mat(k,784) = -rxt(k,411)*y(k,70)
         mat(k,859) = -(rxt(k,416) + rxt(k,417)) * y(k,70)
         mat(k,981) = -rxt(k,421)*y(k,70)
         mat(k,944) = -rxt(k,432)*y(k,70)
         mat(k,714) = -rxt(k,435)*y(k,70)
         mat(k,701) = -rxt(k,447)*y(k,70)
         mat(k,690) = -rxt(k,462)*y(k,70)
         mat(k,926) = -rxt(k,479)*y(k,70)
         mat(k,727) = -rxt(k,482)*y(k,70)
         mat(k,1593) = mat(k,1593) + rxt(k,123)*y(k,69)
         mat(k,1310) = rxt(k,123)*y(k,68) + rxt(k,345)*y(k,106)
         mat(k,396) = rxt(k,345)*y(k,69)
         mat(k,179) = -(rxt(k,272)*y(k,77))
         mat(k,877) = -rxt(k,272)*y(k,71)
         mat(k,1680) = rxt(k,208)*y(k,28)
         mat(k,1786) = rxt(k,208)*y(k,5) + 2.000_r8*rxt(k,178)*y(k,28)
         mat(k,226) = -(rxt(k,261)*y(k,68) + rxt(k,262)*y(k,96))
         mat(k,1514) = -rxt(k,261)*y(k,72)
         mat(k,1323) = -rxt(k,262)*y(k,72)
         mat(k,404) = -(rxt(k,263)*y(k,69) + rxt(k,265)*y(k,70) + rxt(k,268)*y(k,96))
         mat(k,1252) = -rxt(k,263)*y(k,75)
         mat(k,1735) = -rxt(k,265)*y(k,75)
         mat(k,1331) = -rxt(k,268)*y(k,75)
         mat(k,880) = -(rxt(k,266)*y(k,5) + rxt(k,267)*y(k,28) + rxt(k,269)*y(k,62) &
                      + rxt(k,270)*y(k,69) + rxt(k,271)*y(k,70) + rxt(k,272)*y(k,71) &
                      + rxt(k,273)*y(k,96))
         mat(k,1684) = -rxt(k,266)*y(k,77)
         mat(k,1791) = -rxt(k,267)*y(k,77)
         mat(k,1401) = -rxt(k,269)*y(k,77)
         mat(k,1288) = -rxt(k,270)*y(k,77)
         mat(k,1757) = -rxt(k,271)*y(k,77)
         mat(k,181) = -rxt(k,272)*y(k,77)
         mat(k,1341) = -rxt(k,273)*y(k,77)
         mat(k,1568) = rxt(k,261)*y(k,72)
         mat(k,1288) = mat(k,1288) + rxt(k,263)*y(k,75)
         mat(k,1757) = mat(k,1757) + rxt(k,265)*y(k,75)
         mat(k,230) = rxt(k,261)*y(k,68)
         mat(k,405) = rxt(k,263)*y(k,69) + rxt(k,265)*y(k,70) + rxt(k,268)*y(k,96)
         mat(k,1341) = mat(k,1341) + rxt(k,268)*y(k,75)
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
         mat(k,574) = -(rxt(k,264)*y(k,96))
         mat(k,1335) = -rxt(k,264)*y(k,78)
         mat(k,1683) = rxt(k,266)*y(k,77)
         mat(k,1788) = rxt(k,267)*y(k,77)
         mat(k,171) = rxt(k,259)*y(k,63) + (rxt(k,260)+.500_r8*rxt(k,274))*y(k,96)
         mat(k,1398) = rxt(k,269)*y(k,77)
         mat(k,1658) = rxt(k,259)*y(k,35)
         mat(k,1267) = rxt(k,270)*y(k,77)
         mat(k,1743) = rxt(k,271)*y(k,77)
         mat(k,180) = rxt(k,272)*y(k,77)
         mat(k,229) = rxt(k,262)*y(k,96)
         mat(k,879) = rxt(k,266)*y(k,5) + rxt(k,267)*y(k,28) + rxt(k,269)*y(k,62) &
                      + rxt(k,270)*y(k,69) + rxt(k,271)*y(k,70) + rxt(k,272)*y(k,71) &
                      + rxt(k,273)*y(k,96)
         mat(k,1335) = mat(k,1335) + (rxt(k,260)+.500_r8*rxt(k,274))*y(k,35) &
                      + rxt(k,262)*y(k,72) + rxt(k,273)*y(k,77)
         mat(k,116) = -(rxt(k,275)*y(k,176))
         mat(k,1815) = -rxt(k,275)*y(k,79)
         mat(k,573) = rxt(k,264)*y(k,96)
         mat(k,1318) = rxt(k,264)*y(k,78)
         mat(k,1042) = -(rxt(k,300)*y(k,90) + rxt(k,301)*y(k,95) + rxt(k,302)*y(k,89) &
                      + rxt(k,346)*y(k,107) + rxt(k,347)*y(k,108) + rxt(k,348) &
                      *y(k,109) + rxt(k,349)*y(k,106) + rxt(k,350)*y(k,105) + rxt(k,372) &
                      *y(k,117) + rxt(k,373)*y(k,118) + rxt(k,374)*y(k,119) + rxt(k,375) &
                      *y(k,120) + rxt(k,404)*y(k,130) + rxt(k,405)*y(k,131) + rxt(k,406) &
                      *y(k,134) + rxt(k,407)*y(k,135) + rxt(k,408)*y(k,133) + rxt(k,426) &
                      *y(k,142) + rxt(k,427)*y(k,143) + rxt(k,428)*y(k,144) + rxt(k,431) &
                      *y(k,141) + rxt(k,471)*y(k,163) + (rxt(k,472) + rxt(k,476) &
                      ) * y(k,164) + rxt(k,473)*y(k,160) + rxt(k,474)*y(k,162) &
                      + rxt(k,475)*y(k,155) + rxt(k,516)*y(k,171) + rxt(k,517) &
                      *y(k,172) + rxt(k,518)*y(k,174) + rxt(k,519)*y(k,175) + rxt(k,520) &
                      *y(k,173) + rxt(k,521)*y(k,170))
         mat(k,1102) = -rxt(k,300)*y(k,86)
         mat(k,1071) = -rxt(k,301)*y(k,86)
         mat(k,344) = -rxt(k,302)*y(k,86)
         mat(k,635) = -rxt(k,346)*y(k,86)
         mat(k,191) = -rxt(k,347)*y(k,86)
         mat(k,122) = -rxt(k,348)*y(k,86)
         mat(k,393) = -rxt(k,349)*y(k,86)
         mat(k,431) = -rxt(k,350)*y(k,86)
         mat(k,455) = -rxt(k,372)*y(k,86)
         mat(k,375) = -rxt(k,373)*y(k,86)
         mat(k,307) = -rxt(k,374)*y(k,86)
         mat(k,200) = -rxt(k,375)*y(k,86)
         mat(k,673) = -rxt(k,404)*y(k,86)
         mat(k,415) = -rxt(k,405)*y(k,86)
         mat(k,211) = -rxt(k,406)*y(k,86)
         mat(k,142) = -rxt(k,407)*y(k,86)
         mat(k,138) = -rxt(k,408)*y(k,86)
         mat(k,971) = -rxt(k,426)*y(k,86)
         mat(k,318) = -rxt(k,427)*y(k,86)
         mat(k,463) = -rxt(k,428)*y(k,86)
         mat(k,852) = -rxt(k,431)*y(k,86)
         mat(k,480) = -rxt(k,471)*y(k,86)
         mat(k,216) = -(rxt(k,472) + rxt(k,476)) * y(k,86)
         mat(k,471) = -rxt(k,473)*y(k,86)
         mat(k,423) = -rxt(k,474)*y(k,86)
         mat(k,686) = -rxt(k,475)*y(k,86)
         mat(k,736) = -rxt(k,516)*y(k,86)
         mat(k,748) = -rxt(k,517)*y(k,86)
         mat(k,489) = -rxt(k,518)*y(k,86)
         mat(k,330) = -rxt(k,519)*y(k,86)
         mat(k,760) = -rxt(k,520)*y(k,86)
         mat(k,772) = -rxt(k,521)*y(k,86)
         mat(k,1380) = -((rxt(k,126) + rxt(k,127) + rxt(k,128)) * y(k,40) + rxt(k,130) &
                      *y(k,68) + rxt(k,131)*y(k,70) + rxt(k,135)*y(k,96) &
                      + 4._r8*rxt(k,140)*y(k,87) + rxt(k,152)*y(k,63) + rxt(k,157) &
                      *y(k,61) + rxt(k,162)*y(k,62) + (rxt(k,172) + rxt(k,173) &
                      ) * y(k,25) + rxt(k,179)*y(k,28) + rxt(k,205)*y(k,3) + rxt(k,211) &
                      *y(k,5) + rxt(k,251)*y(k,21))
         mat(k,1204) = -(rxt(k,126) + rxt(k,127) + rxt(k,128)) * y(k,87)
         mat(k,1583) = -rxt(k,130)*y(k,87)
         mat(k,1770) = -rxt(k,131)*y(k,87)
         mat(k,1353) = -rxt(k,135)*y(k,87)
         mat(k,1666) = -rxt(k,152)*y(k,87)
         mat(k,1715) = -rxt(k,157)*y(k,87)
         mat(k,1412) = -rxt(k,162)*y(k,87)
         mat(k,1494) = -(rxt(k,172) + rxt(k,173)) * y(k,87)
         mat(k,1799) = -rxt(k,179)*y(k,87)
         mat(k,990) = -rxt(k,205)*y(k,87)
         mat(k,1692) = -rxt(k,211)*y(k,87)
         mat(k,823) = -rxt(k,251)*y(k,87)
         mat(k,990) = mat(k,990) + rxt(k,204)*y(k,17)
         mat(k,1692) = mat(k,1692) + rxt(k,216)*y(k,96)
         mat(k,1605) = rxt(k,204)*y(k,3) + rxt(k,168)*y(k,25) + rxt(k,248)*y(k,63) &
                      + rxt(k,249)*y(k,68)
         mat(k,354) = rxt(k,233)*y(k,25) + rxt(k,234)*y(k,96)
         mat(k,223) = rxt(k,236)*y(k,25) + rxt(k,237)*y(k,96)
         mat(k,823) = mat(k,823) + rxt(k,175)*y(k,28) + rxt(k,252)*y(k,61)
         mat(k,796) = rxt(k,256)*y(k,92)
         mat(k,1494) = mat(k,1494) + rxt(k,168)*y(k,17) + rxt(k,233)*y(k,18) &
                      + rxt(k,236)*y(k,20) + rxt(k,171)*y(k,43)
         mat(k,1799) = mat(k,1799) + rxt(k,175)*y(k,21) + rxt(k,186)*y(k,96)
         mat(k,388) = rxt(k,258)*y(k,96)
         mat(k,174) = .500_r8*rxt(k,274)*y(k,96)
         mat(k,1204) = mat(k,1204) + rxt(k,129)*y(k,69) + rxt(k,486)*y(k,167)
         mat(k,361) = rxt(k,171)*y(k,25) + rxt(k,125)*y(k,68) + rxt(k,134)*y(k,96)
         mat(k,1715) = mat(k,1715) + rxt(k,252)*y(k,21)
         mat(k,1666) = mat(k,1666) + rxt(k,248)*y(k,17) + rxt(k,155)*y(k,96)
         mat(k,1583) = mat(k,1583) + rxt(k,249)*y(k,17) + rxt(k,125)*y(k,43)
         mat(k,1301) = rxt(k,129)*y(k,40)
         mat(k,1770) = mat(k,1770) + rxt(k,137)*y(k,96)
         mat(k,577) = rxt(k,264)*y(k,96)
         mat(k,1463) = rxt(k,256)*y(k,23)
         mat(k,1353) = mat(k,1353) + rxt(k,216)*y(k,5) + rxt(k,234)*y(k,18) &
                      + rxt(k,237)*y(k,20) + rxt(k,186)*y(k,28) + rxt(k,258)*y(k,31) &
                      + .500_r8*rxt(k,274)*y(k,35) + rxt(k,134)*y(k,43) + rxt(k,155) &
                      *y(k,63) + rxt(k,137)*y(k,70) + rxt(k,264)*y(k,78)
         mat(k,518) = rxt(k,486)*y(k,40)
         mat(k,433) = -(rxt(k,142)*y(k,68) + rxt(k,143)*y(k,69) + rxt(k,314)*y(k,97))
         mat(k,1532) = -rxt(k,142)*y(k,88)
         mat(k,1256) = -rxt(k,143)*y(k,88)
         mat(k,653) = -rxt(k,314)*y(k,88)
         mat(k,1532) = mat(k,1532) + rxt(k,304)*y(k,89)
         mat(k,1017) = .900_r8*rxt(k,302)*y(k,89) + .800_r8*rxt(k,300)*y(k,90)
         mat(k,341) = rxt(k,304)*y(k,68) + .900_r8*rxt(k,302)*y(k,86)
         mat(k,1087) = .800_r8*rxt(k,300)*y(k,86)
         mat(k,340) = -(rxt(k,302)*y(k,86) + rxt(k,303)*y(k,69) + (rxt(k,304) &
                      + rxt(k,305)) * y(k,68))
         mat(k,1011) = -rxt(k,302)*y(k,89)
         mat(k,1248) = -rxt(k,303)*y(k,89)
         mat(k,1523) = -(rxt(k,304) + rxt(k,305)) * y(k,89)
         mat(k,1104) = -(rxt(k,300)*y(k,86) + rxt(k,337)*y(k,100) + rxt(k,366) &
                      *y(k,110) + rxt(k,393)*y(k,121) + rxt(k,461)*y(k,145) + rxt(k,491) &
                      *y(k,165))
         mat(k,1044) = -rxt(k,300)*y(k,90)
         mat(k,901) = -rxt(k,337)*y(k,90)
         mat(k,809) = -rxt(k,366)*y(k,90)
         mat(k,835) = -rxt(k,393)*y(k,90)
         mat(k,937) = -rxt(k,461)*y(k,90)
         mat(k,917) = -rxt(k,491)*y(k,90)
         mat(k,867) = rxt(k,309)*y(k,95)
         mat(k,1710) = rxt(k,311)*y(k,95)
         mat(k,1577) = rxt(k,304)*y(k,89)
         mat(k,1296) = rxt(k,308)*y(k,91)
         mat(k,346) = rxt(k,304)*y(k,68)
         mat(k,501) = rxt(k,308)*y(k,69)
         mat(k,1073) = rxt(k,309)*y(k,55) + rxt(k,311)*y(k,61)
         mat(k,497) = -(rxt(k,306)*y(k,68) + (rxt(k,307) + rxt(k,308)) * y(k,69))
         mat(k,1538) = -rxt(k,306)*y(k,91)
         mat(k,1262) = -(rxt(k,307) + rxt(k,308)) * y(k,91)
         mat(k,434) = rxt(k,314)*y(k,97)
         mat(k,654) = rxt(k,314)*y(k,88)
         mat(k,1465) = -(rxt(k,107)*y(k,41) + rxt(k,108)*y(k,176) + (rxt(k,110) &
                      + rxt(k,111)) * y(k,69) + rxt(k,112)*y(k,70) + (rxt(k,160) &
                      + rxt(k,161)) * y(k,56) + rxt(k,193)*y(k,8) + rxt(k,194)*y(k,9) &
                      + rxt(k,195)*y(k,11) + rxt(k,196)*y(k,12) + rxt(k,197)*y(k,13) &
                      + rxt(k,198)*y(k,14) + rxt(k,199)*y(k,15) + (rxt(k,200) &
                      + rxt(k,201)) * y(k,49) + rxt(k,220)*y(k,10) + rxt(k,221) &
                      *y(k,24) + rxt(k,222)*y(k,42) + (rxt(k,223) + rxt(k,224) &
                      ) * y(k,45) + rxt(k,229)*y(k,33) + rxt(k,230)*y(k,34) + rxt(k,243) &
                      *y(k,16) + rxt(k,244)*y(k,18) + rxt(k,245)*y(k,46) + rxt(k,246) &
                      *y(k,47) + rxt(k,247)*y(k,48) + (rxt(k,255) + rxt(k,256) &
                      + rxt(k,257)) * y(k,23))
         mat(k,1644) = -rxt(k,107)*y(k,92)
         mat(k,1876) = -rxt(k,108)*y(k,92)
         mat(k,1303) = -(rxt(k,110) + rxt(k,111)) * y(k,92)
         mat(k,1772) = -rxt(k,112)*y(k,92)
         mat(k,113) = -(rxt(k,160) + rxt(k,161)) * y(k,92)
         mat(k,33) = -rxt(k,193)*y(k,92)
         mat(k,67) = -rxt(k,194)*y(k,92)
         mat(k,47) = -rxt(k,195)*y(k,92)
         mat(k,77) = -rxt(k,196)*y(k,92)
         mat(k,51) = -rxt(k,197)*y(k,92)
         mat(k,82) = -rxt(k,198)*y(k,92)
         mat(k,55) = -rxt(k,199)*y(k,92)
         mat(k,957) = -(rxt(k,200) + rxt(k,201)) * y(k,92)
         mat(k,73) = -rxt(k,220)*y(k,92)
         mat(k,250) = -rxt(k,221)*y(k,92)
         mat(k,44) = -rxt(k,222)*y(k,92)
         mat(k,544) = -(rxt(k,223) + rxt(k,224)) * y(k,92)
         mat(k,98) = -rxt(k,229)*y(k,92)
         mat(k,106) = -rxt(k,230)*y(k,92)
         mat(k,295) = -rxt(k,243)*y(k,92)
         mat(k,355) = -rxt(k,244)*y(k,92)
         mat(k,102) = -rxt(k,245)*y(k,92)
         mat(k,94) = -rxt(k,246)*y(k,92)
         mat(k,167) = -rxt(k,247)*y(k,92)
         mat(k,797) = -(rxt(k,255) + rxt(k,256) + rxt(k,257)) * y(k,92)
         mat(k,1303) = mat(k,1303) + rxt(k,143)*y(k,88)
         mat(k,1051) = .850_r8*rxt(k,301)*y(k,95)
         mat(k,438) = rxt(k,143)*y(k,69)
         mat(k,1080) = .850_r8*rxt(k,301)*y(k,86)
         mat(k,556) = -(rxt(k,114)*y(k,68) + rxt(k,115)*y(k,69) + rxt(k,419)*y(k,141) &
                      + rxt(k,433)*y(k,145))
         mat(k,1545) = -rxt(k,114)*y(k,93)
         mat(k,1266) = -rxt(k,115)*y(k,93)
         mat(k,846) = -rxt(k,419)*y(k,93)
         mat(k,928) = -rxt(k,433)*y(k,93)
         mat(k,1134) = rxt(k,116)*y(k,94)
         mat(k,1545) = mat(k,1545) + rxt(k,118)*y(k,94)
         mat(k,1266) = mat(k,1266) + rxt(k,119)*y(k,94)
         mat(k,1742) = rxt(k,120)*y(k,94)
         mat(k,151) = rxt(k,116)*y(k,32) + rxt(k,118)*y(k,68) + rxt(k,119)*y(k,69) &
                      + rxt(k,120)*y(k,70)
         mat(k,150) = -(rxt(k,116)*y(k,32) + rxt(k,118)*y(k,68) + rxt(k,119)*y(k,69) &
                      + rxt(k,120)*y(k,70))
         mat(k,1119) = -rxt(k,116)*y(k,94)
         mat(k,1511) = -rxt(k,118)*y(k,94)
         mat(k,1238) = -rxt(k,119)*y(k,94)
         mat(k,1728) = -rxt(k,120)*y(k,94)
         mat(k,1238) = mat(k,1238) + rxt(k,110)*y(k,92)
         mat(k,1439) = rxt(k,110)*y(k,69)
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
         mat(k,1072) = -(rxt(k,301)*y(k,86) + rxt(k,309)*y(k,55) + rxt(k,311)*y(k,61) &
                      + rxt(k,336)*y(k,100) + rxt(k,365)*y(k,110) + rxt(k,392) &
                      *y(k,121) + rxt(k,429)*y(k,137) + rxt(k,460)*y(k,145) + rxt(k,492) &
                      *y(k,165))
         mat(k,1043) = -rxt(k,301)*y(k,95)
         mat(k,866) = -rxt(k,309)*y(k,95)
         mat(k,1709) = -rxt(k,311)*y(k,95)
         mat(k,900) = -rxt(k,336)*y(k,95)
         mat(k,808) = -rxt(k,365)*y(k,95)
         mat(k,834) = -rxt(k,392)*y(k,95)
         mat(k,780) = -rxt(k,429)*y(k,95)
         mat(k,936) = -rxt(k,460)*y(k,95)
         mat(k,916) = -rxt(k,492)*y(k,95)
         mat(k,1154) = rxt(k,312)*y(k,97)
         mat(k,1295) = rxt(k,303)*y(k,89) + rxt(k,307)*y(k,91) + rxt(k,315)*y(k,97) &
                      + rxt(k,319)*y(k,98)
         mat(k,1764) = rxt(k,417)*y(k,141) + rxt(k,421)*y(k,142)
         mat(k,345) = rxt(k,303)*y(k,69)
         mat(k,500) = rxt(k,307)*y(k,69)
         mat(k,657) = rxt(k,312)*y(k,32) + rxt(k,315)*y(k,69)
         mat(k,273) = rxt(k,319)*y(k,69)
         mat(k,853) = rxt(k,417)*y(k,70)
         mat(k,972) = rxt(k,421)*y(k,70)
         mat(k,1352) = -(rxt(k,133)*y(k,41) + rxt(k,134)*y(k,43) + rxt(k,135)*y(k,87) &
                      + rxt(k,136)*y(k,68) + rxt(k,137)*y(k,70) + (4._r8*rxt(k,138) &
                      + 4._r8*rxt(k,139)) * y(k,96) + rxt(k,141)*y(k,52) + rxt(k,155) &
                      *y(k,63) + rxt(k,156)*y(k,55) + rxt(k,164)*y(k,62) + rxt(k,165) &
                      *y(k,51) + rxt(k,184)*y(k,29) + (rxt(k,186) + rxt(k,187) &
                      ) * y(k,28) + rxt(k,189)*y(k,49) + rxt(k,192)*y(k,54) + rxt(k,216) &
                      *y(k,5) + rxt(k,218)*y(k,45) + rxt(k,232)*y(k,16) + rxt(k,234) &
                      *y(k,18) + rxt(k,235)*y(k,19) + rxt(k,237)*y(k,20) + rxt(k,239) &
                      *y(k,24) + rxt(k,240)*y(k,46) + rxt(k,241)*y(k,47) + rxt(k,242) &
                      *y(k,48) + rxt(k,250)*y(k,17) + rxt(k,253)*y(k,22) + rxt(k,254) &
                      *y(k,23) + rxt(k,258)*y(k,31) + (rxt(k,260) + rxt(k,274) &
                      ) * y(k,35) + rxt(k,262)*y(k,72) + rxt(k,264)*y(k,78) + rxt(k,268) &
                      *y(k,75) + rxt(k,273)*y(k,77) + rxt(k,414)*y(k,137))
         mat(k,1641) = -rxt(k,133)*y(k,96)
         mat(k,360) = -rxt(k,134)*y(k,96)
         mat(k,1379) = -rxt(k,135)*y(k,96)
         mat(k,1582) = -rxt(k,136)*y(k,96)
         mat(k,1769) = -rxt(k,137)*y(k,96)
         mat(k,278) = -rxt(k,141)*y(k,96)
         mat(k,1665) = -rxt(k,155)*y(k,96)
         mat(k,871) = -rxt(k,156)*y(k,96)
         mat(k,1411) = -rxt(k,164)*y(k,96)
         mat(k,1222) = -rxt(k,165)*y(k,96)
         mat(k,585) = -rxt(k,184)*y(k,96)
         mat(k,1798) = -(rxt(k,186) + rxt(k,187)) * y(k,96)
         mat(k,955) = -rxt(k,189)*y(k,96)
         mat(k,551) = -rxt(k,192)*y(k,96)
         mat(k,1691) = -rxt(k,216)*y(k,96)
         mat(k,543) = -rxt(k,218)*y(k,96)
         mat(k,294) = -rxt(k,232)*y(k,96)
         mat(k,353) = -rxt(k,234)*y(k,96)
         mat(k,58) = -rxt(k,235)*y(k,96)
         mat(k,222) = -rxt(k,237)*y(k,96)
         mat(k,249) = -rxt(k,239)*y(k,96)
         mat(k,101) = -rxt(k,240)*y(k,96)
         mat(k,93) = -rxt(k,241)*y(k,96)
         mat(k,166) = -rxt(k,242)*y(k,96)
         mat(k,1604) = -rxt(k,250)*y(k,96)
         mat(k,237) = -rxt(k,253)*y(k,96)
         mat(k,795) = -rxt(k,254)*y(k,96)
         mat(k,387) = -rxt(k,258)*y(k,96)
         mat(k,173) = -(rxt(k,260) + rxt(k,274)) * y(k,96)
         mat(k,232) = -rxt(k,262)*y(k,96)
         mat(k,576) = -rxt(k,264)*y(k,96)
         mat(k,408) = -rxt(k,268)*y(k,96)
         mat(k,884) = -rxt(k,273)*y(k,96)
         mat(k,783) = -rxt(k,414)*y(k,96)
         mat(k,1604) = mat(k,1604) + rxt(k,249)*y(k,68)
         mat(k,237) = mat(k,237) + .300_r8*rxt(k,253)*y(k,96)
         mat(k,795) = mat(k,795) + rxt(k,255)*y(k,92)
         mat(k,1493) = rxt(k,173)*y(k,87)
         mat(k,569) = rxt(k,227)*y(k,176)
         mat(k,1203) = rxt(k,132)*y(k,70) + 2.000_r8*rxt(k,127)*y(k,87)
         mat(k,1641) = mat(k,1641) + rxt(k,124)*y(k,68) + rxt(k,107)*y(k,92)
         mat(k,360) = mat(k,360) + rxt(k,125)*y(k,68)
         mat(k,543) = mat(k,543) + rxt(k,217)*y(k,68) + rxt(k,223)*y(k,92)
         mat(k,955) = mat(k,955) + rxt(k,188)*y(k,68) + rxt(k,200)*y(k,92)
         mat(k,445) = rxt(k,219)*y(k,68)
         mat(k,551) = mat(k,551) + rxt(k,191)*y(k,68)
         mat(k,1714) = rxt(k,157)*y(k,87)
         mat(k,1665) = mat(k,1665) + rxt(k,152)*y(k,87)
         mat(k,1582) = mat(k,1582) + rxt(k,249)*y(k,17) + rxt(k,124)*y(k,41) &
                      + rxt(k,125)*y(k,43) + rxt(k,217)*y(k,45) + rxt(k,188)*y(k,49) &
                      + rxt(k,219)*y(k,53) + rxt(k,191)*y(k,54) + rxt(k,130)*y(k,87)
         mat(k,1769) = mat(k,1769) + rxt(k,132)*y(k,40) + rxt(k,131)*y(k,87)
         mat(k,1048) = rxt(k,427)*y(k,143)
         mat(k,1379) = mat(k,1379) + rxt(k,173)*y(k,25) + 2.000_r8*rxt(k,127)*y(k,40) &
                      + rxt(k,157)*y(k,61) + rxt(k,152)*y(k,63) + rxt(k,130)*y(k,68) &
                      + rxt(k,131)*y(k,70)
         mat(k,1462) = rxt(k,255)*y(k,23) + rxt(k,107)*y(k,41) + rxt(k,223)*y(k,45) &
                      + rxt(k,200)*y(k,49) + 2.000_r8*rxt(k,108)*y(k,176)
         mat(k,1352) = mat(k,1352) + .300_r8*rxt(k,253)*y(k,22)
         mat(k,598) = rxt(k,329)*y(k,176)
         mat(k,977) = rxt(k,422)*y(k,176)
         mat(k,319) = rxt(k,427)*y(k,86)
         mat(k,723) = rxt(k,484)*y(k,176)
         mat(k,1873) = rxt(k,227)*y(k,39) + 2.000_r8*rxt(k,108)*y(k,92) + rxt(k,329) &
                      *y(k,101) + rxt(k,422)*y(k,142) + rxt(k,484)*y(k,166)
         mat(k,655) = -(rxt(k,312)*y(k,32) + rxt(k,314)*y(k,88) + rxt(k,315)*y(k,69))
         mat(k,1137) = -rxt(k,312)*y(k,97)
         mat(k,435) = -rxt(k,314)*y(k,97)
         mat(k,1273) = -rxt(k,315)*y(k,97)
         mat(k,1553) = rxt(k,305)*y(k,89) + rxt(k,306)*y(k,91) + rxt(k,318)*y(k,98) &
                      + rxt(k,324)*y(k,99)
         mat(k,1026) = rxt(k,316)*y(k,98) + rxt(k,321)*y(k,99)
         mat(k,342) = rxt(k,305)*y(k,68)
         mat(k,498) = rxt(k,306)*y(k,68)
         mat(k,272) = rxt(k,318)*y(k,68) + rxt(k,316)*y(k,86)
         mat(k,267) = rxt(k,324)*y(k,68) + rxt(k,321)*y(k,86)
         mat(k,270) = -(rxt(k,316)*y(k,86) + rxt(k,318)*y(k,68) + rxt(k,319)*y(k,69))
         mat(k,1005) = -rxt(k,316)*y(k,98)
         mat(k,1516) = -rxt(k,318)*y(k,98)
         mat(k,1243) = -rxt(k,319)*y(k,98)
         mat(k,1005) = mat(k,1005) + rxt(k,320)*y(k,99)
         mat(k,264) = rxt(k,320)*y(k,86)
         mat(k,263) = -((rxt(k,320) + rxt(k,321)) * y(k,86) + rxt(k,324)*y(k,68))
         mat(k,1004) = -(rxt(k,320) + rxt(k,321)) * y(k,99)
         mat(k,1515) = -rxt(k,324)*y(k,99)
         mat(k,898) = -(rxt(k,325)*y(k,70) + rxt(k,333)*y(k,69) + rxt(k,336)*y(k,95) &
                      + rxt(k,337)*y(k,90))
         mat(k,1758) = -rxt(k,325)*y(k,100)
         mat(k,1289) = -rxt(k,333)*y(k,100)
         mat(k,1067) = -rxt(k,336)*y(k,100)
         mat(k,1098) = -rxt(k,337)*y(k,100)
         mat(k,1193) = rxt(k,330)*y(k,103) + rxt(k,332)*y(k,104)
         mat(k,1569) = rxt(k,326)*y(k,101)
         mat(k,1758) = mat(k,1758) + rxt(k,327)*y(k,101)
         mat(k,1038) = rxt(k,350)*y(k,105) + rxt(k,349)*y(k,106) + rxt(k,346)*y(k,107) &
                      + rxt(k,347)*y(k,108) + rxt(k,348)*y(k,109)
         mat(k,594) = rxt(k,326)*y(k,68) + rxt(k,327)*y(k,70)
         mat(k,367) = rxt(k,330)*y(k,40)
         mat(k,241) = rxt(k,332)*y(k,40)
         mat(k,430) = rxt(k,350)*y(k,86)
         mat(k,392) = rxt(k,349)*y(k,86)
         mat(k,634) = rxt(k,346)*y(k,86)
         mat(k,190) = rxt(k,347)*y(k,86)
         mat(k,121) = rxt(k,348)*y(k,86)
         mat(k,593) = -(rxt(k,326)*y(k,68) + rxt(k,327)*y(k,70) + rxt(k,328)*y(k,41) &
                      + rxt(k,329)*y(k,176))
         mat(k,1548) = -rxt(k,326)*y(k,101)
         mat(k,1744) = -rxt(k,327)*y(k,101)
         mat(k,1624) = -rxt(k,328)*y(k,101)
         mat(k,1839) = -rxt(k,329)*y(k,101)
         mat(k,1548) = mat(k,1548) + rxt(k,334)*y(k,102)
         mat(k,1744) = mat(k,1744) + rxt(k,325)*y(k,100)
         mat(k,896) = rxt(k,325)*y(k,70)
         mat(k,185) = rxt(k,334)*y(k,68)
         mat(k,184) = -(rxt(k,334)*y(k,68))
         mat(k,1512) = -rxt(k,334)*y(k,102)
         mat(k,1239) = rxt(k,333)*y(k,100)
         mat(k,894) = rxt(k,333)*y(k,69)
         mat(k,366) = -(rxt(k,330)*y(k,40) + rxt(k,331)*y(k,32))
         mat(k,1176) = -rxt(k,330)*y(k,103)
         mat(k,1124) = -rxt(k,331)*y(k,103)
         mat(k,1620) = rxt(k,328)*y(k,101)
         mat(k,592) = rxt(k,328)*y(k,41) + rxt(k,329)*y(k,176)
         mat(k,1826) = rxt(k,329)*y(k,101)
         mat(k,240) = -(rxt(k,332)*y(k,40) + 4._r8*rxt(k,335)*y(k,104))
         mat(k,1171) = -rxt(k,332)*y(k,104)
         mat(k,1121) = rxt(k,331)*y(k,103)
         mat(k,365) = rxt(k,331)*y(k,32)
         mat(k,428) = -(rxt(k,339)*y(k,32) + rxt(k,350)*y(k,86))
         mat(k,1128) = -rxt(k,339)*y(k,105)
         mat(k,1016) = -rxt(k,350)*y(k,105)
         mat(k,1531) = rxt(k,343)*y(k,106)
         mat(k,1255) = rxt(k,345)*y(k,106)
         mat(k,1086) = rxt(k,337)*y(k,100)
         mat(k,1057) = rxt(k,336)*y(k,100)
         mat(k,895) = rxt(k,337)*y(k,90) + rxt(k,336)*y(k,95)
         mat(k,390) = rxt(k,343)*y(k,68) + rxt(k,345)*y(k,69)
         mat(k,389) = -(rxt(k,343)*y(k,68) + rxt(k,345)*y(k,69) + rxt(k,349)*y(k,86))
         mat(k,1527) = -rxt(k,343)*y(k,106)
         mat(k,1250) = -rxt(k,345)*y(k,106)
         mat(k,1013) = -rxt(k,349)*y(k,106)
         mat(k,1527) = mat(k,1527) + rxt(k,342)*y(k,107)
         mat(k,631) = rxt(k,342)*y(k,68)
         mat(k,633) = -(rxt(k,340)*y(k,32) + rxt(k,341)*y(k,176) + rxt(k,342)*y(k,68) &
                      + rxt(k,346)*y(k,86))
         mat(k,1135) = -rxt(k,340)*y(k,107)
         mat(k,1842) = -rxt(k,341)*y(k,107)
         mat(k,1551) = -rxt(k,342)*y(k,107)
         mat(k,1025) = -rxt(k,346)*y(k,107)
         mat(k,189) = -(rxt(k,347)*y(k,86))
         mat(k,1000) = -rxt(k,347)*y(k,108)
         mat(k,1120) = rxt(k,339)*y(k,105) + rxt(k,340)*y(k,107)
         mat(k,427) = rxt(k,339)*y(k,32)
         mat(k,630) = rxt(k,340)*y(k,32)
         mat(k,120) = -(rxt(k,348)*y(k,86))
         mat(k,997) = -rxt(k,348)*y(k,109)
         mat(k,629) = rxt(k,341)*y(k,176)
         mat(k,1816) = rxt(k,341)*y(k,107)
         mat(k,806) = -(rxt(k,351)*y(k,70) + rxt(k,365)*y(k,95) + rxt(k,366)*y(k,90))
         mat(k,1754) = -rxt(k,351)*y(k,110)
         mat(k,1064) = -rxt(k,365)*y(k,110)
         mat(k,1095) = -rxt(k,366)*y(k,110)
         mat(k,1191) = rxt(k,360)*y(k,115)
         mat(k,1564) = rxt(k,352)*y(k,111)
         mat(k,1034) = rxt(k,372)*y(k,117) + rxt(k,373)*y(k,118) + rxt(k,374)*y(k,119) &
                      + rxt(k,375)*y(k,120)
         mat(k,621) = rxt(k,352)*y(k,68)
         mat(k,335) = rxt(k,360)*y(k,40)
         mat(k,454) = rxt(k,372)*y(k,86)
         mat(k,374) = rxt(k,373)*y(k,86)
         mat(k,306) = rxt(k,374)*y(k,86)
         mat(k,199) = rxt(k,375)*y(k,86)
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
         mat(k,620) = -(rxt(k,352)*y(k,68) + rxt(k,353)*y(k,70) + rxt(k,368)*y(k,69) &
                      + rxt(k,369)*y(k,176))
         mat(k,1550) = -rxt(k,352)*y(k,111)
         mat(k,1746) = -rxt(k,353)*y(k,111)
         mat(k,1270) = -rxt(k,368)*y(k,111)
         mat(k,1841) = -rxt(k,369)*y(k,111)
         mat(k,1185) = rxt(k,376)*y(k,115)
         mat(k,1550) = mat(k,1550) + rxt(k,354)*y(k,112)
         mat(k,1746) = mat(k,1746) + rxt(k,351)*y(k,110)
         mat(k,805) = rxt(k,351)*y(k,70)
         mat(k,300) = rxt(k,354)*y(k,68)
         mat(k,334) = rxt(k,376)*y(k,40)
         mat(k,298) = -(rxt(k,354)*y(k,68) + rxt(k,355)*y(k,70))
         mat(k,1519) = -rxt(k,354)*y(k,112)
         mat(k,1731) = -rxt(k,355)*y(k,112)
         mat(k,1519) = mat(k,1519) + rxt(k,356)*y(k,113)
         mat(k,1731) = mat(k,1731) + rxt(k,353)*y(k,111)
         mat(k,617) = rxt(k,353)*y(k,70)
         mat(k,604) = rxt(k,356)*y(k,68)
         mat(k,606) = -(rxt(k,356)*y(k,68) + rxt(k,357)*y(k,176) + rxt(k,359)*y(k,40))
         mat(k,1549) = -rxt(k,356)*y(k,113)
         mat(k,1840) = -rxt(k,357)*y(k,113)
         mat(k,1184) = -rxt(k,359)*y(k,113)
         mat(k,1269) = rxt(k,368)*y(k,111)
         mat(k,1745) = rxt(k,355)*y(k,112)
         mat(k,619) = rxt(k,368)*y(k,69)
         mat(k,299) = rxt(k,355)*y(k,70)
         mat(k,193) = -(rxt(k,358)*y(k,40))
         mat(k,1169) = -rxt(k,358)*y(k,114)
         mat(k,616) = rxt(k,369)*y(k,176)
         mat(k,603) = rxt(k,357)*y(k,176)
         mat(k,1820) = rxt(k,369)*y(k,111) + rxt(k,357)*y(k,113)
         mat(k,333) = -((rxt(k,360) + rxt(k,376)) * y(k,40) + 4._r8*rxt(k,367) &
                      *y(k,115))
         mat(k,1175) = -(rxt(k,360) + rxt(k,376)) * y(k,115)
         mat(k,1175) = mat(k,1175) + rxt(k,359)*y(k,113) + rxt(k,358)*y(k,114)
         mat(k,605) = rxt(k,359)*y(k,40)
         mat(k,194) = rxt(k,358)*y(k,40)
         mat(k,332) = 2.000_r8*rxt(k,367)*y(k,115)
         mat(k,453) = -(rxt(k,361)*y(k,70) + rxt(k,371)*y(k,69) + rxt(k,372)*y(k,86))
         mat(k,1737) = -rxt(k,361)*y(k,117)
         mat(k,1257) = -rxt(k,371)*y(k,117)
         mat(k,1018) = -rxt(k,372)*y(k,117)
         mat(k,1534) = rxt(k,362)*y(k,118)
         mat(k,1088) = rxt(k,366)*y(k,110)
         mat(k,1058) = rxt(k,365)*y(k,110)
         mat(k,804) = rxt(k,366)*y(k,90) + rxt(k,365)*y(k,95)
         mat(k,373) = rxt(k,362)*y(k,68)
         mat(k,372) = -(rxt(k,362)*y(k,68) + rxt(k,373)*y(k,86))
         mat(k,1525) = -rxt(k,362)*y(k,118)
         mat(k,1012) = -rxt(k,373)*y(k,118)
         mat(k,1525) = mat(k,1525) + rxt(k,364)*y(k,119) + rxt(k,363)*y(k,120)
         mat(k,1733) = rxt(k,361)*y(k,117)
         mat(k,452) = rxt(k,361)*y(k,70)
         mat(k,305) = rxt(k,364)*y(k,68)
         mat(k,198) = rxt(k,363)*y(k,68)
         mat(k,304) = -(rxt(k,364)*y(k,68) + rxt(k,374)*y(k,86))
         mat(k,1520) = -rxt(k,364)*y(k,119)
         mat(k,1008) = -rxt(k,374)*y(k,119)
         mat(k,1246) = rxt(k,371)*y(k,117)
         mat(k,451) = rxt(k,371)*y(k,69)
         mat(k,197) = -(rxt(k,363)*y(k,68) + rxt(k,375)*y(k,86))
         mat(k,1513) = -rxt(k,363)*y(k,120)
         mat(k,1001) = -rxt(k,375)*y(k,120)
         mat(k,832) = -(rxt(k,377)*y(k,70) + rxt(k,392)*y(k,95) + rxt(k,393)*y(k,90))
         mat(k,1755) = -rxt(k,377)*y(k,121)
         mat(k,1065) = -rxt(k,392)*y(k,121)
         mat(k,1096) = -rxt(k,393)*y(k,121)
         mat(k,1192) = rxt(k,389)*y(k,129)
         mat(k,1565) = rxt(k,378)*y(k,122)
         mat(k,1035) = rxt(k,404)*y(k,130) + rxt(k,405)*y(k,131) + rxt(k,408)*y(k,133) &
                      + rxt(k,406)*y(k,134) + rxt(k,407)*y(k,135)
         mat(k,533) = rxt(k,378)*y(k,68)
         mat(k,133) = rxt(k,389)*y(k,40)
         mat(k,672) = rxt(k,404)*y(k,86)
         mat(k,414) = rxt(k,405)*y(k,86)
         mat(k,137) = rxt(k,408)*y(k,86)
         mat(k,210) = rxt(k,406)*y(k,86)
         mat(k,141) = rxt(k,407)*y(k,86)
         mat(k,532) = -(rxt(k,378)*y(k,68) + rxt(k,379)*y(k,70) + rxt(k,381)*y(k,176) &
                      + rxt(k,383)*y(k,69) + rxt(k,384)*y(k,32))
         mat(k,1542) = -rxt(k,378)*y(k,122)
         mat(k,1741) = -rxt(k,379)*y(k,122)
         mat(k,1836) = -rxt(k,381)*y(k,122)
         mat(k,1265) = -rxt(k,383)*y(k,122)
         mat(k,1133) = -rxt(k,384)*y(k,122)
         mat(k,1542) = mat(k,1542) + rxt(k,380)*y(k,123)
         mat(k,1741) = mat(k,1741) + rxt(k,377)*y(k,121)
         mat(k,830) = rxt(k,377)*y(k,70)
         mat(k,313) = rxt(k,380)*y(k,68)
         mat(k,311) = -(rxt(k,380)*y(k,68) + rxt(k,386)*y(k,69))
         mat(k,1521) = -rxt(k,380)*y(k,123)
         mat(k,1247) = -rxt(k,386)*y(k,123)
         mat(k,1521) = mat(k,1521) + rxt(k,385)*y(k,127)
         mat(k,1732) = rxt(k,379)*y(k,122)
         mat(k,530) = rxt(k,379)*y(k,70)
         mat(k,125) = rxt(k,385)*y(k,68)
         mat(k,204) = -(rxt(k,388)*y(k,40) + 4._r8*rxt(k,391)*y(k,124))
         mat(k,1170) = -rxt(k,388)*y(k,124)
         mat(k,529) = rxt(k,381)*y(k,176)
         mat(k,398) = rxt(k,382)*y(k,176)
         mat(k,1821) = rxt(k,381)*y(k,122) + rxt(k,382)*y(k,126)
         mat(k,202) = 2.000_r8*rxt(k,391)*y(k,124)
         mat(k,399) = -(rxt(k,382)*y(k,176) + rxt(k,390)*y(k,40))
         mat(k,1828) = -rxt(k,382)*y(k,126)
         mat(k,1178) = -rxt(k,390)*y(k,126)
         mat(k,1528) = rxt(k,387)*y(k,128)
         mat(k,1251) = rxt(k,383)*y(k,122)
         mat(k,531) = rxt(k,383)*y(k,69)
         mat(k,129) = rxt(k,387)*y(k,68)
         mat(k,124) = -(rxt(k,385)*y(k,68))
         mat(k,1508) = -rxt(k,385)*y(k,127)
         mat(k,1117) = rxt(k,384)*y(k,122)
         mat(k,528) = rxt(k,384)*y(k,32)
         mat(k,128) = -(rxt(k,387)*y(k,68))
         mat(k,1509) = -rxt(k,387)*y(k,128)
         mat(k,1236) = rxt(k,386)*y(k,123)
         mat(k,310) = rxt(k,386)*y(k,69)
         mat(k,132) = -(rxt(k,389)*y(k,40))
         mat(k,1168) = -rxt(k,389)*y(k,129)
         mat(k,1168) = mat(k,1168) + rxt(k,388)*y(k,124) + rxt(k,390)*y(k,126)
         mat(k,203) = rxt(k,388)*y(k,40)
         mat(k,397) = rxt(k,390)*y(k,40)
         mat(k,671) = -(rxt(k,394)*y(k,70) + rxt(k,398)*y(k,69) + rxt(k,402)*y(k,176) &
                      + rxt(k,403)*y(k,32) + rxt(k,404)*y(k,86))
         mat(k,1748) = -rxt(k,394)*y(k,130)
         mat(k,1274) = -rxt(k,398)*y(k,130)
         mat(k,1844) = -rxt(k,402)*y(k,130)
         mat(k,1138) = -rxt(k,403)*y(k,130)
         mat(k,1027) = -rxt(k,404)*y(k,130)
         mat(k,1554) = rxt(k,395)*y(k,131)
         mat(k,1748) = mat(k,1748) + rxt(k,400)*y(k,131)
         mat(k,1092) = rxt(k,393)*y(k,121)
         mat(k,1060) = rxt(k,392)*y(k,121)
         mat(k,831) = rxt(k,393)*y(k,90) + rxt(k,392)*y(k,95)
         mat(k,413) = rxt(k,395)*y(k,68) + rxt(k,400)*y(k,70)
         mat(k,412) = -(rxt(k,395)*y(k,68) + (rxt(k,400) + rxt(k,401)) * y(k,70) &
                      + rxt(k,405)*y(k,86))
         mat(k,1530) = -rxt(k,395)*y(k,131)
         mat(k,1736) = -(rxt(k,400) + rxt(k,401)) * y(k,131)
         mat(k,1014) = -rxt(k,405)*y(k,131)
         mat(k,1530) = mat(k,1530) + rxt(k,399)*y(k,132)
         mat(k,1736) = mat(k,1736) + rxt(k,394)*y(k,130)
         mat(k,670) = rxt(k,394)*y(k,70)
         mat(k,284) = rxt(k,399)*y(k,68)
         mat(k,283) = -(rxt(k,399)*y(k,68))
         mat(k,1517) = -rxt(k,399)*y(k,132)
         mat(k,1244) = rxt(k,398)*y(k,130) + rxt(k,397)*y(k,134)
         mat(k,1730) = rxt(k,401)*y(k,131)
         mat(k,669) = rxt(k,398)*y(k,69)
         mat(k,411) = rxt(k,401)*y(k,70)
         mat(k,209) = rxt(k,397)*y(k,69)
         mat(k,136) = -(rxt(k,408)*y(k,86))
         mat(k,998) = -rxt(k,408)*y(k,133)
         mat(k,666) = rxt(k,402)*y(k,176)
         mat(k,1817) = rxt(k,402)*y(k,130)
         mat(k,208) = -(rxt(k,397)*y(k,69) + rxt(k,406)*y(k,86))
         mat(k,1240) = -rxt(k,397)*y(k,134)
         mat(k,1002) = -rxt(k,406)*y(k,134)
         mat(k,140) = -(rxt(k,407)*y(k,86))
         mat(k,999) = -rxt(k,407)*y(k,135)
         mat(k,1118) = rxt(k,403)*y(k,130)
         mat(k,667) = rxt(k,403)*y(k,32)
         mat(k,492) = -(rxt(k,409)*y(k,69) + rxt(k,410)*y(k,70))
         mat(k,1261) = -rxt(k,409)*y(k,136)
         mat(k,1739) = -rxt(k,410)*y(k,136)
         mat(k,1023) = rxt(k,431)*y(k,141) + rxt(k,426)*y(k,142) + rxt(k,427)*y(k,143) &
                      + rxt(k,428)*y(k,144)
         mat(k,845) = rxt(k,431)*y(k,86)
         mat(k,966) = rxt(k,426)*y(k,86)
         mat(k,317) = rxt(k,427)*y(k,86)
         mat(k,461) = rxt(k,428)*y(k,86)
         mat(k,778) = -(rxt(k,411)*y(k,70) + rxt(k,414)*y(k,96) + rxt(k,429)*y(k,95))
         mat(k,1753) = -rxt(k,411)*y(k,137)
         mat(k,1337) = -rxt(k,414)*y(k,137)
         mat(k,1063) = -rxt(k,429)*y(k,137)
         mat(k,1283) = rxt(k,409)*y(k,136)
         mat(k,1753) = mat(k,1753) + rxt(k,410)*y(k,136) + rxt(k,417)*y(k,141)
         mat(k,1337) = mat(k,1337) + rxt(k,415)*y(k,138)
         mat(k,493) = rxt(k,409)*y(k,69) + rxt(k,410)*y(k,70)
         mat(k,161) = rxt(k,415)*y(k,96)
         mat(k,848) = rxt(k,417)*y(k,70)
         mat(k,160) = -(rxt(k,412)*y(k,176) + rxt(k,415)*y(k,96))
         mat(k,1818) = -rxt(k,412)*y(k,138)
         mat(k,1319) = -rxt(k,415)*y(k,138)
         mat(k,1729) = rxt(k,411)*y(k,137) + rxt(k,421)*y(k,142)
         mat(k,1319) = mat(k,1319) + rxt(k,414)*y(k,137)
         mat(k,777) = rxt(k,411)*y(k,70) + rxt(k,414)*y(k,96)
         mat(k,964) = rxt(k,421)*y(k,70)
         mat(k,36) = -(rxt(k,413)*y(k,176))
         mat(k,1813) = -rxt(k,413)*y(k,139)
         mat(k,159) = rxt(k,412)*y(k,176)
         mat(k,1813) = mat(k,1813) + rxt(k,412)*y(k,138)
         mat(k,35) = rxt(k,413)*y(k,176)
         mat(k,1812) = rxt(k,413)*y(k,139)
         mat(k,849) = -((rxt(k,416) + rxt(k,417)) * y(k,70) + rxt(k,418)*y(k,176) &
                      + rxt(k,419)*y(k,93) + rxt(k,420)*y(k,69) + rxt(k,431)*y(k,86))
         mat(k,1756) = -(rxt(k,416) + rxt(k,417)) * y(k,141)
         mat(k,1858) = -rxt(k,418)*y(k,141)
         mat(k,558) = -rxt(k,419)*y(k,141)
         mat(k,1286) = -rxt(k,420)*y(k,141)
         mat(k,1036) = -rxt(k,431)*y(k,141)
         mat(k,1566) = rxt(k,430)*y(k,142)
         mat(k,968) = rxt(k,430)*y(k,68)
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
         mat(k,970) = -(rxt(k,421)*y(k,70) + rxt(k,422)*y(k,176) + rxt(k,423)*y(k,41) &
                      + rxt(k,426)*y(k,86) + rxt(k,430)*y(k,68))
         mat(k,1761) = -rxt(k,421)*y(k,142)
         mat(k,1864) = -rxt(k,422)*y(k,142)
         mat(k,1632) = -rxt(k,423)*y(k,142)
         mat(k,1041) = -rxt(k,426)*y(k,142)
         mat(k,1573) = -rxt(k,430)*y(k,142)
         mat(k,1573) = mat(k,1573) + rxt(k,424)*y(k,144)
         mat(k,1761) = mat(k,1761) + rxt(k,416)*y(k,141)
         mat(k,560) = rxt(k,419)*y(k,141)
         mat(k,1070) = rxt(k,429)*y(k,137)
         mat(k,779) = rxt(k,429)*y(k,95)
         mat(k,851) = rxt(k,416)*y(k,70) + rxt(k,419)*y(k,93)
         mat(k,462) = rxt(k,424)*y(k,68)
         mat(k,316) = -(rxt(k,427)*y(k,86))
         mat(k,1009) = -rxt(k,427)*y(k,143)
         mat(k,1174) = rxt(k,425)*y(k,144)
         mat(k,1619) = rxt(k,423)*y(k,142)
         mat(k,843) = rxt(k,418)*y(k,176)
         mat(k,965) = rxt(k,423)*y(k,41) + rxt(k,422)*y(k,176)
         mat(k,459) = rxt(k,425)*y(k,40)
         mat(k,1823) = rxt(k,418)*y(k,141) + rxt(k,422)*y(k,142)
         mat(k,460) = -(rxt(k,424)*y(k,68) + rxt(k,425)*y(k,40) + rxt(k,428)*y(k,86))
         mat(k,1535) = -rxt(k,424)*y(k,144)
         mat(k,1179) = -rxt(k,425)*y(k,144)
         mat(k,1019) = -rxt(k,428)*y(k,144)
         mat(k,1258) = rxt(k,420)*y(k,141)
         mat(k,844) = rxt(k,420)*y(k,69)
         mat(k,933) = -(rxt(k,432)*y(k,70) + rxt(k,433)*y(k,93) + rxt(k,436)*y(k,69) &
                      + rxt(k,460)*y(k,95) + rxt(k,461)*y(k,90))
         mat(k,1760) = -rxt(k,432)*y(k,145)
         mat(k,559) = -rxt(k,433)*y(k,145)
         mat(k,1291) = -rxt(k,436)*y(k,145)
         mat(k,1069) = -rxt(k,460)*y(k,145)
         mat(k,1100) = -rxt(k,461)*y(k,145)
         mat(k,1195) = rxt(k,451)*y(k,150)
         mat(k,1571) = rxt(k,434)*y(k,146)
         mat(k,1040) = rxt(k,475)*y(k,155) + rxt(k,473)*y(k,160) + rxt(k,474)*y(k,162) &
                      + rxt(k,471)*y(k,163) + (rxt(k,472)+rxt(k,476))*y(k,164)
         mat(k,709) = rxt(k,434)*y(k,68)
         mat(k,523) = rxt(k,451)*y(k,40)
         mat(k,685) = rxt(k,475)*y(k,86)
         mat(k,470) = rxt(k,473)*y(k,86)
         mat(k,422) = rxt(k,474)*y(k,86)
         mat(k,479) = rxt(k,471)*y(k,86)
         mat(k,215) = (rxt(k,472)+rxt(k,476))*y(k,86)
         mat(k,708) = -(rxt(k,434)*y(k,68) + rxt(k,435)*y(k,70) + rxt(k,438)*y(k,176) &
                      + rxt(k,441)*y(k,69) + rxt(k,442)*y(k,32))
         mat(k,1557) = -rxt(k,434)*y(k,146)
         mat(k,1751) = -rxt(k,435)*y(k,146)
         mat(k,1847) = -rxt(k,438)*y(k,146)
         mat(k,1277) = -rxt(k,441)*y(k,146)
         mat(k,1141) = -rxt(k,442)*y(k,146)
         mat(k,1557) = mat(k,1557) + rxt(k,437)*y(k,147)
         mat(k,1751) = mat(k,1751) + rxt(k,432)*y(k,145)
         mat(k,557) = rxt(k,433)*y(k,145)
         mat(k,931) = rxt(k,432)*y(k,70) + rxt(k,433)*y(k,93)
         mat(k,695) = rxt(k,437)*y(k,68)
         mat(k,694) = -(rxt(k,437)*y(k,68) + rxt(k,447)*y(k,70) + rxt(k,450)*y(k,40))
         mat(k,1556) = -rxt(k,437)*y(k,147)
         mat(k,1750) = -rxt(k,447)*y(k,147)
         mat(k,1187) = -rxt(k,450)*y(k,147)
         mat(k,1556) = mat(k,1556) + rxt(k,445)*y(k,148) + rxt(k,443)*y(k,151)
         mat(k,1276) = rxt(k,436)*y(k,145)
         mat(k,1750) = mat(k,1750) + rxt(k,435)*y(k,146)
         mat(k,930) = rxt(k,436)*y(k,69)
         mat(k,707) = rxt(k,435)*y(k,70)
         mat(k,645) = rxt(k,445)*y(k,68)
         mat(k,508) = rxt(k,443)*y(k,68)
         mat(k,644) = -(rxt(k,439)*y(k,176) + rxt(k,440)*y(k,32) + rxt(k,445)*y(k,68) &
                      + rxt(k,446)*y(k,40))
         mat(k,1843) = -rxt(k,439)*y(k,148)
         mat(k,1136) = -rxt(k,440)*y(k,148)
         mat(k,1552) = -rxt(k,445)*y(k,148)
         mat(k,1186) = -rxt(k,446)*y(k,148)
         mat(k,1272) = rxt(k,441)*y(k,146)
         mat(k,1747) = rxt(k,447)*y(k,147)
         mat(k,706) = rxt(k,441)*y(k,69)
         mat(k,693) = rxt(k,447)*y(k,70)
         mat(k,253) = -(rxt(k,449)*y(k,40) + 4._r8*rxt(k,455)*y(k,149))
         mat(k,1172) = -rxt(k,449)*y(k,149)
         mat(k,703) = rxt(k,438)*y(k,176)
         mat(k,641) = rxt(k,439)*y(k,176)
         mat(k,505) = rxt(k,448)*y(k,176)
         mat(k,1822) = rxt(k,438)*y(k,146) + rxt(k,439)*y(k,148) + rxt(k,448)*y(k,151)
         mat(k,522) = -(rxt(k,451)*y(k,40) + rxt(k,452)*y(k,69) + 4._r8*rxt(k,458) &
                      *y(k,150))
         mat(k,1183) = -rxt(k,451)*y(k,150)
         mat(k,1264) = -rxt(k,452)*y(k,150)
         mat(k,1183) = mat(k,1183) + rxt(k,450)*y(k,147) + rxt(k,446)*y(k,148) &
                      + rxt(k,449)*y(k,149) + rxt(k,444)*y(k,151)
         mat(k,1541) = rxt(k,454)*y(k,152)
         mat(k,692) = rxt(k,450)*y(k,40)
         mat(k,643) = rxt(k,446)*y(k,40)
         mat(k,254) = rxt(k,449)*y(k,40)
         mat(k,507) = rxt(k,444)*y(k,40)
         mat(k,88) = rxt(k,454)*y(k,68)
         mat(k,506) = -(rxt(k,443)*y(k,68) + rxt(k,444)*y(k,40) + rxt(k,448)*y(k,176) &
                      + 4._r8*rxt(k,459)*y(k,151))
         mat(k,1539) = -rxt(k,443)*y(k,151)
         mat(k,1181) = -rxt(k,444)*y(k,151)
         mat(k,1834) = -rxt(k,448)*y(k,151)
         mat(k,1131) = rxt(k,442)*y(k,146) + rxt(k,440)*y(k,148)
         mat(k,704) = rxt(k,442)*y(k,32)
         mat(k,642) = rxt(k,440)*y(k,32)
         mat(k,87) = -(rxt(k,454)*y(k,68) + 4._r8*rxt(k,457)*y(k,152))
         mat(k,1507) = -rxt(k,454)*y(k,152)
         mat(k,1507) = mat(k,1507) + rxt(k,453)*y(k,153)
         mat(k,145) = rxt(k,453)*y(k,68)
         mat(k,146) = -(rxt(k,453)*y(k,68) + 4._r8*rxt(k,456)*y(k,153))
         mat(k,1510) = -rxt(k,453)*y(k,153)
         mat(k,1237) = rxt(k,452)*y(k,150)
         mat(k,521) = rxt(k,452)*y(k,69)
         mat(k,252) = 2.000_r8*rxt(k,455)*y(k,149)
         mat(k,684) = -(rxt(k,462)*y(k,70) + rxt(k,464)*y(k,69) + rxt(k,466)*y(k,176) &
                      + rxt(k,467)*y(k,32) + rxt(k,475)*y(k,86))
         mat(k,1749) = -rxt(k,462)*y(k,155)
         mat(k,1275) = -rxt(k,464)*y(k,155)
         mat(k,1845) = -rxt(k,466)*y(k,155)
         mat(k,1139) = -rxt(k,467)*y(k,155)
         mat(k,1028) = -rxt(k,475)*y(k,155)
         mat(k,1555) = rxt(k,463)*y(k,160)
         mat(k,1093) = rxt(k,461)*y(k,145)
         mat(k,1061) = rxt(k,460)*y(k,145)
         mat(k,929) = rxt(k,461)*y(k,90) + rxt(k,460)*y(k,95)
         mat(k,469) = rxt(k,463)*y(k,68)
         mat(k,86) = 2.000_r8*rxt(k,457)*y(k,152)
         mat(k,144) = 2.000_r8*rxt(k,456)*y(k,153)
         mat(k,520) = 2.000_r8*rxt(k,458)*y(k,150)
         mat(k,504) = 2.000_r8*rxt(k,459)*y(k,151)
         mat(k,468) = -(rxt(k,463)*y(k,68) + rxt(k,473)*y(k,86))
         mat(k,1536) = -rxt(k,463)*y(k,160)
         mat(k,1020) = -rxt(k,473)*y(k,160)
         mat(k,1536) = mat(k,1536) + rxt(k,465)*y(k,161)
         mat(k,1738) = rxt(k,462)*y(k,155)
         mat(k,682) = rxt(k,462)*y(k,70)
         mat(k,288) = rxt(k,465)*y(k,68)
         mat(k,287) = -(rxt(k,465)*y(k,68))
         mat(k,1518) = -rxt(k,465)*y(k,161)
         mat(k,1245) = rxt(k,464)*y(k,155) + rxt(k,470)*y(k,162) + rxt(k,468)*y(k,163) &
                      + rxt(k,478)*y(k,164)
         mat(k,680) = rxt(k,464)*y(k,69)
         mat(k,419) = rxt(k,470)*y(k,69)
         mat(k,474) = rxt(k,468)*y(k,69)
         mat(k,214) = rxt(k,478)*y(k,69)
         mat(k,420) = -(rxt(k,470)*y(k,69) + rxt(k,474)*y(k,86))
         mat(k,1254) = -rxt(k,470)*y(k,162)
         mat(k,1015) = -rxt(k,474)*y(k,162)
         mat(k,681) = rxt(k,466)*y(k,176)
         mat(k,475) = rxt(k,469)*y(k,176)
         mat(k,1829) = rxt(k,466)*y(k,155) + rxt(k,469)*y(k,163)
         mat(k,477) = -(rxt(k,468)*y(k,69) + rxt(k,469)*y(k,176) + rxt(k,471)*y(k,86))
         mat(k,1260) = -rxt(k,468)*y(k,163)
         mat(k,1831) = -rxt(k,469)*y(k,163)
         mat(k,1021) = -rxt(k,471)*y(k,163)
         mat(k,1129) = rxt(k,467)*y(k,155)
         mat(k,683) = rxt(k,467)*y(k,32)
         mat(k,213) = -((rxt(k,472) + rxt(k,476)) * y(k,86) + rxt(k,478)*y(k,69))
         mat(k,1003) = -(rxt(k,472) + rxt(k,476)) * y(k,164)
         mat(k,1241) = -rxt(k,478)*y(k,164)
         mat(k,914) = -(rxt(k,479)*y(k,70) + rxt(k,481)*y(k,69) + rxt(k,491)*y(k,90) &
                      + rxt(k,492)*y(k,95))
         mat(k,1759) = -rxt(k,479)*y(k,165)
         mat(k,1290) = -rxt(k,481)*y(k,165)
         mat(k,1099) = -rxt(k,491)*y(k,165)
         mat(k,1068) = -rxt(k,492)*y(k,165)
         mat(k,1194) = rxt(k,486)*y(k,167) + rxt(k,487)*y(k,168) + rxt(k,489)*y(k,169)
         mat(k,1570) = rxt(k,480)*y(k,166)
         mat(k,1039) = rxt(k,521)*y(k,170) + rxt(k,516)*y(k,171) + rxt(k,517)*y(k,172) &
                      + rxt(k,520)*y(k,173) + rxt(k,518)*y(k,174) + rxt(k,519) &
                      *y(k,175)
         mat(k,719) = rxt(k,480)*y(k,68)
         mat(k,515) = rxt(k,486)*y(k,40)
         mat(k,380) = rxt(k,487)*y(k,40)
         mat(k,258) = rxt(k,489)*y(k,40)
         mat(k,771) = rxt(k,521)*y(k,86)
         mat(k,735) = rxt(k,516)*y(k,86)
         mat(k,747) = rxt(k,517)*y(k,86)
         mat(k,759) = rxt(k,520)*y(k,86)
         mat(k,488) = rxt(k,518)*y(k,86)
         mat(k,329) = rxt(k,519)*y(k,86)
         mat(k,718) = -(rxt(k,480)*y(k,68) + rxt(k,482)*y(k,70) + rxt(k,484)*y(k,176) &
                      + rxt(k,485)*y(k,41))
         mat(k,1558) = -rxt(k,480)*y(k,166)
         mat(k,1752) = -rxt(k,482)*y(k,166)
         mat(k,1848) = -rxt(k,484)*y(k,166)
         mat(k,1625) = -rxt(k,485)*y(k,166)
         mat(k,1558) = mat(k,1558) + rxt(k,483)*y(k,167)
         mat(k,1752) = mat(k,1752) + rxt(k,479)*y(k,165)
         mat(k,912) = rxt(k,479)*y(k,70)
         mat(k,514) = rxt(k,483)*y(k,68)
         mat(k,513) = -(rxt(k,483)*y(k,68) + rxt(k,486)*y(k,40))
         mat(k,1540) = -rxt(k,483)*y(k,167)
         mat(k,1182) = -rxt(k,486)*y(k,167)
         mat(k,1263) = rxt(k,481)*y(k,165)
         mat(k,1740) = rxt(k,482)*y(k,166)
         mat(k,911) = rxt(k,481)*y(k,69)
         mat(k,717) = rxt(k,482)*y(k,70)
         mat(k,379) = -(rxt(k,487)*y(k,40) + rxt(k,488)*y(k,32))
         mat(k,1177) = -rxt(k,487)*y(k,168)
         mat(k,1125) = -rxt(k,488)*y(k,168)
         mat(k,1621) = rxt(k,485)*y(k,166)
         mat(k,716) = rxt(k,485)*y(k,41) + rxt(k,484)*y(k,176)
         mat(k,1827) = rxt(k,484)*y(k,166)
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
         mat(k,257) = -(rxt(k,489)*y(k,40) + 4._r8*rxt(k,490)*y(k,169))
         mat(k,1173) = -rxt(k,489)*y(k,169)
         mat(k,1122) = rxt(k,488)*y(k,168)
         mat(k,378) = rxt(k,488)*y(k,32)
         mat(k,770) = -(rxt(k,495)*y(k,69) + rxt(k,497)*y(k,68) + rxt(k,499)*y(k,32) &
                      + rxt(k,500)*y(k,176) + rxt(k,521)*y(k,86))
         mat(k,1282) = -rxt(k,495)*y(k,170)
         mat(k,1562) = -rxt(k,497)*y(k,170)
         mat(k,1145) = -rxt(k,499)*y(k,170)
         mat(k,1852) = -rxt(k,500)*y(k,170)
         mat(k,1032) = -rxt(k,521)*y(k,170)
         mat(k,1094) = rxt(k,491)*y(k,165)
         mat(k,1062) = rxt(k,492)*y(k,165)
         mat(k,913) = rxt(k,491)*y(k,90) + rxt(k,492)*y(k,95)
         mat(k,731) = -(rxt(k,501)*y(k,68) + rxt(k,503)*y(k,32) + rxt(k,511)*y(k,69) &
                      + rxt(k,514)*y(k,176) + rxt(k,516)*y(k,86))
         mat(k,1559) = -rxt(k,501)*y(k,171)
         mat(k,1142) = -rxt(k,503)*y(k,171)
         mat(k,1279) = -rxt(k,511)*y(k,171)
         mat(k,1849) = -rxt(k,514)*y(k,171)
         mat(k,1029) = -rxt(k,516)*y(k,171)
         mat(k,744) = -(rxt(k,505)*y(k,68) + rxt(k,512)*y(k,32) + rxt(k,513)*y(k,176) &
                      + rxt(k,517)*y(k,86))
         mat(k,1560) = -rxt(k,505)*y(k,172)
         mat(k,1143) = -rxt(k,512)*y(k,172)
         mat(k,1850) = -rxt(k,513)*y(k,172)
         mat(k,1030) = -rxt(k,517)*y(k,172)
         mat(k,1280) = rxt(k,495)*y(k,170) + rxt(k,511)*y(k,171) + rxt(k,506)*y(k,173)
         mat(k,768) = rxt(k,495)*y(k,69)
         mat(k,732) = rxt(k,511)*y(k,69)
         mat(k,756) = rxt(k,506)*y(k,69)
         mat(k,757) = -(rxt(k,506)*y(k,69) + rxt(k,507)*y(k,32) + rxt(k,509)*y(k,176) &
                      + rxt(k,520)*y(k,86))
         mat(k,1281) = -rxt(k,506)*y(k,173)
         mat(k,1144) = -rxt(k,507)*y(k,173)
         mat(k,1851) = -rxt(k,509)*y(k,173)
         mat(k,1031) = -rxt(k,520)*y(k,173)
         mat(k,1561) = rxt(k,497)*y(k,170) + rxt(k,501)*y(k,171) + rxt(k,505)*y(k,172)
         mat(k,769) = rxt(k,497)*y(k,68)
         mat(k,733) = rxt(k,501)*y(k,68)
         mat(k,745) = rxt(k,505)*y(k,68)
         mat(k,486) = -(rxt(k,515)*y(k,176) + rxt(k,518)*y(k,86))
         mat(k,1832) = -rxt(k,515)*y(k,174)
         mat(k,1022) = -rxt(k,518)*y(k,174)
         mat(k,1130) = rxt(k,499)*y(k,170) + rxt(k,503)*y(k,171) + rxt(k,512)*y(k,172) &
                      + rxt(k,507)*y(k,173)
         mat(k,766) = rxt(k,499)*y(k,32)
         mat(k,730) = rxt(k,503)*y(k,32)
         mat(k,742) = rxt(k,512)*y(k,32)
         mat(k,754) = rxt(k,507)*y(k,32)
         mat(k,328) = -(rxt(k,519)*y(k,86))
         mat(k,1010) = -rxt(k,519)*y(k,175)
         mat(k,765) = rxt(k,500)*y(k,176)
         mat(k,729) = rxt(k,514)*y(k,176)
         mat(k,741) = rxt(k,513)*y(k,176)
         mat(k,753) = rxt(k,509)*y(k,176)
         mat(k,485) = rxt(k,515)*y(k,176)
         mat(k,1824) = rxt(k,500)*y(k,170) + rxt(k,514)*y(k,171) + rxt(k,513)*y(k,172) &
                      + rxt(k,509)*y(k,173) + rxt(k,515)*y(k,174)
         mat(k,1886) = -(rxt(k,108)*y(k,92) + rxt(k,227)*y(k,39) + rxt(k,275)*y(k,79) &
                      + rxt(k,329)*y(k,101) + rxt(k,341)*y(k,107) + rxt(k,357) &
                      *y(k,113) + rxt(k,369)*y(k,111) + rxt(k,381)*y(k,122) + rxt(k,382) &
                      *y(k,126) + rxt(k,402)*y(k,130) + rxt(k,412)*y(k,138) + rxt(k,413) &
                      *y(k,139) + rxt(k,418)*y(k,141) + rxt(k,422)*y(k,142) + rxt(k,438) &
                      *y(k,146) + rxt(k,439)*y(k,148) + rxt(k,448)*y(k,151) + rxt(k,466) &
                      *y(k,155) + rxt(k,469)*y(k,163) + rxt(k,484)*y(k,166) + rxt(k,500) &
                      *y(k,170) + rxt(k,509)*y(k,173) + rxt(k,513)*y(k,172) + rxt(k,514) &
                      *y(k,171) + rxt(k,515)*y(k,174))
         mat(k,1475) = -rxt(k,108)*y(k,176)
         mat(k,572) = -rxt(k,227)*y(k,176)
         mat(k,119) = -rxt(k,275)*y(k,176)
         mat(k,602) = -rxt(k,329)*y(k,176)
         mat(k,640) = -rxt(k,341)*y(k,176)
         mat(k,615) = -rxt(k,357)*y(k,176)
         mat(k,628) = -rxt(k,369)*y(k,176)
         mat(k,539) = -rxt(k,381)*y(k,176)
         mat(k,403) = -rxt(k,382)*y(k,176)
         mat(k,678) = -rxt(k,402)*y(k,176)
         mat(k,162) = -rxt(k,412)*y(k,176)
         mat(k,37) = -rxt(k,413)*y(k,176)
         mat(k,860) = -rxt(k,418)*y(k,176)
         mat(k,982) = -rxt(k,422)*y(k,176)
         mat(k,715) = -rxt(k,438)*y(k,176)
         mat(k,651) = -rxt(k,439)*y(k,176)
         mat(k,512) = -rxt(k,448)*y(k,176)
         mat(k,691) = -rxt(k,466)*y(k,176)
         mat(k,484) = -rxt(k,469)*y(k,176)
         mat(k,728) = -rxt(k,484)*y(k,176)
         mat(k,776) = -rxt(k,500)*y(k,176)
         mat(k,764) = -rxt(k,509)*y(k,176)
         mat(k,752) = -rxt(k,513)*y(k,176)
         mat(k,740) = -rxt(k,514)*y(k,176)
         mat(k,491) = -rxt(k,515)*y(k,176)
         mat(k,297) = rxt(k,232)*y(k,96)
         mat(k,1617) = rxt(k,250)*y(k,96)
         mat(k,357) = rxt(k,234)*y(k,96)
         mat(k,60) = rxt(k,235)*y(k,96)
         mat(k,225) = rxt(k,237)*y(k,96)
         mat(k,239) = rxt(k,253)*y(k,96)
         mat(k,803) = rxt(k,254)*y(k,96)
         mat(k,1211) = rxt(k,128)*y(k,87) + rxt(k,330)*y(k,103) + rxt(k,332)*y(k,104) &
                      + rxt(k,358)*y(k,114) + rxt(k,360)*y(k,115) + rxt(k,388) &
                      *y(k,124) + rxt(k,389)*y(k,129) + rxt(k,449)*y(k,149) &
                      + rxt(k,451)*y(k,150) + rxt(k,487)*y(k,168) + rxt(k,489) &
                      *y(k,169)
         mat(k,1654) = rxt(k,133)*y(k,96)
         mat(k,364) = rxt(k,134)*y(k,96)
         mat(k,547) = rxt(k,218)*y(k,96)
         mat(k,169) = rxt(k,242)*y(k,96)
         mat(k,963) = (rxt(k,288)+rxt(k,293))*y(k,53) + (rxt(k,281)+rxt(k,287) &
                       +rxt(k,292))*y(k,54) + rxt(k,189)*y(k,96)
         mat(k,1235) = rxt(k,165)*y(k,96)
         mat(k,282) = rxt(k,141)*y(k,96)
         mat(k,449) = (rxt(k,288)+rxt(k,293))*y(k,49)
         mat(k,555) = (rxt(k,281)+rxt(k,287)+rxt(k,292))*y(k,49) + rxt(k,192)*y(k,96)
         mat(k,1312) = rxt(k,470)*y(k,162)
         mat(k,1056) = rxt(k,348)*y(k,109) + rxt(k,408)*y(k,133) + rxt(k,474)*y(k,162) &
                      + rxt(k,519)*y(k,175)
         mat(k,1392) = rxt(k,128)*y(k,40) + rxt(k,135)*y(k,96)
         mat(k,1365) = rxt(k,232)*y(k,16) + rxt(k,250)*y(k,17) + rxt(k,234)*y(k,18) &
                      + rxt(k,235)*y(k,19) + rxt(k,237)*y(k,20) + rxt(k,253)*y(k,22) &
                      + rxt(k,254)*y(k,23) + rxt(k,133)*y(k,41) + rxt(k,134)*y(k,43) &
                      + rxt(k,218)*y(k,45) + rxt(k,242)*y(k,48) + rxt(k,189)*y(k,49) &
                      + rxt(k,165)*y(k,51) + rxt(k,141)*y(k,52) + rxt(k,192)*y(k,54) &
                      + rxt(k,135)*y(k,87) + 2.000_r8*rxt(k,138)*y(k,96)
         mat(k,371) = rxt(k,330)*y(k,40)
         mat(k,245) = rxt(k,332)*y(k,40)
         mat(k,123) = rxt(k,348)*y(k,86)
         mat(k,196) = rxt(k,358)*y(k,40)
         mat(k,339) = rxt(k,360)*y(k,40)
         mat(k,207) = rxt(k,388)*y(k,40)
         mat(k,135) = rxt(k,389)*y(k,40)
         mat(k,139) = rxt(k,408)*y(k,86)
         mat(k,256) = rxt(k,449)*y(k,40)
         mat(k,527) = rxt(k,451)*y(k,40)
         mat(k,426) = rxt(k,470)*y(k,69) + rxt(k,474)*y(k,86)
         mat(k,384) = rxt(k,487)*y(k,40)
         mat(k,262) = rxt(k,489)*y(k,40)
         mat(k,331) = rxt(k,519)*y(k,86)
      end do
      end subroutine nlnmat09
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
         mat(k, 32) = mat(k, 32) + lmat(k, 32)
         mat(k, 34) = mat(k, 34) + lmat(k, 34)
         mat(k, 36) = mat(k, 36) + lmat(k, 36)
         mat(k, 38) = lmat(k, 38)
         mat(k, 39) = lmat(k, 39)
         mat(k, 40) = lmat(k, 40)
         mat(k, 41) = mat(k, 41) + lmat(k, 41)
         mat(k, 42) = mat(k, 42) + lmat(k, 42)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 46) = mat(k, 46) + lmat(k, 46)
         mat(k, 48) = mat(k, 48) + lmat(k, 48)
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 52) = mat(k, 52) + lmat(k, 52)
         mat(k, 53) = mat(k, 53) + lmat(k, 53)
         mat(k, 54) = mat(k, 54) + lmat(k, 54)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 61) = lmat(k, 61)
         mat(k, 62) = lmat(k, 62)
         mat(k, 63) = lmat(k, 63)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 65) = mat(k, 65) + lmat(k, 65)
         mat(k, 66) = mat(k, 66) + lmat(k, 66)
         mat(k, 68) = mat(k, 68) + lmat(k, 68)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 70) = mat(k, 70) + lmat(k, 70)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 84) = lmat(k, 84)
         mat(k, 85) = lmat(k, 85)
         mat(k, 87) = mat(k, 87) + lmat(k, 87)
         mat(k, 91) = mat(k, 91) + lmat(k, 91)
         mat(k, 92) = mat(k, 92) + lmat(k, 92)
         mat(k, 95) = mat(k, 95) + lmat(k, 95)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 99) = mat(k, 99) + lmat(k, 99)
         mat(k, 100) = mat(k, 100) + lmat(k, 100)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = mat(k, 104) + lmat(k, 104)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 108) = lmat(k, 108)
         mat(k, 109) = lmat(k, 109)
         mat(k, 110) = lmat(k, 110)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 117) = lmat(k, 117)
         mat(k, 118) = lmat(k, 118)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 132) = mat(k, 132) + lmat(k, 132)
         mat(k, 136) = mat(k, 136) + lmat(k, 136)
         mat(k, 140) = mat(k, 140) + lmat(k, 140)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = lmat(k, 152)
         mat(k, 153) = lmat(k, 153)
         mat(k, 154) = lmat(k, 154)
         mat(k, 155) = lmat(k, 155)
         mat(k, 156) = lmat(k, 156)
         mat(k, 157) = lmat(k, 157)
         mat(k, 158) = lmat(k, 158)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = lmat(k, 178)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 182) = lmat(k, 182)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 186) = lmat(k, 186)
         mat(k, 187) = mat(k, 187) + lmat(k, 187)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 208) = mat(k, 208) + lmat(k, 208)
         mat(k, 213) = mat(k, 213) + lmat(k, 213)
         mat(k, 218) = mat(k, 218) + lmat(k, 218)
         mat(k, 220) = lmat(k, 220)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 226) = mat(k, 226) + lmat(k, 226)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 228) = lmat(k, 228)
         mat(k, 234) = mat(k, 234) + lmat(k, 234)
         mat(k, 236) = lmat(k, 236)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 240) = mat(k, 240) + lmat(k, 240)
         mat(k, 241) = mat(k, 241) + lmat(k, 241)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 244) = lmat(k, 244)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 261) = lmat(k, 261)
         mat(k, 263) = mat(k, 263) + lmat(k, 263)
         mat(k, 264) = mat(k, 264) + lmat(k, 264)
         mat(k, 265) = lmat(k, 265)
         mat(k, 266) = lmat(k, 266)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 268) = lmat(k, 268)
         mat(k, 269) = lmat(k, 269)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 271) = lmat(k, 271)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 278) = mat(k, 278) + lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = mat(k, 280) + lmat(k, 280)
         mat(k, 281) = lmat(k, 281)
         mat(k, 283) = mat(k, 283) + lmat(k, 283)
         mat(k, 287) = mat(k, 287) + lmat(k, 287)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 293) = mat(k, 293) + lmat(k, 293)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 321) = lmat(k, 321)
         mat(k, 322) = lmat(k, 322)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 326) = mat(k, 326) + lmat(k, 326)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 335) = mat(k, 335) + lmat(k, 335)
         mat(k, 337) = lmat(k, 337)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 350) = lmat(k, 350)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 360) = mat(k, 360) + lmat(k, 360)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 367) = mat(k, 367) + lmat(k, 367)
         mat(k, 370) = lmat(k, 370)
         mat(k, 372) = mat(k, 372) + lmat(k, 372)
         mat(k, 379) = mat(k, 379) + lmat(k, 379)
         mat(k, 380) = mat(k, 380) + lmat(k, 380)
         mat(k, 383) = lmat(k, 383)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 391) = lmat(k, 391)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 399) = mat(k, 399) + lmat(k, 399)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 412) = mat(k, 412) + lmat(k, 412)
         mat(k, 420) = mat(k, 420) + lmat(k, 420)
         mat(k, 428) = mat(k, 428) + lmat(k, 428)
         mat(k, 429) = lmat(k, 429)
         mat(k, 433) = mat(k, 433) + lmat(k, 433)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 444) = lmat(k, 444)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 450) = lmat(k, 450)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 460) = mat(k, 460) + lmat(k, 460)
         mat(k, 468) = mat(k, 468) + lmat(k, 468)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 486) = mat(k, 486) + lmat(k, 486)
         mat(k, 487) = lmat(k, 487)
         mat(k, 490) = mat(k, 490) + lmat(k, 490)
         mat(k, 492) = mat(k, 492) + lmat(k, 492)
         mat(k, 497) = mat(k, 497) + lmat(k, 497)
         mat(k, 506) = mat(k, 506) + lmat(k, 506)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 515) = mat(k, 515) + lmat(k, 515)
         mat(k, 517) = mat(k, 517) + lmat(k, 517)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 540) = mat(k, 540) + lmat(k, 540)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 549) = mat(k, 549) + lmat(k, 549)
         mat(k, 551) = mat(k, 551) + lmat(k, 551)
         mat(k, 552) = mat(k, 552) + lmat(k, 552)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 561) = mat(k, 561) + lmat(k, 561)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 575) = lmat(k, 575)
         mat(k, 578) = lmat(k, 578)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 582) = mat(k, 582) + lmat(k, 582)
         mat(k, 584) = mat(k, 584) + lmat(k, 584)
         mat(k, 586) = lmat(k, 586)
         mat(k, 587) = mat(k, 587) + lmat(k, 587)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 590) = mat(k, 590) + lmat(k, 590)
         mat(k, 593) = mat(k, 593) + lmat(k, 593)
         mat(k, 594) = mat(k, 594) + lmat(k, 594)
         mat(k, 599) = mat(k, 599) + lmat(k, 599)
         mat(k, 606) = mat(k, 606) + lmat(k, 606)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 633) = mat(k, 633) + lmat(k, 633)
         mat(k, 644) = mat(k, 644) + lmat(k, 644)
         mat(k, 655) = mat(k, 655) + lmat(k, 655)
         mat(k, 656) = lmat(k, 656)
         mat(k, 658) = lmat(k, 658)
         mat(k, 668) = lmat(k, 668)
         mat(k, 671) = mat(k, 671) + lmat(k, 671)
         mat(k, 679) = lmat(k, 679)
         mat(k, 684) = mat(k, 684) + lmat(k, 684)
         mat(k, 694) = mat(k, 694) + lmat(k, 694)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 718) = mat(k, 718) + lmat(k, 718)
         mat(k, 731) = mat(k, 731) + lmat(k, 731)
         mat(k, 734) = lmat(k, 734)
         mat(k, 743) = lmat(k, 743)
         mat(k, 744) = mat(k, 744) + lmat(k, 744)
         mat(k, 746) = lmat(k, 746)
         mat(k, 750) = mat(k, 750) + lmat(k, 750)
         mat(k, 755) = lmat(k, 755)
         mat(k, 757) = mat(k, 757) + lmat(k, 757)
         mat(k, 758) = lmat(k, 758)
         mat(k, 763) = mat(k, 763) + lmat(k, 763)
         mat(k, 767) = lmat(k, 767)
         mat(k, 770) = mat(k, 770) + lmat(k, 770)
         mat(k, 778) = mat(k, 778) + lmat(k, 778)
         mat(k, 787) = lmat(k, 787)
         mat(k, 789) = mat(k, 789) + lmat(k, 789)
         mat(k, 790) = mat(k, 790) + lmat(k, 790)
         mat(k, 792) = lmat(k, 792)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 795) = mat(k, 795) + lmat(k, 795)
         mat(k, 799) = lmat(k, 799)
         mat(k, 800) = mat(k, 800) + lmat(k, 800)
         mat(k, 801) = mat(k, 801) + lmat(k, 801)
         mat(k, 803) = mat(k, 803) + lmat(k, 803)
         mat(k, 804) = mat(k, 804) + lmat(k, 804)
         mat(k, 806) = mat(k, 806) + lmat(k, 806)
         mat(k, 807) = lmat(k, 807)
         mat(k, 819) = mat(k, 819) + lmat(k, 819)
         mat(k, 831) = mat(k, 831) + lmat(k, 831)
         mat(k, 832) = mat(k, 832) + lmat(k, 832)
         mat(k, 833) = lmat(k, 833)
         mat(k, 849) = mat(k, 849) + lmat(k, 849)
         mat(k, 862) = lmat(k, 862)
         mat(k, 864) = mat(k, 864) + lmat(k, 864)
         mat(k, 865) = lmat(k, 865)
         mat(k, 878) = lmat(k, 878)
         mat(k, 880) = mat(k, 880) + lmat(k, 880)
         mat(k, 888) = mat(k, 888) + lmat(k, 888)
         mat(k, 895) = mat(k, 895) + lmat(k, 895)
         mat(k, 898) = mat(k, 898) + lmat(k, 898)
         mat(k, 899) = lmat(k, 899)
         mat(k, 913) = mat(k, 913) + lmat(k, 913)
         mat(k, 914) = mat(k, 914) + lmat(k, 914)
         mat(k, 915) = lmat(k, 915)
         mat(k, 929) = mat(k, 929) + lmat(k, 929)
         mat(k, 933) = mat(k, 933) + lmat(k, 933)
         mat(k, 935) = lmat(k, 935)
         mat(k, 951) = mat(k, 951) + lmat(k, 951)
         mat(k, 953) = mat(k, 953) + lmat(k, 953)
         mat(k, 958) = mat(k, 958) + lmat(k, 958)
         mat(k, 970) = mat(k, 970) + lmat(k, 970)
         mat(k, 985) = mat(k, 985) + lmat(k, 985)
         mat(k,1042) = mat(k,1042) + lmat(k,1042)
         mat(k,1072) = mat(k,1072) + lmat(k,1072)
         mat(k,1073) = mat(k,1073) + lmat(k,1073)
         mat(k,1083) = mat(k,1083) + lmat(k,1083)
         mat(k,1104) = mat(k,1104) + lmat(k,1104)
         mat(k,1126) = mat(k,1126) + lmat(k,1126)
         mat(k,1156) = mat(k,1156) + lmat(k,1156)
         mat(k,1163) = mat(k,1163) + lmat(k,1163)
         mat(k,1201) = mat(k,1201) + lmat(k,1201)
         mat(k,1220) = mat(k,1220) + lmat(k,1220)
         mat(k,1222) = mat(k,1222) + lmat(k,1222)
         mat(k,1224) = lmat(k,1224)
         mat(k,1242) = lmat(k,1242)
         mat(k,1243) = mat(k,1243) + lmat(k,1243)
         mat(k,1273) = mat(k,1273) + lmat(k,1273)
         mat(k,1294) = lmat(k,1294)
         mat(k,1295) = mat(k,1295) + lmat(k,1295)
         mat(k,1299) = mat(k,1299) + lmat(k,1299)
         mat(k,1303) = mat(k,1303) + lmat(k,1303)
         mat(k,1305) = mat(k,1305) + lmat(k,1305)
         mat(k,1352) = mat(k,1352) + lmat(k,1352)
         mat(k,1368) = mat(k,1368) + lmat(k,1368)
         mat(k,1380) = mat(k,1380) + lmat(k,1380)
         mat(k,1409) = mat(k,1409) + lmat(k,1409)
         mat(k,1411) = mat(k,1411) + lmat(k,1411)
         mat(k,1413) = mat(k,1413) + lmat(k,1413)
         mat(k,1416) = mat(k,1416) + lmat(k,1416)
         mat(k,1421) = mat(k,1421) + lmat(k,1421)
         mat(k,1465) = mat(k,1465) + lmat(k,1465)
         mat(k,1467) = mat(k,1467) + lmat(k,1467)
         mat(k,1497) = mat(k,1497) + lmat(k,1497)
         mat(k,1515) = mat(k,1515) + lmat(k,1515)
         mat(k,1516) = mat(k,1516) + lmat(k,1516)
         mat(k,1553) = mat(k,1553) + lmat(k,1553)
         mat(k,1575) = lmat(k,1575)
         mat(k,1587) = mat(k,1587) + lmat(k,1587)
         mat(k,1596) = mat(k,1596) + lmat(k,1596)
         mat(k,1601) = mat(k,1601) + lmat(k,1601)
         mat(k,1610) = mat(k,1610) + lmat(k,1610)
         mat(k,1611) = lmat(k,1611)
         mat(k,1648) = mat(k,1648) + lmat(k,1648)
         mat(k,1663) = mat(k,1663) + lmat(k,1663)
         mat(k,1664) = mat(k,1664) + lmat(k,1664)
         mat(k,1667) = mat(k,1667) + lmat(k,1667)
         mat(k,1670) = mat(k,1670) + lmat(k,1670)
         mat(k,1673) = mat(k,1673) + lmat(k,1673)
         mat(k,1675) = mat(k,1675) + lmat(k,1675)
         mat(k,1686) = mat(k,1686) + lmat(k,1686)
         mat(k,1696) = mat(k,1696) + lmat(k,1696)
         mat(k,1700) = mat(k,1700) + lmat(k,1700)
         mat(k,1706) = mat(k,1706) + lmat(k,1706)
         mat(k,1708) = lmat(k,1708)
         mat(k,1710) = mat(k,1710) + lmat(k,1710)
         mat(k,1719) = mat(k,1719) + lmat(k,1719)
         mat(k,1724) = mat(k,1724) + lmat(k,1724)
         mat(k,1742) = mat(k,1742) + lmat(k,1742)
         mat(k,1768) = mat(k,1768) + lmat(k,1768)
         mat(k,1772) = mat(k,1772) + lmat(k,1772)
         mat(k,1774) = mat(k,1774) + lmat(k,1774)
         mat(k,1780) = mat(k,1780) + lmat(k,1780)
         mat(k,1802) = mat(k,1802) + lmat(k,1802)
         mat(k,1803) = mat(k,1803) + lmat(k,1803)
         mat(k,1810) = mat(k,1810) + lmat(k,1810)
         mat(k,1870) = mat(k,1870) + lmat(k,1870)
         mat(k,1873) = mat(k,1873) + lmat(k,1873)
         mat(k,1876) = mat(k,1876) + lmat(k,1876)
         mat(k,1878) = mat(k,1878) + lmat(k,1878)
         mat(k,1880) = lmat(k,1880)
         mat(k,1886) = mat(k,1886) + lmat(k,1886)
         mat(k, 147) = 0._r8
         mat(k, 165) = 0._r8
         mat(k, 205) = 0._r8
         mat(k, 312) = 0._r8
         mat(k, 400) = 0._r8
         mat(k, 421) = 0._r8
         mat(k, 425) = 0._r8
         mat(k, 446) = 0._r8
         mat(k, 457) = 0._r8
         mat(k, 466) = 0._r8
         mat(k, 476) = 0._r8
         mat(k, 478) = 0._r8
         mat(k, 483) = 0._r8
         mat(k, 526) = 0._r8
         mat(k, 535) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 591) = 0._r8
         mat(k, 595) = 0._r8
         mat(k, 607) = 0._r8
         mat(k, 608) = 0._r8
         mat(k, 611) = 0._r8
         mat(k, 613) = 0._r8
         mat(k, 614) = 0._r8
         mat(k, 618) = 0._r8
         mat(k, 622) = 0._r8
         mat(k, 624) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 646) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 662) = 0._r8
         mat(k, 663) = 0._r8
         mat(k, 665) = 0._r8
         mat(k, 676) = 0._r8
         mat(k, 689) = 0._r8
         mat(k, 696) = 0._r8
         mat(k, 697) = 0._r8
         mat(k, 702) = 0._r8
         mat(k, 705) = 0._r8
         mat(k, 711) = 0._r8
         mat(k, 720) = 0._r8
         mat(k, 724) = 0._r8
         mat(k, 785) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 802) = 0._r8
         mat(k, 810) = 0._r8
         mat(k, 812) = 0._r8
         mat(k, 813) = 0._r8
         mat(k, 814) = 0._r8
         mat(k, 817) = 0._r8
         mat(k, 820) = 0._r8
         mat(k, 822) = 0._r8
         mat(k, 829) = 0._r8
         mat(k, 836) = 0._r8
         mat(k, 837) = 0._r8
         mat(k, 839) = 0._r8
         mat(k, 842) = 0._r8
         mat(k, 847) = 0._r8
         mat(k, 850) = 0._r8
         mat(k, 854) = 0._r8
         mat(k, 857) = 0._r8
         mat(k, 863) = 0._r8
         mat(k, 868) = 0._r8
         mat(k, 872) = 0._r8
         mat(k, 874) = 0._r8
         mat(k, 885) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 903) = 0._r8
         mat(k, 905) = 0._r8
         mat(k, 906) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 918) = 0._r8
         mat(k, 919) = 0._r8
         mat(k, 921) = 0._r8
         mat(k, 922) = 0._r8
         mat(k, 923) = 0._r8
         mat(k, 924) = 0._r8
         mat(k, 927) = 0._r8
         mat(k, 932) = 0._r8
         mat(k, 934) = 0._r8
         mat(k, 938) = 0._r8
         mat(k, 939) = 0._r8
         mat(k, 941) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 952) = 0._r8
         mat(k, 956) = 0._r8
         mat(k, 960) = 0._r8
         mat(k, 961) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 969) = 0._r8
         mat(k, 973) = 0._r8
         mat(k, 974) = 0._r8
         mat(k, 980) = 0._r8
         mat(k, 986) = 0._r8
         mat(k, 987) = 0._r8
         mat(k, 989) = 0._r8
         mat(k, 991) = 0._r8
         mat(k, 992) = 0._r8
         mat(k, 996) = 0._r8
         mat(k,1006) = 0._r8
         mat(k,1007) = 0._r8
         mat(k,1024) = 0._r8
         mat(k,1033) = 0._r8
         mat(k,1046) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1050) = 0._r8
         mat(k,1053) = 0._r8
         mat(k,1054) = 0._r8
         mat(k,1055) = 0._r8
         mat(k,1059) = 0._r8
         mat(k,1074) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1077) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1079) = 0._r8
         mat(k,1082) = 0._r8
         mat(k,1084) = 0._r8
         mat(k,1085) = 0._r8
         mat(k,1089) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1105) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1107) = 0._r8
         mat(k,1108) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1110) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1115) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1132) = 0._r8
         mat(k,1140) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1149) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1153) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1164) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1166) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1197) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1205) = 0._r8
         mat(k,1206) = 0._r8
         mat(k,1209) = 0._r8
         mat(k,1214) = 0._r8
         mat(k,1215) = 0._r8
         mat(k,1216) = 0._r8
         mat(k,1217) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1219) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1223) = 0._r8
         mat(k,1225) = 0._r8
         mat(k,1226) = 0._r8
         mat(k,1227) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1229) = 0._r8
         mat(k,1231) = 0._r8
         mat(k,1232) = 0._r8
         mat(k,1233) = 0._r8
         mat(k,1234) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1271) = 0._r8
         mat(k,1278) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1285) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1307) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1334) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1355) = 0._r8
         mat(k,1375) = 0._r8
         mat(k,1377) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1385) = 0._r8
         mat(k,1397) = 0._r8
         mat(k,1402) = 0._r8
         mat(k,1403) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1405) = 0._r8
         mat(k,1406) = 0._r8
         mat(k,1407) = 0._r8
         mat(k,1408) = 0._r8
         mat(k,1414) = 0._r8
         mat(k,1415) = 0._r8
         mat(k,1417) = 0._r8
         mat(k,1418) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1447) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1451) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1457) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1460) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1470) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1502) = 0._r8
         mat(k,1503) = 0._r8
         mat(k,1506) = 0._r8
         mat(k,1529) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1546) = 0._r8
         mat(k,1563) = 0._r8
         mat(k,1576) = 0._r8
         mat(k,1580) = 0._r8
         mat(k,1585) = 0._r8
         mat(k,1595) = 0._r8
         mat(k,1600) = 0._r8
         mat(k,1603) = 0._r8
         mat(k,1606) = 0._r8
         mat(k,1607) = 0._r8
         mat(k,1613) = 0._r8
         mat(k,1614) = 0._r8
         mat(k,1615) = 0._r8
         mat(k,1616) = 0._r8
         mat(k,1622) = 0._r8
         mat(k,1626) = 0._r8
         mat(k,1627) = 0._r8
         mat(k,1628) = 0._r8
         mat(k,1629) = 0._r8
         mat(k,1630) = 0._r8
         mat(k,1633) = 0._r8
         mat(k,1634) = 0._r8
         mat(k,1635) = 0._r8
         mat(k,1636) = 0._r8
         mat(k,1637) = 0._r8
         mat(k,1639) = 0._r8
         mat(k,1640) = 0._r8
         mat(k,1642) = 0._r8
         mat(k,1643) = 0._r8
         mat(k,1647) = 0._r8
         mat(k,1649) = 0._r8
         mat(k,1650) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1652) = 0._r8
         mat(k,1653) = 0._r8
         mat(k,1659) = 0._r8
         mat(k,1660) = 0._r8
         mat(k,1661) = 0._r8
         mat(k,1662) = 0._r8
         mat(k,1668) = 0._r8
         mat(k,1669) = 0._r8
         mat(k,1672) = 0._r8
         mat(k,1674) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1677) = 0._r8
         mat(k,1678) = 0._r8
         mat(k,1685) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1688) = 0._r8
         mat(k,1689) = 0._r8
         mat(k,1694) = 0._r8
         mat(k,1697) = 0._r8
         mat(k,1698) = 0._r8
         mat(k,1699) = 0._r8
         mat(k,1702) = 0._r8
         mat(k,1704) = 0._r8
         mat(k,1711) = 0._r8
         mat(k,1712) = 0._r8
         mat(k,1717) = 0._r8
         mat(k,1721) = 0._r8
         mat(k,1727) = 0._r8
         mat(k,1734) = 0._r8
         mat(k,1763) = 0._r8
         mat(k,1765) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1775) = 0._r8
         mat(k,1776) = 0._r8
         mat(k,1782) = 0._r8
         mat(k,1794) = 0._r8
         mat(k,1795) = 0._r8
         mat(k,1796) = 0._r8
         mat(k,1801) = 0._r8
         mat(k,1805) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1809) = 0._r8
         mat(k,1811) = 0._r8
         mat(k,1825) = 0._r8
         mat(k,1830) = 0._r8
         mat(k,1833) = 0._r8
         mat(k,1835) = 0._r8
         mat(k,1838) = 0._r8
         mat(k,1846) = 0._r8
         mat(k,1853) = 0._r8
         mat(k,1854) = 0._r8
         mat(k,1855) = 0._r8
         mat(k,1856) = 0._r8
         mat(k,1857) = 0._r8
         mat(k,1859) = 0._r8
         mat(k,1860) = 0._r8
         mat(k,1861) = 0._r8
         mat(k,1862) = 0._r8
         mat(k,1863) = 0._r8
         mat(k,1865) = 0._r8
         mat(k,1866) = 0._r8
         mat(k,1867) = 0._r8
         mat(k,1868) = 0._r8
         mat(k,1871) = 0._r8
         mat(k,1874) = 0._r8
         mat(k,1875) = 0._r8
         mat(k,1877) = 0._r8
         mat(k,1879) = 0._r8
         mat(k,1881) = 0._r8
         mat(k,1882) = 0._r8
         mat(k,1883) = 0._r8
         mat(k,1884) = 0._r8
         mat(k,1885) = 0._r8
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
         mat(k, 36) = mat(k, 36) - dti(k)
         mat(k, 38) = mat(k, 38) - dti(k)
         mat(k, 41) = mat(k, 41) - dti(k)
         mat(k, 45) = mat(k, 45) - dti(k)
         mat(k, 49) = mat(k, 49) - dti(k)
         mat(k, 53) = mat(k, 53) - dti(k)
         mat(k, 57) = mat(k, 57) - dti(k)
         mat(k, 61) = mat(k, 61) - dti(k)
         mat(k, 64) = mat(k, 64) - dti(k)
         mat(k, 69) = mat(k, 69) - dti(k)
         mat(k, 74) = mat(k, 74) - dti(k)
         mat(k, 79) = mat(k, 79) - dti(k)
         mat(k, 84) = mat(k, 84) - dti(k)
         mat(k, 87) = mat(k, 87) - dti(k)
         mat(k, 91) = mat(k, 91) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 99) = mat(k, 99) - dti(k)
         mat(k, 104) = mat(k, 104) - dti(k)
         mat(k, 108) = mat(k, 108) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 116) = mat(k, 116) - dti(k)
         mat(k, 120) = mat(k, 120) - dti(k)
         mat(k, 124) = mat(k, 124) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 132) = mat(k, 132) - dti(k)
         mat(k, 136) = mat(k, 136) - dti(k)
         mat(k, 140) = mat(k, 140) - dti(k)
         mat(k, 146) = mat(k, 146) - dti(k)
         mat(k, 150) = mat(k, 150) - dti(k)
         mat(k, 153) = mat(k, 153) - dti(k)
         mat(k, 160) = mat(k, 160) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 170) = mat(k, 170) - dti(k)
         mat(k, 176) = mat(k, 176) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 184) = mat(k, 184) - dti(k)
         mat(k, 189) = mat(k, 189) - dti(k)
         mat(k, 193) = mat(k, 193) - dti(k)
         mat(k, 197) = mat(k, 197) - dti(k)
         mat(k, 204) = mat(k, 204) - dti(k)
         mat(k, 208) = mat(k, 208) - dti(k)
         mat(k, 213) = mat(k, 213) - dti(k)
         mat(k, 218) = mat(k, 218) - dti(k)
         mat(k, 226) = mat(k, 226) - dti(k)
         mat(k, 234) = mat(k, 234) - dti(k)
         mat(k, 240) = mat(k, 240) - dti(k)
         mat(k, 246) = mat(k, 246) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 257) = mat(k, 257) - dti(k)
         mat(k, 263) = mat(k, 263) - dti(k)
         mat(k, 270) = mat(k, 270) - dti(k)
         mat(k, 276) = mat(k, 276) - dti(k)
         mat(k, 283) = mat(k, 283) - dti(k)
         mat(k, 287) = mat(k, 287) - dti(k)
         mat(k, 291) = mat(k, 291) - dti(k)
         mat(k, 298) = mat(k, 298) - dti(k)
         mat(k, 304) = mat(k, 304) - dti(k)
         mat(k, 311) = mat(k, 311) - dti(k)
         mat(k, 316) = mat(k, 316) - dti(k)
         mat(k, 320) = mat(k, 320) - dti(k)
         mat(k, 328) = mat(k, 328) - dti(k)
         mat(k, 333) = mat(k, 333) - dti(k)
         mat(k, 340) = mat(k, 340) - dti(k)
         mat(k, 349) = mat(k, 349) - dti(k)
         mat(k, 358) = mat(k, 358) - dti(k)
         mat(k, 366) = mat(k, 366) - dti(k)
         mat(k, 372) = mat(k, 372) - dti(k)
         mat(k, 379) = mat(k, 379) - dti(k)
         mat(k, 385) = mat(k, 385) - dti(k)
         mat(k, 389) = mat(k, 389) - dti(k)
         mat(k, 399) = mat(k, 399) - dti(k)
         mat(k, 404) = mat(k, 404) - dti(k)
         mat(k, 412) = mat(k, 412) - dti(k)
         mat(k, 420) = mat(k, 420) - dti(k)
         mat(k, 428) = mat(k, 428) - dti(k)
         mat(k, 433) = mat(k, 433) - dti(k)
         mat(k, 442) = mat(k, 442) - dti(k)
         mat(k, 453) = mat(k, 453) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 468) = mat(k, 468) - dti(k)
         mat(k, 477) = mat(k, 477) - dti(k)
         mat(k, 486) = mat(k, 486) - dti(k)
         mat(k, 492) = mat(k, 492) - dti(k)
         mat(k, 497) = mat(k, 497) - dti(k)
         mat(k, 506) = mat(k, 506) - dti(k)
         mat(k, 513) = mat(k, 513) - dti(k)
         mat(k, 522) = mat(k, 522) - dti(k)
         mat(k, 532) = mat(k, 532) - dti(k)
         mat(k, 540) = mat(k, 540) - dti(k)
         mat(k, 549) = mat(k, 549) - dti(k)
         mat(k, 556) = mat(k, 556) - dti(k)
         mat(k, 564) = mat(k, 564) - dti(k)
         mat(k, 574) = mat(k, 574) - dti(k)
         mat(k, 582) = mat(k, 582) - dti(k)
         mat(k, 593) = mat(k, 593) - dti(k)
         mat(k, 606) = mat(k, 606) - dti(k)
         mat(k, 620) = mat(k, 620) - dti(k)
         mat(k, 633) = mat(k, 633) - dti(k)
         mat(k, 644) = mat(k, 644) - dti(k)
         mat(k, 655) = mat(k, 655) - dti(k)
         mat(k, 671) = mat(k, 671) - dti(k)
         mat(k, 684) = mat(k, 684) - dti(k)
         mat(k, 694) = mat(k, 694) - dti(k)
         mat(k, 708) = mat(k, 708) - dti(k)
         mat(k, 718) = mat(k, 718) - dti(k)
         mat(k, 731) = mat(k, 731) - dti(k)
         mat(k, 744) = mat(k, 744) - dti(k)
         mat(k, 757) = mat(k, 757) - dti(k)
         mat(k, 770) = mat(k, 770) - dti(k)
         mat(k, 778) = mat(k, 778) - dti(k)
         mat(k, 789) = mat(k, 789) - dti(k)
         mat(k, 806) = mat(k, 806) - dti(k)
         mat(k, 819) = mat(k, 819) - dti(k)
         mat(k, 832) = mat(k, 832) - dti(k)
         mat(k, 849) = mat(k, 849) - dti(k)
         mat(k, 864) = mat(k, 864) - dti(k)
         mat(k, 880) = mat(k, 880) - dti(k)
         mat(k, 898) = mat(k, 898) - dti(k)
         mat(k, 914) = mat(k, 914) - dti(k)
         mat(k, 933) = mat(k, 933) - dti(k)
         mat(k, 951) = mat(k, 951) - dti(k)
         mat(k, 970) = mat(k, 970) - dti(k)
         mat(k, 985) = mat(k, 985) - dti(k)
         mat(k,1042) = mat(k,1042) - dti(k)
         mat(k,1072) = mat(k,1072) - dti(k)
         mat(k,1104) = mat(k,1104) - dti(k)
         mat(k,1156) = mat(k,1156) - dti(k)
         mat(k,1201) = mat(k,1201) - dti(k)
         mat(k,1220) = mat(k,1220) - dti(k)
         mat(k,1299) = mat(k,1299) - dti(k)
         mat(k,1352) = mat(k,1352) - dti(k)
         mat(k,1380) = mat(k,1380) - dti(k)
         mat(k,1413) = mat(k,1413) - dti(k)
         mat(k,1465) = mat(k,1465) - dti(k)
         mat(k,1497) = mat(k,1497) - dti(k)
         mat(k,1587) = mat(k,1587) - dti(k)
         mat(k,1610) = mat(k,1610) - dti(k)
         mat(k,1648) = mat(k,1648) - dti(k)
         mat(k,1673) = mat(k,1673) - dti(k)
         mat(k,1700) = mat(k,1700) - dti(k)
         mat(k,1724) = mat(k,1724) - dti(k)
         mat(k,1780) = mat(k,1780) - dti(k)
         mat(k,1810) = mat(k,1810) - dti(k)
         mat(k,1886) = mat(k,1886) - dti(k)
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
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
