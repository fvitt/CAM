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
         mat(k,820) = rxt(k,293)*y(k,26)
         mat(k,460) = rxt(k,293)*y(k,4)
         mat(k,382) = (rxt(k,403)+rxt(k,408))*y(k,51)
         mat(k,228) = (rxt(k,403)+rxt(k,408))*y(k,47)
         mat(k,849) = -(4._r8*rxt(k,290)*y(k,4) + (rxt(k,291) + rxt(k,292) + rxt(k,293) &
                      ) * y(k,26) + rxt(k,294)*y(k,100) + rxt(k,295)*y(k,59) + rxt(k,296) &
                      *y(k,60) + rxt(k,299)*y(k,66) + rxt(k,300)*y(k,109) + rxt(k,377) &
                      *y(k,74))
         mat(k,489) = -(rxt(k,291) + rxt(k,292) + rxt(k,293)) * y(k,4)
         mat(k,520) = -rxt(k,294)*y(k,4)
         mat(k,770) = -rxt(k,295)*y(k,4)
         mat(k,885) = -rxt(k,296)*y(k,4)
         mat(k,571) = -rxt(k,299)*y(k,4)
         mat(k,817) = -rxt(k,300)*y(k,4)
         mat(k,352) = -rxt(k,377)*y(k,4)
         mat(k,152) = rxt(k,297)*y(k,66)
         mat(k,302) = rxt(k,314)*y(k,105)
         mat(k,236) = rxt(k,308)*y(k,66)
         mat(k,571) = mat(k,571) + rxt(k,297)*y(k,5) + rxt(k,308)*y(k,51)
         mat(k,719) = rxt(k,289)*y(k,97)
         mat(k,379) = rxt(k,289)*y(k,68)
         mat(k,692) = rxt(k,314)*y(k,43)
         mat(k,146) = -(rxt(k,297)*y(k,66))
         mat(k,535) = -rxt(k,297)*y(k,5)
         mat(k,825) = rxt(k,296)*y(k,60)
         mat(k,859) = rxt(k,296)*y(k,4)
         mat(k,589) = -(rxt(k,232)*y(k,98) + rxt(k,287)*y(k,97) + rxt(k,352)*y(k,61) &
                      + rxt(k,353)*y(k,66) + rxt(k,354)*y(k,109))
         mat(k,620) = -rxt(k,232)*y(k,16)
         mat(k,373) = -rxt(k,287)*y(k,16)
         mat(k,447) = -rxt(k,352)*y(k,16)
         mat(k,563) = -rxt(k,353)*y(k,16)
         mat(k,809) = -rxt(k,354)*y(k,16)
         mat(k,419) = rxt(k,239)*y(k,26) + rxt(k,358)*y(k,59)
         mat(k,115) = .300_r8*rxt(k,360)*y(k,109)
         mat(k,481) = rxt(k,239)*y(k,20)
         mat(k,762) = rxt(k,358)*y(k,20)
         mat(k,809) = mat(k,809) + .300_r8*rxt(k,360)*y(k,21)
         mat(k,415) = -(rxt(k,239)*y(k,26) + rxt(k,357)*y(k,100) + rxt(k,358)*y(k,59))
         mat(k,476) = -rxt(k,239)*y(k,20)
         mat(k,507) = -rxt(k,357)*y(k,20)
         mat(k,757) = -rxt(k,358)*y(k,20)
         mat(k,114) = .700_r8*rxt(k,360)*y(k,109)
         mat(k,804) = .700_r8*rxt(k,360)*y(k,21)
         mat(k,112) = -(rxt(k,360)*y(k,109))
         mat(k,783) = -rxt(k,360)*y(k,21)
         mat(k,411) = rxt(k,357)*y(k,100)
         mat(k,495) = rxt(k,357)*y(k,20)
         mat(k,459) = 2.000_r8*rxt(k,241)*y(k,26)
         mat(k,311) = (rxt(k,401)+rxt(k,406)+rxt(k,411))*y(k,47) + rxt(k,245)*y(k,98)
         mat(k,381) = (rxt(k,401)+rxt(k,406)+rxt(k,411))*y(k,27) + (rxt(k,396) &
                       +rxt(k,402)+rxt(k,407))*y(k,52)
         mat(k,285) = (rxt(k,396)+rxt(k,402)+rxt(k,407))*y(k,47)
         mat(k,600) = rxt(k,245)*y(k,27)
         mat(k,458) = 2.000_r8*rxt(k,273)*y(k,26)
         mat(k,478) = -(rxt(k,239)*y(k,20) + (4._r8*rxt(k,240) + 4._r8*rxt(k,241) &
                      + 4._r8*rxt(k,242) + 4._r8*rxt(k,273)) * y(k,26) + rxt(k,243) &
                      *y(k,100) + rxt(k,244)*y(k,59) + rxt(k,246)*y(k,60) + rxt(k,250) &
                      *y(k,66) + (rxt(k,251) + rxt(k,252)) * y(k,109) + (rxt(k,291) &
                      + rxt(k,292) + rxt(k,293)) * y(k,4) + rxt(k,378)*y(k,74))
         mat(k,416) = -rxt(k,239)*y(k,26)
         mat(k,509) = -rxt(k,243)*y(k,26)
         mat(k,759) = -rxt(k,244)*y(k,26)
         mat(k,874) = -rxt(k,246)*y(k,26)
         mat(k,560) = -rxt(k,250)*y(k,26)
         mat(k,806) = -(rxt(k,251) + rxt(k,252)) * y(k,26)
         mat(k,838) = -(rxt(k,291) + rxt(k,292) + rxt(k,293)) * y(k,26)
         mat(k,344) = -rxt(k,378)*y(k,26)
         mat(k,318) = rxt(k,247)*y(k,66)
         mat(k,398) = rxt(k,272)*y(k,105)
         mat(k,289) = rxt(k,261)*y(k,66) + rxt(k,260)*y(k,98) + rxt(k,262)*y(k,109)
         mat(k,560) = mat(k,560) + rxt(k,247)*y(k,27) + rxt(k,261)*y(k,52)
         mat(k,708) = rxt(k,238)*y(k,98)
         mat(k,89) = rxt(k,383)*y(k,74)
         mat(k,344) = mat(k,344) + rxt(k,383)*y(k,69)
         mat(k,617) = rxt(k,260)*y(k,52) + rxt(k,238)*y(k,68) + rxt(k,237)*y(k,100)
         mat(k,509) = mat(k,509) + rxt(k,237)*y(k,98)
         mat(k,681) = rxt(k,272)*y(k,47)
         mat(k,806) = mat(k,806) + rxt(k,262)*y(k,52)
         mat(k,314) = -(rxt(k,245)*y(k,98) + rxt(k,247)*y(k,66) + rxt(k,248)*y(k,109) &
                      + (rxt(k,401) + rxt(k,406) + rxt(k,411)) * y(k,47))
         mat(k,610) = -rxt(k,245)*y(k,27)
         mat(k,552) = -rxt(k,247)*y(k,27)
         mat(k,798) = -rxt(k,248)*y(k,27)
         mat(k,391) = -(rxt(k,401) + rxt(k,406) + rxt(k,411)) * y(k,27)
         mat(k,470) = rxt(k,246)*y(k,60)
         mat(k,866) = rxt(k,246)*y(k,26)
         mat(k,154) = -((rxt(k,362) + rxt(k,366)) * y(k,109))
         mat(k,787) = -(rxt(k,362) + rxt(k,366)) * y(k,29)
         mat(k,575) = rxt(k,352)*y(k,61) + rxt(k,353)*y(k,66) + rxt(k,287)*y(k,97) &
                      + rxt(k,232)*y(k,98) + rxt(k,354)*y(k,109)
         mat(k,433) = rxt(k,352)*y(k,16)
         mat(k,536) = rxt(k,353)*y(k,16) + rxt(k,373)*y(k,70)
         mat(k,97) = rxt(k,373)*y(k,66) + rxt(k,374)*y(k,109)
         mat(k,367) = rxt(k,287)*y(k,16)
         mat(k,604) = rxt(k,232)*y(k,16)
         mat(k,787) = mat(k,787) + rxt(k,354)*y(k,16) + rxt(k,374)*y(k,70)
         mat(k,40) = -(rxt(k,326)*y(k,105))
         mat(k,668) = -rxt(k,326)*y(k,31)
         mat(k,56) = -(rxt(k,327)*y(k,105))
         mat(k,669) = -rxt(k,327)*y(k,32)
         mat(k,104) = -(rxt(k,370)*y(k,61) + (rxt(k,371) + rxt(k,388)) * y(k,109))
         mat(k,431) = -rxt(k,370)*y(k,33)
         mat(k,782) = -(rxt(k,371) + rxt(k,388)) * y(k,33)
         mat(k,220) = -(rxt(k,322)*y(k,39) + rxt(k,323)*y(k,112) + rxt(k,324)*y(k,49))
         mat(k,725) = -rxt(k,322)*y(k,37)
         mat(k,893) = -rxt(k,323)*y(k,37)
         mat(k,325) = -rxt(k,324)*y(k,37)
         mat(k,42) = 2.000_r8*rxt(k,326)*y(k,105)
         mat(k,58) = rxt(k,327)*y(k,105)
         mat(k,673) = 2.000_r8*rxt(k,326)*y(k,31) + rxt(k,327)*y(k,32)
         mat(k,357) = -((rxt(k,143) + rxt(k,144) + rxt(k,145)) * y(k,100) + rxt(k,146) &
                      *y(k,67) + rxt(k,151)*y(k,68))
         mat(k,504) = -(rxt(k,143) + rxt(k,144) + rxt(k,145)) * y(k,38)
         mat(k,651) = -rxt(k,146)*y(k,38)
         mat(k,705) = -rxt(k,151)*y(k,38)
         mat(k,581) = rxt(k,354)*y(k,109)
         mat(k,155) = rxt(k,366)*y(k,109)
         mat(k,222) = rxt(k,322)*y(k,39)
         mat(k,727) = rxt(k,322)*y(k,37) + rxt(k,139)*y(k,66) + rxt(k,234)*y(k,98) &
                      + rxt(k,111)*y(k,105) + rxt(k,153)*y(k,109)
         mat(k,297) = rxt(k,314)*y(k,105)
         mat(k,393) = rxt(k,272)*y(k,105)
         mat(k,277) = rxt(k,183)*y(k,109)
         mat(k,555) = rxt(k,139)*y(k,39) + rxt(k,156)*y(k,109)
         mat(k,101) = rxt(k,374)*y(k,109)
         mat(k,201) = rxt(k,379)*y(k,109)
         mat(k,342) = rxt(k,384)*y(k,109)
         mat(k,612) = rxt(k,234)*y(k,39)
         mat(k,676) = rxt(k,111)*y(k,39) + rxt(k,314)*y(k,43) + rxt(k,272)*y(k,47)
         mat(k,801) = rxt(k,354)*y(k,16) + rxt(k,366)*y(k,29) + rxt(k,153)*y(k,39) &
                      + rxt(k,183)*y(k,53) + rxt(k,156)*y(k,66) + rxt(k,374)*y(k,70) &
                      + rxt(k,379)*y(k,73) + rxt(k,384)*y(k,74)
         mat(k,739) = -(rxt(k,111)*y(k,105) + rxt(k,139)*y(k,66) + rxt(k,153)*y(k,109) &
                      + rxt(k,234)*y(k,98) + rxt(k,322)*y(k,37))
         mat(k,689) = -rxt(k,111)*y(k,39)
         mat(k,568) = -rxt(k,139)*y(k,39)
         mat(k,814) = -rxt(k,153)*y(k,39)
         mat(k,625) = -rxt(k,234)*y(k,39)
         mat(k,225) = -rxt(k,322)*y(k,39)
         mat(k,362) = rxt(k,143)*y(k,100)
         mat(k,517) = rxt(k,143)*y(k,38)
         mat(k,179) = -(rxt(k,140)*y(k,66) + rxt(k,154)*y(k,109) + rxt(k,235)*y(k,98))
         mat(k,540) = -rxt(k,140)*y(k,41)
         mat(k,791) = -rxt(k,154)*y(k,41)
         mat(k,608) = -rxt(k,235)*y(k,41)
         mat(k,500) = 2.000_r8*rxt(k,162)*y(k,100)
         mat(k,791) = mat(k,791) + 2.000_r8*rxt(k,159)*y(k,109)
         mat(k,69) = rxt(k,390)*y(k,112)
         mat(k,888) = rxt(k,390)*y(k,76)
         mat(k,296) = -(rxt(k,305)*y(k,66) + rxt(k,306)*y(k,109) + (rxt(k,313) &
                      + rxt(k,314)) * y(k,105))
         mat(k,550) = -rxt(k,305)*y(k,43)
         mat(k,796) = -rxt(k,306)*y(k,43)
         mat(k,674) = -(rxt(k,313) + rxt(k,314)) * y(k,43)
         mat(k,579) = rxt(k,287)*y(k,97)
         mat(k,368) = rxt(k,287)*y(k,16) + rxt(k,288)*y(k,100)
         mat(k,503) = rxt(k,288)*y(k,97)
         mat(k,395) = -(rxt(k,257)*y(k,66) + rxt(k,258)*y(k,109) + (rxt(k,271) &
                      + rxt(k,272)) * y(k,105) + (rxt(k,396) + rxt(k,402) + rxt(k,407) &
                      ) * y(k,52) + (rxt(k,401) + rxt(k,406) + rxt(k,411)) * y(k,27) &
                      + (rxt(k,403) + rxt(k,408)) * y(k,51))
         mat(k,557) = -rxt(k,257)*y(k,47)
         mat(k,803) = -rxt(k,258)*y(k,47)
         mat(k,678) = -(rxt(k,271) + rxt(k,272)) * y(k,47)
         mat(k,288) = -(rxt(k,396) + rxt(k,402) + rxt(k,407)) * y(k,47)
         mat(k,316) = -(rxt(k,401) + rxt(k,406) + rxt(k,411)) * y(k,47)
         mat(k,232) = -(rxt(k,403) + rxt(k,408)) * y(k,47)
         mat(k,583) = rxt(k,232)*y(k,98)
         mat(k,475) = rxt(k,252)*y(k,109)
         mat(k,728) = rxt(k,234)*y(k,98)
         mat(k,180) = rxt(k,235)*y(k,98)
         mat(k,288) = mat(k,288) + rxt(k,260)*y(k,98)
         mat(k,614) = rxt(k,232)*y(k,16) + rxt(k,234)*y(k,39) + rxt(k,235)*y(k,41) &
                      + rxt(k,260)*y(k,52) + rxt(k,236)*y(k,100)
         mat(k,506) = rxt(k,236)*y(k,98)
         mat(k,803) = mat(k,803) + rxt(k,252)*y(k,26)
         mat(k,216) = rxt(k,322)*y(k,39) + rxt(k,324)*y(k,49) + rxt(k,323)*y(k,112)
         mat(k,722) = rxt(k,322)*y(k,37)
         mat(k,324) = rxt(k,324)*y(k,37)
         mat(k,889) = rxt(k,323)*y(k,37)
         mat(k,326) = -(rxt(k,203)*y(k,109) + rxt(k,324)*y(k,37))
         mat(k,799) = -rxt(k,203)*y(k,49)
         mat(k,221) = -rxt(k,324)*y(k,49)
         mat(k,580) = rxt(k,352)*y(k,61)
         mat(k,315) = (rxt(k,401)+rxt(k,406)+rxt(k,411))*y(k,47)
         mat(k,108) = rxt(k,370)*y(k,61)
         mat(k,392) = (rxt(k,401)+rxt(k,406)+rxt(k,411))*y(k,27)
         mat(k,867) = rxt(k,197)*y(k,109)
         mat(k,438) = rxt(k,352)*y(k,16) + rxt(k,370)*y(k,33)
         mat(k,799) = mat(k,799) + rxt(k,197)*y(k,60)
         mat(k,126) = -(rxt(k,163)*y(k,109))
         mat(k,786) = -rxt(k,163)*y(k,50)
         mat(k,858) = rxt(k,195)*y(k,100)
         mat(k,497) = rxt(k,195)*y(k,60)
         mat(k,230) = -(rxt(k,308)*y(k,66) + (rxt(k,403) + rxt(k,408)) * y(k,47))
         mat(k,544) = -rxt(k,308)*y(k,51)
         mat(k,389) = -(rxt(k,403) + rxt(k,408)) * y(k,51)
         mat(k,829) = rxt(k,294)*y(k,100)
         mat(k,501) = rxt(k,294)*y(k,4)
         mat(k,287) = -(rxt(k,260)*y(k,98) + rxt(k,261)*y(k,66) + rxt(k,262)*y(k,109) &
                      + (rxt(k,396) + rxt(k,402) + rxt(k,407)) * y(k,47))
         mat(k,609) = -rxt(k,260)*y(k,52)
         mat(k,549) = -rxt(k,261)*y(k,52)
         mat(k,795) = -rxt(k,262)*y(k,52)
         mat(k,390) = -(rxt(k,396) + rxt(k,402) + rxt(k,407)) * y(k,52)
         mat(k,468) = rxt(k,243)*y(k,100)
         mat(k,313) = rxt(k,248)*y(k,109)
         mat(k,502) = rxt(k,243)*y(k,26)
         mat(k,795) = mat(k,795) + rxt(k,248)*y(k,27)
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
         mat(k,276) = -(rxt(k,166)*y(k,59) + (rxt(k,167) + rxt(k,168) + rxt(k,169) &
                      ) * y(k,60) + rxt(k,170)*y(k,67) + rxt(k,183)*y(k,109) + rxt(k,439) &
                      *y(k,108))
         mat(k,754) = -rxt(k,166)*y(k,53)
         mat(k,864) = -(rxt(k,167) + rxt(k,168) + rxt(k,169)) * y(k,53)
         mat(k,648) = -rxt(k,170)*y(k,53)
         mat(k,794) = -rxt(k,183)*y(k,53)
         mat(k,241) = -rxt(k,439)*y(k,53)
         mat(k,548) = rxt(k,164)*y(k,101) + rxt(k,436)*y(k,104)
         mat(k,648) = mat(k,648) + rxt(k,437)*y(k,104)
         mat(k,260) = 1.100_r8*rxt(k,432)*y(k,102) + .200_r8*rxt(k,430)*y(k,103)
         mat(k,141) = rxt(k,164)*y(k,66)
         mat(k,192) = 1.100_r8*rxt(k,432)*y(k,99)
         mat(k,249) = .200_r8*rxt(k,430)*y(k,99)
         mat(k,137) = rxt(k,436)*y(k,66) + rxt(k,437)*y(k,67)
         mat(k,854) = rxt(k,196)*y(k,61)
         mat(k,429) = rxt(k,196)*y(k,60)
         mat(k,768) = -(rxt(k,166)*y(k,53) + rxt(k,178)*y(k,61) + rxt(k,184)*y(k,100) &
                      + rxt(k,185)*y(k,68) + rxt(k,186)*y(k,66) + rxt(k,244)*y(k,26) &
                      + rxt(k,295)*y(k,4) + rxt(k,358)*y(k,20) + rxt(k,441)*y(k,108))
         mat(k,282) = -rxt(k,166)*y(k,59)
         mat(k,453) = -rxt(k,178)*y(k,59)
         mat(k,518) = -rxt(k,184)*y(k,59)
         mat(k,717) = -rxt(k,185)*y(k,59)
         mat(k,569) = -rxt(k,186)*y(k,59)
         mat(k,487) = -rxt(k,244)*y(k,59)
         mat(k,847) = -rxt(k,295)*y(k,59)
         mat(k,424) = -rxt(k,358)*y(k,59)
         mat(k,245) = -rxt(k,441)*y(k,59)
         mat(k,282) = mat(k,282) + 2.000_r8*rxt(k,168)*y(k,60) + rxt(k,170)*y(k,67) &
                      + rxt(k,183)*y(k,109)
         mat(k,883) = 2.000_r8*rxt(k,168)*y(k,53) + rxt(k,171)*y(k,66) + rxt(k,380) &
                      *y(k,74)
         mat(k,569) = mat(k,569) + rxt(k,171)*y(k,60)
         mat(k,662) = rxt(k,170)*y(k,53) + rxt(k,165)*y(k,101)
         mat(k,350) = rxt(k,380)*y(k,60)
         mat(k,144) = rxt(k,165)*y(k,67)
         mat(k,815) = rxt(k,183)*y(k,53)
         mat(k,886) = -((rxt(k,167) + rxt(k,168) + rxt(k,169)) * y(k,53) + (rxt(k,171) &
                      + rxt(k,173)) * y(k,66) + rxt(k,172)*y(k,68) + rxt(k,195) &
                      *y(k,100) + rxt(k,196)*y(k,61) + rxt(k,197)*y(k,109) + rxt(k,246) &
                      *y(k,26) + rxt(k,296)*y(k,4) + rxt(k,380)*y(k,74))
         mat(k,284) = -(rxt(k,167) + rxt(k,168) + rxt(k,169)) * y(k,60)
         mat(k,572) = -(rxt(k,171) + rxt(k,173)) * y(k,60)
         mat(k,720) = -rxt(k,172)*y(k,60)
         mat(k,521) = -rxt(k,195)*y(k,60)
         mat(k,456) = -rxt(k,196)*y(k,60)
         mat(k,818) = -rxt(k,197)*y(k,60)
         mat(k,490) = -rxt(k,246)*y(k,60)
         mat(k,850) = -rxt(k,296)*y(k,60)
         mat(k,353) = -rxt(k,380)*y(k,60)
         mat(k,850) = mat(k,850) + rxt(k,295)*y(k,59)
         mat(k,426) = rxt(k,358)*y(k,59)
         mat(k,490) = mat(k,490) + rxt(k,244)*y(k,59)
         mat(k,131) = rxt(k,163)*y(k,109)
         mat(k,771) = rxt(k,295)*y(k,4) + rxt(k,358)*y(k,20) + rxt(k,244)*y(k,26) &
                      + 2.000_r8*rxt(k,178)*y(k,61) + rxt(k,186)*y(k,66) + rxt(k,185) &
                      *y(k,68) + rxt(k,184)*y(k,100)
         mat(k,456) = mat(k,456) + 2.000_r8*rxt(k,178)*y(k,59) + rxt(k,179)*y(k,66) &
                      + rxt(k,177)*y(k,100) + rxt(k,180)*y(k,109)
         mat(k,572) = mat(k,572) + rxt(k,186)*y(k,59) + rxt(k,179)*y(k,61)
         mat(k,720) = mat(k,720) + rxt(k,185)*y(k,59)
         mat(k,521) = mat(k,521) + rxt(k,184)*y(k,59) + rxt(k,177)*y(k,61)
         mat(k,818) = mat(k,818) + rxt(k,163)*y(k,50) + rxt(k,180)*y(k,61)
         mat(k,443) = -(rxt(k,177)*y(k,100) + rxt(k,178)*y(k,59) + rxt(k,179)*y(k,66) &
                      + rxt(k,180)*y(k,109) + rxt(k,196)*y(k,60) + rxt(k,352)*y(k,16) &
                      + rxt(k,370)*y(k,33))
         mat(k,508) = -rxt(k,177)*y(k,61)
         mat(k,758) = -rxt(k,178)*y(k,61)
         mat(k,559) = -rxt(k,179)*y(k,61)
         mat(k,805) = -rxt(k,180)*y(k,61)
         mat(k,873) = -rxt(k,196)*y(k,61)
         mat(k,585) = -rxt(k,352)*y(k,61)
         mat(k,109) = -rxt(k,370)*y(k,61)
         mat(k,150) = rxt(k,297)*y(k,66)
         mat(k,317) = rxt(k,247)*y(k,66) + rxt(k,245)*y(k,98) + rxt(k,248)*y(k,109)
         mat(k,224) = rxt(k,324)*y(k,49)
         mat(k,329) = rxt(k,324)*y(k,37) + rxt(k,203)*y(k,109)
         mat(k,873) = mat(k,873) + rxt(k,173)*y(k,66) + rxt(k,172)*y(k,68)
         mat(k,559) = mat(k,559) + rxt(k,297)*y(k,5) + rxt(k,247)*y(k,27) + rxt(k,173) &
                      *y(k,60)
         mat(k,707) = rxt(k,172)*y(k,60)
         mat(k,616) = rxt(k,245)*y(k,27)
         mat(k,805) = mat(k,805) + rxt(k,248)*y(k,27) + rxt(k,203)*y(k,49)
         mat(k,562) = -(rxt(k,126)*y(k,68) + 4._r8*rxt(k,128)*y(k,66) + rxt(k,129) &
                      *y(k,67) + rxt(k,139)*y(k,39) + rxt(k,140)*y(k,41) + rxt(k,147) &
                      *y(k,100) + rxt(k,156)*y(k,109) + (rxt(k,171) + rxt(k,173) &
                      ) * y(k,60) + rxt(k,179)*y(k,61) + rxt(k,186)*y(k,59) + rxt(k,247) &
                      *y(k,27) + rxt(k,250)*y(k,26) + rxt(k,257)*y(k,47) + rxt(k,261) &
                      *y(k,52) + rxt(k,297)*y(k,5) + rxt(k,299)*y(k,4) + rxt(k,305) &
                      *y(k,43) + rxt(k,308)*y(k,51) + rxt(k,353)*y(k,16) + rxt(k,373) &
                      *y(k,70) + (rxt(k,434) + rxt(k,435)) * y(k,102) + rxt(k,436) &
                      *y(k,104))
         mat(k,710) = -rxt(k,126)*y(k,66)
         mat(k,655) = -rxt(k,129)*y(k,66)
         mat(k,733) = -rxt(k,139)*y(k,66)
         mat(k,182) = -rxt(k,140)*y(k,66)
         mat(k,511) = -rxt(k,147)*y(k,66)
         mat(k,808) = -rxt(k,156)*y(k,66)
         mat(k,876) = -(rxt(k,171) + rxt(k,173)) * y(k,66)
         mat(k,446) = -rxt(k,179)*y(k,66)
         mat(k,761) = -rxt(k,186)*y(k,66)
         mat(k,319) = -rxt(k,247)*y(k,66)
         mat(k,480) = -rxt(k,250)*y(k,66)
         mat(k,400) = -rxt(k,257)*y(k,66)
         mat(k,290) = -rxt(k,261)*y(k,66)
         mat(k,151) = -rxt(k,297)*y(k,66)
         mat(k,840) = -rxt(k,299)*y(k,66)
         mat(k,299) = -rxt(k,305)*y(k,66)
         mat(k,233) = -rxt(k,308)*y(k,66)
         mat(k,588) = -rxt(k,353)*y(k,66)
         mat(k,102) = -rxt(k,373)*y(k,66)
         mat(k,193) = -(rxt(k,434) + rxt(k,435)) * y(k,66)
         mat(k,138) = -rxt(k,436)*y(k,66)
         mat(k,359) = rxt(k,145)*y(k,100)
         mat(k,279) = rxt(k,166)*y(k,59) + rxt(k,167)*y(k,60) + rxt(k,170)*y(k,67) &
                      + rxt(k,439)*y(k,108)
         mat(k,761) = mat(k,761) + rxt(k,166)*y(k,53)
         mat(k,876) = mat(k,876) + rxt(k,167)*y(k,53)
         mat(k,655) = mat(k,655) + rxt(k,170)*y(k,53) + rxt(k,375)*y(k,73) &
                      + rxt(k,381)*y(k,74) + rxt(k,438)*y(k,104) + (rxt(k,114) &
                       +rxt(k,115))*y(k,105) + rxt(k,444)*y(k,110)
         mat(k,202) = rxt(k,375)*y(k,67)
         mat(k,346) = rxt(k,381)*y(k,67)
         mat(k,263) = rxt(k,430)*y(k,103) + 1.150_r8*rxt(k,431)*y(k,108)
         mat(k,511) = mat(k,511) + rxt(k,145)*y(k,38)
         mat(k,250) = rxt(k,430)*y(k,99)
         mat(k,138) = mat(k,138) + rxt(k,438)*y(k,67)
         mat(k,683) = (rxt(k,114)+rxt(k,115))*y(k,67)
         mat(k,242) = rxt(k,439)*y(k,53) + 1.150_r8*rxt(k,431)*y(k,99)
         mat(k,808) = mat(k,808) + 2.000_r8*rxt(k,158)*y(k,109)
         mat(k,213) = rxt(k,444)*y(k,67)
         mat(k,658) = -(rxt(k,114)*y(k,105) + rxt(k,120)*y(k,106) + rxt(k,129)*y(k,66) &
                      + rxt(k,146)*y(k,38) + rxt(k,165)*y(k,101) + rxt(k,170)*y(k,53) &
                      + rxt(k,375)*y(k,73) + rxt(k,381)*y(k,74) + rxt(k,433)*y(k,102) &
                      + (rxt(k,437) + rxt(k,438)) * y(k,104) + rxt(k,444)*y(k,110))
         mat(k,686) = -rxt(k,114)*y(k,67)
         mat(k,34) = -rxt(k,120)*y(k,67)
         mat(k,565) = -rxt(k,129)*y(k,67)
         mat(k,360) = -rxt(k,146)*y(k,67)
         mat(k,142) = -rxt(k,165)*y(k,67)
         mat(k,280) = -rxt(k,170)*y(k,67)
         mat(k,203) = -rxt(k,375)*y(k,67)
         mat(k,348) = -rxt(k,381)*y(k,67)
         mat(k,194) = -rxt(k,433)*y(k,67)
         mat(k,139) = -(rxt(k,437) + rxt(k,438)) * y(k,67)
         mat(k,214) = -rxt(k,444)*y(k,67)
         mat(k,843) = 2.000_r8*rxt(k,290)*y(k,4) + (rxt(k,292)+rxt(k,293))*y(k,26) &
                      + rxt(k,299)*y(k,66) + rxt(k,294)*y(k,100)
         mat(k,421) = rxt(k,357)*y(k,100)
         mat(k,483) = (rxt(k,292)+rxt(k,293))*y(k,4) + (2.000_r8*rxt(k,240) &
                       +2.000_r8*rxt(k,241))*y(k,26) + rxt(k,250)*y(k,66) + rxt(k,243) &
                      *y(k,100) + rxt(k,252)*y(k,109)
         mat(k,360) = mat(k,360) + rxt(k,151)*y(k,68) + rxt(k,143)*y(k,100)
         mat(k,129) = rxt(k,163)*y(k,109)
         mat(k,280) = mat(k,280) + rxt(k,169)*y(k,60)
         mat(k,764) = rxt(k,185)*y(k,68) + rxt(k,441)*y(k,108)
         mat(k,879) = rxt(k,169)*y(k,53) + rxt(k,171)*y(k,66) + rxt(k,172)*y(k,68)
         mat(k,449) = rxt(k,179)*y(k,66) + rxt(k,177)*y(k,100)
         mat(k,565) = mat(k,565) + rxt(k,299)*y(k,4) + rxt(k,250)*y(k,26) + rxt(k,171) &
                      *y(k,60) + rxt(k,179)*y(k,61) + 2.000_r8*rxt(k,128)*y(k,66) &
                      + 2.000_r8*rxt(k,126)*y(k,68) + rxt(k,147)*y(k,100) + rxt(k,119) &
                      *y(k,106) + rxt(k,156)*y(k,109)
         mat(k,658) = mat(k,658) + 2.000_r8*rxt(k,120)*y(k,106)
         mat(k,713) = rxt(k,151)*y(k,38) + rxt(k,185)*y(k,59) + rxt(k,172)*y(k,60) &
                      + 2.000_r8*rxt(k,126)*y(k,66) + rxt(k,376)*y(k,73) + rxt(k,382) &
                      *y(k,74) + rxt(k,289)*y(k,97) + rxt(k,238)*y(k,98) &
                      + 2.000_r8*rxt(k,148)*y(k,100) + 2.000_r8*rxt(k,116)*y(k,105) &
                      + rxt(k,157)*y(k,109)
         mat(k,203) = mat(k,203) + rxt(k,376)*y(k,68)
         mat(k,348) = mat(k,348) + rxt(k,382)*y(k,68)
         mat(k,374) = rxt(k,289)*y(k,68) + rxt(k,288)*y(k,100)
         mat(k,622) = rxt(k,238)*y(k,68) + rxt(k,236)*y(k,100)
         mat(k,514) = rxt(k,294)*y(k,4) + rxt(k,357)*y(k,20) + rxt(k,243)*y(k,26) &
                      + rxt(k,143)*y(k,38) + rxt(k,177)*y(k,61) + rxt(k,147)*y(k,66) &
                      + 2.000_r8*rxt(k,148)*y(k,68) + rxt(k,288)*y(k,97) + rxt(k,236) &
                      *y(k,98) + 2.000_r8*rxt(k,162)*y(k,100) + rxt(k,155)*y(k,109)
         mat(k,686) = mat(k,686) + 2.000_r8*rxt(k,116)*y(k,68)
         mat(k,34) = mat(k,34) + rxt(k,119)*y(k,66) + 2.000_r8*rxt(k,120)*y(k,67)
         mat(k,243) = rxt(k,441)*y(k,59)
         mat(k,811) = rxt(k,252)*y(k,26) + rxt(k,163)*y(k,50) + rxt(k,156)*y(k,66) &
                      + rxt(k,157)*y(k,68) + rxt(k,155)*y(k,100)
         mat(k,715) = -(rxt(k,116)*y(k,105) + rxt(k,126)*y(k,66) + rxt(k,148)*y(k,100) &
                      + rxt(k,151)*y(k,38) + rxt(k,157)*y(k,109) + rxt(k,172)*y(k,60) &
                      + rxt(k,185)*y(k,59) + rxt(k,238)*y(k,98) + rxt(k,289)*y(k,97) &
                      + rxt(k,376)*y(k,73) + rxt(k,382)*y(k,74))
         mat(k,688) = -rxt(k,116)*y(k,68)
         mat(k,567) = -rxt(k,126)*y(k,68)
         mat(k,516) = -rxt(k,148)*y(k,68)
         mat(k,361) = -rxt(k,151)*y(k,68)
         mat(k,813) = -rxt(k,157)*y(k,68)
         mat(k,881) = -rxt(k,172)*y(k,68)
         mat(k,766) = -rxt(k,185)*y(k,68)
         mat(k,624) = -rxt(k,238)*y(k,68)
         mat(k,376) = -rxt(k,289)*y(k,68)
         mat(k,204) = -rxt(k,376)*y(k,68)
         mat(k,349) = -rxt(k,382)*y(k,68)
         mat(k,567) = mat(k,567) + rxt(k,129)*y(k,67)
         mat(k,660) = rxt(k,129)*y(k,66)
         mat(k,86) = -(rxt(k,383)*y(k,74))
         mat(k,334) = -rxt(k,383)*y(k,69)
         mat(k,822) = rxt(k,291)*y(k,26)
         mat(k,461) = rxt(k,291)*y(k,4) + 2.000_r8*rxt(k,242)*y(k,26)
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
         mat(k,96) = -(rxt(k,373)*y(k,66) + rxt(k,374)*y(k,109))
         mat(k,530) = -rxt(k,373)*y(k,70)
         mat(k,781) = -rxt(k,374)*y(k,70)
         mat(k,199) = -(rxt(k,375)*y(k,67) + rxt(k,376)*y(k,68) + rxt(k,379)*y(k,109))
         mat(k,643) = -rxt(k,375)*y(k,73)
         mat(k,702) = -rxt(k,376)*y(k,73)
         mat(k,792) = -rxt(k,379)*y(k,73)
         mat(k,341) = -(rxt(k,377)*y(k,4) + rxt(k,378)*y(k,26) + rxt(k,380)*y(k,60) &
                      + rxt(k,381)*y(k,67) + rxt(k,382)*y(k,68) + rxt(k,383)*y(k,69) &
                      + rxt(k,384)*y(k,109))
         mat(k,832) = -rxt(k,377)*y(k,74)
         mat(k,472) = -rxt(k,378)*y(k,74)
         mat(k,868) = -rxt(k,380)*y(k,74)
         mat(k,650) = -rxt(k,381)*y(k,74)
         mat(k,704) = -rxt(k,382)*y(k,74)
         mat(k,88) = -rxt(k,383)*y(k,74)
         mat(k,800) = -rxt(k,384)*y(k,74)
         mat(k,554) = rxt(k,373)*y(k,70)
         mat(k,650) = mat(k,650) + rxt(k,375)*y(k,73)
         mat(k,704) = mat(k,704) + rxt(k,376)*y(k,73)
         mat(k,100) = rxt(k,373)*y(k,66)
         mat(k,200) = rxt(k,375)*y(k,67) + rxt(k,376)*y(k,68) + rxt(k,379)*y(k,109)
         mat(k,800) = mat(k,800) + rxt(k,379)*y(k,73)
         mat(k,305) = -(rxt(k,389)*y(k,109))
         mat(k,797) = -rxt(k,389)*y(k,75)
         mat(k,830) = rxt(k,377)*y(k,74)
         mat(k,469) = rxt(k,378)*y(k,74)
         mat(k,107) = rxt(k,370)*y(k,61) + (rxt(k,371)+.500_r8*rxt(k,388))*y(k,109)
         mat(k,865) = rxt(k,380)*y(k,74)
         mat(k,437) = rxt(k,370)*y(k,33)
         mat(k,649) = rxt(k,381)*y(k,74)
         mat(k,703) = rxt(k,382)*y(k,74)
         mat(k,87) = rxt(k,383)*y(k,74)
         mat(k,99) = rxt(k,374)*y(k,109)
         mat(k,340) = rxt(k,377)*y(k,4) + rxt(k,378)*y(k,26) + rxt(k,380)*y(k,60) &
                      + rxt(k,381)*y(k,67) + rxt(k,382)*y(k,68) + rxt(k,383)*y(k,69) &
                      + rxt(k,384)*y(k,109)
         mat(k,797) = mat(k,797) + (rxt(k,371)+.500_r8*rxt(k,388))*y(k,33) &
                      + rxt(k,374)*y(k,70) + rxt(k,384)*y(k,74)
         mat(k,70) = -(rxt(k,390)*y(k,112))
         mat(k,890) = -rxt(k,390)*y(k,76)
         mat(k,304) = rxt(k,389)*y(k,109)
         mat(k,777) = rxt(k,389)*y(k,75)
         mat(k,74) = -(rxt(k,315)*y(k,66))
         mat(k,527) = -rxt(k,315)*y(k,83)
         mat(k,821) = rxt(k,302)*y(k,90)
         mat(k,171) = rxt(k,302)*y(k,4)
         mat(k,122) = -(rxt(k,275)*y(k,98) + rxt(k,276)*y(k,66) + rxt(k,277)*y(k,109) &
                      + (rxt(k,418) + rxt(k,424) + rxt(k,429)) * y(k,47))
         mat(k,603) = -rxt(k,275)*y(k,84)
         mat(k,532) = -rxt(k,276)*y(k,84)
         mat(k,785) = -rxt(k,277)*y(k,84)
         mat(k,385) = -(rxt(k,418) + rxt(k,424) + rxt(k,429)) * y(k,84)
         mat(k,464) = rxt(k,254)*y(k,90)
         mat(k,174) = rxt(k,254)*y(k,26)
         mat(k,158) = -(rxt(k,206)*y(k,109) + rxt(k,325)*y(k,37))
         mat(k,788) = -rxt(k,206)*y(k,85)
         mat(k,217) = -rxt(k,325)*y(k,85)
         mat(k,576) = rxt(k,355)*y(k,93)
         mat(k,105) = rxt(k,372)*y(k,93)
         mat(k,386) = (rxt(k,418)+rxt(k,424)+rxt(k,429))*y(k,84)
         mat(k,123) = (rxt(k,418)+rxt(k,424)+rxt(k,429))*y(k,47)
         mat(k,175) = rxt(k,201)*y(k,109)
         mat(k,165) = rxt(k,355)*y(k,16) + rxt(k,372)*y(k,33)
         mat(k,788) = mat(k,788) + rxt(k,201)*y(k,90)
         mat(k,44) = -(rxt(k,211)*y(k,109))
         mat(k,774) = -rxt(k,211)*y(k,86)
         mat(k,168) = rxt(k,199)*y(k,100)
         mat(k,492) = rxt(k,199)*y(k,90)
         mat(k,47) = -(rxt(k,214)*y(k,59) + (rxt(k,215) + rxt(k,216) + rxt(k,217) &
                      ) * y(k,60) + rxt(k,218)*y(k,67) + rxt(k,226)*y(k,109))
         mat(k,745) = -rxt(k,214)*y(k,87)
         mat(k,853) = -(rxt(k,215) + rxt(k,216) + rxt(k,217)) * y(k,87)
         mat(k,634) = -rxt(k,218)*y(k,87)
         mat(k,775) = -rxt(k,226)*y(k,87)
         mat(k,526) = rxt(k,212)*y(k,88)
         mat(k,28) = rxt(k,212)*y(k,66)
         mat(k,27) = -(rxt(k,212)*y(k,66) + rxt(k,213)*y(k,67))
         mat(k,523) = -rxt(k,212)*y(k,88)
         mat(k,631) = -rxt(k,213)*y(k,88)
         mat(k,120) = -(rxt(k,181)*y(k,61) + rxt(k,189)*y(k,53) + rxt(k,227)*y(k,100) &
                      + rxt(k,228)*y(k,68) + rxt(k,229)*y(k,66) + rxt(k,253)*y(k,26) &
                      + rxt(k,301)*y(k,4) + rxt(k,359)*y(k,20))
         mat(k,432) = -rxt(k,181)*y(k,89)
         mat(k,269) = -rxt(k,189)*y(k,89)
         mat(k,496) = -rxt(k,227)*y(k,89)
         mat(k,699) = -rxt(k,228)*y(k,89)
         mat(k,531) = -rxt(k,229)*y(k,89)
         mat(k,463) = -rxt(k,253)*y(k,89)
         mat(k,824) = -rxt(k,301)*y(k,89)
         mat(k,412) = -rxt(k,359)*y(k,89)
         mat(k,269) = mat(k,269) + rxt(k,191)*y(k,90)
         mat(k,857) = rxt(k,216)*y(k,87)
         mat(k,531) = mat(k,531) + rxt(k,219)*y(k,90)
         mat(k,638) = rxt(k,218)*y(k,87) + rxt(k,213)*y(k,88)
         mat(k,337) = rxt(k,385)*y(k,90)
         mat(k,48) = rxt(k,216)*y(k,60) + rxt(k,218)*y(k,67) + rxt(k,226)*y(k,109)
         mat(k,29) = rxt(k,213)*y(k,67)
         mat(k,173) = rxt(k,191)*y(k,53) + rxt(k,219)*y(k,66) + rxt(k,385)*y(k,74)
         mat(k,784) = rxt(k,226)*y(k,87)
         mat(k,177) = -((rxt(k,190) + rxt(k,191) + rxt(k,192)) * y(k,53) + rxt(k,199) &
                      *y(k,100) + rxt(k,200)*y(k,61) + rxt(k,201)*y(k,109) + rxt(k,202) &
                      *y(k,93) + (rxt(k,219) + rxt(k,221)) * y(k,66) + rxt(k,220) &
                      *y(k,68) + rxt(k,254)*y(k,26) + rxt(k,302)*y(k,4) + rxt(k,385) &
                      *y(k,74))
         mat(k,271) = -(rxt(k,190) + rxt(k,191) + rxt(k,192)) * y(k,90)
         mat(k,499) = -rxt(k,199)*y(k,90)
         mat(k,436) = -rxt(k,200)*y(k,90)
         mat(k,790) = -rxt(k,201)*y(k,90)
         mat(k,167) = -rxt(k,202)*y(k,90)
         mat(k,539) = -(rxt(k,219) + rxt(k,221)) * y(k,90)
         mat(k,701) = -rxt(k,220)*y(k,90)
         mat(k,467) = -rxt(k,254)*y(k,90)
         mat(k,828) = -rxt(k,302)*y(k,90)
         mat(k,338) = -rxt(k,385)*y(k,90)
         mat(k,828) = mat(k,828) + rxt(k,301)*y(k,89)
         mat(k,413) = rxt(k,359)*y(k,89)
         mat(k,467) = mat(k,467) + rxt(k,253)*y(k,89)
         mat(k,750) = rxt(k,223)*y(k,93)
         mat(k,436) = mat(k,436) + rxt(k,181)*y(k,89)
         mat(k,539) = mat(k,539) + rxt(k,229)*y(k,89) + rxt(k,224)*y(k,93)
         mat(k,701) = mat(k,701) + rxt(k,228)*y(k,89)
         mat(k,46) = rxt(k,211)*y(k,109)
         mat(k,121) = rxt(k,301)*y(k,4) + rxt(k,359)*y(k,20) + rxt(k,253)*y(k,26) &
                      + rxt(k,181)*y(k,61) + rxt(k,229)*y(k,66) + rxt(k,228)*y(k,68) &
                      + rxt(k,227)*y(k,100)
         mat(k,167) = mat(k,167) + rxt(k,223)*y(k,59) + rxt(k,224)*y(k,66) &
                      + rxt(k,222)*y(k,100) + rxt(k,225)*y(k,109)
         mat(k,499) = mat(k,499) + rxt(k,227)*y(k,89) + rxt(k,222)*y(k,93)
         mat(k,790) = mat(k,790) + rxt(k,211)*y(k,86) + rxt(k,225)*y(k,93)
         mat(k,428) = rxt(k,200)*y(k,90)
         mat(k,170) = rxt(k,200)*y(k,61)
         mat(k,169) = rxt(k,202)*y(k,93)
         mat(k,162) = rxt(k,202)*y(k,90)
         mat(k,166) = -(rxt(k,198)*y(k,60) + rxt(k,202)*y(k,90) + rxt(k,222)*y(k,100) &
                      + rxt(k,223)*y(k,59) + rxt(k,224)*y(k,66) + rxt(k,225)*y(k,109) &
                      + rxt(k,355)*y(k,16) + rxt(k,372)*y(k,33))
         mat(k,861) = -rxt(k,198)*y(k,93)
         mat(k,176) = -rxt(k,202)*y(k,93)
         mat(k,498) = -rxt(k,222)*y(k,93)
         mat(k,749) = -rxt(k,223)*y(k,93)
         mat(k,538) = -rxt(k,224)*y(k,93)
         mat(k,789) = -rxt(k,225)*y(k,93)
         mat(k,577) = -rxt(k,355)*y(k,93)
         mat(k,106) = -rxt(k,372)*y(k,93)
         mat(k,218) = rxt(k,325)*y(k,85)
         mat(k,538) = mat(k,538) + rxt(k,315)*y(k,83) + rxt(k,276)*y(k,84) &
                      + rxt(k,221)*y(k,90)
         mat(k,700) = rxt(k,220)*y(k,90)
         mat(k,76) = rxt(k,315)*y(k,66)
         mat(k,124) = rxt(k,276)*y(k,66) + rxt(k,275)*y(k,98) + rxt(k,277)*y(k,109)
         mat(k,159) = rxt(k,325)*y(k,37) + rxt(k,206)*y(k,109)
         mat(k,176) = mat(k,176) + rxt(k,221)*y(k,66) + rxt(k,220)*y(k,68)
         mat(k,606) = rxt(k,275)*y(k,84)
         mat(k,789) = mat(k,789) + rxt(k,277)*y(k,84) + rxt(k,206)*y(k,85)
         mat(k,852) = rxt(k,198)*y(k,93)
         mat(k,161) = rxt(k,198)*y(k,60)
         mat(k,95) = -(rxt(k,130)*y(k,66) + rxt(k,131)*y(k,67) + rxt(k,138)*y(k,68) &
                      + rxt(k,141)*y(k,41) + rxt(k,142)*y(k,39) + rxt(k,149)*y(k,100) &
                      + rxt(k,160)*y(k,109) + (rxt(k,174) + rxt(k,176)) * y(k,60) &
                      + rxt(k,182)*y(k,61) + rxt(k,188)*y(k,59) + rxt(k,249)*y(k,27) &
                      + rxt(k,255)*y(k,26) + rxt(k,259)*y(k,47) + rxt(k,263)*y(k,52) &
                      + rxt(k,298)*y(k,5) + rxt(k,303)*y(k,4) + rxt(k,307)*y(k,43) &
                      + rxt(k,309)*y(k,51) + rxt(k,356)*y(k,16))
         mat(k,529) = -rxt(k,130)*y(k,95)
         mat(k,637) = -rxt(k,131)*y(k,95)
         mat(k,698) = -rxt(k,138)*y(k,95)
         mat(k,178) = -rxt(k,141)*y(k,95)
         mat(k,724) = -rxt(k,142)*y(k,95)
         mat(k,494) = -rxt(k,149)*y(k,95)
         mat(k,780) = -rxt(k,160)*y(k,95)
         mat(k,856) = -(rxt(k,174) + rxt(k,176)) * y(k,95)
         mat(k,430) = -rxt(k,182)*y(k,95)
         mat(k,747) = -rxt(k,188)*y(k,95)
         mat(k,312) = -rxt(k,249)*y(k,95)
         mat(k,462) = -rxt(k,255)*y(k,95)
         mat(k,384) = -rxt(k,259)*y(k,95)
         mat(k,286) = -rxt(k,263)*y(k,95)
         mat(k,145) = -rxt(k,298)*y(k,95)
         mat(k,823) = -rxt(k,303)*y(k,95)
         mat(k,295) = -rxt(k,307)*y(k,95)
         mat(k,229) = -rxt(k,309)*y(k,95)
         mat(k,574) = -rxt(k,356)*y(k,95)
         mat(k,268) = rxt(k,190)*y(k,90)
         mat(k,637) = mat(k,637) + (rxt(k,135)+rxt(k,136))*y(k,111)
         mat(k,172) = rxt(k,190)*y(k,53)
         mat(k,79) = (rxt(k,135)+rxt(k,136))*y(k,67)
         mat(k,92) = -(rxt(k,117)*y(k,105) + rxt(k,127)*y(k,66) + rxt(k,150)*y(k,100) &
                      + rxt(k,152)*y(k,38) + rxt(k,161)*y(k,109) + rxt(k,175)*y(k,60) &
                      + rxt(k,187)*y(k,59) + rxt(k,256)*y(k,98) + rxt(k,304)*y(k,97) &
                      + rxt(k,386)*y(k,74) + rxt(k,387)*y(k,73))
         mat(k,671) = -rxt(k,117)*y(k,96)
         mat(k,528) = -rxt(k,127)*y(k,96)
         mat(k,493) = -rxt(k,150)*y(k,96)
         mat(k,355) = -rxt(k,152)*y(k,96)
         mat(k,779) = -rxt(k,161)*y(k,96)
         mat(k,855) = -rxt(k,175)*y(k,96)
         mat(k,746) = -rxt(k,187)*y(k,96)
         mat(k,601) = -rxt(k,256)*y(k,96)
         mat(k,365) = -rxt(k,304)*y(k,96)
         mat(k,335) = -rxt(k,386)*y(k,96)
         mat(k,197) = -rxt(k,387)*y(k,96)
         mat(k,636) = rxt(k,131)*y(k,95)
         mat(k,94) = rxt(k,131)*y(k,67)
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
         mat(k,370) = -(rxt(k,287)*y(k,16) + rxt(k,288)*y(k,100) + rxt(k,289)*y(k,68))
         mat(k,582) = -rxt(k,287)*y(k,97)
         mat(k,505) = -rxt(k,288)*y(k,97)
         mat(k,706) = -rxt(k,289)*y(k,97)
         mat(k,834) = 4.000_r8*rxt(k,290)*y(k,4) + (rxt(k,291)+rxt(k,292))*y(k,26) &
                      + rxt(k,295)*y(k,59) + rxt(k,299)*y(k,66) + rxt(k,377)*y(k,74) &
                      + rxt(k,300)*y(k,109)
         mat(k,474) = (rxt(k,291)+rxt(k,292))*y(k,4)
         mat(k,298) = rxt(k,305)*y(k,66) + rxt(k,313)*y(k,105) + rxt(k,306)*y(k,109)
         mat(k,756) = rxt(k,295)*y(k,4)
         mat(k,556) = rxt(k,299)*y(k,4) + rxt(k,305)*y(k,43)
         mat(k,343) = rxt(k,377)*y(k,4)
         mat(k,677) = rxt(k,313)*y(k,43)
         mat(k,802) = rxt(k,300)*y(k,4) + rxt(k,306)*y(k,43)
         mat(k,621) = -(rxt(k,232)*y(k,16) + rxt(k,234)*y(k,39) + rxt(k,235)*y(k,41) &
                      + (rxt(k,236) + rxt(k,237)) * y(k,100) + rxt(k,238)*y(k,68) &
                      + rxt(k,245)*y(k,27) + rxt(k,260)*y(k,52))
         mat(k,590) = -rxt(k,232)*y(k,98)
         mat(k,735) = -rxt(k,234)*y(k,98)
         mat(k,183) = -rxt(k,235)*y(k,98)
         mat(k,513) = -(rxt(k,236) + rxt(k,237)) * y(k,98)
         mat(k,712) = -rxt(k,238)*y(k,98)
         mat(k,320) = -rxt(k,245)*y(k,98)
         mat(k,291) = -rxt(k,260)*y(k,98)
         mat(k,842) = rxt(k,292)*y(k,26)
         mat(k,420) = rxt(k,239)*y(k,26)
         mat(k,482) = rxt(k,292)*y(k,4) + rxt(k,239)*y(k,20) + (4.000_r8*rxt(k,240) &
                       +2.000_r8*rxt(k,242))*y(k,26) + rxt(k,244)*y(k,59) + rxt(k,250) &
                      *y(k,66) + rxt(k,378)*y(k,74) + rxt(k,251)*y(k,109)
         mat(k,59) = rxt(k,327)*y(k,105)
         mat(k,402) = rxt(k,257)*y(k,66) + rxt(k,271)*y(k,105) + rxt(k,258)*y(k,109)
         mat(k,763) = rxt(k,244)*y(k,26)
         mat(k,564) = rxt(k,250)*y(k,26) + rxt(k,257)*y(k,47)
         mat(k,347) = rxt(k,378)*y(k,26)
         mat(k,685) = rxt(k,327)*y(k,32) + rxt(k,271)*y(k,47)
         mat(k,810) = rxt(k,251)*y(k,26) + rxt(k,258)*y(k,47)
         mat(k,259) = -(rxt(k,430)*y(k,103) + rxt(k,431)*y(k,108) + rxt(k,432) &
                      *y(k,102))
         mat(k,248) = -rxt(k,430)*y(k,99)
         mat(k,240) = -rxt(k,431)*y(k,99)
         mat(k,191) = -rxt(k,432)*y(k,99)
         mat(k,510) = -((rxt(k,143) + rxt(k,144) + rxt(k,145)) * y(k,38) + rxt(k,147) &
                      *y(k,66) + rxt(k,148)*y(k,68) + rxt(k,155)*y(k,109) &
                      + 4._r8*rxt(k,162)*y(k,100) + rxt(k,177)*y(k,61) + rxt(k,184) &
                      *y(k,59) + rxt(k,195)*y(k,60) + (rxt(k,236) + rxt(k,237) &
                      ) * y(k,98) + rxt(k,243)*y(k,26) + rxt(k,288)*y(k,97) + rxt(k,294) &
                      *y(k,4) + rxt(k,357)*y(k,20))
         mat(k,358) = -(rxt(k,143) + rxt(k,144) + rxt(k,145)) * y(k,100)
         mat(k,561) = -rxt(k,147)*y(k,100)
         mat(k,709) = -rxt(k,148)*y(k,100)
         mat(k,807) = -rxt(k,155)*y(k,100)
         mat(k,445) = -rxt(k,177)*y(k,100)
         mat(k,760) = -rxt(k,184)*y(k,100)
         mat(k,875) = -rxt(k,195)*y(k,100)
         mat(k,618) = -(rxt(k,236) + rxt(k,237)) * y(k,100)
         mat(k,479) = -rxt(k,243)*y(k,100)
         mat(k,371) = -rxt(k,288)*y(k,100)
         mat(k,839) = -rxt(k,294)*y(k,100)
         mat(k,417) = -rxt(k,357)*y(k,100)
         mat(k,839) = mat(k,839) + rxt(k,300)*y(k,109)
         mat(k,587) = rxt(k,352)*y(k,61) + rxt(k,353)*y(k,66) + rxt(k,287)*y(k,97) &
                      + rxt(k,232)*y(k,98)
         mat(k,417) = mat(k,417) + rxt(k,239)*y(k,26) + rxt(k,358)*y(k,59)
         mat(k,479) = mat(k,479) + rxt(k,239)*y(k,20) + rxt(k,251)*y(k,109)
         mat(k,156) = rxt(k,362)*y(k,109)
         mat(k,110) = .500_r8*rxt(k,388)*y(k,109)
         mat(k,358) = mat(k,358) + rxt(k,146)*y(k,67)
         mat(k,181) = rxt(k,140)*y(k,66) + rxt(k,235)*y(k,98) + rxt(k,154)*y(k,109)
         mat(k,760) = mat(k,760) + rxt(k,358)*y(k,20)
         mat(k,445) = mat(k,445) + rxt(k,352)*y(k,16) + rxt(k,180)*y(k,109)
         mat(k,561) = mat(k,561) + rxt(k,353)*y(k,16) + rxt(k,140)*y(k,41)
         mat(k,654) = rxt(k,146)*y(k,38)
         mat(k,709) = mat(k,709) + rxt(k,157)*y(k,109)
         mat(k,307) = rxt(k,389)*y(k,109)
         mat(k,371) = mat(k,371) + rxt(k,287)*y(k,16)
         mat(k,618) = mat(k,618) + rxt(k,232)*y(k,16) + rxt(k,235)*y(k,41)
         mat(k,807) = mat(k,807) + rxt(k,300)*y(k,4) + rxt(k,251)*y(k,26) + rxt(k,362) &
                      *y(k,29) + .500_r8*rxt(k,388)*y(k,33) + rxt(k,154)*y(k,41) &
                      + rxt(k,180)*y(k,61) + rxt(k,157)*y(k,68) + rxt(k,389)*y(k,75)
         mat(k,140) = -(rxt(k,164)*y(k,66) + rxt(k,165)*y(k,67))
         mat(k,534) = -rxt(k,164)*y(k,101)
         mat(k,640) = -rxt(k,165)*y(k,101)
         mat(k,534) = mat(k,534) + rxt(k,434)*y(k,102)
         mat(k,254) = .900_r8*rxt(k,432)*y(k,102) + .800_r8*rxt(k,430)*y(k,103)
         mat(k,186) = rxt(k,434)*y(k,66) + .900_r8*rxt(k,432)*y(k,99)
         mat(k,246) = .800_r8*rxt(k,430)*y(k,99)
         mat(k,187) = -(rxt(k,432)*y(k,99) + rxt(k,433)*y(k,67) + (rxt(k,434) &
                      + rxt(k,435)) * y(k,66))
         mat(k,255) = -rxt(k,432)*y(k,102)
         mat(k,642) = -rxt(k,433)*y(k,102)
         mat(k,541) = -(rxt(k,434) + rxt(k,435)) * y(k,102)
         mat(k,247) = -(rxt(k,430)*y(k,99))
         mat(k,258) = -rxt(k,430)*y(k,103)
         mat(k,274) = rxt(k,439)*y(k,108)
         mat(k,752) = rxt(k,441)*y(k,108)
         mat(k,546) = rxt(k,434)*y(k,102)
         mat(k,646) = rxt(k,438)*y(k,104)
         mat(k,190) = rxt(k,434)*y(k,66)
         mat(k,136) = rxt(k,438)*y(k,67)
         mat(k,239) = rxt(k,439)*y(k,53) + rxt(k,441)*y(k,59)
         mat(k,133) = -(rxt(k,436)*y(k,66) + (rxt(k,437) + rxt(k,438)) * y(k,67))
         mat(k,533) = -rxt(k,436)*y(k,104)
         mat(k,639) = -(rxt(k,437) + rxt(k,438)) * y(k,104)
         mat(k,687) = -(rxt(k,111)*y(k,39) + rxt(k,112)*y(k,112) + (rxt(k,114) &
                      + rxt(k,115)) * y(k,67) + rxt(k,116)*y(k,68) + (rxt(k,271) &
                      + rxt(k,272)) * y(k,47) + (rxt(k,313) + rxt(k,314)) * y(k,43) &
                      + rxt(k,326)*y(k,31) + rxt(k,327)*y(k,32))
         mat(k,737) = -rxt(k,111)*y(k,105)
         mat(k,907) = -rxt(k,112)*y(k,105)
         mat(k,659) = -(rxt(k,114) + rxt(k,115)) * y(k,105)
         mat(k,714) = -rxt(k,116)*y(k,105)
         mat(k,404) = -(rxt(k,271) + rxt(k,272)) * y(k,105)
         mat(k,300) = -(rxt(k,313) + rxt(k,314)) * y(k,105)
         mat(k,43) = -rxt(k,326)*y(k,105)
         mat(k,60) = -rxt(k,327)*y(k,105)
         mat(k,659) = mat(k,659) + rxt(k,165)*y(k,101)
         mat(k,265) = .850_r8*rxt(k,431)*y(k,108)
         mat(k,143) = rxt(k,165)*y(k,67)
         mat(k,244) = .850_r8*rxt(k,431)*y(k,99)
         mat(k,33) = -(rxt(k,119)*y(k,66) + rxt(k,120)*y(k,67))
         mat(k,524) = -rxt(k,119)*y(k,106)
         mat(k,632) = -rxt(k,120)*y(k,106)
         mat(k,524) = mat(k,524) + rxt(k,123)*y(k,107)
         mat(k,632) = mat(k,632) + rxt(k,124)*y(k,107)
         mat(k,695) = rxt(k,125)*y(k,107)
         mat(k,35) = rxt(k,123)*y(k,66) + rxt(k,124)*y(k,67) + rxt(k,125)*y(k,68)
         mat(k,36) = -(rxt(k,123)*y(k,66) + rxt(k,124)*y(k,67) + rxt(k,125)*y(k,68))
         mat(k,525) = -rxt(k,123)*y(k,107)
         mat(k,633) = -rxt(k,124)*y(k,107)
         mat(k,696) = -rxt(k,125)*y(k,107)
         mat(k,633) = mat(k,633) + rxt(k,114)*y(k,105)
         mat(k,667) = rxt(k,114)*y(k,67)
         mat(k,238) = -(rxt(k,431)*y(k,99) + rxt(k,439)*y(k,53) + rxt(k,441)*y(k,59))
         mat(k,257) = -rxt(k,431)*y(k,108)
         mat(k,273) = -rxt(k,439)*y(k,108)
         mat(k,751) = -rxt(k,441)*y(k,108)
         mat(k,645) = rxt(k,433)*y(k,102) + rxt(k,437)*y(k,104) + rxt(k,444)*y(k,110)
         mat(k,189) = rxt(k,433)*y(k,67)
         mat(k,135) = rxt(k,437)*y(k,67)
         mat(k,208) = rxt(k,444)*y(k,67)
         mat(k,816) = -(rxt(k,153)*y(k,39) + rxt(k,154)*y(k,41) + rxt(k,155)*y(k,100) &
                      + rxt(k,156)*y(k,66) + rxt(k,157)*y(k,68) + (4._r8*rxt(k,158) &
                      + 4._r8*rxt(k,159)) * y(k,109) + rxt(k,163)*y(k,50) + rxt(k,180) &
                      *y(k,61) + rxt(k,183)*y(k,53) + rxt(k,197)*y(k,60) + rxt(k,203) &
                      *y(k,49) + rxt(k,248)*y(k,27) + (rxt(k,251) + rxt(k,252) &
                      ) * y(k,26) + rxt(k,258)*y(k,47) + rxt(k,262)*y(k,52) + rxt(k,300) &
                      *y(k,4) + rxt(k,306)*y(k,43) + rxt(k,354)*y(k,16) + rxt(k,360) &
                      *y(k,21) + (rxt(k,362) + rxt(k,366)) * y(k,29) + (rxt(k,371) &
                      + rxt(k,388)) * y(k,33) + rxt(k,374)*y(k,70) + rxt(k,379) &
                      *y(k,73) + rxt(k,384)*y(k,74) + rxt(k,389)*y(k,75))
         mat(k,741) = -rxt(k,153)*y(k,109)
         mat(k,184) = -rxt(k,154)*y(k,109)
         mat(k,519) = -rxt(k,155)*y(k,109)
         mat(k,570) = -rxt(k,156)*y(k,109)
         mat(k,718) = -rxt(k,157)*y(k,109)
         mat(k,130) = -rxt(k,163)*y(k,109)
         mat(k,454) = -rxt(k,180)*y(k,109)
         mat(k,283) = -rxt(k,183)*y(k,109)
         mat(k,884) = -rxt(k,197)*y(k,109)
         mat(k,331) = -rxt(k,203)*y(k,109)
         mat(k,321) = -rxt(k,248)*y(k,109)
         mat(k,488) = -(rxt(k,251) + rxt(k,252)) * y(k,109)
         mat(k,407) = -rxt(k,258)*y(k,109)
         mat(k,292) = -rxt(k,262)*y(k,109)
         mat(k,848) = -rxt(k,300)*y(k,109)
         mat(k,301) = -rxt(k,306)*y(k,109)
         mat(k,596) = -rxt(k,354)*y(k,109)
         mat(k,116) = -rxt(k,360)*y(k,109)
         mat(k,157) = -(rxt(k,362) + rxt(k,366)) * y(k,109)
         mat(k,111) = -(rxt(k,371) + rxt(k,388)) * y(k,109)
         mat(k,103) = -rxt(k,374)*y(k,109)
         mat(k,205) = -rxt(k,379)*y(k,109)
         mat(k,351) = -rxt(k,384)*y(k,109)
         mat(k,309) = -rxt(k,389)*y(k,109)
         mat(k,596) = mat(k,596) + rxt(k,353)*y(k,66)
         mat(k,116) = mat(k,116) + .300_r8*rxt(k,360)*y(k,109)
         mat(k,226) = rxt(k,323)*y(k,112)
         mat(k,363) = rxt(k,151)*y(k,68) + 2.000_r8*rxt(k,144)*y(k,100)
         mat(k,741) = mat(k,741) + rxt(k,139)*y(k,66) + rxt(k,111)*y(k,105)
         mat(k,184) = mat(k,184) + rxt(k,140)*y(k,66)
         mat(k,301) = mat(k,301) + rxt(k,305)*y(k,66) + rxt(k,313)*y(k,105)
         mat(k,407) = mat(k,407) + rxt(k,257)*y(k,66) + rxt(k,271)*y(k,105)
         mat(k,235) = rxt(k,308)*y(k,66)
         mat(k,292) = mat(k,292) + rxt(k,261)*y(k,66)
         mat(k,769) = rxt(k,184)*y(k,100)
         mat(k,454) = mat(k,454) + rxt(k,177)*y(k,100)
         mat(k,570) = mat(k,570) + rxt(k,353)*y(k,16) + rxt(k,139)*y(k,39) &
                      + rxt(k,140)*y(k,41) + rxt(k,305)*y(k,43) + rxt(k,257)*y(k,47) &
                      + rxt(k,308)*y(k,51) + rxt(k,261)*y(k,52) + rxt(k,147)*y(k,100)
         mat(k,718) = mat(k,718) + rxt(k,151)*y(k,38) + rxt(k,148)*y(k,100)
         mat(k,627) = rxt(k,237)*y(k,100)
         mat(k,519) = mat(k,519) + 2.000_r8*rxt(k,144)*y(k,38) + rxt(k,184)*y(k,59) &
                      + rxt(k,177)*y(k,61) + rxt(k,147)*y(k,66) + rxt(k,148)*y(k,68) &
                      + rxt(k,237)*y(k,98)
         mat(k,691) = rxt(k,111)*y(k,39) + rxt(k,313)*y(k,43) + rxt(k,271)*y(k,47) &
                      + 2.000_r8*rxt(k,112)*y(k,112)
         mat(k,816) = mat(k,816) + .300_r8*rxt(k,360)*y(k,21)
         mat(k,911) = rxt(k,323)*y(k,37) + 2.000_r8*rxt(k,112)*y(k,105)
         mat(k,207) = -(rxt(k,444)*y(k,67))
         mat(k,644) = -rxt(k,444)*y(k,110)
         mat(k,543) = rxt(k,435)*y(k,102) + rxt(k,436)*y(k,104)
         mat(k,188) = rxt(k,435)*y(k,66)
         mat(k,134) = rxt(k,436)*y(k,66)
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
         mat(k,78) = -(rxt(k,132)*y(k,39) + rxt(k,133)*y(k,112) + (rxt(k,135) &
                      + rxt(k,136)) * y(k,67) + rxt(k,137)*y(k,68) + (rxt(k,285) &
                      + rxt(k,286)) * y(k,47) + (rxt(k,319) + rxt(k,320)) * y(k,43) &
                      + rxt(k,328)*y(k,31) + rxt(k,329)*y(k,32))
         mat(k,723) = -rxt(k,132)*y(k,111)
         mat(k,891) = -rxt(k,133)*y(k,111)
         mat(k,635) = -(rxt(k,135) + rxt(k,136)) * y(k,111)
         mat(k,697) = -rxt(k,137)*y(k,111)
         mat(k,383) = -(rxt(k,285) + rxt(k,286)) * y(k,111)
         mat(k,294) = -(rxt(k,319) + rxt(k,320)) * y(k,111)
         mat(k,41) = -rxt(k,328)*y(k,111)
         mat(k,57) = -rxt(k,329)*y(k,111)
         mat(k,914) = -(rxt(k,112)*y(k,105) + rxt(k,323)*y(k,37) + rxt(k,390)*y(k,76))
         mat(k,694) = -rxt(k,112)*y(k,112)
         mat(k,227) = -rxt(k,323)*y(k,112)
         mat(k,73) = -rxt(k,390)*y(k,112)
         mat(k,599) = rxt(k,354)*y(k,109)
         mat(k,117) = rxt(k,360)*y(k,109)
         mat(k,364) = rxt(k,145)*y(k,100)
         mat(k,744) = rxt(k,153)*y(k,109)
         mat(k,185) = rxt(k,154)*y(k,109)
         mat(k,303) = rxt(k,306)*y(k,109)
         mat(k,410) = (rxt(k,403)+rxt(k,408))*y(k,51) + (rxt(k,396)+rxt(k,402) &
                       +rxt(k,407))*y(k,52) + rxt(k,258)*y(k,109)
         mat(k,333) = rxt(k,203)*y(k,109)
         mat(k,132) = rxt(k,163)*y(k,109)
         mat(k,237) = (rxt(k,403)+rxt(k,408))*y(k,47)
         mat(k,293) = (rxt(k,396)+rxt(k,402)+rxt(k,407))*y(k,47) + rxt(k,262)*y(k,109)
         mat(k,522) = rxt(k,145)*y(k,38) + rxt(k,155)*y(k,109)
         mat(k,819) = rxt(k,354)*y(k,16) + rxt(k,360)*y(k,21) + rxt(k,153)*y(k,39) &
                      + rxt(k,154)*y(k,41) + rxt(k,306)*y(k,43) + rxt(k,258)*y(k,47) &
                      + rxt(k,203)*y(k,49) + rxt(k,163)*y(k,50) + rxt(k,262)*y(k,52) &
                      + rxt(k,155)*y(k,100) + 2.000_r8*rxt(k,158)*y(k,109)
      end do
      end subroutine nlnmat05
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
         mat(k, 27) = mat(k, 27) + lmat(k, 27)
         mat(k, 30) = lmat(k, 30)
         mat(k, 31) = lmat(k, 31)
         mat(k, 32) = lmat(k, 32)
         mat(k, 33) = mat(k, 33) + lmat(k, 33)
         mat(k, 34) = mat(k, 34) + lmat(k, 34)
         mat(k, 35) = mat(k, 35) + lmat(k, 35)
         mat(k, 36) = mat(k, 36) + lmat(k, 36)
         mat(k, 37) = lmat(k, 37)
         mat(k, 38) = lmat(k, 38)
         mat(k, 39) = lmat(k, 39)
         mat(k, 40) = mat(k, 40) + lmat(k, 40)
         mat(k, 42) = mat(k, 42) + lmat(k, 42)
         mat(k, 44) = mat(k, 44) + lmat(k, 44)
         mat(k, 45) = lmat(k, 45)
         mat(k, 46) = mat(k, 46) + lmat(k, 46)
         mat(k, 47) = mat(k, 47) + lmat(k, 47)
         mat(k, 49) = lmat(k, 49)
         mat(k, 50) = lmat(k, 50)
         mat(k, 51) = lmat(k, 51)
         mat(k, 52) = lmat(k, 52)
         mat(k, 53) = lmat(k, 53)
         mat(k, 54) = lmat(k, 54)
         mat(k, 55) = lmat(k, 55)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 61) = lmat(k, 61)
         mat(k, 62) = lmat(k, 62)
         mat(k, 63) = lmat(k, 63)
         mat(k, 64) = lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = lmat(k, 68)
         mat(k, 70) = mat(k, 70) + lmat(k, 70)
         mat(k, 71) = lmat(k, 71)
         mat(k, 72) = lmat(k, 72)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 77) = lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 80) = lmat(k, 80)
         mat(k, 81) = lmat(k, 81)
         mat(k, 82) = lmat(k, 82)
         mat(k, 83) = lmat(k, 83)
         mat(k, 84) = lmat(k, 84)
         mat(k, 85) = lmat(k, 85)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 89) = mat(k, 89) + lmat(k, 89)
         mat(k, 90) = lmat(k, 90)
         mat(k, 91) = lmat(k, 91)
         mat(k, 92) = mat(k, 92) + lmat(k, 92)
         mat(k, 93) = lmat(k, 93)
         mat(k, 95) = mat(k, 95) + lmat(k, 95)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 98) = lmat(k, 98)
         mat(k, 104) = mat(k, 104) + lmat(k, 104)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 113) = lmat(k, 113)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 118) = lmat(k, 118)
         mat(k, 119) = lmat(k, 119)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 125) = lmat(k, 125)
         mat(k, 126) = mat(k, 126) + lmat(k, 126)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = lmat(k, 128)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 133) = mat(k, 133) + lmat(k, 133)
         mat(k, 140) = mat(k, 140) + lmat(k, 140)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 147) = lmat(k, 147)
         mat(k, 148) = lmat(k, 148)
         mat(k, 149) = lmat(k, 149)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 160) = lmat(k, 160)
         mat(k, 163) = lmat(k, 163)
         mat(k, 164) = lmat(k, 164)
         mat(k, 166) = mat(k, 166) + lmat(k, 166)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 172) = mat(k, 172) + lmat(k, 172)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 177) = mat(k, 177) + lmat(k, 177)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 187) = mat(k, 187) + lmat(k, 187)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 206) = lmat(k, 206)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 208) = mat(k, 208) + lmat(k, 208)
         mat(k, 209) = lmat(k, 209)
         mat(k, 210) = lmat(k, 210)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 220) = mat(k, 220) + lmat(k, 220)
         mat(k, 223) = lmat(k, 223)
         mat(k, 230) = mat(k, 230) + lmat(k, 230)
         mat(k, 231) = lmat(k, 231)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 245) = mat(k, 245) + lmat(k, 245)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 270) = lmat(k, 270)
         mat(k, 275) = lmat(k, 275)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 287) = mat(k, 287) + lmat(k, 287)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 306) = lmat(k, 306)
         mat(k, 308) = lmat(k, 308)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 314) = mat(k, 314) + lmat(k, 314)
         mat(k, 315) = mat(k, 315) + lmat(k, 315)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 318) = mat(k, 318) + lmat(k, 318)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 322) = lmat(k, 322)
         mat(k, 326) = mat(k, 326) + lmat(k, 326)
         mat(k, 331) = mat(k, 331) + lmat(k, 331)
         mat(k, 332) = lmat(k, 332)
         mat(k, 339) = lmat(k, 339)
         mat(k, 341) = mat(k, 341) + lmat(k, 341)
         mat(k, 346) = mat(k, 346) + lmat(k, 346)
         mat(k, 357) = mat(k, 357) + lmat(k, 357)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 402) = mat(k, 402) + lmat(k, 402)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 456) = mat(k, 456) + lmat(k, 456)
         mat(k, 478) = mat(k, 478) + lmat(k, 478)
         mat(k, 480) = mat(k, 480) + lmat(k, 480)
         mat(k, 482) = mat(k, 482) + lmat(k, 482)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 543) = mat(k, 543) + lmat(k, 543)
         mat(k, 547) = lmat(k, 547)
         mat(k, 562) = mat(k, 562) + lmat(k, 562)
         mat(k, 575) = mat(k, 575) + lmat(k, 575)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 594) = lmat(k, 594)
         mat(k, 604) = mat(k, 604) + lmat(k, 604)
         mat(k, 613) = lmat(k, 613)
         mat(k, 614) = mat(k, 614) + lmat(k, 614)
         mat(k, 615) = lmat(k, 615)
         mat(k, 618) = mat(k, 618) + lmat(k, 618)
         mat(k, 621) = mat(k, 621) + lmat(k, 621)
         mat(k, 644) = mat(k, 644) + lmat(k, 644)
         mat(k, 645) = mat(k, 645) + lmat(k, 645)
         mat(k, 647) = lmat(k, 647)
         mat(k, 655) = mat(k, 655) + lmat(k, 655)
         mat(k, 658) = mat(k, 658) + lmat(k, 658)
         mat(k, 659) = mat(k, 659) + lmat(k, 659)
         mat(k, 668) = mat(k, 668) + lmat(k, 668)
         mat(k, 669) = mat(k, 669) + lmat(k, 669)
         mat(k, 673) = mat(k, 673) + lmat(k, 673)
         mat(k, 676) = mat(k, 676) + lmat(k, 676)
         mat(k, 677) = mat(k, 677) + lmat(k, 677)
         mat(k, 679) = lmat(k, 679)
         mat(k, 682) = lmat(k, 682)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 684) = lmat(k, 684)
         mat(k, 685) = mat(k, 685) + lmat(k, 685)
         mat(k, 686) = mat(k, 686) + lmat(k, 686)
         mat(k, 687) = mat(k, 687) + lmat(k, 687)
         mat(k, 689) = mat(k, 689) + lmat(k, 689)
         mat(k, 690) = lmat(k, 690)
         mat(k, 691) = mat(k, 691) + lmat(k, 691)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 710) = mat(k, 710) + lmat(k, 710)
         mat(k, 713) = mat(k, 713) + lmat(k, 713)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 715) = mat(k, 715) + lmat(k, 715)
         mat(k, 739) = mat(k, 739) + lmat(k, 739)
         mat(k, 752) = mat(k, 752) + lmat(k, 752)
         mat(k, 753) = lmat(k, 753)
         mat(k, 754) = mat(k, 754) + lmat(k, 754)
         mat(k, 761) = mat(k, 761) + lmat(k, 761)
         mat(k, 768) = mat(k, 768) + lmat(k, 768)
         mat(k, 773) = lmat(k, 773)
         mat(k, 776) = lmat(k, 776)
         mat(k, 802) = mat(k, 802) + lmat(k, 802)
         mat(k, 804) = mat(k, 804) + lmat(k, 804)
         mat(k, 807) = mat(k, 807) + lmat(k, 807)
         mat(k, 810) = mat(k, 810) + lmat(k, 810)
         mat(k, 816) = mat(k, 816) + lmat(k, 816)
         mat(k, 819) = mat(k, 819) + lmat(k, 819)
         mat(k, 834) = mat(k, 834) + lmat(k, 834)
         mat(k, 840) = mat(k, 840) + lmat(k, 840)
         mat(k, 849) = mat(k, 849) + lmat(k, 849)
         mat(k, 867) = mat(k, 867) + lmat(k, 867)
         mat(k, 876) = mat(k, 876) + lmat(k, 876)
         mat(k, 883) = mat(k, 883) + lmat(k, 883)
         mat(k, 884) = mat(k, 884) + lmat(k, 884)
         mat(k, 886) = mat(k, 886) + lmat(k, 886)
         mat(k, 897) = lmat(k, 897)
         mat(k, 903) = lmat(k, 903)
         mat(k, 907) = mat(k, 907) + lmat(k, 907)
         mat(k, 909) = lmat(k, 909)
         mat(k, 911) = mat(k, 911) + lmat(k, 911)
         mat(k, 914) = mat(k, 914) + lmat(k, 914)
         mat(k, 195) = 0._r8
         mat(k, 196) = 0._r8
         mat(k, 198) = 0._r8
         mat(k, 211) = 0._r8
         mat(k, 212) = 0._r8
         mat(k, 215) = 0._r8
         mat(k, 219) = 0._r8
         mat(k, 234) = 0._r8
         mat(k, 251) = 0._r8
         mat(k, 252) = 0._r8
         mat(k, 253) = 0._r8
         mat(k, 256) = 0._r8
         mat(k, 261) = 0._r8
         mat(k, 262) = 0._r8
         mat(k, 264) = 0._r8
         mat(k, 266) = 0._r8
         mat(k, 267) = 0._r8
         mat(k, 272) = 0._r8
         mat(k, 278) = 0._r8
         mat(k, 281) = 0._r8
         mat(k, 310) = 0._r8
         mat(k, 323) = 0._r8
         mat(k, 327) = 0._r8
         mat(k, 328) = 0._r8
         mat(k, 330) = 0._r8
         mat(k, 336) = 0._r8
         mat(k, 345) = 0._r8
         mat(k, 354) = 0._r8
         mat(k, 356) = 0._r8
         mat(k, 366) = 0._r8
         mat(k, 369) = 0._r8
         mat(k, 372) = 0._r8
         mat(k, 375) = 0._r8
         mat(k, 377) = 0._r8
         mat(k, 378) = 0._r8
         mat(k, 380) = 0._r8
         mat(k, 387) = 0._r8
         mat(k, 388) = 0._r8
         mat(k, 394) = 0._r8
         mat(k, 396) = 0._r8
         mat(k, 397) = 0._r8
         mat(k, 399) = 0._r8
         mat(k, 401) = 0._r8
         mat(k, 403) = 0._r8
         mat(k, 405) = 0._r8
         mat(k, 406) = 0._r8
         mat(k, 408) = 0._r8
         mat(k, 409) = 0._r8
         mat(k, 414) = 0._r8
         mat(k, 418) = 0._r8
         mat(k, 422) = 0._r8
         mat(k, 423) = 0._r8
         mat(k, 425) = 0._r8
         mat(k, 427) = 0._r8
         mat(k, 434) = 0._r8
         mat(k, 435) = 0._r8
         mat(k, 439) = 0._r8
         mat(k, 440) = 0._r8
         mat(k, 441) = 0._r8
         mat(k, 442) = 0._r8
         mat(k, 444) = 0._r8
         mat(k, 448) = 0._r8
         mat(k, 450) = 0._r8
         mat(k, 451) = 0._r8
         mat(k, 452) = 0._r8
         mat(k, 455) = 0._r8
         mat(k, 457) = 0._r8
         mat(k, 465) = 0._r8
         mat(k, 466) = 0._r8
         mat(k, 471) = 0._r8
         mat(k, 473) = 0._r8
         mat(k, 477) = 0._r8
         mat(k, 484) = 0._r8
         mat(k, 485) = 0._r8
         mat(k, 486) = 0._r8
         mat(k, 491) = 0._r8
         mat(k, 512) = 0._r8
         mat(k, 515) = 0._r8
         mat(k, 537) = 0._r8
         mat(k, 542) = 0._r8
         mat(k, 545) = 0._r8
         mat(k, 551) = 0._r8
         mat(k, 553) = 0._r8
         mat(k, 558) = 0._r8
         mat(k, 566) = 0._r8
         mat(k, 573) = 0._r8
         mat(k, 578) = 0._r8
         mat(k, 584) = 0._r8
         mat(k, 586) = 0._r8
         mat(k, 591) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 593) = 0._r8
         mat(k, 595) = 0._r8
         mat(k, 597) = 0._r8
         mat(k, 598) = 0._r8
         mat(k, 602) = 0._r8
         mat(k, 605) = 0._r8
         mat(k, 607) = 0._r8
         mat(k, 611) = 0._r8
         mat(k, 619) = 0._r8
         mat(k, 623) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 628) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 630) = 0._r8
         mat(k, 641) = 0._r8
         mat(k, 652) = 0._r8
         mat(k, 653) = 0._r8
         mat(k, 656) = 0._r8
         mat(k, 657) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 663) = 0._r8
         mat(k, 664) = 0._r8
         mat(k, 665) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 670) = 0._r8
         mat(k, 672) = 0._r8
         mat(k, 675) = 0._r8
         mat(k, 680) = 0._r8
         mat(k, 693) = 0._r8
         mat(k, 711) = 0._r8
         mat(k, 716) = 0._r8
         mat(k, 721) = 0._r8
         mat(k, 726) = 0._r8
         mat(k, 729) = 0._r8
         mat(k, 730) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 732) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 738) = 0._r8
         mat(k, 740) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 743) = 0._r8
         mat(k, 748) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 765) = 0._r8
         mat(k, 767) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 778) = 0._r8
         mat(k, 793) = 0._r8
         mat(k, 812) = 0._r8
         mat(k, 826) = 0._r8
         mat(k, 827) = 0._r8
         mat(k, 831) = 0._r8
         mat(k, 833) = 0._r8
         mat(k, 835) = 0._r8
         mat(k, 836) = 0._r8
         mat(k, 837) = 0._r8
         mat(k, 841) = 0._r8
         mat(k, 844) = 0._r8
         mat(k, 845) = 0._r8
         mat(k, 846) = 0._r8
         mat(k, 851) = 0._r8
         mat(k, 860) = 0._r8
         mat(k, 862) = 0._r8
         mat(k, 863) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 870) = 0._r8
         mat(k, 871) = 0._r8
         mat(k, 872) = 0._r8
         mat(k, 877) = 0._r8
         mat(k, 878) = 0._r8
         mat(k, 880) = 0._r8
         mat(k, 882) = 0._r8
         mat(k, 887) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 894) = 0._r8
         mat(k, 895) = 0._r8
         mat(k, 896) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 900) = 0._r8
         mat(k, 901) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 904) = 0._r8
         mat(k, 905) = 0._r8
         mat(k, 906) = 0._r8
         mat(k, 908) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 912) = 0._r8
         mat(k, 913) = 0._r8
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
         mat(k, 24) = mat(k, 24) - dti(k)
         mat(k, 27) = mat(k, 27) - dti(k)
         mat(k, 30) = mat(k, 30) - dti(k)
         mat(k, 33) = mat(k, 33) - dti(k)
         mat(k, 36) = mat(k, 36) - dti(k)
         mat(k, 38) = mat(k, 38) - dti(k)
         mat(k, 40) = mat(k, 40) - dti(k)
         mat(k, 44) = mat(k, 44) - dti(k)
         mat(k, 47) = mat(k, 47) - dti(k)
         mat(k, 49) = mat(k, 49) - dti(k)
         mat(k, 53) = mat(k, 53) - dti(k)
         mat(k, 56) = mat(k, 56) - dti(k)
         mat(k, 61) = mat(k, 61) - dti(k)
         mat(k, 64) = mat(k, 64) - dti(k)
         mat(k, 70) = mat(k, 70) - dti(k)
         mat(k, 74) = mat(k, 74) - dti(k)
         mat(k, 78) = mat(k, 78) - dti(k)
         mat(k, 80) = mat(k, 80) - dti(k)
         mat(k, 86) = mat(k, 86) - dti(k)
         mat(k, 92) = mat(k, 92) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 104) = mat(k, 104) - dti(k)
         mat(k, 112) = mat(k, 112) - dti(k)
         mat(k, 120) = mat(k, 120) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 126) = mat(k, 126) - dti(k)
         mat(k, 133) = mat(k, 133) - dti(k)
         mat(k, 140) = mat(k, 140) - dti(k)
         mat(k, 146) = mat(k, 146) - dti(k)
         mat(k, 154) = mat(k, 154) - dti(k)
         mat(k, 158) = mat(k, 158) - dti(k)
         mat(k, 166) = mat(k, 166) - dti(k)
         mat(k, 177) = mat(k, 177) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 187) = mat(k, 187) - dti(k)
         mat(k, 199) = mat(k, 199) - dti(k)
         mat(k, 207) = mat(k, 207) - dti(k)
         mat(k, 220) = mat(k, 220) - dti(k)
         mat(k, 230) = mat(k, 230) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 247) = mat(k, 247) - dti(k)
         mat(k, 259) = mat(k, 259) - dti(k)
         mat(k, 276) = mat(k, 276) - dti(k)
         mat(k, 287) = mat(k, 287) - dti(k)
         mat(k, 296) = mat(k, 296) - dti(k)
         mat(k, 305) = mat(k, 305) - dti(k)
         mat(k, 314) = mat(k, 314) - dti(k)
         mat(k, 326) = mat(k, 326) - dti(k)
         mat(k, 341) = mat(k, 341) - dti(k)
         mat(k, 357) = mat(k, 357) - dti(k)
         mat(k, 370) = mat(k, 370) - dti(k)
         mat(k, 395) = mat(k, 395) - dti(k)
         mat(k, 415) = mat(k, 415) - dti(k)
         mat(k, 443) = mat(k, 443) - dti(k)
         mat(k, 478) = mat(k, 478) - dti(k)
         mat(k, 510) = mat(k, 510) - dti(k)
         mat(k, 562) = mat(k, 562) - dti(k)
         mat(k, 589) = mat(k, 589) - dti(k)
         mat(k, 621) = mat(k, 621) - dti(k)
         mat(k, 658) = mat(k, 658) - dti(k)
         mat(k, 687) = mat(k, 687) - dti(k)
         mat(k, 715) = mat(k, 715) - dti(k)
         mat(k, 739) = mat(k, 739) - dti(k)
         mat(k, 768) = mat(k, 768) - dti(k)
         mat(k, 816) = mat(k, 816) - dti(k)
         mat(k, 849) = mat(k, 849) - dti(k)
         mat(k, 886) = mat(k, 886) - dti(k)
         mat(k, 914) = mat(k, 914) - dti(k)
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
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
