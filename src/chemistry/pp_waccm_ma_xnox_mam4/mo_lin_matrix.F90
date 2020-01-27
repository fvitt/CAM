      module mo_lin_matrix
      use chem_mods, only: veclen
      private
      public :: linmat
      contains
      subroutine linmat01( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
      do k = 1,avec_len
         mat(k,1) = -( het_rates(k,1) )
         mat(k,2) = -( het_rates(k,2) )
         mat(k,61) = -( rxt(k,40) + het_rates(k,3) )
         mat(k,849) = -( rxt(k,41) + het_rates(k,4) )
         mat(k,152) = rxt(k,42)
         mat(k,146) = -( rxt(k,42) + rxt(k,43) + rxt(k,397) + rxt(k,400) + rxt(k,405) &
                 + het_rates(k,5) )
         mat(k,589) = -( rxt(k,34) + rxt(k,35) + het_rates(k,16) )
         mat(k,115) = rxt(k,36)
         mat(k,684) = rxt(k,364)*y(k,22) + rxt(k,365)*y(k,22)
         mat(k,415) = -( het_rates(k,20) )
         mat(k,615) = rxt(k,233)*y(k,22)
         mat(k,223) = rxt(k,321)*y(k,22)
         mat(k,804) = rxt(k,361)*y(k,22)
         mat(k,679) = rxt(k,363)*y(k,22)
         mat(k,112) = -( rxt(k,36) + het_rates(k,21) )
         mat(k,38) = -( rxt(k,57) + het_rates(k,24) )
         mat(k,21) = -( rxt(k,58) + rxt(k,274) + het_rates(k,25) )
         mat(k,478) = -( rxt(k,59) + het_rates(k,26) )
         mat(k,318) = rxt(k,61)
         mat(k,89) = rxt(k,73)
         mat(k,22) = 2.000_r8*rxt(k,274)
         mat(k,314) = -( rxt(k,60) + rxt(k,61) + rxt(k,399) + rxt(k,404) + rxt(k,410) &
                 + het_rates(k,27) )
         mat(k,154) = -( het_rates(k,29) )
         mat(k,575) = rxt(k,34) + rxt(k,35)
         mat(k,97) = rxt(k,105)
         mat(k,604) = rxt(k,335)*y(k,19)
         mat(k,206) = rxt(k,442)*y(k,30)
         mat(k,40) = -( rxt(k,62) + het_rates(k,31) )
         mat(k,668) = rxt(k,265)*y(k,8) + rxt(k,267)*y(k,11) + 2.000_r8*rxt(k,268)*y(k,12) &
                      + 2.000_r8*rxt(k,269)*y(k,13) + rxt(k,270)*y(k,14) &
                      + rxt(k,310)*y(k,9) + 2.000_r8*rxt(k,312)*y(k,40) &
                      + rxt(k,345)*y(k,45) + rxt(k,346)*y(k,46)
         mat(k,773) = rxt(k,340)*y(k,45) + rxt(k,341)*y(k,46)
         mat(k,56) = -( rxt(k,63) + het_rates(k,32) )
         mat(k,669) = rxt(k,266)*y(k,10) + rxt(k,267)*y(k,11) + rxt(k,344)*y(k,44)
         mat(k,776) = rxt(k,339)*y(k,44)
         mat(k,104) = -( het_rates(k,33) )
         mat(k,3) = -( het_rates(k,34) )
         mat(k,4) = -( het_rates(k,35) )
         mat(k,5) = -( het_rates(k,36) )
         mat(k,220) = -( rxt(k,321)*y(k,22) + het_rates(k,37) )
         mat(k,42) = 2.000_r8*rxt(k,62)
         mat(k,58) = rxt(k,63)
         mat(k,54) = rxt(k,70)
         mat(k,673) = rxt(k,269)*y(k,13) + rxt(k,310)*y(k,9)
         mat(k,357) = -( het_rates(k,38) )
         mat(k,897) = rxt(k,1) + 2.000_r8*rxt(k,3)
         mat(k,581) = 2.000_r8*rxt(k,34)
         mat(k,113) = rxt(k,36)
         mat(k,297) = rxt(k,65)
         mat(k,393) = rxt(k,69)
         mat(k,55) = rxt(k,70)
         mat(k,676) = rxt(k,364)*y(k,22)
         mat(k,739) = -( het_rates(k,39) )
         mat(k,909) = rxt(k,2)
         mat(k,594) = rxt(k,35)
         mat(k,689) = rxt(k,365)*y(k,22)
         mat(k,179) = -( rxt(k,4) + het_rates(k,41) )
         mat(k,500) = .500_r8*rxt(k,391)
         mat(k,24) = -( rxt(k,104) + het_rates(k,42) )
         mat(k,296) = -( rxt(k,65) + het_rates(k,43) )
         mat(k,395) = -( rxt(k,69) + het_rates(k,47) )
         mat(k,614) = rxt(k,233)*y(k,22) + rxt(k,330)*y(k,15) + rxt(k,332)*y(k,17) &
                      + 2.000_r8*rxt(k,335)*y(k,19) + rxt(k,337)*y(k,23)
         mat(k,53) = -( rxt(k,70) + het_rates(k,48) )
         mat(k,216) = rxt(k,321)*y(k,22)
         mat(k,326) = -( rxt(k,11) + het_rates(k,49) )
         mat(k,81) = 2.000_r8*rxt(k,392) + 2.000_r8*rxt(k,395) + 2.000_r8*rxt(k,398) &
                      + 2.000_r8*rxt(k,409)
         mat(k,867) = .500_r8*rxt(k,393)
         mat(k,438) = rxt(k,394)
         mat(k,148) = rxt(k,397) + rxt(k,400) + rxt(k,405)
         mat(k,315) = rxt(k,399) + rxt(k,404) + rxt(k,410)
         mat(k,126) = -( rxt(k,12) + rxt(k,13) + rxt(k,204) + het_rates(k,50) )
         mat(k,230) = -( rxt(k,71) + het_rates(k,51) )
         mat(k,147) = rxt(k,397) + rxt(k,400) + rxt(k,405)
         mat(k,287) = -( rxt(k,72) + het_rates(k,52) )
         mat(k,313) = rxt(k,399) + rxt(k,404) + rxt(k,410)
         mat(k,276) = -( rxt(k,79) + het_rates(k,53) )
         mat(k,754) = rxt(k,17)
         mat(k,210) = rxt(k,443)
         mat(k,80) = -( rxt(k,15) + rxt(k,16) + rxt(k,205) + rxt(k,392) + rxt(k,395) &
                      + rxt(k,398) + rxt(k,409) + het_rates(k,55) )
         mat(k,6) = -( het_rates(k,56) )
         mat(k,7) = -( het_rates(k,57) )
         mat(k,8) = -( het_rates(k,58) )
         mat(k,768) = -( rxt(k,17) + rxt(k,18) + het_rates(k,59) )
         mat(k,84) = rxt(k,16)
         mat(k,883) = rxt(k,19) + .500_r8*rxt(k,393)
         mat(k,453) = rxt(k,21)
         mat(k,245) = rxt(k,440)
         mat(k,690) = 2.000_r8*rxt(k,193)*y(k,54)
         mat(k,886) = -( rxt(k,19) + rxt(k,393) + het_rates(k,60) )
         mat(k,332) = rxt(k,11)
         mat(k,131) = rxt(k,13) + rxt(k,204)
         mat(k,85) = rxt(k,15) + rxt(k,205)
         mat(k,456) = rxt(k,20)
         mat(k,153) = rxt(k,42)
         mat(k,322) = rxt(k,61)
         mat(k,443) = -( rxt(k,20) + rxt(k,21) + rxt(k,394) + het_rates(k,61) )
         mat(k,127) = rxt(k,12)
         mat(k,82) = rxt(k,15) + rxt(k,16) + rxt(k,205)
         mat(k,150) = rxt(k,43)
         mat(k,317) = rxt(k,60)
         mat(k,9) = -( het_rates(k,62) )
         mat(k,10) = -( het_rates(k,63) )
         mat(k,11) = -( het_rates(k,64) )
         mat(k,12) = -( het_rates(k,65) )
         mat(k,562) = -( rxt(k,88) + rxt(k,89) + rxt(k,90) + rxt(k,91) + rxt(k,92) &
                      + rxt(k,93) + het_rates(k,66) )
         mat(k,903) = rxt(k,3)
         mat(k,655) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,94) + rxt(k,96) + rxt(k,98) &
                      + 2.000_r8*rxt(k,99) + 2.000_r8*rxt(k,100) + rxt(k,101) + rxt(k,102) &
                      + rxt(k,103)
         mat(k,710) = rxt(k,7)
         mat(k,83) = rxt(k,16)
         mat(k,761) = rxt(k,17)
         mat(k,876) = rxt(k,19)
         mat(k,446) = rxt(k,20)
         mat(k,840) = rxt(k,41)
         mat(k,480) = rxt(k,59)
         mat(k,90) = rxt(k,73)
         mat(k,346) = rxt(k,106)
         mat(k,308) = rxt(k,107)
         mat(k,72) = rxt(k,108)
         mat(k,683) = rxt(k,113)
         mat(k,658) = -( rxt(k,5) + rxt(k,6) + rxt(k,94) + rxt(k,95) + rxt(k,96) &
                      + rxt(k,97) + rxt(k,98) + rxt(k,99) + rxt(k,100) + rxt(k,101) &
                      + rxt(k,102) + rxt(k,103) + het_rates(k,67) )
         mat(k,713) = rxt(k,7)
         mat(k,449) = rxt(k,21)
         mat(k,34) = rxt(k,109) + rxt(k,118)
         mat(k,37) = rxt(k,110)
         mat(k,686) = rxt(k,194)*y(k,54)
         mat(k,715) = -( rxt(k,7) + rxt(k,8) + het_rates(k,68) )
         mat(k,86) = -( rxt(k,73) + het_rates(k,69) )
         mat(k,96) = -( rxt(k,105) + het_rates(k,70) )
         mat(k,13) = -( het_rates(k,71) )
         mat(k,14) = -( het_rates(k,72) )
         mat(k,199) = -( het_rates(k,73) )
         mat(k,98) = rxt(k,105)
         mat(k,339) = rxt(k,106)
         mat(k,341) = -( rxt(k,106) + het_rates(k,74) )
         mat(k,306) = rxt(k,107)
         mat(k,305) = -( rxt(k,107) + het_rates(k,75) )
         mat(k,71) = rxt(k,108)
         mat(k,70) = -( rxt(k,108) + het_rates(k,76) )
         mat(k,25) = rxt(k,104)
         mat(k,15) = -( het_rates(k,77) )
         mat(k,16) = -( het_rates(k,78) )
         mat(k,17) = -( het_rates(k,79) )
         mat(k,18) = -( het_rates(k,80) )
         mat(k,19) = -( het_rates(k,81) )
         mat(k,20) = -( het_rates(k,82) )
         mat(k,74) = -( rxt(k,74) + rxt(k,75) + rxt(k,412) + rxt(k,417) + rxt(k,423) &
                 + het_rates(k,83) )
         mat(k,122) = -( rxt(k,76) + rxt(k,77) + rxt(k,416) + rxt(k,422) + rxt(k,428) &
                 + het_rates(k,84) )
         mat(k,158) = -( rxt(k,22) + het_rates(k,85) )
         mat(k,75) = rxt(k,412) + rxt(k,417) + rxt(k,423)
         mat(k,67) = rxt(k,413) + rxt(k,419) + rxt(k,425)
         mat(k,31) = rxt(k,414) + rxt(k,420) + rxt(k,426)
         mat(k,50) = 2.000_r8*rxt(k,415) + 2.000_r8*rxt(k,421) + 2.000_r8*rxt(k,427)
         mat(k,123) = rxt(k,416) + rxt(k,422) + rxt(k,428)
         mat(k,44) = -( rxt(k,23) + rxt(k,24) + rxt(k,207) + het_rates(k,86) )
         mat(k,47) = -( het_rates(k,87) )
         mat(k,118) = rxt(k,25)
         mat(k,27) = -( het_rates(k,88) )
         mat(k,120) = -( rxt(k,25) + het_rates(k,89) )
         mat(k,173) = rxt(k,26)
         mat(k,66) = rxt(k,27)
         mat(k,164) = rxt(k,30)
         mat(k,177) = -( rxt(k,26) + het_rates(k,90) )
         mat(k,160) = rxt(k,22)
         mat(k,46) = rxt(k,24) + rxt(k,207)
         mat(k,68) = rxt(k,28) + rxt(k,208)
         mat(k,52) = rxt(k,29) + rxt(k,209)
         mat(k,167) = rxt(k,31)
         mat(k,77) = rxt(k,74)
         mat(k,125) = rxt(k,76)
         mat(k,64) = -( rxt(k,27) + rxt(k,28) + rxt(k,208) + rxt(k,413) + rxt(k,419) &
                      + rxt(k,425) + het_rates(k,91) )
         mat(k,49) = -( rxt(k,29) + rxt(k,209) + rxt(k,415) + rxt(k,421) + rxt(k,427) &
                 + het_rates(k,92) )
         mat(k,166) = -( rxt(k,30) + rxt(k,31) + het_rates(k,93) )
         mat(k,45) = rxt(k,23)
         mat(k,51) = rxt(k,29) + rxt(k,209)
         mat(k,32) = rxt(k,32) + rxt(k,33) + rxt(k,210)
         mat(k,76) = rxt(k,75)
         mat(k,124) = rxt(k,77)
         mat(k,30) = -( rxt(k,32) + rxt(k,33) + rxt(k,210) + rxt(k,414) + rxt(k,420) &
                      + rxt(k,426) + het_rates(k,94) )
         mat(k,95) = -( het_rates(k,95) )
         mat(k,93) = rxt(k,9)
         mat(k,119) = rxt(k,25)
         mat(k,172) = rxt(k,26)
         mat(k,65) = rxt(k,27)
         mat(k,163) = rxt(k,31)
         mat(k,79) = rxt(k,134)
      end do
      end subroutine linmat01
      subroutine linmat02( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
      do k = 1,avec_len
         mat(k,92) = -( rxt(k,9) + rxt(k,10) + het_rates(k,96) )
         mat(k,370) = -( het_rates(k,97) )
         mat(k,62) = rxt(k,40)
         mat(k,834) = rxt(k,41)
         mat(k,149) = rxt(k,43)
         mat(k,298) = rxt(k,65)
         mat(k,231) = rxt(k,71)
         mat(k,677) = rxt(k,265)*y(k,8) + rxt(k,310)*y(k,9) + 3.000_r8*rxt(k,311)*y(k,23) &
                      + 2.000_r8*rxt(k,312)*y(k,40) + 2.000_r8*rxt(k,342)*y(k,15) &
                      + rxt(k,343)*y(k,17)
         mat(k,613) = 2.000_r8*rxt(k,330)*y(k,15) + rxt(k,332)*y(k,17) &
                      + 3.000_r8*rxt(k,337)*y(k,23)
         mat(k,802) = 2.000_r8*rxt(k,331)*y(k,15) + rxt(k,333)*y(k,17) &
                      + 3.000_r8*rxt(k,338)*y(k,23)
         mat(k,621) = -( rxt(k,233)*y(k,22) + rxt(k,330)*y(k,15) + rxt(k,332)*y(k,17) &
                      + rxt(k,335)*y(k,19) + rxt(k,337)*y(k,23) + het_rates(k,98) )
         mat(k,63) = rxt(k,40)
         mat(k,39) = 2.000_r8*rxt(k,57)
         mat(k,23) = 2.000_r8*rxt(k,58)
         mat(k,482) = rxt(k,59)
         mat(k,320) = rxt(k,60)
         mat(k,59) = rxt(k,63)
         mat(k,402) = rxt(k,69)
         mat(k,291) = rxt(k,72)
         mat(k,685) = 4.000_r8*rxt(k,264)*y(k,7) + rxt(k,265)*y(k,8) &
                      + 2.000_r8*rxt(k,266)*y(k,10) + 2.000_r8*rxt(k,267)*y(k,11) &
                      + 2.000_r8*rxt(k,268)*y(k,12) + rxt(k,269)*y(k,13) &
                      + 2.000_r8*rxt(k,270)*y(k,14) + rxt(k,344)*y(k,44) &
                      + rxt(k,345)*y(k,45) + rxt(k,346)*y(k,46)
         mat(k,810) = 3.000_r8*rxt(k,334)*y(k,18) + rxt(k,336)*y(k,19) &
                      + rxt(k,339)*y(k,44) + rxt(k,340)*y(k,45) + rxt(k,341)*y(k,46)
         mat(k,259) = -( het_rates(k,99) )
         mat(k,753) = rxt(k,18)
         mat(k,275) = rxt(k,79)
         mat(k,547) = rxt(k,88) + rxt(k,89) + rxt(k,90) + rxt(k,91) + rxt(k,92) &
                      + rxt(k,93)
         mat(k,647) = rxt(k,94) + rxt(k,95) + rxt(k,96) + rxt(k,97) + rxt(k,98) &
                      + rxt(k,101) + rxt(k,102) + rxt(k,103)
         mat(k,510) = -( rxt(k,391) + het_rates(k,100) )
         mat(k,128) = rxt(k,13) + rxt(k,204)
         mat(k,618) = rxt(k,332)*y(k,17) + rxt(k,335)*y(k,19)
         mat(k,807) = rxt(k,333)*y(k,17) + rxt(k,336)*y(k,19)
         mat(k,682) = rxt(k,364)*y(k,22)
         mat(k,140) = -( het_rates(k,101) )
         mat(k,187) = -( het_rates(k,102) )
         mat(k,247) = -( het_rates(k,103) )
         mat(k,752) = rxt(k,18)
         mat(k,239) = rxt(k,440)
         mat(k,209) = rxt(k,443)
         mat(k,133) = -( het_rates(k,104) )
         mat(k,270) = rxt(k,79)
         mat(k,687) = -( rxt(k,113) + rxt(k,193)*y(k,54) + rxt(k,194)*y(k,54) &
                      + rxt(k,264)*y(k,7) + rxt(k,265)*y(k,8) + rxt(k,266)*y(k,10) &
                      + rxt(k,267)*y(k,11) + rxt(k,268)*y(k,12) + rxt(k,269)*y(k,13) &
                      + rxt(k,270)*y(k,14) + rxt(k,310)*y(k,9) + rxt(k,311)*y(k,23) &
                      + rxt(k,312)*y(k,40) + rxt(k,342)*y(k,15) + rxt(k,343)*y(k,17) &
                      + rxt(k,344)*y(k,44) + rxt(k,345)*y(k,45) + rxt(k,346)*y(k,46) &
                      + rxt(k,363)*y(k,22) + rxt(k,364)*y(k,22) + rxt(k,365)*y(k,22) &
                 + het_rates(k,105) )
         mat(k,907) = rxt(k,2)
         mat(k,659) = rxt(k,6)
         mat(k,714) = rxt(k,8)
         mat(k,33) = -( rxt(k,109) + rxt(k,118) + het_rates(k,106) )
         mat(k,695) = rxt(k,8)
         mat(k,35) = rxt(k,122) + rxt(k,121)*y(k,30)
         mat(k,36) = -( rxt(k,110) + rxt(k,122) + rxt(k,121)*y(k,30) + het_rates(k,107) )
         mat(k,238) = -( rxt(k,440) + het_rates(k,108) )
         mat(k,645) = rxt(k,95) + rxt(k,97)
         mat(k,208) = rxt(k,442)*y(k,30)
         mat(k,816) = -( rxt(k,331)*y(k,15) + rxt(k,333)*y(k,17) + rxt(k,334)*y(k,18) &
                      + rxt(k,336)*y(k,19) + rxt(k,338)*y(k,23) + rxt(k,339)*y(k,44) &
                      + rxt(k,340)*y(k,45) + rxt(k,341)*y(k,46) + rxt(k,361)*y(k,22) &
                 + het_rates(k,109) )
         mat(k,911) = rxt(k,1)
         mat(k,184) = 2.000_r8*rxt(k,4)
         mat(k,331) = rxt(k,11)
         mat(k,130) = rxt(k,12)
         mat(k,116) = rxt(k,36)
         mat(k,235) = rxt(k,71)
         mat(k,292) = rxt(k,72)
         mat(k,884) = .500_r8*rxt(k,393)
         mat(k,691) = rxt(k,363)*y(k,22)
         mat(k,207) = -( rxt(k,443) + rxt(k,442)*y(k,30) + het_rates(k,110) )
         mat(k,543) = rxt(k,88) + rxt(k,89) + rxt(k,90) + rxt(k,91) + rxt(k,92) &
                      + rxt(k,93)
         mat(k,644) = rxt(k,94) + rxt(k,96) + rxt(k,98) + rxt(k,101) + rxt(k,102) &
                      + rxt(k,103)
         mat(k,78) = -( rxt(k,134) + rxt(k,230)*y(k,54) + rxt(k,231)*y(k,54) &
                      + rxt(k,278)*y(k,7) + rxt(k,279)*y(k,8) + rxt(k,280)*y(k,10) &
                      + rxt(k,281)*y(k,11) + rxt(k,282)*y(k,12) + rxt(k,283)*y(k,13) &
                      + rxt(k,284)*y(k,14) + rxt(k,316)*y(k,9) + rxt(k,317)*y(k,23) &
                      + rxt(k,318)*y(k,40) + rxt(k,347)*y(k,15) + rxt(k,348)*y(k,17) &
                      + rxt(k,349)*y(k,44) + rxt(k,350)*y(k,45) + rxt(k,351)*y(k,46) &
                      + rxt(k,367)*y(k,22) + rxt(k,368)*y(k,22) + rxt(k,369)*y(k,22) &
                 + het_rates(k,111) )
         mat(k,91) = rxt(k,10)
         mat(k,914) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,112) )
         mat(k,26) = rxt(k,104)
         mat(k,819) = rxt(k,331)*y(k,15) + rxt(k,333)*y(k,17) + rxt(k,334)*y(k,18) &
                      + rxt(k,336)*y(k,19) + rxt(k,341)*y(k,46) + rxt(k,361)*y(k,22)
      end do
      end subroutine linmat02
      subroutine linmat( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
      call linmat01( avec_len, mat, y, rxt, het_rates )
      call linmat02( avec_len, mat, y, rxt, het_rates )
      end subroutine linmat
      end module mo_lin_matrix
