      module mo_indprd
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: indprd
      contains
      subroutine indprd( class, prod, nprod, y, extfrc, rxt, chnkpnts )
      use chem_mods, only : gas_pcnst, extcnt, rxntot
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: chnkpnts
      integer, intent(in) :: nprod
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(in) :: extfrc(chnkpnts,extcnt)
      real(r8), intent(inout) :: prod(chnkpnts,nprod)
!--------------------------------------------------------------------
! ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,1) = + extfrc(:,17)
         prod(:,2) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,9) = 0._r8
         prod(:,10) = 0._r8
         prod(:,11) = 0._r8
         prod(:,12) = 0._r8
         prod(:,13) = 0._r8
         prod(:,14) = 0._r8
         prod(:,15) =.100_r8*rxt(:,308)*y(:,149)*y(:,29)
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = (rxt(:,265)*y(:,68) +rxt(:,267)*y(:,96) +rxt(:,275)*y(:,68) + &
                 rxt(:,295)*y(:,53) +.500_r8*rxt(:,296)*y(:,54) + &
                 .800_r8*rxt(:,301)*y(:,83) +rxt(:,302)*y(:,84) + &
                 .500_r8*rxt(:,351)*y(:,118) +1.800_r8*rxt(:,461)*y(:,193))*y(:,227) &
                  + (2.000_r8*rxt(:,291)*y(:,210) +.900_r8*rxt(:,292)*y(:,211) + &
                 rxt(:,294)*y(:,136) +2.000_r8*rxt(:,341)*y(:,222) + &
                 rxt(:,365)*y(:,218) +rxt(:,390)*y(:,234))*y(:,210) &
                  + (.200_r8*rxt(:,308)*y(:,29) +.100_r8*rxt(:,352)*y(:,120) + &
                 .270_r8*rxt(:,440)*y(:,6) +.270_r8*rxt(:,443)*y(:,119))*y(:,149) &
                  + (rxt(:,342)*y(:,211) +.450_r8*rxt(:,343)*y(:,216) + &
                 2.000_r8*rxt(:,344)*y(:,222))*y(:,222) &
                  + (.500_r8*rxt(:,450)*y(:,211) +.900_r8*rxt(:,452)*y(:,136)) &
                 *y(:,231) +rxt(:,37)*y(:,54) +.400_r8*rxt(:,60)*y(:,154) +rxt(:,65) &
                 *y(:,189) +.800_r8*rxt(:,69)*y(:,193)
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) = 0._r8
         prod(:,23) = 0._r8
         prod(:,24) =rxt(:,151)*y(:,137)*y(:,121)
         prod(:,25) = 0._r8
         prod(:,26) = 0._r8
         prod(:,27) = 0._r8
         prod(:,28) = 0._r8
         prod(:,29) =rxt(:,479)*y(:,227)*y(:,129) +rxt(:,488)*y(:,130)
         prod(:,30) = (rxt(:,412)*y(:,212) +rxt(:,415)*y(:,221) +rxt(:,418)*y(:,223) + &
                 rxt(:,422)*y(:,156))*y(:,137) +.500_r8*rxt(:,351)*y(:,227)*y(:,118) &
                  +.200_r8*rxt(:,447)*y(:,225)*y(:,136) +.500_r8*rxt(:,459)*y(:,192) &
                 *y(:,138)
         prod(:,31) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,134) = 0._r8
         prod(:,135) = 0._r8
         prod(:,1) = + extfrc(:,10)
         prod(:,2) = + extfrc(:,11)
         prod(:,160) = 0._r8
         prod(:,61) = 0._r8
         prod(:,106) = 0._r8
         prod(:,62) = 0._r8
         prod(:,103) = 0._r8
         prod(:,111) = 0._r8
         prod(:,86) = 0._r8
         prod(:,128) = 0._r8
         prod(:,94) = 0._r8
         prod(:,77) = 0._r8
         prod(:,95) = 0._r8
         prod(:,189) =rxt(:,79)*y(:,37) +rxt(:,80)*y(:,38) +2.000_r8*rxt(:,86)*y(:,44) &
                  +rxt(:,87)*y(:,46) +3.000_r8*rxt(:,90)*y(:,58) +2.000_r8*rxt(:,98) &
                 *y(:,87)
         prod(:,74) = 0._r8
         prod(:,205) = 0._r8
         prod(:,126) = 0._r8
         prod(:,75) = 0._r8
         prod(:,90) = 0._r8
         prod(:,84) = 0._r8
         prod(:,125) = 0._r8
         prod(:,78) = 0._r8
         prod(:,91) = 0._r8
         prod(:,85) = 0._r8
         prod(:,164) = 0._r8
         prod(:,100) = 0._r8
         prod(:,54) = 0._r8
         prod(:,79) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,204) =.180_r8*rxt(:,40)*y(:,57)
         prod(:,176) = 0._r8
         prod(:,52) = 0._r8
         prod(:,162) = 0._r8
         prod(:,181) = 0._r8
         prod(:,124) = 0._r8
         prod(:,119) = 0._r8
         prod(:,148) = 0._r8
         prod(:,104) = 0._r8
         prod(:,194) =4.000_r8*rxt(:,78)*y(:,36) +rxt(:,79)*y(:,37) &
                  +2.000_r8*rxt(:,81)*y(:,39) +2.000_r8*rxt(:,82)*y(:,40) &
                  +2.000_r8*rxt(:,83)*y(:,41) +rxt(:,84)*y(:,42) +2.000_r8*rxt(:,85) &
                 *y(:,43) +3.000_r8*rxt(:,88)*y(:,47) +rxt(:,89)*y(:,49) +rxt(:,100) &
                 *y(:,91) +rxt(:,101)*y(:,92) +rxt(:,102)*y(:,93)
         prod(:,60) = 0._r8
         prod(:,53) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,201) = 0._r8
         prod(:,163) = 0._r8
         prod(:,170) =.380_r8*rxt(:,40)*y(:,57) +rxt(:,41)*y(:,69) + extfrc(:,9)
         prod(:,9) = 0._r8
         prod(:,10) = 0._r8
         prod(:,11) = 0._r8
         prod(:,55) =rxt(:,79)*y(:,37) +rxt(:,80)*y(:,38) +rxt(:,82)*y(:,40) &
                  +2.000_r8*rxt(:,83)*y(:,41) +2.000_r8*rxt(:,84)*y(:,42) +rxt(:,85) &
                 *y(:,43) +2.000_r8*rxt(:,98)*y(:,87) +rxt(:,101)*y(:,92) +rxt(:,102) &
                 *y(:,93)
         prod(:,65) =rxt(:,81)*y(:,39) +rxt(:,82)*y(:,40) +rxt(:,100)*y(:,91)
         prod(:,67) = 0._r8
         prod(:,82) = 0._r8
         prod(:,12) = 0._r8
         prod(:,13) = 0._r8
         prod(:,14) = 0._r8
         prod(:,56) = 0._r8
         prod(:,147) =rxt(:,80)*y(:,38) +rxt(:,84)*y(:,42)
         prod(:,166) = 0._r8
         prod(:,157) = 0._r8
         prod(:,191) = (rxt(:,39) +.330_r8*rxt(:,40))*y(:,57)
         prod(:,177) =1.440_r8*rxt(:,40)*y(:,57)
         prod(:,130) = 0._r8
         prod(:,57) = 0._r8
         prod(:,151) = 0._r8
         prod(:,203) = 0._r8
         prod(:,64) = 0._r8
         prod(:,149) = 0._r8
         prod(:,72) = 0._r8
         prod(:,190) = 0._r8
         prod(:,101) = 0._r8
         prod(:,146) = 0._r8
         prod(:,152) = 0._r8
         prod(:,169) = 0._r8
         prod(:,73) = 0._r8
         prod(:,171) = 0._r8
         prod(:,89) = 0._r8
         prod(:,58) = 0._r8
         prod(:,155) = 0._r8
         prod(:,129) = 0._r8
         prod(:,121) = 0._r8
         prod(:,179) = 0._r8
         prod(:,105) = 0._r8
         prod(:,138) = 0._r8
         prod(:,49) = 0._r8
         prod(:,180) = 0._r8
         prod(:,87) = 0._r8
         prod(:,117) = 0._r8
         prod(:,88) = 0._r8
         prod(:,123) = 0._r8
         prod(:,159) = 0._r8
         prod(:,184) = 0._r8
         prod(:,97) = + extfrc(:,16)
         prod(:,83) = 0._r8
         prod(:,98) = 0._r8
         prod(:,167) = 0._r8
         prod(:,15) = 0._r8
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,51) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,197) = + extfrc(:,2)
         prod(:,202) = + extfrc(:,3)
         prod(:,193) = 0._r8
         prod(:,22) = 0._r8
         prod(:,23) = 0._r8
         prod(:,24) = 0._r8
         prod(:,156) = 0._r8
         prod(:,99) = 0._r8
         prod(:,25) = + extfrc(:,12)
         prod(:,26) = + extfrc(:,13)
         prod(:,27) = 0._r8
         prod(:,28) = + extfrc(:,14)
         prod(:,200) =.180_r8*rxt(:,40)*y(:,57) +rxt(:,41)*y(:,69) + (rxt(:,5) + &
                 2.000_r8*rxt(:,6))
         prod(:,199) = 0._r8
         prod(:,92) = 0._r8
         prod(:,96) = 0._r8
         prod(:,76) = 0._r8
         prod(:,112) = 0._r8
         prod(:,59) = 0._r8
         prod(:,113) = 0._r8
         prod(:,63) = 0._r8
         prod(:,93) = 0._r8
         prod(:,29) = + extfrc(:,6)
         prod(:,30) = + extfrc(:,7)
         prod(:,122) = 0._r8
         prod(:,102) = 0._r8
         prod(:,118) = 0._r8
         prod(:,182) = 0._r8
         prod(:,154) = + extfrc(:,4)
         prod(:,80) = 0._r8
         prod(:,31) = + extfrc(:,8)
         prod(:,32) = + extfrc(:,1)
         prod(:,33) = 0._r8
         prod(:,34) = 0._r8
         prod(:,35) = 0._r8
         prod(:,36) = 0._r8
         prod(:,37) = 0._r8
         prod(:,38) = 0._r8
         prod(:,39) = 0._r8
         prod(:,40) = 0._r8
         prod(:,41) = 0._r8
         prod(:,42) = 0._r8
         prod(:,43) = 0._r8
         prod(:,44) = 0._r8
         prod(:,45) = 0._r8
         prod(:,46) = 0._r8
         prod(:,47) = 0._r8
         prod(:,48) = 0._r8
         prod(:,50) = + extfrc(:,5)
         prod(:,68) = 0._r8
         prod(:,131) = 0._r8
         prod(:,132) = 0._r8
         prod(:,114) = 0._r8
         prod(:,165) = 0._r8
         prod(:,168) = 0._r8
         prod(:,137) = 0._r8
         prod(:,66) = 0._r8
         prod(:,69) = 0._r8
         prod(:,70) = 0._r8
         prod(:,139) = 0._r8
         prod(:,71) = 0._r8
         prod(:,107) = 0._r8
         prod(:,120) = 0._r8
         prod(:,161) = 0._r8
         prod(:,115) = 0._r8
         prod(:,108) = 0._r8
         prod(:,153) = 0._r8
         prod(:,150) = 0._r8
         prod(:,133) = 0._r8
         prod(:,188) = 0._r8
         prod(:,192) =rxt(:,87)*y(:,46) +rxt(:,89)*y(:,49) +rxt(:,39)*y(:,57)
         prod(:,144) = 0._r8
         prod(:,127) = 0._r8
         prod(:,81) = 0._r8
         prod(:,140) = 0._r8
         prod(:,198) = 0._r8
         prod(:,109) = 0._r8
         prod(:,183) = 0._r8
         prod(:,186) = 0._r8
         prod(:,185) = 0._r8
         prod(:,141) = 0._r8
         prod(:,187) = 0._r8
         prod(:,158) = 0._r8
         prod(:,136) = 0._r8
         prod(:,174) = 0._r8
         prod(:,195) =rxt(:,12)*y(:,122) +rxt(:,5)
         prod(:,196) =.330_r8*rxt(:,40)*y(:,57) + extfrc(:,15)
         prod(:,110) = 0._r8
         prod(:,145) = 0._r8
         prod(:,175) = 0._r8
         prod(:,173) = 0._r8
         prod(:,172) = 0._r8
         prod(:,142) = 0._r8
         prod(:,178) = 0._r8
         prod(:,143) = 0._r8
         prod(:,116) = 0._r8
         prod(:,206) =.050_r8*rxt(:,40)*y(:,57)
      end if
      end subroutine indprd
      end module mo_indprd
