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
         prod(:,1) = + extfrc(:,24)
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
         prod(:,15) =.100_r8*rxt(:,493)*y(:,135)*y(:,29)
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = (rxt(:,437)*y(:,62) +rxt(:,439)*y(:,87) +rxt(:,448)*y(:,62) + &
                 rxt(:,476)*y(:,50) +.500_r8*rxt(:,477)*y(:,51) + &
                 .800_r8*rxt(:,483)*y(:,74) +rxt(:,484)*y(:,75) + &
                 .500_r8*rxt(:,553)*y(:,109) +1.800_r8*rxt(:,718)*y(:,179))*y(:,251) &
                  + (2.000_r8*rxt(:,471)*y(:,226) +.900_r8*rxt(:,472)*y(:,227) + &
                 rxt(:,474)*y(:,124) +2.000_r8*rxt(:,540)*y(:,239) + &
                 rxt(:,574)*y(:,235) +rxt(:,618)*y(:,261))*y(:,226) &
                  + (.200_r8*rxt(:,493)*y(:,29) +.100_r8*rxt(:,554)*y(:,111) + &
                 .270_r8*rxt(:,688)*y(:,6) +.270_r8*rxt(:,693)*y(:,110))*y(:,135) &
                  + (rxt(:,541)*y(:,227) +.450_r8*rxt(:,542)*y(:,233) + &
                 2.000_r8*rxt(:,543)*y(:,239))*y(:,239) &
                  + (.500_r8*rxt(:,704)*y(:,227) +.900_r8*rxt(:,706)*y(:,124)) &
                 *y(:,256) +rxt(:,53)*y(:,51) +.400_r8*rxt(:,76)*y(:,140) +rxt(:,81) &
                 *y(:,175) +.800_r8*rxt(:,85)*y(:,179)
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) = 0._r8
         prod(:,23) = 0._r8
         prod(:,24) =rxt(:,238)*y(:,125)*y(:,112)
         prod(:,25) = 0._r8
         prod(:,26) = 0._r8
         prod(:,27) = 0._r8
         prod(:,28) = 0._r8
         prod(:,29) =rxt(:,745)*y(:,251)*y(:,120) +rxt(:,755)*y(:,121)
         prod(:,30) = (rxt(:,646)*y(:,228) +rxt(:,651)*y(:,238) +rxt(:,656)*y(:,240) + &
                 rxt(:,663)*y(:,142))*y(:,125) +.500_r8*rxt(:,553)*y(:,251)*y(:,109) &
                  +.200_r8*rxt(:,699)*y(:,246)*y(:,124) +.500_r8*rxt(:,715)*y(:,178) &
                 *y(:,126)
         prod(:,31) = (rxt(:,648)*y(:,228) +rxt(:,653)*y(:,238) +rxt(:,658)*y(:,240) + &
                 rxt(:,666)*y(:,142))*y(:,197) +.200_r8*rxt(:,721)*y(:,203)*y(:,124) &
                  +.500_r8*rxt(:,717)*y(:,200)*y(:,178) +.500_r8*rxt(:,562)*y(:,251) &
                 *y(:,192)
         prod(:,32) = 0._r8
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,147) = 0._r8
         prod(:,149) = 0._r8
         prod(:,1) = + extfrc(:,2)
         prod(:,2) = + extfrc(:,3)
         prod(:,180) = 0._r8
         prod(:,50) = 0._r8
         prod(:,108) = 0._r8
         prod(:,51) = 0._r8
         prod(:,111) = 0._r8
         prod(:,114) = 0._r8
         prod(:,87) = 0._r8
         prod(:,139) = 0._r8
         prod(:,91) = 0._r8
         prod(:,69) = 0._r8
         prod(:,116) = 0._r8
         prod(:,215) =rxt(:,105)*y(:,34) +rxt(:,106)*y(:,35) +2.000_r8*rxt(:,112) &
                 *y(:,41) +rxt(:,113)*y(:,43) +3.000_r8*rxt(:,116)*y(:,55) &
                  +2.000_r8*rxt(:,124)*y(:,78)
         prod(:,70) = 0._r8
         prod(:,231) = 0._r8
         prod(:,132) = 0._r8
         prod(:,71) = 0._r8
         prod(:,94) = 0._r8
         prod(:,88) = 0._r8
         prod(:,131) = 0._r8
         prod(:,79) = 0._r8
         prod(:,95) = 0._r8
         prod(:,89) = 0._r8
         prod(:,190) = 0._r8
         prod(:,110) = 0._r8
         prod(:,40) = 0._r8
         prod(:,80) = 0._r8
         prod(:,218) =.180_r8*rxt(:,56)*y(:,54)
         prod(:,193) = 0._r8
         prod(:,38) = 0._r8
         prod(:,183) = 0._r8
         prod(:,206) = 0._r8
         prod(:,126) = 0._r8
         prod(:,123) = 0._r8
         prod(:,162) = 0._r8
         prod(:,112) = 0._r8
         prod(:,222) =4.000_r8*rxt(:,104)*y(:,33) +rxt(:,105)*y(:,34) &
                  +2.000_r8*rxt(:,107)*y(:,36) +2.000_r8*rxt(:,108)*y(:,37) &
                  +2.000_r8*rxt(:,109)*y(:,38) +rxt(:,110)*y(:,39) &
                  +2.000_r8*rxt(:,111)*y(:,40) +3.000_r8*rxt(:,114)*y(:,44) &
                  +rxt(:,115)*y(:,46) +rxt(:,126)*y(:,82) +rxt(:,127)*y(:,83) &
                  +rxt(:,128)*y(:,84)
         prod(:,53) = 0._r8
         prod(:,36) = 0._r8
         prod(:,220) = 0._r8
         prod(:,177) = 0._r8
         prod(:,196) = (rxt(:,57) +rxt(:,139))*y(:,63) +.380_r8*rxt(:,56)*y(:,54) &
                  + extfrc(:,1)
         prod(:,62) =rxt(:,105)*y(:,34) +rxt(:,106)*y(:,35) +rxt(:,108)*y(:,37) &
                  +2.000_r8*rxt(:,109)*y(:,38) +2.000_r8*rxt(:,110)*y(:,39) &
                  +rxt(:,111)*y(:,40) +2.000_r8*rxt(:,124)*y(:,78) +rxt(:,127)*y(:,83) &
                  +rxt(:,128)*y(:,84)
         prod(:,77) =rxt(:,107)*y(:,36) +rxt(:,108)*y(:,37) +rxt(:,126)*y(:,82)
         prod(:,55) = 0._r8
         prod(:,104) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,42) = 0._r8
         prod(:,170) =rxt(:,106)*y(:,35) +rxt(:,110)*y(:,39)
         prod(:,192) = 0._r8
         prod(:,174) = 0._r8
         prod(:,230) = (rxt(:,55) +.330_r8*rxt(:,56))*y(:,54)
         prod(:,202) =1.440_r8*rxt(:,56)*y(:,54)
         prod(:,140) = 0._r8
         prod(:,43) = 0._r8
         prod(:,168) = 0._r8
         prod(:,214) = 0._r8
         prod(:,75) = 0._r8
         prod(:,161) = 0._r8
         prod(:,76) = 0._r8
         prod(:,229) = 0._r8
         prod(:,117) = 0._r8
         prod(:,156) = 0._r8
         prod(:,167) = 0._r8
         prod(:,194) = 0._r8
         prod(:,67) = 0._r8
         prod(:,197) = 0._r8
         prod(:,99) = 0._r8
         prod(:,44) = 0._r8
         prod(:,178) = 0._r8
         prod(:,141) = 0._r8
         prod(:,135) = 0._r8
         prod(:,204) = 0._r8
         prod(:,109) = 0._r8
         prod(:,151) = 0._r8
         prod(:,34) = 0._r8
         prod(:,205) = 0._r8
         prod(:,97) = 0._r8
         prod(:,124) = 0._r8
         prod(:,101) = 0._r8
         prod(:,127) = 0._r8
         prod(:,181) = 0._r8
         prod(:,210) = 0._r8
         prod(:,176) = (.800_r8*rxt(:,141) +rxt(:,144) +rxt(:,145) + &
                 .800_r8*rxt(:,147)) + extfrc(:,18)
         prod(:,86) = 0._r8
         prod(:,105) = 0._r8
         prod(:,191) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,37) = 0._r8
         prod(:,9) = 0._r8
         prod(:,227) = + extfrc(:,8)
         prod(:,221) = + extfrc(:,7)
         prod(:,226) = 0._r8
         prod(:,172) = 0._r8
         prod(:,106) = 0._r8
         prod(:,10) = + extfrc(:,9)
         prod(:,11) = + extfrc(:,10)
         prod(:,12) = 0._r8
         prod(:,13) = + extfrc(:,11)
         prod(:,228) = (rxt(:,57) +rxt(:,139))*y(:,63) +.180_r8*rxt(:,56)*y(:,54)
         prod(:,217) = 0._r8
         prod(:,224) = 0._r8
         prod(:,92) = 0._r8
         prod(:,102) = 0._r8
         prod(:,68) = 0._r8
         prod(:,121) = 0._r8
         prod(:,45) = 0._r8
         prod(:,143) = 0._r8
         prod(:,52) = 0._r8
         prod(:,93) = 0._r8
         prod(:,14) = + extfrc(:,12)
         prod(:,15) = + extfrc(:,4)
         prod(:,128) = 0._r8
         prod(:,107) = 0._r8
         prod(:,154) = 0._r8
         prod(:,213) = 0._r8
         prod(:,175) = + extfrc(:,5)
         prod(:,78) = 0._r8
         prod(:,16) = + extfrc(:,13)
         prod(:,17) = + extfrc(:,14)
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) = 0._r8
         prod(:,23) = 0._r8
         prod(:,24) = 0._r8
         prod(:,25) = 0._r8
         prod(:,26) = 0._r8
         prod(:,27) = 0._r8
         prod(:,28) = 0._r8
         prod(:,29) = 0._r8
         prod(:,30) = 0._r8
         prod(:,31) = 0._r8
         prod(:,32) = 0._r8
         prod(:,33) = 0._r8
         prod(:,35) = + extfrc(:,6)
         prod(:,56) = 0._r8
         prod(:,137) = 0._r8
         prod(:,146) = 0._r8
         prod(:,118) = 0._r8
         prod(:,189) = 0._r8
         prod(:,195) = 0._r8
         prod(:,138) = 0._r8
         prod(:,54) = 0._r8
         prod(:,41) = 0._r8
         prod(:,81) = 0._r8
         prod(:,119) = 0._r8
         prod(:,184) = 0._r8
         prod(:,57) = 0._r8
         prod(:,96) = 0._r8
         prod(:,82) = 0._r8
         prod(:,83) = 0._r8
         prod(:,103) = 0._r8
         prod(:,58) = 0._r8
         prod(:,59) = 0._r8
         prod(:,60) = + extfrc(:,20)
         prod(:,61) = 0._r8
         prod(:,98) = 0._r8
         prod(:,185) = 0._r8
         prod(:,187) = 0._r8
         prod(:,72) = 0._r8
         prod(:,90) = 0._r8
         prod(:,186) = 0._r8
         prod(:,46) = 0._r8
         prod(:,73) = 0._r8
         prod(:,120) = 0._r8
         prod(:,84) = 0._r8
         prod(:,130) = 0._r8
         prod(:,129) = 0._r8
         prod(:,74) = 0._r8
         prod(:,63) = 0._r8
         prod(:,64) = 0._r8
         prod(:,39) = 0._r8
         prod(:,100) = 0._r8
         prod(:,65) = 0._r8
         prod(:,152) = 0._r8
         prod(:,66) = 0._r8
         prod(:,113) = 0._r8
         prod(:,150) = 0._r8
         prod(:,182) = 0._r8
         prod(:,142) = 0._r8
         prod(:,133) = 0._r8
         prod(:,188) = 0._r8
         prod(:,173) = 0._r8
         prod(:,157) = 0._r8
         prod(:,212) = 0._r8
         prod(:,216) =rxt(:,113)*y(:,43) +rxt(:,115)*y(:,46) +rxt(:,55)*y(:,54)
         prod(:,166) = 0._r8
         prod(:,160) = (rxt(:,142) +rxt(:,143) +rxt(:,144) +rxt(:,145) +rxt(:,146) + &
                 rxt(:,148)) + extfrc(:,15)
         prod(:,153) = 0._r8
         prod(:,115) = 0._r8
         prod(:,171) = 0._r8
         prod(:,223) = 0._r8
         prod(:,134) = 0._r8
         prod(:,208) = 0._r8
         prod(:,207) = 0._r8
         prod(:,209) = 0._r8
         prod(:,164) = 0._r8
         prod(:,211) = 0._r8
         prod(:,179) = 0._r8
         prod(:,155) = 0._r8
         prod(:,125) = (1.200_r8*rxt(:,141) +rxt(:,142) +rxt(:,146) + &
                 1.200_r8*rxt(:,147)) + extfrc(:,19)
         prod(:,148) = (rxt(:,143) +rxt(:,148)) + extfrc(:,22)
         prod(:,159) = 0._r8
         prod(:,122) = (rxt(:,142) +rxt(:,144) +rxt(:,145) +rxt(:,146))
         prod(:,200) = 0._r8
         prod(:,225) =rxt(:,14)*y(:,113)
         prod(:,47) = 0._r8
         prod(:,48) = 0._r8
         prod(:,158) = + extfrc(:,23)
         prod(:,219) =.330_r8*rxt(:,56)*y(:,54) + extfrc(:,16)
         prod(:,144) = + extfrc(:,17)
         prod(:,136) = 0._r8
         prod(:,169) = 0._r8
         prod(:,201) = 0._r8
         prod(:,199) = 0._r8
         prod(:,198) = 0._r8
         prod(:,163) = 0._r8
         prod(:,49) = + extfrc(:,21)
         prod(:,85) = 0._r8
         prod(:,203) = 0._r8
         prod(:,165) = 0._r8
         prod(:,145) = 0._r8
         prod(:,232) =.050_r8*rxt(:,56)*y(:,54)
      end if
      end subroutine indprd
      end module mo_indprd
