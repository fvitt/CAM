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
         prod(:,1) = 0._r8
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
         prod(:,15) = 0._r8
         prod(:,16) = 0._r8
         prod(:,17) = (rxt(:,362)*y(:,109) +rxt(:,366)*y(:,109))*y(:,29)
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) =rxt(:,167)*y(:,60)*y(:,53)
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,1) = + extfrc(:,3)
         prod(:,2) = + extfrc(:,12)
         prod(:,34) = 0._r8
         prod(:,88) = 0._r8
         prod(:,51) = 0._r8
         prod(:,80) =.180_r8*rxt(:,38)*y(:,22)
         prod(:,75) =rxt(:,53)*y(:,17) +rxt(:,55)*y(:,19) +rxt(:,37)*y(:,22)
         prod(:,45) = 0._r8
         prod(:,27) = 0._r8
         prod(:,21) = 0._r8
         prod(:,77) = 0._r8
         prod(:,69) = 0._r8
         prod(:,52) = (rxt(:,39) +rxt(:,78))*y(:,30) +.380_r8*rxt(:,38)*y(:,22) &
                  + extfrc(:,1)
         prod(:,28) =rxt(:,45)*y(:,8) +rxt(:,46)*y(:,9) +rxt(:,48)*y(:,11) &
                  +2.000_r8*rxt(:,49)*y(:,12) +2.000_r8*rxt(:,50)*y(:,13) +rxt(:,51) &
                 *y(:,14) +2.000_r8*rxt(:,64)*y(:,40) +rxt(:,67)*y(:,45) +rxt(:,68) &
                 *y(:,46)
         prod(:,33) =rxt(:,47)*y(:,10) +rxt(:,48)*y(:,11) +rxt(:,66)*y(:,44)
         prod(:,44) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,60) =rxt(:,46)*y(:,9) +rxt(:,50)*y(:,13)
         prod(:,72) = (rxt(:,37) +.330_r8*rxt(:,38))*y(:,22)
         prod(:,85) =1.440_r8*rxt(:,38)*y(:,22)
         prod(:,56) = 0._r8
         prod(:,22) = 0._r8
         prod(:,67) = 0._r8
         prod(:,74) = 0._r8
         prod(:,32) = 0._r8
         prod(:,70) = 0._r8
         prod(:,48) = 0._r8
         prod(:,61) = 0._r8
         prod(:,66) = 0._r8
         prod(:,65) = (.800_r8*rxt(:,81) +rxt(:,83) +.800_r8*rxt(:,85) +rxt(:,87)) &
                  + extfrc(:,19)
         prod(:,39) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,86) = + extfrc(:,4)
         prod(:,89) = + extfrc(:,2)
         prod(:,76) = 0._r8
         prod(:,9) = + extfrc(:,8)
         prod(:,10) = + extfrc(:,9)
         prod(:,11) = 0._r8
         prod(:,12) = + extfrc(:,10)
         prod(:,79) = (rxt(:,39) +rxt(:,78))*y(:,30) +.180_r8*rxt(:,38)*y(:,22)
         prod(:,82) = 0._r8
         prod(:,84) = 0._r8
         prod(:,40) = 0._r8
         prod(:,43) = 0._r8
         prod(:,13) = + extfrc(:,5)
         prod(:,14) = + extfrc(:,11)
         prod(:,58) = 0._r8
         prod(:,71) = 0._r8
         prod(:,68) = + extfrc(:,13)
         prod(:,36) = 0._r8
         prod(:,15) = + extfrc(:,6)
         prod(:,16) = + extfrc(:,7)
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,37) = 0._r8
         prod(:,47) = 0._r8
         prod(:,53) = 0._r8
         prod(:,29) = 0._r8
         prod(:,30) = 0._r8
         prod(:,23) = 0._r8
         prod(:,46) = 0._r8
         prod(:,55) = 0._r8
         prod(:,35) = 0._r8
         prod(:,31) = 0._r8
         prod(:,54) = 0._r8
         prod(:,24) = 0._r8
         prod(:,42) = 0._r8
         prod(:,41) = 0._r8
         prod(:,73) =rxt(:,45)*y(:,8) +rxt(:,46)*y(:,9) +2.000_r8*rxt(:,52)*y(:,15) &
                  +rxt(:,53)*y(:,17) +3.000_r8*rxt(:,56)*y(:,23) +2.000_r8*rxt(:,64) &
                 *y(:,40)
         prod(:,81) =4.000_r8*rxt(:,44)*y(:,7) +rxt(:,45)*y(:,8) +2.000_r8*rxt(:,47) &
                 *y(:,10) +2.000_r8*rxt(:,48)*y(:,11) +2.000_r8*rxt(:,49)*y(:,12) &
                  +rxt(:,50)*y(:,13) +2.000_r8*rxt(:,51)*y(:,14) +3.000_r8*rxt(:,54) &
                 *y(:,18) +rxt(:,55)*y(:,19) +rxt(:,66)*y(:,44) +rxt(:,67)*y(:,45) &
                  +rxt(:,68)*y(:,46)
         prod(:,64) = (rxt(:,80) +rxt(:,82) +rxt(:,83) +rxt(:,84) +rxt(:,86) + &
                 rxt(:,87)) + extfrc(:,21)
         prod(:,78) = 0._r8
         prod(:,50) = (rxt(:,80) +1.200_r8*rxt(:,81) +1.200_r8*rxt(:,85) +rxt(:,86)) &
                  + extfrc(:,17)
         prod(:,57) = (rxt(:,82) +rxt(:,84)) + extfrc(:,15)
         prod(:,63) = 0._r8
         prod(:,49) = (rxt(:,80) +rxt(:,83) +rxt(:,86) +rxt(:,87)) + extfrc(:,16)
         prod(:,83) =rxt(:,14)*y(:,54)
         prod(:,25) = 0._r8
         prod(:,26) = 0._r8
         prod(:,62) = + extfrc(:,14)
         prod(:,87) =.330_r8*rxt(:,38)*y(:,22) + extfrc(:,18)
         prod(:,59) = + extfrc(:,20)
         prod(:,38) = 0._r8
         prod(:,90) =.050_r8*rxt(:,38)*y(:,22)
      end if
      end subroutine indprd
      end module mo_indprd
