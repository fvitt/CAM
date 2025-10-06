!-----------------------------------------------------------------------------
! High-latitude dynamo inputs
!-----------------------------------------------------------------------------
module high_lat_pot_mod

  use shr_kind_mod, only: r8 => shr_kind_r8 ! 8-byte reals

  real(r8), allocatable, public :: phihm(:,:) ! high-latitude potential
  real(r8), allocatable, public :: fachm(:,:) ! high-latitude field-aligned current

 contains

  !-----------------------------------------------------------------------------
  ! high-latitude electric potential
  !-----------------------------------------------------------------------------
  subroutine edyn3d_highlat_potential_get(hl_pot)

    use mpi_module, only: mlond0, mlond1, mlatd0, mlatd1, mlon0, mlon1
    use params_module, only: nmlat_h, nmlat_T1, nmlon

    real(r8), intent(out) :: hl_pot(2,mlatd0:mlatd1,mlond0:mlond1)

    integer :: h,i,j,jj

    do h = 1,2
       do i = max(1,mlond0), min(mlond1,nmlon)
          do j = max(1,mlatd0),min(mlatd1,nmlat_h)
             if (h==1) then
                jj = j
             else
                jj = nmlat_T1 - j + 1
             end if
             hl_pot(h,j,i) = phihm(i,jj)

             ! wrap around longitude points
             if (i==1) then
                hl_pot(h,j,0) = phihm(nmlon,jj)
             else if (i==nmlon) then
                hl_pot(h,j,nmlon+1) = phihm(1,jj)
             end if

          end do
       end do
    end do

  end subroutine edyn3d_highlat_potential_get

  !-----------------------------------------------------------------------------
  ! high-latitude field-aligned currents
  !-----------------------------------------------------------------------------
  subroutine edyn3d_highlat_currents_get(hl_fac)

    use mpi_module, only: mlond0, mlond1, mlatd0, mlatd1, mlon0, mlon1
    use params_module, only: nmlat_h, nmlat_T1, nmlon

    real(r8), intent(out) :: hl_fac(2,mlatd0:mlatd1,mlond0:mlond1)

    integer :: h,i,j,jj

    do h = 1,2
       do i = max(1,mlond0), min(mlond1,nmlon)
          do j = max(1,mlatd0),min(mlatd1,nmlat_h)
             if (h==1) then
                jj = j
             else
                jj = nmlat_T1 - j + 1
             end if
             hl_fac(h,j,i) = fachm(i,jj)*1.e-6_r8 ! convert uA/m2 to A/m2

             ! wrap around longitude points
             if (i==1) then
                hl_fac(h,j,0) = fachm(nmlon,jj)*1.e-6_r8
             else if (i==nmlon) then
                hl_fac(h,j,nmlon+1) = fachm(1,jj)*1.e-6_r8
             end if

          end do
       end do
    end do

  end subroutine edyn3d_highlat_currents_get

end module
