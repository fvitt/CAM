module zmean_fields
  use shr_kind_mod, only: r8=>SHR_KIND_R8
  use ppgrid,       only: begchunk, endchunk, pcols, pver
  use Zonal_Mean,   only: ZonalMean_t

  implicit none

  type(ZonalMean_t) :: ZMobj

  integer, parameter :: ZonalNbasis = 12

contains

  subroutine zmean_fields_init
     call ZMobj%init(ZonalNbasis)
  end subroutine zmean_fields_init

  function zmean_3d( fld ) result(fldzm)
    real(r8), intent(in) :: fld(pcols,pver,begchunk:endchunk)

    real(r8) :: fldzm(pcols,pver,begchunk:endchunk)

    real(r8) :: Zonal_Bamp3d(ZonalNbasis,pver)

    call ZMobj%calc_amps(fld,Zonal_Bamp3d)
    call ZMobj%eval_grid(Zonal_Bamp3d,fldzm)

  end function zmean_3d

end module zmean_fields
