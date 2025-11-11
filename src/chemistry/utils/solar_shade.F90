module solar_shade
  use shr_kind_mod, only: r8 => shr_kind_r8
  use solar_irrad_data, only: nbins, we
  use physics_types,only : physics_state
  use ppgrid, only : pcols, begchunk, endchunk

  implicit none

  real(r8), public, protected, allocatable :: sun_shade(:,:,:) ! chunk, col, wavelen dependent

contains

  subroutine solar_shade_init()

    allocate(sun_shade(nbins, pcols, begchunk:endchunk))
    sun_shade = 1._r8

  end subroutine solar_shade_init

end module solar_shade
