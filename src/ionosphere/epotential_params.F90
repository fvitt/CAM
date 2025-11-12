module epotential_params
  use shr_kind_mod, only: r8 => shr_kind_r8

  logical,  public :: epot_active = .false.
  real(r8), public :: epot_crit_colats(2) = -huge(1._r8)

  logical, public :: edyn3d_active = .false.
  integer, public :: edyn3d_nmlat_h = 0
  integer, public :: edyn3d_nmlon = 0
  integer, public :: edyn3d_nhgt = 0

end module epotential_params
