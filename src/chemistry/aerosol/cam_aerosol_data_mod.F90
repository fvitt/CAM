module cam_aerosol_data_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physics_buffer, only: physics_buffer_desc
  use physics_types, only: physics_state
  use aerosol_data_mod, only: aerosol_data
  use constituents,   only: pcnst

  implicit none

  type, abstract, extends(aerosol_data) :: cam_aerosol_data
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
     logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
  contains
  end type cam_aerosol_data

contains
end module cam_aerosol_data_mod
