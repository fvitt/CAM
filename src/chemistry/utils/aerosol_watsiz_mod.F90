module aerosol_watsiz_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use aero_model,     only: wetdep_lq
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc

  use modal_aero_calcsize,   only: modal_aero_calcsize_sub
  use modal_aero_wateruptake,only: modal_aero_wateruptake_dr

  implicit none

contains

  subroutine aerosol_watsiz_tend(state, ptend, dt, pbuf)
    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    !print*,'FVDBG.aerosol_watsiz_tend.. wetdep_lq: ',wetdep_lq

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_wetdep', lq=wetdep_lq)

    ! Do calculations of mode radius and water uptake if:
    ! 1) modal aerosols are affecting the climate, or
    ! 2) prognostic modal aerosols are enabled
    call modal_aero_calcsize_sub(state, ptend, dt, pbuf)
    call modal_aero_wateruptake_dr(state, pbuf)

  end subroutine aerosol_watsiz_tend

end module aerosol_watsiz_mod
