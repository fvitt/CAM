module mo_mean_mass

  implicit none

  private
  public :: set_mean_mass

contains

  subroutine set_mean_mass( ncol, lchnk, mmr, mbar )
    !-----------------------------------------------------------------
    !        ... Set the invariant densities (molecules/cm**3)
    !-----------------------------------------------------------------

    use shr_kind_mod,     only : r8 => shr_kind_r8
    use ppgrid,           only : pver, pcols
    use chem_mods,        only : adv_mass, gas_pcnst
    use physconst,        only : mwdry                   ! molecular weight of dry air
    use cam_abortutils,   only : endrun
    use phys_control,     only : waccmx_is               !WACCM-X runtime switch
    use air_composition,  only : mbarv

    !-----------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------
    integer, intent(in)   ::      ncol,lchnk
    real(r8), intent(in)  ::      mmr(:,:,:)           ! species concentrations (kg/kg)
    real(r8), intent(out) ::      mbar(:,:)            ! mean mass (g/mole)

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       mbar(:ncol,:) = mbarv(:ncol,:,lchnk)
    else
       mbar(:ncol,:pver) = mwdry
    endif

  end subroutine set_mean_mass

end module mo_mean_mass
