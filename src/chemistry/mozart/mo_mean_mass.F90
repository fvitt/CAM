
module mo_mean_mass

  implicit none

  private
  public :: set_mean_mass, init_mean_mass

  integer :: id_o2, id_o, id_h, id_n

contains

  subroutine init_mean_mass
    use mo_chem_utls, only : get_spc_ndx

    implicit none

    id_o2 = get_spc_ndx('O2')
    id_o  = get_spc_ndx('O')
    id_h  = get_spc_ndx('H')
    id_n  = get_spc_ndx('N')

  endsubroutine init_mean_mass

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

    implicit none

    !-----------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------
    integer, intent(in)   ::      ncol,lchnk
    real(r8), intent(in)  ::      mmr(:,:,:)           ! species concentrations (kg/kg)
    real(r8), intent(out) ::      mbar(:,:)            ! mean mass (g/mole)

    !-----------------------------------------------------------------
    !        ... Local variables
    !-----------------------------------------------------------------

    !-------------------------------------------
    !  Mean mass not fixed for WACCM-X
    !-------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       mbar(:ncol,:) = mbarv(:ncol,:,lchnk)
    else
       mbar(:ncol,:pver) = mwdry
    endif

  end subroutine set_mean_mass

end module mo_mean_mass
