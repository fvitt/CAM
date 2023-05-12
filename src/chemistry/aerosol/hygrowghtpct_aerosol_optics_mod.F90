module hygrowghtpct_aerosol_optics_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_optics_mod, only: aerosol_optics
  use aerosol_state_mod, only: aerosol_state
  use aerosol_properties_mod, only: aerosol_properties
  use table_interp_mod, only: table_interp

  implicit none

  private
  public :: hygrowghtpct_aerosol_optics

  type, extends(aerosol_optics) :: hygrowghtpct_aerosol_optics

     real(r8), allocatable :: totalmmr(:,:)
     real(r8), allocatable :: wgtpct(:,:)

     real(r8), pointer :: tbl_wgtpct(:)

     real(r8), pointer :: sw_hygro_ext_wtp(:,:)
     real(r8), pointer :: sw_hygro_ssa_wtp(:,:)
     real(r8), pointer :: sw_hygro_asm_wtp(:,:)
     real(r8), pointer :: lw_hygro_abs_wtp(:,:)

   contains

     procedure :: sw_props
     procedure :: lw_props

     final :: destructor

  end type hygrowghtpct_aerosol_optics

  interface hygrowghtpct_aerosol_optics
     procedure :: constructor
  end interface hygrowghtpct_aerosol_optics

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(aero_props, aero_state, ilist, ibin, ncol, nlev) result(newobj)

    class(aerosol_properties),intent(in) :: aero_props
    class(aerosol_state),intent(in) :: aero_state
    integer, intent(in) :: ilist
    integer, intent(in) :: ibin
    integer, intent(in) :: ncol, nlev


    type(hygrowghtpct_aerosol_optics), pointer :: newobj

    integer :: ierr, nspec
    integer :: ispec

    real(r8), pointer :: specmmr(:,:)   ! species mass mixing ratio
    real(r8), pointer :: wgtpct_in(:,:) !  weight precent of H2SO4/H2O solution

    allocate(newobj, stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%totalmmr(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%wgtpct(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    wgtpct_in => aero_state%wgtpct()

    call aero_props%optics_params(ilist, ibin, wgtpct=newobj%tbl_wgtpct)

    newobj%wgtpct(:ncol,:) = wgtpct_in(:ncol,:)
    where ( newobj%wgtpct < minval(newobj%tbl_wgtpct) )
       newobj%wgtpct = minval(newobj%tbl_wgtpct)
    end where
    where ( newobj%wgtpct > maxval(newobj%tbl_wgtpct) )
       newobj%wgtpct = maxval(newobj%tbl_wgtpct)
    end where

    nspec = aero_props%nspecies(ilist, ibin)

    newobj%totalmmr(:,:) = 0._r8

    do ispec = 1,nspec

       call aero_state%get_ambient_mmr(ilist,ispec,ibin,specmmr)
       newobj%totalmmr(:ncol,:nlev) = newobj%totalmmr(:ncol,:nlev) + specmmr(:ncol,:nlev)

    end do

    call aero_props%optics_params(ilist, ibin,  &
         sw_hygro_ext_wtp=newobj%sw_hygro_ext_wtp, &
         sw_hygro_ssa_wtp=newobj%sw_hygro_ssa_wtp, &
         sw_hygro_asm_wtp=newobj%sw_hygro_asm_wtp, &
         lw_hygro_ext_wtp=newobj%lw_hygro_abs_wtp)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine sw_props(self, ncol, ilev, iwav, pext, pabs, palb, pasm)

    class(hygrowghtpct_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol
    integer, intent(in) :: ilev
    integer, intent(in) :: iwav
    real(r8),intent(out) :: pext(ncol)
    real(r8),intent(out) :: pabs(ncol)
    real(r8),intent(out) :: palb(ncol)
    real(r8),intent(out) :: pasm(ncol)

    integer :: icol

    do icol = 1, ncol

       pext(icol) = table_interp( self%tbl_wgtpct, self%sw_hygro_ext_wtp(:,iwav), self%wgtpct(icol,ilev) )
       pabs(icol) = (1._r8 - table_interp( self%tbl_wgtpct, self%sw_hygro_ssa_wtp(:,iwav), self%wgtpct(icol,ilev) ) ) * pext(icol)
       pasm(icol) = table_interp( self%tbl_wgtpct, self%sw_hygro_asm_wtp(:,iwav), self%wgtpct(icol,ilev) )

       pext(icol) = pext(icol)*self%totalmmr(icol,ilev)
       pabs(icol) = pabs(icol)*self%totalmmr(icol,ilev)
       pabs(icol) = max(0._r8,pabs(icol))
       pabs(icol) = min(pext(icol),pabs(icol))

       palb(icol) = 1._r8-pabs(icol)/max(pext(icol),1.e-40_r8)
       palb(icol) = 1._r8-pabs(icol)/max(pext(icol),1.e-40_r8)

    end do

  end subroutine sw_props

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine lw_props(self, ncol, ilev, iwav, pabs)

    class(hygrowghtpct_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol
    integer, intent(in) :: ilev
    integer, intent(in) :: iwav
    real(r8),intent(out) :: pabs(ncol)

    integer :: icol

    do icol = 1, ncol

       pabs(icol) = table_interp( self%tbl_wgtpct, self%lw_hygro_abs_wtp(:,iwav), self%wgtpct(icol,ilev) )
       pabs(icol) = pabs(icol)*self%totalmmr(icol,ilev)
       pabs(icol) = max(0._r8,pabs(icol))

    end do

  end subroutine lw_props

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)

    type(hygrowghtpct_aerosol_optics), intent(inout) :: self

    deallocate(self%totalmmr)
    deallocate(self%wgtpct)

    nullify(self%tbl_wgtpct)
    nullify(self%sw_hygro_ext_wtp)
    nullify(self%sw_hygro_ssa_wtp)
    nullify(self%sw_hygro_asm_wtp)
    nullify(self%lw_hygro_abs_wtp)

  end subroutine destructor

end module hygrowghtpct_aerosol_optics_mod
