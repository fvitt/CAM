module hygrocoreshell_aerosol_optics_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_optics_mod, only: aerosol_optics
  use aerosol_state_mod, only: aerosol_state
  use aerosol_properties_mod, only: aerosol_properties
  use table_interp_mod, only: table_interp

  implicit none

  private
  public :: hygrocoreshell_aerosol_optics

  !> hygrocoreshell_aerosol_optics
  !! Table look up implementation of aerosol_optics to parameterize aerosol
  !! radiative properties in terms of core mass fraction, black carbon/dust fraction,
  !! kappa and relative humidity
  type, extends(aerosol_optics) :: hygrocoreshell_aerosol_optics

     real(r8), allocatable :: totalmmr(:,:) ! total mmr of the aerosol
     real(r8), allocatable :: corefrac(:,:) ! mass fraction that is core
     real(r8), allocatable :: bcdust(:,:)   ! mass fraction of bc vs (bc + dust)
     real(r8), allocatable :: kappa(:,:)    ! hygroscopicity
     real(r8), allocatable :: relh(:,:)     ! relative humidity

     real(r8), pointer :: sw_hygro_coreshell_ext(:,:,:,:,:) => null() ! short wave extinction table
     real(r8), pointer :: sw_hygro_coreshell_ssa(:,:,:,:,:) => null() ! short wave single-scatter albedo table
     real(r8), pointer :: sw_hygro_coreshell_asm(:,:,:,:,:) => null() ! short wave asymmetry table
     real(r8), pointer :: lw_hygro_coreshell_abs(:,:,:,:,:) => null() ! long wave absorption table

     real(r8), pointer :: tbl_corefrac(:) => null() ! core fraction dimension values
     real(r8), pointer :: tbl_bcdust(:) => null()   ! bc/(bc + dust) fraction dimension values
     real(r8), pointer :: tbl_kap(:) => null()      ! hygroscopicity dimension values
     real(r8), pointer :: tbl_relh(:) => null()     ! relative humidity dimension values

   contains

     procedure :: sw_props
     procedure :: lw_props

     final :: destructor

  end type hygrocoreshell_aerosol_optics

  interface hygrocoreshell_aerosol_optics
     procedure :: constructor
  end interface hygrocoreshell_aerosol_optics

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
 function constructor(aero_props, aero_state, ilist, ibin, ncol, nlev, relhum) result(newobj)

    class(aerosol_properties),intent(in) :: aero_props ! aerosol_properties object
    class(aerosol_state),intent(in) :: aero_state      ! aerosol_state object
    integer, intent(in) :: ilist  ! climate or a diagnostic list number
    integer, intent(in) :: ibin   ! bin number
    integer, intent(in) :: ncol   ! number of columns
    integer, intent(in) :: nlev   ! number of levels
    real(r8),intent(in) :: relhum(ncol,nlev) ! relative humidity

    type(hygrocoreshell_aerosol_optics), pointer :: newobj

    integer :: ierr, nspec
    integer :: ilev, ispec, icol

    real(r8), pointer :: kappa_ptr(:,:) ! hygroscopicity
    real(r8), pointer :: specmmr(:,:)   ! species mass mixing ratio

    real(r8) :: coremmr(ncol,nlev)
    real(r8) :: coredustmmr(ncol,nlev)
    real(r8) :: corebcmmr(ncol,nlev)
    real(r8) :: shellmmr(ncol,nlev)
    real(r8) :: bcdustmmr(ncol,nlev)

    character(len=32) :: spectype  ! species type
    character(len=32) :: specmorph
    real(r8)          :: specdens  ! species density (kg/m3)

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

    allocate(newobj%corefrac(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%bcdust(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%kappa(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%relh(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    nspec = aero_props%nspecies(ilist,ibin)

    coremmr(:,:) = 0._r8
    coredustmmr(:,:) = 0._r8
    corebcmmr(:,:) = 0._r8
    shellmmr(:,:) = 0._r8

    do ispec = 1,nspec

       call aero_state%get_ambient_mmr(ilist,ispec,ibin,specmmr)

       call aero_props%get(ibin, ispec, list_ndx=ilist, density=specdens, &
                           spectype=spectype, specmorph=specmorph)

       if (trim(specmorph) == 'core') then
          if (trim(spectype) == 'dust') then
             coredustmmr(:ncol,:nlev) = coredustmmr(:ncol,:nlev) + specmmr(:ncol,:nlev)
          end if
          if (trim(spectype) == 'black-c') then
             corebcmmr(:ncol,:nlev) = corebcmmr(:ncol,:nlev) + specmmr(:ncol,:nlev)
          end if
          coremmr(:ncol,:nlev) = coremmr(:ncol,:nlev) + specmmr(:ncol,:nlev)
       else if (trim(specmorph) == 'shell') then
          shellmmr(:ncol,:nlev) = shellmmr(:ncol,:nlev) + specmmr(:ncol,:nlev)
       else
          nullify(newobj)
          return
       end if

    end do

    newobj%totalmmr(:,:) = coremmr(:,:) + shellmmr(:,:)
    bcdustmmr(:,:) = corebcmmr(:,:) + coredustmmr(:,:)

    do ilev = 1, nlev
       do icol = 1, ncol

          if (newobj%totalmmr(icol,ilev) > 0._r8) then
             newobj%corefrac(icol,ilev) = coremmr(icol,ilev) / newobj%totalmmr(icol,ilev)
          else
             newobj%corefrac(icol,ilev) = 0._r8
          end if
          newobj%corefrac(icol,ilev) = max(0._r8, min(1.0_r8, newobj%corefrac(icol,ilev)))

          if (bcdustmmr(icol,ilev) > 0._r8) then
             newobj%bcdust(icol,ilev) = corebcmmr(icol,ilev) / bcdustmmr(icol,ilev)
          else
             newobj%bcdust(icol,ilev) = 0._r8
          end if
          newobj%bcdust(icol,ilev) = max(0._r8, min(1.0_r8, newobj%bcdust(icol,ilev)))

       end do
    end do

    kappa_ptr => aero_state%hygroscopicity(ilist, ibin)
    newobj%kappa(:ncol,:) = kappa_ptr(:ncol,:)

    call aero_props%optics_params(ilist, ibin, &
         corefrac=newobj%tbl_corefrac, kap=newobj%tbl_kap, &
         bcdust=newobj%tbl_bcdust, relh=newobj%tbl_relh)


    where ( newobj%kappa < minval(newobj%tbl_kap) )
       newobj%kappa = minval(newobj%tbl_kap)
    end where
    where ( newobj%kappa > maxval(newobj%tbl_kap) )
       newobj%kappa = maxval(newobj%tbl_kap)
    end where

    newobj%relh(:ncol,:) = relhum(:ncol,:)

    where ( newobj%relh < minval(newobj%tbl_relh) )
       newobj%relh = minval(newobj%tbl_relh)
    end where
    where ( newobj%relh > maxval(newobj%tbl_relh) )
       newobj%relh = maxval(newobj%tbl_relh)
    end where

    where (newobj%corefrac < minval(newobj%tbl_corefrac))
       newobj%corefrac = minval(newobj%tbl_corefrac)
    end where
    where (newobj%corefrac > maxval(newobj%tbl_corefrac))
       newobj%corefrac = maxval(newobj%tbl_corefrac)
    end where

    where (newobj%bcdust < minval(newobj%tbl_bcdust))
       newobj%bcdust = minval(newobj%tbl_bcdust)
    end where
    where (newobj%bcdust > maxval(newobj%tbl_bcdust))
       newobj%bcdust = maxval(newobj%tbl_bcdust)
    end where

    ! long wave optical properties table
    call aero_props%optics_params(ilist, ibin,  &
         sw_hygro_coreshell_ext=newobj%sw_hygro_coreshell_ext, &
         sw_hygro_coreshell_ssa=newobj%sw_hygro_coreshell_ssa, &
         sw_hygro_coreshell_asm=newobj%sw_hygro_coreshell_asm, &
         lw_hygro_coreshell_ext=newobj%lw_hygro_coreshell_abs)

  end function constructor

  !------------------------------------------------------------------------------
  ! returns short wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine sw_props(self, ncol, ilev, iwav, pext, pabs, palb, pasm)

    class(hygrocoreshell_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pext(ncol) ! parameterized specific extinction (m2/kg)
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)
    real(r8),intent(out) :: palb(ncol) ! parameterized asymmetry factor
    real(r8),intent(out) :: pasm(ncol) ! parameterized single scattering albedo

    integer :: icol

    do icol = 1, ncol

       pext(icol) = table_interp( self%tbl_relh, self%tbl_corefrac, self%tbl_bcdust, self%tbl_kap, &
            self%sw_hygro_coreshell_ext(:,iwav,:,:,:), &
            self%relh(icol,ilev), self%corefrac(icol,ilev), self%bcdust(icol,ilev), self%kappa(icol,ilev) )

       pabs(icol) = (1._r8 - table_interp( self%tbl_relh, self%tbl_corefrac, self%tbl_bcdust, self%tbl_kap, &
            self%sw_hygro_coreshell_ssa(:,iwav,:,:,:), &
            self%relh(icol,ilev), self%corefrac(icol,ilev), self%bcdust(icol,ilev), self%kappa(icol,ilev) ) ) &
            * pext(icol)

       pasm(icol) = table_interp( self%tbl_relh, self%tbl_corefrac, self%tbl_bcdust, self%tbl_kap, &
            self%sw_hygro_coreshell_asm(:,iwav,:,:,:), &
            self%relh(icol,ilev), self%corefrac(icol,ilev), self%bcdust(icol,ilev), self%kappa(icol,ilev) )

       pext(icol) = pext(icol)*self%totalmmr(icol,ilev)
       pabs(icol) = pabs(icol)*self%totalmmr(icol,ilev)
       pabs(icol) = max(0._r8,pabs(icol))
       pabs(icol) = min(pext(icol),pabs(icol))

       palb(icol) = 1._r8-pabs(icol)/max(pext(icol),1.e-40_r8)
       palb(icol) = 1._r8-pabs(icol)/max(pext(icol),1.e-40_r8)

    end do

  end subroutine sw_props

  !------------------------------------------------------------------------------
  ! returns long wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine lw_props(self, ncol, ilev, iwav, pabs)

    class(hygrocoreshell_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)

    integer :: icol

    do icol = 1, ncol

       pabs(icol) = table_interp( self%tbl_relh, self%tbl_corefrac,  self%tbl_bcdust, self%tbl_kap, &
            self%lw_hygro_coreshell_abs(:,iwav,:,:,:), &
            self%relh(icol,ilev), self%corefrac(icol,ilev), self%bcdust(icol,ilev), self%kappa(icol,ilev) )

       pabs(icol) = pabs(icol)*self%totalmmr(icol,ilev)
       pabs(icol) = max(0._r8,pabs(icol))
    end do

  end subroutine lw_props

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)

    type(hygrocoreshell_aerosol_optics), intent(inout) :: self

    deallocate(self%totalmmr)
    deallocate(self%corefrac)
    deallocate(self%bcdust)
    deallocate(self%kappa)
    deallocate(self%relh)

    nullify(self%tbl_corefrac)
    nullify(self%tbl_bcdust)
    nullify(self%tbl_kap)
    nullify(self%tbl_relh)
    nullify(self%sw_hygro_coreshell_ext)
    nullify(self%sw_hygro_coreshell_ssa)
    nullify(self%sw_hygro_coreshell_asm)
    nullify(self%lw_hygro_coreshell_abs)

  end subroutine destructor

end module hygrocoreshell_aerosol_optics_mod
