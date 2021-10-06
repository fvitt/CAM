module modal_spcam_aerosol_data_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use modal_cam_aerosol_data_mod, only: modal_cam_aerosol_data
  use aerosol_data_mod, only: ptr2d_t
  use constituents, only: pcnst, cnst_spec_class_gas, cnst_species_class, cnst_name
  use cam_history, only: outfld
  use ppgrid, only: pcols
  use physconst, only: gravit

  implicit none

  type, extends(modal_cam_aerosol_data) :: modal_spcam_aerosol_data
  contains
     procedure :: initialize => spcam_data_initialize
     procedure :: update => spcam_data_update
     procedure :: write_diags => spcam_write_diags
  end type modal_spcam_aerosol_data


contains

  subroutine spcam_data_initialize(self)
    class(modal_spcam_aerosol_data), intent(inout) :: self

    integer :: m

    call self%modal_cam_aerosol_data%initialize()

    do m=1, pcnst
       if (cnst_species_class(m) == cnst_spec_class_gas) then
          self%lq(m) = .true.
       end if
    end do

  end subroutine spcam_data_initialize



  subroutine spcam_data_update(self, raer, qqcw, raercol, raercol_cw, rgascol, colnum, dtinv)

    class(modal_spcam_aerosol_data), intent(inout) :: self
    type(ptr2d_t), intent(in) :: raer(:)
    type(ptr2d_t), intent(inout) :: qqcw(:)
    real(r8), intent(in) :: raercol(:,:)
    real(r8), intent(in) :: raercol_cw(:,:)
    real(r8), intent(in) :: rgascol(:,:)
    integer,  intent(in) :: colnum
    real(r8), intent(in) :: dtinv

    integer :: m,mm

    call self%cam_aerosol_data%update(raer, qqcw, raercol, raercol_cw, rgascol, colnum, dtinv)

    !
    ! Gas tendency
    !
    mm=0
    do m=1, pcnst
       if (cnst_species_class(m) == cnst_spec_class_gas) then
          mm=mm+1
          self%ptend%q(colnum, :, m) = (rgascol(:,mm)-self%state%q(colnum, :, m)) * dtinv
       end if
    end do

  end subroutine spcam_data_update

  subroutine spcam_write_diags(self)
    class(modal_spcam_aerosol_data), intent(in) :: self

    integer :: i,m
    real(r8) :: coltendgas(pcols)
    character(len=32) :: fieldnamegas

    !
    ! output column-integrated Gas tendency (this should be zero)
    !
    do m=1, pcnst
       if (cnst_species_class(m) == cnst_spec_class_gas) then
          do i=1, self%state%ncol
             coltendgas(i) = sum( self%state%pdel(i,:)*self%ptend%q(i,:,m) )/gravit
          end do
          fieldnamegas = trim(cnst_name(m)) // '_mixnuc1sp'
          call outfld( trim(fieldnamegas), coltendgas, pcols, self%state%lchnk)
       end if
    end do

  end subroutine spcam_write_diags

end module modal_spcam_aerosol_data_mod
