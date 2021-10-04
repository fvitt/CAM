module cam_aerosol_data_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physics_buffer, only: physics_buffer_desc
  use physics_types,  only: physics_state
  use physics_types,  only: physics_ptend
  use constituents,   only: pcnst
  use ppgrid,         only: pver, pcols
  use ref_pres,       only: top_lev => trop_cloud_top_lev
  use physconst,      only: gravit

  use aerosol_data_mod, only: aerosol_data, ptr2d_t

  implicit none

  type, abstract, extends(aerosol_data) :: cam_aerosol_data
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
     type(physics_ptend), pointer :: ptend => null()
     logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
   contains
     procedure :: update => cam_data_update
  end type cam_aerosol_data

contains

  subroutine cam_data_update(self, pdel, raer, qqcw, raercol, raercol_cw, colnum, dtinv, &
                             coltend, coltend_cw)
    class(cam_aerosol_data), intent(inout) :: self
    real(r8), intent(in) :: pdel(:,:)    ! pressure thickess of layer (Pa)
    type(ptr2d_t), intent(in) :: raer(:)
    type(ptr2d_t), intent(inout) :: qqcw(:)
    real(r8), intent(in) :: raercol(:,:)
    real(r8), intent(in) :: raercol_cw(:,:)
    integer,  intent(in) :: colnum
    real(r8), intent(in) :: dtinv

    real(r8), intent(out) :: coltend(:)
    real(r8), intent(out) :: coltend_cw(:)

    integer :: m,l, lptr, mm
    real(r8) :: raertend(pver)
    real(r8) :: qqcwtend(pver)

    coltend = 0._r8
    raertend = 0._r8
    qqcwtend = 0._r8
    coltend_cw = 0._r8

    do m = 1, self%mtotal
       do l = 0, self%nmasses(m)

          mm   = self%indexer(m,l)
          lptr = self%cnstndx(m,l)

          raertend(top_lev:pver) = (raercol(top_lev:pver,mm) - raer(mm)%fld(colnum,top_lev:pver))*dtinv
          qqcwtend(top_lev:pver) = (raercol_cw(top_lev:pver,mm) - qqcw(mm)%fld(colnum,top_lev:pver))*dtinv

          coltend(mm)    = sum( pdel(colnum,:)*raertend )/gravit
          coltend_cw(mm) = sum( pdel(colnum,:)*qqcwtend )/gravit

          ! some interstetial species are might not be advected, check:
          if (lptr>0) then ! adveced species
             self%ptend%q(colnum,:,lptr) = 0.0_r8
             self%ptend%q(colnum,top_lev:pver,lptr) = raertend(top_lev:pver)        ! set tendencies for interstitial aerosol
          else
             raer(mm)%fld(colnum,:) = 0.0_r8
             raer(mm)%fld(colnum,top_lev:pver) = raercol(top_lev:pver,mm) ! update interstitial aerosol
          end if

          qqcw(mm)%fld(colnum,:) = 0.0_r8
          qqcw(mm)%fld(colnum,top_lev:pver) = raercol_cw(top_lev:pver,mm) ! update cloud-borne aerosol

       end do
    end do

  end subroutine cam_data_update

end module cam_aerosol_data_mod
