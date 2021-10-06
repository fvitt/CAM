module aerosol_data_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  ! ptr2d_t is used to create arrays of pointers to 2D fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  type, abstract :: aerosol_data
     real(r8), allocatable :: sigmag_amode(:)! geometric standard deviation for each aerosol mode
     real(r8), allocatable :: dgnumlo_amode(:)
     real(r8), allocatable :: dgnumhi_amode(:)
     character(len=24), allocatable :: fieldname(:)    ! names for drop nuc tendency output fields
     character(len=24), allocatable :: fieldname_cw(:) ! names for drop nuc tendency output fields
     integer, allocatable :: indexer(:,:)
     integer, allocatable :: cnstndx(:,:)
     integer :: ncnst_tot
     integer :: mtotal
     integer, allocatable :: nmasses(:)
     integer, allocatable :: nspec(:)
     real(r8), allocatable :: coltend(:,:)       ! column tendency for diagnostic output
     real(r8), allocatable :: coltend_cw(:,:)    ! column tendency
   contains
     procedure(aero_initialize), deferred :: initialize
     procedure(aero_loadaer),    deferred :: loadaer
     procedure(aero_set_ptrs),   deferred :: set_ptrs
     procedure(aero_update),     deferred :: update
  end type aerosol_data

  interface

     subroutine aero_initialize(self)
       import
       class(aerosol_data), intent(inout) :: self
     end subroutine aero_initialize

     subroutine aero_loadaer( self, istart, istop, k, m, cs, phase, naerosol, vaerosol, hygro )
       import

       class(aerosol_data), intent(in) :: self
       ! inputs
       integer,  intent(in) :: istart      ! start column index (1 <= istart <= istop <= pcols)
       integer,  intent(in) :: istop       ! stop column index
       integer,  intent(in) :: m           ! mode or bin index
       integer,  intent(in) :: k           ! level index
       real(r8), intent(in) :: cs(:,:)     ! air density (kg/m3)
       integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

       ! outputs
       real(r8), intent(out) :: naerosol(:)  ! number conc (1/m3)
       real(r8), intent(out) :: vaerosol(:)  ! volume conc (m3/m3)
       real(r8), intent(out) :: hygro(:)     ! bulk hygroscopicity of mode

     end subroutine aero_loadaer

     subroutine aero_set_ptrs( self, raer, qqcw )
       import

       class(aerosol_data), intent(in) :: self
       type(ptr2d_t), intent(out) :: raer(:)
       type(ptr2d_t), intent(out) :: qqcw(:)

     end subroutine aero_set_ptrs

     subroutine aero_update(self, raer, qqcw, raercol, raercol_cw, rgascol, colnum, dtinv)
       import
       class(aerosol_data), intent(inout) :: self
       type(ptr2d_t), intent(in) :: raer(:)
       type(ptr2d_t), intent(inout) :: qqcw(:)
       real(r8), intent(in) :: raercol(:,:)
       real(r8), intent(in) :: raercol_cw(:,:)
       real(r8), intent(in) :: rgascol(:,:)
       integer,  intent(in) :: colnum
       real(r8), intent(in) :: dtinv

     end subroutine aero_update

  end interface

contains


end module aerosol_data_mod
