! Base Class for aerosol models

module aerosol_model_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physics_buffer, only: physics_buffer_desc
  use physics_types,  only: physics_state
  use cam_abortutils, only: endrun

  implicit none

  type, abstract :: aerosol_model
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
   contains
     procedure :: create => aero_create
     procedure(aero_model_init), deferred :: model_init
     procedure(aero_model_final), deferred :: model_final
     procedure :: destroy => aero_destroy
     procedure(aero_loadaer), deferred :: loadaer
  end type aerosol_model

  interface

     subroutine aero_model_init( self )
       import
       class(aerosol_model), intent(inout) :: self
     end subroutine aero_model_init

     subroutine aero_model_final( self )
       import
       class(aerosol_model), intent(inout) :: self
     end subroutine aero_model_final

     subroutine aero_loadaer( self, istart, istop, k, m, cs, phase, naerosol, vaerosol, hygro)
       import

       class(aerosol_model), intent(in) :: self
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

  end interface

contains

  subroutine aero_create( self, state, pbuf )
    class(aerosol_model), intent(inout) :: self
    type(physics_state), target, intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    self%state => state
    self%pbuf => pbuf

    call self%model_init()

  end subroutine aero_create

  subroutine aero_destroy( self )
    class(aerosol_model), intent(inout) :: self

    nullify(self%state)
    nullify(self%pbuf)

    call self%model_final()

  end subroutine aero_destroy


 end module aerosol_model_mod
