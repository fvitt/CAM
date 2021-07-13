! Base Class for aerosol models

module aerosol_model_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use physics_buffer, only: physics_buffer_desc
  use physics_types,  only: physics_state
  use cam_abortutils, only: endrun
  use ppgrid,         only: pcols, pver
  use ref_pres,       only: top_lev => trop_cloud_top_lev
  use physconst,      only: rhoh2o, mwh2o, r_universal, rh2o, pi

  implicit none

  type, abstract :: aerosol_model
     integer :: mtotal
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
     real(r8), allocatable :: amcubecoef(:)
     real(r8), allocatable :: argfactor(:)
   contains
     procedure :: create => aero_create
     procedure(aero_model_init), deferred :: model_init
     procedure(aero_model_final), deferred :: model_final
     procedure :: destroy => aero_destroy
     procedure(aero_loadaer), deferred :: loadaer
     procedure :: ccncalc => aero_ccncalc
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

  integer,  parameter :: psat=6    ! number of supersaturations to calc ccn concentration
  real(r8), parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn concentration
                       (/ 0.02_r8, 0.05_r8, 0.1_r8, 0.2_r8, 0.5_r8, 1.0_r8 /)

  !     mathematical constants
  real(r8), parameter :: zero     = 0._r8
  real(r8), parameter :: third    = 1._r8/3._r8
  real(r8), parameter :: twothird = 2._r8*third
  real(r8), parameter :: sixth    = 1._r8/6._r8
  real(r8), parameter :: sq2      = sqrt(2._r8)
  real(r8), parameter :: sqpi     = sqrt(pi)

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

  subroutine aero_ccncalc(self, state, cs, ccn)

    ! calculates number concentration of aerosols activated as CCN at
    ! supersaturation supersat.
    ! assumes an internal mixture of a multiple externally-mixed aerosol modes
    ! cgs units

    ! Ghan et al., Atmos. Res., 1993, 198-221.

    ! arguments
    class(aerosol_model), intent(in) :: self

    type(physics_state), target, intent(in)    :: state


    real(r8), intent(in)  :: cs(pcols,pver)       ! air density (kg/m3)
    real(r8), intent(out) :: ccn(pcols,pver,psat) ! number conc of aerosols activated at supersat (#/m3)

    ! local

    integer :: lchnk ! chunk index
    integer :: ncol  ! number of columns
    real(r8), pointer :: tair(:,:)     ! air temperature (K)

    real(r8) :: naerosol(pcols) ! interstit+activated aerosol number conc (/m3)
    real(r8) :: vaerosol(pcols) ! interstit+activated aerosol volume conc (m3/m3)

    real(r8) :: amcube(pcols)
    real(r8) :: surften_coef
    real(r8) :: a(pcols) ! surface tension parameter
    real(r8) :: hygro(pcols)  ! aerosol hygroscopicity
    real(r8) :: sm(pcols)  ! critical supersaturation at mode radius
    real(r8) :: arg(pcols)

    integer :: l,m,n,i,k
    real(r8) :: log,cc
    real(r8) :: smcoef(pcols)
    integer :: phase ! phase of aerosol

    real(r8), parameter :: smcoefcoef=2._r8/sqrt(27._r8)
    real(r8), parameter :: super(psat)=supersat(:psat)*0.01_r8
    real(r8), parameter :: surften=0.076_r8    ! surface tension of water w/respect to air (N/m)

    surften_coef=2._r8*mwh2o*surften/(r_universal*rhoh2o)

    !-------------------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    tair  => state%t


    ccn = 0._r8
    do k=top_lev,pver

       do i=1,ncol
          a(i)=surften_coef/tair(i,k)
          smcoef(i)=smcoefcoef*a(i)*sqrt(a(i))
       end do

       do m=1,self%mtotal

          phase=3 ! interstitial+cloudborne

          call self%loadaer( 1, ncol, k, m, cs, phase, naerosol, vaerosol, hygro)

          where(naerosol(:ncol)>1.e-3_r8 .and. hygro(:ncol)>0._r8)
             amcube(:ncol)=self%amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
             sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol)) ! critical supersaturation
          elsewhere
             sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
          endwhere
          do l=1,psat
             do i=1,ncol
                arg(i)=self%argfactor(m)*log(sm(i)/super(l))
                ccn(i,k,l)=ccn(i,k,l)+naerosol(i)*0.5_r8*(1._r8-erf(arg(i)))
             enddo
          enddo
       enddo
    enddo
    ccn(:ncol,:,:)=ccn(:ncol,:,:)*1.e-6_r8 ! convert from #/m3 to #/cm3

  end subroutine aero_ccncalc


 end module aerosol_model_mod
