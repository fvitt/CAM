module majorsp_diffusion

!--------------------------------------------------------------------------
! This module computes the diffusion of major species (O2 and O) mass mixing
! ratio. This routine computes both the molecular and eddy diffusivity. This
! is adapted from the major species diffusion calculation of TIME-GCM.
!
! Calling sequence:
!   initialization:
!      init
!         call mspd_init
!
!   interfacing:
!      tphysac
!         (after vertical_diffusion_tend)
!         call mspd_intr
!            call mspdiff
!
!---------------------------Code history--------------------------------
! Adapted from TIME-GCM (comp.F): H.-L. Liu, Nov 2003
!--------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver, pverp
  use constituents, only: pcnst, cnst_name, cnst_get_ind, cnst_mw
  use cam_history,  only: outfld
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc
  use infnan,       only: nan, assignment(=)

  implicit none

  private          ! Make default type private to the module
  save
!-----------------------
! Public interfaces
!-----------------------
  public mspd_init   ! Initialization
  public mspd_intr   ! Full routine
!-----------------------
! Private data
!-----------------------

!  real(r8) :: rmass_o2, rmass_o1, rmass_n2               ! molecular weight kg/kmol
  real(r8) :: rmass_o2, rmass_o1, rmass_h, rmass_he, rmass_n2     ! molecular weight kg/kmol
  real(r8) :: rmassinv_o2, rmassinv_o1, rmassinv_n2      ! 1/rmass_o2...
  real(r8) :: phi(2,3)                                   ! mutual diffusion constants of
                                                         ! major constituents
  real(r8) :: delta(2,2)                                 ! unit matrix

  real(r8), parameter :: t00=273._r8                     ! reference temperature
  real(r8), parameter :: ptref=5.e-5_r8                  ! thermosphere reference pressure (Pa)
  real(r8), parameter :: tau=1.86e3_r8                   ! diffusive time constant (sec).
  real(r8), parameter :: protonmass=1.6726e-27_r8        ! Proton mass (kg)
  real(r8), parameter :: mmrMin=1.e-20_r8                ! lower limit of o2 and o mixing ratio
  real(r8), parameter :: N2mmrMin=1.e-6_r8               ! lower limit of n2 mixing ratios
  real(r8), parameter :: HEmmrMin=1.e-7_r8               ! lower limit of he mixing ratios
  real(r8), parameter :: HEmmrMax=0.9_r8                 ! upper limit of he mixing ratios

  integer :: indx_O2                                     ! cnst index for o2
  integer :: indx_O                                      ! cnst index for o
  integer :: indx_H                                      ! cnst index for h
  integer :: indx_HE                                     ! cnst index for he
!  integer, parameter :: io2=1, io1=2                     ! local indices to o2 , o respectively
!  logical :: fixed_ubc(2)                                ! flag for fixed upper boundary condition
  integer, parameter :: io2=1, io1=2, ihe=3              ! local indices to o2 , o, and he respectively
  logical :: fixed_ubc(3)                                ! flag for fixed upper boundary condition

  real(r8) :: o2mmr_ubc(pcols)                           ! MMR of O2 at top boundary (specified)
  real(r8) :: ommr_ubc(pcols)                            ! MMR of O at top boundary
  real(r8) :: hemmr_ubc(pcols)                           ! MMR flux of HE at top boundary

!  character(len=8), private :: mjdiffnam(2)              ! names of v-diff tendencies
  character(len=10), private :: mjdiffnam(5)              ! names of v-diff tendencies

  logical, parameter :: debug = .false.

contains

!===============================================================================
  subroutine mspd_init()

    !-------------------------------------------------------------------------------
    ! Define constants and coeficient matrices, phi and delta, in the initialization.
    !-------------------------------------------------------------------------------
    use constituents, only: cnst_mw, cnst_fixed_ubc
    use cam_history,  only: addfld, add_default
    use phys_control, only: phys_getopts

    !------------------------------Arguments--------------------------------

    !---------------------------Local storage-------------------------------
    logical :: history_waccmx

    call phys_getopts(history_waccmx_out=history_waccmx)

    !-----------------------------------------------------------
    ! Get required molecular weights
    !-----------------------------------------------------------
    call cnst_get_ind('O2', indx_O2, abort=.true.)
    call cnst_get_ind('O',  indx_O, abort=.true.)
    call cnst_get_ind('H',  indx_H, abort=.true.)
    call cnst_get_ind('HE', indx_HE, abort=.true.)

    rmass_o2 = cnst_mw(indx_O2)
    rmass_o1 = cnst_mw(indx_O)
    rmass_h  = cnst_mw(indx_H)
    rmass_he = cnst_mw(indx_HE)
    rmass_n2 = 28._r8

    rmassinv_o2 = 1._r8/rmass_o2
    rmassinv_o1 = 1._r8/rmass_o1
    rmassinv_n2 = 1._r8/rmass_n2

    !--------------------------------------------------------------------
    ! Get fixed upper boundary flags and set vertical range for diffusion
    !--------------------------------------------------------------------
    fixed_ubc(io2) = cnst_fixed_ubc(indx_O2)
    fixed_ubc(io1) = cnst_fixed_ubc(indx_O)
    fixed_ubc(ihe) = cnst_fixed_ubc(indx_HE)

    !------------------------------------------------
    ! Set diffusion constants and setup matrix
    !------------------------------------------------
    phi(:,1)=(/0._r8  ,0.673_r8/)
    phi(:,2)=(/1.35_r8,0._r8   /)
    phi(:,3)=(/1.11_r8,0.769_r8/)
    delta(:,1)=(/1._r8,0._r8/)
    delta(:,2)=(/0._r8,1._r8/)

   ! Set names of major diffusion tendencies and declare them as history variables
    mjdiffnam(1) = 'MD'//cnst_name(indx_O2)
    call addfld (mjdiffnam(1),(/ 'lev' /), 'A','kg/kg/s','Major diffusion of '//cnst_name(indx_O2))
    mjdiffnam(2) = 'MD'//cnst_name(indx_O)
    call addfld (mjdiffnam(2),(/ 'lev' /), 'A','kg/kg/s','Major diffusion of '//cnst_name(indx_O))

   ! Set names of major diffusion tendencies from comp_wx TGCM routine and declare them as history variables
    mjdiffnam(3) = 'comp_wx_O2'
    call addfld (mjdiffnam(3),(/ 'lev' /), 'A','kg/kg','comp_wx major diffusion of '//cnst_name(indx_O2))
    mjdiffnam(4) = 'comp_wx_O'
    call addfld (mjdiffnam(4),(/ 'lev' /), 'A','kg/kg','comp_wx major diffusion of '//cnst_name(indx_O))
    mjdiffnam(5) = 'comp_wx_HE'
    call addfld (mjdiffnam(5),(/ 'lev' /), 'A','kg/kg','comp_wx major diffusion of '//cnst_name(indx_HE))

    call addfld ('MBARV' , (/ 'lev' /),'I','g/mole','Variable Mean Mass')

    if (history_waccmx) then
       call add_default (mjdiffnam(1), 1, ' ')
       call add_default (mjdiffnam(2), 1, ' ')
       call add_default (mjdiffnam(3), 1, ' ')
       call add_default (mjdiffnam(4), 1, ' ')
       call add_default (mjdiffnam(5), 1, ' ')
       call add_default ('MBARV', 1, ' ')
    end if

  end subroutine mspd_init

!===============================================================================
  subroutine mspd_intr(ztodt    ,state    ,ptend)

!-------------------------------------------------------------------------------
! interface routine. output tendency.
!-------------------------------------------------------------------------------
    use physics_types,   only: physics_state, physics_ptend
    use upper_bc,        only: ubc_get_vals
    use upper_bc,        only: ubc_get_flxs
    use air_composition, only: rairv, mbarv
    use ref_pres,        only: nbot_molec
    use physconst,    only: gravit
    use helium_ubc_mod,  only: helium_ubc_fluxes

!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    type(physics_state), intent(in)     :: state   ! Physics state variables
    type(physics_ptend), intent(inout)  :: ptend   ! indivdual parameterization tendencies
!---------------------------Local storage-------------------------------
    real(r8) :: rztodt                             ! 1/ztodt
!    real(r8) :: tendo2o(pcols,pver,2)              ! temporary array for o2 and o tendency
    real(r8) :: tendo2ohe(pcols,pver,3)            ! temporary array for o2 o, and he tendencies
    real(r8) :: ubc_mmr(pcols,pcnst)               ! upper bndy mixing ratios (kg/kg)
    real(r8) :: ubc_t(pcols)                       ! upper bndy temperature (K)
    real(r8) :: ubc_flux(pcols,pcnst)              ! upper bndy mixing ratio flux (kg/kg/s?)
    integer :: lchnk                               ! chunk identifier
    integer :: ncol                                ! number of atmospheric columns
!    integer :: i, k                                ! indexing integers
    integer :: i, k, kk, icol                      ! indexing integers

    ! For comp_wx call
    integer :: nlevp1

    real(r8) :: tlbc,bo2,bo1,bhe,bh,he_ubc,p_ubc     ! For lower boundary
    real(r8) :: step,dfactor,pscaleheight,expzmid,p0
    real(r8),dimension(nbot_molec) :: &
      difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
      o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm,dz,expzm
    real(r8),dimension(3,nbot_molec) :: prod
    real(r8),dimension(3,3,nbot_molec) :: loss
    real(r8),dimension(nbot_molec) :: o2_upd,o1_upd,he_upd
    real(r8),dimension(pcols,pver) :: o2_upd_cols,o1_upd_cols,he_upd_cols,h_upd_cols
    real(r8),dimension(pcols,pver) :: o2_upd_cols_tend,o1_upd_cols_tend,he_upd_cols_tend

    !--------------------------------------------------------------------------------------------
    ! local constants
    !--------------------------------------------------------------------------------------------
    rztodt = 1._r8/ztodt
    lchnk = state%lchnk
    ncol  = state%ncol

!    !----------------------------------------------------------------------------------------------
!    ! Store the o2 and o tendencies calculated from vertical_diffusion (due to eddy diffusion only)
!    !----------------------------------------------------------------------------------------------
!    tendo2o(:ncol,:,io2) = ptend%q(:ncol,:,indx_O2)
!    tendo2o(:ncol,:,io1) = ptend%q(:ncol,:,indx_O)

    !----------------------------------------------------------------------------------------------
    ! Store the o2, o, and he tendencies calculated from vertical_diffusion (due to eddy diffusion only)
    !----------------------------------------------------------------------------------------------
    tendo2ohe(:ncol,:,io2) = ptend%q(:ncol,:,indx_O2)
    tendo2ohe(:ncol,:,io1) = ptend%q(:ncol,:,indx_O)
    tendo2ohe(:ncol,:,ihe) = ptend%q(:ncol,:,indx_HE)

    o2_upd_cols_tend(:ncol,:) = ptend%q(:ncol,:,indx_O2)
    o1_upd_cols_tend(:ncol,:) = ptend%q(:ncol,:,indx_O)
    he_upd_cols_tend(:ncol,:) = ptend%q(:ncol,:,indx_HE)

    !----------------------------------------------------------------------
    ! Operate on copies of the input states, convert to tendencies at end.
    !----------------------------------------------------------------------
    ptend%q(:ncol,:,indx_O2) = state%q(:ncol,:,indx_O2)
    ptend%q(:ncol,:,indx_O) = state%q(:ncol,:,indx_O)
    ptend%q(:ncol,:,indx_HE) = state%q(:ncol,:,indx_HE)

    o2_upd_cols(:ncol,:) = state%q(:ncol,:,indx_O2)
    o1_upd_cols(:ncol,:) = state%q(:ncol,:,indx_O)
    he_upd_cols(:ncol,:) = state%q(:ncol,:,indx_HE)
    h_upd_cols(:ncol,:)  = state%q(:ncol,:,indx_H)

    if (fixed_ubc(io2) .or. fixed_ubc(io1)) then
       !-------------------------------------------
       ! set upper boundary values of O2 and O MMR.
       !-------------------------------------------
       call ubc_get_vals( lchnk, ncol, state%pint, state%zi, ubc_t, ubc_mmr )
       o2mmr_ubc(:ncol) = ubc_mmr(:ncol,indx_O2)
       ommr_ubc(:ncol) = ubc_mmr(:ncol,indx_O)
    endif

       call ubc_get_flxs( state%lchnk, ncol, state%pint, state%zi, state%t, state%q, state%omega, state%phis, ubc_flux )
!       hemmr_ubc(:ncol) = ubc_flux(:ncol,indx_HE)
       hemmr_ubc(:ncol) = helium_ubc_fluxes(:ncol,lchnk)

if (masterproc .and. debug) write(iulog,*) 'comp_wx: lchnk,hemmr_ubc(:ncol) all columns after assignment: ', lchnk,hemmr_ubc(:ncol)

    ! Since this is a combined tendency, retain the old name for output
    ! and debugging purposes.
    ptend%name  = trim(ptend%name)//"+mspd"
    ptend%lq(indx_O2) = .TRUE.
    ptend%lq(indx_O) = .TRUE.
    !---------------------------------------------
    ! Call the major species diffusion subroutine.
    !---------------------------------------------
    call mspdiff (lchnk      ,ncol       ,                                     &
                  state%t    ,ptend%q    ,state%pmid ,state%pint ,             &
                  state%pdel ,ztodt      ,rairv(:,:,lchnk),  mbarv(:,:,lchnk))

    !----------------------------------------------------------------------
    ! Operate on copies of the input states for comp_wx also, convert to tendencies at end.
    !----------------------------------------------------------------------
    ptend%q(:ncol,:,indx_O2) = state%q(:ncol,:,indx_O2)
    ptend%q(:ncol,:,indx_O) = state%q(:ncol,:,indx_O)
    ptend%q(:ncol,:,indx_HE) = state%q(:ncol,:,indx_HE)

    step = ztodt/2._r8

    !
    ! Eddy diffusion set to zero since already calculated in vertical_diffusion
    !
    dfactor = 0._r8
    difk(:) = 0._r8
    wmid(:) = 0._r8
    expzmid = 1._r8
    !
    ! Chemical production/loss set to zero since already done is chemistry
    !
    prod(1:3,1:nbot_molec)   = 0._r8
    loss(1:3,1:3,1:nbot_molec) = 0._r8
    !
    ! Set adv and nm to values to zero
    !
    o2_hadv(1:nbot_molec) = 0._r8
    o1_hadv(1:nbot_molec) = 0._r8
    he_hadv(1:nbot_molec) = 0._r8
    o2_nm(1:nbot_molec)   = 0._r8
    o1_nm(1:nbot_molec)   = 0._r8
    he_nm(1:nbot_molec)   = 0._r8

    tn = nan
    tni = nan

    do iCol = 1,ncol

      tlbc   = state%t(iCol,nbot_molec+1)
      bo2    = state%q(iCol,nbot_molec+1,indx_O2)
      bo1    = state%q(iCol,nbot_molec+1,indx_O)
      bhe    = state%q(iCol,nbot_molec+1,indx_HE)
      bh     = state%q(iCol,nbot_molec+1,indx_H)
      he_ubc = hemmr_ubc(iCol)

      if (masterproc .and. iCol == 1 .and. debug) &
           write(iulog,*) 'comp_wx: iCol,he_ubc,hemmr_ubc(iCol) before ubc calc first column: ', iCol, he_ubc, hemmr_ubc(iCol)

      kk = 0
      do k = nbot_molec,2,-1

        kk = kk + 1
        tn(kk)       = state%t(iCol,k)
        tni(kk)      = .5_r8 * (state%t(iCol,k) + state%t(iCol,k-1))
        o2i(kk)      = .5_r8 * (state%q(iCol,k,indx_O2) + state%q(iCol,k-1,indx_O2))
        o1i(kk)      = .5_r8 * (state%q(iCol,k,indx_O) + state%q(iCol,k-1,indx_O))
        hei(kk)      = .5_r8 * (state%q(iCol,k,indx_HE) + state%q(iCol,k-1,indx_HE))
        mbar(kk)     = mbarv(iCol,k,lchnk)
        barm(kk)     = .5_r8 * (mbarv(iCol,k,lchnk) + mbarv(iCol,k-1,lchnk))
        pScaleHeight = .5_r8*(rairv(iCol,k,lchnk)*state%t(iCol,k) + rairv(iCol,k-1,lchnk)*state%t(iCol,k)) / gravit
!        wmid(kk)     = -state%omega(iCol,k) / (0.5_r8 * (state%pint(iCol,k-1) + state%pint(iCol,k))) * pScaleHeight
        dz(kk)       = (state%pmid(iCol,k) - state%pmid(iCol,k-1)) / state%pint(iCol,k)
        expzm(kk)    = state%pmid(iCol,k) / ptref

      enddo ! kk=1,nbot_molec-1
      !
      ! Top:
      !
      tn(nbot_molec)	 = state%t(iCol,1)
      tni(nbot_molec)	 = 1.5_r8*state%t(iCol,1)-.5_r8*state%t(iCol,2)
      o2i(nbot_molec)	 = 1.5_r8*state%q(iCol,1,indx_O2)-.5_r8*state%q(iCol,2,indx_O2)
      o1i(nbot_molec)	 = 1.5_r8*state%q(iCol,1,indx_O)-.5_r8*state%q(iCol,2,indx_O)
      hei(nbot_molec)	 = 1.5_r8*state%q(iCol,1,indx_HE)-.5_r8*state%q(iCol,2,indx_HE)
      mbar(nbot_molec)   = mbarv(iCol,1,lchnk)
      barm(nbot_molec)   = 1.5_r8*mbarv(iCol,1,lchnk)-.5_r8*mbarv(iCol,2,lchnk)
      pScaleHeight	 = .5_r8*(rairv(iCol,1,lchnk)*state%t(iCol,1) + rairv(iCol,2,lchnk)*state%t(iCol,2)) / gravit
!      wmid(nbot_molec) = -state%omega(iCol,1) / (0.5_r8 * (state%pint(iCol,1) + state%pint(iCol,2))) * pScaleHeight
      p_ubc = state%pmid(iCol,1)*state%pmid(iCol,1)/state%pmid(iCol,2)
      dz(nbot_molec)    = (state%pmid(iCol,1)-p_ubc)/state%pint(iCol,1)
      expzm(nbot_molec) = state%pmid(iCol,1) / ptref

!      he_ubc = 3.0E-12_r8

if (masterproc.and.debug) write(iulog,*) 'mspd_intr: iCol, lchnk, he_ubc before comp_wx: ', iCol, lchnk, he_ubc

      call comp_wx(iCol,lchnk,step,dfactor,tlbc,bo2,bo1,bh,bhe,he_ubc,difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
	       o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm,prod,loss, &
	       nbot_molec,dz,expzm,expzmid,o2_upd,o1_upd,he_upd)

if (debug) then
if (masterproc .and. iCol == 1) write(iulog,*) 'mspd_intr: after comp_wx o2_upd first column all levels ',iCol, o2_upd(:)
if (masterproc .and. iCol == 1) write(iulog,*) 'mspd_intr: after comp_wx o1_upd first column all levels ',iCol, o1_upd(:)
if (masterproc .and. iCol == 1) write(iulog,*) 'mspd_intr: after comp_wx he_upd first column all levels ',iCol, he_upd(:)

!if (masterproc .and. iCol <= 10) write(iulog,*) 'mspd_intr: after comp_wx o2_upd first column all levels ',iCol, o2_upd(:)
!if (masterproc .and. iCol <= 10) write(iulog,*) 'mspd_intr: after comp_wx o1_upd first column all levels ',iCol, o1_upd(:)
!if (masterproc .and. iCol <= 10) write(iulog,*) 'mspd_intr: after comp_wx he_upd first column all levels ',iCol, he_upd(:)
end if

       kk = 0
       do k = 1,nbot_molec

	 kk = nbot_molec - k + 1

         o2_upd_cols(iCol,kk) = o2_upd(k)
         o1_upd_cols(iCol,kk) = o1_upd(k)
         he_upd_cols(iCol,kk) = he_upd(k)

       enddo

     enddo ! iCol loop

!write(iulog,*) 'mspd_intr: MIN/MAX o2_upd_cols,o1_upd_cols,he_upd_cols before output call : ', &
!                     MINVAL(o2_upd_cols(:,:)),MAXVAL(o2_upd_cols(:,:)),  &
!                     MINVAL(o1_upd_cols(:,:)),MAXVAL(o1_upd_cols(:,:)),  &
!                     MINVAL(he_upd_cols(:,:)),MAXVAL(he_upd_cols(:,:))
!
!write(iulog,*) 'mspd_intr: MIN/MAX ptend%q(:,:,indx_O2),ptend%q(:,:,indx_O),ptend(:,:,indx_HE) before output call : ', &
!                     MINVAL(ptend%q(:,:,indx_O2)),MAXVAL(ptend%q(:,:,indx_O2)),  &
!                     MINVAL(ptend%q(:,:,indx_O)),MAXVAL(ptend%q(:,:,indx_O)),  &
!                     MINVAL(ptend%q(:,:,indx_HE)),MAXVAL(ptend%q(:,:,indx_HE))
if (debug) then
if (masterproc) write(iulog,*) 'mspd_intr: after iCol loop o2_upd_cols first column all levels ',o2_upd_cols(1,:)
if (masterproc) write(iulog,*) 'mspd_intr: after iCol loop o1_upd_cols first column all levels ',o1_upd_cols(1,:)
if (masterproc) write(iulog,*) 'mspd_intr: after iCol loop he_upd_cols first column all levels ',he_upd_cols(1,:)
end if
    !---------------------------------------------------------------
    ! Check for N2 greater than one
    !---------------------------------------------------------------
    do i=1,ncol
       do k=1,nbot_molec

	  if(1._r8-mmrMin-o2_upd_cols(i,k)-o1_upd_cols(i,k)-he_upd_cols(i,k)-h_upd_cols(i,k) < 0._r8) then
	     o2_upd_cols(i,k) = o2_upd_cols(i,k)*((1._r8-N2mmrMin-h_upd_cols(i,k))/(o2_upd_cols(i,k)+o1_upd_cols(i,k)+he_upd_cols(i,k)))
	     o1_upd_cols(i,k) = o1_upd_cols(i,k)*((1._r8-N2mmrMin-h_upd_cols(i,k))/(o2_upd_cols(i,k)+o1_upd_cols(i,k)+he_upd_cols(i,k)))
	     he_upd_cols(i,k) = he_upd_cols(i,k)*((1._r8-N2mmrMin-h_upd_cols(i,k))/(o2_upd_cols(i,k)+o1_upd_cols(i,k)+he_upd_cols(i,k)))
	  endif

       enddo
    enddo
if (debug) then
if (masterproc) write(iulog,*) 'mspd_intr: after N2 check loop o2_upd_cols first column all levels ',o2_upd_cols(1,:)
if (masterproc) write(iulog,*) 'mspd_intr: after N2 check loop o1_upd_cols first column all levels ',o1_upd_cols(1,:)
if (masterproc) write(iulog,*) 'mspd_intr: after N2 check loop he_upd_cols first column all levels ',he_upd_cols(1,:)
end if
    call outfld(mjdiffnam(3),o2_upd_cols(:,:),pcols,lchnk)
    call outfld(mjdiffnam(4),o1_upd_cols(:,:),pcols,lchnk)
    call outfld(mjdiffnam(5),he_upd_cols(:,:),pcols,lchnk)

    !---------------------------------------------
    ! Update O2 and O tendencies and output
    !---------------------------------------------
    do k=1,pver
       do i=1,ncol
!          ptend%q(i,k,indx_O2) = (ptend%q(i,k,indx_O2)-state%q(i,k,indx_O2))*rztodt  &
!                                 +tendo2ohe(i,k,io2)
!!                                 +tendo2o(i,k,io2)
!          ptend%q(i,k,indx_O) = (ptend%q(i,k,indx_O)-state%q(i,k,indx_O))*rztodt     &
!                                 +tendo2ohe(i,k,io1)
!!                                 +tendo2o(i,k,io1)

	  o2_upd_cols_tend(i,k) = (o2_upd_cols(i,k) - state%q(i,k,indx_O2)) * rztodt  &
                                  + tendo2ohe(i,k,io2)
	  o1_upd_cols_tend(i,k) = (o1_upd_cols(i,k) - state%q(i,k,indx_O)) * rztodt  &
                                  + tendo2ohe(i,k,io1)
	  he_upd_cols_tend(i,k) = (he_upd_cols(i,k) - state%q(i,k,indx_HE)) * rztodt  &
                                  + tendo2ohe(i,k,ihe)

	  ptend%q(i,k,indx_O2) = o2_upd_cols_tend(i,k)
	  ptend%q(i,k,indx_O)  = o1_upd_cols_tend(i,k)
	  ptend%q(i,k,indx_HE) = he_upd_cols_tend(i,k)

       enddo
    enddo
if (debug) then
if (masterproc) write(iulog,*) 'mspd_intr: after ptend calc o2_upd_cols_tend first column all levels ',o2_upd_cols_tend(1,:)
if (masterproc) write(iulog,*) 'mspd_intr: after ptend calc o1_upd_cols_tend first column all levels ',o1_upd_cols_tend(1,:)
if (masterproc) write(iulog,*) 'mspd_intr: after ptend calc he_upd_cols_tend first column all levels ',he_upd_cols_tend(1,:)
end if
    call outfld(mjdiffnam(1),ptend%q(1,1,indx_O2),pcols,lchnk)
    call outfld(mjdiffnam(2),ptend%q(1,1,indx_O),pcols,lchnk)

  end subroutine mspd_intr

!-----------------------------------------------------------------------
!  pure subroutine comp_wx(step,dfactor,tlbc,bo2,bo1,bh,bhe,he_ubc, &
  subroutine comp_wx(iCol,lchnk,step,dfactor,tlbc,bo2,bo1,bh,bhe,he_ubc, &
    difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
    o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm, &
    prod,loss,nlevp1,dz,expzm,expzmid,o2_upd,o1_upd,he_upd)

! advance major species O2, O, He and N2

!    use params_module,only:nlevp1,dz
!    use cons_module,only:expzm,expzmid,grav,p0,rmass_o2,rmass_o1,rmass_he,rmass_n2
!    use lbc_module,only:fb,b
!    use matutil_module,only:matinv3
!     use physconst,    only: gravit

    integer,intent(in) :: nlevp1, iCol, lchnk

    real(r8),intent(in) :: step,dfactor,tlbc,bo2,bo1,bh,bhe,he_ubc,expzmid
    real(r8),dimension(nlevp1),intent(in) :: &
      difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
      o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm,dz,expzm
    real(r8),dimension(3,nlevp1),intent(inout) :: prod
    real(r8),dimension(3,3,nlevp1),intent(inout) :: loss
    real(r8),dimension(nlevp1),intent(out) :: o2_upd,o1_upd,he_upd

    ! exponent factor for diff_fac
    real(r8),dimension(3),parameter :: ss = (/1.710_r8,1.749_r8,1.718_r8/)

    ! mutual thermal diffusion coefficients among major species
    real(r8),dimension(3,4),parameter :: &
      psi = reshape( &
       (/0.0_r8 ,0.673_r8,0.270_r8, &
        1.35_r8,0.0_r8  ,0.404_r8, &
        2.16_r8,1.616_r8,0.0_r8  , &
        1.11_r8,0.769_r8,0.322_r8/),(/3,4/))

    real(r8),parameter :: tau = 1.86e3_r8, t00 = 273, &
      thdiffalpha = -0.38_r8 ! thermal diffusion coefficient (alpha) for Helium
    integer,dimension(3,3),parameter :: delta = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
    integer :: k,m,n
    real(r8) :: &
      bn2,bmbar, &      ! at midpoint level 0 (not interface level 1)
      flx00,o1_ub,he_ub ! Helium Mass Flux at upper boundary
    real(r8),dimension(3) :: epep, &
      diff_fac ! correction factor for diffusion coefficients between He and O2, O, N2
    real(r8),dimension(3,3) :: invalpha
    real(r8),dimension(nlevp1) :: dtdz,dmdz,wks1, &
      eddyp,eddyq,eddyr,eddyppart,eddyrpart,eddyp1part,eddyr1part
    real(r8),dimension(3,nlevp1) :: &
      ep,fk,upd,dpdt,eddydif,veradv,loss_out,moldif
    real(r8),dimension(3,3,nlevp1) :: &
      alpha,molp,molq,molr,molp1,molr1,pk,qk,rk

    ! lower boundary condition
    real(r8),dimension(3) :: fb
    real(r8),dimension(3,3) :: b

    real(r8), parameter :: grav_cgs = 870._r8 ! waccmx altitudes (cgs units) cm/sec2
    real(r8), parameter :: p0 = 5.e-4_r8 ! cgs units

if (masterproc.and.debug) write(iulog,*) 'comp_wx: top of routine iCol,lchnk,he_ubc : ', iCol,lchnk,he_ubc

! N2, mbar at midpoint level 0 (not interface level 1)
    bn2 = max(1-bo2-bo1-bh-bhe,0.0_r8)
    bmbar = 1/(bo2/rmass_o2+ &
               bo1/rmass_o1+ &
               bh/rmass_h+   &
               bhe/rmass_he+ &
               bn2/rmass_n2)

    dtdz(1) = (tn(1)-tlbc)*2/dz(1)
    dmdz(1) = (mbar(1)-bmbar)/dz(1)
    do k = 2,nlevp1
      dtdz(k) = (tn(k)-tn(k-1))/dz(k)
      dmdz(k) = (mbar(k)-mbar(k-1))/dz(k)
    enddo

!if (masterproc) write(iulog,*) 'comp_wx: dtdz(1:10),dmdz(1:10), after k loop : ',&
!                                   dtdz(1:10),dtdz(1:10),dmdz(1:10),dmdz(1:10)

! WKS1 = MBAR/M4*(T00/T)**0.25/TAU
    wks1 = barm*(t00/tni)**0.25_r8/(tau*rmass_n2)

!if (masterproc) write(iulog,*) 'comp_wx: wks1(1:10), after calc : ',&
!                                   wks1(1:10)

! EP = 1-(M+DMBAR/DZ)/MBAR
    ep(1,:) = 1-(rmass_o2+dmdz)/barm
    ep(2,:) = 1-(rmass_o1+dmdz)/barm
    ep(3,:) = 1-(rmass_he+dmdz)/barm-thdiffalpha*dtdz/tni

!if (masterproc) write(iulog,*) 'comp_wx: ep(1,:),ep(2,:),ep(3,:) after calc : ',&
!                                   ep(1,:),ep(2,:),ep(3,:)

    do k = 1,nlevp1

! correction factors for mutual diffusion between He and O2, O, N2
      do n = 1,3
        diff_fac(n) = (tni(k)/t00)**(1.75_r8-ss(n))
      enddo

! alpha matrix
      alpha(1,1,k) = -psi(1,4)- &
        (psi(1,2)-psi(1,4))*o1i(k)- &
        (diff_fac(1)*psi(1,3)-psi(1,4))*hei(k)
      alpha(2,2,k) = -psi(2,4)- &
        (psi(2,1)-psi(2,4))*o2i(k)- &
        (diff_fac(2)*psi(2,3)-psi(2,4))*hei(k)
      alpha(3,3,k) = -diff_fac(3)*psi(3,4)- &
        (diff_fac(1)*psi(3,1)-diff_fac(3)*psi(3,4))*o2i(k)- &
        (diff_fac(2)*psi(3,2)-diff_fac(3)*psi(3,4))*o1i(k)
      alpha(1,2,k) = (psi(1,2)-psi(1,4))*o2i(k)
      alpha(1,3,k) = (diff_fac(1)*psi(1,3)-psi(1,4))*o2i(k)
      alpha(2,1,k) = (psi(2,1)-psi(2,4))*o1i(k)
      alpha(2,3,k) = (diff_fac(2)*psi(2,3)-psi(2,4))*o1i(k)
      alpha(3,1,k) = (diff_fac(1)*psi(3,1)-diff_fac(3)*psi(3,4))*hei(k)
      alpha(3,2,k) = (diff_fac(2)*psi(3,2)-diff_fac(3)*psi(3,4))*hei(k)

! molecular diffusion coefficients of O2, O, He
      invalpha = matinv3_wx(alpha(:,:,k))
      do n = 1,3
        do m = 1,3
          molp (m,n,k) = invalpha(m,n)*wks1(k)*(1/dz(k)+ep(n,k)/2)
          molr1(m,n,k) = invalpha(m,n)*wks1(k)*(1/dz(k)-ep(n,k)/2)
        enddo
      enddo
    enddo

!if (masterproc) write(iulog,*) 'comp_wx: o2 ep(1,:) after calc : ', ep(1,:)
!if (masterproc) write(iulog,*) 'comp_wx: o1 ep(2,:) after calc : ', ep(2,:)
!if (masterproc) write(iulog,*) 'comp_wx: he ep(3,:) after calc : ', ep(3,:)
!if (masterproc) write(iulog,*) 'comp_wx: o2 alpha(1,1,:) after calc : ', alpha(1,1,:)
!if (masterproc) write(iulog,*) 'comp_wx: o1 alpha(1,2,:) after calc : ', alpha(1,2,:)
!if (masterproc) write(iulog,*) 'comp_wx: he alpha(1,3,:) after calc : ', alpha(1,3,:)

! eddy diffusion coefficients (part) (difk=0)
    eddyppart  = difk*(1/dz-dmdz/(barm*2))
    eddyr1part = difk*(1/dz+dmdz/(barm*2))
!
    do k = 1,nlevp1-1
      molp1(:,:,k) = molp (:,:,k+1)
      molr (:,:,k) = molr1(:,:,k+1)
      eddyp1part(k) = eddyppart (k+1)
      eddyrpart (k) = eddyr1part(k+1)
    enddo
    molp1(:,:,nlevp1) = 2*molp (:,:,nlevp1)-molp (:,:,nlevp1-1)
    molr (:,:,nlevp1) = 2*molr1(:,:,nlevp1)-molr1(:,:,nlevp1-1)
    eddyp1part(nlevp1) =     2*eddyppart (nlevp1)-eddyppart (nlevp1-1)
    eddyrpart (nlevp1) = max(2*eddyr1part(nlevp1)-eddyr1part(nlevp1-1),0.0_r8)

    molq = molp1+molr1

!write(iulog,*) 'comp_wx: dfactor, expzmid, MIN/MAX eddyppart, before eddy : ',&
!   dfactor, expzmid, MINVAL(eddyppart(:)),MAXVAL(eddyppart(:))

! finish the remaining part of eddy diffusion coefficients (all zero)
    eddyp = dfactor*eddyppart/expzmid
    eddyr = dfactor*eddyrpart*expzmid
    eddyq = dfactor*(eddyp1part*expzmid+eddyr1part/expzmid)

!write(iulog,*) 'comp_wx: MIN/MAX expzm,delta,eddyp,wmid,dz before mol loop : ', MINVAL(expzm(:)), &
!                     MAXVAL(expzm(:)), &
!		     MINVAL(delta(:,:)),MAXVAL(delta(:,:)), &
!		     MINVAL(eddyp(:)),MAXVAL(eddyp(:)), &
!		     MINVAL(wmid(:)),MAXVAL(wmid(:)), &
!		     MINVAL(dz(:)),MAXVAL(dz(:))

    do n = 1,3
      do m = 1,3
        pk(m,n,:) = (molp(m,n,:)-expzm*delta(m,n)*(eddyp+wmid/2))/dz
        rk(m,n,:) = (molr(m,n,:)-expzm*delta(m,n)*(eddyr-wmid/2))/dz
        qk(m,n,:) = -molq(m,n,:)/dz+ &
          expzm*(delta(m,n)*(eddyq/dz+1/(2*step))-loss(m,n,:))
      enddo
    enddo

!write(iulog,*) 'comp_wx: MIN/MAX pk,molp after mol loop : ', MINVAL(pk(:,:,nlevp1-10:nlevp1)), &
!                     MAXVAL(pk(:,:,nlevp1-10:nlevp1)), &
!		     MINVAL(molp(:,:,nlevp1-10:nlevp1)),MAXVAL(molp(:,:,nlevp1-10:nlevp1))

! add explicit source terms to fk (no chemical production or advection)
    fk(1,:) = expzm*(prod(1,:)+o2i(:)/(2*step)-o2_hadv)
    fk(2,:) = expzm*(prod(2,:)+o1i(:)/(2*step)-o1_hadv)
    fk(3,:) = expzm*(prod(3,:)+hei(:)/(2*step)-he_hadv)

! lower boundaries

    call init_lbc(dz(1),b,fb)

    qk(:,:,1) = qk(:,:,1)+matmul(pk(:,:,1),b)
    do n = 1,3
      fk(n,1) = fk(n,1)-dot_product(pk(n,:,1),fb)
    enddo
    pk(:,:,1) = 0

! upper boundary
    epep = (2+ep(:,nlevp1)*dz(nlevp1))/(2-ep(:,nlevp1)*dz(nlevp1))
    do n = 1,3
      do m = 1,3
        qk(m,n,nlevp1-1) = qk(m,n,nlevp1-1)+epep(n)*rk(m,n,nlevp1-1)
      enddo
    enddo

if (masterproc.and.debug) write(iulog,*) 'comp_wx: iCol,lchnk,he_ubc before upper boundary calc : ', iCol,lchnk,he_ubc

! Eric Sutton: calculate Helium lateral exospheric transport mass flux at upper boundary
    flx00 = wks1(nlevp1)*p0/grav_cgs
    o1_ub = he_ubc*(alpha(2,3,nlevp1)-alpha(2,2,nlevp1))/(flx00*(1/dz(nlevp1)-ep(2,nlevp1)/2))
    he_ub = he_ubc*(alpha(3,3,nlevp1)-alpha(3,2,nlevp1))/(flx00*(1/dz(nlevp1)-ep(3,nlevp1)/2))
    fk(:,nlevp1-1) = fk(:,nlevp1-1)-rk(:,2,nlevp1-1)*o1_ub-rk(:,3,nlevp1-1)*he_ub
    rk(:,:,nlevp1-1) = 0

!if (masterproc) write(iulog,*) 'comp_wx: pk(1,1,:),qk(1,1,1:10),rk(1,1,:),fk(1,:) before blktri : ', &
!                     pk(1,1,:),qk(1,1,:),rk(1,1,:),fk(1,:)
!if (masterproc) write(iulog,*) 'comp_wx: pk(2,1,:),qk(2,1,1:10),rk(2,1,:),fk(2,:) before blktri : ', &
!                     pk(2,1,:),qk(2,1,:),rk(2,1,:),fk(2,:)
!if (masterproc) write(iulog,*) 'comp_wx: pk(3,1,:),qk(3,1,1:10),rk(3,1,:),fk(3,:) before blktri : ', &
!                     pk(3,1,:),qk(3,1,:),rk(3,1,:),fk(3,:)

!write(iulog,*) 'comp_wx: pk(1,1,1:10),qk(1,1,1:10),rk(1,1,1:10),fk(1:1:10) before blktri : ', &
!                     pk(1,1,1:10),qk(1,1,1:10),rk(1,1,1:10),fk(1,1:10)

if (debug) then
if (masterproc) write(iulog,*) 'comp_wx: pk,qk,rk,fk before blktri : ', &
                     MINVAL(pk(:,:,:)), MAXVAL(pk(:,:,:)), &
		     MINVAL(qk(:,:,:)),MAXVAL(qk(:,:,:)),&
		     MINVAL(rk(:,:,:)), MAXVAL(rk(:,:,:)), &
		     MINVAL(fk(:,:)), MAXVAL(fk(:,:))
end if
    upd = blktri_tgcm(pk,qk,rk,fk,nlevp1)

!if (masterproc) write(iulog,*) 'comp_wx: o2 upd(1,:) all levels after blktri ', upd(1,:)
!if (masterproc) write(iulog,*) 'comp_wx: o1 upd(2,:) all levels after blktri ', upd(2,:)
!if (masterproc) write(iulog,*) 'comp_wx: he upd(3,:) all levels after blktri ', upd(3,:)

! upper boundaries
    upd(:,nlevp1) = epep*upd(:,nlevp1-1)
    upd(2,nlevp1) = upd(2,nlevp1)+o1_ub
    upd(3,nlevp1) = upd(3,nlevp1)+he_ub

!    dpdt(1,:) = (upd(1,:)-o2_nm)/(2*step)
!    dpdt(2,:) = (upd(2,:)-o1_nm)/(2*step)
!    dpdt(3,:) = (upd(3,:)-he_nm)/(2*step)
!    do n = 1,3
!      do k = 1,nlevp1
!        loss_out(n,k) = dot_product(loss(n,:,k),upd(:,k))
!      enddo
!      do k = 2,nlevp1-1
!        moldif(n,k) = &
!          (dot_product(molp(n,:,k),upd(:,k-1))- &
!           dot_product(molq(n,:,k),upd(:,k  ))+ &
!           dot_product(molr(n,:,k),upd(:,k+1)))/dz/expzm(k)
!        eddydif(n,k) = &
!          (eddyp(k)*upd(n,k-1)- &
!           eddyq(k)*upd(n,k  )+ &
!           eddyr(k)*upd(n,k+1))/dz
!        veradv(n,k) = wmid(k)*(upd(n,k+1)-upd(n,k-1))/(2*dz)
!      enddo
!    enddo

! ensure non-negative O2, O, He
    o2_upd = max(upd(1,:),mmrMin)
    o1_upd = max(upd(2,:),mmrMin)
    he_upd = max(upd(3,:),HEmmrMin)
    he_upd = min(he_upd,HEmmrMax)
!! ensure non-negative O2, O, He
!    o2_upd = max(upd(1,:),0.0_r8)
!    o1_upd = max(upd(2,:),0.0_r8)
!    he_upd = max(upd(3,:),0.0_r8)

!if (upd(1,10) <= mmrMin) write(iulog,*) 'comp_wx: zero or negative o2 value at level 10, upd(1,:), o2_upd(10) ', upd(1,:), o2_upd(:)
!if (upd(2,10) <= mmrMin) write(iulog,*) 'comp_wx: zero or negative o1 value at level 10, upd(2,:), o1_upd(10) ', upd(2,:), o1_upd(:)


  endsubroutine comp_wx
!-----------------------------------------------------------------------

!===============================================================================
  subroutine mspdiff (lchnk      ,ncol       ,                                     &
                      t          ,q          ,pmid       ,pint       ,             &
                      pdel       ,ztodt      ,rairv      ,mbarv)
!-----------------------------------------------------------------------
! Driver routine to compute major species diffusion (O2 and O).

! Turbulent diffusivities and boundary layer nonlocal transport terms are
! obtained from the turbulence module.
!---------------------------Arguments------------------------------------
    use ref_pres,     only: lev0 => nbot_molec
    use physconst,    only: gravit

    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns
    real(r8), intent(in) :: t(pcols,pver)          ! temperature input
    real(r8), intent(in) :: pmid(pcols,pver)       ! midpoint pressures
    real(r8), intent(in) :: pint(pcols,pverp)      ! interface pressures
    real(r8), intent(in) :: pdel(pcols,pver)       ! thickness between interfaces
    real(r8), intent(in) :: ztodt                  ! 2 delta-t
    real(r8), intent(in) :: rairv(pcols,pver)                  ! composition dependent gas "constant"
    real(r8), intent(in) :: mbarv(pcols,pver)                  ! composition dependent mean mass

    real(r8), intent(inout) :: q(pcols,pver,pcnst) ! constituents

!---------------------------Local storage-------------------------------
    real(r8) :: o2(pcols,pver), o1(pcols,pver)     ! o2, o1 mixing ratio (kg/kg moist air)
    real(r8) :: h_atom(pcols,pver)                 ! H mixing ratio
    real(r8) :: rztodt                             ! 1/ztodt
    real(r8) :: dz(pcols,pver)                     ! log-pressure interval between interfaces
    real(r8) :: rdz(pcols,pver)                    ! 1./dz, defined on midpoints
    real(r8) :: dzmid(pcols,pverp)                 ! log-pressure interval between midpoints
    real(r8) :: rdzmid(pcols,pverp)                ! 1./dzmid, defined on interfaces
    real(r8) :: ak(pcols,2,2,2)                    ! coefficient matrix "Alfa"
    real(r8) :: ep(pcols,2,2)                      ! coefficient matrix
    real(r8) :: difk(pcols,pverp)                  ! eddy diffusion normalized by scale height (1/sec)
    real(r8) :: expzm(pcols,pver)                  ! exp(-z)=pmid/ptref
    real(r8) :: expzi(pcols,pverp)                 ! exp(-z)=pint/ptref
    real(r8) :: wks1(pcols)
    real(r8) :: wks3(pcols),wks4(pcols)            ! temporary working arrays
    real(r8) :: psclht(pcols,pverp)                ! pressure scale height
    real(r8) :: p_ubc(pcols)                       ! extrapolated pressure at upper boundary level, pmid(1)^2=pmid(2)*p_ubc
    real(r8) :: rair_ubc(pcols)                    ! extrapolated rair at upper boundary level
    real(r8) :: mbar_ubc(pcols)                    ! extrapolated mbar at upper boundary level
    real(r8) :: flb(pcols,2)                       ! lower boundary condition for o2 and o, now
                                                   ! calculated locally.
    real(r8) :: fub(pcols,2)                       ! upper boundary condition for o2 and o, now
                                                   ! calculated locally.
    real(r8) :: fk(pcols,2)                        ! temporary working array for rhs
    real(r8) :: pk(pcols,2,2)                      ! temp array for coefficients on lower diagonal
    real(r8) :: rk(pcols,2,2)                      ! temp array for coefficients on upper diagonal
    real(r8) :: qk(pcols,2,2)                      ! temp array for coefficients on diagonal
    real(r8) :: apk(2,2,pcols,pver)                ! coefficients on lower diagonal
    real(r8) :: ark(2,2,pcols,pver)                ! coefficients on upper diagonal
    real(r8) :: aqk(2,2,pcols,pver)                ! coefficients on diagonal
    real(r8) :: rfk(2,pcols,pver)                  ! rhs of the array equation
    real(r8) :: betawk(2,2,pcols,pver)             ! working arrays for blktri solver.
    real(r8) :: gammawk(2,2,pcols,pver)            ! working arrays for blktri solver.
    real(r8) :: ywk(2,pcols,pver)                  ! working arrays for blktri solver.
    real(r8) :: xwk(2,pcols,pver)                  ! working arrays for blktri solver.
    integer  :: nlevs                              ! number of levels
    real(r8) :: t_ubc(pcols)                       ! Temperature at top boundary
    integer  :: i, k, km, kp, m, ktmp, isp, kk, kr

    !---------------------------------------------------
    ! Set vertical grid and get time step for diffusion
    !---------------------------------------------------
    nlevs = lev0
    rztodt = 1._r8/ztodt

    !------------------------------------------------------
    ! Get species to diffuse and set upper/lower boundaries
    !------------------------------------------------------
    o2(:ncol,:) = q(:ncol,:,indx_O2)
    o1(:ncol,:) = q(:ncol,:,indx_O)
    h_atom(:ncol,:) = q(:ncol,:,indx_H)

    flb(:ncol,1) = o2(:ncol,lev0+1)      ! fixed lower boundary condition
    flb(:ncol,2) = o1(:ncol,lev0+1)
    if(fixed_ubc(io2).or.fixed_ubc(io1)) then
       fub(:ncol,1) = o2mmr_ubc(:ncol)      ! fixed upper boundary condition
       fub(:ncol,2) = ommr_ubc(:ncol)
    endif

    !------------------------------------------------------------------
    ! Get log-pressure intervals between midpoints and interface points
    !------------------------------------------------------------------
    dz(:ncol,:) = pdel(:ncol,:)/pmid(:ncol,:)
    rdz(:ncol,:) = 1._r8/dz(:ncol,:)
    do k=2,pver
       do i=1,ncol
          dzmid(i,k) = (pmid(i,k)-pmid(i,k-1))/pint(i,k)
          rdzmid(i,k) = 1._r8/dzmid(i,k)
       enddo
    enddo
    do i=1,ncol
       p_ubc(i) = pmid(i,1)*pmid(i,1)/pmid(i,2)
       dzmid(i,1) = (pmid(i,1)-p_ubc(i))/pint(i,1)
       rdzmid(i,1) = 1._r8/dzmid(i,1)
    enddo
    do i=1,ncol
       dzmid(i,pverp) = dzmid(i,pver)
       rdzmid(i,pverp) = 1._r8/dzmid(i,pverp)
    enddo

    !------------------------------------------------------------------
    ! Get log-pressure intervals between midpoints and interface points
    !------------------------------------------------------------------
    expzi(:ncol,:) = pint(:ncol,:)/ptref
    expzm(:ncol,:) = pmid(:ncol,:)/ptref

    !------------------------------------------------------------------
    ! Get pressure scale height
    !------------------------------------------------------------------
    do k=2,pver
       do i=1,ncol
          psclht(i,k) = .5_r8*(rairv(i,k)*t(i,k)+rairv(i,k-1)*t(i,k-1))/gravit
       enddo
    enddo
    do i=1,ncol
       rair_ubc(i) = 1.5_r8*rairv(i,1)-.5_r8*rairv(i,2)
       t_ubc(i) = 1.5_r8*t(i,1)-.5_r8*t(i,2)
       psclht(i,1) = .5_r8*(rairv(i,1)*t(i,1)+rair_ubc(i)*t_ubc(i))/gravit
       psclht(i,pverp) = psclht(i,pver)
    enddo

    !------------------------------------------------------------------
    ! Initialize scale height normalized eddy diffusion
    !------------------------------------------------------------------
    do k=1,pverp
       do i=1,ncol
          difk(i,k) = 0._r8          ! eddy diffusion already calculated in vertical_diffusion
       enddo
    enddo

    call outfld ('MBARV', mbarv(:,:), pcols, lchnk)

    !------------------------------------------------------------------
    ! Set up mean mass working array
    !------------------------------------------------------------------
    ! ep, ak at the interface level immediately below midpoint level nbot_molec

    ! WKS4 = .5*(DMBAR/DZ)/MBAR
    do i=1,ncol
       wks4(i) = (mbarv(i,lev0)-mbarv(i,lev0+1))/                               &
                 (dzmid(i,lev0+1)*(mbarv(i,lev0+1)+mbarv(i,lev0)))
    enddo

    !-----------------------------------
    ! Calculate coefficient matrices
    !-----------------------------------
    km = 1
    kp = 2
    do i=1, ncol
       ep(i,io2,kp) = 1._r8-(2._r8/(mbarv(i,lev0+1)+mbarv(i,lev0)))*                  &
                   (rmass_o2+(mbarv(i,lev0)-mbarv(i,lev0+1))*rdzmid(i,lev0+1))
       ep(i,io1,kp) = 1._r8-(2._r8/(mbarv(i,lev0+1)+mbarv(i,lev0)))*                  &
                   (rmass_o1+(mbarv(i,lev0)-mbarv(i,lev0+1))*rdzmid(i,lev0+1))
    enddo


    do m=1,2
      do i=1,ncol
         ak(i,io2,m,kp) =                                           &
            -delta(io2,m)*(phi(io1,3)+(phi(io1,io2)-phi(io1,3))*    &
            .5_r8*(o2(i,lev0+1)+o2(i,lev0)))-(1._r8-delta(io2,m))*        &
            (phi(io2,m)-phi(io2,3))*.5_r8*(o2(i,lev0+1)+o2(i,lev0))
         ak(i,io1,m,kp) =                                           &
            -delta(io1,m)*(phi(io2,3)+(phi(io2,io1)-phi(io2,3))*    &
            .5_r8*(o1(i,lev0+1)+o1(i,lev0)))-(1._r8-delta(io1,m))*        &
            (phi(io1,m)-phi(io1,3))*.5_r8*(o1(i,lev0+1)+o1(i,lev0))
      enddo
    enddo
!
! WKS1=MBAR/M3*(T00/(T0+T))*0.25/(TAU*DET(ak)) ak at the interface level
! immediately below midpoint level nbot_molec.
    do i=1,ncol
      wks1(i) = 0.5_r8*(mbarv(i,lev0+1)+mbarv(i,lev0))*rmassinv_n2*       &
        (2._r8*t00/(t(i,lev0+1)+t(i,lev0)))**0.25_r8/                      &
        (tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-ak(i,1,2,kp)*ak(i,2,1,kp)))
    enddo
!
! Complete calculation of ak at the interface level immediately below midpoint
! level nbot_molec.
    do m=1,2
      do i=1,ncol
        ak(i,io2,m,kp) = ak(i,io2,m,kp)*wks1(i)
        ak(i,io1,m,kp) = ak(i,io1,m,kp)*wks1(i)
      enddo
    enddo

    km = 1
    kp = 2
    do k=lev0,2,-1
       ktmp = km
       km = kp
       kp = ktmp
       do i=1,ncol
          ep(i,io2,kp) = 1._r8-(2._r8/(mbarv(i,k)+mbarv(i,k-1)))*(rmass_o2+      &
                         (mbarv(i,k-1)-mbarv(i,k))*rdzmid(i,k))
          ep(i,io1,kp) = 1._r8-(2._r8/(mbarv(i,k)+mbarv(i,k-1)))*(rmass_o1+      &
                         (mbarv(i,k-1)-mbarv(i,k))*rdzmid(i,k))
       enddo

       do m=1,2
          do i=1,ncol
             ak(i,io2,m,kp) =                                               &
                  -delta(io2,m)*(phi(io1,3)+(phi(io1,io2)-phi(io1,3))*      &
                  .5_r8*(o2(i,k)+o2(i,k-1)))-                                  &
                  (1._r8-delta(io2,m))*(phi(io2,m)-phi(io2,3))*                &
                  .5_r8*(o2(i,k)+o2(i,k-1))

             ak(i,io1,m,kp) =                                               &
                  -delta(io1,m)*(phi(io2,3)+(phi(io2,io1)-phi(io2,3))*      &
                  .5_r8*(o1(i,k)+o1(i,k-1)))-                                  &
                  (1._r8-delta(io1,m))*(phi(io1,m)-phi(io1,3))*                &
                  .5_r8*(o1(i,k)+o1(i,k-1))

          enddo
       enddo

    !---------------------------------------------
    ! Calculate coefficients for diagonals and rhs
    !---------------------------------------------
!
! WKS1=MBAR/M3*(T00/(T0+T))**0.25/(TAU*DET(ALFA))
       do i=1,ncol
          wks1(i) = 0.5_r8*(mbarv(i,k)+mbarv(i,k-1))*rmassinv_n2*              &
               (2._r8*t00/(t(i,k)+t(i,k-1)))**0.25_r8/                          &
               (tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-ak(i,1,2,kp)*ak(i,2,1,kp)))
          wks3(i) = wks4(i)
          wks4(i) = (mbarv(i,k-1)-mbarv(i,k))/                              &
               (dzmid(i,k)*(mbarv(i,k)+mbarv(i,k-1)))
       enddo

!
! FINISH CALCULATING AK(K+1/2) AND GENERATE PK, QK, RK
       do m=1,2
          do isp=io2,io1
             do i=1,ncol
                ak(i,isp,m,kp) = ak(i,isp,m,kp)*wks1(i)

                pk(i,isp,m) = (ak(i,isp,m,km)*(rdzmid(i,k+1)+ep(i,m,km)/2._r8)-   &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)-                  &
                     wks3(i))*delta(isp,m))*rdz(i,k)

                rk(i,isp,m) = (ak(i,isp,m,kp)*(rdzmid(i,k)-ep(i,m,kp)/2._r8)-     &
                     expzi(i,k)*difk(i,k)*(rdzmid(i,k)+                        &
                     wks4(i))*delta(isp,m))*rdz(i,k)

                qk(i,isp,m) = -(ak(i,isp,m,km)*(rdzmid(i,k+1)-ep(i,m,km)/2._r8)+  &
                     ak(i,isp,m,kp)*(rdzmid(i,k)+ep(i,m,kp)/2._r8))*rdz(i,k)+     &
                     ((expzi(i,k)*difk(i,k)*(rdzmid(i,k)-wks4(i))+             &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)+wks3(i)))*        &
                     rdz(i,k)+expzm(i,k)*rztodt)*delta(isp,m)

             enddo
          enddo
       enddo

       do i=1,ncol
          fk(i,io2) = expzm(i,k)*o2(i,k)*rztodt
          fk(i,io1) = expzm(i,k)*o1(i,k)*rztodt
       enddo

       !----------------------------
       ! Lower boundary
       !----------------------------
       if (k==lev0) then
          do m=1,2
            do i=1,ncol
              fk(i,io2) = fk(i,io2)-pk(i,io2,m)*flb(i,m)
              fk(i,io1) = fk(i,io1)-pk(i,io1,m)*flb(i,m)
              pk(i,:,m) = 0._r8
            enddo
          enddo
       endif

       kr = lev0-k+1

       do i=1,ncol
          do m=1,2
             do kk=1,2
                apk(kk,m,i,kr) = pk(i,kk,m)
                aqk(kk,m,i,kr) = qk(i,kk,m)
                ark(kk,m,i,kr) = rk(i,kk,m)
             enddo
          enddo
       enddo

       do i=1,ncol
          rfk(io2,i,kr) = fk(i,io2)
          rfk(io1,i,kr) = fk(i,io1)
       enddo

    enddo

    !----------------------------
    ! Upper boundary
    !----------------------------
    k=1
    ktmp = km
    km = kp
    kp = ktmp
    if(fixed_ubc(io2).or.fixed_ubc(io1)) then
       do i=1,ncol
          mbar_ubc(i) = 1._r8/(o2mmr_ubc(i)*rmassinv_o2+ommr_ubc(i)*rmassinv_o1+ &
               (1._r8-o2mmr_ubc(i)-ommr_ubc(i))*rmassinv_n2)
          ep(i,io2,kp) = 1._r8-(2._r8/(mbarv(i,k)+mbar_ubc(i)))*(rmass_o2+      &
               (mbar_ubc(i)-mbarv(i,k))*rdzmid(i,k))
          ep(i,io1,kp) = 1._r8-(2._r8/(mbarv(i,k)+mbar_ubc(i)))*(rmass_o1+      &
               (mbar_ubc(i)-mbarv(i,k))*rdzmid(i,k))
       enddo

       do m=1,2
          do i=1,ncol
             ak(i,io2,m,kp) =                                               &
                  -delta(io2,m)*(phi(io1,3)+(phi(io1,io2)-phi(io1,3))*      &
                  .5_r8*(o2(i,k)+o2mmr_ubc(i)))-                                  &
                  (1._r8-delta(io2,m))*(phi(io2,m)-phi(io2,3))*                &
                  .5_r8*(o2(i,k)+o2mmr_ubc(i))

             ak(i,io1,m,kp) =                                               &
                  -delta(io1,m)*(phi(io2,3)+(phi(io2,io1)-phi(io2,3))*      &
                  .5_r8*(o1(i,k)+ommr_ubc(i)))-                                  &
                  (1._r8-delta(io1,m))*(phi(io1,m)-phi(io1,3))*                &
                  .5_r8*(o1(i,k)+ommr_ubc(i))

          enddo
       enddo

!
! WKS1=MBAR/M3*(T00/(T0+T))**0.25/(TAU*DET(ALFA))
       do i=1,ncol
          wks1(i) = 0.5_r8*(mbarv(i,k)+mbar_ubc(i))*rmassinv_n2*              &
               (2._r8*t00/(t(i,k)+t_ubc(i)))**0.25_r8/                          &
               (tau*(ak(i,1,1,kp)*ak(i,2,2,kp)-ak(i,1,2,kp)*ak(i,2,1,kp)))
          wks3(i) = wks4(i)
          wks4(i) = (mbar_ubc(i)-mbarv(i,k))/                              &
               (dzmid(i,k)*(mbarv(i,k)+mbar_ubc(i)))
       enddo

!
! FINISH CALCULATING AK(K+1/2) AND GENERATE PK, QK, RK
       do m=1,2
          do isp=io2,io1
             do i=1,ncol
                ak(i,isp,m,kp) = ak(i,isp,m,kp)*wks1(i)

                pk(i,isp,m) = (ak(i,isp,m,km)*(rdzmid(i,k+1)+ep(i,m,km)/2._r8)-   &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)-                  &
                     wks3(i))*delta(isp,m))*rdz(i,k)

                rk(i,isp,m) = (ak(i,isp,m,kp)*(rdzmid(i,k)-ep(i,m,kp)/2._r8)-     &
                     expzi(i,k)*difk(i,k)*(rdzmid(i,k)+                        &
                     wks4(i))*delta(isp,m))*rdz(i,k)

                qk(i,isp,m) = -(ak(i,isp,m,km)*(rdzmid(i,k+1)-ep(i,m,km)/2._r8)+  &
                     ak(i,isp,m,kp)*(rdzmid(i,k)+ep(i,m,kp)/2._r8))*rdz(i,k)+     &
                     ((expzi(i,k)*difk(i,k)*(rdzmid(i,k)-wks4(i))+             &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)+wks3(i)))*        &
                     rdz(i,k)+expzm(i,k)*rztodt)*delta(isp,m)

             enddo
          enddo
       enddo

       do i=1,ncol
          fk(i,io2) = expzm(i,k)*o2(i,k)*rztodt
          fk(i,io1) = expzm(i,k)*o1(i,k)*rztodt
       enddo
       do m=1,2
          do i=1,ncol
             fk(i,io2) = fk(i,io2)-rk(i,io2,m)*fub(i,m)
             fk(i,io1) = fk(i,io1)-rk(i,io1,m)*fub(i,m)
             rk(i,:,m) = 0._r8
          enddo
       enddo

    else

       do i=1,ncol
          wks3(i) = wks4(i)
       enddo
       do m=1,2
          do isp=io2,io1
             do i=1,ncol
                pk(i,isp,m) = (ak(i,isp,m,km)*(rdzmid(i,k+1)+ep(i,m,km)/2._r8)-   &
                     expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)-                  &
                     wks3(i))*delta(isp,m))*rdz(i,k)

                qk(i,isp,m) = -(ak(i,isp,m,km)*(rdzmid(i,k+1)-ep(i,m,km)/2._r8))  &
                     *rdz(i,k)+     &
                     (expzi(i,k+1)*difk(i,k+1)*(rdzmid(i,k+1)+wks3(i))*        &
                     rdz(i,k)+expzm(i,k)*rztodt)*delta(isp,m)

             enddo
          enddo
       enddo

       do i=1,ncol
          fk(i,io2) = expzm(i,k)*o2(i,k)*rztodt
          fk(i,io1) = expzm(i,k)*o1(i,k)*rztodt
       enddo

    endif

    kr = lev0-k+1

    do i=1,ncol
       do m=1,2
          do kk=1,2
             apk(kk,m,i,kr) = pk(i,kk,m)
             aqk(kk,m,i,kr) = qk(i,kk,m)
             ark(kk,m,i,kr) = rk(i,kk,m)
          enddo
       enddo
    enddo

    do i=1,ncol
       rfk(io2,i,kr) = fk(i,io2)
       rfk(io1,i,kr) = fk(i,io1)
    enddo
if (debug) then
if (masterproc) write(iulog,*) 'mspdiff: apk,aqk,ark,rfk before blktri first column: ', &
                     MINVAL(apk(:,:,1,:)), MAXVAL(apk(:,:,1,:)), &
		     MINVAL(aqk(:,:,1,:)),MAXVAL(aqk(:,:,1,:)),&
		     MINVAL(ark(:,:,1,:)), MAXVAL(ark(:,:,1,:)), &
		     MINVAL(rfk(:,1,:)), MAXVAL(rfk(:,1,:))
end if
    !------------------------------------
    ! Call solver to get diffused species
    !------------------------------------
    call blktri(apk,aqk,ark,rfk,pcols,1,ncol,pver,1,nlevs,    &
                betawk, gammawk, ywk, xwk)
if (debug) then
if (masterproc) write(iulog,*) 'mspdiff: o2 xwk(1,1,:) all levels after blktri : ', xwk(1,1,:)
if (masterproc) write(iulog,*) 'mspdiff: o1 xwk(2,1,:) all levels after blktri : ', xwk(2,1,:)
end if
    do k=lev0,1,-1
       kr = lev0-k+1
       do i=1,ncol
          o2(i,k) = xwk(1,i,kr)
          o1(i,k) = xwk(2,i,kr)
       enddo
    enddo

    !---------------------------------------------------------------
    ! Ensure non-negative O2 and O and check for N2 greater than one
    !---------------------------------------------------------------
    do i=1,ncol
       do k=1,lev0

          if (o2(i,k) < mmrMin) o2(i,k) = mmrMin
          if (o1(i,k) < mmrMin) o1(i,k) = mmrMin

          if(1._r8-mmrMin-o2(i,k)-o1(i,k)-h_atom(i,k) < 0._r8) then
             o2(i,k) = o2(i,k)*((1._r8-N2mmrMin-h_atom(i,k))/(o2(i,k)+o1(i,k)))
             o1(i,k) = o1(i,k)*((1._r8-N2mmrMin-h_atom(i,k))/(o2(i,k)+o1(i,k)))
          endif
       enddo
    enddo

    q(:ncol,:,indx_O2) = o2(:ncol,:)
    q(:ncol,:,indx_O)  = o1(:ncol,:)

  end subroutine mspdiff


!===============================================================================
      SUBROUTINE BLKTRI(A,B,C,F,IF,I1,I2,KF,K1,K2,BETA,GAMMA,Y,X)
      implicit none
!     ****
!     ****     This procedure solves (I2-I1+1) tridiagonal block matrix
!     ****     systems in which all blocks are 2 x 2 matrices.
!     ****
!     ****     Each system may be written:
!     ****
!     ****      A(K) * X(K-1) + B(K) * X(K) + Z(K) * X(K+1) = F(K)
!     ****
!     ****      where:
!     ****
!     ****       K = K1,K2,1
!     ****
!     ****       A(K), B(K), C(K) are given (2 x 2) matrices.
!     ****
!     ****       The F(k) are given two componente vectors.
!     ****
!     ****       The system is to be solved for the two component
!     ****       vectors, X(K).
!     ****
!     ****       A(K1) = C(K2) = 0.
!     ****
!     ****      BETA(K), GAMMA(K), (K = K1,K2,1), are work space for
!     ****      (2 x 2) matrices.
!     ****
!     ****      Y(K), (K = K1,K2,1), is work space for two component
!     ****      vectors.
!     ****
!     ****     Algorithm: (See Isaacson and Keller p55)
!     ****
!     ****      Forward sweep from K = K1 to K = K2:
!     ****
!     ****       BETA(K1) = B(K1)**(-1)
!     ****
!     ****       Y(K1) = BETA(K1)*F(K1)
!     ****
!     ****       GAMMA(K) = BETA(K)*C(K),        K = K1,(K2-1),1
!     ****
!     ****       BETA(K) = (B(K) - A(K)*GAMMA(K-1))**(-1),
!     ****                                       K = K1+1,K2,1
!     ****
!     ****       Y(K) = BETA(K)*(F(K) - A(K)*Y(K-1)),
!     ****                                       K = K1+1,K2,1
!     ****
!     ****      Backward sweep, K = K2,K1,-1
!     ****
!     ****       X(K2) = Y(K2)
!     ****
!     ****       X(K) = Y(K) - GAMMA(K)*X(K+1),  K = K2-1,K1,-1
!     ****
!     ****     Dimension statements:
!     ****
!     ****      Block matrices are dimensioned thus:
!     ****
!     ****       MATRIX(2,2,IF,KF)
!     ****
!     ****      Two component vectors are similarly treated:
!     ****
!     ****       VECTOR(2,IF,KF)
!     ****
!     ****     Our block matrix scheme spans the range, (K = K1,K2,1),
!     ****     where (1 .LE. K1 .LT. K2 .LE. KF)
!     ****
!     ****     Similarly, we are solving (I1-I2+1) systems
!     ****     simultaneously as the index, I, spans the range,
!     ****     (I = I1,I2,1), where (1 .LE. I1 .LT. I2 .LE. IF)
!     ****
!     ****
!     ****     Dimension statements:
!     SUBROUTINE BLKTRI(A,B,C,F,IF,I1,I2,KF,K1,K2,BETA,GAMMA,Y,X)
!     ****
! Args:
      integer,intent(in)   :: if,i1,i2,kf,k1,k2
      real(r8),intent(in)  :: a(2,2,if,kf), b(2,2,if,kf), c(2,2,if,kf),      &
                              f(2,if,kf)
      real(r8),intent(out) :: beta(2,2,if,kf), gamma(2,2,if,kf),             &
                              y(2,if,kf), x(2,if,kf)
!
!     DIMENSION A(2,2,IF,KF), B(2,2,IF,KF), C(2,2,IF,KF), F(2,IF,KF),
!    1  BETA(2,2,IF,KF), GAMMA(2,2,IF,KF), Y(2,IF,KF), X(2,IF,KF)
!
! Local:
      integer :: i,k
!     ****
!     ****     Lower boundary at K = K1
!     ****
      DO I = I1,I2
!     ****
!     ****     Y(1,I,K1) = determinant(B(K1))
!     ****
        Y(1,I,K1) = B(1,1,I,K1)*B(2,2,I,K1) - B(1,2,I,K1)*B(2,1,I,K1)
!     ****
!     ****     BETA(K1) = B(K1)**(-1)
!     ****
        BETA(1,1,I,K1) = B(2,2,I,K1)/Y(1,I,K1)
        BETA(1,2,I,K1) = -B(1,2,I,K1)/Y(1,I,K1)
        BETA(2,1,I,K1) = -B(2,1,I,K1)/Y(1,I,K1)
        BETA(2,2,I,K1) = B(1,1,I,K1)/Y(1,I,K1)
!     ****
!     ****     Y(K1) = BETA(K1)*F(K1)
!     ****
        Y(1,I,K1) = BETA(1,1,I,K1)*F(1,I,K1) + BETA(1,2,I,K1)*F(2,I,K1)
        Y(2,I,K1) = BETA(2,1,I,K1)*F(1,I,K1) + BETA(2,2,I,K1)*F(2,I,K1)
      ENDDO
!     ****
!     ****     Now deal with levels (K1+1),K2,1
!     ****
      DO K = K1+1,K2
        DO I = I1,I2
!         ****
!         ****     GAMMA(K-1) = BETA(K-1)*C(K-1)
!         ****
          GAMMA(1,1,I,K-1) = BETA(1,1,I,K-1)*C(1,1,I,K-1) +              &
            BETA(1,2,I,K-1)*C(2,1,I,K-1)
          GAMMA(1,2,I,K-1) = BETA(1,1,I,K-1)*C(1,2,I,K-1) +              &
            BETA(1,2,I,K-1)*C(2,2,I,K-1)
          GAMMA(2,1,I,K-1) = BETA(2,1,I,K-1)*C(1,1,I,K-1) +              &
            BETA(2,2,I,K-1)*C(2,1,I,K-1)
          GAMMA(2,2,I,K-1) = BETA(2,1,I,K-1)*C(1,2,I,K-1) +              &
            BETA(2,2,I,K-1)*C(2,2,I,K-1)
!         ****
!         ****     GAMMA(K) = B(K) - A(K)*GAMMA(K-1)
!         ****
          GAMMA(1,1,I,K) = B(1,1,I,K) - A(1,1,I,K)*GAMMA(1,1,I,K-1) -    &
            A(1,2,I,K)*GAMMA(2,1,I,K-1)
          GAMMA(1,2,I,K) = B(1,2,I,K) - A(1,1,I,K)*GAMMA(1,2,I,K-1) -    &
            A(1,2,I,K)*GAMMA(2,2,I,K-1)
          GAMMA(2,1,I,K) = B(2,1,I,K) - A(2,1,I,K)*GAMMA(1,1,I,K-1) -    &
            A(2,2,I,K)*GAMMA(2,1,I,K-1)
          GAMMA(2,2,I,K) = B(2,2,I,K) - A(2,1,I,K)*GAMMA(1,2,I,K-1) -    &
            A(2,2,I,K)*GAMMA(2,2,I,K-1)
!         ****
!         ****     Y(1,I,K) = determinant(GAMMA(K))
!         ****
          Y(1,I,K) = GAMMA(1,1,I,K)*GAMMA(2,2,I,K) -                     &
            GAMMA(1,2,I,K)*GAMMA(2,1,I,K)
!         ****
!         ****     BETA(K) = GAMMA(K)**(-1)
!         ****
          BETA(1,1,I,K) = GAMMA(2,2,I,K)/Y(1,I,K)
          BETA(1,2,I,K) = -GAMMA(1,2,I,K)/Y(1,I,K)
          BETA(2,1,I,K) = -GAMMA(2,1,I,K)/Y(1,I,K)
          BETA(2,2,I,K) = GAMMA(1,1,I,K)/Y(1,I,K)
!         ****
!         ****     X(K) = F(K) - A(K)*Y(K-1)
!         ****
          X(1,I,K) = F(1,I,K) - A(1,1,I,K)*Y(1,I,K-1) -                  &
            A(1,2,I,K)*Y(2,I,K-1)
          X(2,I,K) = F(2,I,K) - A(2,1,I,K)*Y(1,I,K-1) -                  &
            A(2,2,I,K)*Y(2,I,K-1)
!         ****
!         ****     Y(K) = BETA(K)*X(K)
!         ****
          Y(1,I,K) = BETA(1,1,I,K)*X(1,I,K) + BETA(1,2,I,K)*X(2,I,K)
          Y(2,I,K) = BETA(2,1,I,K)*X(1,I,K) + BETA(2,2,I,K)*X(2,I,K)
        ENDDO
      ENDDO
!     ****
!     ****     Backward sweep to determine final solution, X(K) for
!     ****     K = K2,K1,-1
!     ****
!     ****      X(K2) = Y(K2)
!     ****
      DO I = I1,I2
        X(1,I,K2) = Y(1,I,K2)
        X(2,I,K2) = Y(2,I,K2)
      ENDDO
!     ****
!     ****      X(K) = Y(K) - GAMMA(K)*X(K+1)
!     ****
      DO K = K2-1,K1,-1
        DO I = I1,I2
          X(1,I,K) = Y(1,I,K) - GAMMA(1,1,I,K)*X(1,I,K+1) -             &
            GAMMA(1,2,I,K)*X(2,I,K+1)
          X(2,I,K) = Y(2,I,K) - GAMMA(2,1,I,K)*X(1,I,K+1) -             &
            GAMMA(2,2,I,K)*X(2,I,K+1)

        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE BLKTRI

  !-------------------------------------------------------------------
  pure function matinv3_wx(A) result(B)
  ! Calculate the inverse of the matrix

    real(r8),dimension(3,3),intent(in) :: A
    real(r8),dimension(3,3) :: B

    B = matadj3_wx(A)/matdet3_wx(A)

  endfunction matinv3_wx

!-------------------------------------------------------------------
  pure function matadj3_wx(A) result(B)
! Calculate the adjugate of the matrix

    real(r8),dimension(3,3),intent(in) :: A
    real(r8),dimension(3,3) :: B

    B(1,1) =  (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) =  (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) =  (A(1,1)*A(2,2) - A(1,2)*A(2,1))

  endfunction matadj3_wx

!-------------------------------------------------------------------
  pure function matdet3_wx(A) result(d)
! Calculate the determinant of the matrix

    real(r8),dimension(3,3),intent(in) :: A
    real(r8) :: d

    d = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
      - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
      + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

  endfunction matdet3_wx
!-------------------------------------------------------------------

  pure subroutine init_lbc(dz, b, fb)

!    use matutil_module,only:matinv3

    real(r8),intent(in) :: dz
    ! lower boundary condition out
    real(r8),intent(out),dimension(3) :: fb
    real(r8),intent(out),dimension(3,3) :: b

    real(r8),parameter :: &
      alfa = 0.234_r8, &    ! lower boundary for O2+O (0.22+0.14)
      pshelb = 0.1154e-5_r8 ! lower boundary for Helium (mmr)
    real(r8),dimension(3),parameter :: &
      g = -(/alfa,0.0_r8,pshelb/) ! g = -(O2+O 0 He)
    real(r8),dimension(3,3),parameter :: &
!     |0 0 0|
! e = |0 1 0|
!     |0 0 0|
      e = reshape((/0,0,0,0,1,0,0,0,0/),(/3,3/)), &
!     |1  1  0|
! f = |0 -1  0|
!     |0  0  1|
      f = reshape((/1,0,0,1,-1,0,0,0,1/),(/3,3/))
    integer :: n
    real(r8),dimension(3,3) :: wm1,wm2,wm3

! calculate matrix b(3,3) and vector fb(3)
! representing the lower boundary condition in major,
! where psi = (O2 O He) are calculated as
! psi(k=-1/2) = b * psi(k=1/2) + fb

! first define 3x3 matrices e, f and length-3 vector g
! in the general lower boundary condition
! e * d(psi)/ds + f * psi + g = 0

! then evaluates b and fb from:
! b = (e/ds - f/2)**(-1) * (e/ds + f/2)
! fb = (e/ds - f/2)**(-1) * g

! wm1 = (e/ds - f/2)
! wm2 = (e/ds + f/2)
    wm1 = e/dz - f/2
    wm2 = e/dz + f/2

! now invert wm1 in wm3
    wm3 = matinv3_wx(wm1)

! b = wm3 * wm2
    b = matmul(wm3,wm2)

! fb = wm3 * g
    do n = 1,3
      fb(n) = dot_product(wm3(n,:),g)
    enddo

  endsubroutine init_lbc

!-----------------------------------------------------------------------
  pure function blktri_tgcm(pk,qk,rk,fk,nk) result(upd)

!    use matutil_module,only:matinv3

    integer,intent(in) :: nk
    real(r8),dimension(3,3,nk),intent(in) :: pk,qk,rk
    real(r8),dimension(3,nk),intent(in) :: fk
    real(r8),dimension(3,nk) :: upd

    integer :: n,k
    real(r8),dimension(3) :: wkv1
    real(r8),dimension(3,3) :: wkm1
    real(r8),dimension(3,nk) :: zz
    real(r8),dimension(3,3,nk) :: gama

    zz(:,1) = 0
    gama(:,:,1) = 0

    do k = 1,nk-1
! ALFA = Q(K)-P(K)*GAMA(K)
! ALFA refers to the block diagonal matrices,
!   and GAMA to the upper block diagonal matrices
!   in the Thomas algorithm solution
!   to the block tridiagonal system of equations
! WKM1 = INV(ALFA)
      wkm1 = matinv3_wx(qk(:,:,k)-matmul(pk(:,:,k),gama(:,:,k)))

! WKV1 = F(K)-P(K)*Z(K)
      do n = 1,3
        wkv1(n) = fk(n,k)-dot_product(pk(n,:,k),zz(:,k))
      enddo

! GAMA(K+1) = WKM1*R(K)
      gama(:,:,k+1) = matmul(wkm1,rk(:,:,k))

! Z(K+1) = WKM1*WKV1
      do n = 1,3
        zz(n,k+1) = dot_product(wkm1(n,:),wkv1)
      enddo
    enddo

! set upper boundary to zero
    upd(:,nk) = 0

! downward sweep
    do k = nk-1,1,-1
      do n = 1,3
        upd(n,k) = zz(n,k+1)-dot_product(gama(n,:,k+1),upd(:,k+1))
      enddo
    enddo

  endfunction blktri_tgcm
!-----------------------------------------------------------------------

end module majorsp_diffusion
