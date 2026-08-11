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

!-----------------------
! Public interfaces
!-----------------------
  public mspd_init   ! Initialization
  public mspd_intr   ! Full routine

!-----------------------
! Private data
!-----------------------

  real(r8) :: rmass_o2, rmass_o1, rmass_h, rmass_he, rmass_n2 ! molecular weights (kg/kmol)

  real(r8), parameter :: ptref=5.e-5_r8                  ! thermosphere reference pressure (Pa)
  real(r8), parameter :: mmrMin=1.e-20_r8                ! lower limit of o2 and o mixing ratio
  real(r8), parameter :: N2mmrMin=1.e-6_r8               ! lower limit of n2 mixing ratios
  real(r8), parameter :: HEmmrMin=1.e-7_r8               ! lower limit of he mixing ratios
  real(r8), parameter :: HEmmrMax=0.9_r8                 ! upper limit of he mixing ratios

  integer :: indx_O2                                     ! cnst index for o2
  integer :: indx_O                                      ! cnst index for o
  integer :: indx_H                                      ! cnst index for h
  integer :: indx_HE                                     ! cnst index for he
  integer, parameter :: io2=1, io1=2, ihe=3              ! local indices to o2 , o, and he respectively
  logical :: fixed_ubc(2)                                ! flag for fixed upper boundary condition

  character(len=8) :: mjdiffnam(3)              ! names of v-diff tendencies

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

    !--------------------------------------------------------------------
    ! Get fixed upper boundary flags and set vertical range for diffusion
    !--------------------------------------------------------------------
    fixed_ubc(io2) = cnst_fixed_ubc(indx_O2)
    fixed_ubc(io1) = cnst_fixed_ubc(indx_O)

   ! Set names of major diffusion tendencies and declare them as history variables
    mjdiffnam(1) = 'MD'//trim(cnst_name(indx_O2))
    call addfld (mjdiffnam(1),(/ 'lev' /), 'A','kg/kg/s','Major diffusion of '//cnst_name(indx_O2))
    mjdiffnam(2) = 'MD'//trim(cnst_name(indx_O))
    call addfld (mjdiffnam(2),(/ 'lev' /), 'A','kg/kg/s','Major diffusion of '//cnst_name(indx_O))
    mjdiffnam(3) = 'MD'//trim(cnst_name(indx_HE))
    call addfld (mjdiffnam(3),(/ 'lev' /), 'A','kg/kg/s','Major diffusion of '//cnst_name(indx_HE))

    call addfld('MOLP11', (/ 'lev' /), 'A',' ','major species diffusion ...')
    call addfld('MOLP22', (/ 'lev' /), 'A',' ','major species diffusion ...')
    call addfld('MOLP33', (/ 'lev' /), 'A',' ','major species diffusion ...')
    call addfld('MOLQ11', (/ 'lev' /), 'A',' ','major species diffusion ...')
    call addfld('MOLQ22', (/ 'lev' /), 'A',' ','major species diffusion ...')
    call addfld('MOLQ33', (/ 'lev' /), 'A',' ','major species diffusion ...')
    call addfld('MOLR11', (/ 'lev' /), 'A',' ','major species diffusion ...')
    call addfld('MOLR22', (/ 'lev' /), 'A',' ','major species diffusion ...')
    call addfld('MOLR33', (/ 'lev' /), 'A',' ','major species diffusion ...')

    call addfld ('MBARV' , (/ 'lev' /),'I','g/mole','Variable Mean Mass')

    if (history_waccmx) then
       call add_default (mjdiffnam(1), 1, ' ')
       call add_default (mjdiffnam(2), 1, ' ')
       call add_default (mjdiffnam(3), 1, ' ')
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
    real(r8) :: tendo2ohe(pcols,pver,3)            ! temporary array for o2 o, and he tendencies
    real(r8) :: ubc_mmr(pcols,pcnst)               ! upper bndy mixing ratios (kg/kg)
    real(r8) :: ubc_t(pcols)                       ! upper bndy temperature (K)
    integer :: lchnk                               ! chunk identifier
    integer :: ncol                                ! number of atmospheric columns
    integer :: i, k, kk, icol                      ! indexing integers

    ! For comp_wx call

    real(r8) :: tlbc,bo2,bo1,bhe,bh,he_ubc,p_ubc     ! For lower boundary
    real(r8) :: step,dfactor,pscaleheight,expzmid
    real(r8),dimension(nbot_molec) :: &
      difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
      o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm,dz,expzm
    real(r8),dimension(3,nbot_molec) :: prod
    real(r8),dimension(3,3,nbot_molec) :: loss
    real(r8),dimension(nbot_molec) :: o2_upd,o1_upd,he_upd
    real(r8),dimension(nbot_molec) :: &
         molp11,molp22,molp33, &
         molq11,molq22,molq33, &
         molr11,molr22,molr33
    real(r8),dimension(pcols,pver) :: o2_upd_cols,o1_upd_cols,he_upd_cols,h_upd_cols
    real(r8),dimension(pcols,pver) :: o2_upd_cols_tend,o1_upd_cols_tend,he_upd_cols_tend
    real(r8),dimension(pcols,pver) :: &
         molp11_cols,molp22_cols,molp33_cols, &
         molq11_cols,molq22_cols,molq33_cols, &
         molr11_cols,molr22_cols,molr33_cols

    real(r8) :: o2mmr_ubc(pcols) ! MMR of O2 at top boundary (specified)
    real(r8) :: ommr_ubc(pcols)  ! MMR of O at top boundary
    real(r8) :: heflx_ubc(pcols) ! MMR flux of HE at top boundary

    !--------------------------------------------------------------------------------------------
    ! local constants
    !--------------------------------------------------------------------------------------------
    rztodt = 1._r8/ztodt
    lchnk = state%lchnk
    ncol  = state%ncol

    molp11=0._r8
    molp22=0._r8
    molp33=0._r8
    molq11=0._r8
    molq22=0._r8
    molq33=0._r8
    molr11=0._r8
    molr22=0._r8
    molr33=0._r8

    molp11_cols=0._r8
    molp22_cols=0._r8
    molp33_cols=0._r8
    molq11_cols=0._r8
    molq22_cols=0._r8
    molq33_cols=0._r8
    molr11_cols=0._r8
    molr22_cols=0._r8
    molr33_cols=0._r8

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

    heflx_ubc(:ncol) = helium_ubc_fluxes(:ncol,lchnk)

    ! Since this is a combined tendency, retain the old name for output
    ! and debugging purposes.
    ptend%name  = trim(ptend%name)//"+mspd"
    ptend%lq(indx_O2) = .TRUE.
    ptend%lq(indx_O) = .TRUE.

    step = ztodt

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

    tn = nan
    tni = nan

    do iCol = 1,ncol

      tlbc   = state%t(iCol,nbot_molec+1)
      bo2    = state%q(iCol,nbot_molec+1,indx_O2)
      bo1    = state%q(iCol,nbot_molec+1,indx_O)
      bhe    = state%q(iCol,nbot_molec+1,indx_HE)
      bh     = state%q(iCol,nbot_molec+1,indx_H)
      he_ubc = heflx_ubc(iCol)

      kk = 0
      do k = nbot_molec,2,-1

        kk = kk + 1
        tn(kk)       = state%t(iCol,k)
        tni(kk)      = .5_r8 * (state%t(iCol,k) + state%t(iCol,k-1))
        o2i(kk)      = .5_r8 * (state%q(iCol,k,indx_O2) + state%q(iCol,k-1,indx_O2))
        o1i(kk)      = .5_r8 * (state%q(iCol,k,indx_O) + state%q(iCol,k-1,indx_O))
        hei(kk)      = .5_r8 * (state%q(iCol,k,indx_HE) + state%q(iCol,k-1,indx_HE))

        o2_nm(kk) = state%q(iCol,k,indx_O2)
        o1_nm(kk) = state%q(iCol,k,indx_O)
        he_nm(kk) = state%q(iCol,k,indx_HE)

        mbar(kk)     = mbarv(iCol,k,lchnk)
        barm(kk)     = .5_r8 * (mbarv(iCol,k,lchnk) + mbarv(iCol,k-1,lchnk))
        pScaleHeight = .5_r8*(rairv(iCol,k,lchnk)*state%t(iCol,k) + rairv(iCol,k-1,lchnk)*state%t(iCol,k)) / gravit
        dz(kk)       = (state%pmid(iCol,k) - state%pmid(iCol,k-1)) / state%pint(iCol,k)
        expzm(kk)    = state%pmid(iCol,k) / ptref

      enddo ! kk=1,nbot_molec-1
      !
      ! Top:
      !
      tn(nbot_molec)     = state%t(iCol,1)
      tni(nbot_molec)    = 1.5_r8*state%t(iCol,1)-.5_r8*state%t(iCol,2)
      o2i(nbot_molec)    = 1.5_r8*state%q(iCol,1,indx_O2)-.5_r8*state%q(iCol,2,indx_O2)
      o1i(nbot_molec)    = 1.5_r8*state%q(iCol,1,indx_O)-.5_r8*state%q(iCol,2,indx_O)
      hei(nbot_molec)    = 1.5_r8*state%q(iCol,1,indx_HE)-.5_r8*state%q(iCol,2,indx_HE)
      mbar(nbot_molec)   = mbarv(iCol,1,lchnk)
      barm(nbot_molec)   = 1.5_r8*mbarv(iCol,1,lchnk)-.5_r8*mbarv(iCol,2,lchnk)
      pScaleHeight       = .5_r8*(rairv(iCol,1,lchnk)*state%t(iCol,1) + rairv(iCol,2,lchnk)*state%t(iCol,2)) / gravit
      p_ubc             = state%pmid(iCol,1)*state%pmid(iCol,1)/state%pmid(iCol,2)
      dz(nbot_molec)    = (state%pmid(iCol,1)-p_ubc)/state%pint(iCol,1)
      expzm(nbot_molec) = state%pmid(iCol,1) / ptref

      ! extrapolate to layer about top
      o2_nm(nbot_molec) = 2._r8*state%q(iCol,1,indx_O2) - state%q(iCol,2,indx_O2)
      o1_nm(nbot_molec) = 2._r8*state%q(iCol,1,indx_O ) - state%q(iCol,2,indx_O )
      he_nm(nbot_molec) = 2._r8*state%q(iCol,1,indx_HE) - state%q(iCol,2,indx_HE)

      call comp_wx(step,dfactor,tlbc,bo2,bo1,bh,bhe,he_ubc,difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
               o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm,prod,loss, &
               nbot_molec,dz,expzm,expzmid,o2_upd,o1_upd,he_upd, &
               molp11,molp22,molp33, molq11,molq22,molq33, molr11,molr22,molr33 )

       kk = 0
       do k = 1,nbot_molec

         kk = nbot_molec - k + 1

         o2_upd_cols(iCol,kk) = o2_upd(k)
         o1_upd_cols(iCol,kk) = o1_upd(k)
         he_upd_cols(iCol,kk) = he_upd(k)

         molp11_cols(icol,kk) = molp11(k)
         molp22_cols(icol,kk) = molp22(k)
         molp22_cols(icol,kk) = molp22(k)

         molq11_cols(icol,kk) = molq11(k)
         molq22_cols(icol,kk) = molq22(k)
         molq22_cols(icol,kk) = molq22(k)

         molr11_cols(icol,kk) = molr11(k)
         molr22_cols(icol,kk) = molr22(k)
         molr22_cols(icol,kk) = molr22(k)

       enddo

    enddo ! iCol loop

    call outfld('MOLP11',molp11_cols(:pcols,:),pcols,lchnk)
    call outfld('MOLP22',molp22_cols(:pcols,:),pcols,lchnk)
    call outfld('MOLP33',molp33_cols(:pcols,:),pcols,lchnk)
    call outfld('MOLQ11',molq11_cols(:pcols,:),pcols,lchnk)
    call outfld('MOLQ22',molq22_cols(:pcols,:),pcols,lchnk)
    call outfld('MOLQ33',molq33_cols(:pcols,:),pcols,lchnk)
    call outfld('MOLR11',molr11_cols(:pcols,:),pcols,lchnk)
    call outfld('MOLR22',molr22_cols(:pcols,:),pcols,lchnk)
    call outfld('MOLR33',molr33_cols(:pcols,:),pcols,lchnk)

    !---------------------------------------------------------------
    ! Check for N2 greater than one
    !---------------------------------------------------------------
    do i=1,ncol
       do k=1,nbot_molec

          if(1._r8-mmrMin-o2_upd_cols(i,k)-o1_upd_cols(i,k)-he_upd_cols(i,k)-h_upd_cols(i,k) < 0._r8) then
             o2_upd_cols(i,k) = o2_upd_cols(i,k)*((1._r8-N2mmrMin-h_upd_cols(i,k)) &
                                                 /(o2_upd_cols(i,k)+o1_upd_cols(i,k)+he_upd_cols(i,k)))
             o1_upd_cols(i,k) = o1_upd_cols(i,k)*((1._r8-N2mmrMin-h_upd_cols(i,k)) &
                                                 /(o2_upd_cols(i,k)+o1_upd_cols(i,k)+he_upd_cols(i,k)))
             he_upd_cols(i,k) = he_upd_cols(i,k)*((1._r8-N2mmrMin-h_upd_cols(i,k)) &
                                                 /(o2_upd_cols(i,k)+o1_upd_cols(i,k)+he_upd_cols(i,k)))
          endif

       enddo
    enddo

    !---------------------------------------------
    ! Update O2 and O tendencies and output
    !---------------------------------------------
    do k=1,pver
       do i=1,ncol

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

    call outfld(mjdiffnam(1),ptend%q(1,1,indx_O2),pcols,lchnk)
    call outfld(mjdiffnam(2),ptend%q(1,1,indx_O),pcols,lchnk)
    call outfld(mjdiffnam(3),ptend%q(1,1,indx_HE),pcols,lchnk)

  end subroutine mspd_intr

!-----------------------------------------------------------------------
  subroutine comp_wx(step,dfactor,tlbc,bo2,bo1,bh,bhe,he_ubc, &
    difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
    o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm, &
    prod,loss,nlevp1,dz,expzm,expzmid,o2_upd,o1_upd,he_upd, &
    molp11,molp22,molp33, molq11,molq22,molq33, molr11,molr22,molr33 )

! advance major species O2, O, He and N2

    integer,intent(in) :: nlevp1

    real(r8),intent(in) :: step,dfactor,tlbc,bo2,bo1,bh,bhe,he_ubc,expzmid
    real(r8),dimension(nlevp1),intent(in) :: &
      difk,tn,tni,o2i,o1i,hei,wmid,mbar,barm, &
      o2_hadv,o1_hadv,he_hadv,o2_nm,o1_nm,he_nm,dz,expzm
    real(r8),dimension(3,nlevp1),intent(inout) :: prod
    real(r8),dimension(3,3,nlevp1),intent(inout) :: loss
    real(r8),dimension(nlevp1),intent(out) :: o2_upd,o1_upd,he_upd
    real(r8),dimension(nlevp1),intent(out) :: &
         molp11,molp22,molp33, molq11,molq22,molq33, molr11,molr22,molr33

    ! exponent factor for diff_fac
    real(r8),dimension(3),parameter :: ss = (/1.710_r8,1.749_r8,1.718_r8/)

    ! mutual thermal diffusion coefficients among major species
    real(r8),dimension(3,4),parameter :: &
      psi = reshape( &
       (/0.0_r8,0.673_r8,0.270_r8, &
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

! WKS1 = MBAR/M4*(T00/T)**0.25/TAU
    wks1 = barm*(t00/tni)**0.25_r8/(tau*rmass_n2)

! EP = 1-(M+DMBAR/DZ)/MBAR
    ep(1,:) = 1-(rmass_o2+dmdz)/barm
    ep(2,:) = 1-(rmass_o1+dmdz)/barm
    ep(3,:) = 1-(rmass_he+dmdz)/barm-thdiffalpha*dtdz/tni

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

! finish the remaining part of eddy diffusion coefficients (all zero)
    eddyp = dfactor*eddyppart/expzmid
    eddyr = dfactor*eddyrpart*expzmid
    eddyq = dfactor*(eddyp1part*expzmid+eddyr1part/expzmid)

    do n = 1,3
      do m = 1,3
        pk(m,n,:) = (molp(m,n,:)-expzm*delta(m,n)*(eddyp+wmid/2))/dz
        rk(m,n,:) = (molr(m,n,:)-expzm*delta(m,n)*(eddyr-wmid/2))/dz
        qk(m,n,:) = -molq(m,n,:)/dz+ &
          expzm*(delta(m,n)*(eddyq/dz+1/(step))-loss(m,n,:))
      enddo
    enddo

    molp11(:) = molp(1,1,:)
    molp22(:) = molp(2,2,:)
    molp33(:) = molp(3,3,:)
    molq11(:) = molq(1,1,:)
    molq22(:) = molq(2,2,:)
    molq33(:) = molq(3,3,:)
    molr11(:) = molr(1,1,:)
    molr22(:) = molr(2,2,:)
    molr33(:) = molr(3,3,:)

! add explicit source terms to fk (no chemical production or advection)
    fk(1,:) = expzm*(prod(1,:)+o2_nm(:)/(step)-o2_hadv)
    fk(2,:) = expzm*(prod(2,:)+o1_nm(:)/(step)-o1_hadv)
    fk(3,:) = expzm*(prod(3,:)+he_nm(:)/(step)-he_hadv)

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

! Eric Sutton: calculate Helium lateral exospheric transport mass flux at upper boundary
    flx00 = wks1(nlevp1)*p0/grav_cgs
    o1_ub = he_ubc*(alpha(2,3,nlevp1)-alpha(2,2,nlevp1))/(flx00*(1/dz(nlevp1)-ep(2,nlevp1)/2))
    he_ub = he_ubc*(alpha(3,3,nlevp1)-alpha(3,2,nlevp1))/(flx00*(1/dz(nlevp1)-ep(3,nlevp1)/2))
    fk(:,nlevp1-1) = fk(:,nlevp1-1)-rk(:,2,nlevp1-1)*o1_ub-rk(:,3,nlevp1-1)*he_ub
    rk(:,:,nlevp1-1) = 0

    upd = blktri_tgcm(pk,qk,rk,fk,nlevp1)

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

  end subroutine comp_wx
!-----------------------------------------------------------------------

  !-------------------------------------------------------------------
  pure function matinv3_wx(A) result(B)
  ! Calculate the inverse of the matrix

    real(r8),dimension(3,3),intent(in) :: A
    real(r8),dimension(3,3) :: B

    B = matadj3_wx(A)/matdet3_wx(A)

  end function matinv3_wx

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

  end function matadj3_wx

!-------------------------------------------------------------------
  pure function matdet3_wx(A) result(d)
! Calculate the determinant of the matrix

    real(r8),dimension(3,3),intent(in) :: A
    real(r8) :: d

    d = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
      - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
      + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

  end function matdet3_wx
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

  end function blktri_tgcm
!-----------------------------------------------------------------------

end module majorsp_diffusion
