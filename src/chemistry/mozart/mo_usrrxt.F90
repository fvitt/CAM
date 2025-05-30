module mo_usrrxt

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use cam_logfile,      only : iulog
  use ppgrid,           only : pver, pcols

  implicit none

  private
  public :: usrrxt, usrrxt_inti, usrrxt_hrates

  save

  integer :: usr_O_O2_ndx
  integer :: usr_HO2_HO2_ndx
  integer :: usr_N2O5_M_ndx
  integer :: usr_HNO3_OH_ndx
  integer :: usr_HO2NO2_M_ndx
  integer :: usr_N2O5_aer_ndx
  integer :: usr_NO3_aer_ndx
  integer :: usr_NO2_aer_ndx
  integer :: usr_CO_OH_a_ndx
  integer :: usr_CO_OH_b_ndx
  integer :: usr_PAN_M_ndx
  integer :: usr_CH3COCH3_OH_ndx
  integer :: usr_MCO3_NO2_ndx
  integer :: usr_MPAN_M_ndx
  integer :: usr_XOOH_OH_ndx
  integer :: usr_SO2_OH_ndx
  integer :: usr_DMS_OH_ndx
  integer :: usr_HO2_aer_ndx
  integer :: usr_GLYOXAL_aer_ndx
  
  integer :: tag_NO2_NO3_ndx
  integer :: tag_NO2_OH_ndx
  integer :: tag_NO2_HO2_ndx
  integer :: tag_C2H4_OH_ndx
  integer :: tag_C3H6_OH_ndx
  integer :: tag_CH3CO3_NO2_ndx
  
!lke-TS1
  integer :: usr_PBZNIT_M_ndx
  integer :: tag_ACBZO2_NO2_ndx
  integer :: usr_ISOPNITA_aer_ndx
  integer :: usr_ISOPNITB_aer_ndx
  integer :: usr_ONITR_aer_ndx
  integer :: usr_HONITR_aer_ndx
  integer :: usr_TERPNIT_aer_ndx
  integer :: usr_NTERPOOH_aer_ndx
  integer :: usr_NC4CHO_aer_ndx
  integer :: usr_NC4CH2OH_aer_ndx
!
  integer :: usr_OA_O2_NDX
  integer :: usr_XNO2NO3_M_ndx
  integer :: usr_NO2XNO3_M_ndx
  integer :: usr_XHNO3_OH_ndx
  integer :: usr_XHO2NO2_M_ndx
  integer :: usr_XNO2NO3_aer_ndx
  integer :: usr_NO2XNO3_aer_ndx
  integer :: usr_XNO3_aer_ndx
  integer :: usr_XNO2_aer_ndx
  integer :: usr_XPAN_M_ndx
  integer :: usr_XMPAN_M_ndx
  integer :: usr_MCO3_XNO2_ndx
  
  integer :: usr_C2O3_NO2_ndx
  integer :: usr_C2H4_OH_ndx
  integer :: usr_XO2N_HO2_ndx
  integer :: usr_C2O3_XNO2_ndx
  
  integer :: tag_XO2N_NO_ndx
  integer :: tag_XO2_HO2_ndx
  integer :: tag_XO2_NO_ndx
  
  integer :: usr_O_O_ndx
  integer :: usr_CL2O2_M_ndx
  integer :: usr_SO3_H2O_ndx
  integer :: tag_CLO_CLO_M_ndx
  
  integer :: ion1_ndx, ion2_ndx, ion3_ndx, ion11_ndx
  integer :: elec1_ndx, elec2_ndx, elec3_ndx
  integer :: elec4_ndx, elec5_ndx, elec6_ndx
  integer :: het1_ndx

  integer, parameter :: nean = 3
  integer :: ean_ndx(nean)
  integer, parameter :: nrpe = 5
  integer :: rpe_ndx(nrpe)
  integer, parameter :: npir = 16
  integer :: pir_ndx(npir)
  integer, parameter :: nedn = 2
  integer :: edn_ndx(nedn)
  integer, parameter :: nnir = 13
  integer :: nir_ndx(nnir)  
  integer, parameter :: niira = 112
  integer :: iira_ndx(niira)
  integer, parameter :: niirb = 14
  integer :: iirb_ndx(niirb)

  integer :: usr_clm_h2o_m_ndx, usr_clm_hcl_m_ndx
  integer :: usr_oh_co_ndx, het_no2_h2o_ndx, usr_oh_dms_ndx, aq_so2_h2o2_ndx, aq_so2_o3_ndx

  integer :: h2o_ndx
!
! jfl
!
  integer, parameter :: num_strat_tau = 22
  integer :: usr_strat_tau_ndx(num_strat_tau)
!WUHU FENG, JOHN PLANE,DAN MARSH, 18/09/2018
  integer :: na02_ndx,na08_ndx, na09_ndx,na10_ndx,na12_ndx,na22_ndx
  integer :: na23_ndx,na29a_ndx, na29b_ndx,na29c_ndx,na29d_ndx,na30_ndx
  integer :: fe04_ndx
  integer :: fe09_ndx, fe22_ndx, fe23_ndx
  integer :: fe25a_ndx,fe25b_ndx, fe25c_ndx, fe25d_ndx
  integer :: fe40_ndx, fe41_ndx,fe42_ndx,fe43_ndx,fe44_ndx,fe45_ndx
  integer :: mg02_ndx,mg04_ndx,mg05_ndx,mg07_ndx,mg08_ndx,mg10_ndx
  integer :: mg20_ndx,mg22_ndx,mg27_ndx,mg28_ndx,mg29_ndx
  integer :: mg30_ndx,mg30b_ndx,mg31c_ndx,mg31d_ndx

  integer :: si01_ndx,si04_ndx,si05_ndx,si06_ndx,si07_ndx
  integer :: si12_ndx,si13_ndx,si14_ndx
  integer :: si16_ndx,si20b_ndx,si20c_ndx,si20d_ndx,si21_ndx
  integer :: si23_ndx

  integer :: ca04_ndx,ca06_ndx,ca08_ndx
  integer :: ca09_ndx, ca12_ndx,ca13_ndx
  integer :: ca17b_ndx
  integer :: ca23_ndx, ca25_ndx,ca26_ndx
  integer :: ca31_ndx, ca31b_ndx,ca31c_ndx,ca31d_ndx, ca32_ndx
  integer :: ca34_ndx

   integer :: k03_ndx,k11_ndx
  integer :: k22_ndx,k23_ndx,k24_ndx,k25_ndx
  integer :: k26_ndx,k27_ndx,k27b_ndx,k28_ndx,k29_ndx
  integer :: k30_ndx,k31_ndx,k32_ndx,k33_ndx,k34_ndx
  integer :: k35_ndx,k36_ndx,k37a_ndx,k37b_ndx,k37c_ndx
  integer :: k37d_ndx,k37e_ndx



!WF

!
!lke++
  integer :: usr_COhc_OH_ndx
  integer :: usr_COme_OH_ndx
  integer :: usr_CO01_OH_ndx
  integer :: usr_CO02_OH_ndx
  integer :: usr_CO03_OH_ndx
  integer :: usr_CO04_OH_ndx
  integer :: usr_CO05_OH_ndx
  integer :: usr_CO06_OH_ndx
  integer :: usr_CO07_OH_ndx
  integer :: usr_CO08_OH_ndx
  integer :: usr_CO09_OH_ndx
  integer :: usr_CO10_OH_ndx
  integer :: usr_CO11_OH_ndx
  integer :: usr_CO12_OH_ndx
  integer :: usr_CO13_OH_ndx
  integer :: usr_CO14_OH_ndx
  integer :: usr_CO15_OH_ndx
  integer :: usr_CO16_OH_ndx
  integer :: usr_CO17_OH_ndx
  integer :: usr_CO18_OH_ndx
  integer :: usr_CO19_OH_ndx
  integer :: usr_CO20_OH_ndx
  integer :: usr_CO21_OH_ndx
  integer :: usr_CO22_OH_ndx
  integer :: usr_CO23_OH_ndx
  integer :: usr_CO24_OH_ndx
  integer :: usr_CO25_OH_ndx
  integer :: usr_CO26_OH_ndx
  integer :: usr_CO27_OH_ndx
  integer :: usr_CO28_OH_ndx
  integer :: usr_CO29_OH_ndx
  integer :: usr_CO30_OH_ndx
  integer :: usr_CO31_OH_ndx
  integer :: usr_CO32_OH_ndx
  integer :: usr_CO33_OH_ndx
  integer :: usr_CO34_OH_ndx
  integer :: usr_CO35_OH_ndx
  integer :: usr_CO36_OH_ndx
  integer :: usr_CO37_OH_ndx
  integer :: usr_CO38_OH_ndx
  integer :: usr_CO39_OH_ndx
  integer :: usr_CO40_OH_ndx
  integer :: usr_CO41_OH_ndx
  integer :: usr_CO42_OH_ndx
!lke--

  real(r8), parameter :: t0     = 300._r8                ! K
  real(r8), parameter :: trlim2 = 17._r8/3._r8           ! K
  real(r8), parameter :: trlim3 = 15._r8/3._r8           ! K

  logical :: has_ion_rxts, has_d_chem

contains

  subroutine usrrxt_inti
    !-----------------------------------------------------------------
    !        ... intialize the user reaction constants module
    !-----------------------------------------------------------------

    use mo_chem_utls,   only : get_rxt_ndx, get_spc_ndx
    use spmd_utils,     only : masterproc

    implicit none

    character(len=4) :: xchar
    character(len=32) :: rxtname
    integer :: i

!
! full tropospheric chemistry
!
    usr_O_O2_ndx         = get_rxt_ndx( 'usr_O_O2' )
    usr_HO2_HO2_ndx      = get_rxt_ndx( 'usr_HO2_HO2' )
    usr_N2O5_M_ndx       = get_rxt_ndx( 'usr_N2O5_M' )
    usr_HNO3_OH_ndx      = get_rxt_ndx( 'usr_HNO3_OH' )
    usr_HO2NO2_M_ndx     = get_rxt_ndx( 'usr_HO2NO2_M' )
    usr_N2O5_aer_ndx     = get_rxt_ndx( 'usr_N2O5_aer' )
    usr_NO3_aer_ndx      = get_rxt_ndx( 'usr_NO3_aer' )
    usr_NO2_aer_ndx      = get_rxt_ndx( 'usr_NO2_aer' )
    usr_CO_OH_a_ndx      = get_rxt_ndx( 'usr_CO_OH_a' )
    usr_CO_OH_b_ndx      = get_rxt_ndx( 'usr_CO_OH_b' )
    usr_PAN_M_ndx        = get_rxt_ndx( 'usr_PAN_M' )
    usr_CH3COCH3_OH_ndx  = get_rxt_ndx( 'usr_CH3COCH3_OH' )
    usr_MCO3_NO2_ndx     = get_rxt_ndx( 'usr_MCO3_NO2' )
    usr_MPAN_M_ndx       = get_rxt_ndx( 'usr_MPAN_M' )
    usr_XOOH_OH_ndx      = get_rxt_ndx( 'usr_XOOH_OH' )
    usr_SO2_OH_ndx       = get_rxt_ndx( 'usr_SO2_OH' )
    usr_DMS_OH_ndx       = get_rxt_ndx( 'usr_DMS_OH' )
    usr_HO2_aer_ndx      = get_rxt_ndx( 'usr_HO2_aer' )
    usr_GLYOXAL_aer_ndx  = get_rxt_ndx( 'usr_GLYOXAL_aer' )
 !
    tag_NO2_NO3_ndx      = get_rxt_ndx( 'tag_NO2_NO3' )
    tag_NO2_OH_ndx       = get_rxt_ndx( 'tag_NO2_OH' )
    tag_NO2_HO2_ndx      = get_rxt_ndx( 'tag_NO2_HO2' )
    tag_C2H4_OH_ndx      = get_rxt_ndx( 'tag_C2H4_OH' )
    tag_C3H6_OH_ndx      = get_rxt_ndx( 'tag_C3H6_OH' )
    tag_CH3CO3_NO2_ndx   = get_rxt_ndx( 'tag_CH3CO3_NO2' )     
!lke-TS1
    usr_PBZNIT_M_ndx     = get_rxt_ndx( 'usr_PBZNIT_M' )
    tag_ACBZO2_NO2_ndx   = get_rxt_ndx( 'tag_ACBZO2_NO2' )
    usr_ISOPNITA_aer_ndx = get_rxt_ndx( 'usr_ISOPNITA_aer' )
    usr_ISOPNITB_aer_ndx = get_rxt_ndx( 'usr_ISOPNITB_aer' )
    usr_ONITR_aer_ndx    = get_rxt_ndx( 'usr_ONITR_aer' )
    usr_HONITR_aer_ndx   = get_rxt_ndx( 'usr_HONITR_aer' )
    usr_TERPNIT_aer_ndx  = get_rxt_ndx( 'usr_TERPNIT_aer' )
    usr_NTERPOOH_aer_ndx = get_rxt_ndx( 'usr_NTERPOOH_aer' )
    usr_NC4CHO_aer_ndx   = get_rxt_ndx( 'usr_NC4CHO_aer' )
    usr_NC4CH2OH_aer_ndx = get_rxt_ndx( 'usr_NC4CH2OH_aer' )
 !
 ! additional reactions for O3A/XNO
 !
    usr_OA_O2_ndx        = get_rxt_ndx( 'usr_OA_O2' )
    usr_XNO2NO3_M_ndx    = get_rxt_ndx( 'usr_XNO2NO3_M' )
    usr_NO2XNO3_M_ndx    = get_rxt_ndx( 'usr_NO2XNO3_M' )
    usr_XNO2NO3_aer_ndx  = get_rxt_ndx( 'usr_XNO2NO3_aer' )
    usr_NO2XNO3_aer_ndx  = get_rxt_ndx( 'usr_NO2XNO3_aer' )
    usr_XHNO3_OH_ndx     = get_rxt_ndx( 'usr_XHNO3_OH' )
    usr_XNO3_aer_ndx     = get_rxt_ndx( 'usr_XNO3_aer' )
    usr_XNO2_aer_ndx     = get_rxt_ndx( 'usr_XNO2_aer' )
    usr_MCO3_XNO2_ndx    = get_rxt_ndx( 'usr_MCO3_XNO2' )
    usr_XPAN_M_ndx       = get_rxt_ndx( 'usr_XPAN_M' )
    usr_XMPAN_M_ndx      = get_rxt_ndx( 'usr_XMPAN_M' )
    usr_XHO2NO2_M_ndx    = get_rxt_ndx( 'usr_XHO2NO2_M' )
!
! reduced hydrocarbon chemistry
!
    usr_C2O3_NO2_ndx     = get_rxt_ndx( 'usr_C2O3_NO2' )
    usr_C2H4_OH_ndx      = get_rxt_ndx( 'usr_C2H4_OH' )
    usr_XO2N_HO2_ndx     = get_rxt_ndx( 'usr_XO2N_HO2' )
    usr_C2O3_XNO2_ndx    = get_rxt_ndx( 'usr_C2O3_XNO2' )
!
    tag_XO2N_NO_ndx      = get_rxt_ndx( 'tag_XO2N_NO' )
    tag_XO2_HO2_ndx      = get_rxt_ndx( 'tag_XO2_HO2' )
    tag_XO2_NO_ndx       = get_rxt_ndx( 'tag_XO2_NO' )
!
! stratospheric chemistry
!
    usr_O_O_ndx          = get_rxt_ndx( 'usr_O_O' )
    usr_CL2O2_M_ndx      = get_rxt_ndx( 'usr_CL2O2_M' )
    usr_SO3_H2O_ndx      = get_rxt_ndx( 'usr_SO3_H2O' )
!    
    tag_CLO_CLO_M_ndx      = get_rxt_ndx( 'tag_CLO_CLO_M' )
    if (tag_CLO_CLO_M_ndx<0) then ! for backwards compatibility
       tag_CLO_CLO_M_ndx   = get_rxt_ndx( 'tag_CLO_CLO' )
    endif
!
! reactions to remove BAM aerosols in the stratosphere
!
    usr_strat_tau_ndx( 1) = get_rxt_ndx( 'usr_CB1_strat_tau' )
    usr_strat_tau_ndx( 2) = get_rxt_ndx( 'usr_CB2_strat_tau' )
    usr_strat_tau_ndx( 3) = get_rxt_ndx( 'usr_OC1_strat_tau' )
    usr_strat_tau_ndx( 4) = get_rxt_ndx( 'usr_OC2_strat_tau' )
    usr_strat_tau_ndx( 5) = get_rxt_ndx( 'usr_SO4_strat_tau' )
    usr_strat_tau_ndx( 6) = get_rxt_ndx( 'usr_SOA_strat_tau' )
    usr_strat_tau_ndx( 7) = get_rxt_ndx( 'usr_NH4_strat_tau' )
    usr_strat_tau_ndx( 8) = get_rxt_ndx( 'usr_NH4NO3_strat_tau' )
    usr_strat_tau_ndx( 9) = get_rxt_ndx( 'usr_SSLT01_strat_tau' )
    usr_strat_tau_ndx(10) = get_rxt_ndx( 'usr_SSLT02_strat_tau' )
    usr_strat_tau_ndx(11) = get_rxt_ndx( 'usr_SSLT03_strat_tau' )
    usr_strat_tau_ndx(12) = get_rxt_ndx( 'usr_SSLT04_strat_tau' )
    usr_strat_tau_ndx(13) = get_rxt_ndx( 'usr_DST01_strat_tau' )
    usr_strat_tau_ndx(14) = get_rxt_ndx( 'usr_DST02_strat_tau' )
    usr_strat_tau_ndx(15) = get_rxt_ndx( 'usr_DST03_strat_tau' )
    usr_strat_tau_ndx(16) = get_rxt_ndx( 'usr_DST04_strat_tau' )
    usr_strat_tau_ndx(17) = get_rxt_ndx( 'usr_SO2t_strat_tau' )
    usr_strat_tau_ndx(18) = get_rxt_ndx( 'usr_SOAI_strat_tau' )
    usr_strat_tau_ndx(19) = get_rxt_ndx( 'usr_SOAM_strat_tau' )
    usr_strat_tau_ndx(20) = get_rxt_ndx( 'usr_SOAB_strat_tau' )
    usr_strat_tau_ndx(21) = get_rxt_ndx( 'usr_SOAT_strat_tau' )
    usr_strat_tau_ndx(22) = get_rxt_ndx( 'usr_SOAX_strat_tau' )
!
! stratospheric aerosol chemistry
!
    het1_ndx             = get_rxt_ndx( 'het1' )
!
! ion chemistry
!   
    ion1_ndx  = get_rxt_ndx( 'ion_Op_O2' )
    ion2_ndx  = get_rxt_ndx( 'ion_Op_N2' )
    ion3_ndx  = get_rxt_ndx( 'ion_N2p_Oa' )
    ion11_ndx = get_rxt_ndx( 'ion_N2p_Ob' )

    elec1_ndx  = get_rxt_ndx( 'elec1' )
    elec2_ndx  = get_rxt_ndx( 'elec2' )
    elec3_ndx  = get_rxt_ndx( 'elec3' )

    do i = 1,nean
      write (xchar,'(i4)') i
      rxtname = 'ean'//trim(adjustl(xchar))
      ean_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,nrpe
      write (xchar,'(i4)') i
      rxtname = 'rpe'//trim(adjustl(xchar))
      rpe_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,npir
      write (xchar,'(i4)') i
      rxtname = 'pir'//trim(adjustl(xchar))
      pir_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,nedn
      write (xchar,'(i4)') i
      rxtname = 'edn'//trim(adjustl(xchar))
      edn_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,nnir
      write (xchar,'(i4)') i
      rxtname = 'nir'//trim(adjustl(xchar))
      nir_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,niira
      write (xchar,'(i4)') i
      rxtname = 'iira'//trim(adjustl(xchar))
      iira_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,niirb
      write (xchar,'(i4)') i
      rxtname = 'iirb'//trim(adjustl(xchar))
      iirb_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    usr_clm_h2o_m_ndx = get_rxt_ndx( 'usr_CLm_H2O_M' )
    usr_clm_hcl_m_ndx = get_rxt_ndx( 'usr_CLm_HCL_M' )

    elec4_ndx  = get_rxt_ndx( 'Op2P_ea' )
    elec5_ndx  = get_rxt_ndx( 'Op2P_eb' )
    elec6_ndx  = get_rxt_ndx( 'Op2D_e' )

    has_ion_rxts = ion1_ndx>0 .and. ion2_ndx>0 .and. ion3_ndx>0 .and. elec1_ndx>0 &
                 .and. elec2_ndx>0 .and. elec3_ndx>0

    has_d_chem = &
         all(ean_ndx>0) .and. &
         all(rpe_ndx>0) .and. &
         all(pir_ndx>0) .and. &
         all(edn_ndx>0) .and. &
         all(nir_ndx>0) .and. &
         all(iira_ndx>0) .and. &
         all(iirb_ndx>0) .and. &
         usr_clm_h2o_m_ndx>0 .and. usr_clm_hcl_m_ndx>0

    h2o_ndx    = get_spc_ndx( 'H2O' )

    !
    ! llnl super fast
    !
    usr_oh_co_ndx  = get_rxt_ndx( 'usr_oh_co' )
    het_no2_h2o_ndx  = get_rxt_ndx( 'het_no2_h2o' )
    usr_oh_dms_ndx  = get_rxt_ndx( 'usr_oh_dms' )
    aq_so2_h2o2_ndx  = get_rxt_ndx( 'aq_so2_h2o2' )
    aq_so2_o3_ndx  = get_rxt_ndx( 'aq_so2_o3' )
!WUHU FENG, JOHN PLANE,DAN MARSH, 18/09/2018
    na02_ndx    = get_rxt_ndx( 'na02' )
    na08_ndx    = get_rxt_ndx( 'na08' )
    na09_ndx    = get_rxt_ndx( 'na09' )
    na10_ndx    = get_rxt_ndx( 'na10' )
    na12_ndx    = get_rxt_ndx( 'na12' )
    na22_ndx    = get_rxt_ndx( 'na22' )
    na23_ndx    = get_rxt_ndx( 'na23' )
    na29a_ndx    = get_rxt_ndx( 'na29a' )
    na29b_ndx    = get_rxt_ndx( 'na29b' )
    na29c_ndx    = get_rxt_ndx( 'na29c' )
    na29d_ndx    = get_rxt_ndx( 'na29d' )
    na30_ndx    = get_rxt_ndx( 'na30' )
    fe04_ndx    = get_rxt_ndx( 'fe04' )
    fe09_ndx    = get_rxt_ndx( 'fe09' )
    fe22_ndx    = get_rxt_ndx( 'fe22' )
    fe23_ndx    = get_rxt_ndx( 'fe23' )
    fe25a_ndx    = get_rxt_ndx( 'fe25a' )
    fe25b_ndx    = get_rxt_ndx( 'fe25b' )
    fe25c_ndx    = get_rxt_ndx( 'fe25c' )
    fe25d_ndx    = get_rxt_ndx( 'fe25d' )
    fe40_ndx    = get_rxt_ndx( 'fe40' )
    fe41_ndx    = get_rxt_ndx( 'fe41' )
    fe42_ndx    = get_rxt_ndx( 'fe42' )
    fe43_ndx    = get_rxt_ndx( 'fe43' )
    fe44_ndx    = get_rxt_ndx( 'fe44' )
    fe45_ndx    = get_rxt_ndx( 'fe45' )
    mg02_ndx    = get_rxt_ndx( 'mg02' )
    mg04_ndx    = get_rxt_ndx( 'mg04' )
    mg05_ndx    = get_rxt_ndx( 'mg05' )
    mg07_ndx    = get_rxt_ndx( 'mg07' )
    mg08_ndx    = get_rxt_ndx( 'mg08' )
    mg10_ndx    = get_rxt_ndx( 'mg10' )
    mg20_ndx    = get_rxt_ndx( 'mg20' )
    mg22_ndx    = get_rxt_ndx( 'mg22' )
    mg27_ndx    = get_rxt_ndx( 'mg27' )
    mg28_ndx    = get_rxt_ndx( 'mg28' )
    mg29_ndx    = get_rxt_ndx( 'mg29' )
    mg30_ndx    = get_rxt_ndx( 'mg30' )
    mg30b_ndx    = get_rxt_ndx( 'mg30b' )
    mg31c_ndx    = get_rxt_ndx( 'mg31c' )
    mg31d_ndx    = get_rxt_ndx( 'mg31d' )
    si01_ndx    = get_rxt_ndx( 'si01' )
    si04_ndx    = get_rxt_ndx( 'si04' )
    si05_ndx    = get_rxt_ndx( 'si05' )
    si06_ndx    = get_rxt_ndx( 'si06' )
    si07_ndx    = get_rxt_ndx( 'si07' )
    si12_ndx    = get_rxt_ndx( 'si12' )
    si13_ndx    = get_rxt_ndx( 'si13' )
    si14_ndx    = get_rxt_ndx( 'si14' )
    si16_ndx    = get_rxt_ndx( 'si16' )
    si20b_ndx    = get_rxt_ndx( 'si20b' )
    si20c_ndx    = get_rxt_ndx( 'si20c' )
    si20d_ndx    = get_rxt_ndx( 'si20d' )
    si21_ndx    = get_rxt_ndx( 'si21' )
    si23_ndx    = get_rxt_ndx( 'si23' )

    ca04_ndx    = get_rxt_ndx( 'ca04' )
    ca06_ndx    = get_rxt_ndx( 'ca06' )
    ca08_ndx    = get_rxt_ndx( 'ca08' )
    ca09_ndx    = get_rxt_ndx( 'ca09' )
!   ca12_ndx    = get_rxt_ndx( 'ca12' )
    ca13_ndx    = get_rxt_ndx( 'ca13' )
    ca17b_ndx    = get_rxt_ndx( 'ca17b' )
    ca23_ndx    = get_rxt_ndx( 'ca23' )
    ca25_ndx    = get_rxt_ndx( 'ca25' )
    ca26_ndx    = get_rxt_ndx( 'ca26' )
    ca31_ndx    = get_rxt_ndx( 'ca31' )
    ca31b_ndx    = get_rxt_ndx( 'ca31b' )
    ca31c_ndx    = get_rxt_ndx( 'ca31c' )
    ca31d_ndx    = get_rxt_ndx( 'ca31d' )
    ca32_ndx    = get_rxt_ndx( 'ca32' )
    ca34_ndx    = get_rxt_ndx( 'ca34' )

    k03_ndx    = get_rxt_ndx( 'k03' )
    k11_ndx    = get_rxt_ndx( 'k11' )
    k22_ndx    = get_rxt_ndx( 'k22' )
    k23_ndx    = get_rxt_ndx( 'k23' )
    k24_ndx    = get_rxt_ndx( 'k24' )
    k25_ndx    = get_rxt_ndx( 'k25' )
    k26_ndx    = get_rxt_ndx( 'k26' )
    k27_ndx    = get_rxt_ndx( 'k27' )
    k27b_ndx    = get_rxt_ndx( 'k27b' )
    k28_ndx    = get_rxt_ndx( 'k28' )
    k29_ndx    = get_rxt_ndx( 'k29' )
    k30_ndx    = get_rxt_ndx( 'k30' )
    k31_ndx    = get_rxt_ndx( 'k31' )
    k32_ndx    = get_rxt_ndx( 'k32' )
    k33_ndx    = get_rxt_ndx( 'k33' )
    k34_ndx    = get_rxt_ndx( 'k34' )
    k35_ndx    = get_rxt_ndx( 'k35' )
    k36_ndx    = get_rxt_ndx( 'k36' )
    k37a_ndx    = get_rxt_ndx( 'k37a' )
    k37b_ndx    = get_rxt_ndx( 'k37b' )
    k37c_ndx    = get_rxt_ndx( 'k37c' )
    k37d_ndx    = get_rxt_ndx( 'k37d' )
    k37e_ndx    = get_rxt_ndx( 'k37e' )



!WF

    
!lke++
! CO tags
!
    usr_COhc_OH_ndx      = get_rxt_ndx( 'usr_COhc_OH' )
    usr_COme_OH_ndx      = get_rxt_ndx( 'usr_COme_OH' )
    usr_CO01_OH_ndx      = get_rxt_ndx( 'usr_CO01_OH' )
    usr_CO02_OH_ndx      = get_rxt_ndx( 'usr_CO02_OH' )
    usr_CO03_OH_ndx      = get_rxt_ndx( 'usr_CO03_OH' )
    usr_CO04_OH_ndx      = get_rxt_ndx( 'usr_CO04_OH' )
    usr_CO05_OH_ndx      = get_rxt_ndx( 'usr_CO05_OH' )
    usr_CO06_OH_ndx      = get_rxt_ndx( 'usr_CO06_OH' )
    usr_CO07_OH_ndx      = get_rxt_ndx( 'usr_CO07_OH' )
    usr_CO08_OH_ndx      = get_rxt_ndx( 'usr_CO08_OH' )
    usr_CO09_OH_ndx      = get_rxt_ndx( 'usr_CO09_OH' )
    usr_CO10_OH_ndx      = get_rxt_ndx( 'usr_CO10_OH' )
    usr_CO11_OH_ndx      = get_rxt_ndx( 'usr_CO11_OH' )
    usr_CO12_OH_ndx      = get_rxt_ndx( 'usr_CO12_OH' )
    usr_CO13_OH_ndx      = get_rxt_ndx( 'usr_CO13_OH' )
    usr_CO14_OH_ndx      = get_rxt_ndx( 'usr_CO14_OH' )
    usr_CO15_OH_ndx      = get_rxt_ndx( 'usr_CO15_OH' )
    usr_CO16_OH_ndx      = get_rxt_ndx( 'usr_CO16_OH' )
    usr_CO17_OH_ndx      = get_rxt_ndx( 'usr_CO17_OH' )
    usr_CO18_OH_ndx      = get_rxt_ndx( 'usr_CO18_OH' )
    usr_CO19_OH_ndx      = get_rxt_ndx( 'usr_CO19_OH' )
    usr_CO20_OH_ndx      = get_rxt_ndx( 'usr_CO20_OH' )
    usr_CO21_OH_ndx      = get_rxt_ndx( 'usr_CO21_OH' )
    usr_CO22_OH_ndx      = get_rxt_ndx( 'usr_CO22_OH' )
    usr_CO23_OH_ndx      = get_rxt_ndx( 'usr_CO23_OH' )
    usr_CO24_OH_ndx      = get_rxt_ndx( 'usr_CO24_OH' )
    usr_CO25_OH_ndx      = get_rxt_ndx( 'usr_CO25_OH' )
    usr_CO26_OH_ndx      = get_rxt_ndx( 'usr_CO26_OH' )
    usr_CO27_OH_ndx      = get_rxt_ndx( 'usr_CO27_OH' )
    usr_CO28_OH_ndx      = get_rxt_ndx( 'usr_CO28_OH' )
    usr_CO29_OH_ndx      = get_rxt_ndx( 'usr_CO29_OH' )
    usr_CO30_OH_ndx      = get_rxt_ndx( 'usr_CO30_OH' )
    usr_CO31_OH_ndx      = get_rxt_ndx( 'usr_CO31_OH' )
    usr_CO32_OH_ndx      = get_rxt_ndx( 'usr_CO32_OH' )
    usr_CO33_OH_ndx      = get_rxt_ndx( 'usr_CO33_OH' )
    usr_CO34_OH_ndx      = get_rxt_ndx( 'usr_CO34_OH' )
    usr_CO35_OH_ndx      = get_rxt_ndx( 'usr_CO35_OH' )
    usr_CO36_OH_ndx      = get_rxt_ndx( 'usr_CO36_OH' )
    usr_CO37_OH_ndx      = get_rxt_ndx( 'usr_CO37_OH' )
    usr_CO38_OH_ndx      = get_rxt_ndx( 'usr_CO38_OH' )
    usr_CO39_OH_ndx      = get_rxt_ndx( 'usr_CO39_OH' )
    usr_CO40_OH_ndx      = get_rxt_ndx( 'usr_CO40_OH' )
    usr_CO41_OH_ndx      = get_rxt_ndx( 'usr_CO41_OH' )
    usr_CO42_OH_ndx      = get_rxt_ndx( 'usr_CO42_OH' )
!lke--

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'usrrxt_inti: diagnostics '
       write(iulog,'(10i5)') usr_O_O2_ndx,usr_HO2_HO2_ndx,tag_NO2_NO3_ndx,usr_N2O5_M_ndx,tag_NO2_OH_ndx,usr_HNO3_OH_ndx &
                            ,tag_NO2_HO2_ndx,usr_HO2NO2_M_ndx,usr_N2O5_aer_ndx,usr_NO3_aer_ndx,usr_NO2_aer_ndx &
                            ,usr_CO_OH_b_ndx,tag_C2H4_OH_ndx,tag_C3H6_OH_ndx,tag_CH3CO3_NO2_ndx,usr_PAN_M_ndx,usr_CH3COCH3_OH_ndx &
                            ,usr_MCO3_NO2_ndx,usr_MPAN_M_ndx,usr_XOOH_OH_ndx,usr_SO2_OH_ndx,usr_DMS_OH_ndx,usr_HO2_aer_ndx &
                            ,usr_GLYOXAL_aer_ndx,usr_ISOPNITA_aer_ndx,usr_ISOPNITB_aer_ndx,usr_ONITR_aer_ndx,usr_HONITR_aer_ndx &
                            ,usr_TERPNIT_aer_ndx,usr_NTERPOOH_aer_ndx,usr_NC4CHO_aer_ndx,usr_NC4CH2OH_aer_ndx

    end if

  end subroutine usrrxt_inti

  subroutine usrrxt( rxt, temp, tempi, tempe, invariants, h2ovmr,  &
                     pmid, m, sulfate, mmr, relhum, strato_sad, &
                     tropchemlev, dlat, ncol, sad_trop, reff_trop, cwat, mbar, pbuf )

!-----------------------------------------------------------------
!        ... set the user specified reaction rates
!-----------------------------------------------------------------
    
    use mo_constants,  only : pi, avo => avogadro, boltz_cgs, rgas
    use chem_mods,     only : nfs, rxntot, gas_pcnst, inv_m_ndx=>indexm
    use mo_setinv,     only : inv_o2_ndx=>o2_ndx, inv_h2o_ndx=>h2o_ndx
    use physics_buffer,only : physics_buffer_desc
    use carma_flags_mod, only : carma_hetchem_feedback
    use aero_model,      only : aero_model_surfarea
    use rad_constituents,only : rad_cnst_get_info

    implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    integer, intent(in)     :: tropchemlev(pcols)         ! trop/strat reaction separation vertical index
    real(r8), intent(in)    :: dlat(:)                    ! degrees latitude
    real(r8), intent(in)    :: temp(pcols,pver)           ! temperature (K); neutral temperature
    real(r8), intent(in)    :: tempi(pcols,pver)          ! ionic temperature (K); only used if ion chemistry
    real(r8), intent(in)    :: tempe(pcols,pver)          ! electronic temperature (K); only used if ion chemistry
    real(r8), intent(in)    :: m(ncol,pver)               ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(ncol,pver)         ! sulfate aerosol (mol/mol)
    real(r8), intent(in)    :: strato_sad(pcols,pver)     ! stratospheric aerosol sad (1/cm)
    real(r8), intent(in)    :: h2ovmr(ncol,pver)          ! water vapor (mol/mol)
    real(r8), intent(in)    :: relhum(ncol,pver)          ! relative humidity
    real(r8), intent(in)    :: pmid(pcols,pver)           ! midpoint pressure (Pa)
    real(r8), intent(in)    :: invariants(ncol,pver,nfs)  ! invariants density (/cm^3)
    real(r8), intent(in)    :: mmr(pcols,pver,gas_pcnst)  ! species concentrations (kg/kg)
    real(r8), intent(in)    :: cwat(ncol,pver) !PJC Condensed Water (liquid+ice) (kg/kg)
    real(r8), intent(in)    :: mbar(ncol,pver) !PJC Molar mass of air (g/mol)
    real(r8), intent(inout) :: rxt(ncol,pver,rxntot)      ! gas phase rates
    real(r8), intent(out)   :: sad_trop(pcols,pver)       ! tropospheric surface area density (cm2/cm3)
    real(r8), intent(out)   :: reff_trop(pcols,pver)      ! tropospheric effective radius (cm)
    type(physics_buffer_desc), pointer :: pbuf(:)

!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------
    
    real(r8), parameter :: dg = 0.1_r8            ! mole diffusion =0.1 cm2/s (Dentener, 1993)

!-----------------------------------------------------------------
! 	... reaction probabilities for heterogeneous reactions
!-----------------------------------------------------------------
    real(r8), parameter :: gamma_n2o5 = 0.10_r8         ! from Jacob, Atm Env, 34, 2131, 2000
    real(r8), parameter :: gamma_ho2  = 0.20_r8         ! 
    real(r8), parameter :: gamma_no2  = 0.0001_r8       ! 
    real(r8), parameter :: gamma_no3  = 0.001_r8        ! 
    real(r8), parameter :: gamma_glyoxal  = 2.0e-4_r8   !  Washenfelder et al, JGR, 2011
!TS1 species
    real(r8), parameter :: gamma_isopnita  = 0.005_r8        ! from Fisher et al., ACP, 2016
    real(r8), parameter :: gamma_isopnitb  = 0.005_r8        ! 
    real(r8), parameter :: gamma_onitr     = 0.005_r8        ! 
    real(r8), parameter :: gamma_honitr    = 0.005_r8        ! 
    real(r8), parameter :: gamma_terpnit   = 0.01_r8         ! 
    real(r8), parameter :: gamma_nterpooh  = 0.01_r8         ! 
    real(r8), parameter :: gamma_nc4cho    = 0.005_r8        !
    real(r8), parameter :: gamma_nc4ch2oh  = 0.005_r8        !


    integer  ::  i, k
!ADDED BY WUHU FENG, JOHN PLANE
    integer  ::  ifwh
    real(r8)     :: logt, v8kinf,v8ko,v8fc,v8k8,v8k
    real(r8)     :: v9kinf,v9ko,v9fc,v9k9,v9k
    real(r8)     :: v10kinf,v10ko,v10fc,v10k10,v10k
    real(r8)     :: v12kinf,v12ko,v12fc,v12k12,v12k
  real(r8) :: mg_vk8_0,mg_vk8_inf
  real(r8) :: mg_vk8term1,mg_vk8term2,mg_vk8term3,mg_vk8term4
  real(r8) :: mg_vk9_0,mg_vk9_inf
  real(r8) :: mg_vk9term1,mg_vk9term2,mg_vk9term3,mg_vk9term4
  real(r8) :: mg_vk10_0,mg_vk10_inf
  real(r8) :: mg_vk10term1,mg_vk10term2,mg_vk10term3,mg_vk10term4
  real(r8) :: mg_vk12_0,mg_vk12_inf
  real(r8) :: mg_vk12term1,mg_vk12term2,mg_vk12term3,mg_vk12term4
  real(r8) :: vk4_0,vk4_inf,vk4term1,vk4term2,vk4term3,vk4term4
  real(r8) :: vk5_0,vk5_inf,vk5term1,vk5term2,vk5term3,vk5term4
  real(r8) :: vk6term
  real(r8) :: vk7term
  real(r8) :: ca_vk17b_0,ca_vk17b_inf
  real(r8) :: ca_vk17bterm1,ca_vk17bterm2,ca_vk17bterm3,ca_vk17bterm4
!END ADD

    integer  ::  l
    real(r8) ::  tp(ncol)                       ! 300/t
    real(r8) ::  tinv(ncol)                     ! 1/t
    real(r8) ::  ko(ncol)   
    real(r8) ::  term1(ncol)
    real(r8) ::  term2(ncol)
    real(r8) ::  kinf(ncol)   
    real(r8) ::  fc(ncol)   
    real(r8) ::  xr(ncol)   
    real(r8) ::  sur(ncol)   
    real(r8) ::  sqrt_t(ncol)                   ! sqrt( temp )
    real(r8) ::  sqrt_t_58(ncol)                ! sqrt( temp / 58.)
    real(r8) ::  exp_fac(ncol)                  ! vector exponential
    real(r8) ::  lwc(ncol)   
    real(r8) ::  ko_m(ncol)   
    real(r8) ::  k0(ncol)   
    real(r8) ::  kinf_m(ncol)   
    real(r8) ::  o2(ncol)
    real(r8) ::  c_n2o5, c_ho2, c_no2, c_no3, c_glyoxal
!TS1 species
    real(r8) ::  c_isopnita, c_isopnitb, c_onitr, c_honitr, c_terpnit, c_nterpooh
    real(r8) ::  c_nc4cho, c_nc4ch2oh

    real(r8) ::  amas
    !-----------------------------------------------------------------
    !	... density of sulfate aerosol
    !-----------------------------------------------------------------
    real(r8), parameter :: gam1 = 0.04_r8                 ! N2O5+SUL ->2HNO3
    real(r8), parameter :: wso4 = 98._r8
    real(r8), parameter :: den  = 1.15_r8                 ! each molecule of SO4(aer) density g/cm3
    !-------------------------------------------------
    ! 	... volume of sulfate particles
    !           assuming mean rm 
    !           continient 0.05um  0.07um  0.09um
    !           ocean      0.09um  0.25um  0.37um
    !                      0.16um                  Blake JGR,7195, 1995
    !-------------------------------------------------
    real(r8), parameter :: rm1  = 0.16_r8*1.e-4_r8             ! mean radii in cm
    real(r8), parameter :: fare = 4._r8*pi*rm1*rm1             ! each mean particle(r=0.1u) area   cm2/cm3

    !-----------------------------------------------------------------------
    !        ... Aqueous phase sulfur quantities for SO2 + H2O2 and SO2 + O3
    !-----------------------------------------------------------------------
    real(r8), parameter  :: HENRY298_H2O2 =  7.45e+04_r8
    real(r8), parameter  :: H298_H2O2     = -1.45e+04_r8
    real(r8), parameter  :: HENRY298_SO2  =  1.23e+00_r8
    real(r8), parameter  :: H298_SO2      = -6.25e+03_r8
    real(r8), parameter  :: K298_SO2_HSO3 =  1.3e-02_r8
    real(r8), parameter  :: H298_SO2_HSO3 = -4.16e+03_r8
    real(r8), parameter  :: R_CONC        =  82.05e+00_r8 / avo
    real(r8), parameter  :: R_CAL         =  rgas * 0.239006e+00_r8
    real(r8), parameter  :: K_AQ          =  7.57e+07_r8
    real(r8), parameter  :: ER_AQ         =  4.43e+03_r8

    real(r8), parameter  :: HENRY298_O3   =  1.13e-02_r8
    real(r8), parameter  :: H298_O3       = -5.04e+03_r8
    real(r8), parameter  :: K298_HSO3_SO3 =  6.6e-08_r8
    real(r8), parameter  :: H298_HSO3_SO3 = -2.23e+03_r8
    real(r8), parameter  :: K0_AQ         =  2.4e+04_r8
    real(r8), parameter  :: ER0_AQ        =  0.0e+00_r8
    real(r8), parameter  :: K1_AQ         =  3.7e+05_r8
    real(r8), parameter  :: ER1_AQ        =  5.53e+03_r8
    real(r8), parameter  :: K2_AQ         =  1.5e+09_r8
    real(r8), parameter  :: ER2_AQ        =  5.28e+03_r8

    real(r8), parameter  :: pH            =  4.5e+00_r8
    
    real(r8), pointer :: sfc(:), dm_aer(:)
    integer :: ntot_amode

    real(r8), pointer :: sfc_array(:,:,:), dm_array(:,:,:)
    
    ! get info about the modal aerosols
    ! get ntot_amode
    call rad_cnst_get_info(0, nmodes=ntot_amode)

    if (ntot_amode>0) then
       allocate(sfc_array(pcols,pver,ntot_amode), dm_array(pcols,pver,ntot_amode) )
    else
       allocate(sfc_array(pcols,pver,5), dm_array(pcols,pver,5) )
    endif

    sfc_array(:,:,:) = 0._r8
    dm_array(:,:,:) = 0._r8
    sad_trop(:,:) = 0._r8
    reff_trop(:,:) = 0._r8
        
    if( usr_NO2_aer_ndx > 0 .or. usr_NO3_aer_ndx > 0 .or. usr_N2O5_aer_ndx > 0 .or. usr_HO2_aer_ndx > 0 ) then

! sad_trop should be set outside of usrrxt ?? 
       if( carma_hetchem_feedback ) then
          sad_trop(:ncol,:pver)=strato_sad(:ncol,:pver)
       else

          call aero_model_surfarea( &
               mmr, rm1, relhum, pmid, temp, strato_sad, sulfate, m, tropchemlev, dlat, &
               het1_ndx, pbuf, ncol, sfc_array, dm_array, sad_trop, reff_trop )

       endif
    endif

    level_loop : do k = 1,pver
       tinv(:)           = 1._r8 / temp(:ncol,k)
       tp(:)             = 300._r8 * tinv(:)
       sqrt_t(:)         = sqrt( temp(:ncol,k) )
       sqrt_t_58(:)      = sqrt( temp(:ncol,k) / 58.0_r8 )

!-----------------------------------------------------------------
!	... o + o2 + m --> o3 + m (JPL15-10)
!-----------------------------------------------------------------
       if( usr_O_O2_ndx > 0 ) then
          rxt(:,k,usr_O_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8
       end if
       if( usr_OA_O2_ndx > 0 ) then
          rxt(:,k,usr_OA_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8
       end if

!-----------------------------------------------------------------
!	... o + o + m -> o2 + m
!-----------------------------------------------------------------
       if ( usr_O_O_ndx > 0 ) then
          rxt(:,k,usr_O_O_ndx) = 2.76e-34_r8 * exp( 720.0_r8*tinv(:) )
       end if
         
!-----------------------------------------------------------------
! 	... cl2o2 + m -> 2*clo + m  (JPL15-10)
!-----------------------------------------------------------------
       if ( usr_CL2O2_M_ndx > 0 ) then
          if ( tag_CLO_CLO_M_ndx > 0 ) then
             ko(:)            = 2.16e-27_r8 * exp( 8537.0_r8* tinv(:) )
             rxt(:,k,usr_CL2O2_M_ndx) = rxt(:,k,tag_CLO_CLO_M_ndx)/ko(:)         
          else
             rxt(:,k,usr_CL2O2_M_ndx) = 0._r8
          end if
       end if
       
!-----------------------------------------------------------------
!       ... so3 + 2*h2o --> h2so4 + h2o
!       Note: this reaction proceeds by the 2 intermediate steps below
!           so3 + h2o --> adduct
!           adduct + h2o --> h2so4 + h2o  
!               (Lovejoy et al., JCP, pp. 19911-19916, 1996)
!	The first order rate constant used here is recommended by JPL 2011.
!	This rate involves the water vapor number density.
!-----------------------------------------------------------------

       if ( usr_SO3_H2O_ndx > 0 ) then
          call comp_exp( exp_fac, 6540.0_r8*tinv(:), ncol )
          if( h2o_ndx > 0 ) then
             fc(:) = 8.5e-21_r8 * m(:,k) * h2ovmr(:,k) * exp_fac(:)
          else
             fc(:) = 8.5e-21_r8 * invariants(:,k,inv_h2o_ndx) * exp_fac(:)
          end if
          rxt(:,k,usr_SO3_H2O_ndx) = 1.0e-20_r8 * fc(:)
       end if
         
!-----------------------------------------------------------------
!	... n2o5 + m --> no2 + no3 + m (JPL15-10)
!-----------------------------------------------------------------
       if( usr_N2O5_M_ndx > 0 ) then
          if( tag_NO2_NO3_ndx > 0 ) then
             call comp_exp( exp_fac, -10840.0_r8*tinv, ncol )
             rxt(:,k,usr_N2O5_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) * 1.724138e26_r8 * exp_fac(:)
          else
             rxt(:,k,usr_N2O5_M_ndx) = 0._r8
          end if
       end if
       if( usr_XNO2NO3_M_ndx > 0 ) then
          if( tag_NO2_NO3_ndx > 0 ) then
             call comp_exp( exp_fac, -10840.0_r8*tinv, ncol )
             rxt(:,k,usr_XNO2NO3_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) *1.724138e26_r8 * exp_fac(:)
          else
             rxt(:,k,usr_XNO2NO3_M_ndx) = 0._r8
          end if
       end if
       if( usr_NO2XNO3_M_ndx > 0 ) then
          if( tag_NO2_NO3_ndx > 0 ) then
             call comp_exp( exp_fac, -10840.0_r8*tinv, ncol )
             rxt(:,k,usr_NO2XNO3_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) * 1.734138e26_r8 * exp_fac(:)
          else
             rxt(:,k,usr_NO2XNO3_M_ndx) = 0._r8
          end if
       end if

!-----------------------------------------------------------------
!	set rates for:
! 	... hno3 + oh --> no3 + h2o
!           ho2no2 + m --> ho2 + no2 + m
!-----------------------------------------------------------------
       if( usr_HNO3_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 1335._r8*tinv, ncol )
          ko(:) = m(:,k) * 6.5e-34_r8 * exp_fac(:)
          call comp_exp( exp_fac, 2199._r8*tinv, ncol )
          ko(:) = ko(:) / (1._r8 + ko(:)/(2.7e-17_r8*exp_fac(:)))
          call comp_exp( exp_fac, 460._r8*tinv, ncol )
          rxt(:,k,usr_HNO3_OH_ndx) = ko(:) + 2.4e-14_r8*exp_fac(:)
       end if
       if( usr_XHNO3_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 1335._r8*tinv, ncol )
          ko(:) = m(:,k) * 6.5e-34_r8 * exp_fac(:)
          call comp_exp( exp_fac, 2199._r8*tinv, ncol )
          ko(:) = ko(:) / (1._r8 + ko(:)/(2.7e-17_r8*exp_fac(:)))
          call comp_exp( exp_fac, 460._r8*tinv, ncol )
          rxt(:,k,usr_XHNO3_OH_ndx) = ko(:) + 2.4e-14_r8*exp_fac(:)
       end if
       if( usr_HO2NO2_M_ndx > 0 ) then
          if( tag_NO2_HO2_ndx > 0 ) then
             call comp_exp( exp_fac, -10900._r8*tinv, ncol )
             rxt(:,k,usr_HO2NO2_M_ndx) = rxt(:,k,tag_NO2_HO2_ndx) * exp_fac(:) / 2.1e-27_r8
          else
             rxt(:,k,usr_HO2NO2_M_ndx) = 0._r8
          end if
       end if
       if( usr_XHO2NO2_M_ndx > 0 ) then
          if( tag_NO2_HO2_ndx > 0 ) then
             call comp_exp( exp_fac, -10900._r8*tinv, ncol )
             rxt(:,k,usr_XHO2NO2_M_ndx) = rxt(:,k,tag_NO2_HO2_ndx) * exp_fac(:) / 2.1e-27_r8
          else
             rxt(:,k,usr_XHO2NO2_M_ndx) = 0._r8
          end if
       end if
!-----------------------------------------------------------------
!           co + oh --> co2 + ho2     (combined branches - do not use with CO_OH_b)
!-----------------------------------------------------------------
       if( usr_CO_OH_a_ndx > 0 ) then
          rxt(:,k,usr_CO_OH_a_ndx) = 1.5e-13_r8 * &
               (1._r8 + 6.e-7_r8*boltz_cgs*m(:,k)*temp(:ncol,k))
       end if
!-----------------------------------------------------------------
! 	... co + oh --> co2 + h (second branch JPL15-10, with CO+OH+M)
!-----------------------------------------------------------------
       if( usr_CO_OH_b_ndx > 0 ) then
         kinf(:)  = 2.1e+09_r8 * (temp(:ncol,k)/ t0)**(6.1_r8)
         ko  (:)  = 1.5e-13_r8

         term1(:) = ko(:) / ( (kinf(:) / m(:,k)) )
         term2(:) = ko(:) / (1._r8 + term1(:))

         term1(:) = log10( term1(:) )
         term1(:) = 1.0_r8 / (1.0_r8 + term1(:)*term1(:))

         rxt(:ncol,k,usr_CO_OH_b_ndx) = term2(:) * (0.6_r8)**term1(:)
       end if

!-----------------------------------------------------------------
!       ... ho2 + ho2 --> h2o2
!       note: this rate involves the water vapor number density
!-----------------------------------------------------------------
       if( usr_HO2_HO2_ndx > 0 ) then

          call comp_exp( exp_fac, 460._r8*tinv, ncol )
          ko(:)   = 3.0e-13_r8 * exp_fac(:)
          call comp_exp( exp_fac, 920._r8*tinv, ncol )
          kinf(:) = 2.1e-33_r8 * m(:,k) * exp_fac(:)
          call comp_exp( exp_fac, 2200._r8*tinv, ncol )

          if( h2o_ndx > 0 ) then
             fc(:) = 1._r8 + 1.4e-21_r8 * m(:,k) * h2ovmr(:,k) * exp_fac(:)
          else
             fc(:) = 1._r8 + 1.4e-21_r8 * invariants(:,k,inv_h2o_ndx) * exp_fac(:)
          end if
          rxt(:,k,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)

       end if

!-----------------------------------------------------------------
!    	... mco3 + no2 -> mpan
!-----------------------------------------------------------------
       if( usr_MCO3_NO2_ndx > 0 ) then
          rxt(:,k,usr_MCO3_NO2_ndx) = 1.1e-11_r8 * tp(:) / m(:,k)
       end if
       if( usr_MCO3_XNO2_ndx > 0 ) then
          rxt(:,k,usr_MCO3_XNO2_ndx) = 1.1e-11_r8 * tp(:) / m(:,k)
       end if

!-----------------------------------------------------------------
!	... pan + m --> ch3co3 + no2 + m (JPL15-10)
!-----------------------------------------------------------------
       call comp_exp( exp_fac, -14000._r8*tinv, ncol )
       if( usr_PAN_M_ndx > 0 ) then
          if( tag_CH3CO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_PAN_M_ndx) = rxt(:,k,tag_CH3CO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_PAN_M_ndx) = 0._r8
          end if
       end if
       if( usr_XPAN_M_ndx > 0 ) then
          if( tag_CH3CO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_XPAN_M_ndx) = rxt(:,k,tag_CH3CO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_XPAN_M_ndx) = 0._r8
          end if
       end if

!-----------------------------------------------------------------
!	... mpan + m --> mco3 + no2 + m (JPL15-10)
!-----------------------------------------------------------------
       if( usr_MPAN_M_ndx > 0 ) then
          if( usr_MCO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_MPAN_M_ndx) = rxt(:,k,usr_MCO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_MPAN_M_ndx) = 0._r8
          end if
       end if
       if( usr_XMPAN_M_ndx > 0 ) then
          if( usr_MCO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_XMPAN_M_ndx) = rxt(:,k,usr_MCO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_XMPAN_M_ndx) = 0._r8
          end if
       end if

!lke-TS1
!-----------------------------------------------------------------
!       ... pbznit + m --> acbzo2 + no2 + m
!-----------------------------------------------------------------
       if( usr_PBZNIT_M_ndx > 0 ) then
          if( tag_ACBZO2_NO2_ndx > 0 ) then
             rxt(:,k,usr_PBZNIT_M_ndx) = rxt(:,k,tag_ACBZO2_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_PBZNIT_M_ndx) = 0._r8
          end if
       end if

!-----------------------------------------------------------------
!       ... xooh + oh -> h2o + oh
!-----------------------------------------------------------------
       if( usr_XOOH_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 253._r8*tinv, ncol )
          rxt(:,k,usr_XOOH_OH_ndx) = temp(:ncol,k)**2._r8 * 7.69e-17_r8 * exp_fac(:)
       end if

!-----------------------------------------------------------------
!       ... ch3coch3 + oh -> ro2 + h2o
!-----------------------------------------------------------------
       if( usr_CH3COCH3_OH_ndx > 0 ) then
          call comp_exp( exp_fac, -2000._r8*tinv, ncol )
          rxt(:,k,usr_CH3COCH3_OH_ndx) = 3.82e-11_r8 * exp_fac(:) + 1.33e-13_r8
       end if

!-----------------------------------------------------------------
!       ... DMS + OH  --> .5 * SO2
!-----------------------------------------------------------------
       if( usr_DMS_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 7460._r8*tinv, ncol )
          ko(:) = 1._r8 + 5.5e-31_r8 * exp_fac * m(:,k) * 0.21_r8
          call comp_exp( exp_fac, 7810._r8*tinv, ncol )
          rxt(:,k,usr_DMS_OH_ndx) = 1.7e-42_r8 * exp_fac * m(:,k) * 0.21_r8 / ko(:)
       end if

!-----------------------------------------------------------------
!       ... SO2 + OH  --> SO4  (REFERENCE?? - not Liao)
!-----------------------------------------------------------------
       if( usr_SO2_OH_ndx > 0 ) then
          fc(:) = 3.0e-31_r8 *(300._r8*tinv(:))**3.3_r8
          ko(:) = fc(:)*m(:,k)/(1._r8 + fc(:)*m(:,k)/1.5e-12_r8) 
          rxt(:,k,usr_SO2_OH_ndx) = ko(:)*.6_r8**(1._r8 + (log10(fc(:)*m(:,k)/1.5e-12_r8))**2._r8)**(-1._r8)
       end if
!
! reduced hydrocarbon scheme
!
       if ( usr_C2O3_NO2_ndx > 0 ) then
          ko(:)   = 2.6e-28_r8 * m(:,k)
          kinf(:) = 1.2e-11_r8
          rxt(:,k,usr_C2O3_NO2_ndx) = (ko/(1._r8+ko/kinf)) * 0.6_r8**(1._r8/(1._r8+(log10(ko/kinf))**2))
       end if
       if ( usr_C2O3_XNO2_ndx > 0 ) then
          ko(:)   = 2.6e-28_r8 * m(:,k)
          kinf(:) = 1.2e-11_r8
          rxt(:,k,usr_C2O3_XNO2_ndx) = (ko/(1._r8+ko/kinf)) * 0.6_r8**(1._r8/(1._r8+(log10(ko/kinf))**2))
       end if
       if ( usr_C2H4_OH_ndx > 0 ) then
          ko(:)   = 1.0e-28_r8 * m(:,k)
          kinf(:) = 8.8e-12_r8
          rxt(:,k,usr_C2H4_OH_ndx) = (ko/(1._r8+ko/kinf)) * 0.6_r8**(1._r8/(1._r8+(log(ko/kinf))**2))
       end if
       if ( usr_XO2N_HO2_ndx > 0 ) then
          rxt(:,k,usr_XO2N_HO2_ndx) = rxt(:,k,tag_XO2N_NO_ndx)*rxt(:,k,tag_XO2_HO2_ndx)/(rxt(:,k,tag_XO2_NO_ndx)+1.e-36_r8)
       end if
       
!
! hydrolysis reactions on wetted aerosols
!      
       if( usr_NO2_aer_ndx > 0 .or. usr_NO3_aer_ndx > 0 .or. usr_N2O5_aer_ndx > 0 .or. usr_HO2_aer_ndx > 0 &
         .or. usr_GLYOXAL_aer_ndx > 0 ) then

          long_loop : do i = 1,ncol

             sfc    => sfc_array(i,k,:)
             dm_aer => dm_array(i,k,:)

             c_n2o5 = 1.40e3_r8 * sqrt_t(i)         ! mean molecular speed of n2o5
             c_no3  = 1.85e3_r8 * sqrt_t(i)         ! mean molecular speed of no3
             c_no2  = 2.15e3_r8 * sqrt_t(i)         ! mean molecular speed of no2
             c_ho2  = 2.53e3_r8 * sqrt_t(i)         ! mean molecular speed of ho2
             c_glyoxal = 1.455e4_r8 * sqrt_t_58(i)  ! mean molecular speed of ho2
             c_isopnita = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of isopnita
             c_isopnitb = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of isopnitb
             c_onitr    = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of onitr 
             c_honitr   = 1.26e3_r8 * sqrt_t(i)         ! mean molecular speed of honitr
             c_terpnit  = 0.992e3_r8 * sqrt_t(i)        ! mean molecular speed of terpnit
             c_nterpooh = 0.957e3_r8 * sqrt_t(i)        ! mean molecular speed of nterpooh
             c_nc4cho   = 1.21e3_r8 * sqrt_t(i)         ! mean molecular speed of nc4cho
             c_nc4ch2oh = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of nc4ch2oh

             !-------------------------------------------------------------------------
             !  Heterogeneous reaction rates for uptake of a gas on an aerosol:
             !    rxt = sfc / ( (rad_aer/Dg_gas) + (4/(c_gas*gamma_gas)))
             !-------------------------------------------------------------------------
             !-------------------------------------------------------------------------
             ! 	... n2o5 -> 2 hno3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_N2O5_aer_ndx > 0 ) then
                rxt(i,k,usr_N2O5_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
             end if
             if( usr_XNO2NO3_aer_ndx > 0 ) then
                rxt(i,k,usr_XNO2NO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
             end if
             if( usr_NO2XNO3_aer_ndx > 0 ) then
                rxt(i,k,usr_NO2XNO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
             end if
             !-------------------------------------------------------------------------
             ! 	... no3 -> hno3  (on sulfate, nh4no3, oc, soa)
             !-------------------------------------------------------------------------
             if( usr_NO3_aer_ndx > 0 ) then
                rxt(i,k,usr_NO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no3, gamma_no3 )
             end if
             if( usr_XNO3_aer_ndx > 0 ) then
                rxt(i,k,usr_XNO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no3, gamma_no3 )
             end if
             !-------------------------------------------------------------------------
             ! 	... no2 -> 0.5 * (ho+no+hno3)  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_NO2_aer_ndx > 0 ) then
                rxt(i,k,usr_NO2_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no2, gamma_no2 )
             end if
             if( usr_XNO2_aer_ndx > 0 ) then
                rxt(i,k,usr_XNO2_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no2, gamma_no2 )
             end if

             !-------------------------------------------------------------------------
             ! 	... ho2 -> 0.5 * h2o2  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_HO2_aer_ndx > 0 ) then
                rxt(i,k,usr_HO2_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_ho2, gamma_ho2 )
             end if
             !-------------------------------------------------------------------------
             !  ... glyoxal ->  soag1  (on sulfate, nh4no3, oc2, soa)  
             ! first order uptake, Fuchs and Sutugin, 1971,  dCg = 1/4 * gamma * ! A * |v_mol| * Cg * dt
             !-------------------------------------------------------------------------
             if( usr_GLYOXAL_aer_ndx > 0 ) then
                rxt(i,k,usr_GLYOXAL_aer_ndx) = hetrxtrate_gly( sfc, c_glyoxal, gamma_glyoxal )
             end if
             !-------------------------------------------------------------------------
             ! 	... ISOPNITA -> HNO3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPNITA_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPNITA_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopnita, gamma_isopnita )
             end if
             !-------------------------------------------------------------------------
             ! 	... ISOPNITB -> HNO3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPNITB_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPNITB_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopnitb, gamma_isopnitb )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ONITR -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ONITR_aer_ndx > 0 ) then
                rxt(i,k,usr_ONITR_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_onitr, gamma_onitr )
             end if
             !-------------------------------------------------------------------------
             ! 	... HONITR -> HNO3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_HONITR_aer_ndx > 0 ) then
                rxt(i,k,usr_HONITR_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_honitr, gamma_honitr )
             end if
             !-------------------------------------------------------------------------
             ! 	... TERPNIT -> HNO3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPNIT_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPNIT_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terpnit, gamma_terpnit )
             end if
             !-------------------------------------------------------------------------
             ! 	...  NTERPOOH -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_NTERPOOH_aer_ndx > 0 ) then
                rxt(i,k,usr_NTERPOOH_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_nterpooh, gamma_nterpooh )
             end if
             !-------------------------------------------------------------------------
             ! 	...  NC4CHO -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_NC4CHO_aer_ndx > 0 ) then
                rxt(i,k,usr_NC4CHO_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_nc4cho, gamma_nc4cho )
             end if
             !-------------------------------------------------------------------------
             ! 	...  NC4CH2OH -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_NC4CH2OH_aer_ndx > 0 ) then
                rxt(i,k,usr_NC4CH2OH_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_nc4ch2oh, gamma_nc4ch2oh )
             end if

          end do long_loop
       end if

       ! LLNL super fast chem reaction rates

       !-----------------------------------------------------------------------
       !     ... CO + OH --> CO2 + HO2
       !-----------------------------------------------------------------------
       if ( usr_oh_co_ndx > 0 ) then
          ko(:)     = 5.9e-33_r8 * tp(:)**1.4_r8
          kinf(:)   = 1.1e-12_r8 * (temp(:ncol,k) / 300._r8)**1.3_r8
          ko_m(:)   = ko(:) * m(:,k)
          k0(:)     = 1.5e-13_r8 * (temp(:ncol,k) / 300._r8)**0.6_r8
          kinf_m(:) = (2.1e+09_r8 * (temp(:ncol,k) / 300._r8)**6.1_r8) / m(:,k)
          rxt(:,k,usr_oh_co_ndx) = (ko_m(:)/(1._r8+(ko_m(:)/kinf(:)))) * &
               0.6_r8**(1._r8/(1._r8+(log10(ko_m(:)/kinf(:)))**2._r8)) + &
               (k0(:)/(1._r8+(k0(:)/kinf_m(:)))) * &
               0.6_r8**(1._r8/(1._r8+(log10(k0(:)/kinf_m(:)))**2._r8))
       endif
       !-----------------------------------------------------------------------
       !     ... NO2 + H2O --> 0.5 HONO + 0.5 HNO3
       !-----------------------------------------------------------------------
       if ( het_no2_h2o_ndx > 0 ) then
          rxt(:,k,het_no2_h2o_ndx) = 4.0e-24_r8
       endif
       !-----------------------------------------------------------------------
       !     ... DMS + OH --> 0.75 SO2 + 0.25 MSA
       !-----------------------------------------------------------------------
       if ( usr_oh_dms_ndx > 0 ) then
          o2(:ncol) = invariants(:ncol,k,inv_o2_ndx)
          rxt(:,k,usr_oh_dms_ndx) = 2.000e-10_r8 * exp(5820.0_r8 * tinv(:)) / &
               ((2.000e29_r8 / o2(:)) + exp(6280.0_r8 * tinv(:)))
       endif
       if ( aq_so2_h2o2_ndx > 0 .or. aq_so2_o3_ndx > 0 ) then
          lwc(:) = cwat(:ncol,k) * invariants(:ncol,k,inv_m_ndx) * mbar(:ncol,k) /avo !PJC convert kg/kg to g/cm3
          !-----------------------------------------------------------------------
          !     ... SO2 + H2O2 --> S(VI)
          !-----------------------------------------------------------------------
          if ( aq_so2_h2o2_ndx > 0 ) then
             rxt(:,k,aq_so2_h2o2_ndx) = lwc(:) * 1.0e-03_r8 * avo * &
                  K_AQ * &

                  exp(ER_AQ * ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  HENRY298_SO2 * &
                  K298_SO2_HSO3 * &
                  HENRY298_H2O2 * &
                  exp(((H298_SO2 + H298_SO2_HSO3 + H298_H2O2) / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (R_CONC * temp(:ncol,k))**2.0e+00_r8 / &

                  (1.0e+00_r8 + 13.0e+00_r8 * 10.0e+00_r8**(-pH))
          endif
          !-----------------------------------------------------------------------
          !     ... SO2 + O3 --> S(VI)
          !-----------------------------------------------------------------------
          if (aq_so2_o3_ndx >0) then
             rxt(:,k,aq_so2_o3_ndx) = lwc(:) * 1.0e-03_r8 * avo * &
                  HENRY298_SO2 * exp((H298_SO2 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (K0_AQ * exp(ER0_AQ * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) + &
                  K298_SO2_HSO3 * exp((H298_SO2_HSO3 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (K1_AQ * exp(ER1_AQ * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) / &
                  10.0e+00_r8**(-pH) + K2_AQ * exp(ER2_AQ * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  K298_HSO3_SO3 * exp((H298_HSO3_SO3 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) / &
                  (10.0e+00_r8**(-pH))**2.0e+00_r8) ) * &
                  HENRY298_O3 * exp((H298_O3 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (R_CONC * temp(:ncol,k))**2.0e+00_r8
          endif
       endif

    if ( has_d_chem ) then

        call comp_exp( exp_fac, -600._r8 * tinv, ncol )
        rxt(:,k,ean_ndx(1))  = 1.e-31_r8 * tp(:) * exp_fac(:)
        rxt(:,k,ean_ndx(2))  = 9.1e-12_r8 * tp(:)**(-1.46_r8)
        call comp_exp( exp_fac, -193._r8 * tinv, ncol )
        rxt(:,k,ean_ndx(3))  = (4.e-30_r8 * exp_fac(:)) * 0.21_r8

        rxt(:,k,rpe_ndx(1))  = 4.2e-6_r8 * tp(:)**0.5_r8
        rxt(:,k,rpe_ndx(2))  = 6.3e-7_r8 * tp(:)**0.5_r8
        rxt(:,k,rpe_ndx(3))  = 2.5e-6_r8 * tp(:)**0.1_r8
        rxt(:,k,rpe_ndx(4))  = 2.48e-6_r8 * tp(:)**0.76_r8
        rxt(:,k,rpe_ndx(5))  = 1.4e-6_r8 * tp(:)**0.4_r8

        rxt(:,k,pir_ndx(1)) = 4.e-30_r8 * tp(:)**2.93_r8
        rxt(:,k,pir_ndx(2))  = 4.6e-27_r8 * tp(:)**4._r8

        call comp_exp( exp_fac, -15900._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(3))  = (2.5e-2_r8 * tp(:)**5._r8) * exp_fac(:)
        rxt(:,k,pir_ndx(4))  = 2.3e-27_r8 * tp(:)**7.5_r8

        call comp_exp( exp_fac, -10272._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(5))  = (2.6e-3_r8 * tp(:)**8.5_r8) * exp_fac(:)
        rxt(:,k,pir_ndx(6))  = 3.6e-27_r8 * tp(:)**8.1_r8

        call comp_exp( exp_fac, -9000._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(7))  = (1.5e-1_r8 * tp(:)**9.1_r8) * exp_fac(:)
        rxt(:,k,pir_ndx(8))  = 4.6e-28_r8 * tp(:)**14._r8

        call comp_exp( exp_fac, -6400._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(9))  = (1.7e-3_r8 * tp(:)**15._r8) * exp_fac(:)
        rxt(:,k,pir_ndx(10)) = 1.35e-28_r8 * tp(:)**2.83_r8

        rxt(:,k,pir_ndx(11)) = 1.e-27_r8 * (308._r8 * tinv(:))**4.7_r8
        rxt(:,k,pir_ndx(12)) = rxt(:,k,pir_ndx(11))
        rxt(:,k,pir_ndx(13)) = 1.4e-29_r8 * tp(:)**4._r8

        call comp_exp( exp_fac, -3872._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(14)) = (3.4e-7_r8 * tp(:)**5._r8) * exp_fac(:)

        rxt(:,k,pir_ndx(15)) = 3.0e-31_r8 * tp(:)**4.3_r8
        call comp_exp( exp_fac, -2093._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(16)) = (1.5e-8_r8 * tp(:)**4.3_r8) * exp_fac(:)

        rxt(:,k,edn_ndx(1)) = 3.1e-10_r8 * tp(:)**0.83_r8
        call comp_exp( exp_fac, -4990._r8 * tinv, ncol )
        rxt(:,k,edn_ndx(2)) = (1.9e-12_r8 * tp(:)**(-1.5_r8)) * exp_fac(:)

        rxt(:,k,nir_ndx(1)) = 1.05e-12_r8 * tp(:)**2.15_r8
        rxt(:,k,nir_ndx(2)) = 2.5e-11_r8 * tp(:)**0.79_r8
        rxt(:,k,nir_ndx(3)) = 7.5e-11_r8 * tp(:)**0.79_r8
        rxt(:,k,nir_ndx(4)) = rxt(:,k,nir_ndx(1))
        rxt(:,k,nir_ndx(5)) = 1.3e-11_r8 * tp(:)**1.64_r8
        rxt(:,k,nir_ndx(6)) = 3.3e-11_r8 * tp(:)**2.38_r8

        call comp_exp( exp_fac, -7300_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(7)) = (1.0e-3_r8 * tp(:)) * exp_fac(:)
        call comp_exp( exp_fac, -7050_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(8)) = (7.2e-4_r8 * tp(:)) * exp_fac(:)
        call comp_exp( exp_fac, -6800_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(9)) = (6.5e-3_r8 * tp(:)) * exp_fac(:)
        call comp_exp( exp_fac, -7600_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(10)) = (5.7e-4_r8 * tp(:)) * exp_fac(:)

        call comp_exp( exp_fac, -7150_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(11)) = (1.5e-2_r8 * tp(:)) * exp_fac(:)

        call comp_exp( exp_fac, -13130_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(12)) = (6.0e-3_r8 * tp(:)) * exp_fac(:)
        rxt(:,k,nir_ndx(13)) = 5.22e-28_r8 * tp(:)**2.62_r8

        rxt(:,k,iira_ndx(1)) = 6.0e-8_r8 * tp(:)**.5_r8
        do i = 2,niira
          rxt(:,k,iira_ndx(i)) = rxt(:,k,iira_ndx(i-1))
        enddo

        rxt(:,k,iirb_ndx(1)) = 1.25e-25_r8 * tp(:)**4._r8
        do i = 2,niirb
          rxt(:,k,iirb_ndx(i)) = rxt(:,k,iirb_ndx(i-1))
        enddo

        call comp_exp( exp_fac, -6600._r8 * tinv, ncol )
        rxt(:,k,usr_clm_h2o_m_ndx) = 2.e-8_r8 * exp_fac(:)

        call comp_exp( exp_fac, -11926._r8 * tinv, ncol )
        rxt(:,k,usr_clm_hcl_m_ndx) =  tinv(:) * exp_fac(:)

     endif
    end do level_loop
    
!-----------------------------------------------------------------
! 	... the ionic rates
!-----------------------------------------------------------------
    if ( has_ion_rxts ) then
       level_loop2 : do k = 1,pver
 	   tp(:ncol)         = (2._r8*tempi(:ncol,k) + temp(:ncol,k)) / ( 3._r8 * t0 )
	   tp(:)             = max( min( tp(:),20._r8 ),1._r8 )
	   rxt(:,k,ion1_ndx) = 2.82e-11_r8 + tp(:)*(-7.74e-12_r8 + tp(:)*(1.073e-12_r8  &
			 + tp(:)*(-5.17e-14_r8 + 9.65e-16_r8*tp(:))))
	   tp(:ncol)         = (.6363_r8*tempi(:ncol,k) + .3637_r8*temp(:ncol,k)) / t0
	   tp(:)             = max( min( tp(:),trlim2 ),1._r8 )
	   rxt(:,k,ion2_ndx) = 1.533e-12_r8 + tp(:)*(-5.92e-13_r8 + tp(:)*8.6e-14_r8)
	   tp(:ncol)         = 2._r8 * t0 /(tempi(:ncol,k) + temp(:ncol,k))
	   where( tp(:ncol) < trlim3 )
		  rxt(:,k,ion3_ndx)  = 1.4e-10_r8 * tp(:)**.44_r8
		  rxt(:,k,ion11_ndx) = 1.e-11_r8 * tp(:)**.23_r8
           elsewhere
		  rxt(:,k,ion3_ndx)  = 5.2e-11_r8 / tp(:)**.2_r8
	      rxt(:,k,ion11_ndx) = 3.6e-12_r8 / tp(:)**.41_r8
	   end where
	   tp(:ncol)          = t0 / tempe(:ncol,k)
	   rxt(:,k,elec1_ndx) = 4.e-7_r8 * tp(:)**.85_r8
	   rxt(:,k,elec3_ndx) = 1.8e-7_r8 * tp(:)**.39_r8
	   where( tp(:ncol) < 4._r8 )
	      rxt(:,k,elec2_ndx) = 2.7e-7_r8 * tp(:)**.7_r8
	   elsewhere
	      rxt(:,k,elec2_ndx) = 1.6e-7_r8 * tp(:)**.55_r8
	   end where
	end do level_loop2
     endif

     ! quenching of O+(2P) and O+(2D) by e to produce O+
     ! See TABLE 1 of Roble (1995) 
     ! drm 2015-07-27
     if (elec4_ndx > 0 .and. elec5_ndx > 0 .and. elec6_ndx > 0) then
         do k=1,pver
            tp(:ncol)          = sqrt(300._r8 / tempe(:ncol,k))
            rxt(:,k,elec4_ndx) = 1.5e-7_r8 * tp(:)
            rxt(:,k,elec5_ndx) = 4.0e-8_r8 * tp(:)
            rxt(:,k,elec6_ndx) = 6.6e-8_r8 * tp(:)
         end do 
     endif
!WUHU FENG, JOHN PLANE,DAN MARSH, 18/09/2018
      do k=1,pver
      do i=1,ncol
!-----------------------------------------------------------------
!ADDED BY WUHU FENG AND DAN MARSH, John Plane, Martyn Chipperfield
!-----------------------------------------------------------------
!!-----------------------------------------------------------------
      logt=log10(temp(i,k))

!!-----------------------------------------------------------------
!!       ... Na chemistry reactions
!!-----------------------------------------------------------------
!!       ... NaO + O -> Na + O2
!!-----------------------------------------------------------------
       rxt(i,k,na02_ndx) = 2.2e-10_r8 * (temp(i,k)/ 200.0_r8)**(0.5_r8)
!!-----------------------------------------------------------------
!!       ... NaOH + CO2 + M -> NaHCO3
!!-----------------------------------------------------------------
       rxt(i,k,na08_ndx) = 1.9e-28_r8 * 200.0_r8 /temp(i,k)
!!-----------------------------------------------------------------
!!       ...  NaHCO3 + H ->  Na + H2O + CO2
!!-----------------------------------------------------------------
       rxt(i,k,na09_ndx) = 1.84e-13_r8 * exp (-1014.0_r8/temp(i,k)) &
                            * (temp(i,k))**(0.777_r8)
!!-----------------------------------------------------------------
!!       ... Na + O2 + M -> NaO2
!!-----------------------------------------------------------------
       rxt(i,k,na10_ndx) = 5.0e-30_r8 * (200.0_r8 /temp(i,k))**(1.22_r8)
!!-----------------------------------------------------------------
!!       ... 2NaHCO3 -> (NaHCO3)2
!!-----------------------------------------------------------------
!       rxt(i,k,na12_ndx) = 8.8e-10_r8 * (200.0_r8 /temp(i,k))**(0.23_r8)
        rxt(i,k,na12_ndx) = 3._r8*8.8e-10_r8 * (200.0_r8 /temp(i,k))**(0.23_r8)
!!-----------------------------------------------------------------
!!       ... Nap + N2 + M -> NaN2p
!!-----------------------------------------------------------------
       rxt(i,k,na22_ndx) = 4.8e-30_r8 * (200.0_r8 /temp(i,k))**(2.2_r8)
!!-----------------------------------------------------------------
!!       ... Nap + CO2 + M -> NaCO2p
!!-----------------------------------------------------------------
       rxt(i,k,na23_ndx) = 3.7e-29_r8 * (200.0_r8 /temp(i,k))**(2.9_r8)
!!-----------------------------------------------------------------
!!-----------------------------------------------------------------
!!       ... NaN2p + e -> Na + N2
!!-----------------------------------------------------------------
       rxt(i,k,na29a_ndx) = 1.0e-6_r8 * (200.0_r8 /temp(i,k))**(0.5_r8)
       rxt(i,k,na29b_ndx) = rxt(i,k,na29a_ndx)
       rxt(i,k,na29c_ndx) = rxt(i,k,na29a_ndx)
       rxt(i,k,na29d_ndx) = rxt(i,k,na29a_ndx)
!!-----------------------------------------------------------------
!!       ... Nap + e -> Na +hv
!!-----------------------------------------------------------------
       rxt(i,k,na30_ndx) = 3.9e-12_r8 * (200.0_r8 /temp(i,k))**(0.74_r8)
!!-----------------------------------------------------------------


!-----------------------------------------------------------------
!ADDED BY WUHU FENG AND DAN MARSH, John Plane, Martyn Chipperfield
!-----------------------------------------------------------------
!       ... Fe chemistry reactions
!-----------------------------------------------------------------
!       ... FeO + O2 + M -> FeO3 + M
!-----------------------------------------------------------------
      rxt(i,k,fe04_ndx)=4.41e-30_r8*(temp(i,k)/200.0_r8)**(0.606_r8)
!-----------------------------------------------------------------
!       ... FeO + H2O + M -> FeO2H2+M
!-----------------------------------------------------------------
      rxt(i,k,fe09_ndx)=5.0e-28_r8*(200.0_r8/temp(i,k))**(1.13_r8)
!-----------------------------------------------------------------
!       ... Fep+ N2 + M -> FeN2p
!-----------------------------------------------------------------
       rxt(i,k,fe22_ndx)=4.0e-30_r8*(temp(i,k)/300.0_r8)**(-1.52_r8)
!-----------------------------------------------------------------
!       ... Fep+ O2 + M -> FeO2p
!-----------------------------------------------------------------
       rxt(i,k,fe23_ndx)=8.3e-30_r8*(temp(i,k)/300.0_r8)**(-1.86_r8)
!-----------------------------------------------------------------
!       ... Fep+ e -> Fe (Changed on April 2013 Feng et al. JGR,2013)
!-----------------------------------------------------------------
!      rxt(i,k,fe25a_ndx)=1.0e-12_r8
       rxt(i,k,fe25a_ndx)=6.5e-12_r8*(temp(i,k)/300.0_r8)**(-0.51_r8)
!      rxt(:,k,fe25a_ndx)=1.2e-10_r8*temp(:ncol,k)**(-0.509_r8)
!-----------------------------------------------------------------
!       ... FeOp+ e -> Fe + O
!-----------------------------------------------------------------
!      rxt(i,k,fe25b_ndx) =3.0e-7_r8 *(200.0_r8/temp(i,k))**(0.5_r8)
! Bones et al. (2015)
       rxt(i,k,fe25b_ndx) =5.5e-7_r8 *(298.0_r8/temp(i,k))**(0.5_r8)

!-----------------------------------------------------------------
!       ... FeO2p+ e -> Fe + O2
!-----------------------------------------------------------------
      rxt(i,k,fe25c_ndx) = 3.0e-7_r8 *(200.0_r8/temp(i,k))**(0.5_r8)
!-----------------------------------------------------------------
!       ... FeN2p+ e -> Fe + N2
!-----------------------------------------------------------------
      rxt(i,k,fe25d_ndx) = 3.0e-7_r8 *(200.0_r8/temp(i,k))**(0.5_r8)
!-----------------------------------------------------------------
!       ... FeOH+ FeOH -> (FeOH)2
!-----------------------------------------------------------------
!     rxt(:,k,fe28_ndx) = 9.0e-10_r8
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!Note FOR CESM2 we have not added PMC scheme in
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... Mg chemistry reactions
!-----------------------------------------------------------------
!       ... MgO + O -> Mg + O2
!-----------------------------------------------------------------
      rxt(i,k,mg02_ndx) = 6.2e-10_r8 * (temp(i,k)/ 295.0_r8)**(1.0_r8/6.0_r8)
!-----------------------------------------------------------------
!       ... MgO2 + O -> MgO + O2
!-----------------------------------------------------------------
      rxt(i,k,mg04_ndx) = 7.9e-11_r8 * (temp(i,k)/ 295.0_r8)**(1.0_r8/6.0_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... MgO + H2O +M --> Mg(OH)2
!
!-----------------------------------------------------------------
!JPL.10-6Final15June2011.pdf   (Page, 2-2) equation in the top page
!-----------------------------------------------------------------
!        mg_vk8_0  = 10**(-32.7538_r8+7.984_r8*logt-2.127_r8*logt*logt)
!        mg_vk8_inf = 3.52e-10_r8 * exp(-334.01_r8/temp(i,k))
!        mg_vk8term1 = mg_vk8_0/(mg_vk8_inf/m(i,k))
!        mg_vk8term2 = mg_vk8_0/(1._r8 + mg_vk8term1)
!        mg_vk8term3= log10(mg_vk8term1)
!        mg_vk8term4= 1.0_r8 / (1.0_r8 + mg_vk8term3*mg_vk8term3)
!        rxt(i,k,mg05_ndx) = mg_vk8term2 * (0.28_r8)**mg_vk8term4
      rxt(i,k,mg05_ndx) = 1.1e-26_r8 * (temp(i,k)/ 200.0_r8)**(-1.59_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... MgO+ O2 +M --> MgO3-> Mg(OH)2
!-----------------------------------------------------------------
!        mg_vk9_0  = 10**(-28.046_r8+1.4226_r8*logt-0.683_r8*logt*logt)
!        mg_vk9_inf = 1.157e-10_r8 * exp(-219.27_r8/temp(i,k))
!        mg_vk9term1 = mg_vk9_0/(mg_vk9_inf/m(i,k))
!        mg_vk9term2 = mg_vk9_0/(1._r8 + mg_vk9term1)
!        mg_vk9term3= log10(mg_vk9term1)
!        mg_vk9term4= 1.0_r8 / (1.0_r8 + mg_vk9term3*mg_vk9term3)
!        rxt(i,k,mg07_ndx) = mg_vk9term2 * (0.34_r8)**mg_vk9term4
      rxt(i,k,mg07_ndx) = 3.8e-29_r8 * (temp(i,k)/ 200.0_r8)**(-1.59_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... MgO+ CO2 +M --> MgCO3->Mg(OH)2
!-----------------------------------------------------------------

!        mg_vk10_0  = 10**(-33.6985_r8+5.8268_r8*logt-1.4939_r8*logt*logt)
!        mg_vk10_inf = 6.79e-10_r8 * exp(-309.96_r8/temp(i,k))
!        mg_vk10term1 = mg_vk10_0/(mg_vk10_inf/m(i,k))
!        mg_vk10term2 = mg_vk10_0/(1._r8 + mg_vk10term1)
!        mg_vk10term3= log10(mg_vk10term1)
!        mg_vk10term4= 1.0_r8 / (1.0_r8 + mg_vk10term3*mg_vk10term3)
!        rxt(i,k,mg08_ndx) = mg_vk10term2 * (0.37_r8)**mg_vk10term4
      rxt(i,k,mg08_ndx) = 5.9e-29_r8 * (temp(i,k)/ 200.0_r8)**(-0.86_r8)

!-----------------------------------------------------------------
!       ... MgO2+ O2 +M --> O2MgO2->Mg(OH)2
!LEAVE IT OUT AFTER DISCUSSION WITH JOHN
!-----------------------------------------------------------------
!        mg_vk12_0  = 6.4e-27_r8*(temp(i,k)/300._r8)**(-3.312_r8)
!        mg_vk12_inf = 1.57e-10_r8 * (temp(i,k)/300._r8)**(0.167_r8)
!        mg_vk12term1 = mg_vk12_0/(mg_vk12_inf/m(i,k))
!        mg_vk12term2 = mg_vk12_0/(1._r8 + mg_vk12term1)
!        mg_vk12term3= log10(mg_vk12term1)
!        mg_vk12term4= 1.0_r8 / (1.0_r8 + mg_vk12term3*mg_vk12term3)
!        rxt(i,k,mg10_ndx) = mg_vk12term2 * (0.37_r8)**mg_vk12term4
      rxt(i,k,mg10_ndx) = 1.8e-26_r8 * (temp(i,k)/ 200.0_r8)**(-2.5_r8)
!-----------------------------------------------------------------
!       ... Mgp+ N2 + M--> Mgp.N2 + M
!-----------------------------------------------------------------
       rxt(i,k,mg20_ndx) = 1.8e-30*(temp(i,k)/200._r8)**(-1.72_r8)
!-----------------------------------------------------------------
!       ... Mgp+ O2 + M--> MgO2p + M
!-----------------------------------------------------------------
       rxt(i,k,mg22_ndx) = 2.4e-30*(temp(i,k)/200._r8)**(-1.65_r8)
!-----------------------------------------------------------------
!       ... Mgp+H2O+M --> MgH2Op +M
!-----------------------------------------------------------------
      rxt(i,k,mg27_ndx) = 2.3e-28*(temp(i,k)/200.0_r8)**(-2.53_r8)
!-----------------------------------------------------------------
!       ... Mgp+CO2+M --> MgCO2p +M
!-----------------------------------------------------------------
      rxt(i,k,mg28_ndx) = 4.6e-29*(temp(i,k)/200.0_r8)**(-1.42_r8)
!-----------------------------------------------------------------
!       ... Mgp+e --> Mg +hv
!-----------------------------------------------------------------
      rxt(i,k,mg29_ndx) = 3.3e-12*(temp(i,k)/200.0_r8)**(-0.6_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... MgXp+ e--> Mg + X (Here X can be O, O2, N2, CO2 or H2O)
!-----------------------------------------------------------------
       rxt(i,k,mg30_ndx) = 3.0e-7*sqrt(200._r8/temp(i,k))
       rxt(i,k,mg30b_ndx) = 3.0e-7*sqrt(200._r8/temp(i,k))
       rxt(i,k,mg31c_ndx) = 3.0e-7*sqrt(200._r8/temp(i,k))
       rxt(i,k,mg31d_ndx) = 3.0e-7*sqrt(200._r8/temp(i,k))
!-----------------------------------------------------------------
!!-----------------------------------------------------------------
!       ... Si chemistry reactions
!-----------------------------------------------------------------
!       ... Si + O2 -> SiO + O
!-----------------------------------------------------------------
      rxt(i,k,si01_ndx)=9.5e-11_r8+1.8e-10_r8*exp(-1.0_r8*temp(i,k)/115.0_r8)
!-----------------------------------------------------------------
!       ... SiO2 + H2O + M -> SiO3H2+M
!JPL.10-6Final15June2011.pdf   (Page, 2-2) equation in the top page
!-----------------------------------------------------------------
         vk4_0  = 5.54e-25_r8 * (temp(i,k)/ 200._r8)**(-3.69_r8)
         vk4_inf = 3.51e-10_r8 * exp(-30.18_r8/temp(i,k))
         vk4term1 = vk4_0/(vk4_inf/m(i,k))
         vk4term2 = vk4_0/(1._r8 + vk4term1)
         vk4term3= log10(vk4term1)
         vk4term4= 1.0_r8 / (1.0_r8 + vk4term3*vk4term3)
         rxt(i,k,si04_ndx) = vk4term2 * (0.188_r8)**vk4term4
!-----------------------------------------------------------------
!       ... SiO3H2+ H2O + M -> SiO4H4 +M
!-----------------------------------------------------------------
         vk5_0  = 2.71e-21_r8 * (temp(i,k)/ 200._r8)**(-4.22_r8)
         vk5_inf = 3.51e-10_r8 * exp(-30.18_r8/temp(i,k))
         vk5term1 = vk5_0/(vk5_inf/m(i,k))
         vk5term2 = vk5_0/(1._r8 + vk5term1)
         vk5term3= log10(vk5term1)
         vk5term4= 1.0_r8 / (1.0_r8 + vk5term3*vk5term3)
         rxt(i,k,si05_ndx) = vk5term2 * (0.249_r8)**vk5term4
!-----------------------------------------------------------------
!       ... SiO + OH -> SiO2 +  H
!-----------------------------------------------------------------
        vk6term=-4.92_r8-3.08_r8*log10(temp(i,k))+0.5_r8*log10(temp(i,k))**2_r8
       rxt(i,k,si06_ndx)=10**vk6term
!-----------------------------------------------------------------
!       ... SiO2 + H -> SiO +  OH
!-----------------------------------------------------------------
        vk7term=-19.83_r8+7.32_r8*log10(temp(i,k))-1.37_r8*log10(temp(i,k))**2_r8
       rxt(i,k,si07_ndx)=10**vk7term
!-----------------------------------------------------------------
!       ... Sip + H2O -> SiOHp +  H
!-----------------------------------------------------------------
       rxt(i,k,si12_ndx)=2.3e-10_r8*(temp(i,k)/300.0_r8)**(-0.5_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... Sip + O2_1D -> SiOp +  O
!-----------------------------------------------------------------
       rxt(i,k,si13_ndx)=3.6e-11_r8*(temp(i,k)/200.0_r8)**(-2.2_r8)
!-----------------------------------------------------------------
!       ... Si+ + O2 + M -> SiO2+ + M
!-----------------------------------------------------------------
       rxt(i,k,si14_ndx)=9.0e-29_r8*(temp(i,k)/200.0_r8)**(-2.34_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... SiO+ + H2O -> SiOH+ + OH
!-----------------------------------------------------------------
       rxt(i,k,si16_ndx)=2.9e-9_r8*(temp(i,k)/200.0_r8)**(-0.5_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... SiO+ + e -> Si + O
!       ... SiOH+ + e -> Si + OH
!       ... SiO2+ + e -> Si + O2
!-----------------------------------------------------------------
       rxt(i,k,si20b_ndx)=4.0e-7_r8*(temp(i,k)/300.0_r8)**(-0.5_r8)
       rxt(i,k,si20c_ndx)=4.0e-7_r8*(temp(i,k)/300.0_r8)**(-0.5_r8)
       rxt(i,k,si20d_ndx)=4.0e-7_r8*(temp(i,k)/300.0_r8)**(-0.5_r8)
!-----------------------------------------------------------------
!       ... SiO + O2p -> Si + O2
!-----------------------------------------------------------------
       rxt(i,k,si21_ndx)=4.2e-9_r8*(temp(i,k)/200.0_r8)**(-0.5_r8)
!!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... Sip + e -> Si + hv
!-----------------------------------------------------------------
       rxt(i,k,si23_ndx)=9.2e-12_r8*(temp(i,k)/200.0_r8)**(-0.56_r8)
!!-----------------------------------------------------------------

!!-----------------------------------------------------------------
!!       ... Ca chemistry reactions
!-----------------------------------------------------------------
!       ... CaO + O2 + M -> CaO2 + M
!-----------------------------------------------------------------
      rxt(i,k,ca04_ndx)=10_r8**(-57.10_r8+19.7_r8*logt-3.410_r8*logt*logt)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... CaO + H2O + M -> CaO2H2 + M
!-----------------------------------------------------------------
      rxt(i,k,ca06_ndx)=10_r8**(-23.39_r8+1.41_r8*logt-0.751_r8*logt*logt)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... CaO + O2 + M -> CaO3 + M
!-----------------------------------------------------------------
      rxt(i,k,ca08_ndx)=10_r8**(-42.19_r8+13.15_r8*logt-2.87_r8*logt*logt)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... CaO + CO2 + N2 -> CaCO3 + N2
!-----------------------------------------------------------------
      rxt(i,k,ca09_ndx)=10_r8**(-36.14_r8+9.24_r8*logt-2.191_r8*logt*logt)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... CaO3 + O -> CaO2 + O2
!This rate constant is changed to 2e-11 by John on 12 Jan 2017
!-----------------------------------------------------------------
!     rxt(i,k,ca12_ndx)=4.e-10_r8*(temp(i,k)/200.0_r8)**(0.5_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... CaO2 + O3-> CaO3 + O2
!-----------------------------------------------------------------
      rxt(i,k,ca13_ndx)=1.e-10_r8*(temp(i,k)/200.0_r8)**(0.5_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!ADDED REACTION
!       ... CaOH + O2+M-> OCaOH + M
!John's comment 7 January 2017
!ca_vk17b_0(:)  = 8.9e-26_r8*(300._r8/temp(:ncol,k))**(4.99_r8)
!-----------------------------------------------------------------
        ca_vk17b_0  = 8.9e-26_r8*(300._r8/temp(i,k))**(4.99_r8)
        ca_vk17b_inf = 1.5e-10_r8 * (temp(i,k)/300._r8)**(0.167_r8)
        ca_vk17bterm1 = ca_vk17b_0/(ca_vk17b_inf/m(i,k))
        ca_vk17bterm2 = ca_vk17b_0/(1._r8 + ca_vk17bterm1)
        ca_vk17bterm3= log10(ca_vk17bterm1)
        ca_vk17bterm4= 1.0_r8 / (1.0_r8 + ca_vk17bterm3*ca_vk17bterm3)
      rxt(i,k,ca17b_ndx) = ca_vk17bterm2 * (0.136_r8)**ca_vk17bterm4
!     rxt(i,k,ca17b_ndx) = 0._r8
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... Cap + O2 + M -> CaO2p + M
!John suggested to change on 11 June 2017
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!      rxt(i,k,ca23_ndx)=3_r8*10_r8**(-26.16_r8-1.113_r8*logt-0.0568_r8*logt*logt)
       rxt(i,k,ca23_ndx)=4.1e-29*((200.0_r8/temp(i,k))**2.37_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... Cap + H2O + M-> CaH2Op + M
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      rxt(i,k,ca25_ndx)=3_r8*10_r8**(-23.88_r8-1.823_r8*logt-0.063_r8*logt*logt)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... Cap + CO2 + M-> CaCO2p + M
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      rxt(i,k,ca26_ndx)=3_r8*10_r8**(-27.94_r8+2.204_r8*logt-1.124_r8*logt*logt)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... CaXp + e-> Ca + X
!CHANGE TO divided by 295K, 5 January 2017
!-----------------------------------------------------------------
!     rxt(i,k,ca31_ndx)=3.e-7_r8*sqrt(200.0_r8/temp(i,k))
      rxt(i,k,ca31_ndx)=3.e-7_r8*sqrt(295.0_r8/temp(i,k))
      rxt(i,k,ca31b_ndx)=3.e-7_r8*sqrt(295.0_r8/temp(i,k))
      rxt(i,k,ca31c_ndx)=3.e-7_r8*sqrt(295.0_r8/temp(i,k))
      rxt(i,k,ca31d_ndx)=3.e-7_r8*sqrt(295.0_r8/temp(i,k))
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... Cap + e-> Ca + hv
!-----------------------------------------------------------------
      rxt(i,k,ca32_ndx)=3.8e-12_r8*(temp(i,k)/200.0_r8)**(-0.9_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... Cap + N2 +M-> CaN2+ +M 
!John Suggested to add on 11 June 2017
!-----------------------------------------------------------------
       rxt(i,k,ca34_ndx)=2.3e-30*((200.0_r8/temp(i,k))**2.49_r8)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!       ... K chemistry reactions
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... K + O2 + M -> KO2 + M
!-----------------------------------------------------------------
      rxt(i,k,k03_ndx) = 1.3e-29_r8 * (temp(i,k)/ 200.0_r8)**(-1.23_r8)
!-----------------------------------------------------------------
!       ... KOH + CO2 +M  -> KHCO3 + M
!-----------------------------------------------------------------
       rxt(i,k,k11_ndx) = 7.1e-28*(temp(i,k)/200._r8)**(-4.2_r8)
!-----------------------------------------------------------------
!       ... Kp + N2 +M  -> KpN2 + M
!-----------------------------------------------------------------
       rxt(i,k,k22_ndx) = 2.3e-30*(temp(i,k)/200._r8)**(-2.39_r8)
!-----------------------------------------------------------------
!       ... Kp + O2 +M  -> KpO2 + M
!-----------------------------------------------------------------
       rxt(i,k,k23_ndx) = 1.2e-30*(temp(i,k)/200._r8)**(-2.12_r8)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       ... Kp + O +M  -> KpO + M
!-----------------------------------------------------------------
       rxt(i,k,k24_ndx) = 8.8e-32*(temp(i,k)/200._r8)**(-1.28_r8)
!-----------------------------------------------------------------
!       ... Kp + CO2 +M  -> KpCO2 + M
!-----------------------------------------------------------------
       rxt(i,k,k25_ndx) = 1.3e-29*(temp(i,k)/200._r8)**(-2.43_r8)
!-----------------------------------------------------------------
!       ... Kp + H2O +M  -> KpH2O + M
!-----------------------------------------------------------------
       rxt(i,k,k26_ndx) = 3.0e-29*(temp(i,k)/200._r8)**(-2.22_r8)
!-----------------------------------------------------------------
!       ... KpN2 + O  -> KpO + N2
!-----------------------------------------------------------------
       rxt(i,k,k27_ndx) = 2.9e-10*(temp(i,k)/200._r8)**(-0.17_r8)
!-----------------------------------------------------------------
!       ... KpO + N2  -> KpN2 + O
!-----------------------------------------------------------------
       rxt(i,k,k27b_ndx) = 2.5e-11*(temp(i,k)/200._r8)**(-0.55_r8)
!-----------------------------------------------------------------
!       ... KpN2 + CO2  -> KpCO2 + N2
!-----------------------------------------------------------------
       rxt(i,k,k28_ndx) = 4.8e-10*(temp(i,k)/200._r8)**(-0.88_r8)
!-----------------------------------------------------------------
!       ... KpO2 + O  -> KpO + O2
!-----------------------------------------------------------------
       rxt(i,k,k29_ndx) = 2.8e-10*(temp(i,k)/200._r8)**(-0.42_r8)
!-----------------------------------------------------------------
!       ... KpO + CO2  -> KpCO2 + O
!-----------------------------------------------------------------
       rxt(i,k,k30_ndx) = 7.1e-10*(temp(i,k)/200._r8)**(-0.21_r8)
!-----------------------------------------------------------------
!       ... KpO + H2O  -> KpH2O + O
!-----------------------------------------------------------------
       rxt(i,k,k31_ndx) = 7.1e-10*(temp(i,k)/200._r8)**(-1.90_r8)
!-----------------------------------------------------------------
!       ... KpO2 + N2  -> KpN2 + O2
!-----------------------------------------------------------------
       rxt(i,k,k32_ndx) = 1.6e-10*(temp(i,k)/200._r8)**(-0.90_r8)
!-----------------------------------------------------------------
!       ... KpO2 + CO2  -> KpCO2 + O2
!-----------------------------------------------------------------
       rxt(i,k,k33_ndx) = 1.5e-10*(temp(i,k)/200._r8)**(-1.94_r8)
!-----------------------------------------------------------------
!       ... KpO2 + H2O  -> KpH2O + O2
!-----------------------------------------------------------------
       rxt(i,k,k34_ndx) = 1.8e-9*(temp(i,k)/200._r8)**(-0.99_r8)
!-----------------------------------------------------------------
!       ... KpN2 + H2O  -> KpH2O + N2
!-----------------------------------------------------------------
       rxt(i,k,k35_ndx) = 2.4e-9*(temp(i,k)/200._r8)**(-0.45_r8)
!-----------------------------------------------------------------
!       ... KpCO2 + H2O  -> KpH2O + CO2
!-----------------------------------------------------------------
       rxt(i,k,k36_ndx) = 1.4e-9*(temp(i,k)/200._r8)**(-1.26_r8)
!-----------------------------------------------------------------
!       ... KpX + e  -> K + X (here X=N2, O2, CO2, H2O, O)
!-----------------------------------------------------------------
       rxt(i,k,k37a_ndx) = 1.0e-6*(temp(i,k)/200._r8)**(-0.5_r8)
       rxt(i,k,k37b_ndx) = rxt(i,k,k37a_ndx)
       rxt(i,k,k37c_ndx) = rxt(i,k,k37a_ndx)
       rxt(i,k,k37d_ndx) = rxt(i,k,k37a_ndx)
       rxt(i,k,k37e_ndx) = rxt(i,k,k37a_ndx)
!-----------------------------------------------------------------




        enddo
        enddo
!!END ADD BY WUHU FENG
!END ADD BY WF et al.

!-----------------------------------------------------------------
!	... tropospheric "aerosol" rate constants
!-----------------------------------------------------------------
     if ( het1_ndx > 0 .AND. (.NOT. usr_N2O5_aer_ndx > 0) ) then
         amas = 4._r8*pi*rm1**3*den/3._r8            ! each mean particle(r=0.1u) mass (g)
         do k = 1,pver
!-------------------------------------------------------------------------
! 	... estimate humidity effect on aerosols (from Shettle and Fenn, 1979)
!           xr is a factor of the increase aerosol radii with hum (hum=0., factor=1)
!-------------------------------------------------------------------------
            xr(:)     = .999151_r8 + relhum(:ncol,k)*(1.90445_r8 + relhum(:ncol,k)*(-6.35204_r8 + relhum(:ncol,k)*5.32061_r8))
!-------------------------------------------------------------------------
! 	... estimate sulfate particles surface area (cm2/cm3) in each grid
!-------------------------------------------------------------------------
            if ( carma_hetchem_feedback ) then
               sur(:ncol) = strato_sad(:ncol,k)
            else
               sur(:) = sulfate(:,k)*m(:,k)/avo*wso4 &              ! xform mixing ratio to g/cm3
                        / amas &                                    ! xform g/cm3 to num particels/cm3
                        * fare &                                    ! xform num particels/cm3 to cm2/cm3
                        * xr(:)*xr(:)                               ! humidity factor
            endif
!-----------------------------------------------------------------
!	... compute the "aerosol" reaction rates
!-----------------------------------------------------------------
!             k = gam * A * velo/4
!
!       where velo = sqrt[ 8*bk*T/pi/(w/av) ]
!             bk = 1.381e-16
!             av = 6.02e23
!             w  = 108 (n2o5)  HO2(33)  CH2O (30)  NH3(15)  
!
!       so that velo = 1.40e3*sqrt(T)  (n2o5)   gama=0.1
!       so that velo = 2.53e3*sqrt(T)  (HO2)    gama>0.2
!       so that velo = 2.65e3*sqrt(T)  (CH2O)   gama>0.022
!       so that velo = 3.75e3*sqrt(T)  (NH3)    gama=0.4
!--------------------------------------------------------
!-----------------------------------------------------------------
!	... use this n2o5 -> 2*hno3 only in tropopause
!-----------------------------------------------------------------
	    rxt(:,k,het1_ndx) = rxt(:,k,het1_ndx) &
                                +.25_r8 * gam1 * sur(:) * 1.40e3_r8 * sqrt( temp(:ncol,k) )
         end do
      end if

!lke++
!-----------------------------------------------------------------
!      ... CO tags
!-----------------------------------------------------------------
      if( usr_CO_OH_b_ndx > 0 ) then
         if( usr_COhc_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_COhc_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_COme_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_COme_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO01_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO01_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO02_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO02_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO03_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO03_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO04_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO04_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO05_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO05_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO06_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO06_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO07_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO07_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO08_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO08_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO09_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO09_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO10_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO10_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO11_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO11_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO12_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO12_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO13_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO13_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO14_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO14_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO15_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO15_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO16_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO16_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO17_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO17_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO18_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO18_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO19_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO19_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO20_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO20_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO21_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO21_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO22_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO22_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO23_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO23_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO24_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO24_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO25_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO25_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO26_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO26_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO27_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO27_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO28_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO28_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO29_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO29_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO30_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO30_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO31_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO31_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO32_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO32_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO33_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO33_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO34_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO34_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO35_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO35_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO36_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO36_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO37_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO37_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO38_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO38_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO39_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO39_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO40_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO40_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO41_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO41_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
         if( usr_CO42_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO42_OH_ndx) = rxt(:ncol,:,usr_CO_OH_b_ndx)
         end if
      end if
!lke--
!
! jfl : additional BAM removal reactions.  Zero out below the tropopause
!
      do l=1,num_strat_tau
!
         if ( usr_strat_tau_ndx(l) > 0 ) then
            do i=1,ncol
               rxt(i,tropchemlev(i)+1:pver,usr_strat_tau_ndx(l)) = 0._r8
            end do
         end if
!
      end do
!

      deallocate( sfc_array, dm_array )

  end subroutine usrrxt

      subroutine usrrxt_hrates( rxt, tempn, tempi, tempe, &
				h2ovmr, m, ncol, kbot )
!-----------------------------------------------------------------
!        ... set the user specified reaction rates for heating
!-----------------------------------------------------------------

      use shr_kind_mod,  only : r8 => shr_kind_r8
      use chem_mods,     only : rxntot
      use ppgrid,        only : pver, pcols

      implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
      integer, intent(in)     :: ncol                         ! number columns in chunk
      integer, intent(in)     :: kbot                         ! heating levels
      real(r8), intent(in)    :: tempn(pcols,pver)            ! neutral temperature (K)
      real(r8), intent(in)    :: tempi(pcols,pver)            ! ion temperature (K)
      real(r8), intent(in)    :: tempe(pcols,pver)            ! electron temperature (K)
      real(r8), intent(in)    :: m(ncol,pver)                 ! total atm density (1/cm^3)
      real(r8), intent(in)    :: h2ovmr(ncol,pver)            ! water vapor (vmr)
      real(r8), intent(inout) :: rxt(ncol,pver,rxntot)        ! gas phase rates
      
!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------

      integer  ::  k
      real(r8), dimension(ncol) :: &
                   tp, &
                   tinv, &
                   ko, &
                   kinf, &
                   fc

!-----------------------------------------------------------------
!	... o + o2 + m --> o3 + m
!-----------------------------------------------------------------
      do k = 1,kbot
         tinv(:ncol)       = 1._r8 / tempn(:ncol,k)
         tp(:)             = 300._r8 * tinv(:)
         rxt(:,k,usr_O_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8

!-----------------------------------------------------------------
!	... o + o + m -> o2 + m
!-----------------------------------------------------------------
         rxt(:,k,usr_O_O_ndx) = 2.76e-34_r8 * exp( 720.0_r8*tinv(:) )

!-----------------------------------------------------------------
!	... ho2 + ho2 --> h2o2
!	Note: this rate involves the water vapor number density
!-----------------------------------------------------------------
         ko(:)   = 3.0e-13_r8  * exp( 460._r8*tinv(:) )
         kinf(:) = 2.1e-33_r8 * m(:,k) * exp( 920._r8*tinv(:) )
         fc(:)   = 1._r8 + 1.4e-21_r8 * m(:,k) * h2ovmr(:,k) * exp( 2200._r8*tinv(:) )
         rxt(:,k,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)

      end do

!-----------------------------------------------------------------
! 	... the ionic rates
!-----------------------------------------------------------------
      if ( has_ion_rxts ) then
         level_loop2 :  do k = 1,kbot
            tp(:ncol)         = (2._r8*tempi(:ncol,k) + tempn(:ncol,k)) / ( 3._r8 * t0 )
            tp(:)             = max( min( tp(:),20._r8 ),1._r8 )
            rxt(:,k,ion1_ndx) = 2.82e-11_r8 + tp(:)*(-7.74e-12_r8 + tp(:)*(1.073e-12_r8  &
                 + tp(:)*(-5.17e-14_r8 + 9.65e-16_r8*tp(:))))
            tp(:ncol)         = (.6363_r8*tempi(:ncol,k) + .3637_r8*tempn(:ncol,k)) / t0
            tp(:)             = max( min( tp(:),trlim2 ),1._r8 )
            rxt(:,k,ion2_ndx) = 1.533e-12_r8 + tp(:)*(-5.92e-13_r8 + tp(:)*8.6e-14_r8)
            tp(:ncol)         = 2._r8 * t0 /(tempi(:ncol,k) + tempn(:ncol,k))
            where( tp(:ncol) < trlim3 )
               rxt(:,k,ion3_ndx)  = 1.4e-10_r8 * tp(:)**.44_r8
            elsewhere
               rxt(:,k,ion3_ndx)  = 5.2e-11_r8 / tp(:)**.2_r8
            endwhere
            tp(:ncol)          = t0 / tempe(:ncol,k)
            rxt(:,k,elec1_ndx) = 4.e-7_r8 * tp(:)**.85_r8
            rxt(:,k,elec3_ndx) = 1.8e-7_r8 * tp(:)**.39_r8
            where( tp(:ncol) < 4._r8 )
               rxt(:,k,elec2_ndx) = 2.7e-7_r8 * tp(:)**.7_r8
            elsewhere
               rxt(:,k,elec2_ndx) = 1.6e-7_r8 * tp(:)**.55_r8
            endwhere
         end do level_loop2
      endif
      end subroutine usrrxt_hrates

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine comp_exp( x, y, n )

    implicit none

    real(r8), intent(out) :: x(:)
    real(r8), intent(in)  :: y(:)
    integer,  intent(in)  :: n
    
#ifdef IBM
    call vexp( x, y, n )
#else
    x(:n) = exp( y(:n) )
#endif

  end subroutine comp_exp

  !-------------------------------------------------------------------------
  !  Heterogeneous reaction rates for uptake of a gas on an aerosol:
  !-------------------------------------------------------------------------
  function hetrxtrate( sfc, dm_aer, dg_gas, c_gas, gamma_gas ) result(rate)

    real(r8), intent(in) :: sfc(:)
    real(r8), intent(in) :: dm_aer(:)
    real(r8), intent(in) :: dg_gas
    real(r8), intent(in) :: c_gas
    real(r8), intent(in) :: gamma_gas
    real(r8) :: rate

    real(r8),allocatable :: rxt(:)
    integer :: n, i

    n = size(sfc)

    allocate(rxt(n))
    do i=1,n
       rxt(i) = sfc(i) / (0.5_r8*dm_aer(i)/dg_gas + (4._r8/(c_gas*gamma_gas)))
    enddo

    rate = sum(rxt)

    deallocate(rxt)

  endfunction hetrxtrate

  !-------------------------------------------------------------------------
  !  Heterogeneous reaction rates for uptake of a glyoxal gas on an aerosol:
  !-------------------------------------------------------------------------
  function hetrxtrate_gly( sfc, c_gas, gamma_gas ) result(rate)

    real(r8), intent(in) :: sfc(:)
    real(r8), intent(in) :: c_gas
    real(r8), intent(in) :: gamma_gas
    real(r8) :: rate

    real(r8),allocatable :: rxt(:)
    integer :: n, i

    n = size(sfc)

    allocate(rxt(n))
    do i=1,n
       rxt(i) =  0.25_r8 * c_gas * sfc(i) * gamma_gas
    enddo

    rate = sum(rxt)

    deallocate(rxt)

  endfunction hetrxtrate_gly


end module mo_usrrxt
