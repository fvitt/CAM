! carma_aero_gasaerexch.F90


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
!
! !MODULE: carma_aero_gasaerexch --- does carma aerosol gas-aerosol exchange for SOA
!
! !INTERFACE:
   module carma_aero_gasaerexch

! !USES:
  use shr_kind_mod,    only:  r8 => shr_kind_r8
  use chem_mods,       only:  gas_pcnst
  use ref_pres,        only:  top_lev => clim_modal_aero_top_lev
  use ppgrid,          only:  pcols, pver
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_info_by_bin, rad_cnst_get_bin_props_by_idx, &
                             rad_cnst_get_info_by_bin_spec

  use cam_logfile,      only: iulog
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_add_field, pbuf_get_field, pbuf_set_field, dtype_r8

  implicit none
  private
   ! Misc private data
  character(len=32), allocatable :: fldname(:)    ! names for interstitial output fields
  character(len=32), allocatable :: fldname_cw(:)    ! names for cloud_borne output fields

  save

! !PUBLIC MEMBER FUNCTIONS:
  public carma_aero_gasaerexch_sub, carma_aero_gasaerexch_init

! !PUBLIC DATA MEMBERS:
  integer, parameter :: pcnstxx = gas_pcnst

  ! description of bin aerosols
  integer, public, protected :: nspec_max = 0
  integer, public, protected :: nbins = 0
  integer, public, protected :: nsoa_vbs = 0
  integer, public, protected :: nsoa = 0
  integer, public, protected :: npoa = 0
  integer, public, protected, allocatable :: nspec(:)

  ! local indexing for bins
  integer, allocatable :: bin_idx(:,:) ! table for local indexing of modal aero number and mmr
  integer :: ncnst_tot                  ! total number of mode number conc + mode species

  integer :: mw_soa          = 250.
  integer :: fracvbs_idx          = 0
  integer, allocatable :: dqdtsoa_idx(:,:)
  integer, allocatable  :: cnsoa(:)         ! true if soa gas is a species and carma soa in bin
  integer, allocatable  :: cnpoa(:)         ! true if soa gas is a species and carma soa in bin
  integer, allocatable  :: l_soag(:)         ! true if soa gas is a species and carma soa in bin


  real(r8), protected, allocatable, public :: soa_equivso4_factor(:)
! this factor converts an soa volume to a volume of so4(+nh4)
! having same hygroscopicity as the soa

  real (r8), allocatable :: frac_vbs(:,:,:,:)
  !real (r8), allocatable :: dqdt_soa(:,:)

  real(r8), pointer :: dqdt_soa(:,:)

  logical, allocatable  :: do_soag_any(:)         ! true if soa gas is a species and carma soa in bin
! !DESCRIPTION: This module implements ...
!
! !REVISION HISTORY:
!
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! list private module data here

!EOC
!----------------------------------------------------------------------


  contains

!----------------------------------------------------------------------

      subroutine carma_aero_gasaerexch_init

!-----------------------------------------------------------------------
!
! Purpose:
!      gas-aerosol exchange SOAG <-> soa
!
! Author: Simone Tilmes
!
!-----------------------------------------------------------------------

!use modal_aero_data
!use modal_aero_rename

use cam_abortutils, only: endrun
use cam_history,    only: addfld, add_default, fieldname_len, horiz_only
use constituents,   only: pcnst, cnst_get_ind, cnst_name
use spmd_utils,     only: masterproc
use phys_control,   only: phys_getopts

implicit none

!-----------------------------------------------------------------------
! arguments

!-----------------------------------------------------------------------
! local
   integer  :: ipair, iq, iqfrm, iqfrm_aa, iqtoo, iqtoo_aa
   integer  :: j, jsoa,p
   integer  :: i, ii
   integer  :: l, l1, l2, lsfrm, lstoo, lunout
   integer  :: l_so4g
   integer  :: m, mm, mfrm, mtoo
   integer  :: n, nacc, nait, ns


   real(r8) :: tmp1, tmp2

   character(len=fieldname_len+3) :: fieldname
   character(len=32)              :: spectype
   character(len=32)              :: spec_name
   character(128)                 :: long_name
   character(128)                 :: msg
   character(8)                   :: unit
   character(len=2)               :: outsoa

   logical                        :: history_aerosol      ! Output the MAM aerosol tendencies
   logical                        :: history_aerocom    ! Output the aerocom history
   !-----------------------------------------------------------------------

        call phys_getopts( history_aerosol_out        = history_aerosol   )

!
   ! get info about the bin aerosols
   ! get nbins

    call rad_cnst_get_info( 0, nbins=nbins)

    allocate( nspec(nbins) )
    allocate( cnsoa(nbins) )
    allocate( cnpoa(nbins) )

    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=nspec(m))
    end do
    ! add plus one to include number, total mmr and nspec
    nspec_max = maxval(nspec) + 2

    ncnst_tot = nspec(1) + 2
    do m = 2, nbins
      ncnst_tot = ncnst_tot + nspec(m) + 2
    end do

   allocate(  bin_idx(nbins,nspec_max),               &
              do_soag_any(nbins),                     &
              frac_vbs(nbins, nsoa_vbs, pcols, pver), &
              ! dqdt_soa(pcols, pver),                  &
              fldname_cw(ncnst_tot),                  &
              fldname(ncnst_tot) )

   ! Local indexing compresses the mode and number/mass indicies into one index.
   ! This indexing is used by the pointer arrays used to reference state and pbuf
   ! fields.
   ! for CARMA we add number = 0, total mass = 1, and mass from each constituence into mm.
   ii = 0
   do m = 1, nbins
      do l = 1, nspec(m) + 2    ! do through nspec plus mmr and number
         ii = ii + 1
         bin_idx(m,l) = ii
      end do
   end do

  ! SOAG / SOA / POM information
  ! Define number of VBS bins (nsoa) based on number of SOAG chemistry species

    nsoa_vbs = 0
    do i = 1, pcnst
      if (cnst_name(i)(:4) == 'SOAG') then
             nsoa_vbs= nsoa_vbs + 1
!!$             write(iulog,*) 'soag = ', cnst_name(i)
      end if
    end do
    allocate( l_soag(nsoa_vbs) )
    nsoa_vbs = 0
    do i = 1, pcnst
      if (cnst_name(i)(:4) == 'SOAG') then
             nsoa_vbs= nsoa_vbs + 1
             l_soag(nsoa_vbs) = i
      end if
    end do
!!$    write(iulog,*) 'l_soag = ', l_soag

   fracvbs_idx      = pbuf_get_index('FRACVBS')

   ! identify number of SOA and POA in CARMA code (CARMA number cn)
    do m = 1, nbins
       cnsoa(m) = 0
       cnpoa(m) = 0
       do l = 1, nspec(m)
           mm = bin_idx(m, l)
           call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
           if (trim(spectype) == 's-organic') then
                    cnsoa(m) = cnsoa(m) + 1
      !              write(iulog,*) 'cnsoa,m,l = ', cnsoa(m),m,l
           end if
           if (trim(spectype) == 'p-organic') then
                    cnpoa(m) = cnpoa(m) + 1
      !              write(iulog,*) 'cnpoa = ', cnpoa(m),m,l
           end if
        end do
     end do
     ! some bins don't contain soa or poa
!!$     write(*,*) 'cnsoa = ', cnsoa
     !write(iulog,*) 'cnsoa = ', cnpoa
     nsoa= maxval(cnsoa)
     npoa= maxval(cnpoa)
!!$     write(iulog,*) 'nsoa = ', nsoa
!!$     write(iulog,*) 'npoa = ', npoa
!!$

    allocate( dqdtsoa_idx(nbins,nsoa)                       )
    do m = 1, nbins
       ns = 0
       do l = 1, nspec(m)
           call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
           if (trim(spectype) == 's-organic') then
              call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=spec_name)
              !write(*,*) 'set index DQDT_'//trim(spec_name)
              ns = ns + 1
              dqdtsoa_idx(m,ns) = pbuf_get_index('DQDT_'//trim(spec_name))
           end if
       end do
    end do

     ! set gas species indices and identify if soa carma is in bin
     !do jsoa = 1, nsoa
     !   id_soag(jsoa) = get_spc_ndx( soag(jsoa) )
     !   do_soag(jsoa) = .true.
     !end do

     do m = 1, nbins
       do_soag_any(m) = .false.
       if (cnsoa(m) .gt. 0) then
              do_soag_any(m) = .true.
       end if
     end do

   !if (nsoa_vbs.eq.nsoa) then
   ! nsoag and CARMA soa have same number, simple mapping, else will need to map to volatility bins
   !   nsoa = nsoag
   !   nsoa_vbs = nsoag
   !else
   !! nsoag may be 5 VBS bins, but CARMA SOA is only 1 VBS bin
   !   nsoa = cnsoa_max
   !   nsoa_vbs = nsoag
   !! define pbuf field that includes fraction of volatily bins here or in aero_model
   !end if
   !write(iulog,*) 'nsoa_vbs = ', nsoa_vbs



!   output results
!
!---------define history fields for new cond/evap diagnostics----------------------------------------

      fieldname=trim('qcon_gaex')
      long_name = trim('3D fields for SOA condensation')
      unit = 'kg/kg/s'
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
!!$      if ( masterproc ) write(*,'(3(a,3x))') 'qcon addfld', fieldname, unit

      do j = 1, nsoa_vbs
          write (outsoa, "(I2.2)") j
          fieldname=trim('qcon_gaex')//outsoa
          long_name = trim('3D fields for SOA condensation for VBS bin')//outsoa
          call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
          if ( history_aerosol ) then
            call add_default( fieldname,  1, ' ' )
          endif
      end do


      fieldname=trim('qevap_gaex')
      long_name = trim('3D fields for SOA evaporation')
      call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
      if ( history_aerosol ) then
         call add_default( fieldname,  1, ' ' )
      endif
!!$      if ( masterproc ) write(*,'(3(a,3x))') 'qevap addfld', fieldname, unit


      do j = 1, nsoa_vbs
          write (outsoa, "(I2.2)") j
          fieldname=trim('qevap_gaex')//outsoa
          long_name = trim('3D fields for SOA evaporation for VBS bin')//outsoa
          call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
          if ( history_aerosol ) then
             call add_default( fieldname,  1, ' ' )
          endif
      end do

!------------------------------------------------------------------------------

!  define history fields for basic gas-aer exchange
      do m = 1, nbins
        do l = 1, nspec(m) + 2    ! do through nspec plus mmr and number
         ii = bin_idx(m,l)
         if (l <= nspec(m) ) then   ! species
            call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=fldname(ii) )
            ! only write out SOA exchange here
            call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
            if (trim(spectype) == 's-organic') then
              fieldname= trim(fldname(ii)) // '_sfgaex1'
              long_name = trim(fldname(ii)) // ' gas-aerosol-exchange primary column tendency'
              unit = 'kg/m2/s'
              call addfld( fieldname, horiz_only, 'A', unit, long_name )
              if ( history_aerosol ) then
                call add_default( fieldname, 1, ' ' )
              endif
!!$              if ( masterproc ) write(*,'(3(a,3x))') 'gasaerexch addfld', fieldname, unit
            end if
        ! else if (l == nspec(m) + 1) then   ! mmr
        !    call rad_cnst_get_info_by_bin(0, m,  mmr_name=fieldname(ii), mmr_name_cw=fieldname_cw(ii))
        ! else if (l == nspec(m) + 2) then   !number
        !    call rad_cnst_get_info_by_bin(0, m, num_name=fieldname(ii), num_name_cw=fieldname_cw(ii))
         end if
        end do
      end do

      return

      end subroutine carma_aero_gasaerexch_init


!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  carma_aero_gasaerexch_sub --- ...
!
! !INTERFACE:
subroutine carma_aero_gasaerexch_sub(                            &
                        pbuf, lchnk,    ncol,     nstep,               &
                        loffset,  deltat,   mbar,                &
                        t,        pmid,     pdel,                &
                        qh2o,               troplev,             &
                        q,                  raervmr,             &
                        wetr_n                                   )

! !USES:
use cam_history,       only:  outfld, fieldname_len
use constituents,      only:  pcnst, cnst_name, cnst_get_ind
use mo_tracname,       only:  solsym
use physconst,         only:  gravit, mwdry, rair
use cam_abortutils,    only:  endrun
use spmd_utils,        only:  iam, masterproc
use time_manager,      only: is_first_step


implicit none

! !PARAMETERS:
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   integer,  intent(in)    :: troplev(pcols)       ! tropopause vertical index
   real(r8), intent(in)    :: deltat               ! time step (s)
   real(r8), intent(in)    :: mbar(ncol,pver)      ! mean wet atmospheric mass ( amu )

   real(r8), intent(inout) :: q(ncol,pver,pcnstxx) ! tracer mixing ratio (TMR) array
                                                   ! *** MUST BE  #/kmol-air for number
                                                   ! *** MUST BE mol/mol-air for mass
                                                   ! *** NOTE ncol dimension
   real(r8), intent(inout) :: raervmr (ncol,pver,ncnst_tot) ! aerosol mixing rations (vmr)
   real(r8), intent(in)    :: t(pcols,pver)        ! temperature at model levels (K)
   real(r8), intent(in)    :: pmid(pcols,pver)     ! pressure at model levels (Pa)
   real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: qh2o(pcols,pver)     ! water vapor mixing ratio (kg/kg)
   real(r8), intent(in)    :: wetr_n(pcols,pver,nbins) !wet geo. mean dia. (cm) of number distrib.

! !DESCRIPTION:
! this version does only do condensation for SOA for CARMA
!     method_soa=0 is no uptake
!     method_soa=1 is irreversible uptake done like h2so4 uptake
!     method_soa=2 is reversible uptake using subr carma_aero_soaexch
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1, ldiag4=-1
   integer, parameter :: method_soa = 2

   integer :: i, iq, itmpa
   integer :: idiagss
   integer :: j, jsoa, jpoa
   integer :: k,p
   integer :: l, l2, lb, lsfrm, lstoo
   integer :: l_soa(nsoa_vbs)
   integer :: mm, m, n, nn, niter, niter_max, ntot_soamode

   character(len=fieldname_len+3) :: fieldname
   character(len=100) :: msg !BSINGH - msg string for endrun calls
   character(len=32)              :: spectype
   character(len=32)              :: spec_name
   character(len=2)               :: outsoa

   real (r8) :: avg_uprt_nh4, avg_uprt_so4, avg_uprt_soa(nsoa_vbs)
   real (r8) :: deltatxx
   !real (r8) :: dqdt_soa(pcols,pver)
   real (r8) :: dqdt_soa_vbs(nbins,nsoa_vbs)
   real (r8) :: dqdt_soa_all(nbins,nsoa,pcols,pver)
   real (r8) :: dqdt_soag(nsoa_vbs)
   real (r8) :: fgain_soa(nbins,nsoa_vbs)
   real (r8) :: g0_soa(nsoa_vbs)
   real (r8) :: mw_poa_host              ! molec wght of poa used in host code
   real (r8) :: mw_soa_host              ! molec wght of poa used in host code
   real (r8) :: pdel_fac
   real (r8) :: num_bin(nbins,pcols,pver)          ! SOA from CARMA
   real (r8) :: soa_vbs(nbins,nsoa_vbs,pcols,pver)          ! SOA from CARMA
   real (r8) :: soa_c(nbins,nsoa,pcols,pver)          ! SOA from CARMA
   real (r8) :: poa_c(nbins,npoa,pcols,pver)          ! POA from CARMA
   real (r8) :: qold_poa(nbins,npoa)       ! POA from CARMA old
   real (r8) :: qold_soa(nbins,nsoa_vbs)   ! SOA on VBS bins  old
   real (r8) :: qnew_soa_vbs(nbins,nsoa_vbs)   ! SOA on VBS bins  new
   real (r8) :: qnew_soa(nbins)            ! SOA new for combined VBS bin new for combined VBS binss
   real (r8) :: qold_soag(nsoa_vbs)
   real (r8) :: sum_dqdt_soa(nsoa_vbs)     !   sum_dqdt_soa = soa tendency from soa  gas uptake (mol/mol/s)
   real (r8) :: sum_uprt_soa(nsoa_vbs)     ! total soa uptake rate over all bin, for each soa vbs bin
   real (r8) :: tmp1, tmp2, tmpa
   real (r8) :: tmp_kxt, tmp_pxt
   real (r8) :: uptkrate(nbins,pcols,pver)
   real (r8) :: uptkratebb(nbins)
   real (r8) :: uptkrate_soa(nbins,nsoa_vbs)
                ! gas-to-aerosol mass transfer rates (1/s)

   !integer, parameter :: nsrflx = 2     ! last dimension of qsrflx
   integer, parameter :: nsrflx = 1     ! only one dimension of qsrflx, no renaming or changes in size for CARMA currently
   real(r8) :: dqdt(ncol,pver,pcnstxx)  ! TMR "delta q" array - NOTE dims
   real(r8) :: qsrflx(pcols,nbins,nsoa)
                              ! process-specific column tracer tendencies
                              ! (1=gas condensation)
   real(r8) :: qcon_vbs(nsoa_vbs,pcols,pver)
   real(r8) :: qevap_vbs(nsoa_vbs,pcols,pver)
   real(r8) :: qcon(pcols,pver)
   real(r8) :: qevap(pcols,pver)
   real(r8) :: total_soag
   real(r8) :: soag(nsoa_vbs)

   real(r8), pointer :: frac_vbs(:,:,:,:)   ! fraction of vbs SOA bins to total SOA

!  following only needed for diagnostics
   real(r8) :: qold(ncol,pver,pcnstxx)  ! NOTE dims
   real(r8) :: qnew(ncol,pver,pcnstxx)  ! NOTE dims
   real(r8) :: qdel(ncol,pver,pcnstxx)  ! NOTE dims
   real(r8) :: dqdtsv1(ncol,pver,pcnstxx)

   type(physics_buffer_desc), pointer :: pbuf(:)

!----------------------------------------------------------------------

!  map CARMA soa to working soa(nbins,nsoa)

   call pbuf_get_field(pbuf, fracvbs_idx,  frac_vbs )



   do m = 1, nbins      ! main loop over aerosol bins
      if (do_soag_any(m)) then  ! only bins that contain soa
         n = 0
         nn = 0
         do l = 1, nspec(m)
              mm = bin_idx(m, l)
              call rad_cnst_get_bin_props_by_idx(0, m, l, spectype=spectype)
              if (trim(spectype) == 's-organic') then
                n = n + 1
                soa_c(m,n,:,:) = raervmr(:,:,mm)
              end if
              if (trim(spectype) == 'p-organic') then
                nn = nn + 1
                poa_c(m,nn,:,:) = raervmr(:,:,mm)
              end if
         end do
         if (npoa .gt. 1) then
             call endrun( 'carma_aero_gasaerexch_sub error: CARMA currently only supports 1 POA element' )
         end if


         if (nsoa_vbs.eq.nsoa) then
           soa_vbs(:,:,:,:) = soa_c(:,:,:,:)
         else
           if (nsoa.eq.1) then
             if (is_first_step()) then
             !first time step initialization only
               do k=top_lev,pver
                  do i=1,ncol
                     total_soag = 0.0_r8
                     do j = 1, nsoa_vbs
                        soag(j) = q(i,k,l_soag(j))
                        total_soag = total_soag + soag(j)
                     end do
                     if (total_soag .gt. 0.0_r8) then
                       do j= 1, nsoa_vbs
                           frac_vbs(m,j,i,k) = soag(j)/total_soag
                       end do
                    !   write(*,*) 'frac_vbs first timestep', frac_vbs(m,:,i,k)
                     end if
                   end do
                end do
             !write(iulog,*) 'first time step done'
             end if
             !write(iulog,*) 'second time step done'
             ! end first time step, after that use fraction from previous time step
             do k=top_lev,pver
               do i=1,ncol
                 do j= 1, nsoa_vbs
                    soa_vbs(m,j,i,k) = frac_vbs(m,j,i,k)*soa_c(m,nsoa,i,k)
                    !write(*,*) 'soa_vbs derived ', j, m, minval(soa_vbs(m,j,:ncol,:)),maxval(soa_vbs(m,j,:ncol,:))
                 end do
               end do
             end do
           else
             ! error message this code only works if SOAG and SOA CARMA have the same number of species, or if SOA CARMA has only one species.
             call endrun( 'carma_aero_gasaerexch_sub error in number of SOA species' )
           end if

         end if

       ! get bin number for all aerosols
       l = nspec(m) + 2
       mm = bin_idx(m, l)
       num_bin(m,:,:) = raervmr(:,:,mm)
      end if
  end do


! SOA will be updated in CARMA

! zero out tendencies and other
   dqdt(:,:,:) = 0.0_r8
   qsrflx(:,:,:) = 0.0_r8

!-------Initialize evap/cond diagnostics (ncols x pver)-----------
   qcon_vbs(:,:,:) = 0.0_r8
   qevap_vbs(:,:,:) = 0.0_r8
   qcon(:,:) = 0.0_r8
   qevap(:,:) = 0.0_r8
!---------------------------------------------------
   !write(iulog,*) 'gas_aer_uptkrates start '
! compute gas-to-aerosol mass transfer rates
! check if only number is needed for this calculatuion!
   call gas_aer_uptkrates( ncol,       loffset,                &
                           num_bin,          t,          pmid,       &
                           wetr_n,                   uptkrate    )

   !      write(iulog,*) 'gas_aer_uptkrates done '

! use this for tendency calcs to avoid generating very small negative values
   deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)


         !write(iulog,*) 'start k and i loop '

   dqdt_soa_all(:,:,:,:) = 0.0_r8
   do k=top_lev,pver
      do i=1,ncol
         sum_uprt_soa(:) = 0.0_r8
         do n = 1, nbins
           if (do_soag_any(n)) then  ! only bins that contain soa
                 uptkratebb(n) = uptkrate(n,i,k)
                 if (npoa .gt. 0) then
                   do j = 1, npoa
                     qold_poa(n,j) = poa_c(n,j,i,k)
                   end do
                 else
                   qold_poa(n,j) = 0.0_r8
                 end if
         !        write(iulog,*) 'jsoa '
                 do jsoa = 1, nsoa_vbs
                  ! 0.81 factor is for gas diffusivity (soa/h2so4)
                  ! (differences in fuch-sutugin and accom coef ignored)
                   fgain_soa(n,jsoa) = uptkratebb(n)*0.81_r8
                   uptkrate_soa(n,jsoa) = fgain_soa(n,jsoa)
                   sum_uprt_soa(jsoa) = sum_uprt_soa(jsoa) + fgain_soa(n,jsoa)
                   qold_soa(n,jsoa) = soa_vbs(n,jsoa,i,k)
                 end do
          else
            qold_poa(n,:) = 0.0_r8
            qold_soa(n,:) = 0.0_r8
            fgain_soa(n,:) = 0.0_r8
          end if
          if (maxval(uptkrate_soa(n,:)) .gt. 0.0_r8) then
           !write(*,*) 'uptkrate_soa ', n, uptkrate_soa(n,:)
          end if
         end do ! n
         !write(*,*) 'sum_uprt_soa', sum_uprt_soa
         !write(iulog,*) 'end n '

         do jsoa = 1, nsoa_vbs
            if (sum_uprt_soa(jsoa) > 0.0_r8) then
               do n = 1, nbins
                  if (do_soag_any(n)) then  ! only bins that contain soa
                    fgain_soa(n,jsoa) = fgain_soa(n,jsoa) / sum_uprt_soa(jsoa)
                  end if
               end do
            end if
         end do
        !write(*,*) 'fgain_soa ', minval(fgain_soa(:,:)), maxval(fgain_soa(:,:))

!   uptake amount (fraction of gas uptaken) over deltat
         do jsoa = 1, nsoa_vbs
         !   write(iulog,*) 'sum_uprt_soa(jsoa) ', sum_uprt_soa(jsoa)
         !   write(iulog,*) 'deltatxx ', deltatxx
            avg_uprt_soa(jsoa) = (1.0_r8 - exp(-deltatxx*sum_uprt_soa(jsoa)))/deltatxx
         end do
        !write(*,*) 'avg_uprt_soa ', avg_uprt_soa

!   sum_dqdt_soa = soa_a tendency from soa   gas uptake (mol/mol/s)

         !write(iulog,*) 'sum_dqdt_soa '
         do jsoa = 1, nsoa_vbs
               sum_dqdt_soa(jsoa) = q(i,k,l_soag(jsoa)) * avg_uprt_soa(jsoa)
         end do

         !write(iulog,*) 'method_soa = ', method_soa

         if (method_soa > 1) then
!   compute TMR tendencies for soag and soa interstial aerosol
!   using soa parameterization
            niter_max = 1000
            dqdt_soa_vbs(:,:) = 0.0_r8
            dqdt_soag(:) = 0.0_r8
            do jsoa = 1, nsoa_vbs
               qold_soag(jsoa) = q(i,k,l_soag(jsoa))
            end do

         ! if (maxval(qold_soa) .gt. 0.0_r8 .or. maxval(qold_soag) .gt. 0.0_r8) then
         !   write(*,*) 'qold_soag ',qold_soag
         !   write(*,*) 'qold_soa before soaexch ',  qold_soa
         !   write(*,*) 'qold_poa before soaexch ',  qold_poa

            mw_poa_host = 12.0_r8
            mw_soa_host = 250.0_r8

         !write(iulog,*) 'carma_aero_soaexch '

            call carma_aero_soaexch( deltat, t(i,k), pmid(i,k), &
                 niter, niter_max, nbins, nsoa_vbs, npoa, &
                 mw_poa_host, mw_soa_host, &
                 qold_soag, qold_soa, qold_poa, uptkrate_soa, &
                 dqdt_soag, dqdt_soa_vbs )
            sum_dqdt_soa(:) = -dqdt_soag(:)

            !write(iulog,*) 'dqdt_soag(:)',dqdt_soag(:)
            !write(*,*) 'dqdt_soa_vbs',dqdt_soa_vbs
            !write(*,*) 'sum_dqdt_soa',sum_dqdt_soa
          !  end if

         else if ( method_soa .eq. 1) then
!   compute TMR tendencies for soa interstial aerosol
!   due to simple gas uptake

           do n = 1, nbins
             if (do_soag_any(n) ) then
                 do jsoa = 1, nsoa_vbs
                   dqdt_soa_vbs(n,jsoa) = fgain_soa(n,jsoa)*sum_dqdt_soa(jsoa)
                 end do
              else
                 dqdt_soa_vbs(:,:) = 0.0_r8
              end if
            end do
         else ! method_soa is neither 1 nor 2, no uptake
            dqdt_soa_vbs(:,:) = 0.0_r8
         end if

         !write(iulog,*) 'update soa'
         !  update soa to calcuate fractions (state variables and pbuf is not updated for SOA, will be done in CARMA)
         pdel_fac = pdel(i,k)/gravit
         qnew_soa(:) =0.0_r8
         qnew_soa_vbs(:,:) =0.0_r8

         do n = 1, nbins
            if ( do_soag_any(n) ) then
               if (nsoa.eq.nsoa_vbs) then
                     do jsoa = 1, nsoa_vbs
                        qsrflx(i,n,jsoa) = qsrflx(i,n,jsoa) + dqdt_soa_vbs(n,jsoa)*pdel_fac
                        dqdt_soa_all(n,nsoa,i,k) = dqdt_soa_vbs(n,jsoa) !  sum up for different volatility bins
                     end do
               else if (nsoa.eq.1) then
                 !write(iulog,*) 'gasaer_exch:  nsoa.eq.1, n', n
                  do jsoa = 1, nsoa_vbs
                     dqdt_soa_all(n,nsoa,i,k) = dqdt_soa_all(n,nsoa,i,k) + dqdt_soa_vbs(n,jsoa) !  sum up for different volatility bins
                   !  write(iulog,*) 'dqdt_soa_all, dqdt_soa_vbs :',  n, nsoa, jsoa, dqdt_soa_all(n,nsoa,i,k), dqdt_soa_vbs(n,jsoa)
                  end do
                  do jsoa = 1, nsoa_vbs
                     qsrflx(i,n,nsoa) = qsrflx(i,n,nsoa) + dqdt_soa_vbs(n,jsoa)*pdel_fac
                     qnew_soa_vbs(n,jsoa) = qold_soa(n,jsoa) + dqdt_soa_vbs(n,jsoa)*deltat
                     qnew_soa(n) = qnew_soa(n) + qnew_soa_vbs(n,jsoa) ! derive new fraction of SOA bin contributions
                  end do
                  do jsoa = 1, nsoa_vbs
                     if (qnew_soa(n) .gt. 0.0_r8) then
                       frac_vbs(n,jsoa,i,k) = qnew_soa_vbs(n,jsoa) / qnew_soa(n)
                       !write(iulog,*) 'frac_vbs:  n, jsoa', n, jsoa, frac_vbs(n,jsoa,i,k)
                     end if
                  end do
               else
                   call endrun( 'carma_aero_gasaerexch_sub error' )
               end if

!------- Add code for condensation/evaporation diagnostics sum of all bin---
               do jsoa = 1, nsoa_vbs
                   if (dqdt_soa_vbs(n,jsoa).ge.0.0_r8) then
                           qcon_vbs(jsoa,i,k)=dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                           qcon(i,k)=qcon(i,k)+dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                   else if (dqdt_soa_vbs(n,jsoa).lt.0.0_r8) then
                           qevap_vbs(jsoa,i,k)=dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                           qevap(i,k)=qevap(i,k)+dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                   endif
               end do
!---------------------------------------------------------------------------------------------------------------------
           end if
         end do ! n

         !write(iulog,*) 'update soag tendencies'
!   compute TMR tendencies for SAOG gas
!   due to simple gas uptake
         do jsoa = 1, nsoa
               !write(iulog,*) 'jsoa',jsoa
               !write(iulog,*) ' dqdt', dqdt(i,k,l_soag(jsoa))
               !write(iulog,*) ' sum_dqdt_soa', -sum_dqdt_soa(jsoa)
               dqdt(i,k,l_soag(jsoa)) = -sum_dqdt_soa(jsoa)
! dqdt for gas is negative of the sum of dqdt for aerosol soa species in each mode: Manish
              !ADD OUTFLD for SAOG loss
              ! qsrflx_soag(i,jsoa) = qsrflx(i,jsoa) + dqdt(i,k,l_soag(jsoa))*pdel_fac
              ! call outfld( fieldname, qsrflx_soag(:,j), pcols, lchnk )
         end do

      end do   ! "i = 1, ncol"
   end do     ! "k = top_lev, pver"
         !write(iulog,*) ' i, k loop done'


!  This applies dqdt tendencies for SOAG only , soa is done in CARMA
!  apply the dqdt to update q
!
        ! write(iulog,*) 'soag tendencies'
   do jsoa = 1, nsoa_vbs
         do k = top_lev, pver
            do i = 1, ncol
               q(i,k,l_soag(jsoa)) = q(i,k,l_soag(jsoa)) + dqdt(i,k,l_soag(jsoa))*deltat
            end do
         end do
   end do


!-----Outfld for condensation/evaporation------------------------------
      call outfld(trim('qcon_gaex'), qcon(:,:), pcols, lchnk )
      call outfld(trim('qevap_gaex'), qevap(:,:), pcols, lchnk )
      do jsoa = 1, nsoa_vbs
        write (outsoa, "(I2.2)") jsoa
        call outfld(trim('qcon_gaex')//outsoa, qcon_vbs(jsoa,:,:), pcols, lchnk )
        call outfld(trim('qevap_gaex')//outsoa, qevap_vbs(jsoa,:,:), pcols, lchnk )
      end do
!-----------------------------------------------------------------------
!   do history file of column-tendency fields over SOA fields (as defined in CARMA) and set pointer
   do m = 1, nbins
      if (do_soag_any(m)) then
         j  = 0
         do l = 1, nspec(m)
            mm = bin_idx(m,l)
            call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
            if (trim(spectype) == 's-organic') then
               j = j + 1
               fieldname= trim(fldname(mm)) // '_sfgaex1'
               do i = 1, ncol
                   qsrflx(i,m,j) = qsrflx(i,m,j)*(mw_soa/mwdry)
               end do
               call outfld( fieldname, qsrflx(:,m,j), pcols, lchnk )

               !set pointer field
               call pbuf_get_field(pbuf, dqdtsoa_idx(m,j),  dqdt_soa )

               ! NEEDS TO BE CHECKED! What units are required for CARMA??
               ! soa in vmr, may need to be transformed to  mass mixing ratios using wet mass not dry mass
               ! convert to kg/kg/s currently in vmr
               do k=top_lev,pver
                  dqdt_soa(:ncol,k) =  dqdt_soa_all(m,j,:ncol,k)*(mw_soa/mbar(:ncol,k))
               end do
               ! write(*,*) 'set pbuf dqdtsoa m ', m, minval(dqdt_soa), maxval(dqdt_soa)
               !call pbuf_set_field(pbuf, dqdtsoa_idx(m,j), dqdt_soa(:ncol,:))
            end if
         end do ! l = ...
      end if
   end do ! l = ...
   ! set new frac_vbs field
   call pbuf_set_field(pbuf, fracvbs_idx,  frac_vbs )

   return
   end subroutine carma_aero_gasaerexch_sub


!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine gas_aer_uptkrates( ncol,       loffset,                &
                              num_bin,          t,          pmid,       &
                              wetr,                     uptkrate    )

!
!                         /
!   computes   uptkrate = | dx  dN/dx  gas_conden_rate(Dp(x))
!                         /
!   using Gauss-Hermite quadrature of order nghq=2
!
!       Dp = particle diameter (cm)
!       x = ln(Dp)
!       dN/dx = log-normal particle number density distribution
!       gas_conden_rate(Dp) = 2 * pi * gasdiffus * Dp * F(Kn,ac)
!           F(Kn,ac) = Fuchs-Sutugin correction factor
!           Kn = Knudsen number
!           ac = accomodation coefficient
!

use physconst, only: mwdry, rair

implicit none


   integer,  intent(in) :: ncol                 ! number of atmospheric column
   integer,  intent(in) :: loffset
   real(r8), intent(in) :: t(pcols,pver)        ! Temperature in Kelvin
   real(r8), intent(in) :: pmid(pcols,pver)     ! Air pressure in Pa
   real(r8), intent(in) :: wetr(pcols,pver,nbins)
   real(r8), intent(in) :: num_bin(nbins,pcols,pver)

   real(r8), intent(out) :: uptkrate(nbins,pcols,pver)
                            ! gas-to-aerosol mass transfer rates (1/s)


! local
   integer, parameter :: nghq = 2
   integer :: i, iq, k, l1, l2, la, n

   ! Can use sqrt here once Lahey is gone.
   real(r8), parameter :: tworootpi = 3.5449077_r8
   real(r8), parameter :: root2 = 1.4142135_r8
   real(r8), parameter :: beta = 2.0_r8

   real(r8) :: aircon
   real(r8) :: const
   real(r8) :: dp, dum_m2v
   real(r8) :: dryvol_a(pcols,pver)
   real(r8) :: gasdiffus, gasspeed
   real(r8) :: freepathx2, fuchs_sutugin
   real(r8) :: knudsen
   real(r8) :: lndp, lndpgn, lnsg
   real(r8) :: num_a
   real(r8) :: rhoair
   real(r8) :: sumghq
   real(r8), save :: xghq(nghq), wghq(nghq) ! quadrature abscissae and weights

   data xghq / 0.70710678_r8, -0.70710678_r8 /
   data wghq / 0.88622693_r8,  0.88622693_r8 /


         !write(iulog,*) 'gas uptake start '
! outermost loop over all bins
   do n = 1, nbins

! loops k and i
     do k=top_lev,pver
      do i=1,ncol
       if (wetr(i,k,n) .gt. 0.0_r8) then

         !write(iulog,*) 'calc rair '
         rhoair = pmid(i,k)/(rair*t(i,k))   ! (kg-air/m3)
!        aircon = 1.0e3*rhoair/mwdry        ! (mol-air/m3)

!!   "bounded" number conc. (#/m3)
!        num_a = dryvol_a(i,k)*v2ncur_a(i,k,n)*aircon

!   number conc. (#/m3) -- note q(i,k,numptr) is (#/kmol-air)
!   so need aircon in (kmol-air/m3)
!   num_bin = number per bin
         !write(iulog,*) 'aircon '
         aircon = rhoair/mwdry              ! (kmol-air/m3)
         num_a = num_bin(n,i,k)*aircon
         !num_a = num_bin(n,i,k)

!   gasdiffus = h2so4 gas diffusivity from mosaic code (m^2/s)
!               (pmid must be Pa)
         !write(iulog,*) 'gasdiffus'
         gasdiffus = 0.557e-4_r8 * (t(i,k)**1.75_r8) / pmid(i,k)
!   gasspeed = h2so4 gas mean molecular speed from mosaic code (m/s)
         gasspeed  = 1.470e1_r8 * sqrt(t(i,k))
!   freepathx2 = 2 * (h2so4 mean free path)  (m)
         !write(iulog,*) 'freepathx2'
         freepathx2 = 6.0_r8*gasdiffus/gasspeed

         !lnsg   = log( sigmag_amode(n) )
         !lndpgn = log( dgncur_awet(i,k,n) )   ! (m)
         !const  = tworootpi * num_a * exp(beta*lndpgn + 0.5_r8*(beta*lnsg)**2)

!st      CARMA; diameter: 2*wetr ; convert wetr from cm to m
         ! dp = 2.0_r8 * 0.01_r8 * wetr(i,k,n)
         dp =  0.01_r8 * wetr(i,k,n)
         !dp =  0.01_r8 * wetr(i,k,n)
!st mam assumes wetr in cm
!         dp = 2.0_r8 * wetr(i,k,n)
         !write(iulog,*) 'const'
         const = tworootpi * num_a * 2.0_r8 * dp
         ! gas_conden_rate(Dp) = const *  gasdiffus *  F(Kn,ac)
         !   knudsen number
         knudsen = freepathx2/dp
         fuchs_sutugin = (0.4875_r8*(1._r8 + knudsen)) /   &
                            (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)
         uptkrate(n,i,k) = const * gasdiffus * fuchs_sutugin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!   sum over gauss-hermite quadrature points
!         sumghq = 0.0_r8
!
!         do iq = 1, nghq

!  Changed by Manish Shrivastava on 7/17/2013 to use accom=1; because we do not know better
!   following assumes accomodation coefficient = ac = 1. instead 0.65 ! answer change needs to be tested
!   (Adams & Seinfeld, 2002, JGR, and references therein)
!           fuchs_sutugin = (0.75*ac*(1. + knudsen)) /
!                           (knudsen*(1.0 + knudsen + 0.283*ac) + 0.75*ac)
!            fuchs_sutugin = (0.4875_r8*(1._r8 + knudsen)) /   &
!                            (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)
!            sumghq = sumghq + wghq(iq)*dp*fuchs_sutugin/(dp**beta)
!         end do
!         uptkrate(n,i,k) = const * gasdiffus * sumghq
!         !write(iulog,*) 'uptkrate ', uptkrate(n,i,k)
         else
          uptkrate(n,i,k) = 0.0_r8
          !write(iulog,*) 'uptkrate not calculated, wetr = 0.'
         end if

      end do   ! "do i = 1, ncol"
      end do   ! "do k = 1, pver"

     !write(iulog,*) 'n, uptkrate   = ',n, minval(uptkrate(n,:,:)),maxval(uptkrate(n,:,:)), minval(wetr(:,:,n)),maxval(wetr(:,:,n)), minval(num_bin(:,:,n)), maxval(num_bin(:,:,n))

   end do   ! "do n = 1, nbins"


   return
   end subroutine gas_aer_uptkrates

!----------------------------------------------------------------------
      subroutine carma_aero_soaexch( dtfull, temp, pres, &
          niter, niter_max, nbins, ntot_soaspec, ntot_poaspec,  &
          mw_poa_host, mw_soa_host, &
          g_soa_in, a_soa_in, a_poa_in, xferrate_in, &
          g_soa_tend, a_soa_tend )
!         g_soa_tend, a_soa_tend, g0_soa, idiagss )

!-----------------------------------------------------------------------
!
! Purpose:
!
! calculates condensation/evaporation of "soa gas"
! to/from multiple aerosol modes in 1 grid cell
!
! key assumptions
! (1) ambient equilibrium vapor pressure of soa gas
!     is given by p0_soa_298 and delh_vap_soa
! (2) equilibrium vapor pressure of soa gas at aerosol
!     particle surface is given by raoults law in the form
!     g_star = g0_soa*[a_soa/(a_soa + a_opoa)]
! (3) (oxidized poa)/(total poa) is equal to frac_opoa (constant)
!
!
! Author: R. Easter and R. Zaveri
! Additions to run with multiple BC, SOA and POM's: Shrivastava et al., 2015
!-----------------------------------------------------------------------

      use mo_constants, only: rgas ! Gas constant (J/K/mol)

      implicit none

      real(r8), intent(in)  :: dtfull                     ! full integration time step (s)
      real(r8), intent(in)  :: temp                       ! air temperature (K)
      real(r8), intent(in)  :: pres                       ! air pressure (Pa)
      integer,  intent(out) :: niter                      ! number of iterations performed
      integer,  intent(in)  :: niter_max                  ! max allowed number of iterations
      integer,  intent(in)  :: nbins                      ! number of bins
      !integer,  intent(in)  :: ntot_soamode               ! number of modes having soa (here set equal to number of all bins)
      integer,  intent(in)  :: ntot_poaspec               ! number of poa species
      integer,  intent(in)  :: ntot_soaspec               ! number of soa species
      real(r8), intent(in)  :: mw_poa_host                ! molec wght of poa used in host code
      real(r8), intent(in)  :: mw_soa_host                ! molec wght of soa used in host code
      real(r8), intent(in)  :: g_soa_in(ntot_soaspec)               ! initial soa gas mixrat (mol/mol at host mw)
      !ntot soaspec = nbin, ntot_soamode = 5
      real(r8), intent(in)  :: a_soa_in(nbins,ntot_soaspec)    ! initial soa aerosol mixrat (mol/mol at host mw)
      real(r8), intent(in)  :: a_poa_in(nbins,ntot_poaspec)    ! initial poa aerosol mixrat (mol/mol at host mw)
      real(r8), intent(in)  :: xferrate_in(nbins,ntot_soaspec) ! gas-aerosol mass transfer rate (1/s)
      real(r8), intent(out) :: g_soa_tend(ntot_soaspec)             ! soa gas mixrat tendency (mol/mol/s at host mw)
      real(r8), intent(out) :: a_soa_tend(nbins,ntot_soaspec)  ! soa aerosol mixrat tendency (mol/mol/s at host mw)
!     integer,  intent(in)  :: idiagss

      integer :: ll
      integer :: m,k

      logical :: skip_soamode(nbins)   ! true if this bin does not have soa

      real(r8), parameter :: a_min1 = 1.0e-40_r8
      real(r8), parameter :: g_min1 = 1.0e-40_r8
      real(r8), parameter :: alpha = 0.05_r8     ! parameter used in calc of time step
      real(r8), parameter :: dtsub_fixed = -1.0_r8  ! fixed sub-step for time integration (s)

      real(r8) :: a_ooa_sum_tmp(nbins)          ! total ooa (=soa+opoa) in a bin
      real(r8) :: a_opoa(nbins)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa(nbins,ntot_soaspec)     ! soa aerosol mixrat (mol/mol at actual mw)
      real(r8) :: a_soa_tmp(nbins,ntot_soaspec) ! temporary soa aerosol mixrat (mol/mol)
      real(r8) :: beta(nbins,ntot_soaspec)      ! dtcur*xferrate
      real(r8) :: delh_vap_soa(ntot_soaspec)           ! delh_vap_soa = heat of vaporization for gas soa (J/mol)
      real(r8) :: del_g_soa_tmp(ntot_soaspec)
      real(r8) :: dtcur                                ! current time step (s)
      real(r8) :: dtmax                                ! = (dtfull-tcur)
      real(r8) :: g0_soa(ntot_soaspec)                 ! ambient soa gas equilib mixrat (mol/mol at actual mw)
      real(r8) :: g_soa(ntot_soaspec)                  ! soa gas mixrat (mol/mol at actual mw)
      real(r8) :: g_star(nbins,ntot_soaspec)    ! soa gas mixrat that is in equilib
                                                       ! with each aerosol mode (mol/mol)
      real(r8) :: mw_poa                               ! actual molec wght of poa
      real(r8) :: mw_soa                               ! actual molec wght of soa
      real(r8) :: opoa_frac(ntot_poaspec)              ! fraction of poa that is opoa
      real(r8) :: phi(nbins,ntot_soaspec)       ! "relative driving force"
      real(r8) :: p0_soa(ntot_soaspec)                 ! soa gas equilib vapor presssure (atm)
      real(r8) :: p0_soa_298(ntot_soaspec)             ! p0_soa_298 = soa gas equilib vapor presssure (atm) at 298 k
      real(r8) :: sat(nbins,ntot_soaspec)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
                                                       !    used by the numerical integration scheme -- it is not a saturation rato!
      real(r8) :: tcur                                 ! current integration time (from 0 s)
      real(r8) :: tmpa, tmpb, tmpf
      real(r8) :: tot_soa(ntot_soaspec)                ! g_soa + sum( a_soa(:) )
      real(r8) :: xferrate(nbins,ntot_soaspec)    ! gas-aerosol mass transfer rate (1/s)

! Changed by Manish Shrivastava
      opoa_frac(:) = 0.0_r8 !POA does not form solution with SOA for all runs; set opoa_frac=0.0_r8  by Manish Shrivastava
      mw_poa = 250.0_r8
      mw_soa = 250.0_r8


      !write(*,*) 'ntot_soaspec ',ntot_soaspec
      ! New SOA properties added by Manish Shrivastava on 09/27/2012
      if (ntot_soaspec ==1) then
         p0_soa_298(:) = 1.0e-10_r8
         delh_vap_soa(:) = 156.0e3_r8
         opoa_frac(:) = 0.1_r8
      elseif (ntot_soaspec ==2) then
         ! same for anthropogenic and biomass burning species
         p0_soa_298 (1) = 1.0e-10_r8
         p0_soa_298 (2) = 1.0e-10_r8
         delh_vap_soa(:) = 156.0e3_r8
      elseif(ntot_soaspec ==5) then
         ! 5 volatility bins for each of the a combined SOA classes ( including biomass burning, fossil fuel, biogenic)
         p0_soa_298 (1) = 9.7831E-13_r8 !soaff0 C*=0.01ug/m3
         p0_soa_298 (2) = 9.7831E-12_r8 !soaff1 C*=0.10ug/m3
         p0_soa_298 (3) = 9.7831E-11_r8 !soaff2 C*=1.0ug/m3
         p0_soa_298 (4) = 9.7831E-10_r8 !soaff3 C*=10.0ug/m3
         p0_soa_298 (5) = 9.7831E-9_r8  !soaff4 C*=100.0ug/m3

         delh_vap_soa(1) = 153.0e3_r8
         delh_vap_soa(2) = 142.0e3_r8
         delh_vap_soa(3) = 131.0e3_r8
         delh_vap_soa(4) = 120.0e3_r8
         delh_vap_soa(5) = 109.0e3_r8
      elseif(ntot_soaspec ==15) then
         !
         ! 5 volatility bins for each of the 3 SOA classes ( biomass burning, fossil fuel, biogenic)
         ! SOA species 1-5 are for anthropogenic while 6-10 are for biomass burning SOA
         ! SOA species 11-15 are for biogenic SOA, based on Cappa et al., Reference needs to be updated
         ! For MW=250.0
         p0_soa_298 (1) = 9.7831E-13_r8 !soaff0 C*=0.01ug/m3
         p0_soa_298 (2) = 9.7831E-12_r8 !soaff1 C*=0.10ug/m3
         p0_soa_298 (3) = 9.7831E-11_r8 !soaff2 C*=1.0ug/m3
         p0_soa_298 (4) = 9.7831E-10_r8 !soaff3 C*=10.0ug/m3
         p0_soa_298 (5) = 9.7831E-9_r8  !soaff4 C*=100.0ug/m3
         p0_soa_298 (6) = 9.7831E-13_r8 !soabb0 C*=0.01ug/m3
         p0_soa_298 (7) = 9.7831E-12_r8 !soabb1 C*=0.10ug/m3
         p0_soa_298 (8) = 9.7831E-11_r8 !soabb2 C*=1.0ug/m3
         p0_soa_298 (9) = 9.7831E-10_r8 !soabb3 C*=10.0ug/m3
         p0_soa_298 (10) = 9.7831E-9_r8  !soabb4 C*=100.0ug/m3
         p0_soa_298 (11) = 9.7831E-13_r8 !soabg0 C*=0.01ug/m3
         p0_soa_298 (12) = 9.7831E-12_r8 !soabg1 C*=0.1ug/m3
         p0_soa_298 (13) = 9.7831E-11_r8 !soabg2 C*=1.0ug/m3
         p0_soa_298 (14) = 9.7831E-10_r8 !soabg3 C*=10.0ug/m3
         p0_soa_298 (15) = 9.7831E-9_r8  !soabg4 C*=100.0ug/m3

         !
         ! have to be adjusted to 15 species, following the numbers by Epstein et al., 2012
         !
         delh_vap_soa(1) = 153.0e3_r8
         delh_vap_soa(2) = 142.0e3_r8
         delh_vap_soa(3) = 131.0e3_r8
         delh_vap_soa(4) = 120.0e3_r8
         delh_vap_soa(5) = 109.0e3_r8
         delh_vap_soa(6) = 153.0e3_r8
         delh_vap_soa(7) = 142.0e3_r8
         delh_vap_soa(8) = 131.0e3_r8
         delh_vap_soa(9) = 120.0e3_r8
         delh_vap_soa(10) = 109.0e3_r8
         delh_vap_soa(11) = 153.0e3_r8
         delh_vap_soa(12) = 142.0e3_r8
         delh_vap_soa(13) = 131.0e3_r8
         delh_vap_soa(14) = 120.0e3_r8
         delh_vap_soa(15) = 109.0e3_r8
      endif

      !BSINGH - Initialized g_soa_tend and a_soa_tend to circumvent the undefined behavior (04/16/12)
      g_soa_tend(:)   = 0.0_r8
      a_soa_tend(:,:) = 0.0_r8
      xferrate(:,:) = 0.0_r8

      ! determine which modes have non-zero transfer rates
      !    and are involved in the soa gas-aerosol transfer
      ! for diameter = 1 nm and number = 1 #/cm3, xferrate ~= 1e-9 s-1
      do m = 1, nbins
         if (do_soag_any(m)) then
           skip_soamode(m) = .false.
           do ll = 1, ntot_soaspec
             xferrate(m,ll) = xferrate_in(m,ll)
           end do
           if (maxval(xferrate(m,:)) .gt. 0.0_r8) then
           !write(*,*) 'xferrate(m,ll)', m, xferrate(m,:)
           end if
         else
           skip_soamode(m) = .true.
         end if
      end do
       !write(*,*) 'xferrate(m,ll)', xferrate

      ! convert incoming mixing ratios from mol/mol at the "host-code" molec. weight (12.0 in cam5)
      !    to mol/mol at the "actual" molec. weight (currently assumed to be 250.0)
      ! also
      !    force things to be non-negative
      !    calc tot_soa(ll)
      !    calc a_opoa (always slightly >0)
      do ll = 1, ntot_soaspec
         tmpf = mw_soa_host/mw_soa
         g_soa(ll) = max( g_soa_in(ll), 0.0_r8 ) * tmpf
         tot_soa(ll) = g_soa(ll)
         do m = 1, nbins
            if ( skip_soamode(m) ) cycle
            a_soa(m,ll) = max( a_soa_in(m,ll), 0.0_r8 ) * tmpf
            tot_soa(ll) = tot_soa(ll) + a_soa(m,ll)
         end do
      end do


      tmpf = mw_poa_host/mw_poa
      do m = 1, nbins
         if ( skip_soamode(m) ) cycle
         a_opoa(m) = 0.0_r8
        !check since it seems like in the modal approach there is a bug, not summing up the values for each specie
         do ll = 1, ntot_poaspec
            tmpf = mw_poa_host/mw_poa
            a_opoa(m) = a_opoa(m) + opoa_frac(ll)*a_poa_in(m,ll)
            a_opoa(m) = max( a_opoa(m), 1.0e-40_r8 )  ! force to small non-zero value
         end do
      end do
      !write(*,*) 'a_opoa ',a_opoa

      ! calc ambient equilibrium soa gas
      do ll = 1, ntot_soaspec
         p0_soa(ll) = p0_soa_298(ll) * &
              exp( -(delh_vap_soa(ll)/rgas)*((1.0_r8/temp)-(1.0_r8/298.0_r8)) )
         g0_soa(ll) = 1.01325e5_r8*p0_soa(ll)/pres
      !   write(*,*) 'p0_soa ',ll, p0_soa(ll)
      !   write(*,*) 'g0_soa ',ll,  g0_soa(ll)
      !   write(*,*) 'p0_soa_298 ',ll, p0_soa_298(ll)
      !   write(*,*) 'delh_vap_soa ',ll, delh_vap_soa(ll)
      !   write(*,*) 'rgas ',rgas
      !   write(*,*) 'temp',temp
      end do
      !write(*,*) 'g0_soa ',minval(g0_soa),maxval(g0_soa)
      ! IF mw of soa EQ 12 (as in the MAM3 default case), this has to be in
      ! should actully talk the mw from the chemistry mechanism and substitute with 12.0

      niter = 0
      tcur = 0.0_r8
      dtcur = 0.0_r8
      phi(:,:) = 0.0_r8
      g_star(:,:) = 0.0_r8

!     if (idiagss > 0) then
!        write(luna,'(a,1p,10e11.3)') 'p0, g0_soa', p0_soa, g0_soa
!        write(luna,'(3a)') &
!           'niter, tcur,   dtcur,    phi(:),                       ', &
!           'g_star(:),                    ', &
!           'a_soa(:),                     g_soa'
!        write(luna,'(3a)') &
!           '                         sat(:),                       ', &
!           'sat(:)*a_soa(:)               ', &
!           'a_opoa(:)'
!        write(luna,'(i3,1p,20e10.2)') niter, tcur, dtcur, &
!           phi(:), g_star(:), a_soa(:), g_soa
!      end if


! integration loop -- does multiple substeps to reach dtfull
time_loop: &
      do while (tcur < dtfull-1.0e-3_r8 )

      niter = niter + 1
      if (niter > niter_max) exit

      tmpa = 0.0_r8  ! time integration parameter for all soa species
      do m = 1, nbins
         if ( skip_soamode(m) ) cycle
         a_ooa_sum_tmp(m) = a_opoa(m) + sum( a_soa(m,1:ntot_soaspec) )
      end do
      do ll = 1, ntot_soaspec
         tmpb = 0.0_r8  ! time integration parameter for a single soa species
         do m = 1, nbins
            if ( skip_soamode(m) ) cycle
            sat(m,ll) = g0_soa(ll)/max( a_ooa_sum_tmp(m), a_min1 )
            g_star(m,ll) = sat(m,ll)*a_soa(m,ll)
            phi(m,ll) = (g_soa(ll) - g_star(m,ll))/max( g_soa(ll), g_star(m,ll), g_min1 )
            tmpb = tmpb + xferrate(m,ll)*abs(phi(m,ll))
         end do
         tmpa = max( tmpa, tmpb )
      end do
      !write(*,*) 'g_star ',minval(g_star),maxval(g_star)

      if (dtsub_fixed > 0.0_r8) then
         dtcur = dtsub_fixed
         tcur = tcur + dtcur
      else
         dtmax = dtfull-tcur
         if (dtmax*tmpa <= alpha) then
! here alpha/tmpa >= dtmax, so this is final substep
            dtcur = dtmax
            tcur = dtfull
         else
            dtcur = alpha/tmpa
            tcur = tcur + dtcur
         end if
      end if

! step 1 - for modes where soa is condensing, estimate "new" a_soa(m,ll)
!    using an explicit calculation with "old" g_soa
!    and g_star(m,ll) calculated using "old" a_soa(m,ll)
! do this to get better estimate of "new" a_soa(m,ll) and sat(m,ll)
      do m = 1, nbins
         if ( skip_soamode(m) ) cycle
         do ll = 1, ntot_soaspec
            ! first ll loop calcs a_soa_tmp(m,ll) & a_ooa_sum_tmp
            a_soa_tmp(m,ll) = a_soa(m,ll)
            beta(m,ll) = dtcur*xferrate(m,ll)
            del_g_soa_tmp(ll) = g_soa(ll) - g_star(m,ll)
            if (del_g_soa_tmp(ll) > 0.0_r8) then
               a_soa_tmp(m,ll) = a_soa(m,ll) + beta(m,ll)*del_g_soa_tmp(ll)
            end if
         end do
         a_ooa_sum_tmp(m) = a_opoa(m) + sum( a_soa_tmp(m,1:ntot_soaspec) )
         do ll = 1, ntot_soaspec
            ! second ll loop calcs sat & g_star
            if (del_g_soa_tmp(ll) > 0.0_r8) then
               sat(m,ll) = g0_soa(ll)/max( a_ooa_sum_tmp(m), a_min1 )
               g_star(m,ll) = sat(m,ll)*a_soa_tmp(m,ll)   ! this just needed for diagnostics
            end if
         end do
      end do
      !write(*,*) 'g_star again ',minval(g_star), maxval(g_star)

! step 2 - implicit in g_soa and semi-implicit in a_soa,
!    with g_star(m,ll) calculated semi-implicitly
      do ll = 1, ntot_soaspec
         tmpa = 0.0_r8
         tmpb = 0.0_r8
         do m = 1, nbins
            if ( skip_soamode(m) ) cycle
            tmpa = tmpa + a_soa(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
            tmpb = tmpb + beta(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
         end do

         g_soa(ll) = (tot_soa(ll) - tmpa)/(1.0_r8 + tmpb)
         !write(*,*) 'g_soa  ',g_soa(ll)
         g_soa(ll) = max( 0.0_r8, g_soa(ll) )
         do m = 1, nbins
            if ( skip_soamode(m) ) cycle
            a_soa(m,ll) = (a_soa(m,ll) + beta(m,ll)*g_soa(ll))/   &
                       (1.0_r8 + beta(m,ll)*sat(m,ll))
         end do
      end do
      !write(*,*) 'a_soa, g_soa ',minval(a_soa),maxval(a_soa), g_soa

!     if (idiagss > 0) then
!        write(luna,'(i3,1p,20e10.2)') niter, tcur, dtcur, &
!           phi(:), g_star(:), a_soa(:,), g_soa
!        write(luna,'(23x,1p,20e10.2)') &
!           sat(:), sat(:)*a_soa(:), a_opoa(:)
!    end if

!     if (niter > 9992000) then
!        write(iulog,*) '*** to many iterations'
!         exit
!     end if

      end do time_loop
      !write(*,*) 'a_soa, g_soa time_loop ',minval(g_star), maxval(g_star), minval(a_soa),maxval(a_soa), g_soa


! calculate outgoing tendencies (at the host-code molec. weight)
! (a_soa & g_soa are at actual mw, but a_soa_in & g_soa_in are at host-code mw)
      do ll = 1, ntot_soaspec
         tmpf = mw_soa/mw_soa_host
         !write(*,*) 'tmf',tmf
         !write(*,*) 'g_soa(ll)', g_soa(ll)
         !write(*,*) 'g_soa_in(ll)', g_soa_in(ll)
         g_soa_tend(ll) = (g_soa(ll)*tmpf - g_soa_in(ll))/dtfull
         !write(*,*) 'g_soa_tend(ll)', g_soa_tend(ll)
         do m = 1, nbins
            if ( skip_soamode(m) ) cycle
            a_soa_tend(m,ll) = (a_soa(m,ll)*tmpf - a_soa_in(m,ll))/dtfull
         !   write(*,*) 'a_soa(ll)', a_soa(m,ll)
         !   write(*,*) 'a_soa_in(ll)', a_soa_in(m,ll)
         !   write(*,*) 'a_soa_tend(ll)', a_soa_tend(ll)
         end do
      end do
      !write(*,*) 'a_soa_tend ',minval(a_soa_tend), maxval(a_soa_tend)


      return

      end subroutine carma_aero_soaexch

!----------------------------------------------------------------------


end module carma_aero_gasaerexch
