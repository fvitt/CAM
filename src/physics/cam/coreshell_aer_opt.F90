module coreshell_aer_opt

! Inspired by model_aer_opt and modifications for CARMA from Pengfei Yu,
! this module implements a core/shell representation for aerosol optical
! properties. This reuses the mode definition structure in the namelist to
! allow a mode to be defined as a bin with two (or more) components: a
! shell and one or more cores, which are each separate tracers that are
! internally mixed.
!
! NOTE: The initial implementation is for a smoke particle that has a
! black carbon core and an organic shell that do not take up water. In
! the future this could be generalized to multiple components that are
! hydrophilic. This would be done by mixing the refractive indices for
! components and the shell components as well as determining the
! hydropscopicity of the particles and adding additional dimension for
! core and shell refractive indices to the optical properties lookup
! table.
!
! Chuck Bardeen
! June 2019


use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
use ppgrid,            only: pcols, pver, pverp
use constituents,      only: pcnst
use spmd_utils,        only: masterproc
use physconst,         only: rhoh2o, rga, rair, pi
use radconstants,      only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
use radconstants,      only: ot_length, get_lw_spectral_boundaries
use rad_constituents,  only: n_diag, rad_cnst_get_call_list, rad_cnst_get_info, rad_cnst_get_info_by_bin, &
                             rad_cnst_get_bin_mmr_by_idx, rad_cnst_get_bin_props_by_idx, &
                             rad_cnst_get_bin_props
use physics_types,     only: physics_state
use physics_buffer, only : pbuf_get_index,physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

use pio,               only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
                             pio_get_var, pio_nowrite, pio_closefile
use cam_pio_utils,     only: cam_pio_openfile
use cam_history,       only: addfld, add_default, outfld, horiz_only
use cam_history_support, only: fillvalue
use cam_logfile,       only: iulog
use perf_mod,          only: t_startf, t_stopf
use cam_abortutils,    only: endrun
use wv_saturation,    only: qsat

use table_interp_mod, only: table_interp
implicit none
private
save

public :: coreshell_aer_opt_readnl, coreshell_aer_opt_init, coreshell_aero_sw, coreshell_aero_lw

! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
! in terms of refractive index and wet radius

character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                       '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

integer :: lw10um_indx = -1

! refractive index for water read in read_water_refindex
complex(r8) :: crefwsw(nswbands) = -huge(1._r8) ! complex refractive index for water visible
complex(r8) :: crefwlw(nlwbands) = -huge(1._r8) ! complex refractive index for water infrared

character(shr_kind_cl) :: water_refindex_file = 'none' ! full pathname for water refractive index dataset

!===============================================================================
CONTAINS
!===============================================================================

subroutine coreshell_aer_opt_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'modal_aer_opt_readnl'

   namelist /modal_aer_opt_nl/ water_refindex_file
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'modal_aer_opt_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, modal_aer_opt_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   call mpibcast(water_refindex_file, len(water_refindex_file), mpichar, 0, mpicom)
#endif

    if (masterproc) then
       write(iulog,*) subname,': water_refindex_file = ',trim(water_refindex_file)
    end if


 end subroutine coreshell_aer_opt_readnl

!===============================================================================

subroutine coreshell_aer_opt_init()

   use ioFileMod,        only: getfil
   use phys_control,     only: phys_getopts

   ! Local variables

   integer  :: i, m

   logical :: call_list(0:n_diag)
   integer :: ilist, nbins, m_prefr, m_prefi
   character(len=12) :: fldname
   character(len=128) :: lngname

   character(len=*), parameter :: routine='coreshell_aer_opt_init'
   real(r8) :: lwavlen_lo(nlwbands), lwavlen_hi(nlwbands)

   logical           :: history_amwg            ! output the variables used by the AMWG diag package
   logical           :: history_aero_optics     ! output aerosol optics diagnostics
   logical           :: history_dust            ! output dust diagnostics

   character(len=256) :: locfile

   !----------------------------------------------------------------------------

   ! Check that dimension sizes in the coefficient arrays used to
   ! parameterize aerosol radiative properties are consistent between this
   ! module and the aerosol physprop files.
   call rad_cnst_get_call_list(call_list)

   call rad_cnst_get_info(0, nbins=nbins)

   call phys_getopts(history_amwg_out        = history_amwg, &
                     history_aero_optics_out = history_aero_optics, &
                     history_dust_out        = history_dust )


   ! Add diagnostic fields to history output.
   !
   ! NOTE: If modal aerosols are also being used, then these names will conflict.
   ! Probably need to define unique field names for the sectional aerosols equivalent
   ! fields.

   ! Add diagnostic fields to history output.

   call addfld ('EXTINCT',    (/ 'lev' /), 'A','/m','Aerosol extinction 550 nm, day only',                   &
                flag_xyfill=.true.)
   call addfld ('EXTINCTUV',  (/ 'lev' /), 'A','/m','Aerosol extinction 350 nm, day only',                   &
                flag_xyfill=.true.)
   call addfld ('EXTINCTNIR', (/ 'lev' /), 'A','/m','Aerosol extinction 1020 nm, day only',                  &
                flag_xyfill=.true.)
   call addfld ('ABSORB',     (/ 'lev' /), 'A','/m','Aerosol absorption, day only',                          &
                flag_xyfill=.true.)
   call addfld ('AODTOT', horiz_only, 'A','1','Aerosol optical depth summed over all sw wavelenghts', flag_xyfill=.true.)
   call addfld ('AODTOTdn', horiz_only, 'A','1','Aerosol optical depth summed over all sw wavelenghts')
   call addfld ('AODVIS',     horiz_only,  'A','  ','Aerosol optical depth 550 nm, day only',                &
                flag_xyfill=.true.)
   call addfld ('AODVISst',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 550 nm, day only',  &
                flag_xyfill=.true.)
   call addfld ('AODUV',      horiz_only,  'A','  ','Aerosol optical depth 350 nm, day only',                &
                flag_xyfill=.true.)
   call addfld ('AODUVst',    horiz_only,  'A','  ','Stratospheric aerosol optical depth 350 nm, day only',  &
                flag_xyfill=.true.)
   call addfld ('AODNIR',     horiz_only,  'A','  ','Aerosol optical depth 1020 nm, day only',               &
                flag_xyfill=.true.)
   call addfld ('AODNIRst',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 1020 nm, day only', &
                flag_xyfill=.true.)
   call addfld ('AODABS',     horiz_only,  'A','  ','Aerosol absorption optical depth 550 nm, day only',     &
                flag_xyfill=.true.)
   call addfld ('AODxASYM',   horiz_only,  'A','  ','Aerosol optical depth 550 * asymmetry factor, day only',&
                flag_xyfill=.true.)
   call addfld ('EXTxASYM',   (/ 'lev' /), 'A','  ','extinction 550 nm * asymmetry factor, day only',        &
                flag_xyfill=.true.)

   call addfld ('EXTINCTdn',    (/ 'lev' /), 'A','/m','Aerosol extinction 550 nm, day night',                 &
                  flag_xyfill=.true.)
   call addfld ('EXTINCTUVdn',  (/ 'lev' /), 'A','/m','Aerosol extinction 350 nm, day night',                 &
                  flag_xyfill=.true.)
   call addfld ('EXTINCTNIRdn', (/ 'lev' /), 'A','/m','Aerosol extinction 1020 nm, day night',                &
                  flag_xyfill=.true.)
   call addfld ('ABSORBdn',     (/ 'lev' /), 'A','/m','Aerosol absorption, day night',                        &
                  flag_xyfill=.true.)
   call addfld ('AODVISdn',     horiz_only,  'A','  ','Aerosol optical depth 550 nm, day night',              &
                  flag_xyfill=.true.)
   call addfld ('AODVISstdn',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 550 nm, day night',&
                  flag_xyfill=.true.)
   call addfld ('AODUVdn',      horiz_only,  'A','  ','Aerosol optical depth 350 nm, day night',              &
                  flag_xyfill=.true.)
   call addfld ('AODUVstdn',    horiz_only,  'A','  ','Stratospheric aerosol optical depth 350 nm, day night',&
                  flag_xyfill=.true.)
   call addfld ('AODNIRdn',     horiz_only,  'A','  ','Aerosol optical depth 1020 nm, day night',             &
                  flag_xyfill=.true.)
   call addfld ('AODNIRstdn',   horiz_only,  'A','  ','Stratospheric aerosol optical depth 1020 nm, day night',&
                 flag_xyfill=.true.)
   call addfld ('AODABSdn',     horiz_only,  'A','  ','Aerosol absorption optical depth 550 nm, day night',   &
                  flag_xyfill=.true.)
   call addfld ('AODxASYMdn',   horiz_only,  'A','  ','Aerosol optical depth 550 * asymmetry factor, day night',&
                flag_xyfill=.true.)
   call addfld ('EXTxASYMdn',   (/ 'lev' /), 'A','  ','extinction 550 * asymmetry factor, day night',        &
                flag_xyfill=.true.)



   do m = 1, nbins
      write(fldname,'(a,i2.2)') 'BURDEN', m
      write(lngname,'(a,i2.2)') 'Aerosol burden bin ', m
      call addfld (fldname, horiz_only, 'A', 'kg/m2', lngname, flag_xyfill=.true.)
      if (history_aero_optics) then
         call add_default (fldname, 1, ' ')
      end if

      write(fldname,'(a,i2.2)') 'AODBIN', m
      write(lngname,'(a,i2)') 'Aerosol optical depth, day only, 550 nm bin ', m
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (history_aero_optics) then
         call add_default (fldname, 1, ' ')
      end if

      write(fldname,'(a,i2.2)') 'AODDUST', m
      write(lngname,'(a,i2,a)') 'Aerosol optical depth, day only, 550 nm mode ',m,' from dust'
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (history_aero_optics) then
         call add_default (fldname, 1, ' ')
      end if

      write(fldname,'(a,i2.2)') 'BURDENdn', m
      write(lngname,'(a,i2)') 'Aerosol burden, day night, bin ', m
      call addfld (fldname, horiz_only, 'A', 'kg/m2', lngname, flag_xyfill=.true.)
      if (history_aero_optics) then
         call add_default (fldname, 1, ' ')
      end if

      write(fldname,'(a,i2.2)') 'AODdnBIN', m
      write(lngname,'(a,i2)') 'Aerosol optical depth 550 nm, day night, bin ', m
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (history_aero_optics) then
         call add_default (fldname, 1, ' ')
      end if

      write(fldname,'(a,i2.2)') 'AODdnDUST', m
      write(lngname,'(a,i2,a)') 'Aerosol optical depth 550 nm, day night, bin ',m,' from dust'
      call addfld (fldname, horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
      if (history_aero_optics) then
         call add_default (fldname, 1, ' ')
      end if

   enddo

   call addfld ('AODDUST',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODSO4',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SO4, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODPOM',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from POM, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODSOA',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SOA, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODBC',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from BC, day only',           &
        flag_xyfill=.true.)
   call addfld ('AODSS',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from seasalt, day only',      &
        flag_xyfill=.true.)
   call addfld ('AODABSBC',      horiz_only, 'A','  ',    'Aerosol absorption optical depth 550 nm from BC, day only',&
        flag_xyfill=.true.)
   call addfld ('BURDENDUST',    horiz_only, 'A','kg/m2', 'Dust aerosol burden, day only'        ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSO4',     horiz_only, 'A','kg/m2', 'Sulfate aerosol burden, day only'     ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENPOM',     horiz_only, 'A','kg/m2', 'POM aerosol burden, day only'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSOA',     horiz_only, 'A','kg/m2', 'SOA aerosol burden, day only'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENBC',      horiz_only, 'A','kg/m2', 'Black carbon aerosol burden, day only',                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSEASALT', horiz_only, 'A','kg/m2', 'Seasalt aerosol burden, day only'     ,                    &
        flag_xyfill=.true.)
   call addfld ('SSAVIS',        horiz_only, 'A','  ',    'Aerosol single-scatter albedo, day only',                  &
        flag_xyfill=.true.)

   call addfld ('AODDUSTdn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODSO4dn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SO4, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODPOMdn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from POM, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODSOAdn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SOA, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODBCdn',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from BC, day night',           &
        flag_xyfill=.true.)
   call addfld ('AODSSdn',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from seasalt, day night',      &
        flag_xyfill=.true.)
   call addfld ('AODABSBCdn',      horiz_only, 'A','  ',    'Aerosol absorption optical depth 550 nm from BC, day night',&
        flag_xyfill=.true.)
   call addfld ('BURDENDUSTdn',    horiz_only, 'A','kg/m2', 'Dust aerosol burden, day night'        ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSO4dn',     horiz_only, 'A','kg/m2', 'Sulfate aerosol burden, day night'     ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENPOMdn',     horiz_only, 'A','kg/m2', 'POM aerosol burden, day night'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSOAdn',     horiz_only, 'A','kg/m2', 'SOA aerosol burden, day night'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENBCdn',      horiz_only, 'A','kg/m2', 'Black carbon aerosol burden, day night',                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSEASALTdn', horiz_only, 'A','kg/m2', 'Seasalt aerosol burden, day night'     ,                    &
        flag_xyfill=.true.)
   call addfld ('SSAVISdn',        horiz_only, 'A','  ',    'Aerosol single-scatter albedo, day night',                  &
        flag_xyfill=.true.)


   if (history_amwg) then
      call add_default ('AODDUST'      , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
   end if

   if (history_dust) then
      call add_default ('AODDUST'     , 1, ' ')
   end if

   if (history_aero_optics) then
      call add_default ('ABSORB'       , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
      call add_default ('AODUV'        , 1, ' ')
      call add_default ('AODNIR'       , 1, ' ')
      call add_default ('AODABS'       , 1, ' ')
      call add_default ('AODABSBC'     , 1, ' ')
      call add_default ('AODDUST'      , 1, ' ')
      call add_default ('AODSO4'       , 1, ' ')
      call add_default ('AODPOM'       , 1, ' ')
      call add_default ('AODSOA'       , 1, ' ')
      call add_default ('AODBC'        , 1, ' ')
      call add_default ('AODSS'        , 1, ' ')
      call add_default ('BURDENDUST'   , 1, ' ')
      call add_default ('BURDENSO4'    , 1, ' ')
      call add_default ('BURDENPOM'    , 1, ' ')
      call add_default ('BURDENSOA'    , 1, ' ')
      call add_default ('BURDENBC'     , 1, ' ')
      call add_default ('BURDENSEASALT', 1, ' ')
      call add_default ('SSAVIS'       , 1, ' ')
      call add_default ('EXTINCT'      , 1, ' ')
      call add_default ('AODxASYM'     , 1, ' ')
      call add_default ('EXTxASYM'     , 1, ' ')

      call add_default ('ABSORBdn'       , 1, ' ')
      call add_default ('AODVISdn'       , 1, ' ')
      call add_default ('AODUVdn'        , 1, ' ')
      call add_default ('AODNIRdn'       , 1, ' ')
      call add_default ('AODABSdn'       , 1, ' ')
      call add_default ('AODABSBCdn'     , 1, ' ')
      call add_default ('AODDUSTdn'      , 1, ' ')
      call add_default ('AODSO4dn'       , 1, ' ')
      call add_default ('AODPOMdn'       , 1, ' ')
      call add_default ('AODSOAdn'       , 1, ' ')
      call add_default ('AODBCdn'        , 1, ' ')
      call add_default ('AODSSdn'        , 1, ' ')
      call add_default ('BURDENDUSTdn'   , 1, ' ')
      call add_default ('BURDENSO4dn'    , 1, ' ')
      call add_default ('BURDENPOMdn'    , 1, ' ')
      call add_default ('BURDENSOAdn'    , 1, ' ')
      call add_default ('BURDENBCdn'     , 1, ' ')
      call add_default ('BURDENSEASALTdn', 1, ' ')
      call add_default ('SSAVISdn'       , 1, ' ')
      call add_default ('EXTINCTdn'      , 1, ' ')
      call add_default ('AODxASYMdn'     , 1, ' ')
      call add_default ('EXTxASYMdn'     , 1, ' ')
  end if

   do ilist = 1, n_diag
      if (call_list(ilist)) then

         call addfld ('EXTINCT'//diag(ilist),  (/ 'lev' /), 'A','/m', &
              'Aerosol extinction', flag_xyfill=.true.)
         call addfld ('ABSORB'//diag(ilist),   (/ 'lev' /), 'A','/m', &
              'Aerosol absorption', flag_xyfill=.true.)
         call addfld ('AODVIS'//diag(ilist),   horiz_only,  'A','  ', &
              'Aerosol optical depth 550 nm', flag_xyfill=.true.)
         call addfld ('AODVISst'//diag(ilist), horiz_only,  'A','  ', &
              'Stratospheric aerosol optical depth 550 nm', flag_xyfill=.true.)
         call addfld ('AODABS'//diag(ilist),   horiz_only,  'A','  ', &
              'Aerosol absorption optical depth 550 nm', flag_xyfill=.true.)

         call addfld ('AODTOT'//diag(ilist), horiz_only, 'A','1', &
              'Aerosol optical depth summed over all sw wavelenghts', flag_xyfill=.true.)
         call addfld ('AODTOTdn'//diag(ilist), horiz_only, 'A','1', &
              'Aerosol optical depth summed over all sw wavelenghts')

         call addfld ('EXTINCTdn'//diag(ilist),    (/ 'lev' /), 'A','/m',&
              'Aerosol extinction 550 nm, day night', flag_xyfill=.true.)
         call addfld ('ABSORBdn'//diag(ilist),     (/ 'lev' /), 'A','/m',&
              'Aerosol absorption, day night',  flag_xyfill=.true.)
         call addfld ('AODVISdn'//diag(ilist),     horiz_only,  'A','  ',&
              'Aerosol optical depth 550 nm, day night',  flag_xyfill=.true.)
         call addfld ('AODVISstdn'//diag(ilist),   horiz_only,  'A','  ',&
              'Stratospheric aerosol optical depth 550 nm, day night', flag_xyfill=.true.)
         call addfld ('AODABSdn'//diag(ilist),     horiz_only,  'A','  ',&
              'Aerosol absorption optical depth 550 nm, day night',  flag_xyfill=.true.)
         call addfld ('EXTxASYMdn'//diag(ilist),   (/ 'lev' /), 'A','  ',&
              'extinction 550 * asymmetry factor, day night',  flag_xyfill=.true.)
         call addfld ('EXTxASYM'//diag(ilist),   (/ 'lev' /), 'A','  ',&
              'extinction 550 nm * asymmetry factor, day only',   flag_xyfill=.true.)

         if (history_aero_optics) then
            call add_default ('EXTINCT'//diag(ilist), 1, ' ')
            call add_default ('ABSORB'//diag(ilist),  1, ' ')
            call add_default ('AODVIS'//diag(ilist),  1, ' ')
            call add_default ('AODVISst'//diag(ilist),  1, ' ')
            call add_default ('AODABS'//diag(ilist),  1, ' ')
         end if

      end if
   end do



   call get_lw_spectral_boundaries(lwavlen_lo, lwavlen_hi, units='um')
   do i = 1,nlwbands
      if ((lwavlen_lo(i)<=10._r8) .and. (lwavlen_hi(i)>=10._r8)) then
         lw10um_indx = i
      end if
   end do
   if (lw10um_indx>0) then
      call addfld('AODABSLW', (/ 'lev' /), 'A','/m','Aerosol long-wave absorption optical depth at 10 microns')
   end if
   call addfld ('TOTABSLW', (/ 'lev' /), 'A',' ', 'LW Aero total abs')

   call getfil(water_refindex_file, locfile)
   call read_water_refindex(locfile)

end subroutine coreshell_aer_opt_init

!===============================================================================

subroutine coreshell_aero_sw(list_idx, state, pbuf, nnite, idxnite, &
                         tauxar, wa, ga, fa)

   use tropopause,      only : tropopause_find

   ! calculates aerosol sw radiative properties

   integer,             intent(in) :: list_idx       ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state          ! state variables

   type(physics_buffer_desc), pointer :: pbuf(:)
   integer,             intent(in) :: nnite          ! number of night columns
   integer,             intent(in) :: idxnite(nnite) ! local column indices of night columns

   real(r8), intent(inout) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
   real(r8), intent(inout) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
   real(r8), intent(inout) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
   real(r8), intent(inout) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

   ! Local variables
   integer :: i, ifld, isw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nbins
   integer :: nspec
   integer :: troplev(pcols)

   real(r8) :: mass(pcols,pver)        ! layer mass
   real(r8) :: air_density(pcols,pver) ! (kg/m3)

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index
   character*32         :: spectype            ! species type
   real(r8)             :: hygro_aer           !

   real(r8)    :: refr(pcols)     ! real part of refractive index
   real(r8)    :: refi(pcols)     ! imaginary part of refractive index
   complex(r8) :: crefin(pcols)   ! complex refractive index
   real(r8), pointer :: tbl_corefrac(:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: extpsw(:,:) ! specific extinction
   real(r8), pointer :: abspsw(:,:) ! specific absorption
   real(r8), pointer :: asmpsw(:,:) ! asymmetry factor

   real(r8) :: vol(pcols)      ! volume concentration of aerosol specie (m3/kg)

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: pext(pcols)     ! parameterized specific extinction (m2/kg)
   real(r8) :: specpext(pcols) ! specific extinction (m2/kg)
   real(r8) :: dopaer(pcols)   ! aerosol optical depth in layer
   real(r8) :: pabs(pcols)     ! parameterized specific absorption (m2/kg)
   real(r8) :: pasm(pcols)     ! parameterized asymmetry factor
   real(r8) :: palb(pcols)     ! parameterized single scattering albedo

   real(r8) :: particlemmr(pcols)   ! mmr of aerosol particle
   real(r8) :: particlerfr(pcols)   ! real refractive index of aerosol particle
   real(r8) :: particlerfi(pcols)   ! imaginary refractive index of aerosol particle
   real(r8) :: particlehygro(pcols) ! hygroscopicity of aerosol particle
   real(r8) :: particlevol(pcols)   ! volume of aerosol particle

   real(r8) :: coredustmmr(pcols)
   real(r8) :: corebcmmr(pcols)

   real(r8) :: coremmr(pcols)   ! mmr of aerosol core
   real(r8) :: corerfr(pcols)   ! real refractive index of aerosol core
   real(r8) :: corerfi(pcols)   ! imaginary refractive index of aerosol core
   real(r8) :: corehygro(pcols) ! hygroscopicity of aerosol core
   real(r8) :: corevol(pcols)   ! volume of aerosol core

   real(r8) :: shellmmr(pcols)   ! mmr of aerosol shell
   real(r8) :: shellrfr(pcols)   ! real refractive index of aerosol shell
   real(r8) :: shellrfi(pcols)   ! imaginary refractive index of aerosol shell
   real(r8) :: shellhygro(pcols) ! hygroscopicity of aerosol shell
   real(r8) :: shellvol(pcols)   ! volume of aerosol shell

   real(r8) :: totalmmr(pcols)   ! total mmr of the aerosol
   real(r8) :: bcdustmmr(pcols)  ! total mmr of dust + bc aerosol
   real(r8) :: corefrac(pcols)   ! mass fraction that is core
   real(r8) :: bcdust(pcols)     ! mass fraction of bc vs (bc + dust)
   real(r8) :: dryvol_m3(pcols)  ! dry volume (m3)
   real(r8) :: dryvol(pcols)     ! dry volume density (m3/kg)
   real(r8) :: wet_dry_frac(pcols) ! wet/dry fraction

   ! Diagnostics
   real(r8) :: extinct(pcols,pver)
   real(r8) :: extinctnir(pcols,pver)
   real(r8) :: extinctuv(pcols,pver)
   real(r8) :: absorb(pcols,pver)
   real(r8) :: aodtot(pcols)               ! extinction optical depth
   real(r8) :: aodvis(pcols)               ! extinction optical depth
   real(r8) :: aodvisst(pcols)             ! stratospheric extinction optical depth
   real(r8) :: aodabs(pcols)               ! absorption optical depth
   real(r8) :: ssavis(pcols)
   real(r8) :: specrefr, specrefi
   real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro

   logical :: savaervis ! true if visible wavelength (0.55 micron)
   logical :: savaernir ! true if near ir wavelength (~0.88 micron)
   logical :: savaeruv  ! true if uv wavelength (~0.35 micron)

   real(r8) :: dustvol(pcols)              ! volume concentration of dust in aerosol bin (m3/kg)

   real(r8) :: burden(pcols)

   real(r8) :: burdendust(pcols), burdenso4(pcols), burdenbc(pcols), &
               burdenpom(pcols), burdensoa(pcols), burdenseasalt(pcols)

   real(r8) :: scatdust(pcols), scatso4(pcols), scatbc(pcols), &
               scatpom(pcols), scatsoa(pcols), scatseasalt(pcols)
   real(r8) :: absdust(pcols), absso4(pcols), absbc(pcols), &
               abspom(pcols), abssoa(pcols), absseasalt(pcols)
   real(r8) :: hygrodust(pcols), hygroso4(pcols), hygrobc(pcols), &
               hygropom(pcols), hygrosoa(pcols), hygroseasalt(pcols)

   real(r8) :: aodc                        ! aod of component

   ! total species AOD
   real(r8) :: dustaod(pcols), so4aod(pcols), bcaod(pcols), &
               pomaod(pcols), soaaod(pcols), seasaltaod(pcols)

   real(r8) :: aodbin(pcols)
   real(r8) :: dustaodbin(pcols)          ! dust aod in aerosol bin
   real(r8) :: aodabsbc(pcols)             ! absorption optical depth of BC
   real(r8) :: asymvis(pcols)              ! asymmetry factor * optical depth
   real(r8) :: asymext(pcols,pver)         ! asymmetry factor * extinction

   real(r8) :: wetvol(pcols)     ! wet volume density (m3/kg)
   real(r8) :: wetvol_m3(pcols)  ! wet volume (m3)
   real(r8) :: watervol(pcols)
   real(r8), pointer :: wetr(:,:)
   real(r8), pointer :: dryr(:,:)

   real(r8) :: aoduv(pcols)               ! extinction optical depth in uv
   real(r8) :: aoduvst(pcols)             ! stratospheric extinction optical depth in uv
   real(r8) :: aodnir(pcols)              ! extinction optical depth in nir
   real(r8) :: aodnirst(pcols)            ! stratospheric extinction optical depth in nir

   real(r8) :: sate(pcols,pver)     ! saturation vapor pressure
   real(r8) :: satq(pcols,pver)     ! saturation specific humidity
   real(r8) :: relh(pcols,pver)

   real(r8), pointer, dimension(:,:) :: crkappa   ! kappa
   real(r8), pointer, dimension(:,:) :: wgtpct  ! weight percent
   integer      :: bin
   real(r8), pointer :: h_ext_coreshell(:,:,:,:,:)
   real(r8), pointer :: h_ssa_coreshell(:,:,:,:,:)
   real(r8), pointer :: h_asm_coreshell(:,:,:,:,:)
   real(r8), pointer :: h_ext_wtp(:,:) ! specific extinction
   real(r8), pointer :: h_ssa_wtp(:,:) ! specific absorption
   real(r8), pointer :: h_asm_wtp(:,:) ! asymmetry factor
   real(r8), pointer, dimension(:,:)   :: cldn

   real(r8), pointer :: tbl_wgtpct(:)
   real(r8), pointer :: tbl_bcdust(:)
   real(r8), pointer :: tbl_kap(:)
   real(r8), pointer :: tbl_relh(:)
   integer :: nwtp
   integer :: nbcdust
   integer :: nkap
   integer :: nrelh

   character(len=ot_length) :: opticstype
   character(len=32) :: bin_name
   character(len=32) :: specmorph

   character(len=32) :: outname

   ! debug output
   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf            ! volume fraction of insoluble aerosol
   character(len=*), parameter :: subname = 'coreshell_aero_sw'

   real(r8), parameter :: four_thirds_pi = pi * 4._r8 / 3._r8

   !----------------------------------------------------------------------------

   nullify(h_ext_coreshell)
   nullify(h_ssa_coreshell)
   nullify(h_asm_coreshell)
   nullify(h_ext_wtp)
   nullify(h_ssa_wtp)
   nullify(h_asm_wtp)

   nullify(tbl_wgtpct)
   nullify(tbl_bcdust)
   nullify(tbl_kap)
   nullify(tbl_relh)

   lchnk = state%lchnk
   ncol  = state%ncol

   mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
   air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

   ! diagnostics for visible band summed over bins
   burdendust(:ncol)     = 0.0_r8
   burdenso4(:ncol)      = 0.0_r8
   burdenpom(:ncol)      = 0.0_r8
   burdensoa(:ncol)      = 0.0_r8
   burdenbc(:ncol)       = 0.0_r8
   burdenseasalt(:ncol)  = 0.0_r8

   extinct(1:ncol,:)     = 0.0_r8
   absorb(1:ncol,:)      = 0.0_r8
   aodtot(1:ncol)        = 0.0_r8
   aodvis(1:ncol)        = 0.0_r8
   aodvisst(1:ncol)      = 0.0_r8
   aodabs(1:ncol)        = 0.0_r8
   aodbin(1:ncol)        = 0.0_r8
   ssavis(1:ncol)        = 0.0_r8
   asymvis(1:ncol)       = 0.0_r8
   asymext(1:ncol,:)     = 0.0_r8

   ! diags for other bands
   extinctuv(1:ncol,:)   = 0.0_r8
   extinctnir(1:ncol,:)  = 0.0_r8
   aoduv(:ncol)          = 0.0_r8
   aodnir(:ncol)         = 0.0_r8
   aoduvst(:ncol)        = 0.0_r8
   aodnirst(:ncol)       = 0.0_r8

   aodabsbc(:ncol)       = 0.0_r8
   dustaod(:ncol)        = 0.0_r8
   so4aod(:ncol)         = 0.0_r8
   pomaod(:ncol)         = 0.0_r8
   soaaod(:ncol)         = 0.0_r8
   bcaod(:ncol)          = 0.0_r8
   seasaltaod(:ncol)     = 0.0_r8

   call tropopause_find(state, troplev)

!st
   ! calculate relative humidity for table lookup into rh grid
   call qsat(state%t(1:ncol,1:pver), state%pmid(1:ncol,1:pver), sate(1:ncol,1:pver), satq(1:ncol,1:pver), ncol,pver)

   relh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / satq(1:ncol,1:pver)
   !There is a  different way to diagnose sub-grid rh for clear-sky only that makes physically
   ! more sense but is not used in MAM. This needs to be tested if it makes any differente
   ! Associate pointers with physics buffer fields
   ! itim = pbuf_old_tim_idx()
   ! call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cldn, (/1,1,itim/),(/pcols,pver,1/))
   ! rh(1:ncol,1:pver) = (state%q(1:ncol,1:pver,1)-cldn(1:ncol,1:pver)*qsat(1:ncol,1:pver)) / qsat(1:ncol,1:pver)
   relh(1:ncol,1:pver) = max(1.e-20_r8,relh(1:ncol,1:pver))


   ! loop over all aerosol bins
   call rad_cnst_get_info(list_idx, nbins=nbins)

   do m = 1, nbins

      call rad_cnst_get_info_by_bin(list_idx, m, bin_name=bin_name)
      call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_wetr"),wetr)
      call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_dryr"),dryr)

      burden(:) = 0._r8
      aodbin(:) = 0.0_r8
      dustaodbin(:) = 0.0_r8


     ! get bin info
     call rad_cnst_get_info_by_bin(list_idx, m, nspec=nspec, bin_name=bin_name)

     ! get optics type
     call rad_cnst_get_bin_props(list_idx, m, opticstype=opticstype)

     select case (trim(opticstype))
     case('sectional')
     ! get bin properties
        call rad_cnst_get_bin_props(list_idx, m, &
          corefrac=tbl_corefrac, &
          extpsw=extpsw, abspsw=abspsw, asmpsw=asmpsw)
     case('hygroscopic_coreshell')
     ! get optical properties for hygroscopic aerosols
        call rad_cnst_get_bin_props(list_idx, m, &
             corefrac=tbl_corefrac, &
             kap=tbl_kap, nkap=nkap, bcdust=tbl_bcdust, nbcdust=nbcdust, &
             relh=tbl_relh, nrelh=nrelh, &
             sw_hygro_coreshell_ext=h_ext_coreshell, &
             sw_hygro_coreshell_ssa=h_ssa_coreshell, sw_hygro_coreshell_asm=h_asm_coreshell)
        ! get kappa
        ! need to read in bin_name
        call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_kappa"),crkappa)
        where ( crkappa(:ncol,:) < minval(tbl_kap) )
           crkappa(:ncol,:) = minval(tbl_kap)
        end where
        where ( crkappa(:ncol,:) > maxval(tbl_kap) )
           crkappa(:ncol,:) = maxval(tbl_kap)
        end where
        where ( relh(:ncol,:) < minval(tbl_relh) )
           relh(:ncol,:) = minval(tbl_relh)
        end where
        where ( relh(:ncol,:) > maxval(tbl_relh) )
           relh(:ncol,:) = maxval(tbl_relh)
        end where
     case('hygroscopic_wtp')
      ! get optical properties for hygroscopic aerosols
         call rad_cnst_get_bin_props(list_idx, m, wgtpct=tbl_wgtpct,nwtp=nwtp, &
           sw_hygro_ext_wtp=h_ext_wtp, sw_hygro_ssa_wtp=h_ssa_wtp, sw_hygro_asm_wtp=h_asm_wtp)
       ! determine weight precent of H2SO4/H2O solution
         call pbuf_get_field(pbuf, pbuf_get_index('WTP'),wgtpct)
         where ( wgtpct(:ncol,:) < minval(tbl_wgtpct) )
            wgtpct(:ncol,:) = minval(tbl_wgtpct)
         end where
         where ( wgtpct(:ncol,:) > maxval(tbl_wgtpct) )
            wgtpct(:ncol,:) = maxval(tbl_wgtpct)
         end where
     case('zero')
          ! zero aerosols types have no optical effect, so do nothing.
     case default
          call endrun('coreshell_aero_sw: unsupported opticstype: '//trim(opticstype))
     end select

     do isw = 1, nswbands
        savaervis = (isw .eq. idx_sw_diag)
        savaeruv  = (isw .eq. idx_uv_diag)
        savaernir = (isw .eq. idx_nir_diag)

        do k = 1, pver

           ! Determine the fraction of the particle that is core.
           coremmr(:ncol)       = 0._r8
           corerfr(:ncol)       = 0._r8
           corerfi(:ncol)       = 0._r8
           corehygro(:ncol)     = 0._r8
           corevol(:ncol)       = 0._r8

           coredustmmr(:ncol)       = 0._r8
           corebcmmr(:ncol)       = 0._r8

           shellmmr(:ncol)      = 0._r8
           shellrfr(:ncol)      = 0._r8
           shellrfi(:ncol)      = 0._r8
           shellhygro(:ncol)    = 0._r8
           shellvol(:ncol)      = 0._r8

           ! particlemmr(:ncol)   = 0._r8
           ! particlerfr(:ncol)   = 0._r8
           ! particlerfi(:ncol)   = 0._r8
           ! particlehygro(:ncol) = 0._r8
           ! particlevol(:ncol)   = 0._r8

           corefrac(:ncol)      = 0._r8
           dryvol(:ncol) = 0._r8
           crefin(:ncol) = 0._r8

           dustvol(:ncol) = 0._r8

           scatdust(:ncol)     = 0._r8
           absdust(:ncol)      = 0._r8
           hygrodust(:ncol)    = 0._r8
           scatso4(:ncol)      = 0._r8
           absso4(:ncol)       = 0._r8
           hygroso4(:ncol)     = 0._r8
           scatbc(:ncol)       = 0._r8
           absbc(:ncol)        = 0._r8
           hygrobc(:ncol)      = 0._r8
           scatpom(:ncol)      = 0._r8
           abspom(:ncol)       = 0._r8
           hygropom(:ncol)     = 0._r8
           scatsoa(:ncol)      = 0._r8
           abssoa(:ncol)       = 0._r8
           hygrosoa(:ncol)     = 0._r8
           scatseasalt(:ncol)  = 0._r8
           absseasalt(:ncol)   = 0._r8
           hygroseasalt(:ncol) = 0._r8

           ! aerosol species loop
           do l = 1, nspec
              call rad_cnst_get_bin_mmr_by_idx(list_idx, m, l, 'a', state, pbuf, specmmr)

              ! NOTE: In the future could get more properties like refractive
              ! index and kappa to allow calculation of core and shell properties
              ! from multiple components.
              call rad_cnst_get_bin_props_by_idx(list_idx, m, l, density_aer=specdens, &
                   refindex_aer_sw=specrefindex, spectype=spectype, &
                   specmorph=specmorph, hygro_aer=hygro_aer)

              vol(:ncol)      = specmmr(:ncol,k) / specdens
              dryvol(:ncol)   = dryvol(:ncol) + vol(:ncol)
              crefin(:ncol)   = crefin(:ncol) + vol(:ncol) * specrefindex(isw)

              if (savaervis) then
                 specrefr = real(specrefindex(isw))
                 specrefi = aimag(specrefindex(isw))
                 burden(:ncol) = burden(:ncol) + specmmr(:ncol,k)*mass(:ncol,k)
                 if (trim(spectype) == 'dust') then
                    burdendust(:ncol) = burdendust(:ncol) + specmmr(:ncol,k)*mass(:ncol,k)
                    dustvol(:ncol)    = vol(:ncol)
                    scatdust(:ncol)   = vol(:ncol)*specrefr
                    absdust(:ncol)    =-vol(:ncol)*specrefi
                    hygrodust(:ncol)  = vol(:ncol)*hygro_aer
                 end if
                 if (trim(spectype) == 'sulfate') then
                    burdenso4(:ncol) = burdenso4(:ncol) + specmmr(:ncol,k)*mass(:ncol,k)
                    scatso4(:ncol)   = vol(:ncol)*specrefr
                    absso4(:ncol)    =-vol(:ncol)*specrefi
                    hygroso4(:ncol)  = vol(:ncol)*hygro_aer
                 end if
                 if (trim(spectype) == 'black-c') then
                    burdenbc(:ncol) = burdenbc(:ncol) + specmmr(:ncol,k)*mass(:ncol,k)
                    scatbc(:ncol)   = vol(:ncol)*specrefr
                    absbc(:ncol)    =-vol(:ncol)*specrefi
                    hygrobc(:ncol)  = vol(:ncol)*hygro_aer
                 end if
                 if (trim(spectype) == 'p-organic') then
                    burdenpom(:ncol) = burdenpom(:ncol) + specmmr(:ncol,k)*mass(:ncol,k)
                    scatpom(:ncol)   = vol(:ncol)*specrefr
                    abspom(:ncol)    =-vol(:ncol)*specrefi
                    hygropom(:ncol)  = vol(:ncol)*hygro_aer
                 end if
                 if (trim(spectype) == 's-organic') then
                    burdensoa(:ncol) = burdensoa(:ncol) + specmmr(:ncol,k)*mass(:ncol,k)
                    scatsoa(:ncol)   = vol(:ncol)*specrefr
                    abssoa(:ncol)    =-vol(:ncol)*specrefi
                    hygrosoa(:ncol)  = vol(:ncol)*hygro_aer
                 end if
                 if (trim(spectype) == 'seasalt') then
                    burdenseasalt(:ncol) = burdenseasalt(:ncol) + specmmr(:ncol,k)*mass(:ncol,k)
                    scatseasalt(:ncol)   = vol(:ncol)*specrefr
                    absseasalt(:ncol)    =-vol(:ncol)*specrefi
                    hygroseasalt(:ncol)  = vol(:ncol)*hygro_aer
                 end if
              end if



              if (trim(specmorph) == 'core') then
                 if (trim(spectype) == 'dust') then
                    coredustmmr(:ncol)    = coredustmmr(:ncol) + specmmr(:ncol,k)
                 end if
                 if (trim(spectype) == 'black-c') then
                    corebcmmr(:ncol)    = corebcmmr(:ncol) + specmmr(:ncol,k)
                 end if
                 coremmr(:ncol)    = coremmr(:ncol) + specmmr(:ncol,k)
                 corevol(:ncol)    = corevol(:ncol) + vol(:ncol)
                 corerfr(:ncol)    = corerfr(:ncol) + vol(:ncol) * real(specrefindex(isw))
                 corerfi(:ncol)    = corerfi(:ncol) + vol(:ncol) * aimag(specrefindex(isw))
                 corehygro(:ncol)  = corehygro(:ncol) + vol(:ncol) * hygro_aer
              else if (trim(specmorph) == 'shell') then
                 shellmmr(:ncol)    = shellmmr(:ncol) + specmmr(:ncol,k)
                 shellvol(:ncol)    = shellvol(:ncol) + vol(:ncol)
                 shellrfr(:ncol)    = shellrfr(:ncol) + vol(:ncol) * real(specrefindex(isw))
                 shellrfi(:ncol)    = shellrfi(:ncol) + vol(:ncol) * aimag(specrefindex(isw))
                 shellhygro(:ncol)  = shellhygro(:ncol) + vol(:ncol) * hygro_aer
!               else if (trim(specmorph) == 'particle') then
!                  particlemmr(:ncol)    = specmmr(:ncol,k)
!                  particlevol(:ncol)    = vol(:ncol)
!                  particlerfr(:ncol)    = vol(:ncol) * real(specrefindex(isw))
!                  particlerfi(:ncol)    = vol(:ncol) * aimag(specrefindex(isw))
!                  particlehygro(:ncol)  = vol(:ncol) * hygro_aer
              else
                 call endrun("coreshell_aero_sw: Unknown specmorph "//trim(specmorph))
              end if

           end do ! species loop
           ! particlemmr(:ncol)    = shellmmr(:ncol) + coremmr(:ncol)
           ! particlevol(:ncol)    = shellvol(:ncol) + corevol(:ncol)
           ! particlerfr(:ncol)    = shellrfr(:ncol) + corerfr(:ncol)
           ! particlerfi(:ncol)    = shellrfi(:ncol) + corerfi(:ncol)
           ! particlehygro(:ncol)  = shellhygro(:ncol) + corehygro(:ncol)


           ! If one of the species is a particle, then it represents the mass of the
           ! whole particle, put includes mass that should be part of the shell. So
           ! subtract the mass of the explicit core and shell to determine the
           ! additional contribution.
           !
           ! NOTE: This is consistent with the way that CARMA defined the concentration
           ! element in the group.
           ! Not used currently
           ! do i = 1, ncol
           !   if (particlemmr(i) .gt. (shellmmr(i) + coremmr(i))) then
           !     shellmmr(i)    = shellmmr(i) + (particlemmr(i) - coremmr(i) - shellmmr(i))
           !     shellvol(i)    = shellvol(i) + (particlevol(i) - corevol(i) - shellvol(i))
           !     shellrfr(i)    = shellrfr(i) + (particlerfr(i) - corerfr(i) - shellrfr(i))
           !     shellrfi(i)    = shellrfi(i) + (particlerfi(i) - corerfi(i) - shellrfi(i))
           !     shellhygro(i)  = shellhygro(i) + (particlehygro(i) - corehygro(i) - shellhygro(i))
           !   end if
           ! end do

           ! If the particles were to swell, then water would need to be added to
           ! the shell and the volume and refractive index adjusted accordingly.

           ! Determine the core/shell ratio (by mass).
           totalmmr(:ncol) = (coremmr(:ncol) + shellmmr(:ncol))
           do i = 1, ncol
              if (totalmmr(i) .gt. 0._r8) then
                 corefrac(i) = coremmr(i) / (totalmmr(i))
              end if
              corefrac(i) = max(0._r8, min(1.0_r8, corefrac(i)))
           end do
           ! Determine the bc/dust ratio (by mass).
           bcdustmmr(:ncol) = (corebcmmr(:ncol) + coredustmmr(:ncol))
           do i = 1, ncol
              if (bcdustmmr(i) .gt. 0._r8) then
                 bcdust(i) = corebcmmr(i) / (bcdustmmr(i))
              end if
              bcdust(i) = max(0._r8, min(1.0_r8, bcdust(i)))
           end do

           if (associated(tbl_corefrac)) then
              where ( corefrac(:ncol) < minval(tbl_corefrac) )
                 corefrac(:ncol) = minval(tbl_corefrac)
              end where
              where ( corefrac(:ncol) > maxval(tbl_corefrac) )
                 corefrac(:ncol) = maxval(tbl_corefrac)
              end where
           endif

           if (associated(tbl_bcdust)) then
              where ( bcdust(:ncol) < minval(tbl_bcdust) )
                 bcdust(:ncol) = minval(tbl_bcdust)
              end where
              where ( bcdust(:ncol) > maxval(tbl_bcdust) )
                 bcdust(:ncol) = maxval(tbl_bcdust)
              end where
           endif

           ! Now interplate in the core/shell dimension.
           do i = 1, ncol

              select case (trim(opticstype))
              case('sectional')
                 pext(i) = table_interp( tbl_corefrac, extpsw(:,isw), corefrac(i) )
                 pabs(i) = table_interp( tbl_corefrac, abspsw(:,isw), corefrac(i) )
                 pasm(i) = table_interp( tbl_corefrac, asmpsw(:,isw), corefrac(i) )
              case('hygroscopic_coreshell')
                 pext(i) = table_interp( tbl_relh, tbl_corefrac, tbl_bcdust, tbl_kap, &
                                         h_ext_coreshell(:,isw,:,:,:), &
                                         relh(i,k), corefrac(i), bcdust(i), crkappa(i,k) )
                 pabs(i) = (1._r8 - table_interp( tbl_relh, tbl_corefrac, tbl_bcdust, tbl_kap, &
                                         h_ssa_coreshell(:,isw,:,:,:), &
                                         relh(i,k), corefrac(i), bcdust(i), crkappa(i,k) ) ) * pext(i)
                 pasm(i) = table_interp( tbl_relh, tbl_corefrac, tbl_bcdust, tbl_kap, &
                                         h_asm_coreshell(:,isw,:,:,:), &
                                         relh(i,k), corefrac(i), bcdust(i), crkappa(i,k) )
              case('hygroscopic_wtp')
                 pext(i) = table_interp( tbl_wgtpct, h_ext_wtp(:,isw), wgtpct(i,k) )
                 pabs(i) = (1._r8 - table_interp( tbl_wgtpct, h_ssa_wtp(:,isw), wgtpct(i,k) ) ) * pext(i)
                 pasm(i) = table_interp( tbl_wgtpct, h_asm_wtp(:,isw), wgtpct(i,k) )
              case default
                 call endrun('coreshell_aero_sw: unsupported opticstype: '//trim(opticstype))
              end select

              specpext(i) = pext(i)
              pext(i) = pext(i)*totalmmr(i)
              pabs(i) = pabs(i)*totalmmr(i)
              pabs(i) = max(0._r8,pabs(i))
              pabs(i) = min(pext(i),pabs(i))

              palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)

              dopaer(i) = pext(i)*mass(i,k)

              aodtot(i) = aodtot(i) + dopaer(i)
           end do

           if (savaeruv) then
              do i = 1, ncol
                 extinctuv(i,k) = extinctuv(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                 aoduv(i) = aoduv(i) + dopaer(i)
                 if (k.le.troplev(i)) then
                    aoduvst(i) = aoduvst(i) + dopaer(i)
                 end if
              end do
           end if

           if (savaernir) then
              do i = 1, ncol
                 extinctnir(i,k) = extinctnir(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                 aodnir(i) = aodnir(i) + dopaer(i)
                 if (k.le.troplev(i)) then
                    aodnirst(i) = aodnirst(i) + dopaer(i)
                 end if
              end do
           endif

           ! Save aerosol optical depth at longest visible wavelength
           ! sum over layers
           if (savaervis) then

              dryvol_m3(:ncol) = four_thirds_pi * (dryr(:ncol,k)**3)
              wetvol_m3(:ncol) = four_thirds_pi * (wetr(:ncol,k)**3)
              watervol(:ncol) = wetvol(:ncol)-dryvol(:ncol)
              where (dryvol_m3(:ncol)>0._r8)
                 wet_dry_frac(:ncol) = wetvol_m3(:ncol)/dryvol_m3(:ncol)
              elsewhere
                 wet_dry_frac(:ncol) = 0.0_r8
              end where
              wetvol(:ncol) = dryvol(:ncol)*wet_dry_frac(:ncol)

              ! aerosol extinction (/m)
              do i = 1, ncol

                 watervol(i) = max(0._r8,watervol(i))

                 extinct(i,k) = extinct(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                 absorb(i,k)  = absorb(i,k) + pabs(i)*air_density(i,k)
                 aodvis(i)    = aodvis(i) + dopaer(i)
                 aodabs(i)    = aodabs(i) + pabs(i)*mass(i,k)
                 aodbin(i)    = aodbin(i) + dopaer(i)
                 ssavis(i)    = ssavis(i) + dopaer(i)*palb(i)
                 asymvis(i)   = asymvis(i) + dopaer(i)*pasm(i)
                 asymext(i,k) = asymext(i,k) + dopaer(i)*pasm(i)*air_density(i,k)/mass(i,k)
                 if (k.le.troplev(i)) then
                    aodvisst(i) = aodvisst(i) + dopaer(i)
                 end if

                 if (wetvol(i) > 1.e-40_r8 .and. vol(i)>0._r8) then

                    dustaodbin(i) = dustaodbin(i) + dopaer(i)*dustvol(i)/wetvol(i)

                    ! partition optical depth into contributions from each constituent
                    ! assume contribution is proportional to refractive index X volume

                    scath2o        = watervol(i)*real(crefwsw(isw))
                    absh2o         = -watervol(i)*aimag(crefwsw(isw))
                    sumscat        = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                                     scatdust(i) + scatseasalt(i) + scath2o
                    sumabs         = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                                     absdust(i) + absseasalt(i) + absh2o
                    sumhygro       = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                                     hygrodust(i) + hygroseasalt(i)

                    scatdust(i)    = (scatdust(i) + scath2o*hygrodust(i)/sumhygro)/sumscat
                    absdust(i)     = (absdust(i) + absh2o*hygrodust(i)/sumhygro)/sumabs

                    scatso4(i)     = (scatso4(i) + scath2o*hygroso4(i)/sumhygro)/sumscat
                    absso4(i)      = (absso4(i) + absh2o*hygroso4(i)/sumhygro)/sumabs

                    scatpom(i)     = (scatpom(i) + scath2o*hygropom(i)/sumhygro)/sumscat
                    abspom(i)      = (abspom(i) + absh2o*hygropom(i)/sumhygro)/sumabs

                    scatsoa(i)     = (scatsoa(i) + scath2o*hygrosoa(i)/sumhygro)/sumscat
                    abssoa(i)      = (abssoa(i) + absh2o*hygrosoa(i)/sumhygro)/sumabs

                    scatbc(i)      = (scatbc(i) + scath2o*hygrobc(i)/sumhygro)/sumscat
                    absbc(i)       = (absbc(i) + absh2o*hygrobc(i)/sumhygro)/sumabs

                    scatseasalt(i) = (scatseasalt(i) + scath2o*hygroseasalt(i)/sumhygro)/sumscat
                    absseasalt(i)  = (absseasalt(i) + absh2o*hygroseasalt(i)/sumhygro)/sumabs

                    aodabsbc(i)    = aodabsbc(i) + absbc(i)*dopaer(i)*(1.0_r8-palb(i))

                    aodc           = (absdust(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust(i))*dopaer(i)
                    dustaod(i)     = dustaod(i) + aodc

                    aodc           = (absso4(i)*(1.0_r8 - palb(i)) + palb(i)*scatso4(i))*dopaer(i)
                    so4aod(i)      = so4aod(i) + aodc

                    aodc           = (abspom(i)*(1.0_r8 - palb(i)) + palb(i)*scatpom(i))*dopaer(i)
                    pomaod(i)      = pomaod(i) + aodc

                    aodc           = (abssoa(i)*(1.0_r8 - palb(i)) + palb(i)*scatsoa(i))*dopaer(i)
                    soaaod(i)      = soaaod(i) + aodc

                    aodc           = (absbc(i)*(1.0_r8 - palb(i)) + palb(i)*scatbc(i))*dopaer(i)
                    bcaod(i)       = bcaod(i) + aodc

                    aodc           = (absseasalt(i)*(1.0_r8 - palb(i)) + palb(i)*scatseasalt(i))*dopaer(i)
                    seasaltaod(i)  = seasaltaod(i) + aodc

                 endif


              end do
           endif

           do i = 1, ncol

              if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 300._r8)) then

                 if (dopaer(i) <= -1.e-10_r8) then
                    write(iulog,*) "ERROR: Negative aerosol optical depth &
                         &in this layer."
                 else
                    write(iulog,*) "WARNING: Aerosol optical depth is &
                         &unreasonably high in this layer."
                 end if

                 write(iulog,*) 'dopaer(', i, ',', k, ',', m, ',', lchnk, ')=', dopaer(i)
                 ! write(iulog,*) 'itab,jtab,ttab,utab=',itab(i),jtab(i),ttab(i),utab(i)
                 write(iulog,*) 'k=', k, ' pext=', pext(i), ' specext=', specpext(i)
                 !write(iulog,*) 'wetvol=', wetvol(i), ' dryvol=', dryvol(i), ' watervol=', watervol(i)
                 ! write(iulog,*) 'cext=',(cext(i,l),l=1,ncoef)
                 ! write(iulog,*) 'crefin=',crefin(i)
                 write(iulog,*) 'nspec=', nspec
                 ! write(iulog,*) 'cheb=', (cheb(nc,m,i,k),nc=2,ncoef)
                 do l = 1, nspec
                    call rad_cnst_get_bin_mmr_by_idx(list_idx, m, l, 'a', state, pbuf, specmmr)
                    call rad_cnst_get_bin_props_by_idx(list_idx, m, l, density_aer=specdens, &
                         refindex_aer_sw=specrefindex)
                    volf = specmmr(i,k)/specdens
                    write(iulog,*) 'l=', l, 'vol(l)=', volf
                    write(iulog,*) 'isw=', isw, 'specrefindex(isw)=', specrefindex(isw)
                    write(iulog,*) 'specdens=', specdens
                 end do

                 nerr_dopaer = nerr_dopaer + 1
                 !if (nerr_dopaer >= nerrmax_dopaer) then
                 if (dopaer(i) < -1.e-10_r8) then
                    write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                    call endrun('exit from '//subname)
                 end if

              end if
           end do

           do i=1,ncol
              tauxar(i,k,isw) = tauxar(i,k,isw) + dopaer(i)
              wa(i,k,isw)     = wa(i,k,isw)     + dopaer(i)*palb(i)
              ga(i,k,isw)     = ga(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)
              fa(i,k,isw)     = fa(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
           end do

        end do ! pver
     end do ! sw bands

     ! bin diagnostics
     ! The diagnostics are currently only output for the climate list.  Code mods will
     ! be necessary to provide output for the rad_diag lists.
     if (list_idx == 0) then

        write(outname,'(a,i2.2)') 'BURDENdn', m
        call outfld(trim(outname), burden, pcols, lchnk)

        write(outname,'(a,i2.2)') 'AODdnBIN', m
        call outfld(trim(outname), aodbin, pcols, lchnk)

        write(outname,'(a,i2.2)') 'AODdnDUST', m
        call outfld(trim(outname), dustaodbin, pcols, lchnk)

        do i = 1, nnite
           burden(idxnite(i))  = fillvalue
           aodbin(idxnite(i)) = fillvalue
           dustaodbin(idxnite(i)) = fillvalue
        end do


        write(outname,'(a,i2.2)') 'BURDEN', m
        call outfld(trim(outname), burden, pcols, lchnk)

        write(outname,'(a,i2.2)') 'AODDUST', m
        call outfld(trim(outname), dustaodbin, pcols, lchnk)

        write(outname,'(a,i2.2)') 'AODBIN', m
        call outfld(trim(outname), aodbin, pcols, lchnk)
     end if

   end do ! nbins


   call outfld('EXTINCTdn'//diag(list_idx),  extinct, pcols, lchnk)
   call outfld('ABSORBdn'//diag(list_idx),   absorb,  pcols, lchnk)
   call outfld('AODVISdn'//diag(list_idx),   aodvis,  pcols, lchnk)
   call outfld('AODABSdn'//diag(list_idx),   aodabs,  pcols, lchnk)
   call outfld('AODTOTdn'//diag(list_idx), aodtot,  pcols, lchnk)
   call outfld('AODVISstdn'//diag(list_idx), aodvisst,pcols, lchnk)
   call outfld('EXTxASYMdn'//diag(list_idx), asymext, pcols, lchnk)


   !Output visible band diagnostics for quantities summed over the bins
   !These fields are put out for diagnostic lists as well as the climate list.
   do i = 1, nnite
      extinct(idxnite(i),:) = fillvalue
      absorb(idxnite(i),:)  = fillvalue
      aodvis(idxnite(i))    = fillvalue
      aodtot(idxnite(i))    = fillvalue
      aodabs(idxnite(i))    = fillvalue
      aodvisst(idxnite(i))  = fillvalue
      asymext(idxnite(i),:) = fillvalue
   end do

   call outfld('EXTINCT'//diag(list_idx),  extinct, pcols, lchnk)
   call outfld('ABSORB'//diag(list_idx),   absorb,  pcols, lchnk)
   call outfld('AODVIS'//diag(list_idx),   aodvis,  pcols, lchnk)
   call outfld('AODTOT'//diag(list_idx),   aodtot,  pcols, lchnk)
   call outfld('AODABS'//diag(list_idx),   aodabs,  pcols, lchnk)
   call outfld('AODVISst'//diag(list_idx), aodvisst,pcols, lchnk)
   call outfld('EXTxASYM'//diag(list_idx), asymext, pcols, lchnk)

   ! These diagnostics are output only for climate list
   if (list_idx == 0) then
      do i = 1, ncol
         if (aodvis(i) > 1.e-10_r8) then
            ssavis(i) = ssavis(i)/aodvis(i)
         else
            ssavis(i) = 0.925_r8
         endif
      end do

      call outfld('SSAVISdn',        ssavis,        pcols, lchnk)
      call outfld('AODxASYMdn',      asymvis,       pcols, lchnk)

      call outfld('EXTINCTUVdn',     extinctuv,     pcols, lchnk)
      call outfld('EXTINCTNIRdn',    extinctnir,    pcols, lchnk)
      call outfld('AODUVdn',         aoduv,         pcols, lchnk)
      call outfld('AODNIRdn',        aodnir,        pcols, lchnk)
      call outfld('AODUVstdn',       aoduvst,       pcols, lchnk)
      call outfld('AODNIRstdn',      aodnirst,      pcols, lchnk)

      call outfld('BURDENDUSTdn',    burdendust,    pcols, lchnk)
      call outfld('BURDENSO4dn' ,    burdenso4,     pcols, lchnk)
      call outfld('BURDENPOMdn' ,    burdenpom,     pcols, lchnk)
      call outfld('BURDENSOAdn' ,    burdensoa,     pcols, lchnk)
      call outfld('BURDENBCdn'  ,    burdenbc,      pcols, lchnk)
      call outfld('BURDENSEASALTdn', burdenseasalt, pcols, lchnk)

      call outfld('AODABSBCdn',      aodabsbc,      pcols, lchnk)

      call outfld('AODDUSTdn',       dustaod,       pcols, lchnk)
      call outfld('AODSO4dn',        so4aod,        pcols, lchnk)
      call outfld('AODPOMdn',        pomaod,        pcols, lchnk)
      call outfld('AODSOAdn',        soaaod,        pcols, lchnk)
      call outfld('AODBCdn',         bcaod,         pcols, lchnk)
      call outfld('AODSSdn',         seasaltaod,    pcols, lchnk)


      do i = 1, nnite
         ssavis(idxnite(i))     = fillvalue
         asymvis(idxnite(i))    = fillvalue

         aoduv(idxnite(i))      = fillvalue
         aodnir(idxnite(i))     = fillvalue
         aoduvst(idxnite(i))    = fillvalue
         aodnirst(idxnite(i))   = fillvalue
         extinctuv(idxnite(i),:)  = fillvalue
         extinctnir(idxnite(i),:) = fillvalue

         burdendust(idxnite(i)) = fillvalue
         burdenso4(idxnite(i))  = fillvalue
         burdenpom(idxnite(i))  = fillvalue
         burdensoa(idxnite(i))  = fillvalue
         burdenbc(idxnite(i))   = fillvalue
         burdenseasalt(idxnite(i)) = fillvalue

         aodabsbc(idxnite(i))   = fillvalue

         dustaod(idxnite(i))    = fillvalue
         so4aod(idxnite(i))     = fillvalue
         pomaod(idxnite(i))     = fillvalue
         soaaod(idxnite(i))     = fillvalue
         bcaod(idxnite(i))      = fillvalue
         seasaltaod(idxnite(i)) = fillvalue
       end do

      call outfld('SSAVIS',        ssavis,        pcols, lchnk)
      call outfld('AODxASYM',      asymvis,       pcols, lchnk)

      call outfld('EXTINCTUV',     extinctuv,     pcols, lchnk)
      call outfld('EXTINCTNIR',    extinctnir,    pcols, lchnk)
      call outfld('AODUV',         aoduv,         pcols, lchnk)
      call outfld('AODNIR',        aodnir,        pcols, lchnk)
      call outfld('AODUVst',       aoduvst,       pcols, lchnk)
      call outfld('AODNIRst',      aodnirst,      pcols, lchnk)

      call outfld('BURDENDUST',    burdendust,    pcols, lchnk)
      call outfld('BURDENSO4' ,    burdenso4,     pcols, lchnk)
      call outfld('BURDENPOM' ,    burdenpom,     pcols, lchnk)
      call outfld('BURDENSOA' ,    burdensoa,     pcols, lchnk)
      call outfld('BURDENBC'  ,    burdenbc,      pcols, lchnk)
      call outfld('BURDENSEASALT', burdenseasalt, pcols, lchnk)

      call outfld('AODABSBC',      aodabsbc,      pcols, lchnk)

      call outfld('AODDUST',       dustaod,       pcols, lchnk)
      call outfld('AODSO4',        so4aod,        pcols, lchnk)
      call outfld('AODPOM',        pomaod,        pcols, lchnk)
      call outfld('AODSOA',        soaaod,        pcols, lchnk)
      call outfld('AODBC',         bcaod,         pcols, lchnk)
      call outfld('AODSS',         seasaltaod,    pcols, lchnk)
   end if

end subroutine coreshell_aero_sw

!===============================================================================

subroutine coreshell_aero_lw(list_idx, state, pbuf, tauxar)


   ! calculates aerosol lw radiative properties

   integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state    ! state variables

   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(inout) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

   ! Local variables
   integer :: i, ifld, ilw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nbins
   integer :: nspec
   character(len=32)    :: spectype            ! species type
   real(r8)             :: hygro_aer           !

   real(r8) :: mass(pcols,pver) ! layer mass

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index

   real(r8) :: particlemmr(pcols)   ! mmr of aerosol particle
   real(r8) :: particlerfr(pcols)   ! real refractive index of aerosol particle
   real(r8) :: particlerfi(pcols)   ! imaginary refractive index of aerosol particle
   real(r8) :: particlehygro(pcols) ! hygroscopicity of aerosol particle
   real(r8) :: particlevol(pcols)   ! volume of aerosol particle

   real(r8) :: coredustmmr(pcols)   !
   real(r8) :: corebcmmr(pcols)   !

   real(r8) :: coremmr(pcols)   ! mmr of aerosol core
   real(r8) :: corerfr(pcols)   ! real refractive index of aerosol core
   real(r8) :: corerfi(pcols)   ! imaginary refractive index of aerosol core
   real(r8) :: corehygro(pcols) ! hygroscopicity of aerosol core
   real(r8) :: corevol(pcols)   ! volume of aerosol core

   real(r8) :: shellmmr(pcols)   ! mmr of aerosol shell
   real(r8) :: shellrfr(pcols)   ! real refractive index of aerosol shell
   real(r8) :: shellrfi(pcols)   ! imaginary refractive index of aerosol shell
   real(r8) :: shellhygro(pcols) ! hygroscopicity of aerosol shell
   real(r8) :: shellvol(pcols)   ! volume of aerosol shell

   real(r8) :: totalmmr(pcols)   ! total mmr of the aerosol
   real(r8) :: corefrac(pcols)   ! mass fraction that is core
   real(r8) :: bcdustmmr(pcols)  ! total mmr of dust + bc aerosol
   real(r8) :: bcdust(pcols)     ! mass fraction of bc vs (bc + dust)

   real(r8) :: vol(pcols)       ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)    ! volume concentration of aerosol bin (m3/kg)
   real(r8) :: refr(pcols)      ! real part of refractive index
   real(r8) :: refi(pcols)      ! imaginary part of refractive index
   complex(r8) :: crefin(pcols) ! complex refractive index
   real(r8), pointer :: tbl_corefrac(:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: tbl_kap(:) ! table of kappa of bins
   real(r8), pointer :: tbl_relh(:) ! table of relh of bins
   real(r8), pointer :: tbl_bcdust(:) ! table of bc-dust ration of bins
   real(r8), pointer :: tbl_wgtpct(:)

   real(r8), pointer :: absplw(:,:) ! specific absorption
   real(r8), pointer :: lw_hygro_coreshell_abs(:,:,:,:,:) ! specific absorption
   real(r8), pointer :: lw_hygro_abs_wtp(:,:) ! specific absorption

   real(r8) :: pabs(pcols)      ! parameterized specific absorption (m2/kg)
   real(r8) :: dopaer(pcols)    ! aerosol optical depth in layer

! for table lookup into relh grid
   real(r8) :: sate(pcols,pver)     ! saturation vapor pressure
   real(r8) :: satq(pcols,pver)     ! saturation specific humidity
   real(r8) :: relh(pcols,pver)
   integer :: nrelh, nbcdust, nkap
   integer :: nwtp

   real(r8) :: lwabs(pcols,pver)

   real(r8), pointer, dimension(:,:) :: crkappa   ! kappa
   real(r8), pointer, dimension(:,:) :: wgtpct  ! weight percent
   real(r8), pointer, dimension(:,:) :: cldn

   character(len=ot_length) :: opticstype
   character(len=32) :: bin_name
   character(len=32) :: specmorph

   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf             ! volume fraction of insoluble aerosol

   character(len=*), parameter :: subname = 'modal_aero_lw'
   !----------------------------------------------------------------------------

   nullify(tbl_corefrac)
   nullify(tbl_kap)
   nullify(tbl_relh)
   nullify(tbl_bcdust)
   nullify(tbl_wgtpct)
   nullify(absplw)
   nullify(lw_hygro_coreshell_abs)
   nullify(lw_hygro_abs_wtp)

   nullify(crkappa)
   nullify(wgtpct)
   nullify(cldn)

   lchnk = state%lchnk
   ncol  = state%ncol

   lwabs = 0._r8

   ! initialize output variables

   ! dry mass in each cell
   mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

   ! calculate relative humidity for table lookup into rh grid
   call qsat(state%t(1:ncol,1:pver), state%pmid(1:ncol,1:pver), sate(1:ncol,1:pver), satq(1:ncol,1:pver), ncol,pver)

   relh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / satq(1:ncol,1:pver)
   !There is a  different way to diagnose sub-grid rh for clear-sky only that makes physically
   ! more sense but is not used in MAM. This needs to be tested if it makes any difference
   ! Associate pointers with physics buffer fields
   ! itim = pbuf_old_tim_idx()
   ! call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cldn, (/1,1,itim/),(/pcols,pver,1/))
   ! rh(1:ncol,1:pver) = (state%q(1:ncol,1:pver,1)-cldn(1:ncol,1:pver)*qsat(1:ncol,1:pver)) / qsat(1:ncol,1:pver)
   relh(1:ncol,1:pver) = max(1.e-20_r8,relh(1:ncol,1:pver))

   ! loop over all aerosol bins
   call rad_cnst_get_info(list_idx, nbins=nbins)

   do m = 1, nbins

      ! get bin info
      call rad_cnst_get_info_by_bin(list_idx, m, nspec=nspec, bin_name=bin_name)
      ! get optics type
      call rad_cnst_get_bin_props(list_idx, m, opticstype=opticstype)

      select case (trim(opticstype))
      case('sectional')
         ! get bin properties
         call rad_cnst_get_bin_props(list_idx, m, &
              corefrac=tbl_corefrac, &
              absplw=absplw)
      case('hygroscopic_coreshell')
         ! get optical properties for hygroscopic aerosols
         call rad_cnst_get_bin_props(list_idx, m, &
              corefrac=tbl_corefrac, &
              kap=tbl_kap, nkap=nkap, bcdust=tbl_bcdust, nbcdust=nbcdust, &
              relh=tbl_relh, nrelh=nrelh, &
              lw_hygro_coreshell_ext=lw_hygro_coreshell_abs)
         ! get kappa
         ! need to read in bin_name
         call pbuf_get_field(pbuf, pbuf_get_index(trim(bin_name)//"_kappa"),crkappa)
         where ( crkappa(:ncol,:) < minval(tbl_kap) )
            crkappa(:ncol,:) = minval(tbl_kap)
         end where
         where ( crkappa(:ncol,:) > maxval(tbl_kap) )
            crkappa(:ncol,:) = maxval(tbl_kap)
         end where
         where ( relh(:ncol,:) < minval(tbl_relh) )
            relh(:ncol,:) = minval(tbl_relh)
         end where
         where ( relh(:ncol,:) > maxval(tbl_relh) )
            relh(:ncol,:) = maxval(tbl_relh)
         end where
      case('hygroscopic_wtp')
         ! get optical properties for hygroscopic aerosols
         call rad_cnst_get_bin_props(list_idx, m, &
              wgtpct=tbl_wgtpct,nwtp=nwtp, &
              lw_hygro_ext_wtp=lw_hygro_abs_wtp)
         ! determine weight precent of H2SO4/H2O solution
         call pbuf_get_field(pbuf, pbuf_get_index('WTP'),wgtpct)
         where ( wgtpct(:ncol,:) < minval(tbl_wgtpct) )
            wgtpct(:ncol,:) = minval(tbl_wgtpct)
         end where
         where ( wgtpct(:ncol,:) > maxval(tbl_wgtpct) )
            wgtpct(:ncol,:) = maxval(tbl_wgtpct)
         end where
      case default
         call endrun('coreshell_aero_lw: unsupported opticstype: '//trim(opticstype))
      end select

      do ilw = 1, nlwbands

         do k = 1, pver

            ! form bulk refractive index. Use volume mixing for infrared
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8

            ! Determine the fraction of the particle that is core.
            coremmr(:ncol)       = 0._r8
            corerfr(:ncol)       = 0._r8
            corerfi(:ncol)       = 0._r8
            corehygro(:ncol)     = 0._r8
            corevol(:ncol)       = 0._r8

            coredustmmr(:ncol)   = 0._r8
            corebcmmr(:ncol)     = 0._r8

            shellmmr(:ncol)      = 0._r8
            shellrfr(:ncol)      = 0._r8
            shellrfi(:ncol)      = 0._r8
            shellhygro(:ncol)    = 0._r8
            shellvol(:ncol)      = 0._r8

            ! not used currently
            ! particlemmr(:ncol)   = 0._r8
            ! particlerfr(:ncol)   = 0._r8
            ! particlerfi(:ncol)   = 0._r8
            ! particlehygro(:ncol) = 0._r8
            ! particlevol(:ncol)   = 0._r8

            corefrac(:ncol) = 0._r8
            dryvol(:ncol) = 0._r8
            crefin(:ncol) = 0._r8
            bcdust(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec
               call rad_cnst_get_bin_mmr_by_idx(list_idx, m, l, 'a', state, pbuf, specmmr)

               ! NOTE: In the future could get more properties like refractive
               ! index and kappa to allow calculation of core and shell properties
               ! from multiple components.
               call rad_cnst_get_bin_props_by_idx(list_idx, m, l, density_aer=specdens, &
                                           refindex_aer_lw=specrefindex, spectype=spectype, &
                                           specmorph=specmorph, hygro_aer=hygro_aer)

               vol(:ncol)      = specmmr(:ncol,k) / specdens
               dryvol(:ncol)   = dryvol(:ncol) + vol(:ncol)
               crefin(:ncol)   = crefin(:ncol) + vol(:ncol) * specrefindex(ilw)

               if (trim(specmorph) == 'core') then
                  if (trim(spectype) == 'dust') then
                     coredustmmr(:ncol)    = coredustmmr(:ncol) + specmmr(:ncol,k)
                  end if
                  if (trim(spectype) == 'black-c') then
                     corebcmmr(:ncol)    = corebcmmr(:ncol) + specmmr(:ncol,k)
                  end if
                  coremmr(:ncol)    = coremmr(:ncol) + specmmr(:ncol,k)
                  corevol(:ncol)    = corevol(:ncol) + vol(:ncol)
                  corerfr(:ncol)    = corerfr(:ncol) + vol(:ncol) * real(specrefindex(ilw))
                  corerfi(:ncol)    = corerfi(:ncol) + vol(:ncol) * aimag(specrefindex(ilw))
                  corehygro(:ncol)  = corehygro(:ncol) + vol(:ncol) * hygro_aer
               else if (trim(specmorph) == 'shell') then
                  shellmmr(:ncol)    = shellmmr(:ncol) + specmmr(:ncol,k)
                  shellvol(:ncol)    = shellvol(:ncol) + vol(:ncol)
                  shellrfr(:ncol)    = shellrfr(:ncol) + vol(:ncol) * real(specrefindex(ilw))
                  shellrfi(:ncol)    = shellrfi(:ncol) + vol(:ncol) * aimag(specrefindex(ilw))
                  shellhygro(:ncol)  = shellhygro(:ncol) + vol(:ncol) * hygro_aer
               else
                 call endrun("coreshell_aero_lw: Unknown specmorph "//trim(specmorph))
               end if

               !if (trim(specmorph) == 'particle') then
               !   particlemmr(:ncol)    = specmmr(:ncol,k)
               !   particlevol(:ncol)    = vol(:ncol)
               !   particlerfr(:ncol)    = vol(:ncol) * real(specrefindex(ilw))
               !   particlerfi(:ncol)    = vol(:ncol) * aimag(specrefindex(ilw))
               !   particlehygro(:ncol)  = vol(:ncol) * hygro_aer
               !end if
            end do
            ! particlemmr(:ncol)    = shellmmr(:ncol) + coremmr(:ncol)
            ! particlevol(:ncol)    = shellvol(:ncol) + corevol(:ncol)
            ! particlerfr(:ncol)    = shellrfr(:ncol) + corerfr(:ncol)
            ! particlerfi(:ncol)    = shellrfi(:ncol) + corerfi(:ncol)
            ! particlehygro(:ncol)  = shellhygro(:ncol) + corehygro(:ncol)


            ! If one of the species is a particle, then it represents the mass of the
            ! whole particle, put includes mass that should be part of the shell. So
            ! subtract the mass of the explicit core and shell to determine the
            ! additional contribution.
            !
            ! NOTE: This is consistent with the way that CARMA defined the concentration
            ! element in the group.
            !do i = 1, ncol
            !  if (particlemmr(i) .gt. (shellmmr(i) + coremmr(i))) then
            !    shellmmr(i)    = shellmmr(i) + (particlemmr(i) - coremmr(i) - shellmmr(i))
            !    shellvol(i)    = shellvol(i) + (particlevol(i) - corevol(i) - shellvol(i))
            !    shellrfr(i)    = shellrfr(i) + (particlerfr(i) - corerfr(i) - shellrfr(i))
            !    shellrfi(i)    = shellrfi(i) + (particlerfi(i) - corerfi(i) - shellrfi(i))
            !    shellhygro(i)  = shellhygro(i) + (particlehygro(i) - corehygro(i) - shellhygro(i))
            !  end if
            !end do

            ! If the particles were to swell, then water would need to be added to
            ! the shell and the volume and refractive index adjusted accordingly.

            ! Determine the core/shell ratio (by mass).
            totalmmr(:ncol) = (coremmr(:ncol) + shellmmr(:ncol))
            do i = 1, ncol
               if (totalmmr(i) .gt. 0._r8) then
                  corefrac(i) = coremmr(i) / (totalmmr(i))
               end if
               corefrac(i) = max(0._r8, min(1.0_r8, corefrac(i)))
            end do
            ! Determine the bc/dust ratio (by mass).
            bcdustmmr(:ncol) = (corebcmmr(:ncol) + coredustmmr(:ncol))
            do i = 1, ncol
               if (bcdustmmr(i) .gt. 0._r8) then
                  bcdust(i) = corebcmmr(i) / (bcdustmmr(i))
               end if
               bcdust(i) = max(0._r8, min(1.0_r8, bcdust(i)))
            end do

            if (associated(tbl_corefrac)) then
               where ( corefrac(:ncol) < minval(tbl_corefrac) )
                  corefrac(:ncol) = minval(tbl_corefrac)
               end where
               where ( corefrac(:ncol) > maxval(tbl_corefrac) )
                  corefrac(:ncol) = maxval(tbl_corefrac)
               end where
            endif

            if (associated(tbl_bcdust)) then
               where ( bcdust(:ncol) < minval(tbl_bcdust) )
                  bcdust(:ncol) = minval(tbl_bcdust)
               end where
               where ( bcdust(:ncol) > maxval(tbl_bcdust) )
                  bcdust(:ncol) = maxval(tbl_bcdust)
               end where
            endif

            ! interpolate coefficients linear in refractive index

            do i = 1, ncol
               select case (trim(opticstype))
               case('sectional')
                  pabs(i) = table_interp( tbl_corefrac, absplw(:,ilw), corefrac(i) )
               case('hygroscopic_coreshell')
                  pabs(i) = table_interp( tbl_relh, tbl_corefrac,  tbl_bcdust, tbl_kap, &
                                          lw_hygro_coreshell_abs(:,ilw,:,:,:), &
                                          relh(i,k), corefrac(i), bcdust(i), crkappa(i,k) )
               case('hygroscopic_wtp')
                  pabs(i) =  table_interp( tbl_wgtpct, lw_hygro_abs_wtp(:,ilw), wgtpct(i,k) )
               case default
                  call endrun('coreshell_aero_lw: unsupported opticstype: '//trim(opticstype))
               end select
            end do

            ! parameterized optical properties
            do i = 1, ncol
               pabs(i)   = pabs(i)*totalmmr(i)
               pabs(i)   = max(0._r8,pabs(i))
               dopaer(i) = pabs(i)*mass(i,k)
               lwabs(i,k) = lwabs(i,k) + pabs(i)
            end do

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 300._r8)) then

                  if (dopaer(i) <= -1.e-10_r8) then
                     write(iulog,*) "ERROR: Negative aerosol optical depth &
                          &in this layer."
                  else
                     write(iulog,*) "WARNING: Aerosol optical depth is &
                          &unreasonably high in this layer."
                  end if

                  write(iulog,*) 'dopaer(',i,',',k,',',m,',',lchnk,')=', dopaer(i)
                  write(iulog,*) 'k=',k,' pabs=', pabs(i)
!                  write(iulog,*) 'cabs=', (cabs(i,l),l=1,ncoef)
                  write(iulog,*) 'crefin=', crefin(i)
                  write(iulog,*) 'nspec=', nspec
                  do l = 1,nspec
                      call rad_cnst_get_bin_mmr_by_idx(list_idx, m, l, 'a', state, pbuf, specmmr)
                      call rad_cnst_get_bin_props_by_idx(list_idx, m, l, density_aer=specdens, &
                                                 refindex_aer_lw=specrefindex)
                     volf = specmmr(i,k)/specdens
                     write(iulog,*) 'l=',l,'vol(l)=',volf
                     write(iulog,*) 'ilw=',ilw,' specrefindex(ilw)=',specrefindex(ilw)
                     write(iulog,*) 'specdens=',specdens
                  end do

                  nerr_dopaer = nerr_dopaer + 1
                  if (nerr_dopaer >= nerrmax_dopaer .or. dopaer(i) < -1.e-10_r8) then
                     write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     call endrun()
                  end if

               end if
            end do
            do i = 1, ncol
               tauxar(i,k,ilw) = tauxar(i,k,ilw) + dopaer(i)
            end do

         end do ! k = 1, pver

      end do  ! nlwbands

   end do ! m = 1, nbins

   call outfld('TOTABSLW', lwabs(:,:), pcols, lchnk)

   if (lw10um_indx>0) then
      call outfld('AODABSLW', tauxar(:,:,lw10um_indx), pcols, lchnk)
   end if

end subroutine coreshell_aero_lw

!===============================================================================
! Private routines
!===============================================================================

subroutine read_water_refindex(infilename)

   ! read water refractive index file and set module data

   character*(*), intent(in) :: infilename   ! modal optics filename

   ! Local variables

   integer            :: i, ierr
   type(file_desc_t)  :: ncid              ! pio file handle
   integer            :: did               ! dimension ids
   integer            :: dimlen            ! dimension lengths
   type(var_desc_t)   :: vid               ! variable ids
   real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
   real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared
   !----------------------------------------------------------------------------

   ! open file
   call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

   ! inquire dimensions.  Check that file values match parameter values.

   ierr = pio_inq_dimid(ncid, 'lw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nlwbands) then
      write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
      call endrun('read_modal_optics: bad lw_band value')
   endif

   ierr = pio_inq_dimid(ncid, 'sw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nswbands) then
      write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
      call endrun('read_modal_optics: bad sw_band value')
   endif

   ! read variables
   ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refrwsw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refiwsw)

   ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refrwlw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refiwlw)

   ! set complex representation of refractive indices as module data
   do i = 1, nswbands
      crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)),kind=r8)
   end do
   do i = 1, nlwbands
      crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)),kind=r8)
   end do

   call pio_closefile(ncid)

end subroutine read_water_refindex

!===============================================================================

end module coreshell_aer_opt
