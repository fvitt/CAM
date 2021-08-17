module ndrop

!---------------------------------------------------------------------------------
! Purpose:
!   CAM Interface for droplet activation by modal aerosols
!
! ***N.B.*** This module is currently hardcoded to recognize only the modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use physconst,        only: pi, rhoh2o, mwh2o, r_universal, rh2o, &
                            gravit, latvap, cpair, epsilo, rair
use constituents,     only: pcnst, cnst_get_ind, cnst_name, cnst_spec_class_gas, cnst_species_class
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

use wv_saturation,    only: qsat
use phys_control,     only: phys_getopts
use ref_pres,         only: top_lev => trop_cloud_top_lev
use shr_spfn_mod,     only: erf => shr_spfn_erf
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr, &
                            rad_cnst_get_aer_props, rad_cnst_get_mode_props,                &
                            rad_cnst_get_mam_mmr_idx, rad_cnst_get_mode_num_idx
use cam_history,      only: addfld, add_default, horiz_only, fieldname_len, outfld
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog

use aerosol_model_mod, only: aerosol_model, ptr2d_t

implicit none
private

public ndrop_init

! description of modal aerosols
integer               :: ntot_amode     ! number of aerosol modes
integer,  allocatable :: nspec_amode(:) ! number of chemical species in each aerosol mode

logical :: history_aerosol      ! Output the MAM aerosol tendencies
character(len=fieldname_len), allocatable :: fieldname(:)    ! names for drop nuc tendency output fields
character(len=fieldname_len), allocatable :: fieldname_cw(:) ! names for drop nuc tendency output fields

! local indexing for MAM
integer :: ncnst_tot                  ! total number of mode number conc + mode species

! modal aerosols
logical :: prog_modal_aero     ! true when modal aerosols are prognostic

!===============================================================================
contains
!===============================================================================

subroutine ndrop_init

   integer  :: ii, l, lptr, m, mm
   character(len=32)   :: tmpname
   character(len=32)   :: tmpname_cw
   character(len=128)  :: long_name
   character(len=8)    :: unit
   logical :: history_amwg         ! output the variables used by the AMWG diag package

   !-------------------------------------------------------------------------------

   ! get info about the modal aerosols
   ! get ntot_amode
   call rad_cnst_get_info(0, nmodes=ntot_amode)

   allocate( nspec_amode(ntot_amode) )

   do m = 1, ntot_amode
      call rad_cnst_get_info(0, m, nspec=nspec_amode(m))
   end do

   ! Init the table for local indexing of mam number conc and mmr.
   ! This table uses species index 0 for the number conc.

   ! Find max number of species in all the modes, and the total
   ! number of mode number concentrations + mode species

   ncnst_tot = nspec_amode(1) + 1
   do m = 2, ntot_amode
      ncnst_tot = ncnst_tot + nspec_amode(m) + 1
   end do

   allocate( fieldname(ncnst_tot), &
             fieldname_cw(ncnst_tot) )

   ! Add dropmixnuc tendencies for all modal aerosol species

   call phys_getopts(history_amwg_out = history_amwg, &
                     history_aerosol_out = history_aerosol, &
                     prog_modal_aero_out=prog_modal_aero)


   mm = 0
   do m = 1, ntot_amode
      do l = 0, nspec_amode(m)   ! loop over number + chem constituents

         mm = mm + 1

         unit = 'kg/m2/s'
         if (l == 0) then   ! number
            unit = '#/m2/s'
         end if

         if (l == 0) then   ! number
            call rad_cnst_get_info(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
         else
            call rad_cnst_get_info(0, m, l, spec_name=tmpname, spec_name_cw=tmpname_cw)
         end if

         fieldname(mm)    = trim(tmpname) // '_mixnuc1'
         fieldname_cw(mm) = trim(tmpname_cw) // '_mixnuc1'

         if (prog_modal_aero) then

            ! Add tendency fields to the history only when prognostic MAM is enabled.
            long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency'
            call addfld(fieldname(mm),    horiz_only, 'A', unit, long_name)

            long_name = trim(tmpname_cw) // ' dropmixnuc mixnuc column tendency'
            call addfld(fieldname_cw(mm), horiz_only, 'A', unit, long_name)

            if (history_aerosol) then
               call add_default(fieldname(mm), 1, ' ')
               call add_default(fieldname_cw(mm), 1, ' ')
            end if

         end if

      end do
   end do

   call addfld('CCN1',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.02%')
   call addfld('CCN2',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.05%')
   call addfld('CCN3',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.1%')
   call addfld('CCN4',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.2%')
   call addfld('CCN5',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.5%')
   call addfld('CCN6',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=1.0%')

   call addfld('WTKE',     (/ 'lev' /), 'A', 'm/s', 'Standard deviation of updraft velocity')
   call addfld('NDROPMIX', (/ 'lev' /), 'A', '#/kg/s', 'Droplet number mixing')
   call addfld('NDROPSRC', (/ 'lev' /), 'A', '#/kg/s', 'Droplet number source')
   call addfld('NDROPSNK', (/ 'lev' /), 'A', '#/kg/s', 'Droplet number loss by microphysics')
   call addfld('NDROPCOL', horiz_only,  'A', '#/m2', 'Column droplet number')

   ! set the add_default fields
   if (history_amwg) then
      call add_default('CCN3', 1, ' ')
   endif

end subroutine ndrop_init

end module ndrop
