module ndrop_carma

!---------------------------------------------------------------------------------
! Purpose:
!   CAM Interface for droplet activation by modal aerosols
!
! ***N.B.*** This module is currently hardcoded to recognize only the modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use ppgrid,           only: pcols, pver
use physconst,        only: pi, rhoh2o, mwh2o, r_universal, rh2o, &
                            gravit, latvap, cpair, rair
use constituents,     only: pcnst, cnst_get_ind, cnst_name, cnst_spec_class_gas, cnst_species_class
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

use wv_saturation,    only: qsat
use phys_control,     only: phys_getopts
use ref_pres,         only: top_lev => trop_cloud_top_lev
use shr_spfn_mod,     only: erf => shr_spfn_erf
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_info_by_bin, rad_cnst_get_info_by_bin_spec,       &
                            rad_cnst_get_bin_props_by_idx, rad_cnst_get_bin_num, rad_cnst_get_bin_mmr_by_idx, &
                            rad_cnst_get_bin_mmr
use cam_history,      only: addfld, add_default, horiz_only, fieldname_len, outfld
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog
use aerosol_model_mod, only: aerosol_model, ptr2d_t

implicit none
private

public ndrop_carma_init

! description of bin aerosols
integer, public, protected :: nbins = 0
integer, public, protected, allocatable :: nspec(:)

logical :: history_aerosol      ! Output the aerosol tendencies
character(len=fieldname_len), allocatable :: fieldname(:)    ! names for drop nuc tendency output fields
character(len=fieldname_len), allocatable :: fieldname_cw(:) ! names for drop nuc tendency output fields

! local indexing for bins
integer :: ncnst_tot                  ! total number of mode number conc + mode species

! modal aerosols
logical :: prog_modal_aero     ! true when aerosols are prognostic   !st make sure to check

!===============================================================================
contains
!===============================================================================

subroutine ndrop_carma_init

   integer  :: ii, l, m, mm
   integer :: idxtmp    = -1
   character(len=32)   :: tmpname
   character(len=32)   :: tmpname_cw
   character(len=128)  :: long_name
   character(len=8)    :: unit
   logical :: history_amwg         ! output the variables used by the AMWG diag package

   !-------------------------------------------------------------------------------

   ! get info about the modal aerosols
   ! get nbins

   call rad_cnst_get_info( 0, nbins=nbins)

   allocate( nspec(nbins) )

   do m = 1, nbins
      call rad_cnst_get_info_by_bin(0, m, nspec=nspec(m))
   end do

   ! Init the table for local indexing of bin number conc and mmr.
   ! This table uses species index 0 for the number conc.

   ! Find max number of species in all the bins, and the total
   ! number of bin number concentrations + mode species

   ncnst_tot = nspec(1) + 2
   do m = 2, nbins
      ncnst_tot = ncnst_tot + nspec(m) + 2
   end do

   allocate( fieldname(ncnst_tot), &
             fieldname_cw(ncnst_tot) )

   ! Add dropmixnuc tendencies for all modal aerosol species

   call phys_getopts(history_amwg_out = history_amwg, &
                     history_aerosol_out = history_aerosol)
   prog_modal_aero = .true.

   mm = 0
   do m = 1, nbins
      do l = 0, nspec(m) + 1  ! loop over bin + aerosol constituents

         mm = mm + 1

         unit = 'kg/m2/s'
         if (l == 0) then   ! number
            unit = '#/m2/s'
         end if

         if (l == 0) then   ! number
            call rad_cnst_get_info_by_bin(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
         else if (l == 1) then
            call rad_cnst_get_info_by_bin(0, m,  mmr_name=tmpname, mmr_name_cw=tmpname_cw)
         else
            call rad_cnst_get_info_by_bin_spec(0, m, l-1, spec_name=tmpname, spec_name_cw=tmpname_cw)
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

end subroutine ndrop_carma_init

end module ndrop_carma
