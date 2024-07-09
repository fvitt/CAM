module aero_deposition_cam

  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: cnst_get_ind, pcnst
  use camsrfexch,   only: cam_out_t
  use aerosol_properties_mod, only: aero_name_len
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private
  public :: aero_deposition_cam_init
  public :: aero_deposition_cam_setwet

  integer :: bcarbon_ndx( pcnst ) = -1
  integer :: bcarbon_cnt = 0
  integer :: ocarbon_ndx( pcnst ) = -1
  integer :: ocarbon_cnt = 0

  class(aerosol_properties), pointer :: aero_props=>null()
  integer :: nele_tot=0            ! total number of aerosol elements

  ! bulk dust bins (meters)

  integer, parameter :: n_bulk_dst_bins = 4
  real(r8), parameter :: bulk_dst_edges(n_bulk_dst_bins+1) = &
       (/1.e-15_r8, 1.25e-6_r8, 2.5e-6_r8, 5.0e-6_r8, 1.0e-2_r8/) ! meters
    ! ????
!!$    ! bulk dust sizes (0.05-0.5, 0.5-1.25, 1.25-2.5, and 2.5-5.0 microns) ????

contains

  !==============================================================================
  subroutine aero_deposition_cam_init(aero_props_in)

    class(aerosol_properties),target, intent(in) :: aero_props_in

    integer :: pcnt, scnt
    character(len=*), parameter :: subrname = 'aero_deposition_cam_init'

    aero_props => aero_props_in

    call get_indices( type='black-c',   indices=bcarbon_ndx, count=bcarbon_cnt )
    call get_indices( type='s-organic', indices=ocarbon_ndx, count=pcnt )
    call get_indices( type='p-organic', indices=ocarbon_ndx(pcnt+1:), count=scnt )
    ocarbon_cnt = pcnt+scnt


    nele_tot = aero_props%ncnst_tot()

  contains

    !==============================================================================
    subroutine get_indices( type, indices, count)

      character(len=*), intent(in) :: type
      integer, intent(out) :: indices(:)
      integer, intent(out) :: count

      integer :: ibin,ispc, ndx, nspec
      character(len=aero_name_len) :: spec_type, spec_name
      real(r8) :: minrad
      logical :: getndx

      count = 0
      indices(:) = -1

      do ibin = 1, aero_props%nbins()

         do ispc = 1, aero_props%nspecies(ibin)

            call aero_props%get(ibin,ispc, spectype=spec_type, specname=spec_name )

            if (spec_type==type) then

               getndx = .true.
               if (getndx) then
                  call cnst_get_ind(spec_name, ndx, abort=.false.)
                  if (ndx>0) then
                     count = count+1
                     indices(count) = ndx
                  endif
               endif

            endif

         enddo

      enddo

    end subroutine get_indices

  end subroutine aero_deposition_cam_init

  !==============================================================================
  subroutine aero_deposition_cam_setwet(aerdepwetis, aerdepwetcw, cam_out)


   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec, ibin, mm, ndx
   integer :: ncol                      ! number of columns

   real(r8) :: dep_fluxes(nele_tot)
   real(r8) :: dst_fluxes(n_bulk_dst_bins)
   character(len=aero_name_len) :: specname, name_c

   ncol = cam_out%ncol

   cam_out%bcphiwet(:) = 0._r8
   cam_out%ocphiwet(:) = 0._r8
   cam_out%dstwet1(:) = 0._r8
   cam_out%dstwet2(:) = 0._r8
   cam_out%dstwet3(:) = 0._r8
   cam_out%dstwet4(:) = 0._r8

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface,
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol

      ! black carbon fluxes
      do ispec=1,bcarbon_cnt
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) &
                             - (aerdepwetis(i,bcarbon_ndx(ispec))+aerdepwetcw(i,bcarbon_ndx(ispec)))
      enddo

      ! organic carbon fluxes
      do ispec=1,ocarbon_cnt
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) &
                             - (aerdepwetis(i,ocarbon_ndx(ispec))+aerdepwetcw(i,ocarbon_ndx(ispec)))
      enddo

      dep_fluxes = 0._r8
      dst_fluxes = 0._r8

      do ibin = 1,aero_props%nbins()
         do ispec = 0,aero_props%nspecies(ibin)
            if (ispec==0) then
               call aero_props%num_names(ibin, specname, name_c)
            else
               call aero_props%get(ibin,ispec, specname=specname)
            end if
            call cnst_get_ind(specname, ndx, abort=.false.)
            if (ndx>0) then
               mm = aero_props%indexer(ibin,ispec)
               dep_fluxes(mm) = - (aerdepwetis(i,ndx)+aerdepwetcw(i,ndx))
            end if
         end do
      end do

      ! rebin dust fluxes to bulk dust bins
      call aero_props%rebin_bulk_fluxes('dust', dep_fluxes, bulk_dst_edges, dst_fluxes)
      cam_out%dstwet1(i) = cam_out%dstwet1(i) + dst_fluxes(1)
      cam_out%dstwet2(i) = cam_out%dstwet2(i) + dst_fluxes(2)
      cam_out%dstwet3(i) = cam_out%dstwet3(i) + dst_fluxes(3)
      cam_out%dstwet4(i) = cam_out%dstwet4(i) + dst_fluxes(4)

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(i) .lt. 0._r8) cam_out%bcphiwet(i) = 0._r8
      if (cam_out%ocphiwet(i) .lt. 0._r8) cam_out%ocphiwet(i) = 0._r8
      if (cam_out%dstwet1(i)  .lt. 0._r8) cam_out%dstwet1(i)  = 0._r8
      if (cam_out%dstwet2(i)  .lt. 0._r8) cam_out%dstwet2(i)  = 0._r8
      if (cam_out%dstwet3(i)  .lt. 0._r8) cam_out%dstwet3(i)  = 0._r8
      if (cam_out%dstwet4(i)  .lt. 0._r8) cam_out%dstwet4(i)  = 0._r8

   enddo

  end subroutine aero_deposition_cam_setwet



end module aero_deposition_cam
