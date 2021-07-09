module carma_fixer_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use spmd_utils,     only: masterproc

  implicit none

contains

 subroutine carma_fix_pbuf( state, pbuf )
    use carma_flags_mod,   only: carma_model
    use rad_constituents,  only: rad_cnst_get_info, rad_cnst_get_info_by_bin
    use rad_constituents,  only: rad_cnst_get_bin_num, rad_cnst_get_info_by_bin_spec
    use rad_constituents,  only: rad_cnst_get_bin_mmr, rad_cnst_get_bin_mmr_by_idx
    use ppgrid,            only: pcols, pver
    use constituents,      only: pcnst, cnst_name, cnst_get_ind
    use physics_buffer,    only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
    use physics_types,     only: physics_state

    type(physics_state), intent(in) :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    character(len=32) :: name
    integer :: l,m, nspec, nbins, idxtmp, ncol, lbuf
    real(r8), pointer :: bin_mmr(:,:)
    real(r8), pointer :: bin_num(:,:)
    real(r8), pointer :: spec_mmr(:,:)
    real(r8), pointer :: rmass(:)
    real(r8) :: total_mmr(pcols,pver)

    logical, parameter :: debug=.false.
    character(len=*), parameter :: prefix = 'carma_fix_pbuf: '

    if (carma_model /= 'trop_strat') return

    ncol = state%ncol

    call rad_cnst_get_info(0, nbins=nbins)

    do m = 1, nbins
       call rad_cnst_get_bin_mmr(0, m, 'a', state, pbuf, bin_mmr)
       call rad_cnst_get_bin_num(0, m, 'a', state, pbuf, bin_num)
       call rad_cnst_get_info_by_bin(0, m, nspec=nspec, bin_name=name)
       call pbuf_get_field(pbuf, pbuf_get_index(trim(name)//"_rmass"),rmass)
       bin_num(:ncol,:) = bin_mmr(:ncol,:)/rmass(1)

       if (debug.and.masterproc) write(iulog,*) prefix//'bin mmr name = '//trim(name)
       lbuf = -1
       total_mmr(:,:) = 0._r8

       do l = 1, nspec
          call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=name)
          if (debug.and.masterproc) write(iulog,*) prefix//'  species mmr name = '//trim(name)
          call cnst_get_ind(name, idxtmp, abort=.false.)
          if ( idxtmp>0 ) then
             call rad_cnst_get_bin_mmr_by_idx(0, m, l, 'a', state, pbuf, spec_mmr)
             total_mmr(:ncol,:) = total_mmr(:ncol,:) + spec_mmr(:ncol,:)
          elseif (lbuf<1) then
             if (debug.and.masterproc) write(iulog,*) prefix//'  species '//trim(name)//' is in pbuf'
             lbuf = l
          else
             call endrun(prefix//' multiple bin species in pbuf')
          endif
       enddo
       if (lbuf>0) then
          call rad_cnst_get_info_by_bin_spec(0, m, lbuf, spec_name=name)
          if (debug.and.masterproc) write(iulog,*) prefix//'force '//trim(name)//' mmr = bin - total'
          call rad_cnst_get_bin_mmr_by_idx(0, m, lbuf, 'a', state, pbuf, spec_mmr)
          spec_mmr(:ncol,:) = max( bin_mmr(:ncol,:)-total_mmr(:ncol,:), 0._r8 )
       endif
    enddo

 end subroutine carma_fix_pbuf

end module carma_fixer_mod
