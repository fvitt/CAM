module zmean_phys_fields
  use shr_kind_mod, only: r8=>SHR_KIND_R8
  use ppgrid,       only: begchunk, endchunk, pcols, pver
  use physics_types,only: physics_state
  use cam_history,  only: addfld, outfld
  use infnan,       only: nan, assignment(=)

  use zmean_fields, only: zmean_fields_init, zmean_3d

  implicit none

contains

  subroutine zmean_phys_fields_init

    call zmean_fields_init

    call addfld ( 'Tfld', (/ 'lev' /), 'A', 'K', 'T field')
    call addfld ( 'Tzmn', (/ 'lev' /), 'A', 'K', 'T zonal mean')

  end subroutine zmean_phys_fields_init

  subroutine zmean_phys_fields_timestep_init(phys_state)
    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    real(r8) :: fld(pcols,pver,begchunk:endchunk)
    real(r8) :: zmfld(pcols,pver,begchunk:endchunk)

    real(r8) :: fld_tmp(pcols,pver)
    integer :: lchnk,ncol, i

    fld(:,:,:) = nan

    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol
       do i = 1,ncol
          fld(i,:,lchnk) = phys_state(lchnk)%t(i,:)
       end do
    end do

    zmfld = zmean_3d( fld )

    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol
       fld_tmp(:ncol,:) = fld(:ncol,:,lchnk)
       call outfld( 'Tfld', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = zmfld(:ncol,:,lchnk)
       call outfld( 'Tzmn', fld_tmp(:ncol,:), ncol, lchnk)
    end do

  end subroutine zmean_phys_fields_timestep_init


end module zmean_phys_fields
