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
    call addfld ( 'Ufld', (/ 'lev' /), 'A', 'K', 'T field')
    call addfld ( 'Uzmn', (/ 'lev' /), 'A', 'K', 'T zonal mean')
    call addfld ( 'Vfld', (/ 'lev' /), 'A', 'K', 'T field')
    call addfld ( 'Vzmn', (/ 'lev' /), 'A', 'K', 'T zonal mean')

  end subroutine zmean_phys_fields_init

  subroutine zmean_phys_fields_timestep_init(phys_state)
    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    real(r8) :: Tfld(pcols,pver,begchunk:endchunk)
    real(r8) :: Tzmfld(pcols,pver,begchunk:endchunk)
    real(r8) :: Ufld(pcols,pver,begchunk:endchunk)
    real(r8) :: Uzmfld(pcols,pver,begchunk:endchunk)
    real(r8) :: Vfld(pcols,pver,begchunk:endchunk)
    real(r8) :: Vzmfld(pcols,pver,begchunk:endchunk)

    real(r8) :: fld_tmp(pcols,pver)
    integer :: lchnk,ncol, i

    real(r8), parameter :: t0 = 200._r8

    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol
       do i = 1,ncol
          Tfld(i,:,lchnk) = phys_state(lchnk)%t(i,:) - t0
          Ufld(i,:,lchnk) = phys_state(lchnk)%u(i,:)
          Vfld(i,:,lchnk) = phys_state(lchnk)%v(i,:)
       end do
    end do

    Tzmfld = zmean_3d( Tfld ) + t0
    Uzmfld = zmean_3d( Ufld )
    Vzmfld = zmean_3d( Vfld )

    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol

       fld_tmp(:ncol,:) = Tfld(:ncol,:,lchnk) + t0
       call outfld( 'Tfld', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = Tzmfld(:ncol,:,lchnk)
       call outfld( 'Tzmn', fld_tmp(:ncol,:), ncol, lchnk)

       fld_tmp(:ncol,:) = Ufld(:ncol,:,lchnk)
       call outfld( 'Ufld', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = Uzmfld(:ncol,:,lchnk)
       call outfld( 'Uzmn', fld_tmp(:ncol,:), ncol, lchnk)

       fld_tmp(:ncol,:) = Vfld(:ncol,:,lchnk)
       call outfld( 'Vfld', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = Vzmfld(:ncol,:,lchnk)
       call outfld( 'Vzmn', fld_tmp(:ncol,:), ncol, lchnk)

    end do

  end subroutine zmean_phys_fields_timestep_init


end module zmean_phys_fields
