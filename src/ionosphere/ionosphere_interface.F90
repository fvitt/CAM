module ionosphere_interface

 ! Dummy interface -- actual ionosphere interface exist in src/ionosphere/modelname

  implicit none

  private

  public :: ionosphere_readnl
  public :: ionosphere_init
  public :: ionosphere_run1
  public :: ionosphere_run2
  public :: ionosphere_init_restart
  public :: ionosphere_write_restart
  public :: ionosphere_read_restart
  public :: ionosphere_final

contains

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_readnl( nlfile )

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  end subroutine ionosphere_readnl

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_init()
    use mag_grid_mod, only: mag_grid_mod_reg
    use cam_history, only: addfld

    call mag_grid_mod_reg()

    call addfld( 'TnPhysIn', (/ 'lev' /), 'I', 'K', 'Nuetral Temperature input on phys grid' )
    call addfld( 'TnPhysOut', (/ 'lev' /), 'I', 'K', 'Nuetral Temperature output on phys grid' )

  end subroutine ionosphere_init

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_run1(pbuf2d)
    use physics_buffer, only: physics_buffer_desc

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

  end subroutine ionosphere_run1

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_run2( phys_state, pbuf2d )
    use shr_kind_mod, only: r8 => shr_kind_r8, cl=>shr_kind_cl
    use ppgrid, only: begchunk, endchunk, pcols, pver
    use physics_types, only: physics_state
    use physics_buffer, only: physics_buffer_desc

    use mag_grid_mod, only: mag_grid_timestep
    use cam_history, only: outfld

    ! args
    type(physics_state),    intent(in) :: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: ncol, nphyscols, lchnk, astat
    integer :: i,j,k

    real(r8), pointer :: physalt(:,:)
    real(r8), pointer :: tn_in(:,:)
    real(r8), pointer :: tn_out(:,:)
    real(r8) :: phys_out(pcols,pver)

    nphyscols = 0
    do lchnk = begchunk, endchunk
       nphyscols = nphyscols + phys_state(lchnk)%ncol
    end do

    allocate(physalt(pver,nphyscols), stat=astat)
    allocate(tn_in(pver,nphyscols), stat=astat)
    allocate(tn_out(pver,nphyscols), stat=astat)

    j = 0
    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol
       call outfld( 'TnPhysIn', phys_state(lchnk)%t, pcols, lchnk )
       do i = 1, ncol
          j = j + 1
          do k = 1, pver
             physalt(k,j) = phys_state(lchnk)%zm(i,k) ! meters
             tn_in(k,j) = phys_state(lchnk)%t(i,k)
          end do
       end do
    end do

    call mag_grid_timestep( nphyscols, pver, physalt, tn_in, tn_out )

    j = 0
    do lchnk = begchunk, endchunk
       phys_out = -huge(1._r8)
       ncol = phys_state(lchnk)%ncol
       do i = 1, ncol
          j = j + 1
          do k = 1, pver
             phys_out(i,k) = tn_out(k,j)
          end do
       end do
       call outfld( 'TnPhysOut', phys_out, pcols, lchnk )
    end do

  end subroutine ionosphere_run2

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_init_restart(File)
    use pio, only: file_desc_t

    type(File_desc_t),  intent(inout) :: File

  end subroutine ionosphere_init_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_write_restart(File)
    use pio, only: file_desc_t

    type(File_desc_t), intent(inout) :: File

  end subroutine ionosphere_write_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_read_restart(File)
    use pio, only: file_desc_t

    type(file_desc_t), intent(inout) :: File

  end subroutine ionosphere_read_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_final

  end subroutine ionosphere_final

end module ionosphere_interface
