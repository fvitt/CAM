module edyn3d_esmf_phys_mesh_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc

  use ESMF

  implicit none

  type(ESMF_Mesh)  :: phys_mesh

contains

  subroutine edyn3d_esmf_phys_mesh_init
    use phys_control, only: phys_getopts
    use ppgrid, only: begchunk, endchunk
    use phys_grid, only: get_ncols_p, get_gcol_p

    integer :: total_cols, rc
    integer :: ncols, chnk, col, dindex
    character(len=cl) :: mesh_file
    integer,allocatable :: decomp(:)
    type(ESMF_DistGrid) :: dist_grid_2d

    character(len=*), parameter :: subname = 'edyn3d_esmf_phys_mesh_init'

    ! physics grid / field


    call phys_getopts(physics_grid_out=mesh_file)

    ! Compute the local decomp
    total_cols = 0
    do chnk = begchunk, endchunk
       total_cols = total_cols + get_ncols_p(chnk)
    end do
    allocate(decomp(total_cols))
    dindex = 0
    do chnk = begchunk, endchunk
       ncols = get_ncols_p(chnk)
       do col = 1, ncols
          dindex = dindex + 1
          decomp(dindex) = get_gcol_p(chnk, col)
       end do
    end do

    ! Create a DistGrid based on the physics decomp
    dist_grid_2d = ESMF_DistGridCreate(arbSeqIndexList=decomp, rc=rc)
    call check_error(subname,'ESMF_DistGridCreate phys decomp dist_grid_2d',rc)

    phys_mesh = ESMF_MeshCreate(trim(mesh_file), ESMF_FILEFORMAT_ESMFMESH,  &
                                elementDistgrid=dist_grid_2d, rc=rc)
    call check_error(subname,'ESMF_MeshCreate phys_mesh',rc)

    deallocate(decomp)

  end subroutine edyn3d_esmf_phys_mesh_init

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine check_error(subname, routine, rc)

    character(len=*), intent(in) :: subname
    character(len=*), intent(in) :: routine
    integer,          intent(in) :: rc

    character(len=cl) :: errmsg

    if (rc /= ESMF_SUCCESS) then
       write(errmsg, '(4a,i0)') trim(subname), ': Error return from ', trim(routine), ', rc = ', rc
       if (masterproc) then
          write(iulog, '(2a)') 'ERROR: ', trim(errmsg)
       end if
       call endrun(trim(errmsg))
    end if
  end subroutine check_error

end module edyn3d_esmf_phys_mesh_mod
