module edyn3D_hemislv_mod

  implicit none

contains

  subroutine edyn3D_hemislv_poten
    use edyn3D_calc_coef_fac_const_rhs, only: edyn3D_calc_coef,edyn3D_calc_FAC,edyn3D_add_coef_ns, &
         edyn3D_gather_coef_ns,edyn3D_const_rhs,edyn3D_scatter_poten, &
         edyn3D_solve_sparse
    use edyn3D_mpi, only: mytid

    call edyn3D_calc_coef       ! - calc LHS & RHS

    call edyn3D_calc_FAC        ! - calc high latitude

    call edyn3D_add_coef_ns     ! - add North & South coef

    call edyn3D_gather_coef_ns  ! - gather coef_ns for solver

    if (mytid == 0) then
#ifdef HAS_SUPERLU_SLV
       call edyn3D_solve_sparse
#else
       call edyn3D_const_rhs     ! - solver - solve for rhs (electric potential)
#endif
    endif

    call edyn3D_scatter_poten   ! - Send global potential to each task

  end subroutine edyn3D_hemislv_poten

end module edyn3D_hemislv_mod
