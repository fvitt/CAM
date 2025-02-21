
  subroutine edyn3D_calculate_efield(fline_p,fline_s1,fline_s2)
  !
  ! Calculates electric field Ed1,Ed2 and drift velocity ve1,ve2 at S1 and S2 points
  ! from potential and magnetic field
  ! ve1 = Ed2/Be3   &   ve2 = -Ed1/Be3

    use edyn3D_params,     only: nmlat_h,nmlatS2_h,ylonm,rho,rho_s,r0
    use edyn3D_mpi,        only: mlon0_p,mlon1_p
    use edyn3D_fieldline,  only: fieldline_p,fieldline_s1,fieldline_s2
    use shr_kind_mod,      only: r8 => shr_kind_r8            ! 8-byte reals
    use cam_logfile,       only: iulog
    use spmd_utils,        only: masterproc

    type(fieldline_p),dimension(mlon0_p-1:mlon1_p+1,nmlat_h,2),intent(inout) :: fline_p
    type(fieldline_s1),dimension(mlon0_p-1:mlon1_p+1,nmlat_h,2),intent(inout) :: fline_s1
    type(fieldline_s2),dimension(mlon0_p-1:mlon1_p+1,nmlatS2_h,2),intent(inout) :: fline_s2

    integer  :: i,j,isn
    real(r8) :: fac,facj
    logical  :: isclose

    fac = 1/(r0*(ylonm(2)-ylonm(1))) ! regular spaced in longitude

    do i = mlon0_p,mlon1_p
      do j = 2,nmlat_h-1 ! S1 loop not the pole and equator
        facj = sqrt(1-0.75_r8*rho(j,1)**2)/r0/2/(rho(j+1,1)-rho(j-1,1))
        do isn = 1,2
          fline_s1(i,j,isn)%ed1 = (fline_p(i,j,isn)%pot-fline_p(i+1,j,isn)%pot)*fac/rho(j,1)
          fline_s1(i,j,isn)%ed2 = facj* &
            (fline_p(i,j-1,isn)%pot+fline_p(i+1,j-1,isn)%pot- &
            fline_p(i,j+1,isn)%pot-fline_p(i+1,j+1,isn)%pot)
        enddo
      enddo

      j = 1 ! pole
      facj = sqrt(1-0.75_r8*rho(j,1)**2)/r0/2/(rho(j+1,1)-rho(j,1))
      do isn = 1,2
        fline_s1(i,j,isn)%ed1 = (fline_p(i,j,isn)%pot-fline_p(i+1,j,isn)%pot)*fac/rho(j+1,1)
        fline_s1(i,j,isn)%ed2 = facj* &
          (fline_p(i,j,isn)%pot+fline_p(i+1,j,isn)%pot- &
          fline_p(i,j+1,isn)%pot-fline_p(i+1,j+1,isn)%pot)
      enddo

      j = nmlat_h ! equator
      facj = sqrt(1-0.75_r8*rho(j,1)**2)/r0/2/(rho(j,1)-rho(j-1,1))
      do isn = 1,2
        fline_s1(i,j,isn)%ed1 = (fline_p(i,j,isn)%pot-fline_p(i+1,j,isn)%pot)*fac/rho(j,1)
        fline_s1(i,j,isn)%ed2 = facj* &
          (fline_p(i,j-1,isn)%pot+fline_p(i+1,j-1,isn)%pot- &
          fline_p(i,j,isn)%pot-fline_p(i+1,j,isn)%pot)
      enddo

      do j = 1,nmlat_h
        do isn = 1,2
            fline_s1(i,j,isn)%ve1 = fline_s1(i,j,isn)%ed2/fline_s1(i,j,isn)%be3(1)
            fline_s1(i,j,isn)%ve2 = -fline_s1(i,j,isn)%ed1/fline_s1(i,j,isn)%be3(1)
        enddo
      enddo

      do j = 1,nmlatS2_h ! S2 loop
        facj = sqrt(1-0.75_r8*rho_s(j,1)**2)/r0/(rho(j+1,1)-rho(j,1))
        do isn = 1,2
          fline_s2(i,j,isn)%ed1 = fac/4/rho_s(j,1)* &
            (fline_p(i-1,j,isn)%pot+fline_p(i-1,j+1,isn)%pot- &
            fline_p(i+1,j,isn)%pot-fline_p(i+1,j+1,isn)%pot)
          fline_s2(i,j,isn)%ed2 = (fline_p(i,j,isn)%pot-fline_p(i,j+1,isn)%pot)*facj
          fline_s2(i,j,isn)%ve1 = fline_s2(i,j,isn)%ed2/fline_s2(i,j,isn)%be3(1)
          fline_s2(i,j,isn)%ve2 = -fline_s2(i,j,isn)%ed1/fline_s2(i,j,isn)%be3(1)
        enddo
      enddo
    enddo

  endsubroutine edyn3D_calculate_efield
!-----------------------------------------------------------------------
