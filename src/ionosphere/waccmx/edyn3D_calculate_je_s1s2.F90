
       subroutine edyn3D_calculate_je_s1s2(fline_s1,fline_s2)
!
       use edyn3D_fieldline,only: fieldline_s1,fieldline_s2
       use edyn3D_params, only: Je2Ion_eq
       use edyn3D_params, only: nmlon,nmlat_h,nmlatS2_h,nhgt_fix,hgt_fix,&
                                hgt_fix_r,r0,ylonm_s,rho,rho_s,rearth_m,pi,ylatm,J3LB, &
			        test_pot,Jpg_add,no_wind,nhgt_fix_r,use_lbJ
       use edyn3d_mpi,    only: mlon0_p,mlon1_p
       use shr_kind_mod,  only: r8 => shr_kind_r8            ! 8-byte reals
       use cam_logfile,   only: iulog
       use spmd_utils,    only: masterproc
       !
       implicit none
       !
       type(fieldline_s1),dimension(mlon0_p-1:mlon1_p+1,nmlat_h,2),intent(inout) :: fline_s1
       type(fieldline_s2),dimension(mlon0_p-1:mlon1_p+1,nmlats2_h,2),intent(inout) :: fline_s2
       !
       integer :: isn,i,j,k,nmax,kmax_r
       real(r8):: ue1,ue2,fac
       logical :: test_pot_loc

! diagnostic for  Jf1Dyn &  Jf2Dyn
       real(r8) :: bhat(3),jpar,je2_woD

       test_pot_loc = test_pot
       if(.not.use_lbJ) J3LB=0._r8

        do isn = 1,2 ! loop over both hemisphere
!          do i=1,nmlon ! loop over all longitudes
          do i=mlon0_p,mlon1_p ! loop over task longitudes
!
! at i+0.5,j,k calculate
! Je1D = sigP*d1*d1*ue2*Be3+(sigH*D-sigP*d1*d2)*ue1*Be3+Je1^Ion
! calculate values for S1-points these are at (i+0.5,j,k)
!
           do j=2,nmlat_h ! loop over all latitudes in one hemisphere NOT THE POLE
              nmax = fline_s1(i,j,isn)%npts ! maximum of points on fieldline
              do k=1,nmax
                if(no_wind) then  ! zero winds
                    ue1 = 0._r8
                    ue2 = 0._r8
                else
                  if(test_pot_loc) then
                    ue1 = fline_s1(i,j,isn)%un(k)
                    ue2 = fline_s1(i,j,isn)%vn(k)
                  else
              ! ue1 = (un,vn)dot d1
                    ue1 = fline_s1(i,j,isn)%un(k)*fline_s1(i,j,isn)%d1(1,k)+ &
                     fline_s1(i,j,isn)%vn(k)*fline_s1(i,j,isn)%d1(2,k)
              ! ue2 = (un,vn)dot d2
                    ue2 = fline_s1(i,j,isn)%un(k)*fline_s1(i,j,isn)%d2(1,k)+ &
                     fline_s1(i,j,isn)%vn(k)*fline_s1(i,j,isn)%d2(2,k)
                  endif
                endif
!
                fac = fline_s1(i,j,isn)%sigP(k)*fline_s1(i,j,isn)%d1d1(k)*ue2* &
                  fline_s1(i,j,isn)%Be3(k)
                fac = fac + (fline_s1(i,j,isn)%sigH(k)* &
                  fline_s1(i,j,isn)%D(k)-fline_s1(i,j,isn)%sigP(k)*fline_s1(i,j,isn)%d1d2(k))* &
                  ue1*fline_s1(i,j,isn)%Be3(k)
                ! !!!!! COMMENT OUT ONCE Jpg works!!!!!!!
                if(Jpg_add) then
                  fac = fac + fline_s1(i,j,isn)%Je1Ion(k)
                end if

                fline_s1(i,j,isn)%Je1D(k)= fac

              end do  ! end height loop
            end do  ! end lat/fieldline loop
!
! at i,j+0.5,k calculate
! Je2 = (sigH*D+sigP*d1*d2)*ue2*Be3-sigP*d2*d2*ue1*Be3+Je1^Ion
! calculate values for S2-points these are at (i,j+0.5,k)
!
             do j=1,nmlatS2_h ! loop over all latitudes in one hemisphere
               nmax = fline_s2(i,j,isn)%npts ! maximum of points on fieldline
               do k=1,nmax
                 if(no_wind) then  ! zero winds
                    ue1 = 0._r8
                    ue2 = 0._r8
                 else
                   if(test_pot_loc) then
                    ue1 = fline_s2(i,j,isn)%un(k)
                    ue2 = fline_s2(i,j,isn)%vn(k)
                   else
              !  ue1 = (un,vn)dot d1
                    ue1 = fline_s2(i,j,isn)%un(k)*fline_s2(i,j,isn)%d1(1,k)+ &
                     fline_s2(i,j,isn)%vn(k)*fline_s2(i,j,isn)%d1(2,k)
              ! ue2 = (un,vn)dot d2
                    ue2 = fline_s2(i,j,isn)%un(k)*fline_s2(i,j,isn)%d2(1,k)+ &
                     fline_s2(i,j,isn)%vn(k)*fline_s2(i,j,isn)%d2(2,k)
                   endif
                 endif
!
                fac = fline_s2(i,j,isn)%sigH(k)*fline_s2(i,j,isn)%D(k)+ &
                  fline_s2(i,j,isn)%sigP(k)*fline_s2(i,j,isn)%d1d2(k)
                fac = fac*ue2*fline_s2(i,j,isn)%Be3(k)

                fac = fac - fline_s2(i,j,isn)%sigP(k) *fline_s2(i,j,isn)%d2d2(k)* &
                  ue1*fline_s2(i,j,isn)%Be3(k)

                if(Jpg_add) then
                  fac = fac + fline_s2(i,j,isn)%Je2Ion(k)
                end if

                fline_s2(i,j,isn)%Je2D(k)= fac

               end do  ! end height loop
             end do  ! end lat/fieldline loop
! lowest equatorial volume at i+0.5,j,k calculate
! Je2 = (sigH*D+sigP*d1*d2)*ue2*Be3-sigP*d2*d2*ue1*Be3+Je1^Ion
! calculate for lowest equatorial volume
! Je1D -> Je1S
! with   Je1S = Je1D -(sigH*D-sigP(d1*d2))/sigP/(d2*d2)*(Je2LB-Je2D) at (i+0.5,j,k)
!      Je2LB is the Je2 given by e.g. through coupling with the lower atmosphere
!         we assume that this is  Je2LB(i)= -J3LB(i,nmlat_h) not exactly since half a height level is
!         between but should be close
! could do for only one hemisphere (since the same point) and then copy into other hemisphere
               k = 1
               j = nmlat_h
               if(no_wind) then  ! zero winds
                  ue1 = 0._r8
                  ue2 = 0._r8
               else
                 if(test_pot_loc) then
                   ue1 = fline_s1(i,j,isn)%un(k)
                   ue2 = fline_s1(i,j,isn)%vn(k)
                  else
            ! ue1 = (un,vn)dot d1
                   ue1 = fline_s1(i,j,isn)%un(k)*fline_s1(i,j,isn)%d1(1,k)+ &
                        fline_s1(i,j,isn)%vn(k)*fline_s1(i,j,isn)%d1(2,k)
            ! ue2 = (un,vn)dot d2
                   ue2 = fline_s1(i,j,isn)%un(k)*fline_s1(i,j,isn)%d2(1,k)+ &
                        fline_s1(i,j,isn)%vn(k)*fline_s1(i,j,isn)%d2(2,k)
                  endif
                endif
!
! Je2D at the equator (there is no S2 point therefore needs to be calculated)
               fac = fline_s1(i,j,isn)%sigH(k)*fline_s1(i,j,isn)%D(k)+ &
                  fline_s1(i,j,isn)%sigP(k)*fline_s1(i,j,isn)%d1d2(k)
               fac = fac*ue2*fline_s1(i,j,isn)%Be3(k)
               fac = fac - fline_s1(i,j,isn)%sigP(k) *fline_s1(i,j,isn)%d2d2(k)* &
                  ue1*fline_s1(i,j,isn)%Be3(k)
               fac = fac +  Je2Ion_eq(i) ! Je2D at i+0.5,j,k

!  Je1S = Je1D -(sigH*D-sigP(d1*d2))/sigP/(d2*d2)*(Je2LB-Je2D) at (i+0.5,j,k)
!         we assume that this is  Je2LB(i)= -J3LB(i,nmlat_h)
               fline_s1(i,j,isn)%Je1D(k)= fline_s1(i,j,isn)%Je1D(k) - &
                (fline_s1(i,j,isn)%sigH(k)*fline_s1(i,j,isn)%D(k)-&
                  fline_s1(i,j,isn)%sigP(k)*fline_s1(i,j,isn)%d1d2(k))/ &
                  fline_s1(i,j,isn)%sigP(k)*fline_s1(i,j,isn)%d2d2(k)* &
                  (-J3LB(i,j,isn) - fac)

!
          end do  ! end longitude loop
       end do ! end hemisphere loop

       end subroutine edyn3D_calculate_je_s1s2
