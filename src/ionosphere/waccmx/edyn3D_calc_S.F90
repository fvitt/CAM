!
       subroutine edyn3D_calc_S
!
! S is the wind driven and ionospheric current sources (89) 
! I3^E are the external current sources from the magnetosphere I3^M and from the lower atmosphere
!    see eq. (108) I3^E = -I3^M(NH+SH) + I3_lb(NH) + I3_lb(SH) 
! eq (115) I3_total = Sum_k[S(NH) + S(SH)]  + I3^E     
!        
       use edyn3D_fieldline, only: fieldline_s1,fline_s1,fieldline_s2,fline_s2, &
                                 fieldline_p,fline_p,fieldline_r,fline_r
       use edyn3D_params,    only: nmlon,nmlat_h,val_fill,rtd,nhgt_fix, &
                               hgt_fix,ylonm,ylatm,nmlatS2_h,J3LB,use_stabil				  
       use edyn3d_mpi,       only: mlon0_p,mlon1_p
       use shr_kind_mod,     only: r8 => shr_kind_r8            ! 8-byte reals
       use cam_logfile,      only: iulog
       use spmd_utils,       only: masterproc
!
       implicit none
!
       integer :: isn,i,j,k,nmax,isn1,isn2,im
       real(r8) :: tmp1, tmp2, sum1(nhgt_fix)
!      
! page 12 equation (89) Art's script
! -S(i,j,k) = M1(i-0.5,j,k)*Je1D(i-0.5,j,k)-M1(i+0.5,j,k)*Je1D(i+0.5,j,k) +
!       M2(i,j-0.5,k)*J2D(i,j-0.5,k)-M2(i,j+0.5,k)*Je1D(i,j+0.5,k)
! the -S(i,j,k) denotes that it is on the left hand side, but for the
!   dynamo equation we need it on the right hand side, therefore
!   multiply by -1
! for equatorial volumes at different heights M2(i,j+0.5,k) = 0
!
! the following is taken care of in cal_je_s1s2
! for lowest equatorial volume S = M1(i-0.5,j,k)*Je1S(i-0.5,j,k)-
!       M1(i+0.5,j,k)*Je1S(i+0.5,j,k) +
!       M2(i,j-0.5,k)*J2D(i,j-0.5,k)
!   with Je1S = Je1D -(sigH*D-sigP(d1*d2)/sigP/(d2*d2)*(Je2LB-Je2D) at (i+0.5,j,k)
!      Je2LB is the Je2 given by e.g. through coupling with the lower atmosphere
!         at the moment is is set to zero
!    
!       write(6,*) 'what is the height',hgt_fix
!       
       do isn = 1,2  ! loop over hemisphere
         do i=mlon0_p,mlon1_p ! loop over task longitudes
!         do i=1,nmlon ! loop over all longitudes 

!         do i=1,nmlon ! loop over all longitudes 

            im = i-1

!           if(i == mlon0_p) then
!             im = mlon1_p ! use wrap around point as i-1 -> nmlon 
!           else
!             im = i-1
!           endif

!           if(i == 1) then
!             im = nmlon ! use wrap around point as i-1 -> nmlon 
!           else
!             im = i-1
!           endif
          fline_p(i,1,isn)%S(:) = val_fill  ! put fill value in the polar flux tubes
!              
          sum1 = 0.
          do j=2,nmlat_h  ! loop over all latitudes in one hemisphere; not pole value
             nmax = fline_p(i,j,isn)%npts ! maximum of points on fieldline
             !
             do k=1,nmax
             !
                if((nmlat_h-j+1) == k) then ! top volume at equator M2(i,j+0.5,k) = 0
                  fline_p(i,j,isn)%S(k) = fline_s1(im,j,isn)%M1(k)*fline_s1(im,j,isn)%Je1D(k)- &
                    fline_s1(i,j,isn)%M1(k)*fline_s1(i,j,isn)%Je1D(k)+ &
                    fline_s2(i,j-1,isn)%M2(k)*fline_s2(i,j-1,isn)%Je2D(k)
                  if(use_stabil) fline_p(i,j,isn)%S(k) = fline_p(i,j,isn)%S(k) + &       ! considering cowling conductivity 
                     fline_s1(im,j,isn)%M1Je1D(k) - fline_s1(i,j,isn)%M1Je1D(k)
                  fline_p(i,j,isn)%S(k) = -fline_p(i,j,isn)%S(k)
                  
                  tmp1 =fline_s1(im,j,isn)%M1(k)*fline_s1(im,j,isn)%Je1D(k)- &
                  fline_s1(i,j,isn)%M1(k)*fline_s1(i,j,isn)%Je1D(k)
                  tmp2 =fline_s2(i,j-1,isn)%M2(k)*fline_s2(i,j-1,isn)%Je2D(k)
                  
!               write(88,'(2(f10.5,x),5(x,e17.10))') fline_p(i,j,isn)%mlon_m, &
!               fline_p(i,j,isn)%mlat_m,fline_p(i,j,isn)%hgt_pt(k), &
!               fline_s1(i,j,isn)%M1(k),fline_s2(i,j-1,isn)%M2(k) , &
!               fline_s1(i,j,isn)%Je1D(k),fline_s2(i,j-1,isn)%Je2D(k)
               
               sum1(k) = sum1(k) + fline_s1(i,j,isn)%M1(k)
               
                else 

                  fline_p(i,j,isn)%S(k) = fline_s1(im,j,isn)%M1(k)*fline_s1(im,j,isn)%Je1D(k)- &
                    fline_s1(i,j,isn)%M1(k)*fline_s1(i,j,isn)%Je1D(k)+ &
                    fline_s2(i,j-1,isn)%M2(k)*fline_s2(i,j-1,isn)%Je2D(k)- &
                    fline_s2(i,j,isn)%M2(k)*fline_s2(i,j,isn)%Je2D(k)
                  if(use_stabil) fline_p(i,j,isn)%S(k) = fline_p(i,j,isn)%S(k) + &         ! considering cowling conductivity 
                     fline_s1(im,j,isn)%M1Je1D(k) - fline_s1(i,j,isn)%M1Je1D(k)
                  fline_p(i,j,isn)%S(k) = -fline_p(i,j,isn)%S(k)
                  
                  tmp1 = fline_s1(i,j,isn)%M1(k)*fline_s1(i,j,isn)%Je1D(k)
                  tmp2 = fline_s2(i,j,isn)%M2(k)*fline_s2(i,j,isn)%Je2D(k)
                  
!          write(88,'(2(f10.5,x),5(x,e17.10))') fline_p(i,j,isn)%mlon_m, &
!              fline_p(i,j,isn)%mlat_m,fline_p(i,j,isn)%hgt_pt(k),&
!              fline_s1(i,j,isn)%M1(k),fline_s2(i,j,isn)%M2(k) ,&
!              fline_s1(i,j,isn)%Je1D(k),fline_s2(i,j,isn)%Je2D(k)
               
               sum1(k) = sum1(k) + fline_s1(i,j,isn)%M1(k)
               
                endif
!                write(99,'(4(x,i4),,1(x,e17.10))') i,j,k,isn,fline_p(i,j,isn)%S(k)
             end do  ! end height loop
             
             
!
! add the current from the lower atmosphere  J3_lb*M3 at r point for k=0.5 index=1      
!   I3S_lb = J3LB*M3  
             fline_p(i,j,isn)%S(1) = fline_p(i,j,isn)%S(1) + J3LB(i,j,isn)*fline_r(i,j,isn)%M3(1)
           
           end do  ! end lat/fieldline loop
!         
         end do  ! end longitude loop
       end do ! end hemisphere loop          j=1 ! pole values

       end subroutine edyn3D_calc_S
