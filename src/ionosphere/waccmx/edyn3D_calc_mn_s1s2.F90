!
     subroutine edyn3D_calc_mn_s1s2

       use edyn3D_fieldline, only: fieldline_s1,fline_s1,fieldline_s2,fline_s2,&
                                   fieldline_r,fline_r,fieldline_p,fline_p
       use edyn3D_params,    only: m1f,m2f,m3f,nmlat_h,nmlatS2_h,nhgt_fix,hgt_fix,&
                               hgt_fix_r,nhgt_fix_r,r0,ylonm,ylonm_s,rho,rho_s,pi,ylatm     
       use edyn3d_mpi,       only: mlon0_p,mlon1_p
       use shr_kind_mod,     only: r8 => shr_kind_r8            ! 8-byte reals
       use cam_logfile,      only: iulog
       use spmd_utils,       only: masterproc

       implicit none
       !       
       ! Local variables
       !         
       integer ::  i,j,k,isn,nmax      
       real(r8) :: dlonm,drho,sigC
!       
       dlonm = ylonm_s(2)-ylonm_s(1) ! Assumes equidistant longitudinal gridpoints
!       
       do isn = 1,2 ! loop over both hemispheres
         do i=mlon0_p,mlon1_p ! loop over task longitudes
!	 
! calculate values for S1-points these are at (i+/-0.5,j,k) updated 8/23/2015
!          	 
           do j=1,nmlat_h ! loop over all latitudes in one hemisphere 
	      nmax = fline_s1(i,j,isn)%npts ! maximum of points on fieldline
              do k=1,nmax
!	        	
		  fline_s1(i,j,isn)%M1(k) = m1f(j,k)/fline_s1(i,j,isn)%F(k)
		  ! if(i.eq.10) write(92,*) j,k,isn,fline_s1(i,j,isn)%M1(k)
		  
		  if(j > 1) then ! do not calculate N1p & N1h for the pole (not used)
! N1P(i+0.5) = M1(i+0.5)*[sigP*d1^2](i+0.5)/R/rho(j)/(phi(i+1)-phi(i))		  
		    fline_s1(i,j,isn)%N1p(k) = fline_s1(i,j,isn)%M1(k)*fline_s1(i,j,isn)%sigP(k)* &
		       fline_s1(i,j,isn)%d1d1(k)/r0/rho(j,isn)/dlonm
	   
! N1H(i+0.5) = M1(i+0.5)*[sigH*D-sigP*d1*d2](i+0.5)*sqrt(1-0.75*rho^2(j))/2/R/(rho(j+1)-rho(j-1))
		    if(j.eq.nmlat_h) then
		       drho = 2._r8*(rho(j,isn)-rho(j-1,isn))
		    else
		       drho = rho(j+1,isn)-rho(j-1,isn)
		    endif
		    fline_s1(i,j,isn)%N1h(k) = fline_s1(i,j,isn)%M1(k)*(fline_s1(i,j,isn)%sigH(k)* &
		       fline_s1(i,j,isn)%D(k)-fline_s1(i,j,isn)%sigP(k)*fline_s1(i,j,isn)%d1d2(k))* &
		       sqrt(1-0.75_r8*rho(j,isn)**2)*0.5_r8/r0/drho
		   endif
		   
! lowest volume at equator : overwrite values from above N1P and N1H (page 12 30 Jan 2014 (Art's notes)
! N1H = 0
! N1P -> N1C
! N1C(i+0.5,j,k) = M1(i+0.5,j,k)*sigC(i+0.5,j,k)/R/rho(j)/(phi(i+1)-phi(i))
! sigC = sigP*d1*d1+(sigH*D-sigP*d1*d2)*(sigH*D+sigP*d1*d2)/sigP/(d2*d2)
!        for i+0.5,j,k		
	          if(k == 1 .and. (nmlat_h-j+1) == k) then ! lowest volume at equator j ==  nmlat_h
		    fline_s1(i,j,isn)%N1h(k) = 0.
                    sigC = fline_s1(i,j,isn)%sigP(k)*fline_s1(i,j,isn)%d1d1(k)
                    sigC =sigC + (fline_s1(i,j,isn)%sigH(k)*fline_s1(i,j,isn)%D(k)- &
		       fline_s1(i,j,isn)%sigP(k)*fline_s1(i,j,isn)%d1d2(k))* &
		       (fline_s1(i,j,isn)%sigH(k)*fline_s1(i,j,isn)%D(k)+ &
		       fline_s1(i,j,isn)%sigP(k)*fline_s1(i,j,isn)%d1d2(k))/ &
                       fline_s1(i,j,isn)%sigP(k)/fline_s1(i,j,isn)%d2d2(k)
		    fline_s1(i,j,isn)%N1p(k) = fline_s1(i,j,isn)%M1(k)*sigC/r0/rho(j,isn)/dlonm
                  endif ! end lowest equatorial volume
	  
               end do  ! end height loop
             end do  ! end lat/fieldline loop
	 
! calculate values for S2-points these are at (i,j+0.5,k)
!          		     
	     do j=1,nmlatS2_h                ! loop over all latitudes in one hemisphere  from pole towards equator
	       nmax = fline_s2(i,j,isn)%npts ! maximum of points on fieldline
               do k=1,nmax                   ! loop over all heights
!
! there is no S2 point at j+0.5 therefore no conditional statement is needed
! S2 is inbetween p-points in horizontal, but not at the pole or the equator
		  fline_s2(i,j,isn)%M2(k) = m2f(j,k)/fline_s2(i,j,isn)%F(k)
		  !if(i.eq.10) write(93,*) j,k,isn,fline_s2(i,j,isn)%M2(k)
!		   
! N2H(j+0.5) = M2(j+0.5)*[sigH*D+sigP*d1*d2]](j+0.5)/2/R/rho(j+0.5)/(phi(i+1)-phi(i-1)))		  
		  fline_s2(i,j,isn)%N2h(k) = fline_s2(i,j,isn)%M2(k)*(fline_s2(i,j,isn)%sigH(k)* &
		     fline_s2(i,j,isn)%D(k)+fline_s2(i,j,isn)%sigP(k)*fline_s2(i,j,isn)%d1d2(k))*0.5/ &
		     r0/2/rho_s(j,isn)/dlonm
! N2P(j+0.5) = M2(j+0.5)[sigP*d2^2](j+0.5)*sqrt(1-0.75 rho^2(j+0.5)/R/(rho(j+1)-rho(j))  
		  fline_s2(i,j,isn)%N2p(k) = fline_s2(i,j,isn)%M2(k)*fline_s2(i,j,isn)%sigP(k)* &
		     fline_s2(i,j,isn)%d2d2(k)*sqrt(1-0.75_r8*rho_s(j,isn)**2)/r0/(rho(j+1,isn)-rho(j,isn))
!	  
               end do  ! end height loop
             end do  ! end lat/fieldline loop
	     
! calculate values for R-points these are at (i,j,k+/-1/2) updated 8/23/2015 -> eq (64') Art's notes page 8
!
!   calculate at the pole since needed for mapping to QD coordinates at poles rho(j) = 0
!		     
	     do j=1,nmlat_h   ! loop over all latitudes in one hemisphere 
	       nmax = fline_r(i,j,isn)%npts ! maximum of points on fieldline
               do k=1,nmax  
!write(iulog,*) 'edyn3D_calc_mn_s1s2: isn,i,j,k,fline_r(i,j,isn)%F(k) ', isn,i,j,k,fline_r(i,j,isn)%F(k)
!                 if (fline_r(i,j,isn)%F(k) == 0._r8) fline_r(i,j,isn)%F(k) = fline_r(i,j,isn)%F(k-1) 
		 fline_r(i,j,isn)%M3(k) = m3f(j,k)/fline_r(i,j,isn)%F(k)
		 ! if(i.eq.10) write(94,*) j,k,isn,fline_r(i,j,isn)%M3(k)
 	         
               end do  ! end height loop
             end do  ! end lat/fieldline loop
       
! calculate values for P-points these are at (i,j,k) needed for Jr calculation Art's notes eq. (64') page 8 (updated 8/3/2015)
! M3(i,j,k) = r(k)^2*(phi(i+0.5)-phi(i-0.5))*sqrt(1-r(k)/R*rho(j)^2)*[sqrt(1-r(k)/R*rho(j-0.5)^2)-
!             sqrt(1-r(k)/R*rho(j+0.5)^2)]/F(i,jk-0.5)
! do not calculate at the pole
!		     
	     do j=2,nmlat_h   ! loop over all latitudes in one hemisphere but not pole
	       nmax = fline_p(i,j,isn)%npts ! maximum of points on fieldline
               do k=1,nmax  
		   fline_p(i,j,isn)%M3(k) = m3f(j,k)/fline_p(i,j,isn)%F(k)
               end do  ! end height loop
             end do  ! end lat/fieldline loop
	     
          end do  ! end longitude loop
	  
       end do ! end hemisphere loop

     end subroutine edyn3D_calc_mn_s1s2
