       
       subroutine edyn3D_get_conduct
       !
       !  Calculate Pedersen and Hall conductances for s1 and s2 field line grids
       !              
       use edyn3D_fieldline, only: fline_s1,fline_s2
       use edyn3D_params,    only: nmlon,nmlat_h,nmlatS2_h,nhgt_fix,hgt_fix
!                                   f107,ap,year,doy,ut,m2km,jpg     
       use edyn3d_mpi,       only: mlon0_p,mlon1_p
       use shr_kind_mod,     only: r8 => shr_kind_r8            ! 8-byte reals
       use cam_logfile,      only: iulog
       use spmd_utils,       only: masterproc       
       !       
       implicit none
       !       
       ! local variables       
       !
       integer :: isn,i,j,nmax,k,status
       !
       ! dimension with the largest possible number of points on fieldline       
       !
       real(r8) :: sigH(nhgt_fix),sigP(nhgt_fix)     ! conductivities in [S/m]
       real(r8) :: glon(nhgt_fix),glat(nhgt_fix), &  ! geog. location in [deg]
           alt(nhgt_fix), &                      ! altitude of points [km]
           sumP, &                      ! Pedersen conductance
           sumH                         ! Hall conductance
!!
!! to calculate Jpg
!	real(r8) :: fac     
!	real(r8) :: teiH(nhgt_fix),rhoH(nhgt_fix),neH(nhgt_fix) ! Ti+Te, ion density, electron density on fieldline  
!	real(r8), allocatable :: tei_s1(:,:,:,:),tei_s2(:,:,:,:) ! Ti+Te
!	real(r8), allocatable :: rho_s1(:,:,:,:),rho_s2(:,:,:,:) ! ion density [#/m3*mol]
!	real(r8), allocatable :: ne_s1(:,:,:,:),ne_s2(:,:,:,:)  ! electron density
!	 
       logical,parameter :: debug=.false. 
!         
       if(debug) write(iulog,*) 'in get_conduc'
       !
       ! Get conductance for s2 grid.  
       ! Do not need conductivities from empirical model since they are already calculated
       !
       do isn = 1,2 ! loop over both hemisphere
         do i=mlon0_p,mlon1_p ! loop over task longitudes
           do j=nmlatS2_h,1,-1 ! loop over all latitudes in one hemisphere
              !	   
              ! S2-grid: get the location of the conductivities (this is the s2-grid)
	      !
	      nmax = fline_s2(i,j,isn)%npts ! maximum of points on fieldlinep
!              glat(1:nmax) = fline_s2(i,j,isn)%glat(:)
!              glon(1:nmax) = fline_s2(i,j,isn)%glon(:)
!	      alt(1:nmax)  = fline_s2(i,j,isn)%hgt_pt(:)
!	      alt(1:nmax)  = alt(1:nmax)*m2km     ! convert from m to km
!!
!! calculate conductivities from empirical model
!!
!             if(debug) write(6,*) 'call conduc_empirical' ,i,j,nmax
!              call conduc_empirical(glon,glat,alt,&
!	         nmax,f107,ap,year,doy,ut,sigH,sigP,teiH,rhoH,neH)
!!
!! put conductivity into fieldline array	      
!              fline_s2(i,j,isn)%sigH(:)= sigH(1:nmax)
!              fline_s2(i,j,isn)%sigP(:)= sigP(1:nmax)
!	      do k=1,nmax-1
!	        if(alt(k).lt.140.) then    
!                  fline_s2(i,j,isn)%sigP(k)= 0.5*sigP(k)
! 		endif 
!	      enddo
              !
              ! Calculate s2 conductances from conductivities
	      !	      
	      sumP = 0
	      sumH = 0
	      do k=1,nmax-1
!	       if(alt(k).gt.140.) then
	       sumP = sumP + fline_s2(i,j,isn)%sigP(k)*2*abs(fline_s2(i,j,isn)%Vmp(k+1)-fline_s2(i,j,isn)%Vmp(k))/ &
	         (fline_s2(i,j,isn)%Bmag(k+1)+fline_s2(i,j,isn)%Bmag(k))
	       sumH = sumH + fline_s2(i,j,isn)%sigH(k)*2*abs(fline_s2(i,j,isn)%Vmp(k+1)-fline_s2(i,j,isn)%Vmp(k))/ &
	         (fline_s2(i,j,isn)%Bmag(k+1)+fline_s2(i,j,isn)%Bmag(k))
!		endif 
	      enddo
	      fline_s2(i,j,isn)%zigP = sumP
	      fline_s2(i,j,isn)%zigH = sumH
!!	       
!! for Jpg calculation save Ti+Te, rho_ion, Ne
!	      if(Jpg) then
!		 tei_s2(i,j,:,isn) = teiH(1:nmax)
!		 rho_s2(i,j,:,isn) = rhoH(1:nmax)
!		 ne_s2(i,j,:,isn)  = neH(1:nmax)
!	      endif		       
	       
	   enddo ! end latitude loop
           !
           ! Get conductances for s1 grid.  
           !	   
           do j=nmlat_h,1,-1 ! loop over all latitudes in one hemisphere
              !
	      ! S1-grid: get the location of the conductivities (this is the s1-grid)
	      !
	      nmax = fline_s1(i,j,isn)%npts ! maximum of points on fieldlinep
!              glat(1:nmax) = fline_s1(i,j,isn)%glat(:)
!              glon(1:nmax) = fline_s1(i,j,isn)%glon(:)
!	      alt(1:nmax)  = fline_s1(i,j,isn)%hgt_pt(:)
!	      alt(1:nmax)  = alt(1:nmax)*m2km     ! convert from m to km
!!
!! calculate conductivities from empirical model
!             if(debug) write(iulog,*) 'call conduc_empirical' ,i,j,nmax
!              call conduc_empirical(glon,glat,alt,&
!	         nmax,f107,ap,year,doy,ut,sigH,sigP,teiH,rhoH,neH)		 
!!
!! put conductivity into fieldline array
!              fline_s1(i,j,isn)%sigH(:)= sigH(1:nmax)
!              fline_s1(i,j,isn)%sigP(:)= sigP(1:nmax)
!	      do k=1,nmax-1
!	        if(alt(k).lt.140.) then    
!                  fline_s1(i,j,isn)%sigP(k)= 0.5*sigP(k)
! 		endif 
!	      enddo
	      
! calculate conductances	      
	      sumP = 0.
	      sumH = 0.
	      do k=1,nmax-1
!	       if(alt(k).gt.140.) then
	         sumP = sumP + fline_s1(i,j,isn)%sigP(k)*2*abs(fline_s1(i,j,isn)%Vmp(k+1)-fline_s1(i,j,isn)%Vmp(k))/ &
	          (fline_s1(i,j,isn)%Bmag(k+1)+fline_s1(i,j,isn)%Bmag(k))
	         sumH = sumH + fline_s1(i,j,isn)%sigH(k)*2*abs(fline_s1(i,j,isn)%Vmp(k+1)-fline_s1(i,j,isn)%Vmp(k))/ &
	          (fline_s1(i,j,isn)%Bmag(k+1)+fline_s1(i,j,isn)%Bmag(k))
!	       endif	 
	      enddo
	      fline_s1(i,j,isn)%zigP = sumP
	      fline_s1(i,j,isn)%zigH = sumH
!!	      
!! for Jpg calculation save Ti+Te, rho_ion, Ne
!             if(Jpg) then
!	        tei_s1(i,j,:,isn) = teiH(1:nmax) ! [K]
!	        rho_s1(i,j,:,isn) = rhoH(1:nmax) ! [#/m3*mol]
!	        ne_s1(i,j,:,isn)  = neH(1:nmax)  ! [#/m-3]
!	     endif	      	      
	      
	    enddo ! end loop over all latitudes in one hemisphere
	    
	  enddo  ! end loop over task longitudes
	enddo  ! end loop over both hemisphere
!!	 
!	if(Jpg) then
!	
!	  call calc_Jpg(tei_s1,tei_s2,rho_s1,rho_s2,ne_s1,ne_s2)
!	  
!	  deallocate(tei_s1,STAT=status)
!	  if(status /= 0 ) write(6,*) 'dealloc tei_s1 failed'
!	  deallocate(tei_s2,STAT=status)
!	  if(status /= 0 ) write(6,*) 'dealloc tei_s2 failed'
!	  deallocate(rho_s1,STAT=status)
!	  if(status /= 0 ) write(6,*) 'dealloc rho_s1 failed'
!	  deallocate(rho_s2,STAT=status)
!	  if(status /= 0 ) write(6,*) 'dealloc rho_s2 failed'
!	  deallocate(ne_s1,STAT=status) 
!	  if(status /= 0 ) write(6,*) 'dealloc ne_s1 failed'
!	  deallocate(ne_s2,STAT=status) 
!         if(status /= 0 ) write(6,*) 'dealloc ne_s2 failed'          
!       endif  
!
       end subroutine edyn3D_get_conduct
