     module edyn3D_calc_coef_fac_const_rhs

     use edyn3D_params, only: nmlon,nmlat_h,nhgt_fix,nlonlat,nmlatS2_h,nmlat_T1
     use edyn3d_mpi,       only: mlon0_p,mlon1_p
     use shr_kind_mod,  only: r8 => shr_kind_r8            ! 8-byte reals
     use cam_logfile,    only: iulog
     use spmd_utils,     only: masterproc
     use cam_abortutils, only: endrun

     implicit none
     private

     public :: edyn3D_calc_coef,edyn3D_calc_FAC,edyn3D_add_coef_ns,edyn3D_gather_coef_ns, &
               edyn3D_const_rhs,edyn3D_scatter_poten

     real(r8), allocatable :: coef(:,:,:,:,:)       ! for each P-point each hemisphere and height
     real(r8), allocatable :: coef_ns(:,:,:)        ! lhs+rhs: for each P-point
     real(r8), allocatable :: coef_ns_glb(:,:,:)    ! lhs+rhs:  globally
     real(r8), allocatable :: coef_ns2(:,:,:,:)     ! lhs+rhs: for each P-point
     real(r8), allocatable :: rhs_ns(:)             ! forcing: for each P-point
     real(r8), allocatable :: lhs_ns(:,:)           ! lhs: for each P-point
     real(r8), allocatable :: poten_glb(:,:,:)      ! Potential: globally

     contains

!-----------------------------------------------------------------------------
     subroutine edyn3D_calc_coef
!
     use edyn3D_fieldline, only: fieldline_s1,fline_s1,fieldline_s2,fline_s2,fline_p

     implicit none
!
! coef ordering from tiegcm
! ^   equatorward
! coef(4) (i-1,j+1)      coef(3) (i,j+1)    coef(2) (i+1,j+1)
! coef(5) (i-1,j)        coef(9) (i,j)      coef(1) (i+1,j)
! coef(6) (i-1,j-1)      coef(7) (i,j-1)    coef(8) (i+1,j-1)
! v   poleward
!
! relationship between P,S1, and S2 point for the same index (i,j)
!  P(i,j) then is really S1(i+0.5,j) and S2(i,j+0.5) with j increasing equatorward
!  coefficient is calculated at P points
!

      integer :: isn,i,j,k,im,nmax,status,ic
      real(r8) :: N2p_p,N2h_p

      !
      ! Allocate arrays for solver left hand side coefficients
      !
      allocate(coef(mlon0_p:mlon1_p,nmlat_h,nhgt_fix,9,2),STAT=status) ! for each P-point each hemisphere and height
      if (status /= 0 ) then
        write(iulog,*) 'alloc coef failed'
	call endrun('edyn3d_calc_coef')
      endif
      allocate(coef_ns(mlon0_p:mlon1_p,nmlat_h,10),STAT=status)        ! lhs+rhs: for each P-point
      if(status /= 0 ) then
        write(iulog,*) 'alloc coef_ns failed'
	call endrun('edyn3d_calc_coef')
      endif
      allocate(coef_ns_glb(nmlon,nmlat_h,10),STAT=status)              ! lhs+rhs: globally
      if(status /= 0 ) then
        write(iulog,*) 'alloc coef_ns_glb failed'
	call endrun('edyn3d_calc_coef')
      endif
      allocate(rhs_ns(nlonlat),STAT=status)                            ! forcing: for each P-point globally
      if(status /= 0 ) then
        write(iulog,*) 'alloc rhs_ns failed'
	call endrun('edyn3d_calc_coef')
      endif
      allocate(lhs_ns(nlonlat,nlonlat),STAT=status)                    ! lhs: globally
      if(status /= 0 ) then
        write(iulog,*) 'alloc lhs_ns failed'
	call endrun('edyn3d_calc_coef')
      endif
      allocate(coef_ns2(mlon0_p:mlon1_p,nmlat_h,2,10),STAT=status)     ! lhs+rhs: for each P-point both hemispheres
      if(status /= 0 ) then
        write(iulog,*) 'alloc coef_ns2 failed'
	call endrun('edyn3d_calc_coef')
      endif
      allocate(poten_glb(nmlon,nmlat_h,2),STAT=status)                 ! Potential: globally (both hemispheres)
      if(status /= 0 ) then
        write(iulog,*) 'alloc poten_glb failed'
	call endrun('edyn3d_calc_coef')
      endif

       do isn = 1,2 ! loop over both hemisphere
         do i=mlon0_p,mlon1_p ! loop over task longitudes
!         do i=1,nmlon ! loop over all longitudes

!           im = i-1

           if(i == mlon0_p) then
             im = mlon1_p
           else
             im = i-1
           endif

!           if(i == 1) then
!             im = nmlon
!           else
!             im = i-1
!           endif

           do j=2,nmlat_h ! loop over all latitudes in one hemisphere NOT THE POLE
              nmax = fline_p(i,j,isn)%npts ! maximum of points on fieldline

              do k=1,nmax  ! longitudinal wrap around
                if((nmlat_h-j+1) == k) then ! top volume at equator
                  N2h_p = 0._r8
                  N2p_p = 0._r8
                else
                  N2h_p = fline_s2(i,j,isn)%N2h(k)  ! i,j+0.5
                  N2p_p = fline_s2(i,j,isn)%N2p(k)  ! i,j+0.5
                endif
                coef(i,j,k,2,isn) = -fline_s1(i,j,isn)%N1h(k) +  N2h_p
                coef(i,j,k,3,isn) = fline_s1(im,j,isn)%N1h(k) - &
                 fline_s1(i,j,isn)%N1h(k)+ N2p_p
                coef(i,j,k,4,isn) = fline_s1(im,j,isn)%N1h(k) - N2h_p
!
                coef(i,j,k,5,isn) = fline_s1(im,j,isn)%N1p(k) + &
                   fline_s2(i,j-1,isn)%N2h(k)-N2h_p
                coef(i,j,k,6,isn) = -fline_s1(im,j,isn)%N1h(k) + &
                   fline_s2(i,j-1,isn)%N2h(k)
                coef(i,j,k,7,isn) = -fline_s1(im,j,isn)%N1h(k) + &
                   fline_s1(i,j,isn)%N1h(k)+fline_s2(i,j-1,isn)%N2p(k)
                coef(i,j,k,8,isn) = fline_s1(i,j,isn)%N1h(k) - &
                   fline_s2(i,j-1,isn)%N2h(k)
                coef(i,j,k,1,isn) = fline_s1(i,j,isn)%N1p(k) - &
                   fline_s2(i,j-1,isn)%N2h(k)+N2h_p
                coef(i,j,k,9,isn) = -fline_s1(im,j,isn)%N1p(k) - &
                   fline_s1(i,j,isn)%N1p(k)-fline_s2(i,j-1,isn)%N2p(k) -&
                   N2p_p
              end do  ! end height loop
           end do  ! end lat/fieldline loop
!
         end do  ! end longitude loop
       end do ! end hemisphere loop
!
! add the coefficients in height to get coefficients for each hemisphere
! needed for calculating high latitude FAC
! initialize
     coef_ns2 = 0.
 !
     do i=mlon0_p,mlon1_p ! loop over task longitudes
!     do i=1,nmlon   ! loop over all longitudes
       do isn = 1,2   ! hemisphere loop
        j=1               ! polar values are set
        coef_ns2(i,j,isn,9)   = 0.5_r8  ! later added together hemispheres to get one
        coef_ns2(i,j,isn,1:8) = 0._r8
        coef_ns2(i,j,isn,10)  = 0._r8  ! set potential at pole
        do j=2,nmlat_h ! loop over all latitudes in one hemisphere NOT THE POLE
          nmax = fline_p(i,j,1)%npts ! maximum of points on fieldline
          do k=1,nmax  ! height loop
            do ic = 1,9 ! 9-point stencil
              coef_ns2(i,j,isn,ic) = coef_ns2(i,j,isn,ic)+ coef(i,j,k,ic,isn)
            end do
            coef_ns2(i,j,isn,10) = coef_ns2(i,j,isn,10)+ fline_p(i,j,isn)%S(k)
          end do  ! end height loop
!
        end do  ! end lat/fieldline loop
       end do  ! end hemisphere loop
!
     end do  ! end longitude loop
!

!write(iulog,*) 'edyn3D_calc_coef: MIN/MAX coef_ns2 ', &
!                 MINVAL(coef_ns2(mlon0_p:mlon1_p,:,:,1)), MAXVAL(coef_ns2(mlon0_p:mlon1_p,:,:,1))

     end subroutine edyn3D_calc_coef
!-----------------------------------------------------------------------------
     subroutine edyn3D_calc_FAC
!
     use edyn3D_fieldline,only: fline_p,fieldline_s1,fline_s1,poten_hl
!
     implicit none
!
     real(r8), dimension(nmlon,nmlat_h,2) :: fac_hl
     real(r8) :: sum,sumP,corr,sumn,sums
!
     integer :: isn,i,j,jj,icof,im,ip
!
     fac_hl = 0._r8
     sum    = 0._r8
     sumn   = 0._r8
     sums   = 0._r8
     sumP   = 0._r8
!
!  Set high latitude potential to 0.01 since no input
!
     poten_hl(:,:) = 0.01_r8
!
!     coef_ns2(:,:,:,1:8) = 0._r8

!write(iulog,*) 'edyn3D_calc_FAC: MIN/MAX coef_ns2 ', &
!                 MINVAL(coef_ns2(mlon0_p:mlon1_p,:,:,1)), MAXVAL(coef_ns2(mlon0_p:mlon1_p,:,:,1))

     do isn = 1,2 ! loop over both hemisphere
       do i=mlon0_p,mlon1_p ! loop over task longitudes
!       do i=1,nmlon ! loop over all longitudes

	 if(i == mlon0_p) then ! wrap around in longitude
	   im = mlon1_p
	 else
	   im = i-1
	 endif
	 if(i == mlon1_p) then
	   ip = mlon0_p
	 else
	   ip = i+1
	 endif

!         if(i == 1) then ! wrap around in longitude
!           im = nmlon
!         else
!           im = i-1
!         endif
!         if(i == nmlon) then
!           ip = 1
!         else
!           ip = i+1
!         endif

         do j=2,nmlat_h ! loop over all latitudes in one hemisphere not the pole (potential set later)
           if(isn.eq.1) then  ! jj is latitude index from pole to pole
            jj = j
             fac_hl(i,j,isn) = poten_hl(ip,jj)*coef_ns2(i,j,isn,1)+ &
             poten_hl(ip,jj+1)*coef_ns2(i,j,isn,2)+ &
             poten_hl(i ,jj+1)*coef_ns2(i,j,isn,3)+ &
             poten_hl(im,jj+1)*coef_ns2(i,j,isn,4)+ &
             poten_hl(im,jj  )*coef_ns2(i,j,isn,5)+ &
             poten_hl(im,jj-1)*coef_ns2(i,j,isn,6)+ &
             poten_hl(i ,jj-1)*coef_ns2(i,j,isn,7)+ &
             poten_hl(ip,jj-1)*coef_ns2(i,j,isn,8)+ &
             poten_hl(i  ,jj )*coef_ns2(i,j,isn,9)
           else
            jj = nmlat_T1 - j + 1
            fac_hl(i,j,isn) = poten_hl(ip,jj)*coef_ns2(i,j,isn,1)+ &
             poten_hl(ip,jj-1)*coef_ns2(i,j,isn,2)+ &
             poten_hl(i ,jj-1)*coef_ns2(i,j,isn,3)+ &
             poten_hl(im,jj-1)*coef_ns2(i,j,isn,4)+ &
             poten_hl(im,jj   )*coef_ns2(i,j,isn,5)+ &
             poten_hl(im,jj+1)*coef_ns2(i,j,isn,6)+ &
             poten_hl(i ,jj+1)*coef_ns2(i,j,isn,7)+ &
             poten_hl(ip,jj+1)*coef_ns2(i,j,isn,8)+ &
             poten_hl(i ,jj   )*coef_ns2(i,j,isn,9)
           endif
            if(fline_s1(i,j,isn)%zigP.gt.1.5) then
              sum  = sum  + fac_hl(i,j,isn)
              sumP = sumP + fline_s1(i,j,isn)%zigP*abs(fac_hl(i,j,isn))
              !
              if(isn.eq.1) sums= sums+fac_hl(i,j,isn)
              if(isn.eq.2) sumn= sumn+fac_hl(i,j,isn)
            else
              fac_hl(i,j,isn) = 0._r8
            end if
         end do  ! end lat/fieldline loop
         !
       end do  ! end longitude loop
     end do ! end hemisphere loop
!
! next two lines commented out to make sure high latitude forcing is included
!     fac_hl= 0._r8
!     corr =  0._r8
!
!     write(iulog,*) 'edyn3D_calc_FAC: sum, sumP', sum,sumP
!     write(iulog,*) 'edyn3D_calc_FAC: Before sum/sumP fac_hl max and min', maxval(fac_hl),minval(fac_hl)
!     if (isclose(sumP,0.0_rp)) then
!       corr = 0
!     else
!       corr = sum/sumP
!     endif
     corr = sum/sumP
!     write(iulog,*) 'edyn3D_calc_FAC: fac_hl max and min', maxval(fac_hl),minval(fac_hl)
!     write(iulog,*) 'edyn3D_calc_FAC: correction fac*area/[zigP*|fac|*area]', corr
!     write(iulog,*) 'edyn3D_calc_FAC: sum fac*area', sum,'sumn fac*area ', sumn, 'sums fac*area ',sums
!
     sum = 0.
     do isn = 1,2 ! loop over both hemisphere
       do i=mlon0_p,mlon1_p ! loop over task longitudes
!       do i=1,nmlon ! loop over all longitudes
         do j=2,nmlat_h ! loop over all latitudes in one hemisphere not the pole (potential set later)

           fac_hl(i,j,isn) = fac_hl(i,j,isn)-fline_s1(i,j,isn)%zigP*abs(fac_hl(i,j,isn))*corr
           sum = sum   + fac_hl(i,j,isn)
           ! put in coef-array
           coef_ns2(i,j,isn,10) =   coef_ns2(i,j,isn,10)+fac_hl(i,j,isn)
           !
!           write(iulog,'(4(x,e15.8))') fline_p(i,j,isn)%mlon_qd(1), &
!                fline_p(i,j,isn)%mlat_qd(1),fac_hl(i,j,isn), fline_s1(i,j,isn)%zigP
           if(fline_p(i,j,isn)%M3(1).ne.0) then
!             fline_p(i,j,isn)%fac_hl = coef_ns2(i,j,isn,10)/ fline_p(i,j,isn)%M3(1)
                fline_p(i,j,isn)%fac_hl = fac_hl(i,j,isn)/ fline_p(i,j,isn)%M3(1)
           endif
         end do  ! end lat/fieldline loop
         !
       end do  ! end longitude loop
     end do ! end hemisphere loop
!
!     write(iulog,*) 'edyn3D_calc_FAC: after corr. sum fac_hl', sum
!
     end subroutine edyn3D_calc_FAC
!-----------------------------------------------------------------------------
     subroutine edyn3D_add_coef_ns
! adds the coefficients from both hemispheres together
!  output:  coef_ns
     use edyn3D_fieldline,only: fline_p

     implicit none
!
     integer :: i,j,k,ic,nmax,status
!
! add hemispheres together and set equatorial boundary condition
     do i=mlon0_p,mlon1_p ! loop over task longitudes
!     do i=1,nmlon ! loop over all longitudes
       do j=1,nmlat_h ! loop over all latitudes in one hemisphere
          do ic = 1,9 ! 9-point stencil
             coef_ns(i,j,ic) = coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic)
          end do
          coef_ns(i,j,10) = coef_ns2(i,j,1,10)+ coef_ns2(i,j,2,10)
!
       end do  ! end lat/fieldline loop
       j=nmlat_h                    ! set equatorial values (page 14 Art's notes)
       coef_ns(i,j,2) = 0._r8
       coef_ns(i,j,3) = 0._r8
       coef_ns(i,j,4) = 0._r8
!
     end do  ! end longitude loop
!
     deallocate(coef,STAT=status)
     if(status /= 0) then
       write(iulog,*) 'deallocation of coef not successful'
       call endrun('edyn3D_add_coef_ns: deallocation of coef not successful')
     endif
     deallocate(coef_ns2,STAT=status)
     if(status /= 0) then
       write(iulog,*) 'deallocation of coef_ns2 not successful'
       call endrun('edyn3D_add_coef_ns: deallocation of coef_ns2 not successful')
     endif

     end subroutine edyn3D_add_coef_ns

!--------------------------------------------------------------------------------
     subroutine edyn3D_gather_coef_ns
!
!    Gather needed global arrays to root task, so it can finish non-parallel
!    part of dynamo (beginning after sub rhspde) as in original code

       use edyn3D_mpi,   only: mytid,ntask,mp_gather_edyn3D
!
!    Local:
!    10 fields to gather: coef_ns(:,:,1-10)
!
       integer, parameter :: nf = 10
       real(r8) :: fmsub(mlon0_p:mlon1_p,nmlat_h,nf)
       real(r8) :: fmglb(nmlon,nmlat_h,nf)
       integer :: i,j,jj,status
!
!    This call excludes halo points in coef_ns
!
       fmsub(:,:,:) = coef_ns(:,:,:)

       call mp_gather_edyn3D(fmsub,mlon0_p,mlon1_p,fmglb,nmlon,nmlat_h,nf)
!
!    Now root task can take over, and work with global arrays:
!
       if (mytid==0) then

	 coef_ns_glb(:,:,:) = fmglb(:,:,:)

       endif ! mytid==0

!write(iulog,*) 'edyn3D_gather_coef_ns: MIN/MAX coef_ns, coef_ns_glb ', &
!                 MINVAL(coef_ns), MAXVAL(coef_ns), MINVAL(coef_ns_glb(mlon0_p:mlon1_p,:,:)), MAXVAL(coef_ns_glb(mlon0_p:mlon1_p,:,:))

!
       deallocate(coef_ns,STAT=status)
       if(status /= 0) then
         write(iulog,*) 'deallocation of coef_ns not successful'
         call endrun('edyn3D_gather_coef_ns: deallocation of coef_ns not successful')
       endif
!
     end subroutine edyn3D_gather_coef_ns

!-----------------------------------------------------------------------------
     subroutine edyn3D_const_rhs
     !
! construct matrix lhs &rhs  for solving with matlab

     use edyn3D_fieldline,only: fline_p
!
!     include 'mkl_lapack.fi'

!     implicit none
!
     integer :: i,j,im,ip,ijm,it,ij,jt,status,isn
!  for solver
     integer:: info,ifail
     integer,parameter :: nrhmax = 1
     integer ::  ipiv(nlonlat)
     character, parameter :: trans='N'
! condition number
     real(r8) :: colsum,anorm ,rcond,work(4*nlonlat)
!     real(r8) :: dlange
     integer :: info_c,iwork(nlonlat)
     character, parameter :: norm='1'
! am 1/2015 for testing
     real(r8) :: pot_ns(nlonlat)
     real(r8) :: z_ns(nlonlat)
     !
!     real(r8) :: rhs_ns(nlonlat)
!     real(r8) :: lhs_ns(nlonlat,nlonlat)
     !
     lhs_ns(:,:) = 0._r8
     rhs_ns(:) = 0._r8
     ij = 0
     isn = 1


! am 1/2015 test by putting analytical potential into stencil and compare RHS
!
     do i=1,nmlon       ! loop over all longitudes
!
       if(i == 1) then ! wrap around in longitude
         im = nmlon
       else
         im = i-1
       endif
       if(i == nmlon) then
         ip = 1
       else
         ip = i+1
       endif
       ! pole values
       j=1  ! rhs = 0. c9=1, and c(1:8) = 0 for the pole
       ij = (i-1)*nmlat_h+j
       rhs_ns(ij) = coef_ns_glb(i,j,10)
       it = (ip-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns_glb(i,j,1)
       it = (ip-1)*nmlat_h+j+1
       lhs_ns(ij,it)  = coef_ns_glb(i,j,2)
       it = (i-1)*nmlat_h+j+1
       lhs_ns(ij,it)  = coef_ns_glb(i,j,3)
       it = (im-1)*nmlat_h+j+1
       lhs_ns(ij,it)  = coef_ns_glb(i,j,4)
       it = (im-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns_glb(i,j,5)
       it = (i-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns_glb(i,j,9)
       !
!       pot_ns(ij) = fline_p(i,j,isn)%pot_test ! am 1/2015 test
       !
       do j=2,nmlat_h-1   ! pole to equator
           ij = (i-1)*nmlat_h+j
           rhs_ns(ij) = coef_ns_glb(i,j,10)
           it = (ip-1)*nmlat_h+j
           lhs_ns(ij,it)  = coef_ns_glb(i,j,1)
           it = (ip-1)*nmlat_h+j+1
           lhs_ns(ij,it)  = coef_ns_glb(i,j,2)
           it = (i-1)*nmlat_h+j+1
           lhs_ns(ij,it)  = coef_ns_glb(i,j,3)
           it = (im-1)*nmlat_h+j+1
           lhs_ns(ij,it)  = coef_ns_glb(i,j,4)
           it = (im-1)*nmlat_h+j
           lhs_ns(ij,it)  = coef_ns_glb(i,j,5)
           it = (im-1)*nmlat_h+j-1
           lhs_ns(ij,it)  = coef_ns_glb(i,j,6)
           it = (i-1)*nmlat_h+j-1
           lhs_ns(ij,it)  = coef_ns_glb(i,j,7)
           it = (ip-1)*nmlat_h+j-1
           lhs_ns(ij,it)  = coef_ns_glb(i,j,8)
           it = (i-1)*nmlat_h+j
           lhs_ns(ij,it)  = coef_ns_glb(i,j,9)
       !
!           pot_ns(ij) = fline_p(i,j,isn)%pot_test  ! am 1/2015 test
       !
       end do  ! end lat/fieldline loop
       ! equator values
       j=nmlat_h
       ij = (i-1)*nmlat_h+j
       rhs_ns(ij) = coef_ns_glb(i,j,10)
       it = (ip-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns_glb(i,j,1)
       it = (im-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns_glb(i,j,5)
       it = (im-1)*nmlat_h+j-1
       lhs_ns(ij,it)  = coef_ns_glb(i,j,6)
       it = (i-1)*nmlat_h+j-1
       lhs_ns(ij,it)  = coef_ns_glb(i,j,7)
       it = (ip-1)*nmlat_h+j-1
       lhs_ns(ij,it)  = coef_ns_glb(i,j,8)
       it = (i-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns_glb(i,j,9)
       !
!       pot_ns(ij) = fline_p(i,j,isn)%pot_test  ! am 1/2015 test
       !
     end do  ! end longitude loop
!
! write(iulog,*) 'edyn3D_const_rhs: MIN/MAX lhs_ns(:,:)', &
!                 MINVAL(lhs_ns(:,:)), MAXVAL(lhs_ns(:,:))
!
!     deallocate(coef_ns_glb,STAT=status)
!     if(status /= 0) then
!       write(iulog,*) 'deallocation of coef_ns_glb not successful'
!       call endrun('edyn3D_scatter_poten')
!     endif
!
! put pot_ns into the lhs lhs_ns -> lhs_ns*pot_ns= rhs_test
!
!       z_ns= 0._r8
!        z_ns = matmul( lhs_ns,pot_ns)

! solve the matrix LHS X = RHS
!
! calculate norm
      anorm = 0._r8
      do i=1,nlonlat
        colsum = 0._r8
        do j=1,nlonlat
      ! compute norm-1 of A-> max_j SUM_i abs(a_ij)
           colsum = colsum + abs(lhs_ns(j,i))
        end do
        anorm = max(anorm,colsum)
      end do
!      write(iulog,*) 'anorm ',anorm
!
!       Factorize A
!
!  write(iulog,*) 'edyn3D_const_rhs: MIN/MAX ABS(lhs_ns(:,:))', &
!                 MINVAL(ABS(lhs_ns(:,:))), MAXVAL(ABS(lhs_ns(:,:)))

!      write (iulog,*) 'factorize A'
      call DGETRF(nlonlat,nlonlat,lhs_ns,nlonlat,ipiv,info)
!      write(iulog,*) 'info ' , info
!
      if (info.EQ.0) then
!
! calculate condition number
        CALL DGECON( norm,nlonlat,lhs_ns,nlonlat,anorm,rcond,work,iwork,info_c)
!        write(iulog,*) 'info_c ' , info_c
!        write(iulog,*) 'condit=', 1/rcond
!
!        Compute solution
!
!         write (iulog,*) 'solve'
         call DGETRS(trans,nlonlat,nrhmax,lhs_ns,nlonlat,ipiv,rhs_ns,nlonlat,info)
!         write (iulog,*) 'after solve',info
!
!        Print solution
!
!         ifail = 0
!         call X04CAF('General',' ',nlonlat,nrhmax,rhs_ns,nlonlat,'Solution(s)',ifail)
         !
        it=0
         do i=1,nmlon
           do j=1,nmlat_h
            it = it +1
!            write(iulog,'(2(i4,x),1(x,e15.8))') i,j,rhs_ns(it)
            poten_glb(i,j,1) = rhs_ns(it)
            poten_glb(i,j,2) = rhs_ns(it)
!            fline_p(i,j,1)%pot = rhs_ns(it)
!            fline_p(i,j,2)%pot = rhs_ns(it)
           enddo
         enddo
      else
         write (iulog,*) 'edyn3D_const_rhs: The factor U is singular'
      end if
!
! write(iulog,*) 'edyn3D_const_rhs: MIN/MAX poten_glb(:,:,:)', &
!                 MINVAL(poten_glb(:,:,:)), MAXVAL(poten_glb(:,:,:))

!     deallocate(lhs_ns,STAT=status)
!     if(status /= 0) then
!	write(iulog,*) 'deallocation of lhs_ns not successful'
!	call endrun('edyn3D_scatter_poten')
!     endif
!
!     deallocate(rhs_ns,STAT=status)
!     if(status /= 0) then
!	write(iulog,*) 'deallocation of rhs_ns not successful'
!	call endrun('edyn3D_scatter_poten')
!     endif
!
     end subroutine edyn3D_const_rhs
!
!--------------------------------------------------------------------------------
     subroutine edyn3D_scatter_poten
!
!    Scatter root task to global arrays, after finishing non-parallel
!    part of dynamo (beginning after sub edyn3D_const_rhs) as in original code

       use edyn3D_mpi,   only: mytid,ntask,mp_scatter_edyn3D
       use edyn3D_fieldline,only: fline_p
!
!    Local:
!    Need to scatter 2 hemispheres: pot(:,:,2)
!
       real(r8) :: potensub(mlon0_p:mlon1_p,nmlat_h,2)
       real(r8) :: potenglb1(nmlon,nmlat_h)
       real(r8) :: potenglb2(nmlon,nmlat_h)
       real(r8) :: potensub1(nmlon,nmlat_h)
       real(r8) :: potensub2(nmlon,nmlat_h)

       integer :: i,j,jj,status
       !
       ! Scatter twice, once for each hemisphere
       !
       potenglb1(:,:) = poten_glb(:,:,1)
       call mp_scatter_edyn3D(potenglb1,mlon0_p,mlon1_p,potensub1,nmlon,nmlat_h)

       potenglb2(:,:) = poten_glb(:,:,2)
       call mp_scatter_edyn3D(potenglb2,mlon0_p,mlon1_p,potensub2,nmlon,nmlat_h)
       !
       ! Add potential to field line structure
       !
       do i=mlon0_p,mlon1_p ! loop over task longitudes
         do j=1,nmlat_h ! loop over all latitudes in one hemisphere

           fline_p(i,j,1)%pot = potensub1(i,j)
           fline_p(i,j,2)%pot = potensub2(i,j)

         end do  ! end lat/fieldline loop
       end do  ! end longitude loop
!
!write(iulog,*) 'edyn3D_scatter_poten: MIN/MAX potensub1, potensub2 ', &
!                 MINVAL(potensub1(mlon0_p:mlon1_p,:)), MAXVAL(potensub1(mlon0_p:mlon1_p,:)), &
!		 MINVAL(potensub2(mlon0_p:mlon1_p,:)), MAXVAL(potensub2(mlon0_p:mlon1_p,:))

!
       deallocate(coef_ns_glb,STAT=status)
       if(status /= 0) then
         write(iulog,*) 'deallocation of coef_ns_glb not successful'
         call endrun('edyn3D_scatter_poten')
       endif

       deallocate(poten_glb,STAT=status)
       if(status /= 0) then
          write(iulog,*) 'deallocation of poten_glb not successful'
	  call endrun('edyn3D_scatter_poten')
        endif

       deallocate(lhs_ns,STAT=status)
       if(status /= 0) then
         write(iulog,*) 'deallocation of lhs_ns not successful'
         call endrun('edyn3D_scatter_poten')
       endif
!
       deallocate(rhs_ns,STAT=status)
       if(status /= 0) then
         write(iulog,*) 'deallocation of rhs_ns not successful'
         call endrun('edyn3D_scatter_poten')
       endif
!
     end subroutine edyn3D_scatter_poten

!-----------------------------------------------------------------------
     end module edyn3D_calc_coef_fac_const_rhs
!--------------------------------------------------------------------------------------------
