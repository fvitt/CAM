#define MSIS_DIAGS
module mo_msis_ubc
!---------------------------------------------------------------
!	... msis upper bndy values
!---------------------------------------------------------------

      use shr_kind_mod,     only: r8 => shr_kind_r8
      use constituents,     only: pcnst

      use cam_abortutils,   only: endrun
      use cam_logfile,      only: iulog
      use cam_history,      only: addfld, horiz_only, outfld

      implicit none

      private
      public  :: msis_ubc_inti, get_msis_ubc, msis_timestep_init

      save

      integer                :: ndx_n=-1, ndx_h=-1, ndx_o=-1, ndx_o2=-1 ! n, h, o, o2 spc indicies
      real(r8), allocatable  :: msis_ubc(:,:,:)                       ! module array for msis ub values (kg/kg)

      logical                :: zonal_average         = .false.       ! use zonal averaged tgcm values

      contains

      subroutine msis_ubc_inti( zonal_avg_in, n_ndx_in,h_ndx_in,o_ndx_in,o2_ndx_in )
!------------------------------------------------------------------
!	... initialize upper boundary values
!------------------------------------------------------------------

        use ppgrid, only : pcols, begchunk, endchunk
        use msis_init, only: msisinit

!------------------------------------------------------------------
!	... dummy args
!------------------------------------------------------------------
      logical, intent(in) :: zonal_avg_in  ! zonal averaging switch
      integer, intent(in) :: n_ndx_in,h_ndx_in,o_ndx_in,o2_ndx_in

!------------------------------------------------------------------
!	... local variables
!------------------------------------------------------------------
      integer  :: astat
      real(r8) :: msis_switches(25) = 1._r8

      call msisinit(parmpath='/terminator-data1/home/fvitt/camdev/waccm_msis_update/src/chemistry/utils/nrlmsis2.1/')

      zonal_average = zonal_avg_in

      if (h_ndx_in>0) then
         ndx_h = h_ndx_in
      endif
      if (n_ndx_in>0) then
         ndx_n = n_ndx_in
      endif
      if (o_ndx_in>0) then
         ndx_o = o_ndx_in
      endif
      if (o2_ndx_in>0) then
         ndx_o2 = o2_ndx_in
      endif

!------------------------------------------------------------------
!	... allocate msis ubc array
!------------------------------------------------------------------
      allocate( msis_ubc(pcols,6,begchunk:endchunk),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'msis_ubc_inti: failed to allocate msis_ubc; error = ',astat
         call endrun('msis_ubc_inti: failed to allocate msis_ubc')
      end if

      if( zonal_average ) then
         msis_switches(7:8)   = 0._r8
         msis_switches(10:14) = 0._r8
      end if

!------------------------------------------------------------------
!	... initialize msis switches
!------------------------------------------------------------------
      call tselec( msis_switches )

      call addfld( 'MSIS_T', horiz_only, 'A', 'K',     'T upper boundary condition from MSIS')
      call addfld( 'MSIS_H', horiz_only, 'A', 'kg/kg', 'H upper boundary condition from MSIS')
      call addfld( 'MSIS_N', horiz_only, 'A', 'kg/kg', 'N upper boundary condition from MSIS')
      call addfld( 'MSIS_O', horiz_only, 'A', 'kg/kg', 'O upper boundary condition from MSIS')
      call addfld( 'MSIS_O2',horiz_only, 'A', 'kg/kg', 'O2 upper boundary condition from MSIS')

      end subroutine msis_ubc_inti

      subroutine msis_timestep_init( ap, f107p_in, f107a_in, state )
!--------------------------------------------------------------------
!	... get the upper boundary values for h, n, o, o2 and temp
!--------------------------------------------------------------------

      use ppgrid,       only : pcols, begchunk, endchunk
      use constituents, only : cnst_mw
      use time_manager, only : get_curr_date, get_calday, get_curr_calday
      use phys_grid,    only : get_ncols_p, get_rlon_all_p, get_rlat_all_p
      use ref_pres,     only : ptop_ref
      use spmd_utils,   only : masterproc
      use physconst,    only : pi
      use cam_control_mod,only : lambm0, eccen, mvelpp, obliqr
      use shr_orb_mod,  only : shr_orb_decl

      use physics_types,only : physics_state

      use msis_constants, only : rp
      use msis_calc, only : msiscalc

!--------------------------------------------------------------------
!	... dummy args
!--------------------------------------------------------------------
      real(r8), intent(in)    ::  ap
      real(r8), intent(in)    ::  f107p_in ! previous day
      real(r8), intent(in)    ::  f107a_in
      type(physics_state), intent(in) :: state(begchunk:endchunk)

!--------------------------------------------------------------------
!	... local variables
!--------------------------------------------------------------------
      real(r8), parameter :: pa2mb       = 1.e-2_r8       ! pascal to mb
      real(r8), parameter :: amu_fac     = 1.65979e-24_r8 ! g/amu
      real(r8), parameter :: r2d         = 180._r8/pi
      integer  ::  i, c, ncol
      integer  ::  yr, mon, day, tod
      integer  ::  yrday
      integer  ::  date
      real(r8) ::  doy
      real(r8) ::  rlons(pcols)
      real(r8) ::  rlats(pcols)
      real(r8) ::  dnom(pcols)
      real(r8) ::  calday, delta, esfact
      real(r8) ::  f107p, f107a

      real(kind=rp) :: ms_day
      real(kind=rp) :: ms_utsec
      real(kind=rp) :: ms_z
      real(kind=rp) :: ms_lat
      real(kind=rp) :: ms_lon
      real(kind=rp) :: ms_sfluxavg,ms_sflux,ms_ap(1:7)
      real(kind=rp) :: ms_tn, ms_dn(1:10)

      !--------------------------------------------------------------------
      !	... get values from msis
      !--------------------------------------------------------------------

      call get_curr_date( yr, mon, day, tod )
      date = 10000*yr + 100*mon + day
      doy = get_calday( date, tod )

      calday = get_curr_calday()

      esfact = 1._r8
      call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, delta, esfact )

      f107p = esfact*f107p_in
      f107a = esfact*f107a_in

      ms_day = real(doy,kind=rp)
      ms_utsec = real(tod,kind=rp)
      ms_ap(:) = 0._rp
      ms_ap(1) = real(ap,kind=rp)
      ms_sfluxavg = real(f107a,kind=rp)
      ms_sflux = real(f107p,kind=rp)

#ifdef MSIS_DIAGS
      if( masterproc ) then
         write(iulog,*) '===================================='
         write(iulog,*) 'msis_timestep_init: diagnostics'
         write(iulog,*) 'yr,mon,day,tod,date,doy,esfact = ', yr, mon, day, tod, date, doy, esfact
         write(iulog,*) '===================================='
      end if
#endif
      chunk_loop : do c = begchunk,endchunk
         ncol = get_ncols_p( c )
         call get_rlat_all_p( c, ncol, rlats )
         call get_rlon_all_p( c, ncol, rlons )
         rlons(:ncol) = r2d * rlons(:ncol)
         rlats(:ncol) = r2d * rlats(:ncol)
         column_loop : do i = 1,ncol

            ms_z = real(state(c)%zi(i,1)*1.e-3_r8,kind=rp) ! km
            ms_lat = real(rlats(i),kind=rp)
            ms_lon = real(rlons(i),kind=rp)

            call msiscalc(ms_day,ms_utsec,ms_z,ms_lat,ms_lon,ms_sfluxavg,ms_sflux,ms_ap,ms_tn,ms_dn)

            msis_ubc(i,1,c) = real(ms_tn, kind=r8) ! temp (K)
            msis_ubc(i,2,c) = real(ms_dn(6),kind=r8)*1.e-6_r8  ! h (molec/cm^3)
            msis_ubc(i,3,c) = real(ms_dn(8),kind=r8)*1.e-6_r8  ! n (molec/cm^3)
            msis_ubc(i,4,c) = real(ms_dn(4),kind=r8)*1.e-6_r8  ! o (molec/cm^3)
            msis_ubc(i,5,c) = real(ms_dn(3),kind=r8)*1.e-6_r8  ! o2 (molec/cm^3)
            msis_ubc(i,6,c) = real(ms_dn(1),kind=r8)*1.e-3_r8  ! total atm dens (g/cm^3)

#ifdef MSIS_DIAGS
            if( masterproc ) then
               write(iulog,*) '===================================='
               write(iulog,*) 'msis_timestep_init: diagnostics for col,chnk = ',i,c
               write(iulog,*) 'day, utsec, alt = ',ms_day,ms_utsec,ms_z
               write(iulog,*) 'msis_temp = ', ms_tn
               write(iulog,*) 'new msis day,utset,z,lat,lon,Tn : ',ms_day,ms_utsec,ms_z,ms_lat,ms_lon,ms_tn
               write(iulog,*) 'msis h,n,o,o2,m = ',msis_ubc(i,2:6,c)
               write(iulog,'(a,5e12.3)') 'msis h,n,o,o2,m = ',msis_ubc(i,2:6,c)
               write(iulog,*) '===================================='
            endif
#endif
         end do column_loop

         !--------------------------------------------------------------------
         !	... transform from molecular density to mass mixing ratio
         !--------------------------------------------------------------------
         dnom(:ncol) = amu_fac/msis_ubc(:ncol,6,c)
         if( ndx_h > 0 ) then
            msis_ubc(:ncol,2,c) = cnst_mw(ndx_h)*msis_ubc(:ncol,2,c)*dnom(:ncol)
         end if
         if( ndx_n > 0 ) then
            msis_ubc(:ncol,3,c) = cnst_mw(ndx_n)*msis_ubc(:ncol,3,c)*dnom(:ncol)
         end if
         if( ndx_o > 0 ) then
            msis_ubc(:ncol,4,c) = cnst_mw(ndx_o)*msis_ubc(:ncol,4,c)*dnom(:ncol)
         end if
         if( ndx_o2 > 0 ) then
            msis_ubc(:ncol,5,c) = cnst_mw(ndx_o2)*msis_ubc(:ncol,5,c)*dnom(:ncol)
         end if
      end do chunk_loop

      end subroutine msis_timestep_init

      subroutine get_msis_ubc( lchunk, ncol, temp, mmr )
!--------------------------------------------------------------------
!	... get the upper boundary values for h, n, o, o2 and temp
!--------------------------------------------------------------------

      use ppgrid,       only : pcols

!--------------------------------------------------------------------
!	... dummy args
!--------------------------------------------------------------------
      integer, intent(in)     :: lchunk            ! chunk id
      integer, intent(in)     :: ncol              ! columns in chunk
      real(r8), intent(out)   :: temp(pcols)       ! msis temperature at top interface (K)
      real(r8), intent(inout) :: mmr(pcols,pcnst)  ! msis concentrations at top interface (kg/kg)

!--------------------------------------------------------------------
!	... set model ubc values from msis
!--------------------------------------------------------------------
      temp(:ncol) = msis_ubc(:ncol,1,lchunk)

      call outfld( 'MSIS_T', msis_ubc(:ncol,1,lchunk), ncol, lchunk)
      call outfld( 'MSIS_H', msis_ubc(:ncol,2,lchunk), ncol, lchunk)
      call outfld( 'MSIS_N', msis_ubc(:ncol,3,lchunk), ncol, lchunk)
      call outfld( 'MSIS_O', msis_ubc(:ncol,4,lchunk), ncol, lchunk)
      call outfld( 'MSIS_O2',msis_ubc(:ncol,5,lchunk), ncol, lchunk)

      if( ndx_h > 0 ) then
         mmr(:ncol,ndx_h) = msis_ubc(:ncol,2,lchunk)
      end if
      if( ndx_n > 0 ) then
         mmr(:ncol,ndx_n) = msis_ubc(:ncol,3,lchunk)
      end if
      if( ndx_o > 0 ) then
         mmr(:ncol,ndx_o) = msis_ubc(:ncol,4,lchunk)
      end if
      if( ndx_o2 > 0 ) then
         mmr(:ncol,ndx_o2) = msis_ubc(:ncol,5,lchunk)
      end if

      end subroutine get_msis_ubc

      end module mo_msis_ubc
