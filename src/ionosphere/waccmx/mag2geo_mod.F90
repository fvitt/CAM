module mag2geo_mod
  use shr_kind_mod   ,only: r8 => shr_kind_r8
  use cam_logfile    ,only: iulog
  use cam_abortutils, only: endrun
#ifdef WACCMX_EDYN_ESMF
  use esmf           ,only: ESMF_FIELD
  use edyn_mpi       ,only: mlon0,mlon1,mlat0,mlat1, lon0,lon1,lat0,lat1
  use edyn_params    ,only: finit
#endif

  implicit none
  private

#ifdef WACCMX_EDYN_ESMF

  public :: boxcar_ave, mag2geo, rpt_ncerr
  public :: check_ncerr
  public :: check_alloc
#endif

contains
#ifdef WACCMX_EDYN_ESMF
  !-----------------------------------------------------------------------
  subroutine boxcar_ave(x,y,lon,lat,mtime,itime,ibox)
    !
    ! perform boxcar average
    !
    ! Args:
    integer,  intent(in)  :: lon
    integer,  intent(in)  :: lat
    integer,  intent(in)  :: mtime
    integer,  intent(in)  :: itime
    integer,  intent(in)  :: ibox
    real(r8), intent(in)  :: x(lon,lat,mtime)
    real(r8), intent(out) :: y(lon,lat)
    
    ! Local:
    integer :: i, iset, iset1
    !
    iset = itime - ibox/2
    if (iset < 1) iset = 1
    iset1 = iset + ibox
    if (iset1 > mtime) then
       iset1 = mtime
       iset = iset1 - ibox
    end if
    !     write(iulog,"('boxcar_ave: mtime,itime,ibox',3i5)")
    !    |  mtime,itime,ibox
    !
    y(:,:) = 0._r8
    do i=iset,iset1
       y(:,:) = y(:,:) + x(:,:,i)
    end do
    if (ibox > 0) y(:,:) = y(:,:)/ibox
    !
  end subroutine boxcar_ave
  !-----------------------------------------------------------------------
  subroutine mag2geo(am,ag,im,jm,dim,djm,lg,lm,nlong,nlatg)
    !
    ! Args:
    integer,  intent(in)  :: lg
    integer,  intent(in)  :: lm
    real(r8), intent(in)  :: am(lm,*)
    real(r8), intent(out) :: ag(lg,*)
    integer,  intent(in)  :: im(lg,*)
    integer,  intent(in)  :: jm(lg,*)
    real(r8), intent(in)  :: dim(lg,*)
    real(r8), intent(in)  :: djm(lg,*)
    integer,  intent(in)  :: nlong
    integer,  intent(in)  :: nlatg
    !
    ! Local:
    integer :: ig,jg
    !
    do jg=1,nlatg
       do ig=1,nlong
          ag(ig,jg) =  &
               am(im(ig,jg)  ,jm(ig,jg))  *(1._r8-dim(ig,jg))*(1._r8-djm(ig,jg))+  &
               am(im(ig,jg)+1,jm(ig,jg))  *    dim(ig,jg) *(1._r8-djm(ig,jg))+  &
               am(im(ig,jg)  ,jm(ig,jg)+1)*(1._r8-dim(ig,jg))*djm(ig,jg)+  &
               am(im(ig,jg)+1,jm(ig,jg)+1)*    dim(ig,jg) *djm(ig,jg)
       end do ! ig=1,nlong
    end do ! jg=1,nlatg
  end subroutine mag2geo
!!$  !-----------------------------------------------------------------------
!!$  subroutine mag2geo_2d(fmag,fgeo,ESMF_mag,ESMF_geo,fname)
!!$    !
!!$    ! Convert field on geomagnetic grid fmag to geographic grid in fgeo.
!!$    !
!!$    use edyn_esmf,only: edyn_esmf_set2d_mag,edyn_esmf_regrid,  &
!!$                        edyn_esmf_get_2dfield
!!$    !
!!$    ! Args:
!!$    real(r8),         intent(in)    :: fmag(mlon0:mlon1,mlat0:mlat1)
!!$    real(r8),         intent(out)   :: fgeo(lon0:lon1,lat0:lat1)
!!$    type(ESMF_Field), intent(inout) :: ESMF_mag, ESMF_geo
!!$    character(len=*), intent(in)    :: fname
!!$    !
!!$    ! Local:
!!$    integer :: j
!!$    character (len=8) :: fnames(1)
!!$    type(ESMF_Field) :: magfields(1)
!!$    real(r8),pointer,dimension(:,:) :: fptr
!!$
!!$    fgeo = finit
!!$    fnames(1) = fname
!!$    magfields(1) = ESMF_mag
!!$    !
!!$    ! Put fmag into ESMF mag field on mag source grid:
!!$    call edyn_esmf_set2d_mag(magfields,fnames,fmag,1, &
!!$         mlon0,mlon1,mlat0,mlat1)
!!$    !
!!$    ! Regrid to geographic destination grid, defining ESMF_geo:
!!$    call edyn_esmf_regrid(ESMF_mag,ESMF_geo,'mag2geo',2)
!!$    !
!!$    ! Put regridded geo field into pointer:
!!$    call edyn_esmf_get_2dfield(ESMF_geo,fptr,fname)
!!$    !      write(iulog,*) 'mag2geo: Max,min fptr = ',maxval(fptr),minval(fptr)
!!$    !
!!$    ! Transfer from pointer to output arg:
!!$    do j=lat0,lat1
!!$       fgeo(:,j) = fptr(:,j)
!!$    end do
!!$    !      write(iulog,*) 'mag2geo: max,min fmag = ',maxval(fmag),minval(fmag)
!!$    !      write(iulog,*) 'mag2geo: max,min fgeo = ',maxval(fgeo),minval(fgeo)
!!$  end subroutine mag2geo_2d
  !-----------------------------------------------------------------------
  subroutine rpt_ncerr(istat,msg)
    !
    ! Handle a netcdf lib error:
    !
    integer,         intent(in) :: istat
    character(len=*),intent(in) :: msg
    !
    write(iulog,"(/72('-'))")
    write(iulog,"('>>> Error from netcdf library:')")
    write(iulog,"(a)") trim(msg)

    write(iulog,"('istat=',i5)") istat
    write(iulog,"(72('-')/)")
    return
  end subroutine rpt_ncerr
  !-----------------------------------------------------------------------
  subroutine check_alloc(ierror, subname, varname, lonp1, latp1, ntimes, lw)
    use spmd_utils, only: masterproc
    integer,           intent(in) :: ierror
    character(len=*),  intent(in) :: subname
    character(len=*),  intent(in) :: varname
    integer, optional, intent(in) :: lonp1
    integer, optional, intent(in) :: latp1
    integer, optional, intent(in) :: ntimes
    integer, optional, intent(in) :: lw
    ! Local variable
    character(len=256) :: errmsg

    if (ierror /= 0) then
       write(errmsg, '(">>> ",a,": error allocating ",a)')                   &
            trim(subname), trim(varname)
       if (present(lonp1)) then
          write(errmsg(len_trim(errmsg)+1:), '(", lonp1 = ",i0)') lonp1
       end if
       if (present(latp1)) then
          write(errmsg(len_trim(errmsg)+1:), '(", latp1 = ",i0)') latp1
       end if
       if (present(ntimes)) then
          write(errmsg(len_trim(errmsg)+1:), '(", ntimes = ",i0)') ntimes
       end if
       if (present(lw)) then
          write(errmsg(len_trim(errmsg)+1:), '(", lw = ",i0)') lw
       end if
       if (masterproc) then
          write(iulog, *) trim(errmsg)
       end if
       call endrun(trim(errmsg))
    end if

  end subroutine check_alloc

  !-----------------------------------------------------------------------
  subroutine check_ncerr(istat, subname, msg)
    use pio, only: pio_noerr
    !
    ! Handle a netcdf lib error:
    !
    integer,          intent(in) :: istat
    character(len=*), intent(in) :: subname
    character(len=*), intent(in) :: msg
    !
    ! Local variable
    character(len=256) :: errmsg
    !
    if (istat /= pio_noerr) then
       write(iulog,"(/72('-'))")
       write(iulog,"('>>> Error from netcdf library:')")
       write(iulog,"(a,': Error getting ',a)") trim(subname), trim(msg)

       write(iulog,"('istat=',i5)") istat
       write(iulog,"(72('-')/)")
       write(errmsg, '("NetCDF Error in ",a,": ",2a,", istat = ",i0)')        &
            trim(subname), 'Error getting ', trim(msg), istat
       call endrun(trim(errmsg))
    end if
  end subroutine check_ncerr

#endif

end module mag2geo_mod
