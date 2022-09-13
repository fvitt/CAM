module zmean_fields
  use shr_kind_mod,  only: r8=>SHR_KIND_R8
  use ppgrid,        only: begchunk, endchunk, pcols, pver
  use Zonal_Mean,    only: ZonalMean_t

  use namelist_utils,only: find_group_name
  use spmd_utils,    only: masterproc, mpi_integer, masterprocid, mpicom
  use cam_logfile,   only: iulog
  use cam_abortutils,only: endrun

  use spmd_utils, only: masterproc

  implicit none

  type(ZonalMean_t) :: ZMobj

  integer :: ZonalNbasis = -huge(1)

contains

  subroutine zmean_fields_readnl(nlfile)
    character(len=*), intent(in) :: nlfile
    integer :: ierr, unitn

    character(len=*), parameter :: prefix = 'zmean_fields_readnl: '

    namelist /zmean_fields_opts/ ZonalNbasis

    ! Read in namelist values
    !------------------------
    if(masterproc) then
       open(newunit=unitn, file=trim(nlfile), status='old')
       call find_group_name(unitn, 'zmean_fields_opts', status=ierr)
       if(ierr == 0) then
          read(unitn,zmean_fields_opts,iostat=ierr)
          if(ierr /= 0) then
             call endrun(prefix//'ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    call MPI_bcast(ZonalNbasis, 1, mpi_integer, masterprocid, mpicom, ierr)

    if (masterproc) then
       write(iulog,*) 'zmean_fields_readnl... ZonalNbasis: ',ZonalNbasis
    endif

  end subroutine zmean_fields_readnl

  subroutine zmean_fields_init

    call ZMobj%init(ZonalNbasis)

  end subroutine zmean_fields_init

  function zmean_3d( fld ) result(fldzm)


    real(r8), intent(in) :: fld(pcols,pver,begchunk:endchunk)

    real(r8) :: fldzm(pcols,pver,begchunk:endchunk)

    real(r8) :: Zonal_Bamp3d(ZonalNbasis,pver)

    call ZMobj%calc_amps(fld,Zonal_Bamp3d)
    call ZMobj%eval_grid(Zonal_Bamp3d,fldzm)

  end function zmean_3d

end module zmean_fields
