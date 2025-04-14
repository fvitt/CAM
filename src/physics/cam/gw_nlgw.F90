module gw_nlgw

!
! This module predicts gravity wave forcings via PyTorch NNs trained to include non-local gravity wave effects
!

use gw_utils, only: r8
use ppgrid,   only: pver !vertical levels
use physics_types,  only: physics_state, physics_ptend
use spmd_utils,     only: mpicom, mstrid=>masterprocid, masterproc, mpi_real8
use cam_abortutils, only: endrun
use cam_logfile,    only: iulog
use physconst, only: cappa

use ftorch

implicit none

public :: gw_nlgw_dp_ml, gw_nlgw_dp_init, gw_nlgw_dp_finalize

private

integer, parameter :: p0 = 100000 ! 1000 hPa (Pa)
integer, parameter :: num_levels = 122 ! From WACCM

type(torch_model) :: nlgw_model ! pytorch model

integer :: ncol ! number of vertical columns

real(r8), dimension(:), allocatable :: &
  lat,     &! latitude (radians)
  lon,     &! longitude (radians)
  ps,      &! surface pressure
  phis      ! surface geopotential
real(r8), dimension(:,:), allocatable :: &
  u,       &! zonal wind (m/s)
  v,       &! meridional wind (m/s)
  omega,   &! vertical pressure velocity (Pa/s)
  t,       &! temperature (K)
  theta,   &! potential temperature (K)
  pmid      ! midpoint pressure (Pa)

real(r8), dimension(:,:), allocatable :: &
  utgw,    &! zonal wind tendency
  vtgw      ! meridional wind tendency

real(r8), dimension(:,:), allocatable, target :: net_inputs
real(r8), dimension(:,:), allocatable, target :: net_outputs

! normalisation means and std devs
real(r8) :: u_mean, v_mean, omega_mean, theta_mean, lat_mean, lon_mean
real(r8) :: u_std, v_std, omega_std, theta_std, lat_std, lon_std

real(r8) :: utgw_mean, vtgw_mean
real(r8) :: utgw_std, vtgw_std

real(r8) :: era5_ak = (

contains

!==========================================================================

subroutine gw_nlgw_dp_ml(state_in, ptend)

  ! inputs
  type(physics_state), intent(in) :: state_in
  ! outputs
  type(physics_ptend), intent(inout) :: ptend

  !---------------------------Local storage-------------------------------
  integer :: i

  type(torch_tensor) :: tensor_in(1), tensor_out(1)
  integer :: ninputs = 1, noutputs = 1
  integer, dimension(2) :: layout = [1 , 2]

  ncol = state_in%ncol

  allocate(lat(ncol))
  allocate(lon(ncol))
  allocate(ps(ncol))
  allocate(phis(ncol))
  allocate(u(ncol,pver))
  allocate(v(ncol,pver))
  allocate(t(ncol,pver))
  allocate(pmid(ncol,pver))
  allocate(theta(ncol,pver))
  allocate(omega(ncol,pver))

  allocate(utgw(ncol,pver))
  allocate(vtgw(ncol,pver))

  allocate(net_inputs(ncol, 4*num_levels+3))
  allocate(net_outputs(ncol, 2*num_levels))

  ! dims = (ncol) TODO check ncol size vs gw_drag
  lat = state_in%lat
  lon = state_in%lon
  ps = state_in%ps
  phis = state_in%phis

  ! dims = (ncol, pver)
  u = state_in%u
  v = state_in%v
  t = state_in%t
  pmid = state_in%pmid
  theta = t * (p0 / pmid) ** cappa
  omega = state_in%omega

  ! Normalise and construct the input
  call normalise_data()
  call construct_input()

  ! send all columns from this process
  call torch_tensor_from_array(tensor_in(1), net_inputs, layout, torch_kCUDA)
  call torch_tensor_from_array(tensor_out(1), net_outputs, layout, torch_kCPU)

  ! Run net forward on data
  call torch_model_forward(nlgw_model, tensor_in, tensor_out)

  ! Extract and denormalise outputs
  call extract_output()
  call denormalise_data()

!  utgw = -d(u_flux)/dp
!  vtgw = -d(v_flux)/dp

  ! update the tendencies
  ptend%u(:ncol,:) = ptend%u(:ncol,:) + utgw(:,:)
  ptend%v(:ncol,:) = ptend%v(:ncol,:) + vtgw(:,:)

  ! Clean up the tensors
  call torch_delete(tensor_in)
  call torch_delete(tensor_out)

  deallocate(lat)
  deallocate(lon)
  deallocate(ps)
  deallocate(phis)
  deallocate(u)
  deallocate(v)
  deallocate(t)
  deallocate(pmid)
  deallocate(theta)
  deallocate(omega)

  deallocate(utgw)
  deallocate(vtgw)

  deallocate(net_inputs)
  deallocate(net_outputs)

end subroutine gw_nlgw_dp_ml


subroutine gw_nlgw_dp_init(neural_net_path, norms_path)

  character(len=132), intent(in) :: neural_net_path  ! Filepath to PyTorch Torchscript net
  character(len=132), intent(in) :: norms_path       ! Filepath to NetCDF normalisation weights

  ! Load the convective drag net from TorchScript file
  call torch_model_load(nlgw_model, neural_net_path, device_type=torch_kCUDA, device_index=0)
  ! read in normalisation weights
  call read_norms()

  if (masterproc) then
     write(iulog,*)'gw_convect_net loaded from: ', neural_net_path
     ! write(iulog,*)'Normalisation weights loaded from: ', norms_path
  endif

end subroutine gw_nlgw_dp_init


subroutine gw_nlgw_dp_finalize()

  deallocate(net_inputs)
  deallocate(net_outputs)
  ! free model memory
  call torch_delete(nlgw_model)

end subroutine gw_nlgw_dp_finalize


subroutine read_norms()

  ! use netcdf
  ! use error_messages, only: handle_ncerr

  ! character(len=132), intent(in) :: norms_path  ! Filepath to NetCDF normalisation weights

  ! integer :: ncid, varid, retva, ierr
  ! character(len=*), parameter :: sub = 'gw_nlgw/F90 read_norms: '

  ! Load normalisation weights from file in master process then broadcast
  ! if (masterproc) then
  !   ! Open the NetCDF file
  !   call handle_ncerr( nf90_open(trim(norms_path), NF90_NOWRITE, ncid), &
  !                      "Error opening NetCDF norms file in gw_ml.F90")

  !   ! We do not need to read in dimensions here as we assume inputs match the grid.

  !   ! Read in variables (means and deviations).
  !   call handle_ncerr( nf90_inq_varid(ncid, 'U_mean', varid), &
  !                      "Error getting U_mean varid from NetCDF Norms file in gw_ml.F90")
  !   call handle_ncerr( nf90_get_var(ncid, varid, u_mean), &
  !                      "Error getting U_mean varid from NetCDF Norms file in gw_ml.F90")

  ! endif

  ! Broadcast normalisation variables to other processes
  ! call mpi_bcast(utgw_std, pver, mpi_real8, mstrid, mpicom, ierr)
  ! if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: utgw_std from gw_ml.F90")

  ! TODO
  ! - remove hardcoded means/std devs
  ! - replace with netcdf load and mpi broadcast
  ! - verify with Aman that these are correct
  ! - verify that the denormalization is correctly applied

  lat_mean = 0._r8
  lon_mean = 0._r8
  u_mean = 6.395471175756457_r8
  v_mean = 0.020313991225046_r8
  theta_mean = 0._r8
  omega_mean = 0.0016040905945274022_r8

  lat_std = 90._r8
  lon_std = 360._r8
  u_std = 3._r8 * 22.175504140184618_r8
  v_std = 3._r8 * 9.84143148277375_r8
  theta_std = 1000._r8
  omega_std = 0.017021397434040318_r8

  utgw_mean = -0.0005112474139891424_r8
  vtgw_mean = -0.0002982954242187403_r8
  utgw_std = 0.0050768547492663395_r8
  vtgw_std = 0.003792741148955207_r8

end subroutine read_norms

subroutine normalise_data()

  lat = (lat-lat_mean)/lat_std
  lon = (lon-lon_mean)/lon_std
  phis = phis / 50000._r8
  u = (u-u_mean)/u_std
  v = (v-v_mean)/v_std
  theta = (theta-theta_mean)/theta_std
  omega = (omega-omega_mean)/omega_std
  omega = omega ** (1.0/3.0) ! cube root of omega

  ! TODO :: currently there is no scaling for phis?

  print * , "min/max lat = " , minval(lat)   , " : " , maxval(lat)
  print * , "min/max lon = " , minval(lon)   , " : " , maxval(lon)
  print * , "min/max phi = " , minval(phis)  , " : " , maxval(phis)
  print * , "min/max u   = " , minval(u)     , " : " , maxval(u)
  print * , "min/max v   = " , minval(v)     , " : " , maxval(v)
  print * , "min/max the = " , minval(theta) , " : " , maxval(theta)
  print * , "min/max ome = " , minval(omega) , " : " , maxval(omega)

end subroutine normalise_data

subroutine construct_input()

  integer :: idx_beg, idx_end

  net_inputs(:,1) = lat
  net_inputs(:,2) = lon
  net_inputs(:,3) = phis

  call interp(u, u_interp)

  idx_beg = idx_end + 1
  idx_end = idx_beg + num_levels
  net_inputs(:,idx_beg:idx_end) = u
  idx_beg = idx_end + 1
  idx_end = idx_beg + num_levels
  net_inputs(:,idx_beg:idx_end) = v
  idx_beg = idx_end + 1
  idx_end = idx_beg + num_levels
  net_inputs(:,idx_beg:idx_end) = theta
  idx_beg = idx_end + 1
  idx_end = idx_beg + num_levels
  net_inputs(:,idx_beg:idx_end) = omega

!  u(0km:80km) 93 level -> input(0:50km) 122 levels xi0  = 0 xi121 = 50
! TODO add interpolation (msg Aman about pressure level grids)

end subroutine construct_input

subroutine extract_output()

  u_flux(:, :) = net_outputs(:,:num_levels)
  v_flux(:, :) = net_outputs(:,num_levels+1:)

end subroutine extract_output

subroutine denormalise_data()

  u_flux = u_flux**3 * u_flux_std + u_flux_mean
  v_flux = v_flux**3 * v_flux_std + v_flux_mean

end subroutine denormalise_data

end module gw_nlgw
