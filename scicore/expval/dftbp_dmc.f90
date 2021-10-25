program dmc
  implicit real(8) (a-h,o-z)
  logical :: debug_extended

  damping = 0.01d0 !damping constant for the update of et in every iteration
  no_walkers_start = 500 !initial number of walkers
  delta_t = 10.0d0
  n_steps = 1000 !number of simulation steps
  debug_extended = .FALSE.

  print*, "delta_t = ", delta_t
  print*, "damping = ", damping
  print*, "no_walkers_start = ", no_walkers_start
  print*, "n_steps = ", n_steps

  call do_simulation(no_walkers_start, delta_t, n_steps, damping, debug_extended)

end program dmc

  subroutine do_simulation(no_walkers_start, delta_t, n_steps, damping, debug_extended)
  use OMP_LIB
  use mod_dftbp
  use dftbp_walker
  use parameters
  implicit real(8) (a-h,o-z)
  !dftbp parameters:------------------------------------------------------------
  real*16 e_path
  character(len=20) :: filetype
  character(len=100) :: k_point_cp_string,filename,file_result, file_w_path, folder_path
  !!Integers for system clock timing
  integer*8 :: c,cr,cm,TimeStart,TimeEnd
  !! position step size
  real*8 :: alpha_pos = 1.d0
  !! lattice step size
  real*8 :: alpha_lat = 2.d0
  !! minimal step size
  real*8 :: alpha0 = 1.d-3
  !! Set true when forces should be preconditioned
  logical :: for_pre = .FALSE.!.FALSE. TRUE
  !! force tolerance
  real*8 :: fnrmtol = 1.d-5
  !! lattice derivative tolerance
  real*8 :: lat_tol = 1.d-5
  !! set true if debug information should be written.
  logical :: debug = .TRUE.
  !! set true if every configuration from each iteration should be written to disk
  logical :: write_configurations = .false.
  !! max number of iterations
  integer :: nit = 3000
  !iterations optimizer
  integer :: it_opt
  !simulation variables:--------------------------------------------------------
  integer :: energy_count = 0, a, thread_num !count the number of energy calculations
  type(walker), allocatable, dimension(:) :: walkers, new_walkers
  type(walker) :: half_walker
  real(8), allocatable, dimension(:) :: all_e0, normal_rand, r
  real(8), allocatable, dimension(:,:) :: fxyz, rxyz
  logical, parameter :: adjust_flag = .FALSE. !.TRUE. = with damping adjusting
  logical :: debug_extended, half_walker_flag = .FALSE.

  !Simulation parameters:-------------------------------------------------------
  reduction_thereshold = 0.2 !reduce the damping constant if the error is below this thereshold
  maxWalkers = no_walkers_start * 100 !kill some walker if the number of walkers exeeds this number

  call read_hessian(e0_theory, debug_extended)
  allocate(fxyz(3,nat))
  allocate(rxyz(3,nat))
  allocate(normal_rand(dim))
  allocate(r(dim))

  !initialzations and memory allocations----------------------------------------
  old_no_walkers = no_walkers_start
  no_walkers = no_walkers_start
  allocate(all_e0(n_steps))
  allocate(walkers(no_walkers))
  folder_path = "dftb_temp_folder"
  OMP_NUM_THREADS = OMP_GET_MAX_THREADS()
  print*, "Number of threads = ", OMP_GET_MAX_THREADS()
  do i=0,OMP_GET_MAX_THREADS()-1
    call setup_dftbp(folder_path, i)
  end do
  call energyandforces(nat, rxyz0, fxyz, etot, alat, deralat, atomnames, 0)
  call initalizeWalkers(walkers,no_walkers,etot)
  print*, "Equilibrium position energy: ",etot
  print*, "E_0 from the harmonic approximation: ", e0_theory
  et = etot + e0_theory !inital guess for the energy
  print*, "Initial guess for the energy: ", et
  e0 = et

  !write the progress of the simulation
  open(unit=12, file="expval", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 12"

  !expectation value------------------------------------------------------------
  exp_val_harm = 0.d0
  exp_val_pes = 0.d0
  exp_val_kin = 0.d0
  do n = 1,10000
    call normal(normal_rand,dim)
    r = normal_rand*sqrt(sigma2)
    !exp_val_pes = exp_val_pes + e_pes(r, 1) + el_kin(r)
    !exp_val_harm = exp_val_harm + e_harm(r) + el_kin(r)
    exp_val_pes = exp_val_pes + e_pes(r, 0) - etot
    exp_val_harm = exp_val_harm + e_harm(r)
    exp_val_kin = exp_val_kin + el_kin(r)

    if (modulo(n,10) .eq. 0) then
      print*, n
      write(unit=12, fmt=*) exp_val_harm/n, exp_val_pes/n, exp_val_kin/n
    end if
  end do

  close(unit=12, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 12"

  print*, "all done!!!"

end subroutine do_simulation


real(8) function e_harm(r)
  use parameters
  implicit none
  real(8), dimension(dim), intent(in) :: r
  integer :: i

  e_harm = 0.d0
  do i=1,dim
    e_harm = e_harm  + 0.5*eigenvalues_full(i+6)*r(i)**2
  end do

end function e_harm

real(8) function e_pes(r, thread_num)
  use OMP_LIB
  use mod_dftbp
  use parameters
  implicit none
  real(8), dimension(dim), intent(in) :: r
  real(8), dimension(3,nat) :: rxyz, fxyz
  real(8) :: epot
  integer :: i, thread_num

  rxyz = reshape(matmul(ev,r),shape(rxyz)) + rxyz0 !convert to kartesian corrdinated for energy calculation

  call energyandforces(nat, rxyz, fxyz, epot, alat, deralat, atomnames, thread_num)

  e_pes = epot

end function e_pes
