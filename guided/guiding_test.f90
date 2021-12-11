program dmc
  use mod_dftbp
  use parameters
  implicit real(8) (a-h,o-z)
  !dftbp parameters:------------------------------------------------------------
  real*16 e_path
  character(len=20) :: filetype
  character(len=100) :: k_point_cp_string,filename,file_result, file_w_path
  character(len=100):: dftb_path
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
  integer :: energy_count = 0, a !count the number of energy calculations
  real(8), allocatable, dimension(:,:) :: energies, fxyz, rxyz
  real(8), allocatable, dimension(:) :: r_temp, mu_w

  !Simulation parameters:-------------------------------------------------------
  no_walkers_start = 20 !initial number of walkers
  maxWalkers = 1000
  delta_t = 0.005d0 !does not work with 0.01 and above
  n_steps = 50 !number of simulation steps
  a = 1
  damping = 0.1d0 !damping constant for the update of et in every iteration
  reduction_thereshold = 0.2 !reduce the damping constant if the error is below this thereshold
  et = -5.62 !inital guess for the energy

  print*, "delta_t = ", delta_t

  call read_hessian(e0_theory, .FALSE.)

  !initialzations and memory allocations----------------------------------------
  dftb_path = "dftb_temp_folder"
  call setup_dftbp(dftb_path, 0)

  allocate(fxyz(3,nat))
  allocate(rxyz(3,nat))
  allocate(r_temp(dim))
  allocate(mu_w(dim))
  mu = 0.d0

  allocate(energies(100,dim))
  do k=1,dim
    mu_w = mu
    ener = 0
    mu_w(k) = mu_w(k)-50*0.1/Bohr_Ang
    do j=1,100
      r_temp = matmul(ev,mu_w)
      rxyz = reshape(r_temp,shape(rxyz)) + rxyz0
      call energyandforces(nat, rxyz, fxyz, ener, alat, deralat, atomnames, 0)
      mu_w(k) = mu_w(k) + 0.1/Bohr_Ang
      energies(j,k) = ener
    end do
    print*, "complete dimension ",k
  end do


  !output final positions of the walkers
  open(unit=13, file="ener", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,dim
    do j=1,100
      write(unit=13, fmt=*) energies(j, i)
    end do
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"

  print*, "all done"

end program dmc
