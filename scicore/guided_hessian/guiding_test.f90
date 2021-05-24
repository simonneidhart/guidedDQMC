program dmc
  use mod_dftbp
  use dftbp_walker
  use parameters
  implicit real(8) (a-h,o-z)
  interface
    subroutine initalizeWalkers(walkers,no_walkers,etot)
      use dftbp_walker
      use parameters
      implicit real(8) (a-h,o-z)
      type(walker), dimension(no_walkers) :: walkers
    end subroutine initalizeWalkers
  end interface
  !dftbp parameters:------------------------------------------------------------
  real*16 e_path
  character(len=20) :: filetype
  character(len=100) :: k_point_cp_string,filename,file_result, file_w_path
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
  type(walker), allocatable, dimension(:) :: walkers, new_walkers
  type(walker) :: half_walker
  real(8), allocatable, dimension(:,:) :: energies
  real(8), dimension(3,nat) :: fxyz, rxyz
  real(8), dimension(24) :: r_temp
  real(8), dimension(dim) :: mu_w

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

  call read_hessian()

  !initialzations and memory allocations----------------------------------------

  call setup_dftbp("dftb_temp_folder/")

  allocate(energies(100,dim))
  do k=1,dim
    mu_w = mu
    ener = 0
    mu_w(k) = mu_w(k)-50*0.1
    do j=1,100
      r_temp = matmul(ev,mu_w)
      rxyz = reshape(r_temp,shape(rxyz)) + rxyz0
      call energyandforces(nat, rxyz, fxyz, ener, alat, deralat, atomnames)
      mu_w(k) = mu_w(k) + 0.1
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

  print*, "all complete"

end program dmc

!uniformly initalizes walker positions in the given box.------------------------
subroutine initalizeWalkers(walkers,no_walkers,etot)
  use dftbp_walker
  use parameters
  implicit real(8) (a-h,o-z)
  type(walker), dimension(no_walkers) :: walkers

  do i=1,no_walkers
    call create_walker(walkers(i), mu, etot, 1.d0)
  end do

end subroutine initalizeWalkers

!creates an array of 3 x nat normal distributed random numbers using----
!the Box-Muller algorithm.
subroutine normal(rand, dim)
  implicit real(8) (a-h,o-z)
  integer :: dim
  real(8), dimension(dim) :: rand
  pi = 4.d0*atan(1.d0)

  n = dim
  if (mod(n,2) .eq. 0) then !even number of random numbers wanted
    call RANDOM_NUMBER(rand)
    do i=1,n/2
      x = rand(2*i-1)
      y = rand(2*i)
      rand(2*i-1) = sqrt(-2.d0*log(x))*cos(2.d0*pi*y)
      rand(2*i) = sqrt(-2.d0*log(x))*sin(2.d0*pi*y)
    end do
  else !odd number of random numbers wanted
    call RANDOM_NUMBER(rand)
    do i=1,n/2
      x = rand(2*i-1)
      y = rand(2*i)
      rand(2*i-1) = sqrt(-2.d0*log(x))*cos(2.d0*pi*y)
      rand(2*i) = sqrt(-2.d0*log(x))*sin(2.d0*pi*y)
    end do
      rand(n) = sqrt(-2.0*log(rand(n)))*cos(2*pi*rand(n+1))
  end if

end subroutine normal
