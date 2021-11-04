program dmc
  implicit real(8) (a-h,o-z)
  logical :: debug_extended

  damping = 0.01d0 !damping constant for the update of et in every iteration
  no_walkers_start = 2000 !initial number of walkers
  delta_t = 1.d-4 !discrete time step
  n_steps = 10000 !number of simulation steps
  debug_extended = .FALSE. !TRUE to show extended information such as the hessian and eigenvalues

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
  real(8), allocatable, dimension(:) :: all_e0
  real(8), allocatable, dimension(:,:) :: fxyz, rxyz
  logical, parameter :: adjust_flag = .FALSE. !.TRUE. = with damping adjusting
  logical :: debug_extended, half_walker_flag = .FALSE.

  !Simulation parameters:-------------------------------------------------------
  reduction_thereshold = 0.2 !reduce the damping constant if the error is below this thereshold
  maxWalkers = no_walkers_start * 100 !kill some walker if the number of walkers exeeds this number

  !initialzations and memory allocations----------------------------------------
  old_no_walkers = no_walkers_start
  no_walkers = no_walkers_start
  allocate(all_e0(n_steps))
  allocate(walkers(no_walkers))
  folder_path = "dftb_temp_folder"
  OMP_NUM_THREADS = OMP_GET_MAX_THREADS()
  print*, "Number of threads = ", OMP_GET_MAX_THREADS()
  do i=0,OMP_GET_MAX_THREADS()-1
    call setup_dftbp(folder_path, i) !set up all folders for the DFTB+ calculations
  end do
  call calculate_guiding_wf(e0_theory, debug_extended)
  allocate(fxyz(3,nat))
  allocate(rxyz(3,nat))
  !inital calculation for the energy E_min
  call energyandforces(nat, rxyz0, fxyz, e_min, alat, deralat, atomnames, 0)
  call initalizeWalkers(walkers,no_walkers,e_min)
  print*, "Equilibrium position energy: ", e_min
  print*, "E_0 from the harmonic approximation: ", e0_theory
  et = e_min + e0_theory !inital guess for the energy
  print*, "Initial guess for the energy: ", et
  e0 = et

  !write the progress of the simulation
  open(unit=12, file="et_noWalkers", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 12"
  write(unit=12, fmt=*) et, no_walkers, e0

  !Start of the simulation------------------------------------------------------
  do n=1,n_steps

    !propagate the walkers in parallel
    !$omp parallel do private(thread_num) shared(et, delta_t)
    do j=1,no_walkers
      !move the walker
      thread_num = OMP_GET_THREAD_NUM()
      call propagate(walkers(j), et, delta_t, thread_num, n)
    end do
    !$omp end parallel do

    !calculate the number of walkers after the branching process--------------
    projected_no_walkers = 0
    do j=1,no_walkers
      if (walkers(j)%weight .gt. 2) then !calculate the new number of walkers (multiplied by two)
        projected_no_walkers = projected_no_walkers + int(walkers(j)%weight)
      elseif ( walkers(j)%weight .lt. 0.5) then
        projected_no_walkers = projected_no_walkers + 0.5
      else
        projected_no_walkers = projected_no_walkers + 1
      end if
    end do

    new_no_walkers = floor(projected_no_walkers)

    if (new_no_walkers .lt. maxWalkers) then! normal case
      allocate(new_walkers(new_no_walkers))
      !split-joint algorithm----------------------------------------------------
      k = 1
      do j=1,no_walkers
        if (k .le. new_no_walkers) then
          if (walkers(j)%weight .gt. 2) then !split
            do i=1,int(walkers(j)%weight)
              if (k .le. new_no_walkers) then
                new_walkers(k) = walkers(j)
                new_walkers(k)%weight = walkers(j)%weight/int(walkers(j)%weight)
                k = k + 1
              end if
            end do
          elseif (walkers(j)%weight .lt. 0.5) then !join
            if (half_walker_flag) then
              half_walker_flag = .FALSE.
              call RANDOM_NUMBER(r)
              if (r .lt. walkers(j)%weight/(walkers(j)%weight+half_weight)) then
                new_walkers(k) = walkers(j)
              else
                new_walkers(k) = half_walker
              end if
              new_walkers(k)%weight = walkers(j)%weight + half_weight
              k = k + 1
            else
              half_walker_flag = .TRUE.
              half_walker = walkers(j)
              half_weight = walkers(j)%weight
            end if
          else !normal case
            new_walkers(k) = walkers(j)
            k = k + 1
          end if
        end if
      end do

      deallocate(walkers)
      allocate(walkers(new_no_walkers))
      walkers = new_walkers
      deallocate(new_walkers)

    else !Too many walkers. Choose no_walkers_start walkers randomly to survive
      print*, "Too many walkers!"
      new_no_walkers = no_walkers_start
      allocate(new_walkers(new_no_walkers))
      do j=1,new_no_walkers
        call RANDOM_NUMBER(r)
        index = int(r*new_no_walkers+1)
        new_walkers(j) = walkers(index)
      end do
      deallocate(walkers)
      allocate(walkers(new_no_walkers))
      walkers = new_walkers
      deallocate(new_walkers)
      no_walkers = no_walkers_start
      old_no_walkers = no_walkers_start
    end if

    if (new_no_walkers .eq. 0) then !to keep the simulation going if the population dies out
      print*, "Walkers died out!"
      no_walkers = no_walkers_start
      new_no_walkers = no_walkers_start
      old_no_walkers = no_walkers_start
      deallocate(walkers)
      allocate(walkers(no_walkers))
      call initalizeWalkers(walkers,no_walkers,e_min)
    end if

    !adjust et with formula form the master exam paper (A=1)
    e0 = 0
    sum_weight = 0
    do j=1,new_no_walkers
      e0 = e0 + walkers(j)%etot * walkers(j)%weight
      sum_weight = sum_weight + walkers(j)%weight
    end do
    e0 = e0/sum_weight !best energy guess
    et = e0 - damping*log(sum_weight/real(no_walkers_start))

    all_e0(n) = e0
    no_walkers = new_no_walkers
    energy_count = energy_count + no_walkers

    !adjust the damping term (set flag to TRUE)
    if (adjust_flag) then
      if (n .gt. 10) then
        if (abs(sum(all_e0(n-10:n-1:1))/10 - e0) .le. reduction_thereshold) then
          damping = damping*0.9
          reduction_thereshold = reduction_thereshold*0.1
          print*, "damping reduced!"
        end if
      end if
    end if

    write(unit=12, fmt=*) et, no_walkers, e0

    if (modulo(n,10) .eq. 0) then !write every 100 steps to the out file
      write(*,'(A,I5,A,I5,A)') "Step ", n, " of ", n_steps, "--------------------------------"
      print*, "number of walkers = ", no_walkers
      print*, "et = ", et
      print*, "e0 = ", e0
    end if
  end do
  !end of the simulation--------------------------------------------------------

  !calculation of the ZPE, ignoring the first 500 steps
  zpe = sum(all_e0(n_steps-min(500,n_steps)+1:n_steps))/min(500,n_steps) - e_min
  print*, "ZPE = ", zpe

  print*, "Number of energy calculations = ", energy_count

  close(unit=12, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 12"

  !output final positions of the walkers
  open(unit=13, file="walker_positions_kartesian", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,no_walkers
    call write_walker_kartesian(walkers(i), 13)
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"

  !energy and weight for restarting
  open(unit=14, file="walker_energies", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 14"
  do i=1,no_walkers
    write(unit=14, fmt=*) walkers(i)%etot, walkers(i)%weight
  end do
  close(unit=14, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 14"

  open(unit=15, file="walker_positions_hessian", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 15"
  do i=1,no_walkers
    call write_walker_hessian(walkers(i), 15)
  end do
  close(unit=15, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 15"

end subroutine do_simulation
