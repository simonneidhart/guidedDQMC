program dftbp_dqmc_main
  implicit none
  logical :: debug_extended, restart, grad_descent
  character(len=100) :: label
  integer :: ios, no_walkers_start, n_steps
  real(8) :: delta_t, damping

  !Reading the dqmc.in input file for the simulation parameters.
  open(unit=11, file="dqmc.in", iostat=ios, action="read")
  if ( ios /= 0 ) stop "Error opening file unit 11"
  read(11,*) label, damping
  read(11,*) label, no_walkers_start
  read(11,*) label, delta_t
  read(11,*) label, n_steps
  read(11,*) label, debug_extended
  read(11,*) label, restart
  read(11,*) label, grad_descent
  close(unit=11, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 11"

  !Writing back the parameters
  print*, "delta_t = ", delta_t
  print*, "damping = ", damping
  print*, "no_walkers_start = ", no_walkers_start
  print*, "n_steps = ", n_steps
  print*, "debug_extended = ", debug_extended
  print*, "continue_run = ", restart
  print*, "grad_descent = ", grad_descent

  call do_simulation(no_walkers_start, delta_t, n_steps, damping, debug_extended, restart, grad_descent)

end program dftbp_dqmc_main

  subroutine do_simulation(no_walkers_start, delta_t, n_steps, damping, debug_extended, restart, grad_descent)
  use OMP_LIB
  use mod_dftbp
  use dftbp_walker
  use dftbp_parameters
  implicit none
  character(len=100) :: folder_path
  integer :: energy_count, thread_num, n_steps, n_lines, n_start, i, j, k, n, ios, index
  integer :: max_walkers, no_walkers_start, no_walkers, new_no_walkers, old_no_walkers
  integer, allocatable, dimension(:) :: n_acc
  type(walker), allocatable, dimension(:) :: walkers, new_walkers
  type(walker) :: half_walker
  real(8) :: r, e0, e0_theory, et, delta_t, damping, e_min, zpe, projected_no_walkers
  real(8) :: half_weight, sum_weight, delta_t_eff, p_acc
  real(8), allocatable, dimension(:) :: all_e0
  real(8), allocatable, dimension(:,:) :: fxyz, rxyz, et_no_walkers_old
  logical :: debug_extended, restart, grad_descent, half_walker_flag = .FALSE.

  !kill some walker if the number of walkers exeeds this number
  max_walkers = no_walkers_start * 100
  energy_count = 0 !count the energy calculations

  !create the dftbp folders-----------------------------------------------------
  folder_path = "dftb_temp_folder"
  print*, "Number of threads = ", OMP_GET_MAX_THREADS()
  do i=0,OMP_GET_MAX_THREADS()-1
    call setup_dftbp(folder_path, i) !set up all folders for the DFTB+ calculations
  end do

  call calculate_guiding_wf(e0_theory, debug_extended, grad_descent)

  !inital calculation for the energy E_min
  allocate(fxyz(3,nat))
  allocate(rxyz(3,nat))
  call energyandforces(nat, rxyz0, fxyz, e_min, alat, deralat, atomnames, 0)
  print*, "Minimum energy on the PES E_min: ", e_min
  print*, "E_0 from the harmonic approximation: ", e0_theory

  if (restart .eq. .TRUE.) then !continue a previous run------------------------

    n_lines = 0 !figure out how many iterations were done before
    open(unit=12, file="et_noWalkers_e0.out", iostat=ios, action="read")
    if ( ios /= 0 ) stop "Error opening file unit 12"
    do
      read(12,*,iostat=ios)
      if (ios/=0) exit
      n_lines = n_lines + 1
    end do
    close(unit=12, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 12"

    allocate(et_no_walkers_old(n_lines,3))

    open(unit=12, file="et_noWalkers_e0.out", iostat=ios, action="read")
    if ( ios /= 0 ) stop "Error opening file unit 12"
    print*, n_lines, "old iterations detected"
    n_start = n_lines + 1
    do i=1,n_lines
      read(12, *) ( et_no_walkers_old(i, j), j = 1, 3) !read old data
    end do
    close(unit=12, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 12"

    open(unit=12, file="et_noWalkers_e0.out", iostat=ios, action="write")
    if ( ios /= 0 ) stop "Error opening file unit 12"
    do i=1,n_lines
      write(unit=12, fmt=*) ( et_no_walkers_old(i, j), j = 1, 3) !write old data
    end do

    !set the values to the ending values of the prevous run
    no_walkers = et_no_walkers_old(n_lines,2)
    old_no_walkers = no_walkers
    et = et_no_walkers_old(n_lines,1)
    e0 = et_no_walkers_old(n_lines,3)
    allocate(walkers(no_walkers))
    call read_restart(walkers, no_walkers)! read the prevous postions, energies
    !and weights
    n_steps = n_steps + n_lines
    write(*,*) "Restart complete: ", no_walkers, " walkers, ", et, " et"

  else !start a new run---------------------------------------------------------
    old_no_walkers = no_walkers_start
    no_walkers = no_walkers_start
    allocate(walkers(no_walkers))
    call initalize_walkers(walkers,no_walkers,e_min)
    n_start = 1
    et = e_min + e0_theory !inital guess for the energy
    print*, "Initial guess for the energy: ", et
    e0 = et
    open(unit=12, file="et_noWalkers_e0.out", iostat=ios, action="write")
    if ( ios /= 0 ) stop "Error opening file unit 12"
  end if

  !Start of the simulation------------------------------------------------------
  print*, " "
  allocate(all_e0(n_steps))
  allocate(n_acc(OMP_GET_MAX_THREADS()))
  delta_t_eff = delta_t
  do n=n_start,n_steps

    n_acc = 0
    !propagate the walkers in parallel
    !$omp parallel do private(thread_num) shared(walkers, et, delta_t_eff, n, n_acc)
    do j=1,no_walkers
      thread_num = OMP_GET_THREAD_NUM()
      call propagate(walkers(j), et, delta_t_eff, thread_num, n, n_acc(thread_num+1)) !move the walker
    end do
    !$omp end parallel do

    p_acc = sum(n_acc)/real(no_walkers) !acceptance probability over the last iteration
    delta_t_eff = delta_t*p_acc !effective delta_t

    !calculate the number of walkers after the branching process to allocate----
    projected_no_walkers = 0.d0
    do j=1,no_walkers
      !calculate the new number of walkers
      if (walkers(j)%weight .gt. 2) then
        projected_no_walkers = projected_no_walkers + int(walkers(j)%weight)
      elseif ( walkers(j)%weight .lt. 0.5) then
        projected_no_walkers = projected_no_walkers + 0.5d0
      else
        projected_no_walkers = projected_no_walkers + 1.d0
      end if
    end do
    new_no_walkers = floor(projected_no_walkers)

    !Too many walkers. Choose no_walkers_start walkers randomly to survive
    if (new_no_walkers .gt. max_walkers) then
      print*, "Too many walkers!------------------------------------------------"
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

    elseif (new_no_walkers .eq. 0) then !walkers died out. Initalize new walkers.
      print*, "Walkers died out!------------------------------------------------"
      no_walkers = no_walkers_start
      new_no_walkers = no_walkers_start
      old_no_walkers = no_walkers_start
      deallocate(walkers)
      allocate(walkers(no_walkers))
      call initalize_walkers(walkers,no_walkers,e_min)

    else! normal case (walker population stable)
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
          else !keep the walker
            new_walkers(k) = walkers(j)
            k = k + 1
          end if
        end if
      end do
    end if

    !copy the walker for the next iteration
    deallocate(walkers)
    allocate(walkers(new_no_walkers))
    walkers = new_walkers
    deallocate(new_walkers)

    !adjust et------------------------------------------------------------------
    e0 = 0
    sum_weight = 0
    do j=1,new_no_walkers
      e0 = e0 + walkers(j)%etot * walkers(j)%weight
      sum_weight = sum_weight + walkers(j)%weight
    end do
    e0 = e0/sum_weight !best energy guess
    et = e0 - damping*log(sum_weight/real(no_walkers_start))

    all_e0(n-n_start+1) = e0
    no_walkers = new_no_walkers
    energy_count = energy_count + no_walkers

    write(unit=12, fmt=*) et, real(no_walkers), e0! write the current iteration

    if (modulo(n,100) .eq. 0) then !write every 100 steps to the out file
      write(*,'(A,I6,A,I6,A)') "Step ", n, " of ", n_steps, "--------------------------------"
      print*, "number of walkers = ", no_walkers
      print*, "et = ", et
      print*, "e0 = ", e0
      print*, "p_acc", p_acc
    end if

    if (modulo(n,1000) .eq. 0) then !save walkers every 1000 steps in case of a crash

      open(unit=13, file="walker_positions_cartesian.out", iostat=ios, action="write")
      if ( ios /= 0 ) stop "Error opening file unit 13"
      do i=1,no_walkers
        call write_walker_cartesian(walkers(i), 13)
      end do
      close(unit=13, iostat=ios, status="keep")
      if ( ios /= 0 ) stop "Error closing file unit 13"

      open(unit=14, file="energies_weights.out", iostat=ios, action="write")
      if ( ios /= 0 ) stop "Error opening file unit 14"
      do i=1,no_walkers
        write(unit=14, fmt=*) walkers(i)%etot, walkers(i)%weight
      end do
      close(unit=14, iostat=ios, status="keep")
      if ( ios /= 0 ) stop "Error closing file unit 14"

      open(unit=15, file="walker_positions.out", iostat=ios, action="write")
      if ( ios /= 0 ) stop "Error opening file unit 15"
      do i=1,no_walkers
        call write_walker_modes(walkers(i), 15)
      end do
      close(unit=15, iostat=ios, status="keep")
      if ( ios /= 0 ) stop "Error closing file unit 15"

      print*, "Writing backup files complete."
    end if
  end do !end of the simulation-------------------------------------------------

  close(unit=12, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 12"

  !output final positions of the walkers
  open(unit=13, file="walker_positions_cartesian.out", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,no_walkers
    call write_walker_cartesian(walkers(i), 13)
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"

  open(unit=15, file="walker_positions.out", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 15"
  do i=1,no_walkers
    call write_walker_modes(walkers(i), 15)
  end do
  close(unit=15, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 15"

  !energy and weight for restarting
  open(unit=14, file="energies_weights.out", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 14"
  do i=1,no_walkers
    write(unit=14, fmt=*) walkers(i)%etot, walkers(i)%weight
  end do
  close(unit=14, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 14"


  if (restart) then ! calculate the zpe over all steps of the current run
    zpe = sum(all_e0)/size(all_e0) - e_min
  else !calculation of the ZPE, ignoring the first 500 steps
    zpe = sum(all_e0(n_steps-min(500,n_steps)+1:n_steps))/min(500,n_steps) - e_min
  end if
  print*, "Number of energy calculations: ", energy_count
  print*, "ZPE = ", zpe

end subroutine do_simulation

!read all the info of the walkers from a previous run
subroutine read_restart(walkers, no_walkers)
  use dftbp_walker
  use dftbp_parameters
  implicit none
  integer :: i, j, k, ios, no_walkers
  type(walker), dimension(no_walkers) :: walkers
  real(8) :: weight, energy
  real(8), dimension(dim) :: r

  open(unit=16, file="walker_positions.out", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 16"
  open(unit=17, file="energies_weights.out", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 17"

  do k=1,size(walkers)
    read(17,*) energy, weight
    do i=1,dim
      read(16,*) r(i)
    end do
    call create_walker_restart(walkers(k), r, energy, weight)
  end do

  close(unit=16, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 16"
  close(unit=17, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 17"

end subroutine read_restart

!Initalize walkers in the equilibrium position with weight 1------------------
subroutine initalize_walkers(walkers,no_walkers,etot)
  use dftbp_walker
  implicit real(8) (a-h,o-z)
  type(walker), dimension(no_walkers) :: walkers

  do i=1,no_walkers
    call create_walker(walkers(i), etot, 1.d0)
  end do

end subroutine initalize_walkers
