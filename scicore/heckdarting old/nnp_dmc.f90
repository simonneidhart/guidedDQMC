program dmc
  use precision
  use atomicStructure
  use RuNNerInterface
  implicit real(dp) (a-h,o-z)
  logical :: debug_extended
  character(len=100) :: geom_filepath
  !real(8), allocatable, dimension(:) :: all_e0

  damping = 0.01_dp !damping constant for the update of et in every iteration
  no_walkers_start = 1000 !initial number of walkers
  delta_t = 10.0_dp
  n_steps = 2000 !number of simulation steps
  geom_filepath = "geom/1.ascii"
  debug_extended = .TRUE.

  print*, "delta_t = ", delta_t
  print*, "damping = ", damping
  print*, "no_walkers_start = ", no_walkers_start
  print*, "n_steps = ", n_steps

  call do_simulation(no_walkers_start, delta_t, n_steps, damping, geom_filepath, debug_extended)

  !print*, "Average E_0: ", sum(all_e0(50:n_steps))/(n_steps-50)

end program dmc

  subroutine do_simulation(no_walkers_start, delta_t, n_steps, damping, geom_filepath, debug_extended)
  use OMP_LIB
  use dftbp_walker
  use parameters
  use precision
  use atomicStructure
  use RuNNerInterface
  implicit real(dp) (a-h,o-z)
  !RuNNer parameters:------------------------------------------------------------
  character(len=*), parameter :: RuNNerPath = '/kernph/finjon00/FHMC-periodic/RuNNer-sockets/RuNNer/build/RuNNer.x'
  character(len=*), parameter :: RuNNerParams = '/kernph/finjon00/FHMC-periodic/train-nnp/6/predict-small/*'
  character(len=100) :: geom_filepath
  type(atStruct) :: ats
  type(atStruct), allocatable, dimension(:) :: all_ats
  real(dp) :: energy
  !simulation variables:--------------------------------------------------------
  integer :: energy_count = 0, a, thread_num !count the number of energy calculations
  type(walker), allocatable, dimension(:) :: walkers, new_walkers
  type(walker) :: half_walker
  real(dp), allocatable, dimension(:) :: all_e0
  real(dp), allocatable, dimension(:,:) :: fxyz
  logical, parameter :: adjust_flag = .FALSE. !.TRUE. = with damping adjusting
  logical :: debug_extended, half_walker_flag = .FALSE.

  !Simulation parameters:-------------------------------------------------------
  !damping = 0.01d0 !damping constant for the update of et in every iteration
  !no_walkers_start = 100 !initial number of walkers
  !maxWalkers = 100000
  !delta_t = 10.0d0 !does not work with 0.01 and above
  !n_steps = 200 !number of simulation steps
  reduction_thereshold = 0.2_dp !reduce the damping constant if the error is below this thereshold
  maxWalkers = no_walkers_start * 100

  call as_readAscii(trim(adjustl(geom_filepath)), ats, .TRUE.)
  allocate(fxyz(3,ats%nat))
  OMP_NUM_THREADS = OMP_GET_MAX_THREADS()
  print*, "Number of threads = ", OMP_GET_MAX_THREADS()
  allocate(globalRuNNerHandles(OMP_GET_MAX_THREADS()))
  allocate(all_ats(OMP_GET_MAX_THREADS()))
  do i=1,OMP_GET_MAX_THREADS()
    all_ats(i) = ats
    call startRuNNer(all_ats(i), RuNNerPath, RuNNerParams, globalRuNNerHandles(i))
  end do
  print*, "RuNNers started"
  call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, etot, fxyz)
  print*, "Energy calculation successfull: ", etot
  call gradient_descent(ats)
  call read_hessian(ats, e0_theory, debug_extended)

  !initialzations and memory allocations----------------------------------------
  old_no_walkers = no_walkers_start
  no_walkers = no_walkers_start
  allocate(all_e0(n_steps))
  allocate(walkers(no_walkers))
  call initalizeWalkers(walkers,no_walkers,etot)
  print*, "Equilibrium position energy: ",etot
  print*, "E_0 from the harmonic approximation: ", e0_theory
  et = etot + e0_theory !inital guess for the energy
  !et = etot
  print*, "Initial guess for the energy: ", et
  e0 = et

  open(unit=12, file="et_noWalkers", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 12"
  write(unit=12, fmt=*) et, no_walkers, e0

  !Start of the simulation------------------------------------------------------
  do n=1,n_steps

    !$omp parallel do private(thread_num) shared(et, delta_t)
    do j=1,no_walkers
      !move the walker
      thread_num = OMP_GET_THREAD_NUM() + 1
      call propagate(walkers(j), et, delta_t, thread_num, all_ats(thread_num))
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
      !if (k-1 .ne. new_no_walkers ) then
        !print*, "no_walkers wrong", k-1, new_no_walkers
      !end if
      deallocate(walkers)
      allocate(walkers(new_no_walkers))
      walkers = new_walkers
      deallocate(new_walkers)

    else !Too many walkers. Choose no_walkers_start walkers randomly to survive-------------
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
      call initalizeWalkers(walkers,no_walkers,etot)
    end if

    !v_bar = 0 !alternative for adjusting et (Mc Coy)
    !do j=1,new_no_walkers
    !v_bar = v_bar + walkers(j)%etot
    !end do
    !v_bar = v_bar/new_no_walkers
    !et = v_bar - damping*(new_no_walkers-no_walkers_start)/no_walkers_start
    !e0 = et

    !adjust et with formula form the master exam paper (A=1)
    e0 = 0._dp
    sum_weight = 0._dp
    do j=1,new_no_walkers
      e0 = e0 + walkers(j)%etot * walkers(j)%weight
      sum_weight = sum_weight + walkers(j)%weight
    end do
    e0 = e0/sum_weight !best energy guess
    et = e0 - damping*log(sum_weight/real(no_walkers_start))

    !a = int(n/100) + 1
    !a = 1
    !all_et(n) = et
    !if (modulo(n,a) .eq. 0) then
    !  et = all_et(n-a+1) - damping/(delta_t*a)*log(real(new_no_walkers)/old_no_walkers)
    !  old_no_walkers = new_no_walkers
    !end if

    all_e0(n) = e0
    no_walkers = new_no_walkers
    energy_count = energy_count + no_walkers

    !adjust the damping term (set flag to true)

    if (adjust_flag) then
      if (n .gt. 10) then
        if (abs(sum(all_e0(n-10:n-1:1))/10 - e0) .le. reduction_thereshold) then
          damping = damping*0.9_dp
          reduction_thereshold = reduction_thereshold*0.1_dp
          print*, "damping reduced!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
      end if
    end if

    write(unit=12, fmt=*) et, no_walkers, e0
    if (modulo(n,10) .eq. 0) then
      write(*,'(A,I5,A,I5,A)') "Step ", n, " of ", n_steps, "--------------------------------"
      print*, "number of walkers = ", no_walkers
      print*, "et = ", et
      print*, "e0 = ", e0
    end if
  end do
  !end of the simulation--------------------------------------------------------

  zpe = sum(all_e0(n_steps-min(100,n_steps)+1:n_steps))/min(100,n_steps) - etot
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

subroutine gradient_descent(ats)
  use parameters
  use precision
  use atomicStructure
  use RuNNerInterface
  implicit none
  type(atStruct) :: ats
  real(dp) :: etot,alpha=1d-3,norm=1,vectornorm=0._dp,vectornorm2=0._dp,angle,angle1
  real(dp) ,  dimension(3,ats%nat) :: fxyz, prev_fxyz
  integer :: i, j, k, ios, iter=0

  call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, etot, fxyz)
  write(unit=*, fmt=*) "The energy in the given state is: ", etot

  do while ( norm>1d-6 )
    iter = iter+1 !count iterations
    prev_fxyz = fxyz
    ats%ats = ats%ats + fxyz*alpha !updating the positions
    call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, etot, fxyz) !calculates negative gradient and total energy

    !calculating the gradient norm and the angle between the vectors
    angle = 0._dp
    norm = 0._dp
    do i=1,nat
      vectornorm = 0._dp
      vectornorm2 = 0._dp
      do j=1,3
        vectornorm = vectornorm+fxyz(j,i)**2
        vectornorm2 = vectornorm2+prev_fxyz(j,i)**2
      end do
      norm = norm + sqrt(vectornorm)/nat !average over all the gradient norms
      !angle between the current and previous gradient vector of the atom i
      angle=angle+acos(dot_product(prev_fxyz(:,i),fxyz(:,i))/(sqrt(vectornorm)*sqrt(vectornorm2)))
    end do
    if (angle/nat > 1.0472) then  !average of all the angles, gradient feedback
      alpha = alpha/2._dp
    else
      alpha = alpha*1.05_dp
    endif
  end do

  write(unit=*, fmt=*) "The energy after the minimization is: ", etot
  write(unit=*, fmt='(I5,A)') iter, " iterations were needed to obtain a precision of 1d-6"
  write(unit=*, fmt=*)

end subroutine gradient_descent
