program dmc
  use trial_wf_pot
  implicit real(8) (a-h,o-z)
  integer :: energy_count = 0, a !count the number of energy calculations
  real(8), allocatable, dimension(:) :: all_e0, weight, new_weight
  real(8), allocatable, dimension(:,:) :: walkers, new_walkers, normal_rand
  real(8), dimension(dim) :: walker_prev, half_walker
  real(8), dimension(dim,2) :: box
  real(8), dimension(2,dim) :: mu_sigma
  logical, parameter :: adjust_flag = .TRUE. !.TRUE. = with damping adjusting
  logical :: half_walker_flag = .FALSE.

  !Simulation parameters:-------------------------------------------------------
  bl = 2 !side length of the box
  !box = reshape((/-bl,&
  !                 bl/),shape(box))
  box = reshape((/-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,&
                   bl, bl, bl, bl, bl, bl, bl, bl, bl, bl, bl, bl/),shape(box)) !size of the inital simulation box
  no_walkers_start = 1000 !initial number of walkers
  maxWalkers = 1000000
  delta_t = 0.1
  n_steps = 500 !number of simulation steps
  a = 1
  damping = 0.05 !damping constant for the update of et in every iteration
  reduction_thereshold = 0.2 !reduce the damping constant if the error is below this thereshold
  et = 10 !inital guess for the energy

  !initialzations and memory allocations----------------------------------------
  old_no_walkers = no_walkers_start
  no_walkers = no_walkers_start
  allocate(all_e0(n_steps))
  allocate(walkers(no_walkers,dim))
  allocate(weight(no_walkers))
  call initalizeWalkers(walkers,no_walkers,dim,box)
  weight = 1.0
  !walkers = reshape((/0, 0, 0/),shape(walkers))

  open(unit=12, file="et_noWalkers", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 12"
  write(unit=12, fmt=*) et, no_walkers
  !print*, "number of walkers = ", no_walkers
  !print*, "et = ", et

  !Start of the simulation------------------------------------------------------
  do n=1,n_steps

    allocate(normal_rand(no_walkers,dim))
    call normal(normal_rand, no_walkers, dim)
    projected_no_walkers = 1
    do j=1,no_walkers
      walker_prev = walkers(j,:)
      !move the walkers according to the Gaussian distribution------------------
      walkers(j,:) = normal_rand(j,:)*sqrt(delta_t) + walkers(j,:) + drift(walkers(j,:))*delta_t
      !walkers(j,:) = normal_rand(j,:)*sqrt(delta_t) + walkers(j,:)
      !calculate the weights----------------------------------------------------
      weight(j) = weight(j)*exp(-delta_t*((el(walkers(j,:)) + el(walker_prev))/2 - et))
      !weight(j) = exp(-delta_t*(v(walkers(j,:)) - et))
      !calculate the number of walkers after the branching process--------------
      if (weight(j) .gt. 2) then !calculate the new number of walkers (multiplied by two)
        projected_no_walkers = projected_no_walkers + int(weight(j))
      elseif ( weight(j) .lt. 0.5) then
        projected_no_walkers = projected_no_walkers + 0.5
      else
        projected_no_walkers = projected_no_walkers + 1
      end if
    end do
    deallocate(normal_rand)

    new_no_walkers = floor(projected_no_walkers)

    if (new_no_walkers .lt. maxWalkers) then! normal case
      allocate(new_walkers(new_no_walkers,dim))
      allocate(new_weight(new_no_walkers))
      !split-joint algorithm----------------------------------------------------
      k = 1
      do j=1,no_walkers
        if (weight(j) .gt. 2) then
          do i=1,int(weight(j))
            new_walkers(k,:) = walkers(j,:)
            new_weight(k) = weight(j)/int(weight(j))
            k = k + 1
          end do
        elseif (weight(j) .lt. 0.5) then
          if (half_walker_flag) then
            half_walker_flag = .FALSE.
            call RANDOM_NUMBER(r)
            if (r .lt. weight(j)/(weight(j)+half_weight)) then
              new_walkers(k,:) = walkers(j,:)
            else
              new_walkers(k,:) = half_walker
            end if
            new_weight(k) = weight(j) + half_weight
            k = k + 1
          else
            half_walker_flag = .TRUE.
            half_walker = walkers(j,:)
            half_weight = weight(j)
          end if
        else
          new_walkers(k,:) = walkers(j,:)
          new_weight(k) = weight(j)
          k = k + 1
        end if
      end do
      if (k-1 .ne. new_no_walkers ) then
        !print*, "no_walkers wrong", k-1, new_no_walkers
      end if
      deallocate(walkers)
      allocate(walkers(new_no_walkers,dim))
      walkers = new_walkers
      deallocate(new_walkers)
      deallocate(weight)
      allocate(weight(new_no_walkers))
      weight = new_weight
      deallocate(new_weight)

    else !Too many walkers. Choose 1000 walkers randomly to survive-------------
      print*, "Too many walkers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      new_no_walkers = no_walkers_start
      allocate(new_walkers(new_no_walkers,dim))
      allocate(new_weight(new_no_walkers))
      do j=1,new_no_walkers
        call RANDOM_NUMBER(r)
        index = int(r*new_no_walkers+1)
        new_walkers(j,dim) = walkers(index,dim)
        new_weight(j) = weight(index)
      end do
      deallocate(walkers)
      allocate(walkers(new_no_walkers,dim))
      walkers = new_walkers
      deallocate(new_walkers)
      deallocate(weight)
      allocate(weight(new_no_walkers))
      weight = new_weight
      deallocate(new_weight)
      no_walkers = no_walkers_start
      old_no_walkers = no_walkers_start
    end if

    if (new_no_walkers .eq. 0) then !to keep the simulation going if the population dies out
      print*, "Walkers died out!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      no_walkers = no_walkers_start
      new_no_walkers = no_walkers_start
      old_no_walkers = no_walkers_start
      deallocate(walkers)
      allocate(walkers(no_walkers,dim))
      call initalizeWalkers(walkers,no_walkers,dim,box)
      deallocate(weight)
      allocate(weight(no_walkers))
      weight = 1.0
    end if

    !v_bar = 0 !alternative for adjusting et (Mc Coy)
    !do j=1,new_no_walkers
    !  v_bar = v_bar + v(walkers(j,:))
    !end do
    !v_bar = v_bar/new_no_walkers
    !et = v_bar - 1/delta_t*(new_no_walkers-no_walkers_start)/no_walkers_start

    !adjust et with formula form the master exam paper (A=1)
    e0 = 0
    do j=1,new_no_walkers
      e0 = e0 + el(walkers(j,:))*weight(j)
    end do
    e0 = e0/sum(weight) !best energy guess
    et = e0 - damping*log(sum(weight)/real(no_walkers_start))

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
          damping = damping*0.9
          reduction_thereshold = reduction_thereshold*0.1
          print*, "damping reduced!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
      end if
    end if

    write(unit=12, fmt=*) e0, no_walkers
    if (modulo(n,1) .eq. 0) then
      write(*,'(A,I5,A,I5,A)') "Step ", n, " of ", n_steps, "--------------------------------"
      print*, "number of walkers = ", no_walkers
      print*, "et = ", et
      print*, "e0 = ", e0
    end if
  end do
  !end of the simulation--------------------------------------------------------

  print*, "number of walkers = ", no_walkers
  print*, "average_e0 = ", sum(all_e0(n-101:n-1))/100
  print*, "e0 = ", e0

  call calc_e0_theory(e0_theory)
  print*, "e0_theory = ", e0_theory
  print*, "Number of energy calculations = ", energy_count

  close(unit=12, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 12"

  !output final positions of the walkers
  open(unit=13, file="walker_positions", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,no_walkers
    write(unit=13, fmt=*) (walkers(i,k), k=1,dim), weight(i)
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"

  !output mu and sigma for the final walker distribution
  call calc_mu_sigma(mu_sigma,walkers,no_walkers,dim)
  open(unit=14, file="mu_sigma", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 14"
  do i=1,dim
    write(unit=14, fmt=*) mu_sigma(1,i), mu_sigma(2,i)
  end do
  close(unit=14, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 14"

end program dmc

!uniformly initalizes walker positions in the given box.------------------------
subroutine initalizeWalkers(walkers,no_walkers,dim,box)
  implicit real(8) (a-h,o-z)
  integer :: dim
  real(8), dimension(no_walkers,dim) :: walkers
  real(8), dimension(dim,2) :: box

  call RANDOM_NUMBER(walkers)
  do i=1,no_walkers
    do j=1,dim
      ! transformation from [0,1] to [a,b]: a + x*(b - a)
      walkers(i,j) = box(j,1) + walkers(i,j)*(box(j,2) - box(j,1))
    end do
  end do

end subroutine initalizeWalkers

!creates an array of no_walkers x dim normal distributed random numbers using----
!the Box-Muller algorithm.
subroutine normal(normal_rand, no_walkers, dim)
  implicit real(8) (a-h,o-z)
  integer :: dim
  real(8), allocatable, dimension(:) :: rand
  real(8), dimension(no_walkers,dim) :: normal_rand
  pi = 4.d0*atan(1.d0)

  n = no_walkers*dim
  if (mod(n,2) .eq. 0) then !even number of random numbers wanted
    allocate(rand(n))
    call RANDOM_NUMBER(rand)
    do i=1,n/2
      x = rand(2*i-1)
      y = rand(2*i)
      rand(2*i-1) = sqrt(-2.d0*log(x))*cos(2.d0*pi*y)
      rand(2*i) = sqrt(-2.d0*log(x))*sin(2.d0*pi*y)
    end do
  else !odd number of random numbers wanted
    allocate(rand(n+1))
    call RANDOM_NUMBER(rand)
    do i=1,n/2
      x = rand(2*i-1)
      y = rand(2*i)
      rand(2*i-1) = sqrt(-2.d0*log(x))*cos(2.d0*pi*y)
      rand(2*i) = sqrt(-2.d0*log(x))*sin(2.d0*pi*y)
    end do
      rand(n) = sqrt(-2.0*log(rand(n)))*cos(2*pi*rand(n+1))
  end if
  normal_rand = reshape(rand,shape(normal_rand))
  deallocate(rand)

end subroutine normal

!calculate mu and sigma for the walkers along all dimensions--------------------
subroutine calc_mu_sigma(mu_sigma,walkers,no_walkers,dim)
  implicit real(8) (a-h,o-z)
  integer :: dim
  real(8), dimension(no_walkers,dim) :: walkers
  real(8), dimension(2,dim) :: mu_sigma

  do i=1,dim
    mu_sigma(1,i) = sum(walkers(:,i))/no_walkers
    mu_sigma(2,i) = 0
    do j=1,no_walkers
      mu_sigma(2,i) = mu_sigma(2,i) + (walkers(j,i) - mu_sigma(1,i))**2
    end do
    mu_sigma(2,i) = sqrt(mu_sigma(2,i)/no_walkers)
  end do

end subroutine calc_mu_sigma
