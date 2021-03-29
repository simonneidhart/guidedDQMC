program dmc

  implicit real(8) (a-h,o-z)
  integer, parameter :: dim = 12 !number of dimensions
  integer :: energy_count = 0, a !count the number of energy calculations
  integer, allocatable, dimension(:) :: branching
  real(8), allocatable, dimension(:) :: all_et
  real(8), allocatable, dimension(:,:) :: walkers, new_walkers, normal_rand
  real(8), dimension(dim) :: k_pot
  real(8), dimension(dim,2) :: box
  real(8), dimension(2,dim) :: mu_sigma
  logical, parameter :: adjust_flag = .FALSE. !.TRUE. = with damping adjusting

  !Simulation parameters:-------------------------------------------------------
  k_pot = (/2.0, 1.1, 0.4, 1.8, 3.2, 1.0, 3.4, 1.4, 1.1, 0.6, 1.9, 0.6/) !parameters for the harmonic potential
  bl = 1 !side length of the box
  box = reshape((/-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,-bl,&
                   bl, bl, bl, bl, bl, bl, bl, bl, bl, bl, bl, bl/),shape(box)) !size of the inital simulation box
  no_walkers_start = 10000 !initial number of walkers
  maxWalkers = 1000000
  delta_t = 0.1
  n_steps = 10000 !number of simulation steps
  a = 1
  damping = 0.1 !damping constant for the update of et in every iteration
  reduction_thereshold = 0.2 !reduce the damping constant if the error is below this thereshold
  et = 8 !inital guess for the energy

  !initialzations and memory allocations----------------------------------------
  old_no_walkers = no_walkers_start
  no_walkers = no_walkers_start
  allocate(all_et(n_steps))
  allocate(walkers(no_walkers,dim))
  call initalizeWalkers(walkers,no_walkers,dim,box)

  open(unit=12, file="et_noWalkers", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file "
  write(unit=12, fmt=*) et, no_walkers
  print*, "number of walkers = ", no_walkers
  print*, "et = ", et

  !Start of the simulation------------------------------------------------------
  do n=1,n_steps

    !step 1: move the walkers according to the Gaussian distribution------------
    allocate(normal_rand(no_walkers,dim))
    call normal(normal_rand, no_walkers, dim)
    walkers = normal_rand*sqrt(delta_t) + walkers
    deallocate(normal_rand)

    !step2: branching process---------------------------------------------------
    allocate(branching(no_walkers)) !defines how may new walkers the walker generates
    do j=1,no_walkers
      q = exp(-delta_t*(v(walkers(j,:),k_pot,dim) - et))
      branching(j) = int(q)
      call RANDOM_NUMBER(r)
      if (q - int(q) .gt. r) then
        branching(j) = branching(j) + 1 !walker survives and gives birth to children
      end if
    end do
    energy_count = energy_count + no_walkers
    new_no_walkers = sum(branching)

    if (new_no_walkers .lt. maxWalkers) then !normal case
      allocate(new_walkers(new_no_walkers,dim))
      k = 1 !index for the new_walkers array
      do j=1,no_walkers
        !duplicate the walkers that have children, copy the walkers that survive, skip the walkers that die.
        do while (branching(j) .gt. 0)
          new_walkers(k,:) = walkers(j,:)
          k = k + 1
          branching(j) = branching(j) - 1
        end do
      end do
      deallocate(walkers)
      allocate(walkers(new_no_walkers,dim))
      walkers = new_walkers
      deallocate(new_walkers)
    else !Too many walkers. Choose 1000 walkers randomly to survive.
      print*, "Too many walkers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      new_no_walkers = no_walkers_start
      allocate(new_walkers(new_no_walkers,dim))
      do j=1,new_no_walkers
        call RANDOM_NUMBER(r)
        index = int(r*new_no_walkers+1)
        new_walkers(j,dim) = walkers(index,dim)
      end do
      deallocate(walkers)
      allocate(walkers(new_no_walkers,dim))
      walkers = new_walkers
      deallocate(new_walkers)
      no_walkers = no_walkers_start
    end if
    deallocate(branching)

    if (new_no_walkers .eq. 0) then !to keep the simulation going if the population dies out
      print*, "Walkers died out!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      no_walkers = no_walkers_start
      new_no_walkers = no_walkers_start
      deallocate(walkers)
      allocate(walkers(no_walkers,dim))
      call initalizeWalkers(walkers,no_walkers,dim,box)
    end if

    !v_bar = 0 !alternative for adjusting et (Mc Coy)
    !do j=1,new_no_walkers
    !  v_bar = v_bar + v(walkers(j,:),k_pot,dim)
    !end do
    !v_bar = v_bar/new_no_walkers
    !et = v_bar - 1/delta_t*(new_no_walkers-no_walkers_start)/no_walkers_start

    !adjust et with formula form the master exam paper (A=1)
    !a = int(n/100) + 1
    all_et(n) = et
    if (modulo(n,a) .eq. 0) then
      et = all_et(n-a+1) - damping/(delta_t*a)*log(real(new_no_walkers)/old_no_walkers)
      old_no_walkers = new_no_walkers
    end if

    no_walkers = new_no_walkers

    !adjust the damping term (set flag to true)

    if (adjust_flag) then
      if (n .gt. 10) then
        if (abs(sum(all_et(n-10:n-1:1))/10 - et) .le. reduction_thereshold) then
          damping = damping*0.9
          reduction_thereshold = reduction_thereshold*0.1
          print*, "damping reduced!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
      end if
    end if

    write(unit=12, fmt=*) et, no_walkers
    if (modulo(n,100) .eq. 0) then
      write(*,'(A,I5,A,I5,A)') "Step ", n, " of ", n_steps, "--------------------------------"
      print*, "number of walkers = ", no_walkers
      print*, "et = ", et
    end if
  end do

  print*, "average et = ", sum(all_et(n-8000:n))/8000


  et_theory = 0 !calculate theoretical value for et
  do i=1,dim
    et_theory = et_theory + 0.5*sqrt(k_pot(i))
  end do
  print*, "et_theory = ", et_theory
  print*, "Number of energy calculations = ", energy_count

  !output final positions of the walkers
  open(unit=13, file="walker_positions", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file "
  do i=1,no_walkers
    write(unit=13, fmt=*) (walkers(i,k), k=1,dim)
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"
  close(unit=12, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 12"
  !output mu and sigma for the final walker distribution
  call calc_mu_sigma(mu_sigma,walkers,no_walkers,dim)
  open(unit=14, file="mu_sigma", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file "
  do i=1,dim
    write(unit=14, fmt=*) mu_sigma(1,i), mu_sigma(2,i)
  end do
  close(unit=14, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 14"

end program dmc

!dim-dimensional harmonic potential function------------------------------------
!v(r_1,...,r_dim) = a_1*r_1**2 + ... + a_dim*r_dim**2
real(8) function v(rxyz,k_pot,dim)
  implicit none
  integer :: dim, i
  real(8), dimension(dim), intent(in) :: rxyz, k_pot

  v = 0
  do i=1,size(k_pot)
    v = v + 0.5 * k_pot(i) * rxyz(i)**2
  end do

end function v

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
