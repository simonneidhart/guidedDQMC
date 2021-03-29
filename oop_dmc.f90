program dmc
  use trial_wf_pot
  use dftbp_walker
  implicit real(8) (a-h,o-z)
  interface
    subroutine initalizeWalkers(walkers,no_walkers,nat,atomnames)
      use dftbp_walker
      implicit real(8) (a-h,o-z)
      real(8) :: rxyz(3,nat)
      type(walker), dimension(no_walkers) :: walkers
      character(len=2) :: atomnames(nat)
      character(len=100) :: filename, names, units, boundary_cond
    end subroutine initalizeWalkers
  end interface
  integer, parameter :: nat = 4 !number of atoms
  integer :: energy_count = 0, a !count the number of energy calculations
  type(walker), allocatable, dimension(:) :: walkers, new_walkers
  type(walker) :: half_walker
  real(8), allocatable, dimension(:) :: all_e0
  real(8), dimension(3,nat) :: fxyz
  real(8), dimension(3,3) :: alat, deralat
  character(len=2) :: atomnames(nat)
  logical, parameter :: adjust_flag = .FALSE. !.TRUE. = with damping adjusting
  logical :: half_walker_flag = .FALSE.

  !Simulation parameters:-------------------------------------------------------
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
  allocate(all_e0(n_steps))
  allocate(walkers(no_walkers))
  call initalizeWalkers(walkers,no_walkers,nat,atomnames)
  call energyandforces_dummy(nat, walkers(1)%rxyz, fxyz, etot, alat, deralat, atomnames)
  do j=1,no_walkers
    walkers(j)%etot = etot
    walkers(j)%fxyz = fxyz
  end do

  open(unit=12, file="et_noWalkers", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 12"
  write(unit=12, fmt=*) et, no_walkers

  !Start of the simulation------------------------------------------------------
  do n=1,n_steps

    projected_no_walkers = 1
    do j=1,no_walkers
      !move the walker
      call propagate(walkers(j), et, delta_t, alat, deralat, atomnames)
      !calculate the number of walkers after the branching process--------------
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
        if (walkers(j)%weight .gt. 2) then !split
          do i=1,int(walkers(j)%weight)
            new_walkers(k) = walkers(j)
            new_walkers(k)%weight = walkers(j)%weight/int(walkers(j)%weight)
            k = k + 1
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
      end do
      if (k-1 .ne. new_no_walkers ) then
        !print*, "no_walkers wrong", k-1, new_no_walkers
      end if
      deallocate(walkers)
      allocate(walkers(new_no_walkers))
      walkers = new_walkers
      deallocate(new_walkers)

    else !Too many walkers. Choose 1000 walkers randomly to survive-------------
      print*, "Too many walkers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
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
      print*, "Walkers died out!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      no_walkers = no_walkers_start
      new_no_walkers = no_walkers_start
      old_no_walkers = no_walkers_start
      deallocate(walkers)
      allocate(walkers(no_walkers))
      call initalizeWalkers(walkers,no_walkers,nat,atomnames)
      call energyandforces_dummy(nat, walkers(1)%rxyz, fxyz, etot, alat, deralat, atomnames)
      do j=1,no_walkers
        walkers(j)%etot = etot
        walkers(j)%fxyz = fxyz
      end do
    end if

    !v_bar = 0 !alternative for adjusting et (Mc Coy)
    !do j=1,new_no_walkers
    !  v_bar = v_bar + v(walkers(j,:))
    !end do
    !v_bar = v_bar/new_no_walkers
    !et = v_bar - 1/delta_t*(new_no_walkers-no_walkers_start)/no_walkers_start

    !adjust et with formula form the master exam paper (A=1)
    e0 = 0
    sum_weight = 0
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
          damping = damping*0.9
          reduction_thereshold = reduction_thereshold*0.1
          print*, "damping reduced!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
      end if
    end if

    write(unit=12, fmt=*) e0, no_walkers
    if (modulo(n,100) .eq. 0) then
      write(*,'(A,I5,A,I5,A)') "Step ", n, " of ", n_steps, "--------------------------------"
      print*, "number of walkers = ", no_walkers
      print*, "et = ", et
      print*, "e0 = ", e0
    end if
  end do
  !end of the simulation--------------------------------------------------------

  print*, "number of walkers = ", no_walkers
  print*, "average_e0 = ", sum(all_e0(n-8000:n-1))/8000
  print*, "e0 = ", e0

  print*, "Number of energy calculations = ", energy_count

  ener = 0
  call calc_e0_theory(ener)
  print*, "theoretical value", ener

  close(unit=12, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 12"

  !output final positions of the walkers
  open(unit=13, file="walker_positions", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,no_walkers
    call write_walker(walkers(i), 13)
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"

end program dmc

subroutine energyandforces_dummy(nat, rxyz, fxyz, etot, alat, deralat, atomnames)
  use trial_wf_pot
  implicit real(8) (a-h,o-z)
  !! Number of atoms.
  integer :: nat
  !! Positions of the atoms.
  real*8, dimension(3,nat), intent(in) :: rxyz
  !! Forces acting on the atoms.
  real*8, dimension(3,nat) :: fxyz
  !! Potential energy of the system.
  real*8 :: etot
  !! Lattice vectors of the periodic cell.
  real*8, dimension(3,3) :: alat
  !! Negative derivative of the energy with respect to the lattice vectors.
  real*8, dimension(3,3) :: deralat
  !! array containg the chemical symbols of each atom
  character(len=2) :: atomnames(nat)

  etot = el(rxyz)
  deralat = 0
  fxyz = reshape(drift(rxyz),shape(fxyz))

end subroutine energyandforces_dummy

!uniformly initalizes walker positions in the given box.------------------------
subroutine initalizeWalkers(walkers,no_walkers,nat,atomnames)
  use dftbp_walker
  implicit real(8) (a-h,o-z)
  real(8), dimension(3,nat) :: rxyz, fxyz
  type(walker), dimension(no_walkers) :: walkers
  character(len=2) :: atomnames(nat)
  character(len=100) :: filename, names, units, boundary_cond
  fxyz = 0.d0

  filename='C2H6.ascii'
  open(unit=11, file=filename, iostat=ios, status="old")
  if ( ios /= 0 ) stop "Error opening file unit 11"

  read(11,*) natp, names, units

   if (trim(units).ne.'atomic' .and. trim(units).ne.'bohr')  then
     units='ang'
     write(*,*) "Input in Angstroem"
   end if

  if (natp.ne.nat) write(*,*)"Number of atoms does not equal the defined parameter!!!!"

  read(11,*) boundary_cond
  if (trim(boundary_cond) .eq. "free") write(*,*) "Free boundary conditions"


  do i = 1, nat
     read(11, *) atomnames(i), ( rxyz(j, i), j = 1, 3 )
     rxyz(:,i) = (/0.d0, 0.d0, 0.d0/)
  end do

  do i=1,no_walkers
    call create_walker(walkers(i), rxyz, fxyz, 0.d0, 1.d0)
    !walkers(i) = walker(nat, rxyz, fxyz, 0.0, 1.0) !default constructor for class walker
  end do

end subroutine initalizeWalkers

!creates an array of 3 x nat normal distributed random numbers using----
!the Box-Muller algorithm.
subroutine normal(normal_rand, nat)
  implicit real(8) (a-h,o-z)
  real(8), allocatable, dimension(:) :: rand
  real(8), dimension(3,nat) :: normal_rand
  pi = 4.d0*atan(1.d0)

  n = 3*nat
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
