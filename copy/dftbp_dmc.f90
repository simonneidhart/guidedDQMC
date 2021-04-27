program dmc
  implicit real(8) (a-h,o-z)
  character(len=2), dimension(nat) :: atomnames = (/"C", "C", "H", "H", "H", "H", "H", "H"/)

  et_start = -5.62 !inital guess for the energy
  damping = 0.01d0 !damping constant for the update of et in every iteration
  no_walkers_start = 50 !initial number of walkers
  maxWalkers = 10000
  delta_t = 1.0d0 !does not work with 0.01 and above
  n_steps = 200 !number of simulation steps



end program dmc

subroutine run_simulation(n_steps, no_walkers_start, delta_t, et_start, maxWalkers, damping, atomnames)
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
  real (8) :: Bohr_Ang
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
  real(8), allocatable, dimension(:) :: all_e0
  real(8), dimension(3,nat) :: fxyz, rxyz
  real(8), dimension(3,3) :: alat, deralat
  real(8), dimension(24) :: r_temp
  character(len=2), dimension(nat) :: atomnames
  logical, parameter :: adjust_flag = .FALSE. !.TRUE. = with damping adjusting
  logical :: half_walker_flag = .FALSE.

  Bohr_Ang = 0.529d0

  alat(1,1)=1.d0
  alat(2,1)=0.d0
  alat(3,1)=0.d0

  alat(1,2)=0.d0
  alat(2,2)=1.d0
  alat(3,2)=0.d0

  alat(1,3)=0.d0
  alat(2,3)=0.d0
  alat(3,3)=1.d0

  alat=alat/Bohr_Ang

  !Simulation parameters:-------------------------------------------------------
  a = 1
  et = et_start
  reduction_thereshold = 0.2 !reduce the damping constant if the error is below this thereshold

  print*, "delta_t = ", delta_t

  call read_parameters()

  !initialzations and memory allocations----------------------------------------
  old_no_walkers = no_walkers_start
  no_walkers = no_walkers_start
  allocate(all_e0(n_steps))
  allocate(walkers(no_walkers))
  call setup_dftbp("dftb_temp_folder/")
  r_temp = 0
  do i=1,dim
    r_temp = r_temp + ev(:,i)*mu(i)
  end do
  rxyz = reshape(r_temp,shape(rxyz))
  call energyandforces(nat, rxyz, fxyz, etot, alat, deralat, atomnames)
  call initalizeWalkers(walkers,no_walkers,etot)
  write(*,*) "Input etot: ",etot

  open(unit=12, file="et_noWalkers", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 12"
  e0 = et
  write(unit=12, fmt=*) et, no_walkers, e0

  !Start of the simulation------------------------------------------------------
  do n=1,n_steps

    projected_no_walkers = 0
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
      call initalizeWalkers(walkers,no_walkers,etot)
    end if

    !v_bar = 0 !alternative for adjusting et (Mc Coy)
    !do j=1,new_no_walkers
    !  v_bar = v_bar + walkers(j)%etot
    !end do
    !v_bar = v_bar/new_no_walkers
    !et = v_bar - damping*(new_no_walkers-no_walkers_start)/no_walkers_start

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

    write(unit=12, fmt=*) et, no_walkers, e0
    if (modulo(n,1) .eq. 0) then
      write(*,'(A,I5,A,I5,A)') "Step ", n, " of ", n_steps, "--------------------------------"
      print*, "number of walkers = ", no_walkers
      print*, "et = ", et
      print*, "e0 = ", e0
    end if
  end do
  !end of the simulation--------------------------------------------------------

  print*, "number of walkers = ", no_walkers
  !print*, "average_e0 = ", sum(all_e0(n-101:n-1))/100
  print*, "e0 = ", e0

  print*, "Number of energy calculations = ", energy_count

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

  open(unit=14, file="walker_energies", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 14"
  do i=1,no_walkers
    write(unit=14, fmt=*) walkers(i)%etot
  end do
  close(unit=14, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 14"
end subroutine run_simulation

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
