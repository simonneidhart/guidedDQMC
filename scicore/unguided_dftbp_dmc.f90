!Unguided DQMC according to McCoy
!Anne B. McCoy. Diffusion monte carlo approaches for investigating the structure and
!vibrational spectra of fluxional systems. International Reviews in Physical Chemistry, 25(1-
!2):77–107, 2006.

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
  real(8), allocatable, dimension(:) :: all_et
  real(8), dimension(3,nat) :: fxyz, rxyz_in
  real(8), dimension(3,3) :: alat, deralat
  real(8), dimension(24) :: r_temp
  logical :: half_walker_flag = .FALSE., restart = .FALSE.

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
  no_walkers_start = 50 !initial number of walkers
  no_walkers_restart = 620
  maxWalkers = 10000
  n_steps = 100 !number of simulation steps
  a = 1
  damping = 0.01d0 !damping constant for the update of et in every iteration
  et = -5.62 !inital guess for the energy
  et_restart = -5.62663416392529

  print*, "damping = ", damping

  call read_parameters()

  !initialzations and memory allocations----------------------------------------
  call setup_dftbp("dftb_temp_folder/")
  allocate(all_et(n_steps))

  !start or restat
  if (restart .eq. .TRUE.) then
    old_no_walkers = no_walkers_restart
    no_walkers = no_walkers_restart
    et = et_restart
    allocate(walkers(no_walkers))
    call read_restart(walkers, no_walkers)
    write(*,*) "Reading restart files complete: ", no_walkers, " walkers, ", et, " et"
  else
    old_no_walkers = no_walkers_start
    no_walkers = no_walkers_start
    allocate(walkers(no_walkers))
    call energyandforces(nat, rxyz0, fxyz, etot, alat, deralat, atomnames)
    call initalizeWalkers(walkers,no_walkers,etot)
    write(*,*) "Input etot: ",etot
  end if

  open(unit=12, file="et_noWalkers", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 12"
  write(unit=12, fmt=*) et, no_walkers, et

  !Start of the simulation------------------------------------------------------
  do n=1,n_steps

    do j=1,size(walkers)
      call propagate(walkers(j), et, alat, deralat)   !move the walker
    end do

    projected_no_walkers = 0
    do j=1,size(walkers)
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
      do j=1,size(walkers)
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
          else !walker survives
            new_walkers(k) = walkers(j)
            k = k + 1
          end if
        end if
      end do
      if (k-1 .ne. new_no_walkers ) then
        print*, "no_walkers wrong", k-1, new_no_walkers
        new_no_walkers = k-1
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
      call initalizeWalkers(walkers,no_walkers,etot)
    end if

    !adjust the trial energy
    v_bar = 0
    sum_weight = 0
    do j=1,size(walkers)
      v_bar = v_bar + walkers(j)%etot * walkers(j)%weight
      sum_weight = sum_weight + walkers(j)%weight
    end do
    v_bar = v_bar/sum_weight
    et = v_bar - damping*(new_no_walkers-no_walkers_start)/no_walkers_start

    all_et(n) = et
    no_walkers = new_no_walkers
    energy_count = energy_count + no_walkers

    write(unit=12, fmt=*) et, no_walkers
    if (modulo(n,10) .eq. 0) then
      write(*,'(A,I5,A,I5,A)') "Step ", n, " of ", n_steps, "--------------------------------"
      print*, "number of walkers = ", no_walkers
      print*, "et = ", et
      !write the walker postions for possible restart
      open(unit=13, file="walker_positions", iostat=ios, action="write")
      if ( ios /= 0 ) stop "Error opening file unit 13"
      do i=1,size(walkers)
        call write_walker(walkers(i), 13)
      end do
      close(unit=13, iostat=ios, status="keep")
      if ( ios /= 0 ) stop "Error closing file unit 13"
    end if
  end do
  !end of the simulation--------------------------------------------------------

  print*, "number of walkers = ", no_walkers
  !print*, "average_et = ", sum(all_et(n-101:n-1))/100
  print*, "et = ", et

  print*, "Number of energy calculations = ", energy_count

  close(unit=12, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 12"

  !output final positions of the walkers
  open(unit=13, file="walker_positions", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,size(walkers)
    call write_walker(walkers(i), 13)
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"

end program dmc

!reading the imput parameters for restarting the algorithm----------------------
subroutine read_restart(walkers, no_walkers)
  use dftbp_walker
  use parameters
  implicit none
  integer :: i, j, k, ios, no_walkers
  type(walker), dimension(no_walkers) :: walkers
  real(8) :: weight
  real(8), dimension(3,nat) :: rxyz

  open(unit=16, file="walker_positions", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 16"

  do k=1,size(walkers)
    read(16,*) weight
    do i=1,nat
      read(16,*) ( rxyz(j, i), j = 1, 3 )
    end do
    call create_walker(walkers(k), rxyz, 1.d0, weight)
  end do

  close(unit=16, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 16"



end subroutine read_restart

!inializes all walkers in the geometry opimized position------------------------
subroutine initalizeWalkers(walkers,no_walkers,etot)
  use dftbp_walker
  use parameters
  implicit real(8) (a-h,o-z)
  type(walker), dimension(no_walkers) :: walkers

  do i=1,no_walkers
    call create_walker(walkers(i), rxyz0, etot, 1.d0)
  end do

end subroutine initalizeWalkers

!creates an array of 3 x nat normal distributed random numbers using
!the Box-Muller algorithm.------------------------------------------------------
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
