module dftbp_walker
  use precision
  use atomicStructure
  use RuNNerInterface
  use parameters
  type walker
    real(dp), allocatable, dimension(:) :: r
    real(dp) :: etot, weight
  contains
    procedure :: create_walker
    procedure :: propagate
    procedure :: write_walker_kartesian

  end type walker

contains

  subroutine create_walker(this, etot_, weight_)
    use precision
    implicit none
    class(walker) :: this
    real(dp) :: etot_, weight_

    allocate(this%r(dim))
    this%r = 0.d0
    this%etot = etot_ + el_kin(this%r)
    this%weight = weight_


  end subroutine create_walker

  subroutine propagate(this, et, delta_t, thread_num, ats)
    use precision
    use atomicStructure
    use RuNNerInterface
    use parameters
    implicit none
    integer :: i
    class(walker) :: this
    type(atStruct) :: ats
    real(dp), dimension(3,nat) :: fxyz
    real(dp), dimension(dim) :: normal_rand
    real(dp), dimension(fulldim) :: r_temp
    real(dp) :: delta_t, etot_prev, et, epot
    integer :: thread_num

    call normal(normal_rand,dim)

    etot_prev = this%etot
    !move the walkers according to the Gaussian distribution--------------------
    this%r = this%r + normal_rand*sqrt(delta_t/new_masses) + drift(this%r)*delta_t/new_masses
    !this%r = this%r + normal_rand*sqrt(ndelta_t) + drift(this%r)*ndelta_t
    !this%r = this%r + normal_rand*sqrt(delta_t/average_mass)
    !calculate the local energy-------------------------------------------------
    r_temp = matmul(ev,this%r)
    ats%ats = reshape(r_temp,shape(fxyz)) + rxyz0

    call getRuNNerEnergiesAndForces(globalRuNNerHandles(thread_num), ats, epot, fxyz)
    this%etot = epot + el_kin(this%r)
    !calculate the weights------------------------------------------------------
    if (-delta_t*((this%etot + etot_prev)/2._dp - et) .ge. 10) then
      print*, "weight update error", (-delta_t*((this%etot + etot_prev)/2._dp - et)), this%etot, etot_prev
    else
      this%weight = this%weight*exp(-delta_t*((this%etot + etot_prev)/2._dp - et))
    end if

  end subroutine propagate

  !calculate the drift velocity grad(psi)/psi-----------------------------------
  function drift(r) result(v)
    use precision
    implicit real(dp) (a-h,o-z)
    real(dp), dimension(dim), intent(in) :: r
    real(dp), dimension(dim) :: v

    do i=1,dim
      if (guided_dims(i)) then
        v(i) = -r(i)/sigma2(i)
      else
        v(i) = 0._dp
      end if
    end do

  end function drift

  !calculate the kinetic part of the local energy psi*H/psi-------------------------------------------
  real(dp) function el_kin(r)
    use precision
    implicit none
    real(dp), dimension(dim), intent(in) :: r
    real(dp), dimension(fulldim) :: r_temp, r_plus, r_minus
    real(dp) :: h
    integer :: i

    el_kin = 0
    do i=1,dim
      if (guided_dims(i)) then
        el_kin = el_kin - 0.5_dp/sigma2(i)*(r(i)**2/sigma2(i) - 1._dp)/new_masses(i)
      end if
    end do

  end function el_kin

  real(dp) function psi(r)
    use precision
    implicit none
    real(dp), dimension(dim), intent(in) :: r
    integer :: i

    psi = 0._dp
    do i=1,dim
      if (guided_dims(i)) then
        psi = psi -r(i)**2/(2._dp*sigma2(i))
      end if
    end do
    psi = exp(psi)

  end function psi

  subroutine write_walker_kartesian(this, iounit)
    use precision
    implicit none
    class(walker) :: this
    integer :: iounit, i
    real(dp), dimension(fulldim) :: r_temp
    real(dp), dimension(3,nat) :: rxyz

    r_temp = matmul(ev,this%r)
    rxyz = reshape(r_temp,shape(rxyz)) + rxyz0

    do i=1,nat
      write(unit=iounit, fmt=*) rxyz(1,i), rxyz(2,i), rxyz(3,i)
    end do

  end subroutine write_walker_kartesian

  subroutine write_walker_hessian(this, iounit)
    implicit none
    class(walker) :: this
    integer :: iounit, i

    do i=1,dim
      write(unit=iounit, fmt=*) this%r(i)
    end do

  end subroutine write_walker_hessian

  !uniformly initalizes walker positions in the given box.------------------------
  subroutine initalizeWalkers(walkers,no_walkers,etot)
    use precision
    implicit real(dp) (a-h,o-z)
    type(walker), dimension(no_walkers) :: walkers

    do i=1,no_walkers
      call create_walker(walkers(i), etot, 1._dp)
    end do

  end subroutine initalizeWalkers

end module dftbp_walker

!creates an array of n normal distributed random numbers using----
!the Box-Muller algorithm.
subroutine normal(rand, n)
  use precision
  implicit real(dp) (a-h,o-z)
  real(dp), dimension(n) :: rand
  pi = 4._dp*atan(1._dp)

  if (mod(n,2) .eq. 0) then !even number of random numbers wanted
    call RANDOM_NUMBER(rand)
    do i=1,n/2
      x = rand(2*i-1)
      y = rand(2*i)
      rand(2*i-1) = sqrt(-2._dp*log(x))*cos(2._dp*pi*y)
      rand(2*i) = sqrt(-2._dp*log(x))*sin(2._dp*pi*y)
    end do
  else !odd number of random numbers wanted
    call RANDOM_NUMBER(rand)
    do i=1,n/2
      x = rand(2*i-1)
      y = rand(2*i)
      rand(2*i-1) = sqrt(-2._dp*log(x))*cos(2._dp*pi*y)
      rand(2*i) = sqrt(-2._dp*log(x))*sin(2._dp*pi*y)
    end do
    call RANDOM_NUMBER(r2) !extra random number
    rand(n) = sqrt(-2._dp*log(rand(n)))*cos(2._dp*pi*r2)
  end if

end subroutine normal
