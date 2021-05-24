module dftbp_walker
  use parameters
  type walker
    real(8), dimension(dim) :: r
    real(8) :: etot, weight
  contains
    procedure :: create_walker
    procedure :: propagate
    procedure :: write_walker

  end type walker

contains

  subroutine create_walker(this, r_, etot_, weight_)
    implicit none
    class(walker) :: this
    real(8), dimension(dim) :: r_
    real(8) :: etot_, weight_

    this%r = r_
    this%etot = etot_ + el_kin(mu)
    this%weight = weight_


  end subroutine create_walker

  subroutine propagate(this, et, delta_t)
    use mod_dftbp
    use parameters
    implicit none
    integer :: i
    class(walker) :: this
    real(8), dimension(3,nat) :: rxyz, fxyz
    real(8), dimension(dim) :: normal_rand
    real(8), dimension(fulldim) :: r_temp
    real(8) :: delta_t, etot_prev, et, epot, ndelta_t = 1.d-3

    call normal(normal_rand,dim)

    etot_prev = this%etot
    !move the walkers according to the Gaussian distribution------------------
    !print*, shape(normal_rand), shape(this%rxyz), shape(this%fxyz)
    this%r = this%r + normal_rand*sqrt(delta_t/new_masses) + drift(this%r)*delta_t/new_masses
    !this%r = this%r + normal_rand*sqrt(ndelta_t) + drift(this%r)*ndelta_t
    !this%r = this%r + normal_rand*sqrt(delta_t/average_mass)
    do i=1,dim
      if (abs(this%r(i) - mu(i)) .gt. 2.0) then
        this%r(i) = mu(i)
        print*, "walker moved to much!"
      end if
    end do
    !calculate the local energy------------------------------------------------
    r_temp = matmul(ev,this%r)
    rxyz = reshape(r_temp,shape(rxyz)) + rxyz0
    call energyandforces(nat, rxyz, fxyz, epot, alat, deralat, atomnames)
    this%etot = epot + el_kin(this%r)
    !this%etot = epot
    !calculate the weights----------------------------------------------------
    this%weight = this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))
    !this%weight = this%weight*exp(-ndelta_t*((this%etot + etot_prev)/2 - et))
    !this%weight = this%weight*exp(-delta_t*(this%etot - et))
    !print*, "walker epot, ekin, weight: ", epot, el_kin(this%r), this%weight
    !print*, "Weight, etot ", this%weight, this%etot


  end subroutine propagate

  !calculate the drift velocity grad(psi)/psi-------------------------------------
  function drift(r) result(v)
    implicit real(8) (a-h,o-z)
    real(8), dimension(dim), intent(in) :: r
    real(8), dimension(dim) :: v

    do i=1,dim
      v(i) = -(r(i)-mu(i))/sigma2(i)
    end do

  end function drift

  !calculate the kinetic part of the local energy psi*H/psi-------------------------------------------
  real(8) function el_kin(r)
    implicit none
    real(8), dimension(dim), intent(in) :: r
    real(8), dimension(fulldim) :: r_temp, r_plus, r_minus
    real(8) :: h
    integer :: i

    el_kin = 0
    do i=1,dim
      el_kin = el_kin - 0.5/sigma2(i)*((r(i)-mu(i))**2/sigma2(i) - 1)/new_masses(i)
      !el_kin = el_kin - 0.5/sigma2(i)*((r(i)-mu(i))**2/sigma2(i) - 1)
    end do

  end function el_kin

  real(8) function psi(r)
    implicit none
    real(8), dimension(dim), intent(in) :: r
    integer :: i

    psi = 0
    do i=1,dim
      psi = psi -(r(i)-mu(i))**2/(2*sigma2(i))
    end do
    psi = exp(psi)

  end function psi

  subroutine write_walker(this, iounit)
    implicit none
    class(walker) :: this
    integer :: iounit, i
    real(8), dimension(fulldim) :: r_temp
    real(8), dimension(3,nat) :: rxyz

    r_temp = 0
    do i=1,dim
      r_temp = r_temp + ev(:,i)*this%r(i)
    end do
    rxyz = reshape(r_temp,shape(rxyz))

    !write(unit=iounit, fmt=*) "weight: ", this%weight, "etot: ", this%etot
    !do i=1,nat
    !  write(unit=iounit, fmt=*) rxyz(:,i)
    !end do

    do i=1,dim
      write(unit=iounit, fmt=*) this%r(i)
    end do

  end subroutine write_walker

  !uniformly initalizes walker positions in the given box.------------------------
  subroutine initalizeWalkers(walkers,no_walkers,etot)
    implicit real(8) (a-h,o-z)
    type(walker), dimension(no_walkers) :: walkers

    do i=1,no_walkers
      call create_walker(walkers(i), mu, etot, 1.d0)
    end do

  end subroutine initalizeWalkers

end module dftbp_walker

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
