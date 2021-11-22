module dftbp_walker
  use parameters
  type walker
    real(8), allocatable, dimension(:) :: r !position coordinates
    real(8) :: etot, weight
  end type walker

contains

  subroutine create_walker(this, etot_, weight_)
    implicit none
    class(walker) :: this
    real(8) :: etot_, weight_
    real(8), dimension(dim) :: normal_rand

    call normal(normal_rand,dim)

    allocate(this%r(dim))
    this%r = normal_rand*sqrt(sigma2)
    this%etot = etot_ + el_kin(this%r)
    this%weight = weight_

  end subroutine create_walker

  subroutine create_walker_restart(this, r_, etot_, weight_)
    implicit none
    class(walker) :: this
    real(8) :: etot_, weight_
    real(8), dimension(dim) :: r_

    allocate(this%r(dim))
    this%r = r_
    this%etot = etot_
    this%weight = weight_

  end subroutine create_walker_restart

  subroutine propagate(this, et, delta_t, thread_num, n)
    use mod_dftbp
    use parameters
    implicit none
    integer :: i, thread_num, n
    class(walker) :: this
    real(8), dimension(3,nat) :: rxyz, fxyz
    real(8), dimension(dim) :: normal_rand, r_prop
    real(8), dimension(fulldim) :: r_temp
    real(8) :: delta_t, etot_prev, et, epot, p_acc, rand

    call normal(normal_rand,dim)

    etot_prev = this%etot
    !move the walkers according to the Gaussian distribution and the drift velocity
    r_prop = this%r + normal_rand*sqrt(delta_t) + drift(this%r,delta_t)*delta_t


    if (n .gt. 20) then
      p_acc = min(1.d0,exp(-dot_product(this%r-r_prop-drift(r_prop,delta_t)*&
      delta_t,this%r-r_prop-drift(r_prop,delta_t)*delta_t))/&
      exp(-dot_product(r_prop-this%r-drift(this%r,delta_t)*delta_t,&
      r_prop-this%r-drift(this%r,delta_t)*delta_t))*psi(r_prop)**2/psi(this%r)**2)
      !p_acc = min(1.d0, exp(dot_product(this%r-r_prop-0.5d0*drift(r_prop,delta_t)*delta_t,drift(r_prop,delta_t)) - &
      !dot_product(r_prop-this%r-0.5d0*drift(this%r,delta_t)*delta_t,drift(this%r,delta_t))+ &
      !dot_product(r_prop/sigma2,r_prop/sigma2)-dot_product(this%r/sigma2,this%r/sigma2)))
    else
      p_acc = 1.d0
    end if

    call RANDOM_NUMBER(rand)

    if (rand .le. p_acc) then !accept
      !print*, "accept", p_acc
      this%r = r_prop

      !calculate the local energy-------------------------------------------------
      r_temp = matmul(ev,this%r)
      rxyz = reshape(r_temp,shape(rxyz)) + rxyz0 !convert to kartesian corrdinated for energy calculation

      call energyandforces(nat, rxyz, fxyz, epot, alat, deralat, atomnames, thread_num)
      this%etot = epot + el_kin(this%r)

      !calculate the weights------------------------------------------------------
      if (-delta_t*((this%etot + etot_prev)/2 - et) .ge. 10) then
        print*, "weight over/underflow: ", (-delta_t*((this%etot + etot_prev)/2 - et)), this%etot, thread_num
        this%etot = etot_prev
      else
        this%weight = this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))
      end if
    end if

  end subroutine propagate

  subroutine propagate_unguided(this, et, delta_t, thread_num)
    use mod_dftbp
    use parameters
    implicit none
    integer :: i, thread_num
    class(walker) :: this
    real(8), dimension(3,nat) :: rxyz, fxyz
    real(8), dimension(dim) :: normal_rand
    real(8), dimension(fulldim) :: r_temp
    real(8) :: delta_t, etot_prev, et, epot

    call normal(normal_rand,dim)

    etot_prev = this%etot
    !move the walkers according to the Gaussian distribution and the drift velocity
    !this%r = this%r + normal_rand*sqrt(delta_t/new_masses) + 0.01*drift(this%r)*delta_t/new_masses
    this%r = this%r + normal_rand*sqrt(delta_t/new_masses) !unguided version

    !calculate the local energy-------------------------------------------------
    r_temp = matmul(ev,this%r)
    rxyz = reshape(r_temp,shape(rxyz)) + rxyz0 !convert to kartesian corrdinated for energy calculation

    call energyandforces(nat, rxyz, fxyz, epot, alat, deralat, atomnames, thread_num)
    this%etot = epot

    !calculate the weights------------------------------------------------------
    if (-delta_t*((this%etot + etot_prev)/2 - et) .ge. 10) then
      print*, "weight over/underflow: ", (-delta_t*((this%etot + etot_prev)/2 - et)), this%etot, etot_prev
    else
      this%weight = this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))
    end if
    !print*, "epot, ekin, etot: ", epot, el_kin(this%r), this%etot
    !print*, "weight: ", this%weight

  end subroutine propagate_unguided

  !calculate the drift velocity grad(psi)/psi-----------------------------------
  function drift(r,delta_t) result(v)
    implicit none
    real(8), dimension(dim), intent(in) :: r
    real(8), dimension(dim) :: v
    real(8) :: dot,delta_t
    integer :: i

    do i=1,dim
      v(i) = -r(i)/sigma2(i)
    end do

    !julich.pdf, page 14
    dot = dot_product(v,v)
    if (dot .gt. 1.d-8) then
      v = (-1.d0 + sqrt(1.d0 + 2.d0*dot*delta_t))/(dot*delta_t)*v
    end if

  end function drift

  !calculate the kinetic part of the local energy psi*H/psi---------------------
  real(8) function el_kin(r)
    implicit none
    real(8), dimension(dim), intent(in) :: r
    real(8), dimension(fulldim) :: r_temp, r_plus, r_minus
    real(8) :: h
    integer :: i

    el_kin = 0
    do i=1,dim
      el_kin = el_kin - 0.5/sigma2(i)*(r(i)**2/sigma2(i) - 1)/new_masses(i)
    end do

  end function el_kin

  real(8) function psi(r)
    implicit none
    real(8), dimension(dim), intent(in) :: r
    integer :: i

    psi = 0
    do i=1,dim
      psi = psi - r(i)**2/(2*sigma2(i))
    end do
    psi = exp(psi)

  end function psi

  subroutine write_walker_kartesian(this, iounit)
    implicit none
    class(walker) :: this
    integer :: iounit, i
    real(8), dimension(fulldim) :: r_temp
    real(8), dimension(3,nat) :: rxyz

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

  !Initalize walkers in the equilibrium position with weight 1------------------
  subroutine initalizeWalkers(walkers,no_walkers,etot)
    implicit real(8) (a-h,o-z)
    type(walker), dimension(no_walkers) :: walkers

    do i=1,no_walkers
      call create_walker(walkers(i), etot, 1.d0)
    end do

  end subroutine initalizeWalkers

end module dftbp_walker

!creates an array of n normal distributed random numbers using
!the Box-Muller algorithm.------------------------------------------------------
subroutine normal(rand, n)
  implicit real(8) (a-h,o-z)
  real(8), dimension(n) :: rand
  pi = 4.d0*atan(1.d0)

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
    call RANDOM_NUMBER(r2) !extra random number
    rand(n) = sqrt(-2.0*log(rand(n)))*cos(2*pi*r2)
  end if

end subroutine normal
