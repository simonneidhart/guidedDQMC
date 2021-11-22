module dftbp_walker
  use dftbp_parameters
  type walker
    !displacements from the minimum energy position
    real(8), allocatable, dimension(:) :: r
    real(8) :: etot, weight
  end type walker

contains

  !initalize a new walker with displacements distributed according to the
  !trial wave function and a given weight and energy
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

  !create a walkers with given displacements, energy and weight
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

  !propagate a walker and calculate the weight
  subroutine propagate(this, et, delta_t, thread_num, n, n_acc)
    use mod_dftbp
    use dftbp_parameters
    implicit none
    integer :: i, thread_num, n, n_acc
    class(walker) :: this
    real(8), dimension(3,nat) :: rxyz, fxyz
    real(8), dimension(dim) :: normal_rand, r_prop
    real(8), dimension(fulldim) :: r_temp
    real(8) :: delta_t, etot_prev, et, epot, p_acc, rand

    call normal(normal_rand,dim)

    etot_prev = this%etot !save the energy before the move
    !move the walkers according to the Gaussian distribution and the drift velocity
    r_prop = this%r + normal_rand*sqrt(delta_t) + drift(this%r,delta_t)*delta_t


    if (n .gt. 20) then !no accept- reject step in the first 20 steps
      !calculate the accepance probability
      p_acc = min(1.d0,exp(-dot_product(this%r-r_prop-drift(r_prop,delta_t)*&
      delta_t,this%r-r_prop-drift(r_prop,delta_t)*delta_t))/&
      exp(-dot_product(r_prop-this%r-drift(this%r,delta_t)*delta_t,&
      r_prop-this%r-drift(this%r,delta_t)*delta_t))*psi(r_prop)**2/psi(this%r)**2)
    else
      p_acc = 1.d0 !always accept in the first 20 steps
    end if

    call RANDOM_NUMBER(rand)

    if (rand .le. p_acc) then !accept
      n_acc = n_acc + 1 !count the accepted steps
      this%r = r_prop

      !calculate the local energy-------------------------------------------------
      r_temp = matmul(ev,this%r)
      rxyz = reshape(r_temp,shape(rxyz)) + rxyz0 !convert to kartesian corrdinated for energy calculation

      call energyandforces(nat, rxyz, fxyz, epot, alat, deralat, atomnames, thread_num)
      this%etot = epot + el_kin(this%r)

      !calculate the weights------------------------------------------------------
      if (-delta_t*((this%etot + etot_prev)/2 - et) .ge. 10) then
        print*, "weight overflow: ", this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))
        print*, "etot: ", this%etot
        print*, "thread number: ", thread_num
        print*, "rxyz: "
        call write_matrix(rxyz,3,nat)
        print*, "r: "
        call write_matrix(this%r,dim,1)
        stop
      else
        !update the weight
        this%weight = this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))
      end if
    end if

  end subroutine propagate

  !unguided version of the update step
  subroutine propagate_unguided(this, et, delta_t, thread_num)
    use mod_dftbp
    use dftbp_parameters
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
    this%r = this%r + normal_rand*sqrt(delta_t/new_masses) !unguided version

    !calculate the local energy-------------------------------------------------
    r_temp = matmul(ev,this%r)
    rxyz = reshape(r_temp,shape(rxyz)) + rxyz0 !convert to kartesian corrdinated for energy calculation

    call energyandforces(nat, rxyz, fxyz, epot, alat, deralat, atomnames, thread_num)
    this%etot = epot

    !calculate the weights------------------------------------------------------
    if (-delta_t*((this%etot + etot_prev)/2 - et) .ge. 10) then
      print*, "weight overflow: ", this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))
      print*, "etot: ", this%etot
      print*, "thread number: ", thread_num
      print*, "rxyz: "
      call write_matrix(rxyz,3,nat)
      print*, "r: "
      call write_matrix(this%r,dim,1)
      stop
    else
      !update the weight
      this%weight = this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))
    end if

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
    !transformation of the drift velocity
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

    el_kin = 0.d0
    do i=1,dim
      el_kin = el_kin - 0.5d0/sigma2(i)*(r(i)**2/sigma2(i) - 1.d0)/new_masses(i)
    end do

  end function el_kin

  !evaluate the trial wave function
  real(8) function psi(r)
    implicit none
    real(8), dimension(dim), intent(in) :: r
    integer :: i

    psi = 0.d0
    do i=1,dim
      psi = psi - r(i)**2/(2.d0*sigma2(i))
    end do
    psi = exp(psi)

  end function psi

  !write the cartesian coordinates of a walker to a file------------------------
  subroutine write_walker_cartesian(this, iounit)
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

  end subroutine write_walker_cartesian

  !write the walker displacements to a file-------------------------------------
  subroutine write_walker_modes(this, iounit)
    implicit none
    class(walker) :: this
    integer :: iounit, i

    do i=1,dim
      write(unit=iounit, fmt=*) this%r(i)
    end do

  end subroutine write_walker_modes

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
