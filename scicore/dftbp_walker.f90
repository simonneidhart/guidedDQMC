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

    !allocate(this%rxyz(3,nat_))
    !allocate(this%fxyz(3,nat_))
    !this%nat = nat_
    this%r = r_
    this%etot = etot_ + el_kin(mu)
    this%weight = weight_


  end subroutine create_walker

  subroutine propagate(this, et, delta_t, alat, deralat, atomnames)
    use mod_dftbp
    use parameters
    implicit none
    integer :: i
    class(walker) :: this
    real(8), dimension(3,nat) :: rxyz, fxyz
    real(8), dimension(dim) :: normal_rand
    real(8), dimension(24) :: r_temp
    real(8), dimension(3,3) :: alat, deralat
    real(8) :: delta_t, etot_prev, et, epot
    character(len=2) :: atomnames(nat)

    !print*, nat

    call normal(normal_rand,dim)

    etot_prev = this%etot
    !move the walkers according to the Gaussian distribution------------------
    !print*, shape(normal_rand), shape(this%rxyz), shape(this%fxyz)
    this%r = this%r + normal_rand*sqrt(delta_t) + drift(this%r)*delta_t
    !calculate the local energy------------------------------------------------
    r_temp = 0
    do i=1,dim
      r_temp = r_temp + ev(:,i)*this%r(i)
    end do
    rxyz = reshape(r_temp,shape(rxyz))
    call energyandforces(nat, rxyz, fxyz, epot, alat, deralat, atomnames)
    this%etot = epot + el_kin(this%r)
    !calculate the weights----------------------------------------------------
    print*, "walker epot, ekin, weight: ", epot, el_kin(this%r), this%weight
    this%weight = this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))
    !print*, "Weight: ", this%weight

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
    implicit real(8) (a-h,o-z)
    real(8), dimension(dim), intent(in) :: r
    real(8) :: m = 1836.1d0

    el_kin = 0
    do i=1,dim
      el_kin = el_kin - 0.5/sigma2(i)*((r(i)-mu(i))**2/sigma2(i) - 1)/m
    end do


  end function el_kin

  subroutine write_walker(this, iounit)
    implicit none
    class(walker) :: this
    integer :: iounit, i
    real(8), dimension(24) :: r_temp
    real(8), dimension(3,nat) :: rxyz

    r_temp = 0
    do i=1,dim
      r_temp = r_temp + ev(:,i)*this%r(i)
    end do
    rxyz = reshape(r_temp,shape(rxyz))

    !write(unit=iounit, fmt=*) "weight: ", this%weight, "etot: ", this%etot
    do i=1,nat
      write(unit=iounit, fmt=*) rxyz(:,i)
    end do

  end subroutine write_walker

end module dftbp_walker
