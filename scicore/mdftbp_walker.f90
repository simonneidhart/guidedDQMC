module dftbp_walker
  use parameters
  type walker
    real(8), dimension(3,nat) :: rxyz
    real(8) :: etot, weight
  contains
    procedure :: create_walker
    procedure :: propagate
    procedure :: write_walker

  end type walker

contains

  subroutine create_walker(this, rxyz_, epot_, weight_)
    implicit none
    class(walker) :: this
    real(8), dimension(3,nat) :: rxyz_
    real(8) :: epot_, weight_

    !allocate(this%rxyz(3,nat_))
    !allocate(this%fxyz(3,nat_))
    !this%nat = nat_
    this%rxyz = rxyz_
    this%etot = epot_ + el_kin(rxyz_)
    this%weight = weight_


  end subroutine create_walker

  subroutine propagate(this, et, alat, deralat)
    use mod_dftbp
    use parameters
    implicit none
    integer :: i
    class(walker) :: this
    real(8), dimension(3,nat) :: fxyz, rxyz_in
    real(8), dimension(3,nat) :: normal_rand
    real(8), dimension(3,3) :: alat, deralat
    real(8) :: et, epot, etot_prev

    call normal(normal_rand,nat)

    etot_prev = this%etot
    !move the walkers according to the Gaussian distribution------------------
    this%rxyz = this%rxyz + normal_rand*sqrt(delta_t) + drift(this%rxyz)*delta_t
    !calculate the local energy------------------------------------------------
    call energyandforces(nat, this%rxyz, fxyz, epot, alat, deralat, atomnames)
    this%etot = epot + el_kin(this%rxyz)
    !calculate the weights----------------------------------------------------
    print*, "walker epot, ekin, weight: ", epot, el_kin(this%rxyz), this%weight
    this%weight = this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))
    !print*, "walker epot, weight: ", epot, this%weight

  end subroutine propagate

  !calculate the drift velocity grad(psi)/psi-------------------------------------
  function drift(r) result(v)
    implicit real(8) (a-h,o-z)
    real(8), dimension(3,nat), intent(in) :: r
    real(8), dimension(3,nat) :: v

    do i=1,nat
      do j=1,3
        v(j,i) = -(r(j,i)-rxyz0(j,i))/sigma2(j,i)
      end do
    end do

  end function drift

  !calculate the kinetic part of the local energy psi*H/psi-------------------------------------------
  real(8) function el_kin(r)
    implicit real(8) (a-h,o-z)
    real(8), dimension(3,nat), intent(in) :: r

    el_kin = 0
    do i=1,nat
      do j=1,3
        !el_kin = el_kin - 0.5/sigma2(j,i)*((r(j,i)-rxyz0(j,i))**2/sigma2(j,i) - 1)/masses(i)
      end do
    end do


  end function el_kin

  subroutine write_walker(this, iounit)
    implicit none
    class(walker) :: this
    integer :: iounit, i
    real(8), dimension(3,nat) :: rxyz

    !write(unit=iounit, fmt=*) "weight: ", this%weight, "etot: ", this%etot
    do i=1,nat
      write(unit=iounit, fmt=*) this%rxyz(:,i)
    end do

  end subroutine write_walker

end module dftbp_walker
