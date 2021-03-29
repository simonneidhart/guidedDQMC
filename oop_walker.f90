module dftbp_walker
  integer, private, parameter :: natt = 4
  type walker
    real(8), dimension(3,natt) :: rxyz, fxyz
    real(8) :: etot, weight
  contains
    procedure :: create_walker
    procedure :: propagate
    procedure :: write_walker

  end type walker

contains

  subroutine create_walker(this, rxyz_, fxyz_, etot_, weight_)
    implicit none
    class(walker) :: this
    real(8), dimension(3,natt) :: rxyz_, fxyz_
    real(8) :: etot_, weight_

    this%rxyz = rxyz_
    this%fxyz = fxyz_
    this%etot = etot_
    this%weight = weight_


  end subroutine create_walker

  subroutine propagate(this, et, delta_t, alat, deralat, atomnames)
    implicit none
    class(walker) :: this
    real(8), dimension(3,natt) :: normal_rand
    real(8), dimension(3,3) :: alat, deralat
    real(8) :: delta_t, etot_prev, et
    character(len=2) :: atomnames(natt)

    call normal(normal_rand,natt)

    etot_prev = this%etot
    !move the walkers according to the Gaussian distribution------------------
    this%rxyz = this%rxyz + normal_rand*sqrt(delta_t) + this%fxyz*delta_t
    call energyandforces_dummy(natt, this%rxyz, this%fxyz, this%etot, alat, deralat, atomnames)
    !calculate the weights----------------------------------------------------
    this%weight = this%weight*exp(-delta_t*((this%etot + etot_prev)/2 - et))

  end subroutine propagate

  subroutine write_walker(this, iounit)
    implicit none
    class(walker) :: this
    integer :: iounit, i

    do i=1,natt
      write(unit=iounit, fmt=*) this%rxyz(:,i)
    end do
    write(unit=iounit, fmt=*) this%weight

  end subroutine write_walker

end module dftbp_walker
