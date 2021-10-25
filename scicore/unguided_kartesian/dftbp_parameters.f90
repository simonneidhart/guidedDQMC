module parameters

  implicit none
  integer :: nat
  real(8), allocatable, dimension(:,:) :: rxyz0, sigma2
  character(len=2), allocatable, dimension(:) :: atomnames
  real(8), allocatable, dimension(:) :: masses
  real(8), parameter ::  delta_t = 10.0
  real(8), parameter :: Bohr_Ang = 0.529177d0

contains

  subroutine read_parameters()
    implicit none
    integer :: i,j, ios
    real(8) :: nat_

    print*, "delta_t = ", delta_t

    !the masses file contains the number of atoms and all the nat atom masses
    open(unit=16, file="masses", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    read(16,*) nat
    allocate(atomnames(nat))
    allocate(rxyz0(3,nat))
    allocate(sigma2(3,nat))
    allocate(masses(nat))

    do i=1,nat
      read(16,*) masses(i)
    end do

    masses = masses*1822.888486

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    open(unit=16, file="geom_opt.xyz", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,nat
      read(16,*) atomnames(i), ( rxyz0(j, i), j = 1, 3 )
    end do

    rxyz0=rxyz0/Bohr_Ang

    write(unit=*, fmt=*) "sigma^2: "
    do i=1,nat
      do j=1,3
        sigma2(j,i) = delta_t/masses(i)
      end do
      write(unit=*, fmt=*) atomnames(i), ( sigma2(j, i), j = 1, 3 )
    end do

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    print*, "masses: ", masses

    write(unit=*, fmt=*) "rxyz0: "
    do i=1,nat
      write(unit=*, fmt=*) atomnames(i), ( rxyz0(j, i), j = 1, 3 )
    end do

  end subroutine read_parameters

end module parameters
