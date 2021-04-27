module parameters

  implicit none
  integer, parameter :: nat = 8, dim = 18
  real(8), dimension(dim) :: mu, sigma2, k_pot
  real(8), dimension(24,dim) :: ev
  real(8) :: average_mass = 6900.d0

contains

  subroutine read_parameters()
    implicit none
    integer :: i,j, ios
    real(8) :: e0_theory

    open(unit=16, file="inital", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,dim
      read(16,*) mu(i)
    end do

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    open(unit=16, file="k_pot", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,dim
      read(16,*) k_pot(i)
      sigma2(i) = 1/sqrt(k_pot(i))
    end do
    print*, "sigma2 from k_pot, modifier = 1"
    write(unit=*, fmt=*) "sigma^2= ", sigma2

    e0_theory = 0
    do i=1,dim
      e0_theory = e0_theory + 0.5*sqrt(k_pot(i))
    end do
    print*, "e0 theory = ", e0_theory

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    open(unit=16, file="eigenvectors", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,24
      read(16,*) ( ev(i, j), j = 1, dim )
    end do

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"


  end subroutine read_parameters



end module parameters
