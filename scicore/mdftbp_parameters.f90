module parameters

  implicit none
  integer, parameter :: nat = 8, dim = 18
  real(8), dimension(3,nat) :: rxyz0, sigma2
  character(len=2) :: atomnames(nat)
  real(8), dimension(nat) :: masses
  real(8), parameter ::  delta_t = 10.0
  real(8), parameter :: Bohr_Ang = 0.529177d0


contains

  subroutine read_parameters()
    implicit none
    integer :: i,j, ios

    print*, "delta_t = ", delta_t

    masses = (/12.011, 12.011, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007/)
    !masses = masses*1836.1

    open(unit=16, file="geom_opt.xyz", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,nat
      read(16,*) atomnames(i), ( rxyz0(j, i), j = 1, 3 )
    end do

    rxyz0=rxyz0/Bohr_Ang

    write(unit=*, fmt=*) "sigma^2: "
    do i=1,nat
      do j=1,3
        sigma2(j,i) = 1/masses(i)
      end do
      write(unit=*, fmt=*) atomnames(i), ( sigma2(j, i), j = 1, 3 )
    end do

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    print*, "masses: ", masses


  end subroutine read_parameters



end module parameters
