program pos_conv
  implicit real(8) (a-h,o-z)
  integer, parameter :: nat = 8
  real(8) :: rxyz(3,nat)


  open(unit=16, file="walker_positions", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 16"
  open(unit=17, file="walker_positions_new", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 17"

  no_walkers = 626

  do k=1,no_walkers
    read(16,*) weight
    do i=1,nat
      read(16,*) ( rxyz(j, i), j = 1, 3 )
      write(unit=17, fmt=*) ( rxyz(j, i), j = 1, 3 )
    end do
  end do

  close(unit=16, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 16"

  close(unit=17, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 17"



end program pos_conv
