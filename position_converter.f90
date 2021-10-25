program pos_conv
  implicit real(8) (a-h,o-z)
  integer, parameter :: nat = 8
  real(8) :: rxyz(3,nat)
  character(len=200):: line1, line2, l3, l4


  open(unit=16, file="et_noWalkers", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 16"
  open(unit=17, file="et_noWalkers_new", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 17"

  no_walkers = 626

  !do k=1,no_walkers
  !  read(16,*) weight
  !  do i=1,nat
  !    read(16,*) ( rxyz(j, i), j = 1, 3 )
  !    write(unit=17, fmt=*) ( rxyz(j, i), j = 1, 3 )
  !  end do
  !end do

  read(16,*) line1 ,l3, l4
  write(unit=17, fmt=*) line1, l3, l4
  do i=1,2000
    read(16,*) line1, l3, l4
    write(unit=17, fmt=*) line1, l3, l4
    read(16,*) line2
  end do

  close(unit=16, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 16"

  close(unit=17, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 17"



end program pos_conv
