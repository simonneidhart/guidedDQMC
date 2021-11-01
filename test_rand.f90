program test_rand

  implicit none
  integer :: nat = 100000, i, ios
  real(8) :: rand(100000)


  call normal(rand,nat)
  open(unit=13, file="nrand", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,nat
    write(unit=13, fmt=*) rand(i)
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"

end program test_rand

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
