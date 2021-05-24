program read_hessian
  implicit real(8) (a-h,o-z)
  integer, parameter ::  nat = 8, dim = 18, fulldim = 24
  real(8), dimension(fulldim,fulldim) :: hess1, hess2, dot1, dot2
  real(8), dimension(fulldim) :: eigenvalues1, eigenvalues2
  real(8), dimension(144,4) :: hess_in
  real(8), allocatable, dimension(:) :: work
  integer :: info


  open(unit=16, file="hessian.out", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 16"

  read(16,*) hess_in

  !do i=1,144
  !  read(16,"(F11.10, F11.10, F11.10, F11.10)") ( hess_in(i, j), j = 1, 4)
  !end do

  close(unit=16, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 16"

  hess1 = reshape(hess_in,shape(hess1))
  hess2 = reshape(transpose(hess_in),shape(hess2))

  print*, "hess1: "
  call write_matrix(hess1,fulldim)

  print*, "hess2: "
  call write_matrix(hess2,fulldim)

  lwork = max(1,3*fulldim-1)
  allocate(work(lwork))

  !call dsyev('V','U',fulldim,hess1,fulldim,eigenvalues1,work,lwork,info)
  !if(info .neqv. 0) print*, "hess1 failed ", info

  !call dsyev('V','U',fulldim,hess2,fulldim,eigenvalues2,work,lwork,info)
  !if(info .neqv. 0) print*, "hess2 failed ", info

  do i=1,fulldim
    do j=1,fulldim
      dot1(i,j) = dot_product(hess1(:,i),hess1(:,j))
      dot2(i,j) = dot_product(hess2(:,i),hess2(:,j))
    end do
  end do


  print*, "eigenvalues1: "
  print*, eigenvalues1

  print*, "eigenvectors1: "
  call write_matrix(hess1,fulldim)

  print*, "dot1: "
  call write_matrix(dot1,fulldim)

  print*, "eigenvalues2: "
  print*, eigenvalues2

  print*, "eigenvectors2: "
  call write_matrix(hess2,fulldim)

  print*, "dot2: "
  call write_matrix(dot2,fulldim)



end program read_hessian

subroutine write_matrix(a,n)
   implicit none
   real(8), dimension(n,n) :: a
   integer :: n,i,j
   write(*,*)

   do i = 1,n
      write(*,"(9999(F6.4,:,','))") (a(i,j), j = 1,n)
   end do
   write(*,*)
end subroutine write_matrix
