module parameters
  implicit none
  integer, parameter :: nat = 8, dim = 18, fulldim = 24
  real(8), parameter :: Bohr_Ang = 0.529d0
  real(8), dimension(dim) :: mu = 0.d0, sigma2, new_masses
  real(8), dimension(fulldim,dim) :: ev
  real(8), dimension(3,3) :: alat, deralat
  real(8), dimension(3,nat) :: rxyz0
  character(len=2) :: atomnames(nat)

contains

  !reads the hessian output from DFTB+ and calculates the derrived quantities
  !for the DQMC calculation.
  subroutine read_hessian(e0_theory)
    implicit none
    real(8), dimension(fulldim,fulldim) :: hess, hess_m, dot, eigenvectors_full
    real(8), dimension(fulldim) :: masses, eigenvalues_full, eigenvalues_full_m
    real(8), dimension(144,4) :: hess_in
    real(8), allocatable, dimension(:) :: work
    real(8) :: e0_theory
    integer :: i, j, ios, lwork, info

    masses = (/12.011, 12.011, 12.011, 12.011, 12.011, 12.011, &
    1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, &
    1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007/)
    masses = masses*1836.1

    alat(1,1)=1.d0
    alat(2,1)=0.d0
    alat(3,1)=0.d0

    alat(1,2)=0.d0
    alat(2,2)=1.d0
    alat(3,2)=0.d0

    alat(1,3)=0.d0
    alat(2,3)=0.d0
    alat(3,3)=1.d0

    open(unit=16, file="hessian", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    read(16,*) hess_in

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    hess = reshape(hess_in,shape(hess))

    do i=1,fulldim
      do j=1,fulldim
        hess_m(i,j) = hess(i,j)/sqrt(masses(i)*masses(j)) !scale the hessian with the masses
      end do
    end do

    hess = 0.5d0 * (hess + transpose(hess)) !symmetrize the hessian matrix
    hess_m = 0.5d0 * (hess_m + transpose(hess_m)) !symmetrize the mass-scaled hessian matrix

    print*, "hessian: "
    call write_matrix(hess,fulldim,fulldim)

    lwork = max(1,3*fulldim-1)
    allocate(work(lwork))

    call dsyev('V','U',fulldim,hess_m,fulldim,eigenvalues_full_m,work,lwork,info)
    if(info .neqv. 0) print*, "Mass-scaled hessian diagonalization failed: ", info

    call dsyev('V','U',fulldim,hess,fulldim,eigenvalues_full,work,lwork,info)
    if(info .neqv. 0) print*, "Hessian diagonalization failed: ", info

    deallocate(work)

    do i=1,fulldim
      do j=1,fulldim
        dot(i,j) = dot_product(hess(:,i),hess(:,j))
      end do
    end do

    print*, "Eigenvalues of the hessian matrix: "
    call write_matrix(eigenvalues_full,fulldim,1)

    print*, "Eigenvectors of the hessian matrix: "
    call write_matrix(hess,fulldim,fulldim)
    ev = hess(:,7:fulldim)

    print*, "Orthonormality of the eigenvectors: "
    call write_matrix(dot,fulldim,fulldim)

    open(unit=12, file="geom_opt.xyz", iostat=ios, status="old")
    if ( ios /= 0 ) stop "Error opening file "

    do i = 1, nat
       read(12, *) atomnames(i), ( rxyz0(j, i), j = 1, 3)
    end do

   !Transofrm angstroem input to atomic units
   alat=alat/Bohr_Ang
   rxyz0=rxyz0/Bohr_Ang

   close(unit=12, iostat=ios, status="keep")
   if ( ios /= 0 ) stop "Error closing file unit 16"

   print*, "rxyz0: "
   call write_matrix(rxyz0,3,nat)

   !calculate the masses in the new coordinate system
   do i=1,dim
     new_masses(i) = eigenvalues_full(i+6)/eigenvalues_full_m(i+6)
   end do

   print*, "new_masses: "
   call write_matrix(new_masses,dim,1)

   !calculate the width of the gaussian for the trial wave function
   do i=1,dim
     sigma2(i) = 2.d0/(sqrt(eigenvalues_full(i+6)*new_masses(i)))
   end do
   print*, "sigma2 from ew, modifier = 2.0"
   call write_matrix(sigma2,dim,1)

   !calculate the theoretical value for e0 based only on the harmonic approximation
   e0_theory = 0
   do i=1,dim
     e0_theory = e0_theory + 0.5*sqrt(eigenvalues_full(i+6)/new_masses(i))
   end do

 end subroutine read_hessian

!prints an n x m matrix to the console
  subroutine write_matrix(a,n,m)
     implicit none
     real(8), dimension(n,m) :: a
     integer :: n,m,i,j
     write(*,*)

     do i = 1,n
        write(*,"(9999(G12.4,:,','))") (a(i,j), j = 1,m)
     end do
     write(*,*)
  end subroutine write_matrix

  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  function inv(A) result(Ainv)
    real(8), dimension(:,:), intent(in) :: A
    real(8), dimension(size(A,1),size(A,2)) :: Ainv

    real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function inv

! old subroutine that reads ev, ew, masses form files calculated by MATLAB
  subroutine read_parameters()
    implicit none
    real(8), dimension(fulldim) :: masses
    integer :: i,j, ios
    real(8) :: e0_theory

    masses = (/12.011, 12.011, 12.011, 12.011, 12.011, 12.011, &
    1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, &
    1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007/)
    masses = masses*1836.1

    open(unit=16, file="new_masses", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,dim
      read(16,*) new_masses(i)
    end do

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    open(unit=16, file="inital_mass", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,dim
      read(16,*) mu(i)
    end do

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    open(unit=16, file="k_pot", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,dim
      !read(16,*) k_pot(i)
      !sigma2(i) = 1/(sqrt(k_pot(i)*new_masses(i)))
    end do
    print*, "sigma2 from pot, modifier = 1"
    write(unit=*, fmt=*) "sigma^2= ", sigma2

    e0_theory = 0
    do i=1,dim
      e0_theory = e0_theory + 0.5/(sigma2(i)*sqrt(new_masses(i)))
    end do
    print*, "e0 theory = ", e0_theory

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    open(unit=16, file="eigenvectors_mass", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,fulldim
      read(16,*) ( ev(i, j), j = 1, dim )
    end do

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

  end subroutine read_parameters

end module parameters
