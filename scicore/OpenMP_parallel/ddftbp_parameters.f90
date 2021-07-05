module parameters
  implicit none
  integer, parameter :: nat = 62, dim = 180, fulldim = 186
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
    real(8) :: e0_theory, m_in
    integer :: i, j, ios, lwork, info

    open(unit=16, file="masses", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    do i=1,nat
      read(16,*) m_in
      masses(3*i - 2) = m_in
      masses(3*i - 1) = m_in
      masses(3*i) = m_in
    end do

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

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

    open(unit=16, file="hessian.out", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    read(16,*) hess

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    !print*, "masses: "
    !call write_matrix(masses,fulldim,1)


    do i=1,fulldim
      do j=1,fulldim
        hess_m(i,j) = hess(i,j)/sqrt(masses(i)*masses(j)) !scale the hessian with the masses
      end do
    end do

    hess = 0.5d0 * (hess + transpose(hess)) !symmetrize the hessian matrix
    hess_m = 0.5d0 * (hess_m + transpose(hess_m)) !symmetrize the mass-scaled hessian matrix

    !print*, "hessian: "
    !call write_matrix(hess,fulldim,fulldim)

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
    !call write_matrix(eigenvalues_full,fulldim,1)

    !print*, "Eigenvectors of the hessian matrix: "
    !call write_matrix(hess,fulldim,fulldim)
    ev = hess(:,7:fulldim)

    !print*, "Orthonormality of the eigenvectors: "
    !call write_matrix(dot,fulldim,fulldim)

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

   !print*, "rxyz0: "
   !call write_matrix(rxyz0,3,nat)

   !calculate the masses in the new coordinate system
   do i=1,dim
     new_masses(i) = eigenvalues_full(i+6)/eigenvalues_full_m(i+6)
   end do

   !print*, "new_masses: "
   !call write_matrix(new_masses,dim,1)

   !calculate the width of the gaussian for the trial wave function
   do i=1,dim
     sigma2(i) = 2.d0/(sqrt(eigenvalues_full(i+6)*new_masses(i)))
   end do
   print*, "sigma2 from ew, modifier = 2.0"
   !call write_matrix(sigma2,dim,1)

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

end module parameters
