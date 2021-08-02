module parameters
  implicit none
  integer :: nat, dim, fulldim
  real(8), parameter :: Bohr_Ang = 0.529d0
  real(8), allocatable, dimension(:) :: sigma2, new_masses
  real(8), allocatable, dimension(:,:) :: ev
  real(8), dimension(3,3) :: alat, deralat
  real(8), allocatable, dimension(:,:) :: rxyz0
  character(len=2), allocatable, dimension(:) :: atomnames

contains

  !reads the hessian output from DFTB+ and calculates the derrived quantities
  !for the DQMC calculation.
  subroutine read_hessian(e0_theory, debug_extended)
    implicit none
    real(8), allocatable, dimension(:,:) :: hess, hess_m, dot, eigenvectors_full
    real(8), allocatable, dimension(:) :: masses, eigenvalues_full, eigenvalues_full_m
    real(8), allocatable, dimension(:) :: work
    real(8) :: e0_theory, m_in
    integer :: i, j, ios, lwork, info
    logical :: debug_extended

    open(unit=16, file="masses", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    read(16,*) nat
    fulldim = nat*3
    dim = fulldim - 6
    write(*,'(A,I5,A,I5,A,I5)'),"Nat = ", nat," 3*Nat = ", fulldim, " 3*Nat - 6 = ", dim

    allocate(rxyz0(3,nat))
    allocate(atomnames(nat))
    allocate(ev(fulldim,dim))
    allocate(sigma2(dim))
    allocate(new_masses(dim))
    allocate(hess(fulldim,fulldim))
    allocate(hess_m(fulldim,fulldim))
    allocate(dot(fulldim,fulldim))
    allocate(eigenvectors_full(fulldim,fulldim))
    allocate(masses(fulldim))
    allocate(eigenvalues_full(fulldim))
    allocate(eigenvalues_full_m(fulldim))

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

    if (debug_extended) then
      print*, "masses: "
      call write_matrix(masses,fulldim,1)
    end if


    do i=1,fulldim
      do j=1,fulldim
        hess_m(i,j) = hess(i,j)/sqrt(masses(i)*masses(j)) !scale the hessian with the masses
      end do
    end do

    hess = 0.5d0 * (hess + transpose(hess)) !symmetrize the hessian matrix
    hess_m = 0.5d0 * (hess_m + transpose(hess_m)) !symmetrize the mass-scaled hessian matrix

    if (debug_extended) then
      print*, "hessian: "
      call write_matrix(hess,fulldim,fulldim)
    end if

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

    if (debug_extended) then
      print*, "Eigenvalues of the hessian matrix: "
      call write_matrix(eigenvalues_full,fulldim,1)
      print*, "Eigenvectors of the hessian matrix: "
      call write_matrix(hess,fulldim,fulldim)
      print*, "Orthonormality of the eigenvectors: "
      call write_matrix(dot,fulldim,fulldim)
    end if

    ev = hess(:,7:fulldim)

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

   !calculate the masses in the new coordinate system
   do i=1,dim
     new_masses(i) = eigenvalues_full(i+6)/eigenvalues_full_m(i+6)
   end do

   if (debug_extended) then
      print*, "rxyz0: "
      call write_matrix(rxyz0,3,nat)
      print*, "new_masses: "
      call write_matrix(new_masses,dim,1)
   end if

   !calculate the width of the gaussian for the trial wave function
   do i=1,dim
     sigma2(i) = 1.d0/(sqrt(eigenvalues_full(i+6)*new_masses(i)))
   end do
   print*, "Sigma2 modifier = 1.0"
   if (debug_extended) then
      call write_matrix(sigma2,dim,1)
   end if

   !calculate the theoretical value for e0 based only on the harmonic approximation
   e0_theory = 0
   do i=1,dim
     e0_theory = e0_theory + 0.5*sqrt(eigenvalues_full(i+6)/new_masses(i))
   end do

   !call subroutine write_node_file(eigenvalues_full)

   call write_new_masses_and_eigenvalues(eigenvalues_full)

 end subroutine read_hessian

 subroutine write_new_masses_and_eigenvalues(eigenvalues_full)
   implicit none
   integer :: i, j, ios
   real(8), dimension(fulldim) :: eigenvalues_full

   open(unit=13, file="new_masses", iostat=ios, action="write")
   if ( ios /= 0 ) stop "Error opening file unit 13"
   do i=1,dim
     write(unit=13, fmt=*) new_masses(i)
   end do
   close(unit=13, iostat=ios, status="keep")
   if ( ios /= 0 ) stop "Error closing file unit 13"

   open(unit=13, file="eigenvalues", iostat=ios, action="write")
   if ( ios /= 0 ) stop "Error opening file unit 13"
   do i=1,dim
     write(unit=13, fmt=*) eigenvalues_full(i+6)
   end do
   close(unit=13, iostat=ios, status="keep")
   if ( ios /= 0 ) stop "Error closing file unit 13"

 end subroutine write_new_masses_and_eigenvalues

 subroutine write_node_file(eigenvalues_full)
   implicit none
   integer :: i, j, ios
   real(8), dimension(fulldim) :: eigenvalues_full

   !output mode file
   open(unit=13, file="modes.txt", iostat=ios, action="write")
   if ( ios /= 0 ) stop "Error opening file unit 13"
   write(unit=13, fmt=*) nat, dim
   do i=1,dim
     write(unit=13, fmt=*) ' '
     write(unit=13, fmt=*) "N ", eigenvalues_full(i+6), "Mode ", i
     do j=1,fulldim
       write(unit=13, fmt=*) ev(j,i)
     end do
   end do
   write(unit=13, fmt=*) "END"
   close(unit=13, iostat=ios, status="keep")
   if ( ios /= 0 ) stop "Error closing file unit 13"

 end subroutine write_node_file

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
