module parameters
  use precision
  use atomicStructure
  use RuNNerInterface
  implicit none
  integer :: nat, dim, fulldim
  real(dp), allocatable, dimension(:) :: sigma2, new_masses
  real(dp), allocatable, dimension(:,:) :: ev
  real(dp), allocatable, dimension(:,:) :: rxyz0

contains

  !reads the hessian output from DFTB+ and calculates the derrived quantities
  !for the DQMC calculation.
  subroutine read_hessian(ats, e0_theory, debug_extended)
    use precision
    use atomicStructure
    use RuNNerInterface
    implicit none
    real(dp), allocatable, dimension(:,:) :: hess, hess_m, dot, eigenvectors_full, fxyz
    real(dp), allocatable, dimension(:) :: masses, eigenvalues_full, eigenvalues_full_m
    real(dp), allocatable, dimension(:) :: work
    type(atStruct) :: ats
    real(dp) :: e0_theory, m_in, etot
    integer(4) :: i, j, ios, lwork, info
    logical :: debug_extended

    open(unit=16, file="masses.txt", iostat=ios)
    if ( ios /= 0 ) stop "Error opening file unit 16"

    read(16,*) nat
    fulldim = nat*3
    dim = fulldim - 6
    write(*,'(A,I5,A,I5,A,I5)') "Nat = ", nat," 3*Nat = ", fulldim, " 3*Nat - 6 = ", dim

    allocate(rxyz0(3,nat))
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
    allocate(fxyz(3,nat))

    do i=1,nat
      read(16,*) m_in
      masses(3*i - 2) = m_in
      masses(3*i - 1) = m_in
      masses(3*i) = m_in
    end do
    masses = masses*1836.1

    close(unit=16, iostat=ios, status="keep")
    if ( ios /= 0 ) stop "Error closing file unit 16"

    call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, etot, fxyz)
    print*, "success: ", etot

    call calculate_hessian_forces(ats, hess)

    if (debug_extended) then
      print*, "masses: "
      !call write_matrix(masses,fulldim,1)
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
      call write_matrix(hess,3,3)
    end if

    lwork = max(1,3*fulldim-1)
    allocate(work(lwork))

    call dsyev('V','U',fulldim,hess_m,fulldim,eigenvalues_full_m,work,lwork,info)
    if(info /= 0) print*, "Mass-scaled hessian diagonalization failed: ", info

    call dsyev('V','U',fulldim,hess,fulldim,eigenvalues_full,work,lwork,info)
    if(info /= 0) print*, "Hessian diagonalization failed: ", info

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
      !call write_matrix(hess,fulldim,fulldim)
      print*, "Orthonormality of the eigenvectors: "
      !call write_matrix(dot,fulldim,fulldim)
      !call write_matrix(dot,3,3)
      do i=1,10
        print*, dot(i,i)
      end do
    end if

    ev = hess(:,7:fulldim)

    rxyz0 = ats%ats

   !calculate the masses in the new coordinate system
   do i=1,dim
     new_masses(i) = eigenvalues_full(i+6)/eigenvalues_full_m(i+6)
   end do

   if (debug_extended) then
      print*, "rxyz0: "
      !call write_matrix(rxyz0,3,nat)
      print*, "new_masses: "
      !call write_matrix(new_masses,dim,1)
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

   call write_new_masses_and_eigenvalues(eigenvalues_full)

   call write_node_file(eigenvalues_full)

 end subroutine read_hessian

 subroutine calculate_hessian(ats, hess)
    use precision
    use atomicStructure
    use RuNNerInterface
    implicit none
    real(dp), dimension(3,nat) :: fxyz
    real(dp), dimension(fulldim, fulldim) :: hess
    real(dp), dimension(fulldim) :: r0, r
    type(atStruct) :: ats
    real(dp) :: h, e_a, e_c, e_g, e_i
    integer :: i,j

    r0 = reshape(ats%ats,shape(r))

    h = 1.d-4
    do i=1,fulldim
      do j=1,fulldim
        r = r0
        r(i) = r(i) + h
        r(j) = r(j) + h
        ats%ats = reshape(r,shape(ats%ats))
        call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, e_a, fxyz)

        r = r0
        r(i) = r(i) - h
        r(j) = r(j) + h
        ats%ats = reshape(r,shape(ats%ats))
        call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, e_c, fxyz)

        r = r0
        r(i) = r(i) + h
        r(j) = r(j) - h
        ats%ats = reshape(r,shape(ats%ats))
        call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, e_g, fxyz)

        r = r0
        r(i) = r(i) - h
        r(j) = r(j) - h
        ats%ats = reshape(r,shape(ats%ats))
        call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, e_i, fxyz)

        hess(i,j) = (e_a - e_c - e_g - + e_i)/(4*h**2)

      end do
    end do

 end subroutine calculate_hessian

 subroutine calculate_hessian_forces(ats, hess)
    use precision
    use atomicStructure
    use RuNNerInterface
    implicit none
    real(dp), dimension(3,nat) :: f_plus, f_minus
    real(dp), dimension(fulldim, fulldim) :: hess
    real(dp), dimension(fulldim) :: r0, r, f_p, f_m
    type(atStruct) :: ats
    real(dp) :: h, energy
    integer :: i,j

    r0 = reshape(ats%ats,shape(r))

    h = 1.d-4
    do i=1,fulldim
        r = r0
        r(i) = r(i) + h
        ats%ats = reshape(r,shape(ats%ats))
        call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, energy, f_plus)
        r = r0
        r(i) = r(i) - h
        ats%ats = reshape(r,shape(ats%ats))
        call getRuNNerEnergiesAndForces(globalRuNNerHandles(1), ats, energy, f_minus)
        f_p = - reshape(f_plus,shape(f_p))
        f_m = - reshape(f_minus,shape(f_m))
        hess(:,i) = (f_p - f_m)/(2*h)
    end do

 end subroutine calculate_hessian_forces

 subroutine write_new_masses_and_eigenvalues(eigenvalues_full)
   implicit none
   integer :: i, j, ios
   real(dp), dimension(fulldim) :: eigenvalues_full

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
   real(dp), dimension(fulldim) :: eigenvalues_full

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
     real(dp), dimension(n,m) :: a
     integer :: n,m,i,j
     write(*,*)

     do i = 1,n
        write(*,"(9999(G12.4,:,','))") (a(i,j), j = 1,m)
     end do
     write(*,*)
  end subroutine write_matrix

end module parameters
