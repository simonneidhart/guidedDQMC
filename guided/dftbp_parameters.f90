module dftbp_parameters
  implicit none
  !global parameters
  integer :: nat, dim, fulldim
  real(8), parameter :: Bohr_Ang = 0.5291772109d0 !conversion factor Bohr to Angsroem
  real(8), parameter :: u_me = 1822.888486d0 !conversion factor u to electron mass
  real(8), allocatable, dimension(:) :: sigma2, new_masses
  real(8), allocatable, dimension(:,:) :: ev !eigenvectors of the hessian
  real(8), dimension(3,3) :: alat, deralat
  real(8), allocatable, dimension(:,:) :: rxyz0 !minimal energy position of the atoms
  character(len=2), allocatable, dimension(:) :: atomnames

contains

!calculate the hessian and the guiding wave function for the DQMC calculation.
subroutine calculate_guiding_wf(e0_theory, debug_extended, grad_descent)
  implicit none
  character(len=100) :: line
  real(8), allocatable, dimension(:,:) :: hess, hess_m, dot, eigenvectors_full, rxyz
  real(8), allocatable, dimension(:) :: masses, eigenvalues_full, eigenvalues_full_m
  real(8), allocatable, dimension(:) :: work
  real(8) :: e0_theory, masses_in
  integer :: i, j, ios, lwork, info, nat2
  logical :: debug_extended, grad_descent

  !the masses file contains the number of atoms and all the N_at atom masses
  open(unit=16, file="masses.in", iostat=ios)
  if ( ios /= 0 ) stop "Error opening file unit 16"

  read(16,*) nat
  fulldim = nat*3 !dimensionality of the full 3*nat system
  dim = fulldim - 6 !dimensionality reduced by the rotations and tranlations
  write(*,'(A,I5,A,I5,A,I5)') "Nat = ", nat," 3*Nat = ", fulldim, " 3*Nat - 6 = ", dim

  allocate(rxyz0(3,nat))
  allocate(rxyz(3,nat))
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

  !read the masses.in file
  do i=1,nat
  read(16,*) masses_in
  masses(3*i - 2) = masses_in
  masses(3*i - 1) = masses_in
  masses(3*i) = masses_in
  end do

  close(unit=16, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 16"

  masses = masses*u_me !convert masses from u to electron masses

  !lattice vectors for non periodic systems
  alat(1,1)=1.d0
  alat(2,1)=0.d0
  alat(3,1)=0.d0

  alat(1,2)=0.d0
  alat(2,2)=1.d0
  alat(3,2)=0.d0

  alat(1,3)=0.d0
  alat(2,3)=0.d0
  alat(3,3)=1.d0

  !read the inital atomic positions
  open(unit=12, file="geometry.xyz", iostat=ios, status="old")
  if ( ios /= 0 ) stop "Error opening file "

  read(unit=12, fmt=*) nat2
  if (nat2 .neqv. nat) then
    print*, "Nat in masses.in does not match with geometry.xyz"
    stop
  end if

  read(unit=12, fmt=*) line
  do i = 1, nat
   read(unit=12, fmt=*) atomnames(i), ( rxyz0(j, i), j = 1, 3)
  end do

  !Transform angstroem input to atomic units
  alat=alat/Bohr_Ang
  rxyz0=rxyz0/Bohr_Ang

  close(unit=12, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 16"

  !find the minimum by using gradient descent
  if (grad_descent) then
    call gradient_descent(rxyz0)
  end if

  print*, "rxyz0 (Bohr): "
  call write_matrix(rxyz0,3,nat)

  !calculate the hessian by force derrivatives
  rxyz = rxyz0
  call calculate_hessian_forces(hess, rxyz)

  !reading the hessian from dftb+
  !open(unit=16, file="hessian.out", iostat=ios)
  !if ( ios /= 0 ) stop "Error opening file unit 16"
  !read(16,*) hess
  !close(unit=16, iostat=ios, status="keep")
  !if ( ios /= 0 ) stop "Error closing file unit 16"

  !scale the hessian with the masses
  do i=1,fulldim
    do j=1,fulldim
      hess_m(i,j) = hess(i,j)/sqrt(masses(i)*masses(j))
    end do
  end do

  hess = 0.5d0 * (hess + transpose(hess)) !symmetrize the hessian matrix
  hess_m = 0.5d0 * (hess_m + transpose(hess_m)) !symmetrize the mass-scaled hessian matrix

  !diagonalize the hessian and the mass-scaled hessian
  lwork = max(1,3*fulldim-1)
  allocate(work(lwork))
  call dsyev('V','U',fulldim,hess_m,fulldim,eigenvalues_full_m,work,lwork,info)
  if(info /= 0) print*, "Mass-scaled hessian diagonalization failed: ", info
  call dsyev('V','U',fulldim,hess,fulldim,eigenvalues_full,work,lwork,info)
  if(info /= 0) print*, "Hessian diagonalization failed: ", info
  deallocate(work)

  !we only need the eigenvectors of the non-zero eigenvalues
  ev = hess(:,7:fulldim)

  !calculate the masses in the new coordinate system
  do i=1,dim
    new_masses(i) = eigenvalues_full(i+6)/eigenvalues_full_m(i+6)
  end do

  if (debug_extended) then
    print*, "masses: "
    call write_matrix(masses,fulldim,1)
    print*, "hessian: "
    call write_matrix(hess,fulldim,fulldim)
    print*, "Eigenvalues of the hessian matrix: "
    call write_matrix(eigenvalues_full,fulldim,1)
  end if

  !calculate the width of the gaussian for the trial wave function
  do i=1,dim
    sigma2(i) = 1.d0/(sqrt(eigenvalues_full_m(i+6))*new_masses(i))
  end do

  if (debug_extended) then !write the extended information
    do i=1,fulldim !check if the eigenvectors are orthagonal
      do j=1,fulldim
        dot(i,j) = dot_product(hess(:,i),hess(:,j))
      end do
    end do
    print*, "Eigenvectors of the hessian matrix: "
    call write_matrix(hess,fulldim,fulldim)
    print*, "Orthonormality of the eigenvectors: "
    call write_matrix(dot,fulldim,fulldim)
    print*, "new_masses: "
    call write_matrix(new_masses,dim,1)
    print*, "sigma^2: "
    call write_matrix(sigma2,dim,1)
  end if

  !calculate the theoretical value for e0 based only on the harmonic approximation
  e0_theory = 0
  do i=1,dim
    e0_theory = e0_theory + 0.5*sqrt(eigenvalues_full_m(i+6))
  end do

  call write_mode_file(eigenvalues_full)

  call write_new_masses_and_eigenvalues(eigenvalues_full)

end subroutine calculate_guiding_wf

 subroutine write_new_masses_and_eigenvalues(eigenvalues_full)
   implicit none
   integer :: i, j, ios
   real(8), dimension(fulldim) :: eigenvalues_full

   open(unit=13, file="new_masses.out", iostat=ios, action="write")
   if ( ios /= 0 ) stop "Error opening file unit 13"
   do i=1,dim
     write(unit=13, fmt=*) new_masses(i)
   end do
   close(unit=13, iostat=ios, status="keep")
   if ( ios /= 0 ) stop "Error closing file unit 13"

   open(unit=13, file="eigenvalues.out", iostat=ios, action="write")
   if ( ios /= 0 ) stop "Error opening file unit 13"
   do i=1,dim
     write(unit=13, fmt=*) eigenvalues_full(i+6)
   end do
   close(unit=13, iostat=ios, status="keep")
   if ( ios /= 0 ) stop "Error closing file unit 13"

 end subroutine write_new_masses_and_eigenvalues

 !output mode file for visuatization
 subroutine write_mode_file(eigenvalues_full)
   implicit none
   integer :: i, j, ios
   real(8), dimension(fulldim) :: eigenvalues_full

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

 end subroutine write_mode_file

 !calculate the hessian by the second derrivative of the postions
 subroutine calculate_hessian(hess, rxyz)
    use mod_dftbp
    implicit none
    real(8), dimension(3,nat) :: fxyz, rxyz
    real(8), dimension(fulldim, fulldim) :: hess
    real(8), dimension(fulldim) :: r0, r
    real(8) :: h, e_a, e_c, e_g, e_i
    integer :: i,j

    r0 = reshape(rxyz,shape(r))

    h = 1.d-6
    do i=1,fulldim
      do j=1,fulldim
        r = r0
        r(i) = r(i) + h
        r(j) = r(j) + h
        rxyz = reshape(r,shape(rxyz))
        call energyandforces(nat, rxyz, fxyz, e_a, alat, deralat, atomnames, 0)

        r = r0
        r(i) = r(i) - h
        r(j) = r(j) + h
        rxyz = reshape(r,shape(rxyz))
        call energyandforces(nat, rxyz, fxyz, e_c, alat, deralat, atomnames, 0)

        r = r0
        r(i) = r(i) + h
        r(j) = r(j) - h
        rxyz = reshape(r,shape(rxyz))
        call energyandforces(nat, rxyz, fxyz, e_g, alat, deralat, atomnames, 0)

        r = r0
        r(i) = r(i) - h
        r(j) = r(j) - h
        rxyz = reshape(r,shape(rxyz))
        call energyandforces(nat, rxyz, fxyz, e_i, alat, deralat, atomnames, 0)

        hess(i,j) = (e_a - e_c - e_g - + e_i)/(4*h**2)

      end do
    end do

 end subroutine calculate_hessian

!calculate the hessian by the derrivative of the forces
 subroutine calculate_hessian_forces(hess, rxyz)
    use mod_dftbp
    implicit none
    real(8), dimension(3,nat) :: f_plus, f_minus, rxyz
    real(8), dimension(fulldim, fulldim) :: hess
    real(8), dimension(fulldim) :: r0, r, f_p, f_m
    real(8) :: h, energy
    integer :: i,j

    r0 = reshape(rxyz,shape(r))

    h = 1.d-6
    do i=1,fulldim
        r = r0
        r(i) = r(i) + h
        rxyz = reshape(r,shape(rxyz))
        call energyandforces(nat, rxyz, f_plus, energy, alat, deralat, atomnames, 0)
        r = r0
        r(i) = r(i) - h
        rxyz = reshape(r,shape(rxyz))
        call energyandforces(nat, rxyz, f_minus, energy, alat, deralat, atomnames, 0)
        f_p = - reshape(f_plus,shape(f_p))
        f_m = - reshape(f_minus,shape(f_m))
        hess(:,i) = (f_p - f_m)/(2*h)
    end do

 end subroutine calculate_hessian_forces

 !find the minima of the PES using steepest descent with gradient feedback
 subroutine gradient_descent(rxyz)
   use mod_dftbp
   implicit none
   real(8) :: etot, alpha=1d-3, norm=1.d0, vectornorm, vectornorm2, angle, angle1
   real(8) ,  dimension(3,nat) :: fxyz, prev_fxyz, rxyz
   integer :: i, j, k, ios, iter=0, max_iter = 1000

   call energyandforces(nat, rxyz, fxyz, etot, alat, deralat, atomnames, 0)
   print*, "The energy in the given configuration is: ", etot

   do while ( norm>1d-6 )
     iter = iter+1 !count iterations
     if (iter .gt. max_iter) then
       print*, "Gradient descent did not converge after ", max_iter, " iterations!"
       stop
     end if
     prev_fxyz = fxyz
     rxyz = rxyz + fxyz*alpha !updating the positions
     !calculates negative gradient and total energy
     call energyandforces(nat, rxyz, fxyz, etot, alat, deralat, atomnames, 0)

     !calculating the gradient norm and the angle between the vectors
     angle = 0.d0
     norm = 0.d0
     do i=1,nat
       vectornorm = 0.d0
       vectornorm2 = 0.d0
       do j=1,3
         vectornorm = vectornorm+fxyz(j,i)**2
         vectornorm2 = vectornorm2+prev_fxyz(j,i)**2
       end do
       norm = norm + sqrt(vectornorm)/nat !average over all the gradient norms
       !angle between the current and previous gradient vector of the atom i
       angle=angle+acos(dot_product(prev_fxyz(:,i),fxyz(:,i))/(sqrt(vectornorm)*sqrt(vectornorm2)))
     end do
     if (angle/nat > 1.0472d0) then  !average of all the angles, gradient feedback
       alpha = alpha/2.d0
     else
       alpha = alpha*1.05d0
     endif
   end do

   print*,  "The energy after the minimization is: ", etot
   write(unit=*, fmt='(I5,A)') iter, " iterations were needed to obtain a precision of 1d-6"
   print*, " "

 end subroutine gradient_descent

!write an n x m matrix
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

end module dftbp_parameters
