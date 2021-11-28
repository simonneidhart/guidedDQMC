!evaluates the PES along the normal modes arround the minimal energy position to
!plot and compare to the harmonic approximation.
program guiding_test
  use mod_dftbp
  use dftbp_parameters
  implicit none
  character(len=100) :: dftb_path
  integer :: i, j, k, n_eval, ios
  real(8), allocatable, dimension(:,:) :: energies, fxyz, rxyz
  real(8), allocatable, dimension(:) :: r_temp, mu_w
  real(8) :: e0_theory, mu, ener

  dftb_path = "dftb_temp_folder"
  call setup_dftbp(dftb_path, 0)

  call calculate_guiding_wf(e0_theory, .FALSE., .TRUE.)

  allocate(fxyz(3,nat))
  allocate(rxyz(3,nat))
  allocate(r_temp(dim))
  allocate(mu_w(dim))
  mu = 0.d0
  n_eval = 100 !100 energy evaluations along each vibrational mode
  allocate(energies(n_eval,dim))

  do k=1,dim
    mu_w = mu
    ener = 0.d0
    mu_w(k) = mu_w(k)-n_eval/2*0.1d0/Bohr_Ang !evaluate every 0.1/Bohr_Ang Bohr
    do j=1,n_eval
      r_temp = matmul(ev,mu_w)
      rxyz = reshape(r_temp,shape(rxyz)) + rxyz0
      call energyandforces(nat, rxyz, fxyz, ener, alat, deralat, atomnames, 0)
      mu_w(k) = mu_w(k) + 0.1d0/Bohr_Ang
      energies(j,k) = ener
    end do
    print*, "Dimension ", k, " complete"
  end do


  !output the PES
  open(unit=13, file="pes.out", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,dim
    do j=1,100
      write(unit=13, fmt=*) energies(j, i)
    end do
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"

  print*, "all done"

end program guiding_test
