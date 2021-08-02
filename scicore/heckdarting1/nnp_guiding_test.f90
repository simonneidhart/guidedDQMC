program test
  use precision
  use atomicStructure
  use RuNNerInterface
  use parameters

  implicit real(dp) (a-h,o-z)
  character(len=*), parameter :: RuNNerPath = '/kernph/finjon00/FHMC-periodic/RuNNer-sockets/RuNNer/build/RuNNer.x'
  character(len=*), parameter :: RuNNerParams = '/kernph/finjon00/FHMC-periodic/train-nnp/6/predict-small/*'
  type(atStruct) :: ats
  real(dp), allocatable, dimension(:,:) :: energies, fxyz, rxyz
  real(dp), allocatable, dimension(:) :: r_temp, mu_w, mu

  call as_readAscii('geom/1.ascii', ats, .true.)
  allocate(fxyz(3, ats%nat))
  allocate(rxyz(3, ats%nat))

  call startRuNNer(ats, RuNNerPath, RuNNerParams, globalRuNNerHandle)

  print*, "RuNNer started"

  call getRuNNerEnergiesAndForces(globalRuNNerHandle, ats, ener, fxyz)

  call read_hessian(ats, e0_theory, .FALSE.)

  allocate(r_temp(dim))
  allocate(mu_w(dim))
  allocate(mu(dim))
  mu = 0.d0

  allocate(energies(50,dim))
  energies = 0.d0
  print*, "total dims: ", dim
  do k=1,dim
    mu_w = mu
    ener = 0
    mu_w(k) = mu_w(k)-25*0.1
    do j=1,50
      r_temp = matmul(ev,mu_w)
      ats%ats = reshape(r_temp,shape(rxyz0)) + rxyz0
      call getRuNNerEnergiesAndForces(globalRuNNerHandle, ats, ener, fxyz)
      mu_w(k) = mu_w(k) + 0.1
      energies(j,k) = ener
    end do
    print*, "complete dimension ",k
  end do


  !output final positions of the walkers
  open(unit=13, file="ener", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,dim
    do j=1,50
      write(unit=13, fmt=*) energies(j, i)
    end do
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"

  print*, "all done"

end program
