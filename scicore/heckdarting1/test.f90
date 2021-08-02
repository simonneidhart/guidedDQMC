program test
  use precision
  use atomicStructure
  use RuNNerInterface
  
  implicit none
  print*, "Start"
  call do_all()
  print*, "Done"


end program

subroutine do_all()
  use precision
  use atomicStructure
  use RuNNerInterface

  implicit none

  character(len=*), parameter :: RuNNerPath = '/kernph/finjon00/FHMC-periodic/RuNNer-sockets/RuNNer/build/RuNNer.x'
  character(len=*), parameter :: RuNNerParams = '/kernph/finjon00/FHMC-periodic/train-nnp/6/predict-small/*'
  type(atStruct) :: ats
  real(dp) :: energy
  real(dp), allocatable :: forces(:,:)

  call as_readAscii('1.ascii', ats, .true.)
  allocate(forces(3, ats%nat))

  call startRuNNer(ats, RuNNerPath, RuNNerParams, globalRuNNerHandle)

  print*, "RuNNer started"

  call getRuNNerEnergiesAndForces(globalRuNNerHandle, ats, energy, forces)

  print*, energy
  print*, ats%nat

end subroutine do_all
