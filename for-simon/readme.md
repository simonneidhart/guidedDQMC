

## geometrien
Die ascii files mit -c im Namen haben einfach eine andere Unit-cell. Die mit -t have ich verschoben, sodass die kleinen Moleküle ganz in der Unit-cell sind. Ansonsten sind sie gleich. 
Für manche der Phasen ist auch noch eine kleinere Unit-cell möglich.

## code
So etwas sollte auf scicore in etwa funktionieren. 

```fortran
program test
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

  call getRuNNerEnergiesAndForces(globalRuNNerHandle, ats, energy, forces)

end program
```

