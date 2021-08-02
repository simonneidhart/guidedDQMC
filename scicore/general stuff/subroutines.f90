
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Subroutines for Penalty Function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine SymmetryPenalty(nat,alat,rxyz,width_cutoff,nex_cutoff,ixyzmax,rcov,lengthfp,nat_sphere,dpenaldr, penalty, Trace)
!  Calculates penalty from DimensionalityMatrix which gets calulated using the atomic fingerprints from xyz2devaldr
!  Also calculates the Derivate dpenaldr

!SymmetryPenalty needs rxyz in Bohr

  implicit real*8 (a-h,o-z)
  parameter(natx_sphere=80)  ! maximum number of atoms in the sphere
  parameter(ns=1,np=1)
  parameter(num_dlamda=1) !Number of calculated derivates of Evalues from DimensionalityMatrix
  !                                            starting from largest Evalue
  dimension rxyz(3,nat) !atomic input positions (unchange on output)
  dimension alat(3,3) ! three column vectors specifying the cell vectors
  dimension rcov(nat) ! covalent radii neede for fingerprints
  real*8,allocatable ::  rxyz_sphere(:,:),rcov_sphere(:) ! positions and covalent radii of atoms in sphere
  real*8,allocatable ::  amplitude(:),deramplitude(:) ! amplitude modulationf function for Gaussians and its derivative
  real*8,allocatable ::  fp(:) ! eigenvlues of the overlap matric which give the finperprint vector for the atomic environment
  real*8,allocatable ::  dfpdr(:,:,:) ! dfpdr(ixyz,iat_sp,l)  derivative of l-th eigenvalues with respect
!                                       to the component ixyz  of the position of atom iat_sp in the sphere
  real*8,allocatable ::  fpall(:,:)   ! All the fingerprints for all the atoms
  real*8,allocatable ::  dfpdr0(:,:,:) ! Work array (same as dfpdr) but for a number of atoms in the sphere that is less or
!                                           equal to natx_sphere
  dimension dpenaldr(3,nat) ! derivative of penalty function with respect to all atomic positions
  dimension dlamdadr(3,nat,nat) ! dlamdadr(ixyz,iat,lat) is the derivative of Eval lat from DimensionalityMatrix with respect
  !                                                 to all atomic positions
  real*8,allocatable ::  dFPdrAll(:,:,:,:) !dFPdrAll(iat,ixyz,jat,iatsphere)  derivatice of the jat-th eigenvalue of the
!                                             overlaps mareix centered at atom iatsohere, with respect to the component ixyz of atom iat
  real*8,allocatable ::  DimMatrix(:,:) ! dimensionality matrix nat*nat
  real*8,allocatable ::  evalsDimM(:)   ! eigenvalues of dimensionality matrix
  real*8,allocatable ::  dDimMdr(:,:,:,:)   !dDimMdr(iat,ixyz,lat,jat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
  real*8,allocatable ::  WORK(:)
  real*8,allocatable ::  EigenvecDimM(:,:) ! can be eliminated
  real*8,allocatable ::  tmpResdL(:,:,:) !  work array
  integer,allocatable ::  indat(:) !contains  the indices of the atoms in the sphere
  integer,allocatable ::  nat_sphere_array(:) ! needed to calculate maximum number of atoms in sphere
  real*8 :: trace, penalty
  !integer :: num_threads,OMP_GET_NUM_THREADS, omp_get_thread_num, numthread
  integer :: nat_sphere_current_max
  integer :: iat,jat,kat,l,nex_cutoff,nat_sphere,llat
  integer :: i,iiat
 ! double precision :: omp_get_wtime()
  double precision :: omp_t1,omp_t2,omp_t3

  integer :: LDA,LWMAX, INFO, LWORK, devalNr, dwithTrace
  LDA=2*nat
  LWMAX=10000

  !Which Eval should be calculated + with or without d[trace(A)-Lamda_i]/dR
  devalNr=nat !Nach welchem Eigenvalue abgeleitet wird (nat= größter Eval)
  dwithTrace=1  !If dwithTrace=1 => d[trace(A)-Lamda_i]/dR || If dwithTrace=0 => only dLamda/dR



  if (natx_sphere*(ns+np*3).gt.lengthfp) stop 'increase lengthfp'

  allocate(rxyz_sphere(3,natx_sphere),rcov_sphere(natx_sphere))
  allocate(amplitude(natx_sphere),deramplitude(natx_sphere))
  allocate(fp(natx_sphere*(ns+np*3)),dfpdr(3,natx_sphere,natx_sphere*(ns+np*3)))
  allocate(dfpdr0(3,natx_sphere,natx_sphere*(ns+np*3)))
  allocate(dFPdrAll(nat,3,nat,natx_sphere*(ns+np*3)), dDimMdr(nat,3,nat,nat)) !matrixA nat times nat matrix with A_ij being an 3 times 32 object
  allocate(fpall((ns+3*np)*natx_sphere,nat), DimMatrix(nat,nat), evalsDimM(nat))
  allocate(WORK(LWMAX), EigenvecDimM(nat,nat), tmpResdL(nat,3,nat) )
  allocate(indat(natx_sphere), nat_sphere_array(nat))



dDimMdr=0.d0
dFPdrAll=0.d0
fpall=0.d0
dfpdr0=0.d0
dfpdr=0.d0

call cpu_time(t1)
!$ omp_t1=omp_get_wtime()



!Parallize
numthread=0
num_threads=0
!$omp parallel private(lat,llat,indat,rxyz_sphere,&
!$omp                  amplitude, &
!$omp                  deramplitude,rcov_sphere,&
!$omp                  nat_sphere,i,fp,xl,yl,zl,&
!$omp                  dfpdr0,dfpdr,kat,iiat,&
!$omp                  l,numthread)

! num_threads=OMP_GET_NUM_THREADS()
! numthread = omp_get_thread_num()

!write(*,*) num_threads, numthread !if (numthread==0)

!$omp do schedule(static)

!Iterate over all atoms in Cell----------------------------------------------------------------------------------------------------------
do lat=1,nat

! Calculate all the fingerprints and FP derivatives
     call atoms_sphere(width_cutoff,nex_cutoff,lat,llat,ixyzmax,nat,natx_sphere,nat_sphere,alat,rxyz,rxyz_sphere, &
            rcov,rcov_sphere,indat,amplitude,deramplitude)

! makes sure that entries that are not filled by xyz2devaldr are zero
     do i=1,natx_sphere*(ns+3*np)
     fp(i)=0.d0
     enddo

      xl = rxyz(1,lat)
      yl = rxyz(2,lat)
      zl = rxyz(3,lat)
      call xyz2devaldr(nat_sphere,rxyz_sphere,rcov_sphere,amplitude,deramplitude,llat,xl,yl,zl,ns,np,fp,dfpdr0)
      call reformat_devaldr(natx_sphere,nat_sphere,ns,np,dfpdr0,dfpdr)

fpall(:,lat)=fp

!Remap atoms from sphere to Cell
!dFPdrAll(Derivate from FP of atom lat(of ref Cell, direction(x,y,z), derivate due atom iiat, fp l) (dF_i/dr)
 do kat=1,nat_sphere
  iiat=indat(kat)
  do l=1,nat_sphere*(ns+3*np)
   dFPdrAll(lat,1,iiat,l)=dFPdrAll(lat,1,iiat,l) + dfpdr(1,kat,l)
   dFPdrAll(lat,2,iiat,l)=dFPdrAll(lat,2,iiat,l) + dfpdr(2,kat,l)
   dFPdrAll(lat,3,iiat,l)=dFPdrAll(lat,3,iiat,l) + dfpdr(3,kat,l)
  end do
 end do

 !Save nat_sphere of each iteration in array to determine nat_sphere_max
 nat_sphere_array(lat)=nat_sphere

!write(*,*) "nat_sphere", lat, nat_sphere
!write(*,*) "fp:", lat, ( fpall(j, lat), j = 1, 5)

end do !End do über Atome in Cell---------------------------------------------------------------------------
!$omp end do
!$omp end parallel



call cpu_time(t2)
!$ omp_t2=omp_get_wtime()

!write(*,*)"fp_timing:", omp_t2-omp_t1

nat_sphere_current_max=0
do lat=1,nat
 if(nat_sphere_array(lat) > nat_sphere_current_max) nat_sphere_current_max=nat_sphere_array(lat)
end do

!write(*,*) "nat_sphere_current_max",nat_sphere_current_max

!Calc DimensionalityMatrix
DimMatrix=0.d0

 do jat=1,nat
  do iat=1,nat
   DimMatrix(iat,jat)=ddot((ns+3*np)*natx_sphere,fpall(1,iat),1,fpall(1,jat),1)
  end do
 end do


!Calculation of the derivative of the DimensionalityMatrix dDimMdr with respect to ixyz and lat
!call calc_dDimMdr(nat,ns,np,natx_sphere,nat_sphere_current_max,fpall,dFPdrAll,dDimMdr) ! dDimMdr

!$ omp_t4=omp_get_wtime()
call calc_dDimMdr(nat,ns,np,natx_sphere,fpall,dFPdrAll,dDimMdr) ! dDimMdr
!$ omp_t5=omp_get_wtime()



evalsDimM=0.d0
EigenvecDimM=0.d0
EigenvecDimM=DimMatrix

LWORK = -1
!dsyev if "N"=> Eigenvalues only if "V"&INFO=0 Eigenvalues and Eigenvectors
call dsyev('V', 'U', nat, EigenvecDimM, nat, evalsDimM, WORK, LWORK, INFO)
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
call dsyev('V', 'U', nat, EigenvecDimM, nat, evalsDimM, WORK, LWORK, INFO)



trace=0.d0
 do i = 1, nat
  trace = trace + DimMatrix(i,i)
 end do
penalty=(evalsDimM(devalNr))


! Calculates the derivate form Evals from DimensionalityMatrix to dr: dLamda_i/dr=<v_i|dA/dr|v_i>
! belonging to the greatest num_dlamda(is an integer number) Evalues of DimMatrix
! Returns dlamdadr(3,nat,nat)
dlamdadr=0.d0
call calc_multi_dEvalDimMdr(nat,num_dlamda,EigenvecDimM,dDimMdr,dlamdadr)

!In the moment only one Eval derivate of DimMatrix is needed for penalty
dpenaldr=dlamdadr(:,:,devalNr)

!For calculation of d[trace(A)-Lamda_i]/dR = dA_11/dR + dA_22/dR +... +dLamda/dR
if(dwithTrace==1) then
dpenaldr=dpenaldr*(-1.d0)
 do iat=1,nat
  do jat=1,nat
   dpenaldr(1,iat)=dpenaldr(1,iat) + dDimMdr(jat,1,iat,jat)
   dpenaldr(2,iat)=dpenaldr(2,iat) + dDimMdr(jat,2,iat,jat)
   dpenaldr(3,iat)=dpenaldr(3,iat) + dDimMdr(jat,3,iat,jat)
  end do
 end do

!Penalty=Trace - Lamda(1)
 penalty=trace-penalty

end if



call cpu_time(t3)
!$ omp_t3=omp_get_wtime()


! write(*,'(3(A,e10.3))')"OMP: T1: ",omp_t2-omp_t1,"      T2: ",omp_t3-omp_t2, "      Tall ",omp_t3-omp_t1, &
! "      TdDimdr ",omp_t5-omp_t4

!write(*,'(3(A,e10.3))')"t2-t1",t2-t1,"      t3-t2 ",t3-t2, "      t_all ",t3-t1

  deallocate(tmpResdL,EigenvecDimM,dDimMdr)
  deallocate(DimMatrix,fpall,dFPdrAll)



  deallocate(rxyz_sphere,rcov_sphere)
  deallocate(amplitude,deramplitude)
  deallocate(fp,dfpdr,evalsDimM)
  deallocate(dfpdr0)
  deallocate(indat)



end subroutine SymmetryPenalty

subroutine calc_dDimMdr(nat,ns,np,natx_sphere,fpall,tmpResdL,dDimMdr)
!Calculate derivate of DimMatrix with respect to atomic position
!Return result in dDimMdr
implicit real*8 (a-h,o-z)
integer :: nat, iat ,jat ,kat
real*8 :: dDimMdr(nat,3,nat,nat) !dDimMdr(iat,jat,ixyz,lat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
real*8 :: dDimMdrWork(nat,3,nat,nat) !WorkArray

real*8 :: tmpResdL(nat*nat*3,(ns+3*np)*natx_sphere) !Working Array dim(3*nat*nat,(ns+3*np)*natx_sphere)
real*8 :: fpall((ns+3*np)*natx_sphere,nat)


call calc_dDimMdr_MatMul(nat,ns,np,natx_sphere,fpall,tmpResdL,dDimMdr)

dDimMdrWork=dDimMdr
do kat=1,nat
 do jat=1,nat
  do iat=1,nat

     dDimMdr(iat,1,kat,jat) = dDimMdr(iat,1,kat,jat)+dDimMdrWork(jat,1,kat,iat)
     dDimMdr(iat,2,kat,jat) = dDimMdr(iat,2,kat,jat)+dDimMdrWork(jat,2,kat,iat)
     dDimMdr(iat,3,kat,jat) = dDimMdr(iat,3,kat,jat)+dDimMdrWork(jat,3,kat,iat)

    end do
  end do
 end do

end subroutine calc_dDimMdr

subroutine calc_dDimMdr_MatMul(nat,ns,np,natx_sphere,fpall,tmpResdL,Res)!dFPdrAll,dDimMdr -> tmpResdL,Res
!Calculate derivate of DimMatrix with respect to atomic position
!Return result in dDimMdr
implicit real*8 (a-h,o-z)
integer :: nat

real*8 :: tmpResdL(nat*nat*3,(ns+3*np)*natx_sphere) !Working Array dim(3*nat*nat,(ns+3*np)*natx_sphere)
real*8 :: fpall((ns+3*np)*natx_sphere,nat)
real*8 :: Res(nat*nat*3,nat)

double precision :: alpha,beta
INTEGER          M, K, N

!tmpResdL being an work array to remap the 4 dim matrix dFPdrAll into an 2 dim matrix
!for faster computing time while doing an matrix multiplication
!tmpResdL is structured like:

! (dFPdrAll with xyz=x,lat=1)
! (dFPdrAll with xyz=y,lat=1)
! (dFPdrAll with xyz=z,lat=1)
! (dFPdrAll with xyz=x,lat=2)
! (dFPdrAll with xyz=y,lat=2)
! (dFPdrAll with xyz=z,lat=2)
         ! .
         ! .
         ! .
         ! .
! (dFPdrAll with xyz=x,lat=nat)
! (dFPdrAll with xyz=y,lat=nat)
! (dFPdrAll with xyz=z,lat=nat)
!
!With e.g. (dFPdrAll with xyz=x,lat=1) being the derivative of the FPall
!With respect to the x direction and to lat=1.
!(dFPdrAll with xyz=x,lat=1) has the dimensions (nat,natx_sphere*(ns+np*3))
!tmpResdL has therefore the Dimsions (nat*nat*3,natx_sphere*(ns+np*3))


!Parameter declaration for DGEMM MatMul
alpha=1.d0
beta=0.d0
M=nat*nat*3
K=(ns+3*np)*natx_sphere
N=nat

!DGEMM for optimized matrix multiplication
!DGEMM computes real matrix Res=alpha*tmpResdL*fpall+beta*Res
CALL DGEMM('N','N',M,N,K,alpha,tmpResdL,M,fpall,K,beta,Res,M)

! The Res Matrix is structured like:
! (dDimMdr with xyz=x,lat=1)
! (dDimMdr with xyz=y,lat=1)
! (dDimMdr with xyz=z,lat=1)
! (dDimMdr with xyz=x,lat=2)
! (dDimMdr with xyz=y,lat=2)
! (dDimMdr with xyz=z,lat=2)
         ! .
         ! .
         ! .
         ! .
! (dDimMdr with xyz=x,lat=nat)
! (dDimMdr with xyz=y,lat=nat)
! (dDimMdr with xyz=z,lat=nat)
!
!With e.g. (dDimMdr with xyz=x,lat=1) being the derivative of the DimensionalityMatrix
!With respect to the x direction and to lat=1 (with dimesionions (nat,nat))
!Res has therefore the Dimsions (nat*nat*3,nat)


end subroutine calc_dDimMdr_MatMul

subroutine calc_multi_dEvalDimMdr(nat, num_dlamda,EigenvecDimM,dDimMdr, dlamdadr)
!Calculate derivate of Eval from DimMatrix with respect to atomic position
!Return result in devalDimMdr
implicit real*8 (a-h,o-z)
integer :: nat, iat ,num_dlamda, i
real*8 :: EigenvecDimM(nat, nat) !EigenvecDimM(iat,lat) Evec of DimMatrix belonging to Eval lat for dEval/dr
real*8 :: dDimMdr(nat,3,nat,nat)  !dDimMdr(iat,jat,ixyz,lat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
real*8 :: dlamdadr(3,nat, nat) ! dlamdadr(ixyz,iat,lat) is the derivative of Eval lat from DimensionalityMatrix with respect
  !                                            to all atomic positions

if(num_dlamda>nat) stop "num_dlamda>nat"


!$omp parallel private(iat, i)

!$omp do schedule(static)
do i = 1, num_dlamda
!iat=Evalue that we calulate the derivate of
iat=nat-i +1

! Calc derivate form Eval with devalNr from A to dr: dLamda_i/dr=<v_i|dA/dr|v_i>
! Returns dpenaldr
call calc_dEvalDimMdr(nat,EigenvecDimM(1,iat),dDimMdr,dlamdadr(1,1,iat))


end do !End do über Atome in Cell---------------------------------------------------------------------------
!$omp end do
!$omp end parallel

! Old single derivative call in SymmetryPenalty:
! Calc derivate form Eval with devalNr from A to dr: dLamda_i/dr=<v_i|dA/dr|v_i>
! Returns dpenaldr
!call calc_dEvalDimMdr(nat,EigenvecDimM(:,devalNr),dDimMdr,dpenaldr)


end subroutine calc_multi_dEvalDimMdr

subroutine calc_dEvalDimMdr(nat,Evec,dDimMdr, devalDimMdr)
!Calculate derivate of Eval from DimMatrix with respect to atomic position
!Return result in devalDimMdr
implicit real*8 (a-h,o-z)
integer :: nat, iat ,jat ,kat
real*8 :: Evec(nat) !Evec of DimMatrix belonging to Eval for dEval/dr
real*8 :: dDimMdr(nat,3,nat,nat)  !dDimMdr(iat,jat,ixyz,lat) derivative of the matrix element iat,jat, with respect
!                                                to ixyz,lat
real*8 :: tmpResdL(nat,3,nat) !Working Array
real*8 :: devalDimMdr(3,nat) ! devalDimMdr(ixyz,lat) Derivate of Eval from DimMatrix with respect to ixyz,lat


!Calc derivate form Eval from DimMatrix to dr: dLamda_i/dr=<v_i|dDimM/dr|v_i>
tmpResdL=0.d0
 do iat=1,nat
  do jat=1,nat
   do kat=1,nat
   tmpResdL(iat,1,jat)=tmpResdL(iat,1,jat)+Evec(kat)*dDimMdr(kat,1,jat,iat)
   tmpResdL(iat,2,jat)=tmpResdL(iat,2,jat)+Evec(kat)*dDimMdr(kat,2,jat,iat)
   tmpResdL(iat,3,jat)=tmpResdL(iat,3,jat)+Evec(kat)*dDimMdr(kat,3,jat,iat)
   end do
  end do
 end do

devalDimMdr=0.d0
 do iat=1,nat
  do kat=1,nat
   devalDimMdr(1,iat)=devalDimMdr(1,iat) + tmpResdL(kat,1,iat)*Evec(kat)
   devalDimMdr(2,iat)=devalDimMdr(2,iat) + tmpResdL(kat,2,iat)*Evec(kat)
   devalDimMdr(3,iat)=devalDimMdr(3,iat) + tmpResdL(kat,3,iat)*Evec(kat)
  end do
 end do

end subroutine calc_dEvalDimMdr




  subroutine reformat_devaldr(natx_sphere,nat_sphere,ns,np,devaldr0,devaldr)
      implicit real*8 (a-h,o-z)
      dimension devaldr0(3,nat_sphere,nat_sphere*(ns+np*3))
      dimension  devaldr(3,natx_sphere,natx_sphere*(ns+np*3))

      do j=1,nat_sphere*(ns+np*3)
      do i=1,nat_sphere
      devaldr(1,i,j)=devaldr0(1,i,j)
      devaldr(2,i,j)=devaldr0(2,i,j)
      devaldr(3,i,j)=devaldr0(3,i,j)
      enddo
      enddo
  end subroutine reformat_devaldr


subroutine xyz2devaldr(nat,rxyz,rcov,amplitude,deramplitude,lat,xl,yl,zl,ns,np,eval,devaldr_out)
! Calculqates the derivative of all eigenvalues of an atomic fingerprint with respect to the atomic positions
! Nat=NatSphere da wir eval von sphere davor hatten bzw. für mich sollte es
! Eval = FP
  implicit real*8 (a-h,o-z)
  real*8 :: rxyz(3,nat), rcov(nat), eval(nat*(ns+np*3)),devaldr_out(3,nat,nat*(ns+np*3))
  real*8 :: amplitude(nat), deramplitude(nat)
  real*8,allocatable ::  ovrlp(:,:), evecn(:,:),evec(:,:),work(:),devaldr(:,:,:)
  real*8  alpha(nat), cs(10),cp(10)
  integer :: norb

    nt=3*np+ns
    norb= nat*(ns+np*3)
    allocate(ovrlp(norb,norb),evecn(norb,norb),evec(norb,norb),devaldr(norb,3,nat))
    lwork=100*norb
    allocate(work(lwork))

    do iat=1,nat
       alpha(iat)=.5d0/rcov(iat)**2
    enddo
    ! Specify the width of the Gaussians if several Gaussians per l-channel are used
    do i=1,10
      cs(i)=sqrt(2.d0)**(i-1)
      cp(i)=sqrt(2.d0)**(i-1)
    enddo

    call crtovrlp(nat,rxyz,alpha,cs,cp,ns,np,ovrlp)

    do iat=1,nat
        do iorb=1,norb
          devaldr(iorb,1,iat)=0.d0
          devaldr(iorb,2,iat)=0.d0
          devaldr(iorb,3,iat)=0.d0
        enddo
    enddo

     call multamp(nat,ovrlp,amplitude,norb,ns,np,evecn)
     call dsyev('V','L', norb, evecn, norb, eval, work, lwork, info)
     if(info/=0) stop ' ERROR in dsyev'
  ! eigenvalues in decreasing order
     do i=1,norb/2
        t1=eval(i)
        t2=eval(norb-i+1)
        eval(i)=t2
        eval(norb-i+1)=t1
     enddo

!return    ! if no derivatives are needed
        call rots(norb,norb,evecn,evec)

! Now calculate derivatives
  !  <s|s>
  do jat=1,nat
    do js=1,ns
      jorb=(jat-1)*nt+js
      aj=alpha(jat)/cs(js)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat= 1,nat
        do is=1,ns
          iorb=(iat-1)*nt+is
          ai= alpha(iat)/cs(is)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2
          t1=ai*aj
          t2=ai+aj

          ! derivatives
          tt= -2.d0*t1/t2

          xil = xi - xl; yil = yi - yl; zil = zi - zl
          xjl = xj - xl; yjl = yj - yl; zjl = zj - zl

            pijampl=amplitude(iat)*amplitude(jat)
            dipjampl=deramplitude(iat)*amplitude(jat)
            djpiampl=deramplitude(jat)*amplitude(iat)

            deri1 =  pijampl*(tt*ovrlp(iorb,jorb)*xij) + dipjampl*ovrlp(iorb,jorb)*xil
            derj1 = -pijampl*(tt*ovrlp(iorb,jorb)*xij) + djpiampl*ovrlp(iorb,jorb)*xjl
            derl1 = -xil*dipjampl*ovrlp(iorb,jorb) - djpiampl*ovrlp(iorb,jorb)*xjl

            deri2 =  pijampl*(tt*ovrlp(iorb,jorb)*yij) + dipjampl*ovrlp(iorb,jorb)*yil
            derj2 = -pijampl*(tt*ovrlp(iorb,jorb)*yij) + djpiampl*ovrlp(iorb,jorb)*yjl
            derl2 = -yil*dipjampl*ovrlp(iorb,jorb) - djpiampl*ovrlp(iorb,jorb)*yjl

            deri3 =  pijampl*(tt*ovrlp(iorb,jorb)*zij) + dipjampl*ovrlp(iorb,jorb)*zil
            derj3 = -pijampl*(tt*ovrlp(iorb,jorb)*zij) + djpiampl*ovrlp(iorb,jorb)*zjl
            derl3 = -zil*dipjampl*ovrlp(iorb,jorb) - djpiampl*ovrlp(iorb,jorb)*zjl


          do korb=1,norb
            kkorb=norb-korb+1  ! Also deriavtive in decreasing order of eigenvalues
            pevec=evec(korb,iorb)*evec(korb,jorb)

            devaldr(kkorb,1,iat)=devaldr(kkorb,1,iat)+pevec*deri1
            devaldr(kkorb,1,jat)=devaldr(kkorb,1,jat)+pevec*derj1
            devaldr(kkorb,1,lat)=devaldr(kkorb,1,lat)+pevec*derl1

            devaldr(kkorb,2,iat)=devaldr(kkorb,2,iat)+pevec*deri2
            devaldr(kkorb,2,jat)=devaldr(kkorb,2,jat)+pevec*derj2
            devaldr(kkorb,2,lat)=devaldr(kkorb,2,lat)+pevec*derl2

            devaldr(kkorb,3,iat)=devaldr(kkorb,3,iat)+pevec*deri3
            devaldr(kkorb,3,jat)=devaldr(kkorb,3,jat)+pevec*derj3
            devaldr(kkorb,3,lat)=devaldr(kkorb,3,lat)+pevec*derl3
          enddo
        enddo
      enddo
    enddo
  enddo

  if (np.eq.0) goto 1111

  !  <pi|sj>
  do jat=1,nat ! kat, kat ! 1,nat
    do js=1,ns

      jorb=(jat-1)*nt+js
      aj=alpha(jat)/cs(js)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do ip=1,np

          iorb=(iat-1)*nt+ns+ip
          ai= alpha(iat)/cp(ip)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2

          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
          t3=-2.d0*sqrt(ai)*aj/t2

          ! derivatives
          t5=-2.d0*t1/t2

          xil = xi - xl; yil = yi - yl; zil = zi - zl
          xjl = xj - xl; yjl = yj - yl; zjl = zj - zl

            pijampl=amplitude(iat)*amplitude(jat)
            !dipjampl=deramplitude(iat)*amplitude(jat)
            !djpiampl=deramplitude(jat)*amplitude(iat)

            derx00i = pijampl*(ovrlp(iorb,jorb)*t5*xij+t3*sij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
            derx00j = -pijampl*(ovrlp(iorb,jorb)*t5*xij+t3*sij) + &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl
            derx00l = -xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -          &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl

            dery00i = pijampl*(ovrlp(iorb,jorb)*t5*yij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
            dery00j = -pijampl*(ovrlp(iorb,jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl
            dery00l = -yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl

            derz00i = pijampl*(ovrlp(iorb,jorb)*t5*zij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
            derz00j = -pijampl*(ovrlp(iorb,jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl
            derz00l = -zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl

            derx10i = pijampl*(ovrlp(iorb+1,jorb)*t5*xij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
            derx10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*xjl
            derx10l = -xil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*xjl

            dery10i = pijampl*(ovrlp(iorb+1,jorb)*t5*yij+t3*sij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
            dery10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*yij+t3*sij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*yjl
            dery10l = -yil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -          &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*yjl

            derz10i = pijampl*(ovrlp(iorb+1,jorb)*t5*zij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
            derz10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*zjl
            derz10l = -zil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*zjl

            derx20i = pijampl*(ovrlp(iorb+2,jorb)*t5*xij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
            derx20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*xjl
            derx20l = -xil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*xjl

            dery20i = pijampl*(ovrlp(iorb+2,jorb)*t5*yij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
            dery20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*yjl
            dery20l = -yil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*yjl

            derz20i = pijampl*(ovrlp(iorb+2,jorb)*t5*zij+t3*sij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
            derz20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*zij+t3*sij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*zjl
            derz20l = -zil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -          &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*zjl

          do korb=1,norb
            kkorb=norb-korb+1  ! Also derivative in decreasing order of eigenvalues
            pevec00=evec(korb,iorb+0)*evec(korb,jorb)
            pevec10=evec(korb,iorb+1)*evec(korb,jorb)
            pevec20=evec(korb,iorb+2)*evec(korb,jorb)

            devaldr(kkorb,1,iat)=devaldr(kkorb,1,iat)+pevec00*derx00i+pevec10*derx10i+pevec20*derx20i
            devaldr(kkorb,1,jat)=devaldr(kkorb,1,jat)+pevec00*derx00j+pevec10*derx10j+pevec20*derx20j
            devaldr(kkorb,1,lat)=devaldr(kkorb,1,lat)+pevec00*derx00l+pevec10*derx10l+pevec20*derx20l

            devaldr(kkorb,2,iat)=devaldr(kkorb,2,iat)+pevec00*dery00i+pevec10*dery10i+pevec20*dery20i
            devaldr(kkorb,2,jat)=devaldr(kkorb,2,jat)+pevec00*dery00j+pevec10*dery10j+pevec20*dery20j
            devaldr(kkorb,2,lat)=devaldr(kkorb,2,lat)+pevec00*dery00l+pevec10*dery10l+pevec20*dery20l

            devaldr(kkorb,3,iat)=devaldr(kkorb,3,iat)+pevec00*derz00i+pevec10*derz10i+pevec20*derz20i
            devaldr(kkorb,3,jat)=devaldr(kkorb,3,jat)+pevec00*derz00j+pevec10*derz10j+pevec20*derz20j
            devaldr(kkorb,3,lat)=devaldr(kkorb,3,lat)+pevec00*derz00l+pevec10*derz10l+pevec20*derz20l

          enddo
        enddo
      enddo
    enddo
  enddo


  !  <si|pj>
  do jat=1,nat
    do jp=1,np

      jorb=(jat-1)*nt+ns+jp
      aj=alpha(jat)/cp(jp)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do is=1,ns
          iorb=(iat-1)*nt+is
          ai= alpha(iat)/cs(is)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2

          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
          t3=+2.d0*sqrt(aj)*ai/t2

          ! derivatives
          !tt= -2.d0*t1/t2 * sij
          t5=-2.d0*t1/t2

          xil = xi - xl; yil = yi - yl; zil = zi - zl
          xjl = xj - xl; yjl = yj - yl; zjl = zj - zl

            pijampl=amplitude(iat)*amplitude(jat)
            !dipjampl=deramplitude(iat)*amplitude(jat)
            !djpiampl=deramplitude(jat)*amplitude(iat)

            derx00i = pijampl*(ovrlp(iorb,jorb)*t5*xij+t3*sij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
            derx00j = -pijampl*(ovrlp(iorb,jorb)*t5*xij+t3*sij) + &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl
            derx00l = -xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -          &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl

            dery00i = pijampl*(ovrlp(iorb,jorb)*t5*yij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
            dery00j = -pijampl*(ovrlp(iorb,jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl
            dery00l = -yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl

            derz00i = pijampl*(ovrlp(iorb,jorb)*t5*zij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
            derz00j = -pijampl*(ovrlp(iorb,jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl
            derz00l = -zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl

            derx01i = pijampl*(ovrlp(iorb,jorb+1)*t5*xij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
            derx01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*xjl
            derx01l = -xil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*xjl

            dery01i = pijampl*(ovrlp(iorb,jorb+1)*t5*yij+t3*sij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
            dery01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*yij+t3*sij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*yjl
            dery01l = -yil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -          &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*yjl

            derz01i = pijampl*(ovrlp(iorb,jorb+1)*t5*zij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
            derz01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*zjl
            derz01l = -zil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*zjl

            derx02i = pijampl*(ovrlp(iorb,jorb+2)*t5*xij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
            derx02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*xjl
            derx02l = -xil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*xjl

            dery02i = pijampl*(ovrlp(iorb,jorb+2)*t5*yij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
            dery02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*yjl
            dery02l = -yil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*yjl

            derz02i = pijampl*(ovrlp(iorb,jorb+2)*t5*zij+t3*sij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
            derz02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*zij+t3*sij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*zjl
            derz02l = -zil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -          &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*zjl

          do korb=1,norb
            kkorb=norb-korb+1  ! Also derivative in decreasing order of eigenvalues
            pevec00=evec(korb,iorb)*evec(korb,jorb+0)
            pevec01=evec(korb,iorb)*evec(korb,jorb+1)
            pevec02=evec(korb,iorb)*evec(korb,jorb+2)

            devaldr(kkorb,1,iat)=devaldr(kkorb,1,iat)+pevec00*derx00i+pevec01*derx01i+pevec02*derx02i
            devaldr(kkorb,1,jat)=devaldr(kkorb,1,jat)+pevec00*derx00j+pevec01*derx01j+pevec02*derx02j
            devaldr(kkorb,1,lat)=devaldr(kkorb,1,lat)+pevec00*derx00l+pevec01*derx01l+pevec02*derx02l

            devaldr(kkorb,2,iat)=devaldr(kkorb,2,iat)+pevec00*dery00i+pevec01*dery01i+pevec02*dery02i
            devaldr(kkorb,2,jat)=devaldr(kkorb,2,jat)+pevec00*dery00j+pevec01*dery01j+pevec02*dery02j
            devaldr(kkorb,2,lat)=devaldr(kkorb,2,lat)+pevec00*dery00l+pevec01*dery01l+pevec02*dery02l

            devaldr(kkorb,3,iat)=devaldr(kkorb,3,iat)+pevec00*derz00i+pevec01*derz01i+pevec02*derz02i
            devaldr(kkorb,3,jat)=devaldr(kkorb,3,jat)+pevec00*derz00j+pevec01*derz01j+pevec02*derz02j
            devaldr(kkorb,3,lat)=devaldr(kkorb,3,lat)+pevec00*derz00l+pevec01*derz01l+pevec02*derz02l
          enddo

        enddo
      enddo
    enddo
  enddo


  !  <p|p>
  do jat=1,nat
    do jp=1,np

      jorb=(jat-1)*nt+ns+jp
      aj=alpha(jat)/cp(jp)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do ip=1,np
          iorb=(iat-1)*nt+ns+ip
          ai= alpha(iat)/cp(ip)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2
          t1=ai*aj
          t2=ai+aj

          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
          t4= 2.d0*sqrt(t1)/t2
          t5=-2.d0*t1/t2

          ! derivatives

          xil = xi - xl; yil = yi - yl; zil = zi - zl
          xjl = xj - xl; yjl = yj - yl; zjl = zj - zl
            pijampl=amplitude(iat)*amplitude(jat)
            !dipjampl=deramplitude(iat)*amplitude(jat)
            !djpiampl=deramplitude(jat)*amplitude(iat)

            derx00i = pijampl*(ovrlp(iorb,jorb)*t5*xij+t4*t5*sij*xij*2.d0) +  &
                   xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
            derx00j = -pijampl*(ovrlp(iorb,jorb)*t5*xij+t4*t5*sij*xij*2.d0) + &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl
            derx00l = -xil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -                      &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*xjl

            dery00i = pijampl*(ovrlp(iorb,jorb)*t5*yij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
            dery00j = -pijampl*(ovrlp(iorb,jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl
            dery00l = -yil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -      &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*yjl

            derz00i = pijampl*(ovrlp(iorb,jorb)*t5*zij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat)
            derz00j = -pijampl*(ovrlp(iorb,jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl
            derz00l = -zil*deramplitude(iat)*ovrlp(iorb,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb)*deramplitude(jat)*zjl

            derx10i = pijampl*(ovrlp(iorb+1,jorb)*t5*xij+t4*t5*sij*yij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
            derx10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*xij+t4*t5*sij*yij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*xjl
            derx10l = -xil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -                   &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*xjl

            dery10i = pijampl*(ovrlp(iorb+1,jorb)*t5*yij+t4*t5*sij*xij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
            dery10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*yij+t4*t5*sij*xij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*yjl
            dery10l = -yil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*yjl

            derz10i = pijampl*(ovrlp(iorb+1,jorb)*t5*zij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat)
            derz10j = -pijampl*(ovrlp(iorb+1,jorb)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*zjl
            derz10l = -zil*deramplitude(iat)*ovrlp(iorb+1,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+1,jorb)*deramplitude(jat)*zjl

            derx20i = pijampl*(ovrlp(iorb+2,jorb)*t5*xij+t4*t5*sij*zij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
            derx20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*xij+t4*t5*sij*zij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*xjl
            derx20l = -xil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*xjl

            dery20i = pijampl*(ovrlp(iorb+2,jorb)*t5*yij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
            dery20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*yjl
            dery20l = -yil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*yjl

            derz20i = pijampl*(ovrlp(iorb+2,jorb)*t5*zij+t4*t5*sij*xij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat)
            derz20j = -pijampl*(ovrlp(iorb+2,jorb)*t5*zij+t4*t5*sij*xij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*zjl
            derz20l = -zil*deramplitude(iat)*ovrlp(iorb+2,jorb)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb+2,jorb)*deramplitude(jat)*zjl

            derx01i = pijampl*(ovrlp(iorb,jorb+1)*t5*xij+t4*t5*sij*yij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
            derx01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*xij+t4*t5*sij*yij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*xjl
            derx01l = -xil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*xjl

            dery01i = pijampl*(ovrlp(iorb,jorb+1)*t5*yij+t4*t5*sij*xij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
            dery01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*yij+t4*t5*sij*xij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*yjl
            dery01l = -yil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*yjl

            derz01i = pijampl*(ovrlp(iorb,jorb+1)*t5*zij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat)
            derz01j = -pijampl*(ovrlp(iorb,jorb+1)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*zjl
            derz01l = -zil*deramplitude(iat)*ovrlp(iorb,jorb+1)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb+1)*deramplitude(jat)*zjl

            derx11i = pijampl*(ovrlp(iorb+1,jorb+1)*t5*xij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat)
            derx11j = -pijampl*(ovrlp(iorb+1,jorb+1)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*xjl
            derx11l = -xil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*xjl

            dery11i = pijampl*(ovrlp(iorb+1,jorb+1)*t5*yij+t4*t5*sij*yij*2.d0) +  &
                   yil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat)
            dery11j = -pijampl*(ovrlp(iorb+1,jorb+1)*t5*yij+t4*t5*sij*yij*2.d0) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*yjl
            dery11l = -yil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat) -                          &
                    amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*yjl

            derz11i = pijampl*(ovrlp(iorb+1,jorb+1)*t5*zij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat)
            derz11j = -pijampl*(ovrlp(iorb+1,jorb+1)*t5*zij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*zjl
            derz11l = -zil*deramplitude(iat)*ovrlp(iorb+1,jorb+1)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+1,jorb+1)*deramplitude(jat)*zjl

            derx21i = pijampl*(ovrlp(iorb+2,jorb+1)*t5*xij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat)
            derx21j = -pijampl*(ovrlp(iorb+2,jorb+1)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*xjl
            derx21l = -xil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*xjl

            dery21i = pijampl*(ovrlp(iorb+2,jorb+1)*t5*yij+t4*t5*sij*zij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat)
            dery21j = -pijampl*(ovrlp(iorb+2,jorb+1)*t5*yij+t4*t5*sij*zij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*yjl
            dery21l = -yil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat)-                  &
                    amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*yjl

            derz21i = pijampl*(ovrlp(iorb+2,jorb+1)*t5*zij+t4*t5*sij*yij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat)
            derz21j = -pijampl*(ovrlp(iorb+2,jorb+1)*t5*zij+t4*t5*sij*yij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*zjl
            derz21l = -zil*deramplitude(iat)*ovrlp(iorb+2,jorb+1)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb+2,jorb+1)*deramplitude(jat)*zjl

            derx02i = pijampl*(ovrlp(iorb,jorb+2)*t5*xij+t4*t5*sij*zij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
            derx02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*xij+t4*t5*sij*zij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*xjl
            derx02l = -xil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*xjl

            dery02i = pijampl*(ovrlp(iorb,jorb+2)*t5*yij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
            dery02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*yjl
            dery02l = -yil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*yjl

            derz02i = pijampl*(ovrlp(iorb,jorb+2)*t5*zij+t4*t5*sij*xij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat)
            derz02j = -pijampl*(ovrlp(iorb,jorb+2)*t5*zij+t4*t5*sij*xij) + &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*zjl
            derz02l = -zil*deramplitude(iat)*ovrlp(iorb,jorb+2)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb,jorb+2)*deramplitude(jat)*zjl

            derx12i = pijampl*(ovrlp(iorb+1,jorb+2)*t5*xij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat)
            derx12j = -pijampl*(ovrlp(iorb+1,jorb+2)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*xjl
            derx12l = -xil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*xjl

            dery12i = pijampl*(ovrlp(iorb+1,jorb+2)*t5*yij+t4*t5*sij*zij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat)
            dery12j = -pijampl*(ovrlp(iorb+1,jorb+2)*t5*yij+t4*t5*sij*zij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*yjl
            dery12l = -yil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*yjl

            derz12i = pijampl*(ovrlp(iorb+1,jorb+2)*t5*zij+t4*t5*sij*yij) +  &
                   zil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat)
            derz12j = -pijampl*(ovrlp(iorb+1,jorb+2)*t5*zij+t4*t5*sij*yij) + &
                    amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*zjl
            derz12l = -zil*deramplitude(iat)*ovrlp(iorb+1,jorb+2)*amplitude(jat) -                 &
                    amplitude(iat)*ovrlp(iorb+1,jorb+2)*deramplitude(jat)*zjl

            derx22i = pijampl*(ovrlp(iorb+2,jorb+2)*t5*xij) +  &
                   xil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat)
            derx22j = -pijampl*(ovrlp(iorb+2,jorb+2)*t5*xij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*xjl
            derx22l = -xil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*xjl

            dery22i = pijampl*(ovrlp(iorb+2,jorb+2)*t5*yij) +  &
                   yil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat)
            dery22j = -pijampl*(ovrlp(iorb+2,jorb+2)*t5*yij) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*yjl
            dery22l = -yil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat) -   &
                    amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*yjl

            derz22i = pijampl*(ovrlp(iorb+2,jorb+2)*t5*zij+t4*t5*sij*zij*2.d0) +  &
                   zil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat)
            derz22j = -pijampl*(ovrlp(iorb+2,jorb+2)*t5*zij+t4*t5*sij*zij*2.d0) + &
                    amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*zjl
            derz22l = -zil*deramplitude(iat)*ovrlp(iorb+2,jorb+2)*amplitude(jat) -                      &
                    amplitude(iat)*ovrlp(iorb+2,jorb+2)*deramplitude(jat)*zjl

          do korb=1,norb
            kkorb=norb-korb+1  ! Also derivative in decreasing order of eigenvalues
            pevec00=evec(korb,iorb+0)*evec(korb,jorb+0)
            pevec10=evec(korb,iorb+1)*evec(korb,jorb+0)
            pevec20=evec(korb,iorb+2)*evec(korb,jorb+0)

            pevec01=evec(korb,iorb+0)*evec(korb,jorb+1)
            pevec11=evec(korb,iorb+1)*evec(korb,jorb+1)
            pevec21=evec(korb,iorb+2)*evec(korb,jorb+1)

            pevec02=evec(korb,iorb+0)*evec(korb,jorb+2)
            pevec12=evec(korb,iorb+1)*evec(korb,jorb+2)
            pevec22=evec(korb,iorb+2)*evec(korb,jorb+2)

            devaldr(kkorb,1,iat)=devaldr(kkorb,1,iat)+pevec00*derx00i+pevec10*derx10i+pevec20*derx20i &
                                                     +pevec01*derx01i+pevec11*derx11i+pevec21*derx21i &
                                                     +pevec02*derx02i+pevec12*derx12i+pevec22*derx22i
            devaldr(kkorb,1,jat)=devaldr(kkorb,1,jat)+pevec00*derx00j+pevec10*derx10j+pevec20*derx20j &
                                                     +pevec01*derx01j+pevec11*derx11j+pevec21*derx21j &
                                                     +pevec02*derx02j+pevec12*derx12j+pevec22*derx22j
            devaldr(kkorb,1,lat)=devaldr(kkorb,1,lat)+pevec00*derx00l+pevec10*derx10l+pevec20*derx20l &
                                                     +pevec01*derx01l+pevec11*derx11l+pevec21*derx21l &
                                                     +pevec02*derx02l+pevec12*derx12l+pevec22*derx22l

            devaldr(kkorb,2,iat)=devaldr(kkorb,2,iat)+pevec00*dery00i+pevec10*dery10i+pevec20*dery20i &
                                                     +pevec01*dery01i+pevec11*dery11i+pevec21*dery21i &
                                                     +pevec02*dery02i+pevec12*dery12i+pevec22*dery22i
            devaldr(kkorb,2,jat)=devaldr(kkorb,2,jat)+pevec00*dery00j+pevec10*dery10j+pevec20*dery20j &
                                                     +pevec01*dery01j+pevec11*dery11j+pevec21*dery21j &
                                                     +pevec02*dery02j+pevec12*dery12j+pevec22*dery22j
            devaldr(kkorb,2,lat)=devaldr(kkorb,2,lat)+pevec00*dery00l+pevec10*dery10l+pevec20*dery20l &
                                                     +pevec01*dery01l+pevec11*dery11l+pevec21*dery21l &
                                                     +pevec02*dery02l+pevec12*dery12l+pevec22*dery22l

            devaldr(kkorb,3,iat)=devaldr(kkorb,3,iat)+pevec00*derz00i+pevec10*derz10i+pevec20*derz20i &
                                                     +pevec01*derz01i+pevec11*derz11i+pevec21*derz21i &
                                                     +pevec02*derz02i+pevec12*derz12i+pevec22*derz22i
            devaldr(kkorb,3,jat)=devaldr(kkorb,3,jat)+pevec00*derz00j+pevec10*derz10j+pevec20*derz20j &
                                                     +pevec01*derz01j+pevec11*derz11j+pevec21*derz21j &
                                                     +pevec02*derz02j+pevec12*derz12j+pevec22*derz22j
            devaldr(kkorb,3,lat)=devaldr(kkorb,3,lat)+pevec00*derz00l+pevec10*derz10l+pevec20*derz20l &
                                                     +pevec01*derz01l+pevec11*derz11l+pevec21*derz21l &
                                                     +pevec02*derz02l+pevec12*derz12l+pevec22*derz22l
          enddo
        enddo
      enddo
    enddo
  enddo


1111 continue
        call rots(norb,3*nat,devaldr(1,1,1),devaldr_out(1,1,1))
!         call copy(norb*3*nat,devaldr(1,1,1),devaldr_out(1,1,1))
  deallocate(evecn,evec,devaldr)
  deallocate(ovrlp)
  deallocate(work)

end subroutine xyz2devaldr

        subroutine copy(n,a,b)
        implicit real*8 (a-h,o-z)
        dimension a(n),b(n)

        do i=1,n
          b(i)=a(i)
        enddo

        return
        end




        subroutine rots(n1,n2,a,b)
        implicit real*8 (a-h,o-z)
        parameter(lot=32)
        dimension a(n1,n2),b(n2,n1)

        do jj=1,n2,lot
        do i=1,n1
        do j=jj,min(n2,jj+(lot-1))
          b(j,i)=a(i,j)
        enddo
        enddo
        enddo

        return
        end




        subroutine rot(n1,n2,a,b)
        implicit real*8 (a-h,o-z)
        dimension a(n1,n2),b(n2,n1)

        do i=1,n1
        do j=1,n2
          b(j,i)=a(i,j)
        enddo
        enddo

        return
        end



        subroutine rotb(n,a,b)
        implicit real*8 (a-h,o-z)
        parameter(lot=16)
        dimension a(n,n),b(n,n)

! loop over blocks
        do ii=1,n,lot
        do jj=1,n,lot
! loop over elements in each block
        do i=ii,min(n,ii+(lot-1))
        do j=jj,min(n,jj+(lot-1))
          b(j,i)=a(i,j)
        enddo
        enddo
        enddo
        enddo

        return
        end


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Subroutines for FingerPrint !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine atoms_sphere(width_cutoff,nex_cutoff,lat,llat,ixyzmax,nat,natx_sphere,nat_sphere,alat,rxyz,rxyz_sphere, &
                        rcov,rcov_sphere,indat,amplitude,deramplitude)
  implicit real*8 (a-h,o-z)
  dimension rxyz_sphere(3, natx_sphere),rcov_sphere(natx_sphere)
  dimension amplitude(natx_sphere),deramplitude(natx_sphere)
  dimension rxyz(3,nat),rcov(nat)
  dimension alat(3, 3),indat(natx_sphere)
  integer :: nat_sphere


  radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff
  radius_cutoff2=radius_cutoff**2
  factor_cutoff=1.d0/(2.d0*nex_cutoff*width_cutoff**2)
 ! write(*,*) 'width_cutoff,radius_cutoff',width_cutoff,radius_cutoff


     nat_sphere=0
     do jat = 1, nat
         do ix = -ixyzmax,ixyzmax
           do iy = -ixyzmax,ixyzmax
             do iz = -ixyzmax,ixyzmax
                xj = rxyz(1, jat) + ix*alat(1,1)+iy*alat(1,2)+iz*alat(1,3)
                yj = rxyz(2, jat) + ix*alat(2,1)+iy*alat(2,2)+iz*alat(2,3)
                zj = rxyz(3, jat) + ix*alat(3,1)+iy*alat(3,2)+iz*alat(3,3)
                dist2 = (xj-rxyz(1, lat))**2+(yj-rxyz(2, lat))**2+(zj-rxyz(3, lat))**2
               ! write(*,*) xj,rxyz(1, lat),yj,rxyz(2, lat),zj,rxyz(3,lat)


                if (dist2.le.radius_cutoff2) then
                    nat_sphere=nat_sphere+1
                    if (nat_sphere.gt.natx_sphere) stop 'enlarge natx_sphere'
                    !amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
                    temp = (1.d0 - dist2*factor_cutoff)**(nex_cutoff-1)
                    amplitude(nat_sphere) = temp*(1.d0 - dist2*factor_cutoff)
                    deramplitude(nat_sphere) = -2.d0*factor_cutoff*nex_cutoff*temp

                    rxyz_sphere(1,nat_sphere)=xj
                    rxyz_sphere(2,nat_sphere)=yj
                    rxyz_sphere(3,nat_sphere)=zj
                    rcov_sphere(nat_sphere)=rcov(jat)
                    indat(nat_sphere)=jat
                    if (dist2.eq.0.d0) llat=nat_sphere
                endif
             enddo
           enddo
         enddo
     enddo
end subroutine atoms_sphere


subroutine crtovrlp(nat,rxyz,alpha,cs,cp,ns,np,ovrlp)
  implicit real*8 (a-h,o-z)
  real*8  rxyz(3,nat)
  real*8 ovrlp(nat*(ns+3*np),nat*(ns+3*np))
  real*8  alpha(nat), cs(10),cp(10)


  if(ns>10 .or. np > 10) stop 'ns > 10   .or.  np > 10  !'


 ! 1- setup the overlap matrix

  !  <s|s>
  do jat=1,nat
    do js=1,ns
      jorb=(jat-1)*(ns+3*np)+js
      aj=alpha(jat)/cs(js)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do is=1,ns
          !!iorb=iat+(is-1)*nat
          iorb=(iat-1)*(ns+3*np)+is
          ai= alpha(iat)/cs(is)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2
          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
          ovrlp(iorb,jorb)=sij

        enddo
      enddo
    enddo
  enddo


  !  <pi|sj>
  do jat=1,nat
    do js=1,ns

      jorb=(jat-1)*(ns+3*np)+js
      aj=alpha(jat)/cs(js)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do ip=1,np
          !!iorb=1+(iat-1)*3+ns*nat + (ip-1)*3*nat
          iorb=(iat-1)*(ns+3*np)+ns+ip
          ai= alpha(iat)/cp(ip)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2

          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
          t3=-2.d0*sqrt(ai)*aj/t2
          ovrlp(iorb  ,jorb  )= t3 * xij *sij
          ovrlp(iorb+1,jorb  )= t3 * yij *sij
          ovrlp(iorb+2,jorb  )= t3 * zij *sij

        enddo
      enddo
    enddo
  enddo


  !  <si|pj>
  do jat=1,nat
    do jp=1,np

      jorb=(jat-1)*(ns+3*np)+ns+jp
      aj=alpha(jat)/cp(jp)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do is=1,ns
          !!iorb=iat+(is-1)*nat
          iorb=(iat-1)*(ns+3*np)+is
          ai= alpha(iat)/cs(is)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2

          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
          t3=+2.d0*sqrt(aj)*ai/t2
          ovrlp(iorb,jorb  )= t3 * xij *sij
          ovrlp(iorb,jorb+1)= t3 * yij *sij
          ovrlp(iorb,jorb+2)= t3 * zij *sij

        enddo
      enddo
    enddo
  enddo


  !  <p|p>
  do jat=1,nat
    do jp=1,np

      jorb=(jat-1)*(ns+3*np)+ns+jp
      aj=alpha(jat)/cp(jp)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do ip=1,np
          iorb=(iat-1)*(ns+3*np)+ns+ip
          ai= alpha(iat)/cp(ip)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2
          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2)
          t4= 2.d0*sqrt(t1)/t2
          t5=-2.d0*t1/t2

          ovrlp(iorb  ,jorb  )= t4 *(1.d0 + t5* xij* xij)  * sij
          ovrlp(iorb+1,jorb  )= t4 *(       t5* yij* xij)  * sij
          ovrlp(iorb+2,jorb  )= t4 *(       t5* zij* xij)  * sij
          ovrlp(iorb  ,jorb+1)= t4 *(       t5* xij* yij)  * sij
          ovrlp(iorb+1,jorb+1)= t4 *(1.d0+  t5* yij* yij)  * sij
          ovrlp(iorb+2,jorb+1)= t4 *(       t5* zij* yij)  * sij
          ovrlp(iorb  ,jorb+2)= t4 *(       t5* xij* zij)  * sij
          ovrlp(iorb+1,jorb+2)= t4 *(       t5* yij* zij)  * sij
          ovrlp(iorb+2,jorb+2)= t4 *(1.d0+  t5* zij* zij)  * sij

        enddo
      enddo
    enddo
  enddo

end subroutine crtovrlp

subroutine multampspd(nat,ovrlp,amplitude,norb,ns,np,nd,ovrla)
  implicit real*8 (a-h,o-z)
  real*8 amplitude(nat), ovrlp(norb,norb),ovrla(norb,norb)

  do jat = 1,nat
    do j = 1,(ns+3*np+5*nd)
      jorb = (jat -1)*(ns+3*np+5*nd) + j
      do iat = 1,nat
        do i = 1,(ns+3*np+5*nd)
          iorb = (iat-1)*(ns+3*np+5*nd) +i
          ovrla(iorb,jorb) = ovrlp(iorb,jorb)*amplitude(iat)*amplitude(jat)
        end do
      end do
    end do
  end do

end subroutine multampspd

subroutine multampoff(nat,ovrlp,amplitude,norb,ns,np,ovrla)
  implicit real*8 (a-h,o-z)
  real*8 amplitude(nat), ovrlp(norb,norb),ovrla(norb,norb)

  do jat = 1,nat
    do j = 1,(ns+3*np)
      jorb = (jat -1)*(ns+3*np) + j
      do iat = 1,nat
        do i = 1,(ns+3*np)
          iorb = (iat-1)*(ns+3*np) +i
          if (iat.eq.jat) then
          ovrla(iorb,jorb) = ovrlp(iorb,jorb)
          else
          ovrla(iorb,jorb) = ovrlp(iorb,jorb)*amplitude(iat)*amplitude(jat)
          endif
        end do
      end do
    end do
  end do

end subroutine multampoff


subroutine multamp(nat,ovrlp,amplitude,norb,ns,np,ovrla)
  implicit real*8 (a-h,o-z)
  real*8 amplitude(nat), ovrlp(norb,norb),ovrla(norb,norb)

  do jat = 1,nat
    do j = 1,(ns+3*np)
      jorb = (jat -1)*(ns+3*np) + j
      do iat = 1,nat
        do i = 1,(ns+3*np)
          iorb = (iat-1)*(ns+3*np) +i
          ovrla(iorb,jorb) = ovrlp(iorb,jorb)*amplitude(iat)*amplitude(jat)
        end do
      end do
    end do
  end do

end subroutine multamp



!
!      subroutine frac2cart(nat, alat, xyzred, rxyz)
!      implicit real*8 (a-h,o-z)
!      dimension alat(3,3), xyzred(3,nat), rxyz(3,nat)
!
!      do iat=1,nat
!         do i = 1, 3
!            t = 0.d0
!            do j = 1, 3
!               t = t + xyzred(j,iat) * alat(i, j)
!            end do
!            rxyz(i,iat) = t
!         end do
!      enddo
!
!     !  do j=1,3
!     !  do i=1,3
!     !  alat(i,j)=alat(i,j)
!     !  enddo
!     !  enddo
!
!      end subroutine frac2cart
!
!
!
!
! subroutine cart2frac(nat,alat,rxyz,rxyzred)
!   !This subrouine will convert the redernal coordinates into cartesian coordinates
!   implicit real*8 (a-h,o-z)
!   dimension rxyzred(3,nat),rxyz(3,nat),alat(3,3),alatinv(3,3)
!
!     div=(alat(1,1)*alat(2,2)*alat(3,3)-alat(1,1)*alat(2,3)*alat(3,2)- &
!          alat(1,2)*alat(2,1)*alat(3,3)+alat(1,2)*alat(2,3)*alat(3,1)+ &
!          alat(1,3)*alat(2,1)*alat(3,2)-alat(1,3)*alat(2,2)*alat(3,1))
!     div=1.d0/div
!       alatinv(1,1) = (alat(2,2)*alat(3,3)-alat(2,3)*alat(3,2))*div
!       alatinv(1,2) =-(alat(1,2)*alat(3,3)-alat(1,3)*alat(3,2))*div
!       alatinv(1,3) = (alat(1,2)*alat(2,3)-alat(1,3)*alat(2,2))*div
!       alatinv(2,1) =-(alat(2,1)*alat(3,3)-alat(2,3)*alat(3,1))*div
!       alatinv(2,2) = (alat(1,1)*alat(3,3)-alat(1,3)*alat(3,1))*div
!       alatinv(2,3) =-(alat(1,1)*alat(2,3)-alat(1,3)*alat(2,1))*div
!       alatinv(3,1) = (alat(2,1)*alat(3,2)-alat(2,2)*alat(3,1))*div
!       alatinv(3,2) =-(alat(1,1)*alat(3,2)-alat(1,2)*alat(3,1))*div
!       alatinv(3,3) = (alat(1,1)*alat(2,2)-alat(1,2)*alat(2,1))*div
!
!       do iat=1,nat
!       rxyzred(1,iat)=alatinv(1,1)*rxyz(1,iat)+alatinv(1,2)*rxyz(2,iat)+alatinv(1,3)*rxyz(3,iat)
!       rxyzred(2,iat)=alatinv(2,1)*rxyz(1,iat)+alatinv(2,2)*rxyz(2,iat)+alatinv(2,3)*rxyz(3,iat)
!       rxyzred(3,iat)=alatinv(3,1)*rxyz(1,iat)+alatinv(3,2)*rxyz(2,iat)+alatinv(3,3)*rxyz(3,iat)
!       enddo
!
!     !  do j=1,3
!     !  do i=1,3
!     !  alat(i,j)=alat(i,j)
!     !  enddo
!     !  enddo
!  end subroutine
!
!
!
!
!
!
!      subroutine back2cell(nat,rxyz,alat)
!      implicit real*8 (a-h,o-z)
!      dimension rxyz(3,nat),xyzred(3,nat),alat(3,3)
!
!         call cart2frac(nat,alat,rxyz,xyzred)
!         do iat=1,nat
!         do l=1,3
!         xyzred(l,iat)=modulo(xyzred(l,iat),1.d0)
!         enddo
!         enddo
!         call frac2cart(nat,alat,xyzred,rxyz)
!
!      end subroutine
