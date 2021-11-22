program dmc
  use mod_dftbp
  use parameters
  implicit real(8) (a-h,o-z)
  !dftbp parameters:------------------------------------------------------------
  real*16 e_path
  character(len=20) :: filetype, extra
  character(len=100) :: k_point_cp_string,filename,file_result, file_w_path
  character(len=100):: dftb_path, line2, line
  !!Integers for system clock timing
  integer*8 :: c,cr,cm,TimeStart,TimeEnd
  !! position step size
  real*8 :: alpha_pos = 1.d0
  !! lattice step size
  real*8 :: alpha_lat = 2.d0
  !! minimal step size
  real*8 :: alpha0 = 1.d-3
  !! Set true when forces should be preconditioned
  logical :: for_pre = .FALSE.!.FALSE. TRUE
  !! force tolerance
  real*8 :: fnrmtol = 1.d-5
  !! lattice derivative tolerance
  real*8 :: lat_tol = 1.d-5
  !! set true if debug information should be written.
  logical :: debug = .TRUE.
  !! set true if every configuration from each iteration should be written to disk
  logical :: write_configurations = .false.
  !! max number of iterations
  integer :: nit = 3000
  !iterations optimizer
  integer :: it_opt
  !simulation variables:--------------------------------------------------------
  integer :: energy_count = 0, a !count the number of energy calculations
  real(8), allocatable, dimension(:,:) :: fxyz, rxyz
  real(8), allocatable, dimension(:) :: energies
  character(len=2) :: atomname

  !Simulation parameters:-------------------------------------------------------

  call read_hessian(e0_theory, .FALSE.)

  !initialzations and memory allocations----------------------------------------
  dftb_path = "dftb_temp_folder"
  call setup_dftbp(dftb_path, 0)

  allocate(fxyz(3,nat))
  allocate(rxyz(3,nat))
  allocate(energies(100))

  call energyandforces(nat, rxyz0, fxyz, ener, alat, deralat, atomnames, 0)
  print*, ener
  rxyz = rxyz0

  !gradient descent to find other minima
  !rot = 21.d0
  !call rotate_posinp(rot, 1.d0, 0.d0)

  !open(unit=9,file='rotate_posinp.xyz',status='unknown')
  !read(9,*) line2
  !read(9,*) line2
  !do iat=61,nat
  !   read(9,'(a100)')line
  !   read(line,*,iostat=ierror) atomname,rxyz(1,iat),rxyz(2,iat),rxyz(3,iat) ,extra
  !   if (ierror .ne. 0) then
  !      read(line,*,iostat=ierror) atomname,rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
  !      extra ='  '
  !   end if
  !   rxyz(:,iat)=rxyz(:,iat)/Bohr_Ang
  !end do

  call energyandforces(nat, rxyz, fxyz, ener, alat, deralat, atomnames, 0)
  print*, ener

  call gradient_descent(rxyz)


  !rotate or translate along different directions
  !do i=1,size(energies)
    !trans = i*0.01d0 - 0.01d0
    !rot = i*1.d0
    !call rotate_posinp(rot, 1.d0, 0.d0)
    !rxyz = rxyz0

    !open(unit=9,file='rotate_posinp.xyz',status='unknown')
    !read(9,*) line2
    !read(9,*) line2
    !do iat=61,nat
    !   read(9,'(a100)')line
    !   read(line,*,iostat=ierror) atomname,rxyz(1,iat),rxyz(2,iat),rxyz(3,iat) ,extra
    !   if (ierror .ne. 0) then
    !      read(line,*,iostat=ierror) atomname,rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
    !      extra ='  '
    !   end if
    !   rxyz(:,iat)=rxyz(:,iat)/Bohr_Ang
    !end do
    !rxyz(1,62) = i*0.01d0 - 0.01d0
    !print*, rxyz(1,62)
    !call energyandforces(nat, rxyz, fxyz, ener, alat, deralat, atomnames, 0)


    !energies(i) = ener
  !end do

  !calculate PES in hessian coordinate system along the vibrational modes
  !do k=1,dim
  !  mu_w = mu
  !  ener = 0
  !  mu_w(k) = mu_w(k)-50*0.1
  !  do j=1,100
  !    r_temp = matmul(ev,mu_w)
  !    rxyz = reshape(r_temp,shape(rxyz)) + rxyz0
  !    call energyandforces(nat, rxyz, fxyz, ener, alat, deralat, atomnames, 0)
  !    mu_w(k) = mu_w(k) + 0.1
  !    energies(j,k) = ener
  !  end do
  !  print*, "complete dimension ",k
  !end do


  !output final positions of the walkers

  !rxyz = rxyz*Bohr_Ang
  open(unit=13, file="pos_gradient.xyz", iostat=ios, action="write")
  if ( ios /= 0 ) stop "Error opening file unit 13"
  do i=1,nat
    write(unit=13, fmt=*) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
  end do
  close(unit=13, iostat=ios, status="keep")
  if ( ios /= 0 ) stop "Error closing file unit 13"


  !open(unit=13, file="ener", iostat=ios, action="write")
  !if ( ios /= 0 ) stop "Error opening file unit 13"
  !do i=1,size(energies)
  !  write(unit=13, fmt=*) energies(i)
  !end do
  !close(unit=13, iostat=ios, status="keep")
  !if ( ios /= 0 ) stop "Error closing file unit 13"

  print*, "all done"

end program dmc

subroutine gradient_descent(rxyz)
  use parameters
  use mod_dftbp
  implicit none
  real(8) :: etot,alpha=1d-3,norm=1,vectornorm=0,vectornorm2=0,angle,angle1
  real(8) ,  dimension(3,nat) :: rxyz
  real(8) ,  dimension(3,nat) :: fxyz, prev_fxyz
  integer :: i, j, ios, iter=0, k

  do i=1,nat
    write(unit=*, fmt=*) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
  end do

  call energyandforces(nat, rxyz, fxyz, etot, alat, deralat, atomnames, 0)
  write(unit=*, fmt=*) "The energy in the given state is: ",etot

  do while ( norm>1d-6 )
    iter = iter+1 !count iterations
    prev_fxyz = fxyz
    rxyz = rxyz+  fxyz*alpha !updating the positions
    call energyandforces(nat, rxyz, fxyz, etot, alat, deralat, atomnames, 0) !calculates negative gradient and total energy

    !calculating the gradient norm and the angle between the vectors
    angle = 0
    norm = 0
    do i=1,nat
      vectornorm = 0
      vectornorm2 = 0
      do j=1,3
        vectornorm = vectornorm+fxyz(j,i)**2
        vectornorm2 = vectornorm2+prev_fxyz(j,i)**2
      end do
      norm = norm + sqrt(vectornorm)/nat !average over all the gradient norms
      !angle between the current and previous gradient vector of the atom i
      angle=angle+acos(dot_product(prev_fxyz(1:3,i),fxyz(1:3,i))/(sqrt(vectornorm)*sqrt(vectornorm2)))
    end do
    if (angle/nat > 1.0472) then  !average of all the angles, gradient feedback
      alpha = alpha/2
    else
      alpha = alpha*1.05
    endif
  end do

  write(unit=*, fmt=*) "The energy after the minimization is: ", etot
  write(unit=*, fmt='(I5,A)') iter, " iterations were needed to obtain a precision of 1d-6"
  write(unit=*, fmt=*)

end subroutine gradient_descent


subroutine rotate_posinp(rot, scale, trans)

    implicit none
    integer, parameter :: natx=2000
    character(len=5) :: atomname(natx)
    character(len=12) :: units
    character(len=20) :: extra(natx)
    character(len=120) :: line,line2
    real(kind=8), dimension(3,natx) :: pos
    real(kind=8), dimension(3) :: pos_s
    real(kind=8) :: scale,t1,t2,t3,phi_1,phi_2,phi_3,trans, rot
    real(kind=8), parameter :: pi1=3.141592654d0
    integer :: nat,iat,ierror

    !write(*,*) 'reading atomic positions from file posinp'
    open(unit=9,file='posinp.xyz',status='unknown')
    read(9,*) nat,units
    !write(*,*) nat,units
    if (nat.gt.natx) stop 'increase natx'
    read(9,*) line2
    !write(*,*) line2

    ! put center of mass at origin
    pos_s(1)=0.d0
    pos_s(2)=0.d0
    pos_s(3)=0.d0
    do iat=1,nat
       !write(*,*)  '      ---------------------------------- '
       read(9,'(a100)')line
       !write(*,*) iat, line

       read(line,*,iostat=ierror) atomname(iat),pos(1,iat),pos(2,iat),pos(3,iat) ,extra(iat)
       if (ierror .ne. 0) then
          read(line,*,iostat=ierror) atomname(iat),pos(1,iat),pos(2,iat),pos(3,iat)
          extra(iat)='  '
       end if
       !write(*,*) atomname(iat),pos(1,iat),pos(2,iat),pos(3,iat)

       pos_s(1)=pos_s(1)+pos(1,iat)
       pos_s(2)=pos_s(2)+pos(2,iat)
       pos_s(3)=pos_s(3)+pos(3,iat)
    end do
    close(unit=9)
    pos_s(1)=pos_s(1)/nat
    pos_s(2)=pos_s(2)/nat
    pos_s(3)=pos_s(3)/nat
    do iat=1,nat
       pos(1,iat)=pos(1,iat)-pos_s(1)
       pos(2,iat)=pos(2,iat)-pos_s(2)
       pos(3,iat)=pos(3,iat)-pos_s(3)
    end do

    !write(*,*) 'rotations in degrees (0<= Phi <=360):'
    !write(*,*)
    !write(*,*) 'around z axis / in xy-plane:'
    !read(*,*) phi_1
    !phi_1=2.d0*pi1*rot/360.d0
    phi_1 = 0.d0
    !write(*,*) 'around y axis / in xz-plane:'
    phi_2=2.d0*pi1*rot/360.d0
    !phi_2 = 0.d0
    !write(*,*) 'around x axis / in yz-plane:'
    !read(*,*) phi_3
    phi_3 = 0.d0
    !phi_3=2.d0*pi1*phi_3/360.d0

    do iat=1,nat
       t1=cos(phi_1)*pos(1,iat)+sin(phi_1)*pos(2,iat)
       t2=-sin(phi_1)*pos(1,iat)+cos(phi_1)*pos(2,iat)
       pos(1,iat)=t1
       pos(2,iat)=t2

       t1=cos(phi_2)*pos(1,iat)+sin(phi_2)*pos(3,iat)
       t3=-sin(phi_2)*pos(1,iat)+cos(phi_2)*pos(3,iat)
       pos(1,iat)=t1
       pos(3,iat)=t3

       t2=cos(phi_3)*pos(2,iat)+sin(phi_3)*pos(3,iat)
       t3=-sin(phi_3)*pos(2,iat)+cos(phi_3)*pos(3,iat)
       pos(2,iat)=t2
       pos(3,iat)=t3
    end do

    ! shift back center of mass to original position
    do iat=1,nat
       pos(1,iat)=pos(1,iat)+pos_s(1)
       pos(2,iat)=pos(2,iat)+pos_s(2)
       pos(3,iat)=pos(3,iat)+pos_s(3)
    end do

    !write(*,*) 'scaling factor=?'
    !read(*,*) scale
    !write(*,*) 'translation'

    !write(*,*) 'writing atomic positions to file rot_posinp'
    open(unit=9,file='rotate_posinp.xyz',status='unknown')
    write(9,*) nat, units
    write(9,'(a100)') line2
    do iat=1,nat
       write(9,'(a5,3x,3(1x,e17.10),4x,a)') atomname(iat),   &
pos(1,iat)*scale+trans,pos(2,iat)*scale,pos(3,iat)*scale,extra(iat)
    end do
    close(unit=9)

end subroutine rotate_posinp
