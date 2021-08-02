module mod_dftbp

  implicit none
  !character(len=2), dimension(:), allocatable, public :: atomnames
  character(len=100), private :: path_short
  logical, private :: is_setup = .FALSE.
  logical :: restart_setup= .FALSE.
  !private
  !public :: dftbp_energy_forces, setup_dftbp

contains
  subroutine dftbp_energy_forces(nat, rxyz, fxyz, etot, alat, deralat, atomnames, thread_num)
    implicit none
    !! communicates with dftbp with files.
    !use counter
    integer :: nat
    !! Number of atoms.
    real*8, dimension(3,nat), intent(in) :: rxyz
    !! Positions of the atoms.
    real*8, dimension(3,nat) :: fxyz
    !! Forces acting on the atoms.
    real*8 :: etot
    !! Potential energy of the system.
    real*8, dimension(3,3) :: alat
    !! Lattice vectors of the periodic cell.
    real*8, dimension(3,3) :: deralat
    !! Negative derivative of the energy with respect to the lattice vectors.
    character(len=2) :: atomnames(nat)
    character(len=2), dimension(nat) :: atomtypes
    !! atom types
    integer :: ntypes
    integer :: atomnumbernames(nat)
    character(len=100):: dftb_path
    !! thread number for choosing the correct forlder in dftp+
    integer :: thread_num
    character(len=10) :: thread_str

    write(thread_str, '(I3.3)') thread_num
    dftb_path = trim(adjustl(path_short))//trim(adjustl(thread_str))//"/"

    if ( .not. is_setup ) then
      stop "dftbp is not set up. call setup_dftb."
    end if
    call get_atom_types(atomnames, nat, ntypes, atomtypes)
    call atomnames2number(atomnames, nat, ntypes, atomtypes, atomnumbernames)
    !call system("rm detailed.out")
    call create_input_files(nat, alat, atomnumbernames, atomtypes, ntypes, rxyz, dftb_path)
    call run_dftb(dftb_path)
    call get_results(nat, etot, fxyz, deralat, dftb_path)
    deralat = -deralat
  end subroutine dftbp_energy_forces


  subroutine setup_dftbp(path,thread_num)
    implicit none
    character(len=100) :: command
    integer :: i, j
    character(len=100):: dftb_path, path
    !! thread number for choosing the correct forlder in dftp+
    integer :: thread_num
    character(len=10) :: thread_str

    path_short = path
    write(thread_str, '(I3.3)') thread_num
    dftb_path = trim(adjustl(path_short))//trim(adjustl(thread_str))//"/"

    is_setup = .TRUE.
    call cleanup_dftbp(dftb_path)

    command = "mkdir -p "//trim(adjustl(dftb_path))
    call system(trim(command))

    command = "cp dftb_in.hsd "//trim(adjustl(dftb_path))
    call system(trim(command))

    command = "cp -a rundftbp.sh "//trim(adjustl(dftb_path))
    call system(trim(command))

    command = "cp input_kpoints.gen "//trim(adjustl(dftb_path))
    call system(trim(command))

    command = "cp -r 3ob-3-1 "//trim(adjustl(dftb_path))
    call system(trim(command))
    !command = "cp -r pbc-0-3 "//trim(adjustl(dftb_path))
    !call system(trim(command))
  end subroutine setup_dftbp

  subroutine cleanup_dftbp(dftb_path)
    implicit none
    character(len=250) :: command
    character(len=100):: dftb_path
    !print*, "deleting temporary dftbp files."
    !call system(trim(dftb%dftb_path)//"rm input_restart.gen")
    !call system(trim(dftb%dftb_path)//"rm dftbp.out")
    !call system(trim(dftb%dftb_path)//"rm charges.bin")
    !call system(trim(dftb%dftb_path)//"rm dftb_pin.hsd")
    !print*, "rm -r "//dftb%dftb_path
    command = "rm -rf "//trim(adjustl(dftb_path))//" *"
    !call system(command)
  end subroutine cleanup_dftbp

  subroutine atomnames2number(atomnames_in, nat, ntypes, atomtypes, atomnumbernames)
    implicit none
    integer :: nat
    !! Number of atoms
    character(len=2), dimension(nat) :: atomnames_in
    !! atomnames
    character(len=2), dimension(nat) :: atomtypes
    !! atom types
    integer :: ntypes
    integer :: atomnumbernames(nat)
    integer :: i, j

    do i = 1, nat, 1
      do j = 1, ntypes, 1
        if ( atomnames_in(i) == atomtypes(j) ) then
          atomnumbernames(i) = j
        end if
      end do
    end do
  end subroutine atomnames2number

  subroutine run_dftb(dftb_path)
      implicit none
      character(len=100):: dftb_path
      character(len=100) :: runfile
      integer :: i, j
      runfile = "cd "//trim(adjustl(dftb_path))//" && ./rundftbp.sh"
      call system(runfile)
      !do j = 1, 30, 1
      !  call execute_command_line(trim(adjustl(runfile)), cmdstat=i)
      !  if (i==0) exit
      !    print*, runfile, "failed to execute, error: ", i, "j= ",j
      !end do
      !call execute_command_line(trim(adjustl(runfile)), cmdstat=i)
      !if ( i /= 0 ) then
      !  print*, runfile, "failed to execute, error: ", i
      !  call execute_command_line(trim(adjustl(runfile)), cmdstat=i)
      !end if
  end subroutine

  subroutine create_input_files(nat, alat, atomnumbernames, atomtypes, ntypes, rxyz, dftb_path)
    implicit none
    character(len=100):: dftb_path
    integer :: nat
    !! Number of atoms
    real*8, intent(in), dimension(3,3) :: alat
    !! Lattice vectors
    integer, dimension(nat) :: atomnumbernames
    !! atomnames in number
    character(len=2), dimension(nat) :: atomtypes
    !! atom types
    integer :: ntypes
    !! number of different atom types
    real*8, dimension(3,nat), intent(in) :: rxyz
    character(len=250) :: inputgeo

    integer :: u, ityp, iat

    inputgeo=trim(adjustl(dftb_path))//"input_geometry.gen"
    open(newunit=u,file=inputgeo)
    write(u,'(i5.5,a)') nat, " C"
    !write(u,'(i5.5,a)') nat, " S"
    write(u,*) (trim(adjustl(atomtypes(ityp)))//" ", ityp=1,ntypes)
    do iat = 1, nat
        write(u,'(i5,1x,i5,3(1x,es25.15))') iat, atomnumbernames(iat), rxyz(:,iat)
        !write(u,*) iat, atomnames(iat), rxyz(:,iat)
    end do

    !write(u,'(3(1x,es25.15))') 0.d0,0.d0,0.d0
    !write(u,'(3(1x,es25.15))') alat(:, 1)
    !write(u,'(3(1x,es25.15))') alat(:, 2)
    !write(u,'(3(1x,es25.15))') alat(:, 3)
    close(u)

    open(newunit=u,file=trim(dftb_path)//"input_restart.gen")

    !if(restart_setup==.FALSE.) then
    !  restart_setup=.TRUE.
    !  write(u,'(a)')'No'
    !else
    !  write(u,'(a)')'Yes'
    !end if
    write(u,'(a)')'No'
    !write(u,'(a)')'Yes'
    close(u)

  end subroutine create_input_files

  subroutine get_results(nat, epot, fxyz, deralat, dftb_path)
    implicit none
    character(len=100):: dftb_path
    !parameters
    integer, intent(in) :: nat
    real*8, intent(out) :: fxyz(3,nat)
    real*8, intent(out) :: epot
    !internal
    integer :: u,n,k,m,l,iat,scratch
    real*8 :: str_matrix(3,3)
    real*8 :: deralat(3,3)
    character(len=250) :: all_line
    integer :: ios

    open(newunit=u, file=trim(dftb_path)//"detailed.out", iostat=ios, status="old", action="read")
    if ( ios /= 0 ) stop "Error opening file detailed.out. check dftbp.out for information"
    do
        read(u,'(a250)',end=99)all_line
        n = len_trim(all_line)
        k = index(all_line(1:n),"Total Mermin free energy")
        if(k.ne.0) then
            m = len_trim(all_line)
            l = scan(all_line(1:m),":",.true.)
            read(all_line(l+1:m),*) epot
            cycle
        endif
        k = index(all_line(1:n),"Total Forces")
        if(k.ne.0) then
          do iat=1,nat
              read(u,*) scratch, fxyz(:,iat)
          enddo
          cycle
        endif
        k = index(all_line(1:n),"Total lattice derivs")
        if(k.ne.0) then
            do iat=1,3
              read(u,*) deralat(:,iat)
            enddo
            cycle
        endif
        k = index(all_line(1:n),"Total stress tensor")
        if(k.ne.0) then
            do iat=1,3
              read(u,*) str_matrix(:,iat)
            enddo
            cycle
        endif
    enddo

    99 continue
    close(u)

  end subroutine get_results

  subroutine get_atom_types(atomnames_in, nat, ntypes, atomtypes)
    implicit none
    integer :: nat
    !! Number of atoms
    character(len=2), dimension(nat) :: atomnames_in
    !! atomnames
    character(len=2), dimension(nat) :: atomtypes
    !! atom types
    integer :: ntypes
    integer i,j

    ntypes = 1
    atomtypes(1) = atomnames_in(1)
    do i = 2, nat, 1
      do j = 1, ntypes, 1
        if ( atomnames_in(i) == atomtypes(j) ) then
          goto 123
        end if
      end do
      ntypes = ntypes + 1
      atomtypes(ntypes) = atomnames_in(i)
      123 continue
    end do
  end subroutine get_atom_types

end module mod_dftbp
