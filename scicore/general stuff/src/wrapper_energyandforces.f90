

!! 2 wrappers for the force calculating subroutine.

! This variant of the subroutine calls the Bazant force field
!!  subroutine energyandforces(nat, rxyz, fxyz, etot, alat, deralat, atomnames)
!!  implicit none
!!  !! Number of atoms.
!!  integer :: iat,nat
!!  !! Positions of the atoms.
!!  real*8, dimension(3,nat), intent(in) :: rxyz
!!  !! Forces acting on the atoms.
!!  real*8, dimension(3,nat) :: fxyz
!!  !! Potential energy of the system.
!!  real*8 :: etot
!!  !! Lattice vectors of the periodic cell.
!!  real*8, dimension(3,3) :: alat
!!  !! Negative derivative of the energy with respect to the lattice vectors.
!!  real*8, dimension(3,3) :: deralat
!!  !! array containg the chemical symbols of each atom
!!  character(len=2) :: atomnames(nat)
!!
!!  !! private variable
!!
!!  !! make sure that that the positions are not changed by the force calculating program.
!!  call  energyandforces_bazant(nat,alat,rxyz,etot,fxyz,deralat)
!!  do iat=1,nat
!!  atomnames(iat)="Si"
!!  enddo
!!end subroutine energyandforces


!This variant will call DFTB+
subroutine energyandforces(nat, rxyz, fxyz, etot, alat, deralat, atomnames)
!! input rxyz and alat must be in atomic units (bohr)
  use mod_dftbp
  implicit none
  !! Number of atoms.
  integer :: iat,nat,i,j
  !! Positions of the atoms.
  real*8, dimension(3,nat), intent(in) :: rxyz
  !! Forces acting on the atoms.
  real*8, dimension(3,nat) :: fxyz
  !! Potential energy of the system.
  real*8 :: etot
  !! Lattice vectors of the periodic cell.
  real*8, dimension(3,3) :: alat
  !! Negative derivative of the energy with respect to the lattice vectors.
  real*8, dimension(3,3) :: deralat
  !! array containg the chemical symbols of each atom
  character(len=2) :: atomnames(nat)

  !! private variable
  real*8, dimension(3,nat) :: rxyz_ang
  real*8, dimension(3,3) :: alat_ang

  !! make sure that that the positions are not changed by the force calculating program.
  !! Strange enough the inputs to DFTB+ should be in angstroem (even though all
  !! output variables are in atomic units
  real*8 :: Bohr_Ang=0.529177d0  ! Conversion factor Bohr to angstroem


  do iat=1,nat
  rxyz_ang(1,iat)=rxyz(1,iat)*Bohr_Ang
  rxyz_ang(2,iat)=rxyz(2,iat)*Bohr_Ang
  rxyz_ang(3,iat)=rxyz(3,iat)*Bohr_Ang
  enddo
  do j=1,3
  do i=1,3
  alat_ang(i,j)=alat(i,j)*Bohr_Ang
  enddo
  enddo


  ! call the desired potential potential here.
  call dftbp_energy_forces(nat, rxyz_ang, fxyz, etot, alat_ang, deralat, atomnames)



end subroutine energyandforces
