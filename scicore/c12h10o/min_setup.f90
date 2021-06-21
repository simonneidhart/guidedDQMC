

 program main
  use mod_dftbp
   implicit real*8 (a-h,o-z)
   real*16 e_path
   parameter(nat=180)  ! number of atoms
   character(len=30) :: units, names, path
   !character(len=2) :: symb !exchanged for atomnames
   dimension rxyz0(3,nat),rxyz(3,nat),alat(3,3),fxyz(3,nat),dd(3,nat,2),displ(3,nat),deralat(3,3),xyzred(3,nat)
   dimension dda(3,3,2),displa(3,3)
   dimension alat0(3,3)
   character(len=20) :: filetype
   character(len=2) :: atomnames(nat)
   character(len=100) :: k_point_cp_string,filename,file_result, file_w_path
   real (8) :: Bohr_Ang
   !!Integers für system clock timing
   integer*8 :: c,cr,cm,TimeStart,TimeEnd

  !! position step size
  real*8 :: alpha_pos = 0.5d0
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


  Bohr_Ang = 0.529d0
  pi=4.d0*atan(1.d0)


   !call setup_dftbp("/dev/shm/dftb_temp_folder_c8_new_1/")

   filetype = 'ascii'  ! periodic boundary conditions
   filename='SiC_20to180_mp_displ_03.ascii' !SiC_32_mp_displ_01 !SiC_08_mp_displ_02
   path='rnd_displaced/'

   write(file_w_path,'(a14,a)') path,trim(filename)


   if (trim(filetype).eq."ascii") then !SiC_08_mp_displ_01
        open(10, file=trim(file_w_path), status="old") !posinp.ascii !SiC_8_exp_72  !ALATcircle_055 !SiC_not_gm !SiC_20 !0624 !SiC_20_mp !SiC_20_mp_changed !SiC_8_mp_changed
        !If init gets used change dftb_in.hsd file to C-O-H + H bridges!! (See dftb_in_mo.hsd file)

         read(10,*) natp, names, units

          if (trim(units).ne.'atomic' .and. trim(units).ne.'bohr')  then
            units='ang'
            write(*,*) "Input in Angstroem"
          end if

         if (natp.ne.nat) write(*,*)"Number of atoms does not equal the defined parameter!!!!"
         read(10,*) dxx,dyx,dyy
         read(10,*) dzx,dzy,dzz

         alat0(1,1)=dxx
         alat0(2,1)=0.d0
         alat0(3,1)=0.d0

         alat0(1,2)=dyx
         alat0(2,2)=dyy
         alat0(3,2)=0.d0

         alat0(1,3)=dzx
         alat0(2,3)=dzy
         alat0(3,3)=dzz

         do i = 1, nat
            read(10, *) ( rxyz0(j, i), j = 1, 3) ,atomnames(i)
         end do
  end if  !  ascii

  write(*,*) "Reading input file completed"

  !Transofrm angstroem input to atomic units
  if (units .eq. 'ang') then
    alat0=alat0/Bohr_Ang
    rxyz0=rxyz0/Bohr_Ang
  end if


        alat=alat0
        rxyz=rxyz0
        alpha_pos = 1.d0
        alpha_lat = 2.d0


        !call setup_dftbp("/dev/shm/dftb_temp_folder_r17/")
        call setup_dftbp("dftb_temp_folder/")

        !!energyandforces is wrapper for dftbp
        call energyandforces(nat, rxyz, fxyz, etot, alat, deralat, atomnames)
        write(*,*) "Input etot",etot

        






end program main
