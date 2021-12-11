# guidedDQMC
Diffusion Quantum Monte Carlo With Gaussian Guiding Wave Functions

compile with: 
ifort mod_dftbp.f90 wrapper_energyandforces.f90 dftbp_parameters.f90 dftbp_walker.f90 dftbp_dqmc_main.f90 -mkl -qopenmp -C -g -check all -fpe0 -traceback -debug extended

