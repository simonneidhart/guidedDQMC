#!/bin/sh
#SBATCH --job-name="dftbp_guided_hessian"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --qos=kernph-1day --partition=psi
#SBATCH --mem-per-cpu=2500M


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "SLURM_NTASKS: $SLURM_NTASKS"
echo "SLURM_NTASKS_PER_NODE: $SLURM_NTASKS_PER_NODE"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"


module purge
ml intel

 srun a.out
