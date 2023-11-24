#!/bin/bash
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --qos=debug
#SBATCH --time=00:15:00
cd $SLURM_SUBMIT_DIR

export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=16
export HDF5_DISABLE_VERSION_CHECK=2
export HDF5_USE_FILE_LOCKING='FALSE'

/usr/bin/time --verbose srun --cpu_bind=cores main3d.gnu.x86-milan.MPROF.MPI.OMP.ex inputs
