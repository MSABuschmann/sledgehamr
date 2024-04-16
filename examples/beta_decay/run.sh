#!/bin/bash
#SBATCH --constraint=gpu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-task=1
#SBATCH --qos=debug
#SBATCH --time=00:01:00
#SBATCH -A m3166_g

cd $SLURM_SUBMIT_DIR

export SLURM_CPU_BIND="cores"
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=16
export HDF5_DISABLE_VERSION_CHECK=2
export HDF5_USE_FILE_LOCKING='FALSE'

/usr/bin/time --verbose srun main3d.gnu.MPROF.MPI.OMP.CUDA.ex inputs
