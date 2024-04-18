#!/bin/bash
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --qos=debug
#SBATCH --time=00:02:00
cd $SLURM_SUBMIT_DIR

export SLURM_CPU_BIND="cores"
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=16

srun main3d.gnu.x86-milan.MPROF.MPI.OMP.ex inputs
