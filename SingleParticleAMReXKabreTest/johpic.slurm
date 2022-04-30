#!/bin/sh
#SBATCH --job-name=JohPIC
#SBATCH --partition=nu-all
#SBATCH --ntasks=4
#SBATCH --nodes=4
#--nodelist=nu-3c.cnca,nu-3d.cnca
#SBATCH --exclude=nu-0d.cnca
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=03:00:00

export OMP_NUM_THREADS=64 
export GOMP_CPU_AFFINITY=0-61
mpirun -np 4 ./main inputs



#####################################################################

/home/jovillalobos/home/amrex/installdir