#!/bin/sh
#SBATCH --job-name relax
#SBATCH --error err.relax
#SBATCH --output out.relax

#SBATCH --ntasks 32
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 7000
#SBATCH --qos serial
#SBATCH --time 0-01:00:00


module purge

ml load gcc/11.3.0
ml load openmpi/4.1.3
ml load gsl/2.7.1
ml load gmp/6.2.1
ml load fftw/3.3.10
ml load hdf5/1.12.2-mpi

SIMDIR=$(pwd)

if [[ ! -d $SIMDIR/OUTPUT ]]; then
	mkdir OUTPUT
fi

srun ./Arepo param.txt > Arepo.out
