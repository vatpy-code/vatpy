#!/bin/sh
#SBATCH --job-name relax
#SBATCH --error err.relax
#SBATCH --output out.relax

#SBATCH --ntasks 32
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 7000
#SBATCH --qos serial
#SBATCH --time 0-05:00:00

#SBATCH -p infiniband

#SBATCH --exclude=jst[327-331]

ml purge
. $HOME/Petersson/modules/ml_load_scitas.sh

export UCX_TLS=rc,dc,ud,sysv,posix,self

SIMDIR=$(pwd)

if [[ ! -d $SIMDIR/OUTPUT ]]; then
	mkdir OUTPUT
fi

srun ./Arepo param.txt > Arepo.out
