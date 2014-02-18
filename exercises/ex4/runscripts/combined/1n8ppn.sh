#!/bin/sh

#PBS -N combined_1n8ppn
#PBS -A freecycle
#PBS -q optimist
#PBS -l walltime=00:08:00
#PBS -l nodes=1:ppn=8:default
###PBS -l pmem=1500MB
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine, compact"

for i in $(seq 8 28);
	do	#There will be an empty newline between each For-loop iteration (AKA "\n\n")
		echo 'k: '$i
		mpirun -npernode 8 ../../../release/ex4 $i
	done
