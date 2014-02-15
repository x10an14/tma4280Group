#!/bin/sh

#PBS -N ex4_jobscript
#PBS -A freecycle
#PBS -l walltime=00:08:00
#PBS -l nodes=2:ppn=4:default
#PBS -q optimist
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine, compact"

for i in $(seq 4 14);
	do	#There will be an empty newline between each For-loop iteration (AKA "\n\n")
		echo 'k: '$i
		mpirun -npernode 4 ../release/ex4 $i
	done