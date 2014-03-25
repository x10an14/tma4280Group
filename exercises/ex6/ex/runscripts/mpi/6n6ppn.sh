#!/bin/sh

#PBS -N combined_2n8ppn
#PBS -A freecycle
#PBS -q optimist
#PBS -l walltime=00:08:00
#PBS -l nodes=6:ppn=12:default
###PBS -l pmem=750MB
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine, compact"

for i in $(seq 2 14 2);
	do	#There will be an empty newline between each For-loop iteration (AKA "\n\n")
		n=$(echo "2^$i" | bc)
		echo 'k: '$i
		echo 'n: '$n
		mpirun -npernode 6 ../../../release/ex6 $n
	done