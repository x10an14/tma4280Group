#!/bin/sh

#PBS -N openmp_1n8ppn
#PBS -A freecycle
#PBS -q optimist
#PBS -l walltime=00:08:00
#PBS -l nodes=1:ppn=8:default
#PBS -l pmem=1500MB
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine, compact"

for i in $(seq 3 28);
	do	#There will be an empty newline between each For-loop iteration (AKA "\n\n")
		echo 'k: '$i
		OMP_NUM_THREADS=8 ../../../release_openmp/ex4 $i
	done
