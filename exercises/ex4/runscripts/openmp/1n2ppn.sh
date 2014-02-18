#!/bin/sh

#PBS -N openmp_1n2ppn
#PBS -A freecycle
#PBS -q optimist
#PBS -l walltime=00:08:00
#PBS -l nodes=1:ppn=2:default
#PBS -l pmem=6000MB
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine, compact"

for i in $(seq 3 25);
	do	#There will be an empty newline between each For-loop iteration (AKA "\n\n")
		echo 'k: '$i
		OMP_NUM_THREADS=2 ../../../release_openmp/ex4 $i
	done
