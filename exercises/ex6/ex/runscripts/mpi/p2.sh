#!/bin/sh

#PBS -N mpi2
#PBS -A freecycle
#PBS -q optimist
#PBS -l walltime=00:10:00
#PBS -l nodes=1:ppn=6:default
#PBS -l pmem=4096MB
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine, compact"

for i in $(seq 3 14);
	do	#There will be an empty newline between each For-loop iteration (AKA "\n\n")
		n=$(echo "2^$i" | bc)
		echo 'n: '$n
		mpirun -npernode 2 ../../../release_mpi/ex6 $n
	done
