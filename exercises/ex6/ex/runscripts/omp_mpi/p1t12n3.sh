#!/bin/sh

#PBS -N combined_p1t12n1
#PBS -A freecycle
#PBS -q optimist
#PBS -l walltime=00:10:00
#PBS -l nodes=3:ppn=12:default
###PBS -l pmem=750MB
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine, compact"

for i in $(seq 3 14);
	do	#There will be an empty newline between each For-loop iteration (AKA "\n\n")
		n=$(echo "2^$i" | bc)
		echo 'n: '$n
		OMP_NUM_THREADS=12 mpirun -npernode 1 ../../../release_omp_mpi/ex6 $n
	done
