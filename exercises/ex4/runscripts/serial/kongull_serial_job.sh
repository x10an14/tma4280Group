#!/bin/sh

#PBS -N serial_1n1ppn
#PBS -A freecycle
#PBS -q optimist
#PBS -l walltime=00:08:00
#PBS -l nodes=1:ppn=1:default
#PBS -l pmem=12000MB
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine, compact"

for j in $(seq 1 2);
do
	for i in $(seq 3 28);
		do	#There will be an empty newline between each For-loop iteration (AKA "\n\n")
			echo 'k: '$i
			../../../release_serial/ex4 $i
		done
done
