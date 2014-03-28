#!/bin/sh

#PBS -N bigJob
#PBS -A freecycle
#PBS -q optimist
#PBS -l walltime=01:00:00
#PBS -l nodes=4:ppn=12:default
#PBS -j oe

cd ${PBS_O_WORKDIR}

module load intelcomp
module load openmpi/1.4.3-intel
KMP_AFFINITY="granularity=fine, compact"

#There will be an empty newline between each For-loop iteration (AKA "\n\n")
for (( p = 1; p < 13; p++ )); do
	if [[ "$p" == '5' ]] && [[ "$p" == '7' ]] && [[ "$p" == '9' ]] && [[ "$p" == "11" ]]; then
			continue
	fi
	echo 'p: '$p
	for (( t = 2; t*p < 49; t++ )); do
		if [[ "$t" == '5' ]] && [[ "$t" == '7' ]] && [[ "$t" == '9' ]] && [[ "$t" == "11" ]]; then
			continue
		fi
		echo 't: '$t
		if [[ "$p" -le '3' ]]; then
			x=2
		elif [[ "$p" -gt '3' ]] && [[ "$p" -lt '8' ]]; then
			#echo '		CHANGED X!!!'
			x=3
		else
			#echo '		CHANGED X!!!'
			x=4
		fi
		#echo 'x: '$x
		for (( i = x; i < 15; i++ )); do
			n=$(echo "2^$i" | bc)
			echo 'n: '$n
			OMP_NUM_THREADS=$t mpirun -npernode $p ../../../release_omp_mpi/ex6 $n
			echo
		done
	done
done
