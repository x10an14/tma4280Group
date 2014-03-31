#!/bin/sh

#PBS -N testJob
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
echo 't: 3'
echo 'nodes: 4'
echo 'n: 4096'
OMP_NUM_THREADS=3 mpirun -npernode 4 ../../../release_omp_mpi/ex6 4096
echo
echo 't: 3'
echo 'nodes: 4'
echo 'n: 8192'
OMP_NUM_THREADS=3 mpirun -npernode 4 ../../../release_omp_mpi/ex6 8192
echo
echo 't: 3'
echo 'nodes: 4'
echo 'n: 16384'
OMP_NUM_THREADS=3 mpirun -npernode 4 ../../../release_omp_mpi/ex6 16384
