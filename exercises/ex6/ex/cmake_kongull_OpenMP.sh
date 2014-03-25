#Load necessary modules (this should be done each time you log into shell).

#Create CMAKE sub-folders
mkdir debug_omp_mpi release_omp_mpi -p

#Delete whatever was there before
rm -rf debug_omp_mpi/* release_omp_mpi/*

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake

#Create makefiles for each sub-folder
cd debug_omp_mpi/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug
make clean && make && chmod +x ex6

cd ../release_omp_mpi/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
make clean && make && chmod +x ex6
