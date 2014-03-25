#Create CMAKE sub-folders
mkdir debug_mpi release_mpi -p

#Delete whatever was there before
rm -rf debug_mpi/* release_mpi/*

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake

#Flag to turn off OpenMP, TODO: Make bash ask if you want to use OpenMP
#  -DENABLE_OPENMP=0

#Create makefiles for each sub-folder
cd debug_mpi/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_OPENMP=0
make clean && make && chmod +x ex6

cd ../release_mpi/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_OPENMP=0
make clean && make && chmod +x ex6
