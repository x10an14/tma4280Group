#Create CMAKE sub-folders
mkdir debug release -p

#Delete whatever was there before
rm -rf debug/* release/*

module load intelcomp
module load openmpi/1.4.3-intel
module load cmake

#Flag to turn off OpenMP, TODO: Make bash ask if you want to use OpenMP
#  -DENABLE_OPENMP=0

#Create makefiles for each sub-folder
cd debug/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_OPENMP=0

cd ../release/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_OPENMP=0
