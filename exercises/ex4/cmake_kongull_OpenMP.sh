#Load necessary modules (this should be done each time you log into shell).
module load intelcomp/13.0.1
module load openmpi/1.4.3-intel
module load cmake

#Create CMAKE sub-folders
mkdir debug release -p

#Create makefiles for each sub-folder
cd debug/
CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug

cd ../release/
CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
