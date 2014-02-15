#Load necessary modules (this should be done each time you log into shell).

#Create CMAKE sub-folders
mkdir debug release -p

#Delete whatever was there before
rm -rf debug/* release/*

# module load intelcomp
# module load openmpi/1.4.3-intel
# module load cmake

#Create makefiles for each sub-folder
cd debug/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug

cd ../release/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
