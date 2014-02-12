#Create CMAKE sub-folders
mkdir debug release -p

#Flag to turn off OpenMP, TODO: Make bash ask if you want to use OpenMP
#  -DENABLE_OPENMP=0

#Create makefiles for each sub-folder
cd debug/
CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_OPENMP=0

cd ../release/
CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_OPENMP=0
