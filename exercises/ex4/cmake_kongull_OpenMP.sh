#Create CMAKE sub-folders
mkdir debug release -p

#Create makefiles for each sub-folder
cd debug/
CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug

cd ../release/
CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
