#Create CMAKE sub-folders
mkdir debug release -p

#Create makefiles for each sub-folder
cd debug/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Debug

cd ../release/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
