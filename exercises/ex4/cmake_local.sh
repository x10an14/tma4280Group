#Create CMAKE sub-folders
mkdir debug release -p

#Flag to turn off OpenMP, TODO: Make bash ask if you want to use OpenMP
#  -DENABLE_OPENMP=0

#Create makefiles for each sub-folder
cd debug/
cmake .. -DCMAKE_BUILD_TYPE=Debug 

cd ../release/
cmake .. -DCMAKE_BUILD_TYPE=Release
