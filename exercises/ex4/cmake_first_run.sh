#Create CMAKE sub-folders
mkdir debug release -p

#Create makefiles for each sub-folder
cd debug/
cmake .. -DCMAKE_BUILD_TYPE=Debug

cd ../release/
cmake .. -DCMAKE_BUILD_TYPE=Release
