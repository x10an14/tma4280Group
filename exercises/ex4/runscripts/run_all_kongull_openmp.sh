cd /home/$USER/tma4280Group/exercises/ex4/runscripts/
#Create CMAKE sub-folder
mkdir r../elease_openmp -p

#Delete whatever was there before
rm -rf ../release_openmp/*

cd ../release_openmp/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
make clean && make
chmod +x ex4

cd ../runscripts/openmp/
chmod +x 1n1ppn.sh
chmod +x 1n2ppn.sh
chmod +x 1n4ppn.sh
chmod +x 1n8ppn.sh

mkdir 1n1ppn -p && cd 1n1ppn
qsub ../1n1ppn.sh

mkdir ../1n2ppn -p && cd ../1n2ppn
qsub ../1n2ppn.sh

mkdir ../1n4ppn -p && cd ../1n4ppn
qsub ../1n4ppn.sh

mkdir ../1n8ppn -p && cd ../1n8ppn
qsub ../1n8ppn.sh
