cd /home/$USER/tma4280Group/exercises/ex4
#Create CMAKE sub-folder
mkdir release -p

#Delete whatever was there before
rm -rf release/*

cd release/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
make clean && make

chmod +x ex4

cd ../runscripts/combined/
mkdir 1n1ppn 1n2ppn 1n4ppn 1n8ppn -p

cd 1n1ppn/
chmod +x ../n1ppn.sh
qsub ../n1ppn.sh

cd ../1n2ppn/
chmod +x ../n2ppn.sh
qsub ../n2ppn.sh

cd ../1n4ppn/
chmod +x ../n4ppn.sh
qsub ../n4ppn.sh

cd ../1n8ppn/
chmod +x ../n8ppn.sh
qsub ../n8ppn.sh