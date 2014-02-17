#Create CMAKE sub-folder
mkdir ../release -p

#Delete whatever was there before
rm -rf ../release/*

cd ../release/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release  -DENABLE_OPENMP=0

make clean && make
chmod +x ex4

cd ../runscripts/mpi/
chmod +x 1n1ppn.sh
chmod +x 1n2ppn.sh
chmod +x 1n4ppn.sh
chmod +x 1n8ppn.sh
chmod +x 2n1ppn.sh
chmod +x 2n2ppn.sh
chmod +x 2n4ppn.sh
chmod +x 2n8ppn.sh

cd kongull_outpt
qsub ../1n1ppn.sh
qsub ../1n2ppn.sh
qsub ../1n4ppn.sh
qsub ../1n8ppn.sh
qsub ../2n1ppn.sh
qsub ../2n2ppn.sh
qsub ../2n4ppn.sh
qusb ../2n8ppn.sh
