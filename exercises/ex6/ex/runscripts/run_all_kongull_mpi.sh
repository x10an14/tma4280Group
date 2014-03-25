cd /home/$USER/tma4280Group/exercises/ex6/runscripts/
#Create CMAKE sub-folder
mkdir ../release_mpi -p

#Delete whatever was there before
rm -rf ../release_mpi/*

cd ../release_mpi/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release  -DENABLE_OPENMP=0
make clean && make
chmod +x ex6

cd ../runscripts/mpi/
chmod +x 1n4ppn.sh
chmod +x 1n6ppn.sh
chmod +x 1n8ppn.sh
chmod +x 1n12ppn.sh
chmod +x 3n12ppn.sh
chmod +x 4n9ppn.sh
chmod +x 6n6ppn.sh
chmod +x 9n4ppn.sh

mkdir 1n4ppn -p && cd 1n4ppn
qsub ../1n4ppn.sh

mkdir ../1n6ppn -p && cd ../1n6ppn
qsub ../1n6ppn.sh

mkdir ../1n8ppn -p && cd ../1n8ppn
qsub ../1n8ppn.sh

mkdir ../1n12ppn -p && cd ../1n12ppn
qsub ../1n12ppn.sh

mkdir ../3n12ppn -p && cd ../3n12ppn
qsub ../3n12ppn.sh

mkdir ../4n9ppn -p && cd ../4n9ppn
qsub ../4n9ppn.sh

mkdir ../6n6ppn -p && cd ../6n6ppn
qsub ../6n6ppn.sh

mkdir ../9n4ppn -p && cd ../9n4ppn
qsub ../9n4ppn.sh
