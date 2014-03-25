cd /home/$USER/tma4280Group/exercises/ex6/ex/

#Compile
./cmake_kongull_OpenMP.sh

cd /home/$USER/tma4280Group/exercises/ex6/ex/runscripts/omp_mpi/
chmod +x 1n4ppn.sh
chmod +x 1n6ppn.sh
chmod +x 2n4ppn.sh
chmod +x 2n8ppn.sh
chmod +x 1n12ppn.sh

mkdir 1n4ppn -p && cd 1n4ppn  && qsub ../1n4ppn.sh
mkdir ../1n6ppn -p && cd ../1n6ppn && qsub ../1n6ppn.sh
mkdir ../2n4ppn -p && cd ../2n4ppn && qsub ../2n4ppn.sh
mkdir ../2n8ppn -p && cd ../2n8ppn && qsub ../2n8ppn.sh
mkdir ../1n12ppn -p && cd ../1n12ppn && qsub ../1n12ppn.sh
