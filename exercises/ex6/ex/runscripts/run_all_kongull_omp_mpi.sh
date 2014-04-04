cd /home/$USER/tma4280Group/exercises/ex6/ex/

#Compile
./cmake_kongull_OpenMP.sh

cd /home/$USER/tma4280Group/exercises/ex6/ex/runscripts/omp_mpi/
chmod +x p4t3n1.sh
chmod +x p6t2n1.sh
chmod +x p8t4n2.sh
chmod +x p1t12n1.sh
chmod +x p12t1n1.sh
chmod +x p24t2n3.sh

mkdir p4t3n1 -p && cd p4t3n1  && qsub ../p4t3n1.sh
mkdir ../p6t2n1 -p && cd ../p6t2n1 && qsub ../p6t2n1.sh
mkdir ../p8t4n2 -p && cd ../p8t4n2 && qsub ../p8t4n2.sh
mkdir ../p1t12n1 -p && cd ../p1t12n1 && qsub ../p1t12n1.sh
mkdir ../p12t1n1 -p && cd ../p12t1n1 && qsub ../p12t1n1.sh
mkdir ../p24t2n3 -p && cd ../p24t2n3 && qsub ../p24t2n3.sh
