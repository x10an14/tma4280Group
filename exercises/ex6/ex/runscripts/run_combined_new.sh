cd /home/$USER/tma4280Group/exercises/ex6/ex/

#Compile
./cmake_kongull_OpenMP.sh

cd /home/$USER/tma4280Group/exercises/ex6/ex/runscripts/omp_mpi/
chmod +x p4t3n3.sh
chmod +x p6t2n3.sh
chmod +x p1t12n3.sh

chmod +x p3t3n4.sh

chmod +x p3t2n6.sh
chmod +x p2t3n6.sh

#3 nodes
mkdir p4t3n3 -p && cd p4t3n3  && qsub ../p4t3n3.sh
mkdir ../p6t2n3 -p && cd ../p6t2n3 && qsub ../p6t2n3.sh
mkdir ../p1t12n3 -p && cd ../p1t12n3 && qsub ../p1t12n3.sh

#4 nodes
mkdir ../p3t3n4 -p && cd ../p3t3n4 && qsub ../p3t3n4.sh

#6 nodes
mkdir ../p3t2n6 -p && cd ../p3t2n6 && qsub ../p3t2n6.sh
mkdir ../p2t3n6 -p && cd ../p2t3n6 && qsub ../p2t3n6.sh
