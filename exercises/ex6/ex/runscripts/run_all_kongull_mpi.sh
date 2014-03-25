cd /home/$USER/tma4280Group/exercises/ex6/ex/

#Compile
./cmake_kongull.sh

cd /home/$USER/tma4280Group/exercises/ex6/ex/runscripts/mpi/
chmod +x p4.sh
chmod +x p6.sh
chmod +x p8.sh
chmod +x p36n4.sh
chmod +x p36n6.sh
chmod +x p36n9.sh
chmod +x p36n3.sh
chmod +x p12.sh

mkdir p4 -p && cd p4 && qsub ../p4.sh
mkdir ../p6 -p && cd ../p6 && qsub ../p6.sh
mkdir ../p8 -p && cd ../p8 && qsub ../p8.sh
mkdir ../p12 -p && cd ../p12 && qsub ../p12.sh
mkdir ../p36n4 -p && cd ../p36n4 && qsub ../p36n4.sh
mkdir ../p36n6 -p && cd ../p36n6 && qsub ../p36n6.sh
mkdir ../p36n9 -p && cd ../p36n9 && qsub ../p36n9.sh
mkdir ../p36n3 -p && cd ../p36n3 && qsub ../p36n3.sh
