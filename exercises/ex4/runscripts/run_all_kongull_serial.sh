cd /home/$USER/tma4280Group/exercises/ex4/runscripts/
#Create CMAKE sub-folder
mkdir ../release_serial -p

#Delete whatever was there before
rm -rf ../release_serial/*

cd ../release_serial/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
make clean && make
chmod +x ex4

cd ../runscripts/serial/
chmod +x kongull_serial_job.sh

mkdir 1n1ppn -p && cd 1n1ppn
qsub ../kongull_serial_job.sh

