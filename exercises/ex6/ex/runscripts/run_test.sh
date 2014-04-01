cd /home/$USER/tma4280Group/exercises/ex6/ex/

#Compile
./cmake_kongull_OpenMP.sh

cd /home/$USER/tma4280Group/exercises/ex6/ex/runscripts/

mkdir combined_n_equal_4 -p && cd combined_n_equal_4
chmod +x testJob.sh

mkdir outpt/ -p && cd outpt && qsub ../testJob.sh
