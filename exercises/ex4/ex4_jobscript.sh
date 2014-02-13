#!/bin/sh
#######################################################################################33
#			myjobscript.sh
#############################################
#
#	This is an example script example.sh

# These commands set up the environment for your job:
# "###PBS" should be replaced with "#PBS" if you want to run that specific PBS commandline, otherwise it will remain "commented out"


# Name of the job
###PBS -N myjob
#PBS -N ex4_debug_no_OpenMP

# Account to run job (change the account name: See https://www.hpc.ntnu.no/display/hpc/Kongull+Batch+Scheduling)
###PBS -A freecycle

# max walltime (must not exceed que max walltime)
###PBS -l walltime=00:04:00
#PBS -l walltime=00:08:00

# Specify resources number of nodes:cores per node
###PBS -l nodes=1:ppn=1
#PBS -l nodes=4:ppn=12:default

# Specify queue to submit to: default, bigmem, express or default
###PBS -q express
#PBS -q optimist

# Specify an email address to be notified on begin, abort, and end
#PBS -M chrischa@stud.ntnu.no

# Ask Arne Morten what this below line is for?
###PBS -j oe

###### !!!!!make sure this is correct before you run this script!!!! ######
# Cd to work directory
cd /home/chrischa/Documents/tma4280Group/exercises/ex4/

module load intelcomp/13.0.1
module load openmpi/1.4.3-intel #Should we load this module when we just test compile too?
KMP_AFFINITY="granularity=fine, compact"


#For-loop example in bash
#for i in $(seq x y);
#	do
#		<bash command1>;<bash command2>;etc;
#	done

for i in $(seq 3 14);
	do
		echo $i #There will be an empty newline between these two commands (AKA "\n\n")
		OMP_NUM_THREADS=3 mpirun -npernode 12 ex4 $1
	done