#!/bin/bash
## Job name:
#SBATCH --job-name=bench
## Project:
## Number of nodes/tasks:
##SBATCH --nodes=3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=128
##SBATCH --cpus-per-task=10
#SBATCH --cpus-per-task=40
## Wall clock limit:
#SBATCH --time=0-20:00:00
##SBATCH --partition=scavenger
#SBATCH --partition=batch

#set -o errexit
#set -o nounset

module restore system
#module --quiet purge  # Clear any inherited modules
#module load intel/2021a
#module load intel-compilers/2021.2.0
module load intel-compilers/2023.2.0
#module load intel-mpi/2021.2.0
module list


#time mpirun ./MAIN > ProgramOutput.txt
time ./MAIN > ProgramOutput.txt
exit 0

