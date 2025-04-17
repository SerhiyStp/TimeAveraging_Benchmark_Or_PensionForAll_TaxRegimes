#!/bin/bash
## Job name:
#SBATCH --job-name=TimeAv
## Project:
#SBATCH --account=...
## Number of nodes/tasks:
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
##SBATCH --cpus-per-task=10
## Wall clock limit:
#SBATCH --time=0-10:00:00

#set -o errexit
#set -o nounset

module restore system
module load intel/2022a
module list


time mpirun ./MAIN > ProgramOutput.txt
#exit 0

