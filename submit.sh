#!/bin/bash
#SBATCH --job-name=TimeAv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
## Wall clock limit:
#SBATCH --time=2-10:00:00

#set -o errexit
#set -o nounset

module restore system
module load intel-compilers/2025.0.0
module list

time ./MAIN > ProgramOutput.txt
#exit 0

