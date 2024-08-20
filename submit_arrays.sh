#!/bin/bash
## Job name:
#SBATCH --job-name=b_b
## Project:
#SBATCH --account=nn9487k 
## Number of nodes/tasks:
##SBATCH --nodes=11
#SBATCH --ntasks=1
#SBATCH --array=1-10
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=128
## Wall clock limit:
#SBATCH --time=0-10:00:00

set -o errexit
set -o nounset

module restore system
#module load intel/2021a
#module load intel/2023a
module load intel/2022a
module list

OUTFILE=results.$SLURM_ARRAY_TASK_ID


#time mpirun ./MAIN
./MAIN $SLURM_ARRAY_TASK_ID > $OUTFILE

