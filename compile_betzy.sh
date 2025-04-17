#!/bin/bash
module restore system
#module load intel/2021a
module load intel/2022a
#module load intel/2022b
#module load intel-compilers/2023.1.0
#module load intel/2023a
#module load OpenMPI/4.1.1-intel-compilers-2021.2.0
#module load OpenMPI/3.1.3-GCC-8.2.0-2.31.1
#module load foss/2019a
#module load intel/2018b

#module load intel-compilers/2021.2.0
#module load intel-mpi/2021.2.0
module list

#mpif90 main_v2.f90 -o fox_v2.exe
make
