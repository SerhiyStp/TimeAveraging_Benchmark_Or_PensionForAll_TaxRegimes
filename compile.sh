#!/bin/bash
module restore system
#module load intel-compilers/2021.2.0
module load intel-compilers/2023.2.0
#module load intel-mpi/2021.2.0
module list

make
