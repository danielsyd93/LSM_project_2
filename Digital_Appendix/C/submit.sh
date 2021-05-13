#!/bin/sh

#BSUB -J gaussSeidel

#BSUB -n 1
#BSUB -R "span[ptile=24]"
#BSUB -W 00:59
#BSUB -q hpcintro
#BSUB -R "rusage[mem=8GB]"
#BSUB -oo Output.out
#BSUB -eo Error.err

module purge
module load studio
module load mpi/3.1.3-gcc-8.2.0

mpicc -o main -O2 -ftree-vectorize -march=native main.c

mpirun ./main 1000 