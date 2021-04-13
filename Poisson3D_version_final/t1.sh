#!/bin/bash

#BSUB -J t1-test
#BSUB -o output_t1/j_%J.out
#BSUB -q hpcintro
#BSUB -n 24
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2048]"
#BSUB -W 01:00

module load studio
module load clang/9.0.0

./runbatch_j_t1.sh
./runbatch_gs_t1.sh
