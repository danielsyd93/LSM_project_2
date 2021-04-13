#!/bin/bash

#BSUB -J gauss-test
#BSUB -o output_gs/gs_%J.out
#BSUB -q hpcintro
#BSUB -n 24
#BSUB -R "rusage[mem=2048]"
#BSUB -R "span[hosts=1]"
#BSUB -W 01:00
module load clang/9.0.0
module load studio

./runbatch_gs.sh

