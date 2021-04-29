# Run command before:
# chmod u+x runme.sh
# Parallel Code:
clear
echo "Running parallel script!"
mpifort main.f90 -O2 -ftree-vectorize -march=native -o runme

mpirun -np 8 ./runme


#mpifort ex1.f90 -O2 -ftree-vectorize -march=native -o runme

