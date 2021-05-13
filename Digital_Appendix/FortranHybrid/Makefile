target = runme

.SUFFIXES:
.SUFFIXES: .f90 .o .c

CC = mpicc
FC = mpifort
FFLAGS = -O3 -ftree-vectorize -march=native -ffast-math -funroll-loops -fopenmp

OBJS = main.o write_vtk.o

LIBS =

.PHONY: all
all: $(target)

.PHONY: new
new: clean $(target)

.PHONY: clean realclean
clean:
	@/bin/rm -f $(OBJS) *.mod

realclean: clean
	@/bin/rm -f $(target)

# linking: the target depends on the objects
$(target): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(target)

.f90.o:
	$(FC) -c $(FFLAGS) $<

.c.o:
	$(CC) -c $<

# dependencies:
main.o: write_vtk.o
routines.o: precision.o write_vtk.o


