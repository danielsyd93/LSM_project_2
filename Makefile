defult: Project

MPICC = mpicc
CFLAGS =
INCLUDES = -I.
OBJS=jacobien.o
LIBS = -lm
_DIST_HEADERS = jacobien.h

jacobien.o: jacobien.c $(_DIST_HEADERS)
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)
jacobi_main.o: jacobi_main.c $(_DIST_HEADERS)
	$(MPICC) -c $< $(INCLUDES) $(CFLAGS)

Project: jacobi_main.o $(OBJS) 
	$(MPICC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o Project