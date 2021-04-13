TARGET_J  		= poisson_j		# Jacobi
TARGET_GS 		= poisson_gs	# Gauss-Seidel
TARGET_GS_RB	= poisson_gs_rb	# Gauss-Seidel Red Black 
TARGET_GS_DA	= poisson_gs_da	# Gauss-Seidel Doacross

SOURCES		= main.c print.c alloc3d.c
OBJECTS		= print.o alloc3d.o
MAIN_J		= main_j.o
MAIN_GS 	= main_gs.o
MAIN_GS_RB	= main_gs_rb.o
MAIN_GS_DA	= main_gs_da.o
OBJS_J		= $(MAIN_J) jacobi.o
OBJS_GS		= $(MAIN_GS) gauss_seidel.o
OBJS_GS_RB	= $(MAIN_GS_RB) gauss_seidel.o
OBJS_GS_DA	= $(MAIN_GS_DA) gauss_seidel.o

# options and settings for the SUNCC compiler
CCSUN	= suncc
CC		= clang
DEFS	= 
OPTSUN	= -g -fast -xopenmp -xloopinfo -xvpara
OPT		= -g -Ofast -fopenmp 
IPO		=
ISA		= 
CHIP	= 
ARCH	= 
PARA	= 
DBUG	= -Wall -Wextra
CFLAGS	= $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) $(DBUG) $(XOPTS) 
LDFLAGS = -lm  

all: $(TARGET_J) $(TARGET_GS) $(TARGET_GS_RB) $(TARGET_GS_DA)


$(TARGET_J): $(OBJECTS) $(OBJS_J)
	$(CC) -o $@ $(CFLAGS) $(OBJS_J) $(OBJECTS) $(LDFLAGS)

$(TARGET_GS): $(OBJECTS) $(OBJS_GS)
	$(CC) -o $@ $(CFLAGS) $(OBJS_GS) $(OBJECTS) $(LDFLAGS)

$(TARGET_GS_RB): $(OBJECTS) $(OBJS_GS_RB)	
	$(CC) -o $@ $(CFLAGS) $(OBJS_GS_RB) $(OBJECTS) $(LDFLAGS)

$(TARGET_GS_DA): $(OBJECTS) $(OBJS_GS_DA)	
	$(CC) -o $@ $(CFLAGS) $(OBJS_GS_DA) $(OBJECTS) $(LDFLAGS)


$(MAIN_J):
	$(CC) -o $@ -D_JACOBI $(CFLAGS) -c main.c 

$(MAIN_GS):
	$(CC) -o $@ -D_GAUSS_SEIDEL $(CFLAGS) -c main.c 

$(MAIN_GS_RB):
	$(CC) -o $@ -D_GAUSS_SEIDEL_RB $(CFLAGS) -c main.c 

$(MAIN_GS_DA):
	$(CC) -o $@ -D_GAUSS_SEIDEL_DA $(CFLAGS) -c main.c 

clean:
	@/bin/rm -f core *.o *~

realclean: clean
	@/bin/rm -f $(TARGET_J) $(TARGET_GS) $(TARGET_GS_RB) $(TARGET_GS_DA)

# DO NOT DELETE
main_j.o: main.c print.h jacobi.h
main_gs.o: main.c print.h gauss_seidel.h
main_gs_rb.o: main.c print.h gauss_seidel.h
main_gs_da.o: main.c print.h gauss_seidel.h
print.o: print.h
