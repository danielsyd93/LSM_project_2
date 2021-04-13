EXECUTABLE=poisson_gs_rb

# define the N values for the size of the cube.
#

N="40 80 120 160 200"

# uncomment and set a reasonable BLKSIZE for the blk version
#
TRESH=0.1

# define the max no. of iterations the driver should use - adjust to
# get a reasonable run time.  You can get an estimate by trying this
# on the command line, i.e. "MFLOPS_MAX_IT=10 ./matmult_...." for the
# problem size you want to analyze.
#
MAX_ITER=9999999

# experiment name 
#

# uncomment the HWCOUNT line, if you want to use hardware counters
# define an option string for the harwdware counters (see output of
# 'collect -h' for valid values.  The format is:
# -h cnt1,on,cnt2,on,...  (up to four counters at a time)
#
# the example below is for L1 hits, L1 misses, L2 hits, L2 misses
#
threads="2 4 8 12 16 20"
for n in $N
do
	for t in $threads
	do
		echo "N size = ${n}, Threads = ${t}"
		
		# start the collect command with the above settings
		OMP_NUM_THREADS=${t} time -p ./$EXECUTABLE $n $MAX_ITER $TRESH 0
	done
done
