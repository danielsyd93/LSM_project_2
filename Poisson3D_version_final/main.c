/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "alloc3d.h"
#include "print.h"
#include <math.h>

#ifdef _JACOBI
#include "jacobi.h"
#endif

#if defined _GAUSS_SEIDEL || defined _GAUSS_SEIDEL_RB || defined _GAUSS_SEIDEL_DA
#include "gauss_seidel.h"
#endif

#define N_DEFAULT 100

int
main(int argc, char *argv[]) {

    int 	N = N_DEFAULT;
    int 	iter_max = 1000;
    double	tolerance;
    double	start_T;
    int		output_type = 0;
    char	*output_prefix = "poisson_res";
    char        *output_ext    = "";
    char	output_filename[FILENAME_MAX];
    double 	***u = NULL;


    /* get the paramters from the command line */
    N         = atoi(argv[1]);	// grid size
    iter_max  = atoi(argv[2]);  // max. no. of iterations
    tolerance = atof(argv[3]);  // tolerance
    start_T   = atof(argv[4]);  // start T for all inner grid points
    if (argc == 6) {
	output_type = atoi(argv[5]);  // ouput type
    }

    // allocate memory
    if ( (u = d_malloc_3d(N, N, N)) == NULL ) {
        perror("array u: allocation failed");
        exit(-1);
    }

	// --------------------
	// OUR CODE STARTS HERE
	// --------------------

	// Initilize the boundary conditions

	int i, j, k;
    // Zero fill u
	for (k = 0; k < N; k++)
        for (j = 0; j < N; j++)
            for (i = 0; i < N; i++)
				u[k][j][i] = start_T;
  
    // u(1, y, z) = u(-1, y, z) = 20
    for (k = 0; k < N; k++){
        for (j = 0; j < N; j++){
            u[k][j][N-1] = 20;
            u[k][j][0] = 20;
        }
    }
    // u(x, 1, z) = 20, u(x, -1, z) = 0
    for (k = 0; k < N; k++){
        for (i = 0; i < N; i++){
            u[k][N-1][i] = 20;
            u[k][0][i] = 0;
        }
    }
    // u(x, y, 1) = u(x, y, -1) = 20
    for (j= 0; j < N; j++){
        for (i = 0; i < N; i++){
            u[N-1][j][i] = 20;
            u[0][j][i] = 20;
        }
    }
    // Initilize the f-matrix
    int k_s, k_max;
    int j_max, i_max; 

    // Allocate memory for f matrix
    double ***f = d_malloc_3d(N, N, N);

    // Define the indices 
    k_s = floor((1./3.) * N);
    k_max = floor(N/2.);
    j_max = floor((1./4.) * N);
    i_max = floor((5.*N)/8.);

    for (k = 0; k < N; k++){
        for (j = 0; j < N; j++){
            for (i = 0; i < N; i++){
                if (i <= i_max && j <= j_max && k >= k_s && k <= k_max){
                    f[k][j][i] = 200;
                } 
                else
                {
                    f[k][j][i] = 0;
                }               
            }
        }
    }
    
	printf("Start ");
	#ifdef _JACOBI
	printf("- Jacobi:\n");
	// Allocate u_old matrix and copy u into u_old
	double *** u_old = d_malloc_3d(N,N,N);
	for (k = 0; k < N; k++)
        for (j = 0; j < N; j++)
            for (i = 0; i < N; i++)
				u_old[k][j][i] = u[k][j][i];
	
	// Run jacobi
	jacobi_omp(u, u_old, f, N, iter_max, tolerance);

	// free extra U
	free(u_old);
	#endif

	#ifdef _GAUSS_SEIDEL
	printf("- Gauss Seidel:\n");
	// Run gauss seidel
	gauss_seidel_omp(u, f, N, iter_max, tolerance);
	#endif

	#ifdef _GAUSS_SEIDEL_RB
	printf("- Gauss Seidel Red Black:\n");
	// Run gauss seidel red black
	gauss_seidel_omp_rb(u, f, N, iter_max, tolerance);
	#endif

	#ifdef _GAUSS_SEIDEL_DA
	printf("- Gauss Seidel Doacross:\n");
	// Run gauss seidel
	gauss_seidel_omp_doacross(u, f, N, iter_max, tolerance);
	#endif
	
	printf("Done!\n");

	free(f);

	// ------------------
	// OUR CODE ENDS HERE
	// ------------------

    // dump  results if wanted 
    switch(output_type) {
	case 0:
	    // no output at all
	    break;
	case 3:
	    output_ext = ".bin";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write binary dump to %s: ", output_filename);
	    print_binary(output_filename, N, u);
	    break;
	case 4:
	    output_ext = ".vtk";
	    sprintf(output_filename, "%s_%d%s", output_prefix, N, output_ext);
	    fprintf(stderr, "Write VTK file to %s: \n", output_filename);
	    print_vtk(output_filename, N, u);
	    break;
	default:
	    fprintf(stderr, "Non-supported output type!\n");
	    break;
    }

    // de-allocate memory
    free(u);

    return(0);
}
