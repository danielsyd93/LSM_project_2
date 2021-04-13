/* jacobi.c - Poisson problem in 3d
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define POW2(x) ((x)*(x))

void jacobi_serial(double *** u, double *** u_old, double *** f, int N, int iter_max, double threshold) {
	// init
    int n,i,j,k;
	double ***tmp;
	double difference = -1.0, pow2_gridsize = POW2(1.0/(N-1)), divide=1.0/6.0;
	
	// Iteration loop with early stopping
	for(n=0; n<iter_max; n++) {
		// Difference init
		difference = 0.0;

		// Update Matrix 
		for(k=1; k<N-1; k++)
			for(j=1; j<N-1; j++)
				for(i=1; i<N-1; i++) {
					u_old[k][j][i] = (	u[k-1][j][i] +
										u[k+1][j][i] +
										u[k][j-1][i] +
										u[k][j+1][i] +
										u[k][j][i-1] +
										u[k][j][i+1] +
										pow2_gridsize * f[k][j][i]	) * divide;

					// running norm calculation
					difference += POW2(u[k][j][i] - u_old[k][j][i]);
				}	

		// Difference using L2-norm
		difference = sqrt(difference);
		if(difference<threshold)
			break;

		// Swap pointers
		tmp = u_old;
		u_old = u;	
		u = tmp;	
	}

	printf("Iterations %d - Difference: %f - Threshold: %f\n", n, difference, threshold);
}

void jacobi_omp(double *** u, double *** u_old, double *** f, int N, int iter_max, double threshold) {
	// init
    int n=0,i,j,k;
	double ***tmp;
	double difference = -1.0, pow2_gridsize = POW2(1.0/(N-1)), divide=1.0/6.0;
	
	// Iteration loop with early stopping
	for(n=0; n<iter_max; n++) {
		// Difference init
		difference = 0.0;

		// Update Matrix
		#pragma omp parallel for default(none) shared(n, u, u_old, pow2_gridsize, divide, N, f) private(i,j,k) reduction(+: difference)
		for(i=1; i<N-1; i++) {
			for(j=1; j<N-1; j++){
				for(k=1; k<N-1; k++) {
					u_old[i][j][k] = (	u[i-1][j][k] +
										u[i+1][j][k] +
										u[i][j-1][k] +
										u[i][j+1][k] +
										u[i][j][k-1] +
										u[i][j][k+1] +
										pow2_gridsize * f[i][j][k]	) * divide;

					// running norm calculation
					difference += POW2(u[i][j][k] - u_old[i][j][k]);
				}
			}
		}
	
		// Difference using L2-norm
		difference = sqrt(difference);
		if(difference<threshold)
			break;

		// Swap pointers
		tmp = u_old;
		u_old = u;	
		u = tmp;
	}

	printf("Iterations %d - Difference: %f - Threshold: %f\n", n, difference, threshold);
}

