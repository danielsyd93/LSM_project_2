/* gauss_seidel.c - Poisson problem in 3d
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define POW2(x) ((x)*(x))

void gauss_seidel_serial(double *** u, double *** f, int N, int iter_max, double threshold) { 
	/// init
    int n,i,j,k;
	double temp, difference = -1.0, pow2_gridsize = POW2(1.0/(N-1)), divide=1.0/6.0;
	
	// Iteration loop with early stopping
	for(n=0; n<iter_max; n++) {
		// Difference init
		difference = 0.0;

		// Update Matrix 
		for(k=1; k<N-1; k++)
			for(j=1; j<N-1; j++)
				for(i=1; i<N-1; i++) {
					// save old u[k][j][i]
					temp = u[k][j][i]; 
					
					// calc new u[k][j][i]
					u[k][j][i] = (		u[k-1][j][i] +
										u[k+1][j][i] +
										u[k][j-1][i] +
										u[k][j+1][i] +
										u[k][j][i-1] +
										u[k][j][i+1] +
										pow2_gridsize * f[k][j][i]	) * divide;

					// running norm calculation
					difference += POW2(u[k][j][i] - temp);
				}
		
		// Difference using L2-norm
		difference = sqrt(difference);
		if(difference<threshold)
			break;
	}
	
	printf("Iterations %d - Difference: %f - Threshold: %f\n", n, difference, threshold);
}

void gauss_seidel_omp(double *** u, double *** f, int N, int iter_max, double threshold) {
	/// init
    int n,i,j,k;
	double temp, difference = -1.0, pow2_gridsize = POW2(1.0/(N-1)), divide=1.0/6.0;
	
	// Iteration loop with early stopping
	for(n=0; n<iter_max; n++) {
		// Difference init
		difference = 0.0;

		// Update Matrix
		#pragma omp parallel for default(none) shared(n, u, pow2_gridsize, divide, N, f) private(i,j,k,temp) reduction(+: difference)
		for(k=1; k<N-1; k++)
			for(j=1; j<N-1; j++)
				for(i=1; i<N-1; i++) {
					// save old u[k][j][i]
					temp = u[k][j][i]; 
					
					// calc new u[k][j][i]
					u[k][j][i] = (		u[k-1][j][i] +
										u[k+1][j][i] +
										u[k][j-1][i] +
										u[k][j+1][i] +
										u[k][j][i-1] +
										u[k][j][i+1] +
										pow2_gridsize * f[k][j][i]	) * divide;

					// running norm calculation
					difference += POW2(u[k][j][i] - temp);
				}
		
		// Difference using L2-norm
		difference = sqrt(difference);
		if(difference<threshold)
			break;
	}
	
	printf("Iterations %d - Difference: %f - Threshold: %f\n", n, difference, threshold);
}

void gauss_seidel_omp_rb(double *** u, double *** f, int N, int iter_max, double threshold) {
	/// init
    int n,i,j,k;
	double temp, difference = -1.0, pow2_gridsize = POW2(1.0/(N-1)), divide=1.0/6.0;
	
	// Iteration loop with early stopping
	for(n=0; n<iter_max; n++) {
		// Difference init
		difference = 0.0;
		
		// Update Matrix Red
		#pragma omp parallel for default(none) shared(n, u, pow2_gridsize, divide, N, f) private(temp, i, j, k) reduction(+:difference)
		for(k=1; k<N-1; k++)
			for(j=1; j<N-1; j++)
				for(i=1; i<N-1; i++)
					if ((k+i+j) % 2 == 0) {
						// save old u[k][j][i]
						temp = u[k][j][i]; 
						
						// calc new u[k][j][i]
						u[k][j][i] = (		u[k-1][j][i] +
											u[k+1][j][i] +
											u[k][j-1][i] +
											u[k][j+1][i] +
											u[k][j][i-1] +
											u[k][j][i+1] +
											pow2_gridsize * f[k][j][i]	) * divide;

						// running norm calculation
						difference += POW2(u[k][j][i] - temp);
					}

		// Update Matrix Black
		#pragma omp parallel for default(none) shared(n, u, pow2_gridsize, divide, N, f) private(temp, i, j, k) reduction(+:difference)
		for(k=1; k<N-1; k++)
			for(j=1; j<N-1; j++)
				for(i=1; i<N-1; i++)
					if ((k+i+j) % 2 == 1) {
						// save old u[k][j][i]
						temp = u[k][j][i]; 
						
						// calc new u[k][j][i]
						u[k][j][i] = (		u[k-1][j][i] +
											u[k+1][j][i] +
											u[k][j-1][i] +
											u[k][j+1][i] +
											u[k][j][i-1] +
											u[k][j][i+1] +
											pow2_gridsize * f[k][j][i]	) * divide;

						// running norm calculation
						difference += POW2(u[k][j][i] - temp);
					}
		
		// Difference using L2-norm
		difference = sqrt(difference);
		if(difference<threshold)
			break;
	}
	
	printf("Iterations %d - Difference: %f - Threshold: %f\n", n, difference, threshold);
}

void gauss_seidel_omp_doacross(double *** u, double *** f, int N, int iter_max, double threshold) {
	/// init
    int n,i,j,k, iter=0, done=0;
	double temp, difference = -1.0, diffout = -1.0, pow2_gridsize = POW2(1.0/(N-1)), divide=1.0/6.0;

	// Iteration loop with early stopping
	#pragma omp parallel default(none) private(n, i, j, k, temp) shared(diffout, iter, done, difference, iter_max, threshold, u, pow2_gridsize, divide, N, f)
	for(n=0; n<iter_max; n++) {
		// Update Matrix
		#pragma omp for ordered(2) schedule(static,1) reduction(+:difference)
		for(k=1; k<N-1; k++)
			for(j=1; j<N-1; j++) {
				#pragma omp ordered depend(sink: k-1, j) depend(sink: k, j-1)			
				for(i=1; i<N-1; i++) {
					// save old u[k][j][i]
					temp = u[k][j][i]; 
					
					// calc new u[k][j][i]
					u[k][j][i] = (		u[k-1][j][i] +
										u[k+1][j][i] +
										u[k][j-1][i] +
										u[k][j+1][i] +
										u[k][j][i-1] +
										u[k][j][i+1] +
										pow2_gridsize * f[k][j][i]	) * divide;

					// running norm calculation
					difference += POW2(temp - u[k][j][i]);				
				}
				#pragma omp ordered depend(source)
			}
		
		//Difference using L2-norm
		#pragma omp master
		{
			difference = sqrt(difference);
			done = (difference<threshold) ? 1 : 0;
			diffout = difference;
			difference = 0.0;
			iter++;
		}

		#pragma omp barrier
		if(done)
			break;
	}

	printf("Iterations %d - Difference: %f - Threshold: %f\n", iter, diffout, threshold);
}

