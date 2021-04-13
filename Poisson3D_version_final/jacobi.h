/* jacobi.h - Poisson problem 
 *
 */

#ifndef _JACOBI_H
#define _JACOBI_H

void jacobi_serial(double ***, double ***, double ***, int, int, double);
void jacobi_omp(double ***, double ***, double ***, int, int, double);
void jacobi_omp_v2(double ***, double ***, double ***, int, int, double);

#endif
