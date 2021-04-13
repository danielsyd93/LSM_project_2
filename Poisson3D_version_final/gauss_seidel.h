/* gauss_seidel.h - Poisson problem
 *
 */
#ifndef _GAUSS_SEIDEL_H
#define _GAUSS_SEIDEL_H

void gauss_seidel_serial(double ***, double ***, int, int, double);
void gauss_seidel_omp(double ***, double ***, int, int, double);
void gauss_seidel_omp_rb(double ***, double ***, int, int, double);
void gauss_seidel_omp_doacross(double ***, double ***, int, int, double);

#endif
