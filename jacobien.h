#ifndef __JACOBIEN__
#define __JACOBIEN__

#include<stdbool.h>

#ifdef __cplusplus
extern "C"{
#endif

void Jakobien_seq(double *u, double *u_new, double *f, int N, int iter_max, double threshold);

void Jakobien_MPI(double *u, double *u_new, double *f, int N, int *Local_N, double threshold, int size, int rank);



#ifdef __cplusplus
} 
#endif 

#endif
