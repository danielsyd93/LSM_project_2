#include<stdlib.h>
#include<stdio.h>
#include<math.h>
//#include<mpi.h>

#include "jacobien.h"

void Jakobien_seq(double *u, double *u_new, double *f, int N, int iter_max, double threshold){

    double *temp;
    double div, gridsize, diff;
    div=1.0/6.0;
    gridsize=1/(N-1)*1/(N-1);
    int iter, i, j, k;
    for(iter=0;iter<iter_max;iter++){
        diff=0.0;
        for(i=1;i<N;i++){
            for(j=1;j<N;j++){
                for(k=1;k<N;k++){
                    u_new[i*N*N+j*N+k]=(u[(i-1)*N*N+j*N+k] + u[(i+1)*N*N+j*N+k] +
                                    u[i*N*N+(j-1)*N+k] + u[i*N*N+(j+1)*N+k] +
                                    u[i*N*N+j*N+(k-1)] + u[i*N*N+j*N+(k+1)] + 
                                    gridsize*f[i*N*N+j*N+k])*div;

                    diff+=(u[i*N*N+j*N+k]-u_new[i*N*N+j*N+k])*(u[i*N*N+j*N+k]-u_new[i*N*N+j*N+k]);
                }
            }
        }
    diff=sqrt(diff);
    if( diff<threshold){
        break;
    }
    temp=u_new;
    u_new=u;
    u=temp;
    }
}
void Jakobien_MPI(double *u, double *u_new, double *f, int N, int *Local_N, double threshold, int size, int rank){


    int n=Local_N[rank];
    double *temp;
    double div, gridsize, diff;
    div=1.0/6.0;
    gridsize=1/(N-1)*1/(N-1);

     
        diff=0.0;
        int z, y, x;
        for(z=1;z<n;z++){
            for(y=1;y<n;y++){
                for(x=1;x<n;x++){
                    u_new[z*n*n+y*n+x]=(u[(z-1)*n*n+y*n+x] + u[(z+1)*n*n+y*n+x] +
                                    u[z*n*n+(y-1)*n+x] + u[z*n*n+(y+1)*n+x] +
                                    u[z*n*n+y*N+(x-1)] + u[z*n*n+y*n+(x+1)] + 
                                    gridsize*f[z*n*n+y*n+x])*div;

                    diff+=(u[z*n*n+y*n+x]-u_new[z*n*n+y*n+x])*(u[z*n*n+y*n+x]-u_new[z*n*n+y*n+x]);
                }
            }
        }
    diff=sqrt(diff);
    
    temp=u_new;
    u_new=u;
    u=temp;
    
}