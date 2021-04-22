#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<mpi.h>
#include<inttypes.h>

#include "jacobien.h"

void print_vtk(const char *fname, int n, double *U);
void init_data(int N, double *F,double *U, double *U_new);

static int is_little_endian(void) ;

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);
    int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Request req;
    int N, *local_N;
    N=atoi(argv[1]);
    local_N=malloc(size*sizeof(int));
    if(local_N==NULL){
        MPI_Finalize();
        printf("error didt work\n");
        return -1;
    }


    if(size==1){
        local_N[0]=N;
    }
    else{
        int n;
        for(int i=0;i<size;i++){
            if(i==0){
                n=floor(N/size)+1;
            }
            else if(0<i && i<(size-1)){
                n=floor(N/size)+2;
            }
            else{
                n=N-(size-1)*floor(N/size)+1;
            }
            local_N[i]=n;
        }

    }
    double *U, *U_new, *F;
    U=malloc(N*N*N*sizeof(double));
    U_new=malloc(N*N*N*sizeof(double));
    F=malloc(N*N*N*sizeof(double));
    if(U==NULL||U_new==NULL||F==NULL){
        MPI_Finalize();
        printf("error didt work\n");
        return -1;
    }

    init_data(N, F, U, U_new);

  
    
    double *send_buffer=malloc(local_N[rank]*local_N[rank]*sizeof(double));
    double *recive_buffer=malloc(local_N[rank]*local_N[rank]*sizeof(double));

    int N_buffer=local_N[rank]*local_N[rank];
    for(int iter=0; iter<1000;iter++){
        Jakobien_MPI(U, U_new, F, N, local_N, 0.01,size,rank);
        double *U_boarder_send, *U_boarder_recive;
        if(rank==0){
            U_boarder_send=&U[(local_N[rank]-2) * local_N[rank] * local_N[rank]+0 * local_N[rank]+0];
            U_boarder_recive=&U[(local_N[rank]-1) * local_N[rank] * local_N[rank]+0 * local_N[rank]+0];
            memcpy(send_buffer,U_boarder_send,N_buffer*sizeof(double));
        }
        else{
            U_boarder_send=&U[1 * local_N[rank] * local_N[rank]+0 * local_N[rank]+0];
            U_boarder_recive=&U[0 * local_N[rank] * local_N[rank]+0 * local_N[rank]+0];
            memcpy(send_buffer,U_boarder_send,N_buffer*sizeof(double));
        }
        int destination;
        if(rank==0){
            destination=1;
        }
        else{
            destination=0;
        }
        MPI_Isend(send_buffer,N_buffer,MPI_DOUBLE,destination,iter,MPI_COMM_WORLD,&req);
        MPI_Recv(recive_buffer,N_buffer,MPI_DOUBLE,destination,iter,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        MPI_Barrier(MPI_COMM_WORLD);
        memcpy(U_boarder_recive,recive_buffer,N_buffer*sizeof(double));

    }
    //MPI_Barrier(MPI_COMM_WORLD);
    print_vtk("test.vtk", N, U);
    free(send_buffer);
    free(recive_buffer);

    free(U);
    free(U_new);
    free(F);
    MPI_Finalize();
    
    
    
    return 0;
}

void init_data(int N, double *F, double *U, double *U_new){
      double h=1/(N-1);
    for(int i=0;i<N*N;i++){
        F[i]=0;
        U[i]=0;
        U_new[i]=0;
    }

    for(int z=lround(1/3*N);z<lround(N/2);z++){
        for(int y=0; y<lround(1/4*N);y++){
            for(int x=0;x<lround(5/8*N);x++){
                F[z*N*N+y*N+x]=200;
            }
        }
    }
    int N_end=N-1;
    for(int z=0;z<N;z++){
        for(int y=0;y<N;y++){
            U[z*N*N+y*N+0]=20;
            U[z*N*N+y*N+N_end]=20;
            U_new[z*N*N+y*N+0]=20;
            U_new[z*N*N+y*N+N_end]=20;
        }
    }
    for(int z=0;z<N;z++){
        for(int x=0;x<N;x++){
            U[z*N*N+0*N+x]=0;
            U[z*N*N+N_end*N+x]=20;
            U_new[z*N*N+0*N+x]=0;
            U_new[z*N*N+N_end*N+x]=20;
        }
    }

    for(int y=0;y<N;y++){
        for(int x=0;x<N;x++){
            U[0*N*N+y*N+x]=20;
            U[N_end*N*N+y*N+x]=20;
            U_new[0*N*N+y*N+x]=20;
            U_new[N_end*N*N+y*N+x]=20;
        }
    }
}


static int is_little_endian(void) {
    int num = 1;
    return (*((char *)&num) == 1);
}

void print_vtk(const char *fname, int n, double *U) {

   double ***u = malloc(n * sizeof(double **) +
                               n * n * sizeof(double *) +
                               n * n * n * sizeof(double));
    if (u == NULL) {
        return;
    }

    for(int i = 0; i < n; i++) {
        u[i] = (double **) u + n + i * n ;
    }

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            u[i][j] = (double *) u + n + n * n + i * n * n + n * n;
        }
    }

    for(int z=0; z<n;z++){
        for(int y=0; y<n;y++){
            for(int x=0; x<n;x++){
                u[z][y][x]=U[z*n*n+y*n+x];
            }
        }
    }



    FILE *f_ptr;
    size_t written;
    size_t items = n * n * n;
    size_t i;

    if ( (f_ptr = fopen(fname, "w")) == NULL ) {
       perror("No output! fopen()");
       return;
    }

    // Write VTK file header
    fprintf(f_ptr, "# vtk DataFile Version 3.0\n");
    fprintf(f_ptr, "saved from function print_vtk.\n");
    fprintf(f_ptr, "BINARY\n");
    fprintf(f_ptr, "DATASET STRUCTURED_POINTS\n");
    fprintf(f_ptr, "DIMENSIONS %d %d %d\n", n, n, n);
    fprintf(f_ptr, "ORIGIN %d %d %d\n", 0, 0, 0);
    fprintf(f_ptr, "SPACING %d %d %d\n", 1, 1, 1);
    fprintf(f_ptr, "POINT_DATA %lu\n", items);
    fprintf(f_ptr, "SCALARS %s %s 1\n", "gray", "double");
    fprintf(f_ptr, "LOOKUP_TABLE default\n");

    if ( is_little_endian() ) {
        // System is little endian, so we need to reverse the byte order.
        written = 0;
        for (i = 0; i < items; ++i) {
            uint64_t crnt = *(uint64_t *)(u[0][0] + i); // Get double as int

            // Reverse byte order and write to file
            crnt = (crnt & 0x00000000FFFFFFFF) << 32 | (crnt & 0xFFFFFFFF00000000) >> 32;
            crnt = (crnt & 0x0000FFFF0000FFFF) << 16 | (crnt & 0xFFFF0000FFFF0000) >> 16;
            crnt = (crnt & 0x00FF00FF00FF00FF) << 8  | (crnt & 0xFF00FF00FF00FF00) >> 8;
            written += fwrite(&crnt, sizeof(uint64_t), 1, f_ptr);
        }
    } else {
        // System is big endian, so just dump the data.
        written = fwrite(u[0][0], sizeof(double), items, f_ptr);
    }

    if ( written != items ) {
	    fprintf(stderr, "Writing failed:  only %lu of %lu items saved!\n",
		written, items);
    }

    fclose(f_ptr);
    free(u);
}