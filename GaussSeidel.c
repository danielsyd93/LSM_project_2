#include <mpi.h>
#include <omp.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

int N = 10;

double*** matrix;
double* matrixWin;
double delta;

struct coordinates{
    int x;
    int y;
    int z;
};

int checkForBoundary(int aI, int aJ, int aK)
{
    double x = (double)aI / (N-1);
    double y = (double)aJ / (N-1);
    double z = (double)aK / (N-1);

    if(x == 0 || x == 1 || y == 0 || y == 1 || z == 0 || z == 1)
    {
        return 1;
    }

    return 0;
}

int getBoundaryValue(int aI, int aJ, int aK)
{
    double x = (double)aI / (N-1);
    double y = (double)aJ / (N-1);
    double z = (double)aK / (N-1);

    if(y == 1)
    {
        return 20;
    }

    if(y == 0)
    {
        return 0;
    }
    if(x == 1 || x == 0)
    {
        return 20;
    }

    if(z == 1 || z == 0)
    {
        return 20;
    }

    printf("Error\n");
    return -1;
}

void initMatrix()
{
    for(int x = 0; x < N; ++x)
    {
        for(int y = 0; y < N; ++y)
        {
            for(int z = 0; z < N; ++z)
            {
                if(checkForBoundary(x, y, z))
                {
                    matrix[x][y][z] = getBoundaryValue(x, y, z);
                }
                else
                {
                    matrix[x][y][z] = 0;
                }
                
            }
        }
    }
}

void initMatrixWin()
{
    for(int x = 0; x < N; ++x)
    {
        for(int y = 0; y < N; ++y)
        {
            for(int z = 0; z < N; ++z)
            {
                if(checkForBoundary(x, y, z))
                {
                    matrixWin[x*N*N+y*N+z] = getBoundaryValue(x, y, z);
                }
                else
                {
                    matrixWin[x*N*N+y*N+z] = 0;
                }
                
            }
        }
    }
}

int getRadiatorValue(int aI, int aJ, int aK)
{
    double x = (double)aI / (N-1);
    double y = (double)aJ / (N-1);
    double z = (double)aK / (N-1);
    if(0 <= x && x <= (double)3/16 && 0 <= y && y <= (double)1/4 && (double)2/6 <= z && z <= (double)1/2)
    {
        return 200;
    }

    return 0;
}

double evalPoint(int aI, int aJ, int aK)
{
    matrix[aI][aJ][aK] = (double)1/6 * (matrix[aI - 1][aJ][aK] + matrix[aI+1][aJ][aK] + matrix[aI][aJ-1][aK]
        + matrix[aI][aJ+1][aK] + matrix[aI][aJ][aK-1] + matrix[aI][aJ][aK+1] + delta * delta * getRadiatorValue(aI, aJ, aK));
    return matrix[aI][aJ][aK];
}

void evalPointWin(int aI, int aJ, int aK)
{
    matrixWin[aI*N*N+aJ*N+aK] = (double)1/6 * (matrixWin[(aI-1)*N*N+aJ*N+aK] + matrixWin[(aI+1)*N*N+aJ*N+aK] + matrixWin[aI*N*N+(aJ-1)*N+aK]
        + matrixWin[aI*N*N+(aJ+1)*N+aK] + matrixWin[aI*N*N+aJ*N+aK-1] + matrixWin[aI*N*N+aJ*N+aK-1] + delta * delta * getRadiatorValue(aI, aJ, aK));
}

void evalSequentually()
{
    for(int x = 1; x < N-1; ++x)
    {
        for(int y = 1; y < N-1; ++y)
        {
            for(int z = 1; z < N-1; ++z)
            {
                evalPoint(x, y, z);
            }
        }
    }
}

double* matrixTo1D()
{
    double* matrix1D = (double*)malloc(sizeof(double) * N * N * N);

    int i = 0;

    for(int x = 0; x < N; ++x)
    {
        for(int y = 0; y < N; ++y)
        {
            for(int z = 0; z < N; ++z)
            {
                matrix1D[i++] = matrix[x][y][z];
            }
        }
    }

    return matrix1D;
}

void printMatrixCommandLine()
{
    for(int x = 0; x < N; ++x)
    {
        for(int y = 0; y < N; ++y)
        {
            for(int z = 0; z < N; ++z)
            {
                printf("%.2f\t", matrix[x][y][z]);
            }
            printf("\n");
        }
        printf("\n\n");
    }
}

void printMinMax()
{
    double min = matrix[0][0][0];
    double max = matrix[0][0][0];

    for(int x = 0; x < N; ++x)
    {
        for(int y = 0; y < N; ++y)
        {
            for(int z = 0; z < N; ++z)
            {
                if(min > matrix[x][y][z])
                {
                    min = matrix[x][y][z];
                }

                if(max < matrix[x][y][z])
                {
                    max = matrix[x][y][z];
                }
            }
        }
    }

    printf("Min: %.2f \t Max: %.2f\n", min, max);
}


int getDiagonal2D(int diagonal, struct coordinates* coordinates)
{
    int count = 0;
    if(diagonal < 3)
    {
        return 0;
    }

    for(int x = 1; x < N-1; ++x)
    {
        for(int y = 1; y < N-1; ++y)
        {
            int z = diagonal-x-y;

            if(0 < z && z < N-1 )
            {
               coordinates[count].x = x;
               coordinates[count].y = y;
               coordinates[count].z = z; 
               count++;
            }
        }
    }
    return count;
}

void wavefrontSeq()
{
    for(int diagonal = 3; diagonal <= (N-2)*3; ++diagonal)
    {
        struct coordinates* coordinates = (struct coordinates*)malloc(sizeof(struct coordinates)*((diagonal-1)*(diagonal-2)/2));

        int count = getDiagonal2D(diagonal,coordinates);
        for(int i = 0; i < count; ++i)
        {
            evalPoint(coordinates[i].x, coordinates[i].y, coordinates[i].z);
        }
        free(coordinates);
    }
}

void wavefrontParallelBcast()
{
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(int diagonal = 3; diagonal <= (N-2)*3; ++diagonal)
    {
        struct coordinates* coordinates = (struct coordinates*)malloc(sizeof(struct coordinates)*((diagonal-1)*(diagonal-2)/2));

        int count = getDiagonal2D(diagonal,coordinates);

        for(int i = 0; i < count; ++i)
        {
            if(i % size == rank)
            {
                double data = evalPoint(coordinates[i].x, coordinates[i].y, coordinates[i].z);
                MPI_Bcast(&data, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
            }
            else
            {
                double data;
                MPI_Bcast(&data, 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
                matrix[coordinates[i].x][coordinates[i].y][coordinates[i].z] = data;
            }
        }

        free(coordinates);
    }
}

void setToZero(int* array, int size)
{
    for(int i = 0; i < size; ++i)
    {
        array[i] = 0;
    }
}

void wavefrontParallel2DBcast()
{
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(int diagonal = 3; diagonal <= (N-2)*3; ++diagonal)
    {
        struct coordinates* coordinates = (struct coordinates*)malloc(sizeof(struct coordinates)*((diagonal-1)*(diagonal-2)/2));

        int count = getDiagonal2D(diagonal,coordinates);

        int sendCount[size];

        setToZero(sendCount, size);

        for(int i = 0; i < count; ++i)
        {
           sendCount[i % size] += 1;
        }

        double** data= (double**)malloc(sizeof(double*)*size);

        for(int i = 0; i < size; ++i)
        {
            data[i] = (double*)malloc(sizeof(double)*(sendCount[i % size]));
        }

        int dataCount = 0;

        for(int i = rank; i < count; i+=size)
        {
            data[rank][dataCount++] = evalPoint(coordinates[i].x, coordinates[i].y, coordinates[i].z);
        }

        for(int i = 0; i < size; ++i)
        {
            if(i % size == rank)
            {
                MPI_Bcast(data[rank], sendCount[rank], MPI_DOUBLE, rank, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Bcast(data[i % size], sendCount[i%size], MPI_DOUBLE, i % size, MPI_COMM_WORLD);
            }
        }

        setToZero(sendCount, size);

        for(int i = 0; i < count; ++i)
        {
            matrix[coordinates[i].x][coordinates[i].y][coordinates[i].z] = data[i % size][sendCount[i%size]++];
        }

        for(int i = 0; i < size; ++i)
        {
            free(data[i]);
        }
        free(data);
        free(coordinates);
    }
}

void wavefrontParallelHybrid()
{
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(int diagonal = 3; diagonal <= (N-2)*3; ++diagonal)
    {
        struct coordinates* coordinates = (struct coordinates*)malloc(sizeof(struct coordinates)*((diagonal-1)*(diagonal-2)/2));

        int count = getDiagonal2D(diagonal,coordinates);

        int sendCount[size];

        setToZero(sendCount, size);

        for(int i = 0; i < count; ++i)
        {
           sendCount[i % size] += 1;
        }

        double** data= (double**)malloc(sizeof(double*)*size);

        for(int i = 0; i < size; ++i)
        {
            data[i] = (double*)malloc(sizeof(double)*(sendCount[i % size]));
        }

        int dataCount = 0;


        #pragma omp parallel for
        for(int i = rank; i < count; i+=size)
        {
            data[rank][i/size] = evalPoint(coordinates[i].x, coordinates[i].y, coordinates[i].z);
        }

        for(int i = 0; i < size; ++i)
        {
            if(i % size == rank)
            {
                MPI_Bcast(data[rank], sendCount[rank], MPI_DOUBLE, rank, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Bcast(data[i % size], sendCount[i%size], MPI_DOUBLE, i % size, MPI_COMM_WORLD);
            }
        }

        setToZero(sendCount, size);

        for(int i = 0; i < count; ++i)
        {
            matrix[coordinates[i].x][coordinates[i].y][coordinates[i].z] = data[i % size][sendCount[i%size]++];
        }

        for(int i = 0; i < size; ++i)
        {
            free(data[i]);
        }
        free(data);
        free(coordinates);
    }
}

void wavefrontParallelWin()
{
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(int diagonal = 3; diagonal <= (N-2)*3; ++diagonal)
    {
        struct coordinates* coordinates = (struct coordinates*)malloc(sizeof(struct coordinates)*((diagonal-1)*(diagonal-2)/2));

        int count = getDiagonal2D(diagonal,coordinates);
        for(int i = rank; i < count; i+=size)
        {
            evalPointWin(coordinates[i].x, coordinates[i].y, coordinates[i].z);
        }
        free(coordinates);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void redDotBlackDotWin()
{
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(int x = 1; x < N-1; ++x)
    {
        for(int y = rank + 1; y < N-1; y+=size)
        {
            int zStart;
            if(x+y % 2 == 0)
            {
                zStart = 1;
            }
            else
            {
                zStart = 2;
            }
            for(int z = zStart; z < N-1; z+=2)
            {
                evalPointWin(x, y, z);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for(int x = 1; x < N-1; ++x)
    {
        for(int y = rank + 1; y < N-1; y+=size)
        {
            int zStart;
            if(x+y % 2 == 0)
            {
                zStart = 2;
            }
            else
            {
                zStart = 1;
            }
            for(int z = zStart; z < N-1; z+=2)
            {
                evalPointWin(x, y, z);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void matrixAllocate(int N)
{
    matrix = (double***)malloc(sizeof(double**)*N);

    for(int i = 0; i < N; ++i)
    {
        matrix[i] = (double**)malloc(sizeof(double*)*N);
        for(int j = 0; j < N; ++j)
        {
            matrix[i][j] = (double*)malloc(sizeof(double)*N);
        }
    }
}

void matrixFree()
{
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            free(matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);
}

void write_vtk(int n, double* u);
void writeParaviewCSV();

int main(int argc, char *argv [])
{
    MPI_Init(&argc, &argv);

    int rank;
    int size;

    N = atoi(argv[1]);
    delta = (double)1/(N-1);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);   

    if(rank == 0)
    {
        printf("Program start, N size: %d,  MPI size: %d, openMP: %d\n", N,  size, omp_get_max_threads()); 
    }

    matrixAllocate(N);

    initMatrix();

    double t0;
    double t;

    if(rank == 0)
    {
        t0 = MPI_Wtime();
    }

    for(int i = 0; i < 100; ++i)
    {
        if(rank == 0)
        {
            evalSequentually();
        }
    }

    if(rank == 0)
    {
        t = MPI_Wtime() - t0;
        printf("Serial took %.4f seconds\n", t);
    }

    initMatrix();

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        t0 = MPI_Wtime();
    }

    for(int i = 0; i < 100; ++i)
    {
        if(rank == 0)
        {
            wavefrontSeq();
        }
    }

    if(rank == 0)
    {
        t = MPI_Wtime() - t0;
        printf("wavefrontSeq took %.4f seconds\n", t);
    }

    initMatrix();

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        t0 = MPI_Wtime();
    }

    for(int i = 0; i < 100; ++i)
    {
        wavefrontParallelBcast();
    }

    if(rank == 0)
    {
        t = MPI_Wtime() - t0;
        printf("wavefrontParallelBcast took %.4f seconds\n", t);
    }

    initMatrix();

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        t0 = MPI_Wtime();
    }

    for(int i = 0; i < 100; ++i)
    {
        wavefrontParallel2DBcast();
    }

    if(rank == 0)
    {
        t = MPI_Wtime() - t0;
        printf("wavefrontParallel2DBcast took %.4f seconds\n", t);
    }


    
    matrixFree();

    MPI_Win win;
    MPI_Aint matrixSize = N*N*N*sizeof(double);

    if (rank == 0)
    {
        MPI_Win_allocate_shared(matrixSize, sizeof(double), MPI_INFO_NULL,
                           MPI_COMM_WORLD, &matrixWin, &win);
    }
    else
    {
        int dispUnit =0;
        MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL,
                              MPI_COMM_WORLD, &matrixWin, &win);
        MPI_Win_shared_query(win, 0, &matrixSize, &dispUnit, &matrixWin);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        initMatrixWin();
        t0 = MPI_Wtime();
    }

    for(int i = 0; i < 100; ++i)
    {
            wavefrontParallelWin();
    }

    if(rank == 0)
    {
        t = MPI_Wtime() - t0;
        printf("wavefrontParallelWin took %.4f seconds\n", t);
    }

    if(rank == 0)
    {
        initMatrixWin();
        t0 = MPI_Wtime();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for(int i = 0; i < 100; ++i)
    {
            redDotBlackDotWin();
    }

    if(rank == 0)
    {
        t = MPI_Wtime() - t0;
        printf("redDotBlackDotWin took %.4f seconds\n\n", t);
    }


    //printMinMax();

    //double* matrix1D = matrixTo1D();
    //write_vtk(N, matrix1D);
    //free(matrix1D);


    MPI_Win_free(&win);
    MPI_Finalize();
}

static int is_little_endian(void) {
    int num = 1;
    return (*((char*)&num) == 1);
}

void write_vtk(int n, double* u) {

    FILE* f_ptr;
    size_t written;
    size_t items = n * n * n;
    size_t i;
    int b;
    unsigned char tmp;

    if((f_ptr = fopen("poisson_u.vtk", "w")) == NULL) {
        perror("No output! fopen()");
        return;
    }

    // Write VTK file header
    fprintf(f_ptr, "# vtk DataFile Version 3.0\n");
    fprintf(f_ptr, "saved from function write_vtk.\n");
    fprintf(f_ptr, "BINARY\n");
    fprintf(f_ptr, "DATASET STRUCTURED_POINTS\n");
    fprintf(f_ptr, "DIMENSIONS %d %d %d\n", n, n, n);
    fprintf(f_ptr, "ORIGIN %d %d %d\n", 0, 0, 0);
    fprintf(f_ptr, "SPACING %d %d %d\n", 1, 1, 1);
    fprintf(f_ptr, "POINT_DATA %lu\n", items);
    fprintf(f_ptr, "SCALARS %s %s 1\n", "gray", "double");
    fprintf(f_ptr, "LOOKUP_TABLE default\n");

    if(is_little_endian()) {
        // System is little endian, so we need to reverse the byte order.
        written = 0;
        for(i = 0; i < items; ++i) {
            uint64_t crnt = *(uint64_t*)(u + i); // Get double as int

            // Reverse byte order and write to file
            crnt = (crnt & 0x00000000FFFFFFFF) << 32 | (crnt & 0xFFFFFFFF00000000) >> 32;
            crnt = (crnt & 0x0000FFFF0000FFFF) << 16 | (crnt & 0xFFFF0000FFFF0000) >> 16;
            crnt = (crnt & 0x00FF00FF00FF00FF) << 8 | (crnt & 0xFF00FF00FF00FF00) >> 8;
            written += fwrite(&crnt, sizeof(uint64_t), 1, f_ptr);
        }
    }
    else {
        // System is big endian, so just dump the data.
        written = fwrite(u, sizeof(double), items, f_ptr);
    }

    if(written != items) {
        fprintf(stderr, "Writing failed:  only %lu of %lu items saved!\n",
            written, items);
    }

    fclose(f_ptr);
}

void writeParaviewCSV()
{
    FILE* fpt;

    fpt = fopen("paraview.csv", "w+");

    fprintf(fpt, "x coord,y coord,z coord, scalar\n");
    for(int x = 0; x < N; ++x)
    {
        for(int y = 0; y < N; ++y)
        {
            for(int z = 0; z < N; ++z)
            {
                fprintf(fpt, "%d,%d,%d,%.2f\n", x, y, z, matrix[x][y][z]);
            }
        }
    }
    fclose(fpt);
}
