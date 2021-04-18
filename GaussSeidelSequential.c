#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

int const N = 100;
double matrix[N][N][N];
double delta = (double)1/N;

bool checkForBoundary(int aI, int aJ, int aK)
{
    double x = (double)aI / (N-1);
    double y = (double)aJ / (N-1);
    double z = (double)aK / (N-1);

    if(x == 0 || x == 1 || y == 0 || y == 1 || z == 0 || z == 1)
    {
        return true;
    }

    return false;
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

void evalPoint(int aI, int aJ, int aK)
{
    matrix[aI][aJ][aK] = (double)1/6 * (matrix[aI - 1][aJ][aK] + matrix[aI+1][aJ][aK] + matrix[aI][aJ-1][aK]
        + matrix[aI][aJ+1][aK] + matrix[aI][aJ][aK-1] + matrix[aI][aJ][aK+1] + delta * delta * getRadiatorValue(aI, aJ, aK));
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

void write_vtk(int n, double* u);
void writeParaviewCSV();

int main()
{
    printf("Program start\n");

    initMatrix();

    for(int i = 0; i < 1000; ++i)
    {
        evalSequentually();
    }

    double* matrix1D = matrixTo1D();

    writeParaviewCSV();
    write_vtk(N, matrix1D);

   // printMatrixCommandLine();

    free(matrix1D);

    printMinMax();
}

/* $Id: print.c,v 1.1 2019/12/12 15:03:38 gbarbd Exp gbarbd $ */
#include <stdio.h>
#include <inttypes.h>

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

    if(false) {
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
