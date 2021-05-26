#include <stdio.h>
#include <mpi.h>
#include "alloc.h"

#define ROWS 9
#define COLS 7 

/* A sample program that demonstrates how to time an MPI program, and how to
 * send (parts of) the rows and columns in a two-dimensional array between
 * two processes using MPI.
 * In the matrix representation used here, columns are simple to copy, as
 * column elements are contiguous in memory. Row copying is more complex as
 * consecutive row elements are not contiguous in memory and require the
 * use of derived MPI_Datatypes to perform the copy.
 */

void zeromatrix(float **matrix);
void printmatrix(char *title, float **matrix);

int main(int argc, char **argv)
{
    double t;
    int i, j, n, p;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);

    if (n != 2) {
        printf("Must be run on exactly 2 processors.\n");
        MPI_Finalize();
        return 1;
    }

    t = MPI_Wtime();

    /* Allocate an COLS * ROWS array. */
    float **matrix = alloc_floatmatrix(COLS, ROWS);

    /* Fill processor 1's matrix with numbers */
    for (i = 0; i < COLS; i++) {
        for ( j = 0; j < ROWS; j++) {
            matrix[i][j] = (i * 10) + j;
        }
    }

    /* Define two MPI_Datatypes for rows that we use later */
    MPI_Datatype fullrowtype, partrowtype;
    MPI_Type_vector(ROWS, 1, ROWS, MPI_FLOAT, &fullrowtype);
    MPI_Type_commit(&fullrowtype);
    MPI_Type_vector(3, 1, ROWS, MPI_FLOAT, &partrowtype);
    MPI_Type_commit(&partrowtype);

    if (p == 0) {
        MPI_Status s;

        zeromatrix(matrix);
        MPI_Recv(matrix[4], ROWS, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &s);
        printmatrix("After receiving column 4:", matrix);

        zeromatrix(matrix);
        MPI_Recv(&matrix[6][2], 4, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &s);
        printmatrix("After receiving column 6, rows 3-5:", matrix);

        zeromatrix(matrix);
        MPI_Recv(matrix[3], ROWS*2, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &s);
        printmatrix("After receiving column 3 and 4:", matrix);

        zeromatrix(matrix);
        MPI_Recv(matrix[0], ROWS*COLS, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &s);
        printmatrix("After receiving all columns:", matrix);

        zeromatrix(matrix);
        MPI_Recv(&matrix[0][6], 1, fullrowtype, 1, 0, MPI_COMM_WORLD, &s);
        printmatrix("After receiving row 6:", matrix);
       
        zeromatrix(matrix); 
        MPI_Recv(&matrix[0][1], 1, partrowtype, 1, 0, MPI_COMM_WORLD, &s);
        printmatrix("After receiving row 1 cols 0-2:", matrix);
        
        zeromatrix(matrix); 
        MPI_Recv(&matrix[4][1], 1, partrowtype, 1, 0, MPI_COMM_WORLD, &s);
        printmatrix("After receiving row 1 cols 4-6:", matrix);
    } else {
        /* Send all of column 4 to processor 0 */
        MPI_Send(matrix[4], ROWS, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

        /* Send column 6 rows 2-5 to processor 0 */
        MPI_Send(&matrix[6][2], 4, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

        /* Send columns 3 and 4 to processor 0 */
        MPI_Send(matrix[3], ROWS*2, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

        /* Send the entire matrix (ie all columns) to processor 0 */
        MPI_Send(matrix[0], ROWS*COLS, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

        /* Send row 6 to processor 0 */
        MPI_Send(&matrix[0][6], 1, fullrowtype, 0, 0, MPI_COMM_WORLD);
        
        /* Send row 1 cols 0-2 to processor 0 */
        MPI_Send(&matrix[0][1], 1, partrowtype, 0, 0, MPI_COMM_WORLD);
        
        /* Send row 1 cols 4-6 to processor 0 */
        MPI_Send(&matrix[4][1], 1, partrowtype, 0, 0, MPI_COMM_WORLD);
    }
    if (p == 0) {
        t = MPI_Wtime() - t;
        printf("Program took %f secs to run.\n", t);
    }

    /* Free the matrix we allocated */
    free_matrix(matrix);

    /* Free the derived MPI_Datatypes */
    MPI_Type_free(&fullrowtype);
    MPI_Type_free(&partrowtype);

    MPI_Finalize();
    return 0;
}

void zeromatrix(float **matrix)
{
    int i, j;
    for ( j = 0; j < ROWS; j++) {
        for ( i = 0; i < COLS; i++) {
            matrix[i][j] = 0.0;
        }
    }
}
void printmatrix(char *title, float **matrix)
{
    int i, j;
    printf("%s\n", title);
    for ( j = 0; j < ROWS; j++) {
        for ( i = 0; i < COLS; i++) {
            printf("%02g ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
