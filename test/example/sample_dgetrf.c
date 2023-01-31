/*
    Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*/

/*
    Computes the LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.
*/

#include <stdio.h>
#include <stdlib.h>

typedef int integer;
/* Generate Random Value */
#define DRAND()  ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;

/* External calls to LAPACK API's */
extern int dgetrf_(integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);

/* Local Function Declaration */
void rand_matrix(void *A, integer M, integer N, integer LDA);

int main( int argc, char** argv )
{
    /* Intialize input matrix sizes */
    integer M = 10, N = 10, LDA = 10;
    double *A;
    integer *ipiv;
    integer i, info;

    /* Allocation of memory to matrix*/
    A = (double  *) malloc(LDA * N * sizeof(double));
    ipiv  = (integer *) malloc(N * sizeof(double));

    /* Intialize matrix with random values */
    rand_matrix(A, M, N, LDA);
    
    printf("Started execution of DGETRF API \n");

    /* Call to the DGETRF API */
    dgetrf_(&M, &N, A, &LDA, ipiv, &info);

    if(info == 0)
        printf("DGETRF API execution sucessfully completed\n");
    else
        printf("DGETRF Execution Failed\n");
    return 0;
}

/* Intialize the matrix with random values */
void rand_matrix(void *A, integer M, integer N, integer LDA)
{
    integer i, j;
    for( i = 0; i < N; i++ )
    {
        for( j = 0; j < M; j++ )
        {
            ((double *)A)[i * LDA + j] = DRAND();
        }
    }
}