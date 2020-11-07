/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Sep 24, 2020
*/

#define DRAND()  ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;
#define SRAND()  ( float ) ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;

/* Initialize matrix with random values */
void rand_matrix_s( float *A, int M, int N, int LDA );
void rand_matrix_d( double *A, int M, int N, int LDA );
void rand_matrix_c( scomplex *A, int M, int N, int LDA );
void rand_matrix_z( dcomplex *A, int M, int N, int LDA );

/* Initialize symmetric matrix with random values */
void rand_sym_matrix_s( float *A, int M, int N, int LDA );
void rand_sym_matrix_d( double *A, int M, int N, int LDA );
void rand_sym_matrix_c( scomplex *A, int M, int N, int LDA );
void rand_sym_matrix_z( dcomplex *A, int M, int N, int LDA );

/* Copy a matrix */
void copy_matrix_s( float *A, float *B, int M, int N, int LDA, int LDB );
void copy_matrix_d( double *A, double *B, int M, int N, int LDA, int LDB );
void copy_matrix_c( scomplex *A, scomplex *B, int M, int N, int LDA, int LDB );
void copy_matrix_z( dcomplex *A, dcomplex *B, int M, int N, int LDA, int LDB );

/* Pack a symmetric matrix in column first order */
void pack_matrix_lt_s( float *A, float *B, int N, int LDA );
void pack_matrix_lt_d( double *A, double *B, int N, int LDA );
void pack_matrix_lt_c( scomplex *A, scomplex *B, int N, int LDA );
void pack_matrix_lt_z( dcomplex *A, dcomplex *B, int N, int LDA );

/* Initialize a matrix with zeros */
void reset_matrix_s( float *A, int M, int N, int LDA );
void reset_matrix_d( double *A, int M, int N, int LDA );
void reset_matrix_c( scomplex *A, int M, int N, int LDA );
void reset_matrix_z( dcomplex *A, int M, int N, int LDA );

/* Set a matrix to identity */
void set_identity_s( float *A, int M, int N, int LDA );
void set_identity_d( double *A, int M, int N, int LDA );
void set_identity_c( scomplex *A, int M, int N, int LDA );
void set_identity_z( dcomplex *A, int M, int N, int LDA );

/* Division of complex types */
void c_div_t(scomplex *cp, scomplex *ap, scomplex *bp);
void z_div_t(dcomplex *cp, dcomplex *ap, dcomplex *bp);
