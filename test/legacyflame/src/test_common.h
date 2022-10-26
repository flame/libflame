/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Sep 24, 2020
*/

#define DRAND()  ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;
#define SRAND()  ( float ) ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;

/* Initialize matrix with random values */
void rand_matrix_s( float *A, integer M, integer N, integer LDA );
void rand_matrix_d( double *A, integer M, integer N, integer LDA );
void rand_matrix_c( scomplex *A, integer M, integer N, integer LDA );
void rand_matrix_z( dcomplex *A, integer M, integer N, integer LDA );

/* Initialize symmetric matrix with random values */
void rand_sym_matrix_s( float *A, integer M, integer N, integer LDA );
void rand_sym_matrix_d( double *A, integer M, integer N, integer LDA );
void rand_sym_matrix_c( scomplex *A, integer M, integer N, integer LDA );
void rand_sym_matrix_z( dcomplex *A, integer M, integer N, integer LDA );

/* Copy a matrix */
void copy_matrix_s( float *A, float *B, integer M, integer N, integer LDA, integer LDB );
void copy_matrix_d( double *A, double *B, integer M, integer N, integer LDA, integer LDB );
void copy_matrix_c( scomplex *A, scomplex *B, integer M, integer N, integer LDA, integer LDB );
void copy_matrix_z( dcomplex *A, dcomplex *B, integer M, integer N, integer LDA, integer LDB );

/* Pack a symmetric matrix in column first order */
void pack_matrix_lt_s( float *A, float *B, integer N, integer LDA );
void pack_matrix_lt_d( double *A, double *B, integer N, integer LDA );
void pack_matrix_lt_c( scomplex *A, scomplex *B, integer N, integer LDA );
void pack_matrix_lt_z( dcomplex *A, dcomplex *B, integer N, integer LDA );

/* Initialize a matrix with zeros */
void reset_matrix_s( float *A, integer M, integer N, integer LDA );
void reset_matrix_d( double *A, integer M, integer N, integer LDA );
void reset_matrix_c( scomplex *A, integer M, integer N, integer LDA );
void reset_matrix_z( dcomplex *A, integer M, integer N, integer LDA );

/* Set a matrix to identity */
void set_identity_s( float *A, integer M, integer N, integer LDA );
void set_identity_d( double *A, integer M, integer N, integer LDA );
void set_identity_c( scomplex *A, integer M, integer N, integer LDA );
void set_identity_z( dcomplex *A, integer M, integer N, integer LDA );

/* Division of complex types */
void c_div_t(scomplex *cp, scomplex *ap, scomplex *bp);
void z_div_t(dcomplex *cp, dcomplex *ap, dcomplex *bp);
