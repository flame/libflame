/*
	Copyright (c) 2022 Advanced Micro Devices, Inc.  All rights reserved.
*/


#define DRAND()  ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;
#define SRAND()  ( float ) ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;

#define min( x, y ) ( (x) < (y) ? (x) : (y) )

#define max( x, y ) ( (x) > (y) ? (x) : (y) )

/* matrix functions*/
void create_matrix(integer datatype, integer m, integer n, void **A);
void create_realtype_matrix(integer datatype, integer m, integer n, void **A);
void free_matrix(void *A);

/* Initialize matrix with random values */
void rand_matrix(integer datatype, integer m, integer n, void *A);
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
void copy_matrix(integer datatype, integer m, integer n, void *A, void *B);
void copy_realtype_matrix(integer datatype, integer m, integer n, void *A, void *B);
void copy_matrix_i( integer *A, integer *B, integer M, integer N, integer LDA, integer LDB );
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

/* work value calculation*/
integer get_work_value( integer datatype, void *work );

/* Diagonal Scaling*/
void diagmv( integer datatype, integer m, integer n, void* x, integer incx, void* a, integer a_rs, integer a_cs );
void scalv( integer datatype, integer n, void* x, integer incx, void* y, integer incy );