/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file test_common.h
 *  @brief Defines function declarations to use in APIs of test suite.
 *  */

#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifndef DATATYPES
#define DATATYPES

#if defined(FLA_ENABLE_ILP64)
typedef int64_t integer;
typedef uint64_t uinteger;
#else
typedef int integer;
typedef unsigned int  uinteger;
#endif

typedef struct scomplex
{
    float real, imag;
} scomplex;

typedef struct dcomplex
{
    double real, imag;
} dcomplex;

#endif  // DATATYPES

#include "test_prototype.h"
#include "validate_common.h"

// global variables
extern integer i_zero , i_one , i_n_one;
extern float s_zero, s_one, s_n_one;
extern double d_zero, d_one, d_n_one;
extern scomplex c_zero, c_one, c_n_one;
extern dcomplex z_zero, z_one, z_n_one;

#define DRAND()  ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;
#define SRAND()  ( float ) ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;
#undef min
#undef max
#define min( x, y ) ( (x) < (y) ? (x) : (y) )

#define max( x, y ) ( (x) > (y) ? (x) : (y) )

// Datatype
#define CONSTANT          101
#define INTEGER           102
#define FLOAT             103
#define DOUBLE            104
#define COMPLEX           105
#define DOUBLE_COMPLEX    106
#define INVALID_TYPE     -106



/* vector functions*/
void create_vector(integer datatype, void **A, integer M);
void create_realtype_vector(integer datatype, void **A, integer M);
void free_vector(void *A);
void reset_vector(integer datatype, void *A, integer M, integer incA);
void rand_vector(integer datatype, void *A, integer M, integer LDA);
void copy_vector(integer datatype, integer M, void *A, integer LDA, void *B, integer LDB);
void copy_realtype_vector(integer datatype, integer M, void *A, integer LDA, void *B, integer LDB);

/* matrix functions*/
void create_matrix(integer datatype, void **A, integer M, integer N);
void create_realtype_matrix(integer datatype, void **A, integer M, integer N);
integer get_datatype(char stype);
integer get_realtype(integer datatype);
void create_block_diagonal_matrix(integer datatype,void* wr, void* wi, void* lambda, integer m, integer n, integer lda);
void* get_m_ptr(integer datatype, void *A, integer M, integer N, integer LDA);
void free_matrix(void *A);
void rand_matrix(integer datatype, void *A, integer M, integer N, integer LDA);
void rand_sym_matrix(integer datatype, void *A, integer M, integer N, integer LDA);
void rand_spd_matrix(integer datatype, char *uplo, void **A, integer m,integer lda);
void rand_hermitian_matrix(integer datatype, integer n, void** A, integer lda);
void copy_matrix(integer datatype, char *uplo, integer M, integer N, void *A, integer LDA, void *B, integer LDB);
void copy_realtype_matrix(integer datatype, char *uplo, integer M, integer N, void *A, integer LDA, void *B, integer LDB);
void reset_matrix(integer datatype, integer M, integer N, void *A, integer LDA);
void set_identity_matrix(integer datatype, integer M, integer N, void *A, integer LDA);
void copy_submatrix(integer datatype, void* A, integer m_A, integer n_A, void *B, integer m_B, integer n_B, integer srow, integer scol);
/* orthgonality property of matrix */
double check_orthogonality(integer datatype, void *A, integer m, integer n, integer lda);
void copy_submatrix(integer datatype, void *A, integer m_A, integer n_A, void *B, integer m_B, integer n_B, integer srow, integer scol);
void rand_sym_tridiag_matrix(integer datatype, void *A, integer M, integer N, integer LDA);
void get_diagonal(integer datatype, void *A, integer m, integer n, integer lda, void *Diag);
void get_subdiagonal(integer datatype, void *A, integer m, integer n, integer lda, void *Subdiag);
void copy_sym_tridiag_matrix(integer datatype, void *D, void *E, integer M, integer N, void *B, integer LDA);

/* Division of complex types */
void c_div_t(scomplex *cp, scomplex *ap, scomplex *bp);
void z_div_t(dcomplex *cp, dcomplex *ap, dcomplex *bp);

/* work value calculation */
integer get_work_value( integer datatype, void *work );

/* Diagonal Scaling*/
void diagmv( integer datatype, integer m, integer n, void* x, integer incx, void* a, integer a_rs, integer a_cs );
void scalv( integer datatype, integer n, void* x, integer incx, void* y, integer incy );

/* set Transpose based on uplo */
void set_transpose(integer datatype, char *uplo, char *trans_A, char *trans_B);

/* Create diagonal matrix by copying elements from vector to matrix */
void diagonalize_vector(integer datatype, void* s, void* sigma, integer m, integer n, integer LDA);

/* To calculate matrix multiplication with real and complex datatypes */
void scgemv(char TRANS, integer real_alpha, integer m, integer n, scomplex* alpha, float* a, integer lda, scomplex* v, integer incv, float beta, scomplex* c, integer inc);

/* To find the maximum from the array */
void get_max(integer datatype, void *arr, void *max_val, integer n);
/* To find the minimum from the array */
void get_min(integer datatype, void *arr, void *min_val, integer n);
/* Reading matrix input data from a file */
void init_matrix_from_file(integer datatype, void* A,integer m, integer n, integer lda, FILE* fptr);
/* Reading vector input data from a file */
void init_vector_from_file(integer datatype, void* A, integer m, integer inc, FILE* fptr);

#endif // TEST_COMMON_H
