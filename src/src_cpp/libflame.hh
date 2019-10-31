/******************************************************************************
* Copyright (c) 2019 - present Advanced Micro Devices, Inc. All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*******************************************************************************/

/*! @file libflame.hh
 *  libflame.hh defines all the overloaded CPP functions to be invoked from
 *  template interfaces
 *  */
#ifndef LIBFLAME_HH
#define LIBFLAME_HH

#include  <complex>
extern "C" {
  #include "../../src/base/flamec/blis/include/blis_type_defs.h"
  #include "../../src/base/flamec/include/FLA_type_defs.h"
  #include "../../src/base/flamec/control/FLA_Cntl.h"
  #include "../../src/base/flamec/include/FLA_main_prototypes.h"
  #include "../../src/base/flamec/include/FLA_macro_defs.h"
  #include "../../src/base/flamec/include/FLA_macro_ptr_defs.h"
  #include "../../src/base/flamec/include/FLA_lapack_prototypes.h"
  #include "lapacke.h"
}

namespace libflame{

  // --- Cholesky factorization ---
  inline int potrf( int matrix_layout, char* uplo, int *n, float* a, int* lda )
  {
    return LAPACKE_spotrf( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int potrf( int matrix_layout, char* uplo, int *n, double* a, int* lda )
  {
    return LAPACKE_dpotrf( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int potrf( int matrix_layout, char* uplo, int *n, lapack_complex_float* a, int* lda )
  {
    return LAPACKE_cpotrf( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int potrf( int matrix_layout, char* uplo, int *n, lapack_complex_double* a, int* lda )
  {
    return LAPACKE_zpotrf( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int potf2( int matrix_layout, char* uplo, int *n, float* a, int* lda )
  {
    return LAPACKE_spotrf( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int potf2( int matrix_layout, char* uplo, int *n, double* a, int* lda )
  {
    return LAPACKE_dpotrf( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int potf2( int matrix_layout, char* uplo, int *n, lapack_complex_float* a, int* lda )
  {
    return LAPACKE_cpotrf( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int potf2( int matrix_layout, char* uplo, int *n, lapack_complex_double* a, int* lda )
  {
    return LAPACKE_zpotrf( matrix_layout, *uplo, *n, a, *lda );
  }

  // --- LU factorization with partial pivoting
  inline int getrf( int matrix_layout, int* m, int *n, float* a, int* lda, int* ipiv )
  {
    return LAPACKE_sgetrf( matrix_layout, *m, *n, a, *lda, ipiv );
  }
  inline int getrf( int matrix_layout, int* m, int *n, double*   a, int* lda, int* ipiv )
  {
    return LAPACKE_dgetrf( matrix_layout, *m, *n, a, *lda, ipiv );
  }
  inline int getrf( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, int* ipiv )
  {
    return LAPACKE_cgetrf( matrix_layout, *m, *n, a, *lda, ipiv );
  }
  inline int getrf( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, int* ipiv )
  {
    return LAPACKE_zgetrf( matrix_layout, *m, *n, a, *lda, ipiv );
  }

  inline int getf2( int matrix_layout, int* m, int *n, float* a, int* lda, int* ipiv )
  {
    return LAPACKE_sgetf2( matrix_layout, *m, *n, a, *lda, ipiv );
  }
  inline int getf2( int matrix_layout, int* m, int *n, double* a, int* lda, int* ipiv )
  {
    return LAPACKE_dgetf2( matrix_layout, *m, *n, a, *lda, ipiv );
  }
  inline int getf2( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, int* ipiv )
  {
    return LAPACKE_cgetf2( matrix_layout, *m, *n, a, *lda, ipiv );
  }
  inline int getf2( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, int* ipiv )
  {
    return LAPACKE_zgetf2( matrix_layout, *m, *n, a, *lda, ipiv );
  }
  // --- QR factorization (classic) ---
  inline int geqrf( int matrix_layout, int* m, int *n, float* a, int* lda, float* tau )
  {
    return LAPACKE_sgeqrf( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int geqrf( int matrix_layout, int* m, int *n, double* a, int* lda, double* tau )
  {
    return LAPACKE_dgeqrf( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int geqrf( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cgeqrf( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int geqrf( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zgeqrf( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int geqr2( int matrix_layout, int* m, int *n, float* a, int* lda, float* tau )
  {
    return LAPACKE_sgeqr2( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int geqr2( int matrix_layout, int* m, int *n, double* a, int* lda, double* tau )
  {
    return LAPACKE_dgeqr2( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int geqr2( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cgeqr2( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int geqr2( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zgeqr2( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int geqpf( int matrix_layout, int* m, int *n, float* a, int* lda, int* jpvt, float* tau )
  {
    printf(" Function sgeqpf() is been deprecated and has been replaced by function sgeqp3() \n");
    return LAPACKE_sgeqp3( matrix_layout, *m, *n, a, *lda, jpvt, tau );
  }
  inline int geqpf( int matrix_layout, int* m, int *n, double* a, int* lda, int* jpvt, double* tau )
  {
    printf(" Function dgeqpf() is been deprecated and has been replaced by function dgeqp3() \n");
    return LAPACKE_dgeqp3( matrix_layout, *m, *n, a, *lda, jpvt, tau );
  }
  inline int geqpf( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, int* jpvt, lapack_complex_float* tau )
  {
    printf(" Function cgeqpf() is been deprecated and has been replaced by function cgeqp3() \n");
    return LAPACKE_cgeqp3( matrix_layout, *m, *n, a, *lda, jpvt, tau);
  }
  inline int geqpf( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, int* jpvt, lapack_complex_double* tau)
  {
    printf(" Function zgeqpf() is been deprecated and has been replaced by function zgeqp3() \n");
    return LAPACKE_zgeqp3( matrix_layout, *m, *n, a, *lda, jpvt, tau);
  }

  inline int geqp3( int matrix_layout, int* m, int *n, float* a, int* lda, int* jpvt, float* tau )
  {
    return LAPACKE_sgeqp3( matrix_layout, *m, *n, a, *lda, jpvt, tau );
  }
  inline int geqp3( int matrix_layout, int* m, int *n, double* a, int* lda, int* jpvt, double* tau )
  {
    return LAPACKE_dgeqp3( matrix_layout, *m, *n, a, *lda, jpvt, tau );
  }
  inline int geqp3( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, int* jpvt, lapack_complex_float* tau )
  {
    return LAPACKE_cgeqp3( matrix_layout, *m, *n, a, *lda, jpvt, tau );
  }
  inline int geqp3( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, int* jpvt, lapack_complex_double* tau)
  {
    return LAPACKE_zgeqp3( matrix_layout, *m, *n, a, *lda, jpvt, tau );
  }

  // --- LQ factorization (classic) ---
  inline int gelqf( int matrix_layout, int* m, int *n, float* a, int* lda, float* tau )
  {
    return LAPACKE_sgelqf( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int gelqf( int matrix_layout, int* m, int *n, double* a, int* lda, double* tau )
  {
    return LAPACKE_dgelqf( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int gelqf( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cgelqf( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int gelqf( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zgelqf( matrix_layout, *m, *n, a, *lda, tau );
  }

  inline int gelq2( int matrix_layout, int* m, int *n, float* a, int* lda, float* tau )
  {
    return LAPACKE_sgelq2( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int gelq2( int matrix_layout, int* m, int *n, double* a, int* lda, double* tau )
  {
    return LAPACKE_dgelq2( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int gelq2( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cgelq2( matrix_layout, *m, *n, a, *lda, tau );
  }
  inline int gelq2( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zgelq2( matrix_layout, *m, *n, a, *lda, tau );
  }

  // --- LS solver ---
  inline int gelsd( int matrix_layout, int* m, int *n, int* nrhs, float* a, int* lda, float* b, int* ldb, float*  s, float*  rcond, int* rank )
  {
    return LAPACKE_sgelsd( matrix_layout, *m, *n, *nrhs, a, *lda, b, *ldb, s, *rcond, rank );
  }
  inline int gelsd( int matrix_layout, int* m, int *n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* s, double* rcond, int* rank )
  {
    return LAPACKE_dgelsd( matrix_layout, *m, *n, *nrhs, a, *lda, b, *ldb, s, *rcond, rank );
  }
  inline int gelsd( int matrix_layout, int* m, int *n, int* nrhs, lapack_complex_float* a, int* lda, lapack_complex_float* b, int* ldb, float*  s, float*  rcond, int* rank )
  {
    return LAPACKE_cgelsd( matrix_layout, *m, *n, *nrhs, a, *lda, b, *ldb, s, *rcond, rank );
  }
  inline int gelsd( int matrix_layout, int* m, int *n, int* nrhs, lapack_complex_double* a, int* lda, lapack_complex_double* b, int* ldb, double* s, double* rcond, int* rank )
  {
    return LAPACKE_zgelsd( matrix_layout, *m, *n, *nrhs, a, *lda, b, *ldb, s, *rcond, rank );
  }

  inline int gelss( int matrix_layout, int* m, int *n, int* nrhs, float* a, int* lda, float* b, int* ldb, float*  s, float*  rcond, int* rank )
  {
    return LAPACKE_sgelss( matrix_layout, *m, *n, *nrhs, a, *lda, b, *ldb, s, *rcond, rank );
  }
  inline int gelss( int matrix_layout, int* m, int *n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* s, double* rcond, int* rank )
  {
    return LAPACKE_dgelss( matrix_layout, *m, *n, *nrhs, a, *lda, b, *ldb, s, *rcond, rank );
  }
  inline int gelss( int matrix_layout, int* m, int *n, int* nrhs, lapack_complex_float* a, int* lda, lapack_complex_float* b, int* ldb, float*  s, float*  rcond, int* rank )
  {
    return LAPACKE_cgelss( matrix_layout, *m, *n, *nrhs, a, *lda, b, *ldb, s, *rcond, rank );
  }
  inline int gelss( int matrix_layout, int* m, int *n, int* nrhs, lapack_complex_double* a, int* lda, lapack_complex_double* b, int* ldb, double* s, double* rcond, int* rank)
  {
    return LAPACKE_zgelss( matrix_layout, *m, *n, *nrhs, a, *lda, b, *ldb, s, *rcond, rank );
  }

  // --- Triangular-transpose matrix multiply ---
  inline int lauum( int matrix_layout, char* uplo, int *n, float* a, int* lda )
  {
    return LAPACKE_slauum( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int lauum( int matrix_layout, char* uplo, int *n, double* a, int* lda )
  {
    return LAPACKE_dlauum( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int lauum( int matrix_layout, char* uplo, int *n, lapack_complex_float* a, int* lda )
  {
    return LAPACKE_clauum( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int lauum( int matrix_layout, char* uplo, int *n, lapack_complex_double* a, int* lda )
  {
    return LAPACKE_zlauum( matrix_layout, *uplo, *n, a, *lda );
  }

  inline int lauu2( int matrix_layout, char* uplo, int *n, float* a, int* lda )
  {
    return LAPACKE_slauum( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int lauu2( int matrix_layout, char* uplo, int *n, double* a, int* lda )
  {
    return LAPACKE_dlauum( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int lauu2( int matrix_layout, char* uplo, int *n, lapack_complex_float* a, int* lda )
  {
    return LAPACKE_clauum( matrix_layout, *uplo, *n, a, *lda );
  }
  inline int lauu2( int matrix_layout, char* uplo, int *n, lapack_complex_double* a, int* lda )
  {
    return LAPACKE_zlauum( matrix_layout, *uplo, *n, a, *lda );
  }

  // --- Symmetric (hermitian) positive definite matrix inversion ---
  inline int potri( int matrix_layout, char* uplo, int*  n, float* buff_A, int*  ldim_A )
  {
    return LAPACKE_spotri( matrix_layout, *uplo, *n, buff_A, *ldim_A );
  }
  inline int potri( int matrix_layout, char* uplo, int*  n, double* buff_A, int*  ldim_A )
  {
    return LAPACKE_dpotri( matrix_layout, *uplo, *n, buff_A, *ldim_A );
  }
  inline int potri( int matrix_layout, char* uplo, int*  n, lapack_complex_float* buff_A, int*  ldim_A )
  {
    return LAPACKE_cpotri( matrix_layout, *uplo, *n, buff_A, *ldim_A );
  }
  inline int potri( int matrix_layout, char* uplo, int*  n, lapack_complex_double* buff_A, int*  ldim_A )
  {
    return LAPACKE_zpotri( matrix_layout, *uplo, *n, buff_A, *ldim_A );
  }

  // --- Triangular matrix inversion ---
  inline int trtri( int matrix_layout, char* uplo, char* diag, int *n, float* a, int* lda )
  {
    return LAPACKE_strtri( matrix_layout, *uplo, *diag, *n, a, *lda );
  }
  inline int trtri( int matrix_layout, char* uplo, char* diag, int *n, double* a, int* lda )
  {
    return LAPACKE_dtrtri( matrix_layout, *uplo, *diag, *n, a, *lda );
  }
  inline int trtri( int matrix_layout, char* uplo, char* diag, int *n, lapack_complex_float* a, int* lda )
  {
    return LAPACKE_ctrtri( matrix_layout, *uplo, *diag, *n, a, *lda );
  }
  inline int trtri( int matrix_layout, char* uplo, char* diag, int *n, lapack_complex_double* a, int* lda )
  {
    return LAPACKE_ztrtri( matrix_layout, *uplo, *diag, *n, a, *lda );
  }

  inline int trti2( int matrix_layout, char* uplo, char* diag, int *n, float* a, int* lda )
  {
    return LAPACKE_strtri( matrix_layout, *uplo, *diag, *n, a, *lda );
  }
  inline int trti2( int matrix_layout, char* uplo, char* diag, int *n, double* a, int* lda )
  {
    return LAPACKE_dtrtri( matrix_layout, *uplo, *diag, *n, a, *lda );
  }
  inline int trti2( int matrix_layout, char* uplo, char* diag, int *n, lapack_complex_float* a, int* lda )
  {
    return LAPACKE_ctrtri( matrix_layout, *uplo, *diag, *n, a, *lda );
  }
  inline int trti2( int matrix_layout, char* uplo, char* diag, int *n, lapack_complex_double* a, int* lda )
  {
    return LAPACKE_ztrtri( matrix_layout, *uplo, *diag, *n, a, *lda );
  }

  // --- Triangular Sylvester equation solve ---
  inline int trsyl( int matrix_layout, char* transa, char* transb, int* isgn, int* m, int *n, float* a, int* lda, float* b, int* ldb, float* c, int* ldc, float* scale )
  {
    return LAPACKE_strsyl( matrix_layout, *transa, *transb, *isgn, *m, *n, a, *lda, b, *ldb, c, *ldc, scale );
  }
  inline int trsyl( int matrix_layout, char* transa, char* transb, int* isgn, int* m, int *n, double* a, int* lda, double* b, int* ldb, double* c, int* ldc, double* scale )
  {
    return LAPACKE_dtrsyl( matrix_layout, *transa, *transb, *isgn, *m, *n, a, *lda, b, *ldb, c, *ldc, scale );
  }
  inline int trsyl( int matrix_layout, char* transa, char* transb, int* isgn, int* m, int *n, lapack_complex_float* a, int* lda, lapack_complex_float* b, int* ldb, lapack_complex_float* c, int* ldc, float* scale )
  {
    return LAPACKE_ctrsyl( matrix_layout, *transa, *transb, *isgn, *m, *n, a, *lda, b, *ldb, c, *ldc, scale );
  }
  inline int trsyl( int matrix_layout, char* transa, char* transb, int* isgn, int* m, int *n, lapack_complex_double* a, int* lda, lapack_complex_double* b, int* ldb, lapack_complex_double* c, int* ldc, double* scale )
  {
    return LAPACKE_ztrsyl( matrix_layout, *transa, *transb, *isgn, *m, *n, a, *lda, b, *ldb, c, *ldc, scale );
  }

  // --- Reduction to upper Hessenberg form ---
  inline int gehrd( int matrix_layout, int *n, int* ilo, int* ihi, float* a, int* lda, float* tau )
  {
    return LAPACKE_sgehrd( matrix_layout, *n, *ilo, *ihi, a, *lda, tau );
  }
  inline int gehrd( int matrix_layout, int *n, int* ilo, int* ihi, double* a, int* lda, double* tau )
  {
    return LAPACKE_dgehrd( matrix_layout, *n, *ilo, *ihi, a, *lda, tau );
  }
  inline int gehrd( int matrix_layout, int *n, int* ilo, int* ihi, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cgehrd( matrix_layout, *n, *ilo, *ihi, a, *lda, tau );
  }
  inline int gehrd( int matrix_layout, int *n, int* ilo, int* ihi, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zgehrd( matrix_layout, *n, *ilo, *ihi, a, *lda, tau );
  }

  inline int gehd2( int matrix_layout, int *n, int* ilo, int* ihi, float* a, int* lda, float* tau )
  {
    return LAPACKE_sgehrd( matrix_layout, *n, *ilo, *ihi, a, *lda, tau );
  }
  inline int gehd2( int matrix_layout, int *n, int* ilo, int* ihi, double* a, int* lda, double* tau )
  {
    return LAPACKE_dgehrd( matrix_layout, *n, *ilo, *ihi, a, *lda, tau );
  }
  inline int gehd2( int matrix_layout, int *n, int* ilo, int* ihi, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cgehrd( matrix_layout, *n, *ilo, *ihi, a, *lda, tau );
  }
  inline int gehd2( int matrix_layout, int *n, int* ilo, int* ihi, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zgehrd( matrix_layout, *n, *ilo, *ihi, a, *lda, tau );
  }

  // --- Reduction to tridiagonal form ---
  inline int sytrd( int matrix_layout, char* uplo, int *n, float* a, int* lda, float*  d, float*  e, float* tau )
  {
    return LAPACKE_ssytrd( matrix_layout, *uplo, *n, a, *lda, d, e, tau );
  }
  inline int sytrd( int matrix_layout, char* uplo, int *n, double* a, int* lda, double* d, double* e, double* tau )
  {
    return LAPACKE_dsytrd( matrix_layout, *uplo, *n, a, *lda, d, e, tau );
  }
  inline int hetrd( int matrix_layout, char* uplo, int *n, lapack_complex_float* a, int* lda, float*  d, float*  e, lapack_complex_float* tau )
  {
    return LAPACKE_chetrd( matrix_layout, *uplo, *n, a, *lda, d, e, tau );
  }
  inline int hetrd( int matrix_layout, char* uplo, int *n, lapack_complex_double* a, int* lda, double* d, double* e, lapack_complex_double* tau )
  {
    return LAPACKE_zhetrd( matrix_layout, *uplo, *n, a, *lda, d, e, tau );
  }

  inline int sytd2( int matrix_layout, char* uplo, int *n, float* a, int* lda, float*  d, float*  e, float* tau )
  {
    return LAPACKE_ssytrd( matrix_layout, *uplo, *n, a, *lda, d, e, tau );
  }
  inline int sytd2( int matrix_layout, char* uplo, int *n, double* a, int* lda, double* d, double* e, double* tau )
  {
    return LAPACKE_dsytrd( matrix_layout, *uplo, *n, a, *lda, d, e, tau );
  }
  inline int hetd2( int matrix_layout, char* uplo, int *n, lapack_complex_float* a, int* lda, float*  d, float*  e, lapack_complex_float* tau )
  {
    return LAPACKE_chetrd( matrix_layout, *uplo, *n, a, *lda, d, e, tau );
  }
  inline int hetd2( int matrix_layout, char* uplo, int *n, lapack_complex_double* a, int* lda, double* d, double* e, lapack_complex_double* tau )
  {
    return LAPACKE_zhetrd( matrix_layout, *uplo, *n, a, *lda, d, e, tau );
  }

  // --- Reduction to bidiagonal form ---
  inline int gebrd( int matrix_layout, int* m, int *n, float* a, int* lda, float*  d, float*  e, float* tauq, float* taup )
  {
    return LAPACKE_sgebrd( matrix_layout, *m, *n, a, *lda, d, e, tauq, taup );
  }
  inline int gebrd( int matrix_layout, int* m, int *n, double* a, int* lda, double* d, double* e, double* tauq, double* taup )
  {
    return LAPACKE_dgebrd( matrix_layout, *m, *n, a, *lda, d, e, tauq, taup );
  }
  inline int gebrd( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, float*  d, float*  e, lapack_complex_float* tauq, lapack_complex_float* taup )
  {
    return LAPACKE_cgebrd( matrix_layout, *m, *n, a, *lda, d, e, tauq, taup );
  }
  inline int gebrd( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, double* d, double* e, lapack_complex_double* tauq, lapack_complex_double* taup )
  {
    return LAPACKE_zgebrd( matrix_layout, *m, *n, a, *lda, d, e, tauq, taup );
  }

  inline int gebd2( int matrix_layout, int* m, int *n, float* a, int* lda, float*  d, float*  e, float* tauq, float* taup )
  {
    return LAPACKE_sgebrd( matrix_layout, *m, *n, a, *lda, d, e, tauq, taup );
  }
  inline int gebd2( int matrix_layout, int* m, int *n, double* a, int* lda, double* d, double* e, double* tauq, double* taup )
  {
    return LAPACKE_dgebrd( matrix_layout, *m, *n, a, *lda, d, e, tauq, taup );
  }
  inline int gebd2( int matrix_layout, int* m, int *n, lapack_complex_float* a, int* lda, float*  d, float*  e, lapack_complex_float* tauq, lapack_complex_float* taup )
  {
    return LAPACKE_cgebrd( matrix_layout, *m, *n, a, *lda, d, e, tauq, taup );
  }
  inline int gebd2( int matrix_layout, int* m, int *n, lapack_complex_double* a, int* lda, double* d, double* e, lapack_complex_double* tauq, lapack_complex_double* taup )
  {
    return LAPACKE_zgebrd( matrix_layout, *m, *n, a, *lda, d, e, tauq, taup );
  }

  // --- Reduce Hermitian-definite generalized eigenproblem to standard form ---
  inline int sygst( int matrix_layout, int* itype, char* uplo, int *n, float* a, int* lda, float* b, int* ldb )
  {
    return LAPACKE_ssygst( matrix_layout, *itype, *uplo, *n, a, *lda, b, *ldb );
  }
  inline int sygst( int matrix_layout, int* itype, char* uplo, int *n, double* a, int* lda, double* b, int* ldb )
  {
    return LAPACKE_dsygst( matrix_layout, *itype, *uplo, *n, a, *lda, b, *ldb );
  }
  inline int hegst( int matrix_layout, int* itype, char* uplo, int *n, lapack_complex_float* a, int* lda, lapack_complex_float* b, int* ldb )
  {
    return LAPACKE_chegst( matrix_layout, *itype, *uplo, *n, a, *lda, b, *ldb );
  }
  inline int hegst( int matrix_layout, int* itype, char* uplo, int *n, lapack_complex_double* a, int* lda, lapack_complex_double* b, int* ldb )
  {
    return LAPACKE_zhegst( matrix_layout, *itype, *uplo, *n, a, *lda, b, *ldb );
  }

  inline int sygs2( int matrix_layout, int* itype, char* uplo, int *n, float* a, int* lda, float* b, int* ldb )
  {
    return LAPACKE_ssygst( matrix_layout, *itype, *uplo, *n, a, *lda, b, *ldb );
  }
  inline int sygs2( int matrix_layout, int* itype, char* uplo, int *n, double* a, int* lda, double* b, int* ldb )
  {
    return LAPACKE_dsygst( matrix_layout, *itype, *uplo, *n, a, *lda, b, *ldb );
  }

  inline int hegs2( int matrix_layout, int* itype, char* uplo, int *n, lapack_complex_float* a, int* lda, lapack_complex_float* b, int* ldb )
  {
    return LAPACKE_chegst( matrix_layout, *itype, *uplo, *n, a, *lda, b, *ldb );
  }
  inline int hegs2( int matrix_layout, int* itype, char* uplo, int *n, lapack_complex_double* a, int* lda, lapack_complex_double* b, int* ldb )
  {
    return LAPACKE_zhegst( matrix_layout, *itype, *uplo, *n, a, *lda, b, *ldb );
  }

  // --- Accumulate block Householder matrix T (classic) ---
  inline int larft( int matrix_layout, char* direct, char* storev, int *n, int* k, float* v, int* ldv, float* tau, float* t, int* ldt )
  {
    return LAPACKE_slarft( matrix_layout, *direct, *storev, *n, *k, v, *ldv, tau,  t, *ldt );
  }
  inline int larft( int matrix_layout, char* direct, char* storev, int *n, int* k, double* v, int* ldv, double* tau, double* t, int* ldt )
  {
    return LAPACKE_dlarft( matrix_layout, *direct, *storev, *n, *k, v, *ldv, tau,  t, *ldt );
  }
  inline int larft( int matrix_layout, char* direct, char* storev, int *n, int* k, lapack_complex_float* v, int* ldv, lapack_complex_float* tau, lapack_complex_float* t, int* ldt )
  {
    return LAPACKE_clarft( matrix_layout, *direct, *storev, *n, *k, v, *ldv, tau,  t, *ldt );
  }
  inline int larft( int matrix_layout, char* direct, char* storev, int *n, int* k, lapack_complex_double* v, int* ldv, lapack_complex_double* tau, lapack_complex_double* t, int* ldt )
  {
    return LAPACKE_zlarft( matrix_layout, *direct, *storev, *n, *k, v, *ldv, tau,  t, *ldt );
  }

  // --- Generate a Householder vector (classic) ---
  inline int larfg( int *n, float* alpha, float* x, int* incx, float* tau )
  {
    return LAPACKE_slarfg( *n, alpha, x, *incx, tau );
  }
  inline int larfg( int *n, double* alpha, double* x, int* incx, double* tau )
  {
    return LAPACKE_dlarfg( *n, alpha, x, *incx, tau );
  }
  inline int larfg( int *n, lapack_complex_float* alpha, lapack_complex_float* x, int* incx, lapack_complex_float* tau )
  {
    return LAPACKE_clarfg( *n, alpha, x, *incx, tau );
  }
  inline int larfg( int *n, lapack_complex_double* alpha, lapack_complex_double* x, int* incx, lapack_complex_double* tau )
  {
    return LAPACKE_zlarfg( *n, alpha, x, *incx, tau );
  }

  inline int larfgp( int *n, float* alpha, float* x, int* incx, float* tau )
  {
    return LAPACKE_slarfg( *n, alpha, x, *incx, tau );
  }
  inline int larfgp( int *n, double* alpha, double* x, int* incx, double* tau )
  {
    return LAPACKE_dlarfg( *n, alpha, x, *incx, tau );
  }
  inline int larfgp( int *n, lapack_complex_float* alpha, lapack_complex_float* x, int* incx, lapack_complex_float* tau )
  {
    return LAPACKE_clarfg( *n, alpha, x, *incx, tau );
  }
  inline int larfgp( int *n, lapack_complex_double* alpha, lapack_complex_double* x, int* incx, lapack_complex_double* tau )
  {
    return LAPACKE_zlarfg( *n, alpha, x, *incx, tau );
  }

  // --- Form Q from QR factorization ---
  inline int orgqr( int matrix_layout, int* m, int *n, int* k, float* a, int* lda, float* tau )
  {
     return LAPACKE_sorgqr( matrix_layout, *m, *n, *k, a, *lda, tau );
  }
  inline int orgqr( int matrix_layout, int* m, int *n, int* k, double* a, int* lda, double* tau )
  {
    return LAPACKE_dorgqr( matrix_layout, *m, *n, *k, a, *lda, tau );
  }
  inline int ungqr( int matrix_layout, int* m, int *n, int* k, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cungqr( matrix_layout, *m, *n, *k, a, *lda, tau );
  }
  inline int ungqr( int matrix_layout, int* m, int *n, int* k, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zungqr( matrix_layout, *m, *n, *k, a, *lda, tau );
  }

  // --- Apply Q or Q' from QR factorization ---
  inline int ormqr( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, float* a, int* lda, float* tau, float* c, int* ldc )
  {
    return LAPACKE_sormqr( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int ormqr( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, double* a, int* lda, double* tau, double* c, int* ldc )
  {
    return LAPACKE_dormqr( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int unmqr( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, lapack_complex_float* a, int* lda, lapack_complex_float* tau, lapack_complex_float* c, int* ldc )
  {
    return LAPACKE_cunmqr( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int unmqr( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, lapack_complex_double* a, int* lda, lapack_complex_double* tau, lapack_complex_double* c, int* ldc )
  {
    return LAPACKE_zunmqr( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }

  inline int orm2r( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, float* a, int* lda, float* tau, float* c, int* ldc )
  {
    return LAPACKE_sormqr( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int orm2r( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, double* a, int* lda, double* tau, double* c, int* ldc )
  {
    return LAPACKE_dormqr( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int unm2r( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, lapack_complex_float* a, int* lda, lapack_complex_float* tau, lapack_complex_float* c, int* ldc )
  {
    return LAPACKE_cunmqr( matrix_layout, *side, *trans, *m, *n, *k,  a, *lda, tau, c, *ldc );
  }
  inline int unm2r( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, lapack_complex_double* a, int* lda, lapack_complex_double* tau, lapack_complex_double* c, int* ldc )
  {
    return LAPACKE_zunmqr( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }

  // --- Form Q from LQ factorization ---
  inline int orglq( int matrix_layout, int* m, int *n, int* k, float* a, int* lda, float* tau )
  {
    return LAPACKE_sorglq( matrix_layout, *m, *n, *k, a, *lda, tau );
  }
  inline int orglq( int matrix_layout, int* m, int *n, int* k, double* a, int* lda, double* tau )
  {
    return LAPACKE_dorglq( matrix_layout, *m, *n, *k, a, *lda, tau );
  }
  inline int unglq( int matrix_layout, int* m, int *n, int* k, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cunglq( matrix_layout, *m, *n, *k, a, *lda, tau );
  }
  inline int unglq( int matrix_layout, int* m, int *n, int* k, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zunglq( matrix_layout, *m, *n, *k, a, *lda, tau );
  }

  // --- Apply Q or Q' from LQ factorization ---
  inline int ormlq( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, float* a, int* lda, float* tau, float* c, int* ldc )
  {
    return LAPACKE_sormlq( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int ormlq( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, double* a, int* lda, double* tau, double* c, int* ldc )
  {
    return LAPACKE_dormlq( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int unmlq( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, lapack_complex_float* a, int* lda, lapack_complex_float* tau, lapack_complex_float* c, int* ldc )
  {
    return LAPACKE_cunmlq( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int unmlq( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, lapack_complex_double* a, int* lda, lapack_complex_double* tau, lapack_complex_double* c, int* ldc )
  {
    return LAPACKE_zunmlq( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }

  inline int orml2( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, float* a, int* lda, float* tau, float* c, int* ldc )
  {
    return LAPACKE_sormlq( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int orml2( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, double* a, int* lda, double* tau, double* c, int* ldc )
  {
    return LAPACKE_dormlq( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int unml2( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, lapack_complex_float* a, int* lda, lapack_complex_float* tau, lapack_complex_float* c, int* ldc )
  {
    return LAPACKE_cunmlq( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int unml2( int matrix_layout, char* side, char* trans, int* m, int *n, int* k, lapack_complex_double* a, int* lda, lapack_complex_double* tau, lapack_complex_double* c, int* ldc )
  {
    return LAPACKE_zunmlq( matrix_layout, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }

  // --- Form Q from tridiagonal reduction ---
  inline int orgtr( int matrix_layout, char* uplo, int* m, float* a, int* lda, float* tau )
  {
    return LAPACKE_sorgtr( matrix_layout, *uplo, *m, a, *lda, tau );
  }
  inline int orgtr( int matrix_layout, char* uplo, int* m, double* a, int* lda, double* tau )
  {
    return LAPACKE_dorgtr( matrix_layout, *uplo, *m, a, *lda, tau );
  }
  inline int ungtr( int matrix_layout, char* uplo, int* m, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cungtr( matrix_layout, *uplo, *m, a, *lda, tau );
  }
  inline int ungtr( int matrix_layout, char* uplo, int* m, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zungtr( matrix_layout, *uplo, *m, a, *lda, tau );
  }

  // --- Apply Q or Q' from tridiagonal reduction ---
  inline int ormtr( int matrix_layout, char* side, char* uplo, char* trans, int* m, int *n, float* a, int* lda, float* tau, float* c, int* ldc )
  {
    return LAPACKE_sormtr( matrix_layout, *side, *uplo, *trans, *m, *n, a, *lda, tau, c, *ldc );
  }
  inline int ormtr( int matrix_layout, char* side, char* uplo, char* trans, int* m, int *n, double* a, int* lda, double* tau, double* c, int* ldc )
  {
    return LAPACKE_dormtr( matrix_layout, *side, *uplo, *trans, *m, *n, a, *lda, tau, c, *ldc );
  }
  inline int unmtr( int matrix_layout, char* side, char* uplo, char* trans, int* m, int *n, lapack_complex_float* a, int* lda, lapack_complex_float* tau, lapack_complex_float* c, int* ldc )
  {
    return LAPACKE_cunmtr( matrix_layout, *side, *uplo, *trans, *m, *n, a, *lda, tau, c, *ldc );
  }
  inline int unmtr( int matrix_layout, char* side, char* uplo, char* trans, int* m, int *n, lapack_complex_double* a, int* lda, lapack_complex_double* tau, lapack_complex_double* c, int* ldc )
  {
    return LAPACKE_zunmtr( matrix_layout, *side, *uplo, *trans, *m, *n, a, *lda, tau, c, *ldc );
  }

  // --- Form Q from bidiagonal reduction ---
  inline int orgbr( int matrix_layout, char* vect, int* m, int *n, int* k, float* a, int* lda, float* tau )
  {
    return LAPACKE_sorgbr( matrix_layout, *vect, *m, *n, *k, a, *lda, tau );
  }
  inline int orgbr( int matrix_layout, char* vect, int* m, int *n, int* k, double* a, int* lda, double* tau )
  {
    return LAPACKE_dorgbr( matrix_layout, *vect, *m, *n, *k, a, *lda, tau );
  }
  inline int ungbr( int matrix_layout, char* vect, int* m, int *n, int* k, lapack_complex_float* a, int* lda, lapack_complex_float* tau )
  {
    return LAPACKE_cungbr( matrix_layout, *vect, *m, *n, *k, a, *lda, tau );
  }
  inline int ungbr( int matrix_layout, char* vect, int* m, int *n, int* k, lapack_complex_double* a, int* lda, lapack_complex_double* tau )
  {
    return LAPACKE_zungbr( matrix_layout, *vect, *m, *n, *k, a, *lda, tau );
  }

  // --- Apply Q or Q' from bidiagonal reduction ---
  inline int ormbr( int matrix_layout, char* vect, char* side, char* trans, int* m, int *n, int* k, float* a, int* lda, float* tau, float* c, int* ldc )
  {
    return LAPACKE_sormbr( matrix_layout, *vect, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int ormbr( int matrix_layout, char* vect, char* side, char* trans, int* m, int *n, int* k, double* a, int* lda, double* tau, double* c, int* ldc )
  {
    return LAPACKE_dormbr( matrix_layout, *vect, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int unmbr( int matrix_layout, char* vect, char* side, char* trans, int* m, int *n, int* k, lapack_complex_float* a, int* lda, lapack_complex_float* tau, lapack_complex_float* c, int* ldc )
  {
    return LAPACKE_cunmbr( matrix_layout, *vect, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }
  inline int unmbr( int matrix_layout, char* vect, char* side, char* trans, int* m, int *n, int* k, lapack_complex_double* a, int* lda, lapack_complex_double* tau, lapack_complex_double* c, int* ldc )
  {
    return LAPACKE_zunmbr( matrix_layout, *vect, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc );
  }

  // --- Tridiagonal QR algorithm ---
  inline int steqr( int matrix_layout, char* jobz, int *n, float* d, float* e, float* z, int* ldz )
  {
    return LAPACKE_ssteqr( matrix_layout, *jobz, *n, d, e, z, *ldz );
  }
  inline int steqr( int matrix_layout, char* jobz, int *n, double* d, double* e, double* z, int* ldz )
  {
    return LAPACKE_dsteqr( matrix_layout, *jobz, *n, d, e, z, *ldz );
  }
  inline int steqr( int matrix_layout, char* jobz, int *n, float* d, float* e, lapack_complex_float* z, int* ldz )
  {
    return LAPACKE_csteqr( matrix_layout, *jobz, *n, d, e, z, *ldz );
  }
  inline int steqr( int matrix_layout, char* jobz, int *n, double* d, double* e, lapack_complex_double* z, int* ldz )
  {
    return LAPACKE_zsteqr( matrix_layout, *jobz, *n, d, e, z, *ldz );
  }

  // --- Tridiagonal divide-and-conquer algorithm ---
  inline int stedc( int matrix_layout, char* compz, int *n, float* d, float* e, float* z, int* ldz )
  {
    return LAPACKE_sstedc( matrix_layout, *compz, *n, d, e, z, *ldz );
  }
  inline int stedc( int matrix_layout, char* compz, int *n, double* d, double* e, double* z, int* ldz )
  {
    return LAPACKE_dstedc( matrix_layout, *compz, *n, d, e, z, *ldz );
  }
  inline int stedc( int matrix_layout, char* compz, int *n, float* d, float* e, lapack_complex_float* z, int* ldz )
  {
    return LAPACKE_cstedc( matrix_layout, *compz, *n, d, e, z, *ldz );
  }
  inline int stedc( int matrix_layout, char* compz, int *n, double* d, double* e, lapack_complex_double* z, int* ldz )
  {
    return LAPACKE_zstedc( matrix_layout, *compz, *n, d, e, z, *ldz );
  }

  // --- Tridiagonal MRRR algorithm ---
  inline int stemr( int matrix_layout, char* jobz, char* range, int *n, float*  d, float*  e, int* vl, int* vu, int* il, int* iu, int* m, float*  w, float* z, int* ldz, int* nzc, int* isuppz, int* tryrac )
  {
    return LAPACKE_sstemr( matrix_layout, *jobz, *range, *n, d, e, *vl, *vu, *il, *iu, m, w, z, *ldz, *nzc, isuppz, tryrac );
  }
  inline int stemr( int matrix_layout, char* jobz, char* range, int *n, double* d, double* e, int* vl, int* vu, int* il, int* iu, int* m, double* w, double* z, int* ldz, int* nzc, int* isuppz, int* tryrac )
  {
    return LAPACKE_dstemr( matrix_layout, *jobz, *range, *n, d, e, *vl, *vu, *il, *iu, m, w, z, *ldz, *nzc, isuppz, tryrac );
  }
  inline int stemr( int matrix_layout, char* jobz, char* range, int *n, float*  d, float*  e, int* vl, int* vu, int* il, int* iu, int* m, float*  w, lapack_complex_float* z, int* ldz, int* nzc, int* isuppz, int* tryrac )
  {
    return LAPACKE_cstemr( matrix_layout, *jobz, *range, *n, d, e, *vl, *vu, *il, *iu, m, w, z, *ldz, *nzc, isuppz, tryrac );
  }
  inline int stemr( int matrix_layout, char* jobz, char* range, int *n, double* d, double* e, int* vl, int* vu, int* il, int* iu, int* m, double* w, lapack_complex_double* z, int* ldz, int* nzc, int* isuppz, int* tryrac )
  {
    return LAPACKE_zstemr( matrix_layout, *jobz, *range, *n, d, e, *vl, *vu, *il, *iu, m, w, z, *ldz, *nzc, isuppz, tryrac );
  }

  // --- Hermitian eigenvalue decomposition (QR algorithm) ---
  inline int syev( int matrix_layout, char* jobz, char* uplo, int *n, float* a, int* lda, float*  w )
  {
    return LAPACKE_ssyev( matrix_layout, *jobz, *uplo, *n, a, *lda, w );
  }
  inline int syev( int matrix_layout, char* jobz, char* uplo, int *n, double* a, int* lda, double* w)
  {
    return LAPACKE_dsyev( matrix_layout, *jobz, *uplo, *n, a, *lda, w );
  }
  inline int heev( int matrix_layout, char* jobz, char* uplo, int *n, lapack_complex_float* a, int* lda, float*  w )
  {
    return LAPACKE_cheev( matrix_layout, *jobz, *uplo, *n, a, *lda, w );
  }
  inline int heev( int matrix_layout, char* jobz, char* uplo, int *n, lapack_complex_double* a, int* lda, double* w)
  {
    return LAPACKE_zheev( matrix_layout, *jobz, *uplo, *n, a, *lda, w );
  }

  // --- Hermitian eigenvalue decomposition (divide-and-conquer) ---
  inline int syevd( int matrix_layout, char* jobz, char* uplo, int *n, float* a, int* lda, float*  w)
  {
    return LAPACKE_ssyevd( matrix_layout, *jobz, *uplo, *n, a, *lda, w );
  }
  inline int syevd( int matrix_layout, char* jobz, char* uplo, int *n, double* a, int* lda, double* w)
  {
    return LAPACKE_dsyevd( matrix_layout, *jobz, *uplo, *n, a, *lda, w );
  }
  inline int heevd( int matrix_layout, char* jobz, char* uplo, int *n, lapack_complex_float* a, int* lda, float*  w )
  {
    return LAPACKE_cheevd( matrix_layout, *jobz, *uplo, *n, a, *lda, w );
  }
  inline int heevd( int matrix_layout, char* jobz, char* uplo, int *n, lapack_complex_double* a, int* lda, double* w )
  {
    return LAPACKE_zheevd( matrix_layout, *jobz, *uplo, *n, a, *lda, w );
  }

  // --- Hermitian eigenvalue decomposition (MRRR) ---
  inline int syevr( int matrix_layout, char* jobz, char* range, char* uplo, int *n, float* a, int* lda, float*  vl, float*  vu, int* il, int* iu, float*  abstol, int* m, float*  w, float* z, int* ldz, int* isuppz )
  {
    return LAPACKE_ssyevr( matrix_layout, *jobz, *range, *uplo, *n, a, *lda, *vl, *vu, *il, *iu, *abstol, m, w, z, *ldz, isuppz );
  }
  inline int syevr( int matrix_layout, char* jobz, char* range, char* uplo, int *n, double* a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, double* z, int* ldz, int* isuppz )
  {
    return LAPACKE_dsyevr( matrix_layout, *jobz, *range, *uplo, *n, a, *lda, *vl, *vu, *il, *iu, *abstol, m, w, z, *ldz, isuppz );
  }
  inline int heevr( int matrix_layout, char* jobz, char* range, char* uplo, int *n, lapack_complex_float* a, int* lda, float*  vl, float*  vu, int* il, int* iu, float*  abstol, int* m, float*  w, lapack_complex_float* z, int* ldz, int* isuppz )
  {
    return LAPACKE_cheevr( matrix_layout, *jobz, *range, *uplo, *n, a, *lda, *vl, *vu, *il, *iu, *abstol, m, w, z, *ldz, isuppz );
  }
  inline int heevr( int matrix_layout, char* jobz, char* range, char* uplo, int *n, lapack_complex_double* a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, lapack_complex_double* z, int* ldz, int* isuppz )
  {
    return LAPACKE_zheevr( matrix_layout, *jobz, *range, *uplo, *n, a, *lda, *vl, *vu, *il, *iu, *abstol, m, w, z, *ldz, isuppz );
  }

  // --- Bidiagonal QR algorithm ---
  inline int bdsqr( int matrix_layout, char* uplo, int *n, int* ncvt, int* nru, int* ncc, float* d, float* e, float* vt, int* ldvt, float* u, int* ldu, float* c, int* ldc )
  {
    return LAPACKE_sbdsqr( matrix_layout, *uplo, *n, *ncvt, *nru, *ncc, d, e, vt, *ldvt, u, *ldu, c, *ldc );
  }
  inline int bdsqr( int matrix_layout, char* uplo, int *n, int* ncvt, int* nru, int* ncc, double* d, double* e, double* vt, int* ldvt, double* u, int* ldu, double* c, int* ldc)
  {
    return LAPACKE_dbdsqr( matrix_layout, *uplo, *n, *ncvt, *nru, *ncc, d, e, vt, *ldvt, u, *ldu, c, *ldc );
  }
  inline int bdsqr( int matrix_layout, char* uplo, int *n, int* ncvt, int* nru, int* ncc, float* d, float* e, lapack_complex_float* vt, int* ldvt, lapack_complex_float* u, int* ldu, lapack_complex_float* c, int* ldc )
  {
    return LAPACKE_cbdsqr( matrix_layout, *uplo, *n, *ncvt, *nru, *ncc, d, e, vt, *ldvt, u, *ldu, c, *ldc );
  }
  inline int bdsqr( int matrix_layout, char* uplo, int *n, int* ncvt, int* nru, int* ncc, double* d, double* e, lapack_complex_double* vt, int* ldvt, lapack_complex_double* u, int* ldu, lapack_complex_double* c, int* ldc)
  {
    return LAPACKE_zbdsqr( matrix_layout, *uplo, *n, *ncvt, *nru, *ncc, d, e, vt, *ldvt, u, *ldu, c, *ldc );
  }

  // --- Bidiagonal divide-and-conquor algorithm ---
  inline int bdsdc( int matrix_layout, char* uplo, char* compq, int *n, float*  d, float*  e, float*  u, int* ldu, float*  vt, int* ldvt, float*  q, int*  iq )
  {
    return LAPACKE_sbdsdc( matrix_layout, *uplo, *compq, *n, d, e, u, *ldu, vt, *ldvt, q, iq );
  }
  inline int bdsdc( int matrix_layout, char* uplo, char* compq, int *n, double* d, double* e, double* u, int* ldu, double* vt, int* ldvt, double* q, int* iq )
  {
    return LAPACKE_dbdsdc( matrix_layout, *uplo, *compq, *n, d, e, u, *ldu, vt, *ldvt, q, iq );
  }

  // --- General matrix singular value decomposition (QR algorithm) ---
  inline int gesvd( int matrix_layout, char* jobu, char* jobv, int* m, int *n, float* a, int* lda, float* s, float* u, int* ldu, float* vt, int* ldvt, float *superb )
  {
    return LAPACKE_sgesvd( matrix_layout, *jobu, *jobv, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt, superb );
  }
  inline int gesvd( int matrix_layout, char* jobu, char* jobv, int* m, int *n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double *superb )
  {
    return LAPACKE_dgesvd( matrix_layout, *jobu, *jobv, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt, superb );
  }
  inline int gesvd( int matrix_layout, char* jobu, char* jobv, int* m, int *n, lapack_complex_float* a, int* lda, float*  s, lapack_complex_float* u, int* ldu, lapack_complex_float* vt, int* ldvt, float *superb )
  {
    return LAPACKE_cgesvd( matrix_layout, *jobu, *jobv, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt, superb );
  }
  inline int gesvd( int matrix_layout, char* jobu, char* jobv, int* m, int *n, lapack_complex_double* a, int* lda, double* s, lapack_complex_double* u, int* ldu, lapack_complex_double* vt, int* ldvt, double *superb )
  {
    return LAPACKE_zgesvd( matrix_layout, *jobu, *jobv, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt, superb );
  }

  // --- General matrix singular value decomposition (divide-and-conquer) ---
  inline int gesdd( int matrix_layout, char* jobz, int* m, int *n, float* a, int* lda, float*  s, float* u, int* ldu, float* vt, int* ldvt )
  {
    return LAPACKE_sgesdd( matrix_layout, *jobz, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt );
  }
  inline int gesdd( int matrix_layout, char* jobz, int* m, int *n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt )
  {
    return LAPACKE_dgesdd( matrix_layout, *jobz, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt );
  }
  inline int gesdd( int matrix_layout, char* jobz, int* m, int *n, lapack_complex_float* a, int* lda, float*  s, lapack_complex_float* u, int* ldu, lapack_complex_float* vt, int* ldvt )
  {
    return LAPACKE_cgesdd( matrix_layout, *jobz, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt );
  }
  inline int gesdd( int matrix_layout, char* jobz, int* m, int *n, lapack_complex_double* a, int* lda, double* s, lapack_complex_double* u, int* ldu, lapack_complex_double* vt, int* ldvt )
  {
    return LAPACKE_zgesdd( matrix_layout, *jobz, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt );
  }

  // --- Swap rows ---
  inline int laswp( int matrix_layout, int *n, float* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    return LAPACKE_slaswp( matrix_layout, *n, a, *lda, *k1, *k2, ipiv, *incx );
  }
  inline int laswp( int matrix_layout, int *n, double* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    return LAPACKE_dlaswp( matrix_layout, *n, a, *lda, *k1, *k2, ipiv, *incx );
  }
  inline int laswp( int matrix_layout, int *n, lapack_complex_float* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    return LAPACKE_claswp( matrix_layout, *n, a, *lda, *k1, *k2, ipiv, *incx );
  }
  inline int laswp( int matrix_layout, int *n, lapack_complex_double* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    return LAPACKE_zlaswp( matrix_layout, *n, a, *lda, *k1, *k2, ipiv, *incx );
  }

  // --- Initialize a matrix ---
  inline int laset( int matrix_layout, char* uplo, int* m, int *n, float* alpha, float* beta, float* a, int* lda )
  {
    return LAPACKE_slaset( matrix_layout, *uplo, *m, *n, *alpha, *beta, a, *lda );
  }
  inline int laset( int matrix_layout, char* uplo, int* m, int *n, double* alpha, double* beta, double* a, int* lda )
  {
    return LAPACKE_dlaset( matrix_layout, *uplo, *m, *n, *alpha, *beta, a, *lda );
  }
  inline int laset( int matrix_layout, char* uplo, int* m, int *n, lapack_complex_float* alpha, lapack_complex_float* beta, lapack_complex_float* a, int* lda )
  {
    return LAPACKE_claset( matrix_layout, *uplo, *m, *n, *alpha, *beta, a, *lda );
  }
  inline int laset( int matrix_layout, char* uplo, int* m, int *n, lapack_complex_double* alpha, lapack_complex_double* beta, lapack_complex_double* a, int* lda )
  {
    return LAPACKE_zlaset( matrix_layout, *uplo, *m, *n, *alpha, *beta, a, *lda );
  }

}//namespace libflame

#endif  //  #ifndef LIBFLAME_HH

