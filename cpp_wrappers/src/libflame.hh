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
#include "FLAME.h"
#include "../../src/map/lapack2flamec/FLA_lapack2flame_util_defs.h"
}

namespace libflame{

  // --- Cholesky factorization ---
  inline int potrf( char* uplo, int* n, float* a, int* lda, int* info )
  {
    return F77_spotrf(uplo, n, a, lda, info );
  }
  inline int potrf( char* uplo, int* n, double* a, int* lda, int* info )
  {
    return F77_dpotrf(uplo, n, a, lda, info );
  }
  inline int potrf( char* uplo, int* n, scomplex* a, int* lda, int* info )
  {
    return F77_cpotrf(uplo, n, a, lda, info );
  }
  inline int potrf( char* uplo, int* n, dcomplex* a, int* lda, int* info )
  {
    return F77_zpotrf(uplo, n, a, lda, info );
  }

  inline int potf2( char* uplo, int* n, float* a, int* lda, int* info )
  {
    return F77_spotf2( uplo, n, a, lda, info );
  }
  inline int potf2( char* uplo, int* n, double*    a, int* lda, int* info )
  {
    return F77_dpotf2( uplo, n, a, lda, info );
  }
  inline int potf2( char* uplo, int* n, scomplex* a, int* lda, int* info )
  {
    return F77_cpotf2( uplo, n, a, lda, info );
  }
  inline int potf2( char* uplo, int* n, dcomplex* a, int* lda, int* info )
  {
    return F77_zpotf2( uplo, n, a, lda, info );
  }

  // --- LU factorization with partial pivoting
  inline int getrf( int* m, int* n, float*    a, int* lda, int* ipiv, int* info )
  {
    return F77_sgetrf( m, n, a, lda, ipiv, info );
  }
  inline int getrf( int* m, int* n, double*   a, int* lda, int* ipiv, int* info )
  {
    return F77_dgetrf( m, n, a, lda, ipiv, info );
  }
  inline int getrf( int* m, int* n, scomplex* a, int* lda, int* ipiv, int* info )
  {
    return F77_cgetrf( m, n, a, lda, ipiv, info );
  }
  inline int getrf( int* m, int* n, dcomplex* a, int* lda, int* ipiv, int* info )
  {
    return F77_zgetrf( m, n, a, lda, ipiv, info );
  }

  inline int getf2( int* m, int* n, float* a, int* lda, int* ipiv, int* info )
  {
    return F77_sgetf2( m, n, a, lda, ipiv, info );
  }
  inline int getf2( int* m, int* n, double* a, int* lda, int* ipiv, int* info )
  {
    return F77_dgetf2( m, n, a, lda, ipiv, info );
  }
  inline int getf2( int* m, int* n, scomplex* a, int* lda, int* ipiv, int* info )
  {
    return F77_cgetf2( m, n, a, lda, ipiv, info );
  }
  inline int getf2( int* m, int* n, dcomplex* a, int* lda, int* ipiv, int* info )
  {
    return F77_zgetf2( m, n, a, lda, ipiv, info );
  }

  // --- QR factorization (classic) ---
  inline int geqrf( int* m, int* n, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    return F77_sgeqrf( m, n, a, lda, tau, work, lwork, info );
  }
  inline int geqrf( int* m, int* n, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    return F77_dgeqrf( m, n, a, lda, tau, work, lwork, info );
  }
  inline int geqrf( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    return F77_cgeqrf( m, n, a, lda, tau, work, lwork, info );
  }
  inline int geqrf( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    return F77_zgeqrf( m, n, a, lda, tau, work, lwork, info );
  }

  inline int geqr2( int* m, int* n, float* a, int* lda, float* tau, float* work, int* info )
  {
    return F77_sgeqr2( m, n, a, lda, tau, work, info );
  }
  inline int geqr2( int* m, int* n, double*    a, int* lda, double*    tau, double*    work, int* info )
  {
    return F77_dgeqr2( m, n, a, lda, tau, work, info );
  }
  inline int geqr2( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* info )
  {
    return F77_cgeqr2( m, n, a, lda, tau, work, info );
  }
  inline int geqr2( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* info )
  {
    return F77_zgeqr2( m, n, a, lda, tau, work, info );
  }


  inline int geqpf( int* m, int* n, float* a, int* lda, int* jpvt, float* tau, float* work,     int* info )
  {
    return F77_sgeqpf( m, n, a, lda, jpvt, tau, work, info );
  }
  inline int geqpf( int* m, int* n, double*    a, int* lda, int* jpvt, double*    tau, double*    work,     int* info )
  {
    return F77_dgeqpf( m, n, a, lda, jpvt, tau, work, info );
  }
  inline int geqpf( int* m, int* n, scomplex* a, int* lda, int* jpvt, scomplex* tau, scomplex* work, float*  rwork, int* info )
  {
    return F77_cgeqpf( m, n, a, lda, jpvt, tau, work, rwork, info );
  }
  inline int geqpf( int* m, int* n, dcomplex* a, int* lda, int* jpvt, dcomplex* tau, dcomplex* work, double* rwork, int* info )
  {
    return F77_zgeqpf( m, n, a, lda, jpvt, tau, work, rwork, info );
  }

  inline int geqp3( int* m, int* n, float* a, int* lda, int* jpvt, float* tau, float* work, int* lwork,     int* info )
  {
    return F77_sgeqp3( m, n, a, lda, jpvt, tau, work, lwork, info );
  }
  inline int geqp3( int* m, int* n, double*    a, int* lda, int* jpvt, double*    tau, double*    work, int* lwork,     int* info )
  {
    return F77_dgeqp3( m, n, a, lda, jpvt, tau, work, lwork, info );
  }
  inline int geqp3( int* m, int* n, scomplex* a, int* lda, int* jpvt, scomplex* tau, scomplex* work, int* lwork, float*  rwork, int* info )
  {
    return F77_cgeqp3( m, n, a, lda, jpvt, tau, work, lwork, rwork, info );
  }
  inline int geqp3( int* m, int* n, dcomplex* a, int* lda, int* jpvt, dcomplex* tau, dcomplex* work, int* lwork, double* rwork, int* info )
  {
    return F77_zgeqp3( m, n, a, lda, jpvt, tau, work, lwork, rwork, info );
  }

  // --- LQ factorization (classic) ---
  inline int gelqf( int* m, int* n, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    return F77_sgelqf( m, n, a, lda, tau, work, lwork, info );
  }
  inline int gelqf( int* m, int* n, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    return F77_dgelqf( m, n, a, lda, tau, work, lwork, info );
  }
  inline int gelqf( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    return F77_cgelqf( m, n, a, lda, tau, work, lwork, info );
  }
  inline int gelqf( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    return F77_zgelqf( m, n, a, lda, tau, work, lwork, info );
  }

  inline int gelq2( int* m, int* n, float* a, int* lda, float* tau, float* work, int* info )
  {
    return F77_sgelq2( m, n, a, lda, tau, work, info );
  }
  inline int gelq2( int* m, int* n, double*    a, int* lda, double*    tau, double*    work, int* info )
  {
    return F77_dgelq2( m, n, a, lda, tau, work, info );
  }
  inline int gelq2( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* info )
  {
    return F77_cgelq2( m, n, a, lda, tau, work, info );
  }
  inline int gelq2( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* info )
  {
    return F77_zgelq2( m, n, a, lda, tau, work, info );
  }

  // --- LS solver ---
  inline int gelsd( int* m, int* n, int* nrhs, float* a, int* lda, float* b, int* ldb, float*  s, float*  rcond, int* rank, float* work, int* lwork,     int* iwork, int* info )
  {
    return F77_sgelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info );
  }
  inline int gelsd( int* m, int* n, int* nrhs, double*    a, int* lda, double*    b, int* ldb, double* s, double* rcond, int* rank, double*    work, int* lwork,     int* iwork, int* info )
  {
    return F77_dgelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info );
  }
  inline int gelsd( int* m, int* n, int* nrhs, scomplex* a, int* lda, scomplex* b, int* ldb, float*  s, float*  rcond, int* rank, scomplex* work, int* lwork, float*  rwork, int* iwork, int* info )
  {
    return F77_cgelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info );
  }
  inline int gelsd( int* m, int* n, int* nrhs, dcomplex* a, int* lda, dcomplex* b, int* ldb, double* s, double* rcond, int* rank, dcomplex* work, int* lwork, double* rwork, int* iwork, int* info )
  {
    return F77_zgelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info );
  }

  inline int gelss( int* m, int* n, int* nrhs, float* a, int* lda, float* b, int* ldb, float*  s, float*  rcond, int* rank, float* work, int* lwork,     int* info )
  {
    return F77_sgelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info );
  }
  inline int gelss( int* m, int* n, int* nrhs, double*    a, int* lda, double*    b, int* ldb, double* s, double* rcond, int* rank, double*    work, int* lwork,     int* info )
  {
    return F77_dgelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info );
  }
  inline int gelss( int* m, int* n, int* nrhs, scomplex* a, int* lda, scomplex* b, int* ldb, float*  s, float*  rcond, int* rank, scomplex* work, int* lwork, float*  rwork, int* info )
  {
    return F77_cgelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info );
  }
  inline int gelss( int* m, int* n, int* nrhs, dcomplex* a, int* lda, dcomplex* b, int* ldb, double* s, double* rcond, int* rank, dcomplex* work, int* lwork, double* rwork, int* info )
  {
    return F77_zgelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info );
  }

  // --- Triangular-transpose matrix multiply ---
  inline int lauum( char* uplo, int* n, float* a, int* lda, int* info )
  {
    return F77_slauum( uplo, n, a, lda, info );
  }
  inline int lauum( char* uplo, int* n, double*    a, int* lda, int* info )
  {
    return F77_dlauum( uplo, n, a, lda, info );
  }
  inline int lauum( char* uplo, int* n, scomplex* a, int* lda, int* info )
  {
    return F77_clauum( uplo, n, a, lda, info );
  }
  inline int lauum( char* uplo, int* n, dcomplex* a, int* lda, int* info )
  {
    return F77_zlauum( uplo, n, a, lda, info );
  }

  inline int lauu2( char* uplo, int* n, float* a, int* lda, int* info )
  {
    return F77_slauu2( uplo, n, a, lda, info );
  }
  inline int lauu2( char* uplo, int* n, double*    a, int* lda, int* info )
  {
    return F77_dlauu2( uplo, n, a, lda, info );
  }
  inline int lauu2( char* uplo, int* n, scomplex* a, int* lda, int* info )
  {
    return F77_clauu2( uplo, n, a, lda, info );
  }
  inline int lauu2( char* uplo, int* n, dcomplex* a, int* lda, int* info )
  {
    return F77_zlauu2( uplo, n, a, lda, info );
  }

  // --- Symmetric (hermitian) positive definite matrix inversion ---
  inline int potri( char* uplo, int*  n, float* buff_A, int*  ldim_A, int*  info )
  {
    return F77_spotri( uplo, n, buff_A, ldim_A, info );
  }
  inline int potri( char* uplo, int*  n, double*    buff_A, int*  ldim_A, int*  info )
  {
    return F77_dpotri( uplo, n, buff_A, ldim_A, info );
  }
  inline int potri( char* uplo, int*  n, scomplex* buff_A, int*  ldim_A, int*  info )
  {
    return F77_cpotri( uplo, n, buff_A, ldim_A, info );
  }
  inline int potri( char* uplo, int*  n, dcomplex* buff_A, int*  ldim_A, int*  info )
  {
    return F77_zpotri( uplo, n, buff_A, ldim_A, info );
  }

  // --- Triangular matrix inversion ---
  inline int trtri( char* uplo, char* diag, int* n, float* a, int* lda, int* info )
  {
    return F77_strtri( uplo, diag, n, a, lda, info );
  }
  inline int trtri( char* uplo, char* diag, int* n, double*    a, int* lda, int* info )
  {
    return F77_dtrtri( uplo, diag, n, a, lda, info );
  }
  inline int trtri( char* uplo, char* diag, int* n, scomplex* a, int* lda, int* info )
  {
    return F77_ctrtri( uplo, diag, n, a, lda, info );
  }
  inline int trtri( char* uplo, char* diag, int* n, dcomplex* a, int* lda, int* info )
  {
    return F77_ztrtri( uplo, diag, n, a, lda, info );
  }

  inline int trti2( char* uplo, char* diag, int* n, float* a, int* lda, int* info )
  {
    return F77_strti2( uplo, diag, n, a, lda, info );
  }
  inline int trti2( char* uplo, char* diag, int* n, double*    a, int* lda, int* info )
  {
    return F77_dtrti2( uplo, diag, n, a, lda, info );
  }
  inline int trti2( char* uplo, char* diag, int* n, scomplex* a, int* lda, int* info )
  {
    return F77_ctrti2( uplo, diag, n, a, lda, info );
  }
  inline int trti2( char* uplo, char* diag, int* n, dcomplex* a, int* lda, int* info )
  {
    return F77_ztrti2( uplo, diag, n, a, lda, info );
  }

  // --- Triangular Sylvester equation solve ---
  inline int trsyl( char* transa, char* transb, int* isgn, int* m, int* n, float* a, int* lda, float* b, int* ldb, float* c, int* ldc, float* scale, int* info )
  {
    return F77_strsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }
  inline int trsyl( char* transa, char* transb, int* isgn, int* m, int* n, double*    a, int* lda, double*    b, int* ldb, double*    c, int* ldc, double*    scale, int* info )
  {
    return F77_dtrsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }
  inline int trsyl( char* transa, char* transb, int* isgn, int* m, int* n, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* c, int* ldc, float* scale, int* info )
  {
    return F77_ctrsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }
  inline int trsyl( char* transa, char* transb, int* isgn, int* m, int* n, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* c, int* ldc, double*    scale, int* info )
  {
    return F77_ztrsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }

  // --- Reduction to upper Hessenberg form ---
  inline int gehrd( int* n, int* ilo, int* ihi, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    return F77_sgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
  }
  inline int gehrd( int* n, int* ilo, int* ihi, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    return F77_dgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
  }
  inline int gehrd( int* n, int* ilo, int* ihi, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    return F77_cgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
  }
  inline int gehrd( int* n, int* ilo, int* ihi, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    return F77_zgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
  }

  inline int gehd2( int* n, int* ilo, int* ihi, float* a, int* lda, float* tau, float* work, int* info )
  {
    return F77_sgehd2( n, ilo, ihi, a, lda, tau, work, info );
  }
  inline int gehd2( int* n, int* ilo, int* ihi, double*    a, int* lda, double*    tau, double*    work, int* info )
  {
    return F77_dgehd2( n, ilo, ihi, a, lda, tau, work, info );
  }
  inline int gehd2( int* n, int* ilo, int* ihi, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* info )
  {
    return F77_cgehd2( n, ilo, ihi, a, lda, tau, work, info );
  }
  inline int gehd2( int* n, int* ilo, int* ihi, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* info )
  {
    return F77_zgehd2( n, ilo, ihi, a, lda, tau, work, info );
  }

  // --- Reduction to tridiagonal form ---
  inline int sytrd( char* uplo, int* n, float* a, int* lda, float*  d, float*  e, float* tau, float* work, int* lwork, int* info )
  {
    return F77_ssytrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
  }
  inline int sytrd( char* uplo, int* n, double*    a, int* lda, double* d, double* e, double*    tau, double*    work, int* lwork, int* info )
  {
    return F77_dsytrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
  }
  inline int hetrd( char* uplo, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    return F77_chetrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
  }
  inline int hetrd( char* uplo, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    return F77_zhetrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
  }

  inline int sytd2( char* uplo, int* n, float* a, int* lda, float*  d, float*  e, float* tau, int* info )
  {
    return F77_ssytd2( uplo, n, a, lda, d, e, tau, info );
  }
  inline int sytd2( char* uplo, int* n, double*    a, int* lda, double* d, double* e, double*    tau, int* info )
  {
    return F77_dsytd2( uplo, n, a, lda, d, e, tau, info );
  }
  inline int hetd2( char* uplo, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tau, int* info )
  {
    return F77_chetd2( uplo, n, a, lda, d, e, tau, info );
  }
  inline int hetd2( char* uplo, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tau, int* info )
  {
    return F77_zhetd2( uplo, n, a, lda, d, e, tau, info );
  }

  // --- Reduction to bidiagonal form ---
  inline int gebrd( int* m, int* n, float* a, int* lda, float*  d, float*  e, float* tauq, float* taup, float* work, int* lwork, int* info )
  {
    return F77_sgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }
  inline int gebrd( int* m, int* n, double*    a, int* lda, double* d, double* e, double*    tauq, double*    taup, double*    work, int* lwork, int* info )
  {
    return F77_dgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }
  inline int gebrd( int* m, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, int* lwork, int* info )
  {
    return F77_cgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }
  inline int gebrd( int* m, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, int* lwork, int* info )
  {
    return F77_zgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }

  inline int gebd2( int* m, int* n, float* a, int* lda, float*  d, float*  e, float* tauq, float* taup, float* work, int* info )
  {
    return F77_sgebd2( m, n, a, lda, d, e, tauq, taup, work, info );
  }
  inline int gebd2( int* m, int* n, double*    a, int* lda, double* d, double* e, double*    tauq, double*    taup, double*    work, int* info )
  {
    return F77_dgebd2( m, n, a, lda, d, e, tauq, taup, work, info );
  }
  inline int gebd2( int* m, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, int* info )
  {
    return F77_cgebd2( m, n, a, lda, d, e, tauq, taup, work, info );
  }
  inline int gebd2( int* m, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, int* info )
  {
    return F77_zgebd2( m, n, a, lda, d, e, tauq, taup, work, info );
  }

  // --- Reduce Hermitian-definite generalized eigenproblem to standard form ---
  inline int sygst( int* itype, char* uplo, int* n, float* a, int* lda, float* b, int* ldb, int* info )
  {
    return F77_ssygst( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline int sygst( int* itype, char* uplo, int* n, double*    a, int* lda, double*    b, int* ldb, int* info )
  {
    return F77_dsygst( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline int hegst( int* itype, char* uplo, int* n, scomplex* a, int* lda, scomplex* b, int* ldb, int* info )
  {
    return F77_chegst( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline int hegst( int* itype, char* uplo, int* n, dcomplex* a, int* lda, dcomplex* b, int* ldb, int* info )
  {
    return F77_zhegst( itype, uplo, n, a, lda, b, ldb, info );
  }

  inline int sygs2( int* itype, char* uplo, int* n, float* a, int* lda, float* b, int* ldb, int* info )
  {
    return F77_ssygs2( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline int sygs2( int* itype, char* uplo, int* n, double*    a, int* lda, double*    b, int* ldb, int* info )
  {
    return F77_dsygs2( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline int hegs2( int* itype, char* uplo, int* n, scomplex* a, int* lda, scomplex* b, int* ldb, int* info )
  {
    return F77_chegs2( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline int hegs2( int* itype, char* uplo, int* n, dcomplex* a, int* lda, dcomplex* b, int* ldb, int* info )
  {
    return F77_zhegs2( itype, uplo, n, a, lda, b, ldb, info );
  }

  // --- Accumulate block Householder matrix T (classic) ---
  inline int larft( char* direct, char* storev, int* n, int* k, float* v, int* ldv, float* tau, float* t, int* ldt )
  {
    return F77_slarft( direct, storev, n, k, v, ldv, tau, t, ldt );
  }
  inline int larft( char* direct, char* storev, int* n, int* k, double*    v, int* ldv, double*    tau, double*    t, int* ldt )
  {
    return F77_dlarft( direct, storev, n, k, v, ldv, tau, t, ldt );
  }
  inline int larft( char* direct, char* storev, int* n, int* k, scomplex* v, int* ldv, scomplex* tau, scomplex* t, int* ldt )
  {
    return F77_clarft( direct, storev, n, k, v, ldv, tau, t, ldt );
  }
  inline int larft( char* direct, char* storev, int* n, int* k, dcomplex* v, int* ldv, dcomplex* tau, dcomplex* t, int* ldt )
  {
    return F77_zlarft( direct, storev, n, k, v, ldv, tau, t, ldt );
  }

  // --- Generate a Householder vector (classic) ---
  inline int larfg( int* n, float* alpha, float* x, int* incx, float* tau )
  {
    return F77_slarfg( n, alpha, x, incx, tau );
  }
  inline int larfg( int* n, double*    alpha, double*    x, int* incx, double*    tau )
  {
    return F77_dlarfg( n, alpha, x, incx, tau );
  }
  inline int larfg( int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* tau )
  {
    return F77_clarfg( n, alpha, x, incx, tau );
  }
  inline int larfg( int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* tau )
  {
    return F77_zlarfg( n, alpha, x, incx, tau );
  }

  inline int larfgp( int* n, float* alpha, float* x, int* incx, float* tau )
  {
    return F77_slarfgp( n, alpha, x, incx, tau );
  }
  inline int larfgp( int* n, double*    alpha, double*    x, int* incx, double*    tau )
  {
    return F77_dlarfgp( n, alpha, x, incx, tau );
  }
  inline int larfgp( int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* tau )
  {
    return F77_clarfgp( n, alpha, x, incx, tau );
  }
  inline int larfgp( int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* tau )
  {
    return F77_zlarfgp( n, alpha, x, incx, tau );
  }

  // --- Form Q from QR factorization ---
  inline int orgqr( int* m, int* n, int* k, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
     return F77_sorgqr( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline int orgqr( int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    return F77_dorgqr( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline int ungqr( int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    work, int* lwork, int* info )
  {
    return F77_cungqr( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline int ungqr( int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    work, int* lwork, int* info )
  {
    return F77_zungqr( m, n, k, a, lda, tau, work, lwork, info );
  }


  // --- Apply Q or Q' from QR factorization ---
  inline int ormqr( char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* lwork, int* info )
  {
    return F77_sormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int ormqr( char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* lwork, int* info )
  {
    return F77_dormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int unmqr( char* side, char* trans, int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    c, int* ldc, scomplex*    work, int* lwork, int* info )
  {
    return F77_cunmqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int unmqr( char* side, char* trans, int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    c, int* ldc, dcomplex*    work, int* lwork, int* info )
  {
    return F77_zunmqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }

  inline int orm2r( char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* info )
  {
    return F77_sorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline int orm2r( char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* info )
  {
    return F77_dorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline int unm2r( char* side, char* trans, int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    c, int* ldc, scomplex*    work, int* info )
  {
    return F77_cunm2r( side, trans, m, n, k,  a, lda, tau, c, ldc, work, info );
  }
  inline int unm2r( char* side, char* trans, int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    c, int* ldc, dcomplex*    work, int* info )
  {
    return F77_zunm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }

  // --- Form Q from LQ factorization ---
  inline int orglq( int* m, int* n, int* k, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    return F77_sorglq( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline int orglq( int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    return F77_dorglq( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline int unglq( int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    return F77_cunglq( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline int unglq( int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    return F77_zunglq( m, n, k, a, lda, tau, work, lwork, info );
  }

  // --- Apply Q or Q' from LQ factorization ---
  inline int ormlq( char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* lwork, int* info )
  {
    return F77_sormlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int ormlq( char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* lwork, int* info )
  {
    return F77_dormlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int unmlq( char* side, char* trans, int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    c, int* ldc, scomplex*    work, int* lwork, int* info )
  {
    return F77_cunmlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int unmlq( char* side, char* trans, int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    c, int* ldc, dcomplex*    work, int* lwork, int* info )
  {
    return F77_zunmlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }

  inline int orml2( char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* info )
  {
    return F77_sorml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline int orml2( char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* info )
  {
    return F77_dorml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline int unml2( char* side, char* trans, int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    c, int* ldc, scomplex*    work, int* info )
  {
    return F77_cunml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline int unml2( char* side, char* trans, int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    c, int* ldc, dcomplex*    work, int* info )
  {
    return F77_zunml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }

  // --- Form Q from tridiagonal reduction ---
  inline int orgtr( char* uplo, int* m, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    return F77_sorgtr( uplo, m, a, lda, tau, work, lwork, info );
  }
  inline int orgtr( char* uplo, int* m, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    return F77_dorgtr( uplo, m, a, lda, tau, work, lwork, info );
  }
  inline int oungtr( char* uplo, int* m, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    return F77_cungtr( uplo, m, a, lda, tau, work, lwork, info );
  }
  inline int oungtr( char* uplo, int* m, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    return F77_zungtr( uplo, m, a, lda, tau, work, lwork, info );
  }

  // --- Apply Q or Q' from tridiagonal reduction ---
  inline int ormtr( char* side, char* uplo, char* trans, int* m, int* n, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* lwork, int* info )
  {
    return F77_sormtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int ormtr( char* side, char* uplo, char* trans, int* m, int* n, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* lwork, int* info )
  {
    return F77_dormtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int unmtr( char* side, char* uplo, char* trans, int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* lwork, int* info )
  {
    return F77_cunmtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int unmtr( char* side, char* uplo, char* trans, int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* lwork, int* info )
  {
    return F77_zunmtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
  }

  // --- Form Q from bidiagonal reduction ---
  inline int orgbr( char* vect, int* m, int* n, int* k, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    return F77_sorgbr( vect, m, n, k, a, lda, tau, work, lwork, info );
  }
  inline int orgbr( char* vect, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    return F77_dorgbr( vect, m, n, k, a, lda, tau, work, lwork, info );
  }
  inline int ungbr( char* vect, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    return F77_cungbr( vect, m, n, k, a, lda, tau, work, lwork, info );
  }
  inline int ungbr( char* vect, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    return F77_zungbr( vect, m, n, k, a, lda, tau, work, lwork, info );
  }

  // --- Apply Q or Q' from bidiagonal reduction ---
  inline int ormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* lwork, int* info )
  {
    return F77_sormbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int ormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* lwork, int* info )
  {
    return F77_dormbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int unmbr( char* vect, char* side, char* trans, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* lwork, int* info )
  {
    return F77_cunmbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline int unmbr( char* vect, char* side, char* trans, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* lwork, int* info )
  {
    return F77_zunmbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }

  // --- Tridiagonal QR algorithm ---
  inline int steqr( char* jobz, int* n, float* d, float* e, float* z, int* ldz, float*  work, int* info )
  {
    return F77_ssteqr( jobz, n, d, e, z, ldz, work, info );
  }
  inline int steqr( char* jobz, int* n, double*    d, double*    e, double*    z, int* ldz, double* work, int* info )
  {
    return F77_dsteqr( jobz, n, d, e, z, ldz, work, info );
  }
  inline int steqr( char* jobz, int* n, float* d, float* e, scomplex* z, int* ldz, float*  work, int* info )
  {
    return F77_csteqr( jobz, n, d, e, z, ldz, work, info );
  }
  inline int steqr( char* jobz, int* n, double*    d, double*    e, dcomplex* z, int* ldz, double* work, int* info )
  {
    return F77_zsteqr( jobz, n, d, e, z, ldz, work, info );
  }

  // --- Tridiagonal divide-and-conquer algorithm ---
  inline int stedc( char* compz, int* n, float* d, float* e, float* z, int* ldz, float* work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    return F77_sstedc( compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info );
  }
  inline int stedc( char* compz, int* n, double*    d, double*    e, double*    z, int* ldz, double*    work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    return F77_dstedc( compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info );
  }
  inline int stedc( char* compz, int* n, float* d, float* e, scomplex* z, int* ldz, scomplex* work, int* lwork, float*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    return F77_cstedc( compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info );
  }
  inline int stedc( char* compz, int* n, double*    d, double*    e, dcomplex* z, int* ldz, dcomplex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    return F77_zstedc( compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info );
  }

  // --- Tridiagonal MRRR algorithm ---
  inline int stemr( char* jobz, char* range, int* n, float*  d, float*  e, int* vl, int* vu, int* il, int* iu, int* m, float*  w, float* z, int* ldz, int* nzc, int* isuppz, int* tryrac, float*  work, int* lwork, int* iwork, int* liwork, int* info )
  {
    return F77_sstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }
  inline int stemr( char* jobz, char* range, int* n, double* d, double* e, int* vl, int* vu, int* il, int* iu, int* m, double* w, double*    z, int* ldz, int* nzc, int* isuppz, int* tryrac, double* work, int* lwork, int* iwork, int* liwork, int* info )
  {
    return F77_dstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }
  inline int stemr( char* jobz, char* range, int* n, float*  d, float*  e, int* vl, int* vu, int* il, int* iu, int* m, float*  w, scomplex* z, int* ldz, int* nzc, int* isuppz, int* tryrac, float*  work, int* lwork, int* iwork, int* liwork, int* info )
  {
    return F77_cstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }
  inline int stemr( char* jobz, char* range, int* n, double* d, double* e, int* vl, int* vu, int* il, int* iu, int* m, double* w, dcomplex* z, int* ldz, int* nzc, int* isuppz, int* tryrac, double* work, int* lwork, int* iwork, int* liwork, int* info )
  {
    return F77_zstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }

  // --- Hermitian eigenvalue decomposition (QR algorithm) ---
  inline int syev( char* jobz, char* uplo, int* n, float* a, int* lda, float*  w, float* work, int* lwork, float*  rwork, int* info )
  {
    return F77_ssyev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info );
  }
  inline int syev( char* jobz, char* uplo, int* n, double*    a, int* lda, double* w, double*    work, int* lwork, double* rwork, int* info )
  {
    return F77_dsyev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info );
  }
  inline int heev( char* jobz, char* uplo, int* n, scomplex* a, int* lda, float*  w, scomplex* work, int* lwork, float*  rwork, int* info )
  {
    return F77_cheev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info );
  }
  inline int heev( char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* info )
  {
    return F77_zheev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info );
  }

  // --- Hermitian eigenvalue decomposition (divide-and-conquer) ---
  inline int syevd( char* jobz, char* uplo, int* n, float* a, int* lda, float*  w, float* work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    return F77_ssyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info );
  }
  inline int syevd( char* jobz, char* uplo, int* n, double*    a, int* lda, double* w, double*    work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    return F77_dsyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info );
  }
  inline int heevd( char* jobz, char* uplo, int* n, scomplex* a, int* lda, float*  w, scomplex* work, int* lwork, float*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    return F77_cheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info );
  }
  inline int heevd( char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    return F77_zheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info );
  }

  // --- Hermitian eigenvalue decomposition (MRRR) ---
  inline int syevr( char* jobz, char* range, char* uplo, int* n, float* a, int* lda, float*  vl, float*  vu, int* il, int* iu, float*  abstol, int* m, float*  w, float* z, int* ldz, int* isuppz, float* work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    return F77_ssyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info );
  }
  inline int syevr( char* jobz, char* range, char* uplo, int* n, double*    a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, double*    z, int* ldz, int* isuppz, double*    work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    return F77_dsyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,  iwork, liwork, info );
  }
  inline int heevr( char* jobz, char* range, char* uplo, int* n, scomplex* a, int* lda, float*  vl, float*  vu, int* il, int* iu, float*  abstol, int* m, float*  w, scomplex* z, int* ldz, int* isuppz, scomplex* work, int* lwork, float*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    return F77_cheevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info );
  }
  inline int heevr( char* jobz, char* range, char* uplo, int* n, dcomplex* a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, dcomplex* z, int* ldz, int* isuppz, dcomplex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    return F77_zheevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info );
  }

  // --- Bidiagonal QR algorithm ---
  inline int bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, float* d, float* e, float* vt, int* ldvt, float* u, int* ldu, float* c, int* ldc, float*  rwork, int* info )
  {
    return F77_sbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info );
  }
  inline int bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, double*    d, double*    e, double*    vt, int* ldvt, double*    u, int* ldu, double*    c, int* ldc, double* rwork, int* info )
  {
    return F77_dbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info );
  }
  inline int bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, float* d, float* e, scomplex* vt, int* ldvt, scomplex* u, int* ldu, scomplex* c, int* ldc, float*  rwork, int* info )
  {
    return F77_cbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info );
  }
  inline int bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, double*    d, double*    e, dcomplex* vt, int* ldvt, dcomplex* u, int* ldu, dcomplex* c, int* ldc, double* rwork, int* info )
  {
    return F77_zbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info );
  }

  // --- Bidiagonal divide-and-conquor algorithm ---
  inline int bdsdc( char* uplo, char* compq, int* n, float*  d, float*  e, float*  u, int* ldu, float*  vt, int* ldvt, float*  q, float*  iq, float*  work, int* iwork, int* info )
  {
    return F77_sbdsdc( uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info );
  }
  inline int bdsdc( char* uplo, char* compq, int* n, double* d, double* e, double* u, int* ldu, double* vt, int* ldvt, double* q, double* iq, double* work, int* iwork, int* info )
  {
    return F77_dbdsdc( uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info );
  }

  // --- General matrix singular value decomposition (QR algorithm) ---
  inline int gesvd( char* jobu, char* jobv, int* m, int* n, float* a, int* lda, float*  s, float* u, int* ldu, float* vt, int* ldvt, float* work, int* lwork,     int* info )
  {
    return F77_sgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info );
  }
  inline int gesvd( char* jobu, char* jobv, int* m, int* n, double* a, int* lda, double* s, double*    u, int* ldu, double*    vt, int* ldvt, double*    work, int* lwork,     int* info )
  {
    return F77_dgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info );
  }
  inline int gesvd( char* jobu, char* jobv, int* m, int* n, scomplex* a, int* lda, float*  s, scomplex* u, int* ldu, scomplex* vt, int* ldvt, scomplex* work, int* lwork, float*  rwork, int* info )
  {
    return F77_cgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info );
  }
  inline int gesvd( char* jobu, char* jobv, int* m, int* n, dcomplex* a, int* lda, double* s, dcomplex* u, int* ldu, dcomplex* vt, int* ldvt, dcomplex* work, int* lwork, double* rwork, int* info )
  {
    return F77_zgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info );
  }

  // --- General matrix singular value decomposition (divide-and-conquer) ---
  inline int gesdd( char* jobz, int* m, int* n, float* a, int* lda, float*  s, float* u, int* ldu, float* vt, int* ldvt, float* work, int* lwork,     int* iwork, int* info )
  {
    return F77_sgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info );
  }
  inline int gesdd( char* jobz, int* m, int* n, double*    a, int* lda, double* s, double*    u, int* ldu, double*    vt, int* ldvt, double*    work, int* lwork,     int* iwork, int* info )
  {
    return F77_dgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info );
  }
  inline int gesdd( char* jobz, int* m, int* n, scomplex* a, int* lda, float*  s, scomplex* u, int* ldu, scomplex* vt, int* ldvt, scomplex* work, int* lwork, float*  rwork, int* iwork, int* info )
  {
    return F77_cgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info );
  }
  inline int gesdd( char* jobz, int* m, int* n, dcomplex* a, int* lda, double* s, dcomplex* u, int* ldu, dcomplex* vt, int* ldvt, dcomplex* work, int* lwork, double* rwork, int* iwork, int* info )
  {
    return F77_zgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info );
  }

  // --- Swap rows ---
  inline int laswp( int* n, float* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    return F77_slaswp( n, a, lda, k1, k2, ipiv, incx );
  }
  inline int laswp( int* n, double*    a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    return F77_dlaswp( n, a, lda, k1, k2, ipiv, incx );
  }
  inline int laswp( int* n, scomplex* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    return F77_claswp( n, a, lda, k1, k2, ipiv, incx );
  }
  inline int laswp( int* n, dcomplex* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    return F77_zlaswp( n, a, lda, k1, k2, ipiv, incx );
  }

  // --- Initialize a matrix ---
  inline int laset( char* uplo, int* m, int* n, float* alpha, float* beta, float* a, int* lda )
  {
    return F77_slaset( uplo, m, n, alpha, beta, a, lda );
  }
  inline int laset( char* uplo, int* m, int* n, double*    alpha, double*    beta, double*    a, int* lda )
  {
    return F77_dlaset( uplo, m, n, alpha, beta, a, lda );
  }
  inline int laset( char* uplo, int* m, int* n, scomplex* alpha, scomplex* beta, scomplex* a, int* lda )
  {
    return F77_claset( uplo, m, n, alpha, beta, a, lda );
  }
  inline int laset( char* uplo, int* m, int* n, dcomplex* alpha, dcomplex* beta, dcomplex* a, int* lda )
  {
    return F77_zlaset( uplo, m, n, alpha, beta, a, lda );
  }

}//namespace libflame

#endif  //  #ifndef LIBFLAME_HH

