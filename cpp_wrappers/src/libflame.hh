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
  inline void potrf( char* uplo, int* n, float* a, int* lda, int* info )
  {
    F77_spotrf(uplo, n, a, lda, info );
  }
  inline void potrf( char* uplo, int* n, double* a, int* lda, int* info )
  {
    F77_dpotrf(uplo, n, a, lda, info );
  }
  inline void potrf( char* uplo, int* n, scomplex* a, int* lda, int* info )
  {
    F77_cpotrf(uplo, n, a, lda, info );
  }

  inline void potrf( char* uplo, int* n, dcomplex* a, int* lda, int* info )
  {
    F77_zpotrf(uplo, n, a, lda, info );
  }  
  
  inline void potf2( char* uplo, int* n, float* a, int* lda, int* info )
  {
    F77_spotf2( uplo, n, a, lda, info );
  }
 
  inline void potf2( char* uplo, int* n, double*    a, int* lda, int* info )
  {
    F77_dpotf2( uplo, n, a, lda, info );
  }
  inline void potf2( char* uplo, int* n, scomplex* a, int* lda, int* info )
  {
    F77_cpotf2( uplo, n, a, lda, info );
  }
  inline void potf2( char* uplo, int* n, dcomplex* a, int* lda, int* info )
  {
    F77_zpotf2( uplo, n, a, lda, info );
  }

  // --- LU factorization with partial pivoting 
  inline void getrf( int* m, int* n, lapack_stype * buff_A, int* ldim_A, int* buff_p, int* info )
  {
    F77_sgetrf( m, n, buff_A, ldim_A, buff_p, info );    
  }
  inline void getrf( int* m, int* n, lapack_dtype* buff_A, int* ldim_A, int* buff_p, int* info)
  {
    F77_dgetrf( m, n, buff_A, ldim_A, buff_p, info );
  }
  inline void getrf( int* m, int* n, lapack_ctype* buff_A, int* ldim_A, int* buff_p, int* info )
  {
    F77_cgetrf( m, n, buff_A, ldim_A, buff_p, info );
  }
  inline void getrf( int* m, int* n, lapack_ztype* buff_A, int* ldim_A, int* buff_p, int* info )
  {
    F77_zgetrf( m, n, buff_A, ldim_A, buff_p, info );
  }  
  inline void getf2( int* m, int* n, float* a, int* lda, int* ipiv, int* info )
  {
    F77_sgetf2( m, n, a, lda, ipiv, info );
  }
  inline void getf2( int* m, int* n, double* a, int* lda, int* ipiv, int* info )
  {
    F77_dgetf2( m, n, a, lda, ipiv, info );
  }
  inline void getf2( int* m, int* n, scomplex* a, int* lda, int* ipiv, int* info )
  {
    F77_cgetf2( m, n, a, lda, ipiv, info );
  }  
  inline void getf2( int* m, int* n, dcomplex* a, int* lda, int* ipiv, int* info )
  {
    F77_zgetf2( m, n, a, lda, ipiv, info );
  }
  
  // --- QR factorization (classic) ---
  inline void geqrf( int* m, int* n, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    F77_sgeqrf( m, n, a, lda, tau, work, lwork, info );
  }
  inline void geqrf( int* m, int* n, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    F77_dgeqrf( m, n, a, lda, tau, work, lwork, info );
  }
  inline void geqrf( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    F77_cgeqrf( m, n, a, lda, tau, work, lwork, info );
  }
  inline void geqrf( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    F77_zgeqrf( m, n, a, lda, tau, work, lwork, info );
  }

  inline void geqr2( int* m, int* n, float* a, int* lda, float* tau, float* work, int* info )
  {
    F77_sgeqr2( m, n, a, lda, tau, work, info );
  }
  inline void geqr2( int* m, int* n, double*    a, int* lda, double*    tau, double*    work, int* info )
  {
    F77_dgeqr2( m, n, a, lda, tau, work, info );
  }
  inline void geqr2( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* info )
  {
    F77_cgeqr2( m, n, a, lda, tau, work, info );
  }
  inline void geqr2( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* info )
  {
    F77_zgeqr2( m, n, a, lda, tau, work, info );
  }


  inline void geqpf( int* m, int* n, float* a, int* lda, int* jpvt, float* tau, float* work,     int* info )
  {
    F77_sgeqpf( m, n, a, lda, jpvt, tau, work, info );
  }
  inline void geqpf( int* m, int* n, double*    a, int* lda, int* jpvt, double*    tau, double*    work,     int* info )
  {
    F77_dgeqpf( m, n, a, lda, jpvt, tau, work, info );
  }
  inline void geqpf( int* m, int* n, scomplex* a, int* lda, int* jpvt, scomplex* tau, scomplex* work, float*  rwork, int* info )
  {
    F77_cgeqpf( m, n, a, lda, jpvt, tau, work, rwork, info );
  }
  inline void geqpf( int* m, int* n, dcomplex* a, int* lda, int* jpvt, dcomplex* tau, dcomplex* work, double* rwork, int* info )
  {
    F77_zgeqpf( m, n, a, lda, jpvt, tau, work, rwork, info );
  }
  
  inline void geqp3( int* m, int* n, float* a, int* lda, int* jpvt, float* tau, float* work, int* lwork,     int* info )
  {
    F77_sgeqp3( m, n, a, lda, jpvt, tau, work, lwork, info );
  }
  inline void geqp3( int* m, int* n, double*    a, int* lda, int* jpvt, double*    tau, double*    work, int* lwork,     int* info )
  {
    F77_dgeqp3( m, n, a, lda, jpvt, tau, work, lwork, info );
  }
  inline void geqp3( int* m, int* n, scomplex* a, int* lda, int* jpvt, scomplex* tau, scomplex* work, int* lwork, float*  rwork, int* info )
  {
    F77_cgeqp3( m, n, a, lda, jpvt, tau, work, lwork, rwork, info );
  }
  inline void geqp3( int* m, int* n, dcomplex* a, int* lda, int* jpvt, dcomplex* tau, dcomplex* work, int* lwork, double* rwork, int* info )
  {
    F77_zgeqp3( m, n, a, lda, jpvt, tau, work, lwork, rwork, info );
  }
  
  // --- LQ factorization (classic) ---  
  inline void gelqf( int* m, int* n, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    F77_sgelqf( m, n, a, lda, tau, work, lwork, info );
  }
  inline void gelqf( int* m, int* n, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    F77_dgelqf( m, n, a, lda, tau, work, lwork, info );
  }
  inline void gelqf( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    F77_cgelqf( m, n, a, lda, tau, work, lwork, info );
  }
  inline void gelqf( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    F77_zgelqf( m, n, a, lda, tau, work, lwork, info );
  }
  
  inline void gelq2( int* m, int* n, float* a, int* lda, float* tau, float* work, int* info )
  {
    F77_sgelq2( m, n, a, lda, tau, work, info );
  }
  inline void gelq2( int* m, int* n, double*    a, int* lda, double*    tau, double*    work, int* info )
  {
    F77_dgelq2( m, n, a, lda, tau, work, info );
  }
  inline void gelq2( int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* info )
  {
    F77_cgelq2( m, n, a, lda, tau, work, info );
  }
  inline void gelq2( int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* info )
  {
    F77_zgelq2( m, n, a, lda, tau, work, info );
  }
  
  // --- LS solver ---  
  inline void gelsd( int* m, int* n, int* nrhs, float* a, int* lda, float* b, int* ldb, float*  s, float*  rcond, int* rank, float* work, int* lwork,     int* iwork, int* info )
  {
    F77_sgelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info );
  }
  inline void gelsd( int* m, int* n, int* nrhs, double*    a, int* lda, double*    b, int* ldb, double* s, double* rcond, int* rank, double*    work, int* lwork,     int* iwork, int* info )
  {
    F77_dgelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info );
  }
  inline void gelsd( int* m, int* n, int* nrhs, scomplex* a, int* lda, scomplex* b, int* ldb, float*  s, float*  rcond, int* rank, scomplex* work, int* lwork, float*  rwork, int* iwork, int* info )
  {
    F77_cgelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info );
  }
  inline void gelsd( int* m, int* n, int* nrhs, dcomplex* a, int* lda, dcomplex* b, int* ldb, double* s, double* rcond, int* rank, dcomplex* work, int* lwork, double* rwork, int* iwork, int* info )
  {
    F77_zgelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info );
  }
  
  inline void gelss( int* m, int* n, int* nrhs, float* a, int* lda, float* b, int* ldb, float*  s, float*  rcond, int* rank, float* work, int* lwork,     int* info )
  {
    F77_sgelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info );
  }
  inline void gelss( int* m, int* n, int* nrhs, double*    a, int* lda, double*    b, int* ldb, double* s, double* rcond, int* rank, double*    work, int* lwork,     int* info )
  {
    F77_dgelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info );
  }
  inline void gelss( int* m, int* n, int* nrhs, scomplex* a, int* lda, scomplex* b, int* ldb, float*  s, float*  rcond, int* rank, scomplex* work, int* lwork, float*  rwork, int* info )
  {
    F77_cgelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info );
  }
  inline void gelss( int* m, int* n, int* nrhs, dcomplex* a, int* lda, dcomplex* b, int* ldb, double* s, double* rcond, int* rank, dcomplex* work, int* lwork, double* rwork, int* info )
  {
    F77_zgelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info );
  }
  
  // --- Triangular-transpose matrix multiply ---
  inline void lauum( char* uplo, int* n, float* a, int* lda, int* info )
  {
    F77_slauum( uplo, n, a, lda, info );
  }
  inline void lauum( char* uplo, int* n, double*    a, int* lda, int* info )
  {
    F77_dlauum( uplo, n, a, lda, info );
  }
  inline void lauum( char* uplo, int* n, scomplex* a, int* lda, int* info )
  {
    F77_clauum( uplo, n, a, lda, info );
  }
  inline void lauum( char* uplo, int* n, dcomplex* a, int* lda, int* info )
  {
    F77_zlauum( uplo, n, a, lda, info );
  }
  
  inline void lauu2( char* uplo, int* n, float* a, int* lda, int* info )
  {
    F77_slauu2( uplo, n, a, lda, info );
  }
  inline void lauu2( char* uplo, int* n, double*    a, int* lda, int* info )
  {
    F77_dlauu2( uplo, n, a, lda, info );
  }
  inline void lauu2( char* uplo, int* n, scomplex* a, int* lda, int* info )
  {
    F77_clauu2( uplo, n, a, lda, info );
  }
  inline void lauu2( char* uplo, int* n, dcomplex* a, int* lda, int* info )
  {
    F77_zlauu2( uplo, n, a, lda, info );
  }
  
  // --- Symmetric (hermitian) positive definite matrix inversion ---  
  inline void potri( char* uplo, int*  n, float* buff_A, int*  ldim_A, int*  info )
  {
    F77_spotri( uplo, n, buff_A, ldim_A, info );
  }
  inline void potri( char* uplo, int*  n, double*    buff_A, int*  ldim_A, int*  info )
  {
    F77_dpotri( uplo, n, buff_A, ldim_A, info );
  }
  inline void potri( char* uplo, int*  n, scomplex* buff_A, int*  ldim_A, int*  info )
  {
    F77_cpotri( uplo, n, buff_A, ldim_A, info );
  }
  inline void potri( char* uplo, int*  n, dcomplex* buff_A, int*  ldim_A, int*  info )
  {
    F77_zpotri( uplo, n, buff_A, ldim_A, info );
  }
  
  // --- Triangular matrix inversion ---  
  inline void trtri( char* uplo, char* diag, int* n, float* a, int* lda, int* info )
  {
    F77_strtri( uplo, diag, n, a, lda, info );
  }
  inline void trtri( char* uplo, char* diag, int* n, double*    a, int* lda, int* info )
  {
    F77_dtrtri( uplo, diag, n, a, lda, info );
  }
  inline void trtri( char* uplo, char* diag, int* n, scomplex* a, int* lda, int* info )
  {
    F77_ctrtri( uplo, diag, n, a, lda, info );
  }
  inline void trtri( char* uplo, char* diag, int* n, dcomplex* a, int* lda, int* info )
  {
    F77_ztrtri( uplo, diag, n, a, lda, info );
  }
  
  inline void trti2( char* uplo, char* diag, int* n, float* a, int* lda, int* info )
  {
    F77_strti2( uplo, diag, n, a, lda, info );
  }
  inline void trti2( char* uplo, char* diag, int* n, double*    a, int* lda, int* info )
  {
    F77_dtrti2( uplo, diag, n, a, lda, info );
  }
  inline void trti2( char* uplo, char* diag, int* n, scomplex* a, int* lda, int* info )
  {
    F77_ctrti2( uplo, diag, n, a, lda, info );
  }
  inline void trti2( char* uplo, char* diag, int* n, dcomplex* a, int* lda, int* info )
  {
    F77_ztrti2( uplo, diag, n, a, lda, info );
  }
  
  // --- Triangular Sylvester equation solve ---  
  inline void trsyl( char* transa, char* transb, int* isgn, int* m, int* n, float* a, int* lda, float* b, int* ldb, float* c, int* ldc, float* scale, int* info )
  {
    F77_strsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }
  inline void trsyl( char* transa, char* transb, int* isgn, int* m, int* n, double*    a, int* lda, double*    b, int* ldb, double*    c, int* ldc, double*    scale, int* info )
  {
    F77_dtrsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }
  inline void trsyl( char* transa, char* transb, int* isgn, int* m, int* n, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* c, int* ldc, float* scale, int* info )
  {
    F77_ctrsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }
  inline void trsyl( char* transa, char* transb, int* isgn, int* m, int* n, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* c, int* ldc, double*    scale, int* info )
  {
    F77_ztrsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }
  
  // --- Reduction to upper Hessenberg form ---  
  inline void gehrd( int* n, int* ilo, int* ihi, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    F77_sgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
  }
  inline void gehrd( int* n, int* ilo, int* ihi, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    F77_dgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
  }
  inline void gehrd( int* n, int* ilo, int* ihi, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    F77_cgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
  }
  inline void gehrd( int* n, int* ilo, int* ihi, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    F77_zgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
  }
  
  inline void gehd2( int* n, int* ilo, int* ihi, float* a, int* lda, float* tau, float* work, int* info )
  {
    F77_sgehd2( n, ilo, ihi, a, lda, tau, work, info );
  }
  inline void gehd2( int* n, int* ilo, int* ihi, double*    a, int* lda, double*    tau, double*    work, int* info )
  {
    F77_dgehd2( n, ilo, ihi, a, lda, tau, work, info );
  }
  inline void gehd2( int* n, int* ilo, int* ihi, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* info )
  {
    F77_cgehd2( n, ilo, ihi, a, lda, tau, work, info );
  }
  inline void gehd2( int* n, int* ilo, int* ihi, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* info )
  {
    F77_zgehd2( n, ilo, ihi, a, lda, tau, work, info );
  }
  
  // --- Reduction to tridiagonal form ---  
  inline void sytrd( char* uplo, int* n, float* a, int* lda, float*  d, float*  e, float* tau, float* work, int* lwork, int* info )
  {
    F77_ssytrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
  }
  inline void sytrd( char* uplo, int* n, double*    a, int* lda, double* d, double* e, double*    tau, double*    work, int* lwork, int* info )
  {
    F77_dsytrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
  }
  inline void hetrd( char* uplo, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    F77_chetrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
  }
  inline void hetrd( char* uplo, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    F77_zhetrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
  }
  
  inline void sytd2( char* uplo, int* n, float* a, int* lda, float*  d, float*  e, float* tau, int* info )
  {
    F77_ssytd2( uplo, n, a, lda, d, e, tau, info );
  }
  inline void sytd2( char* uplo, int* n, double*    a, int* lda, double* d, double* e, double*    tau, int* info )
  {
    F77_dsytd2( uplo, n, a, lda, d, e, tau, info );
  }
  inline void hetd2( char* uplo, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tau, int* info )
  {
    F77_chetd2( uplo, n, a, lda, d, e, tau, info );
  }
  inline void hetd2( char* uplo, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tau, int* info )
  {
    F77_zhetd2( uplo, n, a, lda, d, e, tau, info );
  }
  
  // --- Reduction to bidiagonal form ---  
  inline void gebrd( int* m, int* n, float* a, int* lda, float*  d, float*  e, float* tauq, float* taup, float* work, int* lwork, int* info )
  {
    F77_sgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }
  inline void gebrd( int* m, int* n, double*    a, int* lda, double* d, double* e, double*    tauq, double*    taup, double*    work, int* lwork, int* info )
  {
    F77_dgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }
  inline void gebrd( int* m, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, int* lwork, int* info )
  {
    F77_cgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }
  inline void gebrd( int* m, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, int* lwork, int* info )
  {
    F77_zgebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }
  
  inline void gebd2( int* m, int* n, float* a, int* lda, float*  d, float*  e, float* tauq, float* taup, float* work, int* info )
  {
    F77_sgebd2( m, n, a, lda, d, e, tauq, taup, work, info );
  }
  inline void gebd2( int* m, int* n, double*    a, int* lda, double* d, double* e, double*    tauq, double*    taup, double*    work, int* info )
  {
    F77_dgebd2( m, n, a, lda, d, e, tauq, taup, work, info );
  }
  inline void gebd2( int* m, int* n, scomplex* a, int* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, int* info )
  {
    F77_cgebd2( m, n, a, lda, d, e, tauq, taup, work, info );
  }
  inline void gebd2( int* m, int* n, dcomplex* a, int* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, int* info )
  {
    F77_zgebd2( m, n, a, lda, d, e, tauq, taup, work, info );
  }
  
  // --- Reduce Hermitian-definite generalized eigenproblem to standard form ---  
  inline void sygst( int* itype, char* uplo, int* n, float* a, int* lda, float* b, int* ldb, int* info )
  {
    F77_ssygst( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline void sygst( int* itype, char* uplo, int* n, double*    a, int* lda, double*    b, int* ldb, int* info )
  {
    F77_dsygst( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline void sygst( int* itype, char* uplo, int* n, scomplex* a, int* lda, scomplex* b, int* ldb, int* info )
  {
    F77_chegst( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline void sygst( int* itype, char* uplo, int* n, dcomplex* a, int* lda, dcomplex* b, int* ldb, int* info )
  {
    F77_zhegst( itype, uplo, n, a, lda, b, ldb, info );
  }
  
  inline void sygs2( int* itype, char* uplo, int* n, float* a, int* lda, float* b, int* ldb, int* info )
  {
    F77_ssygs2( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline void sygs2( int* itype, char* uplo, int* n, double*    a, int* lda, double*    b, int* ldb, int* info )
  {
    F77_dsygs2( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline void hegs2( int* itype, char* uplo, int* n, scomplex* a, int* lda, scomplex* b, int* ldb, int* info )
  {
    F77_chegs2( itype, uplo, n, a, lda, b, ldb, info );
  }
  inline void hegs2( int* itype, char* uplo, int* n, dcomplex* a, int* lda, dcomplex* b, int* ldb, int* info )
  {
    F77_zhegs2( itype, uplo, n, a, lda, b, ldb, info );
  }
  
  // --- Accumulate block Householder matrix T (classic) ---  
  inline void larft( char* direct, char* storev, int* n, int* k, float* v, int* ldv, float* tau, float* t, int* ldt )
  {
    F77_slarft( direct, storev, n, k, v, ldv, tau, t, ldt );
  }
  inline void larft( char* direct, char* storev, int* n, int* k, double*    v, int* ldv, double*    tau, double*    t, int* ldt )
  {
    F77_dlarft( direct, storev, n, k, v, ldv, tau, t, ldt );
  }
  inline void larft( char* direct, char* storev, int* n, int* k, scomplex* v, int* ldv, scomplex* tau, scomplex* t, int* ldt )
  {
    F77_clarft( direct, storev, n, k, v, ldv, tau, t, ldt );
  }
  inline void larft( char* direct, char* storev, int* n, int* k, dcomplex* v, int* ldv, dcomplex* tau, dcomplex* t, int* ldt )
  {
    F77_zlarft( direct, storev, n, k, v, ldv, tau, t, ldt );
  }
  
  // --- Generate a Householder vector (classic) ---  
  inline void larfg( int* n, float* alpha, float* x, int* incx, float* tau )
  {
    F77_slarfg( n, alpha, x, incx, tau );
  }
  inline void larfg( int* n, double*    alpha, double*    x, int* incx, double*    tau )
  {
    F77_dlarfg( n, alpha, x, incx, tau );
  }
  inline void larfg( int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* tau )
  {
    F77_clarfg( n, alpha, x, incx, tau );
  }
  inline void larfg( int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* tau )
  {
    F77_zlarfg( n, alpha, x, incx, tau );
  }
  
  inline void larfgp( int* n, float* alpha, float* x, int* incx, float* tau )
  {
    F77_slarfgp( n, alpha, x, incx, tau );
  }
  inline void larfgp( int* n, double*    alpha, double*    x, int* incx, double*    tau )
  {
    F77_dlarfgp( n, alpha, x, incx, tau );
  }
  inline void larfgp( int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* tau )
  {
    F77_clarfgp( n, alpha, x, incx, tau );
  }
  inline void larfgp( int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* tau )
  {
    F77_zlarfgp( n, alpha, x, incx, tau );
  }
  
  // --- Form Q from QR factorization ---  
  inline void orgqr( int* m, int* n, int* k, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
     F77_sorgqr( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline void orgqr( int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    F77_dorgqr( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline void orgqr( int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    work, int* lwork, int* info )
  {
    F77_cungqr( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline void orgqr( int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    work, int* lwork, int* info )
  {
    F77_zungqr( m, n, k, a, lda, tau, work, lwork, info );
  }
  
  
  // --- Apply Q or Q' from QR factorization ---  
  inline void sormqr( char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* lwork, int* info )
  {
    F77_sormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void dormqr( char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* lwork, int* info )
  {
    F77_dormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void cunmqr( char* side, char* trans, int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    c, int* ldc, scomplex*    work, int* lwork, int* info )
  {
    F77_cunmqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void zunmqr( char* side, char* trans, int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    c, int* ldc, dcomplex*    work, int* lwork, int* info )
  {
    F77_zunmqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  
  inline void orm2r( char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* info )
  {
    F77_sorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline void orm2r( char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* info )
  {
    F77_dorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline void orm2r( char* side, char* trans, int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    c, int* ldc, scomplex*    work, int* info )
  {
    F77_cunm2r( side, trans, m, n, k,  a, lda, tau, c, ldc, work, info );
  }
  inline void orm2r( char* side, char* trans, int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    c, int* ldc, dcomplex*    work, int* info )
  {
    F77_zunm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  
  // --- Form Q from LQ factorization ---  
  inline void orglq( int* m, int* n, int* k, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    F77_sorglq( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline void orglq( int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    F77_dorglq( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline void orglq( int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    F77_cunglq( m, n, k, a, lda, tau, work, lwork, info );
  }
  inline void orglq( int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    F77_zunglq( m, n, k, a, lda, tau, work, lwork, info );
  }
  
  // --- Apply Q or Q' from LQ factorization ---  
  inline void ormlq( char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* lwork, int* info )
  {
    F77_sormlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void ormlq( char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* lwork, int* info )
  {
    F77_dormlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void ormlq( char* side, char* trans, int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    c, int* ldc, scomplex*    work, int* lwork, int* info )
  {
    F77_cunmlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void ormlq( char* side, char* trans, int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    c, int* ldc, dcomplex*    work, int* lwork, int* info )
  {
    F77_zunmlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  
  inline void orml2( char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* info )
  {
    F77_sorml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline void orml2( char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* info )
  {
    F77_dorml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline void orml2( char* side, char* trans, int* m, int* n, int* k, scomplex*    a, int* lda, scomplex*    tau, scomplex*    c, int* ldc, scomplex*    work, int* info )
  {
    F77_cunml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  inline void orml2( char* side, char* trans, int* m, int* n, int* k, dcomplex*    a, int* lda, dcomplex*    tau, dcomplex*    c, int* ldc, dcomplex*    work, int* info )
  {
    F77_zunml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  
  // --- Form Q from tridiagonal reduction ---  
  inline void orgtr( char* uplo, int* m, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    F77_sorgtr( uplo, m, a, lda, tau, work, lwork, info );
  }
  inline void orgtr( char* uplo, int* m, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    F77_dorgtr( uplo, m, a, lda, tau, work, lwork, info );
  }
  inline void orgtr( char* uplo, int* m, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    F77_cungtr( uplo, m, a, lda, tau, work, lwork, info );
  }
  inline void orgtr( char* uplo, int* m, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    F77_zungtr( uplo, m, a, lda, tau, work, lwork, info );
  }
  
  // --- Apply Q or Q' from tridiagonal reduction ---  
  inline void ormtr( char* side, char* uplo, char* trans, int* m, int* n, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* lwork, int* info )
  {
    F77_sormtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void ormtr( char* side, char* uplo, char* trans, int* m, int* n, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* lwork, int* info )
  {
    F77_dormtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void ormtr( char* side, char* uplo, char* trans, int* m, int* n, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* lwork, int* info )
  {
    F77_cunmtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void ormtr( char* side, char* uplo, char* trans, int* m, int* n, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* lwork, int* info )
  {
    F77_zunmtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
  }
  
  // --- Form Q from bidiagonal reduction ---  
  inline void orgbr( char* vect, int* m, int* n, int* k, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
  {
    F77_sorgbr( vect, m, n, k, a, lda, tau, work, lwork, info );
  }
  inline void orgbr( char* vect, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    work, int* lwork, int* info )
  {
    F77_dorgbr( vect, m, n, k, a, lda, tau, work, lwork, info );
  }
  inline void orgbr( char* vect, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* work, int* lwork, int* info )
  {
    F77_cungbr( vect, m, n, k, a, lda, tau, work, lwork, info );
  }
  inline void orgbr( char* vect, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info )
  {
    F77_zungbr( vect, m, n, k, a, lda, tau, work, lwork, info );
  }
  
  // --- Apply Q or Q' from bidiagonal reduction ---  
  inline void ormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* lwork, int* info )
  {
    F77_sormbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void ormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, double*    a, int* lda, double*    tau, double*    c, int* ldc, double*    work, int* lwork, int* info )
  {
    F77_dormbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void ormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, scomplex* a, int* lda, scomplex* tau, scomplex* c, int* ldc, scomplex* work, int* lwork, int* info )
  {
    F77_cunmbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  inline void ormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, dcomplex* a, int* lda, dcomplex* tau, dcomplex* c, int* ldc, dcomplex* work, int* lwork, int* info )
  {
    F77_zunmbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  
  // --- Tridiagonal QR algorithm ---  
  inline void steqr( char* jobz, int* n, float* d, float* e, float* z, int* ldz, float*  work, int* info )
  {
    F77_ssteqr( jobz, n, d, e, z, ldz, work, info );
  }
  inline void steqr( char* jobz, int* n, double*    d, double*    e, double*    z, int* ldz, double* work, int* info )
  {
    F77_dsteqr( jobz, n, d, e, z, ldz, work, info ); 
  }
  inline void steqr( char* jobz, int* n, float* d, float* e, scomplex* z, int* ldz, float*  work, int* info )
  {
    F77_csteqr( jobz, n, d, e, z, ldz, work, info ); 
  }
  inline void steqr( char* jobz, int* n, double*    d, double*    e, dcomplex* z, int* ldz, double* work, int* info )
  {
    F77_zsteqr( jobz, n, d, e, z, ldz, work, info ); 
  }
  
  // --- Tridiagonal divide-and-conquer algorithm ---
  inline void stedc( char* compz, int* n, float* d, float* e, float* z, int* ldz, float* work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    F77_sstedc( compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info );
  }
  inline void stedc( char* compz, int* n, double*    d, double*    e, double*    z, int* ldz, double*    work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    F77_dstedc( compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info );
  }
  inline void stedc( char* compz, int* n, float* d, float* e, scomplex* z, int* ldz, scomplex* work, int* lwork, float*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    F77_cstedc( compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info );
  }
  inline void stedc( char* compz, int* n, double*    d, double*    e, dcomplex* z, int* ldz, dcomplex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    F77_zstedc( compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info );
  }
  
  // --- Tridiagonal MRRR algorithm ---  
  inline void stemr( char* jobz, char* range, int* n, float*  d, float*  e, int* vl, int* vu, int* il, int* iu, int* m, float*  w, float* z, int* ldz, int* nzc, int* isuppz, int* tryrac, float*  work, int* lwork, int* iwork, int* liwork, int* info )
  {
    F77_sstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }
  inline void stemr( char* jobz, char* range, int* n, double* d, double* e, int* vl, int* vu, int* il, int* iu, int* m, double* w, double*    z, int* ldz, int* nzc, int* isuppz, int* tryrac, double* work, int* lwork, int* iwork, int* liwork, int* info )
  {
    F77_dstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }
  inline void stemr( char* jobz, char* range, int* n, float*  d, float*  e, int* vl, int* vu, int* il, int* iu, int* m, float*  w, scomplex* z, int* ldz, int* nzc, int* isuppz, int* tryrac, float*  work, int* lwork, int* iwork, int* liwork, int* info )
  {
    F77_cstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }
  inline void stemr( char* jobz, char* range, int* n, double* d, double* e, int* vl, int* vu, int* il, int* iu, int* m, double* w, dcomplex* z, int* ldz, int* nzc, int* isuppz, int* tryrac, double* work, int* lwork, int* iwork, int* liwork, int* info )
  {
    F77_zstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }
  
  // --- Hermitian eigenvalue decomposition (QR algorithm) ---  
  inline void syev( char* jobz, char* uplo, int* n, float* a, int* lda, float*  w, float* work, int* lwork, float*  rwork, int* info )
  {
    F77_ssyev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info ); 
  }
  inline void syev( char* jobz, char* uplo, int* n, double*    a, int* lda, double* w, double*    work, int* lwork, double* rwork, int* info )
  {
    F77_dsyev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info ); 
  }
  inline void syev( char* jobz, char* uplo, int* n, scomplex* a, int* lda, float*  w, scomplex* work, int* lwork, float*  rwork, int* info )
  {
    F77_cheev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info ); 
  }
  inline void syev( char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* info )
  {
    F77_zheev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info ); 
  }
  
  // --- Hermitian eigenvalue decomposition (divide-and-conquer) ---  
  inline void syevd( char* jobz, char* uplo, int* n, float* a, int* lda, float*  w, float* work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    F77_ssyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info ); 
  }
  inline void syevd( char* jobz, char* uplo, int* n, double*    a, int* lda, double* w, double*    work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    F77_dsyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info ); 
  }
  inline void syevd( char* jobz, char* uplo, int* n, scomplex* a, int* lda, float*  w, scomplex* work, int* lwork, float*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    F77_cheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info ); 
  }
  inline void syevd( char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    F77_zheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info ); 
  }
  
  // --- Hermitian eigenvalue decomposition (MRRR) ---
  inline void syevr( char* jobz, char* range, char* uplo, int* n, float* a, int* lda, float*  vl, float*  vu, int* il, int* iu, float*  abstol, int* m, float*  w, float* z, int* ldz, int* isuppz, float* work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    F77_ssyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info ); 
  }
  inline void syevr( char* jobz, char* range, char* uplo, int* n, double*    a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, double*    z, int* ldz, int* isuppz, double*    work, int* lwork,          int* iwork, int* liwork, int* info )
  {
    F77_dsyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,  iwork, liwork, info );
  }
  inline void syevr( char* jobz, char* range, char* uplo, int* n, scomplex* a, int* lda, float*  vl, float*  vu, int* il, int* iu, float*  abstol, int* m, float*  w, scomplex* z, int* ldz, int* isuppz, scomplex* work, int* lwork, float*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    F77_cheevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info ); 
  }
  inline void syevr( char* jobz, char* range, char* uplo, int* n, dcomplex* a, int* lda, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, double* w, dcomplex* z, int* ldz, int* isuppz, dcomplex* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    F77_zheevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info ); 
  }
  
  // --- Bidiagonal QR algorithm ---  
  inline void bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, float* d, float* e, float* vt, int* ldvt, float* u, int* ldu, float* c, int* ldc, float*  rwork, int* info )
  {
    F77_sbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info ); 
  }
  inline void bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, double*    d, double*    e, double*    vt, int* ldvt, double*    u, int* ldu, double*    c, int* ldc, double* rwork, int* info )
  {
    F77_dbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info ); 
  }
  inline void bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, float* d, float* e, scomplex* vt, int* ldvt, scomplex* u, int* ldu, scomplex* c, int* ldc, float*  rwork, int* info )
  {
    F77_cbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info ); 
  }
  inline void bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, double*    d, double*    e, dcomplex* vt, int* ldvt, dcomplex* u, int* ldu, dcomplex* c, int* ldc, double* rwork, int* info )
  {
    F77_zbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info ); 
  }
  
  // --- Bidiagonal divide-and-conquor algorithm ---  
  inline void bdsdc( char* uplo, char* compq, int* n, float*  d, float*  e, float*  u, int* ldu, float*  vt, int* ldvt, float*  q, float*  iq, float*  work, int* iwork, int* info )
  {
    F77_sbdsdc( uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info ); 
  }
  inline void bdsdc( char* uplo, char* compq, int* n, double* d, double* e, double* u, int* ldu, double* vt, int* ldvt, double* q, double* iq, double* work, int* iwork, int* info )
  {
    F77_dbdsdc( uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info ); 
  }
  
  // --- General matrix singular value decomposition (QR algorithm) ---  
  inline void gesvd( char* jobu, char* jobv, int* m, int* n, float* a, int* lda, float*  s, float* u, int* ldu, float* vt, int* ldvt, float* work, int* lwork,     int* info )
  {
    F77_sgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info );
  }
  inline void gesvd( char* jobu, char* jobv, int* m, int* n, double* a, int* lda, double* s, double*    u, int* ldu, double*    vt, int* ldvt, double*    work, int* lwork,     int* info )
  {
    F77_dgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info );
  }
  inline void gesvd( char* jobu, char* jobv, int* m, int* n, scomplex* a, int* lda, float*  s, scomplex* u, int* ldu, scomplex* vt, int* ldvt, scomplex* work, int* lwork, float*  rwork, int* info )
  {
    F77_cgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info );
  }
  inline void gesvd( char* jobu, char* jobv, int* m, int* n, dcomplex* a, int* lda, double* s, dcomplex* u, int* ldu, dcomplex* vt, int* ldvt, dcomplex* work, int* lwork, double* rwork, int* info )
  {
    F77_zgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info );
  }
  
  // --- General matrix singular value decomposition (divide-and-conquer) ---  
  inline void gesdd( char* jobz, int* m, int* n, float* a, int* lda, float*  s, float* u, int* ldu, float* vt, int* ldvt, float* work, int* lwork,     int* iwork, int* info )
  {
    F77_sgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info );
  }
  inline void gesdd( char* jobz, int* m, int* n, double*    a, int* lda, double* s, double*    u, int* ldu, double*    vt, int* ldvt, double*    work, int* lwork,     int* iwork, int* info )
  {
    F77_dgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info );
  }
  inline void gesdd( char* jobz, int* m, int* n, scomplex* a, int* lda, float*  s, scomplex* u, int* ldu, scomplex* vt, int* ldvt, scomplex* work, int* lwork, float*  rwork, int* iwork, int* info )
  {
    F77_cgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info );
  }
  inline void gesdd( char* jobz, int* m, int* n, dcomplex* a, int* lda, double* s, dcomplex* u, int* ldu, dcomplex* vt, int* ldvt, dcomplex* work, int* lwork, double* rwork, int* iwork, int* info )
  {
    F77_zgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info );
  }
  
  // --- Swap rows ---  
  inline void laswp( int* n, float* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    F77_slaswp( n, a, lda, k1, k2, ipiv, incx );
  }
  inline void laswp( int* n, double*    a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    F77_dlaswp( n, a, lda, k1, k2, ipiv, incx );
  }
  inline void laswp( int* n, scomplex* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    F77_claswp( n, a, lda, k1, k2, ipiv, incx );
  }
  inline void laswp( int* n, dcomplex* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    F77_zlaswp( n, a, lda, k1, k2, ipiv, incx );
  }
  
  // --- Initialize a matrix ---  
  inline void laset( char* uplo, int* m, int* n, float* alpha, float* beta, float* a, int* lda )
  {
    F77_slaset( uplo, m, n, alpha, beta, a, lda );
  }
  inline void laset( char* uplo, int* m, int* n, double*    alpha, double*    beta, double*    a, int* lda )
  {
    F77_dlaset( uplo, m, n, alpha, beta, a, lda );
  }
  inline void laset( char* uplo, int* m, int* n, scomplex* alpha, scomplex* beta, scomplex* a, int* lda )
  {
    F77_claset( uplo, m, n, alpha, beta, a, lda );
  }
  inline void laset( char* uplo, int* m, int* n, dcomplex* alpha, dcomplex* beta, dcomplex* a, int* lda )
  {
    F77_zlaset( uplo, m, n, alpha, beta, a, lda );
  }

}//namespace libflame

#endif  //  #ifndef LIBFLAME_HH

