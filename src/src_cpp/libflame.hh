/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame.hh
 *  libflame.hh defines all the overloaded CPP functions to be invoked from
 *  template interfaces
 *  */
#ifndef LIBFLAME_HH
#define LIBFLAME_HH

#include  <FLAME.h>

namespace libflame {

// Cholesky factorization of a real symmetric positive definite matrix a
inline integer potrf(char* uplo, integer* n, float* a, integer* lda, integer* info)
{
  return spotrf_(uplo, n, a, lda, info);
}
inline integer potrf(char* uplo, integer* n, double* a, integer* lda, integer* info)
{
  return dpotrf_(uplo, n, a, lda, info);
}
inline integer potrf(char* uplo, integer* n, scomplex* a, integer* lda, integer* info)
{
  return cpotrf_(uplo, n, a, lda, info);
}
inline integer potrf(char* uplo, integer* n, dcomplex* a, integer* lda, integer* info)
{
  return zpotrf_(uplo, n, a, lda, info);
}

// --- computes the Cholesky factorization of a symmetric/Hermitian positive definite matrix (unblocked algorithm) ---
inline integer potf2(char* uplo, integer*n, float* a, integer* lda, integer* info)
{
  return spotf2_(uplo, n, a, lda, info);
}
inline integer potf2(char* uplo, integer*n, double* a, integer* lda, integer* info)
{
  return dpotf2_(uplo, n, a, lda, info);
}
inline integer potf2(char* uplo, integer*n, scomplex* a, integer* lda, integer* info)
{
  return cpotf2_(uplo, n, a, lda, info);
}
inline integer potf2(char* uplo, integer*n, dcomplex* a, integer* lda, integer* info)
{
  return zpotf2_(uplo, n, a, lda, info);
}

// --- LU factorization with partial pivoting
inline integer getrf(integer* m, integer* n, float* a, integer* lda, integer* ipiv, integer* info)
{
  return sgetrf_(m, n, a, lda, ipiv, info);
}
inline integer getrf(integer* m, integer* n, double*   a, integer* lda, integer* ipiv, integer* info)
{
  return dgetrf_(m, n, a, lda, ipiv, info);
}
inline integer getrf(integer* m, integer* n, scomplex* a, integer* lda, integer* ipiv, integer* info)
{
  return cgetrf_(m, n, a, lda, ipiv, info);
}
inline integer getrf(integer* m, integer* n, dcomplex* a, integer* lda, integer* ipiv, integer* info)
{
  return zgetrf_(m, n, a, lda, ipiv, info);
}


//--- computes the LU factorization of a general m-by-n matrix using partial pivoting with row interchanges (unblocked algorithm) ---
inline integer getf2(integer* m, integer* n, float* a, integer* lda, integer* ipiv, integer* info)
{
  return sgetf2_(m, n, a, lda, ipiv, info);
}
inline integer getf2(integer* m, integer* n, double* a, integer* lda, integer* ipiv, integer* info)
{
  return dgetf2_(m, n, a, lda, ipiv, info);
}
inline integer getf2(integer* m, integer* n, scomplex* a, integer* lda, integer* ipiv, integer* info)
{
  return cgetf2_(m, n, a, lda, ipiv, info);
}
inline integer getf2(integer* m, integer* n, dcomplex* a, integer* lda, integer* ipiv, integer* info)
{
  return zgetf2_(m, n, a, lda, ipiv, info);
}

// --- QR factorization (classic) ---
inline integer geqrf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sgeqrf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer geqrf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dgeqrf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer geqrf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cgeqrf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer geqrf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zgeqrf_(m, n, a, lda, tau, work, lwork, info);
}
  
//--- computes the QR factorization of a general rectangular matrix using an unblocked algorithm ---  
inline integer geqr2(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* info)
{
  return sgeqr2_(m, n, a, lda, tau, work, info);
}
inline integer geqr2(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* info)
{
  return dgeqr2_(m, n, a, lda, tau, work, info);
}
inline integer geqr2(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* info)
{
  return cgeqr2_(m, n, a, lda, tau, work, info);
}
inline integer geqr2(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* info)
{
  return zgeqr2_(m, n, a, lda, tau, work, info);
}
  
//--- QR factorization with column pivoting of a real m-by-n matrix a ---  
inline integer geqpf(integer* m, integer* n, float* a, integer* lda, integer* jpvt, float* tau, float* work, integer* info)
{
  printf(" Function sgeqpf() has been deprecated. Please use sgep3() instead.\n");
  return sgeqpf_(m, n, a, lda, jpvt, tau, work, info);
}
inline integer geqpf(integer* m, integer* n, double* a, integer* lda, integer* jpvt, double* tau, double* work, integer* info)
{
  printf(" Function dgeqpf() has been deprecated. Please use dgep3() instead.\n");
  return dgeqpf_(m, n, a, lda, jpvt, tau, work, info);
}
inline integer geqpf(integer* m, integer* n, scomplex* a, integer* lda, integer* jpvt, scomplex* tau, scomplex* work, float* rwork, integer* info)
{
  printf(" Function cgeqpf() has been deprecated. Please use cgep3() instead.\n");
  return cgeqpf_(m, n, a, lda, jpvt, tau, work, rwork, info);
}
inline integer geqpf(integer* m, integer* n, dcomplex* a, integer* lda, integer* jpvt, dcomplex* tau, dcomplex* work, double* rwork, integer* info)
{
  printf(" Function zgeqpf() has been deprecated. Please use zgep3() instead.\n");
  return zgeqpf_(m, n, a, lda, jpvt, tau, work, rwork, info);
}

// --- QR factorization with column pivoting of a real m-by-n matrix ---
inline integer geqp3(integer* m, integer* n, float* a, integer* lda, integer* jpvt, float* tau, float* work, integer* lwork, integer* info)
{
  return sgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
}
inline integer geqp3(integer* m, integer* n, double* a, integer* lda, integer* jpvt, double* tau, double* work, integer* lwork, integer* info)
{
  return dgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
}
inline integer geqp3(integer* m, integer* n, scomplex* a, integer* lda, integer* jpvt, scomplex* tau, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
}
inline integer geqp3(integer* m, integer* n, dcomplex* a, integer* lda, integer* jpvt, dcomplex* tau, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
}

// --- LQ factorization (classic) ---
inline integer gelqf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sgelqf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer gelqf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dgelqf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer gelqf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cgelqf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer gelqf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zgelqf_(m, n, a, lda, tau, work, lwork, info);
}

// --- computes the LQ factorization of a general rectangular matrix using an unblocked algorithm ---
inline integer gelq2(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* info)
{
  return sgelq2_(m, n, a, lda, tau, work, info);
}
inline integer gelq2(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* info)
{
  return dgelq2_(m, n, a, lda, tau, work, info);
}
inline integer gelq2(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* info)
{
  return cgelq2_(m, n, a, lda, tau, work, info);
}
inline integer gelq2(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* info)
{
  return zgelq2_(m, n, a, lda, tau, work, info);
}

// --- LS solver ---
inline integer gelsd(integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, float*  s, float*  rcond, integer* rank, float* work, integer* lwork, integer* iwork, integer* info)
{
  return sgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
}
inline integer gelsd(integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, double* s, double* rcond, integer* rank, double* work, integer* lwork, integer* iwork, integer* info)
{
  return dgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
}
inline integer gelsd(integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, float*  s, float*  rcond, integer* rank, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* info)
{
  return cgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info);
}
inline integer gelsd(integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* s, double* rcond, integer* rank, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* info)
{
  return zgelsd_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info);
}

// ---  minimum-norm solution to a real linear least squares problem ---
inline integer gelss(integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, float*  s, float*  rcond, integer* rank, float* work, integer* lwork, integer* info)
{
  return sgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
}
inline integer gelss(integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, double* s, double* rcond, integer* rank, double* work, integer* lwork, integer* info)
{
  return dgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
}
inline integer gelss(integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, float*  s, float*  rcond, integer* rank, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info);
}
inline integer gelss(integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* s, double* rcond, integer* rank, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info);
}

// --- Triangular-transpose matrix multiply ---
inline integer lauum(char* uplo, integer* n, float* a, integer* lda, integer* info)
{
  return slauum_(uplo, n, a, lda, info);
}
inline integer lauum(char* uplo, integer* n, double* a, integer* lda, integer* info)
{
  return dlauum_(uplo, n, a, lda, info);
}
inline integer lauum(char* uplo, integer* n, scomplex* a, integer* lda, integer* info)
{
  return clauum_(uplo, n, a, lda, info);
}
inline integer lauum(char* uplo, integer* n, dcomplex* a, integer* lda, integer* info)
{
  return zlauum_(uplo, n, a, lda, info);
}

//--- computes the product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked algorithm) ---
inline integer lauu2(char* uplo, integer* n, float* a, integer* lda, integer* info)
{
  return slauu2_(uplo, n, a, lda, info);
}
inline integer lauu2(char* uplo, integer* n, double* a, integer* lda, integer* info)
{
  return dlauu2_(uplo, n, a, lda, info);
}
inline integer lauu2(char* uplo, integer* n, scomplex* a, integer* lda, integer* info)
{
  return clauu2_(uplo, n, a, lda, info);
}
inline integer lauu2(char* uplo, integer* n, dcomplex* a, integer* lda, integer* info)
{
  return zlauu2_(uplo, n, a, lda, info);
}
// --- Symmetric (hermitian) positive definite matrix inversion ---
inline integer potri(char* uplo, integer*  n, float* buff_A, integer*  ldim_A, integer* info)
{
  return spotri_(uplo, n, buff_A, ldim_A, info);
}
inline integer potri(char* uplo, integer*  n, double* buff_A, integer*  ldim_A, integer* info)
{
  return dpotri_(uplo, n, buff_A, ldim_A, info);
}
inline integer potri(char* uplo, integer*  n, scomplex* buff_A, integer*  ldim_A, integer* info)
{
  return cpotri_(uplo, n, buff_A, ldim_A, info);
}
inline integer potri(char* uplo, integer*  n, dcomplex* buff_A, integer*  ldim_A, integer* info)
{
  return zpotri_(uplo, n, buff_A, ldim_A, info);
}

// --- Triangular matrix inversion ---
inline integer trtri(char* uplo, char* diag, integer* n, float* a, integer* lda, integer* info)
{
  return strtri_(uplo, diag, n, a, lda, info);
}
inline integer trtri(char* uplo, char* diag, integer* n, double* a, integer* lda, integer* info)
{
  return dtrtri_(uplo, diag, n, a, lda, info);
}
inline integer trtri(char* uplo, char* diag, integer* n, scomplex* a, integer* lda, integer* info)
{
  return ctrtri_(uplo, diag, n, a, lda, info);
}
inline integer trtri(char* uplo, char* diag, integer* n, dcomplex* a, integer* lda, integer* info)
{
  return ztrtri_(uplo, diag, n, a, lda, info);
}

//--- computes the inverse of a triangular matrix (unblocked algorithm)---
inline integer trti2(char* uplo, char* diag, integer* n, float* a, integer* lda, integer* info)
{
  return strti2_(uplo, diag, n, a, lda, info);
}
inline integer trti2(char* uplo, char* diag, integer* n, double* a, integer* lda, integer* info)
{
  return dtrti2_(uplo, diag, n, a, lda, info);
}
inline integer trti2(char* uplo, char* diag, integer* n, scomplex* a, integer* lda, integer* info)
{
  return ctrti2_(uplo, diag, n, a, lda, info);
}
inline integer trti2(char* uplo, char* diag, integer* n, dcomplex* a, integer* lda, integer* info)
{
  return ztrti2_(uplo, diag, n, a, lda, info);
}

// --- Triangular Sylvester equation solve ---
inline integer trsyl(char* transa, char* transb, integer* isgn, integer* m, integer* n, float* a, integer* lda, float* b, integer* ldb, float* c, integer* ldc, float* scale, integer* info)
{
  return strsyl_(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
}
inline integer trsyl(char* transa, char* transb, integer* isgn, integer* m, integer* n, double* a, integer* lda, double* b, integer* ldb, double* c, integer* ldc, double* scale, integer* info)
{
  return dtrsyl_(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
}
inline integer trsyl(char* transa, char* transb, integer* isgn, integer* m, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* c, integer* ldc, float* scale, integer* info)
{
  return ctrsyl_(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
}
inline integer trsyl(char* transa, char* transb, integer* isgn, integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* c, integer* ldc, double* scale, integer* info)
{
  return ztrsyl_(transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info);
}

// --- Reduction to upper Hessenberg form ---
inline integer gehrd(integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sgehrd_(n, ilo, ihi, a, lda, tau, work, lwork, info);
}
inline integer gehrd(integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dgehrd_(n, ilo, ihi, a, lda, tau, work, lwork, info);
}
inline integer gehrd(integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cgehrd_(n, ilo, ihi, a, lda, tau, work, lwork, info);
}
inline integer gehrd(integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zgehrd_(n, ilo, ihi, a, lda, tau, work, lwork, info);
}

// --- reduces a general square matrix to upper Hessenberg form using an unblocked algorithm ---
inline integer gehd2(integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* tau, float* work, integer* info)
{
  return sgehd2_(n, ilo, ihi, a, lda, tau, work, info);
}
inline integer gehd2(integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* tau, double* work, integer* info)
{
  return dgehd2_(n, ilo, ihi, a, lda, tau, work, info);
}
inline integer gehd2(integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* info)
{
  return cgehd2_(n, ilo, ihi, a, lda, tau, work, info);
}
inline integer gehd2(integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* info)
{
  return zgehd2_(n, ilo, ihi, a, lda, tau, work, info);
}

// --- Reduction to tridiagonal form ---
inline integer sytrd(char* uplo, integer* n, float* a, integer* lda, float*  d, float*  e, float* tau, float* work, integer* lwork, integer* info)
{
  return ssytrd_(uplo, n, a, lda, d, e, tau, work, lwork, info);
}
inline integer sytrd(char* uplo, integer* n, double* a, integer* lda, double* d, double* e, double* tau, double* work, integer* lwork, integer* info)
{
  return dsytrd_(uplo, n, a, lda, d, e, tau, work, lwork, info);
}

//---reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarity transformation ---
inline integer sytd2(char* uplo, integer* n, float* a, integer* lda, float*  d, float*  e, float* tau, integer* info)
{
  return ssytd2_(uplo, n, a, lda, d, e, tau, info);
}
inline integer sytd2(char* uplo, integer* n, double* a, integer* lda, double* d, double* e, double* tau, integer* info)
{
  return dsytd2_(uplo, n, a, lda, d, e, tau, info);
}

//--- reduces a Hermitian matrix to real symmetric tridiagonal form by an unitary similarity transformation ---
inline integer hetd2(char* uplo, integer* n, scomplex* a, integer* lda, float*  d, float*  e, scomplex* tau, integer* info)
{
  return chetd2_(uplo, n, a, lda, d, e, tau, info);
}
inline integer hetd2(char* uplo, integer* n, dcomplex* a, integer* lda, double* d, double* e, dcomplex* tau, integer* info)
{
  return zhetd2_(uplo, n, a, lda, d, e, tau, info);
}

// --- Reduction to bidiagonal form ---
inline integer gebrd(integer* m, integer* n, float* a, integer* lda, float*  d, float*  e, float* tauq, float* taup, float* work, integer* lwork, integer* info)
{
  return sgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
}
inline integer gebrd(integer* m, integer* n, double* a, integer* lda, double* d, double* e, double* tauq, double* taup, double* work, integer* lwork, integer* info)
{
  return dgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
}
inline integer gebrd(integer* m, integer* n, scomplex* a, integer* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, integer* lwork, integer* info)
{
  return cgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
}
inline integer gebrd(integer* m, integer* n, dcomplex* a, integer* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, integer* lwork, integer* info)
{
  return zgebrd_(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
}

// --- reduces a general matrix to bidiagonal form using an unblocked algorithm ---
inline integer gebd2(integer* m, integer* n, float* a, integer* lda, float*  d, float*  e, float* tauq, float* taup, float* work, integer* info)
{
  return sgebd2_(m, n, a, lda, d, e, tauq, taup, work, info);
}
inline integer gebd2(integer* m, integer* n, double* a, integer* lda, double* d, double* e, double* tauq, double* taup, double* work, integer* info)
{
  return dgebd2_(m, n, a, lda, d, e, tauq, taup, work, info);
}
inline integer gebd2(integer* m, integer* n, scomplex* a, integer* lda, float*  d, float*  e, scomplex* tauq, scomplex* taup, scomplex* work, integer* info)
{
  return cgebd2_(m, n, a, lda, d, e, tauq, taup, work, info);
}
inline integer gebd2(integer* m, integer* n, dcomplex* a, integer* lda, double* d, double* e, dcomplex* tauq, dcomplex* taup, dcomplex* work, integer* info)
{
  return zgebd2_(m, n, a, lda, d, e, tauq, taup, work, info);
}

// --- Reduce Hermitian-definite generalized eigenproblem to standard form ---
inline integer sygst(integer* itype, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  return ssygst_(itype, uplo, n, a, lda, b, ldb, info);
}
inline integer sygst(integer* itype, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  return dsygst_(itype, uplo, n, a, lda, b, ldb, info);
}

// --- reduces a complex Hermitian-definite generalized eigenproblem to standard form ---
inline integer hegst(integer* itype, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  return chegst_(itype, uplo, n, a, lda, b, ldb, info);
}
inline integer hegst(integer* itype, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  return zhegst_(itype, uplo, n, a, lda, b, ldb, info);
}

// --- reduces a symmetric definite generalized eigenproblem to standard form ---
inline integer sygs2(integer* itype, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  return ssygs2_(itype, uplo, n, a, lda, b, ldb, info);
}
inline integer sygs2(integer* itype, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  return dsygs2_(itype, uplo, n, a, lda, b, ldb, info);
}

// --- reduces a Hermitian definite generalized eigenproblem to standard form ---
inline integer hegs2(integer* itype, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  return chegs2_(itype, uplo, n, a, lda, b, ldb, info);
}
inline integer hegs2(integer* itype, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  return zhegs2_(itype, uplo, n, a, lda, b, ldb, info);
}

// --- Accumulate block Householder matrix T (classic) ---
inline integer larft(char* direct, char* storev, integer* n, integer* k, float* v, integer* ldv, float* tau, float* t, integer* ldt)
{
  return slarft_(direct, storev, n, k, v, ldv, tau,  t, ldt);
}
inline integer larft(char* direct, char* storev, integer* n, integer* k, double* v, integer* ldv, double* tau, double* t, integer* ldt)
{
  return dlarft_(direct, storev, n, k, v, ldv, tau,  t, ldt);
}
inline integer larft(char* direct, char* storev, integer* n, integer* k, scomplex* v, integer* ldv, scomplex* tau, scomplex* t, integer* ldt)
{
  return clarft_(direct, storev, n, k, v, ldv, tau,  t, ldt);
}
inline integer larft(char* direct, char* storev, integer* n, integer* k, dcomplex* v, integer* ldv, dcomplex* tau, dcomplex* t, integer* ldt)
{
  return zlarft_(direct, storev, n, k, v, ldv, tau,  t, ldt);
}

// --- Generate a Householder vector (classic) ---
inline integer larfg( integer* n, float* alpha, float* x, integer* incx, float* tau)
{
  return slarfg_(n, alpha, x, incx, tau);
}
inline integer larfg( integer* n, double* alpha, double* x, integer* incx, double* tau)
{
  return dlarfg_(n, alpha, x, incx, tau);
}
inline integer larfg( integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* tau)
{
  return clarfg_(n, alpha, x, incx, tau);
}
inline integer larfg( integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* tau)
{
  return zlarfg_(n, alpha, x, incx, tau);
}

// --- generates an elementary reflector (Householder matrix) with non-negative beta ---
inline integer larfgp( integer* n, float* alpha, float* x, integer* incx, float* tau)
{
  return slarfgp_(n, alpha, x, incx, tau);
}
inline integer larfgp( integer* n, double* alpha, double* x, integer* incx, double* tau)
{
  return dlarfgp_(n, alpha, x, incx, tau);
}
inline integer larfgp( integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* tau)
{
  return clarfgp_(n, alpha, x, incx, tau);
}
inline integer larfgp( integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* tau)
{
  return zlarfgp_(n, alpha, x, incx, tau);
}

// --- Form Q from QR factorization ---
inline integer orgqr(integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sorgqr_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer orgqr(integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dorgqr_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- generates an M-by-N complex matrix Q with orthonormal columns ---
inline integer ungqr(integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cungqr_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer ungqr(integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zungqr_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- Apply Q or Q' from QR factorization ---
inline integer ormqr(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  return sormqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer ormqr(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  return dormqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

inline integer unmqr(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cunmqr_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}
inline integer unmqr(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zunmqr_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// ---multiplies a general matrix by the orthogonal matrix from a QR factorization determined by sgeqrf ---  
inline integer orm2r(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* info)
{
  return sorm2r_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline integer orm2r(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* info)
{
  return dorm2r_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

inline integer unm2r(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  return cunm2r_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, info);
}
inline integer unm2r(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  return zunm2r_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

// --- Form Q from LQ factorization ---
inline integer orglq(integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sorglq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer orglq(integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dorglq_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- generates an M-by-N complex matrix Q ---
inline integer unglq(integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cunglq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer unglq(integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zunglq_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- Apply Q or Q' from LQ factorization ---
inline integer ormlq(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  return sormlq_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer ormlq(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  return dormlq_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// --- overwrites the general complex M-by-N matrix C ---
inline integer unmlq(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cunmlq_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer unmlq(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zunmlq_(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// --- multiplies a general matrix by the orthogonal matrix from a LQ factorization determined by sgelqf ---
inline integer orml2(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* info)
{
  return sorml2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline integer orml2(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* info)
{
  return dorml2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

inline integer unml2(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  return cunml2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline integer unml2(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  return zunml2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

// --- Form Q from tridiagonal reduction ---
inline integer orgtr(char* uplo, integer* m, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sorgtr_(uplo, m, a, lda, tau, work, lwork, info);
}
inline integer orgtr(char* uplo, integer* m, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dorgtr_(uplo, m, a, lda, tau, work, lwork, info);
}

// --- generates a complex unitary matrix Q ---
inline integer ungtr(char* uplo, integer* m, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cungtr_(uplo, m, a, lda, tau, work, lwork, info);
}
inline integer ungtr(char* uplo, integer* m, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zungtr_(uplo, m, a, lda, tau, work, lwork, info);
}

// --- Apply Q or Q' from tridiagonal reduction ---
inline integer ormtr(char* side, char* uplo, char* trans, integer* m, integer* n, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  return sormtr_(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer ormtr(char* side, char* uplo, char* trans, integer* m, integer* n, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  return dormtr_(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
}

// --- UNMTR overwrites the general complex M-by-N matrix C ---
inline integer unmtr(char* side, char* uplo, char* trans, integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cunmtr_(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer unmtr(char* side, char* uplo, char* trans, integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zunmtr_(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);
}

// --- Form Q from bidiagonal reduction ---
inline integer orgbr(char* vect, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sorgbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
}
inline integer orgbr(char* vect, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dorgbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
}

// --- generates one of the complex unitary matrices Q ---  
inline integer ungbr(char* vect, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cungbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
}
inline integer ungbr(char* vect, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zungbr_(vect, m, n, k, a, lda, tau, work, lwork, info);
}

// --- Apply Q or Q' from bidiagonal reduction ---
inline integer ormbr(char* vect, char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  return sormbr_(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer ormbr(char* vect, char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  return dormbr_(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// --- overwrites the general complex M-by-N matrix C ---
inline integer unmbr(char* vect, char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cunmbr_(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer unmbr(char* vect, char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zunmbr_(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// --- Tridiagonal QR algorithm ---
inline integer steqr(char* compz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* info)
{
  return ssteqr_(compz, n, d, e, z, ldz, work, info);
}
inline integer steqr(char* compz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* info)
{
  return dsteqr_(compz, n, d, e, z, ldz, work, info);
}
inline integer steqr(char* compz, integer* n, float* d, float* e, scomplex* z, integer* ldz, float* work, integer* info)
{
  return csteqr_(compz, n, d, e, z, ldz, work, info);
}
inline integer steqr(char* compz, integer* n, double* d, double* e, dcomplex* z, integer* ldz, double* work, integer* info)
{
  return zsteqr_(compz, n, d, e, z, ldz, work, info);
}

// --- Tridiagonal divide-and-conquer algorithm ---
inline integer stedc(char* compz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return sstedc_(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer stedc(char* compz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dstedc_(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer stedc(char* compz, integer* n, float* d, float* e, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return cstedc_(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer stedc(char* compz, integer* n, double* d, double* e, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zstedc_(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- Tridiagonal MRRR algorithm ---
inline integer stemr(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, integer* m, float* w, float* z, integer* ldz, integer* nzc, integer* isuppz, logical* tryrac, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return sstemr_(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info);
}
inline integer stemr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, integer* m, double* w, double* z, integer* ldz, integer* nzc, integer* isuppz, integer* tryrac, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dstemr_(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info);
}
inline integer stemr(char* jobz, char* range, integer* n, float*  d, float*  e, float* vl, float* vu, integer* il, integer* iu, integer* m, float*  w, scomplex* z, integer* ldz, integer* nzc, integer* isuppz, integer* tryrac, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return cstemr_(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info);
}
inline integer stemr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, integer* m, double* w, dcomplex* z, integer* ldz, integer* nzc, integer* isuppz, integer* tryrac, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return zstemr_(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info);
}

// --- Hermitian eigenvalue decomposition (QR algorithm) ---
inline integer syev(char* jobz, char* uplo, integer* n, float* a, integer* lda, float*  w, float* work, integer* lwork, float* rwork, integer* info)
{
  return ssyev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}
inline integer syev(char* jobz, char* uplo, integer* n, double* a, integer* lda, double* w, double* work, integer* lwork, double* rwork, integer* info)
{
  return dsyev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

// --- performs the symmetric rank-1 update of a complex symmetric matrix ---
inline integer syr(char* uplo, integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* a, integer* lda)
{
  return csyr_(uplo, n, alpha, x, incx, a, lda);
}
inline integer syr(char* uplo, integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* a, integer* lda)
{
  return zsyr_(uplo, n, alpha, x, incx, a, lda);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline integer heev(char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, float* w, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}
inline integer heev(char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, double* w, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zheev_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline integer heevd(char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, float*  w, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return cheevd_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer heevd(char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, double* w, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zheevd_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline integer heevr(char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, float*  vl, float*  vu, integer* il, integer* iu, float*  abstol, integer* m, float*  w, scomplex* z, integer* ldz, integer* isuppz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return cheevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer heevr(char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, integer* isuppz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zheevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- Hermitian eigenvalue decomposition (divide-and-conquer) ---
inline integer syevd(char* jobz, char* uplo, integer* n, float* a, integer* lda, float*  w, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ssyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
}
inline integer syevd(char* jobz, char* uplo, integer* n, double* a, integer* lda, double* w, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dsyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
}

// --- Hermitian eigenvalue decomposition (MRRR) ---
inline integer syevr(char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ssyevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline integer syevr(char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dsyevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}

// --- Bidiagonal QR algorithm ---
inline integer bdsqr(char* uplo, integer* n, integer* ncvt, integer* nru, integer* ncc, float* d, float* e, float* vt, integer* ldvt, float* u, integer* ldu, float* c, integer* ldc, float* rwork, integer* info)
{
  return sbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
}
inline integer bdsqr(char* uplo, integer* n, integer* ncvt, integer* nru, integer* ncc, double* d, double* e, double* vt, integer* ldvt, double* u, integer* ldu, double* c, integer* ldc, double* rwork, integer* info)
{
  return dbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
}
inline integer bdsqr(char* uplo, integer* n, integer* ncvt, integer* nru, integer* ncc, float* d, float* e, scomplex* vt, integer* ldvt, scomplex* u, integer* ldu, scomplex* c, integer* ldc, float* rwork, integer* info)
{
  return cbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
}
inline integer bdsqr(char* uplo, integer* n, integer* ncvt, integer* nru, integer* ncc, double* d, double* e, dcomplex* vt, integer* ldvt, dcomplex* u, integer* ldu, dcomplex* c, integer* ldc, double* rwork, integer* info)
{
  return zbdsqr_(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
}

// --- Bidiagonal divide-and-conquor algorithm ---
inline integer bdsdc(char* uplo, char* compq, integer* n, float*  d, float*  e, float*  u, integer* ldu, float*  vt, integer* ldvt, float*  q, float*  iq, float* work, integer* iwork, integer* info)
{
  return sbdsdc_(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info);
}
inline integer bdsdc(char* uplo, char* compq, integer* n, double* d, double* e, double* u, integer* ldu, double* vt, integer* ldvt, double* q, double* iq, double* work, integer* iwork, integer* info)
{
  return dbdsdc_(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info);
}

// --- computes the singular value decomposition (SVD) of a real N-by-N (upper or lower) bidiagonal matrix ---
inline integer bdsvdx(char* uplo, char* jobz, char* range, integer* n, float* d, float* e, float *vl, float *vu, integer* il, integer* iu, integer* ns, float* s, float* z, integer* ldz, float* work, integer* iwork, integer* info)
{
  return sbdsvdx_(uplo, jobz, range, n, d, e, vl, vu, il, iu, ns, s, z, ldz, work, iwork, info);
}
inline integer bdsvdx(char* uplo, char* jobz, char* range, integer* n, double* d, double* e, double *vl, double *vu, integer* il, integer* iu, integer* ns, double* s, double* z, integer* ldz, double* work, integer* iwork, integer* info)
{
  return dbdsvdx_(uplo, jobz, range, n, d, e, vl, vu, il, iu, ns, s, z, ldz, work, iwork, info);
}

// --- computes the reciprocal condition numbers for the eigenvectors of a real symmetric or complex Hermitian matrix ---
inline integer disna(char* job, integer* m, integer* n,  float* d, float* sep, integer* info)
{
  return sdisna_(job, m, n,  d, sep, info);
}
inline integer disna(char* job, integer* m, integer* n,  double* d, double* sep, integer* info)
{
  return ddisna_(job, m, n,  d, sep, info);
}

// --- General matrix singular value decomposition (QR algorithm) ---
inline integer gesvd(char* jobu, char* jobv, integer* m, integer* n, float* a, integer* lda, float* s, float* u, integer* ldu, float* vt, integer* ldvt, float* work, integer* lwork, integer* info)
{
  return sgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}
inline integer gesvd(char* jobu, char* jobv, integer* m, integer* n, double* a, integer* lda, double* s, double* u, integer* ldu, double* vt, integer* ldvt, double* work, integer* lwork, integer* info)
{
  return dgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
}
inline integer gesvd(char* jobu, char* jobv, integer* m, integer* n, scomplex* a, integer* lda, float*  s, scomplex* u, integer* ldu, scomplex* vt, integer* ldvt, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
}
inline integer gesvd(char* jobu, char* jobv, integer* m, integer* n, dcomplex* a, integer* lda, double* s, dcomplex* u, integer* ldu, dcomplex* vt, integer* ldvt, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
}

// --- General matrix singular value decomposition (divide-and-conquer) ---
inline integer gesdd(char* jobz, integer* m, integer* n, float* a, integer* lda, float*  s, float* u, integer* ldu, float* vt, integer* ldvt, float* work, integer* lwork, integer* iwork, integer* info)
{
  return sgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
inline integer gesdd(char* jobz, integer* m, integer* n, double* a, integer* lda, double* s, double* u, integer* ldu, double* vt, integer* ldvt, double* work, integer* lwork, integer* iwork, integer* info)
{
  return dgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
inline integer gesdd(char* jobz, integer* m, integer* n, scomplex* a, integer* lda, float*  s, scomplex* u, integer* ldu, scomplex* vt, integer* ldvt, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* info)
{
  return cgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
}
inline integer gesdd(char* jobz, integer* m, integer* n, dcomplex* a, integer* lda, double* s, dcomplex* u, integer* ldu, dcomplex* vt, integer* ldvt, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* info)
{
  return zgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
}

// --- Swap rows ---
inline integer laswp(integer* n, float* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx)
{
  return slaswp_(n, a, lda, k1, k2, ipiv, incx);
}
inline integer laswp(integer* n, double* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx)
{
  return dlaswp_(n, a, lda, k1, k2, ipiv, incx);
}
inline integer laswp(integer* n, scomplex* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx)
{
  return claswp_(n, a, lda, k1, k2, ipiv, incx);
}
inline integer laswp(integer* n, dcomplex* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx)
{
  return zlaswp_(n, a, lda, k1, k2, ipiv, incx);
}

// --- Initialize a matrix ---
inline integer laset(char* uplo, integer* m, integer* n, float* alpha, float* beta, float* a, integer* lda)
{
  return slaset_(uplo, m, n, alpha, beta, a, lda);
}
inline integer laset(char* uplo, integer* m, integer* n, double* alpha, double* beta, double* a, integer* lda)
{
  return dlaset_(uplo, m, n, alpha, beta, a, lda);
}
inline integer laset(char* uplo, integer* m, integer* n, scomplex* alpha, scomplex* beta, scomplex* a, integer* lda)
{
  return claset_(uplo, m, n, alpha, beta, a, lda);
}
inline integer laset(char* uplo, integer* m, integer* n, dcomplex* alpha, dcomplex* beta, dcomplex* a, integer* lda)
{
  return zlaset_(uplo, m, n, alpha, beta, a, lda);
}

  // --- Bidiagonal block cs decomposition of orthogonal/unitary matrix  ---
inline integer bbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, integer* m, integer* p, integer* q, float* theta, float* phi, float* u1, integer* ldu1, float* u2, integer* ldu2, float* v1t, integer* ldv1t, float* v2t, integer* ldv2t, float* b11d, float* b11e, float* b12d, float* b12e, float* b21d, float* b21e, float* b22d, float* b22e, float* rwork, integer* lrwork, integer* info)
{
  return sbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, lrwork, info);
}
inline integer bbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, integer* m, integer* p, integer* q, double* theta, double* phi, double* u1, integer* ldu1, double* u2, integer* ldu2, double* v1t, integer* ldv1t, double* v2t, integer* ldv2t, double* b11d, double* b11e, double* b12d, double* b12e, double* b21d, double* b21e, double* b22d, double* b22e, double* rwork, integer* lrwork, integer* info)
{
  return dbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, lrwork, info);
}
inline integer bbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, integer* m, integer* p, integer* q, float* theta, float* phi, scomplex* u1, integer* ldu1, scomplex* u2, integer* ldu2, scomplex* v1t, integer* ldv1t, scomplex* v2t, integer* ldv2t, float* b11d, float* b11e, float* b12d, float* b12e, float* b21d, float* b21e, float* b22d, float* b22e, float* rwork, integer* lrwork, integer* info)
{
  return cbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, lrwork, info);
}
inline integer bbcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, integer* m, integer* p, integer* q, double* theta, double* phi, dcomplex* u1, integer* ldu1, dcomplex* u2, integer* ldu2, dcomplex* v1t, integer* ldv1t, dcomplex* v2t, integer* ldv2t, double* b11d, double* b11e, double* b12d, double* b12e, double* b21d, double* b21e, double* b22d, double* b22e, double* rwork, integer* lrwork, integer* info)
{
  return zbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, rwork, lrwork, info);
}
  // ---  Reduces a general band matrix to bidiagonal form.  ---
inline integer gbbrd(char* vect, integer* m, integer* n, integer* ncc, integer* kl, integer* ku, float* ab, integer* ldab, float* d, float* e, float* q, integer* ldq, float* pt, integer* ldpt, float* c, integer* ldc, float* work, integer* info)
{
  return sgbbrd_(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, info);
}
inline integer gbbrd(char* vect, integer* m, integer* n, integer* ncc, integer* kl, integer* ku, double* ab, integer* ldab, double* d, double* e, double* q, integer* ldq, double* pt, integer* ldpt, double* c, integer* ldc, double* work, integer* info)
{
  return dgbbrd_(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, info);
}
inline integer gbbrd(char* vect, integer* m, integer* n, integer* ncc, integer* kl, integer* ku, scomplex* ab, integer* ldab, float* d, float* e, scomplex* q, integer* ldq, scomplex* pt, integer* ldpt, scomplex* c, integer* ldc, scomplex* work, float* rwork, integer* info)
{
  return cgbbrd_(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, rwork, info);
}
inline integer gbbrd(char* vect, integer* m, integer* n, integer* ncc, integer* kl, integer* ku, dcomplex* ab, integer* ldab, double* d, double* e, dcomplex* q, integer* ldq, dcomplex* pt, integer* ldpt, dcomplex* c, integer* ldc, dcomplex* work, double* rwork, integer* info)
{
  return zgbbrd_(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, rwork, info);
}

// ---estimates the reciprocal of the condition number of a real general band matrix A, in either the 1-norm or the infinity-norm ---
inline integer gbcon(char* norm, integer* n, integer* kl, integer* ku, float* ab, integer* ldab, integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  return sgbcon_(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, iwork, info);
}
inline integer gbcon(char* norm, integer* n, integer* kl, integer* ku, double* ab, integer* ldab,  integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  return dgbcon_(norm, n, kl, ku, ab, ldab,  ipiv, anorm, rcond, work, iwork, info);
}
inline integer gbcon(char* norm, integer* n, integer* kl, integer* ku, scomplex* ab, integer* ldab,  integer* ipiv, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  return cgbcon_(norm, n, kl, ku, ab, ldab,  ipiv, anorm, rcond, work, rwork, info);
}
inline integer gbcon(char* norm, integer* n, integer* kl, integer* ku, dcomplex* ab, integer* ldab,  integer* ipiv, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  return zgbcon_(norm, n, kl, ku, ab, ldab,  ipiv, anorm, rcond, work, rwork, info);
}

// --- Computes row and column scalings intended to equilibrate  band matrix and reduce its condition number---
inline integer gbequ(integer* m, integer* n, integer* kl, integer* ku,  float* ab, integer* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  return sgbequ_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline integer gbequ(integer* m, integer* n, integer* kl, integer* ku,  double* ab, integer* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  return dgbequ_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline integer gbequ(integer* m, integer* n, integer* kl, integer* ku,  scomplex* ab, integer* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  return cgbequ_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline integer gbequ(integer* m, integer* n, integer* kl, integer* ku,  dcomplex* ab, integer* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  return zgbequ_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}

// --- Computes row and column scalings intended to equilibrate  band matrix and reduce its condition number(scaling factor:power of radix)---
inline integer gbequb(integer* m, integer* n, integer* kl, integer* ku,  float* ab, integer* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  return sgbequb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline integer gbequb(integer* m, integer* n, integer* kl, integer* ku,  double* ab, integer* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  return dgbequb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline integer gbequb(integer* m, integer* n, integer* kl, integer* ku,  scomplex* ab, integer* ldab, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  return cgbequb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}
inline integer gbequb(integer* m, integer* n, integer* kl, integer* ku,  dcomplex* ab, integer* ldab, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  return zgbequb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
}

// ---Improves computed solution to banded matrix and provides error bounds and backward estimates ---
inline integer gbrfs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, float* ab, integer* ldab, float* afb, integer* ldafb, integer* ipiv, float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sgbrfs_(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer gbrfs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  double* ab, integer* ldab, double* afb, integer* ldafb,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dgbrfs_(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer gbrfs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  scomplex* ab, integer* ldab, scomplex* afb, integer* ldafb,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cgbrfs_(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer gbrfs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  dcomplex* ab, integer* ldab, dcomplex* afb, integer* ldafb,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zgbrfs_(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// ---Improves computed solution to banded matrix and provides backward estimates,component wise and normwise error bounds---
inline integer gbrfsx(char* trans, char* equed, integer* n, integer* kl, integer* ku, integer* nrhs,  float* ab, integer* ldab,  float* afb, integer* ldafb,  integer* ipiv,  float* r,  float* c,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  return sgbrfsx_(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer gbrfsx(char* trans, char* equed, integer* n, integer* kl, integer* ku, integer* nrhs,  double* ab, integer* ldab,  double* afb, integer* ldafb,  integer* ipiv,  double* r,  double* c,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  return dgbrfsx_(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer gbrfsx(char* trans, char* equed, integer* n, integer* kl, integer* ku, integer* nrhs,  scomplex* ab, integer* ldab,  scomplex* afb, integer* ldafb,  integer* ipiv,  float* r,  float* c,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return cgbrfsx_(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer gbrfsx(char* trans, char* equed, integer* n, integer* kl, integer* ku, integer* nrhs,  dcomplex* ab, integer* ldab,  dcomplex* afb, integer* ldafb,  integer* ipiv,  double* r,  double* c,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zgbrfsx_(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb,  ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// ---Computes the solution to system of linear equations A * X = B for GB matrices ---
inline integer gbsv(integer* n, integer* kl, integer* ku, integer* nrhs, float* ab, integer* ldab, integer* ipiv, float* b, integer* ldb, integer* info)
{
  return sgbsv_(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
inline integer gbsv(integer* n, integer* kl, integer* ku, integer* nrhs, double* ab, integer* ldab, integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dgbsv_(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
inline integer gbsv(integer* n, integer* kl, integer* ku, integer* nrhs, scomplex* ab, integer* ldab, integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return cgbsv_(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
inline integer gbsv(integer* n, integer* kl, integer* ku, integer* nrhs, dcomplex* ab, integer* ldab, integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zgbsv_(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}

// ---It computes the solution to system of linear equations A * X = B for GB matrices along with error bounds ---
inline integer gbsvx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, float* ab, integer* ldab, float* afb, integer* ldafb, integer* ipiv, char* equed, float* r, float* c, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sgbsvx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer gbsvx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, double* ab, integer* ldab, double* afb, integer* ldafb, integer* ipiv, char* equed, double* r, double* c, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dgbsvx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer gbsvx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, scomplex* ab, integer* ldab, scomplex* afb, integer* ldafb, integer* ipiv, char* equed, float* r, float* c, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cgbsvx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer gbsvx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* afb, integer* ldafb, integer* ipiv, char* equed, double* r, double* c, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zgbsvx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// ---computes the solution to system of linear equations A * X = B for GB matrices with error, normalisations ---
inline integer gbsvxx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, float* ab, integer* ldab, float* afb, integer* ldafb, integer* ipiv, char* equed, float* r, float* c, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  return sgbsvxx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, 
  rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer gbsvxx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, double* ab, integer* ldab, double* afb, integer* ldafb, integer* ipiv, char* equed, double* r, double* c, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  return dgbsvxx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
  ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer gbsvxx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, scomplex* ab, integer* ldab, scomplex* afb, integer* ldafb, integer* ipiv, char* equed, float* r, float* c, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return cgbsvxx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer gbsvxx(char* fact, char* trans, integer* n, integer* kl, integer* ku, integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* afb, integer* ldafb, integer* ipiv, char* equed, double* r, double* c, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zgbsvxx_(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, 
  rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// ---Computes the LU factorization of a general band matrix using partial pivot---
inline integer gbtrf(integer* m, integer* n, integer* kl, integer* ku, float* ab, integer* ldab, integer* ipiv, integer* info)
{
  return sgbtrf_(m, n, kl, ku, ab, ldab, ipiv, info);
}
inline integer gbtrf(integer* m, integer* n, integer* kl, integer* ku, double* ab, integer* ldab, integer* ipiv, integer* info)
{
  return dgbtrf_(m, n, kl, ku, ab, ldab, ipiv, info);
}
inline integer gbtrf(integer* m, integer* n, integer* kl, integer* ku, scomplex* ab, integer* ldab, integer* ipiv, integer* info)
{
  return cgbtrf_(m, n, kl, ku, ab, ldab, ipiv, info);
}
inline integer gbtrf(integer* m, integer* n, integer* kl, integer* ku, dcomplex* ab, integer* ldab, integer* ipiv, integer* info)
{
  return zgbtrf_(m, n, kl, ku, ab, ldab, ipiv, info);
}

// ---solves a system of linear equationsA * X = B  or  A**T * X = B with a general band matrix A using the LU factorization computed by GBTRF ---
inline integer gbtrs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  float* ab, integer* ldab,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  return sgbtrs_(trans, n, kl, ku, nrhs,  ab, ldab,  ipiv, b, ldb, info);
}
inline integer gbtrs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  double* ab, integer* ldab,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dgbtrs_(trans, n, kl, ku, nrhs,  ab, ldab,  ipiv, b, ldb, info);
}
inline integer gbtrs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  scomplex* ab, integer* ldab,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return cgbtrs_(trans, n, kl, ku, nrhs,  ab, ldab,  ipiv, b, ldb, info);
}
inline integer gbtrs(char* trans, integer* n, integer* kl, integer* ku, integer* nrhs,  dcomplex* ab, integer* ldab,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zgbtrs_(trans, n, kl, ku, nrhs,  ab, ldab,  ipiv, b, ldb, info);
}

// ---Forms right left eigen vectors of real general matrix ---
inline integer gebak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  float* scale, integer* m, float* v, integer* ldv, integer* info)
{
  return sgebak_(job, side, n, ilo, ihi,  scale, m, v, ldv, info);
}
inline integer gebak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  double* scale, integer* m, double* v, integer* ldv, integer* info)
{
  return dgebak_(job, side, n, ilo, ihi,  scale, m, v, ldv, info);
}
inline integer gebak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  float* scale, integer* m, scomplex* v, integer* ldv, integer* info)
{
  return cgebak_(job, side, n, ilo, ihi,  scale, m, v, ldv, info);
}
inline integer gebak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  double* scale, integer* m, dcomplex* v, integer* ldv, integer* info)
{
  return zgebak_(job, side, n, ilo, ihi,  scale, m, v, ldv, info);
}

// --- Balances a real matrix ---
inline integer gebal(char* job, integer* n, float* a, integer* lda, integer* ilo, integer* ihi, float* scale, integer* info)
{
  return sgebal_(job, n, a, lda, ilo, ihi, scale, info);
}
inline integer gebal(char* job, integer* n, double* a, integer* lda, integer* ilo, integer* ihi, double* scale, integer* info)
{
  return dgebal_(job, n, a, lda, ilo, ihi, scale, info);
}
inline integer gebal(char* job, integer* n, scomplex* a, integer* lda, integer* ilo, integer* ihi, float* scale, integer* info)
{
  return cgebal_(job, n, a, lda, ilo, ihi, scale, info);
}
inline integer gebal(char* job, integer* n, dcomplex* a, integer* lda, integer* ilo, integer* ihi, double* scale, integer* info)
{
  return zgebal_(job, n, a, lda, ilo, ihi, scale, info);
}

// ---Computes row and column scaling to reduce condition number of matrix ---
inline integer geequb(integer* m, integer* n,  float* a, integer* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  return sgeequb_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline integer geequb(integer* m, integer* n,  double* a, integer* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  return dgeequb_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline integer geequb(integer* m, integer* n,  scomplex* a, integer* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  return cgeequb_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline integer geequb(integer* m, integer* n,  dcomplex* a, integer* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  return zgeequb_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}

// ---Computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices ---
inline integer gees(char* jobvs, char* sort, void* select, integer* n, float* a, integer* lda, integer* sdim, float* wr, float* wi, float* vs, integer* ldvs, float* work, integer* lwork, logical* bwork, integer* info)
{
  return sgees_(jobvs, sort, (L_fp)select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info);
}
inline integer gees(char* jobvs, char* sort, void* select, integer* n, double* a, integer* lda, integer* sdim, double* wr, double* wi, double* vs, integer* ldvs, double* work, integer* lwork, logical* bwork, integer* info)
{
  return dgees_(jobvs, sort, (L_fp)select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info);
}
inline integer gees(char* jobvs, char* sort, void* select, integer* n, scomplex* a, integer* lda, integer* sdim, scomplex* w, scomplex* vs, integer* ldvs, scomplex* work, integer* lwork, float* rwork, logical* bwork, integer* info)
{
  return cgees_(jobvs, sort, (L_fp)select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info);
}
inline integer gees(char* jobvs, char* sort, void* select, integer* n, dcomplex* a, integer* lda, integer* sdim, dcomplex* w, dcomplex* vs, integer* ldvs, dcomplex* work, integer* lwork, double* rwork, logical* bwork, integer* info)
{
  return zgees_(jobvs, sort, (L_fp)select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info);
}

// ---computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices with rconde and rcondv---
inline integer geesx(char* jobvs, char* sort, void* select, char* sense, integer* n, float* a, integer* lda, integer* sdim, float* wr, float* wi, float* vs, integer* ldvs, float* rconde, float* rcondv, float* work, integer* lwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  return sgeesx_(jobvs, sort, (L_fp)select, sense, n, a, lda, sdim, wr, wi, vs, ldvs, rconde, rcondv, work, lwork, iwork, liwork, bwork, info);
}
inline integer geesx(char* jobvs, char* sort, void* select, char* sense, integer* n, double* a, integer* lda, integer* sdim, double* wr, double* wi, double* vs, integer* ldvs, double* rconde, double* rcondv, double* work, integer* lwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  return dgeesx_(jobvs, sort, (L_fp)select, sense, n, a, lda, sdim, wr, wi, vs, ldvs, rconde, rcondv, work, lwork, iwork, liwork, bwork, info);
}
inline integer geesx(char* jobvs, char* sort, void* select, char* sense, integer* n, scomplex* a, integer* lda, integer* sdim, scomplex* w, scomplex* vs, integer* ldvs, float* rconde, float* rcondv, scomplex* work, integer* lwork, float* rwork, logical* bwork, integer* info)
{
  return cgeesx_(jobvs, sort, (L_fp)select, sense, n, a, lda, sdim, w, vs, ldvs, rconde, rcondv, work, lwork, rwork, bwork, info);
}
inline integer geesx(char* jobvs, char* sort, void* select, char* sense, integer* n, dcomplex* a, integer* lda, integer* sdim, dcomplex* w, dcomplex* vs, integer* ldvs, double* rconde, double* rcondv, dcomplex* work, integer* lwork, double* rwork, logical* bwork, integer* info)
{
  return zgeesx_(jobvs, sort, (L_fp)select, sense, n, a, lda, sdim, w, vs, ldvs, rconde, rcondv, work, lwork, rwork, bwork, info);
}

// ---Computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices ---
inline integer geev(char* jobvl, char* jobvr, integer* n, float* a, integer* lda, float* wr, float* wi, float* vl, integer* ldvl, float* vr, integer* ldvr, float* work, integer* lwork, integer* info)
{
  return sgeev_(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline integer geev(char* jobvl, char* jobvr, integer* n, double* a, integer* lda, double* wr, double* wi, double* vl, integer* ldvl, double* vr, integer* ldvr, double* work, integer* lwork, integer* info)
{
  return dgeev_(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline integer geev(char* jobvl, char* jobvr, integer* n, scomplex* a, integer* lda, scomplex* w, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cgeev_(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}
inline integer geev(char* jobvl, char* jobvr, integer* n, dcomplex* a, integer* lda, dcomplex* w, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zgeev_(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices(enabling conditions) ---
inline integer geevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, float* a, integer* lda, float* wr, float* wi, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* ilo, integer* ihi, float* scale, float* abnrm, float* rconde, float* rcondv, float* work, integer* lwork, integer* iwork, integer* info)
{
  return sgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
}
inline integer geevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, double* a, integer* lda, double* wr, double* wi, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* ilo, integer* ihi, double* scale, double* abnrm, double* rconde, double* rcondv, double* work, integer* lwork, integer* iwork, integer* info)
{
  return dgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
}
inline integer geevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, scomplex* a, integer* lda, scomplex* w, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* ilo, integer* ihi, float* scale, float* abnrm, float* rconde, float* rcondv, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
}
inline integer geevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, dcomplex* a, integer* lda, dcomplex* w, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* ilo, integer* ihi, double* scale, double* abnrm, double* rconde, double* rcondv, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zgeevx_(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
}

// ---Computes the singular value decomposition (SVD) of a real matrix ---
inline integer gejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt,
             char* jobp, integer* m, integer* n, float* a, integer* lda, float* sva,
             float* u, integer* ldu, float* v, integer* ldv,
             float* work, integer* lwork, integer* iwork, integer* info)
{
  return sgejsv_(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work, lwork, iwork, info);
}
inline integer gejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt,
             char* jobp, integer* m, integer* n, double* a, integer* lda, double* sva,
             double* u, integer* ldu, double* v, integer* ldv, double* work, integer* lwork, integer* iwork, integer* info)
{
  return dgejsv_(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv,  work, lwork, iwork, info);
}
inline integer gejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt, char* jobp, integer* m,
             integer* n, scomplex* a, integer* lda, float* sva, scomplex* u, integer* ldu,
             scomplex* v, integer* ldv, scomplex* cwork, integer* lwork, float* rwork, integer* lrwork, 
             integer* iwork, integer* info)
{
  return cgejsv_(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, cwork, lwork, rwork, lrwork, iwork, info);
}
inline integer gejsv(char* joba, char* jobu, char* jobv, char* jobr, char* jobt, char* jobp, integer* m, integer* n, dcomplex* a, integer* lda, double* sva, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* cwork, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* info)
{
  return zgejsv_(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, cwork, lwork, rwork, lrwork, iwork, info);
}

// --- Computes LQ factorization of a real matrix ---
inline integer gelq(integer* m, integer* n, float* a, integer* lda, float* t, integer* tsize, float* work, integer* lwork, integer* info)
{
  return sgelq_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline integer gelq(integer* m, integer* n, double* a, integer* lda, double* t, integer* tsize, double* work, integer* lwork, integer* info)
{
  return dgelq_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline integer gelq(integer* m, integer* n, scomplex* a, integer* lda, scomplex* t, integer* tsize, scomplex* work, integer* lwork, integer* info)
{
  return cgelq_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline integer gelq(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* t, integer* tsize, dcomplex* work, integer* lwork, integer* info)
{
  return zgelq_(m, n, a, lda, t, tsize, work, lwork, info);
}

// --- Solves overdetermined or underdetermined systems for GE matrices ---
inline integer gels(char* trans, integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  return sgels_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline integer gels(char* trans, integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  return dgels_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline integer gels(char* trans, integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return cgels_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline integer gels(char* trans, integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zgels_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}

// --- Solves overdetermined or underdetermined systems for GE matrices ---
inline integer gelsy(integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, integer* jpvt, float* rcond, integer* rank, float* work, integer* lwork, integer* info)
{
  return sgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
}
inline integer gelsy(integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, integer* jpvt, double* rcond, integer* rank, double* work, integer* lwork, integer* info)
{
  return dgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
}
inline integer gelsy(integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* jpvt, float* rcond, integer* rank, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info);
}
inline integer gelsy(integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* jpvt, double* rcond, integer* rank, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zgelsy_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info);
}

inline integer gelsx(integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, integer* jpvt, float* rcond, integer* rank, float* work, integer* info)
{
  printf(" Function sgelsx() has been deprecated. Please use sgelsy() instead.\n"); 
  return sgelsx_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, info);
}
inline integer gelsx(integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, integer* jpvt, double* rcond, integer* rank, double* work, integer* info)
{
  printf(" Function dgelsx() has been deprecated. Please use dgelsy() instead.\n"); 
  return dgelsx_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, info);
}
inline integer gelsx(integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* jpvt, float* rcond, integer* rank, scomplex* work, float* rwork, integer* info)
{
  printf(" Function cgelsx() has been deprecated. Please use cgelsy() instead.\n"); 
  return cgelsx_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, rwork, info);
}
inline integer gelsx(integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* jpvt, double* rcond, integer* rank, dcomplex* work, double* rwork, integer* info)
{
  printf(" Function zgelsx() has been deprecated. Please use zgelsy() instead.\n"); 
  return zgelsx_(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, rwork, info);
}

// --- Overwrites general matrix with a form compatible with orthogonal matrix ---
inline integer gemlq(char* side, char* trans, integer* m, integer* n, integer* k,  float* a, integer* lda,  float* t, integer* tsize, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  return sgemlq_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline integer gemlq(char* side, char* trans, integer* m, integer* n, integer* k,  double* a, integer* lda,  double* t, integer* tsize, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  return dgemlq_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline integer gemlq(char* side, char* trans, integer* m, integer* n, integer* k,  scomplex* a, integer* lda,  scomplex* t, integer* tsize, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cgemlq_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline integer gemlq(char* side, char* trans, integer* m, integer* n, integer* k,  dcomplex* a, integer* lda,  dcomplex* t, integer* tsize, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zgemlq_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}

// --- Multiples a matrix C by a real orthogonal or complex unitary matrix Q, as computed by ?geqr---
inline integer gemqr(char* side, char* trans, integer* m, integer* n, integer* k,  float* a, integer* lda,  float* t, integer* tsize, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  return sgemqr_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline integer gemqr(char* side, char* trans, integer* m, integer* n, integer* k,  double* a, integer* lda,  double* t, integer* tsize, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  return dgemqr_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline integer gemqr(char* side, char* trans, integer* m, integer* n, integer* k,  scomplex* a, integer* lda,  scomplex* t, integer* tsize, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cgemqr_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}
inline integer gemqr(char* side, char* trans, integer* m, integer* n, integer* k,  dcomplex* a, integer* lda,  dcomplex* t, integer* tsize, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zgemqr_(side, trans, m, n, k,  a, lda,  t, tsize, c, ldc, work, lwork, info);
}

// ---Multiplies a general matrix by the orthogonal/unitary matrix Q of the QR factorization formed by geqrt. ---
inline integer gemqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* nb,  float* v, integer* ldv,  float* t, integer* ldt, float* c, integer* ldc, float* work, integer* info)
{
  return sgemqrt_(side, trans, m, n, k, nb,  v, ldv,  t, ldt, c, ldc, work, info);
}
inline integer gemqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* nb,  double* v, integer* ldv,  double* t, integer* ldt, double* c, integer* ldc, double* work, integer* info)
{
  return dgemqrt_(side, trans, m, n, k, nb,  v, ldv,  t, ldt, c, ldc, work, info);
}
inline integer gemqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* nb,  scomplex* v, integer* ldv,  scomplex* t, integer* ldt, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  return cgemqrt_(side, trans, m, n, k, nb,  v, ldv,  t, ldt, c, ldc, work, info);
}
inline integer gemqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* nb,  dcomplex* v, integer* ldv,  dcomplex* t, integer* ldt, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  return zgemqrt_(side, trans, m, n, k, nb,  v, ldv,  t, ldt, c, ldc, work, info);
}

// --- computes a QL factorization of M-by-N matrix ---
inline integer geqlf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sgeqlf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer geqlf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dgeqlf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer geqlf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cgeqlf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer geqlf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zgeqlf_(m, n, a, lda, tau, work, lwork, info);
}

// --- computes a QR factorization of M-by-N matrix ---
inline integer geqr(integer* m, integer* n, float* a, integer* lda, float* t, integer* tsize, float* work, integer* lwork, integer* info)
{
  return sgeqr_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline integer geqr(integer* m, integer* n, double* a, integer* lda, double* t, integer* tsize, double* work, integer* lwork, integer* info)
{
  return dgeqr_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline integer geqr(integer* m, integer* n, scomplex* a, integer* lda, scomplex* t, integer* tsize, scomplex* work, integer* lwork, integer* info)
{
  return cgeqr_(m, n, a, lda, t, tsize, work, lwork, info);
}
inline integer geqr(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* t, integer* tsize, dcomplex* work, integer* lwork, integer* info)
{
  return zgeqr_(m, n, a, lda, t, tsize, work, lwork, info);
}

// --- computes a QR factorization of a M-by-N matrix ---
inline integer geqrfp(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sgeqrfp_(m, n, a, lda, tau, work, lwork, info);
}
inline integer geqrfp(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dgeqrfp_(m, n, a, lda, tau, work, lwork, info);
}
inline integer geqrfp(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cgeqrfp_(m, n, a, lda, tau, work, lwork, info);
}
inline integer geqrfp(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zgeqrfp_(m, n, a, lda, tau, work, lwork, info);
}

// --- computes a blocked QR factorization---
inline integer geqrt(integer* m, integer* n, integer* nb, float* a, integer* lda, float* t, integer* ldt, float* work, integer* info)
{
  return sgeqrt_(m, n, nb, a, lda, t, ldt, work, info);
}
inline integer geqrt(integer* m, integer* n, integer* nb, double* a, integer* lda, double* t, integer* ldt, double* work, integer* info)
{
  return dgeqrt_(m, n, nb, a, lda, t, ldt, work, info);
}
inline integer geqrt(integer* m, integer* n, integer* nb, scomplex* a, integer* lda, scomplex* t, integer* ldt, scomplex* work, integer* info)
{
  return cgeqrt_(m, n, nb, a, lda, t, ldt, work, info);
}
inline integer geqrt(integer* m, integer* n, integer* nb, dcomplex* a, integer* lda, dcomplex* t, integer* ldt, dcomplex* work, integer* info)
{
  return zgeqrt_(m, n, nb, a, lda, t, ldt, work, info);
}

// --- computes a blocked QR factorization ---
inline integer geqrt2(integer* m, integer* n, float* a, integer* lda, float* t, integer* ldt, integer* info)
{
  return sgeqrt2_(m, n, a, lda, t, ldt, info);
}
inline integer geqrt2(integer* m, integer* n, double* a, integer* lda, double* t, integer* ldt, integer* info)
{
  return dgeqrt2_(m, n, a, lda, t, ldt, info);
}
inline integer geqrt2(integer* m, integer* n, scomplex* a, integer* lda, scomplex* t, integer* ldt, integer* info)
{
  return cgeqrt2_(m, n, a, lda, t, ldt, info);
}
inline integer geqrt2(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* t, integer* ldt, integer* info)
{
  return zgeqrt2_(m, n, a, lda, t, ldt, info);
}

// --- computes a blocked QR factorization ---
inline integer geqrt3(integer* m, integer* n, float* a, integer* lda, float* t, integer* ldt, integer* info)
{
  return sgeqrt3_(m, n, a, lda, t, ldt, info);
}
inline integer geqrt3(integer* m, integer* n, double* a, integer* lda, double* t, integer* ldt, integer* info)
{
  return dgeqrt3_(m, n, a, lda, t, ldt, info);
}
inline integer geqrt3(integer* m, integer* n, scomplex* a, integer* lda, scomplex* t, integer* ldt, integer* info)
{
  return cgeqrt3_(m, n, a, lda, t, ldt, info);
}
inline integer geqrt3(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* t, integer* ldt, integer* info)
{
  return zgeqrt3_(m, n, a, lda, t, ldt, info);
}

// --- Improves the computed solution to a system of linear equations---
inline integer gerfs(char* trans, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sgerfs_(trans, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer gerfs(char* trans, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dgerfs_(trans, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer gerfs(char* trans, integer* n, integer* nrhs, scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cgerfs_(trans, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer gerfs(char* trans, integer* n, integer* nrhs, dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zgerfs_(trans, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- GERFSX improves the computed solution to a system of linear equations ---
inline integer gerfsx(char* trans, char* equed, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  integer* ipiv,  float* r,  float* c,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  return sgerfsx_(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer gerfsx(char* trans, char* equed, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  integer* ipiv,  double* r,  double* c,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  return dgerfsx_(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer gerfsx(char* trans, char* equed, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  float* r,  float* c,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return cgerfsx_(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer gerfsx(char* trans, char* equed, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  double* r,  double* c,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zgerfsx_(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes a RQ factorization of a M-by-N matrix ---
inline integer gerqf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sgerqf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer gerqf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dgerqf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer gerqf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cgerqf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer gerqf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zgerqf_(m, n, a, lda, tau, work, lwork, info);
}

// --- computes the solution to a real system of linear equations ---
inline integer gesv(integer* n, integer* nrhs, float* a, integer* lda, integer* ipiv, float* b, integer* ldb, integer* info)
{
  return sgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
inline integer gesv(integer* n, integer* nrhs, double* a, integer* lda, integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
inline integer gesv(integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return cgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
inline integer gesv(integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}

// --- computes the singular value decomposition (SVD) ---
inline integer gesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, integer* m, integer* n, float* a, integer* lda, float* s, float* u, integer* ldu, float* v, integer* ldv, integer* numrank, integer* iwork, integer* liwork, float* work, integer* lwork, float* rwork, integer* lrwork, integer* info)
{
  return sgesvdq_(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, 
            iwork, liwork, work, lwork, rwork, lrwork, info);
}
inline integer gesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, integer* m, integer* n, double* a, integer* lda, double* s, double* u, integer* ldu, double* v, integer* ldv, integer* numrank, integer* iwork, integer* liwork, double* work, integer* lwork, double* rwork, integer* lrwork, integer* info)
{
  return dgesvdq_(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, iwork, liwork, work, lwork, rwork, lrwork, info);
}
inline integer gesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, integer* m, 
            integer* n, scomplex* a, integer* lda, float* s, scomplex* u, integer* ldu,
            scomplex* v, integer* ldv, integer* numrank, integer* iwork, integer* liwork,
            scomplex* cwork, integer* lcwork, float* rwork, integer* lrwork, integer* info)
{
  return cgesvdq_(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, iwork, liwork, cwork, lcwork, rwork, lrwork, info);
}
inline integer gesvdq(char* joba, char* jobp, char* jobr, char* jobu, char* jobv, integer* m, integer* n, dcomplex* a, integer* lda, double* s, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, integer* numrank, integer* iwork, integer* liwork, dcomplex* cwork, integer* lcwork, double* rwork, integer* lrwork, integer* info)
{
  return zgesvdq_(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, iwork, liwork, cwork, lcwork, rwork, lrwork, info);
}

// --- computes the singular value decomposition (SVD) for GE matrices ---
inline integer gesvdx(char* jobu, char* jobvt, char* range, integer* m, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, integer* ns, float* s, float* u, integer* ldu, float* vt, integer* ldvt, float* work, integer* lwork, integer* iwork, integer* info)
{
  return sgesvdx_(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
inline integer gesvdx(char* jobu, char* jobvt, char* range, integer* m, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, integer* ns, double* s, double* u, integer* ldu, double* vt, integer* ldvt, double* work, integer* lwork, integer* iwork, integer* info)
{
  return dgesvdx_(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
}
inline integer gesvdx(char* jobu, char* jobvt, char* range, integer* m, integer* n, scomplex* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, integer* ns, float* s, scomplex* u, integer* ldu, scomplex* vt, integer* ldvt, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* info)
{
  return cgesvdx_(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
}
inline integer gesvdx(char* jobu, char* jobvt, char* range, integer* m, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, integer* ns, double* s, dcomplex* u, integer* ldu, dcomplex* vt, integer* ldvt, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* info)
{
  return zgesvdx_(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
}

// --- computes the singular value decomposition (SVD) of a M-by-N matrix ---
inline integer gesvj(char* joba, char* jobu, char* jobv, integer* m, integer* n, float* a, integer* lda, float* sva, integer* mv, float* v, integer* ldv, float* work, integer* lwork, integer* info)
{
  return sgesvj_(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info);
}
inline integer gesvj(char* joba, char* jobu, char* jobv, integer* m, integer* n, double* a, integer* lda, double* sva, integer* mv, double* v, integer* ldv, double* work, integer* lwork, integer* info)
{
  return dgesvj_(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info);
}
inline integer gesvj(char* joba, char* jobu, char* jobv, integer* m, integer* n, scomplex* a, integer* lda, float* sva, integer* mv, scomplex* v, integer* ldv, scomplex* cwork, integer* lwork, float* rwork, integer* lrwork, integer* info)
{
  return cgesvj_(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, cwork, lwork, rwork, lrwork, info);
}
inline integer gesvj(char* joba, char* jobu, char* jobv, integer* m, integer* n, dcomplex* a, integer* lda, double* sva, integer* mv, dcomplex* v, integer* ldv, dcomplex* cwork, integer* lwork, double* rwork, integer* lrwork, integer* info)
{
  return zgesvj_(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, cwork, lwork, rwork, lrwork, info);
}

// --- computes the solution to system of linear equations ---
inline integer gesvx(char* fact, char* trans, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, integer* ipiv, char* equed, float* r, float* c, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sgesvx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer gesvx(char* fact, char* trans, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, integer* ipiv, char* equed, double* r, double* c, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dgesvx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer gesvx(char* fact, char* trans, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv, char* equed, float* r, float* c, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cgesvx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer gesvx(char* fact, char* trans, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv, char* equed, double* r, double* c, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zgesvx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equation ---
inline integer gesvxx(char* fact, char* trans, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, integer* ipiv, char* equed, float* r, float* c, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  return sgesvxx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer gesvxx(char* fact, char* trans, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, integer* ipiv, char* equed, double* r, double* c, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  return dgesvxx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer gesvxx(char* fact, char* trans, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv, char* equed, float* r, float* c, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return cgesvxx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer gesvxx(char* fact, char* trans, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv, char* equed, double* r, double* c, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zgesvxx_(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes an LU factorization ---
inline integer getrf2(integer* m, integer* n, float* a, integer* lda, integer* ipiv, integer* info)
{
  return sgetrf2_(m, n, a, lda, ipiv, info);
}
inline integer getrf2(integer* m, integer* n, double* a, integer* lda, integer* ipiv, integer* info)
{
  return dgetrf2_(m, n, a, lda, ipiv, info);
}
inline integer getrf2(integer* m, integer* n, scomplex* a, integer* lda, integer* ipiv, integer* info)
{
  return cgetrf2_(m, n, a, lda, ipiv, info);
}
inline integer getrf2(integer* m, integer* n, dcomplex* a, integer* lda, integer* ipiv, integer* info)
{
  return zgetrf2_(m, n, a, lda, ipiv, info);
}

// --- computes the inverse of a matrix using the LU factorization ---
inline integer getri(integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* lwork, integer* info)
{
  return sgetri_(n, a, lda, ipiv, work, lwork, info);
}
inline integer getri(integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* lwork, integer* info)
{
  return dgetri_(n, a, lda, ipiv, work, lwork, info);
}
inline integer getri(integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return cgetri_(n, a, lda, ipiv, work, lwork, info);
}
inline integer getri(integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zgetri_(n, a, lda, ipiv, work, lwork, info);
}

// --- solves a system of linear equations ---
inline integer getrs(char* trans, integer* n, integer* nrhs,  float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  return sgetrs_(trans, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer getrs(char* trans, integer* n, integer* nrhs,  double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dgetrs_(trans, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer getrs(char* trans, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return cgetrs_(trans, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer getrs(char* trans, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zgetrs_(trans, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves overdetermined or underdetermined real linear systems ---
inline integer getsls(char* trans, integer* m, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  return sgetsls_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline integer getsls(char* trans, integer* m, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  return dgetsls_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline integer getsls(char* trans, integer* m, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return cgetsls_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}
inline integer getsls(char* trans, integer* m, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zgetsls_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}

// --- forms the right or left eigenvectors of a real generalized eigenvalue problem ---
inline integer ggbak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  float* lscale,  float* rscale, integer* m, float* v, integer* ldv, integer* info)
{
  return sggbak_(job, side, n, ilo, ihi,  lscale,  rscale, m, v, ldv, info);
}
inline integer ggbak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  double* lscale,  double* rscale, integer* m, double* v, integer* ldv, integer* info)
{
  return dggbak_(job, side, n, ilo, ihi,  lscale,  rscale, m, v, ldv, info);
}
inline integer ggbak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  float* lscale,  float* rscale, integer* m, scomplex* v, integer* ldv, integer* info)
{
  return cggbak_(job, side, n, ilo, ihi,  lscale,  rscale, m, v, ldv, info);
}
inline integer ggbak(char* job, char* side, integer* n, integer* ilo, integer* ihi,  double* lscale,  double* rscale, integer* m, dcomplex* v, integer* ldv, integer* info)
{
  return zggbak_(job, side, n, ilo, ihi,  lscale,  rscale, m, v, ldv, info);
}

// --- balances a pair of general real matrices (A,B) ---
inline integer ggbal(char* job, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* ilo, integer* ihi, float* lscale, float* rscale, float* work, integer* info)
{
  return sggbal_(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
}
inline integer ggbal(char* job, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* ilo, integer* ihi, double* lscale, double* rscale, double* work, integer* info)
{
  return dggbal_(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
}
inline integer ggbal(char* job, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* ilo, integer* ihi, float* lscale, float* rscale, float* work, integer* info)
{
  return cggbal_(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
}
inline integer ggbal(char* job, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* ilo, integer* ihi, double* lscale, double* rscale, double* work, integer* info)
{
  return zggbal_(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
}

// --- GGES computes the eigenvalues---
inline integer gges(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* sdim, float* alphar, float* alphai, float* beta, float* vsl, integer* ldvsl, float* vsr, integer* ldvsr, float* work, integer* lwork, logical* bwork, integer* info)
{
  return sgges_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info);
}
inline integer gges(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* sdim, double* alphar, double* alphai, double* beta, double* vsl, integer* ldvsl, double* vsr, integer* ldvsr, double* work, integer* lwork, logical* bwork, integer* info)
{
  return dgges_(jobvsl, jobvsr, sort, (L_fp) selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info);
}
inline integer gges(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* sdim, scomplex* alpha, scomplex* beta, scomplex* vsl, integer* ldvsl, scomplex* vsr, integer* ldvsr, scomplex* work, integer* lwork, float* rwork, logical* bwork, integer* info)
{
  return cgges_(jobvsl, jobvsr, sort, (L_fp) selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info);
}
inline integer gges(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* sdim, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, integer* ldvsl, dcomplex* vsr, integer* ldvsr, dcomplex* work, integer* lwork, double* rwork, logical* bwork, integer* info)
{
  return zgges_(jobvsl, jobvsr, sort, (L_fp) selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info);
}

inline integer gegs(char* jobvsl, char* jobvsr, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vsl, integer* ldvsl, float* vsr, integer* ldvsr, float* work, integer* lwork, integer* info)
{
  printf(" Function sgegs() has been deprecated. Please use sgges() instead.\n");  
  return sgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, info);
}
inline integer gegs(char* jobvsl, char* jobvsr, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vsl, integer* ldvsl, double* vsr, integer* ldvsr, double* work, integer* lwork, integer* info)
{
  printf(" Function dgegs() has been deprecated. Please use dgges() instead.\n");  
  return dgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, info);
}
inline integer gegs(char* jobvsl, char* jobvsr, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vsl, integer* ldvsl, scomplex* vsr, integer* ldvsr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  printf(" Function cgegs() has been deprecated. Please use cgges() instead.\n");  
  return cgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, info);
}
inline integer gegs(char* jobvsl, char* jobvsr, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, integer* ldvsl, dcomplex* vsr, integer* ldvsr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  printf(" Function zgegs() has been deprecated. Please use zgges() instead.\n"); 
  return zgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, info);
}

// --- SGGES3 computes the eigenvalues---
inline integer gges3(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* sdim, float* alphar, float* alphai, float* beta, float* vsl, integer* ldvsl, float* vsr, integer* ldvsr, float* work, integer* lwork, logical* bwork, integer* info)
{
  return sgges3_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info);
}
inline integer gges3(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* sdim, double* alphar, double* alphai, double* beta, double* vsl, integer* ldvsl, double* vsr, integer* ldvsr, double* work, integer* lwork, logical* bwork, integer* info)
{
  return dgges3_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info);
}
inline integer gges3(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* sdim, scomplex* alpha, scomplex* beta, scomplex* vsl, integer* ldvsl, scomplex* vsr, integer* ldvsr, scomplex* work, integer* lwork, float* rwork, logical* bwork, integer* info)
{
  return cgges3_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info);
}
inline integer gges3(char* jobvsl, char* jobvsr, char* sort, void *selctg, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* sdim, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, integer* ldvsl, dcomplex* vsr, integer* ldvsr, dcomplex* work, integer* lwork, double* rwork, logical* bwork, integer* info)
{
  return zgges3_(jobvsl, jobvsr, sort, (L_fp)selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info);
}

// --- computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices ---
inline integer ggesx(char* jobvsl, char* jobvsr, char* sort, void *selctg, char* sense, integer* n, float* a, integer* lda, float* b, integer* ldb, integer* sdim, float* alphar, float* alphai, float* beta, float* vsl, integer* ldvsl, float* vsr, integer* ldvsr, float* rconde, float* rcondv, float* work, integer* lwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  return sggesx_(jobvsl, jobvsr, sort, (L_fp)selctg, sense, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, bwork, info);
}
inline integer ggesx(char* jobvsl, char* jobvsr, char* sort, void *selctg, char* sense, integer* n, double* a, integer* lda, double* b, integer* ldb, integer* sdim, double* alphar, double* alphai, double* beta, double* vsl, integer* ldvsl, double* vsr, integer* ldvsr, double* rconde, double* rcondv, double* work, integer* lwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  return dggesx_(jobvsl, jobvsr, sort, (L_fp)selctg, sense, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, bwork, info);
}
inline integer ggesx(char* jobvsl, char* jobvsr, char* sort, void *selctg, char* sense, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* sdim, scomplex* alpha, scomplex* beta, scomplex* vsl, integer* ldvsl, scomplex* vsr, integer* ldvsr, float* rconde, float* rcondv, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  return cggesx_(jobvsl, jobvsr, sort, (L_fp)selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork, bwork, info);
}
inline integer ggesx(char* jobvsl, char* jobvsr, char* sort, void *selctg, char* sense, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* sdim, dcomplex* alpha, dcomplex* beta, dcomplex* vsl, integer* ldvsl, dcomplex* vsr, integer* ldvsr, double* rconde, double* rcondv, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* liwork, logical* bwork, integer* info)
{
  return zggesx_(jobvsl, jobvsr, sort, (L_fp)selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork, bwork, info);
}

// --- computes the eigenvalues, the left and/or right eigenvectors for GE matrices---
inline integer ggev(char* jobvl, char* jobvr, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vl, integer* ldvl, float* vr, integer* ldvr, float* work, integer* lwork, integer* info)
{
  return sggev_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline integer ggev(char* jobvl, char* jobvr, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vl, integer* ldvl, double* vr, integer* ldvr, double* work, integer* lwork, integer* info)
{
  return dggev_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline integer ggev(char* jobvl, char* jobvr, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}
inline integer ggev(char* jobvl, char* jobvr, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

inline integer gegv(char* jobvl, char* jobvr, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vl, integer* ldvl, float* vr, integer* ldvr, float* work, integer* lwork, integer* info)
{
  printf(" Function sgegv() has been deprecated. Please use sggev() instead.\n"); 
  return sgegv_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline integer gegv(char* jobvl, char* jobvr, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vl, integer* ldvl, double* vr, integer* ldvr, double* work, integer* lwork, integer* info)
{
  printf(" Function dgegv() has been deprecated. Please use dggev() instead.\n");  
  return dgegv_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline integer gegv(char* jobvl, char* jobvr, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  printf(" Function cgegv() has been deprecated. Please use cggev() instead.\n");  
  return cgegv_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}
inline integer gegv(char* jobvl, char* jobvr, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  printf(" Function zgegv() has been deprecated. Please use zggev() instead.\n");  
  return zgegv_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

// --- computes the eigenvalues, the left and/or right eigenvectors for GE matrices (blocked algorithm)---
inline integer ggev3(char* jobvl, char* jobvr, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vl, integer* ldvl, float* vr, integer* ldvr, float* work, integer* lwork, integer* info)
{
  return sggev3_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline integer ggev3(char* jobvl, char* jobvr, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vl, integer* ldvl, double* vr, integer* ldvr, double* work, integer* lwork, integer* info)
{
  return dggev3_(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline integer ggev3(char* jobvl, char* jobvr, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cggev3_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}
inline integer ggev3(char* jobvl, char* jobvr, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zggev3_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

// --- SGGEVX computes the eigenvalues, the left and/or right eigenvectors for GE matrices ---
inline integer ggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* ilo, integer* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* rconde, float* rcondv, float* work, integer* lwork, integer* iwork, logical* bwork, integer* info)
{
  return sggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info);
}
inline integer ggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* ilo, integer* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* rconde, double* rcondv, double* work, integer* lwork, integer* iwork, logical* bwork, integer* info)
{
  return dggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info);
}
inline integer ggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* ilo, integer* ihi, float* lscale, float* rscale, float* abnrm, float* bbnrm, float* rconde, float* rcondv, scomplex* work, integer* lwork, float* rwork, integer* iwork, logical* bwork, integer* info)
{
  return cggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork, info);
}
inline integer ggevx(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* ilo, integer* ihi, double* lscale, double* rscale, double* abnrm, double* bbnrm, double* rconde, double* rcondv, dcomplex* work, integer* lwork, double* rwork, integer* iwork, logical* bwork, integer* info)
{
  return zggevx_(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, rwork, iwork, bwork, info);
}

// --- solves a general Gauss-Markov linear model (GLM) problem ---
inline integer ggglm(integer* n, integer* m, integer* p, float* a, integer* lda, float* b, integer* ldb, float* d, float* x, float* y, float* work, integer* lwork, integer* info)
{
  return sggglm_(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
}
inline integer ggglm(integer* n, integer* m, integer* p, double* a, integer* lda, double* b, integer* ldb, double* d, double* x, double* y, double* work, integer* lwork, integer* info)
{
  return dggglm_(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
}
inline integer ggglm(integer* n, integer* m, integer* p, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* d, scomplex* x, scomplex* y, scomplex* work, integer* lwork, integer* info)
{
  return cggglm_(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
}
inline integer ggglm(integer* n, integer* m, integer* p, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* d, dcomplex* x, dcomplex* y, dcomplex* work, integer* lwork, integer* info)
{
  return zggglm_(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
}

// --- reduces a pair of real matrices (A,B) ---
inline integer gghd3(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* b, integer* ldb, float* q, integer* ldq, float* z, integer* ldz, float* work, integer* lwork, integer* info)
{
  return sgghd3_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info);
}
inline integer gghd3(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* b, integer* ldb, double* q, integer* ldq, double* z, integer* ldz, double* work, integer* lwork, integer* info)
{
  return dgghd3_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info);
}
inline integer gghd3(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* q, integer* ldq, scomplex* z, integer* ldz, scomplex* work, integer* lwork, integer* info)
{
  return cgghd3_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info);
}
inline integer gghd3(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, integer* info)
{
  return zgghd3_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info);
}

// --- reduces a pair of real matrices (A,B) ---
inline integer gghrd(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* b, integer* ldb, float* q, integer* ldq, float* z, integer* ldz, integer* info)
{
  return sgghrd_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
}
inline integer gghrd(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* b, integer* ldb, double* q, integer* ldq, double* z, integer* ldz, integer* info)
{
  return dgghrd_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
}
inline integer gghrd(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* q, integer* ldq, scomplex* z, integer* ldz, integer* info)
{
  return cgghrd_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
}
inline integer gghrd(char* compq, char* compz, integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, integer* info)
{
  return zgghrd_(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
}

// --- solves overdetermined or underdetermined systems for OTHER matrices ---
inline integer gglse(integer* m, integer* n, integer* p, float* a, integer* lda, float* b, integer* ldb, float* c, float* d, float* x, float* work, integer* lwork, integer* info)
{
  return sgglse_(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info);
}
inline integer gglse(integer* m, integer* n, integer* p, double* a, integer* lda, double* b, integer* ldb, double* c, double* d, double* x, double* work, integer* lwork, integer* info)
{
  return dgglse_(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info);
}
inline integer gglse(integer* m, integer* n, integer* p, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* c, scomplex* d, scomplex* x, scomplex* work, integer* lwork, integer* info)
{
  return cgglse_(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info);
}
inline integer gglse(integer* m, integer* n, integer* p, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* c, dcomplex* d, dcomplex* x, dcomplex* work, integer* lwork, integer* info)
{
  return zgglse_(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info);
}

// --- computes a generalized QR factorization ---
inline integer ggqrf(integer* n, integer* m, integer* p, float* a, integer* lda, float* taua, float* b, integer* ldb, float* taub, float* work, integer* lwork, integer* info)
{
  return sggqrf_(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline integer ggqrf(integer* n, integer* m, integer* p, double* a, integer* lda, double* taua, double* b, integer* ldb, double* taub, double* work, integer* lwork, integer* info)
{
  return dggqrf_(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline integer ggqrf(integer* n, integer* m, integer* p, scomplex* a, integer* lda, scomplex* taua, scomplex* b, integer* ldb, scomplex* taub, scomplex* work, integer* lwork, integer* info)
{
  return cggqrf_(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline integer ggqrf(integer* n, integer* m, integer* p, dcomplex* a, integer* lda, dcomplex* taua, dcomplex* b, integer* ldb, dcomplex* taub, dcomplex* work, integer* lwork, integer* info)
{
  return zggqrf_(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info);
}

// --- computes a generalized RQ factorization ---
inline integer ggrqf(integer* m, integer* p, integer* n, float* a, integer* lda, float* taua, float* b, integer* ldb, float* taub, float* work, integer* lwork, integer* info)
{
  return sggrqf_(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline integer ggrqf(integer* m, integer* p, integer* n, double* a, integer* lda, double* taua, double* b, integer* ldb, double* taub, double* work, integer* lwork, integer* info)
{
  return dggrqf_(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline integer ggrqf(integer* m, integer* p, integer* n, scomplex* a, integer* lda, scomplex* taua, scomplex* b, integer* ldb, scomplex* taub, scomplex* work, integer* lwork, integer* info)
{
  return cggrqf_(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info);
}
inline integer ggrqf(integer* m, integer* p, integer* n, dcomplex* a, integer* lda, dcomplex* taua, dcomplex* b, integer* ldb, dcomplex* taub, dcomplex* work, integer* lwork, integer* info)
{
  return zggrqf_(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info);
}

// --- computes the singular value decomposition (SVD) for OTHER matrices ---
inline integer ggsvd3(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, float* a, integer* lda, float* b, integer* ldb, float* alpha, float* beta, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, float* work, integer* lwork, integer* iwork, integer* info)
{
  return sggsvd3_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info);
}
inline integer ggsvd3(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, double* a, integer* lda, double* b, integer* ldb, double* alpha, double* beta, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, double* work, integer* lwork, integer* iwork, integer* info)
{
  return dggsvd3_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info);
}
inline integer ggsvd3(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* alpha, float* beta, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* info)
{
  return cggsvd3_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, rwork, iwork, info);
}
inline integer ggsvd3(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* alpha, double* beta, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* info)
{
  return zggsvd3_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, rwork, iwork, info);
}

// --- computes the singular value decomposition ---
inline integer ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, float* a, integer* lda, float* b, integer* ldb, float* alpha, float* beta, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, float* work, integer* iwork, integer* info)
{
  printf(" Function sggsvd() has been deprecated. Please use sggsvd3() instead.\n"); 
  return sggsvd_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, iwork, info);
}
inline integer ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, double* a, integer* lda, double* b, integer* ldb, double* alpha, double* beta, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, double* work, integer* iwork, integer* info)
{
  printf(" Function dggsvd() has been deprecated. Please use dggsvd3() instead.\n");  
  return dggsvd_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, iwork, info);
}
inline integer ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* alpha, float* beta, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, scomplex* work, float* rwork, integer* iwork, integer* info)
{
  printf(" Function cggsvd() has been deprecated. Please use cggsvd3() instead.\n");  
  return cggsvd_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, rwork, iwork, info);
}
inline integer ggsvd(char* jobu, char* jobv, char* jobq, integer* m, integer* n, integer* p, integer* k, integer* l, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* alpha, double* beta, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, dcomplex* work, double* rwork, integer* iwork, integer* info)
{
  printf(" Function dggsvd() has been deprecated. Please use dggsvd3() instead.\n");  
  return zggsvd_(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, rwork, iwork, info);
}

// --- computes orthogonal matrices U, V and Q ---
inline integer ggsvp3(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, float* a, integer* lda, float* b, integer* ldb, float* tola, float* tolb, integer* k, integer* l, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, integer* iwork, float* tau, float* work, integer* lwork, integer* info)
{
  return sggsvp3_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, lwork, info);
}
inline integer ggsvp3(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, double* a, integer* lda, double* b, integer* ldb, double* tola, double* tolb, integer* k, integer* l, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, integer* iwork, double* tau, double* work, integer* lwork, integer* info)
{
  return dggsvp3_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, lwork, info);
}
inline integer ggsvp3(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* tola, float* tolb, integer* k, integer* l, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, integer* iwork, float* rwork, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cggsvp3_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info);
}
inline integer ggsvp3(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* tola, double* tolb, integer* k, integer* l, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, integer* iwork, double* rwork, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zggsvp3_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info);
}

// --- computes orthogonal matrices U, V and Q ---
inline integer ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, float* a, integer* lda, float* b, integer* ldb, float* tola, float* tolb, integer* k, integer* l, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, integer* iwork, float* tau, float* work, integer* info)
{
  printf(" Function sggsvp() has been deprecated. Please use sggsvp3() instead.\n");  
  return sggsvp_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, info);
}
inline integer ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, double* a, integer* lda, double* b, integer* ldb, double* tola, double* tolb, integer* k, integer* l, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, integer* iwork, double* tau, double* work, integer* info)
{
  printf(" Function dggsvp() has been deprecated. Please use dggsvp3() instead.\n");  
  return dggsvp_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work, info);
}
inline integer ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* tola, float* tolb, integer* k, integer* l, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, integer* iwork, float* rwork, scomplex* tau, scomplex* work, integer* info)
{
  printf(" Function cggsvp() has been deprecated. Please use cggsvp3() instead.\n");  
  return cggsvp_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, info);
}
inline integer ggsvp(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* tola, double* tolb, integer* k, integer* l, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, integer* iwork, double* rwork, dcomplex* tau, dcomplex* work, integer* info)
{
  printf(" Function zggsvp() has been deprecated. Please use zggsvp3() instead.\n"); 
  return zggsvp_(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, info);
}

// --- estimates the reciprocal of the condition number ---
inline integer gtcon(char* norm, integer* n,  float* dl,  float* d, float* du, float* du2, integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  return sgtcon_(norm, n,  dl,  d,  du,  du2,  ipiv, anorm, rcond, work, iwork, info);
}
inline integer gtcon(char* norm, integer* n,  double* dl, double* d, double* du, double* du2, integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  return dgtcon_(norm, n,  dl,  d,  du,  du2,  ipiv, anorm, rcond, work, iwork, info);
}
inline integer gtcon(char* norm, integer* n,  scomplex* dl, scomplex* d, scomplex* du, scomplex* du2, integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  return cgtcon_(norm, n,  dl,  d,  du,  du2,  ipiv, anorm, rcond, work, info);
}
inline integer gtcon(char* norm, integer* n,  dcomplex* dl, dcomplex* d,  dcomplex* du, dcomplex* du2, integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  return zgtcon_(norm, n,  dl,  d,  du,  du2,  ipiv, anorm, rcond, work, info);
}

// --- improves the computed solution to a system of linear equations ---
inline integer gtrfs(char* trans, integer* n, integer* nrhs, float* dl, float* d,  float* du,  float* dlf, float* df, float* duf, float* du2, integer* ipiv, float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sgtrfs_(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer gtrfs(char* trans, integer* n, integer* nrhs,  double* dl,  double* d,  double* du,  double* dlf,  double* df,  double* duf,  double* du2,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dgtrfs_(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer gtrfs(char* trans, integer* n, integer* nrhs,  scomplex* dl,  scomplex* d,  scomplex* du, scomplex* dlf,  scomplex* df,  scomplex* duf,  scomplex* du2,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cgtrfs_(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer gtrfs(char* trans, integer* n, integer* nrhs,  dcomplex* dl,  dcomplex* d,  dcomplex* du,  dcomplex* dlf,  dcomplex* df,  dcomplex* duf,  dcomplex* du2,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zgtrfs_(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equations ---
inline integer gtsv(integer* n, integer* nrhs, float* dl, float* d, float* du, float* b, integer* ldb, integer* info)
{
  return sgtsv_(n, nrhs, dl, d, du, b, ldb, info);
}
inline integer gtsv(integer* n, integer* nrhs, double* dl, double* d, double* du, double* b, integer* ldb, integer* info)
{
  return dgtsv_(n, nrhs, dl, d, du, b, ldb, info);
}
inline integer gtsv(integer* n, integer* nrhs, scomplex* dl, scomplex* d, scomplex* du, scomplex* b, integer* ldb, integer* info)
{
  return cgtsv_(n, nrhs, dl, d, du, b, ldb, info);
}
inline integer gtsv(integer* n, integer* nrhs, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* b, integer* ldb, integer* info)
{
  return zgtsv_(n, nrhs, dl, d, du, b, ldb, info);
}

// --- uses the LU factorization to compute the solution ---
inline integer gtsvx(char* fact, char* trans, integer* n, integer* nrhs,  float* dl,  float* d, float* du, float* dlf, float* df, float* duf, float* du2, integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sgtsvx_(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer gtsvx(char* fact, char* trans, integer* n, integer* nrhs,  double* dl,  double* d,  double* du, double* dlf, double* df, double* duf, double* du2, integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dgtsvx_(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer gtsvx(char* fact, char* trans, integer* n, integer* nrhs,  scomplex* dl,  scomplex* d,  scomplex* du, scomplex* dlf, scomplex* df, scomplex* duf, scomplex* du2, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cgtsvx_(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer gtsvx(char* fact, char* trans, integer* n, integer* nrhs,  dcomplex* dl,  dcomplex* d,  dcomplex* du, dcomplex* dlf, dcomplex* df, dcomplex* duf, dcomplex* du2, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zgtsvx_(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes an LU factorization of a real tridiagonal matrix A ---
inline integer gttrf(integer* n, float* dl, float* d, float* du, float* du2, integer* ipiv, integer* info)
{
  return sgttrf_(n, dl, d, du, du2, ipiv, info);
}
inline integer gttrf(integer* n, double* dl, double* d, double* du, double* du2, integer* ipiv, integer* info)
{
  return dgttrf_(n, dl, d, du, du2, ipiv, info);
}
inline integer gttrf(integer* n, scomplex* dl, scomplex* d, scomplex* du, scomplex* du2, integer* ipiv, integer* info)
{
  return cgttrf_(n, dl, d, du, du2, ipiv, info);
}
inline integer gttrf(integer* n, dcomplex* dl, dcomplex* d, dcomplex* du, dcomplex* du2, integer* ipiv, integer* info)
{
  return zgttrf_(n, dl, d, du, du2, ipiv, info);
}

// --- solves one of the systems of equations ---
inline integer gttrs(char* trans, integer* n, integer* nrhs,  float* dl,  float* d,  float* du,  float* du2,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  return sgttrs_(trans, n, nrhs,  dl,  d,  du,  du2,  ipiv, b, ldb, info);
}
inline integer gttrs(char* trans, integer* n, integer* nrhs,  double* dl,  double* d,  double* du,  double* du2,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dgttrs_(trans, n, nrhs,  dl,  d,  du,  du2,  ipiv, b, ldb, info);
}
inline integer gttrs(char* trans, integer* n, integer* nrhs,  scomplex* dl,  scomplex* d,  scomplex* du,  scomplex* du2,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return cgttrs_(trans, n, nrhs,  dl,  d,  du,  du2,  ipiv, b, ldb, info);
}
inline integer gttrs(char* trans, integer* n, integer* nrhs,  dcomplex* dl,  dcomplex* d,  dcomplex* du,  dcomplex* du2,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zgttrs_(trans, n, nrhs,  dl,  d,  du,  du2,  ipiv, b, ldb, info);
}

// --- estimates the 1-norm of a square matrix ---
inline integer lacn2(integer* n, float* v, float* x, integer* isgn, float* est, integer* kase, integer* isave)
{
  return slacn2_(n, v, x, isgn, est, kase, isave);
}
inline integer lacn2(integer* n, double* v, double* x, integer* isgn, double* est, integer* kase, integer* isave)
{
  return dlacn2_(n, v, x, isgn, est, kase, isave);
}
inline integer lacn2(integer* n, scomplex* v, scomplex* x, float* est, integer* kase, integer* isave)
{
  return clacn2_(n, v, x, est, kase, isave);
}
inline integer lacn2(integer* n, dcomplex* v, dcomplex* x, double* est, integer* kase, integer* isave)
{
  return zlacn2_(n, v, x, est, kase, isave);
}

// --- copies all or part of one two-dimensional array ---
inline integer lacpy(char* uplo, integer* m, integer* n,  float* a, integer* lda, float* b, integer* ldb)
{
  return slacpy_(uplo, m, n,  a, lda, b, ldb);
}
inline integer lacpy(char* uplo, integer* m, integer* n,  double* a, integer* lda, double* b, integer* ldb)
{
  return dlacpy_(uplo, m, n,  a, lda, b, ldb);
}
inline integer lacpy(char* uplo, integer* m, integer* n,  scomplex* a, integer* lda, scomplex* b, integer* ldb)
{
  return clacpy_(uplo, m, n,  a, lda, b, ldb);
}
inline integer lacpy(char* uplo, integer* m, integer* n, dcomplex* a, integer* lda,dcomplex* b, integer* ldb)
{
  return zlacpy_(uplo, m, n,  a, lda, b, ldb);
}

inline integer slag2d(integer* m, integer* n, float* sa, integer* ldsa, double* a, integer* lda, integer* info)
{
  return slag2d_(m, n, sa, ldsa, a, lda, info);
}

// --- converts a double precision matrix to a single precision matrix ---
inline integer dlag2s(integer*m, integer* n, double* a, integer*lda, float* sa, integer* ldsa, integer* info)
{
  return dlag2s_(m, n, a, lda, sa, ldsa, info);
}
// --- converts a complex single precision matrix to a complex double precision matrix. ---
inline integer clag2z(integer*m, integer* n, scomplex* sa, integer*ldsa, dcomplex* a, integer* lda, integer* info)
{
  return clag2z_(m, n, sa, ldsa, a, lda, info);
} 
// --- converts a single precision matrix to a double precision matrix ---
inline integer zlag2c(integer*m, integer* n, dcomplex* a, integer*lda, scomplex* sa, integer* ldsa, integer* info)
{
  return zlag2c_(m, n, a, lda, sa, ldsa, info);
}

// --- returns sqrt(x2+y2) ---
inline float lapy2(float* x, float* y)
{
  return slapy2_(x, y);
}
inline double lapy2(double* x, double* y)
{
  return dlapy2_(x, y);
}

// --- returns sqrt(x2+y2+z2) ---
inline float lapy3(float* x, float* y, float* z)
{
  return slapy3_(x, y, z);
}
inline double lapy3(double* x, double* y, double* z)
{
  return dlapy3_(x, y, z);
}

// --- generates a plane rotation so that the diagonal is nonnegative ---
inline float lartgp(float* f, float* g, float* cs, float* sn, float* r)
{
  return slartgp_(f, g, cs, sn, r);
}
inline double lartgp(double* f, double* g, double* cs, double* sn, double* r)
{
  return dlartgp_(f, g, cs, sn, r);
}

// --- generates a plane rotation designed to introduce a bulge in implicit QR iteration for the bidiagonal SVD problem ---
inline integer lartgs(float* x, float* y, float* sigma, float* cs, float* sn)
{
  return slartgs_(x, y, sigma, cs, sn);
}
inline integer lartgs(double* x, double* y, double* sigma, double* cs, double* sn)
{
  return dlartgs_(x, y, sigma, cs, sn);
}

// ---  returns the value of the largest absolute value of any element of a general rectangular matrix ---
inline float lange(char* norm, integer* m, integer* n,  float* a, integer* lda, float* work)
{
  return slange_(norm, m, n,  a, lda, work);
}
inline double lange(char* norm, integer* m, integer* n,  double* a, integer* lda, double* work)
{
  return dlange_(norm, m, n,  a, lda, work);
}
inline float lange(char* norm, integer* m, integer* n,  scomplex* a, integer* lda, float* work)
{
  return clange_(norm, m, n,  a, lda, work);
}
inline double lange(char* norm, integer* m, integer* n,  dcomplex* a, integer* lda, double* work)
{
  return zlange_(norm, m, n,  a, lda, work);
}

// --- returns the largest absolute value of a real symmetric matrix---
inline float lansy(char* norm, char* uplo, integer* n,  float* a, integer* lda, float* work)
{
  return slansy_(norm, uplo, n,  a, lda, work);
}
inline double lansy(char* norm, char* uplo, integer* n,  double* a, integer* lda, double* work)
{
  return dlansy_(norm, uplo, n,  a, lda, work);
}
inline float lansy(char* norm, char* uplo, integer* n,  scomplex* a, integer* lda, scomplex* work)
{
  return clansy_(norm, uplo, n,  a, lda, work);
}
inline double lansy(char* norm, char* uplo, integer* n,  dcomplex* a, integer* lda, dcomplex* work)
{
  return zlansy_(norm, uplo, n,  a, lda, work);
}

// --- returns the largest absolute value of a trapezoidal or triangular matrix---
inline float lantr(char* norm, char* uplo, char* diag, integer* m, integer* n,  float* a, integer* lda, float* work)
{
  return slantr_(norm, uplo, diag, m, n,  a, lda, work);
}
inline double lantr(char* norm, char* uplo, char* diag, integer* m, integer* n,  double* a, integer* lda, double* work)
{
  return dlantr_(norm, uplo, diag, m, n,  a, lda, work);
}
inline float lantr(char* norm, char* uplo, char* diag, integer* m, integer* n,  scomplex* a, integer* lda, float* work)
{
  return clantr_(norm, uplo, diag, m, n,  a, lda, work);
}
inline double lantr(char* norm, char* uplo, char* diag, integer* m, integer* n,  dcomplex* a, integer* lda, double* work)
{
  return zlantr_(norm, uplo, diag, m, n,  a, lda, work);
}

// --- rearranges rows of a matrix as specified by a permutation vector. ---
inline integer lapmr(logical* forwrd, integer* m, integer* n, float* x, integer* ldx, integer* k)
{
  return slapmr_(forwrd, m, n, x, ldx, k);
}
inline integer lapmr(logical* forwrd, integer* m, integer* n, double* x, integer* ldx, integer* k)
{
  return dlapmr_(forwrd, m, n, x, ldx, k);
}
inline integer lapmr(logical* forwrd, integer* m, integer* n, scomplex* x, integer* ldx, integer* k)
{
  return clapmr_(forwrd, m, n, x, ldx, k);
}
inline integer lapmr(logical* forwrd, integer* m, integer* n, dcomplex* x, integer* ldx, integer* k)
{
  return zlapmr_(forwrd, m, n, x, ldx, k);
}

// --- performs a forward or backward permutation of the columns of a matrix ---
inline integer lapmt(logical* forwrd, integer* m, integer* n, float* x, integer* ldx, integer* k)
{
  return slapmt_(forwrd, m, n, x, ldx, k);
}
inline integer lapmt(logical* forwrd, integer* m, integer* n, double* x, integer* ldx, integer* k)
{
  return dlapmt_(forwrd, m, n, x, ldx, k);
}
inline integer lapmt(logical* forwrd, integer* m, integer* n, scomplex* x, integer* ldx, integer* k)
{
  return clapmt_(forwrd, m, n, x, ldx, k);
}
inline integer lapmt(logical* forwrd, integer* m, integer* n, dcomplex* x, integer* ldx, integer* k)
{
  return zlapmt_(forwrd, m, n, x, ldx, k);
}

// --- applies a block reflector or its transpose to a general rectangular matrix ---
inline integer larfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k,  float* v, integer* ldv,  float* t, integer* ldt, float* c, integer* ldc, float* work, integer* ldwork)
{
  return slarfb_(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline integer larfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k,  double* v, integer* ldv,  double* t, integer* ldt, double* c, integer* ldc, double* work, integer* ldwork)
{
  return dlarfb_(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline integer larfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k,  scomplex* v, integer* ldv,  scomplex* t, integer* ldt, scomplex* c, integer* ldc, scomplex* work, integer* ldwork)
{
  return clarfb_(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline integer larfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k,  dcomplex* v, integer* ldv,  dcomplex* t, integer* ldt, dcomplex* c, integer* ldc, dcomplex* work, integer* ldwork)
{
  return zlarfb_(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork);
}

// --- applies an elementary reflector to a general rectangular matrix ---
inline integer larfx(char* side, integer* m, integer* n, float* v, float* tau, float* c, integer* ldc, float* work)
{
  return slarfx_(side, m, n,  v, tau, c, ldc, work);
}
inline integer larfx(char* side, integer* m, integer* n, double* v, double* tau, double* c, integer* ldc, double* work)
{
  return dlarfx_(side, m, n,  v, tau, c, ldc, work);
}
inline integer larfx(char* side, integer* m, integer* n, scomplex* v, scomplex* tau, scomplex* c, integer* ldc, scomplex* work)
{
  return clarfx_(side, m, n,  v, tau, c, ldc, work);
}
inline integer larfx(char* side, integer* m, integer* n, dcomplex* v, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work)
{
  return zlarfx_(side, m, n,  v, tau, c, ldc, work);
}

// ---  returns a vector of random numbers from a uniform or normal distribution ---
inline integer larnv(integer* idist, integer* iseed, integer* n, float* x)
{
  return slarnv_(idist, iseed, n, x);
}
inline integer larnv(integer* idist, integer* iseed, integer* n, double* x)
{
  return dlarnv_(idist, iseed, n, x);
}
inline integer larnv(integer* idist, integer* iseed, integer* n, scomplex* x)
{
  return clarnv_(idist, iseed, n, x);
}
inline integer larnv(integer* idist, integer* iseed, integer* n, dcomplex* x)
{
  return zlarnv_(idist, iseed, n, x);
}

// --- multiplies a general rectangular matrix by a real scalar defined as cto/cfrom ---
inline integer lascl(char* type, integer* kl, integer* ku, float* cfrom, float* cto, integer* m, integer* n, float* a, integer* lda, integer* info)
{
  return slascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
}
inline integer lascl(char* type, integer* kl, integer* ku, double* cfrom, double* cto, integer* m, integer* n, double* a, integer* lda, integer* info)
{
  return dlascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
}
inline integer lascl(char* type, integer* kl, integer* ku, float* cfrom, float* cto, integer* m, integer* n, scomplex* a, integer* lda, integer* info)
{
  return clascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
}
inline integer lascl(char* type, integer* kl, integer* ku, double* cfrom, double* cto, integer* m, integer* n, dcomplex* a, integer* lda, integer* info)
{
  return zlascl_(type, kl, ku, cfrom, cto, m, n, a, lda, info);
}

// --- sorts numbers in increasing or decreasing order. ---
inline integer lasrt(char* id, integer* n, float* d, integer* info)
{
  return slasrt_(id, n, d, info);
}
inline integer lasrt(char* id, integer* n, double* d, integer* info)
{
  return dlasrt_(id, n, d, info);
}

// --- updates a sum of squares represented in scaled form ---
inline integer lassq(integer* n, float* x, integer* incx, float* scale, float* sumsq)
{
  return slassq_(n, x, incx, scale, sumsq);
}
inline integer lassq(integer* n, double* x, integer* incx, double* scale, double* sumsq)
{
  return dlassq_(n, x, incx, scale, sumsq);
}
inline integer lassq(integer* n, scomplex* x, integer* incx, float* scale, float* sumsq)
{
  return classq_(n, x, incx, scale, sumsq);
}
inline integer lassq(integer* n, dcomplex* x, integer* incx, double* scale, double* sumsq)
{
  return zlassq_(n, x, incx, scale, sumsq);
}

// --- generates a real orthogonal matrix Q ---
inline integer opgtr(char* uplo, integer* n, float* ap, float* tau, float* q, integer* ldq, float *work, integer *info)
{
  return sopgtr_(uplo, n, ap, tau, q, ldq, work, info);
}
inline integer opgtr(char* uplo, integer* n, double* ap, double* tau, double* q, integer* ldq, double *work, integer *info)
{
  return dopgtr_(uplo, n, ap, tau, q, ldq, work, info);
}
inline integer upgtr(char* uplo, integer* n, scomplex* ap, scomplex* tau, scomplex* q, integer* ldq, scomplex* work, integer* info)
{
  return cupgtr_(uplo, n, ap, tau, q, ldq, work, info);
}
inline integer upgtr(char* uplo, integer* n, dcomplex* ap, dcomplex* tau, dcomplex* q, integer* ldq, dcomplex* work, integer* info)
{
  return zupgtr_(uplo, n, ap, tau, q, ldq, work, info);
}

// --- overwrites the general real M-by-N matrix ---
inline integer opmtr(char* side, char* uplo, char* trans, integer* m, integer* n,  float* ap,  float* tau, float* c, integer* ldc, float* work, integer* info)
{
  return sopmtr_(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
}
inline integer opmtr(char* side, char* uplo, char* trans, integer* m, integer* n,  double* ap,  double* tau, double* c, integer* ldc, double* work, integer* info)
{
  return dopmtr_(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
}

inline integer upmtr(char* side, char* uplo, char* trans, integer* m, integer* n,  scomplex* ap, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  return cupmtr_(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
}
inline integer upmtr(char* side, char* uplo, char* trans, integer* m, integer* n,  dcomplex* ap,  dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  return zupmtr_(side, uplo, trans, m, n, ap, tau, c, ldc, work, info);
}

// --- simultaneously bidiagonalizes the blocks of an M-by-M partitioned orthogonal matrix X ---
inline integer orbdb(char* trans, char* signs, integer* m, integer* p, integer* q, float* x11, integer* ldx11, float* x12, integer* ldx12, float* x21, integer* ldx21, float* x22, integer* ldx22, float* theta, float* phi, float* taup1, float* taup2, float* tauq1, float* tauq2, float* work, integer* lwork, integer* info)
{
  return sorbdb_(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info);
}
inline integer orbdb(char* trans, char* signs, integer* m, integer* p, integer* q, double* x11, integer* ldx11, double* x12, integer* ldx12, double* x21, integer* ldx21, double* x22, integer* ldx22, double* theta, double* phi, double* taup1, double* taup2, double* tauq1, double* tauq2, double* work, integer* lwork, integer* info)
{
  return dorbdb_(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info);
}
inline integer unbdb(char* trans, char* signs, integer* m, integer* p, integer* q, scomplex* x11, integer* ldx11, scomplex* x12, integer* ldx12, scomplex* x21, integer* ldx21, scomplex* x22, integer* ldx22, float* theta, float* phi, scomplex* taup1, scomplex* taup2, scomplex* tauq1, scomplex* tauq2, scomplex* work, integer* lwork, integer* info)
{
  return cunbdb_(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info);
}
inline integer unbdb(char* trans, char* signs, integer* m, integer* p, integer* q, dcomplex* x11, integer* ldx11, dcomplex* x12, integer* ldx12, dcomplex* x21, integer* ldx21, dcomplex* x22, integer* ldx22, double* theta, double* phi, dcomplex* taup1, dcomplex* taup2, dcomplex* tauq1, dcomplex* tauq2, dcomplex* work, integer* lwork, integer* info)
{
  return zunbdb_(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info);
}

// --- computes the CS decomposition of an M-by-M partitioned orthogonal matrix X ---
inline integer orcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, integer* m, integer* p, integer* q, float* x11, integer* ldx11, float* x12, integer* ldx12, float* x21, integer* ldx21, float* x22, integer* ldx22, float* theta, float* u1, integer* ldu1, float* u2, integer* ldu2, float* v1t, integer* ldv1t, float* v2t, integer* ldv2t, float *work, integer *lwork, integer *iwork, integer *info)
{
  return sorcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, iwork, info);
}
inline integer orcsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, integer* m, integer* p, integer* q, double* x11, integer* ldx11, double* x12, integer* ldx12, double* x21, integer* ldx21, double* x22, integer* ldx22, double* theta, double* u1, integer* ldu1, double* u2, integer* ldu2, double* v1t, integer* ldv1t, double* v2t, integer* ldv2t, double *work, integer *lwork, integer *iwork, integer *info)
{
  return dorcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, iwork, info);
}
inline integer uncsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, integer* m, integer* p, integer* q, scomplex* x11, integer* ldx11, scomplex* x12, integer* ldx12, scomplex* x21, integer* ldx21, scomplex* x22, integer* ldx22, float* theta, scomplex* u1, integer* ldu1, scomplex* u2, integer* ldu2, scomplex* v1t, integer* ldv1t, scomplex* v2t, integer* ldv2t, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* info)
{
  return cuncsd_(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, rwork, lrwork, iwork, info);
}
inline integer uncsd(char* jobu1, char* jobu2, char* jobv1t, char* jobv2t, char* trans, char* signs, integer* m, integer* p, integer* q, dcomplex* x11, integer* ldx11, dcomplex* x12, integer* ldx12, dcomplex* x21, integer* ldx21, dcomplex* x22, integer* ldx22, double* theta, dcomplex* u1, integer* ldu1, dcomplex* u2, integer* ldu2, dcomplex* v1t, integer* ldv1t, dcomplex* v2t, integer* ldv2t, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* info)
{
  return zuncsd_(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, rwork, lrwork, iwork, info);
}

// --- computes the CS decomposition of an M-by-Q matrix X with orthonormal columns ---
inline integer orcsd2by1(char* jobu1, char* jobu2, char* jobv1t, integer* m, integer* p, integer* q, float* x11, integer* ldx11, float* x21, integer* ldx21, float* theta, float* u1, integer* ldu1, float* u2, integer* ldu2, float* v1t, integer* ldv1t, float *work, integer *lwork, integer *iwork, integer *info)
{
  return sorcsd2by1_(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork, iwork, info);
}
inline integer orcsd2by1(char* jobu1, char* jobu2, char* jobv1t, integer* m, integer* p, integer* q, double* x11, integer* ldx11, double* x21, integer* ldx21, double* theta, double* u1, integer* ldu1, double* u2, integer* ldu2, double* v1t, integer* ldv1t, double *work, integer *lwork, integer *iwork, integer *info)
{
  return dorcsd2by1_(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork, iwork, info);
}
inline integer uncsd2by1(char* jobu1, char* jobu2, char* jobv1t, integer* m, integer* p, integer* q, scomplex* x11, integer* ldx11, scomplex* x21, integer* ldx21, float* theta, scomplex* u1, integer* ldu1, scomplex* u2, integer* ldu2, scomplex* v1t, integer* ldv1t, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* info)
{
  return cuncsd2by1_(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork, rwork, lrwork, iwork, info);
}
inline integer uncsd2by1(char* jobu1, char* jobu2, char* jobv1t, integer* m, integer* p, integer* q, dcomplex* x11, integer* ldx11, dcomplex* x21, integer* ldx21, double* theta, dcomplex* u1, integer* ldu1, dcomplex* u2, integer* ldu2, dcomplex* v1t, integer* ldv1t, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* info)
{
  return zuncsd2by1_(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork, rwork, lrwork, iwork, info);
}

// --- generates an M-by-N real matrix Q with orthonormal columns ---
inline integer orgql(integer* m, integer* n, integer* k, float* a, integer* lda,  float* tau, float* work, integer* lwork, integer* info)
{
  return sorgql_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer orgql(integer* m, integer* n, integer* k, double* a, integer* lda,  double* tau, double* work, integer* lwork, integer* info)
{
  return dorgql_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer ungql(integer* m, integer* n, integer* k, scomplex* a, integer* lda,  scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cungql_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer ungql(integer* m, integer* n, integer* k, dcomplex* a, integer* lda,  dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zungql_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- generates an M-by-N real matrix Q with orthonormal rows ---
inline integer orgrq(integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return sorgrq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer orgrq(integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dorgrq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer ungrq(integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cungrq_(m, n, k, a, lda, tau, work, lwork, info);
}
inline integer ungrq(integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zungrq_(m, n, k, a, lda, tau, work, lwork, info);
}

// --- overwrites the general real M-by-N matrix ---
inline integer ormql(char* side, char* trans, integer* m, integer* n, integer* k, float* a, integer* lda, float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  return sormql_(side, trans, m, n, k,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline integer ormql(char* side, char* trans, integer* m, integer* n, integer* k, double* a, integer* lda, double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  return dormql_(side, trans, m, n, k,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline integer unmql(char* side, char* trans, integer* m, integer* n, integer* k, scomplex* a, integer* lda, scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cunmql_(side, trans, m, n, k,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline integer unmql(char* side, char* trans, integer* m, integer* n, integer* k, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zunmql_(side, trans, m, n, k,  a, lda,  tau, c, ldc, work, lwork, info);
}

// --- overwrites the general real M-by-N matrix ---
inline integer ormrq(char* side, char* trans, integer* m, integer* n, integer* k,  float* a, integer* lda,  float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  return sormrq_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}
inline integer ormrq(char* side, char* trans, integer* m, integer* n, integer* k,  double* a, integer* lda,  double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  return dormrq_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}
inline integer unmrq(char* side, char* trans, integer* m, integer* n, integer* k,  scomplex* a, integer* lda,  scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cunmrq_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}
inline integer unmrq(char* side, char* trans, integer* m, integer* n, integer* k,  dcomplex* a, integer* lda,  dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zunmrq_(side, trans, m, n, k,  a, lda, tau, c, ldc, work, lwork, info);
}

// --- overwrites the general real M-by-N matrix---
inline integer ormrz(char* side, char* trans, integer* m, integer* n, integer* k, integer* l,  float* a, integer* lda,  float* tau, float* c, integer* ldc, float* work, integer* lwork, integer* info)
{
  return sormrz_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer ormrz(char* side, char* trans, integer* m, integer* n, integer* k, integer* l,  double* a, integer* lda,  double* tau, double* c, integer* ldc, double* work, integer* lwork, integer* info)
{
  return dormrz_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer unmrz(char* side, char* trans, integer* m, integer* n, integer* k, integer* l,  scomplex* a, integer* lda,  scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cunmrz_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
}
inline integer unmrz(char* side, char* trans, integer* m, integer* n, integer* k, integer* l,  dcomplex* a, integer* lda,  dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zunmrz_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info);
}

inline integer latzm(char* side, integer* m, integer* n, float* v, integer* incv, float* tau, float* c1, float* c2, integer* ldc, float* work)
{
  return slatzm_(side, m, n, v, incv, tau, c1, c2, ldc, work);
}
inline integer latzm(char* side, integer* m, integer* n, double* v, integer* incv, double* tau, double* c1, double* c2, integer* ldc, double* work)
{
  return dlatzm_(side, m, n, v, incv, tau, c1, c2, ldc, work);
}
inline integer latzm(char* side, integer* m, integer* n, scomplex* v, integer* incv, scomplex* tau, scomplex* c1, scomplex* c2, integer* ldc, scomplex* work)
{
  printf(" Function clatzm() has been deprecated. Please use cunmrz() instead.\n"); 
  return clatzm_(side, m, n, v, incv, tau, c1, c2, ldc, work);
}
inline integer latzm(char* side, integer* m, integer* n, dcomplex* v, integer* incv, dcomplex* tau, dcomplex* c1, dcomplex* c2, integer* ldc, dcomplex* work)
{
  printf(" Function zlatzm() has been deprecated. Please use zunmrz() instead.\n"); 
  return zlatzm_(side, m, n, v, incv, tau, c1, c2, ldc, work);
}

// --- generates a real orthogonal matrix Q ---
inline integer orghr(integer* n, integer* ilo, integer* ihi, float* a, integer* lda, float* tau, float *work, integer *lwork, integer *info)
{
  return sorghr_(n, ilo, ihi, a, lda,  tau, work, lwork, info);
}
inline integer orghr(integer* n, integer* ilo, integer* ihi, double* a, integer* lda, double* tau, double *work, integer *lwork, integer *info)
{
  return dorghr_(n, ilo, ihi, a, lda,  tau, work, lwork, info);
}
inline integer unghr(integer* n, integer* ilo, integer* ihi, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return cunghr_(n, ilo, ihi, a, lda,  tau, work, lwork, info);
}
inline integer unghr(integer* n, integer* ilo, integer* ihi, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zunghr_(n, ilo, ihi, a, lda,  tau, work, lwork, info);
}

// --- overwrites the general complex M-by-N matrix C ---
inline integer ormhr(char* side, char* trans, integer* m, integer* n, integer* ilo, integer* ihi,  float* a, integer* lda,  float* tau, float* c, integer* ldc, float *work, integer *lwork, integer *info)
{
  return sormhr_(side, trans, m, n, ilo, ihi,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline integer ormhr(char* side, char* trans, integer* m, integer* n, integer* ilo, integer* ihi,  double* a, integer* lda,  double* tau, double* c, integer* ldc, double *work, integer *lwork, integer *info)
{
  return dormhr_(side, trans, m, n, ilo, ihi,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline integer unmhr(char* side, char* trans, integer* m, integer* n, integer* ilo, integer* ihi,  scomplex* a, integer* lda,  scomplex* tau, scomplex* c, integer* ldc, scomplex* work, integer* lwork, integer* info)
{
  return cunmhr_(side, trans, m, n, ilo, ihi,  a, lda,  tau, c, ldc, work, lwork, info);
}
inline integer unmhr(char* side, char* trans, integer* m, integer* n, integer* ilo, integer* ihi,  dcomplex* a, integer* lda,  dcomplex* tau, dcomplex* c, integer* ldc, dcomplex* work, integer* lwork, integer* info)
{
  return zunmhr_(side, trans, m, n, ilo, ihi,  a, lda,  tau, c, ldc, work, lwork, info);
}

// --- estimates the reciprocal of the condition number ---
inline integer pbcon(char* uplo, integer* n, integer* kd,  float* ab, integer* ldab, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  return spbcon_(uplo, n, kd,  ab, ldab, anorm, rcond, work, iwork, info);
}
inline integer pbcon(char* uplo, integer* n, integer* kd,  double* ab, integer* ldab, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  return dpbcon_(uplo, n, kd,  ab, ldab, anorm, rcond, work, iwork, info);
}
inline integer pbcon(char* uplo, integer* n, integer* kd,  scomplex* ab, integer* ldab, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  return cpbcon_(uplo, n, kd,  ab, ldab, anorm, rcond, work, rwork, info);
}
inline integer pbcon(char* uplo, integer* n, integer* kd,  dcomplex* ab, integer* ldab, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  return zpbcon_(uplo, n, kd,  ab, ldab, anorm, rcond, work, rwork, info);
}

// --- computes row and column scalings---
inline integer pbequ(char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* s, float* scond, float* amax, integer* info)
{
  return spbequ_(uplo, n, kd,  ab, ldab, s, scond, amax, info);
}
inline integer pbequ(char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* s, double* scond, double* amax, integer* info)
{
  return dpbequ_(uplo, n, kd,  ab, ldab, s, scond, amax, info);
}
inline integer pbequ(char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* s, float* scond, float* amax, integer* info)
{
  return cpbequ_(uplo, n, kd,  ab, ldab, s, scond, amax, info);
}
inline integer pbequ(char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* s, double* scond, double* amax, integer* info)
{
  return zpbequ_(uplo, n, kd,  ab, ldab, s, scond, amax, info);
}

// --- improves the computed solution to a system of linear equations ---
inline integer pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  float* ab, integer* ldab,  float* afb, integer* ldafb,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return spbrfs_(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  double* ab, integer* ldab,  double* afb, integer* ldafb,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dpbrfs_(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  scomplex* ab, integer* ldab,  scomplex* afb, integer* ldafb,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cpbrfs_(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer pbrfs(char* uplo, integer* n, integer* kd, integer* nrhs,  dcomplex* ab, integer* ldab,  dcomplex* afb, integer* ldafb,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zpbrfs_(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- computes a split Cholesky factorization of a real symmetric positive definite band matrix---
inline integer pbstf(char* uplo, integer* n, integer* kb, float* bb, integer* ldbb, integer* info)
{
  return spbstf_(uplo, n, kb, bb, ldbb, info);
}
inline integer pbstf(char* uplo, integer* n, integer* kb, double* bb, integer* ldbb, integer* info)
{
  return dpbstf_(uplo, n, kb, bb, ldbb, info);
}
inline integer pbstf(char* uplo, integer* n, integer* kb, scomplex* bb, integer* ldbb, integer* info)
{
  return cpbstf_(uplo, n, kb, bb, ldbb, info);
}
inline integer pbstf(char* uplo, integer* n, integer* kb, dcomplex* bb, integer* ldbb, integer* info)
{
  return zpbstf_(uplo, n, kb, bb, ldbb, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline integer pbsv(char* uplo, integer* n, integer* kd, integer* nrhs, float* ab, integer* ldab, float* b, integer* ldb, integer* info)
{
  return spbsv_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline integer pbsv(char* uplo, integer* n, integer* kd, integer* nrhs, double* ab, integer* ldab, double* b, integer* ldb, integer* info)
{
  return dpbsv_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline integer pbsv(char* uplo, integer* n, integer* kd, integer* nrhs, scomplex* ab, integer* ldab, scomplex* b, integer* ldb, integer* info)
{
  return cpbsv_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline integer pbsv(char* uplo, integer* n, integer* kd, integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* b, integer* ldb, integer* info)
{
  return zpbsv_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}

// --- computes the solution to system of linear equations ---
inline integer pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, float* ab, integer* ldab, float* afb, integer* ldafb, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return spbsvx_(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, double* ab, integer* ldab, double* afb, integer* ldafb, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dpbsvx_(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, scomplex* ab, integer* ldab, scomplex* afb, integer* ldafb, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cpbsvx_(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer pbsvx(char* fact, char* uplo, integer* n, integer* kd, integer* nrhs, dcomplex* ab, integer* ldab, dcomplex* afb, integer* ldafb, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zpbsvx_(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the Cholesky factorization ---
inline integer pbtrf(char* uplo, integer* n, integer* kd, float* ab, integer*ldab, integer* info)
{
  return spbtrf_(uplo, n, kd, ab, ldab, info);
}
inline integer pbtrf(char* uplo, integer* n, integer* kd, double* ab, integer*ldab, integer* info)
{
  return dpbtrf_(uplo, n, kd, ab, ldab, info);
}
inline integer pbtrf(char* uplo, integer* n, integer* kd, scomplex* ab, integer*ldab, integer* info)
{
  return cpbtrf_(uplo, n, kd, ab, ldab, info);
}
inline integer pbtrf(char* uplo, integer* n, integer* kd, dcomplex* ab, integer*ldab, integer* info)
{
  return zpbtrf_(uplo, n, kd, ab, ldab, info);
}

// --- solves a system of linear equations ---
inline integer pbtrs(char* uplo, integer* n, integer* kd, integer* nrhs,  float* ab, integer* ldab, float* b, integer* ldb, integer* info)
{
  return spbtrs_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline integer pbtrs(char* uplo, integer* n, integer* kd, integer* nrhs,  double* ab, integer* ldab, double* b, integer* ldb, integer* info)
{
  return dpbtrs_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline integer pbtrs(char* uplo, integer* n, integer* kd, integer* nrhs,  scomplex* ab, integer* ldab, scomplex* b, integer* ldb, integer* info)
{
  return cpbtrs_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}
inline integer pbtrs(char* uplo, integer* n, integer* kd, integer* nrhs,  dcomplex* ab, integer* ldab, dcomplex* b, integer* ldb, integer* info)
{
  return zpbtrs_(uplo, n, kd, nrhs, ab, ldab, b, ldb, info);
}

// --- computes the Cholesky factorization of a real symmetric ---
inline integer pftrf(char* transr, char* uplo, integer* n, float* a, integer* info)
{
  return spftrf_(transr, uplo, n, a, info);
}
inline integer pftrf(char* transr, char* uplo, integer* n, double* a, integer* info)
{
  return dpftrf_(transr, uplo, n, a, info);
}
inline integer pftrf(char* transr, char* uplo, integer* n, scomplex* a, integer* info)
{
  return cpftrf_(transr, uplo, n, a, info);
}
inline integer pftrf(char* transr, char* uplo, integer* n, dcomplex* a, integer* info)
{
  return zpftrf_(transr, uplo, n, a, info);
}

// --- computes the inverse of a real (symmetric) positive definite matrix ---
inline integer pftri(char* transr, char* uplo, integer* n, float* a, integer* info)
{
  return spftri_(transr, uplo, n, a, info);
}
inline integer pftri(char* transr, char* uplo, integer* n, double* a, integer* info)
{
  return dpftri_(transr, uplo, n, a, info);
}
inline integer pftri(char* transr, char* uplo, integer* n, scomplex* a, integer* info)
{
  return cpftri_(transr, uplo, n, a, info);
}
inline integer pftri(char* transr, char* uplo, integer* n, dcomplex* a, integer* info)
{
  return zpftri_(transr, uplo, n, a, info);
}

// --- solves a system of linear equations A*X = B with a symmetric matrix ---
inline integer pftrs(char* transr, char* uplo, integer* n, integer* nrhs,  float* a, float* b, integer* ldb, integer* info)
{
  return spftrs_(transr, uplo, n, nrhs,  a, b, ldb, info);
}
inline integer pftrs(char* transr, char* uplo, integer* n, integer* nrhs,  double* a, double* b, integer* ldb, integer* info)
{
  return dpftrs_(transr, uplo, n, nrhs,  a, b, ldb, info);
}
inline integer pftrs(char* transr, char* uplo, integer* n, integer* nrhs, scomplex* a, scomplex* b, integer* ldb, integer* info)
{
  return cpftrs_(transr, uplo, n, nrhs,  a, b, ldb, info);
}
inline integer pftrs(char* transr, char* uplo, integer* n, integer* nrhs, dcomplex* a, dcomplex* b, integer* ldb, integer* info)
{
  return zpftrs_(transr, uplo, n, nrhs,  a, b, ldb, info);
}

// --- estimates the reciprocal of the condition number ---
inline integer pocon(char* uplo, integer* n,  float* a, integer* lda, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  return spocon_(uplo, n,  a, lda, anorm, rcond, work, iwork, info);
}
inline integer pocon(char* uplo, integer* n,  double* a, integer* lda, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  return dpocon_(uplo, n,  a, lda, anorm, rcond, work, iwork, info);
}
inline integer pocon(char* uplo, integer* n,  scomplex* a, integer* lda, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  return cpocon_(uplo, n,  a, lda, anorm, rcond, work, rwork, info);
}
inline integer pocon(char* uplo, integer* n,  dcomplex* a, integer* lda, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  return zpocon_(uplo, n,  a, lda, anorm, rcond, work, rwork, info);
}

// --- computes row and column scalings ---
inline integer poequ(integer* n,  float* a, integer* lda, float* s, float* scond, float* amax, integer* info)
{
  return spoequ_(n, a, lda, s, scond, amax, info);
}
inline integer poequ(integer* n,  double* a, integer* lda, double* s, double* scond, double* amax, integer* info)
{
  return dpoequ_(n, a, lda, s, scond, amax, info);
}
inline integer poequ(integer* n,  scomplex* a, integer* lda, float* s, float* scond, float* amax, integer* info)
{
  return cpoequ_(n, a, lda, s, scond, amax, info);
}
inline integer poequ(integer* n,  dcomplex* a, integer* lda, double* s, double* scond, double* amax, integer* info)
{
  return zpoequ_(n, a, lda, s, scond, amax, info);
}

// --- computes row and column scalings ---
inline integer poequb(integer* n,  float* a, integer* lda, float* s, float* scond, float* amax, integer* info)
{
  return spoequb_(n, a, lda, s, scond, amax, info);
}
inline integer poequb(integer* n,  double* a, integer* lda, double* s, double* scond, double* amax, integer* info)
{
  return dpoequb_(n, a, lda, s, scond, amax, info);
}
inline integer poequb(integer* n,  scomplex* a, integer* lda, float* s, float* scond, float* amax, integer* info)
{
  return cpoequb_(n, a, lda, s, scond, amax, info);
}
inline integer poequb(integer* n,  dcomplex* a, integer* lda, double* s, double* scond, double* amax, integer* info)
{
  return zpoequb_(n, a, lda, s, scond, amax, info);
}

// --- improves the computed solution to a system of linear equations ---
inline integer porfs(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sporfs_(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer porfs(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dporfs_(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer porfs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cporfs_(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer porfs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zporfs_(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- improves the computed solution to a system of linear equations ---
inline integer porfsx(char* uplo, char* equed, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  float* s,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  return sporfsx_(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer porfsx(char* uplo, char* equed, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  double* s,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  return dporfsx_(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer porfsx(char* uplo, char* equed, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  float* s,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return cporfsx_(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer porfsx(char* uplo, char* equed, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  double* s,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zporfsx_(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for PO matrices ---
inline integer posv(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  return sposv_(uplo, n, nrhs, a, lda, b, ldb, info);
}
inline integer posv(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  return dposv_(uplo, n, nrhs, a, lda, b, ldb, info);
}
inline integer posv(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  return cposv_(uplo, n, nrhs, a, lda, b, ldb, info);
}
inline integer posv(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  return zposv_(uplo, n, nrhs, a, lda, b, ldb, info);
}

// --- computes the solution to system of linear equations A * X = B for PO matrices ---
inline integer posvx(char* fact, char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sposvx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer posvx(char* fact, char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dposvx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer posvx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cposvx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer posvx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zposvx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for PO matrices ---
inline integer posvxx(char* fact, char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  return sposvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer posvxx(char* fact, char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  return dposvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer posvxx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return cposvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer posvxx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zposvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// ---computes the Cholesky factorization of a real symmetric matrix A---
inline integer potrf2(char* uplo, integer* n, float* a, integer* lda, integer* info)
{
  return spotrf2_(uplo, n, a, lda, info);
}
inline integer potrf2(char* uplo, integer* n, double* a, integer* lda, integer* info)
{
  return dpotrf2_(uplo, n, a, lda, info);
}
inline integer potrf2(char* uplo, integer* n, scomplex* a, integer* lda, integer* info)
{
  return cpotrf2_(uplo, n, a, lda, info);
}
inline integer potrf2(char* uplo, integer* n, dcomplex* a, integer* lda, integer* info)
{
  return zpotrf2_(uplo, n, a, lda, info);
}

// --- solves a system of linear equations A*X = B ---
inline integer potrs(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  return spotrs_(uplo, n, nrhs,  a, lda, b, ldb, info);
}
inline integer potrs(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  return dpotrs_(uplo, n, nrhs,  a, lda, b, ldb, info);
}
inline integer potrs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  return cpotrs_(uplo, n, nrhs,  a, lda, b, ldb, info);
}
inline integer potrs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  return zpotrs_(uplo, n, nrhs,  a, lda, b, ldb, info);
}

// --- estimates the reciprocal of the condition number ---
inline integer ppcon(char* uplo, integer* n, float* ap, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  return sppcon_(uplo, n,  ap, anorm, rcond, work, iwork, info);
}
inline integer ppcon(char* uplo, integer* n, double* ap, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  return dppcon_(uplo, n,  ap, anorm, rcond, work, iwork, info);
}
inline integer ppcon(char* uplo, integer* n, scomplex* ap, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  return cppcon_(uplo, n,  ap, anorm, rcond, work, rwork, info);
}
inline integer ppcon(char* uplo, integer* n, dcomplex* ap, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  return zppcon_(uplo, n,  ap, anorm, rcond, work, rwork, info);
}

// --- computes row and column scalings ---
inline integer ppequ(char* uplo, integer* n, float* ap, float* s, float* scond, float* amax, integer* info)
{
  return sppequ_(uplo, n, ap, s, scond, amax, info);
}
inline integer ppequ(char* uplo, integer* n, double* ap, double* s, double* scond, double* amax, integer* info)
{
  return dppequ_(uplo, n, ap, s, scond, amax, info);
}
inline integer ppequ(char* uplo, integer* n, scomplex* ap, float* s, float* scond, float* amax, integer* info)
{
  return cppequ_(uplo, n, ap, s, scond, amax, info);
}
inline integer ppequ(char* uplo, integer* n, dcomplex* ap, double* s, double* scond, double* amax, integer* info)
{
  return zppequ_(uplo, n, ap, s, scond, amax, info);
}

// --- improves the computed solution to a system of linear equations ---
inline integer pprfs(char* uplo, integer* n, integer* nrhs,  float* ap,  float* afp,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return spprfs_(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer pprfs(char* uplo, integer* n, integer* nrhs,  double* ap,  double* afp,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dpprfs_(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer pprfs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  scomplex* afp,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cpprfs_(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer pprfs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  dcomplex* afp,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zpprfs_(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline integer ppsv(char* uplo, integer* n, integer* nrhs, float* ap, float* b, integer* ldb, integer* info)
{
  return sppsv_(uplo, n, nrhs, ap, b, ldb, info);
}
inline integer ppsv(char* uplo, integer* n, integer* nrhs, double* ap, double* b, integer* ldb, integer* info)
{
  return dppsv_(uplo, n, nrhs, ap, b, ldb, info);
}
inline integer ppsv(char* uplo, integer* n, integer* nrhs, scomplex* ap, scomplex* b, integer* ldb, integer* info)
{
  return cppsv_(uplo, n, nrhs, ap, b, ldb, info);
}
inline integer ppsv(char* uplo, integer* n, integer* nrhs, dcomplex* ap, dcomplex* b, integer* ldb, integer* info)
{
  return zppsv_(uplo, n, nrhs, ap, b, ldb, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline integer ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, float* ap, float* afp, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sppsvx_(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, double* ap, double* afp, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dppsvx_(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* ap, scomplex* afp, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cppsvx_(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer ppsvx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* ap, dcomplex* afp, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zppsvx_(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the Cholesky factorization of a real symmetric matrix ---
inline integer pptrf(char* uplo, integer* n, float* ap, integer* info)
{
  return spptrf_(uplo, n, ap, info);
}
inline integer pptrf(char* uplo, integer* n, double* ap, integer* info)
{
  return dpptrf_(uplo, n, ap, info);
}
inline integer pptrf(char* uplo, integer* n, scomplex* ap, integer* info)
{
  return cpptrf_(uplo, n, ap, info);
}
inline integer pptrf(char* uplo, integer* n, dcomplex* ap, integer* info)
{
  return zpptrf_(uplo, n, ap, info);
}

// --- computes the inverse of a real symmetric matrix ---
inline integer pptri(char* uplo, integer* n, float* ap, integer* info)
{
  return spptri_(uplo, n, ap, info);
}
inline integer pptri(char* uplo, integer* n, double* ap, integer* info)
{
  return dpptri_(uplo, n, ap, info);
}
inline integer pptri(char* uplo, integer* n, scomplex* ap, integer* info)
{
  return cpptri_(uplo, n, ap, info);
}
inline integer pptri(char* uplo, integer* n, dcomplex* ap, integer* info)
{
  return zpptri_(uplo, n, ap, info);
}

// --- solves a system of linear equations A*X = B with a symmetric matrix ---
inline integer pptrs(char* uplo, integer* n, integer* nrhs,  float* ap, float* b, integer* ldb, integer* info)
{
  return spptrs_(uplo, n, nrhs, ap, b, ldb, info);
}
inline integer pptrs(char* uplo, integer* n, integer* nrhs,  double* ap, double* b, integer* ldb, integer* info)
{
  return dpptrs_(uplo, n, nrhs, ap, b, ldb, info);
}
inline integer pptrs(char* uplo, integer* n, integer* nrhs,  scomplex* ap, scomplex* b, integer* ldb, integer* info)
{
  return cpptrs_(uplo, n, nrhs, ap, b, ldb, info);
}
inline integer pptrs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap, dcomplex* b, integer* ldb, integer* info)
{
  return zpptrs_(uplo, n, nrhs, ap, b, ldb, info);
}

// --- computes the Cholesky factorization ---
inline integer pstrf(char* uplo, integer* n, float* a, integer* lda, integer* piv, integer* rank, float* tol, float* work, integer* info)
{
  return spstrf_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline integer pstrf(char* uplo, integer* n, double* a, integer* lda, integer* piv, integer* rank, double* tol, double* work, integer* info)
{
  return dpstrf_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline integer pstrf(char* uplo, integer* n, scomplex* a, integer* lda, integer* piv, integer* rank, float* tol, float* work, integer* info)
{
  return cpstrf_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline integer pstrf(char* uplo, integer* n, dcomplex* a, integer* lda, integer* piv, integer* rank, double* tol, double* work, integer* info)
{
  return zpstrf_(uplo, n, a, lda, piv, rank, tol, work, info);
}

// --- computes the reciprocal of the condition number ---
inline integer ptcon(integer* n,  float* d,  float* e, float* anorm, float* rcond, float* rwork, integer* info)
{
  return sptcon_(n, d, e, anorm, rcond, rwork, info);
}
inline integer ptcon(integer* n,  double* d,  double* e, double* anorm, double* rcond, double* rwork, integer* info)
{
  return dptcon_(n, d, e, anorm, rcond, rwork, info);
}
inline integer ptcon(integer* n,  float* d,  scomplex* e, float* anorm, float* rcond, float* rwork, integer* info)
{
  return cptcon_(n, d, e, anorm, rcond, rwork, info);
}
inline integer ptcon(integer* n,  double* d,  dcomplex* e, double* anorm, double* rcond, double* rwork, integer* info)
{
  return zptcon_(n, d, e, anorm, rcond, rwork, info);
}

// --- computes all eigenvalues and, optionally, eigenvectors of a symmetric matrix ---
inline integer pteqr(char* compz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* info)
{
  return spteqr_(compz, n, d, e, z, ldz, work, info);
}
inline integer pteqr(char* compz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* info)
{
  return dpteqr_(compz, n, d, e, z, ldz, work, info);
}
inline integer pteqr(char* compz, integer* n, float* d, float* e, scomplex* z, integer* ldz, float* work, integer* info)
{
  return cpteqr_(compz, n, d, e, z, ldz, work, info);
}
inline integer pteqr(char* compz, integer* n, double* d, double* e, dcomplex* z, integer* ldz, double* work, integer* info)
{
  return zpteqr_(compz, n, d, e, z, ldz, work, info);
}

// --- improves the computed solution to a system of linear equations---
inline integer ptrfs(integer* n, integer* nrhs,  float* d,  float* e,  float* df,  float* ef,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* info)
{
  return sptrfs_(n, nrhs,  d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, info);
}
inline integer ptrfs(integer* n, integer* nrhs,  double* d,  double* e,  double* df,  double* ef,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* info)
{
  return dptrfs_(n, nrhs,  d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, info);
}
inline integer ptrfs(char *uplo, integer* n, integer* nrhs,  float* d,  scomplex* e,  float* df,  scomplex* ef,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cptrfs_(uplo, n, nrhs,  d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer ptrfs(char *uplo, integer* n, integer* nrhs,  double* d,  dcomplex* e,  double* df,  dcomplex* ef,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zptrfs_(uplo, n, nrhs,  d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for PT matrices ---
inline integer ptsv(integer* n, integer* nrhs, float* d, float* e, float* b, integer* ldb, integer* info)
{
  return sptsv_(n, nrhs, d, e, b, ldb, info);
}
inline integer ptsv(integer* n, integer* nrhs, double* d, double* e, double* b, integer* ldb, integer* info)
{
  return dptsv_(n, nrhs, d, e, b, ldb, info);
}
inline integer ptsv(integer* n, integer* nrhs, float* d, scomplex* e, scomplex* b, integer* ldb, integer* info)
{
  return cptsv_(n, nrhs, d, e, b, ldb, info);
}
inline integer ptsv(integer* n, integer* nrhs, double* d, dcomplex* e, dcomplex* b, integer* ldb, integer* info)
{
  return zptsv_(n, nrhs, d, e, b, ldb, info);
}

// --- computes the solution to system of linear equations A * X = B for PT matrices ---
inline integer ptsvx(char* fact, integer* n, integer* nrhs,  float* d,  float* e, float* df, float* ef,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* info)
{
  return sptsvx_(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, info);
}
inline integer ptsvx(char* fact, integer* n, integer* nrhs,  double* d,  double* e, double* df, double* ef,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* info)
{
  return dptsvx_(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, info);
}
inline integer ptsvx(char* fact, integer* n, integer* nrhs,  float* d,  scomplex* e, float* df, scomplex* ef,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cptsvx_(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer ptsvx(char* fact, integer* n, integer* nrhs,  double* d,  dcomplex* e, double* df, dcomplex* ef,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zptsvx_(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- computes the L*D*L**T factorization of a real symmetric matrix ---
inline integer pttrf(integer* n, float* d, float* e, integer* info)
{
  return spttrf_(n, d, e, info);
}
inline integer pttrf(integer* n, double* d, double* e, integer* info)
{
  return dpttrf_(n, d, e, info);
}
inline integer pttrf(integer* n, float* d, scomplex* e, integer* info)
{
  return cpttrf_(n, d, e, info);
}
inline integer pttrf(integer* n, double* d, dcomplex* e, integer* info)
{
  return zpttrf_(n, d, e, info);
}

// --- solves a tridiagonal system of the form A * X = B---
inline integer pttrs(integer* n, integer* nrhs,  float* d,  float* e, float* b, integer* ldb, integer* info)
{
  return spttrs_(n, nrhs,  d, e, b, ldb, info);
}
inline integer pttrs(integer* n, integer* nrhs,  double* d,  double* e, double* b, integer* ldb, integer* info)
{
  return dpttrs_(n, nrhs,  d, e, b, ldb, info);
}
inline integer pttrs(char *uplo, integer* n, integer* nrhs,  float* d,  scomplex* e, scomplex* b, integer* ldb, integer* info)
{
  return cpttrs_(uplo, n, nrhs,  d, e, b, ldb, info);
}
inline integer pttrs(char *uplo, integer* n, integer* nrhs,  double* d,  dcomplex* e, dcomplex* b, integer* ldb, integer* info)
{
  return zpttrs_(uplo, n, nrhs,  d, e, b, ldb, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer sbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* info)
{
  return ssbev_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, info);
}
inline integer sbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* info)
{
  return dsbev_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, info);
}
inline integer hbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return chbev_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, info);
}
inline integer hbev_2stage(char* jobz, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zhbev_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer sbev(char* jobz, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* w, float* z, integer* ldz, float* work, integer* info)
{
  return ssbev_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, info);
}
inline integer sbev(char* jobz, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* w, double* z, integer* ldz, double* work, integer* info)
{
  return dsbev_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, info);
}
inline integer hbev(char* jobz, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* info)
{
  return chbev_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, rwork, info);
}
inline integer hbev(char* jobz, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* info)
{
  return zhbev_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer sbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ssbevd_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer sbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dsbevd_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer hbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return chbevd_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer hbevd_2stage(char* jobz, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zhbevd_2stage_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer sbevd(char* jobz, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ssbevd_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer sbevd(char* jobz, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dsbevd_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer hbevd(char* jobz, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return chbevd_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer hbevd(char* jobz, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zhbevd_(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer sbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  return ssbevx_2stage_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline integer sbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  return dsbevx_2stage_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline integer hbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, scomplex* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  return chbevx_2stage_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}
inline integer hbevx_2stage(char* jobz, char* range, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, dcomplex* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  return zhbevx_2stage_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer sbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  return ssbevx_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline integer sbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  return dsbevx_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline integer hbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, scomplex* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  return chbevx_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}
inline integer hbevx(char* jobz, char* range, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, dcomplex* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  return zhbevx_(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}

// --- reduces a real symmetric-definite banded generalized eigenproblem  to standard form ---
inline integer sbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, float* ab, integer* ldab,  float* bb, integer* ldbb, float* x, integer* ldx, float* work, integer* info)
{
  return ssbgst_(vect, uplo, n, ka, kb, ab, ldab,  bb, ldbb, x, ldx, work, info);
}
inline integer sbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, double* ab, integer* ldab,  double* bb, integer* ldbb, double* x, integer* ldx, double* work, integer* info)
{
  return dsbgst_(vect, uplo, n, ka, kb, ab, ldab,  bb, ldbb, x, ldx, work, info);
}
inline integer hbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, scomplex* ab, integer* ldab,  scomplex* bb, integer* ldbb, scomplex* x, integer* ldx, scomplex* work, float* rwork, integer* info)
{
  return chbgst_(vect, uplo, n, ka, kb, ab, ldab,  bb, ldbb, x, ldx, work, rwork, info);
}
inline integer hbgst(char* vect, char* uplo, integer* n, integer* ka, integer* kb, dcomplex* ab, integer* ldab,  dcomplex* bb, integer* ldbb, dcomplex* x, integer* ldx, dcomplex* work, double* rwork, integer* info)
{
  return zhbgst_(vect, uplo, n, ka, kb, ab, ldab,  bb, ldbb, x, ldx, work, rwork, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors ---
inline integer sbgv(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, float* ab, integer* ldab, float* bb, integer* ldbb, float* w, float* z, integer* ldz, float* work, integer* info)
{
  return ssbgv_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, info);
}
inline integer sbgv(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, double* ab, integer* ldab, double* bb, integer* ldbb, double* w, double* z, integer* ldz, double* work, integer* info)
{
  return dsbgv_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, info);
}
inline integer hbgv(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, scomplex* ab, integer* ldab, scomplex* bb, integer* ldbb, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* info)
{
  return chbgv_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, rwork, info);
}
inline integer hbgv(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, dcomplex* ab, integer* ldab, dcomplex* bb, integer* ldbb, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* info)
{
  return zhbgv_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, rwork, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors ---
inline integer sbgvd(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, float* ab, integer* ldab, float* bb, integer* ldbb, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ssbgvd_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer sbgvd(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, double* ab, integer* ldab, double* bb, integer* ldbb, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dsbgvd_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer hbgvd(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, scomplex* ab, integer* ldab, scomplex* bb, integer* ldbb, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return chbgvd_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer hbgvd(char* jobz, char* uplo, integer* n, integer* ka, integer* kb, dcomplex* ab, integer* ldab, dcomplex* bb, integer* ldbb, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zhbgvd_(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes selected eigenvalues, and optionally, eigenvectors ---
inline integer sbgvx(char* jobz, char* range, char* uplo, integer* n, integer* ka, integer* kb, float* ab, integer* ldab, float* bb, integer* ldbb, float* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  return ssbgvx_(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline integer sbgvx(char* jobz, char* range, char* uplo, integer* n, integer* ka, integer* kb, double* ab, integer* ldab, double* bb, integer* ldbb, double* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  return dsbgvx_(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline integer hbgvx(char* jobz, char* range, char* uplo, integer* n, integer* ka, integer* kb, scomplex* ab, integer* ldab, scomplex* bb, integer* ldbb, scomplex* q, integer* ldq, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  return chbgvx_(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}
inline integer hbgvx(char* jobz, char* range, char* uplo, integer* n, integer* ka, integer* kb, dcomplex* ab, integer* ldab, dcomplex* bb, integer* ldbb, dcomplex* q, integer* ldq, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  return zhbgvx_(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}

// --- reduces a real symmetric band matrix A to symmetric tridiagonal form T ---
inline integer sbtrd(char* vect, char* uplo, integer* n, integer* kd, float* ab, integer* ldab, float* d, float* e, float* q, integer* ldq, float* work, integer* info)
{
  return ssbtrd_(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
}
inline integer sbtrd(char* vect, char* uplo, integer* n, integer* kd, double* ab, integer* ldab, double* d, double* e, double* q, integer* ldq, double* work, integer* info)
{
  return dsbtrd_(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
}
inline integer hbtrd(char* vect, char* uplo, integer* n, integer* kd, scomplex* ab, integer* ldab, float* d, float* e, scomplex* q, integer* ldq, scomplex* work, integer* info)
{
  return chbtrd_(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
}
inline integer hbtrd(char* vect, char* uplo, integer* n, integer* kd, dcomplex* ab, integer* ldab, double* d, double* e, dcomplex* q, integer* ldq, dcomplex* work, integer* info)
{
  return zhbtrd_(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info);
}

// --- performs a symmetric rank-k operation for matrix in RFP format ---
inline integer sfrk(char* transr, char* uplo, char* trans, integer* n, integer* k, float* alpha, float* a, integer* lda, float* beta, float* c)
{
  return ssfrk_(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}
inline integer sfrk(char* transr, char* uplo, char* trans, integer* n, integer* k, double* alpha,  double* a, integer* lda, double* beta, double* c)
{
  return dsfrk_(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}
inline integer hfrk(char* transr, char* uplo, char* trans, integer* n, integer* k, float* alpha,  scomplex* a, integer* lda, float* beta, scomplex* c)
{
  return chfrk_(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}
inline integer hfrk(char* transr, char* uplo, char* trans, integer* n, integer* k, double* alpha,  dcomplex* a, integer* lda, double* beta, dcomplex* c)
{
  return zhfrk_(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
}

// --- estimates the reciprocal of the condition number ---
inline integer spcon(char* uplo, integer* n, float* ap, integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  return sspcon_(uplo, n, ap, ipiv, anorm, rcond, work, iwork, info);
}
inline integer spcon(char* uplo, integer* n, double* ap, integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  return dspcon_(uplo, n, ap, ipiv, anorm, rcond, work, iwork, info);
}
inline integer spcon(char* uplo, integer* n, scomplex* ap, integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  return cspcon_(uplo, n, ap, ipiv, anorm, rcond, work, info);
}
inline integer spcon(char* uplo, integer* n,  dcomplex* ap, integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  return zspcon_(uplo, n, ap, ipiv, anorm, rcond, work, info);
}
inline integer hpcon(char* uplo, integer* n, scomplex* ap, integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  return chpcon_(uplo, n, ap, ipiv, anorm, rcond, work, info);
}
inline integer hpcon(char* uplo, integer* n, dcomplex* ap, integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  return zhpcon_(uplo, n, ap, ipiv, anorm, rcond, work, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer spev(char *jobz, char *uplo, integer* n, float* ap, float* w, float* z, integer*ldq, float* work, integer* info)
{
  return sspev_(jobz, uplo, n, ap, w, z, ldq, work, info);
}
inline integer spev(char* jobz, char* uplo, integer* n, double* ap, double* w, double* z, integer* ldq, double* work, integer* info)
{
  return dspev_(jobz, uplo, n, ap, w, z, ldq, work, info);
}
inline integer hpev(char* jobz, char* uplo, integer* n, scomplex* ap, float* w, scomplex* z, integer* ldq, scomplex* work, float* rwork, integer* info)
{
  return chpev_(jobz, uplo, n, ap, w, z, ldq, work, rwork, info);
}
inline integer hpev(char* jobz, char* uplo, integer* n, dcomplex* ap, double* w, dcomplex* z, integer* ldq, dcomplex* work, double* rwork, integer* info)
{
  return zhpev_(jobz, uplo, n, ap, w, z, ldq, work, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer spevd(char* jobz, char* uplo, integer* n, float* ap, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return sspevd_(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer spevd(char* jobz, char* uplo, integer* n, double* ap, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dspevd_(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer hpevd(char* jobz, char* uplo, integer* n, scomplex* ap, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return chpevd_(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer hpevd(char* jobz, char* uplo, integer* n, dcomplex* ap, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zhpevd_(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer spevx(char* jobz, char* range, char* uplo, integer* n, float* ap, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  return sspevx_(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline integer spevx(char* jobz, char* range, char* uplo, integer* n, double* ap, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  return dspevx_(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline integer hpevx(char* jobz, char* range, char* uplo, integer* n, scomplex* ap, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  return chpevx_(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}
inline integer hpevx(char* jobz, char* range, char* uplo, integer* n, dcomplex* ap, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  return zhpevx_(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}

// --- reduces a real symmetric-definite generalized eigenproblem to standard form ---
inline integer spgst(integer* itype, char* uplo, integer* n, float* ap,  float* bp, integer* info)
{
  return sspgst_(itype, uplo, n, ap, bp, info);
}
inline integer spgst(integer* itype, char* uplo, integer* n, double* ap,  double* bp, integer* info)
{
  return dspgst_(itype, uplo, n, ap, bp, info);
}
inline integer hpgst(integer* itype, char* uplo, integer* n, scomplex* ap,  scomplex* bp, integer* info)
{
  return chpgst_(itype, uplo, n, ap, bp, info);
}
inline integer hpgst(integer* itype, char* uplo, integer* n, dcomplex* ap,  dcomplex* bp, integer* info)
{
  return zhpgst_(itype, uplo, n, ap, bp, info);
}

// --- computes all the eigenvalues ---
inline integer spgv(integer* itype, char* jobz, char* uplo, integer* n, float* ap, float* bp, float* w, float* z, integer* ldz, float* work, integer* info)
{
  return sspgv_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info);
}
inline integer spgv(integer* itype, char* jobz, char* uplo, integer* n, double* ap, double* bp, double* w, double* z, integer* ldz, double* work, integer* info)
{
  return dspgv_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info);
}
inline integer hpgv(integer* itype, char* jobz, char* uplo, integer* n, scomplex* ap, scomplex* bp, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* info)
{
  return chpgv_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, info);
}
inline integer hpgv(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* ap, dcomplex* bp, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* info)
{
  return zhpgv_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, info);
}

// --- computes all the eigenvalues ---
inline integer spgvd(integer* itype, char* jobz, char* uplo, integer* n, float* ap, float* bp, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return sspgvd_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer spgvd(integer* itype, char* jobz, char* uplo, integer* n, double* ap, double* bp, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dspgvd_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer hpgvd(integer* itype, char* jobz, char* uplo, integer* n, scomplex* ap, scomplex* bp, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return chpgvd_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer hpgvd(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* ap, dcomplex* bp, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zhpgvd_(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes all the eigenvalues ---
inline integer spgvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, float* ap, float* bp, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  return sspgvx_(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline integer spgvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, double* ap, double* bp, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  return dspgvx_(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline integer hpgvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, scomplex* ap, scomplex* bp, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  return chpgvx_(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}
inline integer hpgvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, dcomplex* ap, dcomplex* bp, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  return zhpgvx_(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline integer spsv(char* uplo, integer* n, integer* nrhs, float* ap, integer* ipiv, float* b, integer* ldb, integer* info)
{
  return sspsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline integer spsv(char* uplo, integer* n, integer* nrhs, double* ap, integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dspsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline integer spsv(char* uplo, integer* n, integer* nrhs, scomplex* ap, integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return cspsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline integer spsv(char* uplo, integer* n, integer* nrhs, dcomplex* ap, integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zspsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline integer hpsv(char* uplo, integer* n, integer* nrhs, scomplex* ap, integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return chpsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}
inline integer hpsv(char* uplo, integer* n, integer* nrhs, dcomplex* ap, integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zhpsv_(uplo, n, nrhs, ap, ipiv, b, ldb, info);
}

// --- computes the solution to system of linear equations A * X = B for OTHER matrices ---
inline integer spsvx(char* fact, char* uplo, integer* n, integer* nrhs,  float* ap, float* afp, integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return sspsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer spsvx(char* fact, char* uplo, integer* n, integer* nrhs,  double* ap, double* afp, integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dspsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info);
}
inline integer spsvx(char* fact, char* uplo, integer* n, integer* nrhs,  scomplex* ap, scomplex* afp, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cspsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer spsvx(char* fact, char* uplo, integer* n, integer* nrhs,  dcomplex* ap, dcomplex* afp, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zspsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer hpsvx(char* fact, char* uplo, integer* n, integer* nrhs,  scomplex* ap, scomplex* afp, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return chpsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}
inline integer hpsvx(char* fact, char* uplo, integer* n, integer* nrhs,  dcomplex* ap, dcomplex* afp, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zhpsvx_(fact, uplo, n, nrhs,  ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
}

// --- reduces a real symmetric matrix A stored in packed form to symmetric tridiagonal form T  ---
inline integer sptrd(char* uplo, integer* n, float* ap, float* d, float* e, float* tau, integer* info)
{
  return ssptrd_(uplo, n, ap, d, e, tau, info);
}
inline integer sptrd(char* uplo, integer* n, double* ap, double* d, double* e, double* tau, integer* info)
{
  return dsptrd_(uplo, n, ap, d, e, tau, info);
}
inline integer hptrd(char* uplo, integer* n, scomplex* ap, float* d, float* e, scomplex* tau, integer* info)
{
  return chptrd_(uplo, n, ap, d, e, tau, info);
}
inline integer hptrd(char* uplo, integer* n, dcomplex* ap, double* d, double* e, dcomplex* tau, integer* info)
{
  return zhptrd_(uplo, n, ap, d, e, tau, info);
}

// --- computes the factorization of a real symmetric matrix  ---
inline integer sptrf(char* uplo, integer* n, float* ap, integer* ipiv, integer* info)
{
  return ssptrf_(uplo, n, ap, ipiv, info);
}
inline integer sptrf(char* uplo, integer* n, double* ap, integer* ipiv, integer* info)
{
  return dsptrf_(uplo, n, ap, ipiv, info);
}
inline integer sptrf(char* uplo, integer* n, scomplex* ap, integer* ipiv, integer* info)
{
  return csptrf_(uplo, n, ap, ipiv, info);
}
inline integer sptrf(char* uplo, integer* n, dcomplex* ap, integer* ipiv, integer* info)
{
  return zsptrf_(uplo, n, ap, ipiv, info);
}
inline integer hptrf(char* uplo, integer* n, scomplex* ap, integer* ipiv, integer* info)
{
  return chptrf_(uplo, n, ap, ipiv, info);
}
inline integer hptrf(char* uplo, integer* n, dcomplex* ap, integer* ipiv, integer* info)
{
  return zhptrf_(uplo, n, ap, ipiv, info);
}

// --- computes the inverse of a real symmetric indefinite matrix ---
inline integer sptri(char* uplo, integer* n, float* ap,  integer* ipiv, float* work, integer* info)
{
  return ssptri_(uplo, n, ap, ipiv, work, info);
}
inline integer sptri(char* uplo, integer* n, double* ap,  integer* ipiv, double* work, integer* info)
{
  return dsptri_(uplo, n, ap, ipiv, work, info);
}
inline integer sptri(char* uplo, integer* n, scomplex* ap,  integer* ipiv, scomplex* work, integer* info)
{
  return csptri_(uplo, n, ap, ipiv, work, info);
}
inline integer sptri(char* uplo, integer* n, dcomplex* ap,  integer* ipiv, dcomplex* work, integer* info)
{
  return zsptri_(uplo, n, ap, ipiv, work, info);
}
inline integer hptri(char* uplo, integer* n, scomplex* ap, integer* ipiv, scomplex* work, integer* info)
{
  return chptri_(uplo, n, ap, ipiv, work, info);
}
inline integer hptri(char* uplo, integer* n, dcomplex* ap, integer* ipiv, dcomplex* work, integer* info)
{
  return zhptri_(uplo, n, ap, ipiv, work, info);
}

// --- solves a system of linear equations A*X = B ---
inline integer sptrs(char* uplo, integer* n, integer* nrhs,  float* ap,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  return ssptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline integer sptrs(char* uplo, integer* n, integer* nrhs,  double* ap,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dsptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline integer sptrs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return csptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline integer sptrs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zsptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline integer hptrs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return chptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}
inline integer hptrs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zhptrs_(uplo, n, nrhs,  ap, ipiv, b, ldb, info);
}

// --- computes the eigenvalues of a symmetric tridiagonal matrix T ---
inline integer stebz(char* range, char* order, integer* n, float* vl, float* vu, integer* il, integer* iu, float* abstol, float* d, float* e, integer* m, integer* nsplit, float* w, integer* iblock, integer* isplit, float* work, integer* iwork, integer* info)
{
  return sstebz_(range, order, n, vl, vu, il, iu, abstol, d,  e, m,  nsplit, w, iblock, isplit, work, iwork, info);
}
inline integer stebz(char* range, char* order, integer* n, double* vl, double* vu, integer* il, integer* iu, double* abstol, double* d, double* e, integer* m, integer* nsplit, double* w, integer* iblock, integer* isplit, double* work, integer* iwork, integer* info)
{
  return dstebz_(range, order, n, vl, vu, il, iu, abstol, d,  e, m,  nsplit, w, iblock, isplit, work, iwork, info);
}

// --- computes selected eigenvalues and, optionally, eigenvectors of a real symmetric tridiagonal matrix T---
inline integer stegr(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return sstegr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline integer stegr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dstegr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline integer stegr(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return cstegr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline integer stegr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return zstegr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvectors of a real symmetric tridiagonal matrix T ---
inline integer stein(integer* n,  float* d,  float* e, integer* m,  float* w,  integer* iblock,  integer* isplit, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  return sstein_(n,  d,  e, m,  w,  iblock,  isplit, z, ldz, work, iwork, ifail, info);
}
inline integer stein(integer* n,  double* d,  double* e, integer* m,  double* w,  integer* iblock,  integer* isplit, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  return dstein_(n,  d,  e, m,  w,  iblock,  isplit, z, ldz, work, iwork, ifail, info);
}
inline integer stein(integer* n,  float* d,  float* e, integer* m,  float* w,  integer* iblock,  integer* isplit, scomplex* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  return cstein_(n,  d,  e, m,  w,  iblock,  isplit, z, ldz, work, iwork, ifail, info);
}
inline integer stein(integer* n,  double* d,  double* e, integer* m,  double* w,  integer* iblock,  integer* isplit, dcomplex* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  return zstein_(n,  d,  e, m,  w,  iblock,  isplit, z, ldz, work, iwork, ifail, info);
}

// --- STERF computes all eigenvalues of a symmetric tridiagonal matrix ---
inline integer sterf(integer* n, float* d, float* e, integer* info)
{
  return ssterf_(n, d, e, info);
}
inline integer sterf(integer* n, double* d, double* e, integer* info)
{
  return dsterf_(n, d, e, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer stev(char* jobz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* info)
{
  return sstev_(jobz, n, d, e, z, ldz, work, info);
}
inline integer stev(char* jobz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* info)
{
  return dstev_(jobz, n, d, e, z, ldz, work, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer stevd(char* jobz, integer* n, float* d, float* e, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return sstevd_(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
}
inline integer stevd(char* jobz, integer* n, double* d, double* e, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dstevd_(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer stevr(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return sstevr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline integer stevr(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dstevr_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices ---
inline integer stevx(char* jobz, char* range, integer* n, float* d, float* e, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* iwork, integer* ifail, integer* info)
{
  return sstevx_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}
inline integer stevx(char* jobz, char* range, integer* n, double* d, double* e, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* iwork, integer* ifail, integer* info)
{
  return dstevx_(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
}

// --- estimates the reciprocal of the condition number of a real symmetric matrix A ---
inline integer sycon_3(char* uplo, integer* n,  float* a, integer* lda,  float* e,  integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  return ssycon_3_(uplo, n,  a, lda,  e, ipiv, anorm, rcond, work, iwork, info);
}
inline integer sycon_3(char* uplo, integer* n,  double* a, integer* lda,  double* e,  integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  return dsycon_3_(uplo, n,  a, lda,  e, ipiv, anorm, rcond, work, iwork, info);
}
inline integer sycon_3(char* uplo, integer* n,  scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  return csycon_3_(uplo, n,  a, lda,  e, ipiv, anorm, rcond, work, info);
}
inline integer sycon_3(char* uplo, integer* n,  dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  return zsycon_3_(uplo, n,  a, lda,  e, ipiv, anorm, rcond, work, info);
}

// --- estimates the reciprocal of the condition number of a real symmetric matrix A ---
inline integer sycon(char* uplo, integer* n,  float* a, integer* lda,  integer* ipiv, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  return ssycon_(uplo, n,  a, lda,  ipiv, anorm, rcond, work, iwork, info);
}
inline integer sycon(char* uplo, integer* n,  double* a, integer* lda,  integer* ipiv, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  return dsycon_(uplo, n,  a, lda,  ipiv, anorm, rcond, work, iwork, info);
}
inline integer sycon(char* uplo, integer* n,  scomplex* a, integer* lda,  integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  return csycon_(uplo, n,  a, lda,  ipiv, anorm, rcond, work, info);
}
inline integer sycon(char* uplo, integer* n,  dcomplex* a, integer* lda,  integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  return zsycon_(uplo, n,  a, lda,  ipiv, anorm, rcond, work, info);
}

// --- convert A given by TRF into L and D and vice-versa ---
inline integer syconv(char* uplo, char* way, integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* info)
{
  return ssyconv_(uplo, way, n, a, lda, ipiv, work, info);
}
inline integer syconv(char* uplo, char* way, integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* info)
{
  return dsyconv_(uplo, way, n, a, lda, ipiv, work, info);
}
inline integer syconv(char* uplo, char* way, integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* info)
{
  return csyconv_(uplo, way, n, a, lda, ipiv, work, info);
}
inline integer syconv(char* uplo, char* way, integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* info)
{
  return zsyconv_(uplo, way, n, a, lda, ipiv, work, info);
}

// --- computes row and column scalings intended to equilibrate a symmetric matrix A ---
inline integer syequb(char* uplo, integer* n,  float* a, integer* lda, float* s, float* scond, float* amax, float* work, integer* info)
{
  return ssyequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}
inline integer syequb(char* uplo, integer* n,  double* a, integer* lda, double* s, double* scond, double* amax, double* work, integer* info)
{
  return dsyequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}
inline integer syequb(char* uplo, integer* n,  scomplex* a, integer* lda, float* s, float* scond, float* amax, scomplex* work, integer* info)
{
  return csyequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}
inline integer syequb(char* uplo, integer* n,  dcomplex* a, integer* lda, double* s, double* scond, double* amax, dcomplex* work, integer* info)
{
  return zsyequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices ---
inline integer syev_2stage(char* jobz, char* uplo, integer* n, float* a, integer* lda, float* w, float* work, integer* lwork, integer* info)
{
  return ssyev_2stage_(jobz, uplo, n, a, lda, w, work, lwork, info);
}
inline integer syev_2stage(char* jobz, char* uplo, integer* n, double* a, integer* lda, double* w, double* work, integer* lwork, integer* info)
{
  return dsyev_2stage_(jobz, uplo, n, a, lda, w, work, lwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices ---
inline integer syevd_2stage(char* jobz, char* uplo, integer* n, float* a, integer* lda, float* w, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ssyevd_2stage_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
}
inline integer syevd_2stage(char* jobz, char* uplo, integer* n, double* a, integer* lda, double* w, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dsyevd_2stage_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices --
inline integer syevr_2stage(char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, integer* isuppz, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ssyevr_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
inline integer syevr_2stage(char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, integer* isuppz, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dsyevr_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices ---
inline integer syevx_2stage(char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  return ssyevx_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline integer syevx_2stage(char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  return dsyevx_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices ---
inline integer syevx(char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  return ssyevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline integer syevx(char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  return dsyevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}

// --- computes all the eigenvalues, the eigenvectors of a real generalized symmetric-definite eigenproblem ---
inline integer sygv_2stage(integer* itype, char* jobz, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, float* w, float* work, integer* lwork, integer* info)
{
  return ssygv_2stage_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
}
inline integer sygv_2stage(integer* itype, char* jobz, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, double* w, double* work, integer* lwork, integer* info)
{
  return dsygv_2stage_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
}

// --- computes all the eigenvalues, the eigenvectors of a real generalized symmetric-definite eigenproblem ---
inline integer sygv(integer* itype, char* jobz, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, float* w, float* work, integer* lwork, integer* info)
{
  return ssygv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
}
inline integer sygv(integer* itype, char* jobz, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, double* w, double* work, integer* lwork, integer* info)
{
  return dsygv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
}

// --- computes all the eigenvalues, the eigenvectors of a real generalized symmetric-definite eigenproblem ---
inline integer sygvd(integer* itype, char* jobz, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, float* w, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ssygvd_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info);
}
inline integer sygvd(integer* itype, char* jobz, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, double* w, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dsygvd_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info);
}

// --- computes selected eigenvalues, eigenvectors of a real generalized symmetric-definite eigenproblem ---
inline integer sygvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, float* a, integer* lda, float* b, integer* ldb, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, float* z, integer* ldz, float* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  return ssygvx_(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}
inline integer sygvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, double* a, integer* lda, double* b, integer* ldb, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, double* z, integer* ldz, double* work, integer* lwork, integer* iwork, integer* ifail, integer* info)
{
  return dsygvx_(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info);
}

// --- improves the computed solution to a system of linear equations---
inline integer syrfs(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return ssyrfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer syrfs(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dsyrfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer syrfs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return csyrfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer syrfs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zsyrfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- improves the computed solution to a system of linear equation---
inline integer syrfsx(char* uplo, char* equed, integer* n, integer* nrhs,  float* a, integer* lda,  float* af, integer* ldaf,  integer* ipiv,  float* s,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  return ssyrfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer syrfsx(char* uplo, char* equed, integer* n, integer* nrhs,  double* a, integer* lda,  double* af, integer* ldaf,  integer* ipiv,  double* s,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  return dsyrfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer syrfsx(char* uplo, char* equed, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  float* s,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return csyrfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer syrfsx(char* uplo, char* equed, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  double* s,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zsyrfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline integer sysv_aa_2stage(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* tb, integer* ltb, integer* ipiv, integer* ipiv2, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  return ssysv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}
inline integer sysv_aa_2stage(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* tb, integer* ltb, integer* ipiv, integer* ipiv2, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  return dsysv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}
inline integer sysv_aa_2stage(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return csysv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}
inline integer sysv_aa_2stage(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zsysv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline integer sysv_aa(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  return ssysv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv_aa(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  return dsysv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv_aa(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return csysv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv_aa(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zsysv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline integer sysv_rk(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* e, integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  return ssysv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv_rk(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* e, integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  return dsysv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv_rk(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* e, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return csysv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv_rk(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* e, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zsysv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline integer sysv_rook(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  return ssysv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv_rook(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  return dsysv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv_rook(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return csysv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv_rook(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zsysv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline integer sysv(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  return ssysv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  return dsysv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return csysv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer sysv(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zsysv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline integer sysvx(char* fact, char* uplo, integer* n, integer* nrhs,  float* a, integer* lda, float* af, integer* ldaf, integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* ferr, float* berr, float* work, integer* lwork, integer* iwork, integer* info)
{
  return ssysvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, iwork, info);
}
inline integer sysvx(char* fact, char* uplo, integer* n, integer* nrhs,  double* a, integer* lda, double* af, integer* ldaf, integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* ferr, double* berr, double* work, integer* lwork, integer* iwork, integer* info)
{
  return dsysvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, iwork, info);
}
inline integer sysvx(char* fact, char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return csysvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info);
}
inline integer sysvx(char* fact, char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zsysvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info);
}

// --- uses the diagonal pivoting factorization to compute the solution to a real system of linear equations ---
inline integer sysvxx(char* fact, char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* af, integer* ldaf, integer* ipiv, char* equed, float* s, float* b, integer* ldb, float* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, float* work, integer* iwork, integer* info)
{
  return ssysvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer sysvxx(char* fact, char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* af, integer* ldaf, integer* ipiv, char* equed, double* s, double* b, integer* ldb, double* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, double* work, integer* iwork, integer* info)
{
  return dsysvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
}
inline integer sysvxx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return csysvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer sysvxx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zsysvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- applies an elementary permutation on the rows and columns of a symmetric matrix ---
inline integer syswapr(char* uplo, integer* n, float* a, integer* lda, integer* i1, integer* i2)
{
  return ssyswapr_(uplo, n, a, lda, i1, i2);
}
inline integer syswapr(char* uplo, integer* n, double* a, integer* lda, integer* i1, integer* i2)
{
  return dsyswapr_(uplo, n, a, lda, i1, i2);
}
inline integer syswapr(char* uplo, integer* n, scomplex* a, integer* lda, integer* i1, integer* i2)
{
  return csyswapr_(uplo, n, a, lda, i1, i2);
}
inline integer syswapr(char* uplo, integer* n, dcomplex* a, integer* lda, integer* i1, integer* i2)
{
  return zsyswapr_(uplo, n, a, lda, i1, i2);
}

// --- computes the factorization of a real symmetric matrix A using the Aasen's algorithm---
inline integer sytrf_aa_2stage(char* uplo, integer* n, float* a, integer* lda, float* tb, integer* ltb, integer* ipiv, integer* ipiv2, float* work, integer* lwork, integer* info)
{
  return ssytrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}
inline integer sytrf_aa_2stage(char* uplo, integer* n, double* a, integer* lda, double* tb, integer* ltb, integer* ipiv, integer* ipiv2, double* work, integer* lwork, integer* info)
{
  return dsytrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}
inline integer sytrf_aa_2stage(char* uplo, integer* n, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* work, integer* lwork, integer* info)
{
  return csytrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}
inline integer sytrf_aa_2stage(char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* work, integer* lwork, integer* info)
{
  return zsytrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}

// --- computes the factorization of a real symmetric matrix A using the Aasen's algorithm ---
inline integer sytrf_aa(char* uplo, integer* n, float* a, integer* lda, integer* ipiv, float* work, integer* lwork, integer* info)
{
  return ssytrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytrf_aa(char* uplo, integer* n, double* a, integer* lda, integer* ipiv, double* work, integer* lwork, integer* info)
{
  return dsytrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytrf_aa(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return csytrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytrf_aa(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zsytrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the factorization of a real symmetric matrix A using the bounded Bunch-Kaufman (rook) diagonal pivoting method ---
inline integer sytrf_rk(char* uplo, integer* n, float* a, integer* lda, float* e, integer* ipiv, float* work, integer* lwork, integer* info)
{
  return ssytrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}
inline integer sytrf_rk(char* uplo, integer* n, double* a, integer* lda, double* e, integer* ipiv, double* work, integer* lwork, integer* info)
{
  return dsytrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}
inline integer sytrf_rk(char* uplo, integer* n, scomplex* a, integer* lda, scomplex* e, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return csytrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}
inline integer sytrf_rk(char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* e, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zsytrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}

// --- computes the factorization of a real symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting method---
inline integer sytrf_rook(char* uplo, integer* n, float* a, integer* lda, integer* ipiv, float* work, integer* lwork, integer* info)
{
  return ssytrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytrf_rook(char* uplo, integer* n, double* a, integer* lda, integer* ipiv, double* work, integer* lwork, integer* info)
{
  return dsytrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytrf_rook(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return csytrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytrf_rook(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zsytrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the factorization of a real symmetric matrix A using the Bunch-Kaufman diagonal pivoting method ---
inline integer sytrf(char* uplo, integer* n, float* a, integer* lda, integer* ipiv, float* work, integer* lwork, integer* info)
{
  return ssytrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytrf(char* uplo, integer* n, double* a, integer* lda, integer* ipiv, double* work, integer* lwork, integer* info)
{
  return dsytrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytrf(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return csytrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytrf(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zsytrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- SYTRI_3 computes the inverse of a real symmetric indefinite matrix A using the factorization computed by SYTRF_RK or SYTRF_BK---
inline integer sytri_3(char* uplo, integer* n, float* a, integer* lda,  float* e,  integer* ipiv, float* work, integer* lwork, integer* info)
{
  return ssytri_3_(uplo, n, a, lda,  e, ipiv, work, lwork, info);
}
inline integer sytri_3(char* uplo, integer* n, double* a, integer* lda,  double* e,  integer* ipiv, double* work, integer* lwork, integer* info)
{
  return dsytri_3_(uplo, n, a, lda,  e, ipiv, work, lwork, info);
}
inline integer sytri_3(char* uplo, integer* n, scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return csytri_3_(uplo, n, a, lda,  e, ipiv, work, lwork, info);
}
inline integer sytri_3(char* uplo, integer* n, dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zsytri_3_(uplo, n, a, lda,  e, ipiv, work, lwork, info);
}

// --- computes the inverse of a real symmetric indefinite matrix A---
inline integer sytri(char* uplo, integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* info)
{
  return ssytri_(uplo, n, a, lda, ipiv, work, info);
}
inline integer sytri(char* uplo, integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* info)
{
  return dsytri_(uplo, n, a, lda, ipiv, work, info);
}
inline integer sytri(char* uplo, integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* info)
{
  return csytri_(uplo, n, a, lda, ipiv, work, info);
}
inline integer sytri(char* uplo, integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* info)
{
  return zsytri_(uplo, n, a, lda, ipiv, work, info);
}

// --- computes the inverse of a REAL symmetric indefinite matrix A using the factorization ---
inline integer sytri2(char* uplo, integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* lwork, integer* info)
{
  return ssytri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytri2(char* uplo, integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* lwork, integer* info)
{
  return dsytri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytri2(char* uplo, integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return csytri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer sytri2(char* uplo, integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zsytri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the inverse of a real symmetric indefinite matrix A using the factorization ---
inline integer sytri2x(char* uplo, integer* n, float* a, integer* lda,  integer* ipiv, float* work, integer* nb, integer* info)
{
  return ssytri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}
inline integer sytri2x(char* uplo, integer* n, double* a, integer* lda,  integer* ipiv, double* work, integer* nb, integer* info)
{
  return dsytri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}
inline integer sytri2x(char* uplo, integer* n, scomplex* a, integer* lda,  integer* ipiv, scomplex* work, integer* nb, integer* info)
{
  return csytri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}
inline integer sytri2x(char* uplo, integer* n, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* work, integer* nb, integer* info)
{
  return zsytri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}

// --- solves a system of linear equations A * X = B with a real symmetric matrix A ---
inline integer sytrs_3(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  float* e,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  return ssytrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv,  b, ldb, info);
}
inline integer sytrs_3(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  double* e,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dsytrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv,  b, ldb, info);
}
inline integer sytrs_3(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return csytrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv,  b, ldb, info);
}
inline integer sytrs_3(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zsytrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv,  b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A---
inline integer sytrs_aa_2stage(char* uplo, integer* n, integer* nrhs, float* a, integer* lda, float* tb, integer* ltb, integer* ipiv, integer* ipiv2, float* b, integer* ldb, integer* info)
{
  return ssytrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}
inline integer sytrs_aa_2stage(char* uplo, integer* n, integer* nrhs, double* a, integer* lda, double* tb, integer* ltb, integer* ipiv, integer* ipiv2, double* b, integer* ldb, integer* info)
{
  return dsytrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}
inline integer sytrs_aa_2stage(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* b, integer* ldb, integer* info)
{
  return csytrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}
inline integer sytrs_aa_2stage(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* b, integer* ldb, integer* info)
{
  return zsytrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A ---
inline integer sytrs_aa(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, float* work, integer* lwork, integer* info)
{
  return ssytrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}
inline integer sytrs_aa(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, double* work, integer* lwork, integer* info)
{
  return dsytrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}
inline integer sytrs_aa(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return csytrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}
inline integer sytrs_aa(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zsytrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A---
inline integer sytrs_rook(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  return ssytrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer sytrs_rook(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dsytrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer sytrs_rook(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return csytrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer sytrs_rook(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zsytrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A ---
inline integer sytrs(char* uplo, integer* n, integer* nrhs,  float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, integer* info)
{
  return ssytrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer sytrs(char* uplo, integer* n, integer* nrhs,  double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, integer* info)
{
  return dsytrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer sytrs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return csytrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer sytrs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zsytrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a real symmetric matrix A ---
inline integer sytrs2(char* uplo, integer* n, integer* nrhs, float* a, integer* lda,  integer* ipiv, float* b, integer* ldb, float* work, integer* info)
{
  return ssytrs2_(uplo, n, nrhs, a, lda,  ipiv, b, ldb, work, info);
}
inline integer sytrs2(char* uplo, integer* n, integer* nrhs, double* a, integer* lda,  integer* ipiv, double* b, integer* ldb, double* work, integer* info)
{
  return dsytrs2_(uplo, n, nrhs, a, lda,  ipiv, b, ldb, work, info);
}
inline integer sytrs2(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* info)
{
  return csytrs2_(uplo, n, nrhs, a, lda,  ipiv, b, ldb, work, info);
}
inline integer sytrs2(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* info)
{
  return zsytrs2_(uplo, n, nrhs, a, lda,  ipiv, b, ldb, work, info);
}

// --- computes the eigenvalues of a real matrix pair (H,T)---
inline integer hgeqz(char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, float* h, integer* ldh, float* t, integer* ldt, float* alphar, float* alphai, float* beta, float* q, integer* ldq, float* z, integer* ldz, float* work, integer* lwork, integer* info)
{
  return shgeqz_(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, info);
}
inline integer hgeqz(char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, double* h, integer* ldh, double* t, integer* ldt, double* alphar, double* alphai, double* beta, double* q, integer* ldq, double* z, integer* ldz, double* work, integer* lwork, integer* info)
{
  return dhgeqz_(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z, ldz, work, lwork, info);
}
inline integer hgeqz(char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, scomplex* h, integer* ldh, scomplex* t, integer* ldt, scomplex* alpha, scomplex* beta, scomplex* q, integer* ldq, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return chgeqz_(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, info);
}
inline integer hgeqz(char* job, char* compq, char* compz, integer* n, integer* ilo, integer* ihi, dcomplex* h, integer* ldh, dcomplex* t, integer* ldt, dcomplex* alpha, dcomplex* beta, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zhgeqz_(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, info);
}

// --- uses inverse iteration to find specified right and/or left eigenvectors of a real upper Hessenberg matrix H ---
inline integer hsein(char* job, char* eigsrc, char* initv, logical* select, integer* n,  float* h, integer* ldh, float* wr,  float* wi, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* mm, integer* m, float* work, integer* ifaill, integer* ifailr, integer* info)
{
  return shsein_(job, eigsrc, initv, select, n,  h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, work, ifaill, ifailr, info);
}
inline integer hsein(char* job, char* eigsrc, char* initv, logical* select, integer* n,  double* h, integer* ldh, double* wr,  double* wi, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* mm, integer* m, double* work, integer* ifaill, integer* ifailr, integer* info)
{
  return dhsein_(job, eigsrc, initv, select, n,  h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, work, ifaill, ifailr, info);
}
inline integer hsein(char* job, char* eigsrc, char* initv, logical* select, integer* n,  scomplex* h, integer* ldh, scomplex* w, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* mm, integer* m, scomplex* work, float* rwork, integer* ifaill, integer* ifailr, integer* info)
{
  return chsein_(job, eigsrc, initv, select, n,  h, ldh, w,  vl, ldvl, vr, ldvr, mm, m, work, rwork, ifaill, ifailr, info);
}
inline integer hsein(char* job, char* eigsrc, char* initv, logical* select, integer* n,  dcomplex* h, integer* ldh, dcomplex* w, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* mm, integer* m, dcomplex* work, double* rwork, integer* ifaill, integer* ifailr, integer* info)
{
  return zhsein_(job, eigsrc, initv, select, n,  h, ldh, w, vl, ldvl, vr, ldvr, mm, m, work, rwork, ifaill, ifailr, info);
}

// --- computes the eigenvalues of a Hessenberg matrix H ---
inline integer hseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, float* h, integer* ldh, float* wr, float* wi, float *z, integer* ldz, float* work, integer* lwork, integer* info)
{
  return shseqr_(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);
}
inline integer hseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, double* h, integer* ldh, double* wr, double* wi, double *z, integer* ldz, double* work, integer* lwork, integer* info)
{
  return dhseqr_(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);
}
inline integer hseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, scomplex* h, integer* ldh, scomplex* w, scomplex*z, integer* ldz, scomplex* work, integer* lwork, integer* info)
{
  return chseqr_(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);
}
inline integer hseqr(char* job, char* compz, integer* n, integer* ilo, integer* ihi, dcomplex* h, integer* ldh, dcomplex* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, integer* info)
{
  return zhseqr_(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);
}

// --- TBCON estimates the reciprocal of the condition number of a triangular band matrix A---
inline integer tbcon(char* norm, char* uplo, char* diag, integer* n, integer* kd,  float* ab, integer* ldab, float* rcond, float* work, integer* iwork, integer* info)
{
  return stbcon_(norm, uplo, diag, n, kd,  ab, ldab, rcond, work, iwork, info);
}
inline integer tbcon(char* norm, char* uplo, char* diag, integer* n, integer* kd,  double* ab, integer* ldab, double* rcond, double* work, integer* iwork, integer* info)
{
  return dtbcon_(norm, uplo, diag, n, kd,  ab, ldab, rcond, work, iwork, info);
}
inline integer tbcon(char* norm, char* uplo, char* diag, integer* n, integer* kd,  scomplex* ab, integer* ldab, float* rcond, scomplex* work, float* rwork, integer* info)
{
  return ctbcon_(norm, uplo, diag, n, kd,  ab, ldab, rcond, work, rwork, info);
}
inline integer tbcon(char* norm, char* uplo, char* diag, integer* n, integer* kd,  dcomplex* ab, integer* ldab, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  return ztbcon_(norm, uplo, diag, n, kd,  ab, ldab, rcond, work, rwork, info);
}

// --- TBRFS provides error bounds and backward error estimates for the solution to a system of linear equations ---
inline integer tbrfs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  float* ab, integer* ldab,  float* b, integer* ldb,  float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return stbrfs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline integer tbrfs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  double* ab, integer* ldab,  double* b, integer* ldb,  double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dtbrfs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline integer tbrfs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  scomplex* ab, integer* ldab,  scomplex* b, integer* ldb,  scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return ctbrfs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}
inline integer tbrfs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  dcomplex* ab, integer* ldab,  dcomplex* b, integer* ldb,  dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return ztbrfs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}

// --- solves a triangular system of the form A * X = B  or  A**T * X = B ---
inline integer tbtrs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  float* ab, integer* ldab, float* b, integer* ldb, integer* info)
{
  return stbtrs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab, b, ldb, info);
}
inline integer tbtrs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  double* ab, integer* ldab, double* b, integer* ldb, integer* info)
{
  return dtbtrs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab, b, ldb, info);
}
inline integer tbtrs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  scomplex* ab, integer* ldab, scomplex* b, integer* ldb, integer* info)
{
  return ctbtrs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab, b, ldb, info);
}
inline integer tbtrs(char* uplo, char* trans, char* diag, integer* n, integer* kd, integer* nrhs,  dcomplex* ab, integer* ldab, dcomplex* b, integer* ldb, integer* info)
{
  return ztbtrs_(uplo, trans, diag, n, kd, nrhs,  ab, ldab, b, ldb, info);
}

// --- solves a matrix equation (one operand is a triangular matrix in RFP format). ---
inline integer tfsm(char* transr, char* side, char* uplo, char* trans, char* diag, integer* m, integer* n, float* alpha, float* a, float* b, integer* ldb)
{
  return stfsm_(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}
inline integer tfsm(char* transr, char* side, char* uplo, char* trans, char* diag, integer* m, integer* n, double* alpha, double* a, double* b, integer* ldb)
{
  return dtfsm_(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}
inline integer tfsm(char* transr, char* side, char* uplo, char* trans, char* diag, integer* m, integer* n, scomplex* alpha, scomplex* a, scomplex* b, integer* ldb)
{
  return ctfsm_(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}
inline integer tfsm(char* transr, char* side, char* uplo, char* trans, char* diag, integer* m, integer* n, dcomplex* alpha, dcomplex* a, dcomplex* b, integer* ldb)
{
  return ztfsm_(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
}

// --- computes the inverse of a triangular matrix A stored in RFP format ---
inline integer tftri(char* transr, char* uplo, char* diag, integer* n, float* a, integer* info)
{
  return stftri_(transr, uplo, diag, n, a, info);
}
inline integer tftri(char* transr, char* uplo, char* diag, integer* n, double* a, integer* info)
{
  return dtftri_(transr, uplo, diag, n, a, info);
}
inline integer tftri(char* transr, char* uplo, char* diag, integer* n, scomplex* a, integer* info)
{
  return ctftri_(transr, uplo, diag, n, a, info);
}
inline integer tftri(char* transr, char* uplo, char* diag, integer* n, dcomplex* a, integer* info)
{
  return ztftri_(transr, uplo, diag, n, a, info);
}

// --- copies a triangular matrix from the rectangular full packed format (TF) to the standard packed format (TP) ---
inline integer tfttp(char* transr, char* uplo, integer* n,  float* arf, float* ap, integer* info)
{
  return stfttp_(transr, uplo, n,  arf, ap, info);
}
inline integer tfttp(char* transr, char* uplo, integer* n,  double* arf, double* ap, integer* info)
{
  return dtfttp_(transr, uplo, n,  arf, ap, info);
}
inline integer tfttp(char* transr, char* uplo, integer* n,  scomplex* arf, scomplex* ap, integer* info)
{
  return ctfttp_(transr, uplo, n,  arf, ap, info);
}
inline integer tfttp(char* transr, char* uplo, integer* n,  dcomplex* arf, dcomplex* ap, integer* info)
{
  return ztfttp_(transr, uplo, n,  arf, ap, info);
}

// --- copies a triangular matrix from the rectangular full packed format (TF) to the standard full format (TR) ---
inline integer tfttr(char* transr, char* uplo, integer* n,  float* arf, float* a, integer* lda, integer* info)
{
  return stfttr_(transr, uplo, n,  arf, a, lda, info);
}
inline integer tfttr(char* transr, char* uplo, integer* n,  double* arf, double* a, integer* lda, integer* info)
{
  return dtfttr_(transr, uplo, n,  arf, a, lda, info);
}
inline integer tfttr(char* transr, char* uplo, integer* n,  scomplex* arf, scomplex* a, integer* lda, integer* info)
{
  return ctfttr_(transr, uplo, n,  arf, a, lda, info);
}
inline integer tfttr(char* transr, char* uplo, integer* n,  dcomplex* arf, dcomplex* a, integer* lda, integer* info)
{
  return ztfttr_(transr, uplo, n,  arf, a, lda, info);
}

// --- computes some or all of the right and/or left eigenvectors of a pair of real matrices ---
inline integer tgevc(char* side, char* howmny,  logical* select, integer* n,  float* s, integer* lds,  float* p, integer* ldp, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* mm, integer* m, float* work, integer* info)
{
  return stgevc_(side, howmny,  select, n,  s, lds,  p, ldp, vl, ldvl, vr, ldvr, mm, m, work, info);
}
inline integer tgevc(char* side, char* howmny,  logical* select, integer* n,  double* s, integer* lds,  double* p, integer* ldp, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* mm, integer* m, double* work, integer* info)
{
  return dtgevc_(side, howmny,  select, n,  s, lds,  p, ldp, vl, ldvl, vr, ldvr, mm, m, work, info);
}
inline integer tgevc(char* side, char* howmny,  logical* select, integer* n,  scomplex* s, integer* lds,  scomplex* p, integer* ldp, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* mm, integer* m, scomplex* work, float* rwork, integer* info)
{
  return ctgevc_(side, howmny,  select, n,  s, lds,  p, ldp, vl, ldvl, vr, ldvr, mm, m, work, rwork, info);
}
inline integer tgevc(char* side, char* howmny,  logical* select, integer* n,  dcomplex* s, integer* lds,  dcomplex* p, integer* ldp, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* mm, integer* m, dcomplex* work, double* rwork, integer* info)
{
  return ztgevc_(side, howmny,  select, n,  s, lds,  p, ldp, vl, ldvl, vr, ldvr, mm, m, work, rwork, info);
}

// --- reorders the generalized real Schur decomposition of a real matrix pair ---
inline integer tgexc(logical* wantq, logical* wantz, integer* n, float* a, integer* lda, float* b, integer* ldb, float* q, integer* ldq, float* z, integer* ldz, integer* ifst, integer* ilst, float* work, integer* lwork, integer* info)
{
  return stgexc_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork, info);
}
inline integer tgexc(logical* wantq, logical* wantz, integer* n, double* a, integer* lda, double* b, integer* ldb, double* q, integer* ldq, double* z, integer* ldz, integer* ifst, integer* ilst, double* work, integer* lwork, integer* info)
{
  return dtgexc_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork, info);
}
inline integer tgexc(logical* wantq, logical* wantz, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* q, integer* ldq, scomplex* z, integer* ldz, integer* ifst, integer* ilst, integer* info)
{
  return ctgexc_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, info);
}
inline integer tgexc(logical* wantq, logical* wantz, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, integer* ifst, integer* ilst, integer* info)
{
  return ztgexc_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, info);
}

// --- reorders the generalized real Schur decomposition of a real matrix pair ---
inline integer tgsen(integer* ijob, logical* wantq, logical* wantz,  logical* select, integer* n, float* a, integer* lda, float* b, integer* ldb, float* alphar, float* alphai, float* beta, float* q, integer* ldq, float* z, integer* ldz, integer* m, float* pl, float* pr, float* dif, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return stgsen_(ijob, wantq, wantz, select, n, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info);
}
inline integer tgsen(integer* ijob, logical* wantq, logical* wantz,  logical* select, integer* n, double* a, integer* lda, double* b, integer* ldb, double* alphar, double* alphai, double* beta, double* q, integer* ldq, double* z, integer* ldz, integer* m, double* pl, double* pr, double* dif, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dtgsen_(ijob, wantq, wantz, select, n, a, lda, b, ldb, alphar, alphai, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info);
}
inline integer tgsen(integer* ijob, logical* wantq, logical* wantz,  logical* select, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* alpha, scomplex* beta, scomplex* q, integer* ldq, scomplex* z, integer* ldz, integer* m, float* pl, float* pr, float* dif, scomplex* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ctgsen_(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info);
}
inline integer tgsen(integer* ijob, logical* wantq, logical* wantz,  logical* select, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* alpha, dcomplex* beta, dcomplex* q, integer* ldq, dcomplex* z, integer* ldz, integer* m, double* pl, double* pr, double* dif, dcomplex* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return ztgsen_(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info);
}

// --- computes the generalized singular value decomposition (GSVD) of two real upper triangular ---
inline integer tgsja(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, integer* k, integer* l, float* a, integer* lda, float* b, integer* ldb, float* tola, float* tolb, float* alpha, float* beta, float* u, integer* ldu, float* v, integer* ldv, float* q, integer* ldq, float* work, integer* ncycle, integer* info)
{
  return stgsja_(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
}
inline integer tgsja(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, integer* k, integer* l, double* a, integer* lda, double* b, integer* ldb, double* tola, double* tolb, double* alpha, double* beta, double* u, integer* ldu, double* v, integer* ldv, double* q, integer* ldq, double* work, integer* ncycle, integer* info)
{
  return dtgsja_(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
}
inline integer tgsja(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, integer* k, integer* l, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* tola, float* tolb, float* alpha, float* beta, scomplex* u, integer* ldu, scomplex* v, integer* ldv, scomplex* q, integer* ldq, scomplex* work, integer* ncycle, integer* info)
{
  return ctgsja_(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
}
inline integer tgsja(char* jobu, char* jobv, char* jobq, integer* m, integer* p, integer* n, integer* k, integer* l, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* tola, double* tolb, double* alpha, double* beta, dcomplex* u, integer* ldu, dcomplex* v, integer* ldv, dcomplex* q, integer* ldq, dcomplex* work, integer* ncycle, integer* info)
{
  return ztgsja_(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
}

// --- estimates reciprocal condition numbers for specified eigenvalues and/or eigenvectors of a matrix pair ---
inline integer tgsna(char* job, char* howmny,  logical* select, integer* n,  float* a, integer* lda,  float* b, integer* ldb,  float* vl, integer* ldvl,  float* vr, integer* ldvr, float* s, float* dif, integer* mm, integer* m, float* work, integer* lwork, integer* iwork, integer* info)
{
  return stgsna_(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info);
}
inline integer tgsna(char* job, char* howmny,  logical* select, integer* n,  double* a, integer* lda,  double* b, integer* ldb,  double* vl, integer* ldvl,  double* vr, integer* ldvr, double* s, double* dif, integer* mm, integer* m, double* work, integer* lwork, integer* iwork, integer* info)
{
  return dtgsna_(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info);
}
inline integer tgsna(char* job, char* howmny,  logical* select, integer* n,  scomplex* a, integer* lda,  scomplex* b, integer* ldb,  scomplex* vl, integer* ldvl,  scomplex* vr, integer* ldvr, float* s, float* dif, integer* mm, integer* m, scomplex* work, integer* lwork, integer* iwork, integer* info)
{
  return ctgsna_(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info);
}
inline integer tgsna(char* job, char* howmny,  logical* select, integer* n,  dcomplex* a, integer* lda,  dcomplex* b, integer* ldb,  dcomplex* vl, integer* ldvl,  dcomplex* vr, integer* ldvr, double* s, double* dif, integer* mm, integer* m, dcomplex* work, integer* lwork, integer* iwork, integer* info)
{
  return ztgsna_(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info);
}

// ---  solves the generalized Sylvester equation ---
inline integer tgsyl(char* trans, integer* ijob, integer* m, integer* n,  float* a, integer* lda,  float* b, integer* ldb, float* c, integer* ldc,  float* d, integer* ldd,  float* e, integer* lde, float* f, integer* ldf, float* scale, float* dif, float* work, integer* lwork, integer* iwork, integer* info)
{
  return stgsyl_(trans, ijob, m, n,  a, lda,  b, ldb, c, ldc,  d, ldd,  e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
}
inline integer tgsyl(char* trans, integer* ijob, integer* m, integer* n,  double* a, integer* lda,  double* b, integer* ldb, double* c, integer* ldc,  double* d, integer* ldd,  double* e, integer* lde, double* f, integer* ldf, double* scale, double* dif, double* work, integer* lwork, integer* iwork, integer* info)
{
  return dtgsyl_(trans, ijob, m, n,  a, lda,  b, ldb, c, ldc,  d, ldd,  e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
}
inline integer tgsyl(char* trans, integer* ijob, integer* m, integer* n,  scomplex* a, integer* lda,  scomplex* b, integer* ldb, scomplex* c, integer* ldc,  scomplex* d, integer* ldd,  scomplex* e, integer* lde, scomplex* f, integer* ldf, float* scale, float* dif, scomplex* work, integer* lwork, integer* iwork, integer* info)
{
  return ctgsyl_(trans, ijob, m, n,  a, lda,  b, ldb, c, ldc,  d, ldd,  e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
}
inline integer tgsyl(char* trans, integer* ijob, integer* m, integer* n,  dcomplex* a, integer* lda,  dcomplex* b, integer* ldb, dcomplex* c, integer* ldc,  dcomplex* d, integer* ldd,  dcomplex* e, integer* lde, dcomplex* f, integer* ldf, double* scale, double* dif, dcomplex* work, integer* lwork, integer* iwork, integer* info)
{
  return ztgsyl_(trans, ijob, m, n,  a, lda,  b, ldb, c, ldc,  d, ldd,  e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
}

// --- estimates the reciprocal of the condition number of a packed triangular matrix A ---
inline integer tpcon(char* norm, char* uplo, char* diag, integer* n,  float* ap, float* rcond, float* work, integer* iwork, integer* info)
{
  return stpcon_(norm, uplo, diag, n,  ap, rcond, work, iwork, info);
}
inline integer tpcon(char* norm, char* uplo, char* diag, integer* n,  double* ap, double* rcond, double* work, integer* iwork, integer* info)
{
  return dtpcon_(norm, uplo, diag, n,  ap, rcond, work, iwork, info);
}
inline integer tpcon(char* norm, char* uplo, char* diag, integer* n,  scomplex* ap, float* rcond, scomplex* work, float* rwork, integer* info)
{
  return ctpcon_(norm, uplo, diag, n,  ap, rcond, work, rwork, info);
}
inline integer tpcon(char* norm, char* uplo, char* diag, integer* n,  dcomplex* ap, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  return ztpcon_(norm, uplo, diag, n,  ap, rcond, work, rwork, info);
}

// --- TPMQRT applies a real orthogonal matrix Q obtained from a triangular-pentagonal ---
inline integer tpmqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* l, integer* nb,  float* v, integer* ldv,  float* t, integer* ldt, float* a, integer* lda, float* b, integer* ldb, float* work, integer* info)
{
  return stpmqrt_(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline integer tpmqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* l, integer* nb,  double* v, integer* ldv,  double* t, integer* ldt, double* a, integer* lda, double* b, integer* ldb, double* work, integer* info)
{
  return dtpmqrt_(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline integer tpmqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* l, integer* nb,  scomplex* v, integer* ldv,  scomplex* t, integer* ldt, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* work, integer* info)
{
  return ctpmqrt_(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline integer tpmqrt(char* side, char* trans, integer* m, integer* n, integer* k, integer* l, integer* nb,  dcomplex* v, integer* ldv,  dcomplex* t, integer* ldt, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* work, integer* info)
{
  return ztpmqrt_(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}

// --- TPQRT computes a blocked QR factorization of a real triangular-pentagonal matrix ---
inline integer tpqrt(integer* m, integer* n, integer* l, integer* nb, float* a, integer* lda, float* b, integer* ldb, float* t, integer* ldt, float* work, integer* info)
{
  return stpqrt_(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info);
}
inline integer tpqrt(integer* m, integer* n, integer* l, integer* nb, double* a, integer* lda, double* b, integer* ldb, double* t, integer* ldt, double* work, integer* info)
{
  return dtpqrt_(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info);
}
inline integer tpqrt(integer* m, integer* n, integer* l, integer* nb, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* t, integer* ldt, scomplex* work, integer* info)
{
  return ctpqrt_(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info);
}
inline integer tpqrt(integer* m, integer* n, integer* l, integer* nb, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* t, integer* ldt, dcomplex* work, integer* info)
{
  return ztpqrt_(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info);
}

// --- computes a QR factorization of a real or complex "triangular-pentagonal" matrix ---
inline integer tpqrt2(integer* m, integer* n, integer* l, float* a, integer* lda, float* b, integer* ldb, float* t, integer* ldt, integer* info)
{
  return stpqrt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline integer tpqrt2(integer* m, integer* n, integer* l, double* a, integer* lda, double* b, integer* ldb, double* t, integer* ldt, integer* info)
{
  return dtpqrt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline integer tpqrt2(integer* m, integer* n, integer* l, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* t, integer* ldt, integer* info)
{
  return ctpqrt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline integer tpqrt2(integer* m, integer* n, integer* l, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* t, integer* ldt, integer* info)
{
  return ztpqrt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}

// --- applies a real or complex "triangular-pentagonal" blocked reflector to a real or complex matrix ---
inline integer tprfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k, integer* l,  float* v, integer* ldv,  float* t, integer* ldt, float* a, integer* lda, float* b, integer* ldb, float* work, integer* ldwork)
{
  return stprfb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork);
}
inline integer tprfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k, integer* l,  double* v, integer* ldv,  double* t, integer* ldt, double* a, integer* lda, double* b, integer* ldb, double* work, integer* ldwork)
{
  return dtprfb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork);
}
inline integer tprfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k, integer* l,  scomplex* v, integer* ldv,  scomplex* t, integer* ldt, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* work, integer* ldwork)
{
  return ctprfb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork);
}
inline integer tprfb(char* side, char* trans, char* direct, char* storev, integer* m, integer* n, integer* k, integer* l,  dcomplex* v, integer* ldv,  dcomplex* t, integer* ldt, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* work, integer* ldwork)
{
  return ztprfb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork);
}

// --- provides error bounds and backward error estimates for the solution to a system of linear equations ---
inline integer tprfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  float* ap,  float* b, integer* ldb,  float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return stprfs_(uplo, trans, diag, n, nrhs,  ap,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline integer tprfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  double* ap,  double* b, integer* ldb,  double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dtprfs_(uplo, trans, diag, n, nrhs,  ap,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline integer tprfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  scomplex* ap,  scomplex* b, integer* ldb,  scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return ctprfs_(uplo, trans, diag, n, nrhs,  ap,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}
inline integer tprfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  dcomplex* ap,  dcomplex* b, integer* ldb,  dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return ztprfs_(uplo, trans, diag, n, nrhs,  ap,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}

// --- computes the inverse of a real upper or lower triangular matrix A stored in packed format ---
inline integer tptri(char* uplo, char* diag, integer* n, float* ap, integer* info)
{
  return stptri_(uplo, diag, n, ap, info);
}
inline integer tptri(char* uplo, char* diag, integer* n, double* ap, integer* info)
{
  return dtptri_(uplo, diag, n, ap, info);
}
inline integer tptri(char* uplo, char* diag, integer* n, scomplex* ap, integer* info)
{
  return ctptri_(uplo, diag, n, ap, info);
}
inline integer tptri(char* uplo, char* diag, integer* n, dcomplex* ap, integer* info)
{
  return ztptri_(uplo, diag, n, ap, info);
}

// --- solves a triangular system of the form A * X = B  or  A**T * X = B ---
inline integer tptrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  float* ap, float* b, integer* ldb, integer* info)
{
  return stptrs_(uplo, trans, diag, n, nrhs,  ap, b, ldb, info);
}
inline integer tptrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  double* ap, double* b, integer* ldb, integer* info)
{
  return dtptrs_(uplo, trans, diag, n, nrhs,  ap, b, ldb, info);
}
inline integer tptrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  scomplex* ap, scomplex* b, integer* ldb, integer* info)
{
  return ctptrs_(uplo, trans, diag, n, nrhs,  ap, b, ldb, info);
}
inline integer tptrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  dcomplex* ap, dcomplex* b, integer* ldb, integer* info)
{
  return ztptrs_(uplo, trans, diag, n, nrhs,  ap, b, ldb, info);
}

// --- copies a triangular matrix from the standard packed format (TP) to the rectangular full packed format (TF) ---
inline integer tpttf(char* transr, char* uplo, integer* n,  float* ap, float* arf, integer* info)
{
  return stpttf_(transr, uplo, n,  ap, arf, info);
}
inline integer tpttf(char* transr, char* uplo, integer* n,  double* ap, double* arf, integer* info)
{
  return dtpttf_(transr, uplo, n,  ap, arf, info);
}
inline integer tpttf(char* transr, char* uplo, integer* n,  scomplex* ap, scomplex* arf, integer* info)
{
  return ctpttf_(transr, uplo, n,  ap, arf, info);
}
inline integer tpttf(char* transr, char* uplo, integer* n,  dcomplex* ap, dcomplex* arf, integer* info)
{
  return ztpttf_(transr, uplo, n,  ap, arf, info);
}

// --- estimates the reciprocal of the condition number of a triangular matrix A ---
inline integer trcon(char* norm, char* uplo, char* diag, integer* n,  float* a, integer* lda, float* rcond, float* work, integer* iwork, integer* info)
{
  return strcon_(norm, uplo, diag, n,  a, lda, rcond, work, iwork, info);
}
inline integer trcon(char* norm, char* uplo, char* diag, integer* n,  double* a, integer* lda, double* rcond, double* work, integer* iwork, integer* info)
{
  return dtrcon_(norm, uplo, diag, n,  a, lda, rcond, work, iwork, info);
}
inline integer trcon(char* norm, char* uplo, char* diag, integer* n,  scomplex* a, integer* lda, float* rcond, scomplex* work, float* rwork, integer* info)
{
  return ctrcon_(norm, uplo, diag, n,  a, lda, rcond, work, rwork, info);
}
inline integer trcon(char* norm, char* uplo, char* diag, integer* n,  dcomplex* a, integer* lda, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  return ztrcon_(norm, uplo, diag, n,  a, lda, rcond, work, rwork, info);
}

// --- computes some or all of the right and/or left eigenvectors of a real upper quasi-triangular matrix T ---
inline integer trevc(char* side, char* howmny, logical* select, integer* n,  float* t, integer* ldt, float* vl, integer* ldvl, float* vr, integer* ldvr, integer* mm, integer* m, float* work, integer* info)
{
  return strevc_(side, howmny, select, n,  t, ldt, vl, ldvl, vr, ldvr, mm, m, work, info);
}
inline integer trevc(char* side, char* howmny, logical* select, integer* n,  double* t, integer* ldt, double* vl, integer* ldvl, double* vr, integer* ldvr, integer* mm, integer* m, double* work, integer* info)
{
  return dtrevc_(side, howmny, select, n,  t, ldt, vl, ldvl, vr, ldvr, mm, m, work, info);
}
inline integer trevc(char* side, char* howmny, logical* select, integer* n,  scomplex* t, integer* ldt, scomplex* vl, integer* ldvl, scomplex* vr, integer* ldvr, integer* mm, integer* m, scomplex* work, float* rwork, integer* info)
{
  return ctrevc_(side, howmny, select, n,  t, ldt, vl, ldvl, vr, ldvr, mm, m, work, rwork, info);
}
inline integer trevc(char* side, char* howmny, logical* select, integer* n,  dcomplex* t, integer* ldt, dcomplex* vl, integer* ldvl, dcomplex* vr, integer* ldvr, integer* mm, integer* m, dcomplex* work, double* rwork, integer* info)
{
  return ztrevc_(side, howmny, select, n,  t, ldt, vl, ldvl, vr, ldvr, mm, m, work, rwork, info);
}

// --- reorders the real Schur factorization of a real matrix A = Q*T*Q**T ---
inline integer trexc(char* compq, integer* n, float* t, integer* ldt, float* q, integer* ldq, integer* ifst, integer* ilst, float* work, integer* info)
{
  return strexc_(compq, n, t, ldt, q, ldq, ifst, ilst, work, info);
}
inline integer trexc(char* compq, integer* n, double* t, integer* ldt, double* q, integer* ldq, integer* ifst, integer* ilst, double* work, integer* info)
{
  return dtrexc_(compq, n, t, ldt, q, ldq, ifst, ilst, work, info);
}
inline integer trexc(char* compq, integer* n, scomplex* t, integer* ldt, scomplex* q, integer* ldq, integer* ifst, integer* ilst, integer* info)
{
  return ctrexc_(compq, n, t, ldt, q, ldq, ifst, ilst, info);
}
inline integer trexc(char* compq, integer* n, dcomplex* t, integer* ldt, dcomplex* q, integer* ldq, integer* ifst, integer* ilst, integer* info)
{
  return ztrexc_(compq, n, t, ldt, q, ldq, ifst, ilst, info);
}

// --- provides error bounds and backward error estimates for the solution to a system of linear equations ---
inline integer trrfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  float* a, integer* lda,  float* b, integer* ldb,  float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return strrfs_(uplo, trans, diag, n, nrhs,  a, lda,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline integer trrfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  double* a, integer* lda,  double* b, integer* ldb,  double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dtrrfs_(uplo, trans, diag, n, nrhs,  a, lda,  b, ldb,  x, ldx, ferr, berr, work, iwork, info);
}
inline integer trrfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* b, integer* ldb,  scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return ctrrfs_(uplo, trans, diag, n, nrhs,  a, lda,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}
inline integer trrfs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* b, integer* ldb,  dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return ztrrfs_(uplo, trans, diag, n, nrhs,  a, lda,  b, ldb,  x, ldx, ferr, berr, work, rwork, info);
}

// --- reorders the real Schur factorization of a real matrix A = Q*T*Q**T ---
inline integer trsen(char* job, char* compq,  logical* select, integer* n, float* t, integer* ldt, float* q, integer* ldq, float* wr, float* wi, integer* m, float* s, float* sep, float* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return strsen_(job, compq,  select, n, t, ldt, q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info);
}
inline integer trsen(char* job, char* compq,  logical* select, integer* n, double* t, integer* ldt, double* q, integer* ldq, double* wr, double* wi, integer* m, double* s, double* sep, double* work, integer* lwork, integer* iwork, integer* liwork, integer* info)
{
  return dtrsen_(job, compq,  select, n, t, ldt, q, ldq, wr, wi, m, s, sep, work, lwork, iwork, liwork, info);
}
inline integer trsen(char* job, char* compq,  logical* select, integer* n, scomplex* t, integer* ldt, scomplex* q, integer* ldq, scomplex* w, integer* m, float* s, float* sep, scomplex* work, integer* lwork, integer* info)
{
  return ctrsen_(job, compq,  select, n, t, ldt, q, ldq, w, m, s, sep, work, lwork, info);
}
inline integer trsen(char* job, char* compq,  logical* select, integer* n, dcomplex* t, integer* ldt, dcomplex* q, integer* ldq, dcomplex* w, integer* m, double* s, double* sep, dcomplex* work, integer* lwork, integer* info)
{
  return ztrsen_(job, compq,  select, n, t, ldt, q, ldq, w, m, s, sep, work, lwork, info);
}

// --- estimates reciprocal condition numbers for specified eigenvalues ---
inline integer trsna(char* job, char* howmny,  logical* select, integer* n,  float* t, integer* ldt,  float* vl, integer* ldvl,  float* vr, integer* ldvr, float* s, float* sep, integer* mm, integer* m, float* work, integer* ldwork, integer* iwork, integer* info)
{
  return strsna_(job, howmny,  select, n,  t, ldt,  vl, ldvl,  vr, ldvr, s, sep, mm,  m, work, ldwork, iwork, info);
}
inline integer trsna(char* job, char* howmny,  logical* select, integer* n,  double* t, integer* ldt,  double* vl, integer* ldvl,  double* vr, integer* ldvr, double* s, double* sep, integer* mm, integer* m, double* work, integer* ldwork, integer* iwork, integer* info)
{
  return dtrsna_(job, howmny,  select, n,  t, ldt,  vl, ldvl,  vr, ldvr, s, sep, mm,  m, work, ldwork, iwork, info);
}
inline integer trsna(char* job, char* howmny,  logical* select, integer* n,  scomplex* t, integer* ldt,  scomplex* vl, integer* ldvl,  scomplex* vr, integer* ldvr, float* s, float* sep, integer* mm, integer* m, scomplex* work, integer* ldwork, float* rwork, integer* info)
{
  return ctrsna_(job, howmny,  select, n,  t, ldt,  vl, ldvl,  vr, ldvr, s, sep, mm,  m, work, ldwork, rwork, info);
}
inline integer trsna(char* job, char* howmny,  logical* select, integer* n,  dcomplex* t, integer* ldt,  dcomplex* vl, integer* ldvl,  dcomplex* vr, integer* ldvr, double* s, double* sep, integer* mm, integer* m, dcomplex* work, integer* ldwork, double* rwork, integer* info)
{
  return ztrsna_(job, howmny,  select, n,  t, ldt,  vl, ldvl,  vr, ldvr, s, sep, mm,  m, work, ldwork, rwork, info);
}

// --- solves a triangular system of the form A * X = B  or  A**T * X = B ---
inline integer trtrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  float* a, integer* lda, float* b, integer* ldb, integer* info)
{
  return strtrs_(uplo, trans, diag, n, nrhs,  a, lda, b, ldb, info);
}
inline integer trtrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  double* a, integer* lda, double* b, integer* ldb, integer* info)
{
  return dtrtrs_(uplo, trans, diag, n, nrhs,  a, lda, b, ldb, info);
}
inline integer trtrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  scomplex* a, integer* lda, scomplex* b, integer* ldb, integer* info)
{
  return ctrtrs_(uplo, trans, diag, n, nrhs,  a, lda, b, ldb, info);
}
inline integer trtrs(char* uplo, char* trans, char* diag, integer* n, integer* nrhs,  dcomplex* a, integer* lda, dcomplex* b, integer* ldb, integer* info)
{
  return ztrtrs_(uplo, trans, diag, n, nrhs,  a, lda, b, ldb, info);
}

// --- copies a triangular matrix from the standard full format (TR) to the rectangular full packed format (TF) ---
inline integer trttf(char* transr, char* uplo, integer* n,  float* a, integer* lda, float* arf, integer* info)
{
  return strttf_(transr, uplo, n,  a, lda, arf, info);
}
inline integer trttf(char* transr, char* uplo, integer* n,  double* a, integer* lda, double* arf, integer* info)
{
  return dtrttf_(transr, uplo, n,  a, lda, arf, info);
}
inline integer trttf(char* transr, char* uplo, integer* n,  scomplex* a, integer* lda, scomplex* arf, integer* info)
{
  return ctrttf_(transr, uplo, n,  a, lda, arf, info);
}
inline integer trttf(char* transr, char* uplo, integer* n,  dcomplex* a, integer* lda, dcomplex* arf, integer* info)
{
  return ztrttf_(transr, uplo, n,  a, lda, arf, info);
}

// --- copies a triangular matrix from the standard full format (TR) to the standard packed format (TP) ---
inline integer trttp(char* uplo, integer* n,  float* a, integer* lda, float* ap, integer* info)
{
  return strttp_(uplo, n,  a, lda, ap, info);
}
inline integer trttp(char* uplo, integer* n,  double* a, integer* lda, double* ap, integer* info)
{
  return dtrttp_(uplo, n,  a, lda, ap, info);
}
inline integer trttp(char* uplo, integer* n,  scomplex* a, integer* lda, scomplex* ap, integer* info)
{
  return ctrttp_(uplo, n,  a, lda, ap, info);
}
inline integer trttp(char* uplo, integer* n,  dcomplex* a, integer* lda, dcomplex* ap, integer* info)
{
  return ztrttp_(uplo, n,  a, lda, ap, info);
}

// --- reduces the M-by-N ( M<=N) real upper trapezoidal matrix A to upper triangular form---
inline integer tzrzf(integer* m, integer* n, float* a, integer* lda, float* tau, float* work, integer* lwork, integer* info)
{
  return stzrzf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer tzrzf(integer* m, integer* n, double* a, integer* lda, double* tau, double* work, integer* lwork, integer* info)
{
  return dtzrzf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer tzrzf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return ctzrzf_(m, n, a, lda, tau, work, lwork, info);
}
inline integer tzrzf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return ztzrzf_(m, n, a, lda, tau, work, lwork, info);
}

inline integer tzrqf(integer* m, integer* n, float* a, integer* lda, float* tau, integer* info)
{
  printf(" Function stzrqf() has been deprecated. Please use stzrzf() instead.\n"); 
  return stzrqf_(m, n, a, lda, tau, info);
}
inline integer tzrqf(integer* m, integer* n, double* a, integer* lda, double* tau, integer* info)
{
  printf(" Function dtzrqf() has been deprecated. Please use dtzrzf() instead.\n"); 
  return dtzrqf_(m, n, a, lda, tau, info);
}
inline integer tzrqf(integer* m, integer* n, scomplex* a, integer* lda, scomplex* tau, integer* info)
{
  printf(" Function ctzrqf() has been deprecated. Please use ctzrzf() instead.\n"); 
  return ctzrqf_(m, n, a, lda, tau, info);
}
inline integer tzrqf(integer* m, integer* n, dcomplex* a, integer* lda, dcomplex* tau, integer* info)
{
  printf(" Function ztzrqf() has been deprecated. Please use ztzrzf() instead.\n"); 
  return ztzrqf_(m, n, a, lda, tau, info);
}

// --- copies a triangular matrix from the standard packed format (TP) to the standard full format (TR) ---
inline integer tpttr(char* uplo, integer* n,  float* ap, float* a, integer* lda, integer* info)
{
  return stpttr_(uplo, n,  ap, a, lda, info);
}
inline integer tpttr(char* uplo, integer* n,  double* ap, double* a, integer* lda, integer* info)
{
  return dtpttr_(uplo, n,  ap, a, lda, info);
}
inline integer tpttr(char* uplo, integer* n,  scomplex* ap, scomplex* a, integer* lda, integer* info)
{
  return ctpttr_(uplo, n,  ap, a, lda, info);
}
inline integer tpttr(char* uplo, integer* n,  dcomplex* ap, dcomplex* a, integer* lda, integer* info)
{
  return ztpttr_(uplo, n,  ap, a, lda, info);
}

// --- improves the computed solution to a system of linear equations---
inline integer sprfs(char* uplo, integer* n, integer* nrhs,  float* ap,  float* afp,  integer* ipiv,  float* b, integer* ldb, float* x, integer* ldx, float* ferr, float* berr, float* work, integer* iwork, integer* info)
{
  return ssprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer sprfs(char* uplo, integer* n, integer* nrhs,  double* ap,  double* afp,  integer* ipiv,  double* b, integer* ldb, double* x, integer* ldx, double* ferr, double* berr, double* work, integer* iwork, integer* info)
{
  return dsprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, iwork, info);
}
inline integer sprfs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  scomplex* afp,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return csprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer sprfs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  dcomplex* afp,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zsprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer hprfs(char* uplo, integer* n, integer* nrhs,  scomplex* ap,  scomplex* afp,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return chprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer hprfs(char* uplo, integer* n, integer* nrhs,  dcomplex* ap,  dcomplex* afp,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zhprfs_(uplo, n, nrhs,  ap,  afp,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- generates random matrices with specified singular values for testing LAPACK programs ---
/* No F2C file for it - need to revisit.
inline integer latms(integer* m, integer* n, char* dist, integer* iseed, char* sym, float* d, integer* mode, float cond, float dmax, integer* kl, integer* ku, char* pack, float* a, integer* lda)
{
  return slatms_(m, n, dist, iseed, sym, d, mode, cond, dmax, kl, ku, pack, a, lda);
}
inline integer latms(integer* m, integer* n, char* dist, integer* iseed, char* sym, double* d, integer* mode, double cond, double dmax, integer* kl, integer* ku, char* pack, double* a, integer* lda)
{
  return dlatms_(m, n, dist, iseed, sym, d, mode, cond, dmax, kl, ku, pack, a, lda);
}
inline integer latms(integer* m, integer* n, char* dist, integer* iseed, char* sym, float* d, integer* mode, float cond, float dmax, integer* kl, integer* ku, char* pack, scomplex* a, integer* lda)
{
  return clatms_(m, n, dist, iseed, sym, d, mode, cond, dmax, kl, ku, pack, a, lda);
}
inline integer latms(integer* m, integer* n, char* dist, integer* iseed, char* sym, double* d, integer* mode, double cond, double dmax, integer* kl, integer* ku, char* pack, dcomplex* a, integer* lda)
{
  return zlatms_(m, n, dist, iseed, sym, d, mode, cond, dmax, kl, ku, pack, a, lda);
}*/

// --- generates a real general m by n matrix A ---
/* No F2C file for it - need to revisit.
inline integer lagge(integer* m, integer* n, integer* kl, integer* ku,  float* d, float* a, integer* lda, integer* iseed)
{
  return slagge_(m, n, kl, ku,  d, a, lda, iseed);
}
inline integer lagge(integer* m, integer* n, integer* kl, integer* ku,  double* d, double* a, integer* lda, integer* iseed)
{
  return dlagge_(m, n, kl, ku,  d, a, lda, iseed);
}
inline integer lagge(integer* m, integer* n, integer* kl, integer* ku,  float* d, scomplex* a, integer* lda, integer* iseed)
{
  return clagge_(m, n, kl, ku,  d, a, lda, iseed);
}
inline integer lagge(integer* m, integer* n, integer* kl, integer* ku,  double* d, dcomplex* a, integer* lda, integer* iseed)
{
  return zlagge_(m, n, kl, ku,  d, a, lda, iseed);
}*/

// ---  conjugates a complex vector ---
inline integer lacgv( integer* n, scomplex* x, integer*incx)
{
  return clacgv_(n, x, incx); 
}
inline integer lacgv( integer* n, dcomplex* x, integer*incx)
{
  return zlacgv_(n, x, incx); 
}

// ---  copies all or part of a real two-dimensional array to a complex array. ---
inline integer lacp2(char *uplo, integer*m, integer* n, float* a, integer*lda, scomplex* b, integer*ldb)
{
  return clacp2_(uplo, m, n, a, lda,  b, ldb); 
}
inline integer lacp2(char *uplo, integer*m, integer* n, double* a, integer*lda, dcomplex* b, integer*ldb)
{
  return zlacp2_(uplo, m, n, a, lda,  b, ldb);
}

// ---  multiplies a complex matrix by a square real matrix ---
inline integer lacrm(integer*m, integer* n, scomplex* a, integer*lda, float* b, integer*ldb, scomplex* c, integer*ldc, float* rwork)
{
  return clacrm_(m, n, a, lda,  b, ldb, c, ldc, rwork); 
}
inline integer lacrm(integer*m, integer* n, dcomplex* a, integer*lda, double* b, integer*ldb, dcomplex* c, integer*ldc, double* rwork)
{
  return zlacrm_(m, n, a, lda,  b, ldb, c, ldc, rwork); 
}

// --- generates a complex hermitian matrix A, by pre- and post- multiplying a real diagonal matrix D ---
/* No F2C file for it - need to revisit.
inline integer laghe(integer* n, integer*k, float* d, scomplex* a, integer*lda, integer* iseed)
{
  return claghe_(n, k, d, a, lda, iseed); 
}
inline integer laghe(integer* n, integer*k, double* d, dcomplex* a, integer*lda, integer* iseed)
{
  return zlaghe_(n, k, d, a, lda, iseed); 
}*/

// --- returns the value of the 1-norm, or the Frobenius norm, or the element of largest absolute value of a complex Hermitian matrix ---
inline float lanhe(char *norm, char *uplo, integer* n, scomplex* a, integer*lda, float* work)
{
  return clanhe_(norm, uplo, n, a, lda, work); 
}
inline double lanhe(char *norm, char *uplo, integer* n, dcomplex* a, integer*lda, double* work)
{
  return zlanhe_(norm, uplo, n, a, lda, work); 
}

// --- returns the value of the 1-norm, or the Frobenius norm, or the element of largest absolute value of a complex Hermitian matrix ---
inline integer larcm(integer*m, integer* n, float* a, integer*lda, scomplex* b, integer*ldb, scomplex* c, integer*ldc, float* rwork)
{
  return clarcm_(m, n, a, lda, b, ldb, c, ldc, rwork); 
}
inline integer larcm(integer*m, integer* n, double* a, integer*lda, dcomplex* b, integer*ldb, dcomplex* c, integer*ldc, double* rwork)
{
  return zlarcm_(m, n, a, lda, b, ldb, c, ldc, rwork); 
}

// --- estimates the reciprocal of the condition number of a general real matrix A ---
inline integer gecon(char* norm, integer* n,  float* a, integer* lda, float* anorm, float* rcond, float* work, integer* iwork, integer* info)
{
  return sgecon_(norm, n,  a, lda, anorm, rcond, work, iwork, info);
}
inline integer gecon(char* norm, integer* n,  double* a, integer* lda, double* anorm, double* rcond, double* work, integer* iwork, integer* info)
{
  return dgecon_(norm, n,  a, lda, anorm, rcond, work, iwork, info);
}
inline integer gecon(char* norm, integer* n,  scomplex* a, integer* lda, float* anorm, float* rcond, scomplex* work, float* rwork, integer* info)
{
  return cgecon_(norm, n,  a, lda, anorm, rcond, work, rwork, info);
}
inline integer gecon(char* norm, integer* n,  dcomplex* a, integer* lda, double* anorm, double* rcond, dcomplex* work, double* rwork, integer* info)
{
  return zgecon_(norm, n,  a, lda, anorm, rcond, work, rwork, info);
}

// --- computes row and column scalings intended to equilibrate an M-by-N matrix A and reduce its condition number ---
inline integer geequ(integer* m, integer* n,  float* a, integer* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  return sgeequ_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline integer geequ(integer* m, integer* n,  double* a, integer* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  return dgeequ_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline integer geequ(integer* m, integer* n,  scomplex* a, integer* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, integer* info)
{
  return cgeequ_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}
inline integer geequ(integer* m, integer* n,  dcomplex* a, integer* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, integer* info)
{
  return zgeequ_(m, n,  a, lda, r, c, rowcnd, colcnd, amax, info);
}

// --- estimates the reciprocal of the condition number of a complex Hermitian matrix A---
inline integer hecon(char* uplo, integer* n,  scomplex* a, integer* lda,  integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  return checon_(uplo, n, a, lda,  ipiv, anorm, rcond, work, info);
}
inline integer hecon(char* uplo, integer* n,  dcomplex* a, integer* lda,  integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  return zhecon_(uplo, n, a, lda,  ipiv, anorm, rcond, work, info);
}

// --- estimates the reciprocal of the condition number---
inline integer hecon_3(char* uplo, integer* n,  scomplex* a, integer* lda,  scomplex* e, integer* ipiv, float* anorm, float* rcond, scomplex* work, integer* info)
{
  return checon_3_(uplo, n, a, lda, e, ipiv, anorm, rcond, work, info);
}
inline integer hecon_3(char* uplo, integer* n,  dcomplex* a, integer* lda,  dcomplex* e, integer* ipiv, double* anorm, double* rcond, dcomplex* work, integer* info)
{
  return zhecon_3_(uplo, n, a, lda,  e, ipiv, anorm, rcond, work, info);
}

// --- computes row and column scalings intended to equilibrate a Hermitian matrix A ---
inline integer heequb(char* uplo, integer* n,  scomplex* a, integer* lda, float* s, float* scond, float* amax, scomplex* work, integer* info)
{
  return cheequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}
inline integer heequb(char* uplo, integer* n,  dcomplex* a, integer* lda, double* s, double* scond, double* amax, dcomplex* work, integer* info)
{
  return zheequb_(uplo, n,  a, lda, s, scond, amax, work, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline integer heev_2stage(char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, float* w, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return cheev_2stage_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}
inline integer heev_2stage(char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, double* w, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zheev_2stage_(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline integer heevd_2stage(char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, float* w, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return cheevd_2stage_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer heevd_2stage(char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, double* w, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zheevd_2stage_(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline integer heevr_2stage(char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, integer* isuppz, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return cheevr_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer heevr_2stage(char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, integer* isuppz, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zheevr_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline integer heevx(char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  return cheevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}
inline integer heevx(char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  return zheevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}

// --- computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices ---
inline integer heevx_2stage(char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  return cheevx_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}
inline integer heevx_2stage(char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  return zheevx_2stage_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors of a complex generalized Hermitian-definite eigenproblem ---
inline integer hegv(integer* itype, char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* w, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return chegv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
}
inline integer hegv(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* w, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zhegv_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors of a complex generalized Hermitian-definite eigenproblem ---
inline integer hegv_2stage(integer* itype, char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* w, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return chegv_2stage_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
}
inline integer hegv_2stage(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* w, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zhegv_2stage_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
}

// --- computes all the eigenvalues, and optionally, the eigenvectors of a complex generalized Hermitian-definite eigenproblem ---
inline integer hegvd(integer* itype, char* jobz, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* w, scomplex* work, integer* lwork, float* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return chegvd_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}
inline integer hegvd(integer* itype, char* jobz, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* w, dcomplex* work, integer* lwork, double* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info)
{
  return zhegvd_(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info);
}

// --- computes selected eigenvalues, and optionally, eigenvectors of a complex generalized Hermitian-definite eigenproblem ---
inline integer hegvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, scomplex* a, integer* lda, scomplex* b, integer* ldb, float* vl, float* vu, integer* il, integer* iu, float* abstol, integer* m, float* w, scomplex* z, integer* ldz, scomplex* work, integer* lwork, float* rwork, integer* iwork, integer* ifail, integer* info)
{
  return chegvx_(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}
inline integer hegvx(integer* itype, char* jobz, char* range, char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* vl, double* vu, integer* il, integer* iu, double* abstol, integer* m, double* w, dcomplex* z, integer* ldz, dcomplex* work, integer* lwork, double* rwork, integer* iwork, integer* ifail, integer* info)
{
  return zhegvx_(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
}

// --- improves the computed solution to a system of linear equations ---
inline integer herfs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* ferr, float* berr, scomplex* work, float* rwork, integer* info)
{
  return cherfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}
inline integer herfs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* ferr, double* berr, dcomplex* work, double* rwork, integer* info)
{
  return zherfs_(uplo, n, nrhs,  a, lda,  af, ldaf,  ipiv,  b, ldb, x, ldx, ferr, berr, work, rwork, info);
}

// --- improves the computed solution to a system of linear equations ---
inline integer herfsx(char* uplo, char* equed, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* af, integer* ldaf,  integer* ipiv,  float* s,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return cherfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer herfsx(char* uplo, char* equed, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* af, integer* ldaf,  integer* ipiv,  double* s,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zherfsx_(uplo, equed, n, nrhs,  a, lda,  af, ldaf,  ipiv,  s,  b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices---
inline integer hesv(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return chesv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer hesv(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zhesv_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices ---
inline integer hesv_aa(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return chesv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer hesv_aa(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zhesv_aa_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices ---
inline integer hesv_aa_2stage(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return chesv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}
inline integer hesv_aa_2stage(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zhesv_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for SY matrices ---
inline integer hesv_rk(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* e, integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return chesv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}
inline integer hesv_rk(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* e, integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zhesv_rk_(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices ---
inline integer hesvx(char* fact, char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv,  scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* ferr, float* berr, scomplex* work, integer* lwork, float* rwork, integer* info)
{
  return chesvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info);
}
inline integer hesvx(char* fact, char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv,  dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* ferr, double* berr, dcomplex* work, integer* lwork, double* rwork, integer* info)
{
  return zhesvx_(fact, uplo, n, nrhs,  a, lda, af, ldaf, ipiv,  b, ldb, x, ldx, rcond, ferr, berr, work, lwork, rwork, info);
}

// --- computes the solution to system of linear equations A * X = B for HE matrices ---
inline integer hesvxx(char* fact, char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* af, integer* ldaf, integer* ipiv, char* equed, float* s, scomplex* b, integer* ldb, scomplex* x, integer* ldx, float* rcond, float* rpvgrw, float* berr, integer* n_err_bnds, float* err_bnds_norm, float* err_bnds_comp, integer* nparams, float* params, scomplex* work, float* rwork, integer* info)
{
  return chesvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}
inline integer hesvxx(char* fact, char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* af, integer* ldaf, integer* ipiv, char* equed, double* s, dcomplex* b, integer* ldb, dcomplex* x, integer* ldx, double* rcond, double* rpvgrw, double* berr, integer* n_err_bnds, double* err_bnds_norm, double* err_bnds_comp, integer* nparams, double* params, dcomplex* work, double* rwork, integer* info)
{
  return zhesvxx_(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
}

// --- applies an elementary permutation on the rows and columns of a Hermitian matrix ---
inline integer heswapr(char* uplo, integer* n, scomplex* a, integer* lda, integer* i1, integer* i2)
{
  return cheswapr_(uplo, n, a, lda, i1, i2);
}
inline integer heswapr(char* uplo, integer* n, dcomplex* a, integer* lda, integer* i1, integer* i2)
{
  return zheswapr_(uplo, n, a, lda, i1, i2);
}

// --- HETRD reduces a complex Hermitian matrix A to real symmetric tridiagonal form T ---
inline integer hetrd(char* uplo, integer* n, scomplex* a, integer* lda, float* d, float* e, scomplex* tau, scomplex* work, integer* lwork, integer* info)
{
  return chetrd_(uplo, n, a, lda, d, e, tau, work, lwork, info);
}
inline integer hetrd(char* uplo, integer* n, dcomplex* a, integer* lda, double* d, double* e, dcomplex* tau, dcomplex* work, integer* lwork, integer* info)
{
  return zhetrd_(uplo, n, a, lda, d, e, tau, work, lwork, info);
}

// --- computes the factorization of a complex Hermitian matrix A ---
inline integer hetrf(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return chetrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer hetrf(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zhetrf_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the factorization of a complex hermitian matrix A ---
inline integer hetrf_aa(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return chetrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer hetrf_aa(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zhetrf_aa_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the factorization of a real hermitian matrix A ---
inline integer hetrf_aa_2stage(char* uplo, integer* n, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* work, integer* lwork, integer* info)
{
  return chetrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}
inline integer hetrf_aa_2stage(char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* work, integer* lwork, integer* info)
{
  return zhetrf_aa_2stage_(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
}

// --- computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman (rook) diagonal pivoting method ---
inline integer hetrf_rk(char* uplo, integer* n, scomplex* a, integer* lda, scomplex* e, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return chetrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}
inline integer hetrf_rk(char* uplo, integer* n, dcomplex* a, integer* lda, dcomplex* e, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zhetrf_rk_(uplo, n, a, lda, e, ipiv, work, lwork, info);
}

//--- computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method ---
inline integer hetrf_rook(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return chetrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer hetrf_rook(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zhetrf_rook_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the inverse of a complex Hermitian indefinite matrix A ---
inline integer hetri(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* info)
{
  return chetri_(uplo, n, a, lda, ipiv, work, info);
}
inline integer hetri(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* info)
{
  return zhetri_(uplo, n, a, lda, ipiv, work, info);
}

// --- computes the inverse of a complex Hermitian indefinite matrix A ---
inline integer hetri_3(char* uplo, integer* n, scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return chetri_3_(uplo, n, a, lda,  e,  ipiv, work, lwork, info);
}
inline integer hetri_3(char* uplo, integer* n, dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zhetri_3_(uplo, n, a, lda,  e,  ipiv, work, lwork, info);
}

// --- computes the inverse of a COMPLEX hermitian indefinite matrix A ---
inline integer hetri2(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* lwork, integer* info)
{
  return chetri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}
inline integer hetri2(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* lwork, integer* info)
{
  return zhetri2_(uplo, n, a, lda, ipiv, work, lwork, info);
}

// --- computes the inverse of a complex Hermitian indefinite matrix A ---
inline integer hetri2x(char* uplo, integer* n, scomplex* a, integer* lda, integer* ipiv, scomplex* work, integer* nb, integer* info)
{
  return chetri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}
inline integer hetri2x(char* uplo, integer* n, dcomplex* a, integer* lda, integer* ipiv, dcomplex* work, integer* nb, integer* info)
{
  return zhetri2x_(uplo, n, a, lda, ipiv, work, nb, info);
}

// --- solves a system of linear equations A*X = B with a complex Hermitian matrix A ---
inline integer hetrs(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return chetrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer hetrs(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zhetrs_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A * X = B with a complex Hermitian matrix A ---
inline integer hetrs_3(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  scomplex* e,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return chetrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv, b, ldb, info);
}
inline integer hetrs_3(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  dcomplex* e,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zhetrs_3_(uplo, n, nrhs,  a, lda,  e,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a complex hermitian matrix A ---
inline integer hetrs_aa(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* lwork, integer* info)
{
  return chetrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}
inline integer hetrs_aa(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* lwork, integer* info)
{
  return zhetrs_aa_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, lwork, info);
}

// --- solves a system of linear equations A*X = B with a real hermitian matrix A ---
inline integer hetrs_aa_2stage(char* uplo, integer* n, integer* nrhs, scomplex* a, integer* lda, scomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, scomplex* b, integer* ldb, integer* info)
{
  return chetrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}
inline integer hetrs_aa_2stage(char* uplo, integer* n, integer* nrhs, dcomplex* a, integer* lda, dcomplex* tb, integer* ltb, integer* ipiv, integer* ipiv2, dcomplex* b, integer* ldb, integer* info)
{
  return zhetrs_aa_2stage_(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
}

// --- computes the solution to a system of linear equations A * X = B for HE matrices ---
inline integer hetrs_rook(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, integer* info)
{
  return chetrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}
inline integer hetrs_rook(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, integer* info)
{
  return zhetrs_rook_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, info);
}

// --- solves a system of linear equations A*X = B with a complex Hermitian matrix A ---
inline integer hetrs2(char* uplo, integer* n, integer* nrhs,  scomplex* a, integer* lda,  integer* ipiv, scomplex* b, integer* ldb, scomplex* work, integer* info)
{
  return chetrs2_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, info);
}
inline integer hetrs2(char* uplo, integer* n, integer* nrhs,  dcomplex* a, integer* lda,  integer* ipiv, dcomplex* b, integer* ldb, dcomplex* work, integer* info)
{
  return zhetrs2_(uplo, n, nrhs,  a, lda,  ipiv, b, ldb, work, info);
}

// --- adds two scaled sum of squares quantities ---
inline integer combssq(float *v1, float *v2)
{
  return scombssq_(v1, v2);
}  
inline integer combssq(double *v1, double *v2)
{
  return dcombssq_(v1, v2);
}

// --- forms the 1-norm of the complex vector using the true absolute value ---
inline float sum1(integer *n, scomplex *cx, integer *incx)
{
  return scsum1_(n, cx, incx);
} 
inline double sum1(integer *n, dcomplex *cx, integer *incx)
{
  return dzsum1_(n, cx, incx);
} 

// --- computes the LU factorization of a general band matrix using the unblocked version of the algorithm ---
inline integer gbtf2(integer *m, integer *n, integer *kl, integer *ku, float *ab, integer *ldab, integer *ipiv, integer *info)
{
  return sgbtf2_(m, n, kl, ku, ab, ldab, ipiv, info); 
} 
inline integer gbtf2(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, integer *ipiv, integer *info)
{
  return dgbtf2_(m, n, kl, ku, ab, ldab, ipiv, info);  
}
inline integer gbtf2(integer *m, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, integer *ipiv, integer *info)
{
  return cgbtf2_(m, n, kl, ku, ab, ldab, ipiv, info);  
}
inline integer gbtf2(integer *m, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, integer *ipiv, integer *info)
{
  return zgbtf2_(m, n, kl, ku, ab, ldab, ipiv, info);  
}

// --- recursively computes a LQ factorization of a float M-by-N matrix A ---
inline integer gelqt3(integer *m, integer *n, float *a, integer *lda, float *t, integer *ldt, integer *info)
{
  return sgelqt3_(m, n, a, lda, t, ldt, info);
}
inline integer gelqt3(integer *m, integer *n, double *a, integer *lda, double *t, integer *ldt, integer *info)
{
  return dgelqt3_(m, n, a, lda, t, ldt, info);
}
inline integer gelqt3(integer *m, integer *n, scomplex *a, integer *lda, scomplex *t, integer *ldt, integer *info)
{
  return cgelqt3_(m, n, a, lda, t, ldt, info);
}
inline integer gelqt3(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, integer *info)
{
  return zgelqt3_(m, n, a, lda, t, ldt, info);
}

// --- computes the QL factorization of a general rectangular matrix using an unblocked algorithm ---
inline integer geql2(integer *m, integer *n, float *a, integer *lda, float *tau, float *work, integer *info)
{
  return sgeql2_(m, n, a, lda, tau, work, info);
}
inline integer geql2(integer *m, integer *n, double *a, integer *lda, double *tau, double *work, integer *info)
{
  return dgeql2_(m, n, a, lda, tau, work, info);
}
inline integer geql2(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  return cgeql2_(m, n, a, lda, tau, work, info);
}
inline integer geql2(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  return zgeql2_(m, n, a, lda, tau, work, info);
}

// --- computes the QR factorization of a general rectangular matrix with non-negative diagonal elements using an unblocked algorithm ---
inline integer geqr2p(integer *m, integer *n, float *a, integer *lda, float *tau, float *work, integer *info)
{
  return sgeqr2p_(m, n, a, lda, tau, work, info);
}
inline integer geqr2p(integer *m, integer *n, double *a, integer *lda, double *tau, double *work, integer *info)
{
  return dgeqr2p_(m, n, a, lda, tau, work, info);
}
inline integer geqr2p(integer *m, integer *n, scomplex *a, integer * lda, scomplex *tau, scomplex *work, integer *info)
{
  return cgeqr2p_(m, n, a, lda, tau, work, info);
}
inline integer geqr2p(integer *m, integer *n, dcomplex *a, integer * lda, dcomplex *tau, dcomplex *work, integer *info)
{
  return zgeqr2p_(m, n, a, lda, tau, work, info);
}

// --- computes the RQ factorization of a general rectangular matrix using an unblocked algorithm ---
inline integer gerq2(integer *m, integer *n, float *a, integer *lda, float *tau, float *work, integer *info)
{
  return sgerq2_(m, n, a, lda, tau, work, info);
}
inline integer gerq2(integer *m, integer *n, double *a, integer *lda, double *tau, double *work, integer *info)
{
  return dgerq2_(m, n, a, lda, tau, work, info);
}
inline integer gerq2(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  return cgerq2_(m, n, a, lda, tau, work, info);
}
inline integer gerq2(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  return zgerq2_(m, n, a, lda, tau, work, info);
}

// --- solves a system of linear equations using the LU factorization with complete pivoting computed by getc2 ---
inline integer gesc2(integer *n, float *a, integer *lda, float *rhs, integer *ipiv, integer *jpiv, float *scale)
{
  return sgesc2_(n, a, lda, rhs, ipiv, jpiv, scale);
}
inline integer gesc2(integer *n, double *a, integer *lda, double *rhs, integer *ipiv, integer *jpiv, double *scale)
{
  return dgesc2_(n, a, lda, rhs, ipiv, jpiv, scale);
}
inline integer gesc2(integer *n, scomplex *a, integer *lda, scomplex *rhs, integer *ipiv, integer *jpiv, float *scale)
{
  return cgesc2_(n, a, lda, rhs, ipiv, jpiv, scale);
}
inline integer gesc2(integer *n, dcomplex *a, integer *lda, dcomplex *rhs, integer *ipiv, integer *jpiv, double *scale)
{
  return zgesc2_(n, a, lda, rhs, ipiv, jpiv, scale);
}

// --- pre-processor for the routine gesvj ---
inline integer gsvj0(char *jobv, integer *m, integer *n, float *a, integer *lda, float *d, float *sva, integer *mv, float *v, integer *ldv, float *eps, float *sfmin, float *tol, integer *nsweep, float *work, integer *lwork, integer *info)
{
  return sgsvj0_(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline integer gsvj0(char *jobv, integer *m, integer *n, double *a, integer *lda, double *d, double *sva, integer *mv, double *v, integer *ldv, double *eps, double *sfmin, double *tol, integer *nsweep, double *work, integer *lwork, integer *info)
{
  return dgsvj0_(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline integer gsvj0(char *jobv, integer *m, integer *n, scomplex *a, integer *lda, scomplex *d, float *sva, integer *mv, scomplex *v, integer *ldv, float *eps, float *sfmin, float *tol, integer *nsweep, scomplex *work, integer *lwork, integer *info)
{
  return cgsvj0_(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline integer gsvj0(char *jobv, integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *d, double *sva, integer *mv, dcomplex *v, integer *ldv, double *eps, double *sfmin, double *tol, integer *nsweep, dcomplex *work, integer *lwork, integer *info)
{
  return zgsvj0_(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}

// --- pre-processor for the routine sgesvj, applies Jacobi rotations targeting only particular pivots ---
inline integer gsvj1(char *jobv, integer *m, integer *n, integer *n1, float *a, integer *lda, float *d, float *sva, integer *mv, float *v, integer *ldv, float *eps, float *sfmin, float *tol, integer *nsweep, float *work, integer *lwork, integer *info)
{
  return sgsvj1_(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline integer gsvj1(char *jobv, integer *m, integer *n, integer *n1, double *a, integer *lda, double *d, double *sva, integer *mv, double *v, integer *ldv, double *eps, double *sfmin, double *tol, integer *nsweep, double *work, integer *lwork, integer *info)
{
  return dgsvj1_(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline integer gsvj1(char *jobv, integer *m, integer *n, integer *n1, scomplex *a, integer *lda, scomplex *d, float *sva, integer *mv, scomplex *v, integer *ldv, float *eps, float *sfmin, float *tol, integer *nsweep, scomplex *work, integer *lwork, integer *info)
{
  return cgsvj1_(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}
inline integer gsvj1(char *jobv, integer *m, integer *n, integer *n1, dcomplex *a, integer *lda, dcomplex *d, double *sva, integer *mv, dcomplex *v, integer *ldv, double *eps, double *sfmin, double *tol, integer *nsweep, dcomplex *work, integer *lwork, integer *info)
{
  return zgsvj1_(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
}

// --- solves a system of linear equations with a tridiagonal matrix using the LU factorization computed by sgttrf ---
inline integer gtts2(integer *itrans, integer *n, integer *nrhs, float *dl, float *d, float *du, float *du2, integer *ipiv, float *b, integer * ldb)
{
  return sgtts2_(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}
inline integer gtts2(integer *itrans, integer *n, integer *nrhs, double *dl, double *d, double *du, double *du2, integer *ipiv, double *b, integer * ldb)
{
  return dgtts2_(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}
inline integer gtts2(integer *itrans, integer *n, integer *nrhs, scomplex *dl, scomplex *d, scomplex *du, scomplex *du2, integer *ipiv, scomplex *b, integer * ldb)
{
  return cgtts2_(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}
inline integer gtts2(integer *itrans, integer *n, integer *nrhs, dcomplex *dl, dcomplex *d, dcomplex *du, dcomplex *du2, integer *ipiv, dcomplex *b, integer * ldb)
{
  return zgtts2_(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}

// --- internal routine used by the HETRD_HB2ST subroutine ---
inline integer hb2st_kernels(char *uplo, logical *wantz, integer * ttype, integer *st, integer *ed, integer *sweep, integer *n, integer * nb, integer *ib, scomplex *a, integer *lda, scomplex *v, scomplex *tau, integer *ldvt, scomplex *work)
{
  return chb2st_kernels_(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
} 
inline integer hb2st_kernels(char *uplo, logical *wantz, integer * ttype, integer *st, integer *ed, integer *sweep, integer *n, integer * nb, integer *ib, dcomplex *a, integer *lda, dcomplex *v, dcomplex *tau, integer *ldvt, dcomplex *work)
{
  return zhb2st_kernels_(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
}

// --- estimates the reciprocal of the condition number fort HE matrices ---
inline integer hecon_rook(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, float *anorm, float *rcond, scomplex *work, integer *info)
{
  return checon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, info); 
} 
inline integer hecon_rook(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, double *anorm, double *rcond, dcomplex *work, integer *info)
{
  return zhecon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, info);  
}
 
// --- computes the solution to a system of linear equations A * X = B for HE matrices ---
inline integer hesv_rook(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv, scomplex *b, integer *ldb, scomplex *work, integer *lwork, integer *info)
{
  return chesv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}
inline integer hesv_rook(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv, dcomplex *b, integer *ldb, dcomplex *work, integer *lwork, integer *info)
{
  return zhesv_rook_(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
}

// --- computes the factorization of a scomplex Hermitian matrix ---
inline integer hetf2(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info)
{
  return chetf2_(uplo, n, a, lda, ipiv, info); 
}
inline integer hetf2(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info)
{
  return zhetf2_(uplo, n, a, lda, ipiv, info); 
}

// --- computes the factorization of a scomplex Hermitian indefinite matrix ---
inline integer hetf2_rk(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *e, integer *ipiv, integer *info)
{
  return chetf2_rk_(uplo, n, a, lda, e, ipiv, info); 
}
inline integer hetf2_rk(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, integer *info)
{
  return zhetf2_rk_(uplo, n, a, lda, e, ipiv, info);  
}

// --- computes the factorization of a scomplex Hermitian indefinite matrix ---
inline integer hetf2_rook(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info)
{
  return chetf2_rook_(uplo, n, a, lda, ipiv, info); 
}
inline integer hetf2_rook(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info)
{
  return zhetf2_rook_(uplo, n, a, lda, ipiv, info); 
}

// --- reduces a scomplex Hermitian matrix A to float symmetric tridiagonal form T  ---
inline integer hetrd_2stage(char *vect, char *uplo, integer *n, scomplex *a, integer *lda, float *d, float *e, scomplex *tau, scomplex *hous2, integer *lhous2, scomplex *work, integer *lwork, integer *info)
{
  return chetrd_2stage_(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);   
}
inline integer hetrd_2stage(char *vect, char *uplo, integer *n, dcomplex *a, integer *lda, double *d, double *e, dcomplex *tau, dcomplex *hous2, integer *lhous2, dcomplex *work, integer *lwork, integer *info)
{
  return zhetrd_2stage_(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);   
} 

// --- reduces a scomplex Hermitian band matrix A to float symmetric tridiagonal form T ---
inline integer hetrd_hb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *d, float * e, scomplex *hous, integer *lhous, scomplex *work, integer *lwork, integer *info)
{
  return chetrd_hb2st_(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
} 
inline integer hetrd_hb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *d, double * e, dcomplex *hous, integer *lhous, dcomplex *work, integer *lwork, integer *info)
{
  return zhetrd_hb2st_(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
}

// --- reduces a scomplex Hermitian matrix A to scomplex Hermitian band-diagonal form AB ---
inline integer hetrd_he2hb(char *uplo, integer *n, integer *kd, scomplex *a, integer *lda, scomplex *ab, integer *ldab, scomplex *tau, scomplex *work, integer *lwork, integer *info)
{
  return chetrd_he2hb_(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
}
inline integer hetrd_he2hb(char *uplo, integer *n, integer *kd, dcomplex *a, integer *lda, dcomplex *ab, integer *ldab, dcomplex *tau, dcomplex *work, integer *lwork, integer *info)
{
  return zhetrd_he2hb_(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
}

// --- computes the inverse of a scomplex Hermitian indefinite matrix A ---
inline integer hetri_3x(char *uplo, integer *n, scomplex *a, integer * lda, scomplex *e, integer *ipiv, scomplex *work, integer *nb, integer * info)
{
  return chetri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}
inline integer hetri_3x(char *uplo, integer *n, dcomplex *a, integer * lda, dcomplex *e, integer *ipiv, dcomplex *work, integer *nb, integer * info)
{
  return zhetri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}

// --- computes the inverse of HE matrix using the factorization obtained with the bounded Bunch-Kaufman ---
inline integer hetri_rook(char *uplo, integer *n, scomplex *a, integer * lda, integer *ipiv, scomplex *work, integer * info)
{
  return chetri_rook_(uplo, n, a, lda, ipiv, work, info);
}
inline integer hetri_rook(char *uplo, integer *n, dcomplex *a, integer * lda, integer *ipiv, dcomplex *work, integer * info)
{
  return zhetri_rook_(uplo, n, a, lda, ipiv, work, info);
}

// ---  tests input for NaN ---
inline integer isnan(float *sin)
{
  return sisnan_(sin);
}
inline integer isnan(double *sin)
{
  return disnan_(sin);
}

// --- computes the Cholesky factorization of a symmetric/Hermitian positive definite band matrix (unblocked algorithm) ---
inline integer pbtf2(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, integer *info)
{
  return spbtf2_(uplo, n, kd, ab, ldab, info);
}
inline integer pbtf2(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, integer *info)
{
  return dpbtf2_(uplo, n, kd, ab, ldab, info);
}
inline integer pbtf2(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, integer *info)
{
  return cpbtf2_(uplo, n, kd, ab, ldab, info);
}
inline integer pbtf2(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, integer *info)
{
  return zpbtf2_(uplo, n, kd, ab, ldab, info);
}

// --- PSTF2 computes the Cholesky factorization with complete pivoting of a float symmetric positive semidefinite matrix ---
inline integer pstf2(char *uplo, integer *n, float *a, integer *lda, integer *piv, integer *rank, float *tol, float *work, integer *info)
{
  return spstf2_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline integer pstf2(char *uplo, integer *n, double *a, integer *lda, integer *piv, integer *rank, double *tol, double *work, integer *info)
{
  return dpstf2_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline integer pstf2(char *uplo, integer *n, scomplex *a, integer *lda, integer *piv, integer *rank, float *tol, float *work, integer *info)
{
  return cpstf2_(uplo, n, a, lda, piv, rank, tol, work, info);
}
inline integer pstf2(char *uplo, integer *n, dcomplex *a, integer *lda, integer *piv, integer *rank, double *tol, double *work, integer *info)
{
  return zpstf2_(uplo, n, a, lda, piv, rank, tol, work, info);
}

// --- solves a tridiagonal system of the form AX=B using the L D LH factorization computed by pttrf ---
inline integer ptts2(integer *n, integer *nrhs, float *d, float *e, float *b, integer *ldb)
{
  return sptts2_(n, nrhs, d, e, b, ldb);
}
inline integer ptts2(integer *n, integer *nrhs, double *d, double *e, double *b, integer *ldb)
{
  return dptts2_(n, nrhs, d, e, b, ldb);
}
inline integer ptts2(integer* iuplo, integer *n, integer *nrhs, float *d, scomplex *e, scomplex *b, integer *ldb)
{
  return cptts2_(iuplo, n, nrhs, d, e, b, ldb);
}
inline integer ptts2(integer* iuplo, integer *n, integer *nrhs, double *d, dcomplex *e, dcomplex *b, integer *ldb)
{
  return zptts2_(iuplo, n, nrhs, d, e, b, ldb);
}

// --- applies a plane rotation with float cosine and scomplex sine to a pair of scomplex vectors ---
inline integer rot(integer *n, scomplex *cx, integer *incx, scomplex * cy, integer *incy, float *c, scomplex *s)
{
  return crot_(n, cx, incx, cy, incy, c, s);
}
inline integer rot(integer *n, dcomplex *cx, integer *incx, dcomplex * cy, integer *incy, double *c, dcomplex *s)
{
  return zrot_(n, cx, incx, cy, incy, c, s);
}

// --- multiplies a vector by the reciprocal of a float scalar ---
inline integer rscl(integer *n, float *sa, float *sx, integer *incx)
{
  return srscl_(n, sa, sx, incx);
}
inline integer rscl(integer *n, double *sa, double *sx, integer *incx)
{
  return drscl_(n, sa, sx, incx);
}
inline integer rscl(integer *n, float *sa, scomplex *sx, integer *incx)
{
  return csrscl_(n, sa, sx, incx);
}
inline integer rscl(integer *n, double *sa, dcomplex *sx, integer *incx)
{
  return zdrscl_(n, sa, sx, incx);
}

// --- internal routine used by the SYTRD_SB2ST subroutine ---
inline integer sb2st_kernels(char *uplo, logical *wantz, integer * ttype, integer *st, integer *ed, integer *sweep, integer *n, integer * nb, integer *ib, float *a, integer *lda, float *v, float *tau, integer *ldvt, float *work)
{
  return ssb2st_kernels_(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
} 
inline integer sb2st_kernels(char *uplo, logical *wantz, integer * ttype, integer *st, integer *ed, integer *sweep, integer *n, integer * nb, integer *ib, double *a, integer *lda, double *v, double *tau, integer *ldvt, double *work)
{
  return dsb2st_kernels_(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work);
}

// --- computes a matrix-vector product for scomplex vectors using a scomplex symmetric packed matrix ---
inline integer spmv(char *uplo, integer *n, scomplex *alpha, scomplex *ap, scomplex *x, integer *incx, scomplex *beta, scomplex *y, integer * incy)
{
  return cspmv_(uplo, n, alpha, ap, x, incx, beta, y, incy);
}
inline integer spmv(char *uplo, integer *n, dcomplex *alpha, dcomplex *ap, dcomplex *x, integer *incx, dcomplex *beta, dcomplex *y, integer * incy)
{
  return zspmv_(uplo, n, alpha, ap, x, incx, beta, y, incy);
}

// --- performs the symmetrical rank-1 update of a scomplex symmetric packed matrix ---
inline integer spr(char *uplo, integer *n, scomplex *alpha, scomplex *x, integer *incx, scomplex *ap)
{
  return cspr_(uplo, n, alpha, x, incx, ap);
}
inline integer spr(char *uplo, integer *n, dcomplex *alpha, dcomplex *x, integer *incx, dcomplex *ap)
{
  return zspr_(uplo, n, alpha, x, incx, ap);
}

// --- computes a matrix-vector product for a scomplex symmetric matrix ---
inline integer symv(char *uplo, integer *n, scomplex *alpha, scomplex * a, integer *lda, scomplex *x, integer *incx, scomplex *beta, scomplex *y, integer *incy)
{
  return csymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy); 
}
inline integer symv(char *uplo, integer *n, dcomplex *alpha, dcomplex * a, integer *lda, dcomplex *x, integer *incx, dcomplex *beta, dcomplex *y, integer *incy)
{
  return zsymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy); 
}

// --- reduces a float symmetric matrix A to float symmetric tridiagonal form T  ---
inline integer sytrd_2stage(char *vect, char *uplo, integer *n, float *a, integer *lda, float *d, float *e, float *tau, float *hous2, integer *lhous2, float *work, integer *lwork, integer *info)
{
  return ssytrd_2stage_(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);   
}
inline integer sytrd_2stage(char *vect, char *uplo, integer *n, double *a, integer *lda, double *d, double *e, double *tau, double *hous2, integer *lhous2, double *work, integer *lwork, integer *info)
{
  return dsytrd_2stage_(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info);   
}

// --- reduces a float symmetric band matrix A to float symmetric tridiagonal form T ---
inline integer sytrd_sb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, float *ab, integer *ldab, float *d, float * e, float *hous, integer *lhous, float *work, integer *lwork, integer *info)
{
  return ssytrd_sb2st_(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
} 
inline integer sytrd_sb2st(char *stage1, char *vect, char *uplo, integer *n, integer *kd, double *ab, integer *ldab, double *d, double * e, double *hous, integer *lhous, double *work, integer *lwork, integer *info)
{
  return dsytrd_sb2st_(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info);
}

// --- reduces a float symmetric matrix A to float symmetric band-diagonal form AB ---
inline integer sytrd_sy2sb(char *uplo, integer *n, integer *kd, float *a, integer *lda, float *ab, integer *ldab, float *tau, float *work, integer *lwork, integer *info)
{
  return ssytrd_sy2sb_(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
}
inline integer sytrd_sy2sb(char *uplo, integer *n, integer *kd, double *a, integer *lda, double *ab, integer *ldab, double *tau, double *work, integer *lwork, integer *info)
{
  return dsytrd_sy2sb_(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info);
}

// --- estimates the reciprocal of the condition number of a float symmetric matrix A ---
inline integer sycon_rook(char *uplo, integer *n, float *a, integer * lda, integer *ipiv, float *anorm, float *rcond, float *work, integer * iwork, integer *info)
{
  return ssycon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, iwork, info);
}
inline integer sycon_rook(char *uplo, integer *n, double *a, integer * lda, integer *ipiv, double *anorm, double *rcond, double *work, integer * iwork, integer *info)
{
  return dsycon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, iwork, info);
}
inline integer sycon_rook(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, float *anorm, float *rcond, scomplex *work, integer *info)
{
  return csycon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, info);
}
inline integer sycon_rook(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, double *anorm, double *rcond, dcomplex *work, integer *info)
{
  return zsycon_rook_(uplo, n, a, lda, ipiv, anorm, rcond, work, info);
}

// --- converts the factorization output format ---
inline integer syconvf(char *uplo, char *way, integer *n, float *a, integer *lda, float *e, integer *ipiv, integer *info)
{
  return ssyconvf_(uplo, way, n, a, lda, e, ipiv, info);
}
inline integer syconvf(char *uplo, char *way, integer *n, double *a, integer *lda, double *e, integer *ipiv, integer *info)
{
  return dsyconvf_(uplo, way, n, a, lda, e, ipiv, info);
}
inline integer syconvf(char *uplo, char *way, integer *n, scomplex *a, integer *lda, scomplex *e, integer *ipiv, integer *info)
{
  return csyconvf_(uplo, way, n, a, lda, e, ipiv, info);
}
inline integer syconvf(char *uplo, char *way, integer *n, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, integer *info)
{
  return zsyconvf_(uplo, way, n, a, lda, e, ipiv, info);
}

// --- converts the factorization output format ---
inline integer syconvf_rook(char *uplo, char *way, integer *n, float *a, integer *lda, float *e, integer *ipiv, integer *info)
{
  return ssyconvf_rook_(uplo, way, n, a, lda, e, ipiv, info);
}
inline integer syconvf_rook(char *uplo, char *way, integer *n, double *a, integer *lda, double *e, integer *ipiv, integer *info)
{
  return dsyconvf_rook_(uplo, way, n, a, lda, e, ipiv, info);
}
inline integer syconvf_rook(char *uplo, char *way, integer *n, scomplex *a, integer *lda, scomplex *e, integer *ipiv, integer *info)
{
  return csyconvf_rook_(uplo, way, n, a, lda, e, ipiv, info);
}
inline integer syconvf_rook(char *uplo, char *way, integer *n, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, integer *info)
{
  return zsyconvf_rook_(uplo, way, n, a, lda, e, ipiv, info);
}

// --- computes the factorization of a float symmetric indefinite matrix, using the diagonal pivoting method ---
inline integer sytf2(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, integer *info)
{
  return ssytf2_(uplo, n, a, lda, ipiv, info);
}
inline integer sytf2(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, integer *info)
{
  return dsytf2_(uplo, n, a, lda, ipiv, info);
}
inline integer sytf2(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info)
{
  return csytf2_(uplo, n, a, lda, ipiv, info);
}
inline integer sytf2(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info)
{
  return zsytf2_(uplo, n, a, lda, ipiv, info);
}

// --- computes the factorization of a float symmetric indefinite matrix ---
inline integer sytf2_rk(char *uplo, integer *n, float *a, integer * lda, float *e, integer *ipiv, integer *info)
{
  return ssytf2_rk_(uplo, n, a, lda, e, ipiv, info);
}
inline integer sytf2_rk(char *uplo, integer *n, double *a, integer * lda, double *e, integer *ipiv, integer *info)
{
  return dsytf2_rk_(uplo, n, a, lda, e, ipiv, info);
}
inline integer sytf2_rk(char *uplo, integer *n, scomplex *a, integer * lda, scomplex *e, integer *ipiv, integer *info)
{
  return csytf2_rk_(uplo, n, a, lda, e, ipiv, info);
}
inline integer sytf2_rk(char *uplo, integer *n, dcomplex *a, integer * lda, dcomplex *e, integer *ipiv, integer *info)
{
  return zsytf2_rk_(uplo, n, a, lda, e, ipiv, info);
}

// --- computes the factorization of a float symmetric indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method ---
inline integer sytf2_rook(char *uplo, integer *n, float *a, integer * lda, integer *ipiv, integer *info)
{
  return ssytf2_rook_(uplo, n, a, lda, ipiv, info);
}
inline integer sytf2_rook(char *uplo, integer *n, double *a, integer * lda, integer *ipiv, integer *info)
{
  return dsytf2_rook_(uplo, n, a, lda, ipiv, info);
}
inline integer sytf2_rook(char *uplo, integer *n, scomplex *a, integer * lda, integer *ipiv, integer *info)
{
  return csytf2_rook_(uplo, n, a, lda, ipiv, info);
}
inline integer sytf2_rook(char *uplo, integer *n, dcomplex *a, integer * lda, integer *ipiv, integer *info)
{
  return zsytf2_rook_(uplo, n, a, lda, ipiv, info);
}

// --- computes the inverse of a float symmetric indefinite matrix A ---
inline integer sytri_3x(char *uplo, integer *n, float *a, integer * lda, float *e, integer *ipiv, float *work, integer *nb, integer *info)
{
  return ssytri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}
inline integer sytri_3x(char *uplo, integer *n, double *a, integer * lda, double *e, integer *ipiv, double *work, integer *nb, integer *info)
{
  return dsytri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}
inline integer sytri_3x(char *uplo, integer *n, scomplex *a, integer * lda, scomplex *e, integer *ipiv, scomplex *work, integer *nb, integer *info)
{
  return csytri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}
inline integer sytri_3x(char *uplo, integer *n, dcomplex *a, integer * lda, dcomplex *e, integer *ipiv, dcomplex *work, integer *nb, integer *info)
{
  return zsytri_3x_(uplo, n, a, lda, e, ipiv, work, nb, info);
}

// --- computes the inverse of a float symmetric matrix A ---
inline integer sytri_rook(char *uplo, integer *n, float *a, integer * lda, integer *ipiv, float *work, integer *info)
{
  return ssytri_rook_(uplo, n, a, lda, ipiv, work, info);
}
inline integer sytri_rook(char *uplo, integer *n, double *a, integer * lda, integer *ipiv, double *work, integer *info)
{
  return dsytri_rook_(uplo, n, a, lda, ipiv, work, info);
}
inline integer sytri_rook(char *uplo, integer *n, scomplex *a, integer * lda, integer *ipiv, scomplex *work, integer *info)
{
  return csytri_rook_(uplo, n, a, lda, ipiv, work, info);
}
inline integer sytri_rook(char *uplo, integer *n, dcomplex *a, integer * lda, integer *ipiv, dcomplex *work, integer *info)
{
  return zsytri_rook_(uplo, n, a, lda, ipiv, work, info);
}

// --- swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an orthogonal equivalence transformation ---
inline integer tgex2(logical *wantq, logical *wantz, integer *n, float *a, integer *lda, float *b, integer *ldb, float *q, integer *ldq, float * z, integer *ldz, integer *j1, integer *n1, integer *n2, float *work, integer *lwork, integer *info)
{
  return stgex2_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork, info);
}
inline integer tgex2(logical *wantq, logical *wantz, integer *n, double *a, integer *lda, double *b, integer *ldb, double *q, integer *ldq, double * z, integer *ldz, integer *j1, integer *n1, integer *n2, double *work, integer *lwork, integer *info)
{
  return dtgex2_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork, info);
}
inline integer tgex2(logical *wantq, logical *wantz, integer *n, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *q, integer *ldq, scomplex *z, integer *ldz, integer *j1, integer *info)
{
  return ctgex2_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, info);
}
inline integer tgex2(logical *wantq, logical *wantz, integer *n, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *q, integer *ldq, dcomplex *z, integer *ldz, integer *j1, integer *info)
{
  return ztgex2_(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, info);
}

// --- solves the generalized Sylvester equation (unblocked algorithm) ---
inline integer tgsy2(char *trans, integer *ijob, integer *m, integer * n, float *a, integer *lda, float *b, integer *ldb, float *c, integer * ldc, float *d, integer *ldd, float *e, integer *lde, float *f, integer *ldf, float *scale, float *rdsum, float *rdscal, integer *iwork, integer *pq, integer *info)
{
  return stgsy2_(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal, iwork, pq, info);
}
inline integer tgsy2(char *trans, integer *ijob, integer *m, integer * n, double *a, integer *lda, double *b, integer *ldb, double *c, integer * ldc, double *d, integer *ldd, double *e, integer *lde, double *f, integer *ldf, double *scale, double *rdsum, double *rdscal, integer *iwork, integer *pq, integer *info)
{
  return dtgsy2_(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal, iwork, pq, info);
}
inline integer tgsy2(char *trans, integer *ijob, integer *m, integer * n, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *c, integer *ldc, scomplex *d, integer *ldd, scomplex *e, integer *lde, scomplex *f, integer *ldf, float *scale, float *rdsum, float *rdscal, integer *info)
{
  return ctgsy2_(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal, info);
}
inline integer tgsy2(char *trans, integer *ijob, integer *m, integer * n, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *c, integer * ldc, dcomplex *d, integer *ldd, dcomplex *e, integer *lde, dcomplex *f, integer *ldf, double *scale, double *rdsum, double *rdscal, integer *info)
{
  return ztgsy2_(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal, info);
}

// --- computes a blocked LQ factorization of a float "triangular-pentagonal" matrix C ---
inline integer tplqt(integer *m, integer *n, integer *l, integer *mb, float *a, integer *lda, float *b, integer *ldb, float *t, integer *ldt, float *work, integer *info)
{
  return stplqt_(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info);
}
inline integer tplqt(integer *m, integer *n, integer *l, integer *mb, double *a, integer *lda, double *b, integer *ldb, double *t, integer *ldt, double *work, integer *info)
{
  return dtplqt_(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info);
}
inline integer tplqt(integer *m, integer *n, integer *l, integer *mb, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *t, integer *ldt, scomplex *work, integer *info)
{
  return ctplqt_(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info);
}
inline integer tplqt(integer *m, integer *n, integer *l, integer *mb, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *t, integer *ldt, dcomplex *work, integer *info)
{
  return ztplqt_(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info);
}

// --- computes a LQ factorization of a float or scomplex "triangular-pentagonal" matrix ---
inline integer tplqt2(integer *m, integer *n, integer *l, float *a, integer *lda, float *b, integer *ldb, float *t, integer *ldt, integer * info)
{
  return stplqt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline integer tplqt2(integer *m, integer *n, integer *l, double *a, integer *lda, double *b, integer *ldb, double *t, integer *ldt, integer * info)
{
  return dtplqt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline integer tplqt2(integer *m, integer *n, integer *l, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *t, integer *ldt, integer * info)
{
  return ctplqt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}
inline integer tplqt2(integer *m, integer *n, integer *l, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *t, integer *ldt, integer * info)
{
  return ztplqt2_(m, n, l, a, lda, b, ldb, t, ldt, info);
}

// --- applies a float orthogonal matrix Q obtained from a "triangular-pentagonal" float block reflector ---
inline integer tpmlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, integer *mb, float *v, integer *ldv, float *t, integer *ldt, float *a, integer *lda, float *b, integer *ldb, float * work, integer *info)
{
  return stpmlqt_(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline integer tpmlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, integer *mb, double *v, integer *ldv, double *t, integer *ldt, double *a, integer *lda, double *b, integer *ldb, double * work, integer *info)
{
  return dtpmlqt_(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline integer tpmlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, integer *mb, scomplex *v, integer *ldv, scomplex *t, integer *ldt, scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex * work, integer *info)
{
  return ctpmlqt_(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}
inline integer tpmlqt(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, integer *mb, dcomplex *v, integer *ldv, dcomplex *t, integer *ldt, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex * work, integer *info)
{
  return ztpmlqt_(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info);
}

// --- computes some or all of the right and/or left eigenvectors of a float upper quasi-triangular matrix T ---
inline integer trevc3(char *side, char *howmny, logical *select, integer *n, float *t, integer *ldt, float *vl, integer *ldvl, float *vr, integer *ldvr, integer *mm, integer *m, float *work, integer *lwork, integer *info)
{
  return strevc3_(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info);
}
inline integer trevc3(char *side, char *howmny, logical *select, integer *n, double *t, integer *ldt, double *vl, integer *ldvl, double *vr, integer *ldvr, integer *mm, integer *m, double *work, integer *lwork, integer *info)
{
  return dtrevc3_(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info);
}
inline integer trevc3(char *side, char *howmny, logical *select, integer *n, scomplex *t, integer *ldt, scomplex *vl, integer *ldvl, scomplex *vr, integer *ldvr, integer *mm, integer *m, scomplex *work, integer *lwork, float *rwork, integer *lrwork, integer *info)
{
  return ctrevc3_(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, rwork, lrwork, info);
}
inline integer trevc3(char *side, char *howmny, logical *select, integer *n, dcomplex *t, integer *ldt, dcomplex *vl, integer *ldvl, dcomplex *vr, integer *ldvr, integer *mm, integer *m, dcomplex *work, integer *lwork, double *rwork, integer *lrwork, integer *info)
{
  return ztrevc3_(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, rwork, lrwork, info);
}

// --- simultaneously bidiagonalizes the blocks of a tall and skinny matrix X with orthonomal columns ---
inline integer orbdb1(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21, integer *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *work, integer *lwork, integer *info)
{
  return sorbdb1_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline integer orbdb1(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21, integer *ldx21, double *theta, double *phi, double *taup1, double *taup2, double *tauq1, double *work, integer *lwork, integer *info)
{
  return dorbdb1_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline integer unbdb1(integer *m, integer *p, integer *q, scomplex * x11, integer *ldx11, scomplex *x21, integer *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *work, integer *lwork, integer *info)
{
  return cunbdb1_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline integer unbdb1(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21, integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info)
{
  return zunbdb1_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}

// --- simultaneously bidiagonalizes the blocks of a tall and skinny matrix X with orthonomal columns ---
inline integer orbdb2(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21, integer *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *work, integer *lwork, integer *info)
{
  return sorbdb2_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline integer orbdb2(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21, integer *ldx21, double *theta, double *phi, double *taup1, double *taup2, double *tauq1, double *work, integer *lwork, integer *info)
{
  return dorbdb2_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline integer unbdb2(integer *m, integer *p, integer *q, scomplex * x11, integer *ldx11, scomplex *x21, integer *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *work, integer *lwork, integer *info)
{
  return cunbdb2_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline integer unbdb2(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21, integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info)
{
  return zunbdb2_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}

// --- simultaneously bidiagonalizes the blocks of a tall and skinny matrix X with orthonomal columns ---
inline integer orbdb3(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21, integer *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *work, integer *lwork, integer *info)
{
  return sorbdb3_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline integer orbdb3(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21, integer *ldx21, double *theta, double *phi, double *taup1, double *taup2, double *tauq1, double *work, integer *lwork, integer *info)
{
  return dorbdb3_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline integer unbdb3(integer *m, integer *p, integer *q, scomplex * x11, integer *ldx11, scomplex *x21, integer *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *work, integer *lwork, integer *info)
{
  return cunbdb3_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}
inline integer unbdb3(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21, integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info)
{
  return zunbdb3_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info);
}

// --- simultaneously bidiagonalizes the blocks of a tall and skinny matrix X with orthonomal columns ---
inline integer orbdb4(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21, integer *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *phantom, float *work, integer *lwork, integer *info)
{
  return sorbdb4_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work, lwork, info);
}
inline integer orbdb4(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21, integer *ldx21, double *theta, double *phi, double *taup1, double *phantom, double *taup2, double *tauq1, double *work, integer *lwork, integer *info)
{
  return dorbdb4_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work, lwork, info);
}
inline integer unbdb4(integer *m, integer *p, integer *q, scomplex * x11, integer *ldx11, scomplex *x21, integer *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *phantom, scomplex *taup2, scomplex *tauq1, scomplex *work, integer *lwork, integer *info)
{
  return cunbdb4_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work, lwork, info);
}
inline integer unbdb4(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21, integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *phantom, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info)
{
  return zunbdb4_(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work, lwork, info);
}

// --- orthogonalizes the column vector ---
inline integer orbdb5(integer *m1, integer *m2, integer *n, float *x1, integer *incx1, float *x2, integer *incx2, float *q1, integer *ldq1, float *q2, integer *ldq2, float *work, integer *lwork, integer *info)
{
  return sorbdb5_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline integer orbdb5(integer *m1, integer *m2, integer *n, double *x1, integer *incx1, double *x2, integer *incx2, double *q1, integer *ldq1, double *q2, integer *ldq2, double *work, integer *lwork, integer *info)
{
  return dorbdb5_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline integer unbdb5(integer *m1, integer *m2, integer *n, scomplex *x1, integer *incx1, scomplex *x2, integer *incx2, scomplex *q1, integer *ldq1, scomplex *q2, integer *ldq2, scomplex *work, integer *lwork, integer *info)
{
  return cunbdb5_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline integer unbdb5(integer *m1, integer *m2, integer *n, dcomplex *x1, integer *incx1, dcomplex *x2, integer *incx2, dcomplex *q1, integer *ldq1, dcomplex *q2, integer *ldq2, dcomplex *work, integer *lwork, integer *info)
{
  return zunbdb5_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}

// --- orthogonalizes the column vector ---
inline integer orbdb6(integer *m1, integer *m2, integer *n, float *x1, integer *incx1, float *x2, integer *incx2, float *q1, integer *ldq1, float *q2, integer *ldq2, float *work, integer *lwork, integer *info)
{
  return sorbdb6_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline integer orbdb6(integer *m1, integer *m2, integer *n, double *x1, integer *incx1, double *x2, integer *incx2, double *q1, integer *ldq1, double *q2, integer *ldq2, double *work, integer *lwork, integer *info)
{
  return dorbdb6_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline integer unbdb6(integer *m1, integer *m2, integer *n, scomplex *x1, integer *incx1, scomplex *x2, integer *incx2, scomplex *q1, integer *ldq1, scomplex *q2, integer *ldq2, scomplex *work, integer *lwork, integer *info)
{
  return cunbdb6_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}
inline integer unbdb6(integer *m1, integer *m2, integer *n, dcomplex *x1, integer *incx1, dcomplex *x2, integer *incx2, dcomplex *q1, integer *ldq1, dcomplex *q2, integer *ldq2, dcomplex *work, integer *lwork, integer *info)
{
  return zunbdb6_(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info);
}

// --- generates all or part of the orthogonal matrix Q from a QL factorization determined by sgeqlf (unblocked algorithm) ---
inline integer org2l(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *work, integer *info)
{
  return sorg2l_(m, n, k, a, lda, tau, work, info);
}
inline integer org2l(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *work, integer *info)
{
  return dorg2l_(m, n, k, a, lda, tau, work, info);
}
inline integer ung2l(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  return cung2l_(m, n, k, a, lda, tau, work, info);
}
inline integer ung2l(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  return zung2l_(m, n, k, a, lda, tau, work, info);
}

// --- generates all or part of the orthogonal matrix Q from a QR factorization determined by sgeqrf ---
inline integer org2r(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *work, integer *info)
{
  return sorg2r_(m, n, k, a, lda, tau, work, info);
}
inline integer org2r(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *work, integer *info)
{
  return dorg2r_(m, n, k, a, lda, tau, work, info);
}
inline integer ung2r(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  return cung2r_(m, n, k, a, lda, tau, work, info);
}
inline integer ung2r(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  return zung2r_(m, n, k, a, lda, tau, work, info);
}

// --- generates an m by n float matrix Q with orthonormal rows ---
inline integer orgl2(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *work, integer *info)
{
  return sorgl2_(m, n, k, a, lda, tau, work, info);
}
inline integer orgl2(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *work, integer *info)
{
  return dorgl2_(m, n, k, a, lda, tau, work, info);
}
inline integer ungl2(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  return cungl2_(m, n, k, a, lda, tau, work, info);
}
inline integer ungl2(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  return zungl2_(m, n, k, a, lda, tau, work, info);
}

// --- generates all or part of the orthogonal matrix Q from an RQ factorization determined by sgerqf ---
inline integer orgr2(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *work, integer *info)
{
  return sorgr2_(m, n, k, a, lda, tau, work, info);
}
inline integer orgr2(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *work, integer *info)
{
  return dorgr2_(m, n, k, a, lda, tau, work, info);
}
inline integer ungr2(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
  return cungr2_(m, n, k, a, lda, tau, work, info);
}
inline integer ungr2(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work, integer *info)
{
  return zungr2_(m, n, k, a, lda, tau, work, info);
}

// --- generates an M-by-N float matrix Q_out with orthonormal columns ---
inline integer orgtsqr(integer *m, integer *n, integer *mb, integer * nb, float *a, integer *lda, float *t, integer *ldt, float *work, integer *lwork, integer *info)
{
  return sorgtsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline integer orgtsqr(integer *m, integer *n, integer *mb, integer * nb, double *a, integer *lda, double *t, integer *ldt, double *work, integer *lwork, integer *info)
{
  return dorgtsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline integer ungtsqr(integer *m, integer *n, integer *mb, integer * nb, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *work, integer *lwork, integer *info)
{
  return cungtsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline integer ungtsqr(integer *m, integer *n, integer *mb, integer * nb, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *work, integer *lwork, integer *info)
{
  return zungtsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}

// --- takes an M-by-N float matrix Q_in with orthonormal columns as input and performs Householder Reconstruction ---
inline integer orhr_col(integer *m, integer *n, integer *nb, float *a, integer *lda, float *t, integer *ldt, float *d, integer *info)
{
  return sorhr_col_(m, n, nb, a, lda, t, ldt, d, info);
}
inline integer orhr_col(integer *m, integer *n, integer *nb, double *a, integer *lda, double *t, integer *ldt, double *d, integer *info)
{
  return dorhr_col_(m, n, nb, a, lda, t, ldt, d, info);
}
inline integer unhr_col(integer *m, integer *n, integer *nb, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *d, integer *info)
{
  return cunhr_col_(m, n, nb, a, lda, t, ldt, d, info);
}
inline integer unhr_col(integer *m, integer *n, integer *nb, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *d, integer *info)
{
  return zunhr_col_(m, n, nb, a, lda, t, ldt, d, info);
}

// --- multiplies a general matrix by the orthogonal matrix from a QL factorization determined by sgeqlf ---
inline integer orm2l(char *side, char *trans, integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *c, integer *ldc, float *work, integer *info)
{
  return sorm2l_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline integer orm2l(char *side, char *trans, integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *c, integer *ldc, double *work, integer *info)
{
  return dorm2l_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline integer unm2l(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *c, integer *ldc, scomplex *work, integer *info)
{
  return cunm2l_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline integer unm2l(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work, integer *info)
{
  return zunm2l_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

// --- multiplies a general matrix by a banded orthogonal matrix ---
inline integer orm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, float *q, integer *ldq, float *c, integer * ldc, float *work, integer *lwork, integer *info)
{
  return sorm22_(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
}
inline integer orm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, double *q, integer *ldq, double *c, integer * ldc, double *work, integer *lwork, integer *info)
{
  return dorm22_(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
}
inline integer unm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, scomplex *q, integer *ldq, scomplex *c, integer * ldc, scomplex *work, integer *lwork, integer *info)
{
  return cunm22_(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
}
inline integer unm22(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, dcomplex *q, integer *ldq, dcomplex *c, integer * ldc, dcomplex *work, integer *lwork, integer *info)
{
  return zunm22_(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info);
}

// --- multiplies a general matrix by the orthogonal matrix from a RQ factorization determined by sgerqf ---
inline integer ormr2(char *side, char *trans, integer *m, integer *n, integer *k, float *a, integer *lda, float *tau, float *c, integer *ldc, float *work, integer *info)
{
  return sormr2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline integer ormr2(char *side, char *trans, integer *m, integer *n, integer *k, double *a, integer *lda, double *tau, double *c, integer *ldc, double *work, integer *info)
{
  return dormr2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline integer unmr2(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *c, integer *ldc, scomplex *work, integer *info)
{
  return cunmr2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}
inline integer unmr2(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work, integer *info)
{
  return zunmr2_(side, trans, m, n, k, a, lda, tau, c, ldc, work, info);
}

// --- multiplies a general matrix by the orthogonal matrix from a RZ factorization determined by stzrzf ---
inline integer ormr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, float *a, integer *lda, float *tau, float *c, integer *ldc, float *work, integer *info)
{
  return sormr3_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
}
inline integer ormr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, double *a, integer *lda, double *tau, double *c, integer *ldc, double *work, integer *info)
{
  return dormr3_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
}
inline integer unmr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, scomplex *a, integer *lda, scomplex *tau, scomplex *c, integer *ldc, scomplex *work, integer *info)
{
  return cunmr3_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
}
inline integer unmr3(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work, integer *info)
{
  return zunmr3_(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info);
}

// --- LA_GBAMV performs a matrix-vector operation to calculate error bounds ---
inline integer la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, float *alpha, float *ab, integer *ldab, float * x, integer *incx, float *beta, float *y, integer *incy)
{
  return sla_gbamv_(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
}
inline integer la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha, double *ab, integer *ldab, double * x, integer *incx, double *beta, double *y, integer *incy)
{
  return dla_gbamv_(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
}
inline integer la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, float *alpha, scomplex *ab, integer *ldab, scomplex * x, integer *incx, float *beta, float *y, integer *incy)
{
  return cla_gbamv_(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
}
inline integer la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha, dcomplex *ab, integer *ldab, dcomplex * x, integer *incx, double *beta, double *y, integer *incy)
{
  return zla_gbamv_(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
}

// --- estimates the Skeel condition number for a general banded matrix ---
inline float la_gbrcond(char *trans, integer *n, integer *kl, integer *ku, float *ab, integer *ldab, float *afb, integer *ldafb, integer *ipiv, integer * cmode, float *c, integer *info, float *work, integer *iwork)
{
  return sla_gbrcond_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, cmode, c, info, work, iwork);
}
inline double la_gbrcond(char *trans, integer *n, integer *kl, integer *ku, double * ab, integer *ldab, double *afb, integer *ldafb, integer *ipiv, integer * cmode, double *c, integer *info, double *work, integer *iwork)
{
  return dla_gbrcond_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, cmode, c, info, work, iwork);
}

// --- computes the infinity norm condition number of op(A)*inv(diag(c)) for general banded matrices ---
inline float la_gbrcond_c(char *trans, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv, float *c, logical *capply, integer *info, scomplex *work, float * rwork)
{
  return cla_gbrcond_c_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, capply, info, work, rwork);
}
inline double la_gbrcond_c(char *trans, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv, double *c, logical *capply, integer *info, dcomplex *work, double * rwork)
{
  return zla_gbrcond_c_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, capply, info, work, rwork);
}

// --- computes the infinity norm condition number of op(A)*diag(x) for general banded matrices ---
inline float la_gbrcond_x(char *trans, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_gbrcond_x_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, x, info, work, rwork);
}
inline double la_gbrcond_x(char *trans, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_gbrcond_x_(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, x, info, work, rwork);
}

// --- improves the computed solution to a system of linear equations for general banded matrices ---
inline integer la_gbrfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *kl, integer *ku, integer *nrhs, float *ab, integer *ldab, float *afb, integer *ldafb, integer *ipiv, logical *colequ, float *c, float *b, integer *ldb, float *y, integer * ldy, float *berr_out, integer *n_norms, float *err_bnds_norm, float *err_bnds_comp, float *res, float *ayb, float *dy, float *y_tail, float *rcond, integer *ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  return sla_gbrfsx_extended_(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_gbrfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *kl, integer *ku, integer *nrhs, double *ab, integer *ldab, double *afb, integer *ldafb, integer *ipiv, logical *colequ, double *c, double *b, integer *ldb, double *y, integer * ldy, double *berr_out, integer *n_norms, double *err_bnds_norm, double *err_bnds_comp, double *res, double *ayb, double *dy, double *y_tail, double *rcond, integer *ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  return dla_gbrfsx_extended_(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_gbrfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *kl, integer *ku, integer *nrhs, scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb, integer * ipiv, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex * y, integer *ldy, float *berr_out, integer *n_norms, float * err_bnds_norm, float *err_bnds_comp, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer *ithresh, float * rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  return cla_gbrfsx_extended_(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_gbrfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *kl, integer *ku, integer *nrhs, dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb, integer * ipiv, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex * y, integer *ldy, double *berr_out, integer *n_norms, double * err_bnds_norm, double *err_bnds_comp, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer *ithresh, double * rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  return zla_gbrfsx_extended_(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- improves the computed solution to a system of linear equations for general matrices  ---
inline integer la_gerfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *nrhs, float *a, integer *lda, float * af, integer *ldaf, integer *ipiv, logical *colequ, float *c, float *b, integer *ldb, float *y, integer *ldy, float *berr_out, integer * n_norms, float *errs_n, float *errs_c, float *res, float *ayb, float *dy, float *y_tail, float *rcond, integer *ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  return sla_gerfsx_extended_(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_gerfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *nrhs, double *a, integer *lda, double * af, integer *ldaf, integer *ipiv, logical *colequ, double *c, double *b, integer *ldb, double *y, integer *ldy, double *berr_out, integer * n_norms, double *errs_n, double *errs_c, double *res, double *ayb, double *dy, double *y_tail, double *rcond, integer *ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  return dla_gerfsx_extended_(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_gerfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex *y, integer *ldy, float *berr_out, integer *n_norms, float *errs_n, float *errs_c, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer * ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  return cla_gerfsx_extended_(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_gerfsx_extended(integer *prec_type, integer * trans_type, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex *y, integer *ldy, double *berr_out, integer *n_norms, double *errs_n, double *errs_c, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer * ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  return zla_gerfsx_extended_(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- improves the computed solution to a system of linear equations for symmetric or Hermitian positive-definite matrices ---
inline integer la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, float *a, integer *lda, float *af, integer * ldaf, logical *colequ, float *c, float *b, integer *ldb, float *y, integer *ldy, float *berr_out, integer *n_norms, float * err_bnds_norm, float *err_bnds_comp, float *res, float *ayb, float *dy, float *y_tail, float *rcond, integer *ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  return sla_porfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, double *a, integer *lda, double *af, integer * ldaf, logical *colequ, double *c, double *b, integer *ldb, double *y, integer *ldy, double *berr_out, integer *n_norms, double * err_bnds_norm, double *err_bnds_comp, double *res, double *ayb, double * dy, double *y_tail, double *rcond, integer *ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  return dla_porfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af, integer *ldaf, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex *y, integer *ldy, float *berr_out, integer *n_norms, float * err_bnds_norm, float *err_bnds_comp, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer *ithresh, float * rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  return cla_porfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_porfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex *y, integer *ldy, double *berr_out, integer *n_norms, double * err_bnds_norm, double *err_bnds_comp, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer *ithresh, double * rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  return zla_porfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- improves the computed solution to a system of linear equations for symmetric indefinite matrices ---
inline integer la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, float *a, integer *lda, float *af, integer * ldaf, integer *ipiv, logical *colequ, float *c, float *b, integer * ldb, float *y, integer *ldy, float *berr_out, integer *n_norms, float *err_bnds_norm, float *err_bnds_comp, float *res, float *ayb, float *dy, float *y_tail, float *rcond, integer *ithresh, float * rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  return sla_syrfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, double *a, integer *lda, double *af, integer * ldaf, integer *ipiv, logical *colequ, double *c, double *b, integer * ldb, double *y, integer *ldy, double *berr_out, integer *n_norms, double *err_bnds_norm, double *err_bnds_comp, double *res, double *ayb, double *dy, double *y_tail, double *rcond, integer *ithresh, double * rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  return dla_syrfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex *y, integer *ldy, float *berr_out, integer * n_norms, float *err_bnds_norm, float *err_bnds_comp, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer * ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  return cla_syrfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_syrfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex *y, integer *ldy, double *berr_out, integer * n_norms, double *err_bnds_norm, double *err_bnds_comp, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer * ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  return zla_syrfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- improves the computed solution to a system of linear equations for Hermitian indefinite matrices ---
inline integer la_herfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, logical *colequ, float *c, scomplex *b, integer *ldb, scomplex *y, integer *ldy, float *berr_out, integer * n_norms, float *err_bnds_norm, float *err_bnds_comp, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail, float *rcond, integer * ithresh, float *rthresh, float *dz_ub, logical *ignore_cwise, integer *info)
{
  return cla_herfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}
inline integer la_herfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, logical *colequ, double *c, dcomplex *b, integer *ldb, dcomplex *y, integer *ldy, double *berr_out, integer * n_norms, double *err_bnds_norm, double *err_bnds_comp, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail, double *rcond, integer * ithresh, double *rthresh, double *dz_ub, logical *ignore_cwise, integer *info)
{
  return zla_herfsx_extended_(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) for a general banded matrix ---
inline float la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, float *ab, integer *ldab, float *afb, integer *ldafb)
{
  return sla_gbrpvgrw_(n, kl, ku, ncols, ab, ldab, afb, ldafb);
}
inline double la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, double *ab, integer *ldab, double *afb, integer *ldafb)
{
  return dla_gbrpvgrw_(n, kl, ku, ncols, ab, ldab, afb, ldafb);
}
inline float la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb)
{
  return cla_gbrpvgrw_(n, kl, ku, ncols, ab, ldab, afb, ldafb);
}
inline double la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb)
{
  return zla_gbrpvgrw_(n, kl, ku, ncols, ab, ldab, afb, ldafb);
}

// --- computes a matrix-vector product using a general matrix to calculate error bounds ---
inline integer la_geamv(integer *trans, integer *m, integer *n, float *alpha, float *a, integer *lda, float *x, integer *incx, float * beta, float *y, integer *incy)
{
  return sla_geamv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline integer la_geamv(integer *trans, integer *m, integer *n, double *alpha, double *a, integer *lda, double *x, integer *incx, double * beta, double *y, integer *incy)
{
  return dla_geamv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline integer la_geamv(integer *trans, integer *m, integer *n, float *alpha, scomplex *a, integer *lda, scomplex *x, integer *incx, float * beta, float *y, integer *incy)
{
  return cla_geamv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline integer la_geamv(integer *trans, integer *m, integer *n, double *alpha, dcomplex *a, integer *lda, dcomplex *x, integer *incx, double * beta, double *y, integer *incy)
{
  return zla_geamv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

// --- estimates the Skeel condition number for a general matrix ---
inline float la_gercond(char *trans, integer *n, float *a, integer *lda, float *af, integer *ldaf, integer *ipiv, integer *cmode, float *c, integer * info, float *work, integer *iwork)
{
  return sla_gercond_(trans, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork); 
}
inline  double la_gercond(char *trans, integer *n, double *a, integer *lda, double *af, integer *ldaf, integer *ipiv, integer *cmode, double *c, integer * info, double *work, integer *iwork)
{
  return dla_gercond_(trans, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork); 
}
inline float la_gercond_x(char *trans, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_gercond_x_(trans, n, a, lda, af, ldaf, ipiv, x, info, work, rwork); 
}
inline double la_gercond_x(char *trans, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_gercond_x_(trans, n, a, lda, af, ldaf, ipiv, x, info, work, rwork); 
}
inline float la_gercond_c(char *trans, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, float *c, logical *capply, integer *info, scomplex *work, float *rwork)
{
  return cla_gercond_c_(trans, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork); 
}
inline double la_gercond_c(char *trans, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, double *c, logical *capply, integer *info, dcomplex *work, double *rwork)
{
  return zla_gercond_c_(trans, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork); 
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) ---
inline float la_gerpvgrw(integer *n, integer *ncols, float *a, integer *lda, float * af, integer *ldaf)
{
  return sla_gerpvgrw_(n, ncols, a, lda, af, ldaf);
}
inline double la_gerpvgrw(integer *n, integer *ncols, double *a, integer *lda, double * af, integer *ldaf)
{
  return dla_gerpvgrw_(n, ncols, a, lda, af, ldaf);
}
inline float la_gerpvgrw(integer *n, integer *ncols, scomplex *a, integer *lda, scomplex *af, integer *ldaf)
{
  return cla_gerpvgrw_(n, ncols, a, lda, af, ldaf);
}
inline double la_gerpvgrw(integer *n, integer *ncols, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf)
{
  return zla_gerpvgrw_(n, ncols, a, lda, af, ldaf);
}

// --- computes a matrix-vector product using a Hermitian indefinite matrix to calculate error bounds ---
inline integer la_heamv(integer *uplo, integer *n, float *alpha, scomplex *a, integer *lda, scomplex *x, integer *incx, float *beta, float *y, integer *incy)
{
  return cla_heamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline integer la_heamv(integer *uplo, integer *n, double *alpha, dcomplex *a, integer *lda, dcomplex *x, integer *incx, double *beta, double *y, integer *incy)
{
  return zla_heamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

// --- computes the infinity norm condition number of op(A)*inv(diag(c)) for Hermitian indefinite matrices ---
inline float la_hercond_c(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, float *c, logical *capply, integer *info, scomplex *work, float *rwork)
{
  return cla_hercond_c_(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork);
}
inline double la_hercond_c(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, double *c, logical *capply, integer *info, dcomplex *work, double *rwork)
{
  return zla_hercond_c_(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork);
}

// --- computes the infinity norm condition number of op(A)*diag(x) for Hermitian indefinite matrices ---
inline float la_hercond_x(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_hercond_x_(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork);
}
inline double la_hercond_x(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_hercond_x_(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork);
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) ---
inline float la_herpvgrw(char *uplo, integer *n, integer *info, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, float *work)
{
  return cla_herpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work); 
}
inline double la_herpvgrw(char *uplo, integer *n, integer *info, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, double *work)
{
  return zla_herpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work); 
}

// --- computes a component-wise relative backward error ---
inline integer la_lin_berr(integer *n, integer *nz, integer *nrhs, float *res, float *ayb, float *berr)
{
  return sla_lin_berr_(n, nz, nrhs, res, ayb, berr);
}
inline integer la_lin_berr(integer *n, integer *nz, integer *nrhs, double *res, double *ayb, double *berr)
{
  return dla_lin_berr_(n, nz, nrhs, res, ayb, berr);
}
inline integer la_lin_berr(integer *n, integer *nz, integer *nrhs, scomplex *res, float *ayb, float *berr)
{
  return cla_lin_berr_(n, nz, nrhs, res, ayb, berr);
}
inline integer la_lin_berr(integer *n, integer *nz, integer *nrhs, dcomplex *res, double *ayb, double *berr)
{
  return zla_lin_berr_(n, nz, nrhs, res, ayb, berr);
}

// --- estimates the Skeel condition number for a symmetric positive-definite matrix ---
inline float la_porcond(char *uplo, integer *n, float *a, integer *lda, float *af, integer *ldaf, integer *cmode, float *c, integer *info, float *work, integer *iwork)
{
  return sla_porcond_(uplo, n, a, lda, af, ldaf, cmode, c, info, work, iwork); 
}
inline double la_porcond(char *uplo, integer *n, double *a, integer *lda, double *af, integer *ldaf, integer *cmode, double *c, integer *info, double *work, integer *iwork)
{
  return dla_porcond_(uplo, n, a, lda, af, ldaf, cmode, c, info, work, iwork); 
}

// --- computes the infinity norm condition number of op(A)*diag(x) for Hermitian positive-definite matrices ---
inline float la_porcond_x(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_porcond_x_(uplo, n, a, lda, af, ldaf, x, info, work, rwork); 
}
inline double la_porcond_x(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_porcond_x_(uplo, n, a, lda, af, ldaf, x, info, work, rwork); 
}

// --- computes the infinity norm condition number of op(A)*inv(diag(c)) for Hermitian positive-definite matrices ---
inline float la_porcond_c(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, float *c, logical *capply, integer *info, scomplex *work, float *rwork)
{
  return cla_porcond_c_(uplo, n, a, lda, af, ldaf, c, capply, info, work, rwork); 
}
inline double la_porcond_c(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, double *c, logical *capply, integer *info, dcomplex *work, double *rwork)
{
  return zla_porcond_c_(uplo, n, a, lda, af, ldaf, c, capply, info, work, rwork); 
}

// --- computes a matrix-vector product using a symmetric indefinite matrix to calculate error bounds ---
inline integer la_syamv(integer *uplo, integer *n, float *alpha, float *a, integer *lda, float *x, integer *incx, float *beta, float *y, integer *incy)
{
  return sla_syamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline integer la_syamv(integer *uplo, integer *n, double *alpha, double *a, integer *lda, double *x, integer *incx, double *beta, double *y, integer *incy)
{
  return dla_syamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline integer la_syamv(integer *uplo, integer *n, float *alpha, scomplex *a, integer *lda, scomplex *x, integer *incx, float *beta, float *y, integer *incy)
{
  return cla_syamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}
inline integer la_syamv(integer *uplo, integer *n, double *alpha, dcomplex *a, integer *lda, dcomplex *x, integer *incx, double *beta, double *y, integer *incy)
{
  return zla_syamv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric or Hermitian positive-definite matrix ---
inline float la_porpvgrw(char *uplo, integer *ncols, float *a, integer *lda, float * af, integer *ldaf, float *work)
{
  return sla_porpvgrw_(uplo, ncols, a, lda, af, ldaf, work);
}
inline double la_porpvgrw(char *uplo, integer *ncols, double *a, integer *lda, double * af, integer *ldaf, double *work)
{
  return dla_porpvgrw_(uplo, ncols, a, lda, af, ldaf, work);
}
inline float la_porpvgrw(char *uplo, integer *ncols, scomplex *a, integer *lda, scomplex * af, integer *ldaf, float *work)
{
  return cla_porpvgrw_(uplo, ncols, a, lda, af, ldaf, work);
}
inline double la_porpvgrw(char *uplo, integer *ncols, dcomplex *a, integer *lda, dcomplex * af, integer *ldaf, double *work)
{
  return zla_porpvgrw_(uplo, ncols, a, lda, af, ldaf, work);
}

// --- estimates the Skeel condition number for a symmetric indefinite matrix ---
inline float la_syrcond(char *uplo, integer *n, float *a, integer *lda, float *af, integer *ldaf, integer *ipiv, integer *cmode, float *c, integer * info, float *work, integer *iwork)
{
  return sla_syrcond_(uplo, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork); 
}
inline double la_syrcond(char *uplo, integer *n, double *a, integer *lda, double *af, integer *ldaf, integer *ipiv, integer *cmode, double *c, integer * info, double *work, integer *iwork)
{
  return dla_syrcond_(uplo, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork); 
}

// --- computes the infinity norm condition number of op(A)*diag(x) for symmetric indefinite matrices ---
inline float la_syrcond_x(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, scomplex *x, integer *info, scomplex *work, float *rwork)
{
  return cla_syrcond_x_(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork); 
}
inline double la_syrcond_x(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, dcomplex *x, integer *info, dcomplex *work, double *rwork)
{
  return zla_syrcond_x_(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork); 
}

// --- LA_SYRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for symmetric indefinite matrices ---
inline float la_syrcond_c(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv, float *c, logical *capply, integer *info, scomplex *work, float *rwork)
{
  return cla_syrcond_c_(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork); 
}
inline double la_syrcond_c(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv, double *c, logical *capply, integer *info, dcomplex *work, double *rwork)
{
  return zla_syrcond_c_(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork); 
}

// --- computes the reciprocal pivot growth factor norm(A)/norm(U) for a symmetric indefinite matrix ---
inline float la_syrpvgrw(char *uplo, integer *n, integer *info, float *a, integer * lda, float *af, integer *ldaf, integer *ipiv, float *work)
{
  return sla_syrpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work);
}
inline double la_syrpvgrw(char *uplo, integer *n, integer *info, double *a, integer * lda, double *af, integer *ldaf, integer *ipiv, double *work)
{
  return dla_syrpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work);
}
inline float la_syrpvgrw(char *uplo, integer *n, integer *info, scomplex *a, integer * lda, scomplex *af, integer *ldaf, integer *ipiv, float *work)
{
  return cla_syrpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work);
}
inline double la_syrpvgrw(char *uplo, integer *n, integer *info, dcomplex *a, integer * lda, dcomplex *af, integer *ldaf, integer *ipiv, double *work)
{
  return zla_syrpvgrw_(uplo, n, info, a, lda, af, ldaf, ipiv, work);
}

// --- adds a vector into a doubled-single vector ---
inline integer la_wwaddw(integer *n, float *x, float *y, float *w)
{
  return sla_wwaddw_(n, x, y, w);
}
inline integer la_wwaddw(integer *n, double *x, double *y, double *w)
{
  return dla_wwaddw_(n, x, y, w);
}
inline integer la_wwaddw(integer *n, scomplex *x, scomplex *y, scomplex *w)
{
  return cla_wwaddw_(n, x, y, w);
}
inline integer la_wwaddw(integer *n, dcomplex *x, dcomplex *y, dcomplex *w)
{
  return zla_wwaddw_(n, x, y, w);
}

// --- takes as input the values computed by LAMCH for underflow and overflow, and returns the square root of each of these values ---
inline integer labad(float *small, float *large)
{
  return slabad_( small, large);
}
inline integer labad(double *small, double *large)
{
  return dlabad_( small, large);
}

// --- reduces the first nb rows and columns of a general matrix to a bidiagonal form ---
inline integer labrd(integer *m, integer *n, integer *nb, float *a, integer *lda, float *d, float *e, float *tauq, float *taup, float *x, integer *ldx, float *y, integer *ldy)
{
  return slabrd_(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
}
inline integer labrd(integer *m, integer *n, integer *nb, double *a, integer *lda, double *d, double *e, double *tauq, double *taup, double *x, integer *ldx, double *y, integer *ldy)
{
  return dlabrd_(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
}
inline integer labrd(integer *m, integer *n, integer *nb, scomplex *a, integer *lda, float *d, float *e, scomplex *tauq, scomplex *taup, scomplex *x, integer *ldx, scomplex *y, integer *ldy)
{
  return clabrd_(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
}
inline integer labrd(integer *m, integer *n, integer *nb, dcomplex *a, integer *lda, double *d, double *e, dcomplex *tauq, dcomplex *taup, dcomplex *x, integer *ldx, dcomplex *y, integer *ldy)
{
  return zlabrd_(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
}

// --- estimates the 1-norm of a square matrix, using reverse communication for evaluating matrix-vector products ---
inline integer lacon(integer *n, float *v, float *x, integer *isgn, float *est, integer *kase)
{
  return slacon_(n, v, x, isgn, est, kase);
}
inline integer lacon(integer *n, double *v, double *x, integer *isgn, double *est, integer *kase)
{
  return dlacon_(n, v, x, isgn, est, kase);
}
inline integer lacon(integer *n, scomplex *v, scomplex *x, float *est, integer *kase)
{
  return clacon_(n, v, x, est, kase);
}
inline integer lacon(integer *n, dcomplex *v, dcomplex *x, double *est, integer *kase)
{
  return zlacon_(n, v, x, est, kase);
}

// --- performs a linear transformation of a pair of scomplex vectors ---
inline integer lacrt(integer *n, scomplex *cx, integer *incx, scomplex * cy, integer *incy, scomplex *c, scomplex *s)
{
  return clacrt_(n, cx, incx, cy, incy, c, s); 
}
inline integer lacrt(integer *n, dcomplex *cx, integer *incx, dcomplex * cy, integer *incy, dcomplex *c, dcomplex *s)
{
  return zlacrt_(n, cx, incx, cy, incy, c, s); 
}

// --- performs scomplex division in float arithmetic, avoiding unnecessary overflow ---
inline integer ladiv(float *a, float *b, float *c, float *d, float *p, float *q)
{
  return sladiv_(a, b, c, d, p, q);
}
inline integer ladiv(double *a, double *b, double *c, double *d, double *p, double *q)
{
  return dladiv_(a, b, c, d, p, q);
}

// --- computes the eigenvalues of a 2-by-2 symmetric matrix ---
inline integer lae2(float *a, float *b, float *c, float *rt1, float *rt2)
{
  return slae2_(a, b, c, rt1, rt2);
}
inline integer lae2(double *a, double *b, double *c, double *rt1, double *rt2)
{
  return dlae2_(a, b, c, rt1, rt2);
}

// --- computes the number of eigenvalues of a float symmetric tridiagonal matrix ---
inline integer laebz(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp, integer *nbmin, float *abstol, float * reltol, float *pivmin, float *d, float *e, float *e2, integer *nval, float *ab, float *c, integer *mout, integer *nab, float *work, integer *iwork, integer *info)
{
  return slaebz_(ijob, nitmax, n, mmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval, ab, c, mout, nab, work, iwork, info);
}
inline integer laebz(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp, integer *nbmin, double *abstol, double * reltol, double *pivmin, double *d, double *e, double *e2, integer *nval, double *ab, double *c, integer *mout, integer *nab, double *work, integer *iwork, integer *info)
{
  return dlaebz_(ijob, nitmax, n, mmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval, ab, c, mout, nab, work, iwork, info);
}

// --- computes all eigenvalues and corresponding eigenvectors of an unreduced symmetric tridiagonal matrix using the divide and conquer method ---
inline integer laed0(integer *icompq, integer *qsiz, integer *n, float *d, float *e, float *q, integer *ldq, float *qstore, integer *ldqs, float *work, integer *iwork, integer *info)
{
  return slaed0_(icompq, qsiz, n, d, e, q, ldq, qstore, ldqs, work, iwork, info);
}
inline integer laed0(integer *icompq, integer *qsiz, integer *n, double *d, double *e, double *q, integer *ldq, double *qstore, integer *ldqs, double *work, integer *iwork, integer *info)
{
  return dlaed0_(icompq, qsiz, n, d, e, q, ldq, qstore, ldqs, work, iwork, info);
}
inline integer laed0(integer *qsiz, integer *n, float *d, float *e, scomplex *q, integer *ldq, scomplex *qstore, integer *ldqs, float *rwork, integer *iwork, integer *info)
{
  return claed0_(qsiz, n, d, e, q, ldq, qstore, ldqs, rwork, iwork, info);
}
inline integer laed0(integer *qsiz, integer *n, double *d, double *e, dcomplex *q, integer *ldq, dcomplex *qstore, integer *ldqs, double *rwork, integer *iwork, integer *info)
{
  return zlaed0_(qsiz, n, d, e, q, ldq, qstore, ldqs, rwork, iwork, info);
}

// --- computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is tridiagonal ---
inline integer laed1(integer *n, float *d, float *q, integer *ldq, integer *indxq, float *rho, integer *cutpnt, float *work, integer * iwork, integer *info)
{
  return slaed1_(n, d, q, ldq, indxq, rho, cutpnt, work, iwork, info); 
}
inline integer laed1(integer *n, double *d, double *q, integer *ldq, integer *indxq, double *rho, integer *cutpnt, double *work, integer * iwork, integer *info)
{
  return dlaed1_(n, d, q, ldq, indxq, rho, cutpnt, work, iwork, info); 
}

// --- Merges eigenvalues and deflates secular equation. Used when the original matrix is tridiagonal ---
inline integer laed2(integer *k, integer *n, integer *n1, float *d, float *q, integer *ldq, integer *indxq, float *rho, float *z, float * dlamda, float *w, float *q2, integer *indx, integer *indxc, integer * indxp, integer *coltyp, integer *info)
{
  return slaed2_(k, n, n1, d, q, ldq, indxq, rho, z, dlamda, w, q2, indx, indxc, indxp, coltyp, info); 
}
inline integer laed2(integer *k, integer *n, integer *n1, double *d, double *q, integer *ldq, integer *indxq, double *rho, double *z, double * dlamda, double *w, double *q2, integer *indx, integer *indxc, integer * indxp, integer *coltyp, integer *info)
{
  return dlaed2_(k, n, n1, d, q, ldq, indxq, rho, z, dlamda, w, q2, indx, indxc, indxp, coltyp, info); 
}

// --- Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is tridiagonal ---
inline integer laed3(integer *k, integer *n, integer *n1, float *d, float *q, integer *ldq, float *rho, float *dlamda, float *q2, integer * indx, integer *ctot, float *w, float *s, integer *info)
{
  return slaed3_(k, n, n1, d, q, ldq, rho, dlamda, q2, indx, ctot, w, s, info); 
}
inline integer laed3(integer *k, integer *n, integer *n1, double *d, double *q, integer *ldq, double *rho, double *dlamda, double *q2, integer * indx, integer *ctot, double *w, double *s, integer *info)
{
  return dlaed3_(k, n, n1, d, q, ldq, rho, dlamda, q2, indx, ctot, w, s, info); 
}

// --- Finds a single root of the secular equation ---
inline integer laed4(integer *n, integer *i, float *d, float *z, float *delta, float *rho, float *dlam, integer *info)
{
  return slaed4_(n, i, d, z, delta, rho, dlam, info); 
}
inline integer laed4(integer *n, integer *i, double *d, double *z, double *delta, double *rho, double *dlam, integer *info)
{
  return dlaed4_(n, i, d, z, delta, rho, dlam, info); 
}

// --- Solves the 2-by-2 secular equation ---
inline integer laed5(integer *i, float *d, float *z, float *delta, float *rho, float *dlam)
{
  return slaed5_(i, d, z, delta, rho, dlam);
}
inline integer laed5(integer *i, double *d, double *z, double *delta, double *rho, double *dlam)
{
  return dlaed5_(i, d, z, delta, rho, dlam);
}

// --- computes one Newton step in solution of the secular equation ---
inline integer laed6(integer *kniter, logical *orgati, float *rho, float *d, float *z, float *finit, float *tau, integer *info)
{
  return slaed6_(kniter, orgati, rho, d, z, finit, tau, info);  
}
inline integer laed6(integer *kniter, logical *orgati, double *rho, double *d, double *z, double *finit, double *tau, integer *info)
{
  return dlaed6_(kniter, orgati, rho, d, z, finit, tau, info);  
}

// --- computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. used when the original matrix is dense ---
inline integer laed7(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, float *d, float *q, integer *ldq, integer *indxq, float *rho, integer *cutpnt, float * qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, float *givnum, float *work, integer *iwork, integer *info)
{
  return slaed7_(icompq, n, qsiz, tlvls, curlvl, curpbm, d, q, ldq, indxq,  rho, cutpnt, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, iwork, info); 
}
inline integer laed7(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, double *d, double *q, integer *ldq, integer *indxq, double *rho, integer *cutpnt, double * qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, double *givnum, double *work, integer *iwork, integer *info)
{
  return dlaed7_(icompq, n, qsiz, tlvls, curlvl, curpbm, d, q, ldq, indxq, rho, cutpnt, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, iwork, info); 
}
inline integer laed7(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, float *d, scomplex *q, integer *ldq, float *rho, integer *indxq, float *qstore, integer * qptr, integer *prmptr, integer *perm, integer *givptr, integer * givcol, float *givnum, scomplex *work, float *rwork, integer *iwork, integer *info)
{
  return claed7_(n, cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, ldq, rho, indxq, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, info); 
}
inline integer laed7(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl, integer *curpbm, double *d, dcomplex *q, integer *ldq, double *rho, integer *indxq, double *qstore, integer * qptr, integer *prmptr, integer *perm, integer *givptr, integer *givcol, double *givnum, dcomplex *work, double *rwork, integer *iwork, integer *info)
{
  return zlaed7_(n, cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, ldq, rho, indxq, qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, info); 
}

// --- Merges eigenvalues and deflates secular equation. Used when the original matrix is dense ---
inline integer laed8(integer *icompq, integer *k, integer *n, integer *qsiz, float *d, float *q, integer *ldq, integer *indxq, float *rho, integer *cutpnt, float *z, float *dlamda, float *q2, integer *ldq2, float *w, integer *perm, integer *givptr, integer *givcol, float * givnum, integer *indxp, integer *indx, integer *info)
{
  return slaed8_(icompq, k, n, qsiz, d, q, ldq, indxq, rho, cutpnt, z, dlamda, q2, ldq2, w, perm, givptr, givcol, givnum, indxp, indx, info);
}
inline integer laed8(integer *icompq, integer *k, integer *n, integer *qsiz, double *d, double *q, integer *ldq, integer *indxq, double *rho, integer *cutpnt, double *z, double *dlamda, double *q2, integer *ldq2, double *w, integer *perm, integer *givptr, integer *givcol, double * givnum, integer *indxp, integer *indx, integer *info)
{
  return dlaed8_(icompq, k, n, qsiz, d, q, ldq, indxq, rho, cutpnt, z, dlamda, q2, ldq2, w, perm, givptr, givcol, givnum, indxp, indx, info);  
}
inline integer laed8(integer *k, integer *n, integer *qsiz, scomplex * q, integer *ldq, float *d, float *rho, integer *cutpnt, float *z, float *dlamda, scomplex *q2, integer *ldq2, float *w, integer *indxp, integer *indx, integer *indxq, integer *perm, integer *givptr, integer *givcol, float *givnum, integer *info)
{
  return claed8_(k, n, qsiz, q, ldq, d, rho, cutpnt, z, dlamda, q2, ldq2, w, indxp, indx, indxq, perm, givptr, givcol, givnum, info);
}
inline integer laed8(integer *k, integer *n, integer *qsiz, dcomplex * q, integer *ldq, double *d, double *rho, integer *cutpnt, double *z, double *dlamda, dcomplex *q2, integer *ldq2, double *w, integer *indxp, integer *indx, integer *indxq, integer *perm, integer *givptr, integer *givcol, double *givnum, integer *info)
{
  return zlaed8_(k, n, qsiz, q, ldq, d, rho, cutpnt, z, dlamda, q2, ldq2, w, indxp, indx, indxq, perm, givptr, givcol, givnum, info);  
}

// --- Finds the roots of the secular equation and updates the eigenvectors. Used when the original matrix is dense ---
inline integer laed9(integer *k, integer *kstart, integer *kstop, integer *n, float *d, float *q, integer *ldq, float *rho, float *dlamda, float *w, float *s, integer *lds, integer *info)
{
  return slaed9_(k, kstart, kstop, n, d, q, ldq, rho, dlamda, w, s, lds, info); 
}
inline integer laed9(integer *k, integer *kstart, integer *kstop, integer *n, double *d, double *q, integer *ldq, double *rho, double *dlamda, double *w, double *s, integer *lds, integer *info)
{
  return dlaed9_(k, kstart, kstop, n, d, q, ldq, rho, dlamda, w, s, lds, info); 
}

// --- Computes the Z vector determining the rank-one modification of the diagonal matrix. Used when the original matrix is dense ---
inline integer laeda(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr, integer *perm, integer *givptr, integer *givcol, float *givnum, float *q, integer *qptr, float *z, float *ztemp, integer *info)
{
  return slaeda_(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z, ztemp, info); 
}
inline integer laeda(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr, integer *perm, integer *givptr, integer *givcol, double *givnum, double *q, integer *qptr, double *z, double *ztemp, integer *info)
{
  return dlaeda_(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z, ztemp, info); 
}

// --- computes a specified right or left eigenvector of an upper Hessenberg matrix by inverse iteration ---
inline integer laein(logical *rightv, logical *noinit, integer *n, float *h, integer *ldh, float *wr, float *wi, float *vr, float *vi, float *b, integer *ldb, float *work, float *eps3, float *smlnum, float *bignum, integer *info)
{
  return slaein_(rightv, noinit, n, h, ldh, wr, wi, vr, vi, b, ldb, work, eps3, smlnum, bignum, info); 
}
inline integer laein(logical *rightv, logical *noinit, integer *n, double *h, integer *ldh, double *wr, double *wi, double *vr, double *vi, double *b, integer *ldb, double *work, double *eps3, double *smlnum, double *bignum, integer *info)
{
  return dlaein_(rightv, noinit, n, h, ldh, wr, wi, vr, vi, b, ldb, work, eps3, smlnum, bignum, info); 
}
inline integer laein(logical *rightv, logical *noinit, integer *n, scomplex *h, integer *ldh, scomplex *w, scomplex *v, scomplex *b, integer *ldb, float *rwork, float *eps3, float *smlnum, integer *info)
{
  return claein_(rightv, noinit, n, h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info); 
}
inline integer laein(logical *rightv, logical *noinit, integer *n, dcomplex *h, integer *ldh, dcomplex *w, dcomplex *v, dcomplex *b, integer *ldb, double *rwork, double *eps3, double *smlnum, integer *info)
{
  return zlaein_(rightv, noinit, n, h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info); 
}

// --- computes the eigenvalues and eigenvectors of a 2-by-2 scomplex symmetric matrix ---
inline integer laesy(scomplex *a, scomplex *b, scomplex *c, scomplex * rt1, scomplex *rt2, scomplex *evscal, scomplex *cs1, scomplex *sn1)
{
  return claesy_(a, b, c, rt1, rt2, evscal, cs1, sn1); 
}
inline integer laesy(dcomplex *a, dcomplex *b, dcomplex *c, dcomplex * rt1, dcomplex *rt2, dcomplex *evscal, dcomplex *cs1, dcomplex *sn1)
{
  return zlaesy_(a, b, c, rt1, rt2, evscal, cs1, sn1); 
}

// --- LAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix ---
inline integer laev2(float *a, float *b, float *c, float *rt1, float * rt2, float *cs1, float *sn1)
{
  return slaev2_(a, b, c, rt1, rt2, cs1, sn1);
}
inline integer laev2(double *a, double *b, double *c, double *rt1, double * rt2, double *cs1, double *sn1)
{
  return dlaev2_(a, b, c, rt1, rt2, cs1, sn1);
}
inline integer laev2(scomplex *a, scomplex *b, scomplex *c, float *rt1, float * rt2, float *cs1, scomplex *sn1)
{
  return claev2_(a, b, c, rt1, rt2, cs1, sn1);
}
inline integer laev2(dcomplex *a, dcomplex *b, dcomplex *c, double *rt1, double * rt2, double *cs1, dcomplex *sn1)
{
  return zlaev2_(a, b, c, rt1, rt2, cs1, sn1);
}

// --- swaps adjacent diagonal blocks of a float upper quasi-triangular matrix in Schur canonical form, by an orthogonal similarity transformation  ---
inline integer laexc(logical *wantq, integer *n, float *t, integer * ldt, float *q, integer *ldq, integer *j1, integer *n1, integer *n2, float *work, integer *info)
{
  return slaexc_(wantq, n, t, ldt, q, ldq, j1, n1, n2, work, info); 
}
inline integer laexc(logical *wantq, integer *n, double *t, integer * ldt, double *q, integer *ldq, integer *j1, integer *n1, integer *n2, double *work, integer *info)
{
  return dlaexc_(wantq, n, t, ldt, q, ldq, j1, n1, n2, work, info); 
}

// --- computes the eigenvalues of a 2-by-2 generalized eigenvalue problem, with scaling as necessary to avoid over-/underflow ---
inline integer lag2(float *a, integer *lda, float *b, integer *ldb, float *safmin, float *scale1, float *scale2, float *wr1, float *wr2, float * wi)
{
  return slag2_( a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi);
}
inline integer lag2(double *a, integer *lda, double *b, integer *ldb, double *safmin, double *scale1, double *scale2, double *wr1, double *wr2, double * wi)
{
  return dlag2_( a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi);
}

// --- computes 2-by-2 orthogonal matrices U, V, and Q, and applies them to matrices A and B such that the rows of the transformed A and B are parallel ---
inline integer lags2(logical *upper, float *a1, float *a2, float *a3, float *b1, float *b2, float *b3, float *csu, float *snu, float *csv, float * snv, float *csq, float *snq)
{
  return slags2_(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
}
inline integer lags2(logical *upper, double *a1, double *a2, double *a3, double *b1, double *b2, double *b3, double *csu, double *snu, double *csv, double * snv, double *csq, double *snq)
{
  return dlags2_(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
}
inline integer lags2(logical *upper, float *a1, scomplex *a2, float *a3, float *b1, scomplex *b2, float *b3, float *csu, scomplex *snu, float *csv, scomplex *snv, float *csq, scomplex *snq)
{
  return clags2_(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
}
inline integer lags2(logical *upper, double *a1, dcomplex *a2, double *a3, double *b1, dcomplex *b2, double *b3, double *csu, dcomplex *snu, double *csv, dcomplex *snv, double *csq, dcomplex *snq)
{
  return zlags2_(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
}

// --- computes an LU factorization of a matrix T-Î»I, where T is a general tridiagonal matrix, and Î» a scalar, using partial pivoting with row interchanges ---
inline integer lagtf(integer *n, float *a, float *lambda, float *b, float *c, float *tol, float *d, integer *in, integer *info)
{
  return slagtf_(n, a, lambda, b, c, tol, d, in, info);
}
inline integer lagtf(integer *n, double *a, double *lambda, double *b, double *c, double *tol, double *d, integer *in, integer *info)
{
  return dlagtf_(n, a, lambda, b, c, tol, d, in, info);
}

// --- performs a matrix-matrix product of the form C = Î±AB+Î²C, where A is a tridiagonal matrix, B and C are rectangular matrices, and Î± and Î² are scalars ---
inline integer lagtm(char *trans, integer *n, integer *nrhs, float *alpha, float *dl, float *d, float *du, float *x, integer *ldx, float * beta, float *b, integer *ldb)
{
  return  slagtm_(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
}
inline integer lagtm(char *trans, integer *n, integer *nrhs, double *alpha, double *dl, double *d, double *du, double *x, integer *ldx, double *beta, double *b, integer *ldb)
{
  return  dlagtm_(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
}
inline integer lagtm(char *trans, integer *n, integer *nrhs, float *alpha, scomplex *dl, scomplex *d, scomplex *du, scomplex *x, integer *ldx, float *beta, scomplex *b, integer *ldb)
{
  return  clagtm_(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
}
inline integer lagtm(char *trans, integer *n, integer *nrhs, double *alpha, dcomplex *dl, dcomplex *d, dcomplex *du, dcomplex *x, integer *ldx, double *beta, dcomplex *b, integer *ldb)
{
  return  zlagtm_(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
}

// --- solves the system of equations (T-Î»I)x = y or (T-Î»I)Tx = y ---
inline integer lagts(integer *job, integer *n, float *a, float *b, float *c, float *d, integer *in, float *y, float *tol, integer *info)
{
  return slagts_(job, n, a, b, c, d, in, y, tol, info);
}
inline integer lagts(integer *job, integer *n, double *a, double *b, double *c, double *d, integer *in, double *y, double *tol, integer *info)
{
  return dlagts_(job, n, a, b, c, d, in, y, tol, info);
}

// --- computes the Generalized Schur factorization of a float 2-by-2 matrix pencil (A,B) where B is upper triangular ---
inline integer lagv2(float *a, integer *lda, float *b, integer *ldb, float *alphar, float *alphai, float *beta, float *csl, float *snl, float * csr, float *snr)
{
  return slagv2_(a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr);
}
inline integer lagv2(double *a, integer *lda, double *b, integer *ldb, double *alphar, double *alphai, double *beta, double *csl, double *snl, double * csr, double *snr)
{
  return dlagv2_(a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr);
}

// --- computes a partial factorization of a scomplex Hermitian indefinite matrix using the Bunch-Kaufman diagonal pivoting method ---
inline integer lahef(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, integer *ipiv, scomplex *w, integer *ldw, integer *info)
{
  return clahef_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline integer lahef(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *w, integer *ldw, integer *info)
{
  return zlahef_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}

// ---  factorizes a panel of a scomplex hermitian matrix A using the Aasen's algorithm ---
inline integer lahef_aa(char *uplo, integer *j1, integer *m, integer *nb, scomplex *a, integer *lda, integer *ipiv, scomplex *h, integer * ldh, scomplex *work) 
{
  return clahef_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work); 
}
inline integer lahef_aa(char *uplo, integer *j1, integer *m, integer *nb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *h, integer * ldh, dcomplex *work) 
{
  return zlahef_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work); 
}

// --- computes a partial factorization of a scomplex Hermitian indefinite matrix using bounded Bunch-Kaufman (rook) diagonal pivoting method ---
inline integer lahef_rk(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, scomplex *e, integer *ipiv, scomplex *w, integer *ldw, integer *info)
{
  return clahef_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}
inline integer lahef_rk(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, dcomplex *w, integer *ldw, integer *info)
{
  return zlahef_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}

// --- computes a partial factorization of a scomplex Hermitian matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting method ---
inline integer lahef_rook(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, integer *ipiv, scomplex *w, integer *ldw, integer *info)
{
  return clahef_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info); 
}
inline integer lahef_rook(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *w, integer *ldw, integer *info)
{
  return zlahef_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info); 
}

// --- computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm ---
inline integer lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, float *h, integer *ldh, float *wr, float * wi, integer *iloz, integer *ihiz, float *z, integer *ldz, integer * info)
{
  return slahqr_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info);
}
inline integer lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, double *h, integer *ldh, double *wr, double * wi, integer *iloz, integer *ihiz, double *z, integer *ldz, integer * info)
{
  return dlahqr_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info);
}
inline integer lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, scomplex *h, integer *ldh, scomplex *w, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, integer * info)
{
  return clahqr_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
}
inline integer lahqr(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, dcomplex *h, integer *ldh, dcomplex *w, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, integer * info)
{
  return zlahqr_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
}

// --- reduces the specified number of first columns of a general rectangular matrix A ---
inline integer lahr2(integer *n, integer *k, integer *nb, float *a, integer *lda, float *tau, float *t, integer *ldt, float *y, integer *ldy)
{
  return slahr2_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline integer lahr2(integer *n, integer *k, integer *nb, double *a, integer *lda, double *tau, double *t, integer *ldt, double *y, integer *ldy)
{
  return dlahr2_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline integer lahr2(integer *n, integer *k, integer *nb, scomplex *a, integer *lda, scomplex *tau, scomplex *t, integer *ldt, scomplex *y, integer *ldy)
{
  return clahr2_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline integer lahr2(integer *n, integer *k, integer *nb, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *t, integer *ldt, dcomplex *y, integer *ldy)
{
  return zlahr2_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}

// --- reduces the first nb columns of a general rectangular matrix A ---
inline integer lahrd(integer *n, integer *k, integer *nb, float *a, integer *lda, float *tau, float *t, integer *ldt, float *y, integer *ldy)
{
  printf(" Function slahrd() has been deprecated. Please use slahr2_() instead.\n"); 
  return slahrd_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline integer lahrd(integer *n, integer *k, integer *nb, double *a, integer *lda, double *tau, double *t, integer *ldt, double *y, integer *ldy)
{
  printf(" Function dlahrd() has been deprecated. Please use dlahr2_() instead.\n"); 
  return dlahrd_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline integer lahrd(integer *n, integer *k, integer *nb, scomplex *a, integer *lda, scomplex *tau, scomplex *t, integer *ldt, scomplex *y, integer *ldy)
{
  printf(" Function clahrd() has been deprecated. Please use clahr2_() instead.\n"); 
  return clahrd_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}
inline integer lahrd(integer *n, integer *k, integer *nb, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *t, integer *ldt, dcomplex *y, integer *ldy)
{
  printf(" Function zlahrd() has been deprecated. Please use zlahr2_() instead.\n"); 
  return zlahrd_(n, k, nb, a, lda, tau, t, ldt, y, ldy);
}

// --- applies one step of incremental condition estimation ---
inline integer laic1(integer *job, integer *j, float *x, float *sest, float *w, float *gamma, float *sestpr, float *s, float *c__)
{
  return slaic1_(job, j, x, sest, w, gamma, sestpr, s, c__);
}
inline integer laic1(integer *job, integer *j, double *x, double *sest, double *w, double *gamma, double *sestpr, double *s, double *c__)
{
  return dlaic1_(job, j, x, sest, w, gamma, sestpr, s, c__);
}
inline integer laic1(integer *job, integer *j, scomplex *x, float *sest, scomplex *w, scomplex *gamma, float *sestpr, scomplex *s, scomplex *c__)
{
  return claic1_(job, j, x, sest, w, gamma, sestpr, s, c__);
}
inline integer laic1(integer *job, integer *j, dcomplex *x, double *sest, dcomplex *w, dcomplex *gamma, double *sestpr, dcomplex *s, dcomplex *c__)
{
  return zlaic1_(job, j, x, sest, w, gamma, sestpr, s, c__);
}

// --- tests input for NaN by comparing two arguments for inequality  ---
inline logical laisnan(float *sin1, float *sin2)
{
  return slaisnan_(sin1, sin2);
}
inline logical laisnan(double *sin1, double *sin2)
{
  return dlaisnan_(sin1, sin2);
}

// --- solves a 1-by-1 or 2-by-2 linear system of equations of the specified form ---
inline integer laln2(logical *ltrans, integer *na, integer *nw, float * smin, float *ca, float *a, integer *lda, float *d1, float *d2, float *b, integer *ldb, float *wr, float *wi, float *x, integer *ldx, float *scale, float *xnorm, integer *info)
{
  return slaln2_(ltrans, na, nw, smin, ca, a, lda, d1, d2, b, ldb, wr, wi, x, ldx, scale, xnorm, info);
}
inline integer laln2(logical *ltrans, integer *na, integer *nw, double * smin, double *ca, double *a, integer *lda, double *d1, double *d2, double *b, integer *ldb, double *wr, double *wi, double *x, integer *ldx, double *scale, double *xnorm, integer *info)
{
  return dlaln2_(ltrans, na, nw, smin, ca, a, lda, d1, d2, b, ldb, wr, wi, x, ldx, scale, xnorm, info);
}

// ---  applies back multiplying factors in solving the least squares problem using divide and conquer SVD approach. Used by sgelsd ---
inline integer lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, float *b, integer *ldb, float *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, float *givnum, integer *ldgnum, float *poles, float *difl, float *difr, float *z, integer *k, float *c, float *s, float *work, integer *info)
{
  return slals0_(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info);
}
inline integer lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, double *b, integer *ldb, double *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, double *givnum, integer *ldgnum, double *poles, double *difl, double *difr, double *z, integer *k, double *c, double *s, double *work, integer *info)
{
  return dlals0_(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info);
}
inline integer lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, scomplex *b, integer *ldb, scomplex *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, float *givnum, integer *ldgnum, float *poles, float * difl, float *difr, float *z, integer *k, float *c, float *s, float * work, integer *info)
{
  return clals0_(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info);
}
inline integer lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, dcomplex *b, integer *ldb, dcomplex *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, double *givnum, integer *ldgnum, double *poles, double * difl, double *difr, double *z, integer *k, double *c, double *s, double * work, integer *info)
{
  return zlals0_(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info);
}

// --- computes the SVD of the coefficient matrix in compact form. Used by sgelsd ---
inline integer lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, float *b, integer *ldb, float *bx, integer *ldbx, float *u, integer *ldu, float *vt, integer *k, float *difl, float *difr, float *z, float *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, float *givnum, float *c, float *s, float *work, integer * iwork, integer *info)
{
  return slalsa_(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
}
inline integer lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, double *b, integer *ldb, double *bx, integer *ldbx, double * u, integer *ldu, double *vt, integer *k, double *difl, double *difr, double * z, double *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, double *givnum, double *c, double *s, double *work, integer * iwork, integer *info)
{
  return dlalsa_(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
}
inline integer lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, scomplex *b, integer *ldb, scomplex *bx, integer *ldbx, float *u, integer *ldu, float *vt, integer *k, float *difl, float *difr, float * z, float *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, float *givnum, float *c, float *s, float *work, integer * iwork, integer *info)
{
  return clalsa_(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
}
inline integer lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, dcomplex *b, integer *ldb, dcomplex *bx, integer *ldbx, double * u, integer *ldu, double *vt, integer *k, double *difl, double *difr, double * z, double *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm, double *givnum, double *c, double *s, double *work, integer * iwork, integer *info)
{
  return zlalsa_(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
}

// --- uses the singular value decomposition of A to solve the least squares problem ---
inline integer lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, float *d, float *e, float *b, integer *ldb, float *rcond, integer *rank, float *work, integer *iwork, integer *info)
{
  return slalsd_(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, iwork, info);
}
inline integer lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, double *d, double *e, double *b, integer *ldb, double *rcond, integer *rank, double *work, integer *iwork, integer *info)
{
  return dlalsd_(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, iwork, info);
}
inline integer lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, float *d, float *e, scomplex *b, integer *ldb, float *rcond, integer *rank, scomplex *work, float* rwork, integer *iwork, integer *info)
{
  return clalsd_(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, rwork, iwork, info);
}
inline integer lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, double *d, double *e, dcomplex *b, integer *ldb, double *rcond, integer *rank, dcomplex *work, double* rwork, integer *iwork, integer *info)
{
  return zlalsd_(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, rwork, iwork, info);
}

// --- creates a permutation list to merge the entries of two independently sorted sets into a single set sorted in ascending order ---
inline integer lamrg(integer *n1, integer *n2, float *a, integer * strd1, integer *strd2, integer *index)
{
  return slamrg_(n1, n2, a, strd1, strd2, index);
}
inline integer lamrg(integer *n1, integer *n2, double *a, integer * strd1, integer *strd2, integer *index)
{
  return dlamrg_(n1, n2, a, strd1, strd2, index);
}

// --- overwrites the general float M-by-N matrix C ---
inline integer lamswlq(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, float *a, integer *lda, float *t, integer *ldt, float *c, integer *ldc, float *work, integer *lwork, integer *info)
{
  return slamswlq_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline integer lamswlq(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, double *a, integer *lda, double * t, integer *ldt, double *c, integer *ldc, double *work, integer *lwork, integer *info)
{
  return dlamswlq_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline integer lamswlq(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, scomplex *a, integer *lda, scomplex * t, integer *ldt, scomplex *c, integer *ldc, scomplex *work, integer *lwork, integer *info)
{
  return clamswlq_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline integer lamswlq(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, dcomplex *a, integer *lda, dcomplex * t, integer *ldt, dcomplex *c, integer *ldc, dcomplex *work, integer *lwork, integer *info)
{
  return zlamswlq_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}

// --- overwrites the general float M-by-N matrix C ---
inline integer lamtsqr(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, float *a, integer *lda, float *t, integer *ldt, float *c, integer *ldc, float *work, integer *lwork, integer *info)
{
  return slamtsqr_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline integer lamtsqr(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, double *a, integer *lda, double * t, integer *ldt, double *c, integer *ldc, double *work, integer *lwork, integer *info)
{
  return dlamtsqr_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline integer lamtsqr(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, scomplex *a, integer *lda, scomplex * t, integer *ldt, scomplex *c, integer *ldc, scomplex *work, integer *lwork, integer *info)
{
  return clamtsqr_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}
inline integer lamtsqr(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, dcomplex *a, integer *lda, dcomplex * t, integer *ldt, dcomplex *c, integer *ldc, dcomplex *work, integer *lwork, integer *info)
{
  return zlamtsqr_(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info);
}

// --- computes the Sturm count ---
inline integer laneg(integer *n, float *d, float *lld, float *sigma, float *pivmin, integer *r__)
{
  return slaneg_(n, d, lld, sigma, pivmin, r__);
}
inline integer laneg(integer *n, double *d, double *lld, double *sigma, double *pivmin, integer *r__)
{
  return dlaneg_(n, d, lld, sigma, pivmin, r__);
}

// --- returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of general band matrix ---
inline float langb(char *norm, integer *n, integer *kl, integer *ku, float *ab, integer *ldab, float *work)
{
  return slangb_(norm, n, kl, ku, ab, ldab, work);
}
inline double langb(char *norm, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, double *work)
{
  return dlangb_(norm, n, kl, ku, ab, ldab, work);
}
inline float langb(char *norm, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, float *work)
{
  return clangb_(norm, n, kl, ku, ab, ldab, work);
}
inline double langb(char *norm, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, double *work)
{
  return zlangb_(norm, n, kl, ku, ab, ldab, work);
}

// --- returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general tridiagonal matrix ---
inline float langt(char *norm, integer *n, float *dl, float *d, float *du)
{
  return slangt_(norm, n, dl, d, du);
}
inline double langt(char *norm, integer *n, double *dl, double *d, double *du)
{
  return dlangt_(norm, n, dl, d, du);
}
inline float langt(char *norm, integer *n, scomplex *dl, scomplex *d, scomplex *du)
{
  return clangt_(norm, n, dl, d, du);
}
inline double langt(char *norm, integer *n, dcomplex *dl, dcomplex *d, dcomplex *du)
{
  return zlangt_(norm, n, dl, d, du);
}

// --- returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a Hermitian band matrix  ---
inline float lanhb(char *norm, char *uplo, integer *n, integer *k, scomplex *ab, integer *ldab, float *work)
{
  return clanhb_(norm, uplo, n, k, ab, ldab, work); 
}
inline double lanhb(char *norm, char *uplo, integer *n, integer *k, dcomplex *ab, integer *ldab, double *work)
{
  return zlanhb_(norm, uplo, n, k, ab, ldab, work); 
}

//--- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a Hermitian matrix in RFP format ---
inline float lanhf(char *norm, char *transr, char *uplo, integer *n, scomplex *a, float *work)
{
  return clanhf_(norm, transr, uplo, n, a, work); 
}
inline double lanhf(char *norm, char *transr, char *uplo, integer *n, dcomplex *a, double *work)
{
  return zlanhf_(norm, transr, uplo, n, a, work); 
}

//--- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a scomplex Hermitian matrix supplied in packed form  ---
inline float lanhp(char *norm, char *uplo, integer *n, scomplex *ap, float *work)
{
  return clanhp_(norm, uplo, n, ap, work); 
}
inline double lanhp(char *norm, char *uplo, integer *n, dcomplex *ap, double *work)
{
  return zlanhp_(norm, uplo, n, ap, work); 
}

// --- returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of an upper Hessenberg matrix ---
inline float lanhs(char *norm, integer *n, float *a, integer *lda, float *work)
{
  return slanhs_(norm, n, a, lda, work);
}
inline double lanhs(char *norm, integer *n, double *a, integer *lda, double *work)
{
  return dlanhs_(norm, n, a, lda, work);
}
inline float lanhs(char *norm, integer *n, scomplex *a, integer *lda, float *work)
{
  return clanhs_(norm, n, a, lda, work);
}
inline double lanhs(char *norm, integer *n, dcomplex *a, integer *lda, double *work)
{
  return zlanhs_(norm, n, a, lda, work);
}

// --- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a scomplex Hermitian tridiagonal matrix  ---
inline float lanht(char *norm, integer *n, float *d, scomplex *e)
{
  return clanht_(norm, n, d, e);
}
inline double lanht(char *norm, integer *n, double *d, dcomplex *e)
{
  return zlanht_(norm, n, d, e);
}

// --- returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a symmetric band matrix  ---
inline float lansb(char *norm, char *uplo, integer *n, integer *k, float *ab, integer *ldab, float *work)
{
  return slansb_(norm, uplo, n, k, ab, ldab, work);
}
inline integer lansb(char *norm, char *uplo, integer *n, integer *k, double *ab, integer *ldab, double *work)
{
  return dlansb_(norm, uplo, n, k, ab, ldab, work);
}
inline float lansb(char *norm, char *uplo, integer *n, integer *k, scomplex *ab, integer *ldab, float *work)
{
  return clansb_(norm, uplo, n, k, ab, ldab, work);
}
inline double lansb(char *norm, char *uplo, integer *n, integer *k, dcomplex *ab, integer *ldab, double *work)
{
  return zlansb_(norm, uplo, n, k, ab, ldab, work);
}

// --- returns the value of the one norm, or the infinity norm, or the element of largest absolute value of a float symmetric matrix A in RFP format  ---
inline float lansf(char *norm, char *transr, char *uplo, integer *n, float *a, float * work)
{
  return slansf_(norm, transr, uplo, n, a, work); 
}
inline double lansf(char *norm, char *transr, char *uplo, integer *n, double *a, double * work)
{
  return dlansf_(norm, transr, uplo, n, a, work); 
}

// --- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a symmetric matrix supplied in packed form ---
inline float lansp(char *norm, char *uplo, integer *n, float *ap, float *work)
{
  return slansp_(norm, uplo, n, ap, work);
}
inline double lansp(char *norm, char *uplo, integer *n, double *ap, double *work)
{
  return dlansp_(norm, uplo, n, ap, work);
}
inline float lansp(char *norm, char *uplo, integer *n, scomplex *ap, float *work)
{
  return clansp_(norm, uplo, n, ap, work);
}
inline double lansp(char *norm, char *uplo, integer *n, dcomplex *ap, double *work)
{
  return zlansp_(norm, uplo, n, ap, work);
}

// --- returns the value of the 1-norm, or the infinity norm, or the element of largest absolute value of a float symmetric tridiagonal matrix  ---
inline float lanst(char *norm, integer *n, float *d, float *e)
{
  return slanst_(norm, n, d, e);
}
inline double lanst(char *norm, integer *n, double *d, double *e)
{
  return dlanst_(norm, n, d, e);
}

// --- computes the Schur factorization of a float 2-by-2 nonsymmetric matrix in standard form  ---
inline integer lanv2(float *a, float *b, float *c, float *d, float * rt1r, float *rt1i, float *rt2r, float *rt2i, float *cs, float *sn)
{
  return slanv2_(a, b, c, d,  rt1r, rt1i, rt2r, rt2i, cs, sn);
}
inline integer lanv2(double *a, double *b, double *c, double *d, double * rt1r, double *rt1i, double *rt2r, double *rt2i, double *cs, double *sn)
{
  return dlanv2_(a, b, c, d,  rt1r, rt1i, rt2r, rt2i, cs, sn);
}

// --- computes the modified LU factorization without pivoting of a float general M-by-N matrix A ---
inline integer laorhr_col_getrfnp(integer *m, integer *n, float *a, integer *lda, float *d, integer *info)
{
  return slaorhr_col_getrfnp_(m, n, a, lda, d, info);
}
inline integer laorhr_col_getrfnp(integer *m, integer *n, double *a, integer *lda, double *d, integer *info)
{
  return dlaorhr_col_getrfnp_(m, n, a, lda, d, info);
}

// --- computes the modified LU factorization without pivoting of a float general M-by-N matrix A ---
inline integer launhr_col_getrfnp(integer *m, integer *n, scomplex *a, integer *lda, scomplex *d, integer *info)
{
  return claunhr_col_getrfnp_(m, n, a, lda, d, info);
}
inline integer launhr_col_getrfnp(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *d, integer *info)
{
  return zlaunhr_col_getrfnp_(m, n, a, lda, d, info);
}

// --- computes the modified LU factorization without pivoting of a float general M-by-N matrix A ---
inline integer laorhr_col_getrfnp2(integer *m, integer *n, float *a, integer *lda, float *d, integer *info)
{
  return slaorhr_col_getrfnp2_(m, n, a, lda, d, info);
}
inline integer laorhr_col_getrfnp2(integer *m, integer *n, double *a, integer *lda, double *d, integer *info)
{
  return dlaorhr_col_getrfnp2_(m, n, a, lda, d, info);
}

// --- computes the modified LU factorization without pivoting of a float general M-by-N matrix A ---
inline integer launhr_col_getrfnp2(integer *m, integer *n, scomplex *a, integer *lda, scomplex *d, integer *info)
{
  return claunhr_col_getrfnp2_(m, n, a, lda, d, info);
}
inline integer launhr_col_getrfnp2(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *d, integer *info)
{
  return zlaunhr_col_getrfnp2_(m, n, a, lda, d, info);
}

// --- measures the linear dependence of two vectors ---
inline integer lapll(integer *n, float *x, integer *incx, float *y, integer *incy, float *ssmin)
{
  return slapll_(n, x, incx, y, incy, ssmin);
}
inline integer lapll(integer *n, double *x, integer *incx, double *y, integer *incy, double *ssmin)
{
  return dlapll_(n, x, incx, y, incy, ssmin);
}
inline integer lapll(integer *n, scomplex *x, integer *incx, scomplex *y, integer *incy, float *ssmin)
{
  return clapll_(n, x, incx, y, incy, ssmin);
}
inline integer lapll(integer *n, dcomplex *x, integer *incx, dcomplex *y, integer *incy, double *ssmin)
{
  return zlapll_(n, x, incx, y, incy, ssmin);
}

// --- scales a general band matrix, using row and column scaling factors computed by sgbequ  ---
inline integer laqgb(integer *m, integer *n, integer *kl, integer *ku, float *ab, integer *ldab, float *r, float *c, float *rowcnd, float *colcnd, float *amax, char *equed)
{
  return slaqgb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
}
inline integer laqgb(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, char *equed)
{
  return dlaqgb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
}
inline integer laqgb(integer *m, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab, float *r, float *c, float *rowcnd, float *colcnd, float *amax, char *equed)
{
  return claqgb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
}
inline integer laqgb(integer *m, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax, char *equed)
{
  return zlaqgb_(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
}

// --- scales a general rectangular matrix, using row and column scaling factors computed by sgeequ ---
inline integer laqge(integer *m, integer *n, float *a, integer *lda, float *r, float *c, float *rowcnd, float *colcnd, float *amax, char * equed)
{
  return slaqge_(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
}
inline integer laqge(integer *m, integer *n, double *a, integer *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, char * equed)
{
  return dlaqge_(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
}
inline integer laqge(integer *m, integer *n, scomplex *a, integer *lda, float *r, float *c, float *rowcnd, float *colcnd, float *amax, char * equed)
{
  return claqge_(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
}
inline integer laqge(integer *m, integer *n, dcomplex *a, integer *lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax, char * equed)
{
  return zlaqge_(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
}

// --- scales a Hermitian band matrix, using scaling factors computed by cpbequ ---
inline integer laqhb(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *s, float *scond, float *amax, char *equed)
{
  return claqhb_(uplo, n, kd, ab, ldab, s, scond, amax, equed); 
}
inline integer laqhb(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *s, double *scond, double *amax, char *equed)
{
  return zlaqhb_(uplo, n, kd, ab, ldab, s, scond, amax, equed); 
}

// --- scales a Hermitian matrix  ---
inline integer laqhe(char *uplo, integer *n, scomplex *a, integer *lda, float *s, float *scond, float *amax, char *equed)
{
  return claqhe_(uplo, n, a, lda, s, scond, amax, equed); 
}
inline integer laqhe(char *uplo, integer *n, dcomplex *a, integer *lda, double *s, double *scond, double *amax, char *equed)
{
  return zlaqhe_(uplo, n, a, lda, s, scond, amax, equed); 
}

// --- scales a Hermitian matrix stored in packed form ---
inline integer laqhp(char *uplo, integer *n, scomplex *ap, float *s, float *scond, float *amax, char *equed)
{
  return claqhp_(uplo, n, ap, s, scond, amax, equed); 
}
inline integer laqhp(char *uplo, integer *n, dcomplex *ap, double *s, double *scond, double *amax, char *equed)
{
  return zlaqhp_(uplo, n, ap, s, scond, amax, equed); 
}

// --- computes a QR factorization with column pivoting of the matrix block ---
inline integer laqp2(integer *m, integer *n, integer *offset, float *a, integer *lda, integer *jpvt, float *tau, float *vn1, float *vn2, float *work)
{
  return slaqp2_(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
}
inline integer laqp2(integer *m, integer *n, integer *offset, double *a, integer *lda, integer *jpvt, double *tau, double *vn1, double *vn2, double *work)
{
  return dlaqp2_(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
}
inline integer laqp2(integer *m, integer *n, integer *offset, scomplex *a, integer *lda, integer *jpvt, scomplex *tau, float *vn1, float *vn2, scomplex *work)
{
  return claqp2_(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
}
inline integer laqp2(integer *m, integer *n, integer *offset, dcomplex *a, integer *lda, integer *jpvt, dcomplex *tau, double *vn1, double *vn2, dcomplex *work)
{
  return zlaqp2_(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work);
}

// --- computes a step of QR factorization with column pivoting of a float m-by-n matrix A by using BLAS level 3 ---
inline integer laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, float *a, integer *lda, integer *jpvt, float *tau, float *vn1, float *vn2, float *auxv, float *f, integer *ldf)
{
  return slaqps_(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
}
inline integer laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, double *a, integer *lda, integer *jpvt, double *tau, double *vn1, double *vn2, double *auxv, double *f, integer *ldf)
{
  return dlaqps_(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
}
inline integer laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, scomplex *a, integer *lda, integer *jpvt, scomplex *tau, float *vn1, float *vn2, scomplex *auxv, scomplex *f, integer *ldf)
{
  return claqps_(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
}
inline integer laqps(integer *m, integer *n, integer *offset, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *jpvt, dcomplex *tau, double *vn1, double *vn2, dcomplex *auxv, dcomplex *f, integer *ldf)
{
  return zlaqps_(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf);
}

// --- computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition ---
inline integer laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, float *h, integer *ldh, float *wr, float *wi, integer *iloz, integer *ihiz, float *z, integer *ldz, float *work, integer *lwork, integer *info)
{
  return slaqr0_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
}
inline integer laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, double *h, integer *ldh, double *wr, double *wi, integer *iloz, integer *ihiz, double *z, integer *ldz, double *work, integer *lwork, integer *info)
{
  return dlaqr0_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
}
inline integer laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, scomplex *h, integer *ldh, scomplex *w, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, scomplex * work, integer *lwork, integer *info)
{
  return claqr0_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
}
inline integer laqr0(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, dcomplex *h, integer *ldh, dcomplex *w, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, dcomplex * work, integer *lwork, integer *info)
{
  return zlaqr0_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
}

// --- sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts  ---
inline integer laqr1(integer *n, float *h, integer *ldh, float *sr1, float *si1, float *sr2, float *si2, float *v)
{
  return slaqr1_(n, h, ldh, sr1, si1, sr2, si2, v);
}
inline integer laqr1(integer *n, double *h, integer *ldh, double *sr1, double *si1, double *sr2, double *si2, double *v)
{
  return dlaqr1_(n, h, ldh, sr1, si1, sr2, si2, v);
}
inline integer laqr1(integer *n, scomplex *h, integer *ldh, scomplex * s1, scomplex *s2, scomplex *v)
{
  return claqr1_(n, h, ldh, s1, s2, v);
}
inline integer laqr1(integer *n, dcomplex *h, integer *ldh, dcomplex * s1, dcomplex *s2, dcomplex *v)
{
  return zlaqr1_(n, h, ldh, s1, s2, v);
}

// --- performs the orthogonal similarity transformation of a Hessenberg matrix ---
inline integer laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, float *h, integer *ldh, integer *iloz, integer *ihiz, float *z, integer *ldz, integer *ns, integer *nd, float *sr, float *si, float *v, integer *ldv, integer *nh, float *t, integer *ldt, integer *nv, float *wv, integer *ldwv, float *work, integer *lwork)
{
  return slaqr2_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline integer laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, double *h, integer *ldh, integer *iloz, integer *ihiz, double *z, integer *ldz, integer *ns, integer *nd, double *sr, double *si, double *v, integer *ldv, integer *nh, double *t, integer *ldt, integer *nv, double *wv, integer *ldwv, double *work, integer *lwork)
{
  return dlaqr2_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline integer laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, scomplex *h, integer *ldh, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, integer *ns, integer *nd, scomplex *sh, scomplex *v, integer *ldv, integer *nh, scomplex *t, integer *ldt, integer *nv, scomplex *wv, integer *ldwv, scomplex *work, integer *lwork)
{
  return claqr2_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline integer laqr2(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, dcomplex *h, integer *ldh, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, integer *ns, integer *nd, dcomplex *sh, dcomplex *v, integer *ldv, integer *nh, dcomplex *t, integer *ldt, integer *nv, dcomplex *wv, integer *ldwv, dcomplex *work, integer *lwork)
{
  return zlaqr2_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}

// --- performs the orthogonal similarity transformation of a Hessenberg matrix ---
inline integer laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, float *h, integer *ldh, integer *iloz, integer *ihiz, float *z, integer *ldz, integer *ns, integer *nd, float *sr, float *si, float *v, integer *ldv, integer *nh, float *t, integer *ldt, integer *nv, float *wv, integer *ldwv, float *work, integer *lwork)
{
  return slaqr3_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline integer laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, double *h, integer *ldh, integer *iloz, integer *ihiz, double *z, integer *ldz, integer *ns, integer *nd, double *sr, double *si, double *v, integer *ldv, integer *nh, double *t, integer *ldt, integer *nv, double *wv, integer *ldwv, double *work, integer *lwork)
{
  return dlaqr3_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline integer laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, scomplex *h, integer *ldh, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, integer *ns, integer *nd, scomplex *sh, scomplex *v, integer *ldv, integer *nh, scomplex *t, integer *ldt, integer *nv, scomplex *wv, integer *ldwv, scomplex *work, integer *lwork)
{
  return claqr3_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}
inline integer laqr3(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw, dcomplex *h, integer *ldh, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, integer *ns, integer *nd, dcomplex *sh, dcomplex *v, integer *ldv, integer *nh, dcomplex *t, integer *ldt, integer *nv, dcomplex *wv, integer *ldwv, dcomplex *work, integer *lwork)
{
  return zlaqr3_(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork);
}

// --- computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition ---
inline integer laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, float *h, integer *ldh, float *wr, float * wi, integer *iloz, integer *ihiz, float *z, integer *ldz, float *work, integer *lwork, integer *info)
{
  return slaqr4_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
}
inline integer laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, double *h, integer *ldh, double *wr, double * wi, integer *iloz, integer *ihiz, double *z, integer *ldz, double *work, integer *lwork, integer *info)
{
  return dlaqr4_(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info);
}
inline integer laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, scomplex *h, integer *ldh, scomplex *w, integer *iloz, integer *ihiz, scomplex *z, integer *ldz, scomplex *work, integer *lwork, integer *info)
{
  return claqr4_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
}
inline integer laqr4(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, dcomplex *h, integer *ldh, dcomplex *w, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, dcomplex *work, integer *lwork, integer *info)
{
  return zlaqr4_(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info);
}

// --- performs a single small-bulge multi-shift QR sweep ---
inline integer laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, float *sr, float *si, float *h, integer *ldh, integer *iloz, integer *ihiz, float *z, integer *ldz, float *v, integer *ldv, float *u, integer *ldu, integer *nv, float *wv, integer *ldwv, integer *nh, float *wh, integer * ldwh)
{
  return slaqr5_(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
}
inline integer laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, double *sr, double *si, double *h, integer *ldh, integer *iloz, integer *ihiz, double *z, integer *ldz, double *v, integer *ldv, double *u, integer *ldu, integer *nv, double *wv, integer *ldwv, integer *nh, double *wh, integer * ldwh)
{
  return dlaqr5_(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
}
inline integer laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, scomplex *s, scomplex *h, integer *ldh, integer *iloz, integer *ihiz, scomplex * z, integer *ldz, scomplex *v, integer *ldv, scomplex *u, integer *ldu, integer *nv, scomplex *wv, integer *ldwv, integer *nh, scomplex *wh, integer *ldwh)
{
  return claqr5_(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
}
inline integer laqr5(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, dcomplex *s, dcomplex *h, integer *ldh, integer *iloz, integer *ihiz, dcomplex *z, integer *ldz, dcomplex *v, integer *ldv, dcomplex *u, integer *ldu, integer *nv, dcomplex *wv, integer *ldwv, integer *nh, dcomplex *wh, integer * ldwh)
{
  return zlaqr5_(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh);
}

// --- scales a symmetric/Hermitian band matrix, using scaling factors computed by spbequ. ---
inline integer laqsb(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, float *s, float *scond, float *amax, char *equed)
{
  return slaqsb_(uplo, n, kd, ab, ldab, s, scond, amax, equed);
}
inline integer laqsb(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, double *s, double *scond, double *amax, char *equed)
{
  return dlaqsb_(uplo, n, kd, ab, ldab, s, scond, amax, equed);
}
inline integer laqsb(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *s, float *scond, float *amax, char *equed)
{
  return claqsb_(uplo, n, kd, ab, ldab, s, scond, amax, equed);
}
inline integer laqsb(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *s, double *scond, double *amax, char *equed)
{
  return zlaqsb_(uplo, n, kd, ab, ldab, s, scond, amax, equed);
}

// --- scales a symmetric/Hermitian matrix in packed storage, using scaling factors computed by sppequ ---
inline integer laqsp(char *uplo, integer *n, float *ap, float *s, float * scond, float *amax, char *equed)
{
  return slaqsp_(uplo, n, ap, s, scond, amax, equed);
}
inline integer laqsp(char *uplo, integer *n, double *ap, double *s, double * scond, double *amax, char *equed)
{
  return dlaqsp_(uplo, n, ap, s, scond, amax, equed);
}
inline integer laqsp(char *uplo, integer *n, scomplex *ap, float *s, float * scond, float *amax, char *equed)
{
  return claqsp_(uplo, n, ap, s, scond, amax, equed);
}
inline integer laqsp(char *uplo, integer *n, dcomplex *ap, double *s, double * scond, double *amax, char *equed)
{
  return zlaqsp_(uplo, n, ap, s, scond, amax, equed);
}

// --- scales a symmetric/Hermitian matrix, using scaling factors computed by spoequ ---
inline integer laqsy(char *uplo, integer *n, float *a, integer *lda, float *s, float *scond, float *amax, char *equed)
{
  return slaqsy_(uplo, n, a, lda, s, scond, amax, equed);
}
inline integer laqsy(char *uplo, integer *n, double *a, integer *lda, double *s, double *scond, double *amax, char *equed)
{
  return dlaqsy_(uplo, n, a, lda, s, scond, amax, equed);
}
inline integer laqsy(char *uplo, integer *n, scomplex *a, integer *lda, float *s, float *scond, float *amax, char *equed)
{
  return claqsy_(uplo, n, a, lda, s, scond, amax, equed);
}
inline integer laqsy(char *uplo, integer *n, dcomplex *a, integer *lda, double *s, double *scond, double *amax, char *equed)
{
  return zlaqsy_(uplo, n, a, lda, s, scond, amax, equed);
}

// --- solves a float quasi-triangular system of equations in float arithmetic ---
inline integer laqtr(logical *ltran, logical *lreal, integer *n, float *t, integer *ldt, float *b, float *w, float *scale, float *x, float *work, integer *info)
{
  return slaqtr_(ltran, lreal, n, t, ldt, b, w, scale, x, work, info);
}
inline integer laqtr(logical *ltran, logical *lreal, integer *n, double *t, integer *ldt, double *b, double *w, double *scale, double *x, double *work, integer *info)
{
  return dlaqtr_(ltran, lreal, n, t, ldt, b, w, scale, x, work, info);
}

// --- computes the (scaled) r-th column of the inverse of the submatrix in rows b1 through bn of the tridiagonal matrix LDLT - Î»I ---
inline integer lar1v(integer *n, integer *b1, integer *bn, float *lambda, float *d, float *l, float *ld, float *lld, float *pivmin, float * gaptol, float *z, logical *wantnc, integer *negcnt, float *ztz, float *mingma, integer *r, integer *isuppz, float *nrminv, float *resid, float *rqcorr, float *work)
{
  return slar1v_(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, rqcorr, work);
}
inline integer lar1v(integer *n, integer *b1, integer *bn, double *lambda, double *d, double *l, double *ld, double *lld, double *pivmin, double * gaptol, double *z, logical *wantnc, integer *negcnt, double *ztz, double * mingma, integer *r, integer *isuppz, double *nrminv, double *resid, double *rqcorr, double *work)
{
  return dlar1v_(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, rqcorr, work);
}
inline integer lar1v(integer *n, integer *b1, integer *bn, float *lambda, float *d, float *l, float *ld, float *lld, float *pivmin, float * gaptol, scomplex *z, logical *wantnc, integer *negcnt, float *ztz, float * mingma, integer *r, integer *isuppz, float *nrminv, float *resid, float *rqcorr, float *work)
{
  return clar1v_(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, rqcorr, work);
}
inline integer lar1v(integer *n, integer *b1, integer *bn, double *lambda, double *d, double *l, double *ld, double *lld, double *pivmin, double * gaptol, dcomplex *z, logical *wantnc, integer *negcnt, double *ztz, double * mingma, integer *r, integer *isuppz, double *nrminv, double *resid, double *rqcorr, double *work)
{
  return zlar1v_(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r, isuppz, nrminv, resid, rqcorr, work);
}

// --- applies a vector of plane rotations with float cosines and float sines from both sides to a sequence of 2-by-2 symmetric/Hermitian matrices ---
inline integer lar2v(integer *n, float *x, float *y, float *z, integer *incx, float *c, float *s, integer *incc)
{
  return slar2v_(n, x, y, z, incx, c, s, incc);
}
inline integer lar2v(integer *n, double *x, double *y, double *z, integer *incx, double *c, double *s, integer *incc)
{
  return dlar2v_(n, x, y, z, incx, c, s, incc);
}
inline integer lar2v(integer *n, scomplex *x, scomplex *y, scomplex *z, integer *incx, float *c, scomplex *s, integer *incc)
{
  return clar2v_(n, x, y, z, incx, c, s, incc);
}
inline integer lar2v(integer *n, dcomplex *x, dcomplex *y, dcomplex *z, integer *incx, double *c, dcomplex *s, integer *incc)
{
  return zlar2v_(n, x, y, z, incx, c, s, incc);
}

// --- applies an elementary reflector to a general rectangular matrix ---
inline integer larf(char *side, integer *m, integer *n, float *v, integer *incv, float *tau, float *c, integer *ldc, float *work)
{
  return slarf_(side, m, n, v, incv, tau, c, ldc, work);
}
inline integer larf(char *side, integer *m, integer *n, double *v, integer *incv, double *tau, double *c, integer *ldc, double *work)
{
  return dlarf_(side, m, n, v, incv, tau, c, ldc, work);
}
inline integer larf(char *side, integer *m, integer *n, scomplex *v, integer *incv, scomplex *tau, scomplex *c, integer *ldc, scomplex *work)
{
  return clarf_(side, m, n, v, incv, tau, c, ldc, work);
}
inline integer larf(char *side, integer *m, integer *n, dcomplex *v, integer *incv, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work)
{
  return zlarf_(side, m, n, v, incv, tau, c, ldc, work);
}

// --- applies an elementary reflector, or Householder matrix, H, to an n x n symmetric matrix C, from both the left and the right ---
inline integer larfy(char *uplo, integer *n, float *v, integer *incv, float *tau, float *c, integer *ldc, float *work)
{
  return slarfy_(uplo, n, v, incv, tau, c, ldc, work);
}
inline integer larfy(char *uplo, integer *n, double *v, integer *incv, double *tau, double *c, integer *ldc, double *work)
{
  return dlarfy_(uplo, n, v, incv, tau, c, ldc, work);
}
inline integer larfy(char *uplo, integer *n, scomplex *v, integer *incv, scomplex *tau, scomplex *c, integer *ldc, scomplex *work)
{
  return clarfy_(uplo, n, v, incv, tau, c, ldc, work);
}
inline integer larfy(char *uplo, integer *n, dcomplex *v, integer *incv, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work)
{
  return zlarfy_(uplo, n, v, incv, tau, c, ldc, work);
}

// --- generates a vector of plane rotations with float cosines and float sines ---
inline integer largv(integer *n, float *x, integer *incx, float *y, integer *incy, float *c, integer *incc)
{
  return slargv_(n, x, incx, y, incy, c, incc);
}
inline integer largv(integer *n, double *x, integer *incx, double *y, integer *incy, double *c, integer *incc)
{
  return dlargv_(n, x, incx, y, incy, c, incc);
}
inline integer largv(integer *n, scomplex *x, integer *incx, scomplex *y, integer *incy, float *c, integer *incc)
{
  return clargv_(n, x, incx, y, incy, c, incc);
}
inline integer largv(integer *n, dcomplex *x, integer *incx, dcomplex *y, integer *incy, double *c, integer *incc)
{
  return zlargv_(n, x, incx, y, incy, c, incc);
}

//--- computes the splitting points with the specified threshold ---
inline integer larra(integer *n, float *d, float *e, float *e2, float * spltol, float *tnrm, integer *nsplit, integer *isplit, integer *info)
{
  return slarra_(n, d, e, e2, spltol, tnrm, nsplit, isplit, info);
}
inline integer larra(integer *n, double *d, double *e, double *e2, double * spltol, double *tnrm, integer *nsplit, integer *isplit, integer *info)
{
  return dlarra_(n, d, e, e2, spltol, tnrm, nsplit, isplit, info);
}

//--- provides limited bisection to locate eigenvalues for more accuracy ---
inline integer larrb(integer *n, float *d, float *lld, integer * ifirst, integer *ilast, float *rtol1, float *rtol2, integer *offset, float *w, float *wgap, float *werr, float *work, integer *iwork, float * pivmin, float *spdiam, integer *twist, integer *info)
{
  return slarrb_(n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, work, iwork, pivmin, spdiam, twist, info);
}
inline integer larrb(integer *n, double *d, double *lld, integer * ifirst, integer *ilast, double *rtol1, double *rtol2, integer *offset, double *w, double *wgap, double *werr, double *work, integer *iwork, double * pivmin, double *spdiam, integer *twist, integer *info)
{
  return dlarrb_(n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, work, iwork, pivmin, spdiam, twist, info);
}

//--- computes the number of eigenvalues of the symmetric tridiagonal matrix ---
inline integer larrc(char *jobt, integer *n, float *vl, float *vu, float *d, float *e, float *pivmin, integer *eigcnt, integer *lcnt, integer * rcnt, integer *info)
{
  return slarrc_(jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info); 
}
inline integer larrc(char *jobt, integer *n, double *vl, double *vu, double *d, double *e, double *pivmin, integer *eigcnt, integer *lcnt, integer * rcnt, integer *info)
{
  return dlarrc_(jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info); 
}

//--- computes the eigenvalues of a symmetric tridiagonal matrix to suitable accuracy ---
inline integer larrd(char *range, char *order, integer *n, float *vl, float *vu, integer *il, integer *iu, float *gers, float *reltol, float * d, float *e, float *e2, float *pivmin, integer *nsplit, integer * isplit, integer *m, float *w, float *werr, float *wl, float *wu, integer * iblock, integer *indexw, float *work, integer *iwork, integer *info)
{
  return slarrd_(range, order, n, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w, werr, wl, wu, iblock, indexw, work, iwork, info);
}
inline integer larrd(char *range, char *order, integer *n, double *vl, double *vu, integer *il, integer *iu, double *gers, double *reltol, double * d, double *e, double *e2, double *pivmin, integer *nsplit, integer * isplit, integer *m, double *w, double *werr, double *wl, double *wu, integer * iblock, integer *indexw, double *work, integer *iwork, integer *info)
{
  return dlarrd_(range, order, n, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w, werr, wl, wu, iblock, indexw, work, iwork, info);
}

//--- given the tridiagonal matrix T, sets small off-diagonal elements to zero and for each unreduced block Ti, finds base representations and eigenvalues ---
inline integer larre(char *range, integer *n, float *vl, float *vu, integer *il, integer *iu, float *d, float *e, float *e2, float *rtol1, float *rtol2, float *spltol, integer *nsplit, integer *isplit, integer * m, float *w, float *werr, float *wgap, integer *iblock, integer *indexw, float *gers, float *pivmin, float *work, integer *iwork, integer *info)
{
  return slarre_(range, n, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr, wgap, iblock, indexw, gers, pivmin, work, iwork, info);
}
inline integer larre(char *range, integer *n, double *vl, double *vu, integer *il, integer *iu, double *d, double *e, double *e2, double *rtol1, double *rtol2, double *spltol, integer *nsplit, integer *isplit, integer * m, double *w, double *werr, double *wgap, integer *iblock, integer *indexw, double *gers, double *pivmin, double *work, integer *iwork, integer *info)
{
  return dlarre_(range, n, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr, wgap, iblock, indexw, gers, pivmin, work, iwork, info);
}

//--- finds a new relatively robust representation such that at least one of the eigenvalues is relatively isolated ---
inline integer larrf(integer *n, float *d, float *l, float *ld, integer *clstrt, integer *clend, float *w, float *wgap, float *werr, float *spdiam, float *clgapl, float *clgapr, float *pivmin, float *sigma, float *dplus, float *lplus, float *work, integer *info)
{
  return slarrf_(n, d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin, sigma, dplus, lplus, work, info);
}
inline integer larrf(integer *n, double *d, double *l, double *ld, integer *clstrt, integer *clend, double *w, double *wgap, double *werr, double *spdiam, double *clgapl, double *clgapr, double *pivmin, double *sigma, double *dplus, double *lplus, double *work, integer *info)
{
  return dlarrf_(n, d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin, sigma, dplus, lplus, work, info);
}

//--- performs refinement of the initial estimates of the eigenvalues of the matrix T ---
inline integer larrj(integer *n, float *d, float *e2, integer *ifirst, integer *ilast, float *rtol, integer *offset, float *w, float *werr, float *work, integer *iwork, float *pivmin, float *spdiam, integer *info)
{
  return slarrj_(n, d, e2, ifirst, ilast, rtol, offset, w, werr, work, iwork, pivmin, spdiam, info);
}
inline integer larrj(integer *n, double *d, double *e2, integer *ifirst, integer *ilast, double *rtol, integer *offset, double *w, double *werr, double *work, integer *iwork, double *pivmin, double *spdiam, integer *info)
{
  return dlarrj_(n, d, e2, ifirst, ilast, rtol, offset, w, werr, work, iwork, pivmin, spdiam, info);
}

//--- computes one eigenvalue of a symmetric tridiagonal matrix T to suitable accuracy ---
inline integer larrk(integer *n, integer *iw, float *gl, float *gu, float *d, float *e2, float *pivmin, float *reltol, float *w, float *werr, integer *info)
{
  return slarrk_(n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info);
}
inline integer larrk(integer *n, integer *iw, double *gl, double *gu, double *d, double *e2, double *pivmin, double *reltol, double *w, double *werr, integer *info)
{
  return dlarrk_(n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info);
}

//--- performs tests to decide whether the symmetric tridiagonal matrix T warrants expensive computations ---
inline integer larrr(integer *n, float *d, float *e, integer *info)
{
  return slarrr_(n, d, e, info); 
}
inline integer larrr(integer *n, double *d, double *e, integer *info)
{
  return dlarrr_(n, d, e, info); 
}

// --- computes the eigenvectors of the tridiagonal matrix T = L D LT given L, D and the eigenvalues of L D LT ---
inline integer larrv(integer *n, float *vl, float *vu, float *d, float *l, float *pivmin, integer *isplit, integer *m, integer *dol, integer *dou, float *minrgp, float *rtol1, float *rtol2, float *w, float *werr, float *wgap, integer *iblock, integer *indexw, float *gers, float *z, integer *ldz, integer *isuppz, float *work, integer *iwork, integer * info)
{
  return slarrv_(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);
}
inline integer larrv(integer *n, double *vl, double *vu, double *d, double *l, double *pivmin, integer *isplit, integer *m, integer *dol, integer *dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, double *wgap, integer *iblock, integer *indexw, double *gers, double *z, integer *ldz, integer *isuppz, double *work, integer *iwork, integer * info)
{
  return dlarrv_(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);
}
inline integer larrv(integer *n, float *vl, float *vu, float *d, float * l, float *pivmin, integer *isplit, integer *m, integer *dol, integer * dou, float *minrgp, float *rtol1, float *rtol2, float *w, float *werr, float *wgap, integer *iblock, integer *indexw, float *gers, scomplex * z, integer *ldz, integer *isuppz, float *work, integer *iwork, integer *info)
{
  return clarrv_(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);;
}
inline integer larrv(integer *n, double *vl, double *vu, double *d, double *l, double *pivmin, integer *isplit, integer *m, integer *dol, integer *dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, double *wgap, integer *iblock, integer *indexw, double *gers, dcomplex *z, integer *ldz, integer *isuppz, double *work, integer *iwork, integer * info)
{
  return zlarrv_(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info);
}

// --- performs reciprocal diagonal scaling on a vector ---
inline integer larscl2(integer *m, integer *n, float *d, float *x, integer *ldx)
{
  return slarscl2_(m, n, d, x, ldx);
}
inline integer larscl2(integer *m, integer *n, double *d, double *x, integer *ldx)
{
  return dlarscl2_(m, n, d, x, ldx);
}
inline integer larscl2(integer *m, integer *n, float *d, scomplex *x, integer *ldx)
{
  return clarscl2_(m, n, d, x, ldx);
}
inline integer larscl2(integer *m, integer *n, double *d, dcomplex *x, integer *ldx)
{
  return zlarscl2_(m, n, d, x, ldx);
}

// --- generates a plane rotation with float cosine and float sine ---
inline integer lartg(float *f, float *g, float *cs, float *sn, float *r__)
{
  return slartg_( f, g, cs, sn, r__);
}
inline integer lartg(double *f, double *g, double *cs, double *sn, double *r__)
{
  return dlartg_( f, g, cs, sn, r__);
}
inline integer lartg(scomplex *f, scomplex *g, float *cs, scomplex *sn, scomplex *r__)
{
  return clartg_( f, g, cs, sn, r__);
}
inline integer lartg(dcomplex *f, dcomplex *g, double *cs, dcomplex *sn, dcomplex *r__)
{
  return zlartg_( f, g, cs, sn, r__);
}

// --- applies a vector of plane rotations with float cosines and float sines to the elements of a pair of vectors ---
inline integer lartv(integer *n, float *x, integer *incx, float *y, integer *incy, float *c, float *s, integer *incc)
{
  return slartv_(n, x, incx, y, incy, c, s, incc);
}
inline integer lartv(integer *n, double *x, integer *incx, double *y, integer *incy, double *c, double *s, integer *incc)
{
  return dlartv_(n, x, incx, y, incy, c, s, incc);
}
inline integer lartv(integer *n, scomplex *x, integer *incx, scomplex *y, integer *incy, float *c, scomplex *s, integer *incc)
{
  return clartv_(n, x, incx, y, incy, c, s, incc);
}
inline integer lartv(integer *n, dcomplex *x, integer *incx, dcomplex *y, integer *incy, double *c, dcomplex *s, integer *incc)
{
  return zlartv_(n, x, incx, y, incy, c, s, incc);
}

// --- returns a vector of n random float numbers from a uniform distribution ---
inline integer laruv(integer *iseed, integer *n, float *x)
{
  return slaruv_(iseed, n, x);  
}
inline integer laruv(integer *iseed, integer *n, double *x)
{
  return dlaruv_(iseed, n, x);  
}

// --- applies an elementary reflector (as returned by stzrzf) to a general matrix ---
inline integer larz(char *side, integer *m, integer *n, integer *l, float *v, integer *incv, float *tau, float *c, integer *ldc, float *work)
{
  return slarz_(side, m, n, l, v, incv, tau, c, ldc, work);
}
inline integer larz(char *side, integer *m, integer *n, integer *l, double *v, integer *incv, double *tau, double *c, integer *ldc, double *work)
{
  return dlarz_(side, m, n, l, v, incv, tau, c, ldc, work);
}
inline integer larz(char *side, integer *m, integer *n, integer *l, scomplex *v, integer *incv, scomplex *tau, scomplex *c, integer *ldc, scomplex *work)
{
  return clarz_(side, m, n, l, v, incv, tau, c, ldc, work);
}
inline integer larz(char *side, integer *m, integer *n, integer *l, dcomplex *v, integer *incv, dcomplex *tau, dcomplex *c, integer *ldc, dcomplex *work)
{
  return zlarz_(side, m, n, l, v, incv, tau, c, ldc, work);
}

// --- applies a block reflector or its transpose to a general matrix ---
inline integer larzb(char *side, char *trans, char *direct, char * storev, integer *m, integer *n, integer *k, integer *l, float *v, integer *ldv, float *t, integer *ldt, float *c, integer *ldc, float *work, integer *ldwork)
{
  return slarzb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline integer larzb(char *side, char *trans, char *direct, char * storev, integer *m, integer *n, integer *k, integer *l, double *v, integer *ldv, double *t, integer *ldt, double *c, integer *ldc, double *work, integer *ldwork)
{
  return dlarzb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline integer larzb(char *side, char *trans, char *direct, char * storev, integer *m, integer *n, integer *k, integer *l, scomplex *v, integer *ldv, scomplex *t, integer *ldt, scomplex *c, integer *ldc, scomplex *work, integer *ldwork)
{
  return clarzb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work, ldwork);
}
inline integer larzb(char *side, char *trans, char *direct, char * storev, integer *m, integer *n, integer *k, integer *l, dcomplex *v, integer *ldv, dcomplex *t, integer *ldt, dcomplex *c, integer *ldc, dcomplex *work, integer *ldwork)
{
  return zlarzb_(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work, ldwork);
}

// --- forms the triangular factor T of a block reflector H = I - vtvH. ---
inline integer larzt(char *direct, char *storev, integer *n, integer * k, float *v, integer *ldv, float *tau, float *t, integer *ldt)
{
  return slarzt_(direct, storev, n, k, v, ldv, tau, t, ldt);
}
inline integer larzt(char *direct, char *storev, integer *n, integer * k, double *v, integer *ldv, double *tau, double *t, integer *ldt)
{
  return dlarzt_(direct, storev, n, k, v, ldv, tau, t, ldt);
}
inline integer larzt(char *direct, char *storev, integer *n, integer * k, scomplex *v, integer *ldv, scomplex *tau, scomplex *t, integer *ldt)
{
  return clarzt_(direct, storev, n, k, v, ldv, tau, t, ldt);
}
inline integer larzt(char *direct, char *storev, integer *n, integer * k, dcomplex *v, integer *ldv, dcomplex *tau, dcomplex *t, integer *ldt)
{
  return zlarzt_(direct, storev, n, k, v, ldv, tau, t, ldt);
}

// --- computes singular values of a 2-by-2 triangular matrix ---
inline integer las2(float *f, float *g, float *h, float *ssmin, float * ssmax)
{
  return slas2_(f, g, h, ssmin, ssmax); 
}
inline integer las2(double *f, double *g, double *h, double *ssmin, double * ssmax)
{
  return dlas2_(f, g, h, ssmin, ssmax); 
}

// --- performs diagonal scaling on a vector ---
inline integer lascl2(integer *m, integer *n, float *d, float *x, integer *ldx)
{
  return slascl2_(m, n, d, x, ldx);
}
inline integer lascl2(integer *m, integer *n, double *d, double *x, integer *ldx)
{
  return dlascl2_(m, n, d, x, ldx);
}
inline integer lascl2(integer *m, integer *n, float *d, scomplex *x, integer *ldx)
{
  return clascl2_(m, n, d, x, ldx);
}
inline integer lascl2(integer *m, integer *n, double *d, dcomplex *x, integer *ldx)
{
  return zlascl2_(m, n, d, x, ldx);
}

// --- computes the singular values of a float upper bidiagonal n-by-m matrix B with diagonal d and off-diagonal e ---
inline integer lasd0(integer *n, integer *sqre, float *d, float *e, float *u, integer *ldu, float *vt, integer *ldvt, integer *smlsiz, integer *iwork, float *work, integer *info)
{
  return slasd0_(n, sqre, d, e, u, ldu, vt, ldvt, smlsiz, iwork, work, info);
}
inline integer lasd0(integer *n, integer *sqre, double *d, double *e, double *u, integer *ldu, double *vt, integer *ldvt, integer *smlsiz, integer *iwork, double *work, integer *info)
{
  return dlasd0_(n, sqre, d, e, u, ldu, vt, ldvt, smlsiz, iwork, work, info); 
}

// --- computes the SVD of an upper bidiagonal matrix B of the specified size ---
inline integer lasd1(integer *nl, integer *nr, integer *sqre, float * d, float *alpha, float *beta, float *u, integer *ldu, float *vt, integer *ldvt, integer *idxq, integer *iwork, float *work, integer * info)
{
  return slasd1_(nl, nr, sqre, d, alpha, beta, u, ldu, vt, ldvt, idxq, iwork, work, info);
}
inline integer lasd1(integer *nl, integer *nr, integer *sqre, double * d, double *alpha, double *beta, double *u, integer *ldu, double *vt, integer *ldvt, integer *idxq, integer *iwork, double *work, integer * info)
{
  return dlasd1_(nl, nr, sqre, d, alpha, beta, u, ldu, vt, ldvt, idxq, iwork, work, info); 
}

// --- merges the two sets of singular values together into a single sorted set ---
inline integer lasd2(integer *nl, integer *nr, integer *sqre, integer *k, float *d, float *z, float *alpha, float *beta, float *u, integer * ldu, float *vt, integer *ldvt, float *dsigma, float *u2, integer *ldu2, float *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info)
{
  return slasd2_(nl, nr, sqre, k, d, z, alpha, beta, u, ldu, vt, ldvt, dsigma, u2, ldu2, vt2, ldvt2, idxp, idx, idxc, idxq, coltyp, info); 
}
inline integer lasd2(integer *nl, integer *nr, integer *sqre, integer *k, double *d, double *z, double *alpha, double *beta, double *u, integer * ldu, double *vt, integer *ldvt, double *dsigma, double *u2, integer *ldu2, double *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info)
{
  return dlasd2_(nl, nr, sqre, k, d, z, alpha, beta, u, ldu, vt, ldvt, dsigma, u2, ldu2, vt2, ldvt2, idxp, idx, idxc, idxq, coltyp, info); 
}

// --- finds all square roots of the roots of the secular equation ---
inline integer lasd3(integer *nl, integer *nr, integer *sqre, integer *k, float *d, float *q, integer *ldq, float *dsigma, float *u, integer * ldu, float *u2, integer *ldu2, float *vt, integer *ldvt, float *vt2, integer *ldvt2, integer *idxc, integer *ctot, float *z, integer * info)
{
  return slasd3_(nl, nr, sqre, k, d, q, ldq, dsigma, u, ldu, u2, ldu2, vt, ldvt, vt2, ldvt2, idxc, ctot, z, info);
}
inline integer lasd3(integer *nl, integer *nr, integer *sqre, integer *k, double *d, double *q, integer *ldq, double *dsigma, double *u, integer * ldu, double *u2, integer *ldu2, double *vt, integer *ldvt, double *vt2, integer *ldvt2, integer *idxc, integer *ctot, double *z, integer * info)
{
  return dlasd3_(nl, nr, sqre, k, d, q, ldq, dsigma, u, ldu, u2, ldu2, vt, ldvt, vt2, ldvt2, idxc, ctot, z, info); 
}

// --- computes the square root of the i-th updated eigenvalue of a positive symmetric rank-one modification to a positive diagonal matrix ---
inline integer lasd4(integer *n, integer *i, float *d, float *z, float *delta, float *rho, float *sigma, float *work, integer *info)
{
  return slasd4_(n, i, d, z, delta, rho, sigma, work, info);
}
inline integer lasd4(integer *n, integer *i, double *d, double *z, double *delta, double *rho, double *sigma, double *work, integer *info)
{
  return dlasd4_(n, i, d, z, delta, rho, sigma, work, info); 
}

// --- computes the square root of the i-th eigenvalue of a positive symmetric rank-one modification of a 2-by-2 diagonal matrix ---
inline integer lasd5(integer *i, float *d, float *z, float *delta, float *rho, float *dsigma, float *work)
{
  return slasd5_(i, d, z, delta, rho, dsigma, work); 
}
inline integer lasd5(integer *i, double *d, double *z, double *delta, double *rho, double *dsigma, double *work)
{
  return dlasd5_(i, d, z, delta, rho, dsigma, work);  
}

// --- computes the SVD of an updated upper bidiagonal matrix obtained by merging two smaller ones by appending a row ---
inline integer lasd6(integer *icompq, integer *nl, integer *nr, integer *sqre, float *d, float *vf, float *vl, float *alpha, float *beta, integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, float *givnum, integer *ldgnum, float *poles, float * difl, float *difr, float *z, integer *k, float *c, float *s, float * work, integer *iwork, integer *info)
{
  return slasd6_(icompq, nl, nr, sqre, d, vf, vl, alpha, beta, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, iwork, info); 
}
inline integer lasd6(integer *icompq, integer *nl, integer *nr, integer *sqre, double *d, double *vf, double *vl, double *alpha, double *beta, integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, double *givnum, integer *ldgnum, double *poles, double * difl, double *difr, double *z, integer *k, double *c, double *s, double * work, integer *iwork, integer *info)
{
  return dlasd6_(icompq, nl, nr, sqre, d, vf, vl, alpha, beta, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, iwork, info); 
}

// --- merges the two sets of singular values together into a single sorted set. Then it tries to deflate the size of the problem ---
inline integer lasd7(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, float *d, float *z, float *zw, float *vf, float *vfw, float *vl, float *vlw, float *alpha, float *beta, float *dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, integer * givptr, integer *givcol, integer *ldgcol, float *givnum, integer * ldgnum, float *c, float *s, integer *info)
{
  return slasd7_(icompq, nl, nr, sqre, k, d, z, zw, vf, vfw, vl, vlw, alpha, beta, dsigma, idx, idxp, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info); 
}
inline integer lasd7(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, double *d, double *z, double *zw, double *vf, double *vfw, double *vl, double *vlw, double *alpha, double *beta, double *dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, integer * givptr, integer *givcol, integer *ldgcol, double *givnum, integer * ldgnum, double *c, double *s, integer *info)
{
  return dlasd7_(icompq, nl, nr, sqre, k, d, z, zw, vf, vfw, vl, vlw, alpha, beta, dsigma, idx, idxp, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info); 
}

// --- finds the square roots of the roots of the secular equation, and stores, for each element in D, the distance to its two nearest poles ---
inline integer lasd8(integer *icompq, integer *k, float *d, float * z, float *vf, float *vl, float *difl, float *difr, integer *lddifr, float *dsigma, float *work, integer *info)
{
  return slasd8_(icompq, k, d, z, vf, vl, difl, difr, lddifr, dsigma, work, info); 
}
inline integer lasd8(integer *icompq, integer *k, double *d, double * z, double *vf, double *vl, double *difl, double *difr, integer *lddifr, double *dsigma, double *work, integer *info)
{
  return dlasd8_(icompq, k, d, z, vf, vl, difl, difr, lddifr, dsigma, work, info); 
}

// --- computes the singular value decomposition (SVD) of a float upper bidiagonal matrix with diagonal d and off-diagonal e ---
inline integer lasda(integer *icompq, integer *smlsiz, integer *n, integer *sqre, float *d, float *e, float *u, integer *ldu, float *vt, integer *k, float *difl, float *difr, float *z, float *poles, integer * givptr, integer *givcol, integer *ldgcol, integer *perm, float *givnum, float *c, float *s, float *work, integer *iwork, integer *info)
{
  return slasda_(icompq, smlsiz, n, sqre, d, e, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info); 
}
inline integer lasda(integer *icompq, integer *smlsiz, integer *n, integer *sqre, double *d, double *e, double *u, integer *ldu, double *vt, integer *k, double *difl, double *difr, double *z, double *poles, integer * givptr, integer *givcol, integer *ldgcol, integer *perm, double *givnum, double *c, double *s, double *work, integer *iwork, integer *info)
{
  return dlasda_(icompq, smlsiz, n, sqre, d, e, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info); 
}

// --- computes the SVD of a float bidiagonal matrix with diagonal d and off-diagonal e ---
inline integer lasdq(char *uplo, integer *sqre, integer *n, integer * ncvt, integer *nru, integer *ncc, float *d, float *e, float *vt, integer *ldvt, float *u, integer *ldu, float *c, integer *ldc, float * work, integer *info)
{
  return slasdq_(uplo, sqre, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info); 
}
inline integer lasdq(char *uplo, integer *sqre, integer *n, integer * ncvt, integer *nru, integer *ncc, double *d, double *e, double *vt, integer *ldvt, double *u, integer *ldu, double *c, integer *ldc, double * work, integer *info)
{
  return dlasdq_(uplo, sqre, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info); 
}

// --- computes the singular values of a float square bidiagonal matrix ---
inline integer lasq1(integer *n, float *d, float *e, float *work, integer *info)
{
  return slasq1_(n, d, e, work, info);
}
inline integer lasq1(integer *n, double *d, double *e, double *work, integer *info)
{
  return dlasq1_(n, d, e, work, info);
}

// --- computes all the eigenvalues of the symmetric positive definite tridiagonal matrix associated with the qd Array Z to high relative accuracy. ---
inline integer lasq2(integer *n, float *z, integer *info)
{
  return slasq2_(n, z, info);
}
inline integer lasq2(integer *n, double *z, integer *info)
{
  return dlasq2_(n, z, info);
}

// --- checks for deflation, computes a shift and calls dqds ---
inline integer lasq3(integer *i0, integer *n0, float *z, integer *pp, float *dmin, float *sigma, float *desig, float *qmax, integer *nfail, integer *iter, integer *ndiv, logical *ieee, integer *ttype, float * dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *g, float * tau)
{
  return slasq3_(i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau);
}
inline integer lasq3(integer *i0, integer *n0, double *z, integer *pp, double *dmin, double *sigma, double *desig, double *qmax, integer *nfail, integer *iter, integer *ndiv, logical *ieee, integer *ttype, double * dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *g, double * tau)
{
  return dlasq3_(i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau);
}

// --- computes an approximation to the smallest eigenvalue using values of d from the previous transform  ---
inline integer lasq4(integer *i0, integer *n0, float *z, integer *pp, integer *n0in, float *dmin, float *dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *tau, integer *ttype, float *g)
{
  return slasq4_(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g);
}
inline integer lasq4(integer *i0, integer *n0, double *z, integer *pp, integer *n0in, double *dmin, double *dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *tau, integer *ttype, double *g)
{
  return dlasq4_(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g);
}

// --- computes one dqds transform in ping-pong form ---
inline integer lasq5(integer *i0, integer *n0, float *z, integer *pp, float *tau, float *sigma, float *dmin, float *dmin1, float *dmin2, float *dn, float *dnm1, float *dnm2, logical *ieee, float *eps)
{
  return slasq5_(i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2, ieee, eps);
}
inline integer lasq5(integer *i0, integer *n0, double *z, integer *pp, double *tau, double *sigma, double *dmin, double *dmin1, double *dmin2, double *dn, double *dnm1, double *dnm2, logical *ieee, double *eps)
{
  return dlasq5_(i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2, ieee, eps);
}

// --- computes one dqd transform in ping-pong form ---
inline integer lasq6(integer *i0, integer *n0, float *z, integer *pp, float *dmin, float *dmin1, float *dmin2, float *dn, float *dnm1, float * dnm2) 
{
  return slasq6_(i0, n0, z, pp, dmin, dmin1, dmin2, dn, dnm1, dnm2);
}
inline integer lasq6(integer *i0, integer *n0, double *z, integer *pp, double *dmin, double *dmin1, double *dmin2, double *dn, double *dnm1, double * dnm2)
{
  return dlasq6_(i0, n0, z, pp, dmin, dmin1, dmin2, dn, dnm1, dnm2);
}

// --- applies a sequence of plane rotations to a general rectangular matrix ---
inline integer lasr(char *side, char *pivot, char *direct, integer *m, integer *n, float *c, float *s, float *a, integer *lda)
{
  return slasr_(side, pivot, direct, m, n, c, s, a, lda);
}
inline integer lasr_(char *side, char *pivot, char *direct, integer *m, integer *n, double *c, double *s, double *a, integer *lda)
{
  return dlasr_(side, pivot, direct, m, n, c, s, a, lda);
}
inline integer lasr(char *side, char *pivot, char *direct, integer *m, integer *n, float *c, float *s, scomplex *a, integer *lda)
{
  return clasr_(side, pivot, direct, m, n, c, s, a, lda);
}
inline integer lasr(char *side, char *pivot, char *direct, integer *m, integer *n, double *c, double *s, dcomplex *a, integer *lda)
{
  return zlasr_(side, pivot, direct, m, n, c, s, a, lda);
}

// --- computes the singular value decomposition of a 2-by-2 triangular matrix ---
inline integer lasv2(float *f, float *g, float *h, float *ssmin, float * ssmax, float *snr, float *csr, float *snl, float *csl)
{
  return slasv2_( f, g, h, ssmin, ssmax, snr, csr, snl, csl);
}
inline integer lasv2(double *f, double *g, double *h, double *ssmin, double * ssmax, double *snr, double *csr, double *snl, double *csl)
{
  return dlasv2_( f, g, h, ssmin, ssmax, snr, csr, snl, csl);
}

// --- computes a blocked Tall-Skinny LQ factorization of a float M-by-N matrix A for M <= N ---
inline integer laswlq(integer *m, integer *n, integer *mb, integer * nb, float *a, integer *lda, float *t, integer *ldt, float *work, integer *lwork, integer *info)
{
  return slaswlq_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline integer laswlq(integer *m, integer *n, integer *mb, integer * nb, double *a, integer *lda, double *t, integer *ldt, double *work, integer *lwork, integer *info)
{
  return dlaswlq_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline integer laswlq(integer *m, integer *n, integer *mb, integer * nb, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *work, integer *lwork, integer *info)
{
  return claswlq_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline integer laswlq(integer *m, integer *n, integer *mb, integer * nb, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *work, integer *lwork, integer *info)
{
  return zlaswlq_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}

// --- solves the Sylvester matrix equation where the matrices are of order 1 or 2 ---
inline integer lasy2(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2, float *tl, integer *ldtl, float *tr, integer * ldtr, float *b, integer *ldb, float *scale, float *x, integer *ldx, float *xnorm, integer *info)
{
  return slasy2_(ltranl, ltranr, isgn, n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx, xnorm, info);
}
inline integer lasy2(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2, double *tl, integer *ldtl, double *tr, integer * ldtr, double *b, integer *ldb, double *scale, double *x, integer *ldx, double *xnorm, integer *info)
{
  return dlasy2_(ltranl, ltranr, isgn, n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx, xnorm, info);
}

// --- computes a partial factorization of a float symmetric matrix using the Bunch-Kaufman diagonal pivoting method ---
inline integer lasyf(char *uplo, integer *n, integer *nb, integer *kb, float *a, integer *lda, integer *ipiv, float *w, integer *ldw, integer *info)
{
  return slasyf_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline integer lasyf(char *uplo, integer *n, integer *nb, integer *kb, double *a, integer *lda, integer *ipiv, double *w, integer *ldw, integer *info)
{
  return dlasyf_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline integer lasyf(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, integer *ipiv, scomplex *w, integer *ldw, integer *info)
{
  return clasyf_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline integer lasyf(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *w, integer *ldw, integer *info)
{
  return zlasyf_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}

// --- factorizes a panel of a float symmetric matrix A using the Aasen's algorithm ---
inline integer lasyf_aa(char *uplo, integer *j1, integer *m, integer *nb, float *a, integer *lda, integer *ipiv, float *h, integer *ldh, float *work)
{
  return slasyf_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
}
inline integer lasyf_aa(char *uplo, integer *j1, integer *m, integer *nb, double *a, integer *lda, integer *ipiv, double *h, integer *ldh, double *work)
{
  return dlasyf_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
}
inline integer lasyf_aa(char *uplo, integer *j1, integer *m, integer *nb, scomplex *a, integer *lda, integer *ipiv, scomplex *h, integer *ldh, scomplex *work)
{
  return clasyf_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
}
inline integer lasyf_aa(char *uplo, integer *j1, integer *m, integer *nb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *h, integer *ldh, dcomplex *work)
{
  return zlasyf_aa_(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
}

// --- computes a partial factorization of a float symmetric indefinite matrix using bounded Bunch-Kaufman (rook) diagonal pivoting method ---
inline integer lasyf_rk(char *uplo, integer *n, integer *nb, integer *kb, float *a, integer *lda, float *e, integer *ipiv, float *w, integer * ldw, integer *info)
{
  return slasyf_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}
inline integer lasyf_rk(char *uplo, integer *n, integer *nb, integer *kb, double *a, integer *lda, double *e, integer *ipiv, double *w, integer * ldw, integer *info)
{
  return dlasyf_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}
inline integer lasyf_rk(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, scomplex *e, integer *ipiv, scomplex *w, integer * ldw, integer *info)
{
  return clasyf_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}
inline integer lasyf_rk(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, dcomplex *e, integer *ipiv, dcomplex *w, integer * ldw, integer *info)
{
  return zlasyf_rk_(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
}

// --- computes a partial factorization of a float symmetric matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method ---
inline integer lasyf_rook(char *uplo, integer *n, integer *nb, integer *kb, float *a, integer *lda, integer *ipiv, float *w, integer * ldw, integer *info)
{
  return slasyf_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline integer lasyf_rook(char *uplo, integer *n, integer *nb, integer *kb, double *a, integer *lda, integer *ipiv, double *w, integer * ldw, integer *info)
{
  return dlasyf_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline integer lasyf_rook(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda, integer *ipiv, scomplex *w, integer * ldw, integer *info)
{
  return clasyf_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}
inline integer lasyf_rook(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda, integer *ipiv, dcomplex *w, integer * ldw, integer *info)
{
  return zlasyf_rook_(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
}

// --- solves a triangular banded system of equations ---
inline integer latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, float *ab, integer *ldab, float *x, float *scale, float *cnorm, integer *info)
{
  return slatbs_(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
}
inline integer latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, double *ab, integer *ldab, double *x, double *scale, double *cnorm, integer *info)
{
  return dlatbs_(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
}
inline integer latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, scomplex *ab, integer *ldab, scomplex *x, float *scale, float *cnorm, integer *info)
{
  return clatbs_(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
}
inline integer latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, dcomplex *ab, integer *ldab, dcomplex *x, double *scale, double *cnorm, integer *info)
{
  return zlatbs_(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
}

// --- uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contribution to the reciprocal Dif-estimate---
inline integer latdf(integer *ijob, integer *n, float *z, integer * ldz, float *rhs, float *rdsum, float *rdscal, integer *ipiv, integer * jpiv)
{
  return slatdf_(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
}
inline integer latdf(integer *ijob, integer *n, double *z, integer * ldz, double *rhs, double *rdsum, double *rdscal, integer *ipiv, integer * jpiv)
{
  return dlatdf_(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
}
inline integer latdf(integer *ijob, integer *n, scomplex *z, integer * ldz, scomplex *rhs, float *rdsum, float *rdscal, integer *ipiv, integer * jpiv)
{
  return clatdf_(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
}
inline integer latdf(integer *ijob, integer *n, dcomplex *z, integer * ldz, dcomplex *rhs, double *rdsum, double *rdscal, integer *ipiv, integer * jpiv)
{
  return zlatdf_(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
}

// --- solves a triangular system of equations with the matrix held in packed storage ---
inline integer latps(char *uplo, char *trans, char *diag, char * normin, integer *n, float *ap, float *x, float *scale, float *cnorm, integer *info)
{
  return slatps_(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
}
inline integer latps(char *uplo, char *trans, char *diag, char * normin, integer *n, double *ap, double *x, double *scale, double *cnorm, integer *info)
{
  return dlatps_(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
}
inline integer latps(char *uplo, char *trans, char *diag, char * normin, integer *n, scomplex *ap, scomplex *x, float *scale, float *cnorm, integer *info)
{
  return clatps_(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
}
inline integer latps(char *uplo, char *trans, char *diag, char * normin, integer *n, dcomplex *ap, dcomplex *x, double *scale, double *cnorm, integer *info)
{
  return zlatps_(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
}

// --- reduces the first nb rows and columns of a symmetric/Hermitian matrix A to float tridiagonal form by an orthogonal similarity transformation ---
inline integer latrd(char *uplo, integer *n, integer *nb, float *a, integer *lda, float *e, float *tau, float *w, integer *ldw)
{
  return slatrd_(uplo, n, nb, a, lda, e, tau, w, ldw);
}
inline integer latrd(char *uplo, integer *n, integer *nb, double *a, integer *lda, double *e, double *tau, double *w, integer *ldw)
{
  return dlatrd_(uplo, n, nb, a, lda, e, tau, w, ldw);
}
inline integer latrd(char *uplo, integer *n, integer *nb, scomplex *a, integer *lda, float *e, scomplex *tau, scomplex *w, integer *ldw)
{
  return clatrd_(uplo, n, nb, a, lda, e, tau, w, ldw);
}
inline integer latrd(char *uplo, integer *n, integer *nb, dcomplex *a, integer *lda, double *e, dcomplex *tau, dcomplex *w, integer *ldw)
{
  return zlatrd_(uplo, n, nb, a, lda, e, tau, w, ldw);
}

// --- solves a triangular system of equations with the scale factor set to prevent overflow ---
inline integer latrs(char *uplo, char *trans, char *diag, char * normin, integer *n, float *a, integer *lda, float *x, float *scale, float *cnorm, integer *info)
{
  return slatrs_(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
}
inline integer latrs(char *uplo, char *trans, char *diag, char * normin, integer *n, double *a, integer *lda, double *x, double *scale, double *cnorm, integer *info)
{
  return dlatrs_(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
}
inline integer latrs(char *uplo, char *trans, char *diag, char * normin, integer *n, scomplex *a, integer *lda, scomplex *x, float *scale, float *cnorm, integer *info)
{
  return clatrs_(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
}
inline integer latrs(char *uplo, char *trans, char *diag, char * normin, integer *n, dcomplex *a, integer *lda, dcomplex *x, double *scale, double *cnorm, integer *info)
{
  return zlatrs_(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
}

// --- factors an upper trapezoidal matrix by means of orthogonal transformations ---
inline integer latrz(integer *m, integer *n, integer *l, float *a, integer *lda, float *tau, float *work)
{
  return slatrz_(m, n, l, a, lda, tau, work);
}
inline integer latrz(integer *m, integer *n, integer *l, double *a, integer *lda, double *tau, double *work)
{
  return dlatrz_(m, n, l, a, lda, tau, work);
}
inline integer latrz(integer *m, integer *n, integer *l, scomplex *a, integer *lda, scomplex *tau, scomplex *work)
{
  return clatrz_(m, n, l, a, lda, tau, work);
}
inline integer latrz(integer *m, integer *n, integer *l, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work)
{
  return zlatrz_(m, n, l, a, lda, tau, work);
}

// --- computes a blocked Tall-Skinny QR factorization of a float M-by-N matrix A for M >= N ---
inline integer latsqr(integer *m, integer *n, integer *mb, integer *nb, float *a, integer *lda, float *t, integer *ldt, float *work, integer *lwork, integer *info)
{
  return slatsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline integer latsqr(integer *m, integer *n, integer *mb, integer *nb, double *a, integer *lda, double *t, integer *ldt, double *work, integer *lwork, integer *info)
{
  return dlatsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline integer latsqr(integer *m, integer *n, integer *mb, integer *nb, scomplex *a, integer *lda, scomplex *t, integer *ldt, scomplex *work, integer *lwork, integer *info)
{
  return clatsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}
inline integer latsqr(integer *m, integer *n, integer *mb, integer *nb, dcomplex *a, integer *lda, dcomplex *t, integer *ldt, dcomplex *work, integer *lwork, integer *info)
{
  return zlatsqr_(m, n, mb, nb, a, lda, t, ldt, work, lwork, info);
}

inline float lantp(char* norm, char* uplo, char* diag, integer* n, float* ap, float* work)
{
  return slantp_(norm, uplo, diag, n, ap, work);
}
inline double lantp(char* norm, char* uplo, char* diag, integer* n, double* ap, double* work)
{
  return dlantp_(norm, uplo, diag, n, ap, work);
}
inline float lantp(char* norm, char* uplo, char* diag, integer* n, scomplex* ap, float* work)
{
  return clantp_(norm, uplo, diag, n, ap, work);
}
inline double lantp(char* norm, char* uplo, char* diag, integer* n, dcomplex* ap, double* work)
{
  return zlantp_(norm, uplo, diag, n, ap, work);
}

inline float lantb(char* norm, char* uplo, char* diag, integer* n, integer* k, float* ab, integer* ldab, float* work)
{
  return slantb_(norm, uplo, diag, n, k, ab, ldab, work);
}
inline double lantb(char* norm, char* uplo, char* diag, integer* n, integer* k, double* ab, integer* ldab, double* work)
{
  return dlantb_(norm, uplo, diag, n, k, ab, ldab, work);
}
inline float lantb(char* norm, char* uplo, char* diag, integer* n, integer* k, scomplex* ab, integer* ldab, float* work)
{
  return clantb_(norm, uplo, diag, n, k, ab, ldab, work);
}
inline double lantb(char* norm, char* uplo, char* diag, integer* n, integer* k, dcomplex* ab, integer* ldab, double* work)
{
  return zlantb_(norm, uplo, diag, n, k, ab, ldab, work);
}

integer gelqt(integer* m, integer* n, integer* mb, float* a, integer* lda, float* t, integer* ldt, float* work, integer* info)
{
  return sgelqt_(m, n, mb, a, lda, t, ldt, work, info);
}
integer gelqt(integer* m, integer* n, integer* mb, double* a, integer* lda, double* t, integer* ldt, double* work, integer* info)
{
  return dgelqt_(m, n, mb, a, lda, t, ldt, work, info);
}
integer gelqt(integer* m, integer* n, integer* mb, scomplex* a, integer* lda, scomplex* t, integer* ldt, scomplex* work, integer* info)
{
  return cgelqt_(m, n, mb, a, lda, t, ldt, work, info);
}
integer gelqt(integer* m, integer* n, integer* mb, dcomplex* a, integer* lda, dcomplex* t, integer* ldt, dcomplex* work, integer* info)
{
  return zgelqt_(m, n, mb, a, lda, t, ldt, work, info);
}

integer gemlqt(char* side, char* trans, integer* m, integer* n, integer* k, integer* mb, float* v, integer* ldv, float* t, integer* ldt, float* c, integer* ldc, float* work, integer* info)
{
  return sgemlqt_(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info);
}
integer gemlqt(char* side, char* trans, integer* m, integer* n, integer* k, integer* mb, double* v, integer* ldv, double* t, integer* ldt, double* c, integer* ldc, double* work, integer* info)
{
  return dgemlqt_(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info);
}
integer gemlqt(char* side, char* trans, integer* m, integer* n, integer* k, integer* mb, scomplex* v, integer* ldv, scomplex* t, integer* ldt, scomplex* c, integer* ldc, scomplex* work, integer* info)
{
  return cgemlqt_(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info);
}
integer gemlqt(char* side, char* trans, integer* m, integer* n, integer* k, integer* mb, dcomplex* v, integer* ldv, dcomplex* t, integer* ldt, dcomplex* c, integer* ldc, dcomplex* work, integer* info)
{
  return zgemlqt_(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info);
}

integer getc2(integer* n, float* a, integer* lda, integer* ipiv, integer* jpiv, integer* info)
{
  return sgetc2_(n, a, lda, ipiv, jpiv, info);
}
integer getc2(integer* n, double* a, integer* lda, integer* ipiv, integer* jpiv, integer* info)
{
  return dgetc2_(n, a, lda, ipiv, jpiv, info);
}
integer getc2(integer* n, scomplex* a, integer* lda, integer* ipiv, integer* jpiv, integer* info)
{
  return cgetc2_(n, a, lda, ipiv, jpiv, info);
}
integer getc2(integer* n, dcomplex* a, integer* lda, integer* ipiv, integer* jpiv, integer* info)
{
  return zgetc2_(n, a, lda, ipiv, jpiv, info);
}

// --- LU factorization without pivoting ---
inline integer getrfnp(integer* m, integer* n, float* a, integer* lda, integer* info)
{
  return sgetrfnp_(m, n, a, lda, info);
}
inline integer getrfnp(integer* m, integer* n, double* a, integer* lda, integer* info)
{
  return dgetrfnp_(m, n, a, lda, info);
}
inline integer getrfnp(integer* m, integer* n, scomplex* a, integer* lda, integer* info)
{
  return cgetrfnp_(m, n, a, lda, info);
}
inline integer getrfnp(integer* m, integer* n, dcomplex* a, integer* lda, integer* info)
{
  return zgetrfnp_(m, n, a, lda, info);
}

inline integer spffrt2(float  *ap, integer *n, integer * ncolm, float  *work, float  *work2)
{
  return sspffrt2_(ap, n, ncolm, work, work2);
}
inline integer spffrt2(double *ap, integer *n, integer * ncolm, double *work, double *work2)
{
  return dspffrt2_(ap, n, ncolm, work, work2);
}
inline integer spffrt2(scomplex *ap, integer *n, integer * ncolm, scomplex *work, scomplex *work2)
{
  return cspffrt2_(ap, n, ncolm, work, work2);
}
inline integer spffrt2(dcomplex *ap, integer *n, integer * ncolm, dcomplex *work, dcomplex *work2)
{
  return zspffrt2_(ap, n, ncolm, work, work2);
}

inline integer spffrtx(float *ap, integer *n, integer * ncolm, float *work, float *work2)
{
  return sspffrtx_(ap, n, ncolm, work, work2);
}
inline integer spffrtx(double *ap, integer *n, integer * ncolm, double *work, double *work2)
{
  return dspffrtx_(ap, n, ncolm, work, work2);
}
inline integer spffrtx(scomplex *ap, integer *n, integer * ncolm, scomplex *work, scomplex *work2)
{
  return cspffrtx_(ap, n, ncolm, work, work2);
}
inline integer spffrtx(dcomplex *ap, integer *n, integer * ncolm, dcomplex *work, dcomplex *work2)
{
  return zspffrtx_(ap, n, ncolm, work, work2);
}

// --- LU factorization complete and incomplete without pivoting ---
inline integer getrfnpi(integer *m, integer *n, integer *nfact, float *a, integer *lda, integer *info)
{
  return sgetrfnpi_(m, n, nfact, a, lda, info);
}
inline integer getrfnpi(integer *m, integer *n, integer *nfact, double *a, integer * lda, integer *info)
{
  return dgetrfnpi_(m, n, nfact, a, lda, info);
}
inline integer getrfnpi(integer *m, integer *n, integer *nfact, scomplex *a, integer *lda, integer *info)
{
  return cgetrfnpi_(m, n, nfact, a, lda, info);
}
inline integer getrfnpi(integer *m, integer *n, integer*nfact, dcomplex *a, integer *lda, integer *info)
{
  return zgetrfnpi_(m, n, nfact, a, lda, info);
}

} //namespace libflame
#endif  //  #ifndef LIBFLAME_HH