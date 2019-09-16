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

/*! @file libflame_interface.hh
 *  libflame_interface.hh defines all the Libflame CPP templated public interfaces
 *  */
#ifndef LIBFLAME_INTERFACE_HH
#define LIBFLAME_INTERFACE_HH

#include "libflame.hh"
namespace libflame{
  
  // --- Cholesky factorization --- 
  
  /*! \brief Cholesky factorization of a real symmetric positive definite matrix A
    
    \b Purpose:
    Cholesky factorization of a real symmetric positive definite matrix A.
    The factorization has the form
        A = U**T * U,  if UPLO = 'U', or
        A = L * L**T,  if UPLO = 'L',
     where U is an upper triangular matrix and L is lower triangular.
     
     This is the block version of the algorithm, calling Level 3 BLAS. 

    \param[in] UPLO
    UPLO is char*. \n
    UPLO specifies output format \n
    = 'U': Output is upper triangular factorization of A \n
    = 'L': Output is lower triangular factorization of A
	
    \param[in] N
    N is int*. \n
    The order of the matrix A. N >= 0
	
    \param[in,out] A
    A is float array, dimension (LDA,N). \n
    On entry, the symmetric matrix A.  If UPLO = 'U', the leading
    N-by-N upper triangular part of A contains the upper
    triangular part of the matrix A, and the strictly lower
    triangular part of A is not referenced.  If UPLO = 'L', the
    leading N-by-N lower triangular part of A contains the lower
    triangular part of the matrix A, and the strictly upper
    triangular part of A is not referenced. \n
    On exit, if INFO = 0, the factor U or L from the Cholesky
    factorization A = U**T*U or A = L*L**T.
	
    \param[in] LDA
    LDA is int*. \n
    The leading dimension of the matrix A, LDA >= max(1,N)
    
	\param[out] INFO
    INFO is int*. \n
    = 0:  successful exit \n
    < 0:  if INFO = -i, the i-th argument had an illegal value \n
    \> 0:  if INFO = i, the leading minor of order i is not \n
          positive definite, and the factorization could not be
          completed.
    */ 
  template< typename T >
  void potrf( char* uplo, int* n, T* a, int* lda, int* info )
  {
    potrf(uplo, n, a, lda, info);
  }

  
  template< typename T >
  void potf2( char* uplo, int* n, T* a, int* lda, int* info )
  {
    potf2(uplo, n, a, lda, info );
  }
  // --- LU factorization with partial pivoting ---  


  template< typename T >
  void getrf(int* m, int* n, T* buff_A, int* ldim_A, int* buff_p, int* info )
  {
    getrf(m, n, buff_A, ldim_A, buff_p, info );
  }

  template< typename T >
  void getf2( int* m, int* n, T* a, int* lda, int* ipiv, int* info )
  {
    dgetf2(m, n, a, lda, ipiv, info );
  }
  // --- QR factorization (classic) ---
  
  template< typename T >
  void geqrf( int* m, int* n, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
  {
    geqrf( m, n, a, lda, tau, work, lwork, info );
  }

  template< typename T >
  void geqr2( int* m, int* n, T* a, int* lda, T* tau, T* work, int* info )
  {
    geqr2( m, n, a, lda, tau, work, info );
  }
  
  template< typename T >
  void geqpf( int* m, int* n, T* a, int* lda, int* jpvt, T* tau,T* work, int* info )
  {
    geqpf( m, n, a, lda, jpvt, tau, work, info );
  }
  
  template< typename Ta , typename Tb>
  void geqpf( int* m, int* n, Ta* a, int* lda, int* jpvt, Ta* tau, Ta* work, Tb*  rwork, int* info )
  {
    geqpf( m, n, a, lda, jpvt, tau, work, rwork, info );
  }
  
  template< typename T >
  void geqp3( int* m, int* n, T* a, int* lda, int* jpvt, T* tau, T* work, int* lwork,    int* info )
  {
    geqp3( m, n, a, lda, jpvt, tau, work, lwork, info );
  }
  
  template< typename Ta, typename Tb >
  void F77_cgeqp3( int* m, int* n, Ta* a, int* lda, int* jpvt, Ta* tau, Ta* work, int* lwork, Tb*  rwork, int* info )
  {
    geqp3( m, n, a, lda, jpvt, tau, work, lwork, rwork, info );
  }
  
  // --- LQ factorization (classic) ---
  template< typename T >
  void gelqf( int* m, int* n, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
  {
    gelqf( m, n, a, lda, tau, work, lwork, info );
  }
  
  template< typename T >
  void gelq2( int* m, int* n, T* a, int* lda, T* tau, T* work, int* info )
  {
    gelq2( m, n, a, lda, tau, work, info );
  }
  
  // --- LS solver ---
  template< typename T >
  void gelsd( int* m, int* n, int* nrhs, T* a, int* lda, T* b, int* ldb, T*  s, T*  rcond, int* rank, T* work, int* lwork,    int* iwork, int* info )
  {
    gelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info );
  }
  template< typename Ta, typename Tb >
  void gelsd( int* m, int* n, int* nrhs, Ta* a, int* lda, Ta* b, int* ldb, Tb*  s, Tb*  rcond, int* rank, Ta* work, int* lwork, Tb*  rwork, int* iwork, int* info )
  {
    gelsd( m,  n, nrhs, a, lda, b, ldb, s, rcond, rank,  work, lwork, rwork, iwork, info );
  }
  
  template< typename T >
  void gelss( int* m, int* n, int* nrhs, T* a, int* lda, T* b, int* ldb, T*  s, T*  rcond, int* rank, T* work, int* lwork,    int* info )
  {
    gelss( m, n, nrhs, a, lda, b,ldb, s, rcond, rank, work, lwork, info );  
  }
  
  template< typename Ta, typename Tb >
  void gelss( int* m, int* n, int* nrhs, Ta* a, int* lda, Ta* b, int* ldb, Tb*  s, Tb*  rcond, int* rank, Ta* work, int* lwork, Tb*  rwork, int* info )
  {
    gelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, work, info );
  }
  
  // --- Triangular-transpose matrix multiply ---
  template< typename T >
  void lauum( char* uplo, int* n, T* a, int* lda, int* info )
  {
    lauum( uplo, n, a, lda, info );
  }
  
  template< typename T >
  void lauu2( char* uplo, int* n, T* a, int* lda, int* info )
  {
    lauu2( uplo, n, a, lda, info );
  }
  
  // --- Symmetric (hermitian) positive definite matrix inversion ---
  template< typename T >
  void potri( char* uplo, int*  n, T* buff_A, int*  ldim_A, int*  info )
  {
    potri( uplo, n, buff_A, ldim_A, info );
  }
  
  // --- Triangular matrix inversion ---
  template< typename T >
  void trtri( char* uplo, char* diag, int* n, T* a, int* lda, int* info )
  {
    trtri( uplo, diag, n, a, lda, info );
  }
  
  template< typename T >
  void trti2( char* uplo, char* diag, int* n, T* a, int* lda, int* info )
  {
    trti2( uplo, diag, n, a, lda, info );
  }
  
  // --- Triangular Sylvester equation solve ---
  template< typename T >
  void trsyl( char* transa, char* transb, int* isgn, int* m, int* n, T* a, int* lda, T* b, int* ldb, T* c, int* ldc, T* scale, int* info )
  {
    trsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }
  
  template< typename Ta, typename Tb >
  void trsyl( char* transa, char* transb, int* isgn, int* m, int* n, Ta* a, int* lda, Ta* b, int* ldb, Ta* c, int* ldc, Tb* scale, int* info )
  {
    trsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
  }
  
  // --- Reduction to upper Hessenberg form ---
  template< typename T >
  void sgehrd( int* n, int* ilo, int* ihi, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
  {
    sgehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
  }
  
  template< typename T >
  void gehd2( int* n, int* ilo, int* ihi, T* a, int* lda, T* tau, T* work, int* info )
  {
    gehd2( n, ilo, ihi, a, lda, tau, work, info );
  }
  
  // --- Reduction to tridiagonal form ---
  template< typename T >
  void sytrd( char* uplo, int* n, T* a, int* lda, T*  d, T*  e, T* tau, T* work, int* lwork, int* info )
  {
    sytrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
  }
  
  template< typename T >
  void sytd2( char* uplo, int* n, T* a, int* lda, T*  d, T*  e, T* tau, int* info )
  {
    sytd2( uplo, n, a, lda, d, e, tau, info );
  }
  
  template< typename Ta, typename Tb >
  void hetd2( char* uplo, int* n, Ta* a, int* lda, Tb*  d, Tb*  e, Ta* tau, int* info )
  {
    hetd2( uplo, n, a, lda, d, e, tau, info );
  }
  
  // --- Reduction to bidiagonal form ---
  template< typename T >
  void gebrd( int* m, int* n, T* a, int* lda, T*  d, T*  e, T* tauq, T* taup, T* work, int* lwork, int* info )
  {
    gebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }
  
  template< typename Ta, typename Tb >
  void gebrd( int* m, int* n, Ta* a, int* lda, Tb*  d, Tb*  e, Ta* tauq, Ta* taup, Ta* work, int* lwork, int* info )
  {
    gebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
  }
  
  template< typename T >
  void gebd2( int* m, int* n, T* a, int* lda, T*  d, T*  e, T* tauq, T* taup, T* work, int* info )
  {
    gebd2( m, n, a, lda, d, e, tauq, taup, work, info );
  }
  
  template< typename Ta, typename Tb >
  void gebd2( int* m, int* n, Ta* a, int* lda, Tb*  d, Tb*  e, Ta* tauq, Ta* taup, Ta* work, int* info )
  {
    gebd2( m, n, a,lda, d, e, tauq, taup, work, info );
  }
  
  // --- Reduce Hermitian-definite generalized eigenproblem to standard form ---
  
  template< typename T >
  void sygst( int* itype, char* uplo, int* n, T* a, int* lda, T* b, int* ldb, int* info )
  {
    sygst( itype, uplo, n, a, lda, b, ldb, info );
  }
  
  template< typename T >
  void sygs2( int* itype, char* uplo, int* n, T* a, int* lda, T* b, int* ldb, int* info )
  {
    sygs2( itype, uplo, n, a, lda, b, ldb, info );
  }
  
  // --- Accumulate block Householder matrix T (classic) ---
  template< typename T >
  void larft( char* direct, char* storev, int* n, int* k, T* v, int* ldv, T* tau, T* t, int* ldt )
  {
    larft( direct, storev, n, k, v, ldv, tau, t, ldt );
  }
  
  // --- Generate a Householder vector (classic) ---
  template< typename T >
  void larfg( int* n, T* alpha, T* x, int* incx, T* tau )
  {
    larfg( n, alpha, x, incx,tau );
  }
  
  template< typename T >
  void larfgp( int* n, T* alpha, T* x, int* incx, T* tau )
  {
    larfgp( n,  alpha, x, incx, tau );
  }
  
  // --- Form Q from QR factorization ---
  template< typename T >
  void orgqr( int* m, int* n, int* k, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
  {
    orgqr(  m, n, k,  a, lda, tau, work, lwork, info );
  }
  
  // --- Apply Q or Q' from QR factorization ---
  template< typename T >
  void ormqr( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
  {
    ormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  
  template< typename T >
  void orm2r( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* info )
  {
    orm2r( side, trans,  m, n, k, a, lda, tau, c, ldc, work, info );  
  }
  
  // --- Form Q from LQ factorization ---
  template< typename T >
  void orglq( int* m, int* n, int* k, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
  {
    orglq( m, n, k, a, lda, tau, work, lwork, info );
  }
  
  // --- Apply Q or Q' from LQ factorization ---
  template< typename T >
  void ormlq( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
  {
    ormlq( side, trans,  m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  
  template< typename T >
  void orml2( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* info )
  {
    orml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
  }
  
  // --- Form Q from tridiagonal reduction ---
  template< typename T >
  void orgtr( char* uplo, int* m, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
  {
    orgtr( uplo, m, a, lda, tau, work, lwork, info );
  }
  
  // --- Apply Q or Q' from tridiagonal reduction ---
  template< typename T >
  void ormtr( char* side, char* uplo, char* trans, int* m, int* n, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
  {
    ormtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
  }
  
  // --- Form Q from bidiagonal reduction ---
  template< typename T >
  void sorgbr( char* vect, int* m, int* n, int* k, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
  {
    sorgbr( vect, m, n, k, a, lda, tau, work, lwork, info );
  }
  
  // --- Apply Q or Q' from bidiagonal reduction ---
  template< typename T >
  void sormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, T* a, int* lda,T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
  {
    sormbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
  }
  
  // --- Tridiagonal QR algorithm ---
  template< typename T >
  void ssteqr( char* jobz, int* n, T* d, T* e, T* z, int* ldz, T*  work, int* info ) 
  {
    ssteqr( jobz, n, d, e, z, ldz, work, info ); 
  }
  
  template< typename Ta, typename Tb >
  void csteqr( char* jobz, int* n, Ta* d, Ta* e, Tb* z, int* ldz, Ta*  work, int* info )
  {
    csteqr( jobz, n, d, e, z, ldz, work, info ); 
  }
  
  // --- Tridiagonal divide-and-conquer algorithm ---
  template< typename T >
  void sstedc( char* compz, int* n, T* d, T* e, T* z, int* ldz, T* work, int* lwork,        int* iwork, int* liwork, int* info )
  {
     sstedc( compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info );
  }
  template< typename Ta, typename Tb >
  void cstedc( char* compz, int* n, Ta* d, Ta* e, Tb* z, int* ldz, Tb* work, int* lwork, Ta*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    cstedc( compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info );
  }
  
  // --- Tridiagonal MRRR algorithm ---
  template< typename T >
  void sstemr( char* jobz, char* range, int* n, T*  d, T*  e, int* vl, int* vu, int* il, int* iu, int* m, T*  w, T* z, int* ldz, int* nzc, int* isuppz, int* tryrac, T*  work, int* lwork, int* iwork, int* liwork, int* info )
  {
    sstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }
  
  template< typename Ta, typename Tb >
  void cstemr( char* jobz, char* range, int* n, Ta*  d, Ta*  e, int* vl, int* vu, int* il, int* iu, int* m, Ta*  w, Tb* z, int* ldz, int* nzc, int* isuppz, int* tryrac, Ta*  work, int* lwork, int* iwork, int* liwork, int* info )
  {
    cstemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
  }
  
  
  // --- Hermitian eigenvalue decomposition (QR algorithm) ---
  template< typename T >
  void ssyev( char* jobz, char* uplo, int* n, T* a, int* lda, T*  w, T* work, int* lwork, T*  rwork, int* info )
  {
    ssyev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info );
  }
  
  template< typename Ta, typename Tb >
  void cheev( char* jobz, char* uplo, int* n, Ta* a, int* lda, Tb*  w, Ta* work, int* lwork, Tb*  rwork, int* info )
  {
     cheev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info ); 
  }
  
  // --- Hermitian eigenvalue decomposition (divide-and-conquer) ---
  template< typename T >
  void ssyevd( char* jobz, char* uplo, int* n, T* a, int* lda, T*  w, T* work, int* lwork,        int* iwork, int* liwork, int* info )
  {
    ssyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info ); 
  }
  
  template< typename Ta, typename Tb >
  void cheevd( char* jobz, char* uplo, int* n, Ta* a, int* lda, Tb*  w, Ta* work, int* lwork, Tb*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    cheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info ); 
  }
  
  // --- Hermitian eigenvalue decomposition (MRRR) ---
  template< typename T >
  void ssyevr( char* jobz, char* range, char* uplo, int* n, T* a, int* lda, T*  vl, T*  vu, int* il, int* iu, T*  abstol, int* m, T*  w, T* z, int* ldz, int* isuppz, T* work, int* lwork,        int* iwork, int* liwork, int* info )
  {
    ssyevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info ); 
  }
  
  template< typename Ta, typename Tb >
  void cheevr( char* jobz, char* range, char* uplo, int* n, Ta* a, int* lda, Tb*  vl, Tb*  vu, int* il, int* iu, Tb*  abstol, int* m, Tb*  w, Ta* z, int* ldz, int* isuppz, scomplex* work, int* lwork, Tb*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
  {
    cheevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info ); 
  }
  
  // --- Bidiagonal QR algorithm ---
  template< typename T >
  void sbdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, T* d, T* e, T* vt, int* ldvt, T* u, int* ldu, T* c, int* ldc, T*  rwork, int* info )
  {
     sbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info ); 
  }
  
  template< typename Ta, typename Tb >
  void cbdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, Ta* d, Ta* e, Tb* vt, int* ldvt, Tb* u, int* ldu, Tb* c, int* ldc, Ta*  rwork, int* info )
  {
     cbdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info ); 
  }
  
  // --- Bidiagonal divide-and-conquor algorithm ---
  template< typename T >
  void sbdsdc( char* uplo, char* compq, int* n, T*  d, T*  e, T*  u, int* ldu, T*  vt, int* ldvt, T*  q, T*  iq, T*  work, int* iwork, int* info )
  {
    sbdsdc( uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info ); 
  }
  
  // --- General matrix singular value decomposition (QR algorithm) ---
  template< typename T >
  void sgesvd( char* jobu, char* jobv, int* m, int* n, T* a, int* lda, T*  s, T* u, int* ldu, T* vt, int* ldvt, T* work, int* lwork,    int* info )
  {
    sgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info );
  }
  
  template< typename Ta, typename Tb >
  void cgesvd( char* jobu, char* jobv, int* m, int* n, Ta* a, int* lda, Tb*  s, Ta* u, int* ldu, Ta* vt, int* ldvt, Ta* work, int* lwork, Tb*  rwork, int* info )
  {
    cgesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info );
  }
  
  // --- General matrix singular value decomposition (divide-and-conquer) ---
  template< typename T >
  void sgesdd( char* jobz, int* m, int* n, T* a, int* lda, T*  s, T* u, int* ldu, T* vt, int* ldvt, T* work, int* lwork,    int* iwork, int* info )
  {
    sgesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info );
  }
  
  template< typename Ta, typename Tb >
  void cgesdd( char* jobz, int* m, int* n, Ta* a, int* lda, Tb*  s, Ta* u, int* ldu, Ta* vt, int* ldvt, Ta* work, int* lwork, Tb*  rwork, int* iwork, int* info )
  {
    cgesdd( jobz, m, n, a, lda, s, u, ldu, vt, vt, work, lwork, rwork, iwork, info );
  }
  
  // --- Swap rows ---
  template< typename T >
  void slaswp( int* n, T* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
  {
    slaswp( n, a, lda, k1, k2, ipiv, incx );
  }
  
  // --- Initialize a matrix ---
  template< typename T >
  void slaset( char* uplo, int* m, int* n, T* alpha, T* beta, T* a, int* lda )
  {
    slaset( uplo, m, n, alpha, beta, a, lda );
  }
}  // namespace libflame
#endif  //  #ifndef INTERFACE_HH


