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

/*! @brief Cholesky factorization of a real symmetric positive definite matrix a
    *
    * @details
    * \b Purpose:
    * \verbatim
        Cholesky factorization of a real symmetric positive definite matrix a.
        The factorization has the form
            A = U**T * U,  if uplo = 'U', or
            A = L * L**T,  if uplo = 'L',
        where U is an upper triangular matrix and L is lower triangular.

        This is the block version of the algorithm, calling Level 3 BLAS.
    \endverbatim

    \param[in] uplo
    uplo is char*. \n
    uplo specifies output format \n
    = 'U': Output is upper triangular factorization of A \n
    = 'L': Output is lower triangular factorization of A

    * @param[in] n
    n is int*. \n
    The order of the matrix a. n >= 0

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n). \n
    On entry, the symmetric matrix a.  If uplo = 'U', the leading
    n-by-n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced.  If uplo = 'L', the
    leading n-by-n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced. \n
    On exit, if info = 0, the factor U or L from the Cholesky
    factorization A = U**T*U or A = L*L**T.

    * @param[in] lda
    lda is int*. \n
    The leading dimension of the matrix a, lda >= max(1,n)

    * @param[out] info
    info is int*. \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  if info = i, the leading minor of order i is not \n
        positive definite, and the factorization could not be
        completed.
    *  */
template< typename T >
int potrf( char* uplo, int* n, T* a, int* lda, int* info )
{
  return potrf(uplo, n, a, lda, info);
}

/*! @brief Cholesky factorization of a real symmetric
    *         positive definite matrix a
    *
    * @details
    * \b Purpose:
    * \verbatim
        Cholesky factorization of a real symmetric positive definite matrix a
        The factorization has the form
            A = U**T * U,  if uplo = 'U', or
            A = L * L**T,  if uplo = 'L',
        where U is an upper triangular matrix and L is lower triangular.

        This is the unblocked version of the algorithm, calling Level 2 BLAS.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    uplo specifies output format \n
    = 'U': Output is upper triangular factorization of A \n
    = 'L': Output is lower triangular factorization of A

    * @param[in] n
    n is int* \n
    The order of the matrix a. n >= 0

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the symmetric matrix a.  If uplo = 'U', the leading
    n-by-n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced.  If uplo = 'L', the
    leading n-by-n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced. \n
    On exit, if info = 0, the factor U or L from the Cholesky
    factorization A = U**T *U  or A = L*L**T.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,n)

    @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -k, the k-th argument had an illegal value \n
    \> 0:  if info = k, the leading minor of order k is not
        positive definite, and the factorization could not be
        completed.
    *  */
template< typename T >
int potf2( char* uplo, int* n, T* a, int* lda, int* info )
{
  return potf2(uplo, n, a, lda, info );
}

/*! @brief LU factorization of a general m-by-n matrix a
    *   using partial pivoting with row interchanges.
    *
    *  @details
    * \b Purpose:
    * \verbatim
        LU factorization of a general m-by-n matrix a using partial pivoting with row interchanges.
        The factorization has the form
            A = P * L * U
        where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower
        trapezoidal if M > N), and U is upper triangular (upper trapezoidal if M < N).

        This is the right-looking Level 3 BLAS version of the algorithm.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix to be factored. \n
    On exit, the factors L and U from the factorization
    A = P*L*U; the unit diagonal elements of L are not stored.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[out] ipiv
    ipiv is int array, dimension (min(m,n)) \n
    The pivot indices; for 1 <= i <= min(m,n), row i of the
    matrix was interchanged with row ipiv(i).

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  if info = i, U(i,i) is exactly zero. The factorization
        has been completed, but the factor U is exactly
        singular, and division by zero will occur if it is used
        to solve a system of equations.
    *  */
template< typename T >
int getrf( int* m, int* n, T* a, int* lda, int* ipiv, int* info )
{
  return getrf(m, n, a, lda, ipiv, info );
}


/*! @brief LU factorization of a general m-by-n matrix a
    *   using partial pivoting with row interchanges.
    *
    * @details
    * \b Purpose:
    * \verbatim
        LU factorization of a general m-by-n matrix a using partial pivoting with row interchanges.
        The factorization has the form
            A = P * L * U
        where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower
        trapezoidal if M > N), and U is upper triangular (upper trapezoidal if M < N).

        This is the right-looking Level 2 BLAS version of the algorithm.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix to be factored. \n
    On exit, the factors L and U from the factorization
    A = P*L*U; the unit diagonal elements of L are not stored.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[out] ipiv
    ipiv is int array, dimension (min(m,n)) \n
    ipiv is int array, dimension (min(m,n)) \n
    The pivot indices; for 1 <= i <= min(m,n), row i of the
    matrix was interchanged with row ipiv(i).

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -k, the k-th argument had an illegal value \n
    \> 0:  if info = k, U(k,k) is exactly zero. The factorization
        has been completed, but the factor U is exactly
        has been completed, but the factor U is exactly
        singular, and division by zero will occur if it is used
        singular, and division by zero will occur if it is used
        to solve a system of equations.
        to solve a system of equations.
    *  */
template< typename T >
int getf2( int* m, int* n, T* a, int* lda, int* ipiv, int* info )
{
  return getf2(m, n, a, lda, ipiv, info );
}

/*! @brief QR factorization of a real m-by-n matrix a
    *
    * @details
    * \b Purpose:
    * \verbatim
        QR factorization of a real m-by-n matrix a
        The factorization has the form
            A = Q * R
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix to be factored. \n
    On exit, the elements on and above the diagonal of the array
    contain the min(m,n)-by-n upper trapezoidal matrix R (R is
    upper triangular if m >= n); the elements below the diagonal,
    with the array tau, represent the orthogonal matrix Q as a
    product of min(m,n) elementary reflectors (see Further
    Details).

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[out] tau
    tau is float/double/scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors (see Further
    Details).

    * @param[out] work
    work is float/double/scomplex/dcomplex array, dimension (min(m,n)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.
    Details).

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work.  lwork >= max(1,n). \n
    For optimum performance lwork >= n*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of elementary reflectors
            Q = H(1) H(2) . . . H(k), where k = min(m,n).

            Each H(i) has the form
            H(i) = I - tau * V * V**T

            where, tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i) = 1;
            V(i+1:M) is stored on exit in A(i+1:M,i), and tau in tau(i).
    \endverbatim
    *  */
template< typename T >
int geqrf( int* m, int* n, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return geqrf( m, n, a, lda, tau, work, lwork, info );
}

/*! @brief QR factorization of a real m-by-n matrix a
    *
    * @details
    * \b Purpose:
    * \verbatim
        QR factorization of a real m-by-n matrix a
        The factorization has the form
            A = Q * R
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix to be factored. \n
    On exit, the elements on and above the diagonal of the array
    contain the min(m,n)-by-n upper trapezoidal matrix R (R is
    upper triangular if m >= n); the elements below the diagonal,
    with the array tau, represent the orthogonal matrix Q as a
    product of elementary reflectors (see Further Details).

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[out] tau
    tau is float/double/scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors (see Further
    Details).

    * @param[out] work
    work is float/double/scomplex/dcomplex array, dimension (n)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of elementary reflectors
            Q = H(1) H(2) . . . H(k), where k = min(m,n).

            Each H(i) has the form
            H(i) = I - tau * V * V**T

            where, tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i) = 1;
            V(i+1:M) is stored on exit in A(i+1:M,i), and tau in tau(i).
    \endverbatim
    *  */
template< typename T >
int geqr2( int* m, int* n, T* a, int* lda, T* tau, T* work, int* info )
{
  return geqr2( m, n, a, lda, tau, work, info );
}


/*! @brief QR factorization of a real m-by-n matrix a
    *   This routine is deprecated and has been
    *   replaced by routine SGEQP3.
    *

    * @details
    * \b Purpose:
    * \verbatim
        QR factorization with column pivoting of a real m-by-n matrix a.
        This routine is deprecated and has been replaced by routine SGEQP3.
        The factorization has the form:
            A*P = Q*R.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the m-by-n matrix a. \n
    On exit, the upper triangle of the array contains the
    min(m,n)-by-n upper triangular matrix R; the elements
    below the diagonal, together with the array tau,
    represent the orthogonal matrix Q as a product of
    min(m,n) elementary reflectors.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[in,out] jpvt
    jpvt is int array, dimension (n) \n
    On entry, if jpvt(i) .ne. 0, the i-th column of A is permuted
    to the front of A*P (a leading column); if jpvt(i) = 0,
    the i-th column of A is a free column. \n
    On exit, if jpvt(i) = k, then the i-th column of A*P
    was the k-th column of A.

    * @param[out] tau
    tau is float/double array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors.

    * @param[out] work
    work is float/double, dimension (3*n)

    * @param[out] info
    info is int*
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of elementary reflectors

            Q = H(1) H(2) . . . H(n)

            Each H(i) has the form

            H = I - tau * V * V**T

            where tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i) = 1;
            V(i+1:m) is stored on exit in A(i+1:m,i).

            The matrix P is represented in jpvt as follows:
            If jpvt(j) = i,
            then the jth column of P is the ith canonical unit vector.

            Partial column norm updating strategy modified by
            Z. Drmac and Z. Bujanovic, Dept. of Mathematics,
            University of Zagreb, Croatia.
            -- April 2011 --
            For more details see LAPACK Working Note 176.
    \endverbatim
    *  */
template< typename T >
int geqpf( int* m, int* n, T* a, int* lda, int* jpvt, T* tau,T* work, int* info )
{
  return geqpf( m, n, a, lda, jpvt, tau, work, info );
}

/*! @brief QR factorization of a real m-by-n matrix a
    *         This routine is deprecated and has been
    *         replaced by routine SGEQP3.
    *

    * @details
    * \b Purpose:
    * \verbatim
        QR factorization with column pivoting of a real m-by-n matrix a:
        This routine is deprecated and has been replaced by routine SGEQP3.
        The factorization has the form:
            A*P = Q*R.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix a. \n
    On exit, the upper triangle of the array contains the
    min(m,n)-by-n upper triangular matrix R; the elements
    below the diagonal, together with the array tau,
    represent the orthogonal matrix Q as a product of
    min(m,n) elementary reflectors.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(0,m)

    * @param[in,out] jpvt
    jpvt is int array, dimension (n) \n
    On entry, if jpvt(i) .ne. 0, the i-th column of A is permuted
    to the front of A*P (a leading column); if jpvt(i) = 0,
    the i-th column of A is a free column. \n
    On exit, if jpvt(i) = k, then the i-th column of A*P
    was the k-th column of A.

    * @param[out] tau
    tau is scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (3*n)

    * @param[out] rwork
    rwork is float/double array, dimension (2*n)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of elementary reflectors

                Q = H(1) H(2) . . . H(n)

            Each H(i) has the form

                H = I - tau * V * V**H

            where tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i) = 1;
            V(i+1:m) is stored on exit in A(i+1:m,i).

            The matrix P is represented in jpvt as follows:
            If jpvt(j) = i,
            then the jth column of P is the ith canonical unit vector.

            Partial column norm updating strategy modified by
                Z. Drmac and Z. Bujanovic, Dept. of Mathematics,
                University of Zagreb, Croatia.
            -- April 2011 --
            For more details see LAPACK Working Note 176.
    \endverbatim
    *  */
template< typename Ta, typename Tb >
int geqpf( int* m, int* n, Ta* a, int* lda, int* jpvt, Ta* tau, Ta* work, Tb* rwork, int* info )
{
  return geqpf( m, n, a, lda, jpvt, tau, work, rwork, info );
}

/*! @brief QR factorization of a real m-by-n matrix a
    *

    * @details
    * \b Purpose:
    * \verbatim
        QR factorization with column pivoting of a real m-by-n matrix a:
            A*P = Q*R
        using Level 3 BLAS.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the m-by-n matrix a. \n
    On exit, the upper triangle of the array contains the
    min(m,n)-by-n upper triangular matrix R; the elements
    below the diagonal, together with the array tau,
    represent the orthogonal matrix Q as a product of
    min(m,n) elementary reflectors.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[in,out] jpvt
    jpvt is int array, dimension (n) \n
    On entry, if jpvt(J).ne.0, the J-th column of A is permuted
    to the front of A*P (a leading column); if jpvt(J) = 0,
    the J-th column of A is a free column. \n
    On exit, if jpvt(J) = K, then the J-th column of A*P
    was the K-th column of A.

    * @param[out] tau
    tau is float/double array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors.

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info=0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= 3*n+1. \n
    For optimal performance lwork >= 2*n+( n+1 )*NB, where NB
    is the optimal blocksize. \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of elementary reflectors

            Q = H(1) H(2) . . . H(k), where k = min(m,n).

            Each H(i) has the form

            H = I - tau * V * V**T

            where tau is a real scalar, and V is a real/complex vector with V(1:i-1) = 0 and V(i) = 1;
            V(i+1:M) is stored on exit in A(i+1:M,i), and tau in tau(i)
    \endverbatim
    *  */
template< typename T >
int geqp3( int* m, int* n, T* a, int* lda, int* jpvt, T* tau, T* work, int* lwork,  int* info )
{
  return geqp3( m, n, a, lda, jpvt, tau, work, lwork, info );
}


/*! @brief QR factorization of a real m-by-n matrix a
    *

    * @details
    * \b Purpose:
    * \verbatim
        QR factorization with column pivoting of a real m-by-n matrix a:
            A*P = Q*R
        using Level 3 BLAS.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix a. \n
    On exit, the upper triangle of the array contains the
    min(m,n)-by-n upper triangular matrix R; the elements
    below the diagonal, together with the array tau,
    represent the orthogonal matrix Q as a product of
    min(m,n) elementary reflectors.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[in,out] jpvt
    jpvt is int array, dimension (n) \n
    On entry, if jpvt(J).ne.0, the J-th column of A is permuted
    to the front of A*P (a leading column); if jpvt(J) = 0,
    the J-th column of A is a free column. \n
    On exit, if jpvt(J) = K, then the J-th column of A*P
    was the K-th column of A.

    * @param[out] tau
    tau is scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info=0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= 3*n+1. \n
    For optimal performance lwork >= 2*n+( n+1 )*NB, where NB
    is the optimal blocksize. \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] rwork
    rwork is float/double array, dimension (2*n)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of elementary reflectors

                Q = H(1) H(2) . . . H(k), where k = min(m,n).

            Each H(i) has the form

                H = I - tau * V * V**H

            where tau is a complex scalar, and V is a real/complex vector with V(1:i-1) = 0 and
        V(i) = 1; V(i+1:M) is stored on exit in A(i+1:M,i), and tau in tau(i)
    \endverbatim
    *  */
template< typename Ta, typename Tb >
int geqp3( int* m, int* n, Ta* a, int* lda, int* jpvt, Ta* tau, Ta* work, int* lwork, Tb* rwork, int* info )
{
  return geqp3( m, n, a, lda, jpvt, tau, work, lwork, rwork, info );
}

/*! @brief LQ factorization of a real m-by-n matrix a
    *

    * @details
    * \b Purpose:
    * \verbatim
        LQ factorization of a real m-by-n matrix a
        The factorization has the form
            A = L * Q
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix. \n
    On exit, the elements on and below the diagonal of the array
    contain the m-by-min(m,n) lower trapezoidal matrix L (L is
    lower triangular if m <= n); the elements above the diagonal,
    with the array tau, represent the orthogonal matrix Q as a
    product of elementary reflectors (see Further Details).

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[out] tau
    tau is float/double/scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors (see Further Details).

    * @param[out] work
    work is float/double/scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int*
    The dimension of the array work.  lwork >= max(1,m). \n
    For optimum performance lwork >= m*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of elementary reflectors
            Q = H(k) . . . H(2) H(1), where k = min(m,n).

            Each H(i) has the form

            H(i) = I - tau * V * V**T

            where tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i) = 1;
            V(i+1:N) is stored on exit in A(i,i+1:N), and tau in tau(i).
    \endverbatim
    *  */
template< typename T >
int gelqf( int* m, int* n, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return gelqf( m, n, a, lda, tau, work, lwork, info );
}

/*! @brief LQ factorization of a real m-by-n matrix a
    *

    * @details
    * \b Purpose:
    * \verbatim
        LQ factorization of a real m-by-n matrix a
        The factorization has the form
        A = Q * R
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix. \n
    On exit, the elements on and below the diagonal of the array
    contain the m by min(m,n) lower trapezoidal matrix L (L is
    lower triangular if m <= n); the elements above the diagonal,
    with the array tau, represent the orthogonal matrix Q as a
    product of elementary reflectors (see Further Details).

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[out] tau
    tau is float/double/scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors (see Further
    Details).

    * @param[out] work
    work is float/double/scomplex/dcomplex array, dimension (m)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of elementary reflectors
            Q = H(k) . . . H(2) H(1), where k = min(m,n).

            Each H(i) has the form
            H(i) = I - tau * V * V**T

            where tau is a real scalar, and V is a real vector with V(1:i-1) = 0 and V(i) = 1;
            V(i+1:N) is stored on exit in A(i,i+1:N), and tau in tau(i).
    \endverbatim
    *  */
template< typename T >
int gelq2( int* m, int* n, T* a, int* lda, T* tau, T* work, int* info )
{
  return gelq2( m, n, a, lda, tau, work, info );
}

/*! @brief The minimum-norm solution to a real linear least squares problem
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of the minimum-norm solution to a real linear least squares problem:
        minimize 2-norm(| B - A*X |)
        using the singular value decomposition (SVD) of A. A is an m-by-n matrix which may be
        rank-deficient.

        Several right hand side vectors B and solution vectors X can be handled in a single call;
        they are stored as the columns of the m-by-nrhs  right hand side matrix b and the n-by-nrhs
        solution matrix X.

        The problem is solved in three steps:
        (1) Reduce the coefficient matrix a to bidiagonal form with Householder transformations,
        reducing the original problem into a "bidiagonal least squares problem" (BLS)
        (2) Solve the BLS using a divide and conquer approach.
        (3) Apply back all the Householder transformations to solve the original least squares
        problem.

        The effective rank of A is determined by treating as zero those singular values which are
        less than rcond times the largest singular value.

        The divide and conquer algorithm makes very mild assumptions about floating point
        arithmetic. It will work on machines with a guard digit in add/subtract, or on those
        binary machines without guard digits which subtract like the Cray X-MP, Cray Y-MP,
        Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without
        guard digits, but we know of none.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in] nrhs
    nrhs is int* \n
    The number of right hand sides, i.e., the number of columns
    of the matrices b and X. nrhs >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n)
    On entry, the m-by-n matrix. \n
    On exit, A has been destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[in,out] b
    b is float/double array, dimension (ldb,nrhs) \n
    On entry, the m-by-n matrix. \n
    On exit, b is overwritten by the n-by-nrhs solution
    matrix X.  If m >= n and rank = n, the residual
    sum-of-squares for the solution in the i-th column is given
    by the sum of squares of elements n+1:m in that column.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the matrix b, ldb >= max(1,max(m,n)).

    * @param[out] s
    s is float/double array, dimension (min(m,n)) \n
    The singular values of A in decreasing order. \n
    The condition number of A in the 2-norm = S(1)/S(min(m,n)).

    * @param[in] rcond
    rcond is float/double* \n
    rcond is used to determine the effective rank of A. \n
    Singular values S(i) <= rcond*S(1) are treated as zero. \n
    If rcond < 0, machine precision is used instead.

    * @param[out] rank
    rank is int* \n
    The effective rank of A, i.e., the number of singular values
    which are greater than rcond*S(1).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork must be at least 1. \n
    The exact minimum amount of workspace needed depends on M,
    N and nrhs. As long as lwork is at least
    12*N + 2*N*SMLSIZ + 8*N*NLVL + N*nrhs + (SMLSIZ+1)**2, \n
    if M is greater than or equal to N or \n
    12*M + 2*M*SMLSIZ + 8*M*NLVL + M*nrhs + (SMLSIZ+1)**2, \n
    if M is less than N, the code will execute correctly. \n
    SMLSIZ is returned by ILAENV and is equal to the maximum
    size of the subproblems at the bottom of the computation
    tree (usually about 25), and \n
    NLVL = max( 0, INT( LOG_2( MIN( m,n )/(SMLSIZ+1) ) ) + 1 ) \n
    For good performance, lwork should generally be larger. \n \n

    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the array work and the
    minimum size of the array iwork, and returns these values as
    the first entries of the work and iwork arrays, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] iwork
    iwork is int array, dimension (max(1,liwork)) \n
    liwork >= max(1, 3*MINMN*NLVL + 11*MINMN), \n
    where MINMN = MIN( m,n ). \n
    On exit, if info = 0, iwork(1) returns the minimum liwork.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  the algorithm for computing the SVD failed to converge;
    if info = i, i off-diagonal elements of an intermediate
    bidiagonal form did not converge to zero.
    *  */
template< typename T >
int gelsd( int* m, int* n, int* nrhs, T* a, int* lda, T* b, int* ldb, T* s, T* rcond, int* rank, T* work, int* lwork, int* iwork, int* info )
{
  return gelsd( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info );
}

/*! @brief The minimum-norm solution to a real linear least squares problem
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of the minimum-norm solution to a real linear least squares problem:
            minimize 2-norm(| B - A*X |)
        using the singular value decomposition (SVD) of A. A is an m-by-n matrix which may be
        rank-deficient.

        Several right hand side vectors B and solution vectors X can be handled in a single call;
        they are stored as the columns of the m-by-nrhs  right hand side matrix b and the n-by-nrhs
        solution matrix X.

        The problem is solved in three steps:
        (1) Reduce the coefficient matrix a to bidiagonal form with Householder transformations,
        reducing the original problem into a "bidiagonal least squares problem" (BLS)
        (2) Solve the BLS using a divide and conquer approach.
        (3) Apply back all the Householder transformations to solve the original least squares
        problem.

        The effective rank of A is determined by treating as zero those singular values which are
        less than rcond times the largest singular value.

        The divide and conquer algorithm makes very mild assumptions about floating point
        arithmetic. It will work on machines with a guard digit in add/subtract, or on those
        binary machines without guard digits which subtract like the Cray X-MP, Cray Y-MP,
        Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without
        guard digits, but we know of none.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in] nrhs
    nrhs is int* \n
    The number of right hand sides, i.e., the number of columns
    of the matrices b and X. nrhs >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix. \n
    On exit, A has been destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[in,out] b
    b is scomplex/dcomplex array, dimension (ldb,nrhs) \n
    On entry, the m-by-n matrix. \n
    On exit, b is overwritten by the n-by-nrhs solution
    matrix X.  If m >= n and rank = n, the residual
    sum-of-squares for the solution in the i-th column is given
    by the sum of squares of elements n+1:m in that column.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the matrix b, ldb >= max(1,max(m,n)).

    * @param[out] s
    s is float/double array, dimension (min(m,n)) \n
    The singular values of A in decreasing order. \n
    The condition number of A in the 2-norm = S(1)/S(min(m,n)).

    * @param[in] rcond
    rcond is float/double \n
    rcond is used to determine the effective rank of A. \n
    Singular values S(i) <= rcond*S(1) are treated as zero. \n
    If rcond < 0, machine precision is used instead.

    * @param[out] rank
    rank is int* \n
    The effective rank of A, i.e., the number of singular values
    which are greater than rcond*S(1).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork must be at least 1.
    The exact minimum amount of workspace needed depends on M,
    N and nrhs. As long as lwork is at least \n
        12*N + 2*N*SMLSIZ + 8*N*NLVL + N*nrhs + (SMLSIZ+1)**2, \n
    if M is greater than or equal to N or \n
        12*M + 2*M*SMLSIZ + 8*M*NLVL + M*nrhs + (SMLSIZ+1)**2, \n
    if M is less than N, the code will execute correctly. \n
    SMLSIZ is returned by ILAENV and is equal to the maximum
    size of the subproblems at the bottom of the computation
    tree (usually about 25), and \n
    NLVL = max( 0, INT( LOG_2( MIN( m,n )/(SMLSIZ+1) ) ) + 1 ) \n
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the array work and the
    minimum size of the array iwork, and returns these values as
    the first entries of the work and iwork arrays, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] rwork
    rwork is float/double array, dimension (max(1,lrwork)) \n
    lrwork >= \n
    10*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*nrhs + \n
    max( (SMLSIZ+1)**2, N*(1+nrhs) + 2*nrhs ) \n
    if M is greater than or equal to N or \n
    10*M + 2*M*SMLSIZ + 8*M*NLVL + 3*SMLSIZ*nrhs + \n
    max( (SMLSIZ+1)**2, N*(1+nrhs) + 2*nrhs ) \n
    if M is less than N, the code will execute correctly. \n
    SMLSIZ is returned by ILAENV and is equal to the maximum
    size of the subproblems at the bottom of the computation
    tree (usually about 25), and \n
    NLVL = max( 0, INT( LOG_2( MIN( m,n )/(SMLSIZ+1) ) ) + 1 ) \n
    On exit, if info = 0, rwork(1) returns the minimum lrwork.

    * @param[out] iwork
    iwork is int array, dimension (max(1,liwork)) \n
    liwork >= max(1, 3*MINMN*NLVL + 11*MINMN), \n
    where MINMN = MIN( m,n ). \n
    On exit, if info = 0, iwork(1) returns the minimum liwork.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  the algorithm for computing the SVD failed to converge;
        if info = i, i off-diagonal elements of an intermediate
        bidiagonal form did not converge to zero.
    *  */
template< typename Ta, typename Tb >
int gelsd( int* m, int* n, int* nrhs, Ta* a, int* lda, Ta* b, int* ldb, Tb* s, Tb*  rcond, int* rank, Ta* work, int* lwork, Tb*  rwork, int* iwork, int* info )
{
  return gelsd( m,  n, nrhs, a, lda, b, ldb, s, rcond, rank,  work, lwork, rwork, iwork, info );
}

/*! @brief The minimum-norm solution to a real linear least squares problem
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of the minimum-norm solution to a real linear least squares problem:
        minimize 2-norm(| B - A*X |)
        using the singular value decomposition (SVD) of A. A is an m-by-n matrix which may be
        rank-deficient.

        Several right hand side vectors B and solution vectors X can be handled in a single call;
        they are stored as the columns of the m-by-nrhs  right hand side matrix b and the n-by-nrhs
        solution matrix X.

        The effective rank of A is determined by treating as zero those singular values which are
        less than rcond times the largest singular value.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in] nrhs
    nrhs is int* \n
    The number of right hand sides, i.e., the number of columns
    of the matrices b and X. nrhs >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the m-by-n matrix. \n
    On exit, the first min(m,n) rows of A are overwritten with
    its right singular vectors, stored rowwise.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[in,out] b
    b is float/double array, dimension (ldb,nrhs) \n
    On entry, the m-by-nrhs  right hand side matrix b. \n
    On exit, b is overwritten by the n-by-nrhs solution
    matrix X.  If m >= n and rank = n, the residual
    sum-of-squares for the solution in the i-th column is given
    by the sum of squares of elements n+1:m in that column.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the matrix b, ldb >= max(1,max(m,n)).

    * @param[out] s
    s is float/double array, dimension (min(m,n)) \n
    The singular values of A in decreasing order.
    The condition number of A in the 2-norm = S(1)/S(min(m,n)).

    * @param[in] rcond
    rcond is float/double* \n
    rcond is used to determine the effective rank of A.
    Singular values S(i) <= rcond*S(1) are treated as zero. \n
    If rcond < 0, machine precision is used instead.

    * @param[out] rank
    rank is int* \n
    The effective rank of A, i.e., the number of singular values
    which are greater than rcond*S(1).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= 1, and also:
    lwork >= 3*min(m,n) + max( 2*min(m,n), max(m,n), nrhs )
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  the algorithm for computing the SVD failed to converge;
    if info = i, i off-diagonal elements of an intermediate
    bidiagonal form did not converge to zero.
    *  */
template< typename T >
int gelss( int* m, int* n, int* nrhs, T* a, int* lda, T* b, int* ldb, T*  s, T* rcond, int* rank, T* work, int* lwork, int* info )
{
  return gelss( m, n, nrhs, a, lda, b,ldb, s, rcond, rank, work, lwork, info );
}

/*! @brief The minimum-norm solution to a real linear least squares problem
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of the minimum-norm solution to a real linear least squares problem:
            minimize 2-norm(| B - A*X |)
        using the singular value decomposition (SVD) of A. A is an m-by-n matrix which may be
        rank-deficient.

        Several right hand side vectors B and solution vectors X can be handled in a single call;
        they are stored as the columns of the m-by-nrhs  right hand side matrix b and the n-by-nrhs
        solution matrix X.

        The effective rank of A is determined by treating as zero those singular values which are
        less than rcond times the largest singular value.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in] nrhs
    nrhs is int* \n
    The number of right hand sides, i.e., the number of columns
    of the matrices b and X. nrhs >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix. \n
    On exit, the first min(m,n) rows of A are overwritten with
    its right singular vectors, stored rowwise.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,m)

    * @param[in,out] b
    b is scomplex/dcomplex array, dimension (ldb,nrhs) \n
    On entry, the m-by-nrhs  right hand side matrix b. \n
    On exit, b is overwritten by the n-by-nrhs solution
    matrix X.  If m >= n and rank = n, the residual
    sum-of-squares for the solution in the i-th column is given
    by the sum of squares of elements n+1:m in that column.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the matrix b, ldb >= max(1,max(m,n)).

    * @param[out] s
    s is float/double array, dimension (min(m,n)) \n
    The singular values of A in decreasing order.
    The condition number of A in the 2-norm = S(1)/S(min(m,n)).

    * @param[in] rcond
    rcond is float/double \n
    rcond is used to determine the effective rank of A.
    Singular values S(i) <= rcond*S(1) are treated as zero. \n
    If rcond < 0, machine precision is used instead.

    * @param[out] rank
    rank is int* \n
    The effective rank of A, i.e., the number of singular values
    which are greater than rcond*S(1).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= 1, and also:
    lwork >= 3*min(m,n) + max( 2*min(m,n), max(m,n), nrhs )
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] rwork
    rwork is float/double array, dimension (5*min(m,n))

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  the algorithm for computing the SVD failed to converge;
        if info = i, i off-diagonal elements of an intermediate
        bidiagonal form did not converge to zero.
    *  */
template< typename Ta, typename Tb >
int gelss( int* m, int* n, int* nrhs, Ta* a, int* lda, Ta* b, int* ldb, Tb*  s, Tb*  rcond, int* rank, Ta* work, int* lwork, Tb* rwork, int* info )
{
  return gelss( m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info );
}

/*! @brief Product UUH or LHL, where U and L are upper or lower triangular matrices (blocked algorithm).
    *

    * @details
    * \b Purpose:
    * \verbatim
        Product UUH or LHL, where U and L are upper or lower triangular matrices (blocked algorithm).
        Computation of the product U * U**T or L**T * L, where the triangular factor U or L is
        stored in the upper or lower triangular part of the array a.

        If uplo = 'U' or 'u' then the upper triangle of the result is stored, overwriting the factor U in A.
        If uplo = 'L' or 'l' then the lower triangle of the result is stored, overwriting the factor L in A.

        This is the blocked form of the algorithm, calling Level 3 BLAS.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    Specifies whether the triangular factor stored in the array a
    is upper or lower triangular: \n
    = 'U':  Upper triangular \n
    = 'L':  Lower triangular

    * @param[in] n
    n is int* \n
    The order of the triangular factor U or L.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the triangular factor U or L. \n
    On exit, if uplo = 'U', the upper triangle of A is
    overwritten with the upper triangle of the product U * U**T; \n
    if uplo = 'L', the lower triangle of A is overwritten with
    the lower triangle of the product L**T * L.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,n)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0: if info = -k, the k-th argument had an illegal value
    *  */
template< typename T >
int lauum( char* uplo, int* n, T* a, int* lda, int* info )
{
  return lauum( uplo, n, a, lda, info );
}

/*! @brief Product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked algorithm).
    *

    * @details
    * \b Purpose:
    * \verbatim
        Product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked algorithm).
        Computation of the product U * U**T or L**T * L, where the triangular factor U or L is
        stored in the upper or lower triangular part of the array a.

        If uplo = 'U' or 'u' then the upper triangle of the result is stored, overwriting the factor U in A.
        If uplo = 'L' or 'l' then the lower triangle of the result is stored, overwriting the factor L in A.

        This is the unblocked form of the algorithm, calling Level 2 BLAS.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    Specifies whether the triangular factor stored in the array a
    is upper or lower triangular: \n
    = 'U':  Upper triangular \n
    = 'L':  Lower triangular

    * @param[in] n
    n is int* \n
    The order of the triangular factor U or L.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the triangular factor U or L. \n
    On exit, if uplo = 'U', the upper triangle of A is
    overwritten with the upper triangle of the product U * U**T; \n
    if uplo = 'L', the lower triangle of A is overwritten with
    the lower triangle of the product L**T * L.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,n)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0: if info = -k, the k-th argument had an illegal value
    *  */
template< typename T >
int lauu2( char* uplo, int* n, T* a, int* lda, int* info )
{
  return lauu2( uplo, n, a, lda, info );
}

/*! @brief Inverse of a real symmetric positive definite matrix.
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of the inverse of a real symmetric positive definite matrix a using the
        Cholesky factorization A = U**T*U or A = L*L**T computed by SPOTRF.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int*  \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] buff_A
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the triangular factor U or L from the Cholesky
    factorization A = U**T*U or A = L*L**T, as computed by
    SPOTRF. \n
    On exit, the upper or lower triangle of the (symmetric)
    inverse of A, overwriting the input factor U or L.

    * @param[in] ldim_A
    ldim_A is int* \n
    The leading dimension of the matrix a, lda >= max(1,n)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  if info = i, the (i,i) element of the factor U or L is
    zero, and the inverse could not be computed.
    *  */
template< typename T >
int potri( char* uplo, int* n, T* buff_A, int*  ldim_A, int* info )
{
  return potri( uplo, n, buff_A, ldim_A, info );
}

/*! @brief Inverse of a real upper or lower triangular matrix.
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of inverse of a real upper or lower triangular matrix

        This is the Level 3 BLAS version of the algorithm.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U':  A is upper triangular; \n
    = 'L':  A is lower triangular.

    * @param[in] diag
    diag is char* \n
    = 'N':  A is non-unit triangular; \n
    = 'U':  A is unit triangular.

    * @param[in] n
    n is int*  \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the triangular matrix a.  If uplo = 'U', the
    leading n-by-n upper triangular part of the array a contains
    the upper triangular matrix, and the strictly lower
    triangular part of a is not referenced.  If uplo = 'L', the
    leading n-by-n lower triangular part of the array a contains
    the lower triangular matrix, and the strictly upper
    triangular part of a is not referenced.  If diag = 'U', the
    diagonal elements of a are also not referenced and are
    assumed to be 1. \n
    On exit, the (triangular) inverse of the original matrix, in
    the same storage format.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,n)

    * @param[out] info
    info is int* \n
    = 0: successful exit \n
    < 0: if info = -i, the i-th argument had an illegal value \n
    \> 0: if info = i, A(i,i) is exactly zero.  The triangular
    matrix is singular and its inverse can not be computed.
    *  */
template< typename T >
int trtri( char* uplo, char* diag, int* n, T* a, int* lda, int* info )
{
  return trtri( uplo, diag, n, a, lda, info );
}


/*! @brief Inverse of a triangular matrix (unblocked algorithm).
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of inverse of a real upper or lower triangular matrix (unblocked algorithm)

        This is the Level 2 BLAS version of the algorithm.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U':  A is upper triangular; \n
    = 'L':  A is lower triangular.

    * @param[in] diag
    diag is char* \n
    = 'N':  A is non-unit triangular; \n
    = 'U':  A is unit triangular.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the triangular matrix a.  If uplo = 'U', the
    leading n by n upper triangular part of the array a contains
    the upper triangular matrix, and the strictly lower
    triangular part of a is not referenced.  If uplo = 'L', the
    leading n by n lower triangular part of the array a contains
    the lower triangular matrix, and the strictly upper
    triangular part of a is not referenced.  If diag = 'U', the
    diagonal elements of a are also not referenced and are
    assumed to be 1. \n
    On exit, the (triangular) inverse of the original matrix, in
    the same storage format.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the matrix a, lda >= max(1,n)

    * @param[out] info
    info is int* \n
    = 0: successful exit \n
    < 0: if info = -k, the k-th argument had an illegal value
    *  */
template< typename T >
int trti2( char* uplo, char* diag, int* n, T* a, int* lda, int* info )
{
  return trti2( uplo, diag, n, a, lda, info );
}

/*! @brief Solving Sylvester matrix equation
    *

    * @details
    * \b Purpose:
    * \verbatim
        Solution for real Sylvester matrix equation:
            op(A)*X + X*op(B) = scale*C or
            op(A)*X - X*op(B) = scale*C,

        where op(A) = A or A**T, and  a and b are both upper quasi- triangular.
        A is M-by-M and B is n-by-n; the right hand side C and the solution X are m-by-n;
        and scale is an output scale factor, set <= 1 to avoid overflow in X.

        a and b must be in Schur canonical form (as returned by SHSEQR), that is, block upper
        triangular with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has its
        diagonal elements equal and its off-diagonal elements of opposite sign.
    \endverbatim

    * @param[in] transa
    transa is char* \n
    Specifies the option op(A): \n
    = 'N': op(A) = A (No transpose) \n
    = 'T': op(A) = A**T (Transpose) \n
    = 'C': op(A) = A**H (Conjugate transpose = Transpose)

    * @param[in] transb
    transb is char* \n
    Specifies the option op(B): \n
    = 'N': op(B) = B (No transpose) \n
    = 'T': op(B) = B**T (Transpose) \n
    = 'C': op(B) = B**H (Conjugate transpose = Transpose)

    * @param[in] isgn
    isgn is int* \n
    Specifies the sign in the equation: \n
    = +1: solve op(A)*X + X*op(B) = scale*C \n
    = -1: solve op(A)*X - X*op(B) = scale*C

    * @param[in] m
    m is int* \n
    The order of the matrix a, and the number of rows in the
    matrices X and C. m >= 0.

    * @param[in] n
    n is int* \n
    The order of the matrix b, and the number of columns in the
    matrices X and C. n >= 0.

    * @param[in] a
    a is float/double array, dimension (lda,m) \n
    The upper quasi-triangular matrix a, in Schur canonical form.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,m).

    * @param[in] b
    b is float/double array, dimension (ldb,n) \n
    The upper quasi-triangular matrix b, in Schur canonical form.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the array b. ldb >= max(1,n).

    * @param[in,out] c
    c is float/double array, dimension (ldc,n) \n
    On entry, the m-by-n right hand side matrix c.
    On exit, c is overwritten by the solution matrix X.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m)

    * @param[out] scale
    scale is float/double* \n
    The scale factor, scale, set <= 1 to avoid overflow in X

    * @param[out] info
    info is int* \n
    = 0: successful exit \n
    < 0: if info = -i, the i-th argument had an illegal value \n
    = 1: a and b have common or very close eigenvalues; perturbed
    values were used to solve the equation (but the matrices
    a and b are unchanged).
    *  */
template< typename T >
int trsyl( char* transa, char* transb, int* isgn, int* m, int* n, T* a, int* lda, T* b, int* ldb, T* c, int* ldc, T* scale, int* info )
{
  return trsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
}

/*! @brief Solving Sylvester matrix equation

* @details
    * \b Purpose:
    * \verbatim
        Solution for complex Sylvester matrix equation:
            op(A)*X + X*op(B) = scale*C or
            op(A)*X - X*op(B) = scale*C,

        where op(A) = A or A**H, and  a and b are both upper quasi- triangular.
        A is M-by-M and B is n-by-n; the right hand side C and the solution X are m-by-n;
        and scale is an output scale factor, set <= 1 to avoid overflow in X.

        a and b must be in Schur canonical form (as returned by SHSEQR), that is, block upper
        triangular with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has its
        diagonal elements equal and its off-diagonal elements of opposite sign.
    \endverbatim

    * @param[in] transa
    transa is char* \n
    Specifies the option op(A): \n
    = 'N': op(A) = A    (No transpose) \n
    = 'C': op(A) = A**H (Conjugate transpose)

    * @param[in] transb
    transb is char* \n
    Specifies the option op(B): \n
    = 'N': op(B) = B    (No transpose) \n
    = 'C': op(B) = B**H (Conjugate transpose)

    * @param[in] isgn
    isgn is int* \n
    Specifies the sign in the equation: \n
    = +1: solve op(A)*X + X*op(B) = scale*C \n
    = -1: solve op(A)*X - X*op(B) = scale*C

    * @param[in] m
    m is int* \n
    The order of the matrix a, and the number of rows in the
    matrices X and C. m >= 0.

    * @param[in] n
    n is int* \n
    The order of the matrix b, and the number of columns in the
    matrices X and C. n >= 0.

    * @param[in] a
    a is scomplex/dcomplex array, dimension (lda,m) \n
    The upper quasi-triangular matrix a, in Schur canonical form.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,m).

    * @param[in] b
    b is scomplex/dcomplex array, dimension (ldb,n) \n
    The upper quasi-triangular matrix b, in Schur canonical form.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the array b. ldb >= max(1,n).

    * @param[in,out] c
    c is scomplex/dcomplex array, dimension (ldc,n) \n
    On entry, the m-by-n right hand side matrix c. \n
    On exit, c is overwritten by the solution matrix X.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m)

    * @param[out] scale
    scale is float\double* \n
    The scale factor, scale, set <= 1 to avoid overflow in X.

    * @param[out] info
    info is int* \n
    = 0: successful exit \n
    < 0: if info = -i, the i-th argument had an illegal value \n
    = 1: a and b have common or very close eigenvalues; perturbed
        values were used to solve the equation (but the matrices
        a and b are unchanged).
    *  */
template< typename Ta, typename Tb >
int trsyl( char* transa, char* transb, int* isgn, int* m, int* n, Ta* a, int* lda, Ta* b, int* ldb, Ta* c, int* ldc, Tb* scale, int* info )
{
  return trsyl( transa, transb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info );
}

/*! @brief Reduction to upper Hessenberg form
    *

    * @details
    * \b Purpose:
    * \verbatim
    Reduction of a real general matrix a to upper Hessenberg form H by an orthogonal
    similarity transformation:
    Q**T * A * Q = H .
    \endverbatim

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in] ilo
    ilo is int* \n

    * @param[in] ihi
    ihi is int* \n
    It is assumed that A is already upper triangular in rows
    and columns 1:ilo-1 and ihi+1:n. \n
    ilo and ihi are normally set by a previous call to SGEBAL;
    otherwise they should be set to 1 and N respectively. See Further Details. \n
    1 <= ilo <= ihi <= N, if N > 0; ilo=1 and ihi=0, if N=0.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the n-by-n general matrix to be reduced. \n
    On exit, the upper triangle and the first subdiagonal of A
    are overwritten with the upper Hessenberg matrix H, and the
    elements below the first subdiagonal, with the array tau,
    represent the orthogonal matrix Q as a product of elementary
    reflectors. See Further Details.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] tau
    tau is float/double/scomplex/dcomplex array, dimension (n-1) \n
    The scalar factors of the elementary reflectors (see Further
    Details). Elements 1:ilo-1 and ihi:N-1 of tau are set to
    zero.

    * @param[out] work
    work is float/double/scomplex/dcomplex array, dimension (lwork) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The length of the array work.  lwork >= max(1,n). \n
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value.
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of (ihi-ilo) elementary reflectors

            Q = H(ilo) H(ilo+1) . . . H(ihi-1).

            Each H(i) has the form

            H(i) = I - tau * V * V**T

            where tau is a real scalar, and V is a real vector with V(1:i) = 0, V(i+1) = 1 and
            V(ihi+1:n) = 0; V(i+2:ihi) is stored on exit in A(i+2:ihi,i), and tau in tau(i).

            The contents of A are illustrated by the following example, with n = 7, ilo = 2 and ihi = 6:
            on entry,   on exit,

            ( a   a   a   a   a   a   a ) (  a   a   h   h   h   h   a )
            (  a   a   a   a   a   a )    (   a   h   h   h   h   a )
            (  a   a   a   a   a   a )    (   h   h   h   h   h   h )
            (  a   a   a   a   a   a )    (   v2  h   h   h   h   h )
            (  a   a   a   a   a   a )    (   v2  v3  h   h   h   h )
            (  a   a   a   a   a   a )    (   v2  v3  v4  h   h   h )
            ( a ) (  a )
            where,
            a denotes an element of the original matrix a,
            h denotes a modified element of the upper Hessenberg matrix H,
            vi denotes an element of the vector defining H(i).
    \endverbatim
    *  */
template< typename T >
int gehrd( int* n, int* ilo, int* ihi, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return gehrd( n, ilo, ihi, a, lda, tau, work, lwork, info );
}

/*! @brief Reduction to upper Hessenberg form using an unblocked algorithm
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real general matrix a to upper Hessenberg form H by an orthogonal
        similarity transformation:
        Q**T * A * Q = H .
    \endverbatim

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in] ilo
    ilo is int*

    * @param[in] ihi
    ihi is int* \n
    It is assumed that A is already upper triangular in rows
    and columns 1:ilo-1 and ihi+1:n. \n
    ilo and ihi are normally set by a previous call to SGEBAL;
    otherwise they should be set to 1 and N respectively. See Further Details. \n
    1 <= ilo <= ihi <= max(1,n).

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the n-by-n general matrix to be reduced. \n
    On exit, the upper triangle and the first subdiagonal of A
    are overwritten with the upper Hessenberg matrix H, and the
    elements below the first subdiagonal, with the array tau,
    represent the orthogonal matrix Q as a product of elementary
    reflectors. See Further Details.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] tau
    tau is float/double/scomplex/dcomplex array, dimension (n-1) \n
    The scalar factors of the elementary reflectors (see Further
    Details).

    * @param[out] work
    work is float/double/scomplex/dcomplex array, dimension (n)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value.
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrix Q is represented as a product of (ihi-ilo) elementary reflectors

            Q = H(ilo) H(ilo+1) . . . H(ihi-1).

            Each H(i) has the form

            H(i) = I - tau * V * V**T

            where tau is a real scalar, and V is a real vector with V(1:i) = 0, V(i+1) = 1 and
            V(ihi+1:n) = 0; V(i+2:ihi) is stored on exit in A(i+2:ihi,i), and tau in tau(i).

            The contents of A are illustrated by the following example, with n = 7, ilo = 2 and ihi = 6:
            on entry,   on exit,

            ( a   a   a   a   a   a   a ) (  a   a   h   h   h   h   a )
            (  a   a   a   a   a   a ) (   a   h   h   h   h   a )
            (  a   a   a   a   a   a ) (   h   h   h   h   h   h )
            (  a   a   a   a   a   a ) (   v2  h   h   h   h   h )
            (  a   a   a   a   a   a ) (   v2  v3  h   h   h   h )
            (  a   a   a   a   a   a ) (   v2  v3  v4  h   h   h )
            ( a ) (  a )
            where,
            a denotes an element of the original matrix a,
            h denotes a modified element of the upper Hessenberg matrix H,
            vi denotes an element of the vector defining H(i).
    \endverbatim
    *  */
template< typename T >
int gehd2( int* n, int* ilo, int* ihi, T* a, int* lda, T* tau, T* work, int* info )
{
  return gehd2( n, ilo, ihi, a, lda, tau, work, info );
}


/*! @brief Reduction of a real symmetric matrix a to real symmetric tridiagonal form
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric matrix a to real symmetric tridiagonal form T by an
        orthogonal similarity transformation:
            Q**T * A * Q = T.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the symmetric matrix a. \n
    If uplo = 'U', the leading
    n-by-n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced. \n
    If uplo = 'L', the leading n-by-n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced. \n
    On exit, if uplo = 'U', the diagonal and first superdiagonal
    of A are overwritten by the corresponding elements of the
    tridiagonal matrix T, and the elements above the first
    superdiagonal, with the array tau, represent the orthogonal
    matrix Q as a product of elementary reflectors; \n
    if uplo = 'L', the diagonal and first subdiagonal of A are over-
    written by the corresponding elements of the tridiagonal
    matrix T, and the elements below the first subdiagonal, with
    the array tau, represent the orthogonal matrix Q as a product
    of elementary reflectors. See Further Details.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] d
    d is float/double array, dimension (n) \n
    The diagonal elements of the tridiagonal matrix T: \n
    D(i) = A(i,i).

    * @param[out] e
    e is float/double array, dimension (n-1) \n
    The off-diagonal elements of the tridiagonal matrix T: \n
    E(i) = A(i,i+1) if uplo = 'U', E(i) = A(i+1,i) if uplo = 'L'.

    * @param[out] tau
    tau is float/double array, dimension (n-1) \n
    The scalar factors of the elementary reflectors (see Further
    Details).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work.  lwork >= 1. \n
    For optimum performance lwork >= n*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            If uplo = 'U', the matrix Q is represented as a product of elementary reflectors

            Q = H(n-1) . . . H(2) H(1).

            Each H(i) has the form

            H(i) = I - tau * V * V**T

            where tau is a real scalar, and V is a real vector with V(i+1:n) = 0 and V(i) = 1;
            V(1:i-1) is stored on exit in A(1:i-1,i+1), and tau in tau(i).

            If uplo = 'L', the matrix Q is represented as a product of elementary reflectors

            Q = H(1) H(2) . . . H(n-1).

            Each H(i) has the form

            H(i) = I - tau * V * V**T

            where tau is a real scalar, and V is a real vector with V(1:i) = 0 and V(i+1) = 1;
            V(i+2:n) is stored on exit in A(i+2:n,i), and tau in tau(i).

            The contents of A on exit are illustrated by the following examples with n = 5:

            if uplo = 'U':  if uplo = 'L':

            (  d   e   v2  v3  v4 )  (  d   )
            (   d   e   v3  v4 )  (  e   d  )
            ( d   e   v4 )  (  v1  e   d )
            (  d   e  )  (  v1  v2  e   d   )
            (   d  )  (  v1  v2  v3  e   d  )

            where d and e denote diagonal and off-diagonal elements of T, and vi denotes an element
            of the vector defining H(i).
    \endverbatim
    *  */
template< typename T >
int sytrd( char* uplo, int* n, T* a, int* lda, T*  d, T*  e, T* tau, T* work, int* lwork, int* info )
{
  return sytrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
}

/*! @brief Reduction of a complex Hermitian matrix a to real symmetric tridiagonal form
    *

    * @details
    * \b Purpose:
    * \verbatim
    Reduction of a complex Hermitian matrix a to real symmetric tridiagonal form T by a
    unitary similarity transformation:
        Q**H * A * Q = T.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the symmetric matrix a. \n
    If uplo = 'U', the leading
    n-by-n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced. \n
    If uplo = 'L', the leading n-by-n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced. \n
    On exit, if uplo = 'U', the diagonal and first superdiagonal
    of A are overwritten by the corresponding elements of the
    tridiagonal matrix T, and the elements above the first
    superdiagonal, with the array tau, represent the orthogonal
    matrix Q as a product of elementary reflectors; \n
    if uplo = 'L', the diagonal and first subdiagonal of A are over-
    written by the corresponding elements of the tridiagonal
    matrix T, and the elements below the first subdiagonal, with
    the array tau, represent the orthogonal matrix Q as a product
    of elementary reflectors. See Further Details.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] d
    d is float/double array, dimension (n) \n
    The diagonal elements of the tridiagonal matrix T: \n
    D(i) = A(i,i).

    * @param[out] e
    e is float/double array, dimension (n-1) \n
    The off-diagonal elements of the tridiagonal matrix T: \n
    E(i) = A(i,i+1) if uplo = 'U', E(i) = A(i+1,i) if uplo = 'L'.

    * @param[out] tau
    tau is scomplex/dcomplex array, dimension (n-1) \n
    The scalar factors of the elementary reflectors (see Further
    Details).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work.  lwork >= 1. \n
    For optimum performance lwork >= n*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            If uplo = 'U', the matrix Q is represented as a product of elementary reflectors

            Q = H(n-1) . . . H(2) H(1).

            Each H(i) has the form

            H(i) = I - tau * V * V**H

            where tau is a complex scalar, and V is a complex vector with V(i+1:n) = 0 and V(i) = 1;
            V(1:i-1) is stored on exit in A(1:i-1,i+1), and tau in tau(i).

            If uplo = 'L', the matrix Q is represented as a product of elementary reflectors

            Q = H(1) H(2) . . . H(n-1).

            Each H(i) has the form

            H(i) = I - tau * V * V**H

            where tau is a complex scalar, and V is a complex vector with V(1:i) = 0 and V(i+1) = 1;
            V(i+2:n) is stored on exit in A(i+2:n,i), and tau in tau(i).

            The contents of A on exit are illustrated by the following examples with n = 5:

            if uplo = 'U':  if uplo = 'L':

            (  d   e   v2  v3  v4 )  (  d   )
            (   d   e   v3  v4 )  (  e   d  )
            ( d   e   v4 )  (  v1  e   d )
            (  d   e  )  (  v1  v2  e   d   )
            (   d  )  (  v1  v2  v3  e   d  )

            where d and e denote diagonal and off-diagonal elements of T, and vi denotes an element
            of the vector defining H(i).
    \endverbatim
    *  */
template< typename Ta, typename Tb >
int hetrd( char* uplo, int* n, Ta* a, int* lda, Tb*  d, Tb*  e, Ta* tau, Ta* work, int* lwork, int* info )
{
  return hetrd( uplo, n, a, lda, d, e, tau, work, lwork, info );
}

/*! @brief Reduction of a real symmetric matrix a to real symmetric tridiagonal form (unblocked algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric matrix a to real symmetric tridiagonal form T by an orthogonal
        similarity transformation(unblocked algorithm):
            Q**T * A * Q = T.
    \endverbatim

    * @param[in] uplo
    uplo is char*
    = 'U':  Upper triangle of A is stored;
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int*
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n)
    On entry, the symmetric matrix a.
    If uplo = 'U', the leading
    n-by-n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced.
    If uplo = 'L', the leading n-by-n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced.
    On exit, if uplo = 'U', the diagonal and first superdiagonal
    of A are overwritten by the corresponding elements of the
    tridiagonal matrix T, and the elements above the first
    superdiagonal, with the array tau, represent the orthogonal
    matrix Q as a product of elementary reflectors; \n
    if uplo = 'L', the diagonal and first subdiagonal of A are over-
    written by the corresponding elements of the tridiagonal
    matrix T, and the elements below the first subdiagonal, with
    the array tau, represent the orthogonal matrix Q as a product
    of elementary reflectors. See Further Details.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] d
    d is float/double array, dimension (n) \n
    The diagonal elements of the tridiagonal matrix T: \n
    D(i) = A(i,i).

    * @param[out] e
    e is float/double array, dimension (n-1) \n
    The off-diagonal elements of the tridiagonal matrix T: \n
    E(i) = A(i,i+1) if uplo = 'U', E(i) = A(i+1,i) if uplo = 'L'.

    * @param[out] tau
    tau is float/double array, dimension (n-1) \n
    The scalar factors of the elementary reflectors (see Further
    Details).

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            If uplo = 'U', the matrix Q is represented as a product of elementary reflectors

            Q = H(n-1) . . . H(2) H(1).

            Each H(i) has the form

            H(i) = I - tau * V * V**T

            where tau is a real scalar, and V is a real vector with V(i+1:n) = 0 and V(i) = 1;
            V(1:i-1) is stored on exit in A(1:i-1,i+1), and tau in tau(i).

            If uplo = 'L', the matrix Q is represented as a product of elementary reflectors

            Q = H(1) H(2) . . . H(n-1).

            Each H(i) has the form

            H(i) = I - tau * V * V**T

            where tau is a real scalar, and V is a real vector with V(1:i) = 0 and V(i+1) = 1;
            V(i+2:n) is stored on exit in A(i+2:n,i), and tau in tau(i).

            The contents of A on exit are illustrated by the following examples with n = 5:

            if uplo = 'U':  if uplo = 'L':

            (  d   e   v2  v3  v4 )  (  d   )
            (   d   e   v3  v4 )  (  e   d  )
            ( d   e   v4 )  (  v1  e   d )
            (  d   e  )  (  v1  v2  e   d   )
            (   d  )  (  v1  v2  v3  e   d  )

            where d and e denote diagonal and off-diagonal elements of T, and vi denotes an element
            of the vector defining H(i).
    \endverbatim
    *  */
template< typename T >
int sytd2( char* uplo, int* n, T* a, int* lda, T*  d, T*  e, T* tau, int* info )
{
  return sytd2( uplo, n, a, lda, d, e, tau, info );
}


/*! @brief Reduction of a Hermitian matrix a to real symmetric tridiagonal form (unblocked algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a Hermitian matrix a to real symmetric tridiagonal form T by an unitary
        similarity transformation(unblocked algorithm):
            Q**T * A * Q = T.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the symmetric matrix a. \n
    If uplo = 'U', the leading
    n-by-n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced. \n
    If uplo = 'L', the leading n-by-n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced. \n
    On exit, if uplo = 'U', the diagonal and first superdiagonal
    of A are overwritten by the corresponding elements of the
    tridiagonal matrix T, and the elements above the first
    superdiagonal, with the array tau, represent the orthogonal
    matrix Q as a product of elementary reflectors; \n
    if uplo = 'L', the diagonal and first subdiagonal of A are over-
    written by the corresponding elements of the tridiagonal
    matrix T, and the elements below the first subdiagonal, with
    the array tau, represent the orthogonal matrix Q as a product
    of elementary reflectors. See Further Details.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] d
    d is float/double array, dimension (n) \n
    The diagonal elements of the tridiagonal matrix T: \n
    D(i) = A(i,i).

    * @param[out] e
    e is float/double array, dimension (n-1) \n
    The off-diagonal elements of the tridiagonal matrix T: \n
    E(i) = A(i,i+1) if uplo = 'U', E(i) = A(i+1,i) if uplo = 'L'.

    * @param[out] tau
    tau is scomplex/dcomplex array, dimension (n-1) \n
    The scalar factors of the elementary reflectors (see Further
    Details).

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *
    * \n
    * **Further Details**
    * \verbatim
            If uplo = 'U', the matrix Q is represented as a product of elementary reflectors

                Q = H(n-1) . . . H(2) H(1).

            Each H(i) has the form

                H(i) = I - tau * V * V**H

            where tau is a complex scalar, and V is a complex vector with V(i+1:n) = 0 and V(i) = 1;
            V(1:i-1) is stored on exit in A(1:i-1,i+1), and tau in tau(i).

            If uplo = 'L', the matrix Q is represented as a product of elementary reflectors

                Q = H(1) H(2) . . . H(n-1).

            Each H(i) has the form

                H(i) = I - tau * V * V**H

            where tau is a complex scalar, and V is a complex complex vector with V(1:i) = 0 and V(i+1) = 1;
            V(i+2:n) is stored on exit in A(i+2:n,i), and tau in tau(i).

            The contents of A on exit are illustrated by the following examples with n = 5:

            if uplo = 'U':                       if uplo = 'L':

                (  d   e   v2  v3  v4 )              (  d                  )
                (      d   e   v3  v4 )              (  e   d              )
                (          d   e   v4 )              (  v1  e   d          )
                (              d   e  )              (  v1  v2  e   d      )
                (                  d  )              (  v1  v2  v3  e   d  )

            where d and e denote diagonal and off-diagonal elements of T, and vi denotes an element
            of the vector defining H(i).
    \endverbatim
    *  */
template< typename Ta, typename Tb >
int hetd2( char* uplo, int* n, Ta* a, int* lda, Tb*  d, Tb*  e, Ta* tau, int* info )
{
  return hetd2( uplo, n, a, lda, d, e, tau, info );
}

/*! @brief Reduction to bidiagonal form
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a general real m-by-n matrix a to upper or lower bidiagonal form B by an
        orthogonal transformation: Q**T * A * P = B.

        If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows in the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns in the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the m-by-n general matrix to be reduced. \n
    On exit, \n
    if m >= n, the diagonal and the first superdiagonal are
    overwritten with the upper bidiagonal matrix b; the
    elements below the diagonal, with the array tauq, represent
    the orthogonal matrix Q as a product of elementary
    reflectors, and the elements above the first superdiagonal,
    with the array taup, represent the orthogonal matrix P as
    a product of elementary reflectors; \n
    if m < n, the diagonal and the first subdiagonal are
    overwritten with the lower bidiagonal matrix b; the
    elements below the first subdiagonal, with the array tauq,
    represent the orthogonal matrix Q as a product of
    elementary reflectors, and the elements above the diagonal,
    with the array taup, represent the orthogonal matrix P as
    a product of elementary reflectors. \n
    See Further Details.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,m).

    * @param[out] d
    d is float/double array, dimension (min(m,n)) \n
    The diagonal elements of the bidiagonal matrix b: \n
    D(i) = A(i,i).

    * @param[out] e
    e is float/double array, dimension (min(m,n)-1) \n
    The off-diagonal elements of the bidiagonal matrix b: \n
    if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; \n
    if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.

    * @param[out] tauq
    tauq is float/double array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors which
    represent the orthogonal matrix Q. \n
    See Further Details.

    * @param[out] taup
    taup is float/double array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors which
    represent the orthogonal matrix P. \n
    See Further Details.

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The length of the array work.  lwork >= max(1,m,n). \n
    For optimum performance lwork >= (M+N)*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value.
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrices Q and P are represented as products of elementary reflectors:

            If m >= n,

            Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)

            Each H(i) and G(i) has the form:

            H(i) = I - tauq * V * V**T  and G(i) = I - taup * u * u**T

            where tauq and taup are real scalars, and V and u are real vectors;
            V(1:i-1) = 0, V(i) = 1, and V(i+1:m) is stored on exit in A(i+1:m,i);
            u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
            tauq is stored in tauq(i) and taup in taup(i).

            If m < n,

            Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)

            Each H(i) and G(i) has the form:

            H(i) = I - tauq * V * V**T  and G(i) = I - taup * u * u**T

            where tauq and taup are real scalars, and V and u are real vectors;
            V(1:i) = 0, V(i+1) = 1, and V(i+2:m) is stored on exit in A(i+2:m,i);
            u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
            tauq is stored in tauq(i) and taup in taup(i).

            The contents of A on exit are illustrated by the following examples:

            m = 6 and n = 5 (m > n): m = 5 and n = 6 (m < n):

            (  d   e   u1  u1  u1 )  (  d   u1  u1  u1  u1  u1 )
            (  v1  d   e   u2  u2 )  (  e   d   u2  u2  u2  u2 )
            (  v1  v2  d   e   u3 )  (  v1  e   d   u3  u3  u3 )
            (  v1  v2  v3  d   e  )  (  v1  v2  e   d   u4  u4 )
            (  v1  v2  v3  v4  d  )  (  v1  v2  v3  e   d   u5 )
            (  v1  v2  v3  v4  v5 )

            where d and e denote diagonal and off-diagonal elements of B, vi denotes an element
            of the vector defining H(i), and ui an element of the vector defining G(i).
    \endverbatim
    *  */
template< typename T >
int gebrd( int* m, int* n, T* a, int* lda, T*  d, T*  e, T* tauq, T* taup, T* work, int* lwork, int* info )
{
  return gebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
}


/*! @brief Reduction to bidiagonal form
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a general complex m-by-n matrix a to upper or lower bidiagonal form B by an
        orthogonal transformation: Q**H * A * P = B.

        If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows in the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns in the matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n general matrix to be reduced. \n
    On exit, \n
    if m >= n, the diagonal and the first superdiagonal are
    overwritten with the upper bidiagonal matrix b; the
    elements below the diagonal, with the array tauq, represent
    the orthogonal matrix Q as a product of elementary
    reflectors, and the elements above the first superdiagonal,
    with the array taup, represent the orthogonal matrix P as
    a product of elementary reflectors; \n
    if m < n, the diagonal and the first subdiagonal are
    overwritten with the lower bidiagonal matrix b; the
    elements below the first subdiagonal, with the array tauq,
    represent the orthogonal matrix Q as a product of
    elementary reflectors, and the elements above the diagonal,
    with the array taup, represent the orthogonal matrix P as
    a product of elementary reflectors. \n
    See Further Details.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,m).

    * @param[out] d
    d is float/double array, dimension (min(m,n)) \n
    The diagonal elements of the bidiagonal matrix b: \n
    D(i) = A(i,i).

    * @param[out] e
    e is float/double array, dimension (min(m,n)-1) \n
    The off-diagonal elements of the bidiagonal matrix b: \n
    if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; \n
    if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.

    * @param[out] tauq
    tauq is scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors which
    represent the orthogonal matrix Q. \n
    See Further Details.

    * @param[out] taup
    taup is scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors which
    represent the orthogonal matrix P.
    \n See Further Details.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The length of the array work.  lwork >= max(1,m,n). \n
    For optimum performance lwork >= (M+N)*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value.
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrices Q and P are represented as products of elementary reflectors:

            If m >= n,

                Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)

            Each H(i) and G(i) has the form:

                H(i) = I - tauq * V * V**H  and G(i) = I - taup * u * u**H

            where tauq and taup are complex scalars, and V and u are complex vectors;
            V(1:i-1) = 0, V(i) = 1, and V(i+1:m) is stored on exit in A(i+1:m,i);
            u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
            tauq is stored in tauq(i) and taup in taup(i).

            If m < n,

                Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)

            Each H(i) and G(i) has the form:

                H(i) = I - tauq * V * V**H  and G(i) = I - taup * u * u**H

            where tauq and taup are complex scalars, and V and u are complex vectors;
            V(1:i) = 0, V(i+1) = 1, and V(i+2:m) is stored on exit in A(i+2:m,i);
            u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
            tauq is stored in tauq(i) and taup in taup(i).

            The contents of A on exit are illustrated by the following examples:

            m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):

                (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
                (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
                (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
                (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
                (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
                (  v1  v2  v3  v4  v5 )

            where d and e denote diagonal and off-diagonal elements of B, vi denotes an element
            of the vector defining H(i), and ui an element of the vector defining G(i).
    \endverbatim
    *  */
template< typename Ta, typename Tb >
int gebrd( int* m, int* n, Ta* a, int* lda, Tb* d, Tb* e, Ta* tauq, Ta* taup, Ta* work, int* lwork, int* info )
{
  return gebrd( m, n, a, lda, d, e, tauq, taup, work, lwork, info );
}


/*! @brief Reduction to bidiagonal form (unblocked algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a general real m-by-n matrix a to upper or lower bidiagonal form B by an
        orthogonal transformation: Q**T * A * P = B.

        If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows in the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns in the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the m-by-n general matrix to be reduced. \n
    On exit, \n
    if m >= n, the diagonal and the first superdiagonal are
    overwritten with the upper bidiagonal matrix b; the
    elements below the diagonal, with the array tauq, represent
    the orthogonal matrix Q as a product of elementary
    reflectors, and the elements above the first superdiagonal,
    with the array taup, represent the orthogonal matrix P as
    a product of elementary reflectors; \n
    if m < n, the diagonal and the first subdiagonal are
    overwritten with the lower bidiagonal matrix b; the
    elements below the first subdiagonal, with the array tauq,
    represent the orthogonal matrix Q as a product of
    elementary reflectors, and the elements above the diagonal,
    with the array taup, represent the orthogonal matrix P as
    a product of elementary reflectors. \n
    See Further Details.

    * @param[in] lda
    lda is int*
    The leading dimension of the array a.  lda >= max(1,m).

    * @param[out] d
    d is float/double array, dimension (min(m,n)) \n
    The diagonal elements of the bidiagonal matrix b: \n
    D(i) = A(i,i).

    * @param[out] e
    e is float/double array, dimension (min(m,n)-1) \n
    The off-diagonal elements of the bidiagonal matrix b: \n
    if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; \n
    if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.

    * @param[out] tauq
    tauq is float/double array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors which
    represent the orthogonal matrix Q. \n
    See Further Details.

    * @param[out] taup
    taup is float/double array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors which
    represent the orthogonal matrix P. \n
    See Further Details.

    * @param[out] work
    work is float/double array, dimension (max(m,n))

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value.
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrices Q and P are represented as products of elementary reflectors:

            If m >= n,

            Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)

            Each H(i) and G(i) has the form:

            H(i) = I - tauq * V * V**T  and G(i) = I - taup * u * u**T

            where tauq and taup are real scalars, and V and u are real vectors;
            V(1:i-1) = 0, V(i) = 1, and V(i+1:m) is stored on exit in A(i+1:m,i);
            u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
            tauq is stored in tauq(i) and taup in taup(i).

            If m < n,

            Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)

            Each H(i) and G(i) has the form:

            H(i) = I - tauq * V * V**T  and G(i) = I - taup * u * u**T

            where tauq and taup are real scalars, and V and u are real vectors;
            V(1:i) = 0, V(i+1) = 1, and V(i+2:m) is stored on exit in A(i+2:m,i);
            u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
            tauq is stored in tauq(i) and taup in taup(i).

            The contents of A on exit are illustrated by the following examples:

            m = 6 and n = 5 (m > n): m = 5 and n = 6 (m < n):

            (  d   e   u1  u1  u1 )  (  d   u1  u1  u1  u1  u1 )
            (  v1  d   e   u2  u2 )  (  e   d   u2  u2  u2  u2 )
            (  v1  v2  d   e   u3 )  (  v1  e   d   u3  u3  u3 )
            (  v1  v2  v3  d   e  )  (  v1  v2  e   d   u4  u4 )
            (  v1  v2  v3  v4  d  )  (  v1  v2  v3  e   d   u5 )
            (  v1  v2  v3  v4  v5 )

            where d and e denote diagonal and off-diagonal elements of B, vi denotes an element of
            the vector defining H(i), and ui an element of the vector defining G(i).
    \endverbatim
    *  */
template< typename T >
int gebd2( int* m, int* n, T* a, int* lda, T*  d, T*  e, T* tauq, T* taup, T* work, int* info )
{
  return gebd2( m, n, a, lda, d, e, tauq, taup, work, info );
}


/*! @brief Reduction to bidiagonal form (unblocked algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a general complex m-by-n matrix a to upper or lower bidiagonal form B by an
        orthogonal transformation: Q**H * A * P = B (unblocked algorithm).

        If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows in the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns in the matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n general matrix to be reduced. \n
    On exit, \n
    if m >= n, the diagonal and the first superdiagonal are
    overwritten with the upper bidiagonal matrix b; the
    elements below the diagonal, with the array tauq, represent
    the orthogonal matrix Q as a product of elementary
    reflectors, and the elements above the first superdiagonal,
    with the array taup, represent the orthogonal matrix P as
    a product of elementary reflectors; \n
    if m < n, the diagonal and the first subdiagonal are
    overwritten with the lower bidiagonal matrix b; the
    elements below the first subdiagonal, with the array tauq,
    represent the orthogonal matrix Q as a product of
    elementary reflectors, and the elements above the diagonal,
    with the array taup, represent the orthogonal matrix P as
    a product of elementary reflectors. \n
    See Further Details.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,m).

    * @param[out] d
    d is float/double array, dimension (min(m,n)) \n
    The diagonal elements of the bidiagonal matrix b: \n
    D(i) = A(i,i).

    * @param[out] e
    e is float/double array, dimension (min(m,n)-1) \n
    The off-diagonal elements of the bidiagonal matrix b: \n
    if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; \n
    if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.

    * @param[out] tauq
    tauq is scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors which
    represent the orthogonal matrix Q. \n
    See Further Details.

    * @param[out] taup
    taup is scomplex/dcomplex array, dimension (min(m,n)) \n
    The scalar factors of the elementary reflectors which
    represent the orthogonal matrix P. \n
    See Further Details.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(m,n))

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value.
    *
    * \n
    * **Further Details**
    * \verbatim
            The matrices Q and P are represented as products of elementary reflectors:

            If m >= n,

                Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)

            Each H(i) and G(i) has the form:

                H(i) = I - tauq * V * V**H  and G(i) = I - taup * u * u**H

            where tauq and taup are complex scalars, and V and u are complex vectors;
            V(1:i-1) = 0, V(i) = 1, and V(i+1:m) is stored on exit in A(i+1:m,i);
            u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
            tauq is stored in tauq(i) and taup in taup(i).

            If m < n,

                Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)

            Each H(i) and G(i) has the form:

                H(i) = I - tauq * V * V**H  and G(i) = I - taup * u * u**H

            where tauq and taup are complex scalars, and V and u are complex vectors;
            V(1:i) = 0, V(i+1) = 1, and V(i+2:m) is stored on exit in A(i+2:m,i);
            u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
            tauq is stored in tauq(i) and taup in taup(i).

            The contents of A on exit are illustrated by the following examples:

            m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):

                (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
                (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
                (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
                (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
                (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
                (  v1  v2  v3  v4  v5 )

            where d and e denote diagonal and off-diagonal elements of B, vi denotes an element of
            the vector defining H(i), and ui an element of the vector defining G(i).
    \endverbatim
    *  */
template< typename Ta, typename Tb >
int gebd2( int* m, int* n, Ta* a, int* lda, Tb*  d, Tb*  e, Ta* tauq, Ta* taup, Ta* work, int* info )
{
  return gebd2( m, n, a,lda, d, e, tauq, taup, work, info );
}

/*! @brief Reduction of a real symmetric-definite generalized eigenproblem to standard form
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric-definite generalized eigenproblem to standard form.

        If itype = 1, the problem is A*X = lambda*B*X,
        and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)

        If itype = 2 or 3, the problem is A*B*X = lambda*X or
        B*A*X = lambda*X, and A is overwritten by U*A*U**T or L**T*A*L.

        B must have been previously factorized as U**T*U or L*L**T by SPOTRF.
    \endverbatim

    * @param[in] itype
    itype is int* \n
    = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); \n
    = 2 or 3: compute U*A*U**T or L**T*A*L.

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored and B is factored as
    U**T*U; \n
    = 'L':  Lower triangle of A is stored and B is factored as
    L*L**T.

    * @param[in] n
    n is int* \n
    The order of the matrices A and B.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the symmetric matrix a.  If uplo = 'U', the leading
    n-by-n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced.  If uplo = 'L', the
    leading n-by-n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced. \n \n
    On exit, if info = 0, the transformed matrix, stored in the
    same format as a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[in] b
    b is float/double array, dimension (ldb,n) \n
    The triangular factor from the Cholesky factorization of B,
    as returned by SPOTRF.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the array b.  ldb >= max(1,n).

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int sygst( int* itype, char* uplo, int* n, T* a, int* lda, T* b, int* ldb, int* info )
{
  return sygst( itype, uplo, n, a, lda, b, ldb, info );
}

/*! @brief Reduction of a complex Hermitian-definite generalized eigenproblem to standard form
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric-definite generalized eigenproblem to standard form.

        If itype = 1, the problem is A*X = lambda*B*X,
        and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)

        If itype = 2 or 3, the problem is A*B*X = lambda*X or
        B*A*X = lambda*X, and A is overwritten by U*A*U**H or L**H*A*L.

        B must have been previously factorized as U**H*U or L*L**H by SPOTRF.
    \endverbatim

    * @param[in] itype
    itype is int* \n
    = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H); \n
    = 2 or 3: compute U*A*U**H or L**H*A*L.

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored and B is factored as
    U**H*U; \n
    = 'L':  Lower triangle of A is stored and B is factored as
    L*L**H.

    * @param[in] n
    n is int* \n
    The order of the matrices A and B.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the symmetric matrix a.  If uplo = 'U', the leading
    n-by-n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced.  If uplo = 'L', the
    leading n-by-n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced. \n \n
    On exit, if info = 0, the transformed matrix, stored in the
    same format as a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[in] b
    b is scomplex/dcomplex array, dimension (ldb,n) \n
    The triangular factor from the Cholesky factorization of B,
    as returned by SPOTRF.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the array b.  ldb >= max(1,n).

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int hegst( int* itype, char* uplo, int* n, T* a, int* lda, T* b, int* ldb, int* info )
{
  return hegst( itype, uplo, n, a, lda, b, ldb, info );
}

/*! @brief Reduction of a symmetric-definite generalized eigenproblem to standard form (unblocked algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a real symmetric definite generalized eigenproblem to standard form, using the
        factorization results obtained from spotrf (unblocked algorithm).
        If itype = 1, the problem is A*X = lambda*B*X,
        and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)

        If itype = 2 or 3, the problem is A*B*X = lambda*X or
        B*A*X = lambda*X, and A is overwritten by U*A*U**T or L**T *A*L.

        B must have been previously factorized as U**T *U or L*L**T by SPOTRF.
    \endverbatim

    * @param[in] itype
    itype is int* \n
    = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); \n
    = 2 or 3: compute U*A*U**T or L**T *A*L.

    * @param[in] uplo
    uplo is char* \n
    Specifies whether the upper or lower triangular part of the
    symmetric matrix a is stored, and how B has been factorized. \n
    = 'U':  Upper triangular \n
    = 'L':  Lower triangular

    * @param[in] n
    n is int* \n
    The order of the matrices A and B.  n >= 0.

    * @param[in,out] a
    a is float/doublearray, dimension (lda,n) \n
    On entry, the symmetric matrix a.  If uplo = 'U', the leading
    n by n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced.  If uplo = 'L', the
    leading n by n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced. \n \n
    On exit, if info = 0, the transformed matrix, stored in the
    same format as a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[in] b
    b is float/double array, dimension (ldb,n) \n
    The triangular factor from the Cholesky factorization of B,
    as returned by SPOTRF.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the array b.  ldb >= max(1,n).

    * @param[out] info
    info is int* \n
    = 0:  successful exit. \n
    < 0:  if info = -i, the i-th argument had an illegal value.
    *  */
template< typename T >
int sygs2( int* itype, char* uplo, int* n, T* a, int* lda, T* b, int* ldb, int* info )
{
  return sygs2( itype, uplo, n, a, lda, b, ldb, info );
}

/*! @brief Reduction of a Hermitian-definite generalized eigenproblem to standard form (unblocked algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a Hermitian definite generalized eigenproblem to standard form, using the
        factorization results obtained from cpotrf (unblocked algorithm).
        If itype = 1, the problem is A*X = lambda*B*X,
        and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)

        If itype = 2 or 3, the problem is A*B*X = lambda*X or
        B*A*X = lambda*X, and A is overwritten by U*A*U**H or L**H *A*L.

        B must have been previously factorized as U**H *U or L*L**H by CPOTRF.
    \endverbatim

    * @param[in] itype
    itype is int* \n
    = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H); \n
    = 2 or 3: compute U*A*U**H or L**H *A*L.

    * @param[in] uplo
    uplo is char* \n
    Specifies whether the upper or lower triangular part of the
    Hermitian matrix a is stored, and how B has been factorized. \n
    = 'U':  Upper triangular \n
    = 'L':  Lower triangular

    * @param[in] n
    n is int* \n
    The order of the matrices A and B.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the symmetric matrix a.  If uplo = 'U', the leading
    n by n upper triangular part of a contains the upper
    triangular part of the matrix a, and the strictly lower
    triangular part of a is not referenced.  If uplo = 'L', the
    leading n by n lower triangular part of a contains the lower
    triangular part of the matrix a, and the strictly upper
    triangular part of a is not referenced. \n \n
    On exit, if info = 0, the transformed matrix, stored in the
    same format as a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[in] b
    b is scomplex/dcomplex array, dimension (ldb,n) \n
    The triangular factor from the Cholesky factorization of B,
    as returned by CPOTRF.

    * @param[in] ldb
    ldb is int* \n
    The leading dimension of the array b.  ldb >= max(1,n).

    * @param[out] info
    info is int* \n
    = 0:  successful exit. \n
    < 0:  if info = -i, the i-th argument had an illegal value.
    *  */
template< typename T >
int hegs2( int* itype, char* uplo, int* n, T* a, int* lda, T* b, int* ldb, int* info )
{
  return hegs2(itype, uplo, n, a, lda, b, ldb, info );
}

/*! @brief Triangular factor of a block reflector
    *

    * @details
    * \b Purpose:
    * \verbatim
        Formation of the triangular factor T of a real block reflector H of order n, which is
        defined as a product of k elementary reflectors.

        If direct = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;

        If direct = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.

        If storev = 'C', the vector which defines the elementary reflector H(i) is stored in the
        i-th column of the array V, and

        H  =  I - V * T * V**T

        If storev = 'R', the vector which defines the elementary reflector H(i) is stored in the
        i-th row of the array V, and

        H  =  I - V**T * T * V
    \endverbatim

    * @param[in] direct
    direct is char* \n
    Specifies the order in which the elementary reflectors are
    multiplied to form the block reflector: \n
    = 'F': H = H(1) H(2) . . . H(k) (Forward) \n
    = 'B': H = H(k) . . . H(2) H(1) (Backward)

    * @param[in] storev
    storev is char* \n
    Specifies how the vectors which define the elementary
    reflectors are stored (see also Further Details): \n
    = 'C': columnwise \n
    = 'R': rowwise

    * @param[in] n
    n is int*
    The order of the block reflector H. n >= 0.

    * @param[in] k
    k is int* \n
    The order of the triangular factor T (= the number of
    elementary reflectors). K >= 1.

    * @param[in] v
    v is float/double/scomplex/dcomplex array, dimension \n
    (ldv,K) if storev = 'C' \n
    (ldv,N) if storev = 'R' \n
    The matrix V. See further details.

    * @param[in] ldv
    ldv is int* \n
    The leading dimension of the array V. \n
    If storev = 'C', ldv >= max(1,N); if storev = 'R', ldv >= K.

    * @param[in] tau
    tau is float/double/scomplex/dcomplex array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i).

    * @param[out] t
    t is float/double/scomplex/dcomplex array, dimension (ldt,K) \n
    The k by k triangular factor T of the block reflector. \n
    If direct = 'F', T is upper triangular; if direct = 'B', T is
    lower triangular. The rest of the array is not used.

    * @param[in] ldt
    ldt is int* \n
    The leading dimension of the array T. ldt >= K.
    *
    * \n
    * **Further Details**
    * \verbatim
            The shape of the matrix V and the storage of the vectors which define the H(i) is best
            illustrated by the following example with n = 5 and k = 3. The elements equal to 1
            are not stored.

            direct = 'F' and storev = 'C':   direct = 'F' and storev = 'R':

            V = (  1 )  V = (  1 v1 v1 v1 v1 )
            ( v1  1 )   (  1 v2 v2 v2 )
            ( v1 v2  1 )   (  1 v3 v3 )
            ( v1 v2 v3 )
            ( v1 v2 v3 )

            direct = 'B' and storev = 'C':   direct = 'B' and storev = 'R':

            V = ( v1 v2 v3 )  V = ( v1 v1  1 )
            ( v1 v2 v3 )   ( v2 v2 v2  1 )
            (  1 v2 v3 )   ( v3 v3 v3 v3  1 )
            (  1 v3 )
            (  1 )
    \endverbatim
    *  */
template< typename T >
int larft( char* direct, char* storev, int* n, int* k, T* v, int* ldv, T* tau, T* t, int* ldt )
{
  return larft( direct, storev, n, k, v, ldv, tau, t, ldt );
}

/*! @brief Generation of an elementary reflector (Householder matrix)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Generation of a real elementary reflector (Householder matrix) H of order n, such that

        H * ( alpha ) = ( beta ),   H**T * H = I.
        (   X   )   (   0  )

        where alpha and beta are scalars, and X is an (n-1)-element real vector.
        H is represented in the form

        H = I - tau * ( 1 ) * ( 1 v**T ) ,
        ( v )

        where tau is a real scalar and v is a real (n-1)-element vector.

        If the elements of X are all zero, then tau = 0 and H is taken to be the unit matrix.

        Otherwise  1 <= tau <= 2.
    \endverbatim

    * @param[in] n
    n is int* \n
    The order of the elementary reflector.

    * @param[in,out] alpha
    alpha is float/double/scomplex/dcomplex* \n
    On entry, the value alpha. \n
    On exit, it is overwritten with the value beta.

    * @param[in,out] x
    x is float/double/scomplex/dcomplex array, dimension (1+(N-2)*abs(incx)) \n
    On entry, the vector X. \n
    On exit, it is overwritten with the vector v.

    * @param[in] incx
    incx is int* \n
    The increment between elements of X. incx > 0.

    * @param[out] tau
    tau is float/double/scomplex/dcomplex* \n
    The value tau.
    *  */
template< typename T >
int larfg( int* n, T* alpha, T* x, int* incx, T* tau )
{
  return larfg( n, alpha, x, incx, tau );
}

/*! @brief Generation of an elementary reflector (Householder matrix) with non-negative beta
    *

    * @details
    * \b Purpose:
    * \verbatim
        Generation of a real elementary reflector H (Householder matrix) with non-negative beta
        of order n, such that

        H * ( alpha ) = ( beta ),   H**T * H = I.
        (   x   )   (   0  )

        where alpha and beta are scalars, beta is non-negative, and x is an (n-1)-element real vector.
        H is represented in the form

        H = I - tau * ( 1 ) * ( 1 v**T ) ,
        ( v )

        where tau is a real scalar and v is a real (n-1)-element vector.

        If the elements of x are all zero, then tau = 0 and H is taken to be the unit matrix.
    \endverbatim

    * @param[in] n
    n is int* \n
    The order of the elementary reflector.

    * @param[in,out] alpha
    alpha is float/double/scomplex/dcomplex* \n
    On entry, the value alpha. \n
    On exit, it is overwritten with the value beta.

    * @param[in,out] x
    x is float/double/scomplex/dcomplex array, dimension (1+(N-2)*abs(incx)) \n
    On entry, the vector x. \n
    On exit, it is overwritten with the vector v.

    * @param[in] incx
    incx is int* \n
    The increment between elements of X. incx > 0.

    * @param[out] tau
    tau is float/double/scomplex/dcomplex* \n
    The value tau.
    *  */
template< typename T >
int larfgp( int* n, T* alpha, T* x, int* incx, T* tau )
{
  return larfgp( n,  alpha, x, incx, tau );
}

/*! @brief Form Q from QR factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Generation of an m-by-n real matrix Q with orthonormal columns, which is defined as the
        first N columns of a product of K elementary reflectors of order M

        Q  =  H(1) H(2) . . . H(k)

        as returned by SGEQRF.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix Q. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix Q. M >= n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines the
    matrix Q. n >= k >= 0. \n

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the i-th column must contain the vector which
    defines the elementary reflector H(i), for i = 1,2,...,k, as
    returned by SGEQRF in the first k columns of its array
    argument a. \n
    On exit, the m-by-n matrix Q.

    * @param[in] lda
    lda is int* \n
    The first dimension of the array a. lda >= max(1,m).

    * @param[in] tau
    tau is float/double array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SGEQRF.

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,n). \n
    For optimum performance lwork >= n*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument has an illegal value
    *  */
template< typename T >
int orgqr( int* m, int* n, int* k, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return orgqr(  m, n, k,  a, lda, tau, work, lwork, info );
}

/*! @brief Form Q from QR factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Generation of an m-by-n complex matrix Q with orthonormal columns, which is defined as the
        first N columns of a product of K elementary reflectors of order M

        Q  =  H(1) H(2) . . . H(k)

        as returned by SGEQRF.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix Q. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix Q. M >= n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines the
    matrix Q. n >= k >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the i-th column must contain the vector which
    defines the elementary reflector H(i), for i = 1,2,...,k, as
    returned by CGEQRF in the first k columns of its array
    argument a. \n
    On exit, the m-by-n matrix Q.

    * @param[in] lda
    lda is int* \n
    The first dimension of the array a. lda >= max(1,m).

    * @param[in] tau
    tau is scomplex/dcomplex array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SGEQRF.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,n).
    For optimum performance lwork >= n*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument has an illegal value
    *  */
template< typename T >
int ungqr( int* m, int* n, int* k, T*   a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return ungqr( m, n, k, a, lda, tau, work, lwork, info );
}

/*! @brief Apply Q or Q' from QR factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from QR factorization
        Overwrite the general real m-by-n matrix c with

        side = 'L'  side = 'R'
        trans = 'N':   Q * C C * Q
        trans = 'T':   Q**T * C C * Q**T

        where Q is a real orthogonal matrix defined as the product of k elementary reflectors

        Q = H(1) H(2) . . . H(k)

        as returned by SGEQRF. Q is of order M if side = 'L' and of order N if side = 'R'.
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**T from the Left; \n
    = 'R': apply Q or Q**T from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'T':  Transpose, apply Q**T.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines
    the matrix Q. \n
    If side = 'L', m >= k >= 0; \n
    if side = 'R', n >= k >= 0.

    * @param[in] a
    a is float/double array, dimension (lda,k) \n
    The i-th column must contain the vector which defines the
    elementary reflector H(i), for i = 1,2,...,k, as returned by
    SGEQRF in the first k columns of its array argument a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. \n
    If side = 'L', lda >= max(1,m); \n
    if side = 'R', lda >= max(1,n).

    * @param[in] tau
    tau is float/double array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SGEQRF.

    * @param[in,out] c
    c is float/double array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If side = 'L', lwork >= max(1,N); \n
    if side = 'R', lwork >= max(1,m). \n
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int ormqr( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
{
  return ormqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
}

/*! @brief Apply Q or Q' from QR factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from QR factorization
        Overwrite the general complex m-by-n matrix c with

        side = 'L'  side = 'R'
        trans = 'N':   Q * C C * Q
        trans = 'H':   Q**H * C C * Q**H

        where Q is a complex unitary matrix defined as the product of k elementary reflectors

        Q = H(1) H(2) . . . H(k)

        as returned by CGEQRF. Q is of order M if side = 'L' and of order N if side = 'R'.
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**H from the Left; \n
    = 'R': apply Q or Q**H from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'H':  Transpose, apply Q**H.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines
    the matrix Q. \n
    If side = 'L', m >= k >= 0; \n
    if side = 'R', n >= k >= 0.

    * @param[in] a
    a is scomplex/dcomplex array, dimension (lda,k) \n
    The i-th column must contain the vector which defines the
    elementary reflector H(i), for i = 1,2,...,k, as returned by
    CGEQRF in the first k columns of its array argument a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. \n
    If side = 'L', lda >= max(1,m); \n
    if side = 'R', lda >= max(1,n).

    * @param[in] tau
    tau is scomplex/dcomplex array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SGEQRF.

    * @param[in,out] c
    c is scomplex/dcomplex array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If side = 'L', lwork >= max(1,n); \n
    if side = 'R', lwork >= max(1,m). \n
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int unmqr( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
{
  return unmqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
}

/*! @brief Multiply a general matrix by the orthogonal matrix from a QR factorization
    *   determined by sgeqrf (unblocked algorithm).
    *

    * @details
    * \b Purpose:
    * \verbatim
        Multiply a general matrix by the orthogonal matrix from a QR factorization determined by
        sgeqrf (unblocked algorithm).
        Overwrite the general real m by n matrix c with

        Q * C  if side = 'L' and trans = 'N', or

        Q**T* C  if side = 'L' and trans = 'T', or

        C * Q  if side = 'R' and trans = 'N', or

        C * Q**T if side = 'R' and trans = 'T',

        where Q is a real orthogonal matrix defined as the product of k elementary reflectors

        Q = H(1) H(2) . . . H(k)

        as returned by SGEQRF. Q is of order m if side = 'L' and of order n if side = 'R'.
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**T from the Left; \n
    = 'R': apply Q or Q**T from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'T':  Transpose, apply Q**T.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines
    the matrix Q. \n
    If side = 'L', m >= k >= 0; \n
    if side = 'R', n >= k >= 0.

    * @param[in] a
    a is float/double array, dimension (lda,k) \n
    The i-th column must contain the vector which defines the
    elementary reflector H(i), for i = 1,2,...,k, as returned by
    SGEQRF in the first k columns of its array argument a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. \n
    If side = 'L', lda >= max(1,m); \n
    if side = 'R', lda >= max(1,n).

    * @param[in] tau
    tau is float/double array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SGEQRF.

    * @param[in,out] c
    c is float/double array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int orm2r( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* info )
{
  return orm2r( side, trans,  m, n, k, a, lda, tau, c, ldc, work, info );
}

/*! @brief Multiply a general matrix by the orthogonal matrix from a QR factorization
    *   determined by sgeqrf (unblocked algorithm).
    *

    * @details
    * \b Purpose:
    * \verbatim
        Multiply a general matrix by the orthogonal matrix from a QR factorization determined by
        sgeqrf (unblocked algorithm).
        Overwrite the general complex m by n matrix c with

        Q * C  if side = 'L' and trans = 'N', or

        Q**H* C  if side = 'L' and trans = 'C', or

        C * Q  if side = 'R' and trans = 'N', or

        C * Q**H if side = 'R' and trans = 'C',

        where Q is a complex unitary matrix defined as the product of k elementary reflectors

        Q = H(1) H(2) . . . H(k)

        as returned by SGEQRF. Q is of order m if side = 'L' and of order n if side = 'R'.
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**H from the Left; \n
    = 'R': apply Q or Q**H from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'C':  Transpose, apply Q**H.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines
    the matrix Q. \n
    If side = 'L', m >= k >= 0; \n
    if side = 'R', n >= k >= 0.

    * @param[in] a
    a is scomplex/dcomplex array, dimension (lda,k) \n
    The i-th column must contain the vector which defines the
    elementary reflector H(i), for i = 1,2,...,k, as returned by
    CGEQRF in the first k columns of its array argument a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. \n
    If side = 'L', lda >= max(1,m); \n
    if side = 'R', lda >= max(1,n).

    * @param[in] tau
    tau is scomplex/dcomplex array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by CGEQRF.

    * @param[in,out] c
    c is scomplex/dcomplex array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int unm2r( char* side, char* trans, int* m, int* n, int* k, T*   a, int* lda, T* tau, T* c, int* ldc, T* work, int* info )
{
  return unm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
}

/*! @brief Form Q from LQ factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Generate an m-by-n real matrix Q with orthonormal rows, which is defined as the first M
        rows of a product of K elementary reflectors of order N

        Q  =  H(k) . . . H(2) H(1)

        as returned by SGELQF.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix Q. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix Q. n >= m.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines the
    matrix Q. m >= k >= 0.

    * @param[in,out] a
    a is float/doublearray, dimension (lda,n) \n
    On entry, the i-th row must contain the vector which defines
    the elementary reflector H(i), for i = 1,2,...,k, as returned
    by SGELQF in the first k rows of its array argument a.
    On exit, the m-by-n matrix Q.

    * @param[in] lda
    lda is int* \n
    The first dimension of the array a. lda >= max(1,m).

    * @param[in] tau
    tau is float/double array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SGELQF.

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,m). \n
    For optimum performance lwork >= m*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument has an illegal value
    *  */
template< typename T >
int orglq( int* m, int* n, int* k, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return orglq( m, n, k, a, lda, tau, work, lwork, info );
}

/*! @brief Form Q from LQ factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Generate an m-by-n complex matrix Q with orthonormal rows, which is defined as the first
        M rows of a product of K elementary reflectors of order N

        Q  =  H(k) . . . H(2) H(1)

        as returned by CGELQF.
    \endverbatim

    * @param[in] m
    m is int* \n
    The number of rows of the matrix Q. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix Q. n >= m.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines the
    matrix Q. m >= k >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the i-th row must contain the vector which defines
    the elementary reflector H(i), for i = 1,2,...,k, as returned
    by CGELQF in the first k rows of its array argument a. \n
    On exit, the m-by-n matrix Q.

    * @param[in] lda
    lda is int* \n
    The first dimension of the array a. lda >= max(1,m).

    * @param[in] tau
    tau is scomplex/dcomplex array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by CGELQF.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,m). \n
    For optimum performance lwork >= m*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument has an illegal value
    *  */
template< typename T >
int unglq( int* m, int* n, int* k, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return unglq( m, n, k, a, lda, tau, work, lwork, info );
}

/*! @brief Apply Q or Q' from LQ factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from LQ factorization. Overwrite the general real m-by-n matrix c with

        side = 'L'  side = 'R'
        trans = 'N':   Q * C C * Q
        trans = 'T':   Q**T * C C * Q**T

        where Q is a real orthogonal matrix defined as the product of k elementary reflectors

        Q = H(k) . . . H(2) H(1)

        as returned by SGELQF. Q is of order M if side = 'L' and of order N if side = 'R'.
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**T from the Left; \n
    = 'R': apply Q or Q**T from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'T':  Transpose, apply Q**T.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines
    the matrix Q. \n
    If side = 'L', m >= k >= 0; \n
    if side = 'R', n >= k >= 0.

    * @param[in] a
    a is float/double/ array, dimension \n
    (lda,m) if side = 'L', \n
    (lda,n) if side = 'R' \n
    The i-th row must contain the vector which defines the
    elementary reflector H(i), for i = 1,2,...,k, as returned by
    SGELQF in the first k rows of its array argument a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,k).

    * @param[in] tau
    tau is float/double array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SGELQF.

    * @param[in,out] c
    c is float/double array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If side = 'L', lwork >= max(1,N); \n
    if side = 'R', lwork >= max(1,m). \n
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int ormlq( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
{
  return ormlq( side, trans,  m, n, k, a, lda, tau, c, ldc, work, lwork, info );
}

/*! @brief Apply Q or Q' from LQ factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from LQ factorization. Overwrite the general complex m-by-n matrix c with

        side = 'L'  side = 'R'
        trans = 'N':   Q * C C * Q
        trans = 'C':   Q**H * C C * Q**H

        where Q is a complex unitary matrix defined as the product of k elementary reflectors

        Q = H(k)**H . . . H(2)**H H(1)**H

        as returned by CGELQF. Q is of order M if side = 'L' and of order N if side = 'R'.
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**H from the Left; \n
    = 'R': apply Q or Q**H from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'C':  Conjugate transpose, apply Q**H.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines
    the matrix Q. \n
    If side = 'L', m >= k >= 0; \n
    if side = 'R', n >= k >= 0.

    * @param[in] a
    a is scomplex/dcomplex array, dimension \n
    (lda,m) if side = 'L', \n
    (lda,n) if side = 'R' \n
    The i-th row must contain the vector which defines the
    elementary reflector H(i), for i = 1,2,...,k, as returned by
    CGELQF in the first k rows of its array argument a.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,k).
    * @param[in] tau
    tau is scomplex/dcomplex array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by CGELQF.

    * @param[in,out] c
    c is scomplex/dcomplex array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If side = 'L', lwork >= max(1,N); \n
    if side = 'R', lwork >= max(1,m). \n
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int unmlq( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
{
  return   unmlq( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
}

/*! @brief Apply Q or Q' from LQ factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from LQ factorization. Overwrite the general real m-by-n matrix c with

        Q * C  if side = 'L' and trans = 'N', or

        Q**T* C  if side = 'L' and trans = 'T', or

        C * Q  if side = 'R' and trans = 'N', or

        C * Q**T if side = 'R' and trans = 'T',

        where Q is a real orthogonal matrix defined as the product of k elementary reflectors

        Q = H(k) . . . H(2) H(1)

        as returned by SGELQF. Q is of order m if side = 'L' and of order n
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**T from the Left; \n
    = 'R': apply Q or Q**T from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'T':  Transpose, apply Q**T.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines
    the matrix Q. \n
    If side = 'L', m >= k >= 0; \n
    if side = 'R', n >= k >= 0. \n

    * @param[in] a
    a is float/double array, dimension \n
    (lda,m) if side = 'L', \n
    (lda,n) if side = 'R' \n
    The i-th row must contain the vector which defines the
    elementary reflector H(i), for i = 1,2,...,k, as returned by
    SGELQF in the first k rows of its array argument a. \n
    a is modified by the routine but restored on exit.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,k).

    * @param[in] tau
    tau is float/double array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SGELQF.

    * @param[in,out] c
    c is float/double array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is float/double array, dimension: \n
    (N) if side = 'L', \n
    (M) if side = 'R' \n

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int orml2( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* info )
{
  return orml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
}

/*! @brief Apply Q or Q' from LQ factorization
    *

    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from LQ factorization. Overwrite the general complex m-by-n matrix c with

        Q * C  if side = 'L' and trans = 'N', or

        Q**H* C  if side = 'L' and trans = 'C', or

        C * Q  if side = 'R' and trans = 'N', or

        C * Q**H if side = 'R' and trans = 'C',

        where Q is a complex unitary matrix defined as the product of k elementary reflectors

        Q = H(k)**H . . . H(2)**H H(1)**H

        as returned by CGELQF. Q is of order m if side = 'L' and of order n
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**T from the Left; \n
    = 'R': apply Q or Q**T from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'T':  Transpose, apply Q**T.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    The number of elementary reflectors whose product defines
    the matrix Q. \n
    If side = 'L', m >= k >= 0; \n
    if side = 'R', n >= k >= 0. \n

    * @param[in] a
    a is scomplex/dcomplex array, dimension \n
    (lda,m) if side = 'L', \n
    (lda,n) if side = 'R' \n
    The i-th row must contain the vector which defines the
    elementary reflector H(i), for i = 1,2,...,k, as returned by
    CGELQF in the first k rows of its array argument a. \n
    a is modified by the routine but restored on exit.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,k).

    * @param[in] tau
    tau is scomplex/dcomplex array, dimension (k) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by CGELQF.

    * @param[in,out] c
    c is scomplex/dcomplex array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is scomplex/dcomplex array, dimension: \n
    (N) if side = 'L', \n
    (M) if side = 'R' \n

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int unml2( char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T*   c, int* ldc, T*   work, int* info )
{
  return unml2( side, trans, m, n, k, a, lda, tau, c, ldc, work, info );
}

/*! @brief Form Q from tridiagonal reduction
    *

    * @details
    * \b Purpose:
    * \verbatim
    Form Q from tridiagonal reduction. Generate a real orthogonal matrix Q which is defined as
    the product of n-1 elementary reflectors of order M, as returned by SSYTRD:

        if uplo = 'U', Q = H(n-1) . . . H(2) H(1),

        if uplo = 'L', Q = H(1) H(2) . . . H(n-1).
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U': Upper triangle of A contains elementary reflectors
    from SSYTRD; \n
    = 'L': Lower triangle of A contains elementary reflectors
    from SSYTRD.

    * @param[in] m
    m is int* \n
    The order of the matrix Q. m >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,m) \n
    On entry, the vectors which define the elementary reflectors,
    as returned by SSYTRD. \n
    On exit, the m-by-m orthogonal matrix Q.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,m).

    * @param[in] tau
    tau is float/double array, dimension (m-1) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SSYTRD.

    * @param[out] work
    work is float/double, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,m-1). \n
    For optimum performance lwork >= (m-1)*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int orgtr( char* uplo, int* m, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return orgtr( uplo, m, a, lda, tau, work, lwork, info );
}

/*! @brief Form Q from tridiagonal reduction
    *

    * @details
    * \b Purpose:
    * \verbatim
    Form Q from tridiagonal reduction. Generate a complex unitary matrix Q which is defined as
    the product of n-1 elementary reflectors of order M, as returned by CHETRD:

        if uplo = 'U', Q = H(n-1) . . . H(2) H(1),

        if uplo = 'L', Q = H(1) H(2) . . . H(n-1).
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U': Upper triangle of A contains elementary reflectors
    from CHETRD; \n
    = 'L': Lower triangle of A contains elementary reflectors
    from CHETRD.

    * @param[in] m
    m is int* \n
    The order of the matrix Q. m >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,m) \n
    On entry, the vectors which define the elementary reflectors,
    as returned by CHETRD. \n
    On exit, the m-by-m orthogonal matrix Q.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,m).

    * @param[in] tau
    tau is scomplex/dcomplex array, dimension (m-1) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by CHETRD.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,m-1). \n
    For optimum performance lwork >= (m-1)*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int ungtr( char* uplo, int* m, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return ungtr( uplo, m, a, lda, tau, work, lwork, info );
}

/*! @brief Apply Q or Q' from tridiagonal reduction
    *

    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from tridiagonal reduction. Overwrite the general real m-by-n matrix c with

        side = 'L'  side = 'R'
        trans = 'N':   Q * C C * Q
        trans = 'T':   Q**T * C C * Q**T

        where Q is a real orthogonal matrix of order nq, with nq = m if side = 'L' and nq = n if
        side = 'R'. Q is defined as the product of nq-1 elementary reflectors, as returned by SSYTRD:

        if uplo = 'U', Q = H(nq-1) . . . H(2) H(1);

        if uplo = 'L', Q = H(1) H(2) . . . H(nq-1).
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**T from the Left; \n
    = 'R': apply Q or Q**T from the Right.

    * @param[in] uplo
    uplo is char* \n
    = 'U': Upper triangle of A contains elementary reflectors
    from SSYTRD; \n
    = 'L': Lower triangle of A contains elementary reflectors
    from SSYTRD.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'T':  Transpose, apply Q**T.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] a
    a is float/double array, dimension \n
    (lda,m) if side = 'L' \n
    (lda,n) if side = 'R' \n
    The vectors which define the elementary reflectors, as
    returned by SSYTRD.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. \n
    lda >= max(1,m) if side = 'L'; lda >= max(1,n) if side = 'R'.

    * @param[in] tau
    tau is float/double array, dimension \n
    (m-1) if side = 'L' \n
    (n-1) if side = 'R' \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by SSYTRD.

    * @param[in,out] c
    c is float/double array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If side = 'L', lwork >= max(1,N); \n
    if side = 'R', lwork >= max(1,m). \n
    For optimum performance lwork >= N*NB if side = 'L', and
    lwork >= m*NB if side = 'R', \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int ormtr( char* side, char* uplo, char* trans, int* m, int* n, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
{
  return ormtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
}

/*! @brief Apply Q or Q' from tridiagonal reduction
    *

    * @details
    * \b Purpose:
    * \verbatim
        Apply Q or Q' from tridiagonal reduction. Overwrite the general complex m-by-n matrix c with

            side = 'L'  side = 'R'
            trans = 'N':   Q * C C * Q
            trans = 'C':   Q**H * C C * Q**H

        where Q is a complex unitary matrix of order nq, with nq = m if side = 'L' and nq = n if
        side = 'R'. Q is defined as the product of nq-1 elementary reflectors, as returned by CHETRD:

            if uplo = 'U', Q = H(nq-1) . . . H(2) H(1);

            if uplo = 'L', Q = H(1) H(2) . . . H(nq-1).
    \endverbatim

    * @param[in] side
    side is char* \n
    = 'L': apply Q or Q**H from the Left; \n
    = 'R': apply Q or Q**H from the Right.

    * @param[in] uplo
    uplo is char* \n
    = 'U': Upper triangle of A contains elementary reflectors
    from CHETRD; \n
    = 'L': Lower triangle of A contains elementary reflectors
    from CHETRD.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q; \n
    = 'C':  Transpose, apply Q**C.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] a
    a is scomplex/dcomplex array, dimension \n
    (lda,m) if side = 'L' \n
    (lda,n) if side = 'R' \n
    The vectors which define the elementary reflectors, as
    returned by CHETRD.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. \n
    lda >= max(1,m) if side = 'L'; lda >= max(1,n) if side = 'R'.

    * @param[in] tau
    tau is scomplex/dcomplex array, dimension \n
    (m-1) if side = 'L' \n
    (n-1) if side = 'R' \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i), as returned by CHETRD.

    * @param[in,out] c
    c is scomplex/dcomplex array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If side = 'L', lwork >= max(1,N); \n
    if side = 'R', lwork >= max(1,m). \n
    For optimum performance lwork >= N*NB if side = 'L', and
    lwork >= m*NB if side = 'R', \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int unmtr( char* side, char* uplo, char* trans, int* m, int* n, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
{
  return unmtr( side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info );
}

/*! @brief Form Q from bidiagonal reduction
    *

    * @details
    * \b Purpose:
    * \verbatim
        Generate one of the real orthogonal matrices Q or P**T determined by SGEBRD when reducing
        a real matrix a to bidiagonal form: A = Q * B * P**T.  Q and P**T are defined as products
        of elementary reflectors H(i) or G(i) respectively.

        If vect = 'Q', A is assumed to have been an M-by-K matrix, and Q is of order M:
        if m >= k, Q = H(1) H(2) . . . H(k) and SORGBR returns the first n columns of Q,
        where m >= n >= k;
        if m < k, Q = H(1) H(2) . . . H(m-1) and SORGBR returns Q as an M-by-M matrix.

        If vect = 'P', A is assumed to have been a K-by-N matrix, and P**T is of order N:
        if k < n, P**T = G(k) . . . G(2) G(1) and SORGBR returns the first m rows of P**T,
        where n >= m >= k;
        if k >= n, P**T = G(n-1) . . . G(2) G(1) and SORGBR returns P**T as an n-by-n matrix.
    \endverbatim

    * @param[in] vect
    vect is char* \n
    Specifies whether the matrix Q or the matrix P**T is
    required, as defined in the transformation applied by SGEBRD:
    = 'Q':  generate Q; \n
    = 'P':  generate P**T.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix Q or P**T to be returned.
    m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix Q or P**T to be returned.
    n >= 0. \n
    If vect = 'Q', m >= n >= min(m,k); \n
    if vect = 'P', n >= m >= min(n,k).

    * @param[in] k
    k is int* \n
    If vect = 'Q', the number of columns in the original M-by-K
    matrix reduced by SGEBRD. \n
    If vect = 'P', the number of rows in the original K-by-N
    matrix reduced by SGEBRD. \n
    k >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the vectors which define the elementary reflectors,
    as returned by SGEBRD. \n
    On exit, the m-by-n matrix Q or P**T.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,m).

    * @param[in] tau
    tau is float/double array, dimension \n
    (min(m,k)) if vect = 'Q' \n
    (min(n,k)) if vect = 'P' \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i) or G(i), which determines Q or P**T, as
    returned by SGEBRD in its array argument tauq or taup.

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,min(m,n)). \n
    For optimum performance lwork >= min(m,n)*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int orgbr( char* vect, int* m, int* n, int* k, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return orgbr( vect, m, n, k, a, lda, tau, work, lwork, info );
}

/*! @brief Form Q from bidiagonal reduction
    *

    * @details
    * \b Purpose:
    * \verbatim
        Generate one of the complex unitary matrices Q or P**H determined by CGEBRD when reducing
        a complex matrix a to bidiagonal form: A = Q * B * P**H.  Q and P**H are defined as
        products of elementary reflectors H(i) or G(i) respectively.

        If vect = 'Q', A is assumed to have been an M-by-K matrix, and Q is of order M:
        if m >= k, Q = H(1) H(2) . . . H(k) and CUNGBR returns the first n columns of Q,
        where m >= n >= k;
        if m < k, Q = H(1) H(2) . . . H(m-1) and CUNGBR returns Q as an M-by-M matrix.

        If vect = 'P', A is assumed to have been a K-by-N matrix, and P**H is of order N:
        if k < n, P**H = G(k) . . . G(2) G(1) and CUNGBR returns the first m rows of P**H,
        where n >= m >= k;
        if k >= n, P**H = G(n-1) . . . G(2) G(1) and CUNGBR returns P**H as an n-by-n matrix.
    \endverbatim

    * @param[in] vect
    vect is char* \n
    Specifies whether the matrix Q or the matrix P**H is
    required, as defined in the transformation applied by CGEBRD: \n
    = 'Q':  generate Q; \n
    = 'P':  generate P**H.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix Q or P**H to be returned.
    m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix Q or P**H to be returned.
    n >= 0. \n
    If vect = 'Q', m >= n >= min(m,k); \n
    if vect = 'P', n >= m >= min(n,k).

    * @param[in] k
    k is int* \n
    If vect = 'Q', the number of columns in the original M-by-K
    matrix reduced by CGEBRD. \n
    If vect = 'P', the number of rows in the original K-by-N
    matrix reduced by CGEBRD.
    k >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the vectors which define the elementary reflectors,
    as returned by CGEBRD. \n
    On exit, the m-by-n matrix Q or P**H.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. lda >= max(1,m).

    * @param[in] tau
    tau is scomplex/dcomplex array, dimension \n
    (min(m,k)) if vect = 'Q' \n
    (min(n,k)) if vect = 'P' \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i) or G(i), which determines Q or P**H, as
    returned by CGEBRD in its array argument tauq or taup.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,min(m,n)). \n
    For optimum performance lwork >= min(m,n)*NB, \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int ungbr( char* vect, int* m, int* n, int* k, T* a, int* lda, T* tau, T* work, int* lwork, int* info )
{
  return ungbr( vect, m, n, k, a, lda, tau, work, lwork, info );
}

/*! @brief Apply Q or Q' from bidiagonal reduction
    *

    * @details
    * \b Purpose:
    * \verbatim
        If vect = 'Q', SORMBR overwrites the general real m-by-n matrix c with
            side = 'L'  side = 'R'
            trans = 'N':   Q * C C * Q
            trans = 'T':   Q**T * C C * Q**T

        If vect = 'P', SORMBR overwrites the general real m-by-n matrix c with
            side = 'L'  side = 'R'
            trans = 'N':   P * C C * P
            trans = 'T':   P**T * C C * P**T

        Here Q and P**T are the orthogonal matrices determined by SGEBRD when reducing a real
        matrix a to bidiagonal form: A = Q * B * P**T. Q and P**T are defined as products of
        elementary reflectors H(i) and G(i) respectively.

        Let nq = m if side = 'L' and nq = n if side = 'R'. Thus nq is the order of the orthogonal
        matrix Q or P**T that is applied.

        If vect = 'Q', A is assumed to have been an NQ-by-K matrix:
            if nq >= k, Q = H(1) H(2) . . . H(k);
            if nq < k, Q = H(1) H(2) . . . H(nq-1).

        If vect = 'P', A is assumed to have been a K-by-NQ matrix:
            if k < nq, P = G(1) G(2) . . . G(k);
            if k >= nq, P = G(1) G(2) . . . G(nq-1).
    \endverbatim

    * @param[in] vect
    vect is char* \n
    = 'Q': apply Q or Q**T; \n
    = 'P': apply P or P**T.

    * @param[in] side
    side is char* \n
    = 'L': apply Q, Q**T, P or P**T from the Left; \n
    = 'R': apply Q, Q**T, P or P**T from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q  or P; \n
    = 'T':  Transpose, apply Q**T or P**T.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    If vect = 'Q', the number of columns in the original
    matrix reduced by SGEBRD. \n
    If vect = 'P', the number of rows in the original
    matrix reduced by SGEBRD. \n
    k >= 0.

    * @param[in] a
    a is float/double array, dimension \n
    (lda,min(nq,k)) if vect = 'Q' \n
    (lda,nq)  if vect = 'P' \n
    The vectors which define the elementary reflectors H(i) and
    G(i), whose products determine the matrices Q and P, as
    returned by SGEBRD.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. \n
    If vect = 'Q', lda >= max(1,nq); \n
    if vect = 'P', lda >= max(1,min(nq,k)).

    * @param[in] tau
    tau is float/double array, dimension (min(nq,k)) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i) or G(i) which determines Q or P, as returned
    by SGEBRD in the array argument tauq or taup.

    * @param[in,out] c
    c is float/double array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
    or P*C or P**T*C or C*P or C*P**T.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If side = 'L', lwork >= max(1,N); \n
    if side = 'R', lwork >= max(1,m). \n
    For optimum performance lwork >= N*NB if side = 'L', and
    lwork >= m*NB if side = 'R', \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int ormbr( char* vect, char* side, char* trans, int* m, int* n, int* k, T* a, int* lda,T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
{
  return ormbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
}

    /*! @brief Apply Q or Q' from bidiagonal reduction
    *

    * @details
    * \b Purpose:
    * \verbatim
        If vect = 'Q', CUNMBR overwrites the general complex m-by-n matrix c with
        side = 'L'  side = 'R'
            trans = 'N':   Q * C C * Q
            trans = 'C':   Q**H * C C * Q**H

        If vect = 'P', CUNMBR overwrites the general complex m-by-n matrix c with
        side = 'L'  side = 'R'
            trans = 'N':   P * C C * P
            trans = 'C':   P**H * C C * P**H

        Here Q and P**H are the orthogonal matrices determined by CGEBRD when reducing a complex
        matrix a to bidiagonal form: A = Q * B * P**H. Q and P**H are defined as products of
        elementary reflectors H(i) and G(i) respectively.

        Let nq = m if side = 'L' and nq = n if side = 'R'. Thus nq is the order of the orthogonal
        matrix Q or P**H that is applied.

        If vect = 'Q', A is assumed to have been an NQ-by-K matrix:
            if nq >= k, Q = H(1) H(2) . . . H(k);
            if nq < k, Q = H(1) H(2) . . . H(nq-1).

        If vect = 'P', A is assumed to have been a K-by-NQ matrix:
            if k < nq, P = G(1) G(2) . . . G(k);
            if k >= nq, P = G(1) G(2) . . . G(nq-1).
    \endverbatim

    * @param[in] vect
    vect is char* \n
    = 'Q': apply Q or Q**H; \n
    = 'P': apply P or P**H.

    * @param[in] side
    side is char* \n
    = 'L': apply Q, Q**H, P or P**H from the Left; \n
    = 'R': apply Q, Q**H, P or P**H from the Right.

    * @param[in] trans
    trans is char* \n
    = 'N':  No transpose, apply Q  or P; \n
    = 'C':  Transpose, apply Q**H or P**H.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix c. m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix c. n >= 0.

    * @param[in] k
    k is int* \n
    If vect = 'Q', the number of columns in the original
    matrix reduced by CGEBRD. \n
    If vect = 'P', the number of rows in the original
    matrix reduced by CGEBRD. \n
    k >= 0.

    * @param[in] a
    a is scomplex/dcomplex array, dimension \n
    (lda,min(nq,k)) if vect = 'Q' \n
    (lda,nq)  if vect = 'P' \n
    The vectors which define the elementary reflectors H(i) and
    G(i), whose products determine the matrices Q and P, as
    returned by CGEBRD.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a. \n
    If vect = 'Q', lda >= max(1,nq); \n
    if vect = 'P', lda >= max(1,min(nq,k)).

    * @param[in] tau
    tau is scomplex/dcomplex array, dimension (min(nq,k)) \n
    tau(i) must contain the scalar factor of the elementary
    reflector H(i) or G(i) which determines Q or P, as returned
    by CGEBRD in the array argument tauq or taup.

    * @param[in,out] c
    c is scomplex/dcomplex array, dimension (ldc,n) \n
    On entry, the m-by-n matrix c. \n
    On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q
    or P*C or P**H*C or C*P or C*P**H.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. ldc >= max(1,m).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If side = 'L', lwork >= max(1,N); \n
    if side = 'R', lwork >= max(1,m). \n
    For optimum performance lwork >= N*NB if side = 'L', and
    lwork >= m*NB if side = 'R', \n
    where NB is the optimal blocksize. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value
    *  */
template< typename T >
int unmbr( char* vect, char* side, char* trans, int* m, int* n, int* k, T* a, int* lda, T* tau, T* c, int* ldc, T* work, int* lwork, int* info )
{
  return unmbr( vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info );
}

/*! @brief Tridiagonal QR algorithm
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of all eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal
        matrix using the implicit QL or QR method. The eigenvectors of a full or band symmetric
        matrix can also be found if SSYTRD or SSPTRD or SSBTRD has been used to reduce this matrix
        to tridiagonal form.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only. \n
    = 'V':  Compute eigenvalues and eigenvectors of the original
    symmetric matrix.  On entry, Z must contain the
    orthogonal matrix used to reduce the original matrix
    to tridiagonal form. \n
    = 'I':  Compute eigenvalues and eigenvectors of the
    tridiagonal matrix.  Z is initialized to the identity
    matrix.

    * @param[in] n
    n is int* \n
    The order of the matrix.  n >= 0.

    * @param[in,out] d
    d is float/double array, dimension (n) \n
    On entry, the diagonal elements of the tridiagonal matrix. \n
    On exit, if info = 0, the eigenvalues in ascending order.

    * @param[in,out] e
    e is float/double array, dimension (n-1) \n
    On entry, the (n-1) subdiagonal elements of the tridiagonal
    matrix. \n
    On exit, E has been destroyed.

    * @param[in,out] z
    z is float/double array, dimension (ldz, n) \n
    On entry, if  jobz = 'V', then Z contains the orthogonal
    matrix used in the reduction to tridiagonal form. \n
    On exit, if info = 0, then if  jobz = 'V', Z contains the
    orthonormal eigenvectors of the original symmetric matrix,
    and if jobz = 'I', Z contains the orthonormal eigenvectors
    of the symmetric tridiagonal matrix. \n
    If jobz = 'N', then Z is not referenced.

    * @param[in] ldz
    ldz is int* \n
    The leading dimension of the array Z.  ldz >= 1, and if
    eigenvectors are desired, then  ldz >= max(1,n).

    * @param[out] work
    work is float/double array, dimension (max(1,2*n-2)) \n
    If jobz = 'N', then work is not referenced.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  the algorithm has failed to find all the eigenvalues in
    a total of 30*N iterations; if info = i, then i
    elements of E have not converged to zero; on exit, D
    and E contain the elements of a symmetric tridiagonal
    matrix which is orthogonally similar to the original
    matrix.
    *  */
template< typename T >
int steqr( char* jobz, int* n, T* d, T* e, T* z, int* ldz, T*  work, int* info )
{
  return steqr( jobz, n, d, e, z, ldz, work, info );
}


/*! @brief Tridiagonal QR algorithm
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of all eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal
        matrix using the implicit QL or QR method. The eigenvectors of a full or band complex
        Hermitian matrix can also be found if CHETRD or CHPTRD or CHBTRD has been used to reduce
        this matrix to tridiagonal form.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only. \n
    = 'V':  Compute eigenvalues and eigenvectors of the original
            Hermitian matrix.  On entry, Z must contain the
            unitary matrix used to reduce the original matrix
            to tridiagonal form. \n
    = 'I':  Compute eigenvalues and eigenvectors of the
            tridiagonal matrix.  Z is initialized to the identity
            matrix.

    * @param[in] n
    n is int* \n
    The order of the matrix.  n >= 0.

    * @param[in,out] d
    d is float/double array, dimension (n) \n
    On entry, the diagonal elements of the tridiagonal matrix. \n
    On exit, if info = 0, the eigenvalues in ascending order.

    * @param[in,out] e
    e is float/double array, dimension (n-1) \n
    On entry, the (n-1) subdiagonal elements of the tridiagonal
    matrix. \n
    On exit, E has been destroyed.

    * @param[in,out] z
    z is scomplex/dcomplex array, dimension (ldz, n) \n
    On entry, if  jobz = 'V', then Z contains the unitary
    matrix used in the reduction to tridiagonal form. \n
    On exit, if info = 0, then if  jobz = 'V', Z contains the
    orthonormal eigenvectors of the original Hermitian matrix,
    and if jobz = 'I', Z contains the orthonormal eigenvectors
    of the symmetric tridiagonal matrix. \n
    If jobz = 'N', then Z is not referenced.

    * @param[in] ldz
    ldz is int* \n
    The leading dimension of the array Z.  ldz >= 1, and if
    eigenvectors are desired, then  ldz >= max(1,n).

    * @param[out] work
    work is float/double array, dimension (max(1,2*n-2)) \n
    If jobz = 'N', then work is not referenced.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  the algorithm has failed to find all the eigenvalues in
        a total of 30*N iterations; if info = i, then i
        elements of E have not converged to zero; on exit, D
        and E contain the elements of a symmetric tridiagonal
        matrix which is orthogonally similar to the original
        matrix.
    *  */
template< typename Ta, typename Tb >
int steqr( char* jobz, int* n, Tb* d, Tb* e, Ta* z, int* ldz, Tb* work, int* info )
{
  return steqr( jobz, n, d, e, z, ldz, work, info );
}

/*! @brief Tridiagonal divide-and-conquer algorithm
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of all eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal
        matrix using the divide and conquer method. The eigenvectors of a full or band real
        symmetric matrix can also be found if SSYTRD or SSPTRD or SSBTRD has been used to reduce
        this matrix to tridiagonal form.

        This code makes very mild assumptions about floating point arithmetic. It will work on
        machines with a guard digit in add/subtract, or on those binary machines without guard
        digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could
        conceivably fail on hexadecimal or decimal machines without guard digits, but we know of
        none.  See SLAED3 for details.
    \endverbatim

    * @param[in] compz
    compz is char* \n
    = 'N':  Compute eigenvalues only. \n
    = 'I':  Compute eigenvectors of tridiagonal matrix also. \n
    = 'V':  Compute eigenvectors of original dense symmetric
    matrix also.  On entry, Z contains the orthogonal
    matrix used to reduce the original matrix to
    tridiagonal form.

    * @param[in] n
    n is int* \n
    The dimension of the symmetric tridiagonal matrix.  n >= 0.

    * @param[in,out] d
    d is float/double array, dimension (n) \n
    On entry, the diagonal elements of the tridiagonal matrix. \n
    On exit, if info = 0, the eigenvalues in ascending order.

    * @param[in,out] e
    e is float/double array, dimension (n-1) \n
    On entry, the subdiagonal elements of the tridiagonal matrix. \n
    On exit, E has been destroyed.

    * @param[in,out] z
    z is float/double array, dimension (ldz,n) \n
    On entry, if compz = 'V', then Z contains the orthogonal
    matrix used in the reduction to tridiagonal form. \n
    On exit, if info = 0, then if compz = 'V', Z contains the
    orthonormal eigenvectors of the original symmetric matrix,
    and if compz = 'I', Z contains the orthonormal eigenvectors
    of the symmetric tridiagonal matrix. \n
    If  compz = 'N', then Z is not referenced.

    * @param[in] ldz
    ldz is int* \n
    The leading dimension of the array Z.  ldz >= 1. \n
    If eigenvectors are desired, then ldz >= max(1,n).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If compz = 'N' or N <= 1 then lwork must be at least 1. \n
    If compz = 'V' and N > 1 then lwork must be at least \n
    ( 1 + 3*N + 2*n*lg N + 4*N**2 ), \n
    where lg( N ) = smallest integer k such
    that 2**k >= N. \n
    If compz = 'I' and N > 1 then lwork must be at least
    ( 1 + 4*N + N**2 ). \n
    Note that for compz = 'I' or 'V', then if N is less than or
    equal to the minimum divide size, usually 25, then lwork need
    only be max(1,2*(N-1)). \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] iwork
    iwork is int array, dimension (max(1,liwork)) \n
    On exit, if info = 0, iwork(1) returns the optimal liwork.

    * @param[in] liwork
    liwork is int* \n
    The dimension of the array iwork. \n
    If compz = 'N' or N <= 1 then liwork must be at least 1. \n
    If compz = 'V' and N > 1 then liwork must be at least \n
    ( 6 + 6*N + 5*N*lg N ). \n
    If compz = 'I' and N > 1 then liwork must be at least \n
    ( 3 + 5*N ). \n
    Note that for compz = 'I' or 'V', then if N is less than or
    equal to the minimum divide size, usually 25, then liwork
    need only be 1. \n \n
    If liwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal size of the iwork array,
    returns this value as the first entry of the iwork array, and
    no error message related to liwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit. \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  The algorithm failed to compute an eigenvalue while
    working on the submatrix lying in rows and columns
    info/(N+1) through mod(info,N+1).
    *  */
template< typename T >
int stedc( char* compz, int* n, T* d, T* e, T* z, int* ldz, T* work, int* lwork, int* iwork, int* liwork, int* info )
{
  return stedc( compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info );
}

/*! @brief Tridiagonal divide-and-conquer algorithm
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of all eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal
        matrix using the divide and conquer method. The eigenvectors of a full or band complex
        Hermitian matrix can also be found if CHETRD or CHPTRD or CHBTRD has been used to reduce
        this matrix to tridiagonal form.

        This code makes very mild assumptions about floating point arithmetic. It will work on
        machines with a guard digit in add/subtract, or on those binary machines without guard
        digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could
        conceivably fail on hexadecimal or decimal machines without guard digits, but we know of
        none.  See SLAED3 for details.
    \endverbatim

    * @param[in] compz
    compz is char* \n
    = 'N':  Compute eigenvalues only. \n
    = 'I':  Compute eigenvectors of tridiagonal matrix also. \n
    = 'V':  Compute eigenvectors of original Hermitian matrix
            also.  On entry, Z contains the unitary matrix used
            to reduce the original matrix to tridiagonal form.

    * @param[in] n
    n is int* \n
    The dimension of the symmetric tridiagonal matrix.  n >= 0.

    * @param[in,out] d
    d is float/double array, dimension (n) \n
    On entry, the diagonal elements of the tridiagonal matrix. \n
    On exit, if info = 0, the eigenvalues in ascending order.

    * @param[in,out] e
    e is float/double array, dimension (n-1) \n
    On entry, the subdiagonal elements of the tridiagonal matrix. \n
    On exit, E has been destroyed.

    * @param[in,out] z
    z is scomplex/dcomplex array, dimension (ldz,n) \n
    On entry, if compz = 'V', then Z contains the unitary
    matrix used in the reduction to tridiagonal form. \n
    On exit, if info = 0, then if compz = 'V', Z contains the
    orthonormal eigenvectors of the original Hermitian matrix,
    and if compz = 'I', Z contains the orthonormal eigenvectors
    of the symmetric tridiagonal matrix. \n
    If  compz = 'N', then Z is not referenced.

    * @param[in] ldz
    ldz is int* \n
    The leading dimension of the array Z.  ldz >= 1. \n
    If eigenvectors are desired, then ldz >= max(1,n).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If compz = 'N' or 'I', or N <= 1, lwork must be at least 1. \n
    If compz = 'V' and N > 1, lwork must be at least N*N. \n
    Note that for compz = 'V', then if N is less than or
    equal to the minimum divide size, usually 25, then lwork need
    only be 1. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal sizes of the work, rwork and
    iwork arrays, returns these values as the first entries of
    the work, rwork and iwork arrays, and no error message
    related to lwork or lrwork or liwork is issued by XERBLA.

    * @param[out] rwork
    rwork is float/double array, dimension (max(1,lrwork)) \n
    On exit, if info = 0, rwork(1) returns the optimal lrwork.

    * @param[in] lrwork
    lrwork is int* \n
    The dimension of the array rwork. \n
    If compz = 'N' or N <= 1, lrwork must be at least 1. \n
    If compz = 'V' and N > 1, lrwork must be at least \n
                1 + 3*N + 2*n*lg N + 4*N**2 , \n
                where lg( N ) = smallest integer k such
                that 2**k >= N. \n
    If compz = 'I' and N > 1, lrwork must be at least \n
                1 + 4*N + 2*n**2 . \n
    Note that for compz = 'I' or 'V', then if N is less than or
    equal to the minimum divide size, usually 25, then lrwork
    need only be max(1,2*(N-1)). \n \n
    If lrwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal sizes of the work, rwork
    and iwork arrays, returns these values as the first entries
    of the work, rwork and iwork arrays, and no error message
    related to lwork or lrwork or liwork is issued by XERBLA.

    @param[out] iwork
    iwork is int array, dimension (max(1,liwork)) \n
    On exit, if info = 0, iwork(1) returns the optimal liwork.

    * @param[in] liwork
    liwork is int* \n
    The dimension of the array iwork. \n
    If compz = 'N' or N <= 1, liwork must be at least 1. \n
    If compz = 'V' or N > 1,  liwork must be at least \n
                            6 + 6*N + 5*N*lg N. \n
    If compz = 'I' or N > 1,  liwork must be at least \n
                            3 + 5*N . \n
    Note that for compz = 'I' or 'V', then if N is less than or
    equal to the minimum divide size, usually 25, then liwork
    need only be 1. \n \n
    If liwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal sizes of the work, rwork
    and iwork arrays, returns these values as the first entries
    of the work, rwork and iwork arrays, and no error message
    related to lwork or lrwork or liwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit. \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  The algorithm failed to compute an eigenvalue while
        working on the submatrix lying in rows and columns
        info/(N+1) through mod(info,N+1).
    *  */
template< typename Ta, typename Tb >
int stedc( char* compz, int* n, Tb* d, Tb* e, Ta* z, int* ldz, Ta* work, int* lwork, Tb*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
{
  return stedc( compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info );
}

/*! @brief Tridiagonal MRRR algorithm
    *

    * @details
    * \b Purpose:
    * \verbatim
        Tridiagonal MRRR algorithm.
        Computation of selected eigenvalues and, optionally, eigenvectors of a real symmetric
        tridiagonal matrix T. Any such unreduced matrix has a well defined set of pairwise
        different real eigenvalues, the corresponding real eigenvectors are pairwise orthogonal.

        The spectrum may be computed either completely or partially by specifying either an
        interval (vl,vu] or a range of indices il:iu for the desired eigenvalues.

        Depending on the number of desired eigenvalues, these are computed either by bisection or
        the dqds algorithm. Numerically orthogonal eigenvectors are computed by the use of various
        suitable L D L^T factorizations near clusters of close eigenvalues (referred to as RRRs,
        Relatively Robust Representations). An informal sketch of the algorithm follows.

        For each unreduced block (submatrix) of T,
        (a) Compute T - sigma I  = L D L^T, so that L and D
        define all the wanted eigenvalues to high relative accuracy.
        This means that small relative changes in the entries of D and L
        cause only small relative changes in the eigenvalues and
        eigenvectors. The standard (unfactored) representation of the
        tridiagonal matrix T does not have this property in general.
        (b) Compute the eigenvalues to suitable accuracy.
        If the eigenvectors are desired, the algorithm attains full
        accuracy of the computed eigenvalues only right before
        the corresponding vectors have to be computed, see steps c) and d).
        (c) For each cluster of close eigenvalues, select a new
        shift close to the cluster, find a new factorization, and refine
        the shifted eigenvalues to suitable accuracy.
        (d) For each eigenvalue with a large enough relative separation compute
        the corresponding eigenvector by forming a rank revealing twisted
        factorization. Go back to (c) for any clusters that remain.

        For more details, see:
        - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations to compute
        orthogonal eigenvectors of symmetric tridiagonal matrices," Linear Algebra and its
        Applications, 387(1), pp. 1-28, August 2004.
        - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and Relative Gaps,"
        SIAM Journal on Matrix Analysis and Applications, Vol. 25,   2004.
        Also LAPACK Working Note 154.
        - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric tridiagonal
        eigenvalue/eigenvector problem", Computer Science Division Technical Report No.
        UCB/CSD-97-971, UC Berkeley, May 1997.

        Further Details
        1.SSTEMR works only on machines which follow IEEE-754 floating-point standard in their
        handling of infinities and NaNs. This permits the use of efficient inner loops avoiding a
        check for zero divisors.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only; \n
    = 'V':  Compute eigenvalues and eigenvectors.

    * @param[in] range
    range is char* \n
    = 'A': all eigenvalues will be found. \n
    = 'V': all eigenvalues in the half-open interval (vl,vu]
    will be found. \n
    = 'I': the il-th through iu-th eigenvalues will be found.

    * @param[in] n
    n is int* \n
    The order of the matrix.  n >= 0.

    * @param[in,out] d
    d is float/double array, dimension (n) \n
    On entry, the N diagonal elements of the tridiagonal matrix
    T. \n
    On exit, D is overwritten.

    * @param[in,out] e
    e is float/double array, dimension (n) \n
    On entry, the (n-1) subdiagonal elements of the tridiagonal
    matrix T in elements 1 to n-1 of E. E(n) need not be set on
    input, but is used internally as workspace. \n
    On exit, E is overwritten.

    * @param[in] vl
    vl is int* \n
    If range='V', the lower bound of the interval to
    be searched for eigenvalues. vl < vu. \n
    Not referenced if range = 'A' or 'I'.

    * @param[in] vu
    vu is T* \n
    If range='V', the upper bound of the interval to
    be searched for eigenvalues. vl < vu. \n
    Not referenced if range = 'A' or 'I'.

    * @param[in] il
    il is int* \n
    If range='I', the index of the
    smallest eigenvalue to be returned. \n
    1 <= il <= iu <= N, if N > 0. \n
    Not referenced if range = 'A' or 'V'.

    * @param[in] iu
    iu is int* \n
    If range='I', the index of the
    largest eigenvalue to be returned. \n
    1 <= il <= iu <= N, if N > 0. \n
    Not referenced if range = 'A' or 'V'.

    * @param[out] m
    m is int* \n
    The total number of eigenvalues found.  0 <= M <= N. \n
    If range = 'A', M = N, and if range = 'I', M = iu-il+1.

    * @param[out] w
    w is float/double array, dimension (n) \n
    The first M elements contain the selected eigenvalues in
    ascending order.

    * @param[out] z
    z is float/double array, dimension (ldz, max(1,m) ) \n
    If jobz = 'V', and if info = 0, then the first M columns of Z
    contain the orthonormal eigenvectors of the matrix T
    corresponding to the selected eigenvalues, with the i-th
    column of Z holding the eigenvector associated with W(i).
    If jobz = 'N', then Z is not referenced. \n
    Note: the user must ensure that at least max(1,m) columns are
    supplied in the array Z; if range = 'V', the exact value of M
    is not known in advance and can be computed with a workspace
    query by setting nzc = -1, see below.

    * @param[in] ldz
    ldz is int* \n
    The leading dimension of the array Z.  ldz >= 1, and if
    jobz = 'V', then ldz >= max(1,n).

    * @param[in] nzc
    nzc is int* \n
    The number of eigenvectors to be held in the array Z. \n
    If range = 'A', then nzc >= max(1,n). \n
    If range = 'V', then nzc >= the number of eigenvalues in (vl,vu]. \n
    If range = 'I', then nzc >= iu-il+1. \n
    If nzc = -1, then a workspace query is assumed; the
    routine calculates the number of columns of the array Z that
    are needed to hold the eigenvectors. \n
    This value is returned as the first entry of the Z array, and
    no error message related to nzc is issued by XERBLA.

    * @param[out] isuppz
    isuppz is int array, dimension ( 2*max(1,m) ) \n
    The support of the eigenvectors in Z, i.e., the indices
    indicating the nonzero elements in Z. The i-th computed eigenvector
    is nonzero only in elements isuppz( 2*i-1 ) through
    isuppz( 2*i ). This is relevant in the case when the matrix
    is split. isuppz is only accessed when jobz is 'V' and N > 0.

    * @param[in,out] tryrac
    tryrac is LOGICAL \n
    If tryrac.EQ..TRUE., indicates that the code should check whether
    the tridiagonal matrix defines its eigenvalues to high relative
    accuracy.  If so, the code uses relative-accuracy preserving
    algorithms that might be (a bit) slower depending on the matrix.
    If the matrix does not define its eigenvalues to high relative
    accuracy, the code can uses possibly faster algorithms.
    If tryrac.EQ..FALSE., the code is not required to guarantee
    relatively accurate eigenvalues and can use the fastest possible
    techniques. \n
    On exit, a .TRUE. tryrac will be set to .FALSE. if the matrix
    does not define its eigenvalues to high relative accuracy.

    * @param[out] work
    work is float/double array, dimension (lwork) \n
    On exit, if info = 0, work(1) returns the optimal
    (and minimal) lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,18*N) \n
    if jobz = 'V', and lwork >= max(1,12*N) if jobz = 'N'. \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] iwork
    iwork is int array, dimension (liwork) \n
    On exit, if info = 0, iwork(1) returns the optimal liwork.

    * @param[in] liwork
    liwork is int* \n
    The dimension of the array iwork.  liwork >= max(1,10*N)
    if the eigenvectors are desired, and liwork >= max(1,8*N)
    if only the eigenvalues are to be computed.
    If liwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal size of the iwork array,
    returns this value as the first entry of the iwork array, and
    no error message related to liwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    On exit, info \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  if info = 1X, internal error in SLARRE,
    if info = 2X, internal error in SLARRV.
    Here, the digit X = ABS( IINFO ) < 10, where IINFO is
    the nonzero error code returned by SLARRE or
    SLARRV, respectively.
    *  */
template< typename T >
int stemr( char* jobz, char* range, int* n, T*  d, T*  e, int* vl, int* vu, int* il, int* iu, int* m, T*  w, T* z, int* ldz, int* nzc, int* isuppz, int* tryrac, T*  work, int* lwork, int* iwork, int* liwork, int* info )
{
  return stemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
}

/*! @brief Tridiagonal MRRR algorithm
    *

    * @details
    * \b Purpose:
    * \verbatim
        Tridiagonal MRRR algorithm.
        Computation of selected eigenvalues and, optionally, eigenvectors of a real symmetric
        tridiagonal matrix T. Any such unreduced matrix has a well defined set of pairwise
        different real eigenvalues, the corresponding real eigenvectors are pairwise orthogonal.

        The spectrum may be computed either completely or partially by specifying either an
        interval (vl,vu] or a range of indices il:iu for the desired eigenvalues.

        Depending on the number of desired eigenvalues, these are computed either by bisection or
        the dqds algorithm. Numerically orthogonal eigenvectors are computed by the use of various
        suitable L D L^T factorizations near clusters of close eigenvalues (referred to as RRRs,
        Relatively Robust Representations). An informal sketch of the algorithm follows.

        For each unreduced block (submatrix) of T,
            (a) Compute T - sigma I  = L D L^T, so that L and D
                define all the wanted eigenvalues to high relative accuracy.
                This means that small relative changes in the entries of D and L
                cause only small relative changes in the eigenvalues and
                eigenvectors. The standard (unfactored) representation of the
                tridiagonal matrix T does not have this property in general.
            (b) Compute the eigenvalues to suitable accuracy.
                If the eigenvectors are desired, the algorithm attains full
                accuracy of the computed eigenvalues only right before
                the corresponding vectors have to be computed, see steps c) and d).
            (c) For each cluster of close eigenvalues, select a new
                shift close to the cluster, find a new factorization, and refine
                the shifted eigenvalues to suitable accuracy.
            (d) For each eigenvalue with a large enough relative separation compute
                the corresponding eigenvector by forming a rank revealing twisted
                factorization. Go back to (c) for any clusters that remain.

        For more details, see:
        - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations to compute
        orthogonal eigenvectors of symmetric tridiagonal matrices," Linear Algebra and its
        Applications, 387(1), pp. 1-28, August 2004.
        - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and Relative Gaps,"
        SIAM Journal on Matrix Analysis and Applications, Vol. 25,   2004.
        Also LAPACK Working Note 154.
        - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric tridiagonal
        eigenvalue/eigenvector problem", Computer Science Division Technical Report No.
        UCB/CSD-97-971, UC Berkeley, May 1997.

        Further Details
        1. CSTEMR works only on machines which follow IEEE-754 floating-point standard in their
        handling of infinities and NaNs. This permits the use of efficient inner loops avoiding a
        check for zero divisors.

        2. LAPACK routines can be used to reduce a complex Hermitean matrix to real symmetric
        tridiagonal form.

        (Any complex Hermitean tridiagonal matrix has real values on its diagonal and potentially
        complex numbers on its off-diagonals. By applying a similarity transform with an
        appropriate diagonal matrix diag(1,e^{i \phy_1}, ... , e^{i \phy_{n-1}}), the complex
        Hermitean matrix can be transformed into a real symmetric matrix and complex arithmetic can
        be entirely avoided.)

        While the eigenvectors of the real symmetric tridiagonal matrix are real, the eigenvectors
        of original complex Hermitean matrix have complex entries in general. Since LAPACK drivers
        overwrite the matrix data with the eigenvectors, CSTEMR accepts complex workspace to
        facilitate interoperability with CUNMTR or CUPMTR.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only; \n
    = 'V':  Compute eigenvalues and eigenvectors.

    * @param[in] range
    range is char* \n
    = 'A': all eigenvalues will be found. \n
    = 'V': all eigenvalues in the half-open interval (vl,vu]
        will be found. \n
    = 'I': the il-th through iu-th eigenvalues will be found.

    * @param[in] n
    n is int* \n
    The order of the matrix.  n >= 0.

    * @param[in,out] d
    d is float/double array, dimension (n) \n
    On entry, the N diagonal elements of the tridiagonal matrix
    T. \n
    On exit, D is overwritten.

    * @param[in,out] e
    e is float/double array, dimension (n) \n
    On entry, the (N-1) subdiagonal elements of the tridiagonal
    matrix T in elements 1 to n-1 of E. E(n) need not be set on
    input, but is used internally as workspace. \n
    On exit, E is overwritten.

    * @param[in] vl
    vl is int* \n
    If range='V', the lower bound of the interval to
    be searched for eigenvalues. vl < vu. \n
    Not referenced if range = 'A' or 'I'.

    * @param[in] vu
    vu is int* \n
    If range='V', the upper bound of the interval to
    be searched for eigenvalues. vl < vu. \n
    Not referenced if range = 'A' or 'I'.

    * @param[in] il
    il is int* \n
    If range='I', the index of the
    smallest eigenvalue to be returned. \n
    1 <= il <= iu <= N, if N > 0. \n
    Not referenced if range = 'A' or 'V'.

    * @param[in] iu
    iu is int* \n
    If range='I', the index of the
    largest eigenvalue to be returned. \n
    1 <= il <= iu <= N, if N > 0. \n
    Not referenced if range = 'A' or 'V'.

    * @param[out] m
    m is int* \n
    The total number of eigenvalues found.  0 <= M <= N. \n
    If range = 'A', M = N, and if range = 'I', M = iu-il+1.

    * @param[out] w
    w is float/double array, dimension (n) \n
    The first M elements contain the selected eigenvalues in
    ascending order.

    * @param[out] z
    z is scomplex/dcomplex array, dimension (ldz, max(1,m) ) \n
    If jobz = 'V', and if info = 0, then the first M columns of Z
    contain the orthonormal eigenvectors of the matrix T
    corresponding to the selected eigenvalues, with the i-th
    column of Z holding the eigenvector associated with W(i). \n
    If jobz = 'N', then Z is not referenced. \n
    Note: the user must ensure that at least max(1,m) columns are
    supplied in the array Z; if range = 'V', the exact value of M
    is not known in advance and can be computed with a workspace
    query by setting nzc = -1, see below.

    * @param[in] ldz
    ldz is int* \n
    The leading dimension of the array Z.  ldz >= 1, and if
    jobz = 'V', then ldz >= max(1,n).

    * @param[in] nzc
    nzc is int* \n
    The number of eigenvectors to be held in the array Z. \n
    If range = 'A', then nzc >= max(1,n). \n
    If range = 'V', then nzc >= the number of eigenvalues in (vl,vu]. \n
    If range = 'I', then nzc >= iu-il+1. \n
    If nzc = -1, then a workspace query is assumed; the
    routine calculates the number of columns of the array Z that
    are needed to hold the eigenvectors. \n
    This value is returned as the first entry of the Z array, and
    no error message related to nzc is issued by XERBLA.

    * @param[out] isuppz
    isuppz is int array, dimension ( 2*max(1,m) ) \n
    The support of the eigenvectors in Z, i.e., the indices
    indicating the nonzero elements in Z. The i-th computed eigenvector
    is nonzero only in elements isuppz( 2*i-1 ) through
    isuppz( 2*i ). This is relevant in the case when the matrix
    is split. isuppz is only accessed when jobz is 'V' and N > 0.

    * @param[in,out] tryrac
    tryrac is LOGICAL \n
    If tryrac.EQ..TRUE., indicates that the code should check whether
    the tridiagonal matrix defines its eigenvalues to high relative
    accuracy.  If so, the code uses relative-accuracy preserving
    algorithms that might be (a bit) slower depending on the matrix.
    If the matrix does not define its eigenvalues to high relative
    accuracy, the code can uses possibly faster algorithms.
    If tryrac.EQ..FALSE., the code is not required to guarantee
    relatively accurate eigenvalues and can use the fastest possible
    techniques. \n
    On exit, a .TRUE. tryrac will be set to .FALSE. if the matrix
    does not define its eigenvalues to high relative accuracy.

    * @param[out] work
    work is float/double array, dimension (lwork) \n
    On exit, if info = 0, work(1) returns the optimal
    (and minimal) lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= max(1,18*N) \n
    if jobz = 'V', and lwork >= max(1,12*N) if jobz = 'N'. \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] iwork
    iwork is int array, dimension (liwork) \n
    On exit, if info = 0, iwork(1) returns the optimal liwork.

    * @param[in] liwork
    liwork is int* \n
    The dimension of the array iwork.  liwork >= max(1,10*N)
    if the eigenvectors are desired, and liwork >= max(1,8*N)
    if only the eigenvalues are to be computed.
    If liwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal size of the iwork array,
    returns this value as the first entry of the iwork array, and
    no error message related to liwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    On exit, info \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  if info = 1X, internal error in SLARRE,
        if info = 2X, internal error in CLARRV.
        Here, the digit X = ABS( IINFO ) < 10, where IINFO is
        the nonzero error code returned by SLARRE or
        CLARRV, respectively.
    *  */
template< typename Ta, typename Tb >
int stemr( char* jobz, char* range, int* n, Tb*  d, Tb*  e, int* vl, int* vu, int* il, int* iu, int* m, Tb*  w, Ta* z, int* ldz, int* nzc, int* isuppz, int* tryrac, Tb*  work, int* lwork, int* iwork, int* liwork, int* info )
{
  return stemr( jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info );
}

/*! @brief Eigenvalue decomposition (QR algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Eigenvalue decomposition (QR algorithm).
        Computation of all eigenvalues and, optionally, eigenvectors of a real symmetric matrix a.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only; \n
    = 'V':  Compute eigenvalues and eigenvectors.

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda, N) \n
    On entry, the symmetric matrix a.  If uplo = 'U', the
    leading n-by-n upper triangular part of a contains the
    upper triangular part of the matrix a.  If uplo = 'L',
    the leading n-by-n lower triangular part of a contains
    the lower triangular part of the matrix a. \n
    On exit, if jobz = 'V', then if info = 0, A contains the
    orthonormal eigenvectors of the matrix a. \n
    If jobz = 'N', then on exit the lower triangle (if uplo='L')
    or the upper triangle (if uplo='U') of A, including the
    diagonal, is destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] w
    w is float/double array, dimension (n) \n
    If info = 0, the eigenvalues in ascending order.

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The length of the array work.  lwork >= max(1,3*N-1). \n
    For optimal efficiency, lwork >= (NB+2)*N, \n
    where NB is the blocksize for SSYTRD returned by ILAENV. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] rwork
    rwork is float/double array, dimension (max(1, 3*N-2))

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  if info = i, the algorithm failed to converge; i
    off-diagonal elements of an intermediate tridiagonal
    form did not converge to zero.
    *  */
template< typename T >
int syev( char* jobz, char* uplo, int* n, T* a, int* lda, T*  w, T* work, int* lwork, T*  rwork, int* info )
{
  return syev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info );
}

/*! @brief Eigenvalue decomposition (QR algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Eigenvalue decomposition (QR algorithm)
        Computation of all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix a.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only; \n
    = 'V':  Compute eigenvalues and eigenvectors.

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda, N) \n
    On entry, the Hermitian matrix a.  If uplo = 'U', the
    leading n-by-n upper triangular part of a contains the
    upper triangular part of the matrix a.  If uplo = 'L',
    the leading n-by-n lower triangular part of a contains
    the lower triangular part of the matrix a. \n
    On exit, if jobz = 'V', then if info = 0, A contains the
    orthonormal eigenvectors of the matrix a. \n
    If jobz = 'N', then on exit the lower triangle (if uplo='L')
    or the upper triangle (if uplo='U') of A, including the
    diagonal, is destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] w
    w is float/double array, dimension (n) \n
    If info = 0, the eigenvalues in ascending order.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The length of the array work.  lwork >= max(1,2*n-1). \n
    For optimal efficiency, lwork >= (NB+1)*n, \n
    where NB is the blocksize for CHETRD returned by ILAENV. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] rwork
    rwork is float/double array, dimension (max(1, 3*N-2)) \n

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  if info = i, the algorithm failed to converge; i
    off-diagonal elements of an intermediate tridiagonal
    form did not converge to zero.
    *  */
template< typename Ta, typename Tb >
int heev( char* jobz, char* uplo, int* n, Ta* a, int* lda, Tb*  w, Ta* work, int* lwork, Tb*  rwork, int* info )
{
  return heev( jobz, uplo, n, a, lda, w, work, lwork, rwork, info );
}

/*! @brief Eigenvalue decomposition (divide-and-conquer)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of the eigenvalues and, optionally, the left and/or right eigenvectors for SY
        matrices. SSYEVD computes all eigenvalues and, optionally, eigenvectors of a real
        symmetric matrix a. If eigenvectors are desired, it uses a divide and conquer algorithm.

        The divide and conquer algorithm makes very mild assumptions about floating point
        arithmetic. It will work on machines with a guard digit in add/subtract, or on those
        binary machines without guard digits which subtract like the Cray X-MP, Cray Y-MP, Cray
        C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without
        guard digits, but we know of none.

        Because of large use of BLAS of level 3, SSYEVD needs N**2 more workspace than SSYEVX.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only; \n
    = 'V':  Compute eigenvalues and eigenvectors.

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda, N) \n
    On entry, the symmetric matrix a.  If uplo = 'U', the
    leading n-by-n upper triangular part of a contains the
    upper triangular part of the matrix a.  If uplo = 'L',
    the leading n-by-n lower triangular part of a contains
    the lower triangular part of the matrix a. \n
    On exit, if jobz = 'V', then if info = 0, A contains the
    orthonormal eigenvectors of the matrix a. \n
    If jobz = 'N', then on exit the lower triangle (if uplo='L')
    or the upper triangle (if uplo='U') of A, including the
    diagonal, is destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] w
    w is float/double array, dimension (n) \n
    If info = 0, the eigenvalues in ascending order.

    * @param[out] work
    work is float/double array, dimension (lwork) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    If N <= 1,   lwork must be at least 1. \n
    If jobz = 'N' and N > 1, lwork must be at least 2*n+1. \n
    If jobz = 'V' and N > 1, lwork must be at least
    1 + 6*N + 2*n**2. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal sizes of the work and iwork
    arrays, returns these values as the first entries of the work
    and iwork arrays, and no error message related to lwork or
    liwork is issued by XERBLA.

    * @param[out] iwork
    iwork is int array, dimension (max(1,liwork)) \n
    On exit, if info = 0, iwork(1) returns the optimal liwork.

    * @param[in] liwork
    liwork is int* \n
    The dimension of the array iwork. \n
    If N <= 1, liwork must be at least 1. \n
    If jobz  = 'N' and N > 1, liwork must be at least 1. \n
    If jobz  = 'V' and N > 1, liwork must be at least 3 + 5*N. \n \n
    If liwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal sizes of the work and
    iwork arrays, returns these values as the first entries of
    the work and iwork arrays, and no error message related to
    lwork or liwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  if info = i and jobz = 'N', then the algorithm failed
    to converge; i off-diagonal elements of an intermediate
    tridiagonal form did not converge to zero;
    if info = i and jobz = 'V', then the algorithm failed
    to compute an eigenvalue while working on the submatrix
    lying in rows and columns info/(N+1) through
    mod(info,N+1).
    *  */
template< typename T >
int syevd( char* jobz, char* uplo, int* n, T* a, int* lda, T*  w, T* work, int* lwork, int* iwork, int* liwork, int* info )
{
  return syevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info );
}

/*! @brief Eigenvalue decomposition (divide-and-conquer)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Eigenvalue decomposition (divide-and-conquer).
        Computation of the eigenvalues and, optionally, the left and/or right eigenvectors for
        HE matrices

        CHEEVD computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian
        matrix a.  If eigenvectors are desired, it uses a divide and conquer algorithm.

        The divide and conquer algorithm makes very mild assumptions about floating point
        arithmetic. It will work on machines with a guard digit in add/subtract, or on those
        binary machines without guard digits which subtract like the Cray X-MP, Cray Y-MP, Cray
        C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without
        guard digits, but we know of none.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only; \n
    = 'V':  Compute eigenvalues and eigenvectors.

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda, N) \n
    On entry, the Hermitian matrix a.  If uplo = 'U', the
    leading n-by-n upper triangular part of a contains the
    upper triangular part of the matrix a.  If uplo = 'L',
    the leading n-by-n lower triangular part of a contains
    the lower triangular part of the matrix a. \n
    On exit, if jobz = 'V', then if info = 0, A contains the
    orthonormal eigenvectors of the matrix a. \n
    If jobz = 'N', then on exit the lower triangle (if uplo='L')
    or the upper triangle (if uplo='U') of A, including the
    diagonal, is destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[out] w
    w is float/double array, dimension (n) \n
    If info = 0, the eigenvalues in ascending order.

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The length of the array work. \n
    If N <= 1, lwork must be at least 1. \n
    If jobz  = 'N' and N > 1, lwork must be at least N + 1. \n
    If jobz  = 'V' and N > 1, lwork must be at least 2*n + N**2. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal sizes of the work, rwork and
    iwork arrays, returns these values as the first entries of
    the work, rwork and iwork arrays, and no error message
    related to lwork or lrwork or liwork is issued by XERBLA.

    * @param[out] rwork
    rwork is float/double array, dimension (lrwork) \n
    On exit, if info = 0, rwork(1) returns the optimal lrwork.

    * @param[in] lrwork
    lrwork is int* \n
    The dimension of the array rwork. \n
    If N <= 1, lrwork must be at least 1. \n
    If jobz  = 'N' and N > 1, lrwork must be at least N. \n
    If jobz  = 'V' and N > 1, lrwork must be at least
    1 + 5*N + 2*n**2. \n \n
    If lrwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal sizes of the work, rwork
    and iwork arrays, returns these values as the first entries
    of the work, rwork and iwork arrays, and no error message
    related to lwork or lrwork or liwork is issued by XERBLA.

    * @param[out] iwork
    iwork is int array, dimension (max(1,liwork)) \n
    On exit, if info = 0, iwork(1) returns the optimal liwork.

    * @param[in] liwork
    liwork is int* \n
    The dimension of the array iwork. \n
    If N <= 1, liwork must be at least 1. \n
    If jobz  = 'N' and N > 1, liwork must be at least 1. \n
    If jobz  = 'V' and N > 1, liwork must be at least 3 + 5*N. \n \n
    If liwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal sizes of the work, rwork
    and iwork arrays, returns these values as the first entries
    of the work, rwork and iwork arrays, and no error message
    related to lwork or lrwork or liwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  if info = i and jobz = 'N', then the algorithm failed
    to converge; i off-diagonal elements of an intermediate
    tridiagonal form did not converge to zero;
    if info = i and jobz = 'V', then the algorithm failed
    to compute an eigenvalue while working on the submatrix
    lying in rows and columns info/(N+1) through
    mod(info,N+1).
    *  */
template< typename Ta, typename Tb >
int heevd( char* jobz, char* uplo, int* n, Ta* a, int* lda, Tb*  w, Ta* work, int* lwork, Tb*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
{
  return heevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info );
}

/*! @brief Hermitian eigenvalue decomposition (MRRR)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Hermitian eigenvalue decomposition (MRRR).
        Computation of eigenvalues and, optionally, the left and/or right eigenvectors for SY
        matrices

        SSYEVR computes selected eigenvalues and, optionally, eigenvectors of a real symmetric
        matrix a. Eigenvalues and eigenvectors can be selected by specifying either a range of
        values or a range of indices for the desired eigenvalues.

        SSYEVR first reduces the matrix a to tridiagonal form T with a call to SSYTRD.  Then,
        whenever possible, SSYEVR calls SSTEMR to compute the eigenspectrum using Relatively
        Robust Representations. SSTEMR computes eigenvalues by the dqds algorithm, while orthogonal
        eigenvectors are computed from various "good" L D L^T representations (also known as
        Relatively Robust Representations). Gram-Schmidt orthogonalization is avoided as far as
        possible. More specifically, the various steps of the algorithm are as follows.

        For each unreduced block (submatrix) of T,
        (a) Compute T - sigma I  = L D L^T, so that L and D
        define all the wanted eigenvalues to high relative accuracy.
        This means that small relative changes in the entries of D and L
        cause only small relative changes in the eigenvalues and
        eigenvectors. The standard (unfactored) representation of the
        tridiagonal matrix T does not have this property in general.
        (b) Compute the eigenvalues to suitable accuracy.
        If the eigenvectors are desired, the algorithm attains full
        accuracy of the computed eigenvalues only right before
        the corresponding vectors have to be computed, see steps c) and d).
        (c) For each cluster of close eigenvalues, select a new
        shift close to the cluster, find a new factorization, and refine
        the shifted eigenvalues to suitable accuracy.
        (d) For each eigenvalue with a large enough relative separation compute
        the corresponding eigenvector by forming a rank revealing twisted
        factorization. Go back to (c) for any clusters that remain.

        The desired accuracy of the output can be specified by the input parameter abstol.

        For more details, see SSTEMR's documentation and:
        - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations to compute
        orthogonal eigenvectors of symmetric tridiagonal matrices," Linear Algebra and its
        Applications, 387(1), pp. 1-28, August 2004.
        - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and   Relative Gaps,"
        SIAM Journal on Matrix Analysis and Applications, Vol. 25,   2004.
        Also LAPACK Working Note 154.
        - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric tridiagonal
        eigenvalue/eigenvector problem",   Computer Science Division Technical Report No.
        UCB/CSD-97-971, UC Berkeley, May 1997.


        Note 1 : SSYEVR calls SSTEMR when the full spectrum is requested on machines which conform
        to the ieee-754 floating point standard. SSYEVR calls SSTEBZ and SSTEIN on non-ieee
        machines and when partial spectrum requests are made.

        Normal execution of SSTEMR may create NaNs and infinities and hence may abort due to a
        floating point exception in environments which do not handle NaNs and infinities in the
        ieee standard default manner.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only; \n
    = 'V':  Compute eigenvalues and eigenvectors.

    * @param[in] range
    range is char* \n
    = 'A': all eigenvalues will be found. \n
    = 'V': all eigenvalues in the half-open interval (vl,vu]
    will be found. \n
    = 'I': the il-th through iu-th eigenvalues will be found. \n
    For range = 'V' or 'I' and iu - il < N - 1, SSTEBZ and
    SSTEIN are called

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda, N) \n
    On entry, the symmetric matrix a.  If uplo = 'U', the
    leading n-by-n upper triangular part of a contains the
    upper triangular part of the matrix a.  If uplo = 'L',
    the leading n-by-n lower triangular part of a contains
    the lower triangular part of the matrix a. \n
    On exit, the lower triangle (if uplo='L') or the upper
    triangle (if uplo='U') of A, including the diagonal, is
    destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[in] vl
    vl is float/double* \n
    If range='V', the lower bound of the interval to
    be searched for eigenvalues. vl < vu. \n
    Not referenced if range = 'A' or 'I'.

    * @param[in] vu
    vu is float/double* \n
    If range='V', the upper bound of the interval to
    be searched for eigenvalues. vl < vu. \n
    Not referenced if range = 'A' or 'I'.

    * @param[in] il
    il is int* \n
    If range='I', the index of the
    smallest eigenvalue to be returned. \n
    1 <= il <= iu <= N, if N > 0; il = 1 and iu = 0 if N = 0. \n
    Not referenced if range = 'A' or 'V'.

    * @param[in] iu
    iu is int* \n
    If range='I', the index of the
    largest eigenvalue to be returned. \n
    1 <= il <= iu <= N, if N > 0; il = 1 and iu = 0 if N = 0. \n
    Not referenced if range = 'A' or 'V'.

    * @param[in] abstol
    abstol is float/double* \n
    The absolute error tolerance for the eigenvalues. \n
    An approximate eigenvalue is accepted as converged
    when it is determined to lie in an interval [a,b]
    of width less than or equal to \n \n
    abstol + EPS *   max( |a|,|b| ) , \n \n
    where EPS is the machine precision.  If abstol is less than
    or equal to zero, then  EPS*|T|  will be used in its place,
    where |T| is the 1-norm of the tridiagonal matrix obtained
    by reducing A to tridiagonal form. \n \n
    See "Computing Small Singular Values of Bidiagonal Matrices
    with Guaranteed High Relative Accuracy," by Demmel and
    Kahan, LAPACK Working Note #3. \n \n
    If high relative accuracy is important, set abstol to
    SLAMCH( 'Safe minimum' ).  Doing so will guarantee that
    eigenvalues are computed to high relative accuracy when
    possible in future releases.  The current code does not
    make any guarantees about high relative accuracy, but
    future releases will. See J. Barlow and J. Demmel,
    "Computing Accurate Eigensystems of Scaled Diagonally
    Dominant Matrices", LAPACK Working Note #7, for a discussion
    of which matrices define their eigenvalues to high relative
    accuracy.

    * @param[out] m
    m is int* \n
    The total number of eigenvalues found.  0 <= M <= N. \n
    If range = 'A', M = N, and if range = 'I', M = iu-il+1.

    * @param[out] w
    w is float/double/scomplex/dcomplex array, dimension (n) \n
    The first M elements contain the selected eigenvalues in
    ascending order.

    * @param[out] z
    z is float/double array, dimension (ldz, max(1,m)) \n
    If jobz = 'V', then if info = 0, the first M columns of Z
    contain the orthonormal eigenvectors of the matrix a
    corresponding to the selected eigenvalues, with the i-th
    column of Z holding the eigenvector associated with W(i). \n
    If jobz = 'N', then Z is not referenced. \n
    Note: the user must ensure that at least max(1,m) columns are
    supplied in the array Z; if range = 'V', the exact value of M
    is not known in advance and an upper bound must be used.
    Supplying N columns is always safe.

    * @param[in] ldz
    ldz is int* \n
    The leading dimension of the array Z.  ldz >= 1, and if
    jobz = 'V', ldz >= max(1,n).

    * @param[out] isuppz
    isuppz is int array, dimension ( 2*max(1,m) ) \n
    The support of the eigenvectors in Z, i.e., the indices
    indicating the nonzero elements in Z. The i-th eigenvector
    is nonzero only in elements isuppz( 2*i-1 ) through
    isuppz( 2*i ). This is an output of SSTEMR (tridiagonal
    matrix). The support of the eigenvectors of A is typically
    1:N because of the orthogonal transformations applied by SORMTR. \n
    Implemented only for range = 'A' or 'I' and iu - il = N - 1

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work.  lwork >= max(1,26*N). \n
    For optimal efficiency, lwork >= (NB+6)*N, \n
    where NB is the max of the blocksize for SSYTRD and SORMTR
    returned by ILAENV. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal sizes of the work and iwork
    arrays, returns these values as the first entries of the work
    and iwork arrays, and no error message related to lwork or
    liwork is issued by XERBLA.

    * @param[out] iwork
    iwork is int array, dimension (max(1,liwork)) \n
    On exit, if info = 0, iwork(1) returns the optimal lwork.

    * @param[in] liwork
    liwork is int* \n
    The dimension of the array iwork.  liwork >= max(1,10*N). \n \n
    If liwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal sizes of the work and
    iwork arrays, returns these values as the first entries of
    the work and iwork arrays, and no error message related to
    lwork or liwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  Internal error
    *  */
template< typename T >
int syevr( char* jobz, char* range, char* uplo, int* n, T* a, int* lda, T* vl, T*  vu, int* il, int* iu, T*  abstol, int* m, T*  w, T* z, int* ldz, int* isuppz, T* work, int* lwork, int* iwork, int* liwork, int* info )
{
  return syevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info );
}

/*! @brief Hermitian eigenvalue decomposition (MRRR)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of eigenvalues and, optionally, the left and/or right eigenvectors for HE
        matrices. CHEEVR computes selected eigenvalues and, optionally, eigenvectors of a complex
        Hermitian matrix a. Eigenvalues and eigenvectors can be selected by specifying either a
        range of values or a range of indices for the desired eigenvalues.

        CHEEVR first reduces the matrix a to tridiagonal form T with a call to CHETRD.  Then,
        whenever possible, CHEEVR calls CSTEMR to compute the eigenspectrum using Relatively
        Robust Representations. CSTEMR computes eigenvalues by the dqds algorithm, while
        orthogonal eigenvectors are computed from various "good" L D L^T representations (also
        known as Relatively Robust Representations). Gram-Schmidt orthogonalization is avoided as
        far as possible. More specifically, the various steps of the algorithm are as follows.

        For each unreduced block (submatrix) of T,
        (a) Compute T - sigma I  = L D L^T, so that L and D
        define all the wanted eigenvalues to high relative accuracy.
        This means that small relative changes in the entries of D and L
        cause only small relative changes in the eigenvalues and
        eigenvectors. The standard (unfactored) representation of the
        tridiagonal matrix T does not have this property in general.
        (b) Compute the eigenvalues to suitable accuracy.
        If the eigenvectors are desired, the algorithm attains full
        accuracy of the computed eigenvalues only right before
        the corresponding vectors have to be computed, see steps c) and d).
        (c) For each cluster of close eigenvalues, select a new
        shift close to the cluster, find a new factorization, and refine
        the shifted eigenvalues to suitable accuracy.
        (d) For each eigenvalue with a large enough relative separation compute
        the corresponding eigenvector by forming a rank revealing twisted
        factorization. Go back to (c) for any clusters that remain.

        The desired accuracy of the output can be specified by the input parameter abstol.

        For more details, see DSTEMR's documentation and:
        - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations to compute
        orthogonal eigenvectors of symmetric tridiagonal matrices," Linear Algebra and its
        Applications, 387(1), pp. 1-28, August 2004.
        - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and   Relative Gaps,"
        SIAM Journal on Matrix Analysis and Applications, Vol. 25,   2004.
        Also LAPACK Working Note 154.
        - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric tridiagonal
        eigenvalue/eigenvector problem",   Computer Science Division Technical Report No.
        UCB/CSD-97-971, UC Berkeley, May 1997.


        Note 1 : CHEEVR calls CSTEMR when the full spectrum is requested on machines which conform
        to the ieee-754 floating point standard. CHEEVR calls SSTEBZ and CSTEIN on non-ieee
        machines and when partial spectrum requests are made.

        Normal execution of CSTEMR may create NaNs and infinities and hence may abort due to a
        floating point exception in environments which do not handle NaNs and infinities in the
        ieee standard default manner.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    = 'N':  Compute eigenvalues only; \n
    = 'V':  Compute eigenvalues and eigenvectors.

    * @param[in] range
    range is char* \n
    = 'A': all eigenvalues will be found. \n
    = 'V': all eigenvalues in the half-open interval (vl,vu]
    will be found. \n
    = 'I': the il-th through iu-th eigenvalues will be found. \n
    For range = 'V' or 'I' and iu - il < N - 1, SSTEBZ and
    CSTEIN are called

    * @param[in] uplo
    uplo is char* \n
    = 'U':  Upper triangle of A is stored; \n
    = 'L':  Lower triangle of A is stored.

    * @param[in] n
    n is int* \n
    The order of the matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda, N) \n
    On entry, the Hermitian matrix a.  If uplo = 'U', the
    leading n-by-n upper triangular part of a contains the
    upper triangular part of the matrix a.  If uplo = 'L',
    the leading n-by-n lower triangular part of a contains
    the lower triangular part of the matrix a. \n
    On exit, the lower triangle (if uplo='L') or the upper
    triangle (if uplo='U') of A, including the diagonal, is
    destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,n).

    * @param[in] vl
    vl is float/double* \n
    If range='V', the lower bound of the interval to
    be searched for eigenvalues. vl < vu. \n
    Not referenced if range = 'A' or 'I'.

    * @param[in] vu
    vu is float/double* \n
    If range='V', the upper bound of the interval to
    be searched for eigenvalues. vl < vu. \n
    Not referenced if range = 'A' or 'I'.

    * @param[in] il
    il is int* \n
    If range='I', the index of the
    smallest eigenvalue to be returned. \n
    1 <= il <= iu <= N, if N > 0; il = 1 and iu = 0 if N = 0. \n
    Not referenced if range = 'A' or 'V'.

    * @param[in] iu
    iu is int* \n
    If range='I', the index of the
    largest eigenvalue to be returned. \n
    1 <= il <= iu <= N, if N > 0; il = 1 and iu = 0 if N = 0. \n
    Not referenced if range = 'A' or 'V'.

    * @param[in] abstol
    abstol is float\double \n
    The absolute error tolerance for the eigenvalues.
    An approximate eigenvalue is accepted as converged
    when it is determined to lie in an interval [a,b]
    of width less than or equal to \n \n
    abstol + EPS *   max( |a|,|b| ) , \n \n
    where EPS is the machine precision.  If abstol is less than
    or equal to zero, then  EPS*|T|  will be used in its place,
    where |T| is the 1-norm of the tridiagonal matrix obtained
    by reducing A to tridiagonal form. \n \n
    See "Computing Small Singular Values of Bidiagonal Matrices
    with Guaranteed High Relative Accuracy," by Demmel and
    Kahan, LAPACK Working Note #3. \n \n
    If high relative accuracy is important, set abstol to
    SLAMCH( 'Safe minimum' ).  Doing so will guarantee that
    eigenvalues are computed to high relative accuracy when
    possible in future releases.  The current code does not
    make any guarantees about high relative accuracy, but
    furutre releases will. See J. Barlow and J. Demmel,
    "Computing Accurate Eigensystems of Scaled Diagonally
    Dominant Matrices", LAPACK Working Note #7, for a discussion
    of which matrices define their eigenvalues to high relative
    accuracy.

    * @param[out] m
    m is int* \n
    The total number of eigenvalues found.  0 <= M <= N. \n
    If range = 'A', M = N, and if range = 'I', M = iu-il+1.

    * @param[out] w
    w is float/double array, dimension (n) \n
    The first M elements contain the selected eigenvalues in
    ascending order.

    * @param[out] z
    z is scomplex/dcomplex array, dimension (ldz, max(1,m)) \n
    If jobz = 'V', then if info = 0, the first M columns of Z
    contain the orthonormal eigenvectors of the matrix a
    corresponding to the selected eigenvalues, with the i-th
    column of Z holding the eigenvector associated with W(i). \n
    If jobz = 'N', then Z is not referenced. \n
    Note: the user must ensure that at least max(1,m) columns are
    supplied in the array Z; if range = 'V', the exact value of M
    is not known in advance and an upper bound must be used.

    * @param[in] ldz
    ldz is int* \n
    The leading dimension of the array Z.  ldz >= 1, and if
    jobz = 'V', ldz >= max(1,n).

    * @param[out] isuppz
    isuppz is int array, dimension ( 2*max(1,m) ) \n
    The support of the eigenvectors in Z, i.e., the indices
    indicating the nonzero elements in Z. The i-th eigenvector
    is nonzero only in elements isuppz( 2*i-1 ) through
    isuppz( 2*i ). This is an output of CSTEMR (tridiagonal
    matrix). The support of the eigenvectors of A is typically
    1:N because of the unitary transformations applied by CUNMTR. \n
    Implemented only for range = 'A' or 'I' and iu - il = N - 1

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The length of the array work.  lwork >= max(1,2*n). \n
    For optimal efficiency, lwork >= (NB+1)*N, \n
    where NB is the max of the blocksize for CHETRD and for
    CUNMTR as returned by ILAENV. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal sizes of the work, rwork and
    iwork arrays, returns these values as the first entries of
    the work, rwork and iwork arrays, and no error message
    related to lwork or lrwork or liwork is issued by XERBLA.

    * @param[out] rwork
    rwork is float/double array, dimension (max(1,lrwork)) \n
    On exit, if info = 0, rwork(1) returns the optimal
    (and minimal) lrwork.

    * @param[in] lrwork
    lrwork is int* \n
    The length of the array rwork.  lrwork >= max(1,24*N). \n \n
    If lrwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal sizes of the work, rwork
    and iwork arrays, returns these values as the first entries
    of the work, rwork and iwork arrays, and no error message
    related to lwork or lrwork or liwork is issued by XERBLA.

    * @param[out] iwork
    iwork is int array, dimension (max(1,liwork)) \n
    On exit, if info = 0, iwork(1) returns the optimal
    (and minimal) liwork.

    * @param[in] liwork
    liwork is int* \n
    The dimension of the array iwork.  liwork >= max(1,10*N). \n \n
    If liwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal sizes of the work, rwork
    and iwork arrays, returns these values as the first entries
    of the work, rwork and iwork arrays, and no error message
    related to lwork or lrwork or liwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  if info = -i, the i-th argument had an illegal value \n
    \> 0:  Internal error
    *  */
template< typename Ta, typename Tb >
int heevr( char* jobz, char* range, char* uplo, int* n, Ta* a, int* lda, Tb*  vl, Tb*  vu, int* il, int* iu, Tb*  abstol, int* m, Tb* w, Ta* z, int* ldz, int* isuppz, Ta* work, int* lwork, Tb*  rwork, int* lrwork, int* iwork, int* liwork, int* info )
{
  return heevr( jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info );
}

/*! @brief Bidiagonal QR algorithm
    *
    * @details
    * \b Purpose:
    * \verbatim
        Bidiagonal QR algorithm.
        Computation of singular values and, optionally, the right and/or left singular vectors
        from the singular value decomposition (SVD) of a real n-by-n (upper or lower) bidiagonal
        matrix b using the implicit zero-shift QR algorithm.  The SVD of B has the form

        B = Q * S * P**T

        where S is the diagonal matrix of singular values, Q is an orthogonal matrix of left
        singular vectors, and P is an orthogonal matrix of right singular vectors.  If left
        singular vectors are requested, this subroutine actually returns U*Q instead of Q, and,
        if right singular vectors are requested, this subroutine returns P**T*vt instead of P**T,
        for given real input matrices U and vt.  When U and vt are the orthogonal matrices that
        reduce a general matrix a to bidiagonal form:

        A = U*B*vt, as computed by SGEBRD, then

        A = (U*Q) * S * (P**T*vt)

        is the SVD of A.  Optionally, the subroutine may also compute Q**T*C for a given real input matrix c.

        See,
        "Computing Small Singular Values of Bidiagonal Matrices With Guaranteed High Relative
        Accuracy," by J. Demmel and W. Kahan, LAPACK Working Note #3 (or SIAM J. Sci. Statist.
        Comput. vol. 11, no. 5, pp. 873-912, Sept 1990) and
        "Accurate singular values and differential qd algorithms," by B. Parlett and V. Fernando,
        Technical Report CPAM-554, Mathematics Department, University of California at Berkeley,
        July 1992 for a detailed description of the algorithm.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U':  B is upper bidiagonal; \n
    = 'L':  B is lower bidiagonal.

    * @param[in] n
    n is int* \n
    The order of the matrix B.  n >= 0.

    * @param[in] ncvt
    ncvt is int* \n
    The number of columns of the matrix vt. ncvt >= 0.

    * @param[in] nru
    nru is int* \n
    The number of rows of the matrix U. nru >= 0.

    * @param[in] ncc
    ncc is int* \n
    The number of columns of the matrix c. ncc >= 0.

    * @param[in,out] d
    d is float/double array, dimension (n) \n
    On entry, the n diagonal elements of the bidiagonal matrix b. \n
    On exit, if info=0, the singular values of B in decreasing
    order.

    * @param[in,out] e
    e is float/double array, dimension (n-1) \n
    On entry, the N-1 offdiagonal elements of the bidiagonal
    matrix b. \n
    On exit, if info = 0, E is destroyed; if info > 0, D and E
    will contain the diagonal and superdiagonal elements of a
    bidiagonal matrix orthogonally equivalent to the one given
    as input.

    * @param[in,out] vt
    vt is float/double array, dimension (ldvt, ncvt) \n
    On entry, an n-by-ncvt matrix vt. \n
    On exit, vt is overwritten by P**T * vt. \n
    Not referenced if ncvt = 0.

    * @param[in] ldvt
    ldvt is int* \n
    The leading dimension of the array vt. \n
    ldvt >= max(1,n) if ncvt > 0; ldvt >= 1 if ncvt = 0.

    * @param[in,out] u
    u is float/double array, dimension (ldu, n) \n
    On entry, an nru-by-n matrix U. \n
    On exit, U is overwritten by U * Q. \n
    Not referenced if nru = 0.

    * @param[in] ldu
    ldu is int* \n
    The leading dimension of the array U.  ldu >= max(1,nru).

    * @param[in,out] c
    c is float/doublearray, dimension (ldc, ncc) \n
    On entry, an n-by-ncc matrix c. \n
    On exit, c is overwritten by Q**T * C. \n
    Not referenced if ncc = 0.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. \n
    ldc >= max(1,n) if ncc > 0; ldc >=1 if ncc = 0.

    * @param[out] rwork
    rwork is float/double array, dimension (4*N)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  If info = -i, the i-th argument had an illegal value \n
    \> 0: \n
    if ncvt = nru = ncc = 0,
    = 1, a split was marked by a positive value in E
    = 2, current block of Z not diagonalized after 30*N
    iterations (in inner while loop)
    = 3, termination criterion of outer while loop not met
    (program created more than N unreduced blocks)
    else ncvt = nru = ncc = 0,
    the algorithm did not converge; D and E contain the
    elements of a bidiagonal matrix which is orthogonally
    similar to the input matrix b;  if info = i, i
    elements of E have not converged to zero.
    * \par
    * \verbatim
        Internal Parameters:
        ====================
        TOLMUL  T*, default = max(10,min(100,EPS**(-1/8)))
        TOLMUL controls the convergence criterion of the QR loop.
        If it is positive, TOLMUL*EPS is the desired relative
        precision in the computed singular values.
        If it is negative, abs(TOLMUL*EPS*sigma_max) is the
        desired absolute accuracy in the computed singular
        values (corresponds to relative accuracy
        abs(TOLMUL*EPS) in the largest singular value.
        abs(TOLMUL) should be between 1 and 1/EPS, and preferably
        between 10 (for fast convergence) and .1/EPS
        (for there to be some accuracy in the results).
        Default is to lose at either one eighth or 2 of the
        available decimal digits in each computed singular value
        (whichever is smaller).

        MAXITR  int*, default = 6
        MAXITR controls the maximum number of passes of the
        algorithm through its inner loop. The algorithms stops
        (and so fails to converge) if the number of passes
        through the inner loop exceeds MAXITR*N**2.
    \endverbatim
    *  */
template< typename T >
int bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, T* d, T* e, T* vt, int* ldvt, T* u, int* ldu, T* c, int* ldc, T* rwork, int* info )
{
  return  bdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info );
}

/*! @brief Bidiagonal QR algorithm
    *

    * @details
    * \b Purpose:
    * \verbatim
        Bidiagonal QR algorithm.
        Computation of singular values and, optionally, the right and/or left singular vectors
        from the singular value decomposition (SVD) of a real n-by-n (upper or lower) bidiagonal
        matrix b using the implicit zero-shift QR algorithm.  The SVD of B has the form

            B = Q * S * P**H

        where S is the diagonal matrix of singular values, Q is an orthogonal matrix of left
        singular vectors,and P is an orthogonal matrix of right singular vectors.  If left
        singular vectors are requested, this subroutine actually returns U*Q instead of Q, and,
        if right singular vectors are requested, this subroutine returns P**H*vt instead of P**H,
        for given complex input matrices U and vt.  When U and vt are the unitary matrices that
        reduce a general matrix a to bidiagonal form:

            A = U*B*vt, as computed by CGEBRD, then

            A = (U*Q) * S * (P**H*vt)

        is the SVD of A. Optionally, the subroutine may also compute Q**H*C for a given complex
        input matrix c.

        See,
        "Computing Small Singular Values of Bidiagonal Matrices With Guaranteed High Relative
        Accuracy," by J. Demmel and W. Kahan, LAPACK Working Note #3 (or SIAM J. Sci. Statist.
        Comput. vol. 11, no. 5, pp. 873-912, Sept 1990) and
        "Accurate singular values and differential qd algorithms," by B. Parlett and V. Fernando,
        Technical Report CPAM-554, Mathematics Department, University of California at Berkeley,
        July 1992 for a detailed description of the algorithm.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U':  B is upper bidiagonal; \n
    = 'L':  B is lower bidiagonal.

    * @param[in] n
    n is int* \n
    The order of the matrix B.  n >= 0.

    * @param[in] ncvt
    ncvt is int* \n
    The number of columns of the matrix vt. ncvt >= 0.

    * @param[in] nru
    nru is int* \n
    The number of rows of the matrix U. nru >= 0.

    * @param[in] ncc
    ncc is int* \n
    The number of columns of the matrix c. ncc >= 0.

    * @param[in,out] d
    d is float/double array, dimension (n) \n
    On entry, the n diagonal elements of the bidiagonal matrix b. \n
    On exit, if info=0, the singular values of B in decreasing
    order.

    * @param[in,out] e
    e is float/double array, dimension (n-1) \n
    On entry, the N-1 offdiagonal elements of the bidiagonal
    matrix b. \n
    On exit, if info = 0, E is destroyed; if info > 0, D and E
    will contain the diagonal and superdiagonal elements of a
    bidiagonal matrix orthogonally equivalent to the one given
    as input.

    * @param[in,out] vt
    vt is scomplex/dcomplex array, dimension (ldvt, ncvt) \n
    On entry, an n-by-ncvt matrix vt. \n
    On exit, vt is overwritten by P**H * vt. \n
    Not referenced if ncvt = 0.

    * @param[in] ldvt
    ldvt is int* \n
    The leading dimension of the array vt. \n
    ldvt >= max(1,n) if ncvt > 0; ldvt >= 1 if ncvt = 0.

    * @param[in,out] u
    u is scomplex/dcomplex array, dimension (ldu, n) \n
    On entry, an nru-by-n matrix U. \n
    On exit, U is overwritten by U * Q. \n
    Not referenced if nru = 0.

    * @param[in] ldu
    ldu is int* \n
    The leading dimension of the array U.  ldu >= max(1,nru).

    * @param[in,out] c
    c is scomplex/dcomplex array, dimension (ldc, ncc) \n
    On entry, an n-by-ncc matrix c. \n
    On exit, c is overwritten by Q**H * C. \n
    Not referenced if ncc = 0.

    * @param[in] ldc
    ldc is int* \n
    The leading dimension of the array c. \n
    ldc >= max(1,n) if ncc > 0; ldc >=1 if ncc = 0.

    * @param[out] rwork
    rwork is float/double array, dimension (4*N)

    * @param[out] info
    info is int* \n
    = 0:  successful exit \n
    < 0:  If info = -i, the i-th argument had an illegal value \n
    \> 0:  the algorithm did not converge; D and E contain the
        elements of a bidiagonal matrix which is orthogonally
        similar to the input matrix b;  if info = i, i
        elements of E have not converged to zero.
    * \par
    * \verbatim
        Internal Parameters:
        ====================
        TOLMUL  float*, default = max(10,min(100,EPS**(-1/8)))
                TOLMUL controls the convergence criterion of the QR loop.
                If it is positive, TOLMUL*EPS is the desired relative
                precision in the computed singular values.
                If it is negative, abs(TOLMUL*EPS*sigma_max) is the
                desired absolute accuracy in the computed singular
                values (corresponds to relative accuracy
                abs(TOLMUL*EPS) in the largest singular value.
                abs(TOLMUL) should be between 1 and 1/EPS, and preferably
                between 10 (for fast convergence) and .1/EPS
                (for there to be some accuracy in the results).
                Default is to lose at either one eighth or 2 of the
                available decimal digits in each computed singular value
                (whichever is smaller).

        MAXITR  int*, default = 6
                MAXITR controls the maximum number of passes of the
                algorithm through its inner loop. The algorithms stops
                (and so fails to converge) if the number of passes
                through the inner loop exceeds MAXITR*N**2.
    \endverbatim
    *  */
template< typename Ta, typename Tb >
int bdsqr( char* uplo, int* n, int* ncvt, int* nru, int* ncc, Tb* d, Tb* e, Ta* vt, int* ldvt, Ta* u, int* ldu, Ta* c, int* ldc, Tb*  rwork, int* info )
{
  return bdsqr( uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info );
}

/*! @brief Bidiagonal divide-and-conquer algorithm
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of singular value decomposition (SVD) of a real n-by-n (upper or lower)
        bidiagonal matrix
        B:  B = U * S * vt,
        using a divide and conquer method, where S is a diagonal matrix with non-negative diagonal
        elements (the singular values of B), and U and vt are orthogonal matrices of left and
        right singular vectors, respectively. SBDSDC can be used to compute all singular values,
        and optionally, singular vectors or singular vectors in compact form.

        This code makes very mild assumptions about floating point arithmetic. It will work on
        machines with a guard digit in add/subtract, or on those binary machines without guard
        digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could
        conceivably fail on hexadecimal or decimal machines without guard digits, but we know of
        none.  See SLASD3 for details.

        The code currently calls SLASDQ if singular values only are desired. However, it can be
        slightly modified to compute singular values using the divide and conquer method.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    = 'U':  B is upper bidiagonal. \n
    = 'L':  B is lower bidiagonal.

    * @param[in] compq
    compq is char* \n
    Specifies whether singular vectors are to be computed
    as follows: \n
    = 'N':  Compute singular values only; \n
    = 'P':  Compute singular values and compute singular
    vectors in compact form; \n
    = 'I':  Compute singular values and singular vectors.

    * @param[in] n
    n is int* \n
    The order of the matrix B.  n >= 0.

    * @param[in,out] d
    d is float/double array, dimension (n) \n
    On entry, the n diagonal elements of the bidiagonal matrix b. \n
    On exit, if info=0, the singular values of B.

    * @param[in,out] e
    e is float/double array, dimension (n-1) \n
    On entry, the elements of E contain the offdiagonal
    elements of the bidiagonal matrix whose SVD is desired. \n
    On exit, E has been destroyed.

    * @param[out] u
    U is float/double array, dimension (ldu,N) \n
    If  compq = 'I', then: \n
    On exit, if info = 0, U contains the left singular vectors
    of the bidiagonal matrix. \n
    For other values of compq, U is not referenced.

    * @param[in] ldu
    ldu is int* \n
    The leading dimension of the array U.  ldu >= 1. \n
    If singular vectors are desired, then ldu >= max( 1, N ).

    * @param[out] vt
    vt is float/double array, dimension (ldvt,N) \n
    If  compq = 'I', then: \n
    On exit, if info = 0, vt**T contains the right singular
    vectors of the bidiagonal matrix. \n
    For other values of compq, vt is not referenced.

    * @param[in] ldvt
    ldvt is int* \n
    The leading dimension of the array vt.  ldvt >= 1. \n
    If singular vectors are desired, then ldvt >= max( 1, N ).

    * @param[out] q
    q is float/double array, dimension (LDQ) \n
    If  compq = 'P', then: \n
    On exit, if info = 0, Q and iq contain the left
    and right singular vectors in a compact form,
    requiring O(N log N) space instead of 2*n**2. \n
    In particular, Q contains all the float data in
    LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))
    words of memory, where SMLSIZ is returned by ILAENV and
    is equal to the maximum size of the subproblems at the
    bottom of the computation tree (usually about 25). \n
    For other values of compq, Q is not referenced.

    * @param[out] iq
    iq is float/double array, dimension (LDIQ) \n
    If  compq = 'P', then: \n
    On exit, if info = 0, Q and iq contain the left
    and right singular vectors in a compact form,
    requiring O(N log N) space instead of 2*n**2. \n
    In particular, iq contains all integer data in
    LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))
    words of memory, where SMLSIZ is returned by ILAENV and
    is equal to the maximum size of the subproblems at the
    bottom of the computation tree (usually about 25). \n
    For other values of compq, iq is not referenced. \n

    * @param[out] work
    work is float/double, dimension (max(1,lwork)) \n
    If compq = 'N' then lwork >= (4 * N). \n
    If compq = 'P' then lwork >= (6 * N). \n
    If compq = 'I' then lwork >= (3 * N**2 + 4 * N).

    * @param[out] iwork
    iwork is int array, dimension (8*N) \n

    * @param[out] info
    info is int* \n
    = 0:  successful exit. \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  The algorithm failed to compute a singular value.
    The update process of divide and conquer failed.
    *  */
template< typename T >
int bdsdc( char* uplo, char* compq, int* n, T*  d, T*  e, T*  u, int* ldu, T*  vt, int* ldvt, T*  q, T*  iq, T*  work, int* iwork, int* info )
{
  return bdsdc( uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info );
}

/*! @brief General matrix singular value decomposition (QR algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of singular value decomposition (SVD) of a real m-by-n matrix a, optionally
        computing the left and/or right singular vectors. The SVD is written

        A = U * SIGMA * transpose(V)

        where SIGMA is an m-by-n matrix which is zero except for its min(m,n) diagonal elements,
        U is an M-by-M orthogonal matrix, and V is an n-by-n orthogonal matrix.  The diagonal
        elements of SIGMA are the singular values of A; they are real and non-negative, and are
        returned in descending order. The first min(m,n) columns of U and V are the left and right
        singular vectors of A.

        Note that the routine returns V**T, not V.
    \endverbatim

    * @param[in] jobu
    jobu is char* \n
    Specifies options for computing all or part of the matrix U: \n
    = 'A':  all M columns of U are returned in array U: \n
    = 'S':  the first min(m,n) columns of U (the left singular
    vectors) are returned in the array U; \n
    = 'O':  the first min(m,n) columns of U (the left singular
    vectors) are overwritten on the array a; \n
    = 'N':  no columns of U (no left singular vectors) are
    computed.

    * @param[in] jobv
    jobv is char* \n
    Specifies options for computing all or part of the matrix
    V**T: \n
    = 'A':  all N rows of V**T are returned in the array vt; \n
    = 'S':  the first min(m,n) rows of V**T (the right singular
    vectors) are returned in the array vt; \n
    = 'O':  the first min(m,n) rows of V**T (the right singular
    vectors) are overwritten on the array a; \n
    = 'N':  no rows of V**T (no right singular vectors) are
    computed. \n \n
    jobv and jobu cannot both be 'O'.

    * @param[in] m
    m is int* \n
    The number of rows of the input matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the input matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the m-by-n matrix a. \n
    On exit, \n
    if jobu = 'O',  A is overwritten with the first min(m,n)
    columns of U (the left singular vectors,
    stored columnwise); \n
    if jobv = 'O', A is overwritten with the first min(m,n)
    rows of V**T (the right singular vectors,
    stored rowwise); \n
    if jobu .ne. 'O' and jobv .ne. 'O', the contents of A
    are destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,m).

    * @param[out] s
    s is float/double array, dimension (min(m,n)) \n
    The singular values of A, sorted so that S(i) >= S(i+1).

    * @param[out] u
    U is float/double array, dimension (ldu,UCOL) \n
    (ldu,M) if jobu = 'A' or (ldu,min(m,n)) if jobu = 'S'. \n
    If jobu = 'A', U contains the m-by-m orthogonal matrix U; \n
    if jobu = 'S', U contains the first min(m,n) columns of U
    (the left singular vectors, stored columnwise); \n
    if jobu = 'N' or 'O', U is not referenced.

    * @param[in] ldu
    ldu is int* \n
    The leading dimension of the array U.  ldu >= 1; if
    jobu = 'S' or 'A', ldu >= M.

    * @param[out] vt
    vt is float/double array, dimension (ldvt,N) \n
    If jobv = 'A', vt contains the n-by-n orthogonal matrix
    V**T; \n
    if jobv = 'S', vt contains the first min(m,n) rows of
    V**T (the right singular vectors, stored rowwise); \n
    if jobv = 'N' or 'O', vt is not referenced.

    * @param[in] ldvt
    ldvt is int* \n
    The leading dimension of the array vt.  ldvt >= 1; if
    jobv = 'A', ldvt >= N; if jobv = 'S', ldvt >= min(m,n).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork; \n
    if info > 0, work(2:MIN(m,n)) contains the unconverged
    superdiagonal elements of an upper bidiagonal matrix b
    whose diagonal is in S (not necessarily sorted). B
    satisfies A = U * B * vt, so it has the same singular values
    as A, and singular vectors related by U and vt.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    lwork >= max(1,5*MIN(m,n)) for the paths (see comments inside code): \n
    - PATH 1  (M much larger than N, jobu='N') \n
    - PATH 1t (N much larger than M, jobv='N') \n
    lwork >= max(1,3*MIN(m,n)+max(m,n),5*MIN(m,n)) for the other paths
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] info
    info is int* \n
    = 0:  successful exit. \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  if SBDSQR did not converge, info specifies how many
    superdiagonals of an intermediate bidiagonal form B
    did not converge to zero. See the description of work
    above for details.
    *  */
template< typename T >
int gesvd( char* jobu, char* jobv, int* m, int* n, T* a, int* lda, T*  s, T* u, int* ldu, T* vt, int* ldvt, T* work, int* lwork, int* info )
{
  return gesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info );
}

/*! @brief General matrix singular value decomposition (QR algorithm)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of singular value decomposition (SVD) of a complex m-by-n matrix a, optionally
        computing the left and/or right singular vectors. The SVD is written

            A = U * SIGMA * conjugate-transpose(V)

        where SIGMA is an m-by-n matrix which is zero except for its min(m,n) diagonal elements,
        U is an M-by-M unitary matrix, and V is an n-by-n unitary matrix.  The diagonal elements
        of SIGMA are the singular values of A; they are real and non-negative, and are returned in
        descending order.  The first min(m,n) columns of U and V are the left and right singular
        vectors of A.

        Note that the routine returns V**H, not V.
    \endverbatim

    * @param[in] jobu
    jobu is char* \n
    Specifies options for computing all or part of the matrix U: \n
    = 'A':  all M columns of U are returned in array U: \n
    = 'S':  the first min(m,n) columns of U (the left singular
            vectors) are returned in the array U; \n
    = 'O':  the first min(m,n) columns of U (the left singular
            vectors) are overwritten on the array a; \n
    = 'N':  no columns of U (no left singular vectors) are
            computed.

    * @param[in] jobv
    jobv is char* \n
    Specifies options for computing all or part of the matrix
    V**H: \n
    = 'A':  all N rows of V**H are returned in the array vt; \n
    = 'S':  the first min(m,n) rows of V**H (the right singular
            vectors) are returned in the array vt; \n
    = 'O':  the first min(m,n) rows of V**H (the right singular
            vectors) are overwritten on the array a; \n
    = 'N':  no rows of V**H (no right singular vectors) are
            computed. \n \n
    jobv and jobu cannot both be 'O'.

    * @param[in] m
    m is int* \n
    The number of rows of the input matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the input matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix a. \n
    On exit, \n
    if jobu = 'O',  A is overwritten with the first min(m,n)
                    columns of U (the left singular vectors,
                    stored columnwise); \n
    if jobv = 'O', A is overwritten with the first min(m,n)
                    rows of V**H (the right singular vectors,
                    stored rowwise); \n
    if jobu .ne. 'O' and jobv .ne. 'O', the contents of A
                    are destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,m).

    * @param[out] s
    s is float/double array, dimension (min(m,n)) \n
    The singular values of A, sorted so that S(i) >= S(i+1).

    * @param[out] u
    U is scomplex/dcomplex array, dimension (ldu,UCOL) \n
    (ldu,M) if jobu = 'A' or (ldu,min(m,n)) if jobu = 'S'. \n
    If jobu = 'A', U contains the m-by-m unitary matrix U; \n
    if jobu = 'S', U contains the first min(m,n) columns of U
    (the left singular vectors, stored columnwise); \n
    if jobu = 'N' or 'O', U is not referenced.

    * @param[in] ldu
    ldu is int* \n
    The leading dimension of the array U.  ldu >= 1; if
    jobu = 'S' or 'A', ldu >= M.

    * @param[out] vt
    vt is scomplex/dcomplex array, dimension (ldvt,N) \n
    If jobv = 'A', vt contains the n-by-n unitary matrix
    V**H; \n
    if jobv = 'S', vt contains the first min(m,n) rows of
    V**H (the right singular vectors, stored rowwise); \n
    if jobv = 'N' or 'O', vt is not referenced.

    * @param[in] ldvt
    ldvt is int* \n
    The leading dimension of the array vt.  ldvt >= 1; if
    jobv = 'A', ldvt >= N; if jobv = 'S', ldvt >= min(m,n).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. \n
    lwork >=  max(1,2*MIN(m,n)+max(m,n)). \n
    For good performance, lwork should generally be larger. \n \n
    If lwork = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of the work array, returns
    this value as the first entry of the work array, and no error
    message related to lwork is issued by XERBLA.

    * @param[out] rwork
    rwork is float/double array, dimension (5*min(m,n)) \n
    On exit, if info > 0, rwork(1:MIN(m,n)-1) contains the
    unconverged superdiagonal elements of an upper bidiagonal
    matrix b whose diagonal is in S (not necessarily sorted).
    B satisfies A = U * B * vt, so it has the same singular
    values as A, and singular vectors related by U and vt.

    * @param[out] info
    info is int* \n
    = 0:  successful exit. \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  if CBDSQR did not converge, info specifies how many
        superdiagonals of an intermediate bidiagonal form B
        did not converge to zero. See the description of rwork
        above for details.
    *  */
template< typename Ta, typename Tb >
int gesvd( char* jobu, char* jobv, int* m, int* n, Ta* a, int* lda, Tb*  s, Ta* u, int* ldu, Ta* vt, int* ldvt, Ta* work, int* lwork, Tb* rwork, int* info )
{
  return gesvd( jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info );
}

/*! @brief General matrix singular value decomposition (divide-and-conquer)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of singular value decomposition (SVD) of a real m-by-n matrix a, optionally
        computing the left and right singular vectors.  If singular vectors are desired, it uses a
        divide-and-conquer algorithm.

        The SVD is written

        A = U * SIGMA * transpose(V)

        where SIGMA is an m-by-n matrix which is zero except for its min(m,n) diagonal elements,
        U is an M-by-M orthogonal matrix, and V is an n-by-n orthogonal matrix.  The diagonal
        elements of SIGMA are the singular values of A; they are real and non-negative, and are
        returned in descending order.  The first min(m,n) columns of U and V are the left and right
        singular vectors of A.

        Note that the routine returns vt = V**T, not V.

        The divide and conquer algorithm makes very mild assumptions about floating point
        arithmetic. It will work on machines with a guard digit in add/subtract, or on those
        binary machines without guard digits which subtract like the Cray X-MP, Cray Y-MP, Cray
        C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without
        guard digits, but we know of none.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    Specifies options for computing all or part of the matrix U: \n
    = 'A':  all M columns of U and all N rows of V**T are
    returned in the arrays U and vt; \n
    = 'S':  the first min(m,n) columns of U and the first
    min(m,n) rows of V**T are returned in the arrays U
    and vt; \n
    = 'O':  If m >= n, the first N columns of U are overwritten
    on the array a and all rows of V**T are returned in
    the array vt; \n
    otherwise, all columns of U are returned in the
    array U and the first M rows of V**T are overwritten
    in the array a; \n
    = 'N':  no columns of U or rows of V**T are computed.

    * @param[in] m
    m is int* \n
    The number of rows of the input matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the input matrix a.  n >= 0.

    * @param[in,out] a
    a is float/double array, dimension (lda,n) \n
    On entry, the m-by-n matrix a. \n
    On exit, \n
    if jobz = 'O',  A is overwritten with the first N columns
    of U (the left singular vectors, stored
    columnwise) if m >= n;
    A is overwritten with the first M rows
    of V**T (the right singular vectors, stored
    rowwise) otherwise. \n
    if jobz .ne. 'O', the contents of A are destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,m).

    * @param[out] s
    s is float/double array, dimension (min(m,n)) \n
    The singular values of A, sorted so that S(i) >= S(i+1).

    * @param[out] u
    U is float/double array, dimension (ldu,UCOL) \n
    UCOL = M if jobz = 'A' or jobz = 'O' and M < N; \n
    UCOL = min(m,n) if jobz = 'S'. \n
    If jobz = 'A' or jobz = 'O' and M < N, U contains the M-by-M
    orthogonal matrix U; \n
    if jobz = 'S', U contains the first min(m,n) columns of U
    (the left singular vectors, stored columnwise); \n
    if jobz = 'O' and m >= n, or jobz = 'N', U is not referenced.

    * @param[in] ldu
    ldu is int* \n
    The leading dimension of the array U.  ldu >= 1; if
    jobz = 'S' or 'A' or jobz = 'O' and M < N, ldu >= M.

    * @param[out] vt
    vt is float/double array, dimension (ldvt,N) \n
    If jobz = 'A' or jobz = 'O' and m >= n, vt contains the
    n-by-n orthogonal matrix V**T; \n
    if jobz = 'S', vt contains the first min(m,n) rows of
    V**T (the right singular vectors, stored rowwise); \n
    if jobz = 'O' and M < N, or jobz = 'N', vt is not referenced.

    * @param[in] ldvt
    ldvt is int* \n
    The leading dimension of the array vt.  ldvt >= 1; \n
    if jobz = 'A' or jobz = 'O' and m >= n, ldvt >= N; \n
    if jobz = 'S', ldvt >= min(m,n).

    * @param[out] work
    work is float/double array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork;

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= 1. \n
    If lwork = -1, a workspace query is assumed.  The optimal
    size for the work array is calculated and stored in work(1),
    and no other work except argument checking is performed. \n \n
    Let mx = max(m,n) and mn = min(m,n). \n
    If jobz = 'N', lwork >= 3*mn + max( mx, 7*mn ). \n
    If jobz = 'O', lwork >= 3*mn + max( mx, 5*mn*mn + 4*mn ). \n
    If jobz = 'S', lwork >= 4*mn*mn + 7*mn. \n
    If jobz = 'A', lwork >= 4*mn*mn + 6*mn + mx. \n
    These are not tight minimums in all cases; see comments inside code.
    For good performance, lwork should generally be larger;
    a query is recommended.

    * @param[out] iwork
    iwork is int array, dimension (8*min(m,n))

    * @param[out] info
    info is int* \n
    = 0:  successful exit. \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  SBDSDC did not converge, updating process failed.
    *  */
template< typename T >
int gesdd( char* jobz, int* m, int* n, T* a, int* lda, T*  s, T* u, int* ldu, T* vt, int* ldvt, T* work, int* lwork, int* iwork, int* info )
{
  return gesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info );
}

/*! @brief General matrix singular value decomposition (divide-and-conquer)
    *

    * @details
    * \b Purpose:
    * \verbatim
        Computation of singular value decomposition (SVD) of a complex m-by-n matrix a, optionally
        computing the left and/or right singular vectors, by using divide-and-conquer method.
        The SVD is written

            A = U * SIGMA * conjugate-transpose(V)

        where SIGMA is an m-by-n matrix which is zero except for its min(m,n) diagonal elements,
        U is an M-by-M unitary matrix, and V is an n-by-n unitary matrix.  The diagonal elements
        of SIGMA are the singular values of A; they are real and non-negative, and are returned in
        descending order.  The first min(m,n) columns of U and V are the left and right singular
        vectors of A.

        Note that the routine returns vt = V**H, not V.

        The divide and conquer algorithm makes very mild assumptions about floating point
        arithmetic. It will work on machines with a guard digit in add/subtract, or on those
        binary machines without guard digits which subtract like the Cray X-MP, Cray Y-MP, Cray
        C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without
        guard digits, but we know of none.
    \endverbatim

    * @param[in] jobz
    jobz is char* \n
    Specifies options for computing all or part of the matrix U: \n
    = 'A':  all M columns of U and all N rows of V**H are
            returned in the arrays U and vt; \n
    = 'S':  the first min(m,n) columns of U and the first
            min(m,n) rows of V**H are returned in the arrays U
            and vt; \n
    = 'O':  If m >= n, the first N columns of U are overwritten
            in the array a and all rows of V**H are returned in
            the array vt; \n
            otherwise, all columns of U are returned in the
            array U and the first M rows of V**H are overwritten
            in the array a; \n
    = 'N':  no columns of U or rows of V**H are computed.

    * @param[in] m
    m is int* \n
    The number of rows of the input matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the input matrix a.  n >= 0.

    * @param[in,out] a
    a is scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the m-by-n matrix a. \n
    On exit, \n
    if jobz = 'O',  A is overwritten with the first N columns
                    of U (the left singular vectors, stored
                    columnwise) if m >= n;
                    A is overwritten with the first M rows
                    of V**H (the right singular vectors, stored
                    rowwise) otherwise. \n
    if jobz .ne. 'O', the contents of A are destroyed.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,m).

    * @param[out] s
    s is float/double array, dimension (min(m,n)) \n
    The singular values of A, sorted so that S(i) >= S(i+1).

    * @param[out] u
    U is scomplex/dcomplex array, dimension (ldu,UCOL) \n
    UCOL = M if jobz = 'A' or jobz = 'O' and M < N; \n
    UCOL = min(m,n) if jobz = 'S'. \n
    If jobz = 'A' or jobz = 'O' and M < N, U contains the M-by-M
    unitary matrix U; \n
    if jobz = 'S', U contains the first min(m,n) columns of U
    (the left singular vectors, stored columnwise); \n
    if jobz = 'O' and m >= n, or jobz = 'N', U is not referenced.

    * @param[in] ldu
    ldu is int* \n
    The leading dimension of the array U.  ldu >= 1; \n
    if jobz = 'S' or 'A' or jobz = 'O' and M < N, ldu >= M.

    * @param[out] vt
    vt is scomplex/dcomplex array, dimension (ldvt,N) \n
    If jobz = 'A' or jobz = 'O' and m >= n, vt contains the
    n-by-n unitary matrix V**H; \n
    if jobz = 'S', vt contains the first min(m,n) rows of
    V**H (the right singular vectors, stored rowwise); \n
    if jobz = 'O' and M < N, or jobz = 'N', vt is not referenced.

    * @param[in] ldvt
    ldvt is int* \n
    The leading dimension of the array vt.  ldvt >= 1; \n
    if jobz = 'A' or jobz = 'O' and m >= n, ldvt >= N; \n
    if jobz = 'S', ldvt >= min(m,n).

    * @param[out] work
    work is scomplex/dcomplex array, dimension (max(1,lwork)) \n
    On exit, if info = 0, work(1) returns the optimal lwork.

    * @param[in] lwork
    lwork is int* \n
    The dimension of the array work. lwork >= 1. \n
    If lwork = -1, a workspace query is assumed.  The optimal
    size for the work array is calculated and stored in work(1),
    and no other work except argument checking is performed. \n \n
    Let mx = max(m,n) and mn = min(m,n). \n
    If jobz = 'N', lwork >= 2*mn + mx. \n
    If jobz = 'O', lwork >= 2*mn*mn + 2*mn + mx. \n
    If jobz = 'S', lwork >=   mn*mn + 3*mn. \n
    If jobz = 'A', lwork >=   mn*mn + 2*mn + mx. \n
    These are not tight minimums in all cases; see comments inside code. \n
    For good performance, lwork should generally be larger;
    a query is recommended.

    * @param[out] rwork
    rwork is float/double array, dimension (max(1,lrwork)) \n
    Let mx = max(m,n) and mn = min(m,n). \n
    If jobz = 'N',    lrwork >= 5*mn (LAPACK <= 3.6 needs 7*mn); \n
    else if mx >> mn, lrwork >= 5*mn*mn + 5*mn; \n
    else              lrwork >= max( 5*mn*mn + 5*mn,
                                    2*mx*mn + 2*mn*mn + mn ).

    * @param[out] iwork
    iwork is int array, dimension (8*min(m,n))

    * @param[out] info
    info is int* \n
    = 0:  successful exit. \n
    < 0:  if info = -i, the i-th argument had an illegal value. \n
    \> 0:  The updating process of SBDSDC did not converge.
    *  */
template< typename Ta, typename Tb >
int gesdd( char* jobz, int* m, int* n, Ta* a, int* lda, Tb*  s, Ta* u, int* ldu, Ta* vt, int* ldvt, Ta* work, int* lwork, Tb*  rwork, int* iwork, int* info )
{
  return gesdd( jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info );
}

/*! @brief Swap rows
    *

    * @details
    * \b Purpose:
    * \verbatim
        Perform a series of row interchanges on the matrix a.
        One row interchange is initiated for each of rows k1 through k2 of A.
    \endverbatim

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.

    * @param[in,out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On entry, the matrix of column dimension N to which the row
    interchanges will be applied. \n
    On exit, the permuted matrix.

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.

    * @param[in] k1
    k1 is int* \n
    The first element of ipiv for which a row interchange will
    be done.

    * @param[in] k2
    k2 is int* \n
    (k2-k1+1) is the number of elements of ipiv for which a row
    interchange will be done.

    * @param[in] ipiv
    ipiv is int array, dimension (k1+(k2-k1)*abs(incx)) \n
    The vector of pivot indices. Only the elements in positions
    k1 through k1+(k2-k1)*abs(incx) of ipiv are accessed. \n
    ipiv(k1+(K-k1)*abs(incx)) = L implies rows K and L are to be
    interchanged.

    * @param[in] incx
    incx is int* \n
    The increment between successive values of ipiv. If incx
    is negative, the pivots are applied in reverse order.
    *  */
template< typename T >
int laswp( int* n, T* a, int* lda, int* k1, int* k2, int* ipiv, int* incx )
{
  return laswp( n, a, lda, k1, k2, ipiv, incx );
}

/*! @brief Initialize the off-diagonal elements and the diagonal elements of a matrix to given values
    *

    * @details
    * \b Purpose:
    * \verbatim
        Initialize an m-by-n matrix a to beta on the diagonal and alpha on the offdiagonals.
    \endverbatim

    * @param[in] uplo
    uplo is char* \n
    Specifies the part of the matrix a to be set. \n
    = 'U':   Upper triangular part is set; the strictly lower
    triangular part of a is not changed. \n
    = 'L':   Lower triangular part is set; the strictly upper
    triangular part of a is not changed. \n
    Otherwise:  All of the matrix a is set.

    * @param[in] m
    m is int* \n
    The number of rows of the matrix a.  m >= 0.

    * @param[in] n
    n is int* \n
    The number of columns of the matrix a.  n >= 0.

    * @param[in] alpha
    alpha is float/double/scomplex/dcomplex* \n
    The constant to which the offdiagonal elements are to be set.

    * @param[in] beta
    beta is float/double/scomplex/dcomplex* \n
    The constant to which the diagonal elements are to be set.

    * @param[out] a
    a is float/double/scomplex/dcomplex array, dimension (lda,n) \n
    On exit, the leading m-by-n submatrix of A is set as follows: \n \n
    if uplo = 'U', A(i,j) = alpha, 1<=i<=j-1, 1<=j<=n, \n
    if uplo = 'L', A(i,j) = alpha, j+1<=i<=m, 1<=j<=n, \n
    otherwise,  A(i,j) = alpha, 1<=i<=m, 1<=j<=n, i.ne.j, \n \n
    and, for all uplo, A(i,i) = beta, 1<=i<=min(m,n).

    * @param[in] lda
    lda is int* \n
    The leading dimension of the array a.  lda >= max(1,m).
    *  */
template< typename T >
int laset( char* uplo, int* m, int* n, T* alpha, T* beta, T* a, int* lda )
{
  return laset( uplo, m, n, alpha, beta, a, lda );
}

}  // namespace libflame
#endif  //  #ifndef INTERFACE_HH


