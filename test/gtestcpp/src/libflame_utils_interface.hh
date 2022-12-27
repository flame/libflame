/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_utils_interface.hh
 *  libflame_utils_interface .hh defines all CPP templated public interfaces
 *  other than Libflame APIs.
 *  */
#ifndef LIBFLAME_UTILS_INTERFACE_HH
#define LIBFLAME_UTILS_INTERFACE_HH

#include "libflame_utils.hh"

namespace libflame_utils {

/*! @brief GEMM performs one of the matrix-matrix operations

 * @details
 * \b Purpose:
    \verbatim
    GEMM performs one of the matrix-matrix operations

    C := alpha*op( A )*op( B ) + beta*C,

    where  op( X ) is one of

      op( X ) = X   or   op( X ) = X**T,

    alpha and beta are scalars, and A, B and C are matrices, with op( A )
    an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
    \endverbatim

 * @param[in] transa
          TRANSA is CHARACTER*1
          On entry, TRANSA specifies the form of op( A ) to be used in
          the matrix multiplication as follows:

              TRANSA = 'N' or 'n',  op( A ) = A.

              TRANSA = 'T' or 't',  op( A ) = A**T.

              TRANSA = 'C' or 'c',  op( A ) = A**T.
 * @param[in] transb
          TRANSB is CHARACTER*1
          On entry, TRANSB specifies the form of op( B ) to be used in
          the matrix multiplication as follows:

              TRANSB = 'N' or 'n',  op( B ) = B.

              TRANSB = 'T' or 't',  op( B ) = B**T.

              TRANSB = 'C' or 'c',  op( B ) = B**T.
 * @param[in] m
          M is INTEGER
          On entry,  M  specifies  the number  of rows  of the  matrix
          op( A )  and of the  matrix  C.  M  must  be at least  zero. \n
 * @param[in] n
          N is INTEGER
           On entry,  N  specifies the number  of columns of the matrix
           op( B ) and the number of columns of the matrix C. N must be
           at least zero. \n
 * @param[in] k
          K is INTEGER
           On entry,  K  specifies  the number of columns of the matrix
           op( A ) and the number of rows of the matrix op( B ). K must
           be at least  zero.
 * @param[in] alpha
           ALPHA is REAL
           On entry, ALPHA specifies the scalar alpha.
 * @param[in] a
           A is REAL array, dimension ( LDA, ka ), where ka is
           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
           part of the array  A  must contain the matrix  A,  otherwise
           the leading  k by m  part of the array  A  must contain  the
           matrix A.
 * @param[in] lda
           LDA is INTEGER
           On entry, LDA specifies the first dimension of A as declared
           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
           LDA must be at least  fla_max( 1, m ), otherwise  LDA must be at
           least  fla_max( 1, k ).
 * @param[in] b
           B is REAL array, dimension ( LDB, kb ), where kb is
           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
           part of the array  B  must contain the matrix  B,  otherwise
           the leading  n by k  part of the array  B  must contain  the
           matrix B.
 * @param[in] ldb
           LDB is INTEGER
           On entry, LDB specifies the first dimension of B as declared
           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
           LDB must be at least  fla_max( 1, k ), otherwise  LDB must be at
           least  fla_max( 1, n ).
 * @param[in] beta
           BETA is REAL
           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
           supplied as zero then C need not be set on input.
 * @param[in, out] c
           C is REAL array, dimension ( LDC, N )
           Before entry, the leading  m by n  part of the array  C must
           contain the matrix  C,  except when  beta  is zero, in which
           case C need not be set on entry.
           On exit, the array  C  is overwritten by the  m by n  matrix
           ( alpha*op( A )*op( B ) + beta*C ).  
 * @param[in] ldc
           LDC is INTEGER
           On entry, LDC specifies the first dimension of C as declared
           in  the  calling  (sub)  program.   LDC  must  be  at  least
           fla_max( 1, m ).
 
 * @return void Nothing.
 * */
template < typename T >
void gemm(char* transa, char* transb, integer* m, integer* n, integer* k, T* alpha, T* a, integer* lda, T* b, integer* ldb, T* beta, T* c, integer* ldc)
{
  gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

/*! @brief COPY copies a vector, x, to a vector, y.

 * @details
 * \b Purpose:
    \verbatim
    COPY copies a vector, x, to a vector, y.
    uses unrolled loops for increments equal to 1.
    \endverbatim

 * @param[in] n
            N is INTEGER
            number of elements in input vector(s)
 * @param[in] x
            X is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
 * @param[in] incx
            INCX is INTEGER
            storage spacing between elements of X
 * @param[out] Y
            Y is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
 * @param[in] incy
            INCY is INTEGER
            storage spacing between elements of SY
 
 * @return void Nothing.
 * */
template < typename T >
void copy(integer* n, T* x, integer* incx, T* y, integer* incy)
{
  copy(n, x, incx, y, incy);
}

} //namespace libflame_utils
#endif  //  #ifndef LIBFLAME_UTILS_INTERFACE_HH