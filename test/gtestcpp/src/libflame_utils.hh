/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_utils.hh
 *  libflame_utils.hh defines all the overloaded CPP functions to be invoked
 *  from template interfaces other than Libflame APIs.
 *  */
#ifndef LIBFLAME_UTILS_HH
#define LIBFLAME_UTILS_HH

#include  <FLAME.h>

namespace libflame_utils {
  
// Level 3 BLAS gemm()
inline void gemm(char* transa, char* transb, integer* m, integer* n, integer* k, float* alpha, float* a, integer* lda, float* b, integer* ldb, float* beta, float* c, integer* ldc)
{
  sgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
inline void gemm( char* transa, char* transb, integer* m, integer* n, integer* k, double* alpha, double* a, integer* lda, double* b, integer* ldb, double* beta, double* c, integer* ldc)
{
  dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
inline void gemm(char* transa, char* transb, integer* m, integer* n, integer* k, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* beta, scomplex* c, integer* ldc)
{
  cgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
inline void gemm(char* transa, char* transb, integer* m, integer* n, integer* k, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* beta, dcomplex* c, integer* ldc)
{
  zgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

// Level 1 BLAS copy()
inline void copy(integer* n, float* x, integer* incx, float* y, integer* incy)
{
  scopy_(n, x, incx, y, incy);
}
inline void copy(integer* n, double* x, integer* incx, double* y, integer* incy)
{
  dcopy_(n, x, incx, y, incy);
}
inline void copy(integer* n, scomplex* x, integer* incx, scomplex* y, integer* incy)
{
  ccopy_(n, x, incx, y, incy);
}
inline void copy(integer* n, dcomplex* x, integer* incx, dcomplex* y, integer* incy)
{
  zcopy_(n, x, incx, y, incy);
}
} //namespace libflame_utils
#endif  //  #ifndef LIBFLAME_UTILS_HH