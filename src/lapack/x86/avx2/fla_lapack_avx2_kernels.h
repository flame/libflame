/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_lapack_avx2_kernels.h
 *  @brief AVX2 Kernel Declarations.
 *  */

#include "immintrin.h"

#ifdef FLA_ENABLE_AMD_OPT
int fla_dhrot3_avx2(integer *n,
                    doublereal *a, integer *lda,
                    doublereal *v, doublereal *tau);
int fla_drot_avx2(integer *n,
                  doublereal *dx, integer *incx,
                  doublereal *dy, integer *incy,
                  doublereal *c__, doublereal *s);
int fla_dscal_ix1_avx2(integer *n, doublereal *da,
                       doublereal *dx, integer *incx);
int fla_dgeqrf_small_avx2(integer *m, integer *n,
                          doublereal *a, integer *lda,
                          doublereal *tau, doublereal *work);
int fla_dhrot3_avx2(integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *tau);
int fla_dscal_ix1_avx2(integer *n, doublereal *da, doublereal *dx, integer *incx);
integer fla_lu_piv_small_d_avx2( integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info);
int fla_sscal_ix1_avx2(integer *n, real *alpha, real *x);
int fla_sger_avx2(integer *m, integer *n, real *alpha, real *x, integer *incx, real *y, integer *incy, real *a, integer *lda);
int fla_zgetrf_small_avx2(integer *m, integer *n,
                          doublecomplex *a, integer *lda,
                          integer *ipiv, integer *info);
int fla_zrot_avx2(integer *n, 
		doublecomplex *cx, integer *incx, 
		doublecomplex *cy, integer *incy, 
		doublereal *c__, doublecomplex *s);
int fla_zscal_avx2(integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx);
int fla_zscal_ix1_avx2(integer *n, doublecomplex *alpha,
                       doublecomplex *x);
#endif

