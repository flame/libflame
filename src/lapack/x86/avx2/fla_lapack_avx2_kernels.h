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
int fla_zscal_ix1_avx2(integer *n, doublecomplex *alpha,
                       doublecomplex *x);
int fla_dgeqrf_small_avx2(integer *m, integer *n,
                          doublereal *a, integer *lda,
                          doublereal *tau, doublereal *work);
#endif

