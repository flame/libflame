/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_lapack_avx2_kernels.h
 *  @brief AVX2 Kernel Declarations.
 *  */

#include "immintrin.h"

#ifdef FLA_ENABLE_AMD_OPT
int fla_dhrot3_avx2(integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *tau);
int fla_drot_avx2(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy, doublereal *c__, doublereal *s);
int fla_zscal_avx2(integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx);
int FLA_LU_piv_small_z_avx2( integer *m, integer *n, doublecomplex *a, integer *lda, integer *ipiv, integer *info);
#endif

