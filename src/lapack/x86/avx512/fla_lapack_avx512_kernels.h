/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_lapack_avx512_kernels.h
 *  @brief AVX512 Kernel Declarations.
 *  */

#include "immintrin.h"

#ifdef FLA_ENABLE_AMD_OPT
int fla_zgetrf_small_avx512(integer *m, integer *n, doublecomplex *a, integer *lda, integer *ipiv, integer *info);
integer fla_dgetrf_small_avx512( integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info);
#endif

