/******************************************************************************
 * * Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_lapack_x86_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#ifdef FLA_ENABLE_AMD_OPT
int fla_dhrot3(integer *n,
               doublereal *a, integer *lda,
               doublereal *v, doublereal *tau);
int fla_drot(integer *n,
             doublereal *dx, integer *incx,
             doublereal *dy, integer *incy,
             doublereal *c__, doublereal *s);
int fla_zscal(integer *n, doublecomplex *alpha,
              doublecomplex *x, integer *incx);
int fla_dgeqrf_small(integer *m, integer *n,
                     doublereal *a, integer *lda,
                     doublereal *tau, doublereal *work);
int fla_sscal(integer *n, real *alpha, 
              real *x, integer *incx);
int fla_sger(integer *m, integer *n, real *alpha,
             real *x, integer *incx, real *y, 
             integer *incy, real *a, integer *lda);
int fla_dgetrf_small_simd(integer *m, integer *n,
                     doublereal *a, integer *lda,
                     integer *ipiv, integer *info);
#endif
