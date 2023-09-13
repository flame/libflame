/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_lapack_avx2_kernels.h
 *  @brief AVX2 Kernel Declarations.
 *  */

#include "immintrin.h"
#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#include "fla_dgeqrf_small_avx2.h"

#if FLA_ENABLE_AMD_OPT
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
int fla_dhrot3_avx2(integer *n,
                    doublereal *a, integer *lda,
                    doublereal *v, doublereal *tau);
int fla_dscal_ix1_avx2(integer *n, doublereal *da, doublereal *dx, integer *incx);
int fla_sscal_ix1_avx2(integer *n, real *alpha, real *x);
int fla_sger_avx2(integer *m, integer *n, real *alpha,
                  real *x, integer *incx,
                  real *y, integer *incy,
                  real *a, integer *lda);
int fla_zgetrf_small_avx2(integer *m, integer *n,
                          dcomplex *a, integer *lda,
                          integer *ipiv, integer *info);
int fla_zrot_avx2(integer *n, 
                  doublecomplex *cx, integer *incx, 
                  doublecomplex *cy, integer *incy, 
                  doublereal *c__, doublecomplex *s);
int fla_zscal_avx2(integer *n, doublecomplex *alpha,
                   doublecomplex *x, integer *incx);
int fla_zscal_ix1_avx2(integer *n, doublecomplex *alpha,
                       doublecomplex *x);
integer fla_dgetrf_small_avx2(integer *m, integer *n,
                              doublereal *a, integer *lda,
                              integer *ipiv,
                              integer *info);
void fla_lu_piv_small_d_update_tr_matrix_avx2(integer i_1, integer mi, integer ni,
                                              doublereal* acur, integer lda_t);
void fla_dgesvd_small6T_avx2(integer *m, integer *n,
                             doublereal *a, integer *lda,
                             doublereal *ql, integer *ldql,
                             doublereal *s,
                             doublereal *u, integer *ldu,
                             doublereal *vt, integer *ldvt,
                             doublereal *work,
                             integer *info);
void fla_dgesvd_nn_small10_avx2(integer *m, integer *n,
                                doublereal *a, integer *lda,
                                doublereal *s,
                                doublereal *work,
                                integer *info);
void fla_dgesvd_nn_small1T_avx2(integer *m, integer *n,
                                doublereal *a, integer *lda,
                                doublereal *s,
                                doublereal *work,
                                integer *info);
integer fla_dgetrf_small_avx2( integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info);
void fla_lu_piv_small_d_update_tr_matrix_avx2(integer i_1, integer mi, integer ni, doublereal* acur, integer lda_t);
int fla_dgetrs_small_trsm_ll_avx2(char *trans, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
#endif

