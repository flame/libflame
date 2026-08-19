/******************************************************************************
 * Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_LAPACK_AVX2_KERNELS_DEFS_H
#define FLA_LAPACK_AVX2_KERNELS_DEFS_H

/*! @file fla_lapack_avx2_kernels.h
 *  @brief AVX2 Kernel Declarations.
 *  */

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#include "fla_dgeqrf_small_avx2.h"
#include "fla_gesvd_small_avx2_kernel.h"
#include "fla_lapack_x86_common.h"
#include "immintrin.h"

#if FLA_ENABLE_AMD_OPT
int fla_dhrot3_avx2(aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *v,
                    doublereal *tau);
int fla_drot_avx2(aocl_int64_t *n, doublereal *dx, aocl_int64_t *incx, doublereal *dy,
                  aocl_int64_t *incy, doublereal *c__, doublereal *s);
int fla_srot_avx2(aocl_int64_t *n, real *dx, aocl_int64_t *incx, real *dy, aocl_int64_t *incy,
                  real *c__, real *s);
int fla_zrot_avx2(aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx, dcomplex *cy,
                  aocl_int64_t *incy, doublereal *c__, dcomplex *s);
int fla_sscal_ix1_avx2(aocl_int64_t *n, real *alpha, real *x);
int fla_dscal_ix1_avx2(aocl_int64_t *n, doublereal *da, doublereal *dx, aocl_int64_t *incx);
int fla_zscal_avx2(aocl_int64_t *n, dcomplex *alpha, dcomplex *x, aocl_int64_t *incx);
int fla_zscal_ix1_avx2(aocl_int64_t *n, dcomplex *alpha, dcomplex *x);
int fla_sger_avx2(aocl_int64_t *m, aocl_int64_t *n, real *alpha, real *x, aocl_int64_t *incx,
                  real *y, aocl_int64_t *incy, real *a, aocl_int64_t *lda);
aocl_int64_t fla_dgetrf_small_avx2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a,
                                   aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
aocl_int64_t fla_sgetrf_small_avx2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda,
                                   aocl_int_t *ipiv, aocl_int64_t *info);
aocl_int64_t fla_sgetrf_4x4_avx2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda,
                                 aocl_int_t *ipiv, aocl_int64_t *info);
void fla_lu_piv_small_d_update_tr_matrix_avx2(aocl_int64_t i_1, aocl_int64_t mi, aocl_int64_t ni,
                                              doublereal *acur, aocl_int64_t lda_t,
                                              doublereal safemin);
void fla_lu_piv_small_s_update_tr_matrix_avx2(aocl_int64_t i_1, aocl_int64_t mi, aocl_int64_t ni,
                                              real *acur, aocl_int64_t lda_t, real safemin);
int fla_dgetrs_small_trsm_ll_avx2(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *a,
                                  aocl_int64_t *lda, aocl_int_t *ipiv, doublereal *b,
                                  aocl_int64_t *ldb, aocl_int64_t *info);
int fla_zgetrf_small_avx2(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                          aocl_int_t *ipiv, aocl_int64_t *info);
int fla_cgetrf_small_avx2(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                          aocl_int_t *ipiv, aocl_int64_t *info);
int fla_dgeqrf_small_avx2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                          doublereal *tau, doublereal *work);
int fla_sgeqrf_small_avx2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau,
                          real *work);
int fla_cgeqrf_small_ker(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                         scomplex *tau, scomplex *work);
int fla_zgeqrf_small_ker(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                         dcomplex *tau, dcomplex *work);
void fla_dgesvd_nn_small1T_avx2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                                doublereal *s, doublereal *work, aocl_int64_t *info);
void fla_dgesvd_small6_avx2(aocl_int64_t wntus, aocl_int64_t wntvs, aocl_int64_t *m,
                            aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *qr,
                            aocl_int64_t *ldqr, doublereal *s, doublereal *u, aocl_int64_t *ldu,
                            doublereal *vt, aocl_int64_t *ldvt, doublereal *work,
                            aocl_int64_t *info);
void fla_dgesvd_small6T_avx2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                             doublereal *ql, aocl_int64_t *ldql, doublereal *s, doublereal *u,
                             aocl_int64_t *ldu, doublereal *vt, aocl_int64_t *ldvt,
                             doublereal *work, aocl_int64_t *info);
void fla_dgesvd_xx_small10_avx2(aocl_int64_t wntu, aocl_int64_t wntv, aocl_int64_t *m,
                                aocl_int64_t *n, aocl_int64_t *ncu, doublereal *a,
                                aocl_int64_t *lda, doublereal *s, doublereal *u, aocl_int64_t *ldu,
                                doublereal *vt, aocl_int64_t *ldvt, doublereal *work,
                                aocl_int64_t *info);
void fla_dgesvd_xs_small10T_avx2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                                 doublereal *s, doublereal *u, aocl_int64_t *ldu, doublereal *vt,
                                 aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *info);
void fla_sgesvd_nn_small1T_avx2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda,
                                real *s, real *work, aocl_int64_t *info);
void fla_sgesvd_small6_avx2(aocl_int64_t wntus, aocl_int64_t wntvs, aocl_int64_t *m,
                            aocl_int64_t *n, real *a, aocl_int64_t *lda, real *qr,
                            aocl_int64_t *ldqr, real *s, real *u, aocl_int64_t *ldu, real *vt,
                            aocl_int64_t *ldvt, real *work, aocl_int64_t *info);
void fla_sgesvd_small6T_avx2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *ql,
                             aocl_int64_t *ldql, real *s, real *u, aocl_int64_t *ldu, real *vt,
                             aocl_int64_t *ldvt, real *work, aocl_int64_t *info);
void fla_sgesvd_xx_small10_avx2(aocl_int64_t wntu, aocl_int64_t wntv, aocl_int64_t *m,
                                aocl_int64_t *n, aocl_int64_t *ncu, real *a, aocl_int64_t *lda,
                                real *s, real *u, aocl_int64_t *ldu, real *vt, aocl_int64_t *ldvt,
                                real *work, aocl_int64_t *info);
void fla_sgesvd_xs_small10T_avx2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda,
                                 real *s, real *u, aocl_int64_t *ldu, real *vt, aocl_int64_t *ldvt,
                                 real *work, aocl_int64_t *info);
doublereal fla_get_max_dabs_element_vector_avx2(aocl_int64_t m, doublereal *a, aocl_int64_t a_dim);
real fla_get_max_sabs_element_vector_avx2(aocl_int64_t m, real *a, aocl_int64_t a_dim);
doublereal fla_get_max_zabs_element_vector_avx2(aocl_int64_t m, dcomplex *a, aocl_int64_t a_dim);
real fla_get_max_cabs_element_vector_avx2(aocl_int64_t m, scomplex *a, aocl_int64_t a_dim);
void fla_dlarf_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                     aocl_int64_t ldr, doublereal *v, doublereal tau,
                                     doublereal *work);
void fla_dlarf1f_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                       aocl_int64_t ldr, doublereal *v, doublereal tau,
                                       doublereal *work);
void fla_dlarf1l_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                       aocl_int64_t ldr, doublereal *v, doublereal tau,
                                       doublereal *work);
void fla_dlarf1_right_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                       aocl_int64_t ldr, doublereal *v, doublereal *c1,
                                       doublereal ntau, doublereal *work);
void fla_dlarf1_right_update_c1_avx2(aocl_int64_t m, doublereal *restrict c1, doublereal ntau,
                                     doublereal *restrict work);
doublereal fla_dnrm2_blas_avx2(aocl_int64_t *sd, doublereal *a, aocl_int64_t *inc);
void fla_zlarf_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                     aocl_int64_t ldr, dcomplex *v, dcomplex *ntau, dcomplex *work);
void fla_zlarf1f_left_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                       aocl_int64_t ldr, dcomplex *v, dcomplex *ntau,
                                       dcomplex *work);
void fla_zlarf1_right_apply_incv1_avx2(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                       aocl_int64_t ldr, dcomplex *v, dcomplex *c1, dcomplex *ntau,
                                       dcomplex *work);
int fla_dpotrf_small_avx2(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                          aocl_int64_t *info);
void fla_dpotri_small_avx2(char *uplo, aocl_int64_t *n, double *A, aocl_int64_t *lda,
                           aocl_int64_t *info);
int fla_dcopy_scal_ix1_avx2(aocl_int64_t n, doublereal da, const doublereal *dx, doublereal *dy);
void fla_dspr_lower_avx2(aocl_int64_t len, doublereal alpha, const doublereal *x, doublereal *ap);
void fla_dlasr_left_pivotv_avx2(logical forward, aocl_int64_t m, aocl_int64_t n, doublereal *c__,
                                doublereal *s, doublereal *a, aocl_int64_t a_dim1);
#endif /* FLA_ENABLE_AMD_OPT */
#endif /* FLA_LAPACK_AVX2_KERNELS_DEFS_H */
