/******************************************************************************
 * * Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

#ifndef FLA_LAPACK_X86_COMMON_H
#define FLA_LAPACK_X86_COMMON_H

/*! @file fla_lapack_x86_common.h
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"

#if FLA_ENABLE_AMD_OPT

void fla_dtranspose(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                    doublereal *b, aocl_int64_t *ldb);
void fla_stranspose(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *b,
                    aocl_int64_t *ldb);
int fla_dhrot3(aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *v, doublereal *tau);
int fla_drot(aocl_int64_t *n, doublereal *dx, aocl_int64_t *incx, doublereal *dy,
             aocl_int64_t *incy, doublereal *c__, doublereal *s);
int fla_zscal(aocl_int64_t *n, dcomplex *alpha, dcomplex *x, aocl_int64_t *incx);
int fla_dgeqrf_small(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                     doublereal *tau, doublereal *work);
int fla_sgeqrf_small(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau,
                     real *work);
int fla_cgeqrf_small(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                     scomplex *work);
int fla_zgeqrf_small(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                     dcomplex *work);
int fla_dgelqf_small(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                     doublereal *tau, doublereal *work);
int fla_sgelqf_small(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau,
                     real *work);
void fla_dscal(aocl_int64_t *n, doublereal *alpha, doublereal *x, aocl_int64_t *incx);
void fla_sscal(aocl_int64_t *n, real *alpha, real *x, aocl_int64_t *incx);
void fla_sger(aocl_int64_t *m, aocl_int64_t *n, real *alpha, real *x, aocl_int64_t *incx, real *y,
              aocl_int64_t *incy, real *a, aocl_int64_t *lda);
int fla_cgetrf_small_simd(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                          aocl_int64_t *info);
int fla_zgetrf_small_simd(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                          aocl_int_t *ipiv, aocl_int64_t *info);
void fla_dgesvd_small6(aocl_int64_t wntus, aocl_int64_t wntvs, aocl_int64_t *m, aocl_int64_t *n,
                       doublereal *a, aocl_int64_t *lda, doublereal *qr, aocl_int64_t *ldqr,
                       doublereal *s, doublereal *u, aocl_int64_t *ldu, doublereal *vt,
                       aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *info);
void fla_dgesvd_nn_small1T(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                           doublereal *s, doublereal *work, aocl_int64_t *info);
void fla_dgesvd_small6T(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                        doublereal *ql, aocl_int64_t *ldql, doublereal *s, doublereal *u,
                        aocl_int64_t *ldu, doublereal *vt, aocl_int64_t *ldvt, doublereal *work,
                        aocl_int64_t *info);
void fla_dgesvd_xx_small10(aocl_int64_t wntu, aocl_int64_t wntv, aocl_int64_t *m, aocl_int64_t *n,
                           aocl_int64_t *ncu, doublereal *a, aocl_int64_t *lda, doublereal *s,
                           doublereal *u, aocl_int64_t *ldu, doublereal *vt, aocl_int64_t *ldvt,
                           doublereal *work, aocl_int64_t *info);
void fla_dgesvd_xs_small10T(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                            doublereal *s, doublereal *u, aocl_int64_t *ldu, doublereal *vt,
                            aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *info);
void fla_sgesvd_small6(aocl_int64_t wntus, aocl_int64_t wntvs, aocl_int64_t *m, aocl_int64_t *n,
                       real *a, aocl_int64_t *lda, real *qr, aocl_int64_t *ldqr, real *s, real *u,
                       aocl_int64_t *ldu, real *vt, aocl_int64_t *ldvt, real *work,
                       aocl_int64_t *info);
void fla_sgesvd_nn_small1T(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *s,
                           real *work, aocl_int64_t *info);
void fla_sgesvd_small6T(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *ql,
                        aocl_int64_t *ldql, real *s, real *u, aocl_int64_t *ldu, real *vt,
                        aocl_int64_t *ldvt, real *work, aocl_int64_t *info);
void fla_sgesvd_xx_small10(aocl_int64_t wntu, aocl_int64_t wntv, aocl_int64_t *m, aocl_int64_t *n,
                           aocl_int64_t *ncu, real *a, aocl_int64_t *lda, real *s, real *u,
                           aocl_int64_t *ldu, real *vt, aocl_int64_t *ldvt, real *work,
                           aocl_int64_t *info);
void fla_sgesvd_xs_small10T(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *s,
                            real *u, aocl_int64_t *ldu, real *vt, aocl_int64_t *ldvt, real *work,
                            aocl_int64_t *info);
int fla_dgetrs_small_notrans(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *a,
                             aocl_int64_t *lda, aocl_int_t *ipiv, doublereal *b, aocl_int64_t *ldb,
                             aocl_int64_t *info);
void lapack_getri_small_d(aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                          doublereal *work, aocl_int64_t *info);
doublereal fla_dnrm2_blas_kernel(aocl_int64_t *sd, doublereal *a, aocl_int64_t *inc);
real fla_get_max_sabs_element_vector(aocl_int64_t m, real *a, aocl_int64_t a_dim);
doublereal fla_get_max_dabs_element_vector(aocl_int64_t m, doublereal *a, aocl_int64_t a_dim);
real fla_get_max_cabs_element_vector(aocl_int64_t m, scomplex *a, aocl_int64_t a_dim);
doublereal fla_get_max_zabs_element_vector(aocl_int64_t m, dcomplex *a, aocl_int64_t a_dim);
void fla_zlarf_left_invc1_opt(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff, aocl_int64_t ldr,
                              dcomplex *v, dcomplex *ntau, dcomplex *work);
void fla_dspr_lower(aocl_int64_t len, doublereal alpha, const doublereal *x, doublereal *ap);
void fla_dcopy_scal(aocl_int64_t n, doublereal alpha, const doublereal *x, aocl_int64_t incx,
                    doublereal *y, aocl_int64_t incy);
void fla_dlarf_small_incv1_simd(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                aocl_int64_t ldr, doublereal *v, doublereal ntau, doublereal *work);
void fla_dlarf1f_small_incv1_simd(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                  aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                  doublereal *work);
void fla_dlarf1l_small_incv1_simd(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                  aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                  doublereal *work);
void fla_dlasr_left_pivotv_simd(logical forward, aocl_int64_t m, aocl_int64_t n, doublereal *c__,
                                doublereal *s, doublereal *a, aocl_int64_t a_dim1);
void fla_dlarf1_small_right_incv1_simd(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                       aocl_int64_t ldr, doublereal *v, doublereal *c1,
                                       doublereal ntau, doublereal *work);
void fla_dlarf1_right_update_c1_simd(aocl_int64_t m, doublereal *c1, doublereal ntau,
                                     doublereal *work);
#endif

#endif /* FLA_LAPACK_X86_COMMON_H */
