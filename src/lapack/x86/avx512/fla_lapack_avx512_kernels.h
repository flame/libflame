/******************************************************************************
 * Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_LAPACK_AVX512_KERNELS_DEFS_H
#define FLA_LAPACK_AVX512_KERNELS_DEFS_H
/*! @file fla_lapack_avx512_kernels.h
 *  @brief AVX512 Kernel Declarations.
 *  */
#include "FLAME.h"
#include "fla_lapack_x86_common.h"
#include "immintrin.h"

#if FLA_ENABLE_AMD_OPT

int fla_drot_avx512(aocl_int64_t *n, doublereal *dx, aocl_int64_t *incx, doublereal *dy,
                    aocl_int64_t *incy, doublereal *c__, doublereal *s);
int fla_zrot_avx512(aocl_int64_t *n, dcomplex *dx, aocl_int64_t *incx, dcomplex *dy,
                    aocl_int64_t *incy, doublereal *c__, dcomplex *s);
int fla_dhrot3_avx512(aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *v,
                      doublereal *tau);
int fla_cgetrf_small_avx512(aocl_int64_t* m, aocl_int64_t* n, scomplex* a, aocl_int64_t* lda, aocl_int_t* ipiv,
                            aocl_int64_t* info);
int fla_zgetrf_small_avx512(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                            aocl_int_t *ipiv, aocl_int64_t *info);
aocl_int64_t fla_dgetrf_small_avx512(aocl_int64_t *m, aocl_int64_t *n, doublereal *a,
                                     aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
aocl_int64_t fla_sgetrf_small_avx512(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda,
                                     aocl_int_t *ipiv, aocl_int64_t *info);
int fla_dscal_ix1_avx512(aocl_int64_t *n, doublereal *da, doublereal *dx, aocl_int64_t *incx);
int fla_sscal_ix1_avx512(aocl_int64_t *n, real *alpha, real *x);
int fla_zscal_ix1_avx512(aocl_int64_t *n, dcomplex *alpha, dcomplex *x);
doublereal fla_get_max_dabs_element_vector_avx512(aocl_int64_t m, doublereal *a,
                                                  aocl_int64_t a_dim);
real fla_get_max_sabs_element_vector_avx512(aocl_int64_t m, real *a, aocl_int64_t a_dim);
doublereal fla_get_max_zabs_element_vector_avx512(aocl_int64_t m, dcomplex *a, aocl_int64_t a_dim);
real fla_get_max_cabs_element_vector_avx512(aocl_int64_t m, scomplex *a, aocl_int64_t a_dim);
void fla_dlarf_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                       aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                       doublereal *work);
void fla_dlarf1f_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                         aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                         doublereal *work);
void fla_dlarf1l_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                         aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                         doublereal *work);
doublereal fla_dnrm2_blas_avx512(aocl_int64_t *sd, doublereal *a, aocl_int64_t *inc);
void fla_zlarf_left_apply_incv1_avx512(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff,
                                       aocl_int64_t ldr, dcomplex *v, dcomplex *ntau,
                                       dcomplex *work);
int fla_dpotrf_small_avx512(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                            aocl_int64_t *info);
void fla_dspr_lower_avx512(aocl_int64_t len, doublereal alpha, const doublereal *x, doublereal *ap);
int fla_dcopy_scal_ix1_avx512(aocl_int64_t n, doublereal da, const doublereal *dx, doublereal *dy);

#endif
#endif
