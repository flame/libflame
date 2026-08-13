/******************************************************************************
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef AOCL_FLA_LAPACK_H
#define AOCL_FLA_LAPACK_H

#include "FLA_type_defs.h"
void cdotc_f2c_(scomplex *r, const aocl_int_t *n, const scomplex *cx, const aocl_int_t *incx, const scomplex *cy,
                const aocl_int_t *incy);
void cdotu_f2c_(scomplex *r, const aocl_int_t *n, const scomplex *cx, const aocl_int_t *incx, const scomplex *cy,
                const aocl_int_t *incy);
void zdotc_f2c_(dcomplex *r, const aocl_int_t *n, const dcomplex *cx, const aocl_int_t *incx, const dcomplex *cy,
                const aocl_int_t *incy);
void zdotu_f2c_(dcomplex *r, const aocl_int_t *n, const dcomplex *cx, const aocl_int_t *incx, const dcomplex *cy,
                const aocl_int_t *incy);
void cspffrtx_fla(scomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, scomplex *work, scomplex *work2);
void dspffrt2_fla(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm, doublereal *work,
                  doublereal *work2);
void zspffrt2_fla(dcomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, dcomplex *work,
                  dcomplex *work2);
void dspffrtx_fla(doublereal *ap, aocl_int64_t *n, aocl_int64_t *ncolm, doublereal *work,
                  doublereal *work2);
void cspffrt2_fla(scomplex *ap, aocl_int64_t *n, aocl_int64_t *ncolm, scomplex *work, scomplex *work2);
void sspffrtx_fla(real *ap, aocl_int64_t *n, aocl_int64_t *ncolm, real *work, real *work2);
void sspffrt2_fla(real *ap, aocl_int64_t *n, aocl_int64_t *ncolm, real *work, real *work2);
void zhegst_fla(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, dcomplex *a,
                aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
void zhegs2_fla(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, dcomplex *a,
                aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
void chegs2_fla(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
void chegst_fla(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
void chetrd_fla(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, real *d__, real *e,
                scomplex *tau, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void zunmtr_fla(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n,
                dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *c__,
                aocl_int64_t *ldc, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void cunmtr_fla(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, scomplex *a,
                aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                aocl_int64_t *lwork, aocl_int64_t *info);
void ssytd2_fla(char *uplo, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *d__, real *e,
                real *tau, aocl_int64_t *info);
void dorgtr_fla(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *tau,
                doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
void sormtr_fla(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, real *a,
                aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                aocl_int64_t *lwork, aocl_int64_t *info);
void dsytrd_fla_parallel(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                         doublereal *d__, doublereal *e, doublereal *tau, doublereal *work,
                         aocl_int64_t *lwork, aocl_int64_t *info);;
void dsytrd_fla(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *d__,
                doublereal *e, doublereal *tau, doublereal *work, aocl_int64_t *lwork,
                aocl_int64_t *info);
void chetd2_fla(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, real *d__, real *e,
                scomplex *tau, aocl_int64_t *info);
void dormtr_fla(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n,
                doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__,
                aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
void ssytrd_fla(char *uplo, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *d__, real *e,
                real *tau, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
void dsytd2_fla(char *uplo, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *d__,
                doublereal *e, doublereal *tau, aocl_int64_t *info);
void zungtr_fla(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                dcomplex *tau, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void zhetd2_fla(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, doublereal *d__,
                doublereal *e, dcomplex *tau, aocl_int64_t *info);
void cungtr_fla(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void zhetrd_fla(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, doublereal *d__,
                doublereal *e, dcomplex *tau, dcomplex *work, aocl_int64_t *lwork,
                aocl_int64_t *info);
void sorgtr_fla(char *uplo, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau, real *work,
                aocl_int64_t *lwork, aocl_int64_t *info);
void sormqr_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a,
                aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                aocl_int64_t *lwork, aocl_int64_t *info);
void zunmqr_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *c__,
                aocl_int64_t *ldc, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void cung2r_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda,
                scomplex *tau, scomplex *work, aocl_int64_t *info);
void sorg2r_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a, aocl_int64_t *lda,
                real *tau, real *work, aocl_int64_t *info);
void cungqr_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda,
                scomplex *tau, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void sgeqr2_fla(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau, real *work,
                aocl_int64_t *info);
void dgeqr2p_fla(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                 doublereal *tau, doublereal *work, aocl_int64_t *info);
void zungqr_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                aocl_int64_t *lda, dcomplex *tau, dcomplex *work, aocl_int64_t *lwork,
                aocl_int64_t *info);
void dgeqp3_fla(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                aocl_int_t *jpvt, doublereal *tau, doublereal *work, aocl_int64_t *lwork,
                aocl_int64_t *info);
void zunm2r_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *c__,
                aocl_int64_t *ldc, dcomplex *work, aocl_int64_t *info);
void cunm2r_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc,
                scomplex *work, aocl_int64_t *info);
void sgeqrfp_fla(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau,
                 real *work, aocl_int64_t *lwork, aocl_int64_t *info);
void dormbr_fla(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n,
                aocl_int64_t *k, doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__,
                aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
void cunmqr_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc,
                scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void dorm2r_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__,
                aocl_int64_t *ldc, doublereal *work, aocl_int64_t *info);
void dgeqrf_fla(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *tau,
                doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
void sormbr_fla(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n,
                aocl_int64_t *k, real *a, aocl_int64_t *lda, real *tau, real *c__,
                aocl_int64_t *ldc, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
void zung2r_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                aocl_int64_t *lda, dcomplex *tau, dcomplex *work, aocl_int64_t *info);
void sgeqpf_fla(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, aocl_int_t *jpvt,
                real *tau, real *work, aocl_int64_t *info);
void dgeqpf_fla(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                aocl_int_t *jpvt, doublereal *tau, doublereal *work, aocl_int64_t *info);
void dgeqr2_fla(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *tau,
                doublereal *work, aocl_int64_t *info);
void dgeqrfp_fla(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                 doublereal *tau, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
void dormqr_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__,
                aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
void sgeqr2p_fla(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau,
                 real *work, aocl_int64_t *info);
void sgeqp3_fla(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, aocl_int_t *jpvt,
                real *tau, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
void sorgqr_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a, aocl_int64_t *lda,
                real *tau, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
void dorg2r_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a, aocl_int64_t *lda,
                doublereal *tau, doublereal *work, aocl_int64_t *info);
void sorm2r_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a,
                aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                aocl_int64_t *info);
void sgeqrf_fla(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau, real *work,
                aocl_int64_t *lwork, aocl_int64_t *info);
void zungl2_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                aocl_int64_t *lda, dcomplex *tau, dcomplex *work, aocl_int64_t *info);
void cunml2_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc,
                scomplex *work, aocl_int64_t *info);
void dorml2_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__,
                aocl_int64_t *ldc, doublereal *work, aocl_int64_t *info);
void cungl2_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda,
                scomplex *tau, scomplex *work, aocl_int64_t *info);
void dorgl2_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a, aocl_int64_t *lda,
                doublereal *tau, doublereal *work, aocl_int64_t *info);
void sorglq_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a, aocl_int64_t *lda,
                real *tau, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
void sorgl2_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a, aocl_int64_t *lda,
                real *tau, real *work, aocl_int64_t *info);
void dorglq_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a, aocl_int64_t *lda,
                doublereal *tau, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
void zunml2_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *c__,
                aocl_int64_t *ldc, dcomplex *work, aocl_int64_t *info);
void cunglq_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda,
                scomplex *tau, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void cunmlq_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc,
                scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void sorml2_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a,
                aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                aocl_int64_t *info);
void dormlq_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__,
                aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
void zunmlq_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *c__,
                aocl_int64_t *ldc, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
void sormlq_fla(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a,
                aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                aocl_int64_t *lwork, aocl_int64_t *info);
void zunglq_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                aocl_int64_t *lda, dcomplex *tau, dcomplex *work, aocl_int64_t *lwork,
                aocl_int64_t *info);
int lapack_dormbr(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n,
                  aocl_int64_t *k, doublereal *a, aocl_int64_t *lda, doublereal *tau,
                  doublereal *c__, aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork,
                  aocl_int64_t *info);
int lapack_dorml2(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                  doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__,
                  aocl_int64_t *ldc, doublereal *work, aocl_int64_t *info);
int lapack_sorgl2(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a, aocl_int64_t *lda,
                  real *tau, real *work, aocl_int64_t *info);
int lapack_sbdsqr(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru,
                  aocl_int64_t *ncc, real *d__, real *e, real *vt, aocl_int64_t *ldvt, real *u,
                  aocl_int64_t *ldu, real *c__, aocl_int64_t *ldc, real *work, aocl_int64_t *info);
int lapack_dorg2r(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a,
                  aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *info);
int lapack_sgesvd(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, real *a,
                  aocl_int64_t *lda, real *s, real *u, aocl_int64_t *ldu, real *vt,
                  aocl_int64_t *ldvt, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dorgqr(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a,
                  aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *lwork,
                  aocl_int64_t *info);
int lapack_dormqr(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                  doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__,
                  aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_sgelqf(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau,
                  real *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_sorglq(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a, aocl_int64_t *lda,
                  real *tau, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dormlq(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                  doublereal *a, aocl_int64_t *lda, doublereal *tau, doublereal *c__,
                  aocl_int64_t *ldc, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dbdsqr(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru,
                  aocl_int64_t *ncc, doublereal *d__, doublereal *e, doublereal *vt,
                  aocl_int64_t *ldvt, doublereal *u, aocl_int64_t *ldu, doublereal *c__,
                  aocl_int64_t *ldc, doublereal *work, aocl_int64_t *info);
int lapack_sorgbr(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a,
                  aocl_int64_t *lda, real *tau, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dorgl2(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a,
                  aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *info);
int lapack_dgebd2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                  doublereal *d__, doublereal *e, doublereal *tauq, doublereal *taup,
                  doublereal *work, aocl_int64_t *info);
int lapack_sgebrd(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *d__, real *e,
                  real *tauq, real *taup, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_sormlq(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                  real *a, aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dorgbr(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a,
                  aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *lwork,
                  aocl_int64_t *info);
int lapack_dgebrd(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                  doublereal *d__, doublereal *e, doublereal *tauq, doublereal *taup,
                  doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dgesvd(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, doublereal *a,
                  aocl_int64_t *lda, doublereal *s, doublereal *u, aocl_int64_t *ldu,
                  doublereal *vt, aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *lwork,
                  aocl_int64_t *info);
int lapack_sgesdd(char *jobz, aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *s,
                  real *u, aocl_int64_t *ldu, real *vt, aocl_int64_t *ldvt, real *work,
                  aocl_int64_t *lwork, aocl_int_t *iwork, aocl_int64_t *info);
int lapack_dgelqf(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                  doublereal *tau, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_dbdsqr_small(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru,
                        doublereal *d__, doublereal *e, doublereal *vt, aocl_int64_t *ldvt,
                        doublereal *u, aocl_int64_t *ldu, aocl_int64_t *info);
int lapack_sbdsqr_small(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru,
                        real *d__, real *e, real *vt, aocl_int64_t *ldvt, real *u,
                        aocl_int64_t *ldu, aocl_int64_t *info);
int lapack_sorml2(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                  real *a, aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                  aocl_int64_t *info);
int lapack_sorgqr(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a, aocl_int64_t *lda,
                  real *tau, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_sormbr(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n,
                  aocl_int64_t *k, real *a, aocl_int64_t *lda, real *tau, real *c__,
                  aocl_int64_t *ldc, real *work, aocl_int64_t *lwork, aocl_int64_t *info);
int lapack_sgelq2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau,
                  real *work, aocl_int64_t *info);
int lapack_dorglq(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a,
                  aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *lwork,
                  aocl_int64_t *info);
int lapack_sorm2r(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                  real *a, aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                  aocl_int64_t *info);
int lapack_dgesdd(char *jobz, aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                  doublereal *s, doublereal *u, aocl_int64_t *ldu, doublereal *vt,
                  aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *lwork, aocl_int_t *iwork,
                  aocl_int64_t *info);
int lapack_dgelq2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                  doublereal *tau, doublereal *work, aocl_int64_t *info);
int lapack_sorg2r(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a, aocl_int64_t *lda,
                  real *tau, real *work, aocl_int64_t *info);
int lapack_sgebd2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *d__, real *e,
                  real *tauq, real *taup, real *work, aocl_int64_t *info);
int lapack_sormqr(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                  real *a, aocl_int64_t *lda, real *tau, real *c__, aocl_int64_t *ldc, real *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
fla_dim_t lapack_cgetf2(fla_dim_t *m, fla_dim_t *n, scomplex *a, fla_dim_t *lda,
	 aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_cgetrf(fla_dim_t *m, fla_dim_t *n, scomplex *a, fla_dim_t *lda,
	 aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_sgetf2(fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda,
	aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_sgetrf(fla_dim_t *m, fla_dim_t *n, real *a, fla_dim_t *lda,
	aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_zgetf2(fla_dim_t *m, fla_dim_t *n, dcomplex *a,
	fla_dim_t *lda, aocl_int_t *ipiv, fla_dim_t *info);
fla_dim_t lapack_zgetrf(fla_dim_t *m, fla_dim_t *n, dcomplex *a,
	fla_dim_t *lda, aocl_int_t *ipiv, fla_dim_t *info);
void fla_dscal(aocl_int64_t *n, doublereal *da, doublereal *dx, aocl_int64_t *incx);
aocl_int64_t fla_i_nint(real *x);

#endif /* defined(AOCL_FLA_LAPACK_H) */