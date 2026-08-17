/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
#ifndef BLAS_INTERFACE_H
#define BLAS_INTERFACE_H

#include "FLA_f2c.h"

void caxpy_(const aocl_int_t *n, const scomplex *ca, const scomplex *cx, const aocl_int_t *incx,
            scomplex *cy, const aocl_int_t *incy);
void ccopy_(const aocl_int_t *n, const scomplex *cx, const aocl_int_t *incx, scomplex *cy,
            const aocl_int_t *incy);
#ifndef FLA_ENABLE_F2C_DOTC
scomplex cdotc_(const aocl_int_t *n, const scomplex *cx, const aocl_int_t *incx, const scomplex *cy,
               const aocl_int_t *incy);
scomplex cdotu_(const aocl_int_t *n, const scomplex *cx, const aocl_int_t *incx, const scomplex *cy,
               const aocl_int_t *incy);
#else
void cdotc_(scomplex *r, const aocl_int_t *n, const scomplex *cx, const aocl_int_t *incx,
            const scomplex *cy, const aocl_int_t *incy);
void cdotu_(scomplex *r, const aocl_int_t *n, const scomplex *cx, const aocl_int_t *incx,
            const scomplex *cy, const aocl_int_t *incy);
#endif
void cgbmv_(const char *trans, const aocl_int_t *m, const aocl_int_t *n, const aocl_int_t *kl,
            const aocl_int_t *ku, const scomplex *alpha, const scomplex *a, const aocl_int_t *lda,
            const scomplex *x, const aocl_int_t *incx, const scomplex *beta, scomplex *y,
            const aocl_int_t *incy);
void cgemm_(const char *transa, const char *transb, const aocl_int_t *m, const aocl_int_t *n,
            const aocl_int_t *k, const scomplex *alpha, const scomplex *a, const aocl_int_t *lda,
            const scomplex *b, const aocl_int_t *ldb, const scomplex *beta, scomplex *c,
            const aocl_int_t *ldc);
void cgemmtr_(const char *uplo, const char *transa, const char *transb, const aocl_int_t *n,
              const aocl_int_t *k, const scomplex *alpha, const scomplex *a, const aocl_int_t *lda,
              const scomplex *b, const aocl_int_t *ldb, const scomplex *beta, scomplex *c,
              const aocl_int_t *ldc);
void cgemv_(const char *trans, const aocl_int_t *m, const aocl_int_t *n, const scomplex *alpha,
            const scomplex *a, const aocl_int_t *lda, const scomplex *x, const aocl_int_t *incx,
            const scomplex *beta, scomplex *y, const aocl_int_t *incy);
void cgerc_(const aocl_int_t *m, const aocl_int_t *n, const scomplex *alpha, const scomplex *x,
            const aocl_int_t *incx, const scomplex *y, const aocl_int_t *incy, scomplex *a,
            const aocl_int_t *lda);
void cgeru_(const aocl_int_t *m, const aocl_int_t *n, const scomplex *alpha, const scomplex *x,
            const aocl_int_t *incx, const scomplex *y, const aocl_int_t *incy, scomplex *a,
            const aocl_int_t *lda);
void chbmv_(const char *uplo, const aocl_int_t *n, const aocl_int_t *k, const scomplex *alpha,
            const scomplex *a, const aocl_int_t *lda, const scomplex *x, const aocl_int_t *incx,
            const scomplex *beta, scomplex *y, const aocl_int_t *incy);
void chemm_(const char *side, const char *uplo, const aocl_int_t *m, const aocl_int_t *n,
            const scomplex *alpha, const scomplex *a, const aocl_int_t *lda, const scomplex *b,
            const aocl_int_t *ldb, const scomplex *beta, scomplex *c, const aocl_int_t *ldc);
void chemv_(const char *uplo, const aocl_int_t *n, const scomplex *alpha, const scomplex *a,
            const aocl_int_t *lda, const scomplex *x, const aocl_int_t *incx, const scomplex *beta,
            scomplex *y, const aocl_int_t *incy);
void cher_(const char *uplo, const aocl_int_t *n, const real *alpha, const scomplex *x,
           const aocl_int_t *incx, scomplex *a, const aocl_int_t *lda);
void cher2_(const char *uplo, const aocl_int_t *n, const scomplex *alpha, const scomplex *x,
            const aocl_int_t *incx, const scomplex *y, const aocl_int_t *incy, scomplex *a,
            const aocl_int_t *lda);
void cher2k_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
             const scomplex *alpha, const scomplex *a, const aocl_int_t *lda, const scomplex *b,
             const aocl_int_t *ldb, const real *beta, scomplex *c, const aocl_int_t *ldc);
void cherk_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
            const real *alpha, const scomplex *a, const aocl_int_t *lda, const real *beta,
            scomplex *c, const aocl_int_t *ldc);
void chpmv_(const char *uplo, const aocl_int_t *n, const scomplex *alpha, const scomplex *ap,
            const scomplex *x, const aocl_int_t *incx, const scomplex *beta, scomplex *y,
            const aocl_int_t *incy);
void chpr_(const char *uplo, const aocl_int_t *n, const real *alpha, const scomplex *x,
           const aocl_int_t *incx, scomplex *ap);
void chpr2_(const char *uplo, const aocl_int_t *n, const scomplex *alpha, const scomplex *x,
            const aocl_int_t *incx, const scomplex *y, const aocl_int_t *incy, scomplex *ap);
void cscal_(const aocl_int_t *n, const scomplex *ca, scomplex *cx, const aocl_int_t *incx);
void csrot_(const aocl_int_t *n, scomplex *cx, const aocl_int_t *incx, scomplex *cy,
            const aocl_int_t *incy, const real *c, const real *s);
void csscal_(const aocl_int_t *n, const real *sa, scomplex *cx, const aocl_int_t *incx);
void cswap_(const aocl_int_t *n, scomplex *cx, const aocl_int_t *incx, scomplex *cy,
            const aocl_int_t *incy);
void csymm_(const char *side, const char *uplo, const aocl_int_t *m, const aocl_int_t *n,
            const scomplex *alpha, const scomplex *a, const aocl_int_t *lda, const scomplex *b,
            const aocl_int_t *ldb, const scomplex *beta, scomplex *c, const aocl_int_t *ldc);
void csyr2k_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
             const scomplex *alpha, const scomplex *a, const aocl_int_t *lda, const scomplex *b,
             const aocl_int_t *ldb, const scomplex *beta, scomplex *c, const aocl_int_t *ldc);
void csyrk_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
            const scomplex *alpha, const scomplex *a, const aocl_int_t *lda, const scomplex *beta,
            scomplex *c, const aocl_int_t *ldc);
void ctbmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const aocl_int_t *k, const scomplex *a, const aocl_int_t *lda, scomplex *x,
            const aocl_int_t *incx);
void ctbsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const aocl_int_t *k, const scomplex *a, const aocl_int_t *lda, scomplex *x,
            const aocl_int_t *incx);
void ctpmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const scomplex *ap, scomplex *x, const aocl_int_t *incx);
void ctpsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const scomplex *ap, scomplex *x, const aocl_int_t *incx);
void ctrmm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const aocl_int_t *m, const aocl_int_t *n, const scomplex *alpha, const scomplex *a,
            const aocl_int_t *lda, scomplex *b, const aocl_int_t *ldb);
void ctrmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const scomplex *a, const aocl_int_t *lda, scomplex *x, const aocl_int_t *incx);
void ctrsm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const aocl_int_t *m, const aocl_int_t *n, const scomplex *alpha, const scomplex *a,
            const aocl_int_t *lda, scomplex *b, const aocl_int_t *ldb);
void ctrsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const scomplex *a, const aocl_int_t *lda, scomplex *x, const aocl_int_t *incx);
doublereal dasum_(const aocl_int_t *n, const doublereal *dx, const aocl_int_t *incx);
void daxpy_(const aocl_int_t *n, const doublereal *da, const doublereal *dx, const aocl_int_t *incx,
            doublereal *dy, const aocl_int_t *incy);
void dcopy_(const aocl_int_t *n, const doublereal *dx, const aocl_int_t *incx, doublereal *dy,
            const aocl_int_t *incy);
doublereal ddot_(const aocl_int_t *n, const doublereal *dx, const aocl_int_t *incx,
                 const doublereal *dy, const aocl_int_t *incy);
void dgbmv_(const char *trans, const aocl_int_t *m, const aocl_int_t *n, const aocl_int_t *kl,
            const aocl_int_t *ku, const doublereal *alpha, const doublereal *a,
            const aocl_int_t *lda, const doublereal *x, const aocl_int_t *incx,
            const doublereal *beta, doublereal *y, const aocl_int_t *incy);
void dgemm_(const char *transa, const char *transb, const aocl_int_t *m, const aocl_int_t *n,
            const aocl_int_t *k, const doublereal *alpha, const doublereal *a,
            const aocl_int_t *lda, const doublereal *b, const aocl_int_t *ldb,
            const doublereal *beta, doublereal *c, const aocl_int_t *ldc);
void dgemmtr_(const char *uplo, const char *transa, const char *transb, const aocl_int_t *n,
              const aocl_int_t *k, const doublereal *alpha, const doublereal *a,
              const aocl_int_t *lda, const doublereal *b, const aocl_int_t *ldb,
              const doublereal *beta, doublereal *c, const aocl_int_t *ldc);
void dgemv_(const char *trans, const aocl_int_t *m, const aocl_int_t *n, const doublereal *alpha,
            const doublereal *a, const aocl_int_t *lda, const doublereal *x, const aocl_int_t *incx,
            const doublereal *beta, doublereal *y, const aocl_int_t *incy);
void dger_(const aocl_int_t *m, const aocl_int_t *n, const doublereal *alpha, const doublereal *x,
           const aocl_int_t *incx, const doublereal *y, const aocl_int_t *incy, doublereal *a,
           const aocl_int_t *lda);
doublereal dnrm2_(const aocl_int_t *n, const doublereal *x, const aocl_int_t *incx);
void drot_(const aocl_int_t *n, doublereal *dx, const aocl_int_t *incx, doublereal *dy,
           const aocl_int_t *incy, const doublereal *c, const doublereal *s);
void drotm_(const aocl_int_t *n, doublereal *dx, const aocl_int_t *incx, doublereal *dy,
            const aocl_int_t *incy, const doublereal *dparam);
void dsbmv_(const char *uplo, const aocl_int_t *n, const aocl_int_t *k, const doublereal *alpha,
            const doublereal *a, const aocl_int_t *lda, const doublereal *x, const aocl_int_t *incx,
            const doublereal *beta, doublereal *y, const aocl_int_t *incy);
void dscal_(const aocl_int_t *n, const doublereal *da, doublereal *dx, const aocl_int_t *incx);
doublereal dsdot_(const aocl_int_t *n, const real *sx, const aocl_int_t *incx, const real *sy,
                  const aocl_int_t *incy);
void dspmv_(const char *uplo, const aocl_int_t *n, const doublereal *alpha, const doublereal *ap,
            const doublereal *x, const aocl_int_t *incx, const doublereal *beta, doublereal *y,
            const aocl_int_t *incy);
void dspr_(const char *uplo, const aocl_int_t *n, const doublereal *alpha, const doublereal *x,
           const aocl_int_t *incx, doublereal *ap);
void dspr2_(const char *uplo, const aocl_int_t *n, const doublereal *alpha, const doublereal *x,
            const aocl_int_t *incx, const doublereal *y, const aocl_int_t *incy, doublereal *ap);
void dswap_(const aocl_int_t *n, doublereal *dx, const aocl_int_t *incx, doublereal *dy,
            const aocl_int_t *incy);
void dsymm_(const char *side, const char *uplo, const aocl_int_t *m, const aocl_int_t *n,
            const doublereal *alpha, const doublereal *a, const aocl_int_t *lda,
            const doublereal *b, const aocl_int_t *ldb, const doublereal *beta, doublereal *c,
            const aocl_int_t *ldc);
void dsymv_(const char *uplo, const aocl_int_t *n, const doublereal *alpha, const doublereal *a,
            const aocl_int_t *lda, const doublereal *x, const aocl_int_t *incx,
            const doublereal *beta, doublereal *y, const aocl_int_t *incy);
void dsyr_(const char *uplo, const aocl_int_t *n, const doublereal *alpha, const doublereal *x,
           const aocl_int_t *incx, doublereal *a, const aocl_int_t *lda);
void dsyr2_(const char *uplo, const aocl_int_t *n, const doublereal *alpha, const doublereal *x,
            const aocl_int_t *incx, const doublereal *y, const aocl_int_t *incy, doublereal *a,
            const aocl_int_t *lda);
void dsyr2k_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
             const doublereal *alpha, const doublereal *a, const aocl_int_t *lda,
             const doublereal *b, const aocl_int_t *ldb, const doublereal *beta, doublereal *c,
             const aocl_int_t *ldc);
void dsyrk_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
            const doublereal *alpha, const doublereal *a, const aocl_int_t *lda,
            const doublereal *beta, doublereal *c, const aocl_int_t *ldc);
void dtbmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const aocl_int_t *k, const doublereal *a, const aocl_int_t *lda, doublereal *x,
            const aocl_int_t *incx);
void dtbsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const aocl_int_t *k, const doublereal *a, const aocl_int_t *lda, doublereal *x,
            const aocl_int_t *incx);
void dtpmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const doublereal *ap, doublereal *x, const aocl_int_t *incx);
void dtpsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const doublereal *ap, doublereal *x, const aocl_int_t *incx);
void dtrmm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const aocl_int_t *m, const aocl_int_t *n, const doublereal *alpha, const doublereal *a,
            const aocl_int_t *lda, doublereal *b, const aocl_int_t *ldb);
void dtrmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const doublereal *a, const aocl_int_t *lda, doublereal *x, const aocl_int_t *incx);
void dtrsm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const aocl_int_t *m, const aocl_int_t *n, const doublereal *alpha, const doublereal *a,
            const aocl_int_t *lda, doublereal *b, const aocl_int_t *ldb);
void dtrsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const doublereal *a, const aocl_int_t *lda, doublereal *x, const aocl_int_t *incx);
doublereal dzasum_(const aocl_int_t *n, dcomplex *zx, const aocl_int_t *incx);
doublereal dznrm2_(const aocl_int_t *n, const dcomplex *x, const aocl_int_t *incx);
aocl_int_t icamax_(const aocl_int_t *n, const scomplex *x, const aocl_int_t *incx);
aocl_int_t idamax_(const aocl_int_t *n, const doublereal *dx, const aocl_int_t *incx);
aocl_int_t isamax_(const aocl_int_t *n, const real *sx, const aocl_int_t *incx);
aocl_int_t izamax_(const aocl_int_t *n, const dcomplex *x, const aocl_int_t *incx);
real sasum_(const aocl_int_t *n, const real *sx, const aocl_int_t *incx);
void saxpy_(const aocl_int_t *n, const real *sa, const real *sx, const aocl_int_t *incx, real *sy,
            const aocl_int_t *incy);
real scasum_(const aocl_int_t *n, scomplex *cx, const aocl_int_t *incx);
real scnrm2_(const aocl_int_t *n, const scomplex *x, const aocl_int_t *incx);
void scopy_(const aocl_int_t *n, const real *sx, const aocl_int_t *incx, real *sy,
            const aocl_int_t *incy);
real sdot_(const aocl_int_t *n, const real *sx, const aocl_int_t *incx, const real *sy,
           const aocl_int_t *incy);
real sdsdot_(const aocl_int_t *n, const real *sb, const real *sx, const aocl_int_t *incx,
             const real *sy, const aocl_int_t *incy);
void sgbmv_(const char *trans, const aocl_int_t *m, const aocl_int_t *n, const aocl_int_t *kl,
            const aocl_int_t *ku, const real *alpha, const real *a, const aocl_int_t *lda,
            const real *x, const aocl_int_t *incx, const real *beta, real *y,
            const aocl_int_t *incy);
void sgemm_(const char *transa, const char *transb, const aocl_int_t *m, const aocl_int_t *n,
            const aocl_int_t *k, const real *alpha, const real *a, const aocl_int_t *lda,
            const real *b, const aocl_int_t *ldb, const real *beta, real *c, const aocl_int_t *ldc);
void sgemmtr_(const char *uplo, const char *transa, const char *transb, const aocl_int_t *n,
              const aocl_int_t *k, const real *alpha, const real *a, const aocl_int_t *lda,
              const real *b, const aocl_int_t *ldb, const real *beta, real *c,
              const aocl_int_t *ldc);
void sgemv_(const char *trans, const aocl_int_t *m, const aocl_int_t *n, const real *alpha,
            const real *a, const aocl_int_t *lda, const real *x, const aocl_int_t *incx,
            const real *beta, real *y, const aocl_int_t *incy);
void sger_(const aocl_int_t *m, const aocl_int_t *n, const real *alpha, const real *x,
           const aocl_int_t *incx, const real *y, const aocl_int_t *incy, real *a,
           const aocl_int_t *lda);
real snrm2_(const aocl_int_t *n, const real *x, const aocl_int_t *incx);
void srot_(const aocl_int_t *n, real *sx, const aocl_int_t *incx, real *sy, const aocl_int_t *incy,
           const real *c, const real *s);
void srotm_(const aocl_int_t *n, real *sx, const aocl_int_t *incx, real *sy, const aocl_int_t *incy,
            const real *sparam);
void ssbmv_(const char *uplo, const aocl_int_t *n, const aocl_int_t *k, const real *alpha,
            const real *a, const aocl_int_t *lda, const real *x, const aocl_int_t *incx,
            const real *beta, real *y, const aocl_int_t *incy);
void sscal_(const aocl_int_t *n, const real *sa, real *sx, const aocl_int_t *incx);
void sspmv_(const char *uplo, const aocl_int_t *n, const real *alpha, const real *ap, const real *x,
            const aocl_int_t *incx, const real *beta, real *y, const aocl_int_t *incy);
void sspr_(const char *uplo, const aocl_int_t *n, const real *alpha, const real *x,
           const aocl_int_t *incx, real *ap);
void sspr2_(const char *uplo, const aocl_int_t *n, const real *alpha, const real *x,
            const aocl_int_t *incx, const real *y, const aocl_int_t *incy, real *ap);
void sswap_(const aocl_int_t *n, real *sx, const aocl_int_t *incx, real *sy,
            const aocl_int_t *incy);
void ssymm_(const char *side, const char *uplo, const aocl_int_t *m, const aocl_int_t *n,
            const real *alpha, const real *a, const aocl_int_t *lda, const real *b,
            const aocl_int_t *ldb, const real *beta, real *c, const aocl_int_t *ldc);
void ssymv_(const char *uplo, const aocl_int_t *n, const real *alpha, const real *a,
            const aocl_int_t *lda, const real *x, const aocl_int_t *incx, const real *beta, real *y,
            const aocl_int_t *incy);
void ssyr_(const char *uplo, const aocl_int_t *n, const real *alpha, const real *x,
           const aocl_int_t *incx, real *a, const aocl_int_t *lda);
void ssyr2_(const char *uplo, const aocl_int_t *n, const real *alpha, const real *x,
            const aocl_int_t *incx, const real *y, const aocl_int_t *incy, real *a,
            const aocl_int_t *lda);
void ssyr2k_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
             const real *alpha, const real *a, const aocl_int_t *lda, const real *b,
             const aocl_int_t *ldb, const real *beta, real *c, const aocl_int_t *ldc);
void ssyrk_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
            const real *alpha, const real *a, const aocl_int_t *lda, const real *beta, real *c,
            const aocl_int_t *ldc);
void stbmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const aocl_int_t *k, const real *a, const aocl_int_t *lda, real *x,
            const aocl_int_t *incx);
void stbsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const aocl_int_t *k, const real *a, const aocl_int_t *lda, real *x,
            const aocl_int_t *incx);
void stpmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const real *ap, real *x, const aocl_int_t *incx);
void stpsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const real *ap, real *x, const aocl_int_t *incx);
void strmm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const aocl_int_t *m, const aocl_int_t *n, const real *alpha, const real *a,
            const aocl_int_t *lda, real *b, const aocl_int_t *ldb);
void strmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const real *a, const aocl_int_t *lda, real *x, const aocl_int_t *incx);
void strsm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const aocl_int_t *m, const aocl_int_t *n, const real *alpha, const real *a,
            const aocl_int_t *lda, real *b, const aocl_int_t *ldb);
void strsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const real *a, const aocl_int_t *lda, real *x, const aocl_int_t *incx);
void xerbla_(const char *srname, const aocl_int_t *info, aocl_int_t srname_len);
void xerbla_array_(const char *srname_array, const aocl_int_t *srname_len, const aocl_int_t *info);
void zaxpy_(const aocl_int_t *n, const dcomplex *za, const dcomplex *zx,
            const aocl_int_t *incx, dcomplex *zy, const aocl_int_t *incy);
void zcopy_(const aocl_int_t *n, const dcomplex *zx, const aocl_int_t *incx, dcomplex *zy,
            const aocl_int_t *incy);
#ifndef FLA_ENABLE_F2C_DOTC
dcomplex zdotc_(const aocl_int_t *n, const dcomplex *zx, const aocl_int_t *incx,
                     const dcomplex *zy, const aocl_int_t *incy);
dcomplex zdotu_(const aocl_int_t *n, const dcomplex *zx, const aocl_int_t *incx,
                     const dcomplex *zy, const aocl_int_t *incy);
#else
void zdotc_(dcomplex *r, const aocl_int_t *n, const dcomplex *zx, const aocl_int_t *incx,
            const dcomplex *zy, const aocl_int_t *incy);
void zdotu_(dcomplex *r, const aocl_int_t *n, const dcomplex *zx, const aocl_int_t *incx,
            const dcomplex *zy, const aocl_int_t *incy);
#endif
void zdrot_(const aocl_int_t *n, dcomplex *zx, const aocl_int_t *incx, dcomplex *zy,
            const aocl_int_t *incy, const doublereal *c, const doublereal *s);
void zdscal_(const aocl_int_t *n, const doublereal *da, dcomplex *zx, const aocl_int_t *incx);
void zgbmv_(const char *trans, const aocl_int_t *m, const aocl_int_t *n, const aocl_int_t *kl,
            const aocl_int_t *ku, const dcomplex *alpha, const dcomplex *a,
            const aocl_int_t *lda, const dcomplex *x, const aocl_int_t *incx,
            const dcomplex *beta, dcomplex *y, const aocl_int_t *incy);
void zgemm_(const char *transa, const char *transb, const aocl_int_t *m, const aocl_int_t *n,
            const aocl_int_t *k, const dcomplex *alpha, const dcomplex *a,
            const aocl_int_t *lda, const dcomplex *b, const aocl_int_t *ldb,
            const dcomplex *beta, dcomplex *c, const aocl_int_t *ldc);
void zgemmtr_(const char *uplo, const char *transa, const char *transb, const aocl_int_t *n,
              const aocl_int_t *k, const dcomplex *alpha, const dcomplex *a,
              const aocl_int_t *lda, const dcomplex *b, const aocl_int_t *ldb,
              const dcomplex *beta, dcomplex *c, const aocl_int_t *ldc);
void zgemv_(const char *trans, const aocl_int_t *m, const aocl_int_t *n, const dcomplex *alpha,
            const dcomplex *a, const aocl_int_t *lda, const dcomplex *x,
            const aocl_int_t *incx, const dcomplex *beta, dcomplex *y,
            const aocl_int_t *incy);
void zgerc_(const aocl_int_t *m, const aocl_int_t *n, const dcomplex *alpha,
            const dcomplex *x, const aocl_int_t *incx, const dcomplex *y,
            const aocl_int_t *incy, dcomplex *a, const aocl_int_t *lda);
void zgeru_(const aocl_int_t *m, const aocl_int_t *n, const dcomplex *alpha,
            const dcomplex *x, const aocl_int_t *incx, const dcomplex *y,
            const aocl_int_t *incy, dcomplex *a, const aocl_int_t *lda);
void zhbmv_(const char *uplo, const aocl_int_t *n, const aocl_int_t *k, const dcomplex *alpha,
            const dcomplex *a, const aocl_int_t *lda, const dcomplex *x,
            const aocl_int_t *incx, const dcomplex *beta, dcomplex *y,
            const aocl_int_t *incy);
void zhemm_(const char *side, const char *uplo, const aocl_int_t *m, const aocl_int_t *n,
            const dcomplex *alpha, const dcomplex *a, const aocl_int_t *lda,
            const dcomplex *b, const aocl_int_t *ldb, const dcomplex *beta,
            dcomplex *c, const aocl_int_t *ldc);
void zhemv_(const char *uplo, const aocl_int_t *n, const dcomplex *alpha,
            const dcomplex *a, const aocl_int_t *lda, const dcomplex *x,
            const aocl_int_t *incx, const dcomplex *beta, dcomplex *y,
            const aocl_int_t *incy);
void zher_(const char *uplo, const aocl_int_t *n, const doublereal *alpha, const dcomplex *x,
           const aocl_int_t *incx, dcomplex *a, const aocl_int_t *lda);
void zher2_(const char *uplo, const aocl_int_t *n, const dcomplex *alpha,
            const dcomplex *x, const aocl_int_t *incx, const dcomplex *y,
            const aocl_int_t *incy, dcomplex *a, const aocl_int_t *lda);
void zher2k_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
             const dcomplex *alpha, const dcomplex *a, const aocl_int_t *lda,
             const dcomplex *b, const aocl_int_t *ldb, const doublereal *beta,
             dcomplex *c, const aocl_int_t *ldc);
void zherk_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
            const doublereal *alpha, const dcomplex *a, const aocl_int_t *lda,
            const doublereal *beta, dcomplex *c, const aocl_int_t *ldc);
void zhpmv_(const char *uplo, const aocl_int_t *n, const dcomplex *alpha,
            const dcomplex *ap, const dcomplex *x, const aocl_int_t *incx,
            const dcomplex *beta, dcomplex *y, const aocl_int_t *incy);
void zhpr_(const char *uplo, const aocl_int_t *n, const doublereal *alpha, const dcomplex *x,
           const aocl_int_t *incx, dcomplex *ap);
void zhpr2_(const char *uplo, const aocl_int_t *n, const dcomplex *alpha,
            const dcomplex *x, const aocl_int_t *incx, const dcomplex *y,
            const aocl_int_t *incy, dcomplex *ap);
void zscal_(const aocl_int_t *n, const dcomplex *za, dcomplex *zx,
            const aocl_int_t *incx);
void zswap_(const aocl_int_t *n, dcomplex *zx, const aocl_int_t *incx, dcomplex *zy,
            const aocl_int_t *incy);
void zsymm_(const char *side, const char *uplo, const aocl_int_t *m, const aocl_int_t *n,
            const dcomplex *alpha, const dcomplex *a, const aocl_int_t *lda,
            const dcomplex *b, const aocl_int_t *ldb, const dcomplex *beta,
            dcomplex *c, const aocl_int_t *ldc);
void zsyr2k_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
             const dcomplex *alpha, const dcomplex *a, const aocl_int_t *lda,
             const dcomplex *b, const aocl_int_t *ldb, const dcomplex *beta,
             dcomplex *c, const aocl_int_t *ldc);
void zsyrk_(const char *uplo, const char *trans, const aocl_int_t *n, const aocl_int_t *k,
            const dcomplex *alpha, const dcomplex *a, const aocl_int_t *lda,
            const dcomplex *beta, dcomplex *c, const aocl_int_t *ldc);
void ztbmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const aocl_int_t *k, const dcomplex *a, const aocl_int_t *lda, dcomplex *x,
            const aocl_int_t *incx);
void ztbsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const aocl_int_t *k, const dcomplex *a, const aocl_int_t *lda, dcomplex *x,
            const aocl_int_t *incx);
void ztpmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const dcomplex *ap, dcomplex *x, const aocl_int_t *incx);
void ztpsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const dcomplex *ap, dcomplex *x, const aocl_int_t *incx);
void ztrmm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const aocl_int_t *m, const aocl_int_t *n, const dcomplex *alpha,
            const dcomplex *a, const aocl_int_t *lda, dcomplex *b, const aocl_int_t *ldb);
void ztrmv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const dcomplex *a, const aocl_int_t *lda, dcomplex *x,
            const aocl_int_t *incx);
void ztrsm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const aocl_int_t *m, const aocl_int_t *n, const dcomplex *alpha,
            const dcomplex *a, const aocl_int_t *lda, dcomplex *b, const aocl_int_t *ldb);
void ztrsv_(const char *uplo, const char *trans, const char *diag, const aocl_int_t *n,
            const dcomplex *a, const aocl_int_t *lda, dcomplex *x,
            const aocl_int_t *incx);
#ifdef FLA_ENABLE_BLAS_EXT_GEMMT
void sgemmt_(const char *uplo, const char *transa, const char *transb, const aocl_int_t *n,
             const aocl_int_t *k, const real *alpha, const real *a, const aocl_int_t *lda,
             const real *b, const aocl_int_t *ldb, const real *beta, real *c,
             const aocl_int_t *ldc);
void dgemmt_(const char *uplo, const char *transa, const char *transb, const aocl_int_t *n,
             const aocl_int_t *k, const doublereal *alpha, const doublereal *a,
             const aocl_int_t *lda, const doublereal *b, const aocl_int_t *ldb,
             const doublereal *beta, doublereal *c, const aocl_int_t *ldc);
void cgemmt_(const char *uplo, const char *transa, const char *transb, const aocl_int_t *n,
             const aocl_int_t *k, const scomplex *alpha, const scomplex *a,
             const aocl_int_t *lda, const scomplex *b, const aocl_int_t *ldb,
             const scomplex *beta, scomplex *c, const aocl_int_t *ldc);
void zgemmt_(const char *uplo, const char *transa, const char *transb, const aocl_int_t *n,
             const aocl_int_t *k, const dcomplex *alpha, const dcomplex *a,
             const aocl_int_t *lda, const dcomplex *b, const aocl_int_t *ldb,
             const dcomplex *beta, dcomplex *c, const aocl_int_t *ldc);
#endif /* FLA_ENABLE_BLAS_EXT_GEMMT */

#endif /* BLAS_INTERFACE_H */