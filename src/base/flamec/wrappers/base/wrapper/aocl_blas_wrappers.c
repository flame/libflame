/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
#include "FLA_f2c.h"
#include "blas_interface.h"

void aocl_blas_caxpy(const aocl_int64_t *n, const scomplex *ca, const scomplex *cx,
                     const aocl_int64_t *incx, scomplex *cy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    caxpy_(n, ca, cx, incx, cy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    caxpy_(&n_lp, ca, cx, &incx_lp, cy, &incy_lp);
#endif
}

void aocl_blas_ccopy(const aocl_int64_t *n, const scomplex *cx, const aocl_int64_t *incx,
                     scomplex *cy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    ccopy_(n, cx, incx, cy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    ccopy_(&n_lp, cx, &incx_lp, cy, &incy_lp);
#endif
}

#ifndef FLA_ENABLE_F2C_DOTC
void cdotc_f2c_(scomplex *r, const aocl_int_t *n, const scomplex *cx, const aocl_int_t *incx, const scomplex *cy,
                const aocl_int_t *incy)
{
    *r = cdotc_(n, cx, incx, cy, incy);
}
void cdotu_f2c_(scomplex *r, const aocl_int_t *n, const scomplex *cx, const aocl_int_t *incx, const scomplex *cy,
                const aocl_int_t *incy)
{
    *r = cdotu_(n, cx, incx, cy, incy);
}
#endif

void aocl_blas_cdotc(scomplex *r, const aocl_int64_t *n, const scomplex *cx,
                     const aocl_int64_t *incx, const scomplex *cy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    cdotc_f2c_(r, n, cx, incx, cy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    cdotc_f2c_(r, &n_lp, cx, &incx_lp, cy, &incy_lp);
#endif
}

void aocl_blas_cdotu(scomplex *r, const aocl_int64_t *n, const scomplex *cx,
                     const aocl_int64_t *incx,
                        const scomplex *cy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    cdotu_f2c_(r, n, cx, incx, cy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    cdotu_f2c_(r, &n_lp, cx, &incx_lp, cy, &incy_lp);
#endif
}

void aocl_blas_cgbmv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const aocl_int64_t *kl, const aocl_int64_t *ku, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, const scomplex *x,
                     const aocl_int64_t *incx, const scomplex *beta, scomplex *y,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    cgbmv_(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t kl_lp = (aocl_int_t)*kl;
    aocl_int_t ku_lp = (aocl_int_t)*ku;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    cgbmv_(trans, &m_lp, &n_lp, &kl_lp, &ku_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_cgemm(const char *transa, const char *transb, const aocl_int64_t *m,
                     const aocl_int64_t *n, const aocl_int64_t *k, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, const scomplex *b,
                     const aocl_int64_t *ldb, const scomplex *beta, scomplex *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    cgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    cgemm_(transa, transb, &m_lp, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_cgemmtr(const char *uplo, const char *transa, const char *transb,
                       const aocl_int64_t *n, const aocl_int64_t *k, const scomplex *alpha,
                       const scomplex *a, const aocl_int64_t *lda, const scomplex *b,
                       const aocl_int64_t *ldb, const scomplex *beta, scomplex *c,
                       const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    cgemmtr_(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    cgemmtr_(uplo, transa, transb, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_cgemv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const scomplex *alpha, const scomplex *a, const aocl_int64_t *lda,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *beta, scomplex *y,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    cgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    cgemv_(trans, &m_lp, &n_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_cgerc(const aocl_int64_t *m, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *y,
                     const aocl_int64_t *incy, scomplex *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    cgerc_(m, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    cgerc_(&m_lp, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

void aocl_blas_cgeru(const aocl_int64_t *m, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *y,
                     const aocl_int64_t *incy, scomplex *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    cgeru_(m, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    cgeru_(&m_lp, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

void aocl_blas_chbmv(const char *uplo, const aocl_int64_t *n, const aocl_int64_t *k,
                     const scomplex *alpha, const scomplex *a, const aocl_int64_t *lda,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *beta, scomplex *y,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    chbmv_(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    chbmv_(uplo, &n_lp, &k_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_chemm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const scomplex *alpha, const scomplex *a,
                     const aocl_int64_t *lda, const scomplex *b, const aocl_int64_t *ldb,
                     const scomplex *beta, scomplex *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    chemm_(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    chemm_(side, uplo, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_chemv(const char *uplo, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, const scomplex *x,
                     const aocl_int64_t *incx, const scomplex *beta, scomplex *y,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    chemv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    chemv_(uplo, &n_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_cher(const char *uplo, const aocl_int64_t *n, const real *alpha, const scomplex *x,
                    const aocl_int64_t *incx, scomplex *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    cher_(uplo, n, alpha, x, incx, a, lda);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    cher_(uplo, &n_lp, alpha, x, &incx_lp, a, &lda_lp);
#endif
}

void aocl_blas_cher2(const char *uplo, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *y,
                     const aocl_int64_t *incy, scomplex *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    cher2_(uplo, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    cher2_(uplo, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

void aocl_blas_cher2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const scomplex *alpha, const scomplex *a,
                      const aocl_int64_t *lda, const scomplex *b, const aocl_int64_t *ldb,
                      const real *beta, scomplex *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    cher2k_(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    cher2k_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_cherk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const real *alpha, const scomplex *a,
                     const aocl_int64_t *lda, const real *beta, scomplex *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    cherk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    cherk_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_chpmv(const char *uplo, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *ap, const scomplex *x, const aocl_int64_t *incx,
                     const scomplex *beta, scomplex *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    chpmv_(uplo, n, alpha, ap, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    chpmv_(uplo, &n_lp, alpha, ap, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_chpr(const char *uplo, const aocl_int64_t *n, const real *alpha, const scomplex *x,
                    const aocl_int64_t *incx, scomplex *ap)
{
#ifdef FLA_ENABLE_ILP64
    chpr_(uplo, n, alpha, x, incx, ap);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    chpr_(uplo, &n_lp, alpha, x, &incx_lp, ap);
#endif
}

void aocl_blas_chpr2(const char *uplo, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *y,
                     const aocl_int64_t *incy, scomplex *ap)
{
#ifdef FLA_ENABLE_ILP64
    chpr2_(uplo, n, alpha, x, incx, y, incy, ap);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    chpr2_(uplo, &n_lp, alpha, x, &incx_lp, y, &incy_lp, ap);
#endif
}

void aocl_blas_cscal(const aocl_int64_t *n, const scomplex *ca, scomplex *cx,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    cscal_(n, ca, cx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    cscal_(&n_lp, ca, cx, &incx_lp);
#endif
}

void aocl_blas_csrot(const aocl_int64_t *n, scomplex *cx, const aocl_int64_t *incx, scomplex *cy,
                     const aocl_int64_t *incy, const real *c, const real *s)
{
#ifdef FLA_ENABLE_ILP64
    csrot_(n, cx, incx, cy, incy, c, s);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    csrot_(&n_lp, cx, &incx_lp, cy, &incy_lp, c, s);
#endif
}

void aocl_blas_csscal(const aocl_int64_t *n, const real *sa, scomplex *cx, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    csscal_(n, sa, cx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    csscal_(&n_lp, sa, cx, &incx_lp);
#endif
}

void aocl_blas_cswap(const aocl_int64_t *n, scomplex *cx, const aocl_int64_t *incx, scomplex *cy,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    cswap_(n, cx, incx, cy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    cswap_(&n_lp, cx, &incx_lp, cy, &incy_lp);
#endif
}

void aocl_blas_csymm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const scomplex *alpha, const scomplex *a,
                     const aocl_int64_t *lda, const scomplex *b, const aocl_int64_t *ldb,
                     const scomplex *beta, scomplex *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    csymm_(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    csymm_(side, uplo, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_csyr2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const scomplex *alpha, const scomplex *a,
                      const aocl_int64_t *lda, const scomplex *b, const aocl_int64_t *ldb,
                      const scomplex *beta, scomplex *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    csyr2k_(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    csyr2k_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_csyrk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const scomplex *alpha, const scomplex *a,
                     const aocl_int64_t *lda, const scomplex *beta, scomplex *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    csyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    csyrk_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_ctbmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const scomplex *a, const aocl_int64_t *lda, scomplex *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ctbmv_(uplo, trans, diag, n, k, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ctbmv_(uplo, trans, diag, &n_lp, &k_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_ctbsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const scomplex *a, const aocl_int64_t *lda, scomplex *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ctbsv_(uplo, trans, diag, n, k, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ctbsv_(uplo, trans, diag, &n_lp, &k_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_ctpmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const scomplex *ap, scomplex *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ctpmv_(uplo, trans, diag, n, ap, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ctpmv_(uplo, trans, diag, &n_lp, ap, x, &incx_lp);
#endif
}

void aocl_blas_ctpsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const scomplex *ap, scomplex *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ctpsv_(uplo, trans, diag, n, ap, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ctpsv_(uplo, trans, diag, &n_lp, ap, x, &incx_lp);
#endif
}

void aocl_blas_ctrmm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, scomplex *b, const aocl_int64_t *ldb)
{
#ifdef FLA_ENABLE_ILP64
    ctrmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;

    ctrmm_(side, uplo, transa, diag, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp);
#endif
}

void aocl_blas_ctrmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const scomplex *a, const aocl_int64_t *lda, scomplex *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ctrmv_(uplo, trans, diag, n, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ctrmv_(uplo, trans, diag, &n_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_ctrsm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, scomplex *b, const aocl_int64_t *ldb)
{
#ifdef FLA_ENABLE_ILP64
    ctrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;

    ctrsm_(side, uplo, transa, diag, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp);
#endif
}

void aocl_blas_ctrsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const scomplex *a, const aocl_int64_t *lda, scomplex *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ctrsv_(uplo, trans, diag, n, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ctrsv_(uplo, trans, diag, &n_lp, a, &lda_lp, x, &incx_lp);
#endif
}

doublereal aocl_blas_dasum(const aocl_int64_t *n, const doublereal *dx, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return dasum_(n, dx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return dasum_(&n_lp, dx, &incx_lp);
#endif
}

void aocl_blas_daxpy(const aocl_int64_t *n, const doublereal *da, const doublereal *dx,
                     const aocl_int64_t *incx, doublereal *dy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    daxpy_(n, da, dx, incx, dy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    daxpy_(&n_lp, da, dx, &incx_lp, dy, &incy_lp);
#endif
}

void aocl_blas_dcopy(const aocl_int64_t *n, const doublereal *dx, const aocl_int64_t *incx,
                     doublereal *dy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    dcopy_(n, dx, incx, dy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    dcopy_(&n_lp, dx, &incx_lp, dy, &incy_lp);
#endif
}

doublereal aocl_blas_ddot(const aocl_int64_t *n, const doublereal *dx, const aocl_int64_t *incx,
                          const doublereal *dy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    return ddot_(n, dx, incx, dy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    return ddot_(&n_lp, dx, &incx_lp, dy, &incy_lp);
#endif
}

void aocl_blas_dgbmv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const aocl_int64_t *kl, const aocl_int64_t *ku, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, const doublereal *x,
                     const aocl_int64_t *incx, const doublereal *beta, doublereal *y,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    dgbmv_(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t kl_lp = (aocl_int_t)*kl;
    aocl_int_t ku_lp = (aocl_int_t)*ku;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    dgbmv_(trans, &m_lp, &n_lp, &kl_lp, &ku_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_dgemm(const char *transa, const char *transb, const aocl_int64_t *m,
                     const aocl_int64_t *n, const aocl_int64_t *k, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, const doublereal *b,
                     const aocl_int64_t *ldb, const doublereal *beta, doublereal *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    dgemm_(transa, transb, &m_lp, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_dgemmtr(const char *uplo, const char *transa, const char *transb,
                       const aocl_int64_t *n, const aocl_int64_t *k, const doublereal *alpha,
                       const doublereal *a, const aocl_int64_t *lda, const doublereal *b,
                       const aocl_int64_t *ldb, const doublereal *beta, doublereal *c,
                       const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    dgemmtr_(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    dgemmtr_(uplo, transa, transb, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_dgemv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const doublereal *alpha, const doublereal *a, const aocl_int64_t *lda,
                     const doublereal *x, const aocl_int64_t *incx, const doublereal *beta,
                     doublereal *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    dgemv_(trans, &m_lp, &n_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_dger(const aocl_int64_t *m, const aocl_int64_t *n, const doublereal *alpha,
                    const doublereal *x, const aocl_int64_t *incx, const doublereal *y,
                    const aocl_int64_t *incy, doublereal *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    dger_(m, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    dger_(&m_lp, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

doublereal aocl_blas_dnrm2(const aocl_int64_t *n, const doublereal *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return dnrm2_(n, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return dnrm2_(&n_lp, x, &incx_lp);
#endif
}

void aocl_blas_drot(const aocl_int64_t *n, doublereal *dx, const aocl_int64_t *incx, doublereal *dy,
                    const aocl_int64_t *incy, const doublereal *c, const doublereal *s)
{
#ifdef FLA_ENABLE_ILP64
    drot_(n, dx, incx, dy, incy, c, s);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    drot_(&n_lp, dx, &incx_lp, dy, &incy_lp, c, s);
#endif
}

void aocl_blas_drotm(const aocl_int64_t *n, doublereal *dx, const aocl_int64_t *incx,
                     doublereal *dy, const aocl_int64_t *incy, const doublereal *dparam)
{
#ifdef FLA_ENABLE_ILP64
    drotm_(n, dx, incx, dy, incy, dparam);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    drotm_(&n_lp, dx, &incx_lp, dy, &incy_lp, dparam);
#endif
}

void aocl_blas_dsbmv(const char *uplo, const aocl_int64_t *n, const aocl_int64_t *k,
                     const doublereal *alpha, const doublereal *a, const aocl_int64_t *lda,
                     const doublereal *x, const aocl_int64_t *incx, const doublereal *beta,
                     doublereal *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    dsbmv_(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    dsbmv_(uplo, &n_lp, &k_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_dscal(const aocl_int64_t *n, const doublereal *da, doublereal *dx,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    dscal_(n, da, dx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    dscal_(&n_lp, da, dx, &incx_lp);
#endif
}

doublereal aocl_blas_dsdot(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx,
                           const real *sy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    return dsdot_(n, sx, incx, sy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    return dsdot_(&n_lp, sx, &incx_lp, sy, &incy_lp);
#endif
}

void aocl_blas_dspmv(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *ap, const doublereal *x, const aocl_int64_t *incx,
                     const doublereal *beta, doublereal *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    dspmv_(uplo, n, alpha, ap, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    dspmv_(uplo, &n_lp, alpha, ap, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_dspr(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                    const doublereal *x, const aocl_int64_t *incx, doublereal *ap)
{
#ifdef FLA_ENABLE_ILP64
    dspr_(uplo, n, alpha, x, incx, ap);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    dspr_(uplo, &n_lp, alpha, x, &incx_lp, ap);
#endif
}

void aocl_blas_dspr2(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *x, const aocl_int64_t *incx, const doublereal *y,
                     const aocl_int64_t *incy, doublereal *ap)
{
#ifdef FLA_ENABLE_ILP64
    dspr2_(uplo, n, alpha, x, incx, y, incy, ap);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    dspr2_(uplo, &n_lp, alpha, x, &incx_lp, y, &incy_lp, ap);
#endif
}

void aocl_blas_dswap(const aocl_int64_t *n, doublereal *dx, const aocl_int64_t *incx,
                     doublereal *dy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    dswap_(n, dx, incx, dy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    dswap_(&n_lp, dx, &incx_lp, dy, &incy_lp);
#endif
}

void aocl_blas_dsymm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const doublereal *alpha, const doublereal *a,
                     const aocl_int64_t *lda, const doublereal *b, const aocl_int64_t *ldb,
                     const doublereal *beta, doublereal *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    dsymm_(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    dsymm_(side, uplo, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_dsymv(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, const doublereal *x,
                     const aocl_int64_t *incx, const doublereal *beta, doublereal *y,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    dsymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    dsymv_(uplo, &n_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_dsyr(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                    const doublereal *x, const aocl_int64_t *incx, doublereal *a,
                    const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    dsyr_(uplo, n, alpha, x, incx, a, lda);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    dsyr_(uplo, &n_lp, alpha, x, &incx_lp, a, &lda_lp);
#endif
}

void aocl_blas_dsyr2(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *x, const aocl_int64_t *incx, const doublereal *y,
                     const aocl_int64_t *incy, doublereal *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    dsyr2_(uplo, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    dsyr2_(uplo, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

void aocl_blas_dsyr2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const doublereal *alpha, const doublereal *a,
                      const aocl_int64_t *lda, const doublereal *b, const aocl_int64_t *ldb,
                      const doublereal *beta, doublereal *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    dsyr2k_(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    dsyr2k_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_dsyrk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const doublereal *alpha, const doublereal *a,
                     const aocl_int64_t *lda, const doublereal *beta, doublereal *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    dsyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    dsyrk_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_dtbmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const doublereal *a, const aocl_int64_t *lda,
                     doublereal *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    dtbmv_(uplo, trans, diag, n, k, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    dtbmv_(uplo, trans, diag, &n_lp, &k_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_dtbsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const doublereal *a, const aocl_int64_t *lda,
                     doublereal *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    dtbsv_(uplo, trans, diag, n, k, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    dtbsv_(uplo, trans, diag, &n_lp, &k_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_dtpmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const doublereal *ap, doublereal *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    dtpmv_(uplo, trans, diag, n, ap, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    dtpmv_(uplo, trans, diag, &n_lp, ap, x, &incx_lp);
#endif
}

void aocl_blas_dtpsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const doublereal *ap, doublereal *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    dtpsv_(uplo, trans, diag, n, ap, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    dtpsv_(uplo, trans, diag, &n_lp, ap, x, &incx_lp);
#endif
}

void aocl_blas_dtrmm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, doublereal *b,
                     const aocl_int64_t *ldb)
{
#ifdef FLA_ENABLE_ILP64
    dtrmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;

    dtrmm_(side, uplo, transa, diag, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp);
#endif
}

void aocl_blas_dtrmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const doublereal *a, const aocl_int64_t *lda, doublereal *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    dtrmv_(uplo, trans, diag, n, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    dtrmv_(uplo, trans, diag, &n_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_dtrsm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, doublereal *b,
                     const aocl_int64_t *ldb)
{
#ifdef FLA_ENABLE_ILP64
    dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;

    dtrsm_(side, uplo, transa, diag, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp);
#endif
}

void aocl_blas_dtrsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const doublereal *a, const aocl_int64_t *lda, doublereal *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    dtrsv_(uplo, trans, diag, n, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    dtrsv_(uplo, trans, diag, &n_lp, a, &lda_lp, x, &incx_lp);
#endif
}

doublereal aocl_blas_dzasum(const aocl_int64_t *n, dcomplex *zx, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return dzasum_(n, zx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return dzasum_(&n_lp, zx, &incx_lp);
#endif
}

doublereal aocl_blas_dznrm2(const aocl_int64_t *n, const dcomplex *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return dznrm2_(n, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return dznrm2_(&n_lp, x, &incx_lp);
#endif
}

aocl_int64_t aocl_blas_icamax(const aocl_int64_t *n, const scomplex *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return icamax_(n, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return icamax_(&n_lp, x, &incx_lp);
#endif
}

aocl_int64_t aocl_blas_idamax(const aocl_int64_t *n, const doublereal *dx, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return idamax_(n, dx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return idamax_(&n_lp, dx, &incx_lp);
#endif
}

aocl_int64_t aocl_blas_isamax(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return isamax_(n, sx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return isamax_(&n_lp, sx, &incx_lp);
#endif
}

aocl_int64_t aocl_blas_izamax(const aocl_int64_t *n, const dcomplex *x,
                              const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return izamax_(n, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return izamax_(&n_lp, x, &incx_lp);
#endif
}

real aocl_blas_sasum(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return sasum_(n, sx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return sasum_(&n_lp, sx, &incx_lp);
#endif
}

void aocl_blas_saxpy(const aocl_int64_t *n, const real *sa, const real *sx,
                     const aocl_int64_t *incx, real *sy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    saxpy_(n, sa, sx, incx, sy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    saxpy_(&n_lp, sa, sx, &incx_lp, sy, &incy_lp);
#endif
}

real aocl_blas_scasum(const aocl_int64_t *n, scomplex *cx, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return scasum_(n, cx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return scasum_(&n_lp, cx, &incx_lp);
#endif
}

real aocl_blas_scnrm2(const aocl_int64_t *n, const scomplex *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return scnrm2_(n, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return scnrm2_(&n_lp, x, &incx_lp);
#endif
}

void aocl_blas_scopy(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx, real *sy,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    scopy_(n, sx, incx, sy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    scopy_(&n_lp, sx, &incx_lp, sy, &incy_lp);
#endif
}

real aocl_blas_sdot(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx, const real *sy,
                    const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    return sdot_(n, sx, incx, sy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    return sdot_(&n_lp, sx, &incx_lp, sy, &incy_lp);
#endif
}

real aocl_blas_sdsdot(const aocl_int64_t *n, const real *sb, const real *sx,
                      const aocl_int64_t *incx, const real *sy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    return sdsdot_(n, sb, sx, incx, sy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    return sdsdot_(&n_lp, sb, sx, &incx_lp, sy, &incy_lp);
#endif
}

void aocl_blas_sgbmv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const aocl_int64_t *kl, const aocl_int64_t *ku, const real *alpha,
                     const real *a, const aocl_int64_t *lda, const real *x,
                     const aocl_int64_t *incx, const real *beta, real *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    sgbmv_(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t kl_lp = (aocl_int_t)*kl;
    aocl_int_t ku_lp = (aocl_int_t)*ku;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    sgbmv_(trans, &m_lp, &n_lp, &kl_lp, &ku_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_sgemm(const char *transa, const char *transb, const aocl_int64_t *m,
                     const aocl_int64_t *n, const aocl_int64_t *k, const real *alpha, const real *a,
                     const aocl_int64_t *lda, const real *b, const aocl_int64_t *ldb,
                     const real *beta, real *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    sgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    sgemm_(transa, transb, &m_lp, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_sgemmtr(const char *uplo, const char *transa, const char *transb,
                       const aocl_int64_t *n, const aocl_int64_t *k, const real *alpha,
                       const real *a, const aocl_int64_t *lda, const real *b,
                       const aocl_int64_t *ldb, const real *beta, real *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    sgemmtr_(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    sgemmtr_(uplo, transa, transb, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_sgemv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const real *alpha, const real *a, const aocl_int64_t *lda, const real *x,
                     const aocl_int64_t *incx, const real *beta, real *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    sgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    sgemv_(trans, &m_lp, &n_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_sger(const aocl_int64_t *m, const aocl_int64_t *n, const real *alpha, const real *x,
                    const aocl_int64_t *incx, const real *y, const aocl_int64_t *incy, real *a,
                    const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    sger_(m, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    sger_(&m_lp, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

real aocl_blas_snrm2(const aocl_int64_t *n, const real *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    return snrm2_(n, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    return snrm2_(&n_lp, x, &incx_lp);
#endif
}

void aocl_blas_srot(const aocl_int64_t *n, real *sx, const aocl_int64_t *incx, real *sy,
                    const aocl_int64_t *incy, const real *c, const real *s)
{
#ifdef FLA_ENABLE_ILP64
    srot_(n, sx, incx, sy, incy, c, s);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    srot_(&n_lp, sx, &incx_lp, sy, &incy_lp, c, s);
#endif
}

void aocl_blas_srotm(const aocl_int64_t *n, real *sx, const aocl_int64_t *incx, real *sy,
                     const aocl_int64_t *incy, const real *sparam)
{
#ifdef FLA_ENABLE_ILP64
    srotm_(n, sx, incx, sy, incy, sparam);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    srotm_(&n_lp, sx, &incx_lp, sy, &incy_lp, sparam);
#endif
}

void aocl_blas_ssbmv(const char *uplo, const aocl_int64_t *n, const aocl_int64_t *k,
                     const real *alpha, const real *a, const aocl_int64_t *lda, const real *x,
                     const aocl_int64_t *incx, const real *beta, real *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    ssbmv_(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    ssbmv_(uplo, &n_lp, &k_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_sscal(const aocl_int64_t *n, const real *sa, real *sx, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    sscal_(n, sa, sx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    sscal_(&n_lp, sa, sx, &incx_lp);
#endif
}

void aocl_blas_sspmv(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *ap,
                     const real *x, const aocl_int64_t *incx, const real *beta, real *y,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    sspmv_(uplo, n, alpha, ap, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    sspmv_(uplo, &n_lp, alpha, ap, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_sspr(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *x,
                    const aocl_int64_t *incx, real *ap)
{
#ifdef FLA_ENABLE_ILP64
    sspr_(uplo, n, alpha, x, incx, ap);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    sspr_(uplo, &n_lp, alpha, x, &incx_lp, ap);
#endif
}

void aocl_blas_sspr2(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *x,
                     const aocl_int64_t *incx, const real *y, const aocl_int64_t *incy, real *ap)
{
#ifdef FLA_ENABLE_ILP64
    sspr2_(uplo, n, alpha, x, incx, y, incy, ap);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    sspr2_(uplo, &n_lp, alpha, x, &incx_lp, y, &incy_lp, ap);
#endif
}

void aocl_blas_sswap(const aocl_int64_t *n, real *sx, const aocl_int64_t *incx, real *sy,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    sswap_(n, sx, incx, sy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    sswap_(&n_lp, sx, &incx_lp, sy, &incy_lp);
#endif
}

void aocl_blas_ssymm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const real *alpha, const real *a,
                     const aocl_int64_t *lda, const real *b, const aocl_int64_t *ldb,
                     const real *beta, real *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    ssymm_(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    ssymm_(side, uplo, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_ssymv(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *a,
                     const aocl_int64_t *lda, const real *x, const aocl_int64_t *incx,
                     const real *beta, real *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    ssymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    ssymv_(uplo, &n_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_ssyr(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *x,
                    const aocl_int64_t *incx, real *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    ssyr_(uplo, n, alpha, x, incx, a, lda);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    ssyr_(uplo, &n_lp, alpha, x, &incx_lp, a, &lda_lp);
#endif
}

void aocl_blas_ssyr2(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *x,
                     const aocl_int64_t *incx, const real *y, const aocl_int64_t *incy, real *a,
                     const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    ssyr2_(uplo, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    ssyr2_(uplo, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

void aocl_blas_ssyr2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const real *alpha, const real *a,
                      const aocl_int64_t *lda, const real *b, const aocl_int64_t *ldb,
                      const real *beta, real *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    ssyr2k_(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    ssyr2k_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_ssyrk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const real *alpha, const real *a,
                     const aocl_int64_t *lda, const real *beta, real *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    ssyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    ssyrk_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_stbmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const real *a, const aocl_int64_t *lda, real *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    stbmv_(uplo, trans, diag, n, k, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    stbmv_(uplo, trans, diag, &n_lp, &k_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_stbsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const real *a, const aocl_int64_t *lda, real *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    stbsv_(uplo, trans, diag, n, k, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    stbsv_(uplo, trans, diag, &n_lp, &k_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_stpmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const real *ap, real *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    stpmv_(uplo, trans, diag, n, ap, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    stpmv_(uplo, trans, diag, &n_lp, ap, x, &incx_lp);
#endif
}

void aocl_blas_stpsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const real *ap, real *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    stpsv_(uplo, trans, diag, n, ap, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    stpsv_(uplo, trans, diag, &n_lp, ap, x, &incx_lp);
#endif
}

void aocl_blas_strmm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const real *alpha, const real *a,
                     const aocl_int64_t *lda, real *b, const aocl_int64_t *ldb)
{
#ifdef FLA_ENABLE_ILP64
    strmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;

    strmm_(side, uplo, transa, diag, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp);
#endif
}

void aocl_blas_strmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const real *a, const aocl_int64_t *lda, real *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    strmv_(uplo, trans, diag, n, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    strmv_(uplo, trans, diag, &n_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_strsm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const real *alpha, const real *a,
                     const aocl_int64_t *lda, real *b, const aocl_int64_t *ldb)
{
#ifdef FLA_ENABLE_ILP64
    strsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;

    strsm_(side, uplo, transa, diag, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp);
#endif
}

void aocl_blas_strsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const real *a, const aocl_int64_t *lda, real *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    strsv_(uplo, trans, diag, n, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    strsv_(uplo, trans, diag, &n_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_xerbla(const char *srname, const aocl_int64_t *info, aocl_int64_t srname_len)
{
#ifdef FLA_ENABLE_ILP64
    xerbla_(srname, info, srname_len);
#else
    aocl_int_t info_lp = (aocl_int_t)*info;

    xerbla_(srname, &info_lp, (aocl_int_t)srname_len);
#endif
}

void aocl_blas_xerbla_array(const char *srname_array, const aocl_int64_t *srname_len,
                            const aocl_int64_t *info)
{
#ifdef FLA_ENABLE_ILP64
    xerbla_array_(srname_array, srname_len, info);
#else
    aocl_int_t srname_len_lp = (aocl_int_t)*srname_len;
    aocl_int_t info_lp = (aocl_int_t)*info;

    xerbla_array_(srname_array, &srname_len_lp, &info_lp);
#endif
}

void aocl_blas_zaxpy(const aocl_int64_t *n, const dcomplex *za, const dcomplex *zx,
                     const aocl_int64_t *incx, dcomplex *zy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zaxpy_(n, za, zx, incx, zy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zaxpy_(&n_lp, za, zx, &incx_lp, zy, &incy_lp);
#endif
}

void aocl_blas_zcopy(const aocl_int64_t *n, const dcomplex *zx, const aocl_int64_t *incx,
                     dcomplex *zy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zcopy_(n, zx, incx, zy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zcopy_(&n_lp, zx, &incx_lp, zy, &incy_lp);
#endif
}

#ifndef FLA_ENABLE_F2C_DOTC
void zdotc_f2c_(dcomplex *r, const aocl_int_t *n, const dcomplex *cx, const aocl_int_t *incx, const dcomplex *cy,
                const aocl_int_t *incy)
{
    *r = zdotc_(n, cx, incx, cy, incy);
}
void zdotu_f2c_(dcomplex *r, const aocl_int_t *n, const dcomplex *cx, const aocl_int_t *incx, const dcomplex *cy,
                const aocl_int_t *incy)
{
    *r = zdotu_(n, cx, incx, cy, incy);
}
#endif

void aocl_blas_zdotc(dcomplex *r, const aocl_int64_t *n, const dcomplex *zx,
                              const aocl_int64_t *incx, const dcomplex *zy,
                              const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zdotc_f2c_(r, n, zx, incx, zy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zdotc_f2c_(r, &n_lp, zx, &incx_lp, zy, &incy_lp);
#endif
}

void aocl_blas_zdotu(dcomplex *r, const aocl_int64_t *n, const dcomplex *zx,
                              const aocl_int64_t *incx, const dcomplex *zy,
                              const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zdotu_f2c_(r, n, zx, incx, zy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zdotu_f2c_(r, &n_lp, zx, &incx_lp, zy, &incy_lp);
#endif
}

void aocl_blas_zdrot(const aocl_int64_t *n, dcomplex *zx, const aocl_int64_t *incx,
                     dcomplex *zy, const aocl_int64_t *incy, const doublereal *c,
                     const doublereal *s)
{
#ifdef FLA_ENABLE_ILP64
    zdrot_(n, zx, incx, zy, incy, c, s);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zdrot_(&n_lp, zx, &incx_lp, zy, &incy_lp, c, s);
#endif
}

void aocl_blas_zdscal(const aocl_int64_t *n, const doublereal *da, dcomplex *zx,
                      const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    zdscal_(n, da, zx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    zdscal_(&n_lp, da, zx, &incx_lp);
#endif
}

void aocl_blas_zgbmv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const aocl_int64_t *kl, const aocl_int64_t *ku, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, const dcomplex *x,
                     const aocl_int64_t *incx, const dcomplex *beta, dcomplex *y,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zgbmv_(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t kl_lp = (aocl_int_t)*kl;
    aocl_int_t ku_lp = (aocl_int_t)*ku;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zgbmv_(trans, &m_lp, &n_lp, &kl_lp, &ku_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_zgemm(const char *transa, const char *transb, const aocl_int64_t *m,
                     const aocl_int64_t *n, const aocl_int64_t *k, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, const dcomplex *b,
                     const aocl_int64_t *ldb, const dcomplex *beta, dcomplex *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    zgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    zgemm_(transa, transb, &m_lp, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_zgemmtr(const char *uplo, const char *transa, const char *transb,
                       const aocl_int64_t *n, const aocl_int64_t *k, const dcomplex *alpha,
                       const dcomplex *a, const aocl_int64_t *lda, const dcomplex *b,
                       const aocl_int64_t *ldb, const dcomplex *beta, dcomplex *c,
                       const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    zgemmtr_(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    zgemmtr_(uplo, transa, transb, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_zgemv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const dcomplex *alpha, const dcomplex *a, const aocl_int64_t *lda,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *beta,
                     dcomplex *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zgemv_(trans, &m_lp, &n_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_zgerc(const aocl_int64_t *m, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *y,
                     const aocl_int64_t *incy, dcomplex *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    zgerc_(m, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    zgerc_(&m_lp, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

void aocl_blas_zgeru(const aocl_int64_t *m, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *y,
                     const aocl_int64_t *incy, dcomplex *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    zgeru_(m, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    zgeru_(&m_lp, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

void aocl_blas_zhbmv(const char *uplo, const aocl_int64_t *n, const aocl_int64_t *k,
                     const dcomplex *alpha, const dcomplex *a, const aocl_int64_t *lda,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *beta,
                     dcomplex *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zhbmv_(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zhbmv_(uplo, &n_lp, &k_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_zhemm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const dcomplex *alpha, const dcomplex *a,
                     const aocl_int64_t *lda, const dcomplex *b, const aocl_int64_t *ldb,
                     const dcomplex *beta, dcomplex *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    zhemm_(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    zhemm_(side, uplo, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_zhemv(const char *uplo, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, const dcomplex *x,
                     const aocl_int64_t *incx, const dcomplex *beta, dcomplex *y,
                     const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zhemv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zhemv_(uplo, &n_lp, alpha, a, &lda_lp, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_zher(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                    const dcomplex *x, const aocl_int64_t *incx, dcomplex *a,
                    const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    zher_(uplo, n, alpha, x, incx, a, lda);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    zher_(uplo, &n_lp, alpha, x, &incx_lp, a, &lda_lp);
#endif
}

void aocl_blas_zher2(const char *uplo, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *y,
                     const aocl_int64_t *incy, dcomplex *a, const aocl_int64_t *lda)
{
#ifdef FLA_ENABLE_ILP64
    zher2_(uplo, n, alpha, x, incx, y, incy, a, lda);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;
    aocl_int_t lda_lp = (aocl_int_t)*lda;

    zher2_(uplo, &n_lp, alpha, x, &incx_lp, y, &incy_lp, a, &lda_lp);
#endif
}

void aocl_blas_zher2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const dcomplex *alpha, const dcomplex *a,
                      const aocl_int64_t *lda, const dcomplex *b, const aocl_int64_t *ldb,
                      const doublereal *beta, dcomplex *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    zher2k_(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    zher2k_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_zherk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const doublereal *alpha, const dcomplex *a,
                     const aocl_int64_t *lda, const doublereal *beta, dcomplex *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    zherk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    zherk_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_zhpmv(const char *uplo, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *ap, const dcomplex *x, const aocl_int64_t *incx,
                     const dcomplex *beta, dcomplex *y, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zhpmv_(uplo, n, alpha, ap, x, incx, beta, y, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zhpmv_(uplo, &n_lp, alpha, ap, x, &incx_lp, beta, y, &incy_lp);
#endif
}

void aocl_blas_zhpr(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                    const dcomplex *x, const aocl_int64_t *incx, dcomplex *ap)
{
#ifdef FLA_ENABLE_ILP64
    zhpr_(uplo, n, alpha, x, incx, ap);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    zhpr_(uplo, &n_lp, alpha, x, &incx_lp, ap);
#endif
}

void aocl_blas_zhpr2(const char *uplo, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *y,
                     const aocl_int64_t *incy, dcomplex *ap)
{
#ifdef FLA_ENABLE_ILP64
    zhpr2_(uplo, n, alpha, x, incx, y, incy, ap);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zhpr2_(uplo, &n_lp, alpha, x, &incx_lp, y, &incy_lp, ap);
#endif
}

void aocl_blas_zscal(const aocl_int64_t *n, const dcomplex *za, dcomplex *zx,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    zscal_(n, za, zx, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    zscal_(&n_lp, za, zx, &incx_lp);
#endif
}

void aocl_blas_zswap(const aocl_int64_t *n, dcomplex *zx, const aocl_int64_t *incx,
                     dcomplex *zy, const aocl_int64_t *incy)
{
#ifdef FLA_ENABLE_ILP64
    zswap_(n, zx, incx, zy, incy);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;
    aocl_int_t incy_lp = (aocl_int_t)*incy;

    zswap_(&n_lp, zx, &incx_lp, zy, &incy_lp);
#endif
}

void aocl_blas_zsymm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const dcomplex *alpha, const dcomplex *a,
                     const aocl_int64_t *lda, const dcomplex *b, const aocl_int64_t *ldb,
                     const dcomplex *beta, dcomplex *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    zsymm_(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    zsymm_(side, uplo, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_zsyr2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const dcomplex *alpha, const dcomplex *a,
                      const aocl_int64_t *lda, const dcomplex *b, const aocl_int64_t *ldb,
                      const dcomplex *beta, dcomplex *c, const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    zsyr2k_(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    zsyr2k_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_zsyrk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const dcomplex *alpha, const dcomplex *a,
                     const aocl_int64_t *lda, const dcomplex *beta, dcomplex *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    zsyrk_(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    zsyrk_(uplo, trans, &n_lp, &k_lp, alpha, a, &lda_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_ztbmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const dcomplex *a, const aocl_int64_t *lda,
                     dcomplex *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ztbmv_(uplo, trans, diag, n, k, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ztbmv_(uplo, trans, diag, &n_lp, &k_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_ztbsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const dcomplex *a, const aocl_int64_t *lda,
                     dcomplex *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ztbsv_(uplo, trans, diag, n, k, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ztbsv_(uplo, trans, diag, &n_lp, &k_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_ztpmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const dcomplex *ap, dcomplex *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ztpmv_(uplo, trans, diag, n, ap, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ztpmv_(uplo, trans, diag, &n_lp, ap, x, &incx_lp);
#endif
}

void aocl_blas_ztpsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const dcomplex *ap, dcomplex *x, const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ztpsv_(uplo, trans, diag, n, ap, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ztpsv_(uplo, trans, diag, &n_lp, ap, x, &incx_lp);
#endif
}

void aocl_blas_ztrmm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, dcomplex *b,
                     const aocl_int64_t *ldb)
{
#ifdef FLA_ENABLE_ILP64
    ztrmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;

    ztrmm_(side, uplo, transa, diag, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp);
#endif
}

void aocl_blas_ztrmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const dcomplex *a, const aocl_int64_t *lda, dcomplex *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ztrmv_(uplo, trans, diag, n, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ztrmv_(uplo, trans, diag, &n_lp, a, &lda_lp, x, &incx_lp);
#endif
}

void aocl_blas_ztrsm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, dcomplex *b,
                     const aocl_int64_t *ldb)
{
#ifdef FLA_ENABLE_ILP64
    ztrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
    aocl_int_t m_lp = (aocl_int_t)*m;
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;

    ztrsm_(side, uplo, transa, diag, &m_lp, &n_lp, alpha, a, &lda_lp, b, &ldb_lp);
#endif
}

void aocl_blas_ztrsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const dcomplex *a, const aocl_int64_t *lda, dcomplex *x,
                     const aocl_int64_t *incx)
{
#ifdef FLA_ENABLE_ILP64
    ztrsv_(uplo, trans, diag, n, a, lda, x, incx);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t incx_lp = (aocl_int_t)*incx;

    ztrsv_(uplo, trans, diag, &n_lp, a, &lda_lp, x, &incx_lp);
#endif
}

#ifdef FLA_ENABLE_BLAS_EXT_GEMMT

void aocl_blas_cgemmt(const char* uplo, const char *transa, const char *transb,
                     const aocl_int64_t *n, const aocl_int64_t *k, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, const scomplex *b,
                     const aocl_int64_t *ldb, const scomplex *beta, scomplex *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    cgemmt_(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    cgemmt_(uplo, transa, transb, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_dgemmt(const char* uplo, const char *transa, const char *transb,
                     const aocl_int64_t *n, const aocl_int64_t *k, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, const doublereal *b,
                     const aocl_int64_t *ldb, const doublereal *beta, doublereal *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    dgemmt_(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    dgemmt_(uplo, transa, transb, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_sgemmt(const char* uplo, const char *transa, const char *transb,
                     const aocl_int64_t *n, const aocl_int64_t *k, const real *alpha,
                     const real *a, const aocl_int64_t *lda, const real *b,
                     const aocl_int64_t *ldb, const real *beta, real *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    sgemmt_(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    sgemmt_(uplo, transa, transb, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}

void aocl_blas_zgemmt(const char* uplo, const char *transa, const char *transb,
                     const aocl_int64_t *n, const aocl_int64_t *k, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, const dcomplex *b,
                     const aocl_int64_t *ldb, const dcomplex *beta, dcomplex *c,
                     const aocl_int64_t *ldc)
{
#ifdef FLA_ENABLE_ILP64
    zgemmt_(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
    aocl_int_t n_lp = (aocl_int_t)*n;
    aocl_int_t k_lp = (aocl_int_t)*k;
    aocl_int_t lda_lp = (aocl_int_t)*lda;
    aocl_int_t ldb_lp = (aocl_int_t)*ldb;
    aocl_int_t ldc_lp = (aocl_int_t)*ldc;

    zgemmt_(uplo, transa, transb, &n_lp, &k_lp, alpha, a, &lda_lp, b, &ldb_lp, beta, c, &ldc_lp);
#endif
}


#endif /* FLA_ENABLE_BLAS_EXT_GEMMT */
