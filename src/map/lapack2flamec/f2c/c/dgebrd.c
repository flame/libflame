/*
    Copyright (C) 2022-2026, Advanced Micro Devices, Inc. All rights reserved.
*/
/* dgebrd.f -- translated by f2c (version 20160102). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c_n1 = -1;
static aocl_int64_t c__3 = 3;
static aocl_int64_t c__2 = 2;
static doublereal c_b21 = -1.;
static doublereal c_b22 = 1.;

extern int fla_thread_get_num_threads(void);

/* Select the number of threads to use for the DLABRD routine based on the size of the matrix*/
static int fla_get_dlabrd_threads(aocl_int64_t m, aocl_int64_t n, int max_threads)
{
    aocl_int64_t min_m_n = fla_min(m, n);
    int chosen_threads;
    if(min_m_n < 1000)
    {
        chosen_threads = 8;
    }
    else if(min_m_n < 4600)
    {
        chosen_threads = 64;
    }
    else if(min_m_n < 6000)
    {
        chosen_threads = 96;
    }
    else
    {
        chosen_threads = 192;
    }

    if(chosen_threads > max_threads)
    {
        chosen_threads = max_threads;
    }
    if(chosen_threads < 1)
    {
        chosen_threads = 1;
    }

    return chosen_threads;
}
/* > \brief \b DGEBRD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGEBRD + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebrd.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebrd.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebrd.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), D( * ), E( * ), TAUP( * ), */
/* $ TAUQ( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEBRD reduces a general real M-by-N matrix A to upper or lower */
/* > bidiagonal form B by an orthogonal transformation: Q**T * A * P = B. */
/* > */
/* > If m >= n, B is upper bidiagonal;
if m < n, B is lower bidiagonal. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows in the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns in the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the M-by-N general matrix to be reduced. */
/* > On exit, */
/* > if m >= n, the diagonal and the first superdiagonal are */
/* > overwritten with the upper bidiagonal matrix B;
the */
/* > elements below the diagonal, with the array TAUQ, represent */
/* > the orthogonal matrix Q as a product of elementary */
/* > reflectors, and the elements above the first superdiagonal, */
/* > with the array TAUP, represent the orthogonal matrix P as */
/* > a product of elementary reflectors;
 */
/* > if m < n, the diagonal and the first subdiagonal are */
/* > overwritten with the lower bidiagonal matrix B;
the */
/* > elements below the first subdiagonal, with the array TAUQ, */
/* > represent the orthogonal matrix Q as a product of */
/* > elementary reflectors, and the elements above the diagonal, */
/* > with the array TAUP, represent the orthogonal matrix P as */
/* > a product of elementary reflectors. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (fla_min(M,N)) */
/* > The diagonal elements of the bidiagonal matrix B: */
/* > D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (fla_min(M,N)-1) */
/* > The off-diagonal elements of the bidiagonal matrix B: */
/* > if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
 */
/* > if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ */
/* > \verbatim */
/* > TAUQ is DOUBLE PRECISION array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors which */
/* > represent the orthogonal matrix Q. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP */
/* > \verbatim */
/* > TAUP is DOUBLE PRECISION array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors which */
/* > represent the orthogonal matrix P. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of the array WORK. LWORK >= fla_max(1,M,N). */
/* > For optimum performance LWORK >= (M+N)*NB, where NB */
/* > is the optimal blocksize. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup doubleGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrices Q and P are represented as products of elementary */
/* > reflectors: */
/* > */
/* > If m >= n, */
/* > */
/* > Q = H(1) H(2) . . . H(n) and P = G(1) G(2) . . . G(n-1) */
/* > */
/* > Each H(i) and G(i) has the form: */
/* > */
/* > H(i) = I - tauq * v * v**T and G(i) = I - taup * u * u**T */
/* > */
/* > where tauq and taup are real scalars, and v and u are real vectors;
 */
/* > v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
 */
/* > u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
 */
/* > tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* > If m < n, */
/* > */
/* > Q = H(1) H(2) . . . H(m-1) and P = G(1) G(2) . . . G(m) */
/* > */
/* > Each H(i) and G(i) has the form: */
/* > */
/* > H(i) = I - tauq * v * v**T and G(i) = I - taup * u * u**T */
/* > */
/* > where tauq and taup are real scalars, and v and u are real vectors;
 */
/* > v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
 */
/* > u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
 */
/* > tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* > The contents of A on exit are illustrated by the following examples: */
/* > */
/* > m = 6 and n = 5 (m > n): m = 5 and n = 6 (m < n): */
/* > */
/* > ( d e u1 u1 u1 ) ( d u1 u1 u1 u1 u1 ) */
/* > ( v1 d e u2 u2 ) ( e d u2 u2 u2 u2 ) */
/* > ( v1 v2 d e u3 ) ( v1 e d u3 u3 u3 ) */
/* > ( v1 v2 v3 d e ) ( v1 v2 e d u4 u4 ) */
/* > ( v1 v2 v3 v4 d ) ( v1 v2 v3 e d u5 ) */
/* > ( v1 v2 v3 v4 v5 ) */
/* > */
/* > where d and e denote diagonal and off-diagonal elements of B, vi */
/* > denotes an element of the vector defining H(i), and ui an element of */
/* > the vector defining G(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int lapack_dgebrd(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                  doublereal *d__, doublereal *e, doublereal *tauq, doublereal *taup,
                  doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    aocl_int64_t i__, j, nb = 0, nx, ws;
    aocl_int64_t nbmin, iinfo, minmn;
    extern /* Subroutine */
        int
        lapack_dgebd2(aocl_int64_t *, aocl_int64_t *, doublereal *, aocl_int64_t *, doublereal *,
                      doublereal *, doublereal *, doublereal *, doublereal *, aocl_int64_t *);
    extern void fla_dlabrd(aocl_int64_t * m, aocl_int64_t * n, aocl_int64_t * nb, doublereal * a,
                           aocl_int64_t * lda, doublereal * d__, doublereal * e, doublereal * tauq,
                           doublereal * taup, doublereal * x, aocl_int64_t * ldx, doublereal * y,
                           aocl_int64_t * ldy);
#ifdef FLA_OPENMP_MULTITHREADING
    extern void fla_dlabrd_var1(
        aocl_int64_t * m, aocl_int64_t * n, aocl_int64_t * nb, doublereal * a, aocl_int64_t * lda,
        doublereal * d__, doublereal * e, doublereal * tauq, doublereal * taup, doublereal * x,
        aocl_int64_t * ldx, doublereal * y, aocl_int64_t * ldy, volatile FLA_BARRIER * barrier,
        doublereal * gemv_a_row_buffer, int requested_num_threads);
#endif
    aocl_int64_t lwkmin, ldwrkx, ldwrky, lwkopt;
    logical lquery;
#ifdef FLA_OPENMP_MULTITHREADING
    volatile FLA_BARRIER barrier;
    int dlabrd_threads = 1;
#endif
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    --work;
    /* Function Body */
    *info = 0;
    minmn = fla_min(*m, *n);
    if(minmn == 0)
    {
        lwkmin = 1;
        lwkopt = 1;
    }
    else
    {
        lwkmin = fla_max(*m, *n);
        /* Computing MAX */
        i__1 = 1;
#ifdef FLA_ENABLE_AMD_OPT
        i__2 = 110;
#else
        static aocl_int64_t c__1 = 1;
        i__2 = aocl_lapack_ilaenv(&c__1, "DGEBRD", " ", m, n, &c_n1, &c_n1); // , expr subst
#endif
        nb = fla_max(i__1, i__2);
        lwkopt = (*m + *n) * nb;
    }
    work[1] = (doublereal)lwkopt;
    lquery = *lwork == -1;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }
    else if(*lwork < lwkmin && !lquery)
    {
        *info = -10;
    }
    if(*info < 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DGEBRD", &i__1, (ftnlen)6);
        return 0;
    }
    else if(lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    minmn = fla_min(*m, *n);
    if(minmn == 0)
    {
        work[1] = 1.;
        return 0;
    }
    ws = fla_max(*m, *n);
    ldwrkx = *m;
    ldwrky = *n;
    if(nb > 1 && nb < minmn)
    {
        /* Set the crossover point NX. */
        /* Computing MAX */
        i__1 = nb;
        i__2 = aocl_lapack_ilaenv(&c__3, "DGEBRD", " ", m, n, &c_n1, &c_n1); // , expr subst
        nx = fla_max(i__1, i__2);
        /* Determine when to switch from blocked to unblocked code. */
        if(nx < minmn)
        {
            ws = lwkopt;
            if(*lwork < ws)
            {
                /* Not enough work space for the optimal NB, consider using */
                /* a smaller block size. */
                nbmin = aocl_lapack_ilaenv(&c__2, "DGEBRD", " ", m, n, &c_n1, &c_n1);
                if(*lwork >= (*m + *n) * nbmin)
                {
                    nb = *lwork / (*m + *n);
                }
                else
                {
                    nb = 1;
                    nx = minmn;
                }
            }
        }
    }
    else
    {
        nx = minmn;
    }
    i__1 = minmn - nx;
    i__2 = nb;
    i__3 = 0;
    i__ = 1;
    i__4 = 0;
#ifdef FLA_OPENMP_MULTITHREADING
    doublereal *gemv_a_row_buffer = NULL;

    if(i__1 > 0 && *m >= *n)
    {
        dlabrd_threads = fla_get_dlabrd_threads(*m, *n, fla_thread_get_num_threads());
        /* fla_dlabrd_var1 (m >= n, OpenMP) uses scratch as: [0..nb-1], [nb..2nb-1],
         * memcpy staging at [2*nb + i__5] (needs up to index 3*nb - 2), and
         * [(*n_sub) + i__5] for X staging (needs up to *n_sub + nb - 1). */
        i__3 = fla_max(3 * nb, fla_max(*m, *n) + nb);
        gemv_a_row_buffer = (doublereal *)malloc((size_t)i__3 * sizeof(doublereal));
    }
#endif

    for(i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
    {
        /* Reduce rows and columns i:i+nb-1 to bidiagonal form and return */
        /* the matrices X and Y which are needed to update the unreduced */
        /* part of the matrix */
        i__3 = *m - i__ + 1;
        i__4 = *n - i__ + 1;
#ifdef FLA_OPENMP_MULTITHREADING
        if(gemv_a_row_buffer != NULL)
        {
            fla_dlabrd_var1(&i__3, &i__4, &nb, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__],
                            &tauq[i__], &taup[i__], &work[1], &ldwrkx, &work[ldwrkx * nb + 1],
                            &ldwrky, &barrier, &gemv_a_row_buffer[0], dlabrd_threads);
        }
        else
#endif
        {
            fla_dlabrd(&i__3, &i__4, &nb, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__],
                       &tauq[i__], &taup[i__], &work[1], &ldwrkx, &work[ldwrkx * nb + 1], &ldwrky);
        }
        /* Update the trailing submatrix A(i+nb:m,i+nb:n), using an update */
        /* of the form A := A - V*Y**T - X*U**T */
        i__3 = *m - i__ - nb + 1;
        i__4 = *n - i__ - nb + 1;
        aocl_blas_dgemm("No transpose", "Transpose", &i__3, &i__4, &nb, &c_b21,
                        &a[i__ + nb + i__ * a_dim1], lda, &work[ldwrkx * nb + nb + 1], &ldwrky,
                        &c_b22, &a[i__ + nb + (i__ + nb) * a_dim1], lda);

        i__3 = *m - i__ - nb + 1;
        i__4 = *n - i__ - nb + 1;

        aocl_blas_dgemm("No transpose", "No transpose", &i__3, &i__4, &nb, &c_b21, &work[nb + 1],
                        &ldwrkx, &a[i__ + (i__ + nb) * a_dim1], lda, &c_b22,
                        &a[i__ + nb + (i__ + nb) * a_dim1], lda);
        /* Copy diagonal and off-diagonal elements of B back into A */
        if(*m >= *n)
        {
            i__3 = i__ + nb - 1;
            for(j = i__; j <= i__3; ++j)
            {
                a[j + j * a_dim1] = d__[j];
                a[j + (j + 1) * a_dim1] = e[j];
                /* L10: */
            }
        }
        else
        {
            i__3 = i__ + nb - 1;
            for(j = i__; j <= i__3; ++j)
            {
                a[j + j * a_dim1] = d__[j];
                a[j + 1 + j * a_dim1] = e[j];
                /* L20: */
            }
        }
        /* L30: */
    }

#ifdef FLA_OPENMP_MULTITHREADING
    if(gemv_a_row_buffer != NULL)
    {
        free(gemv_a_row_buffer);
    }
#endif
    /* Use unblocked code to reduce the remainder of the matrix */
    i__2 = *m - i__ + 1;
    i__1 = *n - i__ + 1;
    lapack_dgebd2(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], &tauq[i__],
                  &taup[i__], &work[1], &iinfo);
    work[1] = (doublereal)ws;
    return 0;
    /* End of DGEBRD */
}
/* lapack_dgebrd */