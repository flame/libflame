/*
 *     Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
 */

/* dorgqr.f -- translated by f2c (version 20160102). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLAME.h"
#include "FLA_f2c.h" /* Table of constant values */
#if FLA_ENABLE_AOCL_BLAS
#include <blis.h>
#endif
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static aocl_int64_t c__3 = 3;
static aocl_int64_t c__2 = 2;

#ifdef FLA_ENABLE_AMD_OPT
extern int fla_thread_get_num_threads(void);

#if defined(FLA_ENABLE_MULTITHREADING) || defined(FLA_OPENMP_MULTITHREADING)
/* Return nonzero when runtime BLIS threading is active, so DORGQR can use the
 * multi-thread tuned path instead of the single-thread tuning. */
static int dorgqr_use_threaded_tuning(void)
{
#if FLA_ENABLE_AOCL_BLAS
    return bli_thread_get_num_threads() > 1;
#else
    return 0;
#endif
}
#endif

/* Select AOCL-tuned DORGQR block sizes for shape and thread-count ranges that
 * benchmark better than the generic ILAENV value. */
static aocl_int64_t dorgqr_tuned_nb(aocl_int64_t nb, aocl_int64_t m, aocl_int64_t n,
                                    aocl_int64_t num_threads)
{
    if(num_threads == 1)
    {
        if(m == n)
        {
            if(n >= 100 && n <= 1000)
            {
                return 32;
            }
            else if(n <= 2000)
            {
                return 64;
            }
            else if(n <= 4000)
            {
                return 96;
            }
            else if(n <= 12000)
            {
                return 128;
            }
            else if(n <= 16000)
            {
                return 160;
            }
        }
        else if(m >= 100 && m <= 500 && n < 100)
        {
            return 32;
        }
    }
    else
    {
        if(m == n)
        {
            if(n > 2500 && n <= 3000)
            {
                return 48;
            }
            else if(n > 3000)
            {
                return 96;
            }
        }
        else if(m > 4000 && m <= 4500 && n > 3500 && n <= 4000)
        {
            return 96;
        }
    }

    return nb;
}

#if FLA_ENABLE_AOCL_BLAS \
    && (defined(FLA_ENABLE_MULTITHREADING) || defined(FLA_OPENMP_MULTITHREADING))
/* Cap BLIS threads for DORGQR shapes where using all requested threads is
 * slower than a smaller tuned thread team. */
static aocl_int64_t dorgqr_tuned_threads(aocl_int64_t m, aocl_int64_t n, aocl_int64_t num_threads)
{
    aocl_int64_t thread_cap = num_threads;

    if(m == n)
    {
        if(n > 6000)
        {
            thread_cap = 24;
        }
        else if(n > 4000)
        {
            thread_cap = 8;
        }
        else if(n > 3200)
        {
            thread_cap = 2;
        }
        else if(n >= 3000)
        {
            thread_cap = 4;
        }
        else if(n > 1000)
        {
            thread_cap = 8;
        }
    }
    else if(m >= 700 && m <= 1300 && n <= 200)
    {
        thread_cap = 1;
    }
    else if(n > 500 && n <= 650 && ((m >= 600 && m <= 750) || (m > 3000 && m <= 4000)))
    {
        thread_cap = 8;
    }
    else if(m > 4000 && m <= 4500 && n > 3500 && n <= 4000)
    {
        thread_cap = 8;
    }

    return thread_cap < num_threads ? thread_cap : num_threads;
}
#endif // #if FLA_ENABLE_AOCL_BLAS && (defined(FLA_ENABLE_MULTITHREADING) ||
       // defined(FLA_OPENMP_MULTITHREADING))
#endif // #ifdef FLA_ENABLE_AMD_OPT

/* > \brief \b DORGQR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DORGQR + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgqr.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgqr.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgqr.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, K, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORGQR generates an M-by-N real matrix Q with orthonormal columns, */
/* > which is defined as the first N columns of a product of K elementary */
/* > reflectors of order M */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by DGEQRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix Q. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix Q. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines the */
/* > matrix Q. N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the i-th column must contain the vector which */
/* > defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* > returned by DGEQRF in the first k columns of its array */
/* > argument A. */
/* > On exit, the M-by-N matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The first dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by DGEQRF. */
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
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
/* > For optimum performance LWORK >= N*NB, where NB is the */
/* > optimal blocksize. */
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
/* > < 0: if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int lapack_dorgqr(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a,
                  aocl_int64_t *lda, doublereal *tau, doublereal *work, aocl_int64_t *lwork,
                  aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    aocl_int64_t i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    extern void dorg2r_fla(aocl_int64_t *, aocl_int64_t *, aocl_int64_t *, doublereal *,
                           aocl_int64_t *, doublereal *, doublereal *, aocl_int64_t *);
    aocl_int64_t ldwork, lwkopt;
    logical lquery;
#if defined(FLA_ENABLE_AMD_OPT)                                                   \
    && (defined(FLA_ENABLE_MULTITHREADING) || defined(FLA_OPENMP_MULTITHREADING)) \
    && FLA_ENABLE_AOCL_BLAS
    aocl_int64_t orig_blis_threads = 0;
    aocl_int64_t tuned_blis_threads = 0;
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
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    nb = 0;
#ifdef FLA_ENABLE_AMD_OPT
    /* precomputed workspace size */
    if(*n == 1)
    {
        work[1] = 32;
    }
    else if(*n <= 6)
    {
        work[1] = 192;
    }
    else
    {
        nb = aocl_lapack_ilaenv(&c__1, "DORGQR", " ", m, n, k, &c_n1);
#if defined(FLA_ENABLE_MULTITHREADING) || defined(FLA_OPENMP_MULTITHREADING)
        if((*n > 64) && (*n == *k) && dorgqr_use_threaded_tuning())
        {
            nb = dorgqr_tuned_nb(nb, *m, *n, 2);
        }
#else
        if(*n == *k)
        {
            nb = dorgqr_tuned_nb(nb, *m, *n, 1);
        }
#endif
        lwkopt = fla_max(1, *n) * nb;
        work[1] = (doublereal)lwkopt;
    }
#else
    nb = aocl_lapack_ilaenv(&c__1, "DORGQR", " ", m, n, k, &c_n1);
    lwkopt = fla_max(1, *n) * nb;
    work[1] = (doublereal)lwkopt;
#endif
    lquery = *lwork == -1;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0 || *n > *m)
    {
        *info = -2;
    }
    else if(*k < 0 || *k > *n)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -5;
    }
    else if(*lwork < fla_max(1, *n) && !lquery)
    {
        *info = -8;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DORGQR", &i__1, (ftnlen)6);
        return 0;
    }
    else if(lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if(*n <= 0)
    {
        work[1] = 1.;
        return 0;
    }
#if defined(FLA_ENABLE_AMD_OPT)                                                   \
    && (defined(FLA_ENABLE_MULTITHREADING) || defined(FLA_OPENMP_MULTITHREADING)) \
    && FLA_ENABLE_AOCL_BLAS
    if((*n > 64) && (*n == *k))
    {
        orig_blis_threads = bli_thread_get_num_threads();
        tuned_blis_threads = dorgqr_tuned_threads(*m, *n, orig_blis_threads);
        if(tuned_blis_threads < orig_blis_threads)
        {
            bli_thread_set_num_threads(tuned_blis_threads);
        }
    }
#endif
    nbmin = 2;
    nx = 0;
    iws = *n;
    if(nb > 1 && nb < *k)
    {
        /* Determine when to cross over from blocked to unblocked code. */
        /* Computing MAX */
        i__1 = 0;
        i__2 = aocl_lapack_ilaenv(&c__3, "DORGQR", " ", m, n, k, &c_n1); // , expr subst
        nx = fla_max(i__1, i__2);
        if(nx < *k)
        {
            /* Determine if workspace is large enough for blocked code. */
            ldwork = *n;
            iws = ldwork * nb;
            if(*lwork < iws)
            {
                /* Not enough workspace to use optimal NB: reduce NB and */
                /* determine the minimum value of NB. */
                nb = *lwork / ldwork;
                /* Computing MAX */
                i__1 = 2;
                i__2 = aocl_lapack_ilaenv(&c__2, "DORGQR", " ", m, n, k, &c_n1); // , expr subst
                nbmin = fla_max(i__1, i__2);
            }
        }
    }
    if(nb >= nbmin && nb < *k && nx < *k)
    {
        /* Use blocked code after the last block. */
        /* The first kk columns are handled by the block method. */
        ki = (*k - nx - 1) / nb * nb;
        /* Computing MIN */
        i__2 = ki + nb; // , expr subst
        kk = fla_min(*k, i__2);
        /* Set A(1:kk,kk+1:n) to zero. */
        i__1 = *n;
        for(j = kk + 1; j <= i__1; ++j)
        {
            for(i__ = 1; i__ <= kk; ++i__)
            {
                a[i__ + j * a_dim1] = 0.;
                /* L10: */
            }
            /* L20: */
        }
    }
    else
    {
        kk = 0;
    }
    /* Use unblocked code for the last or only block. */
    if(kk < *n)
    {
        i__1 = *m - kk;
        i__2 = *n - kk;
        i__3 = *k - kk;
        dorg2r_fla(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &tau[kk + 1], &work[1],
                   &iinfo);
    }
    if(kk > 0)
    {
        /* Use blocked code */
        i__1 = -nb;
        for(i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1)
        {
            /* Computing MIN */
            i__3 = *k - i__ + 1; // , expr subst
            ib = fla_min(nb, i__3);
            if(i__ + ib <= *n)
            {
                /* Form the triangular factor of the block reflector */
                /* H = H(i) H(i+1) . . . H(i+ib-1) */
                i__2 = *m - i__ + 1;
                aocl_lapack_dlarft("Forward", "Columnwise", &i__2, &ib, &a[i__ + i__ * a_dim1], lda,
                                   &tau[i__], &work[1], &ldwork);
                /* Apply H to A(i:m,i+ib:n) from the left */
                i__3 = *n - i__ - ib + 1;
                aocl_lapack_dlarfb("Left", "No transpose", "Forward", "Columnwise", &i__2, &i__3,
                                   &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &ldwork,
                                   &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib + 1], &ldwork);
            }
            /* Apply H to rows i:m of current block */
            i__2 = *m - i__ + 1;
            dorg2r_fla(&i__2, &ib, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);
            /* Set rows 1:i-1 of current block to zero */
            i__2 = i__ + ib - 1;
            for(j = i__; j <= i__2; ++j)
            {
                i__3 = i__ - 1;
                for(l = 1; l <= i__3; ++l)
                {
                    a[l + j * a_dim1] = 0.;
                    /* L30: */
                }
                /* L40: */
            }
            /* L50: */
        }
    }
#if defined(FLA_ENABLE_AMD_OPT)                                                   \
    && (defined(FLA_ENABLE_MULTITHREADING) || defined(FLA_OPENMP_MULTITHREADING)) \
    && FLA_ENABLE_AOCL_BLAS
    if((*n == *k) && (tuned_blis_threads < orig_blis_threads))
    {
        bli_thread_set_num_threads(orig_blis_threads);
    }
#endif
    work[1] = (doublereal)iws;
    return 0;
    /* End of DORGQR */
}
/* lapack_dorgqr */
