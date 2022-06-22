/* ../netlib/v3.9.0/dgelqt.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DGELQT */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGEQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelqt. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelqt. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelqt. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGELQT( M, N, MB, A, LDA, T, LDT, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDT, M, N, MB */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), T( LDT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGELQT computes a blocked LQ factorization of a real M-by-N matrix A */
/* > using the compact WY representation of Q. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* > MB is INTEGER */
/* > The block size to be used in the blocked QR. MIN(M,N) >= MB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the elements on and below the diagonal of the array */
/* > contain the M-by-MIN(M,N) lower trapezoidal matrix L (L is */
/* > lower triangular if M <= N);
the elements above the diagonal */
/* > are the rows of V. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is DOUBLE PRECISION array, dimension (LDT,MIN(M,N)) */
/* > The upper triangular block reflectors stored in compact form */
/* > as a sequence of upper triangular blocks. See below */
/* > for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= MB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MB*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup doubleGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix V stores the elementary reflectors H(i) in the i-th row */
/* > above the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* > V = ( 1 v1 v1 v1 v1 ) */
/* > ( 1 v2 v2 v2 ) */
/* > ( 1 v3 v3 ) */
/* > */
/* > */
/* > where the vi's represent the vectors which define H(i), which are returned */
/* > in the matrix A. The 1's along the diagonal of V are not stored in A. */
/* > Let K=MIN(M,N). The number of blocks is B = ceiling(K/MB), where each */
/* > block is of order MB except for the last block, which is of order */
/* > IB = K - (B-1)*MB. For each of the B blocks, a upper triangular block */
/* > reflector factor is computed: T1, T2, ..., TB. The MB-by-MB (and IB-by-IB */
/* > for the last block) T's are stored in the MB-by-K matrix T as */
/* > */
/* > T = (T1 T2 ... TB). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int dgelqt_(integer *m, integer *n, integer *mb, doublereal * a, integer *lda, doublereal *t, integer *ldt, doublereal *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgelqt inputs: m %" FLA_IS ", n %" FLA_IS ", mb %" FLA_IS ", lda %" FLA_IS ", ldt %" FLA_IS "",*m, *n, *mb, *lda, *ldt);
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    integer i__, k, ib, iinfo;
    extern /* Subroutine */
    int dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *), xerbla_(char *, integer *), dgelqt3_(integer *, integer *, doublereal *, integer *, doublereal *, integer *, integer *);
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;
    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*mb < 1 || *mb > min(*m,*n) && min(*m,*n) > 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    else if (*ldt < *mb)
    {
        *info = -7;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGELQT", &i__1);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    k = min(*m,*n);
    if (k == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Blocked loop of length K */
    i__1 = k;
    i__2 = *mb;
    for (i__ = 1;
            i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
            i__ += i__2)
    {
        /* Computing MIN */
        i__3 = k - i__ + 1;
        ib = min(i__3,*mb);
        /* Compute the LQ factorization of the current block A(I:M,I:I+IB-1) */
        i__3 = *n - i__ + 1;
        dgelqt3_(&ib, &i__3, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &iinfo);
        if (i__ + ib <= *m)
        {
            /* Update by applying H**T to A(I:M,I+IB:N) from the right */
            i__3 = *m - i__ - ib + 1;
            i__4 = *n - i__ + 1;
            i__5 = *m - i__ - ib + 1;
            dlarfb_("R", "N", "F", "R", &i__3, &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &a[i__ + ib + i__ * a_dim1], lda, &work[1], &i__5);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of DGELQT */
}
/* dgelqt_ */
