/* ../netlib/dorgrq.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
/* > \brief \b DORGRQ */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DORGRQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgrq. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgrq. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgrq. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */
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
/* > DORGRQ generates an M-by-N real matrix Q with orthonormal rows, */
/* > which is defined as the last M rows of a product of K elementary */
/* > reflectors of order N */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by DGERQF. */
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
/* > The number of columns of the matrix Q. N >= M. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines the */
/* > matrix Q. M >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the (m-k+i)-th row must contain the vector which */
/* > defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* > returned by DGERQF in the last k rows of its array argument */
/* > A. */
/* > On exit, the M-by-N matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The first dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by DGERQF. */
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
/* > The dimension of the array WORK. LWORK >= max(1,M). */
/* > For optimum performance LWORK >= M*NB, where NB is the */
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
/* > \date November 2011 */
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int dorgrq_(integer *m, integer *n, integer *k, doublereal * a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer i__, j, l, ib, nb, ii, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */
    int dorgr2_(integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *), dlarft_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer ldwork, lwkopt;
    logical lquery;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < *m)
    {
        *info = -2;
    }
    else if (*k < 0 || *k > *m)
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    if (*info == 0)
    {
        if (*m <= 0)
        {
            lwkopt = 1;
        }
        else
        {
            nb = ilaenv_(&c__1, "DORGRQ", " ", m, n, k, &c_n1);
            lwkopt = *m * nb;
        }
        work[1] = (doublereal) lwkopt;
        if (*lwork < max(1,*m) && ! lquery)
        {
            *info = -8;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DORGRQ", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*m <= 0)
    {
        return 0;
    }
    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < *k)
    {
        /* Determine when to cross over from blocked to unblocked code. */
        /* Computing MAX */
        i__1 = 0;
        i__2 = ilaenv_(&c__3, "DORGRQ", " ", m, n, k, &c_n1); // , expr subst
        nx = max(i__1,i__2);
        if (nx < *k)
        {
            /* Determine if workspace is large enough for blocked code. */
            ldwork = *m;
            iws = ldwork * nb;
            if (*lwork < iws)
            {
                /* Not enough workspace to use optimal NB: reduce NB and */
                /* determine the minimum value of NB. */
                nb = *lwork / ldwork;
                /* Computing MAX */
                i__1 = 2;
                i__2 = ilaenv_(&c__2, "DORGRQ", " ", m, n, k, &c_n1); // , expr subst
                nbmin = max(i__1,i__2);
            }
        }
    }
    if (nb >= nbmin && nb < *k && nx < *k)
    {
        /* Use blocked code after the first block. */
        /* The last kk rows are handled by the block method. */
        /* Computing MIN */
        i__1 = *k;
        i__2 = (*k - nx + nb - 1) / nb * nb; // , expr subst
        kk = min(i__1,i__2);
        /* Set A(1:m-kk,n-kk+1:n) to zero. */
        i__1 = *n;
        for (j = *n - kk + 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m - kk;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
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
    /* Use unblocked code for the first or only block. */
    i__1 = *m - kk;
    i__2 = *n - kk;
    i__3 = *k - kk;
    dorgr2_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &iinfo) ;
    if (kk > 0)
    {
        /* Use blocked code */
        i__1 = *k;
        i__2 = nb;
        for (i__ = *k - kk + 1;
                i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
                i__ += i__2)
        {
            /* Computing MIN */
            i__3 = nb;
            i__4 = *k - i__ + 1; // , expr subst
            ib = min(i__3,i__4);
            ii = *m - *k + i__;
            if (ii > 1)
            {
                /* Form the triangular factor of the block reflector */
                /* H = H(i+ib-1) . . . H(i+1) H(i) */
                i__3 = *n - *k + i__ + ib - 1;
                dlarft_("Backward", "Rowwise", &i__3, &ib, &a[ii + a_dim1], lda, &tau[i__], &work[1], &ldwork);
                /* Apply H**T to A(1:m-k+i-1,1:n-k+i+ib-1) from the right */
                i__3 = ii - 1;
                i__4 = *n - *k + i__ + ib - 1;
                dlarfb_("Right", "Transpose", "Backward", "Rowwise", &i__3, & i__4, &ib, &a[ii + a_dim1], lda, &work[1], &ldwork, & a[a_offset], lda, &work[ib + 1], &ldwork);
            }
            /* Apply H**T to columns 1:n-k+i+ib-1 of current block */
            i__3 = *n - *k + i__ + ib - 1;
            dorgr2_(&ib, &i__3, &ib, &a[ii + a_dim1], lda, &tau[i__], &work[1] , &iinfo);
            /* Set columns n-k+i+ib:n of current block to zero */
            i__3 = *n;
            for (l = *n - *k + i__ + ib;
                    l <= i__3;
                    ++l)
            {
                i__4 = ii + ib - 1;
                for (j = ii;
                        j <= i__4;
                        ++j)
                {
                    a[j + l * a_dim1] = 0.;
                    /* L30: */
                }
                /* L40: */
            }
            /* L50: */
        }
    }
    work[1] = (doublereal) iws;
    return 0;
    /* End of DORGRQ */
}
/* dorgrq_ */
