/* ../netlib/zung2r.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZUNG2R */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZUNG2R + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zung2r. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zung2r. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zung2r. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, K, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNG2R generates an m by n complex matrix Q with orthonormal columns, */
/* > which is defined as the first n columns of a product of k elementary */
/* > reflectors of order m */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by ZGEQRF. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the i-th column must contain the vector which */
/* > defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* > returned by ZGEQRF in the first k columns of its array */
/* > argument A. */
/* > On exit, the m by n matrix Q. */
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
/* > TAU is COMPLEX*16 array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by ZGEQRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N) */
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
/* > \ingroup complex16OTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int zung2r_(integer *m, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex * work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;
    /* Local variables */
    integer i__, j, l;
    extern /* Subroutine */
    int zscal_(integer *, doublecomplex *, doublecomplex *, integer *), zlarf_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *), xerbla_(char *, integer *);
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
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0 || *n > *m)
    {
        *info = -2;
    }
    else if (*k < 0 || *k > *n)
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZUNG2R", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n <= 0)
    {
        return 0;
    }
    /* Initialise columns k+1:n to columns of the unit matrix */
    i__1 = *n;
    for (j = *k + 1;
            j <= i__1;
            ++j)
    {
        i__2 = *m;
        for (l = 1;
                l <= i__2;
                ++l)
        {
            i__3 = l + j * a_dim1;
            a[i__3].r = 0.;
            a[i__3].i = 0.; // , expr subst
            /* L10: */
        }
        i__2 = j + j * a_dim1;
        a[i__2].r = 1.;
        a[i__2].i = 0.; // , expr subst
        /* L20: */
    }
    for (i__ = *k;
            i__ >= 1;
            --i__)
    {
        /* Apply H(i) to A(i:m,i:n) from the left */
        if (i__ < *n)
        {
            i__1 = i__ + i__ * a_dim1;
            a[i__1].r = 1.;
            a[i__1].i = 0.; // , expr subst
            i__1 = *m - i__ + 1;
            i__2 = *n - i__;
            zlarf_("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[ i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
        }
        if (i__ < *m)
        {
            i__1 = *m - i__;
            i__2 = i__;
            z__1.r = -tau[i__2].r;
            z__1.i = -tau[i__2].i; // , expr subst
            zscal_(&i__1, &z__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
        }
        i__1 = i__ + i__ * a_dim1;
        i__2 = i__;
        z__1.r = 1. - tau[i__2].r;
        z__1.i = 0. - tau[i__2].i; // , expr subst
        a[i__1].r = z__1.r;
        a[i__1].i = z__1.i; // , expr subst
        /* Set A(1:i-1,i) to zero */
        i__1 = i__ - 1;
        for (l = 1;
                l <= i__1;
                ++l)
        {
            i__2 = l + i__ * a_dim1;
            a[i__2].r = 0.;
            a[i__2].i = 0.; // , expr subst
            /* L30: */
        }
        /* L40: */
    }
    return 0;
    /* End of ZUNG2R */
}
/* zung2r_ */
