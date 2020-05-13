/* ../netlib/dgetc2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b10 = -1.;
/* > \brief \b DGETC2 computes the LU factorization with complete pivoting of the general n-by-n matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGETC2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetc2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetc2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetc2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGETC2( N, A, LDA, IPIV, JPIV, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), JPIV( * ) */
/* DOUBLE PRECISION A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGETC2 computes an LU factorization with complete pivoting of the */
/* > n-by-n matrix A. The factorization has the form A = P * L * U * Q, */
/* > where P and Q are permutation matrices, L is lower triangular with */
/* > unit diagonal elements and U is upper triangular. */
/* > */
/* > This is the Level 2 BLAS algorithm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA, N) */
/* > On entry, the n-by-n matrix A to be factored. */
/* > On exit, the factors L and U from the factorization */
/* > A = P*L*U*Q;
the unit diagonal elements of L are not stored. */
/* > If U(k, k) appears to be less than SMIN, U(k, k) is given the */
/* > value of SMIN, i.e., giving a nonsingular perturbed system. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension(N). */
/* > The pivot indices;
for 1 <= i <= N, row i of the */
/* > matrix has been interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] JPIV */
/* > \verbatim */
/* > JPIV is INTEGER array, dimension(N). */
/* > The pivot indices;
for 1 <= j <= N, column j of the */
/* > matrix has been interchanged with column JPIV(j). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > > 0: if INFO = k, U(k, k) is likely to produce owerflow if */
/* > we try to solve for x in Ax = b. So U is perturbed to */
/* > avoid the overflow. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup doubleGEauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* ===================================================================== */
/* Subroutine */
int dgetc2_(integer *n, doublereal *a, integer *lda, integer *ipiv, integer *jpiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    /* Local variables */
    integer i__, j, ip, jp;
    doublereal eps;
    integer ipv, jpv;
    extern /* Subroutine */
    int dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal smin, xmax;
    extern /* Subroutine */
    int dswap_(integer *, doublereal *, integer *, doublereal *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *);
    doublereal bignum, smlnum;
    /* -- LAPACK auxiliary routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Set constants to control overflow */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --jpiv;
    /* Function Body */
    *info = 0;
    eps = dlamch_("P");
    smlnum = dlamch_("S") / eps;
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    /* Factorize A using complete pivoting. */
    /* Set pivots less than SMIN to SMIN. */
    i__1 = *n - 1;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        /* Find max element in matrix A */
        xmax = 0.;
        i__2 = *n;
        for (ip = i__;
                ip <= i__2;
                ++ip)
        {
            i__3 = *n;
            for (jp = i__;
                    jp <= i__3;
                    ++jp)
            {
                if ((d__1 = a[ip + jp * a_dim1], f2c_abs(d__1)) >= xmax)
                {
                    xmax = (d__1 = a[ip + jp * a_dim1], f2c_abs(d__1));
                    ipv = ip;
                    jpv = jp;
                }
                /* L10: */
            }
            /* L20: */
        }
        if (i__ == 1)
        {
            /* Computing MAX */
            d__1 = eps * xmax;
            smin = max(d__1,smlnum);
        }
        /* Swap rows */
        if (ipv != i__)
        {
            dswap_(n, &a[ipv + a_dim1], lda, &a[i__ + a_dim1], lda);
        }
        ipiv[i__] = ipv;
        /* Swap columns */
        if (jpv != i__)
        {
            dswap_(n, &a[jpv * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], & c__1);
        }
        jpiv[i__] = jpv;
        /* Check for singularity */
        if ((d__1 = a[i__ + i__ * a_dim1], f2c_abs(d__1)) < smin)
        {
            *info = i__;
            a[i__ + i__ * a_dim1] = smin;
        }
        i__2 = *n;
        for (j = i__ + 1;
                j <= i__2;
                ++j)
        {
            a[j + i__ * a_dim1] /= a[i__ + i__ * a_dim1];
            /* L30: */
        }
        i__2 = *n - i__;
        i__3 = *n - i__;
        dger_(&i__2, &i__3, &c_b10, &a[i__ + 1 + i__ * a_dim1], &c__1, &a[i__ + (i__ + 1) * a_dim1], lda, &a[i__ + 1 + (i__ + 1) * a_dim1], lda);
        /* L40: */
    }
    if ((d__1 = a[*n + *n * a_dim1], f2c_abs(d__1)) < smin)
    {
        *info = *n;
        a[*n + *n * a_dim1] = smin;
    }
    /* Set last pivots to N */
    ipiv[*n] = *n;
    jpiv[*n] = *n;
    return 0;
    /* End of DGETC2 */
}
/* dgetc2_ */
