/* ../netlib/cgetc2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static complex c_b10 =
{
    -1.f,-0.f
}
;
/* > \brief \b CGETC2 computes the LU factorization with complete pivoting of the general n-by-n matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGETC2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetc2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetc2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetc2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGETC2( N, A, LDA, IPIV, JPIV, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), JPIV( * ) */
/* COMPLEX A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGETC2 computes an LU factorization, using complete pivoting, of the */
/* > n-by-n matrix A. The factorization has the form A = P * L * U * Q, */
/* > where P and Q are permutation matrices, L is lower triangular with */
/* > unit diagonal elements and U is upper triangular. */
/* > */
/* > This is a level 1 BLAS version of the algorithm. */
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
/* > A is COMPLEX array, dimension (LDA, N) */
/* > On entry, the n-by-n matrix to be factored. */
/* > On exit, the factors L and U from the factorization */
/* > A = P*L*U*Q;
the unit diagonal elements of L are not stored. */
/* > If U(k, k) appears to be less than SMIN, U(k, k) is given the */
/* > value of SMIN, giving a nonsingular perturbed system. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1, N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N). */
/* > The pivot indices;
for 1 <= i <= N, row i of the */
/* > matrix has been interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] JPIV */
/* > \verbatim */
/* > JPIV is INTEGER array, dimension (N). */
/* > The pivot indices;
for 1 <= j <= N, column j of the */
/* > matrix has been interchanged with column JPIV(j). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > > 0: if INFO = k, U(k, k) is likely to produce overflow if */
/* > one tries to solve for x in Ax = b. So U is perturbed */
/* > to avoid the overflow. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup complexGEauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* ===================================================================== */
/* Subroutine */
int cgetc2_(integer *n, complex *a, integer *lda, integer * ipiv, integer *jpiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;
    complex q__1;
    /* Builtin functions */
    double c_abs(complex *);
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    integer i__, j, ip, jp;
    real eps;
    integer ipv, jpv;
    real smin, xmax;
    extern /* Subroutine */
    int cgeru_(integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, integer *), cswap_(integer *, complex *, integer *, complex *, integer *), slabad_(real *, real *);
    extern real slamch_(char *);
    real bignum, smlnum;
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
    eps = slamch_("P");
    smlnum = slamch_("S") / eps;
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);
    /* Factorize A using complete pivoting. */
    /* Set pivots less than SMIN to SMIN */
    i__1 = *n - 1;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        /* Find max element in matrix A */
        xmax = 0.f;
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
                if (c_abs(&a[ip + jp * a_dim1]) >= xmax)
                {
                    xmax = c_abs(&a[ip + jp * a_dim1]);
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
            r__1 = eps * xmax;
            smin = max(r__1,smlnum);
        }
        /* Swap rows */
        if (ipv != i__)
        {
            cswap_(n, &a[ipv + a_dim1], lda, &a[i__ + a_dim1], lda);
        }
        ipiv[i__] = ipv;
        /* Swap columns */
        if (jpv != i__)
        {
            cswap_(n, &a[jpv * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], & c__1);
        }
        jpiv[i__] = jpv;
        /* Check for singularity */
        if (c_abs(&a[i__ + i__ * a_dim1]) < smin)
        {
            *info = i__;
            i__2 = i__ + i__ * a_dim1;
            q__1.r = smin;
            q__1.i = 0.f; // , expr subst
            a[i__2].r = q__1.r;
            a[i__2].i = q__1.i; // , expr subst
        }
        i__2 = *n;
        for (j = i__ + 1;
                j <= i__2;
                ++j)
        {
            i__3 = j + i__ * a_dim1;
            c_div(&q__1, &a[j + i__ * a_dim1], &a[i__ + i__ * a_dim1]);
            a[i__3].r = q__1.r;
            a[i__3].i = q__1.i; // , expr subst
            /* L30: */
        }
        i__2 = *n - i__;
        i__3 = *n - i__;
        cgeru_(&i__2, &i__3, &c_b10, &a[i__ + 1 + i__ * a_dim1], &c__1, &a[ i__ + (i__ + 1) * a_dim1], lda, &a[i__ + 1 + (i__ + 1) * a_dim1], lda);
        /* L40: */
    }
    if (c_abs(&a[*n + *n * a_dim1]) < smin)
    {
        *info = *n;
        i__1 = *n + *n * a_dim1;
        q__1.r = smin;
        q__1.i = 0.f; // , expr subst
        a[i__1].r = q__1.r;
        a[i__1].i = q__1.i; // , expr subst
    }
    /* Set last pivots to N */
    ipiv[*n] = *n;
    jpiv[*n] = *n;
    return 0;
    /* End of CGETC2 */
}
/* cgetc2_ */
