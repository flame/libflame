/* ../netlib/dsyequb.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DSYEQUB */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSYEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyequb .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyequb .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyequb .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, N */
/* DOUBLE PRECISION AMAX, SCOND */
/* CHARACTER UPLO */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), S( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYEQUB computes row and column scalings intended to equilibrate a */
/* > symmetric matrix A and reduce its condition number */
/* > (with respect to the two-norm). S contains the scale factors, */
/* > S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with */
/* > elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal. This */
/* > choice of S puts the condition number of B within a factor N of the */
/* > smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**T;
*/
/* > = 'L': Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > The N-by-N symmetric matrix whose scaling */
/* > factors are to be computed. Only the diagonal elements of A */
/* > are referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is DOUBLE PRECISION array, dimension (N) */
/* > If INFO = 0, S contains the scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] SCOND */
/* > \verbatim */
/* > SCOND is DOUBLE PRECISION */
/* > If INFO = 0, S contains the ratio of the smallest S(i) to */
/* > the largest S(i). If SCOND >= 0.1 and AMAX is neither too */
/* > large nor too small, it is not worth scaling by S. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* > AMAX is DOUBLE PRECISION */
/* > Absolute value of largest matrix element. If AMAX is very */
/* > close to overflow or very close to underflow, the matrix */
/* > should be scaled. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the i-th diagonal element is nonpositive. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup doubleSYcomputational */
/* > \par References: */
/* ================ */
/* > */
/* > Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n */
/* > Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n */
/* > DOI 10.1023/B:NUMA.0000016606.32820.69 \n */
/* > Tech report version: http://ruready.utah.edu/archive/papers/bin.pdf */
/* > */
/* ===================================================================== */
/* Subroutine */
int dsyequb_(char *uplo, integer *n, doublereal *a, integer * lda, doublereal *s, doublereal *scond, doublereal *amax, doublereal * work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;
    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_di(doublereal *, integer *);
    /* Local variables */
    doublereal d__;
    integer i__, j;
    doublereal t, u, c0, c1, c2, si;
    logical up;
    doublereal avg, std, tol, base;
    integer iter;
    doublereal smin, smax, scale;
    extern logical lsame_(char *, char *);
    doublereal sumsq;
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal bignum;
    extern /* Subroutine */
    int dlassq_(integer *, doublereal *, integer *, doublereal *, doublereal *);
    doublereal smlnum;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    --work;
    /* Function Body */
    *info = 0;
    if (! (lsame_(uplo, "U") || lsame_(uplo, "L")))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < max(1,*n))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DSYEQUB", &i__1);
        return 0;
    }
    up = lsame_(uplo, "U");
    *amax = 0.;
    /* Quick return if possible. */
    if (*n == 0)
    {
        *scond = 1.;
        return 0;
    }
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        s[i__] = 0.;
    }
    *amax = 0.;
    if (up)
    {
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j - 1;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                /* Computing MAX */
                d__2 = s[i__];
                d__3 = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)); // , expr subst
                s[i__] = max(d__2,d__3);
                /* Computing MAX */
                d__2 = s[j];
                d__3 = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)); // , expr subst
                s[j] = max(d__2,d__3);
                /* Computing MAX */
                d__2 = *amax;
                d__3 = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)); // , expr subst
                *amax = max(d__2,d__3);
            }
            /* Computing MAX */
            d__2 = s[j];
            d__3 = (d__1 = a[j + j * a_dim1], f2c_abs(d__1)); // , expr subst
            s[j] = max(d__2,d__3);
            /* Computing MAX */
            d__2 = *amax;
            d__3 = (d__1 = a[j + j * a_dim1], f2c_abs(d__1)); // , expr subst
            *amax = max(d__2,d__3);
        }
    }
    else
    {
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Computing MAX */
            d__2 = s[j];
            d__3 = (d__1 = a[j + j * a_dim1], f2c_abs(d__1)); // , expr subst
            s[j] = max(d__2,d__3);
            /* Computing MAX */
            d__2 = *amax;
            d__3 = (d__1 = a[j + j * a_dim1], f2c_abs(d__1)); // , expr subst
            *amax = max(d__2,d__3);
            i__2 = *n;
            for (i__ = j + 1;
                    i__ <= i__2;
                    ++i__)
            {
                /* Computing MAX */
                d__2 = s[i__];
                d__3 = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)); // , expr subst
                s[i__] = max(d__2,d__3);
                /* Computing MAX */
                d__2 = s[j];
                d__3 = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)); // , expr subst
                s[j] = max(d__2,d__3);
                /* Computing MAX */
                d__2 = *amax;
                d__3 = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)); // , expr subst
                *amax = max(d__2,d__3);
            }
        }
    }
    i__1 = *n;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        s[j] = 1. / s[j];
    }
    tol = 1. / sqrt(*n * 2.);
    for (iter = 1;
            iter <= 100;
            ++iter)
    {
        scale = 0.;
        sumsq = 0.;
        /* BETA = |A|S */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            work[i__] = 0.;
        }
        if (up)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    t = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1));
                    work[i__] += (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)) * s[ j];
                    work[j] += (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)) * s[ i__];
                }
                work[j] += (d__1 = a[j + j * a_dim1], f2c_abs(d__1)) * s[j];
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                work[j] += (d__1 = a[j + j * a_dim1], f2c_abs(d__1)) * s[j];
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    t = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1));
                    work[i__] += (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)) * s[ j];
                    work[j] += (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)) * s[ i__];
                }
            }
        }
        /* avg = s^T beta / n */
        avg = 0.;
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            avg += s[i__] * work[i__];
        }
        avg /= *n;
        std = 0.;
        i__1 = *n * 3;
        for (i__ = (*n << 1) + 1;
                i__ <= i__1;
                ++i__)
        {
            work[i__] = s[i__ - (*n << 1)] * work[i__ - (*n << 1)] - avg;
        }
        dlassq_(n, &work[(*n << 1) + 1], &c__1, &scale, &sumsq);
        std = scale * sqrt(sumsq / *n);
        if (std < tol * avg)
        {
            goto L999;
        }
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            t = (d__1 = a[i__ + i__ * a_dim1], f2c_abs(d__1));
            si = s[i__];
            c2 = (*n - 1) * t;
            c1 = (*n - 2) * (work[i__] - t * si);
            c0 = -(t * si) * si + work[i__] * 2 * si - *n * avg;
            d__ = c1 * c1 - c0 * 4 * c2;
            if (d__ <= 0.)
            {
                *info = -1;
                return 0;
            }
            si = c0 * -2 / (c1 + sqrt(d__));
            d__ = si - s[i__];
            u = 0.;
            if (up)
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    t = (d__1 = a[j + i__ * a_dim1], f2c_abs(d__1));
                    u += s[j] * t;
                    work[j] += d__ * t;
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    t = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1));
                    u += s[j] * t;
                    work[j] += d__ * t;
                }
            }
            else
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    t = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1));
                    u += s[j] * t;
                    work[j] += d__ * t;
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    t = (d__1 = a[j + i__ * a_dim1], f2c_abs(d__1));
                    u += s[j] * t;
                    work[j] += d__ * t;
                }
            }
            avg += (u + work[i__]) * d__ / *n;
            s[i__] = si;
        }
    }
L999:
    smlnum = dlamch_("SAFEMIN");
    bignum = 1. / smlnum;
    smin = bignum;
    smax = 0.;
    t = 1. / sqrt(avg);
    base = dlamch_("B");
    u = 1. / log(base);
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = (integer) (u * log(s[i__] * t));
        s[i__] = pow_di(&base, &i__2);
        /* Computing MIN */
        d__1 = smin;
        d__2 = s[i__]; // , expr subst
        smin = min(d__1,d__2);
        /* Computing MAX */
        d__1 = smax;
        d__2 = s[i__]; // , expr subst
        smax = max(d__1,d__2);
    }
    *scond = max(smin,smlnum) / min(smax,bignum);
    return 0;
}
/* dsyequb_ */
