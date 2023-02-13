/* ../netlib/csyequb.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CSYEQUB */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSYEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyequb .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyequb .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyequb .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, N */
/* REAL AMAX, SCOND */
/* CHARACTER UPLO */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), WORK( * ) */
/* REAL S( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYEQUB computes row and column scalings intended to equilibrate a */
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
/* > A is COMPLEX array, dimension (LDA,N) */
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
/* > S is REAL array, dimension (N) */
/* > If INFO = 0, S contains the scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] SCOND */
/* > \verbatim */
/* > SCOND is REAL */
/* > If INFO = 0, S contains the ratio of the smallest S(i) to */
/* > the largest S(i). If SCOND >= 0.1 and AMAX is neither too */
/* > large nor too small, it is not worth scaling by S. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* > AMAX is REAL */
/* > Absolute value of largest matrix element. If AMAX is very */
/* > close to overflow or very close to underflow, the matrix */
/* > should be scaled. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (3*N) */
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
/* > \ingroup complexSYcomputational */
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
int csyequb_(char *uplo, integer *n, complex *a, integer * lda, real *s, real *scond, real *amax, complex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4;
    doublereal d__1;
    complex q__1, q__2, q__3, q__4;
    /* Builtin functions */
    double r_imag(complex *), sqrt(doublereal), log(doublereal), pow_ri(real * , integer *);
    /* Local variables */
    real d__;
    integer i__, j;
    real t, u, c0, c1, c2, si;
    logical up;
    real avg, std, tol, base;
    integer iter;
    real smin, smax, scale;
    extern logical lsame_(char *, char *);
    real sumsq;
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real bignum;
    extern /* Subroutine */
    int classq_(integer *, complex *, integer *, real *, real *);
    real smlnum;
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
    /* .. Statement Functions .. */
    /* .. */
    /* Statement Function Definitions */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
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
        xerbla_("CSYEQUB", &i__1);
        return 0;
    }
    up = lsame_(uplo, "U");
    *amax = 0.f;
    /* Quick return if possible. */
    if (*n == 0)
    {
        *scond = 1.f;
        return 0;
    }
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        s[i__] = 0.f;
    }
    *amax = 0.f;
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
                i__3 = i__ + j * a_dim1;
                r__3 = s[i__];
                r__4 = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)); // , expr subst
                s[i__] = max(r__3,r__4);
                /* Computing MAX */
                i__3 = i__ + j * a_dim1;
                r__3 = s[j];
                r__4 = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)); // , expr subst
                s[j] = max(r__3,r__4);
                /* Computing MAX */
                i__3 = i__ + j * a_dim1;
                r__3 = *amax;
                r__4 = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)); // , expr subst
                *amax = max(r__3,r__4);
            }
            /* Computing MAX */
            i__2 = j + j * a_dim1;
            r__3 = s[j];
            r__4 = (r__1 = a[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[j + j * a_dim1]), f2c_abs(r__2)); // , expr subst
            s[j] = max(r__3,r__4);
            /* Computing MAX */
            i__2 = j + j * a_dim1;
            r__3 = *amax;
            r__4 = (r__1 = a[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[j + j * a_dim1]), f2c_abs(r__2)); // , expr subst
            *amax = max(r__3,r__4);
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
            i__2 = j + j * a_dim1;
            r__3 = s[j];
            r__4 = (r__1 = a[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[j + j * a_dim1]), f2c_abs(r__2)); // , expr subst
            s[j] = max(r__3,r__4);
            /* Computing MAX */
            i__2 = j + j * a_dim1;
            r__3 = *amax;
            r__4 = (r__1 = a[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[j + j * a_dim1]), f2c_abs(r__2)); // , expr subst
            *amax = max(r__3,r__4);
            i__2 = *n;
            for (i__ = j + 1;
                    i__ <= i__2;
                    ++i__)
            {
                /* Computing MAX */
                i__3 = i__ + j * a_dim1;
                r__3 = s[i__];
                r__4 = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)); // , expr subst
                s[i__] = max(r__3,r__4);
                /* Computing MAX */
                i__3 = i__ + j * a_dim1;
                r__3 = s[j];
                r__4 = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)); // , expr subst
                s[j] = max(r__3,r__4);
                /* Computing MAX */
                i__3 = i__ + j * a_dim1;
                r__3 = *amax;
                r__4 = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2)); // , expr subst
                *amax = max(r__3,r__4);
            }
        }
    }
    i__1 = *n;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        s[j] = 1.f / s[j];
    }
    tol = 1.f / sqrt(*n * 2.f);
    for (iter = 1;
            iter <= 100;
            ++iter)
    {
        scale = 0.f;
        sumsq = 0.f;
        /* beta = |A|s */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__;
            work[i__2].r = 0.f;
            work[i__2].i = 0.f; // , expr subst
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
                    i__3 = i__ + j * a_dim1;
                    t = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2));
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    r__3 = ((r__1 = a[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[ i__ + j * a_dim1]), f2c_abs(r__2))) * s[j];
                    q__1.r = work[i__4].r + r__3;
                    q__1.i = work[i__4].i; // , expr subst
                    work[i__3].r = q__1.r;
                    work[i__3].i = q__1.i; // , expr subst
                    i__3 = j;
                    i__4 = j;
                    i__5 = i__ + j * a_dim1;
                    r__3 = ((r__1 = a[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[ i__ + j * a_dim1]), f2c_abs(r__2))) * s[i__];
                    q__1.r = work[i__4].r + r__3;
                    q__1.i = work[i__4].i; // , expr subst
                    work[i__3].r = q__1.r;
                    work[i__3].i = q__1.i; // , expr subst
                }
                i__2 = j;
                i__3 = j;
                i__4 = j + j * a_dim1;
                r__3 = ((r__1 = a[i__4].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[j + j * a_dim1]), f2c_abs(r__2))) * s[j];
                q__1.r = work[i__3].r + r__3;
                q__1.i = work[i__3].i; // , expr subst
                work[i__2].r = q__1.r;
                work[i__2].i = q__1.i; // , expr subst
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                i__3 = j;
                i__4 = j + j * a_dim1;
                r__3 = ((r__1 = a[i__4].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[j + j * a_dim1]), f2c_abs(r__2))) * s[j];
                q__1.r = work[i__3].r + r__3;
                q__1.i = work[i__3].i; // , expr subst
                work[i__2].r = q__1.r;
                work[i__2].i = q__1.i; // , expr subst
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * a_dim1;
                    t = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2));
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    r__3 = ((r__1 = a[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[ i__ + j * a_dim1]), f2c_abs(r__2))) * s[j];
                    q__1.r = work[i__4].r + r__3;
                    q__1.i = work[i__4].i; // , expr subst
                    work[i__3].r = q__1.r;
                    work[i__3].i = q__1.i; // , expr subst
                    i__3 = j;
                    i__4 = j;
                    i__5 = i__ + j * a_dim1;
                    r__3 = ((r__1 = a[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[ i__ + j * a_dim1]), f2c_abs(r__2))) * s[i__];
                    q__1.r = work[i__4].r + r__3;
                    q__1.i = work[i__4].i; // , expr subst
                    work[i__3].r = q__1.r;
                    work[i__3].i = q__1.i; // , expr subst
                }
            }
        }
        /* avg = s^T beta / n */
        avg = 0.f;
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__;
            i__3 = i__;
            q__2.r = s[i__2] * work[i__3].r;
            q__2.i = s[i__2] * work[i__3].i; // , expr subst
            q__1.r = avg + q__2.r;
            q__1.i = q__2.i; // , expr subst
            avg = q__1.r;
        }
        avg /= *n;
        std = 0.f;
        i__1 = *n << 1;
        for (i__ = *n + 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__;
            i__3 = i__ - *n;
            i__4 = i__ - *n;
            q__2.r = s[i__3] * work[i__4].r;
            q__2.i = s[i__3] * work[i__4].i; // , expr subst
            q__1.r = q__2.r - avg;
            q__1.i = q__2.i; // , expr subst
            work[i__2].r = q__1.r;
            work[i__2].i = q__1.i; // , expr subst
        }
        classq_(n, &work[*n + 1], &c__1, &scale, &sumsq);
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
            i__2 = i__ + i__ * a_dim1;
            t = (r__1 = a[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + i__ * a_dim1]), f2c_abs(r__2));
            si = s[i__];
            c2 = (*n - 1) * t;
            i__2 = *n - 2;
            i__3 = i__;
            r__1 = t * si;
            q__2.r = work[i__3].r - r__1;
            q__2.i = work[i__3].i; // , expr subst
            d__1 = (doublereal) i__2;
            q__1.r = d__1 * q__2.r;
            q__1.i = d__1 * q__2.i; // , expr subst
            c1 = q__1.r;
            r__1 = -(t * si) * si;
            i__2 = i__;
            d__1 = 2.;
            q__4.r = d__1 * work[i__2].r;
            q__4.i = d__1 * work[i__2].i; // , expr subst
            q__3.r = si * q__4.r;
            q__3.i = si * q__4.i; // , expr subst
            q__2.r = r__1 + q__3.r;
            q__2.i = q__3.i; // , expr subst
            r__2 = *n * avg;
            q__1.r = q__2.r - r__2;
            q__1.i = q__2.i; // , expr subst
            c0 = q__1.r;
            d__ = c1 * c1 - c0 * 4 * c2;
            if (d__ <= 0.f)
            {
                *info = -1;
                return 0;
            }
            si = c0 * -2 / (c1 + sqrt(d__));
            d__ = si - s[i__];
            u = 0.f;
            if (up)
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = j + i__ * a_dim1;
                    t = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[j + i__ * a_dim1]), f2c_abs(r__2));
                    u += s[j] * t;
                    i__3 = j;
                    i__4 = j;
                    r__1 = d__ * t;
                    q__1.r = work[i__4].r + r__1;
                    q__1.i = work[i__4].i; // , expr subst
                    work[i__3].r = q__1.r;
                    work[i__3].i = q__1.i; // , expr subst
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = i__ + j * a_dim1;
                    t = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2));
                    u += s[j] * t;
                    i__3 = j;
                    i__4 = j;
                    r__1 = d__ * t;
                    q__1.r = work[i__4].r + r__1;
                    q__1.i = work[i__4].i; // , expr subst
                    work[i__3].r = q__1.r;
                    work[i__3].i = q__1.i; // , expr subst
                }
            }
            else
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = i__ + j * a_dim1;
                    t = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2));
                    u += s[j] * t;
                    i__3 = j;
                    i__4 = j;
                    r__1 = d__ * t;
                    q__1.r = work[i__4].r + r__1;
                    q__1.i = work[i__4].i; // , expr subst
                    work[i__3].r = q__1.r;
                    work[i__3].i = q__1.i; // , expr subst
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = j + i__ * a_dim1;
                    t = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[j + i__ * a_dim1]), f2c_abs(r__2));
                    u += s[j] * t;
                    i__3 = j;
                    i__4 = j;
                    r__1 = d__ * t;
                    q__1.r = work[i__4].r + r__1;
                    q__1.i = work[i__4].i; // , expr subst
                    work[i__3].r = q__1.r;
                    work[i__3].i = q__1.i; // , expr subst
                }
            }
            i__2 = i__;
            q__4.r = u + work[i__2].r;
            q__4.i = work[i__2].i; // , expr subst
            q__3.r = d__ * q__4.r;
            q__3.i = d__ * q__4.i; // , expr subst
            d__1 = (doublereal) (*n);
            q__2.r = q__3.r / d__1;
            q__2.i = q__3.i / d__1; // , expr subst
            q__1.r = avg + q__2.r;
            q__1.i = q__2.i; // , expr subst
            avg = q__1.r;
            s[i__] = si;
        }
    }
L999:
    smlnum = slamch_("SAFEMIN");
    bignum = 1.f / smlnum;
    smin = bignum;
    smax = 0.f;
    t = 1.f / sqrt(avg);
    base = slamch_("B");
    u = 1.f / log(base);
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = (integer) (u * log(s[i__] * t));
        s[i__] = pow_ri(&base, &i__2);
        /* Computing MIN */
        r__1 = smin;
        r__2 = s[i__]; // , expr subst
        smin = min(r__1,r__2);
        /* Computing MAX */
        r__1 = smax;
        r__2 = s[i__]; // , expr subst
        smax = max(r__1,r__2);
    }
    *scond = max(smin,smlnum) / min(smax,bignum);
    return 0;
}
/* csyequb_ */
