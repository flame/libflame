/* ../netlib/sdisna.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SDISNA */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SDISNA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sdisna. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sdisna. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sdisna. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SDISNA( JOB, M, N, D, SEP, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOB */
/* INTEGER INFO, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), SEP( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SDISNA computes the reciprocal condition numbers for the eigenvectors */
/* > of a real symmetric or complex Hermitian matrix or for the left or */
/* > right singular vectors of a general m-by-n matrix. The reciprocal */
/* > condition number is the 'gap' between the corresponding eigenvalue or */
/* > singular value and the nearest other one. */
/* > */
/* > The bound on the error, measured by angle in radians, in the I-th */
/* > computed vector is given by */
/* > */
/* > SLAMCH( 'E' ) * ( ANORM / SEP( I ) ) */
/* > */
/* > where ANORM = 2-norm(A) = max( f2c_abs( D(j) ) ). SEP(I) is not allowed */
/* > to be smaller than SLAMCH( 'E' )*ANORM in order to limit the size of */
/* > the error bound. */
/* > */
/* > SDISNA may also be used to compute error bounds for eigenvectors of */
/* > the generalized symmetric definite eigenproblem. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOB */
/* > \verbatim */
/* > JOB is CHARACTER*1 */
/* > Specifies for which problem the reciprocal condition numbers */
/* > should be computed: */
/* > = 'E': the eigenvectors of a symmetric/Hermitian matrix;
*/
/* > = 'L': the left singular vectors of a general matrix;
*/
/* > = 'R': the right singular vectors of a general matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > If JOB = 'L' or 'R', the number of columns of the matrix, */
/* > in which case N >= 0. Ignored if JOB = 'E'. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (M) if JOB = 'E' */
/* > dimension (min(M,N)) if JOB = 'L' or 'R' */
/* > The eigenvalues (if JOB = 'E') or singular values (if JOB = */
/* > 'L' or 'R') of the matrix, in either increasing or decreasing */
/* > order. If singular values, they must be non-negative. */
/* > \endverbatim */
/* > */
/* > \param[out] SEP */
/* > \verbatim */
/* > SEP is REAL array, dimension (M) if JOB = 'E' */
/* > dimension (min(M,N)) if JOB = 'L' or 'R' */
/* > The reciprocal condition numbers of the vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup auxOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int sdisna_(char *job, integer *m, integer *n, real *d__, real *sep, integer *info)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;
    /* Local variables */
    integer i__, k;
    real eps;
    logical decr, left, incr, sing, eigen;
    extern logical lsame_(char *, char *);
    real anorm;
    logical right;
    real oldgap;
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real newgap, thresh;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    --sep;
    --d__;
    /* Function Body */
    *info = 0;
    eigen = lsame_(job, "E");
    left = lsame_(job, "L");
    right = lsame_(job, "R");
    sing = left || right;
    if (eigen)
    {
        k = *m;
    }
    else if (sing)
    {
        k = min(*m,*n);
    }
    if (! eigen && ! sing)
    {
        *info = -1;
    }
    else if (*m < 0)
    {
        *info = -2;
    }
    else if (k < 0)
    {
        *info = -3;
    }
    else
    {
        incr = TRUE_;
        decr = TRUE_;
        i__1 = k - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            if (incr)
            {
                incr = incr && d__[i__] <= d__[i__ + 1];
            }
            if (decr)
            {
                decr = decr && d__[i__] >= d__[i__ + 1];
            }
            /* L10: */
        }
        if (sing && k > 0)
        {
            if (incr)
            {
                incr = incr && 0.f <= d__[1];
            }
            if (decr)
            {
                decr = decr && d__[k] >= 0.f;
            }
        }
        if (! (incr || decr))
        {
            *info = -4;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SDISNA", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (k == 0)
    {
        return 0;
    }
    /* Compute reciprocal condition numbers */
    if (k == 1)
    {
        sep[1] = slamch_("O");
    }
    else
    {
        oldgap = (r__1 = d__[2] - d__[1], f2c_abs(r__1));
        sep[1] = oldgap;
        i__1 = k - 1;
        for (i__ = 2;
                i__ <= i__1;
                ++i__)
        {
            newgap = (r__1 = d__[i__ + 1] - d__[i__], f2c_abs(r__1));
            sep[i__] = min(oldgap,newgap);
            oldgap = newgap;
            /* L20: */
        }
        sep[k] = oldgap;
    }
    if (sing)
    {
        if (left && *m > *n || right && *m < *n)
        {
            if (incr)
            {
                sep[1] = min(sep[1],d__[1]);
            }
            if (decr)
            {
                /* Computing MIN */
                r__1 = sep[k];
                r__2 = d__[k]; // , expr subst
                sep[k] = min(r__1,r__2);
            }
        }
    }
    /* Ensure that reciprocal condition numbers are not less than */
    /* threshold, in order to limit the size of the error bound */
    eps = slamch_("E");
    safmin = slamch_("S");
    /* Computing MAX */
    r__2 = f2c_abs(d__[1]);
    r__3 = (r__1 = d__[k], f2c_abs(r__1)); // , expr subst
    anorm = max(r__2,r__3);
    if (anorm == 0.f)
    {
        thresh = eps;
    }
    else
    {
        /* Computing MAX */
        r__1 = eps * anorm;
        thresh = max(r__1,safmin);
    }
    i__1 = k;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        /* Computing MAX */
        r__1 = sep[i__];
        sep[i__] = max(r__1,thresh);
        /* L30: */
    }
    return 0;
    /* End of SDISNA */
}
/* sdisna_ */
