/* ../netlib/dpteqr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;
static integer c__0 = 0;
static integer c__1 = 1;
/* > \brief \b DPTEQR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DPTEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpteqr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpteqr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpteqr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER COMPZ */
/* INTEGER INFO, LDZ, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPTEQR computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric positive definite tridiagonal matrix by first factoring the */
/* > matrix using DPTTRF, and then calling DBDSQR to compute the singular */
/* > values of the bidiagonal factor. */
/* > */
/* > This routine computes the eigenvalues of the positive definite */
/* > tridiagonal matrix to high relative accuracy. This means that if the */
/* > eigenvalues range over many orders of magnitude in size, then the */
/* > small eigenvalues and corresponding eigenvectors will be computed */
/* > more accurately than, for example, with the standard QR method. */
/* > */
/* > The eigenvectors of a full or band symmetric positive definite matrix */
/* > can also be found if DSYTRD, DSPTRD, or DSBTRD has been used to */
/* > reduce this matrix to tridiagonal form. (The reduction to tridiagonal */
/* > form, however, may preclude the possibility of obtaining high */
/* > relative accuracy in the small eigenvalues of the original matrix, if */
/* > these eigenvalues range over many orders of magnitude.) */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] COMPZ */
/* > \verbatim */
/* > COMPZ is CHARACTER*1 */
/* > = 'N': Compute eigenvalues only. */
/* > = 'V': Compute eigenvectors of original symmetric */
/* > matrix also. Array Z contains the orthogonal */
/* > matrix used to reduce the original matrix to */
/* > tridiagonal form. */
/* > = 'I': Compute eigenvectors of tridiagonal matrix also. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > On entry, the n diagonal elements of the tridiagonal */
/* > matrix. */
/* > On normal exit, D contains the eigenvalues, in descending */
/* > order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N-1) */
/* > On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* > matrix. */
/* > On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* > On entry, if COMPZ = 'V', the orthogonal matrix used in the */
/* > reduction to tridiagonal form. */
/* > On exit, if COMPZ = 'V', the orthonormal eigenvectors of the */
/* > original symmetric matrix;
*/
/* > if COMPZ = 'I', the orthonormal eigenvectors of the */
/* > tridiagonal matrix. */
/* > If INFO > 0 on exit, Z contains the eigenvectors associated */
/* > with only the stored eigenvalues. */
/* > If COMPZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1, and if */
/* > COMPZ = 'V' or 'I', LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = i, and i is: */
/* > <= N the Cholesky factorization of the matrix could */
/* > not be performed because the i-th principal minor */
/* > was not positive definite. */
/* > > N the SVD algorithm failed to converge;
*/
/* > if INFO = N+i, i off-diagonal elements of the */
/* > bidiagonal factor did not converge to zero. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doublePTcomputational */
/* ===================================================================== */
/* Subroutine */
int dpteqr_(char *compz, integer *n, doublereal *d__, doublereal *e, doublereal *z__, integer *ldz, doublereal *work, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    doublereal c__[1] /* was [1][1] */
    ;
    integer i__;
    doublereal vt[1] /* was [1][1] */
    ;
    integer nru;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *), dbdsqr_(char *, integer *, integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
    integer icompz;
    extern /* Subroutine */
    int dpttrf_(integer *, doublereal *, doublereal *, integer *);
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    /* Function Body */
    *info = 0;
    if (lsame_(compz, "N"))
    {
        icompz = 0;
    }
    else if (lsame_(compz, "V"))
    {
        icompz = 1;
    }
    else if (lsame_(compz, "I"))
    {
        icompz = 2;
    }
    else
    {
        icompz = -1;
    }
    if (icompz < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n))
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DPTEQR", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    if (*n == 1)
    {
        if (icompz > 0)
        {
            z__[z_dim1 + 1] = 1.;
        }
        return 0;
    }
    if (icompz == 2)
    {
        dlaset_("Full", n, n, &c_b7, &c_b8, &z__[z_offset], ldz);
    }
    /* Call DPTTRF to factor the matrix. */
    dpttrf_(n, &d__[1], &e[1], info);
    if (*info != 0)
    {
        return 0;
    }
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        d__[i__] = sqrt(d__[i__]);
        /* L10: */
    }
    i__1 = *n - 1;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        e[i__] *= d__[i__];
        /* L20: */
    }
    /* Call DBDSQR to compute the singular values/vectors of the */
    /* bidiagonal factor. */
    if (icompz > 0)
    {
        nru = *n;
    }
    else
    {
        nru = 0;
    }
    dbdsqr_("Lower", n, &c__0, &nru, &c__0, &d__[1], &e[1], vt, &c__1, &z__[ z_offset], ldz, c__, &c__1, &work[1], info);
    /* Square the singular values. */
    if (*info == 0)
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            d__[i__] *= d__[i__];
            /* L30: */
        }
    }
    else
    {
        *info = *n + *info;
    }
    return 0;
    /* End of DPTEQR */
}
/* dpteqr_ */
