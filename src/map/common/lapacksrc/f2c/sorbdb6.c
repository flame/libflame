/* ../netlib/sorbdb6.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b4 = 1.f;
static real c_b5 = 0.f;
static integer c__1 = 1;
static real c_b12 = -1.f;
/* > \brief \b SORBDB6 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORBDB6 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorbdb6 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorbdb6 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorbdb6 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORBDB6( M1, M2, N, X1, INCX1, X2, INCX2, Q1, LDQ1, Q2, */
/* LDQ2, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX1, INCX2, INFO, LDQ1, LDQ2, LWORK, M1, M2, */
/* $ N */
/* .. */
/* .. Array Arguments .. */
/* REAL Q1(LDQ1,*), Q2(LDQ2,*), WORK(*), X1(*), X2(*) */
/* .. */
/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > */
/* > SORBDB6 orthogonalizes the column vector */
/* > X = [ X1 ] */
/* > [ X2 ] */
/* > with respect to the columns of */
/* > Q = [ Q1 ] . */
/* > [ Q2 ] */
/* > The columns of Q must be orthonormal. */
/* > */
/* > If the projection is zero according to Kahan's "twice is enough" */
/* > criterion, then the zero vector is returned. */
/* > */
/* >\endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M1 */
/* > \verbatim */
/* > M1 is INTEGER */
/* > The dimension of X1 and the number of rows in Q1. 0 <= M1. */
/* > \endverbatim */
/* > */
/* > \param[in] M2 */
/* > \verbatim */
/* > M2 is INTEGER */
/* > The dimension of X2 and the number of rows in Q2. 0 <= M2. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns in Q1 and Q2. 0 <= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X1 */
/* > \verbatim */
/* > X1 is REAL array, dimension (M1) */
/* > On entry, the top part of the vector to be orthogonalized. */
/* > On exit, the top part of the projected vector. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX1 */
/* > \verbatim */
/* > INCX1 is INTEGER */
/* > Increment for entries of X1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X2 */
/* > \verbatim */
/* > X2 is REAL array, dimension (M2) */
/* > On entry, the bottom part of the vector to be */
/* > orthogonalized. On exit, the bottom part of the projected */
/* > vector. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX2 */
/* > \verbatim */
/* > INCX2 is INTEGER */
/* > Increment for entries of X2. */
/* > \endverbatim */
/* > */
/* > \param[in] Q1 */
/* > \verbatim */
/* > Q1 is REAL array, dimension (LDQ1, N) */
/* > The top part of the orthonormal basis matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ1 */
/* > \verbatim */
/* > LDQ1 is INTEGER */
/* > The leading dimension of Q1. LDQ1 >= M1. */
/* > \endverbatim */
/* > */
/* > \param[in] Q2 */
/* > \verbatim */
/* > Q2 is REAL array, dimension (LDQ2, N) */
/* > The bottom part of the orthonormal basis matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ2 */
/* > \verbatim */
/* > LDQ2 is INTEGER */
/* > The leading dimension of Q2. LDQ2 >= M2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= N. */
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
/* > \date July 2012 */
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int sorbdb6_(integer *m1, integer *m2, integer *n, real *x1, integer *incx1, real *x2, integer *incx2, real *q1, integer *ldq1, real *q2, integer *ldq2, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1;
    real r__1, r__2;
    /* Local variables */
    integer i__;
    real scl1, scl2, ssq1, ssq2;
    extern /* Subroutine */
    int sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), xerbla_(char *, integer *), slassq_(integer *, real *, integer *, real *, real *);
    real normsq1, normsq2;
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* July 2012 */
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
    /* .. Intrinsic Function .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test input arguments */
    /* Parameter adjustments */
    --x1;
    --x2;
    q1_dim1 = *ldq1;
    q1_offset = 1 + q1_dim1;
    q1 -= q1_offset;
    q2_dim1 = *ldq2;
    q2_offset = 1 + q2_dim1;
    q2 -= q2_offset;
    --work;
    /* Function Body */
    *info = 0;
    if (*m1 < 0)
    {
        *info = -1;
    }
    else if (*m2 < 0)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*incx1 < 1)
    {
        *info = -5;
    }
    else if (*incx2 < 1)
    {
        *info = -7;
    }
    else if (*ldq1 < max(1,*m1))
    {
        *info = -9;
    }
    else if (*ldq2 < max(1,*m2))
    {
        *info = -11;
    }
    else if (*lwork < *n)
    {
        *info = -13;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORBDB6", &i__1);
        return 0;
    }
    /* First, project X onto the orthogonal complement of Q's column */
    /* space */
    scl1 = 0.f;
    ssq1 = 1.f;
    slassq_(m1, &x1[1], incx1, &scl1, &ssq1);
    scl2 = 0.f;
    ssq2 = 1.f;
    slassq_(m2, &x2[1], incx2, &scl2, &ssq2);
    /* Computing 2nd power */
    r__1 = scl1;
    /* Computing 2nd power */
    r__2 = scl2;
    normsq1 = r__1 * r__1 * ssq1 + r__2 * r__2 * ssq2;
    if (*m1 == 0)
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            work[i__] = 0.f;
        }
    }
    else
    {
        sgemv_("C", m1, n, &c_b4, &q1[q1_offset], ldq1, &x1[1], incx1, &c_b5, &work[1], &c__1);
    }
    sgemv_("C", m2, n, &c_b4, &q2[q2_offset], ldq2, &x2[1], incx2, &c_b4, & work[1], &c__1);
    sgemv_("N", m1, n, &c_b12, &q1[q1_offset], ldq1, &work[1], &c__1, &c_b4, & x1[1], incx1);
    sgemv_("N", m2, n, &c_b12, &q2[q2_offset], ldq2, &work[1], &c__1, &c_b4, & x2[1], incx2);
    scl1 = 0.f;
    ssq1 = 1.f;
    slassq_(m1, &x1[1], incx1, &scl1, &ssq1);
    scl2 = 0.f;
    ssq2 = 1.f;
    slassq_(m2, &x2[1], incx2, &scl2, &ssq2);
    /* Computing 2nd power */
    r__1 = scl1;
    /* Computing 2nd power */
    r__2 = scl2;
    normsq2 = r__1 * r__1 * ssq1 + r__2 * r__2 * ssq2;
    /* If projection is sufficiently large in norm, then stop. */
    /* If projection is zero, then stop. */
    /* Otherwise, project again. */
    if (normsq2 >= normsq1 * .01f)
    {
        return 0;
    }
    if (normsq2 == 0.f)
    {
        return 0;
    }
    normsq1 = normsq2;
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        work[i__] = 0.f;
    }
    if (*m1 == 0)
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            work[i__] = 0.f;
        }
    }
    else
    {
        sgemv_("C", m1, n, &c_b4, &q1[q1_offset], ldq1, &x1[1], incx1, &c_b5, &work[1], &c__1);
    }
    sgemv_("C", m2, n, &c_b4, &q2[q2_offset], ldq2, &x2[1], incx2, &c_b4, & work[1], &c__1);
    sgemv_("N", m1, n, &c_b12, &q1[q1_offset], ldq1, &work[1], &c__1, &c_b4, & x1[1], incx1);
    sgemv_("N", m2, n, &c_b12, &q2[q2_offset], ldq2, &work[1], &c__1, &c_b4, & x2[1], incx2);
    scl1 = 0.f;
    ssq1 = 1.f;
    slassq_(m1, &x1[1], incx1, &scl1, &ssq1);
    scl2 = 0.f;
    ssq2 = 1.f;
    slassq_(m1, &x1[1], incx1, &scl1, &ssq1);
    /* Computing 2nd power */
    r__1 = scl1;
    /* Computing 2nd power */
    r__2 = scl2;
    normsq2 = r__1 * r__1 * ssq1 + r__2 * r__2 * ssq2;
    /* If second projection is sufficiently large in norm, then do */
    /* nothing more. Alternatively, if it shrunk significantly, then */
    /* truncate it to zero. */
    if (normsq2 < normsq1 * .01f)
    {
        i__1 = *m1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            x1[i__] = 0.f;
        }
        i__1 = *m2;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            x2[i__] = 0.f;
        }
    }
    return 0;
    /* End of SORBDB6 */
}
/* sorbdb6_ */
