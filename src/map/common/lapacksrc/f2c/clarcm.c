/* ../netlib/clarcm.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b6 = 1.f;
static real c_b7 = 0.f;
/* > \brief \b CLARCM copies all or part of a real two-dimensional array to a complex array. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARCM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarcm. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarcm. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarcm. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK ) */
/* .. Scalar Arguments .. */
/* INTEGER LDA, LDB, LDC, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), RWORK( * ) */
/* COMPLEX B( LDB, * ), C( LDC, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARCM performs a very simple matrix-matrix multiplication: */
/* > C := A * B, */
/* > where A is M by M and real;
B is M by N and complex;
*/
/* > C is M by N and complex. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A and of the matrix C. */
/* > M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns and rows of the matrix B and */
/* > the number of columns of the matrix C. */
/* > N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, M) */
/* > A contains the M by M matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >=max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB, N) */
/* > B contains the M by N matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >=max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is COMPLEX array, dimension (LDC, N) */
/* > C contains the M by N matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >=max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (2*M*N) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int clarcm_(integer *m, integer *n, real *a, integer *lda, complex *b, integer *ldb, complex *c__, integer *ldc, real *rwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    integer i__, j, l;
    extern /* Subroutine */
    int sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
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
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --rwork;
    /* Function Body */
    if (*m == 0 || *n == 0)
    {
        return 0;
    }
    i__1 = *n;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        i__2 = *m;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            i__3 = i__ + j * b_dim1;
            rwork[(j - 1) * *m + i__] = b[i__3].r;
            /* L10: */
        }
        /* L20: */
    }
    l = *m * *n + 1;
    sgemm_("N", "N", m, n, m, &c_b6, &a[a_offset], lda, &rwork[1], m, &c_b7, & rwork[l], m);
    i__1 = *n;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        i__2 = *m;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            i__3 = i__ + j * c_dim1;
            i__4 = l + (j - 1) * *m + i__ - 1;
            c__[i__3].r = rwork[i__4];
            c__[i__3].i = 0.f; // , expr subst
            /* L30: */
        }
        /* L40: */
    }
    i__1 = *n;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        i__2 = *m;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            rwork[(j - 1) * *m + i__] = r_imag(&b[i__ + j * b_dim1]);
            /* L50: */
        }
        /* L60: */
    }
    sgemm_("N", "N", m, n, m, &c_b6, &a[a_offset], lda, &rwork[1], m, &c_b7, & rwork[l], m);
    i__1 = *n;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        i__2 = *m;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            i__3 = i__ + j * c_dim1;
            i__4 = i__ + j * c_dim1;
            r__1 = c__[i__4].r;
            i__5 = l + (j - 1) * *m + i__ - 1;
            q__1.r = r__1;
            q__1.i = rwork[i__5]; // , expr subst
            c__[i__3].r = q__1.r;
            c__[i__3].i = q__1.i; // , expr subst
            /* L70: */
        }
        /* L80: */
    }
    return 0;
    /* End of CLARCM */
}
/* clarcm_ */
