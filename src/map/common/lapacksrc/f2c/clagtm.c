/* ../netlib/clagtm.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLAGTM performs a matrix-matrix product of the form C = αAB+βC, where A is a tridiagonal matr ix, B and C are rectangular matrices, and α and β are scalars, which may be 0, 1, or -1. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAGTM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clagtm. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clagtm. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clagtm. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, */
/* B, LDB ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER LDB, LDX, N, NRHS */
/* REAL ALPHA, BETA */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX B( LDB, * ), D( * ), DL( * ), DU( * ), */
/* $ X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAGTM performs a matrix-vector product of the form */
/* > */
/* > B := alpha * A * X + beta * B */
/* > */
/* > where A is a tridiagonal matrix of order N, B and X are N by NRHS */
/* > matrices, and alpha and beta are real scalars, each of which may be */
/* > 0., 1., or -1. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the operation applied to A. */
/* > = 'N': No transpose, B := alpha * A * X + beta * B */
/* > = 'T': Transpose, B := alpha * A**T * X + beta * B */
/* > = 'C': Conjugate transpose, B := alpha * A**H * X + beta * B */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrices X and B. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is REAL */
/* > The scalar alpha. ALPHA must be 0., 1., or -1.;
otherwise, */
/* > it is assumed to be 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* > DL is COMPLEX array, dimension (N-1) */
/* > The (n-1) sub-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is COMPLEX array, dimension (N) */
/* > The diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is COMPLEX array, dimension (N-1) */
/* > The (n-1) super-diagonal elements of T. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (LDX,NRHS) */
/* > The N by NRHS matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= max(N,1). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* > BETA is REAL */
/* > The scalar beta. BETA must be 0., 1., or -1.;
otherwise, */
/* > it is assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the N by NRHS matrix B. */
/* > On exit, B is overwritten by the matrix expression */
/* > B := alpha * A * X + beta * B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(N,1). */
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
int clagtm_(char *trans, integer *n, integer *nrhs, real * alpha, complex *dl, complex *d__, complex *du, complex *x, integer * ldx, real *beta, complex *b, integer *ldb)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7, q__8, q__9;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    integer i__, j;
    extern logical lsame_(char *, char *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    if (*n == 0)
    {
        return 0;
    }
    /* Multiply B by BETA if BETA.NE.1. */
    if (*beta == 0.f)
    {
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *n;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * b_dim1;
                b[i__3].r = 0.f;
                b[i__3].i = 0.f; // , expr subst
                /* L10: */
            }
            /* L20: */
        }
    }
    else if (*beta == -1.f)
    {
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *n;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * b_dim1;
                q__1.r = -b[i__4].r;
                q__1.i = -b[i__4].i; // , expr subst
                b[i__3].r = q__1.r;
                b[i__3].i = q__1.i; // , expr subst
                /* L30: */
            }
            /* L40: */
        }
    }
    if (*alpha == 1.f)
    {
        if (lsame_(trans, "N"))
        {
            /* Compute B := B + A*X */
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    q__2.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i;
                    q__2.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4] .r; // , expr subst
                    q__1.r = b[i__3].r + q__2.r;
                    q__1.i = b[i__3].i + q__2.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    q__3.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i;
                    q__3.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4] .r; // , expr subst
                    q__2.r = b[i__3].r + q__3.r;
                    q__2.i = b[i__3].i + q__3.i; // , expr subst
                    i__5 = j * x_dim1 + 2;
                    q__4.r = du[1].r * x[i__5].r - du[1].i * x[i__5].i;
                    q__4.i = du[1].r * x[i__5].i + du[1].i * x[i__5] .r; // , expr subst
                    q__1.r = q__2.r + q__4.r;
                    q__1.i = q__2.i + q__4.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n - 1 + j * x_dim1;
                    q__3.r = dl[i__4].r * x[i__5].r - dl[i__4].i * x[i__5].i;
                    q__3.i = dl[i__4].r * x[i__5].i + dl[i__4].i * x[ i__5].r; // , expr subst
                    q__2.r = b[i__3].r + q__3.r;
                    q__2.i = b[i__3].i + q__3.i; // , expr subst
                    i__6 = *n;
                    i__7 = *n + j * x_dim1;
                    q__4.r = d__[i__6].r * x[i__7].r - d__[i__6].i * x[i__7] .i;
                    q__4.i = d__[i__6].r * x[i__7].i + d__[i__6] .i * x[i__7].r; // , expr subst
                    q__1.r = q__2.r + q__4.r;
                    q__1.i = q__2.i + q__4.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n - 1;
                    for (i__ = 2;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        i__5 = i__ - 1;
                        i__6 = i__ - 1 + j * x_dim1;
                        q__4.r = dl[i__5].r * x[i__6].r - dl[i__5].i * x[i__6] .i;
                        q__4.i = dl[i__5].r * x[i__6].i + dl[i__5] .i * x[i__6].r; // , expr subst
                        q__3.r = b[i__4].r + q__4.r;
                        q__3.i = b[i__4].i + q__4.i; // , expr subst
                        i__7 = i__;
                        i__8 = i__ + j * x_dim1;
                        q__5.r = d__[i__7].r * x[i__8].r - d__[i__7].i * x[ i__8].i;
                        q__5.i = d__[i__7].r * x[i__8].i + d__[i__7].i * x[i__8].r; // , expr subst
                        q__2.r = q__3.r + q__5.r;
                        q__2.i = q__3.i + q__5.i; // , expr subst
                        i__9 = i__;
                        i__10 = i__ + 1 + j * x_dim1;
                        q__6.r = du[i__9].r * x[i__10].r - du[i__9].i * x[ i__10].i;
                        q__6.i = du[i__9].r * x[i__10].i + du[i__9].i * x[i__10].r; // , expr subst
                        q__1.r = q__2.r + q__6.r;
                        q__1.i = q__2.i + q__6.i; // , expr subst
                        b[i__3].r = q__1.r;
                        b[i__3].i = q__1.i; // , expr subst
                        /* L50: */
                    }
                }
                /* L60: */
            }
        }
        else if (lsame_(trans, "T"))
        {
            /* Compute B := B + A**T * X */
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    q__2.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i;
                    q__2.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4] .r; // , expr subst
                    q__1.r = b[i__3].r + q__2.r;
                    q__1.i = b[i__3].i + q__2.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    q__3.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i;
                    q__3.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4] .r; // , expr subst
                    q__2.r = b[i__3].r + q__3.r;
                    q__2.i = b[i__3].i + q__3.i; // , expr subst
                    i__5 = j * x_dim1 + 2;
                    q__4.r = dl[1].r * x[i__5].r - dl[1].i * x[i__5].i;
                    q__4.i = dl[1].r * x[i__5].i + dl[1].i * x[i__5] .r; // , expr subst
                    q__1.r = q__2.r + q__4.r;
                    q__1.i = q__2.i + q__4.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n - 1 + j * x_dim1;
                    q__3.r = du[i__4].r * x[i__5].r - du[i__4].i * x[i__5].i;
                    q__3.i = du[i__4].r * x[i__5].i + du[i__4].i * x[ i__5].r; // , expr subst
                    q__2.r = b[i__3].r + q__3.r;
                    q__2.i = b[i__3].i + q__3.i; // , expr subst
                    i__6 = *n;
                    i__7 = *n + j * x_dim1;
                    q__4.r = d__[i__6].r * x[i__7].r - d__[i__6].i * x[i__7] .i;
                    q__4.i = d__[i__6].r * x[i__7].i + d__[i__6] .i * x[i__7].r; // , expr subst
                    q__1.r = q__2.r + q__4.r;
                    q__1.i = q__2.i + q__4.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n - 1;
                    for (i__ = 2;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        i__5 = i__ - 1;
                        i__6 = i__ - 1 + j * x_dim1;
                        q__4.r = du[i__5].r * x[i__6].r - du[i__5].i * x[i__6] .i;
                        q__4.i = du[i__5].r * x[i__6].i + du[i__5] .i * x[i__6].r; // , expr subst
                        q__3.r = b[i__4].r + q__4.r;
                        q__3.i = b[i__4].i + q__4.i; // , expr subst
                        i__7 = i__;
                        i__8 = i__ + j * x_dim1;
                        q__5.r = d__[i__7].r * x[i__8].r - d__[i__7].i * x[ i__8].i;
                        q__5.i = d__[i__7].r * x[i__8].i + d__[i__7].i * x[i__8].r; // , expr subst
                        q__2.r = q__3.r + q__5.r;
                        q__2.i = q__3.i + q__5.i; // , expr subst
                        i__9 = i__;
                        i__10 = i__ + 1 + j * x_dim1;
                        q__6.r = dl[i__9].r * x[i__10].r - dl[i__9].i * x[ i__10].i;
                        q__6.i = dl[i__9].r * x[i__10].i + dl[i__9].i * x[i__10].r; // , expr subst
                        q__1.r = q__2.r + q__6.r;
                        q__1.i = q__2.i + q__6.i; // , expr subst
                        b[i__3].r = q__1.r;
                        b[i__3].i = q__1.i; // , expr subst
                        /* L70: */
                    }
                }
                /* L80: */
            }
        }
        else if (lsame_(trans, "C"))
        {
            /* Compute B := B + A**H * X */
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    r_cnjg(&q__3, &d__[1]);
                    i__4 = j * x_dim1 + 1;
                    q__2.r = q__3.r * x[i__4].r - q__3.i * x[i__4].i;
                    q__2.i = q__3.r * x[i__4].i + q__3.i * x[i__4].r; // , expr subst
                    q__1.r = b[i__3].r + q__2.r;
                    q__1.i = b[i__3].i + q__2.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    r_cnjg(&q__4, &d__[1]);
                    i__4 = j * x_dim1 + 1;
                    q__3.r = q__4.r * x[i__4].r - q__4.i * x[i__4].i;
                    q__3.i = q__4.r * x[i__4].i + q__4.i * x[i__4].r; // , expr subst
                    q__2.r = b[i__3].r + q__3.r;
                    q__2.i = b[i__3].i + q__3.i; // , expr subst
                    r_cnjg(&q__6, &dl[1]);
                    i__5 = j * x_dim1 + 2;
                    q__5.r = q__6.r * x[i__5].r - q__6.i * x[i__5].i;
                    q__5.i = q__6.r * x[i__5].i + q__6.i * x[i__5].r; // , expr subst
                    q__1.r = q__2.r + q__5.r;
                    q__1.i = q__2.i + q__5.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    r_cnjg(&q__4, &du[*n - 1]);
                    i__4 = *n - 1 + j * x_dim1;
                    q__3.r = q__4.r * x[i__4].r - q__4.i * x[i__4].i;
                    q__3.i = q__4.r * x[i__4].i + q__4.i * x[i__4].r; // , expr subst
                    q__2.r = b[i__3].r + q__3.r;
                    q__2.i = b[i__3].i + q__3.i; // , expr subst
                    r_cnjg(&q__6, &d__[*n]);
                    i__5 = *n + j * x_dim1;
                    q__5.r = q__6.r * x[i__5].r - q__6.i * x[i__5].i;
                    q__5.i = q__6.r * x[i__5].i + q__6.i * x[i__5].r; // , expr subst
                    q__1.r = q__2.r + q__5.r;
                    q__1.i = q__2.i + q__5.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n - 1;
                    for (i__ = 2;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        r_cnjg(&q__5, &du[i__ - 1]);
                        i__5 = i__ - 1 + j * x_dim1;
                        q__4.r = q__5.r * x[i__5].r - q__5.i * x[i__5].i;
                        q__4.i = q__5.r * x[i__5].i + q__5.i * x[i__5] .r; // , expr subst
                        q__3.r = b[i__4].r + q__4.r;
                        q__3.i = b[i__4].i + q__4.i; // , expr subst
                        r_cnjg(&q__7, &d__[i__]);
                        i__6 = i__ + j * x_dim1;
                        q__6.r = q__7.r * x[i__6].r - q__7.i * x[i__6].i;
                        q__6.i = q__7.r * x[i__6].i + q__7.i * x[i__6] .r; // , expr subst
                        q__2.r = q__3.r + q__6.r;
                        q__2.i = q__3.i + q__6.i; // , expr subst
                        r_cnjg(&q__9, &dl[i__]);
                        i__7 = i__ + 1 + j * x_dim1;
                        q__8.r = q__9.r * x[i__7].r - q__9.i * x[i__7].i;
                        q__8.i = q__9.r * x[i__7].i + q__9.i * x[i__7] .r; // , expr subst
                        q__1.r = q__2.r + q__8.r;
                        q__1.i = q__2.i + q__8.i; // , expr subst
                        b[i__3].r = q__1.r;
                        b[i__3].i = q__1.i; // , expr subst
                        /* L90: */
                    }
                }
                /* L100: */
            }
        }
    }
    else if (*alpha == -1.f)
    {
        if (lsame_(trans, "N"))
        {
            /* Compute B := B - A*X */
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    q__2.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i;
                    q__2.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4] .r; // , expr subst
                    q__1.r = b[i__3].r - q__2.r;
                    q__1.i = b[i__3].i - q__2.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    q__3.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i;
                    q__3.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4] .r; // , expr subst
                    q__2.r = b[i__3].r - q__3.r;
                    q__2.i = b[i__3].i - q__3.i; // , expr subst
                    i__5 = j * x_dim1 + 2;
                    q__4.r = du[1].r * x[i__5].r - du[1].i * x[i__5].i;
                    q__4.i = du[1].r * x[i__5].i + du[1].i * x[i__5] .r; // , expr subst
                    q__1.r = q__2.r - q__4.r;
                    q__1.i = q__2.i - q__4.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n - 1 + j * x_dim1;
                    q__3.r = dl[i__4].r * x[i__5].r - dl[i__4].i * x[i__5].i;
                    q__3.i = dl[i__4].r * x[i__5].i + dl[i__4].i * x[ i__5].r; // , expr subst
                    q__2.r = b[i__3].r - q__3.r;
                    q__2.i = b[i__3].i - q__3.i; // , expr subst
                    i__6 = *n;
                    i__7 = *n + j * x_dim1;
                    q__4.r = d__[i__6].r * x[i__7].r - d__[i__6].i * x[i__7] .i;
                    q__4.i = d__[i__6].r * x[i__7].i + d__[i__6] .i * x[i__7].r; // , expr subst
                    q__1.r = q__2.r - q__4.r;
                    q__1.i = q__2.i - q__4.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n - 1;
                    for (i__ = 2;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        i__5 = i__ - 1;
                        i__6 = i__ - 1 + j * x_dim1;
                        q__4.r = dl[i__5].r * x[i__6].r - dl[i__5].i * x[i__6] .i;
                        q__4.i = dl[i__5].r * x[i__6].i + dl[i__5] .i * x[i__6].r; // , expr subst
                        q__3.r = b[i__4].r - q__4.r;
                        q__3.i = b[i__4].i - q__4.i; // , expr subst
                        i__7 = i__;
                        i__8 = i__ + j * x_dim1;
                        q__5.r = d__[i__7].r * x[i__8].r - d__[i__7].i * x[ i__8].i;
                        q__5.i = d__[i__7].r * x[i__8].i + d__[i__7].i * x[i__8].r; // , expr subst
                        q__2.r = q__3.r - q__5.r;
                        q__2.i = q__3.i - q__5.i; // , expr subst
                        i__9 = i__;
                        i__10 = i__ + 1 + j * x_dim1;
                        q__6.r = du[i__9].r * x[i__10].r - du[i__9].i * x[ i__10].i;
                        q__6.i = du[i__9].r * x[i__10].i + du[i__9].i * x[i__10].r; // , expr subst
                        q__1.r = q__2.r - q__6.r;
                        q__1.i = q__2.i - q__6.i; // , expr subst
                        b[i__3].r = q__1.r;
                        b[i__3].i = q__1.i; // , expr subst
                        /* L110: */
                    }
                }
                /* L120: */
            }
        }
        else if (lsame_(trans, "T"))
        {
            /* Compute B := B - A**T*X */
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    q__2.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i;
                    q__2.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4] .r; // , expr subst
                    q__1.r = b[i__3].r - q__2.r;
                    q__1.i = b[i__3].i - q__2.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    i__4 = j * x_dim1 + 1;
                    q__3.r = d__[1].r * x[i__4].r - d__[1].i * x[i__4].i;
                    q__3.i = d__[1].r * x[i__4].i + d__[1].i * x[i__4] .r; // , expr subst
                    q__2.r = b[i__3].r - q__3.r;
                    q__2.i = b[i__3].i - q__3.i; // , expr subst
                    i__5 = j * x_dim1 + 2;
                    q__4.r = dl[1].r * x[i__5].r - dl[1].i * x[i__5].i;
                    q__4.i = dl[1].r * x[i__5].i + dl[1].i * x[i__5] .r; // , expr subst
                    q__1.r = q__2.r - q__4.r;
                    q__1.i = q__2.i - q__4.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n - 1 + j * x_dim1;
                    q__3.r = du[i__4].r * x[i__5].r - du[i__4].i * x[i__5].i;
                    q__3.i = du[i__4].r * x[i__5].i + du[i__4].i * x[ i__5].r; // , expr subst
                    q__2.r = b[i__3].r - q__3.r;
                    q__2.i = b[i__3].i - q__3.i; // , expr subst
                    i__6 = *n;
                    i__7 = *n + j * x_dim1;
                    q__4.r = d__[i__6].r * x[i__7].r - d__[i__6].i * x[i__7] .i;
                    q__4.i = d__[i__6].r * x[i__7].i + d__[i__6] .i * x[i__7].r; // , expr subst
                    q__1.r = q__2.r - q__4.r;
                    q__1.i = q__2.i - q__4.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n - 1;
                    for (i__ = 2;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        i__5 = i__ - 1;
                        i__6 = i__ - 1 + j * x_dim1;
                        q__4.r = du[i__5].r * x[i__6].r - du[i__5].i * x[i__6] .i;
                        q__4.i = du[i__5].r * x[i__6].i + du[i__5] .i * x[i__6].r; // , expr subst
                        q__3.r = b[i__4].r - q__4.r;
                        q__3.i = b[i__4].i - q__4.i; // , expr subst
                        i__7 = i__;
                        i__8 = i__ + j * x_dim1;
                        q__5.r = d__[i__7].r * x[i__8].r - d__[i__7].i * x[ i__8].i;
                        q__5.i = d__[i__7].r * x[i__8].i + d__[i__7].i * x[i__8].r; // , expr subst
                        q__2.r = q__3.r - q__5.r;
                        q__2.i = q__3.i - q__5.i; // , expr subst
                        i__9 = i__;
                        i__10 = i__ + 1 + j * x_dim1;
                        q__6.r = dl[i__9].r * x[i__10].r - dl[i__9].i * x[ i__10].i;
                        q__6.i = dl[i__9].r * x[i__10].i + dl[i__9].i * x[i__10].r; // , expr subst
                        q__1.r = q__2.r - q__6.r;
                        q__1.i = q__2.i - q__6.i; // , expr subst
                        b[i__3].r = q__1.r;
                        b[i__3].i = q__1.i; // , expr subst
                        /* L130: */
                    }
                }
                /* L140: */
            }
        }
        else if (lsame_(trans, "C"))
        {
            /* Compute B := B - A**H*X */
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (*n == 1)
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    r_cnjg(&q__3, &d__[1]);
                    i__4 = j * x_dim1 + 1;
                    q__2.r = q__3.r * x[i__4].r - q__3.i * x[i__4].i;
                    q__2.i = q__3.r * x[i__4].i + q__3.i * x[i__4].r; // , expr subst
                    q__1.r = b[i__3].r - q__2.r;
                    q__1.i = b[i__3].i - q__2.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                }
                else
                {
                    i__2 = j * b_dim1 + 1;
                    i__3 = j * b_dim1 + 1;
                    r_cnjg(&q__4, &d__[1]);
                    i__4 = j * x_dim1 + 1;
                    q__3.r = q__4.r * x[i__4].r - q__4.i * x[i__4].i;
                    q__3.i = q__4.r * x[i__4].i + q__4.i * x[i__4].r; // , expr subst
                    q__2.r = b[i__3].r - q__3.r;
                    q__2.i = b[i__3].i - q__3.i; // , expr subst
                    r_cnjg(&q__6, &dl[1]);
                    i__5 = j * x_dim1 + 2;
                    q__5.r = q__6.r * x[i__5].r - q__6.i * x[i__5].i;
                    q__5.i = q__6.r * x[i__5].i + q__6.i * x[i__5].r; // , expr subst
                    q__1.r = q__2.r - q__5.r;
                    q__1.i = q__2.i - q__5.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n + j * b_dim1;
                    i__3 = *n + j * b_dim1;
                    r_cnjg(&q__4, &du[*n - 1]);
                    i__4 = *n - 1 + j * x_dim1;
                    q__3.r = q__4.r * x[i__4].r - q__4.i * x[i__4].i;
                    q__3.i = q__4.r * x[i__4].i + q__4.i * x[i__4].r; // , expr subst
                    q__2.r = b[i__3].r - q__3.r;
                    q__2.i = b[i__3].i - q__3.i; // , expr subst
                    r_cnjg(&q__6, &d__[*n]);
                    i__5 = *n + j * x_dim1;
                    q__5.r = q__6.r * x[i__5].r - q__6.i * x[i__5].i;
                    q__5.i = q__6.r * x[i__5].i + q__6.i * x[i__5].r; // , expr subst
                    q__1.r = q__2.r - q__5.r;
                    q__1.i = q__2.i - q__5.i; // , expr subst
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = *n - 1;
                    for (i__ = 2;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        r_cnjg(&q__5, &du[i__ - 1]);
                        i__5 = i__ - 1 + j * x_dim1;
                        q__4.r = q__5.r * x[i__5].r - q__5.i * x[i__5].i;
                        q__4.i = q__5.r * x[i__5].i + q__5.i * x[i__5] .r; // , expr subst
                        q__3.r = b[i__4].r - q__4.r;
                        q__3.i = b[i__4].i - q__4.i; // , expr subst
                        r_cnjg(&q__7, &d__[i__]);
                        i__6 = i__ + j * x_dim1;
                        q__6.r = q__7.r * x[i__6].r - q__7.i * x[i__6].i;
                        q__6.i = q__7.r * x[i__6].i + q__7.i * x[i__6] .r; // , expr subst
                        q__2.r = q__3.r - q__6.r;
                        q__2.i = q__3.i - q__6.i; // , expr subst
                        r_cnjg(&q__9, &dl[i__]);
                        i__7 = i__ + 1 + j * x_dim1;
                        q__8.r = q__9.r * x[i__7].r - q__9.i * x[i__7].i;
                        q__8.i = q__9.r * x[i__7].i + q__9.i * x[i__7] .r; // , expr subst
                        q__1.r = q__2.r - q__8.r;
                        q__1.i = q__2.i - q__8.i; // , expr subst
                        b[i__3].r = q__1.r;
                        b[i__3].i = q__1.i; // , expr subst
                        /* L150: */
                    }
                }
                /* L160: */
            }
        }
    }
    return 0;
    /* End of CLAGTM */
}
/* clagtm_ */
