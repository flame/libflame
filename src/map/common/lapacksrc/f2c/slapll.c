/* ../netlib/slapll.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLAPLL measures the linear dependence of two vectors. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAPLL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapll. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapll. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapll. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAPLL( N, X, INCX, Y, INCY, SSMIN ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, INCY, N */
/* REAL SSMIN */
/* .. */
/* .. Array Arguments .. */
/* REAL X( * ), Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given two column vectors X and Y, let */
/* > */
/* > A = ( X Y ). */
/* > */
/* > The subroutine first computes the QR factorization of A = Q*R, */
/* > and then computes the SVD of the 2-by-2 upper triangular matrix R. */
/* > The smaller singular value of R is returned in SSMIN, which is used */
/* > as the measurement of the linear dependency of the vectors X and Y. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The length of the vectors X and Y. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is REAL array, */
/* > dimension (1+(N-1)*INCX) */
/* > On entry, X contains the N-vector X. */
/* > On exit, X is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is REAL array, */
/* > dimension (1+(N-1)*INCY) */
/* > On entry, Y contains the N-vector Y. */
/* > On exit, Y is overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > The increment between successive elements of Y. INCY > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMIN */
/* > \verbatim */
/* > SSMIN is REAL */
/* > The smallest singular value of the N-by-2 matrix A = ( X Y ). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int slapll_(integer *n, real *x, integer *incx, real *y, integer *incy, real *ssmin)
{
    /* System generated locals */
    integer i__1;
    /* Local variables */
    real c__, a11, a12, a22, tau;
    extern real sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */
    int slas2_(real *, real *, real *, real *, real *) ;
    real ssmax;
    extern /* Subroutine */
    int saxpy_(integer *, real *, real *, integer *, real *, integer *), slarfg_(integer *, real *, real *, integer *, real *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --y;
    --x;
    /* Function Body */
    if (*n <= 1)
    {
        *ssmin = 0.f;
        return 0;
    }
    /* Compute the QR factorization of the N-by-2 matrix ( X Y ) */
    slarfg_(n, &x[1], &x[*incx + 1], incx, &tau);
    a11 = x[1];
    x[1] = 1.f;
    c__ = -tau * sdot_(n, &x[1], incx, &y[1], incy);
    saxpy_(n, &c__, &x[1], incx, &y[1], incy);
    i__1 = *n - 1;
    slarfg_(&i__1, &y[*incy + 1], &y[(*incy << 1) + 1], incy, &tau);
    a12 = y[1];
    a22 = y[*incy + 1];
    /* Compute the SVD of 2-by-2 Upper triangular matrix. */
    slas2_(&a11, &a12, &a22, ssmin, &ssmax);
    return 0;
    /* End of SLAPLL */
}
/* slapll_ */
