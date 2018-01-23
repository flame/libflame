/* ../netlib/zlapll.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLAPLL measures the linear dependence of two vectors. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAPLL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlapll. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlapll. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlapll. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAPLL( N, X, INCX, Y, INCY, SSMIN ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, INCY, N */
/* DOUBLE PRECISION SSMIN */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 X( * ), Y( * ) */
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
/* > X is COMPLEX*16 array, dimension (1+(N-1)*INCX) */
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
/* > Y is COMPLEX*16 array, dimension (1+(N-1)*INCY) */
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
/* > SSMIN is DOUBLE PRECISION */
/* > The smallest singular value of the N-by-2 matrix A = ( X Y ). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int zlapll_(integer *n, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, doublereal *ssmin)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);
    /* Local variables */
    doublecomplex c__, a11, a12, a22, tau;
    extern /* Subroutine */
    int dlas2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    extern /* Double Complex */
    VOID zdotc_f2c_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    doublereal ssmax;
    extern /* Subroutine */
    int zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), zlarfg_( integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *);
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
        *ssmin = 0.;
        return 0;
    }
    /* Compute the QR factorization of the N-by-2 matrix ( X Y ) */
    zlarfg_(n, &x[1], &x[*incx + 1], incx, &tau);
    a11.r = x[1].r;
    a11.i = x[1].i; // , expr subst
    x[1].r = 1.;
    x[1].i = 0.; // , expr subst
    d_cnjg(&z__3, &tau);
    z__2.r = -z__3.r;
    z__2.i = -z__3.i; // , expr subst
    zdotc_f2c_(&z__4, n, &x[1], incx, &y[1], incy);
    z__1.r = z__2.r * z__4.r - z__2.i * z__4.i;
    z__1.i = z__2.r * z__4.i + z__2.i * z__4.r; // , expr subst
    c__.r = z__1.r;
    c__.i = z__1.i; // , expr subst
    zaxpy_(n, &c__, &x[1], incx, &y[1], incy);
    i__1 = *n - 1;
    zlarfg_(&i__1, &y[*incy + 1], &y[(*incy << 1) + 1], incy, &tau);
    a12.r = y[1].r;
    a12.i = y[1].i; // , expr subst
    i__1 = *incy + 1;
    a22.r = y[i__1].r;
    a22.i = y[i__1].i; // , expr subst
    /* Compute the SVD of 2-by-2 Upper triangular matrix. */
    d__1 = z_abs(&a11);
    d__2 = z_abs(&a12);
    d__3 = z_abs(&a22);
    dlas2_(&d__1, &d__2, &d__3, ssmin, &ssmax);
    return 0;
    /* End of ZLAPLL */
}
/* zlapll_ */
