/* ../netlib/zsymv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZSYMV computes a matrix-vector product for a complex symmetric matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZSYMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsymv.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsymv.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsymv.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INCX, INCY, LDA, N */
/* COMPLEX*16 ALPHA, BETA */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), X( * ), Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYMV performs the matrix-vector operation */
/* > */
/* > y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric matrix. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > On entry, UPLO specifies whether the upper or lower */
/* > triangular part of the array A is to be referenced as */
/* > follows: */
/* > */
/* > UPLO = 'U' or 'u' Only the upper triangular part of A */
/* > is to be referenced. */
/* > */
/* > UPLO = 'L' or 'l' Only the lower triangular part of A */
/* > is to be referenced. */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > On entry, N specifies the order of the matrix A. */
/* > N must be at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is COMPLEX*16 */
/* > On entry, ALPHA specifies the scalar alpha. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension ( LDA, N ) */
/* > Before entry, with UPLO = 'U' or 'u', the leading n by n */
/* > upper triangular part of the array A must contain the upper */
/* > triangular part of the symmetric matrix and the strictly */
/* > lower triangular part of A is not referenced. */
/* > Before entry, with UPLO = 'L' or 'l', the leading n by n */
/* > lower triangular part of the array A must contain the lower */
/* > triangular part of the symmetric matrix and the strictly */
/* > upper triangular part of A is not referenced. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > On entry, LDA specifies the first dimension of A as declared */
/* > in the calling (sub) program. LDA must be at least */
/* > max( 1, N ). */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension at least */
/* > ( 1 + ( N - 1 )*f2c_abs( INCX ) ). */
/* > Before entry, the incremented array X must contain the N- */
/* > element vector x. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > On entry, INCX specifies the increment for the elements of */
/* > X. INCX must not be zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* > BETA is COMPLEX*16 */
/* > On entry, BETA specifies the scalar beta. When BETA is */
/* > supplied as zero then Y need not be set on input. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is COMPLEX*16 array, dimension at least */
/* > ( 1 + ( N - 1 )*f2c_abs( INCY ) ). */
/* > Before entry, the incremented array Y must contain the n */
/* > element vector y. On exit, Y is overwritten by the updated */
/* > vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > On entry, INCY specifies the increment for the elements of */
/* > Y. INCY must not be zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16SYauxiliary */
/* ===================================================================== */
/* Subroutine */
int zsymv_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3, z__4;
    /* Local variables */
    integer i__, j, ix, iy, jx, jy, kx, ky, info;
    doublecomplex temp1, temp2;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;
    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L"))
    {
        info = 1;
    }
    else if (*n < 0)
    {
        info = 2;
    }
    else if (*lda < max(1,*n))
    {
        info = 5;
    }
    else if (*incx == 0)
    {
        info = 7;
    }
    else if (*incy == 0)
    {
        info = 10;
    }
    if (info != 0)
    {
        xerbla_("ZSYMV ", &info);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && beta->i == 0.))
    {
        return 0;
    }
    /* Set up the start points in X and Y. */
    if (*incx > 0)
    {
        kx = 1;
    }
    else
    {
        kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0)
    {
        ky = 1;
    }
    else
    {
        ky = 1 - (*n - 1) * *incy;
    }
    /* Start the operations. In this version the elements of A are */
    /* accessed sequentially with one pass through the triangular part */
    /* of A. */
    /* First form y := beta*y. */
    if (beta->r != 1. || beta->i != 0.)
    {
        if (*incy == 1)
        {
            if (beta->r == 0. && beta->i == 0.)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = i__;
                    y[i__2].r = 0.;
                    y[i__2].i = 0.; // , expr subst
                    /* L10: */
                }
            }
            else
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = i__;
                    i__3 = i__;
                    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i;
                    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3] .r; // , expr subst
                    y[i__2].r = z__1.r;
                    y[i__2].i = z__1.i; // , expr subst
                    /* L20: */
                }
            }
        }
        else
        {
            iy = ky;
            if (beta->r == 0. && beta->i == 0.)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = iy;
                    y[i__2].r = 0.;
                    y[i__2].i = 0.; // , expr subst
                    iy += *incy;
                    /* L30: */
                }
            }
            else
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = iy;
                    i__3 = iy;
                    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i;
                    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3] .r; // , expr subst
                    y[i__2].r = z__1.r;
                    y[i__2].i = z__1.i; // , expr subst
                    iy += *incy;
                    /* L40: */
                }
            }
        }
    }
    if (alpha->r == 0. && alpha->i == 0.)
    {
        return 0;
    }
    if (lsame_(uplo, "U"))
    {
        /* Form y when A is stored in upper triangle. */
        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i;
                z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r; // , expr subst
                temp1.r = z__1.r;
                temp1.i = z__1.i; // , expr subst
                temp2.r = 0.;
                temp2.i = 0.; // , expr subst
                i__2 = j - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i;
                    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5] .r; // , expr subst
                    z__1.r = y[i__4].r + z__2.r;
                    z__1.i = y[i__4].i + z__2.i; // , expr subst
                    y[i__3].r = z__1.r;
                    y[i__3].i = z__1.i; // , expr subst
                    i__3 = i__ + j * a_dim1;
                    i__4 = i__;
                    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i;
                    z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[ i__4].r; // , expr subst
                    z__1.r = temp2.r + z__2.r;
                    z__1.i = temp2.i + z__2.i; // , expr subst
                    temp2.r = z__1.r;
                    temp2.i = z__1.i; // , expr subst
                    /* L50: */
                }
                i__2 = j;
                i__3 = j;
                i__4 = j + j * a_dim1;
                z__3.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i;
                z__3.i = temp1.r * a[i__4].i + temp1.i * a[i__4].r; // , expr subst
                z__2.r = y[i__3].r + z__3.r;
                z__2.i = y[i__3].i + z__3.i; // , expr subst
                z__4.r = alpha->r * temp2.r - alpha->i * temp2.i;
                z__4.i = alpha->r * temp2.i + alpha->i * temp2.r; // , expr subst
                z__1.r = z__2.r + z__4.r;
                z__1.i = z__2.i + z__4.i; // , expr subst
                y[i__2].r = z__1.r;
                y[i__2].i = z__1.i; // , expr subst
                /* L60: */
            }
        }
        else
        {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = jx;
                z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i;
                z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r; // , expr subst
                temp1.r = z__1.r;
                temp1.i = z__1.i; // , expr subst
                temp2.r = 0.;
                temp2.i = 0.; // , expr subst
                ix = kx;
                iy = ky;
                i__2 = j - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = i__ + j * a_dim1;
                    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i;
                    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5] .r; // , expr subst
                    z__1.r = y[i__4].r + z__2.r;
                    z__1.i = y[i__4].i + z__2.i; // , expr subst
                    y[i__3].r = z__1.r;
                    y[i__3].i = z__1.i; // , expr subst
                    i__3 = i__ + j * a_dim1;
                    i__4 = ix;
                    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i;
                    z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[ i__4].r; // , expr subst
                    z__1.r = temp2.r + z__2.r;
                    z__1.i = temp2.i + z__2.i; // , expr subst
                    temp2.r = z__1.r;
                    temp2.i = z__1.i; // , expr subst
                    ix += *incx;
                    iy += *incy;
                    /* L70: */
                }
                i__2 = jy;
                i__3 = jy;
                i__4 = j + j * a_dim1;
                z__3.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i;
                z__3.i = temp1.r * a[i__4].i + temp1.i * a[i__4].r; // , expr subst
                z__2.r = y[i__3].r + z__3.r;
                z__2.i = y[i__3].i + z__3.i; // , expr subst
                z__4.r = alpha->r * temp2.r - alpha->i * temp2.i;
                z__4.i = alpha->r * temp2.i + alpha->i * temp2.r; // , expr subst
                z__1.r = z__2.r + z__4.r;
                z__1.i = z__2.i + z__4.i; // , expr subst
                y[i__2].r = z__1.r;
                y[i__2].i = z__1.i; // , expr subst
                jx += *incx;
                jy += *incy;
                /* L80: */
            }
        }
    }
    else
    {
        /* Form y when A is stored in lower triangle. */
        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i;
                z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r; // , expr subst
                temp1.r = z__1.r;
                temp1.i = z__1.i; // , expr subst
                temp2.r = 0.;
                temp2.i = 0.; // , expr subst
                i__2 = j;
                i__3 = j;
                i__4 = j + j * a_dim1;
                z__2.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i;
                z__2.i = temp1.r * a[i__4].i + temp1.i * a[i__4].r; // , expr subst
                z__1.r = y[i__3].r + z__2.r;
                z__1.i = y[i__3].i + z__2.i; // , expr subst
                y[i__2].r = z__1.r;
                y[i__2].i = z__1.i; // , expr subst
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i;
                    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5] .r; // , expr subst
                    z__1.r = y[i__4].r + z__2.r;
                    z__1.i = y[i__4].i + z__2.i; // , expr subst
                    y[i__3].r = z__1.r;
                    y[i__3].i = z__1.i; // , expr subst
                    i__3 = i__ + j * a_dim1;
                    i__4 = i__;
                    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i;
                    z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[ i__4].r; // , expr subst
                    z__1.r = temp2.r + z__2.r;
                    z__1.i = temp2.i + z__2.i; // , expr subst
                    temp2.r = z__1.r;
                    temp2.i = z__1.i; // , expr subst
                    /* L90: */
                }
                i__2 = j;
                i__3 = j;
                z__2.r = alpha->r * temp2.r - alpha->i * temp2.i;
                z__2.i = alpha->r * temp2.i + alpha->i * temp2.r; // , expr subst
                z__1.r = y[i__3].r + z__2.r;
                z__1.i = y[i__3].i + z__2.i; // , expr subst
                y[i__2].r = z__1.r;
                y[i__2].i = z__1.i; // , expr subst
                /* L100: */
            }
        }
        else
        {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = jx;
                z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i;
                z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r; // , expr subst
                temp1.r = z__1.r;
                temp1.i = z__1.i; // , expr subst
                temp2.r = 0.;
                temp2.i = 0.; // , expr subst
                i__2 = jy;
                i__3 = jy;
                i__4 = j + j * a_dim1;
                z__2.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i;
                z__2.i = temp1.r * a[i__4].i + temp1.i * a[i__4].r; // , expr subst
                z__1.r = y[i__3].r + z__2.r;
                z__1.i = y[i__3].i + z__2.i; // , expr subst
                y[i__2].r = z__1.r;
                y[i__2].i = z__1.i; // , expr subst
                ix = jx;
                iy = jy;
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    ix += *incx;
                    iy += *incy;
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = i__ + j * a_dim1;
                    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i;
                    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5] .r; // , expr subst
                    z__1.r = y[i__4].r + z__2.r;
                    z__1.i = y[i__4].i + z__2.i; // , expr subst
                    y[i__3].r = z__1.r;
                    y[i__3].i = z__1.i; // , expr subst
                    i__3 = i__ + j * a_dim1;
                    i__4 = ix;
                    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i;
                    z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[ i__4].r; // , expr subst
                    z__1.r = temp2.r + z__2.r;
                    z__1.i = temp2.i + z__2.i; // , expr subst
                    temp2.r = z__1.r;
                    temp2.i = z__1.i; // , expr subst
                    /* L110: */
                }
                i__2 = jy;
                i__3 = jy;
                z__2.r = alpha->r * temp2.r - alpha->i * temp2.i;
                z__2.i = alpha->r * temp2.i + alpha->i * temp2.r; // , expr subst
                z__1.r = y[i__3].r + z__2.r;
                z__1.i = y[i__3].i + z__2.i; // , expr subst
                y[i__2].r = z__1.r;
                y[i__2].i = z__1.i; // , expr subst
                jx += *incx;
                jy += *incy;
                /* L120: */
            }
        }
    }
    return 0;
    /* End of ZSYMV */
}
/* zsymv_ */
