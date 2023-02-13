/* ../netlib/cspmv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CSPMV computes a matrix-vector product for complex vectors using a complex symmetric packed mat rix */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSPMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cspmv.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cspmv.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cspmv.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSPMV( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INCX, INCY, N */
/* COMPLEX ALPHA, BETA */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX AP( * ), X( * ), Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSPMV performs the matrix-vector operation */
/* > */
/* > y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric matrix, supplied in packed form. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > On entry, UPLO specifies whether the upper or lower */
/* > triangular part of the matrix A is supplied in the packed */
/* > array AP as follows: */
/* > */
/* > UPLO = 'U' or 'u' The upper triangular part of A is */
/* > supplied in AP. */
/* > */
/* > UPLO = 'L' or 'l' The lower triangular part of A is */
/* > supplied in AP. */
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
/* > ALPHA is COMPLEX */
/* > On entry, ALPHA specifies the scalar alpha. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* > AP is COMPLEX array, dimension at least */
/* > ( ( N*( N + 1 ) )/2 ). */
/* > Before entry, with UPLO = 'U' or 'u', the array AP must */
/* > contain the upper triangular part of the symmetric matrix */
/* > packed sequentially, column by column, so that AP( 1 ) */
/* > contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* > and a( 2, 2 ) respectively, and so on. */
/* > Before entry, with UPLO = 'L' or 'l', the array AP must */
/* > contain the lower triangular part of the symmetric matrix */
/* > packed sequentially, column by column, so that AP( 1 ) */
/* > contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* > and a( 3, 1 ) respectively, and so on. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension at least */
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
/* > BETA is COMPLEX */
/* > On entry, BETA specifies the scalar beta. When BETA is */
/* > supplied as zero then Y need not be set on input. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is COMPLEX array, dimension at least */
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
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int cspmv_(char *uplo, integer *n, complex *alpha, complex * ap, complex *x, integer *incx, complex *beta, complex *y, integer * incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3, q__4;
    /* Local variables */
    integer i__, j, k, kk, ix, iy, jx, jy, kx, ky, info;
    complex temp1, temp2;
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --y;
    --x;
    --ap;
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
    else if (*incx == 0)
    {
        info = 6;
    }
    else if (*incy == 0)
    {
        info = 9;
    }
    if (info != 0)
    {
        xerbla_("CSPMV ", &info);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0 || alpha->r == 0.f && alpha->i == 0.f && (beta->r == 1.f && beta->i == 0.f))
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
    /* Start the operations. In this version the elements of the array AP */
    /* are accessed sequentially with one pass through AP. */
    /* First form y := beta*y. */
    if (beta->r != 1.f || beta->i != 0.f)
    {
        if (*incy == 1)
        {
            if (beta->r == 0.f && beta->i == 0.f)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = i__;
                    y[i__2].r = 0.f;
                    y[i__2].i = 0.f; // , expr subst
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
                    q__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i;
                    q__1.i = beta->r * y[i__3].i + beta->i * y[i__3] .r; // , expr subst
                    y[i__2].r = q__1.r;
                    y[i__2].i = q__1.i; // , expr subst
                    /* L20: */
                }
            }
        }
        else
        {
            iy = ky;
            if (beta->r == 0.f && beta->i == 0.f)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = iy;
                    y[i__2].r = 0.f;
                    y[i__2].i = 0.f; // , expr subst
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
                    q__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i;
                    q__1.i = beta->r * y[i__3].i + beta->i * y[i__3] .r; // , expr subst
                    y[i__2].r = q__1.r;
                    y[i__2].i = q__1.i; // , expr subst
                    iy += *incy;
                    /* L40: */
                }
            }
        }
    }
    if (alpha->r == 0.f && alpha->i == 0.f)
    {
        return 0;
    }
    kk = 1;
    if (lsame_(uplo, "U"))
    {
        /* Form y when AP contains the upper triangle. */
        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                q__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i;
                q__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r; // , expr subst
                temp1.r = q__1.r;
                temp1.i = q__1.i; // , expr subst
                temp2.r = 0.f;
                temp2.i = 0.f; // , expr subst
                k = kk;
                i__2 = j - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = k;
                    q__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i;
                    q__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5] .r; // , expr subst
                    q__1.r = y[i__4].r + q__2.r;
                    q__1.i = y[i__4].i + q__2.i; // , expr subst
                    y[i__3].r = q__1.r;
                    y[i__3].i = q__1.i; // , expr subst
                    i__3 = k;
                    i__4 = i__;
                    q__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i;
                    q__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[ i__4].r; // , expr subst
                    q__1.r = temp2.r + q__2.r;
                    q__1.i = temp2.i + q__2.i; // , expr subst
                    temp2.r = q__1.r;
                    temp2.i = q__1.i; // , expr subst
                    ++k;
                    /* L50: */
                }
                i__2 = j;
                i__3 = j;
                i__4 = kk + j - 1;
                q__3.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i;
                q__3.i = temp1.r * ap[i__4].i + temp1.i * ap[i__4].r; // , expr subst
                q__2.r = y[i__3].r + q__3.r;
                q__2.i = y[i__3].i + q__3.i; // , expr subst
                q__4.r = alpha->r * temp2.r - alpha->i * temp2.i;
                q__4.i = alpha->r * temp2.i + alpha->i * temp2.r; // , expr subst
                q__1.r = q__2.r + q__4.r;
                q__1.i = q__2.i + q__4.i; // , expr subst
                y[i__2].r = q__1.r;
                y[i__2].i = q__1.i; // , expr subst
                kk += j;
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
                q__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i;
                q__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r; // , expr subst
                temp1.r = q__1.r;
                temp1.i = q__1.i; // , expr subst
                temp2.r = 0.f;
                temp2.i = 0.f; // , expr subst
                ix = kx;
                iy = ky;
                i__2 = kk + j - 2;
                for (k = kk;
                        k <= i__2;
                        ++k)
                {
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = k;
                    q__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i;
                    q__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5] .r; // , expr subst
                    q__1.r = y[i__4].r + q__2.r;
                    q__1.i = y[i__4].i + q__2.i; // , expr subst
                    y[i__3].r = q__1.r;
                    y[i__3].i = q__1.i; // , expr subst
                    i__3 = k;
                    i__4 = ix;
                    q__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i;
                    q__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[ i__4].r; // , expr subst
                    q__1.r = temp2.r + q__2.r;
                    q__1.i = temp2.i + q__2.i; // , expr subst
                    temp2.r = q__1.r;
                    temp2.i = q__1.i; // , expr subst
                    ix += *incx;
                    iy += *incy;
                    /* L70: */
                }
                i__2 = jy;
                i__3 = jy;
                i__4 = kk + j - 1;
                q__3.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i;
                q__3.i = temp1.r * ap[i__4].i + temp1.i * ap[i__4].r; // , expr subst
                q__2.r = y[i__3].r + q__3.r;
                q__2.i = y[i__3].i + q__3.i; // , expr subst
                q__4.r = alpha->r * temp2.r - alpha->i * temp2.i;
                q__4.i = alpha->r * temp2.i + alpha->i * temp2.r; // , expr subst
                q__1.r = q__2.r + q__4.r;
                q__1.i = q__2.i + q__4.i; // , expr subst
                y[i__2].r = q__1.r;
                y[i__2].i = q__1.i; // , expr subst
                jx += *incx;
                jy += *incy;
                kk += j;
                /* L80: */
            }
        }
    }
    else
    {
        /* Form y when AP contains the lower triangle. */
        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                q__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i;
                q__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r; // , expr subst
                temp1.r = q__1.r;
                temp1.i = q__1.i; // , expr subst
                temp2.r = 0.f;
                temp2.i = 0.f; // , expr subst
                i__2 = j;
                i__3 = j;
                i__4 = kk;
                q__2.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i;
                q__2.i = temp1.r * ap[i__4].i + temp1.i * ap[i__4].r; // , expr subst
                q__1.r = y[i__3].r + q__2.r;
                q__1.i = y[i__3].i + q__2.i; // , expr subst
                y[i__2].r = q__1.r;
                y[i__2].i = q__1.i; // , expr subst
                k = kk + 1;
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = k;
                    q__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i;
                    q__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5] .r; // , expr subst
                    q__1.r = y[i__4].r + q__2.r;
                    q__1.i = y[i__4].i + q__2.i; // , expr subst
                    y[i__3].r = q__1.r;
                    y[i__3].i = q__1.i; // , expr subst
                    i__3 = k;
                    i__4 = i__;
                    q__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i;
                    q__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[ i__4].r; // , expr subst
                    q__1.r = temp2.r + q__2.r;
                    q__1.i = temp2.i + q__2.i; // , expr subst
                    temp2.r = q__1.r;
                    temp2.i = q__1.i; // , expr subst
                    ++k;
                    /* L90: */
                }
                i__2 = j;
                i__3 = j;
                q__2.r = alpha->r * temp2.r - alpha->i * temp2.i;
                q__2.i = alpha->r * temp2.i + alpha->i * temp2.r; // , expr subst
                q__1.r = y[i__3].r + q__2.r;
                q__1.i = y[i__3].i + q__2.i; // , expr subst
                y[i__2].r = q__1.r;
                y[i__2].i = q__1.i; // , expr subst
                kk += *n - j + 1;
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
                q__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i;
                q__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r; // , expr subst
                temp1.r = q__1.r;
                temp1.i = q__1.i; // , expr subst
                temp2.r = 0.f;
                temp2.i = 0.f; // , expr subst
                i__2 = jy;
                i__3 = jy;
                i__4 = kk;
                q__2.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i;
                q__2.i = temp1.r * ap[i__4].i + temp1.i * ap[i__4].r; // , expr subst
                q__1.r = y[i__3].r + q__2.r;
                q__1.i = y[i__3].i + q__2.i; // , expr subst
                y[i__2].r = q__1.r;
                y[i__2].i = q__1.i; // , expr subst
                ix = jx;
                iy = jy;
                i__2 = kk + *n - j;
                for (k = kk + 1;
                        k <= i__2;
                        ++k)
                {
                    ix += *incx;
                    iy += *incy;
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = k;
                    q__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i;
                    q__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5] .r; // , expr subst
                    q__1.r = y[i__4].r + q__2.r;
                    q__1.i = y[i__4].i + q__2.i; // , expr subst
                    y[i__3].r = q__1.r;
                    y[i__3].i = q__1.i; // , expr subst
                    i__3 = k;
                    i__4 = ix;
                    q__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i;
                    q__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[ i__4].r; // , expr subst
                    q__1.r = temp2.r + q__2.r;
                    q__1.i = temp2.i + q__2.i; // , expr subst
                    temp2.r = q__1.r;
                    temp2.i = q__1.i; // , expr subst
                    /* L110: */
                }
                i__2 = jy;
                i__3 = jy;
                q__2.r = alpha->r * temp2.r - alpha->i * temp2.i;
                q__2.i = alpha->r * temp2.i + alpha->i * temp2.r; // , expr subst
                q__1.r = y[i__3].r + q__2.r;
                q__1.i = y[i__3].i + q__2.i; // , expr subst
                y[i__2].r = q__1.r;
                y[i__2].i = q__1.i; // , expr subst
                jx += *incx;
                jy += *incy;
                kk += *n - j + 1;
                /* L120: */
            }
        }
    }
    return 0;
    /* End of CSPMV */
}
/* cspmv_ */
