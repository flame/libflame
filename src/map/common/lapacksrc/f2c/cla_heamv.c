/* ../netlib/cla_heamv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLA_HEAMV computes a matrix-vector product using a Hermitian indefinite matrix to calculate err or bounds. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLA_HEAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_hea mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_hea mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_hea mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLA_HEAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, */
/* INCY ) */
/* .. Scalar Arguments .. */
/* REAL ALPHA, BETA */
/* INTEGER INCX, INCY, LDA, N, UPLO */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), X( * ) */
/* REAL Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLA_SYAMV performs the matrix-vector operation */
/* > */
/* > y := alpha*f2c_abs(A)*f2c_abs(x) + beta*f2c_abs(y), */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > n by n symmetric matrix. */
/* > */
/* > This function is primarily used in calculating error bounds. */
/* > To protect against underflow during evaluation, components in */
/* > the resulting vector are perturbed away from zero by (N+1) */
/* > times the underflow threshold. To prevent unnecessarily large */
/* > errors for block-structure embedded in general matrices, */
/* > "symbolically" zero components are not perturbed. A zero */
/* > entry is considered "symbolic" if all multiplications involved */
/* > in computing that entry have at least one zero multiplicand. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is INTEGER */
/* > On entry, UPLO specifies whether the upper or lower */
/* > triangular part of the array A is to be referenced as */
/* > follows: */
/* > */
/* > UPLO = BLAS_UPPER Only the upper triangular part of A */
/* > is to be referenced. */
/* > */
/* > UPLO = BLAS_LOWER Only the lower triangular part of A */
/* > is to be referenced. */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > On entry, N specifies the number of columns of the matrix A. */
/* > N must be at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is REAL . */
/* > On entry, ALPHA specifies the scalar alpha. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array of DIMENSION ( LDA, n ). */
/* > Before entry, the leading m by n part of the array A must */
/* > contain the matrix of coefficients. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > On entry, LDA specifies the first dimension of A as declared */
/* > in the calling (sub) program. LDA must be at least */
/* > max( 1, n ). */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension */
/* > ( 1 + ( n - 1 )*f2c_abs( INCX ) ) */
/* > Before entry, the incremented array X must contain the */
/* > vector x. */
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
/* > BETA is REAL . */
/* > On entry, BETA specifies the scalar beta. When BETA is */
/* > supplied as zero then Y need not be set on input. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is REAL array, dimension */
/* > ( 1 + ( n - 1 )*f2c_abs( INCY ) ) */
/* > Before entry with BETA non-zero, the incremented array Y */
/* > must contain the vector y. On exit, Y is overwritten by the */
/* > updated vector y. */
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
/* > \ingroup complexHEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Level 2 Blas routine. */
/* > */
/* > -- Written on 22-October-1986. */
/* > Jack Dongarra, Argonne National Lab. */
/* > Jeremy Du Croz, Nag Central Office. */
/* > Sven Hammarling, Nag Central Office. */
/* > Richard Hanson, Sandia National Labs. */
/* > -- Modified for the absolute-value product, April 2006 */
/* > Jason Riedy, UC Berkeley */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int cla_heamv_(integer *uplo, integer *n, real *alpha, complex *a, integer *lda, complex *x, integer *incx, real *beta, real *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2;
    /* Builtin functions */
    double r_imag(complex *), r_sign(real *, real *);
    /* Local variables */
    integer i__, j;
    logical symb_zero__;
    integer iy, jx, kx, ky, info;
    real temp, safe1;
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilauplo_(char *);
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
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function Definitions .. */
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
    if (*uplo != ilauplo_("U") && *uplo != ilauplo_("L") )
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
        xerbla_("CHEMV ", &info);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0 || *alpha == 0.f && *beta == 1.f)
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
    /* Set SAFE1 essentially to be the underflow threshold times the */
    /* number of additions in each row. */
    safe1 = slamch_("Safe minimum");
    safe1 = (*n + 1) * safe1;
    /* Form y := alpha*f2c_abs(A)*f2c_abs(x) + beta*f2c_abs(y). */
    /* The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to */
    /* the inexact flag. Still doesn't help change the iteration order */
    /* to per-column. */
    iy = ky;
    if (*incx == 1)
    {
        if (*uplo == ilauplo_("U"))
        {
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                if (*beta == 0.f)
                {
                    symb_zero__ = TRUE_;
                    y[iy] = 0.f;
                }
                else if (y[iy] == 0.f)
                {
                    symb_zero__ = TRUE_;
                }
                else
                {
                    symb_zero__ = FALSE_;
                    y[iy] = *beta * (r__1 = y[iy], f2c_abs(r__1));
                }
                if (*alpha != 0.f)
                {
                    i__2 = i__;
                    for (j = 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = j + i__ * a_dim1;
                        temp = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag( &a[j + i__ * a_dim1]), f2c_abs(r__2));
                        i__3 = j;
                        symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[ i__3].i == 0.f || temp == 0.f);
                        i__3 = j;
                        y[iy] += *alpha * ((r__1 = x[i__3].r, f2c_abs(r__1)) + ( r__2 = r_imag(&x[j]), f2c_abs(r__2))) * temp;
                    }
                    i__2 = *n;
                    for (j = i__ + 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = i__ + j * a_dim1;
                        temp = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag( &a[i__ + j * a_dim1]), f2c_abs(r__2));
                        i__3 = j;
                        symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[ i__3].i == 0.f || temp == 0.f);
                        i__3 = j;
                        y[iy] += *alpha * ((r__1 = x[i__3].r, f2c_abs(r__1)) + ( r__2 = r_imag(&x[j]), f2c_abs(r__2))) * temp;
                    }
                }
                if (! symb_zero__)
                {
                    y[iy] += r_sign(&safe1, &y[iy]);
                }
                iy += *incy;
            }
        }
        else
        {
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                if (*beta == 0.f)
                {
                    symb_zero__ = TRUE_;
                    y[iy] = 0.f;
                }
                else if (y[iy] == 0.f)
                {
                    symb_zero__ = TRUE_;
                }
                else
                {
                    symb_zero__ = FALSE_;
                    y[iy] = *beta * (r__1 = y[iy], f2c_abs(r__1));
                }
                if (*alpha != 0.f)
                {
                    i__2 = i__;
                    for (j = 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = i__ + j * a_dim1;
                        temp = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag( &a[i__ + j * a_dim1]), f2c_abs(r__2));
                        i__3 = j;
                        symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[ i__3].i == 0.f || temp == 0.f);
                        i__3 = j;
                        y[iy] += *alpha * ((r__1 = x[i__3].r, f2c_abs(r__1)) + ( r__2 = r_imag(&x[j]), f2c_abs(r__2))) * temp;
                    }
                    i__2 = *n;
                    for (j = i__ + 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = j + i__ * a_dim1;
                        temp = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag( &a[j + i__ * a_dim1]), f2c_abs(r__2));
                        i__3 = j;
                        symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[ i__3].i == 0.f || temp == 0.f);
                        i__3 = j;
                        y[iy] += *alpha * ((r__1 = x[i__3].r, f2c_abs(r__1)) + ( r__2 = r_imag(&x[j]), f2c_abs(r__2))) * temp;
                    }
                }
                if (! symb_zero__)
                {
                    y[iy] += r_sign(&safe1, &y[iy]);
                }
                iy += *incy;
            }
        }
    }
    else
    {
        if (*uplo == ilauplo_("U"))
        {
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                if (*beta == 0.f)
                {
                    symb_zero__ = TRUE_;
                    y[iy] = 0.f;
                }
                else if (y[iy] == 0.f)
                {
                    symb_zero__ = TRUE_;
                }
                else
                {
                    symb_zero__ = FALSE_;
                    y[iy] = *beta * (r__1 = y[iy], f2c_abs(r__1));
                }
                jx = kx;
                if (*alpha != 0.f)
                {
                    i__2 = i__;
                    for (j = 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = j + i__ * a_dim1;
                        temp = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag( &a[j + i__ * a_dim1]), f2c_abs(r__2));
                        i__3 = j;
                        symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[ i__3].i == 0.f || temp == 0.f);
                        i__3 = jx;
                        y[iy] += *alpha * ((r__1 = x[i__3].r, f2c_abs(r__1)) + ( r__2 = r_imag(&x[jx]), f2c_abs(r__2))) * temp;
                        jx += *incx;
                    }
                    i__2 = *n;
                    for (j = i__ + 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = i__ + j * a_dim1;
                        temp = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag( &a[i__ + j * a_dim1]), f2c_abs(r__2));
                        i__3 = j;
                        symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[ i__3].i == 0.f || temp == 0.f);
                        i__3 = jx;
                        y[iy] += *alpha * ((r__1 = x[i__3].r, f2c_abs(r__1)) + ( r__2 = r_imag(&x[jx]), f2c_abs(r__2))) * temp;
                        jx += *incx;
                    }
                }
                if (! symb_zero__)
                {
                    y[iy] += r_sign(&safe1, &y[iy]);
                }
                iy += *incy;
            }
        }
        else
        {
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                if (*beta == 0.f)
                {
                    symb_zero__ = TRUE_;
                    y[iy] = 0.f;
                }
                else if (y[iy] == 0.f)
                {
                    symb_zero__ = TRUE_;
                }
                else
                {
                    symb_zero__ = FALSE_;
                    y[iy] = *beta * (r__1 = y[iy], f2c_abs(r__1));
                }
                jx = kx;
                if (*alpha != 0.f)
                {
                    i__2 = i__;
                    for (j = 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = i__ + j * a_dim1;
                        temp = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag( &a[i__ + j * a_dim1]), f2c_abs(r__2));
                        i__3 = j;
                        symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[ i__3].i == 0.f || temp == 0.f);
                        i__3 = jx;
                        y[iy] += *alpha * ((r__1 = x[i__3].r, f2c_abs(r__1)) + ( r__2 = r_imag(&x[jx]), f2c_abs(r__2))) * temp;
                        jx += *incx;
                    }
                    i__2 = *n;
                    for (j = i__ + 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = j + i__ * a_dim1;
                        temp = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag( &a[j + i__ * a_dim1]), f2c_abs(r__2));
                        i__3 = j;
                        symb_zero__ = symb_zero__ && (x[i__3].r == 0.f && x[ i__3].i == 0.f || temp == 0.f);
                        i__3 = jx;
                        y[iy] += *alpha * ((r__1 = x[i__3].r, f2c_abs(r__1)) + ( r__2 = r_imag(&x[jx]), f2c_abs(r__2))) * temp;
                        jx += *incx;
                    }
                }
                if (! symb_zero__)
                {
                    y[iy] += r_sign(&safe1, &y[iy]);
                }
                iy += *incy;
            }
        }
    }
    return 0;
    /* End of CLA_HEAMV */
}
/* cla_heamv__ */
