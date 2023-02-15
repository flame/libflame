/* ../netlib/zptrfs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublecomplex c_b16 =
{
    1.,0.
}
;
/* > \brief \b ZPTRFS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZPTRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zptrfs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zptrfs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zptrfs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, */
/* FERR, BERR, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION BERR( * ), D( * ), DF( * ), FERR( * ), */
/* $ RWORK( * ) */
/* COMPLEX*16 B( LDB, * ), E( * ), EF( * ), WORK( * ), */
/* $ X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPTRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is Hermitian positive definite */
/* > and tridiagonal, and provides error bounds and backward error */
/* > estimates for the solution. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the superdiagonal or the subdiagonal of the */
/* > tridiagonal matrix A is stored and the form of the */
/* > factorization: */
/* > = 'U': E is the superdiagonal of A, and A = U**H*D*U;
*/
/* > = 'L': E is the subdiagonal of A, and A = L*D*L**H. */
/* > (The two forms are equivalent if A is real.) */
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
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > The n real diagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) off-diagonal elements of the tridiagonal matrix A */
/* > (see UPLO). */
/* > \endverbatim */
/* > */
/* > \param[in] DF */
/* > \verbatim */
/* > DF is DOUBLE PRECISION array, dimension (N) */
/* > The n diagonal elements of the diagonal matrix D from */
/* > the factorization computed by ZPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] EF */
/* > \verbatim */
/* > EF is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) off-diagonal elements of the unit bidiagonal */
/* > factor U or L from the factorization computed by ZPTTRF */
/* > (see UPLO). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > The right hand side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (LDX,NRHS) */
/* > On entry, the solution matrix X, as computed by ZPTTRS. */
/* > On exit, the improved solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* > FERR is DOUBLE PRECISION array, dimension (NRHS) */
/* > The forward error bound for each solution vector */
/* > X(j) (the j-th column of the solution matrix X). */
/* > If XTRUE is the true solution corresponding to X(j), FERR(j) */
/* > is an estimated upper bound for the magnitude of the largest */
/* > element in (X(j) - XTRUE) divided by the magnitude of the */
/* > largest element in X(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* > BERR is DOUBLE PRECISION array, dimension (NRHS) */
/* > The componentwise relative backward error of each solution */
/* > vector X(j) (i.e., the smallest relative change in */
/* > any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > ITMAX is the maximum number of steps of iterative refinement. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16PTcomputational */
/* ===================================================================== */
/* Subroutine */
int zptrfs_(char *uplo, integer *n, integer *nrhs, doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal * rwork, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, d__11, d__12;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *);
    /* Local variables */
    integer i__, j;
    doublereal s;
    doublecomplex bi, cx, dx, ex;
    integer ix, nz;
    doublereal eps, safe1, safe2;
    extern logical lsame_(char *, char *);
    integer count;
    logical upper;
    extern /* Subroutine */
    int zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *);
    extern integer idamax_(integer *, doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal lstres;
    extern /* Subroutine */
    int zpttrs_(char *, integer *, integer *, doublereal *, doublecomplex *, doublecomplex *, integer *, integer *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --e;
    --df;
    --ef;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;
    --rwork;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*nrhs < 0)
    {
        *info = -3;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -9;
    }
    else if (*ldx < max(1,*n))
    {
        *info = -11;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZPTRFS", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *nrhs == 0)
    {
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            ferr[j] = 0.;
            berr[j] = 0.;
            /* L10: */
        }
        return 0;
    }
    /* NZ = maximum number of nonzero elements in each row of A, plus 1 */
    nz = 4;
    eps = dlamch_("Epsilon");
    safmin = dlamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;
    /* Do for each right hand side */
    i__1 = *nrhs;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        count = 1;
        lstres = 3.;
L20: /* Loop until stopping criterion is satisfied. */
        /* Compute residual R = B - A * X. Also compute */
        /* f2c_abs(A)*f2c_abs(x) + f2c_abs(b) for use in the backward error bound. */
        if (upper)
        {
            if (*n == 1)
            {
                i__2 = j * b_dim1 + 1;
                bi.r = b[i__2].r;
                bi.i = b[i__2].i; // , expr subst
                i__2 = j * x_dim1 + 1;
                z__1.r = d__[1] * x[i__2].r;
                z__1.i = d__[1] * x[i__2].i; // , expr subst
                dx.r = z__1.r;
                dx.i = z__1.i; // , expr subst
                z__1.r = bi.r - dx.r;
                z__1.i = bi.i - dx.i; // , expr subst
                work[1].r = z__1.r;
                work[1].i = z__1.i; // , expr subst
                rwork[1] = (d__1 = bi.r, f2c_abs(d__1)) + (d__2 = d_imag(&bi), f2c_abs(d__2)) + ((d__3 = dx.r, f2c_abs(d__3)) + (d__4 = d_imag(&dx), f2c_abs(d__4)));
            }
            else
            {
                i__2 = j * b_dim1 + 1;
                bi.r = b[i__2].r;
                bi.i = b[i__2].i; // , expr subst
                i__2 = j * x_dim1 + 1;
                z__1.r = d__[1] * x[i__2].r;
                z__1.i = d__[1] * x[i__2].i; // , expr subst
                dx.r = z__1.r;
                dx.i = z__1.i; // , expr subst
                i__2 = j * x_dim1 + 2;
                z__1.r = e[1].r * x[i__2].r - e[1].i * x[i__2].i;
                z__1.i = e[ 1].r * x[i__2].i + e[1].i * x[i__2].r; // , expr subst
                ex.r = z__1.r;
                ex.i = z__1.i; // , expr subst
                z__2.r = bi.r - dx.r;
                z__2.i = bi.i - dx.i; // , expr subst
                z__1.r = z__2.r - ex.r;
                z__1.i = z__2.i - ex.i; // , expr subst
                work[1].r = z__1.r;
                work[1].i = z__1.i; // , expr subst
                i__2 = j * x_dim1 + 2;
                rwork[1] = (d__1 = bi.r, f2c_abs(d__1)) + (d__2 = d_imag(&bi), f2c_abs(d__2)) + ((d__3 = dx.r, f2c_abs(d__3)) + (d__4 = d_imag(&dx), f2c_abs(d__4))) + ((d__5 = e[1].r, f2c_abs(d__5)) + (d__6 = d_imag(&e[1]), f2c_abs(d__6))) * ((d__7 = x[ i__2].r, f2c_abs(d__7)) + (d__8 = d_imag(&x[j * x_dim1 + 2]), f2c_abs(d__8)));
                i__2 = *n - 1;
                for (i__ = 2;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    bi.r = b[i__3].r;
                    bi.i = b[i__3].i; // , expr subst
                    d_cnjg(&z__2, &e[i__ - 1]);
                    i__3 = i__ - 1 + j * x_dim1;
                    z__1.r = z__2.r * x[i__3].r - z__2.i * x[i__3].i;
                    z__1.i = z__2.r * x[i__3].i + z__2.i * x[i__3].r; // , expr subst
                    cx.r = z__1.r;
                    cx.i = z__1.i; // , expr subst
                    i__3 = i__;
                    i__4 = i__ + j * x_dim1;
                    z__1.r = d__[i__3] * x[i__4].r;
                    z__1.i = d__[i__3] * x[ i__4].i; // , expr subst
                    dx.r = z__1.r;
                    dx.i = z__1.i; // , expr subst
                    i__3 = i__;
                    i__4 = i__ + 1 + j * x_dim1;
                    z__1.r = e[i__3].r * x[i__4].r - e[i__3].i * x[i__4].i;
                    z__1.i = e[i__3].r * x[i__4].i + e[i__3].i * x[ i__4].r; // , expr subst
                    ex.r = z__1.r;
                    ex.i = z__1.i; // , expr subst
                    i__3 = i__;
                    z__3.r = bi.r - cx.r;
                    z__3.i = bi.i - cx.i; // , expr subst
                    z__2.r = z__3.r - dx.r;
                    z__2.i = z__3.i - dx.i; // , expr subst
                    z__1.r = z__2.r - ex.r;
                    z__1.i = z__2.i - ex.i; // , expr subst
                    work[i__3].r = z__1.r;
                    work[i__3].i = z__1.i; // , expr subst
                    i__3 = i__ - 1;
                    i__4 = i__ - 1 + j * x_dim1;
                    i__5 = i__;
                    i__6 = i__ + 1 + j * x_dim1;
                    rwork[i__] = (d__1 = bi.r, f2c_abs(d__1)) + (d__2 = d_imag(& bi), f2c_abs(d__2)) + ((d__3 = e[i__3].r, f2c_abs(d__3)) + (d__4 = d_imag(&e[i__ - 1]), f2c_abs(d__4))) * (( d__5 = x[i__4].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[ i__ - 1 + j * x_dim1]), f2c_abs(d__6))) + ((d__7 = dx.r, f2c_abs(d__7)) + (d__8 = d_imag(&dx), f2c_abs(d__8)) ) + ((d__9 = e[i__5].r, f2c_abs(d__9)) + (d__10 = d_imag(&e[i__]), f2c_abs(d__10))) * ((d__11 = x[i__6] .r, f2c_abs(d__11)) + (d__12 = d_imag(&x[i__ + 1 + j * x_dim1]), f2c_abs(d__12)));
                    /* L30: */
                }
                i__2 = *n + j * b_dim1;
                bi.r = b[i__2].r;
                bi.i = b[i__2].i; // , expr subst
                d_cnjg(&z__2, &e[*n - 1]);
                i__2 = *n - 1 + j * x_dim1;
                z__1.r = z__2.r * x[i__2].r - z__2.i * x[i__2].i;
                z__1.i = z__2.r * x[i__2].i + z__2.i * x[i__2].r; // , expr subst
                cx.r = z__1.r;
                cx.i = z__1.i; // , expr subst
                i__2 = *n;
                i__3 = *n + j * x_dim1;
                z__1.r = d__[i__2] * x[i__3].r;
                z__1.i = d__[i__2] * x[i__3] .i; // , expr subst
                dx.r = z__1.r;
                dx.i = z__1.i; // , expr subst
                i__2 = *n;
                z__2.r = bi.r - cx.r;
                z__2.i = bi.i - cx.i; // , expr subst
                z__1.r = z__2.r - dx.r;
                z__1.i = z__2.i - dx.i; // , expr subst
                work[i__2].r = z__1.r;
                work[i__2].i = z__1.i; // , expr subst
                i__2 = *n - 1;
                i__3 = *n - 1 + j * x_dim1;
                rwork[*n] = (d__1 = bi.r, f2c_abs(d__1)) + (d__2 = d_imag(&bi), f2c_abs(d__2)) + ((d__3 = e[i__2].r, f2c_abs(d__3)) + (d__4 = d_imag(&e[*n - 1]), f2c_abs(d__4))) * ((d__5 = x[i__3].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[*n - 1 + j * x_dim1]), f2c_abs(d__6))) + ((d__7 = dx.r, f2c_abs(d__7)) + (d__8 = d_imag(&dx), f2c_abs(d__8)));
            }
        }
        else
        {
            if (*n == 1)
            {
                i__2 = j * b_dim1 + 1;
                bi.r = b[i__2].r;
                bi.i = b[i__2].i; // , expr subst
                i__2 = j * x_dim1 + 1;
                z__1.r = d__[1] * x[i__2].r;
                z__1.i = d__[1] * x[i__2].i; // , expr subst
                dx.r = z__1.r;
                dx.i = z__1.i; // , expr subst
                z__1.r = bi.r - dx.r;
                z__1.i = bi.i - dx.i; // , expr subst
                work[1].r = z__1.r;
                work[1].i = z__1.i; // , expr subst
                rwork[1] = (d__1 = bi.r, f2c_abs(d__1)) + (d__2 = d_imag(&bi), f2c_abs(d__2)) + ((d__3 = dx.r, f2c_abs(d__3)) + (d__4 = d_imag(&dx), f2c_abs(d__4)));
            }
            else
            {
                i__2 = j * b_dim1 + 1;
                bi.r = b[i__2].r;
                bi.i = b[i__2].i; // , expr subst
                i__2 = j * x_dim1 + 1;
                z__1.r = d__[1] * x[i__2].r;
                z__1.i = d__[1] * x[i__2].i; // , expr subst
                dx.r = z__1.r;
                dx.i = z__1.i; // , expr subst
                d_cnjg(&z__2, &e[1]);
                i__2 = j * x_dim1 + 2;
                z__1.r = z__2.r * x[i__2].r - z__2.i * x[i__2].i;
                z__1.i = z__2.r * x[i__2].i + z__2.i * x[i__2].r; // , expr subst
                ex.r = z__1.r;
                ex.i = z__1.i; // , expr subst
                z__2.r = bi.r - dx.r;
                z__2.i = bi.i - dx.i; // , expr subst
                z__1.r = z__2.r - ex.r;
                z__1.i = z__2.i - ex.i; // , expr subst
                work[1].r = z__1.r;
                work[1].i = z__1.i; // , expr subst
                i__2 = j * x_dim1 + 2;
                rwork[1] = (d__1 = bi.r, f2c_abs(d__1)) + (d__2 = d_imag(&bi), f2c_abs(d__2)) + ((d__3 = dx.r, f2c_abs(d__3)) + (d__4 = d_imag(&dx), f2c_abs(d__4))) + ((d__5 = e[1].r, f2c_abs(d__5)) + (d__6 = d_imag(&e[1]), f2c_abs(d__6))) * ((d__7 = x[ i__2].r, f2c_abs(d__7)) + (d__8 = d_imag(&x[j * x_dim1 + 2]), f2c_abs(d__8)));
                i__2 = *n - 1;
                for (i__ = 2;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    bi.r = b[i__3].r;
                    bi.i = b[i__3].i; // , expr subst
                    i__3 = i__ - 1;
                    i__4 = i__ - 1 + j * x_dim1;
                    z__1.r = e[i__3].r * x[i__4].r - e[i__3].i * x[i__4].i;
                    z__1.i = e[i__3].r * x[i__4].i + e[i__3].i * x[ i__4].r; // , expr subst
                    cx.r = z__1.r;
                    cx.i = z__1.i; // , expr subst
                    i__3 = i__;
                    i__4 = i__ + j * x_dim1;
                    z__1.r = d__[i__3] * x[i__4].r;
                    z__1.i = d__[i__3] * x[ i__4].i; // , expr subst
                    dx.r = z__1.r;
                    dx.i = z__1.i; // , expr subst
                    d_cnjg(&z__2, &e[i__]);
                    i__3 = i__ + 1 + j * x_dim1;
                    z__1.r = z__2.r * x[i__3].r - z__2.i * x[i__3].i;
                    z__1.i = z__2.r * x[i__3].i + z__2.i * x[i__3].r; // , expr subst
                    ex.r = z__1.r;
                    ex.i = z__1.i; // , expr subst
                    i__3 = i__;
                    z__3.r = bi.r - cx.r;
                    z__3.i = bi.i - cx.i; // , expr subst
                    z__2.r = z__3.r - dx.r;
                    z__2.i = z__3.i - dx.i; // , expr subst
                    z__1.r = z__2.r - ex.r;
                    z__1.i = z__2.i - ex.i; // , expr subst
                    work[i__3].r = z__1.r;
                    work[i__3].i = z__1.i; // , expr subst
                    i__3 = i__ - 1;
                    i__4 = i__ - 1 + j * x_dim1;
                    i__5 = i__;
                    i__6 = i__ + 1 + j * x_dim1;
                    rwork[i__] = (d__1 = bi.r, f2c_abs(d__1)) + (d__2 = d_imag(& bi), f2c_abs(d__2)) + ((d__3 = e[i__3].r, f2c_abs(d__3)) + (d__4 = d_imag(&e[i__ - 1]), f2c_abs(d__4))) * (( d__5 = x[i__4].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[ i__ - 1 + j * x_dim1]), f2c_abs(d__6))) + ((d__7 = dx.r, f2c_abs(d__7)) + (d__8 = d_imag(&dx), f2c_abs(d__8)) ) + ((d__9 = e[i__5].r, f2c_abs(d__9)) + (d__10 = d_imag(&e[i__]), f2c_abs(d__10))) * ((d__11 = x[i__6] .r, f2c_abs(d__11)) + (d__12 = d_imag(&x[i__ + 1 + j * x_dim1]), f2c_abs(d__12)));
                    /* L40: */
                }
                i__2 = *n + j * b_dim1;
                bi.r = b[i__2].r;
                bi.i = b[i__2].i; // , expr subst
                i__2 = *n - 1;
                i__3 = *n - 1 + j * x_dim1;
                z__1.r = e[i__2].r * x[i__3].r - e[i__2].i * x[i__3].i;
                z__1.i = e[i__2].r * x[i__3].i + e[i__2].i * x[i__3] .r; // , expr subst
                cx.r = z__1.r;
                cx.i = z__1.i; // , expr subst
                i__2 = *n;
                i__3 = *n + j * x_dim1;
                z__1.r = d__[i__2] * x[i__3].r;
                z__1.i = d__[i__2] * x[i__3] .i; // , expr subst
                dx.r = z__1.r;
                dx.i = z__1.i; // , expr subst
                i__2 = *n;
                z__2.r = bi.r - cx.r;
                z__2.i = bi.i - cx.i; // , expr subst
                z__1.r = z__2.r - dx.r;
                z__1.i = z__2.i - dx.i; // , expr subst
                work[i__2].r = z__1.r;
                work[i__2].i = z__1.i; // , expr subst
                i__2 = *n - 1;
                i__3 = *n - 1 + j * x_dim1;
                rwork[*n] = (d__1 = bi.r, f2c_abs(d__1)) + (d__2 = d_imag(&bi), f2c_abs(d__2)) + ((d__3 = e[i__2].r, f2c_abs(d__3)) + (d__4 = d_imag(&e[*n - 1]), f2c_abs(d__4))) * ((d__5 = x[i__3].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[*n - 1 + j * x_dim1]), f2c_abs(d__6))) + ((d__7 = dx.r, f2c_abs(d__7)) + (d__8 = d_imag(&dx), f2c_abs(d__8)));
            }
        }
        /* Compute componentwise relative backward error from formula */
        /* max(i) ( f2c_abs(R(i)) / ( f2c_abs(A)*f2c_abs(X) + f2c_abs(B) )(i) ) */
        /* where f2c_abs(Z) is the componentwise absolute value of the matrix */
        /* or vector Z. If the i-th component of the denominator is less */
        /* than SAFE2, then SAFE1 is added to the i-th components of the */
        /* numerator and denominator before dividing. */
        s = 0.;
        i__2 = *n;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            if (rwork[i__] > safe2)
            {
                /* Computing MAX */
                i__3 = i__;
                d__3 = s;
                d__4 = ((d__1 = work[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&work[i__]), f2c_abs(d__2))) / rwork[i__]; // , expr subst
                s = max(d__3,d__4);
            }
            else
            {
                /* Computing MAX */
                i__3 = i__;
                d__3 = s;
                d__4 = ((d__1 = work[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&work[i__]), f2c_abs(d__2)) + safe1) / (rwork[i__] + safe1); // , expr subst
                s = max(d__3,d__4);
            }
            /* L50: */
        }
        berr[j] = s;
        /* Test stopping criterion. Continue iterating if */
        /* 1) The residual BERR(J) is larger than machine epsilon, and */
        /* 2) BERR(J) decreased by at least a factor of 2 during the */
        /* last iteration, and */
        /* 3) At most ITMAX iterations tried. */
        if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5)
        {
            /* Update solution and try again. */
            zpttrs_(uplo, n, &c__1, &df[1], &ef[1], &work[1], n, info);
            zaxpy_(n, &c_b16, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
            lstres = berr[j];
            ++count;
            goto L20;
        }
        /* Bound error from formula */
        /* norm(X - XTRUE) / norm(X) .le. FERR = */
        /* norm( f2c_abs(inv(A))* */
        /* ( f2c_abs(R) + NZ*EPS*( f2c_abs(A)*f2c_abs(X)+f2c_abs(B) ))) / norm(X) */
        /* where */
        /* norm(Z) is the magnitude of the largest component of Z */
        /* inv(A) is the inverse of A */
        /* f2c_abs(Z) is the componentwise absolute value of the matrix or */
        /* vector Z */
        /* NZ is the maximum number of nonzeros in any row of A, plus 1 */
        /* EPS is machine epsilon */
        /* The i-th component of f2c_abs(R)+NZ*EPS*(f2c_abs(A)*f2c_abs(X)+f2c_abs(B)) */
        /* is incremented by SAFE1 if the i-th component of */
        /* f2c_abs(A)*f2c_abs(X) + f2c_abs(B) is less than SAFE2. */
        i__2 = *n;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            if (rwork[i__] > safe2)
            {
                i__3 = i__;
                rwork[i__] = (d__1 = work[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&work[i__]), f2c_abs(d__2)) + nz * eps * rwork[i__] ;
            }
            else
            {
                i__3 = i__;
                rwork[i__] = (d__1 = work[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&work[i__]), f2c_abs(d__2)) + nz * eps * rwork[i__] + safe1;
            }
            /* L60: */
        }
        ix = idamax_(n, &rwork[1], &c__1);
        ferr[j] = rwork[ix];
        /* Estimate the norm of inv(A). */
        /* Solve M(A) * x = e, where M(A) = (m(i,j)) is given by */
        /* m(i,j) = f2c_abs(A(i,j)); i = j; */
        /* m(i,j) = -f2c_abs(A(i,j)), i .ne. j, */
        /* and e = [ 1, 1, ..., 1 ]**T. Note M(A) = M(L)*D*M(L)**H. */
        /* Solve M(L) * x = e. */
        rwork[1] = 1.;
        i__2 = *n;
        for (i__ = 2;
                i__ <= i__2;
                ++i__)
        {
            rwork[i__] = rwork[i__ - 1] * z_abs(&ef[i__ - 1]) + 1.;
            /* L70: */
        }
        /* Solve D * M(L)**H * x = b. */
        rwork[*n] /= df[*n];
        for (i__ = *n - 1;
                i__ >= 1;
                --i__)
        {
            rwork[i__] = rwork[i__] / df[i__] + rwork[i__ + 1] * z_abs(&ef[ i__]);
            /* L80: */
        }
        /* Compute norm(inv(A)) = max(x(i)), 1<=i<=n. */
        ix = idamax_(n, &rwork[1], &c__1);
        ferr[j] *= (d__1 = rwork[ix], f2c_abs(d__1));
        /* Normalize error. */
        lstres = 0.;
        i__2 = *n;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            /* Computing MAX */
            d__1 = lstres;
            d__2 = z_abs(&x[i__ + j * x_dim1]); // , expr subst
            lstres = max(d__1,d__2);
            /* L90: */
        }
        if (lstres != 0.)
        {
            ferr[j] /= lstres;
        }
        /* L100: */
    }
    return 0;
    /* End of ZPTRFS */
}
/* zptrfs_ */
