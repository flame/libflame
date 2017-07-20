/* ../netlib/zgtrfs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b18 = -1.;
static doublereal c_b19 = 1.;
static doublecomplex c_b26 =
{
    1.,0.
}
;
/* > \brief \b ZGTRFS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGTRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgtrfs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgtrfs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgtrfs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, */
/* IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* DOUBLE PRECISION BERR( * ), FERR( * ), RWORK( * ) */
/* COMPLEX*16 B( LDB, * ), D( * ), DF( * ), DL( * ), */
/* $ DLF( * ), DU( * ), DU2( * ), DUF( * ), */
/* $ WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGTRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is tridiagonal, and provides */
/* > error bounds and backward error estimates for the solution. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations: */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate transpose) */
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
/* > \param[in] DL */
/* > \verbatim */
/* > DL is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) subdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is COMPLEX*16 array, dimension (N) */
/* > The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) superdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DLF */
/* > \verbatim */
/* > DLF is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) multipliers that define the matrix L from the */
/* > LU factorization of A as computed by ZGTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] DF */
/* > \verbatim */
/* > DF is COMPLEX*16 array, dimension (N) */
/* > The n diagonal elements of the upper triangular matrix U from */
/* > the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DUF */
/* > \verbatim */
/* > DUF is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* > DU2 is COMPLEX*16 array, dimension (N-2) */
/* > The (n-2) elements of the second superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices;
for 1 <= i <= n, row i of the matrix was */
/* > interchanged with row IPIV(i). IPIV(i) will always be either */
/* > i or i+1;
IPIV(i) = i indicates a row interchange was not */
/* > required. */
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
/* > On entry, the solution matrix X, as computed by ZGTTRS. */
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
/* > The estimated forward error bound for each solution vector */
/* > X(j) (the j-th column of the solution matrix X). */
/* > If XTRUE is the true solution corresponding to X(j), FERR(j) */
/* > is an estimated upper bound for the magnitude of the largest */
/* > element in (X(j) - XTRUE) divided by the magnitude of the */
/* > largest element in X(j). The estimate is as reliable as */
/* > the estimate for RCOND, and is almost always a slight */
/* > overestimate of the true error. */
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
/* > WORK is COMPLEX*16 array, dimension (2*N) */
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
/* > \ingroup complex16GTcomputational */
/* ===================================================================== */
/* Subroutine */
int zgtrfs_(char *trans, integer *n, integer *nrhs, doublecomplex *dl, doublecomplex *d__, doublecomplex *du, doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, d__11, d__12, d__13, d__14;
    doublecomplex z__1;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    integer i__, j;
    doublereal s;
    integer nz;
    doublereal eps;
    integer kase;
    doublereal safe1, safe2;
    extern logical lsame_(char *, char *);
    integer isave[3], count;
    extern /* Subroutine */
    int zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), zlacn2_( integer *, doublecomplex *, doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    doublereal safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *), zlagtm_( char *, integer *, integer *, doublereal *, doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *, integer *, doublereal *, doublecomplex *, integer *);
    logical notran;
    char transn[1], transt[1];
    doublereal lstres;
    extern /* Subroutine */
    int zgttrs_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex * , integer *, doublecomplex *, integer *, integer *);
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    --dlf;
    --df;
    --duf;
    --du2;
    --ipiv;
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
    notran = lsame_(trans, "N");
    if (! notran && ! lsame_(trans, "T") && ! lsame_( trans, "C"))
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
        *info = -13;
    }
    else if (*ldx < max(1,*n))
    {
        *info = -15;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGTRFS", &i__1);
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
    if (notran)
    {
        *(unsigned char *)transn = 'N';
        *(unsigned char *)transt = 'C';
    }
    else
    {
        *(unsigned char *)transn = 'C';
        *(unsigned char *)transt = 'N';
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
        /* Compute residual R = B - op(A) * X, */
        /* where op(A) = A, A**T, or A**H, depending on TRANS. */
        zcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
        zlagtm_(trans, n, &c__1, &c_b18, &dl[1], &d__[1], &du[1], &x[j * x_dim1 + 1], ldx, &c_b19, &work[1], n);
        /* Compute f2c_abs(op(A))*f2c_abs(x) + f2c_abs(b) for use in the backward */
        /* error bound. */
        if (notran)
        {
            if (*n == 1)
            {
                i__2 = j * b_dim1 + 1;
                i__3 = j * x_dim1 + 1;
                rwork[1] = (d__1 = b[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&b[ j * b_dim1 + 1]), f2c_abs(d__2)) + ((d__3 = d__[1].r, f2c_abs( d__3)) + (d__4 = d_imag(&d__[1]), f2c_abs(d__4))) * (( d__5 = x[i__3].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[j * x_dim1 + 1]), f2c_abs(d__6)));
            }
            else
            {
                i__2 = j * b_dim1 + 1;
                i__3 = j * x_dim1 + 1;
                i__4 = j * x_dim1 + 2;
                rwork[1] = (d__1 = b[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&b[ j * b_dim1 + 1]), f2c_abs(d__2)) + ((d__3 = d__[1].r, f2c_abs( d__3)) + (d__4 = d_imag(&d__[1]), f2c_abs(d__4))) * (( d__5 = x[i__3].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[j * x_dim1 + 1]), f2c_abs(d__6))) + ((d__7 = du[1].r, f2c_abs( d__7)) + (d__8 = d_imag(&du[1]), f2c_abs(d__8))) * ((d__9 = x[i__4].r, f2c_abs(d__9)) + (d__10 = d_imag(&x[j * x_dim1 + 2]), f2c_abs(d__10)));
                i__2 = *n - 1;
                for (i__ = 2;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__ - 1;
                    i__5 = i__ - 1 + j * x_dim1;
                    i__6 = i__;
                    i__7 = i__ + j * x_dim1;
                    i__8 = i__;
                    i__9 = i__ + 1 + j * x_dim1;
                    rwork[i__] = (d__1 = b[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&b[i__ + j * b_dim1]), f2c_abs(d__2)) + ((d__3 = dl[i__4].r, f2c_abs(d__3)) + (d__4 = d_imag(&dl[i__ - 1]), f2c_abs(d__4))) * ((d__5 = x[i__5].r, f2c_abs(d__5) ) + (d__6 = d_imag(&x[i__ - 1 + j * x_dim1]), f2c_abs( d__6))) + ((d__7 = d__[i__6].r, f2c_abs(d__7)) + ( d__8 = d_imag(&d__[i__]), f2c_abs(d__8))) * ((d__9 = x[i__7].r, f2c_abs(d__9)) + (d__10 = d_imag(&x[i__ + j * x_dim1]), f2c_abs(d__10))) + ((d__11 = du[i__8].r, f2c_abs(d__11)) + (d__12 = d_imag(&du[i__]), f2c_abs( d__12))) * ((d__13 = x[i__9].r, f2c_abs(d__13)) + ( d__14 = d_imag(&x[i__ + 1 + j * x_dim1]), f2c_abs( d__14)));
                    /* L30: */
                }
                i__2 = *n + j * b_dim1;
                i__3 = *n - 1;
                i__4 = *n - 1 + j * x_dim1;
                i__5 = *n;
                i__6 = *n + j * x_dim1;
                rwork[*n] = (d__1 = b[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&b[ *n + j * b_dim1]), f2c_abs(d__2)) + ((d__3 = dl[i__3].r, f2c_abs(d__3)) + (d__4 = d_imag(&dl[*n - 1]), f2c_abs(d__4))) * ((d__5 = x[i__4].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[* n - 1 + j * x_dim1]), f2c_abs(d__6))) + ((d__7 = d__[i__5] .r, f2c_abs(d__7)) + (d__8 = d_imag(&d__[*n]), f2c_abs(d__8))) * ((d__9 = x[i__6].r, f2c_abs(d__9)) + (d__10 = d_imag(& x[*n + j * x_dim1]), f2c_abs(d__10)));
            }
        }
        else
        {
            if (*n == 1)
            {
                i__2 = j * b_dim1 + 1;
                i__3 = j * x_dim1 + 1;
                rwork[1] = (d__1 = b[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&b[ j * b_dim1 + 1]), f2c_abs(d__2)) + ((d__3 = d__[1].r, f2c_abs( d__3)) + (d__4 = d_imag(&d__[1]), f2c_abs(d__4))) * (( d__5 = x[i__3].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[j * x_dim1 + 1]), f2c_abs(d__6)));
            }
            else
            {
                i__2 = j * b_dim1 + 1;
                i__3 = j * x_dim1 + 1;
                i__4 = j * x_dim1 + 2;
                rwork[1] = (d__1 = b[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&b[ j * b_dim1 + 1]), f2c_abs(d__2)) + ((d__3 = d__[1].r, f2c_abs( d__3)) + (d__4 = d_imag(&d__[1]), f2c_abs(d__4))) * (( d__5 = x[i__3].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[j * x_dim1 + 1]), f2c_abs(d__6))) + ((d__7 = dl[1].r, f2c_abs( d__7)) + (d__8 = d_imag(&dl[1]), f2c_abs(d__8))) * ((d__9 = x[i__4].r, f2c_abs(d__9)) + (d__10 = d_imag(&x[j * x_dim1 + 2]), f2c_abs(d__10)));
                i__2 = *n - 1;
                for (i__ = 2;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__ - 1;
                    i__5 = i__ - 1 + j * x_dim1;
                    i__6 = i__;
                    i__7 = i__ + j * x_dim1;
                    i__8 = i__;
                    i__9 = i__ + 1 + j * x_dim1;
                    rwork[i__] = (d__1 = b[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&b[i__ + j * b_dim1]), f2c_abs(d__2)) + ((d__3 = du[i__4].r, f2c_abs(d__3)) + (d__4 = d_imag(&du[i__ - 1]), f2c_abs(d__4))) * ((d__5 = x[i__5].r, f2c_abs(d__5) ) + (d__6 = d_imag(&x[i__ - 1 + j * x_dim1]), f2c_abs( d__6))) + ((d__7 = d__[i__6].r, f2c_abs(d__7)) + ( d__8 = d_imag(&d__[i__]), f2c_abs(d__8))) * ((d__9 = x[i__7].r, f2c_abs(d__9)) + (d__10 = d_imag(&x[i__ + j * x_dim1]), f2c_abs(d__10))) + ((d__11 = dl[i__8].r, f2c_abs(d__11)) + (d__12 = d_imag(&dl[i__]), f2c_abs( d__12))) * ((d__13 = x[i__9].r, f2c_abs(d__13)) + ( d__14 = d_imag(&x[i__ + 1 + j * x_dim1]), f2c_abs( d__14)));
                    /* L40: */
                }
                i__2 = *n + j * b_dim1;
                i__3 = *n - 1;
                i__4 = *n - 1 + j * x_dim1;
                i__5 = *n;
                i__6 = *n + j * x_dim1;
                rwork[*n] = (d__1 = b[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&b[ *n + j * b_dim1]), f2c_abs(d__2)) + ((d__3 = du[i__3].r, f2c_abs(d__3)) + (d__4 = d_imag(&du[*n - 1]), f2c_abs(d__4))) * ((d__5 = x[i__4].r, f2c_abs(d__5)) + (d__6 = d_imag(&x[* n - 1 + j * x_dim1]), f2c_abs(d__6))) + ((d__7 = d__[i__5] .r, f2c_abs(d__7)) + (d__8 = d_imag(&d__[*n]), f2c_abs(d__8))) * ((d__9 = x[i__6].r, f2c_abs(d__9)) + (d__10 = d_imag(& x[*n + j * x_dim1]), f2c_abs(d__10)));
            }
        }
        /* Compute componentwise relative backward error from formula */
        /* max(i) ( f2c_abs(R(i)) / ( f2c_abs(op(A))*f2c_abs(X) + f2c_abs(B) )(i) ) */
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
            zgttrs_(trans, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[ 1], &work[1], n, info);
            zaxpy_(n, &c_b26, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
            lstres = berr[j];
            ++count;
            goto L20;
        }
        /* Bound error from formula */
        /* norm(X - XTRUE) / norm(X) .le. FERR = */
        /* norm( f2c_abs(inv(op(A)))* */
        /* ( f2c_abs(R) + NZ*EPS*( f2c_abs(op(A))*f2c_abs(X)+f2c_abs(B) ))) / norm(X) */
        /* where */
        /* norm(Z) is the magnitude of the largest component of Z */
        /* inv(op(A)) is the inverse of op(A) */
        /* f2c_abs(Z) is the componentwise absolute value of the matrix or */
        /* vector Z */
        /* NZ is the maximum number of nonzeros in any row of A, plus 1 */
        /* EPS is machine epsilon */
        /* The i-th component of f2c_abs(R)+NZ*EPS*(f2c_abs(op(A))*f2c_abs(X)+f2c_abs(B)) */
        /* is incremented by SAFE1 if the i-th component of */
        /* f2c_abs(op(A))*f2c_abs(X) + f2c_abs(B) is less than SAFE2. */
        /* Use ZLACN2 to estimate the infinity-norm of the matrix */
        /* inv(op(A)) * diag(W), */
        /* where W = f2c_abs(R) + NZ*EPS*( f2c_abs(op(A))*f2c_abs(X)+f2c_abs(B) ))) */
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
        kase = 0;
L70:
        zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
        if (kase != 0)
        {
            if (kase == 1)
            {
                /* Multiply by diag(W)*inv(op(A)**H). */
                zgttrs_(transt, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], & ipiv[1], &work[1], n, info);
                i__2 = *n;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__;
                    z__1.r = rwork[i__4] * work[i__5].r;
                    z__1.i = rwork[i__4] * work[i__5].i; // , expr subst
                    work[i__3].r = z__1.r;
                    work[i__3].i = z__1.i; // , expr subst
                    /* L80: */
                }
            }
            else
            {
                /* Multiply by inv(op(A))*diag(W). */
                i__2 = *n;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__;
                    z__1.r = rwork[i__4] * work[i__5].r;
                    z__1.i = rwork[i__4] * work[i__5].i; // , expr subst
                    work[i__3].r = z__1.r;
                    work[i__3].i = z__1.i; // , expr subst
                    /* L90: */
                }
                zgttrs_(transn, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], & ipiv[1], &work[1], n, info);
            }
            goto L70;
        }
        /* Normalize error. */
        lstres = 0.;
        i__2 = *n;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            /* Computing MAX */
            i__3 = i__ + j * x_dim1;
            d__3 = lstres;
            d__4 = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[i__ + j * x_dim1]), f2c_abs(d__2)); // , expr subst
            lstres = max(d__3,d__4);
            /* L100: */
        }
        if (lstres != 0.)
        {
            ferr[j] /= lstres;
        }
        /* L110: */
    }
    return 0;
    /* End of ZGTRFS */
}
/* zgtrfs_ */
