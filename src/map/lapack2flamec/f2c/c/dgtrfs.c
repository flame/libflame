/* ../netlib/dgtrfs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b18 = -1.;
static doublereal c_b19 = 1.;
/* > \brief \b DGTRFS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGTRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtrfs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtrfs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtrfs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, */
/* IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), IWORK( * ) */
/* DOUBLE PRECISION B( LDB, * ), BERR( * ), D( * ), DF( * ), */
/* $ DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), */
/* $ FERR( * ), WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGTRFS improves the computed solution to a system of linear */
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
/* > = 'C': A**H * X = B (Conjugate transpose = Transpose) */
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
/* > DL is DOUBLE PRECISION array, dimension (N-1) */
/* > The (n-1) subdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is DOUBLE PRECISION array, dimension (N-1) */
/* > The (n-1) superdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DLF */
/* > \verbatim */
/* > DLF is DOUBLE PRECISION array, dimension (N-1) */
/* > The (n-1) multipliers that define the matrix L from the */
/* > LU factorization of A as computed by DGTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] DF */
/* > \verbatim */
/* > DF is DOUBLE PRECISION array, dimension (N) */
/* > The n diagonal elements of the upper triangular matrix U from */
/* > the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DUF */
/* > \verbatim */
/* > DUF is DOUBLE PRECISION array, dimension (N-1) */
/* > The (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* > DU2 is DOUBLE PRECISION array, dimension (N-2) */
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
/* > B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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
/* > X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
/* > On entry, the solution matrix X, as computed by DGTTRS. */
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
/* > WORK is DOUBLE PRECISION array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N) */
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
/* > \ingroup doubleGTcomputational */
/* ===================================================================== */
/* Subroutine */
int dgtrfs_(char *trans, integer *n, integer *nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *dlf, doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal * ferr, doublereal *berr, doublereal *work, integer *iwork, integer * info)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    /* Local variables */
    integer i__, j;
    doublereal s;
    integer nz;
    doublereal eps;
    integer kase;
    doublereal safe1, safe2;
    extern logical lsame_(char *, char *);
    integer isave[3];
    extern /* Subroutine */
    int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *), daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    integer count;
    extern /* Subroutine */
    int dlacn2_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int dlagtm_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    logical notran;
    char transn[1];
    extern /* Subroutine */
    int dgttrs_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
    char transt[1];
    doublereal lstres;
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
    --iwork;
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
        xerbla_("DGTRFS", &i__1);
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
        *(unsigned char *)transt = 'T';
    }
    else
    {
        *(unsigned char *)transn = 'T';
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
        dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
        dlagtm_(trans, n, &c__1, &c_b18, &dl[1], &d__[1], &du[1], &x[j * x_dim1 + 1], ldx, &c_b19, &work[*n + 1], n);
        /* Compute f2c_abs(op(A))*f2c_abs(x) + f2c_abs(b) for use in the backward */
        /* error bound. */
        if (notran)
        {
            if (*n == 1)
            {
                work[1] = (d__1 = b[j * b_dim1 + 1], f2c_abs(d__1)) + (d__2 = d__[ 1] * x[j * x_dim1 + 1], f2c_abs(d__2));
            }
            else
            {
                work[1] = (d__1 = b[j * b_dim1 + 1], f2c_abs(d__1)) + (d__2 = d__[ 1] * x[j * x_dim1 + 1], f2c_abs(d__2)) + (d__3 = du[1] * x[j * x_dim1 + 2], f2c_abs(d__3));
                i__2 = *n - 1;
                for (i__ = 2;
                        i__ <= i__2;
                        ++i__)
                {
                    work[i__] = (d__1 = b[i__ + j * b_dim1], f2c_abs(d__1)) + ( d__2 = dl[i__ - 1] * x[i__ - 1 + j * x_dim1], f2c_abs( d__2)) + (d__3 = d__[i__] * x[i__ + j * x_dim1], f2c_abs(d__3)) + (d__4 = du[i__] * x[i__ + 1 + j * x_dim1], f2c_abs(d__4));
                    /* L30: */
                }
                work[*n] = (d__1 = b[*n + j * b_dim1], f2c_abs(d__1)) + (d__2 = dl[*n - 1] * x[*n - 1 + j * x_dim1], f2c_abs(d__2)) + ( d__3 = d__[*n] * x[*n + j * x_dim1], f2c_abs(d__3));
            }
        }
        else
        {
            if (*n == 1)
            {
                work[1] = (d__1 = b[j * b_dim1 + 1], f2c_abs(d__1)) + (d__2 = d__[ 1] * x[j * x_dim1 + 1], f2c_abs(d__2));
            }
            else
            {
                work[1] = (d__1 = b[j * b_dim1 + 1], f2c_abs(d__1)) + (d__2 = d__[ 1] * x[j * x_dim1 + 1], f2c_abs(d__2)) + (d__3 = dl[1] * x[j * x_dim1 + 2], f2c_abs(d__3));
                i__2 = *n - 1;
                for (i__ = 2;
                        i__ <= i__2;
                        ++i__)
                {
                    work[i__] = (d__1 = b[i__ + j * b_dim1], f2c_abs(d__1)) + ( d__2 = du[i__ - 1] * x[i__ - 1 + j * x_dim1], f2c_abs( d__2)) + (d__3 = d__[i__] * x[i__ + j * x_dim1], f2c_abs(d__3)) + (d__4 = dl[i__] * x[i__ + 1 + j * x_dim1], f2c_abs(d__4));
                    /* L40: */
                }
                work[*n] = (d__1 = b[*n + j * b_dim1], f2c_abs(d__1)) + (d__2 = du[*n - 1] * x[*n - 1 + j * x_dim1], f2c_abs(d__2)) + ( d__3 = d__[*n] * x[*n + j * x_dim1], f2c_abs(d__3));
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
            if (work[i__] > safe2)
            {
                /* Computing MAX */
                d__2 = s;
                d__3 = (d__1 = work[*n + i__], f2c_abs(d__1)) / work[ i__]; // , expr subst
                s = max(d__2,d__3);
            }
            else
            {
                /* Computing MAX */
                d__2 = s;
                d__3 = ((d__1 = work[*n + i__], f2c_abs(d__1)) + safe1) / (work[i__] + safe1); // , expr subst
                s = max(d__2,d__3);
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
            dgttrs_(trans, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[ 1], &work[*n + 1], n, info);
            daxpy_(n, &c_b19, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1) ;
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
        /* Use DLACN2 to estimate the infinity-norm of the matrix */
        /* inv(op(A)) * diag(W), */
        /* where W = f2c_abs(R) + NZ*EPS*( f2c_abs(op(A))*f2c_abs(X)+f2c_abs(B) ))) */
        i__2 = *n;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            if (work[i__] > safe2)
            {
                work[i__] = (d__1 = work[*n + i__], f2c_abs(d__1)) + nz * eps * work[i__];
            }
            else
            {
                work[i__] = (d__1 = work[*n + i__], f2c_abs(d__1)) + nz * eps * work[i__] + safe1;
            }
            /* L60: */
        }
        kase = 0;
L70:
        dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], & kase, isave);
        if (kase != 0)
        {
            if (kase == 1)
            {
                /* Multiply by diag(W)*inv(op(A)**T). */
                dgttrs_(transt, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], & ipiv[1], &work[*n + 1], n, info);
                i__2 = *n;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    work[*n + i__] = work[i__] * work[*n + i__];
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
                    work[*n + i__] = work[i__] * work[*n + i__];
                    /* L90: */
                }
                dgttrs_(transn, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], & ipiv[1], &work[*n + 1], n, info);
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
            d__2 = lstres;
            d__3 = (d__1 = x[i__ + j * x_dim1], f2c_abs(d__1)); // , expr subst
            lstres = max(d__2,d__3);
            /* L100: */
        }
        if (lstres != 0.)
        {
            ferr[j] /= lstres;
        }
        /* L110: */
    }
    return 0;
    /* End of DGTRFS */
}
/* dgtrfs_ */
