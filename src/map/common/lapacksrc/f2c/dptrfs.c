/* ../netlib/dptrfs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b11 = 1.;
/* > \brief \b DPTRFS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DPTRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dptrfs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dptrfs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dptrfs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, */
/* BERR, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION B( LDB, * ), BERR( * ), D( * ), DF( * ), */
/* $ E( * ), EF( * ), FERR( * ), WORK( * ), */
/* $ X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPTRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric positive definite */
/* > and tridiagonal, and provides error bounds and backward error */
/* > estimates for the solution. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > The n diagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N-1) */
/* > The (n-1) subdiagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] DF */
/* > \verbatim */
/* > DF is DOUBLE PRECISION array, dimension (N) */
/* > The n diagonal elements of the diagonal matrix D from the */
/* > factorization computed by DPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] EF */
/* > \verbatim */
/* > EF is DOUBLE PRECISION array, dimension (N-1) */
/* > The (n-1) subdiagonal elements of the unit bidiagonal factor */
/* > L from the factorization computed by DPTTRF. */
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
/* > On entry, the solution matrix X, as computed by DPTTRS. */
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
/* > WORK is DOUBLE PRECISION array, dimension (2*N) */
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
/* > \ingroup doublePTcomputational */
/* ===================================================================== */
/* Subroutine */
int dptrfs_(integer *n, integer *nrhs, doublereal *d__, doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, doublereal *work, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;
    /* Local variables */
    integer i__, j;
    doublereal s, bi, cx, dx, ex;
    integer ix, nz;
    doublereal eps, safe1, safe2;
    extern /* Subroutine */
    int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    integer count;
    extern doublereal dlamch_(char *);
    extern integer idamax_(integer *, doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal lstres;
    extern /* Subroutine */
    int dpttrs_(integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, integer *);
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
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
    /* Function Body */
    *info = 0;
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*nrhs < 0)
    {
        *info = -2;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -8;
    }
    else if (*ldx < max(1,*n))
    {
        *info = -10;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DPTRFS", &i__1);
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
        if (*n == 1)
        {
            bi = b[j * b_dim1 + 1];
            dx = d__[1] * x[j * x_dim1 + 1];
            work[*n + 1] = bi - dx;
            work[1] = f2c_abs(bi) + f2c_abs(dx);
        }
        else
        {
            bi = b[j * b_dim1 + 1];
            dx = d__[1] * x[j * x_dim1 + 1];
            ex = e[1] * x[j * x_dim1 + 2];
            work[*n + 1] = bi - dx - ex;
            work[1] = f2c_abs(bi) + f2c_abs(dx) + f2c_abs(ex);
            i__2 = *n - 1;
            for (i__ = 2;
                    i__ <= i__2;
                    ++i__)
            {
                bi = b[i__ + j * b_dim1];
                cx = e[i__ - 1] * x[i__ - 1 + j * x_dim1];
                dx = d__[i__] * x[i__ + j * x_dim1];
                ex = e[i__] * x[i__ + 1 + j * x_dim1];
                work[*n + i__] = bi - cx - dx - ex;
                work[i__] = f2c_abs(bi) + f2c_abs(cx) + f2c_abs(dx) + f2c_abs(ex);
                /* L30: */
            }
            bi = b[*n + j * b_dim1];
            cx = e[*n - 1] * x[*n - 1 + j * x_dim1];
            dx = d__[*n] * x[*n + j * x_dim1];
            work[*n + *n] = bi - cx - dx;
            work[*n] = f2c_abs(bi) + f2c_abs(cx) + f2c_abs(dx);
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
            /* L40: */
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
            dpttrs_(n, &c__1, &df[1], &ef[1], &work[*n + 1], n, info);
            daxpy_(n, &c_b11, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1) ;
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
            if (work[i__] > safe2)
            {
                work[i__] = (d__1 = work[*n + i__], f2c_abs(d__1)) + nz * eps * work[i__];
            }
            else
            {
                work[i__] = (d__1 = work[*n + i__], f2c_abs(d__1)) + nz * eps * work[i__] + safe1;
            }
            /* L50: */
        }
        ix = idamax_(n, &work[1], &c__1);
        ferr[j] = work[ix];
        /* Estimate the norm of inv(A). */
        /* Solve M(A) * x = e, where M(A) = (m(i,j)) is given by */
        /* m(i,j) = f2c_abs(A(i,j)); i = j; */
        /* m(i,j) = -f2c_abs(A(i,j)), i .ne. j, */
        /* and e = [ 1, 1, ..., 1 ]**T. Note M(A) = M(L)*D*M(L)**T. */
        /* Solve M(L) * x = e. */
        work[1] = 1.;
        i__2 = *n;
        for (i__ = 2;
                i__ <= i__2;
                ++i__)
        {
            work[i__] = work[i__ - 1] * (d__1 = ef[i__ - 1], f2c_abs(d__1)) + 1.;
            /* L60: */
        }
        /* Solve D * M(L)**T * x = b. */
        work[*n] /= df[*n];
        for (i__ = *n - 1;
                i__ >= 1;
                --i__)
        {
            work[i__] = work[i__] / df[i__] + work[i__ + 1] * (d__1 = ef[i__], f2c_abs(d__1));
            /* L70: */
        }
        /* Compute norm(inv(A)) = max(x(i)), 1<=i<=n. */
        ix = idamax_(n, &work[1], &c__1);
        ferr[j] *= (d__1 = work[ix], f2c_abs(d__1));
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
            /* L80: */
        }
        if (lstres != 0.)
        {
            ferr[j] /= lstres;
        }
        /* L90: */
    }
    return 0;
    /* End of DPTRFS */
}
/* dptrfs_ */
