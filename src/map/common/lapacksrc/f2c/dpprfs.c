/* ../netlib/dpprfs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;
/* > \brief \b DPPRFS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DPPRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpprfs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpprfs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpprfs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, */
/* BERR, WORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION AFP( * ), AP( * ), B( LDB, * ), BERR( * ), */
/* $ FERR( * ), WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPPRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric positive definite */
/* > and packed, and provides error bounds and backward error estimates */
/* > for the solution. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
*/
/* > = 'L': Lower triangle of A is stored. */
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
/* > of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* > AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* > The upper or lower triangle of the symmetric matrix A, packed */
/* > columnwise in a linear array. The j-th column of A is stored */
/* > in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] AFP */
/* > \verbatim */
/* > AFP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* > The triangular factor U or L from the Cholesky factorization */
/* > A = U**T*U or A = L*L**T, as computed by DPPTRF/ZPPTRF, */
/* > packed columnwise in a linear array in the same format as A */
/* > (see AP). */
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
/* > On entry, the solution matrix X, as computed by DPPTRS. */
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
/* > \date November 2011 */
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int dpprfs_(char *uplo, integer *n, integer *nrhs, doublereal *ap, doublereal *afp, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;
    /* Local variables */
    integer i__, j, k;
    doublereal s;
    integer ik, kk;
    doublereal xk;
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
    int dspmv_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    logical upper;
    extern /* Subroutine */
    int dlacn2_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    doublereal safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal lstres;
    extern /* Subroutine */
    int dpptrs_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    --ap;
    --afp;
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
        *info = -7;
    }
    else if (*ldx < max(1,*n))
    {
        *info = -9;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DPPRFS", &i__1);
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
    nz = *n + 1;
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
        /* Compute residual R = B - A * X */
        dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
        dspmv_(uplo, n, &c_b12, &ap[1], &x[j * x_dim1 + 1], &c__1, &c_b14, & work[*n + 1], &c__1);
        /* Compute componentwise relative backward error from formula */
        /* max(i) ( f2c_abs(R(i)) / ( f2c_abs(A)*f2c_abs(X) + f2c_abs(B) )(i) ) */
        /* where f2c_abs(Z) is the componentwise absolute value of the matrix */
        /* or vector Z. If the i-th component of the denominator is less */
        /* than SAFE2, then SAFE1 is added to the i-th components of the */
        /* numerator and denominator before dividing. */
        i__2 = *n;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            work[i__] = (d__1 = b[i__ + j * b_dim1], f2c_abs(d__1));
            /* L30: */
        }
        /* Compute f2c_abs(A)*f2c_abs(X) + f2c_abs(B). */
        kk = 1;
        if (upper)
        {
            i__2 = *n;
            for (k = 1;
                    k <= i__2;
                    ++k)
            {
                s = 0.;
                xk = (d__1 = x[k + j * x_dim1], f2c_abs(d__1));
                ik = kk;
                i__3 = k - 1;
                for (i__ = 1;
                        i__ <= i__3;
                        ++i__)
                {
                    work[i__] += (d__1 = ap[ik], f2c_abs(d__1)) * xk;
                    s += (d__1 = ap[ik], f2c_abs(d__1)) * (d__2 = x[i__ + j * x_dim1], f2c_abs(d__2));
                    ++ik;
                    /* L40: */
                }
                work[k] = work[k] + (d__1 = ap[kk + k - 1], f2c_abs(d__1)) * xk + s;
                kk += k;
                /* L50: */
            }
        }
        else
        {
            i__2 = *n;
            for (k = 1;
                    k <= i__2;
                    ++k)
            {
                s = 0.;
                xk = (d__1 = x[k + j * x_dim1], f2c_abs(d__1));
                work[k] += (d__1 = ap[kk], f2c_abs(d__1)) * xk;
                ik = kk + 1;
                i__3 = *n;
                for (i__ = k + 1;
                        i__ <= i__3;
                        ++i__)
                {
                    work[i__] += (d__1 = ap[ik], f2c_abs(d__1)) * xk;
                    s += (d__1 = ap[ik], f2c_abs(d__1)) * (d__2 = x[i__ + j * x_dim1], f2c_abs(d__2));
                    ++ik;
                    /* L60: */
                }
                work[k] += s;
                kk += *n - k + 1;
                /* L70: */
            }
        }
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
            /* L80: */
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
            dpptrs_(uplo, n, &c__1, &afp[1], &work[*n + 1], n, info);
            daxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1) ;
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
        /* Use DLACN2 to estimate the infinity-norm of the matrix */
        /* inv(A) * diag(W), */
        /* where W = f2c_abs(R) + NZ*EPS*( f2c_abs(A)*f2c_abs(X)+f2c_abs(B) ))) */
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
            /* L90: */
        }
        kase = 0;
L100:
        dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], & kase, isave);
        if (kase != 0)
        {
            if (kase == 1)
            {
                /* Multiply by diag(W)*inv(A**T). */
                dpptrs_(uplo, n, &c__1, &afp[1], &work[*n + 1], n, info);
                i__2 = *n;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    work[*n + i__] = work[i__] * work[*n + i__];
                    /* L110: */
                }
            }
            else if (kase == 2)
            {
                /* Multiply by inv(A)*diag(W). */
                i__2 = *n;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    work[*n + i__] = work[i__] * work[*n + i__];
                    /* L120: */
                }
                dpptrs_(uplo, n, &c__1, &afp[1], &work[*n + 1], n, info);
            }
            goto L100;
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
            /* L130: */
        }
        if (lstres != 0.)
        {
            ferr[j] /= lstres;
        }
        /* L140: */
    }
    return 0;
    /* End of DPPRFS */
}
/* dpprfs_ */
