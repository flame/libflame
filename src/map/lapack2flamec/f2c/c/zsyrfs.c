/* ../netlib/zsyrfs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    1.,0.
}
;
static integer c__1 = 1;
/* > \brief \b ZSYRFS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZSYRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyrfs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyrfs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyrfs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, */
/* X, LDX, FERR, BERR, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* DOUBLE PRECISION BERR( * ), FERR( * ), RWORK( * ) */
/* COMPLEX*16 A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/* $ WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric indefinite, and */
/* > provides error bounds and backward error estimates for the solution. */
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
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > The symmetric matrix A. If UPLO = 'U', the leading N-by-N */
/* > upper triangular part of A contains the upper triangular part */
/* > of the matrix A, and the strictly lower triangular part of A */
/* > is not referenced. If UPLO = 'L', the leading N-by-N lower */
/* > triangular part of A contains the lower triangular part of */
/* > the matrix A, and the strictly upper triangular part of A is */
/* > not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* > AF is COMPLEX*16 array, dimension (LDAF,N) */
/* > The factored form of the matrix A. AF contains the block */
/* > diagonal matrix D and the multipliers used to obtain the */
/* > factor U or L from the factorization A = U*D*U**T or */
/* > A = L*D*L**T as computed by ZSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by ZSYTRF. */
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
/* > On entry, the solution matrix X, as computed by ZSYTRS. */
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
/* > \date November 2011 */
/* > \ingroup complex16SYcomputational */
/* ===================================================================== */
/* Subroutine */
int zsyrfs_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    integer i__, j, k;
    doublereal s, xk;
    integer nz;
    doublereal eps;
    integer kase;
    doublereal safe1, safe2;
    extern logical lsame_(char *, char *);
    integer isave[3], count;
    logical upper;
    extern /* Subroutine */
    int zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), zsymv_( char *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *), zlacn2_(integer *, doublecomplex *, doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    doublereal safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal lstres;
    extern /* Subroutine */
    int zsytrs_(char *, integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *, integer *);
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
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
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*ldaf < max(1,*n))
    {
        *info = -7;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -10;
    }
    else if (*ldx < max(1,*n))
    {
        *info = -12;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZSYRFS", &i__1);
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
        zcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
        z__1.r = -1.;
        z__1.i = -0.; // , expr subst
        zsymv_(uplo, n, &z__1, &a[a_offset], lda, &x[j * x_dim1 + 1], &c__1, & c_b1, &work[1], &c__1);
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
            i__3 = i__ + j * b_dim1;
            rwork[i__] = (d__1 = b[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&b[ i__ + j * b_dim1]), f2c_abs(d__2));
            /* L30: */
        }
        /* Compute f2c_abs(A)*f2c_abs(X) + f2c_abs(B). */
        if (upper)
        {
            i__2 = *n;
            for (k = 1;
                    k <= i__2;
                    ++k)
            {
                s = 0.;
                i__3 = k + j * x_dim1;
                xk = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[k + j * x_dim1]), f2c_abs(d__2));
                i__3 = k - 1;
                for (i__ = 1;
                        i__ <= i__3;
                        ++i__)
                {
                    i__4 = i__ + k * a_dim1;
                    rwork[i__] += ((d__1 = a[i__4].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_abs(d__2))) * xk;
                    i__4 = i__ + k * a_dim1;
                    i__5 = i__ + j * x_dim1;
                    s += ((d__1 = a[i__4].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[ i__ + k * a_dim1]), f2c_abs(d__2))) * ((d__3 = x[i__5] .r, f2c_abs(d__3)) + (d__4 = d_imag(&x[i__ + j * x_dim1]), f2c_abs(d__4)));
                    /* L40: */
                }
                i__3 = k + k * a_dim1;
                rwork[k] = rwork[k] + ((d__1 = a[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), f2c_abs(d__2))) * xk + s;
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
                i__3 = k + j * x_dim1;
                xk = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[k + j * x_dim1]), f2c_abs(d__2));
                i__3 = k + k * a_dim1;
                rwork[k] += ((d__1 = a[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(& a[k + k * a_dim1]), f2c_abs(d__2))) * xk;
                i__3 = *n;
                for (i__ = k + 1;
                        i__ <= i__3;
                        ++i__)
                {
                    i__4 = i__ + k * a_dim1;
                    rwork[i__] += ((d__1 = a[i__4].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_abs(d__2))) * xk;
                    i__4 = i__ + k * a_dim1;
                    i__5 = i__ + j * x_dim1;
                    s += ((d__1 = a[i__4].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[ i__ + k * a_dim1]), f2c_abs(d__2))) * ((d__3 = x[i__5] .r, f2c_abs(d__3)) + (d__4 = d_imag(&x[i__ + j * x_dim1]), f2c_abs(d__4)));
                    /* L60: */
                }
                rwork[k] += s;
                /* L70: */
            }
        }
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
            zsytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], n, info);
            zaxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
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
        /* Use ZLACN2 to estimate the infinity-norm of the matrix */
        /* inv(A) * diag(W), */
        /* where W = f2c_abs(R) + NZ*EPS*( f2c_abs(A)*f2c_abs(X)+f2c_abs(B) ))) */
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
            /* L90: */
        }
        kase = 0;
L100:
        zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
        if (kase != 0)
        {
            if (kase == 1)
            {
                /* Multiply by diag(W)*inv(A**T). */
                zsytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[ 1], n, info);
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
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__;
                    z__1.r = rwork[i__4] * work[i__5].r;
                    z__1.i = rwork[i__4] * work[i__5].i; // , expr subst
                    work[i__3].r = z__1.r;
                    work[i__3].i = z__1.i; // , expr subst
                    /* L120: */
                }
                zsytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[ 1], n, info);
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
            i__3 = i__ + j * x_dim1;
            d__3 = lstres;
            d__4 = (d__1 = x[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[i__ + j * x_dim1]), f2c_abs(d__2)); // , expr subst
            lstres = max(d__3,d__4);
            /* L130: */
        }
        if (lstres != 0.)
        {
            ferr[j] /= lstres;
        }
        /* L140: */
    }
    return 0;
    /* End of ZSYRFS */
}
/* zsyrfs_ */
