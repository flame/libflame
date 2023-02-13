/* ../netlib/zcposv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    -1.,0.
}
;
static doublecomplex c_b2 =
{
    1.,0.
}
;
static integer c__1 = 1;
/* > \brief <b> ZCPOSV computes the solution to system of linear equations A * X = B for PO matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZCPOSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zcposv. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zcposv. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zcposv. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZCPOSV( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, WORK, */
/* SWORK, RWORK, ITER, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, ITER, LDA, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION RWORK( * ) */
/* COMPLEX SWORK( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), WORK( N, * ), */
/* $ X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZCPOSV computes the solution to a complex system of linear equations */
/* > A * X = B, */
/* > where A is an N-by-N Hermitian positive definite matrix and X and B */
/* > are N-by-NRHS matrices. */
/* > */
/* > ZCPOSV first attempts to factorize the matrix in COMPLEX and use this */
/* > factorization within an iterative refinement procedure to produce a */
/* > solution with COMPLEX*16 normwise backward error quality (see below). */
/* > If the approach fails the method switches to a COMPLEX*16 */
/* > factorization and solve. */
/* > */
/* > The iterative refinement is not going to be a winning strategy if */
/* > the ratio COMPLEX performance over COMPLEX*16 performance is too */
/* > small. A reasonable strategy should take the number of right-hand */
/* > sides and the size of the matrix into account. This might be done */
/* > with a call to ILAENV in the future. Up to now, we always try */
/* > iterative refinement. */
/* > */
/* > The iterative refinement process is stopped if */
/* > ITER > ITERMAX */
/* > or for all the RHS we have: */
/* > RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX */
/* > where */
/* > o ITER is the number of the current iteration in the iterative */
/* > refinement process */
/* > o RNRM is the infinity-norm of the residual */
/* > o XNRM is the infinity-norm of the solution */
/* > o ANRM is the infinity-operator-norm of the matrix A */
/* > o EPS is the machine epsilon returned by DLAMCH('Epsilon') */
/* > The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00 */
/* > respectively. */
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
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, */
/* > dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > Note that the imaginary parts of the diagonal */
/* > elements need not be set and are assumed to be zero. */
/* > */
/* > On exit, if iterative refinement has been successfully used */
/* > (INFO.EQ.0 and ITER.GE.0, see description below), then A is */
/* > unchanged, if double precision factorization has been used */
/* > (INFO.EQ.0 and ITER.LT.0, see description below), then the */
/* > array A contains the factor U or L from the Cholesky */
/* > factorization A = U**H*U or A = L*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > The N-by-NRHS right hand side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (LDX,NRHS) */
/* > If INFO = 0, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N*NRHS) */
/* > This array is used to hold the residual vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] SWORK */
/* > \verbatim */
/* > SWORK is COMPLEX array, dimension (N*(N+NRHS)) */
/* > This array is used to use the single precision matrix and the */
/* > right-hand sides or solutions in single precision. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] ITER */
/* > \verbatim */
/* > ITER is INTEGER */
/* > < 0: iterative refinement has failed, COMPLEX*16 */
/* > factorization has been performed */
/* > -1 : the routine fell back to full precision for */
/* > implementation- or machine-specific reasons */
/* > -2 : narrowing the precision induced an overflow, */
/* > the routine fell back to full precision */
/* > -3 : failure of CPOTRF */
/* > -31: stop the iterative refinement after the 30th */
/* > iterations */
/* > > 0: iterative refinement has been sucessfully used. */
/* > Returns the number of iterations */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the leading minor of order i of */
/* > (COMPLEX*16) A is not positive definite, so the */
/* > factorization could not be completed, and the solution */
/* > has not been computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complex16POsolve */
/* ===================================================================== */
/* Subroutine */
int zcposv_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublecomplex *work, complex *swork, doublereal *rwork, integer *iter, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, work_dim1, work_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    /* Local variables */
    integer i__;
    doublereal cte, eps, anrm;
    integer ptsa;
    doublereal rnrm, xnrm;
    integer ptsx;
    extern logical lsame_(char *, char *);
    integer iiter;
    extern /* Subroutine */
    int zhemm_(char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), zlag2c_(integer *, integer *, doublecomplex *, integer *, complex *, integer *, integer *), clag2z_(integer *, integer *, complex *, integer *, doublecomplex *, integer *, integer *), zlat2c_(char *, integer *, doublecomplex *, integer *, complex *, integer *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern doublereal zlanhe_(char *, char *, integer *, doublecomplex *, integer *, doublereal *);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */
    int cpotrf_(char *, integer *, complex *, integer *, integer *), zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *), cpotrs_(char *, integer *, integer *, complex *, integer *, complex *, integer *, integer *), zpotrf_(char *, integer *, doublecomplex *, integer *, integer *), zpotrs_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex * , integer *, integer *);
    /* -- LAPACK driver routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    work_dim1 = *n;
    work_offset = 1 + work_dim1;
    work -= work_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --swork;
    --rwork;
    /* Function Body */
    *info = 0;
    *iter = 0;
    /* Test the input parameters. */
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L"))
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
        xerbla_("ZCPOSV", &i__1);
        return 0;
    }
    /* Quick return if (N.EQ.0). */
    if (*n == 0)
    {
        return 0;
    }
    /* Skip single precision iterative refinement if a priori slower */
    /* than double precision factorization. */
    if (FALSE_)
    {
        *iter = -1;
        goto L40;
    }
    /* Compute some constants. */
    anrm = zlanhe_("I", uplo, n, &a[a_offset], lda, &rwork[1]);
    eps = dlamch_("Epsilon");
    cte = anrm * eps * sqrt((doublereal) (*n)) * 1.;
    /* Set the indices PTSA, PTSX for referencing SA and SX in SWORK. */
    ptsa = 1;
    ptsx = ptsa + *n * *n;
    /* Convert B from double precision to single precision and store the */
    /* result in SX. */
    zlag2c_(n, nrhs, &b[b_offset], ldb, &swork[ptsx], n, info);
    if (*info != 0)
    {
        *iter = -2;
        goto L40;
    }
    /* Convert A from double precision to single precision and store the */
    /* result in SA. */
    zlat2c_(uplo, n, &a[a_offset], lda, &swork[ptsa], n, info);
    if (*info != 0)
    {
        *iter = -2;
        goto L40;
    }
    /* Compute the Cholesky factorization of SA. */
    cpotrf_(uplo, n, &swork[ptsa], n, info);
    if (*info != 0)
    {
        *iter = -3;
        goto L40;
    }
    /* Solve the system SA*SX = SB. */
    cpotrs_(uplo, n, nrhs, &swork[ptsa], n, &swork[ptsx], n, info);
    /* Convert SX back to COMPLEX*16 */
    clag2z_(n, nrhs, &swork[ptsx], n, &x[x_offset], ldx, info);
    /* Compute R = B - AX (R is WORK). */
    zlacpy_("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n);
    zhemm_("Left", uplo, n, nrhs, &c_b1, &a[a_offset], lda, &x[x_offset], ldx, &c_b2, &work[work_offset], n);
    /* Check whether the NRHS normwise backward errors satisfy the */
    /* stopping criterion. If yes, set ITER=0 and return. */
    i__1 = *nrhs;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = izamax_(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1;
        xnrm = (d__1 = x[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[izamax_(n, & x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1]), f2c_abs(d__2));
        i__2 = izamax_(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * work_dim1;
        rnrm = (d__1 = work[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&work[ izamax_(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * work_dim1]), f2c_abs(d__2));
        if (rnrm > xnrm * cte)
        {
            goto L10;
        }
    }
    /* If we are here, the NRHS normwise backward errors satisfy the */
    /* stopping criterion. We are good to exit. */
    *iter = 0;
    return 0;
L10:
    for (iiter = 1;
            iiter <= 30;
            ++iiter)
    {
        /* Convert R (in WORK) from double precision to single precision */
        /* and store the result in SX. */
        zlag2c_(n, nrhs, &work[work_offset], n, &swork[ptsx], n, info);
        if (*info != 0)
        {
            *iter = -2;
            goto L40;
        }
        /* Solve the system SA*SX = SR. */
        cpotrs_(uplo, n, nrhs, &swork[ptsa], n, &swork[ptsx], n, info);
        /* Convert SX back to double precision and update the current */
        /* iterate. */
        clag2z_(n, nrhs, &swork[ptsx], n, &work[work_offset], n, info);
        i__1 = *nrhs;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            zaxpy_(n, &c_b2, &work[i__ * work_dim1 + 1], &c__1, &x[i__ * x_dim1 + 1], &c__1);
        }
        /* Compute R = B - AX (R is WORK). */
        zlacpy_("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n);
        zhemm_("L", uplo, n, nrhs, &c_b1, &a[a_offset], lda, &x[x_offset], ldx, &c_b2, &work[work_offset], n);
        /* Check whether the NRHS normwise backward errors satisfy the */
        /* stopping criterion. If yes, set ITER=IITER>0 and return. */
        i__1 = *nrhs;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = izamax_(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1;
            xnrm = (d__1 = x[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&x[izamax_( n, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1]), f2c_abs( d__2));
            i__2 = izamax_(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * work_dim1;
            rnrm = (d__1 = work[i__2].r, f2c_abs(d__1)) + (d__2 = d_imag(&work[ izamax_(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * work_dim1]), f2c_abs(d__2));
            if (rnrm > xnrm * cte)
            {
                goto L20;
            }
        }
        /* If we are here, the NRHS normwise backward errors satisfy the */
        /* stopping criterion, we are good to exit. */
        *iter = iiter;
        return 0;
L20: /* L30: */
        ;
    }
    /* If we are at this place of the code, this is because we have */
    /* performed ITER=ITERMAX iterations and never satisified the */
    /* stopping criterion, set up the ITER flag accordingly and follow */
    /* up on double precision routine. */
    *iter = -31;
L40: /* Single-precision iterative refinement failed to converge to a */
    /* satisfactory solution, so we resort to double precision. */
    zpotrf_(uplo, n, &a[a_offset], lda, info);
    if (*info != 0)
    {
        return 0;
    }
    zlacpy_("All", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    zpotrs_(uplo, n, nrhs, &a[a_offset], lda, &x[x_offset], ldx, info);
    return 0;
    /* End of ZCPOSV. */
}
/* zcposv_ */
