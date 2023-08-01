/* ../netlib/sptsvx.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief <b> SPTSVX computes the solution to system of linear equations A * X = B for PT matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SPTSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sptsvx. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sptsvx. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sptsvx. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, */
/* RCOND, FERR, BERR, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER FACT */
/* INTEGER INFO, LDB, LDX, N, NRHS */
/* REAL RCOND */
/* .. */
/* .. Array Arguments .. */
/* REAL B( LDB, * ), BERR( * ), D( * ), DF( * ), */
/* $ E( * ), EF( * ), FERR( * ), WORK( * ), */
/* $ X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPTSVX uses the factorization A = L*D*L**T to compute the solution */
/* > to a real system of linear equations A*X = B, where A is an N-by-N */
/* > symmetric positive definite tridiagonal matrix and X and B are */
/* > N-by-NRHS matrices. */
/* > */
/* > Error bounds on the solution and a condition estimate are also */
/* > provided. */
/* > \endverbatim */
/* > \par Description: */
/* ================= */
/* > */
/* > \verbatim */
/* > */
/* > The following steps are performed: */
/* > */
/* > 1. If FACT = 'N', the matrix A is factored as A = L*D*L**T, where L */
/* > is a unit lower bidiagonal matrix and D is diagonal. The */
/* > factorization can also be regarded as having the form */
/* > A = U**T*D*U. */
/* > */
/* > 2. If the leading i-by-i principal minor is not positive definite, */
/* > then the routine returns with INFO = i. Otherwise, the factored */
/* > form of A is used to estimate the condition number of the matrix */
/* > A. If the reciprocal of the condition number is less than machine */
/* > precision, INFO = N+1 is returned as a warning, but the routine */
/* > still goes on to solve for X and compute error bounds as */
/* > described below. */
/* > */
/* > 3. The system of equations is solved for X using the factored form */
/* > of A. */
/* > */
/* > 4. Iterative refinement is applied to improve the computed solution */
/* > matrix and calculate error bounds and backward error estimates */
/* > for it. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] FACT */
/* > \verbatim */
/* > FACT is CHARACTER*1 */
/* > Specifies whether or not the factored form of A has been */
/* > supplied on entry. */
/* > = 'F': On entry, DF and EF contain the factored form of A. */
/* > D, E, DF, and EF will not be modified. */
/* > = 'N': The matrix A will be copied to DF and EF and */
/* > factored. */
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
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The n diagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > The (n-1) subdiagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DF */
/* > \verbatim */
/* > DF is REAL array, dimension (N) */
/* > If FACT = 'F', then DF is an input argument and on entry */
/* > contains the n diagonal elements of the diagonal matrix D */
/* > from the L*D*L**T factorization of A. */
/* > If FACT = 'N', then DF is an output argument and on exit */
/* > contains the n diagonal elements of the diagonal matrix D */
/* > from the L*D*L**T factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EF */
/* > \verbatim */
/* > EF is REAL array, dimension (N-1) */
/* > If FACT = 'F', then EF is an input argument and on entry */
/* > contains the (n-1) subdiagonal elements of the unit */
/* > bidiagonal factor L from the L*D*L**T factorization of A. */
/* > If FACT = 'N', then EF is an output argument and on exit */
/* > contains the (n-1) subdiagonal elements of the unit */
/* > bidiagonal factor L from the L*D*L**T factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,NRHS) */
/* > The N-by-NRHS right hand side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* > X is REAL array, dimension (LDX,NRHS) */
/* > If INFO = 0 of INFO = N+1, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > The reciprocal condition number of the matrix A. If RCOND */
/* > is less than the machine precision (in particular, if */
/* > RCOND = 0), the matrix is singular to working precision. */
/* > This condition is indicated by a return code of INFO > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* > FERR is REAL array, dimension (NRHS) */
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
/* > BERR is REAL array, dimension (NRHS) */
/* > The componentwise relative backward error of each solution */
/* > vector X(j) (i.e., the smallest relative change in any */
/* > element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, and i is */
/* > <= N: the leading minor of order i of A is */
/* > not positive definite, so the factorization */
/* > could not be completed, and the solution has not */
/* > been computed. RCOND = 0 is returned. */
/* > = N+1: U is nonsingular, but RCOND is less than machine */
/* > precision, meaning that the matrix is singular */
/* > to working precision. Nevertheless, the */
/* > solution and error bounds are computed because */
/* > there are a number of situations where the */
/* > computed solution can be more accurate than the */
/* > value of RCOND would suggest. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realPTsolve */
/* ===================================================================== */
/* Subroutine */
int sptsvx_(char *fact, integer *n, integer *nrhs, real *d__, real *e, real *df, real *ef, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"sptsvx inputs: fact %c, n %" FLA_IS ", nrhs %" FLA_IS ", ldb %" FLA_IS ", ldx %" FLA_IS "",*fact, *n, *nrhs, *ldb, *ldx);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1;
    /* Local variables */
    extern logical lsame_(char *, char *);
    real anorm;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *);
    extern real slamch_(char *);
    logical nofact;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len), slacpy_( char *, integer *, integer *, real *, integer *, real *, integer * );
    extern real slanst_(char *, integer *, real *, real *);
    extern /* Subroutine */
    int sptcon_(integer *, real *, real *, real *, real *, real *, integer *), sptrfs_(integer *, integer *, real *, real *, real *, real *, real *, integer *, real *, integer *, real *, real *, real *, integer *), spttrf_(integer *, real *, real *, integer *), spttrs_(integer *, integer *, real *, real *, real *, integer *, integer *);
    /* -- LAPACK driver routine (version 3.4.2) -- */
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
    nofact = lsame_(fact, "N");
    if (! nofact && ! lsame_(fact, "F"))
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
    else if (*ldb < fla_max(1,*n))
    {
        *info = -9;
    }
    else if (*ldx < fla_max(1,*n))
    {
        *info = -11;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SPTSVX", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    if (nofact)
    {
        /* Compute the L*D*L**T (or U**T*D*U) factorization of A. */
        scopy_(n, &d__[1], &c__1, &df[1], &c__1);
        if (*n > 1)
        {
            i__1 = *n - 1;
            scopy_(&i__1, &e[1], &c__1, &ef[1], &c__1);
        }
        spttrf_(n, &df[1], &ef[1], info);
        /* Return if INFO is non-zero. */
        if (*info > 0)
        {
            *rcond = 0.f;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return 0;
        }
    }
    /* Compute the norm of the matrix A. */
    anorm = slanst_("1", n, &d__[1], &e[1]);
    /* Compute the reciprocal of the condition number of A. */
    sptcon_(n, &df[1], &ef[1], &anorm, rcond, &work[1], info);
    /* Compute the solution vectors X. */
    slacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    spttrs_(n, nrhs, &df[1], &ef[1], &x[x_offset], ldx, info);
    /* Use iterative refinement to improve the computed solutions and */
    /* compute error bounds and backward error estimates for them. */
    sptrfs_(n, nrhs, &d__[1], &e[1], &df[1], &ef[1], &b[b_offset], ldb, &x[ x_offset], ldx, &ferr[1], &berr[1], &work[1], info);
    /* Set INFO = N+1 if the matrix is singular to working precision. */
    if (*rcond < slamch_("Epsilon"))
    {
        *info = *n + 1;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of SPTSVX */
}
/* sptsvx_ */
