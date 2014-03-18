/* ../netlib/csysvx.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief <b> CSYSVX computes the solution to system of linear equations A * X = B for SY matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSYSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csysvx. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csysvx. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csysvx. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, */
/* LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, */
/* RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER FACT, UPLO */
/* INTEGER INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS */
/* REAL RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL BERR( * ), FERR( * ), RWORK( * ) */
/* COMPLEX A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/* $ WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYSVX uses the diagonal pivoting factorization to compute the */
/* > solution to a complex system of linear equations A * X = B, */
/* > where A is an N-by-N symmetric matrix and X and B are N-by-NRHS */
/* > matrices. */
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
/* > 1. If FACT = 'N', the diagonal pivoting method is used to factor A. */
/* > The form of the factorization is */
/* > A = U * D * U**T, if UPLO = 'U', or */
/* > A = L * D * L**T, if UPLO = 'L', */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is symmetric and block diagonal with */
/* > 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > 2. If some D(i,i)=0, so that D is exactly singular, then the routine */
/* > returns with INFO = i. Otherwise, the factored form of A is used */
/* > to estimate the condition number of the matrix A. If the */
/* > reciprocal of the condition number is less than machine precision, */
/* > INFO = N+1 is returned as a warning, but the routine still goes on */
/* > to solve for X and compute error bounds as described below. */
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
/* > = 'F': On entry, AF and IPIV contain the factored form */
/* > of A. A, AF and IPIV will not be modified. */
/* > = 'N': The matrix A will be copied to AF and factored. */
/* > \endverbatim */
/* > */
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
/* > of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
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
/* > \param[in,out] AF */
/* > \verbatim */
/* > AF is COMPLEX array, dimension (LDAF,N) */
/* > If FACT = 'F', then AF is an input argument and on entry */
/* > contains the block diagonal matrix D and the multipliers used */
/* > to obtain the factor U or L from the factorization */
/* > A = U*D*U**T or A = L*D*L**T as computed by CSYTRF. */
/* > */
/* > If FACT = 'N', then AF is an output argument and on exit */
/* > returns the block diagonal matrix D and the multipliers used */
/* > to obtain the factor U or L from the factorization */
/* > A = U*D*U**T or A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > If FACT = 'F', then IPIV is an input argument and on entry */
/* > contains details of the interchanges and the block structure */
/* > of D, as determined by CSYTRF. */
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* > interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and */
/* > columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* > is a 2-by-2 diagonal block. If UPLO = 'L' and IPIV(k) = */
/* > IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were */
/* > interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > */
/* > If FACT = 'N', then IPIV is an output argument and on exit */
/* > contains details of the interchanges and the block structure */
/* > of D, as determined by CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
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
/* > X is COMPLEX array, dimension (LDX,NRHS) */
/* > If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > The estimate of the reciprocal condition number of the matrix */
/* > A. If RCOND is less than the machine precision (in */
/* > particular, if RCOND = 0), the matrix is singular to working */
/* > precision. This condition is indicated by a return code of */
/* > INFO > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* > FERR is REAL array, dimension (NRHS) */
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
/* > BERR is REAL array, dimension (NRHS) */
/* > The componentwise relative backward error of each solution */
/* > vector X(j) (i.e., the smallest relative change in */
/* > any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of WORK. LWORK >= max(1,2*N), and for best */
/* > performance, when FACT = 'N', LWORK >= max(1,2*N,N*NB), where */
/* > NB is the optimal blocksize for CSYTRF. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, and i is */
/* > <= N: D(i,i) is exactly zero. The factorization */
/* > has been completed but the factor D is exactly */
/* > singular, so the solution and error bounds could */
/* > not be computed. RCOND = 0 is returned. */
/* > = N+1: D is nonsingular, but RCOND is less than machine */
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
/* > \date April 2012 */
/* > \ingroup complexSYsolve */
/* ===================================================================== */
/* Subroutine */
int csysvx_(char *fact, char *uplo, integer *n, integer * nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer * ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, integer *lwork, real *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    /* Local variables */
    integer nb;
    extern logical lsame_(char *, char *);
    real anorm;
    extern real slamch_(char *);
    logical nofact;
    extern /* Subroutine */
    int clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern real clansy_(char *, char *, integer *, complex *, integer *, real *);
    extern /* Subroutine */
    int csycon_(char *, integer *, complex *, integer *, integer *, real *, real *, complex *, integer *), csyrfs_(char *, integer *, integer *, complex *, integer *, complex *, integer *, integer *, complex *, integer *, complex *, integer *, real *, real *, complex *, real *, integer *), csytrf_(char *, integer *, complex *, integer *, integer *, complex *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */
    int csytrs_(char *, integer *, integer *, complex *, integer *, integer *, complex *, integer *, integer *);
    /* -- LAPACK driver routine (version 3.4.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* April 2012 */
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
    nofact = lsame_(fact, "N");
    lquery = *lwork == -1;
    if (! nofact && ! lsame_(fact, "F"))
    {
        *info = -1;
    }
    else if (! lsame_(uplo, "U") && ! lsame_(uplo, "L"))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*nrhs < 0)
    {
        *info = -4;
    }
    else if (*lda < max(1,*n))
    {
        *info = -6;
    }
    else if (*ldaf < max(1,*n))
    {
        *info = -8;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -11;
    }
    else if (*ldx < max(1,*n))
    {
        *info = -13;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 1; // , expr subst
        if (*lwork < max(i__1,i__2) && ! lquery)
        {
            *info = -18;
        }
    }
    if (*info == 0)
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 1; // , expr subst
        lwkopt = max(i__1,i__2);
        if (nofact)
        {
            nb = ilaenv_(&c__1, "CSYTRF", uplo, n, &c_n1, &c_n1, &c_n1);
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = *n * nb; // , expr subst
            lwkopt = max(i__1,i__2);
        }
        work[1].r = (real) lwkopt;
        work[1].i = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CSYSVX", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    if (nofact)
    {
        /* Compute the factorization A = U*D*U**T or A = L*D*L**T. */
        clacpy_(uplo, n, n, &a[a_offset], lda, &af[af_offset], ldaf);
        csytrf_(uplo, n, &af[af_offset], ldaf, &ipiv[1], &work[1], lwork, info);
        /* Return if INFO is non-zero. */
        if (*info > 0)
        {
            *rcond = 0.f;
            return 0;
        }
    }
    /* Compute the norm of the matrix A. */
    anorm = clansy_("I", uplo, n, &a[a_offset], lda, &rwork[1]);
    /* Compute the reciprocal of the condition number of A. */
    csycon_(uplo, n, &af[af_offset], ldaf, &ipiv[1], &anorm, rcond, &work[1], info);
    /* Compute the solution vectors X. */
    clacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    csytrs_(uplo, n, nrhs, &af[af_offset], ldaf, &ipiv[1], &x[x_offset], ldx, info);
    /* Use iterative refinement to improve the computed solutions and */
    /* compute error bounds and backward error estimates for them. */
    csyrfs_(uplo, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1] , &rwork[1], info);
    /* Set INFO = N+1 if the matrix is singular to working precision. */
    if (*rcond < slamch_("Epsilon"))
    {
        *info = *n + 1;
    }
    work[1].r = (real) lwkopt;
    work[1].i = 0.f; // , expr subst
    return 0;
    /* End of CSYSVX */
}
/* csysvx_ */
