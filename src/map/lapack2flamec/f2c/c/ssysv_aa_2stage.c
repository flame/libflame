/* ../netlib/v3.9.0/ssysv_aa_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
/* > \brief <b> SSYSV_AA_2STAGE computes the solution to system of linear equations A * X = B for SY matrices </b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSYSV_AA_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysv_a a_2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysv_a a_2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysv_a a_2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSYSV_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, */
/* IPIV, IPIV2, B, LDB, WORK, LWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER N, NRHS, LDA, LTB, LDB, LWORK, INFO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), IPIV2( * ) */
/* REAL A( LDA, * ), TB( * ), B( LDB, *), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYSV_AA_2STAGE computes the solution to a real system of */
/* > linear equations */
/* > A * X = B, */
/* > where A is an N-by-N symmetric matrix and X and B are N-by-NRHS */
/* > matrices. */
/* > */
/* > Aasen's 2-stage algorithm is used to factor A as */
/* > A = U**T * T * U, if UPLO = 'U', or */
/* > A = L * T * L**T, if UPLO = 'L', */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is symmetric and band. The matrix T is */
/* > then LU-factored with partial pivoting. The factored form of A */
/* > is then used to solve the system of equations A * X = B. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
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
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the symmetric matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, L is stored below (or above) the subdiaonal blocks, */
/* > when UPLO is 'L' (or 'U'). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TB */
/* > \verbatim */
/* > TB is REAL array, dimension (LTB) */
/* > On exit, details of the LU factorization of the band matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LTB */
/* > \verbatim */
/* > LTB is INTEGER */
/* > The size of the array TB. LTB >= 4*N, internally */
/* > used to select NB such that LTB >= (3*NB+1)*N. */
/* > */
/* > If LTB = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal size of LTB, */
/* > returns this value as the first entry of TB, and */
/* > no error message related to LTB is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > On exit, it contains the details of the interchanges, i.e., */
/* > the row and column k of A were interchanged with the */
/* > row and column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV2 */
/* > \verbatim */
/* > IPIV2 is INTEGER array, dimension (N) */
/* > On exit, it contains the details of the interchanges, i.e., */
/* > the row and column k of T were interchanged with the */
/* > row and column IPIV(k). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,NRHS) */
/* > On entry, the right hand side matrix B. */
/* > On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL workspace of size LWORK */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The size of WORK. LWORK >= N, internally used to select NB */
/* > such that LWORK >= N*NB. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal size of the WORK array, */
/* > returns this value as the first entry of the WORK array, and */
/* > no error message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = i, band LU factorization failed on i-th column */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup realSYsolve */
/* ===================================================================== */
/* Subroutine */
int ssysv_aa_2stage_(char *uplo, integer *n, integer *nrhs, real *a, integer *lda, real *tb, integer *ltb, integer *ipiv, integer *ipiv2, real *b, integer *ldb, real *work, integer *lwork, integer * info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"ssysv_aa_2stage inputs: uplo %c, n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS ", ltb %" FLA_IS ", ldb %" FLA_IS "",*uplo, *n, *nrhs, *lda, *ltb, *ldb);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    extern /* Subroutine */
    int ssytrf_aa_2stage_(char *, integer *, real *, integer *, real *, integer *, integer *, integer *, real *, integer *, integer *), ssytrs_aa_2stage_(char *, integer *, integer *, real *, integer *, real *, integer *, integer *, integer *, real *, integer *, integer *);
    extern logical lsame_(char *, char *);
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    integer lwkopt;
    logical tquery, wquery;
    /* -- LAPACK driver routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
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
    --tb;
    --ipiv;
    --ipiv2;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    wquery = *lwork == -1;
    tquery = *ltb == -1;
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
    else if (*lda < fla_max(1,*n))
    {
        *info = -5;
    }
    else if (*ltb < *n << 2 && ! tquery)
    {
        *info = -7;
    }
    else if (*ldb < fla_max(1,*n))
    {
        *info = -11;
    }
    else if (*lwork < *n && ! wquery)
    {
        *info = -13;
    }
    if (*info == 0)
    {
        ssytrf_aa_2stage_(uplo, n, &a[a_offset], lda, &tb[1], &c_n1, &ipiv[1], &ipiv2[1], &work[1], &c_n1, info);
        lwkopt = (integer) work[1];
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSYSV_AA_2STAGE", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (wquery || tquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Compute the factorization A = U**T*T*U or A = L*T*L**T. */
    ssytrf_aa_2stage_(uplo, n, &a[a_offset], lda, &tb[1], ltb, &ipiv[1], & ipiv2[1], &work[1], lwork, info);
    if (*info == 0)
    {
        /* Solve the system A*X = B, overwriting B with X. */
        ssytrs_aa_2stage_(uplo, n, nrhs, &a[a_offset], lda, &tb[1], ltb, & ipiv[1], &ipiv2[1], &b[b_offset], ldb, info);
    }
    work[1] = (real) lwkopt;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of SSYSV_AA_2STAGE */
}
/* ssysv_aa_2stage__ */

