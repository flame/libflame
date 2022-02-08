/* ../netlib/v3.9.0/zhesv_rk.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
/* > \brief <b> ZHESV_RK computes the solution to system of linear equations A * X = B for SY matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHESV_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhesv_r k.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhesv_r k.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhesv_r k.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHESV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, */
/* WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, LWORK, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), E( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > ZHESV_RK computes the solution to a complex system of linear */
/* > equations A * X = B, where A is an N-by-N Hermitian matrix */
/* > and X and B are N-by-NRHS matrices. */
/* > */
/* > The bounded Bunch-Kaufman (rook) diagonal pivoting method is used */
/* > to factor A as */
/* > A = P*U*D*(U**H)*(P**T), if UPLO = 'U', or */
/* > A = P*L*D*(L**H)*(P**T), if UPLO = 'L', */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**H (or L**H) is the conjugate of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is Hermitian and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > ZHETRF_RK is called to compute the factorization of a complex */
/* > Hermitian matrix. The factored form of A is then used to solve */
/* > the system of equations A * X = B by calling BLAS3 routine ZHETRS_3. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > Hermitian matrix A is stored: */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. */
/* > If UPLO = 'U': the leading N-by-N upper triangular part */
/* > of A contains the upper triangular part of the matrix A, */
/* > and the strictly lower triangular part of A is not */
/* > referenced. */
/* > */
/* > If UPLO = 'L': the leading N-by-N lower triangular part */
/* > of A contains the lower triangular part of the matrix A, */
/* > and the strictly upper triangular part of A is not */
/* > referenced. */
/* > */
/* > On exit, if INFO = 0, diagonal of the block diagonal */
/* > matrix D and factors U or L as computed by ZHETRF_RK: */
/* > a) ONLY diagonal elements of the Hermitian block diagonal */
/* > matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*/
/* > (superdiagonal (or subdiagonal) elements of D */
/* > are stored on exit in array E), and */
/* > b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* > If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* > For more info see the description of ZHETRF_RK routine. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is COMPLEX*16 array, dimension (N) */
/* > On exit, contains the output computed by the factorization */
/* > routine ZHETRF_RK, i.e. the superdiagonal (or subdiagonal) */
/* > elements of the Hermitian block diagonal matrix D */
/* > with 1-by-1 or 2-by-2 diagonal blocks, where */
/* > If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
*/
/* > If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* > NOTE: For 1-by-1 diagonal block D(k), where */
/* > 1 <= k <= N, the element E(k) is set to 0 in both */
/* > UPLO = 'U' or UPLO = 'L' cases. */
/* > */
/* > For more info see the description of ZHETRF_RK routine. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D, */
/* > as determined by ZHETRF_RK. */
/* > */
/* > For more info see the description of ZHETRF_RK routine. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the N-by-NRHS right hand side matrix B. */
/* > On exit, if INFO = 0, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension ( MAX(1,LWORK) ). */
/* > Work array used in the factorization stage. */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of WORK. LWORK >= 1. For best performance */
/* > of factorization stage LWORK >= max(1,N*NB), where NB is */
/* > the optimal blocksize for ZHETRF_RK. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
*/
/* > the routine only calculates the optimal size of the WORK */
/* > array for factorization stage, returns this value as */
/* > the first entry of the WORK array, and no error message */
/* > related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > */
/* > < 0: If INFO = -k, the k-th argument had an illegal value */
/* > */
/* > > 0: If INFO = k, the matrix A is singular, because: */
/* > If UPLO = 'U': column k in the upper */
/* > triangular part of A contains all zeros. */
/* > If UPLO = 'L': column k in the lower */
/* > triangular part of A contains all zeros. */
/* > */
/* > Therefore D(k,k) is exactly zero, and superdiagonal */
/* > elements of column k of U (or subdiagonal elements of */
/* > column k of L ) are all zeros. The factorization has */
/* > been completed, but the block diagonal matrix D is */
/* > exactly singular, and division by zero will occur if */
/* > it is used to solve a system of equations. */
/* > */
/* > NOTE: INFO only stores the first occurrence of */
/* > a singularity, any subsequent occurrence of singularity */
/* > is not stored in INFO even though the factorization */
/* > always completes. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup complex16HEsolve */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > December 2016, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* > School of Mathematics, */
/* > University of Manchester */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int zhesv_rk_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"zhesv_rk inputs: uplo %c, n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS ", lwork %" FLA_IS "",*uplo, *n, *nrhs, *lda, *ldb, *lwork);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    extern /* Subroutine */
    int zhetrs_3_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *), zhetrf_rk_(char *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    integer lwkopt;
    logical lquery;
    /* -- LAPACK driver routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
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
    --e;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
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
        *info = -9;
    }
    else if (*lwork < 1 && ! lquery)
    {
        *info = -11;
    }
    if (*info == 0)
    {
        if (*n == 0)
        {
            lwkopt = 1;
        }
        else
        {
            zhetrf_rk_(uplo, n, &a[a_offset], lda, &e[1], &ipiv[1], &work[1], &c_n1, info);
            lwkopt = (integer) work[1].r;
        }
        work[1].r = (doublereal) lwkopt;
        work[1].i = 0.; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHESV_RK ", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Compute the factorization A = P*U*D*(U**H)*(P**T) or */
    /* A = P*U*D*(U**H)*(P**T). */
    zhetrf_rk_(uplo, n, &a[a_offset], lda, &e[1], &ipiv[1], &work[1], lwork, info);
    if (*info == 0)
    {
        /* Solve the system A*X = B with BLAS3 solver, overwriting B with X. */
        zhetrs_3_(uplo, n, nrhs, &a[a_offset], lda, &e[1], &ipiv[1], &b[ b_offset], ldb, info);
    }
    work[1].r = (doublereal) lwkopt;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of ZHESV_RK */
}
/* zhesv_rk__ */

