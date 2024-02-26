/* ../netlib/v3.9.0/zhetri_3.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b ZHETRI_3 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHETRI_3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri_ 3.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri_ 3.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri_ 3.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHETRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), E( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > ZHETRI_3 computes the inverse of a complex Hermitian indefinite */
/* > matrix A using the factorization computed by ZHETRF_RK or ZHETRF_BK: */
/* > */
/* > A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**H (or L**H) is the conjugate of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is Hermitian and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > ZHETRI_3 sets the leading dimension of the workspace before calling */
/* > ZHETRI_3X that actually computes the inverse. This is the blocked */
/* > version of the algorithm, calling Level 3 BLAS. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are */
/* > stored as an upper or lower triangular matrix. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, diagonal of the block diagonal matrix D and */
/* > factors U or L as computed by ZHETRF_RK and ZHETRF_BK: */
/* > a) ONLY diagonal elements of the Hermitian block diagonal */
/* > matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*/
/* > (superdiagonal (or subdiagonal) elements of D */
/* > should be provided on entry in array E), and */
/* > b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* > If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* > On exit, if INFO = 0, the Hermitian inverse of the original */
/* > matrix. */
/* > If UPLO = 'U': the upper triangular part of the inverse */
/* > is formed and the part of A below the diagonal is not */
/* > referenced;
*/
/* > If UPLO = 'L': the lower triangular part of the inverse */
/* > is formed and the part of A above the diagonal is not */
/* > referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is COMPLEX*16 array, dimension (N) */
/* > On entry, contains the superdiagonal (or subdiagonal) */
/* > elements of the Hermitian block diagonal matrix D */
/* > with 1-by-1 or 2-by-2 diagonal blocks, where */
/* > If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*/
/* > If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* > NOTE: For 1-by-1 diagonal block D(k), where */
/* > 1 <= k <= N, the element E(k) is not referenced in both */
/* > UPLO = 'U' or UPLO = 'L' cases. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by ZHETRF_RK or ZHETRF_BK. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N+NB+1)*(NB+3). */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of WORK. LWORK >= (N+NB+1)*(NB+3). */
/* > */
/* > If LDWORK = -1, then a workspace query is assumed;
*/
/* > the routine only calculates the optimal size of the optimal */
/* > size of the WORK array, returns this value as the first */
/* > entry of the WORK array, and no error message related to */
/* > LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, D(i,i) = 0;
the matrix is singular and its */
/* > inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup complex16HEcomputational */
/* > \par Contributors: */
/* ================== */
/* > \verbatim */
/* > */
/* > November 2017, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int zhetri_3_(char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, doublecomplex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhetri_3 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS ", lwork %" FLA_IS "",*uplo, *n, *lda, *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    extern /* Subroutine */
    int zhetri_3x_(char *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    integer nb;
    extern logical lsame_(char *, char *);
    logical upper;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
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
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1;
    /* Determine the block size */
    /* Computing MAX */
    i__1 = 1;
    i__2 = ilaenv_(&c__1, "ZHETRI_3", uplo, n, &c_n1, &c_n1, &c_n1); // , expr subst
    nb = fla_max(i__1,i__2);
    lwkopt = (*n + nb + 1) * (nb + 3);
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*n))
    {
        *info = -4;
    }
    else if (*lwork < lwkopt && ! lquery)
    {
        *info = -8;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHETRI_3", &i__1, (ftnlen)8);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (lquery)
    {
        work[1].r = (doublereal) lwkopt;
        work[1].i = 0.; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    zhetri_3x_(uplo, n, &a[a_offset], lda, &e[1], &ipiv[1], &work[1], &nb, info);
    work[1].r = (doublereal) lwkopt;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZHETRI_3 */
}
/* zhetri_3__ */

