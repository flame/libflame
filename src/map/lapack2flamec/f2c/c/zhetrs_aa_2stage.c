/* ../netlib/v3.9.0/zhetrs_aa_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static doublecomplex c_b1 = {
1.,0.}
;
 static integer c__1 = 1;
 static integer c_n1 = -1;
 /* > \brief \b ZHETRS_AA_2STAGE */
 /* @generated from SRC/dsytrs_aa_2stage.f, fortran d -> c, Mon Oct 30 11:59:02 2017 */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download ZHETRS_AA_2STAGE + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs_ aa_2stage.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs_ aa_2stage.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs_ aa_2stage.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE ZHETRS_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, */
 /* IPIV2, B, LDB, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* INTEGER N, NRHS, LDA, LTB, LDB, INFO */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IPIV( * ), IPIV2( * ) */
 /* COMPLEX*16 A( LDA, * ), TB( * ), B( LDB, * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > ZHETRS_AA_2STAGE solves a system of linear equations A*X = B with a */
 /* > hermitian matrix A using the factorization A = U**H*T*U or */
 /* > A = L*T*L**H computed by ZHETRF_AA_2STAGE. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] UPLO */
 /* > \verbatim */
 /* > UPLO is CHARACTER*1 */
 /* > Specifies whether the details of the factorization are stored */
 /* > as an upper or lower triangular matrix. */
 /* > = 'U': Upper triangular, form is A = U**H*T*U;
 */
 /* > = 'L': Lower triangular, form is A = L*T*L**H. */
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
 /* > \param[in] A */
 /* > \verbatim */
 /* > A is COMPLEX*16 array, dimension (LDA,N) */
 /* > Details of factors computed by ZHETRF_AA_2STAGE. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,N). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] TB */
 /* > \verbatim */
 /* > TB is COMPLEX*16 array, dimension (LTB) */
 /* > Details of factors computed by ZHETRF_AA_2STAGE. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LTB */
 /* > \verbatim */
 /* > LTB is INTEGER */
 /* > The size of the array TB. LTB >= 4*N. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] IPIV */
 /* > \verbatim */
 /* > IPIV is INTEGER array, dimension (N) */
 /* > Details of the interchanges as computed by */
 /* > ZHETRF_AA_2STAGE. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] IPIV2 */
 /* > \verbatim */
 /* > IPIV2 is INTEGER array, dimension (N) */
 /* > Details of the interchanges as computed by */
 /* > ZHETRF_AA_2STAGE. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] B */
 /* > \verbatim */
 /* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
 /* > On entry, the right hand side matrix B. */
 /* > On exit, the solution matrix X. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDB */
 /* > \verbatim */
 /* > LDB is INTEGER */
 /* > The leading dimension of the array B. LDB >= max(1,N). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2017 */
 /* > \ingroup complex16SYcomputational */
 /* ===================================================================== */
 /* Subroutine */
 int zhetrs_aa_2stage_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *tb, integer *ltb, integer *ipiv, integer *ipiv2, doublecomplex *b, integer *ldb, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, b_dim1, b_offset, i__1;
 /* Local variables */
 integer nb, ldtb;
 extern logical lsame_(char *, char *);
 logical upper;
 extern /* Subroutine */
 int ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(char *, integer *), zgbtrs_(char *, integer *, integer *, integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *, integer *), zlaswp_(integer *, doublecomplex *, integer *, integer *, integer *, integer *, integer *);
 /* -- LAPACK computational routine (version 3.8.0) -- */
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
 /* Function Body */
 *info = 0;
 upper = lsame_(uplo, "U");
 if (! upper && ! lsame_(uplo, "L")) {
 *info = -1;
 }
 else if (*n < 0) {
 *info = -2;
 }
 else if (*nrhs < 0) {
 *info = -3;
 }
 else if (*lda < max(1,*n)) {
 *info = -5;
 }
 else if (*ltb < *n << 2) {
 *info = -7;
 }
 else if (*ldb < max(1,*n)) {
 *info = -11;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("ZHETRS_AA_2STAGE", &i__1);
 return 0;
 }
 /* Quick return if possible */
 if (*n == 0 || *nrhs == 0) {
 return 0;
 }
 /* Read NB and compute LDTB */
 nb = (integer) tb[1].r;
 ldtb = *ltb / *n;
 if (upper) {
 /* Solve A*X = B, where A = U**H*T*U. */
 if (*n > nb) {
 /* Pivot, P**T * B -> B */
 i__1 = nb + 1;
 zlaswp_(nrhs, &b[b_offset], ldb, &i__1, n, &ipiv[1], &c__1);
 /* Compute (U**H \ B) -> B [ (U**H \P**T * B) ] */
 i__1 = *n - nb;
 ztrsm_("L", "U", "C", "U", &i__1, nrhs, &c_b1, &a[(nb + 1) * a_dim1 + 1], lda, &b[nb + 1 + b_dim1], ldb);
 }
 /* Compute T \ B -> B [ T \ (U**H \P**T * B) ] */
 zgbtrs_("N", n, &nb, &nb, nrhs, &tb[1], &ldtb, &ipiv2[1], &b[b_offset] , ldb, info);
 if (*n > nb) {
 /* Compute (U \ B) -> B [ U \ (T \ (U**H \P**T * B) ) ] */
 i__1 = *n - nb;
 ztrsm_("L", "U", "N", "U", &i__1, nrhs, &c_b1, &a[(nb + 1) * a_dim1 + 1], lda, &b[nb + 1 + b_dim1], ldb);
 /* Pivot, P * B -> B [ P * (U \ (T \ (U**H \P**T * B) )) ] */
 i__1 = nb + 1;
 zlaswp_(nrhs, &b[b_offset], ldb, &i__1, n, &ipiv[1], &c_n1);
 }
 }
 else {
 /* Solve A*X = B, where A = L*T*L**H. */
 if (*n > nb) {
 /* Pivot, P**T * B -> B */
 i__1 = nb + 1;
 zlaswp_(nrhs, &b[b_offset], ldb, &i__1, n, &ipiv[1], &c__1);
 /* Compute (L \ B) -> B [ (L \P**T * B) ] */
 i__1 = *n - nb;
 ztrsm_("L", "L", "N", "U", &i__1, nrhs, &c_b1, &a[nb + 1 + a_dim1] , lda, &b[nb + 1 + b_dim1], ldb);
 }
 /* Compute T \ B -> B [ T \ (L \P**T * B) ] */
 zgbtrs_("N", n, &nb, &nb, nrhs, &tb[1], &ldtb, &ipiv2[1], &b[b_offset] , ldb, info);
 if (*n > nb) {
 /* Compute (L**H \ B) -> B [ L**H \ (T \ (L \P**T * B) ) ] */
 i__1 = *n - nb;
 ztrsm_("L", "L", "C", "U", &i__1, nrhs, &c_b1, &a[nb + 1 + a_dim1] , lda, &b[nb + 1 + b_dim1], ldb);
 /* Pivot, P * B -> B [ P * (L**H \ (T \ (L \P**T * B) )) ] */
 i__1 = nb + 1;
 zlaswp_(nrhs, &b[b_offset], ldb, &i__1, n, &ipiv[1], &c_n1);
 }
 }
 return 0;
 /* End of ZHETRS_AA_2STAGE */
 }
 /* zhetrs_aa_2stage__ */
 