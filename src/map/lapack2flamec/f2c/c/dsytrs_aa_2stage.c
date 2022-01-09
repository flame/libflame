/* ../netlib/v3.9.0/dsytrs_aa_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 static doublereal c_b10 = 1.;
 static integer c_n1 = -1;
 /* > \brief \b DSYTRS_AA_2STAGE */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DSYTRS_AA_2STAGE + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrs_ aa_2stage.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrs_ aa_2stage.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrs_ aa_2stage.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE DSYTRS_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, */
 /* IPIV2, B, LDB, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* INTEGER N, NRHS, LDA, LTB, LDB, INFO */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IPIV( * ), IPIV2( * ) */
 /* DOUBLE PRECISION A( LDA, * ), TB( * ), B( LDB, * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DSYTRS_AA_2STAGE solves a system of linear equations A*X = B with a real */
 /* > symmetric matrix A using the factorization A = U**T*T*U or */
 /* > A = L*T*L**T computed by DSYTRF_AA_2STAGE. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] UPLO */
 /* > \verbatim */
 /* > UPLO is CHARACTER*1 */
 /* > Specifies whether the details of the factorization are stored */
 /* > as an upper or lower triangular matrix. */
 /* > = 'U': Upper triangular, form is A = U**T*T*U;
 */
 /* > = 'L': Lower triangular, form is A = L*T*L**T. */
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
 /* > A is DOUBLE PRECISION array, dimension (LDA,N) */
 /* > Details of factors computed by DSYTRF_AA_2STAGE. */
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
 /* > TB is DOUBLE PRECISION array, dimension (LTB) */
 /* > Details of factors computed by DSYTRF_AA_2STAGE. */
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
 /* > DSYTRF_AA_2STAGE. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] IPIV2 */
 /* > \verbatim */
 /* > IPIV2 is INTEGER array, dimension (N) */
 /* > Details of the interchanges as computed by */
 /* > DSYTRF_AA_2STAGE. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] B */
 /* > \verbatim */
 /* > B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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
 /* > \ingroup doubleSYcomputational */
 /* ===================================================================== */
 /* Subroutine */
 int dsytrs_aa_2stage_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *tb, integer *ltb, integer * ipiv, integer *ipiv2, doublereal *b, integer *ldb, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, b_dim1, b_offset, i__1;
 /* Local variables */
 integer nb, ldtb;
 extern logical lsame_(char *, char *);
 extern /* Subroutine */
 int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
 logical upper;
 extern /* Subroutine */
 int xerbla_(char *, integer *), dgbtrs_( char *, integer *, integer *, integer *, integer *, doublereal *, integer *, integer *, doublereal *, integer *, integer *), dlaswp_(integer *, doublereal *, integer *, integer *, integer *, integer *, integer *);
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
 xerbla_("DSYTRS_AA_2STAGE", &i__1);
 return 0;
 }
 /* Quick return if possible */
 if (*n == 0 || *nrhs == 0) {
 return 0;
 }
 /* Read NB and compute LDTB */
 nb = (integer) tb[1];
 ldtb = *ltb / *n;
 if (upper) {
 /* Solve A*X = B, where A = U**T*T*U. */
 if (*n > nb) {
 /* Pivot, P**T * B -> B */
 i__1 = nb + 1;
 dlaswp_(nrhs, &b[b_offset], ldb, &i__1, n, &ipiv[1], &c__1);
 /* Compute (U**T \ B) -> B [ (U**T \P**T * B) ] */
 i__1 = *n - nb;
 dtrsm_("L", "U", "T", "U", &i__1, nrhs, &c_b10, &a[(nb + 1) * a_dim1 + 1], lda, &b[nb + 1 + b_dim1], ldb);
 }
 /* Compute T \ B -> B [ T \ (U**T \P**T * B) ] */
 dgbtrs_("N", n, &nb, &nb, nrhs, &tb[1], &ldtb, &ipiv2[1], &b[b_offset] , ldb, info);
 if (*n > nb) {
 /* Compute (U \ B) -> B [ U \ (T \ (U**T \P**T * B) ) ] */
 i__1 = *n - nb;
 dtrsm_("L", "U", "N", "U", &i__1, nrhs, &c_b10, &a[(nb + 1) * a_dim1 + 1], lda, &b[nb + 1 + b_dim1], ldb);
 /* Pivot, P * B -> B [ P * (U \ (T \ (U**T \P**T * B) )) ] */
 i__1 = nb + 1;
 dlaswp_(nrhs, &b[b_offset], ldb, &i__1, n, &ipiv[1], &c_n1);
 }
 }
 else {
 /* Solve A*X = B, where A = L*T*L**T. */
 if (*n > nb) {
 /* Pivot, P**T * B -> B */
 i__1 = nb + 1;
 dlaswp_(nrhs, &b[b_offset], ldb, &i__1, n, &ipiv[1], &c__1);
 /* Compute (L \ B) -> B [ (L \P**T * B) ] */
 i__1 = *n - nb;
 dtrsm_("L", "L", "N", "U", &i__1, nrhs, &c_b10, &a[nb + 1 + a_dim1], lda, &b[nb + 1 + b_dim1], ldb);
 }
 /* Compute T \ B -> B [ T \ (L \P**T * B) ] */
 dgbtrs_("N", n, &nb, &nb, nrhs, &tb[1], &ldtb, &ipiv2[1], &b[b_offset] , ldb, info);
 if (*n > nb) {
 /* Compute (L**T \ B) -> B [ L**T \ (T \ (L \P**T * B) ) ] */
 i__1 = *n - nb;
 dtrsm_("L", "L", "T", "U", &i__1, nrhs, &c_b10, &a[nb + 1 + a_dim1], lda, &b[nb + 1 + b_dim1], ldb);
 /* Pivot, P * B -> B [ P * (L**T \ (T \ (L \P**T * B) )) ] */
 i__1 = nb + 1;
 dlaswp_(nrhs, &b[b_offset], ldb, &i__1, n, &ipiv[1], &c_n1);
 }
 }
 return 0;
 /* End of DSYTRS_AA_2STAGE */
 }
 /* dsytrs_aa_2stage__ */
 