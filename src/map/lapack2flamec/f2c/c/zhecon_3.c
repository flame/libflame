/* ../netlib/v3.9.0/zhecon_3.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 /* > \brief \b ZHECON_3 */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download ZHECON_3 + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhecon_ 3.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhecon_ 3.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhecon_ 3.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE ZHECON_3( UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, */
 /* WORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* INTEGER INFO, LDA, N */
 /* DOUBLE PRECISION ANORM, RCOND */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IPIV( * ) */
 /* COMPLEX*16 A( LDA, * ), E ( * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > ZHECON_3 estimates the reciprocal of the condition number (in the */
 /* > 1-norm) of a complex Hermitian matrix A using the factorization */
 /* > computed by ZHETRF_RK or ZHETRF_BK: */
 /* > */
 /* > A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T), */
 /* > */
 /* > where U (or L) is unit upper (or lower) triangular matrix, */
 /* > U**H (or L**H) is the conjugate of U (or L), P is a permutation */
 /* > matrix, P**T is the transpose of P, and D is Hermitian and block */
 /* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
 /* > */
 /* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
 /* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
 /* > This routine uses BLAS3 solver ZHETRS_3. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] UPLO */
 /* > \verbatim */
 /* > UPLO is CHARACTER*1 */
 /* > Specifies whether the details of the factorization are */
 /* > stored as an upper or lower triangular matrix: */
 /* > = 'U': Upper triangular, form is A = P*U*D*(U**H)*(P**T);
 */
 /* > = 'L': Lower triangular, form is A = P*L*D*(L**H)*(P**T). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The order of the matrix A. N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] A */
 /* > \verbatim */
 /* > A is COMPLEX*16 array, dimension (LDA,N) */
 /* > Diagonal of the block diagonal matrix D and factors U or L */
 /* > as computed by ZHETRF_RK and ZHETRF_BK: */
 /* > a) ONLY diagonal elements of the Hermitian block diagonal */
 /* > matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
 */
 /* > (superdiagonal (or subdiagonal) elements of D */
 /* > should be provided on entry in array E), and */
 /* > b) If UPLO = 'U': factor U in the superdiagonal part of A. */
 /* > If UPLO = 'L': factor L in the subdiagonal part of A. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,N). */
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
 /* > \param[in] ANORM */
 /* > \verbatim */
 /* > ANORM is DOUBLE PRECISION */
 /* > The 1-norm of the original matrix A. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] RCOND */
 /* > \verbatim */
 /* > RCOND is DOUBLE PRECISION */
 /* > The reciprocal of the condition number of the matrix A, */
 /* > computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
 /* > estimate of the 1-norm of inv(A) computed in this routine. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is COMPLEX*16 array, dimension (2*N) */
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
 /* > \date June 2017 */
 /* > \ingroup complex16HEcomputational */
 /* > \par Contributors: */
 /* ================== */
 /* > \verbatim */
 /* > */
 /* > June 2017, Igor Kozachenko, */
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
 int zhecon_3_(char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *work, integer *info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
 char buffer[256];
 snprintf(buffer, 256,"zhecon_3 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS ", ipiv %" FLA_IS "",*uplo, *n, *lda, *ipiv);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, i__1, i__2;
 /* Local variables */
 extern /* Subroutine */
 int zhetrs_3_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
 integer i__, kase;
 extern logical lsame_(char *, char *);
 integer isave[3];
 logical upper;
 extern /* Subroutine */
 int zlacn2_(integer *, doublecomplex *, doublecomplex *, doublereal *, integer *, integer *), xerbla_( char *, integer *);
 doublereal ainvnm;
 /* -- LAPACK computational routine (version 3.7.1) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* June 2017 */
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
 if (! upper && ! lsame_(uplo, "L")) {
 *info = -1;
 }
 else if (*n < 0) {
 *info = -2;
 }
 else if (*lda < max(1,*n)) {
 *info = -4;
 }
 else if (*anorm < 0.) {
 *info = -7;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("ZHECON_3", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Quick return if possible */
 *rcond = 0.;
 if (*n == 0) {
 *rcond = 1.;
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 else if (*anorm <= 0.) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Check that the diagonal matrix D is nonsingular. */
 if (upper) {
 /* Upper triangular storage: examine D from bottom to top */
 for (i__ = *n;
 i__ >= 1;
 --i__) {
 i__1 = i__ + i__ * a_dim1;
 if (ipiv[i__] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.)) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 }
 }
 else {
 /* Lower triangular storage: examine D from top to bottom. */
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = i__ + i__ * a_dim1;
 if (ipiv[i__] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.)) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 }
 }
 /* Estimate the 1-norm of the inverse. */
 kase = 0;
 L30: zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
 if (kase != 0) {
 /* Multiply by inv(L*D*L**H) or inv(U*D*U**H). */
 zhetrs_3_(uplo, n, &c__1, &a[a_offset], lda, &e[1], &ipiv[1], &work[ 1], n, info);
 goto L30;
 }
 /* Compute the estimate of the reciprocal condition number. */
 if (ainvnm != 0.) {
 *rcond = 1. / ainvnm / *anorm;
 }
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of ZHECON_3 */
 }
 /* zhecon_3__ */
 
