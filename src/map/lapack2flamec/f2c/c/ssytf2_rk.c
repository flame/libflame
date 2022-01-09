/* ../netlib/v3.9.0/ssytf2_rk.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 /* > \brief \b SSYTF2_RK computes the factorization of a real symmetric indefinite matrix using the bounded Bu nch-Kaufman (rook) diagonal pivoting method (BLAS2 unblocked algorithm). */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download SSYTF2_RK + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytf2_ rk.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytf2_ rk.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytf2_ rk.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE SSYTF2_RK( UPLO, N, A, LDA, E, IPIV, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* INTEGER INFO, LDA, N */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IPIV( * ) */
 /* REAL A( LDA, * ), E ( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > SSYTF2_RK computes the factorization of a real symmetric matrix A */
 /* > using the bounded Bunch-Kaufman (rook) diagonal pivoting method: */
 /* > */
 /* > A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T), */
 /* > */
 /* > where U (or L) is unit upper (or lower) triangular matrix, */
 /* > U**T (or L**T) is the transpose of U (or L), P is a permutation */
 /* > matrix, P**T is the transpose of P, and D is symmetric and block */
 /* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
 /* > */
 /* > This is the unblocked version of the algorithm, calling Level 2 BLAS. */
 /* > For more information see Further Details section. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] UPLO */
 /* > \verbatim */
 /* > UPLO is CHARACTER*1 */
 /* > Specifies whether the upper or lower triangular part of the */
 /* > symmetric matrix A is stored: */
 /* > = 'U': Upper triangular */
 /* > = 'L': Lower triangular */
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
 /* > A is REAL array, dimension (LDA,N) */
 /* > On entry, the symmetric matrix A. */
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
 /* > On exit, contains: */
 /* > a) ONLY diagonal elements of the symmetric block diagonal */
 /* > matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
 */
 /* > (superdiagonal (or subdiagonal) elements of D */
 /* > are stored on exit in array E), and */
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
 /* > \param[out] E */
 /* > \verbatim */
 /* > E is REAL array, dimension (N) */
 /* > On exit, contains the superdiagonal (or subdiagonal) */
 /* > elements of the symmetric block diagonal matrix D */
 /* > with 1-by-1 or 2-by-2 diagonal blocks, where */
 /* > If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
 */
 /* > If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
 /* > */
 /* > NOTE: For 1-by-1 diagonal block D(k), where */
 /* > 1 <= k <= N, the element E(k) is set to 0 in both */
 /* > UPLO = 'U' or UPLO = 'L' cases. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IPIV */
 /* > \verbatim */
 /* > IPIV is INTEGER array, dimension (N) */
 /* > IPIV describes the permutation matrix P in the factorization */
 /* > of matrix A as follows. The absolute value of IPIV(k) */
 /* > represents the index of row and column that were */
 /* > interchanged with the k-th row and column. The value of UPLO */
 /* > describes the order in which the interchanges were applied. */
 /* > Also, the sign of IPIV represents the block structure of */
 /* > the symmetric block diagonal matrix D with 1-by-1 or 2-by-2 */
 /* > diagonal blocks which correspond to 1 or 2 interchanges */
 /* > at each factorization step. For more info see Further */
 /* > Details section. */
 /* > */
 /* > If UPLO = 'U', */
 /* > ( in factorization order, k decreases from N to 1 ): */
 /* > a) A single positive entry IPIV(k) > 0 means: */
 /* > D(k,k) is a 1-by-1 diagonal block. */
 /* > If IPIV(k) != k, rows and columns k and IPIV(k) were */
 /* > interchanged in the matrix A(1:N,1:N);
 */
 /* > If IPIV(k) = k, no interchange occurred. */
 /* > */
 /* > b) A pair of consecutive negative entries */
 /* > IPIV(k) < 0 and IPIV(k-1) < 0 means: */
 /* > D(k-1:k,k-1:k) is a 2-by-2 diagonal block. */
 /* > (NOTE: negative entries in IPIV appear ONLY in pairs). */
 /* > 1) If -IPIV(k) != k, rows and columns */
 /* > k and -IPIV(k) were interchanged */
 /* > in the matrix A(1:N,1:N). */
 /* > If -IPIV(k) = k, no interchange occurred. */
 /* > 2) If -IPIV(k-1) != k-1, rows and columns */
 /* > k-1 and -IPIV(k-1) were interchanged */
 /* > in the matrix A(1:N,1:N). */
 /* > If -IPIV(k-1) = k-1, no interchange occurred. */
 /* > */
 /* > c) In both cases a) and b), always ABS( IPIV(k) ) <= k. */
 /* > */
 /* > d) NOTE: Any entry IPIV(k) is always NONZERO on output. */
 /* > */
 /* > If UPLO = 'L', */
 /* > ( in factorization order, k increases from 1 to N ): */
 /* > a) A single positive entry IPIV(k) > 0 means: */
 /* > D(k,k) is a 1-by-1 diagonal block. */
 /* > If IPIV(k) != k, rows and columns k and IPIV(k) were */
 /* > interchanged in the matrix A(1:N,1:N). */
 /* > If IPIV(k) = k, no interchange occurred. */
 /* > */
 /* > b) A pair of consecutive negative entries */
 /* > IPIV(k) < 0 and IPIV(k+1) < 0 means: */
 /* > D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
 /* > (NOTE: negative entries in IPIV appear ONLY in pairs). */
 /* > 1) If -IPIV(k) != k, rows and columns */
 /* > k and -IPIV(k) were interchanged */
 /* > in the matrix A(1:N,1:N). */
 /* > If -IPIV(k) = k, no interchange occurred. */
 /* > 2) If -IPIV(k+1) != k+1, rows and columns */
 /* > k-1 and -IPIV(k-1) were interchanged */
 /* > in the matrix A(1:N,1:N). */
 /* > If -IPIV(k+1) = k+1, no interchange occurred. */
 /* > */
 /* > c) In both cases a) and b), always ABS( IPIV(k) ) >= k. */
 /* > */
 /* > d) NOTE: Any entry IPIV(k) is always NONZERO on output. */
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
 /* > \ingroup singleSYcomputational */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > TODO: put further details */
 /* > \endverbatim */
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
 /* > 01-01-96 - Based on modifications by */
 /* > J. Lewis, Boeing Computer Services Company */
 /* > A. Petitet, Computer Science Dept., */
 /* > Univ. of Tenn., Knoxville abd , USA */
 /* > \endverbatim */
 /* ===================================================================== */
 /* Subroutine */
 int ssytf2_rk_(char *uplo, integer *n, real *a, integer * lda, real *e, integer *ipiv, integer *info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
 char buffer[256];
 snprintf(buffer, 256,"ssytf2_rk inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "",*uplo, *n, *lda);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, i__1, i__2;
 real r__1;
 /* Builtin functions */
 double sqrt(doublereal);
 /* Local variables */
 integer i__, j, k, p;
 real t, d11, d12, d21, d22;
 integer ii, kk, kp;
 real wk, wkm1, wkp1;
 logical done;
 integer imax, jmax;
 extern /* Subroutine */
 int ssyr_(char *, integer *, real *, real *, integer *, real *, integer *);
 real alpha;
 extern logical lsame_(char *, char *);
 extern /* Subroutine */
 int sscal_(integer *, real *, real *, integer *);
 real sfmin;
 integer itemp, kstep;
 real stemp;
 logical upper;
 extern /* Subroutine */
 int sswap_(integer *, real *, integer *, real *, integer *);
 real absakk;
 extern real slamch_(char *);
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 extern integer isamax_(integer *, real *, integer *);
 real colmax, rowmax;
 /* -- LAPACK computational routine (version 3.7.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* December 2016 */
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
 --e;
 --ipiv;
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
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("SSYTF2_RK", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Initialize ALPHA for use in choosing pivot block size. */
 alpha = (sqrt(17.f) + 1.f) / 8.f;
 /* Compute machine safe minimum */
 sfmin = slamch_("S");
 if (upper) {
 /* Factorize A as U*D*U**T using the upper triangle of A */
 /* Initialize the first entry of array E, where superdiagonal */
 /* elements of D are stored */
 e[1] = 0.f;
 /* K is the main loop index, decreasing from N to 1 in steps of */
 /* 1 or 2 */
 k = *n;
 L10: /* If K < 1, exit from loop */
 if (k < 1) {
 goto L34;
 }
 kstep = 1;
 p = k;
 /* Determine rows and columns to be interchanged and whether */
 /* a 1-by-1 or 2-by-2 pivot block will be used */
 absakk = (r__1 = a[k + k * a_dim1], f2c_abs(r__1));
 /* IMAX is the row-index of the largest off-diagonal element in */
 /* column K, and COLMAX is its absolute value. */
 /* Determine both COLMAX and IMAX. */
 if (k > 1) {
 i__1 = k - 1;
 imax = isamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
 colmax = (r__1 = a[imax + k * a_dim1], f2c_abs(r__1));
 }
 else {
 colmax = 0.f;
 }
 if (max(absakk,colmax) == 0.f) {
 /* Column K is zero or underflow: set INFO and continue */
 if (*info == 0) {
 *info = k;
 }
 kp = k;
 /* Set E( K ) to zero */
 if (k > 1) {
 e[k] = 0.f;
 }
 }
 else {
 /* Test for interchange */
 /* Equivalent to testing for (used to handle NaN and Inf) */
 /* ABSAKK.GE.ALPHA*COLMAX */
 if (! (absakk < alpha * colmax)) {
 /* no interchange, */
 /* use 1-by-1 pivot block */
 kp = k;
 }
 else {
 done = FALSE_;
 /* Loop until pivot found */
 L12: /* Begin pivot search loop body */
 /* JMAX is the column-index of the largest off-diagonal */
 /* element in row IMAX, and ROWMAX is its absolute value. */
 /* Determine both ROWMAX and JMAX. */
 if (imax != k) {
 i__1 = k - imax;
 jmax = imax + isamax_(&i__1, &a[imax + (imax + 1) * a_dim1], lda);
 rowmax = (r__1 = a[imax + jmax * a_dim1], f2c_abs(r__1));
 }
 else {
 rowmax = 0.f;
 }
 if (imax > 1) {
 i__1 = imax - 1;
 itemp = isamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
 stemp = (r__1 = a[itemp + imax * a_dim1], f2c_abs(r__1));
 if (stemp > rowmax) {
 rowmax = stemp;
 jmax = itemp;
 }
 }
 /* Equivalent to testing for (used to handle NaN and Inf) */
 /* ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */
 if (! ((r__1 = a[imax + imax * a_dim1], f2c_abs(r__1)) < alpha * rowmax)) {
 /* interchange rows and columns K and IMAX, */
 /* use 1-by-1 pivot block */
 kp = imax;
 done = TRUE_;
 /* Equivalent to testing for ROWMAX .EQ. COLMAX, */
 /* used to handle NaN and Inf */
 }
 else if (p == jmax || rowmax <= colmax) {
 /* interchange rows and columns K+1 and IMAX, */
 /* use 2-by-2 pivot block */
 kp = imax;
 kstep = 2;
 done = TRUE_;
 }
 else {
 /* Pivot NOT found, set variables and repeat */
 p = imax;
 colmax = rowmax;
 imax = jmax;
 }
 /* End pivot search loop body */
 if (! done) {
 goto L12;
 }
 }
 /* Swap TWO rows and TWO columns */
 /* First swap */
 if (kstep == 2 && p != k) {
 /* Interchange rows and column K and P in the leading */
 /* submatrix A(1:k,1:k) if we have a 2-by-2 pivot */
 if (p > 1) {
 i__1 = p - 1;
 sswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &c__1);
 }
 if (p < k - 1) {
 i__1 = k - p - 1;
 sswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 1) * a_dim1], lda);
 }
 t = a[k + k * a_dim1];
 a[k + k * a_dim1] = a[p + p * a_dim1];
 a[p + p * a_dim1] = t;
 /* Convert upper triangle of A into U form by applying */
 /* the interchanges in columns k+1:N. */
 if (k < *n) {
 i__1 = *n - k;
 sswap_(&i__1, &a[k + (k + 1) * a_dim1], lda, &a[p + (k + 1) * a_dim1], lda);
 }
 }
 /* Second swap */
 kk = k - kstep + 1;
 if (kp != kk) {
 /* Interchange rows and columns KK and KP in the leading */
 /* submatrix A(1:k,1:k) */
 if (kp > 1) {
 i__1 = kp - 1;
 sswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
 }
 if (kk > 1 && kp < kk - 1) {
 i__1 = kk - kp - 1;
 sswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + ( kp + 1) * a_dim1], lda);
 }
 t = a[kk + kk * a_dim1];
 a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
 a[kp + kp * a_dim1] = t;
 if (kstep == 2) {
 t = a[k - 1 + k * a_dim1];
 a[k - 1 + k * a_dim1] = a[kp + k * a_dim1];
 a[kp + k * a_dim1] = t;
 }
 /* Convert upper triangle of A into U form by applying */
 /* the interchanges in columns k+1:N. */
 if (k < *n) {
 i__1 = *n - k;
 sswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k + 1) * a_dim1], lda);
 }
 }
 /* Update the leading submatrix */
 if (kstep == 1) {
 /* 1-by-1 pivot block D(k): column k now holds */
 /* W(k) = U(k)*D(k) */
 /* where U(k) is the k-th column of U */
 if (k > 1) {
 /* Perform a rank-1 update of A(1:k-1,1:k-1) and */
 /* store U(k) in column k */
 if ((r__1 = a[k + k * a_dim1], f2c_abs(r__1)) >= sfmin) {
 /* Perform a rank-1 update of A(1:k-1,1:k-1) as */
 /* A := A - U(k)*D(k)*U(k)**T */
 /* = A - W(k)*1/D(k)*W(k)**T */
 d11 = 1.f / a[k + k * a_dim1];
 i__1 = k - 1;
 r__1 = -d11;
 ssyr_(uplo, &i__1, &r__1, &a[k * a_dim1 + 1], &c__1, & a[a_offset], lda);
 /* Store U(k) in column k */
 i__1 = k - 1;
 sscal_(&i__1, &d11, &a[k * a_dim1 + 1], &c__1);
 }
 else {
 /* Store L(k) in column K */
 d11 = a[k + k * a_dim1];
 i__1 = k - 1;
 for (ii = 1;
 ii <= i__1;
 ++ii) {
 a[ii + k * a_dim1] /= d11;
 /* L16: */
 }
 /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
 /* A := A - U(k)*D(k)*U(k)**T */
 /* = A - W(k)*(1/D(k))*W(k)**T */
 /* = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */
 i__1 = k - 1;
 r__1 = -d11;
 ssyr_(uplo, &i__1, &r__1, &a[k * a_dim1 + 1], &c__1, & a[a_offset], lda);
 }
 /* Store the superdiagonal element of D in array E */
 e[k] = 0.f;
 }
 }
 else {
 /* 2-by-2 pivot block D(k): columns k and k-1 now hold */
 /* ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */
 /* where U(k) and U(k-1) are the k-th and (k-1)-th columns */
 /* of U */
 /* Perform a rank-2 update of A(1:k-2,1:k-2) as */
 /* A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
 /* = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T */
 /* and store L(k) and L(k+1) in columns k and k+1 */
 if (k > 2) {
 d12 = a[k - 1 + k * a_dim1];
 d22 = a[k - 1 + (k - 1) * a_dim1] / d12;
 d11 = a[k + k * a_dim1] / d12;
 t = 1.f / (d11 * d22 - 1.f);
 for (j = k - 2;
 j >= 1;
 --j) {
 wkm1 = t * (d11 * a[j + (k - 1) * a_dim1] - a[j + k * a_dim1]);
 wk = t * (d22 * a[j + k * a_dim1] - a[j + (k - 1) * a_dim1]);
 for (i__ = j;
 i__ >= 1;
 --i__) {
 a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ + k * a_dim1] / d12 * wk - a[i__ + (k - 1) * a_dim1] / d12 * wkm1;
 /* L20: */
 }
 /* Store U(k) and U(k-1) in cols k and k-1 for row J */
 a[j + k * a_dim1] = wk / d12;
 a[j + (k - 1) * a_dim1] = wkm1 / d12;
 /* L30: */
 }
 }
 /* Copy superdiagonal elements of D(K) to E(K) and */
 /* ZERO out superdiagonal entry of A */
 e[k] = a[k - 1 + k * a_dim1];
 e[k - 1] = 0.f;
 a[k - 1 + k * a_dim1] = 0.f;
 }
 /* End column K is nonsingular */
 }
 /* Store details of the interchanges in IPIV */
 if (kstep == 1) {
 ipiv[k] = kp;
 }
 else {
 ipiv[k] = -p;
 ipiv[k - 1] = -kp;
 }
 /* Decrease K and return to the start of the main loop */
 k -= kstep;
 goto L10;
 L34: ;
 }
 else {
 /* Factorize A as L*D*L**T using the lower triangle of A */
 /* Initialize the unused last entry of the subdiagonal array E. */
 e[*n] = 0.f;
 /* K is the main loop index, increasing from 1 to N in steps of */
 /* 1 or 2 */
 k = 1;
 L40: /* If K > N, exit from loop */
 if (k > *n) {
 goto L64;
 }
 kstep = 1;
 p = k;
 /* Determine rows and columns to be interchanged and whether */
 /* a 1-by-1 or 2-by-2 pivot block will be used */
 absakk = (r__1 = a[k + k * a_dim1], f2c_abs(r__1));
 /* IMAX is the row-index of the largest off-diagonal element in */
 /* column K, and COLMAX is its absolute value. */
 /* Determine both COLMAX and IMAX. */
 if (k < *n) {
 i__1 = *n - k;
 imax = k + isamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
 colmax = (r__1 = a[imax + k * a_dim1], f2c_abs(r__1));
 }
 else {
 colmax = 0.f;
 }
 if (max(absakk,colmax) == 0.f) {
 /* Column K is zero or underflow: set INFO and continue */
 if (*info == 0) {
 *info = k;
 }
 kp = k;
 /* Set E( K ) to zero */
 if (k < *n) {
 e[k] = 0.f;
 }
 }
 else {
 /* Test for interchange */
 /* Equivalent to testing for (used to handle NaN and Inf) */
 /* ABSAKK.GE.ALPHA*COLMAX */
 if (! (absakk < alpha * colmax)) {
 /* no interchange, use 1-by-1 pivot block */
 kp = k;
 }
 else {
 done = FALSE_;
 /* Loop until pivot found */
 L42: /* Begin pivot search loop body */
 /* JMAX is the column-index of the largest off-diagonal */
 /* element in row IMAX, and ROWMAX is its absolute value. */
 /* Determine both ROWMAX and JMAX. */
 if (imax != k) {
 i__1 = imax - k;
 jmax = k - 1 + isamax_(&i__1, &a[imax + k * a_dim1], lda);
 rowmax = (r__1 = a[imax + jmax * a_dim1], f2c_abs(r__1));
 }
 else {
 rowmax = 0.f;
 }
 if (imax < *n) {
 i__1 = *n - imax;
 itemp = imax + isamax_(&i__1, &a[imax + 1 + imax * a_dim1] , &c__1);
 stemp = (r__1 = a[itemp + imax * a_dim1], f2c_abs(r__1));
 if (stemp > rowmax) {
 rowmax = stemp;
 jmax = itemp;
 }
 }
 /* Equivalent to testing for (used to handle NaN and Inf) */
 /* ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX */
 if (! ((r__1 = a[imax + imax * a_dim1], f2c_abs(r__1)) < alpha * rowmax)) {
 /* interchange rows and columns K and IMAX, */
 /* use 1-by-1 pivot block */
 kp = imax;
 done = TRUE_;
 /* Equivalent to testing for ROWMAX .EQ. COLMAX, */
 /* used to handle NaN and Inf */
 }
 else if (p == jmax || rowmax <= colmax) {
 /* interchange rows and columns K+1 and IMAX, */
 /* use 2-by-2 pivot block */
 kp = imax;
 kstep = 2;
 done = TRUE_;
 }
 else {
 /* Pivot NOT found, set variables and repeat */
 p = imax;
 colmax = rowmax;
 imax = jmax;
 }
 /* End pivot search loop body */
 if (! done) {
 goto L42;
 }
 }
 /* Swap TWO rows and TWO columns */
 /* First swap */
 if (kstep == 2 && p != k) {
 /* Interchange rows and column K and P in the trailing */
 /* submatrix A(k:n,k:n) if we have a 2-by-2 pivot */
 if (p < *n) {
 i__1 = *n - p;
 sswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p * a_dim1], &c__1);
 }
 if (p > k + 1) {
 i__1 = p - k - 1;
 sswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[p + (k + 1) * a_dim1], lda);
 }
 t = a[k + k * a_dim1];
 a[k + k * a_dim1] = a[p + p * a_dim1];
 a[p + p * a_dim1] = t;
 /* Convert lower triangle of A into L form by applying */
 /* the interchanges in columns 1:k-1. */
 if (k > 1) {
 i__1 = k - 1;
 sswap_(&i__1, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
 }
 }
 /* Second swap */
 kk = k + kstep - 1;
 if (kp != kk) {
 /* Interchange rows and columns KK and KP in the trailing */
 /* submatrix A(k:n,k:n) */
 if (kp < *n) {
 i__1 = *n - kp;
 sswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
 }
 if (kk < *n && kp > kk + 1) {
 i__1 = kp - kk - 1;
 sswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + ( kk + 1) * a_dim1], lda);
 }
 t = a[kk + kk * a_dim1];
 a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
 a[kp + kp * a_dim1] = t;
 if (kstep == 2) {
 t = a[k + 1 + k * a_dim1];
 a[k + 1 + k * a_dim1] = a[kp + k * a_dim1];
 a[kp + k * a_dim1] = t;
 }
 /* Convert lower triangle of A into L form by applying */
 /* the interchanges in columns 1:k-1. */
 if (k > 1) {
 i__1 = k - 1;
 sswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
 }
 }
 /* Update the trailing submatrix */
 if (kstep == 1) {
 /* 1-by-1 pivot block D(k): column k now holds */
 /* W(k) = L(k)*D(k) */
 /* where L(k) is the k-th column of L */
 if (k < *n) {
 /* Perform a rank-1 update of A(k+1:n,k+1:n) and */
 /* store L(k) in column k */
 if ((r__1 = a[k + k * a_dim1], f2c_abs(r__1)) >= sfmin) {
 /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
 /* A := A - L(k)*D(k)*L(k)**T */
 /* = A - W(k)*(1/D(k))*W(k)**T */
 d11 = 1.f / a[k + k * a_dim1];
 i__1 = *n - k;
 r__1 = -d11;
 ssyr_(uplo, &i__1, &r__1, &a[k + 1 + k * a_dim1], & c__1, &a[k + 1 + (k + 1) * a_dim1], lda);
 /* Store L(k) in column k */
 i__1 = *n - k;
 sscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
 }
 else {
 /* Store L(k) in column k */
 d11 = a[k + k * a_dim1];
 i__1 = *n;
 for (ii = k + 1;
 ii <= i__1;
 ++ii) {
 a[ii + k * a_dim1] /= d11;
 /* L46: */
 }
 /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
 /* A := A - L(k)*D(k)*L(k)**T */
 /* = A - W(k)*(1/D(k))*W(k)**T */
 /* = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */
 i__1 = *n - k;
 r__1 = -d11;
 ssyr_(uplo, &i__1, &r__1, &a[k + 1 + k * a_dim1], & c__1, &a[k + 1 + (k + 1) * a_dim1], lda);
 }
 /* Store the subdiagonal element of D in array E */
 e[k] = 0.f;
 }
 }
 else {
 /* 2-by-2 pivot block D(k): columns k and k+1 now hold */
 /* ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */
 /* where L(k) and L(k+1) are the k-th and (k+1)-th columns */
 /* of L */
 /* Perform a rank-2 update of A(k+2:n,k+2:n) as */
 /* A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T */
 /* = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T */
 /* and store L(k) and L(k+1) in columns k and k+1 */
 if (k < *n - 1) {
 d21 = a[k + 1 + k * a_dim1];
 d11 = a[k + 1 + (k + 1) * a_dim1] / d21;
 d22 = a[k + k * a_dim1] / d21;
 t = 1.f / (d11 * d22 - 1.f);
 i__1 = *n;
 for (j = k + 2;
 j <= i__1;
 ++j) {
 /* Compute D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */
 wk = t * (d11 * a[j + k * a_dim1] - a[j + (k + 1) * a_dim1]);
 wkp1 = t * (d22 * a[j + (k + 1) * a_dim1] - a[j + k * a_dim1]);
 /* Perform a rank-2 update of A(k+2:n,k+2:n) */
 i__2 = *n;
 for (i__ = j;
 i__ <= i__2;
 ++i__) {
 a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ + k * a_dim1] / d21 * wk - a[i__ + (k + 1) * a_dim1] / d21 * wkp1;
 /* L50: */
 }
 /* Store L(k) and L(k+1) in cols k and k+1 for row J */
 a[j + k * a_dim1] = wk / d21;
 a[j + (k + 1) * a_dim1] = wkp1 / d21;
 /* L60: */
 }
 }
 /* Copy subdiagonal elements of D(K) to E(K) and */
 /* ZERO out subdiagonal entry of A */
 e[k] = a[k + 1 + k * a_dim1];
 e[k + 1] = 0.f;
 a[k + 1 + k * a_dim1] = 0.f;
 }
 /* End column K is nonsingular */
 }
 /* Store details of the interchanges in IPIV */
 if (kstep == 1) {
 ipiv[k] = kp;
 }
 else {
 ipiv[k] = -p;
 ipiv[k + 1] = -kp;
 }
 /* Increase K and return to the start of the main loop */
 k += kstep;
 goto L40;
 L64: ;
 }
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of SSYTF2_RK */
 }
 /* ssytf2_rk__ */
 
