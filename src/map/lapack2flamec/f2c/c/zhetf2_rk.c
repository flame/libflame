/* ../netlib/v3.9.0/zhetf2_rk.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 /* > \brief \b ZHETF2_RK computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman (rook) diagonal pivoting method (BLAS2 unblocked algorithm). */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download ZHETF2_RK + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetf2_ rk.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetf2_ rk.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetf2_ rk.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE ZHETF2_RK( UPLO, N, A, LDA, E, IPIV, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* INTEGER INFO, LDA, N */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IPIV( * ) */
 /* COMPLEX*16 A( LDA, * ), E ( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > ZHETF2_RK computes the factorization of a complex Hermitian matrix A */
 /* > using the bounded Bunch-Kaufman (rook) diagonal pivoting method: */
 /* > */
 /* > A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T), */
 /* > */
 /* > where U (or L) is unit upper (or lower) triangular matrix, */
 /* > U**H (or L**H) is the conjugate of U (or L), P is a permutation */
 /* > matrix, P**T is the transpose of P, and D is Hermitian and block */
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
 /* > Hermitian matrix A is stored: */
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
 /* > On exit, contains: */
 /* > a) ONLY diagonal elements of the Hermitian block diagonal */
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
 /* > E is COMPLEX*16 array, dimension (N) */
 /* > On exit, contains the superdiagonal (or subdiagonal) */
 /* > elements of the Hermitian block diagonal matrix D */
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
 /* > the Hermitian block diagonal matrix D with 1-by-1 or 2-by-2 */
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
 /* > \ingroup complex16HEcomputational */
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
 int zhetf2_rk_(char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, integer *info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
 char buffer[256];
 snprintf(buffer, 256,"zhetf2_rk inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "",*uplo, *n, *lda);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
 doublereal d__1, d__2;
 doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;
 /* Builtin functions */
 double sqrt(doublereal), d_imag(doublecomplex *);
 void d_cnjg(doublecomplex *, doublecomplex *);
 /* Local variables */
 doublereal d__;
 integer i__, j, k, p;
 doublecomplex t;
 doublereal r1, d11;
 doublecomplex d12;
 doublereal d22;
 doublecomplex d21;
 integer ii, kk, kp;
 doublecomplex wk;
 doublereal tt;
 doublecomplex wkm1, wkp1;
 logical done;
 integer imax, jmax;
 extern /* Subroutine */
 int zher_(char *, integer *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *);
 doublereal alpha;
 extern logical lsame_(char *, char *);
 doublereal dtemp, sfmin;
 integer itemp, kstep;
 logical upper;
 extern /* Subroutine */
 int zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);
 extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *);
 doublereal absakk;
 extern /* Subroutine */
 int xerbla_(char *, integer *), zdscal_( integer *, doublereal *, doublecomplex *, integer *);
 doublereal colmax;
 extern integer izamax_(integer *, doublecomplex *, integer *);
 doublereal rowmax;
 /* -- LAPACK computational routine (version 3.7.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* December 2016 */
 /* .. Scalar Arguments .. */
 /* .. */
 /* .. Array Arguments .. */
 /* .. */
 /* ====================================================================== */
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
 /* .. Statement Functions .. */
 /* .. */
 /* .. Statement Function definitions .. */
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
 xerbla_("ZHETF2_RK", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Initialize ALPHA for use in choosing pivot block size. */
 alpha = (sqrt(17.) + 1.) / 8.;
 /* Compute machine safe minimum */
 sfmin = dlamch_("S");
 if (upper) {
 /* Factorize A as U*D*U**H using the upper triangle of A */
 /* Initialize the first entry of array E, where superdiagonal */
 /* elements of D are stored */
 e[1].r = 0.; e[1].i = 0.; // , expr subst  
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
 i__1 = k + k * a_dim1;
 absakk = (d__1 = a[i__1].r, f2c_dabs(d__1));
 /* IMAX is the row-index of the largest off-diagonal element in */
 /* column K, and COLMAX is its absolute value. */
 /* Determine both COLMAX and IMAX. */
 if (k > 1) {
 i__1 = k - 1;
 imax = izamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
 i__1 = imax + k * a_dim1;
 colmax = (d__1 = a[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&a[imax + k * a_dim1]), f2c_dabs(d__2));
 }
 else {
 colmax = 0.;
 }
 if (max(absakk,colmax) == 0.) {
 /* Column K is zero or underflow: set INFO and continue */
 if (*info == 0) {
 *info = k;
 }
 kp = k;
 i__1 = k + k * a_dim1;
 i__2 = k + k * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 /* Set E( K ) to zero */
 if (k > 1) {
 i__1 = k;
 e[i__1].r = 0.; e[i__1].i = 0.; // , expr subst  
 }
 }
 else {
 /* ============================================================ */
 /* BEGIN pivot search */
 /* Case(1) */
 /* Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
 /* (used to handle NaN and Inf) */
 if (! (absakk < alpha * colmax)) {
 /* no interchange, use 1-by-1 pivot block */
 kp = k;
 }
 else {
 done = FALSE_;
 /* Loop until pivot found */
 L12: /* BEGIN pivot search loop body */
 /* JMAX is the column-index of the largest off-diagonal */
 /* element in row IMAX, and ROWMAX is its absolute value. */
 /* Determine both ROWMAX and JMAX. */
 if (imax != k) {
 i__1 = k - imax;
 jmax = imax + izamax_(&i__1, &a[imax + (imax + 1) * a_dim1], lda);
 i__1 = imax + jmax * a_dim1;
 rowmax = (d__1 = a[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(& a[imax + jmax * a_dim1]), f2c_dabs(d__2));
 }
 else {
 rowmax = 0.;
 }
 if (imax > 1) {
 i__1 = imax - 1;
 itemp = izamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
 i__1 = itemp + imax * a_dim1;
 dtemp = (d__1 = a[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&a[ itemp + imax * a_dim1]), f2c_dabs(d__2));
 if (dtemp > rowmax) {
 rowmax = dtemp;
 jmax = itemp;
 }
 }
 /* Case(2) */
 /* Equivalent to testing for */
 /* ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX */
 /* (used to handle NaN and Inf) */
 i__1 = imax + imax * a_dim1;
 if (! ((d__1 = a[i__1].r, f2c_dabs(d__1)) < alpha * rowmax)) {
 /* interchange rows and columns K and IMAX, */
 /* use 1-by-1 pivot block */
 kp = imax;
 done = TRUE_;
 /* Case(3) */
 /* Equivalent to testing for ROWMAX.EQ.COLMAX, */
 /* (used to handle NaN and Inf) */
 }
 else if (p == jmax || rowmax <= colmax) {
 /* interchange rows and columns K-1 and IMAX, */
 /* use 2-by-2 pivot block */
 kp = imax;
 kstep = 2;
 done = TRUE_;
 /* Case(4) */
 }
 else {
 /* Pivot not found: set params and repeat */
 p = imax;
 colmax = rowmax;
 imax = jmax;
 }
 /* END pivot search loop body */
 if (! done) {
 goto L12;
 }
 }
 /* END pivot search */
 /* ============================================================ */
 /* KK is the column of A where pivoting step stopped */
 kk = k - kstep + 1;
 /* For only a 2x2 pivot, interchange rows and columns K and P */
 /* in the leading submatrix A(1:k,1:k) */
 if (kstep == 2 && p != k) {
 /* (1) Swap columnar parts */
 if (p > 1) {
 i__1 = p - 1;
 zswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &c__1);
 }
 /* (2) Swap and conjugate middle parts */
 i__1 = k - 1;
 for (j = p + 1;
 j <= i__1;
 ++j) {
 d_cnjg(&z__1, &a[j + k * a_dim1]);
 t.r = z__1.r; t.i = z__1.i; // , expr subst  
 i__2 = j + k * a_dim1;
 d_cnjg(&z__1, &a[p + j * a_dim1]);
 a[i__2].r = z__1.r; a[i__2].i = z__1.i; // , expr subst  
 i__2 = p + j * a_dim1;
 a[i__2].r = t.r; a[i__2].i = t.i; // , expr subst  
 /* L14: */
 }
 /* (3) Swap and conjugate corner elements at row-col interserction */
 i__1 = p + k * a_dim1;
 d_cnjg(&z__1, &a[p + k * a_dim1]);
 a[i__1].r = z__1.r; a[i__1].i = z__1.i; // , expr subst  
 /* (4) Swap diagonal elements at row-col intersection */
 i__1 = k + k * a_dim1;
 r1 = a[i__1].r;
 i__1 = k + k * a_dim1;
 i__2 = p + p * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 i__1 = p + p * a_dim1;
 a[i__1].r = r1; a[i__1].i = 0.; // , expr subst  
 /* Convert upper triangle of A into U form by applying */
 /* the interchanges in columns k+1:N. */
 if (k < *n) {
 i__1 = *n - k;
 zswap_(&i__1, &a[k + (k + 1) * a_dim1], lda, &a[p + (k + 1) * a_dim1], lda);
 }
 }
 /* For both 1x1 and 2x2 pivots, interchange rows and */
 /* columns KK and KP in the leading submatrix A(1:k,1:k) */
 if (kp != kk) {
 /* (1) Swap columnar parts */
 if (kp > 1) {
 i__1 = kp - 1;
 zswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
 }
 /* (2) Swap and conjugate middle parts */
 i__1 = kk - 1;
 for (j = kp + 1;
 j <= i__1;
 ++j) {
 d_cnjg(&z__1, &a[j + kk * a_dim1]);
 t.r = z__1.r; t.i = z__1.i; // , expr subst  
 i__2 = j + kk * a_dim1;
 d_cnjg(&z__1, &a[kp + j * a_dim1]);
 a[i__2].r = z__1.r; a[i__2].i = z__1.i; // , expr subst  
 i__2 = kp + j * a_dim1;
 a[i__2].r = t.r; a[i__2].i = t.i; // , expr subst  
 /* L15: */
 }
 /* (3) Swap and conjugate corner elements at row-col interserction */
 i__1 = kp + kk * a_dim1;
 d_cnjg(&z__1, &a[kp + kk * a_dim1]);
 a[i__1].r = z__1.r; a[i__1].i = z__1.i; // , expr subst  
 /* (4) Swap diagonal elements at row-col intersection */
 i__1 = kk + kk * a_dim1;
 r1 = a[i__1].r;
 i__1 = kk + kk * a_dim1;
 i__2 = kp + kp * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 i__1 = kp + kp * a_dim1;
 a[i__1].r = r1; a[i__1].i = 0.; // , expr subst  
 if (kstep == 2) {
 /* (*) Make sure that diagonal element of pivot is real */
 i__1 = k + k * a_dim1;
 i__2 = k + k * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 /* (5) Swap row elements */
 i__1 = k - 1 + k * a_dim1;
 t.r = a[i__1].r; t.i = a[i__1].i; // , expr subst  
 i__1 = k - 1 + k * a_dim1;
 i__2 = kp + k * a_dim1;
 a[i__1].r = a[i__2].r; a[i__1].i = a[i__2].i; // , expr subst  
 i__1 = kp + k * a_dim1;
 a[i__1].r = t.r; a[i__1].i = t.i; // , expr subst  
 }
 /* Convert upper triangle of A into U form by applying */
 /* the interchanges in columns k+1:N. */
 if (k < *n) {
 i__1 = *n - k;
 zswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k + 1) * a_dim1], lda);
 }
 }
 else {
 /* (*) Make sure that diagonal element of pivot is real */
 i__1 = k + k * a_dim1;
 i__2 = k + k * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 if (kstep == 2) {
 i__1 = k - 1 + (k - 1) * a_dim1;
 i__2 = k - 1 + (k - 1) * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
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
 i__1 = k + k * a_dim1;
 if ((d__1 = a[i__1].r, f2c_dabs(d__1)) >= sfmin) {
 /* Perform a rank-1 update of A(1:k-1,1:k-1) as */
 /* A := A - U(k)*D(k)*U(k)**T */
 /* = A - W(k)*1/D(k)*W(k)**T */
 i__1 = k + k * a_dim1;
 d11 = 1. / a[i__1].r;
 i__1 = k - 1;
 d__1 = -d11;
 zher_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, & a[a_offset], lda);
 /* Store U(k) in column k */
 i__1 = k - 1;
 zdscal_(&i__1, &d11, &a[k * a_dim1 + 1], &c__1);
 }
 else {
 /* Store L(k) in column K */
 i__1 = k + k * a_dim1;
 d11 = a[i__1].r;
 i__1 = k - 1;
 for (ii = 1;
 ii <= i__1;
 ++ii) {
 i__2 = ii + k * a_dim1;
 i__3 = ii + k * a_dim1;
 z__1.r = a[i__3].r / d11; z__1.i = a[i__3].i / d11; // , expr subst  
 a[i__2].r = z__1.r; a[i__2].i = z__1.i; // , expr subst  
 /* L16: */
 }
 /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
 /* A := A - U(k)*D(k)*U(k)**T */
 /* = A - W(k)*(1/D(k))*W(k)**T */
 /* = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */
 i__1 = k - 1;
 d__1 = -d11;
 zher_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, & a[a_offset], lda);
 }
 /* Store the superdiagonal element of D in array E */
 i__1 = k;
 e[i__1].r = 0.; e[i__1].i = 0.; // , expr subst  
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
 /* D = |A12| */
 i__1 = k - 1 + k * a_dim1;
 d__1 = a[i__1].r;
 d__2 = d_imag(&a[k - 1 + k * a_dim1]);
 d__ = dlapy2_(&d__1, &d__2);
 i__1 = k + k * a_dim1;
 z__1.r = a[i__1].r / d__; z__1.i = a[i__1].i / d__; // , expr subst  
 d11 = z__1.r;
 i__1 = k - 1 + (k - 1) * a_dim1;
 z__1.r = a[i__1].r / d__; z__1.i = a[i__1].i / d__; // , expr subst  
 d22 = z__1.r;
 i__1 = k - 1 + k * a_dim1;
 z__1.r = a[i__1].r / d__; z__1.i = a[i__1].i / d__; // , expr subst  
 d12.r = z__1.r; d12.i = z__1.i; // , expr subst  
 tt = 1. / (d11 * d22 - 1.);
 for (j = k - 2;
 j >= 1;
 --j) {
 /* Compute D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */
 i__1 = j + (k - 1) * a_dim1;
 z__3.r = d11 * a[i__1].r; z__3.i = d11 * a[i__1].i; // , expr subst  
 d_cnjg(&z__5, &d12);
 i__2 = j + k * a_dim1;
 z__4.r = z__5.r * a[i__2].r - z__5.i * a[i__2].i; z__4.i = z__5.r * a[i__2].i + z__5.i * a[i__2] .r; // , expr subst  
 z__2.r = z__3.r - z__4.r; z__2.i = z__3.i - z__4.i; // , expr subst  
 z__1.r = tt * z__2.r; z__1.i = tt * z__2.i; // , expr subst  
 wkm1.r = z__1.r; wkm1.i = z__1.i; // , expr subst  
 i__1 = j + k * a_dim1;
 z__3.r = d22 * a[i__1].r; z__3.i = d22 * a[i__1].i; // , expr subst  
 i__2 = j + (k - 1) * a_dim1;
 z__4.r = d12.r * a[i__2].r - d12.i * a[i__2].i; z__4.i = d12.r * a[i__2].i + d12.i * a[i__2] .r; // , expr subst  
 z__2.r = z__3.r - z__4.r; z__2.i = z__3.i - z__4.i; // , expr subst  
 z__1.r = tt * z__2.r; z__1.i = tt * z__2.i; // , expr subst  
 wk.r = z__1.r; wk.i = z__1.i; // , expr subst  
 /* Perform a rank-2 update of A(1:k-2,1:k-2) */
 for (i__ = j;
 i__ >= 1;
 --i__) {
 i__1 = i__ + j * a_dim1;
 i__2 = i__ + j * a_dim1;
 i__3 = i__ + k * a_dim1;
 z__4.r = a[i__3].r / d__; z__4.i = a[i__3].i / d__; // , expr subst  
 d_cnjg(&z__5, &wk);
 z__3.r = z__4.r * z__5.r - z__4.i * z__5.i; z__3.i = z__4.r * z__5.i + z__4.i * z__5.r; // , expr subst  
 z__2.r = a[i__2].r - z__3.r; z__2.i = a[i__2].i - z__3.i; // , expr subst  
 i__4 = i__ + (k - 1) * a_dim1;
 z__7.r = a[i__4].r / d__; z__7.i = a[i__4].i / d__; // , expr subst  
 d_cnjg(&z__8, &wkm1);
 z__6.r = z__7.r * z__8.r - z__7.i * z__8.i; z__6.i = z__7.r * z__8.i + z__7.i * z__8.r; // , expr subst  
 z__1.r = z__2.r - z__6.r; z__1.i = z__2.i - z__6.i; // , expr subst  
 a[i__1].r = z__1.r; a[i__1].i = z__1.i; // , expr subst  
 /* L20: */
 }
 /* Store U(k) and U(k-1) in cols k and k-1 for row J */
 i__1 = j + k * a_dim1;
 z__1.r = wk.r / d__; z__1.i = wk.i / d__; // , expr subst  
 a[i__1].r = z__1.r; a[i__1].i = z__1.i; // , expr subst  
 i__1 = j + (k - 1) * a_dim1;
 z__1.r = wkm1.r / d__; z__1.i = wkm1.i / d__; // , expr subst  
 a[i__1].r = z__1.r; a[i__1].i = z__1.i; // , expr subst  
 /* (*) Make sure that diagonal element of pivot is real */
 i__1 = j + j * a_dim1;
 i__2 = j + j * a_dim1;
 d__1 = a[i__2].r;
 z__1.r = d__1; z__1.i = 0.; // , expr subst  
 a[i__1].r = z__1.r; a[i__1].i = z__1.i; // , expr subst  
 /* L30: */
 }
 }
 /* Copy superdiagonal elements of D(K) to E(K) and */
 /* ZERO out superdiagonal entry of A */
 i__1 = k;
 i__2 = k - 1 + k * a_dim1;
 e[i__1].r = a[i__2].r; e[i__1].i = a[i__2].i; // , expr subst  
 i__1 = k - 1;
 e[i__1].r = 0.; e[i__1].i = 0.; // , expr subst  
 i__1 = k - 1 + k * a_dim1;
 a[i__1].r = 0.; a[i__1].i = 0.; // , expr subst  
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
 /* Factorize A as L*D*L**H using the lower triangle of A */
 /* Initialize the unused last entry of the subdiagonal array E. */
 i__1 = *n;
 e[i__1].r = 0.; e[i__1].i = 0.; // , expr subst  
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
 i__1 = k + k * a_dim1;
 absakk = (d__1 = a[i__1].r, f2c_dabs(d__1));
 /* IMAX is the row-index of the largest off-diagonal element in */
 /* column K, and COLMAX is its absolute value. */
 /* Determine both COLMAX and IMAX. */
 if (k < *n) {
 i__1 = *n - k;
 imax = k + izamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
 i__1 = imax + k * a_dim1;
 colmax = (d__1 = a[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&a[imax + k * a_dim1]), f2c_dabs(d__2));
 }
 else {
 colmax = 0.;
 }
 if (max(absakk,colmax) == 0.) {
 /* Column K is zero or underflow: set INFO and continue */
 if (*info == 0) {
 *info = k;
 }
 kp = k;
 i__1 = k + k * a_dim1;
 i__2 = k + k * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 /* Set E( K ) to zero */
 if (k < *n) {
 i__1 = k;
 e[i__1].r = 0.; e[i__1].i = 0.; // , expr subst  
 }
 }
 else {
 /* ============================================================ */
 /* BEGIN pivot search */
 /* Case(1) */
 /* Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
 /* (used to handle NaN and Inf) */
 if (! (absakk < alpha * colmax)) {
 /* no interchange, use 1-by-1 pivot block */
 kp = k;
 }
 else {
 done = FALSE_;
 /* Loop until pivot found */
 L42: /* BEGIN pivot search loop body */
 /* JMAX is the column-index of the largest off-diagonal */
 /* element in row IMAX, and ROWMAX is its absolute value. */
 /* Determine both ROWMAX and JMAX. */
 if (imax != k) {
 i__1 = imax - k;
 jmax = k - 1 + izamax_(&i__1, &a[imax + k * a_dim1], lda);
 i__1 = imax + jmax * a_dim1;
 rowmax = (d__1 = a[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(& a[imax + jmax * a_dim1]), f2c_dabs(d__2));
 }
 else {
 rowmax = 0.;
 }
 if (imax < *n) {
 i__1 = *n - imax;
 itemp = imax + izamax_(&i__1, &a[imax + 1 + imax * a_dim1] , &c__1);
 i__1 = itemp + imax * a_dim1;
 dtemp = (d__1 = a[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&a[ itemp + imax * a_dim1]), f2c_dabs(d__2));
 if (dtemp > rowmax) {
 rowmax = dtemp;
 jmax = itemp;
 }
 }
 /* Case(2) */
 /* Equivalent to testing for */
 /* ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX */
 /* (used to handle NaN and Inf) */
 i__1 = imax + imax * a_dim1;
 if (! ((d__1 = a[i__1].r, f2c_dabs(d__1)) < alpha * rowmax)) {
 /* interchange rows and columns K and IMAX, */
 /* use 1-by-1 pivot block */
 kp = imax;
 done = TRUE_;
 /* Case(3) */
 /* Equivalent to testing for ROWMAX.EQ.COLMAX, */
 /* (used to handle NaN and Inf) */
 }
 else if (p == jmax || rowmax <= colmax) {
 /* interchange rows and columns K+1 and IMAX, */
 /* use 2-by-2 pivot block */
 kp = imax;
 kstep = 2;
 done = TRUE_;
 /* Case(4) */
 }
 else {
 /* Pivot not found: set params and repeat */
 p = imax;
 colmax = rowmax;
 imax = jmax;
 }
 /* END pivot search loop body */
 if (! done) {
 goto L42;
 }
 }
 /* END pivot search */
 /* ============================================================ */
 /* KK is the column of A where pivoting step stopped */
 kk = k + kstep - 1;
 /* For only a 2x2 pivot, interchange rows and columns K and P */
 /* in the trailing submatrix A(k:n,k:n) */
 if (kstep == 2 && p != k) {
 /* (1) Swap columnar parts */
 if (p < *n) {
 i__1 = *n - p;
 zswap_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + 1 + p * a_dim1], &c__1);
 }
 /* (2) Swap and conjugate middle parts */
 i__1 = p - 1;
 for (j = k + 1;
 j <= i__1;
 ++j) {
 d_cnjg(&z__1, &a[j + k * a_dim1]);
 t.r = z__1.r; t.i = z__1.i; // , expr subst  
 i__2 = j + k * a_dim1;
 d_cnjg(&z__1, &a[p + j * a_dim1]);
 a[i__2].r = z__1.r; a[i__2].i = z__1.i; // , expr subst  
 i__2 = p + j * a_dim1;
 a[i__2].r = t.r; a[i__2].i = t.i; // , expr subst  
 /* L44: */
 }
 /* (3) Swap and conjugate corner elements at row-col interserction */
 i__1 = p + k * a_dim1;
 d_cnjg(&z__1, &a[p + k * a_dim1]);
 a[i__1].r = z__1.r; a[i__1].i = z__1.i; // , expr subst  
 /* (4) Swap diagonal elements at row-col intersection */
 i__1 = k + k * a_dim1;
 r1 = a[i__1].r;
 i__1 = k + k * a_dim1;
 i__2 = p + p * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 i__1 = p + p * a_dim1;
 a[i__1].r = r1; a[i__1].i = 0.; // , expr subst  
 /* Convert lower triangle of A into L form by applying */
 /* the interchanges in columns 1:k-1. */
 if (k > 1) {
 i__1 = k - 1;
 zswap_(&i__1, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
 }
 }
 /* For both 1x1 and 2x2 pivots, interchange rows and */
 /* columns KK and KP in the trailing submatrix A(k:n,k:n) */
 if (kp != kk) {
 /* (1) Swap columnar parts */
 if (kp < *n) {
 i__1 = *n - kp;
 zswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
 }
 /* (2) Swap and conjugate middle parts */
 i__1 = kp - 1;
 for (j = kk + 1;
 j <= i__1;
 ++j) {
 d_cnjg(&z__1, &a[j + kk * a_dim1]);
 t.r = z__1.r; t.i = z__1.i; // , expr subst  
 i__2 = j + kk * a_dim1;
 d_cnjg(&z__1, &a[kp + j * a_dim1]);
 a[i__2].r = z__1.r; a[i__2].i = z__1.i; // , expr subst  
 i__2 = kp + j * a_dim1;
 a[i__2].r = t.r; a[i__2].i = t.i; // , expr subst  
 /* L45: */
 }
 /* (3) Swap and conjugate corner elements at row-col interserction */
 i__1 = kp + kk * a_dim1;
 d_cnjg(&z__1, &a[kp + kk * a_dim1]);
 a[i__1].r = z__1.r; a[i__1].i = z__1.i; // , expr subst  
 /* (4) Swap diagonal elements at row-col intersection */
 i__1 = kk + kk * a_dim1;
 r1 = a[i__1].r;
 i__1 = kk + kk * a_dim1;
 i__2 = kp + kp * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 i__1 = kp + kp * a_dim1;
 a[i__1].r = r1; a[i__1].i = 0.; // , expr subst  
 if (kstep == 2) {
 /* (*) Make sure that diagonal element of pivot is real */
 i__1 = k + k * a_dim1;
 i__2 = k + k * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 /* (5) Swap row elements */
 i__1 = k + 1 + k * a_dim1;
 t.r = a[i__1].r; t.i = a[i__1].i; // , expr subst  
 i__1 = k + 1 + k * a_dim1;
 i__2 = kp + k * a_dim1;
 a[i__1].r = a[i__2].r; a[i__1].i = a[i__2].i; // , expr subst  
 i__1 = kp + k * a_dim1;
 a[i__1].r = t.r; a[i__1].i = t.i; // , expr subst  
 }
 /* Convert lower triangle of A into L form by applying */
 /* the interchanges in columns 1:k-1. */
 if (k > 1) {
 i__1 = k - 1;
 zswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
 }
 }
 else {
 /* (*) Make sure that diagonal element of pivot is real */
 i__1 = k + k * a_dim1;
 i__2 = k + k * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 if (kstep == 2) {
 i__1 = k + 1 + (k + 1) * a_dim1;
 i__2 = k + 1 + (k + 1) * a_dim1;
 d__1 = a[i__2].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 }
 }
 /* Update the trailing submatrix */
 if (kstep == 1) {
 /* 1-by-1 pivot block D(k): column k of A now holds */
 /* W(k) = L(k)*D(k), */
 /* where L(k) is the k-th column of L */
 if (k < *n) {
 /* Perform a rank-1 update of A(k+1:n,k+1:n) and */
 /* store L(k) in column k */
 /* Handle division by a small number */
 i__1 = k + k * a_dim1;
 if ((d__1 = a[i__1].r, f2c_dabs(d__1)) >= sfmin) {
 /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
 /* A := A - L(k)*D(k)*L(k)**T */
 /* = A - W(k)*(1/D(k))*W(k)**T */
 i__1 = k + k * a_dim1;
 d11 = 1. / a[i__1].r;
 i__1 = *n - k;
 d__1 = -d11;
 zher_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], & c__1, &a[k + 1 + (k + 1) * a_dim1], lda);
 /* Store L(k) in column k */
 i__1 = *n - k;
 zdscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
 }
 else {
 /* Store L(k) in column k */
 i__1 = k + k * a_dim1;
 d11 = a[i__1].r;
 i__1 = *n;
 for (ii = k + 1;
 ii <= i__1;
 ++ii) {
 i__2 = ii + k * a_dim1;
 i__3 = ii + k * a_dim1;
 z__1.r = a[i__3].r / d11; z__1.i = a[i__3].i / d11; // , expr subst  
 a[i__2].r = z__1.r; a[i__2].i = z__1.i; // , expr subst  
 /* L46: */
 }
 /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
 /* A := A - L(k)*D(k)*L(k)**T */
 /* = A - W(k)*(1/D(k))*W(k)**T */
 /* = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T */
 i__1 = *n - k;
 d__1 = -d11;
 zher_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], & c__1, &a[k + 1 + (k + 1) * a_dim1], lda);
 }
 /* Store the subdiagonal element of D in array E */
 i__1 = k;
 e[i__1].r = 0.; e[i__1].i = 0.; // , expr subst  
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
 /* D = |A21| */
 i__1 = k + 1 + k * a_dim1;
 d__1 = a[i__1].r;
 d__2 = d_imag(&a[k + 1 + k * a_dim1]);
 d__ = dlapy2_(&d__1, &d__2);
 i__1 = k + 1 + (k + 1) * a_dim1;
 d11 = a[i__1].r / d__;
 i__1 = k + k * a_dim1;
 d22 = a[i__1].r / d__;
 i__1 = k + 1 + k * a_dim1;
 z__1.r = a[i__1].r / d__; z__1.i = a[i__1].i / d__; // , expr subst  
 d21.r = z__1.r; d21.i = z__1.i; // , expr subst  
 tt = 1. / (d11 * d22 - 1.);
 i__1 = *n;
 for (j = k + 2;
 j <= i__1;
 ++j) {
 /* Compute D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J */
 i__2 = j + k * a_dim1;
 z__3.r = d11 * a[i__2].r; z__3.i = d11 * a[i__2].i; // , expr subst  
 i__3 = j + (k + 1) * a_dim1;
 z__4.r = d21.r * a[i__3].r - d21.i * a[i__3].i; z__4.i = d21.r * a[i__3].i + d21.i * a[i__3] .r; // , expr subst  
 z__2.r = z__3.r - z__4.r; z__2.i = z__3.i - z__4.i; // , expr subst  
 z__1.r = tt * z__2.r; z__1.i = tt * z__2.i; // , expr subst  
 wk.r = z__1.r; wk.i = z__1.i; // , expr subst  
 i__2 = j + (k + 1) * a_dim1;
 z__3.r = d22 * a[i__2].r; z__3.i = d22 * a[i__2].i; // , expr subst  
 d_cnjg(&z__5, &d21);
 i__3 = j + k * a_dim1;
 z__4.r = z__5.r * a[i__3].r - z__5.i * a[i__3].i; z__4.i = z__5.r * a[i__3].i + z__5.i * a[i__3] .r; // , expr subst  
 z__2.r = z__3.r - z__4.r; z__2.i = z__3.i - z__4.i; // , expr subst  
 z__1.r = tt * z__2.r; z__1.i = tt * z__2.i; // , expr subst  
 wkp1.r = z__1.r; wkp1.i = z__1.i; // , expr subst  
 /* Perform a rank-2 update of A(k+2:n,k+2:n) */
 i__2 = *n;
 for (i__ = j;
 i__ <= i__2;
 ++i__) {
 i__3 = i__ + j * a_dim1;
 i__4 = i__ + j * a_dim1;
 i__5 = i__ + k * a_dim1;
 z__4.r = a[i__5].r / d__; z__4.i = a[i__5].i / d__; // , expr subst  
 d_cnjg(&z__5, &wk);
 z__3.r = z__4.r * z__5.r - z__4.i * z__5.i; z__3.i = z__4.r * z__5.i + z__4.i * z__5.r; // , expr subst  
 z__2.r = a[i__4].r - z__3.r; z__2.i = a[i__4].i - z__3.i; // , expr subst  
 i__6 = i__ + (k + 1) * a_dim1;
 z__7.r = a[i__6].r / d__; z__7.i = a[i__6].i / d__; // , expr subst  
 d_cnjg(&z__8, &wkp1);
 z__6.r = z__7.r * z__8.r - z__7.i * z__8.i; z__6.i = z__7.r * z__8.i + z__7.i * z__8.r; // , expr subst  
 z__1.r = z__2.r - z__6.r; z__1.i = z__2.i - z__6.i; // , expr subst  
 a[i__3].r = z__1.r; a[i__3].i = z__1.i; // , expr subst  
 /* L50: */
 }
 /* Store L(k) and L(k+1) in cols k and k+1 for row J */
 i__2 = j + k * a_dim1;
 z__1.r = wk.r / d__; z__1.i = wk.i / d__; // , expr subst  
 a[i__2].r = z__1.r; a[i__2].i = z__1.i; // , expr subst  
 i__2 = j + (k + 1) * a_dim1;
 z__1.r = wkp1.r / d__; z__1.i = wkp1.i / d__; // , expr subst  
 a[i__2].r = z__1.r; a[i__2].i = z__1.i; // , expr subst  
 /* (*) Make sure that diagonal element of pivot is real */
 i__2 = j + j * a_dim1;
 i__3 = j + j * a_dim1;
 d__1 = a[i__3].r;
 z__1.r = d__1; z__1.i = 0.; // , expr subst  
 a[i__2].r = z__1.r; a[i__2].i = z__1.i; // , expr subst  
 /* L60: */
 }
 }
 /* Copy subdiagonal elements of D(K) to E(K) and */
 /* ZERO out subdiagonal entry of A */
 i__1 = k;
 i__2 = k + 1 + k * a_dim1;
 e[i__1].r = a[i__2].r; e[i__1].i = a[i__2].i; // , expr subst  
 i__1 = k + 1;
 e[i__1].r = 0.; e[i__1].i = 0.; // , expr subst  
 i__1 = k + 1 + k * a_dim1;
 a[i__1].r = 0.; a[i__1].i = 0.; // , expr subst  
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
 /* End of ZHETF2_RK */
 }
 /* zhetf2_rk__ */
 
