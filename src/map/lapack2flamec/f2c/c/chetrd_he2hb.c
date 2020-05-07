/* ../netlib/v3.9.0/chetrd_he2hb.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static complex c_b1 = {
0.f,0.f}
;
 static complex c_b2 = {
1.f,0.f}
;
 static integer c__4 = 4;
 static integer c_n1 = -1;
 static integer c__1 = 1;
 static real c_b33 = 1.f;
 /* > \brief \b CHETRD_HE2HB */
 /* @generated from zhetrd_he2hb.f, fortran z -> c, Wed Dec 7 08:22:40 2016 */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download CHETRD_HE2HB + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrd. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrd. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrd. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE CHETRD_HE2HB( UPLO, N, KD, A, LDA, AB, LDAB, TAU, */
 /* WORK, LWORK, INFO ) */
 /* IMPLICIT NONE */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* INTEGER INFO, LDA, LDAB, LWORK, N, KD */
 /* .. */
 /* .. Array Arguments .. */
 /* COMPLEX A( LDA, * ), AB( LDAB, * ), */
 /* TAU( * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > CHETRD_HE2HB reduces a complex Hermitian matrix A to complex Hermitian */
 /* > band-diagonal form AB by a unitary similarity transformation: */
 /* > Q**H * A * Q = AB. */
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
 /* > \param[in] KD */
 /* > \verbatim */
 /* > KD is INTEGER */
 /* > The number of superdiagonals of the reduced matrix if UPLO = 'U', */
 /* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
 /* > The reduced matrix is stored in the array AB. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] A */
 /* > \verbatim */
 /* > A is COMPLEX array, dimension (LDA,N) */
 /* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
 /* > N-by-N upper triangular part of A contains the upper */
 /* > triangular part of the matrix A, and the strictly lower */
 /* > triangular part of A is not referenced. If UPLO = 'L', the */
 /* > leading N-by-N lower triangular part of A contains the lower */
 /* > triangular part of the matrix A, and the strictly upper */
 /* > triangular part of A is not referenced. */
 /* > On exit, if UPLO = 'U', the diagonal and first superdiagonal */
 /* > of A are overwritten by the corresponding elements of the */
 /* > tridiagonal matrix T, and the elements above the first */
 /* > superdiagonal, with the array TAU, represent the unitary */
 /* > matrix Q as a product of elementary reflectors;
 if UPLO */
 /* > = 'L', the diagonal and first subdiagonal of A are over- */
 /* > written by the corresponding elements of the tridiagonal */
 /* > matrix T, and the elements below the first subdiagonal, with */
 /* > the array TAU, represent the unitary matrix Q as a product */
 /* > of elementary reflectors. See Further Details. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,N). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] AB */
 /* > \verbatim */
 /* > AB is COMPLEX array, dimension (LDAB,N) */
 /* > On exit, the upper or lower triangle of the Hermitian band */
 /* > matrix A, stored in the first KD+1 rows of the array. The */
 /* > j-th column of A is stored in the j-th column of the array AB */
 /* > as follows: */
 /* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
 */
 /* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDAB */
 /* > \verbatim */
 /* > LDAB is INTEGER */
 /* > The leading dimension of the array AB. LDAB >= KD+1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] TAU */
 /* > \verbatim */
 /* > TAU is COMPLEX array, dimension (N-KD) */
 /* > The scalar factors of the elementary reflectors (see Further */
 /* > Details). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is COMPLEX array, dimension (LWORK) */
 /* > On exit, if INFO = 0, or if LWORK=-1, */
 /* > WORK(1) returns the size of LWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The dimension of the array WORK which should be calculated */
 /* > by a workspace query. LWORK = MAX(1, LWORK_QUERY) */
 /* > If LWORK = -1, then a workspace query is assumed;
 the routine */
 /* > only calculates the optimal size of the WORK array, returns */
 /* > this value as the first entry of the WORK array, and no error */
 /* > message related to LWORK is issued by XERBLA. */
 /* > LWORK_QUERY = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD */
 /* > where FACTOPTNB is the blocking used by the QR or LQ */
 /* > algorithm, usually FACTOPTNB=128 is a good choice otherwise */
 /* > putting LWORK=-1 will provide the size of WORK. */
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
 /* > \ingroup complexHEcomputational */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > Implemented by Azzam Haidar. */
 /* > */
 /* > All details are available on technical report, SC11, SC13 papers. */
 /* > */
 /* > Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
 /* > Parallel reduction to condensed forms for symmetric eigenvalue problems */
 /* > using aggregated fine-grained and memory-aware kernels. In Proceedings */
 /* > of 2011 International Conference for High Performance Computing, */
 /* > Networking, Storage and Analysis (SC '11), New York, NY, USA, */
 /* > Article 8 , 11 pages. */
 /* > http://doi.acm.org/10.1145/2063384.2063394 */
 /* > */
 /* > A. Haidar, J. Kurzak, P. Luszczek, 2013. */
 /* > An improved parallel singular value algorithm and its implementation */
 /* > for multicore hardware, In Proceedings of 2013 International Conference */
 /* > for High Performance Computing, Networking, Storage and Analysis (SC '13). */
 /* > Denver, Colorado, USA, 2013. */
 /* > Article 90, 12 pages. */
 /* > http://doi.acm.org/10.1145/2503210.2503292 */
 /* > */
 /* > A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
 /* > A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
 /* > calculations based on fine-grained memory aware tasks. */
 /* > International Journal of High Performance Computing Applications. */
 /* > Volume 28 Issue 2, Pages 196-209, May 2014. */
 /* > http://hpc.sagepub.com/content/28/2/196 */
 /* > */
 /* > \endverbatim */
 /* > */
 /* > \verbatim */
 /* > */
 /* > If UPLO = 'U', the matrix Q is represented as a product of elementary */
 /* > reflectors */
 /* > */
 /* > Q = H(k)**H . . . H(2)**H H(1)**H, where k = n-kd. */
 /* > */
 /* > Each H(i) has the form */
 /* > */
 /* > H(i) = I - tau * v * v**H */
 /* > */
 /* > where tau is a complex scalar, and v is a complex vector with */
 /* > v(1:i+kd-1) = 0 and v(i+kd) = 1;
 conjg(v(i+kd+1:n)) is stored on exit in */
 /* > A(i,i+kd+1:n), and tau in TAU(i). */
 /* > */
 /* > If UPLO = 'L', the matrix Q is represented as a product of elementary */
 /* > reflectors */
 /* > */
 /* > Q = H(1) H(2) . . . H(k), where k = n-kd. */
 /* > */
 /* > Each H(i) has the form */
 /* > */
 /* > H(i) = I - tau * v * v**H */
 /* > */
 /* > where tau is a complex scalar, and v is a complex vector with */
 /* > v(kd+1:i) = 0 and v(i+kd+1) = 1;
 v(i+kd+2:n) is stored on exit in */
 /* > A(i+kd+2:n,i), and tau in TAU(i). */
 /* > */
 /* > The contents of A on exit are illustrated by the following examples */
 /* > with n = 5: */
 /* > */
 /* > if UPLO = 'U': if UPLO = 'L': */
 /* > */
 /* > ( ab ab/v1 v1 v1 v1 ) ( ab ) */
 /* > ( ab ab/v2 v2 v2 ) ( ab/v1 ab ) */
 /* > ( ab ab/v3 v3 ) ( v1 ab/v2 ab ) */
 /* > ( ab ab/v4 ) ( v1 v2 ab/v3 ab ) */
 /* > ( ab ) ( v1 v2 v3 ab/v4 ab ) */
 /* > */
 /* > where d and e denote diagonal and off-diagonal elements of T, and vi */
 /* > denotes an element of the vector defining H(i). */
 /* > \endverbatim */
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int chetrd_he2hb_(char *uplo, integer *n, integer *kd, complex *a, integer *lda, complex *ab, integer *ldab, complex *tau, complex *work, integer *lwork, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
 complex q__1;
 /* Local variables */
 integer i__, j, lk, pk, pn, lt, lw, ls1, ls2, ldt, ldw, lds1, lds2;
 extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
 integer tpos, wpos, s1pos, s2pos;
 extern /* Subroutine */
 int cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *), chemm_(char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
 extern logical lsame_(char *, char *);
 integer iinfo;
 extern /* Subroutine */
 int ccopy_(integer *, complex *, integer *, complex *, integer *);
 integer lwmin;
 logical upper;
 extern /* Subroutine */
 int cher2k_(char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, real *, complex *, integer *), cgelqf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), cgeqrf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), clarft_(char *, char *, integer *, integer *, complex *, integer *, complex *, complex *, integer *), claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *), xerbla_(char *, integer *);
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
 /* .. Parameters .. */
 /* .. */
 /* .. Local Scalars .. */
 /* .. */
 /* .. External Subroutines .. */
 /* .. */
 /* .. Intrinsic Functions .. */
 /* .. */
 /* .. External Functions .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Determine the minimal workspace size required */
 /* and test the input parameters */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 ab_dim1 = *ldab;
 ab_offset = 1 + ab_dim1;
 ab -= ab_offset;
 --tau;
 --work;
 /* Function Body */
 *info = 0;
 upper = lsame_(uplo, "U");
 lquery = *lwork == -1;
 lwmin = ilaenv2stage_(&c__4, "CHETRD_HE2HB", "", n, kd, &c_n1, &c_n1);
 if (! upper && ! lsame_(uplo, "L")) {
 *info = -1;
 }
 else if (*n < 0) {
 *info = -2;
 }
 else if (*kd < 0) {
 *info = -3;
 }
 else if (*lda < max(1,*n)) {
 *info = -5;
 }
 else /* if(complicated condition) */
 {
 /* Computing MAX */
 i__1 = 1; i__2 = *kd + 1; // , expr subst  
 if (*ldab < max(i__1,i__2)) {
 *info = -7;
 }
 else if (*lwork < lwmin && ! lquery) {
 *info = -10;
 }
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("CHETRD_HE2HB", &i__1);
 return 0;
 }
 else if (lquery) {
 work[1].r = (real) lwmin; work[1].i = 0.f; // , expr subst  
 return 0;
 }
 /* Quick return if possible */
 /* Copy the upper/lower portion of A into AB */
 if (*n <= *kd + 1) {
 if (upper) {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 /* Computing MIN */
 i__2 = *kd + 1;
 lk = min(i__2,i__);
 ccopy_(&lk, &a[i__ - lk + 1 + i__ * a_dim1], &c__1, &ab[*kd + 1 - lk + 1 + i__ * ab_dim1], &c__1);
 /* L100: */
 }
 }
 else {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 /* Computing MIN */
 i__2 = *kd + 1; i__3 = *n - i__ + 1; // , expr subst  
 lk = min(i__2,i__3);
 ccopy_(&lk, &a[i__ + i__ * a_dim1], &c__1, &ab[i__ * ab_dim1 + 1], &c__1);
 /* L110: */
 }
 }
 work[1].r = 1.f; work[1].i = 0.f; // , expr subst  
 return 0;
 }
 /* Determine the pointer position for the workspace */
 ldt = *kd;
 lds1 = *kd;
 lt = ldt * *kd;
 lw = *n * *kd;
 ls1 = lds1 * *kd;
 ls2 = lwmin - lt - lw - ls1;
 /* LS2 = N*MAX(KD,FACTOPTNB) */
 tpos = 1;
 wpos = tpos + lt;
 s1pos = wpos + lw;
 s2pos = s1pos + ls1;
 if (upper) {
 ldw = *kd;
 lds2 = *kd;
 }
 else {
 ldw = *n;
 lds2 = *n;
 }
 /* Set the workspace of the triangular matrix T to zero once such a */
 /* way every time T is generated the upper/lower portion will be always zero */
 claset_("A", &ldt, kd, &c_b1, &c_b1, &work[tpos], &ldt);
 if (upper) {
 i__1 = *n - *kd;
 i__2 = *kd;
 for (i__ = 1;
 i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
 i__ += i__2) {
 pn = *n - i__ - *kd + 1;
 /* Computing MIN */
 i__3 = *n - i__ - *kd + 1;
 pk = min(i__3,*kd);
 /* Compute the LQ factorization of the current block */
 cgelqf_(kd, &pn, &a[i__ + (i__ + *kd) * a_dim1], lda, &tau[i__], & work[s2pos], &ls2, &iinfo);
 /* Copy the upper portion of A into AB */
 i__3 = i__ + pk - 1;
 for (j = i__;
 j <= i__3;
 ++j) {
 /* Computing MIN */
 i__4 = *kd; i__5 = *n - j; // , expr subst  
 lk = min(i__4,i__5) + 1;
 i__4 = *ldab - 1;
 ccopy_(&lk, &a[j + j * a_dim1], lda, &ab[*kd + 1 + j * ab_dim1], &i__4);
 /* L20: */
 }
 claset_("Lower", &pk, &pk, &c_b1, &c_b2, &a[i__ + (i__ + *kd) * a_dim1], lda);
 /* Form the matrix T */
 clarft_("Forward", "Rowwise", &pn, &pk, &a[i__ + (i__ + *kd) * a_dim1], lda, &tau[i__], &work[tpos], &ldt);
 /* Compute W: */
 cgemm_("Conjugate", "No transpose", &pk, &pn, &pk, &c_b2, &work[ tpos], &ldt, &a[i__ + (i__ + *kd) * a_dim1], lda, &c_b1, & work[s2pos], &lds2);
 chemm_("Right", uplo, &pk, &pn, &c_b2, &a[i__ + *kd + (i__ + *kd) * a_dim1], lda, &work[s2pos], &lds2, &c_b1, &work[wpos], & ldw);
 cgemm_("No transpose", "Conjugate", &pk, &pk, &pn, &c_b2, &work[ wpos], &ldw, &work[s2pos], &lds2, &c_b1, &work[s1pos], & lds1);
 q__1.r = -.5f; q__1.i = -0.f; // , expr subst  
 cgemm_("No transpose", "No transpose", &pk, &pn, &pk, &q__1, & work[s1pos], &lds1, &a[i__ + (i__ + *kd) * a_dim1], lda, & c_b2, &work[wpos], &ldw);
 /* Update the unreduced submatrix A(i+kd:n,i+kd:n), using */
 /* an update of the form: A := A - V'*W - W'*V */
 q__1.r = -1.f; q__1.i = -0.f; // , expr subst  
 cher2k_(uplo, "Conjugate", &pn, &pk, &q__1, &a[i__ + (i__ + *kd) * a_dim1], lda, &work[wpos], &ldw, &c_b33, &a[i__ + *kd + ( i__ + *kd) * a_dim1], lda);
 /* L10: */
 }
 /* Copy the upper band to AB which is the band storage matrix */
 i__2 = *n;
 for (j = *n - *kd + 1;
 j <= i__2;
 ++j) {
 /* Computing MIN */
 i__1 = *kd; i__3 = *n - j; // , expr subst  
 lk = min(i__1,i__3) + 1;
 i__1 = *ldab - 1;
 ccopy_(&lk, &a[j + j * a_dim1], lda, &ab[*kd + 1 + j * ab_dim1], & i__1);
 /* L30: */
 }
 }
 else {
 /* Reduce the lower triangle of A to lower band matrix */
 i__2 = *n - *kd;
 i__1 = *kd;
 for (i__ = 1;
 i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
 i__ += i__1) {
 pn = *n - i__ - *kd + 1;
 /* Computing MIN */
 i__3 = *n - i__ - *kd + 1;
 pk = min(i__3,*kd);
 /* Compute the QR factorization of the current block */
 cgeqrf_(&pn, kd, &a[i__ + *kd + i__ * a_dim1], lda, &tau[i__], & work[s2pos], &ls2, &iinfo);
 /* Copy the upper portion of A into AB */
 i__3 = i__ + pk - 1;
 for (j = i__;
 j <= i__3;
 ++j) {
 /* Computing MIN */
 i__4 = *kd; i__5 = *n - j; // , expr subst  
 lk = min(i__4,i__5) + 1;
 ccopy_(&lk, &a[j + j * a_dim1], &c__1, &ab[j * ab_dim1 + 1], & c__1);
 /* L50: */
 }
 claset_("Upper", &pk, &pk, &c_b1, &c_b2, &a[i__ + *kd + i__ * a_dim1], lda);
 /* Form the matrix T */
 clarft_("Forward", "Columnwise", &pn, &pk, &a[i__ + *kd + i__ * a_dim1], lda, &tau[i__], &work[tpos], &ldt);
 /* Compute W: */
 cgemm_("No transpose", "No transpose", &pn, &pk, &pk, &c_b2, &a[ i__ + *kd + i__ * a_dim1], lda, &work[tpos], &ldt, &c_b1, &work[s2pos], &lds2);
 chemm_("Left", uplo, &pn, &pk, &c_b2, &a[i__ + *kd + (i__ + *kd) * a_dim1], lda, &work[s2pos], &lds2, &c_b1, &work[wpos], & ldw);
 cgemm_("Conjugate", "No transpose", &pk, &pk, &pn, &c_b2, &work[ s2pos], &lds2, &work[wpos], &ldw, &c_b1, &work[s1pos], & lds1);
 q__1.r = -.5f; q__1.i = -0.f; // , expr subst  
 cgemm_("No transpose", "No transpose", &pn, &pk, &pk, &q__1, &a[ i__ + *kd + i__ * a_dim1], lda, &work[s1pos], &lds1, & c_b2, &work[wpos], &ldw);
 /* Update the unreduced submatrix A(i+kd:n,i+kd:n), using */
 /* an update of the form: A := A - V*W' - W*V' */
 q__1.r = -1.f; q__1.i = -0.f; // , expr subst  
 cher2k_(uplo, "No transpose", &pn, &pk, &q__1, &a[i__ + *kd + i__ * a_dim1], lda, &work[wpos], &ldw, &c_b33, &a[i__ + *kd + (i__ + *kd) * a_dim1], lda);
 /* ================================================================== */
 /* RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED */
 /* DO 45 J = I, I+PK-1 */
 /* LK = MIN( KD, N-J ) + 1 */
 /* CALL CCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 ) */
 /* 45 CONTINUE */
 /* ================================================================== */
 /* L40: */
 }
 /* Copy the lower band to AB which is the band storage matrix */
 i__1 = *n;
 for (j = *n - *kd + 1;
 j <= i__1;
 ++j) {
 /* Computing MIN */
 i__2 = *kd; i__3 = *n - j; // , expr subst  
 lk = min(i__2,i__3) + 1;
 ccopy_(&lk, &a[j + j * a_dim1], &c__1, &ab[j * ab_dim1 + 1], & c__1);
 /* L60: */
 }
 }
 work[1].r = (real) lwmin; work[1].i = 0.f; // , expr subst  
 return 0;
 /* End of CHETRD_HE2HB */
 }
 /* chetrd_he2hb__ */
 