/* ../netlib/v3.9.0/clamtsqr.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__0 = 0;
 /* > \brief \b CLAMTSQR */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE CLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, */
 /* $ LDT, C, LDC, WORK, LWORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER SIDE, TRANS */
 /* INTEGER INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC */
 /* .. */
 /* .. Array Arguments .. */
 /* COMPLEX A( LDA, * ), WORK( * ), C(LDC, * ), */
 /* $ T( LDT, * ) */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > CLAMTSQR overwrites the general complex M-by-N matrix C with */
 /* > */
 /* > */
 /* > SIDE = 'L' SIDE = 'R' */
 /* > TRANS = 'N': Q * C C * Q */
 /* > TRANS = 'C': Q**H * C C * Q**H */
 /* > where Q is a real orthogonal matrix defined as the product */
 /* > of blocked elementary reflectors computed by tall skinny */
 /* > QR factorization (CLATSQR) */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] SIDE */
 /* > \verbatim */
 /* > SIDE is CHARACTER*1 */
 /* > = 'L': apply Q or Q**H from the Left;
 */
 /* > = 'R': apply Q or Q**H from the Right. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] TRANS */
 /* > \verbatim */
 /* > TRANS is CHARACTER*1 */
 /* > = 'N': No transpose, apply Q;
 */
 /* > = 'C': Conjugate Transpose, apply Q**H. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] M */
 /* > \verbatim */
 /* > M is INTEGER */
 /* > The number of rows of the matrix A. M >=0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The number of columns of the matrix C. M >= N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] K */
 /* > \verbatim */
 /* > K is INTEGER */
 /* > The number of elementary reflectors whose product defines */
 /* > the matrix Q. */
 /* > N >= K >= 0;
 */
 /* > */
 /* > \endverbatim */
 /* > */
 /* > \param[in] MB */
 /* > \verbatim */
 /* > MB is INTEGER */
 /* > The block size to be used in the blocked QR. */
 /* > MB > N. (must be the same as DLATSQR) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] NB */
 /* > \verbatim */
 /* > NB is INTEGER */
 /* > The column block size to be used in the blocked QR. */
 /* > N >= NB >= 1. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] A */
 /* > \verbatim */
 /* > A is COMPLEX array, dimension (LDA,K) */
 /* > The i-th column must contain the vector which defines the */
 /* > blockedelementary reflector H(i), for i = 1,2,...,k, as */
 /* > returned by DLATSQR in the first k columns of */
 /* > its array argument A. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. */
 /* > If SIDE = 'L', LDA >= max(1,M);
 */
 /* > if SIDE = 'R', LDA >= max(1,N). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] T */
 /* > \verbatim */
 /* > T is COMPLEX array, dimension */
 /* > ( N * Number of blocks(CEIL(M-K/MB-K)), */
 /* > The blocked upper triangular block reflectors stored in compact form */
 /* > as a sequence of upper triangular blocks. See below */
 /* > for further details. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDT */
 /* > \verbatim */
 /* > LDT is INTEGER */
 /* > The leading dimension of the array T. LDT >= NB. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] C */
 /* > \verbatim */
 /* > C is COMPLEX array, dimension (LDC,N) */
 /* > On entry, the M-by-N matrix C. */
 /* > On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDC */
 /* > \verbatim */
 /* > LDC is INTEGER */
 /* > The leading dimension of the array C. LDC >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > (workspace) COMPLEX array, dimension (MAX(1,LWORK)) */
 /* > */
 /* > \endverbatim */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The dimension of the array WORK. */
 /* > */
 /* > If SIDE = 'L', LWORK >= max(1,N)*NB;
 */
 /* > if SIDE = 'R', LWORK >= max(1,MB)*NB. */
 /* > If LWORK = -1, then a workspace query is assumed;
 the routine */
 /* > only calculates the optimal size of the WORK array, returns */
 /* > this value as the first entry of the WORK array, and no error */
 /* > message related to LWORK is issued by XERBLA. */
 /* > */
 /* > \endverbatim */
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
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > Tall-Skinny QR (TSQR) performs QR by a sequence of orthogonal transformations, */
 /* > representing Q as a product of other orthogonal matrices */
 /* > Q = Q(1) * Q(2) * . . . * Q(k) */
 /* > where each Q(i) zeros out subdiagonal entries of a block of MB rows of A: */
 /* > Q(1) zeros out the subdiagonal entries of rows 1:MB of A */
 /* > Q(2) zeros out the bottom MB-N rows of rows [1:N,MB+1:2*MB-N] of A */
 /* > Q(3) zeros out the bottom MB-N rows of rows [1:N,2*MB-N+1:3*MB-2*N] of A */
 /* > . . . */
 /* > */
 /* > Q(1) is computed by GEQRT, which represents Q(1) by Householder vectors */
 /* > stored under the diagonal of rows 1:MB of A, and by upper triangular */
 /* > block reflectors, stored in array T(1:LDT,1:N). */
 /* > For more information see Further Details in GEQRT. */
 /* > */
 /* > Q(i) for i>1 is computed by TPQRT, which represents Q(i) by Householder vectors */
 /* > stored in rows [(i-1)*(MB-N)+N+1:i*(MB-N)+N] of A, and by upper triangular */
 /* > block reflectors, stored in array T(1:LDT,(i-1)*N+1:i*N). */
 /* > The last Q(k) may use fewer rows. */
 /* > For more information see Further Details in TPQRT. */
 /* > */
 /* > For more details of the overall algorithm, see the description of */
 /* > Sequential TSQR in Section 2.2 of [1]. */
 /* > */
 /* > [1] “Communication-Optimal Parallel and Sequential QR and LU Factorizations, */
 /* > J. Demmel, L. Grigori, M. Hoemmen, J. Langou, */
 /* > SIAM J. Sci. Comput, vol. 34, no. 1, 2012 */
 /* > \endverbatim */
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int clamtsqr_(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, complex *a, integer *lda, complex *t, integer *ldt, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
 char buffer[256]; 
#if FLA_ENABLE_ILP64 
 snprintf(buffer, 256,"clamtsqr inputs: side %c, trans %c, m %lld, n %lld, k %lld, mb %lld, nb %lld, lda %lld, ldt %lld, ldc %lld, lwork %lld",*side, *trans, *m, *n, *k, *mb, *nb, *lda, *ldt, *ldc, *lwork);
#else 
 snprintf(buffer, 256,"clamtsqr inputs: side %c, trans %c, m %d, n %d, k %d, mb %d, nb %d, lda %d, ldt %d, ldc %d, lwork %d",*side, *trans, *m, *n, *k, *mb, *nb, *lda, *ldt, *ldc, *lwork);
#endif
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, i__3;
 /* Local variables */
 integer i__, ii, kk, lw, ctr;
 logical left, tran;
 extern logical lsame_(char *, char *);
 logical right;
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 logical notran, lquery;
 extern /* Subroutine */
 int cgemqrt_(char *, char *, integer *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *), ctpmqrt_(char *, char *, integer *, integer *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *);
 /* -- LAPACK computational routine (version 3.7.1) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* June 2017 */
 /* .. Scalar Arguments .. */
 /* .. */
 /* .. Array Arguments .. */
 /* .. */
 /* ===================================================================== */
 /* .. */
 /* .. Local Scalars .. */
 /* .. */
 /* .. External Functions .. */
 /* .. External Subroutines .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Test the input arguments */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 t_dim1 = *ldt;
 t_offset = 1 + t_dim1;
 t -= t_offset;
 c_dim1 = *ldc;
 c_offset = 1 + c_dim1;
 c__ -= c_offset;
 --work;
 /* Function Body */
 lquery = *lwork < 0;
 notran = lsame_(trans, "N");
 tran = lsame_(trans, "C");
 left = lsame_(side, "L");
 right = lsame_(side, "R");
 if (left) {
 lw = *n * *nb;
 }
 else {
 lw = *m * *nb;
 }
 *info = 0;
 if (! left && ! right) {
 *info = -1;
 }
 else if (! tran && ! notran) {
 *info = -2;
 }
 else if (*m < 0) {
 *info = -3;
 }
 else if (*n < 0) {
 *info = -4;
 }
 else if (*k < 0) {
 *info = -5;
 }
 else if (*lda < max(1,*k)) {
 *info = -9;
 }
 else if (*ldt < max(1,*nb)) {
 *info = -11;
 }
 else if (*ldc < max(1,*m)) {
 *info = -13;
 }
 else if (*lwork < max(1,lw) && ! lquery) {
 *info = -15;
 }
 /* Determine the block size if it is tall skinny or short and wide */
 if (*info == 0) {
 work[1].r = (real) lw; work[1].i = 0.f; // , expr subst  
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("CLAMTSQR", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 else if (lquery) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Quick return if possible */
 /* Computing MIN */
 i__1 = min(*m,*n);
 if (min(i__1,*k) == 0) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Computing MAX */
 i__1 = max(*m,*n);
 if (*mb <= *k || *mb >= max(i__1,*k)) {
 cgemqrt_(side, trans, m, n, k, nb, &a[a_offset], lda, &t[t_offset], ldt, &c__[c_offset], ldc, &work[1], info);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 if (left && notran) {
 /* Multiply Q to the last block of C */
 kk = (*m - *k) % (*mb - *k);
 ctr = (*m - *k) / (*mb - *k);
 if (kk > 0) {
 ii = *m - kk + 1;
 ctpmqrt_("L", "N", &kk, n, k, &c__0, nb, &a[ii + a_dim1], lda, &t[ (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[ii + c_dim1], ldc, &work[1], info);
 }
 else {
 ii = *m + 1;
 }
 i__1 = *mb + 1;
 i__2 = -(*mb - *k);
 for (i__ = ii - (*mb - *k);
 i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
 i__ += i__2) {
 /* Multiply Q to the current block of C (I:I+MB,1:N) */
 --ctr;
 i__3 = *mb - *k;
 ctpmqrt_("L", "N", &i__3, n, k, &c__0, nb, &a[i__ + a_dim1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info);
 }
 /* Multiply Q to the first block of C (1:MB,1:N) */
 cgemqrt_("L", "N", mb, n, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &c__[c_dim1 + 1], ldc, &work[1], info);
 }
 else if (left && tran) {
 /* Multiply Q to the first block of C */
 kk = (*m - *k) % (*mb - *k);
 ii = *m - kk + 1;
 ctr = 1;
 cgemqrt_("L", "C", mb, n, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &c__[c_dim1 + 1], ldc, &work[1], info);
 i__2 = ii - *mb + *k;
 i__1 = *mb - *k;
 for (i__ = *mb + 1;
 i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
 i__ += i__1) {
 /* Multiply Q to the current block of C (I:I+MB,1:N) */
 i__3 = *mb - *k;
 ctpmqrt_("L", "C", &i__3, n, k, &c__0, nb, &a[i__ + a_dim1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info);
 ++ctr;
 }
 if (ii <= *m) {
 /* Multiply Q to the last block of C */
 ctpmqrt_("L", "C", &kk, n, k, &c__0, nb, &a[ii + a_dim1], lda, &t[ (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[ii + c_dim1], ldc, &work[1], info);
 }
 }
 else if (right && tran) {
 /* Multiply Q to the last block of C */
 kk = (*n - *k) % (*mb - *k);
 ctr = (*n - *k) / (*mb - *k);
 if (kk > 0) {
 ii = *n - kk + 1;
 ctpmqrt_("R", "C", m, &kk, k, &c__0, nb, &a[ii + a_dim1], lda, &t[ (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info);
 }
 else {
 ii = *n + 1;
 }
 i__1 = *mb + 1;
 i__2 = -(*mb - *k);
 for (i__ = ii - (*mb - *k);
 i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
 i__ += i__2) {
 /* Multiply Q to the current block of C (1:M,I:I+MB) */
 --ctr;
 i__3 = *mb - *k;
 ctpmqrt_("R", "C", m, &i__3, k, &c__0, nb, &a[i__ + a_dim1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info);
 }
 /* Multiply Q to the first block of C (1:M,1:MB) */
 cgemqrt_("R", "C", m, mb, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &c__[c_dim1 + 1], ldc, &work[1], info);
 }
 else if (right && notran) {
 /* Multiply Q to the first block of C */
 kk = (*n - *k) % (*mb - *k);
 ii = *n - kk + 1;
 ctr = 1;
 cgemqrt_("R", "N", m, mb, k, nb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &c__[c_dim1 + 1], ldc, &work[1], info);
 i__2 = ii - *mb + *k;
 i__1 = *mb - *k;
 for (i__ = *mb + 1;
 i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
 i__ += i__1) {
 /* Multiply Q to the current block of C (1:M,I:I+MB) */
 i__3 = *mb - *k;
 ctpmqrt_("R", "N", m, &i__3, k, &c__0, nb, &a[i__ + a_dim1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info);
 ++ctr;
 }
 if (ii <= *n) {
 /* Multiply Q to the last block of C */
 ctpmqrt_("R", "N", m, &kk, k, &c__0, nb, &a[ii + a_dim1], lda, &t[ (ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info);
 }
 }
 work[1].r = (real) lw; work[1].i = 0.f; // , expr subst  
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of CLAMTSQR */
 }
 /* clamtsqr_ */
 