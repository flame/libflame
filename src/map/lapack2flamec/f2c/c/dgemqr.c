/* ../netlib/v3.9.0/dgemqr.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* > \brief \b DGEMQR */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE DGEMQR( SIDE, TRANS, M, N, K, A, LDA, T, */
 /* $ TSIZE, C, LDC, WORK, LWORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER SIDE, TRANS */
 /* INTEGER INFO, LDA, M, N, K, LDT, TSIZE, LWORK, LDC */
 /* .. */
 /* .. Array Arguments .. */
 /* DOUBLE PRECISION A( LDA, * ), T( * ), C( LDC, * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DGEMQR overwrites the general real M-by-N matrix C with */
 /* > */
 /* > SIDE = 'L' SIDE = 'R' */
 /* > TRANS = 'N': Q * C C * Q */
 /* > TRANS = 'T': Q**T * C C * Q**T */
 /* > */
 /* > where Q is a real orthogonal matrix defined as the product */
 /* > of blocked elementary reflectors computed by tall skinny */
 /* > QR factorization (DGEQR) */
 /* > */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] SIDE */
 /* > \verbatim */
 /* > SIDE is CHARACTER*1 */
 /* > = 'L': apply Q or Q**T from the Left;
 */
 /* > = 'R': apply Q or Q**T from the Right. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] TRANS */
 /* > \verbatim */
 /* > TRANS is CHARACTER*1 */
 /* > = 'N': No transpose, apply Q;
 */
 /* > = 'T': Transpose, apply Q**T. */
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
 /* > The number of columns of the matrix C. N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] K */
 /* > \verbatim */
 /* > K is INTEGER */
 /* > The number of elementary reflectors whose product defines */
 /* > the matrix Q. */
 /* > If SIDE = 'L', M >= K >= 0;
 */
 /* > if SIDE = 'R', N >= K >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] A */
 /* > \verbatim */
 /* > A is DOUBLE PRECISION array, dimension (LDA,K) */
 /* > Part of the data structure to represent Q as returned by DGEQR. */
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
 /* > T is DOUBLE PRECISION array, dimension (MAX(5,TSIZE)). */
 /* > Part of the data structure to represent Q as returned by DGEQR. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] TSIZE */
 /* > \verbatim */
 /* > TSIZE is INTEGER */
 /* > The dimension of the array T. TSIZE >= 5. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] C */
 /* > \verbatim */
 /* > C is DOUBLE PRECISION array, dimension (LDC,N) */
 /* > On entry, the M-by-N matrix C. */
 /* > On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */
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
 /* > (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The dimension of the array WORK. */
 /* > If LWORK = -1, then a workspace query is assumed. The routine */
 /* > only calculates the size of the WORK array, returns this */
 /* > value as WORK(1), and no error message related to WORK */
 /* > is issued by XERBLA. */
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
 /* > \par Further Details */
 /* ==================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > These details are particular for this LAPACK implementation. Users should not */
 /* > take them for granted. These details may change in the future, and are not likely */
 /* > true for another LAPACK implementation. These details are relevant if one wants */
 /* > to try to understand the code. They are not part of the interface. */
 /* > */
 /* > In this version, */
 /* > */
 /* > T(2): row block size (MB) */
 /* > T(3): column block size (NB) */
 /* > T(6:TSIZE): data structure needed for Q, computed by */
 /* > DLATSQR or DGEQRT */
 /* > */
 /* > Depending on the matrix dimensions M and N, and row and column */
 /* > block sizes MB and NB returned by ILAENV, DGEQR will use either */
 /* > DLATSQR (if the matrix is tall-and-skinny) or DGEQRT to compute */
 /* > the QR factorization. */
 /* > This version of DGEMQR will use either DLAMTSQR or DGEMQRT to */
 /* > multiply matrix Q by another matrix. */
 /* > Further Details in DLATMSQR or DGEMQRT. */
 /* > */
 /* > \endverbatim */
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int dgemqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *t, integer * tsize, doublereal *c__, integer *ldc, doublereal *work, integer * lwork, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, c_dim1, c_offset, i__1;
 /* Local variables */
 extern /* Subroutine */
 int dlamtsqr_(char *, char *, integer *, integer * , integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *);
 integer mb, nb, mn, lw;
 logical left, tran;
 extern logical lsame_(char *, char *);
 logical right;
 integer nblcks;
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 logical notran, lquery;
 extern /* Subroutine */
 int dgemqrt_(char *, char *, integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
 /* -- LAPACK computational routine (version 3.7.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* December 2016 */
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
 /* Test the input arguments */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 --t;
 c_dim1 = *ldc;
 c_offset = 1 + c_dim1;
 c__ -= c_offset;
 --work;
 /* Function Body */
 lquery = *lwork == -1;
 notran = lsame_(trans, "N");
 tran = lsame_(trans, "T");
 left = lsame_(side, "L");
 right = lsame_(side, "R");
 mb = (integer) t[2];
 nb = (integer) t[3];
 if (left) {
 lw = *n * nb;
 mn = *m;
 }
 else {
 lw = mb * nb;
 mn = *n;
 }
 if (mb > *k && mn > *k) {
 if ((mn - *k) % (mb - *k) == 0) {
 nblcks = (mn - *k) / (mb - *k);
 }
 else {
 nblcks = (mn - *k) / (mb - *k) + 1;
 }
 }
 else {
 nblcks = 1;
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
 else if (*k < 0 || *k > mn) {
 *info = -5;
 }
 else if (*lda < max(1,mn)) {
 *info = -7;
 }
 else if (*tsize < 5) {
 *info = -9;
 }
 else if (*ldc < max(1,*m)) {
 *info = -11;
 }
 else if (*lwork < max(1,lw) && ! lquery) {
 *info = -13;
 }
 if (*info == 0) {
 work[1] = (doublereal) lw;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("DGEMQR", &i__1);
 return 0;
 }
 else if (lquery) {
 return 0;
 }
 /* Quick return if possible */
 /* Computing MIN */
 i__1 = min(*m,*n);
 if (min(i__1,*k) == 0) {
 return 0;
 }
 /* Computing MAX */
 i__1 = max(*m,*n);
 if (left && *m <= *k || right && *n <= *k || mb <= *k || mb >= max(i__1,* k)) {
 dgemqrt_(side, trans, m, n, k, &nb, &a[a_offset], lda, &t[6], &nb, & c__[c_offset], ldc, &work[1], info);
 }
 else {
 dlamtsqr_(side, trans, m, n, k, &mb, &nb, &a[a_offset], lda, &t[6], & nb, &c__[c_offset], ldc, &work[1], lwork, info);
 }
 work[1] = (doublereal) lw;
 return 0;
 /* End of DGEMQR */
 }
 /* dgemqr_ */
 