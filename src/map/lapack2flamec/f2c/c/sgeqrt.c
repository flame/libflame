/* ../netlib/v3.9.0/sgeqrt.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* > \brief \b SGEQRT */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download SGEQRT + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrt. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrt. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrt. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE SGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* INTEGER INFO, LDA, LDT, M, N, NB */
 /* .. */
 /* .. Array Arguments .. */
 /* REAL A( LDA, * ), T( LDT, * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > SGEQRT computes a blocked QR factorization of a real M-by-N matrix A */
 /* > using the compact WY representation of Q. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] M */
 /* > \verbatim */
 /* > M is INTEGER */
 /* > The number of rows of the matrix A. M >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The number of columns of the matrix A. N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] NB */
 /* > \verbatim */
 /* > NB is INTEGER */
 /* > The block size to be used in the blocked QR. MIN(M,N) >= NB >= 1. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] A */
 /* > \verbatim */
 /* > A is REAL array, dimension (LDA,N) */
 /* > On entry, the M-by-N matrix A. */
 /* > On exit, the elements on and above the diagonal of the array */
 /* > contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
 /* > upper triangular if M >= N);
 the elements below the diagonal */
 /* > are the columns of V. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] T */
 /* > \verbatim */
 /* > T is REAL array, dimension (LDT,MIN(M,N)) */
 /* > The upper triangular block reflectors stored in compact form */
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
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is REAL array, dimension (NB*N) */
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
 /* > \ingroup realGEcomputational */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > The matrix V stores the elementary reflectors H(i) in the i-th column */
 /* > below the diagonal. For example, if M=5 and N=3, the matrix V is */
 /* > */
 /* > V = ( 1 ) */
 /* > ( v1 1 ) */
 /* > ( v1 v2 1 ) */
 /* > ( v1 v2 v3 ) */
 /* > ( v1 v2 v3 ) */
 /* > */
 /* > where the vi's represent the vectors which define H(i), which are returned */
 /* > in the matrix A. The 1's along the diagonal of V are not stored in A. */
 /* > */
 /* > Let K=MIN(M,N). The number of blocks is B = ceiling(K/NB), where each */
 /* > block is of order NB except for the last block, which is of order */
 /* > IB = K - (B-1)*NB. For each of the B blocks, a upper triangular block */
 /* > reflector factor is computed: T1, T2, ..., TB. The NB-by-NB (and IB-by-IB */
 /* > for the last block) T's are stored in the NB-by-K matrix T as */
 /* > */
 /* > T = (T1 T2 ... TB). */
 /* > \endverbatim */
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int sgeqrt_(integer *m, integer *n, integer *nb, real *a, integer *lda, real *t, integer *ldt, real *work, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;
 /* Local variables */
 integer i__, k, ib, iinfo;
 extern /* Subroutine */
 int slarfb_(char *, char *, char *, char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *), sgeqrt2_( integer *, integer *, real *, integer *, real *, integer *, integer *), sgeqrt3_(integer *, integer *, real *, integer *, real *, integer *, integer *);
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
 --work;
 /* Function Body */
 *info = 0;
 if (*m < 0) {
 *info = -1;
 }
 else if (*n < 0) {
 *info = -2;
 }
 else if (*nb < 1 || *nb > min(*m,*n) && min(*m,*n) > 0) {
 *info = -3;
 }
 else if (*lda < max(1,*m)) {
 *info = -5;
 }
 else if (*ldt < *nb) {
 *info = -7;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("SGEQRT", &i__1);
 return 0;
 }
 /* Quick return if possible */
 k = min(*m,*n);
 if (k == 0) {
 return 0;
 }
 /* Blocked loop of length K */
 i__1 = k;
 i__2 = *nb;
 for (i__ = 1;
 i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
 i__ += i__2) {
 /* Computing MIN */
 i__3 = k - i__ + 1;
 ib = min(i__3,*nb);
 /* Compute the QR factorization of the current block A(I:M,I:I+IB-1) */
 if (TRUE_) {
 i__3 = *m - i__ + 1;
 sgeqrt3_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &iinfo);
 }
 else {
 i__3 = *m - i__ + 1;
 sgeqrt2_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &iinfo);
 }
 if (i__ + ib <= *n) {
 /* Update by applying H**T to A(I:M,I+IB:N) from the left */
 i__3 = *m - i__ + 1;
 i__4 = *n - i__ - ib + 1;
 i__5 = *n - i__ - ib + 1;
 slarfb_("L", "T", "F", "C", &i__3, &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + ib) * a_dim1], lda, &work[1], &i__5);
 }
 }
 return 0;
 /* End of SGEQRT */
 }
 /* sgeqrt_ */
 