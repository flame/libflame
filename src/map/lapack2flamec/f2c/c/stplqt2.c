/* ../netlib/v3.9.0/stplqt2.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static real c_b4 = 1.f;
 static real c_b10 = 0.f;
 /* > \brief \b STPLQT2 computes a LQ factorization of a real or complex "triangular-pentagonal" matrix, which is composed of a triangular block and a pentagonal block, using the compact WY representation for Q. */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download STPLQT2 + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stplqt2 .f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stplqt2 .f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stplqt2 .f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE STPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO ) */
 /* .. Scalar Arguments .. */
 /* INTEGER INFO, LDA, LDB, LDT, N, M, L */
 /* .. */
 /* .. Array Arguments .. */
 /* REAL A( LDA, * ), B( LDB, * ), T( LDT, * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > STPLQT2 computes a LQ a factorization of a real "triangular-pentagonal" */
 /* > matrix C, which is composed of a triangular block A and pentagonal block B, */
 /* > using the compact WY representation for Q. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] M */
 /* > \verbatim */
 /* > M is INTEGER */
 /* > The total number of rows of the matrix B. */
 /* > M >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The number of columns of the matrix B, and the order of */
 /* > the triangular matrix A. */
 /* > N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] L */
 /* > \verbatim */
 /* > L is INTEGER */
 /* > The number of rows of the lower trapezoidal part of B. */
 /* > MIN(M,N) >= L >= 0. See Further Details. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] A */
 /* > \verbatim */
 /* > A is REAL array, dimension (LDA,M) */
 /* > On entry, the lower triangular M-by-M matrix A. */
 /* > On exit, the elements on and below the diagonal of the array */
 /* > contain the lower triangular matrix L. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] B */
 /* > \verbatim */
 /* > B is REAL array, dimension (LDB,N) */
 /* > On entry, the pentagonal M-by-N matrix B. The first N-L columns */
 /* > are rectangular, and the last L columns are lower trapezoidal. */
 /* > On exit, B contains the pentagonal matrix V. See Further Details. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDB */
 /* > \verbatim */
 /* > LDB is INTEGER */
 /* > The leading dimension of the array B. LDB >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] T */
 /* > \verbatim */
 /* > T is REAL array, dimension (LDT,M) */
 /* > The N-by-N upper triangular factor T of the block reflector. */
 /* > See Further Details. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDT */
 /* > \verbatim */
 /* > LDT is INTEGER */
 /* > The leading dimension of the array T. LDT >= max(1,M) */
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
 /* > \ingroup doubleOTHERcomputational */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > The input matrix C is a M-by-(M+N) matrix */
 /* > */
 /* > C = [ A ][ B ] */
 /* > */
 /* > */
 /* > where A is an lower triangular M-by-M matrix, and B is M-by-N pentagonal */
 /* > matrix consisting of a M-by-(N-L) rectangular matrix B1 left of a M-by-L */
 /* > upper trapezoidal matrix B2: */
 /* > */
 /* > B = [ B1 ][ B2 ] */
 /* > [ B1 ] <- M-by-(N-L) rectangular */
 /* > [ B2 ] <- M-by-L lower trapezoidal. */
 /* > */
 /* > The lower trapezoidal matrix B2 consists of the first L columns of a */
 /* > N-by-N lower triangular matrix, where 0 <= L <= MIN(M,N). If L=0, */
 /* > B is rectangular M-by-N;
 if M=L=N, B is lower triangular. */
 /* > */
 /* > The matrix W stores the elementary reflectors H(i) in the i-th row */
 /* > above the diagonal (of A) in the M-by-(M+N) input matrix C */
 /* > */
 /* > C = [ A ][ B ] */
 /* > [ A ] <- lower triangular M-by-M */
 /* > [ B ] <- M-by-N pentagonal */
 /* > */
 /* > so that W can be represented as */
 /* > */
 /* > W = [ I ][ V ] */
 /* > [ I ] <- identity, M-by-M */
 /* > [ V ] <- M-by-N, same form as B. */
 /* > */
 /* > Thus, all of information needed for W is contained on exit in B, which */
 /* > we call V above. Note that V has the same form as B;
 that is, */
 /* > */
 /* > W = [ V1 ][ V2 ] */
 /* > [ V1 ] <- M-by-(N-L) rectangular */
 /* > [ V2 ] <- M-by-L lower trapezoidal. */
 /* > */
 /* > The rows of V represent the vectors which define the H(i)'s. */
 /* > The (M+N)-by-(M+N) block reflector H is then given by */
 /* > */
 /* > H = I - W**T * T * W */
 /* > */
 /* > where W^H is the conjugate transpose of W and T is the upper triangular */
 /* > factor of the block reflector. */
 /* > \endverbatim */
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int stplqt2_(integer *m, integer *n, integer *l, real *a, integer *lda, real *b, integer *ldb, real *t, integer *ldt, integer * info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
 char buffer[256];
 snprintf(buffer, 256,"stplqt2 inputs: m %" FLA_IS ", n %" FLA_IS ", l %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS ", ldt %" FLA_IS "",*m, *n, *l, *lda, *ldb, *ldt);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, i__3;
 /* Local variables */
 integer i__, j, p, mp, np;
 extern /* Subroutine */
 int sger_(integer *, integer *, real *, real *, integer *, real *, integer *, real *, integer *);
 real alpha;
 extern /* Subroutine */
 int sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), strmv_(char *, char *, char *, integer *, real *, integer *, real *, integer *), xerbla_( char *, integer *), slarfg_(integer *, real *, real *, integer *, real *);
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
 b_dim1 = *ldb;
 b_offset = 1 + b_dim1;
 b -= b_offset;
 t_dim1 = *ldt;
 t_offset = 1 + t_dim1;
 t -= t_offset;
 /* Function Body */
 *info = 0;
 if (*m < 0) {
 *info = -1;
 }
 else if (*n < 0) {
 *info = -2;
 }
 else if (*l < 0 || *l > min(*m,*n)) {
 *info = -3;
 }
 else if (*lda < max(1,*m)) {
 *info = -5;
 }
 else if (*ldb < max(1,*m)) {
 *info = -7;
 }
 else if (*ldt < max(1,*m)) {
 *info = -9;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("STPLQT2", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Quick return if possible */
 if (*n == 0 || *m == 0) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 i__1 = *m;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 /* Generate elementary reflector H(I) to annihilate B(I,:) */
 p = *n - *l + min(*l,i__);
 i__2 = p + 1;
 slarfg_(&i__2, &a[i__ + i__ * a_dim1], &b[i__ + b_dim1], ldb, &t[i__ * t_dim1 + 1]);
 if (i__ < *m) {
 /* W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)] */
 i__2 = *m - i__;
 for (j = 1;
 j <= i__2;
 ++j) {
 t[*m + j * t_dim1] = a[i__ + j + i__ * a_dim1];
 }
 i__2 = *m - i__;
 sgemv_("N", &i__2, &p, &c_b4, &b[i__ + 1 + b_dim1], ldb, &b[i__ + b_dim1], ldb, &c_b4, &t[*m + t_dim1], ldt);
 /* C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H */
 alpha = -t[i__ * t_dim1 + 1];
 i__2 = *m - i__;
 for (j = 1;
 j <= i__2;
 ++j) {
 a[i__ + j + i__ * a_dim1] += alpha * t[*m + j * t_dim1];
 }
 i__2 = *m - i__;
 sger_(&i__2, &p, &alpha, &t[*m + t_dim1], ldt, &b[i__ + b_dim1], ldb, &b[i__ + 1 + b_dim1], ldb);
 }
 }
 i__1 = *m;
 for (i__ = 2;
 i__ <= i__1;
 ++i__) {
 /* T(I,1:I-1) := C(I:I-1,1:N) * (alpha * C(I,I:N)^H) */
 alpha = -t[i__ * t_dim1 + 1];
 i__2 = i__ - 1;
 for (j = 1;
 j <= i__2;
 ++j) {
 t[i__ + j * t_dim1] = 0.f;
 }
 /* Computing MIN */
 i__2 = i__ - 1;
 p = min(i__2,*l);
 /* Computing MIN */
 i__2 = *n - *l + 1;
 np = min(i__2,*n);
 /* Computing MIN */
 i__2 = p + 1;
 mp = min(i__2,*m);
 /* Triangular part of B2 */
 i__2 = p;
 for (j = 1;
 j <= i__2;
 ++j) {
 t[i__ + j * t_dim1] = alpha * b[i__ + (*n - *l + j) * b_dim1];
 }
 strmv_("L", "N", "N", &p, &b[np * b_dim1 + 1], ldb, &t[i__ + t_dim1], ldt);
 /* Rectangular part of B2 */
 i__2 = i__ - 1 - p;
 sgemv_("N", &i__2, l, &alpha, &b[mp + np * b_dim1], ldb, &b[i__ + np * b_dim1], ldb, &c_b10, &t[i__ + mp * t_dim1], ldt);
 /* B1 */
 i__2 = i__ - 1;
 i__3 = *n - *l;
 sgemv_("N", &i__2, &i__3, &alpha, &b[b_offset], ldb, &b[i__ + b_dim1], ldb, &c_b4, &t[i__ + t_dim1], ldt);
 /* T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1) */
 i__2 = i__ - 1;
 strmv_("L", "T", "N", &i__2, &t[t_offset], ldt, &t[i__ + t_dim1], ldt);
 /* T(I,I) = tau(I) */
 t[i__ + i__ * t_dim1] = t[i__ * t_dim1 + 1];
 t[i__ * t_dim1 + 1] = 0.f;
 }
 i__1 = *m;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = *m;
 for (j = i__ + 1;
 j <= i__2;
 ++j) {
 t[i__ + j * t_dim1] = t[j + i__ * t_dim1];
 t[j + i__ * t_dim1] = 0.f;
 }
 }
 /* End of STPLQT2 */
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* stplqt2_ */
 
