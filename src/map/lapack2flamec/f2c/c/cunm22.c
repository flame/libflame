/* ../netlib/v3.9.0/cunm22.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static complex c_b1 = {
1.f,0.f}
;
 /* > \brief \b CUNM22 multiplies a general matrix by a banded unitary matrix. */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download CUNM22 + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunm22. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunm22. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunm22. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE CUNM22( SIDE, TRANS, M, N, N1, N2, Q, LDQ, C, LDC, */
 /* $ WORK, LWORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER SIDE, TRANS */
 /* INTEGER M, N, N1, N2, LDQ, LDC, LWORK, INFO */
 /* .. */
 /* .. Array Arguments .. */
 /* COMPLEX Q( LDQ, * ), C( LDC, * ), WORK( * ) */
 /* .. */
 /* > \par Purpose */
 /* ============ */
 /* > */
 /* > \verbatim */
 /* > */
 /* > CUNM22 overwrites the general complex M-by-N matrix C with */
 /* > */
 /* > SIDE = 'L' SIDE = 'R' */
 /* > TRANS = 'N': Q * C C * Q */
 /* > TRANS = 'C': Q**H * C C * Q**H */
 /* > */
 /* > where Q is a complex unitary matrix of order NQ, with NQ = M if */
 /* > SIDE = 'L' and NQ = N if SIDE = 'R'. */
 /* > The unitary matrix Q processes a 2-by-2 block structure */
 /* > */
 /* > [ Q11 Q12 ] */
 /* > Q = [ ] */
 /* > [ Q21 Q22 ], */
 /* > */
 /* > where Q12 is an N1-by-N1 lower triangular matrix and Q21 is an */
 /* > N2-by-N2 upper triangular matrix. */
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
 /* > = 'N': apply Q (No transpose);
 */
 /* > = 'C': apply Q**H (Conjugate transpose). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] M */
 /* > \verbatim */
 /* > M is INTEGER */
 /* > The number of rows of the matrix C. M >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The number of columns of the matrix C. N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N1 */
 /* > \param[in] N2 */
 /* > \verbatim */
 /* > N1 is INTEGER */
 /* > N2 is INTEGER */
 /* > The dimension of Q12 and Q21, respectively. N1, N2 >= 0. */
 /* > The following requirement must be satisfied: */
 /* > N1 + N2 = M if SIDE = 'L' and N1 + N2 = N if SIDE = 'R'. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] Q */
 /* > \verbatim */
 /* > Q is COMPLEX array, dimension */
 /* > (LDQ,M) if SIDE = 'L' */
 /* > (LDQ,N) if SIDE = 'R' */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDQ */
 /* > \verbatim */
 /* > LDQ is INTEGER */
 /* > The leading dimension of the array Q. */
 /* > LDQ >= max(1,M) if SIDE = 'L';
 LDQ >= max(1,N) if SIDE = 'R'. */
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
 /* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
 /* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The dimension of the array WORK. */
 /* > If SIDE = 'L', LWORK >= max(1,N);
 */
 /* > if SIDE = 'R', LWORK >= max(1,M). */
 /* > For optimum performance LWORK >= M*N. */
 /* > */
 /* > If LWORK = -1, then a workspace query is assumed;
 the routine */
 /* > only calculates the optimal size of the WORK array, returns */
 /* > this value as the first entry of the WORK array, and no error */
 /* > message related to LWORK is issued by XERBLA. */
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
 /* > \date January 2015 */
 /* > \ingroup complexOTHERcomputational */
 /* ===================================================================== */
 /* Subroutine */
 int cunm22_(char *side, char *trans, integer *m, integer *n, integer *n1, integer *n2, complex *q, integer *ldq, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
 char buffer[256]; 
 snprintf(buffer, 256,"cunm22 inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", n1 %" FLA_IS ", n2 %" FLA_IS ", ldq %" FLA_IS ", ldc %" FLA_IS ", lwork %" FLA_IS "",*side, *trans, *m, *n, *n1, *n2, *ldq, *ldc, *lwork);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer q_dim1, q_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;
 complex q__1;
 /* Local variables */
 integer i__, nb, nq, nw, len;
 logical left;
 extern /* Subroutine */
 int cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
 extern logical lsame_(char *, char *);
 extern /* Subroutine */
 int ctrmm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *), clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), xerbla_(char *, integer *);
 logical notran;
 integer ldwork, lwkopt;
 logical lquery;
 /* -- LAPACK computational routine (version 3.7.1) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* January 2015 */
 /* .. Scalar Arguments .. */
 /* .. */
 /* .. Array Arguments .. */
 /* .. */
 /* ===================================================================== */
 /* .. Parameters .. */
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
 q_dim1 = *ldq;
 q_offset = 1 + q_dim1;
 q -= q_offset;
 c_dim1 = *ldc;
 c_offset = 1 + c_dim1;
 c__ -= c_offset;
 --work;
 /* Function Body */
 *info = 0;
 left = lsame_(side, "L");
 notran = lsame_(trans, "N");
 lquery = *lwork == -1;
 /* NQ is the order of Q;
 */
 /* NW is the minimum dimension of WORK. */
 if (left) {
 nq = *m;
 }
 else {
 nq = *n;
 }
 nw = nq;
 if (*n1 == 0 || *n2 == 0) {
 nw = 1;
 }
 if (! left && ! lsame_(side, "R")) {
 *info = -1;
 }
 else if (! lsame_(trans, "N") && ! lsame_(trans, "C")) {
 *info = -2;
 }
 else if (*m < 0) {
 *info = -3;
 }
 else if (*n < 0) {
 *info = -4;
 }
 else if (*n1 < 0 || *n1 + *n2 != nq) {
 *info = -5;
 }
 else if (*n2 < 0) {
 *info = -6;
 }
 else if (*ldq < max(1,nq)) {
 *info = -8;
 }
 else if (*ldc < max(1,*m)) {
 *info = -10;
 }
 else if (*lwork < nw && ! lquery) {
 *info = -12;
 }
 if (*info == 0) {
 lwkopt = *m * *n;
 q__1.r = (real) lwkopt; q__1.i = 0.f; // , expr subst  
 work[1].r = q__1.r; work[1].i = q__1.i; // , expr subst  
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("CUNM22", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 else if (lquery) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Quick return if possible */
 if (*m == 0 || *n == 0) {
 work[1].r = 1.f; work[1].i = 0.f; // , expr subst  
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Degenerate cases (N1 = 0 or N2 = 0) are handled using CTRMM. */
 if (*n1 == 0) {
 ctrmm_(side, "Upper", trans, "Non-Unit", m, n, &c_b1, &q[q_offset], ldq, &c__[c_offset], ldc);
 work[1].r = 1.f; work[1].i = 0.f; // , expr subst  
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 else if (*n2 == 0) {
 ctrmm_(side, "Lower", trans, "Non-Unit", m, n, &c_b1, &q[q_offset], ldq, &c__[c_offset], ldc);
 work[1].r = 1.f; work[1].i = 0.f; // , expr subst  
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Compute the largest chunk size available from the workspace. */
 /* Computing MAX */
 i__1 = 1; i__2 = min(*lwork,lwkopt) / nq; // , expr subst  
 nb = max(i__1,i__2);
 if (left) {
 if (notran) {
 i__1 = *n;
 i__2 = nb;
 for (i__ = 1;
 i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
 i__ += i__2) {
 /* Computing MIN */
 i__3 = nb; i__4 = *n - i__ + 1; // , expr subst  
 len = min(i__3,i__4);
 ldwork = *m;
 /* Multiply bottom part of C by Q12. */
 clacpy_("All", n1, &len, &c__[*n2 + 1 + i__ * c_dim1], ldc, & work[1], &ldwork);
 ctrmm_("Left", "Lower", "No Transpose", "Non-Unit", n1, &len, &c_b1, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[1], & ldwork);
 /* Multiply top part of C by Q11. */
 cgemm_("No Transpose", "No Transpose", n1, &len, n2, &c_b1, & q[q_offset], ldq, &c__[i__ * c_dim1 + 1], ldc, &c_b1, &work[1], &ldwork);
 /* Multiply top part of C by Q21. */
 clacpy_("All", n2, &len, &c__[i__ * c_dim1 + 1], ldc, &work[* n1 + 1], &ldwork);
 ctrmm_("Left", "Upper", "No Transpose", "Non-Unit", n2, &len, &c_b1, &q[*n1 + 1 + q_dim1], ldq, &work[*n1 + 1], & ldwork);
 /* Multiply bottom part of C by Q22. */
 cgemm_("No Transpose", "No Transpose", n2, &len, n1, &c_b1, & q[*n1 + 1 + (*n2 + 1) * q_dim1], ldq, &c__[*n2 + 1 + i__ * c_dim1], ldc, &c_b1, &work[*n1 + 1], &ldwork);
 /* Copy everything back. */
 clacpy_("All", m, &len, &work[1], &ldwork, &c__[i__ * c_dim1 + 1], ldc);
 }
 }
 else {
 i__2 = *n;
 i__1 = nb;
 for (i__ = 1;
 i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
 i__ += i__1) {
 /* Computing MIN */
 i__3 = nb; i__4 = *n - i__ + 1; // , expr subst  
 len = min(i__3,i__4);
 ldwork = *m;
 /* Multiply bottom part of C by Q21**H. */
 clacpy_("All", n2, &len, &c__[*n1 + 1 + i__ * c_dim1], ldc, & work[1], &ldwork);
 ctrmm_("Left", "Upper", "Conjugate", "Non-Unit", n2, &len, & c_b1, &q[*n1 + 1 + q_dim1], ldq, &work[1], &ldwork);
 /* Multiply top part of C by Q11**H. */
 cgemm_("Conjugate", "No Transpose", n2, &len, n1, &c_b1, &q[ q_offset], ldq, &c__[i__ * c_dim1 + 1], ldc, &c_b1, & work[1], &ldwork);
 /* Multiply top part of C by Q12**H. */
 clacpy_("All", n1, &len, &c__[i__ * c_dim1 + 1], ldc, &work[* n2 + 1], &ldwork);
 ctrmm_("Left", "Lower", "Conjugate", "Non-Unit", n1, &len, & c_b1, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[*n2 + 1], &ldwork);
 /* Multiply bottom part of C by Q22**H. */
 cgemm_("Conjugate", "No Transpose", n1, &len, n2, &c_b1, &q[* n1 + 1 + (*n2 + 1) * q_dim1], ldq, &c__[*n1 + 1 + i__ * c_dim1], ldc, &c_b1, &work[*n2 + 1], &ldwork);
 /* Copy everything back. */
 clacpy_("All", m, &len, &work[1], &ldwork, &c__[i__ * c_dim1 + 1], ldc);
 }
 }
 }
 else {
 if (notran) {
 i__1 = *m;
 i__2 = nb;
 for (i__ = 1;
 i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
 i__ += i__2) {
 /* Computing MIN */
 i__3 = nb; i__4 = *m - i__ + 1; // , expr subst  
 len = min(i__3,i__4);
 ldwork = len;
 /* Multiply right part of C by Q21. */
 clacpy_("All", &len, n2, &c__[i__ + (*n1 + 1) * c_dim1], ldc, &work[1], &ldwork);
 ctrmm_("Right", "Upper", "No Transpose", "Non-Unit", &len, n2, &c_b1, &q[*n1 + 1 + q_dim1], ldq, &work[1], &ldwork);
 /* Multiply left part of C by Q11. */
 cgemm_("No Transpose", "No Transpose", &len, n2, n1, &c_b1, & c__[i__ + c_dim1], ldc, &q[q_offset], ldq, &c_b1, & work[1], &ldwork);
 /* Multiply left part of C by Q12. */
 clacpy_("All", &len, n1, &c__[i__ + c_dim1], ldc, &work[*n2 * ldwork + 1], &ldwork);
 ctrmm_("Right", "Lower", "No Transpose", "Non-Unit", &len, n1, &c_b1, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[*n2 * ldwork + 1], &ldwork);
 /* Multiply right part of C by Q22. */
 cgemm_("No Transpose", "No Transpose", &len, n1, n2, &c_b1, & c__[i__ + (*n1 + 1) * c_dim1], ldc, &q[*n1 + 1 + (*n2 + 1) * q_dim1], ldq, &c_b1, &work[*n2 * ldwork + 1], & ldwork);
 /* Copy everything back. */
 clacpy_("All", &len, n, &work[1], &ldwork, &c__[i__ + c_dim1], ldc);
 }
 }
 else {
 i__2 = *m;
 i__1 = nb;
 for (i__ = 1;
 i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
 i__ += i__1) {
 /* Computing MIN */
 i__3 = nb; i__4 = *m - i__ + 1; // , expr subst  
 len = min(i__3,i__4);
 ldwork = len;
 /* Multiply right part of C by Q12**H. */
 clacpy_("All", &len, n1, &c__[i__ + (*n2 + 1) * c_dim1], ldc, &work[1], &ldwork);
 ctrmm_("Right", "Lower", "Conjugate", "Non-Unit", &len, n1, & c_b1, &q[(*n2 + 1) * q_dim1 + 1], ldq, &work[1], & ldwork);
 /* Multiply left part of C by Q11**H. */
 cgemm_("No Transpose", "Conjugate", &len, n1, n2, &c_b1, &c__[ i__ + c_dim1], ldc, &q[q_offset], ldq, &c_b1, &work[1] , &ldwork);
 /* Multiply left part of C by Q21**H. */
 clacpy_("All", &len, n2, &c__[i__ + c_dim1], ldc, &work[*n1 * ldwork + 1], &ldwork);
 ctrmm_("Right", "Upper", "Conjugate", "Non-Unit", &len, n2, & c_b1, &q[*n1 + 1 + q_dim1], ldq, &work[*n1 * ldwork + 1], &ldwork);
 /* Multiply right part of C by Q22**H. */
 cgemm_("No Transpose", "Conjugate", &len, n2, n1, &c_b1, &c__[ i__ + (*n2 + 1) * c_dim1], ldc, &q[*n1 + 1 + (*n2 + 1) * q_dim1], ldq, &c_b1, &work[*n1 * ldwork + 1], & ldwork);
 /* Copy everything back. */
 clacpy_("All", &len, n, &work[1], &ldwork, &c__[i__ + c_dim1], ldc);
 }
 }
 }
 q__1.r = (real) lwkopt; q__1.i = 0.f; // , expr subst  
 work[1].r = q__1.r; work[1].i = q__1.i; // , expr subst  
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of CUNM22 */
 }
 /* cunm22_ */
 