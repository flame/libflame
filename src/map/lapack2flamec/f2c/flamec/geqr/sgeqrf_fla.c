/* sgeqrf.f -- translated by f2c (version 20000121). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 static integer c_n1 = -1;
 static integer c__3 = 3;
 static integer c__2 = 2;
 /* > \brief \b SGEQRF */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download SGEQRF + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrf. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrf. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrf. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* INTEGER INFO, LDA, LWORK, M, N */
 /* .. */
 /* .. Array Arguments .. */
 /* REAL A( LDA, * ), TAU( * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > SGEQRF computes a QR factorization of a real M-by-N matrix A: */
 /* > */
 /* > A = Q * ( R ), */
 /* > ( 0 ) */
 /* > */
 /* > where: */
 /* > */
 /* > Q is a M-by-M orthogonal matrix;
 */
 /* > R is an upper-triangular N-by-N matrix;
 */
 /* > 0 is a (M-N)-by-N zero matrix, if M > N. */
 /* > */
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
 /* > \param[in,out] A */
 /* > \verbatim */
 /* > A is REAL array, dimension (LDA,N) */
 /* > On entry, the M-by-N matrix A. */
 /* > On exit, the elements on and above the diagonal of the array */
 /* > contain the min(M,N)-by-N upper trapezoidal matrix R (R is */
 /* > upper triangular if m >= n);
 the elements below the diagonal, */
 /* > with the array TAU, represent the orthogonal matrix Q as a */
 /* > product of min(m,n) elementary reflectors (see Further */
 /* > Details). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] TAU */
 /* > \verbatim */
 /* > TAU is REAL array, dimension (min(M,N)) */
 /* > The scalar factors of the elementary reflectors (see Further */
 /* > Details). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is REAL array, dimension (MAX(1,LWORK)) */
 /* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The dimension of the array WORK. LWORK >= max(1,N). */
 /* > For optimum performance LWORK >= N*NB, where NB is */
 /* > the optimal blocksize. */
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
 /* > \ingroup realGEcomputational */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > The matrix Q is represented as a product of elementary reflectors */
 /* > */
 /* > Q = H(1) H(2) . . . H(k), where k = min(m,n). */
 /* > */
 /* > Each H(i) has the form */
 /* > */
 /* > H(i) = I - tau * v * v**T */
 /* > */
 /* > where tau is a real scalar, and v is a real vector with */
 /* > v(1:i-1) = 0 and v(i) = 1;
 v(i+1:m) is stored on exit in A(i+1:m,i), */
 /* > and tau in TAU(i). */
 /* > \endverbatim */
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int sgeqrf_fla(integer *m, integer *n, real *a, integer *lda, real *tau, real *work, integer *lwork, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
 /* Local variables */
 integer i__, k, nbmin, iinfo;
 extern /* Subroutine */
 int sgeqr2_fla(integer *, integer *, real *, integer *, real *, real *, integer *);
 integer ib, nb, nx;
 extern /* Subroutine */
 int slarfb_(char *, char *, char *, char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *);
 extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
 extern /* Subroutine */
 int slarft_(char *, char *, integer *, integer *, real *, integer *, real *, real *, integer *);
 integer ldwork, lwkopt;
 logical lquery;
 integer iws;
 /* -- LAPACK computational routine -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* .. Scalar Arguments .. */
 /* .. */
 /* .. Array Arguments .. */
 /* .. */
 /* ===================================================================== */
 /* .. Local Scalars .. */
 /* .. */
 /* .. External Subroutines .. */
 /* .. */
 /* .. Intrinsic Functions .. */
 /* .. */
 /* .. External Functions .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Test the input arguments */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1 * 1;
 a -= a_offset;
 #if AOCL_FLA_PROGRESS_H
        AOCL_FLA_PROGRESS_VAR;
 #endif

 --tau;
 --work;
 /* Function Body */
 *info = 0;
 nb = ilaenv_(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1);
 lwkopt = *n * nb;
 work[1] = (real) lwkopt;
 lquery = *lwork == -1;
 if (*m < 0) {
 *info = -1;
 }
 else if (*n < 0) {
 *info = -2;
 }
 else if (*lda < max(1,*m)) {
 *info = -4;
 }
 else if (*lwork < max(1,*n) && ! lquery) {
 *info = -7;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("SGEQRF", &i__1);
 return 0;
 }
 else if (lquery) {
 return 0;
 }
 /* Quick return if possible */
 k = min(*m,*n);
 if (k == 0) {
 work[1] = 1.f;
 return 0;
 }
 nbmin = 2;
 nx = 0;
 iws = *n;

 #if AOCL_FLA_PROGRESS_H
        step_count =0;
    #ifndef FLA_ENABLE_WINDOWS_BUILD
        if(!aocl_fla_progress_ptr)
              aocl_fla_progress_ptr=aocl_fla_progress;
    #endif 
 #endif
 if (nb > 1 && nb < k) {
 /* Determine when to cross over from blocked to unblocked code. */
 /* Computing MAX */
 i__1 = 0; i__2 = ilaenv_(&c__3, "SGEQRF", " ", m, n, &c_n1, &c_n1); // , expr subst  
 nx = max(i__1,i__2);
 if (nx < k) {
 /* Determine if workspace is large enough for blocked code. */
 ldwork = *n;
 iws = ldwork * nb;
 if (*lwork < iws) {
 /* Not enough workspace to use optimal NB: reduce NB and */
 /* determine the minimum value of NB. */
 nb = *lwork / ldwork;
 /* Computing MAX */
 i__1 = 2; i__2 = ilaenv_(&c__2, "SGEQRF", " ", m, n, &c_n1, & c_n1); // , expr subst  
 nbmin = max(i__1,i__2);
 }
 }
 }
 if (nb >= nbmin && nb < k && nx < k) {
 /* Use blocked code initially */
 i__1 = k - nx;
 i__2 = nb;
 for (i__ = 1;
 i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
 i__ += i__2) {
 /* Computing MIN */
 i__3 = k - i__ + 1;
 ib = min(i__3,nb);
 /* Compute the QR factorization of the current block */
 /* A(i:m,i:i+ib-1) */
 i__3 = *m - i__ + 1;

  #if AOCL_FLA_PROGRESS_H

        if(aocl_fla_progress_ptr){
                step_count+=ib;
                AOCL_FLA_PROGRESS_FUNC_PTR("SGEQRF",6,&step_count,&thread_id,&total_threads);
         }

  #endif

 sgeqr2_fla(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[ 1], &iinfo);
 if (i__ + ib <= *n) {
 /* Form the triangular factor of the block reflector */
 /* H = H(i) H(i+1) . . . H(i+ib-1) */
 i__3 = *m - i__ + 1;
 slarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &ldwork);
 /* Apply H**T to A(i:m,i+ib:n) from the left */
 i__3 = *m - i__ + 1;
 i__4 = *n - i__ - ib + 1;
 slarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, & i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], & ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib + 1], &ldwork);
 }
 /* L10: */
 }
 }
 else {
 i__ = 1;
 }
 /* Use unblocked code to factor the last or only block. */
 if (i__ <= k) {
 i__2 = *m - i__ + 1;
 i__1 = *n - i__ + 1;

  #if AOCL_FLA_PROGRESS_H

         if(aocl_fla_progress_ptr){
                  step_count=k;
                   AOCL_FLA_PROGRESS_FUNC_PTR("SGEQRF",6,&step_count,&thread_id,&total_threads);
          }

  #endif

 sgeqr2_fla(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1] , &iinfo);
 }
 work[1] = (real) iws;
 return 0;
 /* End of SGEQRF */
 }
 /* sgeqrf_ */
 