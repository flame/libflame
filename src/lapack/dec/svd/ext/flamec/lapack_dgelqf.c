/* dgelqf.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLAME.h"
 #include "FLA_f2c.h" /* Table of constant values */
#ifdef FLA_ENABLE_AMD_OPT
extern int dgeqrf_fla(integer *m, integer *n,
                      doublereal *a, integer * lda,
                      doublereal *tau,
                      doublereal *work, integer *lwork,
                      integer *info);
int fla_dgelqf_tran(integer *m, integer *n,
                    doublereal *a, integer * lda,
                    doublereal *tau,
                    doublereal *work, integer *lwork,
                    integer *info);
/* Transpose copy from matrix 'a' to 'b'.
 * Memory locations of a and b are assumed to be different
 * without any overlapping memory
 * */
static void fla_dtranspose(integer *m, integer *n,
                           doublereal *a, integer *lda,
                           doublereal *b, integer *ldb)
{
    integer i, j;

    /* Offset adjustments */
    a -= (1 + *lda);
    b -= (1 + *ldb);

    /* Do the transpose copy */
    for (i = 1; i <= *n; i++)
    {
        for (j = 1; j <= *m; j++)
        {
            b[i + j * *ldb] = a[i * *lda + j];
        }
    }
}
#endif
 static integer c__1 = 1;
 static integer c_n1 = -1;
 static integer c__3 = 3;
 static integer c__2 = 2;
 /* > \brief \b DGELQF */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DGELQF + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelqf. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelqf. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelqf. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* INTEGER INFO, LDA, LWORK, M, N */
 /* .. */
 /* .. Array Arguments .. */
 /* DOUBLE PRECISION A( LDA, * ), TAU( * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DGELQF computes an LQ factorization of a real M-by-N matrix A: */
 /* > */
 /* > A = ( L 0 ) * Q */
 /* > */
 /* > where: */
 /* > */
 /* > Q is a N-by-N orthogonal matrix;
 */
 /* > L is a lower-triangular M-by-M matrix;
 */
 /* > 0 is a M-by-(N-M) zero matrix, if M < N. */
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
 /* > A is DOUBLE PRECISION array, dimension (LDA,N) */
 /* > On entry, the M-by-N matrix A. */
 /* > On exit, the elements on and below the diagonal of the array */
 /* > contain the m-by-fla_min(m,n) lower trapezoidal matrix L (L is */
 /* > lower triangular if m <= n);
 the elements above the diagonal, */
 /* > with the array TAU, represent the orthogonal matrix Q as a */
 /* > product of elementary reflectors (see Further Details). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= fla_max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] TAU */
 /* > \verbatim */
 /* > TAU is DOUBLE PRECISION array, dimension (fla_min(M,N)) */
 /* > The scalar factors of the elementary reflectors (see Further */
 /* > Details). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
 /* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The dimension of the array WORK. LWORK >= fla_max(1,M). */
 /* > For optimum performance LWORK >= M*NB, where NB is the */
 /* > optimal blocksize. */
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
 /* > \ingroup doubleGEcomputational */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > The matrix Q is represented as a product of elementary reflectors */
 /* > */
 /* > Q = H(k) . . . H(2) H(1), where k = fla_min(m,n). */
 /* > */
 /* > Each H(i) has the form */
 /* > */
 /* > H(i) = I - tau * v * v**T */
 /* > */
 /* > where tau is a real scalar, and v is a real vector with */
 /* > v(1:i-1) = 0 and v(i) = 1;
 v(i+1:n) is stored on exit in A(i,i+1:n), */
 /* > and tau in TAU(i). */
 /* > \endverbatim */
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int lapack_dgelqf(integer *m, integer *n, doublereal *a, integer * lda, doublereal *tau, doublereal *work, integer *lwork, integer *info) {
#ifdef FLA_ENABLE_AMD_OPT
/* LQ Factorization through QR using below equation:
 *              LQ = (QR)**H
 * */
 if (*m <= FLA_DELQF_TRAN_THRESH && *n <= FLA_DELQF_TRAN_THRESH)
 {
   doublereal *at;

   if (*lwork == -1)
   {
       /* Call QR to get appropriate work buffer size */
       dgeqrf_fla(n, m, a, n, tau, work, lwork, info);
   }
   else
   {
       /* Allocate transpose matrix */
       at = malloc(*n * *m * sizeof(doublereal));

       /* Do transpose and store it in at */
       fla_dtranspose(m, n, a, lda, at, n);

       /* Call QR for the transposed n x m matrix at */
       dgeqrf_fla(n, m, at, n, tau, work, lwork, info);

       /* Transpose at and store back in a */
       fla_dtranspose(n, m, at, n, a, lda);

       /* Free the transpose matrix */
       free(at);
   }
   return 0;
 }
#endif
 /* System generated locals */
 integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
 /* Local variables */
 integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
 extern /* Subroutine */
 int lapack_dgelq2(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *), dlarft_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
 extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
 integer ldwork, lwkopt;
 logical lquery;
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
 a_offset = 1 + a_dim1;
 a -= a_offset;
 --tau;
 --work;
 /* Function Body */
 *info = 0;
 nb = ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1);
 lwkopt = *m * nb;
 work[1] = (doublereal) lwkopt;
 lquery = *lwork == -1;
 if (*m < 0) {
 *info = -1;
 }
 else if (*n < 0) {
 *info = -2;
 }
 else if (*lda < fla_max(1,*m)) {
 *info = -4;
 }
 else if (*lwork < fla_max(1,*m) && ! lquery) {
 *info = -7;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("DGELQF", &i__1);
 return 0;
 }
 else if (lquery) {
 return 0;
 }
 /* Quick return if possible */
 k = fla_min(*m,*n);
 if (k == 0) {
 work[1] = 1.;
 return 0;
 }
 nbmin = 2;
 nx = 0;
 iws = *m;
 if (nb > 1 && nb < k) {
 /* Determine when to cross over from blocked to unblocked code. */
 /* Computing MAX */
 i__1 = 0; i__2 = ilaenv_(&c__3, "DGELQF", " ", m, n, &c_n1, &c_n1); // , expr subst  
 nx = fla_max(i__1,i__2);
 if (nx < k) {
 /* Determine if workspace is large enough for blocked code. */
 ldwork = *m;
 iws = ldwork * nb;
 if (*lwork < iws) {
 /* Not enough workspace to use optimal NB: reduce NB and */
 /* determine the minimum value of NB. */
 nb = *lwork / ldwork;
 /* Computing MAX */
 i__1 = 2; i__2 = ilaenv_(&c__2, "DGELQF", " ", m, n, &c_n1, & c_n1); // , expr subst  
 nbmin = fla_max(i__1,i__2);
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
 ib = fla_min(i__3,nb);
 /* Compute the LQ factorization of the current block */
 /* A(i:i+ib-1,i:n) */
 i__3 = *n - i__ + 1;
 lapack_dgelq2(&ib, &i__3, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[ 1], &iinfo);
 if (i__ + ib <= *m) {
 /* Form the triangular factor of the block reflector */
 /* H = H(i) H(i+1) . . . H(i+ib-1) */
 i__3 = *n - i__ + 1;
 dlarft_("Forward", "Rowwise", &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &ldwork);
 /* Apply H to A(i+ib:m,i:n) from the right */
 i__3 = *m - i__ - ib + 1;
 i__4 = *n - i__ + 1;
 dlarfb_("Right", "No transpose", "Forward", "Rowwise", &i__3, &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], & ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib + 1], &ldwork);
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
 lapack_dgelq2(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1] , &iinfo);
 }
 work[1] = (doublereal) iws;
 return 0;
 /* End of DGELQF */
 }
 /* lapack_dgelqf */
