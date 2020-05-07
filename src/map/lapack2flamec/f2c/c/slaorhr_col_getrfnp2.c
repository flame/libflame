/* ../netlib/v3.9.0/slaorhr_col_getrfnp2.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static real c_b3 = 1.f;
 static integer c__1 = 1;
 static real c_b19 = -1.f;
 /* > \brief \b SLAORHR_COL_GETRFNP2 */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DLAORHR_GETRF2NP + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaorhr _col_getrfnp2.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaorhr _col_getrfnp2.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaorhr _col_getrfnp2.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE SLAORHR_COL_GETRFNP2( M, N, A, LDA, D, INFO ) */
 /* .. Scalar Arguments .. */
 /* INTEGER INFO, LDA, M, N */
 /* .. */
 /* .. Array Arguments .. */
 /* REAL A( LDA, * ), D( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > SLAORHR_COL_GETRFNP2 computes the modified LU factorization without */
 /* > pivoting of a real general M-by-N matrix A. The factorization has */
 /* > the form: */
 /* > */
 /* > A - S = L * U, */
 /* > */
 /* > where: */
 /* > S is a m-by-n diagonal sign matrix with the diagonal D, so that */
 /* > D(i) = S(i,i), 1 <= i <= min(M,N). The diagonal D is constructed */
 /* > as D(i)=-SIGN(A(i,i)), where A(i,i) is the value after performing */
 /* > i-1 steps of Gaussian elimination. This means that the diagonal */
 /* > element at each step of "modified" Gaussian elimination is at */
 /* > least one in absolute value (so that division-by-zero not */
 /* > possible during the division by the diagonal element);
 */
 /* > */
 /* > L is a M-by-N lower triangular matrix with unit diagonal elements */
 /* > (lower trapezoidal if M > N);
 */
 /* > */
 /* > and U is a M-by-N upper triangular matrix */
 /* > (upper trapezoidal if M < N). */
 /* > */
 /* > This routine is an auxiliary routine used in the Householder */
 /* > reconstruction routine SORHR_COL. In SORHR_COL, this routine is */
 /* > applied to an M-by-N matrix A with orthonormal columns, where each */
 /* > element is bounded by one in absolute value. With the choice of */
 /* > the matrix S above, one can show that the diagonal element at each */
 /* > step of Gaussian elimination is the largest (in absolute value) in */
 /* > the column on or below the diagonal, so that no pivoting is required */
 /* > for numerical stability [1]. */
 /* > */
 /* > For more details on the Householder reconstruction algorithm, */
 /* > including the modified LU factorization, see [1]. */
 /* > */
 /* > This is the recursive version of the LU factorization algorithm. */
 /* > Denote A - S by B. The algorithm divides the matrix B into four */
 /* > submatrices: */
 /* > */
 /* > [ B11 | B12 ] where B11 is n1 by n1, */
 /* > B = [ -----|----- ] B21 is (m-n1) by n1, */
 /* > [ B21 | B22 ] B12 is n1 by n2, */
 /* > B22 is (m-n1) by n2, */
 /* > with n1 = min(m,n)/2, n2 = n-n1. */
 /* > */
 /* > */
 /* > The subroutine calls itself to factor B11, solves for B21, */
 /* > solves for B12, updates B22, then calls itself to factor B22. */
 /* > */
 /* > For more details on the recursive LU algorithm, see [2]. */
 /* > */
 /* > SLAORHR_COL_GETRFNP2 is called to factorize a block by the blocked */
 /* > routine SLAORHR_COL_GETRFNP, which uses blocked code calling */
 /* . Level 3 BLAS to update the submatrix. However, SLAORHR_COL_GETRFNP2 */
 /* > is self-sufficient and can be used without SLAORHR_COL_GETRFNP. */
 /* > */
 /* > [1] "Reconstructing Householder vectors from tall-skinny QR", */
 /* > G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen, */
 /* > E. Solomonik, J. Parallel Distrib. Comput., */
 /* > vol. 85, pp. 3-31, 2015. */
 /* > */
 /* > [2] "Recursion leads to automatic variable blocking for dense linear */
 /* > algebra algorithms", F. Gustavson, IBM J. of Res. and Dev., */
 /* > vol. 41, no. 6, pp. 737-755, 1997. */
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
 /* > On entry, the M-by-N matrix to be factored. */
 /* > On exit, the factors L and U from the factorization */
 /* > A-S=L*U;
 the unit diagonal elements of L are not stored. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] D */
 /* > \verbatim */
 /* > D is REAL array, dimension min(M,N) */
 /* > The diagonal elements of the diagonal M-by-N sign matrix S, */
 /* > D(i) = S(i,i), where 1 <= i <= min(M,N). The elements can */
 /* > be only plus or minus one. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value */
 /* > \endverbatim */
 /* > */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2019 */
 /* > \ingroup realGEcomputational */
 /* > \par Contributors: */
 /* ================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > November 2019, Igor Kozachenko, */
 /* > Computer Science Division, */
 /* > University of California, Berkeley */
 /* > */
 /* > \endverbatim */
 /* ===================================================================== */
 /* Subroutine */
 int slaorhr_col_getrfnp2_(integer *m, integer *n, real *a, integer *lda, real *d__, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, i__1;
 real r__1;
 /* Builtin functions */
 double r_sign(real *, real *);
 /* Local variables */
 integer i__, n1, n2, iinfo;
 extern /* Subroutine */
 int sscal_(integer *, real *, real *, integer *), sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
 real sfmin;
 extern /* Subroutine */
 int strsm_(char *, char *, char *, char *, integer *, integer *, real *, real *, integer *, real *, integer * );
 extern real slamch_(char *);
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 /* -- LAPACK computational routine (version 3.9.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* November 2019 */
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
 /* Test the input parameters */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 --d__;
 /* Function Body */
 *info = 0;
 if (*m < 0) {
 *info = -1;
 }
 else if (*n < 0) {
 *info = -2;
 }
 else if (*lda < max(1,*m)) {
 *info = -4;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("SLAORHR_COL_GETRFNP2", &i__1);
 return 0;
 }
 /* Quick return if possible */
 if (min(*m,*n) == 0) {
 return 0;
 }
 if (*m == 1) {
 /* One row case, (also recursion termination case), */
 /* use unblocked code */
 /* Transfer the sign */
 d__[1] = -r_sign(&c_b3, &a[a_dim1 + 1]);
 /* Construct the row of U */
 a[a_dim1 + 1] -= d__[1];
 }
 else if (*n == 1) {
 /* One column case, (also recursion termination case), */
 /* use unblocked code */
 /* Transfer the sign */
 d__[1] = -r_sign(&c_b3, &a[a_dim1 + 1]);
 /* Construct the row of U */
 a[a_dim1 + 1] -= d__[1];
 /* Scale the elements 2:M of the column */
 /* Determine machine safe minimum */
 sfmin = slamch_("S");
 /* Construct the subdiagonal elements of L */
 if ((r__1 = a[a_dim1 + 1], abs(r__1)) >= sfmin) {
 i__1 = *m - 1;
 r__1 = 1.f / a[a_dim1 + 1];
 sscal_(&i__1, &r__1, &a[a_dim1 + 2], &c__1);
 }
 else {
 i__1 = *m;
 for (i__ = 2;
 i__ <= i__1;
 ++i__) {
 a[i__ + a_dim1] /= a[a_dim1 + 1];
 }
 }
 }
 else {
 /* Divide the matrix B into four submatrices */
 n1 = min(*m,*n) / 2;
 n2 = *n - n1;
 /* Factor B11, recursive call */
 slaorhr_col_getrfnp2_(&n1, &n1, &a[a_offset], lda, &d__[1], &iinfo);
 /* Solve for B21 */
 i__1 = *m - n1;
 strsm_("R", "U", "N", "N", &i__1, &n1, &c_b3, &a[a_offset], lda, &a[ n1 + 1 + a_dim1], lda);
 /* Solve for B12 */
 strsm_("L", "L", "N", "U", &n1, &n2, &c_b3, &a[a_offset], lda, &a[(n1 + 1) * a_dim1 + 1], lda);
 /* Update B22, i.e. compute the Schur complement */
 /* B22 := B22 - B21*B12 */
 i__1 = *m - n1;
 sgemm_("N", "N", &i__1, &n2, &n1, &c_b19, &a[n1 + 1 + a_dim1], lda, & a[(n1 + 1) * a_dim1 + 1], lda, &c_b3, &a[n1 + 1 + (n1 + 1) * a_dim1], lda);
 /* Factor B22, recursive call */
 i__1 = *m - n1;
 slaorhr_col_getrfnp2_(&i__1, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &d__[n1 + 1], &iinfo);
 }
 return 0;
 /* End of SLAORHR_COL_GETRFNP2 */
 }
 /* slaorhr_col_getrfnp2__ */
 