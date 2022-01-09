/* ../netlib/v3.9.0/dlaorhr_col_getrfnp2.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static doublereal c_b3 = 1.;
 static integer c__1 = 1;
 static doublereal c_b19 = -1.;
 /* > \brief \b DLAORHR_COL_GETRFNP2 */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DLAORHR_GETRF2NP + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaorhr _col_getrfnp2.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaorhr _col_getrfnp2.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaorhr _col_getrfnp2.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE DLAORHR_COL_GETRFNP2( M, N, A, LDA, D, INFO ) */
 /* .. Scalar Arguments .. */
 /* INTEGER INFO, LDA, M, N */
 /* .. */
 /* .. Array Arguments .. */
 /* DOUBLE PRECISION A( LDA, * ), D( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DLAORHR_COL_GETRFNP2 computes the modified LU factorization without */
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
 /* > reconstruction routine DORHR_COL. In DORHR_COL, this routine is */
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
 /* > DLAORHR_COL_GETRFNP2 is called to factorize a block by the blocked */
 /* > routine DLAORHR_COL_GETRFNP, which uses blocked code calling */
 /* . Level 3 BLAS to update the submatrix. However, DLAORHR_COL_GETRFNP2 */
 /* > is self-sufficient and can be used without DLAORHR_COL_GETRFNP. */
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
 /* > A is DOUBLE PRECISION array, dimension (LDA,N) */
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
 /* > D is DOUBLE PRECISION array, dimension min(M,N) */
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
 /* > \ingroup doubleGEcomputational */
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
 int dlaorhr_col_getrfnp2_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *d__, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, i__1;
 doublereal d__1;
 /* Builtin functions */
 double d_sign(doublereal *, doublereal *);
 /* Local variables */
 integer i__, n1, n2;
 extern /* Subroutine */
 int dscal_(integer *, doublereal *, doublereal *, integer *), dgemm_(char *, char *, integer *, integer *, integer * , doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
 integer iinfo;
 doublereal sfmin;
 extern /* Subroutine */
 int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
 extern doublereal dlamch_(char *);
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
 xerbla_("DLAORHR_COL_GETRFNP2", &i__1);
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
 d__[1] = -d_sign(&c_b3, &a[a_dim1 + 1]);
 /* Construct the row of U */
 a[a_dim1 + 1] -= d__[1];
 }
 else if (*n == 1) {
 /* One column case, (also recursion termination case), */
 /* use unblocked code */
 /* Transfer the sign */
 d__[1] = -d_sign(&c_b3, &a[a_dim1 + 1]);
 /* Construct the row of U */
 a[a_dim1 + 1] -= d__[1];
 /* Scale the elements 2:M of the column */
 /* Determine machine safe minimum */
 sfmin = dlamch_("S");
 /* Construct the subdiagonal elements of L */
 if ((d__1 = a[a_dim1 + 1], f2c_dabs(d__1)) >= sfmin) {
 i__1 = *m - 1;
 d__1 = 1. / a[a_dim1 + 1];
 dscal_(&i__1, &d__1, &a[a_dim1 + 2], &c__1);
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
 dlaorhr_col_getrfnp2_(&n1, &n1, &a[a_offset], lda, &d__[1], &iinfo);
 /* Solve for B21 */
 i__1 = *m - n1;
 dtrsm_("R", "U", "N", "N", &i__1, &n1, &c_b3, &a[a_offset], lda, &a[ n1 + 1 + a_dim1], lda);
 /* Solve for B12 */
 dtrsm_("L", "L", "N", "U", &n1, &n2, &c_b3, &a[a_offset], lda, &a[(n1 + 1) * a_dim1 + 1], lda);
 /* Update B22, i.e. compute the Schur complement */
 /* B22 := B22 - B21*B12 */
 i__1 = *m - n1;
 dgemm_("N", "N", &i__1, &n2, &n1, &c_b19, &a[n1 + 1 + a_dim1], lda, & a[(n1 + 1) * a_dim1 + 1], lda, &c_b3, &a[n1 + 1 + (n1 + 1) * a_dim1], lda);
 /* Factor B22, recursive call */
 i__1 = *m - n1;
 dlaorhr_col_getrfnp2_(&i__1, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &d__[n1 + 1], &iinfo);
 }
 return 0;
 /* End of DLAORHR_COL_GETRFNP2 */
 }
 /* dlaorhr_col_getrfnp2__ */
 