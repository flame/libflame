/* ../netlib/v3.9.0/slangb.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 /* > \brief \b SLANGB returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of general band matrix. */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download SLANGB + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slangb. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slangb. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slangb. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* REAL FUNCTION SLANGB( NORM, N, KL, KU, AB, LDAB, */
 /* WORK ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER NORM */
 /* INTEGER KL, KU, LDAB, N */
 /* .. */
 /* .. Array Arguments .. */
 /* REAL AB( LDAB, * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > SLANGB returns the value of the one norm, or the Frobenius norm, or */
 /* > the infinity norm, or the element of largest absolute value of an */
 /* > n by n band matrix A, with kl sub-diagonals and ku super-diagonals. */
 /* > \endverbatim */
 /* > */
 /* > \return SLANGB */
 /* > \verbatim */
 /* > */
 /* > SLANGB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
 /* > ( */
 /* > ( norm1(A), NORM = '1', 'O' or 'o' */
 /* > ( */
 /* > ( normI(A), NORM = 'I' or 'i' */
 /* > ( */
 /* > ( normF(A), NORM = 'F', 'f', 'E' or 'e' */
 /* > */
 /* > where norm1 denotes the one norm of a matrix (maximum column sum), */
 /* > normI denotes the infinity norm of a matrix (maximum row sum) and */
 /* > normF denotes the Frobenius norm of a matrix (square root of sum of */
 /* > squares). Note that max(abs(A(i,j))) is not a consistent matrix norm. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] NORM */
 /* > \verbatim */
 /* > NORM is CHARACTER*1 */
 /* > Specifies the value to be returned in SLANGB as described */
 /* > above. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The order of the matrix A. N >= 0. When N = 0, SLANGB is */
 /* > set to zero. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] KL */
 /* > \verbatim */
 /* > KL is INTEGER */
 /* > The number of sub-diagonals of the matrix A. KL >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] KU */
 /* > \verbatim */
 /* > KU is INTEGER */
 /* > The number of super-diagonals of the matrix A. KU >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] AB */
 /* > \verbatim */
 /* > AB is REAL array, dimension (LDAB,N) */
 /* > The band matrix A, stored in rows 1 to KL+KU+1. The j-th */
 /* > column of A is stored in the j-th column of the array AB as */
 /* > follows: */
 /* > AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDAB */
 /* > \verbatim */
 /* > LDAB is INTEGER */
 /* > The leading dimension of the array AB. LDAB >= KL+KU+1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is REAL array, dimension (MAX(1,LWORK)), */
 /* > where LWORK >= N when NORM = 'I';
 otherwise, WORK is not */
 /* > referenced. */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date December 2016 */
 /* > \ingroup realGBauxiliary */
 /* ===================================================================== */
 real slangb_(char *norm, integer *n, integer *kl, integer *ku, real *ab, integer *ldab, real *work) {
 /* System generated locals */
 integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
 real ret_val, r__1;
 /* Builtin functions */
 double sqrt(doublereal);
 /* Local variables */
 extern /* Subroutine */
 int scombssq_(real *, real *);
 integer i__, j, k, l;
 real sum, ssq[2], temp;
 extern logical lsame_(char *, char *);
 real value;
 extern logical sisnan_(real *);
 real colssq[2];
 extern /* Subroutine */
 int slassq_(integer *, real *, integer *, real *, real *);
 /* -- LAPACK auxiliary routine (version 3.7.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* December 2016 */
 /* .. Scalar Arguments .. */
 /* .. */
 /* .. Array Arguments .. */
 /* .. */
 /* ===================================================================== */
 /* .. Parameters .. */
 /* .. */
 /* .. Local Scalars .. */
 /* .. */
 /* .. Local Arrays .. */
 /* .. */
 /* .. External Functions .. */
 /* .. */
 /* .. External Subroutines .. */
 /* .. */
 /* .. Intrinsic Functions .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Parameter adjustments */
 ab_dim1 = *ldab;
 ab_offset = 1 + ab_dim1;
 ab -= ab_offset;
 --work;
 /* Function Body */
 if (*n == 0) {
 value = 0.f;
 }
 else if (lsame_(norm, "M")) {
 /* Find max(abs(A(i,j))). */
 value = 0.f;
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 /* Computing MAX */
 i__2 = *ku + 2 - j;
 /* Computing MIN */
 i__4 = *n + *ku + 1 - j; i__5 = *kl + *ku + 1; // , expr subst  
 i__3 = min(i__4,i__5);
 for (i__ = max(i__2,1);
 i__ <= i__3;
 ++i__) {
 temp = (r__1 = ab[i__ + j * ab_dim1], abs(r__1));
 if (value < temp || sisnan_(&temp)) {
 value = temp;
 }
 /* L10: */
 }
 /* L20: */
 }
 }
 else if (lsame_(norm, "O") || *(unsigned char *) norm == '1') {
 /* Find norm1(A). */
 value = 0.f;
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 sum = 0.f;
 /* Computing MAX */
 i__3 = *ku + 2 - j;
 /* Computing MIN */
 i__4 = *n + *ku + 1 - j; i__5 = *kl + *ku + 1; // , expr subst  
 i__2 = min(i__4,i__5);
 for (i__ = max(i__3,1);
 i__ <= i__2;
 ++i__) {
 sum += (r__1 = ab[i__ + j * ab_dim1], abs(r__1));
 /* L30: */
 }
 if (value < sum || sisnan_(&sum)) {
 value = sum;
 }
 /* L40: */
 }
 }
 else if (lsame_(norm, "I")) {
 /* Find normI(A). */
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 work[i__] = 0.f;
 /* L50: */
 }
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 k = *ku + 1 - j;
 /* Computing MAX */
 i__2 = 1; i__3 = j - *ku; // , expr subst  
 /* Computing MIN */
 i__5 = *n; i__6 = j + *kl; // , expr subst  
 i__4 = min(i__5,i__6);
 for (i__ = max(i__2,i__3);
 i__ <= i__4;
 ++i__) {
 work[i__] += (r__1 = ab[k + i__ + j * ab_dim1], abs(r__1));
 /* L60: */
 }
 /* L70: */
 }
 value = 0.f;
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 temp = work[i__];
 if (value < temp || sisnan_(&temp)) {
 value = temp;
 }
 /* L80: */
 }
 }
 else if (lsame_(norm, "F") || lsame_(norm, "E")) {
 /* Find normF(A). */
 /* SSQ(1) is scale */
 /* SSQ(2) is sum-of-squares */
 /* For better accuracy, sum each column separately. */
 ssq[0] = 0.f;
 ssq[1] = 1.f;
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 /* Computing MAX */
 i__4 = 1; i__2 = j - *ku; // , expr subst  
 l = max(i__4,i__2);
 k = *ku + 1 - j + l;
 colssq[0] = 0.f;
 colssq[1] = 1.f;
 /* Computing MIN */
 i__2 = *n; i__3 = j + *kl; // , expr subst  
 i__4 = min(i__2,i__3) - l + 1;
 slassq_(&i__4, &ab[k + j * ab_dim1], &c__1, colssq, &colssq[1]);
 scombssq_(ssq, colssq);
 /* L90: */
 }
 value = ssq[0] * sqrt(ssq[1]);
 }
 ret_val = value;
 return ret_val;
 /* End of SLANGB */
 }
 /* slangb_ */
 