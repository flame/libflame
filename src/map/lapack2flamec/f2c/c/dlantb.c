/* ../netlib/v3.9.0/dlantb.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 /* > \brief \b DLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele ment of largest absolute value of a triangular band matrix. */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DLANTB + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlantb. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlantb. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlantb. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* DOUBLE PRECISION FUNCTION DLANTB( NORM, UPLO, DIAG, N, K, AB, */
 /* LDAB, WORK ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER DIAG, NORM, UPLO */
 /* INTEGER K, LDAB, N */
 /* .. */
 /* .. Array Arguments .. */
 /* DOUBLE PRECISION AB( LDAB, * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DLANTB returns the value of the one norm, or the Frobenius norm, or */
 /* > the infinity norm, or the element of largest absolute value of an */
 /* > n by n triangular band matrix A, with ( k + 1 ) diagonals. */
 /* > \endverbatim */
 /* > */
 /* > \return DLANTB */
 /* > \verbatim */
 /* > */
 /* > DLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
 /* > Specifies the value to be returned in DLANTB as described */
 /* > above. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] UPLO */
 /* > \verbatim */
 /* > UPLO is CHARACTER*1 */
 /* > Specifies whether the matrix A is upper or lower triangular. */
 /* > = 'U': Upper triangular */
 /* > = 'L': Lower triangular */
 /* > \endverbatim */
 /* > */
 /* > \param[in] DIAG */
 /* > \verbatim */
 /* > DIAG is CHARACTER*1 */
 /* > Specifies whether or not the matrix A is unit triangular. */
 /* > = 'N': Non-unit triangular */
 /* > = 'U': Unit triangular */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The order of the matrix A. N >= 0. When N = 0, DLANTB is */
 /* > set to zero. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] K */
 /* > \verbatim */
 /* > K is INTEGER */
 /* > The number of super-diagonals of the matrix A if UPLO = 'U', */
 /* > or the number of sub-diagonals of the matrix A if UPLO = 'L'. */
 /* > K >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] AB */
 /* > \verbatim */
 /* > AB is DOUBLE PRECISION array, dimension (LDAB,N) */
 /* > The upper or lower triangular band matrix A, stored in the */
 /* > first k+1 rows of AB. The j-th column of A is stored */
 /* > in the j-th column of the array AB as follows: */
 /* > if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
 */
 /* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+k). */
 /* > Note that when DIAG = 'U', the elements of the array AB */
 /* > corresponding to the diagonal elements of the matrix A are */
 /* > not referenced, but are assumed to be one. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDAB */
 /* > \verbatim */
 /* > LDAB is INTEGER */
 /* > The leading dimension of the array AB. LDAB >= K+1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
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
 /* > \ingroup doubleOTHERauxiliary */
 /* ===================================================================== */
 doublereal dlantb_(char *norm, char *uplo, char *diag, integer *n, integer *k, doublereal *ab, integer *ldab, doublereal *work) {
 /* System generated locals */
 integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
 doublereal ret_val, d__1;
 /* Builtin functions */
 double sqrt(doublereal);
 /* Local variables */
 extern /* Subroutine */
 int dcombssq_(doublereal *, doublereal *);
 integer i__, j, l;
 doublereal sum, ssq[2];
 logical udiag;
 extern logical lsame_(char *, char *);
 doublereal value;
 extern logical disnan_(doublereal *);
 extern /* Subroutine */
 int dlassq_(integer *, doublereal *, integer *, doublereal *, doublereal *);
 doublereal colssq[2];
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
 value = 0.;
 }
 else if (lsame_(norm, "M")) {
 /* Find max(abs(A(i,j))). */
 if (lsame_(diag, "U")) {
 value = 1.;
 if (lsame_(uplo, "U")) {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 /* Computing MAX */
 i__2 = *k + 2 - j;
 i__3 = *k;
 for (i__ = max(i__2,1);
 i__ <= i__3;
 ++i__) {
 sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L10: */
 }
 /* L20: */
 }
 }
 else {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 /* Computing MIN */
 i__2 = *n + 1 - j; i__4 = *k + 1; // , expr subst  
 i__3 = min(i__2,i__4);
 for (i__ = 2;
 i__ <= i__3;
 ++i__) {
 sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L30: */
 }
 /* L40: */
 }
 }
 }
 else {
 value = 0.;
 if (lsame_(uplo, "U")) {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 /* Computing MAX */
 i__3 = *k + 2 - j;
 i__2 = *k + 1;
 for (i__ = max(i__3,1);
 i__ <= i__2;
 ++i__) {
 sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L50: */
 }
 /* L60: */
 }
 }
 else {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 /* Computing MIN */
 i__3 = *n + 1 - j; i__4 = *k + 1; // , expr subst  
 i__2 = min(i__3,i__4);
 for (i__ = 1;
 i__ <= i__2;
 ++i__) {
 sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L70: */
 }
 /* L80: */
 }
 }
 }
 }
 else if (lsame_(norm, "O") || *(unsigned char *) norm == '1') {
 /* Find norm1(A). */
 value = 0.;
 udiag = lsame_(diag, "U");
 if (lsame_(uplo, "U")) {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 if (udiag) {
 sum = 1.;
 /* Computing MAX */
 i__2 = *k + 2 - j;
 i__3 = *k;
 for (i__ = max(i__2,1);
 i__ <= i__3;
 ++i__) {
 sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
 /* L90: */
 }
 }
 else {
 sum = 0.;
 /* Computing MAX */
 i__3 = *k + 2 - j;
 i__2 = *k + 1;
 for (i__ = max(i__3,1);
 i__ <= i__2;
 ++i__) {
 sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
 /* L100: */
 }
 }
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L110: */
 }
 }
 else {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 if (udiag) {
 sum = 1.;
 /* Computing MIN */
 i__3 = *n + 1 - j; i__4 = *k + 1; // , expr subst  
 i__2 = min(i__3,i__4);
 for (i__ = 2;
 i__ <= i__2;
 ++i__) {
 sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
 /* L120: */
 }
 }
 else {
 sum = 0.;
 /* Computing MIN */
 i__3 = *n + 1 - j; i__4 = *k + 1; // , expr subst  
 i__2 = min(i__3,i__4);
 for (i__ = 1;
 i__ <= i__2;
 ++i__) {
 sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
 /* L130: */
 }
 }
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L140: */
 }
 }
 }
 else if (lsame_(norm, "I")) {
 /* Find normI(A). */
 value = 0.;
 if (lsame_(uplo, "U")) {
 if (lsame_(diag, "U")) {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 work[i__] = 1.;
 /* L150: */
 }
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 l = *k + 1 - j;
 /* Computing MAX */
 i__2 = 1; i__3 = j - *k; // , expr subst  
 i__4 = j - 1;
 for (i__ = max(i__2,i__3);
 i__ <= i__4;
 ++i__) {
 work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs( d__1));
 /* L160: */
 }
 /* L170: */
 }
 }
 else {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 work[i__] = 0.;
 /* L180: */
 }
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 l = *k + 1 - j;
 /* Computing MAX */
 i__4 = 1; i__2 = j - *k; // , expr subst  
 i__3 = j;
 for (i__ = max(i__4,i__2);
 i__ <= i__3;
 ++i__) {
 work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs( d__1));
 /* L190: */
 }
 /* L200: */
 }
 }
 }
 else {
 if (lsame_(diag, "U")) {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 work[i__] = 1.;
 /* L210: */
 }
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 l = 1 - j;
 /* Computing MIN */
 i__4 = *n; i__2 = j + *k; // , expr subst  
 i__3 = min(i__4,i__2);
 for (i__ = j + 1;
 i__ <= i__3;
 ++i__) {
 work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs( d__1));
 /* L220: */
 }
 /* L230: */
 }
 }
 else {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 work[i__] = 0.;
 /* L240: */
 }
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 l = 1 - j;
 /* Computing MIN */
 i__4 = *n; i__2 = j + *k; // , expr subst  
 i__3 = min(i__4,i__2);
 for (i__ = j;
 i__ <= i__3;
 ++i__) {
 work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs( d__1));
 /* L250: */
 }
 /* L260: */
 }
 }
 }
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 sum = work[i__];
 if (value < sum || disnan_(&sum)) {
 value = sum;
 }
 /* L270: */
 }
 }
 else if (lsame_(norm, "F") || lsame_(norm, "E")) {
 /* Find normF(A). */
 /* SSQ(1) is scale */
 /* SSQ(2) is sum-of-squares */
 /* For better accuracy, sum each column separately. */
 if (lsame_(uplo, "U")) {
 if (lsame_(diag, "U")) {
 ssq[0] = 1.;
 ssq[1] = (doublereal) (*n);
 if (*k > 0) {
 i__1 = *n;
 for (j = 2;
 j <= i__1;
 ++j) {
 colssq[0] = 0.;
 colssq[1] = 1.;
 /* Computing MIN */
 i__4 = j - 1;
 i__3 = min(i__4,*k);
 /* Computing MAX */
 i__2 = *k + 2 - j;
 dlassq_(&i__3, &ab[max(i__2,1) + j * ab_dim1], &c__1, colssq, &colssq[1]);
 dcombssq_(ssq, colssq);
 /* L280: */
 }
 }
 }
 else {
 ssq[0] = 0.;
 ssq[1] = 1.;
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 colssq[0] = 0.;
 colssq[1] = 1.;
 /* Computing MIN */
 i__4 = j; i__2 = *k + 1; // , expr subst  
 i__3 = min(i__4,i__2);
 /* Computing MAX */
 i__5 = *k + 2 - j;
 dlassq_(&i__3, &ab[max(i__5,1) + j * ab_dim1], &c__1, colssq, &colssq[1]);
 dcombssq_(ssq, colssq);
 /* L290: */
 }
 }
 }
 else {
 if (lsame_(diag, "U")) {
 ssq[0] = 1.;
 ssq[1] = (doublereal) (*n);
 if (*k > 0) {
 i__1 = *n - 1;
 for (j = 1;
 j <= i__1;
 ++j) {
 colssq[0] = 0.;
 colssq[1] = 1.;
 /* Computing MIN */
 i__4 = *n - j;
 i__3 = min(i__4,*k);
 dlassq_(&i__3, &ab[j * ab_dim1 + 2], &c__1, colssq, & colssq[1]);
 dcombssq_(ssq, colssq);
 /* L300: */
 }
 }
 }
 else {
 ssq[0] = 0.;
 ssq[1] = 1.;
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 colssq[0] = 0.;
 colssq[1] = 1.;
 /* Computing MIN */
 i__4 = *n - j + 1; i__2 = *k + 1; // , expr subst  
 i__3 = min(i__4,i__2);
 dlassq_(&i__3, &ab[j * ab_dim1 + 1], &c__1, colssq, & colssq[1]);
 dcombssq_(ssq, colssq);
 /* L310: */
 }
 }
 }
 value = ssq[0] * sqrt(ssq[1]);
 }
 ret_val = value;
 return ret_val;
 /* End of DLANTB */
 }
 /* dlantb_ */
 