/* ../netlib/v3.9.0/zlahef_aa.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static doublecomplex c_b1 = {
0.,0.}
;
 static doublecomplex c_b2 = {
1.,0.}
;
 static integer c__1 = 1;
 /* > \brief \b ZLAHEF_AA */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download ZLAHEF_AA + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahef_ aa.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahef_ aa.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahef_ aa.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE ZLAHEF_AA( UPLO, J1, M, NB, A, LDA, IPIV, */
 /* H, LDH, WORK ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* INTEGER J1, M, NB, LDA, LDH */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IPIV( * ) */
 /* COMPLEX*16 A( LDA, * ), H( LDH, * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DLAHEF_AA factorizes a panel of a complex hermitian matrix A using */
 /* > the Aasen's algorithm. The panel consists of a set of NB rows of A */
 /* > when UPLO is U, or a set of NB columns when UPLO is L. */
 /* > */
 /* > In order to factorize the panel, the Aasen's algorithm requires the */
 /* > last row, or column, of the previous panel. The first row, or column, */
 /* > of A is set to be the first row, or column, of an identity matrix, */
 /* > which is used to factorize the first panel. */
 /* > */
 /* > The resulting J-th row of U, or J-th column of L, is stored in the */
 /* > (J-1)-th row, or column, of A (without the unit diagonals), while */
 /* > the diagonal and subdiagonal of A are overwritten by those of T. */
 /* > */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] UPLO */
 /* > \verbatim */
 /* > UPLO is CHARACTER*1 */
 /* > = 'U': Upper triangle of A is stored;
 */
 /* > = 'L': Lower triangle of A is stored. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] J1 */
 /* > \verbatim */
 /* > J1 is INTEGER */
 /* > The location of the first row, or column, of the panel */
 /* > within the submatrix of A, passed to this routine, e.g., */
 /* > when called by ZHETRF_AA, for the first panel, J1 is 1, */
 /* > while for the remaining panels, J1 is 2. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] M */
 /* > \verbatim */
 /* > M is INTEGER */
 /* > The dimension of the submatrix. M >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] NB */
 /* > \verbatim */
 /* > NB is INTEGER */
 /* > The dimension of the panel to be facotorized. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] A */
 /* > \verbatim */
 /* > A is COMPLEX*16 array, dimension (LDA,M) for */
 /* > the first panel, while dimension (LDA,M+1) for the */
 /* > remaining panels. */
 /* > */
 /* > On entry, A contains the last row, or column, of */
 /* > the previous panel, and the trailing submatrix of A */
 /* > to be factorized, except for the first panel, only */
 /* > the panel is passed. */
 /* > */
 /* > On exit, the leading panel is factorized. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,N). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IPIV */
 /* > \verbatim */
 /* > IPIV is INTEGER array, dimension (N) */
 /* > Details of the row and column interchanges, */
 /* > the row and column k were interchanged with the row and */
 /* > column IPIV(k). */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] H */
 /* > \verbatim */
 /* > H is COMPLEX*16 workspace, dimension (LDH,NB). */
 /* > */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDH */
 /* > \verbatim */
 /* > LDH is INTEGER */
 /* > The leading dimension of the workspace H. LDH >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is COMPLEX*16 workspace, dimension (M). */
 /* > \endverbatim */
 /* > */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2017 */
 /* > \ingroup complex16HEcomputational */
 /* ===================================================================== */
 /* Subroutine */
 int zlahef_aa_(char *uplo, integer *j1, integer *m, integer *nb, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex * h__, integer *ldh, doublecomplex *work) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
 char buffer[256];
 snprintf(buffer, 256,"zlahef_aa inputs: uplo %c, j1 %" FLA_IS ", m %" FLA_IS ", nb %" FLA_IS ", lda %" FLA_IS ", ldh %" FLA_IS "",*uplo, *j1, *m, *nb, *lda, *ldh);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, h_dim1, h_offset, i__1, i__2;
 doublereal d__1;
 doublecomplex z__1, z__2;
 /* Builtin functions */
 void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, doublecomplex *, doublecomplex *);
 /* Local variables */
 integer j, k, i1, k1, i2, mj;
 doublecomplex piv, alpha;
 extern logical lsame_(char *, char *);
 extern /* Subroutine */
 int zscal_(integer *, doublecomplex *, doublecomplex *, integer *), zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *), zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), zlacgv_( integer *, doublecomplex *, integer *);
 extern integer izamax_(integer *, doublecomplex *, integer *);
 extern /* Subroutine */
 int zlaset_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, integer *);
 /* -- LAPACK computational routine (version 3.8.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* November 2017 */
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
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 --ipiv;
 h_dim1 = *ldh;
 h_offset = 1 + h_dim1;
 h__ -= h_offset;
 --work;
 /* Function Body */
 j = 1;
 /* K1 is the first column of the panel to be factorized */
 /* i.e., K1 is 2 for the first block column, and 1 for the rest of the blocks */
 k1 = 2 - *j1 + 1;
 if (lsame_(uplo, "U")) {
 /* ..................................................... */
 /* Factorize A as U**T*D*U using the upper triangle of A */
 /* ..................................................... */
 L10: if (j > min(*m,*nb)) {
 goto L20;
 }
 /* K is the column to be factorized */
 /* when being called from ZHETRF_AA, */
 /* > for the first block column, J1 is 1, hence J1+J-1 is J, */
 /* > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */
 k = *j1 + j - 1;
 if (j == *m) {
 /* Only need to compute T(J, J) */
 mj = 1;
 }
 else {
 mj = *m - j + 1;
 }
 /* H(J:N, J) := A(J, J:N) - H(J:N, 1:(J-1)) * L(J1:(J-1), J), */
 /* where H(J:N, J) has been initialized to be A(J, J:N) */
 if (k > 2) {
 /* K is the column to be factorized */
 /* > for the first block column, K is J, skipping the first two */
 /* columns */
 /* > for the rest of the columns, K is J+1, skipping only the */
 /* first column */
 i__1 = j - k1;
 zlacgv_(&i__1, &a[j * a_dim1 + 1], &c__1);
 i__1 = j - k1;
 z__1.r = -1.; z__1.i = -0.; // , expr subst  
 zgemv_("No transpose", &mj, &i__1, &z__1, &h__[j + k1 * h_dim1], ldh, &a[j * a_dim1 + 1], &c__1, &c_b2, &h__[j + j * h_dim1], &c__1);
 i__1 = j - k1;
 zlacgv_(&i__1, &a[j * a_dim1 + 1], &c__1);
 }
 /* Copy H(i:n, i) into WORK */
 zcopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);
 if (j > k1) {
 /* Compute WORK := WORK - L(J-1, J:N) * T(J-1,J), */
 /* where A(J-1, J) stores T(J-1, J) and A(J-2, J:N) stores U(J-1, J:N) */
 d_cnjg(&z__2, &a[k - 1 + j * a_dim1]);
 z__1.r = -z__2.r; z__1.i = -z__2.i; // , expr subst  
 alpha.r = z__1.r; alpha.i = z__1.i; // , expr subst  
 zaxpy_(&mj, &alpha, &a[k - 2 + j * a_dim1], lda, &work[1], &c__1);
 }
 /* Set A(J, J) = T(J, J) */
 i__1 = k + j * a_dim1;
 d__1 = work[1].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 if (j < *m) {
 /* Compute WORK(2:N) = T(J, J) L(J, (J+1):N) */
 /* where A(J, J) stores T(J, J) and A(J-1, (J+1):N) stores U(J, (J+1):N) */
 if (k > 1) {
 i__1 = k + j * a_dim1;
 z__1.r = -a[i__1].r; z__1.i = -a[i__1].i; // , expr subst  
 alpha.r = z__1.r; alpha.i = z__1.i; // , expr subst  
 i__1 = *m - j;
 zaxpy_(&i__1, &alpha, &a[k - 1 + (j + 1) * a_dim1], lda, & work[2], &c__1);
 }
 /* Find max(|WORK(2:n)|) */
 i__1 = *m - j;
 i2 = izamax_(&i__1, &work[2], &c__1) + 1;
 i__1 = i2;
 piv.r = work[i__1].r; piv.i = work[i__1].i; // , expr subst  
 /* Apply hermitian pivot */
 if (i2 != 2 && (piv.r != 0. || piv.i != 0.)) {
 /* Swap WORK(I1) and WORK(I2) */
 i1 = 2;
 i__1 = i2;
 i__2 = i1;
 work[i__1].r = work[i__2].r; work[i__1].i = work[i__2].i; // , expr subst  
 i__1 = i1;
 work[i__1].r = piv.r; work[i__1].i = piv.i; // , expr subst  
 /* Swap A(I1, I1+1:N) with A(I1+1:N, I2) */
 i1 = i1 + j - 1;
 i2 = i2 + j - 1;
 i__1 = i2 - i1 - 1;
 zswap_(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda, &a[* j1 + i1 + i2 * a_dim1], &c__1);
 i__1 = i2 - i1;
 zlacgv_(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda);
 i__1 = i2 - i1 - 1;
 zlacgv_(&i__1, &a[*j1 + i1 + i2 * a_dim1], &c__1);
 /* Swap A(I1, I2+1:N) with A(I2, I2+1:N) */
 if (i2 < *m) {
 i__1 = *m - i2;
 zswap_(&i__1, &a[*j1 + i1 - 1 + (i2 + 1) * a_dim1], lda, & a[*j1 + i2 - 1 + (i2 + 1) * a_dim1], lda);
 }
 /* Swap A(I1, I1) with A(I2,I2) */
 i__1 = i1 + *j1 - 1 + i1 * a_dim1;
 piv.r = a[i__1].r; piv.i = a[i__1].i; // , expr subst  
 i__1 = *j1 + i1 - 1 + i1 * a_dim1;
 i__2 = *j1 + i2 - 1 + i2 * a_dim1;
 a[i__1].r = a[i__2].r; a[i__1].i = a[i__2].i; // , expr subst  
 i__1 = *j1 + i2 - 1 + i2 * a_dim1;
 a[i__1].r = piv.r; a[i__1].i = piv.i; // , expr subst  
 /* Swap H(I1, 1:J1) with H(I2, 1:J1) */
 i__1 = i1 - 1;
 zswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
 ipiv[i1] = i2;
 if (i1 > k1 - 1) {
 /* Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
 /* skipping the first column */
 i__1 = i1 - k1 + 1;
 zswap_(&i__1, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * a_dim1 + 1], &c__1);
 }
 }
 else {
 ipiv[j + 1] = j + 1;
 }
 /* Set A(J, J+1) = T(J, J+1) */
 i__1 = k + (j + 1) * a_dim1;
 a[i__1].r = work[2].r; a[i__1].i = work[2].i; // , expr subst  
 if (j < *nb) {
 /* Copy A(J+1:N, J+1) into H(J:N, J), */
 i__1 = *m - j;
 zcopy_(&i__1, &a[k + 1 + (j + 1) * a_dim1], lda, &h__[j + 1 + (j + 1) * h_dim1], &c__1);
 }
 /* Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1), */
 /* where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1) */
 if (j < *m - 1) {
 i__1 = k + (j + 1) * a_dim1;
 if (a[i__1].r != 0. || a[i__1].i != 0.) {
 z_div(&z__1, &c_b2, &a[k + (j + 1) * a_dim1]);
 alpha.r = z__1.r; alpha.i = z__1.i; // , expr subst  
 i__1 = *m - j - 1;
 zcopy_(&i__1, &work[3], &c__1, &a[k + (j + 2) * a_dim1], lda);
 i__1 = *m - j - 1;
 zscal_(&i__1, &alpha, &a[k + (j + 2) * a_dim1], lda);
 }
 else {
 i__1 = *m - j - 1;
 zlaset_("Full", &c__1, &i__1, &c_b1, &c_b1, &a[k + (j + 2) * a_dim1], lda);
 }
 }
 }
 ++j;
 goto L10;
 L20: ;
 }
 else {
 /* ..................................................... */
 /* Factorize A as L*D*L**T using the lower triangle of A */
 /* ..................................................... */
 L30: if (j > min(*m,*nb)) {
 goto L40;
 }
 /* K is the column to be factorized */
 /* when being called from ZHETRF_AA, */
 /* > for the first block column, J1 is 1, hence J1+J-1 is J, */
 /* > for the rest of the columns, J1 is 2, and J1+J-1 is J+1, */
 k = *j1 + j - 1;
 if (j == *m) {
 /* Only need to compute T(J, J) */
 mj = 1;
 }
 else {
 mj = *m - j + 1;
 }
 /* H(J:N, J) := A(J:N, J) - H(J:N, 1:(J-1)) * L(J, J1:(J-1))^T, */
 /* where H(J:N, J) has been initialized to be A(J:N, J) */
 if (k > 2) {
 /* K is the column to be factorized */
 /* > for the first block column, K is J, skipping the first two */
 /* columns */
 /* > for the rest of the columns, K is J+1, skipping only the */
 /* first column */
 i__1 = j - k1;
 zlacgv_(&i__1, &a[j + a_dim1], lda);
 i__1 = j - k1;
 z__1.r = -1.; z__1.i = -0.; // , expr subst  
 zgemv_("No transpose", &mj, &i__1, &z__1, &h__[j + k1 * h_dim1], ldh, &a[j + a_dim1], lda, &c_b2, &h__[j + j * h_dim1], & c__1);
 i__1 = j - k1;
 zlacgv_(&i__1, &a[j + a_dim1], lda);
 }
 /* Copy H(J:N, J) into WORK */
 zcopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);
 if (j > k1) {
 /* Compute WORK := WORK - L(J:N, J-1) * T(J-1,J), */
 /* where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1) */
 d_cnjg(&z__2, &a[j + (k - 1) * a_dim1]);
 z__1.r = -z__2.r; z__1.i = -z__2.i; // , expr subst  
 alpha.r = z__1.r; alpha.i = z__1.i; // , expr subst  
 zaxpy_(&mj, &alpha, &a[j + (k - 2) * a_dim1], &c__1, &work[1], & c__1);
 }
 /* Set A(J, J) = T(J, J) */
 i__1 = j + k * a_dim1;
 d__1 = work[1].r;
 a[i__1].r = d__1; a[i__1].i = 0.; // , expr subst  
 if (j < *m) {
 /* Compute WORK(2:N) = T(J, J) L((J+1):N, J) */
 /* where A(J, J) = T(J, J) and A((J+1):N, J-1) = L((J+1):N, J) */
 if (k > 1) {
 i__1 = j + k * a_dim1;
 z__1.r = -a[i__1].r; z__1.i = -a[i__1].i; // , expr subst  
 alpha.r = z__1.r; alpha.i = z__1.i; // , expr subst  
 i__1 = *m - j;
 zaxpy_(&i__1, &alpha, &a[j + 1 + (k - 1) * a_dim1], &c__1, & work[2], &c__1);
 }
 /* Find max(|WORK(2:n)|) */
 i__1 = *m - j;
 i2 = izamax_(&i__1, &work[2], &c__1) + 1;
 i__1 = i2;
 piv.r = work[i__1].r; piv.i = work[i__1].i; // , expr subst  
 /* Apply hermitian pivot */
 if (i2 != 2 && (piv.r != 0. || piv.i != 0.)) {
 /* Swap WORK(I1) and WORK(I2) */
 i1 = 2;
 i__1 = i2;
 i__2 = i1;
 work[i__1].r = work[i__2].r; work[i__1].i = work[i__2].i; // , expr subst  
 i__1 = i1;
 work[i__1].r = piv.r; work[i__1].i = piv.i; // , expr subst  
 /* Swap A(I1+1:N, I1) with A(I2, I1+1:N) */
 i1 = i1 + j - 1;
 i2 = i2 + j - 1;
 i__1 = i2 - i1 - 1;
 zswap_(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[ i2 + (*j1 + i1) * a_dim1], lda);
 i__1 = i2 - i1;
 zlacgv_(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1);
 i__1 = i2 - i1 - 1;
 zlacgv_(&i__1, &a[i2 + (*j1 + i1) * a_dim1], lda);
 /* Swap A(I2+1:N, I1) with A(I2+1:N, I2) */
 if (i2 < *m) {
 i__1 = *m - i2;
 zswap_(&i__1, &a[i2 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[i2 + 1 + (*j1 + i2 - 1) * a_dim1], &c__1);
 }
 /* Swap A(I1, I1) with A(I2, I2) */
 i__1 = i1 + (*j1 + i1 - 1) * a_dim1;
 piv.r = a[i__1].r; piv.i = a[i__1].i; // , expr subst  
 i__1 = i1 + (*j1 + i1 - 1) * a_dim1;
 i__2 = i2 + (*j1 + i2 - 1) * a_dim1;
 a[i__1].r = a[i__2].r; a[i__1].i = a[i__2].i; // , expr subst  
 i__1 = i2 + (*j1 + i2 - 1) * a_dim1;
 a[i__1].r = piv.r; a[i__1].i = piv.i; // , expr subst  
 /* Swap H(I1, I1:J1) with H(I2, I2:J1) */
 i__1 = i1 - 1;
 zswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
 ipiv[i1] = i2;
 if (i1 > k1 - 1) {
 /* Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
 /* skipping the first column */
 i__1 = i1 - k1 + 1;
 zswap_(&i__1, &a[i1 + a_dim1], lda, &a[i2 + a_dim1], lda);
 }
 }
 else {
 ipiv[j + 1] = j + 1;
 }
 /* Set A(J+1, J) = T(J+1, J) */
 i__1 = j + 1 + k * a_dim1;
 a[i__1].r = work[2].r; a[i__1].i = work[2].i; // , expr subst  
 if (j < *nb) {
 /* Copy A(J+1:N, J+1) into H(J+1:N, J), */
 i__1 = *m - j;
 zcopy_(&i__1, &a[j + 1 + (k + 1) * a_dim1], &c__1, &h__[j + 1 + (j + 1) * h_dim1], &c__1);
 }
 /* Compute L(J+2, J+1) = WORK( 3:N ) / T(J, J+1), */
 /* where A(J, J+1) = T(J, J+1) and A(J+2:N, J) = L(J+2:N, J+1) */
 if (j < *m - 1) {
 i__1 = j + 1 + k * a_dim1;
 if (a[i__1].r != 0. || a[i__1].i != 0.) {
 z_div(&z__1, &c_b2, &a[j + 1 + k * a_dim1]);
 alpha.r = z__1.r; alpha.i = z__1.i; // , expr subst  
 i__1 = *m - j - 1;
 zcopy_(&i__1, &work[3], &c__1, &a[j + 2 + k * a_dim1], & c__1);
 i__1 = *m - j - 1;
 zscal_(&i__1, &alpha, &a[j + 2 + k * a_dim1], &c__1);
 }
 else {
 i__1 = *m - j - 1;
 zlaset_("Full", &i__1, &c__1, &c_b1, &c_b1, &a[j + 2 + k * a_dim1], lda);
 }
 }
 }
 ++j;
 goto L30;
 L40: ;
 }
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of ZLAHEF_AA */
 }
 /* zlahef_aa__ */
 
