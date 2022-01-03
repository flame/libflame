/* ../netlib/v3.9.0/dlasyf_aa.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static doublereal c_b6 = -1.;
 static integer c__1 = 1;
 static doublereal c_b8 = 1.;
 static doublereal c_b22 = 0.;
 /* > \brief \b DLASYF_AA */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DLASYF_AA + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasyf_ aa.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasyf_ aa.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasyf_ aa.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE DLASYF_AA( UPLO, J1, M, NB, A, LDA, IPIV, */
 /* H, LDH, WORK ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* INTEGER J1, M, NB, LDA, LDH */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IPIV( * ) */
 /* DOUBLE PRECISION A( LDA, * ), H( LDH, * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DLATRF_AA factorizes a panel of a real symmetric matrix A using */
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
 /* > when called by DSYTRF_AA, for the first panel, J1 is 1, */
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
 /* > A is DOUBLE PRECISION array, dimension (LDA,M) for */
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
 /* > The leading dimension of the array A. LDA >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IPIV */
 /* > \verbatim */
 /* > IPIV is INTEGER array, dimension (M) */
 /* > Details of the row and column interchanges, */
 /* > the row and column k were interchanged with the row and */
 /* > column IPIV(k). */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] H */
 /* > \verbatim */
 /* > H is DOUBLE PRECISION workspace, dimension (LDH,NB). */
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
 /* > WORK is DOUBLE PRECISION workspace, dimension (M). */
 /* > \endverbatim */
 /* > */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2017 */
 /* > \ingroup doubleSYcomputational */
 /* ===================================================================== */
 /* Subroutine */
 int dlasyf_aa_(char *uplo, integer *j1, integer *m, integer *nb, doublereal *a, integer *lda, integer *ipiv, doublereal *h__, integer *ldh, doublereal *work) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
 char buffer[256]; 
 snprintf(buffer, 256,"dlasyf_aa inputs: uplo %c, j1 %" FLA_IS ", m %" FLA_IS ", nb %" FLA_IS ", lda %" FLA_IS ", ldh %" FLA_IS "",*uplo, *j1, *m, *nb, *lda, *ldh);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, h_dim1, h_offset, i__1;
 /* Local variables */
 integer j, k, i1, k1, i2, mj;
 doublereal piv, alpha;
 extern /* Subroutine */
 int dscal_(integer *, doublereal *, doublereal *, integer *);
 extern logical lsame_(char *, char *);
 extern /* Subroutine */
 int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), dcopy_(integer *, doublereal *, integer *, doublereal *, integer *), dswap_(integer *, doublereal *, integer *, doublereal *, integer *), daxpy_( integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
 extern integer idamax_(integer *, doublereal *, integer *);
 extern /* Subroutine */
 int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *);
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
 /* when being called from DSYTRF_AA, */
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
 /* H(J:M, J) := A(J, J:M) - H(J:M, 1:(J-1)) * L(J1:(J-1), J), */
 /* where H(J:M, J) has been initialized to be A(J, J:M) */
 if (k > 2) {
 /* K is the column to be factorized */
 /* > for the first block column, K is J, skipping the first two */
 /* columns */
 /* > for the rest of the columns, K is J+1, skipping only the */
 /* first column */
 i__1 = j - k1;
 dgemv_("No transpose", &mj, &i__1, &c_b6, &h__[j + k1 * h_dim1], ldh, &a[j * a_dim1 + 1], &c__1, &c_b8, &h__[j + j * h_dim1], &c__1);
 }
 /* Copy H(i:M, i) into WORK */
 dcopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);
 if (j > k1) {
 /* Compute WORK := WORK - L(J-1, J:M) * T(J-1,J), */
 /* where A(J-1, J) stores T(J-1, J) and A(J-2, J:M) stores U(J-1, J:M) */
 alpha = -a[k - 1 + j * a_dim1];
 daxpy_(&mj, &alpha, &a[k - 2 + j * a_dim1], lda, &work[1], &c__1);
 }
 /* Set A(J, J) = T(J, J) */
 a[k + j * a_dim1] = work[1];
 if (j < *m) {
 /* Compute WORK(2:M) = T(J, J) L(J, (J+1):M) */
 /* where A(J, J) stores T(J, J) and A(J-1, (J+1):M) stores U(J, (J+1):M) */
 if (k > 1) {
 alpha = -a[k + j * a_dim1];
 i__1 = *m - j;
 daxpy_(&i__1, &alpha, &a[k - 1 + (j + 1) * a_dim1], lda, & work[2], &c__1);
 }
 /* Find max(|WORK(2:M)|) */
 i__1 = *m - j;
 i2 = idamax_(&i__1, &work[2], &c__1) + 1;
 piv = work[i2];
 /* Apply symmetric pivot */
 if (i2 != 2 && piv != 0.) {
 /* Swap WORK(I1) and WORK(I2) */
 i1 = 2;
 work[i2] = work[i1];
 work[i1] = piv;
 /* Swap A(I1, I1+1:M) with A(I1+1:M, I2) */
 i1 = i1 + j - 1;
 i2 = i2 + j - 1;
 i__1 = i2 - i1 - 1;
 dswap_(&i__1, &a[*j1 + i1 - 1 + (i1 + 1) * a_dim1], lda, &a[* j1 + i1 + i2 * a_dim1], &c__1);
 /* Swap A(I1, I2+1:M) with A(I2, I2+1:M) */
 if (i2 < *m) {
 i__1 = *m - i2;
 dswap_(&i__1, &a[*j1 + i1 - 1 + (i2 + 1) * a_dim1], lda, & a[*j1 + i2 - 1 + (i2 + 1) * a_dim1], lda);
 }
 /* Swap A(I1, I1) with A(I2,I2) */
 piv = a[i1 + *j1 - 1 + i1 * a_dim1];
 a[*j1 + i1 - 1 + i1 * a_dim1] = a[*j1 + i2 - 1 + i2 * a_dim1];
 a[*j1 + i2 - 1 + i2 * a_dim1] = piv;
 /* Swap H(I1, 1:J1) with H(I2, 1:J1) */
 i__1 = i1 - 1;
 dswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
 ipiv[i1] = i2;
 if (i1 > k1 - 1) {
 /* Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
 /* skipping the first column */
 i__1 = i1 - k1 + 1;
 dswap_(&i__1, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * a_dim1 + 1], &c__1);
 }
 }
 else {
 ipiv[j + 1] = j + 1;
 }
 /* Set A(J, J+1) = T(J, J+1) */
 a[k + (j + 1) * a_dim1] = work[2];
 if (j < *nb) {
 /* Copy A(J+1:M, J+1) into H(J:M, J), */
 i__1 = *m - j;
 dcopy_(&i__1, &a[k + 1 + (j + 1) * a_dim1], lda, &h__[j + 1 + (j + 1) * h_dim1], &c__1);
 }
 /* Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1), */
 /* where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1) */
 if (j < *m - 1) {
 if (a[k + (j + 1) * a_dim1] != 0.) {
 alpha = 1. / a[k + (j + 1) * a_dim1];
 i__1 = *m - j - 1;
 dcopy_(&i__1, &work[3], &c__1, &a[k + (j + 2) * a_dim1], lda);
 i__1 = *m - j - 1;
 dscal_(&i__1, &alpha, &a[k + (j + 2) * a_dim1], lda);
 }
 else {
 i__1 = *m - j - 1;
 dlaset_("Full", &c__1, &i__1, &c_b22, &c_b22, &a[k + (j + 2) * a_dim1], lda);
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
 /* when being called from DSYTRF_AA, */
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
 /* H(J:M, J) := A(J:M, J) - H(J:M, 1:(J-1)) * L(J, J1:(J-1))^T, */
 /* where H(J:M, J) has been initialized to be A(J:M, J) */
 if (k > 2) {
 /* K is the column to be factorized */
 /* > for the first block column, K is J, skipping the first two */
 /* columns */
 /* > for the rest of the columns, K is J+1, skipping only the */
 /* first column */
 i__1 = j - k1;
 dgemv_("No transpose", &mj, &i__1, &c_b6, &h__[j + k1 * h_dim1], ldh, &a[j + a_dim1], lda, &c_b8, &h__[j + j * h_dim1], & c__1);
 }
 /* Copy H(J:M, J) into WORK */
 dcopy_(&mj, &h__[j + j * h_dim1], &c__1, &work[1], &c__1);
 if (j > k1) {
 /* Compute WORK := WORK - L(J:M, J-1) * T(J-1,J), */
 /* where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1) */
 alpha = -a[j + (k - 1) * a_dim1];
 daxpy_(&mj, &alpha, &a[j + (k - 2) * a_dim1], &c__1, &work[1], & c__1);
 }
 /* Set A(J, J) = T(J, J) */
 a[j + k * a_dim1] = work[1];
 if (j < *m) {
 /* Compute WORK(2:M) = T(J, J) L((J+1):M, J) */
 /* where A(J, J) = T(J, J) and A((J+1):M, J-1) = L((J+1):M, J) */
 if (k > 1) {
 alpha = -a[j + k * a_dim1];
 i__1 = *m - j;
 daxpy_(&i__1, &alpha, &a[j + 1 + (k - 1) * a_dim1], &c__1, & work[2], &c__1);
 }
 /* Find max(|WORK(2:M)|) */
 i__1 = *m - j;
 i2 = idamax_(&i__1, &work[2], &c__1) + 1;
 piv = work[i2];
 /* Apply symmetric pivot */
 if (i2 != 2 && piv != 0.) {
 /* Swap WORK(I1) and WORK(I2) */
 i1 = 2;
 work[i2] = work[i1];
 work[i1] = piv;
 /* Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
 i1 = i1 + j - 1;
 i2 = i2 + j - 1;
 i__1 = i2 - i1 - 1;
 dswap_(&i__1, &a[i1 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[ i2 + (*j1 + i1) * a_dim1], lda);
 /* Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
 if (i2 < *m) {
 i__1 = *m - i2;
 dswap_(&i__1, &a[i2 + 1 + (*j1 + i1 - 1) * a_dim1], &c__1, &a[i2 + 1 + (*j1 + i2 - 1) * a_dim1], &c__1);
 }
 /* Swap A(I1, I1) with A(I2, I2) */
 piv = a[i1 + (*j1 + i1 - 1) * a_dim1];
 a[i1 + (*j1 + i1 - 1) * a_dim1] = a[i2 + (*j1 + i2 - 1) * a_dim1];
 a[i2 + (*j1 + i2 - 1) * a_dim1] = piv;
 /* Swap H(I1, I1:J1) with H(I2, I2:J1) */
 i__1 = i1 - 1;
 dswap_(&i__1, &h__[i1 + h_dim1], ldh, &h__[i2 + h_dim1], ldh);
 ipiv[i1] = i2;
 if (i1 > k1 - 1) {
 /* Swap L(1:I1-1, I1) with L(1:I1-1, I2), */
 /* skipping the first column */
 i__1 = i1 - k1 + 1;
 dswap_(&i__1, &a[i1 + a_dim1], lda, &a[i2 + a_dim1], lda);
 }
 }
 else {
 ipiv[j + 1] = j + 1;
 }
 /* Set A(J+1, J) = T(J+1, J) */
 a[j + 1 + k * a_dim1] = work[2];
 if (j < *nb) {
 /* Copy A(J+1:M, J+1) into H(J+1:M, J), */
 i__1 = *m - j;
 dcopy_(&i__1, &a[j + 1 + (k + 1) * a_dim1], &c__1, &h__[j + 1 + (j + 1) * h_dim1], &c__1);
 }
 /* Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1), */
 /* where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1) */
 if (j < *m - 1) {
 if (a[j + 1 + k * a_dim1] != 0.) {
 alpha = 1. / a[j + 1 + k * a_dim1];
 i__1 = *m - j - 1;
 dcopy_(&i__1, &work[3], &c__1, &a[j + 2 + k * a_dim1], & c__1);
 i__1 = *m - j - 1;
 dscal_(&i__1, &alpha, &a[j + 2 + k * a_dim1], &c__1);
 }
 else {
 i__1 = *m - j - 1;
 dlaset_("Full", &i__1, &c__1, &c_b22, &c_b22, &a[j + 2 + k * a_dim1], lda);
 }
 }
 }
 ++j;
 goto L30;
 L40: ;
 }
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of DLASYF_AA */
 }
 /* dlasyf_aa__ */
 