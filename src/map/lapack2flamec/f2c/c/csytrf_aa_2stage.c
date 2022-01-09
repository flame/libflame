/* ../netlib/v3.9.0/csytrf_aa_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static complex c_b1 = {
0.f,0.f}
;
 static complex c_b2 = {
1.f,0.f}
;
 static integer c__1 = 1;
 static integer c_n1 = -1;
 /* > \brief \b CSYTRF_AA_2STAGE */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download CSYTRF_AA_2STAGE + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrf_ aa_2stage.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrf_ aa_2stage.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrf_ aa_2stage.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE CSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, */
 /* IPIV2, WORK, LWORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* INTEGER N, LDA, LTB, LWORK, INFO */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IPIV( * ), IPIV2( * ) */
 /* COMPLEX A( LDA, * ), TB( * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > CSYTRF_AA_2STAGE computes the factorization of a complex symmetric matrix A */
 /* > using the Aasen's algorithm. The form of the factorization is */
 /* > */
 /* > A = U**T*T*U or A = L*T*L**T */
 /* > */
 /* > where U (or L) is a product of permutation and unit upper (lower) */
 /* > triangular matrices, and T is a complex symmetric band matrix with the */
 /* > bandwidth of NB (NB is internally selected and stored in TB( 1 ), and T is */
 /* > LU factorized with partial pivoting). */
 /* > */
 /* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
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
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The order of the matrix A. N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] A */
 /* > \verbatim */
 /* > A is COMPLEX array, dimension (LDA,N) */
 /* > On entry, the hermitian matrix A. If UPLO = 'U', the leading */
 /* > N-by-N upper triangular part of A contains the upper */
 /* > triangular part of the matrix A, and the strictly lower */
 /* > triangular part of A is not referenced. If UPLO = 'L', the */
 /* > leading N-by-N lower triangular part of A contains the lower */
 /* > triangular part of the matrix A, and the strictly upper */
 /* > triangular part of A is not referenced. */
 /* > */
 /* > On exit, L is stored below (or above) the subdiaonal blocks, */
 /* > when UPLO is 'L' (or 'U'). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,N). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] TB */
 /* > \verbatim */
 /* > TB is COMPLEX array, dimension (LTB) */
 /* > On exit, details of the LU factorization of the band matrix. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LTB */
 /* > \verbatim */
 /* > LTB is INTEGER */
 /* > The size of the array TB. LTB >= 4*N, internally */
 /* > used to select NB such that LTB >= (3*NB+1)*N. */
 /* > */
 /* > If LTB = -1, then a workspace query is assumed;
 the */
 /* > routine only calculates the optimal size of LTB, */
 /* > returns this value as the first entry of TB, and */
 /* > no error message related to LTB is issued by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IPIV */
 /* > \verbatim */
 /* > IPIV is INTEGER array, dimension (N) */
 /* > On exit, it contains the details of the interchanges, i.e., */
 /* > the row and column k of A were interchanged with the */
 /* > row and column IPIV(k). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IPIV2 */
 /* > \verbatim */
 /* > IPIV2 is INTEGER array, dimension (N) */
 /* > On exit, it contains the details of the interchanges, i.e., */
 /* > the row and column k of T were interchanged with the */
 /* > row and column IPIV(k). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is COMPLEX workspace of size LWORK */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The size of WORK. LWORK >= N, internally used to select NB */
 /* > such that LWORK >= N*NB. */
 /* > */
 /* > If LWORK = -1, then a workspace query is assumed;
 the */
 /* > routine only calculates the optimal size of the WORK array, */
 /* > returns this value as the first entry of the WORK array, and */
 /* > no error message related to LWORK is issued by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value. */
 /* > > 0: if INFO = i, band LU factorization failed on i-th column */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2017 */
 /* > \ingroup complexSYcomputational */
 /* ===================================================================== */
 /* Subroutine */
 int csytrf_aa_2stage_(char *uplo, integer *n, complex *a, integer *lda, complex *tb, integer *ltb, integer *ipiv, integer * ipiv2, complex *work, integer *lwork, integer *info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
 char buffer[256]; 
#if FLA_ENABLE_ILP64 
 snprintf(buffer, 256,"csytrf_aa_2stage inputs: uplo %c, n %lld, lda %lld, ltb %lld, lwork %lld",*uplo, *n, *lda, *ltb, *lwork);
#else 
 snprintf(buffer, 256,"csytrf_aa_2stage inputs: uplo %c, n %d, lda %d, ltb %d, lwork %d",*uplo, *n, *lda, *ltb, *lwork);
#endif
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
 complex q__1;
 /* Local variables */
 integer i__, j, k, i1, i2, jb, kb, nb, td, nt;
 complex piv;
 integer ldtb;
 extern /* Subroutine */
 int cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
 extern logical lsame_(char *, char *);
 integer iinfo;
 extern /* Subroutine */
 int ccopy_(integer *, complex *, integer *, complex *, integer *), cswap_(integer *, complex *, integer *, complex *, integer *), ctrsm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *);
 logical upper;
 extern /* Subroutine */
 int cgbtrf_(integer *, integer *, integer *, integer *, complex *, integer *, integer *, integer *), cgetrf_( integer *, integer *, complex *, integer *, integer *, integer *), clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *), xerbla_( char *, integer *);
 extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
 logical tquery, wquery;
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
 /* Test the input parameters. */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 --tb;
 --ipiv;
 --ipiv2;
 --work;
 /* Function Body */
 *info = 0;
 upper = lsame_(uplo, "U");
 wquery = *lwork == -1;
 tquery = *ltb == -1;
 if (! upper && ! lsame_(uplo, "L")) {
 *info = -1;
 }
 else if (*n < 0) {
 *info = -2;
 }
 else if (*lda < max(1,*n)) {
 *info = -4;
 }
 else if (*ltb < *n << 2 && ! tquery) {
 *info = -6;
 }
 else if (*lwork < *n && ! wquery) {
 *info = -10;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("CSYTRF_AA_2STAGE", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Answer the query */
 nb = ilaenv_(&c__1, "CSYTRF_AA_2STAGE", uplo, n, &c_n1, &c_n1, &c_n1);
 if (*info == 0) {
 if (tquery) {
 i__1 = (nb * 3 + 1) * *n;
 tb[1].r = (real) i__1; tb[1].i = 0.f; // , expr subst  
 }
 if (wquery) {
 i__1 = *n * nb;
 work[1].r = (real) i__1; work[1].i = 0.f; // , expr subst  
 }
 }
 if (tquery || wquery) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Quick return */
 if (*n == 0) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Determine the number of the block size */
 ldtb = *ltb / *n;
 if (ldtb < nb * 3 + 1) {
 nb = (ldtb - 1) / 3;
 }
 if (*lwork < nb * *n) {
 nb = *lwork / *n;
 }
 /* Determine the number of the block columns */
 nt = (*n + nb - 1) / nb;
 td = nb << 1;
 kb = min(nb,*n);
 /* Initialize vectors/matrices */
 i__1 = kb;
 for (j = 1;
 j <= i__1;
 ++j) {
 ipiv[j] = j;
 }
 /* Save NB */
 tb[1].r = (real) nb; tb[1].i = 0.f; // , expr subst  
 if (upper) {
 /* ..................................................... */
 /* Factorize A as U**T*D*U using the upper triangle of A */
 /* ..................................................... */
 i__1 = nt - 1;
 for (j = 0;
 j <= i__1;
 ++j) {
 /* Generate Jth column of W and H */
 /* Computing MIN */
 i__2 = nb; i__3 = *n - j * nb; // , expr subst  
 kb = min(i__2,i__3);
 i__2 = j - 1;
 for (i__ = 1;
 i__ <= i__2;
 ++i__) {
 if (i__ == 1) {
 /* H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J) */
 if (i__ == j - 1) {
 jb = nb + kb;
 }
 else {
 jb = nb << 1;
 }
 i__3 = ldtb - 1;
 cgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2, &tb[td + 1 + i__ * nb * ldtb], &i__3, &a[(i__ - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &c_b1, &work[i__ * nb + 1], n);
 }
 else {
 /* H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J) */
 if (i__ == j - 1) {
 jb = (nb << 1) + kb;
 }
 else {
 jb = nb * 3;
 }
 i__3 = ldtb - 1;
 cgemm_("NoTranspose", "NoTranspose", &nb, &kb, &jb, &c_b2, &tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, &a[(i__ - 2) * nb + 1 + (j * nb + 1) * a_dim1], lda, &c_b1, &work[i__ * nb + 1], n);
 }
 }
 /* Compute T(J,J) */
 i__2 = ldtb - 1;
 clacpy_("Upper", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb * ldtb], &i__2);
 if (j > 1) {
 /* T(J,J) = U(1:J,J)'*H(1:J) */
 i__2 = (j - 1) * nb;
 q__1.r = -1.f; q__1.i = -0.f; // , expr subst  
 i__3 = ldtb - 1;
 cgemm_("Transpose", "NoTranspose", &kb, &kb, &i__2, &q__1, &a[ (j * nb + 1) * a_dim1 + 1], lda, &work[nb + 1], n, & c_b2, &tb[td + 1 + j * nb * ldtb], &i__3);
 /* T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J) */
 i__2 = ldtb - 1;
 cgemm_("Transpose", "NoTranspose", &kb, &nb, &kb, &c_b2, &a[( j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b1, &work[ 1], n);
 q__1.r = -1.f; q__1.i = -0.f; // , expr subst  
 i__2 = ldtb - 1;
 cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &nb, &q__1, & work[1], n, &a[(j - 2) * nb + 1 + (j * nb + 1) * a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], & i__2);
 }
 /* Expand T(J,J) into full format */
 i__2 = kb;
 for (i__ = 1;
 i__ <= i__2;
 ++i__) {
 i__3 = kb;
 for (k = i__ + 1;
 k <= i__3;
 ++k) {
 i__4 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
 i__5 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
 tb[i__4].r = tb[i__5].r; tb[i__4].i = tb[i__5].i; // , expr subst  
 }
 }
 if (j > 0) {
 /* CALL CHEGST( 1, 'Upper', KB, */
 /* $ TB( TD+1 + (J*NB)*LDTB ), LDTB-1, */
 /* $ A( (J-1)*NB+1, J*NB+1 ), LDA, IINFO ) */
 i__2 = ldtb - 1;
 ctrsm_("L", "U", "T", "N", &kb, &kb, &c_b2, &a[(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb * ldtb], &i__2);
 i__2 = ldtb - 1;
 ctrsm_("R", "U", "N", "N", &kb, &kb, &c_b2, &a[(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb * ldtb], &i__2);
 }
 if (j < nt - 1) {
 if (j > 0) {
 /* Compute H(J,J) */
 if (j == 1) {
 i__2 = ldtb - 1;
 cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &kb, & c_b2, &tb[td + 1 + j * nb * ldtb], &i__2, &a[( j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &c_b1, &work[j * nb + 1], n);
 }
 else {
 i__2 = nb + kb;
 i__3 = ldtb - 1;
 cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, &c_b2, &tb[td + nb + 1 + (j - 1) * nb * ldtb], &i__3, &a[(j - 2) * nb + 1 + (j * nb + 1) * a_dim1], lda, &c_b1, &work[j * nb + 1], n);
 }
 /* Update with the previous column */
 i__2 = *n - (j + 1) * nb;
 i__3 = j * nb;
 q__1.r = -1.f; q__1.i = -0.f; // , expr subst  
 cgemm_("Transpose", "NoTranspose", &nb, &i__2, &i__3, & q__1, &work[nb + 1], n, &a[((j + 1) * nb + 1) * a_dim1 + 1], lda, &c_b2, &a[j * nb + 1 + ((j + 1) * nb + 1) * a_dim1], lda);
 }
 /* Copy panel to workspace to call CGETRF */
 i__2 = nb;
 for (k = 1;
 k <= i__2;
 ++k) {
 i__3 = *n - (j + 1) * nb;
 ccopy_(&i__3, &a[j * nb + k + ((j + 1) * nb + 1) * a_dim1] , lda, &work[(k - 1) * *n + 1], &c__1);
 }
 /* Factorize panel */
 i__2 = *n - (j + 1) * nb;
 cgetrf_(&i__2, &nb, &work[1], n, &ipiv[(j + 1) * nb + 1], & iinfo);
 /* IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
 /* INFO = IINFO+(J+1)*NB */
 /* END IF */
 /* Copy panel back */
 i__2 = nb;
 for (k = 1;
 k <= i__2;
 ++k) {
 i__3 = *n - (j + 1) * nb;
 ccopy_(&i__3, &work[(k - 1) * *n + 1], &c__1, &a[j * nb + k + ((j + 1) * nb + 1) * a_dim1], lda);
 }
 /* Compute T(J+1, J), zero out for GEMM update */
 /* Computing MIN */
 i__2 = nb; i__3 = *n - (j + 1) * nb; // , expr subst  
 kb = min(i__2,i__3);
 i__2 = ldtb - 1;
 claset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * nb * ldtb], &i__2);
 i__2 = ldtb - 1;
 clacpy_("Upper", &kb, &nb, &work[1], n, &tb[td + nb + 1 + j * nb * ldtb], &i__2);
 if (j > 0) {
 i__2 = ldtb - 1;
 ctrsm_("R", "U", "N", "U", &kb, &nb, &c_b2, &a[(j - 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + nb + 1 + j * nb * ldtb], &i__2);
 }
 /* Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM */
 /* updates */
 i__2 = nb;
 for (k = 1;
 k <= i__2;
 ++k) {
 i__3 = kb;
 for (i__ = 1;
 i__ <= i__3;
 ++i__) {
 i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1) * ldtb;
 i__5 = td + nb + i__ - k + 1 + (j * nb + k - 1) * ldtb;
 tb[i__4].r = tb[i__5].r; tb[i__4].i = tb[i__5].i; // , expr subst  
 }
 }
 claset_("Lower", &kb, &nb, &c_b1, &c_b2, &a[j * nb + 1 + ((j + 1) * nb + 1) * a_dim1], lda);
 /* Apply pivots to trailing submatrix of A */
 i__2 = kb;
 for (k = 1;
 k <= i__2;
 ++k) {
 /* > Adjust ipiv */
 ipiv[(j + 1) * nb + k] += (j + 1) * nb;
 i1 = (j + 1) * nb + k;
 i2 = ipiv[(j + 1) * nb + k];
 if (i1 != i2) {
 /* > Apply pivots to previous columns of L */
 i__3 = k - 1;
 cswap_(&i__3, &a[(j + 1) * nb + 1 + i1 * a_dim1], & c__1, &a[(j + 1) * nb + 1 + i2 * a_dim1], & c__1);
 /* > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
 if (i2 > i1 + 1) {
 i__3 = i2 - i1 - 1;
 cswap_(&i__3, &a[i1 + (i1 + 1) * a_dim1], lda, &a[ i1 + 1 + i2 * a_dim1], &c__1);
 }
 /* > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
 if (i2 < *n) {
 i__3 = *n - i2;
 cswap_(&i__3, &a[i1 + (i2 + 1) * a_dim1], lda, &a[ i2 + (i2 + 1) * a_dim1], lda);
 }
 /* > Swap A(I1, I1) with A(I2, I2) */
 i__3 = i1 + i1 * a_dim1;
 piv.r = a[i__3].r; piv.i = a[i__3].i; // , expr subst  
 i__3 = i1 + i1 * a_dim1;
 i__4 = i2 + i2 * a_dim1;
 a[i__3].r = a[i__4].r; a[i__3].i = a[i__4].i; // , expr subst  
 i__3 = i2 + i2 * a_dim1;
 a[i__3].r = piv.r; a[i__3].i = piv.i; // , expr subst  
 /* > Apply pivots to previous columns of L */
 if (j > 0) {
 i__3 = j * nb;
 cswap_(&i__3, &a[i1 * a_dim1 + 1], &c__1, &a[i2 * a_dim1 + 1], &c__1);
 }
 }
 }
 }
 }
 }
 else {
 /* ..................................................... */
 /* Factorize A as L*D*L**T using the lower triangle of A */
 /* ..................................................... */
 i__1 = nt - 1;
 for (j = 0;
 j <= i__1;
 ++j) {
 /* Generate Jth column of W and H */
 /* Computing MIN */
 i__2 = nb; i__3 = *n - j * nb; // , expr subst  
 kb = min(i__2,i__3);
 i__2 = j - 1;
 for (i__ = 1;
 i__ <= i__2;
 ++i__) {
 if (i__ == 1) {
 /* H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)' */
 if (i__ == j - 1) {
 jb = nb + kb;
 }
 else {
 jb = nb << 1;
 }
 i__3 = ldtb - 1;
 cgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b2, & tb[td + 1 + i__ * nb * ldtb], &i__3, &a[j * nb + 1 + ((i__ - 1) * nb + 1) * a_dim1], lda, &c_b1, & work[i__ * nb + 1], n);
 }
 else {
 /* H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)' */
 if (i__ == j - 1) {
 jb = (nb << 1) + kb;
 }
 else {
 jb = nb * 3;
 }
 i__3 = ldtb - 1;
 cgemm_("NoTranspose", "Transpose", &nb, &kb, &jb, &c_b2, & tb[td + nb + 1 + (i__ - 1) * nb * ldtb], &i__3, & a[j * nb + 1 + ((i__ - 2) * nb + 1) * a_dim1], lda, &c_b1, &work[i__ * nb + 1], n);
 }
 }
 /* Compute T(J,J) */
 i__2 = ldtb - 1;
 clacpy_("Lower", &kb, &kb, &a[j * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb * ldtb], &i__2);
 if (j > 1) {
 /* T(J,J) = L(J,1:J)*H(1:J) */
 i__2 = (j - 1) * nb;
 q__1.r = -1.f; q__1.i = -0.f; // , expr subst  
 i__3 = ldtb - 1;
 cgemm_("NoTranspose", "NoTranspose", &kb, &kb, &i__2, &q__1, & a[j * nb + 1 + a_dim1], lda, &work[nb + 1], n, &c_b2, &tb[td + 1 + j * nb * ldtb], &i__3);
 /* T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)' */
 i__2 = ldtb - 1;
 cgemm_("NoTranspose", "NoTranspose", &kb, &nb, &kb, &c_b2, &a[ j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[ td + nb + 1 + (j - 1) * nb * ldtb], &i__2, &c_b1, & work[1], n);
 q__1.r = -1.f; q__1.i = -0.f; // , expr subst  
 i__2 = ldtb - 1;
 cgemm_("NoTranspose", "Transpose", &kb, &kb, &nb, &q__1, & work[1], n, &a[j * nb + 1 + ((j - 2) * nb + 1) * a_dim1], lda, &c_b2, &tb[td + 1 + j * nb * ldtb], & i__2);
 }
 /* Expand T(J,J) into full format */
 i__2 = kb;
 for (i__ = 1;
 i__ <= i__2;
 ++i__) {
 i__3 = kb;
 for (k = i__ + 1;
 k <= i__3;
 ++k) {
 i__4 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
 i__5 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
 tb[i__4].r = tb[i__5].r; tb[i__4].i = tb[i__5].i; // , expr subst  
 }
 }
 if (j > 0) {
 /* CALL CHEGST( 1, 'Lower', KB, */
 /* $ TB( TD+1 + (J*NB)*LDTB ), LDTB-1, */
 /* $ A( J*NB+1, (J-1)*NB+1 ), LDA, IINFO ) */
 i__2 = ldtb - 1;
 ctrsm_("L", "L", "N", "N", &kb, &kb, &c_b2, &a[j * nb + 1 + (( j - 1) * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb * ldtb], &i__2);
 i__2 = ldtb - 1;
 ctrsm_("R", "L", "T", "N", &kb, &kb, &c_b2, &a[j * nb + 1 + (( j - 1) * nb + 1) * a_dim1], lda, &tb[td + 1 + j * nb * ldtb], &i__2);
 }
 /* Symmetrize T(J,J) */
 i__2 = kb;
 for (i__ = 1;
 i__ <= i__2;
 ++i__) {
 i__3 = kb;
 for (k = i__ + 1;
 k <= i__3;
 ++k) {
 i__4 = td - (k - (i__ + 1)) + (j * nb + k - 1) * ldtb;
 i__5 = td + (k - i__) + 1 + (j * nb + i__ - 1) * ldtb;
 tb[i__4].r = tb[i__5].r; tb[i__4].i = tb[i__5].i; // , expr subst  
 }
 }
 if (j < nt - 1) {
 if (j > 0) {
 /* Compute H(J,J) */
 if (j == 1) {
 i__2 = ldtb - 1;
 cgemm_("NoTranspose", "Transpose", &kb, &kb, &kb, & c_b2, &tb[td + 1 + j * nb * ldtb], &i__2, &a[ j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], lda, &c_b1, &work[j * nb + 1], n);
 }
 else {
 i__2 = nb + kb;
 i__3 = ldtb - 1;
 cgemm_("NoTranspose", "Transpose", &kb, &kb, &i__2, & c_b2, &tb[td + nb + 1 + (j - 1) * nb * ldtb], &i__3, &a[j * nb + 1 + ((j - 2) * nb + 1) * a_dim1], lda, &c_b1, &work[j * nb + 1], n);
 }
 /* Update with the previous column */
 i__2 = *n - (j + 1) * nb;
 i__3 = j * nb;
 q__1.r = -1.f; q__1.i = -0.f; // , expr subst  
 cgemm_("NoTranspose", "NoTranspose", &i__2, &nb, &i__3, & q__1, &a[(j + 1) * nb + 1 + a_dim1], lda, &work[ nb + 1], n, &c_b2, &a[(j + 1) * nb + 1 + (j * nb + 1) * a_dim1], lda);
 }
 /* Factorize panel */
 i__2 = *n - (j + 1) * nb;
 cgetrf_(&i__2, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &ipiv[(j + 1) * nb + 1], &iinfo);
 /* IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN */
 /* INFO = IINFO+(J+1)*NB */
 /* END IF */
 /* Compute T(J+1, J), zero out for GEMM update */
 /* Computing MIN */
 i__2 = nb; i__3 = *n - (j + 1) * nb; // , expr subst  
 kb = min(i__2,i__3);
 i__2 = ldtb - 1;
 claset_("Full", &kb, &nb, &c_b1, &c_b1, &tb[td + nb + 1 + j * nb * ldtb], &i__2);
 i__2 = ldtb - 1;
 clacpy_("Upper", &kb, &nb, &a[(j + 1) * nb + 1 + (j * nb + 1) * a_dim1], lda, &tb[td + nb + 1 + j * nb * ldtb], & i__2);
 if (j > 0) {
 i__2 = ldtb - 1;
 ctrsm_("R", "L", "T", "U", &kb, &nb, &c_b2, &a[j * nb + 1 + ((j - 1) * nb + 1) * a_dim1], lda, &tb[td + nb + 1 + j * nb * ldtb], &i__2);
 }
 /* Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM */
 /* updates */
 i__2 = nb;
 for (k = 1;
 k <= i__2;
 ++k) {
 i__3 = kb;
 for (i__ = 1;
 i__ <= i__3;
 ++i__) {
 i__4 = td - nb + k - i__ + 1 + (j * nb + nb + i__ - 1) * ldtb;
 i__5 = td + nb + i__ - k + 1 + (j * nb + k - 1) * ldtb;
 tb[i__4].r = tb[i__5].r; tb[i__4].i = tb[i__5].i; // , expr subst  
 }
 }
 claset_("Upper", &kb, &nb, &c_b1, &c_b2, &a[(j + 1) * nb + 1 + (j * nb + 1) * a_dim1], lda);
 /* Apply pivots to trailing submatrix of A */
 i__2 = kb;
 for (k = 1;
 k <= i__2;
 ++k) {
 /* > Adjust ipiv */
 ipiv[(j + 1) * nb + k] += (j + 1) * nb;
 i1 = (j + 1) * nb + k;
 i2 = ipiv[(j + 1) * nb + k];
 if (i1 != i2) {
 /* > Apply pivots to previous columns of L */
 i__3 = k - 1;
 cswap_(&i__3, &a[i1 + ((j + 1) * nb + 1) * a_dim1], lda, &a[i2 + ((j + 1) * nb + 1) * a_dim1], lda);
 /* > Swap A(I1+1:M, I1) with A(I2, I1+1:M) */
 if (i2 > i1 + 1) {
 i__3 = i2 - i1 - 1;
 cswap_(&i__3, &a[i1 + 1 + i1 * a_dim1], &c__1, &a[ i2 + (i1 + 1) * a_dim1], lda);
 }
 /* > Swap A(I2+1:M, I1) with A(I2+1:M, I2) */
 if (i2 < *n) {
 i__3 = *n - i2;
 cswap_(&i__3, &a[i2 + 1 + i1 * a_dim1], &c__1, &a[ i2 + 1 + i2 * a_dim1], &c__1);
 }
 /* > Swap A(I1, I1) with A(I2, I2) */
 i__3 = i1 + i1 * a_dim1;
 piv.r = a[i__3].r; piv.i = a[i__3].i; // , expr subst  
 i__3 = i1 + i1 * a_dim1;
 i__4 = i2 + i2 * a_dim1;
 a[i__3].r = a[i__4].r; a[i__3].i = a[i__4].i; // , expr subst  
 i__3 = i2 + i2 * a_dim1;
 a[i__3].r = piv.r; a[i__3].i = piv.i; // , expr subst  
 /* > Apply pivots to previous columns of L */
 if (j > 0) {
 i__3 = j * nb;
 cswap_(&i__3, &a[i1 + a_dim1], lda, &a[i2 + a_dim1], lda);
 }
 }
 }
 /* Apply pivots to previous columns of L */
 /* CALL CLASWP( J*NB, A( 1, 1 ), LDA, */
 /* $ (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 ) */
 }
 }
 }
 /* Factor the band matrix */
 cgbtrf_(n, n, &nb, &nb, &tb[1], &ldtb, &ipiv2[1], info);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of CSYTRF_AA_2STAGE */
 }
 /* csytrf_aa_2stage__ */
 