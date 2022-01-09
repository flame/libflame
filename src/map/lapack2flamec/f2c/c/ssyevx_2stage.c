/* ../netlib/v3.9.0/ssyevx_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 static integer c_n1 = -1;
 static integer c__2 = 2;
 static integer c__3 = 3;
 static integer c__4 = 4;
 /* > \brief <b> SSYEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b> */
 /* @generated from dsyevx_2stage.f, fortran d -> s, Sat Nov 5 23:55:46 2016 */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download SSYEVX_2STAGE + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevx_ 2stage.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevx_ 2stage.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevx_ 2stage.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE SSYEVX_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, */
 /* IL, IU, ABSTOL, M, W, Z, LDZ, WORK, */
 /* LWORK, IWORK, IFAIL, INFO ) */
 /* IMPLICIT NONE */
 /* .. Scalar Arguments .. */
 /* CHARACTER JOBZ, RANGE, UPLO */
 /* INTEGER IL, INFO, IU, LDA, LDZ, LWORK, M, N */
 /* REAL ABSTOL, VL, VU */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IFAIL( * ), IWORK( * ) */
 /* REAL A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > SSYEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
 /* > of a real symmetric matrix A using the 2stage technique for */
 /* > the reduction to tridiagonal. Eigenvalues and eigenvectors can be */
 /* > selected by specifying either a range of values or a range of indices */
 /* > for the desired eigenvalues. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] JOBZ */
 /* > \verbatim */
 /* > JOBZ is CHARACTER*1 */
 /* > = 'N': Compute eigenvalues only;
 */
 /* > = 'V': Compute eigenvalues and eigenvectors. */
 /* > Not available in this release. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] RANGE */
 /* > \verbatim */
 /* > RANGE is CHARACTER*1 */
 /* > = 'A': all eigenvalues will be found. */
 /* > = 'V': all eigenvalues in the half-open interval (VL,VU] */
 /* > will be found. */
 /* > = 'I': the IL-th through IU-th eigenvalues will be found. */
 /* > \endverbatim */
 /* > */
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
 /* > A is REAL array, dimension (LDA, N) */
 /* > On entry, the symmetric matrix A. If UPLO = 'U', the */
 /* > leading N-by-N upper triangular part of A contains the */
 /* > upper triangular part of the matrix A. If UPLO = 'L', */
 /* > the leading N-by-N lower triangular part of A contains */
 /* > the lower triangular part of the matrix A. */
 /* > On exit, the lower triangle (if UPLO='L') or the upper */
 /* > triangle (if UPLO='U') of A, including the diagonal, is */
 /* > destroyed. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,N). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] VL */
 /* > \verbatim */
 /* > VL is REAL */
 /* > If RANGE='V', the lower bound of the interval to */
 /* > be searched for eigenvalues. VL < VU. */
 /* > Not referenced if RANGE = 'A' or 'I'. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] VU */
 /* > \verbatim */
 /* > VU is REAL */
 /* > If RANGE='V', the upper bound of the interval to */
 /* > be searched for eigenvalues. VL < VU. */
 /* > Not referenced if RANGE = 'A' or 'I'. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] IL */
 /* > \verbatim */
 /* > IL is INTEGER */
 /* > If RANGE='I', the index of the */
 /* > smallest eigenvalue to be returned. */
 /* > 1 <= IL <= IU <= N, if N > 0;
 IL = 1 and IU = 0 if N = 0. */
 /* > Not referenced if RANGE = 'A' or 'V'. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] IU */
 /* > \verbatim */
 /* > IU is INTEGER */
 /* > If RANGE='I', the index of the */
 /* > largest eigenvalue to be returned. */
 /* > 1 <= IL <= IU <= N, if N > 0;
 IL = 1 and IU = 0 if N = 0. */
 /* > Not referenced if RANGE = 'A' or 'V'. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] ABSTOL */
 /* > \verbatim */
 /* > ABSTOL is REAL */
 /* > The absolute error tolerance for the eigenvalues. */
 /* > An approximate eigenvalue is accepted as converged */
 /* > when it is determined to lie in an interval [a,b] */
 /* > of width less than or equal to */
 /* > */
 /* > ABSTOL + EPS * max( |a|,|b| ) , */
 /* > */
 /* > where EPS is the machine precision. If ABSTOL is less than */
 /* > or equal to zero, then EPS*|T| will be used in its place, */
 /* > where |T| is the 1-norm of the tridiagonal matrix obtained */
 /* > by reducing A to tridiagonal form. */
 /* > */
 /* > Eigenvalues will be computed most accurately when ABSTOL is */
 /* > set to twice the underflow threshold 2*SLAMCH('S'), not zero. */
 /* > If this routine returns with INFO>0, indicating that some */
 /* > eigenvectors did not converge, try setting ABSTOL to */
 /* > 2*SLAMCH('S'). */
 /* > */
 /* > See "Computing Small Singular Values of Bidiagonal Matrices */
 /* > with Guaranteed High Relative Accuracy," by Demmel and */
 /* > Kahan, LAPACK Working Note #3. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] M */
 /* > \verbatim */
 /* > M is INTEGER */
 /* > The total number of eigenvalues found. 0 <= M <= N. */
 /* > If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] W */
 /* > \verbatim */
 /* > W is REAL array, dimension (N) */
 /* > On normal exit, the first M elements contain the selected */
 /* > eigenvalues in ascending order. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] Z */
 /* > \verbatim */
 /* > Z is REAL array, dimension (LDZ, max(1,M)) */
 /* > If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
 /* > contain the orthonormal eigenvectors of the matrix A */
 /* > corresponding to the selected eigenvalues, with the i-th */
 /* > column of Z holding the eigenvector associated with W(i). */
 /* > If an eigenvector fails to converge, then that column of Z */
 /* > contains the latest approximation to the eigenvector, and the */
 /* > index of the eigenvector is returned in IFAIL. */
 /* > If JOBZ = 'N', then Z is not referenced. */
 /* > Note: the user must ensure that at least max(1,M) columns are */
 /* > supplied in the array Z;
 if RANGE = 'V', the exact value of M */
 /* > is not known in advance and an upper bound must be used. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDZ */
 /* > \verbatim */
 /* > LDZ is INTEGER */
 /* > The leading dimension of the array Z. LDZ >= 1, and if */
 /* > JOBZ = 'V', LDZ >= max(1,N). */
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
 /* > The length of the array WORK. LWORK >= 1, when N <= 1;
 */
 /* > otherwise */
 /* > If JOBZ = 'N' and N > 1, LWORK must be queried. */
 /* > LWORK = MAX(1, 8*N, dimension) where */
 /* > dimension = max(stage1,stage2) + (KD+1)*N + 3*N */
 /* > = N*KD + N*max(KD+1,FACTOPTNB) */
 /* > + max(2*KD*KD, KD*NTHREADS) */
 /* > + (KD+1)*N + 3*N */
 /* > where KD is the blocking size of the reduction, */
 /* > FACTOPTNB is the blocking used by the QR or LQ */
 /* > algorithm, usually FACTOPTNB=128 is a good choice */
 /* > NTHREADS is the number of threads used when */
 /* > openMP compilation is enabled, otherwise =1. */
 /* > If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available */
 /* > */
 /* > If LWORK = -1, then a workspace query is assumed;
 the routine */
 /* > only calculates the optimal size of the WORK array, returns */
 /* > this value as the first entry of the WORK array, and no error */
 /* > message related to LWORK is issued by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IWORK */
 /* > \verbatim */
 /* > IWORK is INTEGER array, dimension (5*N) */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IFAIL */
 /* > \verbatim */
 /* > IFAIL is INTEGER array, dimension (N) */
 /* > If JOBZ = 'V', then if INFO = 0, the first M elements of */
 /* > IFAIL are zero. If INFO > 0, then IFAIL contains the */
 /* > indices of the eigenvectors that failed to converge. */
 /* > If JOBZ = 'N', then IFAIL is not referenced. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value */
 /* > > 0: if INFO = i, then i eigenvectors failed to converge. */
 /* > Their indices are stored in array IFAIL. */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date June 2016 */
 /* > \ingroup realSYeigen */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > All details about the 2stage techniques are available in: */
 /* > */
 /* > Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
 /* > Parallel reduction to condensed forms for symmetric eigenvalue problems */
 /* > using aggregated fine-grained and memory-aware kernels. In Proceedings */
 /* > of 2011 International Conference for High Performance Computing, */
 /* > Networking, Storage and Analysis (SC '11), New York, NY, USA, */
 /* > Article 8 , 11 pages. */
 /* > http://doi.acm.org/10.1145/2063384.2063394 */
 /* > */
 /* > A. Haidar, J. Kurzak, P. Luszczek, 2013. */
 /* > An improved parallel singular value algorithm and its implementation */
 /* > for multicore hardware, In Proceedings of 2013 International Conference */
 /* > for High Performance Computing, Networking, Storage and Analysis (SC '13). */
 /* > Denver, Colorado, USA, 2013. */
 /* > Article 90, 12 pages. */
 /* > http://doi.acm.org/10.1145/2503210.2503292 */
 /* > */
 /* > A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
 /* > A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
 /* > calculations based on fine-grained memory aware tasks. */
 /* > International Journal of High Performance Computing Applications. */
 /* > Volume 28 Issue 2, Pages 196-209, May 2014. */
 /* > http://hpc.sagepub.com/content/28/2/196 */
 /* > */
 /* > \endverbatim */
 /* ===================================================================== */
 /* Subroutine */
 int ssyevx_2stage_(char *jobz, char *range, char *uplo, integer *n, real *a, integer *lda, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, real *z__, integer * ldz, real *work, integer *lwork, integer *iwork, integer *ifail, integer *info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
 char buffer[256];
 snprintf(buffer, 256,"ssyevx_2stage inputs: jobz %c, range %c, uplo %c, n %" FLA_IS ", lda %" FLA_IS ", il %" FLA_IS ", iu %" FLA_IS ", ldz %" FLA_IS "",*jobz, *range, *uplo, *n, *lda, *il, *iu, *ldz);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
 real r__1, r__2;
 /* Builtin functions */
 double sqrt(doublereal);
 /* Local variables */
 integer i__, j, ib, kd, jj;
 real eps, vll, vuu, tmp1;
 integer indd, inde;
 extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
 real anrm;
 integer imax;
 real rmin, rmax;
 logical test;
 integer itmp1, indee;
 real sigma;
 extern logical lsame_(char *, char *);
 integer iinfo;
 extern /* Subroutine */
 int sscal_(integer *, real *, real *, integer *);
 char order[1];
 integer lhtrd, lwmin;
 logical lower;
 integer lwtrd;
 extern /* Subroutine */
 int scopy_(integer *, real *, integer *, real *, integer *), sswap_(integer *, real *, integer *, real *, integer * ), ssytrd_2stage_(char *, char *, integer *, real *, integer *, real *, real *, real *, real *, integer *, real *, integer *, integer *);
 logical wantz, alleig, indeig;
 integer iscale, indibl;
 logical valeig;
 extern real slamch_(char *);
 real safmin;
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 real abstll, bignum;
 integer indtau, indisp, indiwo, indwkn;
 extern /* Subroutine */
 int slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *);
 integer indwrk;
 extern /* Subroutine */
 int sstein_(integer *, real *, real *, integer *, real *, integer *, integer *, real *, integer *, real *, integer * , integer *, integer *), ssterf_(integer *, real *, real *, integer *);
 integer llwrkn, llwork, nsplit;
 real smlnum;
 extern real slansy_(char *, char *, integer *, real *, integer *, real *);
 extern /* Subroutine */
 int sstebz_(char *, char *, integer *, real *, real *, integer *, integer *, real *, real *, real *, integer *, integer *, real *, integer *, integer *, real *, integer *, integer *), sorgtr_(char *, integer *, real *, integer *, real *, real *, integer *, integer *);
 logical lquery;
 extern /* Subroutine */
 int ssteqr_(char *, integer *, real *, real *, real *, integer *, real *, integer *), sormtr_(char *, char *, char *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *);
 integer indhous;
 /* -- LAPACK driver routine (version 3.8.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* June 2016 */
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
 /* Test the input parameters. */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 --w;
 z_dim1 = *ldz;
 z_offset = 1 + z_dim1;
 z__ -= z_offset;
 --work;
 --iwork;
 --ifail;
 /* Function Body */
 lower = lsame_(uplo, "L");
 wantz = lsame_(jobz, "V");
 alleig = lsame_(range, "A");
 valeig = lsame_(range, "V");
 indeig = lsame_(range, "I");
 lquery = *lwork == -1;
 *info = 0;
 if (! lsame_(jobz, "N")) {
 *info = -1;
 }
 else if (! (alleig || valeig || indeig)) {
 *info = -2;
 }
 else if (! (lower || lsame_(uplo, "U"))) {
 *info = -3;
 }
 else if (*n < 0) {
 *info = -4;
 }
 else if (*lda < max(1,*n)) {
 *info = -6;
 }
 else {
 if (valeig) {
 if (*n > 0 && *vu <= *vl) {
 *info = -8;
 }
 }
 else if (indeig) {
 if (*il < 1 || *il > max(1,*n)) {
 *info = -9;
 }
 else if (*iu < min(*n,*il) || *iu > *n) {
 *info = -10;
 }
 }
 }
 if (*info == 0) {
 if (*ldz < 1 || wantz && *ldz < *n) {
 *info = -15;
 }
 }
 if (*info == 0) {
 if (*n <= 1) {
 lwmin = 1;
 work[1] = (real) lwmin;
 }
 else {
 kd = ilaenv2stage_(&c__1, "SSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1);
 ib = ilaenv2stage_(&c__2, "SSYTRD_2STAGE", jobz, n, &kd, &c_n1, & c_n1);
 lhtrd = ilaenv2stage_(&c__3, "SSYTRD_2STAGE", jobz, n, &kd, &ib, & c_n1);
 lwtrd = ilaenv2stage_(&c__4, "SSYTRD_2STAGE", jobz, n, &kd, &ib, & c_n1);
 /* Computing MAX */
 i__1 = *n << 3; i__2 = *n * 3 + lhtrd + lwtrd; // , expr subst  
 lwmin = max(i__1,i__2);
 work[1] = (real) lwmin;
 }
 if (*lwork < lwmin && ! lquery) {
 *info = -17;
 }
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("SSYEVX_2STAGE", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 else if (lquery) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Quick return if possible */
 *m = 0;
 if (*n == 0) {
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 if (*n == 1) {
 if (alleig || indeig) {
 *m = 1;
 w[1] = a[a_dim1 + 1];
 }
 else {
 if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
 *m = 1;
 w[1] = a[a_dim1 + 1];
 }
 }
 if (wantz) {
 z__[z_dim1 + 1] = 1.f;
 }
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Get machine constants. */
 safmin = slamch_("Safe minimum");
 eps = slamch_("Precision");
 smlnum = safmin / eps;
 bignum = 1.f / smlnum;
 rmin = sqrt(smlnum);
 /* Computing MIN */
 r__1 = sqrt(bignum); r__2 = 1.f / sqrt(sqrt(safmin)); // , expr subst  
 rmax = min(r__1,r__2);
 /* Scale matrix to allowable range, if necessary. */
 iscale = 0;
 abstll = *abstol;
 if (valeig) {
 vll = *vl;
 vuu = *vu;
 }
 anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1]);
 if (anrm > 0.f && anrm < rmin) {
 iscale = 1;
 sigma = rmin / anrm;
 }
 else if (anrm > rmax) {
 iscale = 1;
 sigma = rmax / anrm;
 }
 if (iscale == 1) {
 if (lower) {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 i__2 = *n - j + 1;
 sscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
 /* L10: */
 }
 }
 else {
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 sscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
 /* L20: */
 }
 }
 if (*abstol > 0.f) {
 abstll = *abstol * sigma;
 }
 if (valeig) {
 vll = *vl * sigma;
 vuu = *vu * sigma;
 }
 }
 /* Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */
 indtau = 1;
 inde = indtau + *n;
 indd = inde + *n;
 indhous = indd + *n;
 indwrk = indhous + lhtrd;
 llwork = *lwork - indwrk + 1;
 ssytrd_2stage_(jobz, uplo, n, &a[a_offset], lda, &work[indd], &work[inde] , &work[indtau], &work[indhous], &lhtrd, &work[indwrk], &llwork, & iinfo);
 /* If all eigenvalues are desired and ABSTOL is less than or equal to */
 /* zero, then call SSTERF or SORGTR and SSTEQR. If this fails for */
 /* some eigenvalue, then try SSTEBZ. */
 test = FALSE_;
 if (indeig) {
 if (*il == 1 && *iu == *n) {
 test = TRUE_;
 }
 }
 if ((alleig || test) && *abstol <= 0.f) {
 scopy_(n, &work[indd], &c__1, &w[1], &c__1);
 indee = indwrk + (*n << 1);
 if (! wantz) {
 i__1 = *n - 1;
 scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
 ssterf_(n, &w[1], &work[indee], info);
 }
 else {
 slacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz);
 sorgtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk] , &llwork, &iinfo);
 i__1 = *n - 1;
 scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
 ssteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[ indwrk], info);
 if (*info == 0) {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 ifail[i__] = 0;
 /* L30: */
 }
 }
 }
 if (*info == 0) {
 *m = *n;
 goto L40;
 }
 *info = 0;
 }
 /* Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */
 if (wantz) {
 *(unsigned char *)order = 'B';
 }
 else {
 *(unsigned char *)order = 'E';
 }
 indibl = 1;
 indisp = indibl + *n;
 indiwo = indisp + *n;
 sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[ inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[ indwrk], &iwork[indiwo], info);
 if (wantz) {
 sstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[ indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], & ifail[1], info);
 /* Apply orthogonal matrix used in reduction to tridiagonal */
 /* form to eigenvectors returned by SSTEIN. */
 indwkn = inde;
 llwrkn = *lwork - indwkn + 1;
 sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[ z_offset], ldz, &work[indwkn], &llwrkn, &iinfo);
 }
 /* If matrix was scaled, then rescale eigenvalues appropriately. */
 L40: if (iscale == 1) {
 if (*info == 0) {
 imax = *m;
 }
 else {
 imax = *info - 1;
 }
 r__1 = 1.f / sigma;
 sscal_(&imax, &r__1, &w[1], &c__1);
 }
 /* If eigenvalues are not in order, then sort them, along with */
 /* eigenvectors. */
 if (wantz) {
 i__1 = *m - 1;
 for (j = 1;
 j <= i__1;
 ++j) {
 i__ = 0;
 tmp1 = w[j];
 i__2 = *m;
 for (jj = j + 1;
 jj <= i__2;
 ++jj) {
 if (w[jj] < tmp1) {
 i__ = jj;
 tmp1 = w[jj];
 }
 /* L50: */
 }
 if (i__ != 0) {
 itmp1 = iwork[indibl + i__ - 1];
 w[i__] = w[j];
 iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
 w[j] = tmp1;
 iwork[indibl + j - 1] = itmp1;
 sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
 if (*info != 0) {
 itmp1 = ifail[i__];
 ifail[i__] = ifail[j];
 ifail[j] = itmp1;
 }
 }
 /* L60: */
 }
 }
 /* Set WORK(1) to optimal workspace size. */
 work[1] = (real) lwmin;
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of SSYEVX_2STAGE */
 }
 /* ssyevx_2stage__ */
 
