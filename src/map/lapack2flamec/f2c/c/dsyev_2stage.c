/* ../netlib/v3.9.0/dsyev_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 static integer c_n1 = -1;
 static integer c__2 = 2;
 static integer c__3 = 3;
 static integer c__4 = 4;
 static integer c__0 = 0;
 static doublereal c_b27 = 1.;
 /* > \brief <b> DSYEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b> */
 /* @precisions fortran d -> s */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DSYEV_2STAGE + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyevd_ 2stage.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyevd_ 2stage.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyevd_ 2stage.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE DSYEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, */
 /* INFO ) */
 /* IMPLICIT NONE */
 /* .. Scalar Arguments .. */
 /* CHARACTER JOBZ, UPLO */
 /* INTEGER INFO, LDA, LWORK, N */
 /* .. */
 /* .. Array Arguments .. */
 /* DOUBLE PRECISION A( LDA, * ), W( * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DSYEV_2STAGE computes all eigenvalues and, optionally, eigenvectors of a */
 /* > real symmetric matrix A using the 2stage technique for */
 /* > the reduction to tridiagonal. */
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
 /* > A is DOUBLE PRECISION array, dimension (LDA, N) */
 /* > On entry, the symmetric matrix A. If UPLO = 'U', the */
 /* > leading N-by-N upper triangular part of A contains the */
 /* > upper triangular part of the matrix A. If UPLO = 'L', */
 /* > the leading N-by-N lower triangular part of A contains */
 /* > the lower triangular part of the matrix A. */
 /* > On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
 /* > orthonormal eigenvectors of the matrix A. */
 /* > If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') */
 /* > or the upper triangle (if UPLO='U') of A, including the */
 /* > diagonal, is destroyed. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,N). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] W */
 /* > \verbatim */
 /* > W is DOUBLE PRECISION array, dimension (N) */
 /* > If INFO = 0, the eigenvalues in ascending order. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is DOUBLE PRECISION array, dimension LWORK */
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
 /* > LWORK = MAX(1, dimension) where */
 /* > dimension = max(stage1,stage2) + (KD+1)*N + 2*N */
 /* > = N*KD + N*max(KD+1,FACTOPTNB) */
 /* > + max(2*KD*KD, KD*NTHREADS) */
 /* > + (KD+1)*N + 2*N */
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
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value */
 /* > > 0: if INFO = i, the algorithm failed to converge;
 i */
 /* > off-diagonal elements of an intermediate tridiagonal */
 /* > form did not converge to zero. */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2017 */
 /* > \ingroup doubleSYeigen */
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
 int dsyev_2stage_(char *jobz, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *w, doublereal *work, integer *lwork, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, i__1;
 doublereal d__1;
 /* Builtin functions */
 double sqrt(doublereal);
 /* Local variables */
 integer ib, kd;
 doublereal eps;
 integer inde;
 extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
 integer imax;
 doublereal anrm, rmin, rmax;
 extern /* Subroutine */
 int dscal_(integer *, doublereal *, doublereal *, integer *);
 doublereal sigma;
 extern logical lsame_(char *, char *);
 integer iinfo;
 extern /* Subroutine */
 int dsytrd_2stage_(char *, char *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
 integer lhtrd, lwmin;
 logical lower;
 integer lwtrd;
 logical wantz;
 extern doublereal dlamch_(char *);
 integer iscale;
 extern /* Subroutine */
 int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
 doublereal safmin;
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 doublereal bignum;
 integer indtau;
 extern /* Subroutine */
 int dsterf_(integer *, doublereal *, doublereal *, integer *);
 extern doublereal dlansy_(char *, char *, integer *, doublereal *, integer *, doublereal *);
 integer indwrk;
 extern /* Subroutine */
 int dorgtr_(char *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, integer *), dsteqr_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *);
 integer llwork;
 doublereal smlnum;
 logical lquery;
 integer indhous;
 /* -- LAPACK driver routine (version 3.8.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* November 2017 */
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
 --work;
 /* Function Body */
 wantz = lsame_(jobz, "V");
 lower = lsame_(uplo, "L");
 lquery = *lwork == -1;
 *info = 0;
 if (! lsame_(jobz, "N")) {
 *info = -1;
 }
 else if (! (lower || lsame_(uplo, "U"))) {
 *info = -2;
 }
 else if (*n < 0) {
 *info = -3;
 }
 else if (*lda < max(1,*n)) {
 *info = -5;
 }
 if (*info == 0) {
 kd = ilaenv2stage_(&c__1, "DSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, & c_n1);
 ib = ilaenv2stage_(&c__2, "DSYTRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1);
 lhtrd = ilaenv2stage_(&c__3, "DSYTRD_2STAGE", jobz, n, &kd, &ib, & c_n1);
 lwtrd = ilaenv2stage_(&c__4, "DSYTRD_2STAGE", jobz, n, &kd, &ib, & c_n1);
 lwmin = (*n << 1) + lhtrd + lwtrd;
 work[1] = (doublereal) lwmin;
 if (*lwork < lwmin && ! lquery) {
 *info = -8;
 }
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("DSYEV_2STAGE ", &i__1);
 return 0;
 }
 else if (lquery) {
 return 0;
 }
 /* Quick return if possible */
 if (*n == 0) {
 return 0;
 }
 if (*n == 1) {
 w[1] = a[a_dim1 + 1];
 work[1] = 2.;
 if (wantz) {
 a[a_dim1 + 1] = 1.;
 }
 return 0;
 }
 /* Get machine constants. */
 safmin = dlamch_("Safe minimum");
 eps = dlamch_("Precision");
 smlnum = safmin / eps;
 bignum = 1. / smlnum;
 rmin = sqrt(smlnum);
 rmax = sqrt(bignum);
 /* Scale matrix to allowable range, if necessary. */
 anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1]);
 iscale = 0;
 if (anrm > 0. && anrm < rmin) {
 iscale = 1;
 sigma = rmin / anrm;
 }
 else if (anrm > rmax) {
 iscale = 1;
 sigma = rmax / anrm;
 }
 if (iscale == 1) {
 dlascl_(uplo, &c__0, &c__0, &c_b27, &sigma, n, n, &a[a_offset], lda, info);
 }
 /* Call DSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */
 inde = 1;
 indtau = inde + *n;
 indhous = indtau + *n;
 indwrk = indhous + lhtrd;
 llwork = *lwork - indwrk + 1;
 dsytrd_2stage_(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[inde], & work[indtau], &work[indhous], &lhtrd, &work[indwrk], &llwork, & iinfo);
 /* For eigenvalues only, call DSTERF. For eigenvectors, first call */
 /* DORGTR to generate the orthogonal matrix, then call DSTEQR. */
 if (! wantz) {
 dsterf_(n, &w[1], &work[inde], info);
 }
 else {
 /* Not available in this release, and argument checking should not */
 /* let it getting here */
 return 0;
 dorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], & llwork, &iinfo);
 dsteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau], info);
 }
 /* If matrix was scaled, then rescale eigenvalues appropriately. */
 if (iscale == 1) {
 if (*info == 0) {
 imax = *n;
 }
 else {
 imax = *info - 1;
 }
 d__1 = 1. / sigma;
 dscal_(&imax, &d__1, &w[1], &c__1);
 }
 /* Set WORK(1) to optimal workspace size. */
 work[1] = (doublereal) lwmin;
 return 0;
 /* End of DSYEV_2STAGE */
 }
 /* dsyev_2stage__ */
 