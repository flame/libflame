/* ../netlib/v3.9.0/zheevd_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 static integer c_n1 = -1;
 static integer c__2 = 2;
 static integer c__3 = 3;
 static integer c__4 = 4;
 static integer c__0 = 0;
 static doublereal c_b28 = 1.;
 /* > \brief <b> ZHEEVD_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b> */
 /* @precisions fortran z -> s d c */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download ZHEEVD_2STAGE + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevd_ 2stage.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevd_ 2stage.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevd_ 2stage.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE ZHEEVD_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, */
 /* RWORK, LRWORK, IWORK, LIWORK, INFO ) */
 /* IMPLICIT NONE */
 /* .. Scalar Arguments .. */
 /* CHARACTER JOBZ, UPLO */
 /* INTEGER INFO, LDA, LIWORK, LRWORK, LWORK, N */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IWORK( * ) */
 /* DOUBLE PRECISION RWORK( * ), W( * ) */
 /* COMPLEX*16 A( LDA, * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > ZHEEVD_2STAGE computes all eigenvalues and, optionally, eigenvectors of a */
 /* > complex Hermitian matrix A using the 2stage technique for */
 /* > the reduction to tridiagonal. If eigenvectors are desired, it uses a */
 /* > divide and conquer algorithm. */
 /* > */
 /* > The divide and conquer algorithm makes very mild assumptions about */
 /* > floating point arithmetic. It will work on machines with a guard */
 /* > digit in add/subtract, or on those binary machines without guard */
 /* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
 /* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
 /* > without guard digits, but we know of none. */
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
 /* > A is COMPLEX*16 array, dimension (LDA, N) */
 /* > On entry, the Hermitian matrix A. If UPLO = 'U', the */
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
 /* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
 /* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The dimension of the array WORK. */
 /* > If N <= 1, LWORK must be at least 1. */
 /* > If JOBZ = 'N' and N > 1, LWORK must be queried. */
 /* > LWORK = MAX(1, dimension) where */
 /* > dimension = max(stage1,stage2) + (KD+1)*N + N+1 */
 /* > = N*KD + N*max(KD+1,FACTOPTNB) */
 /* > + max(2*KD*KD, KD*NTHREADS) */
 /* > + (KD+1)*N + N+1 */
 /* > where KD is the blocking size of the reduction, */
 /* > FACTOPTNB is the blocking used by the QR or LQ */
 /* > algorithm, usually FACTOPTNB=128 is a good choice */
 /* > NTHREADS is the number of threads used when */
 /* > openMP compilation is enabled, otherwise =1. */
 /* > If JOBZ = 'V' and N > 1, LWORK must be at least 2*N + N**2 */
 /* > */
 /* > If LWORK = -1, then a workspace query is assumed;
 the routine */
 /* > only calculates the optimal sizes of the WORK, RWORK and */
 /* > IWORK arrays, returns these values as the first entries of */
 /* > the WORK, RWORK and IWORK arrays, and no error message */
 /* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] RWORK */
 /* > \verbatim */
 /* > RWORK is DOUBLE PRECISION array, */
 /* > dimension (LRWORK) */
 /* > On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LRWORK */
 /* > \verbatim */
 /* > LRWORK is INTEGER */
 /* > The dimension of the array RWORK. */
 /* > If N <= 1, LRWORK must be at least 1. */
 /* > If JOBZ = 'N' and N > 1, LRWORK must be at least N. */
 /* > If JOBZ = 'V' and N > 1, LRWORK must be at least */
 /* > 1 + 5*N + 2*N**2. */
 /* > */
 /* > If LRWORK = -1, then a workspace query is assumed;
 the */
 /* > routine only calculates the optimal sizes of the WORK, RWORK */
 /* > and IWORK arrays, returns these values as the first entries */
 /* > of the WORK, RWORK and IWORK arrays, and no error message */
 /* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IWORK */
 /* > \verbatim */
 /* > IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
 /* > On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LIWORK */
 /* > \verbatim */
 /* > LIWORK is INTEGER */
 /* > The dimension of the array IWORK. */
 /* > If N <= 1, LIWORK must be at least 1. */
 /* > If JOBZ = 'N' and N > 1, LIWORK must be at least 1. */
 /* > If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N. */
 /* > */
 /* > If LIWORK = -1, then a workspace query is assumed;
 the */
 /* > routine only calculates the optimal sizes of the WORK, RWORK */
 /* > and IWORK arrays, returns these values as the first entries */
 /* > of the WORK, RWORK and IWORK arrays, and no error message */
 /* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value */
 /* > > 0: if INFO = i and JOBZ = 'N', then the algorithm failed */
 /* > to converge;
 i off-diagonal elements of an intermediate */
 /* > tridiagonal form did not converge to zero;
 */
 /* > if INFO = i and JOBZ = 'V', then the algorithm failed */
 /* > to compute an eigenvalue while working on the submatrix */
 /* > lying in rows and columns INFO/(N+1) through */
 /* > mod(INFO,N+1). */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2017 */
 /* > \ingroup complex16HEeigen */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > Modified description of INFO. Sven, 16 Feb 05. */
 /* > \par Contributors: */
 /* ================== */
 /* > */
 /* > Jeff Rutter, Computer Science Division, University of California */
 /* > at Berkeley, USA */
 /* > */
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
 int zheevd_2stage_(char *jobz, char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info) {
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
 doublereal anrm;
 integer imax;
 doublereal rmin, rmax;
 extern /* Subroutine */
 int dscal_(integer *, doublereal *, doublereal *, integer *);
 doublereal sigma;
 extern logical lsame_(char *, char *);
 extern /* Subroutine */
 int zhetrd_2stage_(char *, char *, integer *, doublecomplex *, integer *, doublereal *, doublereal *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
 integer iinfo, lhtrd, lwmin;
 logical lower;
 integer llrwk, lwtrd;
 logical wantz;
 integer indwk2, llwrk2;
 extern doublereal dlamch_(char *);
 integer iscale;
 doublereal safmin;
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 doublereal bignum;
 extern doublereal zlanhe_(char *, char *, integer *, doublecomplex *, integer *, doublereal *);
 integer indtau;
 extern /* Subroutine */
 int dsterf_(integer *, doublereal *, doublereal *, integer *), zlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublecomplex *, integer *, integer *), zstedc_(char *, integer *, doublereal *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, integer *, integer *, integer *, integer *);
 integer indrwk, indwrk, liwmin;
 extern /* Subroutine */
 int zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
 integer lrwmin, llwork;
 doublereal smlnum;
 logical lquery;
 extern /* Subroutine */
 int zunmtr_(char *, char *, char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
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
 --rwork;
 --iwork;
 /* Function Body */
 wantz = lsame_(jobz, "V");
 lower = lsame_(uplo, "L");
 lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;
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
 if (*n <= 1) {
 lwmin = 1;
 lrwmin = 1;
 liwmin = 1;
 }
 else {
 kd = ilaenv2stage_(&c__1, "ZHETRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1);
 ib = ilaenv2stage_(&c__2, "ZHETRD_2STAGE", jobz, n, &kd, &c_n1, & c_n1);
 lhtrd = ilaenv2stage_(&c__3, "ZHETRD_2STAGE", jobz, n, &kd, &ib, & c_n1);
 lwtrd = ilaenv2stage_(&c__4, "ZHETRD_2STAGE", jobz, n, &kd, &ib, & c_n1);
 if (wantz) {
 lwmin = (*n << 1) + *n * *n;
 /* Computing 2nd power */
 i__1 = *n;
 lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
 liwmin = *n * 5 + 3;
 }
 else {
 lwmin = *n + 1 + lhtrd + lwtrd;
 lrwmin = *n;
 liwmin = 1;
 }
 }
 work[1].r = (doublereal) lwmin; work[1].i = 0.; // , expr subst  
 rwork[1] = (doublereal) lrwmin;
 iwork[1] = liwmin;
 if (*lwork < lwmin && ! lquery) {
 *info = -8;
 }
 else if (*lrwork < lrwmin && ! lquery) {
 *info = -10;
 }
 else if (*liwork < liwmin && ! lquery) {
 *info = -12;
 }
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("ZHEEVD_2STAGE", &i__1);
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
 i__1 = a_dim1 + 1;
 w[1] = a[i__1].r;
 if (wantz) {
 i__1 = a_dim1 + 1;
 a[i__1].r = 1.; a[i__1].i = 0.; // , expr subst  
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
 anrm = zlanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1]);
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
 zlascl_(uplo, &c__0, &c__0, &c_b28, &sigma, n, n, &a[a_offset], lda, info);
 }
 /* Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form. */
 inde = 1;
 indrwk = inde + *n;
 llrwk = *lrwork - indrwk + 1;
 indtau = 1;
 indhous = indtau + *n;
 indwrk = indhous + lhtrd;
 llwork = *lwork - indwrk + 1;
 indwk2 = indwrk + *n * *n;
 llwrk2 = *lwork - indwk2 + 1;
 zhetrd_2stage_(jobz, uplo, n, &a[a_offset], lda, &w[1], &rwork[inde], & work[indtau], &work[indhous], &lhtrd, &work[indwrk], &llwork, & iinfo);
 /* For eigenvalues only, call DSTERF. For eigenvectors, first call */
 /* ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
 /* tridiagonal matrix, then call ZUNMTR to multiply it to the */
 /* Householder transformations represented as Householder vectors in */
 /* A. */
 if (! wantz) {
 dsterf_(n, &w[1], &rwork[inde], info);
 }
 else {
 zstedc_("I", n, &w[1], &rwork[inde], &work[indwrk], n, &work[indwk2], &llwrk2, &rwork[indrwk], &llrwk, &iwork[1], liwork, info);
 zunmtr_("L", uplo, "N", n, n, &a[a_offset], lda, &work[indtau], &work[ indwrk], n, &work[indwk2], &llwrk2, &iinfo);
 zlacpy_("A", n, n, &work[indwrk], n, &a[a_offset], lda);
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
 work[1].r = (doublereal) lwmin; work[1].i = 0.; // , expr subst  
 rwork[1] = (doublereal) lrwmin;
 iwork[1] = liwmin;
 return 0;
 /* End of ZHEEVD_2STAGE */
 }
 /* zheevd_2stage__ */
 