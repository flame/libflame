/* ../netlib/v3.9.0/ssbev_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__2 = 2;
 static integer c_n1 = -1;
 static integer c__3 = 3;
 static integer c__4 = 4;
 static real c_b21 = 1.f;
 static integer c__1 = 1;
 /* > \brief <b> SSBEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b> */
 /* @generated from dsbev_2stage.f, fortran d -> s, Sat Nov 5 23:58:09 2016 */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download SSBEV_2STAGE + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbev_2 stage.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbev_2 stage.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbev_2 stage.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE SSBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, */
 /* WORK, LWORK, INFO ) */
 /* IMPLICIT NONE */
 /* .. Scalar Arguments .. */
 /* CHARACTER JOBZ, UPLO */
 /* INTEGER INFO, KD, LDAB, LDZ, N, LWORK */
 /* .. */
 /* .. Array Arguments .. */
 /* REAL AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > SSBEV_2STAGE computes all the eigenvalues and, optionally, eigenvectors of */
 /* > a real symmetric band matrix A using the 2stage technique for */
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
 /* > \param[in] KD */
 /* > \verbatim */
 /* > KD is INTEGER */
 /* > The number of superdiagonals of the matrix A if UPLO = 'U', */
 /* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] AB */
 /* > \verbatim */
 /* > AB is REAL array, dimension (LDAB, N) */
 /* > On entry, the upper or lower triangle of the symmetric band */
 /* > matrix A, stored in the first KD+1 rows of the array. The */
 /* > j-th column of A is stored in the j-th column of the array AB */
 /* > as follows: */
 /* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
 */
 /* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd). */
 /* > */
 /* > On exit, AB is overwritten by values generated during the */
 /* > reduction to tridiagonal form. If UPLO = 'U', the first */
 /* > superdiagonal and the diagonal of the tridiagonal matrix T */
 /* > are returned in rows KD and KD+1 of AB, and if UPLO = 'L', */
 /* > the diagonal and first subdiagonal of T are returned in the */
 /* > first two rows of AB. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDAB */
 /* > \verbatim */
 /* > LDAB is INTEGER */
 /* > The leading dimension of the array AB. LDAB >= KD + 1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] W */
 /* > \verbatim */
 /* > W is REAL array, dimension (N) */
 /* > If INFO = 0, the eigenvalues in ascending order. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] Z */
 /* > \verbatim */
 /* > Z is REAL array, dimension (LDZ, N) */
 /* > If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
 /* > eigenvectors of the matrix A, with the i-th column of Z */
 /* > holding the eigenvector associated with W(i). */
 /* > If JOBZ = 'N', then Z is not referenced. */
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
 /* > WORK is REAL array, dimension LWORK */
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
 /* > dimension = (2KD+1)*N + KD*NTHREADS + N */
 /* > where KD is the size of the band. */
 /* > NTHREADS is the number of threads used when */
 /* > openMP compilation is enabled, otherwise =1. */
 /* > If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available. */
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
 /* > \ingroup realOTHEReigen */
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
 int ssbev_2stage_(char *jobz, char *uplo, integer *n, integer *kd, real *ab, integer *ldab, real *w, real *z__, integer * ldz, real *work, integer *lwork, integer *info) {
 /* System generated locals */
 integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
 real r__1;
 /* Builtin functions */
 double sqrt(doublereal);
 /* Local variables */
 integer ib;
 real eps;
 integer inde;
 extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
 real anrm;
 integer imax;
 real rmin, rmax;
 extern /* Subroutine */
 int ssytrd_sb2st_(char *, char *, char *, integer *, integer *, real *, integer *, real *, real *, real *, integer *, real *, integer *, integer *);
 real sigma;
 extern logical lsame_(char *, char *);
 integer iinfo;
 extern /* Subroutine */
 int sscal_(integer *, real *, real *, integer *);
 integer lhtrd, lwmin;
 logical lower;
 integer lwtrd;
 logical wantz;
 integer iscale;
 extern real slamch_(char *);
 real safmin;
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 real bignum;
 extern real slansb_(char *, char *, integer *, integer *, real *, integer *, real *);
 extern /* Subroutine */
 int slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *, integer *, integer *);
 integer indwrk;
 extern /* Subroutine */
 int ssterf_(integer *, real *, real *, integer *);
 integer llwork;
 real smlnum;
 logical lquery;
 extern /* Subroutine */
 int ssteqr_(char *, integer *, real *, real *, real *, integer *, real *, integer *);
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
 ab_dim1 = *ldab;
 ab_offset = 1 + ab_dim1;
 ab -= ab_offset;
 --w;
 z_dim1 = *ldz;
 z_offset = 1 + z_dim1;
 z__ -= z_offset;
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
 else if (*kd < 0) {
 *info = -4;
 }
 else if (*ldab < *kd + 1) {
 *info = -6;
 }
 else if (*ldz < 1 || wantz && *ldz < *n) {
 *info = -9;
 }
 if (*info == 0) {
 if (*n <= 1) {
 lwmin = 1;
 work[1] = (real) lwmin;
 }
 else {
 ib = ilaenv2stage_(&c__2, "SSYTRD_SB2ST", jobz, n, kd, &c_n1, & c_n1);
 lhtrd = ilaenv2stage_(&c__3, "SSYTRD_SB2ST", jobz, n, kd, &ib, & c_n1);
 lwtrd = ilaenv2stage_(&c__4, "SSYTRD_SB2ST", jobz, n, kd, &ib, & c_n1);
 lwmin = *n + lhtrd + lwtrd;
 work[1] = (real) lwmin;
 }
 if (*lwork < lwmin && ! lquery) {
 *info = -11;
 }
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("SSBEV_2STAGE ", &i__1);
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
 if (lower) {
 w[1] = ab[ab_dim1 + 1];
 }
 else {
 w[1] = ab[*kd + 1 + ab_dim1];
 }
 if (wantz) {
 z__[z_dim1 + 1] = 1.f;
 }
 return 0;
 }
 /* Get machine constants. */
 safmin = slamch_("Safe minimum");
 eps = slamch_("Precision");
 smlnum = safmin / eps;
 bignum = 1.f / smlnum;
 rmin = sqrt(smlnum);
 rmax = sqrt(bignum);
 /* Scale matrix to allowable range, if necessary. */
 anrm = slansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1]);
 iscale = 0;
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
 slascl_("B", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, info);
 }
 else {
 slascl_("Q", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, info);
 }
 }
 /* Call SSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form. */
 inde = 1;
 indhous = inde + *n;
 indwrk = indhous + lhtrd;
 llwork = *lwork - indwrk + 1;
 ssytrd_sb2st_("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], &work[ inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, &iinfo);
 /* For eigenvalues only, call SSTERF. For eigenvectors, call SSTEQR. */
 if (! wantz) {
 ssterf_(n, &w[1], &work[inde], info);
 }
 else {
 ssteqr_(jobz, n, &w[1], &work[inde], &z__[z_offset], ldz, &work[ indwrk], info);
 }
 /* If matrix was scaled, then rescale eigenvalues appropriately. */
 if (iscale == 1) {
 if (*info == 0) {
 imax = *n;
 }
 else {
 imax = *info - 1;
 }
 r__1 = 1.f / sigma;
 sscal_(&imax, &r__1, &w[1], &c__1);
 }
 /* Set WORK(1) to optimal workspace size. */
 work[1] = (real) lwmin;
 return 0;
 /* End of SSBEV_2STAGE */
 }
 /* ssbev_2stage__ */
 