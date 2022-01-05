/* ../netlib/v3.9.0/dsbevx_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__2 = 2;
 static integer c_n1 = -1;
 static integer c__3 = 3;
 static integer c__4 = 4;
 static doublereal c_b24 = 1.;
 static integer c__1 = 1;
 static doublereal c_b45 = 0.;
 /* > \brief <b> DSBEVX_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b> */
 /* @precisions fortran d -> s */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DSBEVX_2STAGE + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsbevx_ 2stage.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsbevx_ 2stage.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsbevx_ 2stage.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE DSBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, */
 /* LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, */
 /* LDZ, WORK, LWORK, IWORK, IFAIL, INFO ) */
 /* IMPLICIT NONE */
 /* .. Scalar Arguments .. */
 /* CHARACTER JOBZ, RANGE, UPLO */
 /* INTEGER IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK */
 /* DOUBLE PRECISION ABSTOL, VL, VU */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IFAIL( * ), IWORK( * ) */
 /* DOUBLE PRECISION AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ), */
 /* $ Z( LDZ, * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DSBEVX_2STAGE computes selected eigenvalues and, optionally, eigenvectors */
 /* > of a real symmetric band matrix A using the 2stage technique for */
 /* > the reduction to tridiagonal. Eigenvalues and eigenvectors can */
 /* > be selected by specifying either a range of values or a range of */
 /* > indices for the desired eigenvalues. */
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
 /* > = 'A': all eigenvalues will be found;
 */
 /* > = 'V': all eigenvalues in the half-open interval (VL,VU] */
 /* > will be found;
 */
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
 /* > \param[in] KD */
 /* > \verbatim */
 /* > KD is INTEGER */
 /* > The number of superdiagonals of the matrix A if UPLO = 'U', */
 /* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] AB */
 /* > \verbatim */
 /* > AB is DOUBLE PRECISION array, dimension (LDAB, N) */
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
 /* > \param[out] Q */
 /* > \verbatim */
 /* > Q is DOUBLE PRECISION array, dimension (LDQ, N) */
 /* > If JOBZ = 'V', the N-by-N orthogonal matrix used in the */
 /* > reduction to tridiagonal form. */
 /* > If JOBZ = 'N', the array Q is not referenced. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDQ */
 /* > \verbatim */
 /* > LDQ is INTEGER */
 /* > The leading dimension of the array Q. If JOBZ = 'V', then */
 /* > LDQ >= max(1,N). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] VL */
 /* > \verbatim */
 /* > VL is DOUBLE PRECISION */
 /* > If RANGE='V', the lower bound of the interval to */
 /* > be searched for eigenvalues. VL < VU. */
 /* > Not referenced if RANGE = 'A' or 'I'. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] VU */
 /* > \verbatim */
 /* > VU is DOUBLE PRECISION */
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
 /* > ABSTOL is DOUBLE PRECISION */
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
 /* > by reducing AB to tridiagonal form. */
 /* > */
 /* > Eigenvalues will be computed most accurately when ABSTOL is */
 /* > set to twice the underflow threshold 2*DLAMCH('S'), not zero. */
 /* > If this routine returns with INFO>0, indicating that some */
 /* > eigenvectors did not converge, try setting ABSTOL to */
 /* > 2*DLAMCH('S'). */
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
 /* > W is DOUBLE PRECISION array, dimension (N) */
 /* > The first M elements contain the selected eigenvalues in */
 /* > ascending order. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] Z */
 /* > \verbatim */
 /* > Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M)) */
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
 /* > WORK is DOUBLE PRECISION array, dimension (LWORK) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The length of the array WORK. LWORK >= 1, when N <= 1;
 */
 /* > otherwise */
 /* > If JOBZ = 'N' and N > 1, LWORK must be queried. */
 /* > LWORK = MAX(1, 7*N, dimension) where */
 /* > dimension = (2KD+1)*N + KD*NTHREADS + 2*N */
 /* > where KD is the size of the band. */
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
 /* > = 0: successful exit. */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value. */
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
 /* > \ingroup doubleOTHEReigen */
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
 int dsbevx_2stage_(char *jobz, char *range, char *uplo, integer *n, integer *kd, doublereal *ab, integer *ldab, doublereal *q, integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer * iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, integer *lwork, integer *iwork, integer *ifail, integer *info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
 char buffer[256]; 
 snprintf(buffer, 256,"dsbevx_2stage inputs: jobz %c, range %c, uplo %c, n %" FLA_IS ", kd %" FLA_IS ", ldab %" FLA_IS ", ldq %" FLA_IS ", il %" FLA_IS ", iu %" FLA_IS ", ldz %" FLA_IS ", lwork %" FLA_IS "",*jobz, *range, *uplo, *n, *kd, *ldab, *ldq, *il, *iu, *ldz, *lwork);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer ab_dim1, ab_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2;
 doublereal d__1, d__2;
 /* Builtin functions */
 double sqrt(doublereal);
 /* Local variables */
 integer i__, j, ib, jj;
 doublereal eps, vll, vuu, tmp1;
 integer indd, inde;
 extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
 doublereal anrm;
 integer imax;
 doublereal rmin, rmax;
 logical test;
 extern /* Subroutine */
 int dsytrd_sb2st_(char *, char *, char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
 integer itmp1, indee;
 extern /* Subroutine */
 int dscal_(integer *, doublereal *, doublereal *, integer *);
 doublereal sigma;
 extern logical lsame_(char *, char *);
 extern /* Subroutine */
 int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
 integer iinfo;
 char order[1];
 integer lhtrd;
 extern /* Subroutine */
 int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *), dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
 integer lwmin;
 logical lower;
 integer lwtrd;
 logical wantz;
 extern doublereal dlamch_(char *);
 logical alleig, indeig;
 integer iscale, indibl;
 extern /* Subroutine */
 int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
 extern doublereal dlansb_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *);
 logical valeig;
 extern /* Subroutine */
 int dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *);
 doublereal safmin;
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 doublereal abstll, bignum;
 integer indisp;
 extern /* Subroutine */
 int dstein_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, integer *), dsterf_(integer *, doublereal *, doublereal *, integer *);
 integer indiwo;
 extern /* Subroutine */
 int dstebz_(char *, char *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
 integer indwrk;
 extern /* Subroutine */
 int dsteqr_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *);
 integer nsplit, llwork;
 doublereal smlnum;
 logical lquery;
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
 ab_dim1 = *ldab;
 ab_offset = 1 + ab_dim1;
 ab -= ab_offset;
 q_dim1 = *ldq;
 q_offset = 1 + q_dim1;
 q -= q_offset;
 --w;
 z_dim1 = *ldz;
 z_offset = 1 + z_dim1;
 z__ -= z_offset;
 --work;
 --iwork;
 --ifail;
 /* Function Body */
 wantz = lsame_(jobz, "V");
 alleig = lsame_(range, "A");
 valeig = lsame_(range, "V");
 indeig = lsame_(range, "I");
 lower = lsame_(uplo, "L");
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
 else if (*kd < 0) {
 *info = -5;
 }
 else if (*ldab < *kd + 1) {
 *info = -7;
 }
 else if (wantz && *ldq < max(1,*n)) {
 *info = -9;
 }
 else {
 if (valeig) {
 if (*n > 0 && *vu <= *vl) {
 *info = -11;
 }
 }
 else if (indeig) {
 if (*il < 1 || *il > max(1,*n)) {
 *info = -12;
 }
 else if (*iu < min(*n,*il) || *iu > *n) {
 *info = -13;
 }
 }
 }
 if (*info == 0) {
 if (*ldz < 1 || wantz && *ldz < *n) {
 *info = -18;
 }
 }
 if (*info == 0) {
 if (*n <= 1) {
 lwmin = 1;
 work[1] = (doublereal) lwmin;
 }
 else {
 ib = ilaenv2stage_(&c__2, "DSYTRD_SB2ST", jobz, n, kd, &c_n1, & c_n1);
 lhtrd = ilaenv2stage_(&c__3, "DSYTRD_SB2ST", jobz, n, kd, &ib, & c_n1);
 lwtrd = ilaenv2stage_(&c__4, "DSYTRD_SB2ST", jobz, n, kd, &ib, & c_n1);
 lwmin = (*n << 1) + lhtrd + lwtrd;
 work[1] = (doublereal) lwmin;
 }
 if (*lwork < lwmin && ! lquery) {
 *info = -20;
 }
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("DSBEVX_2STAGE ", &i__1);
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
 *m = 1;
 if (lower) {
 tmp1 = ab[ab_dim1 + 1];
 }
 else {
 tmp1 = ab[*kd + 1 + ab_dim1];
 }
 if (valeig) {
 if (! (*vl < tmp1 && *vu >= tmp1)) {
 *m = 0;
 }
 }
 if (*m == 1) {
 w[1] = tmp1;
 if (wantz) {
 z__[z_dim1 + 1] = 1.;
 }
 }
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Get machine constants. */
 safmin = dlamch_("Safe minimum");
 eps = dlamch_("Precision");
 smlnum = safmin / eps;
 bignum = 1. / smlnum;
 rmin = sqrt(smlnum);
 /* Computing MIN */
 d__1 = sqrt(bignum); d__2 = 1. / sqrt(sqrt(safmin)); // , expr subst  
 rmax = min(d__1,d__2);
 /* Scale matrix to allowable range, if necessary. */
 iscale = 0;
 abstll = *abstol;
 if (valeig) {
 vll = *vl;
 vuu = *vu;
 }
 else {
 vll = 0.;
 vuu = 0.;
 }
 anrm = dlansb_("M", uplo, n, kd, &ab[ab_offset], ldab, &work[1]);
 if (anrm > 0. && anrm < rmin) {
 iscale = 1;
 sigma = rmin / anrm;
 }
 else if (anrm > rmax) {
 iscale = 1;
 sigma = rmax / anrm;
 }
 if (iscale == 1) {
 if (lower) {
 dlascl_("B", kd, kd, &c_b24, &sigma, n, n, &ab[ab_offset], ldab, info);
 }
 else {
 dlascl_("Q", kd, kd, &c_b24, &sigma, n, n, &ab[ab_offset], ldab, info);
 }
 if (*abstol > 0.) {
 abstll = *abstol * sigma;
 }
 if (valeig) {
 vll = *vl * sigma;
 vuu = *vu * sigma;
 }
 }
 /* Call DSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form. */
 indd = 1;
 inde = indd + *n;
 indhous = inde + *n;
 indwrk = indhous + lhtrd;
 llwork = *lwork - indwrk + 1;
 dsytrd_sb2st_("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &work[indd], &work[inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, & iinfo);
 /* If all eigenvalues are desired and ABSTOL is less than or equal */
 /* to zero, then call DSTERF or SSTEQR. If this fails for some */
 /* eigenvalue, then try DSTEBZ. */
 test = FALSE_;
 if (indeig) {
 if (*il == 1 && *iu == *n) {
 test = TRUE_;
 }
 }
 if ((alleig || test) && *abstol <= 0.) {
 dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
 indee = indwrk + (*n << 1);
 if (! wantz) {
 i__1 = *n - 1;
 dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
 dsterf_(n, &w[1], &work[indee], info);
 }
 else {
 dlacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz);
 i__1 = *n - 1;
 dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
 dsteqr_(jobz, n, &w[1], &work[indee], &z__[z_offset], ldz, &work[ indwrk], info);
 if (*info == 0) {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 ifail[i__] = 0;
 /* L10: */
 }
 }
 }
 if (*info == 0) {
 *m = *n;
 goto L30;
 }
 *info = 0;
 }
 /* Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */
 if (wantz) {
 *(unsigned char *)order = 'B';
 }
 else {
 *(unsigned char *)order = 'E';
 }
 indibl = 1;
 indisp = indibl + *n;
 indiwo = indisp + *n;
 dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[ inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[ indwrk], &iwork[indiwo], info);
 if (wantz) {
 dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[ indisp], &z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], & ifail[1], info);
 /* Apply orthogonal matrix used in reduction to tridiagonal */
 /* form to eigenvectors returned by DSTEIN. */
 i__1 = *m;
 for (j = 1;
 j <= i__1;
 ++j) {
 dcopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
 dgemv_("N", n, n, &c_b24, &q[q_offset], ldq, &work[1], &c__1, & c_b45, &z__[j * z_dim1 + 1], &c__1);
 /* L20: */
 }
 }
 /* If matrix was scaled, then rescale eigenvalues appropriately. */
 L30: if (iscale == 1) {
 if (*info == 0) {
 imax = *m;
 }
 else {
 imax = *info - 1;
 }
 d__1 = 1. / sigma;
 dscal_(&imax, &d__1, &w[1], &c__1);
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
 /* L40: */
 }
 if (i__ != 0) {
 itmp1 = iwork[indibl + i__ - 1];
 w[i__] = w[j];
 iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
 w[j] = tmp1;
 iwork[indibl + j - 1] = itmp1;
 dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
 if (*info != 0) {
 itmp1 = ifail[i__];
 ifail[i__] = ifail[j];
 ifail[j] = itmp1;
 }
 }
 /* L50: */
 }
 }
 /* Set WORK(1) to optimal workspace size. */
 work[1] = (doublereal) lwmin;
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of DSBEVX_2STAGE */
 }
 /* dsbevx_2stage__ */
 