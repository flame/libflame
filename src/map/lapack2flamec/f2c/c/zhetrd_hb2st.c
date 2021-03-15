/* ../netlib/v3.9.0/zhetrd_hb2st.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static doublecomplex c_b1 = {
0.,0.}
;
 static integer c__2 = 2;
 static integer c_n1 = -1;
 static integer c__3 = 3;
 static integer c__4 = 4;
 /* > \brief \b ZHETRD_HB2ST reduces a complex Hermitian band matrix A to real symmetric tridiagonal form T */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download ZHETRD_HB2ST + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbtrd_ hb2st.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbtrd_ hb2st.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbtrd_ hb2st.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE ZHETRD_HB2ST( STAGE1, VECT, UPLO, N, KD, AB, LDAB, */
 /* D, E, HOUS, LHOUS, WORK, LWORK, INFO ) */
 /* #if defined(_OPENMP) */
 /* use omp_lib */
 /* #endif */
 /* IMPLICIT NONE */
 /* .. Scalar Arguments .. */
 /* CHARACTER STAGE1, UPLO, VECT */
 /* INTEGER N, KD, IB, LDAB, LHOUS, LWORK, INFO */
 /* .. */
 /* .. Array Arguments .. */
 /* DOUBLE PRECISION D( * ), E( * ) */
 /* COMPLEX*16 AB( LDAB, * ), HOUS( * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > ZHETRD_HB2ST reduces a complex Hermitian band matrix A to real symmetric */
 /* > tridiagonal form T by a unitary similarity transformation: */
 /* > Q**H * A * Q = T. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] STAGE1 */
 /* > \verbatim */
 /* > STAGE1 is CHARACTER*1 */
 /* > = 'N': "No": to mention that the stage 1 of the reduction */
 /* > from dense to band using the zhetrd_he2hb routine */
 /* > was not called before this routine to reproduce AB. */
 /* > In other term this routine is called as standalone. */
 /* > = 'Y': "Yes": to mention that the stage 1 of the */
 /* > reduction from dense to band using the zhetrd_he2hb */
 /* > routine has been called to produce AB (e.g., AB is */
 /* > the output of zhetrd_he2hb. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] VECT */
 /* > \verbatim */
 /* > VECT is CHARACTER*1 */
 /* > = 'N': No need for the Housholder representation, */
 /* > and thus LHOUS is of size max(1, 4*N);
 */
 /* > = 'V': the Householder representation is needed to */
 /* > either generate or to apply Q later on, */
 /* > then LHOUS is to be queried and computed. */
 /* > (NOT AVAILABLE IN THIS RELEASE). */
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
 /* > AB is COMPLEX*16 array, dimension (LDAB,N) */
 /* > On entry, the upper or lower triangle of the Hermitian band */
 /* > matrix A, stored in the first KD+1 rows of the array. The */
 /* > j-th column of A is stored in the j-th column of the array AB */
 /* > as follows: */
 /* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
 */
 /* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd). */
 /* > On exit, the diagonal elements of AB are overwritten by the */
 /* > diagonal elements of the tridiagonal matrix T;
 if KD > 0, the */
 /* > elements on the first superdiagonal (if UPLO = 'U') or the */
 /* > first subdiagonal (if UPLO = 'L') are overwritten by the */
 /* > off-diagonal elements of T;
 the rest of AB is overwritten by */
 /* > values generated during the reduction. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDAB */
 /* > \verbatim */
 /* > LDAB is INTEGER */
 /* > The leading dimension of the array AB. LDAB >= KD+1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] D */
 /* > \verbatim */
 /* > D is DOUBLE PRECISION array, dimension (N) */
 /* > The diagonal elements of the tridiagonal matrix T. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] E */
 /* > \verbatim */
 /* > E is DOUBLE PRECISION array, dimension (N-1) */
 /* > The off-diagonal elements of the tridiagonal matrix T: */
 /* > E(i) = T(i,i+1) if UPLO = 'U';
 E(i) = T(i+1,i) if UPLO = 'L'. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] HOUS */
 /* > \verbatim */
 /* > HOUS is COMPLEX*16 array, dimension LHOUS, that */
 /* > store the Householder representation. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LHOUS */
 /* > \verbatim */
 /* > LHOUS is INTEGER */
 /* > The dimension of the array HOUS. LHOUS = MAX(1, dimension) */
 /* > If LWORK = -1, or LHOUS=-1, */
 /* > then a query is assumed;
 the routine */
 /* > only calculates the optimal size of the HOUS array, returns */
 /* > this value as the first entry of the HOUS array, and no error */
 /* > message related to LHOUS is issued by XERBLA. */
 /* > LHOUS = MAX(1, dimension) where */
 /* > dimension = 4*N if VECT='N' */
 /* > not available now if VECT='H' */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is COMPLEX*16 array, dimension LWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The dimension of the array WORK. LWORK = MAX(1, dimension) */
 /* > If LWORK = -1, or LHOUS=-1, */
 /* > then a workspace query is assumed;
 the routine */
 /* > only calculates the optimal size of the WORK array, returns */
 /* > this value as the first entry of the WORK array, and no error */
 /* > message related to LWORK is issued by XERBLA. */
 /* > LWORK = MAX(1, dimension) where */
 /* > dimension = (2KD+1)*N + KD*NTHREADS */
 /* > where KD is the blocking size of the reduction, */
 /* > FACTOPTNB is the blocking used by the QR or LQ */
 /* > algorithm, usually FACTOPTNB=128 is a good choice */
 /* > NTHREADS is the number of threads used when */
 /* > openMP compilation is enabled, otherwise =1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2017 */
 /* > \ingroup complex16OTHERcomputational */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > Implemented by Azzam Haidar. */
 /* > */
 /* > All details are available on technical report, SC11, SC13 papers. */
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
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int zhetrd_hb2st_(char *stage1, char *vect, char *uplo, integer *n, integer *kd, doublecomplex *ab, integer *ldab, doublereal *d__, doublereal *e, doublecomplex *hous, integer *lhous, doublecomplex *work, integer *lwork, integer *info) {
 /* System generated locals */
 integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
 doublecomplex z__1;
 /* Builtin functions */
 double z_abs(doublecomplex *);
 /* Local variables */
 integer abofdpos, nthreads, i__, k, m, stepercol, ed, ib, st, blklastind, lda, tid, ldv;
 doublecomplex tmp;
 integer stt, inda;
 extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
 integer thed, myid, indw, apos, dpos, indv, edind, debug;
 extern logical lsame_(char *, char *);
 integer lhmin, sizea, shift, stind, colpt, lwmin, awpos;
 logical wantq, upper;
 integer grsiz, sizev, ttype, abdpos;
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 integer thgrid, thgrnb, indtau;
 doublereal abstmp;
 integer ofdpos;
 extern /* Subroutine */
 int zhb2st_kernels_(char *, logical *, integer *, integer *, integer *, integer *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *), zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *), zlaset_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, integer *);
 logical lquery, afters1;
 extern /* Subroutine */
 int f90_exit_(void);
 integer ceiltmp, sweepid, nbtiles, sizetau, thgrsiz;
 /* #if defined(_OPENMP) */
 /* use omp_lib */
 /* #endif */
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
 /* .. */
 /* .. Local Scalars .. */
 /* .. */
 /* .. External Subroutines .. */
 /* .. */
 /* .. Intrinsic Functions .. */
 /* .. */
 /* .. External Functions .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Determine the minimal workspace size required. */
 /* Test the input parameters */
 /* Parameter adjustments */
 ab_dim1 = *ldab;
 ab_offset = 1 + ab_dim1;
 ab -= ab_offset;
 --d__;
 --e;
 --hous;
 --work;
 /* Function Body */
 debug = 0;
 *info = 0;
 afters1 = lsame_(stage1, "Y");
 wantq = lsame_(vect, "V");
 upper = lsame_(uplo, "U");
 lquery = *lwork == -1 || *lhous == -1;
 /* Determine the block size, the workspace size and the hous size. */
 ib = ilaenv2stage_(&c__2, "ZHETRD_HB2ST", vect, n, kd, &c_n1, &c_n1);
 lhmin = ilaenv2stage_(&c__3, "ZHETRD_HB2ST", vect, n, kd, &ib, &c_n1);
 lwmin = ilaenv2stage_(&c__4, "ZHETRD_HB2ST", vect, n, kd, &ib, &c_n1);
 if (! afters1 && ! lsame_(stage1, "N")) {
 *info = -1;
 }
 else if (! lsame_(vect, "N")) {
 *info = -2;
 }
 else if (! upper && ! lsame_(uplo, "L")) {
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
 else if (*lhous < lhmin && ! lquery) {
 *info = -11;
 }
 else if (*lwork < lwmin && ! lquery) {
 *info = -13;
 }
 if (*info == 0) {
 hous[1].r = (doublereal) lhmin; hous[1].i = 0.; // , expr subst  
 work[1].r = (doublereal) lwmin; work[1].i = 0.; // , expr subst  
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("ZHETRD_HB2ST", &i__1);
 return 0;
 }
 else if (lquery) {
 return 0;
 }
 /* Quick return if possible */
 if (*n == 0) {
 hous[1].r = 1.; hous[1].i = 0.; // , expr subst  
 work[1].r = 1.; work[1].i = 0.; // , expr subst  
 return 0;
 }
 /* Determine pointer position */
 ldv = *kd + ib;
 sizetau = *n << 1;
 sizev = *n << 1;
 indtau = 1;
 indv = indtau + sizetau;
 lda = (*kd << 1) + 1;
 sizea = lda * *n;
 inda = 1;
 indw = inda + sizea;
 nthreads = 1;
 tid = 0;
 if (upper) {
 apos = inda + *kd;
 awpos = inda;
 dpos = apos + *kd;
 ofdpos = dpos - 1;
 abdpos = *kd + 1;
 abofdpos = *kd;
 }
 else {
 apos = inda;
 awpos = inda + *kd + 1;
 dpos = apos;
 ofdpos = dpos + 1;
 abdpos = 1;
 abofdpos = 2;
 }
 /* Case KD=0: */
 /* The matrix is diagonal. We just copy it (convert to "real" for */
 /* complex because D is double and the imaginary part should be 0) */
 /* and store it in D. A sequential code here is better or */
 /* in a parallel environment it might need two cores for D and E */
 if (*kd == 0) {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = abdpos + i__ * ab_dim1;
 d__[i__] = ab[i__2].r;
 /* L30: */
 }
 i__1 = *n - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 e[i__] = 0.;
 /* L40: */
 }
 hous[1].r = 1.; hous[1].i = 0.; // , expr subst  
 work[1].r = 1.; work[1].i = 0.; // , expr subst  
 return 0;
 }
 /* Case KD=1: */
 /* The matrix is already Tridiagonal. We have to make diagonal */
 /* and offdiagonal elements real, and store them in D and E. */
 /* For that, for real precision just copy the diag and offdiag */
 /* to D and E while for the COMPLEX case the bulge chasing is */
 /* performed to convert the hermetian tridiagonal to symmetric */
 /* tridiagonal. A simpler coversion formula might be used, but then */
 /* updating the Q matrix will be required and based if Q is generated */
 /* or not this might complicate the story. */
 if (*kd == 1) {
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = abdpos + i__ * ab_dim1;
 d__[i__] = ab[i__2].r;
 /* L50: */
 }
 /* make off-diagonal elements real and copy them to E */
 if (upper) {
 i__1 = *n - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = abofdpos + (i__ + 1) * ab_dim1;
 tmp.r = ab[i__2].r; tmp.i = ab[i__2].i; // , expr subst  
 abstmp = z_abs(&tmp);
 i__2 = abofdpos + (i__ + 1) * ab_dim1;
 ab[i__2].r = abstmp; ab[i__2].i = 0.; // , expr subst  
 e[i__] = abstmp;
 if (abstmp != 0.) {
 z__1.r = tmp.r / abstmp; z__1.i = tmp.i / abstmp; // , expr subst  
 tmp.r = z__1.r; tmp.i = z__1.i; // , expr subst  
 }
 else {
 tmp.r = 1.; tmp.i = 0.; // , expr subst  
 }
 if (i__ < *n - 1) {
 i__2 = abofdpos + (i__ + 2) * ab_dim1;
 i__3 = abofdpos + (i__ + 2) * ab_dim1;
 z__1.r = ab[i__3].r * tmp.r - ab[i__3].i * tmp.i; z__1.i = ab[i__3].r * tmp.i + ab[i__3].i * tmp.r; // , expr subst  
 ab[i__2].r = z__1.r; ab[i__2].i = z__1.i; // , expr subst  
 }
 /* IF( WANTZ ) THEN */
 /* CALL ZSCAL( N, DCONJG( TMP ), Q( 1, I+1 ), 1 ) */
 /* END IF */
 /* L60: */
 }
 }
 else {
 i__1 = *n - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = abofdpos + i__ * ab_dim1;
 tmp.r = ab[i__2].r; tmp.i = ab[i__2].i; // , expr subst  
 abstmp = z_abs(&tmp);
 i__2 = abofdpos + i__ * ab_dim1;
 ab[i__2].r = abstmp; ab[i__2].i = 0.; // , expr subst  
 e[i__] = abstmp;
 if (abstmp != 0.) {
 z__1.r = tmp.r / abstmp; z__1.i = tmp.i / abstmp; // , expr subst  
 tmp.r = z__1.r; tmp.i = z__1.i; // , expr subst  
 }
 else {
 tmp.r = 1.; tmp.i = 0.; // , expr subst  
 }
 if (i__ < *n - 1) {
 i__2 = abofdpos + (i__ + 1) * ab_dim1;
 i__3 = abofdpos + (i__ + 1) * ab_dim1;
 z__1.r = ab[i__3].r * tmp.r - ab[i__3].i * tmp.i; z__1.i = ab[i__3].r * tmp.i + ab[i__3].i * tmp.r; // , expr subst  
 ab[i__2].r = z__1.r; ab[i__2].i = z__1.i; // , expr subst  
 }
 /* IF( WANTQ ) THEN */
 /* CALL ZSCAL( N, TMP, Q( 1, I+1 ), 1 ) */
 /* END IF */
 /* L70: */
 }
 }
 hous[1].r = 1.; hous[1].i = 0.; // , expr subst  
 work[1].r = 1.; work[1].i = 0.; // , expr subst  
 return 0;
 }
 /* Main code start here. */
 /* Reduce the hermitian band of A to a tridiagonal matrix. */
 thgrsiz = *n;
 grsiz = 1;
 shift = 3;
 /* NBTILES = CEILING( REAL(N)/REAL(KD) ) */
 nbtiles = *n / *kd;
 ceiltmp = *n - nbtiles * *kd;
 if (ceiltmp != 0) {
 ++nbtiles;
 }
 /* STEPERCOL = CEILING( REAL(SHIFT)/REAL(GRSIZ) ) */
 stepercol = shift / grsiz;
 ceiltmp = shift - stepercol * grsiz;
 if (ceiltmp != 0) {
 ++stepercol;
 }
 /* THGRNB = CEILING( REAL(N-1)/REAL(THGRSIZ) ) */
 thgrnb = (*n - 1) / thgrsiz;
 ceiltmp = *n - 1 - thgrnb * thgrsiz;
 if (ceiltmp != 0) {
 ++thgrnb;
 }
 i__1 = *kd + 1;
 zlacpy_("A", &i__1, n, &ab[ab_offset], ldab, &work[apos], &lda) ;
 zlaset_("A", kd, n, &c_b1, &c_b1, &work[awpos], &lda);
 /* openMP parallelisation start here */
 /* #if defined(_OPENMP) */
 /* !$OMP PARALLEL PRIVATE( TID, THGRID, BLKLASTIND ) */
 /* !$OMP$ PRIVATE( THED, I, M, K, ST, ED, STT, SWEEPID ) */
 /* !$OMP$ PRIVATE( MYID, TTYPE, COLPT, STIND, EDIND ) */
 /* !$OMP$ SHARED ( UPLO, WANTQ, INDV, INDTAU, HOUS, WORK) */
 /* !$OMP$ SHARED ( N, KD, IB, NBTILES, LDA, LDV, INDA ) */
 /* !$OMP$ SHARED ( STEPERCOL, THGRNB, THGRSIZ, GRSIZ, SHIFT ) */
 /* !$OMP MASTER */
 /* #endif */
 /* main bulge chasing loop */
 i__1 = thgrnb;
 for (thgrid = 1;
 thgrid <= i__1;
 ++thgrid) {
 stt = (thgrid - 1) * thgrsiz + 1;
 /* Computing MIN */
 i__2 = stt + thgrsiz - 1; i__3 = *n - 1; // , expr subst  
 thed = min(i__2,i__3);
 i__2 = *n - 1;
 for (i__ = stt;
 i__ <= i__2;
 ++i__) {
 ed = min(i__,thed);
 if (stt > ed) {
 break;
 }
 i__3 = stepercol;
 for (m = 1;
 m <= i__3;
 ++m) {
 st = stt;
 i__4 = ed;
 for (sweepid = st;
 sweepid <= i__4;
 ++sweepid) {
 i__5 = grsiz;
 for (k = 1;
 k <= i__5;
 ++k) {
 myid = (i__ - sweepid) * (stepercol * grsiz) + (m - 1) * grsiz + k;
 if (myid == 1) {
 ttype = 1;
 }
 else {
 ttype = myid % 2 + 2;
 }
 if (ttype == 2) {
 colpt = myid / 2 * *kd + sweepid;
 stind = colpt - *kd + 1;
 edind = min(colpt,*n);
 blklastind = colpt;
 }
 else {
 colpt = (myid + 1) / 2 * *kd + sweepid;
 stind = colpt - *kd + 1;
 edind = min(colpt,*n);
 if (stind >= edind - 1 && edind == *n) {
 blklastind = *n;
 }
 else {
 blklastind = 0;
 }
 }
 /* Call the kernel */
 /* #if defined(_OPENMP) */
 /* IF( TTYPE.NE.1 ) THEN */
 /* !$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1)) */
 /* !$OMP$ DEPEND(in:WORK(MYID-1)) */
 /* !$OMP$ DEPEND(out:WORK(MYID)) */
 /* TID = OMP_GET_THREAD_NUM() */
 /* CALL ZHB2ST_KERNELS( UPLO, WANTQ, TTYPE, */
 /* $ STIND, EDIND, SWEEPID, N, KD, IB, */
 /* $ WORK ( INDA ), LDA, */
 /* $ HOUS( INDV ), HOUS( INDTAU ), LDV, */
 /* $ WORK( INDW + TID*KD ) ) */
 /* !$OMP END TASK */
 /* ELSE */
 /* !$OMP TASK DEPEND(in:WORK(MYID+SHIFT-1)) */
 /* !$OMP$ DEPEND(out:WORK(MYID)) */
 /* TID = OMP_GET_THREAD_NUM() */
 /* CALL ZHB2ST_KERNELS( UPLO, WANTQ, TTYPE, */
 /* $ STIND, EDIND, SWEEPID, N, KD, IB, */
 /* $ WORK ( INDA ), LDA, */
 /* $ HOUS( INDV ), HOUS( INDTAU ), LDV, */
 /* $ WORK( INDW + TID*KD ) ) */
 /* !$OMP END TASK */
 /* ENDIF */
 /* #else */
 zhb2st_kernels_(uplo, &wantq, &ttype, &stind, &edind, &sweepid, n, kd, &ib, &work[inda], &lda, & hous[indv], &hous[indtau], &ldv, &work[indw + tid * *kd]);
 /* #endif */
 if (blklastind >= *n - 1) {
 ++stt;
 break;
 }
 /* L140: */
 }
 /* L130: */
 }
 /* L120: */
 }
 /* L110: */
 }
 /* L100: */
 }
 /* #if defined(_OPENMP) */
 /* !$OMP END MASTER */
 /* !$OMP END PARALLEL */
 /* #endif */
 /* Copy the diagonal from A to D. Note that D is REAL thus only */
 /* the Real part is needed, the imaginary part should be zero. */
 i__1 = *n;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = dpos + (i__ - 1) * lda;
 d__[i__] = work[i__2].r;
 /* L150: */
 }
 /* Copy the off diagonal from A to E. Note that E is REAL thus only */
 /* the Real part is needed, the imaginary part should be zero. */
 if (upper) {
 i__1 = *n - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = ofdpos + i__ * lda;
 e[i__] = work[i__2].r;
 /* L160: */
 }
 }
 else {
 i__1 = *n - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = ofdpos + (i__ - 1) * lda;
 e[i__] = work[i__2].r;
 /* L170: */
 }
 }
 hous[1].r = (doublereal) lhmin; hous[1].i = 0.; // , expr subst  
 work[1].r = (doublereal) lwmin; work[1].i = 0.; // , expr subst  
 return 0;
 /* End of ZHETRD_HB2ST */
 }
 /* zhetrd_hb2st__ */
 