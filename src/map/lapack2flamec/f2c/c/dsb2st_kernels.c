/* ../netlib/v3.9.0/dsb2st_kernels.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static integer c__1 = 1;
 /* > \brief \b DSB2ST_KERNELS */
 /* @generated from zhb2st_kernels.f, fortran z -> d, Wed Dec 7 08:22:39 2016 */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DSB2ST_KERNELS + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsb2st_ kernels.f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsb2st_ kernels.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsb2st_ kernels.f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE DSB2ST_KERNELS( UPLO, WANTZ, TTYPE, */
 /* ST, ED, SWEEP, N, NB, IB, */
 /* A, LDA, V, TAU, LDVT, WORK) */
 /* IMPLICIT NONE */
 /* .. Scalar Arguments .. */
 /* CHARACTER UPLO */
 /* LOGICAL WANTZ */
 /* INTEGER TTYPE, ST, ED, SWEEP, N, NB, IB, LDA, LDVT */
 /* .. */
 /* .. Array Arguments .. */
 /* DOUBLE PRECISION A( LDA, * ), V( * ), */
 /* TAU( * ), WORK( * ) */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DSB2ST_KERNELS is an internal routine used by the DSYTRD_SB2ST */
 /* > subroutine. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] UPLO */
 /* > \verbatim */
 /* > UPLO is CHARACTER*1 */
 /* > \endverbatim */
 /* > */
 /* > \param[in] WANTZ */
 /* > \verbatim */
 /* > WANTZ is LOGICAL which indicate if Eigenvalue are requested or both */
 /* > Eigenvalue/Eigenvectors. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] TTYPE */
 /* > \verbatim */
 /* > TTYPE is INTEGER */
 /* > \endverbatim */
 /* > */
 /* > \param[in] ST */
 /* > \verbatim */
 /* > ST is INTEGER */
 /* > internal parameter for indices. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] ED */
 /* > \verbatim */
 /* > ED is INTEGER */
 /* > internal parameter for indices. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] SWEEP */
 /* > \verbatim */
 /* > SWEEP is INTEGER */
 /* > internal parameter for indices. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER. The order of the matrix A. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] NB */
 /* > \verbatim */
 /* > NB is INTEGER. The size of the band. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] IB */
 /* > \verbatim */
 /* > IB is INTEGER. */
 /* > \endverbatim */
 /* > */
 /* > \param[in, out] A */
 /* > \verbatim */
 /* > A is DOUBLE PRECISION array. A pointer to the matrix A. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER. The leading dimension of the matrix A. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] V */
 /* > \verbatim */
 /* > V is DOUBLE PRECISION array, dimension 2*n if eigenvalues only are */
 /* > requested or to be queried for vectors. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] TAU */
 /* > \verbatim */
 /* > TAU is DOUBLE PRECISION array, dimension (2*n). */
 /* > The scalar factors of the Householder reflectors are stored */
 /* > in this array. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDVT */
 /* > \verbatim */
 /* > LDVT is INTEGER. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is DOUBLE PRECISION array. Workspace of size nb. */
 /* > \endverbatim */
 /* > */
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
 int dsb2st_kernels_(char *uplo, logical *wantz, integer * ttype, integer *st, integer *ed, integer *sweep, integer *n, integer * nb, integer *ib, doublereal *a, integer *lda, doublereal *v, doublereal *tau, integer *ldvt, doublereal *work) {
 /* System generated locals */
 integer a_dim1, a_offset, i__1, i__2;
 doublereal d__1;
 /* Local variables */
 integer i__, j1, j2, lm, ln;
 doublereal ctmp;
 integer dpos, vpos;
 extern logical lsame_(char *, char *);
 logical upper;
 extern /* Subroutine */
 int dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *);
 integer ajeter;
 extern /* Subroutine */
 int dlarfx_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *), dlarfy_(char *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *);
 integer ofdpos, taupos;
 /* -- LAPACK computational routine (version 3.7.1) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* June 2017 */
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
 /* .. External Functions .. */
 /* .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 --v;
 --tau;
 --work;
 /* Function Body */
 ajeter = *ib + *ldvt;
 upper = lsame_(uplo, "U");
 if (upper) {
 dpos = (*nb << 1) + 1;
 ofdpos = *nb << 1;
 }
 else {
 dpos = 1;
 ofdpos = 2;
 }
 /* Upper case */
 if (upper) {
 if (*wantz) {
 vpos = (*sweep - 1) % 2 * *n + *st;
 taupos = (*sweep - 1) % 2 * *n + *st;
 }
 else {
 vpos = (*sweep - 1) % 2 * *n + *st;
 taupos = (*sweep - 1) % 2 * *n + *st;
 }
 if (*ttype == 1) {
 lm = *ed - *st + 1;
 v[vpos] = 1.;
 i__1 = lm - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 v[vpos + i__] = a[ofdpos - i__ + (*st + i__) * a_dim1];
 a[ofdpos - i__ + (*st + i__) * a_dim1] = 0.;
 /* L10: */
 }
 ctmp = a[ofdpos + *st * a_dim1];
 dlarfg_(&lm, &ctmp, &v[vpos + 1], &c__1, &tau[taupos]);
 a[ofdpos + *st * a_dim1] = ctmp;
 lm = *ed - *st + 1;
 d__1 = tau[taupos];
 i__1 = *lda - 1;
 dlarfy_(uplo, &lm, &v[vpos], &c__1, &d__1, &a[dpos + *st * a_dim1] , &i__1, &work[1]);
 }
 if (*ttype == 3) {
 lm = *ed - *st + 1;
 d__1 = tau[taupos];
 i__1 = *lda - 1;
 dlarfy_(uplo, &lm, &v[vpos], &c__1, &d__1, &a[dpos + *st * a_dim1] , &i__1, &work[1]);
 }
 if (*ttype == 2) {
 j1 = *ed + 1;
 /* Computing MIN */
 i__1 = *ed + *nb;
 j2 = min(i__1,*n);
 ln = *ed - *st + 1;
 lm = j2 - j1 + 1;
 if (lm > 0) {
 d__1 = tau[taupos];
 i__1 = *lda - 1;
 dlarfx_("Left", &ln, &lm, &v[vpos], &d__1, &a[dpos - *nb + j1 * a_dim1], &i__1, &work[1]);
 if (*wantz) {
 vpos = (*sweep - 1) % 2 * *n + j1;
 taupos = (*sweep - 1) % 2 * *n + j1;
 }
 else {
 vpos = (*sweep - 1) % 2 * *n + j1;
 taupos = (*sweep - 1) % 2 * *n + j1;
 }
 v[vpos] = 1.;
 i__1 = lm - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 v[vpos + i__] = a[dpos - *nb - i__ + (j1 + i__) * a_dim1];
 a[dpos - *nb - i__ + (j1 + i__) * a_dim1] = 0.;
 /* L30: */
 }
 ctmp = a[dpos - *nb + j1 * a_dim1];
 dlarfg_(&lm, &ctmp, &v[vpos + 1], &c__1, &tau[taupos]);
 a[dpos - *nb + j1 * a_dim1] = ctmp;
 i__1 = ln - 1;
 i__2 = *lda - 1;
 dlarfx_("Right", &i__1, &lm, &v[vpos], &tau[taupos], &a[dpos - *nb + 1 + j1 * a_dim1], &i__2, &work[1]);
 }
 }
 /* Lower case */
 }
 else {
 if (*wantz) {
 vpos = (*sweep - 1) % 2 * *n + *st;
 taupos = (*sweep - 1) % 2 * *n + *st;
 }
 else {
 vpos = (*sweep - 1) % 2 * *n + *st;
 taupos = (*sweep - 1) % 2 * *n + *st;
 }
 if (*ttype == 1) {
 lm = *ed - *st + 1;
 v[vpos] = 1.;
 i__1 = lm - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 v[vpos + i__] = a[ofdpos + i__ + (*st - 1) * a_dim1];
 a[ofdpos + i__ + (*st - 1) * a_dim1] = 0.;
 /* L20: */
 }
 dlarfg_(&lm, &a[ofdpos + (*st - 1) * a_dim1], &v[vpos + 1], &c__1, &tau[taupos]);
 lm = *ed - *st + 1;
 d__1 = tau[taupos];
 i__1 = *lda - 1;
 dlarfy_(uplo, &lm, &v[vpos], &c__1, &d__1, &a[dpos + *st * a_dim1] , &i__1, &work[1]);
 }
 if (*ttype == 3) {
 lm = *ed - *st + 1;
 d__1 = tau[taupos];
 i__1 = *lda - 1;
 dlarfy_(uplo, &lm, &v[vpos], &c__1, &d__1, &a[dpos + *st * a_dim1] , &i__1, &work[1]);
 }
 if (*ttype == 2) {
 j1 = *ed + 1;
 /* Computing MIN */
 i__1 = *ed + *nb;
 j2 = min(i__1,*n);
 ln = *ed - *st + 1;
 lm = j2 - j1 + 1;
 if (lm > 0) {
 i__1 = *lda - 1;
 dlarfx_("Right", &lm, &ln, &v[vpos], &tau[taupos], &a[dpos + * nb + *st * a_dim1], &i__1, &work[1]);
 if (*wantz) {
 vpos = (*sweep - 1) % 2 * *n + j1;
 taupos = (*sweep - 1) % 2 * *n + j1;
 }
 else {
 vpos = (*sweep - 1) % 2 * *n + j1;
 taupos = (*sweep - 1) % 2 * *n + j1;
 }
 v[vpos] = 1.;
 i__1 = lm - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 v[vpos + i__] = a[dpos + *nb + i__ + *st * a_dim1];
 a[dpos + *nb + i__ + *st * a_dim1] = 0.;
 /* L40: */
 }
 dlarfg_(&lm, &a[dpos + *nb + *st * a_dim1], &v[vpos + 1], & c__1, &tau[taupos]);
 i__1 = ln - 1;
 d__1 = tau[taupos];
 i__2 = *lda - 1;
 dlarfx_("Left", &lm, &i__1, &v[vpos], &d__1, &a[dpos + *nb - 1 + (*st + 1) * a_dim1], &i__2, &work[1]);
 }
 }
 }
 return 0;
 /* END OF DSB2ST_KERNELS */
 }
 /* dsb2st_kernels__ */
 