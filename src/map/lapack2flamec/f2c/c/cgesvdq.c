/* ../netlib/v3.9.0/cgesvdq.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static complex c_b1 = {
0.f,0.f}
;
 static complex c_b2 = {
1.f,0.f}
;
 static integer c_n1 = -1;
 static integer c__1 = 1;
 static real c_b74 = 0.f;
 static integer c__0 = 0;
 static real c_b87 = 1.f;
 static logical c_false = FALSE_;
 /* > \brief <b> CGESVDQ computes the singular value decomposition (SVD) with a QR-Preconditioned QR SVD Method for GE matrices</b> */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download CGESVDQ + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesvdq .f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesvdq .f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesvdq .f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE CGESVDQ( JOBA, JOBP, JOBR, JOBU, JOBV, M, N, A, LDA, */
 /* S, U, LDU, V, LDV, NUMRANK, IWORK, LIWORK, */
 /* CWORK, LCWORK, RWORK, LRWORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* IMPLICIT NONE */
 /* CHARACTER JOBA, JOBP, JOBR, JOBU, JOBV */
 /* INTEGER M, N, LDA, LDU, LDV, NUMRANK, LIWORK, LCWORK, LRWORK, */
 /* INFO */
 /* .. */
 /* .. Array Arguments .. */
 /* COMPLEX A( LDA, * ), U( LDU, * ), V( LDV, * ), CWORK( * ) */
 /* REAL S( * ), RWORK( * ) */
 /* INTEGER IWORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > CGESVDQ computes the singular value decomposition (SVD) of a complex */
 /* > M-by-N matrix A, where M >= N. The SVD of A is written as */
 /* > [++] [xx] [x0] [xx] */
 /* > A = U * SIGMA * V^*, [++] = [xx] * [ox] * [xx] */
 /* > [++] [xx] */
 /* > where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal */
 /* > matrix, and V is an N-by-N unitary matrix. The diagonal elements */
 /* > of SIGMA are the singular values of A. The columns of U and V are the */
 /* > left and the right singular vectors of A, respectively. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] JOBA */
 /* > \verbatim */
 /* > JOBA is CHARACTER*1 */
 /* > Specifies the level of accuracy in the computed SVD */
 /* > = 'A' The requested accuracy corresponds to having the backward */
 /* > error bounded by || delta A ||_F <= f(m,n) * EPS * || A ||_F, */
 /* > where EPS = SLAMCH('Epsilon'). This authorises CGESVDQ to */
 /* > truncate the computed triangular factor in a rank revealing */
 /* > QR factorization whenever the truncated part is below the */
 /* > threshold of the order of EPS * ||A||_F. This is aggressive */
 /* > truncation level. */
 /* > = 'M' Similarly as with 'A', but the truncation is more gentle: it */
 /* > is allowed only when there is a drop on the diagonal of the */
 /* > triangular factor in the QR factorization. This is medium */
 /* > truncation level. */
 /* > = 'H' High accuracy requested. No numerical rank determination based */
 /* > on the rank revealing QR factorization is attempted. */
 /* > = 'E' Same as 'H', and in addition the condition number of column */
 /* > scaled A is estimated and returned in RWORK(1). */
 /* > N^(-1/4)*RWORK(1) <= ||pinv(A_scaled)||_2 <= N^(1/4)*RWORK(1) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] JOBP */
 /* > \verbatim */
 /* > JOBP is CHARACTER*1 */
 /* > = 'P' The rows of A are ordered in decreasing order with respect to */
 /* > ||A(i,:)||_\infty. This enhances numerical accuracy at the cost */
 /* > of extra data movement. Recommended for numerical robustness. */
 /* > = 'N' No row pivoting. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] JOBR */
 /* > \verbatim */
 /* > JOBR is CHARACTER*1 */
 /* > = 'T' After the initial pivoted QR factorization, CGESVD is applied to */
 /* > the adjoint R**H of the computed triangular factor R. This involves */
 /* > some extra data movement (matrix transpositions). Useful for */
 /* > experiments, research and development. */
 /* > = 'N' The triangular factor R is given as input to CGESVD. This may be */
 /* > preferred as it involves less data movement. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] JOBU */
 /* > \verbatim */
 /* > JOBU is CHARACTER*1 */
 /* > = 'A' All M left singular vectors are computed and returned in the */
 /* > matrix U. See the description of U. */
 /* > = 'S' or 'U' N = min(M,N) left singular vectors are computed and returned */
 /* > in the matrix U. See the description of U. */
 /* > = 'R' Numerical rank NUMRANK is determined and only NUMRANK left singular */
 /* > vectors are computed and returned in the matrix U. */
 /* > = 'F' The N left singular vectors are returned in factored form as the */
 /* > product of the Q factor from the initial QR factorization and the */
 /* > N left singular vectors of (R**H , 0)**H. If row pivoting is used, */
 /* > then the necessary information on the row pivoting is stored in */
 /* > IWORK(N+1:N+M-1). */
 /* > = 'N' The left singular vectors are not computed. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] JOBV */
 /* > \verbatim */
 /* > JOBV is CHARACTER*1 */
 /* > = 'A', 'V' All N right singular vectors are computed and returned in */
 /* > the matrix V. */
 /* > = 'R' Numerical rank NUMRANK is determined and only NUMRANK right singular */
 /* > vectors are computed and returned in the matrix V. This option is */
 /* > allowed only if JOBU = 'R' or JOBU = 'N';
 otherwise it is illegal. */
 /* > = 'N' The right singular vectors are not computed. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] M */
 /* > \verbatim */
 /* > M is INTEGER */
 /* > The number of rows of the input matrix A. M >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The number of columns of the input matrix A. M >= N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] A */
 /* > \verbatim */
 /* > A is COMPLEX array of dimensions LDA x N */
 /* > On entry, the input matrix A. */
 /* > On exit, if JOBU .NE. 'N' or JOBV .NE. 'N', the lower triangle of A contains */
 /* > the Householder vectors as stored by CGEQP3. If JOBU = 'F', these Householder */
 /* > vectors together with CWORK(1:N) can be used to restore the Q factors from */
 /* > the initial pivoted QR factorization of A. See the description of U. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER. */
 /* > The leading dimension of the array A. LDA >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] S */
 /* > \verbatim */
 /* > S is REAL array of dimension N. */
 /* > The singular values of A, ordered so that S(i) >= S(i+1). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] U */
 /* > \verbatim */
 /* > U is COMPLEX array, dimension */
 /* > LDU x M if JOBU = 'A';
 see the description of LDU. In this case, */
 /* > on exit, U contains the M left singular vectors. */
 /* > LDU x N if JOBU = 'S', 'U', 'R' ;
 see the description of LDU. In this */
 /* > case, U contains the leading N or the leading NUMRANK left singular vectors. */
 /* > LDU x N if JOBU = 'F' ;
 see the description of LDU. In this case U */
 /* > contains N x N unitary matrix that can be used to form the left */
 /* > singular vectors. */
 /* > If JOBU = 'N', U is not referenced. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDU */
 /* > \verbatim */
 /* > LDU is INTEGER. */
 /* > The leading dimension of the array U. */
 /* > If JOBU = 'A', 'S', 'U', 'R', LDU >= max(1,M). */
 /* > If JOBU = 'F', LDU >= max(1,N). */
 /* > Otherwise, LDU >= 1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] V */
 /* > \verbatim */
 /* > V is COMPLEX array, dimension */
 /* > LDV x N if JOBV = 'A', 'V', 'R' or if JOBA = 'E' . */
 /* > If JOBV = 'A', or 'V', V contains the N-by-N unitary matrix V**H;
 */
 /* > If JOBV = 'R', V contains the first NUMRANK rows of V**H (the right */
 /* > singular vectors, stored rowwise, of the NUMRANK largest singular values). */
 /* > If JOBV = 'N' and JOBA = 'E', V is used as a workspace. */
 /* > If JOBV = 'N', and JOBA.NE.'E', V is not referenced. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDV */
 /* > \verbatim */
 /* > LDV is INTEGER */
 /* > The leading dimension of the array V. */
 /* > If JOBV = 'A', 'V', 'R', or JOBA = 'E', LDV >= max(1,N). */
 /* > Otherwise, LDV >= 1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] NUMRANK */
 /* > \verbatim */
 /* > NUMRANK is INTEGER */
 /* > NUMRANK is the numerical rank first determined after the rank */
 /* > revealing QR factorization, following the strategy specified by the */
 /* > value of JOBA. If JOBV = 'R' and JOBU = 'R', only NUMRANK */
 /* > leading singular values and vectors are then requested in the call */
 /* > of CGESVD. The final value of NUMRANK might be further reduced if */
 /* > some singular values are computed as zeros. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IWORK */
 /* > \verbatim */
 /* > IWORK is INTEGER array, dimension (max(1, LIWORK)). */
 /* > On exit, IWORK(1:N) contains column pivoting permutation of the */
 /* > rank revealing QR factorization. */
 /* > If JOBP = 'P', IWORK(N+1:N+M-1) contains the indices of the sequence */
 /* > of row swaps used in row pivoting. These can be used to restore the */
 /* > left singular vectors in the case JOBU = 'F'. */
 /* > */
 /* > If LIWORK, LCWORK, or LRWORK = -1, then on exit, if INFO = 0, */
 /* > LIWORK(1) returns the minimal LIWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LIWORK */
 /* > \verbatim */
 /* > LIWORK is INTEGER */
 /* > The dimension of the array IWORK. */
 /* > LIWORK >= N + M - 1, if JOBP = 'P';
 */
 /* > LIWORK >= N if JOBP = 'N'. */
 /* > */
 /* > If LIWORK = -1, then a workspace query is assumed;
 the routine */
 /* > only calculates and returns the optimal and minimal sizes */
 /* > for the CWORK, IWORK, and RWORK arrays, and no error */
 /* > message related to LCWORK is issued by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] CWORK */
 /* > \verbatim */
 /* > CWORK is COMPLEX array, dimension (max(2, LCWORK)), used as a workspace. */
 /* > On exit, if, on entry, LCWORK.NE.-1, CWORK(1:N) contains parameters */
 /* > needed to recover the Q factor from the QR factorization computed by */
 /* > CGEQP3. */
 /* > */
 /* > If LIWORK, LCWORK, or LRWORK = -1, then on exit, if INFO = 0, */
 /* > CWORK(1) returns the optimal LCWORK, and */
 /* > CWORK(2) returns the minimal LCWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] LCWORK */
 /* > \verbatim */
 /* > LCWORK is INTEGER */
 /* > The dimension of the array CWORK. It is determined as follows: */
 /* > Let LWQP3 = N+1, LWCON = 2*N, and let */
 /* > LWUNQ = {
 MAX( N, 1 ), if JOBU = 'R', 'S', or 'U' */
 /* > {
 MAX( M, 1 ), if JOBU = 'A' */
 /* > LWSVD = MAX( 3*N, 1 ) */
 /* > LWLQF = MAX( N/2, 1 ), LWSVD2 = MAX( 3*(N/2), 1 ), LWUNLQ = MAX( N, 1 ), */
 /* > LWQRF = MAX( N/2, 1 ), LWUNQ2 = MAX( N, 1 ) */
 /* > Then the minimal value of LCWORK is: */
 /* > = MAX( N + LWQP3, LWSVD ) if only the singular values are needed;
 */
 /* > = MAX( N + LWQP3, LWCON, LWSVD ) if only the singular values are needed, */
 /* > and a scaled condition estimate requested;
 */
 /* > */
 /* > = N + MAX( LWQP3, LWSVD, LWUNQ ) if the singular values and the left */
 /* > singular vectors are requested;
 */
 /* > = N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ) if the singular values and the left */
 /* > singular vectors are requested, and also */
 /* > a scaled condition estimate requested;
 */
 /* > */
 /* > = N + MAX( LWQP3, LWSVD ) if the singular values and the right */
 /* > singular vectors are requested;
 */
 /* > = N + MAX( LWQP3, LWCON, LWSVD ) if the singular values and the right */
 /* > singular vectors are requested, and also */
 /* > a scaled condition etimate requested;
 */
 /* > */
 /* > = N + MAX( LWQP3, LWSVD, LWUNQ ) if the full SVD is requested with JOBV = 'R';
 */
 /* > independent of JOBR;
 */
 /* > = N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ) if the full SVD is requested, */
 /* > JOBV = 'R' and, also a scaled condition */
 /* > estimate requested;
 independent of JOBR;
 */
 /* > = MAX( N + MAX( LWQP3, LWSVD, LWUNQ ), */
 /* > N + MAX( LWQP3, N/2+LWLQF, N/2+LWSVD2, N/2+LWUNLQ, LWUNQ) ) if the */
 /* > full SVD is requested with JOBV = 'A' or 'V', and */
 /* > JOBR ='N' */
 /* > = MAX( N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ), */
 /* > N + MAX( LWQP3, LWCON, N/2+LWLQF, N/2+LWSVD2, N/2+LWUNLQ, LWUNQ ) ) */
 /* > if the full SVD is requested with JOBV = 'A' or 'V', and */
 /* > JOBR ='N', and also a scaled condition number estimate */
 /* > requested. */
 /* > = MAX( N + MAX( LWQP3, LWSVD, LWUNQ ), */
 /* > N + MAX( LWQP3, N/2+LWQRF, N/2+LWSVD2, N/2+LWUNQ2, LWUNQ ) ) if the */
 /* > full SVD is requested with JOBV = 'A', 'V', and JOBR ='T' */
 /* > = MAX( N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ), */
 /* > N + MAX( LWQP3, LWCON, N/2+LWQRF, N/2+LWSVD2, N/2+LWUNQ2, LWUNQ ) ) */
 /* > if the full SVD is requested with JOBV = 'A', 'V' and */
 /* > JOBR ='T', and also a scaled condition number estimate */
 /* > requested. */
 /* > Finally, LCWORK must be at least two: LCWORK = MAX( 2, LCWORK ). */
 /* > */
 /* > If LCWORK = -1, then a workspace query is assumed;
 the routine */
 /* > only calculates and returns the optimal and minimal sizes */
 /* > for the CWORK, IWORK, and RWORK arrays, and no error */
 /* > message related to LCWORK is issued by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] RWORK */
 /* > \verbatim */
 /* > RWORK is REAL array, dimension (max(1, LRWORK)). */
 /* > On exit, */
 /* > 1. If JOBA = 'E', RWORK(1) contains an estimate of the condition */
 /* > number of column scaled A. If A = C * D where D is diagonal and C */
 /* > has unit columns in the Euclidean norm, then, assuming full column rank, */
 /* > N^(-1/4) * RWORK(1) <= ||pinv(C)||_2 <= N^(1/4) * RWORK(1). */
 /* > Otherwise, RWORK(1) = -1. */
 /* > 2. RWORK(2) contains the number of singular values computed as */
 /* > exact zeros in CGESVD applied to the upper triangular or trapeziodal */
 /* > R (from the initial QR factorization). In case of early exit (no call to */
 /* > CGESVD, such as in the case of zero matrix) RWORK(2) = -1. */
 /* > */
 /* > If LIWORK, LCWORK, or LRWORK = -1, then on exit, if INFO = 0, */
 /* > RWORK(1) returns the minimal LRWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LRWORK */
 /* > \verbatim */
 /* > LRWORK is INTEGER. */
 /* > The dimension of the array RWORK. */
 /* > If JOBP ='P', then LRWORK >= MAX(2, M, 5*N);
 */
 /* > Otherwise, LRWORK >= MAX(2, 5*N). */
 /* > */
 /* > If LRWORK = -1, then a workspace query is assumed;
 the routine */
 /* > only calculates and returns the optimal and minimal sizes */
 /* > for the CWORK, IWORK, and RWORK arrays, and no error */
 /* > message related to LCWORK is issued by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit. */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value. */
 /* > > 0: if CBDSQR did not converge, INFO specifies how many superdiagonals */
 /* > of an intermediate bidiagonal form B (computed in CGESVD) did not */
 /* > converge to zero. */
 /* > \endverbatim */
 /* > \par Further Details: */
 /* ======================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > 1. The data movement (matrix transpose) is coded using simple nested */
 /* > DO-loops because BLAS and LAPACK do not provide corresponding subroutines. */
 /* > Those DO-loops are easily identified in this source code - by the CONTINUE */
 /* > statements labeled with 11**. In an optimized version of this code, the */
 /* > nested DO loops should be replaced with calls to an optimized subroutine. */
 /* > 2. This code scales A by 1/SQRT(M) if the largest ABS(A(i,j)) could cause */
 /* > column norm overflow. This is the minial precaution and it is left to the */
 /* > SVD routine (CGESVD) to do its own preemptive scaling if potential over- */
 /* > or underflows are detected. To avoid repeated scanning of the array A, */
 /* > an optimal implementation would do all necessary scaling before calling */
 /* > CGESVD and the scaling in CGESVD can be switched off. */
 /* > 3. Other comments related to code optimization are given in comments in the */
 /* > code, enlosed in [[double brackets]]. */
 /* > \endverbatim */
 /* > \par Bugs, examples and comments */
 /* =========================== */
 /* > \verbatim */
 /* > Please report all bugs and send interesting examples and/or comments to */
 /* > drmac@math.hr. Thank you. */
 /* > \endverbatim */
 /* > \par References */
 /* =============== */
 /* > \verbatim */
 /* > [1] Zlatko Drmac, Algorithm 977: A QR-Preconditioned QR SVD Method for */
 /* > Computing the SVD with High Accuracy. ACM Trans. Math. Softw. */
 /* > 44(1): 11:1-11:30 (2017) */
 /* > */
 /* > SIGMA library, xGESVDQ section updated February 2016. */
 /* > Developed and coded by Zlatko Drmac, Department of Mathematics */
 /* > University of Zagreb, Croatia, drmac@math.hr */
 /* > \endverbatim */
 /* > \par Contributors: */
 /* ================== */
 /* > */
 /* > \verbatim */
 /* > Developed and coded by Zlatko Drmac, Department of Mathematics */
 /* > University of Zagreb, Croatia, drmac@math.hr */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2018 */
 /* > \ingroup complexGEsing */
 /* ===================================================================== */
 /* Subroutine */
 int cgesvdq_(char *joba, char *jobp, char *jobr, char *jobu, char *jobv, integer *m, integer *n, complex *a, integer *lda, real *s, complex *u, integer *ldu, complex *v, integer *ldv, integer *numrank, integer *iwork, integer *liwork, complex *cwork, integer *lcwork, real *rwork, integer *lrwork, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, i__3;
 real r__1;
 complex q__1;
 /* Builtin functions */
 double sqrt(doublereal), c_abs(complex *);
 void r_cnjg(complex *, complex *);
 /* Local variables */
 integer optratio, p, q, n1, nr;
 real big;
 integer lwrk_cgeqp3__, lwrk_cgelqf__, lwrk_cgeqrf__, lwrk_cgesvd__, lwrk_cunmlq__, lwrk_cunmqr__, lwrk_cgesvd2__, ierr;
 complex ctmp;
 real rtmp;
 integer lwrk_cunmqr2__;
 logical lsvc0, accla;
 integer lwqp3;
 logical acclh, acclm, conda;
 extern logical lsame_(char *, char *);
 logical lsvec;
 real sfmin, epsln;
 integer lwcon;
 logical rsvec;
 integer lwlqf, lwqrf;
 logical wntua;
 integer lwsvd;
 logical wntva, dntwu, dntwv, wntuf;
 integer lwunq;
 logical wntur, wntus, wntvr;
 extern /* Subroutine */
 int cgeqp3_(integer *, integer *, complex *, integer *, integer *, complex *, complex *, integer *, real *, integer *);
 extern real scnrm2_(integer *, complex *, integer *);
 integer lwsvd2, lwunq2;
 extern real clange_(char *, integer *, integer *, complex *, integer *, real *);
 extern /* Subroutine */
 int cgelqf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), clascl_( char *, integer *, integer *, real *, real *, integer *, integer * , complex *, integer *, integer *);
 real sconda;
 extern /* Subroutine */
 int cgeqrf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), csscal_( integer *, real *, complex *, integer *);
 extern real slamch_(char *);
 extern /* Subroutine */
 int cgesvd_(char *, char *, integer *, integer *, complex *, integer *, real *, complex *, integer *, complex *, integer *, complex *, integer *, real *, integer *), clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *), xerbla_(char *, integer *), clapmt_(logical *, integer *, integer *, complex *, integer *, integer *), slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real * , integer *, integer *), cpocon_(char *, integer *, complex *, integer *, real *, real *, complex *, real *, integer * );
 extern integer isamax_(integer *, real *, integer *);
 extern /* Subroutine */
 int claswp_(integer *, complex *, integer *, integer *, integer *, integer *, integer *), slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
 complex cdummy[1];
 extern /* Subroutine */
 int cunmlq_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *), cunmqr_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *);
 integer minwrk;
 logical rtrans;
 real rdummy[1];
 logical lquery;
 integer lwunlq, optwrk;
 logical rowprm;
 integer minwrk2;
 logical ascaled;
 integer optwrk2, iminwrk, rminwrk;
 /* .. Scalar Arguments .. */
 /* .. */
 /* .. Array Arguments .. */
 /* ===================================================================== */
 /* .. Parameters .. */
 /* .. */
 /* .. Local Scalars .. */
 /* .. */
 /* .. Local Arrays */
 /* .. */
 /* .. External Subroutines (BLAS, LAPACK) */
 /* .. */
 /* .. External Functions (BLAS, LAPACK) */
 /* .. */
 /* .. Intrinsic Functions .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Test the input arguments */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 --s;
 u_dim1 = *ldu;
 u_offset = 1 + u_dim1;
 u -= u_offset;
 v_dim1 = *ldv;
 v_offset = 1 + v_dim1;
 v -= v_offset;
 --iwork;
 --cwork;
 --rwork;
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
 /* Function Body */
#if AOCL_DTL_LOG_ENABLE
 char buffer[256];
 snprintf(buffer, 256, "cgesvdq inputs: joba %c, jobp %c, jobr %c, jobu %c, jobv %c, m %d, n %d, lda %d, ldu %d, ldv %d, numrank %d, liwork %d, lcwork %d, lrwork %d\n", *joba, *jobp, *jobr, *jobu, *jobv, *m, *n, *lda, *ldu, *ldv, *numrank, *liwork, *lcwork, *lrwork);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 wntus = lsame_(jobu, "S") || lsame_(jobu, "U");
 wntur = lsame_(jobu, "R");
 wntua = lsame_(jobu, "A");
 wntuf = lsame_(jobu, "F");
 lsvc0 = wntus || wntur || wntua;
 lsvec = lsvc0 || wntuf;
 dntwu = lsame_(jobu, "N");
 wntvr = lsame_(jobv, "R");
 wntva = lsame_(jobv, "A") || lsame_(jobv, "V");
 rsvec = wntvr || wntva;
 dntwv = lsame_(jobv, "N");
 accla = lsame_(joba, "A");
 acclm = lsame_(joba, "M");
 conda = lsame_(joba, "E");
 acclh = lsame_(joba, "H") || conda;
 rowprm = lsame_(jobp, "P");
 rtrans = lsame_(jobr, "T");
 if (rowprm) {
 /* Computing MAX */
 i__1 = 1; i__2 = *n + *m - 1; // , expr subst  
 iminwrk = max(i__1,i__2);
 /* Computing MAX */
 i__1 = max(2,*m); i__2 = *n * 5; // , expr subst  
 rminwrk = max(i__1,i__2);
 }
 else {
 iminwrk = max(1,*n);
 /* Computing MAX */
 i__1 = 2; i__2 = *n * 5; // , expr subst  
 rminwrk = max(i__1,i__2);
 }
 lquery = *liwork == -1 || *lcwork == -1 || *lrwork == -1;
 *info = 0;
 if (! (accla || acclm || acclh)) {
 *info = -1;
 }
 else if (! (rowprm || lsame_(jobp, "N"))) {
 *info = -2;
 }
 else if (! (rtrans || lsame_(jobr, "N"))) {
 *info = -3;
 }
 else if (! (lsvec || dntwu)) {
 *info = -4;
 }
 else if (wntur && wntva) {
 *info = -5;
 }
 else if (! (rsvec || dntwv)) {
 *info = -5;
 }
 else if (*m < 0) {
 *info = -6;
 }
 else if (*n < 0 || *n > *m) {
 *info = -7;
 }
 else if (*lda < max(1,*m)) {
 *info = -9;
 }
 else if (*ldu < 1 || lsvc0 && *ldu < *m || wntuf && *ldu < *n) {
 *info = -12;
 }
 else if (*ldv < 1 || rsvec && *ldv < *n || conda && *ldv < *n) {
 *info = -14;
 }
 else if (*liwork < iminwrk && ! lquery) {
 *info = -17;
 }
 if (*info == 0) {
 /* Compute workspace */
 /* .. compute the minimal and the optimal workspace lengths */
 /* [[The expressions for computing the minimal and the optimal */
 /* values of LCWORK are written with a lot of redundancy and */
 /* can be simplified. However, this detailed form is easier for */
 /* maintenance and modifications of the code.]] */
 /* .. minimal workspace length for CGEQP3 of an M x N matrix */
 lwqp3 = *n + 1;
 /* .. minimal workspace length for CUNMQR to build left singular vectors */
 if (wntus || wntur) {
 lwunq = max(*n,1);
 }
 else if (wntua) {
 lwunq = max(*m,1);
 }
 /* .. minimal workspace length for CPOCON of an N x N matrix */
 lwcon = *n << 1;
 /* .. CGESVD of an N x N matrix */
 /* Computing MAX */
 i__1 = *n * 3;
 lwsvd = max(i__1,1);
 if (lquery) {
 cgeqp3_(m, n, &a[a_offset], lda, &iwork[1], cdummy, cdummy, &c_n1, rdummy, &ierr);
 lwrk_cgeqp3__ = (integer) cdummy[0].r;
 if (wntus || wntur) {
 cunmqr_("L", "N", m, n, n, &a[a_offset], lda, cdummy, &u[ u_offset], ldu, cdummy, &c_n1, &ierr);
 lwrk_cunmqr__ = (integer) cdummy[0].r;
 }
 else if (wntua) {
 cunmqr_("L", "N", m, m, n, &a[a_offset], lda, cdummy, &u[ u_offset], ldu, cdummy, &c_n1, &ierr);
 lwrk_cunmqr__ = (integer) cdummy[0].r;
 }
 else {
 lwrk_cunmqr__ = 0;
 }
 }
 minwrk = 2;
 optwrk = 2;
 if (! (lsvec || rsvec)) {
 /* .. minimal and optimal sizes of the complex workspace if */
 /* only the singular values are requested */
 if (conda) {
 /* Computing MAX */
 i__1 = *n + lwqp3; i__1 = max(i__1,lwcon); // , expr subst  
 minwrk = max(i__1,lwsvd);
 }
 else {
 /* Computing MAX */
 i__1 = *n + lwqp3;
 minwrk = max(i__1,lwsvd);
 }
 if (lquery) {
 cgesvd_("N", "N", n, n, &a[a_offset], lda, &s[1], &u[u_offset] , ldu, &v[v_offset], ldv, cdummy, &c_n1, rdummy, & ierr);
 lwrk_cgesvd__ = (integer) cdummy[0].r;
 if (conda) {
 /* Computing MAX */
 i__1 = *n + lwrk_cgeqp3__; i__2 = *n + lwcon; i__1 = max( i__1,i__2); // ; expr subst  
 optwrk = max(i__1,lwrk_cgesvd__);
 }
 else {
 /* Computing MAX */
 i__1 = *n + lwrk_cgeqp3__;
 optwrk = max(i__1,lwrk_cgesvd__);
 }
 }
 }
 else if (lsvec && ! rsvec) {
 /* .. minimal and optimal sizes of the complex workspace if the */
 /* singular values and the left singular vectors are requested */
 if (conda) {
 /* Computing MAX */
 i__1 = max(lwqp3,lwcon); i__1 = max(i__1,lwsvd); // , expr subst  
 minwrk = *n + max(i__1,lwunq);
 }
 else {
 /* Computing MAX */
 i__1 = max(lwqp3,lwsvd);
 minwrk = *n + max(i__1,lwunq);
 }
 if (lquery) {
 if (rtrans) {
 cgesvd_("N", "O", n, n, &a[a_offset], lda, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, cdummy, &c_n1, rdummy, &ierr);
 }
 else {
 cgesvd_("O", "N", n, n, &a[a_offset], lda, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, cdummy, &c_n1, rdummy, &ierr);
 }
 lwrk_cgesvd__ = (integer) cdummy[0].r;
 if (conda) {
 /* Computing MAX */
 i__1 = max(lwrk_cgeqp3__,lwcon); i__1 = max(i__1, lwrk_cgesvd__); // , expr subst  
 optwrk = *n + max(i__1,lwrk_cunmqr__);
 }
 else {
 /* Computing MAX */
 i__1 = max(lwrk_cgeqp3__,lwrk_cgesvd__);
 optwrk = *n + max(i__1,lwrk_cunmqr__);
 }
 }
 }
 else if (rsvec && ! lsvec) {
 /* .. minimal and optimal sizes of the complex workspace if the */
 /* singular values and the right singular vectors are requested */
 if (conda) {
 /* Computing MAX */
 i__1 = max(lwqp3,lwcon);
 minwrk = *n + max(i__1,lwsvd);
 }
 else {
 minwrk = *n + max(lwqp3,lwsvd);
 }
 if (lquery) {
 if (rtrans) {
 cgesvd_("O", "N", n, n, &a[a_offset], lda, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, cdummy, &c_n1, rdummy, &ierr);
 }
 else {
 cgesvd_("N", "O", n, n, &a[a_offset], lda, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, cdummy, &c_n1, rdummy, &ierr);
 }
 lwrk_cgesvd__ = (integer) cdummy[0].r;
 if (conda) {
 /* Computing MAX */
 i__1 = max(lwrk_cgeqp3__,lwcon);
 optwrk = *n + max(i__1,lwrk_cgesvd__);
 }
 else {
 optwrk = *n + max(lwrk_cgeqp3__,lwrk_cgesvd__);
 }
 }
 }
 else {
 /* .. minimal and optimal sizes of the complex workspace if the */
 /* full SVD is requested */
 if (rtrans) {
 /* Computing MAX */
 i__1 = max(lwqp3,lwsvd);
 minwrk = max(i__1,lwunq);
 if (conda) {
 minwrk = max(minwrk,lwcon);
 }
 minwrk += *n;
 if (wntva) {
 /* .. minimal workspace length for N x N/2 CGEQRF */
 /* Computing MAX */
 i__1 = *n / 2;
 lwqrf = max(i__1,1);
 /* .. minimal workspace lengt for N/2 x N/2 CGESVD */
 /* Computing MAX */
 i__1 = *n / 2 * 3;
 lwsvd2 = max(i__1,1);
 lwunq2 = max(*n,1);
 /* Computing MAX */
 i__1 = lwqp3, i__2 = *n / 2 + lwqrf, i__1 = max(i__1,i__2) , i__2 = *n / 2 + lwsvd2, i__1 = max(i__1,i__2); i__2 = *n / 2 + lwunq2; i__1 = max(i__1,i__2); // ; expr subst  
 minwrk2 = max(i__1,lwunq);
 if (conda) {
 minwrk2 = max(minwrk2,lwcon);
 }
 minwrk2 = *n + minwrk2;
 minwrk = max(minwrk,minwrk2);
 }
 }
 else {
 /* Computing MAX */
 i__1 = max(lwqp3,lwsvd);
 minwrk = max(i__1,lwunq);
 if (conda) {
 minwrk = max(minwrk,lwcon);
 }
 minwrk += *n;
 if (wntva) {
 /* .. minimal workspace length for N/2 x N CGELQF */
 /* Computing MAX */
 i__1 = *n / 2;
 lwlqf = max(i__1,1);
 /* Computing MAX */
 i__1 = *n / 2 * 3;
 lwsvd2 = max(i__1,1);
 lwunlq = max(*n,1);
 /* Computing MAX */
 i__1 = lwqp3, i__2 = *n / 2 + lwlqf, i__1 = max(i__1,i__2) , i__2 = *n / 2 + lwsvd2, i__1 = max(i__1,i__2); i__2 = *n / 2 + lwunlq; i__1 = max(i__1,i__2); // ; expr subst  
 minwrk2 = max(i__1,lwunq);
 if (conda) {
 minwrk2 = max(minwrk2,lwcon);
 }
 minwrk2 = *n + minwrk2;
 minwrk = max(minwrk,minwrk2);
 }
 }
 if (lquery) {
 if (rtrans) {
 cgesvd_("O", "A", n, n, &a[a_offset], lda, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, cdummy, &c_n1, rdummy, &ierr);
 lwrk_cgesvd__ = (integer) cdummy[0].r;
 /* Computing MAX */
 i__1 = max(lwrk_cgeqp3__,lwrk_cgesvd__);
 optwrk = max(i__1,lwrk_cunmqr__);
 if (conda) {
 optwrk = max(optwrk,lwcon);
 }
 optwrk = *n + optwrk;
 if (wntva) {
 i__1 = *n / 2;
 cgeqrf_(n, &i__1, &u[u_offset], ldu, cdummy, cdummy, & c_n1, &ierr);
 lwrk_cgeqrf__ = (integer) cdummy[0].r;
 i__1 = *n / 2;
 i__2 = *n / 2;
 cgesvd_("S", "O", &i__1, &i__2, &v[v_offset], ldv, &s[ 1], &u[u_offset], ldu, &v[v_offset], ldv, cdummy, &c_n1, rdummy, &ierr);
 lwrk_cgesvd2__ = (integer) cdummy[0].r;
 i__1 = *n / 2;
 cunmqr_("R", "C", n, n, &i__1, &u[u_offset], ldu, cdummy, &v[v_offset], ldv, cdummy, &c_n1, & ierr);
 lwrk_cunmqr2__ = (integer) cdummy[0].r;
 /* Computing MAX */
 i__1 = lwrk_cgeqp3__, i__2 = *n / 2 + lwrk_cgeqrf__, i__1 = max(i__1,i__2), i__2 = *n / 2 + lwrk_cgesvd2__; i__1 = max(i__1,i__2); i__2 = *n / 2 + lwrk_cunmqr2__; // ; expr subst  
 optwrk2 = max(i__1,i__2);
 if (conda) {
 optwrk2 = max(optwrk2,lwcon);
 }
 optwrk2 = *n + optwrk2;
 optwrk = max(optwrk,optwrk2);
 }
 }
 else {
 cgesvd_("S", "O", n, n, &a[a_offset], lda, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, cdummy, &c_n1, rdummy, &ierr);
 lwrk_cgesvd__ = (integer) cdummy[0].r;
 /* Computing MAX */
 i__1 = max(lwrk_cgeqp3__,lwrk_cgesvd__);
 optwrk = max(i__1,lwrk_cunmqr__);
 if (conda) {
 optwrk = max(optwrk,lwcon);
 }
 optwrk = *n + optwrk;
 if (wntva) {
 i__1 = *n / 2;
 cgelqf_(&i__1, n, &u[u_offset], ldu, cdummy, cdummy, & c_n1, &ierr);
 lwrk_cgelqf__ = (integer) cdummy[0].r;
 i__1 = *n / 2;
 i__2 = *n / 2;
 cgesvd_("S", "O", &i__1, &i__2, &v[v_offset], ldv, &s[ 1], &u[u_offset], ldu, &v[v_offset], ldv, cdummy, &c_n1, rdummy, &ierr);
 lwrk_cgesvd2__ = (integer) cdummy[0].r;
 i__1 = *n / 2;
 cunmlq_("R", "N", n, n, &i__1, &u[u_offset], ldu, cdummy, &v[v_offset], ldv, cdummy, &c_n1, & ierr);
 lwrk_cunmlq__ = (integer) cdummy[0].r;
 /* Computing MAX */
 i__1 = lwrk_cgeqp3__, i__2 = *n / 2 + lwrk_cgelqf__, i__1 = max(i__1,i__2), i__2 = *n / 2 + lwrk_cgesvd2__; i__1 = max(i__1,i__2); i__2 = *n / 2 + lwrk_cunmlq__; // ; expr subst  
 optwrk2 = max(i__1,i__2);
 if (conda) {
 optwrk2 = max(optwrk2,lwcon);
 }
 optwrk2 = *n + optwrk2;
 optwrk = max(optwrk,optwrk2);
 }
 }
 }
 }
 minwrk = max(2,minwrk);
 optwrk = max(2,optwrk);
 if (*lcwork < minwrk && ! lquery) {
 *info = -19;
 }
 }
 if (*info == 0 && *lrwork < rminwrk && ! lquery) {
 *info = -21;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("CGESVDQ", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 else if (lquery) {
 /* Return optimal workspace */
 iwork[1] = iminwrk;
 cwork[1].r = (real) optwrk; cwork[1].i = 0.f; // , expr subst  
 cwork[2].r = (real) minwrk; cwork[2].i = 0.f; // , expr subst  
 rwork[1] = (real) rminwrk;
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 /* Quick return if the matrix is void. */
 if (*m == 0 || *n == 0) {
 /* .. all output is void. */
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 big = slamch_("O");
 ascaled = FALSE_;
 if (rowprm) {
 /* .. reordering the rows in decreasing sequence in the */
 /* ell-infinity norm - this enhances numerical robustness in */
 /* the case of differently scaled rows. */
 i__1 = *m;
 for (p = 1;
 p <= i__1;
 ++p) {
 /* RWORK(p) = ABS( A(p,ICAMAX(N,A(p,1),LDA)) ) */
 /* [[CLANGE will return NaN if an entry of the p-th row is Nan]] */
 rwork[p] = clange_("M", &c__1, n, &a[p + a_dim1], lda, rdummy);
 /* .. check for NaN's and Inf's */
 if (rwork[p] != rwork[p] || rwork[p] * 0.f != 0.f) {
 *info = -8;
 i__2 = -(*info);
 xerbla_("CGESVDQ", &i__2);
 return 0;
 }
 /* L1904: */
 }
 i__1 = *m - 1;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = *m - p + 1;
 q = isamax_(&i__2, &rwork[p], &c__1) + p - 1;
 iwork[*n + p] = q;
 if (p != q) {
 rtmp = rwork[p];
 rwork[p] = rwork[q];
 rwork[q] = rtmp;
 }
 /* L1952: */
 }
 if (rwork[1] == 0.f) {
 /* Quick return: A is the M x N zero matrix. */
 *numrank = 0;
 slaset_("G", n, &c__1, &c_b74, &c_b74, &s[1], n);
 if (wntus) {
 claset_("G", m, n, &c_b1, &c_b2, &u[u_offset], ldu) ;
 }
 if (wntua) {
 claset_("G", m, m, &c_b1, &c_b2, &u[u_offset], ldu) ;
 }
 if (wntva) {
 claset_("G", n, n, &c_b1, &c_b2, &v[v_offset], ldv) ;
 }
 if (wntuf) {
 claset_("G", n, &c__1, &c_b1, &c_b1, &cwork[1], n);
 claset_("G", m, n, &c_b1, &c_b2, &u[u_offset], ldu) ;
 }
 i__1 = *n;
 for (p = 1;
 p <= i__1;
 ++p) {
 iwork[p] = p;
 /* L5001: */
 }
 if (rowprm) {
 i__1 = *n + *m - 1;
 for (p = *n + 1;
 p <= i__1;
 ++p) {
 iwork[p] = p - *n;
 /* L5002: */
 }
 }
 if (conda) {
 rwork[1] = -1.f;
 }
 rwork[2] = -1.f;
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 if (rwork[1] > big / sqrt((real) (*m))) {
 /* .. to prevent overflow in the QR factorization, scale the */
 /* matrix by 1/sqrt(M) if too large entry detected */
 r__1 = sqrt((real) (*m));
 clascl_("G", &c__0, &c__0, &r__1, &c_b87, m, n, &a[a_offset], lda, &ierr);
 ascaled = TRUE_;
 }
 i__1 = *m - 1;
 claswp_(n, &a[a_offset], lda, &c__1, &i__1, &iwork[*n + 1], &c__1);
 }
 /* .. At this stage, preemptive scaling is done only to avoid column */
 /* norms overflows during the QR factorization. The SVD procedure should */
 /* have its own scaling to save the singular values from overflows and */
 /* underflows. That depends on the SVD procedure. */
 if (! rowprm) {
 rtmp = clange_("M", m, n, &a[a_offset], lda, &rwork[1]);
 if (rtmp != rtmp || rtmp * 0.f != 0.f) {
 *info = -8;
 i__1 = -(*info);
 xerbla_("CGESVDQ", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 if (rtmp > big / sqrt((real) (*m))) {
 /* .. to prevent overflow in the QR factorization, scale the */
 /* matrix by 1/sqrt(M) if too large entry detected */
 r__1 = sqrt((real) (*m));
 clascl_("G", &c__0, &c__0, &r__1, &c_b87, m, n, &a[a_offset], lda, &ierr);
 ascaled = TRUE_;
 }
 }
 /* .. QR factorization with column pivoting */
 /* A * P = Q * [ R ] */
 /* [ 0 ] */
 i__1 = *n;
 for (p = 1;
 p <= i__1;
 ++p) {
 /* .. all columns are free columns */
 iwork[p] = 0;
 /* L1963: */
 }
 i__1 = *lcwork - *n;
 cgeqp3_(m, n, &a[a_offset], lda, &iwork[1], &cwork[1], &cwork[*n + 1], & i__1, &rwork[1], &ierr);
 /* If the user requested accuracy level allows truncation in the */
 /* computed upper triangular factor, the matrix R is examined and, */
 /* if possible, replaced with its leading upper trapezoidal part. */
 epsln = slamch_("E");
 sfmin = slamch_("S");
 /* SMALL = SFMIN / EPSLN */
 nr = *n;
 if (accla) {
 /* Standard absolute error bound suffices. All sigma_i with */
 /* sigma_i < N*EPS*||A||_F are flushed to zero. This is an */
 /* aggressive enforcement of lower numerical rank by introducing a */
 /* backward error of the order of N*EPS*||A||_F. */
 nr = 1;
 rtmp = sqrt((real) (*n)) * epsln;
 i__1 = *n;
 for (p = 2;
 p <= i__1;
 ++p) {
 if (c_abs(&a[p + p * a_dim1]) < rtmp * c_abs(&a[a_dim1 + 1])) {
 goto L3002;
 }
 ++nr;
 /* L3001: */
 }
 L3002: ;
 }
 else if (acclm) {
 /* .. similarly as above, only slightly more gentle (less aggressive). */
 /* Sudden drop on the diagonal of R is used as the criterion for being */
 /* close-to-rank-deficient. The threshold is set to EPSLN=SLAMCH('E'). */
 /* [[This can be made more flexible by replacing this hard-coded value */
 /* with a user specified threshold.]] Also, the values that underflow */
 /* will be truncated. */
 nr = 1;
 i__1 = *n;
 for (p = 2;
 p <= i__1;
 ++p) {
 if (c_abs(&a[p + p * a_dim1]) < epsln * c_abs(&a[p - 1 + (p - 1) * a_dim1]) || c_abs(&a[p + p * a_dim1]) < sfmin) {
 goto L3402;
 }
 ++nr;
 /* L3401: */
 }
 L3402: ;
 }
 else {
 /* .. RRQR not authorized to determine numerical rank except in the */
 /* obvious case of zero pivots. */
 /* .. inspect R for exact zeros on the diagonal;
 */
 /* R(i,i)=0 => R(i:N,i:N)=0. */
 nr = 1;
 i__1 = *n;
 for (p = 2;
 p <= i__1;
 ++p) {
 if (c_abs(&a[p + p * a_dim1]) == 0.f) {
 goto L3502;
 }
 ++nr;
 /* L3501: */
 }
 L3502: if (conda) {
 /* Estimate the scaled condition number of A. Use the fact that it is */
 /* the same as the scaled condition number of R. */
 /* .. V is used as workspace */
 clacpy_("U", n, n, &a[a_offset], lda, &v[v_offset], ldv);
 /* Only the leading NR x NR submatrix of the triangular factor */
 /* is considered. Only if NR=N will this give a reliable error */
 /* bound. However, even for NR < N, this can be used on an */
 /* expert level and obtain useful information in the sense of */
 /* perturbation theory. */
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 rtmp = scnrm2_(&p, &v[p * v_dim1 + 1], &c__1);
 r__1 = 1.f / rtmp;
 csscal_(&p, &r__1, &v[p * v_dim1 + 1], &c__1);
 /* L3053: */
 }
 if (! (lsvec || rsvec)) {
 cpocon_("U", &nr, &v[v_offset], ldv, &c_b87, &rtmp, &cwork[1], &rwork[1], &ierr);
 }
 else {
 cpocon_("U", &nr, &v[v_offset], ldv, &c_b87, &rtmp, &cwork[*n + 1], &rwork[1], &ierr);
 }
 sconda = 1.f / sqrt(rtmp);
 /* For NR=N, SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1), */
 /* N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA */
 /* See the reference [1] for more details. */
 }
 }
 if (wntur) {
 n1 = nr;
 }
 else if (wntus || wntuf) {
 n1 = *n;
 }
 else if (wntua) {
 n1 = *m;
 }
 if (! (rsvec || lsvec)) {
 /* ....................................................................... */
 /* .. only the singular values are requested */
 /* ....................................................................... */
 if (rtrans) {
 /* .. compute the singular values of R**H = [A](1:NR,1:N)**H */
 /* .. set the lower triangle of [A] to [A](1:NR,1:N)**H and */
 /* the upper triangle of [A] to zero. */
 i__1 = min(*n,nr);
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = p + p * a_dim1;
 r_cnjg(&q__1, &a[p + p * a_dim1]);
 a[i__2].r = q__1.r; a[i__2].i = q__1.i; // , expr subst  
 i__2 = *n;
 for (q = p + 1;
 q <= i__2;
 ++q) {
 i__3 = q + p * a_dim1;
 r_cnjg(&q__1, &a[p + q * a_dim1]);
 a[i__3].r = q__1.r; a[i__3].i = q__1.i; // , expr subst  
 if (q <= nr) {
 i__3 = p + q * a_dim1;
 a[i__3].r = 0.f; a[i__3].i = 0.f; // , expr subst  
 }
 /* L1147: */
 }
 /* L1146: */
 }
 cgesvd_("N", "N", n, &nr, &a[a_offset], lda, &s[1], &u[u_offset], ldu, &v[v_offset], ldv, &cwork[1], lcwork, &rwork[1], info);
 }
 else {
 /* .. compute the singular values of R = [A](1:NR,1:N) */
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("L", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda);
 }
 cgesvd_("N", "N", &nr, n, &a[a_offset], lda, &s[1], &u[u_offset], ldu, &v[v_offset], ldv, &cwork[1], lcwork, &rwork[1], info);
 }
 }
 else if (lsvec && ! rsvec) {
 /* ....................................................................... */
 /* .. the singular values and the left singular vectors requested */
 /* ......................................................................."""""""" */
 if (rtrans) {
 /* .. apply CGESVD to R**H */
 /* .. copy R**H into [U] and overwrite [U] with the right singular */
 /* vectors of R */
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = *n;
 for (q = p;
 q <= i__2;
 ++q) {
 i__3 = q + p * u_dim1;
 r_cnjg(&q__1, &a[p + q * a_dim1]);
 u[i__3].r = q__1.r; u[i__3].i = q__1.i; // , expr subst  
 /* L1193: */
 }
 /* L1192: */
 }
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("U", &i__1, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) + 1] , ldu);
 }
 /* .. the left singular vectors not computed, the NR right singular */
 /* vectors overwrite [U](1:NR,1:NR) as conjugate transposed. These */
 /* will be pre-multiplied by Q to build the left singular vectors of A. */
 i__1 = *lcwork - *n;
 cgesvd_("N", "O", n, &nr, &u[u_offset], ldu, &s[1], &u[u_offset], ldu, &u[u_offset], ldu, &cwork[*n + 1], &i__1, &rwork[1], info);
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = p + p * u_dim1;
 r_cnjg(&q__1, &u[p + p * u_dim1]);
 u[i__2].r = q__1.r; u[i__2].i = q__1.i; // , expr subst  
 i__2 = nr;
 for (q = p + 1;
 q <= i__2;
 ++q) {
 r_cnjg(&q__1, &u[q + p * u_dim1]);
 ctmp.r = q__1.r; ctmp.i = q__1.i; // , expr subst  
 i__3 = q + p * u_dim1;
 r_cnjg(&q__1, &u[p + q * u_dim1]);
 u[i__3].r = q__1.r; u[i__3].i = q__1.i; // , expr subst  
 i__3 = p + q * u_dim1;
 u[i__3].r = ctmp.r; u[i__3].i = ctmp.i; // , expr subst  
 /* L1120: */
 }
 /* L1119: */
 }
 }
 else {
 /* .. apply CGESVD to R */
 /* .. copy R into [U] and overwrite [U] with the left singular vectors */
 clacpy_("U", &nr, n, &a[a_offset], lda, &u[u_offset], ldu);
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("L", &i__1, &i__2, &c_b1, &c_b1, &u[u_dim1 + 2], ldu);
 }
 /* .. the right singular vectors not computed, the NR left singular */
 /* vectors overwrite [U](1:NR,1:NR) */
 i__1 = *lcwork - *n;
 cgesvd_("O", "N", &nr, n, &u[u_offset], ldu, &s[1], &u[u_offset], ldu, &v[v_offset], ldv, &cwork[*n + 1], &i__1, &rwork[1], info);
 /* .. now [U](1:NR,1:NR) contains the NR left singular vectors of */
 /* R. These will be pre-multiplied by Q to build the left singular */
 /* vectors of A. */
 }
 /* .. assemble the left singular vector matrix U of dimensions */
 /* (M x NR) or (M x N) or (M x M). */
 if (nr < *m && ! wntuf) {
 i__1 = *m - nr;
 claset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1], ldu);
 if (nr < n1) {
 i__1 = n1 - nr;
 claset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * u_dim1 + 1], ldu);
 i__1 = *m - nr;
 i__2 = n1 - nr;
 claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + (nr + 1) * u_dim1], ldu);
 }
 }
 /* The Q matrix from the first QRF is built into the left singular */
 /* vectors matrix U. */
 if (! wntuf) {
 i__1 = *lcwork - *n;
 cunmqr_("L", "N", m, &n1, n, &a[a_offset], lda, &cwork[1], &u[ u_offset], ldu, &cwork[*n + 1], &i__1, &ierr);
 }
 if (rowprm && ! wntuf) {
 i__1 = *m - 1;
 claswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[*n + 1], & c_n1);
 }
 }
 else if (rsvec && ! lsvec) {
 /* ....................................................................... */
 /* .. the singular values and the right singular vectors requested */
 /* ....................................................................... */
 if (rtrans) {
 /* .. apply CGESVD to R**H */
 /* .. copy R**H into V and overwrite V with the left singular vectors */
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = *n;
 for (q = p;
 q <= i__2;
 ++q) {
 i__3 = q + p * v_dim1;
 r_cnjg(&q__1, &a[p + q * a_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 /* L1166: */
 }
 /* L1165: */
 }
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1] , ldv);
 }
 /* .. the left singular vectors of R**H overwrite V, the right singular */
 /* vectors not computed */
 if (wntvr || nr == *n) {
 i__1 = *lcwork - *n;
 cgesvd_("O", "N", n, &nr, &v[v_offset], ldv, &s[1], &u[ u_offset], ldu, &u[u_offset], ldu, &cwork[*n + 1], & i__1, &rwork[1], info);
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = p + p * v_dim1;
 r_cnjg(&q__1, &v[p + p * v_dim1]);
 v[i__2].r = q__1.r; v[i__2].i = q__1.i; // , expr subst  
 i__2 = nr;
 for (q = p + 1;
 q <= i__2;
 ++q) {
 r_cnjg(&q__1, &v[q + p * v_dim1]);
 ctmp.r = q__1.r; ctmp.i = q__1.i; // , expr subst  
 i__3 = q + p * v_dim1;
 r_cnjg(&q__1, &v[p + q * v_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 i__3 = p + q * v_dim1;
 v[i__3].r = ctmp.r; v[i__3].i = ctmp.i; // , expr subst  
 /* L1122: */
 }
 /* L1121: */
 }
 if (nr < *n) {
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = *n;
 for (q = nr + 1;
 q <= i__2;
 ++q) {
 i__3 = p + q * v_dim1;
 r_cnjg(&q__1, &v[q + p * v_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 /* L1104: */
 }
 /* L1103: */
 }
 }
 clapmt_(&c_false, &nr, n, &v[v_offset], ldv, &iwork[1]);
 }
 else {
 /* .. need all N right singular vectors and NR < N */
 /* [!] This is simple implementation that augments [V](1:N,1:NR) */
 /* by padding a zero block. In the case NR << N, a more efficient */
 /* way is to first use the QR factorization. For more details */
 /* how to implement this, see the " FULL SVD " branch. */
 i__1 = *n - nr;
 claset_("G", n, &i__1, &c_b1, &c_b1, &v[(nr + 1) * v_dim1 + 1] , ldv);
 i__1 = *lcwork - *n;
 cgesvd_("O", "N", n, n, &v[v_offset], ldv, &s[1], &u[u_offset] , ldu, &u[u_offset], ldu, &cwork[*n + 1], &i__1, & rwork[1], info);
 i__1 = *n;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = p + p * v_dim1;
 r_cnjg(&q__1, &v[p + p * v_dim1]);
 v[i__2].r = q__1.r; v[i__2].i = q__1.i; // , expr subst  
 i__2 = *n;
 for (q = p + 1;
 q <= i__2;
 ++q) {
 r_cnjg(&q__1, &v[q + p * v_dim1]);
 ctmp.r = q__1.r; ctmp.i = q__1.i; // , expr subst  
 i__3 = q + p * v_dim1;
 r_cnjg(&q__1, &v[p + q * v_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 i__3 = p + q * v_dim1;
 v[i__3].r = ctmp.r; v[i__3].i = ctmp.i; // , expr subst  
 /* L1124: */
 }
 /* L1123: */
 }
 clapmt_(&c_false, n, n, &v[v_offset], ldv, &iwork[1]);
 }
 }
 else {
 /* .. aply CGESVD to R */
 /* .. copy R into V and overwrite V with the right singular vectors */
 clacpy_("U", &nr, n, &a[a_offset], lda, &v[v_offset], ldv);
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("L", &i__1, &i__2, &c_b1, &c_b1, &v[v_dim1 + 2], ldv);
 }
 /* .. the right singular vectors overwrite V, the NR left singular */
 /* vectors stored in U(1:NR,1:NR) */
 if (wntvr || nr == *n) {
 i__1 = *lcwork - *n;
 cgesvd_("N", "O", &nr, n, &v[v_offset], ldv, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, &cwork[*n + 1], & i__1, &rwork[1], info);
 clapmt_(&c_false, &nr, n, &v[v_offset], ldv, &iwork[1]);
 /* .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H */
 }
 else {
 /* .. need all N right singular vectors and NR < N */
 /* [!] This is simple implementation that augments [V](1:NR,1:N) */
 /* by padding a zero block. In the case NR << N, a more efficient */
 /* way is to first use the LQ factorization. For more details */
 /* how to implement this, see the " FULL SVD " branch. */
 i__1 = *n - nr;
 claset_("G", &i__1, n, &c_b1, &c_b1, &v[nr + 1 + v_dim1], ldv);
 i__1 = *lcwork - *n;
 cgesvd_("N", "O", n, n, &v[v_offset], ldv, &s[1], &u[u_offset] , ldu, &v[v_offset], ldv, &cwork[*n + 1], &i__1, & rwork[1], info);
 clapmt_(&c_false, n, n, &v[v_offset], ldv, &iwork[1]);
 }
 /* .. now [V] contains the adjoint of the matrix of the right singular */
 /* vectors of A. */
 }
 }
 else {
 /* ....................................................................... */
 /* .. FULL SVD requested */
 /* ....................................................................... */
 if (rtrans) {
 /* .. apply CGESVD to R**H [[this option is left for R&D&T]] */
 if (wntvr || nr == *n) {
 /* .. copy R**H into [V] and overwrite [V] with the left singular */
 /* vectors of R**H */
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = *n;
 for (q = p;
 q <= i__2;
 ++q) {
 i__3 = q + p * v_dim1;
 r_cnjg(&q__1, &a[p + q * a_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 /* L1169: */
 }
 /* L1168: */
 }
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1], ldv);
 }
 /* .. the left singular vectors of R**H overwrite [V], the NR right */
 /* singular vectors of R**H stored in [U](1:NR,1:NR) as conjugate */
 /* transposed */
 i__1 = *lcwork - *n;
 cgesvd_("O", "A", n, &nr, &v[v_offset], ldv, &s[1], &v[ v_offset], ldv, &u[u_offset], ldu, &cwork[*n + 1], & i__1, &rwork[1], info);
 /* .. assemble V */
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = p + p * v_dim1;
 r_cnjg(&q__1, &v[p + p * v_dim1]);
 v[i__2].r = q__1.r; v[i__2].i = q__1.i; // , expr subst  
 i__2 = nr;
 for (q = p + 1;
 q <= i__2;
 ++q) {
 r_cnjg(&q__1, &v[q + p * v_dim1]);
 ctmp.r = q__1.r; ctmp.i = q__1.i; // , expr subst  
 i__3 = q + p * v_dim1;
 r_cnjg(&q__1, &v[p + q * v_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 i__3 = p + q * v_dim1;
 v[i__3].r = ctmp.r; v[i__3].i = ctmp.i; // , expr subst  
 /* L1116: */
 }
 /* L1115: */
 }
 if (nr < *n) {
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = *n;
 for (q = nr + 1;
 q <= i__2;
 ++q) {
 i__3 = p + q * v_dim1;
 r_cnjg(&q__1, &v[q + p * v_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 /* L1102: */
 }
 /* L1101: */
 }
 }
 clapmt_(&c_false, &nr, n, &v[v_offset], ldv, &iwork[1]);
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = p + p * u_dim1;
 r_cnjg(&q__1, &u[p + p * u_dim1]);
 u[i__2].r = q__1.r; u[i__2].i = q__1.i; // , expr subst  
 i__2 = nr;
 for (q = p + 1;
 q <= i__2;
 ++q) {
 r_cnjg(&q__1, &u[q + p * u_dim1]);
 ctmp.r = q__1.r; ctmp.i = q__1.i; // , expr subst  
 i__3 = q + p * u_dim1;
 r_cnjg(&q__1, &u[p + q * u_dim1]);
 u[i__3].r = q__1.r; u[i__3].i = q__1.i; // , expr subst  
 i__3 = p + q * u_dim1;
 u[i__3].r = ctmp.r; u[i__3].i = ctmp.i; // , expr subst  
 /* L1118: */
 }
 /* L1117: */
 }
 if (nr < *m && ! wntuf) {
 i__1 = *m - nr;
 claset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1] , ldu);
 if (nr < n1) {
 i__1 = n1 - nr;
 claset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * u_dim1 + 1], ldu);
 i__1 = *m - nr;
 i__2 = n1 - nr;
 claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + ( nr + 1) * u_dim1], ldu);
 }
 }
 }
 else {
 /* .. need all N right singular vectors and NR < N */
 /* .. copy R**H into [V] and overwrite [V] with the left singular */
 /* vectors of R**H */
 /* [[The optimal ratio N/NR for using QRF instead of padding */
 /* with zeros. Here hard coded to 2;
 it must be at least */
 /* two due to work space constraints.]] */
 /* OPTRATIO = ILAENV(6, 'CGESVD', 'S' // 'O', NR,N,0,0) */
 /* OPTRATIO = MAX( OPTRATIO, 2 ) */
 optratio = 2;
 if (optratio * nr > *n) {
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = *n;
 for (q = p;
 q <= i__2;
 ++q) {
 i__3 = q + p * v_dim1;
 r_cnjg(&q__1, &a[p + q * a_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 /* L1199: */
 }
 /* L1198: */
 }
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1], ldv);
 }
 i__1 = *n - nr;
 claset_("A", n, &i__1, &c_b1, &c_b1, &v[(nr + 1) * v_dim1 + 1], ldv);
 i__1 = *lcwork - *n;
 cgesvd_("O", "A", n, n, &v[v_offset], ldv, &s[1], &v[ v_offset], ldv, &u[u_offset], ldu, &cwork[*n + 1], &i__1, &rwork[1], info);
 i__1 = *n;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = p + p * v_dim1;
 r_cnjg(&q__1, &v[p + p * v_dim1]);
 v[i__2].r = q__1.r; v[i__2].i = q__1.i; // , expr subst  
 i__2 = *n;
 for (q = p + 1;
 q <= i__2;
 ++q) {
 r_cnjg(&q__1, &v[q + p * v_dim1]);
 ctmp.r = q__1.r; ctmp.i = q__1.i; // , expr subst  
 i__3 = q + p * v_dim1;
 r_cnjg(&q__1, &v[p + q * v_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 i__3 = p + q * v_dim1;
 v[i__3].r = ctmp.r; v[i__3].i = ctmp.i; // , expr subst  
 /* L1114: */
 }
 /* L1113: */
 }
 clapmt_(&c_false, n, n, &v[v_offset], ldv, &iwork[1]);
 /* .. assemble the left singular vector matrix U of dimensions */
 /* (M x N1), i.e. (M x N) or (M x M). */
 i__1 = *n;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = p + p * u_dim1;
 r_cnjg(&q__1, &u[p + p * u_dim1]);
 u[i__2].r = q__1.r; u[i__2].i = q__1.i; // , expr subst  
 i__2 = *n;
 for (q = p + 1;
 q <= i__2;
 ++q) {
 r_cnjg(&q__1, &u[q + p * u_dim1]);
 ctmp.r = q__1.r; ctmp.i = q__1.i; // , expr subst  
 i__3 = q + p * u_dim1;
 r_cnjg(&q__1, &u[p + q * u_dim1]);
 u[i__3].r = q__1.r; u[i__3].i = q__1.i; // , expr subst  
 i__3 = p + q * u_dim1;
 u[i__3].r = ctmp.r; u[i__3].i = ctmp.i; // , expr subst  
 /* L1112: */
 }
 /* L1111: */
 }
 if (*n < *m && ! wntuf) {
 i__1 = *m - *n;
 claset_("A", &i__1, n, &c_b1, &c_b1, &u[*n + 1 + u_dim1], ldu);
 if (*n < n1) {
 i__1 = n1 - *n;
 claset_("A", n, &i__1, &c_b1, &c_b1, &u[(*n + 1) * u_dim1 + 1], ldu);
 i__1 = *m - *n;
 i__2 = n1 - *n;
 claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[*n + 1 + (*n + 1) * u_dim1], ldu);
 }
 }
 }
 else {
 /* .. copy R**H into [U] and overwrite [U] with the right */
 /* singular vectors of R */
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = *n;
 for (q = p;
 q <= i__2;
 ++q) {
 i__3 = q + (nr + p) * u_dim1;
 r_cnjg(&q__1, &a[p + q * a_dim1]);
 u[i__3].r = q__1.r; u[i__3].i = q__1.i; // , expr subst  
 /* L1197: */
 }
 /* L1196: */
 }
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("U", &i__1, &i__2, &c_b1, &c_b1, &u[(nr + 2) * u_dim1 + 1], ldu);
 }
 i__1 = *lcwork - *n - nr;
 cgeqrf_(n, &nr, &u[(nr + 1) * u_dim1 + 1], ldu, &cwork[*n + 1], &cwork[*n + nr + 1], &i__1, &ierr);
 i__1 = nr;
 for (p = 1;
 p <= i__1;
 ++p) {
 i__2 = *n;
 for (q = 1;
 q <= i__2;
 ++q) {
 i__3 = q + p * v_dim1;
 r_cnjg(&q__1, &u[p + (nr + q) * u_dim1]);
 v[i__3].r = q__1.r; v[i__3].i = q__1.i; // , expr subst  
 /* L1144: */
 }
 /* L1143: */
 }
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1], ldv);
 i__1 = *lcwork - *n - nr;
 cgesvd_("S", "O", &nr, &nr, &v[v_offset], ldv, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, &cwork[*n + nr + 1], &i__1, &rwork[1], info);
 i__1 = *n - nr;
 claset_("A", &i__1, &nr, &c_b1, &c_b1, &v[nr + 1 + v_dim1] , ldv);
 i__1 = *n - nr;
 claset_("A", &nr, &i__1, &c_b1, &c_b1, &v[(nr + 1) * v_dim1 + 1], ldv);
 i__1 = *n - nr;
 i__2 = *n - nr;
 claset_("A", &i__1, &i__2, &c_b1, &c_b2, &v[nr + 1 + (nr + 1) * v_dim1], ldv);
 i__1 = *lcwork - *n - nr;
 cunmqr_("R", "C", n, n, &nr, &u[(nr + 1) * u_dim1 + 1], ldu, &cwork[*n + 1], &v[v_offset], ldv, &cwork[*n + nr + 1], &i__1, &ierr);
 clapmt_(&c_false, n, n, &v[v_offset], ldv, &iwork[1]);
 /* .. assemble the left singular vector matrix U of dimensions */
 /* (M x NR) or (M x N) or (M x M). */
 if (nr < *m && ! wntuf) {
 i__1 = *m - nr;
 claset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1], ldu);
 if (nr < n1) {
 i__1 = n1 - nr;
 claset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * u_dim1 + 1], ldu);
 i__1 = *m - nr;
 i__2 = n1 - nr;
 claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + (nr + 1) * u_dim1], ldu);
 }
 }
 }
 }
 }
 else {
 /* .. apply CGESVD to R [[this is the recommended option]] */
 if (wntvr || nr == *n) {
 /* .. copy R into [V] and overwrite V with the right singular vectors */
 clacpy_("U", &nr, n, &a[a_offset], lda, &v[v_offset], ldv);
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("L", &i__1, &i__2, &c_b1, &c_b1, &v[v_dim1 + 2], ldv);
 }
 /* .. the right singular vectors of R overwrite [V], the NR left */
 /* singular vectors of R stored in [U](1:NR,1:NR) */
 i__1 = *lcwork - *n;
 cgesvd_("S", "O", &nr, n, &v[v_offset], ldv, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, &cwork[*n + 1], & i__1, &rwork[1], info);
 clapmt_(&c_false, &nr, n, &v[v_offset], ldv, &iwork[1]);
 /* .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H */
 /* .. assemble the left singular vector matrix U of dimensions */
 /* (M x NR) or (M x N) or (M x M). */
 if (nr < *m && ! wntuf) {
 i__1 = *m - nr;
 claset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1] , ldu);
 if (nr < n1) {
 i__1 = n1 - nr;
 claset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * u_dim1 + 1], ldu);
 i__1 = *m - nr;
 i__2 = n1 - nr;
 claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + ( nr + 1) * u_dim1], ldu);
 }
 }
 }
 else {
 /* .. need all N right singular vectors and NR < N */
 /* .. the requested number of the left singular vectors */
 /* is then N1 (N or M) */
 /* [[The optimal ratio N/NR for using LQ instead of padding */
 /* with zeros. Here hard coded to 2;
 it must be at least */
 /* two due to work space constraints.]] */
 /* OPTRATIO = ILAENV(6, 'CGESVD', 'S' // 'O', NR,N,0,0) */
 /* OPTRATIO = MAX( OPTRATIO, 2 ) */
 optratio = 2;
 if (optratio * nr > *n) {
 clacpy_("U", &nr, n, &a[a_offset], lda, &v[v_offset], ldv);
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("L", &i__1, &i__2, &c_b1, &c_b1, &v[v_dim1 + 2], ldv);
 }
 /* .. the right singular vectors of R overwrite [V], the NR left */
 /* singular vectors of R stored in [U](1:NR,1:NR) */
 i__1 = *n - nr;
 claset_("A", &i__1, n, &c_b1, &c_b1, &v[nr + 1 + v_dim1], ldv);
 i__1 = *lcwork - *n;
 cgesvd_("S", "O", n, n, &v[v_offset], ldv, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, &cwork[*n + 1], &i__1, &rwork[1], info);
 clapmt_(&c_false, n, n, &v[v_offset], ldv, &iwork[1]);
 /* .. now [V] contains the adjoint of the matrix of the right */
 /* singular vectors of A. The leading N left singular vectors */
 /* are in [U](1:N,1:N) */
 /* .. assemble the left singular vector matrix U of dimensions */
 /* (M x N1), i.e. (M x N) or (M x M). */
 if (*n < *m && ! wntuf) {
 i__1 = *m - *n;
 claset_("A", &i__1, n, &c_b1, &c_b1, &u[*n + 1 + u_dim1], ldu);
 if (*n < n1) {
 i__1 = n1 - *n;
 claset_("A", n, &i__1, &c_b1, &c_b1, &u[(*n + 1) * u_dim1 + 1], ldu);
 i__1 = *m - *n;
 i__2 = n1 - *n;
 claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[*n + 1 + (*n + 1) * u_dim1], ldu);
 }
 }
 }
 else {
 clacpy_("U", &nr, n, &a[a_offset], lda, &u[nr + 1 + u_dim1], ldu);
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("L", &i__1, &i__2, &c_b1, &c_b1, &u[nr + 2 + u_dim1], ldu);
 }
 i__1 = *lcwork - *n - nr;
 cgelqf_(&nr, n, &u[nr + 1 + u_dim1], ldu, &cwork[*n + 1], &cwork[*n + nr + 1], &i__1, &ierr);
 clacpy_("L", &nr, &nr, &u[nr + 1 + u_dim1], ldu, &v[ v_offset], ldv);
 if (nr > 1) {
 i__1 = nr - 1;
 i__2 = nr - 1;
 claset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1], ldv);
 }
 i__1 = *lcwork - *n - nr;
 cgesvd_("S", "O", &nr, &nr, &v[v_offset], ldv, &s[1], &u[ u_offset], ldu, &v[v_offset], ldv, &cwork[*n + nr + 1], &i__1, &rwork[1], info);
 i__1 = *n - nr;
 claset_("A", &i__1, &nr, &c_b1, &c_b1, &v[nr + 1 + v_dim1] , ldv);
 i__1 = *n - nr;
 claset_("A", &nr, &i__1, &c_b1, &c_b1, &v[(nr + 1) * v_dim1 + 1], ldv);
 i__1 = *n - nr;
 i__2 = *n - nr;
 claset_("A", &i__1, &i__2, &c_b1, &c_b2, &v[nr + 1 + (nr + 1) * v_dim1], ldv);
 i__1 = *lcwork - *n - nr;
 cunmlq_("R", "N", n, n, &nr, &u[nr + 1 + u_dim1], ldu, & cwork[*n + 1], &v[v_offset], ldv, &cwork[*n + nr + 1], &i__1, &ierr);
 clapmt_(&c_false, n, n, &v[v_offset], ldv, &iwork[1]);
 /* .. assemble the left singular vector matrix U of dimensions */
 /* (M x NR) or (M x N) or (M x M). */
 if (nr < *m && ! wntuf) {
 i__1 = *m - nr;
 claset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1], ldu);
 if (nr < n1) {
 i__1 = n1 - nr;
 claset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * u_dim1 + 1], ldu);
 i__1 = *m - nr;
 i__2 = n1 - nr;
 claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + (nr + 1) * u_dim1], ldu);
 }
 }
 }
 }
 /* .. end of the "R**H or R" branch */
 }
 /* The Q matrix from the first QRF is built into the left singular */
 /* vectors matrix U. */
 if (! wntuf) {
 i__1 = *lcwork - *n;
 cunmqr_("L", "N", m, &n1, n, &a[a_offset], lda, &cwork[1], &u[ u_offset], ldu, &cwork[*n + 1], &i__1, &ierr);
 }
 if (rowprm && ! wntuf) {
 i__1 = *m - 1;
 claswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[*n + 1], & c_n1);
 }
 /* ... end of the "full SVD" branch */
 }
 /* Check whether some singular values are returned as zeros, e.g. */
 /* due to underflow, and update the numerical rank. */
 p = nr;
 for (q = p;
 q >= 1;
 --q) {
 if (s[q] > 0.f) {
 goto L4002;
 }
 --nr;
 /* L4001: */
 }
 L4002: /* .. if numerical rank deficiency is detected, the truncated */
 /* singular values are set to zero. */
 if (nr < *n) {
 i__1 = *n - nr;
 slaset_("G", &i__1, &c__1, &c_b74, &c_b74, &s[nr + 1], n);
 }
 /* .. undo scaling;
 this may cause overflow in the largest singular */
 /* values. */
 if (ascaled) {
 r__1 = sqrt((real) (*m));
 slascl_("G", &c__0, &c__0, &c_b87, &r__1, &nr, &c__1, &s[1], n, &ierr);
 }
 if (conda) {
 rwork[1] = sconda;
 }
 rwork[2] = (real) (p - nr);
 /* .. p-NR is the number of singular values that are computed as */
 /* exact zeros in CGESVD() applied to the (possibly truncated) */
 /* full row rank triangular (trapezoidal) factor of A. */
 *numrank = nr;
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of CGESVDQ */
 }
 /* cgesvdq_ */
 