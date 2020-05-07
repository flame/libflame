/* ../netlib/v3.9.0/sbdsvdx.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static real c_b10 = 1.f;
 static doublereal c_b14 = -.125;
 static integer c__1 = 1;
 static real c_b20 = 0.f;
 static integer c__2 = 2;
 /* > \brief \b SBDSVDX */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download SBDSVDX + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sbdsvdx .f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sbdsvdx .f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sbdsvdx .f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE SBDSVDX( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU, */
 /* $ NS, S, Z, LDZ, WORK, IWORK, INFO ) */
 /* .. Scalar Arguments .. */
 /* CHARACTER JOBZ, RANGE, UPLO */
 /* INTEGER IL, INFO, IU, LDZ, N, NS */
 /* REAL VL, VU */
 /* .. */
 /* .. Array Arguments .. */
 /* INTEGER IWORK( * ) */
 /* REAL D( * ), E( * ), S( * ), WORK( * ), */
 /* Z( LDZ, * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > SBDSVDX computes the singular value decomposition (SVD) of a real */
 /* > N-by-N (upper or lower) bidiagonal matrix B, B = U * S * VT, */
 /* > where S is a diagonal matrix with non-negative diagonal elements */
 /* > (the singular values of B), and U and VT are orthogonal matrices */
 /* > of left and right singular vectors, respectively. */
 /* > */
 /* > Given an upper bidiagonal B with diagonal D = [ d_1 d_2 ... d_N ] */
 /* > and superdiagonal E = [ e_1 e_2 ... e_N-1 ], SBDSVDX computes the */
 /* > singular value decompositon of B through the eigenvalues and */
 /* > eigenvectors of the N*2-by-N*2 tridiagonal matrix */
 /* > */
 /* > | 0 d_1 | */
 /* > | d_1 0 e_1 | */
 /* > TGK = | e_1 0 d_2 | */
 /* > | d_2 . . | */
 /* > | . . . | */
 /* > */
 /* > If (s,u,v) is a singular triplet of B with ||u|| = ||v|| = 1, then */
 /* > (+/-s,q), ||q|| = 1, are eigenpairs of TGK, with q = P * ( u' +/-v' ) / */
 /* > sqrt(2) = ( v_1 u_1 v_2 u_2 ... v_n u_n ) / sqrt(2), and */
 /* > P = [ e_{
n+1}
 e_{
1}
 e_{
n+2}
 e_{
2}
 ... ]. */
 /* > */
 /* > Given a TGK matrix, one can either a) compute -s,-v and change signs */
 /* > so that the singular values (and corresponding vectors) are already in */
 /* > descending order (as in SGESVD/SGESDD) or b) compute s,v and reorder */
 /* > the values (and corresponding vectors). SBDSVDX implements a) by */
 /* > calling SSTEVX (bisection plus inverse iteration, to be replaced */
 /* > with a version of the Multiple Relative Robust Representation */
 /* > algorithm. (See P. Willems and B. Lang, A framework for the MR^3 */
 /* > algorithm: theory and implementation, SIAM J. Sci. Comput., */
 /* > 35:740-766, 2013.) */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] UPLO */
 /* > \verbatim */
 /* > UPLO is CHARACTER*1 */
 /* > = 'U': B is upper bidiagonal;
 */
 /* > = 'L': B is lower bidiagonal. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] JOBZ */
 /* > \verbatim */
 /* > JOBZ is CHARACTER*1 */
 /* > = 'N': Compute singular values only;
 */
 /* > = 'V': Compute singular values and singular vectors. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] RANGE */
 /* > \verbatim */
 /* > RANGE is CHARACTER*1 */
 /* > = 'A': all singular values will be found. */
 /* > = 'V': all singular values in the half-open interval [VL,VU) */
 /* > will be found. */
 /* > = 'I': the IL-th through IU-th singular values will be found. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The order of the bidiagonal matrix. N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] D */
 /* > \verbatim */
 /* > D is REAL array, dimension (N) */
 /* > The n diagonal elements of the bidiagonal matrix B. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] E */
 /* > \verbatim */
 /* > E is REAL array, dimension (max(1,N-1)) */
 /* > The (n-1) superdiagonal elements of the bidiagonal matrix */
 /* > B in elements 1 to N-1. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] VL */
 /* > \verbatim */
 /* > VL is REAL */
 /* > If RANGE='V', the lower bound of the interval to */
 /* > be searched for singular values. VU > VL. */
 /* > Not referenced if RANGE = 'A' or 'I'. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] VU */
 /* > \verbatim */
 /* > VU is REAL */
 /* > If RANGE='V', the upper bound of the interval to */
 /* > be searched for singular values. VU > VL. */
 /* > Not referenced if RANGE = 'A' or 'I'. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] IL */
 /* > \verbatim */
 /* > IL is INTEGER */
 /* > If RANGE='I', the index of the */
 /* > smallest singular value to be returned. */
 /* > 1 <= IL <= IU <= min(M,N), if min(M,N) > 0. */
 /* > Not referenced if RANGE = 'A' or 'V'. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] IU */
 /* > \verbatim */
 /* > IU is INTEGER */
 /* > If RANGE='I', the index of the */
 /* > largest singular value to be returned. */
 /* > 1 <= IL <= IU <= min(M,N), if min(M,N) > 0. */
 /* > Not referenced if RANGE = 'A' or 'V'. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] NS */
 /* > \verbatim */
 /* > NS is INTEGER */
 /* > The total number of singular values found. 0 <= NS <= N. */
 /* > If RANGE = 'A', NS = N, and if RANGE = 'I', NS = IU-IL+1. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] S */
 /* > \verbatim */
 /* > S is REAL array, dimension (N) */
 /* > The first NS elements contain the selected singular values in */
 /* > ascending order. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] Z */
 /* > \verbatim */
 /* > Z is REAL array, dimension (2*N,K) */
 /* > If JOBZ = 'V', then if INFO = 0 the first NS columns of Z */
 /* > contain the singular vectors of the matrix B corresponding to */
 /* > the selected singular values, with U in rows 1 to N and V */
 /* > in rows N+1 to N*2, i.e. */
 /* > Z = [ U ] */
 /* > [ V ] */
 /* > If JOBZ = 'N', then Z is not referenced. */
 /* > Note: The user must ensure that at least K = NS+1 columns are */
 /* > supplied in the array Z;
 if RANGE = 'V', the exact value of */
 /* > NS is not known in advance and an upper bound must be used. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDZ */
 /* > \verbatim */
 /* > LDZ is INTEGER */
 /* > The leading dimension of the array Z. LDZ >= 1, and if */
 /* > JOBZ = 'V', LDZ >= max(2,N*2). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is REAL array, dimension (14*N) */
 /* > \endverbatim */
 /* > */
 /* > \param[out] IWORK */
 /* > \verbatim */
 /* > IWORK is INTEGER array, dimension (12*N) */
 /* > If JOBZ = 'V', then if INFO = 0, the first NS elements of */
 /* > IWORK are zero. If INFO > 0, then IWORK contains the indices */
 /* > of the eigenvectors that failed to converge in DSTEVX. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value */
 /* > > 0: if INFO = i, then i eigenvectors failed to converge */
 /* > in SSTEVX. The indices of the eigenvectors */
 /* > (as returned by SSTEVX) are stored in the */
 /* > array IWORK. */
 /* > if INFO = N*2 + 1, an internal error occurred. */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date June 2016 */
 /* > \ingroup realOTHEReigen */
 /* ===================================================================== */
 /* Subroutine */
 int sbdsvdx_(char *uplo, char *jobz, char *range, integer *n, real *d__, real *e, real *vl, real *vu, integer *il, integer *iu, integer *ns, real *s, real *z__, integer *ldz, real *work, integer * iwork, integer *info) {
 /* System generated locals */
 integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
 real r__1, r__2, r__3, r__4;
 doublereal d__1;
 /* Builtin functions */
 double r_sign(real *, real *), sqrt(doublereal), pow_dd(doublereal *, doublereal *);
 /* Local variables */
 integer i__, j, k;
 real d1;
 integer j1, j2;
 real mu, eps;
 integer nsl;
 real tol, ulp;
 integer nru, nrv;
 real emin;
 integer ntgk;
 real smin, smax;
 extern real sdot_(integer *, real *, integer *, real *, integer *);
 real nrmu, nrmv;
 logical sveq0;
 extern real snrm2_(integer *, real *, integer *);
 integer idbeg;
 real sqrt2;
 integer idend, isbeg;
 extern logical lsame_(char *, char *);
 integer idtgk, ietgk;
 extern /* Subroutine */
 int sscal_(integer *, real *, real *, integer *);
 integer iltgk, itemp, icolz;
 logical allsv;
 integer idptr;
 logical indsv;
 integer ieptr, iutgk;
 real vltgk;
 logical lower;
 real zjtji;
 logical split, valsv;
 integer isplt;
 real ortol, vutgk;
 extern /* Subroutine */
 int scopy_(integer *, real *, integer *, real *, integer *), sswap_(integer *, real *, integer *, real *, integer * );
 logical wantz;
 char rngvx[1];
 integer irowu, irowv;
 extern /* Subroutine */
 int saxpy_(integer *, real *, real *, integer *, real *, integer *);
 integer irowz, iifail;
 extern real slamch_(char *);
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 extern integer isamax_(integer *, real *, integer *);
 real abstol;
 extern /* Subroutine */
 int slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
 real thresh;
 integer iiwork;
 extern /* Subroutine */
 int sstevx_(char *, char *, integer *, real *, real *, real *, real *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *, integer * );
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
 --d__;
 --e;
 --s;
 z_dim1 = *ldz;
 z_offset = 1 + z_dim1;
 z__ -= z_offset;
 --work;
 --iwork;
 /* Function Body */
 allsv = lsame_(range, "A");
 valsv = lsame_(range, "V");
 indsv = lsame_(range, "I");
 wantz = lsame_(jobz, "V");
 lower = lsame_(uplo, "L");
 *info = 0;
 if (! lsame_(uplo, "U") && ! lower) {
 *info = -1;
 }
 else if (! (wantz || lsame_(jobz, "N"))) {
 *info = -2;
 }
 else if (! (allsv || valsv || indsv)) {
 *info = -3;
 }
 else if (*n < 0) {
 *info = -4;
 }
 else if (*n > 0) {
 if (valsv) {
 if (*vl < 0.f) {
 *info = -7;
 }
 else if (*vu <= *vl) {
 *info = -8;
 }
 }
 else if (indsv) {
 if (*il < 1 || *il > max(1,*n)) {
 *info = -9;
 }
 else if (*iu < min(*n,*il) || *iu > *n) {
 *info = -10;
 }
 }
 }
 if (*info == 0) {
 if (*ldz < 1 || wantz && *ldz < *n << 1) {
 *info = -14;
 }
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("SBDSVDX", &i__1);
 return 0;
 }
 /* Quick return if possible (N.LE.1) */
 *ns = 0;
 if (*n == 0) {
 return 0;
 }
 if (*n == 1) {
 if (allsv || indsv) {
 *ns = 1;
 s[1] = abs(d__[1]);
 }
 else {
 if (*vl < abs(d__[1]) && *vu >= abs(d__[1])) {
 *ns = 1;
 s[1] = abs(d__[1]);
 }
 }
 if (wantz) {
 z__[z_dim1 + 1] = r_sign(&c_b10, &d__[1]);
 z__[z_dim1 + 2] = 1.f;
 }
 return 0;
 }
 abstol = slamch_("Safe Minimum") * 2;
 ulp = slamch_("Precision");
 eps = slamch_("Epsilon");
 sqrt2 = sqrt(2.f);
 ortol = sqrt(ulp);
 /* Criterion for splitting is taken from SBDSQR when singular */
 /* values are computed to relative accuracy TOL. (See J. Demmel and */
 /* W. Kahan, Accurate singular values of bidiagonal matrices, SIAM */
 /* J. Sci. and Stat. Comput., 11:873–912, 1990.) */
 /* Computing MAX */
 /* Computing MIN */
 d__1 = (doublereal) eps;
 r__3 = 100.f; r__4 = pow_dd(&d__1, &c_b14); // , expr subst  
 r__1 = 10.f; r__2 = min(r__3,r__4); // , expr subst  
 tol = max(r__1,r__2) * eps;
 /* Compute approximate maximum, minimum singular values. */
 i__ = isamax_(n, &d__[1], &c__1);
 smax = (r__1 = d__[i__], abs(r__1));
 i__1 = *n - 1;
 i__ = isamax_(&i__1, &e[1], &c__1);
 /* Computing MAX */
 r__2 = smax; r__3 = (r__1 = e[i__], abs(r__1)); // , expr subst  
 smax = max(r__2,r__3);
 /* Compute threshold for neglecting D's and E's. */
 smin = abs(d__[1]);
 if (smin != 0.f) {
 mu = smin;
 i__1 = *n;
 for (i__ = 2;
 i__ <= i__1;
 ++i__) {
 mu = (r__2 = d__[i__], abs(r__2)) * (mu / (mu + (r__1 = e[i__ - 1] , abs(r__1))));
 smin = min(smin,mu);
 if (smin == 0.f) {
 goto L2;
 }
 }
 L2: ;
 }
 smin /= sqrt((real) (*n));
 thresh = tol * smin;
 /* Check for zeros in D and E (splits), i.e. submatrices. */
 i__1 = *n - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 if ((r__1 = d__[i__], abs(r__1)) <= thresh) {
 d__[i__] = 0.f;
 }
 if ((r__1 = e[i__], abs(r__1)) <= thresh) {
 e[i__] = 0.f;
 }
 }
 if ((r__1 = d__[*n], abs(r__1)) <= thresh) {
 d__[*n] = 0.f;
 }
 /* Pointers for arrays used by SSTEVX. */
 idtgk = 1;
 ietgk = idtgk + (*n << 1);
 itemp = ietgk + (*n << 1);
 iifail = 1;
 iiwork = iifail + (*n << 1);
 /* Set RNGVX, which corresponds to RANGE for SSTEVX in TGK mode. */
 /* VL,VU or IL,IU are redefined to conform to implementation a) */
 /* described in the leading comments. */
 iltgk = 0;
 iutgk = 0;
 vltgk = 0.f;
 vutgk = 0.f;
 if (allsv) {
 /* All singular values will be found. We aim at -s (see */
 /* leading comments) with RNGVX = 'I'. IL and IU are set */
 /* later (as ILTGK and IUTGK) according to the dimension */
 /* of the active submatrix. */
 *(unsigned char *)rngvx = 'I';
 if (wantz) {
 i__1 = *n << 1;
 i__2 = *n + 1;
 slaset_("F", &i__1, &i__2, &c_b20, &c_b20, &z__[z_offset], ldz);
 }
 }
 else if (valsv) {
 /* Find singular values in a half-open interval. We aim */
 /* at -s (see leading comments) and we swap VL and VU */
 /* (as VUTGK and VLTGK), changing their signs. */
 *(unsigned char *)rngvx = 'V';
 vltgk = -(*vu);
 vutgk = -(*vl);
 /* WORK( IDTGK:IDTGK+2*N-1 ) = ZERO */
 i__1 = *n << 1;
 for (j1 = 1;
 j1 <= i__1;
 ++j1) {
 work[idtgk - 1 + j1] = 0.f;
 }
 scopy_(n, &d__[1], &c__1, &work[ietgk], &c__2);
 i__1 = *n - 1;
 scopy_(&i__1, &e[1], &c__1, &work[ietgk + 1], &c__2);
 i__1 = *n << 1;
 sstevx_("N", "V", &i__1, &work[idtgk], &work[ietgk], &vltgk, &vutgk, & iltgk, &iltgk, &abstol, ns, &s[1], &z__[z_offset], ldz, &work[ itemp], &iwork[iiwork], &iwork[iifail], info);
 if (*ns == 0) {
 return 0;
 }
 else {
 if (wantz) {
 i__1 = *n << 1;
 slaset_("F", &i__1, ns, &c_b20, &c_b20, &z__[z_offset], ldz);
 }
 }
 }
 else if (indsv) {
 /* Find the IL-th through the IU-th singular values. We aim */
 /* at -s (see leading comments) and indices are mapped into */
 /* values, therefore mimicking SSTEBZ, where */
 /* GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN */
 /* GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN */
 iltgk = *il;
 iutgk = *iu;
 *(unsigned char *)rngvx = 'V';
 /* WORK( IDTGK:IDTGK+2*N-1 ) = ZERO */
 i__1 = *n << 1;
 for (j1 = 1;
 j1 <= i__1;
 ++j1) {
 work[idtgk - 1 + j1] = 0.f;
 }
 scopy_(n, &d__[1], &c__1, &work[ietgk], &c__2);
 i__1 = *n - 1;
 scopy_(&i__1, &e[1], &c__1, &work[ietgk + 1], &c__2);
 i__1 = *n << 1;
 sstevx_("N", "I", &i__1, &work[idtgk], &work[ietgk], &vltgk, &vltgk, & iltgk, &iltgk, &abstol, ns, &s[1], &z__[z_offset], ldz, &work[ itemp], &iwork[iiwork], &iwork[iifail], info);
 vltgk = s[1] - smax * 2.f * ulp * *n;
 /* WORK( IDTGK:IDTGK+2*N-1 ) = ZERO */
 i__1 = *n << 1;
 for (j1 = 1;
 j1 <= i__1;
 ++j1) {
 work[idtgk - 1 + j1] = 0.f;
 }
 scopy_(n, &d__[1], &c__1, &work[ietgk], &c__2);
 i__1 = *n - 1;
 scopy_(&i__1, &e[1], &c__1, &work[ietgk + 1], &c__2);
 i__1 = *n << 1;
 sstevx_("N", "I", &i__1, &work[idtgk], &work[ietgk], &vutgk, &vutgk, & iutgk, &iutgk, &abstol, ns, &s[1], &z__[z_offset], ldz, &work[ itemp], &iwork[iiwork], &iwork[iifail], info);
 vutgk = s[1] + smax * 2.f * ulp * *n;
 vutgk = min(vutgk,0.f);
 /* If VLTGK=VUTGK, SSTEVX returns an error message, */
 /* so if needed we change VUTGK slightly. */
 if (vltgk == vutgk) {
 vltgk -= tol;
 }
 if (wantz) {
 i__1 = *n << 1;
 i__2 = *iu - *il + 1;
 slaset_("F", &i__1, &i__2, &c_b20, &c_b20, &z__[z_offset], ldz);
 }
 }
 /* Initialize variables and pointers for S, Z, and WORK. */
 /* NRU, NRV: number of rows in U and V for the active submatrix */
 /* IDBEG, ISBEG: offsets for the entries of D and S */
 /* IROWZ, ICOLZ: offsets for the rows and columns of Z */
 /* IROWU, IROWV: offsets for the rows of U and V */
 *ns = 0;
 nru = 0;
 nrv = 0;
 idbeg = 1;
 isbeg = 1;
 irowz = 1;
 icolz = 1;
 irowu = 2;
 irowv = 1;
 split = FALSE_;
 sveq0 = FALSE_;
 /* Form the tridiagonal TGK matrix. */
 /* S( 1:N ) = ZERO */
 i__1 = *n;
 for (j1 = 1;
 j1 <= i__1;
 ++j1) {
 s[j1] = 0.f;
 }
 work[ietgk + (*n << 1) - 1] = 0.f;
 /* WORK( IDTGK:IDTGK+2*N-1 ) = ZERO */
 i__1 = *n << 1;
 for (j1 = 1;
 j1 <= i__1;
 ++j1) {
 work[idtgk - 1 + j1] = 0.f;
 }
 scopy_(n, &d__[1], &c__1, &work[ietgk], &c__2);
 i__1 = *n - 1;
 scopy_(&i__1, &e[1], &c__1, &work[ietgk + 1], &c__2);
 /* Check for splits in two levels, outer level */
 /* in E and inner level in D. */
 i__1 = *n << 1;
 for (ieptr = 2;
 ieptr <= i__1;
 ieptr += 2) {
 if (work[ietgk + ieptr - 1] == 0.f) {
 /* Split in E (this piece of B is square) or bottom */
 /* of the (input bidiagonal) matrix. */
 isplt = idbeg;
 idend = ieptr - 1;
 i__2 = idend;
 for (idptr = idbeg;
 idptr <= i__2;
 idptr += 2) {
 if (work[ietgk + idptr - 1] == 0.f) {
 /* Split in D (rectangular submatrix). Set the number */
 /* of rows in U and V (NRU and NRV) accordingly. */
 if (idptr == idbeg) {
 /* D=0 at the top. */
 sveq0 = TRUE_;
 if (idbeg == idend) {
 nru = 1;
 nrv = 1;
 }
 }
 else if (idptr == idend) {
 /* D=0 at the bottom. */
 sveq0 = TRUE_;
 nru = (idend - isplt) / 2 + 1;
 nrv = nru;
 if (isplt != idbeg) {
 ++nru;
 }
 }
 else {
 if (isplt == idbeg) {
 /* Split: top rectangular submatrix. */
 nru = (idptr - idbeg) / 2;
 nrv = nru + 1;
 }
 else {
 /* Split: middle square submatrix. */
 nru = (idptr - isplt) / 2 + 1;
 nrv = nru;
 }
 }
 }
 else if (idptr == idend) {
 /* Last entry of D in the active submatrix. */
 if (isplt == idbeg) {
 /* No split (trivial case). */
 nru = (idend - idbeg) / 2 + 1;
 nrv = nru;
 }
 else {
 /* Split: bottom rectangular submatrix. */
 nrv = (idend - isplt) / 2 + 1;
 nru = nrv + 1;
 }
 }
 ntgk = nru + nrv;
 if (ntgk > 0) {
 /* Compute eigenvalues/vectors of the active */
 /* submatrix according to RANGE: */
 /* if RANGE='A' (ALLSV) then RNGVX = 'I' */
 /* if RANGE='V' (VALSV) then RNGVX = 'V' */
 /* if RANGE='I' (INDSV) then RNGVX = 'V' */
 iltgk = 1;
 iutgk = ntgk / 2;
 if (allsv || vutgk == 0.f) {
 if (sveq0 || smin < eps || ntgk % 2 > 0) {
 /* Special case: eigenvalue equal to zero or very */
 /* small, additional eigenvector is needed. */
 ++iutgk;
 }
 }
 /* Workspace needed by SSTEVX: */
 /* WORK( ITEMP: ): 2*5*NTGK */
 /* IWORK( 1: ): 2*6*NTGK */
 sstevx_(jobz, rngvx, &ntgk, &work[idtgk + isplt - 1], & work[ietgk + isplt - 1], &vltgk, &vutgk, &iltgk, & iutgk, &abstol, &nsl, &s[isbeg], &z__[irowz + icolz * z_dim1], ldz, &work[itemp], &iwork[iiwork] , &iwork[iifail], info);
 if (*info != 0) {
 /* Exit with the error code from SSTEVX. */
 return 0;
 }
 /* EMIN = ABS( MAXVAL( S( ISBEG:ISBEG+NSL-1 ) ) ) */
 i__3 = nsl;
 for (j1 = 1;
 j1 <= i__3;
 ++j1) {
 /* Computing MAX */
 r__1 = d1; r__2 = s[j1 - 1 + isbeg]; // , expr subst  
 d1 = max(r__1,r__2);
 }
 emin = abs(d1);
 if (nsl > 0 && wantz) {
 /* Normalize u=Z([2,4,...],:) and v=Z([1,3,...],:), */
 /* changing the sign of v as discussed in the leading */
 /* comments. The norms of u and v may be (slightly) */
 /* different from 1/sqrt(2) if the corresponding */
 /* eigenvalues are very small or too close. We check */
 /* those norms and, if needed, reorthogonalize the */
 /* vectors. */
 if (nsl > 1 && vutgk == 0.f && ntgk % 2 == 0 && emin == 0.f && ! split) {
 /* D=0 at the top or bottom of the active submatrix: */
 /* one eigenvalue is equal to zero;
 concatenate the */
 /* eigenvectors corresponding to the two smallest */
 /* eigenvalues. */
 /* Z( IROWZ:IROWZ+NTGK-1,ICOLZ+NSL-2 ) = */
 /* $ Z( IROWZ:IROWZ+NTGK-1,ICOLZ+NSL-2 ) + */
 /* $ Z( IROWZ:IROWZ+NTGK-1,ICOLZ+NSL-1 ) */
 /* Z( IROWZ:IROWZ+NTGK-1,ICOLZ+NSL-1 ) = */
 /* $ ZERO */
 i__3 = ntgk;
 for (j1 = 1;
 j1 <= i__3;
 ++j1) {
 z__[j1 - 1 + irowz + (icolz + nsl - 2) * z_dim1] += z__[j1 - 1 + irowz + ( icolz + nsl - 1) * z_dim1];
 z__[j1 - 1 + irowz + (icolz + nsl - 1) * z_dim1] = 0.f;
 }
 /* IF( IUTGK*2.GT.NTGK ) THEN */
 /* Eigenvalue equal to zero or very small. */
 /* NSL = NSL - 1 */
 /* END IF */
 }
 /* Computing MIN */
 i__4 = nsl - 1; i__5 = nru - 1; // , expr subst  
 i__3 = min(i__4,i__5);
 for (i__ = 0;
 i__ <= i__3;
 ++i__) {
 nrmu = snrm2_(&nru, &z__[irowu + (icolz + i__) * z_dim1], &c__2);
 if (nrmu == 0.f) {
 *info = (*n << 1) + 1;
 return 0;
 }
 r__1 = 1.f / nrmu;
 sscal_(&nru, &r__1, &z__[irowu + (icolz + i__) * z_dim1], &c__2);
 if (nrmu != 1.f && (r__1 = nrmu - ortol, abs(r__1) ) * sqrt2 > 1.f) {
 i__4 = i__ - 1;
 for (j = 0;
 j <= i__4;
 ++j) {
 zjtji = -sdot_(&nru, &z__[irowu + (icolz + j) * z_dim1], &c__2, &z__[irowu + (icolz + i__) * z_dim1], &c__2);
 saxpy_(&nru, &zjtji, &z__[irowu + (icolz + j) * z_dim1], &c__2, &z__[irowu + (icolz + i__) * z_dim1], &c__2);
 }
 nrmu = snrm2_(&nru, &z__[irowu + (icolz + i__) * z_dim1], &c__2);
 r__1 = 1.f / nrmu;
 sscal_(&nru, &r__1, &z__[irowu + (icolz + i__) * z_dim1], &c__2);
 }
 }
 /* Computing MIN */
 i__4 = nsl - 1; i__5 = nrv - 1; // , expr subst  
 i__3 = min(i__4,i__5);
 for (i__ = 0;
 i__ <= i__3;
 ++i__) {
 nrmv = snrm2_(&nrv, &z__[irowv + (icolz + i__) * z_dim1], &c__2);
 if (nrmv == 0.f) {
 *info = (*n << 1) + 1;
 return 0;
 }
 r__1 = -1.f / nrmv;
 sscal_(&nrv, &r__1, &z__[irowv + (icolz + i__) * z_dim1], &c__2);
 if (nrmv != 1.f && (r__1 = nrmv - ortol, abs(r__1) ) * sqrt2 > 1.f) {
 i__4 = i__ - 1;
 for (j = 0;
 j <= i__4;
 ++j) {
 zjtji = -sdot_(&nrv, &z__[irowv + (icolz + j) * z_dim1], &c__2, &z__[irowv + (icolz + i__) * z_dim1], &c__2);
 saxpy_(&nru, &zjtji, &z__[irowv + (icolz + j) * z_dim1], &c__2, &z__[irowv + (icolz + i__) * z_dim1], &c__2);
 }
 nrmv = snrm2_(&nrv, &z__[irowv + (icolz + i__) * z_dim1], &c__2);
 r__1 = 1.f / nrmv;
 sscal_(&nrv, &r__1, &z__[irowv + (icolz + i__) * z_dim1], &c__2);
 }
 }
 if (vutgk == 0.f && idptr < idend && ntgk % 2 > 0) {
 /* D=0 in the middle of the active submatrix (one */
 /* eigenvalue is equal to zero): save the corresponding */
 /* eigenvector for later use (when bottom of the */
 /* active submatrix is reached). */
 split = TRUE_;
 /* Z( IROWZ:IROWZ+NTGK-1,N+1 ) = */
 /* $ Z( IROWZ:IROWZ+NTGK-1,NS+NSL ) */
 /* Z( IROWZ:IROWZ+NTGK-1,NS+NSL ) = */
 /* $ ZERO */
 i__3 = ntgk;
 for (j1 = 1;
 j1 <= i__3;
 ++j1) {
 z__[j1 - 1 + irowz + (*n + 1) * z_dim1] = z__[ j1 - 1 + irowz + (*ns + nsl) * z_dim1] ;
 z__[j1 - 1 + irowz + (*ns + nsl) * z_dim1] = 0.f;
 }
 }
 }
 /* ** WANTZ **! */
 nsl = min(nsl,nru);
 sveq0 = FALSE_;
 /* Absolute values of the eigenvalues of TGK. */
 i__3 = nsl - 1;
 for (i__ = 0;
 i__ <= i__3;
 ++i__) {
 s[isbeg + i__] = (r__1 = s[isbeg + i__], abs(r__1));
 }
 /* Update pointers for TGK, S and Z. */
 isbeg += nsl;
 irowz += ntgk;
 icolz += nsl;
 irowu = irowz;
 irowv = irowz + 1;
 isplt = idptr + 1;
 *ns += nsl;
 nru = 0;
 nrv = 0;
 }
 /* ** NTGK.GT.0 **! */
 if (irowz < *n << 1 && wantz) {
 /* Z( 1:IROWZ-1, ICOLZ ) = ZERO */
 i__3 = irowz - 1;
 for (j1 = 1;
 j1 <= i__3;
 ++j1) {
 z__[j1 + icolz * z_dim1] = 0.f;
 }
 }
 }
 /* ** IDPTR loop **! */
 if (split && wantz) {
 /* Bring back eigenvector corresponding */
 /* to eigenvalue equal to zero. */
 /* Z( IDBEG:IDEND-NTGK+1,ISBEG-1 ) = */
 /* $ Z( IDBEG:IDEND-NTGK+1,ISBEG-1 ) + */
 /* $ Z( IDBEG:IDEND-NTGK+1,N+1 ) */
 /* Z( IDBEG:IDEND-NTGK+1,N+1 ) = 0 */
 i__2 = idend - ntgk + 2 - idbeg;
 for (j1 = 1;
 j1 <= i__2;
 ++j1) {
 z__[j1 - 1 + idbeg + (isbeg - 1) * z_dim1] += z__[j1 - 1 + idbeg + (*n + 1) * z_dim1];
 z__[j1 - 1 + idbeg + (*n + 1) * z_dim1] = 0.f;
 }
 }
 --irowv;
 ++irowu;
 idbeg = ieptr + 1;
 sveq0 = FALSE_;
 split = FALSE_;
 }
 /* ** Check for split in E **! */
 }
 /* Sort the singular values into decreasing order (insertion sort on */
 /* singular values, but only one transposition per singular vector) */
 /* ** IEPTR loop **! */
 i__1 = *ns - 1;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 k = 1;
 smin = s[1];
 i__2 = *ns + 1 - i__;
 for (j = 2;
 j <= i__2;
 ++j) {
 if (s[j] <= smin) {
 k = j;
 smin = s[j];
 }
 }
 if (k != *ns + 1 - i__) {
 s[k] = s[*ns + 1 - i__];
 s[*ns + 1 - i__] = smin;
 if (wantz) {
 i__2 = *n << 1;
 sswap_(&i__2, &z__[k * z_dim1 + 1], &c__1, &z__[(*ns + 1 - i__) * z_dim1 + 1], &c__1);
 }
 }
 }
 /* If RANGE=I, check for singular values/vectors to be discarded. */
 if (indsv) {
 k = *iu - *il + 1;
 if (k < *ns) {
 /* S( K+1:NS ) = ZERO */
 i__1 = *ns - k;
 for (j1 = 1;
 j1 <= i__1;
 ++j1) {
 s[j1 + k] = 0.f;
 }
 /* IF( WANTZ ) Z( 1:N*2,K+1:NS ) = ZERO */
 if (wantz) {
 i__1 = *ns - k;
 for (j2 = 1;
 j2 <= i__1;
 ++j2) {
 i__2 = *n << 1;
 for (j1 = 1;
 j1 <= i__2;
 ++j1) {
 z__[j1 + (j2 + k) * z_dim1] = 0.f;
 }
 }
 }
 *ns = k;
 }
 }
 /* Reorder Z: U = Z( 1:N,1:NS ), V = Z( N+1:N*2,1:NS ). */
 /* If B is a lower diagonal, swap U and V. */
 if (wantz) {
 i__1 = *ns;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = *n << 1;
 scopy_(&i__2, &z__[i__ * z_dim1 + 1], &c__1, &work[1], &c__1);
 if (lower) {
 scopy_(n, &work[2], &c__2, &z__[*n + 1 + i__ * z_dim1], &c__1) ;
 scopy_(n, &work[1], &c__2, &z__[i__ * z_dim1 + 1], &c__1);
 }
 else {
 scopy_(n, &work[2], &c__2, &z__[i__ * z_dim1 + 1], &c__1);
 scopy_(n, &work[1], &c__2, &z__[*n + 1 + i__ * z_dim1], &c__1) ;
 }
 }
 }
 return 0;
 /* End of SBDSVDX */
 }
 /* sbdsvdx_ */
 