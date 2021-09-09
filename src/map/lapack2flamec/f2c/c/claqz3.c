/* claqz3.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static complex c_b1 = {
0.f,0.f}
;
 static complex c_b2 = {
1.f,0.f}
;
 static integer c__1 = 1;
 static logical c_true = TRUE_;
 /* > \brief \b CLAQZ3 */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download CLAQZ3 + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/CLAQZ3. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/CLAQZ3. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/CLAQZ3. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE CLAQZ3( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSHIFTS, */
 /* $ NBLOCK_DESIRED, ALPHA, BETA, A, LDA, B, LDB, Q, LDQ, Z, LDZ, */
 /* $ QC, LDQC, ZC, LDZC, WORK, LWORK, INFO ) */
 /* IMPLICIT NONE */
 /* Function arguments */
 /* LOGICAL, INTENT( IN ) :: ILSCHUR, ILQ, ILZ */
 /* INTEGER, INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, */
 /* $ NSHIFTS, NBLOCK_DESIRED, LDQC, LDZC */
 /* COMPLEX, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
 /* $ Z( LDZ, * ), QC( LDQC, * ), ZC( LDZC, * ), WORK( * ), */
 /* $ ALPHA( * ), BETA( * ) */
 /* INTEGER, INTENT( OUT ) :: INFO */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > CLAQZ3 Executes a single multishift QZ sweep */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] ILSCHUR */
 /* > \verbatim */
 /* > ILSCHUR is LOGICAL */
 /* > Determines whether or not to update the full Schur form */
 /* > \endverbatim */
 /* > \param[in] ILQ */
 /* > \verbatim */
 /* > ILQ is LOGICAL */
 /* > Determines whether or not to update the matrix Q */
 /* > \endverbatim */
 /* > */
 /* > \param[in] ILZ */
 /* > \verbatim */
 /* > ILZ is LOGICAL */
 /* > Determines whether or not to update the matrix Z */
 /* > \endverbatim */
 /* > */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The order of the matrices A, B, Q, and Z. N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] ILO */
 /* > \verbatim */
 /* > ILO is INTEGER */
 /* > \endverbatim */
 /* > */
 /* > \param[in] IHI */
 /* > \verbatim */
 /* > IHI is INTEGER */
 /* > \endverbatim */
 /* > */
 /* > \param[in] NSHIFTS */
 /* > \verbatim */
 /* > NSHIFTS is INTEGER */
 /* > The desired number of shifts to use */
 /* > \endverbatim */
 /* > */
 /* > \param[in] NBLOCK_DESIRED */
 /* > \verbatim */
 /* > NBLOCK_DESIRED is INTEGER */
 /* > The desired size of the computational windows */
 /* > \endverbatim */
 /* > */
 /* > \param[in] ALPHA */
 /* > \verbatim */
 /* > ALPHA is COMPLEX array. SR contains */
 /* > the alpha parts of the shifts to use. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] BETA */
 /* > \verbatim */
 /* > BETA is COMPLEX array. SS contains */
 /* > the scale of the shifts to use. */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] A */
 /* > \verbatim */
 /* > A is COMPLEX array, dimension (LDA, N) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max( 1, N ). */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] B */
 /* > \verbatim */
 /* > B is COMPLEX array, dimension (LDB, N) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDB */
 /* > \verbatim */
 /* > LDB is INTEGER */
 /* > The leading dimension of the array B. LDB >= max( 1, N ). */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] Q */
 /* > \verbatim */
 /* > Q is COMPLEX array, dimension (LDQ, N) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDQ */
 /* > \verbatim */
 /* > LDQ is INTEGER */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] Z */
 /* > \verbatim */
 /* > Z is COMPLEX array, dimension (LDZ, N) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDZ */
 /* > \verbatim */
 /* > LDZ is INTEGER */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] QC */
 /* > \verbatim */
 /* > QC is COMPLEX array, dimension (LDQC, NBLOCK_DESIRED) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDQC */
 /* > \verbatim */
 /* > LDQC is INTEGER */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] ZC */
 /* > \verbatim */
 /* > ZC is COMPLEX array, dimension (LDZC, NBLOCK_DESIRED) */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDZC */
 /* > \verbatim */
 /* > LDZ is INTEGER */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
 /* > On exit, if INFO >= 0, WORK(1) returns the optimal LWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > LWORK is INTEGER */
 /* > The dimension of the array WORK. LWORK >= max(1,N). */
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
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Thijs Steel, KU Leuven */
 /* > \date May 2020 */
 /* > \ingroup complexGEcomputational */
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int claqz3_(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nshifts, integer * nblock_desired__, complex *alpha, complex *beta, complex *a, integer * lda, complex *b, integer *ldb, complex *q, integer *ldq, complex *z__, integer *ldz, complex *qc, integer *ldqc, complex *zc, integer *ldzc, complex *work, integer *lwork, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, qc_dim1, qc_offset, zc_dim1, zc_offset, i__1, i__2, i__3, i__4, i__5;
 complex q__1, q__2, q__3;
 /* Builtin functions */
 double c_abs(complex *), sqrt(doublereal);
 void r_cnjg(complex *, complex *);
 /* Local variables */
 real c__;
 integer i__, j, k;
 complex s;
 integer np, ns;
 complex temp;
 extern /* Subroutine */
 int crot_(integer *, complex *, integer *, complex *, integer *, real *, complex *);
 integer npos;
 complex temp2, temp3;
 real scale;
 extern /* Subroutine */
 int cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *), claqz1_(logical *, logical *, integer *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, integer *, integer *, complex *, integer *, integer *, integer *, complex *, integer *), slabad_(real *, real *);
 extern real slamch_(char *);
 integer nblock;
 extern /* Subroutine */
 int claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *), clartg_(complex *, complex *, real *, complex *, complex *);
 real safmin;
 extern /* Subroutine */
 int xerbla_(char *, integer *);
 real safmax;
 extern /* Subroutine */
 int clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *);
 integer ishift, istopb, swidth, istopm, sheight, istartb, istartm;
 /* Function arguments */
 /* Parameters */
 /* Local scalars */
 /* External Functions */
 /* Parameter adjustments */
 --alpha;
 --beta;
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 b_dim1 = *ldb;
 b_offset = 1 + b_dim1;
 b -= b_offset;
 q_dim1 = *ldq;
 q_offset = 1 + q_dim1;
 q -= q_offset;
 z_dim1 = *ldz;
 z_offset = 1 + z_dim1;
 z__ -= z_offset;
 qc_dim1 = *ldqc;
 qc_offset = 1 + qc_dim1;
 qc -= qc_offset;
 zc_dim1 = *ldzc;
 zc_offset = 1 + zc_dim1;
 zc -= zc_offset;
 --work;
 /* Function Body */
 *info = 0;
 if (*nblock_desired__ < *nshifts + 1) {
 *info = -8;
 }
 if (*lwork == -1) {
 /* workspace query, quick return */
 i__1 = *n * *nblock_desired__;
 work[1].r = (real) i__1; work[1].i = 0.f; // , expr subst  
 return 0;
 }
 else if (*lwork < *n * *nblock_desired__) {
 *info = -25;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("CLAQZ3", &i__1);
 return 0;
 }
 /* Executable statements */
 /* Get machine constants */
 safmin = slamch_("SAFE MINIMUM");
 safmax = 1.f / safmin;
 slabad_(&safmin, &safmax);
 if (*ilo >= *ihi) {
 return 0;
 }
 if (*ilschur) {
 istartm = 1;
 istopm = *n;
 }
 else {
 istartm = *ilo;
 istopm = *ihi;
 }
 ns = *nshifts;
 /* Computing MAX */
 i__1 = *nblock_desired__ - ns;
 npos = max(i__1,1);
 /* The following block introduces the shifts and chases */
 /* them down one by one just enough to make space for */
 /* the other shifts. The near-the-diagonal block is */
 /* of size (ns+1) x ns. */
 i__1 = ns + 1;
 i__2 = ns + 1;
 claset_("FULL", &i__1, &i__2, &c_b1, &c_b2, &qc[qc_offset], ldqc);
 claset_("FULL", &ns, &ns, &c_b1, &c_b2, &zc[zc_offset], ldzc);
 i__1 = ns;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 /* Introduce the shift */
 scale = sqrt(c_abs(&alpha[i__])) * sqrt(c_abs(&beta[i__]));
 if (scale >= safmin && scale <= safmax) {
 i__2 = i__;
 i__3 = i__;
 q__1.r = alpha[i__3].r / scale; q__1.i = alpha[i__3].i / scale; // , expr subst  
 alpha[i__2].r = q__1.r; alpha[i__2].i = q__1.i; // , expr subst  
 i__2 = i__;
 i__3 = i__;
 q__1.r = beta[i__3].r / scale; q__1.i = beta[i__3].i / scale; // , expr subst  
 beta[i__2].r = q__1.r; beta[i__2].i = q__1.i; // , expr subst  
 }
 i__2 = i__;
 i__3 = *ilo + *ilo * a_dim1;
 q__2.r = beta[i__2].r * a[i__3].r - beta[i__2].i * a[i__3].i; q__2.i = beta[i__2].r * a[i__3].i + beta[i__2].i * a[i__3].r; // , expr subst  
 i__4 = i__;
 i__5 = *ilo + *ilo * b_dim1;
 q__3.r = alpha[i__4].r * b[i__5].r - alpha[i__4].i * b[i__5].i; q__3.i = alpha[i__4].r * b[i__5].i + alpha[i__4].i * b[i__5] .r; // , expr subst  
 q__1.r = q__2.r - q__3.r; q__1.i = q__2.i - q__3.i; // , expr subst  
 temp2.r = q__1.r; temp2.i = q__1.i; // , expr subst  
 i__2 = i__;
 i__3 = *ilo + 1 + *ilo * a_dim1;
 q__1.r = beta[i__2].r * a[i__3].r - beta[i__2].i * a[i__3].i; q__1.i = beta[i__2].r * a[i__3].i + beta[i__2].i * a[i__3].r; // , expr subst  
 temp3.r = q__1.r; temp3.i = q__1.i; // , expr subst  
 if (c_abs(&temp2) > safmax || c_abs(&temp3) > safmax) {
 temp2.r = 1.f; temp2.i = 0.f; // , expr subst  
 temp3.r = 0.f; temp3.i = 0.f; // , expr subst  
 }
 clartg_(&temp2, &temp3, &c__, &s, &temp);
 crot_(&ns, &a[*ilo + *ilo * a_dim1], lda, &a[*ilo + 1 + *ilo * a_dim1] , lda, &c__, &s);
 crot_(&ns, &b[*ilo + *ilo * b_dim1], ldb, &b[*ilo + 1 + *ilo * b_dim1] , ldb, &c__, &s);
 i__2 = ns + 1;
 r_cnjg(&q__1, &s);
 crot_(&i__2, &qc[qc_dim1 + 1], &c__1, &qc[(qc_dim1 << 1) + 1], &c__1, &c__, &q__1);
 /* Chase the shift down */
 i__2 = ns - i__;
 for (j = 1;
 j <= i__2;
 ++j) {
 i__3 = *ihi - *ilo + 1;
 i__4 = ns + 1;
 claqz1_(&c_true, &c_true, &j, &c__1, &ns, &i__3, &a[*ilo + *ilo * a_dim1], lda, &b[*ilo + *ilo * b_dim1], ldb, &i__4, &c__1, &qc[qc_offset], ldqc, &ns, &c__1, &zc[zc_offset], ldzc);
 }
 }
 /* Update the rest of the pencil */
 /* Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm) */
 /* from the left with Qc(1:ns+1,1:ns+1)' */
 sheight = ns + 1;
 swidth = istopm - (*ilo + ns) + 1;
 if (swidth > 0) {
 cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[qc_offset], ldqc, &a[*ilo + (*ilo + ns) * a_dim1], lda, &c_b1, &work[1], & sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[*ilo + (*ilo + ns) * a_dim1], lda);
 cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[qc_offset], ldqc, &b[*ilo + (*ilo + ns) * b_dim1], ldb, &c_b1, &work[1], & sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[*ilo + (*ilo + ns) * b_dim1], ldb);
 }
 if (*ilq) {
 cgemm_("N", "N", n, &sheight, &sheight, &c_b2, &q[*ilo * q_dim1 + 1], ldq, &qc[qc_offset], ldqc, &c_b1, &work[1], n);
 clacpy_("ALL", n, &sheight, &work[1], n, &q[*ilo * q_dim1 + 1], ldq);
 }
 /* Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1) */
 /* from the right with Zc(1:ns,1:ns) */
 sheight = *ilo - 1 - istartm + 1;
 swidth = ns;
 if (sheight > 0) {
 cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &a[istartm + *ilo * a_dim1], lda, &zc[zc_offset], ldzc, &c_b1, &work[1], & sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + * ilo * a_dim1], lda);
 cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &b[istartm + *ilo * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b1, &work[1], & sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + * ilo * b_dim1], ldb);
 }
 if (*ilz) {
 cgemm_("N", "N", n, &swidth, &swidth, &c_b2, &z__[*ilo * z_dim1 + 1], ldz, &zc[zc_offset], ldzc, &c_b1, &work[1], n);
 clacpy_("ALL", n, &swidth, &work[1], n, &z__[*ilo * z_dim1 + 1], ldz);
 }
 /* The following block chases the shifts down to the bottom */
 /* right block. If possible, a shift is moved down npos */
 /* positions at a time */
 k = *ilo;
 while(k < *ihi - ns) {
 /* Computing MIN */
 i__1 = *ihi - ns - k;
 np = min(i__1,npos);
 /* Size of the near-the-diagonal block */
 nblock = ns + np;
 /* istartb points to the first row we will be updating */
 istartb = k + 1;
 /* istopb points to the last column we will be updating */
 istopb = k + nblock - 1;
 i__1 = ns + np;
 i__2 = ns + np;
 claset_("FULL", &i__1, &i__2, &c_b1, &c_b2, &qc[qc_offset], ldqc);
 i__1 = ns + np;
 i__2 = ns + np;
 claset_("FULL", &i__1, &i__2, &c_b1, &c_b2, &zc[zc_offset], ldzc);
 /* Near the diagonal shift chase */
 for (i__ = ns - 1;
 i__ >= 0;
 --i__) {
 i__1 = np - 1;
 for (j = 0;
 j <= i__1;
 ++j) {
 /* Move down the block with index k+i+j, updating */
 /* the (ns+np x ns+np) block: */
 /* (k:k+ns+np,k:k+ns+np-1) */
 i__2 = k + i__ + j;
 i__3 = k + 1;
 claqz1_(&c_true, &c_true, &i__2, &istartb, &istopb, ihi, &a[ a_offset], lda, &b[b_offset], ldb, &nblock, &i__3, & qc[qc_offset], ldqc, &nblock, &k, &zc[zc_offset], ldzc);
 }
 }
 /* Update rest of the pencil */
 /* Update A(k+1:k+ns+np, k+ns+np:istopm) and */
 /* B(k+1:k+ns+np, k+ns+np:istopm) */
 /* from the left with Qc(1:ns+np,1:ns+np)' */
 sheight = ns + np;
 swidth = istopm - (k + ns + np) + 1;
 if (swidth > 0) {
 cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[ qc_offset], ldqc, &a[k + 1 + (k + ns + np) * a_dim1], lda, &c_b1, &work[1], &sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[k + 1 + ( k + ns + np) * a_dim1], lda);
 cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[ qc_offset], ldqc, &b[k + 1 + (k + ns + np) * b_dim1], ldb, &c_b1, &work[1], &sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[k + 1 + ( k + ns + np) * b_dim1], ldb);
 }
 if (*ilq) {
 cgemm_("N", "N", n, &nblock, &nblock, &c_b2, &q[(k + 1) * q_dim1 + 1], ldq, &qc[qc_offset], ldqc, &c_b1, &work[1], n);
 clacpy_("ALL", n, &nblock, &work[1], n, &q[(k + 1) * q_dim1 + 1], ldq);
 }
 /* Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1) */
 /* from the right with Zc(1:ns+np,1:ns+np) */
 sheight = k - istartm + 1;
 swidth = nblock;
 if (sheight > 0) {
 cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &a[istartm + k * a_dim1], lda, &zc[zc_offset], ldzc, &c_b1, &work[1], & sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + k * a_dim1], lda);
 cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &b[istartm + k * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b1, &work[1], & sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + k * b_dim1], ldb);
 }
 if (*ilz) {
 cgemm_("N", "N", n, &nblock, &nblock, &c_b2, &z__[k * z_dim1 + 1], ldz, &zc[zc_offset], ldzc, &c_b1, &work[1], n);
 clacpy_("ALL", n, &nblock, &work[1], n, &z__[k * z_dim1 + 1], ldz);
 }
 k += np;
 }
 /* The following block removes the shifts from the bottom right corner */
 /* one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi). */
 claset_("FULL", &ns, &ns, &c_b1, &c_b2, &qc[qc_offset], ldqc);
 i__1 = ns + 1;
 i__2 = ns + 1;
 claset_("FULL", &i__1, &i__2, &c_b1, &c_b2, &zc[zc_offset], ldzc);
 /* istartb points to the first row we will be updating */
 istartb = *ihi - ns + 1;
 /* istopb points to the last column we will be updating */
 istopb = *ihi;
 i__1 = ns;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 /* Chase the shift down to the bottom right corner */
 i__2 = *ihi - 1;
 for (ishift = *ihi - i__;
 ishift <= i__2;
 ++ishift) {
 i__3 = *ihi - ns + 1;
 i__4 = ns + 1;
 i__5 = *ihi - ns;
 claqz1_(&c_true, &c_true, &ishift, &istartb, &istopb, ihi, &a[ a_offset], lda, &b[b_offset], ldb, &ns, &i__3, &qc[ qc_offset], ldqc, &i__4, &i__5, &zc[zc_offset], ldzc);
 }
 }
 /* Update rest of the pencil */
 /* Update A(ihi-ns+1:ihi, ihi+1:istopm) */
 /* from the left with Qc(1:ns,1:ns)' */
 sheight = ns;
 swidth = istopm - (*ihi + 1) + 1;
 if (swidth > 0) {
 cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[qc_offset], ldqc, &a[*ihi - ns + 1 + (*ihi + 1) * a_dim1], lda, &c_b1, & work[1], &sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[*ihi - ns + 1 + (*ihi + 1) * a_dim1], lda);
 cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[qc_offset], ldqc, &b[*ihi - ns + 1 + (*ihi + 1) * b_dim1], ldb, &c_b1, & work[1], &sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[*ihi - ns + 1 + (*ihi + 1) * b_dim1], ldb);
 }
 if (*ilq) {
 cgemm_("N", "N", n, &ns, &ns, &c_b2, &q[(*ihi - ns + 1) * q_dim1 + 1], ldq, &qc[qc_offset], ldqc, &c_b1, &work[1], n);
 clacpy_("ALL", n, &ns, &work[1], n, &q[(*ihi - ns + 1) * q_dim1 + 1], ldq);
 }
 /* Update A(istartm:ihi-ns,ihi-ns:ihi) */
 /* from the right with Zc(1:ns+1,1:ns+1) */
 sheight = *ihi - ns - istartm + 1;
 swidth = ns + 1;
 if (sheight > 0) {
 cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &a[istartm + (* ihi - ns) * a_dim1], lda, &zc[zc_offset], ldzc, &c_b1, &work[ 1], &sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + (* ihi - ns) * a_dim1], lda);
 cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &b[istartm + (* ihi - ns) * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b1, &work[ 1], &sheight);
 clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + (* ihi - ns) * b_dim1], ldb);
 }
 if (*ilz) {
 i__1 = ns + 1;
 i__2 = ns + 1;
 cgemm_("N", "N", n, &i__1, &i__2, &c_b2, &z__[(*ihi - ns) * z_dim1 + 1], ldz, &zc[zc_offset], ldzc, &c_b1, &work[1], n);
 i__1 = ns + 1;
 clacpy_("ALL", n, &i__1, &work[1], n, &z__[(*ihi - ns) * z_dim1 + 1], ldz);
 }
 return 0;
 }
 /* claqz3_ */
 