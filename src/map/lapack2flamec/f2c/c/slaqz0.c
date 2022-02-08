/* slaqz0.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__17 = 17;
static integer c_n1 = -1;
static real c_b25 = 0.f;
static real c_b26 = 1.f;
static integer c__1 = 1;
/* > \brief \b SLAQZ0 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAQZ0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqz0. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqz0. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqz0. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAQZ0( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, */
/* $ LDB, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, REC, */
/* $ INFO ) */
/* IMPLICIT NONE */
/* Arguments */
/* CHARACTER, INTENT( IN ) :: WANTS, WANTQ, WANTZ */
/* INTEGER, INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, */
/* $ REC */
/* INTEGER, INTENT( OUT ) :: INFO */
/* REAL, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ Z( LDZ, * ), ALPHAR( * ), ALPHAI( * ), BETA( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQZ0 computes the eigenvalues of a real matrix pair (H,T), */
/* > where H is an upper Hessenberg matrix and T is upper triangular, */
/* > using the double-shift QZ method. */
/* > Matrix pairs of this type are produced by the reduction to */
/* > generalized upper Hessenberg form of a real matrix pair (A,B): */
/* > */
/* > A = Q1*H*Z1**T, B = Q1*T*Z1**T, */
/* > */
/* > as computed by SGGHRD. */
/* > */
/* > If JOB='S', then the Hessenberg-triangular pair (H,T) is */
/* > also reduced to generalized Schur form, */
/* > */
/* > H = Q*S*Z**T, T = Q*P*Z**T, */
/* > */
/* > where Q and Z are orthogonal matrices, P is an upper triangular */
/* > matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2 */
/* > diagonal blocks. */
/* > */
/* > The 1-by-1 blocks correspond to real eigenvalues of the matrix pair */
/* > (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of */
/* > eigenvalues. */
/* > */
/* > Additionally, the 2-by-2 upper triangular diagonal blocks of P */
/* > corresponding to 2-by-2 blocks of S are reduced to positive diagonal */
/* > form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0, */
/* > P(j,j) > 0, and P(j+1,j+1) > 0. */
/* > */
/* > Optionally, the orthogonal matrix Q from the generalized Schur */
/* > factorization may be postmultiplied into an input matrix Q1, and the */
/* > orthogonal matrix Z may be postmultiplied into an input matrix Z1. */
/* > If Q1 and Z1 are the orthogonal matrices from SGGHRD that reduced */
/* > the matrix pair (A,B) to generalized upper Hessenberg form, then the */
/* > output matrices Q1*Q and Z1*Z are the orthogonal factors from the */
/* > generalized Schur factorization of (A,B): */
/* > */
/* > A = (Q1*Q)*S*(Z1*Z)**T, B = (Q1*Q)*P*(Z1*Z)**T. */
/* > */
/* > To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently, */
/* > of (A,B)) are computed as a pair of values (alpha,beta), where alpha is */
/* > complex and beta real. */
/* > If beta is nonzero, lambda = alpha / beta is an eigenvalue of the */
/* > generalized nonsymmetric eigenvalue problem (GNEP) */
/* > A*x = lambda*B*x */
/* > and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the */
/* > alternate form of the GNEP */
/* > mu*A*y = B*y. */
/* > Real eigenvalues can be read directly from the generalized Schur */
/* > form: */
/* > alpha = S(i,i), beta = P(i,i). */
/* > */
/* > Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix */
/* > Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973), */
/* > pp. 241--256. */
/* > */
/* > Ref: B. Kagstrom, D. Kressner, "Multishift Variants of the QZ */
/* > Algorithm with Aggressive Early Deflation", SIAM J. Numer. */
/* > Anal., 29(2006), pp. 199--227. */
/* > */
/* > Ref: T. Steel, D. Camps, K. Meerbergen, R. Vandebril "A multishift, */
/* > multipole rational QZ method with agressive early deflation" */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTS */
/* > \verbatim */
/* > WANTS is CHARACTER*1 */
/* > = 'E': Compute eigenvalues only;
*/
/* > = 'S': Compute eigenvalues and the Schur form. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTQ */
/* > \verbatim */
/* > WANTQ is CHARACTER*1 */
/* > = 'N': Left Schur vectors (Q) are not computed;
*/
/* > = 'I': Q is initialized to the unit matrix and the matrix Q */
/* > of left Schur vectors of (A,B) is returned;
*/
/* > = 'V': Q must contain an orthogonal matrix Q1 on entry and */
/* > the product Q1*Q is returned. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is CHARACTER*1 */
/* > = 'N': Right Schur vectors (Z) are not computed;
*/
/* > = 'I': Z is initialized to the unit matrix and the matrix Z */
/* > of right Schur vectors of (A,B) is returned;
*/
/* > = 'V': Z must contain an orthogonal matrix Z1 on entry and */
/* > the product Z1*Z is returned. */
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
/* > ILO and IHI mark the rows and columns of A which are in */
/* > Hessenberg form. It is assumed that A is already upper */
/* > triangular in rows and columns 1:ILO-1 and IHI+1:N. */
/* > If N > 0, 1 <= ILO <= IHI <= N;
if N = 0, ILO=1 and IHI=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, N) */
/* > On entry, the N-by-N upper Hessenberg matrix A. */
/* > On exit, if JOB = 'S', A contains the upper quasi-triangular */
/* > matrix S from the generalized Schur factorization. */
/* > If JOB = 'E', the diagonal blocks of A match those of S, but */
/* > the rest of A is unspecified. */
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
/* > B is REAL array, dimension (LDB, N) */
/* > On entry, the N-by-N upper triangular matrix B. */
/* > On exit, if JOB = 'S', B contains the upper triangular */
/* > matrix P from the generalized Schur factorization;
*/
/* > 2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S */
/* > are reduced to positive diagonal form, i.e., if A(j+1,j) is */
/* > non-zero, then B(j+1,j) = B(j,j+1) = 0, B(j,j) > 0, and */
/* > B(j+1,j+1) > 0. */
/* > If JOB = 'E', the diagonal blocks of B match those of P, but */
/* > the rest of B is unspecified. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAR */
/* > \verbatim */
/* > ALPHAR is REAL array, dimension (N) */
/* > The real parts of each scalar alpha defining an eigenvalue */
/* > of GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* > ALPHAI is REAL array, dimension (N) */
/* > The imaginary parts of each scalar alpha defining an */
/* > eigenvalue of GNEP. */
/* > If ALPHAI(j) is zero, then the j-th eigenvalue is real;
if */
/* > positive, then the j-th and (j+1)-st eigenvalues are a */
/* > complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is REAL array, dimension (N) */
/* > The scalars beta that define the eigenvalues of GNEP. */
/* > Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and */
/* > beta = BETA(j) represent the j-th eigenvalue of the matrix */
/* > pair (A,B), in one of the forms lambda = alpha/beta or */
/* > mu = beta/alpha. Since either lambda or mu may overflow, */
/* > they should not, in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDQ, N) */
/* > On entry, if COMPQ = 'V', the orthogonal matrix Q1 used in */
/* > the reduction of (A,B) to generalized Hessenberg form. */
/* > On exit, if COMPQ = 'I', the orthogonal matrix of left Schur */
/* > vectors of (A,B), and if COMPQ = 'V', the orthogonal matrix */
/* > of left Schur vectors of (A,B). */
/* > Not referenced if COMPQ = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= 1. */
/* > If COMPQ='V' or 'I', then LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ, N) */
/* > On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in */
/* > the reduction of (A,B) to generalized Hessenberg form. */
/* > On exit, if COMPZ = 'I', the orthogonal matrix of */
/* > right Schur vectors of (H,T), and if COMPZ = 'V', the */
/* > orthogonal matrix of right Schur vectors of (A,B). */
/* > Not referenced if COMPZ = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1. */
/* > If COMPZ='V' or 'I', then LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
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
/* > \param[in] REC */
/* > \verbatim */
/* > REC is INTEGER */
/* > REC indicates the current recursion level. Should be set */
/* > to 0 on first call. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > = 1,...,N: the QZ iteration did not converge. (A,B) is not */
/* > in Schur form, but ALPHAR(i), ALPHAI(i), and */
/* > BETA(i), i=INFO+1,...,N should be correct. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Thijs Steel, KU Leuven */
/* > \date May 2020 */
/* > \ingroup doubleGEcomputational */
/* > */
/* ===================================================================== */
/* Subroutine */
int slaqz0_(char *wants, char *wantq, char *wantz, integer * n, integer *ilo, integer *ihi, real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *q, integer *ldq, real *z__, integer *ldz, real *work, integer *lwork, integer *rec, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer aed_info__;
    extern /* Subroutine */
    int f90_cycle_(void);
    integer shiftpos, lworkreq, i__, k;
    real c1;
    integer k2;
    real s1;
    integer norm_info__, ld, ns, n_deflated__, nw, sweep_info__, nbr;
    logical ilq, ilz;
    real ulp;
    integer nsr, nwr, nmin;
    real temp, swap;
    integer n_undeflated__;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *);
    extern logical lsame_(char *, char *);
    integer iiter, maxit, rcost, istop, itemp1, itemp2;
    extern /* Subroutine */
    int slaqz3_(logical *, logical *, logical *, integer *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *, integer *, real *, real *, real *, real *, integer *, real *, integer *, real *, integer *, integer *, integer *), slaqz4_( logical *, logical *, logical *, integer *, integer *, integer *, integer *, integer *, real *, real *, real *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *), slabad_(real *, real *);
    integer nibble, nblock;
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real safmax;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    real eshift;
    char jbcmpz[3];
    extern /* Subroutine */
    int slaset_(char *, integer *, integer *, real *, real *, real *, integer *), slartg_(real *, real *, real *, real *, real *), shgeqz_(char *, char *, char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, real *, real *, integer *, real *, integer *, real *, integer *, integer *);
    integer iwantq, iwants, istart;
    real smlnum;
    integer istopm, iwantz, istart2;
    extern /* Subroutine */
    int f90_exit_(void);
    logical ilschur;
    integer nshifts, istartm;
    /* Arguments */
    /* const */
    /* Parameters */
    /* Local scalars */
    /* External Functions */
    /* Decode wantS,wantQ,wantZ */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alphar;
    --alphai;
    --beta;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    /* Function Body */
    if (lsame_(wants, "E"))
    {
        ilschur = FALSE_;
        iwants = 1;
    }
    else if (lsame_(wants, "S"))
    {
        ilschur = TRUE_;
        iwants = 2;
    }
    else
    {
        iwants = 0;
    }
    if (lsame_(wantq, "N"))
    {
        ilq = FALSE_;
        iwantq = 1;
    }
    else if (lsame_(wantq, "V"))
    {
        ilq = TRUE_;
        iwantq = 2;
    }
    else if (lsame_(wantq, "I"))
    {
        ilq = TRUE_;
        iwantq = 3;
    }
    else
    {
        iwantq = 0;
    }
    if (lsame_(wantz, "N"))
    {
        ilz = FALSE_;
        iwantz = 1;
    }
    else if (lsame_(wantz, "V"))
    {
        ilz = TRUE_;
        iwantz = 2;
    }
    else if (lsame_(wantz, "I"))
    {
        ilz = TRUE_;
        iwantz = 3;
    }
    else
    {
        iwantz = 0;
    }
    /* Check Argument Values */
    *info = 0;
    if (iwants == 0)
    {
        *info = -1;
    }
    else if (iwantq == 0)
    {
        *info = -2;
    }
    else if (iwantz == 0)
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*ilo < 1)
    {
        *info = -5;
    }
    else if (*ihi > *n || *ihi < *ilo - 1)
    {
        *info = -6;
    }
    else if (*lda < *n)
    {
        *info = -8;
    }
    else if (*ldb < *n)
    {
        *info = -10;
    }
    else if (*ldq < 1 || ilq && *ldq < *n)
    {
        *info = -15;
    }
    else if (*ldz < 1 || ilz && *ldz < *n)
    {
        *info = -17;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLAQZ0", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n <= 0)
    {
        work[1] = 1.f;
        return 0;
    }
    /* Get the parameters */
    *(unsigned char *)jbcmpz = *(unsigned char *)wants;
    *(unsigned char *)&jbcmpz[1] = *(unsigned char *)wantq;
    *(unsigned char *)&jbcmpz[2] = *(unsigned char *)wantz;
    nmin = ilaenv_(&c__12, "SLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    nwr = ilaenv_(&c__13, "SLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    nwr = max(2,nwr);
    /* Computing MIN */
    i__1 = *ihi - *ilo + 1;
    i__2 = (*n - 1) / 3;
    i__1 = min(i__1,i__2); // ; expr subst
    nwr = min(i__1,nwr);
    nibble = ilaenv_(&c__14, "SLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    nsr = ilaenv_(&c__15, "SLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    /* Computing MIN */
    i__1 = nsr, i__2 = (*n + 6) / 9;
    i__1 = min(i__1,i__2);
    i__2 = *ihi - * ilo; // ; expr subst
    nsr = min(i__1,i__2);
    /* Computing MAX */
    i__1 = 2;
    i__2 = nsr - nsr % 2; // , expr subst
    nsr = max(i__1,i__2);
    rcost = ilaenv_(&c__17, "SLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    itemp1 = (integer) (nsr / sqrt((nsr << 1) / ((real) rcost / 100 * *n) + 1) );
    itemp1 = ((itemp1 - 1) / 4 << 2) + 4;
    nbr = nsr + itemp1;
    if (*n < nmin || *rec >= 2)
    {
        shgeqz_(wants, wantq, wantz, n, ilo, ihi, &a[a_offset], lda, &b[ b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &q[q_offset], ldq, &z__[z_offset], ldz, &work[1], lwork, info);
        return 0;
    }
    /* Find out required workspace */
    /* Workspace query to slaqz3 */
    nw = max(nwr,nmin);
    slaqz3_(&ilschur, &ilq, &ilz, n, ilo, ihi, &nw, &a[a_offset], lda, &b[ b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, & n_undeflated__, &n_deflated__, &alphar[1], &alphai[1], &beta[1], & work[1], &nw, &work[1], &nw, &work[1], &c_n1, rec, &aed_info__);
    itemp1 = (integer) work[1];
    /* Workspace query to slaqz4 */
    slaqz4_(&ilschur, &ilq, &ilz, n, ilo, ihi, &nsr, &nbr, &alphar[1], & alphai[1], &beta[1], &a[a_offset], lda, &b[b_offset], ldb, &q[ q_offset], ldq, &z__[z_offset], ldz, &work[1], &nbr, &work[1], & nbr, &work[1], &c_n1, &sweep_info__);
    itemp2 = (integer) work[1];
    /* Computing MAX */
    /* Computing 2nd power */
    i__3 = nw;
    /* Computing 2nd power */
    i__4 = nbr;
    i__1 = itemp1 + (i__3 * i__3 << 1);
    i__2 = itemp2 + (i__4 * i__4 << 1); // , expr subst
    lworkreq = max(i__1,i__2);
    if (*lwork == -1)
    {
        work[1] = (real) lworkreq;
        return 0;
    }
    else if (*lwork < lworkreq)
    {
        *info = -19;
    }
    if (*info != 0)
    {
        xerbla_("SLAQZ0", info);
        return 0;
    }
    /* Initialize Q and Z */
    if (iwantq == 3)
    {
        slaset_("FULL", n, n, &c_b25, &c_b26, &q[q_offset], ldq);
    }
    if (iwantz == 3)
    {
        slaset_("FULL", n, n, &c_b25, &c_b26, &z__[z_offset], ldz);
    }
    /* Get machine constants */
    safmin = slamch_("SAFE MINIMUM");
    safmax = 1.f / safmin;
    slabad_(&safmin, &safmax);
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real) (*n) / ulp);
    istart = *ilo;
    istop = *ihi;
    maxit = (*ihi - *ilo + 1) * 3;
    ld = 0;
    i__1 = maxit;
    for (iiter = 1;
            iiter <= i__1;
            ++iiter)
    {
        if (iiter >= maxit)
        {
            *info = istop + 1;
            goto L80;
        }
        if (istart + 1 >= istop)
        {
            istop = istart;
            break;
        }
        /* Check deflations at the end */
        /* Computing MAX */
        r__4 = smlnum;
        r__5 = ulp * ((r__1 = a[istop - 1 + (istop - 1) * a_dim1], f2c_abs(r__1)) + (r__2 = a[istop - 2 + (istop - 2) * a_dim1], f2c_abs(r__2))); // , expr subst
        if ((r__3 = a[istop - 1 + (istop - 2) * a_dim1], f2c_abs(r__3)) <= max( r__4,r__5))
        {
            a[istop - 1 + (istop - 2) * a_dim1] = 0.f;
            istop += -2;
            ld = 0;
            eshift = 0.f;
        }
        else /* if(complicated condition) */
        {
            /* Computing MAX */
            r__4 = smlnum;
            r__5 = ulp * ((r__1 = a[istop + istop * a_dim1], f2c_abs(r__1)) + (r__2 = a[istop - 1 + (istop - 1) * a_dim1], f2c_abs(r__2))); // , expr subst
            if ((r__3 = a[istop + (istop - 1) * a_dim1], f2c_abs(r__3)) <= max( r__4,r__5))
            {
                a[istop + (istop - 1) * a_dim1] = 0.f;
                --istop;
                ld = 0;
                eshift = 0.f;
            }
        }
        /* Check deflations at the start */
        /* Computing MAX */
        r__4 = smlnum;
        r__5 = ulp * ((r__1 = a[istart + 1 + (istart + 1) * a_dim1], f2c_abs(r__1)) + (r__2 = a[istart + 2 + (istart + 2) * a_dim1], f2c_abs(r__2))); // , expr subst
        if ((r__3 = a[istart + 2 + (istart + 1) * a_dim1], f2c_abs(r__3)) <= max( r__4,r__5))
        {
            a[istart + 2 + (istart + 1) * a_dim1] = 0.f;
            istart += 2;
            ld = 0;
            eshift = 0.f;
        }
        else /* if(complicated condition) */
        {
            /* Computing MAX */
            r__4 = smlnum;
            r__5 = ulp * ((r__1 = a[istart + istart * a_dim1], f2c_abs(r__1)) + (r__2 = a[istart + 1 + (istart + 1) * a_dim1], f2c_abs(r__2)));  // , expr subst
            if ((r__3 = a[istart + 1 + istart * a_dim1], f2c_abs(r__3)) <= max( r__4,r__5))
            {
                a[istart + 1 + istart * a_dim1] = 0.f;
                ++istart;
                ld = 0;
                eshift = 0.f;
            }
        }
        if (istart + 1 >= istop)
        {
            break;
        }
        /* Check interior deflations */
        istart2 = istart;
        i__2 = istart + 1;
        for (k = istop;
                k >= i__2;
                --k)
        {
            /* Computing MAX */
            r__4 = smlnum;
            r__5 = ulp * ((r__1 = a[k + k * a_dim1], f2c_abs(r__1)) + (r__2 = a[k - 1 + (k - 1) * a_dim1], f2c_abs(r__2))); // , expr subst
            if ((r__3 = a[k + (k - 1) * a_dim1], f2c_abs(r__3)) <= max(r__4,r__5))
            {
                a[k + (k - 1) * a_dim1] = 0.f;
                istart2 = k;
                break;
            }
        }
        /* Get range to apply rotations to */
        if (ilschur)
        {
            istartm = 1;
            istopm = *n;
        }
        else
        {
            istartm = istart2;
            istopm = istop;
        }
        /* Check infinite eigenvalues, this is done without blocking so might */
        /* slow down the method when many infinite eigenvalues are present */
        k = istop;
        while(k >= istart2)
        {
            temp = 0.f;
            if (k < istop)
            {
                temp += (r__1 = b[k + (k + 1) * b_dim1], f2c_abs(r__1));
            }
            if (k > istart2)
            {
                temp += (r__1 = b[k - 1 + k * b_dim1], f2c_abs(r__1));
            }
            /* Computing MAX */
            r__2 = smlnum;
            r__3 = ulp * temp; // , expr subst
            if ((r__1 = b[k + k * b_dim1], f2c_abs(r__1)) < max(r__2,r__3))
            {
                /* A diagonal element of B is negligable, move it */
                /* to the top and deflate it */
                i__2 = istart2 + 1;
                for (k2 = k;
                        k2 >= i__2;
                        --k2)
                {
                    slartg_(&b[k2 - 1 + k2 * b_dim1], &b[k2 - 1 + (k2 - 1) * b_dim1], &c1, &s1, &temp);
                    b[k2 - 1 + k2 * b_dim1] = temp;
                    b[k2 - 1 + (k2 - 1) * b_dim1] = 0.f;
                    i__3 = k2 - 2 - istartm + 1;
                    srot_(&i__3, &b[istartm + k2 * b_dim1], &c__1, &b[istartm + (k2 - 1) * b_dim1], &c__1, &c1, &s1);
                    /* Computing MIN */
                    i__4 = k2 + 1;
                    i__3 = min(i__4,istop) - istartm + 1;
                    srot_(&i__3, &a[istartm + k2 * a_dim1], &c__1, &a[istartm + (k2 - 1) * a_dim1], &c__1, &c1, &s1);
                    if (ilz)
                    {
                        srot_(n, &z__[k2 * z_dim1 + 1], &c__1, &z__[(k2 - 1) * z_dim1 + 1], &c__1, &c1, &s1);
                    }
                    if (k2 < istop)
                    {
                        slartg_(&a[k2 + (k2 - 1) * a_dim1], &a[k2 + 1 + (k2 - 1) * a_dim1], &c1, &s1, &temp);
                        a[k2 + (k2 - 1) * a_dim1] = temp;
                        a[k2 + 1 + (k2 - 1) * a_dim1] = 0.f;
                        i__3 = istopm - k2 + 1;
                        srot_(&i__3, &a[k2 + k2 * a_dim1], lda, &a[k2 + 1 + k2 * a_dim1], lda, &c1, &s1);
                        i__3 = istopm - k2 + 1;
                        srot_(&i__3, &b[k2 + k2 * b_dim1], ldb, &b[k2 + 1 + k2 * b_dim1], ldb, &c1, &s1);
                        if (ilq)
                        {
                            srot_(n, &q[k2 * q_dim1 + 1], &c__1, &q[(k2 + 1) * q_dim1 + 1], &c__1, &c1, &s1);
                        }
                    }
                }
                if (istart2 < istop)
                {
                    slartg_(&a[istart2 + istart2 * a_dim1], &a[istart2 + 1 + istart2 * a_dim1], &c1, &s1, &temp);
                    a[istart2 + istart2 * a_dim1] = temp;
                    a[istart2 + 1 + istart2 * a_dim1] = 0.f;
                    i__2 = istopm - (istart2 + 1) + 1;
                    srot_(&i__2, &a[istart2 + (istart2 + 1) * a_dim1], lda, & a[istart2 + 1 + (istart2 + 1) * a_dim1], lda, &c1, &s1);
                    i__2 = istopm - (istart2 + 1) + 1;
                    srot_(&i__2, &b[istart2 + (istart2 + 1) * b_dim1], ldb, & b[istart2 + 1 + (istart2 + 1) * b_dim1], ldb, &c1, &s1);
                    if (ilq)
                    {
                        srot_(n, &q[istart2 * q_dim1 + 1], &c__1, &q[(istart2 + 1) * q_dim1 + 1], &c__1, &c1, &s1);
                    }
                }
                ++istart2;
            }
            --k;
        }
        /* istart2 now points to the top of the bottom right */
        /* unreduced Hessenberg block */
        if (istart2 >= istop)
        {
            istop = istart2 - 1;
            ld = 0;
            eshift = 0.f;
            continue;
        }
        nw = nwr;
        nshifts = nsr;
        nblock = nbr;
        if (istop - istart2 + 1 < nmin)
        {
            /* Setting nw to the size of the subblock will make AED deflate */
            /* all the eigenvalues. This is slightly more efficient than just */
            /* using qz_small because the off diagonal part gets updated via BLAS. */
            if (istop - istart + 1 < nmin)
            {
                nw = istop - istart + 1;
                istart2 = istart;
            }
            else
            {
                nw = istop - istart2 + 1;
            }
        }
        /* Time for AED */
        /* Computing 2nd power */
        i__2 = nw;
        /* Computing 2nd power */
        i__3 = nw;
        /* Computing 2nd power */
        i__5 = nw;
        i__4 = *lwork - (i__5 * i__5 << 1);
        slaqz3_(&ilschur, &ilq, &ilz, n, &istart2, &istop, &nw, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &n_undeflated__, &n_deflated__, &alphar[1], &alphai[1], & beta[1], &work[1], &nw, &work[i__2 * i__2 + 1], &nw, &work[( i__3 * i__3 << 1) + 1], &i__4, rec, &aed_info__);
        if (n_deflated__ > 0)
        {
            istop -= n_deflated__;
            ld = 0;
            eshift = 0.f;
        }
        if (n_deflated__ * 100 > nibble * (n_deflated__ + n_undeflated__) || istop - istart2 + 1 < nmin)
        {
            /* AED has uncovered many eigenvalues. Skip a QZ sweep and run */
            /* AED again. */
            continue;
        }
        ++ld;
        /* Computing MIN */
        i__2 = nshifts;
        i__3 = istop - istart2; // , expr subst
        ns = min(i__2,i__3);
        ns = min(ns,n_undeflated__);
        shiftpos = istop - n_deflated__ - n_undeflated__ + 1;
        /* Shuffle shifts to put double shifts in front */
        /* This ensures that we don't split up a double shift */
        i__2 = shiftpos + n_undeflated__ - 1;
        for (i__ = shiftpos;
                i__ <= i__2;
                i__ += 2)
        {
            if (alphai[i__] != -alphai[i__ + 1])
            {
                swap = alphar[i__];
                alphar[i__] = alphar[i__ + 1];
                alphar[i__ + 1] = alphar[i__ + 2];
                alphar[i__ + 2] = swap;
                swap = alphai[i__];
                alphai[i__] = alphai[i__ + 1];
                alphai[i__ + 1] = alphai[i__ + 2];
                alphai[i__ + 2] = swap;
                swap = beta[i__];
                beta[i__] = beta[i__ + 1];
                beta[i__ + 1] = beta[i__ + 2];
                beta[i__ + 2] = swap;
            }
        }
        if (ld % 6 == 0)
        {
            /* Exceptional shift. Chosen for no particularly good reason. */
            if ((real) maxit * safmin * (r__1 = a[istop + (istop - 1) * a_dim1], f2c_abs(r__1)) < (r__2 = a[istop - 1 + (istop - 1) * a_dim1], f2c_abs(r__2)))
            {
                eshift = a[istop + (istop - 1) * a_dim1] / b[istop - 1 + ( istop - 1) * b_dim1];
            }
            else
            {
                eshift += 1.f / (safmin * (real) maxit);
            }
            alphar[shiftpos] = 1.f;
            alphar[shiftpos + 1] = 0.f;
            alphai[shiftpos] = 0.f;
            alphai[shiftpos + 1] = 0.f;
            beta[shiftpos] = eshift;
            beta[shiftpos + 1] = eshift;
            ns = 2;
        }
        /* Time for a QZ sweep */
        /* Computing 2nd power */
        i__2 = nblock;
        /* Computing 2nd power */
        i__3 = nblock;
        /* Computing 2nd power */
        i__5 = nblock;
        i__4 = *lwork - (i__5 * i__5 << 1);
        slaqz4_(&ilschur, &ilq, &ilz, n, &istart2, &istop, &ns, &nblock, & alphar[shiftpos], &alphai[shiftpos], &beta[shiftpos], &a[ a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[ z_offset], ldz, &work[1], &nblock, &work[i__2 * i__2 + 1], & nblock, &work[(i__3 * i__3 << 1) + 1], &i__4, &sweep_info__);
    }
    /* Call SHGEQZ to normalize the eigenvalue blocks and set the eigenvalues */
    /* If all the eigenvalues have been found, SHGEQZ will not do any iterations */
    /* and only normalize the blocks. In case of a rare convergence failure, */
    /* the single shift might perform better. */
L80:
    shgeqz_(wants, wantq, wantz, n, ilo, ihi, &a[a_offset], lda, &b[b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &q[q_offset], ldq, &z__[ z_offset], ldz, &work[1], lwork, &norm_info__);
    *info = norm_info__;
    return 0;
}
/* slaqz0_ */
