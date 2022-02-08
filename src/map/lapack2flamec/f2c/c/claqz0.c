/* claqz0.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    0.f,0.f
}
;
static complex c_b2 =
{
    1.f,0.f
}
;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__17 = 17;
static integer c_n1 = -1;
static integer c__1 = 1;
/* > \brief \b CLAQZ0 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQZ0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/CLAQZ0. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/CLAQZ0. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/CLAQZ0. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQZ0( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, */
/* $ LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, RWORK, REC, */
/* $ INFO ) */
/* IMPLICIT NONE */
/* Arguments */
/* CHARACTER, INTENT( IN ) :: WANTS, WANTQ, WANTZ */
/* INTEGER, INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, */
/* $ REC */
/* INTEGER, INTENT( OUT ) :: INFO */
/* COMPLEX, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ Z( LDZ, * ), ALPHA( * ), BETA( * ), WORK( * ) */
/* REAL, INTENT( OUT ) :: RWORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQZ0 computes the eigenvalues of a matrix pair (H,T), */
/* > where H is an upper Hessenberg matrix and T is upper triangular, */
/* > using the double-shift QZ method. */
/* > Matrix pairs of this type are produced by the reduction to */
/* > generalized upper Hessenberg form of a matrix pair (A,B): */
/* > */
/* > A = Q1*H*Z1**H, B = Q1*T*Z1**H, */
/* > */
/* > as computed by CGGHRD. */
/* > */
/* > If JOB='S', then the Hessenberg-triangular pair (H,T) is */
/* > also reduced to generalized Schur form, */
/* > */
/* > H = Q*S*Z**H, T = Q*P*Z**H, */
/* > */
/* > where Q and Z are unitary matrices, P and S are an upper triangular */
/* > matrices. */
/* > */
/* > Optionally, the unitary matrix Q from the generalized Schur */
/* > factorization may be postmultiplied into an input matrix Q1, and the */
/* > unitary matrix Z may be postmultiplied into an input matrix Z1. */
/* > If Q1 and Z1 are the unitary matrices from CGGHRD that reduced */
/* > the matrix pair (A,B) to generalized upper Hessenberg form, then the */
/* > output matrices Q1*Q and Z1*Z are the unitary factors from the */
/* > generalized Schur factorization of (A,B): */
/* > */
/* > A = (Q1*Q)*S*(Z1*Z)**H, B = (Q1*Q)*P*(Z1*Z)**H. */
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
/* > Eigenvalues can be read directly from the generalized Schur */
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
/* > = 'V': Q must contain an unitary matrix Q1 on entry and */
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
/* > = 'V': Z must contain an unitary matrix Z1 on entry and */
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
/* > A is COMPLEX array, dimension (LDA, N) */
/* > On entry, the N-by-N upper Hessenberg matrix A. */
/* > On exit, if JOB = 'S', A contains the upper triangular */
/* > matrix S from the generalized Schur factorization. */
/* > If JOB = 'E', the diagonal of A matches that of S, but */
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
/* > B is COMPLEX array, dimension (LDB, N) */
/* > On entry, the N-by-N upper triangular matrix B. */
/* > On exit, if JOB = 'S', B contains the upper triangular */
/* > matrix P from the generalized Schur factorization. */
/* > If JOB = 'E', the diagonal of B matches that of P, but */
/* > the rest of B is unspecified. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* > ALPHA is COMPLEX array, dimension (N) */
/* > Each scalar alpha defining an eigenvalue */
/* > of GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is COMPLEX array, dimension (N) */
/* > The scalars beta that define the eigenvalues of GNEP. */
/* > Together, the quantities alpha = ALPHA(j) and */
/* > beta = BETA(j) represent the j-th eigenvalue of the matrix */
/* > pair (A,B), in one of the forms lambda = alpha/beta or */
/* > mu = beta/alpha. Since either lambda or mu may overflow, */
/* > they should not, in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDQ, N) */
/* > On entry, if COMPQ = 'V', the unitary matrix Q1 used in */
/* > the reduction of (A,B) to generalized Hessenberg form. */
/* > On exit, if COMPQ = 'I', the unitary matrix of left Schur */
/* > vectors of (A,B), and if COMPQ = 'V', the unitary matrix */
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
/* > Z is COMPLEX array, dimension (LDZ, N) */
/* > On entry, if COMPZ = 'V', the unitary matrix Z1 used in */
/* > the reduction of (A,B) to generalized Hessenberg form. */
/* > On exit, if COMPZ = 'I', the unitary matrix of */
/* > right Schur vectors of (H,T), and if COMPZ = 'V', the */
/* > unitary matrix of right Schur vectors of (A,B). */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (N) */
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
/* > in Schur form, but ALPHA(i) and */
/* > BETA(i), i=INFO+1,...,N should be correct. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Thijs Steel, KU Leuven */
/* > \date May 2020 */
/* > \ingroup complexGEcomputational */
/* > */
/* ===================================================================== */
/* Subroutine */
int claqz0_(char *wants, char *wantq, char *wantz, integer * n, integer *ilo, integer *ihi, complex *a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta, complex *q, integer *ldq, complex *z__, integer *ldz, complex *work, integer *lwork, real * rwork, integer *rec, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2;
    complex q__1, q__2;
    /* Builtin functions */
    double sqrt(doublereal), c_abs(complex *);
    void r_cnjg(complex *, complex *), c_div(complex *, complex *, complex *);
    /* Local variables */
    integer aed_info__;
    extern /* Subroutine */
    int f90_cycle_(void);
    integer shiftpos, lworkreq, k;
    real c1;
    integer k2;
    complex s1;
    integer norm_info__, ld, ns, n_deflated__, nw, sweep_info__, nbr;
    logical ilq, ilz;
    real ulp;
    integer nsr, nwr, nmin;
    complex temp;
    extern /* Subroutine */
    int crot_(integer *, complex *, integer *, complex *, integer *, real *, complex *);
    integer n_undeflated__;
    extern logical lsame_(char *, char *);
    integer iiter, maxit;
    real tempr;
    integer rcost, istop;
    extern /* Subroutine */
    int claqz2_(logical *, logical *, logical *, integer *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, integer *, integer *, complex *, complex *, complex *, integer *, complex *, integer *, complex *, integer *, real *, integer *, integer *), claqz3_(logical *, logical *, logical *, integer *, integer *, integer *, integer *, integer *, complex *, complex *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, integer *);
    integer itemp1, itemp2;
    extern /* Subroutine */
    int slabad_(real *, real *);
    integer nibble;
    extern real slamch_(char *);
    integer nblock;
    extern /* Subroutine */
    int claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *), clartg_(complex *, complex *, real *, complex *, complex *);
    real safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real safmax;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int chgeqz_(char *, char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, complex *, complex *, integer *, complex *, integer *, complex *, integer *, real *, integer *);
    complex eshift;
    char jbcmpz[3];
    integer iwantq, iwants, istart;
    real smlnum;
    integer istopm, iwantz, istart2;
    extern /* Subroutine */
    int f90_exit_(void);
    logical ilschur;
    integer nshifts, istartm;
    /* Arguments */
    /* const for wants,wantq, wantz, n, ilo, ihi, lda, ldb, ldq, ldz, lwork,rec */
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
    --alpha;
    --beta;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
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
        xerbla_("CLAQZ0", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n <= 0)
    {
        work[1].r = 1.f;
        work[1].i = 0.f; // , expr subst
        return 0;
    }
    /* Get the parameters */
    *(unsigned char *)jbcmpz = *(unsigned char *)wants;
    *(unsigned char *)&jbcmpz[1] = *(unsigned char *)wantq;
    *(unsigned char *)&jbcmpz[2] = *(unsigned char *)wantz;
    nmin = ilaenv_(&c__12, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    nwr = ilaenv_(&c__13, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    nwr = max(2,nwr);
    /* Computing MIN */
    i__1 = *ihi - *ilo + 1;
    i__2 = (*n - 1) / 3;
    i__1 = min(i__1,i__2); // ; expr subst
    nwr = min(i__1,nwr);
    nibble = ilaenv_(&c__14, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    nsr = ilaenv_(&c__15, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    /* Computing MIN */
    i__1 = nsr, i__2 = (*n + 6) / 9;
    i__1 = min(i__1,i__2);
    i__2 = *ihi - * ilo; // ; expr subst
    nsr = min(i__1,i__2);
    /* Computing MAX */
    i__1 = 2;
    i__2 = nsr - nsr % 2; // , expr subst
    nsr = max(i__1,i__2);
    rcost = ilaenv_(&c__17, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork);
    itemp1 = (integer) (nsr / sqrt((nsr << 1) / ((real) rcost / 100 * *n) + 1) );
    itemp1 = ((itemp1 - 1) / 4 << 2) + 4;
    nbr = nsr + itemp1;
    if (*n < nmin || *rec >= 2)
    {
        chgeqz_(wants, wantq, wantz, n, ilo, ihi, &a[a_offset], lda, &b[ b_offset], ldb, &alpha[1], &beta[1], &q[q_offset], ldq, &z__[ z_offset], ldz, &work[1], lwork, &rwork[1], info);
        return 0;
    }
    /* Find out required workspace */
    /* Workspace query to CLAQZ2 */
    nw = max(nwr,nmin);
    claqz2_(&ilschur, &ilq, &ilz, n, ilo, ihi, &nw, &a[a_offset], lda, &b[ b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, & n_undeflated__, &n_deflated__, &alpha[1], &beta[1], &work[1], &nw, &work[1], &nw, &work[1], &c_n1, &rwork[1], rec, &aed_info__);
    itemp1 = (integer) work[1].r;
    /* Workspace query to CLAQZ3 */
    claqz3_(&ilschur, &ilq, &ilz, n, ilo, ihi, &nsr, &nbr, &alpha[1], &beta[1], &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[ z_offset], ldz, &work[1], &nbr, &work[1], &nbr, &work[1], &c_n1, & sweep_info__);
    itemp2 = (integer) work[1].r;
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
        r__1 = (real) lworkreq;
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
        return 0;
    }
    else if (*lwork < lworkreq)
    {
        *info = -19;
    }
    if (*info != 0)
    {
        xerbla_("CLAQZ0", info);
        return 0;
    }
    /* Initialize Q and Z */
    if (iwantq == 3)
    {
        claset_("FULL", n, n, &c_b1, &c_b2, &q[q_offset], ldq);
    }
    if (iwantz == 3)
    {
        claset_("FULL", n, n, &c_b1, &c_b2, &z__[z_offset], ldz);
    }
    /* Get machine constants */
    safmin = slamch_("SAFE MINIMUM");
    safmax = 1.f / safmin;
    slabad_(&safmin, &safmax);
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real) (*n) / ulp);
    istart = *ilo;
    istop = *ihi;
    maxit = (*ihi - *ilo + 1) * 30;
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
        r__1 = smlnum;
        r__2 = ulp * (c_abs(&a[istop + istop * a_dim1]) + c_abs(&a[istop - 1 + (istop - 1) * a_dim1])); // , expr subst
        if (c_abs(&a[istop + (istop - 1) * a_dim1]) <= max(r__1,r__2))
        {
            i__2 = istop + (istop - 1) * a_dim1;
            a[i__2].r = 0.f;
            a[i__2].i = 0.f; // , expr subst
            --istop;
            ld = 0;
            eshift.r = 0.f;
            eshift.i = 0.f; // , expr subst
        }
        /* Check deflations at the start */
        /* Computing MAX */
        r__1 = smlnum;
        r__2 = ulp * (c_abs(&a[istart + istart * a_dim1]) + c_abs(&a[istart + 1 + (istart + 1) * a_dim1])); // , expr subst
        if (c_abs(&a[istart + 1 + istart * a_dim1]) <= max(r__1,r__2))
        {
            i__2 = istart + 1 + istart * a_dim1;
            a[i__2].r = 0.f;
            a[i__2].i = 0.f; // , expr subst
            ++istart;
            ld = 0;
            eshift.r = 0.f;
            eshift.i = 0.f; // , expr subst
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
            r__1 = smlnum;
            r__2 = ulp * (c_abs(&a[k + k * a_dim1]) + c_abs(&a[ k - 1 + (k - 1) * a_dim1])); // , expr subst
            if (c_abs(&a[k + (k - 1) * a_dim1]) <= max(r__1,r__2))
            {
                i__3 = k + (k - 1) * a_dim1;
                a[i__3].r = 0.f;
                a[i__3].i = 0.f; // , expr subst
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
            tempr = 0.f;
            if (k < istop)
            {
                tempr += c_abs(&b[k + (k + 1) * b_dim1]);
            }
            if (k > istart2)
            {
                tempr += c_abs(&b[k - 1 + k * b_dim1]);
            }
            /* Computing MAX */
            r__1 = smlnum;
            r__2 = ulp * tempr; // , expr subst
            if (c_abs(&b[k + k * b_dim1]) < max(r__1,r__2))
            {
                /* A diagonal element of B is negligable, move it */
                /* to the top and deflate it */
                i__2 = istart2 + 1;
                for (k2 = k;
                        k2 >= i__2;
                        --k2)
                {
                    clartg_(&b[k2 - 1 + k2 * b_dim1], &b[k2 - 1 + (k2 - 1) * b_dim1], &c1, &s1, &temp);
                    i__3 = k2 - 1 + k2 * b_dim1;
                    b[i__3].r = temp.r;
                    b[i__3].i = temp.i; // , expr subst
                    i__3 = k2 - 1 + (k2 - 1) * b_dim1;
                    b[i__3].r = 0.f;
                    b[i__3].i = 0.f; // , expr subst
                    i__3 = k2 - 2 - istartm + 1;
                    crot_(&i__3, &b[istartm + k2 * b_dim1], &c__1, &b[istartm + (k2 - 1) * b_dim1], &c__1, &c1, &s1);
                    /* Computing MIN */
                    i__4 = k2 + 1;
                    i__3 = min(i__4,istop) - istartm + 1;
                    crot_(&i__3, &a[istartm + k2 * a_dim1], &c__1, &a[istartm + (k2 - 1) * a_dim1], &c__1, &c1, &s1);
                    if (ilz)
                    {
                        crot_(n, &z__[k2 * z_dim1 + 1], &c__1, &z__[(k2 - 1) * z_dim1 + 1], &c__1, &c1, &s1);
                    }
                    if (k2 < istop)
                    {
                        clartg_(&a[k2 + (k2 - 1) * a_dim1], &a[k2 + 1 + (k2 - 1) * a_dim1], &c1, &s1, &temp);
                        i__3 = k2 + (k2 - 1) * a_dim1;
                        a[i__3].r = temp.r;
                        a[i__3].i = temp.i; // , expr subst
                        i__3 = k2 + 1 + (k2 - 1) * a_dim1;
                        a[i__3].r = 0.f;
                        a[i__3].i = 0.f; // , expr subst
                        i__3 = istopm - k2 + 1;
                        crot_(&i__3, &a[k2 + k2 * a_dim1], lda, &a[k2 + 1 + k2 * a_dim1], lda, &c1, &s1);
                        i__3 = istopm - k2 + 1;
                        crot_(&i__3, &b[k2 + k2 * b_dim1], ldb, &b[k2 + 1 + k2 * b_dim1], ldb, &c1, &s1);
                        if (ilq)
                        {
                            r_cnjg(&q__1, &s1);
                            crot_(n, &q[k2 * q_dim1 + 1], &c__1, &q[(k2 + 1) * q_dim1 + 1], &c__1, &c1, &q__1);
                        }
                    }
                }
                if (istart2 < istop)
                {
                    clartg_(&a[istart2 + istart2 * a_dim1], &a[istart2 + 1 + istart2 * a_dim1], &c1, &s1, &temp);
                    i__2 = istart2 + istart2 * a_dim1;
                    a[i__2].r = temp.r;
                    a[i__2].i = temp.i; // , expr subst
                    i__2 = istart2 + 1 + istart2 * a_dim1;
                    a[i__2].r = 0.f;
                    a[i__2].i = 0.f; // , expr subst
                    i__2 = istopm - (istart2 + 1) + 1;
                    crot_(&i__2, &a[istart2 + (istart2 + 1) * a_dim1], lda, & a[istart2 + 1 + (istart2 + 1) * a_dim1], lda, &c1, &s1);
                    i__2 = istopm - (istart2 + 1) + 1;
                    crot_(&i__2, &b[istart2 + (istart2 + 1) * b_dim1], ldb, & b[istart2 + 1 + (istart2 + 1) * b_dim1], ldb, &c1, &s1);
                    if (ilq)
                    {
                        r_cnjg(&q__1, &s1);
                        crot_(n, &q[istart2 * q_dim1 + 1], &c__1, &q[(istart2 + 1) * q_dim1 + 1], &c__1, &c1, &q__1);
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
            eshift.r = 0.f;
            eshift.i = 0.f; // , expr subst
            continue;
        }
        nw = nwr;
        nshifts = nsr;
        nblock = nbr;
        if (istop - istart2 + 1 < nmin)
        {
            /* Setting nw to the size of the subblock will make AED deflate */
            /* all the eigenvalues. This is slightly more efficient than just */
            /* using CHGEQZ because the off diagonal part gets updated via BLAS. */
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
        claqz2_(&ilschur, &ilq, &ilz, n, &istart2, &istop, &nw, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &n_undeflated__, &n_deflated__, &alpha[1], &beta[1], & work[1], &nw, &work[i__2 * i__2 + 1], &nw, &work[(i__3 * i__3 << 1) + 1], &i__4, &rwork[1], rec, &aed_info__);
        if (n_deflated__ > 0)
        {
            istop -= n_deflated__;
            ld = 0;
            eshift.r = 0.f;
            eshift.i = 0.f; // , expr subst
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
        if (ld % 6 == 0)
        {
            /* Exceptional shift. Chosen for no particularly good reason. */
            if ((real) maxit * safmin * c_abs(&a[istop + (istop - 1) * a_dim1] ) < c_abs(&a[istop - 1 + (istop - 1) * a_dim1]))
            {
                c_div(&q__1, &a[istop + (istop - 1) * a_dim1], &b[istop - 1 + (istop - 1) * b_dim1]);
                eshift.r = q__1.r;
                eshift.i = q__1.i; // , expr subst
            }
            else
            {
                r__1 = safmin * (real) maxit;
                q__2.r = 1.f / r__1;
                q__2.i = 0.f / r__1; // , expr subst
                q__1.r = eshift.r + q__2.r;
                q__1.i = eshift.i + q__2.i; // , expr subst
                eshift.r = q__1.r;
                eshift.i = q__1.i; // , expr subst
            }
            i__2 = shiftpos;
            alpha[i__2].r = 1.f;
            alpha[i__2].i = 0.f; // , expr subst
            i__2 = shiftpos;
            beta[i__2].r = eshift.r;
            beta[i__2].i = eshift.i; // , expr subst
            ns = 1;
        }
        /* Time for a QZ sweep */
        /* Computing 2nd power */
        i__2 = nblock;
        /* Computing 2nd power */
        i__3 = nblock;
        /* Computing 2nd power */
        i__5 = nblock;
        i__4 = *lwork - (i__5 * i__5 << 1);
        claqz3_(&ilschur, &ilq, &ilz, n, &istart2, &istop, &ns, &nblock, & alpha[shiftpos], &beta[shiftpos], &a[a_offset], lda, &b[ b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &work[ 1], &nblock, &work[i__2 * i__2 + 1], &nblock, &work[(i__3 * i__3 << 1) + 1], &i__4, &sweep_info__);
    }
    /* Call CHGEQZ to normalize the eigenvalue blocks and set the eigenvalues */
    /* If all the eigenvalues have been found, CHGEQZ will not do any iterations */
    /* and only normalize the blocks. In case of a rare convergence failure, */
    /* the single shift might perform better. */
L80:
    chgeqz_(wants, wantq, wantz, n, ilo, ihi, &a[a_offset], lda, &b[b_offset], ldb, &alpha[1], &beta[1], &q[q_offset], ldq, &z__[z_offset], ldz, &work[1], lwork, &rwork[1], &norm_info__);
    *info = norm_info__;
    return 0;
}
/* claqz0_ */

