/* ../netlib/stgex2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__4 = 4;
static real c_b5 = 0.f;
static integer c__1 = 1;
static integer c__2 = 2;
static real c_b42 = 1.f;
static real c_b48 = -1.f;
static integer c__0 = 0;
/* > \brief \b STGEX2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an orthogon al equivalence transformation. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STGEX2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgex2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgex2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgex2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, */
/* LDZ, J1, N1, N2, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* LOGICAL WANTQ, WANTZ */
/* INTEGER INFO, J1, LDA, LDB, LDQ, LDZ, LWORK, N, N1, N2 */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGEX2 swaps adjacent diagonal blocks (A11, B11) and (A22, B22) */
/* > of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair */
/* > (A, B) by an orthogonal equivalence transformation. */
/* > */
/* > (A, B) must be in generalized real Schur canonical form (as returned */
/* > by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2 */
/* > diagonal blocks. B is upper triangular. */
/* > */
/* > Optionally, the matrices Q and Z of generalized Schur vectors are */
/* > updated. */
/* > */
/* > Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T */
/* > Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTQ */
/* > \verbatim */
/* > WANTQ is LOGICAL */
/* > .TRUE. : update the left transformation matrix Q;
*/
/* > .FALSE.: do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is LOGICAL */
/* > .TRUE. : update the right transformation matrix Z;
*/
/* > .FALSE.: do not update Z. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL arrays, dimensions (LDA,N) */
/* > On entry, the matrix A in the pair (A, B). */
/* > On exit, the updated matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL arrays, dimensions (LDB,N) */
/* > On entry, the matrix B in the pair (A, B). */
/* > On exit, the updated matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDZ,N) */
/* > On entry, if WANTQ = .TRUE., the orthogonal matrix Q. */
/* > On exit, the updated matrix Q. */
/* > Not referenced if WANTQ = .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= 1. */
/* > If WANTQ = .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ,N) */
/* > On entry, if WANTZ =.TRUE., the orthogonal matrix Z. */
/* > On exit, the updated matrix Z. */
/* > Not referenced if WANTZ = .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1. */
/* > If WANTZ = .TRUE., LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] J1 */
/* > \verbatim */
/* > J1 is INTEGER */
/* > The index to the first block (A11, B11). 1 <= J1 <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* > N1 is INTEGER */
/* > The order of the first block (A11, B11). N1 = 0, 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[in] N2 */
/* > \verbatim */
/* > N2 is INTEGER */
/* > The order of the second block (A22, B22). N2 = 0, 1 or 2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)). */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > LWORK >= MAX( N*(N2+N1), (N2+N1)*(N2+N1)*2 ) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > =0: Successful exit */
/* > >0: If INFO = 1, the transformed matrix (A, B) would be */
/* > too far from generalized Schur form;
the blocks are */
/* > not swapped and (A, B) and (Q, Z) are unchanged. */
/* > The problem of swapping is too ill-conditioned. */
/* > <0: If INFO = -16: LWORK is too small. Appropriate value */
/* > for LWORK is returned in WORK(1). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realGEauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > In the current code both weak and strong stability tests are */
/* > performed. The user can omit the strong stability test by changing */
/* > the internal logical parameter WANDS to .FALSE.. See ref. [2] for */
/* > details. */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* > \par References: */
/* ================ */
/* > */
/* > \verbatim */
/* > */
/* > [1] B. Kagstrom;
A Direct Method for Reordering Eigenvalues in the */
/* > Generalized Real Schur Form of a Regular Matrix Pair (A, B), in */
/* > M.S. Moonen et al (eds), Linear Algebra for Large Scale and */
/* > Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218. */
/* > */
/* > [2] B. Kagstrom and P. Poromaa;
Computing Eigenspaces with Specified */
/* > Eigenvalues of a Regular Matrix Pair (A, B) and Condition */
/* > Estimation: Theory, Algorithms and Software, */
/* > Report UMINF - 94.04, Department of Computing Science, Umea */
/* > University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working */
/* > Note 87. To appear in Numerical Algorithms, 1996. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int stgex2_(logical *wantq, logical *wantz, integer *n, real *a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real * z__, integer *ldz, integer *j1, integer *n1, integer *n2, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real f, g;
    integer i__, m;
    real s[16] /* was [4][4] */
    , t[16] /* was [4][4] */
    , be[2], ai[2], ar[2], sa, sb, li[16] /* was [4][4] */
    , ir[16] /* was [4][4] */
    , ss, ws, eps;
    logical weak;
    real ddum;
    integer idum;
    real taul[4], dsum, taur[4], scpy[16] /* was [4][4] */
    , tcpy[16] /* was [4][4] */
    ;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *);
    real scale, bqra21, brqa21;
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    real licop[16] /* was [4][4] */
    ;
    integer linfo;
    extern /* Subroutine */
    int sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    real ircop[16] /* was [4][4] */
    , dnorm;
    integer iwork[4];
    extern /* Subroutine */
    int slagv2_(real *, integer *, real *, integer *, real *, real *, real *, real *, real *, real *, real *), sgeqr2_( integer *, integer *, real *, integer *, real *, real *, integer * ), sgerq2_(integer *, integer *, real *, integer *, real *, real * , integer *), sorg2r_(integer *, integer *, integer *, real *, integer *, real *, real *, integer *), sorgr2_(integer *, integer *, integer *, real *, integer *, real *, real *, integer *), sorm2r_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *), sormr2_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *);
    real dscale;
    extern /* Subroutine */
    int stgsy2_(char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, real * , integer *, real *, integer *, real *, integer *, real *, real *, real *, integer *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *), slartg_(real *, real *, real *, real *, real *);
    real thresh;
    extern /* Subroutine */
    int slaset_(char *, integer *, integer *, real *, real *, real *, integer *), slassq_(integer *, real *, integer *, real *, real *);
    real smlnum;
    logical strong;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* Replaced various illegal calls to SCOPY by calls to SLASET, or by DO */
    /* loops. Sven Hammarling, 1/5/02. */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
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
    --work;
    /* Function Body */
    *info = 0;
    /* Quick return if possible */
    if (*n <= 1 || *n1 <= 0 || *n2 <= 0)
    {
        return 0;
    }
    if (*n1 > *n || *j1 + *n1 > *n)
    {
        return 0;
    }
    m = *n1 + *n2;
    /* Computing MAX */
    i__1 = *n * m;
    i__2 = m * m << 1; // , expr subst
    if (*lwork < max(i__1,i__2))
    {
        *info = -16;
        /* Computing MAX */
        i__1 = *n * m;
        i__2 = m * m << 1; // , expr subst
        work[1] = (real) max(i__1,i__2);
        return 0;
    }
    weak = FALSE_;
    strong = FALSE_;
    /* Make a local copy of selected block */
    slaset_("Full", &c__4, &c__4, &c_b5, &c_b5, li, &c__4);
    slaset_("Full", &c__4, &c__4, &c_b5, &c_b5, ir, &c__4);
    slacpy_("Full", &m, &m, &a[*j1 + *j1 * a_dim1], lda, s, &c__4);
    slacpy_("Full", &m, &m, &b[*j1 + *j1 * b_dim1], ldb, t, &c__4);
    /* Compute threshold for testing acceptance of swapping. */
    eps = slamch_("P");
    smlnum = slamch_("S") / eps;
    dscale = 0.f;
    dsum = 1.f;
    slacpy_("Full", &m, &m, s, &c__4, &work[1], &m);
    i__1 = m * m;
    slassq_(&i__1, &work[1], &c__1, &dscale, &dsum);
    slacpy_("Full", &m, &m, t, &c__4, &work[1], &m);
    i__1 = m * m;
    slassq_(&i__1, &work[1], &c__1, &dscale, &dsum);
    dnorm = dscale * sqrt(dsum);
    /* THRES has been changed from */
    /* THRESH = MAX( TEN*EPS*SA, SMLNUM ) */
    /* to */
    /* THRESH = MAX( TWENTY*EPS*SA, SMLNUM ) */
    /* on 04/01/10. */
    /* "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by */
    /* Jim Demmel and Guillaume Revy. See forum post 1783. */
    /* Computing MAX */
    r__1 = eps * 20.f * dnorm;
    thresh = max(r__1,smlnum);
    if (m == 2)
    {
        /* CASE 1: Swap 1-by-1 and 1-by-1 blocks. */
        /* Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks */
        /* using Givens rotations and perform the swap tentatively. */
        f = s[5] * t[0] - t[5] * s[0];
        g = s[5] * t[4] - t[5] * s[4];
        sb = f2c_abs(t[5]);
        sa = f2c_abs(s[5]);
        slartg_(&f, &g, &ir[4], ir, &ddum);
        ir[1] = -ir[4];
        ir[5] = ir[0];
        srot_(&c__2, s, &c__1, &s[4], &c__1, ir, &ir[1]);
        srot_(&c__2, t, &c__1, &t[4], &c__1, ir, &ir[1]);
        if (sa >= sb)
        {
            slartg_(s, &s[1], li, &li[1], &ddum);
        }
        else
        {
            slartg_(t, &t[1], li, &li[1], &ddum);
        }
        srot_(&c__2, s, &c__4, &s[1], &c__4, li, &li[1]);
        srot_(&c__2, t, &c__4, &t[1], &c__4, li, &li[1]);
        li[5] = li[0];
        li[4] = -li[1];
        /* Weak stability test: */
        /* |S21| + |T21| <= O(EPS * F-norm((S, T))) */
        ws = f2c_abs(s[1]) + f2c_abs(t[1]);
        weak = ws <= thresh;
        if (! weak)
        {
            goto L70;
        }
        if (TRUE_)
        {
            /* Strong stability test: */
            /* F-norm((A-QL**T*S*QR, B-QL**T*T*QR)) <= O(EPS*F-norm((A, B))) */
            slacpy_("Full", &m, &m, &a[*j1 + *j1 * a_dim1], lda, &work[m * m + 1], &m);
            sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, s, &c__4, &c_b5, & work[1], &m);
            sgemm_("N", "T", &m, &m, &m, &c_b48, &work[1], &m, ir, &c__4, & c_b42, &work[m * m + 1], &m);
            dscale = 0.f;
            dsum = 1.f;
            i__1 = m * m;
            slassq_(&i__1, &work[m * m + 1], &c__1, &dscale, &dsum);
            slacpy_("Full", &m, &m, &b[*j1 + *j1 * b_dim1], ldb, &work[m * m + 1], &m);
            sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, t, &c__4, &c_b5, & work[1], &m);
            sgemm_("N", "T", &m, &m, &m, &c_b48, &work[1], &m, ir, &c__4, & c_b42, &work[m * m + 1], &m);
            i__1 = m * m;
            slassq_(&i__1, &work[m * m + 1], &c__1, &dscale, &dsum);
            ss = dscale * sqrt(dsum);
            strong = ss <= thresh;
            if (! strong)
            {
                goto L70;
            }
        }
        /* Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and */
        /* (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)). */
        i__1 = *j1 + 1;
        srot_(&i__1, &a[*j1 * a_dim1 + 1], &c__1, &a[(*j1 + 1) * a_dim1 + 1], &c__1, ir, &ir[1]);
        i__1 = *j1 + 1;
        srot_(&i__1, &b[*j1 * b_dim1 + 1], &c__1, &b[(*j1 + 1) * b_dim1 + 1], &c__1, ir, &ir[1]);
        i__1 = *n - *j1 + 1;
        srot_(&i__1, &a[*j1 + *j1 * a_dim1], lda, &a[*j1 + 1 + *j1 * a_dim1], lda, li, &li[1]);
        i__1 = *n - *j1 + 1;
        srot_(&i__1, &b[*j1 + *j1 * b_dim1], ldb, &b[*j1 + 1 + *j1 * b_dim1], ldb, li, &li[1]);
        /* Set N1-by-N2 (2,1) - blocks to ZERO. */
        a[*j1 + 1 + *j1 * a_dim1] = 0.f;
        b[*j1 + 1 + *j1 * b_dim1] = 0.f;
        /* Accumulate transformations into Q and Z if requested. */
        if (*wantz)
        {
            srot_(n, &z__[*j1 * z_dim1 + 1], &c__1, &z__[(*j1 + 1) * z_dim1 + 1], &c__1, ir, &ir[1]);
        }
        if (*wantq)
        {
            srot_(n, &q[*j1 * q_dim1 + 1], &c__1, &q[(*j1 + 1) * q_dim1 + 1], &c__1, li, &li[1]);
        }
        /* Exit with INFO = 0 if swap was successfully performed. */
        return 0;
    }
    else
    {
        /* CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2 */
        /* and 2-by-2 blocks. */
        /* Solve the generalized Sylvester equation */
        /* S11 * R - L * S22 = SCALE * S12 */
        /* T11 * R - L * T22 = SCALE * T12 */
        /* for R and L. Solutions in LI and IR. */
        slacpy_("Full", n1, n2, &t[(*n1 + 1 << 2) - 4], &c__4, li, &c__4);
        slacpy_("Full", n1, n2, &s[(*n1 + 1 << 2) - 4], &c__4, &ir[*n2 + 1 + ( *n1 + 1 << 2) - 5], &c__4);
        stgsy2_("N", &c__0, n1, n2, s, &c__4, &s[*n1 + 1 + (*n1 + 1 << 2) - 5] , &c__4, &ir[*n2 + 1 + (*n1 + 1 << 2) - 5], &c__4, t, &c__4, & t[*n1 + 1 + (*n1 + 1 << 2) - 5], &c__4, li, &c__4, &scale, & dsum, &dscale, iwork, &idum, &linfo);
        /* Compute orthogonal matrix QL: */
        /* QL**T * LI = [ TL ] */
        /* [ 0 ] */
        /* where */
        /* LI = [ -L ] */
        /* [ SCALE * identity(N2) ] */
        i__1 = *n2;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            sscal_(n1, &c_b48, &li[(i__ << 2) - 4], &c__1);
            li[*n1 + i__ + (i__ << 2) - 5] = scale;
            /* L10: */
        }
        sgeqr2_(&m, n2, li, &c__4, taul, &work[1], &linfo);
        if (linfo != 0)
        {
            goto L70;
        }
        sorg2r_(&m, &m, n2, li, &c__4, taul, &work[1], &linfo);
        if (linfo != 0)
        {
            goto L70;
        }
        /* Compute orthogonal matrix RQ: */
        /* IR * RQ**T = [ 0 TR], */
        /* where IR = [ SCALE * identity(N1), R ] */
        i__1 = *n1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            ir[*n2 + i__ + (i__ << 2) - 5] = scale;
            /* L20: */
        }
        sgerq2_(n1, &m, &ir[*n2], &c__4, taur, &work[1], &linfo);
        if (linfo != 0)
        {
            goto L70;
        }
        sorgr2_(&m, &m, n1, ir, &c__4, taur, &work[1], &linfo);
        if (linfo != 0)
        {
            goto L70;
        }
        /* Perform the swapping tentatively: */
        sgemm_("T", "N", &m, &m, &m, &c_b42, li, &c__4, s, &c__4, &c_b5, & work[1], &m);
        sgemm_("N", "T", &m, &m, &m, &c_b42, &work[1], &m, ir, &c__4, &c_b5, s, &c__4);
        sgemm_("T", "N", &m, &m, &m, &c_b42, li, &c__4, t, &c__4, &c_b5, & work[1], &m);
        sgemm_("N", "T", &m, &m, &m, &c_b42, &work[1], &m, ir, &c__4, &c_b5, t, &c__4);
        slacpy_("F", &m, &m, s, &c__4, scpy, &c__4);
        slacpy_("F", &m, &m, t, &c__4, tcpy, &c__4);
        slacpy_("F", &m, &m, ir, &c__4, ircop, &c__4);
        slacpy_("F", &m, &m, li, &c__4, licop, &c__4);
        /* Triangularize the B-part by an RQ factorization. */
        /* Apply transformation (from left) to A-part, giving S. */
        sgerq2_(&m, &m, t, &c__4, taur, &work[1], &linfo);
        if (linfo != 0)
        {
            goto L70;
        }
        sormr2_("R", "T", &m, &m, &m, t, &c__4, taur, s, &c__4, &work[1], & linfo);
        if (linfo != 0)
        {
            goto L70;
        }
        sormr2_("L", "N", &m, &m, &m, t, &c__4, taur, ir, &c__4, &work[1], & linfo);
        if (linfo != 0)
        {
            goto L70;
        }
        /* Compute F-norm(S21) in BRQA21. (T21 is 0.) */
        dscale = 0.f;
        dsum = 1.f;
        i__1 = *n2;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            slassq_(n1, &s[*n2 + 1 + (i__ << 2) - 5], &c__1, &dscale, &dsum);
            /* L30: */
        }
        brqa21 = dscale * sqrt(dsum);
        /* Triangularize the B-part by a QR factorization. */
        /* Apply transformation (from right) to A-part, giving S. */
        sgeqr2_(&m, &m, tcpy, &c__4, taul, &work[1], &linfo);
        if (linfo != 0)
        {
            goto L70;
        }
        sorm2r_("L", "T", &m, &m, &m, tcpy, &c__4, taul, scpy, &c__4, &work[1] , info);
        sorm2r_("R", "N", &m, &m, &m, tcpy, &c__4, taul, licop, &c__4, &work[ 1], info);
        if (linfo != 0)
        {
            goto L70;
        }
        /* Compute F-norm(S21) in BQRA21. (T21 is 0.) */
        dscale = 0.f;
        dsum = 1.f;
        i__1 = *n2;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            slassq_(n1, &scpy[*n2 + 1 + (i__ << 2) - 5], &c__1, &dscale, & dsum);
            /* L40: */
        }
        bqra21 = dscale * sqrt(dsum);
        /* Decide which method to use. */
        /* Weak stability test: */
        /* F-norm(S21) <= O(EPS * F-norm((S, T))) */
        if (bqra21 <= brqa21 && bqra21 <= thresh)
        {
            slacpy_("F", &m, &m, scpy, &c__4, s, &c__4);
            slacpy_("F", &m, &m, tcpy, &c__4, t, &c__4);
            slacpy_("F", &m, &m, ircop, &c__4, ir, &c__4);
            slacpy_("F", &m, &m, licop, &c__4, li, &c__4);
        }
        else if (brqa21 >= thresh)
        {
            goto L70;
        }
        /* Set lower triangle of B-part to zero */
        i__1 = m - 1;
        i__2 = m - 1;
        slaset_("Lower", &i__1, &i__2, &c_b5, &c_b5, &t[1], &c__4);
        if (TRUE_)
        {
            /* Strong stability test: */
            /* F-norm((A-QL*S*QR**T, B-QL*T*QR**T)) <= O(EPS*F-norm((A,B))) */
            slacpy_("Full", &m, &m, &a[*j1 + *j1 * a_dim1], lda, &work[m * m + 1], &m);
            sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, s, &c__4, &c_b5, & work[1], &m);
            sgemm_("N", "N", &m, &m, &m, &c_b48, &work[1], &m, ir, &c__4, & c_b42, &work[m * m + 1], &m);
            dscale = 0.f;
            dsum = 1.f;
            i__1 = m * m;
            slassq_(&i__1, &work[m * m + 1], &c__1, &dscale, &dsum);
            slacpy_("Full", &m, &m, &b[*j1 + *j1 * b_dim1], ldb, &work[m * m + 1], &m);
            sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, t, &c__4, &c_b5, & work[1], &m);
            sgemm_("N", "N", &m, &m, &m, &c_b48, &work[1], &m, ir, &c__4, & c_b42, &work[m * m + 1], &m);
            i__1 = m * m;
            slassq_(&i__1, &work[m * m + 1], &c__1, &dscale, &dsum);
            ss = dscale * sqrt(dsum);
            strong = ss <= thresh;
            if (! strong)
            {
                goto L70;
            }
        }
        /* If the swap is accepted ("weakly" and "strongly"), apply the */
        /* transformations and set N1-by-N2 (2,1)-block to zero. */
        slaset_("Full", n1, n2, &c_b5, &c_b5, &s[*n2], &c__4);
        /* copy back M-by-M diagonal block starting at index J1 of (A, B) */
        slacpy_("F", &m, &m, s, &c__4, &a[*j1 + *j1 * a_dim1], lda) ;
        slacpy_("F", &m, &m, t, &c__4, &b[*j1 + *j1 * b_dim1], ldb) ;
        slaset_("Full", &c__4, &c__4, &c_b5, &c_b5, t, &c__4);
        /* Standardize existing 2-by-2 blocks. */
        i__1 = m * m;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            work[i__] = 0.f;
            /* L50: */
        }
        work[1] = 1.f;
        t[0] = 1.f;
        idum = *lwork - m * m - 2;
        if (*n2 > 1)
        {
            slagv2_(&a[*j1 + *j1 * a_dim1], lda, &b[*j1 + *j1 * b_dim1], ldb, ar, ai, be, &work[1], &work[2], t, &t[1]);
            work[m + 1] = -work[2];
            work[m + 2] = work[1];
            t[*n2 + (*n2 << 2) - 5] = t[0];
            t[4] = -t[1];
        }
        work[m * m] = 1.f;
        t[m + (m << 2) - 5] = 1.f;
        if (*n1 > 1)
        {
            slagv2_(&a[*j1 + *n2 + (*j1 + *n2) * a_dim1], lda, &b[*j1 + *n2 + (*j1 + *n2) * b_dim1], ldb, taur, taul, &work[m * m + 1], &work[*n2 * m + *n2 + 1], &work[*n2 * m + *n2 + 2], &t[* n2 + 1 + (*n2 + 1 << 2) - 5], &t[m + (m - 1 << 2) - 5]);
            work[m * m] = work[*n2 * m + *n2 + 1];
            work[m * m - 1] = -work[*n2 * m + *n2 + 2];
            t[m + (m << 2) - 5] = t[*n2 + 1 + (*n2 + 1 << 2) - 5];
            t[m - 1 + (m << 2) - 5] = -t[m + (m - 1 << 2) - 5];
        }
        sgemm_("T", "N", n2, n1, n2, &c_b42, &work[1], &m, &a[*j1 + (*j1 + * n2) * a_dim1], lda, &c_b5, &work[m * m + 1], n2);
        slacpy_("Full", n2, n1, &work[m * m + 1], n2, &a[*j1 + (*j1 + *n2) * a_dim1], lda);
        sgemm_("T", "N", n2, n1, n2, &c_b42, &work[1], &m, &b[*j1 + (*j1 + * n2) * b_dim1], ldb, &c_b5, &work[m * m + 1], n2);
        slacpy_("Full", n2, n1, &work[m * m + 1], n2, &b[*j1 + (*j1 + *n2) * b_dim1], ldb);
        sgemm_("N", "N", &m, &m, &m, &c_b42, li, &c__4, &work[1], &m, &c_b5, & work[m * m + 1], &m);
        slacpy_("Full", &m, &m, &work[m * m + 1], &m, li, &c__4);
        sgemm_("N", "N", n2, n1, n1, &c_b42, &a[*j1 + (*j1 + *n2) * a_dim1], lda, &t[*n2 + 1 + (*n2 + 1 << 2) - 5], &c__4, &c_b5, &work[1], n2);
        slacpy_("Full", n2, n1, &work[1], n2, &a[*j1 + (*j1 + *n2) * a_dim1], lda);
        sgemm_("N", "N", n2, n1, n1, &c_b42, &b[*j1 + (*j1 + *n2) * b_dim1], ldb, &t[*n2 + 1 + (*n2 + 1 << 2) - 5], &c__4, &c_b5, &work[1], n2);
        slacpy_("Full", n2, n1, &work[1], n2, &b[*j1 + (*j1 + *n2) * b_dim1], ldb);
        sgemm_("T", "N", &m, &m, &m, &c_b42, ir, &c__4, t, &c__4, &c_b5, & work[1], &m);
        slacpy_("Full", &m, &m, &work[1], &m, ir, &c__4);
        /* Accumulate transformations into Q and Z if requested. */
        if (*wantq)
        {
            sgemm_("N", "N", n, &m, &m, &c_b42, &q[*j1 * q_dim1 + 1], ldq, li, &c__4, &c_b5, &work[1], n);
            slacpy_("Full", n, &m, &work[1], n, &q[*j1 * q_dim1 + 1], ldq);
        }
        if (*wantz)
        {
            sgemm_("N", "N", n, &m, &m, &c_b42, &z__[*j1 * z_dim1 + 1], ldz, ir, &c__4, &c_b5, &work[1], n);
            slacpy_("Full", n, &m, &work[1], n, &z__[*j1 * z_dim1 + 1], ldz);
        }
        /* Update (A(J1:J1+M-1, M+J1:N), B(J1:J1+M-1, M+J1:N)) and */
        /* (A(1:J1-1, J1:J1+M), B(1:J1-1, J1:J1+M)). */
        i__ = *j1 + m;
        if (i__ <= *n)
        {
            i__1 = *n - i__ + 1;
            sgemm_("T", "N", &m, &i__1, &m, &c_b42, li, &c__4, &a[*j1 + i__ * a_dim1], lda, &c_b5, &work[1], &m);
            i__1 = *n - i__ + 1;
            slacpy_("Full", &m, &i__1, &work[1], &m, &a[*j1 + i__ * a_dim1], lda);
            i__1 = *n - i__ + 1;
            sgemm_("T", "N", &m, &i__1, &m, &c_b42, li, &c__4, &b[*j1 + i__ * b_dim1], ldb, &c_b5, &work[1], &m);
            i__1 = *n - i__ + 1;
            slacpy_("Full", &m, &i__1, &work[1], &m, &b[*j1 + i__ * b_dim1], ldb);
        }
        i__ = *j1 - 1;
        if (i__ > 0)
        {
            sgemm_("N", "N", &i__, &m, &m, &c_b42, &a[*j1 * a_dim1 + 1], lda, ir, &c__4, &c_b5, &work[1], &i__);
            slacpy_("Full", &i__, &m, &work[1], &i__, &a[*j1 * a_dim1 + 1], lda);
            sgemm_("N", "N", &i__, &m, &m, &c_b42, &b[*j1 * b_dim1 + 1], ldb, ir, &c__4, &c_b5, &work[1], &i__);
            slacpy_("Full", &i__, &m, &work[1], &i__, &b[*j1 * b_dim1 + 1], ldb);
        }
        /* Exit with INFO = 0 if swap was successfully performed. */
        return 0;
    }
    /* Exit with INFO = 1 if swap was rejected. */
L70:
    *info = 1;
    return 0;
    /* End of STGEX2 */
}
/* stgex2_ */
