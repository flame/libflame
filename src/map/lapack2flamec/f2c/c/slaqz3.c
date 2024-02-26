/* slaqz3.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static logical c_true = TRUE_;
static integer c_n1 = -1;
static integer c__1 = 1;
static real c_b16 = 0.f;
static real c_b17 = 1.f;
/* > \brief \b SLAQZ3 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAQZ3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqz3. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqz3. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqz3. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAQZ3( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NW, A, LDA, B, */
/* $ LDB, Q, LDQ, Z, LDZ, NS, ND, ALPHAR, ALPHAI, BETA, QC, LDQC, */
/* $ ZC, LDZC, WORK, LWORK, REC, INFO ) */
/* IMPLICIT NONE */
/* Arguments */
/* LOGICAL, INTENT( IN ) :: ILSCHUR, ILQ, ILZ */
/* INTEGER, INTENT( IN ) :: N, ILO, IHI, NW, LDA, LDB, LDQ, LDZ, */
/* $ LDQC, LDZC, LWORK, REC */
/* REAL, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ Z( LDZ, * ), ALPHAR( * ), ALPHAI( * ), BETA( * ) */
/* INTEGER, INTENT( OUT ) :: NS, ND, INFO */
/* REAL :: QC( LDQC, * ), ZC( LDZC, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQZ3 performs AED */
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
/* > ILO and IHI mark the rows and columns of (A,B) which */
/* > are to be normalized */
/* > \endverbatim */
/* > */
/* > \param[in] NW */
/* > \verbatim */
/* > NW is INTEGER */
/* > The desired size of the deflation window. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDQ, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] NS */
/* > \verbatim */
/* > NS is INTEGER */
/* > The number of unconverged eigenvalues available to */
/* > use as shifts. */
/* > \endverbatim */
/* > */
/* > \param[out] ND */
/* > \verbatim */
/* > ND is INTEGER */
/* > The number of converged eigenvalues found. */
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
/* > \param[in,out] QC */
/* > \verbatim */
/* > QC is REAL array, dimension (LDQC, NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQC */
/* > \verbatim */
/* > LDQC is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] ZC */
/* > \verbatim */
/* > ZC is REAL array, dimension (LDZC, NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZC */
/* > \verbatim */
/* > LDZ is INTEGER */
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
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
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
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Thijs Steel, KU Leuven */
/* > \date May 2020 */
/* > \ingroup doubleGEcomputational */
/* > */
/* ===================================================================== */
/* Subroutine */
int slaqz3_(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nw, real *a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *z__, integer *ldz, integer *ns, integer *nd, real *alphar, real *alphai, real *beta, real *qc, integer *ldqc, real *zc, integer *ldzc, real * work, integer *lwork, integer *rec, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slaqz3 inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", nw %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS ", ldq %" FLA_IS ", ldz %" FLA_IS ", ns %" FLA_IS ", nd %" FLA_IS ", ldqc %" FLA_IS ", ldzc %" FLA_IS ", rec %" FLA_IS "",*n, *ilo, *ihi, *nw, *lda, *ldb, *ldq, *ldz, *ns, *nd, *ldqc, *ldzc, *rec);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, qc_dim1, qc_offset, zc_dim1, zc_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer lworkreq, i__, j, k;
    real s, c1;
    integer k2;
    real s1;
    integer jw;
    real ulp;
    integer stgexc_info__, ifst;
    real temp;
    integer ilst;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *), slag2_(real *, integer *, real *, integer *, real *, real *, real *, real *, real *, real *);
    logical bulge;
    real atemp;
    extern /* Subroutine */
    int sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    integer kwbot, kwtop, qz_small_info__;
    extern /* Subroutine */
    int slaqz0_(char *, char *, char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, real *, real *, real *, integer *, real *, integer *, real *, integer *, integer *, integer *), slaqz2_( logical *, logical *, integer *, integer *, integer *, integer *, real *, integer *, real *, integer *, integer *, integer *, real *, integer *, integer *, integer *, real *, integer *), slabad_( real *, real *);
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real safmax;
    extern /* Subroutine */
    int slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *), slaset_(char *, integer *, integer *, real *, real *, real *, integer *), stgexc_( logical *, logical *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *, integer *, real *, integer *, integer *), slartg_(real *, real *, real *, real *, real *);
    integer istopm;
    real smlnum;
    integer istartm;
    /* Arguments */
    /* Parameters */
    /* Local Scalars */
    /* External Functions */
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
    --alphar;
    --alphai;
    --beta;
    qc_dim1 = *ldqc;
    qc_offset = 1 + qc_dim1;
    qc -= qc_offset;
    zc_dim1 = *ldzc;
    zc_offset = 1 + zc_dim1;
    zc -= zc_offset;
    --work;
    /* Function Body */
    *info = 0;
    /* Set up deflation window */
    /* Computing MIN */
    i__1 = *nw;
    i__2 = *ihi - *ilo + 1; // , expr subst
    jw = fla_min(i__1,i__2);
    kwtop = *ihi - jw + 1;
    if (kwtop == *ilo)
    {
        s = 0.f;
    }
    else
    {
        s = a[kwtop + (kwtop - 1) * a_dim1];
    }
    /* Determine required workspace */
    ifst = 1;
    ilst = jw;
    stgexc_(&c_true, &c_true, &jw, &a[a_offset], lda, &b[b_offset], ldb, &qc[ qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, &ilst, &work[1], & c_n1, &stgexc_info__);
    lworkreq = (integer) work[1];
    i__1 = *rec + 1;
    slaqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, & b[kwtop + kwtop * b_dim1], ldb, &alphar[1], &alphai[1], &beta[1], &qc[qc_offset], ldqc, &zc[zc_offset], ldzc, &work[1], &c_n1, & i__1, &qz_small_info__);
    /* Computing MAX */
    /* Computing 2nd power */
    i__3 = jw;
    i__1 = lworkreq;
    i__2 = (integer) work[1] + (i__3 * i__3 << 1); // , expr subst
    lworkreq = fla_max(i__1,i__2);
    /* Computing MAX */
    /* Computing 2nd power */
    i__3 = *nw;
    i__1 = lworkreq, i__2 = *n * *nw;
    i__1 = fla_max(i__1,i__2);
    i__2 = (i__3 * i__3 << 1) + *n; // ; expr subst
    lworkreq = fla_max(i__1,i__2);
    if (*lwork == -1)
    {
        /* workspace query, quick return */
        work[1] = (real) lworkreq;
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (*lwork < lworkreq)
    {
        *info = -26;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLAQZ3", &i__1, (ftnlen)6);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Get machine constants */
    safmin = slamch_("SAFE MINIMUM");
    safmax = 1.f / safmin;
    slabad_(&safmin, &safmax);
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real) (*n) / ulp);
    if (*ihi == kwtop)
    {
        /* 1 by 1 deflation window, just try a regular deflation */
        alphar[kwtop] = a[kwtop + kwtop * a_dim1];
        alphai[kwtop] = 0.f;
        beta[kwtop] = b[kwtop + kwtop * b_dim1];
        *ns = 1;
        *nd = 0;
        /* Computing MAX */
        r__2 = smlnum;
        r__3 = ulp * (r__1 = a[kwtop + kwtop * a_dim1], f2c_abs( r__1)); // , expr subst
        if (f2c_abs(s) <= fla_max(r__2,r__3))
        {
            *ns = 0;
            *nd = 1;
            if (kwtop > *ilo)
            {
                a[kwtop + (kwtop - 1) * a_dim1] = 0.f;
            }
        }
    }
    /* Store window in case of convergence failure */
    slacpy_("ALL", &jw, &jw, &a[kwtop + kwtop * a_dim1], lda, &work[1], &jw);
    /* Computing 2nd power */
    i__1 = jw;
    slacpy_("ALL", &jw, &jw, &b[kwtop + kwtop * b_dim1], ldb, &work[i__1 * i__1 + 1], &jw);
    /* Transform window to real schur form */
    slaset_("FULL", &jw, &jw, &c_b16, &c_b17, &qc[qc_offset], ldqc) ;
    slaset_("FULL", &jw, &jw, &c_b16, &c_b17, &zc[zc_offset], ldzc) ;
    /* Computing 2nd power */
    i__1 = jw;
    /* Computing 2nd power */
    i__3 = jw;
    i__2 = *lwork - (i__3 * i__3 << 1);
    i__4 = *rec + 1;
    slaqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, & b[kwtop + kwtop * b_dim1], ldb, &alphar[1], &alphai[1], &beta[1], &qc[qc_offset], ldqc, &zc[zc_offset], ldzc, &work[(i__1 * i__1 << 1) + 1], &i__2, &i__4, &qz_small_info__);
    if (qz_small_info__ != 0)
    {
        /* Convergence failure, restore the window and exit */
        *nd = 0;
        *ns = jw - qz_small_info__;
        slacpy_("ALL", &jw, &jw, &work[1], &jw, &a[kwtop + kwtop * a_dim1], lda);
        /* Computing 2nd power */
        i__1 = jw;
        slacpy_("ALL", &jw, &jw, &work[i__1 * i__1 + 1], &jw, &b[kwtop + kwtop * b_dim1], ldb);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Deflation detection loop */
    if (kwtop == *ilo || s == 0.f)
    {
        kwbot = kwtop - 1;
    }
    else
    {
        kwbot = *ihi;
        k = 1;
        k2 = 1;
        while(k <= jw)
        {
            bulge = FALSE_;
            if (kwbot - kwtop + 1 >= 2)
            {
                bulge = a[kwbot + (kwbot - 1) * a_dim1] != 0.f;
            }
            if (bulge)
            {
                /* Try to deflate complex conjugate eigenvalue pair */
                temp = (r__3 = a[kwbot + kwbot * a_dim1], f2c_abs(r__3)) + sqrt(( r__1 = a[kwbot + (kwbot - 1) * a_dim1], f2c_abs(r__1))) * sqrt((r__2 = a[kwbot - 1 + kwbot * a_dim1], f2c_abs(r__2)) );
                if (temp == 0.f)
                {
                    temp = f2c_abs(s);
                }
                /* Computing MAX */
                r__3 = (r__1 = s * qc[(kwbot - kwtop) * qc_dim1 + 1], f2c_abs( r__1));
                r__4 = (r__2 = s * qc[(kwbot - kwtop + 1) * qc_dim1 + 1], f2c_abs(r__2)); // , expr subst
                /* Computing MAX */
                r__5 = smlnum;
                r__6 = ulp * temp; // , expr subst
                if (fla_max(r__3,r__4) <= fla_max(r__5,r__6))
                {
                    /* Deflatable */
                    kwbot += -2;
                }
                else
                {
                    /* Not deflatable, move out of the way */
                    ifst = kwbot - kwtop + 1;
                    ilst = k2;
                    stgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1], lda, &b[kwtop + kwtop * b_dim1], ldb, &qc[ qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, & ilst, &work[1], lwork, &stgexc_info__);
                    k2 += 2;
                }
                k += 2;
            }
            else
            {
                /* Try to deflate real eigenvalue */
                temp = (r__1 = a[kwbot + kwbot * a_dim1], f2c_abs(r__1));
                if (temp == 0.f)
                {
                    temp = f2c_abs(s);
                }
                /* Computing MAX */
                r__2 = ulp * temp;
                if ((r__1 = s * qc[(kwbot - kwtop + 1) * qc_dim1 + 1], f2c_abs( r__1)) <= fla_max(r__2,smlnum))
                {
                    /* Deflatable */
                    --kwbot;
                }
                else
                {
                    /* Not deflatable, move out of the way */
                    ifst = kwbot - kwtop + 1;
                    ilst = k2;
                    stgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1], lda, &b[kwtop + kwtop * b_dim1], ldb, &qc[ qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, & ilst, &work[1], lwork, &stgexc_info__);
                    ++k2;
                }
                ++k;
            }
        }
    }
    /* Store eigenvalues */
    *nd = *ihi - kwbot;
    *ns = jw - *nd;
    k = kwtop;
    while(k <= *ihi)
    {
        bulge = FALSE_;
        if (k < *ihi)
        {
            if (a[k + 1 + k * a_dim1] != 0.f)
            {
                bulge = TRUE_;
            }
        }
        if (bulge)
        {
            /* 2x2 eigenvalue block */
            slag2_(&a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb, &safmin, &beta[k], &beta[k + 1], &alphar[k], &alphar[k + 1], & alphai[k]);
            alphai[k + 1] = -alphai[k];
            k += 2;
        }
        else
        {
            /* 1x1 eigenvalue block */
            alphar[k] = a[k + k * a_dim1];
            alphai[k] = 0.f;
            beta[k] = b[k + k * b_dim1];
            ++k;
        }
    }
    if (kwtop != *ilo && s != 0.f)
    {
        /* Reflect spike back, this will create optimally packed bulges */
        /* A( KWTOP:KWBOT, KWTOP-1 ) = A( KWTOP, KWTOP-1 )*QC( 1, */
        /* $ 1:JW-ND ) */
        atemp = a[kwtop + (kwtop - 1) * a_dim1];
        j = 1;
        i__1 = kwbot;
        for (i__ = kwtop;
                i__ <= i__1;
                ++i__)
        {
            a[i__ + (kwtop - 1) * a_dim1] = atemp * qc[j * qc_dim1 + 1];
            ++j;
        }
        i__1 = kwtop;
        for (k = kwbot - 1;
                k >= i__1;
                --k)
        {
            slartg_(&a[k + (kwtop - 1) * a_dim1], &a[k + 1 + (kwtop - 1) * a_dim1], &c1, &s1, &temp);
            a[k + (kwtop - 1) * a_dim1] = temp;
            a[k + 1 + (kwtop - 1) * a_dim1] = 0.f;
            /* Computing MAX */
            i__2 = kwtop;
            i__3 = k - 1; // , expr subst
            k2 = fla_max(i__2,i__3);
            i__2 = *ihi - k2 + 1;
            srot_(&i__2, &a[k + k2 * a_dim1], lda, &a[k + 1 + k2 * a_dim1], lda, &c1, &s1);
            i__2 = *ihi - (k - 1) + 1;
            srot_(&i__2, &b[k + (k - 1) * b_dim1], ldb, &b[k + 1 + (k - 1) * b_dim1], ldb, &c1, &s1);
            srot_(&jw, &qc[(k - kwtop + 1) * qc_dim1 + 1], &c__1, &qc[(k + 1 - kwtop + 1) * qc_dim1 + 1], &c__1, &c1, &s1);
        }
        /* Chase bulges down */
        istartm = kwtop;
        istopm = *ihi;
        k = kwbot - 1;
        while(k >= kwtop)
        {
            if (k >= kwtop + 1 && a[k + 1 + (k - 1) * a_dim1] != 0.f)
            {
                /* Move double pole block down and remove it */
                i__1 = kwbot - 2;
                for (k2 = k - 1;
                        k2 <= i__1;
                        ++k2)
                {
                    i__2 = kwtop + jw - 1;
                    slaqz2_(&c_true, &c_true, &k2, &kwtop, &i__2, &kwbot, &a[ a_offset], lda, &b[b_offset], ldb, &jw, &kwtop, & qc[qc_offset], ldqc, &jw, &kwtop, &zc[zc_offset], ldzc);
                }
                k += -2;
            }
            else
            {
                /* k points to single shift */
                i__1 = kwbot - 2;
                for (k2 = k;
                        k2 <= i__1;
                        ++k2)
                {
                    /* Move shift down */
                    slartg_(&b[k2 + 1 + (k2 + 1) * b_dim1], &b[k2 + 1 + k2 * b_dim1], &c1, &s1, &temp);
                    b[k2 + 1 + (k2 + 1) * b_dim1] = temp;
                    b[k2 + 1 + k2 * b_dim1] = 0.f;
                    i__2 = k2 + 2 - istartm + 1;
                    srot_(&i__2, &a[istartm + (k2 + 1) * a_dim1], &c__1, &a[ istartm + k2 * a_dim1], &c__1, &c1, &s1);
                    i__2 = k2 - istartm + 1;
                    srot_(&i__2, &b[istartm + (k2 + 1) * b_dim1], &c__1, &b[ istartm + k2 * b_dim1], &c__1, &c1, &s1);
                    srot_(&jw, &zc[(k2 + 1 - kwtop + 1) * zc_dim1 + 1], &c__1, &zc[(k2 - kwtop + 1) * zc_dim1 + 1], &c__1, &c1, &s1);
                    slartg_(&a[k2 + 1 + k2 * a_dim1], &a[k2 + 2 + k2 * a_dim1], &c1, &s1, &temp);
                    a[k2 + 1 + k2 * a_dim1] = temp;
                    a[k2 + 2 + k2 * a_dim1] = 0.f;
                    i__2 = istopm - k2;
                    srot_(&i__2, &a[k2 + 1 + (k2 + 1) * a_dim1], lda, &a[k2 + 2 + (k2 + 1) * a_dim1], lda, &c1, &s1);
                    i__2 = istopm - k2;
                    srot_(&i__2, &b[k2 + 1 + (k2 + 1) * b_dim1], ldb, &b[k2 + 2 + (k2 + 1) * b_dim1], ldb, &c1, &s1);
                    srot_(&jw, &qc[(k2 + 1 - kwtop + 1) * qc_dim1 + 1], &c__1, &qc[(k2 + 2 - kwtop + 1) * qc_dim1 + 1], &c__1, & c1, &s1);
                }
                /* Remove the shift */
                slartg_(&b[kwbot + kwbot * b_dim1], &b[kwbot + (kwbot - 1) * b_dim1], &c1, &s1, &temp);
                b[kwbot + kwbot * b_dim1] = temp;
                b[kwbot + (kwbot - 1) * b_dim1] = 0.f;
                i__1 = kwbot - istartm;
                srot_(&i__1, &b[istartm + kwbot * b_dim1], &c__1, &b[istartm + (kwbot - 1) * b_dim1], &c__1, &c1, &s1);
                i__1 = kwbot - istartm + 1;
                srot_(&i__1, &a[istartm + kwbot * a_dim1], &c__1, &a[istartm + (kwbot - 1) * a_dim1], &c__1, &c1, &s1);
                srot_(&jw, &zc[(kwbot - kwtop + 1) * zc_dim1 + 1], &c__1, &zc[ (kwbot - 1 - kwtop + 1) * zc_dim1 + 1], &c__1, &c1, & s1);
                --k;
            }
        }
    }
    /* Apply Qc and Zc to rest of the matrix */
    if (*ilschur)
    {
        istartm = 1;
        istopm = *n;
    }
    else
    {
        istartm = *ilo;
        istopm = *ihi;
    }
    if (istopm - *ihi > 0)
    {
        i__1 = istopm - *ihi;
        sgemm_("T", "N", &jw, &i__1, &jw, &c_b17, &qc[qc_offset], ldqc, &a[ kwtop + (*ihi + 1) * a_dim1], lda, &c_b16, &work[1], &jw);
        i__1 = istopm - *ihi;
        slacpy_("ALL", &jw, &i__1, &work[1], &jw, &a[kwtop + (*ihi + 1) * a_dim1], lda);
        i__1 = istopm - *ihi;
        sgemm_("T", "N", &jw, &i__1, &jw, &c_b17, &qc[qc_offset], ldqc, &b[ kwtop + (*ihi + 1) * b_dim1], ldb, &c_b16, &work[1], &jw);
        i__1 = istopm - *ihi;
        slacpy_("ALL", &jw, &i__1, &work[1], &jw, &b[kwtop + (*ihi + 1) * b_dim1], ldb);
    }
    if (*ilq)
    {
        sgemm_("N", "N", n, &jw, &jw, &c_b17, &q[kwtop * q_dim1 + 1], ldq, & qc[qc_offset], ldqc, &c_b16, &work[1], n);
        slacpy_("ALL", n, &jw, &work[1], n, &q[kwtop * q_dim1 + 1], ldq);
    }
    if (kwtop - 1 - istartm + 1 > 0)
    {
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        sgemm_("N", "N", &i__1, &jw, &jw, &c_b17, &a[istartm + kwtop * a_dim1], lda, &zc[zc_offset], ldzc, &c_b16, &work[1], &i__2);
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        slacpy_("ALL", &i__1, &jw, &work[1], &i__2, &a[istartm + kwtop * a_dim1], lda);
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        sgemm_("N", "N", &i__1, &jw, &jw, &c_b17, &b[istartm + kwtop * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b16, &work[1], &i__2);
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        slacpy_("ALL", &i__1, &jw, &work[1], &i__2, &b[istartm + kwtop * b_dim1], ldb);
    }
    if (*ilz)
    {
        sgemm_("N", "N", n, &jw, &jw, &c_b17, &z__[kwtop * z_dim1 + 1], ldz, & zc[zc_offset], ldzc, &c_b16, &work[1], n);
        slacpy_("ALL", n, &jw, &work[1], n, &z__[kwtop * z_dim1 + 1], ldz);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
}
/* slaqz3_ */
