/* dlaqz3.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static logical c_true = TRUE_;
static integer c_n1 = -1;
static integer c__1 = 1;
static doublereal c_b16 = 0.;
static doublereal c_b17 = 1.;
/* > \brief \b DLAQZ3 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQZ3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqz3. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqz3. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqz3. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQZ3( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NW, A, LDA, B, */
/* $ LDB, Q, LDQ, Z, LDZ, NS, ND, ALPHAR, ALPHAI, BETA, QC, LDQC, */
/* $ ZC, LDZC, WORK, LWORK, REC, INFO ) */
/* IMPLICIT NONE */
/* Arguments */
/* LOGICAL, INTENT( IN ) :: ILSCHUR, ILQ, ILZ */
/* INTEGER, INTENT( IN ) :: N, ILO, IHI, NW, LDA, LDB, LDQ, LDZ, */
/* $ LDQC, LDZC, LWORK, REC */
/* DOUBLE PRECISION, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), */
/* $ Q( LDQ, * ), Z( LDZ, * ), ALPHAR( * ), ALPHAI( * ), BETA( * ) */
/* INTEGER, INTENT( OUT ) :: NS, ND, INFO */
/* DOUBLE PRECISION :: QC( LDQC, * ), ZC( LDZC, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQZ3 performs AED */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ILSCHUR */
/* > \verbatim */
/* > ILSCHUR is LOGICAL */
/* > Determines whether or not to update the full Schur form */
/* > \endverbatim */
/* > */
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
/* > A is DOUBLE PRECISION array, dimension (LDA, N) */
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
/* > B is DOUBLE PRECISION array, dimension (LDB, N) */
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
/* > Q is DOUBLE PRECISION array, dimension (LDQ, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension (LDZ, N) */
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
/* > ALPHAR is DOUBLE PRECISION array, dimension (N) */
/* > The real parts of each scalar alpha defining an eigenvalue */
/* > of GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* > ALPHAI is DOUBLE PRECISION array, dimension (N) */
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
/* > BETA is DOUBLE PRECISION array, dimension (N) */
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
/* > QC is DOUBLE PRECISION array, dimension (LDQC, NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQC */
/* > \verbatim */
/* > LDQC is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] ZC */
/* > \verbatim */
/* > ZC is DOUBLE PRECISION array, dimension (LDZC, NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZC */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Thijs Steel, KU Leuven */
/* > \date May 2020 */
/* > \ingroup doubleGEcomputational */
/* > */
/* ===================================================================== */
/* Subroutine */
int dlaqz3_(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nw, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *q, integer * ldq, doublereal *z__, integer *ldz, integer *ns, integer *nd, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal * qc, integer *ldqc, doublereal *zc, integer *ldzc, doublereal *work, integer *lwork, integer *rec, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"dlaqz3 inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", nw %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS ", ldq %" FLA_IS ", ldz %" FLA_IS ", ldqc %" FLA_IS ", ldzc %" FLA_IS ", lwork %" FLA_IS ", rec %" FLA_IS "",*n, *ilo, *ihi, *nw, *lda, *ldb, *ldq, *ldz, *ldqc, *ldzc, *lwork, *rec);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, qc_dim1, qc_offset, zc_dim1, zc_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    complex q__1, q__2;
    doublecomplex z__1;
    /* Builtin functions */
    double sqrt(doublereal);
    void r_cnjg(complex *, complex *);
    /* Local variables */
    integer lworkreq, k;
    doublereal s, c1;
    integer k2;
    doublereal s1;
    integer jw, jki, jli;
    doublereal ulp;
    integer dtgexc_info__, ifst;
    doublereal temp;
    extern /* Subroutine */
    int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *);
    integer ilst;
    extern /* Subroutine */
    int dlag2_(doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *), dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    logical bulge;
    integer kwbot, kwtop, qz_small_info__;
    extern /* Subroutine */
    int dlaqz0_(char *, char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, integer *), dlaqz2_(logical *, logical *, integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal safmax;
    extern /* Subroutine */
    int dtgexc_(logical *, logical *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, integer *, integer *), dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *), dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    integer istopm;
    doublereal smlnum;
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
    jw = min(i__1,i__2);
    kwtop = *ihi - jw + 1;
    if (kwtop == *ilo)
    {
        s = 0.;
    }
    else
    {
        s = a[kwtop + (kwtop - 1) * a_dim1];
    }
    /* Determine required workspace */
    ifst = 1;
    ilst = jw;
    dtgexc_(&c_true, &c_true, &jw, &a[a_offset], lda, &b[b_offset], ldb, &qc[ qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, &ilst, &work[1], & c_n1, &dtgexc_info__);
    lworkreq = (integer) work[1];
    i__1 = *rec + 1;
    dlaqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, & b[kwtop + kwtop * b_dim1], ldb, &alphar[1], &alphai[1], &beta[1], &qc[qc_offset], ldqc, &zc[zc_offset], ldzc, &work[1], &c_n1, & i__1, &qz_small_info__);
    /* Computing MAX */
    /* Computing 2nd power */
    i__3 = jw;
    i__1 = lworkreq;
    i__2 = (integer) work[1] + (i__3 * i__3 << 1); // , expr subst
    lworkreq = max(i__1,i__2);
    /* Computing MAX */
    /* Computing 2nd power */
    i__3 = *nw;
    i__1 = lworkreq, i__2 = *n * *nw;
    i__1 = max(i__1,i__2);
    i__2 = (i__3 * i__3 << 1) + *n; // ; expr subst
    lworkreq = max(i__1,i__2);
    if (*lwork == -1)
    {
        /* workspace query, quick return */
        work[1] = (doublereal) lworkreq;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (*lwork < lworkreq)
    {
        *info = -26;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DLAQZ3", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Get machine constants */
    safmin = dlamch_("SAFE MINIMUM");
    safmax = 1. / safmin;
    dlabad_(&safmin, &safmax);
    ulp = dlamch_("PRECISION");
    smlnum = safmin * ((doublereal) (*n) / ulp);
    if (*ihi == kwtop)
    {
        /* 1 by 1 deflation window, just try a regular deflation */
        alphar[kwtop] = a[kwtop + kwtop * a_dim1];
        alphai[kwtop] = 0.;
        beta[kwtop] = b[kwtop + kwtop * b_dim1];
        *ns = 1;
        *nd = 0;
        /* Computing MAX */
        d__2 = smlnum;
        d__3 = ulp * (d__1 = a[kwtop + kwtop * a_dim1], f2c_abs( d__1)); // , expr subst
        if (f2c_abs(s) <= max(d__2,d__3))
        {
            *ns = 0;
            *nd = 1;
            if (kwtop > *ilo)
            {
                a[kwtop + (kwtop - 1) * a_dim1] = 0.;
            }
        }
    }
    /* Store window in case of convergence failure */
    dlacpy_("ALL", &jw, &jw, &a[kwtop + kwtop * a_dim1], lda, &work[1], &jw);
    /* Computing 2nd power */
    i__1 = jw;
    dlacpy_("ALL", &jw, &jw, &b[kwtop + kwtop * b_dim1], ldb, &work[i__1 * i__1 + 1], &jw);
    /* Transform window to real schur form */
    dlaset_("FULL", &jw, &jw, &c_b16, &c_b17, &qc[qc_offset], ldqc) ;
    dlaset_("FULL", &jw, &jw, &c_b16, &c_b17, &zc[zc_offset], ldzc) ;
    /* Computing 2nd power */
    i__1 = jw;
    /* Computing 2nd power */
    i__3 = jw;
    i__2 = *lwork - (i__3 * i__3 << 1);
    i__4 = *rec + 1;
    dlaqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, & b[kwtop + kwtop * b_dim1], ldb, &alphar[1], &alphai[1], &beta[1], &qc[qc_offset], ldqc, &zc[zc_offset], ldzc, &work[(i__1 * i__1 << 1) + 1], &i__2, &i__4, &qz_small_info__);
    if (qz_small_info__ != 0)
    {
        /* Convergence failure, restore the window and exit */
        *nd = 0;
        *ns = jw - qz_small_info__;
        dlacpy_("ALL", &jw, &jw, &work[1], &jw, &a[kwtop + kwtop * a_dim1], lda);
        /* Computing 2nd power */
        i__1 = jw;
        dlacpy_("ALL", &jw, &jw, &work[i__1 * i__1 + 1], &jw, &b[kwtop + kwtop * b_dim1], ldb);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Deflation detection loop */
    if (kwtop == *ilo || s == 0.)
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
                bulge = a[kwbot + (kwbot - 1) * a_dim1] != 0.;
            }
            if (bulge)
            {
                /* Try to deflate complex conjugate eigenvalue pair */
                temp = (d__3 = a[kwbot + kwbot * a_dim1], f2c_abs(d__3)) + sqrt(( d__1 = a[kwbot + (kwbot - 1) * a_dim1], f2c_abs(d__1))) * sqrt((d__2 = a[kwbot - 1 + kwbot * a_dim1], f2c_abs(d__2)) );
                if (temp == 0.)
                {
                    temp = f2c_abs(s);
                }
                /* Computing MAX */
                d__3 = (d__1 = s * qc[(kwbot - kwtop) * qc_dim1 + 1], f2c_abs( d__1));
                d__4 = (d__2 = s * qc[(kwbot - kwtop + 1) * qc_dim1 + 1], f2c_abs(d__2)); // , expr subst
                /* Computing MAX */
                d__5 = smlnum;
                d__6 = ulp * temp; // , expr subst
                if (max(d__3,d__4) <= max(d__5,d__6))
                {
                    /* Deflatable */
                    kwbot += -2;
                }
                else
                {
                    /* Not deflatable, move out of the way */
                    ifst = kwbot - kwtop + 1;
                    ilst = k2;
                    dtgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1], lda, &b[kwtop + kwtop * b_dim1], ldb, &qc[ qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, & ilst, &work[1], lwork, &dtgexc_info__);
                    k2 += 2;
                }
                k += 2;
            }
            else
            {
                /* Try to deflate real eigenvalue */
                temp = (d__1 = a[kwbot + kwbot * a_dim1], f2c_abs(d__1));
                if (temp == 0.)
                {
                    temp = f2c_abs(s);
                }
                /* Computing MAX */
                d__2 = ulp * temp;
                if ((d__1 = s * qc[(kwbot - kwtop + 1) * qc_dim1 + 1], f2c_abs( d__1)) <= max(d__2,smlnum))
                {
                    /* Deflatable */
                    --kwbot;
                }
                else
                {
                    /* Not deflatable, move out of the way */
                    ifst = kwbot - kwtop + 1;
                    ilst = k2;
                    dtgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1], lda, &b[kwtop + kwtop * b_dim1], ldb, &qc[ qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, & ilst, &work[1], lwork, &dtgexc_info__);
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
            if (a[k + 1 + k * a_dim1] != 0.)
            {
                bulge = TRUE_;
            }
        }
        if (bulge)
        {
            /* 2x2 eigenvalue block */
            dlag2_(&a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb, &safmin, &beta[k], &beta[k + 1], &alphar[k], &alphar[k + 1], & alphai[k]);
            alphai[k + 1] = -alphai[k];
            k += 2;
        }
        else
        {
            /* 1x1 eigenvalue block */
            alphar[k] = a[k + k * a_dim1];
            alphai[k] = 0.;
            beta[k] = b[k + k * b_dim1];
            ++k;
        }
    }
    if (kwtop != *ilo && s != 0.)
    {
        /* Reflect spike back, this will create optimally packed bulges */
        /* A( KWTOP:KWBOT, KWTOP-1 ) = A( KWTOP, KWTOP-1 )*QC( 1, */
        /* $ 1:JW-ND ) */
        i__1 = kwbot;
        for (jki = kwtop;
                jki <= i__1;
                ++jki)
        {
            i__2 = jw - *nd;
            for (jli = 1;
                    jli <= i__2;
                    ++jli)
            {
                i__3 = jli * qc_dim1 + 1;
                q__1.r = qc[i__3];
                q__1.i = 0.f; // , expr subst
                i__4 = jki + (kwtop - 1) * a_dim1;
                i__5 = kwtop + (kwtop - 1) * a_dim1;
                r_cnjg(&q__2, &q__1);
                z__1.r = a[i__5] * q__2.r;
                z__1.i = a[i__5] * q__2.i; // , expr subst
                a[i__4] = z__1.r;
            }
        }
        i__1 = kwtop;
        for (k = kwbot - 1;
                k >= i__1;
                --k)
        {
            dlartg_(&a[k + (kwtop - 1) * a_dim1], &a[k + 1 + (kwtop - 1) * a_dim1], &c1, &s1, &temp);
            a[k + (kwtop - 1) * a_dim1] = temp;
            a[k + 1 + (kwtop - 1) * a_dim1] = 0.;
            /* Computing MAX */
            i__2 = kwtop;
            i__3 = k - 1; // , expr subst
            k2 = max(i__2,i__3);
            i__2 = *ihi - k2 + 1;
            drot_(&i__2, &a[k + k2 * a_dim1], lda, &a[k + 1 + k2 * a_dim1], lda, &c1, &s1);
            i__2 = *ihi - (k - 1) + 1;
            drot_(&i__2, &b[k + (k - 1) * b_dim1], ldb, &b[k + 1 + (k - 1) * b_dim1], ldb, &c1, &s1);
            drot_(&jw, &qc[(k - kwtop + 1) * qc_dim1 + 1], &c__1, &qc[(k + 1 - kwtop + 1) * qc_dim1 + 1], &c__1, &c1, &s1);
        }
        /* Chase bulges down */
        istartm = kwtop;
        istopm = *ihi;
        k = kwbot - 1;
        while(k >= kwtop)
        {
            if (k >= kwtop + 1 && a[k + 1 + (k - 1) * a_dim1] != 0.)
            {
                /* Move double pole block down and remove it */
                i__1 = kwbot - 2;
                for (k2 = k - 1;
                        k2 <= i__1;
                        ++k2)
                {
                    i__2 = kwtop + jw - 1;
                    dlaqz2_(&c_true, &c_true, &k2, &kwtop, &i__2, &kwbot, &a[ a_offset], lda, &b[b_offset], ldb, &jw, &kwtop, & qc[qc_offset], ldqc, &jw, &kwtop, &zc[zc_offset], ldzc);
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
                    dlartg_(&b[k2 + 1 + (k2 + 1) * b_dim1], &b[k2 + 1 + k2 * b_dim1], &c1, &s1, &temp);
                    b[k2 + 1 + (k2 + 1) * b_dim1] = temp;
                    b[k2 + 1 + k2 * b_dim1] = 0.;
                    i__2 = k2 + 2 - istartm + 1;
                    drot_(&i__2, &a[istartm + (k2 + 1) * a_dim1], &c__1, &a[ istartm + k2 * a_dim1], &c__1, &c1, &s1);
                    i__2 = k2 - istartm + 1;
                    drot_(&i__2, &b[istartm + (k2 + 1) * b_dim1], &c__1, &b[ istartm + k2 * b_dim1], &c__1, &c1, &s1);
                    drot_(&jw, &zc[(k2 + 1 - kwtop + 1) * zc_dim1 + 1], &c__1, &zc[(k2 - kwtop + 1) * zc_dim1 + 1], &c__1, &c1, &s1);
                    dlartg_(&a[k2 + 1 + k2 * a_dim1], &a[k2 + 2 + k2 * a_dim1], &c1, &s1, &temp);
                    a[k2 + 1 + k2 * a_dim1] = temp;
                    a[k2 + 2 + k2 * a_dim1] = 0.;
                    i__2 = istopm - k2;
                    drot_(&i__2, &a[k2 + 1 + (k2 + 1) * a_dim1], lda, &a[k2 + 2 + (k2 + 1) * a_dim1], lda, &c1, &s1);
                    i__2 = istopm - k2;
                    drot_(&i__2, &b[k2 + 1 + (k2 + 1) * b_dim1], ldb, &b[k2 + 2 + (k2 + 1) * b_dim1], ldb, &c1, &s1);
                    drot_(&jw, &qc[(k2 + 1 - kwtop + 1) * qc_dim1 + 1], &c__1, &qc[(k2 + 2 - kwtop + 1) * qc_dim1 + 1], &c__1, & c1, &s1);
                }
                /* Remove the shift */
                dlartg_(&b[kwbot + kwbot * b_dim1], &b[kwbot + (kwbot - 1) * b_dim1], &c1, &s1, &temp);
                b[kwbot + kwbot * b_dim1] = temp;
                b[kwbot + (kwbot - 1) * b_dim1] = 0.;
                i__1 = kwbot - istartm;
                drot_(&i__1, &b[istartm + kwbot * b_dim1], &c__1, &b[istartm + (kwbot - 1) * b_dim1], &c__1, &c1, &s1);
                i__1 = kwbot - istartm + 1;
                drot_(&i__1, &a[istartm + kwbot * a_dim1], &c__1, &a[istartm + (kwbot - 1) * a_dim1], &c__1, &c1, &s1);
                drot_(&jw, &zc[(kwbot - kwtop + 1) * zc_dim1 + 1], &c__1, &zc[ (kwbot - 1 - kwtop + 1) * zc_dim1 + 1], &c__1, &c1, & s1);
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
        dgemm_("T", "N", &jw, &i__1, &jw, &c_b17, &qc[qc_offset], ldqc, &a[ kwtop + (*ihi + 1) * a_dim1], lda, &c_b16, &work[1], &jw);
        i__1 = istopm - *ihi;
        dlacpy_("ALL", &jw, &i__1, &work[1], &jw, &a[kwtop + (*ihi + 1) * a_dim1], lda);
        i__1 = istopm - *ihi;
        dgemm_("T", "N", &jw, &i__1, &jw, &c_b17, &qc[qc_offset], ldqc, &b[ kwtop + (*ihi + 1) * b_dim1], ldb, &c_b16, &work[1], &jw);
        i__1 = istopm - *ihi;
        dlacpy_("ALL", &jw, &i__1, &work[1], &jw, &b[kwtop + (*ihi + 1) * b_dim1], ldb);
    }
    if (*ilq)
    {
        dgemm_("N", "N", n, &jw, &jw, &c_b17, &q[kwtop * q_dim1 + 1], ldq, & qc[qc_offset], ldqc, &c_b16, &work[1], n);
        dlacpy_("ALL", n, &jw, &work[1], n, &q[kwtop * q_dim1 + 1], ldq);
    }
    if (kwtop - 1 - istartm + 1 > 0)
    {
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        dgemm_("N", "N", &i__1, &jw, &jw, &c_b17, &a[istartm + kwtop * a_dim1], lda, &zc[zc_offset], ldzc, &c_b16, &work[1], &i__2);
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        dlacpy_("ALL", &i__1, &jw, &work[1], &i__2, &a[istartm + kwtop * a_dim1], lda);
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        dgemm_("N", "N", &i__1, &jw, &jw, &c_b17, &b[istartm + kwtop * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b16, &work[1], &i__2);
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        dlacpy_("ALL", &i__1, &jw, &work[1], &i__2, &b[istartm + kwtop * b_dim1], ldb);
    }
    if (*ilz)
    {
        dgemm_("N", "N", n, &jw, &jw, &c_b17, &z__[kwtop * z_dim1 + 1], ldz, & zc[zc_offset], ldzc, &c_b16, &work[1], n);
        dlacpy_("ALL", n, &jw, &work[1], n, &z__[kwtop * z_dim1 + 1], ldz);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
}
/* dlaqz3_ */

