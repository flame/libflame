/* ../netlib/cuncsd2by1.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__0 = 0;
static integer c_n1 = -1;
static integer c__1 = 1;
static logical c_false = FALSE_;
/* > \brief \b CUNCSD2BY1 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CUNCSD2BY1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cuncsd2 by1.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cuncsd2 by1.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cuncsd2 by1.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CUNCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11, */
/* X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T, */
/* LDV1T, WORK, LWORK, RWORK, LRWORK, IWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBU1, JOBU2, JOBV1T */
/* INTEGER INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21, */
/* $ M, P, Q */
/* INTEGER LRWORK, LRWORKMIN, LRWORKOPT */
/* .. */
/* .. Array Arguments .. */
/* REAL RWORK(*) */
/* REAL THETA(*) */
/* COMPLEX U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*), */
/* $ X11(LDX11,*), X21(LDX21,*) */
/* INTEGER IWORK(*) */
/* .. */
/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > */
/* > CUNCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with */
/* > orthonormal columns that has been partitioned into a 2-by-1 block */
/* > structure: */
/* > */
/* > [ I 0 0 ] */
/* > [ 0 C 0 ] */
/* > [ X11 ] [ U1 | ] [ 0 0 0 ] */
/* > X = [-----] = [---------] [----------] V1**T . */
/* > [ X21 ] [ | U2 ] [ 0 0 0 ] */
/* > [ 0 S 0 ] */
/* > [ 0 0 I ] */
/* > */
/* > X11 is P-by-Q. The unitary matrices U1, U2, V1, and V2 are P-by-P, */
/* > (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are */
/* > R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in */
/* > which R = MIN(P,M-P,Q,M-Q). */
/* > */
/* >\endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBU1 */
/* > \verbatim */
/* > JOBU1 is CHARACTER */
/* > = 'Y': U1 is computed;
*/
/* > otherwise: U1 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU2 */
/* > \verbatim */
/* > JOBU2 is CHARACTER */
/* > = 'Y': U2 is computed;
*/
/* > otherwise: U2 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV1T */
/* > \verbatim */
/* > JOBV1T is CHARACTER */
/* > = 'Y': V1T is computed;
*/
/* > otherwise: V1T is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows and columns in X. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* > P is INTEGER */
/* > The number of rows in X11 and X12. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* > Q is INTEGER */
/* > The number of columns in X11 and X21. 0 <= Q <= M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* > X11 is COMPLEX array, dimension (LDX11,Q) */
/* > On entry, part of the unitary matrix whose CSD is */
/* > desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* > LDX11 is INTEGER */
/* > The leading dimension of X11. LDX11 >= MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* > X21 is COMPLEX array, dimension (LDX21,Q) */
/* > On entry, part of the unitary matrix whose CSD is */
/* > desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* > LDX21 is INTEGER */
/* > The leading dimension of X21. LDX21 >= MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* > THETA is COMPLEX array, dimension (R), in which R = */
/* > MIN(P,M-P,Q,M-Q). */
/* > C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and */
/* > S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ). */
/* > \endverbatim */
/* > */
/* > \param[out] U1 */
/* > \verbatim */
/* > U1 is COMPLEX array, dimension (P) */
/* > If JOBU1 = 'Y', U1 contains the P-by-P unitary matrix U1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU1 */
/* > \verbatim */
/* > LDU1 is INTEGER */
/* > The leading dimension of U1. If JOBU1 = 'Y', LDU1 >= */
/* > MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[out] U2 */
/* > \verbatim */
/* > U2 is COMPLEX array, dimension (M-P) */
/* > If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) unitary */
/* > matrix U2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU2 */
/* > \verbatim */
/* > LDU2 is INTEGER */
/* > The leading dimension of U2. If JOBU2 = 'Y', LDU2 >= */
/* > MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] V1T */
/* > \verbatim */
/* > V1T is COMPLEX array, dimension (Q) */
/* > If JOBV1T = 'Y', V1T contains the Q-by-Q matrix unitary */
/* > matrix V1**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV1T */
/* > \verbatim */
/* > LDV1T is INTEGER */
/* > The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >= */
/* > MAX(1,Q). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > If INFO > 0 on exit, WORK(2:R) contains the values PHI(1), */
/* > ..., PHI(R-1) that, together with THETA(1), ..., THETA(R), */
/* > define the matrix in intermediate bidiagonal-block form */
/* > remaining after nonconvergence. INFO specifies the number */
/* > of nonzero PHI's. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > \endverbatim */
/* > \verbatim */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the work array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* > On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > If INFO > 0 on exit, RWORK(2:R) contains the values PHI(1), */
/* > ..., PHI(R-1) that, together with THETA(1), ..., THETA(R), */
/* > define the matrix in intermediate bidiagonal-block form */
/* > remaining after nonconvergence. INFO specifies the number */
/* > of nonzero PHI's. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* > LRWORK is INTEGER */
/* > The dimension of the array RWORK. */
/* > */
/* > If LRWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the RWORK array, returns */
/* > this value as the first entry of the work array, and no error */
/* > message related to LRWORK is issued by XERBLA. */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (M-MIN(P,M-P,Q,M-Q)) */
/* > \endverbatim */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: CBBCSD did not converge. See the description of WORK */
/* > above for details. */
/* > \endverbatim */
/* > \par References: */
/* > ================ */
/* > */
/* > [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* > Algorithms, 50(1):33-65, 2009. */
/* > */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date July 2012 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int cuncsd2by1_(char *jobu1, char *jobu2, char *jobv1t, integer *m, integer *p, integer *q, complex *x11, integer *ldx11, complex *x21, integer *ldx21, real *theta, complex *u1, integer *ldu1, complex *u2, integer *ldu2, complex *v1t, integer *ldv1t, complex * work, integer *lwork, real *rwork, integer *lrwork, integer *iwork, integer *info)
{
    /* System generated locals */
    integer u1_dim1, u1_offset, u2_dim1, u2_offset, v1t_dim1, v1t_offset, x11_dim1, x11_offset, x21_dim1, x21_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9;
    /* Local variables */
    integer lworkmin, lworkopt, i__, j, r__, childinfo, lorglqmin, lorgqrmin, lorglqopt, lrworkmin, lorgqropt, lrworkopt, ib11d, ib11e, ib12d, ib12e, ib21d, ib21e, ib22d, ib22e, iphi;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int ccopy_(integer *, complex *, integer *, complex *, integer *);
    integer itaup1, itaup2, itauq1;
    logical wantu1, wantu2;
    extern /* Subroutine */
    int cbbcsd_();
    integer ibbcsd, lbbcsd, iorbdb, lorbdb;
    extern /* Subroutine */
    int clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), xerbla_(char *, integer *), clapmr_(logical *, integer *, integer *, complex *, integer *, integer *), clapmt_(logical *, integer *, integer *, complex *, integer *, integer *), cunglq_();
    integer iorglq;
    extern /* Subroutine */
    int cungqr_();
    integer lorglq, iorgqr, lorgqr;
    extern /* Subroutine */
    int cunbdb1_(), cunbdb2_();
    logical lquery;
    extern /* Subroutine */
    int cunbdb3_(), cunbdb4_();
    logical wantv1t;
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* July 2012 */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Function .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test input arguments */
    /* Parameter adjustments */
    x11_dim1 = *ldx11;
    x11_offset = 1 + x11_dim1;
    x11 -= x11_offset;
    x21_dim1 = *ldx21;
    x21_offset = 1 + x21_dim1;
    x21 -= x21_offset;
    --theta;
    u1_dim1 = *ldu1;
    u1_offset = 1 + u1_dim1;
    u1 -= u1_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    v1t_dim1 = *ldv1t;
    v1t_offset = 1 + v1t_dim1;
    v1t -= v1t_offset;
    --work;
    --rwork;
    --iwork;
    /* Function Body */
    *info = 0;
    wantu1 = lsame_(jobu1, "Y");
    wantu2 = lsame_(jobu2, "Y");
    wantv1t = lsame_(jobv1t, "Y");
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -4;
    }
    else if (*p < 0 || *p > *m)
    {
        *info = -5;
    }
    else if (*q < 0 || *q > *m)
    {
        *info = -6;
    }
    else if (*ldx11 < max(1,*p))
    {
        *info = -8;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *p; // , expr subst
        if (*ldx21 < max(i__1,i__2))
        {
            *info = -10;
        }
        else if (wantu1 && *ldu1 < *p)
        {
            *info = -13;
        }
        else if (wantu2 && *ldu2 < *m - *p)
        {
            *info = -15;
        }
        else if (wantv1t && *ldv1t < *q)
        {
            *info = -17;
        }
    }
    /* Computing MIN */
    i__1 = *p, i__2 = *m - *p, i__1 = min(i__1,i__2);
    i__1 = min(i__1,*q);
    i__2 = *m - *q; // ; expr subst
    r__ = min(i__1,i__2);
    /* Compute workspace */
    /* WORK layout: */
    /* |-----------------------------------------| */
    /* | LWORKOPT (1) | */
    /* |-----------------------------------------| */
    /* | TAUP1 (MAX(1,P)) | */
    /* | TAUP2 (MAX(1,M-P)) | */
    /* | TAUQ1 (MAX(1,Q)) | */
    /* |-----------------------------------------| */
    /* | CUNBDB WORK | CUNGQR WORK | CUNGLQ WORK | */
    /* | | | | */
    /* | | | | */
    /* | | | | */
    /* | | | | */
    /* |-----------------------------------------| */
    /* RWORK layout: */
    /* |------------------| */
    /* | LRWORKOPT (1) | */
    /* |------------------| */
    /* | PHI (MAX(1,R-1)) | */
    /* |------------------| */
    /* | B11D (R) | */
    /* | B11E (R-1) | */
    /* | B12D (R) | */
    /* | B12E (R-1) | */
    /* | B21D (R) | */
    /* | B21E (R-1) | */
    /* | B22D (R) | */
    /* | B22E (R-1) | */
    /* | CBBCSD RWORK | */
    /* |------------------| */
    if (*info == 0)
    {
        iphi = 2;
        /* Computing MAX */
        i__1 = 1;
        i__2 = r__ - 1; // , expr subst
        ib11d = iphi + max(i__1,i__2);
        ib11e = ib11d + max(1,r__);
        /* Computing MAX */
        i__1 = 1;
        i__2 = r__ - 1; // , expr subst
        ib12d = ib11e + max(i__1,i__2);
        ib12e = ib12d + max(1,r__);
        /* Computing MAX */
        i__1 = 1;
        i__2 = r__ - 1; // , expr subst
        ib21d = ib12e + max(i__1,i__2);
        ib21e = ib21d + max(1,r__);
        /* Computing MAX */
        i__1 = 1;
        i__2 = r__ - 1; // , expr subst
        ib22d = ib21e + max(i__1,i__2);
        ib22e = ib22d + max(1,r__);
        /* Computing MAX */
        i__1 = 1;
        i__2 = r__ - 1; // , expr subst
        ibbcsd = ib22e + max(i__1,i__2);
        itaup1 = 2;
        itaup2 = itaup1 + max(1,*p);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *p; // , expr subst
        itauq1 = itaup2 + max(i__1,i__2);
        iorbdb = itauq1 + max(1,*q);
        iorgqr = itauq1 + max(1,*q);
        iorglq = itauq1 + max(1,*q);
        if (r__ == *q)
        {
            cunbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], & c_n1, &childinfo);
            lorbdb = (integer) work[1].r;
            if (*p >= *m - *p)
            {
                cungqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1, &childinfo);
                lorgqrmin = max(1,*p);
                lorgqropt = (integer) work[1].r;
            }
            else
            {
                i__1 = *m - *p;
                i__2 = *m - *p;
                cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1] , &c_n1, &childinfo);
                /* Computing MAX */
                i__1 = 1;
                i__2 = *m - *p; // , expr subst
                lorgqrmin = max(i__1,i__2);
                lorgqropt = (integer) work[1].r;
            }
            /* Computing MAX */
            i__2 = 0;
            i__3 = *q - 1; // , expr subst
            i__1 = max(i__2,i__3);
            /* Computing MAX */
            i__5 = 0;
            i__6 = *q - 1; // , expr subst
            i__4 = max(i__5,i__6);
            /* Computing MAX */
            i__8 = 0;
            i__9 = *q - 1; // , expr subst
            i__7 = max(i__8,i__9);
            cunglq_(&i__1, &i__4, &i__7, &v1t[v1t_offset], ldv1t, &c__0, & work[1], &c_n1, &childinfo);
            /* Computing MAX */
            i__1 = 1;
            i__2 = *q - 1; // , expr subst
            lorglqmin = max(i__1,i__2);
            lorglqopt = (integer) work[1].r;
            cbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &c__0, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[ v1t_offset], ldv1t, &c__0, &c__1, &c__0, &c__0, &c__0, & c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, & childinfo);
            lbbcsd = (integer) rwork[1];
        }
        else if (r__ == *p)
        {
            cunbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], & c_n1, &childinfo);
            lorbdb = (integer) work[1].r;
            if (*p - 1 >= *m - *p)
            {
                i__1 = *p - 1;
                i__2 = *p - 1;
                i__3 = *p - 1;
                cungqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, & c__0, &work[1], &c_n1, &childinfo);
                /* Computing MAX */
                i__1 = 1;
                i__2 = *p - 1; // , expr subst
                lorgqrmin = max(i__1,i__2);
                lorgqropt = (integer) work[1].r;
            }
            else
            {
                i__1 = *m - *p;
                i__2 = *m - *p;
                cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &c__0, &work[1] , &c_n1, &childinfo);
                /* Computing MAX */
                i__1 = 1;
                i__2 = *m - *p; // , expr subst
                lorgqrmin = max(i__1,i__2);
                lorgqropt = (integer) work[1].r;
            }
            cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], & c_n1, &childinfo);
            lorglqmin = max(1,*q);
            lorglqopt = (integer) work[1].r;
            cbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &c__0, &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &c__0, &c__0, &c__0, &c__0, & c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &childinfo);
            lbbcsd = (integer) rwork[1];
        }
        else if (r__ == *m - *p)
        {
            cunbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &work[1], & c_n1, &childinfo);
            lorbdb = (integer) work[1].r;
            if (*p >= *m - *p - 1)
            {
                cungqr_(p, p, q, &u1[u1_offset], ldu1, &c__0, &work[1], &c_n1, &childinfo);
                lorgqrmin = max(1,*p);
                lorgqropt = (integer) work[1].r;
            }
            else
            {
                i__1 = *m - *p - 1;
                i__2 = *m - *p - 1;
                i__3 = *m - *p - 1;
                cungqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, & c__0, &work[1], &c_n1, &childinfo);
                /* Computing MAX */
                i__1 = 1;
                i__2 = *m - *p - 1; // , expr subst
                lorgqrmin = max(i__1,i__2);
                lorgqropt = (integer) work[1].r;
            }
            cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &c__0, &work[1], & c_n1, &childinfo);
            lorglqmin = max(1,*q);
            lorglqopt = (integer) work[1].r;
            i__1 = *m - *q;
            i__2 = *m - *p;
            cbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1] , &c__0, &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[ u2_offset], ldu2, &u1[u1_offset], ldu1, &c__0, &c__0, & c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, &childinfo);
            lbbcsd = (integer) rwork[1];
        }
        else
        {
            cunbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &theta[1], &c__0, &c__0, &c__0, &c__0, &c__0, & work[1], &c_n1, &childinfo);
            lorbdb = *m + (integer) work[1].r;
            if (*p >= *m - *p)
            {
                i__1 = *m - *q;
                cungqr_(p, p, &i__1, &u1[u1_offset], ldu1, &c__0, &work[1], & c_n1, &childinfo);
                lorgqrmin = max(1,*p);
                lorgqropt = (integer) work[1].r;
            }
            else
            {
                i__1 = *m - *p;
                i__2 = *m - *p;
                i__3 = *m - *q;
                cungqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &c__0, & work[1], &c_n1, &childinfo);
                /* Computing MAX */
                i__1 = 1;
                i__2 = *m - *p; // , expr subst
                lorgqrmin = max(i__1,i__2);
                lorgqropt = (integer) work[1].r;
            }
            cunglq_(q, q, q, &v1t[v1t_offset], ldv1t, &c__0, &work[1], &c_n1, &childinfo);
            lorglqmin = max(1,*q);
            lorglqopt = (integer) work[1].r;
            i__1 = *m - *p;
            i__2 = *m - *q;
            cbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1] , &c__0, &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, & c__0, &c__1, &v1t[v1t_offset], ldv1t, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &c__0, &rwork[1], &c_n1, & childinfo);
            lbbcsd = (integer) rwork[1];
        }
        lrworkmin = ibbcsd + lbbcsd - 1;
        lrworkopt = lrworkmin;
        rwork[1] = (real) lrworkopt;
        /* Computing MAX */
        i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqrmin - 1;
        i__1 = max( i__1,i__2);
        i__2 = iorglq + lorglqmin - 1; // ; expr subst
        lworkmin = max(i__1,i__2);
        /* Computing MAX */
        i__1 = iorbdb + lorbdb - 1, i__2 = iorgqr + lorgqropt - 1;
        i__1 = max( i__1,i__2);
        i__2 = iorglq + lorglqopt - 1; // ; expr subst
        lworkopt = max(i__1,i__2);
        work[1].r = (real) lworkopt;
        work[1].i = 0.f; // , expr subst
        if (*lwork < lworkmin && ! lquery)
        {
            *info = -19;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CUNCSD2BY1", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    lorgqr = *lwork - iorgqr + 1;
    lorglq = *lwork - iorglq + 1;
    /* Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q, */
    /* in which R = MIN(P,M-P,Q,M-Q) */
    if (r__ == *q)
    {
        /* Case 1: R = Q */
        /* Simultaneously bidiagonalize X11 and X21 */
        cunbdb1_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, & theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[ itauq1], &work[iorbdb], &lorbdb, &childinfo);
        /* Accumulate Householder reflectors */
        if (wantu1 && *p > 0)
        {
            clacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1);
            cungqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[ iorgqr], &lorgqr, &childinfo);
        }
        if (wantu2 && *m - *p > 0)
        {
            i__1 = *m - *p;
            clacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], ldu2);
            i__1 = *m - *p;
            i__2 = *m - *p;
            cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], & work[iorgqr], &lorgqr, &childinfo);
        }
        if (wantv1t && *q > 0)
        {
            i__1 = v1t_dim1 + 1;
            v1t[i__1].r = 1.f;
            v1t[i__1].i = 0.f; // , expr subst
            i__1 = *q;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                i__2 = j * v1t_dim1 + 1;
                v1t[i__2].r = 0.f;
                v1t[i__2].i = 0.f; // , expr subst
                i__2 = j + v1t_dim1;
                v1t[i__2].r = 0.f;
                v1t[i__2].i = 0.f; // , expr subst
            }
            i__1 = *q - 1;
            i__2 = *q - 1;
            clacpy_("U", &i__1, &i__2, &x21[(x21_dim1 << 1) + 1], ldx21, &v1t[ (v1t_dim1 << 1) + 2], ldv1t);
            i__1 = *q - 1;
            i__2 = *q - 1;
            i__3 = *q - 1;
            cunglq_(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, & work[itauq1], &work[iorglq], &lorglq, &childinfo);
        }
        /* Simultaneously diagonalize X11 and X21. */
        cbbcsd_(jobu1, jobu2, jobv1t, "N", "N", m, p, q, &theta[1], &rwork[ iphi], &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[ v1t_offset], ldv1t, &c__0, &c__1, &rwork[ib11d], &rwork[ib11e] , &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, & childinfo);
        /* Permute rows and columns to place zero submatrices in */
        /* preferred positions */
        if (*q > 0 && wantu2)
        {
            i__1 = *q;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                iwork[i__] = *m - *p - *q + i__;
            }
            i__1 = *m - *p;
            for (i__ = *q + 1;
                    i__ <= i__1;
                    ++i__)
            {
                iwork[i__] = i__ - *q;
            }
            i__1 = *m - *p;
            i__2 = *m - *p;
            clapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
        }
    }
    else if (r__ == *p)
    {
        /* Case 2: R = P */
        /* Simultaneously bidiagonalize X11 and X21 */
        cunbdb2_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, & theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[ itauq1], &work[iorbdb], &lorbdb, &childinfo);
        /* Accumulate Householder reflectors */
        if (wantu1 && *p > 0)
        {
            i__1 = u1_dim1 + 1;
            u1[i__1].r = 1.f;
            u1[i__1].i = 0.f; // , expr subst
            i__1 = *p;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                i__2 = j * u1_dim1 + 1;
                u1[i__2].r = 0.f;
                u1[i__2].i = 0.f; // , expr subst
                i__2 = j + u1_dim1;
                u1[i__2].r = 0.f;
                u1[i__2].i = 0.f; // , expr subst
            }
            i__1 = *p - 1;
            i__2 = *p - 1;
            clacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[( u1_dim1 << 1) + 2], ldu1);
            i__1 = *p - 1;
            i__2 = *p - 1;
            i__3 = *p - 1;
            cungqr_(&i__1, &i__2, &i__3, &u1[(u1_dim1 << 1) + 2], ldu1, &work[ itaup1], &work[iorgqr], &lorgqr, &childinfo);
        }
        if (wantu2 && *m - *p > 0)
        {
            i__1 = *m - *p;
            clacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], ldu2);
            i__1 = *m - *p;
            i__2 = *m - *p;
            cungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], & work[iorgqr], &lorgqr, &childinfo);
        }
        if (wantv1t && *q > 0)
        {
            clacpy_("U", p, q, &x11[x11_offset], ldx11, &v1t[v1t_offset], ldv1t);
            cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[ iorglq], &lorglq, &childinfo);
        }
        /* Simultaneously diagonalize X11 and X21. */
        cbbcsd_(jobv1t, "N", jobu1, jobu2, "T", m, q, p, &theta[1], &rwork[ iphi], &v1t[v1t_offset], ldv1t, &c__0, &c__1, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &rwork[ib11d], &rwork[ib11e], & rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ib21e], & rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, & childinfo);
        /* Permute rows and columns to place identity submatrices in */
        /* preferred positions */
        if (*q > 0 && wantu2)
        {
            i__1 = *q;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                iwork[i__] = *m - *p - *q + i__;
            }
            i__1 = *m - *p;
            for (i__ = *q + 1;
                    i__ <= i__1;
                    ++i__)
            {
                iwork[i__] = i__ - *q;
            }
            i__1 = *m - *p;
            i__2 = *m - *p;
            clapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
        }
    }
    else if (r__ == *m - *p)
    {
        /* Case 3: R = M-P */
        /* Simultaneously bidiagonalize X11 and X21 */
        cunbdb3_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, & theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[ itauq1], &work[iorbdb], &lorbdb, &childinfo);
        /* Accumulate Householder reflectors */
        if (wantu1 && *p > 0)
        {
            clacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1);
            cungqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[ iorgqr], &lorgqr, &childinfo);
        }
        if (wantu2 && *m - *p > 0)
        {
            i__1 = u2_dim1 + 1;
            u2[i__1].r = 1.f;
            u2[i__1].i = 0.f; // , expr subst
            i__1 = *m - *p;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                i__2 = j * u2_dim1 + 1;
                u2[i__2].r = 0.f;
                u2[i__2].i = 0.f; // , expr subst
                i__2 = j + u2_dim1;
                u2[i__2].r = 0.f;
                u2[i__2].i = 0.f; // , expr subst
            }
            i__1 = *m - *p - 1;
            i__2 = *m - *p - 1;
            clacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[( u2_dim1 << 1) + 2], ldu2);
            i__1 = *m - *p - 1;
            i__2 = *m - *p - 1;
            i__3 = *m - *p - 1;
            cungqr_(&i__1, &i__2, &i__3, &u2[(u2_dim1 << 1) + 2], ldu2, &work[ itaup2], &work[iorgqr], &lorgqr, &childinfo);
        }
        if (wantv1t && *q > 0)
        {
            i__1 = *m - *p;
            clacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], ldv1t);
            cunglq_(q, q, &r__, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[ iorglq], &lorglq, &childinfo);
        }
        /* Simultaneously diagonalize X11 and X21. */
        i__1 = *m - *q;
        i__2 = *m - *p;
        cbbcsd_("N", jobv1t, jobu2, jobu1, "T", m, &i__1, &i__2, &theta[1], & rwork[iphi], &c__0, &c__1, &v1t[v1t_offset], ldv1t, &u2[ u2_offset], ldu2, &u1[u1_offset], ldu1, &rwork[ib11d], &rwork[ ib11e], &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, &childinfo);
        /* Permute rows and columns to place identity submatrices in */
        /* preferred positions */
        if (*q > r__)
        {
            i__1 = r__;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                iwork[i__] = *q - r__ + i__;
            }
            i__1 = *q;
            for (i__ = r__ + 1;
                    i__ <= i__1;
                    ++i__)
            {
                iwork[i__] = i__ - r__;
            }
            if (wantu1)
            {
                clapmt_(&c_false, p, q, &u1[u1_offset], ldu1, &iwork[1]);
            }
            if (wantv1t)
            {
                clapmr_(&c_false, q, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
            }
        }
    }
    else
    {
        /* Case 4: R = M-Q */
        /* Simultaneously bidiagonalize X11 and X21 */
        i__1 = lorbdb - *m;
        cunbdb4_(m, p, q, &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, & theta[1], &rwork[iphi], &work[itaup1], &work[itaup2], &work[ itauq1], &work[iorbdb], &work[iorbdb + *m], &i__1, &childinfo) ;
        /* Accumulate Householder reflectors */
        if (wantu1 && *p > 0)
        {
            ccopy_(p, &work[iorbdb], &c__1, &u1[u1_offset], &c__1);
            i__1 = *p;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                i__2 = j * u1_dim1 + 1;
                u1[i__2].r = 0.f;
                u1[i__2].i = 0.f; // , expr subst
            }
            i__1 = *p - 1;
            i__2 = *m - *q - 1;
            clacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &u1[( u1_dim1 << 1) + 2], ldu1);
            i__1 = *m - *q;
            cungqr_(p, p, &i__1, &u1[u1_offset], ldu1, &work[itaup1], &work[ iorgqr], &lorgqr, &childinfo);
        }
        if (wantu2 && *m - *p > 0)
        {
            i__1 = *m - *p;
            ccopy_(&i__1, &work[iorbdb + *p], &c__1, &u2[u2_offset], &c__1);
            i__1 = *m - *p;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                i__2 = j * u2_dim1 + 1;
                u2[i__2].r = 0.f;
                u2[i__2].i = 0.f; // , expr subst
            }
            i__1 = *m - *p - 1;
            i__2 = *m - *q - 1;
            clacpy_("L", &i__1, &i__2, &x21[x21_dim1 + 2], ldx21, &u2[( u2_dim1 << 1) + 2], ldu2);
            i__1 = *m - *p;
            i__2 = *m - *p;
            i__3 = *m - *q;
            cungqr_(&i__1, &i__2, &i__3, &u2[u2_offset], ldu2, &work[itaup2], &work[iorgqr], &lorgqr, &childinfo);
        }
        if (wantv1t && *q > 0)
        {
            i__1 = *m - *q;
            clacpy_("U", &i__1, q, &x21[x21_offset], ldx21, &v1t[v1t_offset], ldv1t);
            i__1 = *p - (*m - *q);
            i__2 = *q - (*m - *q);
            clacpy_("U", &i__1, &i__2, &x11[*m - *q + 1 + (*m - *q + 1) * x11_dim1], ldx11, &v1t[*m - *q + 1 + (*m - *q + 1) * v1t_dim1], ldv1t);
            i__1 = -(*p) + *q;
            i__2 = *q - *p;
            clacpy_("U", &i__1, &i__2, &x21[*m - *q + 1 + (*p + 1) * x21_dim1] , ldx21, &v1t[*p + 1 + (*p + 1) * v1t_dim1], ldv1t);
            cunglq_(q, q, q, &v1t[v1t_offset], ldv1t, &work[itauq1], &work[ iorglq], &lorglq, &childinfo);
        }
        /* Simultaneously diagonalize X11 and X21. */
        i__1 = *m - *p;
        i__2 = *m - *q;
        cbbcsd_(jobu2, jobu1, "N", jobv1t, "N", m, &i__1, &i__2, &theta[1], & rwork[iphi], &u2[u2_offset], ldu2, &u1[u1_offset], ldu1, & c__0, &c__1, &v1t[v1t_offset], ldv1t, &rwork[ib11d], &rwork[ ib11e], &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsd, &childinfo);
        /* Permute rows and columns to place identity submatrices in */
        /* preferred positions */
        if (*p > r__)
        {
            i__1 = r__;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                iwork[i__] = *p - r__ + i__;
            }
            i__1 = *p;
            for (i__ = r__ + 1;
                    i__ <= i__1;
                    ++i__)
            {
                iwork[i__] = i__ - r__;
            }
            if (wantu1)
            {
                clapmt_(&c_false, p, p, &u1[u1_offset], ldu1, &iwork[1]);
            }
            if (wantv1t)
            {
                clapmr_(&c_false, p, q, &v1t[v1t_offset], ldv1t, &iwork[1]);
            }
        }
    }
    return 0;
    /* End of CUNCSD2BY1 */
}
/* cuncsd2by1_ */
