/* ../netlib/dorbdb1.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DORBDB1 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DORBDB1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorbdb1 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorbdb1 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorbdb1 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DORBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, */
/* TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LWORK, M, P, Q, LDX11, LDX21 */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION PHI(*), THETA(*) */
/* DOUBLE PRECISION TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), */
/* $ X11(LDX11,*), X21(LDX21,*) */
/* .. */
/* > \par Purpose: */
/* > ============= */
/* > */
/* >\verbatim */
/* > */
/* > DORBDB1 simultaneously bidiagonalizes the blocks of a tall and skinny */
/* > matrix X with orthonomal columns: */
/* > */
/* > [ B11 ] */
/* > [ X11 ] [ P1 | ] [ 0 ] */
/* > [-----] = [---------] [-----] Q1**T . */
/* > [ X21 ] [ | P2 ] [ B21 ] */
/* > [ 0 ] */
/* > */
/* > X11 is P-by-Q, and X21 is (M-P)-by-Q. Q must be no larger than P, */
/* > M-P, or M-Q. Routines DORBDB2, DORBDB3, and DORBDB4 handle cases in */
/* > which Q is not the minimum dimension. */
/* > */
/* > The orthogonal matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P), */
/* > and (M-Q)-by-(M-Q), respectively. They are represented implicitly by */
/* > Householder vectors. */
/* > */
/* > B11 and B12 are Q-by-Q bidiagonal matrices represented implicitly by */
/* > angles THETA, PHI. */
/* > */
/* >\endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows X11 plus the number of rows in X21. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* > P is INTEGER */
/* > The number of rows in X11. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* > Q is INTEGER */
/* > The number of columns in X11 and X21. 0 <= Q <= */
/* > MIN(P,M-P,M-Q). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* > X11 is DOUBLE PRECISION array, dimension (LDX11,Q) */
/* > On entry, the top block of the matrix X to be reduced. On */
/* > exit, the columns of tril(X11) specify reflectors for P1 and */
/* > the rows of triu(X11,1) specify reflectors for Q1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* > LDX11 is INTEGER */
/* > The leading dimension of X11. LDX11 >= P. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* > X21 is DOUBLE PRECISION array, dimension (LDX21,Q) */
/* > On entry, the bottom block of the matrix X to be reduced. On */
/* > exit, the columns of tril(X21) specify reflectors for P2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* > LDX21 is INTEGER */
/* > The leading dimension of X21. LDX21 >= M-P. */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* > THETA is DOUBLE PRECISION array, dimension (Q) */
/* > The entries of the bidiagonal blocks B11, B21 are defined by */
/* > THETA and PHI. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] PHI */
/* > \verbatim */
/* > PHI is DOUBLE PRECISION array, dimension (Q-1) */
/* > The entries of the bidiagonal blocks B11, B21 are defined by */
/* > THETA and PHI. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP1 */
/* > \verbatim */
/* > TAUP1 is DOUBLE PRECISION array, dimension (P) */
/* > The scalar factors of the elementary reflectors that define */
/* > P1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP2 */
/* > \verbatim */
/* > TAUP2 is DOUBLE PRECISION array, dimension (M-P) */
/* > The scalar factors of the elementary reflectors that define */
/* > P2. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ1 */
/* > \verbatim */
/* > TAUQ1 is DOUBLE PRECISION array, dimension (Q) */
/* > The scalar factors of the elementary reflectors that define */
/* > Q1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= M-Q. */
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
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* > */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date July 2012 */
/* > \ingroup doubleOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The upper-bidiagonal blocks B11, B21 are represented implicitly by */
/* > angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry */
/* > in each bidiagonal band is a product of a sine or cosine of a THETA */
/* > with a sine or cosine of a PHI. See [1] or DORCSD for details. */
/* > */
/* > P1, P2, and Q1 are represented as products of elementary reflectors. */
/* > See DORCSD2BY1 for details on generating P1, P2, and Q1 using DORGQR */
/* > and DORGLQ. */
/* > \endverbatim */
/* > \par References: */
/* ================ */
/* > */
/* > [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* > Algorithms, 50(1):33-65, 2009. */
/* > */
/* ===================================================================== */
/* Subroutine */
int dorbdb1_(integer *m, integer *p, integer *q, doublereal * x11, integer *ldx11, doublereal *x21, integer *ldx21, doublereal * theta, doublereal *phi, doublereal *taup1, doublereal *taup2, doublereal *tauq1, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer x11_dim1, x11_offset, x21_dim1, x21_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    /* Builtin functions */
    double atan2(doublereal, doublereal), cos(doublereal), sin(doublereal), sqrt(doublereal);
    /* Local variables */
    integer lworkmin, lworkopt;
    doublereal c__;
    integer i__;
    doublereal s;
    integer childinfo;
    extern /* Subroutine */
    int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *);
    extern doublereal dnrm2_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern /* Subroutine */
    int dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *);
    integer ilarf, llarf;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    logical lquery;
    extern /* Subroutine */
    int dorbdb5_(integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *);
    integer iorbdb5, lorbdb5;
    extern /* Subroutine */
    int dlarfgp_(integer *, doublereal *, doublereal * , integer *, doublereal *);
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* July 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ==================================================================== */
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
    --phi;
    --taup1;
    --taup2;
    --tauq1;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*p < *q || *m - *p < *q)
    {
        *info = -2;
    }
    else if (*q < 0 || *m - *q < *q)
    {
        *info = -3;
    }
    else if (*ldx11 < max(1,*p))
    {
        *info = -5;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *p; // , expr subst
        if (*ldx21 < max(i__1,i__2))
        {
            *info = -7;
        }
    }
    /* Compute workspace */
    if (*info == 0)
    {
        ilarf = 2;
        /* Computing MAX */
        i__1 = *p - 1, i__2 = *m - *p - 1;
        i__1 = max(i__1,i__2);
        i__2 = *q - 1; // ; expr subst
        llarf = max(i__1,i__2);
        iorbdb5 = 2;
        lorbdb5 = *q - 2;
        /* Computing MAX */
        i__1 = ilarf + llarf - 1;
        i__2 = iorbdb5 + lorbdb5 - 1; // , expr subst
        lworkopt = max(i__1,i__2);
        lworkmin = lworkopt;
        work[1] = (doublereal) lworkopt;
        if (*lwork < lworkmin && ! lquery)
        {
            *info = -14;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DORBDB1", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Reduce columns 1, ..., Q of X11 and X21 */
    i__1 = *q;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = *p - i__ + 1;
        dlarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + 1 + i__ * x11_dim1], &c__1, &taup1[i__]);
        i__2 = *m - *p - i__ + 1;
        dlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + 1 + i__ * x21_dim1], &c__1, &taup2[i__]);
        theta[i__] = atan2(x21[i__ + i__ * x21_dim1], x11[i__ + i__ * x11_dim1]);
        c__ = cos(theta[i__]);
        s = sin(theta[i__]);
        x11[i__ + i__ * x11_dim1] = 1.;
        x21[i__ + i__ * x21_dim1] = 1.;
        i__2 = *p - i__ + 1;
        i__3 = *q - i__;
        dlarf_("L", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], &c__1, &taup1[ i__], &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &work[ilarf]);
        i__2 = *m - *p - i__ + 1;
        i__3 = *q - i__;
        dlarf_("L", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], &c__1, &taup2[ i__], &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &work[ilarf]);
        if (i__ < *q)
        {
            i__2 = *q - i__;
            drot_(&i__2, &x11[i__ + (i__ + 1) * x11_dim1], ldx11, &x21[i__ + ( i__ + 1) * x21_dim1], ldx21, &c__, &s);
            i__2 = *q - i__;
            dlarfgp_(&i__2, &x21[i__ + (i__ + 1) * x21_dim1], &x21[i__ + (i__ + 2) * x21_dim1], ldx21, &tauq1[i__]);
            s = x21[i__ + (i__ + 1) * x21_dim1];
            x21[i__ + (i__ + 1) * x21_dim1] = 1.;
            i__2 = *p - i__;
            i__3 = *q - i__;
            dlarf_("R", &i__2, &i__3, &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &tauq1[i__], &x11[i__ + 1 + (i__ + 1) * x11_dim1], ldx11, &work[ilarf]);
            i__2 = *m - *p - i__;
            i__3 = *q - i__;
            dlarf_("R", &i__2, &i__3, &x21[i__ + (i__ + 1) * x21_dim1], ldx21, &tauq1[i__], &x21[i__ + 1 + (i__ + 1) * x21_dim1], ldx21, &work[ilarf]);
            i__2 = *p - i__;
            /* Computing 2nd power */
            d__1 = dnrm2_(&i__2, &x11[i__ + 1 + (i__ + 1) * x11_dim1], &c__1, &x11[i__ + 1 + (i__ + 1) * x11_dim1], &c__1);
            i__3 = *m - *p - i__;
            /* Computing 2nd power */
            d__2 = dnrm2_(&i__3, &x21[i__ + 1 + (i__ + 1) * x21_dim1], &c__1, &x21[i__ + 1 + (i__ + 1) * x21_dim1], &c__1);
            c__ = sqrt(d__1 * d__1 + d__2 * d__2);
            phi[i__] = atan2(s, c__);
            i__2 = *p - i__;
            i__3 = *m - *p - i__;
            i__4 = *q - i__ - 1;
            dorbdb5_(&i__2, &i__3, &i__4, &x11[i__ + 1 + (i__ + 1) * x11_dim1] , &c__1, &x21[i__ + 1 + (i__ + 1) * x21_dim1], &c__1, & x11[i__ + 1 + (i__ + 2) * x11_dim1], ldx11, &x21[i__ + 1 + (i__ + 2) * x21_dim1], ldx21, &work[iorbdb5], &lorbdb5, &childinfo);
        }
    }
    return 0;
    /* End of DORBDB1 */
}
/* dorbdb1_ */
