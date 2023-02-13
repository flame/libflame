/* ../netlib/claqr5.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
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
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__2 = 2;
/* > \brief \b CLAQR5 performs a single small-bulge multi-shift QR sweep. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQR5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr5. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr5. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr5. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S, */
/* H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV, */
/* WV, LDWV, NH, WH, LDWH ) */
/* .. Scalar Arguments .. */
/* INTEGER IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, */
/* $ LDWH, LDWV, LDZ, N, NH, NSHFTS, NV */
/* LOGICAL WANTT, WANTZ */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ), */
/* $ WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQR5 called by CLAQR0 performs a */
/* > single small-bulge multi-shift QR sweep. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTT */
/* > \verbatim */
/* > WANTT is logical scalar */
/* > WANTT = .true. if the triangular Schur factor */
/* > is being computed. WANTT is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is logical scalar */
/* > WANTZ = .true. if the unitary Schur factor is being */
/* > computed. WANTZ is set to .false. otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in] KACC22 */
/* > \verbatim */
/* > KACC22 is integer with value 0, 1, or 2. */
/* > Specifies the computation mode of far-from-diagonal */
/* > orthogonal updates. */
/* > = 0: CLAQR5 does not accumulate reflections and does not */
/* > use matrix-matrix multiply to update far-from-diagonal */
/* > matrix entries. */
/* > = 1: CLAQR5 accumulates reflections and uses matrix-matrix */
/* > multiply to update the far-from-diagonal matrix entries. */
/* > = 2: CLAQR5 accumulates reflections, uses matrix-matrix */
/* > multiply to update the far-from-diagonal matrix entries, */
/* > and takes advantage of 2-by-2 block structure during */
/* > matrix multiplies. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is integer scalar */
/* > N is the order of the Hessenberg matrix H upon which this */
/* > subroutine operates. */
/* > \endverbatim */
/* > */
/* > \param[in] KTOP */
/* > \verbatim */
/* > KTOP is integer scalar */
/* > \endverbatim */
/* > */
/* > \param[in] KBOT */
/* > \verbatim */
/* > KBOT is integer scalar */
/* > These are the first and last rows and columns of an */
/* > isolated diagonal block upon which the QR sweep is to be */
/* > applied. It is assumed without a check that */
/* > either KTOP = 1 or H(KTOP,KTOP-1) = 0 */
/* > and */
/* > either KBOT = N or H(KBOT+1,KBOT) = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NSHFTS */
/* > \verbatim */
/* > NSHFTS is integer scalar */
/* > NSHFTS gives the number of simultaneous shifts. NSHFTS */
/* > must be positive and even. */
/* > \endverbatim */
/* > */
/* > \param[in,out] S */
/* > \verbatim */
/* > S is COMPLEX array of size (NSHFTS) */
/* > S contains the shifts of origin that define the multi- */
/* > shift QR sweep. On output S may be reordered. */
/* > \endverbatim */
/* > */
/* > \param[in,out] H */
/* > \verbatim */
/* > H is COMPLEX array of size (LDH,N) */
/* > On input H contains a Hessenberg matrix. On output a */
/* > multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied */
/* > to the isolated diagonal block in rows and columns KTOP */
/* > through KBOT. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is integer scalar */
/* > LDH is the leading dimension of H just as declared in the */
/* > calling procedure. LDH.GE.MAX(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] ILOZ */
/* > \verbatim */
/* > ILOZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHIZ */
/* > \verbatim */
/* > IHIZ is INTEGER */
/* > Specify the rows of Z to which transformations must be */
/* > applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is COMPLEX array of size (LDZ,IHI) */
/* > If WANTZ = .TRUE., then the QR Sweep unitary */
/* > similarity transformation is accumulated into */
/* > Z(ILOZ:IHIZ,ILO:IHI) from the right. */
/* > If WANTZ = .FALSE., then Z is unreferenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is integer scalar */
/* > LDA is the leading dimension of Z just as declared in */
/* > the calling procedure. LDZ.GE.N. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX array of size (LDV,NSHFTS/2) */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is integer scalar */
/* > LDV is the leading dimension of V as declared in the */
/* > calling procedure. LDV.GE.3. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* > U is COMPLEX array of size */
/* > (LDU,3*NSHFTS-3) */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is integer scalar */
/* > LDU is the leading dimension of U just as declared in the */
/* > in the calling subroutine. LDU.GE.3*NSHFTS-3. */
/* > \endverbatim */
/* > */
/* > \param[in] NH */
/* > \verbatim */
/* > NH is integer scalar */
/* > NH is the number of columns in array WH available for */
/* > workspace. NH.GE.1. */
/* > \endverbatim */
/* > */
/* > \param[out] WH */
/* > \verbatim */
/* > WH is COMPLEX array of size (LDWH,NH) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWH */
/* > \verbatim */
/* > LDWH is integer scalar */
/* > Leading dimension of WH just as declared in the */
/* > calling procedure. LDWH.GE.3*NSHFTS-3. */
/* > \endverbatim */
/* > */
/* > \param[in] NV */
/* > \verbatim */
/* > NV is integer scalar */
/* > NV is the number of rows in WV agailable for workspace. */
/* > NV.GE.1. */
/* > \endverbatim */
/* > */
/* > \param[out] WV */
/* > \verbatim */
/* > WV is COMPLEX array of size */
/* > (LDWV,3*NSHFTS-3) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWV */
/* > \verbatim */
/* > LDWV is integer scalar */
/* > LDWV is the leading dimension of WV as declared in the */
/* > in the calling subroutine. LDWV.GE.NV. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Karen Braman and Ralph Byers, Department of Mathematics, */
/* > University of Kansas, USA */
/* > \par References: */
/* ================ */
/* > */
/* > K. Braman, R. Byers and R. Mathias, The Multi-Shift QR */
/* > Algorithm Part I: Maintaining Well Focused Shifts, and Level 3 */
/* > Performance, SIAM Journal of Matrix Analysis, volume 23, pages */
/* > 929--947, 2002. */
/* > */
/* ===================================================================== */
/* Subroutine */
int claqr5_(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop, integer *kbot, integer *nshfts, complex *s, complex *h__, integer *ldh, integer *iloz, integer *ihiz, complex * z__, integer *ldz, complex *v, integer *ldv, complex *u, integer *ldu, integer *nv, complex *wv, integer *ldwv, integer *nh, complex *wh, integer *ldwh)
{
    /* System generated locals */
    integer h_dim1, h_offset, u_dim1, u_offset, v_dim1, v_offset, wh_dim1, wh_offset, wv_dim1, wv_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, i__11;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7, q__8;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    double r_imag(complex *);
    /* Local variables */
    integer j, k, m, i2, j2, i4, j4, k1;
    real h11, h12, h21, h22;
    integer m22, ns, nu;
    complex vt[3];
    real scl;
    integer kdu, kms;
    real ulp;
    integer knz, kzs;
    real tst1, tst2;
    complex beta;
    logical blk22, bmp22;
    integer mend, jcol, jlen, jbot, mbot, jtop, jrow, mtop;
    complex alpha;
    logical accum;
    extern /* Subroutine */
    int cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    integer ndcol, incol, krcol, nbmps;
    extern /* Subroutine */
    int ctrmm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *), claqr1_(integer *, complex *, integer *, complex *, complex *, complex *), slabad_( real *, real *), clarfg_(integer *, complex *, complex *, integer *, complex *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *);
    real safmin, safmax;
    complex refsum;
    integer mstart;
    real smlnum;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ================================================================ */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* ==== If there are no shifts, then there is nothing to do. ==== */
    /* Parameter adjustments */
    --s;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    wv_dim1 = *ldwv;
    wv_offset = 1 + wv_dim1;
    wv -= wv_offset;
    wh_dim1 = *ldwh;
    wh_offset = 1 + wh_dim1;
    wh -= wh_offset;
    /* Function Body */
    if (*nshfts < 2)
    {
        return 0;
    }
    /* ==== If the active block is empty or 1-by-1, then there */
    /* . is nothing to do. ==== */
    if (*ktop >= *kbot)
    {
        return 0;
    }
    /* ==== NSHFTS is supposed to be even, but if it is odd, */
    /* . then simply reduce it by one. ==== */
    ns = *nshfts - *nshfts % 2;
    /* ==== Machine constants for deflation ==== */
    safmin = slamch_("SAFE MINIMUM");
    safmax = 1.f / safmin;
    slabad_(&safmin, &safmax);
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real) (*n) / ulp);
    /* ==== Use accumulated reflections to update far-from-diagonal */
    /* . entries ? ==== */
    accum = *kacc22 == 1 || *kacc22 == 2;
    /* ==== If so, exploit the 2-by-2 block structure? ==== */
    blk22 = ns > 2 && *kacc22 == 2;
    /* ==== clear trash ==== */
    if (*ktop + 2 <= *kbot)
    {
        i__1 = *ktop + 2 + *ktop * h_dim1;
        h__[i__1].r = 0.f;
        h__[i__1].i = 0.f; // , expr subst
    }
    /* ==== NBMPS = number of 2-shift bulges in the chain ==== */
    nbmps = ns / 2;
    /* ==== KDU = width of slab ==== */
    kdu = nbmps * 6 - 3;
    /* ==== Create and chase chains of NBMPS bulges ==== */
    i__1 = *kbot - 2;
    i__2 = nbmps * 3 - 2;
    for (incol = (1 - nbmps) * 3 + *ktop - 1;
            i__2 < 0 ? incol >= i__1 : incol <= i__1;
            incol += i__2)
    {
        ndcol = incol + kdu;
        if (accum)
        {
            claset_("ALL", &kdu, &kdu, &c_b1, &c_b2, &u[u_offset], ldu);
        }
        /* ==== Near-the-diagonal bulge chase. The following loop */
        /* . performs the near-the-diagonal part of a small bulge */
        /* . multi-shift QR sweep. Each 6*NBMPS-2 column diagonal */
        /* . chunk extends from column INCOL to column NDCOL */
        /* . (including both column INCOL and column NDCOL). The */
        /* . following loop chases a 3*NBMPS column long chain of */
        /* . NBMPS bulges 3*NBMPS-2 columns to the right. (INCOL */
        /* . may be less than KTOP and and NDCOL may be greater than */
        /* . KBOT indicating phantom columns from which to chase */
        /* . bulges before they are actually introduced or to which */
        /* . to chase bulges beyond column KBOT.) ==== */
        /* Computing MIN */
        i__4 = incol + nbmps * 3 - 3;
        i__5 = *kbot - 2; // , expr subst
        i__3 = min(i__4,i__5);
        for (krcol = incol;
                krcol <= i__3;
                ++krcol)
        {
            /* ==== Bulges number MTOP to MBOT are active double implicit */
            /* . shift bulges. There may or may not also be small */
            /* . 2-by-2 bulge, if there is room. The inactive bulges */
            /* . (if any) must wait until the active bulges have moved */
            /* . down the diagonal to make room. The phantom matrix */
            /* . paradigm described above helps keep track. ==== */
            /* Computing MAX */
            i__4 = 1;
            i__5 = (*ktop - 1 - krcol + 2) / 3 + 1; // , expr subst
            mtop = max(i__4,i__5);
            /* Computing MIN */
            i__4 = nbmps;
            i__5 = (*kbot - krcol) / 3; // , expr subst
            mbot = min(i__4,i__5);
            m22 = mbot + 1;
            bmp22 = mbot < nbmps && krcol + (m22 - 1) * 3 == *kbot - 2;
            /* ==== Generate reflections to chase the chain right */
            /* . one column. (The minimum value of K is KTOP-1.) ==== */
            i__4 = mbot;
            for (m = mtop;
                    m <= i__4;
                    ++m)
            {
                k = krcol + (m - 1) * 3;
                if (k == *ktop - 1)
                {
                    claqr1_(&c__3, &h__[*ktop + *ktop * h_dim1], ldh, &s[(m << 1) - 1], &s[m * 2], &v[m * v_dim1 + 1]);
                    i__5 = m * v_dim1 + 1;
                    alpha.r = v[i__5].r;
                    alpha.i = v[i__5].i; // , expr subst
                    clarfg_(&c__3, &alpha, &v[m * v_dim1 + 2], &c__1, &v[m * v_dim1 + 1]);
                }
                else
                {
                    i__5 = k + 1 + k * h_dim1;
                    beta.r = h__[i__5].r;
                    beta.i = h__[i__5].i; // , expr subst
                    i__5 = m * v_dim1 + 2;
                    i__6 = k + 2 + k * h_dim1;
                    v[i__5].r = h__[i__6].r;
                    v[i__5].i = h__[i__6].i; // , expr subst
                    i__5 = m * v_dim1 + 3;
                    i__6 = k + 3 + k * h_dim1;
                    v[i__5].r = h__[i__6].r;
                    v[i__5].i = h__[i__6].i; // , expr subst
                    clarfg_(&c__3, &beta, &v[m * v_dim1 + 2], &c__1, &v[m * v_dim1 + 1]);
                    /* ==== A Bulge may collapse because of vigilant */
                    /* . deflation or destructive underflow. In the */
                    /* . underflow case, try the two-small-subdiagonals */
                    /* . trick to try to reinflate the bulge. ==== */
                    i__5 = k + 3 + k * h_dim1;
                    i__6 = k + 3 + (k + 1) * h_dim1;
                    i__7 = k + 3 + (k + 2) * h_dim1;
                    if (h__[i__5].r != 0.f || h__[i__5].i != 0.f || (h__[i__6] .r != 0.f || h__[i__6].i != 0.f) || h__[i__7].r == 0.f && h__[i__7].i == 0.f)
                    {
                        /* ==== Typical case: not collapsed (yet). ==== */
                        i__5 = k + 1 + k * h_dim1;
                        h__[i__5].r = beta.r;
                        h__[i__5].i = beta.i; // , expr subst
                        i__5 = k + 2 + k * h_dim1;
                        h__[i__5].r = 0.f;
                        h__[i__5].i = 0.f; // , expr subst
                        i__5 = k + 3 + k * h_dim1;
                        h__[i__5].r = 0.f;
                        h__[i__5].i = 0.f; // , expr subst
                    }
                    else
                    {
                        /* ==== Atypical case: collapsed. Attempt to */
                        /* . reintroduce ignoring H(K+1,K) and H(K+2,K). */
                        /* . If the fill resulting from the new */
                        /* . reflector is too large, then abandon it. */
                        /* . Otherwise, use the new one. ==== */
                        claqr1_(&c__3, &h__[k + 1 + (k + 1) * h_dim1], ldh, & s[(m << 1) - 1], &s[m * 2], vt);
                        alpha.r = vt[0].r;
                        alpha.i = vt[0].i; // , expr subst
                        clarfg_(&c__3, &alpha, &vt[1], &c__1, vt);
                        r_cnjg(&q__2, vt);
                        i__5 = k + 1 + k * h_dim1;
                        r_cnjg(&q__5, &vt[1]);
                        i__6 = k + 2 + k * h_dim1;
                        q__4.r = q__5.r * h__[i__6].r - q__5.i * h__[i__6].i;
                        q__4.i = q__5.r * h__[i__6].i + q__5.i * h__[ i__6].r; // , expr subst
                        q__3.r = h__[i__5].r + q__4.r;
                        q__3.i = h__[i__5].i + q__4.i; // , expr subst
                        q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                        q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                        refsum.r = q__1.r;
                        refsum.i = q__1.i; // , expr subst
                        i__5 = k + 2 + k * h_dim1;
                        q__3.r = refsum.r * vt[1].r - refsum.i * vt[1].i;
                        q__3.i = refsum.r * vt[1].i + refsum.i * vt[1] .r; // , expr subst
                        q__2.r = h__[i__5].r - q__3.r;
                        q__2.i = h__[i__5].i - q__3.i; // , expr subst
                        q__1.r = q__2.r;
                        q__1.i = q__2.i; // , expr subst
                        q__5.r = refsum.r * vt[2].r - refsum.i * vt[2].i;
                        q__5.i = refsum.r * vt[2].i + refsum.i * vt[2] .r; // , expr subst
                        q__4.r = q__5.r;
                        q__4.i = q__5.i; // , expr subst
                        i__6 = k + k * h_dim1;
                        i__7 = k + 1 + (k + 1) * h_dim1;
                        i__8 = k + 2 + (k + 2) * h_dim1;
                        if ((r__1 = q__1.r, f2c_abs(r__1)) + (r__2 = r_imag(&q__1) , f2c_abs(r__2)) + ((r__3 = q__4.r, f2c_abs(r__3)) + ( r__4 = r_imag(&q__4), f2c_abs(r__4))) > ulp * (( r__5 = h__[i__6].r, f2c_abs(r__5)) + (r__6 = r_imag(&h__[k + k * h_dim1]), f2c_abs(r__6)) + (( r__7 = h__[i__7].r, f2c_abs(r__7)) + (r__8 = r_imag(&h__[k + 1 + (k + 1) * h_dim1]), f2c_abs( r__8))) + ((r__9 = h__[i__8].r, f2c_abs(r__9)) + ( r__10 = r_imag(&h__[k + 2 + (k + 2) * h_dim1]) , f2c_abs(r__10)))))
                        {
                            /* ==== Starting a new bulge here would */
                            /* . create non-negligible fill. Use */
                            /* . the old one with trepidation. ==== */
                            i__5 = k + 1 + k * h_dim1;
                            h__[i__5].r = beta.r;
                            h__[i__5].i = beta.i; // , expr subst
                            i__5 = k + 2 + k * h_dim1;
                            h__[i__5].r = 0.f;
                            h__[i__5].i = 0.f; // , expr subst
                            i__5 = k + 3 + k * h_dim1;
                            h__[i__5].r = 0.f;
                            h__[i__5].i = 0.f; // , expr subst
                        }
                        else
                        {
                            /* ==== Stating a new bulge here would */
                            /* . create only negligible fill. */
                            /* . Replace the old reflector with */
                            /* . the new one. ==== */
                            i__5 = k + 1 + k * h_dim1;
                            i__6 = k + 1 + k * h_dim1;
                            q__1.r = h__[i__6].r - refsum.r;
                            q__1.i = h__[ i__6].i - refsum.i; // , expr subst
                            h__[i__5].r = q__1.r;
                            h__[i__5].i = q__1.i; // , expr subst
                            i__5 = k + 2 + k * h_dim1;
                            h__[i__5].r = 0.f;
                            h__[i__5].i = 0.f; // , expr subst
                            i__5 = k + 3 + k * h_dim1;
                            h__[i__5].r = 0.f;
                            h__[i__5].i = 0.f; // , expr subst
                            i__5 = m * v_dim1 + 1;
                            v[i__5].r = vt[0].r;
                            v[i__5].i = vt[0].i; // , expr subst
                            i__5 = m * v_dim1 + 2;
                            v[i__5].r = vt[1].r;
                            v[i__5].i = vt[1].i; // , expr subst
                            i__5 = m * v_dim1 + 3;
                            v[i__5].r = vt[2].r;
                            v[i__5].i = vt[2].i; // , expr subst
                        }
                    }
                }
                /* L10: */
            }
            /* ==== Generate a 2-by-2 reflection, if needed. ==== */
            k = krcol + (m22 - 1) * 3;
            if (bmp22)
            {
                if (k == *ktop - 1)
                {
                    claqr1_(&c__2, &h__[k + 1 + (k + 1) * h_dim1], ldh, &s[( m22 << 1) - 1], &s[m22 * 2], &v[m22 * v_dim1 + 1]) ;
                    i__4 = m22 * v_dim1 + 1;
                    beta.r = v[i__4].r;
                    beta.i = v[i__4].i; // , expr subst
                    clarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 * v_dim1 + 1]);
                }
                else
                {
                    i__4 = k + 1 + k * h_dim1;
                    beta.r = h__[i__4].r;
                    beta.i = h__[i__4].i; // , expr subst
                    i__4 = m22 * v_dim1 + 2;
                    i__5 = k + 2 + k * h_dim1;
                    v[i__4].r = h__[i__5].r;
                    v[i__4].i = h__[i__5].i; // , expr subst
                    clarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 * v_dim1 + 1]);
                    i__4 = k + 1 + k * h_dim1;
                    h__[i__4].r = beta.r;
                    h__[i__4].i = beta.i; // , expr subst
                    i__4 = k + 2 + k * h_dim1;
                    h__[i__4].r = 0.f;
                    h__[i__4].i = 0.f; // , expr subst
                }
            }
            /* ==== Multiply H by reflections from the left ==== */
            if (accum)
            {
                jbot = min(ndcol,*kbot);
            }
            else if (*wantt)
            {
                jbot = *n;
            }
            else
            {
                jbot = *kbot;
            }
            i__4 = jbot;
            for (j = max(*ktop,krcol);
                    j <= i__4;
                    ++j)
            {
                /* Computing MIN */
                i__5 = mbot;
                i__6 = (j - krcol + 2) / 3; // , expr subst
                mend = min(i__5,i__6);
                i__5 = mend;
                for (m = mtop;
                        m <= i__5;
                        ++m)
                {
                    k = krcol + (m - 1) * 3;
                    r_cnjg(&q__2, &v[m * v_dim1 + 1]);
                    i__6 = k + 1 + j * h_dim1;
                    r_cnjg(&q__6, &v[m * v_dim1 + 2]);
                    i__7 = k + 2 + j * h_dim1;
                    q__5.r = q__6.r * h__[i__7].r - q__6.i * h__[i__7].i;
                    q__5.i = q__6.r * h__[i__7].i + q__6.i * h__[i__7] .r; // , expr subst
                    q__4.r = h__[i__6].r + q__5.r;
                    q__4.i = h__[i__6].i + q__5.i; // , expr subst
                    r_cnjg(&q__8, &v[m * v_dim1 + 3]);
                    i__8 = k + 3 + j * h_dim1;
                    q__7.r = q__8.r * h__[i__8].r - q__8.i * h__[i__8].i;
                    q__7.i = q__8.r * h__[i__8].i + q__8.i * h__[i__8] .r; // , expr subst
                    q__3.r = q__4.r + q__7.r;
                    q__3.i = q__4.i + q__7.i; // , expr subst
                    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                    q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                    refsum.r = q__1.r;
                    refsum.i = q__1.i; // , expr subst
                    i__6 = k + 1 + j * h_dim1;
                    i__7 = k + 1 + j * h_dim1;
                    q__1.r = h__[i__7].r - refsum.r;
                    q__1.i = h__[i__7].i - refsum.i; // , expr subst
                    h__[i__6].r = q__1.r;
                    h__[i__6].i = q__1.i; // , expr subst
                    i__6 = k + 2 + j * h_dim1;
                    i__7 = k + 2 + j * h_dim1;
                    i__8 = m * v_dim1 + 2;
                    q__2.r = refsum.r * v[i__8].r - refsum.i * v[i__8].i;
                    q__2.i = refsum.r * v[i__8].i + refsum.i * v[i__8] .r; // , expr subst
                    q__1.r = h__[i__7].r - q__2.r;
                    q__1.i = h__[i__7].i - q__2.i; // , expr subst
                    h__[i__6].r = q__1.r;
                    h__[i__6].i = q__1.i; // , expr subst
                    i__6 = k + 3 + j * h_dim1;
                    i__7 = k + 3 + j * h_dim1;
                    i__8 = m * v_dim1 + 3;
                    q__2.r = refsum.r * v[i__8].r - refsum.i * v[i__8].i;
                    q__2.i = refsum.r * v[i__8].i + refsum.i * v[i__8] .r; // , expr subst
                    q__1.r = h__[i__7].r - q__2.r;
                    q__1.i = h__[i__7].i - q__2.i; // , expr subst
                    h__[i__6].r = q__1.r;
                    h__[i__6].i = q__1.i; // , expr subst
                    /* L20: */
                }
                /* L30: */
            }
            if (bmp22)
            {
                k = krcol + (m22 - 1) * 3;
                /* Computing MAX */
                i__4 = k + 1;
                i__5 = jbot;
                for (j = max(i__4,*ktop);
                        j <= i__5;
                        ++j)
                {
                    r_cnjg(&q__2, &v[m22 * v_dim1 + 1]);
                    i__4 = k + 1 + j * h_dim1;
                    r_cnjg(&q__5, &v[m22 * v_dim1 + 2]);
                    i__6 = k + 2 + j * h_dim1;
                    q__4.r = q__5.r * h__[i__6].r - q__5.i * h__[i__6].i;
                    q__4.i = q__5.r * h__[i__6].i + q__5.i * h__[i__6] .r; // , expr subst
                    q__3.r = h__[i__4].r + q__4.r;
                    q__3.i = h__[i__4].i + q__4.i; // , expr subst
                    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                    q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                    refsum.r = q__1.r;
                    refsum.i = q__1.i; // , expr subst
                    i__4 = k + 1 + j * h_dim1;
                    i__6 = k + 1 + j * h_dim1;
                    q__1.r = h__[i__6].r - refsum.r;
                    q__1.i = h__[i__6].i - refsum.i; // , expr subst
                    h__[i__4].r = q__1.r;
                    h__[i__4].i = q__1.i; // , expr subst
                    i__4 = k + 2 + j * h_dim1;
                    i__6 = k + 2 + j * h_dim1;
                    i__7 = m22 * v_dim1 + 2;
                    q__2.r = refsum.r * v[i__7].r - refsum.i * v[i__7].i;
                    q__2.i = refsum.r * v[i__7].i + refsum.i * v[i__7] .r; // , expr subst
                    q__1.r = h__[i__6].r - q__2.r;
                    q__1.i = h__[i__6].i - q__2.i; // , expr subst
                    h__[i__4].r = q__1.r;
                    h__[i__4].i = q__1.i; // , expr subst
                    /* L40: */
                }
            }
            /* ==== Multiply H by reflections from the right. */
            /* . Delay filling in the last row until the */
            /* . vigilant deflation check is complete. ==== */
            if (accum)
            {
                jtop = max(*ktop,incol);
            }
            else if (*wantt)
            {
                jtop = 1;
            }
            else
            {
                jtop = *ktop;
            }
            i__5 = mbot;
            for (m = mtop;
                    m <= i__5;
                    ++m)
            {
                i__4 = m * v_dim1 + 1;
                if (v[i__4].r != 0.f || v[i__4].i != 0.f)
                {
                    k = krcol + (m - 1) * 3;
                    /* Computing MIN */
                    i__6 = *kbot;
                    i__7 = k + 3; // , expr subst
                    i__4 = min(i__6,i__7);
                    for (j = jtop;
                            j <= i__4;
                            ++j)
                    {
                        i__6 = m * v_dim1 + 1;
                        i__7 = j + (k + 1) * h_dim1;
                        i__8 = m * v_dim1 + 2;
                        i__9 = j + (k + 2) * h_dim1;
                        q__4.r = v[i__8].r * h__[i__9].r - v[i__8].i * h__[ i__9].i;
                        q__4.i = v[i__8].r * h__[i__9].i + v[ i__8].i * h__[i__9].r; // , expr subst
                        q__3.r = h__[i__7].r + q__4.r;
                        q__3.i = h__[i__7].i + q__4.i; // , expr subst
                        i__10 = m * v_dim1 + 3;
                        i__11 = j + (k + 3) * h_dim1;
                        q__5.r = v[i__10].r * h__[i__11].r - v[i__10].i * h__[ i__11].i;
                        q__5.i = v[i__10].r * h__[i__11].i + v[i__10].i * h__[i__11].r; // , expr subst
                        q__2.r = q__3.r + q__5.r;
                        q__2.i = q__3.i + q__5.i; // , expr subst
                        q__1.r = v[i__6].r * q__2.r - v[i__6].i * q__2.i;
                        q__1.i = v[i__6].r * q__2.i + v[i__6].i * q__2.r; // , expr subst
                        refsum.r = q__1.r;
                        refsum.i = q__1.i; // , expr subst
                        i__6 = j + (k + 1) * h_dim1;
                        i__7 = j + (k + 1) * h_dim1;
                        q__1.r = h__[i__7].r - refsum.r;
                        q__1.i = h__[i__7].i - refsum.i; // , expr subst
                        h__[i__6].r = q__1.r;
                        h__[i__6].i = q__1.i; // , expr subst
                        i__6 = j + (k + 2) * h_dim1;
                        i__7 = j + (k + 2) * h_dim1;
                        r_cnjg(&q__3, &v[m * v_dim1 + 2]);
                        q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                        q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                        q__1.r = h__[i__7].r - q__2.r;
                        q__1.i = h__[i__7].i - q__2.i; // , expr subst
                        h__[i__6].r = q__1.r;
                        h__[i__6].i = q__1.i; // , expr subst
                        i__6 = j + (k + 3) * h_dim1;
                        i__7 = j + (k + 3) * h_dim1;
                        r_cnjg(&q__3, &v[m * v_dim1 + 3]);
                        q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                        q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                        q__1.r = h__[i__7].r - q__2.r;
                        q__1.i = h__[i__7].i - q__2.i; // , expr subst
                        h__[i__6].r = q__1.r;
                        h__[i__6].i = q__1.i; // , expr subst
                        /* L50: */
                    }
                    if (accum)
                    {
                        /* ==== Accumulate U. (If necessary, update Z later */
                        /* . with with an efficient matrix-matrix */
                        /* . multiply.) ==== */
                        kms = k - incol;
                        /* Computing MAX */
                        i__4 = 1;
                        i__6 = *ktop - incol; // , expr subst
                        i__7 = kdu;
                        for (j = max(i__4,i__6);
                                j <= i__7;
                                ++j)
                        {
                            i__4 = m * v_dim1 + 1;
                            i__6 = j + (kms + 1) * u_dim1;
                            i__8 = m * v_dim1 + 2;
                            i__9 = j + (kms + 2) * u_dim1;
                            q__4.r = v[i__8].r * u[i__9].r - v[i__8].i * u[ i__9].i;
                            q__4.i = v[i__8].r * u[i__9].i + v[i__8].i * u[i__9].r; // , expr subst
                            q__3.r = u[i__6].r + q__4.r;
                            q__3.i = u[i__6].i + q__4.i; // , expr subst
                            i__10 = m * v_dim1 + 3;
                            i__11 = j + (kms + 3) * u_dim1;
                            q__5.r = v[i__10].r * u[i__11].r - v[i__10].i * u[ i__11].i;
                            q__5.i = v[i__10].r * u[i__11] .i + v[i__10].i * u[i__11].r; // , expr subst
                            q__2.r = q__3.r + q__5.r;
                            q__2.i = q__3.i + q__5.i; // , expr subst
                            q__1.r = v[i__4].r * q__2.r - v[i__4].i * q__2.i;
                            q__1.i = v[i__4].r * q__2.i + v[i__4].i * q__2.r; // , expr subst
                            refsum.r = q__1.r;
                            refsum.i = q__1.i; // , expr subst
                            i__4 = j + (kms + 1) * u_dim1;
                            i__6 = j + (kms + 1) * u_dim1;
                            q__1.r = u[i__6].r - refsum.r;
                            q__1.i = u[i__6].i - refsum.i; // , expr subst
                            u[i__4].r = q__1.r;
                            u[i__4].i = q__1.i; // , expr subst
                            i__4 = j + (kms + 2) * u_dim1;
                            i__6 = j + (kms + 2) * u_dim1;
                            r_cnjg(&q__3, &v[m * v_dim1 + 2]);
                            q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                            q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                            q__1.r = u[i__6].r - q__2.r;
                            q__1.i = u[i__6].i - q__2.i; // , expr subst
                            u[i__4].r = q__1.r;
                            u[i__4].i = q__1.i; // , expr subst
                            i__4 = j + (kms + 3) * u_dim1;
                            i__6 = j + (kms + 3) * u_dim1;
                            r_cnjg(&q__3, &v[m * v_dim1 + 3]);
                            q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                            q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                            q__1.r = u[i__6].r - q__2.r;
                            q__1.i = u[i__6].i - q__2.i; // , expr subst
                            u[i__4].r = q__1.r;
                            u[i__4].i = q__1.i; // , expr subst
                            /* L60: */
                        }
                    }
                    else if (*wantz)
                    {
                        /* ==== U is not accumulated, so update Z */
                        /* . now by multiplying by reflections */
                        /* . from the right. ==== */
                        i__7 = *ihiz;
                        for (j = *iloz;
                                j <= i__7;
                                ++j)
                        {
                            i__4 = m * v_dim1 + 1;
                            i__6 = j + (k + 1) * z_dim1;
                            i__8 = m * v_dim1 + 2;
                            i__9 = j + (k + 2) * z_dim1;
                            q__4.r = v[i__8].r * z__[i__9].r - v[i__8].i * z__[i__9].i;
                            q__4.i = v[i__8].r * z__[ i__9].i + v[i__8].i * z__[i__9].r; // , expr subst
                            q__3.r = z__[i__6].r + q__4.r;
                            q__3.i = z__[i__6] .i + q__4.i; // , expr subst
                            i__10 = m * v_dim1 + 3;
                            i__11 = j + (k + 3) * z_dim1;
                            q__5.r = v[i__10].r * z__[i__11].r - v[i__10].i * z__[i__11].i;
                            q__5.i = v[i__10].r * z__[ i__11].i + v[i__10].i * z__[i__11].r; // , expr subst
                            q__2.r = q__3.r + q__5.r;
                            q__2.i = q__3.i + q__5.i; // , expr subst
                            q__1.r = v[i__4].r * q__2.r - v[i__4].i * q__2.i;
                            q__1.i = v[i__4].r * q__2.i + v[i__4].i * q__2.r; // , expr subst
                            refsum.r = q__1.r;
                            refsum.i = q__1.i; // , expr subst
                            i__4 = j + (k + 1) * z_dim1;
                            i__6 = j + (k + 1) * z_dim1;
                            q__1.r = z__[i__6].r - refsum.r;
                            q__1.i = z__[ i__6].i - refsum.i; // , expr subst
                            z__[i__4].r = q__1.r;
                            z__[i__4].i = q__1.i; // , expr subst
                            i__4 = j + (k + 2) * z_dim1;
                            i__6 = j + (k + 2) * z_dim1;
                            r_cnjg(&q__3, &v[m * v_dim1 + 2]);
                            q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                            q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                            q__1.r = z__[i__6].r - q__2.r;
                            q__1.i = z__[i__6] .i - q__2.i; // , expr subst
                            z__[i__4].r = q__1.r;
                            z__[i__4].i = q__1.i; // , expr subst
                            i__4 = j + (k + 3) * z_dim1;
                            i__6 = j + (k + 3) * z_dim1;
                            r_cnjg(&q__3, &v[m * v_dim1 + 3]);
                            q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                            q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                            q__1.r = z__[i__6].r - q__2.r;
                            q__1.i = z__[i__6] .i - q__2.i; // , expr subst
                            z__[i__4].r = q__1.r;
                            z__[i__4].i = q__1.i; // , expr subst
                            /* L70: */
                        }
                    }
                }
                /* L80: */
            }
            /* ==== Special case: 2-by-2 reflection (if needed) ==== */
            k = krcol + (m22 - 1) * 3;
            if (bmp22)
            {
                i__5 = m22 * v_dim1 + 1;
                if (v[i__5].r != 0.f || v[i__5].i != 0.f)
                {
                    /* Computing MIN */
                    i__7 = *kbot;
                    i__4 = k + 3; // , expr subst
                    i__5 = min(i__7,i__4);
                    for (j = jtop;
                            j <= i__5;
                            ++j)
                    {
                        i__7 = m22 * v_dim1 + 1;
                        i__4 = j + (k + 1) * h_dim1;
                        i__6 = m22 * v_dim1 + 2;
                        i__8 = j + (k + 2) * h_dim1;
                        q__3.r = v[i__6].r * h__[i__8].r - v[i__6].i * h__[ i__8].i;
                        q__3.i = v[i__6].r * h__[i__8].i + v[ i__6].i * h__[i__8].r; // , expr subst
                        q__2.r = h__[i__4].r + q__3.r;
                        q__2.i = h__[i__4].i + q__3.i; // , expr subst
                        q__1.r = v[i__7].r * q__2.r - v[i__7].i * q__2.i;
                        q__1.i = v[i__7].r * q__2.i + v[i__7].i * q__2.r; // , expr subst
                        refsum.r = q__1.r;
                        refsum.i = q__1.i; // , expr subst
                        i__7 = j + (k + 1) * h_dim1;
                        i__4 = j + (k + 1) * h_dim1;
                        q__1.r = h__[i__4].r - refsum.r;
                        q__1.i = h__[i__4].i - refsum.i; // , expr subst
                        h__[i__7].r = q__1.r;
                        h__[i__7].i = q__1.i; // , expr subst
                        i__7 = j + (k + 2) * h_dim1;
                        i__4 = j + (k + 2) * h_dim1;
                        r_cnjg(&q__3, &v[m22 * v_dim1 + 2]);
                        q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                        q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                        q__1.r = h__[i__4].r - q__2.r;
                        q__1.i = h__[i__4].i - q__2.i; // , expr subst
                        h__[i__7].r = q__1.r;
                        h__[i__7].i = q__1.i; // , expr subst
                        /* L90: */
                    }
                    if (accum)
                    {
                        kms = k - incol;
                        /* Computing MAX */
                        i__5 = 1;
                        i__7 = *ktop - incol; // , expr subst
                        i__4 = kdu;
                        for (j = max(i__5,i__7);
                                j <= i__4;
                                ++j)
                        {
                            i__5 = m22 * v_dim1 + 1;
                            i__7 = j + (kms + 1) * u_dim1;
                            i__6 = m22 * v_dim1 + 2;
                            i__8 = j + (kms + 2) * u_dim1;
                            q__3.r = v[i__6].r * u[i__8].r - v[i__6].i * u[ i__8].i;
                            q__3.i = v[i__6].r * u[i__8].i + v[i__6].i * u[i__8].r; // , expr subst
                            q__2.r = u[i__7].r + q__3.r;
                            q__2.i = u[i__7].i + q__3.i; // , expr subst
                            q__1.r = v[i__5].r * q__2.r - v[i__5].i * q__2.i;
                            q__1.i = v[i__5].r * q__2.i + v[i__5].i * q__2.r; // , expr subst
                            refsum.r = q__1.r;
                            refsum.i = q__1.i; // , expr subst
                            i__5 = j + (kms + 1) * u_dim1;
                            i__7 = j + (kms + 1) * u_dim1;
                            q__1.r = u[i__7].r - refsum.r;
                            q__1.i = u[i__7].i - refsum.i; // , expr subst
                            u[i__5].r = q__1.r;
                            u[i__5].i = q__1.i; // , expr subst
                            i__5 = j + (kms + 2) * u_dim1;
                            i__7 = j + (kms + 2) * u_dim1;
                            r_cnjg(&q__3, &v[m22 * v_dim1 + 2]);
                            q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                            q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                            q__1.r = u[i__7].r - q__2.r;
                            q__1.i = u[i__7].i - q__2.i; // , expr subst
                            u[i__5].r = q__1.r;
                            u[i__5].i = q__1.i; // , expr subst
                            /* L100: */
                        }
                    }
                    else if (*wantz)
                    {
                        i__4 = *ihiz;
                        for (j = *iloz;
                                j <= i__4;
                                ++j)
                        {
                            i__5 = m22 * v_dim1 + 1;
                            i__7 = j + (k + 1) * z_dim1;
                            i__6 = m22 * v_dim1 + 2;
                            i__8 = j + (k + 2) * z_dim1;
                            q__3.r = v[i__6].r * z__[i__8].r - v[i__6].i * z__[i__8].i;
                            q__3.i = v[i__6].r * z__[ i__8].i + v[i__6].i * z__[i__8].r; // , expr subst
                            q__2.r = z__[i__7].r + q__3.r;
                            q__2.i = z__[i__7] .i + q__3.i; // , expr subst
                            q__1.r = v[i__5].r * q__2.r - v[i__5].i * q__2.i;
                            q__1.i = v[i__5].r * q__2.i + v[i__5].i * q__2.r; // , expr subst
                            refsum.r = q__1.r;
                            refsum.i = q__1.i; // , expr subst
                            i__5 = j + (k + 1) * z_dim1;
                            i__7 = j + (k + 1) * z_dim1;
                            q__1.r = z__[i__7].r - refsum.r;
                            q__1.i = z__[ i__7].i - refsum.i; // , expr subst
                            z__[i__5].r = q__1.r;
                            z__[i__5].i = q__1.i; // , expr subst
                            i__5 = j + (k + 2) * z_dim1;
                            i__7 = j + (k + 2) * z_dim1;
                            r_cnjg(&q__3, &v[m22 * v_dim1 + 2]);
                            q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                            q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                            q__1.r = z__[i__7].r - q__2.r;
                            q__1.i = z__[i__7] .i - q__2.i; // , expr subst
                            z__[i__5].r = q__1.r;
                            z__[i__5].i = q__1.i; // , expr subst
                            /* L110: */
                        }
                    }
                }
            }
            /* ==== Vigilant deflation check ==== */
            mstart = mtop;
            if (krcol + (mstart - 1) * 3 < *ktop)
            {
                ++mstart;
            }
            mend = mbot;
            if (bmp22)
            {
                ++mend;
            }
            if (krcol == *kbot - 2)
            {
                ++mend;
            }
            i__4 = mend;
            for (m = mstart;
                    m <= i__4;
                    ++m)
            {
                /* Computing MIN */
                i__5 = *kbot - 1;
                i__7 = krcol + (m - 1) * 3; // , expr subst
                k = min(i__5,i__7);
                /* ==== The following convergence test requires that */
                /* . the tradition small-compared-to-nearby-diagonals */
                /* . criterion and the Ahues & Tisseur (LAWN 122, 1997) */
                /* . criteria both be satisfied. The latter improves */
                /* . accuracy in some examples. Falling back on an */
                /* . alternate convergence criterion when TST1 or TST2 */
                /* . is zero (as done here) is traditional but probably */
                /* . unnecessary. ==== */
                i__5 = k + 1 + k * h_dim1;
                if (h__[i__5].r != 0.f || h__[i__5].i != 0.f)
                {
                    i__5 = k + k * h_dim1;
                    i__7 = k + 1 + (k + 1) * h_dim1;
                    tst1 = (r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(& h__[k + k * h_dim1]), f2c_abs(r__2)) + ((r__3 = h__[ i__7].r, f2c_abs(r__3)) + (r__4 = r_imag(&h__[k + 1 + (k + 1) * h_dim1]), f2c_abs(r__4)));
                    if (tst1 == 0.f)
                    {
                        if (k >= *ktop + 1)
                        {
                            i__5 = k + (k - 1) * h_dim1;
                            tst1 += (r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + (k - 1) * h_dim1]), f2c_abs( r__2));
                        }
                        if (k >= *ktop + 2)
                        {
                            i__5 = k + (k - 2) * h_dim1;
                            tst1 += (r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + (k - 2) * h_dim1]), f2c_abs( r__2));
                        }
                        if (k >= *ktop + 3)
                        {
                            i__5 = k + (k - 3) * h_dim1;
                            tst1 += (r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + (k - 3) * h_dim1]), f2c_abs( r__2));
                        }
                        if (k <= *kbot - 2)
                        {
                            i__5 = k + 2 + (k + 1) * h_dim1;
                            tst1 += (r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + 2 + (k + 1) * h_dim1]), f2c_abs(r__2));
                        }
                        if (k <= *kbot - 3)
                        {
                            i__5 = k + 3 + (k + 1) * h_dim1;
                            tst1 += (r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + 3 + (k + 1) * h_dim1]), f2c_abs(r__2));
                        }
                        if (k <= *kbot - 4)
                        {
                            i__5 = k + 4 + (k + 1) * h_dim1;
                            tst1 += (r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + 4 + (k + 1) * h_dim1]), f2c_abs(r__2));
                        }
                    }
                    i__5 = k + 1 + k * h_dim1;
                    /* Computing MAX */
                    r__3 = smlnum;
                    r__4 = ulp * tst1; // , expr subst
                    if ((r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[ k + 1 + k * h_dim1]), f2c_abs(r__2)) <= max(r__3,r__4) )
                    {
                        /* Computing MAX */
                        i__5 = k + 1 + k * h_dim1;
                        i__7 = k + (k + 1) * h_dim1;
                        r__5 = (r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + 1 + k * h_dim1]), f2c_abs(r__2));
                        r__6 = (r__3 = h__[i__7].r, f2c_abs(r__3)) + ( r__4 = r_imag(&h__[k + (k + 1) * h_dim1]), f2c_abs(r__4)); // , expr subst
                        h12 = max(r__5,r__6);
                        /* Computing MIN */
                        i__5 = k + 1 + k * h_dim1;
                        i__7 = k + (k + 1) * h_dim1;
                        r__5 = (r__1 = h__[i__5].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + 1 + k * h_dim1]), f2c_abs(r__2));
                        r__6 = (r__3 = h__[i__7].r, f2c_abs(r__3)) + ( r__4 = r_imag(&h__[k + (k + 1) * h_dim1]), f2c_abs(r__4)); // , expr subst
                        h21 = min(r__5,r__6);
                        i__5 = k + k * h_dim1;
                        i__7 = k + 1 + (k + 1) * h_dim1;
                        q__2.r = h__[i__5].r - h__[i__7].r;
                        q__2.i = h__[i__5] .i - h__[i__7].i; // , expr subst
                        q__1.r = q__2.r;
                        q__1.i = q__2.i; // , expr subst
                        /* Computing MAX */
                        i__6 = k + 1 + (k + 1) * h_dim1;
                        r__5 = (r__1 = h__[i__6].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + 1 + (k + 1) * h_dim1]), f2c_abs( r__2));
                        r__6 = (r__3 = q__1.r, f2c_abs(r__3)) + ( r__4 = r_imag(&q__1), f2c_abs(r__4)); // , expr subst
                        h11 = max(r__5,r__6);
                        i__5 = k + k * h_dim1;
                        i__7 = k + 1 + (k + 1) * h_dim1;
                        q__2.r = h__[i__5].r - h__[i__7].r;
                        q__2.i = h__[i__5] .i - h__[i__7].i; // , expr subst
                        q__1.r = q__2.r;
                        q__1.i = q__2.i; // , expr subst
                        /* Computing MIN */
                        i__6 = k + 1 + (k + 1) * h_dim1;
                        r__5 = (r__1 = h__[i__6].r, f2c_abs(r__1)) + (r__2 = r_imag(&h__[k + 1 + (k + 1) * h_dim1]), f2c_abs( r__2));
                        r__6 = (r__3 = q__1.r, f2c_abs(r__3)) + ( r__4 = r_imag(&q__1), f2c_abs(r__4)); // , expr subst
                        h22 = min(r__5,r__6);
                        scl = h11 + h12;
                        tst2 = h22 * (h11 / scl);
                        /* Computing MAX */
                        r__1 = smlnum;
                        r__2 = ulp * tst2; // , expr subst
                        if (tst2 == 0.f || h21 * (h12 / scl) <= max(r__1,r__2) )
                        {
                            i__5 = k + 1 + k * h_dim1;
                            h__[i__5].r = 0.f;
                            h__[i__5].i = 0.f; // , expr subst
                        }
                    }
                }
                /* L120: */
            }
            /* ==== Fill in the last row of each bulge. ==== */
            /* Computing MIN */
            i__4 = nbmps;
            i__5 = (*kbot - krcol - 1) / 3; // , expr subst
            mend = min(i__4,i__5);
            i__4 = mend;
            for (m = mtop;
                    m <= i__4;
                    ++m)
            {
                k = krcol + (m - 1) * 3;
                i__5 = m * v_dim1 + 1;
                i__7 = m * v_dim1 + 3;
                q__2.r = v[i__5].r * v[i__7].r - v[i__5].i * v[i__7].i;
                q__2.i = v[i__5].r * v[i__7].i + v[i__5].i * v[i__7] .r; // , expr subst
                i__6 = k + 4 + (k + 3) * h_dim1;
                q__1.r = q__2.r * h__[i__6].r - q__2.i * h__[i__6].i;
                q__1.i = q__2.r * h__[i__6].i + q__2.i * h__[i__6].r; // , expr subst
                refsum.r = q__1.r;
                refsum.i = q__1.i; // , expr subst
                i__5 = k + 4 + (k + 1) * h_dim1;
                q__1.r = -refsum.r;
                q__1.i = -refsum.i; // , expr subst
                h__[i__5].r = q__1.r;
                h__[i__5].i = q__1.i; // , expr subst
                i__5 = k + 4 + (k + 2) * h_dim1;
                q__2.r = -refsum.r;
                q__2.i = -refsum.i; // , expr subst
                r_cnjg(&q__3, &v[m * v_dim1 + 2]);
                q__1.r = q__2.r * q__3.r - q__2.i * q__3.i;
                q__1.i = q__2.r * q__3.i + q__2.i * q__3.r; // , expr subst
                h__[i__5].r = q__1.r;
                h__[i__5].i = q__1.i; // , expr subst
                i__5 = k + 4 + (k + 3) * h_dim1;
                i__7 = k + 4 + (k + 3) * h_dim1;
                r_cnjg(&q__3, &v[m * v_dim1 + 3]);
                q__2.r = refsum.r * q__3.r - refsum.i * q__3.i;
                q__2.i = refsum.r * q__3.i + refsum.i * q__3.r; // , expr subst
                q__1.r = h__[i__7].r - q__2.r;
                q__1.i = h__[i__7].i - q__2.i; // , expr subst
                h__[i__5].r = q__1.r;
                h__[i__5].i = q__1.i; // , expr subst
                /* L130: */
            }
            /* ==== End of near-the-diagonal bulge chase. ==== */
            /* L140: */
        }
        /* ==== Use U (if accumulated) to update far-from-diagonal */
        /* . entries in H. If required, use U to update Z as */
        /* . well. ==== */
        if (accum)
        {
            if (*wantt)
            {
                jtop = 1;
                jbot = *n;
            }
            else
            {
                jtop = *ktop;
                jbot = *kbot;
            }
            if (! blk22 || incol < *ktop || ndcol > *kbot || ns <= 2)
            {
                /* ==== Updates not exploiting the 2-by-2 block */
                /* . structure of U. K1 and NU keep track of */
                /* . the location and size of U in the special */
                /* . cases of introducing bulges and chasing */
                /* . bulges off the bottom. In these special */
                /* . cases and in case the number of shifts */
                /* . is NS = 2, there is no 2-by-2 block */
                /* . structure to exploit. ==== */
                /* Computing MAX */
                i__3 = 1;
                i__4 = *ktop - incol; // , expr subst
                k1 = max(i__3,i__4);
                /* Computing MAX */
                i__3 = 0;
                i__4 = ndcol - *kbot; // , expr subst
                nu = kdu - max(i__3,i__4) - k1 + 1;
                /* ==== Horizontal Multiply ==== */
                i__3 = jbot;
                i__4 = *nh;
                for (jcol = min(ndcol,*kbot) + 1;
                        i__4 < 0 ? jcol >= i__3 : jcol <= i__3;
                        jcol += i__4)
                {
                    /* Computing MIN */
                    i__5 = *nh;
                    i__7 = jbot - jcol + 1; // , expr subst
                    jlen = min(i__5,i__7);
                    cgemm_("C", "N", &nu, &jlen, &nu, &c_b2, &u[k1 + k1 * u_dim1], ldu, &h__[incol + k1 + jcol * h_dim1], ldh, &c_b1, &wh[wh_offset], ldwh);
                    clacpy_("ALL", &nu, &jlen, &wh[wh_offset], ldwh, &h__[ incol + k1 + jcol * h_dim1], ldh);
                    /* L150: */
                }
                /* ==== Vertical multiply ==== */
                i__4 = max(*ktop,incol) - 1;
                i__3 = *nv;
                for (jrow = jtop;
                        i__3 < 0 ? jrow >= i__4 : jrow <= i__4;
                        jrow += i__3)
                {
                    /* Computing MIN */
                    i__5 = *nv;
                    i__7 = max(*ktop,incol) - jrow; // , expr subst
                    jlen = min(i__5,i__7);
                    cgemm_("N", "N", &jlen, &nu, &nu, &c_b2, &h__[jrow + ( incol + k1) * h_dim1], ldh, &u[k1 + k1 * u_dim1], ldu, &c_b1, &wv[wv_offset], ldwv);
                    clacpy_("ALL", &jlen, &nu, &wv[wv_offset], ldwv, &h__[ jrow + (incol + k1) * h_dim1], ldh);
                    /* L160: */
                }
                /* ==== Z multiply (also vertical) ==== */
                if (*wantz)
                {
                    i__3 = *ihiz;
                    i__4 = *nv;
                    for (jrow = *iloz;
                            i__4 < 0 ? jrow >= i__3 : jrow <= i__3;
                            jrow += i__4)
                    {
                        /* Computing MIN */
                        i__5 = *nv;
                        i__7 = *ihiz - jrow + 1; // , expr subst
                        jlen = min(i__5,i__7);
                        cgemm_("N", "N", &jlen, &nu, &nu, &c_b2, &z__[jrow + ( incol + k1) * z_dim1], ldz, &u[k1 + k1 * u_dim1], ldu, &c_b1, &wv[wv_offset], ldwv);
                        clacpy_("ALL", &jlen, &nu, &wv[wv_offset], ldwv, &z__[ jrow + (incol + k1) * z_dim1], ldz) ;
                        /* L170: */
                    }
                }
            }
            else
            {
                /* ==== Updates exploiting U's 2-by-2 block structure. */
                /* . (I2, I4, J2, J4 are the last rows and columns */
                /* . of the blocks.) ==== */
                i2 = (kdu + 1) / 2;
                i4 = kdu;
                j2 = i4 - i2;
                j4 = kdu;
                /* ==== KZS and KNZ deal with the band of zeros */
                /* . along the diagonal of one of the triangular */
                /* . blocks. ==== */
                kzs = j4 - j2 - (ns + 1);
                knz = ns + 1;
                /* ==== Horizontal multiply ==== */
                i__4 = jbot;
                i__3 = *nh;
                for (jcol = min(ndcol,*kbot) + 1;
                        i__3 < 0 ? jcol >= i__4 : jcol <= i__4;
                        jcol += i__3)
                {
                    /* Computing MIN */
                    i__5 = *nh;
                    i__7 = jbot - jcol + 1; // , expr subst
                    jlen = min(i__5,i__7);
                    /* ==== Copy bottom of H to top+KZS of scratch ==== */
                    /* (The first KZS rows get multiplied by zero.) ==== */
                    clacpy_("ALL", &knz, &jlen, &h__[incol + 1 + j2 + jcol * h_dim1], ldh, &wh[kzs + 1 + wh_dim1], ldwh);
                    /* ==== Multiply by U21**H ==== */
                    claset_("ALL", &kzs, &jlen, &c_b1, &c_b1, &wh[wh_offset], ldwh);
                    ctrmm_("L", "U", "C", "N", &knz, &jlen, &c_b2, &u[j2 + 1 + (kzs + 1) * u_dim1], ldu, &wh[kzs + 1 + wh_dim1] , ldwh);
                    /* ==== Multiply top of H by U11**H ==== */
                    cgemm_("C", "N", &i2, &jlen, &j2, &c_b2, &u[u_offset], ldu, &h__[incol + 1 + jcol * h_dim1], ldh, &c_b2, &wh[wh_offset], ldwh);
                    /* ==== Copy top of H to bottom of WH ==== */
                    clacpy_("ALL", &j2, &jlen, &h__[incol + 1 + jcol * h_dim1] , ldh, &wh[i2 + 1 + wh_dim1], ldwh);
                    /* ==== Multiply by U21**H ==== */
                    ctrmm_("L", "L", "C", "N", &j2, &jlen, &c_b2, &u[(i2 + 1) * u_dim1 + 1], ldu, &wh[i2 + 1 + wh_dim1], ldwh);
                    /* ==== Multiply by U22 ==== */
                    i__5 = i4 - i2;
                    i__7 = j4 - j2;
                    cgemm_("C", "N", &i__5, &jlen, &i__7, &c_b2, &u[j2 + 1 + ( i2 + 1) * u_dim1], ldu, &h__[incol + 1 + j2 + jcol * h_dim1], ldh, &c_b2, &wh[i2 + 1 + wh_dim1], ldwh);
                    /* ==== Copy it back ==== */
                    clacpy_("ALL", &kdu, &jlen, &wh[wh_offset], ldwh, &h__[ incol + 1 + jcol * h_dim1], ldh);
                    /* L180: */
                }
                /* ==== Vertical multiply ==== */
                i__3 = max(incol,*ktop) - 1;
                i__4 = *nv;
                for (jrow = jtop;
                        i__4 < 0 ? jrow >= i__3 : jrow <= i__3;
                        jrow += i__4)
                {
                    /* Computing MIN */
                    i__5 = *nv;
                    i__7 = max(incol,*ktop) - jrow; // , expr subst
                    jlen = min(i__5,i__7);
                    /* ==== Copy right of H to scratch (the first KZS */
                    /* . columns get multiplied by zero) ==== */
                    clacpy_("ALL", &jlen, &knz, &h__[jrow + (incol + 1 + j2) * h_dim1], ldh, &wv[(kzs + 1) * wv_dim1 + 1], ldwv);
                    /* ==== Multiply by U21 ==== */
                    claset_("ALL", &jlen, &kzs, &c_b1, &c_b1, &wv[wv_offset], ldwv);
                    ctrmm_("R", "U", "N", "N", &jlen, &knz, &c_b2, &u[j2 + 1 + (kzs + 1) * u_dim1], ldu, &wv[(kzs + 1) * wv_dim1 + 1], ldwv);
                    /* ==== Multiply by U11 ==== */
                    cgemm_("N", "N", &jlen, &i2, &j2, &c_b2, &h__[jrow + ( incol + 1) * h_dim1], ldh, &u[u_offset], ldu, & c_b2, &wv[wv_offset], ldwv);
                    /* ==== Copy left of H to right of scratch ==== */
                    clacpy_("ALL", &jlen, &j2, &h__[jrow + (incol + 1) * h_dim1], ldh, &wv[(i2 + 1) * wv_dim1 + 1], ldwv);
                    /* ==== Multiply by U21 ==== */
                    i__5 = i4 - i2;
                    ctrmm_("R", "L", "N", "N", &jlen, &i__5, &c_b2, &u[(i2 + 1) * u_dim1 + 1], ldu, &wv[(i2 + 1) * wv_dim1 + 1] , ldwv);
                    /* ==== Multiply by U22 ==== */
                    i__5 = i4 - i2;
                    i__7 = j4 - j2;
                    cgemm_("N", "N", &jlen, &i__5, &i__7, &c_b2, &h__[jrow + ( incol + 1 + j2) * h_dim1], ldh, &u[j2 + 1 + (i2 + 1) * u_dim1], ldu, &c_b2, &wv[(i2 + 1) * wv_dim1 + 1], ldwv);
                    /* ==== Copy it back ==== */
                    clacpy_("ALL", &jlen, &kdu, &wv[wv_offset], ldwv, &h__[ jrow + (incol + 1) * h_dim1], ldh);
                    /* L190: */
                }
                /* ==== Multiply Z (also vertical) ==== */
                if (*wantz)
                {
                    i__4 = *ihiz;
                    i__3 = *nv;
                    for (jrow = *iloz;
                            i__3 < 0 ? jrow >= i__4 : jrow <= i__4;
                            jrow += i__3)
                    {
                        /* Computing MIN */
                        i__5 = *nv;
                        i__7 = *ihiz - jrow + 1; // , expr subst
                        jlen = min(i__5,i__7);
                        /* ==== Copy right of Z to left of scratch (first */
                        /* . KZS columns get multiplied by zero) ==== */
                        clacpy_("ALL", &jlen, &knz, &z__[jrow + (incol + 1 + j2) * z_dim1], ldz, &wv[(kzs + 1) * wv_dim1 + 1], ldwv);
                        /* ==== Multiply by U12 ==== */
                        claset_("ALL", &jlen, &kzs, &c_b1, &c_b1, &wv[ wv_offset], ldwv);
                        ctrmm_("R", "U", "N", "N", &jlen, &knz, &c_b2, &u[j2 + 1 + (kzs + 1) * u_dim1], ldu, &wv[(kzs + 1) * wv_dim1 + 1], ldwv);
                        /* ==== Multiply by U11 ==== */
                        cgemm_("N", "N", &jlen, &i2, &j2, &c_b2, &z__[jrow + ( incol + 1) * z_dim1], ldz, &u[u_offset], ldu, &c_b2, &wv[wv_offset], ldwv);
                        /* ==== Copy left of Z to right of scratch ==== */
                        clacpy_("ALL", &jlen, &j2, &z__[jrow + (incol + 1) * z_dim1], ldz, &wv[(i2 + 1) * wv_dim1 + 1], ldwv);
                        /* ==== Multiply by U21 ==== */
                        i__5 = i4 - i2;
                        ctrmm_("R", "L", "N", "N", &jlen, &i__5, &c_b2, &u[( i2 + 1) * u_dim1 + 1], ldu, &wv[(i2 + 1) * wv_dim1 + 1], ldwv);
                        /* ==== Multiply by U22 ==== */
                        i__5 = i4 - i2;
                        i__7 = j4 - j2;
                        cgemm_("N", "N", &jlen, &i__5, &i__7, &c_b2, &z__[ jrow + (incol + 1 + j2) * z_dim1], ldz, &u[j2 + 1 + (i2 + 1) * u_dim1], ldu, &c_b2, &wv[(i2 + 1) * wv_dim1 + 1], ldwv);
                        /* ==== Copy the result back to Z ==== */
                        clacpy_("ALL", &jlen, &kdu, &wv[wv_offset], ldwv, & z__[jrow + (incol + 1) * z_dim1], ldz);
                        /* L200: */
                    }
                }
            }
        }
        /* L210: */
    }
    /* ==== End of CLAQR5 ==== */
    return 0;
}
/* claqr5_ */
