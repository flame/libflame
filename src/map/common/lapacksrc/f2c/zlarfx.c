/* ../netlib/zlarfx.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling whe n the reflector has order ≤ 10. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLARFX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfx. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfx. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfx. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLARFX( SIDE, M, N, V, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER LDC, M, N */
/* COMPLEX*16 TAU */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARFX applies a complex elementary reflector H to a complex m by n */
/* > matrix C, from either the left or the right. H is represented in the */
/* > form */
/* > */
/* > H = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar and v is a complex vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix */
/* > */
/* > This version uses inline code if H has order < 11. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': form H * C */
/* > = 'R': form C * H */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is COMPLEX*16 array, dimension (M) if SIDE = 'L' */
/* > or (N) if SIDE = 'R' */
/* > The vector v in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX*16 array, dimension (LDC,N) */
/* > On entry, the m by n matrix C. */
/* > On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* > or C * H if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N) if SIDE = 'L' */
/* > or (M) if SIDE = 'R' */
/* > WORK is not referenced if H has order < 11. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int zlarfx_(char *side, integer *m, integer *n, doublecomplex *v, doublecomplex *tau, doublecomplex *c__, integer * ldc, doublecomplex *work)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, i__11;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8, z__9, z__10, z__11, z__12, z__13, z__14, z__15, z__16, z__17, z__18, z__19;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    integer j;
    doublecomplex t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, v6, v7, v8, v9, t10, v10, sum;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int zlarf_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    if (tau->r == 0. && tau->i == 0.)
    {
        return 0;
    }
    if (lsame_(side, "L"))
    {
        /* Form H * C, where H has order m. */
        switch (*m)
        {
        case 1:
            goto L10;
        case 2:
            goto L30;
        case 3:
            goto L50;
        case 4:
            goto L70;
        case 5:
            goto L90;
        case 6:
            goto L110;
        case 7:
            goto L130;
        case 8:
            goto L150;
        case 9:
            goto L170;
        case 10:
            goto L190;
        }
        /* Code for general M */
        zlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1]);
        goto L410;
L10: /* Special code for 1 x 1 Householder */
        z__3.r = tau->r * v[1].r - tau->i * v[1].i;
        z__3.i = tau->r * v[1].i + tau->i * v[1].r; // , expr subst
        d_cnjg(&z__4, &v[1]);
        z__2.r = z__3.r * z__4.r - z__3.i * z__4.i;
        z__2.i = z__3.r * z__4.i + z__3.i * z__4.r; // , expr subst
        z__1.r = 1. - z__2.r;
        z__1.i = 0. - z__2.i; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__1.r = t1.r * c__[i__3].r - t1.i * c__[i__3].i;
            z__1.i = t1.r * c__[i__3].i + t1.i * c__[i__3].r; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L20: */
        }
        goto L410;
L30: /* Special code for 2 x 2 Householder */
        d_cnjg(&z__1, &v[1]);
        v1.r = z__1.r;
        v1.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[2]);
        v2.r = z__1.r;
        v2.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            z__2.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__2.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j * c_dim1 + 2;
            z__3.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__3.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__1.r = z__2.r + z__3.r;
            z__1.i = z__2.i + z__3.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L40: */
        }
        goto L410;
L50: /* Special code for 3 x 3 Householder */
        d_cnjg(&z__1, &v[1]);
        v1.r = z__1.r;
        v1.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[2]);
        v2.r = z__1.r;
        v2.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[3]);
        v3.r = z__1.r;
        v3.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            z__3.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__3.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j * c_dim1 + 2;
            z__4.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__4.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__2.r = z__3.r + z__4.r;
            z__2.i = z__3.i + z__4.i; // , expr subst
            i__4 = j * c_dim1 + 3;
            z__5.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__5.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__1.r = z__2.r + z__5.r;
            z__1.i = z__2.i + z__5.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L60: */
        }
        goto L410;
L70: /* Special code for 4 x 4 Householder */
        d_cnjg(&z__1, &v[1]);
        v1.r = z__1.r;
        v1.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[2]);
        v2.r = z__1.r;
        v2.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[3]);
        v3.r = z__1.r;
        v3.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[4]);
        v4.r = z__1.r;
        v4.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            z__4.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__4.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j * c_dim1 + 2;
            z__5.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__5.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__3.r = z__4.r + z__5.r;
            z__3.i = z__4.i + z__5.i; // , expr subst
            i__4 = j * c_dim1 + 3;
            z__6.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__6.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__2.r = z__3.r + z__6.r;
            z__2.i = z__3.i + z__6.i; // , expr subst
            i__5 = j * c_dim1 + 4;
            z__7.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__7.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__1.r = z__2.r + z__7.r;
            z__1.i = z__2.i + z__7.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L80: */
        }
        goto L410;
L90: /* Special code for 5 x 5 Householder */
        d_cnjg(&z__1, &v[1]);
        v1.r = z__1.r;
        v1.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[2]);
        v2.r = z__1.r;
        v2.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[3]);
        v3.r = z__1.r;
        v3.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[4]);
        v4.r = z__1.r;
        v4.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[5]);
        v5.r = z__1.r;
        v5.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            z__5.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__5.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j * c_dim1 + 2;
            z__6.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__6.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__4.r = z__5.r + z__6.r;
            z__4.i = z__5.i + z__6.i; // , expr subst
            i__4 = j * c_dim1 + 3;
            z__7.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__7.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__3.r = z__4.r + z__7.r;
            z__3.i = z__4.i + z__7.i; // , expr subst
            i__5 = j * c_dim1 + 4;
            z__8.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__8.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__2.r = z__3.r + z__8.r;
            z__2.i = z__3.i + z__8.i; // , expr subst
            i__6 = j * c_dim1 + 5;
            z__9.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__9.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__1.r = z__2.r + z__9.r;
            z__1.i = z__2.i + z__9.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L100: */
        }
        goto L410;
L110: /* Special code for 6 x 6 Householder */
        d_cnjg(&z__1, &v[1]);
        v1.r = z__1.r;
        v1.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[2]);
        v2.r = z__1.r;
        v2.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[3]);
        v3.r = z__1.r;
        v3.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[4]);
        v4.r = z__1.r;
        v4.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[5]);
        v5.r = z__1.r;
        v5.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[6]);
        v6.r = z__1.r;
        v6.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            z__6.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__6.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j * c_dim1 + 2;
            z__7.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__7.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__5.r = z__6.r + z__7.r;
            z__5.i = z__6.i + z__7.i; // , expr subst
            i__4 = j * c_dim1 + 3;
            z__8.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__8.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__4.r = z__5.r + z__8.r;
            z__4.i = z__5.i + z__8.i; // , expr subst
            i__5 = j * c_dim1 + 4;
            z__9.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__9.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__3.r = z__4.r + z__9.r;
            z__3.i = z__4.i + z__9.i; // , expr subst
            i__6 = j * c_dim1 + 5;
            z__10.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__10.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__2.r = z__3.r + z__10.r;
            z__2.i = z__3.i + z__10.i; // , expr subst
            i__7 = j * c_dim1 + 6;
            z__11.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__11.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__1.r = z__2.r + z__11.r;
            z__1.i = z__2.i + z__11.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L120: */
        }
        goto L410;
L130: /* Special code for 7 x 7 Householder */
        d_cnjg(&z__1, &v[1]);
        v1.r = z__1.r;
        v1.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[2]);
        v2.r = z__1.r;
        v2.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[3]);
        v3.r = z__1.r;
        v3.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[4]);
        v4.r = z__1.r;
        v4.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[5]);
        v5.r = z__1.r;
        v5.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[6]);
        v6.r = z__1.r;
        v6.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[7]);
        v7.r = z__1.r;
        v7.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v7);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t7.r = z__1.r;
        t7.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            z__7.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__7.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j * c_dim1 + 2;
            z__8.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__8.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__6.r = z__7.r + z__8.r;
            z__6.i = z__7.i + z__8.i; // , expr subst
            i__4 = j * c_dim1 + 3;
            z__9.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__9.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__5.r = z__6.r + z__9.r;
            z__5.i = z__6.i + z__9.i; // , expr subst
            i__5 = j * c_dim1 + 4;
            z__10.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__10.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__4.r = z__5.r + z__10.r;
            z__4.i = z__5.i + z__10.i; // , expr subst
            i__6 = j * c_dim1 + 5;
            z__11.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__11.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__3.r = z__4.r + z__11.r;
            z__3.i = z__4.i + z__11.i; // , expr subst
            i__7 = j * c_dim1 + 6;
            z__12.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__12.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__2.r = z__3.r + z__12.r;
            z__2.i = z__3.i + z__12.i; // , expr subst
            i__8 = j * c_dim1 + 7;
            z__13.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i;
            z__13.i = v7.r * c__[i__8].i + v7.i * c__[i__8].r; // , expr subst
            z__1.r = z__2.r + z__13.r;
            z__1.i = z__2.i + z__13.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 7;
            i__3 = j * c_dim1 + 7;
            z__2.r = sum.r * t7.r - sum.i * t7.i;
            z__2.i = sum.r * t7.i + sum.i * t7.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L140: */
        }
        goto L410;
L150: /* Special code for 8 x 8 Householder */
        d_cnjg(&z__1, &v[1]);
        v1.r = z__1.r;
        v1.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[2]);
        v2.r = z__1.r;
        v2.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[3]);
        v3.r = z__1.r;
        v3.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[4]);
        v4.r = z__1.r;
        v4.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[5]);
        v5.r = z__1.r;
        v5.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[6]);
        v6.r = z__1.r;
        v6.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[7]);
        v7.r = z__1.r;
        v7.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v7);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t7.r = z__1.r;
        t7.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[8]);
        v8.r = z__1.r;
        v8.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v8);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t8.r = z__1.r;
        t8.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            z__8.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__8.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j * c_dim1 + 2;
            z__9.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__9.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__7.r = z__8.r + z__9.r;
            z__7.i = z__8.i + z__9.i; // , expr subst
            i__4 = j * c_dim1 + 3;
            z__10.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__10.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__6.r = z__7.r + z__10.r;
            z__6.i = z__7.i + z__10.i; // , expr subst
            i__5 = j * c_dim1 + 4;
            z__11.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__11.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__5.r = z__6.r + z__11.r;
            z__5.i = z__6.i + z__11.i; // , expr subst
            i__6 = j * c_dim1 + 5;
            z__12.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__12.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__4.r = z__5.r + z__12.r;
            z__4.i = z__5.i + z__12.i; // , expr subst
            i__7 = j * c_dim1 + 6;
            z__13.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__13.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__3.r = z__4.r + z__13.r;
            z__3.i = z__4.i + z__13.i; // , expr subst
            i__8 = j * c_dim1 + 7;
            z__14.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i;
            z__14.i = v7.r * c__[i__8].i + v7.i * c__[i__8].r; // , expr subst
            z__2.r = z__3.r + z__14.r;
            z__2.i = z__3.i + z__14.i; // , expr subst
            i__9 = j * c_dim1 + 8;
            z__15.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i;
            z__15.i = v8.r * c__[i__9].i + v8.i * c__[i__9].r; // , expr subst
            z__1.r = z__2.r + z__15.r;
            z__1.i = z__2.i + z__15.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 7;
            i__3 = j * c_dim1 + 7;
            z__2.r = sum.r * t7.r - sum.i * t7.i;
            z__2.i = sum.r * t7.i + sum.i * t7.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 8;
            i__3 = j * c_dim1 + 8;
            z__2.r = sum.r * t8.r - sum.i * t8.i;
            z__2.i = sum.r * t8.i + sum.i * t8.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L160: */
        }
        goto L410;
L170: /* Special code for 9 x 9 Householder */
        d_cnjg(&z__1, &v[1]);
        v1.r = z__1.r;
        v1.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[2]);
        v2.r = z__1.r;
        v2.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[3]);
        v3.r = z__1.r;
        v3.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[4]);
        v4.r = z__1.r;
        v4.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[5]);
        v5.r = z__1.r;
        v5.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[6]);
        v6.r = z__1.r;
        v6.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[7]);
        v7.r = z__1.r;
        v7.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v7);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t7.r = z__1.r;
        t7.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[8]);
        v8.r = z__1.r;
        v8.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v8);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t8.r = z__1.r;
        t8.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[9]);
        v9.r = z__1.r;
        v9.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v9);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t9.r = z__1.r;
        t9.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            z__9.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__9.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j * c_dim1 + 2;
            z__10.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__10.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__8.r = z__9.r + z__10.r;
            z__8.i = z__9.i + z__10.i; // , expr subst
            i__4 = j * c_dim1 + 3;
            z__11.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__11.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__7.r = z__8.r + z__11.r;
            z__7.i = z__8.i + z__11.i; // , expr subst
            i__5 = j * c_dim1 + 4;
            z__12.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__12.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__6.r = z__7.r + z__12.r;
            z__6.i = z__7.i + z__12.i; // , expr subst
            i__6 = j * c_dim1 + 5;
            z__13.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__13.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__5.r = z__6.r + z__13.r;
            z__5.i = z__6.i + z__13.i; // , expr subst
            i__7 = j * c_dim1 + 6;
            z__14.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__14.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__4.r = z__5.r + z__14.r;
            z__4.i = z__5.i + z__14.i; // , expr subst
            i__8 = j * c_dim1 + 7;
            z__15.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i;
            z__15.i = v7.r * c__[i__8].i + v7.i * c__[i__8].r; // , expr subst
            z__3.r = z__4.r + z__15.r;
            z__3.i = z__4.i + z__15.i; // , expr subst
            i__9 = j * c_dim1 + 8;
            z__16.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i;
            z__16.i = v8.r * c__[i__9].i + v8.i * c__[i__9].r; // , expr subst
            z__2.r = z__3.r + z__16.r;
            z__2.i = z__3.i + z__16.i; // , expr subst
            i__10 = j * c_dim1 + 9;
            z__17.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i;
            z__17.i = v9.r * c__[i__10].i + v9.i * c__[i__10].r; // , expr subst
            z__1.r = z__2.r + z__17.r;
            z__1.i = z__2.i + z__17.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 7;
            i__3 = j * c_dim1 + 7;
            z__2.r = sum.r * t7.r - sum.i * t7.i;
            z__2.i = sum.r * t7.i + sum.i * t7.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 8;
            i__3 = j * c_dim1 + 8;
            z__2.r = sum.r * t8.r - sum.i * t8.i;
            z__2.i = sum.r * t8.i + sum.i * t8.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 9;
            i__3 = j * c_dim1 + 9;
            z__2.r = sum.r * t9.r - sum.i * t9.i;
            z__2.i = sum.r * t9.i + sum.i * t9.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L180: */
        }
        goto L410;
L190: /* Special code for 10 x 10 Householder */
        d_cnjg(&z__1, &v[1]);
        v1.r = z__1.r;
        v1.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[2]);
        v2.r = z__1.r;
        v2.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[3]);
        v3.r = z__1.r;
        v3.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[4]);
        v4.r = z__1.r;
        v4.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[5]);
        v5.r = z__1.r;
        v5.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[6]);
        v6.r = z__1.r;
        v6.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[7]);
        v7.r = z__1.r;
        v7.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v7);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t7.r = z__1.r;
        t7.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[8]);
        v8.r = z__1.r;
        v8.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v8);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t8.r = z__1.r;
        t8.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[9]);
        v9.r = z__1.r;
        v9.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v9);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t9.r = z__1.r;
        t9.i = z__1.i; // , expr subst
        d_cnjg(&z__1, &v[10]);
        v10.r = z__1.r;
        v10.i = z__1.i; // , expr subst
        d_cnjg(&z__2, &v10);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t10.r = z__1.r;
        t10.i = z__1.i; // , expr subst
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j * c_dim1 + 1;
            z__10.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__10.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j * c_dim1 + 2;
            z__11.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__11.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__9.r = z__10.r + z__11.r;
            z__9.i = z__10.i + z__11.i; // , expr subst
            i__4 = j * c_dim1 + 3;
            z__12.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__12.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__8.r = z__9.r + z__12.r;
            z__8.i = z__9.i + z__12.i; // , expr subst
            i__5 = j * c_dim1 + 4;
            z__13.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__13.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__7.r = z__8.r + z__13.r;
            z__7.i = z__8.i + z__13.i; // , expr subst
            i__6 = j * c_dim1 + 5;
            z__14.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__14.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__6.r = z__7.r + z__14.r;
            z__6.i = z__7.i + z__14.i; // , expr subst
            i__7 = j * c_dim1 + 6;
            z__15.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__15.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__5.r = z__6.r + z__15.r;
            z__5.i = z__6.i + z__15.i; // , expr subst
            i__8 = j * c_dim1 + 7;
            z__16.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i;
            z__16.i = v7.r * c__[i__8].i + v7.i * c__[i__8].r; // , expr subst
            z__4.r = z__5.r + z__16.r;
            z__4.i = z__5.i + z__16.i; // , expr subst
            i__9 = j * c_dim1 + 8;
            z__17.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i;
            z__17.i = v8.r * c__[i__9].i + v8.i * c__[i__9].r; // , expr subst
            z__3.r = z__4.r + z__17.r;
            z__3.i = z__4.i + z__17.i; // , expr subst
            i__10 = j * c_dim1 + 9;
            z__18.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i;
            z__18.i = v9.r * c__[i__10].i + v9.i * c__[i__10].r; // , expr subst
            z__2.r = z__3.r + z__18.r;
            z__2.i = z__3.i + z__18.i; // , expr subst
            i__11 = j * c_dim1 + 10;
            z__19.r = v10.r * c__[i__11].r - v10.i * c__[i__11].i;
            z__19.i = v10.r * c__[i__11].i + v10.i * c__[i__11].r; // , expr subst
            z__1.r = z__2.r + z__19.r;
            z__1.i = z__2.i + z__19.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 1;
            i__3 = j * c_dim1 + 1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 2;
            i__3 = j * c_dim1 + 2;
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 3;
            i__3 = j * c_dim1 + 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 4;
            i__3 = j * c_dim1 + 4;
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 5;
            i__3 = j * c_dim1 + 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 6;
            i__3 = j * c_dim1 + 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 7;
            i__3 = j * c_dim1 + 7;
            z__2.r = sum.r * t7.r - sum.i * t7.i;
            z__2.i = sum.r * t7.i + sum.i * t7.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 8;
            i__3 = j * c_dim1 + 8;
            z__2.r = sum.r * t8.r - sum.i * t8.i;
            z__2.i = sum.r * t8.i + sum.i * t8.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 9;
            i__3 = j * c_dim1 + 9;
            z__2.r = sum.r * t9.r - sum.i * t9.i;
            z__2.i = sum.r * t9.i + sum.i * t9.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j * c_dim1 + 10;
            i__3 = j * c_dim1 + 10;
            z__2.r = sum.r * t10.r - sum.i * t10.i;
            z__2.i = sum.r * t10.i + sum.i * t10.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L200: */
        }
        goto L410;
    }
    else
    {
        /* Form C * H, where H has order n. */
        switch (*n)
        {
        case 1:
            goto L210;
        case 2:
            goto L230;
        case 3:
            goto L250;
        case 4:
            goto L270;
        case 5:
            goto L290;
        case 6:
            goto L310;
        case 7:
            goto L330;
        case 8:
            goto L350;
        case 9:
            goto L370;
        case 10:
            goto L390;
        }
        /* Code for general N */
        zlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1]);
        goto L410;
L210: /* Special code for 1 x 1 Householder */
        z__3.r = tau->r * v[1].r - tau->i * v[1].i;
        z__3.i = tau->r * v[1].i + tau->i * v[1].r; // , expr subst
        d_cnjg(&z__4, &v[1]);
        z__2.r = z__3.r * z__4.r - z__3.i * z__4.i;
        z__2.i = z__3.r * z__4.i + z__3.i * z__4.r; // , expr subst
        z__1.r = 1. - z__2.r;
        z__1.i = 0. - z__2.i; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__1.r = t1.r * c__[i__3].r - t1.i * c__[i__3].i;
            z__1.i = t1.r * c__[i__3].i + t1.i * c__[i__3].r; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L220: */
        }
        goto L410;
L230: /* Special code for 2 x 2 Householder */
        v1.r = v[1].r;
        v1.i = v[1].i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        v2.r = v[2].r;
        v2.i = v[2].i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            z__2.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__2.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j + (c_dim1 << 1);
            z__3.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__3.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__1.r = z__2.r + z__3.r;
            z__1.i = z__2.i + z__3.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L240: */
        }
        goto L410;
L250: /* Special code for 3 x 3 Householder */
        v1.r = v[1].r;
        v1.i = v[1].i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        v2.r = v[2].r;
        v2.i = v[2].i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        v3.r = v[3].r;
        v3.i = v[3].i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            z__3.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__3.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j + (c_dim1 << 1);
            z__4.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__4.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__2.r = z__3.r + z__4.r;
            z__2.i = z__3.i + z__4.i; // , expr subst
            i__4 = j + c_dim1 * 3;
            z__5.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__5.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__1.r = z__2.r + z__5.r;
            z__1.i = z__2.i + z__5.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L260: */
        }
        goto L410;
L270: /* Special code for 4 x 4 Householder */
        v1.r = v[1].r;
        v1.i = v[1].i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        v2.r = v[2].r;
        v2.i = v[2].i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        v3.r = v[3].r;
        v3.i = v[3].i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        v4.r = v[4].r;
        v4.i = v[4].i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            z__4.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__4.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j + (c_dim1 << 1);
            z__5.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__5.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__3.r = z__4.r + z__5.r;
            z__3.i = z__4.i + z__5.i; // , expr subst
            i__4 = j + c_dim1 * 3;
            z__6.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__6.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__2.r = z__3.r + z__6.r;
            z__2.i = z__3.i + z__6.i; // , expr subst
            i__5 = j + (c_dim1 << 2);
            z__7.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__7.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__1.r = z__2.r + z__7.r;
            z__1.i = z__2.i + z__7.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L280: */
        }
        goto L410;
L290: /* Special code for 5 x 5 Householder */
        v1.r = v[1].r;
        v1.i = v[1].i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        v2.r = v[2].r;
        v2.i = v[2].i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        v3.r = v[3].r;
        v3.i = v[3].i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        v4.r = v[4].r;
        v4.i = v[4].i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        v5.r = v[5].r;
        v5.i = v[5].i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            z__5.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__5.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j + (c_dim1 << 1);
            z__6.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__6.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__4.r = z__5.r + z__6.r;
            z__4.i = z__5.i + z__6.i; // , expr subst
            i__4 = j + c_dim1 * 3;
            z__7.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__7.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__3.r = z__4.r + z__7.r;
            z__3.i = z__4.i + z__7.i; // , expr subst
            i__5 = j + (c_dim1 << 2);
            z__8.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__8.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__2.r = z__3.r + z__8.r;
            z__2.i = z__3.i + z__8.i; // , expr subst
            i__6 = j + c_dim1 * 5;
            z__9.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__9.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__1.r = z__2.r + z__9.r;
            z__1.i = z__2.i + z__9.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L300: */
        }
        goto L410;
L310: /* Special code for 6 x 6 Householder */
        v1.r = v[1].r;
        v1.i = v[1].i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        v2.r = v[2].r;
        v2.i = v[2].i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        v3.r = v[3].r;
        v3.i = v[3].i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        v4.r = v[4].r;
        v4.i = v[4].i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        v5.r = v[5].r;
        v5.i = v[5].i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        v6.r = v[6].r;
        v6.i = v[6].i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            z__6.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__6.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j + (c_dim1 << 1);
            z__7.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__7.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__5.r = z__6.r + z__7.r;
            z__5.i = z__6.i + z__7.i; // , expr subst
            i__4 = j + c_dim1 * 3;
            z__8.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__8.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__4.r = z__5.r + z__8.r;
            z__4.i = z__5.i + z__8.i; // , expr subst
            i__5 = j + (c_dim1 << 2);
            z__9.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__9.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__3.r = z__4.r + z__9.r;
            z__3.i = z__4.i + z__9.i; // , expr subst
            i__6 = j + c_dim1 * 5;
            z__10.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__10.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__2.r = z__3.r + z__10.r;
            z__2.i = z__3.i + z__10.i; // , expr subst
            i__7 = j + c_dim1 * 6;
            z__11.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__11.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__1.r = z__2.r + z__11.r;
            z__1.i = z__2.i + z__11.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L320: */
        }
        goto L410;
L330: /* Special code for 7 x 7 Householder */
        v1.r = v[1].r;
        v1.i = v[1].i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        v2.r = v[2].r;
        v2.i = v[2].i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        v3.r = v[3].r;
        v3.i = v[3].i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        v4.r = v[4].r;
        v4.i = v[4].i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        v5.r = v[5].r;
        v5.i = v[5].i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        v6.r = v[6].r;
        v6.i = v[6].i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        v7.r = v[7].r;
        v7.i = v[7].i; // , expr subst
        d_cnjg(&z__2, &v7);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t7.r = z__1.r;
        t7.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            z__7.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__7.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j + (c_dim1 << 1);
            z__8.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__8.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__6.r = z__7.r + z__8.r;
            z__6.i = z__7.i + z__8.i; // , expr subst
            i__4 = j + c_dim1 * 3;
            z__9.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__9.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__5.r = z__6.r + z__9.r;
            z__5.i = z__6.i + z__9.i; // , expr subst
            i__5 = j + (c_dim1 << 2);
            z__10.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__10.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__4.r = z__5.r + z__10.r;
            z__4.i = z__5.i + z__10.i; // , expr subst
            i__6 = j + c_dim1 * 5;
            z__11.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__11.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__3.r = z__4.r + z__11.r;
            z__3.i = z__4.i + z__11.i; // , expr subst
            i__7 = j + c_dim1 * 6;
            z__12.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__12.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__2.r = z__3.r + z__12.r;
            z__2.i = z__3.i + z__12.i; // , expr subst
            i__8 = j + c_dim1 * 7;
            z__13.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i;
            z__13.i = v7.r * c__[i__8].i + v7.i * c__[i__8].r; // , expr subst
            z__1.r = z__2.r + z__13.r;
            z__1.i = z__2.i + z__13.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 7;
            i__3 = j + c_dim1 * 7;
            z__2.r = sum.r * t7.r - sum.i * t7.i;
            z__2.i = sum.r * t7.i + sum.i * t7.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L340: */
        }
        goto L410;
L350: /* Special code for 8 x 8 Householder */
        v1.r = v[1].r;
        v1.i = v[1].i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        v2.r = v[2].r;
        v2.i = v[2].i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        v3.r = v[3].r;
        v3.i = v[3].i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        v4.r = v[4].r;
        v4.i = v[4].i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        v5.r = v[5].r;
        v5.i = v[5].i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        v6.r = v[6].r;
        v6.i = v[6].i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        v7.r = v[7].r;
        v7.i = v[7].i; // , expr subst
        d_cnjg(&z__2, &v7);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t7.r = z__1.r;
        t7.i = z__1.i; // , expr subst
        v8.r = v[8].r;
        v8.i = v[8].i; // , expr subst
        d_cnjg(&z__2, &v8);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t8.r = z__1.r;
        t8.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            z__8.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__8.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j + (c_dim1 << 1);
            z__9.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__9.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__7.r = z__8.r + z__9.r;
            z__7.i = z__8.i + z__9.i; // , expr subst
            i__4 = j + c_dim1 * 3;
            z__10.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__10.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__6.r = z__7.r + z__10.r;
            z__6.i = z__7.i + z__10.i; // , expr subst
            i__5 = j + (c_dim1 << 2);
            z__11.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__11.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__5.r = z__6.r + z__11.r;
            z__5.i = z__6.i + z__11.i; // , expr subst
            i__6 = j + c_dim1 * 5;
            z__12.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__12.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__4.r = z__5.r + z__12.r;
            z__4.i = z__5.i + z__12.i; // , expr subst
            i__7 = j + c_dim1 * 6;
            z__13.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__13.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__3.r = z__4.r + z__13.r;
            z__3.i = z__4.i + z__13.i; // , expr subst
            i__8 = j + c_dim1 * 7;
            z__14.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i;
            z__14.i = v7.r * c__[i__8].i + v7.i * c__[i__8].r; // , expr subst
            z__2.r = z__3.r + z__14.r;
            z__2.i = z__3.i + z__14.i; // , expr subst
            i__9 = j + (c_dim1 << 3);
            z__15.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i;
            z__15.i = v8.r * c__[i__9].i + v8.i * c__[i__9].r; // , expr subst
            z__1.r = z__2.r + z__15.r;
            z__1.i = z__2.i + z__15.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 7;
            i__3 = j + c_dim1 * 7;
            z__2.r = sum.r * t7.r - sum.i * t7.i;
            z__2.i = sum.r * t7.i + sum.i * t7.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 3);
            i__3 = j + (c_dim1 << 3);
            z__2.r = sum.r * t8.r - sum.i * t8.i;
            z__2.i = sum.r * t8.i + sum.i * t8.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L360: */
        }
        goto L410;
L370: /* Special code for 9 x 9 Householder */
        v1.r = v[1].r;
        v1.i = v[1].i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        v2.r = v[2].r;
        v2.i = v[2].i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        v3.r = v[3].r;
        v3.i = v[3].i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        v4.r = v[4].r;
        v4.i = v[4].i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        v5.r = v[5].r;
        v5.i = v[5].i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        v6.r = v[6].r;
        v6.i = v[6].i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        v7.r = v[7].r;
        v7.i = v[7].i; // , expr subst
        d_cnjg(&z__2, &v7);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t7.r = z__1.r;
        t7.i = z__1.i; // , expr subst
        v8.r = v[8].r;
        v8.i = v[8].i; // , expr subst
        d_cnjg(&z__2, &v8);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t8.r = z__1.r;
        t8.i = z__1.i; // , expr subst
        v9.r = v[9].r;
        v9.i = v[9].i; // , expr subst
        d_cnjg(&z__2, &v9);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t9.r = z__1.r;
        t9.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            z__9.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__9.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j + (c_dim1 << 1);
            z__10.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__10.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__8.r = z__9.r + z__10.r;
            z__8.i = z__9.i + z__10.i; // , expr subst
            i__4 = j + c_dim1 * 3;
            z__11.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__11.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__7.r = z__8.r + z__11.r;
            z__7.i = z__8.i + z__11.i; // , expr subst
            i__5 = j + (c_dim1 << 2);
            z__12.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__12.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__6.r = z__7.r + z__12.r;
            z__6.i = z__7.i + z__12.i; // , expr subst
            i__6 = j + c_dim1 * 5;
            z__13.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__13.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__5.r = z__6.r + z__13.r;
            z__5.i = z__6.i + z__13.i; // , expr subst
            i__7 = j + c_dim1 * 6;
            z__14.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__14.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__4.r = z__5.r + z__14.r;
            z__4.i = z__5.i + z__14.i; // , expr subst
            i__8 = j + c_dim1 * 7;
            z__15.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i;
            z__15.i = v7.r * c__[i__8].i + v7.i * c__[i__8].r; // , expr subst
            z__3.r = z__4.r + z__15.r;
            z__3.i = z__4.i + z__15.i; // , expr subst
            i__9 = j + (c_dim1 << 3);
            z__16.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i;
            z__16.i = v8.r * c__[i__9].i + v8.i * c__[i__9].r; // , expr subst
            z__2.r = z__3.r + z__16.r;
            z__2.i = z__3.i + z__16.i; // , expr subst
            i__10 = j + c_dim1 * 9;
            z__17.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i;
            z__17.i = v9.r * c__[i__10].i + v9.i * c__[i__10].r; // , expr subst
            z__1.r = z__2.r + z__17.r;
            z__1.i = z__2.i + z__17.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 7;
            i__3 = j + c_dim1 * 7;
            z__2.r = sum.r * t7.r - sum.i * t7.i;
            z__2.i = sum.r * t7.i + sum.i * t7.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 3);
            i__3 = j + (c_dim1 << 3);
            z__2.r = sum.r * t8.r - sum.i * t8.i;
            z__2.i = sum.r * t8.i + sum.i * t8.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 9;
            i__3 = j + c_dim1 * 9;
            z__2.r = sum.r * t9.r - sum.i * t9.i;
            z__2.i = sum.r * t9.i + sum.i * t9.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L380: */
        }
        goto L410;
L390: /* Special code for 10 x 10 Householder */
        v1.r = v[1].r;
        v1.i = v[1].i; // , expr subst
        d_cnjg(&z__2, &v1);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t1.r = z__1.r;
        t1.i = z__1.i; // , expr subst
        v2.r = v[2].r;
        v2.i = v[2].i; // , expr subst
        d_cnjg(&z__2, &v2);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t2.r = z__1.r;
        t2.i = z__1.i; // , expr subst
        v3.r = v[3].r;
        v3.i = v[3].i; // , expr subst
        d_cnjg(&z__2, &v3);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t3.r = z__1.r;
        t3.i = z__1.i; // , expr subst
        v4.r = v[4].r;
        v4.i = v[4].i; // , expr subst
        d_cnjg(&z__2, &v4);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t4.r = z__1.r;
        t4.i = z__1.i; // , expr subst
        v5.r = v[5].r;
        v5.i = v[5].i; // , expr subst
        d_cnjg(&z__2, &v5);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t5.r = z__1.r;
        t5.i = z__1.i; // , expr subst
        v6.r = v[6].r;
        v6.i = v[6].i; // , expr subst
        d_cnjg(&z__2, &v6);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t6.r = z__1.r;
        t6.i = z__1.i; // , expr subst
        v7.r = v[7].r;
        v7.i = v[7].i; // , expr subst
        d_cnjg(&z__2, &v7);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t7.r = z__1.r;
        t7.i = z__1.i; // , expr subst
        v8.r = v[8].r;
        v8.i = v[8].i; // , expr subst
        d_cnjg(&z__2, &v8);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t8.r = z__1.r;
        t8.i = z__1.i; // , expr subst
        v9.r = v[9].r;
        v9.i = v[9].i; // , expr subst
        d_cnjg(&z__2, &v9);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t9.r = z__1.r;
        t9.i = z__1.i; // , expr subst
        v10.r = v[10].r;
        v10.i = v[10].i; // , expr subst
        d_cnjg(&z__2, &v10);
        z__1.r = tau->r * z__2.r - tau->i * z__2.i;
        z__1.i = tau->r * z__2.i + tau->i * z__2.r; // , expr subst
        t10.r = z__1.r;
        t10.i = z__1.i; // , expr subst
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + c_dim1;
            z__10.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i;
            z__10.i = v1.r * c__[i__2].i + v1.i * c__[i__2].r; // , expr subst
            i__3 = j + (c_dim1 << 1);
            z__11.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i;
            z__11.i = v2.r * c__[i__3].i + v2.i * c__[i__3].r; // , expr subst
            z__9.r = z__10.r + z__11.r;
            z__9.i = z__10.i + z__11.i; // , expr subst
            i__4 = j + c_dim1 * 3;
            z__12.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i;
            z__12.i = v3.r * c__[i__4].i + v3.i * c__[i__4].r; // , expr subst
            z__8.r = z__9.r + z__12.r;
            z__8.i = z__9.i + z__12.i; // , expr subst
            i__5 = j + (c_dim1 << 2);
            z__13.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i;
            z__13.i = v4.r * c__[i__5].i + v4.i * c__[i__5].r; // , expr subst
            z__7.r = z__8.r + z__13.r;
            z__7.i = z__8.i + z__13.i; // , expr subst
            i__6 = j + c_dim1 * 5;
            z__14.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i;
            z__14.i = v5.r * c__[i__6].i + v5.i * c__[i__6].r; // , expr subst
            z__6.r = z__7.r + z__14.r;
            z__6.i = z__7.i + z__14.i; // , expr subst
            i__7 = j + c_dim1 * 6;
            z__15.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i;
            z__15.i = v6.r * c__[i__7].i + v6.i * c__[i__7].r; // , expr subst
            z__5.r = z__6.r + z__15.r;
            z__5.i = z__6.i + z__15.i; // , expr subst
            i__8 = j + c_dim1 * 7;
            z__16.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i;
            z__16.i = v7.r * c__[i__8].i + v7.i * c__[i__8].r; // , expr subst
            z__4.r = z__5.r + z__16.r;
            z__4.i = z__5.i + z__16.i; // , expr subst
            i__9 = j + (c_dim1 << 3);
            z__17.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i;
            z__17.i = v8.r * c__[i__9].i + v8.i * c__[i__9].r; // , expr subst
            z__3.r = z__4.r + z__17.r;
            z__3.i = z__4.i + z__17.i; // , expr subst
            i__10 = j + c_dim1 * 9;
            z__18.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i;
            z__18.i = v9.r * c__[i__10].i + v9.i * c__[i__10].r; // , expr subst
            z__2.r = z__3.r + z__18.r;
            z__2.i = z__3.i + z__18.i; // , expr subst
            i__11 = j + c_dim1 * 10;
            z__19.r = v10.r * c__[i__11].r - v10.i * c__[i__11].i;
            z__19.i = v10.r * c__[i__11].i + v10.i * c__[i__11].r; // , expr subst
            z__1.r = z__2.r + z__19.r;
            z__1.i = z__2.i + z__19.i; // , expr subst
            sum.r = z__1.r;
            sum.i = z__1.i; // , expr subst
            i__2 = j + c_dim1;
            i__3 = j + c_dim1;
            z__2.r = sum.r * t1.r - sum.i * t1.i;
            z__2.i = sum.r * t1.i + sum.i * t1.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 1);
            i__3 = j + (c_dim1 << 1);
            z__2.r = sum.r * t2.r - sum.i * t2.i;
            z__2.i = sum.r * t2.i + sum.i * t2.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 3;
            i__3 = j + c_dim1 * 3;
            z__2.r = sum.r * t3.r - sum.i * t3.i;
            z__2.i = sum.r * t3.i + sum.i * t3.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 2);
            i__3 = j + (c_dim1 << 2);
            z__2.r = sum.r * t4.r - sum.i * t4.i;
            z__2.i = sum.r * t4.i + sum.i * t4.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 5;
            i__3 = j + c_dim1 * 5;
            z__2.r = sum.r * t5.r - sum.i * t5.i;
            z__2.i = sum.r * t5.i + sum.i * t5.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 6;
            i__3 = j + c_dim1 * 6;
            z__2.r = sum.r * t6.r - sum.i * t6.i;
            z__2.i = sum.r * t6.i + sum.i * t6.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 7;
            i__3 = j + c_dim1 * 7;
            z__2.r = sum.r * t7.r - sum.i * t7.i;
            z__2.i = sum.r * t7.i + sum.i * t7.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + (c_dim1 << 3);
            i__3 = j + (c_dim1 << 3);
            z__2.r = sum.r * t8.r - sum.i * t8.i;
            z__2.i = sum.r * t8.i + sum.i * t8.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 9;
            i__3 = j + c_dim1 * 9;
            z__2.r = sum.r * t9.r - sum.i * t9.i;
            z__2.i = sum.r * t9.i + sum.i * t9.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            i__2 = j + c_dim1 * 10;
            i__3 = j + c_dim1 * 10;
            z__2.r = sum.r * t10.r - sum.i * t10.i;
            z__2.i = sum.r * t10.i + sum.i * t10.r; // , expr subst
            z__1.r = c__[i__3].r - z__2.r;
            z__1.i = c__[i__3].i - z__2.i; // , expr subst
            c__[i__2].r = z__1.r;
            c__[i__2].i = z__1.i; // , expr subst
            /* L400: */
        }
        goto L410;
    }
L410:
    return 0;
    /* End of ZLARFX */
}
/* zlarfx_ */
