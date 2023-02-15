/* ../netlib/zlaqr1.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H a nd specified shifts. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAQR1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr1. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr1. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr1. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAQR1( N, H, LDH, S1, S2, V ) */
/* .. Scalar Arguments .. */
/* COMPLEX*16 S1, S2 */
/* INTEGER LDH, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 H( LDH, * ), V( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given a 2-by-2 or 3-by-3 matrix H, ZLAQR1 sets v to a */
/* > scalar multiple of the first column of the product */
/* > */
/* > (*) K = (H - s1*I)*(H - s2*I) */
/* > */
/* > scaling to avoid overflows and most underflows. */
/* > */
/* > This is useful for starting double implicit shift bulges */
/* > in the QR algorithm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is integer */
/* > Order of the matrix H. N must be either 2 or 3. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* > H is COMPLEX*16 array of dimension (LDH,N) */
/* > The 2-by-2 or 3-by-3 matrix H in (*). */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is integer */
/* > The leading dimension of H as declared in */
/* > the calling procedure. LDH.GE.N */
/* > \endverbatim */
/* > */
/* > \param[in] S1 */
/* > \verbatim */
/* > S1 is COMPLEX*16 */
/* > \endverbatim */
/* > */
/* > \param[in] S2 */
/* > \verbatim */
/* > S2 is COMPLEX*16 */
/* > */
/* > S1 and S2 are the shifts defining K in (*) above. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX*16 array of dimension N */
/* > A scalar multiple of the first column of the */
/* > matrix K in (*). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Karen Braman and Ralph Byers, Department of Mathematics, */
/* > University of Kansas, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
int zlaqr1_(integer *n, doublecomplex *h__, integer *ldh, doublecomplex *s1, doublecomplex *s2, doublecomplex *v)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    doublereal s;
    doublecomplex h21s, h31s;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --v;
    /* Function Body */
    if (*n == 2)
    {
        i__1 = h_dim1 + 1;
        z__2.r = h__[i__1].r - s2->r;
        z__2.i = h__[i__1].i - s2->i; // , expr subst
        z__1.r = z__2.r;
        z__1.i = z__2.i; // , expr subst
        i__2 = h_dim1 + 2;
        s = (d__1 = z__1.r, f2c_abs(d__1)) + (d__2 = d_imag(&z__1), f2c_abs(d__2)) + ( (d__3 = h__[i__2].r, f2c_abs(d__3)) + (d__4 = d_imag(&h__[h_dim1 + 2]), f2c_abs(d__4)));
        if (s == 0.)
        {
            v[1].r = 0.;
            v[1].i = 0.; // , expr subst
            v[2].r = 0.;
            v[2].i = 0.; // , expr subst
        }
        else
        {
            i__1 = h_dim1 + 2;
            z__1.r = h__[i__1].r / s;
            z__1.i = h__[i__1].i / s; // , expr subst
            h21s.r = z__1.r;
            h21s.i = z__1.i; // , expr subst
            i__1 = (h_dim1 << 1) + 1;
            z__2.r = h21s.r * h__[i__1].r - h21s.i * h__[i__1].i;
            z__2.i = h21s.r * h__[i__1].i + h21s.i * h__[i__1].r; // , expr subst
            i__2 = h_dim1 + 1;
            z__4.r = h__[i__2].r - s1->r;
            z__4.i = h__[i__2].i - s1->i; // , expr subst
            i__3 = h_dim1 + 1;
            z__6.r = h__[i__3].r - s2->r;
            z__6.i = h__[i__3].i - s2->i; // , expr subst
            z__5.r = z__6.r / s;
            z__5.i = z__6.i / s; // , expr subst
            z__3.r = z__4.r * z__5.r - z__4.i * z__5.i;
            z__3.i = z__4.r * z__5.i + z__4.i * z__5.r; // , expr subst
            z__1.r = z__2.r + z__3.r;
            z__1.i = z__2.i + z__3.i; // , expr subst
            v[1].r = z__1.r;
            v[1].i = z__1.i; // , expr subst
            i__1 = h_dim1 + 1;
            i__2 = (h_dim1 << 1) + 2;
            z__4.r = h__[i__1].r + h__[i__2].r;
            z__4.i = h__[i__1].i + h__[ i__2].i; // , expr subst
            z__3.r = z__4.r - s1->r;
            z__3.i = z__4.i - s1->i; // , expr subst
            z__2.r = z__3.r - s2->r;
            z__2.i = z__3.i - s2->i; // , expr subst
            z__1.r = h21s.r * z__2.r - h21s.i * z__2.i;
            z__1.i = h21s.r * z__2.i + h21s.i * z__2.r; // , expr subst
            v[2].r = z__1.r;
            v[2].i = z__1.i; // , expr subst
        }
    }
    else
    {
        i__1 = h_dim1 + 1;
        z__2.r = h__[i__1].r - s2->r;
        z__2.i = h__[i__1].i - s2->i; // , expr subst
        z__1.r = z__2.r;
        z__1.i = z__2.i; // , expr subst
        i__2 = h_dim1 + 2;
        i__3 = h_dim1 + 3;
        s = (d__1 = z__1.r, f2c_abs(d__1)) + (d__2 = d_imag(&z__1), f2c_abs(d__2)) + ( (d__3 = h__[i__2].r, f2c_abs(d__3)) + (d__4 = d_imag(&h__[h_dim1 + 2]), f2c_abs(d__4))) + ((d__5 = h__[i__3].r, f2c_abs(d__5)) + (d__6 = d_imag(&h__[h_dim1 + 3]), f2c_abs(d__6)));
        if (s == 0.)
        {
            v[1].r = 0.;
            v[1].i = 0.; // , expr subst
            v[2].r = 0.;
            v[2].i = 0.; // , expr subst
            v[3].r = 0.;
            v[3].i = 0.; // , expr subst
        }
        else
        {
            i__1 = h_dim1 + 2;
            z__1.r = h__[i__1].r / s;
            z__1.i = h__[i__1].i / s; // , expr subst
            h21s.r = z__1.r;
            h21s.i = z__1.i; // , expr subst
            i__1 = h_dim1 + 3;
            z__1.r = h__[i__1].r / s;
            z__1.i = h__[i__1].i / s; // , expr subst
            h31s.r = z__1.r;
            h31s.i = z__1.i; // , expr subst
            i__1 = h_dim1 + 1;
            z__4.r = h__[i__1].r - s1->r;
            z__4.i = h__[i__1].i - s1->i; // , expr subst
            i__2 = h_dim1 + 1;
            z__6.r = h__[i__2].r - s2->r;
            z__6.i = h__[i__2].i - s2->i; // , expr subst
            z__5.r = z__6.r / s;
            z__5.i = z__6.i / s; // , expr subst
            z__3.r = z__4.r * z__5.r - z__4.i * z__5.i;
            z__3.i = z__4.r * z__5.i + z__4.i * z__5.r; // , expr subst
            i__3 = (h_dim1 << 1) + 1;
            z__7.r = h__[i__3].r * h21s.r - h__[i__3].i * h21s.i;
            z__7.i = h__[i__3].r * h21s.i + h__[i__3].i * h21s.r; // , expr subst
            z__2.r = z__3.r + z__7.r;
            z__2.i = z__3.i + z__7.i; // , expr subst
            i__4 = h_dim1 * 3 + 1;
            z__8.r = h__[i__4].r * h31s.r - h__[i__4].i * h31s.i;
            z__8.i = h__[i__4].r * h31s.i + h__[i__4].i * h31s.r; // , expr subst
            z__1.r = z__2.r + z__8.r;
            z__1.i = z__2.i + z__8.i; // , expr subst
            v[1].r = z__1.r;
            v[1].i = z__1.i; // , expr subst
            i__1 = h_dim1 + 1;
            i__2 = (h_dim1 << 1) + 2;
            z__5.r = h__[i__1].r + h__[i__2].r;
            z__5.i = h__[i__1].i + h__[ i__2].i; // , expr subst
            z__4.r = z__5.r - s1->r;
            z__4.i = z__5.i - s1->i; // , expr subst
            z__3.r = z__4.r - s2->r;
            z__3.i = z__4.i - s2->i; // , expr subst
            z__2.r = h21s.r * z__3.r - h21s.i * z__3.i;
            z__2.i = h21s.r * z__3.i + h21s.i * z__3.r; // , expr subst
            i__3 = h_dim1 * 3 + 2;
            z__6.r = h__[i__3].r * h31s.r - h__[i__3].i * h31s.i;
            z__6.i = h__[i__3].r * h31s.i + h__[i__3].i * h31s.r; // , expr subst
            z__1.r = z__2.r + z__6.r;
            z__1.i = z__2.i + z__6.i; // , expr subst
            v[2].r = z__1.r;
            v[2].i = z__1.i; // , expr subst
            i__1 = h_dim1 + 1;
            i__2 = h_dim1 * 3 + 3;
            z__5.r = h__[i__1].r + h__[i__2].r;
            z__5.i = h__[i__1].i + h__[ i__2].i; // , expr subst
            z__4.r = z__5.r - s1->r;
            z__4.i = z__5.i - s1->i; // , expr subst
            z__3.r = z__4.r - s2->r;
            z__3.i = z__4.i - s2->i; // , expr subst
            z__2.r = h31s.r * z__3.r - h31s.i * z__3.i;
            z__2.i = h31s.r * z__3.i + h31s.i * z__3.r; // , expr subst
            i__3 = (h_dim1 << 1) + 3;
            z__6.r = h21s.r * h__[i__3].r - h21s.i * h__[i__3].i;
            z__6.i = h21s.r * h__[i__3].i + h21s.i * h__[i__3].r; // , expr subst
            z__1.r = z__2.r + z__6.r;
            z__1.i = z__2.i + z__6.i; // , expr subst
            v[3].r = z__1.r;
            v[3].i = z__1.i; // , expr subst
        }
    }
    return 0;
}
/* zlaqr1_ */
