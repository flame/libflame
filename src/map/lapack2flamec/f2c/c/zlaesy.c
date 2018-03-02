/* ../netlib/zlaesy.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    1.,0.
}
;
static integer c__2 = 2;
/* > \brief \b ZLAESY computes the eigenvalues and eigenvectors of a 2-by-2 complex symmetric matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAESY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaesy. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaesy. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaesy. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 ) */
/* .. Scalar Arguments .. */
/* COMPLEX*16 A, B, C, CS1, EVSCAL, RT1, RT2, SN1 */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAESY computes the eigendecomposition of a 2-by-2 symmetric matrix */
/* > ( ( A, B );
( B, C ) ) */
/* > provided the norm of the matrix of eigenvectors is larger than */
/* > some threshold value. */
/* > */
/* > RT1 is the eigenvalue of larger absolute value, and RT2 of */
/* > smaller absolute value. If the eigenvectors are computed, then */
/* > on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence */
/* > */
/* > [ CS1 SN1 ] . [ A B ] . [ CS1 -SN1 ] = [ RT1 0 ] */
/* > [ -SN1 CS1 ] [ B C ] [ SN1 CS1 ] [ 0 RT2 ] */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 */
/* > The ( 1, 1 ) element of input matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX*16 */
/* > The ( 1, 2 ) element of input matrix. The ( 2, 1 ) element */
/* > is also given by B, since the 2-by-2 matrix is symmetric. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is COMPLEX*16 */
/* > The ( 2, 2 ) element of input matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1 */
/* > \verbatim */
/* > RT1 is COMPLEX*16 */
/* > The eigenvalue of larger modulus. */
/* > \endverbatim */
/* > */
/* > \param[out] RT2 */
/* > \verbatim */
/* > RT2 is COMPLEX*16 */
/* > The eigenvalue of smaller modulus. */
/* > \endverbatim */
/* > */
/* > \param[out] EVSCAL */
/* > \verbatim */
/* > EVSCAL is COMPLEX*16 */
/* > The complex value by which the eigenvector matrix was scaled */
/* > to make it orthonormal. If EVSCAL is zero, the eigenvectors */
/* > were not computed. This means one of two things: the 2-by-2 */
/* > matrix could not be diagonalized, or the norm of the matrix */
/* > of eigenvectors before scaling was larger than the threshold */
/* > value THRESH (set below). */
/* > \endverbatim */
/* > */
/* > \param[out] CS1 */
/* > \verbatim */
/* > CS1 is COMPLEX*16 */
/* > \endverbatim */
/* > */
/* > \param[out] SN1 */
/* > \verbatim */
/* > SN1 is COMPLEX*16 */
/* > If EVSCAL .NE. 0, ( CS1, SN1 ) is the unit right eigenvector */
/* > for RT1. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16SYauxiliary */
/* ===================================================================== */
/* Subroutine */
int zlaesy_(doublecomplex *a, doublecomplex *b, doublecomplex *c__, doublecomplex *rt1, doublecomplex *rt2, doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1)
{
    /* System generated locals */
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;
    /* Builtin functions */
    double z_abs(doublecomplex *);
    void pow_zi(doublecomplex *, doublecomplex *, integer *), z_sqrt( doublecomplex *, doublecomplex *), z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    /* Local variables */
    doublecomplex s, t;
    doublereal z__;
    doublecomplex tmp;
    doublereal babs, tabs, evnorm;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Special case: The matrix is actually diagonal. */
    /* To avoid divide by zero later, we treat this case separately. */
    if (z_abs(b) == 0.)
    {
        rt1->r = a->r, rt1->i = a->i;
        rt2->r = c__->r, rt2->i = c__->i;
        if (z_abs(rt1) < z_abs(rt2))
        {
            tmp.r = rt1->r;
            tmp.i = rt1->i; // , expr subst
            rt1->r = rt2->r, rt1->i = rt2->i;
            rt2->r = tmp.r, rt2->i = tmp.i;
            cs1->r = 0., cs1->i = 0.;
            sn1->r = 1., sn1->i = 0.;
        }
        else
        {
            cs1->r = 1., cs1->i = 0.;
            sn1->r = 0., sn1->i = 0.;
        }
    }
    else
    {
        /* Compute the eigenvalues and eigenvectors. */
        /* The characteristic equation is */
        /* lambda **2 - (A+C) lambda + (A*C - B*B) */
        /* and we solve it using the quadratic formula. */
        z__2.r = a->r + c__->r;
        z__2.i = a->i + c__->i; // , expr subst
        z__1.r = z__2.r * .5;
        z__1.i = z__2.i * .5; // , expr subst
        s.r = z__1.r;
        s.i = z__1.i; // , expr subst
        z__2.r = a->r - c__->r;
        z__2.i = a->i - c__->i; // , expr subst
        z__1.r = z__2.r * .5;
        z__1.i = z__2.i * .5; // , expr subst
        t.r = z__1.r;
        t.i = z__1.i; // , expr subst
        /* Take the square root carefully to avoid over/under flow. */
        babs = z_abs(b);
        tabs = z_abs(&t);
        z__ = max(babs,tabs);
        if (z__ > 0.)
        {
            z__5.r = t.r / z__;
            z__5.i = t.i / z__; // , expr subst
            pow_zi(&z__4, &z__5, &c__2);
            z__7.r = b->r / z__;
            z__7.i = b->i / z__; // , expr subst
            pow_zi(&z__6, &z__7, &c__2);
            z__3.r = z__4.r + z__6.r;
            z__3.i = z__4.i + z__6.i; // , expr subst
            z_sqrt(&z__2, &z__3);
            z__1.r = z__ * z__2.r;
            z__1.i = z__ * z__2.i; // , expr subst
            t.r = z__1.r;
            t.i = z__1.i; // , expr subst
        }
        /* Compute the two eigenvalues. RT1 and RT2 are exchanged */
        /* if necessary so that RT1 will have the greater magnitude. */
        z__1.r = s.r + t.r;
        z__1.i = s.i + t.i; // , expr subst
        rt1->r = z__1.r, rt1->i = z__1.i;
        z__1.r = s.r - t.r;
        z__1.i = s.i - t.i; // , expr subst
        rt2->r = z__1.r, rt2->i = z__1.i;
        if (z_abs(rt1) < z_abs(rt2))
        {
            tmp.r = rt1->r;
            tmp.i = rt1->i; // , expr subst
            rt1->r = rt2->r, rt1->i = rt2->i;
            rt2->r = tmp.r, rt2->i = tmp.i;
        }
        /* Choose CS1 = 1 and SN1 to satisfy the first equation, then */
        /* scale the components of this eigenvector so that the matrix */
        /* of eigenvectors X satisfies X * X**T = I . (No scaling is */
        /* done if the norm of the eigenvalue matrix is less than THRESH.) */
        z__2.r = rt1->r - a->r;
        z__2.i = rt1->i - a->i; // , expr subst
        z_div(&z__1, &z__2, b);
        sn1->r = z__1.r, sn1->i = z__1.i;
        tabs = z_abs(sn1);
        if (tabs > 1.)
        {
            /* Computing 2nd power */
            d__2 = 1. / tabs;
            d__1 = d__2 * d__2;
            z__5.r = sn1->r / tabs;
            z__5.i = sn1->i / tabs; // , expr subst
            pow_zi(&z__4, &z__5, &c__2);
            z__3.r = d__1 + z__4.r;
            z__3.i = z__4.i; // , expr subst
            z_sqrt(&z__2, &z__3);
            z__1.r = tabs * z__2.r;
            z__1.i = tabs * z__2.i; // , expr subst
            t.r = z__1.r;
            t.i = z__1.i; // , expr subst
        }
        else
        {
            z__3.r = sn1->r * sn1->r - sn1->i * sn1->i;
            z__3.i = sn1->r * sn1->i + sn1->i * sn1->r; // , expr subst
            z__2.r = z__3.r + 1.;
            z__2.i = z__3.i + 0.; // , expr subst
            z_sqrt(&z__1, &z__2);
            t.r = z__1.r;
            t.i = z__1.i; // , expr subst
        }
        evnorm = z_abs(&t);
        if (evnorm >= .1)
        {
            z_div(&z__1, &c_b1, &t);
            evscal->r = z__1.r, evscal->i = z__1.i;
            cs1->r = evscal->r, cs1->i = evscal->i;
            z__1.r = sn1->r * evscal->r - sn1->i * evscal->i;
            z__1.i = sn1->r * evscal->i + sn1->i * evscal->r; // , expr subst
            sn1->r = z__1.r, sn1->i = z__1.i;
        }
        else
        {
            evscal->r = 0., evscal->i = 0.;
        }
    }
    return 0;
    /* End of ZLAESY */
}
/* zlaesy_ */
