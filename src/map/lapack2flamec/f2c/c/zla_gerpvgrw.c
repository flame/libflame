/* ../netlib/zla_gerpvgrw.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLA_GERPVGRW multiplies a square real matrix by a complex matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLA_GERPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_ger pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_ger pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_ger pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* DOUBLE PRECISION FUNCTION ZLA_GERPVGRW( N, NCOLS, A, LDA, AF, */
/* LDAF ) */
/* .. Scalar Arguments .. */
/* INTEGER N, NCOLS, LDA, LDAF */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), AF( LDAF, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > ZLA_GERPVGRW computes the reciprocal pivot growth factor */
/* > norm(A)/norm(U). The "max absolute element" norm is used. If this is */
/* > much less than 1, the stability of the LU factorization of the */
/* > (equilibrated) matrix A could be poor. This also means that the */
/* > solution X, estimated condition numbers, and error bounds could be */
/* > unreliable. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCOLS */
/* > \verbatim */
/* > NCOLS is INTEGER */
/* > The number of columns of the matrix A. NCOLS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* > AF is DOUBLE PRECISION array, dimension (LDAF,N) */
/* > The factors L and U from the factorization */
/* > A = P*L*U as computed by ZGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= max(1,N). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16GEcomputational */
/* ===================================================================== */
doublereal zla_gerpvgrw_(integer *n, integer *ncols, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3;
    doublereal ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    integer i__, j;
    doublereal amax, umax, rpvgrw;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function Definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    /* Function Body */
    rpvgrw = 1.;
    i__1 = *ncols;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        amax = 0.;
        umax = 0.;
        i__2 = *n;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            /* Computing MAX */
            i__3 = i__ + j * a_dim1;
            d__3 = (d__1 = a[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[i__ + j * a_dim1]), f2c_abs(d__2));
            amax = max(d__3,amax);
        }
        i__2 = j;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            /* Computing MAX */
            i__3 = i__ + j * af_dim1;
            d__3 = (d__1 = af[i__3].r, f2c_abs(d__1)) + (d__2 = d_imag(&af[i__ + j * af_dim1]), f2c_abs(d__2));
            umax = max(d__3,umax);
        }
        if (umax != 0.)
        {
            /* Computing MIN */
            d__1 = amax / umax;
            rpvgrw = min(d__1,rpvgrw);
        }
    }
    ret_val = rpvgrw;
    return ret_val;
}
/* zla_gerpvgrw__ */
