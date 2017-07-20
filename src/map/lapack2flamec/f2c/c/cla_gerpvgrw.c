/* ../netlib/cla_gerpvgrw.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLA_GERPVGRW multiplies a square real matrix by a complex matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLA_GERPVGRW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_ger pvgrw.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_ger pvgrw.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_ger pvgrw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLA_GERPVGRW( N, NCOLS, A, LDA, AF, LDAF ) */
/* .. Scalar Arguments .. */
/* INTEGER N, NCOLS, LDA, LDAF */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), AF( LDAF, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > CLA_GERPVGRW computes the reciprocal pivot growth factor */
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
/* > A is COMPLEX array, dimension (LDA,N) */
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
/* > AF is COMPLEX array, dimension (LDAF,N) */
/* > The factors L and U from the factorization */
/* > A = P*L*U as computed by CGETRF. */
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
/* > \ingroup complexGEcomputational */
/* ===================================================================== */
real cla_gerpvgrw_(integer *n, integer *ncols, complex *a, integer *lda, complex *af, integer *ldaf)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3;
    real ret_val, r__1, r__2, r__3;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    integer i__, j;
    real amax, umax, rpvgrw;
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
    rpvgrw = 1.f;
    i__1 = *ncols;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        amax = 0.f;
        umax = 0.f;
        i__2 = *n;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            /* Computing MAX */
            i__3 = i__ + j * a_dim1;
            r__3 = (r__1 = a[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[i__ + j * a_dim1]), f2c_abs(r__2));
            amax = max(r__3,amax);
        }
        i__2 = j;
        for (i__ = 1;
                i__ <= i__2;
                ++i__)
        {
            /* Computing MAX */
            i__3 = i__ + j * af_dim1;
            r__3 = (r__1 = af[i__3].r, f2c_abs(r__1)) + (r__2 = r_imag(&af[i__ + j * af_dim1]), f2c_abs(r__2));
            umax = max(r__3,umax);
        }
        if (umax != 0.f)
        {
            /* Computing MIN */
            r__1 = amax / umax;
            rpvgrw = min(r__1,rpvgrw);
        }
    }
    ret_val = rpvgrw;
    return ret_val;
}
/* cla_gerpvgrw__ */
