/* ../netlib/dlarz.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b5 = 1.;
/* > \brief \b DLARZ applies an elementary reflector (as returned by stzrzf) to a general matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarz.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarz.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarz.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER INCV, L, LDC, M, N */
/* DOUBLE PRECISION TAU */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARZ applies a real elementary reflector H to a real M-by-N */
/* > matrix C, from either the left or the right. H is represented in the */
/* > form */
/* > */
/* > H = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar and v is a real vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > */
/* > */
/* > H is a product of k elementary reflectors as returned by DTZRZF. */
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
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > The number of entries of the vector V containing */
/* > the meaningful part of the Householder vectors. */
/* > If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is DOUBLE PRECISION array, dimension (1+(L-1)*f2c_abs(INCV)) */
/* > The vector v in the representation of H as returned by */
/* > DTZRZF. V is not used if TAU = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* > INCV is INTEGER */
/* > The increment between elements of v. INCV <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* > or C * H if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension */
/* > (N) if SIDE = 'L' */
/* > or (M) if SIDE = 'R' */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int dlarz_(char *side, integer *m, integer *n, integer *l, doublereal *v, integer *incv, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work)
{
    /* System generated locals */
    integer c_dim1, c_offset;
    doublereal d__1;
    /* Local variables */
    extern /* Subroutine */
    int dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), dcopy_(integer *, doublereal *, integer *, doublereal *, integer *), daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *) ;
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    if (lsame_(side, "L"))
    {
        /* Form H * C */
        if (*tau != 0.)
        {
            /* w( 1:n ) = C( 1, 1:n ) */
            dcopy_(n, &c__[c_offset], ldc, &work[1], &c__1);
            /* w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )**T * v( 1:l ) */
            dgemv_("Transpose", l, n, &c_b5, &c__[*m - *l + 1 + c_dim1], ldc, &v[1], incv, &c_b5, &work[1], &c__1);
            /* C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n ) */
            d__1 = -(*tau);
            daxpy_(n, &d__1, &work[1], &c__1, &c__[c_offset], ldc);
            /* C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ... */
            /* tau * v( 1:l ) * w( 1:n )**T */
            d__1 = -(*tau);
            dger_(l, n, &d__1, &v[1], incv, &work[1], &c__1, &c__[*m - *l + 1 + c_dim1], ldc);
        }
    }
    else
    {
        /* Form C * H */
        if (*tau != 0.)
        {
            /* w( 1:m ) = C( 1:m, 1 ) */
            dcopy_(m, &c__[c_offset], &c__1, &work[1], &c__1);
            /* w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l ) */
            dgemv_("No transpose", m, l, &c_b5, &c__[(*n - *l + 1) * c_dim1 + 1], ldc, &v[1], incv, &c_b5, &work[1], &c__1);
            /* C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m ) */
            d__1 = -(*tau);
            daxpy_(m, &d__1, &work[1], &c__1, &c__[c_offset], &c__1);
            /* C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ... */
            /* tau * w( 1:m ) * v( 1:l )**T */
            d__1 = -(*tau);
            dger_(m, l, &d__1, &work[1], &c__1, &v[1], incv, &c__[(*n - *l + 1) * c_dim1 + 1], ldc);
        }
    }
    return 0;
    /* End of DLARZ */
}
/* dlarz_ */
