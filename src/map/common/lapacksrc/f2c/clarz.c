/* ../netlib/clarz.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    1.f,0.f
}
;
static integer c__1 = 1;
/* > \brief \b CLARZ applies an elementary reflector (as returned by stzrzf) to a general matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARZ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarz.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarz.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarz.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE */
/* INTEGER INCV, L, LDC, M, N */
/* COMPLEX TAU */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX C( LDC, * ), V( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARZ applies a complex elementary reflector H to a complex */
/* > M-by-N matrix C, from either the left or the right. H is represented */
/* > in the form */
/* > */
/* > H = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar and v is a complex vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > */
/* > To apply H**H (the conjugate transpose of H), supply conjg(tau) instead */
/* > tau. */
/* > */
/* > H is a product of k elementary reflectors as returned by CTZRZF. */
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
/* > V is COMPLEX array, dimension (1+(L-1)*f2c_abs(INCV)) */
/* > The vector v in the representation of H as returned by */
/* > CTZRZF. V is not used if TAU = 0. */
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
/* > TAU is COMPLEX */
/* > The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX array, dimension (LDC,N) */
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
/* > WORK is COMPLEX array, dimension */
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
/* > \ingroup complexOTHERcomputational */
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
int clarz_(char *side, integer *m, integer *n, integer *l, complex *v, integer *incv, complex *tau, complex *c__, integer *ldc, complex *work)
{
    /* System generated locals */
    integer c_dim1, c_offset;
    complex q__1;
    /* Local variables */
    extern /* Subroutine */
    int cgerc_(integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, integer *), cgemv_(char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int cgeru_(integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, integer *), ccopy_(integer *, complex *, integer *, complex *, integer *), caxpy_(integer *, complex *, complex *, integer *, complex *, integer *), clacgv_(integer *, complex *, integer *);
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
        if (tau->r != 0.f || tau->i != 0.f)
        {
            /* w( 1:n ) = conjg( C( 1, 1:n ) ) */
            ccopy_(n, &c__[c_offset], ldc, &work[1], &c__1);
            clacgv_(n, &work[1], &c__1);
            /* w( 1:n ) = conjg( w( 1:n ) + C( m-l+1:m, 1:n )**H * v( 1:l ) ) */
            cgemv_("Conjugate transpose", l, n, &c_b1, &c__[*m - *l + 1 + c_dim1], ldc, &v[1], incv, &c_b1, &work[1], &c__1);
            clacgv_(n, &work[1], &c__1);
            /* C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n ) */
            q__1.r = -tau->r;
            q__1.i = -tau->i; // , expr subst
            caxpy_(n, &q__1, &work[1], &c__1, &c__[c_offset], ldc);
            /* C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ... */
            /* tau * v( 1:l ) * w( 1:n )**H */
            q__1.r = -tau->r;
            q__1.i = -tau->i; // , expr subst
            cgeru_(l, n, &q__1, &v[1], incv, &work[1], &c__1, &c__[*m - *l + 1 + c_dim1], ldc);
        }
    }
    else
    {
        /* Form C * H */
        if (tau->r != 0.f || tau->i != 0.f)
        {
            /* w( 1:m ) = C( 1:m, 1 ) */
            ccopy_(m, &c__[c_offset], &c__1, &work[1], &c__1);
            /* w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l ) */
            cgemv_("No transpose", m, l, &c_b1, &c__[(*n - *l + 1) * c_dim1 + 1], ldc, &v[1], incv, &c_b1, &work[1], &c__1);
            /* C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m ) */
            q__1.r = -tau->r;
            q__1.i = -tau->i; // , expr subst
            caxpy_(m, &q__1, &work[1], &c__1, &c__[c_offset], &c__1);
            /* C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ... */
            /* tau * w( 1:m ) * v( 1:l )**H */
            q__1.r = -tau->r;
            q__1.i = -tau->i; // , expr subst
            cgerc_(m, l, &q__1, &work[1], &c__1, &v[1], incv, &c__[(*n - *l + 1) * c_dim1 + 1], ldc);
        }
    }
    return 0;
    /* End of CLARZ */
}
/* clarz_ */
