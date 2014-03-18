/* ../netlib/zgehd2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZGEHD2 reduces a general square matrix to upper Hessenberg form using an unblocked algorithm. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGEHD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgehd2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgehd2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgehd2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER IHI, ILO, INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEHD2 reduces a complex general matrix A to upper Hessenberg form H */
/* > by a unitary similarity transformation: Q**H * A * Q = H . */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > */
/* > It is assumed that A is already upper triangular in rows */
/* > and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
/* > set by a previous call to ZGEBAL;
otherwise they should be */
/* > set to 1 and N respectively. See Further Details. */
/* > 1 <= ILO <= IHI <= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the n by n general matrix to be reduced. */
/* > On exit, the upper triangle and the first subdiagonal of A */
/* > are overwritten with the upper Hessenberg matrix H, and the */
/* > elements below the first subdiagonal, with the array TAU, */
/* > represent the unitary matrix Q as a product of elementary */
/* > reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (N-1) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16GEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of (ihi-ilo) elementary */
/* > reflectors */
/* > */
/* > Q = H(ilo) H(ilo+1) . . . H(ihi-1). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a complex vector with */
/* > v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0;
v(i+2:ihi) is stored on */
/* > exit in A(i+2:ihi,i), and tau in TAU(i). */
/* > */
/* > The contents of A are illustrated by the following example, with */
/* > n = 7, ilo = 2 and ihi = 6: */
/* > */
/* > on entry, on exit, */
/* > */
/* > ( a a a a a a a ) ( a a h h h h a ) */
/* > ( a a a a a a ) ( a h h h h a ) */
/* > ( a a a a a a ) ( h h h h h h ) */
/* > ( a a a a a a ) ( v2 h h h h h ) */
/* > ( a a a a a a ) ( v2 v3 h h h h ) */
/* > ( a a a a a a ) ( v2 v3 v4 h h h ) */
/* > ( a ) ( a ) */
/* > */
/* > where a denotes an element of the original matrix A, h denotes a */
/* > modified element of the upper Hessenberg matrix H, and vi denotes an */
/* > element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int zgehd2_(integer *n, integer *ilo, integer *ihi, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex * work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    integer i__;
    doublecomplex alpha;
    extern /* Subroutine */
    int zlarf_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *), xerbla_(char *, integer *), zlarfg_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *);
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
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*ilo < 1 || *ilo > max(1,*n))
    {
        *info = -2;
    }
    else if (*ihi < min(*ilo,*n) || *ihi > *n)
    {
        *info = -3;
    }
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGEHD2", &i__1);
        return 0;
    }
    i__1 = *ihi - 1;
    for (i__ = *ilo;
            i__ <= i__1;
            ++i__)
    {
        /* Compute elementary reflector H(i) to annihilate A(i+2:ihi,i) */
        i__2 = i__ + 1 + i__ * a_dim1;
        alpha.r = a[i__2].r;
        alpha.i = a[i__2].i; // , expr subst
        i__2 = *ihi - i__;
        /* Computing MIN */
        i__3 = i__ + 2;
        zlarfg_(&i__2, &alpha, &a[min(i__3,*n) + i__ * a_dim1], &c__1, &tau[ i__]);
        i__2 = i__ + 1 + i__ * a_dim1;
        a[i__2].r = 1.;
        a[i__2].i = 0.; // , expr subst
        /* Apply H(i) to A(1:ihi,i+1:ihi) from the right */
        i__2 = *ihi - i__;
        zlarf_("Right", ihi, &i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[ i__], &a[(i__ + 1) * a_dim1 + 1], lda, &work[1]);
        /* Apply H(i)**H to A(i+1:ihi,i+1:n) from the left */
        i__2 = *ihi - i__;
        i__3 = *n - i__;
        d_cnjg(&z__1, &tau[i__]);
        zlarf_("Left", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &c__1, &z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &work[1]);
        i__2 = i__ + 1 + i__ * a_dim1;
        a[i__2].r = alpha.r;
        a[i__2].i = alpha.i; // , expr subst
        /* L10: */
    }
    return 0;
    /* End of ZGEHD2 */
}
/* zgehd2_ */
