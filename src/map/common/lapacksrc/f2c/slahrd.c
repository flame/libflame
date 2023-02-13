/* ../netlib/slahrd.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b4 = -1.f;
static real c_b5 = 1.f;
static integer c__1 = 1;
static real c_b38 = 0.f;
/* > \brief \b SLAHRD reduces the first nb columns of a general rectangular matrix A so that elements below th e k-th subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformati on to the unreduced part of A. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slahrd. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slahrd. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slahrd. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY ) */
/* .. Scalar Arguments .. */
/* INTEGER K, LDA, LDT, LDY, N, NB */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), T( LDT, NB ), TAU( NB ), */
/* $ Y( LDY, NB ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAHRD reduces the first NB columns of a real general n-by-(n-k+1) */
/* > matrix A so that elements below the k-th subdiagonal are zero. The */
/* > reduction is performed by an orthogonal similarity transformation */
/* > Q**T * A * Q. The routine returns the matrices V and T which determine */
/* > Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T. */
/* > */
/* > This is an OBSOLETE auxiliary routine. */
/* > This routine will be 'deprecated' in a future release. */
/* > Please use the new routine SLAHR2 instead. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The offset for the reduction. Elements below the k-th */
/* > subdiagonal in the first NB columns are reduced to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The number of columns to be reduced. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N-K+1) */
/* > On entry, the n-by-(n-k+1) general matrix A. */
/* > On exit, the elements on and above the k-th subdiagonal in */
/* > the first NB columns are overwritten with the corresponding */
/* > elements of the reduced matrix;
the elements below the k-th */
/* > subdiagonal, with the array TAU, represent the matrix Q as a */
/* > product of elementary reflectors. The other columns of A are */
/* > unchanged. See Further Details. */
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
/* > TAU is REAL array, dimension (NB) */
/* > The scalar factors of the elementary reflectors. See Further */
/* > Details. */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is REAL array, dimension (LDT,NB) */
/* > The upper triangular matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] Y */
/* > \verbatim */
/* > Y is REAL array, dimension (LDY,NB) */
/* > The n-by-nb matrix Y. */
/* > \endverbatim */
/* > */
/* > \param[in] LDY */
/* > \verbatim */
/* > LDY is INTEGER */
/* > The leading dimension of the array Y. LDY >= N. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of nb elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(nb). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar, and v is a real vector with */
/* > v(1:i+k-1) = 0, v(i+k) = 1;
v(i+k+1:n) is stored on exit in */
/* > A(i+k+1:n,i), and tau in TAU(i). */
/* > */
/* > The elements of the vectors v together form the (n-k+1)-by-nb matrix */
/* > V which is needed, with T and Y, to apply the transformation to the */
/* > unreduced part of the matrix, using an update of the form: */
/* > A := (I - V*T*V**T) * (A - Y*V**T). */
/* > */
/* > The contents of A on exit are illustrated by the following example */
/* > with n = 7, k = 3 and nb = 2: */
/* > */
/* > ( a h a a a ) */
/* > ( a h a a a ) */
/* > ( a h a a a ) */
/* > ( h h a a a ) */
/* > ( v1 h a a a ) */
/* > ( v1 v2 a a a ) */
/* > ( v1 v2 a a a ) */
/* > */
/* > where a denotes an element of the original matrix A, h denotes a */
/* > modified element of the upper Hessenberg matrix H, and vi denotes an */
/* > element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int slahrd_(integer *n, integer *k, integer *nb, real *a, integer *lda, real *tau, real *t, integer *ldt, real *y, integer *ldy)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, y_dim1, y_offset, i__1, i__2, i__3;
    real r__1;
    /* Local variables */
    integer i__;
    real ei;
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *), sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), scopy_( integer *, real *, integer *, real *, integer *), saxpy_(integer * , real *, real *, integer *, real *, integer *), strmv_(char *, char *, char *, integer *, real *, integer *, real *, integer *), slarfg_(integer *, real *, real *, integer *, real *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --tau;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    /* Function Body */
    if (*n <= 1)
    {
        return 0;
    }
    i__1 = *nb;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if (i__ > 1)
        {
            /* Update A(1:n,i) */
            /* Compute i-th column of A - Y * V**T */
            i__2 = i__ - 1;
            sgemv_("No transpose", n, &i__2, &c_b4, &y[y_offset], ldy, &a[*k + i__ - 1 + a_dim1], lda, &c_b5, &a[i__ * a_dim1 + 1], & c__1);
            /* Apply I - V * T**T * V**T to this column (call it b) from the */
            /* left, using the last column of T as workspace */
            /* Let V = ( V1 ) and b = ( b1 ) (first I-1 rows) */
            /* ( V2 ) ( b2 ) */
            /* where V1 is unit lower triangular */
            /* w := V1**T * b1 */
            i__2 = i__ - 1;
            scopy_(&i__2, &a[*k + 1 + i__ * a_dim1], &c__1, &t[*nb * t_dim1 + 1], &c__1);
            i__2 = i__ - 1;
            strmv_("Lower", "Transpose", "Unit", &i__2, &a[*k + 1 + a_dim1], lda, &t[*nb * t_dim1 + 1], &c__1);
            /* w := w + V2**T *b2 */
            i__2 = *n - *k - i__ + 1;
            i__3 = i__ - 1;
            sgemv_("Transpose", &i__2, &i__3, &c_b5, &a[*k + i__ + a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b5, &t[*nb * t_dim1 + 1], &c__1);
            /* w := T**T *w */
            i__2 = i__ - 1;
            strmv_("Upper", "Transpose", "Non-unit", &i__2, &t[t_offset], ldt, &t[*nb * t_dim1 + 1], &c__1);
            /* b2 := b2 - V2*w */
            i__2 = *n - *k - i__ + 1;
            i__3 = i__ - 1;
            sgemv_("No transpose", &i__2, &i__3, &c_b4, &a[*k + i__ + a_dim1], lda, &t[*nb * t_dim1 + 1], &c__1, &c_b5, &a[*k + i__ + i__ * a_dim1], &c__1);
            /* b1 := b1 - V1*w */
            i__2 = i__ - 1;
            strmv_("Lower", "No transpose", "Unit", &i__2, &a[*k + 1 + a_dim1] , lda, &t[*nb * t_dim1 + 1], &c__1);
            i__2 = i__ - 1;
            saxpy_(&i__2, &c_b4, &t[*nb * t_dim1 + 1], &c__1, &a[*k + 1 + i__ * a_dim1], &c__1);
            a[*k + i__ - 1 + (i__ - 1) * a_dim1] = ei;
        }
        /* Generate the elementary reflector H(i) to annihilate */
        /* A(k+i+1:n,i) */
        i__2 = *n - *k - i__ + 1;
        /* Computing MIN */
        i__3 = *k + i__ + 1;
        slarfg_(&i__2, &a[*k + i__ + i__ * a_dim1], &a[min(i__3,*n) + i__ * a_dim1], &c__1, &tau[i__]);
        ei = a[*k + i__ + i__ * a_dim1];
        a[*k + i__ + i__ * a_dim1] = 1.f;
        /* Compute Y(1:n,i) */
        i__2 = *n - *k - i__ + 1;
        sgemv_("No transpose", n, &i__2, &c_b5, &a[(i__ + 1) * a_dim1 + 1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b38, &y[i__ * y_dim1 + 1], &c__1);
        i__2 = *n - *k - i__ + 1;
        i__3 = i__ - 1;
        sgemv_("Transpose", &i__2, &i__3, &c_b5, &a[*k + i__ + a_dim1], lda, & a[*k + i__ + i__ * a_dim1], &c__1, &c_b38, &t[i__ * t_dim1 + 1], &c__1);
        i__2 = i__ - 1;
        sgemv_("No transpose", n, &i__2, &c_b4, &y[y_offset], ldy, &t[i__ * t_dim1 + 1], &c__1, &c_b5, &y[i__ * y_dim1 + 1], &c__1);
        sscal_(n, &tau[i__], &y[i__ * y_dim1 + 1], &c__1);
        /* Compute T(1:i,i) */
        i__2 = i__ - 1;
        r__1 = -tau[i__];
        sscal_(&i__2, &r__1, &t[i__ * t_dim1 + 1], &c__1);
        i__2 = i__ - 1;
        strmv_("Upper", "No transpose", "Non-unit", &i__2, &t[t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1) ;
        t[i__ + i__ * t_dim1] = tau[i__];
        /* L10: */
    }
    a[*k + *nb + *nb * a_dim1] = ei;
    return 0;
    /* End of SLAHRD */
}
/* slahrd_ */
