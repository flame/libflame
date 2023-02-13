/* ../netlib/chetd2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b2 =
{
    0.f,0.f
}
;
static integer c__1 = 1;
/* > \brief \b CHETD2 reduces a Hermitian matrix to real symmetric tridiagonal form by an unitary similarity t ransformation (unblocked algorithm). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHETD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetd2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetd2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetd2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHETD2( UPLO, N, A, LDA, D, E, TAU, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E( * ) */
/* COMPLEX A( LDA, * ), TAU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHETD2 reduces a complex Hermitian matrix A to real symmetric */
/* > tridiagonal form T by a unitary similarity transformation: */
/* > Q**H * A * Q = T. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > Hermitian matrix A is stored: */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
/* > n-by-n upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading n-by-n lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > On exit, if UPLO = 'U', the diagonal and first superdiagonal */
/* > of A are overwritten by the corresponding elements of the */
/* > tridiagonal matrix T, and the elements above the first */
/* > superdiagonal, with the array TAU, represent the unitary */
/* > matrix Q as a product of elementary reflectors;
if UPLO */
/* > = 'L', the diagonal and first subdiagonal of A are over- */
/* > written by the corresponding elements of the tridiagonal */
/* > matrix T, and the elements below the first subdiagonal, with */
/* > the array TAU, represent the unitary matrix Q as a product */
/* > of elementary reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The diagonal elements of the tridiagonal matrix T: */
/* > D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > The off-diagonal elements of the tridiagonal matrix T: */
/* > E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (N-1) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). */
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
/* > \ingroup complexHEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > If UPLO = 'U', the matrix Q is represented as a product of elementary */
/* > reflectors */
/* > */
/* > Q = H(n-1) . . . H(2) H(1). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a complex vector with */
/* > v(i+1:n) = 0 and v(i) = 1;
v(1:i-1) is stored on exit in */
/* > A(1:i-1,i+1), and tau in TAU(i). */
/* > */
/* > If UPLO = 'L', the matrix Q is represented as a product of elementary */
/* > reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(n-1). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a complex vector with */
/* > v(1:i) = 0 and v(i+1) = 1;
v(i+2:n) is stored on exit in A(i+2:n,i), */
/* > and tau in TAU(i). */
/* > */
/* > The contents of A on exit are illustrated by the following examples */
/* > with n = 5: */
/* > */
/* > if UPLO = 'U': if UPLO = 'L': */
/* > */
/* > ( d e v2 v3 v4 ) ( d ) */
/* > ( d e v3 v4 ) ( e d ) */
/* > ( d e v4 ) ( v1 e d ) */
/* > ( d e ) ( v1 v2 e d ) */
/* > ( d ) ( v1 v2 v3 e d ) */
/* > */
/* > where d and e denote diagonal and off-diagonal elements of T, and vi */
/* > denotes an element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int chetd2_fla(char *uplo, integer *n, complex *a, integer *lda, real *d__, real *e, complex *tau, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;
    complex q__1, q__2, q__3, q__4;
    /* Local variables */
    integer i__;
    complex taui;
    extern /* Subroutine */
    int cher2_(char *, integer *, complex *, complex * , integer *, complex *, integer *, complex *, integer *);
    complex alpha;
    extern /* Complex */
    VOID cdotc_f2c_(complex *, integer *, complex *, integer *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int chemv_(char *, integer *, complex *, complex * , integer *, complex *, integer *, complex *, complex *, integer * ), caxpy_(integer *, complex *, complex *, integer *, complex *, integer *);
    logical upper;
    extern /* Subroutine */
    int clarfg_(integer *, complex *, complex *, integer *, complex *), xerbla_(char *, integer *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < max(1,*n))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHETD2", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n <= 0)
    {
        return 0;
    }
    if (upper)
    {
        /* Reduce the upper triangle of A */
        i__1 = *n + *n * a_dim1;
        i__2 = *n + *n * a_dim1;
        r__1 = a[i__2].r;
        a[i__1].r = r__1;
        a[i__1].i = 0.f; // , expr subst
        for (i__ = *n - 1;
                i__ >= 1;
                --i__)
        {
            /* Generate elementary reflector H(i) = I - tau * v * v**H */
            /* to annihilate A(1:i-1,i+1) */
            i__1 = i__ + (i__ + 1) * a_dim1;
            alpha.r = a[i__1].r;
            alpha.i = a[i__1].i; // , expr subst
            clarfg_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &taui);
            i__1 = i__;
            e[i__1] = alpha.r;
            if (taui.r != 0.f || taui.i != 0.f)
            {
                /* Apply H(i) from both sides to A(1:i,1:i) */
                i__1 = i__ + (i__ + 1) * a_dim1;
                a[i__1].r = 1.f;
                a[i__1].i = 0.f; // , expr subst
                /* Compute x := tau * A * v storing x in TAU(1:i) */
                chemv_(uplo, &i__, &taui, &a[a_offset], lda, &a[(i__ + 1) * a_dim1 + 1], &c__1, &c_b2, &tau[1], &c__1);
                /* Compute w := x - 1/2 * tau * (x**H * v) * v */
                q__3.r = -.5f;
                q__3.i = -0.f; // , expr subst
                q__2.r = q__3.r * taui.r - q__3.i * taui.i;
                q__2.i = q__3.r * taui.i + q__3.i * taui.r; // , expr subst
                cdotc_f2c_(&q__4, &i__, &tau[1], &c__1, &a[(i__ + 1) * a_dim1 + 1] , &c__1);
                q__1.r = q__2.r * q__4.r - q__2.i * q__4.i;
                q__1.i = q__2.r * q__4.i + q__2.i * q__4.r; // , expr subst
                alpha.r = q__1.r;
                alpha.i = q__1.i; // , expr subst
                caxpy_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[ 1], &c__1);
                /* Apply the transformation as a rank-2 update: */
                /* A := A - v * w**H - w * v**H */
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cher2_(uplo, &i__, &q__1, &a[(i__ + 1) * a_dim1 + 1], &c__1, & tau[1], &c__1, &a[a_offset], lda);
            }
            else
            {
                i__1 = i__ + i__ * a_dim1;
                i__2 = i__ + i__ * a_dim1;
                r__1 = a[i__2].r;
                a[i__1].r = r__1;
                a[i__1].i = 0.f; // , expr subst
            }
            i__1 = i__ + (i__ + 1) * a_dim1;
            i__2 = i__;
            a[i__1].r = e[i__2];
            a[i__1].i = 0.f; // , expr subst
            i__1 = i__ + 1;
            i__2 = i__ + 1 + (i__ + 1) * a_dim1;
            d__[i__1] = a[i__2].r;
            i__1 = i__;
            tau[i__1].r = taui.r;
            tau[i__1].i = taui.i; // , expr subst
            /* L10: */
        }
        i__1 = a_dim1 + 1;
        d__[1] = a[i__1].r;
    }
    else
    {
        /* Reduce the lower triangle of A */
        i__1 = a_dim1 + 1;
        i__2 = a_dim1 + 1;
        r__1 = a[i__2].r;
        a[i__1].r = r__1;
        a[i__1].i = 0.f; // , expr subst
        i__1 = *n - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            /* Generate elementary reflector H(i) = I - tau * v * v**H */
            /* to annihilate A(i+2:n,i) */
            i__2 = i__ + 1 + i__ * a_dim1;
            alpha.r = a[i__2].r;
            alpha.i = a[i__2].i; // , expr subst
            i__2 = *n - i__;
            /* Computing MIN */
            i__3 = i__ + 2;
            clarfg_(&i__2, &alpha, &a[min(i__3,*n) + i__ * a_dim1], &c__1, & taui);
            i__2 = i__;
            e[i__2] = alpha.r;
            if (taui.r != 0.f || taui.i != 0.f)
            {
                /* Apply H(i) from both sides to A(i+1:n,i+1:n) */
                i__2 = i__ + 1 + i__ * a_dim1;
                a[i__2].r = 1.f;
                a[i__2].i = 0.f; // , expr subst
                /* Compute x := tau * A * v storing y in TAU(i:n-1) */
                i__2 = *n - i__;
                chemv_(uplo, &i__2, &taui, &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b2, &tau[ i__], &c__1);
                /* Compute w := x - 1/2 * tau * (x**H * v) * v */
                q__3.r = -.5f;
                q__3.i = -0.f; // , expr subst
                q__2.r = q__3.r * taui.r - q__3.i * taui.i;
                q__2.i = q__3.r * taui.i + q__3.i * taui.r; // , expr subst
                i__2 = *n - i__;
                cdotc_f2c_(&q__4, &i__2, &tau[i__], &c__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
                q__1.r = q__2.r * q__4.r - q__2.i * q__4.i;
                q__1.i = q__2.r * q__4.i + q__2.i * q__4.r; // , expr subst
                alpha.r = q__1.r;
                alpha.i = q__1.i; // , expr subst
                i__2 = *n - i__;
                caxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[ i__], &c__1);
                /* Apply the transformation as a rank-2 update: */
                /* A := A - v * w**H - w * v**H */
                i__2 = *n - i__;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cher2_(uplo, &i__2, &q__1, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__], &c__1, &a[i__ + 1 + (i__ + 1) * a_dim1], lda);
            }
            else
            {
                i__2 = i__ + 1 + (i__ + 1) * a_dim1;
                i__3 = i__ + 1 + (i__ + 1) * a_dim1;
                r__1 = a[i__3].r;
                a[i__2].r = r__1;
                a[i__2].i = 0.f; // , expr subst
            }
            i__2 = i__ + 1 + i__ * a_dim1;
            i__3 = i__;
            a[i__2].r = e[i__3];
            a[i__2].i = 0.f; // , expr subst
            i__2 = i__;
            i__3 = i__ + i__ * a_dim1;
            d__[i__2] = a[i__3].r;
            i__2 = i__;
            tau[i__2].r = taui.r;
            tau[i__2].i = taui.i; // , expr subst
            /* L20: */
        }
        i__1 = *n;
        i__2 = *n + *n * a_dim1;
        d__[i__1] = a[i__2].r;
    }
    return 0;
    /* End of CHETD2 */
}
/* chetd2_ */
