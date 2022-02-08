/* ../netlib/v3.9.0/ztplqt2.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    0.,0.
}
;
static doublecomplex c_b2 =
{
    1.,0.
}
;
/* > \brief \b ZTPLQT2 computes a LQ factorization of a real or complex "triangular-pentagonal" matrix, which is composed of a triangular block and a pentagonal block, using the compact WY representation for Q. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZTPLQT2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztplqt2 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztplqt2 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztplqt2 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDB, LDT, N, M, L */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), T( LDT, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTPLQT2 computes a LQ a factorization of a complex "triangular-pentagonal" */
/* > matrix C, which is composed of a triangular block A and pentagonal block B, */
/* > using the compact WY representation for Q. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The total number of rows of the matrix B. */
/* > M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix B, and the order of */
/* > the triangular matrix A. */
/* > N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > The number of rows of the lower trapezoidal part of B. */
/* > MIN(M,N) >= L >= 0. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,M) */
/* > On entry, the lower triangular M-by-M matrix A. */
/* > On exit, the elements on and below the diagonal of the array */
/* > contain the lower triangular matrix L. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,N) */
/* > On entry, the pentagonal M-by-N matrix B. The first N-L columns */
/* > are rectangular, and the last L columns are lower trapezoidal. */
/* > On exit, B contains the pentagonal matrix V. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is COMPLEX*16 array, dimension (LDT,M) */
/* > The N-by-N upper triangular factor T of the block reflector. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= max(1,M) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2017 */
/* > \ingroup doubleOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The input matrix C is a M-by-(M+N) matrix */
/* > */
/* > C = [ A ][ B ] */
/* > */
/* > */
/* > where A is an lower triangular M-by-M matrix, and B is M-by-N pentagonal */
/* > matrix consisting of a M-by-(N-L) rectangular matrix B1 left of a M-by-L */
/* > upper trapezoidal matrix B2: */
/* > */
/* > B = [ B1 ][ B2 ] */
/* > [ B1 ] <- M-by-(N-L) rectangular */
/* > [ B2 ] <- M-by-L lower trapezoidal. */
/* > */
/* > The lower trapezoidal matrix B2 consists of the first L columns of a */
/* > N-by-N lower triangular matrix, where 0 <= L <= MIN(M,N). If L=0, */
/* > B is rectangular M-by-N;
if M=L=N, B is lower triangular. */
/* > */
/* > The matrix W stores the elementary reflectors H(i) in the i-th row */
/* > above the diagonal (of A) in the M-by-(M+N) input matrix C */
/* > */
/* > C = [ A ][ B ] */
/* > [ A ] <- lower triangular M-by-M */
/* > [ B ] <- M-by-N pentagonal */
/* > */
/* > so that W can be represented as */
/* > */
/* > W = [ I ][ V ] */
/* > [ I ] <- identity, M-by-M */
/* > [ V ] <- M-by-N, same form as B. */
/* > */
/* > Thus, all of information needed for W is contained on exit in B, which */
/* > we call V above. Note that V has the same form as B;
that is, */
/* > */
/* > W = [ V1 ][ V2 ] */
/* > [ V1 ] <- M-by-(N-L) rectangular */
/* > [ V2 ] <- M-by-L lower trapezoidal. */
/* > */
/* > The rows of V represent the vectors which define the H(i)'s. */
/* > The (M+N)-by-(M+N) block reflector H is then given by */
/* > */
/* > H = I - W**T * T * W */
/* > */
/* > where W^H is the conjugate transpose of W and T is the upper triangular */
/* > factor of the block reflector. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int ztplqt2_(integer *m, integer *n, integer *l, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *t, integer *ldt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    integer i__, j, p, mp, np;
    doublecomplex alpha;
    extern /* Subroutine */
    int zgerc_(integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *), zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *), ztrmv_(char *, char *, char *, integer *, doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(char *, integer *), zlarfg_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *);
    /* -- LAPACK computational routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
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
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*l < 0 || *l > min(*m,*n))
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    else if (*ldb < max(1,*m))
    {
        *info = -7;
    }
    else if (*ldt < max(1,*m))
    {
        *info = -9;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZTPLQT2", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *m == 0)
    {
        return 0;
    }
    i__1 = *m;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        /* Generate elementary reflector H(I) to annihilate B(I,:) */
        p = *n - *l + min(*l,i__);
        i__2 = p + 1;
        zlarfg_(&i__2, &a[i__ + i__ * a_dim1], &b[i__ + b_dim1], ldb, &t[i__ * t_dim1 + 1]);
        i__2 = i__ * t_dim1 + 1;
        d_cnjg(&z__1, &t[i__ * t_dim1 + 1]);
        t[i__2].r = z__1.r;
        t[i__2].i = z__1.i; // , expr subst
        if (i__ < *m)
        {
            i__2 = p;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ + j * b_dim1;
                d_cnjg(&z__1, &b[i__ + j * b_dim1]);
                b[i__3].r = z__1.r;
                b[i__3].i = z__1.i; // , expr subst
            }
            /* W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)] */
            i__2 = *m - i__;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = *m + j * t_dim1;
                i__4 = i__ + j + i__ * a_dim1;
                t[i__3].r = a[i__4].r;
                t[i__3].i = a[i__4].i; // , expr subst
            }
            i__2 = *m - i__;
            zgemv_("N", &i__2, &p, &c_b2, &b[i__ + 1 + b_dim1], ldb, &b[i__ + b_dim1], ldb, &c_b2, &t[*m + t_dim1], ldt);
            /* C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H */
            i__2 = i__ * t_dim1 + 1;
            z__1.r = -t[i__2].r;
            z__1.i = -t[i__2].i; // , expr subst
            alpha.r = z__1.r;
            alpha.i = z__1.i; // , expr subst
            i__2 = *m - i__;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ + j + i__ * a_dim1;
                i__4 = i__ + j + i__ * a_dim1;
                i__5 = *m + j * t_dim1;
                z__2.r = alpha.r * t[i__5].r - alpha.i * t[i__5].i;
                z__2.i = alpha.r * t[i__5].i + alpha.i * t[i__5].r; // , expr subst
                z__1.r = a[i__4].r + z__2.r;
                z__1.i = a[i__4].i + z__2.i; // , expr subst
                a[i__3].r = z__1.r;
                a[i__3].i = z__1.i; // , expr subst
            }
            i__2 = *m - i__;
            z__1.r = alpha.r;
            z__1.i = alpha.i; // , expr subst
            zgerc_(&i__2, &p, &z__1, &t[*m + t_dim1], ldt, &b[i__ + b_dim1], ldb, &b[i__ + 1 + b_dim1], ldb);
            i__2 = p;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ + j * b_dim1;
                d_cnjg(&z__1, &b[i__ + j * b_dim1]);
                b[i__3].r = z__1.r;
                b[i__3].i = z__1.i; // , expr subst
            }
        }
    }
    i__1 = *m;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        /* T(I,1:I-1) := C(I:I-1,1:N)**H * (alpha * C(I,I:N)) */
        i__2 = i__ * t_dim1 + 1;
        z__1.r = -t[i__2].r;
        z__1.i = -t[i__2].i; // , expr subst
        alpha.r = z__1.r;
        alpha.i = z__1.i; // , expr subst
        i__2 = i__ - 1;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * t_dim1;
            t[i__3].r = 0.;
            t[i__3].i = 0.; // , expr subst
        }
        /* Computing MIN */
        i__2 = i__ - 1;
        p = min(i__2,*l);
        /* Computing MIN */
        i__2 = *n - *l + 1;
        np = min(i__2,*n);
        /* Computing MIN */
        i__2 = p + 1;
        mp = min(i__2,*m);
        i__2 = *n - *l + p;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * b_dim1;
            d_cnjg(&z__1, &b[i__ + j * b_dim1]);
            b[i__3].r = z__1.r;
            b[i__3].i = z__1.i; // , expr subst
        }
        /* Triangular part of B2 */
        i__2 = p;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * t_dim1;
            i__4 = i__ + (*n - *l + j) * b_dim1;
            z__1.r = alpha.r * b[i__4].r - alpha.i * b[i__4].i;
            z__1.i = alpha.r * b[i__4].i + alpha.i * b[i__4].r; // , expr subst
            t[i__3].r = z__1.r;
            t[i__3].i = z__1.i; // , expr subst
        }
        ztrmv_("L", "N", "N", &p, &b[np * b_dim1 + 1], ldb, &t[i__ + t_dim1], ldt);
        /* Rectangular part of B2 */
        i__2 = i__ - 1 - p;
        zgemv_("N", &i__2, l, &alpha, &b[mp + np * b_dim1], ldb, &b[i__ + np * b_dim1], ldb, &c_b1, &t[i__ + mp * t_dim1], ldt);
        /* B1 */
        i__2 = i__ - 1;
        i__3 = *n - *l;
        zgemv_("N", &i__2, &i__3, &alpha, &b[b_offset], ldb, &b[i__ + b_dim1], ldb, &c_b2, &t[i__ + t_dim1], ldt);
        /* T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1) */
        i__2 = i__ - 1;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * t_dim1;
            d_cnjg(&z__1, &t[i__ + j * t_dim1]);
            t[i__3].r = z__1.r;
            t[i__3].i = z__1.i; // , expr subst
        }
        i__2 = i__ - 1;
        ztrmv_("L", "C", "N", &i__2, &t[t_offset], ldt, &t[i__ + t_dim1], ldt);
        i__2 = i__ - 1;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * t_dim1;
            d_cnjg(&z__1, &t[i__ + j * t_dim1]);
            t[i__3].r = z__1.r;
            t[i__3].i = z__1.i; // , expr subst
        }
        i__2 = *n - *l + p;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * b_dim1;
            d_cnjg(&z__1, &b[i__ + j * b_dim1]);
            b[i__3].r = z__1.r;
            b[i__3].i = z__1.i; // , expr subst
        }
        /* T(I,I) = tau(I) */
        i__2 = i__ + i__ * t_dim1;
        i__3 = i__ * t_dim1 + 1;
        t[i__2].r = t[i__3].r;
        t[i__2].i = t[i__3].i; // , expr subst
        i__2 = i__ * t_dim1 + 1;
        t[i__2].r = 0.;
        t[i__2].i = 0.; // , expr subst
    }
    i__1 = *m;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = *m;
        for (j = i__ + 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * t_dim1;
            i__4 = j + i__ * t_dim1;
            t[i__3].r = t[i__4].r;
            t[i__3].i = t[i__4].i; // , expr subst
            i__3 = j + i__ * t_dim1;
            t[i__3].r = 0.;
            t[i__3].i = 0.; // , expr subst
        }
    }
    /* End of ZTPLQT2 */
    return 0;
}
/* ztplqt2_ */
