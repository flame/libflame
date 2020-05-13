/* ../netlib/ctpqrt2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    1.f,0.f
}
;
static complex c_b2 =
{
    0.f,0.f
}
;
static integer c__1 = 1;
/* > \brief \b CTPQRT2 computes a QR factorization of a real or complex "triangular-pentagonal" matrix, which is composed of a triangular block and a pentagonal block, using the compact WY representation for Q. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CTPQRT2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctpqrt2 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctpqrt2 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctpqrt2 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CTPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDB, LDT, N, M, L */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), B( LDB, * ), T( LDT, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPQRT2 computes a QR factorization of a complex "triangular-pentagonal" */
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
/* > The number of rows of the upper trapezoidal part of B. */
/* > MIN(M,N) >= L >= 0. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the upper triangular N-by-N matrix A. */
/* > On exit, the elements on and above the diagonal of the array */
/* > contain the upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,N) */
/* > On entry, the pentagonal M-by-N matrix B. The first M-L rows */
/* > are rectangular, and the last L rows are upper trapezoidal. */
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
/* > T is COMPLEX array, dimension (LDT,N) */
/* > The N-by-N upper triangular factor T of the block reflector. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= max(1,N) */
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
/* > \date September 2012 */
/* > \ingroup complexOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The input matrix C is a (N+M)-by-N matrix */
/* > */
/* > C = [ A ] */
/* > [ B ] */
/* > */
/* > where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal */
/* > matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N */
/* > upper trapezoidal matrix B2: */
/* > */
/* > B = [ B1 ] <- (M-L)-by-N rectangular */
/* > [ B2 ] <- L-by-N upper trapezoidal. */
/* > */
/* > The upper trapezoidal matrix B2 consists of the first L rows of a */
/* > N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N). If L=0, */
/* > B is rectangular M-by-N;
if M=L=N, B is upper triangular. */
/* > */
/* > The matrix W stores the elementary reflectors H(i) in the i-th column */
/* > below the diagonal (of A) in the (N+M)-by-N input matrix C */
/* > */
/* > C = [ A ] <- upper triangular N-by-N */
/* > [ B ] <- M-by-N pentagonal */
/* > */
/* > so that W can be represented as */
/* > */
/* > W = [ I ] <- identity, N-by-N */
/* > [ V ] <- M-by-N, same form as B. */
/* > */
/* > Thus, all of information needed for W is contained on exit in B, which */
/* > we call V above. Note that V has the same form as B;
that is, */
/* > */
/* > V = [ V1 ] <- (M-L)-by-N rectangular */
/* > [ V2 ] <- L-by-N upper trapezoidal. */
/* > */
/* > The columns of V represent the vectors which define the H(i)'s. */
/* > The (M+N)-by-(M+N) block reflector H is then given by */
/* > */
/* > H = I - W * T * W**H */
/* > */
/* > where W**H is the conjugate transpose of W and T is the upper triangular */
/* > factor of the block reflector. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int ctpqrt2_(integer *m, integer *n, integer *l, complex *a, integer *lda, complex *b, integer *ldb, complex *t, integer *ldt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    integer i__, j, p, mp, np;
    extern /* Subroutine */
    int cgerc_(integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, integer *);
    complex alpha;
    extern /* Subroutine */
    int cgemv_(char *, integer *, integer *, complex * , complex *, integer *, complex *, integer *, complex *, complex * , integer *), ctrmv_(char *, char *, char *, integer *, complex *, integer *, complex *, integer *), clarfg_(integer *, complex *, complex *, integer *, complex *), xerbla_(char *, integer *);
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
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*ldb < max(1,*m))
    {
        *info = -7;
    }
    else if (*ldt < max(1,*n))
    {
        *info = -9;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CTPQRT2", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *m == 0)
    {
        return 0;
    }
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        /* Generate elementary reflector H(I) to annihilate B(:,I) */
        p = *m - *l + min(*l,i__);
        i__2 = p + 1;
        clarfg_(&i__2, &a[i__ + i__ * a_dim1], &b[i__ * b_dim1 + 1], &c__1, & t[i__ + t_dim1]);
        if (i__ < *n)
        {
            /* W(1:N-I) := C(I:M,I+1:N)**H * C(I:M,I) [use W = T(:,N)] */
            i__2 = *n - i__;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = j + *n * t_dim1;
                r_cnjg(&q__1, &a[i__ + (i__ + j) * a_dim1]);
                t[i__3].r = q__1.r;
                t[i__3].i = q__1.i; // , expr subst
            }
            i__2 = *n - i__;
            cgemv_("C", &p, &i__2, &c_b1, &b[(i__ + 1) * b_dim1 + 1], ldb, &b[ i__ * b_dim1 + 1], &c__1, &c_b1, &t[*n * t_dim1 + 1], & c__1);
            /* C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)**H */
            r_cnjg(&q__2, &t[i__ + t_dim1]);
            q__1.r = -q__2.r;
            q__1.i = -q__2.i; // , expr subst
            alpha.r = q__1.r;
            alpha.i = q__1.i; // , expr subst
            i__2 = *n - i__;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ + (i__ + j) * a_dim1;
                i__4 = i__ + (i__ + j) * a_dim1;
                r_cnjg(&q__3, &t[j + *n * t_dim1]);
                q__2.r = alpha.r * q__3.r - alpha.i * q__3.i;
                q__2.i = alpha.r * q__3.i + alpha.i * q__3.r; // , expr subst
                q__1.r = a[i__4].r + q__2.r;
                q__1.i = a[i__4].i + q__2.i; // , expr subst
                a[i__3].r = q__1.r;
                a[i__3].i = q__1.i; // , expr subst
            }
            i__2 = *n - i__;
            cgerc_(&p, &i__2, &alpha, &b[i__ * b_dim1 + 1], &c__1, &t[*n * t_dim1 + 1], &c__1, &b[(i__ + 1) * b_dim1 + 1], ldb);
        }
    }
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        /* T(1:I-1,I) := C(I:M,1:I-1)**H * (alpha * C(I:M,I)) */
        i__2 = i__ + t_dim1;
        q__1.r = -t[i__2].r;
        q__1.i = -t[i__2].i; // , expr subst
        alpha.r = q__1.r;
        alpha.i = q__1.i; // , expr subst
        i__2 = i__ - 1;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = j + i__ * t_dim1;
            t[i__3].r = 0.f;
            t[i__3].i = 0.f; // , expr subst
        }
        /* Computing MIN */
        i__2 = i__ - 1;
        p = min(i__2,*l);
        /* Computing MIN */
        i__2 = *m - *l + 1;
        mp = min(i__2,*m);
        /* Computing MIN */
        i__2 = p + 1;
        np = min(i__2,*n);
        /* Triangular part of B2 */
        i__2 = p;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = j + i__ * t_dim1;
            i__4 = *m - *l + j + i__ * b_dim1;
            q__1.r = alpha.r * b[i__4].r - alpha.i * b[i__4].i;
            q__1.i = alpha.r * b[i__4].i + alpha.i * b[i__4].r; // , expr subst
            t[i__3].r = q__1.r;
            t[i__3].i = q__1.i; // , expr subst
        }
        ctrmv_("U", "C", "N", &p, &b[mp + b_dim1], ldb, &t[i__ * t_dim1 + 1], &c__1);
        /* Rectangular part of B2 */
        i__2 = i__ - 1 - p;
        cgemv_("C", l, &i__2, &alpha, &b[mp + np * b_dim1], ldb, &b[mp + i__ * b_dim1], &c__1, &c_b2, &t[np + i__ * t_dim1], &c__1);
        /* B1 */
        i__2 = *m - *l;
        i__3 = i__ - 1;
        cgemv_("C", &i__2, &i__3, &alpha, &b[b_offset], ldb, &b[i__ * b_dim1 + 1], &c__1, &c_b1, &t[i__ * t_dim1 + 1], &c__1);
        /* T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I) */
        i__2 = i__ - 1;
        ctrmv_("U", "N", "N", &i__2, &t[t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1);
        /* T(I,I) = tau(I) */
        i__2 = i__ + i__ * t_dim1;
        i__3 = i__ + t_dim1;
        t[i__2].r = t[i__3].r;
        t[i__2].i = t[i__3].i; // , expr subst
        i__2 = i__ + t_dim1;
        t[i__2].r = 0.f;
        t[i__2].i = 0.f; // , expr subst
    }
    /* End of CTPQRT2 */
    return 0;
}
/* ctpqrt2_ */
