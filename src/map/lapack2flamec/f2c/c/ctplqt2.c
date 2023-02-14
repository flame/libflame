/* ../netlib/v3.9.0/ctplqt2.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    0.f,0.f
}
;
static complex c_b2 =
{
    1.f,0.f
}
;
/* > \brief \b CTPLQT2 */
/* Definition: */
/* =========== */
/* SUBROUTINE CTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO ) */
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
/* > CTPLQT2 computes a LQ a factorization of a complex "triangular-pentagonal" */
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
/* > A is COMPLEX array, dimension (LDA,M) */
/* > On entry, the lower triangular M-by-M matrix A. */
/* > On exit, the elements on and below the diagonal of the array */
/* > contain the lower triangular matrix L. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,N) */
/* > On entry, the pentagonal M-by-N matrix B. The first N-L columns */
/* > are rectangular, and the last L columns are lower trapezoidal. */
/* > On exit, B contains the pentagonal matrix V. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is COMPLEX array, dimension (LDT,M) */
/* > The N-by-N upper triangular factor T of the block reflector. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= fla_max(1,M) */
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
int ctplqt2_(integer *m, integer *n, integer *l, complex *a, integer *lda, complex *b, integer *ldb, complex *t, integer *ldt, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"ctplqt2 inputs: m %lld, n %lld, l %lld, lda %lld, ldb %lld, ldt %lld",*m, *n, *l, *lda, *ldb, *ldt);
#else
    snprintf(buffer, 256,"ctplqt2 inputs: m %d, n %d, l %d, lda %d, ldb %d, ldt %d",*m, *n, *l, *lda, *ldb, *ldt);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    integer i__, j, p, mp, np;
    extern /* Subroutine */
    int cgerc_(integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, integer *);
    complex alpha;
    extern /* Subroutine */
    int cgemv_(char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *), ctrmv_(char *, char *, char *, integer *, complex *, integer *, complex *, integer *), clarfg_(integer *, complex *, complex *, integer *, complex *), xerbla_(char *, integer *);
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
    else if (*l < 0 || *l > fla_min(*m,*n))
    {
        *info = -3;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -5;
    }
    else if (*ldb < fla_max(1,*m))
    {
        *info = -7;
    }
    else if (*ldt < fla_max(1,*m))
    {
        *info = -9;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CTPLQT2", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *m == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    i__1 = *m;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        /* Generate elementary reflector H(I) to annihilate B(I,:) */
        p = *n - *l + fla_min(*l,i__);
        i__2 = p + 1;
        clarfg_(&i__2, &a[i__ + i__ * a_dim1], &b[i__ + b_dim1], ldb, &t[i__ * t_dim1 + 1]);
        i__2 = i__ * t_dim1 + 1;
        r_cnjg(&q__1, &t[i__ * t_dim1 + 1]);
        t[i__2].r = q__1.r;
        t[i__2].i = q__1.i; // , expr subst
        if (i__ < *m)
        {
            i__2 = p;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ + j * b_dim1;
                r_cnjg(&q__1, &b[i__ + j * b_dim1]);
                b[i__3].r = q__1.r;
                b[i__3].i = q__1.i; // , expr subst
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
            cgemv_("N", &i__2, &p, &c_b2, &b[i__ + 1 + b_dim1], ldb, &b[i__ + b_dim1], ldb, &c_b2, &t[*m + t_dim1], ldt);
            /* C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H */
            i__2 = i__ * t_dim1 + 1;
            q__1.r = -t[i__2].r;
            q__1.i = -t[i__2].i; // , expr subst
            alpha.r = q__1.r;
            alpha.i = q__1.i; // , expr subst
            i__2 = *m - i__;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ + j + i__ * a_dim1;
                i__4 = i__ + j + i__ * a_dim1;
                i__5 = *m + j * t_dim1;
                q__2.r = alpha.r * t[i__5].r - alpha.i * t[i__5].i;
                q__2.i = alpha.r * t[i__5].i + alpha.i * t[i__5].r; // , expr subst
                q__1.r = a[i__4].r + q__2.r;
                q__1.i = a[i__4].i + q__2.i; // , expr subst
                a[i__3].r = q__1.r;
                a[i__3].i = q__1.i; // , expr subst
            }
            i__2 = *m - i__;
            q__1.r = alpha.r;
            q__1.i = alpha.i; // , expr subst
            cgerc_(&i__2, &p, &q__1, &t[*m + t_dim1], ldt, &b[i__ + b_dim1], ldb, &b[i__ + 1 + b_dim1], ldb);
            i__2 = p;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ + j * b_dim1;
                r_cnjg(&q__1, &b[i__ + j * b_dim1]);
                b[i__3].r = q__1.r;
                b[i__3].i = q__1.i; // , expr subst
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
        q__1.r = -t[i__2].r;
        q__1.i = -t[i__2].i; // , expr subst
        alpha.r = q__1.r;
        alpha.i = q__1.i; // , expr subst
        i__2 = i__ - 1;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * t_dim1;
            t[i__3].r = 0.f;
            t[i__3].i = 0.f; // , expr subst
        }
        /* Computing MIN */
        i__2 = i__ - 1;
        p = fla_min(i__2,*l);
        /* Computing MIN */
        i__2 = *n - *l + 1;
        np = fla_min(i__2,*n);
        /* Computing MIN */
        i__2 = p + 1;
        mp = fla_min(i__2,*m);
        i__2 = *n - *l + p;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * b_dim1;
            r_cnjg(&q__1, &b[i__ + j * b_dim1]);
            b[i__3].r = q__1.r;
            b[i__3].i = q__1.i; // , expr subst
        }
        /* Triangular part of B2 */
        i__2 = p;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * t_dim1;
            i__4 = i__ + (*n - *l + j) * b_dim1;
            q__1.r = alpha.r * b[i__4].r - alpha.i * b[i__4].i;
            q__1.i = alpha.r * b[i__4].i + alpha.i * b[i__4].r; // , expr subst
            t[i__3].r = q__1.r;
            t[i__3].i = q__1.i; // , expr subst
        }
        ctrmv_("L", "N", "N", &p, &b[np * b_dim1 + 1], ldb, &t[i__ + t_dim1], ldt);
        /* Rectangular part of B2 */
        i__2 = i__ - 1 - p;
        cgemv_("N", &i__2, l, &alpha, &b[mp + np * b_dim1], ldb, &b[i__ + np * b_dim1], ldb, &c_b1, &t[i__ + mp * t_dim1], ldt);
        /* B1 */
        i__2 = i__ - 1;
        i__3 = *n - *l;
        cgemv_("N", &i__2, &i__3, &alpha, &b[b_offset], ldb, &b[i__ + b_dim1], ldb, &c_b2, &t[i__ + t_dim1], ldt);
        /* T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1) */
        i__2 = i__ - 1;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * t_dim1;
            r_cnjg(&q__1, &t[i__ + j * t_dim1]);
            t[i__3].r = q__1.r;
            t[i__3].i = q__1.i; // , expr subst
        }
        i__2 = i__ - 1;
        ctrmv_("L", "C", "N", &i__2, &t[t_offset], ldt, &t[i__ + t_dim1], ldt);
        i__2 = i__ - 1;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * t_dim1;
            r_cnjg(&q__1, &t[i__ + j * t_dim1]);
            t[i__3].r = q__1.r;
            t[i__3].i = q__1.i; // , expr subst
        }
        i__2 = *n - *l + p;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = i__ + j * b_dim1;
            r_cnjg(&q__1, &b[i__ + j * b_dim1]);
            b[i__3].r = q__1.r;
            b[i__3].i = q__1.i; // , expr subst
        }
        /* T(I,I) = tau(I) */
        i__2 = i__ + i__ * t_dim1;
        i__3 = i__ * t_dim1 + 1;
        t[i__2].r = t[i__3].r;
        t[i__2].i = t[i__3].i; // , expr subst
        i__2 = i__ * t_dim1 + 1;
        t[i__2].r = 0.f;
        t[i__2].i = 0.f; // , expr subst
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
            t[i__3].r = 0.f;
            t[i__3].i = 0.f; // , expr subst
        }
    }
    /* End of CTPLQT2 */
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
}
/* ctplqt2_ */
