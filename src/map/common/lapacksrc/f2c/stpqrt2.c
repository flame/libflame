/* ../netlib/stpqrt2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b5 = 1.f;
static real c_b17 = 0.f;
/* > \brief \b STPQRT2 computes a QR factorization of a real or complex "triangular-pentagonal" matrix, which is composed of a triangular block and a pentagonal block, using the compact WY representation for Q. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STPQRT2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpqrt2 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpqrt2 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpqrt2 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDB, LDT, N, M, L */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), T( LDT, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPQRT2 computes a QR factorization of a real "triangular-pentagonal" */
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
/* > A is REAL array, dimension (LDA,N) */
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
/* > B is REAL array, dimension (LDB,N) */
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
/* > T is REAL array, dimension (LDT,N) */
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
/* > \ingroup realOTHERcomputational */
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
/* > H = I - W * T * W^H */
/* > */
/* > where W^H is the conjugate transpose of W and T is the upper triangular */
/* > factor of the block reflector. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int stpqrt2_(integer *m, integer *n, integer *l, real *a, integer *lda, real *b, integer *ldb, real *t, integer *ldt, integer * info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, j, p, mp, np;
    extern /* Subroutine */
    int sger_(integer *, integer *, real *, real *, integer *, real *, integer *, real *, integer *);
    real alpha;
    extern /* Subroutine */
    int sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), strmv_(char *, char *, char *, integer *, real *, integer *, real *, integer *), xerbla_( char *, integer *), slarfg_(integer *, real *, real *, integer *, real *);
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
        xerbla_("STPQRT2", &i__1);
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
        slarfg_(&i__2, &a[i__ + i__ * a_dim1], &b[i__ * b_dim1 + 1], &c__1, & t[i__ + t_dim1]);
        if (i__ < *n)
        {
            /* W(1:N-I) := C(I:M,I+1:N)^H * C(I:M,I) [use W = T(:,N)] */
            i__2 = *n - i__;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                t[j + *n * t_dim1] = a[i__ + (i__ + j) * a_dim1];
            }
            i__2 = *n - i__;
            sgemv_("T", &p, &i__2, &c_b5, &b[(i__ + 1) * b_dim1 + 1], ldb, &b[ i__ * b_dim1 + 1], &c__1, &c_b5, &t[*n * t_dim1 + 1], & c__1);
            /* C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)^H */
            alpha = -t[i__ + t_dim1];
            i__2 = *n - i__;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                a[i__ + (i__ + j) * a_dim1] += alpha * t[j + *n * t_dim1];
            }
            i__2 = *n - i__;
            sger_(&p, &i__2, &alpha, &b[i__ * b_dim1 + 1], &c__1, &t[*n * t_dim1 + 1], &c__1, &b[(i__ + 1) * b_dim1 + 1], ldb);
        }
    }
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        /* T(1:I-1,I) := C(I:M,1:I-1)^H * (alpha * C(I:M,I)) */
        alpha = -t[i__ + t_dim1];
        i__2 = i__ - 1;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            t[j + i__ * t_dim1] = 0.f;
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
            t[j + i__ * t_dim1] = alpha * b[*m - *l + j + i__ * b_dim1];
        }
        strmv_("U", "T", "N", &p, &b[mp + b_dim1], ldb, &t[i__ * t_dim1 + 1], &c__1);
        /* Rectangular part of B2 */
        i__2 = i__ - 1 - p;
        sgemv_("T", l, &i__2, &alpha, &b[mp + np * b_dim1], ldb, &b[mp + i__ * b_dim1], &c__1, &c_b17, &t[np + i__ * t_dim1], &c__1);
        /* B1 */
        i__2 = *m - *l;
        i__3 = i__ - 1;
        sgemv_("T", &i__2, &i__3, &alpha, &b[b_offset], ldb, &b[i__ * b_dim1 + 1], &c__1, &c_b5, &t[i__ * t_dim1 + 1], &c__1);
        /* T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I) */
        i__2 = i__ - 1;
        strmv_("U", "N", "N", &i__2, &t[t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1);
        /* T(I,I) = tau(I) */
        t[i__ + i__ * t_dim1] = t[i__ + t_dim1];
        t[i__ + t_dim1] = 0.f;
    }
    /* End of STPQRT2 */
    return 0;
}
/* stpqrt2_ */
