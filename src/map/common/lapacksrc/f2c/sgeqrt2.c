/* ../netlib/sgeqrt2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b5 = 1.f;
static real c_b7 = 0.f;
/* > \brief \b SGEQRT2 computes a QR factorization of a general real or complex matrix using the compact WY re presentation of Q. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGEQRT2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrt2 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrt2 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrt2 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGEQRT2( M, N, A, LDA, T, LDT, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDT, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), T( LDT, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEQRT2 computes a QR factorization of a real M-by-N matrix A, */
/* > using the compact WY representation of Q. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the real M-by-N matrix A. On exit, the elements on and */
/* > above the diagonal contain the N-by-N upper triangular matrix R;
the */
/* > elements below the diagonal are the columns of V. See below for */
/* > further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is REAL array, dimension (LDT,N) */
/* > The N-by-N upper triangular factor of the block reflector. */
/* > The elements on and above the diagonal contain the block */
/* > reflector T;
the elements below the diagonal are not used. */
/* > See below for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= max(1,N). */
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
/* > \ingroup realGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix V stores the elementary reflectors H(i) in the i-th column */
/* > below the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* > V = ( 1 ) */
/* > ( v1 1 ) */
/* > ( v1 v2 1 ) */
/* > ( v1 v2 v3 ) */
/* > ( v1 v2 v3 ) */
/* > */
/* > where the vi's represent the vectors which define H(i), which are returned */
/* > in the matrix A. The 1's along the diagonal of V are not stored in A. The */
/* > block reflector H is then given by */
/* > */
/* > H = I - V * T * V**T */
/* > */
/* > where V**T is the transpose of V. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int sgeqrt2_(integer *m, integer *n, real *a, integer *lda, real *t, integer *ldt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, k;
    real aii;
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
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
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
    else if (*lda < max(1,*m))
    {
        *info = -4;
    }
    else if (*ldt < max(1,*n))
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGEQRT2", &i__1);
        return 0;
    }
    k = min(*m,*n);
    i__1 = k;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        /* Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1) */
        i__2 = *m - i__ + 1;
        /* Computing MIN */
        i__3 = i__ + 1;
        slarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3,*m) + i__ * a_dim1] , &c__1, &t[i__ + t_dim1]);
        if (i__ < *n)
        {
            /* Apply H(i) to A(I:M,I+1:N) from the left */
            aii = a[i__ + i__ * a_dim1];
            a[i__ + i__ * a_dim1] = 1.f;
            /* W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)] */
            i__2 = *m - i__ + 1;
            i__3 = *n - i__;
            sgemv_("T", &i__2, &i__3, &c_b5, &a[i__ + (i__ + 1) * a_dim1], lda, &a[i__ + i__ * a_dim1], &c__1, &c_b7, &t[*n * t_dim1 + 1], &c__1);
            /* A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H */
            alpha = -t[i__ + t_dim1];
            i__2 = *m - i__ + 1;
            i__3 = *n - i__;
            sger_(&i__2, &i__3, &alpha, &a[i__ + i__ * a_dim1], &c__1, &t[*n * t_dim1 + 1], &c__1, &a[i__ + (i__ + 1) * a_dim1], lda);
            a[i__ + i__ * a_dim1] = aii;
        }
    }
    i__1 = *n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        aii = a[i__ + i__ * a_dim1];
        a[i__ + i__ * a_dim1] = 1.f;
        /* T(1:I-1,I) := alpha * A(I:M,1:I-1)**T * A(I:M,I) */
        alpha = -t[i__ + t_dim1];
        i__2 = *m - i__ + 1;
        i__3 = i__ - 1;
        sgemv_("T", &i__2, &i__3, &alpha, &a[i__ + a_dim1], lda, &a[i__ + i__ * a_dim1], &c__1, &c_b7, &t[i__ * t_dim1 + 1], &c__1);
        a[i__ + i__ * a_dim1] = aii;
        /* T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I) */
        i__2 = i__ - 1;
        strmv_("U", "N", "N", &i__2, &t[t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1);
        /* T(I,I) = tau(I) */
        t[i__ + i__ * t_dim1] = t[i__ + t_dim1];
        t[i__ + t_dim1] = 0.f;
    }
    /* End of SGEQRT2 */
    return 0;
}
/* sgeqrt2_ */
