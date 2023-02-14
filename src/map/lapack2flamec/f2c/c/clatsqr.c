/* clatsqr.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__0 = 0;
/* > \brief \b CLATSQR */
/* Definition: */
/* =========== */
/* SUBROUTINE CLATSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, */
/* LWORK, INFO) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N, MB, NB, LDT, LWORK */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), T( LDT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLATSQR computes a blocked Tall-Skinny QR factorization of */
/* > a complex M-by-N matrix A for M >= N: */
/* > */
/* > A = Q * ( R ), */
/* > ( 0 ) */
/* > */
/* > where: */
/* > */
/* > Q is a M-by-M orthogonal matrix, stored on exit in an implicit */
/* > form in the elements below the diagonal of the array A and in */
/* > the elements of the array T;
*/
/* > */
/* > R is an upper-triangular N-by-N matrix, stored on exit in */
/* > the elements on and above the diagonal of the array A. */
/* > */
/* > 0 is a (M-N)-by-N zero matrix, and is not stored. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* > MB is INTEGER */
/* > The row block size to be used in the blocked QR. */
/* > MB > N. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The column block size to be used in the blocked QR. */
/* > N >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the elements on and above the diagonal */
/* > of the array contain the N-by-N upper triangular matrix R;
*/
/* > the elements below the diagonal represent Q by the columns */
/* > of blocked V (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is COMPLEX array, */
/* > dimension (LDT, N * Number_of_row_blocks) */
/* > where Number_of_row_blocks = CEIL((M-N)/(MB-N)) */
/* > The blocked upper triangular block reflectors stored in compact form */
/* > as a sequence of upper triangular blocks. */
/* > See Further Details below. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > (workspace) COMPLEX array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > The dimension of the array WORK. LWORK >= NB*N. */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
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
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > Tall-Skinny QR (TSQR) performs QR by a sequence of orthogonal transformations, */
/* > representing Q as a product of other orthogonal matrices */
/* > Q = Q(1) * Q(2) * . . . * Q(k) */
/* > where each Q(i) zeros out subdiagonal entries of a block of MB rows of A: */
/* > Q(1) zeros out the subdiagonal entries of rows 1:MB of A */
/* > Q(2) zeros out the bottom MB-N rows of rows [1:N,MB+1:2*MB-N] of A */
/* > Q(3) zeros out the bottom MB-N rows of rows [1:N,2*MB-N+1:3*MB-2*N] of A */
/* > . . . */
/* > */
/* > Q(1) is computed by GEQRT, which represents Q(1) by Householder vectors */
/* > stored under the diagonal of rows 1:MB of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,1:N). */
/* > For more information see Further Details in GEQRT. */
/* > */
/* > Q(i) for i>1 is computed by TPQRT, which represents Q(i) by Householder vectors */
/* > stored in rows [(i-1)*(MB-N)+N+1:i*(MB-N)+N] of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,(i-1)*N+1:i*N). */
/* > The last Q(k) may use fewer rows. */
/* > For more information see Further Details in TPQRT. */
/* > */
/* > For more details of the overall algorithm, see the description of */
/* > Sequential TSQR in Section 2.2 of [1]. */
/* > */
/* > [1] “Communication-Optimal Parallel and Sequential QR and LU Factorizations, */
/* > J. Demmel, L. Grigori, M. Hoemmen, J. Langou, */
/* > SIAM J. Sci. Comput, vol. 34, no. 1, 2012 */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int clatsqr_(integer *m, integer *n, integer *mb, integer * nb, complex *a, integer *lda, complex *t, integer *ldt, complex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"clatsqr inputs: m %lld, n %lld, mb %lld, nb %lld, lda %lld, ldt %lld, lwork %lld",*m, *n, *mb, *nb, *lda, *ldt, *lwork);
#else
    snprintf(buffer, 256,"clatsqr inputs: m %d, n %d, mb %d, nb %d, lda %d, ldt %d, lwork %d",*m, *n, *mb, *nb, *lda, *ldt, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, ii, kk, ctr;
    extern /* Subroutine */
    int xerbla_(char *, integer *), cgeqrt_( integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *), ctpqrt_(integer *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *);
    logical lquery;
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. -- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. EXTERNAL FUNCTIONS .. */
    /* .. EXTERNAL SUBROUTINES .. */
    /* .. INTRINSIC FUNCTIONS .. */
    /* .. */
    /* .. EXECUTABLE STATEMENTS .. */
    /* TEST THE INPUT ARGUMENTS */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0 || *m < *n)
    {
        *info = -2;
    }
    else if (*mb < 1)
    {
        *info = -3;
    }
    else if (*nb < 1 || *nb > *n && *n > 0)
    {
        *info = -4;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -6;
    }
    else if (*ldt < *nb)
    {
        *info = -8;
    }
    else if (*lwork < *n * *nb && ! lquery)
    {
        *info = -10;
    }
    if (*info == 0)
    {
        i__1 = *nb * *n;
        work[1].r = (real) i__1;
        work[1].i = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CLATSQR", &i__1);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (fla_min(*m,*n) == 0)
    {
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* The QR Decomposition */
    if (*mb <= *n || *mb >= *m)
    {
        cgeqrt_(m, n, nb, &a[a_offset], lda, &t[t_offset], ldt, &work[1], info);
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    kk = (*m - *n) % (*mb - *n);
    ii = *m - kk + 1;
    /* Compute the QR factorization of the first block A(1:MB,1:N) */
    cgeqrt_(mb, n, nb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &work[1], info) ;
    ctr = 1;
    i__1 = ii - *mb + *n;
    i__2 = *mb - *n;
    for (i__ = *mb + 1;
            i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
            i__ += i__2)
    {
        /* Compute the QR factorization of the current block A(I:I+MB-N,1:N) */
        i__3 = *mb - *n;
        ctpqrt_(&i__3, n, &c__0, nb, &a[a_dim1 + 1], lda, &a[i__ + a_dim1], lda, &t[(ctr * *n + 1) * t_dim1 + 1], ldt, &work[1], info);
        ++ctr;
    }
    /* Compute the QR factorization of the last block A(II:M,1:N) */
    if (ii <= *m)
    {
        ctpqrt_(&kk, n, &c__0, nb, &a[a_dim1 + 1], lda, &a[ii + a_dim1], lda, &t[(ctr * *n + 1) * t_dim1 + 1], ldt, &work[1], info);
    }
    i__2 = *n * *nb;
    work[1].r = (real) i__2;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CLATSQR */
}
/* clatsqr_ */
