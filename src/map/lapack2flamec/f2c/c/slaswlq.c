/* slaswlq.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__0 = 0;
/* > \brief \b SLASWLQ */
/* Definition: */
/* =========== */
/* SUBROUTINE SLASWLQ( M, N, MB, NB, A, LDA, T, LDT, WORK, */
/* LWORK, INFO) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N, MB, NB, LDT, LWORK */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), T( LDT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASWLQ computes a blocked Tall-Skinny LQ factorization of */
/* > a real M-by-N matrix A for M <= N: */
/* > */
/* > A = ( L 0 ) * Q, */
/* > */
/* > where: */
/* > */
/* > Q is a n-by-N orthogonal matrix, stored on exit in an implicit */
/* > form in the elements above the diagonal of the array A and in */
/* > the elements of the array T;
*/
/* > L is a lower-triangular M-by-M matrix stored on exit in */
/* > the elements on and below the diagonal of the array A. */
/* > 0 is a M-by-(N-M) zero matrix, if M < N, and is not stored. */
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
/* > The number of columns of the matrix A. N >= M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* > MB is INTEGER */
/* > The row block size to be used in the blocked QR. */
/* > M >= MB >= 1 */
/* > \endverbatim */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The column block size to be used in the blocked QR. */
/* > NB > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the elements on and below the diagonal */
/* > of the array contain the N-by-N lower triangular matrix L;
*/
/* > the elements above the diagonal represent Q by the rows */
/* > of blocked V (see Further Details). */
/* > */
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
/* > T is REAL array, */
/* > dimension (LDT, N * Number_of_row_blocks) */
/* > where Number_of_row_blocks = CEIL((N-M)/(NB-M)) */
/* > The blocked upper triangular block reflectors stored in compact form */
/* > as a sequence of upper triangular blocks. */
/* > See Further Details below. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= MB. */
/* > \endverbatim */
/* > */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > (workspace) REAL array, dimension (MAX(1,LWORK)) */
/* > */
/* > \endverbatim */
/* > \param[in] LWORK */
/* > \verbatim */
/* > The dimension of the array WORK. LWORK >= MB * M. */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > */
/* > \endverbatim */
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
/* > Short-Wide LQ (SWLQ) performs LQ by a sequence of orthogonal transformations, */
/* > representing Q as a product of other orthogonal matrices */
/* > Q = Q(1) * Q(2) * . . . * Q(k) */
/* > where each Q(i) zeros out upper diagonal entries of a block of NB rows of A: */
/* > Q(1) zeros out the upper diagonal entries of rows 1:NB of A */
/* > Q(2) zeros out the bottom MB-N rows of rows [1:M,NB+1:2*NB-M] of A */
/* > Q(3) zeros out the bottom MB-N rows of rows [1:M,2*NB-M+1:3*NB-2*M] of A */
/* > . . . */
/* > */
/* > Q(1) is computed by GELQT, which represents Q(1) by Householder vectors */
/* > stored under the diagonal of rows 1:MB of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,1:N). */
/* > For more information see Further Details in GELQT. */
/* > */
/* > Q(i) for i>1 is computed by TPLQT, which represents Q(i) by Householder vectors */
/* > stored in columns [(i-1)*(NB-M)+M+1:i*(NB-M)+M] of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,(i-1)*M+1:i*M). */
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
int slaswlq_(integer *m, integer *n, integer *mb, integer * nb, real *a, integer *lda, real *t, integer *ldt, real *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slaswlq inputs: m %" FLA_IS ", n %" FLA_IS ", mb %" FLA_IS ", nb %" FLA_IS ", lda %" FLA_IS ", ldt %" FLA_IS ", lwork %" FLA_IS "",*m, *n, *mb, *nb, *lda, *ldt, *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, ii, kk, ctr;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len), sgelqt_( integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *);
    logical lquery;
    extern /* Subroutine */
    int stplqt_(integer *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *);
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
    else if (*n < 0 || *n < *m)
    {
        *info = -2;
    }
    else if (*mb < 1 || *mb > *m && *m > 0)
    {
        *info = -3;
    }
    else if (*nb <= 0)
    {
        *info = -4;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -6;
    }
    else if (*ldt < *mb)
    {
        *info = -8;
    }
    else if (*lwork < *m * *mb && ! lquery)
    {
        *info = -10;
    }
    if (*info == 0)
    {
        work[1] = (real) (*mb * *m);
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLASWLQ", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    if (fla_min(*m,*n) == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* The LQ Decomposition */
    if (*m >= *n || *nb <= *m || *nb >= *n)
    {
        sgelqt_(m, n, mb, &a[a_offset], lda, &t[t_offset], ldt, &work[1], info);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    kk = (*n - *m) % (*nb - *m);
    ii = *n - kk + 1;
    /* Compute the LQ factorization of the first block A(1:M,1:NB) */
    sgelqt_(m, nb, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &work[1], info) ;
    ctr = 1;
    i__1 = ii - *nb + *m;
    i__2 = *nb - *m;
    for (i__ = *nb + 1;
            i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
            i__ += i__2)
    {
        /* Compute the QR factorization of the current block A(1:M,I:I+NB-M) */
        i__3 = *nb - *m;
        stplqt_(m, &i__3, &c__0, mb, &a[a_dim1 + 1], lda, &a[i__ * a_dim1 + 1], lda, &t[(ctr * *m + 1) * t_dim1 + 1], ldt, &work[1], info);
        ++ctr;
    }
    /* Compute the QR factorization of the last block A(1:M,II:N) */
    if (ii <= *n)
    {
        stplqt_(m, &kk, &c__0, mb, &a[a_dim1 + 1], lda, &a[ii * a_dim1 + 1], lda, &t[(ctr * *m + 1) * t_dim1 + 1], ldt, &work[1], info);
    }
    work[1] = (real) (*m * *mb);
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of SLASWLQ */
}
/* slaswlq_ */
