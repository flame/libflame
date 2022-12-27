/* zlamswlq.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__0 = 0;
/* > \brief \b ZLAMSWLQ */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, */
/* $ LDT, C, LDC, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
/* INTEGER INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), WORK( * ), C(LDC, * ), */
/* $ T( LDT, * ) */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAMSWLQ overwrites the general complex M-by-N matrix C with */
/* > */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'C': Q**H * C C * Q**H */
/* > where Q is a complex unitary matrix defined as the product of blocked */
/* > elementary reflectors computed by short wide LQ */
/* > factorization (ZLASWLQ) */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': apply Q or Q**H from the Left;
*/
/* > = 'R': apply Q or Q**H from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': No transpose, apply Q;
*/
/* > = 'C': Conjugate Transpose, apply Q**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. M >=0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines */
/* > the matrix Q. */
/* > M >= K >= 0;
*/
/* > */
/* > \endverbatim */
/* > \param[in] MB */
/* > \verbatim */
/* > MB is INTEGER */
/* > The row block size to be used in the blocked LQ. */
/* > M >= MB >= 1 */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The column block size to be used in the blocked LQ. */
/* > NB > M. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension */
/* > (LDA,M) if SIDE = 'L', */
/* > (LDA,N) if SIDE = 'R' */
/* > The i-th row must contain the vector which defines the blocked */
/* > elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* > ZLASWLQ in the first k rows of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= MAX(1,K). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* > T is COMPLEX*16 array, dimension */
/* > ( M * Number of blocks(CEIL(N-K/NB-K)), */
/* > The blocked upper triangular block reflectors stored in compact form */
/* > as a sequence of upper triangular blocks. See below */
/* > for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= MB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX*16 array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If SIDE = 'L', LWORK >= fla_max(1,NB) * MB;
*/
/* > if SIDE = 'R', LWORK >= fla_max(1,M) * MB. */
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
/* > Short-Wide LQ (SWLQ) performs LQ by a sequence of unitary transformations, */
/* > representing Q as a product of other unitary matrices */
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
/* > For more information see Further Details in TPLQT. */
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
int zlamswlq_(char *side, char *trans, integer *m, integer * n, integer *k, integer *mb, integer *nb, doublecomplex *a, integer * lda, doublecomplex *t, integer *ldt, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlamswlq inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", mb %" FLA_IS ", nb %" FLA_IS ", lda %" FLA_IS ", ldt %" FLA_IS ", ldc %" FLA_IS ", lwork %" FLA_IS "",*side, *trans, *m, *n, *k, *mb, *nb, *lda, *ldt, *ldc, *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, ii, kk, lw, ctr;
    logical left, tran;
    extern logical lsame_(char *, char *);
    logical right;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    logical notran, lquery;
    extern /* Subroutine */
    int zgemlqt_(char *, char *, integer *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *), ztpmlqt_(char *, char *, integer *, integer *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
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
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    lquery = *lwork < 0;
    notran = lsame_(trans, "N");
    tran = lsame_(trans, "C");
    left = lsame_(side, "L");
    right = lsame_(side, "R");
    if (left)
    {
        lw = *n * *mb;
    }
    else
    {
        lw = *m * *mb;
    }
    *info = 0;
    if (! left && ! right)
    {
        *info = -1;
    }
    else if (! tran && ! notran)
    {
        *info = -2;
    }
    else if (*k < 0)
    {
        *info = -5;
    }
    else if (*m < *k)
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*k < *mb || *mb < 1)
    {
        *info = -6;
    }
    else if (*lda < fla_max(1,*k))
    {
        *info = -9;
    }
    else if (*ldt < fla_max(1,*mb))
    {
        *info = -11;
    }
    else if (*ldc < fla_max(1,*m))
    {
        *info = -13;
    }
    else if (*lwork < fla_max(1,lw) && ! lquery)
    {
        *info = -15;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZLAMSWLQ", &i__1);
        work[1].r = (doublereal) lw;
        work[1].i = 0.; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (lquery)
    {
        work[1].r = (doublereal) lw;
        work[1].i = 0.; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    /* Computing MIN */
    i__1 = fla_min(*m,*n);
    if (fla_min(i__1,*k) == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Computing MAX */
    i__1 = fla_max(*m,*n);
    if (*nb <= *k || *nb >= fla_max(i__1,*k))
    {
        zgemlqt_(side, trans, m, n, k, mb, &a[a_offset], lda, &t[t_offset], ldt, &c__[c_offset], ldc, &work[1], info);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (left && tran)
    {
        /* Multiply Q to the last block of C */
        kk = (*m - *k) % (*nb - *k);
        ctr = (*m - *k) / (*nb - *k);
        if (kk > 0)
        {
            ii = *m - kk + 1;
            ztpmlqt_("L", "C", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[ii + c_dim1], ldc, &work[1], info);
        }
        else
        {
            ii = *m + 1;
        }
        i__1 = *nb + 1;
        i__2 = -(*nb - *k);
        for (i__ = ii - (*nb - *k);
                i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
                i__ += i__2)
        {
            /* Multiply Q to the current block of C (1:M,I:I+NB) */
            --ctr;
            i__3 = *nb - *k;
            ztpmlqt_("L", "C", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info);
        }
        /* Multiply Q to the first block of C (1:M,1:NB) */
        zgemlqt_("L", "C", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &c__[c_dim1 + 1], ldc, &work[1], info);
    }
    else if (left && notran)
    {
        /* Multiply Q to the first block of C */
        kk = (*m - *k) % (*nb - *k);
        ii = *m - kk + 1;
        ctr = 1;
        zgemlqt_("L", "N", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &c__[c_dim1 + 1], ldc, &work[1], info);
        i__2 = ii - *nb + *k;
        i__1 = *nb - *k;
        for (i__ = *nb + 1;
                i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
                i__ += i__1)
        {
            /* Multiply Q to the current block of C (I:I+NB,1:N) */
            i__3 = *nb - *k;
            ztpmlqt_("L", "N", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[i__ + c_dim1], ldc, &work[1], info);
            ++ctr;
        }
        if (ii <= *m)
        {
            /* Multiply Q to the last block of C */
            ztpmlqt_("L", "N", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[ii + c_dim1], ldc, &work[1], info);
        }
    }
    else if (right && notran)
    {
        /* Multiply Q to the last block of C */
        kk = (*n - *k) % (*nb - *k);
        ctr = (*n - *k) / (*nb - *k);
        if (kk > 0)
        {
            ii = *n - kk + 1;
            ztpmlqt_("R", "N", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info);
        }
        else
        {
            ii = *n + 1;
        }
        i__1 = *nb + 1;
        i__2 = -(*nb - *k);
        for (i__ = ii - (*nb - *k);
                i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
                i__ += i__2)
        {
            /* Multiply Q to the current block of C (1:M,I:I+MB) */
            --ctr;
            i__3 = *nb - *k;
            ztpmlqt_("R", "N", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info);
        }
        /* Multiply Q to the first block of C (1:M,1:MB) */
        zgemlqt_("R", "N", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &c__[c_dim1 + 1], ldc, &work[1], info);
    }
    else if (right && tran)
    {
        /* Multiply Q to the first block of C */
        kk = (*n - *k) % (*nb - *k);
        ii = *n - kk + 1;
        zgemlqt_("R", "C", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt, &c__[c_dim1 + 1], ldc, &work[1], info);
        ctr = 1;
        i__2 = ii - *nb + *k;
        i__1 = *nb - *k;
        for (i__ = *nb + 1;
                i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
                i__ += i__1)
        {
            /* Multiply Q to the current block of C (1:M,I:I+MB) */
            i__3 = *nb - *k;
            ztpmlqt_("R", "C", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[i__ * c_dim1 + 1], ldc, &work[1], info);
            ++ctr;
        }
        if (ii <= *n)
        {
            /* Multiply Q to the last block of C */
            ztpmlqt_("R", "C", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda, &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc, &c__[ii * c_dim1 + 1], ldc, &work[1], info);
        }
    }
    work[1].r = (doublereal) lw;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZLAMSWLQ */
}
/* zlamswlq_ */
