/* ../netlib/cunmql.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;
/* > \brief \b CUNMQL */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CUNMQL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmql. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmql. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmql. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/* WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
/* INTEGER INFO, K, LDA, LDC, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), C( LDC, * ), TAU( * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNMQL overwrites the general complex M-by-N matrix C with */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'C': Q**H * C C * Q**H */
/* > */
/* > where Q is a complex unitary matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* > Q = H(k) . . . H(2) H(1) */
/* > */
/* > as returned by CGEQLF. Q is of order M if SIDE = 'L' and of order N */
/* > if SIDE = 'R'. */
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
/* > = 'C': Transpose, apply Q**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix C. M >= 0. */
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
/* > If SIDE = 'L', M >= K >= 0;
*/
/* > if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,K) */
/* > The i-th column must contain the vector which defines the */
/* > elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* > CGEQLF in the last k columns of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. */
/* > If SIDE = 'L', LDA >= fla_max(1,M);
*/
/* > if SIDE = 'R', LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by CGEQLF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX array, dimension (LDC,N) */
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
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If SIDE = 'L', LWORK >= fla_max(1,N);
*/
/* > if SIDE = 'R', LWORK >= fla_max(1,M). */
/* > For good performance, LWORK should generally be larger. */
/* > */
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
/* > \date December 2016 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int cunmql_(char *side, char *trans, integer *m, integer *n, integer *k, complex *a, integer *lda, complex *tau, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"cunmql inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS ", lwork %" FLA_IS "",*side, *trans, *m, *n, *k, *lda, *ldc, *lwork);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2], i__4, i__5;
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */

    /* Local variables */
    integer i__, i1, i2, i3, ib, nb, mi, ni, nq, nw, iwt;
    logical left;
    extern logical lsame_(char *, char *);
    integer nbmin, iinfo;
    extern /* Subroutine */
    int cunm2l_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *), clarfb_(char *, char *, char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *), clarft_(char *, char *, integer *, integer *, complex *, integer *, complex *, complex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    logical notran;
    integer ldwork, lwkopt;
    logical lquery;
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
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
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");
    lquery = *lwork == -1;
    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left)
    {
        nq = *m;
        nw = fla_max(1,*n);
    }
    else
    {
        nq = *n;
        nw = fla_max(1,*m);
    }
    if (! left && ! lsame_(side, "R"))
    {
        *info = -1;
    }
    else if (! notran && ! lsame_(trans, "C"))
    {
        *info = -2;
    }
    else if (*m < 0)
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*k < 0 || *k > nq)
    {
        *info = -5;
    }
    else if (*lda < fla_max(1,nq))
    {
        *info = -7;
    }
    else if (*ldc < fla_max(1,*m))
    {
        *info = -10;
    }
    else if (*lwork < nw && ! lquery)
    {
        *info = -12;
    }
    if (*info == 0)
    {
        /* Compute the workspace requirements */
        if (*m == 0 || *n == 0)
        {
            lwkopt = 1;
        }
        else
        {
            /* Computing MIN */
            i__1 = 64;
            i__2 = ilaenv_(&c__1, "CUNMQL", ch__1, m, n, k, &c_n1); // , expr subst
            nb = fla_min(i__1,i__2);
            lwkopt = nw * nb + 4160;
        }
        work[1].r = (real) lwkopt;
        work[1].i = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CUNMQL", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Determine the block size */
    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k)
    {
        if (*lwork < nw * nb + 4160)
        {
            nb = (*lwork - 4160) / ldwork;
            /* Computing MAX */
            i__1 = 2;
            i__2 = ilaenv_(&c__2, "CUNMQL", ch__1, m, n, k, &c_n1); // , expr subst
            nbmin = fla_max(i__1,i__2);
        }
    }
    if (nb < nbmin || nb >= *k)
    {
        /* Use unblocked code */
        cunm2l_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[ c_offset], ldc, &work[1], &iinfo);
    }
    else
    {
        /* Use blocked code */
        iwt = nw * nb + 1;
        if (left && notran || ! left && ! notran)
        {
            i1 = 1;
            i2 = *k;
            i3 = nb;
        }
        else
        {
            i1 = (*k - 1) / nb * nb + 1;
            i2 = 1;
            i3 = -nb;
        }
        if (left)
        {
            ni = *n;
        }
        else
        {
            mi = *m;
        }
        i__1 = i2;
        i__2 = i3;
        for (i__ = i1;
                i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
                i__ += i__2)
        {
            /* Computing MIN */
            i__4 = nb;
            i__5 = *k - i__ + 1; // , expr subst
            ib = fla_min(i__4,i__5);
            /* Form the triangular factor of the block reflector */
            /* H = H(i+ib-1) . . . H(i+1) H(i) */
            i__4 = nq - *k + i__ + ib - 1;
            clarft_("Backward", "Columnwise", &i__4, &ib, &a[i__ * a_dim1 + 1], lda, &tau[i__], &work[iwt], &c__65);
            if (left)
            {
                /* H or H**H is applied to C(1:m-k+i+ib-1,1:n) */
                mi = *m - *k + i__ + ib - 1;
            }
            else
            {
                /* H or H**H is applied to C(1:m,1:n-k+i+ib-1) */
                ni = *n - *k + i__ + ib - 1;
            }
            /* Apply H or H**H */
            clarfb_(side, trans, "Backward", "Columnwise", &mi, &ni, &ib, &a[ i__ * a_dim1 + 1], lda, &work[iwt], &c__65, &c__[c_offset], ldc, &work[1], &ldwork);
            /* L10: */
        }
    }
    work[1].r = (real) lwkopt;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CUNMQL */
}
/* cunmql_ */
