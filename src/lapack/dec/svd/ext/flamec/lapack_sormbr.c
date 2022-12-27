/* ../netlib/sormbr.f -- translated by f2c (version 20000121). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
/* > \brief \b SORMBR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORMBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormbr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormbr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormbr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, */
/* LDC, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS, VECT */
/* INTEGER INFO, K, LDA, LDC, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), C( LDC, * ), TAU( * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > If VECT = 'Q', SORMBR overwrites the general real M-by-N matrix C */
/* > with */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'T': Q**T * C C * Q**T */
/* > */
/* > If VECT = 'P', SORMBR overwrites the general real M-by-N matrix C */
/* > with */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': P * C C * P */
/* > TRANS = 'T': P**T * C C * P**T */
/* > */
/* > Here Q and P**T are the orthogonal matrices determined by SGEBRD when */
/* > reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and */
/* > P**T are defined as products of elementary reflectors H(i) and G(i) */
/* > respectively. */
/* > */
/* > Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the */
/* > order of the orthogonal matrix Q or P**T that is applied. */
/* > */
/* > If VECT = 'Q', A is assumed to have been an NQ-by-K matrix: */
/* > if nq >= k, Q = H(1) H(2) . . . H(k);
*/
/* > if nq < k, Q = H(1) H(2) . . . H(nq-1). */
/* > */
/* > If VECT = 'P', A is assumed to have been a K-by-NQ matrix: */
/* > if k < nq, P = G(1) G(2) . . . G(k);
*/
/* > if k >= nq, P = G(1) G(2) . . . G(nq-1). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] VECT */
/* > \verbatim */
/* > VECT is CHARACTER*1 */
/* > = 'Q': apply Q or Q**T;
*/
/* > = 'P': apply P or P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': apply Q, Q**T, P or P**T from the Left;
*/
/* > = 'R': apply Q, Q**T, P or P**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': No transpose, apply Q or P;
*/
/* > = 'T': Transpose, apply Q**T or P**T. */
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
/* > If VECT = 'Q', the number of columns in the original */
/* > matrix reduced by SGEBRD. */
/* > If VECT = 'P', the number of rows in the original */
/* > matrix reduced by SGEBRD. */
/* > K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension */
/* > (LDA,fla_min(nq,K)) if VECT = 'Q' */
/* > (LDA,nq) if VECT = 'P' */
/* > The vectors which define the elementary reflectors H(i) and */
/* > G(i), whose products determine the matrices Q and P, as */
/* > returned by SGEBRD. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. */
/* > If VECT = 'Q', LDA >= fla_max(1,nq);
*/
/* > if VECT = 'P', LDA >= fla_max(1,fla_min(nq,K)). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is REAL array, dimension (fla_min(nq,K)) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i) or G(i) which determines Q or P, as returned */
/* > by SGEBRD in the array argument TAUQ or TAUP. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is REAL array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q */
/* > or P*C or P**T*C or C*P or C*P**T. */
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
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
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
/* > For optimum performance LWORK >= N*NB if SIDE = 'L', and */
/* > LWORK >= M*NB if SIDE = 'R', where NB is the optimal */
/* > blocksize. */
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
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int lapack_sormbr(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2];
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */

    /* Local variables */
    logical left;
    extern logical lsame_(char *, char *);
    integer iinfo, i1, i2, nb, mi, ni, nq, nw;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    logical notran, applyq;
    char transt[1];
    extern /* Subroutine */
    int lapack_sormlq_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */
    int lapack_sormqr(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *);
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
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
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    applyq = lsame_(vect, "Q");
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");
    lquery = *lwork == -1;
    /* NQ is the order of Q or P and NW is the minimum dimension of WORK */
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
    if (! applyq && ! lsame_(vect, "P"))
    {
        *info = -1;
    }
    else if (! left && ! lsame_(side, "R"))
    {
        *info = -2;
    }
    else if (! notran && ! lsame_(trans, "T"))
    {
        *info = -3;
    }
    else if (*m < 0)
    {
        *info = -4;
    }
    else if (*n < 0)
    {
        *info = -5;
    }
    else if (*k < 0)
    {
        *info = -6;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = fla_min(nq,*k); // , expr subst
        if (applyq && *lda < fla_max(1,nq) || ! applyq && *lda < fla_max(i__1,i__2))
        {
            *info = -8;
        }
        else if (*ldc < fla_max(1,*m))
        {
            *info = -11;
        }
        else if (*lwork < nw && ! lquery)
        {
            *info = -13;
        }
    }
    if (*info == 0)
    {
        if (applyq)
        {
            if (left)
            {
                i__1 = *m - 1;
                i__2 = *m - 1;
                nb = ilaenv_(&c__1, "SORMQR", ch__1, &i__1, n, &i__2, &c_n1);
            }
            else
            {
                i__1 = *n - 1;
                i__2 = *n - 1;
                nb = ilaenv_(&c__1, "SORMQR", ch__1, m, &i__1, &i__2, &c_n1);
            }
        }
        else
        {
            if (left)
            {
                i__1 = *m - 1;
                i__2 = *m - 1;
                nb = ilaenv_(&c__1, "SORMLQ", ch__1, &i__1, n, &i__2, &c_n1);
            }
            else
            {
                i__1 = *n - 1;
                i__2 = *n - 1;
                nb = ilaenv_(&c__1, "SORMLQ", ch__1, m, &i__1, &i__2, &c_n1);
            }
        }
        lwkopt = nw * nb;
        work[1] = (real) lwkopt;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORMBR", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    work[1] = 1.f;
    if (*m == 0 || *n == 0)
    {
        return 0;
    }
    if (applyq)
    {
        /* Apply Q */
        if (nq >= *k)
        {
            /* Q was determined by a call to SGEBRD with nq >= k */
            lapack_sormqr(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[ c_offset], ldc, &work[1], lwork, &iinfo);
        }
        else if (nq > 1)
        {
            /* Q was determined by a call to SGEBRD with nq < k */
            if (left)
            {
                mi = *m - 1;
                ni = *n;
                i1 = 2;
                i2 = 1;
            }
            else
            {
                mi = *m;
                ni = *n - 1;
                i1 = 1;
                i2 = 2;
            }
            i__1 = nq - 1;
            lapack_sormqr(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo);
        }
    }
    else
    {
        /* Apply P */
        if (notran)
        {
            *(unsigned char *)transt = 'T';
        }
        else
        {
            *(unsigned char *)transt = 'N';
        }
        if (nq > *k)
        {
            /* P was determined by a call to SGEBRD with nq > k */
            lapack_sormlq(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[ c_offset], ldc, &work[1], lwork, &iinfo);
        }
        else if (nq > 1)
        {
            /* P was determined by a call to SGEBRD with nq <= k */
            if (left)
            {
                mi = *m - 1;
                ni = *n;
                i1 = 2;
                i2 = 1;
            }
            else
            {
                mi = *m;
                ni = *n - 1;
                i1 = 1;
                i2 = 2;
            }
            i__1 = nq - 1;
            lapack_sormlq(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda, &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, & iinfo);
        }
    }
    work[1] = (real) lwkopt;
    return 0;
    /* End of SORMBR */
}
/* lapack_sormbr */
