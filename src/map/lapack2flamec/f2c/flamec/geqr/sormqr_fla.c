/* ../netlib/sormqr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;
/* > \brief \b SORMQR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORMQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormqr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormqr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormqr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/* WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
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
/* > SORMQR overwrites the general real M-by-N matrix C with */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'T': Q**T * C C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by SGEQRF. Q is of order M if SIDE = 'L' and of order N */
/* > if SIDE = 'R'. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': apply Q or Q**T from the Left;
*/
/* > = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': No transpose, apply Q;
*/
/* > = 'T': Transpose, apply Q**T. */
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
/* > A is REAL array, dimension (LDA,K) */
/* > The i-th column must contain the vector which defines the */
/* > elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* > SGEQRF in the first k columns of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. */
/* > If SIDE = 'L', LDA >= max(1,M);
*/
/* > if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is REAL array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by SGEQRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is REAL array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= max(1,M). */
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
/* > If SIDE = 'L', LWORK >= max(1,N);
*/
/* > if SIDE = 'R', LWORK >= max(1,M). */
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
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int sormqr_fla(char *side, char *trans, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__4, i__5;
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */

    /* Local variables */
    integer i__;
    real t[4160] /* was [65][64] */
    ;
    integer i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iws;
    logical left;
    extern logical lsame_(char *, char *);
    integer nbmin, iinfo;
    extern /* Subroutine */
    int sorm2r_fla(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *), slarfb_(char *, char *, char *, char * , integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int slarft_(char *, char *, integer *, integer *, real *, integer *, real *, real *, integer *);
    logical notran;
    integer ldwork, lwkopt;
    logical lquery;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Local Arrays .. */
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
        nw = *n;
    }
    else
    {
        nq = *n;
        nw = *m;
    }
    if (! left && ! lsame_(side, "R"))
    {
        *info = -1;
    }
    else if (! notran && ! lsame_(trans, "T"))
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
    else if (*lda < max(1,nq))
    {
        *info = -7;
    }
    else if (*ldc < max(1,*m))
    {
        *info = -10;
    }
    else if (*lwork < max(1,nw) && ! lquery)
    {
        *info = -12;
    }
    if (*info == 0)
    {
        /* Determine the block size. NB may be at most NBMAX, where NBMAX */
        /* is used to define the local array T. */
        /* Computing MIN */
        i__1 = 64;
        i__2 = ilaenv_(&c__1, "SORMQR", ch__1, m, n, k, &c_n1); // , expr subst
        nb = min(i__1,i__2);
        lwkopt = max(1,nw) * nb;
        work[1] = (real) lwkopt;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORMQR", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0 || *k == 0)
    {
        work[1] = 1.f;
        return 0;
    }
    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k)
    {
        iws = nw * nb;
        if (*lwork < iws)
        {
            nb = *lwork / ldwork;
            /* Computing MAX */
            i__1 = 2;
            i__2 = ilaenv_(&c__2, "SORMQR", ch__1, m, n, k, &c_n1); // , expr subst
            nbmin = max(i__1,i__2);
        }
    }
    else
    {
        iws = nw;
    }
    if (nb < nbmin || nb >= *k)
    {
        /* Use unblocked code */
        sorm2r_fla(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[ c_offset], ldc, &work[1], &iinfo);
    }
    else
    {
        /* Use blocked code */
        if (left && ! notran || ! left && notran)
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
            jc = 1;
        }
        else
        {
            mi = *m;
            ic = 1;
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
            ib = min(i__4,i__5);
            /* Form the triangular factor of the block reflector */
            /* H = H(i) H(i+1) . . . H(i+ib-1) */
            i__4 = nq - i__ + 1;
            slarft_("Forward", "Columnwise", &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], t, &c__65) ;
            if (left)
            {
                /* H or H**T is applied to C(i:m,1:n) */
                mi = *m - i__ + 1;
                ic = i__;
            }
            else
            {
                /* H or H**T is applied to C(1:m,i:n) */
                ni = *n - i__ + 1;
                jc = i__;
            }
            /* Apply H or H**T */
            slarfb_(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &a[ i__ + i__ * a_dim1], lda, t, &c__65, &c__[ic + jc * c_dim1], ldc, &work[1], &ldwork);
            /* L10: */
        }
    }
    work[1] = (real) lwkopt;
    return 0;
    /* End of SORMQR */
}
/* sormqr_ */
