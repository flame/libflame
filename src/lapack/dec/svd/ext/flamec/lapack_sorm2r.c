/* ../netlib/sorm2r.f -- translated by f2c (version 20000121). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SORM2R multiplies a general matrix by the orthogonal matrix from a QR factorization determined by sgeqrf (unblocked algorithm). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORM2R + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorm2r. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorm2r. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorm2r. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/* WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
/* INTEGER INFO, K, LDA, LDC, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORM2R overwrites the general real m by n matrix C with */
/* > */
/* > Q * C if SIDE = 'L' and TRANS = 'N', or */
/* > */
/* > Q**T* C if SIDE = 'L' and TRANS = 'T', or */
/* > */
/* > C * Q if SIDE = 'R' and TRANS = 'N', or */
/* > */
/* > C * Q**T if SIDE = 'R' and TRANS = 'T', */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by SGEQRF. Q is of order m if SIDE = 'L' and of order n */
/* > if SIDE = 'R'. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': apply Q or Q**T from the Left */
/* > = 'R': apply Q or Q**T from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': apply Q (No transpose) */
/* > = 'T': apply Q**T (Transpose) */
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
/* > A is modified by the routine but restored on exit. */
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
/* > On entry, the m by n matrix C. */
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
/* > WORK is REAL array, dimension */
/* > (N) if SIDE = 'L', */
/* > (M) if SIDE = 'R' */
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
int lapack_sorm2r(char *side, char *trans, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, real *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    /* Local variables */
    logical left;
    integer i__;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int slarf_(char *, integer *, integer *, real *, integer *, real *, real *, integer *, real *);
    integer i1, i2, i3, ic, jc, mi, ni, nq;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    logical notran;
    real aii;
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");
    /* NQ is the order of Q */
    if (left)
    {
        nq = *m;
    }
    else
    {
        nq = *n;
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
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORM2R", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0 || *k == 0)
    {
        return 0;
    }
    if (left && ! notran || ! left && notran)
    {
        i1 = 1;
        i2 = *k;
        i3 = 1;
    }
    else
    {
        i1 = *k;
        i2 = 1;
        i3 = -1;
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
        if (left)
        {
            /* H(i) is applied to C(i:m,1:n) */
            mi = *m - i__ + 1;
            ic = i__;
        }
        else
        {
            /* H(i) is applied to C(1:m,i:n) */
            ni = *n - i__ + 1;
            jc = i__;
        }
        /* Apply H(i) */
        aii = a[i__ + i__ * a_dim1];
        a[i__ + i__ * a_dim1] = 1.f;
        slarf_(side, &mi, &ni, &a[i__ + i__ * a_dim1], &c__1, &tau[i__], &c__[ ic + jc * c_dim1], ldc, &work[1]);
        a[i__ + i__ * a_dim1] = aii;
        /* L10: */
    }
    return 0;
    /* End of SORM2R */
}
/* sorm2r_ */
