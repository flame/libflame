/* ../netlib/cunmtr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b CUNMTR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CUNMTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmtr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmtr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmtr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, */
/* WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS, UPLO */
/* INTEGER INFO, LDA, LDC, LWORK, M, N */
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
/* > CUNMTR overwrites the general complex M-by-N matrix C with */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'C': Q**H * C C * Q**H */
/* > */
/* > where Q is a complex unitary matrix of order nq, with nq = m if */
/* > SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of */
/* > nq-1 elementary reflectors, as returned by CHETRD: */
/* > */
/* > if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
*/
/* > */
/* > if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1). */
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
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A contains elementary reflectors */
/* > from CHETRD;
*/
/* > = 'L': Lower triangle of A contains elementary reflectors */
/* > from CHETRD. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': No transpose, apply Q;
*/
/* > = 'C': Conjugate transpose, apply Q**H. */
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
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension */
/* > (LDA,M) if SIDE = 'L' */
/* > (LDA,N) if SIDE = 'R' */
/* > The vectors which define the elementary reflectors, as */
/* > returned by CHETRD. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. */
/* > LDA >= fla_max(1,M) if SIDE = 'L';
LDA >= fla_max(1,N) if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension */
/* > (M-1) if SIDE = 'L' */
/* > (N-1) if SIDE = 'R' */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by CHETRD. */
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
/* > For optimum performance LWORK >= N*NB if SIDE = 'L', and */
/* > LWORK >=M*NB if SIDE = 'R', where NB is the optimal */
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
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int cunmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, complex *a, integer *lda, complex *tau, complex *c__, integer *ldc, complex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"cunmtr inputs: side %c, uplo %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS ", lwork %" FLA_IS "",*side, *uplo, *trans, *m, *n, *lda, *ldc, *lwork);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__2, i__3;
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */

    /* Local variables */
    integer i1, i2, nb, mi, ni, nq, nw;
    logical left;
    extern logical lsame_(char *, char *);
    integer iinfo;
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int cunmql_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *), cunmqr_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *);
    integer lwkopt;
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
    upper = lsame_(uplo, "U");
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
    else if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -2;
    }
    else if (! lsame_(trans, "N") && ! lsame_(trans, "C"))
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
    else if (*lda < fla_max(1,nq))
    {
        *info = -7;
    }
    else if (*ldc < fla_max(1,*m))
    {
        *info = -10;
    }
    else if (*lwork < fla_max(1,nw) && ! lquery)
    {
        *info = -12;
    }
    if (*info == 0)
    {
        if (upper)
        {
            if (left)
            {
                i__2 = *m - 1;
                i__3 = *m - 1;
                nb = ilaenv_(&c__1, "CUNMQL", ch__1, &i__2, n, &i__3, &c_n1);
            }
            else
            {
                i__2 = *n - 1;
                i__3 = *n - 1;
                nb = ilaenv_(&c__1, "CUNMQL", ch__1, m, &i__2, &i__3, &c_n1);
            }
        }
        else
        {
            if (left)
            {
                i__2 = *m - 1;
                i__3 = *m - 1;
                nb = ilaenv_(&c__1, "CUNMQR", ch__1, &i__2, n, &i__3, &c_n1);
            }
            else
            {
                i__2 = *n - 1;
                i__3 = *n - 1;
                nb = ilaenv_(&c__1, "CUNMQR", ch__1, m, &i__2, &i__3, &c_n1);
            }
        }
        lwkopt = fla_max(1,nw) * nb;
        work[1].r = (real) lwkopt;
        work[1].i = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__2 = -(*info);
        xerbla_("CUNMTR", &i__2);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0 || nq == 1)
    {
        work[1].r = 1.f;
        work[1].i = 0.f; // , expr subst
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    if (left)
    {
        mi = *m - 1;
        ni = *n;
    }
    else
    {
        mi = *m;
        ni = *n - 1;
    }
    if (upper)
    {
        /* Q was determined by a call to CHETRD with UPLO = 'U' */
        i__2 = nq - 1;
        cunmql_(side, trans, &mi, &ni, &i__2, &a[(a_dim1 << 1) + 1], lda, & tau[1], &c__[c_offset], ldc, &work[1], lwork, &iinfo);
    }
    else
    {
        /* Q was determined by a call to CHETRD with UPLO = 'L' */
        if (left)
        {
            i1 = 2;
            i2 = 1;
        }
        else
        {
            i1 = 1;
            i2 = 2;
        }
        i__2 = nq - 1;
        cunmqr_(side, trans, &mi, &ni, &i__2, &a[a_dim1 + 2], lda, &tau[1], & c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo);
    }
    work[1].r = (real) lwkopt;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CUNMTR */
}
/* cunmtr_ */
