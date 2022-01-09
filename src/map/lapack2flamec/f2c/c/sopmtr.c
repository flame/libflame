/* ../netlib/sopmtr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SOPMTR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SOPMTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sopmtr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sopmtr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sopmtr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS, UPLO */
/* INTEGER INFO, LDC, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL AP( * ), C( LDC, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SOPMTR overwrites the general real M-by-N matrix C with */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'T': Q**T * C C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix of order nq, with nq = m if */
/* > SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of */
/* > nq-1 elementary reflectors, as returned by SSPTRD using packed */
/* > storage: */
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
/* > = 'L': apply Q or Q**T from the Left;
*/
/* > = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangular packed storage used in previous */
/* > call to SSPTRD;
*/
/* > = 'L': Lower triangular packed storage used in previous */
/* > call to SSPTRD. */
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
/* > \param[in] AP */
/* > \verbatim */
/* > AP is REAL array, dimension */
/* > (M*(M+1)/2) if SIDE = 'L' */
/* > (N*(N+1)/2) if SIDE = 'R' */
/* > The vectors which define the elementary reflectors, as */
/* > returned by SSPTRD. AP is modified by the routine but */
/* > restored on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is REAL array, dimension (M-1) if SIDE = 'L' */
/* > or (N-1) if SIDE = 'R' */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by SSPTRD. */
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
/* > WORK is REAL array, dimension */
/* > (N) if SIDE = 'L' */
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
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int sopmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, real *ap, real *tau, real *c__, integer *ldc, real *work, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
    char buffer[256]; 
    snprintf(buffer, 256,"sopmtr inputs: side %c, uplo %c, trans %c, m %d, n %d, ldc %d",*side, *uplo, *trans, *m, *n, *ldc);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2;
    /* Local variables */
    integer i__, i1, i2, i3, ic, jc, ii, mi, ni, nq;
    real aii;
    logical left;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int slarf_(char *, integer *, integer *, real *, integer *, real *, real *, integer *, real *);
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    logical notran, forwrd;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    --ap;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");
    upper = lsame_(uplo, "U");
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
    else if (! upper && ! lsame_(uplo, "L"))
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
    else if (*ldc < max(1,*m))
    {
        *info = -9;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SOPMTR", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    if (upper)
    {
        /* Q was determined by a call to SSPTRD with UPLO = 'U' */
        forwrd = left && notran || ! left && ! notran;
        if (forwrd)
        {
            i1 = 1;
            i2 = nq - 1;
            i3 = 1;
            ii = 2;
        }
        else
        {
            i1 = nq - 1;
            i2 = 1;
            i3 = -1;
            ii = nq * (nq + 1) / 2 - 1;
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
            if (left)
            {
                /* H(i) is applied to C(1:i,1:n) */
                mi = i__;
            }
            else
            {
                /* H(i) is applied to C(1:m,1:i) */
                ni = i__;
            }
            /* Apply H(i) */
            aii = ap[ii];
            ap[ii] = 1.f;
            slarf_(side, &mi, &ni, &ap[ii - i__ + 1], &c__1, &tau[i__], &c__[ c_offset], ldc, &work[1]);
            ap[ii] = aii;
            if (forwrd)
            {
                ii = ii + i__ + 2;
            }
            else
            {
                ii = ii - i__ - 1;
            }
            /* L10: */
        }
    }
    else
    {
        /* Q was determined by a call to SSPTRD with UPLO = 'L'. */
        forwrd = left && ! notran || ! left && notran;
        if (forwrd)
        {
            i1 = 1;
            i2 = nq - 1;
            i3 = 1;
            ii = 2;
        }
        else
        {
            i1 = nq - 1;
            i2 = 1;
            i3 = -1;
            ii = nq * (nq + 1) / 2 - 1;
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
        i__2 = i2;
        i__1 = i3;
        for (i__ = i1;
                i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
                i__ += i__1)
        {
            aii = ap[ii];
            ap[ii] = 1.f;
            if (left)
            {
                /* H(i) is applied to C(i+1:m,1:n) */
                mi = *m - i__;
                ic = i__ + 1;
            }
            else
            {
                /* H(i) is applied to C(1:m,i+1:n) */
                ni = *n - i__;
                jc = i__ + 1;
            }
            /* Apply H(i) */
            slarf_(side, &mi, &ni, &ap[ii], &c__1, &tau[i__], &c__[ic + jc * c_dim1], ldc, &work[1]);
            ap[ii] = aii;
            if (forwrd)
            {
                ii = ii + nq - i__ + 1;
            }
            else
            {
                ii = ii - nq + i__ - 2;
            }
            /* L20: */
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of SOPMTR */
}
/* sopmtr_ */
