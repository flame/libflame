/* dgemlqt.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DGEMLQT */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGEMLQT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgemlqt .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgemlqt .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgemlqt .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGEMLQT( SIDE, TRANS, M, N, K, MB, V, LDV, T, LDT, */
/* C, LDC, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
/* INTEGER INFO, K, LDV, LDC, M, N, MB, LDT */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEMLQT overwrites the general real M-by-N matrix C with */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q C C Q */
/* > TRANS = 'T': Q**T C C Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of K */
/* > elementary reflectors: */
/* > */
/* > Q = H(1) H(2) . . . H(K) = I - V T V**T */
/* > */
/* > generated using the compact WY representation as returned by DGELQT. */
/* > */
/* > Q is of order M if SIDE = 'L' and of order N if SIDE = 'R'. */
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
/* > = 'C': Transpose, apply Q**T. */
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
/* > \param[in] MB */
/* > \verbatim */
/* > MB is INTEGER */
/* > The block size used for the storage of T. K >= MB >= 1. */
/* > This must be the same value of MB used to generate T */
/* > in DGELQT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is DOUBLE PRECISION array, dimension */
/* > (LDV,M) if SIDE = 'L', */
/* > (LDV,N) if SIDE = 'R' */
/* > The i-th row must contain the vector which defines the */
/* > elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* > DGELQT in the first K rows of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V. LDV >= fla_max(1,K). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* > T is DOUBLE PRECISION array, dimension (LDT,K) */
/* > The upper triangular factors of the block reflectors */
/* > as returned by DGELQT, stored as a MB-by-K matrix. */
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
/* > C is DOUBLE PRECISION array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q. */
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
/* > WORK is DOUBLE PRECISION array. The dimension of */
/* > WORK is N*MB if SIDE = 'L', or M*MB if SIDE = 'R'. */
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
/* > \ingroup doubleGEcomputational */
/* ===================================================================== */
/* Subroutine */
int dgemlqt_(char *side, char *trans, integer *m, integer *n, integer *k, integer *mb, doublereal *v, integer *ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc, doublereal *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgemlqt inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", mb %" FLA_IS ", ldv %" FLA_IS ", ldt %" FLA_IS ", ldc %" FLA_IS "",*side, *trans, *m, *n, *k, *mb, *ldv, *ldt, *ldc);
    /* System generated locals */
    integer v_dim1, v_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer i__, q, ib, kf;
    logical left, tran;
    extern logical lsame_(char *, char *);
    logical right;
    extern /* Subroutine */
    int dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *), xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical notran;
    integer ldwork;
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
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* .. Test the input arguments .. */
    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    left = lsame_(side, "L");
    right = lsame_(side, "R");
    tran = lsame_(trans, "T");
    notran = lsame_(trans, "N");
    if (left)
    {
        ldwork = fla_max(1,*n);
        q = *m;
    }
    else if (right)
    {
        ldwork = fla_max(1,*m);
        q = *n;
    }
    if (! left && ! right)
    {
        *info = -1;
    }
    else if (! tran && ! notran)
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
    else if (*k < 0 || *k > q)
    {
        *info = -5;
    }
    else if (*mb < 1 || *mb > *k && *k > 0)
    {
        *info = -6;
    }
    else if (*ldv < fla_max(1,*k))
    {
        *info = -8;
    }
    else if (*ldt < *mb)
    {
        *info = -10;
    }
    else if (*ldc < fla_max(1,*m))
    {
        *info = -12;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGEMLQT", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* .. Quick return if possible .. */
    if (*m == 0 || *n == 0 || *k == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (left && notran)
    {
        i__1 = *k;
        i__2 = *mb;
        for (i__ = 1;
                i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
                i__ += i__2)
        {
            /* Computing MIN */
            i__3 = *mb;
            i__4 = *k - i__ + 1; // , expr subst
            ib = fla_min(i__3,i__4);
            i__3 = *m - i__ + 1;
            dlarfb_("L", "T", "F", "R", &i__3, n, &ib, &v[i__ + i__ * v_dim1], ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, &work[1], &ldwork);
        }
    }
    else if (right && tran)
    {
        i__2 = *k;
        i__1 = *mb;
        for (i__ = 1;
                i__1 < 0 ? i__ >= i__2 : i__ <= i__2;
                i__ += i__1)
        {
            /* Computing MIN */
            i__3 = *mb;
            i__4 = *k - i__ + 1; // , expr subst
            ib = fla_min(i__3,i__4);
            i__3 = *n - i__ + 1;
            dlarfb_("R", "N", "F", "R", m, &i__3, &ib, &v[i__ + i__ * v_dim1], ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], ldc, &work[1], &ldwork);
        }
    }
    else if (left && tran)
    {
        kf = (*k - 1) / *mb * *mb + 1;
        i__1 = -(*mb);
        for (i__ = kf;
                i__1 < 0 ? i__ >= 1 : i__ <= 1;
                i__ += i__1)
        {
            /* Computing MIN */
            i__2 = *mb;
            i__3 = *k - i__ + 1; // , expr subst
            ib = fla_min(i__2,i__3);
            i__2 = *m - i__ + 1;
            dlarfb_("L", "N", "F", "R", &i__2, n, &ib, &v[i__ + i__ * v_dim1], ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ + c_dim1], ldc, &work[1], &ldwork);
        }
    }
    else if (right && notran)
    {
        kf = (*k - 1) / *mb * *mb + 1;
        i__1 = -(*mb);
        for (i__ = kf;
                i__1 < 0 ? i__ >= 1 : i__ <= 1;
                i__ += i__1)
        {
            /* Computing MIN */
            i__2 = *mb;
            i__3 = *k - i__ + 1; // , expr subst
            ib = fla_min(i__2,i__3);
            i__2 = *n - i__ + 1;
            dlarfb_("R", "T", "F", "R", m, &i__2, &ib, &v[i__ + i__ * v_dim1], ldv, &t[i__ * t_dim1 + 1], ldt, &c__[i__ * c_dim1 + 1], ldc, &work[1], &ldwork);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of DGEMLQT */
}
/* dgemlqt_ */
