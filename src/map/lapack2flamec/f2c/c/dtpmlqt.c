/* ../netlib/v3.9.0/dtpmlqt.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DTPMLQT */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DTPMQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtpmlqt .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtpmlqt .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtpmlqt .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DTPMLQT( SIDE, TRANS, M, N, K, L, MB, V, LDV, T, LDT, */
/* A, LDA, B, LDB, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
/* INTEGER INFO, K, LDV, LDA, LDB, M, N, L, MB, LDT */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION V( LDV, * ), A( LDA, * ), B( LDB, * ), */
/* $ T( LDT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPMQRT applies a real orthogonal matrix Q obtained from a */
/* > "triangular-pentagonal" real block reflector H to a general */
/* > real matrix C, which consists of two blocks A and B. */
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
/* > The number of rows of the matrix B. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines */
/* > the matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > The order of the trapezoidal part of V. */
/* > K >= L >= 0. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* > MB is INTEGER */
/* > The block size used for the storage of T. K >= MB >= 1. */
/* > This must be the same value of MB used to generate T */
/* > in DTPLQT. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* > V is DOUBLE PRECISION array, dimension (LDV,K) */
/* > The i-th row must contain the vector which defines the */
/* > elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* > DTPLQT in B. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V. */
/* > If SIDE = 'L', LDV >= max(1,M);
*/
/* > if SIDE = 'R', LDV >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* > T is DOUBLE PRECISION array, dimension (LDT,K) */
/* > The upper triangular factors of the block reflectors */
/* > as returned by DTPLQT, stored as a MB-by-K matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= MB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension */
/* > (LDA,N) if SIDE = 'L' or */
/* > (LDA,K) if SIDE = 'R' */
/* > On entry, the K-by-N or M-by-K matrix A. */
/* > On exit, A is overwritten by the corresponding block of */
/* > Q*C or Q**T*C or C*Q or C*Q**T. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. */
/* > If SIDE = 'L', LDC >= max(1,K);
*/
/* > If SIDE = 'R', LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,N) */
/* > On entry, the M-by-N matrix B. */
/* > On exit, B is overwritten by the corresponding block of */
/* > Q*C or Q**T*C or C*Q or C*Q**T. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. */
/* > LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array. The dimension of WORK is */
/* > N*MB if SIDE = 'L', or M*MB if SIDE = 'R'. */
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
/* > \date November 2017 */
/* > \ingroup doubleOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The columns of the pentagonal matrix V contain the elementary reflectors */
/* > H(1), H(2), ..., H(K);
V is composed of a rectangular block V1 and a */
/* > trapezoidal block V2: */
/* > */
/* > V = [V1] [V2]. */
/* > */
/* > */
/* > The size of the trapezoidal block V2 is determined by the parameter L, */
/* > where 0 <= L <= K;
V2 is lower trapezoidal, consisting of the first L */
/* > rows of a K-by-K upper triangular matrix. If L=K, V2 is lower triangular;
*/
/* > if L=0, there is no trapezoidal block, hence V = V1 is rectangular. */
/* > */
/* > If SIDE = 'L': C = [A] where A is K-by-N, B is M-by-N and V is K-by-M. */
/* > [B] */
/* > */
/* > If SIDE = 'R': C = [A B] where A is M-by-K, B is M-by-N and V is K-by-N. */
/* > */
/* > The real orthogonal matrix Q is formed from V and T. */
/* > */
/* > If TRANS='N' and SIDE='L', C is on exit replaced with Q * C. */
/* > */
/* > If TRANS='T' and SIDE='L', C is on exit replaced with Q**T * C. */
/* > */
/* > If TRANS='N' and SIDE='R', C is on exit replaced with C * Q. */
/* > */
/* > If TRANS='T' and SIDE='R', C is on exit replaced with C * Q**T. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int dtpmlqt_(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, integer *mb, doublereal *v, integer *ldv, doublereal *t, integer *ldt, doublereal *a, integer *lda, doublereal * b, integer *ldb, doublereal *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dtpmlqt inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", l %" FLA_IS ", mb %" FLA_IS ", ldv %" FLA_IS ", ldt %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "",*side, *trans, *m, *n, *k, *l, *mb, *ldv, *ldt, *lda, *ldb);
    /* System generated locals */
    integer v_dim1, v_offset, a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer i__, ib, lb, nb, kf, ldaq;
    logical left, tran;
    extern logical lsame_(char *, char *);
    logical right;
    extern /* Subroutine */
    int xerbla_(char *, integer *), dtprfb_( char *, char *, char *, char *, integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
    logical notran;
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    /* Function Body */
    *info = 0;
    left = lsame_(side, "L");
    right = lsame_(side, "R");
    tran = lsame_(trans, "T");
    notran = lsame_(trans, "N");
    if (left)
    {
        ldaq = max(1,*k);
    }
    else if (right)
    {
        ldaq = max(1,*m);
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
    else if (*k < 0)
    {
        *info = -5;
    }
    else if (*l < 0 || *l > *k)
    {
        *info = -6;
    }
    else if (*mb < 1 || *mb > *k && *k > 0)
    {
        *info = -7;
    }
    else if (*ldv < *k)
    {
        *info = -9;
    }
    else if (*ldt < *mb)
    {
        *info = -11;
    }
    else if (*lda < ldaq)
    {
        *info = -13;
    }
    else if (*ldb < max(1,*m))
    {
        *info = -15;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DTPMLQT", &i__1);
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
            ib = min(i__3,i__4);
            /* Computing MIN */
            i__3 = *m - *l + i__ + ib - 1;
            nb = min(i__3,*m);
            if (i__ >= *l)
            {
                lb = 0;
            }
            else
            {
                lb = 0;
            }
            dtprfb_("L", "T", "F", "R", &nb, n, &ib, &lb, &v[i__ + v_dim1], ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &b[ b_offset], ldb, &work[1], &ib);
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
            ib = min(i__3,i__4);
            /* Computing MIN */
            i__3 = *n - *l + i__ + ib - 1;
            nb = min(i__3,*n);
            if (i__ >= *l)
            {
                lb = 0;
            }
            else
            {
                lb = nb - *n + *l - i__ + 1;
            }
            dtprfb_("R", "N", "F", "R", m, &nb, &ib, &lb, &v[i__ + v_dim1], ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], lda, &b[b_offset], ldb, &work[1], m);
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
            ib = min(i__2,i__3);
            /* Computing MIN */
            i__2 = *m - *l + i__ + ib - 1;
            nb = min(i__2,*m);
            if (i__ >= *l)
            {
                lb = 0;
            }
            else
            {
                lb = 0;
            }
            dtprfb_("L", "N", "F", "R", &nb, n, &ib, &lb, &v[i__ + v_dim1], ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ + a_dim1], lda, &b[ b_offset], ldb, &work[1], &ib);
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
            ib = min(i__2,i__3);
            /* Computing MIN */
            i__2 = *n - *l + i__ + ib - 1;
            nb = min(i__2,*n);
            if (i__ >= *l)
            {
                lb = 0;
            }
            else
            {
                lb = nb - *n + *l - i__ + 1;
            }
            dtprfb_("R", "T", "F", "R", m, &nb, &ib, &lb, &v[i__ + v_dim1], ldv, &t[i__ * t_dim1 + 1], ldt, &a[i__ * a_dim1 + 1], lda, &b[b_offset], ldb, &work[1], m);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of DTPMLQT */
}
/* dtpmlqt_ */

