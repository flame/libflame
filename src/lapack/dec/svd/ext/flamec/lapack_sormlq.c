/* ../netlib/sormlq.f -- translated by f2c (version 20000121). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLAME.h"
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;
/* > \brief \b SORMLQ */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORMLQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormlq. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormlq. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormlq. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
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
/* > SORMLQ overwrites the general real M-by-N matrix C with */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'T': Q**T * C C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* > Q = H(k) . . . H(2) H(1) */
/* > */
/* > as returned by SGELQF. Q is of order M if SIDE = 'L' and of order N */
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
/* > A is REAL array, dimension */
/* > (LDA,M) if SIDE = 'L', */
/* > (LDA,N) if SIDE = 'R' */
/* > The i-th row must contain the vector which defines the */
/* > elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* > SGELQF in the first k rows of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,K). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is REAL array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by SGELQF. */
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
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int lapack_sormlq(char *side, char *trans, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2], i__4, i__5;
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */

    /* Local variables */
    logical left;
    integer i__;
    extern logical lsame_(char *, char *);
    integer nbmin, iinfo, i1, i2, i3, ib, ic, jc;
    extern /* Subroutine */
    int lapack_sorml2(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *);
    integer nb, mi, ni, nq, nw;
    extern /* Subroutine */
    int slarfb_(char *, char *, char *, char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int slarft_(char *, char *, integer *, integer *, real *, integer *, real *, real *, integer *);
    logical notran;
    integer ldwork;
    char transt[1];
    integer lwkopt;
    logical lquery;
    integer iwt;
    int thread_id, actual_num_threads;
    integer index, mi_sub, ni_sub;
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
    lquery = *lwork == -1;
    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left)
    {
        nq = *m;
        nw = max(1,*n);
    }
    else
    {
        nq = *n;
        nw = max(1,*m);
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
    else if (*lda < max(1,*k))
    {
        *info = -7;
    }
    else if (*ldc < max(1,*m))
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
        /* Computing MIN */
        i__1 = 64;
        i__2 = ilaenv_(&c__1, "SORMLQ", ch__1, m, n, k, &c_n1); // , expr subst
        nb = min(i__1,i__2);
        lwkopt = nw * nb + 4160;
        work[1] = (real) lwkopt;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORMLQ", &i__1);
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
        if (*lwork < lwkopt)
        {
            nb = (*lwork - 4160) / ldwork;
            /* Computing MAX */
            i__1 = 2;
            i__2 = ilaenv_(&c__2, "SORMLQ", ch__1, m, n, k, &c_n1); // , expr subst
            nbmin = max(i__1,i__2);
        }
    }
    if (nb < nbmin || nb >= *k)
    {
        /* Use unblocked code */
        lapack_sorml2(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[ c_offset], ldc, &work[1], &iinfo);
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
            jc = 1;
        }
        else
        {
            mi = *m;
            ic = 1;
        }
        if (notran)
        {
            *(unsigned char *)transt = 'T';
        }
        else
        {
            *(unsigned char *)transt = 'N';
        }
        i__1 = i2;
        i__2 = i3;

#ifdef FLA_OPENMP_MULTITHREADING
        /* Get optimum thread number for DORMLQ*/
        FLA_Thread_optimum( FLA_ORMLQ, &actual_num_threads);
        #pragma omp parallel num_threads(actual_num_threads) private(i__, thread_id, mi_sub, ni_sub, index)
        {
            thread_id = omp_get_thread_num();
#else
        {
#endif
            for (i__ = i1;
                    i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
                    i__ += i__2)
            {
                /* Computing MIN */
#ifdef FLA_OPENMP_MULTITHREADING
                /* Compute triangular factor of the block reflector in a single thread */
                #pragma omp single
#endif
                {
                    i__4 = nb;
                    i__5 = *k - i__ + 1; // , expr subst
                    ib = min(i__4,i__5);
                    /* Form the triangular factor of the block reflector */
                    /* H = H(i) H(i+1) . . . H(i+ib-1) */
                    i__4 = nq - i__ + 1;
                    slarft_("Forward", "Rowwise", &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[iwt], &c__65);
                }
                if (left)
                {
                    /* H or H**T is applied to C(i:m,1:n) */
                    mi = *m - i__ + 1;
                    ic = i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, ni, &ni_sub, &index);
                    mi_sub = mi;
#endif
                }
                else
                {
                    /* H or H**T is applied to C(1:m,i:n) */
                    ni = *n - i__ + 1;
                    jc = i__;
#ifdef FLA_OPENMP_MULTITHREADING
                    /* Determine the sub partition range of current thread */
                    FLA_Thread_get_subrange(thread_id, actual_num_threads, mi, &mi_sub, &index);
                    ni_sub = ni;
#endif
                }

                /* Apply H or H**T */
#ifdef FLA_OPENMP_MULTITHREADING
                if(left)
                    slarfb_(side, transt, "Forward", "Rowwise", &mi_sub, &ni_sub, &ib, &a[i__ + i__ * a_dim1], lda, &work[iwt], &c__65, &c__[ic + (index + jc) * c_dim1], ldc, &work[1 + index], &ldwork);
                else
                    slarfb_(side, transt, "Forward", "Rowwise", &mi_sub, &ni_sub, &ib, &a[i__ + i__ * a_dim1], lda, &work[iwt], &c__65, &c__[index + ic + jc * c_dim1], ldc, &work[1 + index], &ldwork);
                #pragma omp barrier
#else
                slarfb_(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &a[i__ + i__ * a_dim1], lda, &work[iwt], &c__65, &c__[ic + jc * c_dim1], ldc, &work[1], &ldwork);
#endif
                /* L10: */
            }
        }
    }
    work[1] = (real) lwkopt;
    return 0;
    /* End of SORMLQ */
}
/* lapack_sormlq */