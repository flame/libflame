/* ../netlib/dtpqrt.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DTPQRT */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DTPQRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtpqrt. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtpqrt. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtpqrt. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DTPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDB, LDT, N, M, L, NB */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPQRT computes a blocked QR factorization of a real */
/* > "triangular-pentagonal" matrix C, which is composed of a */
/* > triangular block A and pentagonal block B, using the compact */
/* > WY representation for Q. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix B. */
/* > M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix B, and the order of the */
/* > triangular matrix A. */
/* > N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > The number of rows of the upper trapezoidal part of B. */
/* > MIN(M,N) >= L >= 0. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The block size to be used in the blocked QR. N >= NB >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the upper triangular N-by-N matrix A. */
/* > On exit, the elements on and above the diagonal of the array */
/* > contain the upper triangular matrix R. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,N) */
/* > On entry, the pentagonal M-by-N matrix B. The first M-L rows */
/* > are rectangular, and the last L rows are upper trapezoidal. */
/* > On exit, B contains the pentagonal matrix V. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is DOUBLE PRECISION array, dimension (LDT,N) */
/* > The upper triangular block reflectors stored in compact form */
/* > as a sequence of upper triangular blocks. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= NB. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (NB*N) */
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
/* > \date November 2013 */
/* > \ingroup doubleOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The input matrix C is a (N+M)-by-N matrix */
/* > */
/* > C = [ A ] */
/* > [ B ] */
/* > */
/* > where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal */
/* > matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N */
/* > upper trapezoidal matrix B2: */
/* > */
/* > B = [ B1 ] <- (M-L)-by-N rectangular */
/* > [ B2 ] <- L-by-N upper trapezoidal. */
/* > */
/* > The upper trapezoidal matrix B2 consists of the first L rows of a */
/* > N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N). If L=0, */
/* > B is rectangular M-by-N;
if M=L=N, B is upper triangular. */
/* > */
/* > The matrix W stores the elementary reflectors H(i) in the i-th column */
/* > below the diagonal (of A) in the (N+M)-by-N input matrix C */
/* > */
/* > C = [ A ] <- upper triangular N-by-N */
/* > [ B ] <- M-by-N pentagonal */
/* > */
/* > so that W can be represented as */
/* > */
/* > W = [ I ] <- identity, N-by-N */
/* > [ V ] <- M-by-N, same form as B. */
/* > */
/* > Thus, all of information needed for W is contained on exit in B, which */
/* > we call V above. Note that V has the same form as B;
that is, */
/* > */
/* > V = [ V1 ] <- (M-L)-by-N rectangular */
/* > [ V2 ] <- L-by-N upper trapezoidal. */
/* > */
/* > The columns of V represent the vectors which define the H(i)'s. */
/* > */
/* > The number of blocks is B = ceiling(N/NB), where each */
/* > block is of order NB except for the last block, which is of order */
/* > IB = N - (B-1)*NB. For each of the B blocks, a upper triangular block */
/* > reflector factor is computed: T1, T2, ..., TB. The NB-by-NB (and IB-by-IB */
/* > for the last block) T's are stored in the NB-by-N matrix T as */
/* > */
/* > T = [T1 T2 ... TB]. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int dtpqrt_(integer *m, integer *n, integer *l, integer *nb, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal * t, integer *ldt, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, ib, lb, mb, iinfo;
    extern /* Subroutine */
    int xerbla_(char *, integer *), dtprfb_( char *, char *, char *, char *, integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *), dtpqrt2_(integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *);
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;
    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*l < 0 || *l > min(*m,*n) && min(*m,*n) >= 0)
    {
        *info = -3;
    }
    else if (*nb < 1 || *nb > *n && *n > 0)
    {
        *info = -4;
    }
    else if (*lda < max(1,*n))
    {
        *info = -6;
    }
    else if (*ldb < max(1,*m))
    {
        *info = -8;
    }
    else if (*ldt < *nb)
    {
        *info = -10;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DTPQRT", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        return 0;
    }
    i__1 = *n;
    i__2 = *nb;
    for (i__ = 1;
            i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
            i__ += i__2)
    {
        /* Compute the QR factorization of the current block */
        /* Computing MIN */
        i__3 = *n - i__ + 1;
        ib = min(i__3,*nb);
        /* Computing MIN */
        i__3 = *m - *l + i__ + ib - 1;
        mb = min(i__3,*m);
        if (i__ >= *l)
        {
            lb = 0;
        }
        else
        {
            lb = mb - *m + *l - i__ + 1;
        }
        dtpqrt2_(&mb, &ib, &lb, &a[i__ + i__ * a_dim1], lda, &b[i__ * b_dim1 + 1], ldb, &t[i__ * t_dim1 + 1], ldt, &iinfo);
        /* Update by applying H**T to B(:,I+IB:N) from the left */
        if (i__ + ib <= *n)
        {
            i__3 = *n - i__ - ib + 1;
            dtprfb_("L", "T", "F", "C", &mb, &i__3, &ib, &lb, &b[i__ * b_dim1 + 1], ldb, &t[i__ * t_dim1 + 1], ldt, &a[i__ + (i__ + ib) * a_dim1], lda, &b[(i__ + ib) * b_dim1 + 1], ldb, &work[1] , &ib);
        }
    }
    return 0;
    /* End of DTPQRT */
}
/* dtpqrt_ */
