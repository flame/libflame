/* ../netlib/v3.9.0/zgeqr.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
/* > \brief \b ZGEQR */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N, TSIZE, LWORK */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), T( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEQR computes a QR factorization of a complex M-by-N matrix A: */
/* > */
/* > A = Q * ( R ), */
/* > ( 0 ) */
/* > */
/* > where: */
/* > */
/* > Q is a M-by-M orthogonal matrix;
*/
/* > R is an upper-triangular N-by-N matrix;
*/
/* > 0 is a (M-N)-by-N zero matrix, if M > N. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the elements on and above the diagonal of the array */
/* > contain the min(M,N)-by-N upper trapezoidal matrix R */
/* > (R is upper triangular if M >= N);
*/
/* > the elements below the diagonal are used to store part of the */
/* > data structure to represent Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is COMPLEX*16 array, dimension (MAX(5,TSIZE)) */
/* > On exit, if INFO = 0, T(1) returns optimal (or either minimal */
/* > or optimal, if query is assumed) TSIZE. See TSIZE for details. */
/* > Remaining T contains part of the data structure used to represent Q. */
/* > If one wants to apply or construct Q, then one needs to keep T */
/* > (in addition to A) and pass it to further subroutines. */
/* > \endverbatim */
/* > */
/* > \param[in] TSIZE */
/* > \verbatim */
/* > TSIZE is INTEGER */
/* > If TSIZE >= 5, the dimension of the array T. */
/* > If TSIZE = -1 or -2, then a workspace query is assumed. The routine */
/* > only calculates the sizes of the T and WORK arrays, returns these */
/* > values as the first entries of the T and WORK arrays, and no error */
/* > message related to T or WORK is issued by XERBLA. */
/* > If TSIZE = -1, the routine calculates optimal size of T for the */
/* > optimum performance and returns this value in T(1). */
/* > If TSIZE = -2, the routine calculates minimal size of T and */
/* > returns this value in T(1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) contains optimal (or either minimal */
/* > or optimal, if query was assumed) LWORK. */
/* > See LWORK for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If LWORK = -1 or -2, then a workspace query is assumed. The routine */
/* > only calculates the sizes of the T and WORK arrays, returns these */
/* > values as the first entries of the T and WORK arrays, and no error */
/* > message related to T or WORK is issued by XERBLA. */
/* > If LWORK = -1, the routine calculates optimal size of WORK for the */
/* > optimal performance and returns this value in WORK(1). */
/* > If LWORK = -2, the routine calculates minimal size of WORK and */
/* > returns this value in WORK(1). */
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
/* > \par Further Details */
/* ==================== */
/* > */
/* > \verbatim */
/* > */
/* > The goal of the interface is to give maximum freedom to the developers for */
/* > creating any QR factorization algorithm they wish. The triangular */
/* > (trapezoidal) R has to be stored in the upper part of A. The lower part of A */
/* > and the array T can be used to store any relevant information for applying or */
/* > constructing the Q factor. The WORK array can safely be discarded after exit. */
/* > */
/* > Caution: One should not expect the sizes of T and WORK to be the same from one */
/* > LAPACK implementation to the other, or even from one execution to the other. */
/* > A workspace query (for T and WORK) is needed at each execution. However, */
/* > for a given execution, the size of T and WORK are fixed and will not change */
/* > from one query to the next. */
/* > */
/* > \endverbatim */
/* > */
/* > \par Further Details particular to this LAPACK implementation: */
/* ============================================================== */
/* > */
/* > \verbatim */
/* > */
/* > These details are particular for this LAPACK implementation. Users should not */
/* > take them for granted. These details may change in the future, and are not likely */
/* > true for another LAPACK implementation. These details are relevant if one wants */
/* > to try to understand the code. They are not part of the interface. */
/* > */
/* > In this version, */
/* > */
/* > T(2): row block size (MB) */
/* > T(3): column block size (NB) */
/* > T(6:TSIZE): data structure needed for Q, computed by */
/* > ZLATSQR or ZGEQRT */
/* > */
/* > Depending on the matrix dimensions M and N, and row and column */
/* > block sizes MB and NB returned by ILAENV, ZGEQR will use either */
/* > ZLATSQR (if the matrix is tall-and-skinny) or ZGEQRT to compute */
/* > the QR factorization. */
/* > */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int zgeqr_(integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *t, integer *tsize, doublecomplex *work, integer * lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"zgeqr inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", tsize %" FLA_IS ", lwork %" FLA_IS "",*m, *n, *lda, *tsize, *lwork);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer mb, nb;
    logical mint, minw;
    integer nblcks;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    logical lminws;
    extern /* Subroutine */
    int zgeqrt_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    logical lquery;
    integer mintsz;
    extern /* Subroutine */
    int zlatsqr_(integer *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    /* -- LAPACK computational routine (version 3.9.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. -- */
    /* November 2019 */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --t;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *tsize == -1 || *tsize == -2 || *lwork == -1 || *lwork == -2;
    mint = FALSE_;
    minw = FALSE_;
    if (*tsize == -2 || *lwork == -2)
    {
        if (*tsize != -1)
        {
            mint = TRUE_;
        }
        if (*lwork != -1)
        {
            minw = TRUE_;
        }
    }
    /* Determine the block size */
    if (min(*m,*n) > 0)
    {
        mb = ilaenv_(&c__1, "ZGEQR ", " ", m, n, &c__1, &c_n1);
        nb = ilaenv_(&c__1, "ZGEQR ", " ", m, n, &c__2, &c_n1);
    }
    else
    {
        mb = *m;
        nb = 1;
    }
    if (mb > *m || mb <= *n)
    {
        mb = *m;
    }
    if (nb > min(*m,*n) || nb < 1)
    {
        nb = 1;
    }
    mintsz = *n + 5;
    if (mb > *n && *m > *n)
    {
        if ((*m - *n) % (mb - *n) == 0)
        {
            nblcks = (*m - *n) / (mb - *n);
        }
        else
        {
            nblcks = (*m - *n) / (mb - *n) + 1;
        }
    }
    else
    {
        nblcks = 1;
    }
    /* Determine if the workspace size satisfies minimal size */
    lminws = FALSE_;
    /* Computing MAX */
    i__1 = 1;
    i__2 = nb * *n * nblcks + 5; // , expr subst
    if ((*tsize < max(i__1,i__2) || *lwork < nb * *n) && *lwork >= *n && * tsize >= mintsz && ! lquery)
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = nb * *n * nblcks + 5; // , expr subst
        if (*tsize < max(i__1,i__2))
        {
            lminws = TRUE_;
            nb = 1;
            mb = *m;
        }
        if (*lwork < nb * *n)
        {
            lminws = TRUE_;
            nb = 1;
        }
    }
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < max(1,*m))
    {
        *info = -4;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = nb * *n * nblcks + 5; // , expr subst
        if (*tsize < max(i__1,i__2) && ! lquery && ! lminws)
        {
            *info = -6;
        }
        else /* if(complicated condition) */
        {
            /* Computing MAX */
            i__1 = 1;
            i__2 = *n * nb; // , expr subst
            if (*lwork < max(i__1,i__2) && ! lquery && ! lminws)
            {
                *info = -8;
            }
        }
    }
    if (*info == 0)
    {
        if (mint)
        {
            t[1].r = (doublereal) mintsz;
            t[1].i = 0.; // , expr subst
        }
        else
        {
            i__1 = nb * *n * nblcks + 5;
            t[1].r = (doublereal) i__1;
            t[1].i = 0.; // , expr subst
        }
        t[2].r = (doublereal) mb;
        t[2].i = 0.; // , expr subst
        t[3].r = (doublereal) nb;
        t[3].i = 0.; // , expr subst
        if (minw)
        {
            i__1 = max(1,*n);
            work[1].r = (doublereal) i__1;
            work[1].i = 0.; // , expr subst
        }
        else
        {
            /* Computing MAX */
            i__2 = 1;
            i__3 = nb * *n; // , expr subst
            i__1 = max(i__2,i__3);
            work[1].r = (doublereal) i__1;
            work[1].i = 0.; // , expr subst
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGEQR", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (min(*m,*n) == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* The QR Decomposition */
    if (*m <= *n || mb <= *n || mb >= *m)
    {
        zgeqrt_(m, n, &nb, &a[a_offset], lda, &t[6], &nb, &work[1], info);
    }
    else
    {
        zlatsqr_(m, n, &mb, &nb, &a[a_offset], lda, &t[6], &nb, &work[1], lwork, info);
    }
    /* Computing MAX */
    i__2 = 1;
    i__3 = nb * *n; // , expr subst
    i__1 = max(i__2,i__3);
    work[1].r = (doublereal) i__1;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of ZGEQR */
}
/* zgeqr_ */

