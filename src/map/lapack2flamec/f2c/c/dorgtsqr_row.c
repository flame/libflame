/* dorgtsqr_row.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b4 = 0.;
static doublereal c_b5 = 1.;
static integer c__0 = 0;
static integer c__1 = 1;
/* > \brief \b DORGTSQR_ROW */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DORGTSQR_ROW + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgtsq r_row.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgtsq r_row.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgtsq r_row.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DORGTSQR_ROW( M, N, MB, NB, A, LDA, T, LDT, WORK, */
/* $ LWORK, INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDT, LWORK, M, N, MB, NB */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), T( LDT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORGTSQR_ROW generates an M-by-N real matrix Q_out with */
/* > orthonormal columns from the output of DLATSQR. These N orthonormal */
/* > columns are the first N columns of a product of complex unitary */
/* > matrices Q(k)_in of order M, which are returned by DLATSQR in */
/* > a special format. */
/* > */
/* > Q_out = first_N_columns_of( Q(1)_in * Q(2)_in * ... * Q(k)_in ). */
/* > */
/* > The input matrices Q(k)_in are stored in row and column blocks in A. */
/* > See the documentation of DLATSQR for more details on the format of */
/* > Q(k)_in, where each Q(k)_in is represented by block Householder */
/* > transformations. This routine calls an auxiliary routine DLARFB_GETT, */
/* > where the computation is performed on each individual block. The */
/* > algorithm first sweeps NB-sized column blocks from the right to left */
/* > starting in the bottom row block and continues to the top row block */
/* > (hence _ROW in the routine name). This sweep is in reverse order of */
/* > the order in which DLATSQR generates the output blocks. */
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
/* > The number of columns of the matrix A. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] MB */
/* > \verbatim */
/* > MB is INTEGER */
/* > The row block size used by DLATSQR to return */
/* > arrays A and T. MB > N. */
/* > (Note that if MB > M, then M is used instead of MB */
/* > as the row block size). */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The column block size used by DLATSQR to return */
/* > arrays A and T. NB >= 1. */
/* > (Note that if NB > N, then N is used instead of NB */
/* > as the column block size). */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > */
/* > On entry: */
/* > */
/* > The elements on and above the diagonal are not used as */
/* > input. The elements below the diagonal represent the unit */
/* > lower-trapezoidal blocked matrix V computed by DLATSQR */
/* > that defines the input matrices Q_in(k) (ones on the */
/* > diagonal are not stored). See DLATSQR for more details. */
/* > */
/* > On exit: */
/* > */
/* > The array A contains an M-by-N orthonormal matrix Q_out, */
/* > i.e the columns of A are orthogonal unit vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* > T is DOUBLE PRECISION array, */
/* > dimension (LDT, N * NIRB) */
/* > where NIRB = Number_of_input_row_blocks */
/* > = MAX( 1, CEIL((M-N)/(MB-N)) ) */
/* > Let NICB = Number_of_input_col_blocks */
/* > = CEIL(N/NB) */
/* > */
/* > The upper-triangular block reflectors used to define the */
/* > input matrices Q_in(k), k=(1:NIRB*NICB). The block */
/* > reflectors are stored in compact form in NIRB block */
/* > reflector sequences. Each of the NIRB block reflector */
/* > sequences is stored in a larger NB-by-N column block of T */
/* > and consists of NICB smaller NB-by-NB upper-triangular */
/* > column blocks. See DLATSQR for more details on the format */
/* > of T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. */
/* > LDT >= fla_max(1,fla_min(NB,N)). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > The dimension of the array WORK. */
/* > LWORK >= NBLOCAL * MAX(NBLOCAL,(N-NBLOCAL)), */
/* > where NBLOCAL=MIN(NB,N). */
/* > If LWORK = -1, then a workspace query is assumed. */
/* > The routine only calculates the optimal size of the WORK */
/* > array, returns this value as the first entry of the WORK */
/* > array, and no error message related to LWORK is issued */
/* > by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* > */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup doubleOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2020, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int dorgtsqr_row_(integer *m, integer *n, integer *mb, integer *nb, doublereal *a, integer *lda, doublereal *t, integer *ldt, doublereal *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dorgtsqr_row inputs: m %" FLA_IS ", n %" FLA_IS ", mb %" FLA_IS ", nb %" FLA_IS ", lda %" FLA_IS ", ldt %" FLA_IS ", lwork  %" FLA_IS "",*m, *n, *mb, *nb, *lda, *ldt, *lwork);
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    integer lworkopt, ib_bottom__, ib, kb, mb1, mb2, m_plus_one__, num_all_row_blocks__, imb, knb;
    extern /* Subroutine */
    int dlarfb_gett_(char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
    integer jb_t__, itmp;
    doublereal dummy[1] /* was [1][1] */
    ;
    extern /* Subroutine */
    int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    logical lquery;
    integer nblocal, kb_last__;
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0 || *m < *n)
    {
        *info = -2;
    }
    else if (*mb <= *n)
    {
        *info = -3;
    }
    else if (*nb < 1)
    {
        *info = -4;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -6;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = fla_min(*nb,*n); // , expr subst
        if (*ldt < fla_max(i__1,i__2))
        {
            *info = -8;
        }
        else if (*lwork < 1 && ! lquery)
        {
            *info = -10;
        }
    }
    nblocal = fla_min(*nb,*n);
    /* Determine the workspace size. */
    if (*info == 0)
    {
        /* Computing MAX */
        i__1 = nblocal;
        i__2 = *n - nblocal; // , expr subst
        lworkopt = nblocal * fla_max(i__1,i__2);
    }
    /* Handle error in the input parameters and handle the workspace query. */
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DORGTSQR_ROW", &i__1);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (lquery)
    {
        work[1] = (doublereal) lworkopt;
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    if (fla_min(*m,*n) == 0)
    {
        work[1] = (doublereal) lworkopt;
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* (0) Set the upper-triangular part of the matrix A to zero and */
    /* its diagonal elements to one. */
    dlaset_("U", m, n, &c_b4, &c_b5, &a[a_offset], lda);
    /* KB_LAST is the column index of the last column block reflector */
    /* in the matrices T and V. */
    kb_last__ = (*n - 1) / nblocal * nblocal + 1;
    /* (1) Bottom-up loop over row blocks of A, except the top row block. */
    /* NOTE: If MB>=M, then the loop is never executed. */
    if (*mb < *m)
    {
        /* MB2 is the row blocking size for the row blocks before the */
        /* first top row block in the matrix A. IB is the row index for */
        /* the row blocks in the matrix A before the first top row block. */
        /* IB_BOTTOM is the row index for the last bottom row block */
        /* in the matrix A. JB_T is the column index of the corresponding */
        /* column block in the matrix T. */
        /* Initialize variables. */
        /* NUM_ALL_ROW_BLOCKS is the number of row blocks in the matrix A */
        /* including the first row block. */
        mb2 = *mb - *n;
        m_plus_one__ = *m + 1;
        itmp = (*m - *mb - 1) / mb2;
        ib_bottom__ = itmp * mb2 + *mb + 1;
        num_all_row_blocks__ = itmp + 2;
        jb_t__ = num_all_row_blocks__ * *n + 1;
        i__1 = *mb + 1;
        i__2 = -mb2;
        for (ib = ib_bottom__;
                i__2 < 0 ? ib >= i__1 : ib <= i__1;
                ib += i__2)
        {
            /* Determine the block size IMB for the current row block */
            /* in the matrix A. */
            /* Computing MIN */
            i__3 = m_plus_one__ - ib;
            imb = fla_min(i__3,mb2);
            /* Determine the column index JB_T for the current column block */
            /* in the matrix T. */
            jb_t__ -= *n;
            /* Apply column blocks of H in the row block from right to left. */
            /* KB is the column index of the current column block reflector */
            /* in the matrices T and V. */
            i__3 = -nblocal;
            for (kb = kb_last__;
                    i__3 < 0 ? kb >= 1 : kb <= 1;
                    kb += i__3)
            {
                /* Determine the size of the current column block KNB in */
                /* the matrices T and V. */
                /* Computing MIN */
                i__4 = nblocal;
                i__5 = *n - kb + 1; // , expr subst
                knb = fla_min(i__4,i__5);
                i__4 = *n - kb + 1;
                dlarfb_gett_("I", &imb, &i__4, &knb, &t[(jb_t__ + kb - 1) * t_dim1 + 1], ldt, &a[kb + kb * a_dim1], lda, &a[ib + kb * a_dim1], lda, &work[1], &knb);
            }
        }
    }
    /* (2) Top row block of A. */
    /* NOTE: If MB>=M, then we have only one row block of A of size M */
    /* and we work on the entire matrix A. */
    mb1 = fla_min(*mb,*m);
    /* Apply column blocks of H in the top row block from right to left. */
    /* KB is the column index of the current block reflector in */
    /* the matrices T and V. */
    i__2 = -nblocal;
    for (kb = kb_last__;
            i__2 < 0 ? kb >= 1 : kb <= 1;
            kb += i__2)
    {
        /* Determine the size of the current column block KNB in */
        /* the matrices T and V. */
        /* Computing MIN */
        i__1 = nblocal;
        i__3 = *n - kb + 1; // , expr subst
        knb = fla_min(i__1,i__3);
        if (mb1 - kb - knb + 1 == 0)
        {
            /* In SLARFB_GETT parameters, when M=0, then the matrix B */
            /* does not exist, hence we need to pass a dummy array */
            /* reference DUMMY(1,1) to B with LDDUMMY=1. */
            i__1 = *n - kb + 1;
            dlarfb_gett_("N", &c__0, &i__1, &knb, &t[kb * t_dim1 + 1], ldt, & a[kb + kb * a_dim1], lda, dummy, &c__1, &work[1], &knb);
        }
        else
        {
            i__1 = mb1 - kb - knb + 1;
            i__3 = *n - kb + 1;
            dlarfb_gett_("N", &i__1, &i__3, &knb, &t[kb * t_dim1 + 1], ldt, & a[kb + kb * a_dim1], lda, &a[kb + knb + kb * a_dim1], lda, &work[1], &knb);
        }
    }
    work[1] = (doublereal) lworkopt;
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of DORGTSQR_ROW */
}
/* dorgtsqr_row__ */

