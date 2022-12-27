/* sgetsqrhrt.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SGETSQRHRT */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGETSQRHRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetsqr hrt.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetsqr hrt.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetsqr hrt.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK, */
/* $ LWORK, INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDT, LWORK, M, N, NB1, NB2, MB1 */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), T( LDT, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGETSQRHRT computes a NB2-sized column blocked QR-factorization */
/* > of a complex M-by-N matrix A with M >= N, */
/* > */
/* > A = Q * R. */
/* > */
/* > The routine uses internally a NB1-sized column blocked and MB1-sized */
/* > row blocked TSQR-factorization and perfors the reconstruction */
/* > of the Householder vectors from the TSQR output. The routine also */
/* > converts the R_tsqr factor from the TSQR-factorization output into */
/* > the R factor that corresponds to the Householder QR-factorization, */
/* > */
/* > A = Q_tsqr * R_tsqr = Q * R. */
/* > */
/* > The output Q and R factors are stored in the same format as in SGEQRT */
/* > (Q is in blocked compact WY-representation). See the documentation */
/* > of SGEQRT for more details on the format. */
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
/* > \param[in] MB1 */
/* > \verbatim */
/* > MB1 is INTEGER */
/* > The row block size to be used in the blocked TSQR. */
/* > MB1 > N. */
/* > \endverbatim */
/* > */
/* > \param[in] NB1 */
/* > \verbatim */
/* > NB1 is INTEGER */
/* > The column block size to be used in the blocked TSQR. */
/* > N >= NB1 >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] NB2 */
/* > \verbatim */
/* > NB2 is INTEGER */
/* > The block size to be used in the blocked QR that is */
/* > output. NB2 >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > */
/* > On entry: an M-by-N matrix A. */
/* > */
/* > On exit: */
/* > a) the elements on and above the diagonal */
/* > of the array contain the N-by-N upper-triangular */
/* > matrix R corresponding to the Householder QR;
*/
/* > b) the elements below the diagonal represent Q by */
/* > the columns of blocked V (compact WY-representation). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* > T is REAL array, dimension (LDT,N)) */
/* > The upper triangular block reflectors stored in compact form */
/* > as a sequence of upper triangular blocks. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= NB2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > (workspace) REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > The dimension of the array WORK. */
/* > LWORK >= MAX( LWT + LW1, MAX( LWT+N*N+LW2, LWT+N*N+N ) ), */
/* > where */
/* > NUM_ALL_ROW_BLOCKS = CEIL((M-N)/(MB1-N)), */
/* > NB1LOCAL = MIN(NB1,N). */
/* > LWT = NUM_ALL_ROW_BLOCKS * N * NB1LOCAL, */
/* > LW1 = NB1LOCAL * N, */
/* > LW2 = NB1LOCAL * MAX( NB1LOCAL, ( N - NB1LOCAL ) ), */
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
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup singleOTHERcomputational */
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
int sgetsqrhrt_(integer *m, integer *n, integer *mb1, integer *nb1, integer *nb2, real *a, integer *lda, real *t, integer * ldt, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    /* Local variables */
    integer lworkopt, i__, j;
    extern /* Subroutine */
    int sorhr_col_(integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *);
    integer lw1, lw2, num_all_row_blocks__, lwt, ldwt, iinfo;
    extern /* Subroutine */
    int sorgtsqr_row_(integer *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *), scopy_(integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *);
    logical lquery;
    extern /* Subroutine */
    int slatsqr_(integer *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *);
    integer nb1local, nb2local;
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
    else if (*mb1 <= *n)
    {
        *info = -3;
    }
    else if (*nb1 < 1)
    {
        *info = -4;
    }
    else if (*nb2 < 1)
    {
        *info = -5;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -7;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = fla_min(*nb2,*n); // , expr subst
        if (*ldt < fla_max(i__1,i__2))
        {
            *info = -9;
        }
        else
        {
            /* Test the input LWORK for the dimension of the array WORK. */
            /* This workspace is used to store array: */
            /* a) Matrix T and WORK for SLATSQR;
            */
            /* b) N-by-N upper-triangular factor R_tsqr;
            */
            /* c) Matrix T and array WORK for SORGTSQR_ROW;
            */
            /* d) Diagonal D for SORHR_COL. */
            if (*lwork < *n * *n + 1 && ! lquery)
            {
                *info = -11;
            }
            else
            {
                /* Set block size for column blocks */
                nb1local = fla_min(*nb1,*n);
                /* Computing MAX */
                r__1 = 1.f;
                r__2 = ceiling_f90_((real) (*m - *n) / (real) (*mb1 - *n)); // , expr subst
                num_all_row_blocks__ = fla_max(r__1,r__2);
                /* Length and leading dimension of WORK array to place */
                /* T array in TSQR. */
                lwt = num_all_row_blocks__ * *n * nb1local;
                ldwt = nb1local;
                /* Length of TSQR work array */
                lw1 = nb1local * *n;
                /* Length of SORGTSQR_ROW work array. */
                /* Computing MAX */
                i__1 = nb1local;
                i__2 = *n - nb1local; // , expr subst
                lw2 = nb1local * fla_max(i__1,i__2);
                /* Computing MAX */
                /* Computing MAX */
                i__3 = lwt + *n * *n + lw2;
                i__4 = lwt + *n * *n + *n; // , expr subst
                i__1 = lwt + lw1;
                i__2 = fla_max(i__3,i__4); // , expr subst
                lworkopt = fla_max(i__1,i__2);
                if (*lwork < fla_max(1,lworkopt) && ! lquery)
                {
                    *info = -11;
                }
            }
        }
    }
    /* Handle error in the input parameters and return workspace query. */
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGETSQRHRT", &i__1);
        return 0;
    }
    else if (lquery)
    {
        work[1] = (real) lworkopt;
        return 0;
    }
    /* Quick return if possible */
    if (fla_min(*m,*n) == 0)
    {
        work[1] = (real) lworkopt;
        return 0;
    }
    nb2local = fla_min(*nb2,*n);
    /* (1) Perform TSQR-factorization of the M-by-N matrix A. */
    slatsqr_(m, n, mb1, &nb1local, &a[a_offset], lda, &work[1], &ldwt, &work[ lwt + 1], &lw1, &iinfo);
    /* (2) Copy the factor R_tsqr stored in the upper-triangular part */
    /* of A into the square matrix in the work array */
    /* WORK(LWT+1:LWT+N*N) column-by-column. */
    i__1 = *n;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        scopy_(&j, &a[j * a_dim1 + 1], &c__1, &work[lwt + *n * (j - 1) + 1], & c__1);
    }
    /* (3) Generate a M-by-N matrix Q with orthonormal columns from */
    /* the result stored below the diagonal in the array A in place. */
    sorgtsqr_row_(m, n, mb1, &nb1local, &a[a_offset], lda, &work[1], &ldwt, & work[lwt + *n * *n + 1], &lw2, &iinfo);
    /* (4) Perform the reconstruction of Householder vectors from */
    /* the matrix Q (stored in A) in place. */
    sorhr_col_(m, n, &nb2local, &a[a_offset], lda, &t[t_offset], ldt, &work[ lwt + *n * *n + 1], &iinfo);
    /* (5) Copy the factor R_tsqr stored in the square matrix in the */
    /* work array WORK(LWT+1:LWT+N*N) into the upper-triangular */
    /* part of A. */
    /* (6) Compute from R_tsqr the factor R_hr corresponding to */
    /* the reconstructed Householder vectors, i.e. R_hr = S * R_tsqr. */
    /* This multiplication by the sign matrix S on the left means */
    /* changing the sign of I-th row of the matrix R_tsqr according */
    /* to sign of the I-th diagonal element DIAG(I) of the matrix S. */
    /* DIAG is stored in WORK( LWT+N*N+1 ) from the SORHR_COL output. */
    /* (5) and (6) can be combined in a single loop, so the rows in A */
    /* are accessed only once. */
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if (work[lwt + *n * *n + i__] == -1.f)
        {
            i__2 = *n;
            for (j = i__;
                    j <= i__2;
                    ++j)
            {
                a[i__ + j * a_dim1] = work[lwt + *n * (j - 1) + i__] * -1.f;
            }
        }
        else
        {
            i__2 = *n - i__ + 1;
            scopy_(&i__2, &work[lwt + *n * (i__ - 1) + i__], n, &a[i__ + i__ * a_dim1], lda);
        }
    }
    work[1] = (real) lworkopt;
    return 0;
    /* End of SGETSQRHRT */
}
/* sgetsqrhrt_ */

