/* ../netlib/v3.9.0/zlaunhr_col_getrfnp.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    1.,0.
}
;
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b ZLAUNHR_COL_GETRFNP */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAUNHR_COL_GETRFNP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaunhr _col_getrfnp.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaunhr _col_getrfnp.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaunhr _col_getrfnp.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAUNHR_COL_GETRFNP( M, N, A, LDA, D, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), D( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAUNHR_COL_GETRFNP computes the modified LU factorization without */
/* > pivoting of a complex general M-by-N matrix A. The factorization has */
/* > the form: */
/* > */
/* > A - S = L * U, */
/* > */
/* > where: */
/* > S is a m-by-n diagonal sign matrix with the diagonal D, so that */
/* > D(i) = S(i,i), 1 <= i <= fla_min(M,N). The diagonal D is constructed */
/* > as D(i)=-SIGN(A(i,i)), where A(i,i) is the value after performing */
/* > i-1 steps of Gaussian elimination. This means that the diagonal */
/* > element at each step of "modified" Gaussian elimination is */
/* > at least one in absolute value (so that division-by-zero not */
/* > not possible during the division by the diagonal element);
*/
/* > */
/* > L is a M-by-N lower triangular matrix with unit diagonal elements */
/* > (lower trapezoidal if M > N);
*/
/* > */
/* > and U is a M-by-N upper triangular matrix */
/* > (upper trapezoidal if M < N). */
/* > */
/* > This routine is an auxiliary routine used in the Householder */
/* > reconstruction routine ZUNHR_COL. In ZUNHR_COL, this routine is */
/* > applied to an M-by-N matrix A with orthonormal columns, where each */
/* > element is bounded by one in absolute value. With the choice of */
/* > the matrix S above, one can show that the diagonal element at each */
/* > step of Gaussian elimination is the largest (in absolute value) in */
/* > the column on or below the diagonal, so that no pivoting is required */
/* > for numerical stability [1]. */
/* > */
/* > For more details on the Householder reconstruction algorithm, */
/* > including the modified LU factorization, see [1]. */
/* > */
/* > This is the blocked right-looking version of the algorithm, */
/* > calling Level 3 BLAS to update the submatrix. To factorize a block, */
/* > this routine calls the recursive routine ZLAUNHR_COL_GETRFNP2. */
/* > */
/* > [1] "Reconstructing Householder vectors from tall-skinny QR", */
/* > G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen, */
/* > E. Solomonik, J. Parallel Distrib. Comput., */
/* > vol. 85, pp. 3-31, 2015. */
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
/* > On entry, the M-by-N matrix to be factored. */
/* > On exit, the factors L and U from the factorization */
/* > A-S=L*U;
the unit diagonal elements of L are not stored. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is COMPLEX*16 array, dimension fla_min(M,N) */
/* > The diagonal elements of the diagonal M-by-N sign matrix S, */
/* > D(i) = S(i,i), where 1 <= i <= fla_min(M,N). The elements can be */
/* > only ( +1.0, 0.0 ) or (-1.0, 0.0 ). */
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
/* > \date November 2019 */
/* > \ingroup complex16GEcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2019, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int zlaunhr_col_getrfnp_(integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *d__, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlaunhr_col_getrfnp inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "",*m, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;
    /* Local variables */
    integer j, jb, nb;
    extern /* Subroutine */
    int zlaunhr_col_getrfnp2_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    integer iinfo;
    extern /* Subroutine */
    int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *), ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    /* -- LAPACK computational routine (version 3.9.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2019 */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
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
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZLAUNHR_COL_GETRFNP", &i__1);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    if (fla_min(*m,*n) == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Determine the block size for this environment. */
    nb = ilaenv_(&c__1, "ZLAUNHR_COL_GETRFNP", " ", m, n, &c_n1, &c_n1);
    if (nb <= 1 || nb >= fla_min(*m,*n))
    {
        /* Use unblocked code. */
        zlaunhr_col_getrfnp2_(m, n, &a[a_offset], lda, &d__[1], info);
    }
    else
    {
        /* Use blocked code. */
        i__1 = fla_min(*m,*n);
        i__2 = nb;
        for (j = 1;
                i__2 < 0 ? j >= i__1 : j <= i__1;
                j += i__2)
        {
            /* Computing MIN */
            i__3 = fla_min(*m,*n) - j + 1;
            jb = fla_min(i__3,nb);
            /* Factor diagonal and subdiagonal blocks. */
            i__3 = *m - j + 1;
            zlaunhr_col_getrfnp2_(&i__3, &jb, &a[j + j * a_dim1], lda, &d__[ j], &iinfo);
            if (j + jb <= *n)
            {
                /* Compute block row of U. */
                i__3 = *n - j - jb + 1;
                ztrsm_("Left", "Lower", "No transpose", "Unit", &jb, &i__3, & c_b1, &a[j + j * a_dim1], lda, &a[j + (j + jb) * a_dim1], lda);
                if (j + jb <= *m)
                {
                    /* Update trailing submatrix. */
                    i__3 = *m - j - jb + 1;
                    i__4 = *n - j - jb + 1;
                    z__1.r = -1.;
                    z__1.i = -0.; // , expr subst
                    zgemm_("No transpose", "No transpose", &i__3, &i__4, &jb, &z__1, &a[j + jb + j * a_dim1], lda, &a[j + (j + jb) * a_dim1], lda, &c_b1, &a[j + jb + (j + jb) * a_dim1], lda);
                }
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZLAUNHR_COL_GETRFNP */
}
/* zlaunhr_col_getrfnp__ */

