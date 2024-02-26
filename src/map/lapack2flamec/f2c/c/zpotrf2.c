/* ../netlib/v3.9.0/zpotrf2.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    1.,0.
}
;
static doublereal c_b11 = -1.;
static doublereal c_b12 = 1.;
/* > \brief \b ZPOTRF2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* Definition: */
/* =========== */
/* SUBROUTINE ZPOTRF2( UPLO, N, A, LDA, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPOTRF2 computes the Cholesky factorization of a Hermitian */
/* > positive definite matrix A using the recursive algorithm. */
/* > */
/* > The factorization has the form */
/* > A = U**H * U, if UPLO = 'U', or */
/* > A = L * L**H, if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular. */
/* > */
/* > This is the recursive version of the algorithm. It divides */
/* > the matrix into four submatrices: */
/* > */
/* > [ A11 | A12 ] where A11 is n1 by n1 and A22 is n2 by n2 */
/* > A = [ -----|----- ] with n1 = n/2 */
/* > [ A21 | A22 ] n2 = n-n1 */
/* > */
/* > The subroutine calls itself to factor A11. Update and scale A21 */
/* > or A12, update A22 then call itself to factor A22. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
*/
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, if INFO = 0, the factor U or L from the Cholesky */
/* > factorization A = U**H*U or A = L*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the leading minor of order i is not */
/* > positive definite, and the factorization could not be */
/* > completed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup complex16POcomputational */
/* ===================================================================== */
/* Subroutine */
int zpotrf2_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zpotrf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "",*uplo, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer n1, n2;
    doublereal ajj;
    extern logical lsame_(char *, char *);
    integer iinfo;
    extern /* Subroutine */
    int zherk_(char *, char *, integer *, integer *, doublereal *, doublecomplex *, integer *, doublereal *, doublecomplex *, integer *);
    logical upper;
    extern /* Subroutine */
    int ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical disnan_(doublereal *);
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
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
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*n))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZPOTRF2", &i__1, (ftnlen)7);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* N=1 case */
    if (*n == 1)
    {
        /* Test for non-positive-definiteness */
        i__1 = a_dim1 + 1;
        ajj = a[i__1].r;
        if (ajj <= 0. || disnan_(&ajj))
        {
            *info = 1;
    AOCL_DTL_TRACE_LOG_EXIT
            return 0;
        }
        /* Factor */
        i__1 = a_dim1 + 1;
        d__1 = sqrt(ajj);
        a[i__1].r = d__1;
        a[i__1].i = 0.; // , expr subst
        /* Use recursive code */
    }
    else
    {
        n1 = *n / 2;
        n2 = *n - n1;
        /* Factor A11 */
        zpotrf2_(uplo, &n1, &a[a_dim1 + 1], lda, &iinfo);
        if (iinfo != 0)
        {
            *info = iinfo;
    AOCL_DTL_TRACE_LOG_EXIT
            return 0;
        }
        /* Compute the Cholesky factorization A = U**H*U */
        if (upper)
        {
            /* Update and scale A12 */
            ztrsm_("L", "U", "C", "N", &n1, &n2, &c_b1, &a[a_dim1 + 1], lda, & a[(n1 + 1) * a_dim1 + 1], lda);
            /* Update and factor A22 */
            zherk_(uplo, "C", &n2, &n1, &c_b11, &a[(n1 + 1) * a_dim1 + 1], lda, &c_b12, &a[n1 + 1 + (n1 + 1) * a_dim1], lda);
            zpotrf2_(uplo, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &iinfo);
            if (iinfo != 0)
            {
                *info = iinfo + n1;
    AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
            /* Compute the Cholesky factorization A = L*L**H */
        }
        else
        {
            /* Update and scale A21 */
            ztrsm_("R", "L", "C", "N", &n2, &n1, &c_b1, &a[a_dim1 + 1], lda, & a[n1 + 1 + a_dim1], lda);
            /* Update and factor A22 */
            zherk_(uplo, "N", &n2, &n1, &c_b11, &a[n1 + 1 + a_dim1], lda, & c_b12, &a[n1 + 1 + (n1 + 1) * a_dim1], lda);
            zpotrf2_(uplo, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &iinfo);
            if (iinfo != 0)
            {
                *info = iinfo + n1;
    AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZPOTRF2 */
}
/* zpotrf2_ */
