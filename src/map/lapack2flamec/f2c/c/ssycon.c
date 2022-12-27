/* ../netlib/ssycon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SSYCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSYCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssycon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssycon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssycon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, */
/* IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N */
/* REAL ANORM, RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), IWORK( * ) */
/* REAL A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric matrix A using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by SSYTRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**T;
*/
/* > = 'L': Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The block diagonal matrix D and the multipliers used to */
/* > obtain the factor U or L as computed by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* > ANORM is REAL */
/* > The 1-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > The reciprocal of the condition number of the matrix A, */
/* > computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
/* > estimate of the 1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N) */
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
/* > \ingroup realSYcomputational */
/* ===================================================================== */
/* Subroutine */
int ssycon_(char *uplo, integer *n, real *a, integer *lda, integer *ipiv, real *anorm, real *rcond, real *work, integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"ssycon inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "",*uplo, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    integer i__, kase;
    extern logical lsame_(char *, char *);
    integer isave[3];
    logical upper;
    extern /* Subroutine */
    int slacn2_(integer *, real *, real *, integer *, real *, integer *, integer *), xerbla_(char *, integer *);
    real ainvnm;
    extern /* Subroutine */
    int ssytrs_(char *, integer *, integer *, real *, integer *, integer *, real *, integer *, integer *);
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    --iwork;
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
    else if (*anorm < 0.f)
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSYCON", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    *rcond = 0.f;
    if (*n == 0)
    {
        *rcond = 1.f;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (*anorm <= 0.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Check that the diagonal matrix D is nonsingular. */
    if (upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for (i__ = *n;
                i__ >= 1;
                --i__)
        {
            if (ipiv[i__] > 0 && a[i__ + i__ * a_dim1] == 0.f)
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return 0;
            }
            /* L10: */
        }
    }
    else
    {
        /* Lower triangular storage: examine D from top to bottom. */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            if (ipiv[i__] > 0 && a[i__ + i__ * a_dim1] == 0.f)
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return 0;
            }
            /* L20: */
        }
    }
    /* Estimate the 1-norm of the inverse. */
    kase = 0;
L30:
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if (kase != 0)
    {
        /* Multiply by inv(L*D*L**T) or inv(U*D*U**T). */
        ssytrs_(uplo, n, &c__1, &a[a_offset], lda, &ipiv[1], &work[1], n, info);
        goto L30;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if (ainvnm != 0.f)
    {
        *rcond = 1.f / ainvnm / *anorm;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of SSYCON */
}
/* ssycon_ */
