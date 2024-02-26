/* ../netlib/zhpcon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZHPCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpcon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpcon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpcon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, N */
/* DOUBLE PRECISION ANORM, RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 AP( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPCON estimates the reciprocal of the condition number of a complex */
/* > Hermitian packed matrix A using the factorization A = U*D*U**H or */
/* > A = L*D*L**H computed by ZHPTRF. */
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
/* > = 'U': Upper triangular, form is A = U*D*U**H;
*/
/* > = 'L': Lower triangular, form is A = L*D*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > The block diagonal matrix D and the multipliers used to */
/* > obtain the factor U or L as computed by ZHPTRF, stored as a */
/* > packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by ZHPTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* > ANORM is DOUBLE PRECISION */
/* > The 1-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is DOUBLE PRECISION */
/* > The reciprocal of the condition number of the matrix A, */
/* > computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an */
/* > estimate of the 1-norm of inv(A) computed in this routine. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (2*N) */
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
/* > \ingroup complex16OTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int zhpcon_(char *uplo, integer *n, doublecomplex *ap, integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex * work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhpcon inputs: uplo %c, n %" FLA_IS "",*uplo, *n);
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    integer i__, ip, kase;
    extern logical lsame_(char *, char *);
    integer isave[3];
    logical upper;
    extern /* Subroutine */
    int zlacn2_(integer *, doublecomplex *, doublecomplex *, doublereal *, integer *, integer *), xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    doublereal ainvnm;
    extern /* Subroutine */
    int zhptrs_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --work;
    --ipiv;
    --ap;
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
    else if (*anorm < 0.)
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHPCON", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    *rcond = 0.;
    if (*n == 0)
    {
        *rcond = 1.;
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (*anorm <= 0.)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Check that the diagonal matrix D is nonsingular. */
    if (upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        ip = *n * (*n + 1) / 2;
        for (i__ = *n;
                i__ >= 1;
                --i__)
        {
            i__1 = ip;
            if (ipiv[i__] > 0 && (ap[i__1].r == 0. && ap[i__1].i == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
            ip -= i__;
            /* L10: */
        }
    }
    else
    {
        /* Lower triangular storage: examine D from top to bottom. */
        ip = 1;
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = ip;
            if (ipiv[i__] > 0 && (ap[i__2].r == 0. && ap[i__2].i == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
            ip = ip + *n - i__ + 1;
            /* L20: */
        }
    }
    /* Estimate the 1-norm of the inverse. */
    kase = 0;
L30:
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
    if (kase != 0)
    {
        /* Multiply by inv(L*D*L**H) or inv(U*D*U**H). */
        zhptrs_(uplo, n, &c__1, &ap[1], &ipiv[1], &work[1], n, info);
        goto L30;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if (ainvnm != 0.)
    {
        *rcond = 1. / ainvnm / *anorm;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZHPCON */
}
/* zhpcon_ */
