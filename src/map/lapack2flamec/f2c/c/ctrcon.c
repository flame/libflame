/* ../netlib/ctrcon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CTRCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CTRCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrcon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrcon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrcon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, */
/* RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORM, UPLO */
/* INTEGER INFO, LDA, N */
/* REAL RCOND */
/* .. */
/* .. Array Arguments .. */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTRCON estimates the reciprocal of the condition number of a */
/* > triangular matrix A, in either the 1-norm or the infinity-norm. */
/* > */
/* > The norm of A is computed and an estimate is obtained for */
/* > norm(inv(A)), then the reciprocal of the condition number is */
/* > computed as */
/* > RCOND = 1 / ( norm(A) * norm(inv(A)) ). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] NORM */
/* > \verbatim */
/* > NORM is CHARACTER*1 */
/* > Specifies whether the 1-norm condition number or the */
/* > infinity-norm condition number is required: */
/* > = '1' or 'O': 1-norm;
*/
/* > = 'I': Infinity-norm. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': A is upper triangular;
*/
/* > = 'L': A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* > DIAG is CHARACTER*1 */
/* > = 'N': A is non-unit triangular;
*/
/* > = 'U': A is unit triangular. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > The triangular matrix A. If UPLO = 'U', the leading N-by-N */
/* > upper triangular part of the array A contains the upper */
/* > triangular matrix, and the strictly lower triangular part of */
/* > A is not referenced. If UPLO = 'L', the leading N-by-N lower */
/* > triangular part of the array A contains the lower triangular */
/* > matrix, and the strictly upper triangular part of A is not */
/* > referenced. If DIAG = 'U', the diagonal elements of A are */
/* > also not referenced and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > The reciprocal of the condition number of the matrix A, */
/* > computed as RCOND = 1/(norm(A) * norm(inv(A))). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (N) */
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
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int ctrcon_(char *norm, char *uplo, char *diag, integer *n, complex *a, integer *lda, real *rcond, complex *work, real *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    real r__1, r__2;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    integer ix, kase, kase1;
    real scale;
    extern logical lsame_(char *, char *);
    integer isave[3];
    real anorm;
    logical upper;
    extern /* Subroutine */
    int clacn2_(integer *, complex *, complex *, real *, integer *, integer *);
    real xnorm;
    extern integer icamax_(integer *, complex *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern real clantr_(char *, char *, char *, integer *, integer *, complex *, integer *, real *);
    real ainvnm;
    extern /* Subroutine */
    int clatrs_(char *, char *, char *, char *, integer *, complex *, integer *, complex *, real *, real *, integer *), csrscl_(integer *, real *, complex *, integer *);
    logical onenrm;
    char normin[1];
    real smlnum;
    logical nounit;
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --rwork;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    nounit = lsame_(diag, "N");
    if (! onenrm && ! lsame_(norm, "I"))
    {
        *info = -1;
    }
    else if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -2;
    }
    else if (! nounit && ! lsame_(diag, "U"))
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*lda < max(1,*n))
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CTRCON", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        *rcond = 1.f;
        return 0;
    }
    *rcond = 0.f;
    smlnum = slamch_("Safe minimum") * (real) max(1,*n);
    /* Compute the norm of the triangular matrix A. */
    anorm = clantr_(norm, uplo, diag, n, n, &a[a_offset], lda, &rwork[1]);
    /* Continue only if ANORM > 0. */
    if (anorm > 0.f)
    {
        /* Estimate the norm of the inverse of A. */
        ainvnm = 0.f;
        *(unsigned char *)normin = 'N';
        if (onenrm)
        {
            kase1 = 1;
        }
        else
        {
            kase1 = 2;
        }
        kase = 0;
L10:
        clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
        if (kase != 0)
        {
            if (kase == kase1)
            {
                /* Multiply by inv(A). */
                clatrs_(uplo, "No transpose", diag, normin, n, &a[a_offset], lda, &work[1], &scale, &rwork[1], info);
            }
            else
            {
                /* Multiply by inv(A**H). */
                clatrs_(uplo, "Conjugate transpose", diag, normin, n, &a[ a_offset], lda, &work[1], &scale, &rwork[1], info);
            }
            *(unsigned char *)normin = 'Y';
            /* Multiply by 1/SCALE if doing so will not cause overflow. */
            if (scale != 1.f)
            {
                ix = icamax_(n, &work[1], &c__1);
                i__1 = ix;
                xnorm = (r__1 = work[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(& work[ix]), f2c_abs(r__2));
                if (scale < xnorm * smlnum || scale == 0.f)
                {
                    goto L20;
                }
                csrscl_(n, &scale, &work[1], &c__1);
            }
            goto L10;
        }
        /* Compute the estimate of the reciprocal condition number. */
        if (ainvnm != 0.f)
        {
            *rcond = 1.f / anorm / ainvnm;
        }
    }
L20:
    return 0;
    /* End of CTRCON */
}
/* ctrcon_ */
