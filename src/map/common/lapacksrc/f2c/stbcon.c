/* ../netlib/stbcon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b STBCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stbcon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stbcon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stbcon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, */
/* IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORM, UPLO */
/* INTEGER INFO, KD, LDAB, N */
/* REAL RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL AB( LDAB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STBCON estimates the reciprocal of the condition number of a */
/* > triangular band matrix A, in either the 1-norm or the infinity-norm. */
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
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals or subdiagonals of the */
/* > triangular band matrix A. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is REAL array, dimension (LDAB,N) */
/* > The upper or lower triangular band matrix A, stored in the */
/* > first kd+1 rows of the array. The j-th column of A is stored */
/* > in the j-th column of the array AB as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*/
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd). */
/* > If DIAG = 'U', the diagonal elements of A are not referenced */
/* > and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
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
/* > WORK is REAL array, dimension (3*N) */
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
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int stbcon_(char *norm, char *uplo, char *diag, integer *n, integer *kd, real *ab, integer *ldab, real *rcond, real *work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1;
    real r__1;
    /* Local variables */
    integer ix, kase, kase1;
    real scale;
    extern logical lsame_(char *, char *);
    integer isave[3];
    real anorm;
    extern /* Subroutine */
    int srscl_(integer *, real *, real *, integer *);
    logical upper;
    real xnorm;
    extern /* Subroutine */
    int slacn2_(integer *, real *, real *, integer *, real *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer isamax_(integer *, real *, integer *);
    extern real slantb_(char *, char *, char *, integer *, integer *, real *, integer *, real *);
    real ainvnm;
    extern /* Subroutine */
    int slatbs_(char *, char *, char *, char *, integer *, integer *, real *, integer *, real *, real *, real *, integer *);
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --work;
    --iwork;
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
    else if (*kd < 0)
    {
        *info = -5;
    }
    else if (*ldab < *kd + 1)
    {
        *info = -7;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STBCON", &i__1);
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
    anorm = slantb_(norm, uplo, diag, n, kd, &ab[ab_offset], ldab, &work[1]);
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
        slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
        if (kase != 0)
        {
            if (kase == kase1)
            {
                /* Multiply by inv(A). */
                slatbs_(uplo, "No transpose", diag, normin, n, kd, &ab[ ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1], info) ;
            }
            else
            {
                /* Multiply by inv(A**T). */
                slatbs_(uplo, "Transpose", diag, normin, n, kd, &ab[ab_offset] , ldab, &work[1], &scale, &work[(*n << 1) + 1], info);
            }
            *(unsigned char *)normin = 'Y';
            /* Multiply by 1/SCALE if doing so will not cause overflow. */
            if (scale != 1.f)
            {
                ix = isamax_(n, &work[1], &c__1);
                xnorm = (r__1 = work[ix], f2c_abs(r__1));
                if (scale < xnorm * smlnum || scale == 0.f)
                {
                    goto L20;
                }
                srscl_(n, &scale, &work[1], &c__1);
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
    /* End of STBCON */
}
/* stbcon_ */
