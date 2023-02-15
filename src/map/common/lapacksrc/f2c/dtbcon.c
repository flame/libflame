/* ../netlib/dtbcon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DTBCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DTBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtbcon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtbcon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtbcon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, */
/* IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORM, UPLO */
/* INTEGER INFO, KD, LDAB, N */
/* DOUBLE PRECISION RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION AB( LDAB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTBCON estimates the reciprocal of the condition number of a */
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
/* > AB is DOUBLE PRECISION array, dimension (LDAB,N) */
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
/* > RCOND is DOUBLE PRECISION */
/* > The reciprocal of the condition number of the matrix A, */
/* > computed as RCOND = 1/(norm(A) * norm(inv(A))). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (3*N) */
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
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int dtbcon_(char *norm, char *uplo, char *diag, integer *n, integer *kd, doublereal *ab, integer *ldab, doublereal *rcond, doublereal *work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1;
    doublereal d__1;
    /* Local variables */
    integer ix, kase, kase1;
    doublereal scale;
    extern logical lsame_(char *, char *);
    integer isave[3];
    extern /* Subroutine */
    int drscl_(integer *, doublereal *, doublereal *, integer *);
    doublereal anorm;
    logical upper;
    doublereal xnorm;
    extern /* Subroutine */
    int dlacn2_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern doublereal dlantb_(char *, char *, char *, integer *, integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */
    int dlatbs_(char *, char *, char *, char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    doublereal ainvnm;
    logical onenrm;
    char normin[1];
    doublereal smlnum;
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
        xerbla_("DTBCON", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        *rcond = 1.;
        return 0;
    }
    *rcond = 0.;
    smlnum = dlamch_("Safe minimum") * (doublereal) max(1,*n);
    /* Compute the norm of the triangular matrix A. */
    anorm = dlantb_(norm, uplo, diag, n, kd, &ab[ab_offset], ldab, &work[1]);
    /* Continue only if ANORM > 0. */
    if (anorm > 0.)
    {
        /* Estimate the norm of the inverse of A. */
        ainvnm = 0.;
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
        dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
        if (kase != 0)
        {
            if (kase == kase1)
            {
                /* Multiply by inv(A). */
                dlatbs_(uplo, "No transpose", diag, normin, n, kd, &ab[ ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1], info) ;
            }
            else
            {
                /* Multiply by inv(A**T). */
                dlatbs_(uplo, "Transpose", diag, normin, n, kd, &ab[ab_offset] , ldab, &work[1], &scale, &work[(*n << 1) + 1], info);
            }
            *(unsigned char *)normin = 'Y';
            /* Multiply by 1/SCALE if doing so will not cause overflow. */
            if (scale != 1.)
            {
                ix = idamax_(n, &work[1], &c__1);
                xnorm = (d__1 = work[ix], f2c_abs(d__1));
                if (scale < xnorm * smlnum || scale == 0.)
                {
                    goto L20;
                }
                drscl_(n, &scale, &work[1], &c__1);
            }
            goto L10;
        }
        /* Compute the estimate of the reciprocal condition number. */
        if (ainvnm != 0.)
        {
            *rcond = 1. / anorm / ainvnm;
        }
    }
L20:
    return 0;
    /* End of DTBCON */
}
/* dtbcon_ */
