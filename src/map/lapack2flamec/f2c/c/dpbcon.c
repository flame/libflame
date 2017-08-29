/* ../netlib/dpbcon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DPBCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DPBCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpbcon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpbcon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpbcon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, */
/* IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, KD, LDAB, N */
/* DOUBLE PRECISION ANORM, RCOND */
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
/* > DPBCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a real symmetric positive definite band matrix using the */
/* > Cholesky factorization A = U**T*U or A = L*L**T computed by DPBTRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangular factor stored in AB;
*/
/* > = 'L': Lower triangular factor stored in AB. */
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
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* > The triangular factor U or L from the Cholesky factorization */
/* > A = U**T*U or A = L*L**T of the band matrix A, stored in the */
/* > first KD+1 rows of the array. The j-th column of U or L is */
/* > stored in the j-th column of the array AB as follows: */
/* > if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
*/
/* > if UPLO ='L', AB(1+i-j,j) = L(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* > ANORM is DOUBLE PRECISION */
/* > The 1-norm (or infinity-norm) of the symmetric band matrix A. */
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
int dpbcon_(char *uplo, integer *n, integer *kd, doublereal * ab, integer *ldab, doublereal *anorm, doublereal *rcond, doublereal * work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1;
    doublereal d__1;
    /* Local variables */
    integer ix, kase;
    doublereal scale;
    extern logical lsame_(char *, char *);
    integer isave[3];
    extern /* Subroutine */
    int drscl_(integer *, doublereal *, doublereal *, integer *);
    logical upper;
    extern /* Subroutine */
    int dlacn2_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    doublereal scalel;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */
    int dlatbs_(char *, char *, char *, char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, integer *);
    doublereal scaleu;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal ainvnm;
    char normin[1];
    doublereal smlnum;
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
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*kd < 0)
    {
        *info = -3;
    }
    else if (*ldab < *kd + 1)
    {
        *info = -5;
    }
    else if (*anorm < 0.)
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DPBCON", &i__1);
        return 0;
    }
    /* Quick return if possible */
    *rcond = 0.;
    if (*n == 0)
    {
        *rcond = 1.;
        return 0;
    }
    else if (*anorm == 0.)
    {
        return 0;
    }
    smlnum = dlamch_("Safe minimum");
    /* Estimate the 1-norm of the inverse. */
    kase = 0;
    *(unsigned char *)normin = 'N';
L10:
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if (kase != 0)
    {
        if (upper)
        {
            /* Multiply by inv(U**T). */
            dlatbs_("Upper", "Transpose", "Non-unit", normin, n, kd, &ab[ ab_offset], ldab, &work[1], &scalel, &work[(*n << 1) + 1], info);
            *(unsigned char *)normin = 'Y';
            /* Multiply by inv(U). */
            dlatbs_("Upper", "No transpose", "Non-unit", normin, n, kd, &ab[ ab_offset], ldab, &work[1], &scaleu, &work[(*n << 1) + 1], info);
        }
        else
        {
            /* Multiply by inv(L). */
            dlatbs_("Lower", "No transpose", "Non-unit", normin, n, kd, &ab[ ab_offset], ldab, &work[1], &scalel, &work[(*n << 1) + 1], info);
            *(unsigned char *)normin = 'Y';
            /* Multiply by inv(L**T). */
            dlatbs_("Lower", "Transpose", "Non-unit", normin, n, kd, &ab[ ab_offset], ldab, &work[1], &scaleu, &work[(*n << 1) + 1], info);
        }
        /* Multiply by 1/SCALE if doing so will not cause overflow. */
        scale = scalel * scaleu;
        if (scale != 1.)
        {
            ix = idamax_(n, &work[1], &c__1);
            if (scale < (d__1 = work[ix], f2c_abs(d__1)) * smlnum || scale == 0.)
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
        *rcond = 1. / ainvnm / *anorm;
    }
L20:
    return 0;
    /* End of DPBCON */
}
/* dpbcon_ */
