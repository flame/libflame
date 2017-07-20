/* ../netlib/zppcon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZPPCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZPPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zppcon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zppcon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zppcon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZPPCON( UPLO, N, AP, ANORM, RCOND, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, N */
/* DOUBLE PRECISION ANORM, RCOND */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION RWORK( * ) */
/* COMPLEX*16 AP( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPPCON estimates the reciprocal of the condition number (in the */
/* > 1-norm) of a complex Hermitian positive definite packed matrix using */
/* > the Cholesky factorization A = U**H*U or A = L*L**H computed by */
/* > ZPPTRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))). */
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
/* > \param[in] AP */
/* > \verbatim */
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > The triangular factor U or L from the Cholesky factorization */
/* > A = U**H*U or A = L*L**H, packed columnwise in a linear */
/* > array. The j-th column of U or L is stored in the array AP */
/* > as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* > ANORM is DOUBLE PRECISION */
/* > The 1-norm (or infinity-norm) of the Hermitian matrix A. */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension (N) */
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
int zppcon_(char *uplo, integer *n, doublecomplex *ap, doublereal *anorm, doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    integer ix, kase;
    doublereal scale;
    extern logical lsame_(char *, char *);
    integer isave[3];
    logical upper;
    extern /* Subroutine */
    int zlacn2_(integer *, doublecomplex *, doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    doublereal scalel, scaleu;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal ainvnm;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */
    int zdrscl_(integer *, doublereal *, doublecomplex *, integer *);
    char normin[1];
    doublereal smlnum;
    extern /* Subroutine */
    int zlatps_(char *, char *, char *, char *, integer *, doublecomplex *, doublecomplex *, doublereal *, doublereal *, integer *);
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
    --rwork;
    --work;
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
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZPPCON", &i__1);
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
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
    if (kase != 0)
    {
        if (upper)
        {
            /* Multiply by inv(U**H). */
            zlatps_("Upper", "Conjugate transpose", "Non-unit", normin, n, & ap[1], &work[1], &scalel, &rwork[1], info);
            *(unsigned char *)normin = 'Y';
            /* Multiply by inv(U). */
            zlatps_("Upper", "No transpose", "Non-unit", normin, n, &ap[1], & work[1], &scaleu, &rwork[1], info);
        }
        else
        {
            /* Multiply by inv(L). */
            zlatps_("Lower", "No transpose", "Non-unit", normin, n, &ap[1], & work[1], &scalel, &rwork[1], info);
            *(unsigned char *)normin = 'Y';
            /* Multiply by inv(L**H). */
            zlatps_("Lower", "Conjugate transpose", "Non-unit", normin, n, & ap[1], &work[1], &scaleu, &rwork[1], info);
        }
        /* Multiply by 1/SCALE if doing so will not cause overflow. */
        scale = scalel * scaleu;
        if (scale != 1.)
        {
            ix = izamax_(n, &work[1], &c__1);
            i__1 = ix;
            if (scale < ((d__1 = work[i__1].r, f2c_abs(d__1)) + (d__2 = d_imag(& work[ix]), f2c_abs(d__2))) * smlnum || scale == 0.)
            {
                goto L20;
            }
            zdrscl_(n, &scale, &work[1], &c__1);
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
    /* End of ZPPCON */
}
/* zppcon_ */
