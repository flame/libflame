/* ../netlib/ztpcon.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZTPCON */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZTPCON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpcon. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpcon. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpcon. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZTPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, RWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORM, UPLO */
/* INTEGER INFO, N */
/* DOUBLE PRECISION RCOND */
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
/* > ZTPCON estimates the reciprocal of the condition number of a packed */
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
/* > \param[in] AP */
/* > \verbatim */
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > The upper or lower triangular matrix A, packed columnwise in */
/* > a linear array. The j-th column of A is stored in the array */
/* > AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > If DIAG = 'U', the diagonal elements of A are not referenced */
/* > and are assumed to be 1. */
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
int ztpcon_(char *norm, char *uplo, char *diag, integer *n, doublecomplex *ap, doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    integer ix, kase, kase1;
    doublereal scale;
    extern logical lsame_(char *, char *);
    integer isave[3];
    doublereal anorm;
    logical upper;
    doublereal xnorm;
    extern /* Subroutine */
    int zlacn2_(integer *, doublecomplex *, doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal ainvnm;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    logical onenrm;
    extern /* Subroutine */
    int zdrscl_(integer *, doublereal *, doublecomplex *, integer *);
    char normin[1];
    extern doublereal zlantp_(char *, char *, char *, integer *, doublecomplex *, doublereal *);
    doublereal smlnum;
    logical nounit;
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
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZTPCON", &i__1);
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
    anorm = zlantp_(norm, uplo, diag, n, &ap[1], &rwork[1]);
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
        zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
        if (kase != 0)
        {
            if (kase == kase1)
            {
                /* Multiply by inv(A). */
                zlatps_(uplo, "No transpose", diag, normin, n, &ap[1], &work[ 1], &scale, &rwork[1], info);
            }
            else
            {
                /* Multiply by inv(A**H). */
                zlatps_(uplo, "Conjugate transpose", diag, normin, n, &ap[1], &work[1], &scale, &rwork[1], info);
            }
            *(unsigned char *)normin = 'Y';
            /* Multiply by 1/SCALE if doing so will not cause overflow. */
            if (scale != 1.)
            {
                ix = izamax_(n, &work[1], &c__1);
                i__1 = ix;
                xnorm = (d__1 = work[i__1].r, f2c_abs(d__1)) + (d__2 = d_imag(& work[ix]), f2c_abs(d__2));
                if (scale < xnorm * smlnum || scale == 0.)
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
            *rcond = 1. / anorm / ainvnm;
        }
    }
L20:
    return 0;
    /* End of ZTPCON */
}
/* ztpcon_ */
