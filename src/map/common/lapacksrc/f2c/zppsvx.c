/* ../netlib/zppsvx.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief <b> ZPPSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZPPSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zppsvx. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zppsvx. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zppsvx. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, */
/* X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED, FACT, UPLO */
/* INTEGER INFO, LDB, LDX, N, NRHS */
/* DOUBLE PRECISION RCOND */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION BERR( * ), FERR( * ), RWORK( * ), S( * ) */
/* COMPLEX*16 AFP( * ), AP( * ), B( LDB, * ), WORK( * ), */
/* $ X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPPSVX uses the Cholesky factorization A = U**H * U or A = L * L**H to */
/* > compute the solution to a complex system of linear equations */
/* > A * X = B, */
/* > where A is an N-by-N Hermitian positive definite matrix stored in */
/* > packed format and X and B are N-by-NRHS matrices. */
/* > */
/* > Error bounds on the solution and a condition estimate are also */
/* > provided. */
/* > \endverbatim */
/* > \par Description: */
/* ================= */
/* > */
/* > \verbatim */
/* > */
/* > The following steps are performed: */
/* > */
/* > 1. If FACT = 'E', real scaling factors are computed to equilibrate */
/* > the system: */
/* > diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B */
/* > Whether or not the system will be equilibrated depends on the */
/* > scaling of the matrix A, but if equilibration is used, A is */
/* > overwritten by diag(S)*A*diag(S) and B by diag(S)*B. */
/* > */
/* > 2. If FACT = 'N' or 'E', the Cholesky decomposition is used to */
/* > factor the matrix A (after equilibration if FACT = 'E') as */
/* > A = U**H * U , if UPLO = 'U', or */
/* > A = L * L**H, if UPLO = 'L', */
/* > where U is an upper triangular matrix, L is a lower triangular */
/* > matrix, and **H indicates conjugate transpose. */
/* > */
/* > 3. If the leading i-by-i principal minor is not positive definite, */
/* > then the routine returns with INFO = i. Otherwise, the factored */
/* > form of A is used to estimate the condition number of the matrix */
/* > A. If the reciprocal of the condition number is less than machine */
/* > precision, INFO = N+1 is returned as a warning, but the routine */
/* > still goes on to solve for X and compute error bounds as */
/* > described below. */
/* > */
/* > 4. The system of equations is solved for X using the factored form */
/* > of A. */
/* > */
/* > 5. Iterative refinement is applied to improve the computed solution */
/* > matrix and calculate error bounds and backward error estimates */
/* > for it. */
/* > */
/* > 6. If equilibration was used, the matrix X is premultiplied by */
/* > diag(S) so that it solves the original system before */
/* > equilibration. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] FACT */
/* > \verbatim */
/* > FACT is CHARACTER*1 */
/* > Specifies whether or not the factored form of the matrix A is */
/* > supplied on entry, and if not, whether the matrix A should be */
/* > equilibrated before it is factored. */
/* > = 'F': On entry, AFP contains the factored form of A. */
/* > If EQUED = 'Y', the matrix A has been equilibrated */
/* > with scaling factors given by S. AP and AFP will not */
/* > be modified. */
/* > = 'N': The matrix A will be copied to AFP and factored. */
/* > = 'E': The matrix A will be equilibrated if necessary, then */
/* > copied to AFP and factored. */
/* > \endverbatim */
/* > */
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
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the Hermitian matrix */
/* > A, packed columnwise in a linear array, except if FACT = 'F' */
/* > and EQUED = 'Y', then A must contain the equilibrated matrix */
/* > diag(S)*A*diag(S). The j-th column of A is stored in the */
/* > array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > See below for further details. A is not modified if */
/* > FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit. */
/* > */
/* > On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by */
/* > diag(S)*A*diag(S). */
/* > \endverbatim */
/* > */
/* > \param[in,out] AFP */
/* > \verbatim */
/* > AFP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > If FACT = 'F', then AFP is an input argument and on entry */
/* > contains the triangular factor U or L from the Cholesky */
/* > factorization A = U**H*U or A = L*L**H, in the same storage */
/* > format as A. If EQUED .ne. 'N', then AFP is the factored */
/* > form of the equilibrated matrix A. */
/* > */
/* > If FACT = 'N', then AFP is an output argument and on exit */
/* > returns the triangular factor U or L from the Cholesky */
/* > factorization A = U**H * U or A = L * L**H of the original */
/* > matrix A. */
/* > */
/* > If FACT = 'E', then AFP is an output argument and on exit */
/* > returns the triangular factor U or L from the Cholesky */
/* > factorization A = U**H * U or A = L * L**H of the equilibrated */
/* > matrix A (see the description of AP for the form of the */
/* > equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in,out] EQUED */
/* > \verbatim */
/* > EQUED is CHARACTER*1 */
/* > Specifies the form of equilibration that was done. */
/* > = 'N': No equilibration (always true if FACT = 'N'). */
/* > = 'Y': Equilibration was done, i.e., A has been replaced by */
/* > diag(S) * A * diag(S). */
/* > EQUED is an input argument if FACT = 'F';
otherwise, it is an */
/* > output argument. */
/* > \endverbatim */
/* > */
/* > \param[in,out] S */
/* > \verbatim */
/* > S is DOUBLE PRECISION array, dimension (N) */
/* > The scale factors for A;
not accessed if EQUED = 'N'. S is */
/* > an input argument if FACT = 'F';
otherwise, S is an output */
/* > argument. If FACT = 'F' and EQUED = 'Y', each element of S */
/* > must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the N-by-NRHS right hand side matrix B. */
/* > On exit, if EQUED = 'N', B is not modified;
if EQUED = 'Y', */
/* > B is overwritten by diag(S) * B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (LDX,NRHS) */
/* > If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to */
/* > the original system of equations. Note that if EQUED = 'Y', */
/* > A and B are modified on exit, and the solution to the */
/* > equilibrated system is inv(diag(S))*X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is DOUBLE PRECISION */
/* > The estimate of the reciprocal condition number of the matrix */
/* > A after equilibration (if done). If RCOND is less than the */
/* > machine precision (in particular, if RCOND = 0), the matrix */
/* > is singular to working precision. This condition is */
/* > indicated by a return code of INFO > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* > FERR is DOUBLE PRECISION array, dimension (NRHS) */
/* > The estimated forward error bound for each solution vector */
/* > X(j) (the j-th column of the solution matrix X). */
/* > If XTRUE is the true solution corresponding to X(j), FERR(j) */
/* > is an estimated upper bound for the magnitude of the largest */
/* > element in (X(j) - XTRUE) divided by the magnitude of the */
/* > largest element in X(j). The estimate is as reliable as */
/* > the estimate for RCOND, and is almost always a slight */
/* > overestimate of the true error. */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* > BERR is DOUBLE PRECISION array, dimension (NRHS) */
/* > The componentwise relative backward error of each solution */
/* > vector X(j) (i.e., the smallest relative change in */
/* > any element of A or B that makes X(j) an exact solution). */
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
/* > > 0: if INFO = i, and i is */
/* > <= N: the leading minor of order i of A is */
/* > not positive definite, so the factorization */
/* > could not be completed, and the solution has not */
/* > been computed. RCOND = 0 is returned. */
/* > = N+1: U is nonsingular, but RCOND is less than machine */
/* > precision, meaning that the matrix is singular */
/* > to working precision. Nevertheless, the */
/* > solution and error bounds are computed because */
/* > there are a number of situations where the */
/* > computed solution can be more accurate than the */
/* > value of RCOND would suggest. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date April 2012 */
/* > \ingroup complex16OTHERsolve */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The packed storage scheme is illustrated by the following example */
/* > when N = 4, UPLO = 'U': */
/* > */
/* > Two-dimensional storage of the Hermitian matrix A: */
/* > */
/* > a11 a12 a13 a14 */
/* > a22 a23 a24 */
/* > a33 a34 (aij = conjg(aji)) */
/* > a44 */
/* > */
/* > Packed storage of the upper triangle of A: */
/* > */
/* > AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ] */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int zppsvx_(char *fact, char *uplo, integer *n, integer * nrhs, doublecomplex *ap, doublecomplex *afp, char *equed, doublereal * s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex * work, doublereal *rwork, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1;
    /* Local variables */
    integer i__, j;
    doublereal amax, smin, smax;
    extern logical lsame_(char *, char *);
    doublereal scond, anorm;
    logical equil, rcequ;
    extern /* Subroutine */
    int zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *);
    logical nofact;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal bignum;
    integer infequ;
    extern doublereal zlanhp_(char *, char *, integer *, doublecomplex *, doublereal *);
    extern /* Subroutine */
    int zlaqhp_(char *, integer *, doublecomplex *, doublereal *, doublereal *, doublereal *, char *), zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *), zppcon_(char *, integer *, doublecomplex *, doublereal *, doublereal *, doublecomplex *, doublereal *, integer *);
    doublereal smlnum;
    extern /* Subroutine */
    int zppequ_(char *, integer *, doublecomplex *, doublereal *, doublereal *, doublereal *, integer *), zpprfs_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, doublereal *, doublecomplex *, doublereal *, integer *), zpptrf_(char *, integer *, doublecomplex *, integer *), zpptrs_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, integer *);
    /* -- LAPACK driver routine (version 3.4.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* April 2012 */
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
    /* Parameter adjustments */
    --ap;
    --afp;
    --s;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;
    --rwork;
    /* Function Body */
    *info = 0;
    nofact = lsame_(fact, "N");
    equil = lsame_(fact, "E");
    if (nofact || equil)
    {
        *(unsigned char *)equed = 'N';
        rcequ = FALSE_;
    }
    else
    {
        rcequ = lsame_(equed, "Y");
        smlnum = dlamch_("Safe minimum");
        bignum = 1. / smlnum;
    }
    /* Test the input parameters. */
    if (! nofact && ! equil && ! lsame_(fact, "F"))
    {
        *info = -1;
    }
    else if (! lsame_(uplo, "U") && ! lsame_(uplo, "L"))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*nrhs < 0)
    {
        *info = -4;
    }
    else if (lsame_(fact, "F") && ! (rcequ || lsame_( equed, "N")))
    {
        *info = -7;
    }
    else
    {
        if (rcequ)
        {
            smin = bignum;
            smax = 0.;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Computing MIN */
                d__1 = smin;
                d__2 = s[j]; // , expr subst
                smin = min(d__1,d__2);
                /* Computing MAX */
                d__1 = smax;
                d__2 = s[j]; // , expr subst
                smax = max(d__1,d__2);
                /* L10: */
            }
            if (smin <= 0.)
            {
                *info = -8;
            }
            else if (*n > 0)
            {
                scond = max(smin,smlnum) / min(smax,bignum);
            }
            else
            {
                scond = 1.;
            }
        }
        if (*info == 0)
        {
            if (*ldb < max(1,*n))
            {
                *info = -10;
            }
            else if (*ldx < max(1,*n))
            {
                *info = -12;
            }
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZPPSVX", &i__1);
        return 0;
    }
    if (equil)
    {
        /* Compute row and column scalings to equilibrate the matrix A. */
        zppequ_(uplo, n, &ap[1], &s[1], &scond, &amax, &infequ);
        if (infequ == 0)
        {
            /* Equilibrate the matrix. */
            zlaqhp_(uplo, n, &ap[1], &s[1], &scond, &amax, equed);
            rcequ = lsame_(equed, "Y");
        }
    }
    /* Scale the right-hand side. */
    if (rcequ)
    {
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *n;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * b_dim1;
                i__4 = i__;
                i__5 = i__ + j * b_dim1;
                z__1.r = s[i__4] * b[i__5].r;
                z__1.i = s[i__4] * b[i__5].i; // , expr subst
                b[i__3].r = z__1.r;
                b[i__3].i = z__1.i; // , expr subst
                /* L20: */
            }
            /* L30: */
        }
    }
    if (nofact || equil)
    {
        /* Compute the Cholesky factorization A = U**H * U or A = L * L**H. */
        i__1 = *n * (*n + 1) / 2;
        zcopy_(&i__1, &ap[1], &c__1, &afp[1], &c__1);
        zpptrf_(uplo, n, &afp[1], info);
        /* Return if INFO is non-zero. */
        if (*info > 0)
        {
            *rcond = 0.;
            return 0;
        }
    }
    /* Compute the norm of the matrix A. */
    anorm = zlanhp_("I", uplo, n, &ap[1], &rwork[1]);
    /* Compute the reciprocal of the condition number of A. */
    zppcon_(uplo, n, &afp[1], &anorm, rcond, &work[1], &rwork[1], info);
    /* Compute the solution matrix X. */
    zlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    zpptrs_(uplo, n, nrhs, &afp[1], &x[x_offset], ldx, info);
    /* Use iterative refinement to improve the computed solution and */
    /* compute error bounds and backward error estimates for it. */
    zpprfs_(uplo, n, nrhs, &ap[1], &afp[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1], &rwork[1], info);
    /* Transform the solution matrix X to a solution of the original */
    /* system. */
    if (rcequ)
    {
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *n;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * x_dim1;
                i__4 = i__;
                i__5 = i__ + j * x_dim1;
                z__1.r = s[i__4] * x[i__5].r;
                z__1.i = s[i__4] * x[i__5].i; // , expr subst
                x[i__3].r = z__1.r;
                x[i__3].i = z__1.i; // , expr subst
                /* L40: */
            }
            /* L50: */
        }
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            ferr[j] /= scond;
            /* L60: */
        }
    }
    /* Set INFO = N+1 if the matrix is singular to working precision. */
    if (*rcond < dlamch_("Epsilon"))
    {
        *info = *n + 1;
    }
    return 0;
    /* End of ZPPSVX */
}
/* zppsvx_ */
