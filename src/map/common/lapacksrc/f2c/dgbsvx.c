/* ../netlib/dgbsvx.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief <b> DGBSVX computes the solution to system of linear equations A * X = B for GB matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGBSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbsvx. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbsvx. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbsvx. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, */
/* LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, */
/* RCOND, FERR, BERR, WORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED, FACT, TRANS */
/* INTEGER INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS */
/* DOUBLE PRECISION RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), IWORK( * ) */
/* DOUBLE PRECISION AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/* $ BERR( * ), C( * ), FERR( * ), R( * ), */
/* $ WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGBSVX uses the LU factorization to compute the solution to a real */
/* > system of linear equations A * X = B, A**T * X = B, or A**H * X = B, */
/* > where A is a band matrix of order N with KL subdiagonals and KU */
/* > superdiagonals, and X and B are N-by-NRHS matrices. */
/* > */
/* > Error bounds on the solution and a condition estimate are also */
/* > provided. */
/* > \endverbatim */
/* > \par Description: */
/* ================= */
/* > */
/* > \verbatim */
/* > */
/* > The following steps are performed by this subroutine: */
/* > */
/* > 1. If FACT = 'E', real scaling factors are computed to equilibrate */
/* > the system: */
/* > TRANS = 'N': diag(R)*A*diag(C) *inv(diag(C))*X = diag(R)*B */
/* > TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B */
/* > TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B */
/* > Whether or not the system will be equilibrated depends on the */
/* > scaling of the matrix A, but if equilibration is used, A is */
/* > overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N') */
/* > or diag(C)*B (if TRANS = 'T' or 'C'). */
/* > */
/* > 2. If FACT = 'N' or 'E', the LU decomposition is used to factor the */
/* > matrix A (after equilibration if FACT = 'E') as */
/* > A = L * U, */
/* > where L is a product of permutation and unit lower triangular */
/* > matrices with KL subdiagonals, and U is upper triangular with */
/* > KL+KU superdiagonals. */
/* > */
/* > 3. If some U(i,i)=0, so that U is exactly singular, then the routine */
/* > returns with INFO = i. Otherwise, the factored form of A is used */
/* > to estimate the condition number of the matrix A. If the */
/* > reciprocal of the condition number is less than machine precision, */
/* > INFO = N+1 is returned as a warning, but the routine still goes on */
/* > to solve for X and compute error bounds as described below. */
/* > */
/* > 4. The system of equations is solved for X using the factored form */
/* > of A. */
/* > */
/* > 5. Iterative refinement is applied to improve the computed solution */
/* > matrix and calculate error bounds and backward error estimates */
/* > for it. */
/* > */
/* > 6. If equilibration was used, the matrix X is premultiplied by */
/* > diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so */
/* > that it solves the original system before equilibration. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] FACT */
/* > \verbatim */
/* > FACT is CHARACTER*1 */
/* > Specifies whether or not the factored form of the matrix A is */
/* > supplied on entry, and if not, whether the matrix A should be */
/* > equilibrated before it is factored. */
/* > = 'F': On entry, AFB and IPIV contain the factored form of */
/* > A. If EQUED is not 'N', the matrix A has been */
/* > equilibrated with scaling factors given by R and C. */
/* > AB, AFB, and IPIV are not modified. */
/* > = 'N': The matrix A will be copied to AFB and factored. */
/* > = 'E': The matrix A will be equilibrated if necessary, then */
/* > copied to AFB and factored. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations. */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* > KL is INTEGER */
/* > The number of subdiagonals within the band of A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* > KU is INTEGER */
/* > The number of superdiagonals within the band of A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* > On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* > The j-th column of A is stored in the j-th column of the */
/* > array AB as follows: */
/* > AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl) */
/* > */
/* > If FACT = 'F' and EQUED is not 'N', then A must have been */
/* > equilibrated by the scaling factors in R and/or C. AB is not */
/* > modified if FACT = 'F' or 'N', or if FACT = 'E' and */
/* > EQUED = 'N' on exit. */
/* > */
/* > On exit, if EQUED .ne. 'N', A is scaled as follows: */
/* > EQUED = 'R': A := diag(R) * A */
/* > EQUED = 'C': A := A * diag(C) */
/* > EQUED = 'B': A := diag(R) * A * diag(C). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AFB */
/* > \verbatim */
/* > AFB is DOUBLE PRECISION array, dimension (LDAFB,N) */
/* > If FACT = 'F', then AFB is an input argument and on entry */
/* > contains details of the LU factorization of the band matrix */
/* > A, as computed by DGBTRF. U is stored as an upper triangular */
/* > band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, */
/* > and the multipliers used during the factorization are stored */
/* > in rows KL+KU+2 to 2*KL+KU+1. If EQUED .ne. 'N', then AFB is */
/* > the factored form of the equilibrated matrix A. */
/* > */
/* > If FACT = 'N', then AFB is an output argument and on exit */
/* > returns details of the LU factorization of A. */
/* > */
/* > If FACT = 'E', then AFB is an output argument and on exit */
/* > returns details of the LU factorization of the equilibrated */
/* > matrix A (see the description of AB for the form of the */
/* > equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* > LDAFB is INTEGER */
/* > The leading dimension of the array AFB. LDAFB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > If FACT = 'F', then IPIV is an input argument and on entry */
/* > contains the pivot indices from the factorization A = L*U */
/* > as computed by DGBTRF;
row i of the matrix was interchanged */
/* > with row IPIV(i). */
/* > */
/* > If FACT = 'N', then IPIV is an output argument and on exit */
/* > contains the pivot indices from the factorization A = L*U */
/* > of the original matrix A. */
/* > */
/* > If FACT = 'E', then IPIV is an output argument and on exit */
/* > contains the pivot indices from the factorization A = L*U */
/* > of the equilibrated matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EQUED */
/* > \verbatim */
/* > EQUED is CHARACTER*1 */
/* > Specifies the form of equilibration that was done. */
/* > = 'N': No equilibration (always true if FACT = 'N'). */
/* > = 'R': Row equilibration, i.e., A has been premultiplied by */
/* > diag(R). */
/* > = 'C': Column equilibration, i.e., A has been postmultiplied */
/* > by diag(C). */
/* > = 'B': Both row and column equilibration, i.e., A has been */
/* > replaced by diag(R) * A * diag(C). */
/* > EQUED is an input argument if FACT = 'F';
otherwise, it is an */
/* > output argument. */
/* > \endverbatim */
/* > */
/* > \param[in,out] R */
/* > \verbatim */
/* > R is DOUBLE PRECISION array, dimension (N) */
/* > The row scale factors for A. If EQUED = 'R' or 'B', A is */
/* > multiplied on the left by diag(R);
if EQUED = 'N' or 'C', R */
/* > is not accessed. R is an input argument if FACT = 'F';
*/
/* > otherwise, R is an output argument. If FACT = 'F' and */
/* > EQUED = 'R' or 'B', each element of R must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (N) */
/* > The column scale factors for A. If EQUED = 'C' or 'B', A is */
/* > multiplied on the right by diag(C);
if EQUED = 'N' or 'R', C */
/* > is not accessed. C is an input argument if FACT = 'F';
*/
/* > otherwise, C is an output argument. If FACT = 'F' and */
/* > EQUED = 'C' or 'B', each element of C must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* > On entry, the right hand side matrix B. */
/* > On exit, */
/* > if EQUED = 'N', B is not modified;
*/
/* > if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by */
/* > diag(R)*B;
*/
/* > if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is */
/* > overwritten by diag(C)*B. */
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
/* > X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
/* > If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X */
/* > to the original system of equations. Note that A and B are */
/* > modified on exit if EQUED .ne. 'N', and the solution to the */
/* > equilibrated system is inv(diag(C))*X if TRANS = 'N' and */
/* > EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C' */
/* > and EQUED = 'R' or 'B'. */
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
/* > WORK is DOUBLE PRECISION array, dimension (3*N) */
/* > On exit, WORK(1) contains the reciprocal pivot growth */
/* > factor norm(A)/norm(U). The "max absolute element" norm is */
/* > used. If WORK(1) is much less than 1, then the stability */
/* > of the LU factorization of the (equilibrated) matrix A */
/* > could be poor. This also means that the solution X, condition */
/* > estimator RCOND, and forward error bound FERR could be */
/* > unreliable. If factorization fails with 0<INFO<=N, then */
/* > WORK(1) contains the reciprocal pivot growth factor for the */
/* > leading INFO columns of A. */
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
/* > > 0: if INFO = i, and i is */
/* > <= N: U(i,i) is exactly zero. The factorization */
/* > has been completed, but the factor U is exactly */
/* > singular, so the solution and error bounds */
/* > could not be computed. RCOND = 0 is returned. */
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
/* > \ingroup doubleGBsolve */
/* ===================================================================== */
/* Subroutine */
int dgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, integer *ipiv, char *equed, doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;
    /* Local variables */
    integer i__, j, j1, j2;
    doublereal amax;
    char norm[1];
    extern logical lsame_(char *, char *);
    doublereal rcmin, rcmax, anorm;
    extern /* Subroutine */
    int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical equil;
    extern doublereal dlangb_(char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *), dlamch_(char *);
    extern /* Subroutine */
    int dlaqgb_(integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, char *), dgbcon_(char *, integer *, integer *, integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, integer *);
    doublereal colcnd;
    extern doublereal dlantb_(char *, char *, char *, integer *, integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */
    int dgbequ_(integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, integer *), dgbrfs_( char *, integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, integer *, integer *), dgbtrf_(integer *, integer *, integer *, integer *, doublereal *, integer *, integer *, integer *);
    logical nofact;
    extern /* Subroutine */
    int dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *), xerbla_(char *, integer *);
    doublereal bignum;
    extern /* Subroutine */
    int dgbtrs_(char *, integer *, integer *, integer *, integer *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
    integer infequ;
    logical colequ;
    doublereal rowcnd;
    logical notran;
    doublereal smlnum;
    logical rowequ;
    doublereal rpvgrw;
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
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    afb_dim1 = *ldafb;
    afb_offset = 1 + afb_dim1;
    afb -= afb_offset;
    --ipiv;
    --r__;
    --c__;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    nofact = lsame_(fact, "N");
    equil = lsame_(fact, "E");
    notran = lsame_(trans, "N");
    if (nofact || equil)
    {
        *(unsigned char *)equed = 'N';
        rowequ = FALSE_;
        colequ = FALSE_;
    }
    else
    {
        rowequ = lsame_(equed, "R") || lsame_(equed, "B");
        colequ = lsame_(equed, "C") || lsame_(equed, "B");
        smlnum = dlamch_("Safe minimum");
        bignum = 1. / smlnum;
    }
    /* Test the input parameters. */
    if (! nofact && ! equil && ! lsame_(fact, "F"))
    {
        *info = -1;
    }
    else if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, "C"))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*kl < 0)
    {
        *info = -4;
    }
    else if (*ku < 0)
    {
        *info = -5;
    }
    else if (*nrhs < 0)
    {
        *info = -6;
    }
    else if (*ldab < *kl + *ku + 1)
    {
        *info = -8;
    }
    else if (*ldafb < (*kl << 1) + *ku + 1)
    {
        *info = -10;
    }
    else if (lsame_(fact, "F") && ! (rowequ || colequ || lsame_(equed, "N")))
    {
        *info = -12;
    }
    else
    {
        if (rowequ)
        {
            rcmin = bignum;
            rcmax = 0.;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Computing MIN */
                d__1 = rcmin;
                d__2 = r__[j]; // , expr subst
                rcmin = min(d__1,d__2);
                /* Computing MAX */
                d__1 = rcmax;
                d__2 = r__[j]; // , expr subst
                rcmax = max(d__1,d__2);
                /* L10: */
            }
            if (rcmin <= 0.)
            {
                *info = -13;
            }
            else if (*n > 0)
            {
                rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
            }
            else
            {
                rowcnd = 1.;
            }
        }
        if (colequ && *info == 0)
        {
            rcmin = bignum;
            rcmax = 0.;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Computing MIN */
                d__1 = rcmin;
                d__2 = c__[j]; // , expr subst
                rcmin = min(d__1,d__2);
                /* Computing MAX */
                d__1 = rcmax;
                d__2 = c__[j]; // , expr subst
                rcmax = max(d__1,d__2);
                /* L20: */
            }
            if (rcmin <= 0.)
            {
                *info = -14;
            }
            else if (*n > 0)
            {
                colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
            }
            else
            {
                colcnd = 1.;
            }
        }
        if (*info == 0)
        {
            if (*ldb < max(1,*n))
            {
                *info = -16;
            }
            else if (*ldx < max(1,*n))
            {
                *info = -18;
            }
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGBSVX", &i__1);
        return 0;
    }
    if (equil)
    {
        /* Compute row and column scalings to equilibrate the matrix A. */
        dgbequ_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &rowcnd, &colcnd, &amax, &infequ);
        if (infequ == 0)
        {
            /* Equilibrate the matrix. */
            dlaqgb_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], & rowcnd, &colcnd, &amax, equed);
            rowequ = lsame_(equed, "R") || lsame_(equed, "B");
            colequ = lsame_(equed, "C") || lsame_(equed, "B");
        }
    }
    /* Scale the right hand side. */
    if (notran)
    {
        if (rowequ)
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
                    b[i__ + j * b_dim1] = r__[i__] * b[i__ + j * b_dim1];
                    /* L30: */
                }
                /* L40: */
            }
        }
    }
    else if (colequ)
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
                b[i__ + j * b_dim1] = c__[i__] * b[i__ + j * b_dim1];
                /* L50: */
            }
            /* L60: */
        }
    }
    if (nofact || equil)
    {
        /* Compute the LU factorization of the band matrix A. */
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Computing MAX */
            i__2 = j - *ku;
            j1 = max(i__2,1);
            /* Computing MIN */
            i__2 = j + *kl;
            j2 = min(i__2,*n);
            i__2 = j2 - j1 + 1;
            dcopy_(&i__2, &ab[*ku + 1 - j + j1 + j * ab_dim1], &c__1, &afb[* kl + *ku + 1 - j + j1 + j * afb_dim1], &c__1);
            /* L70: */
        }
        dgbtrf_(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);
        /* Return if INFO is non-zero. */
        if (*info > 0)
        {
            /* Compute the reciprocal pivot growth factor of the */
            /* leading rank-deficient INFO columns of A. */
            anorm = 0.;
            i__1 = *info;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Computing MAX */
                i__2 = *ku + 2 - j;
                /* Computing MIN */
                i__4 = *n + *ku + 1 - j;
                i__5 = *kl + *ku + 1; // , expr subst
                i__3 = min(i__4,i__5);
                for (i__ = max(i__2,1);
                        i__ <= i__3;
                        ++i__)
                {
                    /* Computing MAX */
                    d__2 = anorm;
                    d__3 = (d__1 = ab[i__ + j * ab_dim1], f2c_abs( d__1)); // , expr subst
                    anorm = max(d__2,d__3);
                    /* L80: */
                }
                /* L90: */
            }
            /* Computing MIN */
            i__3 = *info - 1;
            i__2 = *kl + *ku; // , expr subst
            i__1 = min(i__3,i__2);
            /* Computing MAX */
            i__4 = 1;
            i__5 = *kl + *ku + 2 - *info; // , expr subst
            rpvgrw = dlantb_("M", "U", "N", info, &i__1, &afb[max(i__4,i__5) + afb_dim1], ldafb, &work[1]);
            if (rpvgrw == 0.)
            {
                rpvgrw = 1.;
            }
            else
            {
                rpvgrw = anorm / rpvgrw;
            }
            work[1] = rpvgrw;
            *rcond = 0.;
            return 0;
        }
    }
    /* Compute the norm of the matrix A and the */
    /* reciprocal pivot growth factor RPVGRW. */
    if (notran)
    {
        *(unsigned char *)norm = '1';
    }
    else
    {
        *(unsigned char *)norm = 'I';
    }
    anorm = dlangb_(norm, n, kl, ku, &ab[ab_offset], ldab, &work[1]);
    i__1 = *kl + *ku;
    rpvgrw = dlantb_("M", "U", "N", n, &i__1, &afb[afb_offset], ldafb, &work[ 1]);
    if (rpvgrw == 0.)
    {
        rpvgrw = 1.;
    }
    else
    {
        rpvgrw = dlangb_("M", n, kl, ku, &ab[ab_offset], ldab, &work[1]) / rpvgrw;
    }
    /* Compute the reciprocal of the condition number of A. */
    dgbcon_(norm, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], &anorm, rcond, &work[1], &iwork[1], info);
    /* Compute the solution matrix X. */
    dlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    dgbtrs_(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[ x_offset], ldx, info);
    /* Use iterative refinement to improve the computed solution and */
    /* compute error bounds and backward error estimates for it. */
    dgbrfs_(trans, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[afb_offset], ldafb, &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], & berr[1], &work[1], &iwork[1], info);
    /* Transform the solution matrix X to a solution of the original */
    /* system. */
    if (notran)
    {
        if (colequ)
        {
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__3 = *n;
                for (i__ = 1;
                        i__ <= i__3;
                        ++i__)
                {
                    x[i__ + j * x_dim1] = c__[i__] * x[i__ + j * x_dim1];
                    /* L100: */
                }
                /* L110: */
            }
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                ferr[j] /= colcnd;
                /* L120: */
            }
        }
    }
    else if (rowequ)
    {
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__3 = *n;
            for (i__ = 1;
                    i__ <= i__3;
                    ++i__)
            {
                x[i__ + j * x_dim1] = r__[i__] * x[i__ + j * x_dim1];
                /* L130: */
            }
            /* L140: */
        }
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            ferr[j] /= rowcnd;
            /* L150: */
        }
    }
    /* Set INFO = N+1 if the matrix is singular to working precision. */
    if (*rcond < dlamch_("Epsilon"))
    {
        *info = *n + 1;
    }
    work[1] = rpvgrw;
    return 0;
    /* End of DGBSVX */
}
/* dgbsvx_ */
