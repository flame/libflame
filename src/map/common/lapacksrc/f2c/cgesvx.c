/* ../netlib/cgesvx.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief <b> CGESVX computes the solution to system of linear equations A * X = B for GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGESVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesvx. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesvx. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesvx. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, */
/* EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, */
/* WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED, FACT, TRANS */
/* INTEGER INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/* REAL RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL BERR( * ), C( * ), FERR( * ), R( * ), */
/* $ RWORK( * ) */
/* COMPLEX A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/* $ WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGESVX uses the LU factorization to compute the solution to a complex */
/* > system of linear equations */
/* > A * X = B, */
/* > where A is an N-by-N matrix and X and B are N-by-NRHS matrices. */
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
/* > A = P * L * U, */
/* > where P is a permutation matrix, L is a unit lower triangular */
/* > matrix, and U is upper triangular. */
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
/* > = 'F': On entry, AF and IPIV contain the factored form of A. */
/* > If EQUED is not 'N', the matrix A has been */
/* > equilibrated with scaling factors given by R and C. */
/* > A, AF, and IPIV are not modified. */
/* > = 'N': The matrix A will be copied to AF and factored. */
/* > = 'E': The matrix A will be equilibrated if necessary, then */
/* > copied to AF and factored. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations: */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate transpose) */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the N-by-N matrix A. If FACT = 'F' and EQUED is */
/* > not 'N', then A must have been equilibrated by the scaling */
/* > factors in R and/or C. A is not modified if FACT = 'F' or */
/* > 'N', or if FACT = 'E' and EQUED = 'N' on exit. */
/* > */
/* > On exit, if EQUED .ne. 'N', A is scaled as follows: */
/* > EQUED = 'R': A := diag(R) * A */
/* > EQUED = 'C': A := A * diag(C) */
/* > EQUED = 'B': A := diag(R) * A * diag(C). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] AF */
/* > \verbatim */
/* > AF is COMPLEX array, dimension (LDAF,N) */
/* > If FACT = 'F', then AF is an input argument and on entry */
/* > contains the factors L and U from the factorization */
/* > A = P*L*U as computed by CGETRF. If EQUED .ne. 'N', then */
/* > AF is the factored form of the equilibrated matrix A. */
/* > */
/* > If FACT = 'N', then AF is an output argument and on exit */
/* > returns the factors L and U from the factorization A = P*L*U */
/* > of the original matrix A. */
/* > */
/* > If FACT = 'E', then AF is an output argument and on exit */
/* > returns the factors L and U from the factorization A = P*L*U */
/* > of the equilibrated matrix A (see the description of A for */
/* > the form of the equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > If FACT = 'F', then IPIV is an input argument and on entry */
/* > contains the pivot indices from the factorization A = P*L*U */
/* > as computed by CGETRF;
row i of the matrix was interchanged */
/* > with row IPIV(i). */
/* > */
/* > If FACT = 'N', then IPIV is an output argument and on exit */
/* > contains the pivot indices from the factorization A = P*L*U */
/* > of the original matrix A. */
/* > */
/* > If FACT = 'E', then IPIV is an output argument and on exit */
/* > contains the pivot indices from the factorization A = P*L*U */
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
/* > R is REAL array, dimension (N) */
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
/* > C is REAL array, dimension (N) */
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
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the N-by-NRHS right hand side matrix B. */
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
/* > X is COMPLEX array, dimension (LDX,NRHS) */
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
/* > RCOND is REAL */
/* > The estimate of the reciprocal condition number of the matrix */
/* > A after equilibration (if done). If RCOND is less than the */
/* > machine precision (in particular, if RCOND = 0), the matrix */
/* > is singular to working precision. This condition is */
/* > indicated by a return code of INFO > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* > FERR is REAL array, dimension (NRHS) */
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
/* > BERR is REAL array, dimension (NRHS) */
/* > The componentwise relative backward error of each solution */
/* > vector X(j) (i.e., the smallest relative change in */
/* > any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (2*N) */
/* > On exit, RWORK(1) contains the reciprocal pivot growth */
/* > factor norm(A)/norm(U). The "max absolute element" norm is */
/* > used. If RWORK(1) is much less than 1, then the stability */
/* > of the LU factorization of the (equilibrated) matrix A */
/* > could be poor. This also means that the solution X, condition */
/* > estimator RCOND, and forward error bound FERR could be */
/* > unreliable. If factorization fails with 0<INFO<=N, then */
/* > RWORK(1) contains the reciprocal pivot growth factor for the */
/* > leading INFO columns of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, and i is */
/* > <= N: U(i,i) is exactly zero. The factorization has */
/* > been completed, but the factor U is exactly */
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
/* > \ingroup complexGEsolve */
/* ===================================================================== */
/* Subroutine */
int cgesvx_(char *fact, char *trans, integer *n, integer * nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer * ipiv, char *equed, real *r__, real *c__, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, real *rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2;
    complex q__1;
    /* Local variables */
    integer i__, j;
    real amax;
    char norm[1];
    extern logical lsame_(char *, char *);
    real rcmin, rcmax, anorm;
    logical equil;
    extern real clange_(char *, integer *, integer *, complex *, integer *, real *);
    extern /* Subroutine */
    int claqge_(integer *, integer *, complex *, integer *, real *, real *, real *, real *, real *, char *) , cgecon_(char *, integer *, complex *, integer *, real *, real *, complex *, real *, integer *);
    real colcnd;
    extern real slamch_(char *);
    extern /* Subroutine */
    int cgeequ_(integer *, integer *, complex *, integer *, real *, real *, real *, real *, real *, integer *);
    logical nofact;
    extern /* Subroutine */
    int cgerfs_(char *, integer *, integer *, complex *, integer *, complex *, integer *, integer *, complex *, integer *, complex *, integer *, real *, real *, complex *, real *, integer *), cgetrf_(integer *, integer *, complex *, integer *, integer *, integer *), clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), xerbla_(char *, integer *);
    real bignum;
    extern real clantr_(char *, char *, char *, integer *, integer *, complex *, integer *, real *);
    integer infequ;
    logical colequ;
    extern /* Subroutine */
    int cgetrs_(char *, integer *, integer *, complex *, integer *, integer *, complex *, integer *, integer *);
    real rowcnd;
    logical notran;
    real smlnum;
    logical rowequ;
    real rpvgrw;
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
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
    --rwork;
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
        smlnum = slamch_("Safe minimum");
        bignum = 1.f / smlnum;
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
    else if (*nrhs < 0)
    {
        *info = -4;
    }
    else if (*lda < max(1,*n))
    {
        *info = -6;
    }
    else if (*ldaf < max(1,*n))
    {
        *info = -8;
    }
    else if (lsame_(fact, "F") && ! (rowequ || colequ || lsame_(equed, "N")))
    {
        *info = -10;
    }
    else
    {
        if (rowequ)
        {
            rcmin = bignum;
            rcmax = 0.f;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Computing MIN */
                r__1 = rcmin;
                r__2 = r__[j]; // , expr subst
                rcmin = min(r__1,r__2);
                /* Computing MAX */
                r__1 = rcmax;
                r__2 = r__[j]; // , expr subst
                rcmax = max(r__1,r__2);
                /* L10: */
            }
            if (rcmin <= 0.f)
            {
                *info = -11;
            }
            else if (*n > 0)
            {
                rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
            }
            else
            {
                rowcnd = 1.f;
            }
        }
        if (colequ && *info == 0)
        {
            rcmin = bignum;
            rcmax = 0.f;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Computing MIN */
                r__1 = rcmin;
                r__2 = c__[j]; // , expr subst
                rcmin = min(r__1,r__2);
                /* Computing MAX */
                r__1 = rcmax;
                r__2 = c__[j]; // , expr subst
                rcmax = max(r__1,r__2);
                /* L20: */
            }
            if (rcmin <= 0.f)
            {
                *info = -12;
            }
            else if (*n > 0)
            {
                colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
            }
            else
            {
                colcnd = 1.f;
            }
        }
        if (*info == 0)
        {
            if (*ldb < max(1,*n))
            {
                *info = -14;
            }
            else if (*ldx < max(1,*n))
            {
                *info = -16;
            }
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGESVX", &i__1);
        return 0;
    }
    if (equil)
    {
        /* Compute row and column scalings to equilibrate the matrix A. */
        cgeequ_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &colcnd, & amax, &infequ);
        if (infequ == 0)
        {
            /* Equilibrate the matrix. */
            claqge_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, & colcnd, &amax, equed);
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
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__;
                    i__5 = i__ + j * b_dim1;
                    q__1.r = r__[i__4] * b[i__5].r;
                    q__1.i = r__[i__4] * b[ i__5].i; // , expr subst
                    b[i__3].r = q__1.r;
                    b[i__3].i = q__1.i; // , expr subst
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
                i__3 = i__ + j * b_dim1;
                i__4 = i__;
                i__5 = i__ + j * b_dim1;
                q__1.r = c__[i__4] * b[i__5].r;
                q__1.i = c__[i__4] * b[i__5] .i; // , expr subst
                b[i__3].r = q__1.r;
                b[i__3].i = q__1.i; // , expr subst
                /* L50: */
            }
            /* L60: */
        }
    }
    if (nofact || equil)
    {
        /* Compute the LU factorization of A. */
        clacpy_("Full", n, n, &a[a_offset], lda, &af[af_offset], ldaf);
        cgetrf_(n, n, &af[af_offset], ldaf, &ipiv[1], info);
        /* Return if INFO is non-zero. */
        if (*info > 0)
        {
            /* Compute the reciprocal pivot growth factor of the */
            /* leading rank-deficient INFO columns of A. */
            rpvgrw = clantr_("M", "U", "N", info, info, &af[af_offset], ldaf, &rwork[1]);
            if (rpvgrw == 0.f)
            {
                rpvgrw = 1.f;
            }
            else
            {
                rpvgrw = clange_("M", n, info, &a[a_offset], lda, &rwork[1]) / rpvgrw;
            }
            rwork[1] = rpvgrw;
            *rcond = 0.f;
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
    anorm = clange_(norm, n, n, &a[a_offset], lda, &rwork[1]);
    rpvgrw = clantr_("M", "U", "N", n, n, &af[af_offset], ldaf, &rwork[1]);
    if (rpvgrw == 0.f)
    {
        rpvgrw = 1.f;
    }
    else
    {
        rpvgrw = clange_("M", n, n, &a[a_offset], lda, &rwork[1]) / rpvgrw;
    }
    /* Compute the reciprocal of the condition number of A. */
    cgecon_(norm, n, &af[af_offset], ldaf, &anorm, rcond, &work[1], &rwork[1], info);
    /* Compute the solution matrix X. */
    clacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    cgetrs_(trans, n, nrhs, &af[af_offset], ldaf, &ipiv[1], &x[x_offset], ldx, info);
    /* Use iterative refinement to improve the computed solution and */
    /* compute error bounds and backward error estimates for it. */
    cgerfs_(trans, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[ 1], &rwork[1], info);
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
                i__2 = *n;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * x_dim1;
                    i__4 = i__;
                    i__5 = i__ + j * x_dim1;
                    q__1.r = c__[i__4] * x[i__5].r;
                    q__1.i = c__[i__4] * x[ i__5].i; // , expr subst
                    x[i__3].r = q__1.r;
                    x[i__3].i = q__1.i; // , expr subst
                    /* L70: */
                }
                /* L80: */
            }
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                ferr[j] /= colcnd;
                /* L90: */
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
            i__2 = *n;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * x_dim1;
                i__4 = i__;
                i__5 = i__ + j * x_dim1;
                q__1.r = r__[i__4] * x[i__5].r;
                q__1.i = r__[i__4] * x[i__5] .i; // , expr subst
                x[i__3].r = q__1.r;
                x[i__3].i = q__1.i; // , expr subst
                /* L100: */
            }
            /* L110: */
        }
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            ferr[j] /= rowcnd;
            /* L120: */
        }
    }
    /* Set INFO = N+1 if the matrix is singular to working precision. */
    if (*rcond < slamch_("Epsilon"))
    {
        *info = *n + 1;
    }
    rwork[1] = rpvgrw;
    return 0;
    /* End of CGESVX */
}
/* cgesvx_ */
