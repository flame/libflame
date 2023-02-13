/* ../netlib/cbbcsd.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    -1.f,0.f
}
;
static doublereal c_b11 = -.125;
static integer c__1 = 1;
/* > \brief \b CBBCSD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CBBCSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cbbcsd. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cbbcsd. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cbbcsd. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, */
/* THETA, PHI, U1, LDU1, U2, LDU2, V1T, LDV1T, */
/* V2T, LDV2T, B11D, B11E, B12D, B12E, B21D, B21E, */
/* B22D, B22E, RWORK, LRWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS */
/* INTEGER INFO, LDU1, LDU2, LDV1T, LDV2T, LRWORK, M, P, Q */
/* .. */
/* .. Array Arguments .. */
/* REAL B11D( * ), B11E( * ), B12D( * ), B12E( * ), */
/* $ B21D( * ), B21E( * ), B22D( * ), B22E( * ), */
/* $ PHI( * ), THETA( * ), RWORK( * ) */
/* COMPLEX U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), */
/* $ V2T( LDV2T, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CBBCSD computes the CS decomposition of a unitary matrix in */
/* > bidiagonal-block form, */
/* > */
/* > */
/* > [ B11 | B12 0 0 ] */
/* > [ 0 | 0 -I 0 ] */
/* > X = [----------------] */
/* > [ B21 | B22 0 0 ] */
/* > [ 0 | 0 0 I ] */
/* > */
/* > [ C | -S 0 0 ] */
/* > [ U1 | ] [ 0 | 0 -I 0 ] [ V1 | ]**H */
/* > = [---------] [---------------] [---------] . */
/* > [ | U2 ] [ S | C 0 0 ] [ | V2 ] */
/* > [ 0 | 0 0 I ] */
/* > */
/* > X is M-by-M, its top-left block is P-by-Q, and Q must be no larger */
/* > than P, M-P, or M-Q. (If Q is not the smallest index, then X must be */
/* > transposed and/or permuted. This can be done in constant time using */
/* > the TRANS and SIGNS options. See CUNCSD for details.) */
/* > */
/* > The bidiagonal matrices B11, B12, B21, and B22 are represented */
/* > implicitly by angles THETA(1:Q) and PHI(1:Q-1). */
/* > */
/* > The unitary matrices U1, U2, V1T, and V2T are input/output. */
/* > The input matrices are pre- or post-multiplied by the appropriate */
/* > singular vector matrices. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBU1 */
/* > \verbatim */
/* > JOBU1 is CHARACTER */
/* > = 'Y': U1 is updated;
*/
/* > otherwise: U1 is not updated. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU2 */
/* > \verbatim */
/* > JOBU2 is CHARACTER */
/* > = 'Y': U2 is updated;
*/
/* > otherwise: U2 is not updated. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV1T */
/* > \verbatim */
/* > JOBV1T is CHARACTER */
/* > = 'Y': V1T is updated;
*/
/* > otherwise: V1T is not updated. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV2T */
/* > \verbatim */
/* > JOBV2T is CHARACTER */
/* > = 'Y': V2T is updated;
*/
/* > otherwise: V2T is not updated. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER */
/* > = 'T': X, U1, U2, V1T, and V2T are stored in row-major */
/* > order;
*/
/* > otherwise: X, U1, U2, V1T, and V2T are stored in column- */
/* > major order. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows and columns in X, the unitary matrix in */
/* > bidiagonal-block form. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* > P is INTEGER */
/* > The number of rows in the top-left block of X. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* > Q is INTEGER */
/* > The number of columns in the top-left block of X. */
/* > 0 <= Q <= MIN(P,M-P,M-Q). */
/* > \endverbatim */
/* > */
/* > \param[in,out] THETA */
/* > \verbatim */
/* > THETA is REAL array, dimension (Q) */
/* > On entry, the angles THETA(1),...,THETA(Q) that, along with */
/* > PHI(1), ...,PHI(Q-1), define the matrix in bidiagonal-block */
/* > form. On exit, the angles whose cosines and sines define the */
/* > diagonal blocks in the CS decomposition. */
/* > \endverbatim */
/* > */
/* > \param[in,out] PHI */
/* > \verbatim */
/* > PHI is REAL array, dimension (Q-1) */
/* > The angles PHI(1),...,PHI(Q-1) that, along with THETA(1),..., */
/* > THETA(Q), define the matrix in bidiagonal-block form. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U1 */
/* > \verbatim */
/* > U1 is COMPLEX array, dimension (LDU1,P) */
/* > On entry, an LDU1-by-P matrix. On exit, U1 is postmultiplied */
/* > by the left singular vector matrix common to [ B11 ;
0 ] and */
/* > [ B12 0 0 ;
0 -I 0 0 ]. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU1 */
/* > \verbatim */
/* > LDU1 is INTEGER */
/* > The leading dimension of the array U1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U2 */
/* > \verbatim */
/* > U2 is COMPLEX array, dimension (LDU2,M-P) */
/* > On entry, an LDU2-by-(M-P) matrix. On exit, U2 is */
/* > postmultiplied by the left singular vector matrix common to */
/* > [ B21 ;
0 ] and [ B22 0 0 ;
0 0 I ]. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU2 */
/* > \verbatim */
/* > LDU2 is INTEGER */
/* > The leading dimension of the array U2. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V1T */
/* > \verbatim */
/* > V1T is COMPLEX array, dimension (LDV1T,Q) */
/* > On entry, a LDV1T-by-Q matrix. On exit, V1T is premultiplied */
/* > by the conjugate transpose of the right singular vector */
/* > matrix common to [ B11 ;
0 ] and [ B21 ;
0 ]. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV1T */
/* > \verbatim */
/* > LDV1T is INTEGER */
/* > The leading dimension of the array V1T. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V2T */
/* > \verbatim */
/* > V2T is COMPLEX array, dimenison (LDV2T,M-Q) */
/* > On entry, a LDV2T-by-(M-Q) matrix. On exit, V2T is */
/* > premultiplied by the conjugate transpose of the right */
/* > singular vector matrix common to [ B12 0 0 ;
0 -I 0 ] and */
/* > [ B22 0 0 ;
0 0 I ]. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV2T */
/* > \verbatim */
/* > LDV2T is INTEGER */
/* > The leading dimension of the array V2T. */
/* > \endverbatim */
/* > */
/* > \param[out] B11D */
/* > \verbatim */
/* > B11D is REAL array, dimension (Q) */
/* > When CBBCSD converges, B11D contains the cosines of THETA(1), */
/* > ..., THETA(Q). If CBBCSD fails to converge, then B11D */
/* > contains the diagonal of the partially reduced top-left */
/* > block. */
/* > \endverbatim */
/* > */
/* > \param[out] B11E */
/* > \verbatim */
/* > B11E is REAL array, dimension (Q-1) */
/* > When CBBCSD converges, B11E contains zeros. If CBBCSD fails */
/* > to converge, then B11E contains the superdiagonal of the */
/* > partially reduced top-left block. */
/* > \endverbatim */
/* > */
/* > \param[out] B12D */
/* > \verbatim */
/* > B12D is REAL array, dimension (Q) */
/* > When CBBCSD converges, B12D contains the negative sines of */
/* > THETA(1), ..., THETA(Q). If CBBCSD fails to converge, then */
/* > B12D contains the diagonal of the partially reduced top-right */
/* > block. */
/* > \endverbatim */
/* > */
/* > \param[out] B12E */
/* > \verbatim */
/* > B12E is REAL array, dimension (Q-1) */
/* > When CBBCSD converges, B12E contains zeros. If CBBCSD fails */
/* > to converge, then B12E contains the subdiagonal of the */
/* > partially reduced top-right block. */
/* > \endverbatim */
/* > */
/* > \param[out] B21D */
/* > \verbatim */
/* > B21D is REAL array, dimension (Q) */
/* > When CBBCSD converges, B21D contains the negative sines of */
/* > THETA(1), ..., THETA(Q). If CBBCSD fails to converge, then */
/* > B21D contains the diagonal of the partially reduced bottom-left */
/* > block. */
/* > \endverbatim */
/* > */
/* > \param[out] B21E */
/* > \verbatim */
/* > B21E is REAL array, dimension (Q-1) */
/* > When CBBCSD converges, B21E contains zeros. If CBBCSD fails */
/* > to converge, then B21E contains the subdiagonal of the */
/* > partially reduced bottom-left block. */
/* > \endverbatim */
/* > */
/* > \param[out] B22D */
/* > \verbatim */
/* > B22D is REAL array, dimension (Q) */
/* > When CBBCSD converges, B22D contains the negative sines of */
/* > THETA(1), ..., THETA(Q). If CBBCSD fails to converge, then */
/* > B22D contains the diagonal of the partially reduced bottom-right */
/* > block. */
/* > \endverbatim */
/* > */
/* > \param[out] B22E */
/* > \verbatim */
/* > B22E is REAL array, dimension (Q-1) */
/* > When CBBCSD converges, B22E contains zeros. If CBBCSD fails */
/* > to converge, then B22E contains the subdiagonal of the */
/* > partially reduced bottom-right block. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* > LRWORK is INTEGER */
/* > The dimension of the array RWORK. LRWORK >= MAX(1,8*Q). */
/* > */
/* > If LRWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal size of the RWORK array, */
/* > returns this value as the first entry of the work array, and */
/* > no error message related to LRWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if CBBCSD did not converge, INFO specifies the number */
/* > of nonzero entries in PHI, and B11D, B11E, etc., */
/* > contain the partially reduced matrix. */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > TOLMUL REAL, default = MAX(10,MIN(100,EPS**(-1/8))) */
/* > TOLMUL controls the convergence criterion of the QR loop. */
/* > Angles THETA(i), PHI(i) are rounded to 0 or PI/2 when they */
/* > are within TOLMUL*EPS of either bound. */
/* > \endverbatim */
/* > \par References: */
/* ================ */
/* > */
/* > [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* > Algorithms, 50(1):33-65, 2009. */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int cbbcsd_(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, integer *m, integer *p, integer *q, real *theta, real *phi, complex *u1, integer *ldu1, complex *u2, integer *ldu2, complex *v1t, integer *ldv1t, complex *v2t, integer *ldv2t, real * b11d, real *b11e, real *b12d, real *b12e, real *b21d, real *b21e, real *b22d, real *b22e, real *rwork, integer *lrwork, integer *info)
{
    /* System generated locals */
    integer u1_dim1, u1_offset, u2_dim1, u2_offset, v1t_dim1, v1t_offset, v2t_dim1, v2t_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4;
    doublereal d__1;
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), cos(doublereal), sin( doublereal), sqrt(doublereal), atan2(doublereal, doublereal);
    /* Local variables */
    logical colmajor;
    real thetamin, thetamax;
    logical restart11, restart12, restart21, restart22;
    integer i__, j;
    real r__, x1, x2, y1, y2;
    integer lrworkmin, lrworkopt;
    real mu, nu, eps, tol;
    integer imin, mini, imax, iter;
    real unfl, temp;
    integer iu1cs, iu2cs;
    extern /* Subroutine */
    int slas2_(real *, real *, real *, real *, real *) ;
    integer iu1sn, iu2sn;
    extern /* Subroutine */
    int cscal_(integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int clasr_(char *, char *, char *, integer *, integer *, real *, real *, complex *, integer *), cswap_(integer *, complex *, integer *, complex *, integer *);
    integer maxit;
    real dummy;
    integer iv1tcs, iv2tcs;
    logical wantu1, wantu2;
    integer iv1tsn, iv2tsn;
    real sigma11, sigma21;
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real thresh, tolmul;
    logical lquery;
    real b11bulge, b12bulge;
    logical wantv1t, wantv2t;
    real b21bulge, b22bulge;
    extern /* Subroutine */
    int slartgp_(real *, real *, real *, real *, real *), slartgs_(real *, real *, real *, real *, real *);
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* =================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test input arguments */
    /* Parameter adjustments */
    --theta;
    --phi;
    u1_dim1 = *ldu1;
    u1_offset = 1 + u1_dim1;
    u1 -= u1_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    v1t_dim1 = *ldv1t;
    v1t_offset = 1 + v1t_dim1;
    v1t -= v1t_offset;
    v2t_dim1 = *ldv2t;
    v2t_offset = 1 + v2t_dim1;
    v2t -= v2t_offset;
    --b11d;
    --b11e;
    --b12d;
    --b12e;
    --b21d;
    --b21e;
    --b22d;
    --b22e;
    --rwork;
    /* Function Body */
    *info = 0;
    lquery = *lrwork == -1;
    wantu1 = lsame_(jobu1, "Y");
    wantu2 = lsame_(jobu2, "Y");
    wantv1t = lsame_(jobv1t, "Y");
    wantv2t = lsame_(jobv2t, "Y");
    colmajor = ! lsame_(trans, "T");
    if (*m < 0)
    {
        *info = -6;
    }
    else if (*p < 0 || *p > *m)
    {
        *info = -7;
    }
    else if (*q < 0 || *q > *m)
    {
        *info = -8;
    }
    else if (*q > *p || *q > *m - *p || *q > *m - *q)
    {
        *info = -8;
    }
    else if (wantu1 && *ldu1 < *p)
    {
        *info = -12;
    }
    else if (wantu2 && *ldu2 < *m - *p)
    {
        *info = -14;
    }
    else if (wantv1t && *ldv1t < *q)
    {
        *info = -16;
    }
    else if (wantv2t && *ldv2t < *m - *q)
    {
        *info = -18;
    }
    /* Quick return if Q = 0 */
    if (*info == 0 && *q == 0)
    {
        lrworkmin = 1;
        rwork[1] = (real) lrworkmin;
        return 0;
    }
    /* Compute workspace */
    if (*info == 0)
    {
        iu1cs = 1;
        iu1sn = iu1cs + *q;
        iu2cs = iu1sn + *q;
        iu2sn = iu2cs + *q;
        iv1tcs = iu2sn + *q;
        iv1tsn = iv1tcs + *q;
        iv2tcs = iv1tsn + *q;
        iv2tsn = iv2tcs + *q;
        lrworkopt = iv2tsn + *q - 1;
        lrworkmin = lrworkopt;
        rwork[1] = (real) lrworkopt;
        if (*lrwork < lrworkmin && ! lquery)
        {
            *info = -28;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CBBCSD", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Get machine constants */
    eps = slamch_("Epsilon");
    unfl = slamch_("Safe minimum");
    /* Computing MAX */
    /* Computing MIN */
    d__1 = (doublereal) eps;
    r__3 = 100.f;
    r__4 = pow_dd(&d__1, &c_b11); // , expr subst
    r__1 = 10.f;
    r__2 = min(r__3,r__4); // , expr subst
    tolmul = max(r__1,r__2);
    tol = tolmul * eps;
    /* Computing MAX */
    r__1 = tol;
    r__2 = *q * 6 * *q * unfl; // , expr subst
    thresh = max(r__1,r__2);
    /* Test for negligible sines or cosines */
    i__1 = *q;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if (theta[i__] < thresh)
        {
            theta[i__] = 0.f;
        }
        else if (theta[i__] > 1.57079632679489662f - thresh)
        {
            theta[i__] = 1.57079632679489662f;
        }
    }
    i__1 = *q - 1;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if (phi[i__] < thresh)
        {
            phi[i__] = 0.f;
        }
        else if (phi[i__] > 1.57079632679489662f - thresh)
        {
            phi[i__] = 1.57079632679489662f;
        }
    }
    /* Initial deflation */
    imax = *q;
    while(imax > 1)
    {
        if (phi[imax - 1] != 0.f)
        {
            break;
        }
        --imax;
    }
    imin = imax - 1;
    if (imin > 1)
    {
        while(phi[imin - 1] != 0.f)
        {
            --imin;
            if (imin <= 1)
            {
                break;
            }
        }
    }
    /* Initialize iteration counter */
    maxit = *q * 6 * *q;
    iter = 0;
    /* Begin main iteration loop */
    while(imax > 1)
    {
        /* Compute the matrix entries */
        b11d[imin] = cos(theta[imin]);
        b21d[imin] = -sin(theta[imin]);
        i__1 = imax - 1;
        for (i__ = imin;
                i__ <= i__1;
                ++i__)
        {
            b11e[i__] = -sin(theta[i__]) * sin(phi[i__]);
            b11d[i__ + 1] = cos(theta[i__ + 1]) * cos(phi[i__]);
            b12d[i__] = sin(theta[i__]) * cos(phi[i__]);
            b12e[i__] = cos(theta[i__ + 1]) * sin(phi[i__]);
            b21e[i__] = -cos(theta[i__]) * sin(phi[i__]);
            b21d[i__ + 1] = -sin(theta[i__ + 1]) * cos(phi[i__]);
            b22d[i__] = cos(theta[i__]) * cos(phi[i__]);
            b22e[i__] = -sin(theta[i__ + 1]) * sin(phi[i__]);
        }
        b12d[imax] = sin(theta[imax]);
        b22d[imax] = cos(theta[imax]);
        /* Abort if not converging;
        otherwise, increment ITER */
        if (iter > maxit)
        {
            *info = 0;
            i__1 = *q;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                if (phi[i__] != 0.f)
                {
                    ++(*info);
                }
            }
            return 0;
        }
        iter = iter + imax - imin;
        /* Compute shifts */
        thetamax = theta[imin];
        thetamin = theta[imin];
        i__1 = imax;
        for (i__ = imin + 1;
                i__ <= i__1;
                ++i__)
        {
            if (theta[i__] > thetamax)
            {
                thetamax = theta[i__];
            }
            if (theta[i__] < thetamin)
            {
                thetamin = theta[i__];
            }
        }
        if (thetamax > 1.57079632679489662f - thresh)
        {
            /* Zero on diagonals of B11 and B22;
            induce deflation with a */
            /* zero shift */
            mu = 0.f;
            nu = 1.f;
        }
        else if (thetamin < thresh)
        {
            /* Zero on diagonals of B12 and B22;
            induce deflation with a */
            /* zero shift */
            mu = 1.f;
            nu = 0.f;
        }
        else
        {
            /* Compute shifts for B11 and B21 and use the lesser */
            slas2_(&b11d[imax - 1], &b11e[imax - 1], &b11d[imax], &sigma11, & dummy);
            slas2_(&b21d[imax - 1], &b21e[imax - 1], &b21d[imax], &sigma21, & dummy);
            if (sigma11 <= sigma21)
            {
                mu = sigma11;
                /* Computing 2nd power */
                r__1 = mu;
                nu = sqrt(1.f - r__1 * r__1);
                if (mu < thresh)
                {
                    mu = 0.f;
                    nu = 1.f;
                }
            }
            else
            {
                nu = sigma21;
                /* Computing 2nd power */
                r__1 = nu;
                mu = sqrt(1.f - r__1 * r__1);
                if (nu < thresh)
                {
                    mu = 1.f;
                    nu = 0.f;
                }
            }
        }
        /* Rotate to produce bulges in B11 and B21 */
        if (mu <= nu)
        {
            slartgs_(&b11d[imin], &b11e[imin], &mu, &rwork[iv1tcs + imin - 1], &rwork[iv1tsn + imin - 1]);
        }
        else
        {
            slartgs_(&b21d[imin], &b21e[imin], &nu, &rwork[iv1tcs + imin - 1], &rwork[iv1tsn + imin - 1]);
        }
        temp = rwork[iv1tcs + imin - 1] * b11d[imin] + rwork[iv1tsn + imin - 1] * b11e[imin];
        b11e[imin] = rwork[iv1tcs + imin - 1] * b11e[imin] - rwork[iv1tsn + imin - 1] * b11d[imin];
        b11d[imin] = temp;
        b11bulge = rwork[iv1tsn + imin - 1] * b11d[imin + 1];
        b11d[imin + 1] = rwork[iv1tcs + imin - 1] * b11d[imin + 1];
        temp = rwork[iv1tcs + imin - 1] * b21d[imin] + rwork[iv1tsn + imin - 1] * b21e[imin];
        b21e[imin] = rwork[iv1tcs + imin - 1] * b21e[imin] - rwork[iv1tsn + imin - 1] * b21d[imin];
        b21d[imin] = temp;
        b21bulge = rwork[iv1tsn + imin - 1] * b21d[imin + 1];
        b21d[imin + 1] = rwork[iv1tcs + imin - 1] * b21d[imin + 1];
        /* Compute THETA(IMIN) */
        /* Computing 2nd power */
        r__1 = b21d[imin];
        /* Computing 2nd power */
        r__2 = b21bulge;
        /* Computing 2nd power */
        r__3 = b11d[imin];
        /* Computing 2nd power */
        r__4 = b11bulge;
        theta[imin] = atan2(sqrt(r__1 * r__1 + r__2 * r__2), sqrt(r__3 * r__3 + r__4 * r__4));
        /* Chase the bulges in B11(IMIN+1,IMIN) and B21(IMIN+1,IMIN) */
        /* Computing 2nd power */
        r__1 = b11d[imin];
        /* Computing 2nd power */
        r__2 = b11bulge;
        /* Computing 2nd power */
        r__3 = thresh;
        if (r__1 * r__1 + r__2 * r__2 > r__3 * r__3)
        {
            slartgp_(&b11bulge, &b11d[imin], &rwork[iu1sn + imin - 1], &rwork[ iu1cs + imin - 1], &r__);
        }
        else if (mu <= nu)
        {
            slartgs_(&b11e[imin], &b11d[imin + 1], &mu, &rwork[iu1cs + imin - 1], &rwork[iu1sn + imin - 1]);
        }
        else
        {
            slartgs_(&b12d[imin], &b12e[imin], &nu, &rwork[iu1cs + imin - 1], &rwork[iu1sn + imin - 1]);
        }
        /* Computing 2nd power */
        r__1 = b21d[imin];
        /* Computing 2nd power */
        r__2 = b21bulge;
        /* Computing 2nd power */
        r__3 = thresh;
        if (r__1 * r__1 + r__2 * r__2 > r__3 * r__3)
        {
            slartgp_(&b21bulge, &b21d[imin], &rwork[iu2sn + imin - 1], &rwork[ iu2cs + imin - 1], &r__);
        }
        else if (nu < mu)
        {
            slartgs_(&b21e[imin], &b21d[imin + 1], &nu, &rwork[iu2cs + imin - 1], &rwork[iu2sn + imin - 1]);
        }
        else
        {
            slartgs_(&b22d[imin], &b22e[imin], &mu, &rwork[iu2cs + imin - 1], &rwork[iu2sn + imin - 1]);
        }
        rwork[iu2cs + imin - 1] = -rwork[iu2cs + imin - 1];
        rwork[iu2sn + imin - 1] = -rwork[iu2sn + imin - 1];
        temp = rwork[iu1cs + imin - 1] * b11e[imin] + rwork[iu1sn + imin - 1] * b11d[imin + 1];
        b11d[imin + 1] = rwork[iu1cs + imin - 1] * b11d[imin + 1] - rwork[ iu1sn + imin - 1] * b11e[imin];
        b11e[imin] = temp;
        if (imax > imin + 1)
        {
            b11bulge = rwork[iu1sn + imin - 1] * b11e[imin + 1];
            b11e[imin + 1] = rwork[iu1cs + imin - 1] * b11e[imin + 1];
        }
        temp = rwork[iu1cs + imin - 1] * b12d[imin] + rwork[iu1sn + imin - 1] * b12e[imin];
        b12e[imin] = rwork[iu1cs + imin - 1] * b12e[imin] - rwork[iu1sn + imin - 1] * b12d[imin];
        b12d[imin] = temp;
        b12bulge = rwork[iu1sn + imin - 1] * b12d[imin + 1];
        b12d[imin + 1] = rwork[iu1cs + imin - 1] * b12d[imin + 1];
        temp = rwork[iu2cs + imin - 1] * b21e[imin] + rwork[iu2sn + imin - 1] * b21d[imin + 1];
        b21d[imin + 1] = rwork[iu2cs + imin - 1] * b21d[imin + 1] - rwork[ iu2sn + imin - 1] * b21e[imin];
        b21e[imin] = temp;
        if (imax > imin + 1)
        {
            b21bulge = rwork[iu2sn + imin - 1] * b21e[imin + 1];
            b21e[imin + 1] = rwork[iu2cs + imin - 1] * b21e[imin + 1];
        }
        temp = rwork[iu2cs + imin - 1] * b22d[imin] + rwork[iu2sn + imin - 1] * b22e[imin];
        b22e[imin] = rwork[iu2cs + imin - 1] * b22e[imin] - rwork[iu2sn + imin - 1] * b22d[imin];
        b22d[imin] = temp;
        b22bulge = rwork[iu2sn + imin - 1] * b22d[imin + 1];
        b22d[imin + 1] = rwork[iu2cs + imin - 1] * b22d[imin + 1];
        /* Inner loop: chase bulges from B11(IMIN,IMIN+2), */
        /* B12(IMIN,IMIN+1), B21(IMIN,IMIN+2), and B22(IMIN,IMIN+1) to */
        /* bottom-right */
        i__1 = imax - 1;
        for (i__ = imin + 1;
                i__ <= i__1;
                ++i__)
        {
            /* Compute PHI(I-1) */
            x1 = sin(theta[i__ - 1]) * b11e[i__ - 1] + cos(theta[i__ - 1]) * b21e[i__ - 1];
            x2 = sin(theta[i__ - 1]) * b11bulge + cos(theta[i__ - 1]) * b21bulge;
            y1 = sin(theta[i__ - 1]) * b12d[i__ - 1] + cos(theta[i__ - 1]) * b22d[i__ - 1];
            y2 = sin(theta[i__ - 1]) * b12bulge + cos(theta[i__ - 1]) * b22bulge;
            /* Computing 2nd power */
            r__1 = x1;
            /* Computing 2nd power */
            r__2 = x2;
            /* Computing 2nd power */
            r__3 = y1;
            /* Computing 2nd power */
            r__4 = y2;
            phi[i__ - 1] = atan2(sqrt(r__1 * r__1 + r__2 * r__2), sqrt(r__3 * r__3 + r__4 * r__4));
            /* Determine if there are bulges to chase or if a new direct */
            /* summand has been reached */
            /* Computing 2nd power */
            r__1 = b11e[i__ - 1];
            /* Computing 2nd power */
            r__2 = b11bulge;
            /* Computing 2nd power */
            r__3 = thresh;
            restart11 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
            /* Computing 2nd power */
            r__1 = b21e[i__ - 1];
            /* Computing 2nd power */
            r__2 = b21bulge;
            /* Computing 2nd power */
            r__3 = thresh;
            restart21 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
            /* Computing 2nd power */
            r__1 = b12d[i__ - 1];
            /* Computing 2nd power */
            r__2 = b12bulge;
            /* Computing 2nd power */
            r__3 = thresh;
            restart12 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
            /* Computing 2nd power */
            r__1 = b22d[i__ - 1];
            /* Computing 2nd power */
            r__2 = b22bulge;
            /* Computing 2nd power */
            r__3 = thresh;
            restart22 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
            /* If possible, chase bulges from B11(I-1,I+1), B12(I-1,I), */
            /* B21(I-1,I+1), and B22(I-1,I). If necessary, restart bulge- */
            /* chasing by applying the original shift again. */
            if (! restart11 && ! restart21)
            {
                slartgp_(&x2, &x1, &rwork[iv1tsn + i__ - 1], &rwork[iv1tcs + i__ - 1], &r__);
            }
            else if (! restart11 && restart21)
            {
                slartgp_(&b11bulge, &b11e[i__ - 1], &rwork[iv1tsn + i__ - 1], &rwork[iv1tcs + i__ - 1], &r__);
            }
            else if (restart11 && ! restart21)
            {
                slartgp_(&b21bulge, &b21e[i__ - 1], &rwork[iv1tsn + i__ - 1], &rwork[iv1tcs + i__ - 1], &r__);
            }
            else if (mu <= nu)
            {
                slartgs_(&b11d[i__], &b11e[i__], &mu, &rwork[iv1tcs + i__ - 1] , &rwork[iv1tsn + i__ - 1]);
            }
            else
            {
                slartgs_(&b21d[i__], &b21e[i__], &nu, &rwork[iv1tcs + i__ - 1] , &rwork[iv1tsn + i__ - 1]);
            }
            rwork[iv1tcs + i__ - 1] = -rwork[iv1tcs + i__ - 1];
            rwork[iv1tsn + i__ - 1] = -rwork[iv1tsn + i__ - 1];
            if (! restart12 && ! restart22)
            {
                slartgp_(&y2, &y1, &rwork[iv2tsn + i__ - 2], &rwork[iv2tcs + i__ - 2], &r__);
            }
            else if (! restart12 && restart22)
            {
                slartgp_(&b12bulge, &b12d[i__ - 1], &rwork[iv2tsn + i__ - 2], &rwork[iv2tcs + i__ - 2], &r__);
            }
            else if (restart12 && ! restart22)
            {
                slartgp_(&b22bulge, &b22d[i__ - 1], &rwork[iv2tsn + i__ - 2], &rwork[iv2tcs + i__ - 2], &r__);
            }
            else if (nu < mu)
            {
                slartgs_(&b12e[i__ - 1], &b12d[i__], &nu, &rwork[iv2tcs + i__ - 2], &rwork[iv2tsn + i__ - 2]);
            }
            else
            {
                slartgs_(&b22e[i__ - 1], &b22d[i__], &mu, &rwork[iv2tcs + i__ - 2], &rwork[iv2tsn + i__ - 2]);
            }
            temp = rwork[iv1tcs + i__ - 1] * b11d[i__] + rwork[iv1tsn + i__ - 1] * b11e[i__];
            b11e[i__] = rwork[iv1tcs + i__ - 1] * b11e[i__] - rwork[iv1tsn + i__ - 1] * b11d[i__];
            b11d[i__] = temp;
            b11bulge = rwork[iv1tsn + i__ - 1] * b11d[i__ + 1];
            b11d[i__ + 1] = rwork[iv1tcs + i__ - 1] * b11d[i__ + 1];
            temp = rwork[iv1tcs + i__ - 1] * b21d[i__] + rwork[iv1tsn + i__ - 1] * b21e[i__];
            b21e[i__] = rwork[iv1tcs + i__ - 1] * b21e[i__] - rwork[iv1tsn + i__ - 1] * b21d[i__];
            b21d[i__] = temp;
            b21bulge = rwork[iv1tsn + i__ - 1] * b21d[i__ + 1];
            b21d[i__ + 1] = rwork[iv1tcs + i__ - 1] * b21d[i__ + 1];
            temp = rwork[iv2tcs + i__ - 2] * b12e[i__ - 1] + rwork[iv2tsn + i__ - 2] * b12d[i__];
            b12d[i__] = rwork[iv2tcs + i__ - 2] * b12d[i__] - rwork[iv2tsn + i__ - 2] * b12e[i__ - 1];
            b12e[i__ - 1] = temp;
            b12bulge = rwork[iv2tsn + i__ - 2] * b12e[i__];
            b12e[i__] = rwork[iv2tcs + i__ - 2] * b12e[i__];
            temp = rwork[iv2tcs + i__ - 2] * b22e[i__ - 1] + rwork[iv2tsn + i__ - 2] * b22d[i__];
            b22d[i__] = rwork[iv2tcs + i__ - 2] * b22d[i__] - rwork[iv2tsn + i__ - 2] * b22e[i__ - 1];
            b22e[i__ - 1] = temp;
            b22bulge = rwork[iv2tsn + i__ - 2] * b22e[i__];
            b22e[i__] = rwork[iv2tcs + i__ - 2] * b22e[i__];
            /* Compute THETA(I) */
            x1 = cos(phi[i__ - 1]) * b11d[i__] + sin(phi[i__ - 1]) * b12e[i__ - 1];
            x2 = cos(phi[i__ - 1]) * b11bulge + sin(phi[i__ - 1]) * b12bulge;
            y1 = cos(phi[i__ - 1]) * b21d[i__] + sin(phi[i__ - 1]) * b22e[i__ - 1];
            y2 = cos(phi[i__ - 1]) * b21bulge + sin(phi[i__ - 1]) * b22bulge;
            /* Computing 2nd power */
            r__1 = y1;
            /* Computing 2nd power */
            r__2 = y2;
            /* Computing 2nd power */
            r__3 = x1;
            /* Computing 2nd power */
            r__4 = x2;
            theta[i__] = atan2(sqrt(r__1 * r__1 + r__2 * r__2), sqrt(r__3 * r__3 + r__4 * r__4));
            /* Determine if there are bulges to chase or if a new direct */
            /* summand has been reached */
            /* Computing 2nd power */
            r__1 = b11d[i__];
            /* Computing 2nd power */
            r__2 = b11bulge;
            /* Computing 2nd power */
            r__3 = thresh;
            restart11 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
            /* Computing 2nd power */
            r__1 = b12e[i__ - 1];
            /* Computing 2nd power */
            r__2 = b12bulge;
            /* Computing 2nd power */
            r__3 = thresh;
            restart12 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
            /* Computing 2nd power */
            r__1 = b21d[i__];
            /* Computing 2nd power */
            r__2 = b21bulge;
            /* Computing 2nd power */
            r__3 = thresh;
            restart21 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
            /* Computing 2nd power */
            r__1 = b22e[i__ - 1];
            /* Computing 2nd power */
            r__2 = b22bulge;
            /* Computing 2nd power */
            r__3 = thresh;
            restart22 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
            /* If possible, chase bulges from B11(I+1,I), B12(I+1,I-1), */
            /* B21(I+1,I), and B22(I+1,I-1). If necessary, restart bulge- */
            /* chasing by applying the original shift again. */
            if (! restart11 && ! restart12)
            {
                slartgp_(&x2, &x1, &rwork[iu1sn + i__ - 1], &rwork[iu1cs + i__ - 1], &r__);
            }
            else if (! restart11 && restart12)
            {
                slartgp_(&b11bulge, &b11d[i__], &rwork[iu1sn + i__ - 1], & rwork[iu1cs + i__ - 1], &r__);
            }
            else if (restart11 && ! restart12)
            {
                slartgp_(&b12bulge, &b12e[i__ - 1], &rwork[iu1sn + i__ - 1], & rwork[iu1cs + i__ - 1], &r__);
            }
            else if (mu <= nu)
            {
                slartgs_(&b11e[i__], &b11d[i__ + 1], &mu, &rwork[iu1cs + i__ - 1], &rwork[iu1sn + i__ - 1]);
            }
            else
            {
                slartgs_(&b12d[i__], &b12e[i__], &nu, &rwork[iu1cs + i__ - 1], &rwork[iu1sn + i__ - 1]);
            }
            if (! restart21 && ! restart22)
            {
                slartgp_(&y2, &y1, &rwork[iu2sn + i__ - 1], &rwork[iu2cs + i__ - 1], &r__);
            }
            else if (! restart21 && restart22)
            {
                slartgp_(&b21bulge, &b21d[i__], &rwork[iu2sn + i__ - 1], & rwork[iu2cs + i__ - 1], &r__);
            }
            else if (restart21 && ! restart22)
            {
                slartgp_(&b22bulge, &b22e[i__ - 1], &rwork[iu2sn + i__ - 1], & rwork[iu2cs + i__ - 1], &r__);
            }
            else if (nu < mu)
            {
                slartgs_(&b21e[i__], &b21e[i__ + 1], &nu, &rwork[iu2cs + i__ - 1], &rwork[iu2sn + i__ - 1]);
            }
            else
            {
                slartgs_(&b22d[i__], &b22e[i__], &mu, &rwork[iu2cs + i__ - 1], &rwork[iu2sn + i__ - 1]);
            }
            rwork[iu2cs + i__ - 1] = -rwork[iu2cs + i__ - 1];
            rwork[iu2sn + i__ - 1] = -rwork[iu2sn + i__ - 1];
            temp = rwork[iu1cs + i__ - 1] * b11e[i__] + rwork[iu1sn + i__ - 1] * b11d[i__ + 1];
            b11d[i__ + 1] = rwork[iu1cs + i__ - 1] * b11d[i__ + 1] - rwork[ iu1sn + i__ - 1] * b11e[i__];
            b11e[i__] = temp;
            if (i__ < imax - 1)
            {
                b11bulge = rwork[iu1sn + i__ - 1] * b11e[i__ + 1];
                b11e[i__ + 1] = rwork[iu1cs + i__ - 1] * b11e[i__ + 1];
            }
            temp = rwork[iu2cs + i__ - 1] * b21e[i__] + rwork[iu2sn + i__ - 1] * b21d[i__ + 1];
            b21d[i__ + 1] = rwork[iu2cs + i__ - 1] * b21d[i__ + 1] - rwork[ iu2sn + i__ - 1] * b21e[i__];
            b21e[i__] = temp;
            if (i__ < imax - 1)
            {
                b21bulge = rwork[iu2sn + i__ - 1] * b21e[i__ + 1];
                b21e[i__ + 1] = rwork[iu2cs + i__ - 1] * b21e[i__ + 1];
            }
            temp = rwork[iu1cs + i__ - 1] * b12d[i__] + rwork[iu1sn + i__ - 1] * b12e[i__];
            b12e[i__] = rwork[iu1cs + i__ - 1] * b12e[i__] - rwork[iu1sn + i__ - 1] * b12d[i__];
            b12d[i__] = temp;
            b12bulge = rwork[iu1sn + i__ - 1] * b12d[i__ + 1];
            b12d[i__ + 1] = rwork[iu1cs + i__ - 1] * b12d[i__ + 1];
            temp = rwork[iu2cs + i__ - 1] * b22d[i__] + rwork[iu2sn + i__ - 1] * b22e[i__];
            b22e[i__] = rwork[iu2cs + i__ - 1] * b22e[i__] - rwork[iu2sn + i__ - 1] * b22d[i__];
            b22d[i__] = temp;
            b22bulge = rwork[iu2sn + i__ - 1] * b22d[i__ + 1];
            b22d[i__ + 1] = rwork[iu2cs + i__ - 1] * b22d[i__ + 1];
        }
        /* Compute PHI(IMAX-1) */
        x1 = sin(theta[imax - 1]) * b11e[imax - 1] + cos(theta[imax - 1]) * b21e[imax - 1];
        y1 = sin(theta[imax - 1]) * b12d[imax - 1] + cos(theta[imax - 1]) * b22d[imax - 1];
        y2 = sin(theta[imax - 1]) * b12bulge + cos(theta[imax - 1]) * b22bulge;
        /* Computing 2nd power */
        r__1 = y1;
        /* Computing 2nd power */
        r__2 = y2;
        phi[imax - 1] = atan2((f2c_abs(x1)), sqrt(r__1 * r__1 + r__2 * r__2));
        /* Chase bulges from B12(IMAX-1,IMAX) and B22(IMAX-1,IMAX) */
        /* Computing 2nd power */
        r__1 = b12d[imax - 1];
        /* Computing 2nd power */
        r__2 = b12bulge;
        /* Computing 2nd power */
        r__3 = thresh;
        restart12 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
        /* Computing 2nd power */
        r__1 = b22d[imax - 1];
        /* Computing 2nd power */
        r__2 = b22bulge;
        /* Computing 2nd power */
        r__3 = thresh;
        restart22 = r__1 * r__1 + r__2 * r__2 <= r__3 * r__3;
        if (! restart12 && ! restart22)
        {
            slartgp_(&y2, &y1, &rwork[iv2tsn + imax - 2], &rwork[iv2tcs + imax - 2], &r__);
        }
        else if (! restart12 && restart22)
        {
            slartgp_(&b12bulge, &b12d[imax - 1], &rwork[iv2tsn + imax - 2], & rwork[iv2tcs + imax - 2], &r__);
        }
        else if (restart12 && ! restart22)
        {
            slartgp_(&b22bulge, &b22d[imax - 1], &rwork[iv2tsn + imax - 2], & rwork[iv2tcs + imax - 2], &r__);
        }
        else if (nu < mu)
        {
            slartgs_(&b12e[imax - 1], &b12d[imax], &nu, &rwork[iv2tcs + imax - 2], &rwork[iv2tsn + imax - 2]);
        }
        else
        {
            slartgs_(&b22e[imax - 1], &b22d[imax], &mu, &rwork[iv2tcs + imax - 2], &rwork[iv2tsn + imax - 2]);
        }
        temp = rwork[iv2tcs + imax - 2] * b12e[imax - 1] + rwork[iv2tsn + imax - 2] * b12d[imax];
        b12d[imax] = rwork[iv2tcs + imax - 2] * b12d[imax] - rwork[iv2tsn + imax - 2] * b12e[imax - 1];
        b12e[imax - 1] = temp;
        temp = rwork[iv2tcs + imax - 2] * b22e[imax - 1] + rwork[iv2tsn + imax - 2] * b22d[imax];
        b22d[imax] = rwork[iv2tcs + imax - 2] * b22d[imax] - rwork[iv2tsn + imax - 2] * b22e[imax - 1];
        b22e[imax - 1] = temp;
        /* Update singular vectors */
        if (wantu1)
        {
            if (colmajor)
            {
                i__1 = imax - imin + 1;
                clasr_("R", "V", "F", p, &i__1, &rwork[iu1cs + imin - 1], & rwork[iu1sn + imin - 1], &u1[imin * u1_dim1 + 1], ldu1);
            }
            else
            {
                i__1 = imax - imin + 1;
                clasr_("L", "V", "F", &i__1, p, &rwork[iu1cs + imin - 1], & rwork[iu1sn + imin - 1], &u1[imin + u1_dim1], ldu1);
            }
        }
        if (wantu2)
        {
            if (colmajor)
            {
                i__1 = *m - *p;
                i__2 = imax - imin + 1;
                clasr_("R", "V", "F", &i__1, &i__2, &rwork[iu2cs + imin - 1], &rwork[iu2sn + imin - 1], &u2[imin * u2_dim1 + 1], ldu2);
            }
            else
            {
                i__1 = imax - imin + 1;
                i__2 = *m - *p;
                clasr_("L", "V", "F", &i__1, &i__2, &rwork[iu2cs + imin - 1], &rwork[iu2sn + imin - 1], &u2[imin + u2_dim1], ldu2);
            }
        }
        if (wantv1t)
        {
            if (colmajor)
            {
                i__1 = imax - imin + 1;
                clasr_("L", "V", "F", &i__1, q, &rwork[iv1tcs + imin - 1], & rwork[iv1tsn + imin - 1], &v1t[imin + v1t_dim1], ldv1t);
            }
            else
            {
                i__1 = imax - imin + 1;
                clasr_("R", "V", "F", q, &i__1, &rwork[iv1tcs + imin - 1], & rwork[iv1tsn + imin - 1], &v1t[imin * v1t_dim1 + 1], ldv1t);
            }
        }
        if (wantv2t)
        {
            if (colmajor)
            {
                i__1 = imax - imin + 1;
                i__2 = *m - *q;
                clasr_("L", "V", "F", &i__1, &i__2, &rwork[iv2tcs + imin - 1], &rwork[iv2tsn + imin - 1], &v2t[imin + v2t_dim1], ldv2t);
            }
            else
            {
                i__1 = *m - *q;
                i__2 = imax - imin + 1;
                clasr_("R", "V", "F", &i__1, &i__2, &rwork[iv2tcs + imin - 1], &rwork[iv2tsn + imin - 1], &v2t[imin * v2t_dim1 + 1], ldv2t);
            }
        }
        /* Fix signs on B11(IMAX-1,IMAX) and B21(IMAX-1,IMAX) */
        if (b11e[imax - 1] + b21e[imax - 1] > 0.f)
        {
            b11d[imax] = -b11d[imax];
            b21d[imax] = -b21d[imax];
            if (wantv1t)
            {
                if (colmajor)
                {
                    cscal_(q, &c_b1, &v1t[imax + v1t_dim1], ldv1t);
                }
                else
                {
                    cscal_(q, &c_b1, &v1t[imax * v1t_dim1 + 1], &c__1);
                }
            }
        }
        /* Compute THETA(IMAX) */
        x1 = cos(phi[imax - 1]) * b11d[imax] + sin(phi[imax - 1]) * b12e[imax - 1];
        y1 = cos(phi[imax - 1]) * b21d[imax] + sin(phi[imax - 1]) * b22e[imax - 1];
        theta[imax] = atan2((f2c_abs(y1)), (f2c_abs(x1)));
        /* Fix signs on B11(IMAX,IMAX), B12(IMAX,IMAX-1), B21(IMAX,IMAX), */
        /* and B22(IMAX,IMAX-1) */
        if (b11d[imax] + b12e[imax - 1] < 0.f)
        {
            b12d[imax] = -b12d[imax];
            if (wantu1)
            {
                if (colmajor)
                {
                    cscal_(p, &c_b1, &u1[imax * u1_dim1 + 1], &c__1);
                }
                else
                {
                    cscal_(p, &c_b1, &u1[imax + u1_dim1], ldu1);
                }
            }
        }
        if (b21d[imax] + b22e[imax - 1] > 0.f)
        {
            b22d[imax] = -b22d[imax];
            if (wantu2)
            {
                if (colmajor)
                {
                    i__1 = *m - *p;
                    cscal_(&i__1, &c_b1, &u2[imax * u2_dim1 + 1], &c__1);
                }
                else
                {
                    i__1 = *m - *p;
                    cscal_(&i__1, &c_b1, &u2[imax + u2_dim1], ldu2);
                }
            }
        }
        /* Fix signs on B12(IMAX,IMAX) and B22(IMAX,IMAX) */
        if (b12d[imax] + b22d[imax] < 0.f)
        {
            if (wantv2t)
            {
                if (colmajor)
                {
                    i__1 = *m - *q;
                    cscal_(&i__1, &c_b1, &v2t[imax + v2t_dim1], ldv2t);
                }
                else
                {
                    i__1 = *m - *q;
                    cscal_(&i__1, &c_b1, &v2t[imax * v2t_dim1 + 1], &c__1);
                }
            }
        }
        /* Test for negligible sines or cosines */
        i__1 = imax;
        for (i__ = imin;
                i__ <= i__1;
                ++i__)
        {
            if (theta[i__] < thresh)
            {
                theta[i__] = 0.f;
            }
            else if (theta[i__] > 1.57079632679489662f - thresh)
            {
                theta[i__] = 1.57079632679489662f;
            }
        }
        i__1 = imax - 1;
        for (i__ = imin;
                i__ <= i__1;
                ++i__)
        {
            if (phi[i__] < thresh)
            {
                phi[i__] = 0.f;
            }
            else if (phi[i__] > 1.57079632679489662f - thresh)
            {
                phi[i__] = 1.57079632679489662f;
            }
        }
        /* Deflate */
        if (imax > 1)
        {
            while(phi[imax - 1] == 0.f)
            {
                --imax;
                if (imax <= 1)
                {
                    break;
                }
            }
        }
        if (imin > imax - 1)
        {
            imin = imax - 1;
        }
        if (imin > 1)
        {
            while(phi[imin - 1] != 0.f)
            {
                --imin;
                if (imin <= 1)
                {
                    break;
                }
            }
        }
        /* Repeat main iteration loop */
    }
    /* Postprocessing: order THETA from least to greatest */
    i__1 = *q;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        mini = i__;
        thetamin = theta[i__];
        i__2 = *q;
        for (j = i__ + 1;
                j <= i__2;
                ++j)
        {
            if (theta[j] < thetamin)
            {
                mini = j;
                thetamin = theta[j];
            }
        }
        if (mini != i__)
        {
            theta[mini] = theta[i__];
            theta[i__] = thetamin;
            if (colmajor)
            {
                if (wantu1)
                {
                    cswap_(p, &u1[i__ * u1_dim1 + 1], &c__1, &u1[mini * u1_dim1 + 1], &c__1);
                }
                if (wantu2)
                {
                    i__2 = *m - *p;
                    cswap_(&i__2, &u2[i__ * u2_dim1 + 1], &c__1, &u2[mini * u2_dim1 + 1], &c__1);
                }
                if (wantv1t)
                {
                    cswap_(q, &v1t[i__ + v1t_dim1], ldv1t, &v1t[mini + v1t_dim1], ldv1t);
                }
                if (wantv2t)
                {
                    i__2 = *m - *q;
                    cswap_(&i__2, &v2t[i__ + v2t_dim1], ldv2t, &v2t[mini + v2t_dim1], ldv2t);
                }
            }
            else
            {
                if (wantu1)
                {
                    cswap_(p, &u1[i__ + u1_dim1], ldu1, &u1[mini + u1_dim1], ldu1);
                }
                if (wantu2)
                {
                    i__2 = *m - *p;
                    cswap_(&i__2, &u2[i__ + u2_dim1], ldu2, &u2[mini + u2_dim1], ldu2);
                }
                if (wantv1t)
                {
                    cswap_(q, &v1t[i__ * v1t_dim1 + 1], &c__1, &v1t[mini * v1t_dim1 + 1], &c__1);
                }
                if (wantv2t)
                {
                    i__2 = *m - *q;
                    cswap_(&i__2, &v2t[i__ * v2t_dim1 + 1], &c__1, &v2t[mini * v2t_dim1 + 1], &c__1);
                }
            }
        }
    }
    return 0;
    /* End of CBBCSD */
}
/* cbbcsd_ */
