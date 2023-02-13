/* ../netlib/dbbcsd.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b10 = -.125;
static doublereal c_b35 = -1.;
static integer c__1 = 1;
/* > \brief \b DBBCSD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DBBCSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbbcsd. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dbbcsd. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dbbcsd. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, */
/* THETA, PHI, U1, LDU1, U2, LDU2, V1T, LDV1T, */
/* V2T, LDV2T, B11D, B11E, B12D, B12E, B21D, B21E, */
/* B22D, B22E, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS */
/* INTEGER INFO, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION B11D( * ), B11E( * ), B12D( * ), B12E( * ), */
/* $ B21D( * ), B21E( * ), B22D( * ), B22E( * ), */
/* $ PHI( * ), THETA( * ), WORK( * ) */
/* DOUBLE PRECISION U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), */
/* $ V2T( LDV2T, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DBBCSD computes the CS decomposition of an orthogonal matrix in */
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
/* > [ U1 | ] [ 0 | 0 -I 0 ] [ V1 | ]**T */
/* > = [---------] [---------------] [---------] . */
/* > [ | U2 ] [ S | C 0 0 ] [ | V2 ] */
/* > [ 0 | 0 0 I ] */
/* > */
/* > X is M-by-M, its top-left block is P-by-Q, and Q must be no larger */
/* > than P, M-P, or M-Q. (If Q is not the smallest index, then X must be */
/* > transposed and/or permuted. This can be done in constant time using */
/* > the TRANS and SIGNS options. See DORCSD for details.) */
/* > */
/* > The bidiagonal matrices B11, B12, B21, and B22 are represented */
/* > implicitly by angles THETA(1:Q) and PHI(1:Q-1). */
/* > */
/* > The orthogonal matrices U1, U2, V1T, and V2T are input/output. */
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
/* > The number of rows and columns in X, the orthogonal matrix in */
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
/* > THETA is DOUBLE PRECISION array, dimension (Q) */
/* > On entry, the angles THETA(1),...,THETA(Q) that, along with */
/* > PHI(1), ...,PHI(Q-1), define the matrix in bidiagonal-block */
/* > form. On exit, the angles whose cosines and sines define the */
/* > diagonal blocks in the CS decomposition. */
/* > \endverbatim */
/* > */
/* > \param[in,out] PHI */
/* > \verbatim */
/* > PHI is DOUBLE PRECISION array, dimension (Q-1) */
/* > The angles PHI(1),...,PHI(Q-1) that, along with THETA(1),..., */
/* > THETA(Q), define the matrix in bidiagonal-block form. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U1 */
/* > \verbatim */
/* > U1 is DOUBLE PRECISION array, dimension (LDU1,P) */
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
/* > U2 is DOUBLE PRECISION array, dimension (LDU2,M-P) */
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
/* > V1T is DOUBLE PRECISION array, dimension (LDV1T,Q) */
/* > On entry, a LDV1T-by-Q matrix. On exit, V1T is premultiplied */
/* > by the transpose of the right singular vector */
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
/* > V2T is DOUBLE PRECISION array, dimenison (LDV2T,M-Q) */
/* > On entry, a LDV2T-by-(M-Q) matrix. On exit, V2T is */
/* > premultiplied by the transpose of the right */
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
/* > B11D is DOUBLE PRECISION array, dimension (Q) */
/* > When DBBCSD converges, B11D contains the cosines of THETA(1), */
/* > ..., THETA(Q). If DBBCSD fails to converge, then B11D */
/* > contains the diagonal of the partially reduced top-left */
/* > block. */
/* > \endverbatim */
/* > */
/* > \param[out] B11E */
/* > \verbatim */
/* > B11E is DOUBLE PRECISION array, dimension (Q-1) */
/* > When DBBCSD converges, B11E contains zeros. If DBBCSD fails */
/* > to converge, then B11E contains the superdiagonal of the */
/* > partially reduced top-left block. */
/* > \endverbatim */
/* > */
/* > \param[out] B12D */
/* > \verbatim */
/* > B12D is DOUBLE PRECISION array, dimension (Q) */
/* > When DBBCSD converges, B12D contains the negative sines of */
/* > THETA(1), ..., THETA(Q). If DBBCSD fails to converge, then */
/* > B12D contains the diagonal of the partially reduced top-right */
/* > block. */
/* > \endverbatim */
/* > */
/* > \param[out] B12E */
/* > \verbatim */
/* > B12E is DOUBLE PRECISION array, dimension (Q-1) */
/* > When DBBCSD converges, B12E contains zeros. If DBBCSD fails */
/* > to converge, then B12E contains the subdiagonal of the */
/* > partially reduced top-right block. */
/* > \endverbatim */
/* > */
/* > \param[out] B21D */
/* > \verbatim */
/* > B21D is DOUBLE PRECISION array, dimension (Q) */
/* > When CBBCSD converges, B21D contains the negative sines of */
/* > THETA(1), ..., THETA(Q). If CBBCSD fails to converge, then */
/* > B21D contains the diagonal of the partially reduced bottom-left */
/* > block. */
/* > \endverbatim */
/* > */
/* > \param[out] B21E */
/* > \verbatim */
/* > B21E is DOUBLE PRECISION array, dimension (Q-1) */
/* > When CBBCSD converges, B21E contains zeros. If CBBCSD fails */
/* > to converge, then B21E contains the subdiagonal of the */
/* > partially reduced bottom-left block. */
/* > \endverbatim */
/* > */
/* > \param[out] B22D */
/* > \verbatim */
/* > B22D is DOUBLE PRECISION array, dimension (Q) */
/* > When CBBCSD converges, B22D contains the negative sines of */
/* > THETA(1), ..., THETA(Q). If CBBCSD fails to converge, then */
/* > B22D contains the diagonal of the partially reduced bottom-right */
/* > block. */
/* > \endverbatim */
/* > */
/* > \param[out] B22E */
/* > \verbatim */
/* > B22E is DOUBLE PRECISION array, dimension (Q-1) */
/* > When CBBCSD converges, B22E contains zeros. If CBBCSD fails */
/* > to converge, then B22E contains the subdiagonal of the */
/* > partially reduced bottom-right block. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= MAX(1,8*Q). */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal size of the WORK array, */
/* > returns this value as the first entry of the work array, and */
/* > no error message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if DBBCSD did not converge, INFO specifies the number */
/* > of nonzero entries in PHI, and B11D, B11E, etc., */
/* > contain the partially reduced matrix. */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > TOLMUL DOUBLE PRECISION, default = MAX(10,MIN(100,EPS**(-1/8))) */
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
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int dbbcsd_(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, integer *m, integer *p, integer *q, doublereal * theta, doublereal *phi, doublereal *u1, integer *ldu1, doublereal *u2, integer *ldu2, doublereal *v1t, integer *ldv1t, doublereal *v2t, integer *ldv2t, doublereal *b11d, doublereal *b11e, doublereal *b12d, doublereal *b12e, doublereal *b21d, doublereal *b21e, doublereal * b22d, doublereal *b22e, doublereal *work, integer *lwork, integer * info)
{
    /* System generated locals */
    integer u1_dim1, u1_offset, u2_dim1, u2_offset, v1t_dim1, v1t_offset, v2t_dim1, v2t_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), cos(doublereal), sin( doublereal), sqrt(doublereal), atan2(doublereal, doublereal);
    /* Local variables */
    logical colmajor;
    doublereal thetamin, thetamax;
    logical restart11, restart12, restart21, restart22;
    integer lworkmin, lworkopt, i__, j;
    doublereal r__, x1, x2, y1, y2, mu, nu, eps, tol;
    integer imin, mini, imax, iter;
    doublereal unfl, temp;
    extern /* Subroutine */
    int dlas2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    integer iu1cs, iu2cs, iu1sn, iu2sn;
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int dlasr_(char *, char *, char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *), dswap_(integer *, doublereal *, integer * , doublereal *, integer *);
    integer maxit;
    doublereal dummy;
    integer iv1tcs, iv2tcs;
    logical wantu1, wantu2;
    integer iv1tsn, iv2tsn;
    extern doublereal dlamch_(char *);
    doublereal sigma11, sigma21;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal thresh, tolmul;
    logical lquery;
    doublereal b11bulge, b12bulge;
    logical wantv1t, wantv2t;
    doublereal b21bulge, b22bulge;
    extern /* Subroutine */
    int dlartgp_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *), dlartgs_( doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
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
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
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
        lworkmin = 1;
        work[1] = (doublereal) lworkmin;
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
        lworkopt = iv2tsn + *q - 1;
        lworkmin = lworkopt;
        work[1] = (doublereal) lworkopt;
        if (*lwork < lworkmin && ! lquery)
        {
            *info = -28;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DBBCSD", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Get machine constants */
    eps = dlamch_("Epsilon");
    unfl = dlamch_("Safe minimum");
    /* Computing MAX */
    /* Computing MIN */
    d__3 = 100.;
    d__4 = pow_dd(&eps, &c_b10); // , expr subst
    d__1 = 10.;
    d__2 = min(d__3,d__4); // , expr subst
    tolmul = max(d__1,d__2);
    tol = tolmul * eps;
    /* Computing MAX */
    d__1 = tol;
    d__2 = *q * 6 * *q * unfl; // , expr subst
    thresh = max(d__1,d__2);
    /* Test for negligible sines or cosines */
    i__1 = *q;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if (theta[i__] < thresh)
        {
            theta[i__] = 0.;
        }
        else if (theta[i__] > 1.57079632679489662 - thresh)
        {
            theta[i__] = 1.57079632679489662;
        }
    }
    i__1 = *q - 1;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if (phi[i__] < thresh)
        {
            phi[i__] = 0.;
        }
        else if (phi[i__] > 1.57079632679489662 - thresh)
        {
            phi[i__] = 1.57079632679489662;
        }
    }
    /* Initial deflation */
    imax = *q;
    while(imax > 1)
    {
        if (phi[imax - 1] != 0.)
        {
            break;
        }
        --imax;
    }
    imin = imax - 1;
    if (imin > 1)
    {
        while(phi[imin - 1] != 0.)
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
                if (phi[i__] != 0.)
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
        if (thetamax > 1.57079632679489662 - thresh)
        {
            /* Zero on diagonals of B11 and B22;
            induce deflation with a */
            /* zero shift */
            mu = 0.;
            nu = 1.;
        }
        else if (thetamin < thresh)
        {
            /* Zero on diagonals of B12 and B22;
            induce deflation with a */
            /* zero shift */
            mu = 1.;
            nu = 0.;
        }
        else
        {
            /* Compute shifts for B11 and B21 and use the lesser */
            dlas2_(&b11d[imax - 1], &b11e[imax - 1], &b11d[imax], &sigma11, & dummy);
            dlas2_(&b21d[imax - 1], &b21e[imax - 1], &b21d[imax], &sigma21, & dummy);
            if (sigma11 <= sigma21)
            {
                mu = sigma11;
                /* Computing 2nd power */
                d__1 = mu;
                nu = sqrt(1. - d__1 * d__1);
                if (mu < thresh)
                {
                    mu = 0.;
                    nu = 1.;
                }
            }
            else
            {
                nu = sigma21;
                /* Computing 2nd power */
                d__1 = nu;
                mu = sqrt(1.f - d__1 * d__1);
                if (nu < thresh)
                {
                    mu = 1.;
                    nu = 0.;
                }
            }
        }
        /* Rotate to produce bulges in B11 and B21 */
        if (mu <= nu)
        {
            dlartgs_(&b11d[imin], &b11e[imin], &mu, &work[iv1tcs + imin - 1], &work[iv1tsn + imin - 1]);
        }
        else
        {
            dlartgs_(&b21d[imin], &b21e[imin], &nu, &work[iv1tcs + imin - 1], &work[iv1tsn + imin - 1]);
        }
        temp = work[iv1tcs + imin - 1] * b11d[imin] + work[iv1tsn + imin - 1] * b11e[imin];
        b11e[imin] = work[iv1tcs + imin - 1] * b11e[imin] - work[iv1tsn + imin - 1] * b11d[imin];
        b11d[imin] = temp;
        b11bulge = work[iv1tsn + imin - 1] * b11d[imin + 1];
        b11d[imin + 1] = work[iv1tcs + imin - 1] * b11d[imin + 1];
        temp = work[iv1tcs + imin - 1] * b21d[imin] + work[iv1tsn + imin - 1] * b21e[imin];
        b21e[imin] = work[iv1tcs + imin - 1] * b21e[imin] - work[iv1tsn + imin - 1] * b21d[imin];
        b21d[imin] = temp;
        b21bulge = work[iv1tsn + imin - 1] * b21d[imin + 1];
        b21d[imin + 1] = work[iv1tcs + imin - 1] * b21d[imin + 1];
        /* Compute THETA(IMIN) */
        /* Computing 2nd power */
        d__1 = b21d[imin];
        /* Computing 2nd power */
        d__2 = b21bulge;
        /* Computing 2nd power */
        d__3 = b11d[imin];
        /* Computing 2nd power */
        d__4 = b11bulge;
        theta[imin] = atan2(sqrt(d__1 * d__1 + d__2 * d__2), sqrt(d__3 * d__3 + d__4 * d__4));
        /* Chase the bulges in B11(IMIN+1,IMIN) and B21(IMIN+1,IMIN) */
        /* Computing 2nd power */
        d__1 = b11d[imin];
        /* Computing 2nd power */
        d__2 = b11bulge;
        /* Computing 2nd power */
        d__3 = thresh;
        if (d__1 * d__1 + d__2 * d__2 > d__3 * d__3)
        {
            dlartgp_(&b11bulge, &b11d[imin], &work[iu1sn + imin - 1], &work[ iu1cs + imin - 1], &r__);
        }
        else if (mu <= nu)
        {
            dlartgs_(&b11e[imin], &b11d[imin + 1], &mu, &work[iu1cs + imin - 1], &work[iu1sn + imin - 1]);
        }
        else
        {
            dlartgs_(&b12d[imin], &b12e[imin], &nu, &work[iu1cs + imin - 1], & work[iu1sn + imin - 1]);
        }
        /* Computing 2nd power */
        d__1 = b21d[imin];
        /* Computing 2nd power */
        d__2 = b21bulge;
        /* Computing 2nd power */
        d__3 = thresh;
        if (d__1 * d__1 + d__2 * d__2 > d__3 * d__3)
        {
            dlartgp_(&b21bulge, &b21d[imin], &work[iu2sn + imin - 1], &work[ iu2cs + imin - 1], &r__);
        }
        else if (nu < mu)
        {
            dlartgs_(&b21e[imin], &b21d[imin + 1], &nu, &work[iu2cs + imin - 1], &work[iu2sn + imin - 1]);
        }
        else
        {
            dlartgs_(&b22d[imin], &b22e[imin], &mu, &work[iu2cs + imin - 1], & work[iu2sn + imin - 1]);
        }
        work[iu2cs + imin - 1] = -work[iu2cs + imin - 1];
        work[iu2sn + imin - 1] = -work[iu2sn + imin - 1];
        temp = work[iu1cs + imin - 1] * b11e[imin] + work[iu1sn + imin - 1] * b11d[imin + 1];
        b11d[imin + 1] = work[iu1cs + imin - 1] * b11d[imin + 1] - work[iu1sn + imin - 1] * b11e[imin];
        b11e[imin] = temp;
        if (imax > imin + 1)
        {
            b11bulge = work[iu1sn + imin - 1] * b11e[imin + 1];
            b11e[imin + 1] = work[iu1cs + imin - 1] * b11e[imin + 1];
        }
        temp = work[iu1cs + imin - 1] * b12d[imin] + work[iu1sn + imin - 1] * b12e[imin];
        b12e[imin] = work[iu1cs + imin - 1] * b12e[imin] - work[iu1sn + imin - 1] * b12d[imin];
        b12d[imin] = temp;
        b12bulge = work[iu1sn + imin - 1] * b12d[imin + 1];
        b12d[imin + 1] = work[iu1cs + imin - 1] * b12d[imin + 1];
        temp = work[iu2cs + imin - 1] * b21e[imin] + work[iu2sn + imin - 1] * b21d[imin + 1];
        b21d[imin + 1] = work[iu2cs + imin - 1] * b21d[imin + 1] - work[iu2sn + imin - 1] * b21e[imin];
        b21e[imin] = temp;
        if (imax > imin + 1)
        {
            b21bulge = work[iu2sn + imin - 1] * b21e[imin + 1];
            b21e[imin + 1] = work[iu2cs + imin - 1] * b21e[imin + 1];
        }
        temp = work[iu2cs + imin - 1] * b22d[imin] + work[iu2sn + imin - 1] * b22e[imin];
        b22e[imin] = work[iu2cs + imin - 1] * b22e[imin] - work[iu2sn + imin - 1] * b22d[imin];
        b22d[imin] = temp;
        b22bulge = work[iu2sn + imin - 1] * b22d[imin + 1];
        b22d[imin + 1] = work[iu2cs + imin - 1] * b22d[imin + 1];
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
            d__1 = x1;
            /* Computing 2nd power */
            d__2 = x2;
            /* Computing 2nd power */
            d__3 = y1;
            /* Computing 2nd power */
            d__4 = y2;
            phi[i__ - 1] = atan2(sqrt(d__1 * d__1 + d__2 * d__2), sqrt(d__3 * d__3 + d__4 * d__4));
            /* Determine if there are bulges to chase or if a new direct */
            /* summand has been reached */
            /* Computing 2nd power */
            d__1 = b11e[i__ - 1];
            /* Computing 2nd power */
            d__2 = b11bulge;
            /* Computing 2nd power */
            d__3 = thresh;
            restart11 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
            /* Computing 2nd power */
            d__1 = b21e[i__ - 1];
            /* Computing 2nd power */
            d__2 = b21bulge;
            /* Computing 2nd power */
            d__3 = thresh;
            restart21 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
            /* Computing 2nd power */
            d__1 = b12d[i__ - 1];
            /* Computing 2nd power */
            d__2 = b12bulge;
            /* Computing 2nd power */
            d__3 = thresh;
            restart12 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
            /* Computing 2nd power */
            d__1 = b22d[i__ - 1];
            /* Computing 2nd power */
            d__2 = b22bulge;
            /* Computing 2nd power */
            d__3 = thresh;
            restart22 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
            /* If possible, chase bulges from B11(I-1,I+1), B12(I-1,I), */
            /* B21(I-1,I+1), and B22(I-1,I). If necessary, restart bulge- */
            /* chasing by applying the original shift again. */
            if (! restart11 && ! restart21)
            {
                dlartgp_(&x2, &x1, &work[iv1tsn + i__ - 1], &work[iv1tcs + i__ - 1], &r__);
            }
            else if (! restart11 && restart21)
            {
                dlartgp_(&b11bulge, &b11e[i__ - 1], &work[iv1tsn + i__ - 1], & work[iv1tcs + i__ - 1], &r__);
            }
            else if (restart11 && ! restart21)
            {
                dlartgp_(&b21bulge, &b21e[i__ - 1], &work[iv1tsn + i__ - 1], & work[iv1tcs + i__ - 1], &r__);
            }
            else if (mu <= nu)
            {
                dlartgs_(&b11d[i__], &b11e[i__], &mu, &work[iv1tcs + i__ - 1], &work[iv1tsn + i__ - 1]);
            }
            else
            {
                dlartgs_(&b21d[i__], &b21e[i__], &nu, &work[iv1tcs + i__ - 1], &work[iv1tsn + i__ - 1]);
            }
            work[iv1tcs + i__ - 1] = -work[iv1tcs + i__ - 1];
            work[iv1tsn + i__ - 1] = -work[iv1tsn + i__ - 1];
            if (! restart12 && ! restart22)
            {
                dlartgp_(&y2, &y1, &work[iv2tsn + i__ - 2], &work[iv2tcs + i__ - 2], &r__);
            }
            else if (! restart12 && restart22)
            {
                dlartgp_(&b12bulge, &b12d[i__ - 1], &work[iv2tsn + i__ - 2], & work[iv2tcs + i__ - 2], &r__);
            }
            else if (restart12 && ! restart22)
            {
                dlartgp_(&b22bulge, &b22d[i__ - 1], &work[iv2tsn + i__ - 2], & work[iv2tcs + i__ - 2], &r__);
            }
            else if (nu < mu)
            {
                dlartgs_(&b12e[i__ - 1], &b12d[i__], &nu, &work[iv2tcs + i__ - 2], &work[iv2tsn + i__ - 2]);
            }
            else
            {
                dlartgs_(&b22e[i__ - 1], &b22d[i__], &mu, &work[iv2tcs + i__ - 2], &work[iv2tsn + i__ - 2]);
            }
            temp = work[iv1tcs + i__ - 1] * b11d[i__] + work[iv1tsn + i__ - 1] * b11e[i__];
            b11e[i__] = work[iv1tcs + i__ - 1] * b11e[i__] - work[iv1tsn + i__ - 1] * b11d[i__];
            b11d[i__] = temp;
            b11bulge = work[iv1tsn + i__ - 1] * b11d[i__ + 1];
            b11d[i__ + 1] = work[iv1tcs + i__ - 1] * b11d[i__ + 1];
            temp = work[iv1tcs + i__ - 1] * b21d[i__] + work[iv1tsn + i__ - 1] * b21e[i__];
            b21e[i__] = work[iv1tcs + i__ - 1] * b21e[i__] - work[iv1tsn + i__ - 1] * b21d[i__];
            b21d[i__] = temp;
            b21bulge = work[iv1tsn + i__ - 1] * b21d[i__ + 1];
            b21d[i__ + 1] = work[iv1tcs + i__ - 1] * b21d[i__ + 1];
            temp = work[iv2tcs + i__ - 2] * b12e[i__ - 1] + work[iv2tsn + i__ - 2] * b12d[i__];
            b12d[i__] = work[iv2tcs + i__ - 2] * b12d[i__] - work[iv2tsn + i__ - 2] * b12e[i__ - 1];
            b12e[i__ - 1] = temp;
            b12bulge = work[iv2tsn + i__ - 2] * b12e[i__];
            b12e[i__] = work[iv2tcs + i__ - 2] * b12e[i__];
            temp = work[iv2tcs + i__ - 2] * b22e[i__ - 1] + work[iv2tsn + i__ - 2] * b22d[i__];
            b22d[i__] = work[iv2tcs + i__ - 2] * b22d[i__] - work[iv2tsn + i__ - 2] * b22e[i__ - 1];
            b22e[i__ - 1] = temp;
            b22bulge = work[iv2tsn + i__ - 2] * b22e[i__];
            b22e[i__] = work[iv2tcs + i__ - 2] * b22e[i__];
            /* Compute THETA(I) */
            x1 = cos(phi[i__ - 1]) * b11d[i__] + sin(phi[i__ - 1]) * b12e[i__ - 1];
            x2 = cos(phi[i__ - 1]) * b11bulge + sin(phi[i__ - 1]) * b12bulge;
            y1 = cos(phi[i__ - 1]) * b21d[i__] + sin(phi[i__ - 1]) * b22e[i__ - 1];
            y2 = cos(phi[i__ - 1]) * b21bulge + sin(phi[i__ - 1]) * b22bulge;
            /* Computing 2nd power */
            d__1 = y1;
            /* Computing 2nd power */
            d__2 = y2;
            /* Computing 2nd power */
            d__3 = x1;
            /* Computing 2nd power */
            d__4 = x2;
            theta[i__] = atan2(sqrt(d__1 * d__1 + d__2 * d__2), sqrt(d__3 * d__3 + d__4 * d__4));
            /* Determine if there are bulges to chase or if a new direct */
            /* summand has been reached */
            /* Computing 2nd power */
            d__1 = b11d[i__];
            /* Computing 2nd power */
            d__2 = b11bulge;
            /* Computing 2nd power */
            d__3 = thresh;
            restart11 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
            /* Computing 2nd power */
            d__1 = b12e[i__ - 1];
            /* Computing 2nd power */
            d__2 = b12bulge;
            /* Computing 2nd power */
            d__3 = thresh;
            restart12 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
            /* Computing 2nd power */
            d__1 = b21d[i__];
            /* Computing 2nd power */
            d__2 = b21bulge;
            /* Computing 2nd power */
            d__3 = thresh;
            restart21 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
            /* Computing 2nd power */
            d__1 = b22e[i__ - 1];
            /* Computing 2nd power */
            d__2 = b22bulge;
            /* Computing 2nd power */
            d__3 = thresh;
            restart22 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
            /* If possible, chase bulges from B11(I+1,I), B12(I+1,I-1), */
            /* B21(I+1,I), and B22(I+1,I-1). If necessary, restart bulge- */
            /* chasing by applying the original shift again. */
            if (! restart11 && ! restart12)
            {
                dlartgp_(&x2, &x1, &work[iu1sn + i__ - 1], &work[iu1cs + i__ - 1], &r__);
            }
            else if (! restart11 && restart12)
            {
                dlartgp_(&b11bulge, &b11d[i__], &work[iu1sn + i__ - 1], &work[ iu1cs + i__ - 1], &r__);
            }
            else if (restart11 && ! restart12)
            {
                dlartgp_(&b12bulge, &b12e[i__ - 1], &work[iu1sn + i__ - 1], & work[iu1cs + i__ - 1], &r__);
            }
            else if (mu <= nu)
            {
                dlartgs_(&b11e[i__], &b11d[i__ + 1], &mu, &work[iu1cs + i__ - 1], &work[iu1sn + i__ - 1]);
            }
            else
            {
                dlartgs_(&b12d[i__], &b12e[i__], &nu, &work[iu1cs + i__ - 1], &work[iu1sn + i__ - 1]);
            }
            if (! restart21 && ! restart22)
            {
                dlartgp_(&y2, &y1, &work[iu2sn + i__ - 1], &work[iu2cs + i__ - 1], &r__);
            }
            else if (! restart21 && restart22)
            {
                dlartgp_(&b21bulge, &b21d[i__], &work[iu2sn + i__ - 1], &work[ iu2cs + i__ - 1], &r__);
            }
            else if (restart21 && ! restart22)
            {
                dlartgp_(&b22bulge, &b22e[i__ - 1], &work[iu2sn + i__ - 1], & work[iu2cs + i__ - 1], &r__);
            }
            else if (nu < mu)
            {
                dlartgs_(&b21e[i__], &b21e[i__ + 1], &nu, &work[iu2cs + i__ - 1], &work[iu2sn + i__ - 1]);
            }
            else
            {
                dlartgs_(&b22d[i__], &b22e[i__], &mu, &work[iu2cs + i__ - 1], &work[iu2sn + i__ - 1]);
            }
            work[iu2cs + i__ - 1] = -work[iu2cs + i__ - 1];
            work[iu2sn + i__ - 1] = -work[iu2sn + i__ - 1];
            temp = work[iu1cs + i__ - 1] * b11e[i__] + work[iu1sn + i__ - 1] * b11d[i__ + 1];
            b11d[i__ + 1] = work[iu1cs + i__ - 1] * b11d[i__ + 1] - work[ iu1sn + i__ - 1] * b11e[i__];
            b11e[i__] = temp;
            if (i__ < imax - 1)
            {
                b11bulge = work[iu1sn + i__ - 1] * b11e[i__ + 1];
                b11e[i__ + 1] = work[iu1cs + i__ - 1] * b11e[i__ + 1];
            }
            temp = work[iu2cs + i__ - 1] * b21e[i__] + work[iu2sn + i__ - 1] * b21d[i__ + 1];
            b21d[i__ + 1] = work[iu2cs + i__ - 1] * b21d[i__ + 1] - work[ iu2sn + i__ - 1] * b21e[i__];
            b21e[i__] = temp;
            if (i__ < imax - 1)
            {
                b21bulge = work[iu2sn + i__ - 1] * b21e[i__ + 1];
                b21e[i__ + 1] = work[iu2cs + i__ - 1] * b21e[i__ + 1];
            }
            temp = work[iu1cs + i__ - 1] * b12d[i__] + work[iu1sn + i__ - 1] * b12e[i__];
            b12e[i__] = work[iu1cs + i__ - 1] * b12e[i__] - work[iu1sn + i__ - 1] * b12d[i__];
            b12d[i__] = temp;
            b12bulge = work[iu1sn + i__ - 1] * b12d[i__ + 1];
            b12d[i__ + 1] = work[iu1cs + i__ - 1] * b12d[i__ + 1];
            temp = work[iu2cs + i__ - 1] * b22d[i__] + work[iu2sn + i__ - 1] * b22e[i__];
            b22e[i__] = work[iu2cs + i__ - 1] * b22e[i__] - work[iu2sn + i__ - 1] * b22d[i__];
            b22d[i__] = temp;
            b22bulge = work[iu2sn + i__ - 1] * b22d[i__ + 1];
            b22d[i__ + 1] = work[iu2cs + i__ - 1] * b22d[i__ + 1];
        }
        /* Compute PHI(IMAX-1) */
        x1 = sin(theta[imax - 1]) * b11e[imax - 1] + cos(theta[imax - 1]) * b21e[imax - 1];
        y1 = sin(theta[imax - 1]) * b12d[imax - 1] + cos(theta[imax - 1]) * b22d[imax - 1];
        y2 = sin(theta[imax - 1]) * b12bulge + cos(theta[imax - 1]) * b22bulge;
        /* Computing 2nd power */
        d__1 = y1;
        /* Computing 2nd power */
        d__2 = y2;
        phi[imax - 1] = atan2((f2c_abs(x1)), sqrt(d__1 * d__1 + d__2 * d__2));
        /* Chase bulges from B12(IMAX-1,IMAX) and B22(IMAX-1,IMAX) */
        /* Computing 2nd power */
        d__1 = b12d[imax - 1];
        /* Computing 2nd power */
        d__2 = b12bulge;
        /* Computing 2nd power */
        d__3 = thresh;
        restart12 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
        /* Computing 2nd power */
        d__1 = b22d[imax - 1];
        /* Computing 2nd power */
        d__2 = b22bulge;
        /* Computing 2nd power */
        d__3 = thresh;
        restart22 = d__1 * d__1 + d__2 * d__2 <= d__3 * d__3;
        if (! restart12 && ! restart22)
        {
            dlartgp_(&y2, &y1, &work[iv2tsn + imax - 2], &work[iv2tcs + imax - 2], &r__);
        }
        else if (! restart12 && restart22)
        {
            dlartgp_(&b12bulge, &b12d[imax - 1], &work[iv2tsn + imax - 2], & work[iv2tcs + imax - 2], &r__);
        }
        else if (restart12 && ! restart22)
        {
            dlartgp_(&b22bulge, &b22d[imax - 1], &work[iv2tsn + imax - 2], & work[iv2tcs + imax - 2], &r__);
        }
        else if (nu < mu)
        {
            dlartgs_(&b12e[imax - 1], &b12d[imax], &nu, &work[iv2tcs + imax - 2], &work[iv2tsn + imax - 2]);
        }
        else
        {
            dlartgs_(&b22e[imax - 1], &b22d[imax], &mu, &work[iv2tcs + imax - 2], &work[iv2tsn + imax - 2]);
        }
        temp = work[iv2tcs + imax - 2] * b12e[imax - 1] + work[iv2tsn + imax - 2] * b12d[imax];
        b12d[imax] = work[iv2tcs + imax - 2] * b12d[imax] - work[iv2tsn + imax - 2] * b12e[imax - 1];
        b12e[imax - 1] = temp;
        temp = work[iv2tcs + imax - 2] * b22e[imax - 1] + work[iv2tsn + imax - 2] * b22d[imax];
        b22d[imax] = work[iv2tcs + imax - 2] * b22d[imax] - work[iv2tsn + imax - 2] * b22e[imax - 1];
        b22e[imax - 1] = temp;
        /* Update singular vectors */
        if (wantu1)
        {
            if (colmajor)
            {
                i__1 = imax - imin + 1;
                dlasr_("R", "V", "F", p, &i__1, &work[iu1cs + imin - 1], & work[iu1sn + imin - 1], &u1[imin * u1_dim1 + 1], ldu1);
            }
            else
            {
                i__1 = imax - imin + 1;
                dlasr_("L", "V", "F", &i__1, p, &work[iu1cs + imin - 1], & work[iu1sn + imin - 1], &u1[imin + u1_dim1], ldu1);
            }
        }
        if (wantu2)
        {
            if (colmajor)
            {
                i__1 = *m - *p;
                i__2 = imax - imin + 1;
                dlasr_("R", "V", "F", &i__1, &i__2, &work[iu2cs + imin - 1], & work[iu2sn + imin - 1], &u2[imin * u2_dim1 + 1], ldu2);
            }
            else
            {
                i__1 = imax - imin + 1;
                i__2 = *m - *p;
                dlasr_("L", "V", "F", &i__1, &i__2, &work[iu2cs + imin - 1], & work[iu2sn + imin - 1], &u2[imin + u2_dim1], ldu2);
            }
        }
        if (wantv1t)
        {
            if (colmajor)
            {
                i__1 = imax - imin + 1;
                dlasr_("L", "V", "F", &i__1, q, &work[iv1tcs + imin - 1], & work[iv1tsn + imin - 1], &v1t[imin + v1t_dim1], ldv1t);
            }
            else
            {
                i__1 = imax - imin + 1;
                dlasr_("R", "V", "F", q, &i__1, &work[iv1tcs + imin - 1], & work[iv1tsn + imin - 1], &v1t[imin * v1t_dim1 + 1], ldv1t);
            }
        }
        if (wantv2t)
        {
            if (colmajor)
            {
                i__1 = imax - imin + 1;
                i__2 = *m - *q;
                dlasr_("L", "V", "F", &i__1, &i__2, &work[iv2tcs + imin - 1], &work[iv2tsn + imin - 1], &v2t[imin + v2t_dim1], ldv2t);
            }
            else
            {
                i__1 = *m - *q;
                i__2 = imax - imin + 1;
                dlasr_("R", "V", "F", &i__1, &i__2, &work[iv2tcs + imin - 1], &work[iv2tsn + imin - 1], &v2t[imin * v2t_dim1 + 1], ldv2t);
            }
        }
        /* Fix signs on B11(IMAX-1,IMAX) and B21(IMAX-1,IMAX) */
        if (b11e[imax - 1] + b21e[imax - 1] > 0.)
        {
            b11d[imax] = -b11d[imax];
            b21d[imax] = -b21d[imax];
            if (wantv1t)
            {
                if (colmajor)
                {
                    dscal_(q, &c_b35, &v1t[imax + v1t_dim1], ldv1t);
                }
                else
                {
                    dscal_(q, &c_b35, &v1t[imax * v1t_dim1 + 1], &c__1);
                }
            }
        }
        /* Compute THETA(IMAX) */
        x1 = cos(phi[imax - 1]) * b11d[imax] + sin(phi[imax - 1]) * b12e[imax - 1];
        y1 = cos(phi[imax - 1]) * b21d[imax] + sin(phi[imax - 1]) * b22e[imax - 1];
        theta[imax] = atan2((f2c_abs(y1)), (f2c_abs(x1)));
        /* Fix signs on B11(IMAX,IMAX), B12(IMAX,IMAX-1), B21(IMAX,IMAX), */
        /* and B22(IMAX,IMAX-1) */
        if (b11d[imax] + b12e[imax - 1] < 0.)
        {
            b12d[imax] = -b12d[imax];
            if (wantu1)
            {
                if (colmajor)
                {
                    dscal_(p, &c_b35, &u1[imax * u1_dim1 + 1], &c__1);
                }
                else
                {
                    dscal_(p, &c_b35, &u1[imax + u1_dim1], ldu1);
                }
            }
        }
        if (b21d[imax] + b22e[imax - 1] > 0.)
        {
            b22d[imax] = -b22d[imax];
            if (wantu2)
            {
                if (colmajor)
                {
                    i__1 = *m - *p;
                    dscal_(&i__1, &c_b35, &u2[imax * u2_dim1 + 1], &c__1);
                }
                else
                {
                    i__1 = *m - *p;
                    dscal_(&i__1, &c_b35, &u2[imax + u2_dim1], ldu2);
                }
            }
        }
        /* Fix signs on B12(IMAX,IMAX) and B22(IMAX,IMAX) */
        if (b12d[imax] + b22d[imax] < 0.)
        {
            if (wantv2t)
            {
                if (colmajor)
                {
                    i__1 = *m - *q;
                    dscal_(&i__1, &c_b35, &v2t[imax + v2t_dim1], ldv2t);
                }
                else
                {
                    i__1 = *m - *q;
                    dscal_(&i__1, &c_b35, &v2t[imax * v2t_dim1 + 1], &c__1);
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
                theta[i__] = 0.;
            }
            else if (theta[i__] > 1.57079632679489662 - thresh)
            {
                theta[i__] = 1.57079632679489662;
            }
        }
        i__1 = imax - 1;
        for (i__ = imin;
                i__ <= i__1;
                ++i__)
        {
            if (phi[i__] < thresh)
            {
                phi[i__] = 0.;
            }
            else if (phi[i__] > 1.57079632679489662 - thresh)
            {
                phi[i__] = 1.57079632679489662;
            }
        }
        /* Deflate */
        if (imax > 1)
        {
            while(phi[imax - 1] == 0.)
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
            while(phi[imin - 1] != 0.)
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
                    dswap_(p, &u1[i__ * u1_dim1 + 1], &c__1, &u1[mini * u1_dim1 + 1], &c__1);
                }
                if (wantu2)
                {
                    i__2 = *m - *p;
                    dswap_(&i__2, &u2[i__ * u2_dim1 + 1], &c__1, &u2[mini * u2_dim1 + 1], &c__1);
                }
                if (wantv1t)
                {
                    dswap_(q, &v1t[i__ + v1t_dim1], ldv1t, &v1t[mini + v1t_dim1], ldv1t);
                }
                if (wantv2t)
                {
                    i__2 = *m - *q;
                    dswap_(&i__2, &v2t[i__ + v2t_dim1], ldv2t, &v2t[mini + v2t_dim1], ldv2t);
                }
            }
            else
            {
                if (wantu1)
                {
                    dswap_(p, &u1[i__ + u1_dim1], ldu1, &u1[mini + u1_dim1], ldu1);
                }
                if (wantu2)
                {
                    i__2 = *m - *p;
                    dswap_(&i__2, &u2[i__ + u2_dim1], ldu2, &u2[mini + u2_dim1], ldu2);
                }
                if (wantv1t)
                {
                    dswap_(q, &v1t[i__ * v1t_dim1 + 1], &c__1, &v1t[mini * v1t_dim1 + 1], &c__1);
                }
                if (wantv2t)
                {
                    i__2 = *m - *q;
                    dswap_(&i__2, &v2t[i__ * v2t_dim1 + 1], &c__1, &v2t[mini * v2t_dim1 + 1], &c__1);
                }
            }
        }
    }
    return 0;
    /* End of DBBCSD */
}
/* dbbcsd_ */
