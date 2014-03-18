/* ../netlib/zuncsd.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
static logical c_false = FALSE_;
/* > \brief \b ZUNCSD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZUNCSD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zuncsd. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zuncsd. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zuncsd. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZUNCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, */
/* SIGNS, M, P, Q, X11, LDX11, X12, */
/* LDX12, X21, LDX21, X22, LDX22, THETA, */
/* U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, */
/* LDV2T, WORK, LWORK, RWORK, LRWORK, */
/* IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS */
/* INTEGER INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12, */
/* $ LDX21, LDX22, LRWORK, LWORK, M, P, Q */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION THETA( * ) */
/* DOUBLE PRECISION RWORK( * ) */
/* COMPLEX*16 U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), */
/* $ V2T( LDV2T, * ), WORK( * ), X11( LDX11, * ), */
/* $ X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, */
/* $ * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNCSD computes the CS decomposition of an M-by-M partitioned */
/* > unitary matrix X: */
/* > */
/* > [ I 0 0 | 0 0 0 ] */
/* > [ 0 C 0 | 0 -S 0 ] */
/* > [ X11 | X12 ] [ U1 | ] [ 0 0 0 | 0 0 -I ] [ V1 | ]**H */
/* > X = [-----------] = [---------] [---------------------] [---------] . */
/* > [ X21 | X22 ] [ | U2 ] [ 0 0 0 | I 0 0 ] [ | V2 ] */
/* > [ 0 S 0 | 0 C 0 ] */
/* > [ 0 0 I | 0 0 0 ] */
/* > */
/* > X11 is P-by-Q. The unitary matrices U1, U2, V1, and V2 are P-by-P, */
/* > (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are */
/* > R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in */
/* > which R = MIN(P,M-P,Q,M-Q). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBU1 */
/* > \verbatim */
/* > JOBU1 is CHARACTER */
/* > = 'Y': U1 is computed;
*/
/* > otherwise: U1 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU2 */
/* > \verbatim */
/* > JOBU2 is CHARACTER */
/* > = 'Y': U2 is computed;
*/
/* > otherwise: U2 is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV1T */
/* > \verbatim */
/* > JOBV1T is CHARACTER */
/* > = 'Y': V1T is computed;
*/
/* > otherwise: V1T is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV2T */
/* > \verbatim */
/* > JOBV2T is CHARACTER */
/* > = 'Y': V2T is computed;
*/
/* > otherwise: V2T is not computed. */
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
/* > \param[in] SIGNS */
/* > \verbatim */
/* > SIGNS is CHARACTER */
/* > = 'O': The lower-left block is made nonpositive (the */
/* > "other" convention);
*/
/* > otherwise: The upper-right block is made nonpositive (the */
/* > "default" convention). */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows and columns in X. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* > P is INTEGER */
/* > The number of rows in X11 and X12. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* > Q is INTEGER */
/* > The number of columns in X11 and X21. 0 <= Q <= M. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* > X11 is COMPLEX*16 array, dimension (LDX11,Q) */
/* > On entry, part of the unitary matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* > LDX11 is INTEGER */
/* > The leading dimension of X11. LDX11 >= MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X12 */
/* > \verbatim */
/* > X12 is COMPLEX*16 array, dimension (LDX12,M-Q) */
/* > On entry, part of the unitary matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX12 */
/* > \verbatim */
/* > LDX12 is INTEGER */
/* > The leading dimension of X12. LDX12 >= MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* > X21 is COMPLEX*16 array, dimension (LDX21,Q) */
/* > On entry, part of the unitary matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* > LDX21 is INTEGER */
/* > The leading dimension of X11. LDX21 >= MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X22 */
/* > \verbatim */
/* > X22 is COMPLEX*16 array, dimension (LDX22,M-Q) */
/* > On entry, part of the unitary matrix whose CSD is desired. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX22 */
/* > \verbatim */
/* > LDX22 is INTEGER */
/* > The leading dimension of X11. LDX22 >= MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* > THETA is DOUBLE PRECISION array, dimension (R), in which R = */
/* > MIN(P,M-P,Q,M-Q). */
/* > C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and */
/* > S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ). */
/* > \endverbatim */
/* > */
/* > \param[out] U1 */
/* > \verbatim */
/* > U1 is COMPLEX*16 array, dimension (P) */
/* > If JOBU1 = 'Y', U1 contains the P-by-P unitary matrix U1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU1 */
/* > \verbatim */
/* > LDU1 is INTEGER */
/* > The leading dimension of U1. If JOBU1 = 'Y', LDU1 >= */
/* > MAX(1,P). */
/* > \endverbatim */
/* > */
/* > \param[out] U2 */
/* > \verbatim */
/* > U2 is COMPLEX*16 array, dimension (M-P) */
/* > If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) unitary */
/* > matrix U2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU2 */
/* > \verbatim */
/* > LDU2 is INTEGER */
/* > The leading dimension of U2. If JOBU2 = 'Y', LDU2 >= */
/* > MAX(1,M-P). */
/* > \endverbatim */
/* > */
/* > \param[out] V1T */
/* > \verbatim */
/* > V1T is COMPLEX*16 array, dimension (Q) */
/* > If JOBV1T = 'Y', V1T contains the Q-by-Q matrix unitary */
/* > matrix V1**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV1T */
/* > \verbatim */
/* > LDV1T is INTEGER */
/* > The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >= */
/* > MAX(1,Q). */
/* > \endverbatim */
/* > */
/* > \param[out] V2T */
/* > \verbatim */
/* > V2T is COMPLEX*16 array, dimension (M-Q) */
/* > If JOBV2T = 'Y', V2T contains the (M-Q)-by-(M-Q) unitary */
/* > matrix V2**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV2T */
/* > \verbatim */
/* > LDV2T is INTEGER */
/* > The leading dimension of V2T. If JOBV2T = 'Y', LDV2T >= */
/* > MAX(1,M-Q). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the work array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension MAX(1,LRWORK) */
/* > On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > If INFO > 0 on exit, RWORK(2:R) contains the values PHI(1), */
/* > ..., PHI(R-1) that, together with THETA(1), ..., THETA(R), */
/* > define the matrix in intermediate bidiagonal-block form */
/* > remaining after nonconvergence. INFO specifies the number */
/* > of nonzero PHI's. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* > LRWORK is INTEGER */
/* > The dimension of the array RWORK. */
/* > */
/* > If LRWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the RWORK array, returns */
/* > this value as the first entry of the work array, and no error */
/* > message related to LRWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (M-MIN(P,M-P,Q,M-Q)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: ZBBCSD did not converge. See the description of RWORK */
/* > above for details. */
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
/* > \ingroup complex16OTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int zuncsd_(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, char *signs, integer *m, integer *p, integer *q, doublecomplex *x11, integer *ldx11, doublecomplex *x12, integer * ldx12, doublecomplex *x21, integer *ldx21, doublecomplex *x22, integer *ldx22, doublereal *theta, doublecomplex *u1, integer *ldu1, doublecomplex *u2, integer *ldu2, doublecomplex *v1t, integer *ldv1t, doublecomplex *v2t, integer *ldv2t, doublecomplex *work, integer * lwork, doublereal *rwork, integer *lrwork, integer *iwork, integer * info)
{
    /* System generated locals */
    integer u1_dim1, u1_offset, u2_dim1, u2_offset, v1t_dim1, v1t_offset, v2t_dim1, v2t_offset, x11_dim1, x11_offset, x12_dim1, x12_offset, x21_dim1, x21_offset, x22_dim1, x22_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    /* Local variables */
    logical colmajor;
    integer lworkmin, lworkopt, i__, j, childinfo, p1, q1, lrworkmin, lrworkopt, lbbcsdwork, lorbdbwork, lorglqwork, lorgqrwork, ib11d, ib11e, ib12d, ib12e, ib21d, ib21e, ib22d, ib22e, iphi;
    logical defaultsigns;
    extern logical lsame_(char *, char *);
    integer lbbcsdworkmin, itaup1, itaup2, itauq1, itauq2, lorbdbworkmin, lbbcsdworkopt;
    logical wantu1, wantu2;
    integer ibbcsd, lorbdbworkopt, iorbdb, lorglqworkmin;
    extern /* Subroutine */
    int zbbcsd_(char *, char *, char *, char *, char * , integer *, integer *, integer *, doublereal *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, integer *);
    integer lorgqrworkmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    integer lorglqworkopt;
    extern /* Subroutine */
    int zunbdb_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, doublereal *, doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *, integer *, integer *);
    integer lorgqrworkopt, iorglq;
    extern /* Subroutine */
    int zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    integer iorgqr;
    extern /* Subroutine */
    int zlapmr_(logical *, integer *, integer *, doublecomplex *, integer *, integer *);
    char signst[1];
    extern /* Subroutine */
    int zlapmt_(logical *, integer *, integer *, doublecomplex *, integer *, integer *);
    char transt[1];
    logical lquery;
    extern /* Subroutine */
    int zunglq_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, integer *), zungqr_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, integer *);
    logical wantv1t, wantv2t, lrquery;
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
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions */
    /* .. */
    /* .. Executable Statements .. */
    /* Test input arguments */
    /* Parameter adjustments */
    x11_dim1 = *ldx11;
    x11_offset = 1 + x11_dim1;
    x11 -= x11_offset;
    x12_dim1 = *ldx12;
    x12_offset = 1 + x12_dim1;
    x12 -= x12_offset;
    x21_dim1 = *ldx21;
    x21_offset = 1 + x21_dim1;
    x21 -= x21_offset;
    x22_dim1 = *ldx22;
    x22_offset = 1 + x22_dim1;
    x22 -= x22_offset;
    --theta;
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
    --work;
    --rwork;
    --iwork;
    /* Function Body */
    *info = 0;
    wantu1 = lsame_(jobu1, "Y");
    wantu2 = lsame_(jobu2, "Y");
    wantv1t = lsame_(jobv1t, "Y");
    wantv2t = lsame_(jobv2t, "Y");
    colmajor = ! lsame_(trans, "T");
    defaultsigns = ! lsame_(signs, "O");
    lquery = *lwork == -1;
    lrquery = *lrwork == -1;
    if (*m < 0)
    {
        *info = -7;
    }
    else if (*p < 0 || *p > *m)
    {
        *info = -8;
    }
    else if (*q < 0 || *q > *m)
    {
        *info = -9;
    }
    else if (colmajor && *ldx11 < max(1,*p))
    {
        *info = -11;
    }
    else if (! colmajor && *ldx11 < max(1,*q))
    {
        *info = -11;
    }
    else if (colmajor && *ldx12 < max(1,*p))
    {
        *info = -13;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        if (! colmajor && *ldx12 < max(i__1,i__2))
        {
            *info = -13;
        }
        else /* if(complicated condition) */
        {
            /* Computing MAX */
            i__1 = 1;
            i__2 = *m - *p; // , expr subst
            if (colmajor && *ldx21 < max(i__1,i__2))
            {
                *info = -15;
            }
            else if (! colmajor && *ldx21 < max(1,*q))
            {
                *info = -15;
            }
            else /* if(complicated condition) */
            {
                /* Computing MAX */
                i__1 = 1;
                i__2 = *m - *p; // , expr subst
                if (colmajor && *ldx22 < max(i__1,i__2))
                {
                    *info = -17;
                }
                else /* if(complicated condition) */
                {
                    /* Computing MAX */
                    i__1 = 1;
                    i__2 = *m - *q; // , expr subst
                    if (! colmajor && *ldx22 < max(i__1,i__2))
                    {
                        *info = -17;
                    }
                    else if (wantu1 && *ldu1 < *p)
                    {
                        *info = -20;
                    }
                    else if (wantu2 && *ldu2 < *m - *p)
                    {
                        *info = -22;
                    }
                    else if (wantv1t && *ldv1t < *q)
                    {
                        *info = -24;
                    }
                    else if (wantv2t && *ldv2t < *m - *q)
                    {
                        *info = -26;
                    }
                }
            }
        }
    }
    /* Work with transpose if convenient */
    /* Computing MIN */
    i__1 = *p;
    i__2 = *m - *p; // , expr subst
    /* Computing MIN */
    i__3 = *q;
    i__4 = *m - *q; // , expr subst
    if (*info == 0 && min(i__1,i__2) < min(i__3,i__4))
    {
        if (colmajor)
        {
            *(unsigned char *)transt = 'T';
        }
        else
        {
            *(unsigned char *)transt = 'N';
        }
        if (defaultsigns)
        {
            *(unsigned char *)signst = 'O';
        }
        else
        {
            *(unsigned char *)signst = 'D';
        }
        zuncsd_(jobv1t, jobv2t, jobu1, jobu2, transt, signst, m, q, p, &x11[ x11_offset], ldx11, &x21[x21_offset], ldx21, &x12[x12_offset], ldx12, &x22[x22_offset], ldx22, &theta[1], &v1t[v1t_offset], ldv1t, &v2t[v2t_offset], ldv2t, &u1[u1_offset], ldu1, &u2[ u2_offset], ldu2, &work[1], lwork, &rwork[1], lrwork, &iwork[ 1], info);
        return 0;
    }
    /* Work with permutation [ 0 I;
    I 0 ] * X * [ 0 I;
    I 0 ] if */
    /* convenient */
    if (*info == 0 && *m - *q < *q)
    {
        if (defaultsigns)
        {
            *(unsigned char *)signst = 'O';
        }
        else
        {
            *(unsigned char *)signst = 'D';
        }
        i__1 = *m - *p;
        i__2 = *m - *q;
        zuncsd_(jobu2, jobu1, jobv2t, jobv1t, trans, signst, m, &i__1, &i__2, &x22[x22_offset], ldx22, &x21[x21_offset], ldx21, &x12[ x12_offset], ldx12, &x11[x11_offset], ldx11, &theta[1], &u2[ u2_offset], ldu2, &u1[u1_offset], ldu1, &v2t[v2t_offset], ldv2t, &v1t[v1t_offset], ldv1t, &work[1], lwork, &rwork[1], lrwork, &iwork[1], info);
        return 0;
    }
    /* Compute workspace */
    if (*info == 0)
    {
        /* Real workspace */
        iphi = 2;
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        ib11d = iphi + max(i__1,i__2);
        ib11e = ib11d + max(1,*q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        ib12d = ib11e + max(i__1,i__2);
        ib12e = ib12d + max(1,*q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        ib21d = ib12e + max(i__1,i__2);
        ib21e = ib21d + max(1,*q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        ib22d = ib21e + max(i__1,i__2);
        ib22e = ib22d + max(1,*q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *q - 1; // , expr subst
        ibbcsd = ib22e + max(i__1,i__2);
        zbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, &theta[1], & theta[1], &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[ v1t_offset], ldv1t, &v2t[v2t_offset], ldv2t, &theta[1], & theta[1], &theta[1], &theta[1], &theta[1], &theta[1], &theta[ 1], &theta[1], &rwork[1], &c_n1, &childinfo);
        lbbcsdworkopt = (integer) rwork[1];
        lbbcsdworkmin = lbbcsdworkopt;
        lrworkopt = ibbcsd + lbbcsdworkopt - 1;
        lrworkmin = ibbcsd + lbbcsdworkmin - 1;
        rwork[1] = (doublereal) lrworkopt;
        /* Complex workspace */
        itaup1 = 2;
        itaup2 = itaup1 + max(1,*p);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *p; // , expr subst
        itauq1 = itaup2 + max(i__1,i__2);
        itauq2 = itauq1 + max(1,*q);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        iorgqr = itauq2 + max(i__1,i__2);
        i__1 = *m - *q;
        i__2 = *m - *q;
        i__3 = *m - *q;
        /* Computing MAX */
        i__5 = 1;
        i__6 = *m - *q; // , expr subst
        i__4 = max(i__5,i__6);
        zungqr_(&i__1, &i__2, &i__3, &u1[u1_offset], &i__4, &u1[u1_offset], & work[1], &c_n1, &childinfo);
        lorgqrworkopt = (integer) work[1].r;
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        lorgqrworkmin = max(i__1,i__2);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        iorglq = itauq2 + max(i__1,i__2);
        i__1 = *m - *q;
        i__2 = *m - *q;
        i__3 = *m - *q;
        /* Computing MAX */
        i__5 = 1;
        i__6 = *m - *q; // , expr subst
        i__4 = max(i__5,i__6);
        zunglq_(&i__1, &i__2, &i__3, &u1[u1_offset], &i__4, &u1[u1_offset], & work[1], &c_n1, &childinfo);
        lorglqworkopt = (integer) work[1].r;
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        lorglqworkmin = max(i__1,i__2);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *m - *q; // , expr subst
        iorbdb = itauq2 + max(i__1,i__2);
        zunbdb_(trans, signs, m, p, q, &x11[x11_offset], ldx11, &x12[ x12_offset], ldx12, &x21[x21_offset], ldx21, &x22[x22_offset], ldx22, &theta[1], &theta[1], &u1[u1_offset], &u2[u2_offset], &v1t[v1t_offset], &v2t[v2t_offset], &work[1], &c_n1, & childinfo);
        lorbdbworkopt = (integer) work[1].r;
        lorbdbworkmin = lorbdbworkopt;
        /* Computing MAX */
        i__1 = iorgqr + lorgqrworkopt, i__2 = iorglq + lorglqworkopt;
        i__1 = max(i__1,i__2);
        i__2 = iorbdb + lorbdbworkopt; // ; expr subst
        lworkopt = max(i__1,i__2) - 1;
        /* Computing MAX */
        i__1 = iorgqr + lorgqrworkmin, i__2 = iorglq + lorglqworkmin;
        i__1 = max(i__1,i__2);
        i__2 = iorbdb + lorbdbworkmin; // ; expr subst
        lworkmin = max(i__1,i__2) - 1;
        i__1 = max(lworkopt,lworkmin);
        work[1].r = (doublereal) i__1;
        work[1].i = 0.; // , expr subst
        if (*lwork < lworkmin && ! (lquery || lrquery))
        {
            *info = -22;
        }
        else if (*lrwork < lrworkmin && ! (lquery || lrquery))
        {
            *info = -24;
        }
        else
        {
            lorgqrwork = *lwork - iorgqr + 1;
            lorglqwork = *lwork - iorglq + 1;
            lorbdbwork = *lwork - iorbdb + 1;
            lbbcsdwork = *lrwork - ibbcsd + 1;
        }
    }
    /* Abort if any illegal arguments */
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZUNCSD", &i__1);
        return 0;
    }
    else if (lquery || lrquery)
    {
        return 0;
    }
    /* Transform to bidiagonal block form */
    zunbdb_(trans, signs, m, p, q, &x11[x11_offset], ldx11, &x12[x12_offset], ldx12, &x21[x21_offset], ldx21, &x22[x22_offset], ldx22, &theta[1] , &rwork[iphi], &work[itaup1], &work[itaup2], &work[itauq1], & work[itauq2], &work[iorbdb], &lorbdbwork, &childinfo);
    /* Accumulate Householder reflectors */
    if (colmajor)
    {
        if (wantu1 && *p > 0)
        {
            zlacpy_("L", p, q, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1);
            zungqr_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[ iorgqr], &lorgqrwork, info);
        }
        if (wantu2 && *m - *p > 0)
        {
            i__1 = *m - *p;
            zlacpy_("L", &i__1, q, &x21[x21_offset], ldx21, &u2[u2_offset], ldu2);
            i__1 = *m - *p;
            i__2 = *m - *p;
            zungqr_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], & work[iorgqr], &lorgqrwork, info);
        }
        if (wantv1t && *q > 0)
        {
            i__1 = *q - 1;
            i__2 = *q - 1;
            zlacpy_("U", &i__1, &i__2, &x11[(x11_dim1 << 1) + 1], ldx11, &v1t[ (v1t_dim1 << 1) + 2], ldv1t);
            i__1 = v1t_dim1 + 1;
            v1t[i__1].r = 1.;
            v1t[i__1].i = 0.; // , expr subst
            i__1 = *q;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                i__2 = j * v1t_dim1 + 1;
                v1t[i__2].r = 0.;
                v1t[i__2].i = 0.; // , expr subst
                i__2 = j + v1t_dim1;
                v1t[i__2].r = 0.;
                v1t[i__2].i = 0.; // , expr subst
            }
            i__1 = *q - 1;
            i__2 = *q - 1;
            i__3 = *q - 1;
            zunglq_(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, & work[itauq1], &work[iorglq], &lorglqwork, info);
        }
        if (wantv2t && *m - *q > 0)
        {
            i__1 = *m - *q;
            zlacpy_("U", p, &i__1, &x12[x12_offset], ldx12, &v2t[v2t_offset], ldv2t);
            if (*m - *p > *q)
            {
                i__1 = *m - *p - *q;
                i__2 = *m - *p - *q;
                zlacpy_("U", &i__1, &i__2, &x22[*q + 1 + (*p + 1) * x22_dim1], ldx22, &v2t[*p + 1 + (*p + 1) * v2t_dim1], ldv2t);
            }
            if (*m > *q)
            {
                i__1 = *m - *q;
                i__2 = *m - *q;
                i__3 = *m - *q;
                zunglq_(&i__1, &i__2, &i__3, &v2t[v2t_offset], ldv2t, &work[ itauq2], &work[iorglq], &lorglqwork, info);
            }
        }
    }
    else
    {
        if (wantu1 && *p > 0)
        {
            zlacpy_("U", q, p, &x11[x11_offset], ldx11, &u1[u1_offset], ldu1);
            zunglq_(p, p, q, &u1[u1_offset], ldu1, &work[itaup1], &work[ iorglq], &lorglqwork, info);
        }
        if (wantu2 && *m - *p > 0)
        {
            i__1 = *m - *p;
            zlacpy_("U", q, &i__1, &x21[x21_offset], ldx21, &u2[u2_offset], ldu2);
            i__1 = *m - *p;
            i__2 = *m - *p;
            zunglq_(&i__1, &i__2, q, &u2[u2_offset], ldu2, &work[itaup2], & work[iorglq], &lorglqwork, info);
        }
        if (wantv1t && *q > 0)
        {
            i__1 = *q - 1;
            i__2 = *q - 1;
            zlacpy_("L", &i__1, &i__2, &x11[x11_dim1 + 2], ldx11, &v1t[( v1t_dim1 << 1) + 2], ldv1t);
            i__1 = v1t_dim1 + 1;
            v1t[i__1].r = 1.;
            v1t[i__1].i = 0.; // , expr subst
            i__1 = *q;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                i__2 = j * v1t_dim1 + 1;
                v1t[i__2].r = 0.;
                v1t[i__2].i = 0.; // , expr subst
                i__2 = j + v1t_dim1;
                v1t[i__2].r = 0.;
                v1t[i__2].i = 0.; // , expr subst
            }
            i__1 = *q - 1;
            i__2 = *q - 1;
            i__3 = *q - 1;
            zungqr_(&i__1, &i__2, &i__3, &v1t[(v1t_dim1 << 1) + 2], ldv1t, & work[itauq1], &work[iorgqr], &lorgqrwork, info);
        }
        if (wantv2t && *m - *q > 0)
        {
            /* Computing MIN */
            i__1 = *p + 1;
            p1 = min(i__1,*m);
            /* Computing MIN */
            i__1 = *q + 1;
            q1 = min(i__1,*m);
            i__1 = *m - *q;
            zlacpy_("L", &i__1, p, &x12[x12_offset], ldx12, &v2t[v2t_offset], ldv2t);
            if (*m > *p + *q)
            {
                i__1 = *m - *p - *q;
                i__2 = *m - *p - *q;
                zlacpy_("L", &i__1, &i__2, &x22[p1 + q1 * x22_dim1], ldx22, & v2t[*p + 1 + (*p + 1) * v2t_dim1], ldv2t);
            }
            i__1 = *m - *q;
            i__2 = *m - *q;
            i__3 = *m - *q;
            zungqr_(&i__1, &i__2, &i__3, &v2t[v2t_offset], ldv2t, &work[ itauq2], &work[iorgqr], &lorgqrwork, info);
        }
    }
    /* Compute the CSD of the matrix in bidiagonal-block form */
    zbbcsd_(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, &theta[1], &rwork[ iphi], &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1t[ v1t_offset], ldv1t, &v2t[v2t_offset], ldv2t, &rwork[ib11d], & rwork[ib11e], &rwork[ib12d], &rwork[ib12e], &rwork[ib21d], &rwork[ ib21e], &rwork[ib22d], &rwork[ib22e], &rwork[ibbcsd], &lbbcsdwork, info);
    /* Permute rows and columns to place identity submatrices in top- */
    /* left corner of (1,1)-block and/or bottom-right corner of (1,2)- */
    /* block and/or bottom-right corner of (2,1)-block and/or top-left */
    /* corner of (2,2)-block */
    if (*q > 0 && wantu2)
    {
        i__1 = *q;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            iwork[i__] = *m - *p - *q + i__;
        }
        i__1 = *m - *p;
        for (i__ = *q + 1;
                i__ <= i__1;
                ++i__)
        {
            iwork[i__] = i__ - *q;
        }
        if (colmajor)
        {
            i__1 = *m - *p;
            i__2 = *m - *p;
            zlapmt_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
        }
        else
        {
            i__1 = *m - *p;
            i__2 = *m - *p;
            zlapmr_(&c_false, &i__1, &i__2, &u2[u2_offset], ldu2, &iwork[1]);
        }
    }
    if (*m > 0 && wantv2t)
    {
        i__1 = *p;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            iwork[i__] = *m - *p - *q + i__;
        }
        i__1 = *m - *q;
        for (i__ = *p + 1;
                i__ <= i__1;
                ++i__)
        {
            iwork[i__] = i__ - *p;
        }
        if (! colmajor)
        {
            i__1 = *m - *q;
            i__2 = *m - *q;
            zlapmt_(&c_false, &i__1, &i__2, &v2t[v2t_offset], ldv2t, &iwork[1] );
        }
        else
        {
            i__1 = *m - *q;
            i__2 = *m - *q;
            zlapmr_(&c_false, &i__1, &i__2, &v2t[v2t_offset], ldv2t, &iwork[1] );
        }
    }
    return 0;
    /* End ZUNCSD */
}
/* zuncsd_ */
