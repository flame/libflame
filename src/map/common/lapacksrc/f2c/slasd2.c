/* ../netlib/slasd2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b30 = 0.f;
/* > \brief \b SLASD2 merges the two sets of singular values together into a single sorted set. Used by sbdsdc . */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLASD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLASD2( NL, NR, SQRE, K, D, Z, ALPHA, BETA, U, LDU, VT, */
/* LDVT, DSIGMA, U2, LDU2, VT2, LDVT2, IDXP, IDX, */
/* IDXC, IDXQ, COLTYP, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, K, LDU, LDU2, LDVT, LDVT2, NL, NR, SQRE */
/* REAL ALPHA, BETA */
/* .. */
/* .. Array Arguments .. */
/* INTEGER COLTYP( * ), IDX( * ), IDXC( * ), IDXP( * ), */
/* $ IDXQ( * ) */
/* REAL D( * ), DSIGMA( * ), U( LDU, * ), */
/* $ U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ), */
/* $ Z( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASD2 merges the two sets of singular values together into a single */
/* > sorted set. Then it tries to deflate the size of the problem. */
/* > There are two ways in which deflation can occur: when two or more */
/* > singular values are close together or if there is a tiny entry in the */
/* > Z vector. For each such occurrence the order of the related secular */
/* > equation problem is reduced by one. */
/* > */
/* > SLASD2 is called from SLASD1. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] NL */
/* > \verbatim */
/* > NL is INTEGER */
/* > The row dimension of the upper block. NL >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] NR */
/* > \verbatim */
/* > NR is INTEGER */
/* > The row dimension of the lower block. NR >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] SQRE */
/* > \verbatim */
/* > SQRE is INTEGER */
/* > = 0: the lower block is an NR-by-NR square matrix. */
/* > = 1: the lower block is an NR-by-(NR+1) rectangular matrix. */
/* > */
/* > The bidiagonal matrix has N = NL + NR + 1 rows and */
/* > M = N + SQRE >= N columns. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* > K is INTEGER */
/* > Contains the dimension of the non-deflated matrix, */
/* > This is the order of the related secular equation. 1 <= K <=N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry D contains the singular values of the two submatrices */
/* > to be combined. On exit D contains the trailing (N-K) updated */
/* > singular values (those which were deflated) sorted into */
/* > increasing order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (N) */
/* > On exit Z contains the updating row vector in the secular */
/* > equation. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is REAL */
/* > Contains the diagonal element associated with the added row. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* > BETA is REAL */
/* > Contains the off-diagonal element associated with the added */
/* > row. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* > U is REAL array, dimension (LDU,N) */
/* > On entry U contains the left singular vectors of two */
/* > submatrices in the two square blocks with corners at (1,1), */
/* > (NL, NL), and (NL+2, NL+2), (N,N). */
/* > On exit U contains the trailing (N-K) updated left singular */
/* > vectors (those which were deflated) in its last N-K columns. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER */
/* > The leading dimension of the array U. LDU >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VT */
/* > \verbatim */
/* > VT is REAL array, dimension (LDVT,M) */
/* > On entry VT**T contains the right singular vectors of two */
/* > submatrices in the two square blocks with corners at (1,1), */
/* > (NL+1, NL+1), and (NL+2, NL+2), (M,M). */
/* > On exit VT**T contains the trailing (N-K) updated right singular */
/* > vectors (those which were deflated) in its last N-K columns. */
/* > In case SQRE =1, the last row of VT spans the right null */
/* > space. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* > LDVT is INTEGER */
/* > The leading dimension of the array VT. LDVT >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] DSIGMA */
/* > \verbatim */
/* > DSIGMA is REAL array, dimension (N) */
/* > Contains a copy of the diagonal elements (K-1 singular values */
/* > and one zero) in the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] U2 */
/* > \verbatim */
/* > U2 is REAL array, dimension (LDU2,N) */
/* > Contains a copy of the first K-1 left singular vectors which */
/* > will be used by SLASD3 in a matrix multiply (SGEMM) to solve */
/* > for the new left singular vectors. U2 is arranged into four */
/* > blocks. The first block contains a column with 1 at NL+1 and */
/* > zero everywhere else;
the second block contains non-zero */
/* > entries only at and above NL;
the third contains non-zero */
/* > entries only below NL+1;
and the fourth is dense. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU2 */
/* > \verbatim */
/* > LDU2 is INTEGER */
/* > The leading dimension of the array U2. LDU2 >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VT2 */
/* > \verbatim */
/* > VT2 is REAL array, dimension (LDVT2,N) */
/* > VT2**T contains a copy of the first K right singular vectors */
/* > which will be used by SLASD3 in a matrix multiply (SGEMM) to */
/* > solve for the new right singular vectors. VT2 is arranged into */
/* > three blocks. The first block contains a row that corresponds */
/* > to the special 0 diagonal element in SIGMA;
the second block */
/* > contains non-zeros only at and before NL +1;
the third block */
/* > contains non-zeros only at and after NL +2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT2 */
/* > \verbatim */
/* > LDVT2 is INTEGER */
/* > The leading dimension of the array VT2. LDVT2 >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] IDXP */
/* > \verbatim */
/* > IDXP is INTEGER array, dimension (N) */
/* > This will contain the permutation used to place deflated */
/* > values of D at the end of the array. On output IDXP(2:K) */
/* > points to the nondeflated D-values and IDXP(K+1:N) */
/* > points to the deflated singular values. */
/* > \endverbatim */
/* > */
/* > \param[out] IDX */
/* > \verbatim */
/* > IDX is INTEGER array, dimension (N) */
/* > This will contain the permutation used to sort the contents of */
/* > D into ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] IDXC */
/* > \verbatim */
/* > IDXC is INTEGER array, dimension (N) */
/* > This will contain the permutation used to arrange the columns */
/* > of the deflated U matrix into three groups: the first group */
/* > contains non-zero entries only at and above NL, the second */
/* > contains non-zero entries only below NL+2, and the third is */
/* > dense. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IDXQ */
/* > \verbatim */
/* > IDXQ is INTEGER array, dimension (N) */
/* > This contains the permutation which separately sorts the two */
/* > sub-problems in D into ascending order. Note that entries in */
/* > the first hlaf of this permutation must first be moved one */
/* > position backward;
and entries in the second half */
/* > must first have NL+1 added to their values. */
/* > \endverbatim */
/* > */
/* > \param[out] COLTYP */
/* > \verbatim */
/* > COLTYP is INTEGER array, dimension (N) */
/* > As workspace, this will contain a label which will indicate */
/* > which of the following types a column in the U2 matrix or a */
/* > row in the VT2 matrix is: */
/* > 1 : non-zero in the upper half only */
/* > 2 : non-zero in the lower half only */
/* > 3 : dense */
/* > 4 : deflated */
/* > */
/* > On exit, it is an array of dimension 4, with COLTYP(I) being */
/* > the dimension of the I-th type columns. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Ming Gu and Huan Ren, Computer Science Division, University of */
/* > California at Berkeley, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
int slasd2_(integer *nl, integer *nr, integer *sqre, integer *k, real *d__, real *z__, real *alpha, real *beta, real *u, integer * ldu, real *vt, integer *ldvt, real *dsigma, real *u2, integer *ldu2, real *vt2, integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info)
{
    /* System generated locals */
    integer u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, vt_offset, vt2_dim1, vt2_offset, i__1;
    real r__1, r__2;
    /* Local variables */
    real c__;
    integer i__, j, m, n;
    real s;
    integer k2;
    real z1;
    integer ct, jp;
    real eps, tau, tol;
    integer psm[4], nlp1, nlp2, idxi, idxj, ctot[4];
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *);
    integer idxjp, jprev;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *);
    extern real slapy2_(real *, real *), slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *), slamrg_( integer *, integer *, real *, integer *, integer *, integer *);
    real hlftol;
    extern /* Subroutine */
    int slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *), slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Arrays .. */
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --z__;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --dsigma;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    vt2_dim1 = *ldvt2;
    vt2_offset = 1 + vt2_dim1;
    vt2 -= vt2_offset;
    --idxp;
    --idx;
    --idxc;
    --idxq;
    --coltyp;
    /* Function Body */
    *info = 0;
    if (*nl < 1)
    {
        *info = -1;
    }
    else if (*nr < 1)
    {
        *info = -2;
    }
    else if (*sqre != 1 && *sqre != 0)
    {
        *info = -3;
    }
    n = *nl + *nr + 1;
    m = n + *sqre;
    if (*ldu < n)
    {
        *info = -10;
    }
    else if (*ldvt < m)
    {
        *info = -12;
    }
    else if (*ldu2 < n)
    {
        *info = -15;
    }
    else if (*ldvt2 < m)
    {
        *info = -17;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLASD2", &i__1);
        return 0;
    }
    nlp1 = *nl + 1;
    nlp2 = *nl + 2;
    /* Generate the first part of the vector Z;
    and move the singular */
    /* values in the first part of D one position backward. */
    z1 = *alpha * vt[nlp1 + nlp1 * vt_dim1];
    z__[1] = z1;
    for (i__ = *nl;
            i__ >= 1;
            --i__)
    {
        z__[i__ + 1] = *alpha * vt[i__ + nlp1 * vt_dim1];
        d__[i__ + 1] = d__[i__];
        idxq[i__ + 1] = idxq[i__] + 1;
        /* L10: */
    }
    /* Generate the second part of the vector Z. */
    i__1 = m;
    for (i__ = nlp2;
            i__ <= i__1;
            ++i__)
    {
        z__[i__] = *beta * vt[i__ + nlp2 * vt_dim1];
        /* L20: */
    }
    /* Initialize some reference arrays. */
    i__1 = nlp1;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        coltyp[i__] = 1;
        /* L30: */
    }
    i__1 = n;
    for (i__ = nlp2;
            i__ <= i__1;
            ++i__)
    {
        coltyp[i__] = 2;
        /* L40: */
    }
    /* Sort the singular values into increasing order */
    i__1 = n;
    for (i__ = nlp2;
            i__ <= i__1;
            ++i__)
    {
        idxq[i__] += nlp1;
        /* L50: */
    }
    /* DSIGMA, IDXC, IDXC, and the first column of U2 */
    /* are used as storage space. */
    i__1 = n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        dsigma[i__] = d__[idxq[i__]];
        u2[i__ + u2_dim1] = z__[idxq[i__]];
        idxc[i__] = coltyp[idxq[i__]];
        /* L60: */
    }
    slamrg_(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);
    i__1 = n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        idxi = idx[i__] + 1;
        d__[i__] = dsigma[idxi];
        z__[i__] = u2[idxi + u2_dim1];
        coltyp[i__] = idxc[idxi];
        /* L70: */
    }
    /* Calculate the allowable deflation tolerance */
    eps = slamch_("Epsilon");
    /* Computing MAX */
    r__1 = f2c_abs(*alpha);
    r__2 = f2c_abs(*beta); // , expr subst
    tol = max(r__1,r__2);
    /* Computing MAX */
    r__2 = (r__1 = d__[n], f2c_abs(r__1));
    tol = eps * 8.f * max(r__2,tol);
    /* There are 2 kinds of deflation -- first a value in the z-vector */
    /* is small, second two (or more) singular values are very close */
    /* together (their difference is small). */
    /* If the value in the z-vector is small, we simply permute the */
    /* array so that the corresponding singular value is moved to the */
    /* end. */
    /* If two values in the D-vector are close, we perform a two-sided */
    /* rotation designed to make one of the corresponding z-vector */
    /* entries zero, and then permute the array so that the deflated */
    /* singular value is moved to the end. */
    /* If there are multiple singular values then the problem deflates. */
    /* Here the number of equal singular values are found. As each equal */
    /* singular value is found, an elementary reflector is computed to */
    /* rotate the corresponding singular subspace so that the */
    /* corresponding components of Z are zero in this new basis. */
    *k = 1;
    k2 = n + 1;
    i__1 = n;
    for (j = 2;
            j <= i__1;
            ++j)
    {
        if ((r__1 = z__[j], f2c_abs(r__1)) <= tol)
        {
            /* Deflate due to small z component. */
            --k2;
            idxp[k2] = j;
            coltyp[j] = 4;
            if (j == n)
            {
                goto L120;
            }
        }
        else
        {
            jprev = j;
            goto L90;
        }
        /* L80: */
    }
L90:
    j = jprev;
L100:
    ++j;
    if (j > n)
    {
        goto L110;
    }
    if ((r__1 = z__[j], f2c_abs(r__1)) <= tol)
    {
        /* Deflate due to small z component. */
        --k2;
        idxp[k2] = j;
        coltyp[j] = 4;
    }
    else
    {
        /* Check if singular values are close enough to allow deflation. */
        if ((r__1 = d__[j] - d__[jprev], f2c_abs(r__1)) <= tol)
        {
            /* Deflation is possible. */
            s = z__[jprev];
            c__ = z__[j];
            /* Find sqrt(a**2+b**2) without overflow or */
            /* destructive underflow. */
            tau = slapy2_(&c__, &s);
            c__ /= tau;
            s = -s / tau;
            z__[j] = tau;
            z__[jprev] = 0.f;
            /* Apply back the Givens rotation to the left and right */
            /* singular vector matrices. */
            idxjp = idxq[idx[jprev] + 1];
            idxj = idxq[idx[j] + 1];
            if (idxjp <= nlp1)
            {
                --idxjp;
            }
            if (idxj <= nlp1)
            {
                --idxj;
            }
            srot_(&n, &u[idxjp * u_dim1 + 1], &c__1, &u[idxj * u_dim1 + 1], & c__1, &c__, &s);
            srot_(&m, &vt[idxjp + vt_dim1], ldvt, &vt[idxj + vt_dim1], ldvt, & c__, &s);
            if (coltyp[j] != coltyp[jprev])
            {
                coltyp[j] = 3;
            }
            coltyp[jprev] = 4;
            --k2;
            idxp[k2] = jprev;
            jprev = j;
        }
        else
        {
            ++(*k);
            u2[*k + u2_dim1] = z__[jprev];
            dsigma[*k] = d__[jprev];
            idxp[*k] = jprev;
            jprev = j;
        }
    }
    goto L100;
L110: /* Record the last singular value. */
    ++(*k);
    u2[*k + u2_dim1] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;
L120: /* Count up the total number of the various types of columns, then */
    /* form a permutation which positions the four column types into */
    /* four groups of uniform structure (although one or more of these */
    /* groups may be empty). */
    for (j = 1;
            j <= 4;
            ++j)
    {
        ctot[j - 1] = 0;
        /* L130: */
    }
    i__1 = n;
    for (j = 2;
            j <= i__1;
            ++j)
    {
        ct = coltyp[j];
        ++ctot[ct - 1];
        /* L140: */
    }
    /* PSM(*) = Position in SubMatrix (of types 1 through 4) */
    psm[0] = 2;
    psm[1] = ctot[0] + 2;
    psm[2] = psm[1] + ctot[1];
    psm[3] = psm[2] + ctot[2];
    /* Fill out the IDXC array so that the permutation which it induces */
    /* will place all type-1 columns first, all type-2 columns next, */
    /* then all type-3's, and finally all type-4's, starting from the */
    /* second column. This applies similarly to the rows of VT. */
    i__1 = n;
    for (j = 2;
            j <= i__1;
            ++j)
    {
        jp = idxp[j];
        ct = coltyp[jp];
        idxc[psm[ct - 1]] = j;
        ++psm[ct - 1];
        /* L150: */
    }
    /* Sort the singular values and corresponding singular vectors into */
    /* DSIGMA, U2, and VT2 respectively. The singular values/vectors */
    /* which were not deflated go into the first K slots of DSIGMA, U2, */
    /* and VT2 respectively, while those which were deflated go into the */
    /* last N - K slots, except that the first column/row will be treated */
    /* separately. */
    i__1 = n;
    for (j = 2;
            j <= i__1;
            ++j)
    {
        jp = idxp[j];
        dsigma[j] = d__[jp];
        idxj = idxq[idx[idxp[idxc[j]]] + 1];
        if (idxj <= nlp1)
        {
            --idxj;
        }
        scopy_(&n, &u[idxj * u_dim1 + 1], &c__1, &u2[j * u2_dim1 + 1], &c__1);
        scopy_(&m, &vt[idxj + vt_dim1], ldvt, &vt2[j + vt2_dim1], ldvt2);
        /* L160: */
    }
    /* Determine DSIGMA(1), DSIGMA(2) and Z(1) */
    dsigma[1] = 0.f;
    hlftol = tol / 2.f;
    if (f2c_abs(dsigma[2]) <= hlftol)
    {
        dsigma[2] = hlftol;
    }
    if (m > n)
    {
        z__[1] = slapy2_(&z1, &z__[m]);
        if (z__[1] <= tol)
        {
            c__ = 1.f;
            s = 0.f;
            z__[1] = tol;
        }
        else
        {
            c__ = z1 / z__[1];
            s = z__[m] / z__[1];
        }
    }
    else
    {
        if (f2c_abs(z1) <= tol)
        {
            z__[1] = tol;
        }
        else
        {
            z__[1] = z1;
        }
    }
    /* Move the rest of the updating row to Z. */
    i__1 = *k - 1;
    scopy_(&i__1, &u2[u2_dim1 + 2], &c__1, &z__[2], &c__1);
    /* Determine the first column of U2, the first row of VT2 and the */
    /* last row of VT. */
    slaset_("A", &n, &c__1, &c_b30, &c_b30, &u2[u2_offset], ldu2);
    u2[nlp1 + u2_dim1] = 1.f;
    if (m > n)
    {
        i__1 = nlp1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            vt[m + i__ * vt_dim1] = -s * vt[nlp1 + i__ * vt_dim1];
            vt2[i__ * vt2_dim1 + 1] = c__ * vt[nlp1 + i__ * vt_dim1];
            /* L170: */
        }
        i__1 = m;
        for (i__ = nlp2;
                i__ <= i__1;
                ++i__)
        {
            vt2[i__ * vt2_dim1 + 1] = s * vt[m + i__ * vt_dim1];
            vt[m + i__ * vt_dim1] = c__ * vt[m + i__ * vt_dim1];
            /* L180: */
        }
    }
    else
    {
        scopy_(&m, &vt[nlp1 + vt_dim1], ldvt, &vt2[vt2_dim1 + 1], ldvt2);
    }
    if (m > n)
    {
        scopy_(&m, &vt[m + vt_dim1], ldvt, &vt2[m + vt2_dim1], ldvt2);
    }
    /* The deflated singular values and their corresponding vectors go */
    /* into the back of D, U, and V respectively. */
    if (n > *k)
    {
        i__1 = n - *k;
        scopy_(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);
        i__1 = n - *k;
        slacpy_("A", &n, &i__1, &u2[(*k + 1) * u2_dim1 + 1], ldu2, &u[(*k + 1) * u_dim1 + 1], ldu);
        i__1 = n - *k;
        slacpy_("A", &i__1, &m, &vt2[*k + 1 + vt2_dim1], ldvt2, &vt[*k + 1 + vt_dim1], ldvt);
    }
    /* Copy CTOT into COLTYP for referencing in SLASD3. */
    for (j = 1;
            j <= 4;
            ++j)
    {
        coltyp[j] = ctot[j - 1];
        /* L190: */
    }
    return 0;
    /* End of SLASD2 */
}
/* slasd2_ */
