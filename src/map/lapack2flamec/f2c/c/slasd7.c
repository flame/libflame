/* ../netlib/slasd7.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SLASD7 merges the two sets of singular values together into a single sorted set. Then it tries to deflate the size of the problem. Used by sbdsdc. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLASD7 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd7. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd7. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd7. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLASD7( ICOMPQ, NL, NR, SQRE, K, D, Z, ZW, VF, VFW, VL, */
/* VLW, ALPHA, BETA, DSIGMA, IDX, IDXP, IDXQ, */
/* PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM, */
/* C, S, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL, */
/* $ NR, SQRE */
/* REAL ALPHA, BETA, C, S */
/* .. */
/* .. Array Arguments .. */
/* INTEGER GIVCOL( LDGCOL, * ), IDX( * ), IDXP( * ), */
/* $ IDXQ( * ), PERM( * ) */
/* REAL D( * ), DSIGMA( * ), GIVNUM( LDGNUM, * ), */
/* $ VF( * ), VFW( * ), VL( * ), VLW( * ), Z( * ), */
/* $ ZW( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASD7 merges the two sets of singular values together into a single */
/* > sorted set. Then it tries to deflate the size of the problem. There */
/* > are two ways in which deflation can occur: when two or more singular */
/* > values are close together or if there is a tiny entry in the Z */
/* > vector. For each such occurrence the order of the related */
/* > secular equation problem is reduced by one. */
/* > */
/* > SLASD7 is called from SLASD6. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ICOMPQ */
/* > \verbatim */
/* > ICOMPQ is INTEGER */
/* > Specifies whether singular vectors are to be computed */
/* > in compact form, as follows: */
/* > = 0: Compute singular values only. */
/* > = 1: Compute singular vectors of upper */
/* > bidiagonal matrix in compact form. */
/* > \endverbatim */
/* > */
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
/* > The bidiagonal matrix has */
/* > N = NL + NR + 1 rows and */
/* > M = N + SQRE >= N columns. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* > K is INTEGER */
/* > Contains the dimension of the non-deflated matrix, this is */
/* > the order of the related secular equation. 1 <= K <=N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension ( N ) */
/* > On entry D contains the singular values of the two submatrices */
/* > to be combined. On exit D contains the trailing (N-K) updated */
/* > singular values (those which were deflated) sorted into */
/* > increasing order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension ( M ) */
/* > On exit Z contains the updating row vector in the secular */
/* > equation. */
/* > \endverbatim */
/* > */
/* > \param[out] ZW */
/* > \verbatim */
/* > ZW is REAL array, dimension ( M ) */
/* > Workspace for Z. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VF */
/* > \verbatim */
/* > VF is REAL array, dimension ( M ) */
/* > On entry, VF(1:NL+1) contains the first components of all */
/* > right singular vectors of the upper block;
and VF(NL+2:M) */
/* > contains the first components of all right singular vectors */
/* > of the lower block. On exit, VF contains the first components */
/* > of all right singular vectors of the bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] VFW */
/* > \verbatim */
/* > VFW is REAL array, dimension ( M ) */
/* > Workspace for VF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* > VL is REAL array, dimension ( M ) */
/* > On entry, VL(1:NL+1) contains the last components of all */
/* > right singular vectors of the upper block;
and VL(NL+2:M) */
/* > contains the last components of all right singular vectors */
/* > of the lower block. On exit, VL contains the last components */
/* > of all right singular vectors of the bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] VLW */
/* > \verbatim */
/* > VLW is REAL array, dimension ( M ) */
/* > Workspace for VL. */
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
/* > \param[out] DSIGMA */
/* > \verbatim */
/* > DSIGMA is REAL array, dimension ( N ) */
/* > Contains a copy of the diagonal elements (K-1 singular values */
/* > and one zero) in the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] IDX */
/* > \verbatim */
/* > IDX is INTEGER array, dimension ( N ) */
/* > This will contain the permutation used to sort the contents of */
/* > D into ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] IDXP */
/* > \verbatim */
/* > IDXP is INTEGER array, dimension ( N ) */
/* > This will contain the permutation used to place deflated */
/* > values of D at the end of the array. On output IDXP(2:K) */
/* > points to the nondeflated D-values and IDXP(K+1:N) */
/* > points to the deflated singular values. */
/* > \endverbatim */
/* > */
/* > \param[in] IDXQ */
/* > \verbatim */
/* > IDXQ is INTEGER array, dimension ( N ) */
/* > This contains the permutation which separately sorts the two */
/* > sub-problems in D into ascending order. Note that entries in */
/* > the first half of this permutation must first be moved one */
/* > position backward;
and entries in the second half */
/* > must first have NL+1 added to their values. */
/* > \endverbatim */
/* > */
/* > \param[out] PERM */
/* > \verbatim */
/* > PERM is INTEGER array, dimension ( N ) */
/* > The permutations (from deflation and sorting) to be applied */
/* > to each singular block. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVPTR */
/* > \verbatim */
/* > GIVPTR is INTEGER */
/* > The number of Givens rotations which took place in this */
/* > subproblem. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVCOL */
/* > \verbatim */
/* > GIVCOL is INTEGER array, dimension ( LDGCOL, 2 ) */
/* > Each pair of numbers indicates a pair of columns to take place */
/* > in a Givens rotation. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGCOL */
/* > \verbatim */
/* > LDGCOL is INTEGER */
/* > The leading dimension of GIVCOL, must be at least N. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVNUM */
/* > \verbatim */
/* > GIVNUM is REAL array, dimension ( LDGNUM, 2 ) */
/* > Each number indicates the C or S value to be used in the */
/* > corresponding Givens rotation. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGNUM */
/* > \verbatim */
/* > LDGNUM is INTEGER */
/* > The leading dimension of GIVNUM, must be at least N. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is REAL */
/* > C contains garbage if SQRE =0 and the C-value of a Givens */
/* > rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is REAL */
/* > S contains garbage if SQRE =0 and the S-value of a Givens */
/* > rotation related to the right null space if SQRE = 1. */
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
int slasd7_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, real *d__, real *z__, real *zw, real *vf, real *vfw, real *vl, real *vlw, real *alpha, real *beta, real *dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm, integer * givptr, integer *givcol, integer *ldgcol, real *givnum, integer * ldgnum, real *c__, real *s, integer *info)
{
    /* System generated locals */
    integer givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, i__1;
    real r__1, r__2;
    /* Local variables */
    integer i__, j, m, n, k2;
    real z1;
    integer jp;
    real eps, tau, tol;
    integer nlp1, nlp2, idxi, idxj;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *);
    integer idxjp, jprev;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *);
    extern real slapy2_(real *, real *), slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *), slamrg_( integer *, integer *, real *, integer *, integer *, integer *);
    real hlftol;
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
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --z__;
    --zw;
    --vf;
    --vfw;
    --vl;
    --vlw;
    --dsigma;
    --idx;
    --idxp;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    /* Function Body */
    *info = 0;
    n = *nl + *nr + 1;
    m = n + *sqre;
    if (*icompq < 0 || *icompq > 1)
    {
        *info = -1;
    }
    else if (*nl < 1)
    {
        *info = -2;
    }
    else if (*nr < 1)
    {
        *info = -3;
    }
    else if (*sqre < 0 || *sqre > 1)
    {
        *info = -4;
    }
    else if (*ldgcol < n)
    {
        *info = -22;
    }
    else if (*ldgnum < n)
    {
        *info = -24;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLASD7", &i__1);
        return 0;
    }
    nlp1 = *nl + 1;
    nlp2 = *nl + 2;
    if (*icompq == 1)
    {
        *givptr = 0;
    }
    /* Generate the first part of the vector Z and move the singular */
    /* values in the first part of D one position backward. */
    z1 = *alpha * vl[nlp1];
    vl[nlp1] = 0.f;
    tau = vf[nlp1];
    for (i__ = *nl;
            i__ >= 1;
            --i__)
    {
        z__[i__ + 1] = *alpha * vl[i__];
        vl[i__] = 0.f;
        vf[i__ + 1] = vf[i__];
        d__[i__ + 1] = d__[i__];
        idxq[i__ + 1] = idxq[i__] + 1;
        /* L10: */
    }
    vf[1] = tau;
    /* Generate the second part of the vector Z. */
    i__1 = m;
    for (i__ = nlp2;
            i__ <= i__1;
            ++i__)
    {
        z__[i__] = *beta * vf[i__];
        vf[i__] = 0.f;
        /* L20: */
    }
    /* Sort the singular values into increasing order */
    i__1 = n;
    for (i__ = nlp2;
            i__ <= i__1;
            ++i__)
    {
        idxq[i__] += nlp1;
        /* L30: */
    }
    /* DSIGMA, IDXC, IDXC, and ZW are used as storage space. */
    i__1 = n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        dsigma[i__] = d__[idxq[i__]];
        zw[i__] = z__[idxq[i__]];
        vfw[i__] = vf[idxq[i__]];
        vlw[i__] = vl[idxq[i__]];
        /* L40: */
    }
    slamrg_(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);
    i__1 = n;
    for (i__ = 2;
            i__ <= i__1;
            ++i__)
    {
        idxi = idx[i__] + 1;
        d__[i__] = dsigma[idxi];
        z__[i__] = zw[idxi];
        vf[i__] = vfw[idxi];
        vl[i__] = vlw[idxi];
        /* L50: */
    }
    /* Calculate the allowable deflation tolerence */
    eps = slamch_("Epsilon");
    /* Computing MAX */
    r__1 = f2c_abs(*alpha);
    r__2 = f2c_abs(*beta); // , expr subst
    tol = max(r__1,r__2);
    /* Computing MAX */
    r__2 = (r__1 = d__[n], f2c_abs(r__1));
    tol = eps * 64.f * max(r__2,tol);
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
            if (j == n)
            {
                goto L100;
            }
        }
        else
        {
            jprev = j;
            goto L70;
        }
        /* L60: */
    }
L70:
    j = jprev;
L80:
    ++j;
    if (j > n)
    {
        goto L90;
    }
    if ((r__1 = z__[j], f2c_abs(r__1)) <= tol)
    {
        /* Deflate due to small z component. */
        --k2;
        idxp[k2] = j;
    }
    else
    {
        /* Check if singular values are close enough to allow deflation. */
        if ((r__1 = d__[j] - d__[jprev], f2c_abs(r__1)) <= tol)
        {
            /* Deflation is possible. */
            *s = z__[jprev];
            *c__ = z__[j];
            /* Find sqrt(a**2+b**2) without overflow or */
            /* destructive underflow. */
            tau = slapy2_(c__, s);
            z__[j] = tau;
            z__[jprev] = 0.f;
            *c__ /= tau;
            *s = -(*s) / tau;
            /* Record the appropriate Givens rotation */
            if (*icompq == 1)
            {
                ++(*givptr);
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
                givcol[*givptr + (givcol_dim1 << 1)] = idxjp;
                givcol[*givptr + givcol_dim1] = idxj;
                givnum[*givptr + (givnum_dim1 << 1)] = *c__;
                givnum[*givptr + givnum_dim1] = *s;
            }
            srot_(&c__1, &vf[jprev], &c__1, &vf[j], &c__1, c__, s);
            srot_(&c__1, &vl[jprev], &c__1, &vl[j], &c__1, c__, s);
            --k2;
            idxp[k2] = jprev;
            jprev = j;
        }
        else
        {
            ++(*k);
            zw[*k] = z__[jprev];
            dsigma[*k] = d__[jprev];
            idxp[*k] = jprev;
            jprev = j;
        }
    }
    goto L80;
L90: /* Record the last singular value. */
    ++(*k);
    zw[*k] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;
L100: /* Sort the singular values into DSIGMA. The singular values which */
    /* were not deflated go into the first K slots of DSIGMA, except */
    /* that DSIGMA(1) is treated separately. */
    i__1 = n;
    for (j = 2;
            j <= i__1;
            ++j)
    {
        jp = idxp[j];
        dsigma[j] = d__[jp];
        vfw[j] = vf[jp];
        vlw[j] = vl[jp];
        /* L110: */
    }
    if (*icompq == 1)
    {
        i__1 = n;
        for (j = 2;
                j <= i__1;
                ++j)
        {
            jp = idxp[j];
            perm[j] = idxq[idx[jp] + 1];
            if (perm[j] <= nlp1)
            {
                --perm[j];
            }
            /* L120: */
        }
    }
    /* The deflated singular values go back into the last N - K slots of */
    /* D. */
    i__1 = n - *k;
    scopy_(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);
    /* Determine DSIGMA(1), DSIGMA(2), Z(1), VF(1), VL(1), VF(M), and */
    /* VL(M). */
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
            *c__ = 1.f;
            *s = 0.f;
            z__[1] = tol;
        }
        else
        {
            *c__ = z1 / z__[1];
            *s = -z__[m] / z__[1];
        }
        srot_(&c__1, &vf[m], &c__1, &vf[1], &c__1, c__, s);
        srot_(&c__1, &vl[m], &c__1, &vl[1], &c__1, c__, s);
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
    /* Restore Z, VF, and VL. */
    i__1 = *k - 1;
    scopy_(&i__1, &zw[2], &c__1, &z__[2], &c__1);
    i__1 = n - 1;
    scopy_(&i__1, &vfw[2], &c__1, &vf[2], &c__1);
    i__1 = n - 1;
    scopy_(&i__1, &vlw[2], &c__1, &vl[2], &c__1);
    return 0;
    /* End of SLASD7 */
}
/* slasd7_ */
