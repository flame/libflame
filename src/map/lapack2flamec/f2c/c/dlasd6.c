/* ../netlib/dlasd6.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__0 = 0;
static doublereal c_b7 = 1.;
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b DLASD6 computes the SVD of an updated upper bidiagonal matrix obtained by merging two smaller o nes by appending a row. Used by sbdsdc. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLASD6 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd6. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd6. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd6. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLASD6( ICOMPQ, NL, NR, SQRE, D, VF, VL, ALPHA, BETA, */
/* IDXQ, PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, */
/* LDGNUM, POLES, DIFL, DIFR, Z, K, C, S, WORK, */
/* IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER GIVPTR, ICOMPQ, INFO, K, LDGCOL, LDGNUM, NL, */
/* $ NR, SQRE */
/* DOUBLE PRECISION ALPHA, BETA, C, S */
/* .. */
/* .. Array Arguments .. */
/* INTEGER GIVCOL( LDGCOL, * ), IDXQ( * ), IWORK( * ), */
/* $ PERM( * ) */
/* DOUBLE PRECISION D( * ), DIFL( * ), DIFR( * ), */
/* $ GIVNUM( LDGNUM, * ), POLES( LDGNUM, * ), */
/* $ VF( * ), VL( * ), WORK( * ), Z( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASD6 computes the SVD of an updated upper bidiagonal matrix B */
/* > obtained by merging two smaller ones by appending a row. This */
/* > routine is used only for the problem which requires all singular */
/* > values and optionally singular vector matrices in factored form. */
/* > B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE. */
/* > A related subroutine, DLASD1, handles the case in which all singular */
/* > values and singular vectors of the bidiagonal matrix are desired. */
/* > */
/* > DLASD6 computes the SVD as follows: */
/* > */
/* > ( D1(in) 0 0 0 ) */
/* > B = U(in) * ( Z1**T a Z2**T b ) * VT(in) */
/* > ( 0 0 D2(in) 0 ) */
/* > */
/* > = U(out) * ( D(out) 0) * VT(out) */
/* > */
/* > where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M */
/* > with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros */
/* > elsewhere;
and the entry b is empty if SQRE = 0. */
/* > */
/* > The singular values of B can be computed using D1, D2, the first */
/* > components of all the right singular vectors of the lower block, and */
/* > the last components of all the right singular vectors of the upper */
/* > block. These components are stored and updated in VF and VL, */
/* > respectively, in DLASD6. Hence U and VT are not explicitly */
/* > referenced. */
/* > */
/* > The singular values are stored in D. The algorithm consists of two */
/* > stages: */
/* > */
/* > The first stage consists of deflating the size of the problem */
/* > when there are multiple singular values or if there is a zero */
/* > in the Z vector. For each such occurence the dimension of the */
/* > secular equation problem is reduced by one. This stage is */
/* > performed by the routine DLASD7. */
/* > */
/* > The second stage consists of calculating the updated */
/* > singular values. This is done by finding the roots of the */
/* > secular equation via the routine DLASD4 (as called by DLASD8). */
/* > This routine also updates VF and VL and computes the distances */
/* > between the updated singular values and the old singular */
/* > values. */
/* > */
/* > DLASD6 is called from DLASDA. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ICOMPQ */
/* > \verbatim */
/* > ICOMPQ is INTEGER */
/* > Specifies whether singular vectors are to be computed in */
/* > factored form: */
/* > = 0: Compute singular values only. */
/* > = 1: Compute singular vectors in factored form as well. */
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
/* > The bidiagonal matrix has row dimension N = NL + NR + 1, */
/* > and column dimension M = N + SQRE. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension ( NL+NR+1 ). */
/* > On entry D(1:NL,1:NL) contains the singular values of the */
/* > upper block, and D(NL+2:N) contains the singular values */
/* > of the lower block. On exit D(1:N) contains the singular */
/* > values of the modified matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VF */
/* > \verbatim */
/* > VF is DOUBLE PRECISION array, dimension ( M ) */
/* > On entry, VF(1:NL+1) contains the first components of all */
/* > right singular vectors of the upper block;
and VF(NL+2:M) */
/* > contains the first components of all right singular vectors */
/* > of the lower block. On exit, VF contains the first components */
/* > of all right singular vectors of the bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* > VL is DOUBLE PRECISION array, dimension ( M ) */
/* > On entry, VL(1:NL+1) contains the last components of all */
/* > right singular vectors of the upper block;
and VL(NL+2:M) */
/* > contains the last components of all right singular vectors of */
/* > the lower block. On exit, VL contains the last components of */
/* > all right singular vectors of the bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ALPHA */
/* > \verbatim */
/* > ALPHA is DOUBLE PRECISION */
/* > Contains the diagonal element associated with the added row. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BETA */
/* > \verbatim */
/* > BETA is DOUBLE PRECISION */
/* > Contains the off-diagonal element associated with the added */
/* > row. */
/* > \endverbatim */
/* > */
/* > \param[out] IDXQ */
/* > \verbatim */
/* > IDXQ is INTEGER array, dimension ( N ) */
/* > This contains the permutation which will reintegrate the */
/* > subproblem just solved back into sorted order, i.e. */
/* > D( IDXQ( I = 1, N ) ) will be in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] PERM */
/* > \verbatim */
/* > PERM is INTEGER array, dimension ( N ) */
/* > The permutations (from deflation and sorting) to be applied */
/* > to each block. Not referenced if ICOMPQ = 0. */
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
/* > leading dimension of GIVCOL, must be at least N. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVNUM */
/* > \verbatim */
/* > GIVNUM is DOUBLE PRECISION array, dimension ( LDGNUM, 2 ) */
/* > Each number indicates the C or S value to be used in the */
/* > corresponding Givens rotation. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGNUM */
/* > \verbatim */
/* > LDGNUM is INTEGER */
/* > The leading dimension of GIVNUM and POLES, must be at least N. */
/* > \endverbatim */
/* > */
/* > \param[out] POLES */
/* > \verbatim */
/* > POLES is DOUBLE PRECISION array, dimension ( LDGNUM, 2 ) */
/* > On exit, POLES(1,*) is an array containing the new singular */
/* > values obtained from solving the secular equation, and */
/* > POLES(2,*) is an array containing the poles in the secular */
/* > equation. Not referenced if ICOMPQ = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] DIFL */
/* > \verbatim */
/* > DIFL is DOUBLE PRECISION array, dimension ( N ) */
/* > On exit, DIFL(I) is the distance between I-th updated */
/* > (undeflated) singular value and the I-th (undeflated) old */
/* > singular value. */
/* > \endverbatim */
/* > */
/* > \param[out] DIFR */
/* > \verbatim */
/* > DIFR is DOUBLE PRECISION array, */
/* > dimension ( LDGNUM, 2 ) if ICOMPQ = 1 and */
/* > dimension ( N ) if ICOMPQ = 0. */
/* > On exit, DIFR(I, 1) is the distance between I-th updated */
/* > (undeflated) singular value and the I+1-th (undeflated) old */
/* > singular value. */
/* > */
/* > If ICOMPQ = 1, DIFR(1:K,2) is an array containing the */
/* > normalizing factors for the right singular vector matrix. */
/* > */
/* > See DLASD8 for details on DIFL and DIFR. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension ( M ) */
/* > The first elements of this array contain the components */
/* > of the deflation-adjusted updating row vector. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* > K is INTEGER */
/* > Contains the dimension of the non-deflated matrix, */
/* > This is the order of the related secular equation. 1 <= K <=N. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION */
/* > C contains garbage if SQRE =0 and the C-value of a Givens */
/* > rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is DOUBLE PRECISION */
/* > S contains garbage if SQRE =0 and the S-value of a Givens */
/* > rotation related to the right null space if SQRE = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension ( 4 * M ) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension ( 3 * N ) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = 1, a singular value did not converge */
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
int dlasd6_(integer *icompq, integer *nl, integer *nr, integer *sqre, doublereal *d__, doublereal *vf, doublereal *vl, doublereal *alpha, doublereal *beta, integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *poles, doublereal *difl, doublereal * difr, doublereal *z__, integer *k, doublereal *c__, doublereal *s, doublereal *work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, poles_dim1, poles_offset, i__1;
    doublereal d__1, d__2;
    /* Local variables */
    integer i__, m, n, n1, n2, iw, idx, idxc, idxp, ivfw, ivlw;
    extern /* Subroutine */
    int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *), dlasd7_(integer *, integer *, integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, integer *, integer *, integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), dlasd8_( integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, integer *), dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *), dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *);
    integer isigma;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal orgnrm;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --vf;
    --vl;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    poles_dim1 = *ldgnum;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    --difl;
    --difr;
    --z__;
    --work;
    --iwork;
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
        *info = -14;
    }
    else if (*ldgnum < n)
    {
        *info = -16;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DLASD6", &i__1);
        return 0;
    }
    /* The following values are for bookkeeping purposes only. They are */
    /* integer pointers which indicate the portion of the workspace */
    /* used by a particular array in DLASD7 and DLASD8. */
    isigma = 1;
    iw = isigma + n;
    ivfw = iw + m;
    ivlw = ivfw + m;
    idx = 1;
    idxc = idx + n;
    idxp = idxc + n;
    /* Scale. */
    /* Computing MAX */
    d__1 = f2c_dabs(*alpha);
    d__2 = f2c_dabs(*beta); // , expr subst
    orgnrm = max(d__1,d__2);
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if ((d__1 = d__[i__], f2c_dabs(d__1)) > orgnrm)
        {
            orgnrm = (d__1 = d__[i__], f2c_dabs(d__1));
        }
        /* L10: */
    }
    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b7, &n, &c__1, &d__[1], &n, info);
    *alpha /= orgnrm;
    *beta /= orgnrm;
    /* Sort and Deflate singular values. */
    dlasd7_(icompq, nl, nr, sqre, k, &d__[1], &z__[1], &work[iw], &vf[1], & work[ivfw], &vl[1], &work[ivlw], alpha, beta, &work[isigma], & iwork[idx], &iwork[idxp], &idxq[1], &perm[1], givptr, &givcol[ givcol_offset], ldgcol, &givnum[givnum_offset], ldgnum, c__, s, info);
    /* Solve Secular Equation, compute DIFL, DIFR, and update VF, VL. */
    dlasd8_(icompq, k, &d__[1], &z__[1], &vf[1], &vl[1], &difl[1], &difr[1], ldgnum, &work[isigma], &work[iw], info);
    /* Handle error returned */
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DLASD8", &i__1);
        return 0;
    }
    /* Save the poles if ICOMPQ = 1. */
    if (*icompq == 1)
    {
        dcopy_(k, &d__[1], &c__1, &poles[poles_dim1 + 1], &c__1);
        dcopy_(k, &work[isigma], &c__1, &poles[(poles_dim1 << 1) + 1], &c__1);
    }
    /* Unscale. */
    dlascl_("G", &c__0, &c__0, &c_b7, &orgnrm, &n, &c__1, &d__[1], &n, info);
    /* Prepare the IDXQ sorting permutation. */
    n1 = *k;
    n2 = n - *k;
    dlamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);
    return 0;
    /* End of DLASD6 */
}
/* dlasd6_ */
