/* ../netlib/slaeda.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c__1 = 1;
static real c_b24 = 1.f;
static real c_b26 = 0.f;
/* > \brief \b SLAEDA used by sstedc. Computes the Z vector determining the rank-one modification of the diago nal matrix. Used when the original matrix is dense. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAEDA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaeda. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaeda. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaeda. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, */
/* GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER CURLVL, CURPBM, INFO, N, TLVLS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER GIVCOL( 2, * ), GIVPTR( * ), PERM( * ), */
/* $ PRMPTR( * ), QPTR( * ) */
/* REAL GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAEDA computes the Z vector corresponding to the merge step in the */
/* > CURLVLth step of the merge process with TLVLS steps for the CURPBMth */
/* > problem. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The dimension of the symmetric tridiagonal matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TLVLS */
/* > \verbatim */
/* > TLVLS is INTEGER */
/* > The total number of merging levels in the overall divide and */
/* > conquer tree. */
/* > \endverbatim */
/* > */
/* > \param[in] CURLVL */
/* > \verbatim */
/* > CURLVL is INTEGER */
/* > The current level in the overall merge routine, */
/* > 0 <= curlvl <= tlvls. */
/* > \endverbatim */
/* > */
/* > \param[in] CURPBM */
/* > \verbatim */
/* > CURPBM is INTEGER */
/* > The current problem in the current level in the overall */
/* > merge routine (counting from upper left to lower right). */
/* > \endverbatim */
/* > */
/* > \param[in] PRMPTR */
/* > \verbatim */
/* > PRMPTR is INTEGER array, dimension (N lg N) */
/* > Contains a list of pointers which indicate where in PERM a */
/* > level's permutation is stored. PRMPTR(i+1) - PRMPTR(i) */
/* > indicates the size of the permutation and incidentally the */
/* > size of the full, non-deflated problem. */
/* > \endverbatim */
/* > */
/* > \param[in] PERM */
/* > \verbatim */
/* > PERM is INTEGER array, dimension (N lg N) */
/* > Contains the permutations (from deflation and sorting) to be */
/* > applied to each eigenblock. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVPTR */
/* > \verbatim */
/* > GIVPTR is INTEGER array, dimension (N lg N) */
/* > Contains a list of pointers which indicate where in GIVCOL a */
/* > level's Givens rotations are stored. GIVPTR(i+1) - GIVPTR(i) */
/* > indicates the number of Givens rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVCOL */
/* > \verbatim */
/* > GIVCOL is INTEGER array, dimension (2, N lg N) */
/* > Each pair of numbers indicates a pair of columns to take place */
/* > in a Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVNUM */
/* > \verbatim */
/* > GIVNUM is REAL array, dimension (2, N lg N) */
/* > Each number indicates the S value to be used in the */
/* > corresponding Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (N**2) */
/* > Contains the square eigenblocks from previous levels, the */
/* > starting positions for blocks are given by QPTR. */
/* > \endverbatim */
/* > */
/* > \param[in] QPTR */
/* > \verbatim */
/* > QPTR is INTEGER array, dimension (N+2) */
/* > Contains a list of pointers which indicate where in Q an */
/* > eigenblock is stored. SQRT( QPTR(i+1) - QPTR(i) ) indicates */
/* > the size of the block. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (N) */
/* > On output this vector contains the updating vector (the last */
/* > row of the first sub-eigenvector matrix and the first row of */
/* > the second sub-eigenvector matrix). */
/* > \endverbatim */
/* > */
/* > \param[out] ZTEMP */
/* > \verbatim */
/* > ZTEMP is REAL array, dimension (N) */
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
/* > \ingroup auxOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */
/* ===================================================================== */
/* Subroutine */
int slaeda_(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr, integer *perm, integer *givptr, integer *givcol, real *givnum, real *q, integer *qptr, real *z__, real *ztemp, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);
    /* Local variables */
    integer i__, k, mid, ptr, curr;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *);
    integer bsiz1, bsiz2, psiz1, psiz2, zptr1;
    extern /* Subroutine */
    int sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), scopy_(integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *);
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    --ztemp;
    --z__;
    --qptr;
    --q;
    givnum -= 3;
    givcol -= 3;
    --givptr;
    --perm;
    --prmptr;
    /* Function Body */
    *info = 0;
    if (*n < 0)
    {
        *info = -1;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLAEDA", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Determine location of first number in second half. */
    mid = *n / 2 + 1;
    /* Gather last/first rows of appropriate eigenblocks into center of Z */
    ptr = 1;
    /* Determine location of lowest level subproblem in the full storage */
    /* scheme */
    i__1 = *curlvl - 1;
    curr = ptr + *curpbm * pow_ii(&c__2, curlvl) + pow_ii(&c__2, &i__1) - 1;
    /* Determine size of these matrices. We add HALF to the value of */
    /* the SQRT in case the machine underestimates one of these square */
    /* roots. */
    bsiz1 = (integer) (sqrt((real) (qptr[curr + 1] - qptr[curr])) + .5f);
    bsiz2 = (integer) (sqrt((real) (qptr[curr + 2] - qptr[curr + 1])) + .5f);
    i__1 = mid - bsiz1 - 1;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        z__[k] = 0.f;
        /* L10: */
    }
    scopy_(&bsiz1, &q[qptr[curr] + bsiz1 - 1], &bsiz1, &z__[mid - bsiz1], & c__1);
    scopy_(&bsiz2, &q[qptr[curr + 1]], &bsiz2, &z__[mid], &c__1);
    i__1 = *n;
    for (k = mid + bsiz2;
            k <= i__1;
            ++k)
    {
        z__[k] = 0.f;
        /* L20: */
    }
    /* Loop through remaining levels 1 -> CURLVL applying the Givens */
    /* rotations and permutation and then multiplying the center matrices */
    /* against the current Z. */
    ptr = pow_ii(&c__2, tlvls) + 1;
    i__1 = *curlvl - 1;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        i__2 = *curlvl - k;
        i__3 = *curlvl - k - 1;
        curr = ptr + *curpbm * pow_ii(&c__2, &i__2) + pow_ii(&c__2, &i__3) - 1;
        psiz1 = prmptr[curr + 1] - prmptr[curr];
        psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
        zptr1 = mid - psiz1;
        /* Apply Givens at CURR and CURR+1 */
        i__2 = givptr[curr + 1] - 1;
        for (i__ = givptr[curr];
                i__ <= i__2;
                ++i__)
        {
            srot_(&c__1, &z__[zptr1 + givcol[(i__ << 1) + 1] - 1], &c__1, & z__[zptr1 + givcol[(i__ << 1) + 2] - 1], &c__1, &givnum[( i__ << 1) + 1], &givnum[(i__ << 1) + 2]);
            /* L30: */
        }
        i__2 = givptr[curr + 2] - 1;
        for (i__ = givptr[curr + 1];
                i__ <= i__2;
                ++i__)
        {
            srot_(&c__1, &z__[mid - 1 + givcol[(i__ << 1) + 1]], &c__1, &z__[ mid - 1 + givcol[(i__ << 1) + 2]], &c__1, &givnum[(i__ << 1) + 1], &givnum[(i__ << 1) + 2]);
            /* L40: */
        }
        psiz1 = prmptr[curr + 1] - prmptr[curr];
        psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
        i__2 = psiz1 - 1;
        for (i__ = 0;
                i__ <= i__2;
                ++i__)
        {
            ztemp[i__ + 1] = z__[zptr1 + perm[prmptr[curr] + i__] - 1];
            /* L50: */
        }
        i__2 = psiz2 - 1;
        for (i__ = 0;
                i__ <= i__2;
                ++i__)
        {
            ztemp[psiz1 + i__ + 1] = z__[mid + perm[prmptr[curr + 1] + i__] - 1];
            /* L60: */
        }
        /* Multiply Blocks at CURR and CURR+1 */
        /* Determine size of these matrices. We add HALF to the value of */
        /* the SQRT in case the machine underestimates one of these */
        /* square roots. */
        bsiz1 = (integer) (sqrt((real) (qptr[curr + 1] - qptr[curr])) + .5f);
        bsiz2 = (integer) (sqrt((real) (qptr[curr + 2] - qptr[curr + 1])) + .5f);
        if (bsiz1 > 0)
        {
            sgemv_("T", &bsiz1, &bsiz1, &c_b24, &q[qptr[curr]], &bsiz1, & ztemp[1], &c__1, &c_b26, &z__[zptr1], &c__1);
        }
        i__2 = psiz1 - bsiz1;
        scopy_(&i__2, &ztemp[bsiz1 + 1], &c__1, &z__[zptr1 + bsiz1], &c__1);
        if (bsiz2 > 0)
        {
            sgemv_("T", &bsiz2, &bsiz2, &c_b24, &q[qptr[curr + 1]], &bsiz2, & ztemp[psiz1 + 1], &c__1, &c_b26, &z__[mid], &c__1);
        }
        i__2 = psiz2 - bsiz2;
        scopy_(&i__2, &ztemp[psiz1 + bsiz2 + 1], &c__1, &z__[mid + bsiz2], & c__1);
        i__2 = *tlvls - k;
        ptr += pow_ii(&c__2, &i__2);
        /* L70: */
    }
    return 0;
    /* End of SLAEDA */
}
/* slaeda_ */
