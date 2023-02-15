/* ../netlib/clalsa.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b9 = 1.f;
static real c_b10 = 0.f;
static integer c__2 = 2;
/* > \brief \b CLALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLALSA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clalsa. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clalsa. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clalsa. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U, */
/* LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR, */
/* GIVCOL, LDGCOL, PERM, GIVNUM, C, S, RWORK, */
/* IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS, */
/* $ SMLSIZ */
/* .. */
/* .. Array Arguments .. */
/* INTEGER GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ), */
/* $ K( * ), PERM( LDGCOL, * ) */
/* REAL C( * ), DIFL( LDU, * ), DIFR( LDU, * ), */
/* $ GIVNUM( LDU, * ), POLES( LDU, * ), RWORK( * ), */
/* $ S( * ), U( LDU, * ), VT( LDU, * ), Z( LDU, * ) */
/* COMPLEX B( LDB, * ), BX( LDBX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLALSA is an itermediate step in solving the least squares problem */
/* > by computing the SVD of the coefficient matrix in compact form (The */
/* > singular vectors are computed as products of simple orthorgonal */
/* > matrices.). */
/* > */
/* > If ICOMPQ = 0, CLALSA applies the inverse of the left singular vector */
/* > matrix of an upper bidiagonal matrix to the right hand side;
and if */
/* > ICOMPQ = 1, CLALSA applies the right singular vector matrix to the */
/* > right hand side. The singular vector matrices were generated in */
/* > compact form by CLALSA. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ICOMPQ */
/* > \verbatim */
/* > ICOMPQ is INTEGER */
/* > Specifies whether the left or the right singular vector */
/* > matrix is involved. */
/* > = 0: Left singular vector matrix */
/* > = 1: Right singular vector matrix */
/* > \endverbatim */
/* > */
/* > \param[in] SMLSIZ */
/* > \verbatim */
/* > SMLSIZ is INTEGER */
/* > The maximum size of the subproblems at the bottom of the */
/* > computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The row and column dimensions of the upper bidiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of columns of B and BX. NRHS must be at least 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension ( LDB, NRHS ) */
/* > On input, B contains the right hand sides of the least */
/* > squares problem in rows 1 through M. */
/* > On output, B contains the solution X in rows 1 through N. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B in the calling subprogram. */
/* > LDB must be at least max(1,MAX( M, N ) ). */
/* > \endverbatim */
/* > */
/* > \param[out] BX */
/* > \verbatim */
/* > BX is COMPLEX array, dimension ( LDBX, NRHS ) */
/* > On exit, the result of applying the left or right singular */
/* > vector matrix to B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDBX */
/* > \verbatim */
/* > LDBX is INTEGER */
/* > The leading dimension of BX. */
/* > \endverbatim */
/* > */
/* > \param[in] U */
/* > \verbatim */
/* > U is REAL array, dimension ( LDU, SMLSIZ ). */
/* > On entry, U contains the left singular vector matrices of all */
/* > subproblems at the bottom level. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER, LDU = > N. */
/* > The leading dimension of arrays U, VT, DIFL, DIFR, */
/* > POLES, GIVNUM, and Z. */
/* > \endverbatim */
/* > */
/* > \param[in] VT */
/* > \verbatim */
/* > VT is REAL array, dimension ( LDU, SMLSIZ+1 ). */
/* > On entry, VT**H contains the right singular vector matrices of */
/* > all subproblems at the bottom level. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER array, dimension ( N ). */
/* > \endverbatim */
/* > */
/* > \param[in] DIFL */
/* > \verbatim */
/* > DIFL is REAL array, dimension ( LDU, NLVL ). */
/* > where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1. */
/* > \endverbatim */
/* > */
/* > \param[in] DIFR */
/* > \verbatim */
/* > DIFR is REAL array, dimension ( LDU, 2 * NLVL ). */
/* > On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record */
/* > distances between singular values on the I-th level and */
/* > singular values on the (I -1)-th level, and DIFR(*, 2 * I) */
/* > record the normalizing factors of the right singular vectors */
/* > matrices of subproblems on I-th level. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* > Z is REAL array, dimension ( LDU, NLVL ). */
/* > On entry, Z(1, I) contains the components of the deflation- */
/* > adjusted updating row vector for subproblems on the I-th */
/* > level. */
/* > \endverbatim */
/* > */
/* > \param[in] POLES */
/* > \verbatim */
/* > POLES is REAL array, dimension ( LDU, 2 * NLVL ). */
/* > On entry, POLES(*, 2 * I -1: 2 * I) contains the new and old */
/* > singular values involved in the secular equations on the I-th */
/* > level. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVPTR */
/* > \verbatim */
/* > GIVPTR is INTEGER array, dimension ( N ). */
/* > On entry, GIVPTR( I ) records the number of Givens */
/* > rotations performed on the I-th problem on the computation */
/* > tree. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVCOL */
/* > \verbatim */
/* > GIVCOL is INTEGER array, dimension ( LDGCOL, 2 * NLVL ). */
/* > On entry, for each I, GIVCOL(*, 2 * I - 1: 2 * I) records the */
/* > locations of Givens rotations performed on the I-th level on */
/* > the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] LDGCOL */
/* > \verbatim */
/* > LDGCOL is INTEGER, LDGCOL = > N. */
/* > The leading dimension of arrays GIVCOL and PERM. */
/* > \endverbatim */
/* > */
/* > \param[in] PERM */
/* > \verbatim */
/* > PERM is INTEGER array, dimension ( LDGCOL, NLVL ). */
/* > On entry, PERM(*, I) records permutations done on the I-th */
/* > level of the computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVNUM */
/* > \verbatim */
/* > GIVNUM is REAL array, dimension ( LDU, 2 * NLVL ). */
/* > On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S- */
/* > values of Givens rotations performed on the I-th level on the */
/* > computation tree. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is REAL array, dimension ( N ). */
/* > On entry, if the I-th subproblem is not square, */
/* > C( I ) contains the C-value of a Givens rotation related to */
/* > the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is REAL array, dimension ( N ). */
/* > On entry, if the I-th subproblem is not square, */
/* > S( I ) contains the S-value of a Givens rotation related to */
/* > the right null space of the I-th subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension at least */
/* > MAX( (SMLSZ+1)*NRHS*3, N*(1+NRHS) + 2*NRHS ). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array. */
/* > The dimension must be at least 3 * N */
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
/* > \ingroup complexOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Ming Gu and Ren-Cang Li, Computer Science Division, University of */
/* > California at Berkeley, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* ===================================================================== */
/* Subroutine */
int clalsa_(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, complex *b, integer *ldb, complex *bx, integer *ldbx, real *u, integer *ldu, real *vt, integer *k, real *difl, real *difr, real *z__, real *poles, integer *givptr, integer *givcol, integer * ldgcol, integer *perm, real *givnum, real *c__, real *s, real *rwork, integer *iwork, integer *info)
{
    /* System generated locals */
    integer givcol_dim1, givcol_offset, perm_dim1, perm_offset, difl_dim1, difl_offset, difr_dim1, difr_offset, givnum_dim1, givnum_offset, poles_dim1, poles_offset, u_dim1, u_offset, vt_dim1, vt_offset, z_dim1, z_offset, b_dim1, b_offset, bx_dim1, bx_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    integer pow_ii(integer *, integer *);
    /* Local variables */
    integer i__, j, i1, ic, lf, nd, ll, nl, nr, im1, nlf, nrf, lvl, ndb1, nlp1, lvl2, nrp1, jcol, nlvl, sqre, jrow, jimag, jreal, inode, ndiml;
    extern /* Subroutine */
    int sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    integer ndimr;
    extern /* Subroutine */
    int ccopy_(integer *, complex *, integer *, complex *, integer *), clals0_(integer *, integer *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, integer *, integer *, integer *, integer *, real *, integer *, real *, real *, real *, real *, integer *, real *, real *, real *, integer *), xerbla_(char *, integer *), slasdt_(integer * , integer *, integer *, integer *, integer *, integer *, integer * );
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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    bx_dim1 = *ldbx;
    bx_offset = 1 + bx_dim1;
    bx -= bx_offset;
    givnum_dim1 = *ldu;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    poles_dim1 = *ldu;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    z_dim1 = *ldu;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    difr_dim1 = *ldu;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    difl_dim1 = *ldu;
    difl_offset = 1 + difl_dim1;
    difl -= difl_offset;
    vt_dim1 = *ldu;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --k;
    --givptr;
    perm_dim1 = *ldgcol;
    perm_offset = 1 + perm_dim1;
    perm -= perm_offset;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    --c__;
    --s;
    --rwork;
    --iwork;
    /* Function Body */
    *info = 0;
    if (*icompq < 0 || *icompq > 1)
    {
        *info = -1;
    }
    else if (*smlsiz < 3)
    {
        *info = -2;
    }
    else if (*n < *smlsiz)
    {
        *info = -3;
    }
    else if (*nrhs < 1)
    {
        *info = -4;
    }
    else if (*ldb < *n)
    {
        *info = -6;
    }
    else if (*ldbx < *n)
    {
        *info = -8;
    }
    else if (*ldu < *n)
    {
        *info = -10;
    }
    else if (*ldgcol < *n)
    {
        *info = -19;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CLALSA", &i__1);
        return 0;
    }
    /* Book-keeping and setting up the computation tree. */
    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    slasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], smlsiz);
    /* The following code applies back the left singular vector factors. */
    /* For applying back the right singular vector factors, go to 170. */
    if (*icompq == 1)
    {
        goto L170;
    }
    /* The nodes on the bottom level of the tree were solved */
    /* by SLASDQ. The corresponding left and right singular vector */
    /* matrices are in explicit form. First apply back the left */
    /* singular vector matrices. */
    ndb1 = (nd + 1) / 2;
    i__1 = nd;
    for (i__ = ndb1;
            i__ <= i__1;
            ++i__)
    {
        /* IC : center row of each node */
        /* NL : number of rows of left subproblem */
        /* NR : number of rows of right subproblem */
        /* NLF: starting row of the left subproblem */
        /* NRF: starting row of the right subproblem */
        i1 = i__ - 1;
        ic = iwork[inode + i1];
        nl = iwork[ndiml + i1];
        nr = iwork[ndimr + i1];
        nlf = ic - nl;
        nrf = ic + 1;
        /* Since B and BX are complex, the following call to SGEMM */
        /* is performed in two steps (real and imaginary parts). */
        /* CALL SGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU, */
        /* $ B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX ) */
        j = nl * *nrhs << 1;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nlf + nl - 1;
            for (jrow = nlf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++j;
                i__4 = jrow + jcol * b_dim1;
                rwork[j] = b[i__4].r;
                /* L10: */
            }
            /* L20: */
        }
        sgemm_("T", "N", &nl, nrhs, &nl, &c_b9, &u[nlf + u_dim1], ldu, &rwork[ (nl * *nrhs << 1) + 1], &nl, &c_b10, &rwork[1], &nl);
        j = nl * *nrhs << 1;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nlf + nl - 1;
            for (jrow = nlf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++j;
                rwork[j] = r_imag(&b[jrow + jcol * b_dim1]);
                /* L30: */
            }
            /* L40: */
        }
        sgemm_("T", "N", &nl, nrhs, &nl, &c_b9, &u[nlf + u_dim1], ldu, &rwork[ (nl * *nrhs << 1) + 1], &nl, &c_b10, &rwork[nl * *nrhs + 1], & nl);
        jreal = 0;
        jimag = nl * *nrhs;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nlf + nl - 1;
            for (jrow = nlf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++jreal;
                ++jimag;
                i__4 = jrow + jcol * bx_dim1;
                i__5 = jreal;
                i__6 = jimag;
                q__1.r = rwork[i__5];
                q__1.i = rwork[i__6]; // , expr subst
                bx[i__4].r = q__1.r;
                bx[i__4].i = q__1.i; // , expr subst
                /* L50: */
            }
            /* L60: */
        }
        /* Since B and BX are complex, the following call to SGEMM */
        /* is performed in two steps (real and imaginary parts). */
        /* CALL SGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU, */
        /* $ B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX ) */
        j = nr * *nrhs << 1;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nrf + nr - 1;
            for (jrow = nrf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++j;
                i__4 = jrow + jcol * b_dim1;
                rwork[j] = b[i__4].r;
                /* L70: */
            }
            /* L80: */
        }
        sgemm_("T", "N", &nr, nrhs, &nr, &c_b9, &u[nrf + u_dim1], ldu, &rwork[ (nr * *nrhs << 1) + 1], &nr, &c_b10, &rwork[1], &nr);
        j = nr * *nrhs << 1;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nrf + nr - 1;
            for (jrow = nrf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++j;
                rwork[j] = r_imag(&b[jrow + jcol * b_dim1]);
                /* L90: */
            }
            /* L100: */
        }
        sgemm_("T", "N", &nr, nrhs, &nr, &c_b9, &u[nrf + u_dim1], ldu, &rwork[ (nr * *nrhs << 1) + 1], &nr, &c_b10, &rwork[nr * *nrhs + 1], & nr);
        jreal = 0;
        jimag = nr * *nrhs;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nrf + nr - 1;
            for (jrow = nrf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++jreal;
                ++jimag;
                i__4 = jrow + jcol * bx_dim1;
                i__5 = jreal;
                i__6 = jimag;
                q__1.r = rwork[i__5];
                q__1.i = rwork[i__6]; // , expr subst
                bx[i__4].r = q__1.r;
                bx[i__4].i = q__1.i; // , expr subst
                /* L110: */
            }
            /* L120: */
        }
        /* L130: */
    }
    /* Next copy the rows of B that correspond to unchanged rows */
    /* in the bidiagonal matrix to BX. */
    i__1 = nd;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        ic = iwork[inode + i__ - 1];
        ccopy_(nrhs, &b[ic + b_dim1], ldb, &bx[ic + bx_dim1], ldbx);
        /* L140: */
    }
    /* Finally go through the left singular vector matrices of all */
    /* the other subproblems bottom-up on the tree. */
    j = pow_ii(&c__2, &nlvl);
    sqre = 0;
    for (lvl = nlvl;
            lvl >= 1;
            --lvl)
    {
        lvl2 = (lvl << 1) - 1;
        /* find the first node LF and last node LL on */
        /* the current level LVL */
        if (lvl == 1)
        {
            lf = 1;
            ll = 1;
        }
        else
        {
            i__1 = lvl - 1;
            lf = pow_ii(&c__2, &i__1);
            ll = (lf << 1) - 1;
        }
        i__1 = ll;
        for (i__ = lf;
                i__ <= i__1;
                ++i__)
        {
            im1 = i__ - 1;
            ic = iwork[inode + im1];
            nl = iwork[ndiml + im1];
            nr = iwork[ndimr + im1];
            nlf = ic - nl;
            nrf = ic + 1;
            --j;
            clals0_(icompq, &nl, &nr, &sqre, nrhs, &bx[nlf + bx_dim1], ldbx, & b[nlf + b_dim1], ldb, &perm[nlf + lvl * perm_dim1], & givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, & givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[ j], &s[j], &rwork[1], info);
            /* L150: */
        }
        /* L160: */
    }
    goto L330;
    /* ICOMPQ = 1: applying back the right singular vector factors. */
L170: /* First now go through the right singular vector matrices of all */
    /* the tree nodes top-down. */
    j = 0;
    i__1 = nlvl;
    for (lvl = 1;
            lvl <= i__1;
            ++lvl)
    {
        lvl2 = (lvl << 1) - 1;
        /* Find the first node LF and last node LL on */
        /* the current level LVL. */
        if (lvl == 1)
        {
            lf = 1;
            ll = 1;
        }
        else
        {
            i__2 = lvl - 1;
            lf = pow_ii(&c__2, &i__2);
            ll = (lf << 1) - 1;
        }
        i__2 = lf;
        for (i__ = ll;
                i__ >= i__2;
                --i__)
        {
            im1 = i__ - 1;
            ic = iwork[inode + im1];
            nl = iwork[ndiml + im1];
            nr = iwork[ndimr + im1];
            nlf = ic - nl;
            nrf = ic + 1;
            if (i__ == ll)
            {
                sqre = 0;
            }
            else
            {
                sqre = 1;
            }
            ++j;
            clals0_(icompq, &nl, &nr, &sqre, nrhs, &b[nlf + b_dim1], ldb, &bx[ nlf + bx_dim1], ldbx, &perm[nlf + lvl * perm_dim1], & givptr[j], &givcol[nlf + lvl2 * givcol_dim1], ldgcol, & givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1], &difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[ j], &s[j], &rwork[1], info);
            /* L180: */
        }
        /* L190: */
    }
    /* The nodes on the bottom level of the tree were solved */
    /* by SLASDQ. The corresponding right singular vector */
    /* matrices are in explicit form. Apply them back. */
    ndb1 = (nd + 1) / 2;
    i__1 = nd;
    for (i__ = ndb1;
            i__ <= i__1;
            ++i__)
    {
        i1 = i__ - 1;
        ic = iwork[inode + i1];
        nl = iwork[ndiml + i1];
        nr = iwork[ndimr + i1];
        nlp1 = nl + 1;
        if (i__ == nd)
        {
            nrp1 = nr;
        }
        else
        {
            nrp1 = nr + 1;
        }
        nlf = ic - nl;
        nrf = ic + 1;
        /* Since B and BX are complex, the following call to SGEMM is */
        /* performed in two steps (real and imaginary parts). */
        /* CALL SGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU, */
        /* $ B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX ) */
        j = nlp1 * *nrhs << 1;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nlf + nlp1 - 1;
            for (jrow = nlf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++j;
                i__4 = jrow + jcol * b_dim1;
                rwork[j] = b[i__4].r;
                /* L200: */
            }
            /* L210: */
        }
        sgemm_("T", "N", &nlp1, nrhs, &nlp1, &c_b9, &vt[nlf + vt_dim1], ldu, & rwork[(nlp1 * *nrhs << 1) + 1], &nlp1, &c_b10, &rwork[1], & nlp1);
        j = nlp1 * *nrhs << 1;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nlf + nlp1 - 1;
            for (jrow = nlf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++j;
                rwork[j] = r_imag(&b[jrow + jcol * b_dim1]);
                /* L220: */
            }
            /* L230: */
        }
        sgemm_("T", "N", &nlp1, nrhs, &nlp1, &c_b9, &vt[nlf + vt_dim1], ldu, & rwork[(nlp1 * *nrhs << 1) + 1], &nlp1, &c_b10, &rwork[nlp1 * * nrhs + 1], &nlp1);
        jreal = 0;
        jimag = nlp1 * *nrhs;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nlf + nlp1 - 1;
            for (jrow = nlf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++jreal;
                ++jimag;
                i__4 = jrow + jcol * bx_dim1;
                i__5 = jreal;
                i__6 = jimag;
                q__1.r = rwork[i__5];
                q__1.i = rwork[i__6]; // , expr subst
                bx[i__4].r = q__1.r;
                bx[i__4].i = q__1.i; // , expr subst
                /* L240: */
            }
            /* L250: */
        }
        /* Since B and BX are complex, the following call to SGEMM is */
        /* performed in two steps (real and imaginary parts). */
        /* CALL SGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU, */
        /* $ B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX ) */
        j = nrp1 * *nrhs << 1;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nrf + nrp1 - 1;
            for (jrow = nrf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++j;
                i__4 = jrow + jcol * b_dim1;
                rwork[j] = b[i__4].r;
                /* L260: */
            }
            /* L270: */
        }
        sgemm_("T", "N", &nrp1, nrhs, &nrp1, &c_b9, &vt[nrf + vt_dim1], ldu, & rwork[(nrp1 * *nrhs << 1) + 1], &nrp1, &c_b10, &rwork[1], & nrp1);
        j = nrp1 * *nrhs << 1;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nrf + nrp1 - 1;
            for (jrow = nrf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++j;
                rwork[j] = r_imag(&b[jrow + jcol * b_dim1]);
                /* L280: */
            }
            /* L290: */
        }
        sgemm_("T", "N", &nrp1, nrhs, &nrp1, &c_b9, &vt[nrf + vt_dim1], ldu, & rwork[(nrp1 * *nrhs << 1) + 1], &nrp1, &c_b10, &rwork[nrp1 * * nrhs + 1], &nrp1);
        jreal = 0;
        jimag = nrp1 * *nrhs;
        i__2 = *nrhs;
        for (jcol = 1;
                jcol <= i__2;
                ++jcol)
        {
            i__3 = nrf + nrp1 - 1;
            for (jrow = nrf;
                    jrow <= i__3;
                    ++jrow)
            {
                ++jreal;
                ++jimag;
                i__4 = jrow + jcol * bx_dim1;
                i__5 = jreal;
                i__6 = jimag;
                q__1.r = rwork[i__5];
                q__1.i = rwork[i__6]; // , expr subst
                bx[i__4].r = q__1.r;
                bx[i__4].i = q__1.i; // , expr subst
                /* L300: */
            }
            /* L310: */
        }
        /* L320: */
    }
L330:
    return 0;
    /* End of CLALSA */
}
/* clalsa_ */
