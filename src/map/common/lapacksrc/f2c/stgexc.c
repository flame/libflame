/* ../netlib/stgexc.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c__2 = 2;
/* > \brief \b STGEXC */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STGEXC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgexc. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgexc. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgexc. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, */
/* LDZ, IFST, ILST, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* LOGICAL WANTQ, WANTZ */
/* INTEGER IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGEXC reorders the generalized real Schur decomposition of a real */
/* > matrix pair (A,B) using an orthogonal equivalence transformation */
/* > */
/* > (A, B) = Q * (A, B) * Z**T, */
/* > */
/* > so that the diagonal block of (A, B) with row index IFST is moved */
/* > to row ILST. */
/* > */
/* > (A, B) must be in generalized real Schur canonical form (as returned */
/* > by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2 */
/* > diagonal blocks. B is upper triangular. */
/* > */
/* > Optionally, the matrices Q and Z of generalized Schur vectors are */
/* > updated. */
/* > */
/* > Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T */
/* > Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] WANTQ */
/* > \verbatim */
/* > WANTQ is LOGICAL */
/* > .TRUE. : update the left transformation matrix Q;
*/
/* > .FALSE.: do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* > WANTZ is LOGICAL */
/* > .TRUE. : update the right transformation matrix Z;
*/
/* > .FALSE.: do not update Z. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the matrix A in generalized real Schur canonical */
/* > form. */
/* > On exit, the updated matrix A, again in generalized */
/* > real Schur canonical form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,N) */
/* > On entry, the matrix B in generalized real Schur canonical */
/* > form (A,B). */
/* > On exit, the updated matrix B, again in generalized */
/* > real Schur canonical form (A,B). */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDZ,N) */
/* > On entry, if WANTQ = .TRUE., the orthogonal matrix Q. */
/* > On exit, the updated matrix Q. */
/* > If WANTQ = .FALSE., Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= 1. */
/* > If WANTQ = .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ,N) */
/* > On entry, if WANTZ = .TRUE., the orthogonal matrix Z. */
/* > On exit, the updated matrix Z. */
/* > If WANTZ = .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1. */
/* > If WANTZ = .TRUE., LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IFST */
/* > \verbatim */
/* > IFST is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] ILST */
/* > \verbatim */
/* > ILST is INTEGER */
/* > Specify the reordering of the diagonal blocks of (A, B). */
/* > The block with row index IFST is moved to row ILST, by a */
/* > sequence of swapping between adjacent blocks. */
/* > On exit, if IFST pointed on entry to the second row of */
/* > a 2-by-2 block, it is changed to point to the first row;
*/
/* > ILST always points to the first row of the block in its */
/* > final position (which may differ from its input value by */
/* > +1 or -1). 1 <= IFST, ILST <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > LWORK >= 1 when N <= 1, otherwise LWORK >= 4*N + 16. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > =0: successful exit. */
/* > <0: if INFO = -i, the i-th argument had an illegal value. */
/* > =1: The transformed matrix pair (A, B) would be too far */
/* > from generalized Schur form;
the problem is ill- */
/* > conditioned. (A, B) may have been partially reordered, */
/* > and ILST points to the first row of the current */
/* > position of the block being moved. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realGEcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* > Umea University, S-901 87 Umea, Sweden. */
/* > \par References: */
/* ================ */
/* > */
/* > \verbatim */
/* > */
/* > [1] B. Kagstrom;
A Direct Method for Reordering Eigenvalues in the */
/* > Generalized Real Schur Form of a Regular Matrix Pair (A, B), in */
/* > M.S. Moonen et al (eds), Linear Algebra for Large Scale and */
/* > Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int stgexc_(logical *wantq, logical *wantz, integer *n, real *a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real * z__, integer *ldz, integer *ifst, integer *ilst, real *work, integer * lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, i__1;
    /* Local variables */
    integer nbf, nbl, here, lwmin;
    extern /* Subroutine */
    int stgex2_(logical *, logical *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *, integer *, integer *, real *, integer *, integer *), xerbla_(char *, integer *);
    integer nbnext;
    logical lquery;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    /* Decode and test input arguments. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if (*n < 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -7;
    }
    else if (*ldq < 1 || *wantq && *ldq < max(1,*n))
    {
        *info = -9;
    }
    else if (*ldz < 1 || *wantz && *ldz < max(1,*n))
    {
        *info = -11;
    }
    else if (*ifst < 1 || *ifst > *n)
    {
        *info = -12;
    }
    else if (*ilst < 1 || *ilst > *n)
    {
        *info = -13;
    }
    if (*info == 0)
    {
        if (*n <= 1)
        {
            lwmin = 1;
        }
        else
        {
            lwmin = (*n << 2) + 16;
        }
        work[1] = (real) lwmin;
        if (*lwork < lwmin && ! lquery)
        {
            *info = -15;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STGEXC", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*n <= 1)
    {
        return 0;
    }
    /* Determine the first row of the specified block and find out */
    /* if it is 1-by-1 or 2-by-2. */
    if (*ifst > 1)
    {
        if (a[*ifst + (*ifst - 1) * a_dim1] != 0.f)
        {
            --(*ifst);
        }
    }
    nbf = 1;
    if (*ifst < *n)
    {
        if (a[*ifst + 1 + *ifst * a_dim1] != 0.f)
        {
            nbf = 2;
        }
    }
    /* Determine the first row of the final block */
    /* and find out if it is 1-by-1 or 2-by-2. */
    if (*ilst > 1)
    {
        if (a[*ilst + (*ilst - 1) * a_dim1] != 0.f)
        {
            --(*ilst);
        }
    }
    nbl = 1;
    if (*ilst < *n)
    {
        if (a[*ilst + 1 + *ilst * a_dim1] != 0.f)
        {
            nbl = 2;
        }
    }
    if (*ifst == *ilst)
    {
        return 0;
    }
    if (*ifst < *ilst)
    {
        /* Update ILST. */
        if (nbf == 2 && nbl == 1)
        {
            --(*ilst);
        }
        if (nbf == 1 && nbl == 2)
        {
            ++(*ilst);
        }
        here = *ifst;
L10: /* Swap with next one below. */
        if (nbf == 1 || nbf == 2)
        {
            /* Current block either 1-by-1 or 2-by-2. */
            nbnext = 1;
            if (here + nbf + 1 <= *n)
            {
                if (a[here + nbf + 1 + (here + nbf) * a_dim1] != 0.f)
                {
                    nbnext = 2;
                }
            }
            stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[ q_offset], ldq, &z__[z_offset], ldz, &here, &nbf, &nbnext, &work[1], lwork, info);
            if (*info != 0)
            {
                *ilst = here;
                return 0;
            }
            here += nbnext;
            /* Test if 2-by-2 block breaks into two 1-by-1 blocks. */
            if (nbf == 2)
            {
                if (a[here + 1 + here * a_dim1] == 0.f)
                {
                    nbf = 3;
                }
            }
        }
        else
        {
            /* Current block consists of two 1-by-1 blocks, each of which */
            /* must be swapped individually. */
            nbnext = 1;
            if (here + 3 <= *n)
            {
                if (a[here + 3 + (here + 2) * a_dim1] != 0.f)
                {
                    nbnext = 2;
                }
            }
            i__1 = here + 1;
            stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[ q_offset], ldq, &z__[z_offset], ldz, &i__1, &c__1, & nbnext, &work[1], lwork, info);
            if (*info != 0)
            {
                *ilst = here;
                return 0;
            }
            if (nbnext == 1)
            {
                /* Swap two 1-by-1 blocks. */
                stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &here, &c__1, &c__1, &work[1], lwork, info);
                if (*info != 0)
                {
                    *ilst = here;
                    return 0;
                }
                ++here;
            }
            else
            {
                /* Recompute NBNEXT in case of 2-by-2 split. */
                if (a[here + 2 + (here + 1) * a_dim1] == 0.f)
                {
                    nbnext = 1;
                }
                if (nbnext == 2)
                {
                    /* 2-by-2 block did not split. */
                    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, & here, &c__1, &nbnext, &work[1], lwork, info);
                    if (*info != 0)
                    {
                        *ilst = here;
                        return 0;
                    }
                    here += 2;
                }
                else
                {
                    /* 2-by-2 block did split. */
                    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, & here, &c__1, &c__1, &work[1], lwork, info);
                    if (*info != 0)
                    {
                        *ilst = here;
                        return 0;
                    }
                    ++here;
                    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, & here, &c__1, &c__1, &work[1], lwork, info);
                    if (*info != 0)
                    {
                        *ilst = here;
                        return 0;
                    }
                    ++here;
                }
            }
        }
        if (here < *ilst)
        {
            goto L10;
        }
    }
    else
    {
        here = *ifst;
L20: /* Swap with next one below. */
        if (nbf == 1 || nbf == 2)
        {
            /* Current block either 1-by-1 or 2-by-2. */
            nbnext = 1;
            if (here >= 3)
            {
                if (a[here - 1 + (here - 2) * a_dim1] != 0.f)
                {
                    nbnext = 2;
                }
            }
            i__1 = here - nbnext;
            stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[ q_offset], ldq, &z__[z_offset], ldz, &i__1, &nbnext, &nbf, &work[1], lwork, info);
            if (*info != 0)
            {
                *ilst = here;
                return 0;
            }
            here -= nbnext;
            /* Test if 2-by-2 block breaks into two 1-by-1 blocks. */
            if (nbf == 2)
            {
                if (a[here + 1 + here * a_dim1] == 0.f)
                {
                    nbf = 3;
                }
            }
        }
        else
        {
            /* Current block consists of two 1-by-1 blocks, each of which */
            /* must be swapped individually. */
            nbnext = 1;
            if (here >= 3)
            {
                if (a[here - 1 + (here - 2) * a_dim1] != 0.f)
                {
                    nbnext = 2;
                }
            }
            i__1 = here - nbnext;
            stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[ q_offset], ldq, &z__[z_offset], ldz, &i__1, &nbnext, & c__1, &work[1], lwork, info);
            if (*info != 0)
            {
                *ilst = here;
                return 0;
            }
            if (nbnext == 1)
            {
                /* Swap two 1-by-1 blocks. */
                stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &here, & nbnext, &c__1, &work[1], lwork, info);
                if (*info != 0)
                {
                    *ilst = here;
                    return 0;
                }
                --here;
            }
            else
            {
                /* Recompute NBNEXT in case of 2-by-2 split. */
                if (a[here + (here - 1) * a_dim1] == 0.f)
                {
                    nbnext = 1;
                }
                if (nbnext == 2)
                {
                    /* 2-by-2 block did not split. */
                    i__1 = here - 1;
                    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, & i__1, &c__2, &c__1, &work[1], lwork, info);
                    if (*info != 0)
                    {
                        *ilst = here;
                        return 0;
                    }
                    here += -2;
                }
                else
                {
                    /* 2-by-2 block did split. */
                    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, & here, &c__1, &c__1, &work[1], lwork, info);
                    if (*info != 0)
                    {
                        *ilst = here;
                        return 0;
                    }
                    --here;
                    stgex2_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, & here, &c__1, &c__1, &work[1], lwork, info);
                    if (*info != 0)
                    {
                        *ilst = here;
                        return 0;
                    }
                    --here;
                }
            }
        }
        if (here > *ilst)
        {
            goto L20;
        }
    }
    *ilst = here;
    work[1] = (real) lwmin;
    return 0;
    /* End of STGEXC */
}
/* stgexc_ */
