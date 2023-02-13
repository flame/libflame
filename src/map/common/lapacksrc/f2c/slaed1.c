/* ../netlib/slaed1.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b SLAED1 used by sstedc. Computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix. Used when the original matrix is tridiagonal. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAED1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed1. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed1. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed1. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER CUTPNT, INFO, LDQ, N */
/* REAL RHO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER INDXQ( * ), IWORK( * ) */
/* REAL D( * ), Q( LDQ, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAED1 computes the updated eigensystem of a diagonal */
/* > matrix after modification by a rank-one symmetric matrix. This */
/* > routine is used only for the eigenproblem which requires all */
/* > eigenvalues and eigenvectors of a tridiagonal matrix. SLAED7 handles */
/* > the case in which eigenvalues only or eigenvalues and eigenvectors */
/* > of a full symmetric matrix (which was reduced to tridiagonal form) */
/* > are desired. */
/* > */
/* > T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out) */
/* > */
/* > where Z = Q**T*u, u is a vector of length N with ones in the */
/* > CUTPNT and CUTPNT + 1 th elements and zeros elsewhere. */
/* > */
/* > The eigenvectors of the original matrix are stored in Q, and the */
/* > eigenvalues are in D. The algorithm consists of three stages: */
/* > */
/* > The first stage consists of deflating the size of the problem */
/* > when there are multiple eigenvalues or if there is a zero in */
/* > the Z vector. For each such occurence the dimension of the */
/* > secular equation problem is reduced by one. This stage is */
/* > performed by the routine SLAED2. */
/* > */
/* > The second stage consists of calculating the updated */
/* > eigenvalues. This is done by finding the roots of the secular */
/* > equation via the routine SLAED4 (as called by SLAED3). */
/* > This routine also calculates the eigenvectors of the current */
/* > problem. */
/* > */
/* > The final stage consists of computing the updated eigenvectors */
/* > directly using the updated eigenvalues. The eigenvectors for */
/* > the current problem are multiplied with the eigenvectors from */
/* > the overall problem. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The dimension of the symmetric tridiagonal matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, the eigenvalues of the rank-1-perturbed matrix. */
/* > On exit, the eigenvalues of the repaired matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDQ,N) */
/* > On entry, the eigenvectors of the rank-1-perturbed matrix. */
/* > On exit, the eigenvectors of the repaired tridiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] INDXQ */
/* > \verbatim */
/* > INDXQ is INTEGER array, dimension (N) */
/* > On entry, the permutation which separately sorts the two */
/* > subproblems in D into ascending order. */
/* > On exit, the permutation which will reintegrate the */
/* > subproblems back into sorted order, */
/* > i.e. D( INDXQ( I = 1, N ) ) will be in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* > RHO is REAL */
/* > The subdiagonal entry used to create the rank-1 modification. */
/* > \endverbatim */
/* > */
/* > \param[in] CUTPNT */
/* > \verbatim */
/* > CUTPNT is INTEGER */
/* > The location of the last eigenvalue in the leading sub-matrix. */
/* > min(1,N) <= CUTPNT <= N/2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (4*N + N**2) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = 1, an eigenvalue did not converge */
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
/* > at Berkeley, USA \n */
/* > Modified by Francoise Tisseur, University of Tennessee */
/* > */
/* ===================================================================== */
/* Subroutine */
int slaed1_(integer *n, real *d__, real *q, integer *ldq, integer *indxq, real *rho, integer *cutpnt, real *work, integer * iwork, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2;
    /* Local variables */
    integer i__, k, n1, n2, is, iw, iz, iq2, cpp1, indx, indxc, indxp;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *), slaed2_(integer *, integer *, integer *, real *, real *, integer *, integer *, real *, real *, real *, real *, real *, integer *, integer *, integer *, integer *, integer *), slaed3_( integer *, integer *, integer *, real *, real *, integer *, real * , real *, real *, integer *, integer *, real *, real *, integer *) ;
    integer idlmda;
    extern /* Subroutine */
    int xerbla_(char *, integer *), slamrg_( integer *, integer *, real *, integer *, integer *, integer *);
    integer coltyp;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
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
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --indxq;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*ldq < max(1,*n))
    {
        *info = -4;
    }
    else /* if(complicated condition) */
    {
        /* Computing MIN */
        i__1 = 1;
        i__2 = *n / 2; // , expr subst
        if (min(i__1,i__2) > *cutpnt || *n / 2 < *cutpnt)
        {
            *info = -7;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLAED1", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* The following values are integer pointers which indicate */
    /* the portion of the workspace */
    /* used by a particular array in SLAED2 and SLAED3. */
    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq2 = iw + *n;
    indx = 1;
    indxc = indx + *n;
    coltyp = indxc + *n;
    indxp = coltyp + *n;
    /* Form the z-vector which consists of the last row of Q_1 and the */
    /* first row of Q_2. */
    scopy_(cutpnt, &q[*cutpnt + q_dim1], ldq, &work[iz], &c__1);
    cpp1 = *cutpnt + 1;
    i__1 = *n - *cutpnt;
    scopy_(&i__1, &q[cpp1 + cpp1 * q_dim1], ldq, &work[iz + *cutpnt], &c__1);
    /* Deflate eigenvalues. */
    slaed2_(&k, n, cutpnt, &d__[1], &q[q_offset], ldq, &indxq[1], rho, &work[ iz], &work[idlmda], &work[iw], &work[iq2], &iwork[indx], &iwork[ indxc], &iwork[indxp], &iwork[coltyp], info);
    if (*info != 0)
    {
        goto L20;
    }
    /* Solve Secular Equation. */
    if (k != 0)
    {
        is = (iwork[coltyp] + iwork[coltyp + 1]) * *cutpnt + (iwork[coltyp + 1] + iwork[coltyp + 2]) * (*n - *cutpnt) + iq2;
        slaed3_(&k, n, cutpnt, &d__[1], &q[q_offset], ldq, rho, &work[idlmda], &work[iq2], &iwork[indxc], &iwork[coltyp], &work[iw], &work[ is], info);
        if (*info != 0)
        {
            goto L20;
        }
        /* Prepare the INDXQ sorting permutation. */
        n1 = k;
        n2 = *n - k;
        slamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &indxq[1]);
    }
    else
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            indxq[i__] = i__;
            /* L10: */
        }
    }
L20:
    return 0;
    /* End of SLAED1 */
}
/* slaed1_ */
