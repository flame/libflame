/* ../netlib/sspgv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SSPGST */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSPGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspgv.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspgv.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspgv.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, ITYPE, LDZ, N */
/* .. */
/* .. Array Arguments .. */
/* REAL AP( * ), BP( * ), W( * ), WORK( * ), */
/* $ Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPGV computes all the eigenvalues and, optionally, the eigenvectors */
/* > of a real generalized symmetric-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x. */
/* > Here A and B are assumed to be symmetric, stored in packed format, */
/* > and B is also positive definite. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ITYPE */
/* > \verbatim */
/* > ITYPE is INTEGER */
/* > Specifies the problem type to be solved: */
/* > = 1: A*x = (lambda)*B*x */
/* > = 2: A*B*x = (lambda)*x */
/* > = 3: B*A*x = (lambda)*x */
/* > \endverbatim */
/* > */
/* > \param[in] JOBZ */
/* > \verbatim */
/* > JOBZ is CHARACTER*1 */
/* > = 'N': Compute eigenvalues only;
*/
/* > = 'V': Compute eigenvalues and eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangles of A and B are stored;
*/
/* > = 'L': Lower triangles of A and B are stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is REAL array, dimension */
/* > (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the symmetric matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* > On exit, the contents of AP are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] BP */
/* > \verbatim */
/* > BP is REAL array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the symmetric matrix */
/* > B, packed columnwise in a linear array. The j-th column of B */
/* > is stored in the array BP as follows: */
/* > if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n. */
/* > */
/* > On exit, the triangular factor U or L from the Cholesky */
/* > factorization B = U**T*U or B = L*L**T, in the same storage */
/* > format as B. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ, N) */
/* > If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of */
/* > eigenvectors. The eigenvectors are normalized as follows: */
/* > if ITYPE = 1 or 2, Z**T*B*Z = I;
*/
/* > if ITYPE = 3, Z**T*inv(B)*Z = I. */
/* > If JOBZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1, and if */
/* > JOBZ = 'V', LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: SPPTRF or SSPEV returned an error code: */
/* > <= N: if INFO = i, SSPEV failed to converge;
*/
/* > i off-diagonal elements of an intermediate */
/* > tridiagonal form did not converge to zero. */
/* > > N: if INFO = n + i, for 1 <= i <= n, then the leading */
/* > minor of order i of B is not positive definite. */
/* > The factorization of B could not be completed and */
/* > no eigenvalues or eigenvectors were computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realOTHEReigen */
/* ===================================================================== */
/* Subroutine */
int sspgv_(integer *itype, char *jobz, char *uplo, integer * n, real *ap, real *bp, real *w, real *z__, integer *ldz, real *work, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    /* Local variables */
    integer j, neig;
    extern logical lsame_(char *, char *);
    char trans[1];
    logical upper;
    extern /* Subroutine */
    int sspev_(char *, char *, integer *, real *, real *, real *, integer *, real *, integer *);
    logical wantz;
    extern /* Subroutine */
    int stpmv_(char *, char *, char *, integer *, real *, real *, integer *), stpsv_(char *, char *, char *, integer *, real *, real *, integer *), xerbla_(char *, integer *), spptrf_(char *, integer *, real *, integer *), sspgst_(integer *, char *, integer *, real *, real *, integer *);
    /* -- LAPACK driver routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --ap;
    --bp;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    /* Function Body */
    wantz = lsame_(jobz, "V");
    upper = lsame_(uplo, "U");
    *info = 0;
    if (*itype < 1 || *itype > 3)
    {
        *info = -1;
    }
    else if (! (wantz || lsame_(jobz, "N")))
    {
        *info = -2;
    }
    else if (! (upper || lsame_(uplo, "L")))
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -9;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSPGV ", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Form a Cholesky factorization of B. */
    spptrf_(uplo, n, &bp[1], info);
    if (*info != 0)
    {
        *info = *n + *info;
        return 0;
    }
    /* Transform problem to standard eigenvalue problem and solve. */
    sspgst_(itype, uplo, n, &ap[1], &bp[1], info);
    sspev_(jobz, uplo, n, &ap[1], &w[1], &z__[z_offset], ldz, &work[1], info);
    if (wantz)
    {
        /* Backtransform eigenvectors to the original problem. */
        neig = *n;
        if (*info > 0)
        {
            neig = *info - 1;
        }
        if (*itype == 1 || *itype == 2)
        {
            /* For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            */
            /* backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y */
            if (upper)
            {
                *(unsigned char *)trans = 'N';
            }
            else
            {
                *(unsigned char *)trans = 'T';
            }
            i__1 = neig;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                stpsv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 1], &c__1);
                /* L10: */
            }
        }
        else if (*itype == 3)
        {
            /* For B*A*x=(lambda)*x;
            */
            /* backtransform eigenvectors: x = L*y or U**T*y */
            if (upper)
            {
                *(unsigned char *)trans = 'T';
            }
            else
            {
                *(unsigned char *)trans = 'N';
            }
            i__1 = neig;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                stpmv_(uplo, trans, "Non-unit", n, &bp[1], &z__[j * z_dim1 + 1], &c__1);
                /* L20: */
            }
        }
    }
    return 0;
    /* End of SSPGV */
}
/* sspgv_ */
