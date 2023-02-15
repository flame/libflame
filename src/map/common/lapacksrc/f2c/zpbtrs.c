/* ../netlib/zpbtrs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZPBTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZPBTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpbtrs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpbtrs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpbtrs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, KD, LDAB, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 AB( LDAB, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPBTRS solves a system of linear equations A*X = B with a Hermitian */
/* > positive definite band matrix A using the Cholesky factorization */
/* > A = U**H *U or A = L*L**H computed by ZPBTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangular factor stored in AB;
*/
/* > = 'L': Lower triangular factor stored in AB. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is COMPLEX*16 array, dimension (LDAB,N) */
/* > The triangular factor U or L from the Cholesky factorization */
/* > A = U**H *U or A = L*L**H of the band matrix A, stored in the */
/* > first KD+1 rows of the array. The j-th column of U or L is */
/* > stored in the j-th column of the array AB as follows: */
/* > if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
*/
/* > if UPLO ='L', AB(1+i-j,j) = L(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the right hand side matrix B. */
/* > On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complex16OTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int zpbtrs_(char *uplo, integer *n, integer *kd, integer * nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer * ldb, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    integer j;
    extern logical lsame_(char *, char *);
    logical upper;
    extern /* Subroutine */
    int ztbsv_(char *, char *, char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(char *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*kd < 0)
    {
        *info = -3;
    }
    else if (*nrhs < 0)
    {
        *info = -4;
    }
    else if (*ldab < *kd + 1)
    {
        *info = -6;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -8;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZPBTRS", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *nrhs == 0)
    {
        return 0;
    }
    if (upper)
    {
        /* Solve A*X = B where A = U**H *U. */
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Solve U**H *X = B, overwriting B with X. */
            ztbsv_("Upper", "Conjugate transpose", "Non-unit", n, kd, &ab[ ab_offset], ldab, &b[j * b_dim1 + 1], &c__1);
            /* Solve U*X = B, overwriting B with X. */
            ztbsv_("Upper", "No transpose", "Non-unit", n, kd, &ab[ab_offset], ldab, &b[j * b_dim1 + 1], &c__1);
            /* L10: */
        }
    }
    else
    {
        /* Solve A*X = B where A = L*L**H. */
        i__1 = *nrhs;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Solve L*X = B, overwriting B with X. */
            ztbsv_("Lower", "No transpose", "Non-unit", n, kd, &ab[ab_offset], ldab, &b[j * b_dim1 + 1], &c__1);
            /* Solve L**H *X = B, overwriting B with X. */
            ztbsv_("Lower", "Conjugate transpose", "Non-unit", n, kd, &ab[ ab_offset], ldab, &b[j * b_dim1 + 1], &c__1);
            /* L20: */
        }
    }
    return 0;
    /* End of ZPBTRS */
}
/* zpbtrs_ */
