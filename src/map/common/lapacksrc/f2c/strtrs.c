/* ../netlib/strtrs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b12 = 1.f;
/* > \brief \b STRTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STRTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strtrs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strtrs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strtrs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, TRANS, UPLO */
/* INTEGER INFO, LDA, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRTRS solves a triangular system of the form */
/* > */
/* > A * X = B or A**T * X = B, */
/* > */
/* > where A is a triangular matrix of order N, and B is an N-by-NRHS */
/* > matrix. A check is made to verify that A is nonsingular. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': A is upper triangular;
*/
/* > = 'L': A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations: */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* > DIAG is CHARACTER*1 */
/* > = 'N': A is non-unit triangular;
*/
/* > = 'U': A is unit triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The triangular matrix A. If UPLO = 'U', the leading N-by-N */
/* > upper triangular part of the array A contains the upper */
/* > triangular matrix, and the strictly lower triangular part of */
/* > A is not referenced. If UPLO = 'L', the leading N-by-N lower */
/* > triangular part of the array A contains the lower triangular */
/* > matrix, and the strictly upper triangular part of A is not */
/* > referenced. If DIAG = 'U', the diagonal elements of A are */
/* > also not referenced and are assumed to be 1. */
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
/* > B is REAL array, dimension (LDB,NRHS) */
/* > On entry, the right hand side matrix B. */
/* > On exit, if INFO = 0, the solution matrix X. */
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
/* > > 0: if INFO = i, the i-th diagonal element of A is zero, */
/* > indicating that the matrix is singular and the solutions */
/* > X have not been computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int strtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, real *a, integer *lda, real *b, integer *ldb, integer * info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int strsm_(char *, char *, char *, char *, integer *, integer *, real *, real *, integer *, real *, integer * ), xerbla_(char *, integer *);
    logical nounit;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    nounit = lsame_(diag, "N");
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! lsame_(trans, "C"))
    {
        *info = -2;
    }
    else if (! nounit && ! lsame_(diag, "U"))
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*nrhs < 0)
    {
        *info = -5;
    }
    else if (*lda < max(1,*n))
    {
        *info = -7;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -9;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STRTRS", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Check for singularity. */
    if (nounit)
    {
        i__1 = *n;
        for (*info = 1;
                *info <= i__1;
                ++(*info))
        {
            if (a[*info + *info * a_dim1] == 0.f)
            {
                return 0;
            }
            /* L10: */
        }
    }
    *info = 0;
    /* Solve A * x = b or A**T * x = b. */
    strsm_("Left", uplo, trans, diag, n, nrhs, &c_b12, &a[a_offset], lda, &b[ b_offset], ldb);
    return 0;
    /* End of STRTRS */
}
/* strtrs_ */
