/* ../netlib/chesv_rook.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b CHESV_ROOK computes the solution to a system of linear equations A * X = B for HE matrices usin g the bounded Bunch-Kaufman ("rook") diagonal pivoting method */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHESV_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chesv_r ook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chesv_r ook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chesv_r ook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHESV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, */
/* LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, LWORK, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHESV_ROOK computes the solution to a complex system of linear equations */
/* > A * X = B, */
/* > where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS */
/* > matrices. */
/* > */
/* > The bounded Bunch-Kaufman ("rook") diagonal pivoting method is used */
/* > to factor A as */
/* > A = U * D * U**T, if UPLO = 'U', or */
/* > A = L * D * L**T, if UPLO = 'L', */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is Hermitian and block diagonal with */
/* > 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > CHETRF_ROOK is called to compute the factorization of a complex */
/* > Hermition matrix A using the bounded Bunch-Kaufman ("rook") diagonal */
/* > pivoting method. */
/* > */
/* > The factored form of A is then used to solve the system */
/* > of equations A * X = B by calling CHETRS_ROOK (uses BLAS 2). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
*/
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, if INFO = 0, the block diagonal matrix D and the */
/* > multipliers used to obtain the factor U or L from the */
/* > factorization A = U*D*U**H or A = L*D*L**H as computed by */
/* > CHETRF_ROOK. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D. */
/* > */
/* > If UPLO = 'U': */
/* > Only the last KB elements of IPIV are set. */
/* > */
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* > interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* > If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and */
/* > columns k and -IPIV(k) were interchanged and rows and */
/* > columns k-1 and -IPIV(k-1) were inerchaged, */
/* > D(k-1:k,k-1:k) is a 2-by-2 diagonal block. */
/* > */
/* > If UPLO = 'L': */
/* > Only the first KB elements of IPIV are set. */
/* > */
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) */
/* > were interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* > If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and */
/* > columns k and -IPIV(k) were interchanged and rows and */
/* > columns k+1 and -IPIV(k+1) were inerchaged, */
/* > D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the N-by-NRHS right hand side matrix B. */
/* > On exit, if INFO = 0, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The length of WORK. LWORK >= 1, and for best performance */
/* > LWORK >= max(1,N*NB), where NB is the optimal blocksize for */
/* > CHETRF_ROOK. */
/* > for LWORK < N, TRS will be done with Level BLAS 2 */
/* > for LWORK >= N, TRS will be done with Level BLAS 3 */
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
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, D(i,i) is exactly zero. The factorization */
/* > has been completed, but the block diagonal matrix D is */
/* > exactly singular, so the solution could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup complexHEsolve */
/* > */
/* > \verbatim */
/* > */
/* > November 2013, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* > School of Mathematics, */
/* > University of Manchester */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int chesv_rook_(char *uplo, integer *n, integer *nrhs, complex *a, integer *lda, integer *ipiv, complex *b, integer *ldb, complex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
    char buffer[256]; 
#if FLA_ENABLE_ILP64 
    snprintf(buffer, 256,"chesv_rook inputs: uplo %c, n %lld, nrhs %lld, lda %lld, ldb %lld, lwork %lld",*uplo, *n, *nrhs, *lda, *ldb, *lwork);
#else 
    snprintf(buffer, 256,"chesv_rook inputs: uplo %c, n %d, nrhs %d, lda %d, ldb %d, lwork %d",*uplo, *n, *nrhs, *lda, *ldb, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    integer nb;
    extern /* Subroutine */
    int chetrf_rook_(char *, integer *, complex *, integer *, integer *, complex *, integer *, integer *), chetrs_rook_(char *, integer *, integer *, complex *, integer *, integer *, complex *, integer *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    /* -- LAPACK driver routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*nrhs < 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -8;
    }
    else if (*lwork < 1 && ! lquery)
    {
        *info = -10;
    }
    if (*info == 0)
    {
        if (*n == 0)
        {
            lwkopt = 1;
        }
        else
        {
            nb = ilaenv_(&c__1, "CHETRF_ROOK", uplo, n, &c_n1, &c_n1, &c_n1);
            lwkopt = *n * nb;
        }
        work[1].r = (real) lwkopt;
        work[1].i = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHESV_ROOK ", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Compute the factorization A = U*D*U**H or A = L*D*L**H. */
    chetrf_rook_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], lwork, info);
    if (*info == 0)
    {
        /* Solve the system A*X = B, overwriting B with X. */
        /* Solve with TRS ( Use Level BLAS 2) */
        chetrs_rook_(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset] , ldb, info);
    }
    work[1].r = (real) lwkopt;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CHESV_ROOK */
}
/* chesv_rook__ */
