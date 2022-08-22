/* ../netlib/v3.9.0/zsytrs_3.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    1.,0.
}
;
/* > \brief \b ZSYTRS_3 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZSYTRS_3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs_ 3.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs_ 3.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs_ 3.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), E( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > ZSYTRS_3 solves a system of linear equations A * X = B with a complex */
/* > symmetric matrix A using the factorization computed */
/* > by ZSYTRF_RK or ZSYTRF_BK: */
/* > */
/* > A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**T (or L**T) is the transpose of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is symmetric and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This algorithm is using Level 3 BLAS. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are */
/* > stored as an upper or lower triangular matrix: */
/* > = 'U': Upper triangular, form is A = P*U*D*(U**T)*(P**T);
*/
/* > = 'L': Lower triangular, form is A = P*L*D*(L**T)*(P**T). */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > Diagonal of the block diagonal matrix D and factors U or L */
/* > as computed by ZSYTRF_RK and ZSYTRF_BK: */
/* > a) ONLY diagonal elements of the symmetric block diagonal */
/* > matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*/
/* > (superdiagonal (or subdiagonal) elements of D */
/* > should be provided on entry in array E), and */
/* > b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* > If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is COMPLEX*16 array, dimension (N) */
/* > On entry, contains the superdiagonal (or subdiagonal) */
/* > elements of the symmetric block diagonal matrix D */
/* > with 1-by-1 or 2-by-2 diagonal blocks, where */
/* > If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*/
/* > If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. */
/* > */
/* > NOTE: For 1-by-1 diagonal block D(k), where */
/* > 1 <= k <= N, the element E(k) is not referenced in both */
/* > UPLO = 'U' or UPLO = 'L' cases. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by ZSYTRF_RK or ZSYTRF_BK. */
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
/* > \date June 2017 */
/* > \ingroup complex16SYcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > June 2017, Igor Kozachenko, */
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
int zsytrs_3_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *e, integer *ipiv, doublecomplex *b, integer *ldb, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zsytrs_3 inputs: uplo %c, n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "", *uplo, *n, *nrhs, *lda, *ldb);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    /* Local variables */
    integer i__, j, k;
    doublecomplex ak, bk;
    integer kp;
    doublecomplex akm1, bkm1, akm1k;
    extern logical lsame_(char *, char *);
    doublecomplex denom;
    extern /* Subroutine */
    int zscal_(integer *, doublecomplex *, doublecomplex *, integer *);
    logical upper;
    extern /* Subroutine */
    int zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(char *, integer *);
    /* -- LAPACK computational routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
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
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --ipiv;
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
        *info = -9;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZSYTRS_3", &i__1);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *nrhs == 0)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (upper)
    {
        /* Begin Upper */
        /* Solve A*X = B, where A = U*D*U**T. */
        /* P**T * B */
        /* Interchange rows K and IPIV(K) of matrix B in the same order */
        /* that the formation order of IPIV(I) vector for Upper case. */
        /* (We can do the simple loop over IPIV with decrement -1, */
        /* since the ABS value of IPIV(I) represents the row index */
        /* of the interchange with row i in both 1x1 and 2x2 pivot cases) */
        for (k = *n;
                k >= 1;
                --k)
        {
            kp = (i__1 = ipiv[k], f2c_dabs(i__1));
            if (kp != k)
            {
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
        }
        /* Compute (U \P**T * B) -> B [ (U \P**T * B) ] */
        ztrsm_("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[ b_offset], ldb);
        /* Compute D \ B -> B [ D \ (U \P**T * B) ] */
        i__ = *n;
        while(i__ >= 1)
        {
            if (ipiv[i__] > 0)
            {
                z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
                zscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
            }
            else if (i__ > 1)
            {
                i__1 = i__;
                akm1k.r = e[i__1].r;
                akm1k.i = e[i__1].i; // , expr subst
                z_div(&z__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
                akm1.r = z__1.r;
                akm1.i = z__1.i; // , expr subst
                z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
                ak.r = z__1.r;
                ak.i = z__1.i; // , expr subst
                z__2.r = akm1.r * ak.r - akm1.i * ak.i;
                z__2.i = akm1.r * ak.i + akm1.i * ak.r; // , expr subst
                z__1.r = z__2.r - 1.;
                z__1.i = z__2.i - 0.; // , expr subst
                denom.r = z__1.r;
                denom.i = z__1.i; // , expr subst
                i__1 = *nrhs;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    z_div(&z__1, &b[i__ - 1 + j * b_dim1], &akm1k);
                    bkm1.r = z__1.r;
                    bkm1.i = z__1.i; // , expr subst
                    z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
                    bk.r = z__1.r;
                    bk.i = z__1.i; // , expr subst
                    i__2 = i__ - 1 + j * b_dim1;
                    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                    z__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                    z__2.r = z__3.r - bk.r;
                    z__2.i = z__3.i - bk.i; // , expr subst
                    z_div(&z__1, &z__2, &denom);
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                    i__2 = i__ + j * b_dim1;
                    z__3.r = akm1.r * bk.r - akm1.i * bk.i;
                    z__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                    z__2.r = z__3.r - bkm1.r;
                    z__2.i = z__3.i - bkm1.i; // , expr subst
                    z_div(&z__1, &z__2, &denom);
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                }
                --i__;
            }
            --i__;
        }
        /* Compute (U**T \ B) -> B [ U**T \ (D \ (U \P**T * B) ) ] */
        ztrsm_("L", "U", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[ b_offset], ldb);
        /* P * B [ P * (U**T \ (D \ (U \P**T * B) )) ] */
        /* Interchange rows K and IPIV(K) of matrix B in reverse order */
        /* from the formation order of IPIV(I) vector for Upper case. */
        /* (We can do the simple loop over IPIV with increment 1, */
        /* since the ABS value of IPIV(I) represents the row index */
        /* of the interchange with row i in both 1x1 and 2x2 pivot cases) */
        i__1 = *n;
        for (k = 1;
                k <= i__1;
                ++k)
        {
            kp = (i__2 = ipiv[k], f2c_dabs(i__2));
            if (kp != k)
            {
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
        }
    }
    else
    {
        /* Begin Lower */
        /* Solve A*X = B, where A = L*D*L**T. */
        /* P**T * B */
        /* Interchange rows K and IPIV(K) of matrix B in the same order */
        /* that the formation order of IPIV(I) vector for Lower case. */
        /* (We can do the simple loop over IPIV with increment 1, */
        /* since the ABS value of IPIV(I) represents the row index */
        /* of the interchange with row i in both 1x1 and 2x2 pivot cases) */
        i__1 = *n;
        for (k = 1;
                k <= i__1;
                ++k)
        {
            kp = (i__2 = ipiv[k], f2c_dabs(i__2));
            if (kp != k)
            {
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
        }
        /* Compute (L \P**T * B) -> B [ (L \P**T * B) ] */
        ztrsm_("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[ b_offset], ldb);
        /* Compute D \ B -> B [ D \ (L \P**T * B) ] */
        i__ = 1;
        while(i__ <= *n)
        {
            if (ipiv[i__] > 0)
            {
                z_div(&z__1, &c_b1, &a[i__ + i__ * a_dim1]);
                zscal_(nrhs, &z__1, &b[i__ + b_dim1], ldb);
            }
            else if (i__ < *n)
            {
                i__1 = i__;
                akm1k.r = e[i__1].r;
                akm1k.i = e[i__1].i; // , expr subst
                z_div(&z__1, &a[i__ + i__ * a_dim1], &akm1k);
                akm1.r = z__1.r;
                akm1.i = z__1.i; // , expr subst
                z_div(&z__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
                ak.r = z__1.r;
                ak.i = z__1.i; // , expr subst
                z__2.r = akm1.r * ak.r - akm1.i * ak.i;
                z__2.i = akm1.r * ak.i + akm1.i * ak.r; // , expr subst
                z__1.r = z__2.r - 1.;
                z__1.i = z__2.i - 0.; // , expr subst
                denom.r = z__1.r;
                denom.i = z__1.i; // , expr subst
                i__1 = *nrhs;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    z_div(&z__1, &b[i__ + j * b_dim1], &akm1k);
                    bkm1.r = z__1.r;
                    bkm1.i = z__1.i; // , expr subst
                    z_div(&z__1, &b[i__ + 1 + j * b_dim1], &akm1k);
                    bk.r = z__1.r;
                    bk.i = z__1.i; // , expr subst
                    i__2 = i__ + j * b_dim1;
                    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                    z__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                    z__2.r = z__3.r - bk.r;
                    z__2.i = z__3.i - bk.i; // , expr subst
                    z_div(&z__1, &z__2, &denom);
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                    i__2 = i__ + 1 + j * b_dim1;
                    z__3.r = akm1.r * bk.r - akm1.i * bk.i;
                    z__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                    z__2.r = z__3.r - bkm1.r;
                    z__2.i = z__3.i - bkm1.i; // , expr subst
                    z_div(&z__1, &z__2, &denom);
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                }
                ++i__;
            }
            ++i__;
        }
        /* Compute (L**T \ B) -> B [ L**T \ (D \ (L \P**T * B) ) ] */
        ztrsm_("L", "L", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[ b_offset], ldb);
        /* P * B [ P * (L**T \ (D \ (L \P**T * B) )) ] */
        /* Interchange rows K and IPIV(K) of matrix B in reverse order */
        /* from the formation order of IPIV(I) vector for Lower case. */
        /* (We can do the simple loop over IPIV with decrement -1, */
        /* since the ABS value of IPIV(I) represents the row index */
        /* of the interchange with row i in both 1x1 and 2x2 pivot cases) */
        for (k = *n;
                k >= 1;
                --k)
        {
            kp = (i__1 = ipiv[k], f2c_dabs(i__1));
            if (kp != k)
            {
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
        }
        /* END Lower */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZSYTRS_3 */
}
/* zsytrs_3__ */

