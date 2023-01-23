/* chetrf_aa.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b2 =
{
    1.f,0.f
}
;
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b CHETRF_AA */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHETRF_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrf_ aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrf_ aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrf_ aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER N, LDA, LWORK, INFO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHETRF_AA computes the factorization of a complex hermitian matrix A */
/* > using the Aasen's algorithm. The form of the factorization is */
/* > */
/* > A = U**H*T*U or A = L*T*L**H */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and T is a hermitian tridiagonal matrix. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
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
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the hermitian matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, the tridiagonal matrix is stored in the diagonals */
/* > and the subdiagonals of A just below (or above) the diagonals, */
/* > and L is stored below (or above) the subdiaonals, when UPLO */
/* > is 'L' (or 'U'). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > On exit, it contains the details of the interchanges, i.e., */
/* > the row and column k of A were interchanged with the */
/* > row and column IPIV(k). */
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
/* > The length of WORK. LWORK >= 2*N. For optimum performance */
/* > LWORK >= N*(1+NB), where NB is the optimal blocksize. */
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
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup complexHEcomputational */
/* ===================================================================== */
/* Subroutine */
int chetrf_aa_(char *uplo, integer *n, complex *a, integer * lda, integer *ipiv, complex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("chetrf inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "",*uplo, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1;
    complex q__1;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    extern /* Subroutine */
    int clahef_aa_(char *, integer *, integer *, integer *, complex *, integer *, integer *, complex *, integer *, complex *);
    integer j, j1, k1, k2, j2, j3, jb, nb, mj, nj;
    complex alpha;
    extern /* Subroutine */
    int cscal_(integer *, complex *, complex *, integer *), cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int ccopy_(integer *, complex *, integer *, complex *, integer *), cswap_(integer *, complex *, integer *, complex *, integer *);
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Determine the block size */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    /* Function Body */
    nb = ilaenv_(&c__1, "CHETRF_AA", uplo, n, &c_n1, &c_n1, &c_n1);
    /* Test the input parameters. */
    *info = 0;
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1;
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*n))
    {
        *info = -4;
    }
    else if (*lwork < *n << 1 && ! lquery)
    {
        *info = -7;
    }
    if (*info == 0)
    {
        lwkopt = (nb + 1) * *n;
        work[1].r = (real) lwkopt;
        work[1].i = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHETRF_AA", &i__1);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (lquery)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return */
    if (*n == 0)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    ipiv[1] = 1;
    if (*n == 1)
    {
        i__1 = a_dim1 + 1;
        i__2 = a_dim1 + 1;
        r__1 = a[i__2].r;
        a[i__1].r = r__1;
        a[i__1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Adjust block size based on the workspace size */
    if (*lwork < (nb + 1) * *n)
    {
        nb = (*lwork - *n) / *n;
    }
    if (upper)
    {
        /* ..................................................... */
        /* Factorize A as U**H*D*U using the upper triangle of A */
        /* ..................................................... */
        /* copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N)) */
        ccopy_(n, &a[a_dim1 + 1], lda, &work[1], &c__1);
        /* J is the main loop index, increasing from 1 to N in steps of */
        /* JB, where JB is the number of columns factorized by CLAHEF;
        */
        /* JB is either NB, or N-J+1 for the last block */
        j = 0;
L10:
        if (j >= *n)
        {
            goto L20;
        }
        /* each step of the main loop */
        /* J is the last column of the previous panel */
        /* J1 is the first column of the current panel */
        /* K1 identifies if the previous column of the panel has been */
        /* explicitly stored, e.g., K1=1 for the first panel, and */
        /* K1=0 for the rest */
        j1 = j + 1;
        /* Computing MIN */
        i__1 = *n - j1 + 1;
        jb = fla_min(i__1,nb);
        k1 = fla_max(1,j) - j;
        /* Panel factorization */
        i__1 = 2 - k1;
        i__2 = *n - j;
        clahef_aa_(uplo, &i__1, &i__2, &jb, &a[fla_max(1,j) + (j + 1) * a_dim1], lda, &ipiv[j + 1], &work[1], n, &work[*n * nb + 1]) ;
        /* Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot) */
        /* Computing MIN */
        i__2 = *n;
        i__3 = j + jb + 1; // , expr subst
        i__1 = fla_min(i__2,i__3);
        for (j2 = j + 2;
                j2 <= i__1;
                ++j2)
        {
            ipiv[j2] += j;
            if (j2 != ipiv[j2] && j1 - k1 > 2)
            {
                i__2 = j1 - k1 - 2;
                cswap_(&i__2, &a[j2 * a_dim1 + 1], &c__1, &a[ipiv[j2] * a_dim1 + 1], &c__1);
            }
        }
        j += jb;
        /* Trailing submatrix update, where */
        /* the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and */
        /* WORK stores the current block of the auxiriarly matrix H */
        if (j < *n)
        {
            /* if the first panel and JB=1 (NB=1), then nothing to do */
            if (j1 > 1 || jb > 1)
            {
                /* Merge rank-1 update with BLAS-3 update */
                r_cnjg(&q__1, &a[j + (j + 1) * a_dim1]);
                alpha.r = q__1.r;
                alpha.i = q__1.i; // , expr subst
                i__1 = j + (j + 1) * a_dim1;
                a[i__1].r = 1.f;
                a[i__1].i = 0.f; // , expr subst
                i__1 = *n - j;
                ccopy_(&i__1, &a[j - 1 + (j + 1) * a_dim1], lda, &work[j + 1 - j1 + 1 + jb * *n], &c__1);
                i__1 = *n - j;
                cscal_(&i__1, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);
                /* K1 identifies if the previous column of the panel has been */
                /* explicitly stored, e.g., K1=0 and K2=1 for the first panel, */
                /* and K1=1 and K2=0 for the rest */
                if (j1 > 1)
                {
                    /* Not first panel */
                    k2 = 1;
                }
                else
                {
                    /* First panel */
                    k2 = 0;
                    /* First update skips the first column */
                    --jb;
                }
                i__1 = *n;
                i__2 = nb;
                for (j2 = j + 1;
                        i__2 < 0 ? j2 >= i__1 : j2 <= i__1;
                        j2 += i__2)
                {
                    /* Computing MIN */
                    i__3 = nb;
                    i__4 = *n - j2 + 1; // , expr subst
                    nj = fla_min(i__3,i__4);
                    /* Update (J2, J2) diagonal block with CGEMV */
                    j3 = j2;
                    for (mj = nj - 1;
                            mj >= 1;
                            --mj)
                    {
                        i__3 = jb + 1;
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("Conjugate transpose", "Transpose", &c__1, &mj, &i__3, &q__1, &a[j1 - k2 + j3 * a_dim1], lda, &work[j3 - j1 + 1 + k1 * *n], n, &c_b2, &a[ j3 + j3 * a_dim1], lda) ;
                        ++j3;
                    }
                    /* Update off-diagonal block of J2-th block row with CGEMM */
                    i__3 = *n - j3 + 1;
                    i__4 = jb + 1;
                    q__1.r = -1.f;
                    q__1.i = -0.f; // , expr subst
                    cgemm_("Conjugate transpose", "Transpose", &nj, &i__3, & i__4, &q__1, &a[j1 - k2 + j2 * a_dim1], lda, & work[j3 - j1 + 1 + k1 * *n], n, &c_b2, &a[j2 + j3 * a_dim1], lda);
                }
                /* Recover T( J, J+1 ) */
                i__2 = j + (j + 1) * a_dim1;
                r_cnjg(&q__1, &alpha);
                a[i__2].r = q__1.r;
                a[i__2].i = q__1.i; // , expr subst
            }
            /* WORK(J+1, 1) stores H(J+1, 1) */
            i__2 = *n - j;
            ccopy_(&i__2, &a[j + 1 + (j + 1) * a_dim1], lda, &work[1], &c__1);
        }
        goto L10;
    }
    else
    {
        /* ..................................................... */
        /* Factorize A as L*D*L**H using the lower triangle of A */
        /* ..................................................... */
        /* copy first column A(1:N, 1) into H(1:N, 1) */
        /* (stored in WORK(1:N)) */
        ccopy_(n, &a[a_dim1 + 1], &c__1, &work[1], &c__1);
        /* J is the main loop index, increasing from 1 to N in steps of */
        /* JB, where JB is the number of columns factorized by CLAHEF;
        */
        /* JB is either NB, or N-J+1 for the last block */
        j = 0;
L11:
        if (j >= *n)
        {
            goto L20;
        }
        /* each step of the main loop */
        /* J is the last column of the previous panel */
        /* J1 is the first column of the current panel */
        /* K1 identifies if the previous column of the panel has been */
        /* explicitly stored, e.g., K1=1 for the first panel, and */
        /* K1=0 for the rest */
        j1 = j + 1;
        /* Computing MIN */
        i__2 = *n - j1 + 1;
        jb = fla_min(i__2,nb);
        k1 = fla_max(1,j) - j;
        /* Panel factorization */
        i__2 = 2 - k1;
        i__1 = *n - j;
        clahef_aa_(uplo, &i__2, &i__1, &jb, &a[j + 1 + fla_max(1,j) * a_dim1], lda, &ipiv[j + 1], &work[1], n, &work[*n * nb + 1]) ;
        /* Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot) */
        /* Computing MIN */
        i__1 = *n;
        i__3 = j + jb + 1; // , expr subst
        i__2 = fla_min(i__1,i__3);
        for (j2 = j + 2;
                j2 <= i__2;
                ++j2)
        {
            ipiv[j2] += j;
            if (j2 != ipiv[j2] && j1 - k1 > 2)
            {
                i__1 = j1 - k1 - 2;
                cswap_(&i__1, &a[j2 + a_dim1], lda, &a[ipiv[j2] + a_dim1], lda);
            }
        }
        j += jb;
        /* Trailing submatrix update, where */
        /* A(J2+1, J1-1) stores L(J2+1, J1) and */
        /* WORK(J2+1, 1) stores H(J2+1, 1) */
        if (j < *n)
        {
            /* if the first panel and JB=1 (NB=1), then nothing to do */
            if (j1 > 1 || jb > 1)
            {
                /* Merge rank-1 update with BLAS-3 update */
                r_cnjg(&q__1, &a[j + 1 + j * a_dim1]);
                alpha.r = q__1.r;
                alpha.i = q__1.i; // , expr subst
                i__2 = j + 1 + j * a_dim1;
                a[i__2].r = 1.f;
                a[i__2].i = 0.f; // , expr subst
                i__2 = *n - j;
                ccopy_(&i__2, &a[j + 1 + (j - 1) * a_dim1], &c__1, &work[j + 1 - j1 + 1 + jb * *n], &c__1);
                i__2 = *n - j;
                cscal_(&i__2, &alpha, &work[j + 1 - j1 + 1 + jb * *n], &c__1);
                /* K1 identifies if the previous column of the panel has been */
                /* explicitly stored, e.g., K1=0 and K2=1 for the first panel, */
                /* and K1=1 and K2=0 for the rest */
                if (j1 > 1)
                {
                    /* Not first panel */
                    k2 = 1;
                }
                else
                {
                    /* First panel */
                    k2 = 0;
                    /* First update skips the first column */
                    --jb;
                }
                i__2 = *n;
                i__1 = nb;
                for (j2 = j + 1;
                        i__1 < 0 ? j2 >= i__2 : j2 <= i__2;
                        j2 += i__1)
                {
                    /* Computing MIN */
                    i__3 = nb;
                    i__4 = *n - j2 + 1; // , expr subst
                    nj = fla_min(i__3,i__4);
                    /* Update (J2, J2) diagonal block with CGEMV */
                    j3 = j2;
                    for (mj = nj - 1;
                            mj >= 1;
                            --mj)
                    {
                        i__3 = jb + 1;
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("No transpose", "Conjugate transpose", &mj, & c__1, &i__3, &q__1, &work[j3 - j1 + 1 + k1 * * n], n, &a[j3 + (j1 - k2) * a_dim1], lda, & c_b2, &a[j3 + j3 * a_dim1], lda);
                        ++j3;
                    }
                    /* Update off-diagonal block of J2-th block column with CGEMM */
                    i__3 = *n - j3 + 1;
                    i__4 = jb + 1;
                    q__1.r = -1.f;
                    q__1.i = -0.f; // , expr subst
                    cgemm_("No transpose", "Conjugate transpose", &i__3, &nj, &i__4, &q__1, &work[j3 - j1 + 1 + k1 * *n], n, &a[ j2 + (j1 - k2) * a_dim1], lda, &c_b2, &a[j3 + j2 * a_dim1], lda);
                }
                /* Recover T( J+1, J ) */
                i__1 = j + 1 + j * a_dim1;
                r_cnjg(&q__1, &alpha);
                a[i__1].r = q__1.r;
                a[i__1].i = q__1.i; // , expr subst
            }
            /* WORK(J+1, 1) stores H(J+1, 1) */
            i__1 = *n - j;
            ccopy_(&i__1, &a[j + 1 + (j + 1) * a_dim1], &c__1, &work[1], & c__1);
        }
        goto L11;
    }
L20:
    work[1].r = (real) lwkopt;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of CHETRF_AA */
}
/* chetrf_aa__ */
