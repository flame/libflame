/* ../netlib/v3.9.0/chegv_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    1.f,0.f
}
;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
/* > \brief \b CHEGV_2STAGE */
/* @generated from zhegv_2stage.f, fortran z -> c, Sun Nov 6 13:09:52 2016 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHEGV_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegv_2 stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegv_2 stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegv_2 stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHEGV_2STAGE( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, */
/* WORK, LWORK, RWORK, INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, ITYPE, LDA, LDB, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* REAL RWORK( * ), W( * ) */
/* COMPLEX A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEGV_2STAGE computes all the eigenvalues, and optionally, the eigenvectors */
/* > of a complex generalized Hermitian-definite eigenproblem, of the form */
/* > A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x. */
/* > Here A and B are assumed to be Hermitian and B is also */
/* > positive definite. */
/* > This routine use the 2stage technique for the reduction to tridiagonal */
/* > which showed higher performance on recent architecture and for large */
/* > sizes N>2000. */
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
/* > Not available in this release. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA, N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the */
/* > leading N-by-N upper triangular part of A contains the */
/* > upper triangular part of the matrix A. If UPLO = 'L', */
/* > the leading N-by-N lower triangular part of A contains */
/* > the lower triangular part of the matrix A. */
/* > */
/* > On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
/* > matrix Z of eigenvectors. The eigenvectors are normalized */
/* > as follows: */
/* > if ITYPE = 1 or 2, Z**H*B*Z = I;
*/
/* > if ITYPE = 3, Z**H*inv(B)*Z = I. */
/* > If JOBZ = 'N', then on exit the upper triangle (if UPLO='U') */
/* > or the lower triangle (if UPLO='L') of A, including the */
/* > diagonal, is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB, N) */
/* > On entry, the Hermitian positive definite matrix B. */
/* > If UPLO = 'U', the leading N-by-N upper triangular part of B */
/* > contains the upper triangular part of the matrix B. */
/* > If UPLO = 'L', the leading N-by-N lower triangular part of B */
/* > contains the lower triangular part of the matrix B. */
/* > */
/* > On exit, if INFO <= N, the part of B containing the matrix is */
/* > overwritten by the triangular factor U or L from the Cholesky */
/* > factorization B = U**H*U or B = L*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
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
/* > The length of the array WORK. LWORK >= 1, when N <= 1;
*/
/* > otherwise */
/* > If JOBZ = 'N' and N > 1, LWORK must be queried. */
/* > LWORK = MAX(1, dimension) where */
/* > dimension = fla_max(stage1,stage2) + (KD+1)*N + N */
/* > = N*KD + N*max(KD+1,FACTOPTNB) */
/* > + fla_max(2*KD*KD, KD*NTHREADS) */
/* > + (KD+1)*N + N */
/* > where KD is the blocking size of the reduction, */
/* > FACTOPTNB is the blocking used by the QR or LQ */
/* > algorithm, usually FACTOPTNB=128 is a good choice */
/* > NTHREADS is the number of threads used when */
/* > openMP compilation is enabled, otherwise =1. */
/* > If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (fla_max(1, 3*N-2)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: CPOTRF or CHEEV returned an error code: */
/* > <= N: if INFO = i, CHEEV failed to converge;
*/
/* > i off-diagonal elements of an intermediate */
/* > tridiagonal form did not converge to zero;
*/
/* > > N: if INFO = N + i, for 1 <= i <= N, then the leading */
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
/* > \date November 2017 */
/* > \ingroup complexHEeigen */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > All details about the 2stage techniques are available in: */
/* > */
/* > Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* > Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* > using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* > of 2011 International Conference for High Performance Computing, */
/* > Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* > Article 8 , 11 pages. */
/* > http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* > A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* > An improved parallel singular value algorithm and its implementation */
/* > for multicore hardware, In Proceedings of 2013 International Conference */
/* > for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* > Denver, Colorado, USA, 2013. */
/* > Article 90, 12 pages. */
/* > http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* > A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* > A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* > calculations based on fine-grained memory aware tasks. */
/* > International Journal of High Performance Computing Applications. */
/* > Volume 28 Issue 2, Pages 196-209, May 2014. */
/* > http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int chegv_2stage_(integer *itype, char *jobz, char *uplo, integer *n, complex *a, integer *lda, complex *b, integer *ldb, real * w, complex *work, integer *lwork, real *rwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"chegv_2stage inputs: itype %lld, jobz %c, uplo %c, n %lld, lda %lld, ldb %lld, lwork %lld",*itype, *jobz, *uplo, *n, *lda, *ldb, *lwork);
#else
    snprintf(buffer, 256,"chegv_2stage inputs: itype %d, jobz %c, uplo %c, n %d, lda %d, ldb %d, lwork %d",*itype, *jobz, *uplo, *n, *lda, *ldb, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    integer ib, kd;
    extern /* Subroutine */
    int cheev_2stage_(char *, char *, integer *, complex *, integer *, real *, complex *, integer *, real *, integer *);
    integer neig;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern logical lsame_(char *, char *);
    integer lhtrd;
    extern /* Subroutine */
    int ctrmm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *);
    integer lwmin;
    char trans[1];
    extern /* Subroutine */
    int ctrsm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *);
    logical upper;
    integer lwtrd;
    logical wantz;
    extern /* Subroutine */
    int xerbla_(char *, integer *), chegst_( integer *, char *, integer *, complex *, integer *, complex *, integer *, integer *), cpotrf_(char *, integer *, complex *, integer *, integer *);
    logical lquery;
    /* -- LAPACK driver routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
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
    --w;
    --work;
    --rwork;
    /* Function Body */
    wantz = lsame_(jobz, "V");
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1;
    *info = 0;
    if (*itype < 1 || *itype > 3)
    {
        *info = -1;
    }
    else if (! lsame_(jobz, "N"))
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
    else if (*lda < fla_max(1,*n))
    {
        *info = -6;
    }
    else if (*ldb < fla_max(1,*n))
    {
        *info = -8;
    }
    if (*info == 0)
    {
        kd = ilaenv2stage_(&c__1, "CHETRD_2STAGE", jobz, n, &c_n1, &c_n1, & c_n1);
        ib = ilaenv2stage_(&c__2, "CHETRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1);
        lhtrd = ilaenv2stage_(&c__3, "CHETRD_2STAGE", jobz, n, &kd, &ib, & c_n1);
        lwtrd = ilaenv2stage_(&c__4, "CHETRD_2STAGE", jobz, n, &kd, &ib, & c_n1);
        lwmin = *n + lhtrd + lwtrd;
        work[1].r = (real) lwmin;
        work[1].i = 0.f; // , expr subst
        if (*lwork < lwmin && ! lquery)
        {
            *info = -11;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHEGV_2STAGE ", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Form a Cholesky factorization of B. */
    cpotrf_(uplo, n, &b[b_offset], ldb, info);
    if (*info != 0)
    {
        *info = *n + *info;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Transform problem to standard eigenvalue problem and solve. */
    chegst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info);
    cheev_2stage_(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, & rwork[1], info);
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
            /* backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y */
            if (upper)
            {
                *(unsigned char *)trans = 'N';
            }
            else
            {
                *(unsigned char *)trans = 'C';
            }
            ctrsm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b1, &b[ b_offset], ldb, &a[a_offset], lda);
        }
        else if (*itype == 3)
        {
            /* For B*A*x=(lambda)*x;
            */
            /* backtransform eigenvectors: x = L*y or U**H *y */
            if (upper)
            {
                *(unsigned char *)trans = 'C';
            }
            else
            {
                *(unsigned char *)trans = 'N';
            }
            ctrmm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b1, &b[ b_offset], ldb, &a[a_offset], lda);
        }
    }
    work[1].r = (real) lwmin;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CHEGV_2STAGE */
}
/* chegv_2stage__ */
