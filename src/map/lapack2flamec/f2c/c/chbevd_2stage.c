/* ../netlib/v3.9.0/chbevd_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    0.f,0.f
}
;
static complex c_b2 =
{
    1.f,0.f
}
;
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__4 = 4;
static real c_b23 = 1.f;
static integer c__1 = 1;
/* > \brief <b> CHBEVD_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b> */
/* @generated from zhbevd_2stage.f, fortran z -> c, Sat Nov 5 23:18:17 2016 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHBEVD_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbevd_ 2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbevd_ 2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbevd_ 2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHBEVD_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, */
/* WORK, LWORK, RWORK, LRWORK, IWORK, */
/* LIWORK, INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL RWORK( * ), W( * ) */
/* COMPLEX AB( LDAB, * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBEVD_2STAGE computes all the eigenvalues and, optionally, eigenvectors of */
/* > a complex Hermitian band matrix A using the 2stage technique for */
/* > the reduction to tridiagonal. If eigenvectors are desired, it */
/* > uses a divide and conquer algorithm. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB, N) */
/* > On entry, the upper or lower triangle of the Hermitian band */
/* > matrix A, stored in the first KD+1 rows of the array. The */
/* > j-th column of A is stored in the j-th column of the array AB */
/* > as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*/
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* > On exit, AB is overwritten by values generated during the */
/* > reduction to tridiagonal form. If UPLO = 'U', the first */
/* > superdiagonal and the diagonal of the tridiagonal matrix T */
/* > are returned in rows KD and KD+1 of AB, and if UPLO = 'L', */
/* > the diagonal and first subdiagonal of T are returned in the */
/* > first two rows of AB. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD + 1. */
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
/* > Z is COMPLEX array, dimension (LDZ, N) */
/* > If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal */
/* > eigenvectors of the matrix A, with the i-th column of Z */
/* > holding the eigenvector associated with W(i). */
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
/* > dimension = (2KD+1)*N + KD*NTHREADS */
/* > where KD is the size of the band. */
/* > NTHREADS is the number of threads used when */
/* > openMP compilation is enabled, otherwise =1. */
/* > If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal sizes of the WORK, RWORK and */
/* > IWORK arrays, returns these values as the first entries of */
/* > the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, */
/* > dimension (LRWORK) */
/* > On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* > LRWORK is INTEGER */
/* > The dimension of array RWORK. */
/* > If N <= 1, LRWORK must be at least 1. */
/* > If JOBZ = 'N' and N > 1, LRWORK must be at least N. */
/* > If JOBZ = 'V' and N > 1, LRWORK must be at least */
/* > 1 + 5*N + 2*N**2. */
/* > */
/* > If LRWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* > On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of array IWORK. */
/* > If JOBZ = 'N' or N <= 1, LIWORK must be at least 1. */
/* > If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N . */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = i, the algorithm failed to converge;
i */
/* > off-diagonal elements of an intermediate tridiagonal */
/* > form did not converge to zero. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup complexOTHEReigen */
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
int chbevd_2stage_(char *jobz, char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, real *w, complex *z__, integer *ldz, complex *work, integer *lwork, real *rwork, integer * lrwork, integer *iwork, integer *liwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"chbevd_2stage inputs: jobz %c, uplo %c, n %lld, kd %lld, ldab %lld, ldz %lld, lwork %lld, lrwork %lld, liwork %lld",*jobz, *uplo, *n, *kd, *ldab, *ldz, *lwork, *lrwork, *liwork);
#else
    snprintf(buffer, 256,"chbevd_2stage inputs: jobz %c, uplo %c, n %d, kd %d, ldab %d, ldz %d, lwork %d, lrwork %d, liwork %d",*jobz, *uplo, *n, *kd, *ldab, *ldz, *lwork, *lrwork, *liwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1, i__2;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer ib;
    real eps;
    extern /* Subroutine */
    int chetrd_hb2st_(char *, char *, char *, integer *, integer *, complex *, integer *, real *, real *, complex *, integer *, complex *, integer *, integer *);
    integer inde;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    real anrm;
    integer imax;
    real rmin, rmax;
    integer llwk2;
    extern /* Subroutine */
    int cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    real sigma;
    extern logical lsame_(char *, char *);
    integer iinfo;
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    integer indwk, lhtrd, lwmin;
    logical lower;
    integer lwtrd, llrwk;
    logical wantz;
    integer indwk2;
    extern real clanhb_(char *, char *, integer *, integer *, complex *, integer *, real *);
    integer iscale;
    extern /* Subroutine */
    int clascl_(char *, integer *, integer *, real *, real *, integer *, integer *, complex *, integer *, integer *), cstedc_(char *, integer *, real *, real *, complex *, integer *, complex *, integer *, real *, integer *, integer *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *);
    real safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real bignum;
    integer indrwk, liwmin;
    extern /* Subroutine */
    int ssterf_(integer *, real *, real *, integer *);
    integer lrwmin, llwork;
    real smlnum;
    logical lquery;
    integer indhous;
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
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    --iwork;
    /* Function Body */
    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");
    lquery = *lwork == -1 || *liwork == -1 || *lrwork == -1;
    *info = 0;
    if (*n <= 1)
    {
        lwmin = 1;
        lrwmin = 1;
        liwmin = 1;
    }
    else
    {
        ib = ilaenv2stage_(&c__2, "CHETRD_HB2ST", jobz, n, kd, &c_n1, &c_n1);
        lhtrd = ilaenv2stage_(&c__3, "CHETRD_HB2ST", jobz, n, kd, &ib, &c_n1);
        lwtrd = ilaenv2stage_(&c__4, "CHETRD_HB2ST", jobz, n, kd, &ib, &c_n1);
        if (wantz)
        {
            /* Computing 2nd power */
            i__1 = *n;
            lwmin = i__1 * i__1 << 1;
            /* Computing 2nd power */
            i__1 = *n;
            lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
            liwmin = *n * 5 + 3;
        }
        else
        {
            /* Computing MAX */
            i__1 = *n;
            i__2 = lhtrd + lwtrd; // , expr subst
            lwmin = max(i__1,i__2);
            lrwmin = *n;
            liwmin = 1;
        }
    }
    if (! lsame_(jobz, "N"))
    {
        *info = -1;
    }
    else if (! (lower || lsame_(uplo, "U")))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*kd < 0)
    {
        *info = -4;
    }
    else if (*ldab < *kd + 1)
    {
        *info = -6;
    }
    else if (*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -9;
    }
    if (*info == 0)
    {
        work[1].r = (real) lwmin;
        work[1].i = 0.f; // , expr subst
        rwork[1] = (real) lrwmin;
        iwork[1] = liwmin;
        if (*lwork < lwmin && ! lquery)
        {
            *info = -11;
        }
        else if (*lrwork < lrwmin && ! lquery)
        {
            *info = -13;
        }
        else if (*liwork < liwmin && ! lquery)
        {
            *info = -15;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHBEVD_2STAGE", &i__1);
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
    if (*n == 1)
    {
        i__1 = ab_dim1 + 1;
        w[1] = ab[i__1].r;
        if (wantz)
        {
            i__1 = z_dim1 + 1;
            z__[i__1].r = 1.f;
            z__[i__1].i = 0.f; // , expr subst
        }
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Get machine constants. */
    safmin = slamch_("Safe minimum");
    eps = slamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1.f / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
    /* Scale matrix to allowable range, if necessary. */
    anrm = clanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1]);
    iscale = 0;
    if (anrm > 0.f && anrm < rmin)
    {
        iscale = 1;
        sigma = rmin / anrm;
    }
    else if (anrm > rmax)
    {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1)
    {
        if (lower)
        {
            clascl_("B", kd, kd, &c_b23, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
        else
        {
            clascl_("Q", kd, kd, &c_b23, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
    }
    /* Call CHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form. */
    inde = 1;
    indrwk = inde + *n;
    llrwk = *lrwork - indrwk + 1;
    indhous = 1;
    indwk = indhous + lhtrd;
    llwork = *lwork - indwk + 1;
    indwk2 = indwk + *n * *n;
    llwk2 = *lwork - indwk2 + 1;
    chetrd_hb2st_("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], & rwork[inde], &work[indhous], &lhtrd, &work[indwk], &llwork, & iinfo);
    /* For eigenvalues only, call SSTERF. For eigenvectors, call CSTEDC. */
    if (! wantz)
    {
        ssterf_(n, &w[1], &rwork[inde], info);
    }
    else
    {
        cstedc_("I", n, &w[1], &rwork[inde], &work[1], n, &work[indwk2], & llwk2, &rwork[indrwk], &llrwk, &iwork[1], liwork, info);
        cgemm_("N", "N", n, n, n, &c_b2, &z__[z_offset], ldz, &work[1], n, & c_b1, &work[indwk2], n);
        clacpy_("A", n, n, &work[indwk2], n, &z__[z_offset], ldz);
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if (iscale == 1)
    {
        if (*info == 0)
        {
            imax = *n;
        }
        else
        {
            imax = *info - 1;
        }
        r__1 = 1.f / sigma;
        sscal_(&imax, &r__1, &w[1], &c__1);
    }
    work[1].r = (real) lwmin;
    work[1].i = 0.f; // , expr subst
    rwork[1] = (real) lrwmin;
    iwork[1] = liwmin;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CHBEVD_2STAGE */
}
/* chbevd_2stage__ */
