/* ../netlib/v3.9.0/zhbev_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__4 = 4;
static doublereal c_b21 = 1.;
static integer c__1 = 1;
/* > \brief <b> ZHBEV_2STAGE computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b> */
/* @precisions fortran z -> s d c */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHBEV_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbev_2 stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbev_2 stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbev_2 stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHBEV_2STAGE( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, */
/* WORK, LWORK, RWORK, INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, KD, LDAB, LDZ, N, LWORK */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION RWORK( * ), W( * ) */
/* COMPLEX*16 AB( LDAB, * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHBEV_2STAGE computes all the eigenvalues and, optionally, eigenvectors of */
/* > a complex Hermitian band matrix A using the 2stage technique for */
/* > the reduction to tridiagonal. */
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
/* > AB is COMPLEX*16 array, dimension (LDAB, N) */
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
/* > W is DOUBLE PRECISION array, dimension (N) */
/* > If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is COMPLEX*16 array, dimension (LDZ, N) */
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
/* > WORK is COMPLEX*16 array, dimension LWORK */
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
/* > RWORK is DOUBLE PRECISION array, dimension (max(1,3*N-2)) */
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
/* > \ingroup complex16OTHEReigen */
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
int zhbev_2stage_(char *jobz, char *uplo, integer *n, integer *kd, doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhbev inputs: jobz %c, uplo %c, n %" FLA_IS ", kd %" FLA_IS ", ldab %" FLA_IS ", ldz %" FLA_IS "",*jobz, *uplo, *n, *kd, *ldab, *ldz);
    /* System generated locals */
    integer ab_dim1, ab_offset, z_dim1, z_offset, i__1;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer ib;
    doublereal eps;
    integer inde;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    doublereal anrm;
    integer imax;
    extern /* Subroutine */
    int zhetrd_hb2st_(char *, char *, char *, integer *, integer *, doublecomplex *, integer *, doublereal *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    doublereal rmin, rmax;
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *);
    integer iinfo, lhtrd, lwmin;
    logical lower;
    integer lwtrd;
    logical wantz;
    extern doublereal dlamch_(char *);
    integer iscale;
    doublereal safmin;
    extern doublereal zlanhb_(char *, char *, integer *, integer *, doublecomplex *, integer *, doublereal *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal bignum;
    extern /* Subroutine */
    int dsterf_(integer *, doublereal *, doublereal *, integer *), zlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublecomplex *, integer *, integer *);
    integer indwrk, indrwk, llwork;
    doublereal smlnum;
    logical lquery;
    extern /* Subroutine */
    int zsteqr_(char *, integer *, doublereal *, doublereal *, doublecomplex *, integer *, doublereal *, integer *);
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
    /* Function Body */
    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");
    lquery = *lwork == -1;
    *info = 0;
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
        if (*n <= 1)
        {
            lwmin = 1;
            work[1].r = (doublereal) lwmin;
            work[1].i = 0.; // , expr subst
        }
        else
        {
            ib = ilaenv2stage_(&c__2, "ZHETRD_HB2ST", jobz, n, kd, &c_n1, & c_n1);
            lhtrd = ilaenv2stage_(&c__3, "ZHETRD_HB2ST", jobz, n, kd, &ib, & c_n1);
            lwtrd = ilaenv2stage_(&c__4, "ZHETRD_HB2ST", jobz, n, kd, &ib, & c_n1);
            lwmin = lhtrd + lwtrd;
            work[1].r = (doublereal) lwmin;
            work[1].i = 0.; // , expr subst
        }
        if (*lwork < lwmin && ! lquery)
        {
            *info = -11;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHBEV_2STAGE ", &i__1);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (lquery)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (*n == 1)
    {
        if (lower)
        {
            i__1 = ab_dim1 + 1;
            w[1] = ab[i__1].r;
        }
        else
        {
            i__1 = *kd + 1 + ab_dim1;
            w[1] = ab[i__1].r;
        }
        if (wantz)
        {
            i__1 = z_dim1 + 1;
            z__[i__1].r = 1.;
            z__[i__1].i = 0.; // , expr subst
        }
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Get machine constants. */
    safmin = dlamch_("Safe minimum");
    eps = dlamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
    /* Scale matrix to allowable range, if necessary. */
    anrm = zlanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1]);
    iscale = 0;
    if (anrm > 0. && anrm < rmin)
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
            zlascl_("B", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
        else
        {
            zlascl_("Q", kd, kd, &c_b21, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
    }
    /* Call ZHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form. */
    inde = 1;
    indhous = 1;
    indwrk = indhous + lhtrd;
    llwork = *lwork - indwrk + 1;
    zhetrd_hb2st_("N", jobz, uplo, n, kd, &ab[ab_offset], ldab, &w[1], & rwork[inde], &work[indhous], &lhtrd, &work[indwrk], &llwork, & iinfo);
    /* For eigenvalues only, call DSTERF. For eigenvectors, call ZSTEQR. */
    if (! wantz)
    {
        dsterf_(n, &w[1], &rwork[inde], info);
    }
    else
    {
        indrwk = inde + *n;
        zsteqr_(jobz, n, &w[1], &rwork[inde], &z__[z_offset], ldz, &rwork[ indrwk], info);
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
        d__1 = 1. / sigma;
        dscal_(&imax, &d__1, &w[1], &c__1);
    }
    /* Set WORK(1) to optimal workspace size. */
    work[1].r = (doublereal) lwmin;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZHBEV_2STAGE */
}
/* zhbev_2stage__ */
