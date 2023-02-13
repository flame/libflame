/* ../netlib/chpevd.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief <b> CHPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHPEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chpevd. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chpevd. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chpevd. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, */
/* RWORK, LRWORK, IWORK, LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, UPLO */
/* INTEGER INFO, LDZ, LIWORK, LRWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL RWORK( * ), W( * ) */
/* COMPLEX AP( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPEVD computes all the eigenvalues and, optionally, eigenvectors of */
/* > a complex Hermitian matrix A in packed storage. If eigenvectors are */
/* > desired, it uses a divide and conquer algorithm. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is COMPLEX array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the Hermitian matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* > On exit, AP is overwritten by values generated during the */
/* > reduction to tridiagonal form. If UPLO = 'U', the diagonal */
/* > and first superdiagonal of the tridiagonal matrix T overwrite */
/* > the corresponding elements of A, and if UPLO = 'L', the */
/* > diagonal and first subdiagonal of T overwrite the */
/* > corresponding elements of A. */
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
/* > On exit, if INFO = 0, WORK(1) returns the required LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of array WORK. */
/* > If N <= 1, LWORK must be at least 1. */
/* > If JOBZ = 'N' and N > 1, LWORK must be at least N. */
/* > If JOBZ = 'V' and N > 1, LWORK must be at least 2*N. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the required sizes of the WORK, RWORK and */
/* > IWORK arrays, returns these values as the first entries of */
/* > the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* > On exit, if INFO = 0, RWORK(1) returns the required LRWORK. */
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
/* > routine only calculates the required sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* > On exit, if INFO = 0, IWORK(1) returns the required LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK is INTEGER */
/* > The dimension of array IWORK. */
/* > If JOBZ = 'N' or N <= 1, LIWORK must be at least 1. */
/* > If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N. */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the required sizes of the WORK, RWORK */
/* > and IWORK arrays, returns these values as the first entries */
/* > of the WORK, RWORK and IWORK arrays, and no error message */
/* > related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
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
/* > \date November 2011 */
/* > \ingroup complexOTHEReigen */
/* ===================================================================== */
/* Subroutine */
int chpevd_(char *jobz, char *uplo, integer *n, complex *ap, real *w, complex *z__, integer *ldz, complex *work, integer *lwork, real *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    real eps;
    integer inde;
    real anrm;
    integer imax;
    real rmin, rmax, sigma;
    extern logical lsame_(char *, char *);
    integer iinfo;
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    integer lwmin, llrwk, llwrk;
    logical wantz;
    integer iscale;
    extern real clanhp_(char *, char *, integer *, complex *, real *);
    extern /* Subroutine */
    int cstedc_(char *, integer *, real *, real *, complex *, integer *, complex *, integer *, real *, integer *, integer *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int csscal_(integer *, real *, complex *, integer *);
    real safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real bignum;
    integer indtau;
    extern /* Subroutine */
    int chptrd_(char *, integer *, complex *, real *, real *, complex *, integer *);
    integer indrwk, indwrk, liwmin;
    extern /* Subroutine */
    int ssterf_(integer *, real *, real *, integer *);
    integer lrwmin;
    extern /* Subroutine */
    int cupmtr_(char *, char *, char *, integer *, integer *, complex *, complex *, complex *, integer *, complex *, integer *);
    real smlnum;
    logical lquery;
    /* -- LAPACK driver routine (version 3.4.0) -- */
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
    --ap;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    --iwork;
    /* Function Body */
    wantz = lsame_(jobz, "V");
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;
    *info = 0;
    if (! (wantz || lsame_(jobz, "N")))
    {
        *info = -1;
    }
    else if (! (lsame_(uplo, "L") || lsame_(uplo, "U")))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*ldz < 1 || wantz && *ldz < *n)
    {
        *info = -7;
    }
    if (*info == 0)
    {
        if (*n <= 1)
        {
            lwmin = 1;
            liwmin = 1;
            lrwmin = 1;
        }
        else
        {
            if (wantz)
            {
                lwmin = *n << 1;
                /* Computing 2nd power */
                i__1 = *n;
                lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
                liwmin = *n * 5 + 3;
            }
            else
            {
                lwmin = *n;
                lrwmin = *n;
                liwmin = 1;
            }
        }
        work[1].r = (real) lwmin;
        work[1].i = 0.f; // , expr subst
        rwork[1] = (real) lrwmin;
        iwork[1] = liwmin;
        if (*lwork < lwmin && ! lquery)
        {
            *info = -9;
        }
        else if (*lrwork < lrwmin && ! lquery)
        {
            *info = -11;
        }
        else if (*liwork < liwmin && ! lquery)
        {
            *info = -13;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHPEVD", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    if (*n == 1)
    {
        w[1] = ap[1].r;
        if (wantz)
        {
            i__1 = z_dim1 + 1;
            z__[i__1].r = 1.f;
            z__[i__1].i = 0.f; // , expr subst
        }
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
    anrm = clanhp_("M", uplo, n, &ap[1], &rwork[1]);
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
        i__1 = *n * (*n + 1) / 2;
        csscal_(&i__1, &sigma, &ap[1], &c__1);
    }
    /* Call CHPTRD to reduce Hermitian packed matrix to tridiagonal form. */
    inde = 1;
    indtau = 1;
    indrwk = inde + *n;
    indwrk = indtau + *n;
    llwrk = *lwork - indwrk + 1;
    llrwk = *lrwork - indrwk + 1;
    chptrd_(uplo, n, &ap[1], &w[1], &rwork[inde], &work[indtau], &iinfo);
    /* For eigenvalues only, call SSTERF. For eigenvectors, first call */
    /* CUPGTR to generate the orthogonal matrix, then call CSTEDC. */
    if (! wantz)
    {
        ssterf_(n, &w[1], &rwork[inde], info);
    }
    else
    {
        cstedc_("I", n, &w[1], &rwork[inde], &z__[z_offset], ldz, &work[ indwrk], &llwrk, &rwork[indrwk], &llrwk, &iwork[1], liwork, info);
        cupmtr_("L", uplo, "N", n, n, &ap[1], &work[indtau], &z__[z_offset], ldz, &work[indwrk], &iinfo);
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
    return 0;
    /* End of CHPEVD */
}
/* chpevd_ */
