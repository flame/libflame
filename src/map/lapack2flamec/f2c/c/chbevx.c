/* ../netlib/chbevx.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
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
static real c_b16 = 1.f;
static integer c__1 = 1;
/* > \brief <b> CHBEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHBEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbevx. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbevx. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbevx. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, */
/* VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, */
/* IWORK, IFAIL, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, RANGE, UPLO */
/* INTEGER IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N */
/* REAL ABSTOL, VL, VU */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IFAIL( * ), IWORK( * ) */
/* REAL RWORK( * ), W( * ) */
/* COMPLEX AB( LDAB, * ), Q( LDQ, * ), WORK( * ), */
/* $ Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHBEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a complex Hermitian band matrix A. Eigenvalues and eigenvectors */
/* > can be selected by specifying either a range of values or a range of */
/* > indices for the desired eigenvalues. */
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
/* > \param[in] RANGE */
/* > \verbatim */
/* > RANGE is CHARACTER*1 */
/* > = 'A': all eigenvalues will be found;
*/
/* > = 'V': all eigenvalues in the half-open interval (VL,VU] */
/* > will be found;
*/
/* > = 'I': the IL-th through IU-th eigenvalues will be found. */
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
/* > reduction to tridiagonal form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD + 1. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDQ, N) */
/* > If JOBZ = 'V', the N-by-N unitary matrix used in the */
/* > reduction to tridiagonal form. */
/* > If JOBZ = 'N', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. If JOBZ = 'V', then */
/* > LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* > VL is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* > VU is REAL */
/* > If RANGE='V', the lower and upper bounds of the interval to */
/* > be searched for eigenvalues. VL < VU. */
/* > Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* > IL is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* > IU is INTEGER */
/* > If RANGE='I', the indices (in ascending order) of the */
/* > smallest and largest eigenvalues to be returned. */
/* > 1 <= IL <= IU <= N, if N > 0;
IL = 1 and IU = 0 if N = 0. */
/* > Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] ABSTOL */
/* > \verbatim */
/* > ABSTOL is REAL */
/* > The absolute error tolerance for the eigenvalues. */
/* > An approximate eigenvalue is accepted as converged */
/* > when it is determined to lie in an interval [a,b] */
/* > of width less than or equal to */
/* > */
/* > ABSTOL + EPS * max( |a|,|b| ) , */
/* > */
/* > where EPS is the machine precision. If ABSTOL is less than */
/* > or equal to zero, then EPS*|T| will be used in its place, */
/* > where |T| is the 1-norm of the tridiagonal matrix obtained */
/* > by reducing AB to tridiagonal form. */
/* > */
/* > Eigenvalues will be computed most accurately when ABSTOL is */
/* > set to twice the underflow threshold 2*SLAMCH('S'), not zero. */
/* > If this routine returns with INFO>0, indicating that some */
/* > eigenvectors did not converge, try setting ABSTOL to */
/* > 2*SLAMCH('S'). */
/* > */
/* > See "Computing Small Singular Values of Bidiagonal Matrices */
/* > with Guaranteed High Relative Accuracy," by Demmel and */
/* > Kahan, LAPACK Working Note #3. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The total number of eigenvalues found. 0 <= M <= N. */
/* > If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > The first M elements contain the selected eigenvalues in */
/* > ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ, max(1,M)) */
/* > If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* > contain the orthonormal eigenvectors of the matrix A */
/* > corresponding to the selected eigenvalues, with the i-th */
/* > column of Z holding the eigenvector associated with W(i). */
/* > If an eigenvector fails to converge, then that column of Z */
/* > contains the latest approximation to the eigenvector, and the */
/* > index of the eigenvector is returned in IFAIL. */
/* > If JOBZ = 'N', then Z is not referenced. */
/* > Note: the user must ensure that at least max(1,M) columns are */
/* > supplied in the array Z;
if RANGE = 'V', the exact value of M */
/* > is not known in advance and an upper bound must be used. */
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
/* > WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (7*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (5*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAIL */
/* > \verbatim */
/* > IFAIL is INTEGER array, dimension (N) */
/* > If JOBZ = 'V', then if INFO = 0, the first M elements of */
/* > IFAIL are zero. If INFO > 0, then IFAIL contains the */
/* > indices of the eigenvectors that failed to converge. */
/* > If JOBZ = 'N', then IFAIL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, then i eigenvectors failed to converge. */
/* > Their indices are stored in array IFAIL. */
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
int chbevx_(char *jobz, char *range, char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, complex *q, integer *ldq, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer * m, real *w, complex *z__, integer *ldz, complex *work, real *rwork, integer *iwork, integer *ifail, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"chbevx inputs: jobz %c, range %c, uplo %c, n %lld, kd %lld, ldab %lld, ldq %lld, il %lld, iu %lld, m %lld, ldz %lld",*jobz, *range, *uplo, *n, *kd, *ldab, *ldq, *il, *iu, *m, *ldz);
#else
    snprintf(buffer, 256,"chbevx inputs: jobz %c, range %c, uplo %c, n %d, kd %d, ldab %d, ldq %d, il %d, iu %d, m %d, ldz %d",*jobz, *range, *uplo, *n, *kd, *ldab, *ldq, *il, *iu, *m, *ldz);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer ab_dim1, ab_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, jj;
    real eps, vll, vuu, tmp1;
    integer indd, inde;
    real anrm;
    integer imax;
    real rmin, rmax;
    logical test;
    complex ctmp1;
    integer itmp1, indee;
    real sigma;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int cgemv_(char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    integer iinfo;
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    char order[1];
    extern /* Subroutine */
    int ccopy_(integer *, complex *, integer *, complex *, integer *), cswap_(integer *, complex *, integer *, complex *, integer *);
    logical lower;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *);
    logical wantz;
    extern real clanhb_(char *, char *, integer *, integer *, complex *, integer *, real *);
    logical alleig, indeig;
    integer iscale, indibl;
    extern /* Subroutine */
    int clascl_(char *, integer *, integer *, real *, real *, integer *, integer *, complex *, integer *, integer *), chbtrd_(char *, char *, integer *, integer *, complex *, integer *, real *, real *, complex *, integer *, complex *, integer *);
    logical valeig;
    extern real slamch_(char *);
    extern /* Subroutine */
    int clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *);
    real safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    real abstll, bignum;
    integer indiwk, indisp;
    extern /* Subroutine */
    int cstein_(integer *, real *, real *, integer *, real *, integer *, integer *, complex *, integer *, real *, integer *, integer *, integer *);
    integer indrwk, indwrk;
    extern /* Subroutine */
    int csteqr_(char *, integer *, real *, real *, complex *, integer *, real *, integer *), ssterf_(integer *, real *, real *, integer *);
    integer nsplit;
    extern /* Subroutine */
    int sstebz_(char *, char *, integer *, real *, real *, integer *, integer *, real *, real *, real *, integer *, integer *, real *, integer *, integer *, real *, integer *, integer *);
    real smlnum;
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
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    --iwork;
    --ifail;
    /* Function Body */
    wantz = lsame_(jobz, "V");
    alleig = lsame_(range, "A");
    valeig = lsame_(range, "V");
    indeig = lsame_(range, "I");
    lower = lsame_(uplo, "L");
    *info = 0;
    if (! (wantz || lsame_(jobz, "N")))
    {
        *info = -1;
    }
    else if (! (alleig || valeig || indeig))
    {
        *info = -2;
    }
    else if (! (lower || lsame_(uplo, "U")))
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*kd < 0)
    {
        *info = -5;
    }
    else if (*ldab < *kd + 1)
    {
        *info = -7;
    }
    else if (wantz && *ldq < max(1,*n))
    {
        *info = -9;
    }
    else
    {
        if (valeig)
        {
            if (*n > 0 && *vu <= *vl)
            {
                *info = -11;
            }
        }
        else if (indeig)
        {
            if (*il < 1 || *il > max(1,*n))
            {
                *info = -12;
            }
            else if (*iu < min(*n,*il) || *iu > *n)
            {
                *info = -13;
            }
        }
    }
    if (*info == 0)
    {
        if (*ldz < 1 || wantz && *ldz < *n)
        {
            *info = -18;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHBEVX", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    *m = 0;
    if (*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    if (*n == 1)
    {
        *m = 1;
        if (lower)
        {
            i__1 = ab_dim1 + 1;
            ctmp1.r = ab[i__1].r;
            ctmp1.i = ab[i__1].i; // , expr subst
        }
        else
        {
            i__1 = *kd + 1 + ab_dim1;
            ctmp1.r = ab[i__1].r;
            ctmp1.i = ab[i__1].i; // , expr subst
        }
        tmp1 = ctmp1.r;
        if (valeig)
        {
            if (! (*vl < tmp1 && *vu >= tmp1))
            {
                *m = 0;
            }
        }
        if (*m == 1)
        {
            w[1] = ctmp1.r;
            if (wantz)
            {
                i__1 = z_dim1 + 1;
                z__[i__1].r = 1.f;
                z__[i__1].i = 0.f; // , expr subst
            }
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
    /* Computing MIN */
    r__1 = sqrt(bignum);
    r__2 = 1.f / sqrt(sqrt(safmin)); // , expr subst
    rmax = min(r__1,r__2);
    /* Scale matrix to allowable range, if necessary. */
    iscale = 0;
    abstll = *abstol;
    if (valeig)
    {
        vll = *vl;
        vuu = *vu;
    }
    else
    {
        vll = 0.f;
        vuu = 0.f;
    }
    anrm = clanhb_("M", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1]);
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
            clascl_("B", kd, kd, &c_b16, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
        else
        {
            clascl_("Q", kd, kd, &c_b16, &sigma, n, n, &ab[ab_offset], ldab, info);
        }
        if (*abstol > 0.f)
        {
            abstll = *abstol * sigma;
        }
        if (valeig)
        {
            vll = *vl * sigma;
            vuu = *vu * sigma;
        }
    }
    /* Call CHBTRD to reduce Hermitian band matrix to tridiagonal form. */
    indd = 1;
    inde = indd + *n;
    indrwk = inde + *n;
    indwrk = 1;
    chbtrd_(jobz, uplo, n, kd, &ab[ab_offset], ldab, &rwork[indd], &rwork[ inde], &q[q_offset], ldq, &work[indwrk], &iinfo);
    /* If all eigenvalues are desired and ABSTOL is less than or equal */
    /* to zero, then call SSTERF or CSTEQR. If this fails for some */
    /* eigenvalue, then try SSTEBZ. */
    test = FALSE_;
    if (indeig)
    {
        if (*il == 1 && *iu == *n)
        {
            test = TRUE_;
        }
    }
    if ((alleig || test) && *abstol <= 0.f)
    {
        scopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
        indee = indrwk + (*n << 1);
        if (! wantz)
        {
            i__1 = *n - 1;
            scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
            ssterf_(n, &w[1], &rwork[indee], info);
        }
        else
        {
            clacpy_("A", n, n, &q[q_offset], ldq, &z__[z_offset], ldz);
            i__1 = *n - 1;
            scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
            csteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, & rwork[indrwk], info);
            if (*info == 0)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    ifail[i__] = 0;
                    /* L10: */
                }
            }
        }
        if (*info == 0)
        {
            *m = *n;
            goto L30;
        }
        *info = 0;
    }
    /* Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN. */
    if (wantz)
    {
        *(unsigned char *)order = 'B';
    }
    else
    {
        *(unsigned char *)order = 'E';
    }
    indibl = 1;
    indisp = indibl + *n;
    indiwk = indisp + *n;
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], & rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], & rwork[indrwk], &iwork[indiwk], info);
    if (wantz)
    {
        cstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], & iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[ indiwk], &ifail[1], info);
        /* Apply unitary matrix used in reduction to tridiagonal */
        /* form to eigenvectors returned by CSTEIN. */
        i__1 = *m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            ccopy_(n, &z__[j * z_dim1 + 1], &c__1, &work[1], &c__1);
            cgemv_("N", n, n, &c_b2, &q[q_offset], ldq, &work[1], &c__1, & c_b1, &z__[j * z_dim1 + 1], &c__1);
            /* L20: */
        }
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
L30:
    if (iscale == 1)
    {
        if (*info == 0)
        {
            imax = *m;
        }
        else
        {
            imax = *info - 1;
        }
        r__1 = 1.f / sigma;
        sscal_(&imax, &r__1, &w[1], &c__1);
    }
    /* If eigenvalues are not in order, then sort them, along with */
    /* eigenvectors. */
    if (wantz)
    {
        i__1 = *m - 1;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__ = 0;
            tmp1 = w[j];
            i__2 = *m;
            for (jj = j + 1;
                    jj <= i__2;
                    ++jj)
            {
                if (w[jj] < tmp1)
                {
                    i__ = jj;
                    tmp1 = w[jj];
                }
                /* L40: */
            }
            if (i__ != 0)
            {
                itmp1 = iwork[indibl + i__ - 1];
                w[i__] = w[j];
                iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
                w[j] = tmp1;
                iwork[indibl + j - 1] = itmp1;
                cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
                if (*info != 0)
                {
                    itmp1 = ifail[i__];
                    ifail[i__] = ifail[j];
                    ifail[j] = itmp1;
                }
            }
            /* L50: */
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CHBEVX */
}
/* chbevx_ */
