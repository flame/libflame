#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b17 = 1.;
/* > \brief <b> dsteqr_helper computes the eigenvalues and  right eigenvectors for SY matrices</b> */
/* =========== DOCUMENTATION =========== */
/* ===================================================================== */
/* Subroutine */
int dsteqr_helper_(char *jobz, char *uplo, integer *n, doublereal * a, integer *lda, doublereal *w, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dsteqr_helper inputs: jobz %c, uplo %c, n %" FLA_IS ", lda %" FLA_IS ", lwork %" FLA_IS ", liwork %" FLA_IS "",*jobz, *uplo, *n, *lda, *lwork, *liwork);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    doublereal eps;
    integer inde;
    doublereal anrm, rmin, rmax;
    integer lopt;
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *);
    integer iinfo, lwmin, liopt;
    logical lower, wantz;
    integer indwk2, llwrk2;
    extern doublereal dlamch_(char *);
    integer iscale;
    extern /* Subroutine */
    int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *), dstedc_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, integer *, integer *), dlacpy_( char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal bignum;
    integer indtau;
    extern /* Subroutine */
    int dsterf_(integer *, doublereal *, doublereal *, integer *);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, integer *, doublereal *);
    integer indwrk, liwmin;
    extern /* Subroutine */
    int dormtr_(char *, char *, char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
    integer llwork;
    doublereal smlnum;
    logical lquery;
    /* -- LAPACK driver routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    --w;
    --work;
    --iwork;
    /* Function Body */
    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");
    lquery = *lwork == -1 || *liwork == -1;
    *info = 0;

    if (*info == 0)
    {
        if (*n <= 1)
        {
            liwmin = 1;
            lwmin = 1;
            lopt = lwmin;
            liopt = liwmin;
        }
        else
        {
            if (wantz)
            {
                liwmin = *n * 5 + 3;
                /* Computing 2nd power */
                i__1 = *n;
                lwmin = *n * 6 + 1 + (i__1 * i__1 << 1);
            }
            else
            {
                liwmin = 1;
                lwmin = (*n << 1) + 1;
            }
            /* Computing MAX */
            i__1 = lwmin;
            i__2 = (*n << 1) + ilaenv_(&c__1, "DSYEVD", uplo, n, &c_n1, &c_n1, &c_n1); // , expr subst
            lopt = max(i__1,i__2);
            liopt = liwmin;
        }
        work[1] = (doublereal) lopt;
        iwork[1] = liopt;
        if (*lwork < lwmin && ! lquery)
        {
            *info = -8;
        }
        else if (*liwork < liwmin && ! lquery)
        {
            *info = -10;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DSYEVD", &i__1);
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
        w[1] = a[a_dim1 + 1];
        if (wantz)
        {
            a[a_dim1 + 1] = 1.;
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
    anrm = dlansy_("M", uplo, n, &a[a_offset], lda, &work[1]);
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
        dlascl_(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, info);
    }
    // printf("reaching after lascl\n");
    /* Call DSYTRD to reduce symmetric matrix to tridiagonal form. */
    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    indwk2 = indwrk + *n * *n;
    llwrk2 = *lwork - indwk2 + 1;
    /* For eigenvalues only, call DSTERF. For eigenvectors, first call */
    /* DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
    /* tridiagonal matrix, then call DORMTR to multiply it by the */
    /* Householder transformations stored in A. */
    dsytrd_(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &
            work[indwrk], &llwork, &iinfo);
    lopt = (integer) ((*n << 1) + work[indwrk]);
    if (! wantz)
    {
        dsterf_(n, &w[1], &work[inde], info);
    }
    else
    {
        dstedc_("I", n, &w[1], &work[inde], &work[indwrk], n, &work[indwk2], & llwrk2, &iwork[1], liwork, info);
        dormtr_("L", uplo, "N", n, n, &a[a_offset], lda, &work[indtau], &work[ indwrk], n, &work[indwk2], &llwrk2, &iinfo);
        dlacpy_("A", n, n, &work[indwrk], n, &a[a_offset], lda);
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
    if (iscale == 1)
    {
        d__1 = 1. / sigma;
        dscal_(n, &d__1, &w[1], &c__1);
    }
    work[1] = (doublereal) lopt;
    iwork[1] = liopt;

    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of DSTEQR_HELPER */
}
/* dsteqr_helper_ */
