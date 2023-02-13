/* ../netlib/dstevx.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief <b> DSTEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSTEVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstevx. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstevx. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstevx. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, */
/* M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBZ, RANGE */
/* INTEGER IL, INFO, IU, LDZ, M, N */
/* DOUBLE PRECISION ABSTOL, VL, VU */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IFAIL( * ), IWORK( * ) */
/* DOUBLE PRECISION D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEVX computes selected eigenvalues and, optionally, eigenvectors */
/* > of a real symmetric tridiagonal matrix A. Eigenvalues and */
/* > eigenvectors can be selected by specifying either a range of values */
/* > or a range of indices for the desired eigenvalues. */
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
/* > = 'A': all eigenvalues will be found. */
/* > = 'V': all eigenvalues in the half-open interval (VL,VU] */
/* > will be found. */
/* > = 'I': the IL-th through IU-th eigenvalues will be found. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > On entry, the n diagonal elements of the tridiagonal matrix */
/* > A. */
/* > On exit, D may be multiplied by a constant factor chosen */
/* > to avoid over/underflow in computing the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (max(1,N-1)) */
/* > On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* > matrix A in elements 1 to N-1 of E. */
/* > On exit, E may be multiplied by a constant factor chosen */
/* > to avoid over/underflow in computing the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* > VL is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* > VU is DOUBLE PRECISION */
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
/* > ABSTOL is DOUBLE PRECISION */
/* > The absolute error tolerance for the eigenvalues. */
/* > An approximate eigenvalue is accepted as converged */
/* > when it is determined to lie in an interval [a,b] */
/* > of width less than or equal to */
/* > */
/* > ABSTOL + EPS * max( |a|,|b| ) , */
/* > */
/* > where EPS is the machine precision. If ABSTOL is less */
/* > than or equal to zero, then EPS*|T| will be used in */
/* > its place, where |T| is the 1-norm of the tridiagonal */
/* > matrix. */
/* > */
/* > Eigenvalues will be computed most accurately when ABSTOL is */
/* > set to twice the underflow threshold 2*DLAMCH('S'), not zero. */
/* > If this routine returns with INFO>0, indicating that some */
/* > eigenvectors did not converge, try setting ABSTOL to */
/* > 2*DLAMCH('S'). */
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
/* > W is DOUBLE PRECISION array, dimension (N) */
/* > The first M elements contain the selected eigenvalues in */
/* > ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M) ) */
/* > If JOBZ = 'V', then if INFO = 0, the first M columns of Z */
/* > contain the orthonormal eigenvectors of the matrix A */
/* > corresponding to the selected eigenvalues, with the i-th */
/* > column of Z holding the eigenvector associated with W(i). */
/* > If an eigenvector fails to converge (INFO > 0), then that */
/* > column of Z contains the latest approximation to the */
/* > eigenvector, and the index of the eigenvector is returned */
/* > in IFAIL. If JOBZ = 'N', then Z is not referenced. */
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
/* > WORK is DOUBLE PRECISION array, dimension (5*N) */
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
/* > \ingroup doubleOTHEReigen */
/* ===================================================================== */
/* Subroutine */
int dstevx_(char *jobz, char *range, integer *n, doublereal * d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work, integer *iwork, integer *ifail, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, jj;
    doublereal eps, vll, vuu, tmp1;
    integer imax;
    doublereal rmin, rmax;
    logical test;
    doublereal tnrm;
    integer itmp1;
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *);
    char order[1];
    extern /* Subroutine */
    int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *), dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical wantz;
    extern doublereal dlamch_(char *);
    logical alleig, indeig;
    integer iscale, indibl;
    logical valeig;
    doublereal safmin;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    doublereal bignum;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *);
    integer indisp;
    extern /* Subroutine */
    int dstein_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, integer *), dsterf_(integer *, doublereal *, doublereal *, integer *);
    integer indiwo;
    extern /* Subroutine */
    int dstebz_(char *, char *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
    integer indwrk;
    extern /* Subroutine */
    int dsteqr_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    integer nsplit;
    doublereal smlnum;
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
    --d__;
    --e;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;
    /* Function Body */
    wantz = lsame_(jobz, "V");
    alleig = lsame_(range, "A");
    valeig = lsame_(range, "V");
    indeig = lsame_(range, "I");
    *info = 0;
    if (! (wantz || lsame_(jobz, "N")))
    {
        *info = -1;
    }
    else if (! (alleig || valeig || indeig))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else
    {
        if (valeig)
        {
            if (*n > 0 && *vu <= *vl)
            {
                *info = -7;
            }
        }
        else if (indeig)
        {
            if (*il < 1 || *il > max(1,*n))
            {
                *info = -8;
            }
            else if (*iu < min(*n,*il) || *iu > *n)
            {
                *info = -9;
            }
        }
    }
    if (*info == 0)
    {
        if (*ldz < 1 || wantz && *ldz < *n)
        {
            *info = -14;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DSTEVX", &i__1);
        return 0;
    }
    /* Quick return if possible */
    *m = 0;
    if (*n == 0)
    {
        return 0;
    }
    if (*n == 1)
    {
        if (alleig || indeig)
        {
            *m = 1;
            w[1] = d__[1];
        }
        else
        {
            if (*vl < d__[1] && *vu >= d__[1])
            {
                *m = 1;
                w[1] = d__[1];
            }
        }
        if (wantz)
        {
            z__[z_dim1 + 1] = 1.;
        }
        return 0;
    }
    /* Get machine constants. */
    safmin = dlamch_("Safe minimum");
    eps = dlamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    /* Computing MIN */
    d__1 = sqrt(bignum);
    d__2 = 1. / sqrt(sqrt(safmin)); // , expr subst
    rmax = min(d__1,d__2);
    /* Scale matrix to allowable range, if necessary. */
    iscale = 0;
    if (valeig)
    {
        vll = *vl;
        vuu = *vu;
    }
    else
    {
        vll = 0.;
        vuu = 0.;
    }
    tnrm = dlanst_("M", n, &d__[1], &e[1]);
    if (tnrm > 0. && tnrm < rmin)
    {
        iscale = 1;
        sigma = rmin / tnrm;
    }
    else if (tnrm > rmax)
    {
        iscale = 1;
        sigma = rmax / tnrm;
    }
    if (iscale == 1)
    {
        dscal_(n, &sigma, &d__[1], &c__1);
        i__1 = *n - 1;
        dscal_(&i__1, &sigma, &e[1], &c__1);
        if (valeig)
        {
            vll = *vl * sigma;
            vuu = *vu * sigma;
        }
    }
    /* If all eigenvalues are desired and ABSTOL is less than zero, then */
    /* call DSTERF or SSTEQR. If this fails for some eigenvalue, then */
    /* try DSTEBZ. */
    test = FALSE_;
    if (indeig)
    {
        if (*il == 1 && *iu == *n)
        {
            test = TRUE_;
        }
    }
    if ((alleig || test) && *abstol <= 0.)
    {
        dcopy_(n, &d__[1], &c__1, &w[1], &c__1);
        i__1 = *n - 1;
        dcopy_(&i__1, &e[1], &c__1, &work[1], &c__1);
        indwrk = *n + 1;
        if (! wantz)
        {
            dsterf_(n, &w[1], &work[1], info);
        }
        else
        {
            dsteqr_("I", n, &w[1], &work[1], &z__[z_offset], ldz, &work[ indwrk], info);
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
            goto L20;
        }
        *info = 0;
    }
    /* Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */
    if (wantz)
    {
        *(unsigned char *)order = 'B';
    }
    else
    {
        *(unsigned char *)order = 'E';
    }
    indwrk = 1;
    indibl = 1;
    indisp = indibl + *n;
    indiwo = indisp + *n;
    dstebz_(range, order, n, &vll, &vuu, il, iu, abstol, &d__[1], &e[1], m, & nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[indwrk], & iwork[indiwo], info);
    if (wantz)
    {
        dstein_(n, &d__[1], &e[1], m, &w[1], &iwork[indibl], &iwork[indisp], & z__[z_offset], ldz, &work[indwrk], &iwork[indiwo], &ifail[1], info);
    }
    /* If matrix was scaled, then rescale eigenvalues appropriately. */
L20:
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
        d__1 = 1. / sigma;
        dscal_(&imax, &d__1, &w[1], &c__1);
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
                /* L30: */
            }
            if (i__ != 0)
            {
                itmp1 = iwork[indibl + i__ - 1];
                w[i__] = w[j];
                iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
                w[j] = tmp1;
                iwork[indibl + j - 1] = itmp1;
                dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
                if (*info != 0)
                {
                    itmp1 = ifail[i__];
                    ifail[i__] = ifail[j];
                    ifail[j] = itmp1;
                }
            }
            /* L40: */
        }
    }
    return 0;
    /* End of DSTEVX */
}
/* dstevx_ */
