/* ../netlib/cstein.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b CSTEIN */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSTEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cstein. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cstein. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cstein. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, */
/* IWORK, IFAIL, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDZ, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IBLOCK( * ), IFAIL( * ), ISPLIT( * ), */
/* $ IWORK( * ) */
/* REAL D( * ), E( * ), W( * ), WORK( * ) */
/* COMPLEX Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSTEIN computes the eigenvectors of a real symmetric tridiagonal */
/* > matrix T corresponding to specified eigenvalues, using inverse */
/* > iteration. */
/* > */
/* > The maximum number of iterations allowed for each eigenvector is */
/* > specified by an internal parameter MAXITS (currently set to 5). */
/* > */
/* > Although the eigenvectors are real, they are stored in a complex */
/* > array, which may be passed to CUNMTR or CUPMTR for back */
/* > transformation to the eigenvectors of a complex Hermitian matrix */
/* > which was reduced to tridiagonal form. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The n diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > The (n-1) subdiagonal elements of the tridiagonal matrix */
/* > T, stored in elements 1 to N-1. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of eigenvectors to be found. 0 <= M <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > The first M elements of W contain the eigenvalues for */
/* > which eigenvectors are to be computed. The eigenvalues */
/* > should be grouped by split-off block and ordered from */
/* > smallest to largest within the block. ( The output array */
/* > W from SSTEBZ with ORDER = 'B' is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[in] IBLOCK */
/* > \verbatim */
/* > IBLOCK is INTEGER array, dimension (N) */
/* > The submatrix indices associated with the corresponding */
/* > eigenvalues in W;
IBLOCK(i)=1 if eigenvalue W(i) belongs to */
/* > the first submatrix from the top, =2 if W(i) belongs to */
/* > the second submatrix, etc. ( The output array IBLOCK */
/* > from SSTEBZ is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[in] ISPLIT */
/* > \verbatim */
/* > ISPLIT is INTEGER array, dimension (N) */
/* > The splitting points, at which T breaks up into submatrices. */
/* > The first submatrix consists of rows/columns 1 to */
/* > ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1 */
/* > through ISPLIT( 2 ), etc. */
/* > ( The output array ISPLIT from SSTEBZ is expected here. ) */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ, M) */
/* > The computed eigenvectors. The eigenvector associated */
/* > with the eigenvalue W(i) is stored in the i-th column of */
/* > Z. Any vector which fails to converge is set to its current */
/* > iterate after MAXITS iterations. */
/* > The imaginary parts of the eigenvectors are set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (5*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAIL */
/* > \verbatim */
/* > IFAIL is INTEGER array, dimension (M) */
/* > On normal exit, all elements of IFAIL are zero. */
/* > If one or more eigenvectors fail to converge after */
/* > MAXITS iterations, then their indices are stored in */
/* > array IFAIL. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, then i eigenvectors failed to converge */
/* > in MAXITS iterations. Their indices are stored in */
/* > array IFAIL. */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > MAXITS INTEGER, default = 5 */
/* > The maximum number of iterations performed. */
/* > */
/* > EXTRA INTEGER, default = 2 */
/* > The number of iterations performed after norm growth */
/* > criterion is satisfied, should be at least 1. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int cstein_(integer *n, real *d__, real *e, integer *m, real *w, integer *iblock, integer *isplit, complex *z__, integer *ldz, real *work, integer *iwork, integer *ifail, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"cstein inputs: n %lld, m %lld, iblock %lld, isplit %lld, ldz %lld",*n, *m, *iblock, *isplit, *ldz);
#else
    snprintf(buffer, 256,"cstein inputs: n %d, m %d, iblock %d, isplit %d, ldz %d",*n, *m, *iblock, *isplit, *ldz);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5;
    complex q__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, b1, j1, bn, jr;
    real xj, scl, eps, ctr, sep, nrm, tol;
    integer its;
    real xjm, eps1;
    integer jblk, nblk, jmax;
    extern real snrm2_(integer *, real *, integer *);
    integer iseed[4], gpind, iinfo;
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *), scopy_(integer *, real *, integer *, real *, integer *);
    real ortol;
    integer indrv1, indrv2, indrv3, indrv4, indrv5;
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len), slagtf_( integer *, real *, real *, real *, real *, real *, real *, integer *, integer *);
    integer nrmchk;
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */
    int slagts_(integer *, integer *, real *, real *, real *, real *, integer *, real *, real *, integer *);
    integer blksiz;
    real onenrm, pertol;
    extern /* Subroutine */
    int slarnv_(integer *, integer *, integer *, real *);
    real stpcrt;
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Local Arrays .. */
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
    --iblock;
    --isplit;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;
    /* Function Body */
    *info = 0;
    stpcrt = 0.f;
    onenrm = 0.f;
    ortol = 0.f;
    gpind= 0;
    xjm = 0.f;
    i__1 = *m;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        ifail[i__] = 0;
        /* L10: */
    }
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*m < 0 || *m > *n)
    {
        *info = -4;
    }
    else if (*ldz < fla_max(1,*n))
    {
        *info = -9;
    }
    else
    {
        i__1 = *m;
        for (j = 2;
                j <= i__1;
                ++j)
        {
            if (iblock[j] < iblock[j - 1])
            {
                *info = -6;
                goto L30;
            }
            if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1])
            {
                *info = -5;
                goto L30;
            }
            /* L20: */
        }
L30:
        ;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CSTEIN", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *m == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (*n == 1)
    {
        i__1 = z_dim1 + 1;
        z__[i__1].r = 1.f;
        z__[i__1].i = 0.f; // , expr subst
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Get machine constants. */
    eps = slamch_("Precision");
    /* Initialize seed for random number generator SLARNV. */
    for (i__ = 1;
            i__ <= 4;
            ++i__)
    {
        iseed[i__ - 1] = 1;
        /* L40: */
    }
    /* Initialize pointers. */
    indrv1 = 0;
    indrv2 = indrv1 + *n;
    indrv3 = indrv2 + *n;
    indrv4 = indrv3 + *n;
    indrv5 = indrv4 + *n;
    /* Compute eigenvectors of matrix blocks. */
    j1 = 1;
    i__1 = iblock[*m];
    for (nblk = 1;
            nblk <= i__1;
            ++nblk)
    {
        /* Find starting and ending indices of block nblk. */
        if (nblk == 1)
        {
            b1 = 1;
        }
        else
        {
            b1 = isplit[nblk - 1] + 1;
        }
        bn = isplit[nblk];
        blksiz = bn - b1 + 1;
        if (blksiz == 1)
        {
            goto L60;
        }
        gpind = j1;
        /* Compute reorthogonalization criterion and stopping criterion. */
        onenrm = (r__1 = d__[b1], f2c_abs(r__1)) + (r__2 = e[b1], f2c_abs(r__2));
        /* Computing MAX */
        r__3 = onenrm;
        r__4 = (r__1 = d__[bn], f2c_abs(r__1)) + (r__2 = e[bn - 1], f2c_abs(r__2)); // , expr subst
        onenrm = fla_max(r__3,r__4);
        i__2 = bn - 1;
        for (i__ = b1 + 1;
                i__ <= i__2;
                ++i__)
        {
            /* Computing MAX */
            r__4 = onenrm;
            r__5 = (r__1 = d__[i__], f2c_abs(r__1)) + (r__2 = e[ i__ - 1], f2c_abs(r__2)) + (r__3 = e[i__], f2c_abs(r__3)); // , expr subst
            onenrm = fla_max(r__4,r__5);
            /* L50: */
        }
        ortol = onenrm * .001f;
        stpcrt = sqrt(.1f / blksiz);
        /* Loop through eigenvalues of block nblk. */
L60:
        jblk = 0;
        i__2 = *m;
        for (j = j1;
                j <= i__2;
                ++j)
        {
            if (iblock[j] != nblk)
            {
                j1 = j;
                goto L180;
            }
            ++jblk;
            xj = w[j];
            /* Skip all the work if the block size is one. */
            if (blksiz == 1)
            {
                work[indrv1 + 1] = 1.f;
                goto L140;
            }
            /* If eigenvalues j and j-1 are too close, add a relatively */
            /* small perturbation. */
            if (jblk > 1)
            {
                eps1 = (r__1 = eps * xj, f2c_abs(r__1));
                pertol = eps1 * 10.f;
                sep = xj - xjm;
                if (sep < pertol)
                {
                    xj = xjm + pertol;
                }
            }
            its = 0;
            nrmchk = 0;
            /* Get random starting vector. */
            slarnv_(&c__2, iseed, &blksiz, &work[indrv1 + 1]);
            /* Copy the matrix T so it won't be destroyed in factorization. */
            scopy_(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
            i__3 = blksiz - 1;
            scopy_(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
            i__3 = blksiz - 1;
            scopy_(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);
            /* Compute LU factors with partial pivoting ( PT = LU ) */
            tol = 0.f;
            slagtf_(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[ indrv3 + 1], &tol, &work[indrv5 + 1], &iwork[1], &iinfo);
            /* Update iteration count. */
L70:
            ++its;
            if (its > 5)
            {
                goto L120;
            }
            /* Normalize and scale the righthand side vector Pb. */
            jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
            /* Computing MAX */
            r__3 = eps;
            r__4 = (r__1 = work[indrv4 + blksiz], f2c_abs(r__1)); // , expr subst
            scl = blksiz * onenrm * fla_max(r__3,r__4) / (r__2 = work[indrv1 + jmax], f2c_abs(r__2));
            sscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
            /* Solve the system LU = Pb. */
            slagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], & work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[ indrv1 + 1], &tol, &iinfo);
            /* Reorthogonalize by modified Gram-Schmidt if eigenvalues are */
            /* close enough. */
            if (jblk == 1)
            {
                goto L110;
            }
            if ((r__1 = xj - xjm, f2c_abs(r__1)) > ortol)
            {
                gpind = j;
            }
            if (gpind != j)
            {
                i__3 = j - 1;
                for (i__ = gpind;
                        i__ <= i__3;
                        ++i__)
                {
                    ctr = 0.f;
                    i__4 = blksiz;
                    for (jr = 1;
                            jr <= i__4;
                            ++jr)
                    {
                        i__5 = b1 - 1 + jr + i__ * z_dim1;
                        ctr += work[indrv1 + jr] * z__[i__5].r;
                        /* L80: */
                    }
                    i__4 = blksiz;
                    for (jr = 1;
                            jr <= i__4;
                            ++jr)
                    {
                        i__5 = b1 - 1 + jr + i__ * z_dim1;
                        work[indrv1 + jr] -= ctr * z__[i__5].r;
                        /* L90: */
                    }
                    /* L100: */
                }
            }
            /* Check the infinity norm of the iterate. */
L110:
            jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
            nrm = (r__1 = work[indrv1 + jmax], f2c_abs(r__1));
            /* Continue for additional iterations after norm reaches */
            /* stopping criterion. */
            if (nrm < stpcrt)
            {
                goto L70;
            }
            ++nrmchk;
            if (nrmchk < 3)
            {
                goto L70;
            }
            goto L130;
            /* If stopping criterion was not satisfied, update info and */
            /* store eigenvector number in array ifail. */
L120:
            ++(*info);
            ifail[*info] = j;
            /* Accept iterate as jth eigenvector. */
L130:
            scl = 1.f / snrm2_(&blksiz, &work[indrv1 + 1], &c__1);
            jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
            if (work[indrv1 + jmax] < 0.f)
            {
                scl = -scl;
            }
            sscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
L140:
            i__3 = *n;
            for (i__ = 1;
                    i__ <= i__3;
                    ++i__)
            {
                i__4 = i__ + j * z_dim1;
                z__[i__4].r = 0.f;
                z__[i__4].i = 0.f; // , expr subst
                /* L150: */
            }
            i__3 = blksiz;
            for (i__ = 1;
                    i__ <= i__3;
                    ++i__)
            {
                i__4 = b1 + i__ - 1 + j * z_dim1;
                i__5 = indrv1 + i__;
                q__1.r = work[i__5];
                q__1.i = 0.f; // , expr subst
                z__[i__4].r = q__1.r;
                z__[i__4].i = q__1.i; // , expr subst
                /* L160: */
            }
            /* Save the shift to check eigenvalue spacing at next */
            /* iteration. */
            xjm = xj;
            /* L170: */
        }
L180:
        ;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CSTEIN */
}
/* cstein_ */
