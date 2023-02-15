/* ../netlib/sstein.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b SSTEIN */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSTEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstein. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstein. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstein. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, */
/* IWORK, IFAIL, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDZ, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IBLOCK( * ), IFAIL( * ), ISPLIT( * ), */
/* $ IWORK( * ) */
/* REAL D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEIN computes the eigenvectors of a real symmetric tridiagonal */
/* > matrix T corresponding to specified eigenvalues, using inverse */
/* > iteration. */
/* > */
/* > The maximum number of iterations allowed for each eigenvector is */
/* > specified by an internal parameter MAXITS (currently set to 5). */
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
/* > T, in elements 1 to N-1. */
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
/* > Z is REAL array, dimension (LDZ, M) */
/* > The computed eigenvectors. The eigenvector associated */
/* > with the eigenvalue W(i) is stored in the i-th column of */
/* > Z. Any vector which fails to converge is set to its current */
/* > iterate after MAXITS iterations. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= max(1,N). */
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
/* > = 0: successful exit. */
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
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int sstein_(integer *n, real *d__, real *e, integer *m, real *w, integer *iblock, integer *isplit, real *z__, integer *ldz, real * work, integer *iwork, integer *ifail, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4, r__5;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, b1, j1, bn;
    real xj, scl, eps, ctr, sep, nrm, tol;
    integer its;
    real xjm, eps1;
    integer jblk, nblk, jmax;
    extern real sdot_(integer *, real *, integer *, real *, integer *), snrm2_(integer *, real *, integer *);
    integer iseed[4], gpind, iinfo;
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    extern real sasum_(integer *, real *, integer *);
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *);
    real ortol;
    extern /* Subroutine */
    int saxpy_(integer *, real *, real *, integer *, real *, integer *);
    integer indrv1, indrv2, indrv3, indrv4, indrv5;
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *), slagtf_( integer *, real *, real *, real *, real *, real *, real *, integer *, integer *);
    integer nrmchk;
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */
    int slagts_(integer *, integer *, real *, real *, real *, real *, integer *, real *, real *, integer *);
    integer blksiz;
    real onenrm, pertol;
    extern /* Subroutine */
    int slarnv_(integer *, integer *, integer *, real *);
    real stpcrt;
    /* -- LAPACK computational routine (version 3.4.0) -- */
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
    else if (*ldz < max(1,*n))
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
        xerbla_("SSTEIN", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *m == 0)
    {
        return 0;
    }
    else if (*n == 1)
    {
        z__[z_dim1 + 1] = 1.f;
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
        gpind = b1;
        /* Compute reorthogonalization criterion and stopping criterion. */
        onenrm = (r__1 = d__[b1], f2c_abs(r__1)) + (r__2 = e[b1], f2c_abs(r__2));
        /* Computing MAX */
        r__3 = onenrm;
        r__4 = (r__1 = d__[bn], f2c_abs(r__1)) + (r__2 = e[bn - 1], f2c_abs(r__2)); // , expr subst
        onenrm = max(r__3,r__4);
        i__2 = bn - 1;
        for (i__ = b1 + 1;
                i__ <= i__2;
                ++i__)
        {
            /* Computing MAX */
            r__4 = onenrm;
            r__5 = (r__1 = d__[i__], f2c_abs(r__1)) + (r__2 = e[ i__ - 1], f2c_abs(r__2)) + (r__3 = e[i__], f2c_abs(r__3)); // , expr subst
            onenrm = max(r__4,r__5);
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
                goto L160;
            }
            ++jblk;
            xj = w[j];
            /* Skip all the work if the block size is one. */
            if (blksiz == 1)
            {
                work[indrv1 + 1] = 1.f;
                goto L120;
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
                goto L100;
            }
            /* Normalize and scale the righthand side vector Pb. */
            /* Computing MAX */
            r__2 = eps;
            r__3 = (r__1 = work[indrv4 + blksiz], f2c_abs(r__1)); // , expr subst
            scl = blksiz * onenrm * max(r__2,r__3) / sasum_(&blksiz, &work[ indrv1 + 1], &c__1);
            sscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
            /* Solve the system LU = Pb. */
            slagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], & work[indrv3 + 1], &work[indrv5 + 1], &iwork[1], &work[ indrv1 + 1], &tol, &iinfo);
            /* Reorthogonalize by modified Gram-Schmidt if eigenvalues are */
            /* close enough. */
            if (jblk == 1)
            {
                goto L90;
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
                    ctr = -sdot_(&blksiz, &work[indrv1 + 1], &c__1, &z__[b1 + i__ * z_dim1], &c__1);
                    saxpy_(&blksiz, &ctr, &z__[b1 + i__ * z_dim1], &c__1, & work[indrv1 + 1], &c__1);
                    /* L80: */
                }
            }
            /* Check the infinity norm of the iterate. */
L90:
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
            goto L110;
            /* If stopping criterion was not satisfied, update info and */
            /* store eigenvector number in array ifail. */
L100:
            ++(*info);
            ifail[*info] = j;
            /* Accept iterate as jth eigenvector. */
L110:
            scl = 1.f / snrm2_(&blksiz, &work[indrv1 + 1], &c__1);
            jmax = isamax_(&blksiz, &work[indrv1 + 1], &c__1);
            if (work[indrv1 + jmax] < 0.f)
            {
                scl = -scl;
            }
            sscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
L120:
            i__3 = *n;
            for (i__ = 1;
                    i__ <= i__3;
                    ++i__)
            {
                z__[i__ + j * z_dim1] = 0.f;
                /* L130: */
            }
            i__3 = blksiz;
            for (i__ = 1;
                    i__ <= i__3;
                    ++i__)
            {
                z__[b1 + i__ - 1 + j * z_dim1] = work[indrv1 + i__];
                /* L140: */
            }
            /* Save the shift to check eigenvalue spacing at next */
            /* iteration. */
            xjm = xj;
            /* L150: */
        }
L160:
        ;
    }
    return 0;
    /* End of SSTEIN */
}
/* sstein_ */
