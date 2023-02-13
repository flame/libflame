/* ../netlib/sstedc.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static real c_b17 = 0.f;
static real c_b18 = 1.f;
static integer c__1 = 1;
/* > \brief \b SSTEBZ */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSTEDC + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstedc. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstedc. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstedc. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, */
/* LIWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER COMPZ */
/* INTEGER INFO, LDZ, LIWORK, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSTEDC computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the divide and conquer method. */
/* > The eigenvectors of a full or band real symmetric matrix can also be */
/* > found if SSYTRD or SSPTRD or SSBTRD has been used to reduce this */
/* > matrix to tridiagonal form. */
/* > */
/* > This code makes very mild assumptions about floating point */
/* > arithmetic. It will work on machines with a guard digit in */
/* > add/subtract, or on those binary machines without guard digits */
/* > which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. */
/* > It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. See SLAED3 for details. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] COMPZ */
/* > \verbatim */
/* > COMPZ is CHARACTER*1 */
/* > = 'N': Compute eigenvalues only. */
/* > = 'I': Compute eigenvectors of tridiagonal matrix also. */
/* > = 'V': Compute eigenvectors of original dense symmetric */
/* > matrix also. On entry, Z contains the orthogonal */
/* > matrix used to reduce the original matrix to */
/* > tridiagonal form. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The dimension of the symmetric tridiagonal matrix. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, the diagonal elements of the tridiagonal matrix. */
/* > On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > On entry, the subdiagonal elements of the tridiagonal matrix. */
/* > On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ,N) */
/* > On entry, if COMPZ = 'V', then Z contains the orthogonal */
/* > matrix used in the reduction to tridiagonal form. */
/* > On exit, if INFO = 0, then if COMPZ = 'V', Z contains the */
/* > orthonormal eigenvectors of the original symmetric matrix, */
/* > and if COMPZ = 'I', Z contains the orthonormal eigenvectors */
/* > of the symmetric tridiagonal matrix. */
/* > If COMPZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. LDZ >= 1. */
/* > If eigenvectors are desired, then LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If COMPZ = 'N' or N <= 1 then LWORK must be at least 1. */
/* > If COMPZ = 'V' and N > 1 then LWORK must be at least */
/* > ( 1 + 3*N + 2*N*lg N + 4*N**2 ), */
/* > where lg( N ) = smallest integer k such */
/* > that 2**k >= N. */
/* > If COMPZ = 'I' and N > 1 then LWORK must be at least */
/* > ( 1 + 4*N + N**2 ). */
/* > Note that for COMPZ = 'I' or 'V', then if N is less than or */
/* > equal to the minimum divide size, usually 25, then LWORK need */
/* > only be max(1,2*(N-1)). */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
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
/* > The dimension of the array IWORK. */
/* > If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1. */
/* > If COMPZ = 'V' and N > 1 then LIWORK must be at least */
/* > ( 6 + 6*N + 5*N*lg N ). */
/* > If COMPZ = 'I' and N > 1 then LIWORK must be at least */
/* > ( 3 + 5*N ). */
/* > Note that for COMPZ = 'I' or 'V', then if N is less than or */
/* > equal to the minimum divide size, usually 25, then LIWORK */
/* > need only be 1. */
/* > */
/* > If LIWORK = -1, then a workspace query is assumed;
the */
/* > routine only calculates the optimal size of the IWORK array, */
/* > returns this value as the first entry of the IWORK array, and */
/* > no error message related to LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: The algorithm failed to compute an eigenvalue while */
/* > working on the submatrix lying in rows and columns */
/* > INFO/(N+1) through mod(INFO,N+1). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup auxOTHERcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA \n */
/* > Modified by Francoise Tisseur, University of Tennessee */
/* > */
/* ===================================================================== */
/* Subroutine */
int sstedc_(char *compz, integer *n, real *d__, real *e, real *z__, integer *ldz, real *work, integer *lwork, integer *iwork, integer *liwork, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k, m;
    real p;
    integer ii, lgn;
    real eps, tiny;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    integer lwmin, start;
    extern /* Subroutine */
    int sswap_(integer *, real *, integer *, real *, integer *), slaed0_(integer *, integer *, integer *, real *, real *, real *, integer *, real *, integer *, real *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer finish;
    extern /* Subroutine */
    int slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *, integer *, integer *), slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *), slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
    integer liwmin, icompz;
    real orgnrm;
    extern real slanst_(char *, integer *, real *, real *);
    extern /* Subroutine */
    int ssterf_(integer *, real *, real *, integer *), slasrt_(char *, integer *, real *, integer *);
    logical lquery;
    integer smlsiz;
    extern /* Subroutine */
    int ssteqr_(char *, integer *, real *, real *, real *, integer *, real *, integer *);
    integer storez, strtrw;
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
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1 || *liwork == -1;
    if (lsame_(compz, "N"))
    {
        icompz = 0;
    }
    else if (lsame_(compz, "V"))
    {
        icompz = 1;
    }
    else if (lsame_(compz, "I"))
    {
        icompz = 2;
    }
    else
    {
        icompz = -1;
    }
    if (icompz < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n))
    {
        *info = -6;
    }
    if (*info == 0)
    {
        /* Compute the workspace requirements */
        smlsiz = ilaenv_(&c__9, "SSTEDC", " ", &c__0, &c__0, &c__0, &c__0);
        if (*n <= 1 || icompz == 0)
        {
            liwmin = 1;
            lwmin = 1;
        }
        else if (*n <= smlsiz)
        {
            liwmin = 1;
            lwmin = *n - 1 << 1;
        }
        else
        {
            lgn = (integer) (log((real) (*n)) / log(2.f));
            if (pow_ii(&c__2, &lgn) < *n)
            {
                ++lgn;
            }
            if (pow_ii(&c__2, &lgn) < *n)
            {
                ++lgn;
            }
            if (icompz == 1)
            {
                /* Computing 2nd power */
                i__1 = *n;
                lwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 2);
                liwmin = *n * 6 + 6 + *n * 5 * lgn;
            }
            else if (icompz == 2)
            {
                /* Computing 2nd power */
                i__1 = *n;
                lwmin = (*n << 2) + 1 + i__1 * i__1;
                liwmin = *n * 5 + 3;
            }
        }
        work[1] = (real) lwmin;
        iwork[1] = liwmin;
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
        xerbla_("SSTEDC", &i__1);
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
        if (icompz != 0)
        {
            z__[z_dim1 + 1] = 1.f;
        }
        return 0;
    }
    /* If the following conditional clause is removed, then the routine */
    /* will use the Divide and Conquer routine to compute only the */
    /* eigenvalues, which requires (3N + 3N**2) real workspace and */
    /* (2 + 5N + 2N lg(N)) integer workspace. */
    /* Since on many architectures SSTERF is much faster than any other */
    /* algorithm for finding eigenvalues only, it is used here */
    /* as the default. If the conditional clause is removed, then */
    /* information on the size of workspace needs to be changed. */
    /* If COMPZ = 'N', use SSTERF to compute the eigenvalues. */
    if (icompz == 0)
    {
        ssterf_(n, &d__[1], &e[1], info);
        goto L50;
    }
    /* If N is smaller than the minimum divide size (SMLSIZ+1), then */
    /* solve the problem with another solver. */
    if (*n <= smlsiz)
    {
        ssteqr_(compz, n, &d__[1], &e[1], &z__[z_offset], ldz, &work[1], info);
    }
    else
    {
        /* If COMPZ = 'V', the Z matrix must be stored elsewhere for later */
        /* use. */
        if (icompz == 1)
        {
            storez = *n * *n + 1;
        }
        else
        {
            storez = 1;
        }
        if (icompz == 2)
        {
            slaset_("Full", n, n, &c_b17, &c_b18, &z__[z_offset], ldz);
        }
        /* Scale. */
        orgnrm = slanst_("M", n, &d__[1], &e[1]);
        if (orgnrm == 0.f)
        {
            goto L50;
        }
        eps = slamch_("Epsilon");
        start = 1;
        /* while ( START <= N ) */
L10:
        if (start <= *n)
        {
            /* Let FINISH be the position of the next subdiagonal entry */
            /* such that E( FINISH ) <= TINY or FINISH = N if no such */
            /* subdiagonal exists. The matrix identified by the elements */
            /* between START and FINISH constitutes an independent */
            /* sub-problem. */
            finish = start;
L20:
            if (finish < *n)
            {
                tiny = eps * sqrt((r__1 = d__[finish], f2c_abs(r__1))) * sqrt(( r__2 = d__[finish + 1], f2c_abs(r__2)));
                if ((r__1 = e[finish], f2c_abs(r__1)) > tiny)
                {
                    ++finish;
                    goto L20;
                }
            }
            /* (Sub) Problem determined. Compute its size and solve it. */
            m = finish - start + 1;
            if (m == 1)
            {
                start = finish + 1;
                goto L10;
            }
            if (m > smlsiz)
            {
                /* Scale. */
                orgnrm = slanst_("M", &m, &d__[start], &e[start]);
                slascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &m, &c__1, &d__[ start], &m, info);
                i__1 = m - 1;
                i__2 = m - 1;
                slascl_("G", &c__0, &c__0, &orgnrm, &c_b18, &i__1, &c__1, &e[ start], &i__2, info);
                if (icompz == 1)
                {
                    strtrw = 1;
                }
                else
                {
                    strtrw = start;
                }
                slaed0_(&icompz, n, &m, &d__[start], &e[start], &z__[strtrw + start * z_dim1], ldz, &work[1], n, &work[storez], & iwork[1], info);
                if (*info != 0)
                {
                    *info = (*info / (m + 1) + start - 1) * (*n + 1) + *info % (m + 1) + start - 1;
                    goto L50;
                }
                /* Scale back. */
                slascl_("G", &c__0, &c__0, &c_b18, &orgnrm, &m, &c__1, &d__[ start], &m, info);
            }
            else
            {
                if (icompz == 1)
                {
                    /* Since QR won't update a Z matrix which is larger than */
                    /* the length of D, we must solve the sub-problem in a */
                    /* workspace and then multiply back into Z. */
                    ssteqr_("I", &m, &d__[start], &e[start], &work[1], &m, & work[m * m + 1], info);
                    slacpy_("A", n, &m, &z__[start * z_dim1 + 1], ldz, &work[ storez], n);
                    sgemm_("N", "N", n, &m, &m, &c_b18, &work[storez], n, & work[1], &m, &c_b17, &z__[start * z_dim1 + 1], ldz);
                }
                else if (icompz == 2)
                {
                    ssteqr_("I", &m, &d__[start], &e[start], &z__[start + start * z_dim1], ldz, &work[1], info);
                }
                else
                {
                    ssterf_(&m, &d__[start], &e[start], info);
                }
                if (*info != 0)
                {
                    *info = start * (*n + 1) + finish;
                    goto L50;
                }
            }
            start = finish + 1;
            goto L10;
        }
        /* endwhile */
        /* If the problem split any number of times, then the eigenvalues */
        /* will not be properly ordered. Here we permute the eigenvalues */
        /* (and the associated eigenvectors) into ascending order. */
        if (m != *n)
        {
            if (icompz == 0)
            {
                /* Use Quick Sort */
                slasrt_("I", n, &d__[1], info);
            }
            else
            {
                /* Use Selection Sort to minimize swaps of eigenvectors */
                i__1 = *n;
                for (ii = 2;
                        ii <= i__1;
                        ++ii)
                {
                    i__ = ii - 1;
                    k = i__;
                    p = d__[i__];
                    i__2 = *n;
                    for (j = ii;
                            j <= i__2;
                            ++j)
                    {
                        if (d__[j] < p)
                        {
                            k = j;
                            p = d__[j];
                        }
                        /* L30: */
                    }
                    if (k != i__)
                    {
                        d__[k] = d__[i__];
                        d__[i__] = p;
                        sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1], &c__1);
                    }
                    /* L40: */
                }
            }
        }
    }
L50:
    work[1] = (real) lwmin;
    iwork[1] = liwmin;
    return 0;
    /* End of SSTEDC */
}
/* sstedc_ */
