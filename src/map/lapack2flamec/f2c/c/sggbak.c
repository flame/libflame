/* sggbak.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SGGBAK */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGGBAK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggbak. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggbak. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggbak. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, */
/* LDV, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOB, SIDE */
/* INTEGER IHI, ILO, INFO, LDV, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL LSCALE( * ), RSCALE( * ), V( LDV, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGBAK forms the right or left eigenvectors of a real generalized */
/* > eigenvalue problem A*x = lambda*B*x, by backward transformation on */
/* > the computed eigenvectors of the balanced pair of matrices output by */
/* > SGGBAL. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOB */
/* > \verbatim */
/* > JOB is CHARACTER*1 */
/* > Specifies the type of backward transformation required: */
/* > = 'N': do nothing, return immediately;
*/
/* > = 'P': do backward transformation for permutation only;
*/
/* > = 'S': do backward transformation for scaling only;
*/
/* > = 'B': do backward transformations for both permutation and */
/* > scaling. */
/* > JOB must be the same as the argument JOB supplied to SGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'R': V contains right eigenvectors;
*/
/* > = 'L': V contains left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of rows of the matrix V. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > The integers ILO and IHI determined by SGGBAL. */
/* > 1 <= ILO <= IHI <= N, if N > 0;
ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in] LSCALE */
/* > \verbatim */
/* > LSCALE is REAL array, dimension (N) */
/* > Details of the permutations and/or scaling factors applied */
/* > to the left side of A and B, as returned by SGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] RSCALE */
/* > \verbatim */
/* > RSCALE is REAL array, dimension (N) */
/* > Details of the permutations and/or scaling factors applied */
/* > to the right side of A and B, as returned by SGGBAL. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of columns of the matrix V. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* > V is REAL array, dimension (LDV,M) */
/* > On entry, the matrix of right or left eigenvectors to be */
/* > transformed, as returned by STGEVC. */
/* > On exit, V is overwritten by the transformed eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the matrix V. LDV >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup realGBcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > See R.C. Ward, Balancing the generalized eigenvalue problem, */
/* > SIAM J. Sci. Stat. Comp. 2 (1981), 141-152. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int sggbak_(char *job, char *side, integer *n, integer *ilo, integer *ihi, real *lscale, real *rscale, integer *m, real *v, integer *ldv, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sggbak inputs: job %c, side %c, n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", m %" FLA_IS ", ldv %" FLA_IS "",*job, *side, *n, *ilo, *ihi, *m, *ldv);
    /* System generated locals */
    integer v_dim1, v_offset, i__1;
    /* Local variables */
    integer i__, k;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    logical leftv;
    extern /* Subroutine */
    int sswap_(integer *, real *, integer *, real *, integer *), xerbla_(char *, integer *);
    logical rightv;
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    --lscale;
    --rscale;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    /* Function Body */
    rightv = lsame_(side, "R");
    leftv = lsame_(side, "L");
    *info = 0;
    if (! lsame_(job, "N") && ! lsame_(job, "P") && ! lsame_(job, "S") && ! lsame_(job, "B"))
    {
        *info = -1;
    }
    else if (! rightv && ! leftv)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*ilo < 1)
    {
        *info = -4;
    }
    else if (*n == 0 && *ihi == 0 && *ilo != 1)
    {
        *info = -4;
    }
    else if (*n > 0 && (*ihi < *ilo || *ihi > fla_max(1,*n)))
    {
        *info = -5;
    }
    else if (*n == 0 && *ilo == 1 && *ihi != 0)
    {
        *info = -5;
    }
    else if (*m < 0)
    {
        *info = -8;
    }
    else if (*ldv < fla_max(1,*n))
    {
        *info = -10;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGGBAK", &i__1);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (*m == 0)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (lsame_(job, "N"))
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (*ilo == *ihi)
    {
        goto L30;
    }
    /* Backward balance */
    if (lsame_(job, "S") || lsame_(job, "B"))
    {
        /* Backward transformation on right eigenvectors */
        if (rightv)
        {
            i__1 = *ihi;
            for (i__ = *ilo;
                    i__ <= i__1;
                    ++i__)
            {
                sscal_(m, &rscale[i__], &v[i__ + v_dim1], ldv);
                /* L10: */
            }
        }
        /* Backward transformation on left eigenvectors */
        if (leftv)
        {
            i__1 = *ihi;
            for (i__ = *ilo;
                    i__ <= i__1;
                    ++i__)
            {
                sscal_(m, &lscale[i__], &v[i__ + v_dim1], ldv);
                /* L20: */
            }
        }
    }
    /* Backward permutation */
L30:
    if (lsame_(job, "P") || lsame_(job, "B"))
    {
        /* Backward permutation on right eigenvectors */
        if (rightv)
        {
            if (*ilo == 1)
            {
                goto L50;
            }
            for (i__ = *ilo - 1;
                    i__ >= 1;
                    --i__)
            {
                k = (integer) rscale[i__];
                if (k == i__)
                {
                    goto L40;
                }
                sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
L40:
                ;
            }
L50:
            if (*ihi == *n)
            {
                goto L70;
            }
            i__1 = *n;
            for (i__ = *ihi + 1;
                    i__ <= i__1;
                    ++i__)
            {
                k = (integer) rscale[i__];
                if (k == i__)
                {
                    goto L60;
                }
                sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
L60:
                ;
            }
        }
        /* Backward permutation on left eigenvectors */
L70:
        if (leftv)
        {
            if (*ilo == 1)
            {
                goto L90;
            }
            for (i__ = *ilo - 1;
                    i__ >= 1;
                    --i__)
            {
                k = (integer) lscale[i__];
                if (k == i__)
                {
                    goto L80;
                }
                sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
L80:
                ;
            }
L90:
            if (*ihi == *n)
            {
                goto L110;
            }
            i__1 = *n;
            for (i__ = *ihi + 1;
                    i__ <= i__1;
                    ++i__)
            {
                k = (integer) lscale[i__];
                if (k == i__)
                {
                    goto L100;
                }
                sswap_(m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
L100:
                ;
            }
        }
    }
L110:
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of SGGBAK */
}
/* sggbak_ */
