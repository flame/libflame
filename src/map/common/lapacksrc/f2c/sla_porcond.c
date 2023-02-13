/* ../netlib/sla_porcond.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SLA_PORCOND estimates the Skeel condition number for a symmetric positive-definite matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLA_PORCOND + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_por cond.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_por cond.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_por cond.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION SLA_PORCOND( UPLO, N, A, LDA, AF, LDAF, CMODE, C, */
/* INFO, WORK, IWORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER N, LDA, LDAF, INFO, CMODE */
/* REAL A( LDA, * ), AF( LDAF, * ), WORK( * ), */
/* $ C( * ) */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLA_PORCOND Estimates the Skeel condition number of op(A) * op2(C) */
/* > where op2 is determined by CMODE as follows */
/* > CMODE = 1 op2(C) = C */
/* > CMODE = 0 op2(C) = I */
/* > CMODE = -1 op2(C) = inv(C) */
/* > The Skeel condition number cond(A) = norminf( |inv(A)||A| ) */
/* > is computed by computing scaling factors R such that */
/* > diag(R)*A*op2(C) is row equilibrated and computing the standard */
/* > infinity-norm condition number. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* > AF is REAL array, dimension (LDAF,N) */
/* > The triangular factor U or L from the Cholesky factorization */
/* > A = U**T*U or A = L*L**T, as computed by SPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] CMODE */
/* > \verbatim */
/* > CMODE is INTEGER */
/* > Determines op2(C) in the formula op(A) * op2(C) as follows: */
/* > CMODE = 1 op2(C) = C */
/* > CMODE = 0 op2(C) = I */
/* > CMODE = -1 op2(C) = inv(C) */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is REAL array, dimension (N) */
/* > The vector C in the formula op(A) * op2(C). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: Successful exit. */
/* > i > 0: The ith argument is invalid. */
/* > \endverbatim */
/* > */
/* > \param[in] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (3*N). */
/* > Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N). */
/* > Workspace. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realPOcomputational */
/* ===================================================================== */
real sla_porcond_(char *uplo, integer *n, real *a, integer *lda, real *af, integer *ldaf, integer *cmode, real *c__, integer *info, real *work, integer *iwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2;
    real ret_val, r__1;
    /* Local variables */
    integer i__, j;
    logical up;
    real tmp;
    integer kase;
    extern logical lsame_(char *, char *);
    integer isave[3];
    extern /* Subroutine */
    int slacn2_(integer *, real *, real *, integer *, real *, integer *, integer *), xerbla_(char *, integer *);
    real ainvnm;
    extern /* Subroutine */
    int spotrs_(char *, integer *, integer *, real *, integer *, real *, integer *, integer *);
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    --c__;
    --work;
    --iwork;
    /* Function Body */
    ret_val = 0.f;
    *info = 0;
    if (*n < 0)
    {
        *info = -2;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLA_PORCOND", &i__1);
        return ret_val;
    }
    if (*n == 0)
    {
        ret_val = 1.f;
        return ret_val;
    }
    up = FALSE_;
    if (lsame_(uplo, "U"))
    {
        up = TRUE_;
    }
    /* Compute the equilibration matrix R such that */
    /* inv(R)*A*C has unit 1-norm. */
    if (up)
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            tmp = 0.f;
            if (*cmode == 1)
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[j + i__ * a_dim1] * c__[j], f2c_abs(r__1));
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[i__ + j * a_dim1] * c__[j], f2c_abs(r__1));
                }
            }
            else if (*cmode == 0)
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[j + i__ * a_dim1], f2c_abs(r__1));
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                }
            }
            else
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[j + i__ * a_dim1] / c__[j], f2c_abs(r__1));
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[i__ + j * a_dim1] / c__[j], f2c_abs(r__1));
                }
            }
            work[(*n << 1) + i__] = tmp;
        }
    }
    else
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            tmp = 0.f;
            if (*cmode == 1)
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[i__ + j * a_dim1] * c__[j], f2c_abs(r__1));
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[j + i__ * a_dim1] * c__[j], f2c_abs(r__1));
                }
            }
            else if (*cmode == 0)
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[j + i__ * a_dim1], f2c_abs(r__1));
                }
            }
            else
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[i__ + j * a_dim1] / c__[j], f2c_abs(r__1));
                }
                i__2 = *n;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    tmp += (r__1 = a[j + i__ * a_dim1] / c__[j], f2c_abs(r__1));
                }
            }
            work[(*n << 1) + i__] = tmp;
        }
    }
    /* Estimate the norm of inv(op(A)). */
    ainvnm = 0.f;
    kase = 0;
L10:
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if (kase != 0)
    {
        if (kase == 2)
        {
            /* Multiply by R. */
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                work[i__] *= work[(*n << 1) + i__];
            }
            if (up)
            {
                spotrs_("Upper", n, &c__1, &af[af_offset], ldaf, &work[1], n, info);
            }
            else
            {
                spotrs_("Lower", n, &c__1, &af[af_offset], ldaf, &work[1], n, info);
            }
            /* Multiply by inv(C). */
            if (*cmode == 1)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] /= c__[i__];
                }
            }
            else if (*cmode == -1)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] *= c__[i__];
                }
            }
        }
        else
        {
            /* Multiply by inv(C**T). */
            if (*cmode == 1)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] /= c__[i__];
                }
            }
            else if (*cmode == -1)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] *= c__[i__];
                }
            }
            if (up)
            {
                spotrs_("Upper", n, &c__1, &af[af_offset], ldaf, &work[1], n, info);
            }
            else
            {
                spotrs_("Lower", n, &c__1, &af[af_offset], ldaf, &work[1], n, info);
            }
            /* Multiply by R. */
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                work[i__] *= work[(*n << 1) + i__];
            }
        }
        goto L10;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if (ainvnm != 0.f)
    {
        ret_val = 1.f / ainvnm;
    }
    return ret_val;
}
/* sla_porcond__ */
