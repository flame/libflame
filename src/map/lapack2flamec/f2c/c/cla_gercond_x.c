/* ../netlib/cla_gercond_x.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLA_GERCOND_X computes the infinity norm condition number of op(A)*diag(x) for general matrices . */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLA_GERCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_ger cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_ger cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_ger cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLA_GERCOND_X( TRANS, N, A, LDA, AF, LDAF, IPIV, X, */
/* INFO, WORK, RWORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER N, LDA, LDAF, INFO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * ) */
/* REAL RWORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > CLA_GERCOND_X computes the infinity norm condition number of */
/* > op(A) * diag(X) where X is a COMPLEX vector. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations: */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate Transpose = Transpose) */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* > AF is COMPLEX array, dimension (LDAF,N) */
/* > The factors L and U from the factorization */
/* > A = P*L*U as computed by CGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices from the factorization A = P*L*U */
/* > as computed by CGETRF;
row i of the matrix was interchanged */
/* > with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (N) */
/* > The vector X in the formula op(A) * diag(X). */
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
/* > WORK is COMPLEX array, dimension (2*N). */
/* > Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (N). */
/* > Workspace. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexGEcomputational */
/* ===================================================================== */
real cla_gercond_x_(char *trans, integer *n, complex *a, integer *lda, complex *af, integer *ldaf, integer *ipiv, complex *x, integer *info, complex *work, real *rwork)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"cla_gercond_x inputs: trans %c, n %lld, lda %lld, ldaf %lld",*trans, *n, *lda, *ldaf);
#else
    snprintf(buffer, 256,"cla_gercond_x inputs: trans %c, n %d, lda %d, ldaf %d",*trans, *n, *lda, *ldaf);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3, i__4;
    real ret_val, r__1, r__2;
    complex q__1, q__2;
    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    integer i__, j;
    real tmp;
    integer kase;
    extern logical lsame_(char *, char *);
    integer isave[3];
    real anorm;
    extern /* Subroutine */
    int clacn2_(integer *, complex *, complex *, real *, integer *, integer *), xerbla_(char *, integer *), cgetrs_(char *, integer *, integer *, complex *, integer *, integer *, complex *, integer *, integer *);
    real ainvnm;
    logical notrans;
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function Definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    --ipiv;
    --x;
    --work;
    --rwork;
    /* Function Body */
    ret_val = 0.f;
    *info = 0;
    notrans = lsame_(trans, "N");
    if (! notrans && ! lsame_(trans, "T") && ! lsame_( trans, "C"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*n))
    {
        *info = -4;
    }
    else if (*ldaf < fla_max(1,*n))
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CLA_GERCOND_X", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return ret_val;
    }
    /* Compute norm of op(A)*op2(C). */
    anorm = 0.f;
    if (notrans)
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            tmp = 0.f;
            i__2 = *n;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ + j * a_dim1;
                i__4 = j;
                q__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i;
                q__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4] .r; // , expr subst
                q__1.r = q__2.r;
                q__1.i = q__2.i; // , expr subst
                tmp += (r__1 = q__1.r, f2c_abs(r__1)) + (r__2 = r_imag(&q__1), f2c_abs(r__2));
            }
            rwork[i__] = tmp;
            anorm = fla_max(anorm,tmp);
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
            i__2 = *n;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = j + i__ * a_dim1;
                i__4 = j;
                q__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i;
                q__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4] .r; // , expr subst
                q__1.r = q__2.r;
                q__1.i = q__2.i; // , expr subst
                tmp += (r__1 = q__1.r, f2c_abs(r__1)) + (r__2 = r_imag(&q__1), f2c_abs(r__2));
            }
            rwork[i__] = tmp;
            anorm = fla_max(anorm,tmp);
        }
    }
    /* Quick return if possible. */
    if (*n == 0)
    {
        ret_val = 1.f;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return ret_val;
    }
    else if (anorm == 0.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return ret_val;
    }
    /* Estimate the norm of inv(op(A)). */
    ainvnm = 0.f;
    kase = 0;
L10:
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
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
                i__2 = i__;
                i__3 = i__;
                i__4 = i__;
                q__1.r = rwork[i__4] * work[i__3].r;
                q__1.i = rwork[i__4] * work[i__3].i; // , expr subst
                work[i__2].r = q__1.r;
                work[i__2].i = q__1.i; // , expr subst
            }
            if (notrans)
            {
                cgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[ 1], &work[1], n, info);
            }
            else
            {
                cgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], n, info);
            }
            /* Multiply by inv(X). */
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = i__;
                c_div(&q__1, &work[i__], &x[i__]);
                work[i__2].r = q__1.r;
                work[i__2].i = q__1.i; // , expr subst
            }
        }
        else
        {
            /* Multiply by inv(X**H). */
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = i__;
                c_div(&q__1, &work[i__], &x[i__]);
                work[i__2].r = q__1.r;
                work[i__2].i = q__1.i; // , expr subst
            }
            if (notrans)
            {
                cgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], n, info);
            }
            else
            {
                cgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[ 1], &work[1], n, info);
            }
            /* Multiply by R. */
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = i__;
                i__3 = i__;
                i__4 = i__;
                q__1.r = rwork[i__4] * work[i__3].r;
                q__1.i = rwork[i__4] * work[i__3].i; // , expr subst
                work[i__2].r = q__1.r;
                work[i__2].i = q__1.i; // , expr subst
            }
        }
        goto L10;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if (ainvnm != 0.f)
    {
        ret_val = 1.f / ainvnm;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return ret_val;
}
/* cla_gercond_x__ */
