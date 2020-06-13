/* ../netlib/zla_gercond_x.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZLA_GERCOND_X computes the infinity norm condition number of op(A)*diag(x) for general matrices . */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLA_GERCOND_X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_ger cond_x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_ger cond_x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_ger cond_x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* DOUBLE PRECISION FUNCTION ZLA_GERCOND_X( TRANS, N, A, LDA, AF, */
/* LDAF, IPIV, X, INFO, */
/* WORK, RWORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER N, LDA, LDAF, INFO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * ) */
/* DOUBLE PRECISION RWORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLA_GERCOND_X computes the infinity norm condition number of */
/* > op(A) * diag(X) where X is a COMPLEX*16 vector. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
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
/* > AF is COMPLEX*16 array, dimension (LDAF,N) */
/* > The factors L and U from the factorization */
/* > A = P*L*U as computed by ZGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices from the factorization A = P*L*U */
/* > as computed by ZGETRF;
row i of the matrix was interchanged */
/* > with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (N) */
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
/* > WORK is COMPLEX*16 array, dimension (2*N). */
/* > Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension (N). */
/* > Workspace. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16GEcomputational */
/* ===================================================================== */
doublereal zla_gercond_x_(char *trans, integer *n, doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex * x, integer *info, doublecomplex *work, doublereal *rwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1, z__2;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    /* Local variables */
    integer i__, j;
    doublereal tmp;
    integer kase;
    extern logical lsame_(char *, char *);
    integer isave[3];
    doublereal anorm;
    extern /* Subroutine */
    int zlacn2_(integer *, doublecomplex *, doublecomplex *, doublereal *, integer *, integer *), xerbla_( char *, integer *);
    doublereal ainvnm;
    extern /* Subroutine */
    int zgetrs_(char *, integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *, integer *);
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
    ret_val = 0.;
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
    else if (*lda < max(1,*n))
    {
        *info = -4;
    }
    else if (*ldaf < max(1,*n))
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZLA_GERCOND_X", &i__1);
        return ret_val;
    }
    /* Compute norm of op(A)*op2(C). */
    anorm = 0.;
    if (notrans)
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            tmp = 0.;
            i__2 = *n;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = i__ + j * a_dim1;
                i__4 = j;
                z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i;
                z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4] .r; // , expr subst
                z__1.r = z__2.r;
                z__1.i = z__2.i; // , expr subst
                tmp += (d__1 = z__1.r, f2c_dabs(d__1)) + (d__2 = d_imag(&z__1), f2c_dabs(d__2));
            }
            rwork[i__] = tmp;
            anorm = max(anorm,tmp);
        }
    }
    else
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            tmp = 0.;
            i__2 = *n;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                i__3 = j + i__ * a_dim1;
                i__4 = j;
                z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i;
                z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4] .r; // , expr subst
                z__1.r = z__2.r;
                z__1.i = z__2.i; // , expr subst
                tmp += (d__1 = z__1.r, f2c_dabs(d__1)) + (d__2 = d_imag(&z__1), f2c_dabs(d__2));
            }
            rwork[i__] = tmp;
            anorm = max(anorm,tmp);
        }
    }
    /* Quick return if possible. */
    if (*n == 0)
    {
        ret_val = 1.;
        return ret_val;
    }
    else if (anorm == 0.)
    {
        return ret_val;
    }
    /* Estimate the norm of inv(op(A)). */
    ainvnm = 0.;
    kase = 0;
L10:
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
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
                z__1.r = rwork[i__4] * work[i__3].r;
                z__1.i = rwork[i__4] * work[i__3].i; // , expr subst
                work[i__2].r = z__1.r;
                work[i__2].i = z__1.i; // , expr subst
            }
            if (notrans)
            {
                zgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[ 1], &work[1], n, info);
            }
            else
            {
                zgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], n, info);
            }
            /* Multiply by inv(X). */
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = i__;
                z_div(&z__1, &work[i__], &x[i__]);
                work[i__2].r = z__1.r;
                work[i__2].i = z__1.i; // , expr subst
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
                z_div(&z__1, &work[i__], &x[i__]);
                work[i__2].r = z__1.r;
                work[i__2].i = z__1.i; // , expr subst
            }
            if (notrans)
            {
                zgetrs_("Conjugate transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], n, info);
            }
            else
            {
                zgetrs_("No transpose", n, &c__1, &af[af_offset], ldaf, &ipiv[ 1], &work[1], n, info);
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
                z__1.r = rwork[i__4] * work[i__3].r;
                z__1.i = rwork[i__4] * work[i__3].i; // , expr subst
                work[i__2].r = z__1.r;
                work[i__2].i = z__1.i; // , expr subst
            }
        }
        goto L10;
    }
    /* Compute the estimate of the reciprocal condition number. */
    if (ainvnm != 0.)
    {
        ret_val = 1. / ainvnm;
    }
    return ret_val;
}
/* zla_gercond_x__ */
