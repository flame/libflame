/* dlaqz1.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLAQZ1 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQZ1 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqz1. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqz1. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqz1. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQZ1( A, LDA, B, LDB, SR1, SR2, SI, BETA1, BETA2, */
/* $ V ) */
/* IMPLICIT NONE */
/* Arguments */
/* INTEGER, INTENT( IN ) :: LDA, LDB */
/* DOUBLE PRECISION, INTENT( IN ) :: A( LDA, * ), B( LDB, * ), SR1, */
/* $ SR2, SI, BETA1, BETA2 */
/* DOUBLE PRECISION, INTENT( OUT ) :: V( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given a 3-by-3 matrix pencil (A,B), DLAQZ1 sets v to a */
/* > scalar multiple of the first column of the product */
/* > */
/* > (*) K = (A - (beta2*sr2 - i*si)*B)*B^(-1)*(beta1*A - (sr2 + i*si2)*B)*B^(-1). */
/* > */
/* > It is assumed that either */
/* > */
/* > 1) sr1 = sr2 */
/* > or */
/* > 2) si = 0. */
/* > */
/* > This is useful for starting double implicit shift bulges */
/* > in the QZ algorithm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > The 3-by-3 matrix A in (*). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of A as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* > \param[in] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,N) */
/* > The 3-by-3 matrix B in (*). */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* > */
/* > \param[in] SR1 */
/* > \verbatim */
/* > SR1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] SR2 */
/* > \verbatim */
/* > SR2 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] SI */
/* > \verbatim */
/* > SI is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] BETA1 */
/* > \verbatim */
/* > BETA1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] BETA2 */
/* > \verbatim */
/* > BETA2 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is DOUBLE PRECISION array, dimension (N) */
/* > A scalar multiple of the first column of the */
/* > matrix K in (*). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Thijs Steel, KU Leuven */
/* > \date May 2020 */
/* > \ingroup doubleGEcomputational */
/* > */
/* ===================================================================== */
/* Subroutine */
int dlaqz1_(doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *sr1, doublereal *sr2, doublereal *si, doublereal *beta1, doublereal *beta2, doublereal *v)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlaqz1 inputs: lda %" FLA_IS ", ldb %" FLA_IS "",*lda, *ldb);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    doublereal w[2], scale1, scale2;
    extern doublereal dlamch_(char *);
    extern logical disnan_(doublereal *);
    doublereal safmin, safmax;
    /* Arguments */
    /* Parameters */
    /* Local scalars */
    /* External Functions */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --v;
    /* Function Body */
    safmin = dlamch_("SAFE MINIMUM");
    safmax = 1. / safmin;
    /* Calculate first shifted vector */
    w[0] = *beta1 * a[a_dim1 + 1] - *sr1 * b[b_dim1 + 1];
    w[1] = *beta1 * a[a_dim1 + 2] - *sr1 * b[b_dim1 + 2];
    scale1 = sqrt((f2c_abs(w[0]))) * sqrt((f2c_abs(w[1])));
    if (scale1 >= safmin && scale1 <= safmax)
    {
        w[0] /= scale1;
        w[1] /= scale1;
    }
    /* Solve linear system */
    w[1] /= b[(b_dim1 << 1) + 2];
    w[0] = (w[0] - b[(b_dim1 << 1) + 1] * w[1]) / b[b_dim1 + 1];
    scale2 = sqrt((f2c_abs(w[0]))) * sqrt((f2c_abs(w[1])));
    if (scale2 >= safmin && scale2 <= safmax)
    {
        w[0] /= scale2;
        w[1] /= scale2;
    }
    /* Apply second shift */
    v[1] = *beta2 * (a[a_dim1 + 1] * w[0] + a[(a_dim1 << 1) + 1] * w[1]) - * sr2 * (b[b_dim1 + 1] * w[0] + b[(b_dim1 << 1) + 1] * w[1]);
    v[2] = *beta2 * (a[a_dim1 + 2] * w[0] + a[(a_dim1 << 1) + 2] * w[1]) - * sr2 * (b[b_dim1 + 2] * w[0] + b[(b_dim1 << 1) + 2] * w[1]);
    v[3] = *beta2 * (a[a_dim1 + 3] * w[0] + a[(a_dim1 << 1) + 3] * w[1]) - * sr2 * (b[b_dim1 + 3] * w[0] + b[(b_dim1 << 1) + 3] * w[1]);
    /* Account for imaginary part */
    v[1] += *si * *si * b[b_dim1 + 1] / scale1 / scale2;
    /* Check for overflow */
    if (f2c_abs(v[1]) > safmax || f2c_abs(v[2]) > safmax || f2c_abs(v[3]) > safmax || disnan_(&v[1]) || disnan_(&v[2]) || disnan_(&v[3]))
    {
        v[1] = 0.;
        v[2] = 0.;
        v[3] = 0.;
    }
    /* End of DLAQZ1 */
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
}
/* dlaqz1_ */

