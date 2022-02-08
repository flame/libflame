/* dlaqz2.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c__1 = 1;
/* > \brief \b DLAQZ2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQZ2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqz2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqz2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqz2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQZ2( ILQ, ILZ, K, ISTARTM, ISTOPM, IHI, A, LDA, B, */
/* $ LDB, NQ, QSTART, Q, LDQ, NZ, ZSTART, Z, LDZ ) */
/* IMPLICIT NONE */
/* Arguments */
/* LOGICAL, INTENT( IN ) :: ILQ, ILZ */
/* INTEGER, INTENT( IN ) :: K, LDA, LDB, LDQ, LDZ, ISTARTM, ISTOPM, */
/* $ NQ, NZ, QSTART, ZSTART, IHI */
/* DOUBLE PRECISION :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, */
/* $ * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQZ2 chases a 2x2 shift bulge in a matrix pencil down a single position */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > */
/* > \param[in] ILQ */
/* > \verbatim */
/* > ILQ is LOGICAL */
/* > Determines whether or not to update the matrix Q */
/* > \endverbatim */
/* > */
/* > \param[in] ILZ */
/* > \verbatim */
/* > ILZ is LOGICAL */
/* > Determines whether or not to update the matrix Z */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > Index indicating the position of the bulge. */
/* > On entry, the bulge is located in */
/* > (A(k+1:k+2,k:k+1),B(k+1:k+2,k:k+1)). */
/* > On exit, the bulge is located in */
/* > (A(k+2:k+3,k+1:k+2),B(k+2:k+3,k+1:k+2)). */
/* > \endverbatim */
/* > */
/* > \param[in] ISTARTM */
/* > \verbatim */
/* > ISTARTM is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] ISTOPM */
/* > \verbatim */
/* > ISTOPM is INTEGER */
/* > Updates to (A,B) are restricted to */
/* > (istartm:k+3,k:istopm). It is assumed */
/* > without checking that istartm <= k+1 and */
/* > k+2 <= istopm */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[inout] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of A as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* > \param[inout] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* > */
/* > \param[in] NQ */
/* > \verbatim */
/* > NQ is INTEGER */
/* > The order of the matrix Q */
/* > \endverbatim */
/* > */
/* > \param[in] QSTART */
/* > \verbatim */
/* > QSTART is INTEGER */
/* > Start index of the matrix Q. Rotations are applied */
/* > To columns k+2-qStart:k+4-qStart of Q. */
/* > \endverbatim */
/* > \param[inout] Q */
/* > \verbatim */
/* > Q is DOUBLE PRECISION array, dimension (LDQ,NQ) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of Q as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* > */
/* > \param[in] NZ */
/* > \verbatim */
/* > NZ is INTEGER */
/* > The order of the matrix Z */
/* > \endverbatim */
/* > */
/* > \param[in] ZSTART */
/* > \verbatim */
/* > ZSTART is INTEGER */
/* > Start index of the matrix Z. Rotations are applied */
/* > To columns k+1-qStart:k+3-qStart of Z. */
/* > \endverbatim */
/* > \param[inout] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension (LDZ,NZ) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of Q as declared in */
/* > the calling procedure. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Thijs Steel, KU Leuven */
/* > \date May 2020 */
/* > \ingroup doubleGEcomputational */
/* > */
/* ===================================================================== */
/* Subroutine */
int dlaqz2_(logical *ilq, logical *ilz, integer *k, integer * istartm, integer *istopm, integer *ihi, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *nq, integer *qstart, doublereal *q, integer *ldq, integer *nz, integer *zstart, doublereal *z__, integer *ldz)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"dlaqz2 inputs: k %" FLA_IS ", istartm %" FLA_IS ", istopm %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS ", nq %" FLA_IS ", qstart %" FLA_IS ", ldq %" FLA_IS ", nz %" FLA_IS ", zstart %" FLA_IS ", ldz %" FLA_IS "",*k, *istartm, *istopm, *ihi, *lda, *ldb, *nq, *qstart, *ldq, *nz, *zstart, *ldz);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2;
    /* Local variables */
    doublereal h__[6] /* was [2][3] */
    , c1, c2, s1, s2;
    integer indi, indj;
    doublereal temp;
    extern /* Subroutine */
    int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *), dlartg_( doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    /* Arguments */
    /* Parameters */
    /* Local variables */
    /* External functions */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    /* Function Body */
    if (*k + 2 == *ihi)
    {
        /* Shift is located on the edge of the matrix, remove it */
        /* H = B( IHI-1:IHI, IHI-2:IHI ) */
        i__1 = *ihi;
        for (indi = *ihi - 1;
                indi <= i__1;
                ++indi)
        {
            i__2 = *ihi;
            for (indj = *ihi - 2;
                    indj <= i__2;
                    ++indj)
            {
                h__[indi + (indj << 1) - 3] = b[indi + indj * b_dim1];
            }
        }
        /* Make H upper triangular */
        dlartg_(h__, &h__[1], &c1, &s1, &temp);
        h__[1] = 0.;
        h__[0] = temp;
        drot_(&c__2, &h__[2], &c__2, &h__[3], &c__2, &c1, &s1);
        dlartg_(&h__[5], &h__[3], &c1, &s1, &temp);
        drot_(&c__1, &h__[4], &c__1, &h__[2], &c__1, &c1, &s1);
        dlartg_(&h__[2], h__, &c2, &s2, &temp);
        i__1 = *ihi - *istartm + 1;
        drot_(&i__1, &b[*istartm + *ihi * b_dim1], &c__1, &b[*istartm + (*ihi - 1) * b_dim1], &c__1, &c1, &s1);
        i__1 = *ihi - *istartm + 1;
        drot_(&i__1, &b[*istartm + (*ihi - 1) * b_dim1], &c__1, &b[*istartm + (*ihi - 2) * b_dim1], &c__1, &c2, &s2);
        b[*ihi - 1 + (*ihi - 2) * b_dim1] = 0.;
        b[*ihi + (*ihi - 2) * b_dim1] = 0.;
        i__1 = *ihi - *istartm + 1;
        drot_(&i__1, &a[*istartm + *ihi * a_dim1], &c__1, &a[*istartm + (*ihi - 1) * a_dim1], &c__1, &c1, &s1);
        i__1 = *ihi - *istartm + 1;
        drot_(&i__1, &a[*istartm + (*ihi - 1) * a_dim1], &c__1, &a[*istartm + (*ihi - 2) * a_dim1], &c__1, &c2, &s2);
        if (*ilz)
        {
            drot_(nz, &z__[(*ihi - *zstart + 1) * z_dim1 + 1], &c__1, &z__[(* ihi - 1 - *zstart + 1) * z_dim1 + 1], &c__1, &c1, &s1);
            drot_(nz, &z__[(*ihi - 1 - *zstart + 1) * z_dim1 + 1], &c__1, & z__[(*ihi - 2 - *zstart + 1) * z_dim1 + 1], &c__1, &c2, & s2);
        }
        dlartg_(&a[*ihi - 1 + (*ihi - 2) * a_dim1], &a[*ihi + (*ihi - 2) * a_dim1], &c1, &s1, &temp);
        a[*ihi - 1 + (*ihi - 2) * a_dim1] = temp;
        a[*ihi + (*ihi - 2) * a_dim1] = 0.;
        i__1 = *istopm - *ihi + 2;
        drot_(&i__1, &a[*ihi - 1 + (*ihi - 1) * a_dim1], lda, &a[*ihi + (*ihi - 1) * a_dim1], lda, &c1, &s1);
        i__1 = *istopm - *ihi + 2;
        drot_(&i__1, &b[*ihi - 1 + (*ihi - 1) * b_dim1], ldb, &b[*ihi + (*ihi - 1) * b_dim1], ldb, &c1, &s1);
        if (*ilq)
        {
            drot_(nq, &q[(*ihi - 1 - *qstart + 1) * q_dim1 + 1], &c__1, &q[(* ihi - *qstart + 1) * q_dim1 + 1], &c__1, &c1, &s1);
        }
        dlartg_(&b[*ihi + *ihi * b_dim1], &b[*ihi + (*ihi - 1) * b_dim1], &c1, &s1, &temp);
        b[*ihi + *ihi * b_dim1] = temp;
        b[*ihi + (*ihi - 1) * b_dim1] = 0.;
        i__1 = *ihi - *istartm;
        drot_(&i__1, &b[*istartm + *ihi * b_dim1], &c__1, &b[*istartm + (*ihi - 1) * b_dim1], &c__1, &c1, &s1);
        i__1 = *ihi - *istartm + 1;
        drot_(&i__1, &a[*istartm + *ihi * a_dim1], &c__1, &a[*istartm + (*ihi - 1) * a_dim1], &c__1, &c1, &s1);
        if (*ilz)
        {
            drot_(nz, &z__[(*ihi - *zstart + 1) * z_dim1 + 1], &c__1, &z__[(* ihi - 1 - *zstart + 1) * z_dim1 + 1], &c__1, &c1, &s1);
        }
    }
    else
    {
        /* Normal operation, move bulge down */
        /* H = B( K+1:K+2, K:K+2 ) */
        i__1 = *k + 2;
        for (indi = *k + 1;
                indi <= i__1;
                ++indi)
        {
            i__2 = *k + 2;
            for (indj = *k;
                    indj <= i__2;
                    ++indj)
            {
                a[indi + indj * a_dim1] = b[indi + indj * b_dim1];
            }
        }
        /* Make H upper triangular */
        dlartg_(h__, &h__[1], &c1, &s1, &temp);
        h__[1] = 0.;
        h__[0] = temp;
        drot_(&c__2, &h__[2], &c__2, &h__[3], &c__2, &c1, &s1);
        /* Calculate Z1 and Z2 */
        dlartg_(&h__[5], &h__[3], &c1, &s1, &temp);
        drot_(&c__1, &h__[4], &c__1, &h__[2], &c__1, &c1, &s1);
        dlartg_(&h__[2], h__, &c2, &s2, &temp);
        /* Apply transformations from the right */
        i__1 = *k + 3 - *istartm + 1;
        drot_(&i__1, &a[*istartm + (*k + 2) * a_dim1], &c__1, &a[*istartm + (* k + 1) * a_dim1], &c__1, &c1, &s1);
        i__1 = *k + 3 - *istartm + 1;
        drot_(&i__1, &a[*istartm + (*k + 1) * a_dim1], &c__1, &a[*istartm + * k * a_dim1], &c__1, &c2, &s2);
        i__1 = *k + 2 - *istartm + 1;
        drot_(&i__1, &b[*istartm + (*k + 2) * b_dim1], &c__1, &b[*istartm + (* k + 1) * b_dim1], &c__1, &c1, &s1);
        i__1 = *k + 2 - *istartm + 1;
        drot_(&i__1, &b[*istartm + (*k + 1) * b_dim1], &c__1, &b[*istartm + * k * b_dim1], &c__1, &c2, &s2);
        if (*ilz)
        {
            drot_(nz, &z__[(*k + 2 - *zstart + 1) * z_dim1 + 1], &c__1, &z__[( *k + 1 - *zstart + 1) * z_dim1 + 1], &c__1, &c1, &s1);
            drot_(nz, &z__[(*k + 1 - *zstart + 1) * z_dim1 + 1], &c__1, &z__[( *k - *zstart + 1) * z_dim1 + 1], &c__1, &c2, &s2);
        }
        b[*k + 1 + *k * b_dim1] = 0.;
        b[*k + 2 + *k * b_dim1] = 0.;
        /* Calculate Q1 and Q2 */
        dlartg_(&a[*k + 2 + *k * a_dim1], &a[*k + 3 + *k * a_dim1], &c1, &s1, &temp);
        a[*k + 2 + *k * a_dim1] = temp;
        a[*k + 3 + *k * a_dim1] = 0.;
        dlartg_(&a[*k + 1 + *k * a_dim1], &a[*k + 2 + *k * a_dim1], &c2, &s2, &temp);
        a[*k + 1 + *k * a_dim1] = temp;
        a[*k + 2 + *k * a_dim1] = 0.;
        /* Apply transformations from the left */
        i__1 = *istopm - *k;
        drot_(&i__1, &a[*k + 2 + (*k + 1) * a_dim1], lda, &a[*k + 3 + (*k + 1) * a_dim1], lda, &c1, &s1);
        i__1 = *istopm - *k;
        drot_(&i__1, &a[*k + 1 + (*k + 1) * a_dim1], lda, &a[*k + 2 + (*k + 1) * a_dim1], lda, &c2, &s2);
        i__1 = *istopm - *k;
        drot_(&i__1, &b[*k + 2 + (*k + 1) * b_dim1], ldb, &b[*k + 3 + (*k + 1) * b_dim1], ldb, &c1, &s1);
        i__1 = *istopm - *k;
        drot_(&i__1, &b[*k + 1 + (*k + 1) * b_dim1], ldb, &b[*k + 2 + (*k + 1) * b_dim1], ldb, &c2, &s2);
        if (*ilq)
        {
            drot_(nq, &q[(*k + 2 - *qstart + 1) * q_dim1 + 1], &c__1, &q[(*k + 3 - *qstart + 1) * q_dim1 + 1], &c__1, &c1, &s1);
            drot_(nq, &q[(*k + 1 - *qstart + 1) * q_dim1 + 1], &c__1, &q[(*k + 2 - *qstart + 1) * q_dim1 + 1], &c__1, &c2, &s2);
        }
    }
    /* End of DLAQZ2 */
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
}
/* dlaqz2_ */

