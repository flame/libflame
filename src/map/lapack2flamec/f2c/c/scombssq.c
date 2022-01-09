/* ../netlib/v3.9.0/scombssq.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* > \brief \b SCOMBSSQ adds two scaled sum of squares quantities */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE SCOMBSSQ( V1, V2 ) */
 /* .. Array Arguments .. */
 /* REAL V1( 2 ), V2( 2 ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > SCOMBSSQ adds two scaled sum of squares quantities, V1 := V1 + V2. */
 /* > That is, */
 /* > */
 /* > V1_scale**2 * V1_sumsq := V1_scale**2 * V1_sumsq */
 /* > + V2_scale**2 * V2_sumsq */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in,out] V1 */
 /* > \verbatim */
 /* > V1 is REAL array, dimension (2). */
 /* > The first scaled sum. */
 /* > V1(1) = V1_scale, V1(2) = V1_sumsq. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] V2 */
 /* > \verbatim */
 /* > V2 is REAL array, dimension (2). */
 /* > The second scaled sum. */
 /* > V2(1) = V2_scale, V2(2) = V2_sumsq. */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2018 */
 /* > \ingroup OTHERauxiliary */
 /* ===================================================================== */
 /* Subroutine */
 int scombssq_(real *v1, real *v2) {
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
 /* System generated locals */
 real r__1;
 /* -- LAPACK auxiliary routine (version 3.7.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* November 2018 */
 /* .. Array Arguments .. */
 /* .. */
 /* ===================================================================== */
 /* .. Parameters .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Parameter adjustments */
 --v2;
 --v1;
 /* Function Body */
 if (v1[1] >= v2[1]) {
 if (v1[1] != 0.f) {
 /* Computing 2nd power */
 r__1 = v2[1] / v1[1];
 v1[2] += r__1 * r__1 * v2[2];
 }
 }
 else {
 /* Computing 2nd power */
 r__1 = v1[1] / v2[1];
 v1[2] = v2[2] + r__1 * r__1 * v1[2];
 v1[1] = v2[1];
 }
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of SCOMBSSQ */
 }
 /* scombssq_ */
 
