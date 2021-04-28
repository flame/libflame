/* ../netlib/sisnan.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* > \brief \b SISNAN tests input for NaN. */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download SISNAN + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sisnan. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sisnan. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sisnan. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* LOGICAL FUNCTION SISNAN( SIN ) */
 /* .. Scalar Arguments .. */
 /* REAL, INTENT(IN) :: SIN */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > SISNAN returns .TRUE. if its argument is NaN, and .FALSE. */
 /* > otherwise. To be replaced by the Fortran 2003 intrinsic in the */
 /* > future. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] SIN */
 /* > \verbatim */
 /* > SIN is REAL */
 /* > Input to test for NaN. */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date June 2017 */
 /* > \ingroup OTHERauxiliary */
 /* ===================================================================== */
 logical sisnan_(real *sin__) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
 /* System generated locals */
 logical ret_val;
 /* Local variables */
 extern logical slaisnan_(real *, real *);
 /* -- LAPACK auxiliary routine (version 3.7.1) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* June 2017 */
 /* .. Scalar Arguments .. */
 /* .. */
 /* ===================================================================== */
 /* .. External Functions .. */
 /* .. */
 /* .. Executable Statements .. */
 ret_val = slaisnan_(sin__, sin__);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return ret_val;
 }
 /* sisnan_ */
 
