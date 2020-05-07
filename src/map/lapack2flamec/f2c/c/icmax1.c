/* ../netlib/v3.9.0/icmax1.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* > \brief \b ICMAX1 finds the index of the first vector element of maximum absolute value. */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download ICMAX1 + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/icmax1. f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/icmax1. f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/icmax1. f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* INTEGER FUNCTION ICMAX1( N, CX, INCX ) */
 /* .. Scalar Arguments .. */
 /* INTEGER INCX, N */
 /* .. */
 /* .. Array Arguments .. */
 /* COMPLEX CX( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > ICMAX1 finds the index of the first vector element of maximum absolute value. */
 /* > */
 /* > Based on ICAMAX from Level 1 BLAS. */
 /* > The change is to use the 'genuine' absolute value. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] N */
 /* > \verbatim */
 /* > N is INTEGER */
 /* > The number of elements in the vector CX. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] CX */
 /* > \verbatim */
 /* > CX is COMPLEX array, dimension (N) */
 /* > The vector CX. The ICMAX1 function returns the index of its first */
 /* > element of maximum absolute value. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] INCX */
 /* > \verbatim */
 /* > INCX is INTEGER */
 /* > The spacing between successive values of CX. INCX >= 1. */
 /* > \endverbatim */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date February 2014 */
 /* > \ingroup complexOTHERauxiliary */
 /* > \par Contributors: */
 /* ================== */
 /* > */
 /* > Nick Higham for use with CLACON. */
 /* ===================================================================== */
 integer icmax1_(integer *n, complex *cx, integer *incx) {
 /* System generated locals */
 integer ret_val, i__1;
 /* Builtin functions */
 double c_abs(complex *);
 /* Local variables */
 integer i__, ix;
 real smax;
 /* -- LAPACK auxiliary routine (version 3.7.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* February 2014 */
 /* .. Scalar Arguments .. */
 /* .. */
 /* .. Array Arguments .. */
 /* .. */
 /* ===================================================================== */
 /* .. Local Scalars .. */
 /* .. */
 /* .. Intrinsic Functions .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Parameter adjustments */
 --cx;
 /* Function Body */
 ret_val = 0;
 if (*n < 1 || *incx <= 0) {
 return ret_val;
 }
 ret_val = 1;
 if (*n == 1) {
 return ret_val;
 }
 if (*incx == 1) {
 /* code for increment equal to 1 */
 smax = c_abs(&cx[1]);
 i__1 = *n;
 for (i__ = 2;
 i__ <= i__1;
 ++i__) {
 if (c_abs(&cx[i__]) > smax) {
 ret_val = i__;
 smax = c_abs(&cx[i__]);
 }
 }
 }
 else {
 /* code for increment not equal to 1 */
 ix = 1;
 smax = c_abs(&cx[1]);
 ix += *incx;
 i__1 = *n;
 for (i__ = 2;
 i__ <= i__1;
 ++i__) {
 if (c_abs(&cx[ix]) > smax) {
 ret_val = i__;
 smax = c_abs(&cx[ix]);
 }
 ix += *incx;
 }
 }
 return ret_val;
 /* End of ICMAX1 */
 }
 /* icmax1_ */
 