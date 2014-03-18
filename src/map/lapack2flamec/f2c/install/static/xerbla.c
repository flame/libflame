/* xerbla.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "FLA_f2c.h"
#include "stdio.h"

/* Table of constant values */

/* Subroutine */ int xerbla_(char *srname, integer *info)
{
    /*  -- LAPACK auxiliary routine (preliminary version) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
    /*     Courant Institute, Argonne National Lab, and Rice University */
    /*     February 29, 1992 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  XERBLA  is an error handler for the LAPACK routines. */
    /*  It is called by an LAPACK routine if an input parameter has an */
    /*  invalid value.  A message is printed and execution stops. */

    /*  Installers may consider modifying the STOP statement in order to */
    /*  call system-specific exception-handling facilities. */

    /*  Arguments */
    /*  ========= */

    /*  SRNAME  (input) CHARACTER*6 */
    /*          The name of the routine which called XERBLA. */

    /*  INFO    (input) INTEGER */
    /*          The position of the invalid parameter in the parameter list */
    /*          of the calling routine. */

    printf("lapack2flame: On entry to %6s, parameter number %2i had an illegal value\n",
           srname, (int)*info);

    /*     End of XERBLA */

    return 0;
} /* xerbla_ */


