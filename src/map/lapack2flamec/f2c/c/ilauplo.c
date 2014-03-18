/* ../netlib/ilauplo.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ILAUPLO */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ILAUPLO + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilauplo .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilauplo .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilauplo .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* INTEGER FUNCTION ILAUPLO( UPLO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine translated from a character string specifying a */
/* > upper- or lower-triangular matrix to the relevant BLAST-specified */
/* > integer constant. */
/* > */
/* > ILAUPLO returns an INTEGER. If ILAUPLO < 0, then the input is not */
/* > a character indicating an upper- or lower-triangular matrix. */
/* > Otherwise ILAUPLO returns the constant value corresponding to UPLO. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup auxOTHERcomputational */
/* ===================================================================== */
integer ilauplo_(char *uplo)
{
    /* System generated locals */
    integer ret_val;
    /* Local variables */
    extern logical lsame_(char *, char *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    if (lsame_(uplo, "U"))
    {
        ret_val = 121;
    }
    else if (lsame_(uplo, "L"))
    {
        ret_val = 122;
    }
    else
    {
        ret_val = -1;
    }
    return ret_val;
    /* End of ILAUPLO */
}
/* ilauplo_ */
