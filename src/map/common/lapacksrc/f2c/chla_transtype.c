/* ../netlib/chla_transtype.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CHLA_TRANSTYPE */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHLA_TRANSTYPE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chla_tr anstype.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chla_tr anstype.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chla_tr anstype.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* CHARACTER*1 FUNCTION CHLA_TRANSTYPE( TRANS ) */
/* .. Scalar Arguments .. */
/* INTEGER TRANS */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine translates from a BLAST-specified integer constant to */
/* > the character string specifying a transposition operation. */
/* > */
/* > CHLA_TRANSTYPE returns an CHARACTER*1. If CHLA_TRANSTYPE is 'X', */
/* > then input is not an integer indicating a transposition operator. */
/* > Otherwise CHLA_TRANSTYPE returns the constant value corresponding to */
/* > TRANS. */
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
/* Character */
VOID chla_transtype_(char *ret_val, integer *trans)
{
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Executable Statements .. */
    if (*trans == 111)
    {
        *(unsigned char *)ret_val = 'N';
    }
    else if (*trans == 112)
    {
        *(unsigned char *)ret_val = 'T';
    }
    else if (*trans == 113)
    {
        *(unsigned char *)ret_val = 'C';
    }
    else
    {
        *(unsigned char *)ret_val = 'X';
    }
    return ;
    /* End of CHLA_TRANSTYPE */
}
/* chla_transtype__ */
