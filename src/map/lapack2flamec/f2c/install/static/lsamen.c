/* netlib/lsamen.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "FLA_f2c.h"

/* > \brief \b LSAMEN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download LSAMEN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/lsamen.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/lsamen.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/lsamen.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       LOGICAL          FUNCTION LSAMEN( N, CA, CB ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER*( * )    CA, CB */
/*       INTEGER            N */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > LSAMEN  tests if the first N letters of CA are the same as the */
/* > first N letters of CB, regardless of case. */
/* > LSAMEN returns .TRUE. if CA and CB are equivalent except for case */
/* > and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA ) */
/* > or LEN( CB ) is less than N. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of characters in CA and CB to be compared. */
/* > \endverbatim */
/* > */
/* > \param[in] CA */
/* > \verbatim */
/* >          CA is CHARACTER*(*) */
/* > \endverbatim */
/* > */
/* > \param[in] CB */
/* > \verbatim */
/* >          CB is CHARACTER*(*) */
/* >          CA and CB specify two character strings of length at least N. */
/* >          Only the first N characters of each string will be accessed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup auxOTHERauxiliary */

/*  ===================================================================== */
logical lsamen_(integer *n, char *ca, char *cb)
{
    /* System generated locals */
    integer i__1;
    logical ret_val;

    /* Builtin functions */
    integer i_len(char *, ftnlen);

    /* Local variables */
    integer i__;
    extern logical lsame_(char *, char *);

    ftnlen ca_len, cb_len;

    ca_len = strlen (ca);
    cb_len = strlen (cb);

    /*  -- LAPACK auxiliary routine (version 3.4.0) -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2011 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /* ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. Executable Statements .. */

    ret_val = FALSE_;
    if (i_len(ca, ca_len) < *n || i_len(cb, cb_len) < *n)
    {
        goto L20;
    }

    /*     Do for each character in the two strings. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {

        /*        Test if the characters are equal using LSAME. */

        if (! lsame_(ca + (i__ - 1), cb + (i__ - 1)))
        {
            goto L20;
        }

        /* L10: */
    }
    ret_val = TRUE_;

L20:
    return ret_val;

    /*     End of LSAMEN */

} /* lsamen_ */

