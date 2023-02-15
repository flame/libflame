/* netlib/xerbla_array.f -- translated by f2c (version 20100827).
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

/* > \brief \b XERBLA_ARRAY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download XERBLA_ARRAY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla_
array.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/xerbla_
array.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/xerbla_
array.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE XERBLA_ARRAY( SRNAME_ARRAY, SRNAME_LEN, INFO) */

/*       .. Scalar Arguments .. */
/*       INTEGER SRNAME_LEN, INFO */
/*       .. */
/*       .. Array Arguments .. */
/*       CHARACTER    SRNAME_ARRAY(SRNAME_LEN) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > XERBLA_ARRAY assists other languages in calling XERBLA, the LAPACK */
/* > and BLAS error handler.  Rather than taking a Fortran string argument */
/* > as the function's name, XERBLA_ARRAY takes an array of single */
/* > characters along with the array's length.  XERBLA_ARRAY then copies */
/* > up to 32 characters of that array into a Fortran string and passes */
/* > that to XERBLA.  If called with a non-positive SRNAME_LEN, */
/* > XERBLA_ARRAY will call XERBLA with a string of all blank characters. */
/* > */
/* > Say some macro or other device makes XERBLA_ARRAY available to C99 */
/* > by a name lapack_xerbla and with a common Fortran calling convention. */
/* > Then a C99 program could invoke XERBLA via: */
/* >    { */
/* >      int flen = strlen(__func__); */
/* >      lapack_xerbla(__func__, &flen, &info); */
/* >    } */
/* > */
/* > Providing XERBLA_ARRAY is not necessary for intercepting LAPACK */
/* > errors.  XERBLA_ARRAY calls XERBLA. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SRNAME_ARRAY */
/* > \verbatim */
/* >          SRNAME_ARRAY is CHARACTER    array, dimension (SRNAME_LEN) */
/* >          The name of the routine which called XERBLA_ARRAY. */
/* > \endverbatim */
/* > */
/* > \param[in] SRNAME_LEN */
/* > \verbatim */
/* >          SRNAME_LEN is INTEGER */
/* >          The length of the name in SRNAME_ARRAY. */
/* > \endverbatim */
/* > */
/* > \param[in] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          The position of the invalid parameter in the parameter list */
/* >          of the calling routine. */
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
/* Subroutine */ int xerbla_array__(char *srname_array__, integer *
                                    srname_len__, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    /* Subroutine */
    int s_copy(char *, char *);
    integer i_len(char *, ftnlen);

    /* Local variables */
    integer i__;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    char srname[32];


    /*  -- LAPACK auxiliary routine (version 3.4.0) -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2011 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /* ===================================================================== */

    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Local Arrays .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. Executable Statements .. */
    /* Parameter adjustments */
    --srname_array__;

    /* Function Body */
    s_copy(srname, "");
    /* Computing MIN */
    i__2 = *srname_len__, i__3 = i_len(srname, (ftnlen)32);
    i__1 = min(i__2,i__3);
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        *(unsigned char *)&srname[i__ - 1] = *(unsigned char *)&
                                             srname_array__[i__];
    }
    xerbla_(srname, info);
    return 0;
} /* xerbla_array__ */

