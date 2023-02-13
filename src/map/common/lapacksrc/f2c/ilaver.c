/* ../netlib/ilaver.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ILAVER returns the LAPACK version. */
/* * */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* Definition: */
/* =========== */
/* SUBROUTINE ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH ) */
/* INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine returns the LAPACK version. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[out] VERS_MAJOR */
/* > return the lapack major version */
/* > */
/* > \param[out] VERS_MINOR */
/* > return the lapack minor version from the major version */
/* > */
/* > \param[out] VERS_PATCH */
/* > return the lapack patch version from the minor version */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int ilaver_(integer *vers_major__, integer *vers_minor__, integer *vers_patch__)
{
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
    /* ===================================================================== */
    /* ===================================================================== */
    *vers_major__ = 3;
    *vers_minor__ = 5;
    *vers_patch__ = 0;
    /* ===================================================================== */
    return 0;
}
/* ilaver_ */
