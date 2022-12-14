/* slarmm.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLARMM */
/* Definition: */
/* =========== */
/* REAL FUNCTION SLARMM( ANORM, BNORM, CNORM ) */
/* .. Scalar Arguments .. */
/* REAL ANORM, BNORM, CNORM */
/* .. */
/* > \par Purpose: */
/* ======= */
/* > */
/* > \verbatim */
/* > */
/* > SLARMM returns a factor s in (0, 1] such that the linear updates */
/* > */
/* > (s * C) - A * (s * B) and (s * C) - (s * A) * B */
/* > */
/* > cannot overflow, where A, B, and C are matrices of conforming */
/* > dimensions. */
/* > */
/* > This is an auxiliary routine so there is no argument checking. */
/* > \endverbatim */
/* Arguments: */
/* ========= */
/* > \param[in] ANORM */
/* > \verbatim */
/* > ANORM is REAL */
/* > The infinity norm of A. ANORM >= 0. */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] BNORM */
/* > \verbatim */
/* > BNORM is REAL */
/* > The infinity norm of B. BNORM >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] CNORM */
/* > \verbatim */
/* > CNORM is REAL */
/* > The infinity norm of C. CNORM >= 0. */
/* > \endverbatim */
/* > */
/* > */
/* ===================================================================== */
/* > References: */
/* > C. C. Kjelgaard Mikkelsen and L. Karlsson, Blocked Algorithms for */
/* > Robust Solution of Triangular Linear Systems. In: International */
/* > Conference on Parallel Processing and Applied Mathematics, pages */
/* > 68--78. Springer, 2017. */
/* > */
/* > \ingroup OTHERauxiliary */
/* ===================================================================== */
real slarmm_(real *anorm, real *bnorm, real *cnorm)
{
    AOCL_DTL_TRACE_LOG_INIT
    /* System generated locals */
    doublereal ret_val;
    /* Local variables */
    extern real slamch_(char *);
    real bignum, smlnum;
    /* .. Scalar Arguments .. */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Determine machine dependent parameters to control overflow. */
    smlnum = slamch_("Safe minimum") / slamch_("Precision");
    bignum = 1.f / smlnum / 4.f;
    /* Compute a scale factor. */
    ret_val = 1.f;
    if (*bnorm <= 1.f)
    {
        if (*anorm * *bnorm > bignum - *cnorm)
        {
            ret_val = .5f;
        }
    }
    else
    {
        if (*anorm > (bignum - *cnorm) / *bnorm)
        {
            ret_val = .5f / *bnorm;
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return ret_val;
    /* ==== End of SLARMM ==== */
}
/* slarmm_ */
