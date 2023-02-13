/* ../netlib/dopgtr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DOPGTR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DOPGTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dopgtr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dopgtr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dopgtr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDQ, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION AP( * ), Q( LDQ, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DOPGTR generates a real orthogonal matrix Q which is defined as the */
/* > product of n-1 elementary reflectors H(i) of order n, as returned by */
/* > DSPTRD using packed storage: */
/* > */
/* > if UPLO = 'U', Q = H(n-1) . . . H(2) H(1), */
/* > */
/* > if UPLO = 'L', Q = H(1) H(2) . . . H(n-1). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangular packed storage used in previous */
/* > call to DSPTRD;
*/
/* > = 'L': Lower triangular packed storage used in previous */
/* > call to DSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix Q. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* > AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* > The vectors which define the elementary reflectors, as */
/* > returned by DSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (N-1) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by DSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* > Q is DOUBLE PRECISION array, dimension (LDQ,N) */
/* > The N-by-N orthogonal matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (N-1) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int dopgtr_(char *uplo, integer *n, doublereal *ap, doublereal *tau, doublereal *q, integer *ldq, doublereal *work, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, j, ij;
    extern logical lsame_(char *, char *);
    integer iinfo;
    logical upper;
    extern /* Subroutine */
    int dorg2l_(integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), dorg2r_fla(integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    --ap;
    --tau;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*ldq < max(1,*n))
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DOPGTR", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    if (upper)
    {
        /* Q was determined by a call to DSPTRD with UPLO = 'U' */
        /* Unpack the vectors which define the elementary reflectors and */
        /* set the last row and column of Q equal to those of the unit */
        /* matrix */
        ij = 2;
        i__1 = *n - 1;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j - 1;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                q[i__ + j * q_dim1] = ap[ij];
                ++ij;
                /* L10: */
            }
            ij += 2;
            q[*n + j * q_dim1] = 0.;
            /* L20: */
        }
        i__1 = *n - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            q[i__ + *n * q_dim1] = 0.;
            /* L30: */
        }
        q[*n + *n * q_dim1] = 1.;
        /* Generate Q(1:n-1,1:n-1) */
        i__1 = *n - 1;
        i__2 = *n - 1;
        i__3 = *n - 1;
        dorg2l_(&i__1, &i__2, &i__3, &q[q_offset], ldq, &tau[1], &work[1], & iinfo);
    }
    else
    {
        /* Q was determined by a call to DSPTRD with UPLO = 'L'. */
        /* Unpack the vectors which define the elementary reflectors and */
        /* set the first row and column of Q equal to those of the unit */
        /* matrix */
        q[q_dim1 + 1] = 1.;
        i__1 = *n;
        for (i__ = 2;
                i__ <= i__1;
                ++i__)
        {
            q[i__ + q_dim1] = 0.;
            /* L40: */
        }
        ij = 3;
        i__1 = *n;
        for (j = 2;
                j <= i__1;
                ++j)
        {
            q[j * q_dim1 + 1] = 0.;
            i__2 = *n;
            for (i__ = j + 1;
                    i__ <= i__2;
                    ++i__)
            {
                q[i__ + j * q_dim1] = ap[ij];
                ++ij;
                /* L50: */
            }
            ij += 2;
            /* L60: */
        }
        if (*n > 1)
        {
            /* Generate Q(2:n,2:n) */
            i__1 = *n - 1;
            i__2 = *n - 1;
            i__3 = *n - 1;
            dorg2r_fla(&i__1, &i__2, &i__3, &q[(q_dim1 << 1) + 2], ldq, &tau[1], &work[1], &iinfo);
        }
    }
    return 0;
    /* End of DOPGTR */
}
/* dopgtr_ */
