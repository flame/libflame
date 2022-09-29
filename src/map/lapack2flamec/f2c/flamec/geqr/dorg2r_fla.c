/* ../netlib/dorg2r.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
#ifdef FLA_ENABLE_AMD_OPT
#include "immintrin.h"
#endif

static integer c__1 = 1;
/* > \brief \b DORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by s geqrf (unblocked algorithm). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DORG2R + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2r. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2r. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2r. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, K, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORG2R generates an m by n real matrix Q with orthonormal columns, */
/* > which is defined as the first n columns of a product of k elementary */
/* > reflectors of order m */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by DGEQRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix Q. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix Q. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines the */
/* > matrix Q. N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the i-th column must contain the vector which */
/* > defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* > returned by DGEQRF in the first k columns of its array */
/* > argument A. */
/* > On exit, the m-by-n matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The first dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by DGEQRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int dorg2r_fla(integer *m, integer *n, integer *k, doublereal * a, integer *lda, doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
#ifdef FLA_ENABLE_AMD_OPT
    integer i;
    doublereal *dx;
    __m256d alphav, x0v;
#endif
    /* Local variables */
    integer i__, j, l;
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *), dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *), xerbla_(char *, integer *);
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0 || *n > *m)
    {
        *info = -2;
    }
    else if (*k < 0 || *k > *n)
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DORG2R", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n <= 0)
    {
        return 0;
    }
    /* Initialise columns k+1:n to columns of the unit matrix */
    i__1 = *n;
    for (j = *k + 1;
            j <= i__1;
            ++j)
    {
        i__2 = *m;
        for (l = 1;
                l <= i__2;
                ++l)
        {
            a[l + j * a_dim1] = 0.;
            /* L10: */
        }
        a[j + j * a_dim1] = 1.;
        /* L20: */
    }

    for (i__ = *k;
            i__ >= 1;
            --i__)
    {
        /* Apply H(i) to A(i:m,i:n) from the left */
        if (i__ < *n)
        {
            a[i__ + i__ * a_dim1] = 1.;
            i__1 = *m - i__ + 1;
            i__2 = *n - i__;
            dlarf_("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[ i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
        }

#ifdef FLA_ENABLE_AMD_OPT
        /* Inline DSCAL for small size */
        if (i__ < *m && *m <= FLA_DSCAL_INLINE_SMALL)
        {
            i__1 = *m - i__;
            d__1 = -tau[i__];
            dx = &a[i__ + i__ * a_dim1];

            /* Load scaling factor */
            alphav = _mm256_set1_pd(d__1);

            /* Scaling with 0 */
            if(d__1 == 0.0)
            {
                for (i = 1; i <= (i__1-3); i+=4 )
                {
                    dx[i] = 0.0;
                    dx[i+1] = 0.0;
                    dx[i+2] = 0.0;
                    dx[i+3] = 0.0;
                }
                for (; i <= i__1; ++i)
                {
                    dx[i] = 0.0;
                }
            }
            /* Scaling factor other than 0 */
            else
            {
                for ( i = 1; i <= (i__1 - 3); i += 4)
                {
                    /* Load the input values */
                    x0v = _mm256_loadu_pd((double const *) &dx[i]);

                    /* perform alpha * x  */
                    x0v = _mm256_mul_pd( alphav, x0v );  

                    /* Store the output */
                    _mm256_storeu_pd((double *) &dx[i], x0v);
                }

                /* Remainder iterations */
                if((i__1-i) >= 2)
                {
                    for ( ; i <= (i__1-1); i += 2 )
                    {
                        dx[i] *= d__1;
                        dx[i+1] *= d__1;
                    }
                    for ( ; i <= i__1; ++i )
                    {
                        dx[i] *= d__1;
                    }
                }
                else
                {
                    for ( ; i <= i__1; ++i )
                    {
                        dx[i] *= d__1;
                    }
                }
            }
        }
        else
        {
            i__1 = *m - i__;
            d__1 = -tau[i__];
            dscal_(&i__1, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
        }
#else
        if (i__ < *m)
        {
            i__1 = *m - i__;
            d__1 = -tau[i__];
            dscal_(&i__1, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
        }
#endif
        a[i__ + i__ * a_dim1] = 1. - tau[i__];
        /* Set A(1:i-1,i) to zero */
        i__1 = i__ - 1;
        for (l = 1;
                l <= i__1;
                ++l)
        {
            a[l + i__ * a_dim1] = 0.;
            /* L30: */
        }
        /* L40: */
    }
    return 0;
    /* End of DORG2R */
}
/* dorg2r_ */
