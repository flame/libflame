/* ../netlib/ssyconv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SSYCONV */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSYCONV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconv .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconv .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconv .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSYCONV( UPLO, WAY, N, A, LDA, IPIV, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO, WAY */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYCONV convert A given by TRF into L and D and vice-versa. */
/* > Get Non-diag elements of D (returned in workspace) and */
/* > apply or reverse permutation done in TRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**T;
*/
/* > = 'L': Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] WAY */
/* > \verbatim */
/* > WAY is CHARACTER*1 */
/* > = 'C': Convert */
/* > = 'R': Revert */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The block diagonal matrix D and the multipliers used to */
/* > obtain the factor U or L as computed by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (N) */
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
/* > \ingroup realSYcomputational */
/* ===================================================================== */
/* Subroutine */
int ssyconv_(char *uplo, char *way, integer *n, real *a, integer *lda, integer *ipiv, real *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    integer i__, j, ip;
    real temp;
    extern logical lsame_(char *, char *);
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    logical convert;
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
    /* .. External Functions .. */
    /* .. External Subroutines .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    convert = lsame_(way, "C");
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (! convert && ! lsame_(way, "R"))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSYCONV", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    if (upper)
    {
        /* A is UPPER */
        /* Convert A (A is upper) */
        /* Convert VALUE */
        if (convert)
        {
            i__ = *n;
            work[1] = 0.f;
            while(i__ > 1)
            {
                if (ipiv[i__] < 0)
                {
                    work[i__] = a[i__ - 1 + i__ * a_dim1];
                    a[i__ - 1 + i__ * a_dim1] = 0.f;
                    --i__;
                }
                else
                {
                    work[i__] = 0.f;
                }
                --i__;
            }
            /* Convert PERMUTATIONS */
            i__ = *n;
            while(i__ >= 1)
            {
                if (ipiv[i__] > 0)
                {
                    ip = ipiv[i__];
                    if (i__ < *n)
                    {
                        i__1 = *n;
                        for (j = i__ + 1;
                                j <= i__1;
                                ++j)
                        {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = temp;
                            /* L12: */
                        }
                    }
                }
                else
                {
                    ip = -ipiv[i__];
                    if (i__ < *n)
                    {
                        i__1 = *n;
                        for (j = i__ + 1;
                                j <= i__1;
                                ++j)
                        {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
                            a[i__ - 1 + j * a_dim1] = temp;
                            /* L13: */
                        }
                    }
                    --i__;
                }
                --i__;
            }
        }
        else
        {
            /* Revert A (A is upper) */
            /* Revert PERMUTATIONS */
            i__ = 1;
            while(i__ <= *n)
            {
                if (ipiv[i__] > 0)
                {
                    ip = ipiv[i__];
                    if (i__ < *n)
                    {
                        i__1 = *n;
                        for (j = i__ + 1;
                                j <= i__1;
                                ++j)
                        {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = temp;
                        }
                    }
                }
                else
                {
                    ip = -ipiv[i__];
                    ++i__;
                    if (i__ < *n)
                    {
                        i__1 = *n;
                        for (j = i__ + 1;
                                j <= i__1;
                                ++j)
                        {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
                            a[i__ - 1 + j * a_dim1] = temp;
                        }
                    }
                }
                ++i__;
            }
            /* Revert VALUE */
            i__ = *n;
            while(i__ > 1)
            {
                if (ipiv[i__] < 0)
                {
                    a[i__ - 1 + i__ * a_dim1] = work[i__];
                    --i__;
                }
                --i__;
            }
        }
    }
    else
    {
        /* A is LOWER */
        if (convert)
        {
            /* Convert A (A is lower) */
            /* Convert VALUE */
            i__ = 1;
            work[*n] = 0.f;
            while(i__ <= *n)
            {
                if (i__ < *n && ipiv[i__] < 0)
                {
                    work[i__] = a[i__ + 1 + i__ * a_dim1];
                    a[i__ + 1 + i__ * a_dim1] = 0.f;
                    ++i__;
                }
                else
                {
                    work[i__] = 0.f;
                }
                ++i__;
            }
            /* Convert PERMUTATIONS */
            i__ = 1;
            while(i__ <= *n)
            {
                if (ipiv[i__] > 0)
                {
                    ip = ipiv[i__];
                    if (i__ > 1)
                    {
                        i__1 = i__ - 1;
                        for (j = 1;
                                j <= i__1;
                                ++j)
                        {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = temp;
                            /* L22: */
                        }
                    }
                }
                else
                {
                    ip = -ipiv[i__];
                    if (i__ > 1)
                    {
                        i__1 = i__ - 1;
                        for (j = 1;
                                j <= i__1;
                                ++j)
                        {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ + 1 + j * a_dim1];
                            a[i__ + 1 + j * a_dim1] = temp;
                            /* L23: */
                        }
                    }
                    ++i__;
                }
                ++i__;
            }
        }
        else
        {
            /* Revert A (A is lower) */
            /* Revert PERMUTATIONS */
            i__ = *n;
            while(i__ >= 1)
            {
                if (ipiv[i__] > 0)
                {
                    ip = ipiv[i__];
                    if (i__ > 1)
                    {
                        i__1 = i__ - 1;
                        for (j = 1;
                                j <= i__1;
                                ++j)
                        {
                            temp = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = temp;
                        }
                    }
                }
                else
                {
                    ip = -ipiv[i__];
                    --i__;
                    if (i__ > 1)
                    {
                        i__1 = i__ - 1;
                        for (j = 1;
                                j <= i__1;
                                ++j)
                        {
                            temp = a[i__ + 1 + j * a_dim1];
                            a[i__ + 1 + j * a_dim1] = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = temp;
                        }
                    }
                }
                --i__;
            }
            /* Revert VALUE */
            i__ = 1;
            while(i__ <= *n - 1)
            {
                if (ipiv[i__] < 0)
                {
                    a[i__ + 1 + i__ * a_dim1] = work[i__];
                    ++i__;
                }
                ++i__;
            }
        }
    }
    return 0;
    /* End of SSYCONV */
}
/* ssyconv_ */
