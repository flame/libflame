/* ../netlib/zlanhp.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZLANHP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele ment of largest absolute value of a complex Hermitian matrix supplied in packed form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLANHP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanhp. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanhp. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanhp. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* DOUBLE PRECISION FUNCTION ZLANHP( NORM, UPLO, N, AP, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM, UPLO */
/* INTEGER N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION WORK( * ) */
/* COMPLEX*16 AP( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLANHP returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > complex hermitian matrix A, supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return ZLANHP */
/* > \verbatim */
/* > */
/* > ZLANHP = ( max(f2c_abs(A(i,j))), NORM = 'M' or 'm' */
/* > ( */
/* > ( norm1(A), NORM = '1', 'O' or 'o' */
/* > ( */
/* > ( normI(A), NORM = 'I' or 'i' */
/* > ( */
/* > ( normF(A), NORM = 'F', 'f', 'E' or 'e' */
/* > */
/* > where norm1 denotes the one norm of a matrix (maximum column sum), */
/* > normI denotes the infinity norm of a matrix (maximum row sum) and */
/* > normF denotes the Frobenius norm of a matrix (square root of sum of */
/* > squares). Note that max(f2c_abs(A(i,j))) is not a consistent matrix norm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] NORM */
/* > \verbatim */
/* > NORM is CHARACTER*1 */
/* > Specifies the value to be returned in ZLANHP as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > hermitian matrix A is supplied. */
/* > = 'U': Upper triangular part of A is supplied */
/* > = 'L': Lower triangular part of A is supplied */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. When N = 0, ZLANHP is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > The upper or lower triangle of the hermitian matrix A, packed */
/* > columnwise in a linear array. The j-th column of A is stored */
/* > in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > Note that the imaginary parts of the diagonal elements need */
/* > not be set and are assumed to be zero. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
/* > where LWORK >= N when NORM = 'I' or '1' or 'O';
otherwise, */
/* > WORK is not referenced. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
doublereal zlanhp_(char *norm, char *uplo, integer *n, doublecomplex *ap, doublereal *work)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;
    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);
    /* Local variables */
    integer i__, j, k;
    doublereal sum, absa, scale;
    extern logical lsame_(char *, char *);
    doublereal value;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */
    int zlassq_(integer *, doublecomplex *, integer *, doublereal *, doublereal *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --work;
    --ap;
    /* Function Body */
    if (*n == 0)
    {
        value = 0.;
    }
    else if (lsame_(norm, "M"))
    {
        /* Find max(f2c_abs(A(i,j))). */
        value = 0.;
        if (lsame_(uplo, "U"))
        {
            k = 0;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = k + j - 1;
                for (i__ = k + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    sum = z_abs(&ap[i__]);
                    if (value < sum || disnan_(&sum))
                    {
                        value = sum;
                    }
                    /* L10: */
                }
                k += j;
                i__2 = k;
                sum = (d__1 = ap[i__2].r, f2c_abs(d__1));
                if (value < sum || disnan_(&sum))
                {
                    value = sum;
                }
                /* L20: */
            }
        }
        else
        {
            k = 1;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = k;
                sum = (d__1 = ap[i__2].r, f2c_abs(d__1));
                if (value < sum || disnan_(&sum))
                {
                    value = sum;
                }
                i__2 = k + *n - j;
                for (i__ = k + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    sum = z_abs(&ap[i__]);
                    if (value < sum || disnan_(&sum))
                    {
                        value = sum;
                    }
                    /* L30: */
                }
                k = k + *n - j + 1;
                /* L40: */
            }
        }
    }
    else if (lsame_(norm, "I") || lsame_(norm, "O") || *(unsigned char *)norm == '1')
    {
        /* Find normI(A) ( = norm1(A), since A is hermitian). */
        value = 0.;
        k = 1;
        if (lsame_(uplo, "U"))
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                sum = 0.;
                i__2 = j - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    absa = z_abs(&ap[k]);
                    sum += absa;
                    work[i__] += absa;
                    ++k;
                    /* L50: */
                }
                i__2 = k;
                work[j] = sum + (d__1 = ap[i__2].r, f2c_abs(d__1));
                ++k;
                /* L60: */
            }
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                sum = work[i__];
                if (value < sum || disnan_(&sum))
                {
                    value = sum;
                }
                /* L70: */
            }
        }
        else
        {
            i__1 = *n;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                work[i__] = 0.;
                /* L80: */
            }
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = k;
                sum = work[j] + (d__1 = ap[i__2].r, f2c_abs(d__1));
                ++k;
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    absa = z_abs(&ap[k]);
                    sum += absa;
                    work[i__] += absa;
                    ++k;
                    /* L90: */
                }
                if (value < sum || disnan_(&sum))
                {
                    value = sum;
                }
                /* L100: */
            }
        }
    }
    else if (lsame_(norm, "F") || lsame_(norm, "E"))
    {
        /* Find normF(A). */
        scale = 0.;
        sum = 1.;
        k = 2;
        if (lsame_(uplo, "U"))
        {
            i__1 = *n;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                i__2 = j - 1;
                zlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
                k += j;
                /* L110: */
            }
        }
        else
        {
            i__1 = *n - 1;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = *n - j;
                zlassq_(&i__2, &ap[k], &c__1, &scale, &sum);
                k = k + *n - j + 1;
                /* L120: */
            }
        }
        sum *= 2;
        k = 1;
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = k;
            if (ap[i__2].r != 0.)
            {
                i__2 = k;
                absa = (d__1 = ap[i__2].r, f2c_abs(d__1));
                if (scale < absa)
                {
                    /* Computing 2nd power */
                    d__1 = scale / absa;
                    sum = sum * (d__1 * d__1) + 1.;
                    scale = absa;
                }
                else
                {
                    /* Computing 2nd power */
                    d__1 = absa / scale;
                    sum += d__1 * d__1;
                }
            }
            if (lsame_(uplo, "U"))
            {
                k = k + i__ + 1;
            }
            else
            {
                k = k + *n - i__ + 1;
            }
            /* L130: */
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    return ret_val;
    /* End of ZLANHP */
}
/* zlanhp_ */
