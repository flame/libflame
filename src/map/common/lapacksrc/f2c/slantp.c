/* ../netlib/slantp.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SLANTP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele ment of largest absolute value of a triangular matrix supplied in packed form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLANTP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantp. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantp. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantp. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION SLANTP( NORM, UPLO, DIAG, N, AP, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORM, UPLO */
/* INTEGER N */
/* .. */
/* .. Array Arguments .. */
/* REAL AP( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANTP returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > triangular matrix A, supplied in packed form. */
/* > \endverbatim */
/* > */
/* > \return SLANTP */
/* > \verbatim */
/* > */
/* > SLANTP = ( max(f2c_abs(A(i,j))), NORM = 'M' or 'm' */
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
/* > Specifies the value to be returned in SLANTP as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the matrix A is upper or lower triangular. */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* > DIAG is CHARACTER*1 */
/* > Specifies whether or not the matrix A is unit triangular. */
/* > = 'N': Non-unit triangular */
/* > = 'U': Unit triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. When N = 0, SLANTP is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* > AP is REAL array, dimension (N*(N+1)/2) */
/* > The upper or lower triangular matrix A, packed columnwise in */
/* > a linear array. The j-th column of A is stored in the array */
/* > AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > Note that when DIAG = 'U', the elements of the array AP */
/* > corresponding to the diagonal elements of the matrix A are */
/* > not referenced, but are assumed to be one. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)), */
/* > where LWORK >= N when NORM = 'I';
otherwise, WORK is not */
/* > referenced. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* ===================================================================== */
real slantp_(char *norm, char *uplo, char *diag, integer *n, real *ap, real * work)
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k;
    real sum, scale;
    logical udiag;
    extern logical lsame_(char *, char *);
    real value;
    extern logical sisnan_(real *);
    extern /* Subroutine */
    int slassq_(integer *, real *, integer *, real *, real *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
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
        value = 0.f;
    }
    else if (lsame_(norm, "M"))
    {
        /* Find max(f2c_abs(A(i,j))). */
        k = 1;
        if (lsame_(diag, "U"))
        {
            value = 1.f;
            if (lsame_(uplo, "U"))
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = k + j - 2;
                    for (i__ = k;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum = (r__1 = ap[i__], f2c_abs(r__1));
                        if (value < sum || sisnan_(&sum))
                        {
                            value = sum;
                        }
                        /* L10: */
                    }
                    k += j;
                    /* L20: */
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = k + *n - j;
                    for (i__ = k + 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum = (r__1 = ap[i__], f2c_abs(r__1));
                        if (value < sum || sisnan_(&sum))
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
        else
        {
            value = 0.f;
            if (lsame_(uplo, "U"))
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = k + j - 1;
                    for (i__ = k;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum = (r__1 = ap[i__], f2c_abs(r__1));
                        if (value < sum || sisnan_(&sum))
                        {
                            value = sum;
                        }
                        /* L50: */
                    }
                    k += j;
                    /* L60: */
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = k + *n - j;
                    for (i__ = k;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum = (r__1 = ap[i__], f2c_abs(r__1));
                        if (value < sum || sisnan_(&sum))
                        {
                            value = sum;
                        }
                        /* L70: */
                    }
                    k = k + *n - j + 1;
                    /* L80: */
                }
            }
        }
    }
    else if (lsame_(norm, "O") || *(unsigned char *) norm == '1')
    {
        /* Find norm1(A). */
        value = 0.f;
        k = 1;
        udiag = lsame_(diag, "U");
        if (lsame_(uplo, "U"))
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (udiag)
                {
                    sum = 1.f;
                    i__2 = k + j - 2;
                    for (i__ = k;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += (r__1 = ap[i__], f2c_abs(r__1));
                        /* L90: */
                    }
                }
                else
                {
                    sum = 0.f;
                    i__2 = k + j - 1;
                    for (i__ = k;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += (r__1 = ap[i__], f2c_abs(r__1));
                        /* L100: */
                    }
                }
                k += j;
                if (value < sum || sisnan_(&sum))
                {
                    value = sum;
                }
                /* L110: */
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (udiag)
                {
                    sum = 1.f;
                    i__2 = k + *n - j;
                    for (i__ = k + 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += (r__1 = ap[i__], f2c_abs(r__1));
                        /* L120: */
                    }
                }
                else
                {
                    sum = 0.f;
                    i__2 = k + *n - j;
                    for (i__ = k;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += (r__1 = ap[i__], f2c_abs(r__1));
                        /* L130: */
                    }
                }
                k = k + *n - j + 1;
                if (value < sum || sisnan_(&sum))
                {
                    value = sum;
                }
                /* L140: */
            }
        }
    }
    else if (lsame_(norm, "I"))
    {
        /* Find normI(A). */
        k = 1;
        if (lsame_(uplo, "U"))
        {
            if (lsame_(diag, "U"))
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] = 1.f;
                    /* L150: */
                }
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = j - 1;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        work[i__] += (r__1 = ap[k], f2c_abs(r__1));
                        ++k;
                        /* L160: */
                    }
                    ++k;
                    /* L170: */
                }
            }
            else
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] = 0.f;
                    /* L180: */
                }
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = j;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        work[i__] += (r__1 = ap[k], f2c_abs(r__1));
                        ++k;
                        /* L190: */
                    }
                    /* L200: */
                }
            }
        }
        else
        {
            if (lsame_(diag, "U"))
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] = 1.f;
                    /* L210: */
                }
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    ++k;
                    i__2 = *n;
                    for (i__ = j + 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        work[i__] += (r__1 = ap[k], f2c_abs(r__1));
                        ++k;
                        /* L220: */
                    }
                    /* L230: */
                }
            }
            else
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] = 0.f;
                    /* L240: */
                }
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *n;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        work[i__] += (r__1 = ap[k], f2c_abs(r__1));
                        ++k;
                        /* L250: */
                    }
                    /* L260: */
                }
            }
        }
        value = 0.f;
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            sum = work[i__];
            if (value < sum || sisnan_(&sum))
            {
                value = sum;
            }
            /* L270: */
        }
    }
    else if (lsame_(norm, "F") || lsame_(norm, "E"))
    {
        /* Find normF(A). */
        if (lsame_(uplo, "U"))
        {
            if (lsame_(diag, "U"))
            {
                scale = 1.f;
                sum = (real) (*n);
                k = 2;
                i__1 = *n;
                for (j = 2;
                        j <= i__1;
                        ++j)
                {
                    i__2 = j - 1;
                    slassq_(&i__2, &ap[k], &c__1, &scale, &sum);
                    k += j;
                    /* L280: */
                }
            }
            else
            {
                scale = 0.f;
                sum = 1.f;
                k = 1;
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    slassq_(&j, &ap[k], &c__1, &scale, &sum);
                    k += j;
                    /* L290: */
                }
            }
        }
        else
        {
            if (lsame_(diag, "U"))
            {
                scale = 1.f;
                sum = (real) (*n);
                k = 2;
                i__1 = *n - 1;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *n - j;
                    slassq_(&i__2, &ap[k], &c__1, &scale, &sum);
                    k = k + *n - j + 1;
                    /* L300: */
                }
            }
            else
            {
                scale = 0.f;
                sum = 1.f;
                k = 1;
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *n - j + 1;
                    slassq_(&i__2, &ap[k], &c__1, &scale, &sum);
                    k = k + *n - j + 1;
                    /* L310: */
                }
            }
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    return ret_val;
    /* End of SLANTP */
}
/* slantp_ */
