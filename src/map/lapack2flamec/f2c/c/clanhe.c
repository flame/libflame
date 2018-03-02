/* ../netlib/clanhe.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLANHE returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele ment of largest absolute value of a complex Hermitian matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLANHE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhe. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhe. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhe. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLANHE( NORM, UPLO, N, A, LDA, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM, UPLO */
/* INTEGER LDA, N */
/* .. */
/* .. Array Arguments .. */
/* REAL WORK( * ) */
/* COMPLEX A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANHE returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > complex hermitian matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANHE */
/* > \verbatim */
/* > */
/* > CLANHE = ( max(f2c_abs(A(i,j))), NORM = 'M' or 'm' */
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
/* > Specifies the value to be returned in CLANHE as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > hermitian matrix A is to be referenced. */
/* > = 'U': Upper triangular part of A is referenced */
/* > = 'L': Lower triangular part of A is referenced */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. When N = 0, CLANHE is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > The hermitian matrix A. If UPLO = 'U', the leading n by n */
/* > upper triangular part of A contains the upper triangular part */
/* > of the matrix A, and the strictly lower triangular part of A */
/* > is not referenced. If UPLO = 'L', the leading n by n lower */
/* > triangular part of A contains the lower triangular part of */
/* > the matrix A, and the strictly upper triangular part of A is */
/* > not referenced. Note that the imaginary parts of the diagonal */
/* > elements need not be set and are assumed to be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(N,1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)), */
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
/* > \ingroup complexHEauxiliary */
/* ===================================================================== */
real clanhe_(char *norm, char *uplo, integer *n, complex *a, integer *lda, real *work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real ret_val, r__1;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    integer i__, j;
    real sum, absa, scale;
    extern logical lsame_(char *, char *);
    real value;
    extern /* Subroutine */
    int classq_(integer *, complex *, integer *, real *, real *);
    extern logical sisnan_(real *);
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    /* Function Body */
    if (*n == 0)
    {
        value = 0.f;
    }
    else if (lsame_(norm, "M"))
    {
        /* Find max(f2c_abs(A(i,j))). */
        value = 0.f;
        if (lsame_(uplo, "U"))
        {
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
                    sum = c_abs(&a[i__ + j * a_dim1]);
                    if (value < sum || sisnan_(&sum))
                    {
                        value = sum;
                    }
                    /* L10: */
                }
                i__2 = j + j * a_dim1;
                sum = (r__1 = a[i__2].r, f2c_abs(r__1));
                if (value < sum || sisnan_(&sum))
                {
                    value = sum;
                }
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
                i__2 = j + j * a_dim1;
                sum = (r__1 = a[i__2].r, f2c_abs(r__1));
                if (value < sum || sisnan_(&sum))
                {
                    value = sum;
                }
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    sum = c_abs(&a[i__ + j * a_dim1]);
                    if (value < sum || sisnan_(&sum))
                    {
                        value = sum;
                    }
                    /* L30: */
                }
                /* L40: */
            }
        }
    }
    else if (lsame_(norm, "I") || lsame_(norm, "O") || *(unsigned char *)norm == '1')
    {
        /* Find normI(A) ( = norm1(A), since A is hermitian). */
        value = 0.f;
        if (lsame_(uplo, "U"))
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                sum = 0.f;
                i__2 = j - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    absa = c_abs(&a[i__ + j * a_dim1]);
                    sum += absa;
                    work[i__] += absa;
                    /* L50: */
                }
                i__2 = j + j * a_dim1;
                work[j] = sum + (r__1 = a[i__2].r, f2c_abs(r__1));
                /* L60: */
            }
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
                work[i__] = 0.f;
                /* L80: */
            }
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j + j * a_dim1;
                sum = work[j] + (r__1 = a[i__2].r, f2c_abs(r__1));
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    absa = c_abs(&a[i__ + j * a_dim1]);
                    sum += absa;
                    work[i__] += absa;
                    /* L90: */
                }
                if (value < sum || sisnan_(&sum))
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
        scale = 0.f;
        sum = 1.f;
        if (lsame_(uplo, "U"))
        {
            i__1 = *n;
            for (j = 2;
                    j <= i__1;
                    ++j)
            {
                i__2 = j - 1;
                classq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
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
                classq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
                /* L120: */
            }
        }
        sum *= 2;
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__ + i__ * a_dim1;
            if (a[i__2].r != 0.f)
            {
                i__2 = i__ + i__ * a_dim1;
                absa = (r__1 = a[i__2].r, f2c_abs(r__1));
                if (scale < absa)
                {
                    /* Computing 2nd power */
                    r__1 = scale / absa;
                    sum = sum * (r__1 * r__1) + 1.f;
                    scale = absa;
                }
                else
                {
                    /* Computing 2nd power */
                    r__1 = absa / scale;
                    sum += r__1 * r__1;
                }
            }
            /* L130: */
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    return ret_val;
    /* End of CLANHE */
}
/* clanhe_ */
