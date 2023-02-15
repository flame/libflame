/* ../netlib/slantb.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SLANTB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele ment of largest absolute value of a triangular band matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLANTB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantb. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantb. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantb. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION SLANTB( NORM, UPLO, DIAG, N, K, AB, */
/* LDAB, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORM, UPLO */
/* INTEGER K, LDAB, N */
/* .. */
/* .. Array Arguments .. */
/* REAL AB( LDAB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANTB returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of an */
/* > n by n triangular band matrix A, with ( k + 1 ) diagonals. */
/* > \endverbatim */
/* > */
/* > \return SLANTB */
/* > \verbatim */
/* > */
/* > SLANTB = ( max(f2c_abs(A(i,j))), NORM = 'M' or 'm' */
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
/* > Specifies the value to be returned in SLANTB as described */
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
/* > The order of the matrix A. N >= 0. When N = 0, SLANTB is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of super-diagonals of the matrix A if UPLO = 'U', */
/* > or the number of sub-diagonals of the matrix A if UPLO = 'L'. */
/* > K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is REAL array, dimension (LDAB,N) */
/* > The upper or lower triangular band matrix A, stored in the */
/* > first k+1 rows of AB. The j-th column of A is stored */
/* > in the j-th column of the array AB as follows: */
/* > if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
*/
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+k). */
/* > Note that when DIAG = 'U', the elements of the array AB */
/* > corresponding to the diagonal elements of the matrix A are */
/* > not referenced, but are assumed to be one. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= K+1. */
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
real slantb_(char *norm, char *uplo, char *diag, integer *n, integer *k, real *ab, integer *ldab, real *work)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
    real ret_val, r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, l;
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
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --work;
    /* Function Body */
    if (*n == 0)
    {
        value = 0.f;
    }
    else if (lsame_(norm, "M"))
    {
        /* Find max(f2c_abs(A(i,j))). */
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
                    /* Computing MAX */
                    i__2 = *k + 2 - j;
                    i__3 = *k;
                    for (i__ = max(i__2,1);
                            i__ <= i__3;
                            ++i__)
                    {
                        sum = (r__1 = ab[i__ + j * ab_dim1], f2c_abs(r__1));
                        if (value < sum || sisnan_(&sum))
                        {
                            value = sum;
                        }
                        /* L10: */
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
                    /* Computing MIN */
                    i__2 = *n + 1 - j;
                    i__4 = *k + 1; // , expr subst
                    i__3 = min(i__2,i__4);
                    for (i__ = 2;
                            i__ <= i__3;
                            ++i__)
                    {
                        sum = (r__1 = ab[i__ + j * ab_dim1], f2c_abs(r__1));
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
                    /* Computing MAX */
                    i__3 = *k + 2 - j;
                    i__2 = *k + 1;
                    for (i__ = max(i__3,1);
                            i__ <= i__2;
                            ++i__)
                    {
                        sum = (r__1 = ab[i__ + j * ab_dim1], f2c_abs(r__1));
                        if (value < sum || sisnan_(&sum))
                        {
                            value = sum;
                        }
                        /* L50: */
                    }
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
                    /* Computing MIN */
                    i__3 = *n + 1 - j;
                    i__4 = *k + 1; // , expr subst
                    i__2 = min(i__3,i__4);
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum = (r__1 = ab[i__ + j * ab_dim1], f2c_abs(r__1));
                        if (value < sum || sisnan_(&sum))
                        {
                            value = sum;
                        }
                        /* L70: */
                    }
                    /* L80: */
                }
            }
        }
    }
    else if (lsame_(norm, "O") || *(unsigned char *) norm == '1')
    {
        /* Find norm1(A). */
        value = 0.f;
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
                    /* Computing MAX */
                    i__2 = *k + 2 - j;
                    i__3 = *k;
                    for (i__ = max(i__2,1);
                            i__ <= i__3;
                            ++i__)
                    {
                        sum += (r__1 = ab[i__ + j * ab_dim1], f2c_abs(r__1));
                        /* L90: */
                    }
                }
                else
                {
                    sum = 0.f;
                    /* Computing MAX */
                    i__3 = *k + 2 - j;
                    i__2 = *k + 1;
                    for (i__ = max(i__3,1);
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += (r__1 = ab[i__ + j * ab_dim1], f2c_abs(r__1));
                        /* L100: */
                    }
                }
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
                    /* Computing MIN */
                    i__3 = *n + 1 - j;
                    i__4 = *k + 1; // , expr subst
                    i__2 = min(i__3,i__4);
                    for (i__ = 2;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += (r__1 = ab[i__ + j * ab_dim1], f2c_abs(r__1));
                        /* L120: */
                    }
                }
                else
                {
                    sum = 0.f;
                    /* Computing MIN */
                    i__3 = *n + 1 - j;
                    i__4 = *k + 1; // , expr subst
                    i__2 = min(i__3,i__4);
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += (r__1 = ab[i__ + j * ab_dim1], f2c_abs(r__1));
                        /* L130: */
                    }
                }
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
        value = 0.f;
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
                    l = *k + 1 - j;
                    /* Computing MAX */
                    i__2 = 1;
                    i__3 = j - *k; // , expr subst
                    i__4 = j - 1;
                    for (i__ = max(i__2,i__3);
                            i__ <= i__4;
                            ++i__)
                    {
                        work[i__] += (r__1 = ab[l + i__ + j * ab_dim1], f2c_abs( r__1));
                        /* L160: */
                    }
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
                    l = *k + 1 - j;
                    /* Computing MAX */
                    i__4 = 1;
                    i__2 = j - *k; // , expr subst
                    i__3 = j;
                    for (i__ = max(i__4,i__2);
                            i__ <= i__3;
                            ++i__)
                    {
                        work[i__] += (r__1 = ab[l + i__ + j * ab_dim1], f2c_abs( r__1));
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
                    l = 1 - j;
                    /* Computing MIN */
                    i__4 = *n;
                    i__2 = j + *k; // , expr subst
                    i__3 = min(i__4,i__2);
                    for (i__ = j + 1;
                            i__ <= i__3;
                            ++i__)
                    {
                        work[i__] += (r__1 = ab[l + i__ + j * ab_dim1], f2c_abs( r__1));
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
                    l = 1 - j;
                    /* Computing MIN */
                    i__4 = *n;
                    i__2 = j + *k; // , expr subst
                    i__3 = min(i__4,i__2);
                    for (i__ = j;
                            i__ <= i__3;
                            ++i__)
                    {
                        work[i__] += (r__1 = ab[l + i__ + j * ab_dim1], f2c_abs( r__1));
                        /* L250: */
                    }
                    /* L260: */
                }
            }
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
                if (*k > 0)
                {
                    i__1 = *n;
                    for (j = 2;
                            j <= i__1;
                            ++j)
                    {
                        /* Computing MIN */
                        i__4 = j - 1;
                        i__3 = min(i__4,*k);
                        /* Computing MAX */
                        i__2 = *k + 2 - j;
                        slassq_(&i__3, &ab[max(i__2,1) + j * ab_dim1], &c__1, &scale, &sum);
                        /* L280: */
                    }
                }
            }
            else
            {
                scale = 0.f;
                sum = 1.f;
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    /* Computing MIN */
                    i__4 = j;
                    i__2 = *k + 1; // , expr subst
                    i__3 = min(i__4,i__2);
                    /* Computing MAX */
                    i__5 = *k + 2 - j;
                    slassq_(&i__3, &ab[max(i__5,1) + j * ab_dim1], &c__1, & scale, &sum);
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
                if (*k > 0)
                {
                    i__1 = *n - 1;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        /* Computing MIN */
                        i__4 = *n - j;
                        i__3 = min(i__4,*k);
                        slassq_(&i__3, &ab[j * ab_dim1 + 2], &c__1, &scale, & sum);
                        /* L300: */
                    }
                }
            }
            else
            {
                scale = 0.f;
                sum = 1.f;
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    /* Computing MIN */
                    i__4 = *n - j + 1;
                    i__2 = *k + 1; // , expr subst
                    i__3 = min(i__4,i__2);
                    slassq_(&i__3, &ab[j * ab_dim1 + 1], &c__1, &scale, &sum);
                    /* L310: */
                }
            }
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    return ret_val;
    /* End of SLANTB */
}
/* slantb_ */
