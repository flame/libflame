/* ../netlib/clantr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele ment of largest absolute value of a trapezoidal or triangular matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLANTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clantr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clantr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clantr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLANTR( NORM, UPLO, DIAG, M, N, A, LDA, */
/* WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, NORM, UPLO */
/* INTEGER LDA, M, N */
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
/* > CLANTR returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > trapezoidal or triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANTR */
/* > \verbatim */
/* > */
/* > CLANTR = ( max(f2c_abs(A(i,j))), NORM = 'M' or 'm' */
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
/* > Specifies the value to be returned in CLANTR as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the matrix A is upper or lower trapezoidal. */
/* > = 'U': Upper trapezoidal */
/* > = 'L': Lower trapezoidal */
/* > Note that A is triangular instead of trapezoidal if M = N. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* > DIAG is CHARACTER*1 */
/* > Specifies whether or not the matrix A has unit diagonal. */
/* > = 'N': Non-unit diagonal */
/* > = 'U': Unit diagonal */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0, and if */
/* > UPLO = 'U', M <= N. When M = 0, CLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0, and if */
/* > UPLO = 'L', N <= M. When N = 0, CLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > The trapezoidal matrix A (A is triangular if M = N). */
/* > If UPLO = 'U', the leading m by n upper trapezoidal part of */
/* > the array A contains the upper trapezoidal matrix, and the */
/* > strictly lower triangular part of A is not referenced. */
/* > If UPLO = 'L', the leading m by n lower trapezoidal part of */
/* > the array A contains the lower trapezoidal matrix, and the */
/* > strictly upper triangular part of A is not referenced. Note */
/* > that when DIAG = 'U', the diagonal elements of A are not */
/* > referenced and are assumed to be one. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(M,1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)), */
/* > where LWORK >= M when NORM = 'I';
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
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
real clantr_(char *norm, char *uplo, char *diag, integer *m, integer *n, complex *a, integer *lda, real *work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real ret_val;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    integer i__, j;
    real sum, scale;
    logical udiag;
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
    if (min(*m,*n) == 0)
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
                    /* Computing MIN */
                    i__3 = *m;
                    i__4 = j - 1; // , expr subst
                    i__2 = min(i__3,i__4);
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
                    i__2 = *m;
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
                    i__2 = min(*m,j);
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum = c_abs(&a[i__ + j * a_dim1]);
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
                    i__2 = *m;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum = c_abs(&a[i__ + j * a_dim1]);
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
                if (udiag && j <= *m)
                {
                    sum = 1.f;
                    i__2 = j - 1;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += c_abs(&a[i__ + j * a_dim1]);
                        /* L90: */
                    }
                }
                else
                {
                    sum = 0.f;
                    i__2 = min(*m,j);
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += c_abs(&a[i__ + j * a_dim1]);
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
                    i__2 = *m;
                    for (i__ = j + 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += c_abs(&a[i__ + j * a_dim1]);
                        /* L120: */
                    }
                }
                else
                {
                    sum = 0.f;
                    i__2 = *m;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        sum += c_abs(&a[i__ + j * a_dim1]);
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
        if (lsame_(uplo, "U"))
        {
            if (lsame_(diag, "U"))
            {
                i__1 = *m;
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
                    /* Computing MIN */
                    i__3 = *m;
                    i__4 = j - 1; // , expr subst
                    i__2 = min(i__3,i__4);
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        work[i__] += c_abs(&a[i__ + j * a_dim1]);
                        /* L160: */
                    }
                    /* L170: */
                }
            }
            else
            {
                i__1 = *m;
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
                    i__2 = min(*m,j);
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        work[i__] += c_abs(&a[i__ + j * a_dim1]);
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
                i__1 = *m;
                for (i__ = *n + 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] = 0.f;
                    /* L220: */
                }
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *m;
                    for (i__ = j + 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        work[i__] += c_abs(&a[i__ + j * a_dim1]);
                        /* L230: */
                    }
                    /* L240: */
                }
            }
            else
            {
                i__1 = *m;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    work[i__] = 0.f;
                    /* L250: */
                }
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *m;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        work[i__] += c_abs(&a[i__ + j * a_dim1]);
                        /* L260: */
                    }
                    /* L270: */
                }
            }
        }
        value = 0.f;
        i__1 = *m;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            sum = work[i__];
            if (value < sum || sisnan_(&sum))
            {
                value = sum;
            }
            /* L280: */
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
                sum = (real) min(*m,*n);
                i__1 = *n;
                for (j = 2;
                        j <= i__1;
                        ++j)
                {
                    /* Computing MIN */
                    i__3 = *m;
                    i__4 = j - 1; // , expr subst
                    i__2 = min(i__3,i__4);
                    classq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
                    /* L290: */
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
                    i__2 = min(*m,j);
                    classq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
                    /* L300: */
                }
            }
        }
        else
        {
            if (lsame_(diag, "U"))
            {
                scale = 1.f;
                sum = (real) min(*m,*n);
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *m - j;
                    /* Computing MIN */
                    i__3 = *m;
                    i__4 = j + 1; // , expr subst
                    classq_(&i__2, &a[min(i__3,i__4) + j * a_dim1], &c__1, & scale, &sum);
                    /* L310: */
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
                    i__2 = *m - j + 1;
                    classq_(&i__2, &a[j + j * a_dim1], &c__1, &scale, &sum);
                    /* L320: */
                }
            }
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    return ret_val;
    /* End of CLANTR */
}
/* clantr_ */
