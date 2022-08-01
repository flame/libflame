/* clanhs.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of an upper Hessenberg matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLANHS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanhs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanhs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanhs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLANHS( NORM, N, A, LDA, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM */
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
/* > CLANHS returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > Hessenberg matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANHS */
/* > \verbatim */
/* > */
/* > CLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
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
/* > squares). Note that max(abs(A(i,j))) is not a consistent matrix norm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] NORM */
/* > \verbatim */
/* > NORM is CHARACTER*1 */
/* > Specifies the value to be returned in CLANHS as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. When N = 0, CLANHS is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > The n by n upper Hessenberg matrix A;
the part of A below the */
/* > first sub-diagonal is not referenced. */
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
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
real clanhs_(char *norm, integer *n, complex *a, integer *lda, real *work)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"clanhs inputs: norm %c, n %lld, lda %lld",*norm, *n, *lda);
#else
    snprintf(buffer, 256,"clanhs inputs: norm %c, n %d, lda %d",*norm, *n, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real ret_val;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    integer i__, j;
    real sum, scale;
    extern logical lsame_(char *, char *);
    real value;
    extern /* Subroutine */
    int classq_(integer *, complex *, integer *, real *, real *);
    extern logical sisnan_(real *);
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
        /* Find max(abs(A(i,j))). */
        value = 0.f;
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Computing MIN */
            i__3 = *n;
            i__4 = j + 1; // , expr subst
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
    else if (lsame_(norm, "O") || *(unsigned char *) norm == '1')
    {
        /* Find norm1(A). */
        value = 0.f;
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            sum = 0.f;
            /* Computing MIN */
            i__3 = *n;
            i__4 = j + 1; // , expr subst
            i__2 = min(i__3,i__4);
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                sum += c_abs(&a[i__ + j * a_dim1]);
                /* L30: */
            }
            if (value < sum || sisnan_(&sum))
            {
                value = sum;
            }
            /* L40: */
        }
    }
    else if (lsame_(norm, "I"))
    {
        /* Find normI(A). */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            work[i__] = 0.f;
            /* L50: */
        }
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Computing MIN */
            i__3 = *n;
            i__4 = j + 1; // , expr subst
            i__2 = min(i__3,i__4);
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                work[i__] += c_abs(&a[i__ + j * a_dim1]);
                /* L60: */
            }
            /* L70: */
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
            /* L80: */
        }
    }
    else if (lsame_(norm, "F") || lsame_(norm, "E"))
    {
        /* Find normF(A). */
        scale = 0.f;
        sum = 1.f;
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Computing MIN */
            i__3 = *n;
            i__4 = j + 1; // , expr subst
            i__2 = min(i__3,i__4);
            classq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
            /* L90: */
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return ret_val;
    /* End of CLANHS */
}
/* clanhs_ */
