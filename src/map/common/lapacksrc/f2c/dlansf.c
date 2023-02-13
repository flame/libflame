/* ../netlib/dlansf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DLANSF returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele ment of largest absolute value of a symmetric matrix in RFP format. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLANSF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlansf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlansf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlansf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* DOUBLE PRECISION FUNCTION DLANSF( NORM, TRANSR, UPLO, N, A, WORK ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM, TRANSR, UPLO */
/* INTEGER N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( 0: * ), WORK( 0: * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANSF returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > real symmetric matrix A in RFP format. */
/* > \endverbatim */
/* > */
/* > \return DLANSF */
/* > \verbatim */
/* > */
/* > DLANSF = ( max(f2c_abs(A(i,j))), NORM = 'M' or 'm' */
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
/* > squares). Note that max(f2c_abs(A(i,j))) is not a matrix norm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] NORM */
/* > \verbatim */
/* > NORM is CHARACTER*1 */
/* > Specifies the value to be returned in DLANSF as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSR */
/* > \verbatim */
/* > TRANSR is CHARACTER*1 */
/* > Specifies whether the RFP format of A is normal or */
/* > transposed format. */
/* > = 'N': RFP format is Normal;
*/
/* > = 'T': RFP format is Transpose. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > On entry, UPLO specifies whether the RFP matrix A came from */
/* > an upper or lower triangular matrix as follows: */
/* > = 'U': RFP A came from an upper triangular matrix;
*/
/* > = 'L': RFP A came from a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. When N = 0, DLANSF is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension ( N*(N+1)/2 );
*/
/* > On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') */
/* > part of the symmetric matrix A stored in RFP format. See the */
/* > "Notes" below for more details. */
/* > Unchanged on exit. */
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
/* > \ingroup doubleOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > We first consider Rectangular Full Packed (RFP) Format when N is */
/* > even. We give an example where N = 6. */
/* > */
/* > AP is Upper AP is Lower */
/* > */
/* > 00 01 02 03 04 05 00 */
/* > 11 12 13 14 15 10 11 */
/* > 22 23 24 25 20 21 22 */
/* > 33 34 35 30 31 32 33 */
/* > 44 45 40 41 42 43 44 */
/* > 55 50 51 52 53 54 55 */
/* > */
/* > */
/* > Let TRANSR = 'N'. RFP holds AP as follows: */
/* > For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last */
/* > three columns of AP upper. The lower triangle A(4:6,0:2) consists of */
/* > the transpose of the first three columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:2,0:2) consists of */
/* > the transpose of the last three columns of AP lower. */
/* > This covers the case N even and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > 03 04 05 33 43 53 */
/* > 13 14 15 00 44 54 */
/* > 23 24 25 10 11 55 */
/* > 33 34 35 20 21 22 */
/* > 00 44 45 30 31 32 */
/* > 01 11 55 40 41 42 */
/* > 02 12 22 50 51 52 */
/* > */
/* > Now let TRANSR = 'T'. RFP A in both UPLO cases is just the */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* > RFP A RFP A */
/* > */
/* > 03 13 23 33 00 01 02 33 00 10 20 30 40 50 */
/* > 04 14 24 34 44 11 12 43 44 11 21 31 41 51 */
/* > 05 15 25 35 45 55 22 53 54 55 22 32 42 52 */
/* > */
/* > */
/* > We then consider Rectangular Full Packed (RFP) Format when N is */
/* > odd. We give an example where N = 5. */
/* > */
/* > AP is Upper AP is Lower */
/* > */
/* > 00 01 02 03 04 00 */
/* > 11 12 13 14 10 11 */
/* > 22 23 24 20 21 22 */
/* > 33 34 30 31 32 33 */
/* > 44 40 41 42 43 44 */
/* > */
/* > */
/* > Let TRANSR = 'N'. RFP holds AP as follows: */
/* > For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last */
/* > three columns of AP upper. The lower triangle A(3:4,0:1) consists of */
/* > the transpose of the first two columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:1,1:2) consists of */
/* > the transpose of the last two columns of AP lower. */
/* > This covers the case N odd and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > 02 03 04 00 33 43 */
/* > 12 13 14 10 11 44 */
/* > 22 23 24 20 21 22 */
/* > 00 33 34 30 31 32 */
/* > 01 11 44 40 41 42 */
/* > */
/* > Now let TRANSR = 'T'. RFP A in both UPLO cases is just the */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > RFP A RFP A */
/* > */
/* > 02 12 22 00 01 00 10 20 30 40 50 */
/* > 03 13 23 33 11 33 11 21 31 41 51 */
/* > 04 14 24 34 44 43 44 22 32 42 52 */
/* > \endverbatim */
/* ===================================================================== */
doublereal dlansf_(char *norm, char *transr, char *uplo, integer *n, doublereal *a, doublereal *work)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k, l;
    doublereal s;
    integer n1;
    doublereal aa;
    integer lda, ifm, noe, ilu;
    doublereal temp, scale;
    extern logical lsame_(char *, char *);
    doublereal value;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */
    int dlassq_(integer *, doublereal *, integer *, doublereal *, doublereal *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    if (*n == 0)
    {
        ret_val = 0.;
        return ret_val;
    }
    else if (*n == 1)
    {
        ret_val = f2c_abs(a[0]);
        return ret_val;
    }
    /* set noe = 1 if n is odd. if n is even set noe=0 */
    noe = 1;
    if (*n % 2 == 0)
    {
        noe = 0;
    }
    /* set ifm = 0 when form='T or 't' and 1 otherwise */
    ifm = 1;
    if (lsame_(transr, "T"))
    {
        ifm = 0;
    }
    /* set ilu = 0 when uplo='U or 'u' and 1 otherwise */
    ilu = 1;
    if (lsame_(uplo, "U"))
    {
        ilu = 0;
    }
    /* set lda = (n+1)/2 when ifm = 0 */
    /* set lda = n when ifm = 1 and noe = 1 */
    /* set lda = n+1 when ifm = 1 and noe = 0 */
    if (ifm == 1)
    {
        if (noe == 1)
        {
            lda = *n;
        }
        else
        {
            /* noe=0 */
            lda = *n + 1;
        }
    }
    else
    {
        /* ifm=0 */
        lda = (*n + 1) / 2;
    }
    if (lsame_(norm, "M"))
    {
        /* Find max(f2c_abs(A(i,j))). */
        k = (*n + 1) / 2;
        value = 0.;
        if (noe == 1)
        {
            /* n is odd */
            if (ifm == 1)
            {
                /* A is n by k */
                i__1 = k - 1;
                for (j = 0;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *n - 1;
                    for (i__ = 0;
                            i__ <= i__2;
                            ++i__)
                    {
                        temp = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
            }
            else
            {
                /* xpose case;
                A is k by n */
                i__1 = *n - 1;
                for (j = 0;
                        j <= i__1;
                        ++j)
                {
                    i__2 = k - 1;
                    for (i__ = 0;
                            i__ <= i__2;
                            ++i__)
                    {
                        temp = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
            }
        }
        else
        {
            /* n is even */
            if (ifm == 1)
            {
                /* A is n+1 by k */
                i__1 = k - 1;
                for (j = 0;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *n;
                    for (i__ = 0;
                            i__ <= i__2;
                            ++i__)
                    {
                        temp = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
            }
            else
            {
                /* xpose case;
                A is k by n+1 */
                i__1 = *n;
                for (j = 0;
                        j <= i__1;
                        ++j)
                {
                    i__2 = k - 1;
                    for (i__ = 0;
                            i__ <= i__2;
                            ++i__)
                    {
                        temp = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
            }
        }
    }
    else if (lsame_(norm, "I") || lsame_(norm, "O") || *(unsigned char *)norm == '1')
    {
        /* Find normI(A) ( = norm1(A), since A is symmetric). */
        if (ifm == 1)
        {
            k = *n / 2;
            if (noe == 1)
            {
                /* n is odd */
                if (ilu == 0)
                {
                    i__1 = k - 1;
                    for (i__ = 0;
                            i__ <= i__1;
                            ++i__)
                    {
                        work[i__] = 0.;
                    }
                    i__1 = k;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        s = 0.;
                        i__2 = k + j - 1;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* -> A(i,j+k) */
                            s += aa;
                            work[i__] += aa;
                        }
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* -> A(j+k,j+k) */
                        work[j + k] = s + aa;
                        if (i__ == k + k)
                        {
                            goto L10;
                        }
                        ++i__;
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* -> A(j,j) */
                        work[j] += aa;
                        s = 0.;
                        i__2 = k - 1;
                        for (l = j + 1;
                                l <= i__2;
                                ++l)
                        {
                            ++i__;
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* -> A(l,j) */
                            s += aa;
                            work[l] += aa;
                        }
                        work[j] += s;
                    }
L10:
                    value = work[0];
                    i__1 = *n - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        temp = work[i__];
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
                else
                {
                    /* ilu = 1 */
                    ++k;
                    /* k=(n+1)/2 for n odd and ilu=1 */
                    i__1 = *n - 1;
                    for (i__ = k;
                            i__ <= i__1;
                            ++i__)
                    {
                        work[i__] = 0.;
                    }
                    for (j = k - 1;
                            j >= 0;
                            --j)
                    {
                        s = 0.;
                        i__1 = j - 2;
                        for (i__ = 0;
                                i__ <= i__1;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* -> A(j+k,i+k) */
                            s += aa;
                            work[i__ + k] += aa;
                        }
                        if (j > 0)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* -> A(j+k,j+k) */
                            s += aa;
                            work[i__ + k] += s;
                            /* i=j */
                            ++i__;
                        }
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* -> A(j,j) */
                        work[j] = aa;
                        s = 0.;
                        i__1 = *n - 1;
                        for (l = j + 1;
                                l <= i__1;
                                ++l)
                        {
                            ++i__;
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* -> A(l,j) */
                            s += aa;
                            work[l] += aa;
                        }
                        work[j] += s;
                    }
                    value = work[0];
                    i__1 = *n - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        temp = work[i__];
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
            }
            else
            {
                /* n is even */
                if (ilu == 0)
                {
                    i__1 = k - 1;
                    for (i__ = 0;
                            i__ <= i__1;
                            ++i__)
                    {
                        work[i__] = 0.;
                    }
                    i__1 = k - 1;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        s = 0.;
                        i__2 = k + j - 1;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* -> A(i,j+k) */
                            s += aa;
                            work[i__] += aa;
                        }
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* -> A(j+k,j+k) */
                        work[j + k] = s + aa;
                        ++i__;
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* -> A(j,j) */
                        work[j] += aa;
                        s = 0.;
                        i__2 = k - 1;
                        for (l = j + 1;
                                l <= i__2;
                                ++l)
                        {
                            ++i__;
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* -> A(l,j) */
                            s += aa;
                            work[l] += aa;
                        }
                        work[j] += s;
                    }
                    value = work[0];
                    i__1 = *n - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        temp = work[i__];
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
                else
                {
                    /* ilu = 1 */
                    i__1 = *n - 1;
                    for (i__ = k;
                            i__ <= i__1;
                            ++i__)
                    {
                        work[i__] = 0.;
                    }
                    for (j = k - 1;
                            j >= 0;
                            --j)
                    {
                        s = 0.;
                        i__1 = j - 1;
                        for (i__ = 0;
                                i__ <= i__1;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* -> A(j+k,i+k) */
                            s += aa;
                            work[i__ + k] += aa;
                        }
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* -> A(j+k,j+k) */
                        s += aa;
                        work[i__ + k] += s;
                        /* i=j */
                        ++i__;
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* -> A(j,j) */
                        work[j] = aa;
                        s = 0.;
                        i__1 = *n - 1;
                        for (l = j + 1;
                                l <= i__1;
                                ++l)
                        {
                            ++i__;
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* -> A(l,j) */
                            s += aa;
                            work[l] += aa;
                        }
                        work[j] += s;
                    }
                    value = work[0];
                    i__1 = *n - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        temp = work[i__];
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
            }
        }
        else
        {
            /* ifm=0 */
            k = *n / 2;
            if (noe == 1)
            {
                /* n is odd */
                if (ilu == 0)
                {
                    n1 = k;
                    /* n/2 */
                    ++k;
                    /* k is the row size and lda */
                    i__1 = *n - 1;
                    for (i__ = n1;
                            i__ <= i__1;
                            ++i__)
                    {
                        work[i__] = 0.;
                    }
                    i__1 = n1 - 1;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        s = 0.;
                        i__2 = k - 1;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(j,n1+i) */
                            work[i__ + n1] += aa;
                            s += aa;
                        }
                        work[j] = s;
                    }
                    /* j=n1=k-1 is special */
                    s = (d__1 = a[j * lda], f2c_abs(d__1));
                    /* A(k-1,k-1) */
                    i__1 = k - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* A(k-1,i+n1) */
                        work[i__ + n1] += aa;
                        s += aa;
                    }
                    work[j] += s;
                    i__1 = *n - 1;
                    for (j = k;
                            j <= i__1;
                            ++j)
                    {
                        s = 0.;
                        i__2 = j - k - 1;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(i,j-k) */
                            work[i__] += aa;
                            s += aa;
                        }
                        /* i=j-k */
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* A(j-k,j-k) */
                        s += aa;
                        work[j - k] += s;
                        ++i__;
                        s = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* A(j,j) */
                        i__2 = *n - 1;
                        for (l = j + 1;
                                l <= i__2;
                                ++l)
                        {
                            ++i__;
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(j,l) */
                            work[l] += aa;
                            s += aa;
                        }
                        work[j] += s;
                    }
                    value = work[0];
                    i__1 = *n - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        temp = work[i__];
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
                else
                {
                    /* ilu=1 */
                    ++k;
                    /* k=(n+1)/2 for n odd and ilu=1 */
                    i__1 = *n - 1;
                    for (i__ = k;
                            i__ <= i__1;
                            ++i__)
                    {
                        work[i__] = 0.;
                    }
                    i__1 = k - 2;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        /* process */
                        s = 0.;
                        i__2 = j - 1;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(j,i) */
                            work[i__] += aa;
                            s += aa;
                        }
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* i=j so process of A(j,j) */
                        s += aa;
                        work[j] = s;
                        /* is initialised here */
                        ++i__;
                        /* i=j process A(j+k,j+k) */
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        s = aa;
                        i__2 = *n - 1;
                        for (l = k + j + 1;
                                l <= i__2;
                                ++l)
                        {
                            ++i__;
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(l,k+j) */
                            s += aa;
                            work[l] += aa;
                        }
                        work[k + j] += s;
                    }
                    /* j=k-1 is special :process col A(k-1,0:k-1) */
                    s = 0.;
                    i__1 = k - 2;
                    for (i__ = 0;
                            i__ <= i__1;
                            ++i__)
                    {
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* A(k,i) */
                        work[i__] += aa;
                        s += aa;
                    }
                    /* i=k-1 */
                    aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                    /* A(k-1,k-1) */
                    s += aa;
                    work[i__] = s;
                    /* done with col j=k+1 */
                    i__1 = *n - 1;
                    for (j = k;
                            j <= i__1;
                            ++j)
                    {
                        /* process col j of A = A(j,0:k-1) */
                        s = 0.;
                        i__2 = k - 1;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(j,i) */
                            work[i__] += aa;
                            s += aa;
                        }
                        work[j] += s;
                    }
                    value = work[0];
                    i__1 = *n - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        temp = work[i__];
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
            }
            else
            {
                /* n is even */
                if (ilu == 0)
                {
                    i__1 = *n - 1;
                    for (i__ = k;
                            i__ <= i__1;
                            ++i__)
                    {
                        work[i__] = 0.;
                    }
                    i__1 = k - 1;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        s = 0.;
                        i__2 = k - 1;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(j,i+k) */
                            work[i__ + k] += aa;
                            s += aa;
                        }
                        work[j] = s;
                    }
                    /* j=k */
                    aa = (d__1 = a[j * lda], f2c_abs(d__1));
                    /* A(k,k) */
                    s = aa;
                    i__1 = k - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* A(k,k+i) */
                        work[i__ + k] += aa;
                        s += aa;
                    }
                    work[j] += s;
                    i__1 = *n - 1;
                    for (j = k + 1;
                            j <= i__1;
                            ++j)
                    {
                        s = 0.;
                        i__2 = j - 2 - k;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(i,j-k-1) */
                            work[i__] += aa;
                            s += aa;
                        }
                        /* i=j-1-k */
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* A(j-k-1,j-k-1) */
                        s += aa;
                        work[j - k - 1] += s;
                        ++i__;
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* A(j,j) */
                        s = aa;
                        i__2 = *n - 1;
                        for (l = j + 1;
                                l <= i__2;
                                ++l)
                        {
                            ++i__;
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(j,l) */
                            work[l] += aa;
                            s += aa;
                        }
                        work[j] += s;
                    }
                    /* j=n */
                    s = 0.;
                    i__1 = k - 2;
                    for (i__ = 0;
                            i__ <= i__1;
                            ++i__)
                    {
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* A(i,k-1) */
                        work[i__] += aa;
                        s += aa;
                    }
                    /* i=k-1 */
                    aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                    /* A(k-1,k-1) */
                    s += aa;
                    work[i__] += s;
                    value = work[0];
                    i__1 = *n - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        temp = work[i__];
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
                else
                {
                    /* ilu=1 */
                    i__1 = *n - 1;
                    for (i__ = k;
                            i__ <= i__1;
                            ++i__)
                    {
                        work[i__] = 0.;
                    }
                    /* j=0 is special :process col A(k:n-1,k) */
                    s = f2c_abs(a[0]);
                    /* A(k,k) */
                    i__1 = k - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        aa = (d__1 = a[i__], f2c_abs(d__1));
                        /* A(k+i,k) */
                        work[i__ + k] += aa;
                        s += aa;
                    }
                    work[k] += s;
                    i__1 = k - 1;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        /* process */
                        s = 0.;
                        i__2 = j - 2;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(j-1,i) */
                            work[i__] += aa;
                            s += aa;
                        }
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* i=j-1 so process of A(j-1,j-1) */
                        s += aa;
                        work[j - 1] = s;
                        /* is initialised here */
                        ++i__;
                        /* i=j process A(j+k,j+k) */
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        s = aa;
                        i__2 = *n - 1;
                        for (l = k + j + 1;
                                l <= i__2;
                                ++l)
                        {
                            ++i__;
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(l,k+j) */
                            s += aa;
                            work[l] += aa;
                        }
                        work[k + j] += s;
                    }
                    /* j=k is special :process col A(k,0:k-1) */
                    s = 0.;
                    i__1 = k - 2;
                    for (i__ = 0;
                            i__ <= i__1;
                            ++i__)
                    {
                        aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                        /* A(k,i) */
                        work[i__] += aa;
                        s += aa;
                    }
                    /* i=k-1 */
                    aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                    /* A(k-1,k-1) */
                    s += aa;
                    work[i__] = s;
                    /* done with col j=k+1 */
                    i__1 = *n;
                    for (j = k + 1;
                            j <= i__1;
                            ++j)
                    {
                        /* process col j-1 of A = A(j-1,0:k-1) */
                        s = 0.;
                        i__2 = k - 1;
                        for (i__ = 0;
                                i__ <= i__2;
                                ++i__)
                        {
                            aa = (d__1 = a[i__ + j * lda], f2c_abs(d__1));
                            /* A(j-1,i) */
                            work[i__] += aa;
                            s += aa;
                        }
                        work[j - 1] += s;
                    }
                    value = work[0];
                    i__1 = *n - 1;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        temp = work[i__];
                        if (value < temp || disnan_(&temp))
                        {
                            value = temp;
                        }
                    }
                }
            }
        }
    }
    else if (lsame_(norm, "F") || lsame_(norm, "E"))
    {
        /* Find normF(A). */
        k = (*n + 1) / 2;
        scale = 0.;
        s = 1.;
        if (noe == 1)
        {
            /* n is odd */
            if (ifm == 1)
            {
                /* A is normal */
                if (ilu == 0)
                {
                    /* A is upper */
                    i__1 = k - 3;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = k - j - 2;
                        dlassq_(&i__2, &a[k + j + 1 + j * lda], &c__1, &scale, &s);
                        /* L at A(k,0) */
                    }
                    i__1 = k - 1;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = k + j - 1;
                        dlassq_(&i__2, &a[j * lda], &c__1, &scale, &s);
                        /* trap U at A(0,0) */
                    }
                    s += s;
                    /* double s for the off diagonal elements */
                    i__1 = k - 1;
                    i__2 = lda + 1;
                    dlassq_(&i__1, &a[k], &i__2, &scale, &s);
                    /* tri L at A(k,0) */
                    i__1 = lda + 1;
                    dlassq_(&k, &a[k - 1], &i__1, &scale, &s);
                    /* tri U at A(k-1,0) */
                }
                else
                {
                    /* ilu=1 & A is lower */
                    i__1 = k - 1;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = *n - j - 1;
                        dlassq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s) ;
                        /* trap L at A(0,0) */
                    }
                    i__1 = k - 2;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
                        /* U at A(0,1) */
                    }
                    s += s;
                    /* double s for the off diagonal elements */
                    i__1 = lda + 1;
                    dlassq_(&k, a, &i__1, &scale, &s);
                    /* tri L at A(0,0) */
                    i__1 = k - 1;
                    i__2 = lda + 1;
                    dlassq_(&i__1, &a[lda], &i__2, &scale, &s);
                    /* tri U at A(0,1) */
                }
            }
            else
            {
                /* A is xpose */
                if (ilu == 0)
                {
                    /* A**T is upper */
                    i__1 = k - 2;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&j, &a[(k + j) * lda], &c__1, &scale, &s);
                        /* U at A(0,k) */
                    }
                    i__1 = k - 2;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&k, &a[j * lda], &c__1, &scale, &s);
                        /* k by k-1 rect. at A(0,0) */
                    }
                    i__1 = k - 2;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = k - j - 1;
                        dlassq_(&i__2, &a[j + 1 + (j + k - 1) * lda], &c__1, & scale, &s);
                        /* L at A(0,k-1) */
                    }
                    s += s;
                    /* double s for the off diagonal elements */
                    i__1 = k - 1;
                    i__2 = lda + 1;
                    dlassq_(&i__1, &a[k * lda], &i__2, &scale, &s);
                    /* tri U at A(0,k) */
                    i__1 = lda + 1;
                    dlassq_(&k, &a[(k - 1) * lda], &i__1, &scale, &s);
                    /* tri L at A(0,k-1) */
                }
                else
                {
                    /* A**T is lower */
                    i__1 = k - 1;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&j, &a[j * lda], &c__1, &scale, &s);
                        /* U at A(0,0) */
                    }
                    i__1 = *n - 1;
                    for (j = k;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&k, &a[j * lda], &c__1, &scale, &s);
                        /* k by k-1 rect. at A(0,k) */
                    }
                    i__1 = k - 3;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = k - j - 2;
                        dlassq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s) ;
                        /* L at A(1,0) */
                    }
                    s += s;
                    /* double s for the off diagonal elements */
                    i__1 = lda + 1;
                    dlassq_(&k, a, &i__1, &scale, &s);
                    /* tri U at A(0,0) */
                    i__1 = k - 1;
                    i__2 = lda + 1;
                    dlassq_(&i__1, &a[1], &i__2, &scale, &s);
                    /* tri L at A(1,0) */
                }
            }
        }
        else
        {
            /* n is even */
            if (ifm == 1)
            {
                /* A is normal */
                if (ilu == 0)
                {
                    /* A is upper */
                    i__1 = k - 2;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = k - j - 1;
                        dlassq_(&i__2, &a[k + j + 2 + j * lda], &c__1, &scale, &s);
                        /* L at A(k+1,0) */
                    }
                    i__1 = k - 1;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = k + j;
                        dlassq_(&i__2, &a[j * lda], &c__1, &scale, &s);
                        /* trap U at A(0,0) */
                    }
                    s += s;
                    /* double s for the off diagonal elements */
                    i__1 = lda + 1;
                    dlassq_(&k, &a[k + 1], &i__1, &scale, &s);
                    /* tri L at A(k+1,0) */
                    i__1 = lda + 1;
                    dlassq_(&k, &a[k], &i__1, &scale, &s);
                    /* tri U at A(k,0) */
                }
                else
                {
                    /* ilu=1 & A is lower */
                    i__1 = k - 1;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = *n - j - 1;
                        dlassq_(&i__2, &a[j + 2 + j * lda], &c__1, &scale, &s) ;
                        /* trap L at A(1,0) */
                    }
                    i__1 = k - 1;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&j, &a[j * lda], &c__1, &scale, &s);
                        /* U at A(0,0) */
                    }
                    s += s;
                    /* double s for the off diagonal elements */
                    i__1 = lda + 1;
                    dlassq_(&k, &a[1], &i__1, &scale, &s);
                    /* tri L at A(1,0) */
                    i__1 = lda + 1;
                    dlassq_(&k, a, &i__1, &scale, &s);
                    /* tri U at A(0,0) */
                }
            }
            else
            {
                /* A is xpose */
                if (ilu == 0)
                {
                    /* A**T is upper */
                    i__1 = k - 1;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&j, &a[(k + 1 + j) * lda], &c__1, &scale, &s);
                        /* U at A(0,k+1) */
                    }
                    i__1 = k - 1;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&k, &a[j * lda], &c__1, &scale, &s);
                        /* k by k rect. at A(0,0) */
                    }
                    i__1 = k - 2;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = k - j - 1;
                        dlassq_(&i__2, &a[j + 1 + (j + k) * lda], &c__1, & scale, &s);
                        /* L at A(0,k) */
                    }
                    s += s;
                    /* double s for the off diagonal elements */
                    i__1 = lda + 1;
                    dlassq_(&k, &a[(k + 1) * lda], &i__1, &scale, &s);
                    /* tri U at A(0,k+1) */
                    i__1 = lda + 1;
                    dlassq_(&k, &a[k * lda], &i__1, &scale, &s);
                    /* tri L at A(0,k) */
                }
                else
                {
                    /* A**T is lower */
                    i__1 = k - 1;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&j, &a[(j + 1) * lda], &c__1, &scale, &s);
                        /* U at A(0,1) */
                    }
                    i__1 = *n;
                    for (j = k + 1;
                            j <= i__1;
                            ++j)
                    {
                        dlassq_(&k, &a[j * lda], &c__1, &scale, &s);
                        /* k by k rect. at A(0,k+1) */
                    }
                    i__1 = k - 2;
                    for (j = 0;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = k - j - 1;
                        dlassq_(&i__2, &a[j + 1 + j * lda], &c__1, &scale, &s) ;
                        /* L at A(0,0) */
                    }
                    s += s;
                    /* double s for the off diagonal elements */
                    i__1 = lda + 1;
                    dlassq_(&k, &a[lda], &i__1, &scale, &s);
                    /* tri L at A(0,1) */
                    i__1 = lda + 1;
                    dlassq_(&k, a, &i__1, &scale, &s);
                    /* tri U at A(0,0) */
                }
            }
        }
        value = scale * sqrt(s);
    }
    ret_val = value;
    return ret_val;
    /* End of DLANSF */
}
/* dlansf_ */
